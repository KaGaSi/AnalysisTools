#include "../src/AnalysisTools.h"

// Uses Identification of the Truly Interfacial Molecules (ITIM) from
// https://doi.org/10.1002/jcc.20852

// TODO: inconsistencies in -wd/-w: should both be allowed separately? If so,
//       all needs checking as sometimes -wd is necessary for -w
// TODO: ...assumes probe + bead distance is 1 (i.e., r_c for DPD bead), that is
//       probe radius and bead radius are the same, specifically 0.5 - really?
//       what about the -r option?
// TODO: find proper test system, then implement the other stuff

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Surface utility determines the first bead (either going from the box's "
  "centre - defaault behaviour - or from its edges - if --in option is used) "
  "in each square prism defined by the given <width> parameter, thus defining "
  "a surface of, e.g., polymer brush or lipid bilayer. The <width> 'slices' "
  "the box into square prisms along the chosen axis (i.e., if z is the "
  "chosen axis, the xy plane is chopped into squares, creating "
  "<width>*<width>*<box length in z> prisms). In each such prism, two beads "
  "are found corresponding to the two surfaces (e.g., polymer brush on both "
  "box edges or the two surfaces of a lipid bilayer inside the box). The "
  "surface coordinates are saved into the <surf.txt> output file, while their "
  "area (calculated from triangles defined by the surface coordinates) are "
  "written into the <area.txt>. What bead types are considered as possible "
  "surface beads is controlled via --bonded and -bt options.",

  "Usage: Surface <input> <width> <surf.txt> <axis> [options]",
  .args = 4, // number of mandatory arguments
  .all = 21, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<width>", nullptr, "width of a single bin", OPT_ARG},
  {"<surf.txt>", nullptr, "average surface", OPT_ARG},
  {"<axis>", nullptr, "calculate along x, y, or z axis", OPT_ARG},
  {"--in", nullptr, "start from the box's edges instead of its centre",
    OPT_EXTRA},
  {"--bonded", nullptr, "use only beads in molecules", OPT_EXTRA},
  {"-bt", "<name(s)>", "bead type(s) to use", OPT_EXTRA},
  {"-wd", "<file> <w>", "calculate distribution of widths with given single "
    "bin width", OPT_EXTRA},
  {"-w", "<file>", "save per-timestep width", OPT_EXTRA},
  {"-a", "<area.txt>", "per-timestep areas", OPT_EXTRA},
  {"-b", "<file>", "save per-timestep surface beads to a coordinate file",
    OPT_EXTRA},
  {"-r", "<float>", "radius of the ITIM probe (default 0.5)", OPT_EXTRA},
  // {"-m", "<mol(s)>", "molecule type(s) to use", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int bt_number, *bt;     // -bt (number of types; list of the types)
  double distr_width,     // -wd (width of distribution bin)
         probe;           // -r (probe radius)
  char width_distr[LINE], // -wd (filename)
       width_avg[LINE],   // -w (filename)
       area_file[LINE];   // -a (filename)
  bool in,                // --in
       bonded;            // --bonded
  FILE_TYPE bead_file;    // -b (filename)
}; //}}}

// calculate area of a triangle given three points (Heron's formula) //{{{
double calc_area(const double A[3], const double B[3], const double C[3]) {
  // triangle's sides vectors
  vec3d AB, AC, BC;
  for (int dd = 0; dd < 3; dd++) {
    AB.v[dd] = B[dd] - A[dd];
    AC.v[dd] = C[dd] - A[dd];
    BC.v[dd] = C[dd] - B[dd];
  }
  double a = 0, b = 0, c = 0; // triangle's sidelengths
  for (int dd = 0; dd < 3; dd++) {
    a = VectLength(BC);
    b = VectLength(AC);
    c = VectLength(AB);
  }
  double s = (a + b + c) / 2;
  double area_sq = s * (s - a) * (s - b) * (s - c);
  if (area_sq < 0) { // near-degenerate (collinear) triangle: fp round-off < 0
    area_sq = 0;
  }
  return sqrt(area_sq);
} //}}}
// trajectory-averaged surface coordinate in a bin (assumes values > 0) //{{{
static double AvgSurf(const ArrNDd *sum_surf, const ArrNDi *values,
                      int i, int j, int aa) {
  return GetArr3D(sum_surf, i, j, aa) / GetArr3D(values, i, j, aa);
} //}}}
// calculate areas of two triangles (if appropriate points are defined) //{{{
void calc_4points(double A[3], double B[3], double C[3], double D[3],
                  double *sum_area, int *triangles) {
  if (A[0] != -1 && B[0] != -1 && D[0] != -1) {
    *sum_area += calc_area(A, B, D);
    (*triangles)++;
  }
  if (A[0] != -1 && C[0] != -1 && D[0] != -1) {
    *sum_area += calc_area(A, C, D);
    (*triangles)++;
  }
} //}}}

// SurfacePoint() //{{{
/*
 * For each grid point (i.e, each probe) find a bead that's close enough
 * and has maximum/minimum (i.e., top/bottom for a bilayer or vice versa
 * for a brush) coordinate normal to the surface.
 * For each bead:
 * 1) find each in-surface circular bin the bead is in (i.e., index for
 *    which the bead is in-surface-plane at most probe radius + bead
 *    radius away from the grid point); may be multiple grid points for
 *    each bead
 * 2) find if the 3D distance is within the distance probe radius + bead
 *    radius
 * Find a bead with maximum/minimum normal coordinate for each bin. So,
 * it's a loop over all considered beads and nested within is a loop over
 * the grid of probes...
 */
void AddPoint(ArrNDd *surf_step, ArrNDb *bin_use, ArrNDi *surf_bead_ids,
              int i, int j, int k, int id, double coor) {
  SetArr3D(surf_step, i, j, k, coor);
  SetArr3D(bin_use, i, j, k, true);
  SetArr3D(surf_bead_ids, i, j, k, id);
}
void SurfacePoint(SYSTEM System, int id, const int map[2], int axis,
                  double width, OPT opt, const int *bins_step, ArrNDb *bin_use,
                  ArrNDd *surf_step, ArrNDi *surf_bead_ids) {
  BEAD *bead = &System.Bead[id];
  double coor[3]; // coor[0] & [1] are in the surface plane
  coor[0] = bead->Position.v[map[0]];
  coor[1] = bead->Position.v[map[1]];
  coor[2] = bead->Position.v[axis];
  // maximum 3D distance between the probe and the bead...
  double dist = System.BeadType[bead->Type].Radius + opt.probe;
  double max_dist = Square(dist); // ...well, square of
  // minimum and maximum possible grid point for specified in-surface coordinate
  // (note it must use the distance itself, not its square)
  int min[2], max[2];
  for (int aa = 0; aa < 2; aa++) {
    min[aa] = (coor[aa] - dist) / width - 1;
    max[aa] = (coor[aa] + dist) / width + 1;
  }
  // go over all those grid ponts
  int tmp[2];
  for (tmp[0] = min[0]; tmp[0] <= max[0]; tmp[0]++) {
    for (tmp[1] = min[1]; tmp[1] <= max[1]; tmp[1]++) {
      // account for pbc in the grid point: if its coordinate is too high/low,
      // use one from the box's other side
      int grid[2];
      for (int aa = 0; aa < 2; aa++) {
        grid[aa] = tmp[aa];
        if (tmp[aa] >= bins_step[aa]) {
          grid[aa] -= bins_step[aa];
        } else if (tmp[aa] < 0) {
          grid[aa] += bins_step[aa];
        }
      }
      // 1) bead's distance from the grid point in the surface plane
      double d[2]; // the two axes distance, accounting for pbc
      for (int aa = 0; aa < 2; aa++) {
        d[aa] = fabs(coor[aa] - grid[aa] * width);
        // account for pbc
        while (fabs(d[aa]) > (System.Box.Length.v[map[aa]] / 2)) {
          d[aa] -= System.Box.Length.v[map[aa]];
        }
      }
      // actual in-surface-plane disance (well, square of)
      d[0] = Square(d[0]) + Square(d[1]);
      // 2) only use beads close enough to the probe (in surface plane)
      if (d[0] <= max_dist) {
        // 'bottom' surface for bilayers or 'top' surface for brushes
        // -sqrt because we need lower intersection of line and sphere
        double axis_coor = coor[2] - sqrt(max_dist - d[0]);
        double old = GetArr3D(surf_step, grid[0], grid[1], 0);
        if ((opt.in && axis_coor <= old) || (!opt.in && axis_coor >= old)) {
          AddPoint(surf_step, bin_use, surf_bead_ids,
                   grid[0], grid[1], 0, id, axis_coor);
        }
        // 'top' surface for bilayers or 'bottom' surface for brushes
        // +sqrt because we need upper intersection of line and sphere
        axis_coor = coor[2] + sqrt(max_dist - d[0]);
        old = GetArr3D(surf_step, grid[0], grid[1], 1);
        if ((opt.in && axis_coor >= old) || (!opt.in && axis_coor <= old)) {
          AddPoint(surf_step, bin_use, surf_bead_ids,
                   grid[0], grid[1], 1, id, axis_coor);
        }
      }
    }
  }
} //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}
  // <width> - distance between probes //{{{
  double width = 0;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <surf.txt> - output file with averaged surface coordinates
  char file_surf[LINE] = "";
  s_strcpy(file_surf, argv[++count], LINE);
  // <axis> - x, y, or z //{{{
  int axis, // which axis? 0=x, 1=y, 2=z
      map[2]; // map the remaining two axes based on the 'axis' variable
  if (argv[++count][0] == 'x' ) {
    axis = 0;
    map[0] = 1;
    map[1] = 2;
  } else if (argv[count][0] == 'y') {
    axis = 1;
    map[0] = 0;
    map[1] = 2;
  } else if (argv[count][0] == 'z') {
    axis = 2;
    map[0] = 0;
    map[1] = 1;
  } else {
    err_msg("must be 'x', 'y', or 'z'");
    PrintErrorOption("<axis>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  //}}}
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // -wd option
  double distr_width = 0;
  int vals[2];
  FileNumbersOption(argc, argv, 0, 1, "-wd", &distr_width,
                    vals, opt.width_distr, 'd');
  FileOption(argc, argv, "-w", opt.width_avg);
  // --in option
  opt.in = BoolOption(argc, argv, "--in");
  // --bonded option
  opt.bonded = BoolOption(argc, argv, "--bonded");
  // -r option (probe size)
  if (!OneNumberOption(argc, argv, "-r", &opt.probe, 'd')) {
    opt.probe = 0.5;
  }
  // -a option
  FileOption(argc, argv, "-a", opt.area_file);
  // -b option
  opt.bead_file = InitFile;
  if (FileOption(argc, argv, "-b", opt.bead_file.name)) {
    opt.bead_file.type = CoordinateFileType(opt.bead_file.name);
  }
  //}}}

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // -bt option //{{{
  opt.bt_number = 0;
  opt.bt = calloc(Count->BeadType, sizeof *opt.bt);
  bool *flag = calloc(Count->BeadType, sizeof *flag);
  if (TypeOption(argc, argv, "-bt", 'b', true, flag, System)) {
    for (int i = 0; i < Count->BeadType; i++) {
      if (flag[i]) {
        opt.bt[opt.bt_number] = i;
        opt.bt_number++;
      }
    }
    if (opt.bonded) {
      err_msg("when both are used, --bonded takes precedence");
      PrintWarnOption("--bonded/-bt");
    }
  }
  free(flag); //}}}

  // specify radius for beads that have none //{{{
  bool warn = false;
  if (opt.bonded) { // use all beads in moleculs (--bonded)
    for (int i = 0; i < Count->Bonded; i++) {
      int btype = System.Bead[System.Bonded[i]].Type;
      if (System.BeadType[btype].Radius == RADIUS) {
        System.BeadType[btype].Radius = 0.5;
        if (!warn) {
          warn = true;
          err_msg("unspecified bead radius (using 0.5)");
          PrintWarning();
        }
      }
    }
  } else if (opt.bt_number > 0) { // use specified bead types (-bt option)
    for (int i = 0; i < opt.bt_number; i++) {
      BEADTYPE *btype = &System.BeadType[opt.bt[i]];
      if (btype->Radius == RADIUS) {
        btype->Radius = 0.5;
        if (!warn) {
          warn = true;
          err_msg("unspecified bead radius (using 0.5)");
          PrintWarning();
        }
      }
    }
  } else { // use all beads (no option specified)
    for (int i = 0; i < Count->Bead; i++) {
      int btype = System.Bead[i].Type;
      if (System.BeadType[btype].Radius == RADIUS) {
        System.BeadType[btype].Radius = 0.5;
        if (!warn) {
          warn = true;
          err_msg("unspecified bead radius (using 0.5)");
          PrintWarning();
        }
      }
    }
  } //}}}

  // read the first timestep from a coordinate file to get box dimensions //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int line_count = 0;
  if (!ReadTimestep(in, fr, &System, &line_count)) {
    exit(1);
  }
  fclose(fr);
  if (System.Box.Volume == -1) {
    err_msg("missing box dimensions");
    PrintError();
    exit(1);
  }
  double sidelength[3];
  sidelength[0] = System.Box.Length.v[map[0]];
  sidelength[1] = System.Box.Length.v[map[1]];
  sidelength[2] = System.Box.Length.v[axis];
  //}}}

  // number of grid points (i.e., of bins), guarding against box enlargement
  int bin_alloc[2];
  bin_alloc[0] = sidelength[0] / width * 10;
  bin_alloc[1] = sidelength[1] / width * 10;

  if (commons.verbose) {
    VerboseOutput(System);
  }

  /*
   * Number of values in each bin each surface: should be equal to the number of
   * steps, but if there's no bead that falls into the given bin (i.e., <width>
   * is too low), than it may be lower
   */
  ArrNDi *values = CreateArr3Di(bin_alloc[0], bin_alloc[1], 2);
  /*
   * Sum of points for each surface (i.e., the 'proper' coordinates in <axis>
   * direction) in each bin
   *
   * sum_surf/values gives average coordinate for the two surfaces
   */
  ArrNDd *sum_surf = CreateArr3Dd(bin_alloc[0], bin_alloc[1], 2);
  if (!values || !sum_surf) {
    ErrorAlloc("values/sum_surf");
  }
  // distribution of widths (i.e., top-bottom surface distances)
  long int *distr = nullptr;
  int distr_bins = sidelength[2] / distr_width * 10;
  double avg_thickness = 0;
  if (distr_width > 0) {
    distr = calloc(distr_bins, sizeof *distr);
    if (!distr) {
      ErrorAlloc("distr");
    }
  }

  // open input coordinate file
  fr = OpenFile(in.coor.name, "r");

  // write initial stuff to the per-timestep area //{{{
  if (opt.area_file[0] != '\0') {
    FILE *out = PrintBylineOpenFile(opt.area_file, argc, argv);
    count = 1;
    fprintf(out, "# (%d) timestep", count++);
    fprintf(out, "; (%d) surface 1", count++);
    fprintf(out, "; (%d) surface 2", count++);
    fprintf(out, "; (%d) middle surface", count++);
    putc('\n', out);
  } //}}}

  if (opt.width_avg[0] != '\0') {
    FILE *fout = PrintBylineOpenFile(opt.width_avg, argc, argv);
    fprintf(fout, "# (1) step; (2) thickness\n");
    fclose(fout);
  }

  // array for writing surface beads (if -b option is used) //{{{
  bool *write = nullptr;
  if (opt.bead_file.name[0] != '\0') {
    InitOutputCoorFile(opt.bead_file, System, argc, argv);
    write = calloc(Count->Bead, sizeof *write);
  } //}}}

  // main loop //{{{
  int count_coor = 0, // count calculated timesteps
      count_used = 0; // count timesteps from the beginning
  line_count = 0; // count lines in the coor file
  bool warn_box_change = false;
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // decide whether to use this timestep (based on -st/-sk/-e) //{{{
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    } //}}}
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);

      // warn once if box size changed //{{{
      if (!warn_box_change &&
          (fabs(sidelength[0] - System.Box.Length.v[map[0]]) > 0.00001 ||
           fabs(sidelength[1] - System.Box.Length.v[map[1]]) > 0.00001 ||
           fabs(sidelength[2] - System.Box.Length.v[axis]) > 0.00001)) {
        err_msg("box size changed; only coordinates inside the "
                "original box are used for surface averaging");
        PrintWarning();
        warn_box_change = true;
      } //}}}
      sidelength[0] = System.Box.Length.v[map[0]];
      sidelength[1] = System.Box.Length.v[map[1]];
      sidelength[2] = System.Box.Length.v[axis];

      // per-timestp grid size //{{{
      /*
       * Recalculate grid size for this timestep, considering that
       *   1) Box.Length might have changed
       *   2) area calculation requires one more point because for the highest
       *      coordinate, 0th's value for surface is used
       */
      // TODO: check for grid size too large (unlikely, but possible)
      int bins_step[2];
      for (int aa = 0; aa < 2; aa++) {
        bins_step[aa] = sidelength[aa] / width + 1;
      } //}}}

      // allocate memory for temporary arrays //{{{
      // surfaces' coordinates in this step
      ArrNDd *surf_step = CreateArr3Dd(bins_step[0], bins_step[1], 2);
      ArrNDi *surf_bead_ids = CreateArr3Di(bins_step[0], bins_step[1], 2);
      // is a bin used in this step? (akin to the values array)
      ArrNDb *bin_use = CreateArr3Db(bins_step[0], bins_step[1], 2);
      if (!surf_step || !surf_bead_ids || !bin_use) {
        ErrorAlloc("surf_step/surf_bead_ids/bin_use");
      }
      FillArrND(surf_bead_ids, -1);
      /*
       * Seed each surface with the coordinate it searches away from, so that
       * the first bead found always replaces it: without --in, surface 1
       * maximises from the box's centre downwards and surface 2 minimises
       * from the centre upwards; with --in, the two are swapped.
       */
      double seed[2];
      if (opt.in) {
        seed[0] = sidelength[2];
        seed[1] = 0;
      } else {
        seed[0] = 0;
        seed[1] = sidelength[2];
      }
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          SetArr3D(surf_step, i, j, 0, seed[0]);
          SetArr3D(surf_step, i, j, 1, seed[1]);
        }
      } //}}}

      // calculate surface //{{{
      if (opt.bonded) { // use all beads in moleculs (--bonded)
        for (int i = 0; i < Count->BondedCoor; i++) {
          int id = System.BondedCoor[i];
          SurfacePoint(System, id, map, axis, width, opt,
                       bins_step, bin_use, surf_step, surf_bead_ids);
        }
      } else if (opt.bt_number > 0) { // use specified bead types (-bt option)
        for (int i = 0; i < opt.bt_number; i++) {
          BEADTYPE *btype = &System.BeadType[opt.bt[i]];
          for (int j = 0; j < btype->Number; j++) {
            int id = btype->Index[j];
            if (System.Bead[id].InTimestep) {
              SurfacePoint(System, id, map, axis, width, opt,
                           bins_step, bin_use, surf_step, surf_bead_ids);
            }
          }
        }
      } else { // use all beads (no option specified)
        for (int i = 0; i < Count->BeadCoor; i++) {
          int id = System.BeadCoor[i];
          SurfacePoint(System, id, map, axis, width, opt,
                       bins_step, bin_use, surf_step, surf_bead_ids);
        }
      }
      //}}}

      // add to sums //{{{
      double avg_thickness_step = 0;
      int avg_thickness_count = 0;
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          for (int aa = 0; aa < 2; aa++) {
            if (GetArr3D(bin_use, i, j, aa)) {
              AddArr3D(sum_surf, i, j, aa, GetArr3D(surf_step, i, j, aa));
              AddArr3D(values, i, j, aa, 1);
            }
          }
          if ((distr_width > 0 || opt.width_avg[0] != '\0') &&
              GetArr3D(bin_use, i, j, 0) && GetArr3D(bin_use, i, j, 1)) {
            double w = fabs(GetArr3D(surf_step, i, j, 0) -
                            GetArr3D(surf_step, i, j, 1));
            if (distr_width > 0) { // distr is allocated only for -wd
              int bin = w / distr_width;
              distr[bin]++;
            }
            avg_thickness += w;
            avg_thickness_step += w;
            avg_thickness_count++;
          }
        }
      } //}}}
      if (opt.width_avg[0] != '\0') {
        FILE *fout = OpenFile(opt.width_avg, "a");
        fprintf(fout, "%5d ", count_coor);
        if (avg_thickness_count > 0) {
          fprintf(fout, " %lf\n", avg_thickness_step / avg_thickness_count);
        } else {
          fprintf(fout, " %lf\n", 0.0);
        }
        fclose(fout);
      }

      // calculate total area as a sum of areas of triangles //{{{
      if (opt.area_file[0] != '\0') {
        // close the surface periodically: the last grid point in each
        // direction takes the 0th point's value (and its 'is it used' flag)
        for (int i = 0; i < bins_step[0]; i++) {
          for (int aa = 0; aa < 2; aa++) {
            SetArr3D(surf_step, i, bins_step[1]-1, aa,
                     GetArr3D(surf_step, i, 0, aa));
            SetArr3D(bin_use, i, bins_step[1]-1, aa,
                     GetArr3D(bin_use, i, 0, aa));
          }
        }
        for (int j = 0; j < bins_step[1]; j++) {
          for (int aa = 0; aa < 2; aa++) {
            SetArr3D(surf_step, bins_step[0]-1, j, aa,
                     GetArr3D(surf_step, 0, j, aa));
            SetArr3D(bin_use, bins_step[0]-1, j, aa,
                     GetArr3D(bin_use, 0, j, aa));
          }
        }
        double area[3] = {0, 0, 0};
        int triangles[3] = {0, 0, 0}; // number of valid triangles per area
        for (int i = 0; i < (bins_step[0] - 1); i++) {
          for (int j = 0; j < (bins_step[1] - 1); j++) {
            // four points defining the two triangles
            double A[3]; //= {-1, -1, -1};
            double B[3]; //= {-1, -1, -1};
            double C[3]; //= {-1, -1, -1};
            double D[3]; //= {-1, -1, -1};
            // top and bottom surfaces //{{{
            for (int aa = 0; aa < 2; aa++) {
              B[0] = -1;
              C[0] = -1;
              if (GetArr3D(bin_use, i, j, aa) &&
                  GetArr3D(bin_use, i + 1, j + 1, aa)) {
                A[0] = 0;
                A[1] = 0;
                A[2] = GetArr3D(surf_step, i, j, aa);
                D[0] = width;
                D[1] = width;
                D[2] = GetArr3D(surf_step, i + 1, j + 1, aa);
                if (GetArr3D(bin_use, i+1, j, aa)) {
                  B[0] = width;
                  B[1] = 0;
                  B[2] = GetArr3D(surf_step, i + 1, j, aa);
                }
                if (GetArr3D(bin_use, i, j + 1, aa)) {
                  C[0] = 0;
                  C[1] = width;
                  C[2] = GetArr3D(surf_step, i, j + 1, aa);
                }
                calc_4points(A, B, C, D, &area[aa], &triangles[aa]);
              }
            } //}}}
            // 'middle' surface //{{{
            B[0] = -1;
            C[0] = -1;
            if (GetArr3D(bin_use, i, j, 0) && GetArr3D(bin_use, i, j, 1) &&
                GetArr3D(bin_use, i + 1, j + 1, 0) &&
                GetArr3D(bin_use, i + 1, j + 1, 1)) {
              A[0] = 0;
              A[1] = 0;
              A[2] = (GetArr3D(surf_step, i, j, 0) +
                      GetArr3D(surf_step, i, j, 1)) / 2;
              D[0] = width;
              D[1] = width;
              D[2] = (GetArr3D(surf_step, i + 1, j + 1, 0) +
                      GetArr3D(surf_step, i + 1, j + 1, 1)) / 2;

              if (GetArr3D(bin_use, i + 1, j, 0) &&
                  GetArr3D(bin_use, i + 1, j, 1)) {
                B[0] = width;
                B[1] = 0;
                B[2] = (GetArr3D(surf_step, i + 1, j, 0) +
                        GetArr3D(surf_step, i + 1, j, 1)) / 2;
              }
              if (GetArr3D(bin_use, i, j + 1, 0) &&
                  GetArr3D(bin_use, i, j + 1, 1)) {
                C[0] = 0;
                C[1] = width;
                C[2] = (GetArr3D(surf_step, i, j + 1, 0) +
                        GetArr3D(surf_step, i, j + 1, 1)) / 2;
              }
              calc_4points(A, B, C, D, &area[2], &triangles[2]);
            }
            //}}}
          }
        }
        // add average triangle areas to the total area if not enough triangles
        int n_triangles = (bins_step[0] - 1) * (bins_step[1] - 1) * 2;
        for (int dd = 0; dd < 3; dd++) {
          if (triangles[dd] > 0) {
            double avg_triangle = area[dd] / triangles[dd];
            area[dd] += avg_triangle * (n_triangles - triangles[dd]);
          }
        }
        double Length_area = sidelength[0] * sidelength[1];
        double width_area = (bins_step[0] - 1) * (bins_step[1] - 1) *
                            Square(width);
        FILE *out = OpenFile(opt.area_file, "a");
        fprintf(out, "%d %lf %lf %lf\n", count_coor,
                                         area[0] * Length_area / width_area,
                                         area[1] * Length_area / width_area,
                                         area[2] * Length_area / width_area);
        fclose(out);
      }
      //}}}

      // save interfacial (surface) beads to a coordinate file //{{{
      if (opt.bead_file.name[0] != '\0') {
        // find which beads were assigned as surface
        InitBoolArray(write, Count->Bead, false);
        // bounds come from the array itself, as it's allocated per timestep
        // (i.e., bins_step, not bin_alloc, points in each direction)
        for (size_t i = 0; i < surf_bead_ids->shape[0]; i++) {
          for (size_t j = 0; j < surf_bead_ids->shape[1]; j++) {
            for (size_t aa = 0; aa < surf_bead_ids->shape[2]; aa++) {
              int id = GetArr3D(surf_bead_ids, i, j, aa);
              if (id > -1) {
                write[id] = true;
              }
            }
          }
        }
        // give [0,0,0] coordinates to beads in the timestep but not on the
        // surface so they don't interfere in vmd visualization
        for (int i = 0; i < Count->BeadCoor; i++) {
          int id = System.BeadCoor[i];
          if (!write[id]) {
            write[id] = true;
            for (int dd = 0; dd < 3; dd++) {
              System.Bead[id].Position.v[dd] = 0;
            }
          }
        }
        WriteTimestep(opt.bead_file, System, count_coor, write, argc, argv);
      } //}}}

      FreeArrND(surf_step);
      FreeArrND(surf_bead_ids);
      FreeArrND(bin_use); //}}}
    } else { //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, commons.silent); //}}}

  // find highest grid point with non-zero surface values //{{{
  /*
   * The two directions must be searched independently - taking both from the
   * first (i,j) found truncates the grid whenever that row happens to be
   * shorter than the others. Also, use the number of values, not their sum,
   * as a surface coordinate can legitimately be 0.
   */
  int max[2] = {0, 0};
  for (size_t i = 0; i < values->shape[0]; i++) {
    for (size_t j = 0; j < values->shape[1]; j++) {
      if (GetArr3D(values, i, j, 0) > 0 || GetArr3D(values, i, j, 1) > 0) {
        if ((int)(i + 1) > max[0]) { // add 1 to go from 0 to max-1
          max[0] = i + 1;
        }
        if ((int)(j + 1) > max[1]) {
          max[1] = j + 1;
        }
      }
    }
  } //}}}

  // write surface to output file //{{{
  // print legend
  FILE *out = PrintBylineOpenFile(file_surf, argc, argv);
  char a[3] = {'x', 'y', 'z'};
  fprintf(out, "# (1) %c coordinate; (2) %c coordinate;", a[map[0]], a[map[1]]);
  fprintf(out, " (3) surface 1; (4) surface 2; (5) average surface\n");
  for (int i = 0; i < max[0]; i++) {
    for (int j = 0; j < max[1]; j++) {
      double surface[3];
      surface[0] = GetArr3D(sum_surf, i, j, 0) / GetArr3D(values, i, j, 0);
      surface[1] = GetArr3D(sum_surf, i, j, 1) / GetArr3D(values, i, j, 1);
      surface[2] = (surface[0] + surface[1]) / 2;

      fprintf(out, "%10.4f %10.4f %10.4f %10.4f %10.4f\n", i * width, j * width,
              surface[0], surface[1], surface[2]);
    }
    putc('\n', out);
  }
  fclose(out); //}}}

  // calculate total area as a sum of areas of triangles //{{{
  if (opt.area_file[0] != '\0' && max[0] > 1 && max[1] > 1) {
    for (int i = 0; i < max[0]; i++) {
      for (int aa = 0; aa < 2; aa++) {
        SetArr3D(values, i, max[1]-1, aa, GetArr3D(values, i, 0, aa));
        SetArr3D(sum_surf, i, max[1]-1, aa, GetArr3D(sum_surf, i, 0, aa));
      }
    }
    for (int j = 0; j < max[1]; j++) {
      for (int aa = 0; aa < 2; aa++) {
        SetArr3D(values, max[0]-1, j, aa, GetArr3D(values, 0, j, aa));
        SetArr3D(sum_surf, max[0]-1, j, aa, GetArr3D(sum_surf, 0, j, aa));
      }
    }
    double area[3] = {0, 0, 0};
    int triangles[3] = {0, 0, 0};
    for (int i = 0; i < (max[0] - 1); i++) {
      for (int j = 0; j < (max[1] - 1); j++) {
        // top and bottom surfaces //{{{
        double A[3]; //= {-1};
        double B[3]; //= {-1};
        double C[3]; //= {-1};
        double D[3]; //= {-1};
        for (int aa = 0; aa < 2; aa++) {
          B[0] = -1;
          C[0] = -1;
          if (GetArr3D(values, i, j, aa) > 0 &&
              GetArr3D(values, i + 1, j + 1, aa) > 0) {
            A[0] = 0;
            A[1] = 0;
            A[2] = AvgSurf(sum_surf, values, i, j, aa);
            D[0] = width;
            D[1] = width;
            D[2] = AvgSurf(sum_surf, values, i + 1, j + 1, aa);
            if (GetArr3D(values, i + 1, j, aa) > 0) {
              B[0] = width;
              B[1] = 0;
              B[2] = AvgSurf(sum_surf, values, i + 1, j, aa);
            }
            if (GetArr3D(values, i, j + 1, aa) > 0) {
              C[0] = 0;
              C[1] = width;
              C[2] = AvgSurf(sum_surf, values, i, j + 1, aa);
            }
            calc_4points(A, B, C, D, &area[aa], &triangles[aa]);
          }
        } //}}}
        // 'middle' surface //{{{
        B[0] = -1;
        C[0] = -1;
        if (GetArr3D(values, i, j, 0) > 0 && GetArr3D(values, i, j, 1) > 0 &&
            GetArr3D(values, i + 1, j + 1, 0) > 0 &&
            GetArr3D(values, i + 1, j + 1, 1) > 0) {
          A[0] = 0;
          A[1] = 0;
          A[2] = (AvgSurf(sum_surf, values, i, j, 0) +
                  AvgSurf(sum_surf, values, i, j, 1)) / 2;
          D[0] = width;
          D[1] = width;
          D[2] = (AvgSurf(sum_surf, values, i + 1, j + 1, 0) +
                  AvgSurf(sum_surf, values, i + 1, j + 1, 1)) / 2;

          if (GetArr3D(values, i + 1, j, 0) > 0 &&
              GetArr3D(values, i + 1, j, 1) > 0) {
            B[0] = width;
            B[1] = 0;
            B[2] = (AvgSurf(sum_surf, values, i + 1, j, 0) +
                    AvgSurf(sum_surf, values, i + 1, j, 1)) / 2;
          }
          if (GetArr3D(values, i, j + 1, 0) > 0 &&
              GetArr3D(values, i, j + 1, 1) > 0) {
            C[0] = 0;
            C[1] = width;
            C[2] = (AvgSurf(sum_surf, values, i, j + 1, 0) +
                    AvgSurf(sum_surf, values, i, j + 1, 1)) / 2;
          }
          calc_4points(A, B, C, D, &area[2], &triangles[2]);
        }
        //}}}
      }
    }

  // add average triangle areas to the total area if not enough triangles
    int n_triangles = (max[0] - 1) * (max[1] - 1) * 2;
    for (int dd = 0; dd < 3; dd++) {
      if (triangles[dd] > 0) {
        double avg_triangle = area[dd] / triangles[dd];
        area[dd] += avg_triangle * (n_triangles - triangles[dd]);
      }
    }
    double Length_area = System.Box.Length.v[map[0]] *
                         System.Box.Length.v[map[1]];
    double width_area = (max[0] - 1) * (max[1] - 1) * Square(width);
    FILE *out = OpenFile(opt.area_file, "a");
    fprintf(out, "# average: (1) surface 1");
    fprintf(out, "; (2) surface 2");
    fprintf(out, "; (3) middle surface\n");
    fprintf(out, "# %lf %lf %lf\n", area[0] * Length_area / width_area,
                                    area[1] * Length_area / width_area,
                                    area[2] * Length_area / width_area);
    fclose(out);
  } //}}}

  // write distribution of widths (-wd option) //{{{
  if (distr_width > 0) {
    long int norm = 0; // normalization factor
    int min = 0, max = distr_bins; // lowest/highest bin
    for (int i = 0; i < distr_bins; i++) {
      norm += distr[i];
      if (distr[i] > 0 && min == 0) {
        min = i;
      }
      if (distr[distr_bins-1-i] > 0 && max == distr_bins) {
        max = distr_bins - i;
      }
    }
    // write data to the file
    out = PrintBylineOpenFile(opt.width_distr, argc, argv);
    fprintf(out, "# (1) distance; (2) distribution\n");
    // for (int i = min; i < max; i++) {
    for (int i = 0; i < distr_bins; i++) {
      double dist = distr_width * (2 * i + 1) / 2;
      fprintf(out, "%10.5f", dist);
      if (distr[i] > 0) {
        fprintf(out, " %lf", (double)(distr[i]) / norm);
      } else {
        fprintf(out, " %lf", 0.0);
      }
      putc('\n', out);
    }
    fprintf(out, "# average thickness: %lf\n", avg_thickness / norm);
    fclose(out);
    // write avg thickness to per-timestep thickness
    if (opt.width_avg[0] != '\0') {
      out = OpenFile(opt.width_avg, "a");
      fprintf(out, "# average thickness: %lf\n", avg_thickness / norm);
      fclose(out);
    }
  } //}}}

  // free arrays //{{{
  FreeSystem(&System);
  FreeArrND(sum_surf);
  FreeArrND(values);
  if (distr_width > 0) {
    free(distr);
  }
  free(opt.bt);
  if (opt.bead_file.name[0] != '\0') {
    free(write);
  }
  //}}}

  return 0;
}
