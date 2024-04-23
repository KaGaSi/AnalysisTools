#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
TODO: CHECK & REWRITE\n\n\
Surface utility determines the first bead (either from the box's centre or \
edges) in each square prism defined by the given <width> parameter, thus \
defining a surface of, e.g., polymer brush or lipid bilayer. The <width> \
slices the box into square prisms along the chosen axis (i.e., if z is the \
chosen axis, the xy plane is chopped into squares, creating \
<width>*<width>*<box length in z> prisms). In each such prism, two beads \
are found corresponding to the two surfaces (e.g., polymer brush on both box \
edges or the two surfaces of a lipid bilayer inside the box).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <surf.txt> <area.txt> "
          "<axis> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<surf.txt>          average surface\n");
  fprintf(ptr, "<area.txt>          per-timestep areas\n");
  fprintf(ptr, "<axis>              calculate along x, y, or z axis\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --in              start from the box's edges "
          "instead of the centre\n");
  fprintf(ptr, "  --bonded          use only beads in molecules\n");
  fprintf(ptr, "  -bt <name(s)>     bead type(s) to use\n");
  // fprintf(ptr, "  -m <mol(s)>       molecule type(s) to use\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int bt_number, *bt; // -bt (number of types; list of the types)
  bool in,            // --in
       bonded;        // --bonded
  FILE_TYPE fout;     // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

/*
 * Number of bins allocated in the two dimensions; it's 10 times larger than the
 * first box dimensions to guard against box size changes (enlargement).
 */
int bin_alloc[2];

// TODO: overall averages to <surf.txt>?

// calculate index in an 1D array that simulates a 3D one /{{{
int id3D(int x, int y, int z) {
  return (z * bin_alloc[0] * bin_alloc[1]) + (y * bin_alloc[0]) + x;
} //}}}

// calculate area of a triangle given three points (Heron's formula) //{{{
double calc_area(double A[3], double B[3], double C[3]) {
  // triangle's sides vectors
  double AB[3], AC[3], BC[3];
  for (int dd = 0; dd < 3; dd++) {
    AB[dd] = B[dd] - A[dd];
    AC[dd] = C[dd] - A[dd];
    BC[dd] = C[dd] - B[dd];
  }
  double a = 0, b = 0, c = 0; // triangle's sidelengths
  for (int dd = 0; dd < 3; dd++) {
    a = VectorLength(BC);
    b = VectorLength(AC);
    c = VectorLength(AB);
  }
  double s = (a + b + c) / 2;
  return sqrt(s * (s - a) * (s - b) * (s - c));
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

// find maximum bin in the two dimensions //{{{
void MaxBin(int bin[2], int max_bin[2]) {
  for (int aa = 0; aa < 2; aa++) {
    if (bin[aa] > max_bin[aa]) {
      max_bin[aa] = bin[aa];
    }
  }
} //}}}
// assign bin to a coordate for going from the middle of the box //{{{
void AddSurfacePointOut(double coor, bool *bin_use, double *surf_step,
                        double sidelength[3], int bin[2]) {
  if (coor <= (sidelength[2] / 2)) {
    if (coor >= surf_step[id3D(bin[0], bin[1], 0)]) {
      surf_step[id3D(bin[0], bin[1], 0)] = coor;
      bin_use[id3D(bin[0], bin[1], 0)] = true;
    }
  } else {
    if (coor <= surf_step[id3D(bin[0], bin[1], 1)]) {
      surf_step[id3D(bin[0], bin[1], 1)] = coor;
      bin_use[id3D(bin[0], bin[1], 1)] = true;
    }
  }
} //}}}
// assign bin to a coordate for going from the sides of the box //{{{
void AddSurfacePointIn(double coor, bool *bin_use, double *surf_step,
                        double sidelength[3], int bin[2]) {
  if (coor <= surf_step[id3D(bin[0], bin[1], 0)]) {
    surf_step[id3D(bin[0], bin[1], 0)] = coor;
    bin_use[id3D(bin[0], bin[1], 0)] = true;
  }
  if (coor >= surf_step[id3D(bin[0], bin[1], 1)]) {
    surf_step[id3D(bin[0], bin[1], 1)] = coor;
    bin_use[id3D(bin[0], bin[1], 1)] = true;
  }
} //}}}
// based on mode, pick a proper bin-assigning function //{{{
void AddSurfacePoint(double coor, bool *bin_use, double *surf_step,
                     double sidelength[3], int bin[2], bool in) {
  if (!in) {
    AddSurfacePointOut(coor, bin_use, surf_step, sidelength, bin);
  } else {
    AddSurfacePointIn(coor, bin_use, surf_step, sidelength, bin);
  }
} //}}}
/*
 * function to assign bins to a single coordinate (including of possible
 * periodic images)
 */ //{{{
void SurfacePoint(SYSTEM System, int id, int map[2], int axis, bool in,
                  double sidelength[3], double width, double *surf_step, bool
                  *bin_use, double lo[2], double hi[2], int max_bin[2]) {
  BEAD *b = &System.Bead[id];
  double coor[3]; // coor[0] & [1] are in the surface plane
  coor[0] = b->Position[map[0]];
  coor[1] = b->Position[map[1]];
  coor[2] = b->Position[axis];
  // fnd bead's bins & test for highest bins ever
  double pbc[2];
  int bin[2];
  for (int aa = 0; aa < 2; aa++) {
    bin[aa] = coor[aa] / width;
    pbc[aa] = coor[aa] + sidelength[aa];
  }
  MaxBin(bin, max_bin);
  AddSurfacePoint(coor[2], bin_use, surf_step, sidelength, bin, in);
  /*
   * use a periodic image if the highest bin is narrower than <width>, so the
   * highest bins don't have too few beads in them
   */
  // in one dimension of the surface plane
  if (pbc[0] > lo[0] && pbc[0] < hi[0]) {
    bin[0] = pbc[0] / width;
    bin[1] = coor[1] / width;
    MaxBin(bin, max_bin);
    AddSurfacePoint(coor[2], bin_use, surf_step, sidelength, bin, in);
  }
  // in second dimension of the surface plane
  if (pbc[1] > lo[1] && pbc[1] < hi[1]) {
    bin[0] = coor[0] / width;
    bin[1] = pbc[1] / width;
    MaxBin(bin, max_bin);
    AddSurfacePoint(coor[2], bin_use, surf_step, sidelength, bin, in);
  }
  // in both dimensions of the surface plane
  if (pbc[0] > lo[0] && pbc[0] < hi[0] &&
      pbc[1] > lo[1] && pbc[1] < hi[1]) {
    bin[0] = pbc[0] / width;
    bin[1] = pbc[1] / width;
    MaxBin(bin, max_bin);
    AddSurfacePoint(coor[2], bin_use, surf_step, sidelength, bin, in);
  }
} //}}}

// TODO: reinstate molecular condition

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 9, all = common + 3, count = 0,
      req_arg = 5;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "--in");
  strcpy(option[count++], "--bonded");
  // strcpy(option[count++], "-m");
  strcpy(option[count++], "-bt");
  OptionCheck(argc, argv, count, req_arg, common, all, option, true); //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // arguments & options before reading system data //{{{
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  snprintf(in.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <width> - width of a single bin
  double width = 0;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  // <surf.txt> - output file with averaged surface coordinates
  char file_surf[LINE] = "";
  snprintf(file_surf, LINE, "%s", argv[++count]);
  // <area.txt> - output file with average area
  char file_area[LINE] = "";
  snprintf(file_area, LINE, "%s", argv[++count]);
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
    strcpy(ERROR_MSG, "must be 'x', 'y', or 'z'");
    PrintErrorOption("<axis>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  //}}}
  opt->c = CommonOptions(argc, argv, LINE, in);
  opt->in = BoolOption(argc, argv, "--in");
  opt->bonded = BoolOption(argc, argv, "--bonded"); //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, opt->c.detailed);
  COUNT *Count = &System.Count;

  // -bt option
  opt->bt_number = 0;
  opt->bt = calloc(Count->BeadType, sizeof *opt->bt);
  bool *flag = calloc(Count->BeadType, sizeof *flag);
  if (BeadTypeOption(argc, argv, "-bt", true, flag, System)) {
    for (int i = 0; i < Count->BeadType; i++) {
      if (flag[i]) {
        opt->bt[opt->bt_number] = i;
        opt->bt_number++;
      }
    }
  }

  // read the first timestep from a coordinate file to get box dimensions
  FILE *fr = OpenFile(in.coor.name, "r");
  int line_count = 0;
  if (!ReadTimestep(in, fr, &System, &line_count)) {
    exit(1);
  }
  fclose(fr);
  double sidelength[3];
  sidelength[0] = System.Box.Length[map[0]];
  sidelength[1] = System.Box.Length[map[1]];
  sidelength[2] = System.Box.Length[axis];

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // number of bins, guarding against box enlargement
  bin_alloc[0] = sidelength[0] / width * 10;
  bin_alloc[1] = sidelength[1] / width * 10;

  /*
   * Number of values in each bin each surface: should be equal to the number of
   * steps, but if there's no bead that falls into the given bin (i.e., <width>
   * is too low), than it may be lower
   */
  int *values = calloc(bin_alloc[0] * bin_alloc[1] * 2, sizeof *values);
  /*
   * Sum of points for each surface (i.e., the 'proper' coordinates in <axis>
   * direction) in each bin
   *
   * sum_surf/values gives average coordinate for the two surfaces
   */
  double *sum_surf = calloc(bin_alloc[0] * bin_alloc[1] * 2, sizeof *sum_surf);

  // open input coordinate file
  fr = OpenFile(in.coor.name, "r");

  // write initial stuff to the per-timestep area //{{{
  PrintByline(file_area, argc, argv);
  FILE *out = OpenFile(file_area, "a");
  count = 1;
  fprintf(out, "# (%d) timestep", count++);
  fprintf(out, "; (%d) surface 1", count++);
  fprintf(out, "; (%d) surface 2", count++);
  fprintf(out, "; (%d) middle surface", count++);
  putc('\n', out); //}}}

  // main loop //{{{
  int count_coor = 0, // count calculated timesteps
      count_used = 0, // count timesteps from the beginning
      max_bin[2] = {0, 0}; // highest bin ids (i.e., number of bins to write)
  line_count = 0; // count lines in the coor file
  bool warn_box_change = false;
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // decide whether to use this timestep (based on -st/-sk/-e) //{{{
    bool use = false;
    if (count_coor >= opt->c.start &&
        (count_coor <= opt->c.end || opt->c.end == -1) &&
        ((count_coor - opt->c.start) % opt->c.skip) == 0) {
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
          (fabs(sidelength[0] - System.Box.Length[map[0]]) > 0.00001 ||
           fabs(sidelength[1] - System.Box.Length[map[1]]) > 0.00001 ||
           fabs(sidelength[2] - System.Box.Length[axis]) > 0.00001)) {
        strcpy(ERROR_MSG, "box size changed; only coordinates inside the "
               "original box are used for surface averaging");
        PrintWarning();
        warn_box_change = true;
      } //}}}
      // number of bins //{{{
      sidelength[0] = System.Box.Length[map[0]];
      sidelength[1] = System.Box.Length[map[1]];
      sidelength[2] = System.Box.Length[axis];

      /*
       * Recalculate number of bins for this timestep, considering that
       *   1) Box.Length might have changed
       *   2) area calculation requires one more because for the highest bin id,
       *      0th's value for surface is used
       */
      int bins_step[2];
      for (int aa = 0; aa < 2; aa++) {
        bins_step[aa] = sidelength[aa] / width + 1;
        // *1000 to ensure width is integer as otherwise it can wreak havoc
        if ((fmod(sidelength[aa] * 1000, width * 1000) / 1000) > 0.001) {
          bins_step[aa]++;
        }
      }
      /*
      * Find lowest and highest coordinates for the topmost bins; these are used
      * to use beads from a periodic image should the highest bin be less wide
      * than 'width'
      */
      double lo[2], hi[2];
      for (int aa = 0; aa < 2; aa++) {
        lo[aa] = (int)(sidelength[aa] / width) * width;
        if (fabs(lo[aa]-sidelength[aa]) > 0.00001) {
          hi[aa] = lo[aa] + width;
        } else {
          hi[aa] = lo[aa];
        }
      } //}}}

      // allocate memory for temporary arrays //{{{
      // surfaces' coordinates in this step
      double *surf_step = calloc(bin_alloc[0] * bin_alloc[1] * 2,
                                 sizeof *surf_step);
      // is a bin used in this step? (akin to the values array)
      bool *bin_use = calloc(bin_alloc[0] * bin_alloc[1] * 2, sizeof *bin_use);
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          if (!opt->in) {
            surf_step[id3D(i, j, 0)] = 0;
            surf_step[id3D(i, j, 1)] = sidelength[2];
          } else {
            surf_step[id3D(i, j, 0)] = sidelength[2];
            surf_step[id3D(i, j, 1)] = 0;
          }
        }
      } //}}}

      // calculate surface //{{{
      if (opt->bonded) { // use all beads in moleculs
        for (int i = 0; i < Count->Bonded; i++) {
          int id = System.Bonded[i];
          SurfacePoint(System, id, map, axis, opt->in, sidelength, width,
                       surf_step, bin_use, lo, hi, max_bin);
        }
      } else if (opt->bt_number > 0) { // use specified bead types
        for (int i = 0; i < opt->bt_number; i++) {
          BEADTYPE *btype = &System.BeadType[opt->bt[i]];
          for (int j = 0; j < btype->Number; j++) {
            int id = btype->Index[j];
            SurfacePoint(System, id, map, axis, opt->in, sidelength, width,
                         surf_step, bin_use, lo, hi, max_bin);
          }
        }
      } else {
        for (int i = 0; i < Count->BeadCoor; i++) {
          int id = System.BeadCoor[i];
          SurfacePoint(System, id, map, axis, opt->in, sidelength, width,
                       surf_step, bin_use, lo, hi, max_bin);
        }
      } //}}}

      // add to sums //{{{
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          if (bin_use[id3D(i, j, 0)]) {
            sum_surf[id3D(i, j, 0)] += surf_step[id3D(i, j, 0)];
            values[id3D(i, j, 0)]++;
          }
          if (bin_use[id3D(i, j, 1)]) {
            sum_surf[id3D(i, j, 1)] += surf_step[id3D(i, j, 1)];
            values[id3D(i, j, 1)]++;
          }
        }
      } //}}}

      // calculate total area as a sum of areas of triangles //{{{
      for (int i = 0; i < bins_step[0]; i++) {
        surf_step[id3D(i, bins_step[1]-1, 0)] = surf_step[id3D(i, 0, 0)];
        surf_step[id3D(i, bins_step[1]-1, 1)] = surf_step[id3D(i, 0, 1)];
        bin_use[id3D(i, bins_step[1]-1, 0)] = bin_use[id3D(i, 0, 0)];
        bin_use[id3D(i, bins_step[1]-1, 1)] = bin_use[id3D(i, 0, 1)];
      }
      for (int j = 0; j < bins_step[1]; j++) {
        surf_step[id3D(bins_step[0]-1, j, 0)] = surf_step[id3D(0, j, 0)];
        surf_step[id3D(bins_step[0]-1, j, 1)] = surf_step[id3D(0, j, 1)];
        bin_use[id3D(bins_step[0]-1, j, 0)] = bin_use[id3D(0, j, 0)];
        bin_use[id3D(bins_step[0]-1, j, 1)] = bin_use[id3D(0, j, 1)];
      }
      double area[3] = {0, 0, 0};
      int triangles[3] = {0, 0, 0}; // number of valid triangles per area
      for (int i = 0; i < (bins_step[0] - 1); i++) {
        for (int j = 0; j < (bins_step[1] - 1); j++) {
          // four points defining the two triangles
          double A[3] = {-1}, B[3] = {-1}, C[3] = {-1}, D[3] = {-1};
          // top and bottom surfaces //{{{
          for (int aa = 0; aa < 2; aa++) {
            if (bin_use[id3D(i, j, aa)] && bin_use[id3D(i+1, j+1, aa)]) {
              A[0] = A[1] = 0;
              A[2] = surf_step[id3D(i, j, aa)];
              D[0] = D[1] = width;
              D[2] = surf_step[id3D(i+1, j+1, aa)];
              if (bin_use[id3D(i+1, j, aa)]) {
                // B[0] = rest[0];
                B[0] = width;
                B[1] = 0;
                B[2] = surf_step[id3D(i+1, j, aa)];
              }
              if (bin_use[id3D(i, j+1, aa)]) {
                C[0] = 0;
                C[1] = width;
                C[2] = surf_step[id3D(i, j+1, aa)];
              }
              calc_4points(A, B, C, D, &area[aa], &triangles[aa]);
            }
          } //}}}
          // 'middle' surface //{{{
          A[0] = B[0] = C[0] = D[0] = -1;
          if (bin_use[id3D(i, j, 0)] && bin_use[id3D(i, j, 1)] &&
              bin_use[id3D(i+1, j+1, 0)] && bin_use[id3D(i+1, j+1, 1)]) {
            A[0] = A[1] = 0;
            A[2] = (surf_step[id3D(i, j, 0)] + surf_step[id3D(i, j, 1)]) / 2;
            D[0] = D[1] = width;
            D[2] = (surf_step[id3D(i+1, j+1, 0)] + surf_step[id3D(i+1, j+1, 1)]) / 2;

            if (bin_use[id3D(i+1, j, 0)] && bin_use[id3D(i+1, j, 1)]) {
              B[0] = width;
              B[1] = 0;
              B[2] = (surf_step[id3D(i+1, j, 0)] + surf_step[id3D(i+1, j, 1)]) / 2;
            }
            if (bin_use[id3D(i, j+1, 0)] && bin_use[id3D(i, j+1, 1)]) {
              // C[0] = rest[1];
              C[0] = width;
              C[1] = 0;
              C[2] = (surf_step[id3D(i, j+1, 0)] + surf_step[id3D(i, j+1, 1)]) / 2;
            }
            calc_4points(A, B, C, D, &area[2], &triangles[2]);
          }
          //}}}
        }
      }
      // add average triangle areas to the total area if not enough triangles
      int n_triangles = (bins_step[0] - 1) * (bins_step[1] - 1) * 2;
      double avg_triangle[3];
      for (int dd = 0; dd < 3; dd++) {
        avg_triangle[dd] = area[dd] / triangles[dd];
        area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
      }
      double Length_area = sidelength[0] * sidelength[1];
      double width_area = (bins_step[0] - 1) * (bins_step[1] - 1) * SQR(width);
      fprintf(out, "%d %lf %lf %lf\n", count_coor,
                                       area[0] * Length_area / width_area,
                                       area[1] * Length_area / width_area,
                                       area[2] * Length_area / width_area);
      //}}}
      free(surf_step);
      free(bin_use); //}}}
    } else { //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  fclose(out);
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // write surface to output file //{{{
  PrintByline(file_surf, argc, argv);
  // print legend
  out = OpenFile(file_surf, "a");
  char a[3] = {'x', 'y', 'z'};
  fprintf(out, "# (1) %c coordinate; (2) %c coordinate;", a[map[0]], a[map[1]]);
  fprintf(out, " (3) surface 1; (4) surface 2; (5) average surface\n");

  // max_bin[] are bin ids, so +1 to get count which includes id=0
  max_bin[0]++;
  max_bin[1]++;
  for (int i = 0; i < max_bin[0]; i++) {
    count = 0;
    for (int j = 0; j < max_bin[1]; j++) {
      if (values[id3D(i, j, 0)] > 0 ||
          values[id3D(i, j, 1)] > 0) {
        count++;
        fprintf(out, "%7.4f", width*(2*i+1)/2);
        fprintf(out, " %7.4f", width*(2*j+1)/2);
        double surface[2];
        for (int aa = 0; aa < 2; aa++) {
          // if (values[i][j][aa] > 0) {
          if (values[id3D(i, j, aa)] > 0) {
            // surface[aa] = surf[i][j][aa] / values[i][j][aa];
            surface[aa] = sum_surf[id3D(i, j, aa)] / values[id3D(i, j, aa)];
            fprintf(out, " %7.4f", surface[aa]);
          } else {
            fprintf(out, "       ?");
          }
        }
        if (values[id3D(i, j, 0)] > 0 ||
            values[id3D(i, j, 1)] > 0) {
          fprintf(out, " %7.4f", (surface[0]+surface[1])/2);
        } else {
          fprintf(out, "       ?");
        }
        fprintf(out, "\n");
      }
    }
    if (count != 0) {
      fprintf(out, "\n");
    }
  }
  fclose(out); //}}}

  // calculate total area as a sum of areas of triangles //{{{
  for (int i = 0; i < max_bin[0]; i++) {
    values[id3D(i, max_bin[1]-1, 0)] = values[id3D(i, 0, 0)];
    values[id3D(i, max_bin[1]-1, 1)] = values[id3D(i, 0, 1)];
    sum_surf[id3D(i, max_bin[1]-1, 0)] = sum_surf[id3D(i, 0, 0)];
    sum_surf[id3D(i, max_bin[1]-1, 1)] = sum_surf[id3D(i, 0, 1)];
  }
  for (int j = 0; j < max_bin[1]; j++) {
    values[id3D(max_bin[0]-1, j, 0)] = values[id3D(0, j, 0)];
    values[id3D(max_bin[0]-1, j, 1)] = values[id3D(0, j, 1)];
    sum_surf[id3D(max_bin[0]-1, j, 0)] = sum_surf[id3D(0, j, 0)];
    sum_surf[id3D(max_bin[0]-1, j, 1)] = sum_surf[id3D(0, j, 1)];
  }
  double area[3] = {0, 0, 0};
  int triangles[3] = {0, 0, 0};
  for (int i = 0; i < (max_bin[0] - 1); i++) {
    for (int j = 0; j < (max_bin[1] - 1); j++) {
      // top and bottom surfaces //{{{
      double A[3] = {-1}, B[3] = {-1}, C[3] = {-1}, D[3] = {-1};
      for (int aa = 0; aa < 2; aa++) {
        if (values[id3D(i, j, aa)] > 0 && values[id3D(i+1, j+1, aa)] > 0) {
          A[0] = A[1] = 0;
          A[2] = sum_surf[id3D(i, j, aa)] / values[id3D(i, j, aa)];
          D[0] = D[1] = width;
          D[2] = sum_surf[id3D(i+1, j+1, aa)] / values[id3D(i+1, j+1, aa)];
          if (values[id3D(i+1, j, aa)] > 0) {
            B[0] = width;
            B[1] = 0;
            B[2] = sum_surf[id3D(i+1, j, aa)] / values[id3D(i+1, j, aa)];
          }
          if (values[id3D(i, j+1, aa)] > 0) {
            C[0] = 0;
            C[1] = width;
            C[2] = sum_surf[id3D(i, j+1, aa)] / values[id3D(i, j+1, aa)];
          }
          calc_4points(A, B, C, D, &area[aa], &triangles[aa]);
        }
      } //}}}
      // 'middle' surface //{{{
      A[0] = B[0] = C[0] = D[0] = -1;
      if (values[id3D(i, j, 0)] > 0 && values[id3D(i, j, 1)] > 0 &&
          values[id3D(i+1, j+1, 0)] > 0 && values[id3D(i+1, j+1, 1)] > 0) {
        A[0] = A[1] = 0;
        A[2] = (sum_surf[id3D(i, j, 0)] / values[id3D(i, j, 0)] +
                sum_surf[id3D(i, j, 1)] / values[id3D(i, j, 1)]) / 2;
        D[0] = D[1] = width;
        D[2] = (sum_surf[id3D(i+1, j+1, 0)] / values[id3D(i+1, j+1, 0)] +
                sum_surf[id3D(i+1, j+1, 1)] / values[id3D(i+1, j+1, 1)]) / 2;

        if (values[id3D(i+1, j, 0)] > 0 &&
            values[id3D(i+1, j, 1)] > 0) {
          B[0] = width;
          B[1] = 0;
          B[2] = (sum_surf[id3D(i+1, j, 0)] / values[id3D(i+1, j, 0)] +
                  sum_surf[id3D(i+1, j, 1)] / values[id3D(i+1, j, 1)]) / 2;
        }
        if (values[id3D(i, j+1, 0)] > 0 &&
            values[id3D(i, j+1, 1)] > 0) {
          C[0] = width;
          C[1] = 0;
          C[2] = (sum_surf[id3D(i, j+1, 0)] / values[id3D(i, j+1, 0)] +
                  sum_surf[id3D(i, j+1, 1)] / values[id3D(i, j+1, 1)]) / 2;
        }
        calc_4points(A, B, C, D, &area[2], &triangles[2]);
      }
      //}}}
    }
  }

  // add average triangle areas to the total area if not enough triangles
  int n_triangles = (max_bin[0] - 1) * (max_bin[1] - 1) * 2;
  double avg_triangle[3];
  for (int dd = 0; dd < 3; dd++) {
    avg_triangle[dd] = area[dd] / triangles[dd];
    area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
  }
  double Length_area = System.Box.Length[0] * System.Box.Length[1];
  double width_area = (max_bin[0] - 1) * (max_bin[1] - 1) * SQR(width);
  out = OpenFile(file_area, "a");
  fprintf(out, "# average: (1) surface 1; (2) surface 2; (3) middle surface\n");
  fprintf(out, "# %lf %lf %lf\n", area[0] * Length_area / width_area,
                                  area[1] * Length_area / width_area,
                                  area[2] * Length_area / width_area);
  fclose(out); //}}}

  FreeSystem(&System);
  free(sum_surf);
  free(values);
  free(opt->bt);
  free(opt);

  return 0;
}
