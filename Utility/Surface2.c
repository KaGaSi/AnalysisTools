#include "../AnalysisTools.h"
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

// Uses Identification of the Truly Interfacial Molecules (ITIM) from
// https://doi.org/10.1002/jcc.20852
// ...assumes probe + bead distance is 1 (i.e., r_c for DPD bead), that is probe
// radius and bead radius are the same, specifically 0.5

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Surface utility determines the first bead (either going from the box's centre \
- defaault behaviour - \
or from its edges - if --in option is used) \
in each square prism defined by the given <width> parameter, \
thus defining a surface of, e.g., polymer brush or lipid bilayer. The <width> \
'slices' the box into square prisms along the chosen axis (i.e., if z is the \
chosen axis, the xy plane is chopped into squares, creating \
<width>*<width>*<box length in z> prisms). In each such prism, two beads \
are found corresponding to the two surfaces (e.g., polymer brush on both box \
edges or the two surfaces of a lipid bilayer inside the box). The surface \
coordinates are saved into the <surf.txt> output file, while their area \
(calculated from triangles defined by the surface coordinates) are written \
into the <area.txt>. What bead types are considered as possible surface \
beads is controlled via --bonded and -bt options.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <surf.txt> "
          "<axis> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<surf.txt>          average surface\n");
  fprintf(ptr, "<axis>              calculate along x, y, or z axis\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --in              start from the box's edges "
          "instead of the centre\n");
  fprintf(ptr, "  --bonded          use only beads in molecules\n");
  fprintf(ptr, "  -bt <name(s)>     bead type(s) to use\n");
  fprintf(ptr, "  -w <file> <width> calculate distribution of widths"
          "with given single bin width\n");
  fprintf(ptr, "  -a <area.txt>     per-timestep areas\n");
  fprintf(ptr, "  -b <file>         save per-timestep surface beads"
          "to a coordinate file\n");
  // fprintf(ptr, "  -m <mol(s)>       molecule type(s) to use\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int bt_number, *bt;    // -bt (number of types; list of the types)
  double distr_width,    // -w (width of distribution bin)
         probe;          // -r (probe radius)
  char distr_file[LINE], // -w (filename)
       area_file[LINE];  // -a (filename)
  bool in,               // --in
       bonded;           // --bonded
  FILE_TYPE bead_file;   // -b (filename)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

/*
 * grid size for probes' coordinates in surface plane; it's 10 times larger than
 * the first box dimensions to guard against box size changes (enlargement).
 */
int bin_alloc[2];

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

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 8, all = common + 7, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "--in");
  strcpy(option[count++], "--bonded");
  // strcpy(option[count++], "-m");
  strcpy(option[count++], "-bt");
  strcpy(option[count++], "-w");
  strcpy(option[count++], "-a");
  strcpy(option[count++], "-b");
  strcpy(option[count++], "-r");
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
  // <width> - distance between probes
  double width = 0;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  // <surf.txt> - output file with averaged surface coordinates
  char file_surf[LINE] = "";
  snprintf(file_surf, LINE, "%s", argv[++count]);
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
  // -w option
  double distr_width = 0;
  int vals[2];
  FileDoubleOption(argc, argv, 1, "-w", &distr_width, vals, opt->distr_file);
  opt->c = CommonOptions(argc, argv, LINE, in);
  opt->in = BoolOption(argc, argv, "--in");
  opt->bonded = BoolOption(argc, argv, "--bonded");
  if (!DoubleOption1(argc, argv, "-r", &opt->probe)) {
    opt->probe = 0.5;
  }
  // -a option
  FileOption(argc, argv, "-a", opt->area_file);
  // -b option
  opt->bead_file = InitFile;
  if (FileOption(argc, argv, "-b", opt->bead_file.name)) {
    opt->bead_file.type = CoordinateFileType(opt->bead_file.name);
  }
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  double max_rad = 0;
  for (int i = 0; i < Count->BeadType; i++) {
    double *r = &System.BeadType[i].Radius;
    if (*r == RADIUS) {
      *r = 0.5;
    }
    if (*r > max_rad) {
      max_rad = *r;
    }
  }

  // -bt option //{{{
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
    if (opt->bonded) {
      strcpy(ERROR_MSG, "when both are used, --bonded takes precedence");
      PrintWarnOption("--bonded/-bt");
    }
  }
  free(flag); //}}}

  // read the first timestep from a coordinate file to get box dimensions //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int line_count = 0;
  if (!ReadTimestep(in, fr, &System, &line_count)) {
    exit(1);
  }
  fclose(fr);
  if (System.Box.Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintError();
    exit(1);
  }
  double sidelength[3];
  sidelength[0] = System.Box.Length[map[0]];
  sidelength[1] = System.Box.Length[map[1]];
  sidelength[2] = System.Box.Length[axis];
  //}}}

  // number of grid points (i.e., of bins), guarding against box enlargement
  bin_alloc[0] = sidelength[0] / width * 10;
  bin_alloc[1] = sidelength[1] / width * 10;

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

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
  // distribution of widths (i.e., distances between top and bottom surface)
  long int *distr = NULL;
  int distr_bins = sidelength[2] / distr_width * 10;
  if (distr_width > 0) {
    distr = calloc(distr_bins, sizeof *distr);
  }

  // open input coordinate file
  fr = OpenFile(in.coor.name, "r");

  // write initial stuff to the per-timestep area //{{{
  if (opt->area_file[0] != '\0') {
    PrintByline(opt->area_file, argc, argv);
    FILE *out = OpenFile(opt->area_file, "a");
    count = 1;
    fprintf(out, "# (%d) timestep", count++);
    fprintf(out, "; (%d) surface 1", count++);
    fprintf(out, "; (%d) surface 2", count++);
    fprintf(out, "; (%d) middle surface", count++);
    putc('\n', out);
  } //}}}

  bool *write = NULL;
  if (opt->bead_file.name[0] != '\0') {
    InitCoorFile(opt->bead_file, System, argc, argv);
    write = calloc(Count->Bead, sizeof *write);
  }

  // main loop //{{{
  int count_coor = 0, // count calculated timesteps
      count_used = 0; // count timesteps from the beginning
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
      sidelength[0] = System.Box.Length[map[0]];
      sidelength[1] = System.Box.Length[map[1]];
      sidelength[2] = System.Box.Length[axis];

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
      double *surf_step = calloc(bin_alloc[0] * bin_alloc[1] * 2,
                                 sizeof *surf_step);
      int *surf_bead_ids = calloc(bin_alloc[0] * bin_alloc[1] * 2,
                                  sizeof *surf_bead_ids);
      InitIntArray(surf_bead_ids, bin_alloc[0] * bin_alloc[1] * 2, -1);
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
      if (opt->bonded) { // use all beads in moleculs (--bonded)
        for (int i = 0; i < Count->BondedCoor; i++) {
          int id = System.BondedCoor[i];
          BEAD *bead = &System.Bead[id];
          double coor[3]; // coor[0] & [1] are in the surface plane
          coor[0] = bead->Position[map[0]];
          coor[1] = bead->Position[map[1]];
          coor[2] = bead->Position[axis];
          // maximum 3D distance between the probe and the bead
          double probe_bead_max_dist =
            SQR(System.BeadType[bead->Type].Radius + opt->probe);
          // go over all the relevant grid points
          int min[2], max[2];
          for (int aa = 0; aa < 2; aa++) {
            min[aa] = (coor[aa] - probe_bead_max_dist) / width - 1;
            max[aa] = (coor[aa] + probe_bead_max_dist) / width + 1;
          }
          int tmp[2];
          for (tmp[0] = min[0]; tmp[0] <= max[0]; tmp[0]++) {
            for (tmp[1] = min[1]; tmp[1] <= max[1]; tmp[1]++) {
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
                while (fabs(d[aa]) > (System.Box.Length[map[aa]] / 2)) {
                  d[aa] -= System.Box.Length[map[aa]];
                }
              }
              // actual in-surface-plane disance
              d[0] = SQR(d[0]) + SQR(d[1]);
              // 2) only use beads at most the maximum distance from probe
              if (d[0] <= probe_bead_max_dist) {
                // +sqrt because for top surface, we need upper intersection of
                // line and sphere
                // 'bottom' surface for bilayers, 'top' for brushes
                int the_bin = id3D(grid[0], grid[1], 0);
                double axis_coor = coor[2] - sqrt(probe_bead_max_dist - d[0]);
                if (axis_coor <= surf_step[the_bin]) {
                  surf_step[the_bin] = axis_coor;
                  bin_use[the_bin] = true;
                  surf_bead_ids[the_bin] = id;
                }
                // 'top' surface for bilayers, 'bottom' for brushes
                the_bin = id3D(grid[0], grid[1], 1);
                axis_coor = coor[2] + sqrt(probe_bead_max_dist - d[0]);
                if (axis_coor >= surf_step[the_bin]) {
                  surf_step[the_bin] = axis_coor;
                  bin_use[the_bin] = true;
                  surf_bead_ids[the_bin] = id;
                }
              }
            }
          }
        }
      } //}}}

      // add to sums //{{{
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          for (int aa = 0; aa < 2; aa++) {
            int id = id3D(i, j, aa);
            if (bin_use[id]) {
              // printf("OK %d %lf %lf\n", aa, sum_surf[id], surf_step[id]);
              sum_surf[id] += surf_step[id];
              values[id]++;
            }
          }
          int id0 = id3D(i, j, 0),
              id1 = id3D(i, j, 1);
          if (distr_width > 0 && surf_step[id0] != -1 && surf_step[id1] != -1) {
            double w = fabs(surf_step[id0] - surf_step[id1]);
            int bin = w / distr_width;
            distr[bin]++;
          }
        }
      } //}}}

      // calculate total area as a sum of areas of triangles //{{{
      if (opt->area_file[0] != '\0') {
        for (int i = 0; i < bins_step[0]; i++) {
          surf_step[id3D(i, bins_step[1]-1, 0)] = surf_step[id3D(i, 0, 0)];
          surf_step[id3D(i, bins_step[1]-1, 1)] = surf_step[id3D(i, 0, 1)];
        }
        for (int j = 0; j < bins_step[1]; j++) {
          surf_step[id3D(bins_step[0]-1, j, 0)] = surf_step[id3D(0, j, 0)];
          surf_step[id3D(bins_step[0]-1, j, 1)] = surf_step[id3D(0, j, 1)];
        }
        double area[3] = {0, 0, 0};
        int triangles[3] = {0, 0, 0}; // number of valid triangles per area
        for (int i = 0; i < (bins_step[0] - 1); i++) {
          for (int j = 0; j < (bins_step[1] - 1); j++) {
            // four points defining the two triangles
            double A[3] = {-1}, B[3] = {-1}, C[3] = {-1}, D[3] = {-1};
            // top and bottom surfaces //{{{
            for (int aa = 0; aa < 2; aa++) {
              if (surf_step[id3D(i, j, aa)] != -1 &&
                  surf_step[id3D(i+1, j+1, aa)] != -1) {
                A[0] = A[1] = 0;
                A[2] = surf_step[id3D(i, j, aa)];
                D[0] = D[1] = width;
                D[2] = surf_step[id3D(i+1, j+1, aa)];
                if (surf_step[id3D(i+1, j, aa)] != -1) {
                  // B[0] = rest[0];
                  B[0] = width;
                  B[1] = 0;
                  B[2] = surf_step[id3D(i+1, j, aa)];
                }
                if (surf_step[id3D(i, j+1, aa)] != -1) {
                  C[0] = 0;
                  C[1] = width;
                  C[2] = surf_step[id3D(i, j+1, aa)];
                }
                calc_4points(A, B, C, D, &area[aa], &triangles[aa]);
              }
            } //}}}
            // 'middle' surface //{{{
            A[0] = B[0] = C[0] = D[0] = -1;
            if (surf_step[id3D(i, j, 0)] != -1 &&
                surf_step[id3D(i, j, 1)] != -1 &&
                surf_step[id3D(i+1, j+1, 0)] != -1 &&
                surf_step[id3D(i+1, j+1, 1)] != -1) {
              A[0] = A[1] = 0;
              A[2] = (surf_step[id3D(i, j, 0)] + surf_step[id3D(i, j, 1)]) / 2;
              D[0] = D[1] = width;
              D[2] = (surf_step[id3D(i+1, j+1, 0)] + surf_step[id3D(i+1, j+1, 1)]) / 2;

              if (surf_step[id3D(i+1, j, 0)] != -1 &&
                  surf_step[id3D(i+1, j, 1)] != -1) {
                B[0] = width;
                B[1] = 0;
                B[2] = (surf_step[id3D(i+1, j, 0)] + surf_step[id3D(i+1, j, 1)]) / 2;
              }
              if (surf_step[id3D(i, j+1, 0)] != -1 &&
                  surf_step[id3D(i, j+1, 1)] != -1) {
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
        FILE *out = OpenFile(opt->area_file, "a");
        fprintf(out, "%d %lf %lf %lf\n", count_coor,
                                         area[0] * Length_area / width_area,
                                         area[1] * Length_area / width_area,
                                         area[2] * Length_area / width_area);
        fclose(out);
      }
      //}}}

      // save interfacial (surface) beads to a coordinate file
      if (opt->bead_file.name[0] != '\0') {
        // find which beads were assigned as surface
        InitBoolArray(write, Count->Bead, false);
        for (int i = 0; i < bin_alloc[0]; i++) {
          for (int j = 0; j < bin_alloc[1]; j++) {
            for (int aa = 0; aa < 2; aa++) {
              int id = surf_bead_ids[id3D(i, j, aa)];
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
              System.Bead[id].Position[dd] = 0;
            }
          }
        }
        WriteTimestep(opt->bead_file, System, count_coor, write, argc, argv);
      }

      free(surf_step);
      free(bin_use);
      free(surf_bead_ids); //}}}
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
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // find highest grid point with non-zero surface values //{{{
  int max[2] = {-1, -1};
  for (int i = (bin_alloc[0] - 1); i >= 0; i--) {
    for (int j = (bin_alloc[1] - 1); j >= 0; j--) {
      if (sum_surf[id3D(i, j, 0)] > 0 && sum_surf[id3D(i, j, 1)] > 0) {
        max[0] = i + 1; // highest point, so add 1 to go from 0 to max-1
        max[1] = j + 1; //
        break;
      }
    }
    if (max[0] > -1) {
      break;
    }
  } //}}}

  // write surface to output file //{{{
  PrintByline(file_surf, argc, argv);
  // print legend
  FILE *out = OpenFile(file_surf, "a");
  char a[3] = {'x', 'y', 'z'};
  fprintf(out, "# (1) %c coordinate; (2) %c coordinate;", a[map[0]], a[map[1]]);
  fprintf(out, " (3) surface 1; (4) surface 2; (5) average surface\n");
  for (int i = 0; i < max[0]; i++) {
    for (int j = 0; j < max[1]; j++) {
      double surface[3];
      surface[0] = sum_surf[id3D(i, j, 0)] / values[id3D(i, j, 0)];
      surface[1] = sum_surf[id3D(i, j, 1)] / values[id3D(i, j, 1)];
      surface[2] = (surface[0] + surface[1]) / 2;

      fprintf(out, "%10.4f %10.4f %10.4f %10.4f %10.4f\n", i*width, j*width,
              surface[0], surface[1], surface[2]);
    }
    putc('\n', out);
  }
  fclose(out); //}}}

  // calculate total area as a sum of areas of triangles //{{{
  if (opt->area_file[0] != '\0') {
    for (int i = 0; i < max[0]; i++) {
      values[id3D(i, max[1]-1, 0)] = values[id3D(i, 0, 0)];
      values[id3D(i, max[1]-1, 1)] = values[id3D(i, 0, 1)];
      sum_surf[id3D(i, max[1]-1, 0)] = sum_surf[id3D(i, 0, 0)];
      sum_surf[id3D(i, max[1]-1, 1)] = sum_surf[id3D(i, 0, 1)];
    }
    for (int j = 0; j < max[1]; j++) {
      values[id3D(max[0]-1, j, 0)] = values[id3D(0, j, 0)];
      values[id3D(max[0]-1, j, 1)] = values[id3D(0, j, 1)];
      sum_surf[id3D(max[0]-1, j, 0)] = sum_surf[id3D(0, j, 0)];
      sum_surf[id3D(max[0]-1, j, 1)] = sum_surf[id3D(0, j, 1)];
    }
    double area[3] = {0, 0, 0};
    int triangles[3] = {0, 0, 0};
    for (int i = 0; i < (max[0] - 1); i++) {
      for (int j = 0; j < (max[1] - 1); j++) {
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
    int n_triangles = (max[0] - 1) * (max[1] - 1) * 2;
    double avg_triangle[3];
    for (int dd = 0; dd < 3; dd++) {
      avg_triangle[dd] = area[dd] / triangles[dd];
      area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
    }
    double Length_area = System.Box.Length[0] * System.Box.Length[1];
    double width_area = (max[0] - 1) * (max[1] - 1) * SQR(width);
    FILE *out = OpenFile(opt->area_file, "a");
    fprintf(out, "# average: (1) surface 1");
    fprintf(out, "; (2) surface 2");
    fprintf(out, "; (3) middle surface\n");
    fprintf(out, "# %lf %lf %lf\n", area[0] * Length_area / width_area,
                                    area[1] * Length_area / width_area,
                                    area[2] * Length_area / width_area);
    fclose(out);
  } //}}}

  // write distribution of widths (-w option)
  if (distr_width > 0) {
    long int norm = 0;
    int min = 0, max = distr_bins;
    for (int i = 0; i < distr_bins; i++) {
      norm += distr[i];
      if (distr[i] > 0 && min == 0) {
        min = i;
      }
      if (distr[distr_bins-1-i] > 0 && max == distr_bins) {
        max = distr_bins - i;
      }
    }

    PrintByline(opt->distr_file, argc, argv);
    out = OpenFile(opt->distr_file, "a");
    fprintf(out, "# (1) distance; (2) distribution\n");
    for (int i = min; i < max; i++) {
      double dist = distr_width * (2 * i + 1) / 2;
      fprintf(out, "%10.5f", dist);
      if (distr[i] > 0) {
        fprintf(out, " %lf", (double)(distr[i]) / norm);
      } else {
        fprintf(out, " %lf", 0.0);
      }
      putc('\n', out);
    }
    fclose(out);
  }

  FreeSystem(&System);
  free(sum_surf);
  free(values);
  if (distr_width > 0) {
    free(distr);
  }
  free(opt->bt);
  if (opt->bead_file.name[0] != '\0') {
    free(write);
  }
  free(opt);

  return 0;
}