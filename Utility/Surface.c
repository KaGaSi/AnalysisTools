#include "../AnalysisTools.h"

// TODO: 'central' average surface
// TODO: overall averages to <surf.txt>?

double calc_area(double A[3], double B[3], double C[3]) {
  double AB[3], AC[3], BC[3];
  for (int dd = 0; dd < 3; dd++) {
    AB[dd] = B[dd] - A[dd];
    AC[dd] = C[dd] - A[dd];
    BC[dd] = C[dd] - B[dd];
  }
  double a = 0, b = 0, c = 0;
  for (int dd = 0; dd < 3; dd++) {
    a += SQR(BC[dd]);
    b += SQR(AC[dd]);
    c += SQR(AB[dd]);
  }
  a = sqrt(a);
  b = sqrt(b);
  c = sqrt(c);
  double s = (a + b + c) / 2;
  return sqrt(s * (s - a) * (s - b) * (s - c));
}

void calc_4points(double A[3], double B[3], double C[3], double D[3],
                  double *sum_area, int *triangles) {
  if (A[0] != -1 && B[0] != -1 && D[0] != -1) {
    // printf("area: %lf\n", calc_area(A, B, D));
    *sum_area += calc_area(A, B, D);
    (*triangles)++;
  }
  if (A[0] != -1 && C[0] != -1 && D[0] != -1) {
    // printf("area: %lf\n", calc_area(A, C, D));
    *sum_area += calc_area(A, C, D);
    (*triangles)++;
  }
}

// TODO: reinstate molecular condition
// TODO: impelment -m & -bt options
// TODO: calculate proper per-timestep area even when box size changes

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
  // fprintf(ptr, "  -m <mol(s)>       molecule type(s) to use\n");
  // fprintf(ptr, "  -bt <name(s)>     bead type(s) to use\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool in;        // --in
  FILE_TYPE fout; // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 9, all = common + 1, count = 0,
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
  // strcpy(option[count++], "-m");
  // strcpy(option[count++], "-bt");
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
  opt->in = BoolOption(argc, argv, "--in"); //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, opt->c.detailed);
  COUNT *Count = &System.Count;
  double sidelength[3];
  sidelength[0] = System.Box.Length[map[0]];
  sidelength[1] = System.Box.Length[map[1]];
  sidelength[2] = System.Box.Length[axis];


  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // set maximum/minimum as box length in the given direction
  double range[2];
  range[0] = 0;
  range[1] = sidelength[2];

  // number of bins //{{{
  int bins_alloc[2];
  bins_alloc[0] = sidelength[0] / width * 10;
  bins_alloc[1] = sidelength[1] / width * 10;
  // // for calculation of area - TODO: change
  // int bins_true[2];
  // bins_true[0] = sidelength[0] / width + 1;
  // bins_true[1] = sidelength[0] / width + 1;
  // // *1000 to ensure width is integer as otherwise it can wreak havoc
  // if ((fmod(sidelength[0] * 1000, width * 1000) / 1000) > 0.001) {
  //   bins_true[0]++;
  // }
  // if ((fmod(sidelength[1] * 1000, width * 1000) / 1000) > 0.001) {
  //   bins_true[1]++;
  // }
  // find lowest and highest coordinates for the topmost bins
  double lo[2], hi[2];
  lo[0] = (int)(sidelength[0] / width) * width;
  lo[1] = (int)(sidelength[1] / width) * width;
  if (fabs(lo[0]-sidelength[0]) > 0.00001) {
    hi[0] = lo[0] + width;
  } else {
    hi[0] = lo[0];
  }
  if (fabs(lo[1]-sidelength[1]) > 0.00001) {
    hi[1] = lo[1] + width;
  } else {
    hi[1] = lo[1];
  } //}}}

  // allocate memory for density arrays //{{{
  double ***surf = malloc(bins_alloc[0] * sizeof(double **)); // sum
  int ***values = malloc(bins_alloc[0] * sizeof(int **)); // number of values
  for (int i = 0; i < bins_alloc[0]; i++) {
    surf[i] = malloc(bins_alloc[1] * sizeof(double *));
    values[i] = malloc(bins_alloc[1] * sizeof(int *));
    for (int j = 0; j < bins_alloc[1]; j++) {
      surf[i][j] = calloc(2, sizeof(double));
      values[i][j] = calloc(2, sizeof(int));
    }
  } //}}}

  // open input coordinate file
  FILE *fr = OpenFile(in.coor.name, "r");

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
      line_count = 0, // count lines in the coor file
      max_bin[2] = {0, 0}; // highest bin ids (i.e., number of bins to write)
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

      // TODO: why range?
      range[0] = 0;
      range[1] = sidelength[2];

      /*
       * Recalculate number of bins for this timestep, considering that
       *   1) Box.Length might have changed
       *   2) area calculation requires one more because for the highest bin id,
       *      0th's value for surface is used
       */
      int bins_step[2];
      bins_step[0] = sidelength[0] / width + 1;
      bins_step[1] = sidelength[1] / width + 1;
      // *1000 to ensure width is integer as otherwise it can wreak havoc
      if ((fmod(sidelength[0] * 1000, width * 1000) / 1000) > 0.001) {
        bins_step[0]++;
      }
      if ((fmod(sidelength[1] * 1000, width * 1000) / 1000) > 0.001) {
        bins_step[1]++;
      }
      /*
      * Find lowest and highest coordinates for the topmost bins; these are used
      * to use beads from a periodic image should the highest bin be less wide
      * than 'width'
      */
      lo[0] = (int)(sidelength[0] / width) * width;
      lo[1] = (int)(sidelength[1] / width) * width;
      for (int aa = 0; aa < 2; aa++) {
        if (fabs(lo[aa]-sidelength[aa]) > 0.00001) {
          hi[aa] = lo[aa] + width;
        } else {
          hi[aa] = lo[aa];
        }
      }
      //}}}

      // allocate memory for temporary arrays //{{{
      double ***temp = calloc(bins_step[0], sizeof(double **));
      bool ***temp2 = calloc(bins_step[0], sizeof(bool **));
      for (int i = 0; i < bins_step[0]; i++) {
        temp[i] = calloc(bins_step[1], sizeof(double *));
        temp2[i] = calloc(bins_step[1], sizeof(bool *));
        for (int j = 0; j < bins_step[1]; j++) {
          temp[i][j] = calloc(2, sizeof(double));
          temp2[i][j] = calloc(2, sizeof(bool));
          if (!opt->in) {
            temp[i][j][0] = 0;
            temp[i][j][1] = sidelength[2];
          } else {
            temp[i][j][0] = sidelength[2];
            temp[i][j][1] = 0;
          }
        }
      } //}}}

      // calculate surface //{{{
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        BEAD *b = &System.Bead[id];
        // int mol = b->Molecule;
        // if (mol == -1) { // consider only beads in molecules
        //   continue;
        // }
        // BEADTYPE *btype = &System.BeadType[b->Type]; // for -bt
        // int mtype = System.Molecule[mol].Type;
        // MOLECULETYPE *mol_type = &System.MoleculeType[mtype]; // for -m
        // TODO: -m & -bt options: add codition like this
        //       if (b->Molecule != -1 && mol_type->Flag && btype->Flag)
        double coor[3];
        coor[0] = b->Position[map[0]];
        coor[1] = b->Position[map[1]];
        coor[2] = b->Position[axis];
        int bin[2];
        bin[0] = coor[0] / width;
        bin[1] = coor[1] / width;
        if (bin[0] > max_bin[0]) {
          max_bin[0] = bin[0];
        }
        if (bin[1] > max_bin[1]) {
          max_bin[1] = bin[1];
        }
        // printf("%d (%d) %d (%d)\n", bin[0], bins[0], bin[1], bins[1]);
        // error if box got bigger
        double pbc[2];
        pbc[0] = coor[0] + sidelength[0];
        pbc[1] = coor[1] + sidelength[1];

        if (!opt->in) { // go from box centre to edges
          if (coor[2] >= temp[bin[0]][bin[1]][0] &&
              coor[2] >= range[0] &&
              coor[2] <= ((range[0] + range[1])/2)) {
            temp[bin[0]][bin[1]][0] = coor[2];
            temp2[bin[0]][bin[1]][0] = true;
          }
          if (coor[2] <= temp[bin[0]][bin[1]][1] &&
              coor[2] >= ((range[0] + range[1])/2) &&
              coor[2] <= range[1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
            temp2[bin[0]][bin[1]][1] = true;
          }
          if (pbc[0] > lo[0] && pbc[0] < hi[0] &&
              pbc[1] > lo[1] && pbc[1] < hi[1]) {
            bin[0] = pbc[0] / width;
            bin[1] = pbc[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] >= temp[bin[0]][bin[1]][0] &&
                coor[2] >= range[0] &&
                coor[2] <= ((range[0] + range[1])/2)) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] <= temp[bin[0]][bin[1]][1] &&
                coor[2] >= ((range[0] + range[1])/2) &&
                coor[2] <= range[1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
          if (pbc[0] > lo[0] && pbc[0] < hi[0]) {
            bin[0] = pbc[0] / width;
            bin[1] = coor[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] >= temp[bin[0]][bin[1]][0] &&
                coor[2] >= range[0] &&
                coor[2] <= ((range[0] + range[1])/2)) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] <= temp[bin[0]][bin[1]][1] &&
                coor[2] >= ((range[0] + range[1])/2) &&
                coor[2] <= range[1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
          if (pbc[1] > lo[1] && pbc[1] < hi[1]) {
            bin[0] = coor[0] / width;
            bin[1] = pbc[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] >= temp[bin[0]][bin[1]][0] &&
                coor[2] >= range[0] &&
                coor[2] <= ((range[0] + range[1])/2)) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] <= temp[bin[0]][bin[1]][1] &&
                coor[2] >= ((range[0] + range[1])/2) &&
                coor[2] <= range[1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
        } else if (coor[2] >= range[0] &&
                   coor[2] <= range[1]) { // go from box edges to centre
          if (coor[2] <= temp[bin[0]][bin[1]][0]) {
            temp[bin[0]][bin[1]][0] = coor[2];
            temp2[bin[0]][bin[1]][0] = true;
          }
          if (coor[2] >= temp[bin[0]][bin[1]][1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
            temp2[bin[0]][bin[1]][1] = true;
          }
          if (pbc[0] > lo[0] && pbc[0] < hi[0] &&
              pbc[1] > lo[1] && pbc[1] < hi[1]) {
            bin[0] = pbc[0] / width;
            bin[1] = pbc[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] <= temp[bin[0]][bin[1]][0]) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] >= temp[bin[0]][bin[1]][1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
          if (pbc[0] > lo[0] && pbc[0] < hi[0]) {
            bin[0] = pbc[0] / width;
            bin[1] = coor[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] <= temp[bin[0]][bin[1]][0]) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] >= temp[bin[0]][bin[1]][1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
          if (pbc[1] > lo[1] && pbc[1] < hi[1]) {
            bin[0] = coor[0] / width;
            bin[1] = pbc[1] / width;
            if (bin[0] > max_bin[0]) {
              max_bin[0] = bin[0];
            }
            if (bin[1] > max_bin[1]) {
              max_bin[1] = bin[1];
            }
            if (coor[2] <= temp[bin[0]][bin[1]][0]) {
              temp[bin[0]][bin[1]][0] = coor[2];
              temp2[bin[0]][bin[1]][0] = true;
            }
            if (coor[2] >= temp[bin[0]][bin[1]][1]) {
              temp[bin[0]][bin[1]][1] = coor[2];
              temp2[bin[0]][bin[1]][1] = true;
            }
          }
        }
      } //}}}

      // add to sums //{{{
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          if (temp2[i][j][0]) {
            surf[i][j][0] += temp[i][j][0];
            values[i][j][0]++;
          }
          if (temp2[i][j][1]) {
            surf[i][j][1] += temp[i][j][1];
            values[i][j][1]++;
          }
        }
      } //}}}

      // calculate total area as a sum of areas of triangles //{{{
      double rest[2];
      for (int aa = 0; aa < 2; aa++) {
        rest[aa] = sidelength[aa] - width * (bins_step[aa] - 2);
        if (fabs(rest[0]) < 0.00001) {
          rest[aa] = width;
        }
      }
      for (int i = 0; i < bins_step[0]; i++) {
        temp[i][bins_step[1]-1][0] = temp[i][bins_step[0]-2][0] + (temp[i][0][0] - temp[i][bins_step[1]-2][0]) * rest[0] / width;
        temp[i][bins_step[1]-1][1] = temp[i][bins_step[0]-2][1] + (temp[i][0][1] - temp[i][bins_step[1]-2][1]) * rest[0] / width;
        temp2[i][bins_step[1]-1][0] = temp2[i][0][0];
        temp2[i][bins_step[1]-1][1] = temp2[i][0][1];
      }
      for (int j = 0; j < bins_step[1]; j++) {
        temp[bins_step[0]-1][j][0] = temp[bins_step[1]-2][j][0] + (temp[0][j][0] - temp[bins_step[1]-2][j][0]) * rest[1] / width;
        temp[bins_step[0]-1][j][1] = temp[bins_step[1]-2][j][1] + (temp[0][j][1] - temp[bins_step[1]-2][j][1]) * rest[1] / width;
        temp2[bins_step[0]-1][j][0] = temp2[0][j][0];
        temp2[bins_step[0]-1][j][1] = temp2[0][j][1];
      }
      double sum_area[3] = {0, 0, 0};
      int triangles[3] = {0, 0, 0};
      for (int i = 0; i < (bins_step[0] - 1); i++) {
        for (int j = 0; j < (bins_step[1] - 1); j++) {
          // find the true width of the bin (last one may be narrower)
          rest[0] = sidelength[0] - width * (bins_step[0] - 2);
          if (fabs(rest[0]) < 0.00001 || i != (bins_step[0] - 2)) {
            rest[0] = width;
          }
          rest[1] = sidelength[1] - width * (bins_step[1] - 2);
          if (fabs(rest[1]) < 0.00001 || j != (bins_step[1] - 2)) {
            rest[1] = width;
          }
          // four points defining the two triangles
          double A[3] = {-1}, B[3] = {-1}, C[3] = {-1}, D[3] = {-1};
          // top and bottom surfaces //{{{
          for (int aa = 0; aa < 2; aa++) {
            if (temp2[i][j][aa] && temp2[i+1][j+1][aa]) {
              A[0] = A[1] = 0;
              A[2] = temp[i][j][aa];
              D[0] = rest[0];
              D[1] = rest[1];
              D[2] = temp[i+1][j+1][aa];
              if (temp2[i+1][j][aa]) {
                B[0] = rest[0];
                B[1] = 0;
                B[2] = temp[i+1][j][aa];
              }
              if (temp2[i][j+1][aa]) {
                C[0] = 0;
                C[1] = rest[1];
                C[2] = temp[i][j+1][aa];
              }
              calc_4points(A, B, C, D, &sum_area[aa], &triangles[aa]);
            }
          } //}}}
          // 'middle' surface //{{{
          A[0] = B[0] = C[0] = D[0] = -1;
          if (temp2[i][j][0] && temp2[i][j][1] &&
              temp2[i+1][j+1][0] && temp2[i+1][j+1][1]) {
            A[0] = A[1] = 0;
            A[2] = (temp[i][j][0] + temp[i][j][1]) / 2;
            D[0] = rest[0];
            D[1] = rest[1];
            D[2] = (temp[i+1][j+1][0] + temp[i+1][j+1][1]) / 2;

            if (temp2[i+1][j][0] && temp2[i+1][j][1]) {
              B[0] = rest[0];
              B[1] = 0;
              B[2] = (temp[i+1][j][0] + temp[i+1][j][1]) / 2;
            }
            if (temp2[i][j+1][0] && temp2[i][j+1][1]) {
              C[0] = rest[1];
              C[1] = 0;
              C[2] = (temp[i][j+1][0] + temp[i][j+1][1]) / 2;
            }
            calc_4points(A, B, C, D, &sum_area[2], &triangles[2]);
          }
          //}}}
        }
      }

      int n_triangles = (bins_step[0] - 1) * (bins_step[1] - 1) * 2;
      double avg_triangle[3];
      for (int dd = 0; dd < 3; dd++) {
        avg_triangle[dd] = sum_area[dd] / triangles[dd];
        sum_area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
      }
      fprintf(out, "%d %lf %lf %lf\n", count_coor,
              sum_area[0], sum_area[1], sum_area[2]);
      //}}}

      // free the temporary array //{{{
      for (int i = 0; i < bins_step[0]; i++) {
        for (int j = 0; j < bins_step[1]; j++) {
          free(temp[i][j]);
          free(temp2[i][j]);
        }
        free(temp[i]);
        free(temp2[i]);
      }
      free(temp);
      free(temp2); //}}}
      //}}}
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
      if (values[i][j][0] > 0 || values[i][j][1] > 0) {
        count++;
        fprintf(out, "%7.4f", width*(2*i+1)/2);
        fprintf(out, " %7.4f", width*(2*j+1)/2);
        double surface[2];
        for (int aa = 0; aa < 2; aa++) {
          if (values[i][j][aa] > 0) {
            surface[aa] = surf[i][j][aa] / values[i][j][aa];
            fprintf(out, " %7.4f", surface[aa]);
          } else {
            fprintf(out, "       ?");
          }
        }
        if (values[i][j][0] > 0 && values[i][j][1] > 0) {
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

  // calculate total area as a sum of areas of triangles
  double sum_area[3] = {0, 0, 0};
  int triangles[3] = {0, 0, 0};
  count = 0;
  for (int i = 0; i < max_bin[0]; i++) {
    values[i][max_bin[1]-1][0] = values[i][0][0];
    values[i][max_bin[1]-1][1] = values[i][0][1];
    surf[i][max_bin[1]-1][0] = surf[i][0][0];
    surf[i][max_bin[1]-1][1] = surf[i][0][1];
  }
  for (int j = 0; j < max_bin[1]; j++) {
    values[max_bin[0]-1][j][0] = values[0][j][0];
    values[max_bin[0]-1][j][1] = values[0][j][1];
    surf[max_bin[0]-1][j][0] = surf[0][j][0];
    surf[max_bin[0]-1][j][1] = surf[0][j][1];
  }
  for (int i = 0; i < (max_bin[0] - 1); i++) {
    for (int j = 0; j < (max_bin[1] - 1); j++) {
      count++;
      // top and bottom surfaces //{{{
      double A[3] = {-1}, B[3] = {-1}, C[3] = {-1}, D[3] = {-1};
      for (int aa = 0; aa < 2; aa++) {
        if (values[i][j][aa] > 0 && values[i+1][j+1][aa] > 0) {
          A[0] = A[1] = 0;
          A[2] = surf[i][j][aa] / values[i][j][aa];
          D[0] = D[1] = width;
          D[2] = surf[i+1][j+1][aa] / values[i][j][aa];
          if (values[i+1][j][aa] > 0) {
            B[0] = width;
            B[1] = 0;
            B[2] = surf[i+1][j][aa] / values[i+1][j][aa];
          }
          if (values[i][j+1][aa] > 0) {
            C[0] = 0;
            C[1] = width;
            C[2] = surf[i][j+1][aa] / values[i][j+1][aa];
          }
          calc_4points(A, B, C, D, &sum_area[aa], &triangles[aa]);
        }
      } //}}}
      // 'middle' surface //{{{
      A[0] = B[0] = C[0] = D[0] = -1;
      if (values[i][j][0] > 0 && values[i][j][1] > 0 &&
          values[i+1][j+1][0] > 0 && values[i+1][j+1][1] > 0) {
          // values[i+1][j][0] > 0 &&
          // values[i+1][j][1] > 0 &&
        A[0] = A[1] = 0;
        A[2] = (surf[i][j][0] / values[i][j][0] +
                surf[i][j][1] / values[i][j][1]) / 2;
        D[0] = D[1] = width;
        D[2] = (surf[i+1][j+1][0] / values[i+1][j+1][0] +
                surf[i+1][j+1][1] / values[i+1][j+1][1]) / 2;

        if (values[i+1][j][0] > 0 && values[i+1][j][1] > 0) {
          B[0] = width;
          B[1] = 0;
          B[2] = (surf[i+1][j][0] / values[i+1][j][0] +
                  surf[i+1][j][1] / values[i+1][j][1]) / 2;
        }
        if (values[i][j+1][0] > 0 && values[i][j+1][1] > 0) {
          C[0] = width;
          C[1] = 0;
          C[2] = (surf[i][j+1][0] / values[i][j+1][0] +
                  surf[i][j+1][1] / values[i][j+1][1]) / 2;
        }
        calc_4points(A, B, C, D, &sum_area[2], &triangles[2]);
      }
      //}}}
    }
  }

  // int n_triangles = (bins[0] - 1) * (bins[1] - 1)  * 2;
  int n_triangles = (max_bin[0] - 1) * (max_bin[1] - 1) * 2;
  double avg_triangle[3];
  for (int dd = 0; dd < 3; dd++) {
    avg_triangle[dd] = sum_area[dd] / triangles[dd];
  }
  for (int dd = 0; dd < 3; dd++) {
    sum_area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
  }
  double Length_area = System.Box.Length[0] * System.Box.Length[1];
  double width_area = (max_bin[0] - 1) * (max_bin[1] - 1) * SQR(width);
  out = OpenFile(file_area, "a");
  fprintf(out, "# average: (1) surface 1; (2) surface 2; (3) middle surface\n");
  fprintf(out, "# %lf %lf %lf\n", sum_area[0] * Length_area / width_area,
                                  sum_area[1] * Length_area / width_area,
                                  sum_area[2] * Length_area / width_area);
  fclose(out);

  FreeSystem(&System);
  for (int i = 0; i < bins_alloc[0]; i++) {
    for (int j = 0; j < bins_alloc[1]; j++) {
      free(surf[i][j]);
      free(values[i][j]);
    }
    free(surf[i]);
    free(values[i]);
  }
  free(surf);
  free(values);
  free(opt);

  return 0;
}
