#include "../AnalysisTools.h"

// TODO: reinstate molecular condition
// TODO: impelment -m & -bt options

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
  fprintf(ptr, "   %s <input> <width> <output> <axis> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "<axis>              calculate along x, y, or z axis\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --in              start from the box's centre "
          "instead of edges\n");
  // fprintf(ptr, "  -m <mol(s)>       molecule type(s) to use\n");
  // fprintf(ptr, "  -bt <name(s)>     bead type(s) to use\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 9, all = common + 1, count = 0,
      req_arg = 4;
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
  OptionCheck(argc, argv, count, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // arguments & options before reading system data //{{{
  // <input> - input coordinate (and structure) file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}
  // <width> - width of a single bin //{{{
  double width = 0;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}
  // <output> - output filename //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]); //}}}
  // <axis> - x, y, or z //{{{
  char axis = 'x';
  while (++count < argc && argv[count][0] != '-') {
    axis = argv[count][0];
    // Error - not x/y/z
    if (axis != 'x' && axis != 'y' && axis != 'z') {
      strcpy(ERROR_MSG, "must be 'x', 'y', or 'z'");
      PrintErrorOption("<width>");
      Help(argv[0], true, common, option);
      exit(1);
    }
  } //}}}
  bool silent, verbose, detailed;
  int start = 1, end = -1, skip = 0;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &start, &end, &skip);
  bool in = BoolOption(argc, argv, "--in"); //}}}

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed);
  COUNT *Count = &System.Count;
  BOX *Box = &System.Box;

  if (verbose) {
    VerboseOutput(System);
  }

  // set maximum/minimum as half a box length in the given direction //{{{
  double range[2];
  switch(axis) {
    case 'x':
      range[0] = 0;
      range[1] = Box->Length[0];
      break;
    case 'y':
      range[0] = 0;
      range[1] = Box->Length[1];
      break;
    case 'z':
      range[0] = 0;
      range[1] = Box->Length[2];
      break;
  } //}}}

  // number of bins //{{{
  int bins[2];
  switch(axis) {
    case 'x':
      bins[0] = Box->Length[1] / width + 1;
      bins[1] = Box->Length[2] / width + 1;
      break;
    case 'y':
      bins[0] = Box->Length[0] / width + 1;
      bins[1] = Box->Length[2] / width + 1;
      break;
    case 'z':
      bins[0] = Box->Length[0] / width + 1;
      bins[1] = Box->Length[1] / width + 1;
      break;
  }
  // TODO: this ensures enough bins when, e.g., vsf file contains pbc smaller
  //       than the coordinate file - think about how to make it safe...
  bins[0] += 10;
  bins[1] += 10;
  //}}}

// TODO: sizeof ...argh!
  // allocate memory for density arrays //{{{
  double ***surf = malloc(bins[0] * sizeof(double **)); // sum
  int ***values = malloc(bins[0] * sizeof(int **)); // number of values
  for (int i = 0; i < bins[0]; i++) {
    surf[i] = malloc(bins[1] * sizeof(double *));
    values[i] = malloc(bins[0] * sizeof(int *));
    for (int j = 0; j < bins[1]; j++) {
      surf[i][j] = calloc(2, sizeof(double));
      values[i][j] = calloc(2, sizeof(int));
    }
  } //}}}

  // open input coordinate file
  FILE *fr = OpenFile(coor_file, "r");

  // main loop //{{{
  int count_coor = 0, // count calculated timesteps
      count_used = 0, // count timesteps from the beginning
      line_count = 0; // count lines in the coor file
  while (true) {
    PrintStep(&count_coor, start, silent);
    // decide whether to use this timestep (based on -st/-sk/-e) //{{{
    bool use = false;
    if (count_coor >= start &&
        (count_coor <= end || end == -1) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    } //}}}

    if (use) { //{{{
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
    // TODO: sizeof ...argh!
      // allocate memory for temporary arrays //{{{
      double ***temp = calloc(bins[0], sizeof(double **));
      for (int i = 0; i < bins[0]; i++) {
        temp[i] = calloc(bins[1], sizeof(double *));
        for (int j = 0; j < bins[1]; j++) {
          temp[i][j] = calloc(2, sizeof(double));
          if (!in) {
            temp[i][j][0] = 0;
            switch(axis) {
              case 'x':
                temp[i][j][1] = Box->Length[0];
                break;
              case 'y':
                temp[i][j][1] = Box->Length[1];
                break;
              case 'z':
                temp[i][j][1] = Box->Length[2];
                break;
            }
          } else {
            switch(axis) {
              case 'x':
                temp[i][j][0] = Box->Length[0];
                break;
              case 'y':
                temp[i][j][0] = Box->Length[1];
                break;
              case 'z':
                temp[i][j][0] = Box->Length[2];
                break;
            }
            temp[i][j][1] = 0;
          }
        }
      } //}}}

    // TODO: check
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
        switch(axis) {
          case 'x':
            coor[0] = b->Position[1];
            coor[1] = b->Position[2];
            coor[2] = b->Position[0];
            break;
          case 'y':
            coor[0] = b->Position[0];
            coor[1] = b->Position[2];
            coor[2] = b->Position[1];
            break;
          case 'z':
            coor[0] = b->Position[0];
            coor[1] = b->Position[1];
            coor[2] = b->Position[2];
            break;
        }
        int bin[2];
        bin[0] = coor[0] / width;
        bin[1] = coor[1] / width;
        if (bin[0] > bins[0] || bin[1] > bins[1]) {
          printf("%lf %lf %lf\n", b->Position[0], b->Position[1], b->Position[2]);
          printf("%d of %d\n", bin[0], bins[0]);
          printf("%d of %d\n", bin[1], bins[1]);
          PrintBox(*Box);
          PrintBox(System.Box);
        }
        if (!in) { // go from box centre to edges
          if (coor[2] >= temp[bin[0]][bin[1]][0] &&
              coor[2] >= range[0] &&
              coor[2] <= ((range[0]+range[1])/2)) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] <= temp[bin[0]][bin[1]][1] &&
              coor[2] >= ((range[0]+range[1])/2) &&
              coor[2] <= range[1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        } else if (coor[2] >= range[0] &&
                   coor[2] <= range[1]) { // go from box edges to centre
          if (coor[2] <= temp[bin[0]][bin[1]][0]) {
            temp[bin[0]][bin[1]][0] = coor[2];
          }
          if (coor[2] >= temp[bin[0]][bin[1]][1]) {
            temp[bin[0]][bin[1]][1] = coor[2];
          }
        }
      } //}}}

      // add to sums //{{{
      for (int i = 0; i < bins[0]; i++) {
        for (int j = 0; j < bins[1]; j++) {
          if (!in) {
            if (temp[i][j][0] > 0) {
              surf[i][j][0] += temp[i][j][0];
              values[i][j][0]++;
            }
            switch(axis) {
              case 'x':
                if (temp[i][j][1] < Box->Length[0]) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
              case 'y':
                if (temp[i][j][1] < Box->Length[1]) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
              case 'z':
                if (temp[i][j][1] < Box->Length[2]) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
            }
          } else {
            switch(axis) {
              case 'x':
                if (temp[i][j][0] != (Box->Length[0]/2)) {
                  surf[i][j][0] += temp[i][j][0];
                  values[i][j][0]++;
                }
                if (temp[i][j][1] != (Box->Length[0]/2)) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
              case 'y':
                if (temp[i][j][0] != (Box->Length[1]/2)) {
                  surf[i][j][0] += temp[i][j][0];
                  values[i][j][0]++;
                }
                if (temp[i][j][1] != (Box->Length[1]/2)) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
              case 'z':
                if (temp[i][j][0] < Box->Length[2]) {
                  surf[i][j][0] += temp[i][j][0];
                  values[i][j][0]++;
                }
                if (temp[i][j][1] > 0) {
                  surf[i][j][1] += temp[i][j][1];
                  values[i][j][1]++;
                }
                break;
            }
          }
        }
      } //}}}

      // free the temporary array //{{{
      for (int i = 0; i < bins[0]; i++) {
        for (int j = 0; j < bins[1]; j++) {
          free(temp[i][j]);
        }
        free(temp[i]);
      }
      free(temp); //}}}
      //}}}
    } else { //{{{
      if (!SkipTimestep(coor_type, fr, coor_file,
                        struct_file, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // write surface to output file //{{{
  PrintByline(output, argc, argv);
  // print legend
  FILE *out = OpenFile(output, "w");
  switch(axis) {
    case 'x':
      fprintf(out, "# (1) y coordinate; (2) z coordinate;");
      break;
    case 'y':
      fprintf(out, "# (1) x coordinate; (2) z coordinate;");
      break;
    case 'z':
      fprintf(out, "# (1) x coordinate; (2) y coordinate;");
      break;
  }
  fprintf(out, " (3) surface 1; (4) surface 2\n");

  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      fprintf(out, "%7.4f", width*(2*i+1)/2);
      fprintf(out, " %7.4f", width*(2*j+1)/2);
      if (values[i][j][0] > 0) {
        fprintf(out, " %7.4f", surf[i][j][0]/values[i][j][0]);
      } else {
        fprintf(out, "       ?");
      }
      if (values[i][j][1] > 0) {
        fprintf(out, " %7.4f", surf[i][j][1]/values[i][j][1]);
      } else {
        fprintf(out, "       ?");
      }
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
  fclose(out); //}}}

  // calculate total area as a sum of areas of triangles
  double sum_area[2] = {0, 0};
  int triangles[2] = {0, 0};
  count = 0;
  int bins_true[2];
  bins_true[0] = System.Box.Length[0] / width + 1;
  bins_true[1] = System.Box.Length[1] / width + 1;
  for (int i = 0; i < bins_true[0]; i++) {
    values[i][bins_true[1]-1][0] = values[i][0][0];
    values[i][bins_true[1]-1][1] = values[i][0][1];
    surf[i][bins_true[1]-1][0] = surf[i][0][0];
    surf[i][bins_true[1]-1][1] = surf[i][0][1];
  }
  for (int j = 0; j < bins_true[1]; j++) {
    values[bins_true[0]-1][j][0] = values[0][j][0];
    values[bins_true[0]-1][j][1] = values[0][j][1];
    surf[bins_true[0]-1][j][0] = surf[0][j][0];
    surf[bins_true[0]-1][j][1] = surf[0][j][1];
  }
  for (int i = 0; i < (bins_true[0] - 1); i++) {
    for (int j = 0; j < (bins_true[1] - 1); j++) {
      count++;
      double remainder[2];
      if (i == (bins_true[0] - 2)) {
        remainder[0] = System.Box.Length[0] - width * (i + 1);
      } else {
        remainder[0] = width;
      }
      if (j == (bins_true[1] - 2)) {
        remainder[1] = System.Box.Length[1] - width * (j + 1);
      } else {
        remainder[1] = width;
      }
      // printf("%lf %lf\n", remainder[0], remainder[1]);
      // first surface //{{{
      // first triangle
      if (values[i][j][0] > 0 && values[i+1][j][0] > 0 && values[i+1][j+1][0] > 0) {
        double A[3], B[3], C[3];
        A[0] = 0;
        A[1] = 0;
        A[2] = surf[i][j][0] / values[i][j][0];
        B[0] = remainder[0];
        B[1] = 0;
        B[2] = surf[i+1][j][0] / values[i+1][j][0];
        C[0] = remainder[0];
        C[1] = remainder[1];
        C[2] = surf[i+1][j+1][0] / values[i+1][j+1][0];
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
        double area = sqrt(s * (s - a) * (s - b) * (s - c));
        sum_area[0] += area;
        triangles[0]++;
      }
      // second triangle
      if (values[i][j][0] > 0 && values[i][j+1][0] > 0 && values[i+1][j+1][0] > 0) {
        double A[3], B[3], C[3];
        A[0] = 0;
        A[1] = 0;
        A[2] = surf[i][j][0] / values[i][j][0];
        B[0] = 0;
        B[1] = remainder[1];
        B[2] = surf[i][j+1][0] / values[i][j+1][0];
        C[0] = remainder[0];
        C[1] = remainder[1];
        C[2] = surf[i+1][j+1][0] / values[i+1][j+1][0];
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
        double area = sqrt(s * (s - a) * (s - b) * (s - c));
        sum_area[0] += area;
        triangles[0]++;
      }
      //}}}
      // second surface //{{{
      // first triangle
      if (values[i][j][1] > 0 && values[i+1][j][1] > 0 && values[i+1][j+1][1] > 0) {
        double A[3], B[3], C[3];
        A[0] = 0;
        A[1] = 0;
        A[2] = surf[i][j][1] / values[i][j][1];
        B[0] = remainder[0];
        B[1] = 0;
        B[2] = surf[i+1][j][1] / values[i+1][j][1];
        C[0] = remainder[0];
        C[1] = remainder[1];
        C[2] = surf[i+1][j+1][1] / values[i+1][j+1][1];
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
        double area = sqrt(s * (s - a) * (s - b) * (s - c));
        sum_area[1] += area;
        triangles[1]++;
      }
      // second triangle
      if (values[i][j][1] > 0 && values[i][j+1][1] > 0 && values[i+1][j+1][1] > 0) {
        double A[3], B[3], C[3];
        A[0] = 0;
        A[1] = 0;
        A[2] = surf[i][j][1] / values[i][j][1];
        B[0] = 0;
        B[1] = remainder[1];
        B[2] = surf[i][j+1][1] / values[i][j+1][1];
        C[0] = remainder[0];
        C[1] = remainder[1];
        C[2] = surf[i+1][j+1][1] / values[i+1][j+1][1];
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
        double area = sqrt(s * (s - a) * (s - b) * (s - c));
        sum_area[1] += area;
        triangles[1]++;
      }
      //}}}
    }
  }

  // int n_triangles = (bins[0] - 1) * (bins[1] - 1)  * 2;
  int n_triangles = (bins_true[0] - 1) * (bins_true[1] -1) * 2;
  double avg_triangle[2];
  avg_triangle[0] = sum_area[0] / triangles[0];
  avg_triangle[1] = sum_area[1] / triangles[1];
  sum_area[0] += avg_triangle[0] * (n_triangles - triangles[0]);
  sum_area[1] += avg_triangle[1] * (n_triangles - triangles[0]);
  // printf("area = %lf (from %d triangles out of %d)\n", sum_area[0], triangles[0], n_triangles);
  // printf("       %lf (from %d triangles out of %d)\n", sum_area[1], triangles[1], n_triangles);
  printf("%lf %lf\n", sum_area[0], sum_area[1]);

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  for (int i = 0; i < bins[0]; i++) {
    for (int j = 0; j < bins[1]; j++) {
      free(surf[i][j]);
      free(values[i][j]);
    }
    free(surf[i]);
    free(values[i]);
  }
  free(surf);
  free(values); //}}}

  return 0;
}
