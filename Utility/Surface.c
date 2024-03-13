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
  BOX *Box = &System.Box;

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // set maximum/minimum as box length in the given direction
  double range[2];
  range[0] = 0;
  range[1] = Box->Length[axis];

  // number of bins //{{{
  int bins[2];
  bins[0] = Box->Length[map[0]] / width + 1;
  bins[1] = Box->Length[map[1]] / width + 1;
  // TODO: this ensures enough bins when, e.g., vsf file may contain pbc smaller
  //       than the coordinate file - think about how to make it safe...
  bins[0] += 10;
  bins[1] += 10;
  //}}}

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
  FILE *fr = OpenFile(in.coor.name, "r");

  // write initial stuff to the per-timestep area
  PrintByline(file_area, argc, argv);
  FILE *out = OpenFile(file_area, "a");
  fprintf(out, "# (1) timestep; (2) surface 1; (3) surface 2; "
          "(4) middle surface\n");

  // main loop //{{{
  int count_coor = 0, // count calculated timesteps
      count_used = 0, // count timesteps from the beginning
      line_count = 0; // count lines in the coor file
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

      // copy some beads to have enough beads for the (width×bin)^2 area //{{{
      double lo[2], hi[2];
      lo[0] = (int)(System.Box.Length[0] / width) * width;
      lo[1] = (int)(System.Box.Length[1] / width) * width;
      hi[0] = lo[0] + width;
      hi[1] = lo[1] + width;
      SYSTEM extra_xy = CopySystem(System);
      SYSTEM extra_x = CopySystem(System);
      SYSTEM extra_y = CopySystem(System);
      extra_xy.Count.BeadCoor = 0;
      extra_x.Count.BeadCoor = 0;
      extra_y.Count.BeadCoor = 0;
      for (int i = 0; i < Count->BeadCoor; i++) {
        int id = System.BeadCoor[i];
        extra_xy.Bead[id].InTimestep = false;
        extra_x.Bead[id].InTimestep = false;
        extra_y.Bead[id].InTimestep = false;
        double (*pos)[3] = &System.Bead[id].Position;
        double pbc[2];
        pbc[0] = (*pos)[0] + System.Box.Length[0];
        pbc[1] = (*pos)[1] + System.Box.Length[1];
        if (pbc[0] > lo[0] && pbc[0] < hi[0] &&
            pbc[1] > lo[1] && pbc[1] < hi[1]) {
          extra_xy.BeadCoor[extra_xy.Count.BeadCoor] = id;
          extra_xy.Bead[id].Position[0] = pbc[0];
          extra_xy.Bead[id].Position[1] = pbc[1];
          extra_xy.Bead[id].Position[2] = (*pos)[2];
          extra_xy.Count.BeadCoor++;
        }
        if (pbc[0] > lo[0] && pbc[0] < hi[0]) {
          extra_x.BeadCoor[extra_x.Count.BeadCoor] = id;
          extra_x.Bead[id].Position[0] = pbc[0];
          extra_x.Bead[id].Position[1] = (*pos)[1];
          extra_x.Bead[id].Position[2] = (*pos)[2];
          extra_x.Count.BeadCoor++;
        }
        if (pbc[1] > lo[1] && pbc[1] < hi[1]) {
          extra_y.BeadCoor[extra_y.Count.BeadCoor] = id;
          extra_y.Bead[id].Position[0] = (*pos)[0];
          extra_y.Bead[id].Position[1] = pbc[1];
          extra_y.Bead[id].Position[2] = (*pos)[2];
          extra_y.Count.BeadCoor++;
        }
      }
      FillInCoor(&extra_xy);
      FillInCoor(&extra_x);
      FillInCoor(&extra_y);
      SYSTEM full = CopySystem(System);
      ConcatenateSystems(&full, extra_xy, System.Box);
      ConcatenateSystems(&full, extra_x, System.Box);
      ConcatenateSystems(&full, extra_y, System.Box);
      FreeSystem(&extra_xy);
      FreeSystem(&extra_x);
      FreeSystem(&extra_y); //}}}

      // allocate memory for temporary arrays //{{{
      double ***temp = calloc(bins[0], sizeof(double **));
      for (int i = 0; i < bins[0]; i++) {
        temp[i] = calloc(bins[1], sizeof(double *));
        for (int j = 0; j < bins[1]; j++) {
          temp[i][j] = calloc(2, sizeof(double));
          if (!opt->in) {
            temp[i][j][0] = 0;
            temp[i][j][1] = Box->Length[axis];
          } else {
            temp[i][j][0] = Box->Length[axis];
            temp[i][j][1] = 0;
          }
        }
      } //}}}

    // TODO: check
    // TODO: get rid of the SYSTEM full and just add some other coordinates as
    //       well: if a bead has 'proper' coordinates, use it several times
      // calculate surface //{{{
      for (int i = 0; i < full.Count.BeadCoor; i++) {
        int id = full.BeadCoor[i];
        BEAD *b = &full.Bead[id];
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
        // TODO: useless check - would only trigger for serious pbc change
        if (bin[0] > bins[0] || bin[1] > bins[1]) {
          printf("%lf %lf %lf\n", b->Position[0], b->Position[1], b->Position[2]);
          printf("%d of %d\n", bin[0], bins[0]);
          printf("%d of %d\n", bin[1], bins[1]);
          PrintBox(*Box);
          PrintBox(System.Box);
        }
        if (!opt->in) { // go from box centre to edges
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
          // if (bin[1] == 32 && bin[0] == 32) {
          //   printf("bin 36: %d %lf %lf %lf\n", id, coor[0], coor[1], coor[2]);
          // }
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
          if (!opt->in) {
            if (temp[i][j][0] > 0) {
              surf[i][j][0] += temp[i][j][0];
              values[i][j][0]++;
            }
            if (temp[i][j][1] < Box->Length[axis]) {
              surf[i][j][1] += temp[i][j][1];
              values[i][j][1]++;
            }
          } else {
            if (temp[i][j][0] < Box->Length[axis]) {
              surf[i][j][0] += temp[i][j][0];
              values[i][j][0]++;
            }
            if (temp[i][j][1] > 0) {
              surf[i][j][1] += temp[i][j][1];
              values[i][j][1]++;
            }
          }
        }
      } //}}}

      // calculate total area as a sum of areas of triangles //{{{
      double sum_area[3] = {0, 0, 0};
      int triangles[3] = {0, 0, 0};
      count = 0;
      // TODO: is the bins_true correct? E.g., Box 10, width a) = 1.0, b) = 1.1
      //       a) bins_true = 11
      //       b) bins_true = 12
      //       a) last temp with data: temp[9]
      //       b) last temp with data: temp[10]
      //       ...is that so & is that right?
      //       ...consider we're using the SYSTEM full; so I thinks it's right
      int bins_true[2];
      bins_true[0] = System.Box.Length[0] / width + 1;
      bins_true[1] = System.Box.Length[1] / width + 1;
      // printf("%d\n", bins_true[0]);
      double remainder[2];
      remainder[0] = fmod(System.Box.Length[0],width);
      // TODO: no idea why fmod() sometimes returns 'width' (i.e., the divisor)
      if (fabs(remainder[0]) > 0.001 && fabs(remainder[0]-width) > 0.001) {
        bins_true[0]++;
      }
      // printf("%lf %lf\n", fabs(remainder[0]), fabs(remainder[0]-width));
      remainder[1] = fmod(System.Box.Length[1],width);
      if (fabs(remainder[1]) > 0.001 && fabs(remainder[1]-width) > 0.001) {
        bins_true[1]++;
      }
      // printf("%d %lf %lf %d %lf\n", bins_true[0], System.Box.Length[0], width,
      //        (int)(System.Box.Length[0]/width), System.Box.Length[0]/width);
      for (int i = 0; i < bins_true[0]; i++) {
        temp[i][bins_true[1]-1][0] = temp[i][0][0];
        temp[i][bins_true[1]-1][1] = temp[i][0][1];
      }
      for (int j = 0; j < bins_true[1]; j++) {
        temp[bins_true[0]-1][j][0] = temp[0][j][0];
        temp[bins_true[0]-1][j][1] = temp[0][j][1];
      }
      for (int i = 0; i < (bins_true[0] - 1); i++) {
        for (int j = 0; j < (bins_true[1] - 1); j++) {
          count++;
          // TODO: assumes in=true !!!
          // first surface //{{{
          // first triangle
          if (temp[i][j][0] < Box->Length[axis] &&
              temp[i+1][j][0] < Box->Length[axis] &&
              temp[i+1][j+1][0] < Box->Length[axis]) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = temp[i][j][0];
            // A[2] = 0;
            B[0] = width;
            B[2] = temp[i+1][j][0];
            // B[2] = 0;
            C[0] = width;
            C[1] = width;
            C[2] = temp[i+1][j+1][0];
            // C[2] = 0;
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
          if (temp[i][j][0] < Box->Length[axis] &&
              temp[i][j+1][0] < Box->Length[axis] &&
              temp[i+1][j+1][0] < Box->Length[axis]) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = temp[i][j][0];
            // A[2] = 0;
            B[1] = width;
            B[2] = temp[i][j+1][0];
            // B[2] = 0;
            C[0] = width;
            C[1] = width;
            C[2] = temp[i+1][j+1][0];
            // C[2] = 0;
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
          if (temp[i][j][1] > 0 &&
              temp[i+1][j][1] > 0 &&
              temp[i+1][j+1][1] > 0) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = temp[i][j][1];
            B[0] = width;
            B[2] = temp[i+1][j][1];
            C[0] = width;
            C[1] = width;
            C[2] = temp[i+1][j+1][1];
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
          if (temp[i][j][1] > 0 &&
              temp[i][j+1][1] > 0 &&
              temp[i+1][j+1][1] > 0) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = temp[i][j][1];
            B[1] = width;
            B[2] = temp[i][j+1][1];
            C[0] = width;
            C[1] = width;
            C[2] = temp[i+1][j+1][1];
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
          // 'middle' surface //{{{
          // first triangle
          if (temp[i][j][0] < Box->Length[axis] &&
              temp[i+1][j][0] < Box->Length[axis] &&
              temp[i+1][j+1][0] < Box->Length[axis] &&
              temp[i][j][1] > 0 &&
              temp[i+1][j][1] > 0 &&
              temp[i+1][j+1][1] > 0) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = (temp[i][j][0] + temp[i][j][1]) / 2;
            B[0] = width;
            B[2] = (temp[i+1][j][0] + temp[i+1][j][1]) / 2;
            C[0] = width;
            C[1] = width;
            C[2] = (temp[i+1][j+1][0] + temp[i+1][j+1][1]) / 2;
            // C[2] = 0;
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
            sum_area[2] += area;
            triangles[2]++;
          }
          // second triangle
          if (temp[i][j][0] < Box->Length[axis] &&
              temp[i][j+1][0] < Box->Length[axis] &&
              temp[i+1][j+1][0] < Box->Length[axis] &&
              temp[i][j][1] > 0 &&
              temp[i][j+1][1] > 0 &&
              temp[i+1][j+1][1] > 0) {
            double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
            A[2] = (temp[i][j][0] + temp[i][j][1]) / 2;
            B[1] = width;
            B[2] = (temp[i][j+1][0] + temp[i][j+1][1]) / 2;
            C[0] = width;
            C[1] = width;
            C[2] = (temp[i+1][j+1][0] + temp[i+1][j+1][1]) / 2;
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
            sum_area[2] += area;
            triangles[2]++;
          }
          //}}}
        }
      }

      // int n_triangles = (bins[0] - 1) * (bins[1] - 1)  * 2;
      int n_triangles = (bins_true[0] - 1) * (bins_true[1] - 1) * 2;
      double avg_triangle[3];
      avg_triangle[0] = sum_area[0] / triangles[0];
      avg_triangle[1] = sum_area[1] / triangles[1];
      avg_triangle[2] = sum_area[2] / triangles[2];
      sum_area[0] += avg_triangle[0] * (n_triangles - triangles[0]);
      sum_area[1] += avg_triangle[1] * (n_triangles - triangles[1]);
      sum_area[2] += avg_triangle[2] * (n_triangles - triangles[2]);
      double Length_area = System.Box.Length[0] * System.Box.Length[1];
      double width_area = (bins_true[0] - 1) * (bins_true[1] - 1) * SQR(width);
      fprintf(out, "%d %lf %lf %lf\n", count_coor,
                                       sum_area[0]*Length_area/width_area,
                                       sum_area[1]*Length_area/width_area,
                                       sum_area[2]*Length_area/width_area);
      //}}}

      // free the temporary array //{{{
      for (int i = 0; i < bins[0]; i++) {
        for (int j = 0; j < bins[1]; j++) {
          free(temp[i][j]);
        }
        free(temp[i]);
      }
      free(temp);
      FreeSystem(&full); //}}}
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
  double sum_area[3] = {0, 0, 0};
  int triangles[3] = {0, 0, 0};
  count = 0;
  // TODO: is the bins_true correct? E.g., Box 36, width a) = 1.0, b) = 1.1
  //       a) bins_true = 36 / 1 + 1 = 37
  //       b) bins_true = 36 / 1.1 + 1 = 33, than + 1 from fmod = 34
  //       a) last temp with data: coor 35-36, i.e. temp[36]
  //       b) last temp with data: coor 35.2-36, i.e. temp[32]
  //       ...is that so & is that right?
  //       ...consider we're using the SYSTEM full; so I think it's right
  //       ...although temp[33] is empty
  //       ...from ~/FA-bilayer/Hair-C16/36_to_32/3 &
  // Surface test.lammpstrj 1.1 surf.txt area.txt z --in -st 2 -e 6 -sk 20 --silent
  int bins_true[2];
  bins_true[0] = System.Box.Length[0] / width + 1;
  bins_true[1] = System.Box.Length[1] / width + 1;
  // TODO: no idea why fmod() sometimes returns 'width' (i.e., the divisor)
  double remainder[2];
  remainder[0] = fmod(System.Box.Length[0],width);
  if (fabs(remainder[0]) > 0.001 && fabs(remainder[0]-width) > 0.001) {
    bins_true[0]++;
  }
  remainder[1] = fmod(System.Box.Length[1],width);
  if (fabs(remainder[1]) > 0.001 && fabs(remainder[1]-width) > 0.001) {
    bins_true[1]++;
  }
  for (int i = 0; i < bins_true[0]; i++) {
    values[i][bins_true[1]-1][0] = values[i][0][0];
    values[i][bins_true[1]-1][1] = values[i][0][1];
    surf[i][bins_true[1]-1][0] = surf[i][0][0];
    surf[i][bins_true[1]-1][1] = surf[i][0][1];
    // printf("%d %d\n", i, values[i][bins_true[1]-1][0]);
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
      // first surface //{{{
      // first triangle
      if (values[i][j][0] > 0 &&
          values[i+1][j][0] > 0 &&
          values[i+1][j+1][0] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] = surf[i][j][0] / values[i][j][0];
        B[0] = width;
        B[2] = surf[i+1][j][0] / values[i+1][j][0];
        C[0] = width;
        C[1] = width;
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
      if (values[i][j][0] > 0 &&
          values[i][j+1][0] > 0 &&
          values[i+1][j+1][0] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] = surf[i][j][0] / values[i][j][0];
        B[0] = 0;
        B[1] = width;
        B[2] = surf[i][j+1][0] / values[i][j+1][0];
        C[0] = width;
        C[1] = width;
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
      if (values[i][j][1] > 0 &&
          values[i+1][j][1] > 0 &&
          values[i+1][j+1][1] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] = surf[i][j][1] / values[i][j][1];
        B[0] = width;
        B[2] = surf[i+1][j][1] / values[i+1][j][1];
        C[0] = width;
        C[1] = width;
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
      if (values[i][j][1] > 0 &&
          values[i][j+1][1] > 0 &&
          values[i+1][j+1][1] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] = surf[i][j][1] / values[i][j][1];
        B[1] = width;
        B[2] = surf[i][j+1][1] / values[i][j+1][1];
        C[0] = width;
        C[1] = width;
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
      // 'middle' surface //{{{
      // first triangle
      if (values[i][j][0] > 0 &&
          values[i+1][j][0] > 0 &&
          values[i+1][j+1][0] > 0 &&
          values[i][j][1] > 0 &&
          values[i+1][j][1] > 0 &&
          values[i+1][j+1][1] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] =  surf[i][j][0] / values[i][j][0];
        A[2] += surf[i][j][1] / values[i][j][1];
        A[2] /= 2;
        B[0] = width;
        B[2] =  surf[i+1][j][0] / values[i+1][j][0];
        B[2] += surf[i+1][j][1] / values[i+1][j][1];
        B[2] /= 2;
        C[0] = width;
        C[1] = width;
        C[2] =  surf[i+1][j+1][0] / values[i+1][j+1][0];
        C[2] += surf[i+1][j+1][1] / values[i+1][j+1][1];
        C[2] /= 2;
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
        sum_area[2] += area;
        triangles[2]++;
      }
      // second triangle
      if (values[i][j][0] > 0 &&
          values[i][j+1][0] > 0 &&
          values[i+1][j+1][0] > 0 &&
          values[i][j][1] > 0 &&
          values[i][j+1][1] > 0 &&
          values[i+1][j+1][1] > 0) {
        double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, C[3] = {0, 0, 0};
        A[2] =  surf[i][j][0] / values[i][j][0];
        A[2] += surf[i][j][1] / values[i][j][1];
        A[2] /= 2;
        B[1] = width;
        B[2] =  surf[i][j+1][0] / values[i][j+1][0];
        B[2] += surf[i][j+1][1] / values[i][j+1][1];
        B[2] /= 2;
        C[0] = width;
        C[1] = width;
        C[2] =  surf[i+1][j+1][0] / values[i+1][j+1][0];
        C[2] += surf[i+1][j+1][1] / values[i+1][j+1][1];
        C[2] /= 2;
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
        sum_area[2] += area;
        triangles[2]++;
      }
      //}}}
    }
  }

  // int n_triangles = (bins[0] - 1) * (bins[1] - 1)  * 2;
  int n_triangles = (bins_true[0] - 1) * (bins_true[1] - 1) * 2;
  double avg_triangle[3];
  for (int dd = 0; dd < 3; dd++) {
    avg_triangle[dd] = sum_area[dd] / triangles[dd];
  }
  for (int dd = 0; dd < 3; dd++) {
    sum_area[dd] += avg_triangle[dd] * (n_triangles - triangles[dd]);
  }
  double Length_area = System.Box.Length[0] * System.Box.Length[1];
  double width_area = (bins_true[0] - 1) * (bins_true[1] - 1) * SQR(width);
  out = OpenFile(file_area, "a");
  fprintf(out, "# average: (1) surface 1; (2) surface 2; (3) middle surface\n");
  fprintf(out, "# %lf %lf %lf\n", sum_area[0]*Length_area/width_area,
                                  sum_area[1]*Length_area/width_area,
                                  sum_area[2]*Length_area/width_area);
  fclose(out);

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
  free(values);
  free(opt);

  return 0;
}
