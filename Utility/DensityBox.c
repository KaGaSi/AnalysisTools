#include "../AnalysisTools.h"

// TODO: possible changing box size: make bins' width variable, keeping their
//       number, and work in relative coordinates (relative to instantaneous
//       dimensions) throughout the code

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
DensityBox utility calculates number \
density for all bead types in the direction of all axes (x, y, and z). \
The utility works properly only for orthogonal boxes that do not change size.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single bin\n");
  fprintf(ptr, "<output>            output density files (automatic ending "
          "-x.rho, -y.rho, and -z.rho)\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -x <name(s)>      exclude specified molecule(s)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 9, all = common + 1, count = 0,
      req_arg = 3;
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
  strcpy(option[count++], "-x");
  OptionCheck(argc, argv, count, req_arg, common, all, option);
  //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  snprintf(in.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <outputt> - filename
  char fout_rho[LINE] = "";
  snprintf(fout_rho, LINE, "%s", argv[++count]);
  fout_rho[LINE-7] = '\0'; // for adding -<axis>.rho

  // options before reading system data
  opt->c = CommonOptions(argc, argv, LINE, in);

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, opt->c.detailed);
  COUNT *Count = &System.Count;
  BOX *box = &System.Box;

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, &System)) {
    exit(1);
  } //}}}

  // number of bins //{{{
  if (box->Volume == -1) {
    strcpy(ERROR_MSG, "missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }
  double bin[3];
  // TODO: *3 to assume box change of at most thrice as big
  //       probably change from width to number of bins per box?
  for (int dd = 0; dd < 3; dd++) {
    bin[dd] = ceil(box->Length[dd] / width) * 3;
  } //}}}

  // allocate memory for arrays //{{{
  // just check if the bead type is at all present in the calculation
  bool *n_beads = calloc(Count->BeadType, sizeof *n_beads);
  long int **rho[3]; // = malloc(3 * sizeof ***rho);
  for (int dd = 0; dd < 3; dd++) {
    rho[dd] = malloc(Count->BeadType * sizeof **rho);
  }
  for (int j = 0; j < Count->BeadType; j++) {
    for (int dd = 0; dd < 3; dd++) {
      rho[dd][j] = calloc(bin[dd], sizeof *rho[dd][j]);
    }
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (count_coor >= opt->c.start &&
        (count_coor <= opt->c.end || opt->c.end == -1) &&
        ((count_coor - opt->c.start) % opt->c.skip) == 0) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      // allocate memory for temporary density arrays
      int **temp_rho[3]; // = malloc(3 * sizeof(int **));
      for (int dd = 0; dd < 3; dd++) {
        temp_rho[dd] = malloc(Count->BeadType * sizeof *temp_rho[dd]);
        for (int i = 0; i < Count->BeadType; i++) {
          temp_rho[dd][i] = calloc(bin[dd], sizeof *temp_rho[dd][i]);
        }
      }

      // calculate densities //{{{
      for (int i = 0; i < Count->BeadCoor; i++) {
        use = true;
        int id = System.BeadCoor[i];
        BEAD *bead = &System.Bead[id];
        int mol = bead->Molecule;
        if (mol != -1) { // do not use excluded molecules (-x option)
          int mtype = System.Molecule[mol].Type;
          use = System.MoleculeType[mtype].Flag;
        }
        if (use) {
          n_beads[bead->Type] = true;
          for (int dd = 0; dd < 3; dd++) {
            int j = bead->Position[dd] / width;
            temp_rho[dd][bead->Type][j]++;
          }
        }
      } //}}}
      // add from temporary density arrays to global density arrays
      for (int j = 0; j < Count->BeadType; j++) {
        for (int dd = 0; dd < 3; dd++) {
          for (int k = 0; k < bin[dd]-1; k++) {
            rho[dd][j][k] += temp_rho[dd][j][k];
          }
        }
      }
      // free temporary density array
      for (int dd = 0; dd < 3; dd++) {
        for (int j = 0; j < Count->BeadType; j++) {
          free(temp_rho[dd][j]);
        }
        free(temp_rho[dd]);
      }
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
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

  // write densities to output file(s) //{{{
  for (int ax = 0; ax < 3; ax++) {
    // axis-based variables
    double volume = width, size = -1;
    int n = 0; // number of bins
    char axis;
    if (ax == 0) {
      axis = 'x';
      size = box->Length[0];
      volume *= box->Length[1] * box->Length[2];
      n = bin[0];
    } else if (ax == 1) {
      axis = 'y';
      size = box->Length[1];
      volume *= box->Length[0] * box->Length[2];
      n = bin[1];
    } else {
      axis = 'z';
      size = box->Length[2];
      volume *= box->Length[0] * box->Length[1];
      n = bin[2];
    }
    char file[LINE]; // filename <output>-<axis>.rho
    if (snprintf(file, LINE, "%s-%c.rho", fout_rho, axis) < 0) {
      ErrorSnprintf();
    }
    // write initial stuff to output density file
    PrintByline(file, argc, argv);
    FILE *fw = OpenFile(file, "a");
    // print bead type names to output file
    fprintf(fw, "# columns: (1) distance");
    count = 1;
    for (int i = 0; i < Count->BeadType; i++) {
      if (n_beads[i]) {
        count++;
        fprintf(fw, "; (%d) %s", count, System.BeadType[i].Name);
      }
    }
    putc('\n', fw);
    // write rdf
    for (int i = 0; i < (n - 1); i++) {
      double dist = width * (2 * i + 1) / 2;
      if (dist > size) { // write only til the max box size
        break;
      }
      fprintf(fw, "%7.3f", dist); // absolute distance
      for (int j = 0; j < Count->BeadType; j++) {
        if (n_beads[j]){
          double temp_rho = rho[ax][j][i] / (volume * count_used);
          fprintf(fw, " %10f", temp_rho);
        }
      }
      putc('\n',fw);
    }
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  for (int dd = 0; dd < 3; dd++) {
    for (int j = 0; j < Count->BeadType; j++) {
      free(rho[dd][j]);
    }
    free(rho[dd]);
  }
  free(n_beads);
  free(opt); //}}}

  return 0;
}
