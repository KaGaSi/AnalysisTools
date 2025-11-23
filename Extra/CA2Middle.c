#include "../src/AnalysisTools.h"
#include <stdbool.h>

char *name = "CA2";

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Take coordinate file with molecule called CA2 and returns \
a file with geometric centre of CA2.\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <max> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output coordinate file\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  int dim[4]; // -d: [0]...1D/2D/3D; rest axes (0-2)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 0, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);
  //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // seed random number generator
  srand(time(0));

  InitOutputCoorFile(fout, System, argc, argv);

  int mtype = FindMoleculeName(name, System);
  MOLECULETYPE *mt = NULL;
  if (mtype != -1) {
    mt = &System.MoleculeType[mtype];
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
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (fout.type == LDATA_FILE && count_used == 1) {
        err_msg("only one timestep can be saved to lammps data file");
        PrintWarnFile(fout.name, "\0", "\0");
        count_coor--;
        break;
      }
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, false, true);
      bool *write = malloc(Count->BeadCoor * sizeof *write);
      InitBoolArray(write, Count->BeadCoor, true);
      if (mtype != -1) {
        for (int i = 0; i < mt->Number; i++) {
          int id = mt->Index[i];
          int beads[mt->nBeads];
          count = 0;
          for (int j = 0; j < mt->nBeads; j++) {
            int b = System.Molecule[id].Bead[j];
            if (System.Bead[b].InTimestep) {
              beads[count] = b;
              count++;
            }
            if (j > 0) {
              write[b] = false;
            }
          }
          double gc[3];
          GeomCentre(count, beads, System.Bead, gc);
          for (int dd = 0; dd < 3; dd++) {
            System.Bead[System.Molecule[id].Bead[0]].Position.v[dd] = gc[dd];
          }
          for (int j = 1; j < mt->nBeads; j++) {
          }
        }
      }
      WriteTimestep(fout, System, count_coor, write, argc, argv);
      free(write);
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
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

  // free memory - to make valgrind happy //{{{
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
