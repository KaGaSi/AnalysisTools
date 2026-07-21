#include "../src/AnalysisTools.h"

char *name = "CA2";

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Take coordinate file with molecule called CA2 and returns "
  "a file with geometric centre of CA2.",

  "Usage: CA2Middle <input> <max> <output> [options]",
  .args = 2, // number of mandatory arguments
  .all = 10, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<output>", nullptr, "output coordinate file", OPT_ARG},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  int dim[4]; // -d: [0]...1D/2D/3D; rest axes (0-2)
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  int count = 0;
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
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // seed random number generator
  srand(time(0));

  InitOutputCoorFile(fout, System, argc, argv);

  int mtype = FindMoleculeName(name, System);
  MOLECULETYPE *mt = nullptr;
  if (mtype != -1) {
    mt = &System.MoleculeType[mtype];
  }

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_coor)) {
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
          vec3d gc = GeomCentre(count, beads, System.Bead);
          for (int dd = 0; dd < 3; dd++) {
            System.Bead[System.Molecule[id].Bead[0]].Position.v[dd] = gc.v[dd];
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
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  //}}}

  return 0;
}
