#include "../src/AnalysisTools.h"

static void ReduceSystem(const SYSTEM System, SYSTEM *Sys, const bool write[],
                         int **b_full_to_red, const OPT opt,
                         const COMMON_OPT commons, const FILE_TYPE fout,
                         const int argc, char *argv[]);
static void ReducedWrite(SYSTEM *Sys, const SYSTEM System,
                         const int b_full_to_red[],
                         const bool *write, bool **write2);
static void SaveReduced(SYSTEM *Sys, const SYSTEM System, const int count_saved,
                        const int count_coor, int **b_full_to_red,
                        const bool write[], const OPT opt,
                        const COMMON_OPT commons, const FILE_TYPE fout,
                        const int argc, char *argv[]);
static void ScaleCoordinates(SYSTEM *System, const double scale);
static void MoveCoordinates(SYSTEM *System, const vec3d move);
static void CopyWrite(const int num, bool *new, bool *old);
static void ConstrainCoordinates(SYSTEM *System, const OPT opt,
                                 bool *write_new, bool **write_orig);
static void TransformAndSave(SYSTEM *System, const OPT opt, const FILE_TYPE f,
                             bool *write, const int count_coor,
                             const int argc, char *argv[]);
static void SaveTimestep(SYSTEM *Sys, SYSTEM *System,
                         const int count_saved, const int count_coor,
                         int **b_full_to_red, bool *write, const OPT opt,
                         const COMMON_OPT commons,
                         const FILE_TYPE fout, const int argc, char *argv[]);
static bool ReadTimestepSilent(const SYS_FILES in, FILE *fr,
                               SYSTEM *System, int *line_count);

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Selected creates new coordinate file in the extension-specified format "
  "that contains specified bead and/or molecule types. Periodic boundary "
  "conditions can be either stripped away or applied (which happens first if "
  "both '--join' and '--wrap' options are used).",

  "Usage: Selected <input> <output> [options]",
  .args = 2, // number of mandatory arguments
  .all = 26, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<output>", NULL, "output coordinate file", OPT_ARG},
  {"-bt", "<bead type>", "bead types to exclude", OPT_EXTRA},
  {"-mt", "<mol type>", "molecule types to exclude", OPT_EXTRA},
  {"--keep", NULL, "save only the specified types instead of excluding them", OPT_EXTRA},
  {"--join", NULL, "join molecules (remove pbc)", OPT_EXTRA},
  {"--wrap", NULL, "wrap coordinates (i.e., apply pbc)", OPT_EXTRA},
  {"-n", "<int(s)>", "save only specified timesteps (--last overrides this option)", OPT_EXTRA},
  {"--last", NULL, "use only the last step (-st/-e/-n options are ignored)", OPT_EXTRA},
  {"-sc", "<float>", "divide all coordinates by given value", OPT_EXTRA},
  {"-m", "3x<float>", "move all coordinates by given vector (-sc option is applied first)", OPT_EXTRA},
  {"-cx", "2x<float>", "constrain x-coordinate to specified dimensions (in fraction of output box); multiple pairs possible", OPT_EXTRA},
  {"-cy", "2x<float>", "constrain y-coordinate to specified dimensions (in fraction of output box); multiple pairs possible", OPT_EXTRA},
  {"-cz", "2x<float>", "constrain z-coordinate to specified dimensions (in fraction of output box); multiple pairs possible", OPT_EXTRA},
  {"--real", NULL, "use real coordinates for -cx/-cy/-cz options instead of box fractions", OPT_EXTRA},
  {"--reduce", NULL, "reduce the structure to contaion only beads in the coordinate file", OPT_EXTRA},
  {"-b", "3x<float>", "set box size for all timesteps", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool bt, mt,               // -bt/-mt
       keep,                 // --keep
       join, wrap, last,     // --join --wrap --last
       real,                 // --real
       reduce;               // --reduce
  int n_save[100], n_number; // -n
  double scale;              // -sc
  vec3d move,                // -m
        box;                 // -b
  double ca[3][100];         // -cx/y/z ... slice(s)' coordinates
  int ca_count[3];           // -cx/y/z ... number of slices per axis
}; //}}}

// reduce system based on the beds present in the timestep //{{{
static void ReduceSystem(const SYSTEM System, SYSTEM *Sys, const bool *write,
                         int **b_full_to_red, const OPT opt,
                         const COMMON_OPT commons, const FILE_TYPE fout,
                         const int argc, char *argv[]) {
  *Sys = CopySystem(System);
  // fill Sys.BeadCoor for beads to be saved
  int count = -1;
  for (int i = 0; i < Sys->Count.BeadCoor; i++) {
    int id = Sys->BeadCoor[i];
    if (write[id]) {
      Sys->BeadCoor[++count] = id;
    } else {
      Sys->Bead[id].InTimestep = false;
    }
  }
  Sys->Count.BeadCoor = count;
  // prune the system and generate full->reduced bead ids transformation
  *b_full_to_red = calloc(System.Count.Bead, sizeof *b_full_to_red);
  if (!b_full_to_red) {
    ErrorAlloc("b_full_to_red");
  }
  PruneSystem2(Sys, *b_full_to_red);
  // print the system to save if required
  if (commons.verbose) {
    fprintf(stdout, "\n################\n");
    fprintf(stdout, "# Saved System #\n");
    fprintf(stdout, "################\n");
    VerboseOutput(*Sys);
  }
  // if output is vtf, write the structure part
  if (fout.type == VTF_FILE) {
    WriteStructure(fout, *Sys, -1, false, argc, argv);
  // if output is vcf, write a new structure file
  } else if (fout.type == VCF_FILE) {
    FILE_TYPE f_vsf;
    s_strcpy(f_vsf.name, fout.name, LINE);
    f_vsf.name[strlen(f_vsf.name)-2] = 's';
    f_vsf.type = StructureFileType(f_vsf.name);
    WriteStructure(f_vsf, *Sys, -1, false, argc, argv);
  }
} //}}}
// populate reduced BeadCoor for a given timestep //{{{
static void ReducedWrite(SYSTEM *Sys, const SYSTEM System,
                         const int b_full_to_red[],
                         const bool *write, bool **write2) {
  Sys->Count.BeadCoor = 0;
  *write2 = calloc(Sys->Count.Bead, sizeof *write2);
  if (!write2) {
    ErrorAlloc("write2");
  }
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id_orig = System.BeadCoor[i];
    int id_new = b_full_to_red[id_orig];
    if (id_new == -1) {
      continue;
    }
    for (int dd = 0; dd < 3; dd++) {
      Sys->Bead[id_new].Position.v[dd] = System.Bead[id_orig].Position.v[dd];
    }
    Sys->BeadCoor[Sys->Count.BeadCoor] = id_new;
    Sys->Box = System.Box;
    (*write2)[id_new] = write[id_orig];
    Sys->Count.BeadCoor++;
  }
} //}}}
// save reduced coordinate (and structure part for the first saved step) //{{{
static void SaveReduced(SYSTEM *Sys, const SYSTEM System, const int count_saved,
                        const int count_coor, int **b_full_to_red,
                        const bool write[], const OPT opt,
                        const COMMON_OPT commons, const FILE_TYPE fout,
                        const int argc, char *argv[]) {
  if (count_saved == 1) {
    ReduceSystem(System, Sys, write, b_full_to_red, opt, commons,
                 fout, argc, argv);
  }
  bool *write2 = NULL;
  ReducedWrite(Sys, System, *b_full_to_red, write, &write2);
  TransformAndSave(Sys, opt, fout, write2, count_coor, argc, argv);
  free(write2);
}
//}}}
static void ScaleCoordinates(SYSTEM *System, const double scale) { //{{{
  if (scale != 1) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position.v[dd] /= scale;
      }
    }
    for (int dd = 0; dd < 3; dd++) {
      System->Box.Length.v[dd] /= scale;
    }
    CalculateBoxData(&System->Box, 0);
  }
} //}}}
static void MoveCoordinates(SYSTEM *System, const vec3d move) { //{{{
  if (move.x != 0 || move.y != 0 || move.z != 0) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[id].Position.v[dd] += move.v[dd];
      }
    }
  }
} //}}}
static void CopyWrite(const int num, bool *new, bool *old) { //{{{
  for (int i = 0; i < num; i++) {
    new[i] = old[i];
  }
} //}}}
// ConstrainCoordinates() //{{{
static void ConstrainCoordinates(SYSTEM *System, const OPT opt,
                                 bool *write_new, bool **write_orig) {
  if (opt.ca_count[0] == 0 && opt.ca_count[1] == 0 && opt.ca_count[2] == 0) {
    return;
  }
  // recalculate constraints if --real was not used
  double con[3][100];
  bool init[3] = {true, true, true};
  for (int dd = 0; dd < 3; dd++) {
    for (int i = 0; i < opt.ca_count[dd]; i++) {
      con[dd][i] = opt.ca[dd][i];
      init[dd] = false;
      if (!opt.real) {
        con[dd][i] *= System->Box.Length.v[dd];
      }
    }
  }
  COUNT *Count = &System->Count;
  // allocate array to store the original 'write'
  *write_orig = malloc(Count->Bead * sizeof *write_orig);
  if (!write_orig) {
    ErrorAlloc("write_orig");
  }
  // save the original 'write' array
  CopyWrite(Count->Bead, *write_orig, write_new);
  // go over beads in the coordinate file
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System->BeadCoor[i];
    // skip constraints for beads not to be written
    if (!(*write_orig)[id]) {
      continue;
    }
    // check -cx/-cy/-cz constraint
    vec3d *pos = &System->Bead[id].Position;
    bool save[3] = {init[0], init[1], init[2]};
    for (int dd = 0; dd < 3; dd++) {
      for (int j = 0; j < opt.ca_count[dd]; j+=2) {
        if (pos->v[dd] >= con[dd][j] && pos->v[dd] <= con[dd][j+1]) {
          save[dd] = true;
          break;
        }
      }
    }
    // if at least one axis constraint isn't met, don't save the bead
    if (!save[0] || !save[1] || !save[2]) {
      write_new[id] = false;
    }
  }
} //}}}
// make all the alterations and saves into one function //{{{
static void TransformAndSave(SYSTEM *System, const OPT opt, const FILE_TYPE f,
                             bool *write, const int count_coor,
                             const int argc, char *argv[]) {
  bool *write2 = NULL; // used in case of -cx/-cy/-cz constraints
  ConstrainCoordinates(System, opt, write, &write2);
  ScaleCoordinates(System, opt.scale);
  MoveCoordinates(System, opt.move);
  WrapJoinCoordinates(System, opt.wrap, opt.join);
  WriteTimestep(f, *System, count_coor, write, argc, argv);
  if (opt.ca_count[0] > 0 ||
      opt.ca_count[1] > 0 ||
      opt.ca_count[2] > 0) {
    // restore the original 'write' array
    CopyWrite(System->Count.Bead, write, write2);
    // free the array alloc'd in ConstrainCoordinates()
    free(write2);
  }
} //}}}
// save timestep, including --reduce and transformation options and all //{{{
static void SaveTimestep(SYSTEM *Sys, SYSTEM *System,
                         const int count_saved, const int count_coor,
                         int **b_full_to_red, bool *write, const OPT opt,
                         const COMMON_OPT commons,
                         const FILE_TYPE fout, const int argc, char *argv[]) {
  if (opt.reduce) {
    SaveReduced(Sys, *System, count_saved, count_coor, b_full_to_red,
                write, opt, commons, fout, argc, argv);
  } else {
    TransformAndSave(System, opt, fout, write, count_coor, argc, argv);
  }
} //}}}

// ReadTimestep wrapper that suppresses all stderr output //{{{
static bool ReadTimestepSilent(const SYS_FILES in, FILE *fr,
                               SYSTEM *System, int *line_count) {
  FILE *devnull = fopen("/dev/null", "w");
  int saved_fd = dup(STDERR_FILENO);
  dup2(fileno(devnull), STDERR_FILENO);
  fclose(devnull);
  bool ok = ReadTimestep(in, fr, System, line_count);
  fflush(stderr);
  dup2(saved_fd, STDERR_FILENO);
  close(saved_fd);
  return ok;
} //}}}
// check if the current words/split[] is the start of a timestep //{{{
static bool IsTimestepStartLine(int type, const SYSTEM *System) {
  if (type == XYZ_FILE) {
    // first line of an XYZ timestep is just the bead count
    long val;
    return words == 1 && IsNaturalNumber(split[0], &val) &&
           val > 0 && val <= System->Count.Bead;
  } else if (type == VTF_FILE || type == VCF_FILE) {
    // matches VtfCheckTimestepLine() logic
    return (words == 1 && split[0][0] == 't') ||
           (words > 1 && split[0][0] == 't' && split[1][0] == 'o') ||
           split[0][0] == 'o';
  } else if (type == LTRJ_FILE) {
    return words >= 2 && strcmp(split[0], "ITEM:") == 0 &&
           strcmp(split[1], "TIMESTEP") == 0;
  }
  return false;
} //}}}
// scan backward from end_pos to find the start of the preceding timestep //{{{
static bool FindPrevTimestepStart(FILE *fr, int type, const SYSTEM *System,
                                  long end_pos, long *ts_start) {
  long pos = end_pos - 1;
  // skip any trailing newlines before end_pos
  while (pos >= 0) {
    fseek(fr, pos, SEEK_SET);
    int c = fgetc(fr);
    if (c != '\n' && c != '\r') break;
    pos--;
  }
  while (pos >= 0) {
    // find the start of the line containing pos
    long line_start = pos;
    while (line_start > 0) {
      fseek(fr, line_start - 1, SEEK_SET);
      if (fgetc(fr) == '\n') break;
      line_start--;
    }
    // read the line and test if it is a timestep start
    fseek(fr, line_start, SEEK_SET);
    ReadAndSplitLine(fr, SPL_STR, " \t\n");
    if (IsTimestepStartLine(type, System)) {
      *ts_start = line_start;
      return true;
    }
    // advance backward past this line
    pos = line_start - 1;
    while (pos >= 0) {
      fseek(fr, pos, SEEK_SET);
      int c = fgetc(fr);
      if (c != '\n' && c != '\r') break;
      pos--;
    }
  }
  return false;
} //}}}

// structure for the callback function //{{{
struct user_data {
  OPT opt;
  COMMON_OPT commons;
  SYSTEM *Sys;
  int **b_full_to_red;
  bool *write;
  FILE_TYPE fout;
  int count_saved;
  int argc;
  char **argv;
}; //}}}
// adaptor for the SaveTimestep() function //{{{
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data *)userdata;
  // handle -n option: only save specified timesteps
  if (p->opt.n_number != -1) {
    bool in_list = false;
    for (int i = 0; i < p->opt.n_number; i++) {
      if (p->opt.n_save[i] == step->coor) {
        in_list = true;
        break;
      }
    }
    if (!in_list) {
      return;
    }
  }
  // warn and skip if lammps data file already has one saved timestep
  if (p->fout.type == LDATA_FILE && p->count_saved == 1) {
    err_msg("only one timestep can be saved to lammps data file");
    PrintWarnFile(p->fout.name, "\0", "\0");
    return;
  }
  // apply -b override
  if (p->opt.box.x != -1) {
    for (int dd = 0; dd < 3; dd++) {
      System->Box.Length.v[dd] = p->opt.box.v[dd];
    }
    CalculateBoxData(&System->Box, 0);
  }
  p->count_saved++;
  SaveTimestep(p->Sys, System, p->count_saved, step->coor, p->b_full_to_red,
               p->write, p->opt, p->commons, p->fout, p->argc, p->argv);
}; //}}}

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
  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }
  opt.keep = BoolOption(argc, argv, "--keep");
  opt.join = BoolOption(argc, argv, "--join");
  opt.wrap = BoolOption(argc, argv, "--wrap");
  opt.last = BoolOption(argc, argv, "--last");
  if (!OneNumberOption(argc, argv, "-sc", &opt.scale, 'd')) {
    opt.scale = 1;
  }
  opt.real = BoolOption(argc, argv, "--real");
  if (!ThreeNumbersOption(argc, argv, "-m", opt.move.v, 'd')) {
    for (int dd = 0; dd < 3; dd++) {
      opt.move.v[dd] = 0;
    }
  }
  opt.reduce = BoolOption(argc, argv, "--reduce");
  // constraints (-cx/-cy/-cz) //{{{
  for (int dd = 0; dd < 3; dd++) {
    char option[10];
    if (dd == 0) {
      s_strcpy(option, "-cx", 10);
    } else if (dd == 1) {
      s_strcpy(option, "-cy", 10);
    } else {
      s_strcpy(option, "-cz", 10);
    }
    if (!NumbersOption(argc, argv, 100, option,
                       &opt.ca_count[dd], opt.ca[dd], 'd')) {
      opt.ca_count[dd] = 0;
    } else if ((opt.ca_count[dd] % 2) != 0) { // not even number of numbers
      goto err_constraint;
    }
    for (int i = 0; i < opt.ca_count[dd]; i+=2 ) {
      if (fabs(opt.ca[0][i] - opt.ca[dd][i+1]) < 0.0001) { // same numbers in a pair
        goto err_constraint;
      } else if (!opt.real && (opt.ca[dd][i] < 0 || opt.ca[dd][i] > 1 ||
                                opt.ca[dd][i+1] < 0 || opt.ca[dd][i+1] > 1)) {
        goto err_constraint;
      } else if (opt.ca[dd][i] > opt.ca[dd][i+1]) { // switch so [i] < [i+1]
        SwapDouble(&opt.ca[dd][i], &opt.ca[dd][i+1]);
      }
    }
  }
  // // -cy
  // if (!NumbersOption(argc, argv, 100, "-cy", &opt.c_count[1], opt.cy, 'd')) {
  //   opt.c_count[1] = 0;
  // } else if ((opt.c_count[1] % 2) != 0) { // not even number of numbers
  //   goto err_constraint;
  // }
  // for (int i = 0; i < opt.c_count[1]; i+=2 ) { // same numbers in a pair
  //   if (fabs(opt.cy[i] - opt.cy[i+1]) < 0.0001) {
  //     goto err_constraint;
  //   } else if (opt.cy[i] > opt.cy[i+1]) { // switch so [i] < [i+1]
  //     SwapDouble(&opt.cy[i], &opt.cy[i+1]);
  //   }
  // }
  // // -cz
  // if (!NumbersOption(argc, argv, 100, "-cz", &opt.c_count[2], opt.cz, 'd')) {
  //   opt.c_count[2] = 0;
  // } else if ((opt.c_count[2] % 2) != 0) { // not even number of numbers
  //   goto err_constraint;
  // }
  // for (int i = 0; i < opt.c_count[2]; i+=2 ) { // same numbers in a pair
  //   if (fabs(opt.cz[i] - opt.cz[i+1]) < 0.0001) {
  //     goto err_constraint;
  //   } else if (opt.cz[i] > opt.cz[i+1]) { // switch so [i] < [i+1]
  //     SwapDouble(&opt.cz[i], &opt.cz[i+1]);
  //   }
  // }
  // error - triggers only via goto command
  if (false) {
    err_constraint:
      err_msg("requires coordinate pairs (with two distinct numbers per pair) "
              "in units of box fraction, i.e., <0,1> (unless --real is used)");
      PrintErrorOption("-cx/-cy/-cz");
      exit(1);
  } //}}}
  if (!ThreeNumbersOption(argc, argv, "-b", opt.box.v, 'd')) {
    for (int dd = 0; dd < 3; dd++) {
      opt.box.v[dd] = -1;
    }
  } //}}}

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  if (opt.box.v[0] != -1) {
    System.Box.Length = opt.box;
    System.Box.alpha = 90;
    System.Box.beta = 90;
    System.Box.gamma = 90;
    if (!CalculateBoxData(&System.Box, 0)) {
      err_msg("CalculateBoxData() function - should not happen!");
      PrintError();
      exit(1);
    }
  }
  if (!opt.real) {
    for (int dd = 0; dd < 3; dd++) {
      opt.move.v[dd] *= System.Box.Length.v[dd];
    }
  }

  if (opt.join && Count->Molecule == 0) {
    err_msg("no molecules to join");
    PrintWarning();
  }

  // specify beads to save (possibly using -bt and/or -mt options) //{{{
  /*
   * reverse=true ... save only the specified species
   * reverse=false ... exclude the specified species
   *
   * First setting all bead/molecule types to !reverse and then adjusting this
   * if -bt/-mt options are present correctly specifies which
   * bead/molecule types to save
   */
  // auxiliary arrays holding which bead/molecule types to save
  bool *write_bt = malloc(Count->BeadType * sizeof *write_bt),
       *write_mt = malloc(Count->MoleculeType * sizeof *write_mt);
  if (!write_bt || !write_mt) {
    ErrorAlloc("write_bt/write_mt");
  }
  // first assume all bead types are saved/excluded based on --keep option...
  InitBoolArray(write_bt, Count->BeadType, !opt.keep);
  // ... then adjust if -bt option is present
  opt.bt = TypeOption(argc, argv, "-bt", 'b', opt.keep, write_bt, System);
  // first assume all molecule types are saved/excluded...
  InitBoolArray(write_mt, Count->MoleculeType, !opt.keep);
  // ... then adjust if -mt option is present
  opt.mt = TypeOption(argc, argv, "-mt", 'm', opt.keep, write_mt, System);
  // array for holding which beads to save/exclude
  bool *write = malloc(Count->Bead * sizeof *write);
  if (!write) {
    ErrorAlloc("write");
  }
  // first assume all are saved/excluded...
  InitBoolArray(write, Count->Bead, !opt.keep);
  // then check possible -mt option...
  if (opt.mt) {
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < mt->Number; j++) {
        int mol = mt->Index[j];
        for (int k = 0; k < mt->nBeads; k++) {
          int id = System.Molecule[mol].Bead[k];
          if (write_mt[i] == opt.keep) {
            write[id] = opt.keep; // save/exclude based on --keep
          }
        }
      }
    }
  }
  // ... and possible -bt option
  if (opt.bt) {
    for (int i = 0; i < Count->Bead; i++) {
      int type = System.Bead[i].Type;
      if (write_bt[type] == opt.keep) {
        write[i] = opt.keep;
      }
    }
  }
  // free the auxiliary arrays
  free(write_mt);
  free(write_bt); //}}}

  // '-n' option - specify timestep ids //{{{
  opt.n_number = -1;
  InitIntArray(opt.n_save, 100, 0);
  NumbersOption(argc, argv, 100, "-n", &opt.n_number, opt.n_save, 'i');
  // ignore -st/-e/-sk when -n is used
  if (opt.n_number != -1) {
    commons.start = 1;
    commons.end = -1;
    commons.skip = 1;
  }
  // SortArrayInt(opt.n_save, opt.n_number, 0);
  SortArray(opt.n_save, opt.n_number, 0, 'i'); //}}}

  if (commons.verbose) {
    if (opt.reduce) {
      fprintf(stdout, "\n##################\n");
      fprintf(stdout, "# Initial System #\n");
      fprintf(stdout, "##################\n");
    }
    VerboseOutput(System);
  }

  InitOutputCoorFile(fout, System, argc, argv);

  // helper variables for --reduce option
  SYSTEM Sys; // the reduced system
  int *b_full_to_red = NULL; // full-system bead ids to reduced-system ids

  if (opt.last) { // read from end of file to find and save the last step //{{{
    FILE *fr = OpenFile(in.coor.name, "r");
    fseek(fr, 0, SEEK_END);
    long pos = ftell(fr); // start scanning from the end
    bool found = false;
    long ts_start;
    while (FindPrevTimestepStart(fr, in.coor.type, &System, pos, &ts_start)) {
      int line_count = 0;
      fseek(fr, ts_start, SEEK_SET);
      if (ReadTimestepSilent(in, fr, &System, &line_count)) {
        SaveTimestep(&Sys, &System, 1, 0, &b_full_to_red,
                     write, opt, commons, fout, argc, argv);
        if (!commons.silent) {
          fprintf(stdout, "Saved last step\n");
          fflush(stdout);
        }
        found = true;
        break;
      } else {
        snprintf(ERROR_MSG, LINE, "disregarding corrupt last step");
        PrintWarnFile(in.coor.name, "\0", "\0");
        pos = ts_start; // try the preceding timestep
      }
    }
    if (!found) {
      remove(fout.name);
      err_msg("no valid timestep found");
      PrintErrorFile(in.coor.name, "\0", "\0");
    }
    fclose(fr);
    //}}}
  } else { // normal: use MainLoopCoor() //{{{
    STEP step = InitStep;
    struct user_data ud = { opt, commons, &Sys, &b_full_to_red,
                            write, fout, 0, argc, argv };
    MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);
    if (step.coor == 0) { // error - input file without a valid timestep
      remove(fout.name);
      err_msg("no valid timestep found");
      PrintErrorFile(in.coor.name, "\0", "\0");
    } else if (commons.start > step.coor) { // warn if no timesteps were written
      remove(fout.name);
      err_msg("no coordinates written (starting timestep is higher "
              "than the total number of timesteps)");
      PrintWarning();
    }
  } //}}}

  // free memory
  FreeSystem(&System);
  free(write);
  if (opt.reduce) {
    FreeSystem(&Sys);
    free(b_full_to_red);
  }

  return 0;
}
