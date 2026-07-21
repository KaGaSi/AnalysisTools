#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AggResidence calculates distribution of residence times for molecules "
  "in any aggregate of at least given size. Note that 1) the utility does not "
  "care about which aggregate (i.e., no labeling or anything like that) and "
  "2) simulation sampling may significantly skew the results (e.g., "
  "a molecule may have come in and out of an aggregate between two saved "
  "configurations.",

  "Usage: AggResidence <in.stru> <in.agg> <agg size> <output> [options]",
  .args = 4, // number of mandatory arguments
  .all = 11, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_VERSION],
  {"<in.stru>", nullptr, "input structure file", OPT_ARG},
  {"<in.agg>", nullptr, "input agg file", OPT_ARG},
  {"<agg size>", nullptr, "minimum aggregate size", OPT_ARG},
  {"<output>", nullptr, "output filemane with the distribution", OPT_ARG},
  {nullptr},
}; //}}}

// structure for options //{{{
struct OPT {
}; //}}}

// all state shared between main() and the per-timestep callback //{{{
struct user_data {
  AGGREGATE *Aggregate;
  long min_agg_size;
  int max_residence_time,
      *residence_time;
  ArrNDi *residence_distr;
}; //}}}

// per-timestep residence time bookkeeping //{{{
static void Calculation(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *ud = userdata;
  AGGREGATE *Aggregate = ud->Aggregate;
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Molecule; i++) { //{{{
    MOLECULE *mol = &System->Molecule[i];
    if (Aggregate[mol->Aggregate].nMolecules >= ud->min_agg_size) {
      ud->residence_time[i]++;
    } else {
      if (ud->residence_time[i] > 0) {
        if (ud->residence_time[i] == ud->max_residence_time) {
          // TODO: huh? ...I guess better define maximum time (see below)?
          err_msg("TOO LONG IN AN AGGREGATE; USING max_residence_time "
                  "...WILL BE REWORKED!");
          PrintWarning();
          AddArr2D(ud->residence_distr, i, ud->max_residence_time - 1, 1);
        } else {
          AddArr2D(ud->residence_distr, i, ud->residence_time[i] - 1, 1);
        }
      }
      ud->residence_time[i] = 0;
    }
  } //}}}
  ReInitAggregate(*System, Aggregate);
} //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  // OPT opt;
  int count = 0;
  // <input> - input structure file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);
  // <in.agg> - input aggregate file
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  // <agg size> - minimum size
  long min_agg_size;
  if (!IsWholeNumber(argv[++count], &min_agg_size)) {
    ErrorNaN("<agg size>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  // <output> - distribution file
  char output[LINE];
  s_strcpy(output, argv[++count], LINE);

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AGGREGATE *Aggregate = nullptr;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // TODO: make it better - not just a random constant
  int max_residence_time = 1000;
  ArrNDi *residence_distr = CreateArr2Di(Count->Molecule, max_residence_time);
  int *residence_time = calloc(Count->Molecule, sizeof *residence_time);
  if (!residence_time) {
    ErrorAlloc("residence_time");
  }

  // main loop //{{{
  struct user_data ud = {
    .Aggregate = Aggregate,
    .min_agg_size = min_agg_size,
    .max_residence_time = max_residence_time,
    .residence_time = residence_time,
    .residence_distr = residence_distr,
  };
  STEP step = InitStep;
  MainLoopAgg(&System, input_agg, commons, &step, Aggregate, Calculation, &ud);
  int count_used = step.used; //}}}

  for (int i = 0; i < Count->Molecule; i++) {
    if (residence_time[i] > 0) {
      if (residence_time[i] == max_residence_time) {
        // TODO: huh? ...I guess 'better' maximum time?
        err_msg("TOO LONG IN AN AGGREGATE; USING max_residence_time "
            "...WILL BE REWORKED!");
        PrintWarning();
        AddArr2D(residence_distr, i, max_residence_time - 1, 1);
      } else {
        AddArr2D(residence_distr, i, residence_time[i] - 1, 1);
      }
    }
  }

  // per molecule normalization factors + overall normalization factor
  double *distr_norm = calloc(Count->Molecule + 1, sizeof *distr_norm);
  for (int i = 0; i < Count->Molecule; i++) {
    for (int j = 0; j < count_used; j++) {
      distr_norm[i] += GetArr2D(residence_distr, i, j);
      distr_norm[Count->Molecule] += GetArr2D(residence_distr, i, j);
    }
  }

  // sum up distributions from all nMolecules
  // count_used is the maximum theoretical residence time
  double *distr = calloc(count_used, sizeof *distr);
  for (int i = 0; i < count_used; i++) {
    for (int j = 0; j < Count->Molecule; j++) {
      distr[i] += GetArr2D(residence_distr, j, i);
    }
  }

  FILE *fw = PrintBylineOpenFile(output, argc, argv);
  count = 1;
  fprintf(fw, "# Column: ");
  fprintf(fw, "(%d) residence time, ", count++);
  fprintf(fw, "(%d) overall distribution, ", count++);
  fprintf(fw, "(%d+) distribution for individual molecules", count++);
  putc('\n', fw);
  int print_mols = 100;
  int ncols = 2 + print_mols;
  int nrows = count_used;
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    ErrorAlloc("data");
  }
  for (int i = 0; i < nrows; i++) {
    count = 0;
    SetArr2D(data, i, count++, i + 1);
    SetArr2D(data, i, count++, distr[i] / distr_norm[Count->Molecule]);
    for (int j = 0; j < print_mols; j++) {
      double val = GetArr2D(residence_distr, j, i);
      if (distr_norm[j] == 0 && val != 0) {
        err_msg("Huh? Something wrong!");
        PrintError();
      }
      if (distr_norm[j] > 0) {
        val /= distr_norm[j];
      }
      SetArr2D(data, i, count++, val);
    }
  }
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data);
  fclose(fw);

  // free memory
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeArrND(residence_distr);
  free(distr);

  return 0;
}
