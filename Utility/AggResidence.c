#include "../src/AnalysisTools.h"

//TODO: warning if 0 mass leading to nan in output file

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AggResidence calculates distribution of residence times for molecules "
  "in any aggregate of at least given size. Note that 1) the utility does not "
  "care about which aggregate (i.e., no labeling or anything like that) and "
  "2) simulation sampling may significantly skew the results (e.g., "
  "a molecule may come in and out of an aggregate between two saved "
  "configurations.",

  "Usage: AggResidence <in.stru> <in.agg> <agg size> <output> [options]",
  .args = 4, // number of mandatory arguments
  .all = 11, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_VERSION],
  {"<in.stru>", NULL, "input structure file", OPT_ARG},
  {"<in.agg>", NULL, "input agg file", OPT_ARG},
  {"<agg size>", NULL, "minimum aggregate size", OPT_ARG},
  {"<output>", NULL, "output filemane with the distribution", OPT_ARG},
  {NULL},
}; //}}}

// structure for options //{{{
struct OPT {
}; //}}}

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

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // open <in.agg> and skip the first two lines //{{{
  FILE *fr = OpenFile(input_agg, "r");
  SkipLine(fr);
  SkipLine(fr);
  //}}}

  // TODO: make it better - not just a random constant
  int max_residence_time = 1000;
  ArrNDi *residence_distr = CreateArr2Di(Count->Molecule, max_residence_time);
  int *residence_time = calloc(Count->Molecule, sizeof *residence_time);

  // main loop //{{{
  int count_step = 0,
      count_used = 0,
      agg_lines = 2; // first two lines already read (skipped)
  while (true) { // cycle ends with 'Last Step' line in agg file
    PrintStep(&count_step, commons.start, commons.silent);

    // decide whether this timestep is to be used
    bool use = false;
    if (UseStep(commons, count_step)) {
      use = true;
    }
    if (use) { //{{{
      if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
        count_step--;
        break;
      }
      count_used++; // just to print at the end
      for (int i = 0; i < Count->Molecule; i++) { //{{{
        MOLECULE *mol = &System.Molecule[i];
        // printf("Size: Aggregate[mol->Aggregate].nMolecules = %d\n",
        //        Aggregate[mol->Aggregate].nMolecules);
        if (Aggregate[mol->Aggregate].nMolecules >= min_agg_size) {
          residence_time[i]++;
        } else {
          if (residence_time[i] > 0) {
            if (residence_time[i] == max_residence_time) {
              err_msg("TOO LONG IN AN AGGREGATE; USING max_residence_time ...WILL BE REWORKED!");
              PrintWarning();
              AddArr2D(residence_distr, i, max_residence_time - 1, 1);
            } else {
              // printf("OK ... %d %d", residence_time[i],
              //        GetArr2D(residence_distr, i, residence_time[i] - 1));
              AddArr2D(residence_distr, i, residence_time[i] - 1, 1);
              // printf(" ... %d\n",
              //        GetArr2D(residence_distr, i, residence_time[i] - 1));
            }
          }
          residence_time[i] = 0;
        }
      } //}}}
      ReInitAggregate(System, Aggregate);
      //}}}
    } else {
      if (!SkipAggregates(fr, input_agg, &agg_lines)) {
        count_step--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_step == commons.end) {
      break;
    }
  }
  fclose(fr);

  // print last step //{{{
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d", count_step);
    fprintf(stdout, " (%d used)\n", count_used);
  } //}}}
  //}}}

  for (int i = 0; i < Count->Molecule; i++) {
    if (residence_time[i] > 0) {
      if (residence_time[i] == max_residence_time) {
        err_msg("TOO LONG IN AN AGGREGATE; USING max_residence_time ...WILL BE REWORKED!");
        PrintWarning();
        AddArr2D(residence_distr, i, max_residence_time - 1, 1);
      } else {
        // printf("OK ... %d %d", residence_time[i],
        //        GetArr2D(residence_distr, i, residence_time[i] - 1));
        AddArr2D(residence_distr, i, residence_time[i] - 1, 1);
        // printf(" ... %d\n",
        //        GetArr2D(residence_distr, i, residence_time[i] - 1));
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
