#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "For each step, print the size of the largest aggregate.",

  "Usage: LargestAgg <input> <in.agg> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 10, // number of valid lines OptSpec (not counting last {nullptr})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input structure file", OPT_ARG},
  {"<in.agg>", nullptr, "input agg file", OPT_ARG},
  {"<output>", nullptr, "output file with largest aggregates", OPT_ARG},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  int count = 0;
  // <input> - input structure file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);
  // <in.agg> - input aggregate file
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  // <output> - file largest aggregates
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
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

  // print the initial stuff to output file //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  count = 1;
  fprintf(fw, "# column: (%d) step, ", count++);
  fprintf(fw, "(%d) largest aggregate size; ", count++);
  fprintf(fw, "number of: (%d) molecules, ", count++);
  fprintf(fw, "(%d) aggreages", count++);
  putc('\n', fw);
  fclose(fw); //}}}

  // open <in.agg> and skip the first two lines //{{{
  FILE *fr = OpenFile(input_agg, "r");
  while (getc(fr) != '\n')
    ;
  while (getc(fr) != '\n')
    ; //}}}

  // main loop //{{{
  int count_step = 0,
      count_used = 0,
      agg_lines = 2; // first two lines already read (skipped)
  while (true) { // cycle ends with 'Last Step' line in agg file
    PrintStep(&count_step, commons.start, commons.silent);
    if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
      count_step--;
      break;
    }

    // decide whether this timestep is to be used for averages and distributions
    bool use = false;
    if (UseStep(commons, count_step)) {
      use = true;
    }
    if (use) { //{{{
      count_used++; // just to print at the end
      int largest = 0;
      for (int i = 0; i < Count->Aggregate; i++) {
        if (Aggregate[i].nMolecules > largest) {
          largest = Aggregate[i].nMolecules;
        }
      }
      // print averages to output file
      fw = OpenFile(fout, "a");
      fprintf(fw, "%5d", count_step);
      fprintf(fw, " %7d", largest);
      fprintf(fw, " %7d", Count->Molecule);
      fprintf(fw, " %7d", Count->Aggregate);
      putc('\n', fw);
      fclose(fw);
      ReInitAggregate(System, Aggregate);
    } //}}}

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
    fprintf(stdout, " (%d used for distributions and overall averages)\n",
            count_used);
  } //}}}
  //}}}

  // free memory
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);

  return 0;
}
