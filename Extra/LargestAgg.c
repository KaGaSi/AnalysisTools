#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
For each step, print the size of the largest aggregate.\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <in.agg> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input structure file\n");
  fprintf(ptr, "<in.agg>            input agg file\n");
  fprintf(ptr, "<output>            output file with largest aggregates\n");
  fprintf(ptr, "[options]\n");
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

  // define options & check their validity
  int common = 7, all = common + 0, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "--verbose", "--silent", "--help",
               "--version");

  // commad line arguments before reading the structure //{{{
  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
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
  //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // print the initial stuff to output file //{{{
  PrintByline(fout, argc, argv);
  FILE *fw = OpenFile(fout, "a");
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
    PrintStep(&count_step, opt->c.start, opt->c.silent);
    if (ReadAggregates(fr, input_agg, &System, Aggregate, &agg_lines) < 0) {
      count_step--;
      break;
    }

    // decide whether this timestep is to be used for averages and distributions
    bool use = false;
    if (UseStep(opt->c, count_step)) {
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
    if (count_step == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  // print last step //{{{
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d", count_step);
    fprintf(stdout, " (%d used for distributions and overall averages)\n",
            count_used);
  } //}}}
  //}}}

  // free memory - to make valgrind happy
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  free(opt);

  return 0;
}
