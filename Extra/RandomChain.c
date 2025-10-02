#include "../AnalysisTools.h"
#include <stdlib.h>

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "Utility description\n");
  }
  fprintf(ptr, "Usage: %s <length> <output> [options]\n\n", cmd);

  fprintf(ptr, "<length>            chain length\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -r <int>          number of generated chains "
          "(default: 1000)\n");
  fprintf(ptr, "  -alpha <double>   ionization fraction "
          "(default: 0.3)\n");
  fprintf(ptr, "  -s <int>          random number generator seed\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int repeat;        // -r
  double alpha;      // -alpha
  int seed;
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 3, all = common + 3, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "--silent", "--help", "--version", "-r", "-alpha", "-s");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
  // chain length
  long length = 0;
  if (!IsNaturalNumber(argv[++count], &length)) {
    ErrorNaN("<length>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  //}}}

  // options //{{{
  SYS_FILES trash = InitSysFiles; // not used
  opt->c = CommonOptions(argc, argv, trash);
  opt->repeat = 1000;
  OneNumberOption(argc, argv, "-r", &opt->repeat, 'i');
  opt->alpha = 0.3;
  OneNumberOption(argc, argv, "-alpha", &opt->alpha, 'd');
  // reasonable default seed: time & process id
  opt->seed = time(0) * getpid();
  OneNumberOption(argc, argv, "-s", &opt->seed, 'i');
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  bool *chain = calloc(length, sizeof *chain);
  size_t shape2D[2] = {length, 2};
  ArrND distr = NewArrND(2, shape2D);

  pcg32_random_t rng;
  pcg32Seed(&rng, (uint64_t)opt->seed);

  int num = length * opt->alpha;
  for (int i = 0; i < opt->repeat; i++) {
    InitBoolArray(chain, length, false);
    for (int j = 0; j < num; j++) {
      int bead = pcg32Rand0Int(&rng, (uint64_t)length);
      if (!chain[bead]) { // not chosen yet
        chain[bead] = true;
      } else { // was already chosen
        j--;
      }
    }
    // pro-forma check that correct number of beads was picked //{{{
    count = 0;
    for (int j = 0; j < length; j++) {
      if (chain[j]) {
        count++;
      }
    }
    if (count != num) {
      snprintf(ERROR_MSG, LINE, "incorrect number of randomly chosen beads; "
               "%s%d%s instead of %s%d%s (rng seed: %s%u%s)",
               ErrYellow(), count, ErrRed(), ErrYellow(), num, ErrRed(),
               ErrYellow(), opt->seed, ErrRed());
      PrintError();
      exit(1);
    } //}}}
    for (int j = 0; j < length; j++) {
      printf("%d ...\n", chain[j]);
    }
    putchar('\n');
    // sequence lengths
    int seq = 1;
    bool type = chain[0];
    for (int j = 1; j < length; j++) {
      if (chain[j] == chain[j-1]) {
        seq++;
      } else {
        distr.data[idx2d(distr, seq - 1, type)]++;
        seq = 1;
        type = !type;
      }
    }
    distr.data[idx2d(distr, seq - 1, type)]++;
    for (int i = 0; i < length; i++) {
      int count[2] = {(int)distr.data[idx2d(distr, i, 0)],
                      (int)distr.data[idx2d(distr, i, 1)]};
      if (count[0] != 0 || count[1] != 0) {
        printf("%d: %d ... %d\n", i + 1, count[0], count[1]);
      }
    }
  }
  // TODO: ...it should correctly capture numbers of sequences, so now normalize
  //       and print distributions to file

  // free memory - to make valgrind happy //{{{
  free(opt);
  free(chain);
  FreeArrND(distr);
  //}}}

  return 0;
}
