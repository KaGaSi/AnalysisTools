#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "RandnomChain generates linear chains of given length with randomly "
  "distributed charged and neutral beads according to given dissociation "
  "<alpha> (0 to 1), also printing distribution of length of charged segments",

  "Usage: RandomChain <length> <out.txt> <out.xyz> [options]",
  .args = 3, // number of mandatory arguments
  .all = 14, // number of valid lines OptSpec (not counting last {nullptr})
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
  {"<length>", nullptr, "chain length", OPT_ARG},
  {"<out.txt>", nullptr, "output distribution file", OPT_ARG},
  {"<out.xyz>", nullptr, "output coordinate (-<n>.xyz ending)", OPT_ARG},
  {"-r", "<int>", "number of generated chains (default: 1000)", OPT_EXTRA},
  {"-alpha", "<double>", "ionization fraction (default: 0.3)", OPT_EXTRA},
  {"-s", "<int>", "random number generator seed", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int repeat;        // -r
  double alpha;      // -alpha
  int seed;
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // chain length
  long length = 0;
  if (!IsNaturalNumber(argv[++count], &length)) {
    ErrorNaN("<length>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  // <out.txt> - output distribution
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <out.xyz> - output xyz file
  char coor_file[LINE] = "";
  s_strcpy(coor_file, argv[++count], LINE);
  //}}}

  // options //{{{
  SYS_FILES trash = InitSysFiles; // not used
  COMMON_OPT commons = CommonOptions(argc, argv, trash);
  opt.repeat = 1000;
  OneNumberOption(argc, argv, "-r", &opt.repeat, 'i');
  opt.alpha = 0.3;
  OneNumberOption(argc, argv, "-alpha", &opt.alpha, 'd');
  // reasonable default seed: time & process id
  opt.seed = time(0) * getpid();
  OneNumberOption(argc, argv, "-s", &opt.seed, 'i');
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  bool *chain = calloc(length, sizeof *chain);
  ArrNDd *distr = CreateArr2Dd(length, 2);
  // size_t shape2D[2] = {length, 2};
  // ArrND distr = NewArrND(2, shape2D);

  pcg32_random_t rng;
  pcg32Seed(&rng, (uint64_t)opt.seed);

  int num[2];
  num[1] = length * opt.alpha;
  num[0] = length - num[1];
  for (int i = 0; i < opt.repeat; i++) {
    InitBoolArray(chain, length, false);
    for (int j = 0; j < num[1]; j++) {
      int bead = pcg32Rand0Int(&rng, (uint64_t)length);
      if (!chain[bead]) { // not chosen yet
        chain[bead] = true;
      } else { // was already chosen, rerun for this 'j'
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
    if (count != num[1]) {
      snprintf(ERROR_MSG, LINE, "incorrect number of randomly chosen beads; "
               "%s%d%s instead of %s%d%s (rng seed: %s%u%s)",
               ErrYellow(), count, ErrRed(), ErrYellow(), num[1], ErrRed(),
               ErrYellow(), opt.seed, ErrRed());
      PrintError();
      exit(1);
    } //}}}
    int seq = 1; // sequence length
    bool type = chain[0]; // true for one type, false for the other
    for (int j = 1; j < length; j++) {
      // if 'j' and 'j-1' beads are the same, continue sequence...
      if (chain[j] == chain[j-1]) {
        seq++;
      // ...if not, save sequence length and restart counting
      } else {
        AddArr2D(distr, seq - 1, type, 1);
        // distr.d[idx2d(distr, seq - 1, type)]++;
        // restart sequence with 'j' as the first bead of the sequence
        seq = 1;
        // change type of bead to the other one
        type = !type;
      }
    }
    AddArr2D(distr, seq - 1, type, 1);
    // distr.d[idx2d(distr, seq - 1, type)]++;
  }

  // normalilzation //{{{
  int norm[2] = {0, 0}; // normalization factor
  for (int i = 0; i < length; i++) {
    norm[0] += GetArr2D(distr, i, 0);
    norm[1] += GetArr2D(distr, i, 1);
    // norm[0] += distr.d[idx2d(distr, i, 0)];
    // norm[1] += distr.d[idx2d(distr, i, 1)];
  }
  // normalize
  for (int i = 0; i < length; i++) {
    SetArr2D(distr, i, 0, GetArr2D(distr, i, 0) / norm[0]);
    SetArr2D(distr, i, 1, GetArr2D(distr, i, 1) / norm[1]);
    // distr.d[idx2d(distr, i, 0)] /= norm[0];
    // distr.d[idx2d(distr, i, 1)] /= norm[1];
  } //}}}

  // calculated integrated distributions //{{{
  ArrNDd *integ = CreateArr2Dd(length, 2);
  SetArr2D(integ, 0, 0, GetArr2D(distr, 0, 0));
  SetArr2D(integ, 0, 1, GetArr2D(distr, 0, 1));
  // ArrND integ = NewArrND(2, shape2D);
  // integ.d[idx2d(integ, 0, 0)] = distr.d[idx2d(integ, 0, 0)];
  // integ.d[idx2d(integ, 0, 1)] = distr.d[idx2d(integ, 0, 1)];
  for (int i = 1; i < length; i++) {
    for (int dd = 0; dd < 2; dd++) {
      double val = GetArr2D(integ, i - 1, dd) + GetArr2D(distr, i, dd);
      SetArr2D(integ, i - 1, dd, val);
    }
    // int id0 = idx2d(integ, i, 0);
    // int id1 = idx2d(integ, i, 1);
    // integ.d[id0] = integ.d[idx2d(integ, i - 1, 0)] + distr.d[id0];
    // integ.d[id1] = integ.d[idx2d(integ, i - 1, 1)] + distr.d[id1];
  } //}}}

  // print distribution to output file //{{{
  // determine width of each column & collate data //{{{
  int columns = 5;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double **data = calloc(length, sizeof *data); // array for data
  for (int i = 0; i < length; i++) {
    data[i] = calloc(columns, sizeof *data[i]);
    count = -1;
    data[i][++count] = i + 1; // sequence length starts from 1
    data[i][++count] = GetArr2D(distr, i, 0);
    data[i][++count] = GetArr2D(distr, i, 1);
    data[i][++count] = GetArr2D(integ, i, 0);
    data[i][++count] = GetArr2D(integ, i, 1);
    // data[i][++count] = distr.d[idx2d(distr, i, 0)];
    // data[i][++count] = distr.d[idx2d(distr, i, 1)];
    // data[i][++count] = integ.d[idx2d(distr, i, 0)];
    // data[i][++count] = integ.d[idx2d(distr, i, 1)];
  }
  FillMaxDigits(columns, length, data, digits); //}}}
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  // print column headers
  fprintf(fw, "# Columns:");
  count = 0;
  fprintf(fw, " (%d) sequence length", ++count);
  fprintf(fw, ", (%d) bead0", ++count);
  fprintf(fw, ", (%d) bead1", ++count);
  fprintf(fw, ", (%d) bead0 (integrated)", ++count);
  fprintf(fw, ", (%d) bead1 (integrated)", ++count);
  putc('\n', fw);
  // print the data
  for (int i = 0; i < length; i++) {
    FprintfRow(fw, columns, data[i], digits);
    free(data[i]);
  }
  free(data);
  fclose(fw);
  //}}}

  // create random chain(s) from the distribution and print sequences //{{{
  char bname[2] = {'O', 'C'};
  int chains = 10;
  for (int xxx = 0; xxx < chains; xxx++) {
    // randomly pick starting type
    bool type = (bool)(pcg32Rand0Int(&rng, (uint64_t)2));
    int start = type;
    // numbers of beads of both types
    int beads[2] = {0, 0};
    // array of possible sequences
    int sequences[length];
    for (int i = 0; i < length; i++) {
      sequences[i] = -1;
    }
    // counter for number of sequences of each type
    int count_seq = 0;
    // total number of beads in the chain
    int total_beads = 0;
    while (true) {
      // random number between 0 and 1
      double percent = (double)(pcg32Rand0Int(&rng, (uint64_t)100)) / 100;
      // find sequence length it corresponds to
      for (int i = 0; i < length; i++) {
        // if (percent < integ.d[idx2d(integ, i, type)]) {
        if (percent < GetArr2D(integ, i, type)) {
          sequences[count_seq] = i + 1;
          beads[type] += i + 1;
          count_seq++;
          total_beads += i + 1;
          break;
        }
      }
      type = !type;
      if (beads[0] > num[0] ||
          beads[1] > num[1] ||
          total_beads > length) {
        for (int dd = 0; dd < 2; dd++) {
          beads[dd] = 0;
        }
        count_seq = 0;
        total_beads = 0;
        type = start;
      } else if (beads[0] == num[0] &&
                 beads[1] == num[1] &&
                 total_beads == length) {
        break;
      }
    }
    // pro-forma check //{{{
    if (beads[0] != num[0] || beads[1] != num[1] || total_beads != length) {
      snprintf(ERROR_MSG, LINE, "incorrectly generated chain; should be vs. is\n"
               "  bead0: %s%d%s vs %s%d%s\n"
               "  bead1: %s%d%s vs %s%d%s\n"
               "  total beads: %s%ld%s vs %s%d%s\n"
               "  (rng seed: %s%d%s)",
               ErrYellow(), num[0], ErrRed(), ErrYellow(), beads[0], ErrRed(),
               ErrYellow(), num[1], ErrRed(), ErrYellow(), beads[1], ErrRed(),
               ErrYellow(), length, ErrRed(), ErrYellow(), total_beads, ErrRed(),
               ErrYellow(), opt.seed, ErrRed());
      PrintError();
      exit(1);
    } //}}}

    // print sequence
    count = 0;
    type = !start;
    int count_beads = 0;
    vec3d coor = { .v = {0, 0, 0} };
    int sig = 1;
    char file[LINE];
    snprintf(file, LINE, "%s-%d.xyz", coor_file, xxx);
    fw = OpenFile(file, "w");
    fprintf(fw, "%d\n", (int)(length));
    fprintf(fw, "%d\n", (int)(length));
    int coor_length = 30;
    for (int i = 0; i < count_seq; i++) {
      for (int j = 0; j < sequences[count]; j++) {
        if ((count_beads % (coor_length + 2)) < (coor_length - 1)) {
          coor.v[0] += sig;
        } else if ((count_beads % (coor_length + 2)) < (coor_length + 1)) {
          coor.v[1]--;
        } else {
          coor.v[1]--;
          sig = -sig;
        }
        fprintf(fw, "%c %lf %lf %lf\n", bname[type], coor.v[0], coor.v[1], coor.v[2]);
        count_beads++;
      }
      type = !type;
      count++;
    }
    fclose(fw);
  }
  //}}}

  // free memory - to make valgrind happy //{{{
  free(chain);
  FreeArrND(distr);
  FreeArrND(integ);
  //}}}

  return 0;
}
