#include "../AnalysisTools.h"

// TODO: very messy!!!

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
PersistenceLength calculates correlation of bond vectors and its standard \
deviation. To get the peristence length, the data must be fitted via, \
typically, an exponential function. Note the utility expects a linear chain \
with ordered beads ids (e.g., for a 4-bead chain, the order must be 1-2-3-4, \
leading to ordered bonds 1-2, 2-3, 3-4; connectivity like 1-4-2-3 with \
bonds 1-4, 2-4, 2-3 could give unexpected results).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file with the persistence length\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)\n");
  fprintf(ptr, "  --joined          specify that <input> contains joined "
          "coordinates\n");
  fprintf(ptr, "  -ns <int>         start with <int>-th bead in a molecule\n");
  fprintf(ptr, "  -ne <int>         end with <int>-th bead in a molecule\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join,  // --joined
       *mt;   // -m
  int ns, ne; // -ns/-ne; first bead and bond and last bead and bond
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  int common = 8, all = common + 4, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--joined", "-m", "-ns", "-ne");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output> - file name with persistence lengths
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // --joined option
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  }
  if (!OneNumberOption(argc, argv, "-ns", &opt->ns, 'i')) {
    opt->ns = 1;
  }
  opt->ns--; // indexing starts from 0
  if (!OneNumberOption(argc, argv, "-ne", &opt->ne, 'i')) {
    opt->ne = HIGHNUM;
  } else {
    opt->ne--; // indexing starts from 0
  }
  if (opt->ns != HIGHNUM && opt->ne != HIGHNUM &&
      (opt->ne - opt->ns) < 2) {
    err_msg("at least three beads are necessary, i.e., <-ns> - <-ne> > 1; "
            "(note that the calculation is meaningful only for longer chains)");
    PrintErrorOption("-ns/-ne");
    exit(1);
  }
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-m <name(s)>' option
  opt->mt = calloc(Count->MoleculeType, sizeof *opt->mt);
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // maximum number of bonds & beads //{{{
  int max_bonds = 0;
  if (opt->ne == HIGHNUM) { // no -ne option -> find longest molecule
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->mt[i] && System.MoleculeType[i].nBonds > max_bonds) {
        max_bonds = System.MoleculeType[i].nBonds;
      }
    }
  } else { // -ne option -> cannot be longer than the specified length
    max_bonds = opt->ne;
  }
  if (max_bonds < opt->ns) {
    err_msg("starting bead is larger than the length of any molecule");
    ErrorOption("-ns");
    exit(1);
  } //}}}

  // TODO: define S1 throuh S3 via kp's text //{{{
  double (**S1)[2] = calloc(Count->MoleculeType, sizeof *S1);
  double **S2 = calloc(Count->MoleculeType, sizeof *S2);
  int **count_S2 = calloc(Count->MoleculeType, sizeof *S2);
  double (**S3)[2] = calloc(Count->MoleculeType, sizeof *S3);
  int **count_S3 = calloc(Count->MoleculeType, sizeof *S3);
  double *bondlength = calloc(Count->MoleculeType, sizeof *bondlength);
  int *count_bonds = calloc(Count->MoleculeType, sizeof *count_bonds);
  for (int i = 0; i < Count->MoleculeType; i++) {
    S1[i] = calloc(max_bonds, sizeof *S1[i]);
    S2[i] = calloc(max_bonds, sizeof *S2[i]);
    count_S2[i] = calloc(max_bonds, sizeof *count_S2[i]);
    S3[i] = calloc(max_bonds, sizeof *S3[i]);
    count_S3[i] = calloc(max_bonds, sizeof *count_S3[i]);
    for (int j = 0; j < Count->MoleculeType; j++) {
      S1[i][j][0] = 0;
      S1[i][j][1] = 0;
      S3[i][j][0] = 0;
      S3[i][j][1] = 0;
    }
  } //}}}

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
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, false, opt->join);
      // go through all molecules //{{{
      for (int i = 0; i < Count->MoleculeType; i++) {
        MOLECULETYPE *mt = &System.MoleculeType[i];
        // last bond id
        int last_bond = mt->nBonds;
        if (opt->ne != HIGHNUM) {
          last_bond = opt->ne;
        }
        int first_bond = 0;
        if (opt->ns > 0) {
          first_bond = opt->ns;
        }
        // use only specified molecule types that are long enough
        if (!opt->mt[i] || mt->nBonds < opt->ns) {
          continue;
        }

        for (int j = 0; j < mt->Number; j++) {
          MOLECULE *mol = &System.Molecule[mt->Index[j]];
          // calculate average bond length in the chain (for S3)
          double l2_sum = 0;
          int nb_sum = 0;
          for (int b = 0; b < mt->nBonds; b++) {
            int a = mol->Bead[mt->Bond[b][0]];
            int c = mol->Bead[mt->Bond[b][1]];
            // double v[3];
            // Vector(System.Bead[a].Position, System.Bead[c].Position, v);
            vec3 v = Vector(System.Bead[a].Position, System.Bead[c].Position);
            l2_sum += Dot(v, v);
            nb_sum += 1;
          }
          l2_sum /= nb_sum;
          // S1 function
          // first bond vector (for S1)
          int b1 = mol->Bead[mt->Bond[first_bond][0]],
              b2 = mol->Bead[mt->Bond[first_bond][1]];
          vec3 bond1 = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
          // last bond vector (for reversed S1)
          b1 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][0]],
          b2 = mol->Bead[mt->Bond[mt->nBonds-first_bond-1][1]];
          vec3 bondN = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
          for (int k = first_bond; k < last_bond; k++) {

            // S1 function & bondlengths //{{{
            // TODO: bondj -> bondk
            // 1->N S1
            b1 = mol->Bead[mt->Bond[k][0]];
            b2 = mol->Bead[mt->Bond[k][1]];
            vec3 bondj = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
            S1[mol->Type][k-first_bond][0] += CosAngle(bondj, bond1);
            // bondlength & count bonds
            bondlength[mol->Type] += VectLength(bondj);
            count_bonds[mol->Type]++;
            // reverse S1
            int bond_id = mt->nBonds - k - 1;
            int bin_id = k - first_bond;
            b1 = mol->Bead[mt->Bond[bond_id][0]];
            b2 = mol->Bead[mt->Bond[bond_id][1]];
            bondj = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
            S1[mol->Type][bin_id][1] += CosAngle(bondN, bondj); //}}}
            for (int l = k; l < last_bond; l++) {
              int lag = l - k;
              // S2 function (classic bond correlation) //{{{
              // first bond vector
              int b1 = mol->Bead[mt->Bond[k][0]],
                  b2 = mol->Bead[mt->Bond[k][1]];
              vec3 bondk = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
              // second bond vector
              b1 = mol->Bead[mt->Bond[l][0]];
              b2 = mol->Bead[mt->Bond[l][1]];
              vec3 bondl = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
              // autocorrelation
              S2[mol->Type][lag] += CosAngle(bondk, bondl);
              count_S2[mol->Type][lag]++; //}}}
              // S3 function (end-to-end distances) //{{{
              b1 = mol->Bead[mt->Bond[k][0]];
              b2 = mol->Bead[mt->Bond[l][1]];
              vec3 Re = Vector(System.Bead[b1].Position, System.Bead[b2].Position);
              double denom = Square(lag + 1) * l2_sum;
              S3[mol->Type][lag][0] += Dot(Re, Re) / denom;
              S3[mol->Type][lag][1] += Dot(Re, Re);
              count_S3[mol->Type][lag]++; //}}}
            }
          }
        }
      } //}}}
    //}}}
    } else { //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // write data //{{{
  // determine width of each column & collate data //{{{
  int datalines = max_bonds - opt->ns;
  // count used molecule types
  count = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      count++;
    }
  }
  int data_per_mtype = 9;
  int columns = count * data_per_mtype + 1;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double *data[datalines];
  // arrays for integrated functions
  double sum_S1[Count->MoleculeType][2]; // [0] ... 1->N; [1] ... reverse
  double sum_S2[Count->MoleculeType];
  double sum_S3[Count->MoleculeType][2]; // [0] ... normalized; [1] .. raw
  // average bond length
  for (int i = 0; i < Count->MoleculeType; i++) {
    sum_S1[i][0] = 0;
    sum_S1[i][1] = 0;
    sum_S2[i] = 0;
    sum_S3[i][0] = 0;
    sum_S3[i][1] = 0;
  }
  for (int lag = 0; lag < datalines; lag++) {
    data[lag] = calloc(columns, sizeof *data[lag]);
    count = -1;
    data[lag][++count] = lag;
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt->mt[j]) {
        // S1 function (from either end)
        double avg[2] = {S1[j][lag][0] / count_used,
                         S1[j][lag][1] / count_used};
        for (int dd = 0; dd < 2; dd++) {
          sum_S1[j][dd] += avg[dd];
          data[lag][++count] = avg[dd];
          data[lag][++count] = sum_S1[j][dd];
        }
        // S2 (autocorrelation function)
        avg[0] = S2[j][lag] / count_S2[j][lag];
        sum_S2[j] += avg[0];
        data[lag][++count] = avg[0];
        data[lag][++count] = sum_S2[j];
        // S3 (end-to-end distances)
        for (int dd = 0; dd < 2; dd++) {
          avg[dd] = S3[j][lag][dd] / count_S3[j][lag];
          sum_S3[j][dd] += avg[dd];
          data[lag][++count] = avg[dd];
          data[lag][++count] = sum_S3[j][dd];
        }
      }
    }
  }
  FillMaxDigits(columns, datalines, data, digits); //}}}
  // print the data //{{{
  // headers
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# for each molecule type: ");
  count = 1;
  fprintf(fw, "(%d) S1, ", count++);
  fprintf(fw, "(%d) int S1, ", count++);
  fprintf(fw, "(%d) S1 (rev), ", count++);
  fprintf(fw, "(%d) int S1 (rev), ", count++);
  fprintf(fw, "(%d) S2, ", count++);
  fprintf(fw, "(%d) int S2", count++);
  fprintf(fw, "(%d) S3 (normalised), ", count++);
  fprintf(fw, "(%d) int S3 (normalised)", count++);
  fprintf(fw, "(%d) S3 (raw)", count++);
  putc('\n', fw);
  fprintf(fw, "# ");
  count = 1;
  fprintf(fw, "(%d) lag; ", count++);
  fprintf(fw, "molecule types: ");
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt->mt[i]) {
      fprintf(fw, "(%d)-(%d) %s", count, count+data_per_mtype-1, mt->Name);
      count += data_per_mtype;
      if (i != (Count->MoleculeType - 1)) {
        fprintf(fw, ", ");
      }
    }
  }
  putc('\n', fw);
  // datalines
  for (int lag = 0; lag < datalines; lag++) {
    FprintfRow(fw, columns, data[lag], digits);
    free(data[lag]);
  } //}}}
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  free(opt->mt);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(S1[i]);
    free(S2[i]);
    free(count_S2[i]);
    free(S3[i]);
    free(count_S3[i]);
  }
  free(S1);
  free(S2);
  free(count_S2);
  free(S3);
  free(count_S3);
  free(bondlength);
  free(count_bonds);
  FreeSystem(&System);
  free(opt);
  //}}}

  return 0;
}
