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

double Dot(double u[3], double v[3]) {
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}
void Vector(double coor1[3], double coor2[3], double v[3]) {
  for (int dd = 0; dd < 3; dd++) {
    v[dd] = coor1[dd] - coor2[dd];
  }
}
// angle between vectors i-j and k-l
static inline double CosAngle(double u[3], double v[3]) {
  return Dot(u, v) / (VectLength(u) * VectLength(v));
}

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
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // maximum number of bonds & beads //{{{
  int max_bonds = 0;
  if (opt->ne == HIGHNUM) { // no -ne option -> find longes molecule
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->mt[i] && System.MoleculeType[i].nBonds > max_bonds) {
        max_bonds = System.MoleculeType[i].nBonds;
      }
    }
  } else { // -ne option -> cannot be longer than the specified length
    max_bonds = opt->ne;
  } //}}}

  // arrays for sums //{{{
  double **lag_cos_phi = calloc(Count->MoleculeType, sizeof *lag_cos_phi),
         **lag_cos_phi2 = calloc(Count->MoleculeType, sizeof *lag_cos_phi2),
         (*avg_bond)[2] = calloc(Count->MoleculeType, sizeof avg_bond[2]);
  long int **count_stuff = calloc(Count->MoleculeType, sizeof *count_stuff);
  for (int i = 0; i < Count->MoleculeType; i++) {
    lag_cos_phi[i] = calloc(max_bonds, sizeof *lag_cos_phi[i]);
    lag_cos_phi2[i] = calloc(max_bonds, sizeof *lag_cos_phi2[i]);
    count_stuff[i] = calloc(max_bonds, sizeof *count_stuff[i]);
  }
  double *Re2 = calloc(Count->MoleculeType, sizeof *Re2);
  double *contour = calloc(Count->MoleculeType, sizeof *contour);
  double *ang_corr_avg = calloc(Count->MoleculeType, sizeof *ang_corr_avg);
  double **ang_corr_lag = calloc(Count->MoleculeType, sizeof *ang_corr_lag);
  double **ang_corr_lag2 = calloc(Count->MoleculeType, sizeof *ang_corr_lag2);
  double **cos_first_j = calloc(Count->MoleculeType, sizeof *cos_first_j);
  double **cos_last_j = calloc(Count->MoleculeType, sizeof *cos_last_j);
  long int *count_ang = calloc(Count->MoleculeType, sizeof *count_ang);
  long int **count_ang_lag = calloc(Count->MoleculeType, sizeof *count_ang_lag);
  for (int i = 0; i < Count->MoleculeType; i++) {
    // TODO: why +1 everywhere?
    ang_corr_lag[i] = calloc(max_bonds + 1, sizeof *ang_corr_lag[i]);
    ang_corr_lag2[i] = calloc(max_bonds + 1, sizeof *ang_corr_lag[i]);
    cos_first_j[i] = calloc(max_bonds + 1, sizeof *cos_first_j[i]);
    cos_last_j[i] = calloc(max_bonds + 1, sizeof *cos_last_j[i]);
    count_ang_lag[i] = calloc(max_bonds + 1, sizeof *count_ang_lag[i]);
  }
  int *count_mols = calloc(Count->MoleculeType, sizeof *count_mols);
  //}}}

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
        MOLECULETYPE *mt_i = &System.MoleculeType[i];
        // last bead id; corresponds to the highest bond with id last_b-1
        int last_b = mt_i->nBeads - 1;
        if (opt->ne != HIGHNUM) {
          last_b = opt->ne;
        }
        // use only specified molecule types that are long enough
        if (!opt->mt[i] &&
            mt_i->nBeads > (opt->ns + 1) &&
            mt_i->nBeads > last_b) {
          continue;
        }
        for (int mm = 0; mm < mt_i->Number; mm++) {
          int mol = mt_i->Index[mm];
          MOLECULE *mol_i = &System.Molecule[mol];
          if (!mol_i->InTimestep) {
            continue;
          }
          // vectors r_1,2 and r_1,N
          double first[3], last[3];
          Vector(System.Bead[mol_i->Bead[opt->ns]].Position,
                 System.Bead[mol_i->Bead[opt->ns+1]].Position, first);
          Vector(System.Bead[mol_i->Bead[opt->ns]].Position,
                 System.Bead[mol_i->Bead[last_b]].Position, last);
          // end-to-end distance squared
          Re2[mol_i->Type] += SqVectLength(last);
          // count molecules of every type
          count_mols[mol_i->Type]++;
          int bonds = last_b - opt->ns; // number of bonds to consider
          for (int lag = 0; lag < bonds; lag++) {
            // go from first bond (same id as first bead) to the highest
            // possible bond (its bond_id is last_bead-1)
            for (int j = opt->ns; j < (last_b - lag); j++) {
              int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                  id2 = mol_i->Bead[mt_i->Bond[j][1]],
                  id3 = mol_i->Bead[mt_i->Bond[j+lag][0]],
                  id4 = mol_i->Bead[mt_i->Bond[j+lag][1]];
              double u[3], v[3];
              Vector(System.Bead[id1].Position, System.Bead[id2].Position, u);
              Vector(System.Bead[id3].Position, System.Bead[id4].Position, v);
              // sum of cos(\phi)
              lag_cos_phi[mol_i->Type][lag] += CosAngle(u, v); // column 2
              // count values - easier than figuring out their number
              count_stuff[mol_i->Type][lag]++;
              if (lag == 0) { // for end-of-file average
                // contour length for end-to-end lp estimate
                contour[mol_i->Type] += VectLength(u);
                // bond lengths
                avg_bond[mol_i->Type][0] += VectLength(u);
                avg_bond[mol_i->Type][1]++;
              }
              if (j == opt->ns) { // for end-of-file average
                // angular correlation method - angle between r_1,2 & r_1,lag
                Vector(System.Bead[id1].Position, System.Bead[id4].Position, v);
                ang_corr_avg[mol_i->Type] += VectLength(v) * CosAngle(u, v);
                Vector(System.Bead[id3].Position, System.Bead[id4].Position, v);
                count_ang[mol_i->Type]++;
              }
              // bond correlation from the other end of the chain
              id1 = mol_i->Bead[mt_i->Bond[mt_i->nBonds-1-j][0]];
              id2 = mol_i->Bead[mt_i->Bond[mt_i->nBonds-1-j][1]];
              id3 = mol_i->Bead[mt_i->Bond[mt_i->nBonds-1-(j+lag)][0]];
              id4 = mol_i->Bead[mt_i->Bond[mt_i->nBonds-1-(j+lag)][1]];
              Vector(System.Bead[id1].Position, System.Bead[id2].Position, u);
              Vector(System.Bead[id3].Position, System.Bead[id4].Position, v);
              lag_cos_phi2[mol_i->Type][lag] += CosAngle(u, v); // column 3
            }
          }
          // the last allowed bead has id last_b, so +1 for number of beads
          for (int j = opt->ns; j < (last_b + 1 - 2); j++) {
            double u[3], v[3];
            int bin_id = j + 2 - opt->ns;
            // from the beginning of the chain
            int id11 = mol_i->Bead[opt->ns];
            int id12 = mol_i->Bead[opt->ns+1];
            int id21 = id11;
            int id22 = mol_i->Bead[j+2];
            Vector(System.Bead[id11].Position, System.Bead[id12].Position, u);
            Vector(System.Bead[id21].Position, System.Bead[id22].Position, v);
            // column 3:
            ang_corr_lag[mol_i->Type][bin_id] += VectLength(v) * CosAngle(u, v);
            // from the end of the chain
            id11 = mol_i->Bead[mt_i->nBeads-opt->ns-1];
            id12 = mol_i->Bead[mt_i->nBeads-opt->ns-2];
            id21 = id11;
            id22 = mol_i->Bead[mt_i->nBeads-j-3];
            Vector(System.Bead[id11].Position, System.Bead[id12].Position, u);
            Vector(System.Bead[id12].Position, System.Bead[id22].Position, v);
            // column 4:
            ang_corr_lag2[mol_i->Type][bin_id] += VectLength(v) * CosAngle(u, v);
            count_ang_lag[mol_i->Type][bin_id]++;
          }
          // go over all bonds; last_b is last bead's id and number of bonds
          // column 5 - from the beginning of the chain
          for (int j = opt->ns; j < last_b; j++) {
            double u[3], v[3];
            int bin_id = j - opt->ns;
            int id11 = mol_i->Bead[opt->ns],
                id12 = mol_i->Bead[opt->ns+1],
                id21 = mol_i->Bead[mt_i->Bond[j][0]],
                id22 = mol_i->Bead[mt_i->Bond[j][1]];
            Vector(System.Bead[id11].Position, System.Bead[id12].Position, u);
            Vector(System.Bead[id21].Position, System.Bead[id22].Position, v);
            // column 5:
            cos_first_j[mol_i->Type][bin_id] += CosAngle(u, v);
          }
          // column 5 - from the end of the chain
          int last_bond_id = mt_i->nBonds - 1 - opt->ns;
          for (int j = last_bond_id; j >= (mt_i->nBonds - last_b); j--) {
            double u[3], v[3];
            int bin_id = last_bond_id - j;
            int id11 = mol_i->Bead[mt_i->Bond[j][0]],
                id12 = mol_i->Bead[mt_i->Bond[j][1]],
                id21 = mol_i->Bead[mt_i->Bond[last_bond_id][0]],
                id22 = mol_i->Bead[mt_i->Bond[last_bond_id][1]];
            Vector(System.Bead[id11].Position, System.Bead[id12].Position, u);
            Vector(System.Bead[id21].Position, System.Bead[id22].Position, v);
            // column 6:
            cos_last_j[mol_i->Type][bin_id] += CosAngle(u, v);
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

  // average all arrays //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      for (int lag = 0; lag < (max_bonds - 1); lag++) {
        lag_cos_phi[i][lag] /= count_stuff[i][lag];
        lag_cos_phi2[i][lag] /= count_stuff[i][lag];
      }
    }
  } //}}}

  // write data to ouptut file //{{{
  // print first lines of output file //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# for each molecule type: bond correlation; ");
  fprintf(fw, "angle r_1,2 and r_1,N multiplied by r_1,N distance; ");
  fprintf(fw, "angle r_1,2 and r_N-1,N; ");
  fprintf(fw, "integrated angle r_1,2 and r_N-1,N\n");
  fprintf(fw, "# (1) distance between bonds/beads;");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      count++;
      fprintf(fw, " (%d-%d) %s", count, count + 9, mt_i->Name);
    }
  }
  putc('\n', fw); //}}}
  // determine width of each column & collate data //{{{
  int columns = Count->MoleculeType * 10 + 1;
  int digits[columns][2];
  InitInt2DArray((int *)digits, columns, 2, 0);
  double *data[max_bonds+1];
  for (int lag = 0; lag < (max_bonds + 1); lag++) {
    data[lag] = calloc(columns, sizeof data[lag]);
    count = -1;
    data[lag][++count] = lag;
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (!opt->mt[j]) {
        continue;
      }
      MOLECULETYPE *mt_j = &System.MoleculeType[j];
      int last_b = mt_j->nBeads - 1;
      if (opt->ne != HIGHNUM) {
        last_b = opt->ne;
      }
      if (lag < mt_j->nBonds && lag < last_b) { // columns 2 & 3
        data[lag][++count] = lag_cos_phi[j][lag];
        data[lag][++count] = lag_cos_phi2[j][lag];
      } else {
        count += 2;
      }
      if (lag < 2) { // columns 4 & 5
        count += 2;
      } else if (lag < mt_j->nBeads && lag < (last_b + 1)) {
        data[lag][++count] = ang_corr_lag[j][lag] / count_ang_lag[j][lag];
        data[lag][++count] = ang_corr_lag2[j][lag] / count_ang_lag[j][lag];
      }
      if (lag < mt_j->nBonds && lag < last_b) { // columns 6 & 7
        data[lag][++count] = cos_first_j[j][lag] / count_mols[j];
        data[lag][++count] = cos_last_j[j][lag] / count_mols[j];
      } else {
        count += 2;
      }
      if (lag == 0) { // columns 8 & 9
        data[lag][++count] = 1;
        data[lag][++count] = 1;
      } else if (lag < mt_j->nBonds && lag < last_b) {
        count++;
        data[lag][count] = data[lag-1][count] + cos_first_j[j][lag] / count_mols[j];
        count++;
        data[lag][count] = data[lag-1][count] + cos_last_j[j][lag] / count_mols[j];
      } else {
        count += 2;
      }
      if (lag == 0) { // columns 10 & 11
        data[lag][++count] = 1;
        data[lag][++count] = 1;
      } else if (lag < mt_j->nBonds && lag < last_b) {
        count++;
        data[lag][count] = data[lag-1][count] + lag_cos_phi[j][lag];
        count++;
        data[lag][count] = data[lag-1][count] + lag_cos_phi2[j][lag];
      } else {
        count += 2;
      }
    }
  }
  FillMaxDigits(columns, max_bonds + 1, data, digits); //}}}
  // <l_p>\approx\sum_{i=1}^N<b1.bi> ... N is number of bonds
  double lp_kp[Count->MoleculeType];
  double lp_kp2[Count->MoleculeType];
  InitDoubleArray(lp_kp, Count->MoleculeType, 0);
  InitDoubleArray(lp_kp2, Count->MoleculeType, 0);
  for (int lag = opt->ns; lag < (max_bonds + 1); lag++) {
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt->mt[j]) {
        lp_kp[j] += cos_first_j[j][lag];
        lp_kp2[j] += cos_last_j[j][lag];
      }
    }
  }
  for (int lag = 0; lag < (max_bonds + 1 - opt->ns); lag++) {
    count = 0;
    Fprintf1(fw, data[lag][count], digits[count]);
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (!opt->mt[j]) {
        continue;
      }
      MOLECULETYPE *mt_j = &System.MoleculeType[j];
      int last_b = mt_j->nBeads - 1;
      if (opt->ne != HIGHNUM) {
        last_b = opt->ne;
      }
      count += 2; // columns 2 & 3
      if (lag < (last_b - opt->ns)) {
        Fprintf1(fw, data[lag][count-1], digits[count-1]);
        Fprintf1(fw, data[lag][count], digits[count]);
      } else {
        fprintf(fw, " %*s", digits[count-1][0] + digits[count-1][1], "?");
        fprintf(fw, " %*s", digits[count][0] + digits[count][1], "?");
      }
      count += 2; // columns 4 & 5
      if (lag < 2 || lag > (last_b - opt->ns)) {
        fprintf(fw, " %*s", digits[count-1][0] + digits[count-1][1], "?");
        fprintf(fw, " %*s", digits[count][0] + digits[count][1], "?");
      } else {
        Fprintf1(fw, data[lag][count-1], digits[count-1]);
        Fprintf1(fw, data[lag][count], digits[count]);
      }
      count += 2; // columns 6 & 7
      if (lag < (last_b - opt->ns)) {
        Fprintf1(fw, data[lag][count-1], digits[count-1]);
        Fprintf1(fw, data[lag][count], digits[count]);
      } else {
        fprintf(fw, " %*s", digits[count-1][0] + digits[count-1][1], "?");
        fprintf(fw, " %*s", digits[count][0] + digits[count][1], "?");
      }
      count += 2; // columns 8 & 9
      if (lag < (last_b - opt->ns)) {
        Fprintf1(fw, data[lag][count-1], digits[count-1]);
        Fprintf1(fw, data[lag][count], digits[count]);
      } else {
        fprintf(fw, " %*s", digits[count-1][0] + digits[count-1][1], "?");
        fprintf(fw, " %*s", digits[count][0] + digits[count][1], "?");
      }
      count += 2; // columns 10 & 11
      if (lag < (last_b - opt->ns)) {
        Fprintf1(fw, data[lag][count-1], digits[count-1]);
        Fprintf1(fw, data[lag][count], digits[count]);
      } else {
        fprintf(fw, " %*s", digits[count-1][0] + digits[count-1][1], "?");
        fprintf(fw, " %*s", digits[count][0] + digits[count][1], "?");
      }
    }
    putc('\n', fw);
    free(data[lag]);
  }
  for (int lag = (max_bonds + 1 - opt->ns); lag < (max_bonds + 1); lag++) {
    free(data[lag]);
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      fprintf(fw, "# l_p from: (1) R_e,");
      fprintf(fw, " (2) angle correlation-all angles; ");
      fprintf(fw, " (3) angle correlation-first to i-th angle");
      double lp = Re2[i] / (2 * contour[i]);
      fprintf(fw, "\n# %lf", lp);
      int beads = System.MoleculeType[i].nBeads - 1 - opt->ns;
      lp = ang_corr_avg[i] / (count_mols[i] * beads);
      fprintf(fw, " %lf", lp);
      fprintf(fw, " %lf", lp_kp[i] / count_mols[i]);
      fprintf(fw, " %lf", lp_kp2[i] / count_mols[i]);
      fprintf(fw, "\n# average bond length: %lf\n",
              avg_bond[i][0] / avg_bond[i][1]);
    }
  }
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  free(opt->mt);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(lag_cos_phi[i]);
    free(lag_cos_phi2[i]);
    free(count_stuff[i]);
  }
  free(lag_cos_phi);
  free(lag_cos_phi2);
  free(count_stuff);
  free(avg_bond);
  free(Re2);
  free(ang_corr_avg);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(ang_corr_lag[i]);
    free(ang_corr_lag2[i]);
    free(cos_first_j[i]);
    free(cos_last_j[i]);
    free(count_ang_lag[i]);
  }
  free(ang_corr_lag);
  free(ang_corr_lag2);
  free(cos_first_j);
  free(cos_last_j);
  free(count_ang_lag);
  free(contour);
  free(count_ang);
  free(count_mols);
  FreeSystem(&System);
  free(opt);
  //}}}

  return 0;
}
