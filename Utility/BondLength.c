#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths in specified \
molecule type(s) for all bonds, printing results per bond type or adding \
all molecules' bonds as well (--all option). \
Note that input structure file with defined bonds must be used. \
The utility can also calculate distribution of \
distances between any two beads in those molecule types (-d option).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a distribution bin\n");
  fprintf(ptr, "<output>            output file with the distribution of "
          "bond lengths\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)");
  fprintf(ptr, "  --joined          specify that <input> contains joined "
          "coordinates\n");
  fprintf(ptr, "  --all             calculate distribution for each bond in "
          "the molecule type(s)\n");
  fprintf(ptr, "  -d <file> [ints]  write distribution of distances "
          "between specified bead pair(s) to <file> (if no [ints] are "
          "provided, the molecule's first and last beads are used)\n");
  fprintf(ptr, "  -w <float>        warn if the length exceeds <float> \n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join,         // --joined
       *mt,          // -m
       all;          // --all
  int n_list[100],   // -d (list of bead id pairs)
      n_number;      // -d (total number of beads in the pairs)
  double warn;       // -w option
  char d_file[LINE]; // -d (output file)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 5, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--joined", "--all", "-d", "-m", "-w");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output> - file name with bond length distribution
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  } //}}}
  opt->all = BoolOption(argc, argv, "--all");
  // '-d' option - specify bead ids to calculate distance between //{{{
  // FileIntegerOption(argc, argv, 0, 100, "-d",
  //                   opt->n_list, &opt->n_number, opt->d_file);
  FileNumbersOption(argc, argv, 0, 100, "-d", opt->n_list,
                    &opt->n_number, opt->d_file, 'i');
  // if '-d' is present without numbers, use first and last for each molecule
  int d_per_set = 2; // it's a bond, so there two beads in each
  if (opt->d_file[0] != '\0' && opt->n_number == 0) {
    opt->n_number = d_per_set;
    opt->n_list[0] = 1;
    opt->n_list[1] = HIGHNUM; // large number to specify last bead
  }
  int d_pair_n = opt->n_number / d_per_set;
  // Error: wrong number of integers //{{{
  if (opt->d_file[0] != '\0' && (opt->n_number % d_per_set) != 0) {
    err_msg("number of bead indexes must be even");
    PrintErrorOption("-d");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < opt->n_number; i += d_per_set) {
    if (opt->n_list[i] == opt->n_list[i+1] ||
        opt->n_list[i] == 0 || opt->n_list[i+1] == 0) {
      err_msg("each pair of bead ids must be non-zero and different");
      PrintErrorOption("-d");
      exit(1);
    }
  } //}}}
  //}}}
  // '-w' option - bond length warning //{{{
  if (!OneNumberOption(argc, argv, "-w", &opt->warn, 'd')) {
    opt->warn = HIGHNUM;
  } //}}}
  //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *box = System.Box.Length;

  // '-m <name(s)>' option
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  /*
   * number of bins: *10 because of the -d option; bondlength should be at most
   * half boxsize, but distance between any two beads in a molecule can be large
  */
  int bins = Max3(box[0], box[1], box[2]) / width * 10;

  // arrays for BeadType-BeadType bonds //{{{
  double ****bond_bt = calloc(Count->MoleculeType, sizeof *bond_bt),
         (***bond_bt_mma)[3] = calloc(Count->MoleculeType,
                                      sizeof (***bond_bt_mma)[3]);
  for (int i = 0; i < Count->MoleculeType; i++) {
    bond_bt[i] = calloc(Count->BeadType, sizeof *bond_bt);
    bond_bt_mma[i] = calloc(Count->BeadType, sizeof ***bond_bt_mma);
    for (int j = 0; j < Count->BeadType; j++) {
      bond_bt[i][j] = calloc(Count->BeadType, sizeof *bond_bt);
      bond_bt_mma[i][j] = calloc(Count->BeadType, sizeof ***bond_bt_mma);
      for (int k = 0; k < Count->BeadType; k++) {
        bond_bt[i][j][k] = calloc(bins, sizeof *bond_bt);
        // set high number as initial minium bond length
        bond_bt_mma[i][j][k][0] = HIGHNUM;
      }
    }
  } //}}}
  // arrays for all bonds in molecules //{{{
  double ***bond_all = NULL;
  double (**bond_all_mma)[3] = NULL;
  if (opt->all) {
    bond_all = calloc(Count->MoleculeType, sizeof *bond_all);
    bond_all_mma = calloc(Count->MoleculeType, sizeof (**bond_all_mma)[3]);
    for (int i = 0; i < Count->MoleculeType; i++) {
      bond_all[i] = calloc(System.MoleculeType[i].nBonds, sizeof *bond_all);
      bond_all_mma[i] = calloc(System.MoleculeType[i].nBonds,
                               sizeof **bond_all_mma);
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        bond_all[i][j] = calloc(bins, sizeof *bond_all);
        // set high number as initial minium bond length
        bond_all_mma[i][j][0] = HIGHNUM;
      }
    }
  } //}}}
  // extra arrays for -d option //{{{
  double ***bond_d = NULL;
  double (**bond_d_mma)[3] = NULL;
  if (opt->d_file[0] != '\0') {
    bond_d = calloc(Count->MoleculeType, sizeof *bond_d),
    bond_d_mma = calloc(Count->MoleculeType, sizeof (**bond_d_mma)[3]);
    for (int i = 0; i < Count->MoleculeType; i++) {
      bond_d[i] = calloc(d_pair_n, sizeof *bond_d);
      bond_d_mma[i] = calloc(d_pair_n, sizeof **bond_d_mma);
      for (int j = 0; j < d_pair_n; j++) {
        bond_d[i][j] = calloc(bins, sizeof *bond_d);
        // set high number as initial minimum
        bond_d_mma[i][j][0] = HIGHNUM;
      }
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
      WrapJoinCoordinates(&System, true, opt->join);
      // calculate bond lengths //{{{
      // go through all molecules
      // TODO: make into for (mtype); for (mtype.index)
      for (int i = 0; i < Count->Molecule; i++) {
        MOLECULE *mol_i = &System.Molecule[i];
        MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
        if (opt->mt[mol_i->Type]) { // use only specified molecule types
          for (int j = 0; j < mt_i->nBonds; j++) {
            // bead ids in the bond
            // TODO: use BondIndices
            int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                id2 = mol_i->Bead[mt_i->Bond[j][1]];
            BEAD *b_1 = &System.Bead[id1],
                 *b_2 = &System.Bead[id2];
            // bond length
            double bond[3];
            for (int dd = 0; dd < 3; dd++) {
              bond[dd] = b_1->Position[dd] - b_2->Position[dd];
            }
            bond[0] = VectLength(bond);
            // warn if bond is too long //{{{
            if (opt->warn != HIGHNUM && bond[0] > opt->warn) {
              snprintf(ERROR_MSG, LINE, "-w option; too long a bond between "
                       "beads %s%d%s and %s%d%s (%s%lf%s)",
                       ErrYellow(), id1, ErrCyan(), ErrYellow(), id2, ErrCyan(),
                       ErrYellow(), bond[0], ErrCyan());
              PrintWarning();
            } //}}}
            // btype1 must be lower than btype2
            int *id_lo, *id_hi;
            if (b_1->Type < b_2->Type) {
              id_lo = &b_1->Type;
              id_hi = &b_2->Type;
            } else {
              id_lo = &b_2->Type;
              id_hi = &b_1->Type;
            }

            // mins & maxes & averages //{{{
            if (bond[0] < bond_bt_mma[mol_i->Type][*id_lo][*id_hi][0]) {
              bond_bt_mma[mol_i->Type][*id_lo][*id_hi][0] = bond[0];
            } else if (bond[0] > bond_bt_mma[mol_i->Type][*id_lo][*id_hi][1]) {
              bond_bt_mma[mol_i->Type][*id_lo][*id_hi][1] = bond[0];
            }
            bond_bt_mma[mol_i->Type][*id_lo][*id_hi][2] += bond[0];
            if (opt->all) {
              if (bond[0] < bond_all_mma[mol_i->Type][j][0]) {
                bond_all_mma[mol_i->Type][j][0] = bond[0];
              } else if (bond[0] > bond_all_mma[mol_i->Type][j][1]) {
                bond_all_mma[mol_i->Type][j][1] = bond[0];
              }
              bond_all_mma[mol_i->Type][j][2] += bond[0];
            }
            //}}}

            int k = bond[0] / width;
            if (k < bins) {
              bond_bt[mol_i->Type][*id_lo][*id_hi][k]++;
              if (opt->all) {
                bond_all[mol_i->Type][j][k]++;
              }
            }
          }
        }
      } //}}}
      // calculate distance (-d option) //{{{
      if (opt->d_file[0] != '\0') {
        for (int i = 0; i < Count->Molecule; i++) {
          MOLECULE *mol_i = &System.Molecule[i];
          MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
          if (opt->mt[mol_i->Type]) { // use only specified molecule types
            for (int j = 0; j < opt->n_number; j += d_per_set) {
              // bead ids to use //{{{
              int id1, id2;
              // use first molecule bead if bead index too high or -1
              if (opt->n_list[j] >= mt_i->nBeads) {
                id1 = mol_i->Bead[mt_i->nBeads-1];
              } else { // use specified index otherwise
                id1 = mol_i->Bead[opt->n_list[j]-1];
              }
              // use last molecule bead if bead index too high or -1
              if (opt->n_list[j+1] >= mt_i->nBeads) {
                id2 = mol_i->Bead[mt_i->nBeads-1];
              } else { // use specified index otherwise
                id2 = mol_i->Bead[opt->n_list[j+1]-1];
              } //}}}
              BEAD *b_1 = &System.Bead[id1], *b_2 = &System.Bead[id2];
              // distance calculation
              double dist[3];
              for (int dd = 0; dd < 3; dd++) {
                dist[dd] = b_1->Position[dd] - b_2->Position[dd];
              }
              dist[0] = VectLength(dist);

              // distance mins & maxes & averages //{{{
              if (dist[0] < bond_d_mma[mol_i->Type][j/2][0]) { // minimum
                bond_d_mma[mol_i->Type][j/2][0] = dist[0];
              } else if (dist[0] > bond_d_mma[mol_i->Type][j/2][1]) { // maximum
                bond_d_mma[mol_i->Type][j/2][1] = dist[0];
              }
              bond_d_mma[mol_i->Type][j/2][2] += dist[0]; // average
              //}}}

              int k = dist[0] / width; // distribution 'bin'
              if (k < bins) {
                bond_d[mol_i->Type][j/2][k]++;
              }
            }
          }
        }
      } //}}}
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // sum up all bonds in molecules (normalization factor) //{{{
  // BeadType-BeadType bonds
  int ***bond_bt_norm = calloc(Count->MoleculeType, sizeof *bond_bt_norm),
      range[2] = {bins, 0}; // min/max distance
  for (int i = 0; i < Count->MoleculeType; i++) {
    bond_bt_norm[i] = calloc(Count->BeadType, sizeof *bond_bt_norm);
    for (int j = 0; j < Count->BeadType; j++) {
      bond_bt_norm[i][j] = calloc(Count->BeadType, sizeof *bond_bt_norm);
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = 0; l < bins; l++) {
          if (bond_bt[i][j][k][l] > 0) {
            bond_bt_norm[i][j][k] += bond_bt[i][j][k][l];
            if (l < range[0]) {
              range[0] = l;
            }
            if (l > range[1]) {
              range[1] = l;
            }
          }
        }
      }
    }
  }
  // all molecules' bonds
  int **bond_all_norm = NULL;
  if (opt->all) {
    bond_all_norm = calloc(Count->MoleculeType, sizeof *bond_all_norm);
    for (int i = 0; i < Count->MoleculeType; i++) {
      bond_all_norm[i] = calloc(System.MoleculeType[i].nBonds,
                                sizeof *bond_all_norm);
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        for (int k = 0; k < bins; k++) {
          bond_all_norm[i][j] += bond_all[i][j][k];
          if (bond_all[i][j][k] && k < range[0]) {
            range[0] = k;
          }
          if (bond_all[i][j][k] && k > range[1]) {
            range[1] = k;
          }
        }
      }
    }
  }
  // include nearest 0 values in the range of distances
  if (range[0] > 0) {
    range[0]--;
  }
  if (range[1] < (bins - 1)) {
    range[1] += 2; // +2 as for loop is range[0]...range[1]-1
  } //}}}

  // write distribution of bond lengths //{{{
  PrintByline(fout, argc, argv);
  // print first lines of output file - molecule names and beadtype pairs //{{{
  FILE *fw = OpenFile(fout, "a");
  fprintf(fw, "# (1) distance\n");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt->mt[i] && mt->nBonds > 0) {
      fprintf(fw, "# %s molecule:", mt->Name);
      for (int j = 0; j < mt->nBTypes; j++) {
        for (int k = j; k < mt->nBTypes; k++) {
          int btype1 = mt->BType[j],
              btype2 = mt->BType[k];
          if (bond_bt_norm[i][btype1][btype2] > 0) {
            count++;
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
          }
        }
      }
      if (opt->all) {
        count++;
        fprintf(fw, " (%d)-(%d) individual bonds", count,
                count + mt->nBonds - 1);
        count += mt->nBonds - 1;
      }
      putc('\n', fw);
    }
  } //}}}
  // write distribution to output file //{{{
  for (int i = range[0]; i < range[1]; i++) {
    fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
    for (int j = 0; j < Count->MoleculeType; j++) {
      MOLECULETYPE *mt = &System.MoleculeType[j];
      if (opt->mt[j] && mt->nBonds > 0) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < mt->nBTypes; k++) {
          for (int l = k; l < mt->nBTypes; l++) {
            int btype1 = mt->BType[k],
                btype2 = mt->BType[l];
            // btype1 must be lower than btype2
            if (btype1 > btype2) {
              SwapInt(&btype1, &btype2);
            }
            if (bond_bt_norm[j][btype1][btype2] > 0) {
              double value = bond_bt[j][btype1][btype2][i] / bond_bt_norm[j][btype1][btype2];
              fprintf(fw, "%10f", value);
            }
          }
        }
        if (opt->all) {
          for (int k = 0; k < mt->nBonds; k++) {
            if (bond_all_norm[j][k] > 0) {
              double value = bond_all[j][k][i] / bond_all_norm[j][k];
              fprintf(fw, "%10f", value);
            } else {
              fprintf(fw, "%10s", "?");
            }
          }
        }
      }
    }
    putc('\n', fw);
  } //}}}
  // write mins, maxes, and averages //{{{
  // legend line
  fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns)\n");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System.MoleculeType[i];
    if (opt->mt[i] && mt->nBonds > 0) {
      fprintf(fw, "# %s molecule:", mt->Name);
      for (int j = 0; j < mt->nBTypes; j++) {
        for (int k = j; k < mt->nBTypes; k++) {
          int btype1 = mt->BType[j],
              btype2 = mt->BType[k];
          if (bond_bt_norm[i][btype1][btype2] > 0) {
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
            count += 3;
          }
        }
      }
      if (opt->all) {
        fprintf(fw, " (%d)-(%d) individual bonds", count,
                count + 3 * mt->nBonds - 1);
        count += 3 * mt->nBonds;
      }
      putc('\n', fw);
    }
  }
  // data line
  putc('#', fw);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        // if this bin is filled, its max must be larger than 0
        if (bond_bt_mma[i][j][k][1] > 0) {
          fprintf(fw, " %lf", bond_bt_mma[i][j][k][0]);
          fprintf(fw, " %lf", bond_bt_mma[i][j][k][1]);
          fprintf(fw, " %lf", bond_bt_mma[i][j][k][2] / bond_bt_norm[i][j][k]);
        }
      }
    }
    if (opt->all) {
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        fprintf(fw, " %lf", bond_all_mma[i][j][0]);
        fprintf(fw, " %lf", bond_all_mma[i][j][1]);
        fprintf(fw, " %lf", bond_all_mma[i][j][2] / bond_all_norm[i][j]);
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // write distribution of distances from '-d' option //{{{
  if (opt->d_file[0] != '\0') {
    PrintByline(opt->d_file, argc, argv);
    // sum up all calculated distances (normalization factors) //{{{
    int d_norm[Count->MoleculeType][d_pair_n];
    range[0] = bins;
    range[1] = 0;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->mt[i]) {
        for (int j = 0; j < d_pair_n; j++) {
          d_norm[i][j] = 0;
          for (int k = 0; k < bins; k++) {
            if (bond_d[i][j][k] > 0) {
              d_norm[i][j] += bond_d[i][j][k];
              if (k < range[0]) {
                range[0] = k;
              }
              if (k > range[1]) {
                range[1] = k;
              }
            }
          }
        }
      }
    }
    // include nearest 0 values in the range of distances
    if (range[0] > 0) {
      range[0]--;
    }
    if (range[1] < (bins - 1)) {
      range[1] += 2; // +2 as for loop is range[0]...range[1]-1
    } //}}}
    fw = OpenFile(opt->d_file, "a");
    // print the first line - molecule names with bead order //{{{
    fprintf(fw, "# bead order in molecule(s) -");
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt->mt[i]) {
        fprintf(fw, " %s:", mt_i->Name);
        for (int j = 0; j < mt_i->nBeads; j++) {
          int btype = mt_i->Bead[j];
          fprintf(fw, " %s", System.BeadType[btype].Name);
        }
        putc(';', fw);
      }
    }
    putc('\n', fw); //}}}
    // print the second line - molecule names and ids with column numbers //{{{
    fprintf(fw, "# (1) distance\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt->mt[i]) {
        fprintf(fw, "# %s:", mt_i->Name);

        for (int j = 0; j < opt->n_number; j += d_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          if (opt->n_list[j] >= mt_i->nBeads &&
              opt->n_list[j+1] >= mt_i->nBeads) {
            continue;
          } //}}}
          // bead ids for the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt->n_list[j] >= mt_i->nBeads) {
            id1 = mt_i->nBeads;
          } else { // use specified index otherwise
            id1 = opt->n_list[j];
          }
          // use last molecule bead if bead index too high or -1
          if (opt->n_list[j+1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->n_list[j+1];
          }
          // write the numbers so the first is higher
          if (id1 > id2) {
            SwapInt(&id1, &id2);
          } //}}}
          fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
        }
      }
      putc('\n', fw);
    } //}}}
    // write the distribution to output file //{{{
    for (int i = range[0]; i < range[1]; i++) {
      fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
      for (int j = 0; j < Count->MoleculeType; j++) {
        if (opt->mt[j]) {
          for (int k = 0; k < opt->n_number; k += d_per_set) {
            // skip id pairs if both are too high for the molecule //{{{
            if (opt->n_list[k] >= System.MoleculeType[j].nBeads &&
                opt->n_list[k+1] >= System.MoleculeType[j].nBeads) {
              continue;
            } //}}}
            double value = (double)(bond_d[j][k/d_per_set][i]) /
                           d_norm[j][(int)(k/d_per_set)];
            fprintf(fw, " %10f", value);
          }
        }
      }
      putc('\n', fw);
    } //}}}
    // write mins and maxes
    // legend line //{{{
    fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns)\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt->mt[i]) {
        fprintf(fw, "# %s:", mt_i->Name);
        for (int j = 0; j < opt->n_number; j += d_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          if (opt->n_list[j] >= mt_i->nBeads &&
              opt->n_list[j+1] >= mt_i->nBeads) {
            continue;
          } //}}}
          // bead ids the distance //{{{
          int id1, id2;
          // use last molecule bead if bead index is too high
          if (opt->n_list[j] >= mt_i->nBeads) {
            id1 = mt_i->nBeads;
          } else { // use specified index otherwise
            id1 = opt->n_list[j];
          }
          // use last molecule bead if bead index is too high
          if (opt->n_list[j+1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->n_list[j+1];
          }
          // write the numbers so the first is lower
          if (id1 > id2) {
            SwapInt(&id1, &id2);
          } //}}}
          fprintf(fw, " (%d) %d-%d", count, id1, id2);
          count += 3;
        }
      }
      putc('\n', fw);
    } //}}}
    // data line //{{{
    putc('#', fw);
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->mt[i]) {
        for (int j = 0; j < opt->n_number; j += d_per_set) {
          // skip id pairs if both are too high for the molecule //{{{
          if (opt->n_list[j] >= System.MoleculeType[i].nBeads &&
              opt->n_list[j+1] >= System.MoleculeType[i].nBeads) {
            continue;
          } //}}}
          // if this bin is filled, its max must be larger than 0
          if (bond_d_mma[i][j/2][1] > 0) {
            fprintf(fw, " %lf", bond_d_mma[i][j/2][0]);
            fprintf(fw, " %lf", bond_d_mma[i][j/2][1]);
            fprintf(fw, " %lf", bond_d_mma[i][j/2][2] / d_norm[i][j]);
          }
        }
      }
    }
    putc('\n', fw); //}}}
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(opt->mt);
  // free arrays for BeadType-BeadType bonds
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        free(bond_bt[i][j][k]);
      }
      free(bond_bt[i][j]);
      free(bond_bt_mma[i][j]);
      free(bond_bt_norm[i][j]);
    }
    free(bond_bt[i]);
    free(bond_bt_mma[i]);
    free(bond_bt_norm[i]);
  }
  free(bond_bt);
  free(bond_bt_mma);
  free(bond_bt_norm);
  // free arrays for all bonds
  if (opt->all) {
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        free(bond_all[i][j]);
      }
      free(bond_all[i]);
      free(bond_all_mma[i]);
      free(bond_all_norm[i]);
    }
    free(bond_all);
    free(bond_all_mma);
    free(bond_all_norm);
  }
  // free arrays for -d option
  if (opt->d_file[0] != '\0') {
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < d_pair_n; j++) {
        free(bond_d[i][j]);
      }
      free(bond_d[i]);
      free(bond_d_mma[i]);
    }
    free(bond_d_mma);
    free(bond_d);
  }
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
