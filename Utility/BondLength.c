#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths in specified \
molecule type(s) for all bonds. \
Note that input structure file with defined bonds must be used. \
The utility can also calculate distribution of \
distances between any two beads in those molecule types (-d option).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a single distribution bin\n");
  fprintf(ptr, "<output>            output file with the distribution of "
          "bond lengths\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)");
  fprintf(ptr, "  --joined          specify that <input> contains joined "
          "coordinates\n");
  fprintf(ptr, "  -d <file> [ints]  write distribution of distances "
          "between specified bead pair(s) to <file> (if no [ints] are "
          "provided, the molecule's first and last beads are used)\n");
  fprintf(ptr, "  -w <float>        warn if the length exceeds <float> \n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join,         // --joined
       *mt;          // -m
  int d_pair[100],   // -d (list of bead id pairs)
      d_number;      // -d (total number of beads in the pairs)
  double warn;       // -w option
  char d_file[LINE]; // -d (output file)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// calculate index in an 1D array that simulates a 2D one /{{{
int id2D(int i1, int i2, int size) {
  return (i2 * size + i1);
} //}}}
// calculate index in an 1D array that simulates a 3D one /{{{
int id3D(int i1, int i2, int i3, int size[2]) {
  return (i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
}
// return three indices (assumes 3D array size 3*size[0]*size[1])
int * bins_id3D(int i1, int i2, int size[2]) {
  static int bin[3];
  for (int dd = 0; dd < 3; dd++) {
    bin[dd] = id3D(i1, i2, dd, size);
  }
  return bin;
} //}}}
// calculate index in an 1D array that simulates a 4D one /{{{
int id4D(int i1, int i2, int i3, int i4, int size[3]) {
  return (i4 * size[0] * size[1] * size[2] +
          i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
}
// return three indices (assumes 4D array size 3*size[0]*size[1]*size[2])
int * bins_id4D(int i1, int i2, int i3, int size[3]) {
  static int bin[3];
  for (int dd = 0; dd < 3; dd++) {
    bin[dd] = id4D(i1, i2, i3, dd, size);
  }
  return bin;
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 8, all = common + 4, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "--joined");
  strcpy(option[count++], "-d");
  strcpy(option[count++], "-m");
  strcpy(option[count++], "-w");
  OptionCheck(argc, argv, count, req_arg, common, all, option, true); //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  snprintf(in.coor.name, LINE, "%s", argv[++count]);
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
  snprintf(fout, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, LINE, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  } //}}}
  // '-d' option - specify bead ids to calculate distance between //{{{
  FileIntegerOption(argc, argv, 0, 100, "-d",
                    opt->d_pair, &opt->d_number, opt->d_file);
  // if '-d' is present without numbers, use first and last for each molecule
  int d_per_set = 2; // it's a bond, so there two beads in each
  if (opt->d_file[0] != '\0' && opt->d_number == 0) {
    opt->d_number = d_per_set;
    opt->d_pair[0] = 1;
    opt->d_pair[1] = HIGHNUM; // large number to specify last bead
  }
  int d_pair_n = opt->d_number / d_per_set;
  // Error: wrong number of integers //{{{
  if (opt->d_file[0] != '\0' && (opt->d_number % d_per_set) != 0) {
    strcpy(ERROR_MSG, "number of bead indexes must be even");
    PrintErrorOption("-d");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < opt->d_number; i += d_per_set) {
    if (opt->d_pair[i] == opt->d_pair[i+1] ||
        opt->d_pair[i] == 0 || opt->d_pair[i+1] == 0) {
      strcpy(ERROR_MSG, "each pair of bead ids must be non-zero and different");
      PrintErrorOption("-d");
      exit(1);
    }
  }
  //}}}
  for (int i = 0; i < opt->d_number; i++) {
    if (opt->d_pair[i] != HIGHNUM) {
      opt->d_pair[i]--; // ids should start with zero (or is -1 if none specified)
    }
  } //}}}
  // '-w' option - bond length warning //{{{
  if (!DoubleOption1(argc, argv, "-w", &opt->warn)) {
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
  if (!MoleculeTypeOption(argc, argv, "-m", true, opt->mt, System)) {
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
  int arr_bt[3];
  arr_bt[0] = Count->MoleculeType;
  arr_bt[1] = Count->BeadType;
  arr_bt[2] = Count->BeadType;
  int arr_bt_size = arr_bt[0] * arr_bt[1] * arr_bt[2];
  double *length_bt = calloc(arr_bt_size * bins, sizeof *length_bt),
         *length_bt_mma = calloc(arr_bt_size * 3, sizeof *length_bt_mma);
  for (int i = 0; i < arr_bt[0]; i++) {
    for (int j = 0; j < arr_bt[1]; j++) {
      for (int k = 0; k < arr_bt[2]; k++) {
        int id = id4D(i, j, k, 0, arr_bt);
        length_bt_mma[id] = 10 * Max3(box[0], box[1], box[2]);
      }
    }
  } //}}}
  // arrays for all bonds in molecules //{{{
  int max_bonds = 0; // find MoleculeType with the most bonds
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (System.MoleculeType[i].nBonds > max_bonds) {
      max_bonds = System.MoleculeType[i].nBonds;
    }
  }
  int arr_bonds[2];
  arr_bonds[0] = Count->MoleculeType;
  arr_bonds[1] = max_bonds;
  int arr_bonds_size = arr_bonds[0] * arr_bonds[1];
  double *length_bonds = calloc(arr_bonds_size * bins, sizeof *length_bt),
         *length_bonds_mma = calloc(arr_bonds_size * 3, sizeof *length_bt_mma);
  for (int i = 0; i < arr_bonds[0]; i++) {
    for (int j = 0; j < arr_bonds[1]; j++) {
      int id = id3D(i, j, 0, arr_bonds);
      length_bonds_mma[id] = 10 * Max3(box[0], box[1], box[2]);
    }
  } //}}}
  // extra arrays for -d option //{{{
  int arr_d[2];
  arr_d[0] = Count->MoleculeType;
  arr_d[1] = d_pair_n;
  int arr_d_size = arr_d[0] * arr_d[1];
  double *d_dist = calloc(arr_d_size * bins, sizeof *d_dist),
         *d_mma = calloc(arr_d_size * 3, sizeof *d_mma);
  for (int i = 0; i < arr_d[0]; i++) {
    for (int j = 0; j < arr_d[1]; j++) {
      int id = id3D(i, j, 0, arr_d);
      d_mma[id] = 10 * Max3(box[0], box[1], box[2]);
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
    if (count_coor >= opt->c.start &&
        (count_coor <= opt->c.end || opt->c.end == -1) &&
        ((count_coor - opt->c.start) % opt->c.skip) == 0) {
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
            bond[0] = VECTORLENGTH(bond);
            // warn if bond is too long //{{{
            if (opt->warn != HIGHNUM && bond[0] > opt->warn) {
              snprintf(ERROR_MSG, LINE, "-w option; too long a bond between "
                       "beads %s%d%s and %s%d%s (%s%lf%s)",
                       ErrYellow(), id1, ErrCyan(), ErrYellow(), id2, ErrCyan(),
                       ErrYellow(), bond[0], ErrCyan());
              PrintWarning();
            } //}}}
            // btype1 must be lower then btype2
            int *id_lo, *id_hi;
            if (b_1->Type < b_2->Type) {
              id_lo = &b_1->Type;
              id_hi = &b_2->Type;
            } else {
              id_lo = &b_2->Type;
              id_hi = &b_1->Type;
            }

            // mins & maxes & averages //{{{
            int *bin_id = bins_id4D(mol_i->Type, *id_lo, *id_hi, arr_bt);
            if (bond[0] < length_bt_mma[bin_id[0]]) {
              length_bt_mma[bin_id[0]] = bond[0];
            } else if (bond[0] > length_bt_mma[bin_id[1]]) {
              length_bt_mma[bin_id[1]] = bond[0];
            }
            length_bt_mma[bin_id[2]] += bond[0];
            //}}}

            int k = bond[0] / width;
            if (k < bins) {
              length_bt[id4D(mol_i->Type, *id_lo, *id_hi, k, arr_bt)]++;
              length_bonds[id3D(mol_i->Type, j, k, arr_bonds)]++;
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
            for (int j = 0; j < opt->d_number; j += d_per_set) {
              // bead ids to use //{{{
              int id1, id2;
              // use first molecule bead if bead index too high or -1
              if (opt->d_pair[j] == -1 || opt->d_pair[j] == HIGHNUM) {
                id1 = mol_i->Bead[0];
              } else { // use specified index otherwise
                id1 = mol_i->Bead[opt->d_pair[j]];
              }
              // use last molecule bead if bead index too high or -1
              if (opt->d_pair[j+1] == -1 || opt->d_pair[j+1] == HIGHNUM) {
                id2 = mol_i->Bead[mt_i->nBeads - 1];
              } else { // use specified index otherwise
                id2 = mol_i->Bead[opt->d_pair[j+1]];
              } //}}}
              BEAD *b_1 = &System.Bead[id1], *b_2 = &System.Bead[id2];
              // distance calculation
              double dist[3];
              for (int dd = 0; dd < 3; dd++) {
                dist[dd] = b_1->Position[dd] - b_2->Position[dd];
              }
              dist[0] = VECTORLENGTH(dist);

              // distance mins & maxes & averages //{{{
              int *bin_id = bins_id3D(mol_i->Type, j / 2, arr_d);
              if (dist[0] < d_mma[bin_id[0]]) { // minimum
                d_mma[bin_id[0]] = dist[0];
              } else if (dist[0] > d_mma[bin_id[1]]) { // maximum
                d_mma[bin_id[1]] = dist[0];
              }
              d_mma[bin_id[2]] += dist[0]; // average
              //}}}

              int k = dist[0] / width; // distribution 'bin'
              if (k < bins) {
                d_dist[id3D(mol_i->Type, j / 2, k, arr_d)]++;
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
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // sum up all bonds in molecules (normalization factor) //{{{
    // BeadType-BeadType bonds
  int *bt_norm = calloc(arr_bt_size, sizeof *bt_norm);
  for (int i = 0; i < arr_bt[0]; i++) {
    for (int j = 0; j < arr_bt[1]; j++) {
      for (int k = j; k < arr_bt[2]; k++) {
        for (int l = 0; l < bins; l++) {
          int bin = id4D(i, j, k, l, arr_bt);
          if (length_bt[bin] > 0) {
            int id = id3D(i, j, k, arr_bt);
            bt_norm[id] += length_bt[bin];
          }
        }
      }
    }
  }
  // all molecules' bonds
  int *bonds_norm = calloc(arr_bonds_size, sizeof *bonds_norm);
  for (int i = 0; i < arr_bt[0]; i++) {
    for (int j = 0; j < arr_bonds[1]; j++) {
      for (int k = 0; k < bins; k++) {
        int bin = id3D(i, j, k, arr_bonds),
            id_norm = id2D(i, j, arr_bonds[0]);
        bonds_norm[id_norm] += length_bonds[bin];
      }
    }
  }
  //}}}

  // write distribution of bond lengths //{{{
  PrintByline(fout, argc, argv);
  // print first line of output file - molecule names and beadtype pairs //{{{
  FILE *fw = OpenFile(fout, "a");
  fprintf(fw, "# (1) distance\n");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      fprintf(fw, "# %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          int btype1 = mt_i->BType[j],
              btype2 = mt_i->BType[k];
          int id = id3D(i, btype1, btype2, arr_bt);
          if (bt_norm[id] > 0) {
            count++;
            fprintf(fw, " (%d) %s-%s,", count, System.BeadType[btype1].Name,
                    System.BeadType[btype2].Name);
          }
        }
      }
      count++;
      fprintf(fw, " (%d)-(%d) individual bonds\n", count,
              count + mt_i->nBonds - 1);
    }
  } //}}}
  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
    for (int j = 0; j < Count->MoleculeType; j++) {
      MOLECULETYPE *mt_j = &System.MoleculeType[j];
      if (opt->mt[j]) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < mt_j->nBTypes; k++) {
          for (int l = k; l < mt_j->nBTypes; l++) {
            int btype1 = mt_j->BType[k],
                btype2 = mt_j->BType[l];
            // btype1 must be lower than btype2
            if (btype1 > btype2) {
              SwapInt(&btype1, &btype2);
            }
            int id_norm = id3D(j, btype1, btype2, arr_bt);
            if (bt_norm[id_norm] > 0) {
              int id = id4D(j, btype1, btype2, i, arr_bt);
              double value = length_bt[id] / bt_norm[id_norm];
              fprintf(fw, "%10f", value);
            }
          }
        }
        for (int k = 0; k < mt_j->nBonds; k++) {
          int id_norm = id2D(j, k, arr_bonds[0]);
          if (bonds_norm[id_norm] > 0) {
            int id = id3D(j, k, i, arr_bonds);
            double value = length_bonds[id] / bonds_norm[id_norm];
            fprintf(fw, "%10f", value);
          } else {
            fprintf(fw, "%10s", "?");
          }
        }
      }
    }
    putc('\n', fw);
  } //}}}
  // write mins, maxes, and averages //{{{
  // legend line
  fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns) -");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      fprintf(fw, " %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          int btype1 = mt_i->BType[j],
              btype2 = mt_i->BType[k];
          int id = id3D(i, btype1, btype2, arr_bt);
          if (bt_norm[id] > 0) {
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
            count += 3;
            if (k == (mt_i->nBTypes - 1)) {
              if (i != (Count->MoleculeType - 1)) {
                putc(';', fw); // molecule's last bead type pair
              }
            } else {
              putc(',', fw); // not the last pair
            }
          }
        }
      }
    }
  }
  putc('\n', fw);
  // data line
  putc('#', fw);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        int *bin_id = bins_id4D(i, j, k, arr_bt);
        // if this bin is filled, its max must be larger than 0
        if (length_bt_mma[bin_id[1]] > 0) {
          int id = id3D(i, j, k, arr_bt);
          fprintf(fw, " %lf", length_bt_mma[bin_id[0]]);
          fprintf(fw, " %lf", length_bt_mma[bin_id[1]]);
          fprintf(fw, " %lf", length_bt_mma[bin_id[2]] / bt_norm[id]);
        }
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // write distribution of distances from '-d' option //{{{
  if (opt->d_file[0] != '\0') {
    // sum up all calculated distances (normalization factors) //{{{
    int d_norm[arr_d[0]][arr_d[1]];
    for (int i = 0; i < arr_d[0]; i++) {
      if (opt->mt[i]) {
        for (int j = 0; j < arr_d[1]; j++) {
          d_norm[i][j] = 0;
          for (int k = 0; k < bins; k++) {
            int bin = id3D(i, j, k, arr_d);
            if (d_dist[bin] > 0) {
              d_norm[i][j] += d_dist[bin];
            }
          }
        }
      }
    } //}}}
    fw = OpenFile(opt->d_file, "w");
    // print command to output file
    putc('#', fw);
    PrintCommand(fw, argc, argv);
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
    fprintf(fw, "# (1) distance;");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt->mt[i]) {
        fprintf(fw, " %s:", mt_i->Name);

        for (int j = 0; j < opt->d_number; j += d_per_set) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt->d_pair[j] == -1 || opt->d_pair[j] == HIGHNUM) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = opt->d_pair[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (opt->d_pair[j+1] == -1 || opt->d_pair[j+1] == HIGHNUM) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->d_pair[j+1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
          // add semicolon for the molecule's last pair, add comma otherwise
          if (j == (opt->d_number - d_per_set)) {
            putc(';', fw);
          } else {
            putc(',', fw);
          }
        }
      }
    }
    putc('\n', fw); //}}}
    // write the distribution to output file //{{{
    for (int i = 0; i < bins; i++) {
      fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
      for (int j = 0; j < Count->MoleculeType; j++) {
        if (opt->mt[j]) {
          for (int k = 0; k < opt->d_number; k += d_per_set) {
            int id = id3D(j, k / d_per_set, i, arr_d);
            double value = (double)(d_dist[id]) / d_norm[j][(int)(k/d_per_set)];
            fprintf(fw, " %10f", value);
          }
        }
      }
      putc('\n', fw);
    } //}}}
    // write mins and maxes //{{{
    // legend line
    fprintf(fw, "# min(1st columns)/max(2nd columns)/average(3rd columns) -");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt->mt[i]) {
        fprintf(fw, " %s:", mt_i->Name);
        for (int j = 0; j < opt->d_number; j += d_per_set) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt->d_pair[j] == -1 || opt->d_pair[j] == HIGHNUM) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = opt->d_pair[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (opt->d_pair[j+1] == -1 || opt->d_pair[j+1] == HIGHNUM) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->d_pair[j+1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", count, id1, id2);
          count += 3;
          // add semicolon if this is the last pair for this molecule, add comma
          // otherwise
          if (j == (opt->d_number - 2)) {
            putc(';', fw);
          } else {
            putc(',', fw);
          }
        }
      }
    }
    putc('\n', fw);
    // data line
    putc('#', fw);
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->mt[i]) {
        for (int j = 0; j < d_pair_n; j++) {
          int *bin_id = bins_id3D(i, j, arr_d);
          // if this bin is filled, its max must be larger than 0
          if (length_bt_mma[bin_id[1]] > 0) {
            fprintf(fw, " %lf", d_mma[bin_id[0]]);
            fprintf(fw, " %lf", d_mma[bin_id[1]]);
            fprintf(fw, " %lf", d_mma[bin_id[2]] / d_norm[i][j]);
          }
        }
      }
    }
    putc('\n', fw); //}}}
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  free(opt->mt);
  free(opt);
  free(length_bt);
  free(length_bt_mma);
  free(length_bonds);
  free(length_bonds_mma);
  free(d_mma);
  free(d_dist);
  free(bt_norm);
  free(bonds_norm);
  //}}}

  return 0;
}
