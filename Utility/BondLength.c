#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
BondLength utility calculates distribution of bond lengths for all bead type \
pairs in specified molecule type(s). It can also calculate distribution of \
distances between any two beads in those molecule type(s).\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>              input coordinate file\n");
  fprintf(ptr, "   <width>              width of a single distribution bin\n");
  fprintf(ptr, "   <output>             output file with the distribution of "
               "bond lengths\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -m <name(s)>      molecules to calculate distances for"
               "(if not present, use all molecule types)");
  fprintf(ptr, "      --joined          specify that <input> contains joined "
               "coordinates\n");
  fprintf(ptr, "      -d <file> [ints]  write distribution of distances "
               "between specified beads to file <out> (if no [ints] are "
               "provided, the molecule's first and last bead are used)\n");
  fprintf(ptr, "      -w <float>        warn if the length exceeds <float> \n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool joined,       // --joined
       *mt;           // -m
  int dist_pair[100], // -d (list of bead id pairs)
      n_dist,         // -d (number of bead id pairs)
      dist_per_set;   // -d (number of beads in a set, i.e., 2 for the bond)
  double warn;        // -w option
  char d_file[LINE];  // -d (output file)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// calculate index in an 1D array that simulates a 4D one /{{{
int id4D(int i1, int i2, int i3, int i4, int size[3]) {
  return (i4 * size[0] * size[1] * size[2] +
          i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
} //}}}
// calculate index in an 1D array that simulates a 3D one /{{{
int id3D(int i1, int i2, int i3, int size[2]) {
  return (i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
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
  // --joined option
  opt->joined = BoolOption(argc, argv, "--joined"); //{{{
  if (opt->joined) { // joined coordinates supplied, so no need to join
    opt->joined = false;
  } else { // molecules need to be joined
    opt->joined = true;
  } //}}}
  // '-d' option - specify bead ids to calculate distance between //{{{
  FileIntegerOption(argc, argv, 0, 100, "-d",
                    opt->dist_pair, &opt->n_dist, opt->d_file);
  // if '-d' is present without numbers - use first and last for each molecule
  opt->dist_per_set = 2;
  if (opt->d_file[0] != '\0' && opt->n_dist == 0) {
    opt->n_dist = 2;
    opt->dist_pair[0] = 1;
    opt->dist_pair[1] = NOT; // large number to specify last bead
  }
  // Error: wrong number of integers //{{{
  if (opt->d_file[0] != '\0' && (opt->n_dist % opt->dist_per_set) != 0) {
    strcpy(ERROR_MSG, "number of bead indexes must be even");
    PrintErrorOption("-d");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < opt->n_dist; i += 2) {
    if (opt->dist_pair[i] == opt->dist_pair[i+1] ||
        opt->dist_pair[i] == 0 || opt->dist_pair[i+1] == 0) {
      strcpy(ERROR_MSG, "each pair of bead ids must be non-zero and different");
      PrintErrorOption("-d");
      exit(1);
    }
  }
  // number bead pairs for -d option
  opt->n_dist = opt->n_dist / opt->dist_per_set;
  //}}}
  for (int i = 0; i < opt->n_dist; i++) {
    opt->dist_pair[i]--; // ids should start with zero (or is -1 if none specified)
  }            //}}}
  // '-w' option - bond length warning //{{{
  if (!DoubleOption1(argc, argv, "-w", &opt->warn)) {
    opt->warn = NOT;
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

  // arrays for distributions and averages //{{{
  int arr_id[3];
  arr_id[0] = Count->MoleculeType;
  arr_id[1] = Count->BeadType;
  arr_id[2] = Count->BeadType;
  double *length = calloc(arr_id[0] * arr_id[1] * arr_id[2] * bins,
                          sizeof *length);
  double *min_max_avg = calloc(arr_id[0] * arr_id[1] * arr_id[2] * 3,
                           sizeof *min_max_avg);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        int id = id4D(i, j, k, 0, arr_id);
        min_max_avg[id] = 10 * Max3(box[0], box[1], box[2]);
      }
    }
  }
  // extra arrays for -d option
  // double *distance[Count->MoleculeType][opt->n_dist];
  int arr_id_d[2];
  arr_id_d[0] = Count->MoleculeType;
  arr_id_d[1] = opt->n_dist;
  double *dist_d = calloc(arr_id_d[0] * arr_id_d[1] * bins, sizeof *dist_d);
  // double min_max_d_option[Count->MoleculeType][opt->n_dist][2];
  double *min_max_d = calloc(arr_id_d[0] * arr_id_d[1] * 2, sizeof *min_max_d);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < opt->n_dist; j++) {
      min_max_d[id3D(i, j, 0, arr_id_d)] = 10 * Max3(box[0], box[1], box[2]);
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
    if (use) {
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, opt->joined);
      // calculate bond lengths //{{{
      // go through all molecules
      for (int i = 0; i < Count->Molecule; i++) {
        MOLECULE *mol_i = &System.Molecule[i];
        MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
        if (opt->mt[mol_i->Type]) { // use only specified molecule types
          for (int j = 0; j < mt_i->nBonds; j++) {
            // bead ids in the bond
            int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                id2 = mol_i->Bead[mt_i->Bond[j][1]];
            BEAD *b_1 = &System.Bead[id1], *b_2 = &System.Bead[id2];
            // bond length
            double bond[3];
            for (int dd = 0; dd < 3; dd++) {
              bond[dd] = b_1->Position[dd] - b_2->Position[dd];
            }
            bond[0] = VectorLength(bond);
            // warn if bond is too long //{{{
            if (opt->warn != NOT && bond[0] > opt->warn) {
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
            int bin_id[3];
            bin_id[0] = id4D(mol_i->Type, *id_lo, *id_hi, 0, arr_id);
            bin_id[1] = id4D(mol_i->Type, *id_lo, *id_hi, 1, arr_id);
            bin_id[2] = id4D(mol_i->Type, *id_lo, *id_hi, 2, arr_id);
            if (bond[0] < min_max_avg[bin_id[0]]) {
              min_max_avg[bin_id[0]] = bond[0];
            } else if (bond[0] > min_max_avg[bin_id[1]]) {
              min_max_avg[bin_id[1]] = bond[0];
            }
            min_max_avg[bin_id[2]] += bond[0];
            //}}}

            int k = bond[0] / width;
            if (k < bins) {
              length[id4D(mol_i->Type, *id_lo, *id_hi, k, arr_id)]++;
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
            for (int j = 0; j < opt->n_dist; j += opt->dist_per_set) {
              // bead ids to use //{{{
              int id1, id2;
              // use first molecule bead if bead index too high or -1
              if (opt->dist_pair[j] == -1 || opt->dist_pair[j] == NOT) {
                id1 = mol_i->Bead[0];
              } else { // use specified index otherwise
                id1 = mol_i->Bead[opt->dist_pair[j]];
              }
              // use last molecule bead if bead index too high or -1
              if (opt->dist_pair[j+1] == -1 || opt->dist_pair[j+1] == NOT) {
                id2 = mol_i->Bead[mt_i->nBeads - 1];
              } else { // use specified index otherwise
                id2 = mol_i->Bead[opt->dist_pair[j+1]];
              } //}}}
              BEAD *b_1 = &System.Bead[id1], *b_2 = &System.Bead[id2];
              // distance calculation
              double dist[3];
              for (int dd = 0; dd < 3; dd++) {
                dist[dd] = b_1->Position[dd] - b_2->Position[dd];
              }
              dist[0] = VectorLength(dist);
              // distance mins & maxes
              int pair = j / 2;
              int id[2];
              id[0] = id3D(mol_i->Type, pair, 0,arr_id_d);
              id[1] = id3D(mol_i->Type, pair, 1, arr_id_d);
              if (dist[0] < min_max_d[id[0]]) {
                min_max_d[id[0]] = dist[0];
              } else if (dist[0] > min_max_d[id[1]]) {
                min_max_d[id[1]] = dist[0];
              }
              int k = dist[0] / width; // distribution 'bin'
              if (k < bins) {
                dist_d[id3D(mol_i->Type, pair, k, arr_id_d)]++;
              }
            }
          }
        }
      } //}}}
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

  // count total number of bonds in molecules //{{{
  int bonds[Count->MoleculeType][Count->BeadType][Count->BeadType];
  // zeroize the array
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        bonds[i][j][k] = 0;
      }
    }
  }
  // count the bonds
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] += length[id4D(i, j, k, l, arr_id)];
        }
      }
    }
  } //}}}

  // write distribution of bond lengths //{{{
  PrintByline(fout, argc, argv);
  // print first line of output file - molecule names and beadtype pairs //{{{
  FILE *fw = OpenFile(fout, "w");
  fprintf(fw, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      fprintf(fw, " %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          int btype1 = mt_i->BType[j], btype2 = mt_i->BType[k];
          if (bonds[i][btype1][btype2] > 0) {
            count++;
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                    System.BeadType[btype2].Name);
            // add semicolon for molecule's last pair; add comma otherwise
            if (k == (mt_i->nBTypes - 1)) {
              putc(';', fw);
            } else {
              putc(',', fw);
            }
          }
        }
      }
    }
  }
  putc('\n', fw); //}}}
  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
    for (int j = 0; j < Count->MoleculeType; j++) {
      MOLECULETYPE *mt_j = &System.MoleculeType[j];
      if (opt->mt[j]) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < mt_j->nBTypes; k++) {
          for (int l = k; l < mt_j->nBTypes; l++) {
            int btype1 = mt_j->BType[k];
            int btype2 = mt_j->BType[l];
            // btype1 must be lower than btype2
            if (btype1 > btype2) {
              SwapInt(&btype1, &btype2);
            }
            if (bonds[j][btype1][btype2] > 0) {
              double value = length[id4D(j, btype1, btype2, i, arr_id)] /
                             bonds[j][btype1][btype2];
              fprintf(fw, "%10f", value);
            }
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
          int btype1 = mt_i->BType[j];
          int btype2 = mt_i->BType[k];
          if (bonds[i][btype1][btype2] > 0) {
            fprintf(fw, " (%d) %s-%s", count, System.BeadType[btype1].Name,
                                              System.BeadType[btype2].Name);
            count += 2;
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
        int bin_id[3];
        bin_id[0] = id4D(i, j, k, 0, arr_id);
        bin_id[1] = id4D(i, j, k, 1, arr_id);
        bin_id[2] = id4D(i, j, k, 2, arr_id);
        if (min_max_avg[bin_id[1]] > 0) {
          fprintf(fw, " %lf", min_max_avg[bin_id[0]]);
          fprintf(fw, " %lf", min_max_avg[bin_id[1]]);
          fprintf(fw, " %lf", min_max_avg[bin_id[2]] / bonds[i][j][k]);
        }
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // write distribution of distances from '-d' option //{{{
  if (opt->d_file[0] != '\0') {
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

        for (int j = 0; j < opt->n_dist; j += opt->dist_per_set) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt->dist_pair[j] == -1 || opt->dist_pair[j] == NOT) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = opt->dist_pair[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (opt->dist_pair[j+1] == -1 || opt->dist_pair[j+1] == NOT) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->dist_pair[j+1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
          // add semicolon for the molecule's last pair, add comma otherwise
          if (j == (opt->n_dist - opt->dist_per_set)) {
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
        MOLECULETYPE *mt_j = &System.MoleculeType[j];
        if (opt->mt[j]) {
          for (int k = 0; k < opt->n_dist; k += opt->dist_per_set) {
            double value = (double)(dist_d[id3D(j, k / 2, i, arr_id_d)]) /
                           (count_used * mt_j->Number);
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
        for (int j = 0; j < opt->n_dist; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (opt->dist_pair[j] == -1 || opt->dist_pair[j] == NOT) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = opt->dist_pair[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (opt->dist_pair[j+1] == -1 || opt->dist_pair[j+1] == NOT) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = opt->dist_pair[j+1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", count, id1, id2);
          count += 2;
          // add semicolon if this is the last pair for this molecule, add comma
          // otherwise
          if (j == (opt->n_dist - 2)) {
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
        for (int j = 0; j < opt->n_dist; j += opt->dist_per_set) {
          int pair = j / 2;
          fprintf(fw, " %lf %lf", min_max_d[id3D(i, pair, 0, arr_id_d)],
                                  min_max_d[id3D(i, pair, 1, arr_id_d)]);
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
  free(length);
  free(min_max_avg);
  free(min_max_d);
  free(dist_d);
  //}}}

  return 0;
}
