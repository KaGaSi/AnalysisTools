#include "../AnalysisTools.h"

void Help(char cmd[50], bool error) { //{{{
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
  fprintf(ptr, "   <mol(s)>             molecule name(s) to calculate bond "
               "lengths for (optional and ignored if '--all' is used)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -m <name(s)>      molecules to calculate distances for"
               "(if not present, use all molecule types)");
  fprintf(ptr, "      --joined          specify that <input> contains joined "
               "coordinates\n");
  fprintf(ptr, "      -d <file> [ints]  write distribution of distances "
               "between specified beads to file <out> (if no [ints] are "
               "provided, the molecule's first and last bead are used)\n");
  fprintf(ptr, "      -w <float>        warn if the length exceeds <float> \n");
  int common = 11;
  char option[common][OPT_LENGTH];
  strcpy(option[0], "-st");
  strcpy(option[1], "-e");
  strcpy(option[2], "-sk");
  strcpy(option[3], "-i");
  strcpy(option[4], "--variable");
  strcpy(option[5], "-pbc");
  strcpy(option[6], "--detailed");
  strcpy(option[7], "-v");
  strcpy(option[8], "--silent");
  strcpy(option[9], "-h");
  strcpy(option[10], "--version");
  CommonHelp(error, common, option);
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 3; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }

  if (count < (req_args - 1) || count == (req_args - 1)) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "--all") != 0 &&
        strcmp(argv[i], "-m") != 0 && strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-d") != 0 && strcmp(argv[i], "-w") != 0 &&
        strcmp(argv[i], "-st") != 0 && strcmp(argv[i], "-e") != 0 &&
        strcmp(argv[i], "-sk") != 0 && strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "--variable") != 0 && strcmp(argv[i], "-pbc") != 0 &&
        strcmp(argv[i], "--detailed") != 0 && strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 && strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type, struct_file,
                       &struct_type)) {
    exit(1);
  } //}}}

  // <width> - width of a single bin //{{{
  // Error - non-numeric argument
  double width = -1;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - file name with bond length distribution //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]); //}}}

  // options before reading system data
  bool silent, verbose, detailed, vtf_var;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
  if (joined) { // --joined means it is joined, but the opposite is needed
    joined = false;
  } else {
    joined = true;
  }

  if (!silent) {
    PrintCommand(stdout, argc, argv);
  }

  int ltrj_start_id = -1; // for lammpstrj structure file, start ids from 0 or 1
  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz, &ltrj_start_id);

  // '-m <name(s)>' option
  bool *use_moltype = calloc(System.Count.MoleculeType, sizeof *use_moltype);
  if (MoleculeTypeOption(argc, argv, "-m", true, use_moltype, System)) {
    exit(1);
  }

  // '-d' option - specify bead ids to calculate distance between //{{{
  int bead[100] = {0},
      number_of_beads,   // number of parameters to -d
      beads_per_set = 2; // the numbers must come in pairs
  char output_d[LINE] = "";
  if (FileIntegerOption(argc, argv, 100, "-d", bead,
                        &number_of_beads, output_d)) {
    exit(1);
  }
  // if '-d' is present, but without numbers - use first and last for each
  // molecule
  if (output_d[0] != '\0' && number_of_beads == 0) {
    number_of_beads = 2;
    bead[0] = 1;
    bead[1] = 1000000;
  }
  // Error: wrong number of integers //{{{
  if (output_d[0] != '\0' && (number_of_beads % beads_per_set) != 0) {
    strcpy(ERROR_MSG, "number of bead indexes must be even");
    PrintErrorOption("-d");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < number_of_beads; i += 2) {
    if (bead[i] == bead[i + 1] || bead[i] == 0 || bead[i + 1] == 0) {
      strcpy(ERROR_MSG, "each pair of bead ids must be non-zero and different");
      PrintErrorOption("-d");
      exit(1);
    }
  } //}}}
  for (int i = 0; i < number_of_beads; i++) {
    bead[i]--; // ids should start with zero (or is -1 if none specified)
  }            //}}}

  int number_of_pairs = number_of_beads / beads_per_set;

  if (verbose) {
    VerboseOutput(System);
  }

  // '-w' option - bond length warning //{{{
  double warn[1] = {-1};
  int trash;
  if (DoubleOption(argc, argv, 1, "-w", &trash, warn)) {
    exit(1);
  } //}}}

  // number of bins
  VECTOR *box = &System.Box.Length;
  int bins = Max3(box->x, box->y, box->z) / (2 * width);

  // arrays for distributions //{{{
  COUNT *Count = &System.Count;
  double *length[Count->MoleculeType][Count->BeadType][Count->BeadType];
  double min_max[Count->MoleculeType][Count->BeadType][Count->BeadType][2];
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        length[i][j][k] = calloc(bins, sizeof *length[i][j][k]);
        min_max[i][j][k][0] = 10 * Max3(box->x, box->y, box->z);
        min_max[i][j][k][1] = 0;
      }
    }
  }
  // extra arrays for -d option
  double *distance[Count->MoleculeType][number_of_pairs];
  double min_max_d_option[Count->MoleculeType][number_of_pairs][2];
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      distance[i][j] = calloc(bins, sizeof *distance[i][j]);
      min_max_d_option[i][j][0] = 10 * Max3(box->x, box->y, box->z);
      min_max_d_option[i][j][1] = 0;
    }
  } //}}}

  // open input coordinate file
  FILE *coor = OpenFile(coor_file, "r");

  // main loop //{{{
  int count_coor = 0,
      count_saved = 0,
      line_count = 0;
  while (true) {
    count_coor++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (count_coor < start) {
        fprintf(stdout, "\rDiscarding step: %d", count_coor);
      } else {
        if (count_coor == start) {
          fprintf(stdout, "\rStarting step: %d    \n", start);
        }
        fprintf(stdout, "\rStep: %d", count_coor);
      }
      fflush(stdout);
    } //}}}
    bool use = false;
    // use every skip-th timestep between start and end
    if (count_coor >= start && (count_coor <= end || end == -1) &&
        ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(coor_type, coor, coor_file, &System, &line_count,
                        ltrj_start_id, vtf_var)) {
        count_coor--;
        break;
      }
      count_saved++;
      WrapJoinCoordinates(&System, false, joined);
      // calculate bond lengths //{{{
      // go through all molecules
      for (int i = 0; i < Count->Molecule; i++) {
        MOLECULE *mol_i = &System.Molecule[i];
        MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
        if (use_moltype[mol_i->Type]) { // use only specified molecule types
          for (int j = 0; j < mt_i->nBonds; j++) {

            // bead ids in the bond
            int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                id2 = mol_i->Bead[mt_i->Bond[j][1]];
            BEAD *b_1 = &System.Bead[id1],
                 *b_2 = &System.Bead[id2];
            // bond length //{{{
            VECTOR bond;
            bond.x = b_1->Position.x - b_2->Position.x;
            bond.y = b_1->Position.y - b_2->Position.y;
            bond.z = b_1->Position.z - b_2->Position.z;
            bond.x = VectorLength(bond); //}}}
            // warn if bond is too long //{{{
            if (warn[0] > -1 && bond.x > warn[0]) {
              snprintf(ERROR_MSG, LINE, "-w option; too long a bond between "
                       "beads %s%d%s and %s%d%s (%s%lf%s)", ErrYellow(), id1,
                       ErrCyan(), ErrYellow(), id2, ErrCyan(), ErrYellow(),
                       bond.x, ErrCyan());
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

            // mins & maxes //{{{
            if (bond.x < min_max[mol_i->Type][*id_lo][*id_hi][0]) {
              min_max[mol_i->Type][*id_lo][*id_hi][0] = bond.x;
            } else if (bond.x > min_max[mol_i->Type][*id_lo][*id_hi][1]) {
              min_max[mol_i->Type][*id_lo][*id_hi][1] = bond.x;
            } //}}}

            int k = bond.x / width;
            if (k < bins) {
              length[mol_i->Type][*id_lo][*id_hi][k]++;
            }
          }
        }
      } //}}}
      // calculate distance (-d option) //{{{
      if (output_d[0] != '\0') {
        for (int i = 0; i < Count->Molecule; i++) {
          MOLECULE *mol_i = &System.Molecule[i];
          MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
          if (use_moltype[mol_i->Type]) { // use only specified molecule types
            for (int j = 0; j < number_of_beads; j += beads_per_set) {
              // bead ids to use //{{{
              int id1, id2;
              // use first molecule bead if bead index too high or -1
              if (bead[j] == -1 || bead[j] >= mt_i->nBeads) {
                id1 = mol_i->Bead[0];
              } else { // use specified index otherwise
                id1 = mol_i->Bead[bead[j]];
              }
              // use last molecule bead if bead index too high or -1
              if (bead[j+1] == -1 ||
                  bead[j+1] >= mt_i->nBeads) {
                id2 = mol_i->Bead[mt_i->nBeads-1];
              } else { // use specified index otherwise
                id2 = mol_i->Bead[bead[j+1]];
              } //}}}
              BEAD *b_1 = &System.Bead[id1],
                   *b_2 = &System.Bead[id2];
              // distance calculation //{{{
              VECTOR dist;
              dist.x = b_1->Position.x - b_2->Position.x;
              dist.y = b_1->Position.y - b_2->Position.y;
              dist.z = b_1->Position.z - b_2->Position.z;
              dist.x = VectorLength(dist); //}}}
              // distance mins & maxes //{{{
              int pair = j / 2;
              if (dist.x < min_max_d_option[mol_i->Type][pair][0]) {
                min_max_d_option[mol_i->Type][pair][0] = dist.x;
              } else if (dist.x > min_max_d_option[mol_i->Type][pair][1]) {
                min_max_d_option[mol_i->Type][pair][1] = dist.x;
              } //}}}
              int k = dist.x / width; // distribution 'bin'
              if (k < bins) {
                distance[mol_i->Type][pair][k]++;
              }
            }
          }
        }
      } //}}}
    }
    if (count_coor == end) {
      break;
    }
  }
  fclose(coor);
  // print last step?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_coor);
  } //}}}

  // count total number of bonds in molecules //{{{
  int bonds[Count->MoleculeType][Count->BeadType][Count->BeadType];
  // zeroize the array
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] = 0;
        }
      }
    }
  }
  // count the bonds
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = 0; l < bins; l++) {
          bonds[i][j][k] += length[i][j][k][l];
        }
      }
    }
  } //}}}

  // write distribution of bond lengths
  FILE *fw = OpenFile(output, "w");
  PrintByline(fw, argc, argv);
  // print first line of output file - molecule names and beadtype pairs //{{{
  fprintf(fw, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (use_moltype[i]) {
      fprintf(fw, " %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          int btype1 = mt_i->BType[j],
              btype2 = mt_i->BType[k];
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
      if (use_moltype[j]) {
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
              double value = length[j][btype1][btype2][i] /
                             bonds[j][btype1][btype2];
              fprintf(fw, "%10f", value);
            }
          }
        }
      }
    }
    putc('\n', fw);
  } //}}}
  // write mins and maxes //{{{
  // legend line
  fprintf(fw, "# mins(odd columns)/maxes(even columns) -");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (use_moltype[i]) {
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
        if (min_max[i][j][k][1] > 0) {
          fprintf(fw, " %lf %lf", min_max[i][j][k][0], min_max[i][j][k][1]);
        }
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw);

  // write distribution of distances from '-d' option
  if (output_d[0] != '\0') {
    fw = OpenFile(output_d, "w");
    // print command to output file
    putc('#', fw);
    PrintCommand(fw, argc, argv);
    // print the first line - molecule names with bead order //{{{
    fprintf(fw, "# bead order in molecule(s) -");
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (use_moltype[i]) {
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
      if (use_moltype[i]) {
        fprintf(fw, " %s:", mt_i->Name);

        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= mt_i->nBeads) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = bead[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j + 1] == -1 || bead[j + 1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = bead[j + 1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", ++count, id1, id2);
          // add semicolon for the molecule's last pair, add comma otherwise
          if (j == (number_of_beads - beads_per_set)) {
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
        if (use_moltype[j]) {
          for (int k = 0; k < number_of_beads; k += beads_per_set) {
            double value = (double)(distance[j][k / 2][i]) /
                           (count_saved * mt_j->Number);
            fprintf(fw, " %10f", value);
          }
        }
      }
      putc('\n', fw);
    } //}}}
    // write mins and maxes //{{{
    // legend line
    fprintf(fw, "# mins(odd columns)/maxes(even columns) -");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (use_moltype[i]) {
        fprintf(fw, " %s:", mt_i->Name);
        for (int j = 0; j < number_of_beads; j += 2) {
          // bead ids the distance //{{{
          int id1, id2;
          // use first molecule bead if bead index too high or -1
          if (bead[j] == -1 || bead[j] >= mt_i->nBeads) {
            id1 = 1;
          } else { // use specified index otherwise
            id1 = bead[j] + 1;
          }
          // use last molecule bead if bead index too high or -1
          if (bead[j+1] == -1 || bead[j+1] >= mt_i->nBeads) {
            id2 = mt_i->nBeads;
          } else { // use specified index otherwise
            id2 = bead[j+1] + 1;
          } //}}}

          fprintf(fw, " (%d) %d-%d", count, id1, id2);
          count += 2;
          // add semicolon if this is the last pair for this molecule, add comma
          // otherwise
          if (j == (number_of_beads - 2)) {
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
      if (use_moltype[i]) {
        for (int j = 0; j < number_of_beads; j += beads_per_set) {
          int pair = j / 2;
          fprintf(fw, " %lf %lf", min_max_d_option[i][pair][0],
                                   min_max_d_option[i][pair][1]);
        }
      }
    }
    putc('\n', fw); //}}}
    fclose(fw);
  }

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        free(length[i][j][k]);
      }
    }
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < number_of_pairs; j++) {
      free(distance[i][j]);
    }
  } //}}}

  return 0;
}
