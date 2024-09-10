#include "../AnalysisTools.h"

// TODO: -a option to not just calculate angles between same-type angles,
//       but between all angles in the molecule; e.g., should a molecule contain
//       2 A-A-A angles, it should produce just one distribution without -a
//       (i.e., average of all A-A-A angles) but three distributions if -a is
//       present (i.e., for each of the two A-A-A angles as well as the average)
// TODO: -n option akin to BondLength's -d option

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
AngleMolecules utility calculates distribution of angles in specified \
molecule type(s) for all angles, dividing them according to different bead \
types. Note that input structure file with defined angles must be used. \
The utility can also calculate distribution of \
angles between any three beads in those molecule types (-n option).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <width> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<width>             width of a distribution bin in degrees\n");
  fprintf(ptr, "<output>            output file with the distribution of angles\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)");
  fprintf(ptr, "      --joined    specify that <input> contains joined "
          "coordinates\n");
  // see BondLength -d
  // fprintf(ptr, "  -n <ints>         bead indices (multiple of 3 <ints>)"
  //         "for angle calculation (default: 1 2 3)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join, // --joined
       *mt;  // -m
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// calculate index in an 1D array that simulates a 5D one /{{{
int id5D(int i1, int i2, int i3, int i4, int i5, int size[4]) {
  return (i5 * size[0] * size[1] * size[2] * size[3]  +
          i4 * size[0] * size[1] * size[2] +
          i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
}
// return three indices (needs 5D array size 3*size[0]*size[1]*size[2]*size[3])
int * bins_id5D(int i1, int i2, int i3, int i4, int size[4]) {
  static int bin[3];
  for (int dd = 0; dd < 3; dd++) {
    bin[dd] = id5D(i1, i2, i3, i4, dd, size);
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
  int common = 8, all = common + 2, count = 0,
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
  strcpy(option[count++], "-m");
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
  int bins = 180 / width;

  // <output> - file name with angle distribution
  char fout[LINE] = "";
  snprintf(fout, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, LINE, in);
  // --joined option
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  } //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-m <name(s)>' option
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!MoleculeTypeOption(argc, argv, "-m", true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // arrays for distributions //{{{
  int arr_id[4];
  arr_id[0] = Count->MoleculeType;
  arr_id[1] = Count->BeadType;
  arr_id[2] = Count->BeadType;
  arr_id[3] = Count->BeadType;
  int arr_id_size = arr_id[0] * arr_id[1] * arr_id[2] * arr_id[3];
  double *angle = calloc(arr_id_size * bins, sizeof *angle);
  double *min_max_avg = calloc(arr_id_size * 3, sizeof *min_max_avg);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        for (int l = 0; l < Count->BeadType; l++) {
          int id = id5D(i, j, k, l, 0, arr_id);
          min_max_avg[id] = 180; // maximum angle in degrees
        }
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
      for (int i = 0; i < Count->Molecule; i++) {
        MOLECULE *mol_i = &System.Molecule[i];
        MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
        if (opt->mt[mol_i->Type]) { // use only specified molecule types
          for (int j = 0; j < mt_i->nAngles; j++) {
            // bead ids in the bond
            int id1 = mol_i->Bead[mt_i->Angle[j][0]],
                id2 = mol_i->Bead[mt_i->Angle[j][1]],
                id3 = mol_i->Bead[mt_i->Angle[j][2]];
            BEAD *b_1 = &System.Bead[id1],
                 *b_2 = &System.Bead[id2],
                 *b_3 = &System.Bead[id3];
            // two vectors to calculate the angle for
            double u[3], v[3];
            for (int dd = 0; dd < 3; dd++) {
              u[dd] = b_1->Position[dd] - b_2->Position[dd];
              v[dd] = b_3->Position[dd] - b_2->Position[dd];
            }
            // calculate angle between the two vectors
            double size[2];
            size[0] = VectorLength(u);
            size[1] = VectorLength(v);
            double scalar = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            double ang = acos(scalar / (size[0] * size[1])); // in rad
            ang *= 180 / PI; // in degrees
            // btype1 must be lower then btype3
            int *id_lo, *id_hi;
            if (b_1->Type < b_3->Type) {
              id_lo = &b_1->Type;
              id_hi = &b_3->Type;
            } else {
              id_lo = &b_3->Type;
              id_hi = &b_1->Type;
            }

            // mins & maxes //{{{
            int *bin_id = bins_id5D(mol_i->Type, *id_lo, b_2->Type, *id_hi, arr_id);
            if (ang < min_max_avg[bin_id[0]]) {
              min_max_avg[bin_id[0]] = ang;
            } else if (ang > min_max_avg[bin_id[1]]) {
              min_max_avg[bin_id[1]] = ang;
            }
            min_max_avg[bin_id[2]] += ang;
            //}}}

            int k = ang / width;
            if (k < bins) {
              angle[id5D(mol_i->Type, *id_lo, b_2->Type, *id_hi, k, arr_id)]++;
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

  // sum up all angles in molecules (normalization factor) //{{{
  int angles[Count->MoleculeType][Count->BeadType][Count->BeadType][Count->BeadType];
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = j; l < Count->BeadType; l++) {
          angles[i][j][k][l] = 0;
          for (int m = 0; m < bins; m++) {
            int id = id5D(i, j, k, l, m, arr_id);
            if (angle[id] > 0) {
              angles[i][j][k][l] += angle[id];
            }
          }
        }
      }
    }
  } //}}}

  // write distribution of bond lengths //{{{
  PrintByline(fout, argc, argv);
  // print first line of output file - molecule names and beadtype trios //{{{
  FILE *fw = OpenFile(fout, "w");
  fprintf(fw, "# (1) distance;");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      fprintf(fw, " %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          for (int l = j; l < mt_i->nBTypes; l++) {
            int btype1 = mt_i->BType[j],
                btype2 = mt_i->BType[k],
                btype3 = mt_i->BType[l];
            if (angles[i][btype1][btype2][btype3] > 0) {
              count++;
              fprintf(fw, " (%d) %s-%s-%s", count,
                      System.BeadType[btype1].Name,
                      System.BeadType[btype2].Name,
                      System.BeadType[btype3].Name);
              // add semicolon for molecule's last trio; add comma otherwise
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
            for (int m = k; m < mt_j->nBTypes; m++) {
              int btype1 = mt_j->BType[k],
                  btype2 = mt_j->BType[l],
                  btype3 = mt_j->BType[m];
              // btype1 must be lower than btype2
              if (btype1 > btype3) {
                SwapInt(&btype1, &btype3);
              }
              if (angles[j][btype1][btype2][btype3] > 0) {
                int id = id5D(j, btype1, btype2, btype3, i, arr_id);
                double value = angle[id] / angles[j][btype1][btype2][btype3];
                fprintf(fw, "%10f", value);
              }
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
          for (int l = j; l < mt_i->nBTypes; l++) {
            int btype1 = mt_i->BType[j],
                btype2 = mt_i->BType[k],
                btype3 = mt_i->BType[l];
            if (angles[i][btype1][btype2][btype3] > 0) {
              fprintf(fw, " (%d) %s-%s-%s", count,
                      System.BeadType[btype1].Name,
                      System.BeadType[btype1].Name,
                      System.BeadType[btype3].Name);
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
  }
  putc('\n', fw);
  // data line
  putc('#', fw);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = j; l < Count->BeadType; l++) {
          int *bin_id = bins_id5D(i, j, k, l, arr_id);
          if (min_max_avg[bin_id[1]] > 0) {
            fprintf(fw, " %lf", min_max_avg[bin_id[0]]);
            fprintf(fw, " %lf", min_max_avg[bin_id[1]]);
            fprintf(fw, " %lf", min_max_avg[bin_id[2]] / angles[i][j][k][l]);
          }
        }
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  free(opt->mt);
  free(opt);
  free(angle);
  free(min_max_avg);
  //}}}

  return 0;
}
