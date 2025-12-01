#include "../src/AnalysisTools.h"
// TODO: error if no angles in the system and no -n option
// TODO: -n option - why extra output file?

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AngleMolecules utility calculates distribution of angles in specified "
  "molecule type(s), dividing the angles according to different bead "
  "types. Note that input structure file with defined angles must be used. "
  "It can also add the distributions for all angles in the molecule type(s) "
  "(--all option). Finally, the utility can also calculate distribution of "
  "angles between any three beads in those molecule types (-n option).",

  "Usage: AngleMolecules <input> <width> <output> [options]\n\n",
  .args = 3,
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
  {"<input>", NULL, "input coordinate file", OPT_ARG},
  {"<width>", NULL, "width of a distribution bin in degrees", OPT_ARG},
  {"<output>", NULL, "output file with the distribution of angles", OPT_ARG},
  {"-m", "<name(s)>", "molecule types to use (default: all)", OPT_EXTRA},
  {"--joined", NULL, "specify that <input> contains joined coordinates", OPT_EXTRA},
  {"--all", NULL, "calculate distribution for all angles", OPT_EXTRA},
  {"-n", "<file> <ints>", "calculate distribution of angles between given bead trios", OPT_EXTRA},
  {NULL}
}; //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {

} //}}}

// structure for options //{{{
struct OPT {
  bool join,       // --joined
       *mt,          // -m
       all;          // --all
  int n_list[100],   // -n (list of bead id pairs)
      n_number;      // -n (total number of beads in the pairs)
  char n_file[LINE]; // -n (output file)
}; //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
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
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <output> - file name with angle distribution
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option
  if (BoolOption(argc, argv, "--joined")) {
    opt.join = false; // joined coordinates supplied, so no need to join
  } else {
    opt.join = true; // molecules need to be joined
  }
  opt.all = BoolOption(argc, argv, "--all");
  // '-n' option - specify bead ids to calculate angles between
  FileNumbersOption(argc, argv, 0, 100, "-n",
                    opt.n_list, &opt.n_number, opt.n_file, 'i');
  // if '-n' is present without numbers, use first and last for each molecule
  int n_per_set = 3; // it's an angle, so there three beads in each
  if (opt.n_file[0] != '\0' && opt.n_number == 0) {
    err_msg("missing bead indices");
    PrintErrorOption("-n");
    exit(1);
  }
  int n_pair_num = opt.n_number / n_per_set;
  // Error: wrong number of integers //{{{
  if (opt.n_file[0] != '\0' && (opt.n_number % n_per_set) != 0) {
    err_msg("number of bead indexes must a multiple of three");
    PrintErrorOption("-n");
    exit(1);
  } //}}}
  // Error: same bead ids //{{{
  for (int i = 0; i < opt.n_number; i += n_per_set) {
    if (opt.n_list[i] == opt.n_list[i+1] ||
        opt.n_list[i+1] == opt.n_list[i+2] ||
        opt.n_list[i] == opt.n_list[i+2] ||
        opt.n_list[i] == 0 ||
        opt.n_list[i+1] == 0 ||
        opt.n_list[i+2] == 0) {
      err_msg("each trio of bead ids must be non-zero and different");
      PrintErrorOption("-n");
      exit(1);
    }
  } //}}}
  //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  int bins = 180 / width;

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-m <name(s)>' option
  opt.mt = calloc(System.Count.MoleculeType, sizeof *opt.mt);
  if (!TypeOption(argc, argv, "-m", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // arrays for distributions //{{{
  // double *****ang = calloc(Count->MoleculeType, sizeof *ang);
  size_t shape[5] = {Count->MoleculeType,
                     Count->BeadType,
                     Count->BeadType,
                     Count->BeadType,
                     bins};
  ArrNDd *ang = CreateArrNDd(5, shape);
  shape[4] = 3;
  // set maximum possible angle as initial minimum
  ArrNDd *ang_mma = CreateArrNDd(5, shape);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        for (int l = 0; l < Count->BeadType; l++) {
          size_t shape5D[5] = {i, j, k, l, 0};
          SetArrND(ang_mma, shape5D, 180);
        }
      }
    }
  } //}}}
  // arrays for all angles in molecules //{{{
  double ***ang_all = NULL, (**ang_all_mma)[3] = NULL;
  if (opt.all) {
    ang_all = calloc(Count->MoleculeType, sizeof *ang_all);
    ang_all_mma = calloc(Count->MoleculeType, sizeof (**ang_all_mma)[3]);
    for (int i = 0; i < Count->MoleculeType; i++) {
      ang_all[i] = calloc(System.MoleculeType[i].nAngles, sizeof *ang_all);
      ang_all_mma[i] = calloc(System.MoleculeType[i].nAngles,
                              sizeof **ang_all_mma);
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        ang_all[i][j] = calloc(bins, sizeof *ang_all);
        // set maximum possible angle as initial minimum
        ang_all_mma[i][j][0] = 180;
      }
    }
  } //}}}
  // extra arrays for -n option //{{{
  double ***ang_n = NULL;
  double (**ang_n_mma)[3] = NULL;
  if (opt.n_file[0] != '\0') {
    ang_n = calloc(Count->MoleculeType, sizeof *ang_n),
    ang_n_mma = calloc(Count->MoleculeType, sizeof (**ang_n_mma)[3]);
    for (int i = 0; i < Count->MoleculeType; i++) {
      ang_n[i] = calloc(n_pair_num, sizeof *ang_n);
      ang_n_mma[i] = calloc(n_pair_num, sizeof **ang_n_mma);
      for (int j = 0; j < n_pair_num; j++) {
        ang_n[i][j] = calloc(bins, sizeof *ang_n);
        // set maximum possible angle as initial minimum
        ang_n_mma[i][j][0] = HIGHNUM;
      }
    }
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, opt.join);
      // calculate angles //{{{
      // go through all molecules
      // TODO: make into for (mtype); for (mtype.index)
      for (int i = 0; i < Count->Molecule; i++) {
        MOLECULE *mol_i = &System.Molecule[i];
        MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
        if (opt.mt[mol_i->Type]) { // use only specified molecule types
          for (int j = 0; j < mt_i->nAngles; j++) {
            // bead ids in the angle
            int id1 = mol_i->Bead[mt_i->Angle[j][0]],
                id2 = mol_i->Bead[mt_i->Angle[j][1]],
                id3 = mol_i->Bead[mt_i->Angle[j][2]];
            BEAD *b_1 = &System.Bead[id1],
                 *b_2 = &System.Bead[id2],
                 *b_3 = &System.Bead[id3];
            // calculate angle between the two vectors in degrees
            vec3d u = Vector(b_1->Position, b_2->Position);
            vec3d v = Vector(b_3->Position, b_2->Position);
            double angle = AngleDegrees(u, v);
            // btype1 must be lower than btype3
            int *id_lo, *id_hi;
            if (b_1->Type < b_3->Type) {
              id_lo = &b_1->Type;
              id_hi = &b_3->Type;
            } else {
              id_lo = &b_3->Type;
              id_hi = &b_1->Type;
            }

            // mins & maxes & averages //{{{
            size_t shape5D_0[5] = {mol_i->Type, *id_lo, b_2->Type, *id_hi, 0};
            size_t shape5D_1[5] = {mol_i->Type, *id_lo, b_2->Type, *id_hi, 1};
            size_t shape5D_2[5] = {mol_i->Type, *id_lo, b_2->Type, *id_hi, 2};
            if (angle < GetArrND(ang_mma, shape5D_0)) {
              SetArrND(ang_mma, shape5D_0, angle);
            } else if (angle > GetArrND(ang_mma, shape5D_1)) {
              SetArrND(ang_mma, shape5D_1, angle);
            }
            AddArrND(ang_mma, shape5D_2, angle);
            if (opt.all) {
              if (angle < ang_all_mma[mol_i->Type][j][0]) {
                ang_all_mma[mol_i->Type][j][0] = angle;
              } else if (angle > ang_all_mma[mol_i->Type][j][1]) {
                ang_all_mma[mol_i->Type][j][1] = angle;
              }
              ang_all_mma[mol_i->Type][j][2] += angle;
            }
            //}}}

            int k = angle / width;
            if (k < bins) {
              size_t shape5D[5] = {mol_i->Type, *id_lo, b_2->Type, *id_hi, k};
              AddArrND(ang, shape5D, 1);
              if (opt.all) {
                ang_all[mol_i->Type][j][k]++;
              }
            }
          }
        }
      } //}}}
      // calculate extra angle (-n option) //{{{
      if (opt.n_file[0] != '\0') {
        for (int i = 0; i < Count->Molecule; i++) {
          MOLECULE *mol_i = &System.Molecule[i];
          MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
          if (opt.mt[mol_i->Type]) { // use only specified molecule types
            for (int j = 0; j < opt.n_number; j += n_per_set) {
              // ignore the angle if any index is too high
              if (opt.n_list[j] > mt_i->nBeads ||
                  opt.n_list[j+1] > mt_i->nBeads ||
                  opt.n_list[j+2] > mt_i->nBeads) {
                continue;
              }
              // bead ids in the angle
              int id1 = mol_i->Bead[opt.n_list[j]-1],
                  id2 = mol_i->Bead[opt.n_list[j+1]-1],
                  id3 = mol_i->Bead[opt.n_list[j+2]-1];
              BEAD *b_1 = &System.Bead[id1],
                   *b_2 = &System.Bead[id2],
                   *b_3 = &System.Bead[id3];
              // calculate angle between the two vectors in degrees
              vec3d u = Vector(b_1->Position, b_2->Position);
              vec3d v = Vector(b_3->Position, b_2->Position);
              double angle = AngleDegrees(u, v);
              // mins & maxes & averages //{{{
              if (angle < ang_n_mma[mol_i->Type][j/n_per_set][0]) {
                ang_n_mma[mol_i->Type][j/n_per_set][0] = angle;
              } else if (angle > ang_n_mma[mol_i->Type][j/n_per_set][1]) {
                ang_n_mma[mol_i->Type][j/n_per_set][1] = angle;
              }
              ang_n_mma[mol_i->Type][j/n_per_set][2] += angle;
              //}}}
              int k = angle / width;
              if (k < bins) {
                ang_n[mol_i->Type][j/n_per_set][k]++;
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
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, commons.silent); //}}}

  // sum up all angles in molecules (normalization factor) //{{{
  // BeadType-BeadType-BeadType angles
  size_t shape4D[4] = {Count->MoleculeType,
                       Count->BeadType,
                       Count->BeadType,
                       Count->BeadType};
  ArrNDi *ang_norm = CreateArrNDi(4, shape4D);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        for (int l = 0; l < Count->BeadType; l++) {
          for (int m = 0; m < bins; m++) {
            size_t shape5D[5] = {i, j, k, l, m};
            if (GetArrND(ang, shape5D) > 0) {
              size_t id[4] = {i, j, k, l};
              SetArrND(ang_norm, id, GetArrND(ang, shape5D));
            }
          }
        }
      }
    }
  }
  // all molecules' angles
  int **ang_all_norm = NULL;
  if (opt.all) {
    ang_all_norm = calloc(Count->MoleculeType, sizeof *ang_all_norm);
    for (int i = 0; i < Count->MoleculeType; i++) {
      ang_all_norm[i] = calloc(System.MoleculeType[i].nAngles,
                               sizeof *ang_all_norm);
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        for (int k = 0; k < bins; k++) {
          ang_all_norm[i][j] += ang_all[i][j][k];
        }
      }
    }
  } //}}}

  // write distribution of angles //{{{
  // print first line of output file - molecule names and beadtype trios //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  fprintf(fw, "# (1) angle\n");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt.mt[i] && mt_i->nAngles > 0) {
      fprintf(fw, "# %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          for (int l = j; l < mt_i->nBTypes; l++) {
            int btype1 = mt_i->BType[j],
                btype2 = mt_i->BType[k],
                btype3 = mt_i->BType[l];
            size_t id[4] = {i, btype1, btype2, btype3};
            if (GetArrND(ang_norm, id) > 0) {
              count++;
              fprintf(fw, " (%d) %s-%s-%s", count,
                      System.BeadType[btype1].Name,
                      System.BeadType[btype2].Name,
                      System.BeadType[btype3].Name);
            }
          }
        }
      }
      if (opt.all) {
        count++;
        fprintf(fw, " (%d)-(%d) individual angles", count,
                count + mt_i->nAngles - 1);
        count += mt_i->nAngles - 1;
      }
      putc('\n', fw);
    }
  } //}}}
  // write distribution to output file //{{{
  for (int i = 0; i < bins; i++) {
    fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
    for (int j = 0; j < Count->MoleculeType; j++) {
      MOLECULETYPE *mt_j = &System.MoleculeType[j];
      if (opt.mt[j]) {
        // go over all beadtype pairs in molecule type 'j'
        for (int k = 0; k < mt_j->nBTypes; k++) {
          for (int l = k; l < mt_j->nBTypes; l++) {
            for (int m = k; m < mt_j->nBTypes; m++) {
              int btype1 = mt_j->BType[k],
                  btype2 = mt_j->BType[l],
                  btype3 = mt_j->BType[m];
              // btype1 must be lower than btype3
              if (btype1 > btype3) {
                SwapInt(&btype1, &btype3);
              }
              size_t id4D[4] = {j, btype1, btype2, btype3};
              if (GetArrND(ang_norm, id4D) > 0) {
                size_t id5D[5] = {j, btype1, btype2, btype3, i};
                double value = GetArrND(ang, id5D) / GetArrND(ang_norm, id4D);
                fprintf(fw, "%10f", value);
              }
            }
          }
        }
        if (opt.all) {
          for (int k = 0; k < mt_j->nAngles; k++) {
            if (ang_all_norm[j][k] > 0) {
              double value = ang_all[j][k][i] / ang_all_norm[j][k];
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
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt.mt[i] && mt_i->nAngles > 0) {
      fprintf(fw, "# %s molecule:", mt_i->Name);
      for (int j = 0; j < mt_i->nBTypes; j++) {
        for (int k = j; k < mt_i->nBTypes; k++) {
          for (int l = j; l < mt_i->nBTypes; l++) {
            int btype1 = mt_i->BType[j],
                btype2 = mt_i->BType[k],
                btype3 = mt_i->BType[l];
            size_t id4D[4] = {i, btype1, btype2, btype3};
            if (GetArrND(ang_norm, id4D) > 0) {
              fprintf(fw, " (%d) %s-%s-%s", count,
                      System.BeadType[btype1].Name,
                      System.BeadType[btype2].Name,
                      System.BeadType[btype3].Name);
              count += 3;
            }
          }
        }
      }
      if (opt.all) {
        fprintf(fw, " (%d)-(%d) individual angles", count,
                count + 3 * mt_i->nAngles - 1);
        count += 3 * mt_i->nAngles;
      }
      putc('\n', fw);
    }
  }
  // data line
  putc('#', fw);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = j; k < Count->BeadType; k++) {
        for (int l = j; l < Count->BeadType; l++) {
          size_t id5D_0[5] = {i, j, k, l, 0};
          size_t id5D_1[5] = {i, j, k, l, 1};
          size_t id5D_2[5] = {i, j, k, l, 2};
          if (GetArrND(ang_mma, id5D_1)) {
            fprintf(fw, " %lf", GetArrND(ang_mma, id5D_0));
            fprintf(fw, " %lf", GetArrND(ang_mma, id5D_1));
            size_t id4d[4] = {i, j, k, l};
            double value = GetArrND(ang_mma, id5D_2) / GetArrND(ang_norm, id4d);
            fprintf(fw, " %lf", value);
          }
        }
      }
    }
    if (opt.all) {
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        fprintf(fw, " %lf", ang_all_mma[i][j][0]);
        fprintf(fw, " %lf", ang_all_mma[i][j][1]);
        fprintf(fw, " %lf", ang_all_mma[i][j][2] / ang_all_norm[i][j]);
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // write distribution of angles from '-n' option //{{{
  if (opt.n_file[0] != '\0') {
    // sum up all calculated angles (normalization factors) //{{{
    int n_norm[Count->MoleculeType][n_pair_num];
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < n_pair_num; j++) {
        n_norm[i][j] = 0;
        for (int k = 0; k < bins; k++) {
          n_norm[i][j] += ang_n[i][j][k];
        }
      }
    } //}}}
    fw = PrintBylineOpenFile(opt.n_file, argc, argv);
    // print molecule names and ids with column numbers //{{{
    fprintf(fw, "# (1) angle\n");
    count = 1;
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System.MoleculeType[i];
      if (opt.mt[i] && System.MoleculeType[i].nBeads >= n_per_set) {
        fprintf(fw, "# %s:", mt_i->Name);

        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // bead ids for the angle
          int id1 = opt.n_list[j],
              id2 = opt.n_list[j+1],
              id3 = opt.n_list[j+2];
          // ignore possible angles with too high numbers
          if (id1 > mt_i->nBeads ||
              id2 > mt_i->nBeads ||
              id3 > mt_i->nBeads) {
            continue;
          }
          // the first number should be lower
          if (id1 > id3) {
            SwapInt(&id1, &id3);
          }
          fprintf(fw, " (%d) %d-%d-%d", ++count, id1, id2, id3);
        }
        putc('\n', fw);
      }
    } //}}}
    // write the distribution to output file //{{{
    for (int i = 0; i < 180; i++) {
      fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
      for (int j = 0; j < Count->MoleculeType; j++) {
        if (opt.mt[j] && System.MoleculeType[j].nBeads >= n_per_set) {
          for (int k = 0; k < n_pair_num; k++) {
            if (n_norm[j][k] > 0) {
              double value = (double)(ang_n[j][k][i]) / n_norm[j][k];
              fprintf(fw, " %10f", value);
            }
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
      if (opt.mt[i] && System.MoleculeType[i].nBeads >= n_per_set) {
        fprintf(fw, "# %s:", mt_i->Name);
        for (int j = 0; j < opt.n_number; j += n_per_set) {
          int id1 = opt.n_list[j],
              id2 = opt.n_list[j+1],
              id3 = opt.n_list[j+2];
          // skip id trios if all are too high for the molecule
          if (id1 > mt_i->nBeads ||
              id2 > mt_i->nBeads ||
              id3 > mt_i->nBeads) {
            continue;
          }
          // first number should be lower than the last one
          if (id1 > id3) {
            SwapInt(&id1, &id3);
          }
          fprintf(fw, " (%d) %d-%d-%d", count, id1, id2, id3);
          count += 3;
        }
        putc('\n', fw);
      }
    } //}}}
    // data line //{{{
    putc('#', fw);
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt.mt[i] && System.MoleculeType[i].nBeads >= n_per_set) {
        for (int j = 0; j < opt.n_number; j += n_per_set) {
          // skip id trios with all ids too high for the molecule
          if (opt.n_list[j] > System.MoleculeType[i].nBeads ||
              opt.n_list[j+1] > System.MoleculeType[i].nBeads ||
              opt.n_list[j+2] > System.MoleculeType[i].nBeads) {
            continue;
          }
          int ang_id = j / n_per_set;
          // if this bin is filled, its max must be larger than 0
          if (n_norm[i][j/n_per_set] > 0) {
            fprintf(fw, " %lf", ang_n_mma[i][ang_id][0]);
            fprintf(fw, " %lf", ang_n_mma[i][ang_id][1]);
            fprintf(fw, " %lf", ang_n_mma[i][ang_id][2] / n_norm[i][ang_id]);
          }
        }
      }
    } //}}}
    putc('\n', fw);
    fclose(fw);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(opt.mt);
  FreeArrND(ang_mma);
  // free arrays for all angles
  if (opt.all) {
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        free(ang_all[i][j]);
      }
      free(ang_all[i]);
      free(ang_all_mma[i]);
      free(ang_all_norm[i]);
    }
    free(ang_all);
    free(ang_all_mma);
    free(ang_all_norm);
  }
  FreeArrND(ang);
  FreeArrND(ang_norm);
  // free arrays for -n option
  if (opt.n_file[0] != '\0') {
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < n_pair_num; j++) {
        free(ang_n[i][j]);
      }
      free(ang_n[i]);
      free(ang_n_mma[i]);
    }
    free(ang_n);
    free(ang_n_mma);
  }
  FreeSystem(&System);
  //}}}

  return 0;
}
