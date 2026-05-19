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
  .args = 3, // number of mandatory arguments
  .all = 16, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
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

// structure for options //{{{
struct OPT {
  bool join,       // --joined
       *mt,          // -m
       all;          // --all
  int n_list[100],   // -n (list of bead id pairs)
      n_number;      // -n (total number of beads in the pairs)
  char n_file[LINE]; // -n (output file)
}; //}}}

void Calculation(SYSTEM *System, OPT opt, ArrNDd *ang, ArrNDd *ang_mma,
                 ArrNDd *ang_all, ArrNDd *ang_all_mma,
                 ArrNDd *ang_n, ArrNDd *ang_n_mma,
                 double width, int bins, int n_per_set) {
  COUNT *Count = &System->Count;
  WrapJoinCoordinates(System, true, opt.join);
  // calculate angles //{{{
  // go through all molecules
  // TODO: make into for (mtype); for (mtype.index)
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol_i = &System->Molecule[i];
    MOLECULETYPE *mt_i = &System->MoleculeType[mol_i->Type];
    if (opt.mt[mol_i->Type]) { // use only specified molecule types
      for (int j = 0; j < mt_i->nAngles; j++) {
        // bead ids in the angle
        int id1 = mol_i->Bead[mt_i->Angle[j][0]],
            id2 = mol_i->Bead[mt_i->Angle[j][1]],
            id3 = mol_i->Bead[mt_i->Angle[j][2]];
        BEAD *b_1 = &System->Bead[id1],
             *b_2 = &System->Bead[id2],
             *b_3 = &System->Bead[id3];
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
          // if (angle < ang_all_mma[mol_i->Type][j][0]) {
          if (angle < GetArr3D(ang_all_mma, mol_i->Type, j, 0)) {
            // ang_all_mma[mol_i->Type][j][0] = angle;
            SetArr3D(ang_all_mma, mol_i->Type, j, 0, angle);
          // } else if (angle > ang_all_mma[mol_i->Type][j][1]) {
          } else if (angle > GetArr3D(ang_all_mma, mol_i->Type, j, 1)) {
            // ang_all_mma[mol_i->Type][j][1] = angle;
            SetArr3D(ang_all_mma, mol_i->Type, j, 1, angle);
          }
          // ang_all_mma[mol_i->Type][j][2] += angle;
          AddArr3D(ang_all_mma, mol_i->Type, j, 2, angle);
        }
        //}}}

        int k = angle / width;
        if (k < bins) {
          size_t shape5D[5] = {mol_i->Type, *id_lo, b_2->Type, *id_hi, k};
          AddArrND(ang, shape5D, 1);
          if (opt.all) {
            // ang_all[mol_i->Type][j][k]++;
            AddArr3D(ang_all, mol_i->Type, j, k, 1);
          }
        }
      }
    }
  } //}}}
  // calculate extra angle (-n option) //{{{
  if (opt.n_file[0] != '\0') {
    for (int i = 0; i < Count->Molecule; i++) {
      MOLECULE *mol_i = &System->Molecule[i];
      MOLECULETYPE *mt_i = &System->MoleculeType[mol_i->Type];
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
          BEAD *b_1 = &System->Bead[id1],
               *b_2 = &System->Bead[id2],
               *b_3 = &System->Bead[id3];
          // calculate angle between the two vectors in degrees
          vec3d u = Vector(b_1->Position, b_2->Position);
          vec3d v = Vector(b_3->Position, b_2->Position);
          double angle = AngleDegrees(u, v);
          // mins & maxes & averages //{{{
          // if (angle < ang_n_mma[mol_i->Type][j/n_per_set][0]) {
          if (angle < GetArr3D(ang_n_mma, mol_i->Type, j / n_per_set, 0)) {
            // ang_n_mma[mol_i->Type][j/n_per_set][0] = angle;
            SetArr3D(ang_n_mma, mol_i->Type, j / n_per_set, 0, angle);
          // } else if (angle > ang_n_mma[mol_i->Type][j/n_per_set][1]) {
          } else if (angle > GetArr3D(ang_n_mma, mol_i->Type,
                                      j / n_per_set, 1)) {
            // ang_n_mma[mol_i->Type][j/n_per_set][1] = angle;
            SetArr3D(ang_n_mma, mol_i->Type, j / n_per_set, 1, angle);
          }
          // ang_n_mma[mol_i->Type][j/n_per_set][2] += angle;
          AddArr3D(ang_n_mma, mol_i->Type, j / n_per_set, 2, angle);
          //}}}
          int k = angle / width;
          if (k < bins) {
            // ang_n[mol_i->Type][j/n_per_set][k]++;
            AddArr3D(ang_n, mol_i->Type, j / n_per_set, k, 1);
          }
        }
      }
    }
  } //}}}
}
// structure for the callback function
struct user_data {
  OPT opt;
  ArrNDd *ang, *ang_mma,
         *ang_all, *ang_all_mma,
         *ang_n, *ang_n_mma;
  double width;
  int bins, n_per_set;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, p->opt, p->ang, p->ang_mma, p->ang_all, p->ang_all_mma,
              p->ang_n, p->ang_n_mma, p->width, p->bins, p->n_per_set);
};

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
  // --joined option (opt.join == true -> needs joining)
  opt.join = !BoolOption(argc, argv, "--joined");
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
  if (!(opt.mt = calloc(System.Count.MoleculeType, sizeof *opt.mt))) {
    ErrorAlloc("opt.mt");
  }
  if (!TypeOption(argc, argv, "-m", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // arrays for distributions //{{{
  size_t shape[5] = {Count->MoleculeType,
                     Count->BeadType,
                     Count->BeadType,
                     Count->BeadType,
                     bins};
  ArrNDd *ang = CreateArrNDd(5, shape);
  if (!ang) {
    ErrorAlloc("ang");
  }
  shape[4] = 3;
  // set maximum possible angle as initial minimum
  ArrNDd *ang_mma = CreateArrNDd(5, shape);
  if (!ang_mma) {
    ErrorAlloc("ang_mma");
  }
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
  // double ***ang_all = NULL, (**ang_all_mma)[3] = NULL;
  ArrNDd *ang_all = NULL;
  ArrNDd *ang_all_mma = NULL;
  // maximum number of angles in all molecules
  int max_angs = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (System.MoleculeType[i].nBonds > max_angs) {
      max_angs = System.MoleculeType[i].nBonds;
    }
  }
  if (opt.all) {
    if (!(ang_all = CreateArr3Dd(Count->MoleculeType, max_angs, bins)) ||
        !(ang_all_mma = CreateArr3Dd(Count->MoleculeType, max_angs, 3))) {
      ErrorAlloc("ang_all/ang_all_mma");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        // set maximum possible angle as initial minimum
        SetArr3D(ang_all_mma, i, j, 0, 180);
      }
    }
  } //}}}
  // extra arrays for -n option //{{{
  // double ***ang_n = NULL;
  // double (**ang_n_mma)[3] = NULL;
  ArrNDd *ang_n = NULL;
  ArrNDd *ang_n_mma = NULL;
  if (opt.n_file[0] != '\0') {
    // ang_n = calloc(Count->MoleculeType, sizeof *ang_n),
    // ang_n_mma = calloc(Count->MoleculeType, sizeof (**ang_n_mma)[3]);
    if (!(ang_n = CreateArr3Dd(Count->MoleculeType, n_pair_num, bins)) ||
        !(ang_n_mma = CreateArr3Dd(Count->MoleculeType, n_pair_num, 3))) {
      ErrorAlloc("ang_all/ang_all_mma");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      // ang_n[i] = calloc(n_pair_num, sizeof *ang_n);
      // ang_n_mma[i] = calloc(n_pair_num, sizeof **ang_n_mma);
      for (int j = 0; j < n_pair_num; j++) {
        // ang_n[i][j] = calloc(bins, sizeof *ang_n);
        // set maximum possible angle as initial minimum
        // ang_n_mma[i][j][0] = HIGHNUM;
        SetArr3D(ang_n_mma, i, j, 0, HIGHNUM);
      }
    }
  } //}}}

  STEP step = InitStep;
  struct user_data ud = { opt, ang, ang_mma,
                          ang_all, ang_all_mma,
                          ang_n, ang_n_mma,
                          width, bins, n_per_set };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // sum up all angles in molecules (normalization factor) //{{{
  // BeadType-BeadType-BeadType angles
  size_t shape4D[4] = {Count->MoleculeType,
                       Count->BeadType,
                       Count->BeadType,
                       Count->BeadType};
  ArrNDi *ang_norm = CreateArrNDi(4, shape4D);
  if (!ang_norm) {
    ErrorAlloc("ang_norm");
  }
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
  ArrNDi *ang_all_norm = NULL;
  if (opt.all) {
    if (!(ang_all_norm = CreateArr2Di(Count->MoleculeType, max_angs))) {
      ErrorAlloc("ang_all_norm");
    }
    // ang_all_norm = calloc(Count->MoleculeType, sizeof *ang_all_norm);
    for (int i = 0; i < Count->MoleculeType; i++) {
      // ang_all_norm[i] = calloc(System.MoleculeType[i].nAngles,
      //                          sizeof *ang_all_norm);
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        for (int k = 0; k < bins; k++) {
          // ang_all_norm[i][j] += ang_all[i][j][k];
          AddArr2D(ang_all_norm, i, j, GetArr3D(ang_all, i, j, k));
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
  // collate data //{{{
  int ncols = count;
  int nrows = bins;
  ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
  if (!data) {
    ErrorAlloc("data");
  }
  for (int i = 0; i < nrows; i++) {
    count = 0;
    SetArr2D(data, i, count++, width * (2 * i + 1) / 2);
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
                SetArr2D(data, i, count++, value);
                // fprintf(fw, "%10f", value);
              }
            }
          }
        }
        if (opt.all) {
          for (int k = 0; k < mt_j->nAngles; k++) {
            // double val = ang_all[j][k][i] / ang_all_norm[j][k];
            double val = GetArr3D(ang_all, j, k, i) /
                         GetArr2D(ang_all_norm, j, k);
            // fprintf(fw, "%10f", val);
            SetArr2D(data, i, count++, val);
          }
        }
      }
    }
  } //}}}
  ComputeColumnWidths(nrows, ncols, data, 6);
  PrintDataAll(fw, nrows, ncols, data);
  FreeArrND(data);
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
        // fprintf(fw, " %lf", ang_all_mma[i][j][0]);
        // fprintf(fw, " %lf", ang_all_mma[i][j][1]);
        // fprintf(fw, " %lf", ang_all_mma[i][j][2] / ang_all_norm[i][j]);
        fprintf(fw, " %lf", GetArr3D(ang_all_mma, i, j, 0));
        fprintf(fw, " %lf", GetArr3D(ang_all_mma, i, j, 1));
        double val = GetArr3D(ang_all_mma, i, j, 2) /
                     GetArr2D(ang_all_mma, i, j);
        fprintf(fw, " %lf", val);
      }
    }
  }
  putc('\n', fw); //}}}
  fclose(fw); //}}}

  // write distribution of angles from '-n' option //{{{
  if (opt.n_file[0] != '\0') {
    // sum up all calculated angles (normalization factors) //{{{
    // int n_norm[Count->MoleculeType][n_pair_num];
    ArrNDi *n_norm = CreateArr2Di(Count->MoleculeType, n_pair_num);
    if (!n_norm) {
      ErrorAlloc("n_norm");
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      for (int j = 0; j < n_pair_num; j++) {
        for (int k = 0; k < bins; k++) {
          // n_norm[i][j] += ang_n[i][j][k];
          AddArr2D(n_norm, i, j, GetArr3D(ang_n, i, j, k));
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
    // collate data //{{{
    int ncols = count;
    int nrows = 180;
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      ErrorAlloc("data");
    }
    for (int i = 0; i < 180; i++) {
      // fprintf(fw, "%7.4f", width * (2 * i + 1) / 2);
      count = 0;
      SetArr2D(data, i, count++, width * (2 * i + 1) / 2);
      for (int j = 0; j < Count->MoleculeType; j++) {
        if (opt.mt[j] && System.MoleculeType[j].nBeads >= n_per_set) {
          for (int k = 0; k < n_pair_num; k++) {
            // if (n_norm[j][k] > 0) {
            if (GetArr2D(n_norm, j, k) > 0) {
              // double val = (double)(ang_n[j][k][i]) / n_norm[j][k];
              double val = (double)(GetArr3D(ang_n, j, k, i)) /
                           GetArr2D(n_norm, j, k);
              // fprintf(fw, " %10f", val);
              SetArr2D(data, i, count++, val);
            }
          }
        }
      }
    } //}}}
    ComputeColumnWidths(nrows, ncols, data, 6);
    PrintDataAll(fw, nrows, ncols, data);
    FreeArrND(data);
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
          if (GetArr2D(n_norm, i, ang_id) > 0) {
            fprintf(fw, " %lf", GetArr3D(ang_n_mma, i, ang_id, 0));
            fprintf(fw, " %lf", GetArr3D(ang_n_mma, i, ang_id, 1));
            double val = GetArr3D(ang_n_mma, i, ang_id, 2) /
                         GetArr2D(n_norm, i, ang_id);
            fprintf(fw, " %lf", val);
          }
        }
      }
    } //}}}
    putc('\n', fw);
    fclose(fw);
    FreeArrND(n_norm);
  } //}}}

  // free memory - to make valgrind happy //{{{
  free(opt.mt);
  FreeArrND(ang_mma);
  // free arrays for all angles
  if (opt.all) {
    FreeArrND(ang_all);
    FreeArrND(ang_all_mma);
    FreeArrND(ang_all_norm);
  }
  FreeArrND(ang);
  FreeArrND(ang_norm);
  // free arrays for -n option
  if (opt.n_file[0] != '\0') {
    FreeArrND(ang_n);
    FreeArrND(ang_n_mma);
  }
  FreeSystem(&System);
  //}}}

  return 0;
}
