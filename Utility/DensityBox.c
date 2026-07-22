#include "../src/AnalysisTools.h"
// TODO: possible changing box size: make bins' width variable, keeping their
//       number, and work in relative coordinates (relative to instantaneous
//       dimensions) throughout the code

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "DensityBox utility calculates number density for all bead types in the "
  "direction of all axes (x, y, and z). The utility works properly only for "
  "orthogonal boxes that do not change size.",

  "Usage: DensityBox <input> <width> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 14, // number of valid lines OptSpec (not counting last {nullptr})
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
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<width>", nullptr, "width of a single bin", OPT_ARG},
  {"<output>", nullptr, "3 output files (appends -<axis>.txt)", OPT_ARG},
  {"-x", "<name(s)>", "exclude specified molecule type(s)", OPT_EXTRA},
  {"--per-bead", nullptr, "calculate density per bead position within each "
    "molecule type (instead of per bead type)", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool *x;       // -x option
  bool per_bead; // --per-bead option
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
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <outputt> - filename
  char fout_rho[LINE] = "";
  s_strcpy(fout_rho, argv[++count], LINE);
  fout_rho[LINE-7] = '\0'; // for adding -<axis>.txt
  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  opt.per_bead = BoolOption(argc, argv, "--per-bead");

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  BOX *box = &System.Box;

  // -x option
  if (!(opt.x = calloc(Count->MoleculeType, sizeof *opt.x))) {
    ErrorAlloc("opt.x");
  }
  InitBoolArray(opt.x, Count->MoleculeType, true);
  TypeOption(argc, argv, "-x", 'm', false, opt.x, System);

  // --per-bead: precompute per-molecule-type offsets and per-bead position index //{{{
  int total_pos = 0;
  int *offset = nullptr;   // offset[mt] = first column index for molecule type mt
  int *bead_pos = nullptr; // bead_pos[bead_id] = index of bead within its molecule
  if (opt.per_bead) {
    if (!(offset = calloc(Count->MoleculeType + 1, sizeof *offset))) {
      ErrorAlloc("offset");
    }
    for (int mt = 0; mt < Count->MoleculeType; mt++) {
      offset[mt] = total_pos;
      total_pos += System.MoleculeType[mt].nBeads;
    }
    offset[Count->MoleculeType] = total_pos;
    if (!(bead_pos = calloc(Count->Bead, sizeof *bead_pos))) {
      ErrorAlloc("bead_pos");
    }
    for (int m = 0; m < Count->Molecule; m++) {
      int mtype = System.Molecule[m].Type;
      for (int b = 0; b < System.MoleculeType[mtype].nBeads; b++) {
        bead_pos[System.Molecule[m].Bead[b]] = b;
      }
    }
  }
  int n_cols = Count->BeadType;
  if (opt.per_bead) {
    n_cols = total_pos;
  } //}}}

  // number of bins //{{{
  if (box->Volume == -1) {
    err_msg("missing box dimensions");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }
  double bin[3];
  // TODO: *3 to assume box change of at most thrice as big
  //       probably change from width to number of bins per box?
  for (int dd = 0; dd < 3; dd++) {
    bin[dd] = ceil(box->Length.v[dd] / width) * 3;
  } //}}}
  int bin_max = Max3(bin[0], bin[1], bin[2]);

  bool *n_used = calloc(n_cols, sizeof *n_used);
  ArrNDli *rho = CreateArr3Dli(3, n_cols, bin_max);
  if (!n_used || !rho) {
    ErrorAlloc("n_used/rho");
  }

  if (commons.verbose) {
    VerboseOutput(System);
  }

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
    if (use) {
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);

      ArrNDi *rho_temp = CreateArr3Di(3, n_cols, bin_max);
      if (!rho_temp) {
        ErrorAlloc("rho_temp");
      }

      // calculate densities //{{{
      for (int i = 0; i < Count->BeadCoor; i++) {
        use = true;
        int id = System.BeadCoor[i];
        BEAD *bead = &System.Bead[id];
        int mol = bead->Molecule;
        if (mol != -1) { // do not use excluded molecules (-x option)
          int mtype = System.Molecule[mol].Type;
          use = opt.x[mtype];
        }
        if (use) {
          int col;
          if (opt.per_bead) {
            if (mol == -1) {
              continue; // skip unbonded beads in --per-bead mode
            }
            int mtype = System.Molecule[mol].Type;
            col = offset[mtype] + bead_pos[id];
          } else {
            col = bead->Type;
          }
          n_used[col] = true;
          for (int dd = 0; dd < 3; dd++) {
            int j = bead->Position.v[dd] / width;
            AddArr3D(rho_temp, dd, col, j, 1);
          }
        }
      } //}}}
      // add from temporary density arrays to global density arrays
      for (int j = 0; j < n_cols; j++) {
        for (int dd = 0; dd < 3; dd++) {
          for (int k = 0; k < bin[dd]-1; k++) {
            AddArr3D(rho, dd, j, k, GetArr3D(rho_temp, dd, j, k));
          }
        }
      }

      FreeArrND(rho_temp);
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, commons.silent); //}}}

  // write densities to output file(s) //{{{
  for (int ax = 0; ax < 3; ax++) {
    // axis-based variables
    double volume = width;
    // TODO: redo when box-changing is reasonably dealt with
    int bins = 0, // number of bins
        size = -1;
    char axis;
    if (ax == 0) {
      axis = 'x';
      size = box->Length.x;
      volume *= box->Length.y * box->Length.z;
      bins = bin[0];
    } else if (ax == 1) {
      axis = 'y';
      size = box->Length.y;
      volume *= box->Length.x * box->Length.z;
      bins = bin[1];
    } else {
      axis = 'z';
      size = box->Length.z;
      volume *= box->Length.x * box->Length.y;
      bins = bin[2];
    }
    char file[LINE]; // filename <output>-<axis>.rho
    if (snprintf(file, LINE, "%s-%c.txt", fout_rho, axis) < 0) {
      ErrorSnprintf();
    }
    // write initial stuff to output density file
    FILE *fw = PrintBylineOpenFile(file, argc, argv);
    // print column headers to output file
    fprintf(fw, "# columns: (1) distance");
    count = 1;
    if (opt.per_bead) {
      for (int mt = 0; mt < Count->MoleculeType; mt++) {
        for (int b = 0; b < System.MoleculeType[mt].nBeads; b++) {
          int col = offset[mt] + b;
          if (n_used[col]) {
            count++;
            fprintf(fw, "; (%d) %s:%d", count, System.MoleculeType[mt].Name, b + 1);
          }
        }
      }
    } else {
      for (int i = 0; i < Count->BeadType; i++) {
        if (n_used[i]) {
          count++;
          fprintf(fw, "; (%d) %s", count, System.BeadType[i].Name);
        }
      }
    }
    putc('\n', fw);
    // collate data
    int ncols = count;
    int nrows;
    // TODO: redo when box-changing is reasonably dealt with
    for (nrows = 0; nrows < bins; nrows++) {
      double dist = width * (2 * nrows + 1) / 2;
      if (dist > size) { // write only til the max box size
        break;
      }
    }
    ArrNDd *data = CreateArr2Dd(nrows + 2, ncols);
    if (!data) {
      err_msg("ArrNDd constructor failed (data)");
      PrintError();
      exit(1);
    }
    for (int i = 0; i < nrows; i++) {
      double dist = width * (2 * i + 1) / 2;
      count = 0;
      SetArr2D(data, i, count++, dist);
      for (int j = 0; j < n_cols; j++) {
        if (n_used[j]) {
          double rho_val = GetArr3D(rho, ax, j, i) / (volume * count_used);
          SetArr2D(data, i, count++, rho_val);
        }
      }
    }
    ComputeColumnWidths(nrows, ncols, data, 6);
    PrintDataAll(fw, nrows, ncols, data);
    FreeArrND(data);
    fclose(fw);
  } //}}}

  FreeSystem(&System);
  FreeArrND(rho);
  free(n_used);
  free(opt.x);
  free(offset);
  free(bead_pos);

  return 0;
}
