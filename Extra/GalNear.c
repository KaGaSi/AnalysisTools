#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "For each timestep, calculate the fraction of 'C' beads and 'CA2' "
  "molecules that are within <dist> of at least one bead of the specified "
  "types (-bt) in the specified molecule types (-mt). The denominator is the "
  "number of 'C' beads / 'CA2' molecules present in the given timestep.",

  "Usage: GalNear <input> <output> <dist> [options]",
  .args = 3,
  .all = 13,
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
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<output>", nullptr, "output file", OPT_ARG},
  {"<dist>", nullptr, "maximum contact distance", OPT_EXTRA},
  {"-mt", "<name(s)>","use specified molecule type(s)", OPT_EXTRA},
  {"-bt", "<name(s)>","use specified bead type(s)", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  bool *mt; // -mt
  bool *bt; // -bt
}; //}}}

// name for the monovalent cation bead
static char *name = "C";
// name for the divalent cation molecule
static char *name_mol = "CA2";

int main(int argc, char *argv[]) {

  // command line arguments //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input>
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <output>
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <dist>
  double dist_check = 0;
  if (!IsPosRealNumber(argv[++count], &dist_check)) {
    ErrorNaN("<dist>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  COMMON_OPT commons = CommonOptions(argc, argv, in); //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  const BOX *boxlength = &System.Box;

  // define variables for mono- and divalent counterions //{{{
  const int bt_name = FindBeadType(name, System);
  const int mt_name = FindMoleculeName(name_mol, System);
  MOLECULETYPE *MolType_name = nullptr;
  BEADTYPE *BType_name = nullptr;
  if (mt_name != -1) {
    MolType_name = &System.MoleculeType[mt_name];
  }
  if (bt_name != -1) {
    BType_name = &System.BeadType[bt_name];
  } //}}}

  // molecule/bead type options //{{{
  opt.mt = calloc(Count->MoleculeType, sizeof *opt.mt);
  if (!TypeOption(argc, argv, "-mt", 'm', true, opt.mt, System)) {
    InitBoolArray(opt.mt, Count->MoleculeType, true);
  }
  opt.bt = calloc(Count->BeadType, sizeof *opt.bt);
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
    InitBoolArray(opt.bt, Count->BeadType, true);
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  // print initial stuff to output file //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  int col = 1;
  fprintf(fw, "# (%d) step", col++);
  if (bt_name != -1) {
    fprintf(fw, ", (%d) fraction of '%s' beads near -mt/-bt beads", col++, name);
  }
  if (mt_name != -1) {
    fprintf(fw, ", (%d) fraction of '%s' molecules near -mt/-bt beads",
            col++, name_mol);
  }
  putc('\n', fw);
  fclose(fw); //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0,
      count_used = 0,
      line_count = 0;
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    bool use = UseStep(commons, count_coor);
    if (use) {
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;

      // build list of target beads: -bt types in -mt molecules //{{{
      int *target = calloc(Count->BeadCoor, sizeof *target);
      int n_target = 0;
      for (int i = 0; i < Count->BondedCoor; i++) {
        int id_i = System.Bonded[i];
        BEAD *b_i = &System.Bead[id_i];
        MOLECULE *m_i = &System.Molecule[b_i->Molecule];
        if (opt.mt[m_i->Type] && opt.bt[b_i->Type]) {
          target[n_target++] = id_i;
        }
      } //}}}

      // fraction of 'name' beads near any target bead //{{{
      double frac_name = 0;
      if (bt_name != -1 && BType_name->InCoor > 0) {
        int near = 0;
        for (int i = 0; i < BType_name->InCoor; i++) {
          int id_i = BType_name->Index[i];
          BEAD *b_i = &System.Bead[id_i];
          for (int j = 0; j < n_target; j++) {
            BEAD *b_j = &System.Bead[target[j]];
            vec3d d = DistancePBC(b_i->Position, b_j->Position, boxlength);
            if (VectLength(d) <= dist_check) {
              near++;
              break;
            }
          }
        }
        frac_name = (double)near / BType_name->InCoor;
      } //}}}

      // fraction of 'name_mol' molecules near any target bead //{{{
      double frac_name_mol = 0;
      if (mt_name != -1) {
        int near = 0, in_ts = 0;
        for (int i = 0; i < MolType_name->Number; i++) {
          MOLECULE *mol = &System.Molecule[MolType_name->Index[i]];
          if (!mol->InTimestep) {
            continue;
          }
          in_ts++;
          bool mol_near = false;
          for (int k = 0; k < MolType_name->nBeads && !mol_near; k++) {
            BEAD *b_k = &System.Bead[mol->Bead[k]];
            if (!b_k->InTimestep) {
              continue;
            }
            for (int j = 0; j < n_target && !mol_near; j++) {
              BEAD *b_j = &System.Bead[target[j]];
              vec3d d = DistancePBC(b_k->Position, b_j->Position, boxlength);
              if (VectLength(d) <= dist_check) {
                mol_near = true;
              }
            }
          }
          if (mol_near) {
            near++;
          }
        }
        if (in_ts > 0) {
          frac_name_mol = (double)near / in_ts;
        }
      } //}}}

      // write to output //{{{
      fw = OpenFile(fout, "a");
      fprintf(fw, "%5d", count_used);
      if (bt_name != -1) {
        fprintf(fw, " %lf", frac_name);
      }
      if (mt_name != -1) {
        fprintf(fw, " %lf", frac_name_mol);
      }
      putc('\n', fw);
      fclose(fw); //}}}

      free(target);
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
  if (!commons.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // free memory //{{{
  free(opt.mt);
  free(opt.bt);
  FreeSystem(&System);
  //}}}

  return 0;
}
