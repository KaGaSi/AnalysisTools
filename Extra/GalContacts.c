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
I DON'THINK THIS IS CURRENT, RIGHT? ...SEEMS TO COUNT ONLY 3-BODY CONTACTS... \
Count contacts intra- and intermolecular contacts between specified bead types \
in specified molecule types. Also counts 3-body contacts, specifically, \
two of the specified bead types with 'C' bead type or 'CA2' molecule type's \
geometric centre (these names are hardcoded).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> <skip> <dist> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "<skip>              number of in-between beads to skip\n");
  fprintf(ptr, "<dist>              minimum distance contact check\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -mt <name(s)>     use specified molecule type(s)\n");
  fprintf(ptr, "  -bt <name(s)>     use specified bead type(s)\n");
  fprintf(ptr, "  ---multi          allow multiple trios with the same ion\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  bool *mt,          // -mt
       *bt,          // -bt
       multi;        // --multi
  int int1;          // -int1
  int int2[2];       // -int2
  char f_file[LINE]; // -f (filename)
  int f_list[100],   // -f (list of numbers)
      f_num;         // -f (number of those numbers)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// name for the monovalent cation
static char *name = "C";
// name for the divalent cation (molecule with two connected beads)
static char *name_mol = "CA2";
// calculate distance between two points, accounting for pbc //{{{
inline static double DistLength(const double v1[3], const double v2[3],
                                const double box[3]) {
  vec3 dist = Distance(v1, v2, box);
  return VectLength(dist);
} //}}}
// get two molecule/bead types ordered so the first one < second one //{{{
typedef struct {
  int a, b;
} Types;
inline static Types SortTypes(const int type1, const int type2) {
  Types tp;
  tp.a = type1;
  tp.b = type2;
  if (tp.a > tp.b) {
    SwapInt(&tp.a, &tp.b);
  }
  return tp;
} //}}}

// get bead position in molecule //{{{
int PosInMol(const int b_id, const SYSTEM System) {
  int m_id = System.Bead[b_id].Molecule;
  int mtype = System.Molecule[m_id].Type;
  for (int i = 0; i < System.MoleculeType[mtype].nBeads; i++) {
    if (b_id == System.Molecule[m_id].Bead[i]) {
      return i;
    }
  }
  err_msg("bead not in the molecule!");
  PrintError();
  exit(1);
} //}}}

// heavily ChatGPT-advised code to not have to repeat the loops //{{{
// Context for print_bead_types
typedef struct {
  int column;
} PrintHeader_ctx;
// Context for compute_averages
typedef struct {
  ArrNDi *intra_step;
  int mt_name, // molecule type id of the divalent counterion 'name_mol'
      bt_name, // bead type id of the monovalent counterion 'name'
      mt; // molecule name (the 0_3G500 or some such)
} PrintAvgContacts_ctx;
// Define a function pointer type for loop body functions
typedef void (*BodyFunc)(int, int, FILE *, SYSTEM, OPT *, void *);
// Function to iterate over bead types
void iterate_btypes(int i, FILE *fw, SYSTEM Sys, OPT *opt,
                    BodyFunc func, void *context) {
  for (int j = 0; j < Sys.MoleculeType[i].nBTypes; j++) {
    for (int k = j; k < Sys.MoleculeType[i].nBTypes; k++) {
      int bt_j = Sys.MoleculeType[i].BType[j];
      int bt_k = Sys.MoleculeType[i].BType[k];
      if (opt->bt[bt_j] && opt->bt[bt_k]) {
        func(bt_j, bt_k, fw, Sys, opt, context);
      }
    }
  }
}
// Function to print bead types
void PrintHeader(int bt_j, int bt_k, FILE *fw,
                 SYSTEM Sys, OPT *opt, void *context) {
  PrintHeader_ctx *ctx = (PrintHeader_ctx *)context;
  char *name_j = Sys.BeadType[bt_j].Name;
  char *name_k = Sys.BeadType[bt_k].Name;
  if (FindBeadType(name, Sys) != -1) {
    fprintf(fw, " (%d) %s-%s-%s;", ctx->column++, name_j, name_k, name);
  }
  if (FindMoleculeName(name_mol, Sys) != -1) {
    fprintf(fw, " (%d) %s-%s-%s", ctx->column++, name_j, name_k, name_mol);
  }
}
// Function to compute averages
void PrintAvgContacts(int bt_j, int bt_k, FILE *fw,
                      SYSTEM Sys, OPT *opt, void *context) {
  PrintAvgContacts_ctx *ctx = (PrintAvgContacts_ctx *)context;
  if (ctx->intra_step) {
    int num_mol = Sys.MoleculeType[ctx->mt].Number;
    if (ctx->mt_name != -1) {
      size_t id4[4] = {ctx->mt, bt_j, bt_k, 0};
      double avg = (double)(GetArrND(ctx->intra_step, id4)) / num_mol;
      fprintf(fw, " %lf", avg);
    }
    if (ctx->bt_name != -1) {
      size_t id4[4] = {ctx->mt, bt_j, bt_k, 1};
      double avg = (double)(GetArrND(ctx->intra_step, id4)) / num_mol;
      fprintf(fw, " %lf", avg);
    }
  }
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity //{{{
  int common = 8, all = common + 3, column = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "-mt", "-bt", "--multi"); //}}}

  column = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++column], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++column], LINE);
  // <skip> - how many beads to skip at least betweem contact-able beads
  long skip = 0;
  if (!IsWholeNumber(argv[++column], &skip)) {
    ErrorNaN("<skip>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  }
  double dist_check = 0;
  if (!IsPosRealNumber(argv[++column], &dist_check)) {
    ErrorNaN("<dist>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  } //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);
  opt->multi = BoolOption(argc, argv, "--multi");

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  double *boxlength = System.Box.Length;

  // define variables for mono- and divalent counterions //{{{
  const int bt_name = FindBeadType(name, System);
  const int mt_name = FindMoleculeName(name_mol, System);
  // int bt_name_mol = -1;
  MOLECULETYPE *MolType_name = NULL;
  BEADTYPE *BType_name = NULL;
  if (mt_name != -1) {
    // bt_name_mol = System.MoleculeType[mt_name].BType[0];
    MolType_name = &System.MoleculeType[mt_name];
  }
  if (bt_name != -1) {
    BType_name = &System.BeadType[bt_name];
  } //}}}

  // molecule/bead type options //{{{
  // molecule types to calculate contacts for
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!TypeOption(argc, argv, "-mt", 'm', true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }
  // bead types to calculate contacts for
  opt->bt = calloc(System.Count.BeadType, sizeof *opt->bt);
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt->bt, System)) {
    InitBoolArray(opt->bt, Count->BeadType, true);
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // arrays for all the necessary stuff //{{{
  // count molecules of each type
  int *c_mtype = calloc(Count->MoleculeType, sizeof *c_mtype);
  // // count per moltype intramolecular contacts //{{{
  // ArrNDli *intra_mol = CreateArr3Dli(Count->MoleculeType,
  //                                    Count->BeadType, Count->BeadType);
  // size_t shape_intra_3body[4] = {Count->MoleculeType, Count->BeadType,
  //                                Count->BeadType, Count->BeadType};
  // ArrNDli *intra_3body = CreateArrNDli(4, shape_intra_3body);
  // count molecules of each type
  // ArrNDi *c_mtype_mtype = CreateArr2Di(Count->MoleculeType,
  //                                      Count->MoleculeType);
  // count per moltype-moltype pair intermolecular contacts
  // size_t shape_inter_mol[4] = {Count->MoleculeType,
  //                              Count->MoleculeType,
  //                              Count->BeadType,
  //                              Count->BeadType};
  // ArrNDli *inter_mol = CreateArrNDli(4, shape_inter_mol);
  // size_t shape_inter_3body[5] = {Count->MoleculeType,
  //                                Count->MoleculeType,
  //                                Count->BeadType,
  //                                Count->BeadType,
  //                                Count->BeadType};
  // ArrNDli *inter_3body = CreateArrNDli(5, shape_inter_3body); //}}}
  //}}}

  // print initial stuff to output file //{{{
  FILE *fw = PrintBylineOpenFile(fout, argc, argv);
  column = 1;
  fprintf(fw, "# (%d) step\n", column++);
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      fprintf(fw, "# molecule %s: ", System.MoleculeType[i].Name);
      PrintHeader_ctx print_ctx = { column };
      iterate_btypes(i, fw, System, opt, PrintHeader, &print_ctx);
      putc('\n', fw);
    }
  }
  fclose(fw); //}}}

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
      // printf("\nmax_contacts: %d\n", max_contacts);
      size_t shape[4] = {Count->MoleculeType,
                         Count->BeadType,
                         Count->BeadType,
                         2}; // ion: 0 - name_mol (CA2), 1 - name (C)
      ArrNDi *count_3body_step = CreateArrNDi(4, shape);
      // is the mono-/divalent counterion already in a trio?
      bool *used_name = NULL;
      if (bt_name != -1) {
        used_name = calloc(BType_name->Number, sizeof *used_name);
      }
      bool *used_name_mol = NULL;
      if (mt_name != -1) {
        used_name_mol = calloc(MolType_name->Number, sizeof *used_name_mol);
      }
      // array of ids of selected beads in selected molecules
      int *mt_beads = calloc(Count->BeadCoor, sizeof *mt_beads);
      int count_mt_beads = 0;
      for (int i = 0; i < Count->BondedCoor; i++) {
        int id_i = System.Bonded[i];
        BEAD *b_i = &System.Bead[id_i];
        MOLECULE *m_i = &System.Molecule[b_i->Molecule];
        if (opt->mt[m_i->Type] && opt->bt[b_i->Type]) {
          mt_beads[count_mt_beads] = id_i;
          count_mt_beads++;
        }
      }

      // TODO: some vmd tcl print
      char tcl[LINE] = "";
      snprintf(tcl, LINE, "contacts-%04d.tcl", count_coor - 1);
      FILE *out_vmd = OpenFile(tcl, "w");
      for (int i = 0; i < count_mt_beads; i++) {
        int id_i = mt_beads[i];
        BEAD *b_i = &System.Bead[id_i];
        MOLECULE *m_i = &System.Molecule[b_i->Molecule];
        // 1) 'name_mol' molecules
        // ... mt_name = id of CA2 name_mol
        // ... MolType_name = MOLECULETYPE thingy of the mt_name
        for (int j = 0; mt_name != -1 && j < MolType_name->Number; j++) {
          MOLECULE *mol = &System.Molecule[MolType_name->Index[j]];
          if ((opt->multi && used_name_mol[j]) || !mol->InTimestep) {
            continue;
          }
          for (int k = 0; k < MolType_name->nBeads; k++) {
            BEAD *b_k = &System.Bead[mol->Bead[k]];
            double d = DistLength(b_i->Position.v, b_k->Position.v,
                                  boxlength);
            if (d > dist_check) {
              continue;
            }
            for (int l = (i + 1); l < count_mt_beads; l++) {
              int id_l = mt_beads[l];
              BEAD *b_l = &System.Bead[id_l];
              MOLECULE *m_l = &System.Molecule[b_l->Molecule];
              double dist = DistLength(b_l->Position.v, b_k->Position.v,
                                       boxlength);
              if (dist > dist_check) {
                continue;
              }
              // for beads in one polymer, ignore beads too close to each other
              if (b_i->Molecule == b_l->Molecule &&
                  abs(PosInMol(id_i, System) -
                      PosInMol(id_l, System)) <= skip) {
                continue;
              }
              used_name_mol[j] = true;
              Types mt = SortTypes(m_i->Type, m_l->Type);
              Types bt = SortTypes(b_i->Type, b_l->Type);
              size_t id4[4] = {mt.a, bt.a, bt.b, 0};
              AddArrND(count_3body_step, id4, 1);
              // TODO: some vmd tcl print
              fprintf(out_vmd, "set rep [expr $rep + 1]\n");
              fprintf(out_vmd, "mol addrep ${mol}\n");
              fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
              fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d or resid %d\n",
                      id_i, id_l, mol->Index);
              if (!opt->multi && used_name_mol[j]) {
                break;
              }
            }
            if (!opt->multi && used_name_mol[j]) {
              break;
            }
          }
        }
        // 2) 'name' beads
        // ... bt_name = id of C bead type name
        // ... BType_name = BEADTYPE thingy of the bt_name
        for (int j = 0; bt_name != -1 && j < BType_name->InCoor; j++) {
          int id_j = BType_name->Index[j];
          if (used_name[j]) { // bead j already in a trio
            continue;
          }
          BEAD *b_j = &System.Bead[id_j];
          double d = DistLength(b_i->Position.v, b_j->Position.v, boxlength);
          if (d > dist_check) {
            continue;
          }
          for (int l = (i + 1); l < count_mt_beads; l++) {
            int id_l = mt_beads[l];
            BEAD *b_l = &System.Bead[id_l];
            MOLECULE *m_l = &System.Molecule[b_l->Molecule];
            double dist = DistLength(b_l->Position.v, b_j->Position.v,
                                     boxlength);
            if (dist > dist_check) {
              continue;
            }
            if (b_i->Molecule == b_l->Molecule &&
                abs(PosInMol(id_i, System) -
                    PosInMol(id_l, System)) <= skip) {
              continue;
            }
            used_name[j] = true;
            Types mt = SortTypes(m_i->Type, m_l->Type);
            Types bt = SortTypes(b_i->Type, b_l->Type);
            size_t id4[4] = {mt.a, bt.a, bt.b, 1};
            AddArrND(count_3body_step, id4, 1);
            // TODO: some vmd tcl print
            fprintf(out_vmd, "set rep [expr $rep + 1]\n");
            fprintf(out_vmd, "mol addrep ${mol}\n");
            fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
            fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                    id_i, id_l, id_j);
            if (!opt->multi && used_name[j]) {
              break;
            }
          }
          if (!opt->multi && used_name[j]) {
            break;
          }
        }
      }
      // TODO: some vmd tcl print
      // check wheter the vmd file is empty; i.e., no contact trios in this step
      fseek(out_vmd, 0, SEEK_END); // Move to the end of the file
      long fileSize = ftell(out_vmd); // Get the current position (file size)
      // close the vmd file
      fclose(out_vmd);
      // remove the vmd file if it's empty
      if (fileSize == 0) {
        remove(tcl);
      }
      // write average number of contacts to a file //{{{
      fw = OpenFile(fout, "a");
      fprintf(fw, "%5d", count_used);
      for (int i = 0; i < Count->MoleculeType; i++) {
        if (opt->mt[i]) {
          PrintAvgContacts_ctx avg_ctx = { count_3body_step,
                                           mt_name, // mt id for name_mol (CA2)
                                           bt_name, // bt id for name (C)
                                           i };
          iterate_btypes(i, fw, System, opt, PrintAvgContacts, &avg_ctx);
        }
      }
      putc('\n', fw);
      fclose(fw); //}}}
      // free temp arrays //{{{
      free(mt_beads);
      if (bt_name != -1) {
        free(used_name);
      }
      if (mt_name != -1) {
        free(used_name_mol);
      }
      FreeArrND(count_3body_step); //}}}
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

  // free memory - to make valgrind happy //{{{
  // FreeArrND(intra_mol);
  // FreeArrND(intra_3body);
  free(c_mtype);
  // FreeArrND(inter_mol);
  // FreeArrND(inter_3body);
  // FreeArrND(c_mtype_mtype);
  free(opt->mt);
  free(opt->bt);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}

// ...just a backup of the main many loop //{{{
      // char tcl[LINE] = "";
      // snprintf(tcl, LINE, "contacts-%04d.tcl", count_coor - 1);
      // FILE *out_vmd = OpenFile(tcl, "w");
      // for (int i = 0; i < Count->BondedCoor; i++) {
      //   int id_i = System.Bonded[i];
      // // for (int i = 0; i < count; i++) {
      // //   int id_i = beads[i];
      //   BEAD *b_i = &System.Bead[id_i];
      //   MOLECULE *m_i = &System.Molecule[b_i->Molecule];
      //   if (!opt->mt[m_i->Type] || !opt->bt[b_i->Type]) {
      //     continue;
      //   }
      //   // 1) 'name_mol' molecules
      //   // ... mt_name = id of CA2 name_mol
      //   // ... MolType_name = MOLECULETYPE thingy of the mt_name
      //   for (int j = 0; mt_name != -1 && j < MolType_name->Number; j++) {
      //     MOLECULE *mol = &System.Molecule[MolType_name->Index[j]];
      //     if (used_name_mol[j] || !mol->InTimestep) {
      //       continue;
      //     }
      //     for (int k = 0; k < MolType_name->nBeads; k++) {
      //       BEAD *b_k = &System.Bead[mol->Bead[k]];
      //       double d = DistLength(b_i->Position.v, b_k->Position.v,
      //                             boxlength);
      //       if (d > dist_check) {
      //         continue;
      //       }
      //       // for (int l = (i + 1); l < Count->BondedCoor; l++) {
      //       //   int id_l = System.Bonded[l];
      //       for (int l = (i + 1); l < count; l++) {
      //         int id_l = beads[l];
      //         BEAD *b_l = &System.Bead[id_l];
      //         MOLECULE *m_l = &System.Molecule[b_l->Molecule];
      //         if (!opt->mt[m_l->Type] || !opt->bt[b_l->Type]) {
      //           continue;
      //         }
      //         double dist = DistLength(b_l->Position.v, b_k->Position.v,
      //                                  boxlength);
      //         if (dist > dist_check) {
      //           continue;
      //         }
      //         // // a) same molecule
      //         // if (b_i->Molecule == b_l->Molecule) {
      //         //   size_t id4[4] = {mt.a, bt.a, bt.b, bt_name_mol};
      //         //   AddArrND(intra_3body, id4, 1);
      //         // // b) different molecule
      //         // } else {
      //         //   size_t id5[5] = {mt.a, mt.b, bt2.a, bt2.b, bt_name_mol};
      //         //   AddArrND(inter_3body, id5, 1);
      //         // }
      //
      //         // for beads in one polymer, ignore beads too close to each other
      //         if (b_i->Molecule == b_l->Molecule &&
      //             abs(PosInMol(id_i, System) -
      //                 PosInMol(id_l, System)) <= skip) {
      //           continue;
      //         }
      //         used_name_mol[j] = true;
      //         Types mt = SortTypes(m_i->Type, m_l->Type);
      //         Types bt = SortTypes(b_i->Type, b_l->Type);
      //         size_t id4[4] = {mt.a, bt.a, bt.b, 0};
      //         AddArrND(count_3body_step, id4, 1);
      //         // TODO: some vmd tcl print
      //         fprintf(out_vmd, "set rep [expr $rep + 1]\n");
      //         fprintf(out_vmd, "mol addrep ${mol}\n");
      //         fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
      //         // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
      //         fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d or resid %d\n",
      //                 id_i, id_l, mol->Index);
      //         if (!opt->multi && used_name_mol[j]) {
      //           // printf("%s (%d) %s (%d) %s (%d)\n",
      //           //        System.BeadType[b_i->Type].Name, id_i,
      //           //        System.BeadType[b_l->Type].Name, id_l,
      //           //        System.BeadType[bt_name].Name, mol->Index);
      //           break;
      //         }
      //       }
      //       if (!opt->multi && used_name_mol[j]) {
      //         break;
      //       }
      //     }
      //   }
      //   // 2) 'name' beads
      //   // ... bt_name = id of C bead type name
      //   // ... BType_name = BEADTYPE thingy of the bt_name
      //   for (int j = 0; bt_name != -1 && j < BType_name->InCoor; j++) {
      //     int id_j = BType_name->Index[j];
      //     if (used_name[j]) { // bead j already in a trio
      //       continue;
      //     }
      //     BEAD *b_j = &System.Bead[id_j];
      //     double d = DistLength(b_i->Position.v, b_j->Position.v, boxlength);
      //     if (d > dist_check) {
      //       continue;
      //     }
      //     // for (int l = (i + 1); l < Count->Bonded; l++) {
      //     //   int id_l = System.Bonded[l];
      //     for (int l = (i + 1); l < count; l++) {
      //       int id_l = beads[l];
      //       BEAD *b_l = &System.Bead[id_l];
      //       MOLECULE *m_l = &System.Molecule[b_l->Molecule];
      //       if (!opt->mt[m_l->Type] || !opt->bt[b_l->Type]) {
      //         continue;
      //       }
      //       double dist = DistLength(b_l->Position.v, b_j->Position.v,
      //                                boxlength);
      //       if (dist > dist_check) {
      //         continue;
      //       }
      //       // // a) same molecule
      //       // if (b_i->Molecule == b_l->Molecule) {
      //       //   // are they far enough in terms of bonds?
      //       //   if (abs(PosInMol(id_i, System) -
      //       //           PosInMol(id_l, System)) <= skip) {
      //       //     continue;
      //       //   }
      //       //   // size_t id4[4] = {mt.a, bt.a, bt.b, bt_name};
      //       //   // AddArrND(intra_3body, id4, 1);
      //       // // b) different molecule
      //       // } else {
      //       //   // size_t id5[5] = {mt.a, mt.b, bt.a, bt.b,
      //       //   //                  bt_name};
      //       //   // AddArrND(inter_3body, id5, 1);
      //       // }
      //       if (b_i->Molecule == b_l->Molecule &&
      //           abs(PosInMol(id_i, System) -
      //               PosInMol(id_l, System)) <= skip) {
      //         continue;
      //       }
      //       used_name[j] = true;
      //       Types mt = SortTypes(m_i->Type, m_l->Type);
      //       Types bt = SortTypes(b_i->Type, b_l->Type);
      //       size_t id4[4] = {mt.a, bt.a, bt.b, 1};
      //       AddArrND(count_3body_step, id4, 1);
      //       // TODO: some vmd tcl print
      //       fprintf(out_vmd, "set rep [expr $rep + 1]\n");
      //       fprintf(out_vmd, "mol addrep ${mol}\n");
      //       fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
      //       // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
      //       fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
      //               id_i, id_l, id_j);
      //       if (!opt->multi && used_name[j]) {
      //         // printf("%s (%d) %s (%d) %s (%d)\n",
      //         //        System.BeadType[b_i->Type].Name, id_i,
      //         //        System.BeadType[b_l->Type].Name, id_l,
      //         //        System.BeadType[b_j->Type].Name, id_j);
      //         break;
      //       }
      //     }
      //     if (!opt->multi && used_name[j]) {
      //       break;
      //     }
      //   }
      // } //}}}
