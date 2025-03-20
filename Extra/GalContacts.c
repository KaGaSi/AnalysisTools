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
  double dist[3];
  Distance(v1, v2, box, dist);
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
  int ****intra_step, // counts of contacts in the one step
      mt_name, // molecule type id of the divalent counterion 'name_mol'
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
      double avg = (double)(ctx->intra_step[ctx->mt][bt_j][bt_k][0]) / num_mol;
      fprintf(fw, " %lf", avg);
    }
    if (ctx->bt_name != -1) {
      double avg = (double)(ctx->intra_step[ctx->mt][bt_j][bt_k][1]) / num_mol;
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
    Help(argv[0], true, common, option);
    exit(1);
  }
  double dist_check = 0;
  if (!IsPosRealNumber(argv[++column], &dist_check)) {
    ErrorNaN("<dist>");
    Help(argv[0], true, common, option);
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
  int bt_name_mol = -1;
  MOLECULETYPE *MolType_name = NULL;
  BEADTYPE *BType_name = NULL;
  if (mt_name != -1) {
    bt_name_mol = System.MoleculeType[mt_name].BType[0];
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
  // count per moltype intramolecular contacts
  long int ***intra_mol = calloc(Count->MoleculeType, sizeof *intra_mol);
  long int ****intra_3body = calloc(Count->MoleculeType, sizeof *intra_3body);
  for (int i = 0; i < Count->MoleculeType; i++) {
    intra_mol[i] = calloc(Count->BeadType, sizeof *intra_mol[i]);
    intra_3body[i] = calloc(Count->BeadType, sizeof *intra_3body[i]);
    for (int j = 0; j < Count->BeadType; j++) {
      intra_mol[i][j] = calloc(Count->BeadType, sizeof *intra_mol[i][j]);
      intra_3body[i][j] = calloc(Count->BeadType, sizeof *intra_3body[i][j]);
      for (int k = 0; k < Count->BeadType; k++) {
        intra_3body[i][j][k] = calloc(Count->BeadType,
                                      sizeof *intra_3body[i][j][k]);
      }
    }
  }
  // count molecules of each type
  int **c_mtype_mtype = calloc(Count->MoleculeType, sizeof *c_mtype_mtype);
  // count per moltype-moltype pair intermolecular contacts
  long int ****inter_mol = calloc(Count->MoleculeType, sizeof *inter_mol);
  long int *****inter_3body = calloc(Count->MoleculeType, sizeof *inter_3body);
  for (int i = 0; i < Count->MoleculeType; i++) {
    c_mtype_mtype[i] = calloc(Count->MoleculeType, sizeof *c_mtype_mtype[i]);
    inter_mol[i] = calloc(Count->MoleculeType, sizeof *inter_mol[i]);
    inter_3body[i] = calloc(Count->MoleculeType, sizeof *inter_3body[i]);
    for (int j = 0; j < Count->MoleculeType; j++) {
      inter_mol[i][j] = calloc(Count->BeadType, sizeof *inter_mol[i][j]);
      inter_3body[i][j] = calloc(Count->BeadType, sizeof *inter_3body[i][j]);
      for (int k = 0; k < Count->BeadType; k++) {
        inter_mol[i][j][k] = calloc(Count->BeadType,
                                    sizeof *inter_mol[i][j][k]);
        inter_3body[i][j][k] = calloc(Count->BeadType,
                                      sizeof *inter_3body[i][j][k]);
        for (int l = 0; l < Count->BeadType; l++) {
          inter_3body[i][j][k][l] = calloc(Count->BeadType,
                                           sizeof *inter_3body[i][j][k][l]);
        }
      }
    }
  } //}}}

  // print initial stuff to output file //{{{
  PrintByline(fout, argc, argv);
  FILE *fw = OpenFile(fout, "a");
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
      int ****intra_step = calloc(Count->MoleculeType, sizeof *intra_step);
      for (int i = 0; i < Count->MoleculeType; i++) {
        intra_step[i] = calloc(Count->BeadType, sizeof *intra_step[i]);
        for (int j = 0; j < Count->BeadType; j++) {
          intra_step[i][j] = calloc(Count->BeadType, sizeof *intra_step[i][j]);
          for (int k = 0; k < Count->BeadType; k++) {
            intra_step[i][j][k] = calloc(2, sizeof *intra_step[i][j][k]);
          }
        }
      }
      // is the mono-/divalent counterion already in a trie?
      bool *used_name = NULL;
      if (bt_name != -1) {
        used_name = calloc(BType_name->Number, sizeof *used_name);
      }
      bool *used_name_mol = NULL;
      if (mt_name != -1) {
        used_name_mol = calloc(MolType_name->Number, sizeof *used_name_mol);
      }
      char tcl[LINE] = "";
      snprintf(tcl, LINE, "contacts-%04d.tcl", count_used);
      FILE *out_vmd = OpenFile(tcl, "w");
      for (int i = 0; i < Count->Bonded; i++) {
        int id_i = System.Bonded[i];
        BEAD *b_i = &System.Bead[id_i];
        MOLECULE *m_i = &System.Molecule[b_i->Molecule];
        if (opt->mt[m_i->Type] && opt->bt[b_i->Type]) {
          // 1) 'name_mol' molecules
          if (mt_name != -1) {
            for (int j = 0; j < MolType_name->Number; j++) {
              MOLECULE *mol = &System.Molecule[MolType_name->Index[j]];
              if (!used_name_mol[j] && mol->InTimestep) {
                for (int k = 0; k < MolType_name->nBeads; k++) {
                  BEAD *b_k = &System.Bead[mol->Bead[k]];
                  double d = DistLength(b_i->Position, b_k->Position,
                                        boxlength);
                  if (d < dist_check) {
                    for (int l = (i + 1); l < Count->Bonded; l++) {
                      int id_l = System.Bonded[l];
                      BEAD *b_l = &System.Bead[id_l];
                      MOLECULE *m_l = &System.Molecule[b_l->Molecule];
                      if (opt->mt[m_l->Type] && opt->bt[b_l->Type]) {
                        double dist = DistLength(b_l->Position, b_k->Position,
                                                 boxlength);
                        if (dist < dist_check) {
                          Types mt = SortTypes(m_i->Type, m_l->Type);
                          Types bt = SortTypes(b_i->Type, b_l->Type);
                          // a) same molecule
                          if (b_i->Molecule == b_l->Molecule) {
                            // are they far enough in terms of bonds?
                            if (abs(PosInMol(id_i, System) -
                                    PosInMol(id_l, System)) > skip) {
                              used_name_mol[j] = true;
                              intra_3body[mt.a][bt.a][bt.b][bt_name_mol]++;
                              intra_step[mt.a][bt.a][bt.b][0]++;
                              fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                              fprintf(out_vmd, "mol addrep ${mol}\n");
                              fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                              // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                              fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d or resid %d\n",
                                      id_i, id_l, mol->Index);
                            }
                          // b) different molecule
                          } else {
                            used_name_mol[j] = true;
                            Types bt = SortTypes(b_i->Type, b_l->Type);
                            inter_3body[mt.a][mt.b][bt.a][bt.b][bt_name_mol]++;
                            intra_step[mt.a][bt.a][bt.b][0]++;
                            fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                            fprintf(out_vmd, "mol addrep ${mol}\n");
                            fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                            // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                            fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d or resid %d\n",
                                    id_i, id_l, mol->Index);
                          }
                        }
                      }
                      if (!opt->multi && used_name_mol[j]) {
                        // printf("%s (%d) %s (%d) %s (%d)\n",
                        //        System.BeadType[b_i->Type].Name, id_i,
                        //        System.BeadType[b_l->Type].Name, id_l,
                        //        System.BeadType[bt_name].Name, mol->Index);
                        break;
                      }
                    }
                  }
                  if (!opt->multi && used_name_mol[j]) {
                    break;
                  }
                }
              }
            }
          }
          // 2) 'name' beads
          if (bt_name != -1) {
            for (int j = 0; j < BType_name->InCoor; j++) {
              int id_j = BType_name->Index[j];
              if (!used_name[j]) {
                BEAD *b_j = &System.Bead[id_j];
                double d = DistLength(b_i->Position, b_j->Position, boxlength);
                if (d < dist_check) {
                  for (int l = (i + 1); l < Count->Bonded; l++) {
                    int id_l = System.Bonded[l];
                    BEAD *b_l = &System.Bead[id_l];
                    MOLECULE *m_l = &System.Molecule[b_l->Molecule];
                    if (opt->mt[m_l->Type] && opt->bt[b_l->Type]) {
                      double dist = DistLength(b_l->Position, b_j->Position,
                                               boxlength);
                      if (dist < dist_check) {
                        Types mt = SortTypes(m_i->Type, m_l->Type);
                        Types bt = SortTypes(b_i->Type, b_l->Type);
                        // a) same molecule
                        if (b_i->Molecule == b_l->Molecule) {
                          // are they far enough in terms of bonds?
                          if (abs(PosInMol(id_i, System) -
                                  PosInMol(id_l, System)) > skip) {
                            used_name[j] = true;
                            intra_3body[mt.a][bt.a][bt.b][bt_name]++;
                            intra_step[mt.a][bt.a][bt.b][1]++;
                            fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                            fprintf(out_vmd, "mol addrep ${mol}\n");
                            fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                            // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                            fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                                    id_i, id_l, id_j);
                          }
                        // b) different molecule
                        } else {
                          used_name[j] = true;
                          inter_3body[mt.a][mt.b][bt.a][bt.b][bt_name]++;
                          intra_step[mt.a][bt.a][bt.b][1]++;
                          fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                          fprintf(out_vmd, "mol addrep ${mol}\n");
                          fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                          // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                          fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                                  id_i, id_l, id_j);
                        }
                      }
                    }
                    if (!opt->multi && used_name[j]) {
                      // printf("%s (%d) %s (%d) %s (%d)\n",
                      //        System.BeadType[b_i->Type].Name, id_i,
                      //        System.BeadType[b_l->Type].Name, id_l,
                      //        System.BeadType[b_j->Type].Name, id_j);
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
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
          PrintAvgContacts_ctx avg_ctx = { intra_step, mt_name, bt_name, i };
          iterate_btypes(i, fw, System, opt, PrintAvgContacts, &avg_ctx);
        }
      }
      putc('\n', fw);
      fclose(fw); //}}}
      // free temp arrays //{{{
      if (bt_name != -1) {
        free(used_name);
      }
      if (mt_name != -1) {
        free(used_name_mol);
      }
      for (int i = 0; i < Count->MoleculeType; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          for (int k = 0; k < Count->BeadType; k++) {
            free(intra_step[i][j][k]);
          }
          free(intra_step[i][j]);
        }
        free(intra_step[i]);
      }
      free(intra_step); //}}}
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
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->BeadType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        free(intra_3body[i][j][k]);
      }
      free(intra_mol[i][j]);
      free(intra_3body[i][j]);
    }
    free(intra_mol[i]);
    free(intra_3body[i]);
  }
  free(intra_mol);
  free(intra_3body);
  free(c_mtype);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < Count->MoleculeType; j++) {
      for (int k = 0; k < Count->BeadType; k++) {
        for (int l = 0; l < Count->BeadType; l++) {
          free(inter_3body[i][j][k][l]);
        }
        free(inter_mol[i][j][k]);
        free(inter_3body[i][j][k]);
      }
      free(inter_mol[i][j]);
      free(inter_3body[i][j]);
    }
    free(inter_mol[i]);
    free(inter_3body[i]);
    free(c_mtype_mtype[i]);
  }
  free(inter_mol);
  free(inter_3body);
  free(c_mtype_mtype);
  free(opt->mt);
  free(opt->bt);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
