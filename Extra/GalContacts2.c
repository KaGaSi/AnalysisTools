#include "../AnalysisTools.h"

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
  fprintf(ptr, "Usage: %s <input> <output> <double> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file\n");
  fprintf(ptr, "<skip>              number of in-between beads to skip\n");
  fprintf(ptr, "<dist>              minimum distance contact check\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -mt <name(s)>     use specified molecule type(s)\n");
  fprintf(ptr, "  -bt <name(s)>     use specified bead type(s)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  // here com option variables
  bool bool_opt;     // --bool
  bool *mt;          // -mt
  bool *bt;          // -bt
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

int main(int argc, char *argv[]) {

  // define options & check their validity //{{{
  int common = 8, all = common + 2, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent", "--help",
              "--version", "-mt", "-bt"); //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // mandatory options //{{{
  // <input> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <output> - output file name
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);
  // <skip> - how many beads to skip at least betweem contact-able beads
  long skip = 0;
  if (!IsWholeNumber(argv[++count], &skip)) {
    ErrorNaN("<skip>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  }
  double dist_check = 0;
  if (!IsPosRealNumber(argv[++count], &dist_check)) {
    ErrorNaN("<dist>");
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  } //}}}

  // options before reading system data
  opt->c = CommonOptions(argc, argv, in);

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
  if (!MoleculeTypeOption(argc, argv, "-mt", true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }
  // bead types to calculate contacts for
  opt->bt = calloc(System.Count.BeadType, sizeof *opt->bt);
  if (!BeadTypeOption(argc, argv, "-bt", true, opt->bt, System)) {
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
  }
  int *contacts1 = calloc(Count->Bonded * Count->Bonded, sizeof *contacts1);
  int *contacts2 = calloc(Count->Bonded * Count->Bonded, sizeof *contacts2);
  //}}}

  // print initial stuff to output file //{{{
  PrintByline(fout, argc, argv);
  FILE *fw = OpenFile(fout, "a");
  count = 1;
  fprintf(fw, "# (%d) step\n", count++);
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      fprintf(fw, "# molecule %s: ", System.MoleculeType[i].Name);
      for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
        for (int k = j; k < System.MoleculeType[i].nBTypes; k++) {
          int btj = System.MoleculeType[i].BType[j];
          int btk = System.MoleculeType[i].BType[k];
          if (opt->bt[btj] && opt->bt[btk]) {
            if (FindBeadType(name, System) != -1) {
              fprintf(fw, " (%d) %s-%s-%s;", count++, System.BeadType[btj].Name,
                      System.BeadType[btk].Name, name);
            }
            if (FindMoleculeName(name_mol, System) != -1) {
              fprintf(fw, " (%d) %s-%s-%s", count++, System.BeadType[btj].Name,
                      System.BeadType[btk].Name, name_mol);
            }
          }
        }
      }
      putc('\n', fw);
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
        used_name = calloc(Count->BeadCoor, sizeof *used_name);
      }
      bool *used_name_mol = NULL;
      if (mt_name != -1) {
        used_name_mol = calloc(Count->BeadCoor, sizeof *used_name_mol);
      }
      // create cell-linked list
      double cell_size = sqrt(dist_check);
      int n_cells[3], *Head, *Link, Dc[14][3];
      LinkedList(System, &Head, &Link, cell_size, n_cells, Dc);
      char tcl[LINE] = "";
      snprintf(tcl, LINE, "contacts2-%04d.tcl", count_used);
      FILE *out_vmd = OpenFile(tcl, "w");
      // find the trios //{{{
      int c1[3];
      for (c1[2] = 0; c1[2] < n_cells[2]; c1[2]++) {
        for (c1[1] = 0; c1[1] < n_cells[1]; c1[1]++) {
          for (c1[0] = 0; c1[0] < n_cells[0]; c1[0]++) {
            int cell1 = SelectCell1(c1, n_cells);
            // select first bead in the cell 'cell1'
            int i = Head[cell1];
            while (i != -1) {
              BEAD *b_i = &System.Bead[System.BeadCoor[i]];
              int mol_i = b_i->Molecule;
              if (opt->bt[b_i->Type] && mol_i != -1 &&
                  opt->mt[System.Molecule[mol_i].Type]) {
                for (int k = 0; k < 14; k++) {
                  int cell2 = SelectCell2(c1, n_cells, Dc, k);
                  // select bead in the cell 'cell2' //{{{
                  int j;
                  if (cell1 == cell2) { // next bead in 'cell1'
                    j = Link[i];
                  } else { // first bead in 'cell2'
                    j = Head[cell2];
                  } //}}}

                  while (j != -1) {
                    BEAD *b_j = &System.Bead[System.BeadCoor[j]];
                    int bt_j = FindBeadType(System.BeadType[b_j->Type].Name,
                                            System);
                    if (bt_j == bt_name) {
                      double dist = DistLength(b_i->Position, b_j->Position,
                                               boxlength);
                      if (dist < dist_check) {
                        for (int kk = 0; kk < 14; kk++) {
                          int cell3 = SelectCell2(c1, n_cells, Dc, kk);
                          // select bead in the cell 'cell3' //{{{
                          int jj;
                          if (cell2 == cell3) { // next bead in 'cell1'
                            jj = Link[j];
                          } else { // first bead in 'cell2'
                            jj = Head[cell3];
                          } //}}}

                          while (jj != -1) {
                            BEAD *b_jj = &System.Bead[System.BeadCoor[jj]];
                            int mol_jj = b_jj->Molecule;
                            if (opt->bt[b_jj->Type] && mol_jj != -1 &&
                                opt->mt[System.Molecule[mol_jj].Type]) {
                              double distt = DistLength(b_jj->Position, b_j->Position,
                                                        boxlength);
                              if (distt < dist_check) {
                                used_name[j] = true;
                                MOLECULE *m_i = &System.Molecule[b_i->Molecule];
                                MOLECULE *m_jj = &System.Molecule[b_jj->Molecule];
                                Types mt = SortTypes(m_i->Type, m_jj->Type);
                                Types bt = SortTypes(b_i->Type, b_jj->Type);
                                // a) same molecule
                                if (b_i->Molecule == b_jj->Molecule) {
                                  // are they far enough in terms of bonds?
                                  if (abs(PosInMol(System.BeadCoor[i], System) -
                                          PosInMol(System.BeadCoor[jj], System)) > skip) {
                                    intra_3body[mt.a][bt.a][bt.b][bt_name]++;
                                    intra_step[mt.a][bt.a][bt.b][1]++;
                                    fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                                    fprintf(out_vmd, "mol addrep ${mol}\n");
                                    fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                                    // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                                    fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                                            System.BeadCoor[i], System.BeadCoor[jj], System.BeadCoor[j]);
                                  }
                                // b) different molecule
                                } else {
                                  inter_3body[mt.a][mt.b][bt.a][bt.b][bt_name]++;
                                  intra_step[mt.a][bt.a][bt.b][1]++;
                                  fprintf(out_vmd, "set rep [expr $rep + 1]\n");
                                  fprintf(out_vmd, "mol addrep ${mol}\n");
                                  fprintf(out_vmd, "mol modstyle  ${rep} ${mol} cpk 1.0 0.0\n");
                                  // fprintf(out_vmd, "mol modcolor  ${rep} ${mol} ColorID 0\n");
                                  fprintf(out_vmd, "mol modselect ${rep} ${mol} index %d %d %d\n",
                                          System.BeadCoor[i], System.BeadCoor[jj], System.BeadCoor[j]);
                                }
                              }
                            }
                            jj = Link[jj];
                            if (used_name[j]) {
                              break;
                            }
                          }
                          if (used_name[j]) {
                            break;
                          }
                        }
                      }
                    }
                    j = Link[j];
                  }
                }
              }
              i = Link[i];
            }
          }
        }
      } //}}}
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
      fprintf(fw, "%5d", count_used);
      for (int i = 0; i < Count->MoleculeType; i++) {
        if (opt->mt[i]) {
          for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
            for (int k = j; k < System.MoleculeType[i].nBTypes; k++) {
              int bt_j = System.MoleculeType[i].BType[j];
              int bt_k = System.MoleculeType[i].BType[k];
              if (opt->bt[bt_j] && opt->bt[bt_k]) {
                double avg = (double)(intra_step[i][bt_j][bt_k][0]) /
                             System.MoleculeType[i].Number;
                fprintf(fw, " %lf", avg);
                avg = (double)(intra_step[i][bt_j][bt_k][1]) /
                      System.MoleculeType[i].Number;
                fprintf(fw, " %lf", avg);
              }
            }
          }
        }
      }
      putc('\n', fw); //}}}
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
  free(contacts1);
  free(contacts2);
  free(opt->mt);
  free(opt->bt);
  free(opt);
  FreeSystem(&System);
  //}}}

  return 0;
}
