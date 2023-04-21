#include "Aggregates.h"
#include "../AnalysisTools.h"

// TODO: -sk/-e/-st options
// TODO: <distance> & <contacts> as options
// TODO: fractional - read coor; to fraction; from fractional distance while
//       calculating

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
Aggregates utility determines which molecules belong to which aggregate on \
the basis of given parameters - the minimum distance at which a pair of beads \
from different molecules is considered in contact and the minimum number of \
such contacts between two molecules to consider them as belonging to the same \
aggregate. Only distances between specified bead types are considered. \
Information about aggregates in each timestep is written to '.agg' file (see \
documentation for the format of this file). Coordinates of joined aggregates \
can be written to an output '.vcf' file (with indexed timesteps).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <distance> <contacts> <out.agg> <bead(s)> \
[options]\n\n",
          cmd);

  fprintf(ptr, "   <input>      input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <distance>   minimum contact distance for bead pairs\n");
  fprintf(ptr, "   <contacts>   minimum number of contacts for aggregate\n");
  fprintf(ptr, "   <out.agg>    output aggregate file\n");
  fprintf(ptr, "   <bead(s)>    bead names for closeness calculation\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -j <out.vcf>   output file with joined coordinates\n");
  CommonHelp(error, n, opt);
} //}}}

// CalculateAggregates() //{{{
void CalculateAggregates(AGGREGATE **Aggregate, SYSTEM *System, double sqdist,
                         int contacts) {
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i]->nMolecules = 0;
    Aggregate[i]->nBeads = 0;
  }

  // allocate & zeroize contact[][] (triangular matrix) and moved array //{{{
  int **contact = malloc(sizeof *contact * Count->Molecule);
  int *moved = malloc(sizeof *moved * Count->Molecule);
  for (int i = 0; i < Count->Molecule; i++) {
    contact[i] = malloc(sizeof *contact[i] * (i + 1));
  }

  // zeroize arrays
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].Aggregate = -1;
    moved[i] = 0;
    for (int j = 0; j < i; j++)
      contact[i][j] = 0;
  } //}}}

  // create cell-linked list
  double cell_size = sqrt(sqdist);
  INTVECTOR n_cells;
  int *Head, *Link;
  int Dcx[14], Dcy[14], Dcz[14];
  LinkedList(*System, &Head, &Link, cell_size, &n_cells, Dcx, Dcy, Dcz);

  for (int c1z = 0; c1z < n_cells.z; c1z++) { //{{{
    for (int c1y = 0; c1y < n_cells.y; c1y++) {
      for (int c1x = 0; c1x < n_cells.x; c1x++) {
        // select first cell
        int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;
        // select first bead in the cell 'cell1'
        int i = Head[cell1];
        while (i != -1) {
          for (int k = 0; k < 14; k++) {
            int c2x = c1x + Dcx[k]; //{{{
            int c2y = c1y + Dcy[k];
            int c2z = c1z + Dcz[k];

            // periodic boundary conditions for cells
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z;

            // select second cell
            int cell2 =
                c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}
            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}
            while (j != -1) {
              if ((*Bead)[i].Molecule != -1 &&
                  (*Bead)[j].Molecule !=
                      -1) { // both i and j must be in molecule)
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mol_i = (*Bead)[i].Molecule;
                int mol_j = (*Bead)[j].Molecule;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                // one must be used and the other not
                if (BeadType[btype_i].Flag && BeadType[btype_j].Flag) {

                  // TODO: fractionals?
                  VECTOR rij = Distance((*Bead)[i].Position,
                                        (*Bead)[j].Position, Box.Length);
                  rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

                  if (rij.x <= sqdist) {
                    if ((*xm_use_mol)[mol_i]) {
                      (*xm_use_mol)[mol_i] = false;
                    } else {
                      (*xm_use_mol)[mol_j] = false;
                    }
                    break;
                  }
                }
              }
              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  } //}}}

  // count contacts between all molecules pairs (using cell linked list) //{{{
  for (int c1z = 0; c1z < n_cells.z; c1z++) {
    for (int c1y = 0; c1y < n_cells.y; c1y++) {
      for (int c1x = 0; c1x < n_cells.x; c1x++) {

        // select first cell
        int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;

        // select first bead in the cell 'cell1'
        int i = Head[cell1];

        while (i != -1) {
          for (int k = 0; k < 14; k++) {
            int c2x = c1x + Dcx[k]; //{{{
            int c2y = c1y + Dcy[k];
            int c2z = c1z + Dcz[k];

            // periodic boundary conditions for cells
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z;

            // select second cell
            int cell2 =
                c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}

            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}

            while (j != -1) {
              int mol_i = (*Bead)[i].Molecule;
              int mol_j = (*Bead)[j].Molecule;
              if (mol_i != -1 &&
                  mol_j != -1) { // both i and j must be in molecule
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                // TODO
                // if (BeadType[btype_i].Use && BeadType[btype_j].Use) {
                if (true) {

                  // TODO: fractionals?
                  // calculate distance between i and j beads
                  VECTOR rij = Distance((*Bead)[i].Position,
                                        (*Bead)[j].Position, Box.Length);
                  rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

                  // are 'i' and 'j' close enough?
                  if (mol_i != mol_j && rij.x <= sqdist) {
                    // xm option
                    if (mol_i > mol_j) {
                      contact[mol_i][mol_j]++;
                    } else {
                      contact[mol_j][mol_i]++;
                    }
                  }
                }
              }
              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  } //}}}

  EvaluateContacts(Counts, Aggregate, Molecule, contacts, contact);

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    SortArray((*Aggregate)[i].Molecule, (*Aggregate)[i].nMolecules, 0);
  } //}}}

  // assign bonded beads to Aggregate struct //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    // go through all molecules in aggregate 'i'
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      // copy all bead in molecule 'mol' to Aggregate struct
      int mtype = (*Molecule)[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
        int beads = (*Aggregate)[i].nBeads;
        (*Aggregate)[i].Bead[beads] = (*Molecule)[mol].Bead[k];
        (*Aggregate)[i].nBeads++;
      }
    }
  } //}}}

  // assign aggregate id to every bonded bead in the aggregate //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      int mtype = (*Molecule)[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
        int id = (*Molecule)[mol].Bead[k];
        (*Bead)[id].Aggregate = i;
      }
    }
  } //}}}

/*
  // find monomeric beads close to aggregates (using cell linked list) //{{{
  for (int c1z = 0; c1z < n_cells.z; c1z++) {
    for (int c1y = 0; c1y < n_cells.y; c1y++) {
      for (int c1x = 0; c1x < n_cells.x; c1x++) {

        // select first cell
        int cell1 = c1x + c1y * n_cells.x + c1z * n_cells.x * n_cells.y;

        int i = Head[cell1];

        while (i != -1) {
          for (int k = 0; k < 14; k++) {
            int c2x = c1x + Dcx[k];
            int c2y = c1y + Dcy[k];
            int c2z = c1z + Dcz[k];

            // cell periodic boundary condition //{{{
            if (c2x >= n_cells.x)
              c2x -= n_cells.x;
            else if (c2x < 0)
              c2x += n_cells.x;

            if (c2y >= n_cells.y)
              c2y -= n_cells.y;
            else if (c2y < 0)
              c2y += n_cells.y;

            if (c2z >= n_cells.z)
              c2z -= n_cells.z; //}}}

            // select second cell
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y;

            int j;
            if (cell1 == cell2) {
              j = Link[i];
            } else {
              j = Head[cell2];
            }

            while (j != -1) {
              // test if the monmeric bead is near aggregate
              if ((*Bead)[i].Molecule == -1 && // monomeric 'i'
                  (*Bead)[j].Molecule != -1) { // 'j' in molecule //{{{

                int agg_j = (*Bead)[j].Aggregate;
                int beads_j = (*Aggregate)[agg_j].nMonomers;

                // test if 'i' is already in 'j''s aggregate //{{{
                bool in_agg = false;
                for (int l = 0; l < (*Bead)[i].nAggregates; l++) {
                  if ((*Bead)[i].Aggregatexxx[l] == agg_j) {
                    in_agg = true;
                    break;
                  }
                } //}}}

                if (!in_agg) {
                  // TODO: fractionals?
                  // calculate distance between i and j beads
                  VECTOR rij = Distance((*Bead)[i].Position,
                                        (*Bead)[j].Position, Box.Length);

                  // test if 'i' is near 'j''s aggregate
                  if ((SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_j].Monomer[beads_j] = i;
                    (*Aggregate)[agg_j].nMonomers++;

                    int aggs = (*Bead)[i].nAggregates;
                    (*Bead)[i].nAggregates++;
                    (*Bead)[i].Aggregatexxx =
                        realloc((*Bead)[i].Aggregatexxx,
                                sizeof *(*Bead)[i].Aggregatexxx *
                                    (*Bead)[i].nAggregates);
                    (*Bead)[i].Aggregatexxx[aggs] = agg_j;
                  }
                }                                     //}}}
              } else if ((*Bead)[j].Molecule == -1 && // monomeric 'j'
                         (*Bead)[i].Molecule != -1) { // 'i' in molecule //{{{

                int agg_i = (*Bead)[i].Aggregatexxx[0];
                int mono_i = (*Aggregate)[agg_i].nMonomers;

                // test if 'j' is already in 'i''s aggregate //{{{
                bool in_agg = false;
                for (int l = 0; l < (*Bead)[j].nAggregates; l++) {
                  if ((*Bead)[j].Aggregatexxx[l] == agg_i) {
                    in_agg = true;
                    break;
                  }
                } //}}}

                if (!in_agg) {
                  // TODO: fractionals?
                  // calculate distance between i and j beads
                  VECTOR rij = Distance((*Bead)[i].Position,
                                        (*Bead)[j].Position, Box.Length);

                  // test if 'j' is near 'i''s aggregate
                  if ((SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_i].Monomer[mono_i] = j;
                    (*Aggregate)[agg_i].nMonomers++;

                    int aggs = (*Bead)[j].nAggregates;
                    (*Bead)[j].nAggregates++;
                    (*Bead)[j].Aggregatexxx =
                        realloc((*Bead)[j].Aggregatexxx,
                                sizeof *(*Bead)[j].Aggregatexxx *
                                    (*Bead)[j].nAggregates);
                    (*Bead)[j].Aggregatexxx[aggs] = agg_i;
                  }
                }
              } //}}}

              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  } //}}}
*/

  // sort monomers in aggregates according to ascending ids //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    SortArray((*Aggregate)[i].Monomer, (*Aggregate)[i].nMonomers, 0);
  } //}}}

  // sort aggregates according to ascending ids of first molecules
  SortAggStruct(Aggregate, *Counts, *Molecule, MoleculeType, Bead, BeadType);

  // free memory //{{{
  free(Head);
  free(Link);
  for (int i = 0; i < (*Counts).Molecules; i++)
    free(contact[i]);
  free(contact);
  free(moved); //}}}
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 11, all = common + 1, count = 0,
      req_arg = 5;
  char option[all][OPT_LENGTH];
  // common options
  // strcpy(option[count++], "-st");
  // strcpy(option[count++], "-e");
  // strcpy(option[count++], "-sk");

  strcpy(option[count++], "-i");
  strcpy(option[count++], "--variable");
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "--reverse");
  strcpy(option[count++], "--join");
  strcpy(option[count++], "--wrap");
  strcpy(option[count++], "-n");
  strcpy(option[count++], "--last");
  OptionCheck(argc, argv, req_arg, common, all, option); //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // parameters for aggregate check //{{{
  double distance;
  if (!IsPosRealNumber(argv[++count], &distance)) {
    ErrorNaN("<distance>");
    Help(argv[0], true, common, option);
    exit(1);
  }
  long contacts;
  if (!IsNaturalNumber(argv[++count], &contacts)) {
    ErrorNaN("<contacts>");
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // <output.agg> - filename of output agg file (must end with .agg) //{{{
  char agg_file[LINE] = "";
  snprintf(agg_file, LINE, "%s", argv[++count]);
  // test if <output.agg> ends with '.agg'
  int ext = 1;
  char extension[1][EXTENSION];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(agg_file, ext, extension) == -1) {
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // options before reading system data
  bool silent, verbose, detailed, vtf_var;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &start, &end, &skip);

  // save coordinates of joined aggregates //{{{
  char join_file[LINE];
  int join_type;
  if (JoinCoorOption(argc, argv, &join_type, join_file)) {
    exit(1);
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  int ltrj_start_id = -1; // for lammpstrj structure file, start ids from 0 or 1
  SYSTEM System = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);

  COUNT *Count = &System.Count;
  // <bead names> - names of bead types to use for closeness calculation //{{{
  // TODO: necessary to assign false?
  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Flag = false;
  }
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], System);
    if (type == -1) {
      ErrorBeadType(argv[count], System);
      exit(1);
    }
    if (System.BeadType[type].Flag) {
      snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
               ErrYellow(), argv[count], ErrCyan());
      PrintWarning();
    }
    System.BeadType[type].Flag = true;
  } //}}}

  // print command to output .agg file //{{{
  FILE *fw_agg = OpenFile(agg_file, "w");
  // TODO: print byline - requires change in reading agg file
  // PrintByline(out, argc, argv);
  PrintCommand(fw_agg, argc, argv);
  fclose(fw_agg); //}}}

  // open input coordinate file
  FILE *fr = OpenFile(coor_file, "r");

  // write bead type names and pbc to <joined.vcf> if '-j' option was used //{{{
  bool *write = calloc(Count->BeadType, sizeof *write);
  InitBoolArray(write, Count->BeadType, true);
  if (join_file[0] != '\0') {
    // open <joined.vcf>
    FILE *fw_coor = OpenFile(join_file, "w");
    PrintByline(fw_coor, argc, argv);
    fclose(fw_coor);
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Count->Molecule, sizeof *Aggregate);
  for (int i = 0; i < Count->Molecule; i++) {
    // assumes all monomeric beads can be near one aggregate (memory-heavy)
    Aggregate[i].Monomer = malloc(sizeof *Aggregate[i].Monomer *
                                  Count->Unbonded);
    // assumes all bonded beads can be in one aggregate (memory-heavy)
    Aggregate[i].Bead = malloc(sizeof *Aggregate[i].Bead * Count->Bonded);
    // assume all molecules can be in one aggregate
    Aggregate[i].Molecule = malloc(sizeof *Aggregate[i].Molecule *
                                   Count->Molecule);
  } //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  // main loop //{{{
  int count_coor = 0,
      count_used = 0,
      line_count = 0;
  while (true) {
    PrintStep(&count_coor, 1, silent);

    if (!ReadTimestep(coor_type, fr, coor_file, &System,
                      &line_count, vtf_var)) {
      count_coor--;
      break;
    }
    count_used++;
    WrapJoinCoordinates(&System, true, false);
    CalculateAggregates(&Aggregate, &System, SQR(distance), contacts);

    // calculate & write joined coordinatest to <out.vcf> if '-j' option is used
    // //{{{
    if (join_file[0] != '\0') {
      // TODO: fractionals!
      RemovePBCMolecules2(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
      // TODO: we're in fractionals, so maybe no need for anything new?
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length, BeadType,
                          &Bead, MoleculeType, Molecule);
      FILE *joined = OpenFile(join_file, "a");
      FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      WriteCoorIndexed(joined, Counts, BeadType, Bead, MoleculeType, Molecule,
                       stuff, Box);
      fclose(joined);
    } //}}}

    // are all molecules accounted for? //{{{
    // TODO: change to Warning + colours
    if (test_count != Counts.Molecules) {
      ErrorPrintError_old();
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr,
              "\nError: not all molecules were assigned to aggregates\n");
      fprintf(stderr, "       Counts.Molecules = \033[1;33m%d\033[1;31m;",
              Counts.Molecules);
      fprintf(stderr, " Molecules in aggregates: \033[1;33m%d\033[1;31m\n\n",
              test_count);
      fprintf(stderr, "\033[0m");
      exit(1);
    } //}}}

    WriteAggregates(count, agg_file, Counts, MoleculeType, Bead, Aggregate);

    // if there's no additional timestep, exit the while loop
    if (LastStep(fr, NULL)) {
      break;
    }
  }
  fclose(fr);
  // print last step count?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count);
  } //}}}

  // print last step number to <output.agg> //{{{
  // open output .agg file for appending
  fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "\nLast Step: %d\n", count);
  fclose(fw_agg); //}}}

  // free memory - to make valgrind happy //{{{
  free(xm_use_mol);
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  free(write);
  //}}}

  return 0;
}
