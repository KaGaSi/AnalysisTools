#include "../AnalysisTools.h"
#include "Aggregates.h"
int *InFile;

// TODO: <distance> & <contacts> as options
// TODO: do -x/-xm options make sense?
// TODO: fractional - read coor; to fraction; from fractional distance while
//       calculating
// TODO: -st/-e options? Only for -j, obviously

void Help(char cmd[50], bool error) { //{{{
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
[options]\n\n", cmd);

  fprintf(ptr, "   <input>      input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <distance>   minimum contact distance for bead pairs\n");
  fprintf(ptr, "   <contacts>   minimum number of contacts for aggregate\n");
  fprintf(ptr, "   <out.agg>    output aggregate file\n");
  fprintf(ptr, "   <bead(s)>    bead names for closeness calculation\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -x <mol(s)>    exclude specified molecule(s)\n");
  fprintf(ptr, "      -xm <mol(s)>   exclude molecules close to \
specified molecule(s)\n");
  fprintf(ptr, "      -j <out.vcf>   output vcf file with joined \
coordinates\n");
  CommonHelp(error);
} //}}}

// CalculateAggregates() //{{{
/**
 * Function to determine distribution of molecules in aggregates.
 */
void CalculateAggregates(AGGREGATE **Aggregate, COUNTS *Counts, double sqdist,
                         int contacts, int *xm_mols, bool **xm_use_mol,
                         BOX Box, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE **Molecule) {

  // zeroize //{{{
  (*Counts).Aggregates = 0;
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Aggregate)[i].nMolecules = 0;
    (*Aggregate)[i].nBeads = 0;
    (*Aggregate)[i].nMonomers = 0;
  }

  for (int i = 0; i < (*Counts).BeadsCoor; i++) {
    (*Bead)[i].nAggregates = 0;
  } //}}}

  // allocate & zeroize contact[][] (triangular matrix) and moved array //{{{
  int **contact = malloc(sizeof *contact * (*Counts).Molecules);
  int *moved = malloc(sizeof *moved * (*Counts).Molecules);
  for (int i = 0; i < (*Counts).Molecules; i++) {
    contact[i] = malloc(sizeof *contact[i] * (i + 1));
  }

  // zeroize arrays
  for (int i = 0; i < (*Counts).Molecules; i++) {
    (*Molecule)[i].Aggregate = -1;
    moved[i] = 0;
    for (int j = 0; j < i; j++)
      contact[i][j] = 0;
  } //}}}

  // create cell-linked list //{{{
  double cell_size = sqrt(sqdist);
  INTVECTOR n_cells;
  int *Head, *Link;
  int Dcx[14], Dcy[14], Dcz[14];
  LinkedList(Box.Length, *Counts, *Bead, &Head, &Link,
             cell_size, &n_cells, Dcx, Dcy, Dcz); //}}}

  // disqualify '-xm'ed molecules //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    if (xm_mols[(*Molecule)[i].Type]) {
      (*xm_use_mol)[i] = false;
    }
  }

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
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}
            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}
            while (j != -1) {
              if ((*Bead)[i].Molecule != -1 && (*Bead)[j].Molecule != -1) { // both i and j must be in molecule)
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mol_i = (*Bead)[i].Molecule;
                int mol_j = (*Bead)[j].Molecule;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                // one must be used and the other not
                if ((BeadType[btype_i].Use &&
                     (*xm_use_mol)[mol_i] && xm_mols[mtype_j]) ||
                    (BeadType[btype_j].Use &&
                     (*xm_use_mol)[mol_j] && xm_mols[mtype_i])) {

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
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y; //}}}

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
              if (mol_i != -1 && mol_j != -1) { // both i and j must be in molecule
                int btype_i = (*Bead)[i].Type;
                int btype_j = (*Bead)[j].Type;
                int mtype_i = (*Molecule)[mol_i].Type;
                int mtype_j = (*Molecule)[mol_j].Type;

                if (BeadType[btype_i].Use && BeadType[btype_j].Use && // beads must be of specified type
                    MoleculeType[mtype_i].Use && MoleculeType[mtype_j].Use && // molecules can't be excluded via -x option
                    (*xm_use_mol)[mol_i] && (*xm_use_mol)[mol_j]) { // molecules can't be excluded via -xm option

                  // TODO: fractionals?
                  // calculate distance between i and j beads
                  VECTOR rij = Distance((*Bead)[i].Position,
                                        (*Bead)[j].Position, Box.Length);
                  rij.x = SQR(rij.x) + SQR(rij.y) + SQR(rij.z);

                  // are 'i' and 'j' close enough?
                  if (mol_i != mol_j && rij.x <= sqdist) {
                    // xm option
                    if (i > j) {
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
        (*Bead)[id].nAggregates = 1;
        (*Bead)[id].Aggregate[0] = i;
      }
    }
  } //}}}

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

                int agg_j = (*Bead)[j].Aggregate[0];
                int beads_j = (*Aggregate)[agg_j].nMonomers;

                // test if 'i' is already in 'j''s aggregate //{{{
                bool in_agg = false;
                for (int l = 0; l < (*Bead)[i].nAggregates; l++) {
                  if ((*Bead)[i].Aggregate[l] == agg_j) {
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
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_j].Monomer[beads_j] = i;
                    (*Aggregate)[agg_j].nMonomers++;

                    int aggs = (*Bead)[i].nAggregates;
                    (*Bead)[i].nAggregates++;
                    (*Bead)[i].Aggregate =
                      realloc((*Bead)[i].Aggregate,
                              sizeof *(*Bead)[i].Aggregate *
                              (*Bead)[i].nAggregates);
                    (*Bead)[i].Aggregate[aggs] = agg_j;
                  }
                } //}}}
              } else if ((*Bead)[j].Molecule == -1 && // monomeric 'j'
                         (*Bead)[i].Molecule != -1) { // 'i' in molecule //{{{

                int agg_i = (*Bead)[i].Aggregate[0];
                int mono_i = (*Aggregate)[agg_i].nMonomers;

                // test if 'j' is already in 'i''s aggregate //{{{
                bool in_agg = false;
                for (int l = 0; l < (*Bead)[j].nAggregates; l++) {
                  if ((*Bead)[j].Aggregate[l] == agg_i) {
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
                  if ((SQR(rij.x)+SQR(rij.y)+SQR(rij.z)) <= sqdist) {
                    (*Aggregate)[agg_i].Monomer[mono_i] = j;
                    (*Aggregate)[agg_i].nMonomers++;

                    int aggs = (*Bead)[j].nAggregates;
                    (*Bead)[j].nAggregates++;
                    (*Bead)[j].Aggregate =
                      realloc((*Bead)[j].Aggregate,
                              sizeof *(*Bead)[j].Aggregate *
                              (*Bead)[j].nAggregates);
                    (*Bead)[j].Aggregate[aggs] = agg_i;
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
  int req_args = 5; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-xm") != 0 &&
        strcmp(argv[i], "-j") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", input_vsf[LINE] = "";
  snprintf(input_coor, LINE, "%s", argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <distance> - number of starting timestep //{{{
  // Error - non-numeric argument
  if (!IsPosReal(argv[++count])) {
    ErrorNaN("<distance>");
    Help(argv[0], true);
    exit(1);
  }
  double distance = atof(argv[count]); //}}}

  // <contacts> - number of steps to skip per one used //{{{
  // Error - non-numeric argument
  if (!IsInteger(argv[++count])) {
    ErrorNaN("<contacts>");
    Help(argv[0], true);
    exit(1);
  }
  int contacts = atoi(argv[count]); //}}}

  // <output.agg> - filename of output agg file (must end with .agg) //{{{
  char output_agg[LINE] = "";
  snprintf(output_agg, LINE, "%s", argv[++count]);

  // test if <output.agg> ends with '.agg'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(output_agg, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);

  // save coordinates of joined aggregates //{{{
  char joined_vcf[LINE];
  if (JoinCoorOption(argc, argv, joined_vcf)) {
    exit(1);
  }
  // test if <out.vcf> filename ends with '.vcf'
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (joined_vcf[0] != '\0') {
    if (ErrorExtension(joined_vcf, ext, extension)) {
      Help(argv[0], true);
      exit(1);
    }
  } //}}}
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // <bead names> - names of bead types to use for closeness calculation //{{{
  // TODO: necessary to assign false?
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Use = false;
  }
  while (++count < argc && argv[count][0] != '-') {
    int type = FindBeadType(argv[count], Counts, BeadType);
    if (type == -1) {
      ErrorPrintError_old();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent bead name ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s\n", argv[count]);
      ResetColour(STDERR_FILENO);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    if (BeadType[type].Use) {
      WarnPrintWarning();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "bead name ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%s", argv[count]);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " specified more than once\n");
      ResetColour(STDERR_FILENO);
    }
    BeadType[type].Use = true;
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }

  // used molecule type = write molecule type -- for now
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  } //}}}

  // '-xm' option //{{{
  int xm_mols[Counts.TypesOfMolecules];
  if (MoleculeTypeOption2(argc, argv, "-xm", xm_mols, Counts, &MoleculeType)) {
    exit(1);
  }

  // set all individual molecules to be used - changes in CalculateAggregates()
  bool *xm_use_mol = calloc(Counts.Molecules, sizeof *xm_use_mol);
  // TODO: Huh? Why?
//// is -xm in use?
//bool xm = false;
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  if (xm_mols[i]) {
//    xm = true;
//    break;
//  }
//} //}}}

  // print command to output .agg file //{{{
  FILE *out;
  if ((out = fopen(output_agg, "w")) == NULL) {
    ErrorFileOpen(output_agg, 'w');
    exit(1);
  }
  // TODO: print byline - requires change in reading agg file
//PrintByline(out, argc, argv);
  PrintCommand(out, argc, argv);
  fclose(out); //}}}

  // open input coordinate file //{{{
  FILE *vcf;
  if ((vcf = fopen(input_coor, "r")) == NULL) {
    ErrorFileOpen(input_coor, 'r');
    exit(1);
  }
  SkipVtfStructure(vcf, struct_lines); //}}}

  // write bead type names and pbc to <joined.vcf> if '-j' option was used //{{{
  if (joined_vcf[0] != '\0') {
    // TODO: huh?
    // bead types are to be written in joined.vcf -- probably will never change,
    // since the beadtypes should correspond to those in .agg file
    for (int i = 0; i < Counts.TypesOfBeads; i++) {
      BeadType[i].Write = true;
    }
    // open <joined.vcf>
    FILE *joined;
    if ((joined = fopen(joined_vcf, "w")) == NULL) {
      ErrorFileOpen(joined_vcf, 'w');
      exit(1);
    }
    PrintByline(joined, argc, argv);
    fclose(joined);
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules,sizeof (AGGREGATE));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate (memory-heavy)
    Aggregate[i].Monomer = malloc(sizeof *Aggregate[i].Monomer *
                                  Counts.Unbonded);
    // assumes all bonded beads can be in one aggregate (memory-heavy)
    Aggregate[i].Bead = malloc(sizeof *Aggregate[i].Bead * Counts.Bonded);
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = malloc(sizeof *Aggregate[i].Molecule *
                                   Counts.Molecules);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // main loop //{{{
  count = 0; // count timesteps
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count);
    } //}}}

    for (int i = 0; i < Counts.Molecules; i++) {
      xm_use_mol[i] = true;
    }

    // read coordinates
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    // wrap box
    RestorePBC(Counts.BeadsCoor, Box, &Bead);
    // TODO: fractionals!
    CalculateAggregates(&Aggregate, &Counts, SQR(distance), contacts,
                        xm_mols, &xm_use_mol, Box, BeadType, &Bead,
                        MoleculeType, &Molecule);

    // calculate & write joined coordinatest to <out.vcf> if '-j' option is used //{{{
    if (joined_vcf[0] != '\0') {
      // TODO: fractionals!
      RemovePBCMolecules(Counts, Box, BeadType, &Bead, MoleculeType, Molecule);
      // TODO: we're in fractionals, so maybe no need for anything new?
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length,
                          BeadType, &Bead, MoleculeType, Molecule);
      // open <joined.vcf> file
      FILE *joined;
      if ((joined = fopen(joined_vcf, "a")) == NULL) {
        ErrorFileOpen(joined_vcf, 'a');
        exit(1);
      }
      FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      WriteCoorIndexed(joined, Counts, BeadType, Bead,
                       MoleculeType, Molecule, stuff, Box);
      fclose(joined);
    } //}}}

    // find the number of aggregates - remove aggregates only of excluded mols //{{{
    int no_excluded_aggs = 0;
    int test_count = 0; // to test that every molecule is in an aggregate
    for (int i = 0; i < Counts.Aggregates; i++) {
      Aggregate[i].Use = false;

      test_count += Aggregate[i].nMolecules;

      if (Aggregate[i].nMolecules != 1 || xm_use_mol[Aggregate[i].Molecule[0]]) {
        for (int j = 0; j < Aggregate[i].nMolecules; j++) {
          int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
          if (MoleculeType[moltype].Use) {
            Aggregate[i].Use = true;
            no_excluded_aggs++;
            break;
          }
        }
      }
    } //}}}

    // are all molecules accounted for? //{{{
  // TODO: change to Warning + colours
    if (test_count != Counts.Molecules) {
      ErrorPrintError_old();
      fprintf(stderr, "\033[1;31m");
      fprintf(stderr, "\nError: not all molecules were assigned to aggregates\n");
      fprintf(stderr, "       Counts.Molecules = \033[1;33m%d\033[1;31m;", Counts.Molecules);
      fprintf(stderr, " Molecules in aggregates: \033[1;33m%d\033[1;31m\n\n", test_count);
      fprintf(stderr, "\033[0m");
      exit(1);
    } //}}}

    WriteAggregates(count, output_agg, Counts, MoleculeType, Bead, Aggregate);

    // if there's no additional timestep, exit the while loop
    if (LastStep(vcf, NULL)) {
      break;
    }
  }
  fclose(vcf);
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
  if ((out = fopen(output_agg, "a")) == NULL) {
    ErrorFileOpen(output_agg, 'a');
    exit(1);
  }

  fprintf(out, "\nLast Step: %d\n", count);

  fclose(out); //}}}

  // free memory - to make valgrind happy //{{{
  free(xm_use_mol);
  FreeAggregate(Counts, &Aggregate);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff);
  //}}}

  return 0;
}
