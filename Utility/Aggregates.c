#include "Aggregates.h"
#include "../AnalysisTools.h"

// TODO: fractional

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Aggregates utility determines which molecules belong to which aggregate on \
the basis of given parameters - the maximum distance at which a pair of beads \
from different molecules is considered in contact and the minimum number of \
such contacts between two molecules to consider them as belonging to the same \
aggregate. Only distances between specified bead types are considered. \
Information about aggregates in each timestep is written to '.agg' file (see \
documentation for the format of this file), and Cartesian coordinates of \
joined aggregates can be written to an output coordinate file (to be used \
for visualization or further analysis by other utilities).\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <out.agg> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<out.agg>           output aggregate file\n");
  fprintf(ptr, "<bead(s)>           bead names for closeness calculation\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -d                maximum distance for contact "
          "(default: 1)\n");
  fprintf(ptr, "  -c                minimum number of contacts (default: 1)\n");
  fprintf(ptr, "  -j <output>       output file with joined coordinates\n");
  CommonHelp(error, n, opt);
} //}}}

// CalculateAggregates() //{{{
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System,
                         double sqdist, int contacts) {
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
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
            int cell2 = c2x + c2y * n_cells.x + c2z * n_cells.x * n_cells.y;
            //}}}

            // select bead in the cell 'cell2' //{{{
            int j;
            if (cell1 == cell2) { // next bead in 'cell1'
              j = Link[i];
            } else { // first bead in 'cell2'
              j = Head[cell2];
            } //}}}

            while (j != -1) {
              BEAD *b_i = &System->Bead[System->BeadCoor[i]],
                   *b_j = &System->Bead[System->BeadCoor[j]];
              int mol_i = b_i->Molecule,
                  mol_j = b_j->Molecule;
              if (mol_i != -1 &&
                  mol_j != -1) { // both i and j must be in molecule
                if (System->BeadType[b_i->Type].Flag &&
                    System->BeadType[b_j->Type].Flag) {
                  // TODO: fractionals?
                  // calculate distance between i and j beads
                  VECTOR rij = Distance(b_i->Position,
                                        b_j->Position, System->Box.Length);
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

  EvaluateContacts(Aggregate, System, contacts, contact);

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    SortArray(Aggregate[i].Molecule, Aggregate[i].nMolecules, 0);
  } //}}}

  // assign bonded beads to Aggregate struct //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    // go through all molecules in aggregate 'i'
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      // copy all bead in molecule 'mol' to Aggregate struct
      int mtype = System->Molecule[mol].Type;
      for (int k = 0; k < System->MoleculeType[mtype].nBeads; k++) {
        int beads = Aggregate[i].nBeads;
        Aggregate[i].Bead[beads] = System->Molecule[mol].Bead[k];
        Aggregate[i].nBeads++;
      }
    }
  } //}}}

  // assign aggregate id to every bonded bead in the aggregate //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      int mtype = System->Molecule[mol].Type;
      for (int k = 0; k < System->MoleculeType[mtype].nBeads; k++) {
        int id = System->Molecule[mol].Bead[k];
        System->Bead[id].Aggregate = i;
      }
    }
  } //}}}

  SortAggStruct(Aggregate, *System);

  // free memory //{{{
  free(Head);
  free(Link);
  for (int i = 0; i < System->Count.Molecule; i++)
    free(contact[i]);
  free(contact);
  free(moved); //}}}
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 10, all = common + 3, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-d");
  strcpy(option[count++], "-c");
  strcpy(option[count++], "-j");
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
  bool silent, verbose, detailed;
  int start = 1, end = -1, skip = 0, pbc_xyz = -1;
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed,
                &pbc_xyz, &start, &end, &skip);

  // -j option - save coordinates of joined aggregates //{{{
  char join_file[LINE];
  int join_type;
  if (JoinCoorOption(argc, argv, &join_type, join_file)) {
    exit(1);
  } //}}}

  // parameters for aggregate check (-d and -c options)
  double distance = 1;
  DoubleOption1(argc, argv, "-d", &distance);
  int contacts = 1;
  IntegerOption1(argc, argv, "-c", &contacts);

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  SYSTEM System = ReadStructure(struct_type, struct_file,
                                coor_type, coor_file, detailed, pbc_xyz);
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

  // set all beads to write to output coordinate file
  bool *write = calloc(Count->Bead, sizeof *write);
  InitBoolArray(write, Count->Bead, true);

  // print command to output .agg file //{{{
  FILE *fw_agg = OpenFile(agg_file, "w");
  PrintByline(fw_agg, argc, argv);
  fclose(fw_agg); //}}}

  if (join_file[0] != '\0') {
    FILE *fw_coor = OpenFile(join_file, "w");
    PrintByline(fw_coor, argc, argv);
    fclose(fw_coor);
  }

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = malloc(Count->Molecule * sizeof *Aggregate);
  for (int i = 0; i < Count->Molecule; i++) {
    // assume all molecules can be in one aggregate
    Aggregate[i].Molecule = malloc(Count->Molecule *
                                   sizeof *Aggregate[i].Molecule);
    // assumes all bonded beads can be in one aggregate (memory-heavy)
    Aggregate[i].Bead = malloc(Count->Bonded * sizeof *Aggregate[i].Bead);
  } //}}}

  if (verbose) {
    VerboseOutput(System);
  }

  FILE *fr = OpenFile(coor_file, "r");
  // main loop //{{{
  int count_coor = 0,
      count_used = 0,
      line_count = 0;
  while (true) {
    PrintStep(&count_coor, 1, silent);
    // decide whether this timestep is to be saved
    bool use = false;
    if (count_coor >= start && (count_coor <= end || end == -1) &&
       ((count_coor - start) % skip) == 0) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(coor_type, fr, coor_file, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      CalculateAggregates(Aggregate, &System, SQR(distance), contacts);
      // calculate & write joined coordinatest to <out.vcf> if '-j' option is used
      if (join_file[0] != '\0') {
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(distance, Aggregate, &System);
        WriteTimestep(join_type, join_file, System, count_coor, write);
      }

      for (int i = 0; i < Count->Aggregate; i++) {
        Aggregate[i].Flag = true;
      }
      WriteAggregates(count_coor, agg_file, System, Aggregate);
      //}}}
    } else { //{{{
      if (!SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == end) {
      break;
    }
  }
  // print last step count?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}
  fclose(fr);

  // print last step number to <output.agg> //{{{
  // open output .agg file for appending
  fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "Last Step: %d\n", count_coor);
  fclose(fw_agg); //}}}

  // free memory //{{{
  free(write);
  FreeAggregate(System.Count, Aggregate);
  FreeSystem(&System);
  //}}}

  return 0;
}
