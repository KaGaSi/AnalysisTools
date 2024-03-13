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

// structure for options //{{{
struct OPT {
  double distance; // -d
  int contacts;    // -c
  FILE_TYPE fout;  // -j
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// CalculateAggregates() //{{{
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System, OPT opt) {
  double sqdist = SQR(opt.distance);
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
  }

  // allocate & zeroize contact[][] (triangular matrix) and moved array //{{{
  int **contact = malloc(sizeof *contact * Count->Molecule);
  int *moved = calloc(Count->Molecule, sizeof *moved);
  for (int i = 0; i < Count->Molecule; i++) {
    contact[i] = calloc(i + 1, sizeof *contact[i]);
  }
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].Aggregate = -1;
  } //}}}

  // create cell-linked list
  double cell_size = sqrt(sqdist);
  int n_cells[3], *Head, *Link, Dc[14][3];
  LinkedList(*System, &Head, &Link, cell_size, n_cells, Dc);

  // count contacts between all molecules pairs (using cell linked list) //{{{
  int c1[3];
  for (c1[2] = 0; c1[2] < n_cells[2]; c1[2]++) {
    for (c1[1] = 0; c1[1] < n_cells[1]; c1[1]++) {
      for (c1[0] = 0; c1[0] < n_cells[0]; c1[0]++) {
        int cell1 = SelectCell1(c1, n_cells);
        // select first bead in the cell 'cell1'
        int i = Head[cell1];
        while (i != -1) {
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
                  double rij[3];
                  Distance(b_i->Position, b_j->Position,
                           System->Box.Length, rij);
                  rij[0] = SQR(rij[0]) + SQR(rij[1]) + SQR(rij[2]);
                  // are 'i' and 'j' close enough?
                  if (mol_i != mol_j && rij[0] <= sqdist) {
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

  EvaluateContacts(Aggregate, System, opt.contacts, contact);

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    SortArrayInt(Aggregate[i].Molecule, Aggregate[i].nMolecules, 0);
  } //}}}

  // assign bonded beads to Aggregate struct //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    AGGREGATE *agg = &Aggregate[i];
    // go through all molecules in aggregate 'i'
    for (int j = 0; j < agg->nMolecules; j++) {
      int mol = agg->Molecule[j];
      // copy all bead in molecule 'mol' to Aggregate struct
      int mtype = System->Molecule[mol].Type;
      for (int k = 0; k < System->MoleculeType[mtype].nBeads; k++) {
        int beads = agg->nBeads;
        agg->Bead[beads] = System->Molecule[mol].Bead[k];
        agg->nBeads++;
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
  int common = 9, all = common + 3, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-e");
  strcpy(option[count++], "-sk");
  strcpy(option[count++], "-i");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-d");
  strcpy(option[count++], "-c");
  strcpy(option[count++], "-j");
  OptionCheck(argc, argv, count, req_arg, common, all, option, false); //}}}

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input coordinate file //{{{
  SYS_FILES in = InitSysFiles;
  snprintf(in.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in)) {
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
  opt->c = CommonOptions(argc, argv, LINE, in);

  // -j option - save coordinates of joined aggregates //{{{
  opt->fout = InitFile;
  if (FileOption(argc, argv, "-j", opt->fout.name)) {
    opt->fout.type = CoordinateFileType(opt->fout.name);
  } //}}}

  // parameters for aggregate check (-d and -c options)
  opt->distance = 1;
  DoubleOption1(argc, argv, "-d", &opt->distance);
  opt->contacts = 1;
  IntegerOption1(argc, argv, "-c", &opt->contacts);

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, opt->c.detailed);
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

  // print command to output .agg (and, possibly, coordinate) file
  PrintByline(agg_file, argc, argv);
  if (opt->fout.name[0] != '\0') {
    if (opt->fout.type == VCF_FILE) {
      PrintByline(opt->fout.name, argc, argv);
    } else if (opt->fout.type == VTF_FILE) {
      WriteStructure(opt->fout, System, -1, false, argc, argv);
    } else {
      FILE *out = OpenFile(opt->fout.name, "w");
      fclose(out);
    }
  }

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = malloc(Count->Molecule * sizeof *Aggregate);
  for (int i = 0; i < Count->Molecule; i++) {
    AGGREGATE *agg = &Aggregate[i];
    // assume all molecules can be in one aggregate
    agg->Molecule = malloc(Count->Molecule * sizeof *agg->Molecule);
    // assumes all bonded beads can be in one aggregate (memory-heavy)
    agg->Bead = malloc(Count->Bonded * sizeof *agg->Bead);
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  FILE *fr = OpenFile(in.coor.name, "r");
  // main loop //{{{
  int count_coor = 0,
      count_used = 0,
      line_count = 0;
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // decide whether this timestep is to be saved
    bool use = false;
    if (count_coor >= opt->c.start &&
        (count_coor <= opt->c.end || opt->c.end == -1) &&
        ((count_coor - opt->c.start) % opt->c.skip) == 0) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      CalculateAggregates(Aggregate, &System, *opt);
      // calculate & write joined coordinatest to <out.vcf> if '-j' option is used
      if (opt->fout.name[0] != '\0') {
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(opt->distance, Aggregate, &System);
        WriteTimestep(opt->fout, System, count_coor, write);
      }

      for (int i = 0; i < Count->Aggregate; i++) {
        Aggregate[i].Flag = true;
      }
      WriteAggregates(count_coor, agg_file, System, Aggregate);
      //}}}
    } else { //{{{
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  // print last step count?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}
  fclose(fr);

  // print last step number to <output.agg>
  // open output .agg file for appending
  FILE *fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "Last Step: %d\n", count_coor);
  fclose(fw_agg);

  free(write);
  FreeAggregate(System.Count, Aggregate);
  FreeSystem(&System);
  free(opt);

  return 0;
}
