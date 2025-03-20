#include "Aggregates.h"
#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Aggregates-NotSameBeads utility works in the same way as Aggregates utility, \
but it does not use contacts between beads of the same type; i.e., if bead \
types 'A' and 'B' are given, it considers only 'A-B' pairs.\n\n");
  }

  fprintf(ptr, "Usage: %s <input> <out.agg> <bead(s)> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<out.agg>           output aggregate file\n");
  fprintf(ptr, "<bead(s)>           bead names for closeness calculation "
          "(at least two are required)\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -d                maximum distance for contact "
          "(default: 1)\n");
  fprintf(ptr, "  -c                minimum number of contacts (default: 1, "
          "max: 255)");
  fprintf(ptr, "  -j <output>       output file with joined coordinates\n");
  fprintf(ptr, "  -w <a> <float(s)> position of wall perpendicular to "
          "given axis <a> at the axis' coordinate(s)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  double distance; // -d
  int contacts;    // -c
  FILE_TYPE fout;  // -j
  double wall[100];     // -w
  int w_count, axis;    // -w
  char w_file[2][LINE]; // -w
  FILE_TYPE j_file[2];  // -w (if -j)
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// CalculateAggregates() //{{{
// note the function doesn't fill in Aggregate[].Bead[] as it's not used here
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System, OPT opt) {
  double sqdist = Square(opt.distance);
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize - just pro-forma; is done in (Re)InitAggregate()
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
  }

  // allocate 2D triangular contact array
  uint8_t **contact = malloc(sizeof *contact * Count->Molecule);
  for (int i = 0; i < Count->Molecule; i++) {
    contact[i] = calloc(i + 1, sizeof *contact[i]);
  }

  // assign in no aggregate to each molecule
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].Aggregate = -1;
  }

  // count contacts between all molecule pairs (using cell linked list) //{{{
  // create cell-linked list
  double cell_size = sqrt(sqdist);
  int n_cells[3], *Head, *Link, Dc[14][3];
  LinkedList(*System, &Head, &Link, cell_size, n_cells, Dc);
  // go over all cells (and beads inside)
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
              if (mol_i != -1 && mol_j != -1) { // i and j must be in molecule
                if (System->BeadType[b_i->Type].Flag &&
                    System->BeadType[b_j->Type].Flag &&
                    b_i->Type != b_j->Type) { // the NotSameBeads bit
                  // calculate distance between i and j beads
                  double rij[3];
                  Distance(b_i->Position, b_j->Position,
                           System->Box.Length, rij);
                  rij[0] = Square(rij[0]) + Square(rij[1]) + Square(rij[2]);
                  // are 'i' and 'j' close enough?
                  if (mol_i != mol_j && rij[0] <= sqdist) {
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
  }
  free(Head);
  free(Link); //}}}

  EvaluateContacts(Aggregate, System, opt.contacts, contact);

  // free contact array //{{{
  for (int i = 0; i < System->Count.Molecule; i++) {
    free(contact[i]);
  }
  free(contact); //}}}

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    SortArray(Aggregate[i].Molecule, Aggregate[i].nMolecules, 0, 'i');
  } //}}}

  SortAggStruct(Aggregate, *System);
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 4, req_arg = 4, count = 0;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, false, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "-d", "-c", "-j", "-w");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // <input> - input coordinate file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output.agg> - filename of output agg file (must end with .agg) //{{{
  char agg_file[LINE] = "";
  s_strcpy(agg_file, argv[++count], LINE);
  // test if <output.agg> ends with '.agg'
  int ext = 1;
  char extension[1][EXTENSION];
  s_strcpy(extension[0], ".agg", EXTENSION);
  if (ErrorExtension(agg_file, ext, extension) == -1) {
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // -j option - save coordinates of joined aggregates
  opt->fout = InitFile;
  if (FileOption(argc, argv, "-j", opt->fout.name)) {
    opt->fout.type = CoordinateFileType(opt->fout.name);
    opt->j_file[0].type = opt->fout.type;
    opt->j_file[1].type = opt->fout.type;
  }
  // parameters for aggregate check (-d and -c options)
  opt->distance = 1;
  OneNumberOption(argc, argv, "-d", &opt->distance, 'd');
  opt->contacts = 1;
  OneNumberOption(argc, argv, "-c", &opt->contacts, 'i');
  if (opt->contacts > 255 || opt->contacts <= 0) {
    err_msg("requires whole number between 0 and 255");
    PrintErrorOption("-c");
    exit(1);
  }
  // wall options //{{{
  // -w option - wall axis and wall coordinates
  opt->w_count = 0;
  char str[LINE];
  if (FileNumbersOption(argc, argv, 1, 100, "-w", opt->wall,
                        &opt->w_count, str, 'd')) {
    if (str[0] == 'x') {
      opt->axis = 0;
    } else if (str[0] == 'y') {
      opt->axis = 1;
    } else if (str[0] == 'z') {
      opt->axis = 2;
    } else {
      err_msg("<a> requires argument 'x', 'y', or 'z'");
      PrintErrorOption("-w");
      exit(1);
    }
    // copy the agg_file to a new string, ending before '.agg'
    s_strcpy(str, agg_file, LINE);
    str[strlen(str)-4] = '\0';
    // append proper endings to the new string
    if (snprintf(opt->w_file[0], LINE, "%s-w.agg", str) < 0) {
      ErrorSnprintf();
    }
    if (snprintf(opt->w_file[1], LINE, "%s-no_w.agg", str) < 0) {
      ErrorSnprintf();
    }
    if (opt->fout.name[0] != '\0') {
      char *last_dot = strrchr(opt->fout.name, '.');
      size_t len_before_dot = last_dot - opt->fout.name;
      s_strcpy(str, opt->fout.name, len_before_dot + 1);
      if (snprintf(opt->j_file[0].name, LINE, "%s-w%s", str, last_dot) < 0) {
        ErrorSnprintf();
      }
      if (snprintf(opt->j_file[1].name, LINE, "%s-no_w%s", str, last_dot) < 0) {
        ErrorSnprintf();
      }
    }
  } //}}}
 //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;
  if (Count->Molecule == 0) {
    err_msg("No molecules in the system");
    PrintErrorFile(in.coor.name, in.stru.name, "\0");
    exit(1);
  }

  // <bead(s)> - names of bead types to use for closeness calculation //{{{
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

  // print command to output .agg (and, possibly, coordinate) file
  PrintByline(agg_file, argc, argv);
  if (opt->fout.name[0] != '\0') {
    InitOutputCoorFile(opt->fout, System, argc, argv);
    if (opt->w_count > 0) {
      InitOutputCoorFile(opt->j_file[0], System, argc, argv);
      InitOutputCoorFile(opt->j_file[1], System, argc, argv);
    }
  }
  if (opt->w_count > 0) {
    PrintByline(opt->w_file[0], argc, argv);
    PrintByline(opt->w_file[1], argc, argv);
  }

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

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
    if (UseStep(opt->c, count_coor)) {
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
      // calculate & write joined coordinatest (-j option)
      if (opt->fout.name[0] != '\0') {
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(opt->distance, Aggregate, &System);
        bool *write = calloc(Count->Bead, sizeof *write);
        InitBoolArray(write, Count->Bead, true);
        WriteTimestep(opt->fout, System, count_coor, write, argc, argv);
        free(write);
      }

      for (int i = 0; i < Count->Aggregate; i++) {
        Aggregate[i].Flag = true;
      }
      WriteAggregates(count_coor, agg_file, System, Aggregate);

      // are there walls (-w option)? //{{{
      if (opt->w_count > 0) {
        // find aggregates touching a wall
        for (int i = 0; i < Count->Aggregate; i++) {
          Aggregate[i].Flag = false;
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol_id = Aggregate[i].Molecule[j];
            MOLECULE *mol = &System.Molecule[mol_id];
            for (int k = 0; k < System.MoleculeType[mol->Type].nBeads; k++) {
              BEAD *b = &System.Bead[mol->Bead[k]];
              if (System.BeadType[b->Type].Flag) {
                for (int l = 0; l < opt->w_count; l++) {
                  double dist = b->Position[opt->axis] - opt->wall[l];
                  if (fabs(dist) < opt->distance) {
                    Aggregate[i].Flag = true;
                    goto next;
                  }
                }
              }
            }
          }
          next:
          ;
        }
        // write the aggregates to *-w.agg file
        WriteAggregates(count_coor, opt->w_file[0], System, Aggregate);
        // write joined coordinates to *-w* file (-j option)?
        if (opt->fout.name[0] != '\0') {
          bool *write = calloc(Count->Bead, sizeof *write);
          for (int i = 0; i < Count->Aggregate; i++) {
            if (Aggregate[i].Flag) {
              for (int j = 0; j < Aggregate[i].nBeads; j++) {
                write[Aggregate[i].Bead[j]] = true;
              }
            }
          }
          WriteTimestep(opt->j_file[0], System, count_coor, write, argc, argv);
          free(write);
        }
        // reverse the Aggregate[].Flag to select aggregates in bulk
        for (int i = 0; i < Count->Aggregate; i++) {
          Aggregate[i].Flag = !Aggregate[i].Flag;
        }
        // write the aggregates to *-no_w.agg file
        WriteAggregates(count_coor, opt->w_file[1], System, Aggregate);
        // write joined coordinates to *-no_w* file (-j option)?
        if (opt->fout.name[0] != '\0') {
          bool *write = calloc(Count->Bead, sizeof *write);
          for (int i = 0; i < Count->Aggregate; i++) {
            if (Aggregate[i].Flag) {
              for (int j = 0; j < Aggregate[i].nBeads; j++) {
                write[Aggregate[i].Bead[j]] = true;
              }
            }
          }
          WriteTimestep(opt->j_file[1], System, count_coor, write, argc, argv);
          free(write);
        }
      } //}}}

      ReInitAggregate(System, Aggregate);
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
  fclose(fr);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // print last step number to <output.agg>
  // open output .agg file for appending
  FILE *fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "Last Step: %d\n", count_coor);
  fclose(fw_agg);
  if (opt->w_count > 0) {
    for (int i = 0; i < 2; i++) {
      fw_agg = OpenFile(opt->w_file[i], "a");
      fprintf(fw_agg, "Last Step: %d\n", count_coor);
      fclose(fw_agg);
    }
  }

  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  free(opt);

  return 0;
}
