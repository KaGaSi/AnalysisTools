#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Aggregates utility determines which molecules belong to which aggregate on"
  "the basis of given parameters - the maximum distance at which a pair of"
  "beads from different molecules is considered in contact and the minimum "
  "number of such contacts between two molecules to consider them as belonging "
  "to the same aggregate. Only distances between specified bead types are "
  "considered. Information about aggregates in each timestep is written to "
  "'.agg' file (see documentation for the format of this file), and Cartesian "
  "coordinates of joined aggregates can be written to an output coordinate "
  "file (to be used for visualization or further analysis by other utilities). "
  "The utility can also differentiate between aggregates near a wall and those "
  "in bulk (-w option); when the axis perpendicular to the wall(s) and the"
  "coordinate(s) of the wall(s) along that axis are provided, any aggregate "
  "containing a <bead(s)> that is at most the contact distance from the wall "
  "is considered near the wall. Aggregates near the wall(s)/in bulk are saved "
  "into two files whose names are based on <out.agg> ('_w' and '_b' is "
  "prepended to the .agg extension)",

  "Usage: Aggregates <coor> <out.agg> <bead(s)>/--all [options]\n\n",
  .args = 2,
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
  {"<coor>", NULL, "input coordinate file", OPT_ARG},
  {"<out.agg>", NULL, "output aggregate file", OPT_ARG},
  {"<bead(s)>/--all", NULL, "bead names for closeness calculation", OPT_EXTRA},
  {"--all", NULL, "use all types (overwrites <bead(s)>)", OPT_EXTRA},
  {"-d", "<float>", "maximum distance for contact (default: 1)", OPT_EXTRA},
  {"-c", "<float>", "minimum number of contacts (default: 1, max: 255)", OPT_EXTRA},
  {"-j", "<output>", "output file with joined coordinates", OPT_EXTRA},
  {"-w", "<a> <float(s)>", "position of wall perpendicular to given axis <a> at the axis' coordinate(s)", OPT_EXTRA},
  {NULL}
}; //}}}

// Help() //{{{
void Help_old(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
} //}}}

// structure for options //{{{
struct OPT {
  double cutoff;      // -d
  int contacts;         // -c
  FILE_TYPE fout;       // -j
  double wall[100];     // -w
  int w_count, axis;    // -w
  char w_file[2][LINE]; // -w
  FILE_TYPE j_file[2];  // -w (if -j)
  bool all;             // --all
}; //}}}

// detect possible contact between two beads //{{{
void CalculateContacts(const int id_i, const int id_j, SYSTEM System,
                       const double dist, ArrNDi *contact) {
  int i = System.BeadCoor[id_i];
  int j = System.BeadCoor[id_j];
  int mol_i = System.Bead[i].Molecule;
  int mol_j = System.Bead[j].Molecule;
  vec3d *pos_i = &System.Bead[i].Position;
  vec3d *pos_j = &System.Bead[j].Position;
  vec3d rij = Distance(pos_i->v, pos_j->v, System.Box.Length);
  rij.v[0] = VectLength(rij);
  // are 'i' and 'j' close enough?
  if (System.Bead[i].Molecule != System.Bead[j].Molecule &&
      rij.v[0] <= dist) {
    if (mol_i > mol_j) {
      AddArr2D(contact, mol_i, mol_j, 1);
    } else {
      AddArr2D(contact, mol_j, mol_i, 1);
    }
  }
}
// structure for the callback function
struct contacts_args {
  double dist;
  ArrNDi *contact;
};
// adaptor for the CalculateContacts() function
static void CalculateContacts_adaptor(int id_i, int id_j,
                                      const SYSTEM System, void *ud) {
  struct contacts_args *p = (struct contacts_args*)ud;
  CalculateContacts(id_i, id_j, System, p->dist, p->contact);
} //}}}
// condition for using specified beads //{{{
static bool CheckBead(int id, SYSTEM System) {
  int i = System.BeadCoor[id];
  int btype = System.Bead[i].Type;
  if (!System.BeadType[btype].Flag ||
      System.Bead[i].Molecule == -1) {
    return false;
  } else {
    return true;
  }
}
// structure for the callback function (empty as only System is needed here)
struct check_args {
};
// adaptor for the CalculatePCF() function
static bool CheckBead_adaptor(int type, SYSTEM System, void *ud) {
  return CheckBead(type, System);
} //}}}

// CalculateAggregates() //{{{
// note the function doesn't fill in Aggregate[].Bead[] as it's not used here
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System, OPT opt) {
  double sqdist = Square(opt.cutoff);
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize - just pro-forma; is done in (Re)InitAggregate()
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
  }

  // array for number of contacts between molecules
  ArrNDi *contact = CreateArr2Di(Count->Molecule, Count->Molecule);

  // assign in no aggregate to each molecule
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].Aggregate = -1;
  }
  // calculate contact pairs
  double cell_size = sqrt(sqdist);
  struct contacts_args args = { sqrt(sqdist), contact};
  struct check_args check = { };
  TraversePairs(*System, cell_size, CalculateContacts_adaptor, &args,
                CheckBead_adaptor, &check);

  EvaluateContacts(Aggregate, System, opt.contacts, contact);

  FreeArrND(contact);

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    SortArray(Aggregate[i].Molecule, Aggregate[i].nMolecules, 0, 'i');
  } //}}}

  SortAggStruct(Aggregate, *System);
} //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
  OPT opt;
  int count = 0;
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
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // options before reading system data //{{{
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // -j option - save coordinates of joined aggregates
  opt.fout = InitFile;
  if (FileOption(argc, argv, "-j", opt.fout.name)) {
    opt.fout.type = CoordinateFileType(opt.fout.name);
    opt.j_file[0].type = opt.fout.type;
    opt.j_file[1].type = opt.fout.type;
  }
  // parameters for aggregate check (-d and -c options)
  opt.cutoff = 1;
  OneNumberOption(argc, argv, "-d", &opt.cutoff, 'd');
  opt.contacts = 1;
  OneNumberOption(argc, argv, "-c", &opt.contacts, 'i');
  if (opt.contacts > 255 || opt.contacts <= 0) {
    err_msg("requires whole number between 0 and 255");
    PrintErrorOption("-c");
    exit(1);
  }
  // wall options //{{{
  // -w option - wall axis and wall coordinates
  opt.w_count = 0;
  char str[LINE];
  if (FileNumbersOption(argc, argv, 1, 100, "-w", opt.wall,
                        &opt.w_count, str, 'd')) {
    if (str[0] == 'x') {
      opt.axis = 0;
    } else if (str[0] == 'y') {
      opt.axis = 1;
    } else if (str[0] == 'z') {
      opt.axis = 2;
    } else {
      err_msg("<a> requires argument 'x', 'y', or 'z'");
      PrintErrorOption("-w");
      exit(1);
    }
    // copy the agg_file to a new string, ending before '.agg'
    s_strcpy(str, agg_file, LINE);
    str[strlen(str)-4] = '\0';
    // append proper endings to the new string
    if (snprintf(opt.w_file[0], LINE, "%s_w.agg", str) < 0) {
      ErrorSnprintf();
    }
    if (snprintf(opt.w_file[1], LINE, "%s_b.agg", str) < 0) {
      ErrorSnprintf();
    }
    if (opt.fout.name[0] != '\0') {
      char *last_dot = strrchr(opt.fout.name, '.');
      size_t len_before_dot = last_dot - opt.fout.name;
      s_strcpy(str, opt.fout.name, len_before_dot + 1);
      if (snprintf(opt.j_file[0].name, LINE, "%s_w%s", str, last_dot) < 0) {
        ErrorSnprintf();
      }
      if (snprintf(opt.j_file[1].name, LINE, "%s_b%s", str, last_dot) < 0) {
        ErrorSnprintf();
      }
    }
  } //}}}
  opt.all = BoolOption(argc, argv, "--all"); //}}}
  //}}}

  if (!commons.silent) {
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
  if (opt.all) {
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = true;
    }
  } else {
    // missing --all as well as any bead type(s)
    // TODO: necessary to assign false?
    for (int i = 0; i < Count->BeadType; i++) {
      System.BeadType[i].Flag = false;
    }
    while (++count < argc && argv[count][0] != '-') {
      int type = FindBeadType(argv[count], System);
      if (type == -1) {
        err_msg("non-existent bead name");
        PrintError();
        ErrorBeadType(argv[count], System);
        exit(1);
      }
      if (System.BeadType[type].Flag) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
                 ErrYellow(), argv[count], ErrCyan());
        PrintWarning();
      }
      System.BeadType[type].Flag = true;
    }
    count--; // while always increments count at least once
    if (count < (HelpDesc.args + 1)) {
      err_msg("missing <bead(s)> or --all option");
      PrintError();
      PrintCommand(stderr, argc, argv);
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}

  // print command to output .agg (and, possibly, coordinate) file
  PrintByline(agg_file, argc, argv);
  if (opt.fout.name[0] != '\0') {
    InitOutputCoorFile(opt.fout, System, argc, argv);
    if (opt.w_count > 0) {
      InitOutputCoorFile(opt.j_file[0], System, argc, argv);
      InitOutputCoorFile(opt.j_file[1], System, argc, argv);
    }
  }
  if (opt.w_count > 0) {
    PrintByline(opt.w_file[0], argc, argv);
    PrintByline(opt.w_file[1], argc, argv);
  }

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  FILE *fr = OpenFile(in.coor.name, "r");
  // main loop //{{{
  int count_coor = 0,
      count_used = 0,
      line_count = 0;
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // decide whether this timestep is to be saved
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    }
    if (use) { //{{{
      if (!ReadTimestep(in, fr, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      WrapJoinCoordinates(&System, true, false);
      CalculateAggregates(Aggregate, &System, opt);
      // calculate & write joined coordinatest (-j option)
      if (opt.fout.name[0] != '\0') {
        FillAggregateBeads(Aggregate, System);
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(opt.cutoff, Aggregate, &System);
        bool *write = calloc(Count->Bead, sizeof *write);
        InitBoolArray(write, Count->Bead, true);
        WriteTimestep(opt.fout, System, count_coor, write, argc, argv);
        free(write);
      }

      for (int i = 0; i < Count->Aggregate; i++) {
        Aggregate[i].Flag = true;
      }
      WriteAggregates(count_coor, agg_file, System, Aggregate);

      // are there walls (-w option)? //{{{
      if (opt.w_count > 0) {
        // find aggregates touching a wall
        for (int i = 0; i < Count->Aggregate; i++) {
          Aggregate[i].Flag = false;
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol_id = Aggregate[i].Molecule[j];
            MOLECULE *mol = &System.Molecule[mol_id];
            for (int k = 0; k < System.MoleculeType[mol->Type].nBeads; k++) {
              BEAD *b = &System.Bead[mol->Bead[k]];
              if (System.BeadType[b->Type].Flag) {
                for (int l = 0; l < opt.w_count; l++) {
                  double dist = b->Position.v[opt.axis] - opt.wall[l];
                  if (fabs(dist) < opt.cutoff) {
                    Aggregate[i].Flag = true; // aggregate i is touching a wall
                    goto next;
                  }
                }
              }
            }
          }
          next:
          ;
        }
        // write the aggregates to *_w.agg file
        WriteAggregates(count_coor, opt.w_file[0], System, Aggregate);
        // reverse the Aggregate[].Flag to select aggregates in bulk
        for (int i = 0; i < Count->Aggregate; i++) {
          Aggregate[i].Flag = !Aggregate[i].Flag;
        }
        // write the aggregates to *_b.agg file
        WriteAggregates(count_coor, opt.w_file[1], System, Aggregate);

        // write joined coordinates to _b/_w files (-j option)?
        if (opt.fout.name[0] != '\0') {
          bool *write = calloc(Count->Bead, sizeof *write);
          // assume all beads are saved (to save unbonded beads)
          InitBoolArray(write, Count->Bead, true);

          // exclude from saving all aggregate beads in bulk
          for (int i = 0; i < Count->Aggregate; i++) {
            // is aggregate in the bulk?
            if (Aggregate[i].Flag) {
              for (int j = 0; j < Aggregate[i].nMolecules; j++) {
                int mol = Aggregate[i].Molecule[j];
                int mtype = System.Molecule[mol].Type;
                for (int k = 0; k < System.MoleculeType[mtype].nBeads; k++) {
                  int id = System.Molecule[mol].Bead[k];
                  if (System.Bead[id].InTimestep) {
                    write[id] = false;
                  }
                }
              }
            }
          }
          // write joined coordinates for wall-touching aggregates to _w file
          WriteTimestep(opt.j_file[0], System, count_coor, write, argc, argv);

          // flip write flag for all bonded beads
          for (int i = 0; i < Count->Bonded; i++) {
            int id = System.Bonded[i];
            write[id] = !write[id];
          }
          // write joined coordinates for bulk aggregates to _b file
          WriteTimestep(opt.j_file[1], System, count_coor, write, argc, argv);
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
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, commons.silent); //}}}

  // print last step number to <output.agg>
  // open output .agg file for appending
  FILE *fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "Last Step: %d\n", count_coor);
  fclose(fw_agg);
  if (opt.w_count > 0) {
    for (int i = 0; i < 2; i++) {
      fw_agg = OpenFile(opt.w_file[i], "a");
      fprintf(fw_agg, "Last Step: %d\n", count_coor);
      fclose(fw_agg);
    }
  }

  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);

  return 0;
}
