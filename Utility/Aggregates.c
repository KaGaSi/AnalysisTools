#include "../src/AnalysisTools.h"
// TODO: well, the contact array - should be spare array or some such

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Aggregates utility determines which molecules belong to which aggregate on"
  "the basis of given parameters - the maximum distance at which a pair of"
  "beads from different molecules is considered in contact and the minimum "
  "number of such contacts between two molecules to consider them as belonging "
  "to the same aggregate. Only distances between specified bead type pairs are "
  "considered; these pairs may be defined either explicitly (via --pairs and "
  "-bt options) or as all possible bead pairs (either between all bead types "
  "or bead types specified by -bt option). "
  "Information about aggregates in each timestep is written to "
  "'.agg' file (see documentation for the format of this file), and Cartesian "
  "coordinates of joined aggregates can be written to an output coordinate "
  "file (to be used for visualization or further analysis by other utilities). "
  "The utility can also differentiate between aggregates near a wall and those "
  "in bulk (-w option); when the axis perpendicular to the wall(s) and the"
  "coordinate(s) of the wall(s) along that axis are provided, any aggregate "
  "containing a <bead(s)> that is at most the contact distance from the wall "
  "is considered near the wall. Aggregates near the wall(s)/in bulk are saved "
  "into two files whose names are based on <out.agg> ('_w' and '_b' is "
  "prepended in front of the .agg extension)",

  "Usage: Aggregates <coor> <out.agg> [options]",
  .args = 2, // number of mandatory arguments
  .all = 18, // number of valid lines OptSpec (not counting last {NULL})
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<coor>", NULL, "input coordinate file", OPT_ARG},
  {"<out.agg>", NULL, "output aggregate file", OPT_ARG},
  {"-bt", "<bead(s)>", "bead types to use for closeness calculation (default: all)", OPT_EXTRA},
  {"--pairs", NULL, "-bt option specifies bead pairs instead (default: all possible pairs)", OPT_EXTRA},
  {"-d", "<float>", "maximum distance for contact (default: 1)", OPT_EXTRA},
  {"-c", "<float>", "minimum number of contacts (default: 1, max: 255)", OPT_EXTRA},
  {"-j", "<coor>", "output file with joined coordinates", OPT_EXTRA},
  {"--no_pbc", NULL, "ignore periodic boundary conditions", OPT_EXTRA},
  {"-w", "<a> <float(s)>", "coordinate(s) on <a> axis of wall(s) perpendicular to the axis", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  double cutoff;        // -d
  int contacts;         // -c
  FILE_TYPE fout;       // -j
  double wall[100];     // -w
  int w_count, axis;    // -w
  char w_file[2][LINE]; // -w
  FILE_TYPE j_file[2];  // -w (if -j)
  bool pairs;           // --pairs
  bool no_pbc;          // --no_pbc
}; //}}}

// detect possible contact between two beads //{{{
void CalculateContacts(const int id_i, const int id_j, SYSTEM System,
                       OPT opt, const double dist, int **contact,
                       const ArrNDb *use_bt_pair, SYSTEM sys_copy) {
  int i = sys_copy.BeadCoor[id_i];
  int j = sys_copy.BeadCoor[id_j];
  BEAD *b_i = &sys_copy.Bead[i];
  BEAD *b_j = &sys_copy.Bead[j];
  // skip if the the pair isn't to be used
  if (!GetArr2D(use_bt_pair, b_i->Type, b_j->Type)) {
    return;
  }
  int mol_i = b_i->Molecule;
  int mol_j = b_j->Molecule;
  vec3d *pos_i = &b_i->Position;
  vec3d *pos_j = &b_j->Position;
  vec3d rij;
  if (!opt.no_pbc) {
    // vec3d *pos_i = &b_i->Position;
    // vec3d *pos_j = &b_j->Position;
    rij = Distance(pos_i->v, pos_j->v, sys_copy.Box.Length);
  } else {
    // RemovePBCMolecule(mol_i, &System);
    // RemovePBCMolecule(mol_j, &System);
    // vec3d *pos_i = &b_i->Position;
    // vec3d *pos_j = &b_j->Position;
    for (int dd = 0; dd < 3; dd++) {
      rij.v[dd] = pos_i->v[dd] - pos_j->v[dd];
    }
  }
  rij.v[0] = VectLength(rij);
  // are 'i' and 'j' close enough?
  if (mol_i != mol_j && rij.v[0] <= dist) {
    if (mol_i > mol_j) {
      contact[mol_i][mol_j]++;
    } else {
      contact[mol_j][mol_i]++;
    }
  }
}
// structure for the callback function
struct contacts_args {
  OPT opt;
  double dist;
  int **contact;
  ArrNDb *use_bt_pair;
  SYSTEM sys_copy;
};
// adaptor for the CalculateContacts() function
static void CalculateContacts_adaptor(int id_i, int id_j,
                                      const SYSTEM System, void *ud) {
  struct contacts_args *p = (struct contacts_args*)ud;
  CalculateContacts(id_i, id_j, System, p->opt, p->dist,
                    p->contact, p->use_bt_pair, p->sys_copy);
} //}}}
// condition for using specified beads //{{{
static bool CheckBead(int id, SYSTEM System, bool *use_bt) {
  int i = System.BeadCoor[id];
  int btype = System.Bead[i].Type;
  if (!use_bt[btype] || System.Bead[i].Molecule == -1) {
    return false;
  } else {
    return true;
  }
}
// structure for the callback function (empty as only System is needed here)
struct check_args {
  bool *use_bt;
};
// adaptor for the CalculatePCF() function
static bool CheckBead_adaptor(int type, SYSTEM System, void *userdata) {
  struct check_args *p = (struct check_args*)userdata;
  return CheckBead(type, System, p->use_bt);
} //}}}

// CalculateAggregates() //{{{
// note the function doesn't fill in Aggregate[].Bead[] as it's not used here
void CalculateAggregates(AGGREGATE *Aggregate, SYSTEM *System,
                         OPT opt, bool *use_bt, ArrNDb *use_bt_pair) {
  double sqdist = Square(opt.cutoff);
  COUNT *Count = &System->Count;
  Count->Aggregate = 0;
  // zeroize - just pro-forma; is done in (Re)InitAggregate()
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
  }

  // // array for number of contacts between molecules
  // ArrNDi *contact = CreateArr2Di(Count->Molecule, Count->Molecule);
  // if (!contact) {
  //   ErrorAlloc("contact");
  // }
  // allocate & zeroize contact[][] (triangular matrix)
  int **contact = calloc(Count->Molecule, sizeof *contact);
  for (int i = 0; i < Count->Molecule; i++) {
    contact[i] = calloc(i + 1, sizeof *contact[i]);
  }

  // assign in no aggregate to each molecule
  for (int i = 0; i < Count->Molecule; i++) {
    System->Molecule[i].Aggregate = -1;
  }
  // calculate contact pairs
  double cell_size = sqrt(sqdist);
  SYSTEM copy = CopySystem(*System);
  WrapJoinCoordinates(System, true, false);
  struct contacts_args args = { opt, sqrt(sqdist), contact, use_bt_pair, copy };
  struct check_args check = { use_bt };
  TraversePairs(*System, cell_size, CalculateContacts_adaptor, &args,
                CheckBead_adaptor, &check);
  FreeSystem(&copy);

  EvaluateContacts(Aggregate, System, opt.contacts, contact);

  // FreeArrND(contact);
  // free memory
  for (int i = 0; i < Count->Molecule; i++) {
    free(contact[i]);
  }
  free(contact);

  // sort molecules in aggregates according to ascending ids //{{{
  for (int i = 0; i < System->Count.Aggregate; i++) {
    SortArray(Aggregate[i].Molecule, Aggregate[i].nMolecules, 0, 'i');
  } //}}}

  SortAggStruct(Aggregate, *System);
} //}}}

// aggregate calculation //{{{
void Calculation(SYSTEM *System, STEP step, OPT opt, COMMON_OPT commons,
                 AGGREGATE *Aggregate, char agg_file[LINE], bool *use_bt,
                 ArrNDb *use_bt_pair, const int argc, char **argv) {
  COUNT *Count = &System->Count;
  CalculateAggregates(Aggregate, System, opt, use_bt, use_bt_pair);
  // calculate & write joined coordinatest (-j option)
  if (opt.fout.name[0] != '\0') {
    FillAggregateBeads(Aggregate, *System);
    WrapJoinCoordinates(System, false, true);
    RemovePBCAggregates(opt.cutoff, Aggregate, System, use_bt);
    bool *write = calloc(Count->Bead, sizeof *write);
    if (!write) {
      ErrorAlloc("write");
    }
    InitBoolArray(write, Count->Bead, true);
    WriteTimestep(opt.fout, *System, step.coor, write, argc, argv);
    free(write);
  }

  for (int i = 0; i < Count->Aggregate; i++) {
    Aggregate[i].Flag = true;
  }
  WriteAggregates(step.coor, agg_file, *System, Aggregate);

  // are there walls (-w option)? //{{{
  if (opt.w_count > 0) {
    // find aggregates touching a wall
    for (int i = 0; i < Count->Aggregate; i++) {
      Aggregate[i].Flag = false;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_id = Aggregate[i].Molecule[j];
        MOLECULE *mol = &System->Molecule[mol_id];
        for (int k = 0; k < System->MoleculeType[mol->Type].nBeads; k++) {
          BEAD *b = &System->Bead[mol->Bead[k]];
          if (System->BeadType[b->Type].Flag) {
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
    WriteAggregates(step.coor, opt.w_file[0], *System, Aggregate);
    // reverse the Aggregate[].Flag to select aggregates in bulk
    for (int i = 0; i < Count->Aggregate; i++) {
      Aggregate[i].Flag = !Aggregate[i].Flag;
    }
    // write the aggregates to *_b.agg file
    WriteAggregates(step.coor, opt.w_file[1], *System, Aggregate);

    // write joined coordinates to _b/_w files (-j option)?
    if (opt.fout.name[0] != '\0') {
      bool *write = calloc(Count->Bead, sizeof *write);
      if (!write) {
        ErrorAlloc("write");
      }
      // assume all beads are saved (to save unbonded beads)
      InitBoolArray(write, Count->Bead, true);

      // exclude from saving all aggregate beads in bulk
      for (int i = 0; i < Count->Aggregate; i++) {
        // is aggregate in the bulk?
        if (Aggregate[i].Flag) {
          for (int j = 0; j < Aggregate[i].nMolecules; j++) {
            int mol = Aggregate[i].Molecule[j];
            int mtype = System->Molecule[mol].Type;
            for (int k = 0; k < System->MoleculeType[mtype].nBeads; k++) {
              int id = System->Molecule[mol].Bead[k];
              if (System->Bead[id].InTimestep) {
                write[id] = false;
              }
            }
          }
        }
      }
      // write joined coordinates for wall-touching aggregates to _w file
      WriteTimestep(opt.j_file[0], *System, step.coor, write, argc, argv);

      // flip write flag for all bonded beads
      for (int i = 0; i < Count->Bonded; i++) {
        int id = System->Bonded[i];
        write[id] = !write[id];
      }
      // write joined coordinates for bulk aggregates to _b file
      WriteTimestep(opt.j_file[1], *System, step.coor, write, argc, argv);
      free(write);
    }
  } //}}}

  ReInitAggregate(*System, Aggregate);
} //}}}
// structure for the callback function
struct user_data {
  OPT opt;
  COMMON_OPT commons;
  AGGREGATE *Aggregate;
  char *agg_file;
  int argc;
  char **argv;
  bool *use_bt;
  ArrNDb *use_bt_pair;
};
// adaptor for the Calculation() function
static void Calculation_adaptor(SYSTEM *System, STEP *step, void *userdata) {
  struct user_data *p = (struct user_data*)userdata;
  Calculation(System, *step, p->opt, p->commons, p->Aggregate,
              p->agg_file, p->use_bt, p->use_bt_pair, p->argc, p->argv);
};

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
  // options before reading system data
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
  opt.pairs = BoolOption(argc, argv, "--pairs");
  opt.no_pbc = BoolOption(argc, argv, "--no_pbc");
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

  // theoretically possible to use bead type?
  // Used in the CheckBead() as that checks only one bead
  bool *use_bt = calloc(Count->BeadType, sizeof *use_bt);
  // use bead type pair?
  // Used in the CalculateAggregates() to cross out unwanted pairs
  ArrNDb *use_bt_pair = CreateArr2Db(Count->BeadType, Count->BeadType);
  if (!use_bt || !use_bt_pair) {
    ErrorAlloc("use_bt/use_bt_pair");
  }
  if (opt.pairs) {
    if (!TypeOptionPair(argc, argv, "-bt", 'b', true, use_bt_pair, System)) {
      err_msg("option -bt is mandatory in this case");
      PrintErrorOption("--pairs");
      exit(1);
    }
    TypeOption(argc, argv, "-bt", 'b', true, use_bt, System);
  } else {
    if (TypeOption(argc, argv, "-bt", 'b', true, use_bt, System)) {
      for (int i = 0; i < Count->BeadType; i++) {
        for (int j = 0; j < Count->BeadType; j++) {
          if (use_bt[i] && use_bt[j]) {
            SetArr2D(use_bt_pair, i, j, true);
          } else {
            SetArr2D(use_bt_pair, i, j, false);
          }
        }
      }
    } else {
      FillArrND(use_bt_pair, true);
    }
  }
  // <bead(s)> - names of bead types to use for closeness calculation //{{{
  // if (opt.all) {
  //   for (int i = 0; i < Count->BeadType; i++) {
  //     System.BeadType[i].Flag = true;
  //   }
  // } else {
  //   bool *use_bt = calloc(Count->BeadType, sizeof *use_bt);
  //   TypeOption(argc, argv, "-bt", 'b', true, use_bt, System);
  //   for (int i = 0; i < Count->BeadType; i++) {
  //     if (use_bt[i]) {
  //       System.BeadType[i].Flag = true;
  //     } else {
  //       System.BeadType[i].Flag = false;
  //     }
  //   }
  //   free(use_bt);
  //   // // missing --all as well as any bead type(s)
  //   // // TODO: necessary to assign false? Well, Flag will not be used!
  //   // for (int i = 0; i < Count->BeadType; i++) {
  //   //   System.BeadType[i].Flag = false;
  //   // }
  //   // while (++count < argc && argv[count][0] != '-') {
  //   //   int type = FindBeadType(argv[count], System);
  //   //   if (type == -1) {
  //   //     err_msg("non-existent bead name");
  //   //     PrintError();
  //   //     ErrorBeadType(argv[count], System);
  //   //     exit(1);
  //   //   }
  //   //   if (System.BeadType[type].Flag) {
  //   //     snprintf(ERROR_MSG, LINE, "bead type %s%s%s specified more than once",
  //   //              ErrYellow(), argv[count], ErrCyan());
  //   //     PrintWarning();
  //   //   }
  //   //   System.BeadType[type].Flag = true;
  //   // }
  //   // count--; // while() always increments count at least once
  //   // if (count < (HelpDesc.args + 1)) {
  //   //   err_msg("missing <bead(s)> or --all option");
  //   //   PrintError();
  //   //   PrintCommand(stderr, argc, argv);
  //   //   Help(true, HelpDesc, opts);
  //   //   exit(1);
  //   // }
  // } //}}}

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

  STEP step = InitStep;
  struct user_data ud = { opt, commons, Aggregate, agg_file, argc, argv,
                          use_bt, use_bt_pair };
  MainLoopCoor(&System, in, commons, &step, Calculation_adaptor, &ud);

  // print last step number to <output.agg>
  // open output .agg file for appending
  FILE *fw_agg = OpenFile(agg_file, "a");
  fprintf(fw_agg, "Last Step: %d\n", step.coor);
  fclose(fw_agg);
  if (opt.w_count > 0) {
    for (int i = 0; i < 2; i++) {
      fw_agg = OpenFile(opt.w_file[i], "a");
      fprintf(fw_agg, "Last Step: %d\n", step.coor);
      fclose(fw_agg);
    }
  }

  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  free(use_bt);
  FreeArrND(use_bt_pair);

  return 0;
}
