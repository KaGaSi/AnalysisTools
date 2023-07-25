#include "../AnalysisTools.h"

// TODO: the '--' for input coordinate file to generate new system?
//       ...that should be a separate GenSystem, no?
// TODO: --bonded option as a 'shorthand' for -bt <list> <all> <bonded> <btypes>
// TODO: -cx|y|z <min> <max> to constrain coordinates for adding beads
// TODO: rotate the molecules randomly (+ switch --no-rot)

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
AddToSystem either creates a system from scratch or adds unbonded beads \
and/or molecules to an existing system. The new components are defined either \
by a FIELD-like file (an input file for DL_MESO simulation program) or by \
vsf/vcf files (-vtf option). In the first case, new components are placed \
randomly (with several possible constraints), while in the second case, the \
provided coordinates are used as is.\n\n");
  }
  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <in.field> <output> [options]\n\n", cmd);

  fprintf(ptr, "      <input>        input coordinate file (vcf format)\n");
  fprintf(ptr, "      <in.field>     input FIELD file with species to add\n");
  fprintf(ptr, "      <output>       output structure and coordinate file "
          "(format: xyz, lammpstrj, or vtf)\n");
  fprintf(ptr, "   [general options]\n");
  fprintf(ptr, "      -ld <float>    specify lowest distance from "
          "chosen bead types (default: none)\n");
  fprintf(ptr, "      -hd <float>    specify highest distance from "
          "chosen bead types (default: none)\n");
  fprintf(ptr, "      -bt <name(s)>  specify bead types new beads "
          "should be far from/near to (default: none)\n");
  fprintf(ptr, "      --bonded       use bonded beads for the distance "
          "condition (overwrites -bt option)\n");
  fprintf(ptr, "      --switch       exchange beads instead of adding them "
          "(by default, bead type with the most beads is used)\n");
  fprintf(ptr, "      -cx 2×<float>  constrain x-coordinate to specified "
          "dimensions\n");
  fprintf(ptr, "      -cy 2×<float>  constrain y-coordinate to specified "
          "dimensions\n");
  fprintf(ptr, "      -cz 2×<float>  constrain z-coordinate to specified "
          "dimensions\n");
  fprintf(ptr, "      -s <int>       seed for random number generator\n");
  CommonHelp(error, n, opt);
} //}}}

// generate random point in a cube (0,length)^3 //{{{
VECTOR RandomCoordinate(BOX box) {
  VECTOR random;
  double number = (double)rand() / ((double)RAND_MAX + 1);
  random.x = number * box.Length.x + box.Low.x;
  number = (double)rand() / ((double)RAND_MAX + 1);
  random.y = number * box.Length.y + box.Low.y;
  number = (double)rand() / ((double)RAND_MAX + 1);
  random.z = number * box.Length.z + box.Low.z;
  return random;
} //}}}

// generate random point constrained by distance from other beads //{{{
VECTOR RandomConstrainedCoor(BOX constrain, SYSTEM S_orig,
                             bool bt_use_orig[], bool bonded, VECTOR box,
                             double lowest_dist, double highest_dist) {
  VECTOR random;
  double min_dist = 0;
  do {
    random = RandomCoordinate(constrain);
    min_dist = SQR(constrain.Length.x + constrain.Low.x) +
               SQR(constrain.Length.y + constrain.Low.y) +
               SQR(constrain.Length.z + constrain.Low.z); // a high number
    if (bonded) { // use all bonded beads
      for (int i = 0; i < S_orig.Count.BondedCoor; i++) {
        int id = S_orig.BondedCoor[i];
        VECTOR dist = Distance(S_orig.Bead[id].Position, random, box);
        dist.x = VectorLength(dist);
        if (dist.x < min_dist) {
          min_dist = dist.x;
        }
      }
    } else { // use specified bead types
      for (int i = 0; i < S_orig.Count.BeadType; i++) {
        if (bt_use_orig[i]) {
          for (int j = 0; j < S_orig.BeadType[i].Number; j++) {
            int id = S_orig.BeadType[i].Index[j];
            if (S_orig.Bead[id].InTimestep) {
              VECTOR dist = Distance(S_orig.Bead[id].Position, random, box);
              dist.x = VectorLength(dist);
              if (dist.x < min_dist) {
                min_dist = dist.x;
              }
            }
          }
        }
      }
    }
  } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
           (highest_dist != -1 && highest_dist <= min_dist));
  return random;
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 8, all = common + 9, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  // common options
  strcpy(option[count++], "-st");
  strcpy(option[count++], "-i");
  // strcpy(option[count++], "--variable"); // TODO: makes no sense, I think
  strcpy(option[count++], "-pbc");
  strcpy(option[count++], "--detailed");
  strcpy(option[count++], "--verbose");
  strcpy(option[count++], "--silent");
  strcpy(option[count++], "--help");
  strcpy(option[count++], "--version");
  // extra options
  strcpy(option[count++], "-ld");
  strcpy(option[count++], "-hd");
  strcpy(option[count++], "-bt");
  strcpy(option[count++], "--bonded");
  // strcpy(option[count++], "-xb");
  strcpy(option[count++], "--switch");
  strcpy(option[count++], "-cx");
  strcpy(option[count++], "-cy");
  strcpy(option[count++], "-cz");
  strcpy(option[count++], "-s");
  OptionCheck(argc, argv, req_arg, common, all, option);
  //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate (and structure) file //{{{
  char coor_file[LINE] = "", struct_file[LINE] = "";
  int coor_type = -1, struct_type = 0;
  snprintf(coor_file, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, coor_file, &coor_type,
                       struct_file, &struct_type)) {
    exit(1);
  } //}}}

  // <in.field> - FIELD file with specis to add //{{{
  char field_add_file[LINE] = "";
  snprintf(field_add_file, LINE, "%s", argv[++count]);
  int field_add_type = StructureFileType(field_add_file, 0);
  if (field_add_type != FIELD_FILE) {
    strcpy(ERROR_MSG, "input FIELD file required");
    PrintErrorFile(field_add_file, "\0", "\0");
    exit(1);
  } //}}}

  // <output> - coordinate and structure output file //{{{
  char out_file[LINE] = "";
  snprintf(out_file, LINE, "%s", argv[++count]);
  int out_type = FullFileType(out_file, 1);
  if (out_type == -1) {
    strcpy(ERROR_MSG, "output file must be lammpstrj, xyz, or vtf");
    PrintErrorFile(out_file, "\0", "\0");
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int timestep = 1, pbc_xyz = -1;
  int trash[1]; // some stuff for unused things in options
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &timestep, trash, trash);
  // lowest and/or highest distance from specified beads //{{{
  double lowest_dist = -1, highest_dist = -1;
  DoubleOption1(argc, argv, "-ld", &lowest_dist);
  DoubleOption1(argc, argv, "-hd", &highest_dist);
  if (lowest_dist != -1 && highest_dist != -1 && lowest_dist >= highest_dist) {
    strcpy(ERROR_MSG, "highest distance must be higher than lowest distance");
    PrintErrorOption("-ld/-hd");
    Help(argv[0], true, common, option);
    exit(1);
  }
  // error: missing '-bt'/'--bonded' when '-ld' and/or '-hd' are used //{{{
  if (highest_dist != -1 || lowest_dist != -1) {
    bool bt = false;
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0 || strcmp(argv[i], "--bonded") == 0) {
        bt = true;
      }
    }
    if (!bt) {
      strcpy(ERROR_MSG, "missing mandatory -bt or --bonded options");
      PrintErrorOption("-ld/-hd");
      Help(argv[0], true, common, option);
      exit(1);
    }
  } //}}}
  //}}}
  // axes constraints (-cx/y/z options) //{{{
  double cx[2] = {-1, -1}, cy[2] = {-1, -1}, cz[2] = {-1, -1};
  if (DoubleOption2(argc, argv, "-cx", cx)) {
    if (cx[0] == cx[1]) {
      strcpy(ERROR_MSG, "two different distance values required");
      PrintErrorOption("-cx");
      exit(1);
    } else if (cx[0] > cx[1]) {
      SwapDouble(&cx[0], &cx[1]);
    }
  }
  if (DoubleOption2(argc, argv, "-cy", cy)) {
    if (cy[0] == cy[1]) {
      strcpy(ERROR_MSG, "two different distance values required");
      PrintErrorOption("-cy");
      exit(1);
    } else if (cy[0] > cy[1]) {
      SwapDouble(&cy[0], &cy[1]);
    }
  }
  if (DoubleOption2(argc, argv, "-cz", cz)) {
    if (cz[0] == cz[1]) {
      strcpy(ERROR_MSG, "two different distance values required");
      PrintErrorOption("-cz");
      exit(1);
    } else if (cz[0] > cz[1]) {
      SwapDouble(&cz[0], &cz[1]);
    }
  } //}}}

  // seed for random number generator (-s option)
  int seed = -1;
  if (IntegerOption(argc, argv, 1, "-s", trash, &seed)) {
    exit(1);
  }
  // exchange beads instead of appending them?
  bool sw = BoolOption(argc, argv, "--switch");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  SYSTEM S_orig = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);

  // find bead type to switch (the most numerous one; solvent, presumably) //{{{
  int sw_type = -1; // type to switch (if --switch is used)
  if (sw) {
    count = 0;
    for (int i = 0; i < S_orig.Count.BeadType; i++) {
      if (S_orig.BeadType[i].Number > count) {
        count = S_orig.BeadType[i].Number;
        sw_type = i;
      }
    }
  } //}}}
  // -bt <name(s)>/--bonded - specify what bead types to use //{{{
  bool *bt_use_orig = calloc(S_orig.Count.BeadType, sizeof *bt_use_orig);
  if (BeadTypeOption(argc, argv, "-bt", false, bt_use_orig, &S_orig)) {
    exit(0);
  }
  bool bonded = BoolOption(argc, argv, "--bonded"); //}}}

  // seed random number generator //{{{
  if (seed != -1) {
    srand(seed);
  } else {
    srand(time(0));
  } //}}}

  // read input coordinates //{{{
  FILE *fr = OpenFile(coor_file, "r");
  int line_count = 0;
  if (coor_type != LDATA_FILE) {
    for (int i = 1; i < timestep; i++) { // from 1 as timestep=1 is the first
      SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count);
    }
  }
  ReadTimestep(coor_type, fr, coor_file, &S_orig, &line_count, vtf_var);
  fclose(fr); //}}}

  // print original system (if present) //{{{
  if (verbose && strlen(coor_file) > 0) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(S_orig);
    if (timestep > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", timestep);
    }
  } //}}}

  SYSTEM S_add = ReadStructure(6, field_add_file, -1, "\0",
                               detailed, vtf_var, pbc_xyz);
  S_add.Count.BeadCoor = S_add.Count.Bead;
  for (int i = 0; i < S_add.Count.Bead; i++) {
    S_add.Bead[i].InTimestep = true;
    S_add.BeadCoor[i] = i;
  }
  // minimize initial coordinates of added molecules //{{{
  for (int i = 0; i < S_add.Count.Molecule; i++) {
    int type = S_add.Molecule[i].Type;
    // int id0 = S_add.Molecule[i].Bead[0];
    // VECTOR zero = S_add.Bead[id0].Position;
    VECTOR zero = GeomCentre(S_add.MoleculeType[type].nBeads, S_add.Molecule[i].Bead, S_add.Bead);
    for (int j = 0; j < S_add.MoleculeType[type].nBeads; j++) {
      int id = S_add.Molecule[i].Bead[j];
      S_add.Bead[id].Position.x -= zero.x;
      S_add.Bead[id].Position.y -= zero.y;
      S_add.Bead[id].Position.z -= zero.z;
    }
  } //}}}

  // set final box size //{{{
  BOX box_new,
      *box_add = &S_add.Box;
  box_new = InitBox;
  if (box_add->Length.x > S_orig.Box.Length.x) {
    box_new.Length.x = box_add->Length.x;
    box_new.Low.x = box_add->Low.x;
  } else {
    box_new.Length.x = S_orig.Box.Length.x;
    box_new.Low.x = S_orig.Box.Low.x;
  }
  if (box_add->Length.y > S_orig.Box.Length.y) {
    box_new.Length.y = box_add->Length.y;
    box_new.Low.y = box_add->Low.y;
  } else {
    box_new.Length.y = S_orig.Box.Length.y;
    box_new.Low.y = S_orig.Box.Low.y;
  }
  if (box_add->Length.z > S_orig.Box.Length.z) {
    box_new.Length.z = box_add->Length.z;
    box_new.Low.z = box_add->Low.z;
  } else {
    box_new.Length.z = S_orig.Box.Length.z;
    box_new.Low.z = S_orig.Box.Low.z;
  }
  CalculateBoxData(&box_new, 0); //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(S_add);
  } //}}}

  // create output System //{{{
  SYSTEM S_new;
  COUNT *C_orig = &S_orig.Count;
  COUNT *C_new = &S_new.Count;
  COUNT *C_add = &S_add.Count;
  int *add_b_id_to_new_b_id; // array of orig_ids to switch (in case of --switch)
  // if not switched, concatenate the new (i.e., original) and the added systems
  if (!sw) { // do not switch, append the new system
    S_new = CopySystem(S_orig);
    ConcatenateSystems(&S_new, S_add, box_new);
  } else { // switch, so transform the system
    // error - too few beads to switch //{{{
    if (C_add->Bead > S_orig.BeadType[sw_type].Number) {
      if (snprintf(ERROR_MSG, LINE, "not enough %s%s%s beads to switch",
                   ErrYellow(), S_orig.BeadType[sw_type].Name, ErrRed()) < 0) {
        strcpy(ERROR_MSG, "something wrong with snprintf()");
        PrintError();
        exit(1);
      }
      PrintError();
      exit(1);
    } //}}}
    add_b_id_to_new_b_id = calloc(C_add->Bead, sizeof *add_b_id_to_new_b_id);
    count = S_orig.BeadType[sw_type].Number - 1;
    for (int i = 0; i < C_add->Bead; i++) {
      int id = S_orig.BeadType[sw_type].Index[count];
      S_orig.Bead[id].InTimestep = false;
      count--;
    }
    PruneSystem(&S_orig);
    S_new = CopySystem(S_orig);
    ConcatenateSystems(&S_new, S_add, box_new);
  } //}}}

  // define constrained box for adding beads (-cx/y/z and/or -hd options) //{{{
  BOX box_constrain = InitBox;
  box_constrain.Length = S_new.Box.Length;
  if (cx[0] != -1) {
    box_constrain.Low.x = cx[0];
    box_constrain.Length.x = cx[1] - cx[0];
  }
  if (cy[0] != -1) {
    box_constrain.Low.y = cy[0];
    box_constrain.Length.y = cy[1] - cy[0];
  }
  if (cz[0] != -1) {
    box_constrain.Low.z = cz[0];
    box_constrain.Length.z = cz[1] - cz[0];
  }
  CalculateBoxData(&box_constrain, 0);
  // minimize box if -hd is used
  if (highest_dist != -1) {
    // find minimum/maximum coordinates of beads for distance check //{{{
    VECTOR max = {0, 0, 0}, min = S_orig.Box.Length;
    if (bonded) { // use all bonded beads
      for (int i = 0; i < S_orig.Count.BondedCoor; i++) {
        int id = S_orig.BondedCoor[i];
        BEAD *b = &S_orig.Bead[id];
        if (b->InTimestep) {
          if (b->Position.x < min.x) {
            min.x = b->Position.x;
          } else if (b->Position.x > max.x) {
            max.x = b->Position.x;
          }
          if (b->Position.y < min.y) {
            min.y = b->Position.y;
          } else if (b->Position.y > max.y) {
            max.y = b->Position.y;
          }
          if (b->Position.z < min.z) {
            min.z = b->Position.z;
          } else if (b->Position.z > max.z) {
            max.z = b->Position.z;
          }
        }
      }
    } else { // use bead types specified by -bt
      for (int i = 0; i < S_orig.Count.BeadType; i++) {
        if (bt_use_orig[i]) {
          for (int j = 0; j < S_orig.BeadType[i].Number; j++) {
            int id = S_orig.BeadType[i].Index[j];
            BEAD *b = &S_orig.Bead[id];
            if (b->InTimestep) {
              if (b->Position.x < min.x) {
                min.x = b->Position.x;
              } else if (b->Position.x > max.x) {
                max.x = b->Position.x;
              }
              if (b->Position.y < min.y) {
                min.y = b->Position.y;
              } else if (b->Position.y > max.y) {
                max.y = b->Position.y;
              }
              if (b->Position.z < min.z) {
                min.z = b->Position.z;
              } else if (b->Position.z > max.z) {
                max.z = b->Position.z;
              }
            }
          }
        }
      }
    } //}}}
    // get the maximum possible coordinate of any added bead //{{{
    max.x += highest_dist;
    max.y += highest_dist;
    max.z += highest_dist;
    if (max.x > (box_constrain.Low.x + box_constrain.Length.x)) {
      max.x = box_constrain.Low.x + box_constrain.Length.x;
    }
    if (max.y > (box_constrain.Low.y + box_constrain.Length.y)) {
      max.y = box_constrain.Low.y + box_constrain.Length.y;
    }
    if (max.z > (box_constrain.Low.z + box_constrain.Length.z)) {
      max.z = box_constrain.Low.z + box_constrain.Length.z;
    } //}}}
    // get the minimum possible coordinate of any added bead //{{{
    min.x -= highest_dist;
    min.y -= highest_dist;
    min.z -= highest_dist;
    if (min.x < box_constrain.Low.x) {
      min.x = box_constrain.Low.x;
    }
    if (min.y < box_constrain.Low.y) {
      min.y = box_constrain.Low.y;
    }
    if (min.z < box_constrain.Low.z) {
      min.z = box_constrain.Low.z;
    } //}}}
    // define the box
    box_constrain.Length.x = max.x - min.x;
    box_constrain.Length.y = max.y - min.y;
    box_constrain.Length.z = max.z - min.z;
    box_constrain.Low = min;
    CalculateBoxData(&box_constrain, 0);
  } //}}}

  // add monomeric beads //{{{
  int id = S_orig.Count.Bead;
  count = 0; // counts beads in sw_type's Index array (check for bonded beads)
  for (int i = 0; i < S_add.Count.Unbonded; i++) {
    VECTOR random;
    if (lowest_dist != -1 || highest_dist != -1) {
      random = RandomConstrainedCoor(box_constrain, S_orig,
                                     bt_use_orig, bonded, S_new.Box.Length,
                                     lowest_dist, highest_dist);
    } else {
      random = RandomCoordinate(box_constrain);
    }

    S_new.Bead[id].Position.x = random.x;
    S_new.Bead[id].Position.y = random.y;
    S_new.Bead[id].Position.z = random.z;
    id++;

    // print number of placed beads? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rMonomers placed: %d", i+1);
    } //}}}
  } //}}}
  // print total number of placed beads? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                           \r");
    }
    fprintf(stdout, "\rMonomer placed: %d\n", S_add.Count.Unbonded);
  } //}}}

  // add molecules //{{{
  for (int i = C_orig->Molecule; i < C_new->Molecule; i++) {
    int mtype = S_new.Molecule[i].Type;

    VECTOR random;
    if (lowest_dist != -1 || highest_dist != -1) {
      random = RandomConstrainedCoor(box_constrain, S_orig,
                                     bt_use_orig, bonded, S_new.Box.Length,
                                     lowest_dist, highest_dist);
    } else {
      random = RandomCoordinate(box_constrain);
    }

    for (int j = 0; j < S_new.MoleculeType[mtype].nBeads; j++) {
      int id_add = S_add.Molecule[i-C_orig->Molecule].Bead[j];
      id = S_new.Molecule[i].Bead[j];
      S_new.Bead[id].Position.x = S_add.Bead[id_add].Position.x + random.x;
      S_new.Bead[id].Position.y = S_add.Bead[id_add].Position.y + random.y;
      S_new.Bead[id].Position.z = S_add.Bead[id_add].Position.z + random.z;
    }

    // print number of placed molecules? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rMolecules placed: %d", i-S_orig.Count.Molecule+1);
    } //}}}
  } //}}}
  // print total number of placed molecules? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                           \r");
    }
    fprintf(stdout, "\rMolecules placed: %d\n", S_add.Count.Molecule);
  } //}}}

  // remove possible duplicities from the system - a fiddly way that works //{{{
  C_new->BeadCoor = C_new->Bead;
  for (int i = 0; i < C_new->BeadType; i++) {
    free(S_new.BeadType[i].Index); // MergeBeadTypes() requires it
  }
  MergeBeadTypes(&S_new, true);
  FillBeadTypeIndex(&S_new); // re-fill after MergeBeadTypes()
  FreeSystem(&S_orig); // to copy to for pruning of various parameters
  S_orig = CopySystem(S_new);
  PruneBondTypes(S_orig, &S_new);
  PruneAngleTypes(S_orig, &S_new);
  PruneDihedralTypes(S_orig, &S_new);
  PruneImproperTypes(S_orig, &S_new);
  // free stuff that's required unalloc'd by MergeMoleculeTypes()
  for (int i = 0; i < C_new->MoleculeType; i++) {
    free(S_new.MoleculeType[i].BType);
    free(S_new.MoleculeType[i].Index);
  }
  MergeMoleculeTypes(&S_new);
  // re-fill the freed stuff
  for (int i = 0; i < C_new->MoleculeType; i++) {
    FillMoleculeTypeBType(&S_new.MoleculeType[i]);
  }
  FillMoleculeTypeIndex(&S_new);
  // prune the new system just for sure - I'm not sure what it does
  PruneSystem(&S_new); //}}}

  // minimize Molecule indices //{{{
  for (int i = 0; i < C_new->Molecule; i++) {
    S_new.Molecule[i].Index = i;
  }
  C_new->HighestResid = C_new->Molecule; //}}}

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(S_new);
  } //}}}

  // write data to output files //{{{
  // vsf file
  char open[2] = "w";
  if (out_type == VTF_FILE) {
    int vsf_def_type = -1;
    count = 0;
    for (int i = 0; i < S_new.Count.BeadType; i++) {
      if (S_new.BeadType[i].Number > count) {
        count = S_new.BeadType[i].Number;
        vsf_def_type = i;
      }
    }
    WriteStructure(out_type, out_file, S_new, vsf_def_type, false);
    open[0] = 'a';
  }
  FILE *out = OpenFile(out_file, open);
  // set all beads to be written
  bool *write = malloc(sizeof *write * S_new.Count.Bead);
  for (int i = 0; i < S_new.Count.Bead; i++) {
    S_new.BeadCoor[i] = i;
    write[i] = true;
  }
  WriteTimestep(out_type, out_file, S_new, 0, write);
  fclose(out);
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_new);
  free(bt_use_orig);
  free(write);
  if (sw) {
    free(add_b_id_to_new_b_id);
  }
  //}}}

  return 0;
}
