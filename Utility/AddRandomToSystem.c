#include "../AnalysisTools.h"

// TODO: output file - any struct/output file
//       Don't forget: if separate vsf/vcf files, one more argument is mandatory
//       ...actually, what if we have FIELD struct + something else?
//       ...well, don't allow FIELD since it does generate coordinates; then
//       again, what if someone wants only that FIELD? Well, let them write it
//       themselves; it's a simple format!
//       Write structure only if struct type data, FIELD, vsf/vtf present
// TODO: not just FIELD file for specifying what to add; actually, I want FIELD
//       as it simply defines molecule 'prototypes', don't I?
//       ...switch --coordinates to use provided coordinates? That should be a
//       separate AddSystemToSystem, no? A separate AddSystemToSystem that
//       accepts any struct/coordinate file
// TODO: the '--' for input coordinate file to generate new system?
//       ...that should be a separate GenSystem, no?

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
  fprintf(ptr, "      --switch       exchange beads instead of adding them "
          "(by default, bead type with the most beads is used)\n");
  // fprintf(ptr, "      -xb <name>     unbondedbead(s) to switch new beads for "
  //         "(instead of appending to the system)\n");
  fprintf(ptr, "      -s <int>       seed for random number generator\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 8, all = common + 5, count = 0,
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
  // strcpy(option[count++], "-xb");
  strcpy(option[count++], "--switch");
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
  // lowest and/or highest distance from beads of type specified by '-bt'
  double lowest_dist = -1;
  if (DoubleOption(argc, argv, 1, "-ld", trash, &lowest_dist)) {
    exit(1);
  }
  double highest_dist = -1;
  if (DoubleOption(argc, argv, 1, "-hd", trash, &highest_dist)) {
    exit(1);
  }
  // error: if '-ld' and/or '-hd' are present, '-bt' must be too //{{{
  if (highest_dist != -1 || lowest_dist != -1) {
    bool bt = false;
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0) {
        bt = true;
      }
    }
    if (!bt) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, "if '-ld' and/or '-hd' is used,");
      fprintf(stderr, "'-bt' must be specified as well\n\n");
      // ColourReset(STDERR_FILENO);
      exit(1);
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
  // -bt <name(s)> - specify what bead types to use //{{{
  bool *bt_use_orig = calloc(S_orig.Count.BeadType, sizeof *bt_use_orig);
  if (BeadTypeOption(argc, argv, "-bt", false, bt_use_orig, &S_orig)) {
    exit(0);
  } //}}}

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
    int id0 = S_add.Molecule[i].Bead[0],
        type = S_add.Molecule[i].Type;
    VECTOR zero = S_add.Bead[id0].Position;
    for (int j = 0; j < S_add.MoleculeType[type].nBeads; j++) {
      int id = S_add.Molecule[i].Bead[j];
      S_add.Bead[id].Position.x -= zero.x;
      S_add.Bead[id].Position.y -= zero.y;
      S_add.Bead[id].Position.z -= zero.z;
    }
  } //}}}

  // set final box size //{{{
  BOX box_n,
      *box_a = &S_add.Box;
  box_n = InitBox;
  if (box_a->Length.x > S_orig.Box.Length.x) {
    box_n.Length.x = box_a->Length.x;
    box_n.Low.x = box_a->Low.x;
  } else {
    box_n.Length.x = S_orig.Box.Length.x;
    box_n.Low.x = S_orig.Box.Low.x;
  }
  if (box_a->Length.y > S_orig.Box.Length.y) {
    box_n.Length.y = box_a->Length.y;
    box_n.Low.y = box_a->Low.y;
  } else {
    box_n.Length.y = S_orig.Box.Length.y;
    box_n.Low.y = S_orig.Box.Low.y;
  }
  if (box_a->Length.z > S_orig.Box.Length.z) {
    box_n.Length.z = box_a->Length.z;
    box_n.Low.z = box_a->Low.z;
  } else {
    box_n.Length.z = S_orig.Box.Length.z;
    box_n.Low.z = S_orig.Box.Low.z;
  }
  CalculateBoxData(&box_n, 0); //}}}

  // minimize box size for adding beads if -hd is used //{{{
  BOX box_r = InitBox;
  if (highest_dist != -1) {
    VECTOR max = {0, 0, 0}, min = box_n.Length;
    for (int i = 0; i < S_orig.Count.BeadCoor; i++) {
      int id = S_orig.BeadCoor[i];
      int btype = S_orig.Bead[id].Type;
      if (bt_use_orig[btype]) {
        BEAD *bead = &S_orig.Bead[id];
        if (bead->Position.x < min.x) {
          min.x = bead->Position.x;
        }
        if (bead->Position.y < min.y) {
          min.y = bead->Position.y;
        }
        if (bead->Position.z < min.z) {
          min.z = bead->Position.z;
        }
        if (bead->Position.x > max.x) {
          max.x = bead->Position.x;
        }
        if (bead->Position.y > max.y) {
          max.y = bead->Position.y;
        }
        if (bead->Position.z > max.z) {
          max.z = bead->Position.z;
        }
      }
    }
    min.x -= highest_dist;
    min.y -= highest_dist;
    min.z -= highest_dist;
    if (min.x < 0) {
      min.x = 0;
    }
    if (min.y < 0) {
      min.y = 0;
    }
    if (min.z < 0) {
      min.z = 0;
    }
    max.x += highest_dist;
    max.y += highest_dist;
    max.z += highest_dist;
    if (max.x > box_n.Length.x) {
      max.x = box_n.Length.x;
    }
    if (max.y > box_n.Length.y) {
      max.y = box_n.Length.y;
    }
    if (max.z > box_n.Length.z) {
      max.z = box_n.Length.z;
    }
    box_r.Low.x = min.x;
    box_r.Low.y = min.y;
    box_r.Low.z = min.z;
    box_r.Length.x = max.x - min.x;
    box_r.Length.y = max.y - min.y;
    box_r.Length.z = max.z - min.z;
    CalculateBoxData(&box_r, 0);
  } else {
    box_r = box_n;
  }
  // error - no box size (should never trigger) //{{{
  if (box_n.Volume == -1) {
    strcpy(ERROR_MSG, "zero box size for the new system");
    PrintError();
    exit(1);
  } //}}}
  //}}}

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
    ConcatenateSystems(&S_new, S_add, box_n);
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
    ConcatenateSystems(&S_new, S_add, box_n);
  } //}}}

  // add monomeric beads //{{{
  int id = S_orig.Count.Bead;
  count = 0; // counts beads in sw_type's Index array (check for bonded beads)
  for (int i = 0; i < S_add.Count.Unbonded; i++) {
    VECTOR random;
    if (lowest_dist != -1 || highest_dist != -1) {
      double min_dist = 0;
      int tries = 0;
      do {
        tries++;
        double number = (double)rand() / ((double)RAND_MAX + 1);
        random.x = number * box_n.Length.x;
        // random.x = number * box_r.Length.x + box_r.Low.x;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.y = number * box_n.Length.y;
        // random.y = number * box_r.Length.y + box_r.Low.y;
        number = (double)rand() / ((double)RAND_MAX + 1);
        random.z = number * box_n.Length.z;
        // random.z = number * box_r.Length.z + box_r.Low.z;

        min_dist = SQR(box_n.Length.x) +
                   SQR(box_n.Length.y) +
                   SQR(box_n.Length.z); // some high number
        for (int j = 0; j < S_orig.Count.BeadType; j++) {
          if (bt_use_orig[j]) {
            for (int k = 0; k < S_orig.BeadType[j].Number; k++) {
              int id = S_orig.BeadType[j].Index[k];
              if (S_orig.Bead[id].InTimestep) {
                VECTOR dist;
                dist = Distance(S_orig.Bead[id].Position,
                                random, S_new.Box.Length);
                dist.x = VectorLength(dist);
                if (dist.x < min_dist) {
                  min_dist = dist.x;
                }
              }
            }
          }
        }
      } while ((lowest_dist != -1 && lowest_dist >= min_dist) ||
               (highest_dist != -1 && highest_dist <= min_dist));
    } else {
      double number = (double)rand() / ((double)RAND_MAX + 1);
      random.x = number * box_n.Length.x;
      number = (double)rand() / ((double)RAND_MAX + 1);
      random.y = number * box_n.Length.y;
      number = (double)rand() / ((double)RAND_MAX + 1);
      random.z = number * box_n.Length.z;
    }
      // printf("\nxxx %d %lf %lf %lf\n", i, random.x, random.y, random.z);

    // if (sw) {
    //   id = add_b_id_to_new_b_id[i];
    // }
    // add the new coordinate
    S_new.Bead[id].Position.x = random.x;
    S_new.Bead[id].Position.y = random.y;
    S_new.Bead[id].Position.z = random.z;
    // printf("\n%d %s\n", id, S_new.BeadType[S_new.Bead[id].Type].Name);
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
    double number = (double)rand() / ((double)RAND_MAX + 1);
    random.x = number * box_n.Length.x;
    number = (double)rand() / ((double)RAND_MAX + 1);
    random.y = number * box_n.Length.y;
    number = (double)rand() / ((double)RAND_MAX + 1);
    random.z = number * box_n.Length.z;

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
