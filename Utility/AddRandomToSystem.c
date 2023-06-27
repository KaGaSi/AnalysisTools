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
  fprintf(ptr, "   %s <input> <out.vsf> <out.vcf> [options]\n\n", cmd);

  fprintf(ptr, "      <input>        input coordinate file (vcf format)\n");
  fprintf(ptr, "      <out.vsf>      output structure file (vsf format)\n");
  fprintf(ptr, "      <out.vcf>      output coordinate file (vcf format)\n");
  fprintf(ptr, "   [general options]\n");
  fprintf(ptr, "      -f <name>            FIELD-like file with molecules "
          "to add (default: FIELD)\n");
  fprintf(ptr, "      -ld <float>          specify lowest distance from "
          "chosen bead types (default: none)\n");
  fprintf(ptr, "      -hd <float>          specify highest distance from "
          "chosen bead types (default: none)\n");
  fprintf(ptr, "      -bt <name(s)>        specify bead types new beads "
          "should be far from/near to (default: none)\n");
  CommonHelp(error, n, opt);
} //}}}

int main(int argc, char *argv[]) {

  // define options //{{{
  int common = 8, all = common + 4, count = 0,
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
  strcpy(option[count++], "-f");
  strcpy(option[count++], "-ld");
  strcpy(option[count++], "-hd");
  strcpy(option[count++], "-bt");
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

  // <output> - output vsf file //{{{
  char file_out_struct[LINE] = "";
  snprintf(file_out_struct, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  int ext = 1;
  char extension[1][EXTENSION];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(file_out_struct, ext, extension) == -1) {
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}
  // <output> - output vcf file //{{{
  char file_out_coor[LINE] = "";
  snprintf(file_out_coor, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(file_out_coor, ext, extension) == -1) {
    Help(argv[0], true, common, option);
    exit(1);
  } //}}}

  // options before reading system data //{{{
  bool silent, verbose, detailed, vtf_var;
  int timestep = 1, pbc_xyz = -1;
  int trash[1]; // some stuff for unused things in options
  CommonOptions(argc, argv, LINE, &verbose, &silent, &detailed, &vtf_var,
                &pbc_xyz, &timestep, trash, trash);
  // -f <add> - FIELD-like file with molecules to add
  char file_add_field[LINE] = "";
  if (FileIntegerOption(argc, argv, 0, "-f", trash, trash, file_add_field)) {
    exit(1);
  }
  if (file_add_field[0] == '\0') {
    strcpy(file_add_field, "FIELD");
  }
  // lowest and/or highest distance from beads of type specified by '-bt'
  double value[1] = {-1};
  if (DoubleOption(argc, argv, 1, "-ld", trash, value)) {
    exit(1);
  }
  double lowest_dist = value[0];
  value[0] = -1;
  if (DoubleOption(argc, argv, 1, "-hd", trash, value)) {
    exit(1);
  }
  double highest_dist = value[0];
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
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  SYSTEM S_orig = ReadStructure(struct_type, struct_file, coor_type, coor_file,
                                detailed, vtf_var, pbc_xyz);

  // -bt <name(s)> - specify what bead types to use //{{{
  bool *bt_use_orig = calloc(S_orig.Count.BeadType, sizeof *bt_use_orig);
  if (BeadTypeOption(argc, argv, "-bt", false, bt_use_orig, &S_orig)) {
    exit(0);
  } //}}}

  // seed random number generator
  srand(time(0));

  FILE *fr = OpenFile(coor_file, "r");
  int line_count = 0;
  for (int i = 1; i < timestep; i++) { // from 1 as timestep=1 is the first
    SkipTimestep(coor_type, fr, coor_file, struct_file, &line_count);
  }
  ReadTimestep(coor_type, fr, coor_file, &S_orig, &line_count, vtf_var);
  fclose(fr);

  // print original system (if present) //{{{
  if (verbose && strlen(coor_file) > 0) {
    fprintf(stdout, "\nORIGINAL SYSTEM\n");
    VerboseOutput(S_orig);
    if (timestep > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", timestep);
    }
  } //}}}

  SYSTEM S_add = ReadStructure(6, file_add_field, -1, "\0",
                               detailed, vtf_var, pbc_xyz);
  S_add.Count.BeadCoor = S_add.Count.Bead;
  for (int i = 0; i < S_add.Count.Bead; i++) {
    S_add.Bead[i].InTimestep = true;
    S_add.BeadCoor[i] = i;
  }

  // SYSTEM S_new;
  // InitSystem(&S_new);

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

  // minimize box size for adding beads if -hd is used
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
  PrintBox(box_n);
  PrintBox(box_r);

  // error - no box size (should never trigger) //{{{
  if (box_n.Volume == -1) {
    strcpy(ERROR_MSG, "zero box size for the new system");
    PrintError();
    exit(1);
  } //}}}

  // print what is to be added //{{{
  if (verbose) {
    fprintf(stdout, "\nBEADS AND MOLECULES TO ADD\n");
    VerboseOutput(S_add);
  } //}}}

  // join original and added systems
  // BOX Box_new = *box_n; // TODO box set up previously; maybe make it a different BOX
  SYSTEM S_new = CopySystem(S_orig);
  ConcatenateSystems(&S_new, S_add, box_n);
  PruneSystem(&S_new);
  //}}}

  // print new system //{{{
  if (verbose) {
    fprintf(stdout, "\nNEW SYSTEM\n");
    VerboseOutput(S_new);
//  PrintBondTypes2(S_new.TypesOfBonds, bond_type);
  } //}}}

  // add monomeric beads //{{{
  int id = S_orig.Count.Bead;
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
                dist = Distance(S_orig.Bead[id].Position, random, S_new.Box.Length);
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
      // printf("\n%d %lf %lf %lf\n", i, random.x, random.y, random.z);
    } else {
      double number = (double)rand() / ((double)RAND_MAX + 1);
      random.x = number * box_n.Length.x;
      number = (double)rand() / ((double)RAND_MAX + 1);
      random.y = number * box_n.Length.y;
      number = (double)rand() / ((double)RAND_MAX + 1);
      random.z = number * box_n.Length.z;
    }
      // printf("\nxxx %d %lf %lf %lf\n", i, random.x, random.y, random.z);

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
  for (int i = S_orig.Count.Molecule; i < S_new.Count.Molecule; i++) {
    int mtype = S_new.Molecule[i].Type;

    VECTOR random;
    double number = (double)rand() / ((double)RAND_MAX + 1);
    random.x = number * box_n.Length.x;
    number = (double)rand() / ((double)RAND_MAX + 1);
    random.y = number * box_n.Length.y;
    number = (double)rand() / ((double)RAND_MAX + 1);
    random.z = number * box_n.Length.z;

    for (int j = 0; j < S_new.MoleculeType[mtype].nBeads; j++) {
      id = S_new.Molecule[i].Bead[j];
      S_new.Bead[id].Position.x += random.x;
      S_new.Bead[id].Position.y += random.y;
      S_new.Bead[id].Position.z += random.z;
    }

    // print number of placed beads? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rMolecules placed: %d", i-S_orig.Count.Molecule+1);
    } //}}}
  } //}}}
  // print total number of placed beads? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                           \r");
    }
    fprintf(stdout, "\rMolecules placed: %d\n", S_add.Count.Molecule);
  } //}}}

  // write data to output files //{{{
  // vsf file
  int vsf_def_type = -1;
  count = 0;
  for (int i = 0; i < S_new.Count.BeadType; i++) {
    if (S_new.BeadType[i].Number > count) {
      count = S_new.BeadType[i].Number;
      vsf_def_type = i;
    }
  }
  WriteStructure(VSF_FILE, file_out_struct, S_new, vsf_def_type, false);

  // .vcf file
  FILE *out = OpenFile(file_out_coor, "w");
  PrintByline(out, argc, argv);
  fclose(out);
  bool *write = malloc(sizeof *write * S_new.Count.Bead);
  for (int i = 0; i < S_new.Count.Bead; i++) {
    write[i] = true; // TODO change somewhere (use different flag for sw)
  }
  // VtfWriteCoorIndexed(out, stuff, write, S_new);
  WriteTimestep(VCF_FILE, file_out_coor, S_new, 0, write);
  //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_new);
  free(bt_use_orig);
  free(write);
  //}}}

  return 0;
}
