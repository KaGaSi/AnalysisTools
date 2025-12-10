#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AddToSystem either creates a system from scratch or adds unbonded beads "
  "and/or molecules to an existing system. The new components are defined "
  "by a FIELD-like file and are placed either randomly or according to "
  "several possible constraints. The new species can either be appended "
  "to the system, or specified beads can be exchanged for the new ones.",

  "Usage: AddToSystem <input> <in.field> <output> [options]",
  .args = 3, // number of mandatory arguments
  .all = 29, // number of valid lines OptSpec (not counting last {NULL})
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
  {"<input>", NULL, "input coordinate file or '-' to generate new system", OPT_ARG},
  {"<in.field>", NULL, "input FIELD file with beads to add", OPT_ARG},
  {"<output>", NULL, "output coordinate file", OPT_ARG},
  {"-o", "<filename>", "output extra structure file", OPT_EXTRA},
  {"-ld", "<float>", "lowest distance from chosen beads (default: none)", OPT_EXTRA},
  {"-hd", "<float>", "highest distance from chosen beads (default: none)", OPT_EXTRA},
  {"-bt", "<name(s)>", "bead types for -hd/-ld (default: none)", OPT_EXTRA},
  {"--bonded", NULL, "use bonded beads for -hd/-ld (overwrites -bt option)", OPT_EXTRA},
  {"-xb", "<bead type>", "what bead type to exchange", OPT_EXTRA},
  {"--add", NULL, "add beads instead of exchanging (overwrites -xb)", OPT_EXTRA},
  {"--no-rotate", NULL, "do not randomly rotate molecules", OPT_EXTRA},
  {"-a", "3x<angle>", "rotate molecules by yaw, pitch, and roll in degrees (overrides --no-rotate)", OPT_EXTRA},
  {"-cx", "2x<float>", "constrain x-coordinate (in fraction of output box)", OPT_EXTRA},
  {"-cy", "2x<float>", "constrain y-coordinate (in fraction of output box)", OPT_EXTRA},
  {"-cz", "2x<float>", "constrain z-coordinate (in fraction of output box)", OPT_EXTRA},
  {"--tail", NULL, "use molecule's last bead for constraint checks (default: molecule's geometric centre)", OPT_EXTRA},
  {"--head", NULL, "use molecule's first bead for constraint checks (overrides --tail)", OPT_EXTRA},
  {"--real", NULL, "use real coordinates for-cx/-cy/-cz/-off options", OPT_EXTRA},
  {"-b", "<x> <y> <z>", "new box dimensions (in real units)", OPT_EXTRA},
  {"-off", "3x<float>", "orignial system's offset (in fractions of the output box)", OPT_EXTRA},
  {"-s", "<int>", "seed for random number generator", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  bool ld, hd;             // -ld/-hd
  double ldist, hdist,     //
         axis[3][2];       // -cx/-cy/-cz
  vec3d angle[3],          // -a
        off[3];            // -off
  bool *bt_use_orig,       // -bt
       *sw_type,           // -xb
       new,                // generate new system from scratch?
       real, add, no_rot,  // --real/--add/--no-rotate
       bonded, head, tail; // --bonded/--head/--tail
  BOX box;                 // -b (then constrained 'box' via -cx/-cy/-cz)
  int seed;                // -s
  FILE_TYPE fout;          // -o
}; //}}}

// generate random point in a cube (0,length)^3 //{{{
vec3d RandomCoordinate(BOX box) {
  vec3d random;
  for (int dd = 0; dd < 3; dd++) {
    double number = (double)(rand()) / ((double)(RAND_MAX) + 1);
    random.v[dd] = number * box.Length.v[dd] + box.Low.v[dd];
  }
  return random;
} //}}}

// generate random point constrained by distance from other beads //{{{
/* What beads to use for distance check (int mode)
 *   0...no checks
 *   1...all bonded beads
 *   2...specified bead types,
 */
void GetMinDist(BEAD bead, vec3d random, vec3d box, double *min_dist) {
  vec3d dist = Distance(bead.Position.v, random.v, box);
  dist.v[0] = VectLength(dist);
  if (dist.v[0] < *min_dist) {
    *min_dist = dist.v[0];
  }
}
vec3d RandomConstrainedCoor(SYSTEM S_orig, int mode, vec3d box, OPT opt) {
  vec3d random;
  if (mode == 0) { // no distance check
    for (int dd = 0; dd < 3; dd++) {
      random = RandomCoordinate(opt.box);
    }
    return random;
  }
  COUNT *C_orig = &S_orig.Count;
  double min_dist = 0;
  do {
    random = RandomCoordinate(opt.box);
    min_dist = 1e6;  // simply a high number
    if (mode == 1) { // use all bonded beads
      for (int i = 0; i < C_orig->BondedCoor; i++) {
        int id = S_orig.BondedCoor[i];
        GetMinDist(S_orig.Bead[id], random, box, &min_dist);
      }
    } else if (mode == 2) { // use specified bead types
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (opt.bt_use_orig[i]) {
          for (int j = 0; j < S_orig.BeadType[i].Number; j++) {
            int id = S_orig.BeadType[i].Index[j];
            if (S_orig.Bead[id].InTimestep) {
              GetMinDist(S_orig.Bead[id], random, box, &min_dist);
            }
          }
        }
      }
    // TODO: this is not yet implemented - I think...
    } else if (mode == 3) { // use first bead of each molecule
      for (int i = 0; i < C_orig->Molecule; i++) {
        int id = S_orig.Molecule[i].Bead[0];
        if (S_orig.Bead[id].InTimestep) {
          GetMinDist(S_orig.Bead[id], random, box, &min_dist);
        }
      }
    } else {
      err_msg("RandomConstrainedCoor(): mode must be 0 to 3");
      PrintError();
      exit(1);
    }
  } while ((opt.ld && opt.ldist >= min_dist) ||
           (opt.hd && opt.hdist <= min_dist));
  return random;
} //}}}

// rotate randomly given collection of beads (e.g., a molecule) //{{{
void Rotate(SYSTEM System, int number, const int *list,
            const double rot_angle[3], double (*new)[3]) {
  // rotation angles around x-, y-, and z-axes
  double alpha, beta, gamma;
  // specified by -a option...
  if (rot_angle[0] != 0 || rot_angle[1] != 0 || rot_angle[2] != 0) {
    alpha = rot_angle[0] / 180 * PI;
    beta  = rot_angle[1] / 180 * PI;
    gamma = rot_angle[2] / 180 * PI;
  // ...or random
  } else {
    alpha = (double)(rand()) / (double)(RAND_MAX) * PI;
    beta  = (double)(rand()) / (double)(RAND_MAX) * PI;
    gamma = (double)(rand()) / (double)(RAND_MAX) * PI;
  }
  double rot[3][3];
  rot[0][0] = cos(alpha) * cos(beta);
  rot[1][0] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
  rot[2][0] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);

  rot[0][1] = sin(alpha) * cos(beta);
  rot[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
  rot[2][1] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);

  rot[0][2] = -sin(beta);
  rot[1][2] = cos(beta) * sin(gamma);
  rot[2][2] = cos(beta) * cos(gamma);
  // generate the rotated coordinates
  for (int i = 0; i < number; i++) {
    for (int dd = 0; dd < 3; dd++) {
      new[i][dd] = rot[dd][0] * System.Bead[list[i]].Position.v[0] +
                   rot[dd][1] * System.Bead[list[i]].Position.v[1] +
                   rot[dd][2] * System.Bead[list[i]].Position.v[2];
    }
  }
} //}}}

int main(int argc, char *argv[]) {

  // commad line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  opt.new = true; // create new system from scratch?
  if (argv[++count][0] != '-') {
    s_strcpy(in.coor.name, argv[count], LINE);
    opt.new = false;
    if (!InputCoorStruct(argc, argv, &in)) {
      exit(1);
    }
  } //}}}
  // <in.field> - FIELD file with specis to add //{{{
  SYS_FILES field = InitSysFiles;
  s_strcpy(field.stru.name, argv[++count], LINE);
  field.stru.type = StructureFileType(field.stru.name);
  if (field.stru.type != FIELD_FILE) {
    err_msg("input FIELD file required");
    PrintErrorFile(field.stru.name, "\0", "\0");
    exit(1);
  } //}}}
  // <output> - coordinate and structure output file //{{{
  FILE_TYPE fout = InitFile;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name); //}}}

  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // output structure file (-o option)
  opt.fout.name[0] = '\0';
  if (FileOption(argc, argv, "-o", opt.fout.name)) {
    opt.fout.type = FileType(opt.fout.name);
  }
  // lowest and/or highest distance from specified beads //{{{
  opt.ld = false;
  opt.hd = false;
  if (!opt.new) { // only if not generating system from scratch
    opt.ld = OneNumberOption(argc, argv, "-ld", &opt.ldist, 'd');
    opt.hd = OneNumberOption(argc, argv, "-hd", &opt.hdist, 'd');
  }
  // errors for -ld/-hd options //{{{
  if ((opt.ld && opt.ldist <= 0) || (opt.hd && opt.hdist <= 0)) {
    err_msg("highest/lowest distance must be positive real number");
    PrintErrorOption("-ld/-hd");
    PrintCommand(stderr, argc, argv);
    Help(true, HelpDesc, opts);
    exit(1);
  }
  if (opt.ld && opt.hd && opt.ldist >= opt.hdist) {
    err_msg("highest distance must be higher than lowest distance");
    PrintErrorOption("-ld/-hd");
    PrintCommand(stderr, argc, argv);
    Help(true, HelpDesc, opts);
    exit(1);
  }
  if (opt.hd || opt.ld) {
    bool bt = false;
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0 || strcmp(argv[i], "--bonded") == 0) {
        bt = true;
        break;
      }
    }
    if (!bt) {
      err_msg("missing mandatory -bt or --bonded options");
      PrintErrorOption("-ld/-hd");
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}
  //}}}
  opt.real = BoolOption(argc, argv, "--real");
  // axes constraints (-cx/y/z options) //{{{
  for (int dd = 0; dd < 3; dd++) {
    opt.axis[dd][0] = -1;
    opt.axis[dd][1] = -1;
    char str[4];
    switch (dd) {
      case 0:
        s_strcpy(str, "-cx", 4);
        break;
      case 1:
        s_strcpy(str, "-cy", 4);
        break;
      case 2:
        s_strcpy(str, "-cz", 4);
        break;
    }
    // TODO: should be able to be negative, right? Consider, e.g., ltrj BOX...
    if (TwoNumbersOption(argc, argv, str, opt.axis[dd], 'd')) {
      if (opt.axis[dd][0] < 0 || opt.axis[dd][1] < 0) {
        err_msg("two non-negative numbers required");
        PrintErrorOption("-cx/-cy/-cz");
        exit(1);
      } else if (opt.axis[dd][0] == opt.axis[dd][1]) {
        err_msg("two different distance values required");
        PrintErrorOption("-cx/-cy/-cz");
        exit(1);
      } else if (opt.axis[dd][0] > opt.axis[dd][1]) {
        SwapDouble(&opt.axis[dd][0], &opt.axis[dd][1]);
      }
    }
    if (!opt.real) {
      if ((opt.axis[dd][0] != -1 && opt.axis[dd][0] > 1) ||
          (opt.axis[dd][1] != -1 && opt.axis[dd][1] > 1)) {
        err_msg("unless --real is used, -cx/y/z must be between 0 and 1");
        PrintErrorOption(str);
        exit(1);
      }
    } //}}}
  }
  // exchange beads instead of appending them?
  opt.add = BoolOption(argc, argv, "--add");
  // always add, if generating the system from scratch
  if (opt.new) {
    opt.add = true;
  }
  // do not rotate molecules?
  opt.no_rot = BoolOption(argc, argv, "--no-rotate");
  // output box dimensions //{{{
  InitDoubleArray(opt.angle->v, 3, 0);
  if (ThreeNumbersOption(argc, argv, "-a", opt.angle->v, 'd')) {
    opt.no_rot = false;
  } //}}}
  opt.head = BoolOption(argc, argv, "--head");
  opt.tail = BoolOption(argc, argv, "--tail");
  // output box dimensions //{{{
  opt.box = InitBox;
  vec3d temp = { .v = {0, 0, 0}};
  if (ThreeNumbersOption(argc, argv, "-b", temp.v, 'd')) {
    opt.box.Length = temp;
    if (count != 3 ||
        opt.box.Length.x <= 0 ||
        opt.box.Length.y <= 0 ||
        opt.box.Length.z <= 0) {
      err_msg("three positive numbers required");
      PrintErrorOption("-b");
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}
  // -off option
  InitDoubleArray(opt.off->v, 3, 0);
  ThreeNumbersOption(argc, argv, "-off", opt.off->v, 'd');
  // seed for random number generator (-s option)
  opt.seed = -1;
  OneNumberOption(argc, argv, "-s", &opt.seed, 'i');
  // warn about options with no effect //{{{
  if (opt.new) {
    for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-bt") == 0 ||
          strcmp(argv[i], "-ld") == 0 ||
          strcmp(argv[i], "-hd") == 0 ||
          strcmp(argv[i], "--bonded") == 0 ||
          strcmp(argv[i], "--add") == 0 ||
          strcmp(argv[i], "-xb") == 0 ||
          strcmp(argv[i], "-st") == 0) {
        err_msg("ignored when creating new system from scratch");
        PrintWarnOption("-bt/-ld/-hd/--bonded/-xb/-st");
        break;
      }
    }
  } //}}}
  //}}}

  SYSTEM S_orig;
  BOX *box = &S_orig.Box;
  if (opt.new) {
    InitSystem(&S_orig);
  } else {
    S_orig = ReadStructure(in, false);
  }
  COUNT *C_orig = &S_orig.Count;

  // find bead type to switch (the most numerous one; solvent, probably) //{{{
  opt.sw_type = NULL;
  if (!opt.add) {
    if (!(opt.sw_type = calloc(C_orig->BeadType, sizeof *opt.sw_type))) {
      ErrorAlloc("opt.sw_type");
    }
    // if -xb option not present, take the most numerous bead type
    if (!TypeOption(argc, argv, "-xb", 'b', true, opt.sw_type, S_orig)) {
      count = 0;
      int bt = 0;
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (S_orig.BeadType[i].Number > count) {
          count = S_orig.BeadType[i].Number;
          bt = i;
        }
      }
      opt.sw_type[bt] = true;
    }
  } //}}}

  // -bt <name(s)>/--bonded - specify what bead types to use //{{{
  opt.bt_use_orig = NULL;
  opt.bonded = false;
  if (!opt.new) {
    if (!(opt.bt_use_orig = calloc(C_orig->BeadType,
                                   sizeof *opt.bt_use_orig))) {
      ErrorAlloc("opt.bt_use_orig");
    }
    opt.bonded = BoolOption(argc, argv, "--bonded");
    TypeOption(argc, argv, "-bt", 'b', true, opt.bt_use_orig, S_orig);
  } //}}}

  // seed random number generator //{{{
  if (opt.seed != -1) {
    srand(opt.seed);
  } else {
    srand(time(0));
  } //}}}

  // read input coordinates //{{{
  if (in.coor.name[0] != '\0') {
    FILE *fr = OpenFile(in.coor.name, "r");
    int line_count = 0;
    for (int i = 1; i < commons.start; i++) { // from 1 as start=1 is the first
      if (!SkipTimestep(in, fr, &line_count)) {
        err_msg("couldn't skip");
        PrintError();
        exit(1);
      }
    }
    if (!ReadTimestep(in, fr, &S_orig, &line_count)) {
      err_msg("no coordinate data (starting step may be too high)");
      PrintErrorFile(in.coor.name, "\0", "\0");
      exit(1);
    }
    fclose(fr);
  } //}}}

  // read input FIELD file defining what to add //{{{
  SYSTEM S_add = ReadStructure(field, false);
  COUNT *C_add = &S_add.Count;
  C_add->BeadCoor = C_add->Bead;
  for (int i = 0; i < C_add->Bead; i++) {
    S_add.Bead[i].InTimestep = true;
    S_add.BeadCoor[i] = i;
  }
  if (opt.new) {
    S_orig.Box = S_add.Box;
  } //}}}

  // print original system (if there is any) //{{{
  if (commons.verbose && !opt.new) {
    fprintf(stdout, "\n==================================================");
    fprintf(stdout, "\nOriginal system");
    fprintf(stdout, "\n==================================================\n");
    VerboseOutput(S_orig);
    if (commons.start > 1) {
      fprintf(stdout, "\n   Using %d. timestep\n", commons.start);
    }
  } //}}}

  // new box if exists //{{{
  if (opt.box.Length.v[0] != -1) {
    if (!opt.new) {
      for (int dd = 0; dd < 3; dd++) {
        opt.box.Low.v[dd] += box->Low.v[dd] +
                             0.5 * (box->Length.v[dd] - opt.box.Length.v[dd]);
      }
    }
    opt.box.alpha = 90;
    opt.box.beta = 90;
    opt.box.gamma = 90;
    CalculateBoxData(&opt.box, 0);
    if (commons.verbose) {
      fprintf(stdout, "\n==================================================");
      printf("\nNew box");
      fprintf(stdout, "\n==================================================\n");
      PrintBox(opt.box);
    }
    *box = opt.box;
  } //}}}

  // move the beads (-off option) //{{{
  if (!opt.new) {
    if (!opt.real) { // transform offset to 'real' units if necessary
      for (int dd = 0; dd < 3; dd++) {
        opt.off->v[dd] *= S_orig.Box.Length.v[dd];
      }
    }
    for (int i = 0; i < C_orig->Bead; i++) {
      int id = S_orig.BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        S_orig.Bead[id].Position.v[dd] += opt.off->v[dd];
      }
    }
  } //}}}

  // minimize initial coordinates of added molecules //{{{
  for (int i = 0; i < C_add->Molecule; i++) {
    int type = S_add.Molecule[i].Type;
    double zero[3];
    if (opt.head) {
      int id0 = S_add.Molecule[i].Bead[0];
      for (int dd = 0; dd < 3; dd++) {
        zero[dd] = S_add.Bead[id0].Position.v[dd];
      }
    } else if (opt.tail) {
      int n = S_add.MoleculeType[S_add.Molecule[i].Type].nBeads;
      int id0 = S_add.Molecule[i].Bead[n-1];
      for (int dd = 0; dd < 3; dd++) {
        zero[dd] = S_add.Bead[id0].Position.v[dd];
      }
    } else {
      GeomCentre(S_add.MoleculeType[type].nBeads,
                 S_add.Molecule[i].Bead, S_add.Bead, zero);
    }
    for (int j = 0; j < S_add.MoleculeType[type].nBeads; j++) {
      int id = S_add.Molecule[i].Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        S_add.Bead[id].Position.v[dd] -= zero[dd];
      }
    }
  } //}}}

  // recalculate possible fractional constraints into true dimensions //{{{
  if (!opt.real) {
    for (int dd = 0; dd < 3; dd++) {
      for (int i = 0; i < 2; i++) {
        if (opt.axis[dd][i] != -1) {
          opt.axis[dd][i] *= box->Length.v[dd];
        }
      }
    }
  } //}}}

  // print what is to be added //{{{
  if (commons.verbose) {
    fprintf(stdout, "\n==================================================");
    fprintf(stdout, "\nBeads and molecules to add");
    fprintf(stdout, "\n==================================================\n");
    VerboseOutput(S_add);
  } //}}}

  // create output System //{{{
  SYSTEM S_out;
  COUNT *C_out = &S_out.Count;
  // if not switched, concatenate the new (i.e., original) and the added systems
  if (!opt.add) { // beads are to be switched, so transform the system
    // error - too few beads to switch //{{{
    // first, count number of beads that can be exchanged
    count = 0;
    for (int i = 0; i < C_orig->BeadType; i++) {
      if (opt.sw_type[i]) {
        count += S_orig.BeadType[i].Number;
      }
    }
    // second, the error?
    if (C_add->Bead > count) {
      err_msg("not enough beads to switch");
      PrintError();
      exit(1);
    } //}}}
    for (int i = 0; i < C_add->Bead; i++) {
      for (int j = 0; j < C_orig->BeadType; j++) {
        if (opt.sw_type[j] && S_orig.BeadType[j].InCoor > 0) {
          count = S_orig.BeadType[j].InCoor - 1;
          int id = S_orig.BeadType[j].Index[count];
          S_orig.Bead[id].InTimestep = false;
          S_orig.BeadType[j].InCoor--;
          break;
        }
      }
    }
    PruneSystem(&S_orig);
  }

  S_out = CopySystem(S_orig);
  ConcatenateSystems(&S_out, S_add, S_orig.Box, false); //}}}

  // define constrained box for adding beads (-cx/y/z and/or -hd options) //{{{
  opt.box = InitBox;
  for (int dd = 0; dd < 3; dd++) {
    opt.box.Length.v[dd] = S_out.Box.Length.v[dd];
  }
  // minimize box if -hd is used
  if (opt.hd) {
    // find minimum/maximum coordinates of beads for distance check //{{{
    double max[3] = {0, 0, 0}, min[3];
    for (int dd = 0; dd < 3; dd++) {
      min[dd] = S_orig.Box.Length.v[dd];
    }
    if (opt.bonded) { // use all bonded beads
      for (int i = 0; i < C_orig->BondedCoor; i++) {
        int id = S_orig.BondedCoor[i];
        BEAD *b = &S_orig.Bead[id];
        for (int dd = 0; dd < 3; dd++) {
          if (b->Position.v[dd] < min[dd]) {
            min[dd] = b->Position.v[dd];
          } else if (b->Position.v[dd] > max[dd]) {
            max[dd] = b->Position.v[dd];
          }
        }
      }
    } else { // use bead types specified by -bt
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (opt.bt_use_orig[i]) {
          for (int j = 0; j < S_orig.BeadType[i].Number; j++) {
            int id = S_orig.BeadType[i].Index[j];
            BEAD *b = &S_orig.Bead[id];
            if (b->InTimestep) {
              for (int dd = 0; dd < 3; dd++) {
                if (b->Position.v[dd] < min[dd]) {
                  min[dd] = b->Position.v[dd];
                }
                if (b->Position.v[dd] > max[dd]) {
                  max[dd] = b->Position.v[dd];
                }
              }
            }
          }
        }
      }
    } //}}}
    // the maximum/minimum possible coordinate of any added bead
    for (int dd = 0; dd < 3; dd++) {
      max[dd] += opt.hdist;
      min[dd] -= opt.hdist;
    }
    // define the box
    for (int dd = 0; dd < 3; dd++) {
      opt.box.Length.v[dd] = max[dd] - min[dd];
      opt.box.Low.v[dd] = min[dd];
    }
    CalculateBoxData(&opt.box, 0);
  }
  for (int dd = 0; dd < 3; dd++) {
    if (opt.axis[dd][0] != -1) {
      opt.box.Low.v[dd] = opt.axis[dd][0];
      opt.box.Length.v[dd] = opt.axis[dd][1] - opt.axis[dd][0];
    }
  }
  CalculateBoxData(&opt.box, 0);
  //}}}

  // what beads to check distance from for placing? //{{{
  int mode = 0; // no check
  if (!opt.new) {
    if (opt.bonded) { // all bonded beads
      mode = 1;
    } else { // possibly some speficied bead type(s)
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (opt.bt_use_orig[i]) { // yes, some specified bead type(s)
          mode = 2;
          break;
        }
      }
    }
  } //}}}

  // add monomeric beads //{{{
  for (int i = 0; i < C_add->Unbonded; i++) {
    vec3d random = RandomConstrainedCoor(S_orig, mode, S_out.Box.Length, opt);
    int id = C_orig->Bead + i;
    for (int dd = 0; dd < 3; dd++) {
      S_out.Bead[id].Position.v[dd] = random.v[dd];
    }
    // print number of placed beads?
    if (!commons.silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rMonomers placed: %d", i + 1);
    }
  } //}}}
  // print total number of placed beads? //{{{
  if (!commons.silent && C_add->Unbonded > 0) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                           \r");
    }
    fprintf(stdout, "\rMonomers placed: %d\n", C_add->Unbonded);
  } //}}}

  // add molecules //{{{
  for (int i = 0; i < C_add->Molecule; i++) {
    int mtype = S_out.Molecule[C_orig->Molecule+i].Type;
    double (*rot)[3];
    if (!(rot = calloc(S_out.MoleculeType[mtype].nBeads, sizeof *rot))) {
      ErrorAlloc("rot");
    }
    if (opt.no_rot) {
      for (int j = 0; j < S_out.MoleculeType[mtype].nBeads; j++) {
        int id_add = S_add.Molecule[i].Bead[j];
        for (int dd = 0; dd < 3; dd++) {
          rot[j][dd] = S_add.Bead[id_add].Position.v[dd];
        }
      }
    } else {
      Rotate(S_add, S_out.MoleculeType[mtype].nBeads,
             S_add.Molecule[i].Bead, opt.angle->v, rot);
    }
    vec3d random = RandomConstrainedCoor(S_orig, mode, S_out.Box.Length, opt);
    for (int j = 0; j < S_out.MoleculeType[mtype].nBeads; j++) {
      int id = S_out.Molecule[C_orig->Molecule+i].Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        S_out.Bead[id].Position.v[dd] = rot[j][dd] + random.v[dd];
      }
    }
    free(rot);
    // print number of placed molecules?
    if (!commons.silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rMolecules placed: %d", i + 1);
    }
  } //}}}
  // print total number of placed molecules? //{{{
  if (!commons.silent && C_add->Molecule > 0) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                           \r");
    }
    fprintf(stdout, "\rMolecules placed: %d\n", C_add->Molecule);
  } //}}}

  // TODO: the whole VtfSystem() and S_out2 needed for it doesn't seem necessry;
  //       who cares if some molecule is called differently because it contains
  //       data unsaveable to vtf? Actually, it might better as it shows the
  //       input(s) contain different molecules...
  SYSTEM S_out2;
  if (opt.fout.name[0] != '\0') {
    S_out2 = CopySystem(S_out);
    if (opt.fout.type == VCF_FILE ||
        opt.fout.type == VSF_FILE ||
        opt.fout.type == VTF_FILE) {
      VtfSystem(&S_out2);
    }
    PruneSystem(&S_out2);
  }
  if (fout.type == VCF_FILE ||
      fout.type == VSF_FILE ||
      fout.type == VTF_FILE) {
    VtfSystem(&S_out);
  }
  PruneSystem(&S_out);

  // print information about new system //{{{
  if (commons.verbose) {
    fprintf(stdout, "\n==================================================");
    fprintf(stdout, "\nNew system");
    fprintf(stdout, "\n==================================================\n");
    VerboseOutput(S_out);
  } //}}}

  // write data to output file(s) //{{{
  bool *write = malloc(sizeof *write * C_out->Bead);
  if (!write) {
    ErrorAlloc("write");
  }
  InitBoolArray(write, C_out->Bead, true); // save all beads
  WriteOutput(S_out, write, fout, false, -1, argc, argv);
  if (opt.fout.name[0] != '\0') {
    WriteOutput(S_out2, write, opt.fout, false, -1, argc, argv);
  } //}}}

  // free memory //{{{
  FreeSystem(&S_orig);
  FreeSystem(&S_add);
  FreeSystem(&S_out);
  if (opt.fout.name[0] != '\0') {
    FreeSystem(&S_out2);
  }
  if (!opt.new) {
    free(opt.bt_use_orig);
    free(opt.sw_type);
  }
  free(write);
  //}}}

  return 0;
}
