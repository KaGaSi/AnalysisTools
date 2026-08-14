#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AddToSystem either creates a system from scratch or adds unbonded beads "
  "and/or molecules to an existing system. The new components are defined "
  "by a FIELD-like file or via -lib and -mol options and are placed either "
  "randomly or according to several possible constraints. By default, the "
  "new species are switched with existing beads, but using the --add flag, "
  "they can be appended to the system, increasing total bead count. Use -ntot "
  "to fill to a target total bead count with given molecules (requires -lib).",

  "Usage: AddToSystem <input> [<in.field>] <output> [options]",
  .args = 2, // minimum: <input> and <output>; <in.field> optional with -lib
  .all = 36, // number of valid lines in OptSpec (not counting last {nullptr})
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
  {"<input>", nullptr, "input coordinate file or '-' to generate new system",
    OPT_ARG},
  {"<in.field>", nullptr, "FIELD file with beads to add (unused with -lib)",
    OPT_ARG},
  {"<output>", nullptr, "output coordinate file", OPT_ARG},
  {"-o", "<filename>", "output extra structure file", OPT_EXTRA},
  {"-lib", "<dir>", "library directory: load molecules via -mol", OPT_EXTRA},
  {"-mol", "<mol> <%/n> ...", "molecule name(s) with percentage (e.g. 5%) "
    "or count; multiple pairs can be specified", OPT_EXTRA},
  {"-sys", "<input> [output]", "system file: assign names to existing molecule "
    "types; optionally output updated info", OPT_EXTRA},
  {"-ld", "<float>", "lowest distance from chosen beads (default: none)",
    OPT_EXTRA},
  {"-hd", "<float>", "highest distance from chosen beads (default: none)",
    OPT_EXTRA},
  {"-bt", "<name(s)>", "bead types for -hd/-ld (default: none)", OPT_EXTRA},
  {"--bonded", nullptr, "use bonded beads for -hd/-ld (overwrites -bt option)",
    OPT_EXTRA},
  {"-xb", "<bead type>", "what bead type to exchange", OPT_EXTRA},
  {"--add", nullptr, "add instead of exchanging beads (overwrites -xb)",
    OPT_EXTRA},
  {"--no-rotate", nullptr, "do not randomly rotate molecules", OPT_EXTRA},
  {"-a", "3x<angle>", "rotate molecules around <x>, <y>, <z> axes by "
    "given degrees (overrides --no-rotate)", OPT_EXTRA},
  {"-cx", "<lo> <hi> ...", "constrain x-coordinate (in fraction of output "
    "box); multiple pairs give a union of ranges", OPT_EXTRA},
  {"-cy", "<lo> <hi> ...", "constrain y-coordinate (in fraction of output "
    "box); multiple pairs give a union of ranges", OPT_EXTRA},
  {"-cz", "<lo> <hi> ...", "constrain z-coordinate (in fraction of output "
    "box); multiple pairs give a union of ranges", OPT_EXTRA},
  {"--no-cion", nullptr, "add only the named molecules, not the counterions they declare (place those yourself)", OPT_EXTRA},
  {"--tail", nullptr, "use molecule's last bead for constraint checks "
    "(default: molecule's geometric centre)", OPT_EXTRA},
  {"--head", nullptr, "use molecule's first bead for constraint checks "
    "(overrides --tail)", OPT_EXTRA},
  {"--real", nullptr, "use real coordinates for-cx/-cy/-cz/-off options",
    OPT_EXTRA},
  {"-b", "<x> <y> <z> [3x<angle>]", "new box dimensions (in real units), "
    "optionally followed by the alpha, beta, and gamma angles", OPT_EXTRA},
  {"-off", "3x<float>", "original system's offset "
    "(in fractions of the output box)", OPT_EXTRA},
  {"-s", "<int>", "seed for random number generator", OPT_EXTRA},
  {"-ntot", "<n> <type>", "fill to n total beads using <type> from library "
    "(requires -lib)", OPT_EXTRA},
  {"-ebt", "<int>", "number of extra bead types (output lammps data file only)",
    OPT_EXTRA},
  {nullptr}
}; //}}}

// maximum number of <lo> <hi> pairs a single -cx/-cy/-cz option can take
enum { MAX_RANGE = 10 };

// molecule specs from -mol option //{{{
typedef struct {
  char name[MOL_NAME];
  double value; // percentage (if is_frac=true), or integer count
  bool is_frac;
} ADD_SPEC; //}}}

// structure for options //{{{
struct OPT {
  bool ld, hd;                        // -ld/-hd
  double ldist, hdist;                //
  /*
   * -cx/-cy/-cz ranges: [range][0] holds the lower and [range][1] the upper
   * bound of each range, with the vector component picking the axis; n_axis
   * gives the number of ranges per axis (0 for an unconstrained one)
   */
  vec3d axis[MAX_RANGE][2],           // -cx/-cy/-cz
        frac[MAX_RANGE][2];           // the ranges as fractions of opt.box
  vec3i n_axis;                       // number of ranges per axis
  bool frac_axis;                     // are the ranges fractions already?
  vec3d angle,              // -a
        off;                // -off
  bool *bt_use_orig,        // -bt
       *sw_type,            // -xb
       new,                 // generate new system from scratch?
       real, add, no_rot,   // --real/--add/--no-rotate
       bonded, head, tail,  // --bonded/--head/--tail
       no_cion;             // --no-cion
  char lib_dir[LINE],       // -lib
       sys_in[LINE],        // -sys
       sys_out[LINE];       //
  BOX box;                  // constrained placement box (via -cx/-cy/-cz/-hd)
  int seed,                 // -s
      ebt;                  // -ebt
  FILE_TYPE fout;           // -o
  ADD_SPEC *lib_mol_add;    // -mol
  int n_lib_mol_add;        //
  int ntot;                 // -ntot
  char ntot_name[MOL_NAME]; // -ntot
}; //}}}

// generate random point uniformly distributed inside the box //{{{
/*
 * The point is drawn in fractional coordinates and transformed to Cartesian
 * ones, so it is uniform inside a triclinic cell as well; for an orthogonal
 * box, the transformation matrix is diagonal and this reduces to
 * s[dd] * OrthoLength[dd] + Low[dd].
 *
 * frac[dd] restricts the fractional coordinate to n[dd] sub-ranges of <0,1),
 * carving out a smaller region of the same shape (see the -cx/-cy/-cz
 * options). With several sub-ranges, one is drawn with a probability
 * proportional to its width, keeping the density over their union uniform.
 */
vec3d RandomCoordinate(BOX box, vec3i n, const vec3d frac[MAX_RANGE][2]) {
  // fractional coordinates, each in one of the requested <lo,hi) sub-ranges
  vec3d s;
  for (int dd = 0; dd < 3; dd++) {
    int r = 0;
    if (n.v[dd] > 1) { // pick a sub-range
      double total = 0;
      for (int i = 0; i < n.v[dd]; i++) {
        total += frac[i][1].v[dd] - frac[i][0].v[dd];
      }
      double pick = (double)(rand()) / ((double)(RAND_MAX) + 1) * total;
      for (r = 0; r < (n.v[dd] - 1); r++) {
        pick -= frac[r][1].v[dd] - frac[r][0].v[dd];
        if (pick < 0) {
          break;
        }
      }
    }
    double number = (double)(rand()) / ((double)(RAND_MAX) + 1);
    s.v[dd] = frac[r][0].v[dd] +
              number * (frac[r][1].v[dd] - frac[r][0].v[dd]);
  }
  // fractional -> Cartesian
  vec3d random;
  for (int dd = 0; dd < 3; dd++) {
    random.v[dd] = box.transform[dd][0] * s.v[0] +
                   box.transform[dd][1] * s.v[1] +
                   box.transform[dd][2] * s.v[2] +
                   box.Low.v[dd];
  }
  return random;
} //}}}

// wrap a point into the simulation box //{{{
/*
 * Bead positions are stored relative to Box.Low (the readers subtract it and
 * the writers add it back), so the cell spans <0,OrthoLength) for an
 * orthogonal box and <0,1) in fractional coordinates for a tilted one.
 */
vec3d WrapIntoBox(const vec3d coor, const BOX *box) {
  if (fabs(box->alpha - 90) < 1e-5 &&
      fabs(box->beta - 90) < 1e-5 &&
      fabs(box->gamma - 90) < 1e-5) {
    return RestorePBC(coor, box->OrthoLength);
  }
  double s[3];
  for (int dd = 0; dd < 3; dd++) {
    s[dd] = box->inverse[dd][0] * coor.v[0] +
            box->inverse[dd][1] * coor.v[1] +
            box->inverse[dd][2] * coor.v[2];
  }
  for (int dd = 0; dd < 3; dd++) {
    s[dd] -= floor(s[dd]);
  }
  vec3d out;
  for (int dd = 0; dd < 3; dd++) {
    out.v[dd] = box->transform[dd][0] * s[0] +
                box->transform[dd][1] * s[1] +
                box->transform[dd][2] * s[2];
  }
  return out;
} //}}}

// generate random point constrained by distance from other beads //{{{
// helper function calculating and updating minium coordinate
void GetMinDist(BEAD bead, vec3d random, const BOX *box, double *min_dist) {
  vec3d dist = DistancePBC(bead.Position, random, box);
  dist.v[0] = VectLength(dist);
  if (dist.v[0] < *min_dist) {
    *min_dist = dist.v[0];
  }
}
/*
 * What beads to use for distance check (int mode)
 *   0...no checks
 *   1...all bonded beads
 *   2...specified bead types,
 *
 * The placement box may stick out of the simulation box (the -hd option
 * defines it as a region around the existing beads), so the generated point
 * must be wrapped back into the simulation box.
 */
vec3d RandomConstrainedCoor(SYSTEM S_orig, int mode, const BOX *box, OPT opt) {
  if (mode == 0) { // no distance check
    vec3d random = RandomCoordinate(opt.box, opt.n_axis, opt.frac);
    return WrapIntoBox(random, box);
  }
  COUNT *C_orig = &S_orig.Count;
  double min_dist = 0;
  int tries = 0;
  vec3d random;
  const int max_tries = 10000000;
  do {
    tries++;
    if (tries > max_tries) {
      err_msg("could not place bead: constraints may be unsatisfiable");
      PrintError();
      exit(1);
    }
    random = RandomCoordinate(opt.box, opt.n_axis, opt.frac);
    min_dist = HUGE_VAL;
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
      err_msg("RandomConstrainedCoor(): mode must be 0 to 2");
      PrintError();
      exit(1);
    }
  } while ((opt.ld && opt.ldist >= min_dist) ||
           (opt.hd && opt.hdist <= min_dist));
  return WrapIntoBox(random, box);
} //}}}

// rotate randomly given collection of beads (e.g., a molecule) //{{{
void Rotate(SYSTEM System, int number, const int *list,
            vec3d rot_angle, vec3d *new) {
  // rotation angles around x-, y-, and z-axes
  double alpha, beta, gamma,
         rot[3][3];
  // specified by -a option...
  if (rot_angle.v[0] != 0 || rot_angle.v[1] != 0 || rot_angle.v[2] != 0) {
    gamma = rot_angle.v[0] / 180 * PI; // around x
    beta  = rot_angle.v[1] / 180 * PI; // around y
    alpha = rot_angle.v[2] / 180 * PI; // around z
    // ZYX (yaw=alpha/Z, pitch=beta/Y, roll=gamma/X)
    rot[0][0] = cos(alpha) * cos(beta);
    rot[0][1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
    rot[0][2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);

    rot[1][0] = sin(alpha) * cos(beta);
    rot[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
    rot[1][2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);

    rot[2][0] = -sin(beta);
    rot[2][1] = cos(beta) * sin(gamma);
    rot[2][2] = cos(beta) * cos(gamma);
  // ...or random
  } else {
    alpha = (double)(rand()) / ((double)(RAND_MAX) + 1) * 2 * PI;
    beta  = acos(1.0 - 2.0 * (double)(rand()) / ((double)(RAND_MAX) + 1));
    gamma = (double)(rand()) / ((double)(RAND_MAX) + 1) * 2 * PI;
    // ZYZ (Rz(alpha)*Ry(beta)*Rz(gamma)) to give properly random angles
    rot[0][0] = cos(alpha) * cos(beta) * cos(gamma) - sin(alpha) * sin(gamma);
    rot[1][0] = sin(alpha) * cos(beta) * cos(gamma) + cos(alpha) * sin(gamma);
    rot[2][0] = -sin(beta) * cos(gamma);

    rot[0][1] = -cos(alpha) * cos(beta) * sin(gamma) - sin(alpha) * cos(gamma);
    rot[1][1] = -sin(alpha) * cos(beta) * sin(gamma) + cos(alpha) * cos(gamma);
    rot[2][1] = sin(beta) * sin(gamma);

    rot[0][2] = cos(alpha) * sin(beta);
    rot[1][2] = sin(alpha) * sin(beta);
    rot[2][2] = cos(beta);
  }
  // generate the rotated coordinates
  for (int i = 0; i < number; i++) {
    vec3d *pos = &System.Bead[list[i]].Position;
    for (int dd = 0; dd < 3; dd++) {
      new[i].v[dd] = rot[dd][0] * pos->v[0] +
                     rot[dd][1] * pos->v[1] +
                     rot[dd][2] * pos->v[2];
    }
  }
} //}}}

int main(int argc, char *argv[]) {

  // command line arguments before reading the structure //{{{
  OptionCheck(argc, argv, false, HelpDesc, opts);
  OPT opt;
  int count = 0;

  // -lib option (determines if [in.field] is necessary)
  FileOption(argc, argv, "-lib", opt.lib_dir);

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
  // <in.field> - FIELD file with species to add (optional when -lib used) //{{{
  SYS_FILES field = InitSysFiles;
  if (opt.lib_dir[0] == '\0') {
    s_strcpy(field.stru.name, argv[++count], LINE);
    field.stru.type = StructureFileType(field.stru.name);
    if (field.stru.type != FIELD_FILE) {
      err_msg("<in.FIELD> file required when -lib is not used");
      PrintError();
      exit(1);
    }
  } //}}}
  // <output> - coordinate and structure output file //{{{
  FILE_TYPE fout = InitFile;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name); //}}}

  COMMON_OPT commons = CommonOptions(argc, argv, in);
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // -o - extra output structure file
  opt.fout.name[0] = '\0';
  if (FileOption(argc, argv, "-o", opt.fout.name)) {
    opt.fout.type = FileType(opt.fout.name);
  }
  // -sys <input> [output] //{{{
  opt.sys_in[0] = '\0';
  opt.sys_out[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-sys") == 0) {
      if ((i + 1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG, "missing file name", LINE);
        PrintErrorOption("-sys");
        exit(1);
      }
      s_strcpy(opt.sys_in, argv[i+1], LINE);
      if ((i + 2) < argc && argv[i+2][0] != '-') {
        s_strcpy(opt.sys_out, argv[i+2], LINE);
      }
      break;
    }
  }
  // When building from scratch, a single -sys arg is output-only
  if (opt.new && opt.sys_in[0] != '\0' && opt.sys_out[0] == '\0') {
    s_strcpy(opt.sys_out, opt.sys_in, LINE);
    opt.sys_in[0] = '\0';
  } //}}}
  // -mol <mol> <%/n> [<mol> <%/n> ...] //{{{
  opt.lib_mol_add = nullptr;
  opt.n_lib_mol_add = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-mol") != 0) {
      continue;
    }
    // consume all name/value pairs until the next flag or end of argv
    int j = i + 1;
    while ((j + 1) < argc && argv[j][0] != '-') {
      opt.lib_mol_add = s_realloc(opt.lib_mol_add, (opt.n_lib_mol_add + 1) *
                                  sizeof *opt.lib_mol_add);
      s_strcpy(opt.lib_mol_add[opt.n_lib_mol_add].name, argv[j], MOL_NAME);
      char *val = argv[j+1];
      int vlen = strlen(val);
      if (val[vlen-1] == '%') {
        char tmp[32] = {0};
        strncpy(tmp, val, vlen - 1);
        double frac;
        if (!IsPosRealNumber(tmp, &frac)) {
          err_msg("invalid fraction (must be <num>%)");
          PrintErrorOption("-mol");
          exit(1);
        }
        opt.lib_mol_add[opt.n_lib_mol_add].value = frac;
        opt.lib_mol_add[opt.n_lib_mol_add].is_frac = true;
      } else {
        long n;
        if (!IsWholeNumber(val, &n) || n <= 0) {
          err_msg("invalid count (use positive integer, optionally with%)");
          PrintErrorOption("-mol");
          exit(1);
        }
        opt.lib_mol_add[opt.n_lib_mol_add].value = n;
        opt.lib_mol_add[opt.n_lib_mol_add].is_frac = false;
      }
      opt.n_lib_mol_add++;
      j += 2;
    }
    break;
  } //}}}
  // -ntot <n> <type>: fill to target total bead count //{{{
  opt.ntot = 0;
  opt.ntot_name[0] = '\0';
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-ntot") == 0) {
      if ((i + 2) >= argc) {
        s_strcpy(ERROR_MSG, "requires <n> and <type>", LINE);
        PrintErrorOption("-ntot");
        exit(1);
      }
      long n;
      if (!IsWholeNumber(argv[i+1], &n) || n <= 0) {
        s_strcpy(ERROR_MSG, "<n> must be a positive integer", LINE);
        PrintErrorOption("-ntot");
        exit(1);
      }
      opt.ntot = (int)n;
      s_strcpy(opt.ntot_name, argv[i+2], MOL_NAME);
      break;
    }
  }
  if (opt.ntot > 0 && opt.lib_dir[0] == '\0') {
    s_strcpy(ERROR_MSG, "missing mandatory -lib option", LINE);
    PrintErrorOption("-ntot");
    exit(1);
  } //}}}

  // lowest and/or highest distance from specified beads //{{{
  opt.ld = false;
  opt.hd = false;
  if (!opt.new) { // only if not generating system from scratch
    opt.ld = OneNumberOption(argc, argv, "-ld", &opt.ldist, 'd');
    opt.hd = OneNumberOption(argc, argv, "-hd", &opt.hdist, 'd');
  }
  // errors for -ld/-hd options //{{{
  if ((opt.ld && opt.ldist < 0) || (opt.hd && opt.hdist <= 0)) {
    err_msg("highest/lowest distance must be positive real number");
    PrintErrorOption("-ld/-hd");
    exit(1);
  }
  if (opt.ld && opt.hd && opt.ldist >= opt.hdist) {
    err_msg("highest distance must be higher than lowest distance");
    PrintErrorOption("-ld/-hd");
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
    }
  } //}}}
  //}}}
  opt.real = BoolOption(argc, argv, "--real");
  // axes constraints (-cx/y/z options) //{{{
  for (int dd = 0; dd < 3; dd++) {
    opt.n_axis.v[dd] = 0;
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
    double val[2*MAX_RANGE];
    int n = 0;
    if (!NumbersOption(argc, argv, 2 * MAX_RANGE, str, &n, val, 'd')) {
      continue;
    }
    if (n % 2 != 0) {
      err_msg("<lo> <hi> pairs required, i.e., an even number of arguments");
      PrintErrorOption(str);
      exit(1);
    }
    opt.n_axis.v[dd] = n / 2;
    for (int i = 0; i < opt.n_axis.v[dd]; i++) {
      opt.axis[i][0].v[dd] = val[2*i];
      opt.axis[i][1].v[dd] = val[2*i+1];
      if (opt.axis[i][0].v[dd] == opt.axis[i][1].v[dd]) {
        err_msg("two different distance values required");
        PrintErrorOption(str);
        exit(1);
      } else if (opt.axis[i][0].v[dd] > opt.axis[i][1].v[dd]) {
        SwapDouble(&opt.axis[i][0].v[dd], &opt.axis[i][1].v[dd]);
      }
      if (!opt.real &&
          (opt.axis[i][0].v[dd] > 1 || opt.axis[i][1].v[dd] > 1)) {
        err_msg("unless --real is used, -cx/y/z must be between 0 and 1");
        PrintErrorOption(str);
        exit(1);
      }
    }
    // sort the ranges and merge the overlapping ones, so that every point of
    // the region is equally likely
    for (int i = 1; i < opt.n_axis.v[dd]; i++) {
      for (int j = i;
           j > 0 && opt.axis[j-1][0].v[dd] > opt.axis[j][0].v[dd]; j--) {
        SwapDouble(&opt.axis[j-1][0].v[dd], &opt.axis[j][0].v[dd]);
        SwapDouble(&opt.axis[j-1][1].v[dd], &opt.axis[j][1].v[dd]);
      }
    }
    int kept = 0;
    bool overlap = false;
    for (int i = 1; i < opt.n_axis.v[dd]; i++) {
      if (opt.axis[i][0].v[dd] <= opt.axis[kept][1].v[dd]) { // overlaps
        overlap = true;
        if (opt.axis[i][1].v[dd] > opt.axis[kept][1].v[dd]) {
          opt.axis[kept][1].v[dd] = opt.axis[i][1].v[dd];
        }
      } else {
        kept++;
        opt.axis[kept][0].v[dd] = opt.axis[i][0].v[dd];
        opt.axis[kept][1].v[dd] = opt.axis[i][1].v[dd];
      }
    }
    if (overlap) {
      snprintf(ERROR_MSG, LINE, "overlapping or touching ranges merged; %d of "
               "the %d given remain", kept + 1, opt.n_axis.v[dd]);
      PrintWarnOption(str);
    }
    opt.n_axis.v[dd] = kept + 1; //}}}
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
  InitDoubleArray(opt.angle.v, 3, 0);
  if (ThreeNumbersOption(argc, argv, "-a", opt.angle.v, 'd')) {
    opt.no_rot = false;
  } //}}}
  opt.head = BoolOption(argc, argv, "--head");
  opt.tail = BoolOption(argc, argv, "--tail");
  opt.no_cion = BoolOption(argc, argv, "--no-cion");
  // new box dimensions (-b option) //{{{
  /*
   * Three numbers define an orthogonal box; three more are the alpha, beta,
   * and gamma angles of a triclinic one.
   */
  BOX newbox = InitBox;
  double temp[6] = {0};
  int n_temp = 0;
  if (NumbersOption(argc, argv, 6, "-b", &n_temp, temp, 'd')) {
    if (n_temp != 3 && n_temp != 6) {
      err_msg("three or six numeric arguments required");
      PrintErrorOption("-b");
      Help(true, HelpDesc, opts);
      exit(1);
    }
    for (int dd = 0; dd < 3; dd++) {
      newbox.Length.v[dd] = temp[dd];
    }
    if (newbox.Length.x <= 0 || newbox.Length.y <= 0 || newbox.Length.z <= 0) {
      err_msg("three positive box dimensions required");
      PrintErrorOption("-b");
      Help(true, HelpDesc, opts);
      exit(1);
    }
    if (n_temp == 6) {
      newbox.alpha = temp[3];
      newbox.beta = temp[4];
      newbox.gamma = temp[5];
      if (newbox.alpha <= 0 || newbox.alpha >= 180 ||
          newbox.beta <= 0 || newbox.beta >= 180 ||
          newbox.gamma <= 0 || newbox.gamma >= 180) {
        err_msg("box angles must be between 0 and 180 degrees");
        PrintErrorOption("-b");
        Help(true, HelpDesc, opts);
        exit(1);
      }
    }
  } //}}}
  // -off option
  InitDoubleArray(opt.off.v, 3, 0);
  ThreeNumbersOption(argc, argv, "-off", opt.off.v, 'd');
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
        PrintWarnOption("-bt/-ld/-hd/--bonded/--add/-xb/-st");
        break;
      }
    }
  } //}}}
  // extra bead types for data output (-ebt option)
  opt.ebt = 0;
  OneNumberOption(argc, argv, "-ebt", &opt.ebt, 'i');
  //}}}

  SYSTEM S_orig;
  BOX *box = &S_orig.Box;
  if (opt.new) {
    InitSystem(&S_orig);
  } else {
    S_orig = ReadStructure(in, false);
  }
  COUNT *C_orig = &S_orig.Count;

  // apply -sys by assigning molecule type names from system_info file
  if (opt.sys_in[0] != '\0') {
    ReadSysInfo(opt.sys_in, &S_orig);
  }

  // find bead type to switch (the most numerous one; solvent, probably) //{{{
  opt.sw_type = nullptr;
  if (!opt.add) {
    if (!(opt.sw_type = calloc(C_orig->BeadType, sizeof *opt.sw_type))) {
      ErrorAlloc("opt.sw_type");
    }
    // if -xb option not present, take the most numerous bead type
    if (!TypeOption(argc, argv, "-xb", 'b', true, opt.sw_type, S_orig)) {
      int max_count = 0;
      int bt = 0;
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (S_orig.BeadType[i].Number > max_count) {
          max_count = S_orig.BeadType[i].Number;
          bt = i;
        }
      }
      opt.sw_type[bt] = true;
    }
  } //}}}

  // -bt <name(s)>/--bonded - specify what bead types to use //{{{
  // TypeOption for -bt is deferred until after RenameBeadTypesFromLibrary, so
  // that library bead type names can be used when -sys is present.
  opt.bt_use_orig = nullptr;
  opt.bonded = false;
  if (!opt.new) {
    if (!(opt.bt_use_orig = calloc(C_orig->BeadType,
                                   sizeof *opt.bt_use_orig))) {
      ErrorAlloc("opt.bt_use_orig");
    }
    opt.bonded = BoolOption(argc, argv, "--bonded");
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

  // build S_add: either from FIELD file or from library //{{{
  SYSTEM S_add;
  LIBRARY lib = {0};
  if (opt.lib_dir[0] != '\0') {
    lib = ReadLibrary(opt.lib_dir);
    if (opt.sys_in[0] == '\0' && !opt.new) {
      s_strcpy(ERROR_MSG, "-lib requires -sys", LINE);
      PrintError();
      exit(1);
    }
    RenameBeadTypesFromLibrary(&S_orig, &lib, opt.lib_dir);
    if (opt.bt_use_orig) {
      TypeOption(argc, argv, "-bt", 'b', true, opt.bt_use_orig, S_orig);
    }
    // -ntot //{{{
    // precompute -mol bead totals, add fill molecules first so that
    // unbonded fill beads precede bonded -mol molecules in lib.System.Bead[]
    if (opt.ntot > 0) {
      int mol_beads = 0;
      for (int i = 0; i < opt.n_lib_mol_add; i++) {
        LIB_MOL_INFO inf = LibraryMoleculeInfo(opt.lib_dir,
                                               opt.lib_mol_add[i].name);
        if (inf.n_beads_total <= 0) {
          FreeMolInfo(&inf);
          continue;
        }
        int nm;
        if (opt.lib_mol_add[i].is_frac) {
          nm = round(opt.lib_mol_add[i].value / 100.0 *
                     C_orig->Bead / inf.n_beads_total);
        } else {
          nm = (int)opt.lib_mol_add[i].value;
        }
        if (nm > 0) {
          mol_beads += nm * inf.n_beads_total;
        }
        FreeMolInfo(&inf);
      }
      int n_fill_beads = opt.ntot;
      if (opt.new) {
        n_fill_beads -= mol_beads;
      } else {
        n_fill_beads -= C_orig->Bead + mol_beads;
      }
      if (n_fill_beads <= 0) {
        int n = C_orig->Bead;
        if (opt.new) {
          n = 0;
        }
        if (snprintf(ERROR_MSG, LINE, "already at or above given target "
                     "(%s%d + %d >= %d%s)",
                     ErrYellow(), n, mol_beads, opt.ntot, ErrCyan()) < 0) {
          ErrorSnprintf();
        }
        PrintWarnOption("-ntot");
      } else {
        LIB_MOL_INFO fill_info = LibraryMoleculeInfo(opt.lib_dir,
                                                     opt.ntot_name);
        if (fill_info.n_beads_total <= 0) {
          if (snprintf(ERROR_MSG, LINE, "no library file for '%s%s%s'",
                       ErrYellow(), opt.ntot_name, ErrRed()) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption("-ntot");
          exit(1);
        }
        int n_mols_fill = n_fill_beads / fill_info.n_beads_total;
        if (n_mols_fill <= 0) {
          if (snprintf(ERROR_MSG, LINE, "fill count rounds to zero "
                       "(%s%d%s beads needed, %s%d%s per %s%s%s)",
                       ErrYellow(), n_fill_beads, ErrRed(),
                       ErrYellow(), fill_info.n_beads_total, ErrRed(),
                       ErrYellow(), opt.ntot_name, ErrRed()) < 0) {
            ErrorSnprintf();
          }
          PrintWarnOption("-ntot");
        } else {
          ReadLibraryMolecule(opt.lib_dir, opt.ntot_name, n_mols_fill,
                              !opt.no_cion, &lib);
        }
        FreeMolInfo(&fill_info);
      }
    } //}}}
    // total beads to use for fraction->count calculation
    int ntot = C_orig->Bead;
    for (int i = 0; i < opt.n_lib_mol_add; i++) {
      LIB_MOL_INFO info = LibraryMoleculeInfo(opt.lib_dir,
                                              opt.lib_mol_add[i].name);
      if (info.n_beads_total <= 0) {
        if (snprintf(ERROR_MSG, LINE, "no library file for %s%s%s",
                     ErrYellow(), opt.lib_mol_add[i].name, ErrRed()) < 0) {
          ErrorSnprintf();
        }
        PrintErrorOption("-lib");
        exit(1);
      }
      int n_mols;
      if (opt.lib_mol_add[i].is_frac) {
        n_mols = round(opt.lib_mol_add[i].value / 100.0 * ntot /
                       info.n_beads_total);
      } else {
        n_mols = opt.lib_mol_add[i].value;
      }
      if (n_mols > 0) {
        ReadLibraryMolecule(opt.lib_dir, opt.lib_mol_add[i].name, n_mols,
                            !opt.no_cion, &lib);
      }
      FreeMolInfo(&info);
    }
    FillSystemNonessentials(&lib.System, true);
    S_add = lib.System; // shallow copy - lib.System arrays now owned by S_add
    S_add.BeadCoor = s_realloc(S_add.BeadCoor,
                               S_add.Count.Bead * sizeof *S_add.BeadCoor);
    S_add.Count.BeadCoor = S_add.Count.Bead;
    for (int i = 0; i < S_add.Count.Bead; i++) {
      S_add.Bead[i].InTimestep = true;
      S_add.BeadCoor[i] = i;
    }
  } else {
    // no -lib: process -bt now (bead types keep their original names)
    if (opt.bt_use_orig) {
      TypeOption(argc, argv, "-bt", 'b', true, opt.bt_use_orig, S_orig);
    }
    S_add = ReadStructure(field, false);
    S_add.Count.BeadCoor = S_add.Count.Bead;
    for (int i = 0; i < S_add.Count.Bead; i++) {
      S_add.Bead[i].InTimestep = true;
      S_add.BeadCoor[i] = i;
    }
    if (opt.new) {
      S_orig.Box = S_add.Box;
    }
  }
  COUNT *C_add = &S_add.Count;
  //}}}

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

  // a box is required to place the new species into //{{{
  /*
   * When adding to an existing system, the box must come from the input files,
   * as -b is defined relative to it; from scratch, it comes from <in.field>
   * (or the library) or from -b.
   */
  if (!opt.new && box->Volume == -1) {
    err_msg("input file(s) contain no box dimensions");
    PrintError();
    exit(1);
  } //}}}

  // new box if exists //{{{
  if (newbox.Length.v[0] != -1) {
    // the angles are either the InitBox default of 90 or given via -b
    if (!CalculateBoxData(&newbox, 0)) {
      PrintErrorOption("-b");
      exit(1);
    }
    // centre the original box in the new one
    if (!opt.new) {
      for (int dd = 0; dd < 3; dd++) {
        newbox.Low.v[dd] += box->Low.v[dd] + 0.5 * (box->OrthoLength.v[dd] -
                                                    newbox.OrthoLength.v[dd]);
      }
    }
    if (commons.verbose) {
      fprintf(stdout, "\n==================================================");
      printf("\nNew box");
      fprintf(stdout, "\n==================================================\n");
      PrintBox(newbox);
    }
    *box = newbox;
  }
  if (box->Volume == -1) {
    err_msg("no box dimensions; specify them via -b or in <in.field>");
    PrintError();
    exit(1);
  } //}}}

  // move the beads (-off option) //{{{
  if (!opt.new) {
    if (!opt.real) { // transform offset to 'real' units if necessary
      // a fraction of each cell vector; the same as multiplying by
      // OrthoLength for an orthogonal box
      const BOX *b = &S_orig.Box;
      vec3d frac = opt.off;
      for (int dd = 0; dd < 3; dd++) {
        opt.off.v[dd] = b->transform[dd][0] * frac.v[0] +
                        b->transform[dd][1] * frac.v[1] +
                        b->transform[dd][2] * frac.v[2];
      }
    }
    for (int i = 0; i < C_orig->Bead; i++) {
      int id = S_orig.BeadCoor[i];
      for (int dd = 0; dd < 3; dd++) {
        S_orig.Bead[id].Position.v[dd] += opt.off.v[dd];
      }
    }
  } //}}}

  // minimize initial coordinates of added molecules //{{{
  for (int i = 0; i < C_add->Molecule; i++) {
    MOLECULE *mol_add = &S_add.Molecule[i];
    MOLECULETYPE *mtype_add = &S_add.MoleculeType[mol_add->Type];
    vec3d zero;
    // specify where is [0,0,0] coordinate
    if (opt.head) { // the first bead
      zero = S_add.Bead[mol_add->Bead[0]].Position;
    } else if (opt.tail) { // the last bead
      int n = mtype_add->nBeads;
      zero = S_add.Bead[mol_add->Bead[n-1]].Position;
    } else { // the molecule's geometric centre
      zero = GeomCentre(mtype_add->nBeads, mol_add->Bead, S_add.Bead);
    }
    for (int j = 0; j < mtype_add->nBeads; j++) {
      int id = mol_add->Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        S_add.Bead[id].Position.v[dd] -= zero.v[dd];
      }
    }
  } //}}}

  // recalculate possible fractional constraints into true dimensions //{{{
  /*
   * Bead positions are stored relative to Box.Low (the readers subtract it and
   * the writers add it back), making Box.Low only metadata for the output. The
   * two branches are therefore not symmetric on purpose: fractions of the box
   * are already relative, while --real takes the coordinates as they appear in
   * the input file, i.e., absolute, and must have Box.Low subtracted.
   *
   * A tilted cell has no Cartesian sub-box: a slab perpendicular to an axis
   * cannot be periodic unless the cell is orthogonal in that plane. The
   * fractions are therefore kept as fractions of the cell vectors, carving out
   * a smaller cell of the same shape. This does not apply to -hd, which always
   * defines an axis-aligned region around the existing beads.
   */
  bool tilted = fabs(box->alpha - 90) > 1e-5 ||
                fabs(box->beta - 90) > 1e-5 ||
                fabs(box->gamma - 90) > 1e-5;
  opt.frac_axis = tilted && !opt.hd;
  if (opt.frac_axis) { // the ranges are already fractions of the cell vectors
    for (int dd = 0; dd < 3; dd++) {
      if (opt.n_axis.v[dd] > 0 && opt.real) {
        err_msg("a tilted box has no Cartesian sub-box; drop --real and "
                "use fractions of the cell instead");
        PrintErrorOption("--real");
        exit(1);
      }
    }
  } else if (!opt.real) {
    for (int dd = 0; dd < 3; dd++) {
      for (int i = 0; i < opt.n_axis.v[dd]; i++) {
        opt.axis[i][0].v[dd] *= box->OrthoLength.v[dd];
        opt.axis[i][1].v[dd] *= box->OrthoLength.v[dd];
      }
    }
  } else {
    for (int dd = 0; dd < 3; dd++) {
      for (int i = 0; i < opt.n_axis.v[dd]; i++) {
        opt.axis[i][0].v[dd] -= box->Low.v[dd];
        opt.axis[i][1].v[dd] -= box->Low.v[dd];
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
        BEADTYPE *btype = &S_orig.BeadType[j];
        if (opt.sw_type[j] && btype->InCoor > 0) {
          count = btype->InCoor - 1;
          int id = btype->Index[count];
          S_orig.Bead[id].InTimestep = false;
          btype->InCoor--;
          break;
        }
      }
    }
    // compact BeadCoor: remove exchanged beads (InTimestep=false)
    int new_coor = 0;
    for (int i = 0; i < C_orig->BeadCoor; i++) {
      int id = S_orig.BeadCoor[i];
      if (S_orig.Bead[id].InTimestep) {
        S_orig.BeadCoor[new_coor] = id;
        new_coor++;
      }
    }
    C_orig->BeadCoor = new_coor;
    PruneSystem(&S_orig, nullptr);
  }

  S_out = CopySystem(S_orig);
  ConcatenateSystems(&S_out, S_add, S_orig.Box, false); //}}}

  // define constrained box for adding beads (-cx/y/z and/or -hd options) //{{{
  /*
   * The placement box copies the shape (i.e., the angles) of the system box,
   * as -cx/-cy/-cz only restrict the fractional coordinates inside it. The -hd
   * option, on the other hand, defines an axis-aligned region around the
   * existing beads, that is, an orthogonal placement box.
   */
  opt.box = InitBox;
  if (opt.hd) {
    for (int dd = 0; dd < 3; dd++) {
      opt.box.Length.v[dd] = S_out.Box.OrthoLength.v[dd];
    }
  } else {
    for (int dd = 0; dd < 3; dd++) {
      opt.box.Length.v[dd] = S_out.Box.Length.v[dd];
    }
    opt.box.alpha = S_out.Box.alpha;
    opt.box.beta = S_out.Box.beta;
    opt.box.gamma = S_out.Box.gamma;
  }
  // minimize box if -hd is used
  if (opt.hd) {
    // find minimum/maximum coordinates of beads for distance check //{{{
    vec3d max, min;
    for (int dd = 0; dd < 3; dd++) {
      min.v[dd] = HUGE_VAL;
      max.v[dd] = -HUGE_VAL;
    }
    if (opt.bonded) { // use all bonded beads
      for (int i = 0; i < C_orig->BondedCoor; i++) {
        int id = S_orig.BondedCoor[i];
        BEAD *b = &S_orig.Bead[id];
        for (int dd = 0; dd < 3; dd++) {
          if (b->Position.v[dd] < min.v[dd]) {
            min.v[dd] = b->Position.v[dd];
          }
          if (b->Position.v[dd] > max.v[dd]) {
            max.v[dd] = b->Position.v[dd];
          }
        }
      }
    } else { // use bead types specified by -bt
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (opt.bt_use_orig[i]) {
          BEADTYPE *bt = &S_orig.BeadType[i];
          for (int j = 0; j < bt->Number; j++) {
            int id = bt->Index[j];
            BEAD *b = &S_orig.Bead[id];
            if (b->InTimestep) {
              for (int dd = 0; dd < 3; dd++) {
                if (b->Position.v[dd] < min.v[dd]) {
                  min.v[dd] = b->Position.v[dd];
                }
                if (b->Position.v[dd] > max.v[dd]) {
                  max.v[dd] = b->Position.v[dd];
                }
              }
            }
          }
        }
      }
    } //}}}
    // no bead to measure the distance from means an undefined box
    if (min.v[0] > max.v[0]) {
      err_msg("no beads to measure the distance from");
      PrintErrorOption("-hd");
      exit(1);
    }
    // the maximum/minimum possible coordinate of any added bead
    for (int dd = 0; dd < 3; dd++) {
      max.v[dd] += opt.hdist;
      min.v[dd] -= opt.hdist;
    }
    /*
     * Define the box; opt.box.Low is where the placement region starts, not
     * the system box's offset (which the bead positions are relative to). The
     * region may stick out of the simulation box, with the protruding part
     * wrapped back inside when a point is generated; if it is longer than the
     * simulation box, however, it would wrap onto itself, oversampling the
     * overlap, so it is capped at the box length in that case. Only an
     * orthogonal box can be capped this way, as an axis-aligned region tiles
     * the space only if the cell is not tilted.
     */
    for (int dd = 0; dd < 3; dd++) {
      double length = max.v[dd] - min.v[dd];
      if (!tilted && length > S_out.Box.OrthoLength.v[dd]) {
        opt.box.Length.v[dd] = S_out.Box.OrthoLength.v[dd];
        opt.box.Low.v[dd] = 0;
      } else {
        opt.box.Length.v[dd] = length;
        opt.box.Low.v[dd] = min.v[dd];
      }
    }
    CalculateBoxData(&opt.box, 0);
  }
  CalculateBoxData(&opt.box, 0);
  // express the axis constraints as fractions of the placement box
  for (int dd = 0; dd < 3; dd++) {
    if (opt.n_axis.v[dd] == 0) { // unconstrained axis spans the whole box
      opt.n_axis.v[dd] = 1;
      opt.frac[0][0].v[dd] = 0;
      opt.frac[0][1].v[dd] = 1;
      continue;
    }
    for (int i = 0; i < opt.n_axis.v[dd]; i++) {
      for (int b = 0; b < 2; b++) {
        if (opt.frac_axis) {
          opt.frac[i][b].v[dd] = opt.axis[i][b].v[dd];
        } else {
          opt.frac[i][b].v[dd] = (opt.axis[i][b].v[dd] - opt.box.Low.v[dd]) /
                                 opt.box.OrthoLength.v[dd];
        }
      }
    }
  }
  //}}}

  // what beads to check distance from for placing? //{{{
  int mode = 0; // no check
  if (!opt.new) {
    if (opt.bonded) { // all bonded beads
      if (C_orig->BondedCoor == 0) {
        err_msg("no bonded beads in the system");
        PrintErrorOption("--bonded");
        exit(1);
      }
      mode = 1;
    } else { // possibly some specified bead type(s)
      for (int i = 0; i < C_orig->BeadType; i++) {
        if (opt.bt_use_orig[i]) { // yes, some specified bead type(s)
          mode = -1;
          if (S_orig.BeadType[i].InCoor > 0) {
            mode = 2;
            break;
          }
        }
      }
      if (mode == -1) {
        err_msg("no beads of specified type(s) present");
        PrintErrorOption("-bt");
        exit(1);
      }
    }
  } //}}}

  // add monomeric beads //{{{
  /*
   * Indexed through S_add.Unbonded[] rather than assuming the free beads are
   * the first C_add->Unbonded of the added block. They are not, as soon as
   * anything unbonded follows a molecule - a counterion after its parent, or
   * water after '-mol CTAC 10 water 200'. The prefix assumption left those
   * beads at the origin while overwriting the molecules' coordinates, which
   * the molecule loop below then quietly put back.
   */
  for (int i = 0; i < C_add->Unbonded; i++) {
    vec3d random = RandomConstrainedCoor(S_orig, mode, &S_out.Box, opt);
    int id = C_orig->Bead + S_add.Unbonded[i];
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
    MOLECULE *s_out_mol = &S_out.Molecule[C_orig->Molecule+i];
    MOLECULETYPE *s_out_mt = &S_out.MoleculeType[s_out_mol->Type];
    vec3d *rot;
    if (!(rot = calloc(s_out_mt->nBeads, sizeof *rot))) {
      ErrorAlloc("rot");
    }
    MOLECULE *s_add_mol = &S_add.Molecule[i];
    if (opt.no_rot) {
      for (int j = 0; j < s_out_mt->nBeads; j++) {
        int id_add = s_add_mol->Bead[j];
        for (int dd = 0; dd < 3; dd++) {
          rot[j].v[dd] = S_add.Bead[id_add].Position.v[dd];
        }
      }
    } else {
      Rotate(S_add, s_out_mt->nBeads, s_add_mol->Bead, opt.angle, rot);
    }
    vec3d random = RandomConstrainedCoor(S_orig, mode, &S_out.Box, opt);
    for (int j = 0; j < s_out_mt->nBeads; j++) {
      int id = s_out_mol->Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        S_out.Bead[id].Position.v[dd] = rot[j].v[dd] + random.v[dd];
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

  // S_out2: separate copy for the secondary output (-o), which may need
  // VtfSystem() stripping independent of the primary output format
  SYSTEM S_out2;
  if (opt.fout.name[0] != '\0') {
    S_out2 = CopySystem(S_out);
    if (opt.fout.type == VCF_FILE ||
        opt.fout.type == VSF_FILE ||
        opt.fout.type == VTF_FILE) {
      VtfSystem(&S_out2);
    }
    PruneSystem(&S_out2, nullptr);
  }
  if (fout.type == VCF_FILE ||
      fout.type == VSF_FILE ||
      fout.type == VTF_FILE) {
    VtfSystem(&S_out);
  }
  PruneSystem(&S_out, nullptr);

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
  if (fout.type == LDATA_FILE && opt.ebt > 0) {
    for (int i = 0; i < opt.ebt; i++) {
      NewBeadType(&S_out.BeadType, &S_out.Count.BeadType, "extra", 0, 1, 1);
    }
  }
  WriteOutput(S_out, write, fout, false, -1, argc, argv);
  if (opt.fout.name[0] != '\0') {
    if (opt.fout.type == LDATA_FILE && opt.ebt > 0) {
      for (int i = 0; i < opt.ebt; i++) {
        NewBeadType(&S_out2.BeadType, &S_out2.Count.BeadType, "extra", 0, 1, 1);
      }
    }
    WriteOutput(S_out2, write, opt.fout, false, -1, argc, argv);
    // the interactions block lives in the library, not in the system, so
    // WriteOutput() cannot produce it (Info does the same for its -o)
    if (opt.lib_dir[0] != '\0' && opt.fout.type == FIELD_FILE) {
      AppendFieldInteractions(opt.fout.name, &S_out2, &lib);
    }
  } //}}}

  // write system_info if -sysout specified //{{{
  if (opt.sys_out[0] != '\0') {
    WriteSysInfo(opt.sys_out, &S_out);
  } //}}}

  // free memory //{{{
  FreeSystem(&S_orig);
  if (opt.lib_dir[0] != '\0') {
    // S_add owns lib.System's arrays, so free everything but that
    free(lib.inter);
    free(lib.bond_id);
    free(lib.angle_id);
  }
  FreeSystem(&S_add);
  FreeSystem(&S_out);
  if (opt.fout.name[0] != '\0') {
    FreeSystem(&S_out2);
  }
  if (!opt.new) {
    free(opt.bt_use_orig);
    if (!opt.add) {
      free(opt.sw_type);
    }
  }
  free(opt.lib_mol_add);
  free(write);
  //}}}

  return 0;
}
