#include "../AnalysisTools.h"

void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
PersistenceLength calculates persistence length for a whole chain or \
a part of it (SOME OPTION). Results are printed as \
<r_i.r_i,i+j>/<cos\\phi_i,i+j> vs. j. Probably...\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file with the persistence length\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)");
  fprintf(ptr, "      --joined    specify that <input> contains joined "
          "coordinates\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join, // --joined
       *mt;  // -m
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// calculate index in an 1D array that simulates a 2D one /{{{
int id2D(int i1, int i2, int size[2]) {
  return (i2 * size[0] +
          i1);
} //}}}
// calculate index in an 1D array that simulates a 3D one /{{{
int id3D(int i1, int i2, int i3, int size[3]) {
  return (i3 * size[0] * size[1] +
          i2 * size[0] +
          i1);
} //}}}

int main(int argc, char *argv[]) {

  int common = 8, all = common + 2, count = 0,
      req_arg = 2;
  char option[all][OPT_LENGTH];
  OptionCheck2(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--joined", "-m");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  snprintf(in.coor.name, LINE, "%s", argv[++count]);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output> - file name with persistence lengths
  char fout[LINE] = "";
  snprintf(fout, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, LINE, in);
  // --joined option
  if (BoolOption(argc, argv, "--joined")) {
    opt->join = false; // joined coordinates supplied, so no need to join
  } else {
    opt->join = true; // molecules need to be joined
  } //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // '-m <name(s)>' option
  opt->mt = calloc(System.Count.MoleculeType, sizeof *opt->mt);
  if (!MoleculeTypeOption(argc, argv, "-m", true, opt->mt, System)) {
    InitBoolArray(opt->mt, Count->MoleculeType, true);
  }

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // maximum number of bonds //{{{
  int max_bonds = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i] && System.MoleculeType[i].nBonds > max_bonds) {
      max_bonds = System.MoleculeType[i].nBonds;
    }
  } //}}}

  // arrays for sums //{{{
  int arr[3];
  arr[0] = Count->MoleculeType;
  arr[1] = max_bonds;
  arr[2] = max_bonds;
  int arr_size = arr[0] * arr[1] * arr[2];
  double *dot = calloc(arr_size, sizeof *dot),
         *cos_phi = calloc(arr_size, sizeof *cos_phi), // TODO: useless?
         *bond = calloc(arr_size, sizeof *bond);
  long int *count_stuff = calloc(arr_size, sizeof *count_stuff);
  int arr_avg_bond[2];
  arr_avg_bond[0] = Count->MoleculeType;
  arr_avg_bond[1] = 2; // [0]...bond length, [1]...number of mols
  int arr_avg_bond_size = arr_avg_bond[0] * arr_avg_bond[1];
  double *avg_bond = calloc(arr_avg_bond_size, sizeof *avg_bond);
  //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
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
      WrapJoinCoordinates(&System, true, opt->join);
      // go through all molecules //{{{
      for (int i = 0; i < Count->MoleculeType; i++) {
        if (opt->mt[i]) {
          for (int mm = 0; mm < System.MoleculeType[i].Number; mm++) {
            int mol = System.MoleculeType[i].Index[mm];
            MOLECULE *mol_i = &System.Molecule[mol];
            MOLECULETYPE *mt_i = &System.MoleculeType[mol_i->Type];
            // use only specified molecule types
            if (mol_i->InTimestep) {
              for (int j = 0; j < mt_i->nBonds; j++) {
                int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                    id2 = mol_i->Bead[mt_i->Bond[j][1]];
                BEAD *bj1 = &System.Bead[id1],
                     *bj2 = &System.Bead[id2];
                double u[3]; // the first bond vector
                for (int dd = 0; dd < 3; dd++) {
                  u[dd] = bj1->Position[dd] - bj2->Position[dd];
                }
                // bond lengths
                double size[2];
                size[0] = VECTORLENGTH(u);
                int id = id2D(mol_i->Type, 0, arr_avg_bond);
                avg_bond[id] += size[0];
                id = id2D(mol_i->Type, 1, arr_avg_bond);
                avg_bond[id]++;
                for (int k = (j + 1); k < mt_i->nBonds; k++) {
                  int idk[2] = {mol_i->Bead[mt_i->Bond[k][0]],
                                mol_i->Bead[mt_i->Bond[k][1]]};
                  BEAD *bk1 = &System.Bead[idk[0]],
                       *bk2 = &System.Bead[idk[1]];
                  double v[3]; // the second bond vector
                  for (int dd = 0; dd < 3; dd++) {
                    v[dd] = bk1->Position[dd] - bk2->Position[dd];
                  }
                  // bond lengths
                  size[1] = VECTORLENGTH(v);
                  // size[1] = sqrt(SQR(v[0]) + SQR(v[1]) * SQR(v[2]));
                  double scalar = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
                  int id = id3D(mol_i->Type, j, k, arr);
                  // sum of dot products
                  dot[id] += scalar;
                  // sum of cos(\phi)
                  cos_phi[id] += scalar / (size[0] * size[1]);
                  // sum of bond lengths
                  bond[id] += size[0] * size[1];
                  // count values - easier than figuring out their number
                  count_stuff[id]++;
                }
              }
            }
          }
        }
      } //}}}
      //}}}
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(fr);
  // print last step?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d (used %d)\n", count_coor, count_used);
  } //}}}

  // average all arrays //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < (mt->nBonds - 1); j++) {
        for (int k = (j + 1); k < mt->nBonds; k++) {
          int id = id3D(i, j, k, arr);
          cos_phi[id] /= count_stuff[id];
          dot[id] /= count_stuff[id];
          bond[id] /= count_stuff[id];
        }
      }
    }
  } //}}}

  // recalculate arrays to be dependent on just distance between bonds //{{{
  int arr_id_dist[2];
  arr_id_dist[0] = Count->MoleculeType;
  arr_id_dist[1] = max_bonds;
  int arr_id_dist_size = arr_id_dist[0] * arr_id_dist[1];
  double *cos_phi_dist = calloc(arr_id_dist_size, sizeof *cos_phi_dist),
         *dot_dist = calloc(arr_id_dist_size, sizeof *dot_dist),
         *bond_dist = calloc(arr_id_dist_size, sizeof *bond_dist);
  int *count_stuff_dist = calloc(arr_id_dist_size, sizeof *count_stuff_dist);
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < (mt->nBonds - 1); j++) {
        count = 0;
        for (int k = (j + 1); k < mt->nBonds; k++) {
          int id2 = id2D(i, k - j - 1, arr_id_dist),
              id3 = id3D(i, j, k, arr);
          cos_phi_dist[id2] += cos_phi[id3];
          dot_dist[id2] += dot[id3];
          bond_dist[id2] += bond[id3];
          count_stuff_dist[id2]++;
        }
      }
    }
  }
  // average the arrays
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < (mt->nBonds - 1); j++) {
        int id2 = id2D(i, j, arr_id_dist);
        cos_phi_dist[id2] /= count_stuff_dist[id2];
        dot_dist[id2] /= count_stuff_dist[id2];
        bond_dist[id2] /= count_stuff_dist[id2];
      }
    }
  } //}}}

  // write data to ouptut file //{{{
  PrintByline(fout, argc, argv);
  // print first line of output file - molecule names and beadtype trios //{{{
  FILE *fw = OpenFile(fout, "a");
  fprintf(fw, "# (1) distance between bonds;");
  fprintf(fw, " molecules:");
  count = 1;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    if (opt->mt[i]) {
      count++;
      fprintf(fw, " (%d) %s molecule:", count, mt_i->Name);
    }
  }
  putc('\n', fw); //}}}
  // write data //{{{
  for (int i = 0; i < arr_id_dist[1]; i++) {
    fprintf(fw, "%d", i + 1);
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt->mt[j]) {
        int id = id2D(j, i, arr_id_dist);
        fprintf(fw, " %lf %lf", cos_phi_dist[id], dot_dist[id] / bond_dist[id]);
      }
    }
    putc('\n', fw);
  } //}}}
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      int id1 = id2D(i, 0, arr_avg_bond);
      int id2 = id2D(i, 1, arr_avg_bond);
      fprintf(fw, "# average bond length: %lf\n", avg_bond[id1] / avg_bond[id2]);
    }
  }
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystem(&System);
  free(opt->mt);
  free(opt);
  free(cos_phi);
  free(cos_phi_dist);
  free(dot);
  free(dot_dist);
  free(count_stuff);
  free(count_stuff_dist);
  free(bond);
  free(bond_dist);
  free(avg_bond);
  //}}}

  return 0;
}
