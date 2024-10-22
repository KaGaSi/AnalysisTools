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
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <output> - file name with persistence lengths
  char fout[LINE] = "";
  s_strcpy(fout, argv[++count], LINE);

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
  double ***cos_phi = calloc(Count->MoleculeType, sizeof *cos_phi), // TODO: useless?
         (*avg_bond)[2] = calloc(Count->MoleculeType, sizeof avg_bond[2]);
  long int ***count_stuff = calloc(Count->MoleculeType, sizeof *count_stuff);
  for (int i = 0; i < Count->MoleculeType; i++) {
    cos_phi[i] = calloc(max_bonds, sizeof *cos_phi);
    count_stuff[i] = calloc(max_bonds, sizeof *count_stuff);
    for (int j = 0; j < max_bonds; j++) {
      cos_phi[i][j] = calloc(max_bonds, sizeof *cos_phi);
      count_stuff[i][j] = calloc(max_bonds, sizeof *count_stuff);
    }
  } //}}}

  // main loop //{{{
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);
    // use every skip-th timestep between start and end
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
                avg_bond[mol_i->Type][0] += size[0];
                avg_bond[mol_i->Type][1]++;
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
                  double scalar = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
                  // sum of cos(\phi)
                  cos_phi[mol_i->Type][j][k] += scalar / (size[0] * size[1]);
                  // count values - easier than figuring out their number
                  count_stuff[mol_i->Type][j][k]++;
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
          cos_phi[i][j][k] /= count_stuff[i][j][k];
        }
      }
    }
  } //}}}

  // recalculate arrays to be dependent on just distance between bonds //{{{
  double **cos_phi_dist = calloc(Count->MoleculeType, sizeof *cos_phi_dist);
  int **count_stuff_dist = calloc(Count->MoleculeType,
                                  sizeof *count_stuff_dist);
  for (int i = 0; i < Count->MoleculeType; i++) {
    cos_phi_dist[i] = calloc(max_bonds, sizeof *cos_phi_dist);
    count_stuff_dist[i] = calloc(max_bonds, sizeof *count_stuff_dist);
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < (mt->nBonds - 1); j++) {
        count = 0;
        for (int k = (j + 1); k < mt->nBonds; k++) {
          cos_phi_dist[i][k-j-1] += cos_phi[i][j][k];
          count_stuff_dist[i][k-j-1]++;
        }
      }
    }
  }
  // average the arrays
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int j = 0; j < (mt->nBonds - 1); j++) {
        cos_phi_dist[i][j] /= count_stuff_dist[i][j];
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
      fprintf(fw, " (%d) %s", count, mt_i->Name);
    }
  }
  putc('\n', fw); //}}}
  // write data //{{{
  for (int i = 0; i < max_bonds; i++) {
    fprintf(fw, "%d", i + 1);
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt->mt[j] && i < System.MoleculeType[j].nBonds) {
        fprintf(fw, " %lf", cos_phi_dist[j][i]);
      }
    }
    putc('\n', fw);
  } //}}}
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      fprintf(fw, "# average bond length: %lf\n",
              avg_bond[i][0] / avg_bond[i][1]);
    }
  }
  fclose(fw); //}}}

  // free memory - to make valgrind happy //{{{
  free(opt->mt);
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < max_bonds; j++) {
      free(cos_phi[i][j]);
      free(count_stuff[i][j]);
    }
    free(cos_phi[i]);
    free(count_stuff[i]);
    free(cos_phi_dist[i]);
    free(count_stuff_dist[i]);
  }
  free(cos_phi);
  free(count_stuff);
  free(cos_phi_dist);
  free(count_stuff_dist);
  free(avg_bond);
  FreeSystem(&System);
  free(opt);
  //}}}

  return 0;
}
