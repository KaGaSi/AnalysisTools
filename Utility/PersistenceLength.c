#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
PersistenceLength calculates correlation of bond vectors and its standard \
deviation. To get the peristence length, the data must be fitted via, \
typically, an exponential function. Note the utility expects a linear chain \
with ordered beads ids (e.g., for a 4-bead chain, the order must be 1-2-3-4, \
leading to ordered bonds 1-2, 2-3, 3-4; connectivity like 1-4-2-3 with \
bonds 1-4, 2-4, 2-3 could give unexpected results).\n\n");
  }
  fprintf(ptr, "Usage: %s <input> <output> [options]\n\n", cmd);

  fprintf(ptr, "<input>             input coordinate file\n");
  fprintf(ptr, "<output>            output file with the persistence length\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -m <name(s)>      molecule types to calculate bond lengths "
          "for (if not present, use all molecule types)\n");
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
  OptionCheck(argc, argv, req_arg, common, all, true, option,
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
  opt->c = CommonOptions(argc, argv, in);
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
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->mt, System)) {
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
  double **cos_phi = calloc(Count->MoleculeType, sizeof *cos_phi),
         **cos2_phi = calloc(Count->MoleculeType, sizeof *cos2_phi),
         (*avg_bond)[2] = calloc(Count->MoleculeType, sizeof avg_bond[2]);
  long int **count_stuff = calloc(Count->MoleculeType, sizeof *count_stuff);
  for (int i = 0; i < Count->MoleculeType; i++) {
    cos_phi[i] = calloc(max_bonds, sizeof *cos_phi[i]);
    cos2_phi[i] = calloc(max_bonds, sizeof *cos2_phi[i]);
    count_stuff[i] = calloc(max_bonds, sizeof *count_stuff[i]);
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
              for (int lag = 0; lag < (mt_i->nBonds - 1); lag++) {
                for (int j = 0; j < (mt_i->nBonds - lag - 1); j++) {
                  int id1 = mol_i->Bead[mt_i->Bond[j][0]],
                      id2 = mol_i->Bead[mt_i->Bond[j][1]];
                  BEAD *bj1 = &System.Bead[id1],
                       *bj2 = &System.Bead[id2];

                  int idk[2] = {mol_i->Bead[mt_i->Bond[j+lag][0]],
                                mol_i->Bead[mt_i->Bond[j+lag][1]]};
                  BEAD *bk1 = &System.Bead[idk[0]],
                       *bk2 = &System.Bead[idk[1]];

                  double u[3]; // the first bond vector
                  double v[3]; // the second bond vector
                  for (int dd = 0; dd < 3; dd++) {
                    u[dd] = bj1->Position[dd] - bj2->Position[dd];
                    v[dd] = bk1->Position[dd] - bk2->Position[dd];
                  }
                  // bond lengths
                  double size[2];
                  size[0] = VectLength(u);
                  size[1] = VectLength(v);
                  avg_bond[mol_i->Type][0] += size[0];
                  avg_bond[mol_i->Type][1]++;
                  double res = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
                  res /= (size[0] * size[1]);
                  // sum of cos(\phi)
                  cos_phi[mol_i->Type][lag] += res;
                  // sum of cos(\phi)^2 for error
                  cos2_phi[mol_i->Type][lag] += Square(res);
                  // count values - easier than figuring out their number
                  count_stuff[mol_i->Type][lag]++;
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
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // average all arrays //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (opt->mt[i]) {
      MOLECULETYPE *mt = &System.MoleculeType[i];
      for (int lag = 0; lag < (mt->nBonds - 1); lag++) {
        cos_phi[i][lag] /= count_stuff[i][lag];
        cos2_phi[i][lag] /= count_stuff[i][lag];
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
  // start at 1 as lag=0 is always 1, end at max_bonds-1 as the last has nothing
  for (int lag = 1; lag < (max_bonds - 1); lag++) {
    fprintf(fw, "%d", lag);
    for (int j = 0; j < Count->MoleculeType; j++) {
      if (opt->mt[j] && lag < System.MoleculeType[j].nBonds) {
        fprintf(fw, " %lf", cos_phi[j][lag]);
        double err = sqrt((cos2_phi[j][lag] - Square(cos_phi[j][lag])) /
                          count_stuff[j][lag]);
        fprintf(fw, " %lf", err);
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
    free(cos_phi[i]);
    free(cos2_phi[i]);
    free(count_stuff[i]);
  }
  free(cos_phi);
  free(cos2_phi);
  free(count_stuff);
  free(avg_bond);
  FreeSystem(&System);
  free(opt);
  //}}}

  return 0;
}
