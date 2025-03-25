#include "../AnalysisTools.h"
#include <gsl/gsl_sf_legendre.h>  // GSL for spherical harmonics
#include <gsl/gsl_sf_coupling.h>  // GSL for Wigner 3-j symbols

#define LMAX 6  // Maximum order for q_l
#define CUTOFF 1.5  // Neighbour cutoff distance

// TODO: the -x option

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Write about BOOP!\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <in.agg> <output> [options]\n\n", cmd);

  fprintf(ptr, "   <input>    input coordinate file\n");
  fprintf(ptr, "   <in.agg>   input agg file\n");
  fprintf(ptr, "   <output>   file with BOOP data for the largest cluster\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined          specify that <input> contains joined "
          "coordinates\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool joined; // --joined
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

// Particle structure
typedef struct {
  double x, y, z;
} Particle;

// Compute the square distance between two particles //{{{
double dist2(Particle a, Particle b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;
  return dx*dx + dy*dy + dz*dz;
} //}}}

// Compute the spherical harmonics Y_{lm} using GSL //{{{
complex double Y_lm(int l, int m, double theta, double phi) {
  if (abs(m) > l) return 0.0;
  double P_lm = gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
  double factor = sqrt((2.0 * l + 1) / (4.0 * PI));
  return factor * P_lm * cexp(I * m * phi);
} //}}}

double compute_q_l(Particle *particles, int npart, int l) { //{{{
  complex double Q_lm_global[2 * l + 1]; // Global order parameter
  for (int m = -l; m <= l; m++) Q_lm_global[m + l] = 0.0;

  int total_neighbours = 0;
  for (int i = 0; i < npart; i++) {
    int count = 0;
    for (int j = 0; j < npart; j++) {
      if (i != j && dist2(particles[i], particles[j]) < Square(CUTOFF)) {
        count++;
      }
    }
    total_neighbours += count;
    // printf("Particle %d has %d neighbours\n", i, count); // Optional per-particle check
  }
  // printf("Average neighbours per particle: %.2f\n", total_neighbours / (double)npart);


  for (int i = 0; i < npart; i++) {
    complex double Q_lm[2 * l + 1]; // Local order parameter
    for (int m = -l; m <= l; m++) Q_lm[m + l] = 0.0;

    int count = 0;
    for (int j = 0; j < npart; j++) {
      if (i != j && dist2(particles[i], particles[j]) < Square(CUTOFF)) {
        count++;
        double dx = particles[j].x - particles[i].x;
        double dy = particles[j].y - particles[i].y;
        double dz = particles[j].z - particles[i].z;
        double r = sqrt(dx * dx + dy * dy + dz * dz);
        double theta = acos(fmax(-1.0, fmin(1.0, dz / r))); // Avoid NaN
        double phi = atan2(dy, dx);
        for (int m = -l; m <= l; m++) {
          Q_lm[m + l] += Y_lm(l, m, theta, phi);
        }
      }
    }
    if (count > 0) {
      for (int m = -l; m <= l; m++) Q_lm[m + l] /= count;
    }

    for (int m = -l; m <= l; m++) {
      Q_lm_global[m + l] += Q_lm[m + l]; // Sum over all particles
    }
  }

  // Final normalisation
  for (int m = -l; m <= l; m++) Q_lm_global[m + l] /= npart;

  double sum = 0.0;
  for (int m = -l; m <= l; m++) sum += cabs(Q_lm_global[m + l]) * cabs(Q_lm_global[m + l]);
  return sqrt(4.0 * PI / (2 * l + 1) * sum);
} //}}}

double compute_W_l(Particle *particles, int npart, int l) { //{{{
  complex double Q_lm_global[2 * l + 1];
  for (int m = -l; m <= l; m++) Q_lm_global[m + l] = 0.0;

  for (int i = 0; i < npart; i++) {
    complex double Q_lm[2 * l + 1];
    for (int m = -l; m <= l; m++) Q_lm[m + l] = 0.0;

    int count = 0;
    for (int j = 0; j < npart; j++) {
      if (i != j && dist2(particles[i], particles[j]) < CUTOFF * CUTOFF) {
        count++;
        double dx = particles[j].x - particles[i].x;
        double dy = particles[j].y - particles[i].y;
        double dz = particles[j].z - particles[i].z;
        double r = sqrt(dx * dx + dy * dy + dz * dz);
        double theta = acos(fmax(-1.0, fmin(1.0, dz / r))); // Avoid NaN
        double phi = atan2(dy, dx);
        for (int m = -l; m <= l; m++) {
          Q_lm[m + l] += Y_lm(l, m, theta, phi);
        }
      }
    }
    if (count > 0) {
      for (int m = -l; m <= l; m++) Q_lm[m + l] /= count;
    }

    for (int m = -l; m <= l; m++) {
      Q_lm_global[m + l] += Q_lm[m + l];
    }
  }

  for (int m = -l; m <= l; m++) Q_lm_global[m + l] /= npart;

  double sum = 0.0;
  for (int m1 = -l; m1 <= l; m1++) {
    for (int m2 = -l; m2 <= l; m2++) {
      for (int m3 = -l; m3 <= l; m3++) {
        double w3j = gsl_sf_coupling_3j(2 * l, 2 * l, 2 * l, 2 * m1, 2 * m2, 2 * m3);
        sum += w3j * creal(Q_lm_global[m1 + l] * Q_lm_global[m2 + l] * Q_lm_global[m3 + l]);
      }
    }
  }

  double q_l_value = compute_q_l(particles, npart, l);
  return (q_l_value > 1e-8) ? sum / pow(q_l_value, 3) : 0.0;
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 8, all = common + 1, count = 0,
      req_arg = 3;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
              "-st", "-e", "-sk", "-i", "--verbose", "--silent",
              "--help", "--version", "--joined");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}

  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  // test if <in.agg> ends with '.agg'
  int ext = 1;
  char extension[2][EXTENSION];
  s_strcpy(extension[0], ".agg", EXTENSION);
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(StripPath(argv[0]), true, common, option);
    exit(1);
  } //}}}

  // <output> - filename with data during simulation run
  char output[LINE];
  s_strcpy(output, argv[++count], LINE);

  // options before reading system data //{{{
  opt->c = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt->joined = true; // joined coordinates supplied, so no need to join
  } else {
    opt->joined = false; // molecules need to be joined
  } //}}}
  //}}}

  // print command to stdout
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // write initial stuff to output file //{{{
  PrintByline(output, argc, argv);
  FILE *out = OpenFile(output, "a");
  // print legend line to output file
  count = 1;
  fprintf(out, "# (%d) timestep", count++);
  fprintf(out, ", (%d) size", count++);
  fprintf(out, ", (%d) q_4", count++);
  fprintf(out, ", (%d) q_6", count++);
  fprintf(out, ", (%d) W_4", count++);
  fprintf(out, ", (%d) W_6", count++);
  putc('\n', out);
  fclose(out); //}}}

  // open input aggregate file and skip the first lines (Aggregate command & blank line) //{{{
  double distance = 1.6; // TODO: read from agg file
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO go for while(fgets()); treatment
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  AGGREGATE *Aggregate = NULL;
  InitAggregate(System, &Aggregate);

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  // FILE_TYPE fcoor;
  // s_strcpy(fcoor.name, "join.vtf", LINE);
  // fcoor.type = CoordinateFileType(fcoor.name);
  // if (fcoor.type == VTF_FILE) {
  //   WriteStructure(fcoor, System, -1, false, argc, argv);
  // }

  for (int i = 0; i < Count->BeadType; i++) {
    System.BeadType[i].Flag = true;
  }

  // main loop //{{{
  FILE *vcf = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0, // count lines in the vcf file
      line_count_agg = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_coor, opt->c.start, opt->c.silent);

    bool use = false;
    if (UseStep(opt->c, count_coor)) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(in, vcf, &System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      ReadAggregates(agg, input_agg, &System, Aggregate, &line_count_agg);
      if (!opt->joined) {
        WrapJoinCoordinates(&System, false, true);
        RemovePBCAggregates(distance, Aggregate, &System);
      }
      // bool *write = malloc(Count->Bead * sizeof *write);
      // InitBoolArray(write, Count->Bead, true);
      // WriteTimestep(fcoor, System, count_coor, write, argc, argv);
      // free(write);

      int largest[2] = {0, 0};
      for (int i = 0; i < Count->Aggregate; i++) {
        if (Aggregate[i].nMolecules > largest[0]) {
          largest[0] = Aggregate[i].nMolecules;
          largest[1] = i;
        }
      }

      Particle *particles = calloc(Aggregate[largest[1]].nBeads,
                                   sizeof *particles);
      for (int i = 0; i < Aggregate[largest[1]].nBeads; i++) {
        int id = Aggregate[largest[1]].Bead[i];
        particles[i].x = System.Bead[id].Position[0];
        particles[i].y = System.Bead[id].Position[1];
        particles[i].z = System.Bead[id].Position[2];
      }

      double q4 = compute_q_l(particles, Aggregate[largest[1]].nBeads, 4);
      double q6 = compute_q_l(particles, Aggregate[largest[1]].nBeads, 6);
      double W4 = compute_W_l(particles, Aggregate[largest[1]].nBeads, 4);
      double W6 = compute_W_l(particles, Aggregate[largest[1]].nBeads, 6);

      // print data to output file //{{{
      out = OpenFile(output, "a");
      fprintf(out, "%5d", count_used); // timestep
      fprintf(out, " %5d", Aggregate[largest[1]].nBeads);
      fprintf(out, " %lf", q4);
      fprintf(out, " %lf", q6);
      fprintf(out, " %lf", W4);
      fprintf(out, " %lf", W6);
      putc('\n', out);
      fclose(out); //}}}

      free(particles);
    }

    // exit the main loop if reached user-specied end timestep
    if (count_coor == opt->c.end) {
      break;
    }
  }
  fclose(vcf);
  fclose(agg);
  PrintLastStep(count_coor, count_used, opt->c.silent); //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  free(opt); //}}}

  return 0;
}
