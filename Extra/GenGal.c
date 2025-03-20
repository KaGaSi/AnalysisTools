#include "../AnalysisTools.h"

int length = 10000; // polysaccharide length
double alpha = 0.4; // dissociation degree
PARAMS bond = {.a = 30, .b = 0.6};
PARAMS angle = {.a = 10, .b = 160};
double box = 50;
char name[LINE] = "0_4G500";

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "Create 1000 Gal_500 molecule; hardcoded stuff...\n\n");
  }

  fprintf(ptr, "Usage: %s <output> [options]\n\n", cmd);

  fprintf(ptr, "<output>            output coordinate file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -o <struct>       extra output file printed after the last "
          "molecule is generated\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  FILE_TYPE fw; // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  srand(time(0));

  // define options
  int common = 4, all = common + 1, count = 0,
      req_arg = 1;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--joined", "--all", "-d", "-m", "-w");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();

  // <output> - output coordinate file
  FILE_TYPE fw_coor;
  snprintf(fw_coor.name, LINE, "%s", argv[++count]);
  fw_coor.type = CoordinateFileType(fw_coor.name);

  SYS_FILES in = InitSysFiles;
  opt->c = CommonOptions(argc, argv, in);
  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // output file (-o option) //{{{
  opt->fw = InitFile;
  if (FileOption(argc, argv, "-o", opt->fw.name)) {
    opt->fw.type = FileType(opt->fw.name);
  } //}}}

  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  // fill COUNT struct //{{{
  Count->Bead = length + length * alpha; // G,GH + C
  Count->BeadCoor = length + length * alpha;
  Count->Bonded = length; // G,GH
  Count->Unbonded = length * alpha; // C
  Count->Molecule = 1;
  Count->HighestResid = 0;
  Count->BondType = 1;
  Count->AngleType = 1; //}}}
  // box size
  System.Box.Length[0] = box;
  System.Box.Length[1] = box;
  System.Box.Length[2] = box;
  // allocate necessary arrays
  System.Bead = realloc(System.Bead, Count->Bead * sizeof *System.Bead);
  System.BeadCoor = realloc(System.BeadCoor,
                            Count->BeadCoor * sizeof *System.BeadCoor);
  // create bead types //{{{
  // dissociated polysaccharide bead
  NewBeadType(&System.BeadType, &Count->BeadType, "G", -1.038200, 2.6, RADIUS);
  System.BeadType[0].Number = length * alpha;
  // non-dissociated polysaccharide bead
  NewBeadType(&System.BeadType, &Count->BeadType, "GH", 0, 2.6, RADIUS);
  System.BeadType[1].Number = length - System.BeadType[0].Number;
  // cation counterion
  NewBeadType(&System.BeadType, &Count->BeadType, "C", 1.038200, 1, RADIUS);
  System.BeadType[2].Number = length * alpha; //}}}
  // create bond and angle types
  System.BondType[0] = bond;
  System.AngleType[0] = angle;
  // create new molecule type, fill bonds and angles //{{{
  NewMolType(&System.MoleculeType, &Count->MoleculeType, name, length,
             length - 1, length - 2, 0, 0);
  MOLECULETYPE *mt0 = &System.MoleculeType[0];
  mt0->Number = 1;
  for (int i = 0; i < mt0->nBonds; i++) {
    mt0->Bond[i][0] = i;
    mt0->Bond[i][1] = i + 1;
    mt0->Bond[i][2] = 0;
  }
  for (int i = 0; i < mt0->nAngles; i++) {
    mt0->Angle[i][0] = i;
    mt0->Angle[i][1] = i + 1;
    mt0->Angle[i][2] = i + 2;
    mt0->Angle[i][3] = 0;
  } //}}}
  // create Molecule struct //{{{
  MOLECULE *mol0 = &System.Molecule[0];
  InitMolecule(mol0);
  mol0->Type = 0;
  mol0->Index = 0;
  mol0->Bead = calloc(mt0->nBeads, sizeof *mol0->Bead); //}}}
  // fill the System (no coordinates yet) //{{{
  double every = 1 / alpha;
  int first = every / 2;
  printf("%d\n", first);
  // bonded beads
  int count1 = 0, count2 = 0;
  int bt_G = FindBeadType("G", System);
  int bt_GH = FindBeadType("GH", System);
  for (int i = 0; i < Count->Bonded; i++) {
    BEAD *b = &System.Bead[i];
    InitBead(b);
    if (i > 0 && ((double)(count1) / i) < alpha) {
      b->Type = bt_G;
      count1++;
    } else { // neutral bead
      b->Type = bt_GH;
      count2++;
    }
    System.BeadCoor[i] = i;
    b->Molecule = 0;
    mt0->Bead[i] = b->Type;
    mol0->Bead[i] = i;
  }
  // unbonded beads (counterions)
  for (int i = Count->Bonded; i < Count->Bead; i++) {
    BEAD *b = &System.Bead[i];
    InitBead(b);
    b->Type = FindBeadType("C", System);
    b->Molecule = -1;
    System.BeadCoor[i] = i;
  }
  FinishSystem(&System); //}}}
  if ((count1 - System.BeadType[bt_G].Number) > 0) {
    for (int i = 0; i < (count1 - System.BeadType[bt_G].Number); i++) {
      printf("Huh? %d\n", i);
    }
  } else if ((count1 - System.BeadType[bt_G].Number) < 0) {
    for (int i = 0; i < (System.BeadType[bt_G].Number - count1); i++) {
      printf("Ah? %d\n", i);
    }
  }

  printf("%d %d %d\n", System.BeadType[0].Number, System.BeadType[1].Number, System.BeadType[2].Number);
  printf("%d %d\n", count1, count2);

  // print initial stuff to output coordinate file //{{{
  if (fw_coor.type == VCF_FILE) {
    PrintByline(fw_coor.name, argc, argv);
  } else if (fw_coor.type == VTF_FILE) {
    WriteStructure(fw_coor, System, -1, false, argc, argv);
  } else { // ensure it's a new file
    FILE *out = OpenFile(fw_coor.name, "w");
    fclose(out);
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  for (int step = 1; step <= 1000; step++) {
    fprintf(stdout, "\rChain %4d", step);
    fflush(stdout);
    // generate chain bead coordinates //{{{
    int id1, id2, id3;
    id1 = mol0->Bead[0];
    id2 = mol0->Bead[1];
    double dist = System.BondType[0].b;
    for (int dd = 0; dd < 3; dd++) {
      System.Bead[id1].Position[dd] = 0;
      System.Bead[id2].Position[dd] = 0;
    }
    System.Bead[id2].Position[2] = dist;
    for (int i = 2; i < (mt0->nBeads); i++) {
      id1 = mol0->Bead[i-2];
      id2 = mol0->Bead[i-1];
      id3 = mol0->Bead[i];
      BEAD *b1 = &System.Bead[id1],
           *b2 = &System.Bead[id2],
           *b3 = &System.Bead[id3];
      // vector from two known coordinates
      double u[3], ul = 0;
      for (int dd = 0; dd < 3; dd++) {
        u[dd] = b1->Position[dd] - b2->Position[dd];
        ul += Square(u[dd]);
      }
      ul = sqrt(ul);
      // create random vector with given length (well, dist) and angle
      double angle, v[3];
      do {
        // random x component: <-dist,dist>
        v[0] = (double)(rand()) / ((double)(RAND_MAX) + 1);
        v[0] = dist * (2 * v[0] - 1);
        // random y component <-(dist^2-c1^2)^0.5,(dist^2-c1^2)^0.5>
        double v1max = sqrt(Square(dist) - Square(v[0]));
        v[1] = (double)(rand()) / ((double)(RAND_MAX) + 1);
        v[1] = v1max * (2 * v[1] - 1);
        // z component to have |v|=dist
        v[2] = sqrt(Square(dist) - (Square(v[0]) + Square(v[1])));
        if (rand() > (RAND_MAX / 2)) {
          v[2] *= -1;
        }
        double vl = 0;
        for (int dd = 0; dd < 3; dd++) {
          vl += Square(v[dd]);
        }
        vl = sqrt(vl);

        angle = acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (ul * vl));
        angle *= 180 / PI;
        // printf(" angle = %lf (%lf)\n", angle, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (ul * vl));
      } while (angle < 155 || angle > 165);
      // generate new bead coordinate
      for (int dd = 0; dd < 3; dd++) {
        b3->Position[dd] = b2->Position[dd] + v[dd];
      }
    } //}}}
    // place centre of mass into the simulation box middle //{{{
    double com[3];
    CentreOfMass(mt0->nBeads, mol0->Bead, System, com);
    for (int i = 0; i < mt0->nBeads; i++) {
      int id = mol0->Bead[i];
      for (int dd = 0; dd < 3; dd++) {
        System.Bead[id].Position[dd] -= com[dd];
        System.Bead[id].Position[dd] += System.Box.Length[dd] / 2;
      }
    } //}}}
    // generate randomly counterion coordinates //{{{
    for (int i = Count->Bonded; i < Count->Bead; i++) {
      BEAD *b = &System.Bead[i];
      for (int dd = 0; dd < 3; dd++) {
        b->Position[dd] = (double)(rand()) / ((double)(RAND_MAX) + 1);
        b->Position[dd] *= System.Box.Length[dd];
      }
    } //}}}
    WriteTimestepAll(fw_coor, System, step, argc, argv);
  }
  fprintf(stdout, "\nDone\n");

  if (opt->fw.name[0] != '\0') {
    WriteOutputAll(System, opt->fw, false, -1, argc, argv);
  }

  // free memory
  FreeSystem(&System);
  free(opt);

  return 0;
}
