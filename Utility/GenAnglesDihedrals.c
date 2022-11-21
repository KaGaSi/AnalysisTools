#include "../AnalysisTools.h"
int *InFile;

// TODO: add reasonable output formatting such as - "type" <ids> "params" for
//       FIELD; <number> "type" <ids> for lammps data file
// TODO: speaking of FIELD, sum up all dihedrals in that case
// TODO: option to print to output file instead of screen

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
The utility calculates all possible angles and dihedrals (proper and \
improper) for a molecule from its bonds and prints their bead indices. \
For now, the molecule of interest must be the first one in the provided \
structure file; other molecules are ignored.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> [options]\n\n", cmd);
  fprintf(ptr, "   <input>   input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -v          verbose output\n");
  fprintf(ptr, "      -h          print this help and exit\n");
  fprintf(ptr, "      --version   print version number and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 1; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <input> - input structure file (must end with .vsf or .vtf) //{{{
  char input_vsf[LINE] = "";
  snprintf(input_vsf, LINE, "%s", argv[++count]);
  // test if <input> ends with '.vsf' or '.vtf'
  int ext = 2;
  char extension[2][EXTENSION];
  strcpy(extension[0], ".vsf");
  strcpy(extension[1], ".vtf");
  if (ErrorExtension(input_vsf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data //{{{
  // -c option - use a coordinate file //{{{
  char input_coor[LINE] = "";
  if (FileOption(argc, argv, "-c", input_coor, LINE)) {
    exit(1);
  }
  bool vtf = false;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  if (input_coor[0] != '\0') {
    int test;
    if ((test=ErrorExtension(input_coor, ext, extension)) == -1) {
      Help(argv[0], true);
      exit(1);
    } else if (test == 1) {
      vtf = true;
    }
  } //}}}
  bool verbose = BoolOption(argc, argv, "-v");
  //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // determine angles and dihedrals and impropers for polyacene //{{{
  // angles //{{{
  MoleculeType[0].nAngles = 0;
  MoleculeType[0].Angle = malloc(sizeof *MoleculeType[0].Angle * 1);
  Counts.TypesOfAngles = 0;
  // go through all beads, considering every one as a central one for an angle
  for (int i = 0; i < MoleculeType[0].nBeads; i++) {
    int id1 = -1,
        id2 = Molecule[0].Bead[i],
        id3 = -1;
    // find a bond that shares a bead with the central one for the angle
    for (int j = 0; j < MoleculeType[0].nBonds; j++) {
      int idx = MoleculeType[0].Bond[j][0];
      idx = Molecule[0].Bead[idx];
      int idy = MoleculeType[0].Bond[j][1];
      idy = Molecule[0].Bead[idy];
      if (id2 == idx || id2 == idy) {
        if (id2 == idx) {
          id1 = idy;
        } else {
          id1 = idx;
        }
        // find another bond sharing a bead with the central one for the angle
        for (int k = 0; k < MoleculeType[0].nBonds; k++) {
          if (j != k) {
            int idx = MoleculeType[0].Bond[k][0];
            idx = Molecule[0].Bead[idx];
            int idy = MoleculeType[0].Bond[k][1];
            idy = Molecule[0].Bead[idy];
            if (id2 == idx || id2 == idy) {
              if (id2 == idx) {
                id3 = idy;
              } else {
                id3 = idx;
              }
              // if this angle is already there, ignore it
              for (int l = 0; l < MoleculeType[0].nAngles; l++) {
                int ida = MoleculeType[0].Angle[l][0],
                    idb = MoleculeType[0].Angle[l][1],
                    idc = MoleculeType[0].Angle[l][2];
                if ((id1 == ida && id2 == idb && id3 == idc) ||
                    (id1 == idc && id2 == idb && id3 == ida)) {
                  id3 = -1;
                  break;
                }
              }
              // new angle
              if (id3 != -1) {
                int n = MoleculeType[0].nAngles;
                MoleculeType[0].nAngles++;
                MoleculeType[0].Angle = realloc(MoleculeType[0].Angle,
                                                sizeof *MoleculeType[0].Angle *
                                                MoleculeType[0].nAngles);
                MoleculeType[0].Angle[n][0] = id1;
                MoleculeType[0].Angle[n][1] = id2;
                MoleculeType[0].Angle[n][2] = id3;
                MoleculeType[0].Angle[n][3] = -1;
                id3 = -1;
              }
            }
          }
        }
      }
    }
  }
  SortAngles(MoleculeType[0].Angle, MoleculeType[0].nAngles); //}}}
  // dihedrals //{{{
  MoleculeType[0].nDihedrals = 0;
  MoleculeType[0].Dihedral = malloc(sizeof *MoleculeType[0].Dihedral * 1);
  Counts.TypesOfDihedrals = 0;
  // go through all bonds, considering every one as a central one for a dihedral
  for (int i = 0; i < MoleculeType[0].nBonds; i++) {
    int id1 = -1;
    int id2 = MoleculeType[0].Bond[i][0];
    id2 = Molecule[0].Bead[id2];
    int id3 = MoleculeType[0].Bond[i][1];
    id3 = Molecule[0].Bead[id3];
    int id4 = -1;
    // find a bond that contains the first bead of the initial bond
    for (int j = 0; j < MoleculeType[0].nBonds; j++) {
      if (i != j) {
        int idx = MoleculeType[0].Bond[j][0];
        idx = Molecule[0].Bead[idx];
        int idy = MoleculeType[0].Bond[j][1];
        idy = Molecule[0].Bead[idy];
        if (id2 == idx || id2 == idy) {
          if (id2 == idx) {
            id1 = idy;
          } else {
            id1 = idx;
          }
          // find a bond that contains the second bead of the initial bond
          for (int k = 0; k < MoleculeType[0].nBonds; k++) {
            if (i != k && j != k) {
              idx = MoleculeType[0].Bond[k][0];
              idx = Molecule[0].Bead[idx];
              idy = MoleculeType[0].Bond[k][1];
              idy = Molecule[0].Bead[idy];
              if (id3 == idx) {
                id4 = idy;
              } else if (id3 == idy) {
                id4 = idx;
              }
              // new dihedral -- no need to check whether it already exists
              if (id4 != -1) {
                int n = MoleculeType[0].nDihedrals;
                MoleculeType[0].nDihedrals++;
                MoleculeType[0].Dihedral = realloc(MoleculeType[0].Dihedral,
                                                   sizeof *MoleculeType[0].Dihedral *
                                                   MoleculeType[0].nDihedrals);
                MoleculeType[0].Dihedral[n][0] = id1;
                MoleculeType[0].Dihedral[n][1] = id2;
                MoleculeType[0].Dihedral[n][2] = id3;
                MoleculeType[0].Dihedral[n][3] = id4;
                MoleculeType[0].Dihedral[n][4] = -1;
                id4 = -1;
              }
            }
          }
        }
      }
    }
  }
  SortDihImp(MoleculeType[0].Dihedral, MoleculeType[0].nDihedrals); //}}}
  // impropers //{{{
  int nImpropers = 0;
  int (*Improper)[5] = malloc(sizeof *Improper * 1);
  // pick a central bead for an improper dihedral
  for (int i = 0; i < MoleculeType[0].nBeads; i++) {
    int id1 = -1,
        id2 = Molecule[0].Bead[i],
        id3 = -1,
        id4 = -1;
    // find a bond containing the central bead
    for (int j = 0; j < MoleculeType[0].nBonds; j++) {
      int idx = MoleculeType[0].Bond[j][0];
      idx = Molecule[0].Bead[idx];
      int idy = MoleculeType[0].Bond[j][1];
      idy = Molecule[0].Bead[idy];
      if (id2 == idx || id2 == idy) {
        if (id2 == idx) {
          id1 = idy;
        } else {
          id1 = idx;
        }
        // find a second bond containing the central bead
        for (int k = (j+1); k < MoleculeType[0].nBonds; k++) {
          int idx = MoleculeType[0].Bond[k][0];
          idx = Molecule[0].Bead[idx];
          int idy = MoleculeType[0].Bond[k][1];
          idy = Molecule[0].Bead[idy];
          if (id2 == idx || id2 == idy) {
            if (id2 == idx) {
              id3 = idy;
            } else {
              id3 = idx;
            }
            // find a second bond containing the central bead
            for (int l = (k+1); l < MoleculeType[0].nBonds; l++) {
              int idx = MoleculeType[0].Bond[l][0];
              idx = Molecule[0].Bead[idx];
              int idy = MoleculeType[0].Bond[l][1];
              idy = Molecule[0].Bead[idy];
              if (id2 == idx || id2 == idy) {
                if (id2 == idx) {
                  id4 = idy;
                } else {
                  id4 = idx;
                }
                /*
            // if this angle is already there, ignore it
            for (int l = 0; l < MoleculeType[0].nAngles; l++) {
              int ida = MoleculeType[0].Angle[l][0],
                  idb = MoleculeType[0].Angle[l][1],
                  idc = MoleculeType[0].Angle[l][2];
              if ((id1 == ida && id2 == idb && id3 == idc) ||
                  (id1 == idc && id2 == idb && id3 == ida)) {
                id3 = -1;
                break;
              }
            }
            */
                // new improper
                if (id4 != -1) {
                  int n = nImpropers;
                  nImpropers++;
                  Improper = realloc(Improper, sizeof *Improper * nImpropers);
                  Improper[n][0] = id1;
                  Improper[n][1] = id2;
                  Improper[n][2] = id3;
                  Improper[n][3] = id4;
                  Improper[n][4] = -1;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  SortDihImp(Improper, nImpropers); //}}}
  // add impropers to dihedrals //{{{
  MoleculeType[0].nDihedrals += nImpropers;
  MoleculeType[0].Dihedral = realloc(MoleculeType[0].Dihedral,
                                     sizeof *MoleculeType[0].Dihedral *
                                     MoleculeType[0].nDihedrals);
  for (int i = 0; i < nImpropers; i++) {
    int id = MoleculeType[0].nDihedrals - nImpropers + i;
    MoleculeType[0].Dihedral[id][0] = Improper[i][0];
    MoleculeType[0].Dihedral[id][1] = Improper[i][1];
    MoleculeType[0].Dihedral[id][2] = Improper[i][2];
    MoleculeType[0].Dihedral[id][3] = Improper[i][3];
    MoleculeType[0].Dihedral[id][4] = Improper[i][4];
  } //}}}
  //}}}

  // print information /{{{
  if (verbose) {
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // print angles & dihedrals //{{{
  if (MoleculeType[0].nAngles > 0) {
    printf("angles %d\n", MoleculeType[0].nAngles);
    for (int i = 0; i < MoleculeType[0].nAngles; i++) {
      printf("harm %4d %4d %4d   0 0\n", MoleculeType[0].Angle[i][0]+1,
                                         MoleculeType[0].Angle[i][1]+1,
                                         MoleculeType[0].Angle[i][2]+1);
    }
  }
  if ((MoleculeType[0].nDihedrals-nImpropers) > 0) {
    printf("dihedrals %d\n", MoleculeType[0].nDihedrals-nImpropers);
    for (int i = 0; i < (MoleculeType[0].nDihedrals-nImpropers); i++) {
      printf("harm %4d %4d %4d %4d   0 0\n", MoleculeType[0].Dihedral[i][0]+1,
                                             MoleculeType[0].Dihedral[i][1]+1,
                                             MoleculeType[0].Dihedral[i][2]+1,
                                             MoleculeType[0].Dihedral[i][3]+1);
    }
  }
  if (nImpropers > 0) {
    printf("impropers %d\n", nImpropers);
    for (int i = 0; i < nImpropers; i++) {
      int id = MoleculeType[0].nDihedrals - nImpropers + i;
      printf("harm %4d %4d %4d %4d   0 0\n", MoleculeType[0].Dihedral[id][0]+1,
                                             MoleculeType[0].Dihedral[id][1]+1,
                                             MoleculeType[0].Dihedral[id][2]+1,
                                             MoleculeType[0].Dihedral[id][3]+1);
    }
  } //}}}

  // free memory - to make valgrind happy
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(Improper);

  return 0;
}
