#include "AnalysisTools.h"

// VerboseOutput() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput(char *input_vcf, COUNTS Counts, VECTOR BoxLength,
                   BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {

  putchar('\n');
  if (BoxLength.x != -1) {
    fprintf(stdout, "Box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } else {
    fprintf(stdout, "Unknown box size because no coordinate file is provided.\n\n");
  }
  PrintCounts(Counts);
  PrintBeadType(Counts, BeadType);
  PrintMoleculeType(Counts, BeadType, MoleculeType);
  putchar('\n');
} //}}}

// PrintCounts()  //{{{
/**
 * Function printing Counts structure.
 */
void PrintCounts(COUNTS Counts) {
  fprintf(stdout, "Counts = {\n");
  fprintf(stdout, "  .TypesOfBeads     = %d,\n", Counts.TypesOfBeads);
  fprintf(stdout, "  .Bonded           = %d,\n", Counts.Bonded);
  fprintf(stdout, "  .Unbonded         = %d,\n", Counts.Unbonded);
  fprintf(stdout, "  .Beads            = %d,\n", Counts.Beads);
  if (Counts.BeadsInVsf != Counts.Beads) {
    fprintf(stdout, "  .BeadsInVsf       = %d,\n", Counts.BeadsInVsf);
  }
  fprintf(stdout, "  .TypesOfMolecules = %d,\n", Counts.TypesOfMolecules);
  fprintf(stdout, "  .Molecules        = %d,\n", Counts.Molecules);
  if (Counts.TypesOfBonds != -1) {
    fprintf(stdout, "  .TypesOfBonds     = %d,\n", Counts.TypesOfBonds);
  }
  if (Counts.TypesOfAngles != -1) {
    fprintf(stdout, "  .TypesOfAngles    = %d,\n", Counts.TypesOfAngles);
  }
  fprintf(stdout, "}\n\n");
} //}}}

// PrintBeadType()  //{{{
/**
 * Function printing BeadType structure.
 */
void PrintBeadType(COUNTS Counts, BEADTYPE *BeadType) {
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(stdout, "BeadType[%2d] = {", i);
    fprintf(stdout, ".Name =%10s, ", BeadType[i].Name);
    fprintf(stdout, ".Number =%7d, ", BeadType[i].Number);
    fprintf(stdout, ".Charge =%9.5f, ", BeadType[i].Charge);
    fprintf(stdout, ".Mass = %.5f}\n", BeadType[i].Mass);
//  fprintf(stdout, "Use = %3s, ", BeadType[i].Use? "Yes":"No");
//  fprintf(stdout, "Write = %3s}\n", BeadType[i].Write? "Yes":"No");
  }
  putchar('\n');
} //}}}

// PrintMoleculeType()  //{{{
/**
 * Function printing MoleculeType structure.
 */
void PrintMoleculeType(COUNTS Counts, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(stdout, "MoleculeType[%2d] = {\n", i);
    fprintf(stdout, "  .Name    = %s,\n", MoleculeType[i].Name);
    fprintf(stdout, "  .Number  = %d,\n", MoleculeType[i].Number);
    // print bead types (list all beads) //{{{
    fprintf(stdout, "  .nBeads  = %d,\n", MoleculeType[i].nBeads);
    fprintf(stdout, "  .Bead    = {");
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].Bead[j]].Name);
    }
    fprintf(stdout, "},\n"); //}}}
    // print bonds if there are any //{{{
    if (MoleculeType[i].nBonds > 0) {
      fprintf(stdout, "  .nBonds  = %d,\n", MoleculeType[i].nBonds);
      fprintf(stdout, "  .Bond    = {");
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d", MoleculeType[i].Bond[j][0]+1, MoleculeType[i].Bond[j][1]+1);
        if (MoleculeType[i].Bond[j][2] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Bond[j][2]+1);
        }
      }
    } //}}}
    // print angles if there are any //{{{
    if (MoleculeType[i].nAngles > 0) {
      fprintf(stdout, "},\n  .nAngles = %d,\n  .Angle   = {", MoleculeType[i].nAngles);
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d", MoleculeType[i].Angle[j][0]+1,
                                    MoleculeType[i].Angle[j][1]+1,
                                    MoleculeType[i].Angle[j][2]+1);
        if (MoleculeType[i].Angle[j][3] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Angle[j][3]+1);
        }
      }
    } //}}}
    // print bead types (just the which are present) //{{{
    fprintf(stdout, "},\n  .nBTypes = %d,\n  .BType   = {", MoleculeType[i].nBTypes);
    for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].BType[j]].Name);
    } //}}}
    fprintf(stdout, "},\n  .Mass    = %.5f,\n}\n", MoleculeType[i].Mass);
  }
} //}}}

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead) {
  fprintf(stdout, "Beads - <i> (<Bead[i].Index>; <Index[i]>)\n");
  for (int i = 0; i < Counts.Beads; i++) {
    int type = Bead[i].Type;
    fprintf(stdout, "   %6d (%6d; %6d) %8s molecule: ", i, Bead[i].Index, Index[i], BeadType[type].Name);
    if (Bead[i].Molecule == -1) {
      fprintf(stdout, "None\n");
    } else {
      fprintf(stdout, "%6d\n", Bead[i].Molecule+1);
    }
  }
} //}}}

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(COUNTS Counts, int *Index,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                   BEADTYPE *BeadType, BEAD *Bead) {
  fprintf(stdout, "Molecules\n");
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    fprintf(stdout, "Molecule %3d (%s):\n", i+1, MoleculeType[type].Name);
    fprintf(stdout, " BEAD INDICES (%d): internal (in vsf)\n", MoleculeType[type].nBeads);
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      fprintf(stdout, "   %d (%d)\n", Molecule[i].Bead[j], Bead[Molecule[i].Bead[j]].Index);
    }
    fprintf(stdout, " BONDS (%d): internal (in vsf)\n", MoleculeType[type].nBonds);
    for (int j = 0; j < MoleculeType[type].nBonds; j++) {
      int bead1 = MoleculeType[type].Bond[j][0];
      int bead1_1 = Molecule[i].Bead[bead1];
      int bead2 = MoleculeType[type].Bond[j][1];
      int bead2_1 = Molecule[i].Bead[bead2];
      fprintf(stdout, "   %d-%d (%d-%d)\n", bead1_1, bead2_1, Bead[bead1_1].Index, Bead[bead2_1].Index);
    }
  }
  fprintf(stdout, "\n");
} //}}}

// PrintAggregate() //{{{
/**
 * Function printing Aggregate structure.
 */
void PrintAggregate(COUNTS Counts, int *Index,
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                    BEAD *Bead, BEADTYPE *BeadType, AGGREGATE *Aggregate) {
  fprintf(stdout, "Aggregates: %d\n", Counts.Aggregates);
  for (int i = 0; i < Counts.Aggregates; i++) {
    // print molecules
    fprintf(stdout, " %d mols:", Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      int type = Molecule[mol].Type;
      fprintf(stdout, " %d (%d)", mol, type);
      if (j != (Aggregate[i].nMolecules-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print bonded beads
    fprintf(stdout, " %d bonded beads:", Aggregate[i].nBeads);
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
      if (j != (Aggregate[i].nBeads-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print monomeric beads
    fprintf(stdout, " %d free beads:", Aggregate[i].nMonomers);
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int bead = Aggregate[i].Monomer[j];
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
      if (j != (Aggregate[i].nMonomers-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
  }
} //}}}

// PrintBondTypes() //{{{
void PrintBondTypes(COUNTS Counts, PARAMS *bond_type) {
  for (int i = 0; i < Counts.TypesOfBonds; i++) {
    fprintf(stdout, "bond %2d: k = %lf, r_0 = %lf\n", i+1, bond_type[i].a, bond_type[i].b);
  }
  putc('\n', stdout);
} //}}}

// PrintAngleTypes() //{{{
void PrintAngleTypes(COUNTS Counts, PARAMS *angle_type) {
  for (int i = 0; i < Counts.TypesOfAngles; i++) {
    fprintf(stdout, "angle %2d: k = %lf, r_0 = %lf\n", i+1, angle_type[i].a, angle_type[i].b);
  }
  putc('\n', stdout);
} //}}}

// FindBeadType() //{{{
/* Function to identify type of bead from its name; returns -1 on non-existent
 * bead name.
 */
int FindBeadType(char *name, COUNTS Counts, BEADTYPE *BeadType) {
  int type;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in BeadType struct
  return -1;
} //}}}

// FindMoleculeType() //{{{
int FindMoleculeType(char *name, COUNTS Counts, MOLECULETYPE *MoleculeType) {
  int type;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return(-1);
} //}}}

// Distance() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength) {
  // distance vector
  VECTOR rij;
  rij.x = id1.x - id2.x;
  rij.y = id1.y - id2.y;
  rij.z = id1.z - id2.z;
  // remove periodic boundary conditions in x-direction
  while (rij.x >= (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x < -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y >= (BoxLength.y/2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y < -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z >= (BoxLength.z/2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z < -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;
  return rij;
} //}}}

// RemovePBCMolecules() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them.
 */
void RemovePBCMolecules(COUNTS Counts, VECTOR BoxLength,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      (*Bead)[Molecule[i].Bead[j]].Flag = false; // no beads moved yet
    }
    // first bead in the first bond is considered moved
    if (MoleculeType[type].nBonds > 0) {
      (*Bead)[Molecule[i].Bead[MoleculeType[type].Bond[0][0]]].Flag = true;
      bool done = false;
      int test = 0; // if too many loops, just leave the loop with error
      while (!done && test < 1000) {
        for (int j = 0; j < MoleculeType[type].nBonds; j++) {
          int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
          int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
          // move id1, if id2 is moved already
          if (!(*Bead)[id1].Flag && (*Bead)[id2].Flag) {
            VECTOR dist = Distance((*Bead)[id2].Position, (*Bead)[id1].Position, BoxLength);
            (*Bead)[id1].Position.x = (*Bead)[id2].Position.x - dist.x;
            (*Bead)[id1].Position.y = (*Bead)[id2].Position.y - dist.y;
            (*Bead)[id1].Position.z = (*Bead)[id2].Position.z - dist.z;
            (*Bead)[id1].Flag = true;
          // move id2, if id1 was moved already
          } else if ((*Bead)[id1].Flag && !(*Bead)[id2].Flag) {
            VECTOR dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position, BoxLength);

            (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
            (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
            (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
            (*Bead)[id2].Flag = true;
          }
        }

        // exit if all beads have moved
        done = true;
        for (int j = 1; j < MoleculeType[type].nBeads; j++) {
          if (!(*Bead)[Molecule[i].Bead[j]].Flag) {
            done = false;
            break;
          }
        }
        test++;
      }
      if (test == 1000) {
        fprintf(stderr, "\033[1;33m");
        fprintf(stderr, "\nWarning: unable connect molecule %s (resid %d)\n\n", MoleculeType[type].Name, i+1);
        fprintf(stderr, "\033[0m");
      }

      // put molecule's centre of mass into the simulation box //{{{
      VECTOR com = CentreOfMass(MoleculeType[type].nBeads, Molecule[i].Bead, *Bead, BeadType);
      // by how many BoxLength's should com be moved?
      // for distant molecules - it shouldn't happen, but better safe than sorry
      INTVECTOR move;
      move.x = com.x / BoxLength.x;
      move.y = com.y / BoxLength.y;
      move.z = com.z / BoxLength.z;
      if (com.x < 0) {
        move.x--;
      }
      if (com.y < 0) {
        move.y--;
      }
      if (com.z < 0) {
        move.z--;
      }
      for (int j = 0; j < MoleculeType[type].nBeads; j++) {
        int bead = Molecule[i].Bead[j];
        (*Bead)[bead].Position.x -= move.x * BoxLength.x;
        (*Bead)[bead].Position.y -= move.y * BoxLength.y;
        (*Bead)[bead].Position.z -= move.z * BoxLength.z;
      } //}}}
    }
  }
} //}}}

// RemovePBCAggregates() //{{{
/**
 * Function to remove periodic boundary conditions from all aggregates,
 * thus joining them.
 */
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {

  bool *moved = malloc(Counts.Molecules*sizeof(bool));

  // go through all aggregates larger than unimers and put all molecules together //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    // negate moved array, but the first molecule is not to move //{{{
    for (int j = 1; j < Counts.Molecules; j++) {
      moved[j] = false;
    }
    moved[0] = true; //}}}

    bool done = false;
    int test = 0; // if too many loops, just exit with error
    while (!done && test < 1000) {

      // go through all molecule pairs
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        for (int k = 0; k < Aggregate[i].nMolecules; k++) {
          int mol1 = Aggregate[i].Molecule[j],
              mol2 = Aggregate[i].Molecule[k],
              mol1_type = Molecule[mol1].Type,
              mol2_type = Molecule[mol2].Type;

          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          if (moved[j] && !moved[k]) { // automatically follows that j != k
            // go through all bead pairs in the two molecules
            for (int l = 0; l < MoleculeType[mol1_type].nBeads; l++) {
              for (int m = 0; m < MoleculeType[mol2_type].nBeads; m++) {
                int bead1 = Molecule[mol1].Bead[l];
                int bead2 = Molecule[mol2].Bead[m];

                // use only bead types that were used to assign molecules to aggregates
                if (BeadType[(*Bead)[bead1].Type].Use &&
                    BeadType[(*Bead)[bead2].Type].Use) {

                  // calculate distance between 'bead1' and 'bead2'
                  VECTOR dist = Distance((*Bead)[bead1].Position, (*Bead)[bead2].Position, BoxLength);
                  dist.x = Length(dist);

                  // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
                  if (dist.x <= distance) {

                    // distance vector between 'bead1' and 'bead2' //{{{
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z; //}}}
                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    while (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    }
                    while (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    while (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    }
                    while (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    while (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    }
                    while (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    } //}}}

                    moved[k] = true;

                    // skip remainder of 'mol2' (or 'k')
                    break;
                  }
                }
              }
              // if molekule 'k' (or 'mol2') has been moved, skip also remainder of molecules 'mol1'
              if (moved[k]) {
                break;
              }
            }
          }
        }
      }

      // check if all molecules have moved //{{{
      done = true;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (!moved[j]) {
          done = false;
          break;
        }
      } //}}}
      test++;
    }
    if (test == 1000) {
      fprintf(stderr, "\033[1;33m");
      fprintf(stderr, "\nWarning: unable connect aggregate containing resids:\n");
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        fprintf(stderr, " %d", Aggregate[i].Molecule[j]);
        if (j != (Aggregate[i].nMolecules)) {
          fprintf(stderr, ",");
        } else {
          fprintf(stderr, "\n");
        }
      }
      fprintf(stderr, "\033[0m");
    }
  }
  free(moved); //}}}

  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    VECTOR com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead,
                              *Bead, BeadType);
    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    INTVECTOR move;
    move.x = com.x / BoxLength.x;
    move.y = com.y / BoxLength.y;
    move.z = com.z / BoxLength.z;
    if (com.x < 0) {
      move.x--;
    }
    if (com.y < 0) {
      move.y--;
    }
    if (com.z < 0) {
      move.z--;
    }
    // move all the beads
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * BoxLength.x;
      (*Bead)[bead].Position.y -= move.y * BoxLength.y;
      (*Bead)[bead].Position.z -= move.z * BoxLength.z;
    }
  } //}}}
} //}}}

// RestorePBC() //{{{
/**
 * Function to restore removed periodic boundary conditions. Used also in case
 * of cell linked lists, because they need coordinates <0, BoxLength>.
 */
void RestorePBC(COUNTS Counts, VECTOR BoxLength, BEAD **Bead) {

  for (int i = 0; i < (Counts.Unbonded+Counts.Bonded); i++) {
    // x direction
    while ((*Bead)[i].Position.x >= BoxLength.x) {
      (*Bead)[i].Position.x -= BoxLength.x;
    }
    while ((*Bead)[i].Position.x < 0) {
      (*Bead)[i].Position.x += BoxLength.x;
    }
    // y direction
    while ((*Bead)[i].Position.y >= BoxLength.y) {
      (*Bead)[i].Position.y -= BoxLength.y;
    }
    while ((*Bead)[i].Position.y < 0) {
      (*Bead)[i].Position.y += BoxLength.y;
    }
    // z direction
    while ((*Bead)[i].Position.z >= BoxLength.z) {
      (*Bead)[i].Position.z -= BoxLength.z;
    }
    while ((*Bead)[i].Position.z < 0) {
      (*Bead)[i].Position.z += BoxLength.z;
    }
  }
} //}}}

// CentreOfMass() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType) {
  VECTOR com = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    com.x += Bead[list[i]].Position.x * BeadType[Bead[list[i]].Type].Mass;
    com.y += Bead[list[i]].Position.y * BeadType[Bead[list[i]].Type].Mass;
    com.z += Bead[list[i]].Position.z * BeadType[Bead[list[i]].Type].Mass;
  }
  com.x /= n;
  com.y /= n;
  com.z /= n;
  return com;
} //}}}

// GeomCentre() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
VECTOR GeomCentre(int n, int *list, BEAD *Bead) {
  VECTOR cog = {0, 0, 0};
  for (int i = 0; i < n; i++) {
    cog.x += Bead[list[i]].Position.x;
    cog.y += Bead[list[i]].Position.y;
    cog.z += Bead[list[i]].Position.z;
  }
  cog.x /= n;
  cog.y /= n;
  cog.z /= n;
  return cog;
} //}}}

// Gyration() //{{{
/**
 * Function to calculate the principle moments of the gyration tensor.
 */
// // jacobi - numerical recipes //{{{
// #define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
// a[k][l]=h+s*(g-h*tau);
// void jacobi(float **a, int n, float d[], float **v, int *nrot) {
// /* Computes all eignvalues and eignvectors of a real symmetric matrix a[1..n][1..n]. On output 
// elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a v[1..n][1..n] is a matrix whose columns contain, on output, the normalised eigenvectors of a. nrot 
// returns the number of JUacobi rotations that were required.
// */
//   int j, iq, ip, i;
//   float tresh, theta, tau, t, sm, s, h, g, c;
//   float b[n], z[n];
//   /* Initialise to the identity matrix */
//   for (ip = 0; ip < n; ip++) {
//     for (iq = 0; iq < n; iq++) {
//       v[ip][iq] = 0.0;
//     }
//     v[ip][ip] = 1.0;
//   }
//   for (ip = 0; ip < n; ip++) {
//     b[ip] = d[ip] = a[ip][ip];
//     z[ip] = 0.0;
//   }
//
//   *nrot = 0;
//   for (i = 0; i < 50; i++) {
//     sm = 0;
//     for (ip = 0; ip < (n-1); ip++) {
//       for (ip = (ip+1); iq < n; iq++) {
//         sm += fabs(a[ip][iq]);
//       }
//     }
//     if (sm == 0.0) {
//       return;
//     }
//     if (i < 4) {
//       tresh = 0.2 * sm / SQR(n);
//     } else {
//       tresh = 0.0;
//     }
//     for (ip = 0; ip < (n-1); ip++) {
//       for (iq = (ip+1); iq < n; iq++) {
//         g = 100.0 * fabs(a[ip][iq]);
//         if (i > 4 && (float)(fabs(d[ip]) + g) == (float)fabs(d[ip])
//           && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) {
//           a[ip][iq] = 0.0;
//         } else if (fabs(a[ip][iq]) > tresh) {
//           h = d[iq] - d[ip];
//           if ((float)(fabs(h)+g) == (float)fabs(h))
//             t = a[ip][iq] / h;
//           else {
//             theta = 0.5 * h / a[ip][iq];
//             t = 1.0 / (fabs(theta) + sqrt(1.0 + SQR(theta)));
//             if (theta < 0) {
//               t = -t;
//             }
//           }
//           c = 1.0 / sqrt(1 + SQR(t));
//           s = t * c;
//           tau = s / (1.0 + c);
//           h = t * a[ip][iq];
//           z[ip] -= h;
//           z[iq] += h;
//           d[ip] -= h;
//           d[iq] += h;
//           a[ip][iq] = 0.0;
//           for (j = 0; j < ip; j++) {
//             ROTATE(a, j, ip, j, iq);
//           }
//           for (j = (ip+1); j < iq; j++) {
//             ROTATE(a, ip, j, j, iq);
//           }
//           for (j = (iq+1); j < n; j++) {
//             ROTATE(a, ip, j, iq, j);
//           }
//           for (j = 0; j < n; j++) {
//             ROTATE(v, j, ip, j, iq);
//           }
//           (*nrot)++;
//         }
//       }
//     }
//     for (ip = 0; ip < n; ip++) {
//       b[ip] += z[ip];
//       d[ip] = b[ip];
//       z[ip] = 0.0;
//     }
//   }
//   printf("Too many iterations in routine Jacobi\n");
// } //}}}
// jacobi - ml //{{{
void jacobi2(double **a, int n, double d[], double **v) {
  int maxit = 1000;
  double eps = 1e-6;
  int r, s;
  double qt[n][n], amax, tsum, c, sgn_ars, z, omg, rmu, rnu, sgn_mu,
         temp1, temp2, temp3;
  for (int it = 0; it < maxit; it++) {
    amax = 0;
    tsum = 0;

    for (int i = 0; i < (n-1); i++) {
      for (int j = (i+1); j < n; j++) {
        tsum += SQR(a[i][j]);
        if (fabs(a[i][j]) > amax) {
          amax = fabs(a[i][j]);
          r = i;
          s = j;
        }
      }
    }

    if (tsum < eps) {
      for (int i = 0; i < 3; i++) {
        d[i] = a[i][i];
      }
      return;
    }

    if (a[r][r] == a[s][s]) {
      c = sqrt(2) / 2;
      sgn_ars = 1;
      if (a[r][s] < 0) {
        sgn_ars = -1;
      }
      z = sgn_ars * c;
    } else {
      omg = -1 * a[r][s];
      rmu = (a[r][r] - a[s][s]) / 2;
      rnu = sqrt(SQR(omg) + SQR(rmu));
      c = sqrt((rnu + fabs(rmu)) / (2 * rnu));
      sgn_mu = 1;
      if (rmu < 0) {
        sgn_mu = -1;
      }
      z = sgn_mu * omg / (2 * rnu * c);
    }

    for (int i = 0; i < n; i++) {
      if (i == r || i == s) {
        continue;
      }
      temp1 = c * a[i][r] - z * a[i][s];
      temp2 = c * a[i][s] + z * a[i][r];
      a[i][r] = temp1;
      a[r][i] = temp1;
      a[i][s] = temp2;
      a[s][i] = temp2;
    }
    temp3 = SQR(c) * a[r][r] + SQR(z) * a[s][s] - 2 * c * z * a[r][s];
    a[s][s] = SQR(c) * a[s][s] + SQR(z) * a[r][r] + 2 * c * z * a[r][s];
    a[r][r] = temp3;
    a[r][s] = 0;
    a[s][r] = 0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        qt[i][j] = 0;
      }
      qt[i][i] = 1;
    }
    qt[r][r] = c;
    qt[s][s] = c;
    qt[r][s] = z;
    qt[s][r] = -z;
    if (it == 0) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          v[i][j] = qt[i][j];
        }
      }
    } else {
      double mul[n][n];
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          double sum = 0;
          for (int k = 0; k < n; k++) {
            sum += v[i][k] * qt[k][j];
          }
          mul[i][j] = sum;
        }
      }
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          v[i][j] = mul[i][j];
        }
      }
    }
  }
} //}}}
void Gyration(int n, int *list, COUNTS Counts, VECTOR BoxLength,
              BEADTYPE *BeadType, BEAD **Bead, double **tensor,
              double eigenvalue[3], double **eigenvector) {
  // gyration tensor (3x3 array)
  // use long double to ensure precision -- previous problem with truncation in short chains
  struct Tensor {
    LONGVECTOR x, y, z;
  } GyrationTensor = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

  VECTOR com = GeomCentre(n, list, *Bead);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    GyrationTensor.x.x += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
    GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
    GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
    GyrationTensor.y.y += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
    GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
    GyrationTensor.z.z += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n; //}}}

  // tensor to pass back
  tensor[0][0] = GyrationTensor.x.x;
  tensor[0][1] = GyrationTensor.x.y;
  tensor[0][2] = GyrationTensor.x.z;
  tensor[1][0] = tensor[0][1];
  tensor[1][1] = GyrationTensor.y.y;
  tensor[1][2] = GyrationTensor.y.z;
  tensor[2][0] = tensor[0][2];
  tensor[2][1] = tensor[1][2];
  tensor[2][2] = GyrationTensor.z.z;

  // Jacobi
  double **a = malloc(3 * sizeof *a);
  for (int i = 0; i < 3; i++) {
    a[i] = calloc(3, sizeof *a);
  }
  a[0][0] = GyrationTensor.x.x;
  a[0][1] = GyrationTensor.x.y;
  a[0][2] = GyrationTensor.x.z;
  a[1][0] = a[0][1];
  a[1][1] = GyrationTensor.y.y;
  a[1][2] = GyrationTensor.y.z;
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  a[2][2] = GyrationTensor.z.z;

  jacobi2(a, 3, eigenvalue, eigenvector);
  VECTOR eigen;
  eigen.x = eigenvalue[0];
  eigen.y = eigenvalue[1];
  eigen.z = eigenvalue[2];
  eigen = Sort3(eigen);
  eigenvalue[0] = eigen.x;
  eigenvalue[1] = eigen.y;
  eigenvalue[2] = eigen.z;

  // sort eigenvalues - TODO

  for (int i = 0; i < 3; i++) {
    free(a[i]);
  }
  free(a);
} //}}}
#if 0
// Gyration - backup() //{{{
/**
 * Function to calculate the principle moments of the gyration tensor.
 */
// jacobi - numerical recipes //{{{
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);
void jacobi(float **a, int n, float d[], float **v, int *nrot) {
/* Computes all eignvalues and eignvectors of a real symmetric matrix a[1..n][1..n]. On output 
elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a v[1..n][1..n] is a matrix whose columns contain, on output, the normalised eigenvectors of a. nrot 
returns the number of JUacobi rotations that were required.
*/
  int j, iq, ip, i;
  float tresh, theta, tau, t, sm, s, h, g, c;
  float b[n], z[n];
  /* Initialise to the identity matrix */
  for (ip = 0; ip < n; ip++) {
    for (iq = 0; iq < n; iq++) {
      v[ip][iq] = 0.0;
    }
    v[ip][ip] = 1.0;
  }
  for (ip = 0; ip < n; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  *nrot = 0;
  for (i = 0; i < 50; i++) {
    sm = 0;
    for (ip = 0; ip < (n-1); ip++) {
      for (ip = (ip+1); iq < n; iq++) {
        sm += fabs(a[ip][iq]);
      }
    }
    if (sm == 0.0) {
      return;
    }
    if (i < 4) {
      tresh = 0.2 * sm / SQR(n);
    } else {
      tresh = 0.0;
    }
    for (ip = 0; ip < (n-1); ip++) {
      for (iq = (ip+1); iq < n; iq++) {
        g = 100.0 * fabs(a[ip][iq]);
        if (i > 4 && (float)(fabs(d[ip]) + g) == (float)fabs(d[ip])
          && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) {
          a[ip][iq] = 0.0;
        } else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
          if ((float)(fabs(h)+g) == (float)fabs(h))
            t = a[ip][iq] / h;
          else {
            theta = 0.5 * h / a[ip][iq];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + SQR(theta)));
            if (theta < 0) {
              t = -t;
            }
          }
          c = 1.0 / sqrt(1 + SQR(t));
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0.0;
          for (j = 0; j < ip; j++) {
            ROTATE(a, j, ip, j, iq);
          }
          for (j = (ip+1); j < iq; j++) {
            ROTATE(a, ip, j, j, iq);
          }
          for (j = (iq+1); j < n; j++) {
            ROTATE(a, ip, j, iq, j);
          }
          for (j = 0; j < n; j++) {
            ROTATE(v, j, ip, j, iq);
          }
          (*nrot)++;
        }
      }
    }
    for (ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  printf("Too many iterations in routine Jacobi\n");
} //}}}
// jacobi - ml //{{{
void jacobi2(double **a, int n, double d[], double **v) {
  int maxit = 1000;
  double eps = 1e-6;
  int r, s;
  double qt[n][n], amax, tsum, c, sgn_ars, z, omg, rmu, rnu, sgn_mu,
         temp1, temp2, temp3;
  for (int it = 0; it < maxit; it++) {
    amax = 0;
    tsum = 0;

    for (int i = 0; i < (n-1); i++) {
      for (int j = (i+1); j < n; j++) {
        tsum += SQR(a[i][j]);
        if (fabs(a[i][j]) > amax) {
          amax = fabs(a[i][j]);
          r = i;
          s = j;
        }
      }
    }

    if (tsum < eps) {
      for (int i = 0; i < 3; i++) {
        d[i] = a[i][i];
      }
      return;
    }

    if (a[r][r] == a[s][s]) {
      c = sqrt(2) / 2;
      sgn_ars = 1;
      if (a[r][s] < 0) {
        sgn_ars = -1;
      }
      z = sgn_ars * c;
    } else {
      omg = -1 * a[r][s];
      rmu = (a[r][r] - a[s][s]) / 2;
      rnu = sqrt(SQR(omg) + SQR(rmu));
      c = sqrt((rnu + fabs(rmu)) / (2 * rnu));
      sgn_mu = 1;
      if (rmu < 0) {
        sgn_mu = -1;
      }
      z = sgn_mu * omg / (2 * rnu * c);
    }

    for (int i = 0; i < n; i++) {
      if (i == r || i == s) {
        continue;
      }
      temp1 = c * a[i][r] - z * a[i][s];
      temp2 = c * a[i][s] + z * a[i][r];
      a[i][r] = temp1;
      a[r][i] = temp1;
      a[i][s] = temp2;
      a[s][i] = temp2;
    }
    temp3 = SQR(c) * a[r][r] + SQR(z) * a[s][s] - 2 * c * z * a[r][s];
    a[s][s] = SQR(c) * a[s][s] + SQR(z) * a[r][r] + 2 * c * z * a[r][s];
    a[r][r] = temp3;
    a[r][s] = 0;
    a[s][r] = 0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        qt[i][j] = 0;
      }
      qt[i][i] = 1;
    }
    qt[r][r] = c;
    qt[s][s] = c;
    qt[r][s] = z;
    qt[s][r] = -z;
    if (it == 0) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          v[i][j] = qt[i][j];
        }
      }
    } else {
      double mul[n][n];
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          double sum = 0;
          for (int k = 0; k < n; k++) {
            sum += v[i][k] * qt[k][j];
          }
          mul[i][j] = sum;
        }
      }
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          v[i][j] = mul[i][j];
        }
      }
    }
    printf("%lf %lf %lf\n",   v[0][0], v[0][1], v[0][2]);
    printf("%lf %lf %lf\n",   v[1][0], v[1][1], v[1][2]);
    printf("%lf %lf %lf\n\n", v[2][0], v[2][1], v[2][2]);
  }
} //}}}
// VECTOR Gyration(int n, int *list, COUNTS Counts, VECTOR BoxLength,
//                 BEADTYPE *BeadType, BEAD **Bead, double tensor[6]) {
//   // gyration tensor (3x3 array)
//   // use long double to ensure precision -- previous problem with truncation in short chains
//   struct Tensor {
//     LONGVECTOR x, y, z;
//   } GyrationTensor = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
//
//   VECTOR com = GeomCentre(n, list, *Bead);
//
//   // move centre of mass to [0,0,0] //{{{
//   for (int i = 0; i < n; i++) {
//     (*Bead)[list[i]].Position.x -= com.x;
//     (*Bead)[list[i]].Position.y -= com.y;
//     (*Bead)[list[i]].Position.z -= com.z;
//   } //}}}
//
//   // calculate gyration tensor //{{{
//   for (int i = 0; i < n; i++) {
//     GyrationTensor.x.x += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.x;
//     GyrationTensor.x.y += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.y;
//     GyrationTensor.x.z += (*Bead)[list[i]].Position.x * (*Bead)[list[i]].Position.z;
//     GyrationTensor.y.y += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.y;
//     GyrationTensor.y.z += (*Bead)[list[i]].Position.y * (*Bead)[list[i]].Position.z;
//     GyrationTensor.z.z += (*Bead)[list[i]].Position.z * (*Bead)[list[i]].Position.z;
//   }
//   GyrationTensor.x.x /= n;
//   GyrationTensor.x.y /= n;
//   GyrationTensor.x.z /= n;
//   GyrationTensor.y.y /= n;
//   GyrationTensor.y.z /= n;
//   GyrationTensor.z.z /= n; //}}}
//
//   tensor[0] = GyrationTensor.x.x;
//   tensor[1] = GyrationTensor.x.y;
//   tensor[2] = GyrationTensor.x.z;
//   tensor[3] = GyrationTensor.y.y;
//   tensor[4] = GyrationTensor.y.z;
//   tensor[5] = GyrationTensor.z.z;
//
//   // char polynomial: a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0 //{{{
//   long double a_cube = -1;
//   long double b_cube = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
//   long double c_cube = - GyrationTensor.x.x * GyrationTensor.y.y
//                        - GyrationTensor.x.x * GyrationTensor.z.z
//                        - GyrationTensor.y.y * GyrationTensor.z.z
//                        + SQR(GyrationTensor.y.z)
//                        + SQR(GyrationTensor.x.y)
//                        + SQR(GyrationTensor.x.z);
//   long double d_cube = + GyrationTensor.x.x * GyrationTensor.y.y * GyrationTensor.z.z
//                        + 2 * GyrationTensor.x.y * GyrationTensor.y.z * GyrationTensor.x.z
//                        - SQR(GyrationTensor.x.z) * GyrationTensor.y.y
//                        - SQR(GyrationTensor.x.y) * GyrationTensor.z.z
//                        - SQR(GyrationTensor.y.z) * GyrationTensor.x.x; //}}}
//
//   // first root: either 0 or Newton's iterative method to get it //{{{
//   long double root0 = 0;
//   if (fabs(d_cube) > 0.0000000001L) {
//     // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
//     long double a_deriv = 3 * a_cube;
//     long double b_deriv = 2 * b_cube;
//     long double c_deriv = c_cube;
//
//     long double root1 = 1;
//
//     while (fabs(root0-root1) > 0.0000000001L) {
//       long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) + c_cube * root0 + d_cube);
//       long double f_deriv_root0 = (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
//       root1 = root0 - f_root0 / f_deriv_root0;
//
//       // swap root0 and root1 for the next iteration
//       long double tmp = root0;
//       root0 = root1;
//       root1 = tmp;
//     }
//   } //}}}
//
//   // determine paremeters of quadratic equation a_quad * x^2 + b_quad * x + c_quad = 0 //{{{
//   // derived by division: (x^3 + (b_cube/a_cube) * x^2 + (c_cube/a_cube) * x + (d_cube/a_cube)):(x - root0)
//   long double a_quad = 1;
//   long double b_quad = b_cube / a_cube + root0;
//   long double c_quad = SQR(root0) + b_cube / a_cube * root0 + c_cube/a_cube; //}}}
//
//   // calculate & sort eigenvalues //{{{
//   LONGVECTOR eigen;
//   eigen.x = root0; // found out by Newton's method
//   // roots of the quadratic equation
//   eigen.y = (-b_quad + sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
//   eigen.z = (-b_quad - sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
//
//   VECTOR eigen2; // change LONGVECTOR to VECTOR
//   eigen2.x = eigen.x;
//   eigen2.y = eigen.y;
//   eigen2.z = eigen.z;
//   eigen2 = Sort3(eigen2); //}}}
//
//   return eigen2;
// } //}}}
#endif

// EvaluateContacts() //{{{
/**
 * Function evaluating contacts for aggregate detection for Aggregate and
 * Aggregate-NotSameBeads utilities.
 */
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule,
                      int contacts, int **contact) {
  // first molecule
  for (int i = 1; i < (*Counts).Molecules; i++) {
    // second molecule
    for (int j = 0; j < i; j++) {
      int agg_i = (*Molecule)[i].Aggregate;
      int agg_j = (*Molecule)[j].Aggregate;
      // molecules 'i' and 'j' are in contact //{{{
      if (contact[i][j] >= contacts) {
        // create new aggregate if 'j' isn'it in any //{{{
        if (agg_j == -1) {
          agg_j = (*Counts).Aggregates;
          (*Molecule)[j].Aggregate = agg_j;

          (*Aggregate)[agg_j].nMolecules = 1;
          (*Aggregate)[agg_j].Molecule[0] = j;

          (*Counts).Aggregates++;
        } //}}}

        // add 'mol_i' to 'agg_j' aggregate (that contains 'mol_j' molecule) if 'i' isn't in an agg //{{{
        if (agg_i == -1) {
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].Molecule[mols] = i;
          (*Aggregate)[agg_j].nMolecules++;

          (*Molecule)[i].Aggregate = agg_j;
        } //}}}

        // 'mol_i' and 'mol_j' molecules are in different aggregate => unite aggregates
        if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {

          // add molecules from aggregate 'agg_i' to 'agg_j' //{{{
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].nMolecules += (*Aggregate)[agg_i].nMolecules;

          // copy molecule ids from Aggregate[agg_i-1] to Aggregate[agg_j-1]
          int id1 = 0;
          for (int k = mols; k < (*Aggregate)[agg_j].nMolecules; k++) {
            int mol = (*Aggregate)[agg_i].Molecule[id1];
            (*Aggregate)[agg_j].Molecule[k] = mol;
            (*Molecule)[mol].Aggregate = agg_j;
            id1++;
          } //}}}

          // move aggregates with id greater then agg_i to id-1 //{{{
          for (int k = (agg_i+1); k < (*Counts).Aggregates; k++) {

            (*Aggregate)[k-1].nMolecules = (*Aggregate)[k].nMolecules;

            // move every molecule from aggregate 'k' to aggregate 'k-1'
            for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
              int mol = (*Aggregate)[k].Molecule[l];
              (*Aggregate)[k-1].Molecule[l] = mol;
              (*Molecule)[mol].Aggregate = k - 1;
            }
          } //}}}

          // reduce number of aggregates (two aggregates were merged)
          (*Counts).Aggregates--;
        } //}}}
      } else if (agg_j == -1) { // or 'i' and 'j' aren't in contact and 'j' isn't in any aggregate =>  new aggregate for 'j' */ //{{{
        int agg_j = (*Counts).Aggregates;
        (*Molecule)[j].Aggregate = agg_j;

        (*Aggregate)[agg_j].nMolecules = 1;
        (*Aggregate)[agg_j].Molecule[0] = j;

        (*Counts).Aggregates++;
      } //}}}
    }
  }

  // check if highest id residue is in aggregate //{{{
  bool test = false;
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    for (int j = 1; j < (*Aggregate)[i].nMolecules; j++) {
      if ((*Aggregate)[i].Molecule[j] == ((*Counts).Molecules-1)) {
        test = 1;
      }
    }
  } //}}}
  // if highest id residue isn't in any aggregate, create separate one //{{{
  if (!test) {
    int aggs = (*Counts).Aggregates;
    (*Aggregate)[aggs].nMolecules = 1;
    (*Aggregate)[aggs].Molecule[0] = (*Counts).Molecules - 1;

    (*Counts).Aggregates++;
  } //}}}
} //}}}

// SortAggStruct() //{{{
/**
 * Sort an Aggregate struct using the bubble sort algorithm. The resulting
 * struct is arranged so that aggregates with the first molecule's lower id
 * come first.
 */
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts,
                   MOLECULE *Molecule, MOLECULETYPE *MoleculeType,
                   BEAD **Bead, BEADTYPE *BeadType) {
  for (int i = 0; i < (Counts.Aggregates-1); i++) {
    bool done = true;
    for (int j = 0; j < (Counts.Aggregates-i-1); j++) {
      if ((*Aggregate)[j].Molecule[0] > (*Aggregate)[j+1].Molecule[0]) {
        SwapInt(&(*Aggregate)[j].nMolecules, &(*Aggregate)[j+1].nMolecules);
        // switch the whole Aggregate[].Molecule array
        int mols; // number of molecules in the larger aggregate
        if ((*Aggregate)[j].nMolecules > (*Aggregate)[j+1].nMolecules) {
          mols = (*Aggregate)[j].nMolecules;
        } else {
          mols = (*Aggregate)[j+1].nMolecules;
        }
        for (int k = 0; k < mols; k++) {
          SwapInt(&(*Aggregate)[j].Molecule[k], &(*Aggregate)[j+1].Molecule[k]);
        }
        // switch bonded beads array
        SwapInt(&(*Aggregate)[j].nBeads, &(*Aggregate)[j+1].nBeads);
        int beads; // number of beads in the larger aggregate
        if ((*Aggregate)[j].nBeads > (*Aggregate)[j+1].nBeads) {
          beads = (*Aggregate)[j].nBeads;
        } else {
          beads = (*Aggregate)[j+1].nBeads;
        }
        for (int k = 0; k < beads; k++) {
          SwapInt(&(*Aggregate)[j].Bead[k], &(*Aggregate)[j+1].Bead[k]);
        }
        // switch monomer beads array
        SwapInt(&(*Aggregate)[j].nMonomers, &(*Aggregate)[j+1].nMonomers);
        int mons; // larger number of monomers of the two aggregates
        if ((*Aggregate)[j].nMonomers > (*Aggregate)[j+1].nMonomers) {
          mons = (*Aggregate)[j].nMonomers;
        } else {
          mons = (*Aggregate)[j+1].nMonomers;
        }
        for (int k = 0; k < mons; k++) {
          SwapInt(&(*Aggregate)[j].Monomer[k], &(*Aggregate)[j+1].Monomer[k]);
        }
        done = false;
      }
    }
    if (done)
      break;
  }

  // re-assign aggregate id to every bonded bead in the aggregate, correcting after sorting //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      int mtype = Molecule[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
        int id = Molecule[mol].Bead[k];
        (*Bead)[id].nAggregates = 1;
        (*Bead)[id].Aggregate[0] = i;
      }
    }
  } //}}}
} //}}}

// LinkedList() //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead,
                int **Head, int **Link, double cell_size, INTVECTOR *n_cells,
                int *Dcx, int *Dcy, int *Dcz) {

  (*n_cells).x = ceil(BoxLength.x/cell_size);
  (*n_cells).y = ceil(BoxLength.y/cell_size);
  (*n_cells).z = ceil(BoxLength.z/cell_size);

  // allocate arrays
  *Head = malloc(((*n_cells).x*(*n_cells).y*(*n_cells).z)*sizeof(int));
  *Link = malloc(Counts.Beads*sizeof(int));
  for (int i = 0; i < ((*n_cells).x*(*n_cells).y*(*n_cells).z); i++) {
    (*Head)[i] = -1;
  }

  // sort beads into cells //{{{
  for (int i = 0; i < Counts.Beads; i++) {
    int cell = (int)(Bead[i].Position.x / cell_size)
             + (int)(Bead[i].Position.y / cell_size) * (*n_cells).x
             + (int)(Bead[i].Position.z / cell_size) * (*n_cells).x * (*n_cells).y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  } //}}}

  // coordinates of adjoining cells //{{{
  int x[14] = {0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
  int y[14] = {0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
  int z[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 14; i++) {
    Dcx[i] = x[i];
    Dcy[i] = y[i];
    Dcz[i] = z[i];
  } //}}}
} //}}}

// SortBonds() //{{{
/**
 * Function to sort a bond array. As each bond contains a 2-member array, first
 * check that first bead id is lower than the second. Then sort the bonds
 * according to the index of the first bead in each bond.
 */
void SortBonds(int **bond, int length) {
  // first, check order in every bond
  for (int j = 0; j < length; j++) {
    if (bond[j][0] > bond[j][1]) {
      SwapInt(&bond[j][0], &bond[j][1]);
    }
  }
  // second, bubble sort bonds
  for (int j = 0; j < (length-1); j++) {
    bool swap = false;
    for (int k = 0; k < (length-j-1); k++) {
      if ((bond[k][0] > bond[k+1][0]) || // swap if first beads are in wrong order
          (bond[k][0] == bond[k+1][0] && // or if they're the same, but second ones are in wrong order
          bond[k][1] > bond[k+1][1])) {
        SwapInt(&bond[k][0], &bond[k+1][0]);
        SwapInt(&bond[k][1], &bond[k+1][1]);
        SwapInt(&bond[k][2], &bond[k+1][2]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// SortAngles() //{{{
/**
 * Function to sort an angle array. As each angle contains a 3-member array
 * with the middle number being the 'centre' of the angle (i.e., it must remain
 * in the middle). Sort it so that the first index is lower than the third one
 * and then ascendingly according to the first indices.
 */
void SortAngles(int **angle, int length) {
  // first, check order of the 1st and 3rd id in every angle
  for (int j = 0; j < length; j++) {
    if (angle[j][0] > angle[j][2]) {
      SwapInt(&angle[j][0], &angle[j][2]);
    }
  }
  // second, bubble sort angles
  for (int j = 0; j < (length-1); j++) {
    bool swap = false;
    for (int k = 0; k < (length-j-1); k++) {
      if ((angle[k][0] > angle[k+1][0]) || // swap if first beads are in wrong order
          (angle[k][0] == angle[k+1][0] &&
           angle[k][2] > angle[k+1][2])) { // or if they're the same, but 3rd ones are in wrong order
        SwapInt(&angle[k][0], &angle[k+1][0]);
        SwapInt(&angle[k][1], &angle[k+1][1]);
        SwapInt(&angle[k][2], &angle[k+1][2]);
        SwapInt(&angle[k][3], &angle[k+1][3]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(COUNTS Counts, BEAD **Bead) {
  for (int i = 0; i < Counts.BeadsInVsf; i++) {
    free((*Bead)[i].Aggregate);
  }
  free(*Bead);
} //}}}

// FreeMolecule() //{{{
/**
 * Free memory allocated for Molecule struct array. This function makes it
 * easier other arrays to the Molecule struct in the future
 */
void FreeMolecule(COUNTS Counts, MOLECULE **Molecule) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Molecule)[i].Bead);
  }
  free(*Molecule);
} //}}}

// FreeMoleculeType() //{{{
/**
 * Free memory allocated for MoleculeType struct array. This function makes
 * it easier other arrays to the MoleculeType struct in the future
 */
void FreeMoleculeType(COUNTS Counts, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      free((*MoleculeType)[i].Bond[j]);
    }
    for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
      free((*MoleculeType)[i].Angle[j]);
    }
    free((*MoleculeType)[i].Bead);
    free((*MoleculeType)[i].Bond);
    free((*MoleculeType)[i].Angle);
    free((*MoleculeType)[i].BType);
  }
  free(*MoleculeType);
} //}}}

// FreeAggregate() //{{{
/**
 * Free memory allocated for Aggregate struct array. This function makes it
 * easier other arrays to the Aggregate struct in the future
 */
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Aggregate)[i].Molecule);
    free((*Aggregate)[i].Bead);
    free((*Aggregate)[i].Monomer);
  }
  free(*Aggregate);
} //}}}
