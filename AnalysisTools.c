#include "AnalysisTools.h"
#include "Errors.h"
#include <string.h>

// TODO: consider BeadType[].Index, System.Bonded, etc. arrays - shouldn't they
//       be filled based on whether the beads are in the timestep? Plus a
//       function that says 'everything's there' for utilities not using
//       coordinates (e.g., DistrAgg). Needs adding IndexCoor count, I guess, so
//       we know how many are in that timestep.
//       ...would speed some cases, and it's easier to use since it does not
//       have to be combined with Bead[].InTimestep anymore.
//       TODO: Done already? ...via System.{Bead,Bonded,Unbonded}Coor arrays

// TODO: consider the names - can't I have different bead types with the same
//       name? ...I'd only need to ensure that going over named beads would
//       encompass all bead types with that name.

// TODO: check Molecule[].InTimestep stuff - is it implemented for all input
//       files?

// TODO: System.Index_mol[] array - what's it for? Useless, right? But it needs
//       something like MolCoor array akin to BeadCoor to specify which
//       molecules are in a timestep. And corresponding Count.MolCoor, I guess.
//       Oh, there's already Count.MoleculeCoor defined - but not used

// STATIC DEFINITIONS
// functions to transform to/from fractional coordinates
static void ToFractional(double coor[3], BOX Box);
// static void ToFractionalCoor(SYSTEM *System);
static void FromFractional(double coor[3], BOX Box);
// static void FromFractionalCoor(SYSTEM *System);
// test whether all bonds in a molecule are connected
static bool ConnectedMolecule(SYSTEM System, int n);
// remove pbc for molecules by joining the molecules
static void RemovePBCMolecules(SYSTEM *System);
// restore pbc by wrapping all coordinates inside the simulation box
static void RestorePBC(SYSTEM *System);

// STATIC IMPLEMENTATIONS
// transform to/from fractional coordinates
/* HOW TO CALCULATE DISTANCE IN TRICLINIC SYSTEM //{{{
//double dist[3];
//dist[0] = (*Bead)[0].Position[0] - (*Bead)[10].Position[0];
//dist[1] = (*Bead)[0].Position[1] - (*Bead)[10].Position[1];
//dist[2] = (*Bead)[0].Position[2] - (*Bead)[10].Position[2];
//printf("dist1 = (%lf, %lf, %lf) = %lf\n", dist[0], dist[1], dist[2],
sqrt(SQR(dist[0])+SQR(dist[1])+SQR(dist[2])));

//double new[3];
//new[0] = Box.transform[0][0] * dist[0] +
//         Box.transform[0][1] * dist[1] +
//         Box.transform[0][2] * dist[2];
//new[1] = Box.transform[1][0] * dist[0] +
//         Box.transform[1][1] * dist[1] +
//         Box.transform[1][2] * dist[2];
//new[2] = Box.transform[2][0] * dist[0] +
//         Box.transform[2][1] * dist[1] +
//         Box.transform[2][2] * dist[2];
//dist[0] = new[0] / a;
//dist[1] = new[1] / b;
//dist[2] = new[2] / c;
//printf("dist2 = (%lf, %lf, %lf) = %lf\n", dist[0], dist[1], dist[2],
sqrt(SQR(dist[0])+SQR(dist[1])+SQR(dist[2])));
*/ //}}}
static void ToFractional(double coor[3], BOX Box) { //{{{
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    double new[3] = {0, 0, 0};
    for (int dd = 0; dd < 3; dd++) {
      new[dd] = Box.inverse[dd][0] * coor[0] +
                Box.inverse[dd][1] * coor[1] +
                Box.inverse[dd][2] * coor[2];
    }
    for (int dd = 0; dd < 3; dd++) {
      coor[dd] = new[dd] * Box.Length[dd];
    }
  }
} //}}}
void ToFractionalCoor(SYSTEM *System) { //{{{
  if (System->Box.alpha != 90 ||
      System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      BEAD *bead = &System->Bead[id];
      ToFractional(bead->Position, System->Box);
    }
  }
} //}}}
static void FromFractional(double coor[3], BOX Box) { //{{{
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    double new[3] = {0, 0, 0};
    for (int dd = 0; dd < 3; dd++) {
      new[dd] = Box.transform[dd][0] * coor[0] / Box.Length[dd] +
                Box.transform[dd][1] * coor[1] / Box.Length[dd] +
                Box.transform[dd][2] * coor[2] / Box.Length[dd];
    }
    for (int dd = 0; dd < 3; dd++) {
      coor[dd] = new[dd];
    }
  }
} //}}}
void FromFractionalCoor(SYSTEM *System) { //{{{
  if (System->Box.alpha != 90 ||
      System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      BEAD *bead = &System->Bead[id];
      FromFractional(bead->Position, System->Box);
    }
  }
} //}}}
// test whether all bonds in a molecule are connected //{{{
static bool ConnectedMolecule(SYSTEM System, int n) {
  MOLECULE *mol = &System.Molecule[n];
  MOLECULETYPE *mt = &System.MoleculeType[mol->Type];
  if (mt->nBonds == 0) {
    return false;
  }
  int *connected = calloc(mt->nBonds, sizeof *connected),
      *unconnected = calloc(mt->nBonds, sizeof *unconnected),
      count_connected = 0, count_unconnected = 0;
  for (int i = 0; i < mt->nBonds; i++) {
    int id1 = mt->Bond[i][0], id2 = mt->Bond[i][1];
    BEAD *b_1 = &System.Bead[mol->Bead[id1]],
         *b_2 = &System.Bead[mol->Bead[id2]];
    if (b_1->InTimestep && b_2->InTimestep) {
      if (count_connected == 0) {
        connected[count_connected] = i;
        count_connected++;
      } else {
        unconnected[count_unconnected] = i;
        count_unconnected++;
      }
    }
  }
  if (count_unconnected == 0 && count_connected == 0) {
    // TODO: this means the molecule isn't in the step - change into int ouptut?
    free(connected);
    free(unconnected);
    return true;
  }

  for (int i = 0; i < count_connected; i++) {
    for (int j = 0; j < count_unconnected; j++) {
      int bond[2], con[2], uncon[2];
      bond[0] = connected[i];
      bond[1] = unconnected[j];
      con[0] = mt->Bond[bond[0]][0];
      con[1] = mt->Bond[bond[0]][1];
      uncon[0] = mt->Bond[bond[1]][0];
      uncon[1] = mt->Bond[bond[1]][1];
      if (con[0] == uncon[0] || con[0] == uncon[1] ||
          con[1] == uncon[0] || con[1] == uncon[1]) {
        connected[count_connected] = bond[1];
        count_connected++;
        count_unconnected--;
        for (int k = j; k < count_unconnected; k++) {
          unconnected[k] = unconnected[k+1];
        }
        j--;
      }
    }
  }
  if (count_unconnected > 0) {
    return false;
  }
  free(connected);
  free(unconnected);
  return true;
} //}}}
// remove pbc for molecules by joining the molecules //{{{
/*
 * Create a list of all bonds ('unconnected' array) with beads that are in the
 * timestep. Then create a connectivity array by going through the list,
 * transferring used bonds into 'connected' array, and finally, use the
 * 'connected' array to join the molecule. As long as there are bonds in the
 * 'unconnected' array, continue creating a new 'connected' array and joining
 * the molecule. This procedure should be able to join molecule of any
 * complexity as well as molecule where some beads ar not connected.
 */
/*
 * TODO: maybe split molecule to connected pieces, connect those and place
 *       their centres of mass nearest each other
 */
static void RemovePBCMolecules(SYSTEM *System) {
  BOX *box = &System->Box;
  // go through all molecules
  for (int mm = 0; mm < System->Count.Molecule; mm++) {
    MOLECULE *mol = &System->Molecule[mm];
    if (!mol->InTimestep) {
      continue;
    }
    MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
    // skip molecule if it is bond-less
    if (mt->nBonds == 0) {
      snprintf(ERROR_MSG, LINE, "molecule %s%d%s (%s%s%s) has no bonds ",
               ErrYellow(), mol->Index, ErrCyan(),
               ErrYellow(), mt->Name, ErrCyan());
      PrintWarning();
      continue;
    }
    // arrays holding bonds already connected and yet unconnected
    int *connected = calloc(mt->nBonds, sizeof *connected),
        *unconnected = calloc(mt->nBonds, sizeof *unconnected),
        count_unconnected = 0;
    // 1)
    for (int i = 0; i < mt->nBonds; i++) {
      int id[2] = {mol->Bead[mt->Bond[i][0]], mol->Bead[mt->Bond[i][1]]};
      BEAD *b_1 = &System->Bead[id[0]],
           *b_2 = &System->Bead[id[1]];
      // printf("%d %d\n", id[0], id[1]);
      if (b_1->InTimestep && b_2->InTimestep) {
        unconnected[count_unconnected] = i;
        count_unconnected++;
      }
    }
    // skip molecule if there is no valid bond in the coordinate file //{{{
    if (count_unconnected == 0) {
      snprintf(ERROR_MSG, LINE, "no bonded beads in the timestep "
               " for molecule %s%d%s (%s%s%s) has no bonds ", ErrYellow(),
               mol->Index, ErrCyan(), ErrYellow(), mt->Name, ErrCyan());
      PrintWarning();
      free(connected);
      free(unconnected);
      continue;
    } //}}}
    while (count_unconnected > 0) {
      int count_connected = 0;
      connected[count_connected] = unconnected[0];
      count_connected++;
      count_unconnected--;
      for (int i = 0; i < count_unconnected; i++) {
        unconnected[i] = unconnected[i+1];
      }
      // 2)
      for (int i = 0; i < count_connected; i++) {
        for (int j = 0; j < count_unconnected; j++) {
          int bond[2], // the connected and unconnected bonds
              con[2], // beads in the already connected bond (bond[0])
              uncon[2]; // beads in the yet unconneced bond (bond[1])
          bond[0] = connected[i];
          bond[1] = unconnected[j];
          con[0] = mt->Bond[bond[0]][0];
          con[1] = mt->Bond[bond[0]][1];
          uncon[0] = mt->Bond[bond[1]][0];
          uncon[1] = mt->Bond[bond[1]][1];
          // if a bead is in both bonds, the unconnected bond becomes connected
          if (con[0] == uncon[0] || con[0] == uncon[1] ||
              con[1] == uncon[0] || con[1] == uncon[1]) {
            connected[count_connected] = bond[1];
            count_connected++;
            count_unconnected--;
            // move unconnected bonds to retain continuous array
            for (int k = j; k < count_unconnected; k++) {
              unconnected[k] = unconnected[k+1];
            }
            // unconnected[j] is again unconnected, so decremenet 'j'
            j--;
          }
        }
      }
      // connect the molecule by going through the list of connected bonds 
      bool *moved = calloc(mt->nBeads, sizeof *moved);
      int first = mt->Bond[connected[0]][0];
      moved[first] = true;
      for (int i = 0; i < count_connected; i++) {
        int bond = connected[i],
            id[2] = {mt->Bond[bond][0], mt->Bond[bond][1]};
        BEAD *b_1 = &System->Bead[mol->Bead[id[0]]],
             *b_2 = &System->Bead[mol->Bead[id[1]]];
        double dist[3];
        if (!moved[id[0]] && moved[id[1]]) {
          Distance(b_2->Position, b_1->Position, box->OrthoLength, dist);
          for (int dd = 0; dd < 3; dd++) {
            b_1->Position[dd] = b_2->Position[dd] - dist[dd];
          }
          moved[id[0]] = true;
        } else if (moved[id[0]] && !moved[id[1]]) {
          Distance(b_1->Position, b_2->Position, box->OrthoLength, dist);
          for (int dd = 0; dd < 3; dd++) {
            b_2->Position[dd] = b_1->Position[dd] - dist[dd];
          }
          moved[id[1]] = true;
        }
      }
      // TODO: CENTRE OF MASS
      free(moved);
    }
    free(connected);
    free(unconnected);
    // put molecule's geometric centre into the simulation box //{{{
    double cog[3];
    GeomCentre(mt->nBeads, mol->Bead, System->Bead, cog);
    // by how many BoxLength's should cog be moved?
    int move[3];
    for (int dd = 0; dd < 3; dd++) {
      move[dd] = cog[dd] / box->OrthoLength[dd];
      if (cog[dd] < 0) {
        move[dd]--;
      }
    }
    for (int j = 0; j < mt->nBeads; j++) {
      int bead = mol->Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[bead].Position[dd] -= move[dd] * box->OrthoLength[dd];
      }
    } //}}}
  }
} //}}}
// remove pbc for molecules by joining the molecules //{{{
/*
 * Create a list of all bonds ('unconnected' array) with beads that are in the
 * timestep. Then create a connectivity array by going through the list,
 * transferring used bonds into 'connected' array, and finally, use the
 * 'connected' array to join the molecule. As long as there are bonds in the
 * 'unconnected' array, continue creating a new 'connected' array and joining
 * the molecule. This procedure should be able to join molecule of any
 * complexity as well as molecule where some beads ar not connected.
 */
/*
 * TODO: maybe split molecule to connected pieces, connect those and place
 *       their centres of mass nearest each other
 */
static void RemovePBCMolecules_old(SYSTEM *System) {
  BOX *box = &System->Box;
  // go through all molecules
  for (int mm = 0; mm < System->Count.Molecule; mm++) {
    MOLECULE *mol = &System->Molecule[mm];
    if (!mol->InTimestep) {
      continue;
    }
    MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
    // skip molecule if it is bond-less
    if (mt->nBonds == 0) {
      snprintf(ERROR_MSG, LINE, "molecule %s%d%s (%s%s%s) has no bonds ",
               ErrYellow(), mol->Index, ErrCyan(),
               ErrYellow(), mt->Name, ErrCyan());
      PrintWarning();
      continue;
    }
    // arrays holding bonds already connected and yet unconnected
    int *connected = calloc(mt->nBonds, sizeof *connected),
        *unconnected = calloc(mt->nBonds, sizeof *unconnected),
        count_unconnected = 0;
    // 1)
    for (int i = 0; i < mt->nBonds; i++) {
      int id[2] = {mol->Bead[mt->Bond[i][0]], mol->Bead[mt->Bond[i][1]]};
      BEAD *b_1 = &System->Bead[id[0]],
           *b_2 = &System->Bead[id[1]];
      // printf("%d %d\n", id[0], id[1]);
      if (b_1->InTimestep && b_2->InTimestep) {
        unconnected[count_unconnected] = i;
        count_unconnected++;
      }
    }
    // skip molecule if there is no valid bond in the coordinate file //{{{
    if (count_unconnected == 0) {
      snprintf(ERROR_MSG, LINE, "no bonded beads in the timestep "
               " for molecule %s%d%s (%s%s%s) has no bonds ", ErrYellow(),
               mol->Index, ErrCyan(), ErrYellow(), mt->Name, ErrCyan());
      PrintWarning();
      free(connected);
      free(unconnected);
      continue;
    } //}}}
    while (count_unconnected > 0) {
      int count_connected = 0;
      connected[count_connected] = unconnected[0];
      count_connected++;
      count_unconnected--;
      for (int i = 0; i < count_unconnected; i++) {
        unconnected[i] = unconnected[i+1];
      }
      // 2)
      for (int i = 0; i < count_connected; i++) {
        for (int j = 0; j < count_unconnected; j++) {
          int bond[2], // the connected and unconnected bonds
              con[2], // beads in the already connected bond (bond[0])
              uncon[2]; // beads in the yet unconneced bond (bond[1])
          bond[0] = connected[i];
          bond[1] = unconnected[j];
          con[0] = mt->Bond[bond[0]][0];
          con[1] = mt->Bond[bond[0]][1];
          uncon[0] = mt->Bond[bond[1]][0];
          uncon[1] = mt->Bond[bond[1]][1];
          // if a bead is in both bonds, the unconnected bond becomes connected
          if (con[0] == uncon[0] || con[0] == uncon[1] ||
              con[1] == uncon[0] || con[1] == uncon[1]) {
            connected[count_connected] = bond[1];
            count_connected++;
            count_unconnected--;
            // move unconnected bonds to retain continuous array
            for (int k = j; k < count_unconnected; k++) {
              unconnected[k] = unconnected[k+1];
            }
            // unconnected[j] is again unconnected, so decremenet 'j'
            j--;
          }
        }
      }
      // connect the molecule by going through the list of connected bonds 
      bool *moved = calloc(mt->nBeads, sizeof *moved);
      int first = mt->Bond[connected[0]][0];
      moved[first] = true;
      for (int i = 0; i < count_connected; i++) {
        int bond = connected[i],
            id[2] = {mt->Bond[bond][0], mt->Bond[bond][1]};
        BEAD *b_1 = &System->Bead[mol->Bead[id[0]]],
             *b_2 = &System->Bead[mol->Bead[id[1]]];
        double dist[3];
        if (!moved[id[0]] && moved[id[1]]) {
          Distance(b_2->Position, b_1->Position, box->OrthoLength, dist);
          for (int dd = 0; dd < 3; dd++) {
            b_1->Position[dd] = b_2->Position[dd] - dist[dd];
          }
          moved[id[0]] = true;
        } else if (moved[id[0]] && !moved[id[1]]) {
          Distance(b_1->Position, b_2->Position, box->OrthoLength, dist);
          for (int dd = 0; dd < 3; dd++) {
            b_2->Position[dd] = b_1->Position[dd] - dist[dd];
          }
          moved[id[1]] = true;
        }
      }
      // TODO: CENTRE OF MASS
      free(moved);
    }
    free(connected);
    free(unconnected);
    // put molecule's geometric centre into the simulation box //{{{
    double cog[3];
    GeomCentre(mt->nBeads, mol->Bead, System->Bead, cog);
    // by how many BoxLength's should cog be moved?
    int move[3];
    for (int dd = 0; dd < 3; dd++) {
      move[dd] = cog[dd] / box->OrthoLength[dd];
      if (cog[dd] < 0) {
        move[dd]--;
      }
    }
    for (int j = 0; j < mt->nBeads; j++) {
      int bead = mol->Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[bead].Position[dd] -= move[dd] * box->OrthoLength[dd];
      }
    } //}}}
  }
} //}}}
// restore pbc by wrapping all coordinates inside the simulation box //{{{
static void RestorePBC(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    int id = System->BeadCoor[i];
    BEAD *bead = &System->Bead[id];
    BOX *box = &System->Box;
    for (int dd = 0; dd < 3; dd++) {
      while (bead->Position[dd] >= box->OrthoLength[dd]) {
        bead->Position[dd] -= box->OrthoLength[dd];
      }
      while (bead->Position[dd] < 0) {
        bead->Position[dd] += box->OrthoLength[dd];
      }
    }
  }
} //}}}

// Helper functions for dealing with SYSTEM structure
// get bead indices from bonds/angles/dihedrals/impropers //{{{
/*
 * They all return an array of size 2*n where n is the number of bead indices
 * (2 for bonds, 3 for angles, 4 for dihedrals and impropers); elements 0..n-1
 * contain internal bead indices for the corresponding array in MoleculeType
 * (Bond, Angle, etc.), while elements n..2*n-1 contain bead indices for in
 * Bead.
 */
extern inline int *BondIndices(SYSTEM System, int mol, int bond) { //{{{
  static int index[4]; // first are 'real' ids, then the intramolecular ones
  int n = 2, type = System.Molecule[mol].Type;
  // intramolecular id
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Bond[bond][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in a bond; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in bond %s%d%s:%s", ErrYellow(),
            bond, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in a bond; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in bond %s%d%s:%s", ErrYellow(), bond, ErrRed(),
            ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *AngleIndices(SYSTEM System, int mol, int angle) { //{{{
  static int index[6]; // first are 'real' ids, then the intramolecular ones
  int n = 3, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Angle[angle][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in an angle; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in angle %s%d%s:%s", ErrYellow(),
            angle, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in an angle; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in angle %s%d%s:%s", ErrYellow(), angle, ErrRed(),
            ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *DihedralIndices(SYSTEM System, int mol, int dihed) { //{{{
  static int index[8]; // first are 'real' ids, then the intramolecular ones
  int n = 4, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Dihedral[dihed][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in an dihedral; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in dihedral %s%d%s:%s",
            ErrYellow(), dihed, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in a dihedral; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in dihedral %s%d%s:%s", ErrYellow(), dihed,
            ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *ImproperIndices(SYSTEM System, int mol, int improper) { //{{{
  static int index[8]; // first are 'real' ids, then the intramolecular ones
  int n = 4, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Improper[improper][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in an improper; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in improper %s%d%s:%s",
            ErrYellow(), improper, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    err_msg("wrong index in a improper; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in improper %s%d%s:%s", ErrYellow(), improper,
            ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
  //}}}
// identify bead type based on name //{{{
int FindBeadType(char name[], SYSTEM System) {
  int type;
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (strcmp(name, System.BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in BeadType struct
  return -1;
} //}}}
// identify molecule type based on name //{{{
int FindMoleculeName(char name[], SYSTEM System) {
  int type = -1;
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    if (strcmp(name, System.MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return type;
} //}}}
// identify molecule type based on name only or on other parameters too //{{{
/*
 * If name==true, then find molecule type according to mode:
 *   0 - check nothing else
 *   1 - check number of beads
 *   2 - check number and types of beads
 *   3 - check everything
 * name = true/false for checking/ignoring molecule name
 */
int FindMoleculeType(SYSTEM Sys1, MOLECULETYPE mt, SYSTEM Sys2, int mode) {
  // just to be sure the function's mode parameter is correct //{{{
  if (mode < 0 || mode > 3) {
    err_msg("FindMoleculeType() - mode parameter must be <0,3>\n");
    PrintError();
    exit(1);
  } //}}}
  MOLECULETYPE *mt_1 = &mt;
  // find the same name
  for (int i = 0; i < Sys2.Count.MoleculeType; i++) {
    MOLECULETYPE *mt_2 = &Sys2.MoleculeType[i];
    if (strcmp(mt_1->Name, mt_2->Name) == 0) {
      if (mode == 0) { // only name checked
        return i;
      }
      // check number of beads //{{{
      if (mt_1->nBeads != mt_2->nBeads) {
        goto end_loop; // not 'continue;' to avoid break/continue combo later
      } else if (mode == 1) {
        return i;
      } //}}}
      // check bead order //{{{
      for (int j = 0; j < mt_1->nBeads; j++) {
        int bt_1 = mt_1->Bead[j],
            bt_2 = mt_2->Bead[j];
        if (!SameBeadType(Sys1.BeadType[bt_1], Sys2.BeadType[bt_2])) {
          goto end_loop;
        }
      }
      if (mode == 2) {
        return i;
      } //}}}
      // check bonds //{{{
      if (mt_1->nBonds != mt_2->nBonds) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nBonds; j++) {
        if (!SameArrayInt(mt_1->Bond[j], mt_2->Bond[j], 2)) {
          goto end_loop;
        }
        if (mt_1->Bond[j][2] == mt_2->Bond[j][2]) {
          if (mt_1->Bond[j][2] != -1) {
            PARAMS *tbond_1 = &Sys1.BondType[mt_1->Bond[j][2]],
                   *tbond_2 = &Sys2.BondType[mt_2->Bond[j][2]];
            if (tbond_1->a != tbond_2->a ||
                tbond_1->b != tbond_2->b ||
                tbond_1->c != tbond_2->c ||
                tbond_1->d != tbond_2->d) {
              goto end_loop;
            }
          }
        } else {
          goto end_loop;
        }
      } //}}}
      // check angles //{{{
      if (mt_1->nAngles != mt_2->nAngles) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nAngles; j++) {
        if (!SameArrayInt(mt_1->Angle[j], mt_2->Angle[j], 3)) {
          goto end_loop;
        }
        if (mt_1->Angle[j][3] != -1 && mt_2->Angle[j][3] != -1) {
          PARAMS *tangle_1 = &Sys1.AngleType[mt_1->Angle[j][3]],
                 *tangle_2 = &Sys2.AngleType[mt_2->Angle[j][3]];
          if (tangle_1->a != tangle_2->a ||
              tangle_1->b != tangle_2->b ||
              tangle_1->c != tangle_2->c ||
              tangle_1->d != tangle_2->d) {
            goto end_loop;
          }
        }
      } //}}}
      // check dihedrals //{{{
      if (mt_1->nDihedrals != mt_2->nDihedrals) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nDihedrals; j++) {
        if (!SameArrayInt(mt_1->Dihedral[j], mt_2->Dihedral[j], 4)) {
          goto end_loop;
        }
        if (mt_1->Dihedral[j][4] != -1 && mt_2->Dihedral[j][4] != -1) {
          PARAMS *tdihed_1 = &Sys1.DihedralType[mt_1->Dihedral[j][4]],
                 *tdihed_2 = &Sys2.DihedralType[mt_2->Dihedral[j][4]];
          if (tdihed_1->a != tdihed_2->a ||
              tdihed_1->b != tdihed_2->b ||
              tdihed_1->c != tdihed_2->c ||
              tdihed_1->d != tdihed_2->d) {
            goto end_loop;
          }
        }
      } //}}}
      // check impropers //{{{
      if (mt_1->nImpropers != mt_2->nImpropers) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nImpropers; j++) {
        if (!SameArrayInt(mt_1->Improper[j], mt_2->Improper[j], 4)) {
          goto end_loop;
        }
        if (mt_1->Improper[j][4] != -1 && mt_1->Improper[j][4] != -1) {
          PARAMS *timpro_1 = &Sys1.ImproperType[mt_1->Improper[j][4]],
                 *timpro_2 = &Sys2.ImproperType[mt_2->Improper[j][4]];
          if (timpro_1->a != timpro_2->a ||
              timpro_1->b != timpro_2->b ||
              timpro_1->c != timpro_2->c ||
              timpro_1->d != timpro_2->d) {
            goto end_loop;
          }
        }
      }         //}}}
      return i; // assumes mode=3, obviously
    }
    end_loop:;
  }
  return -1;
} //}}}

// Helper functions for manipulating coordinates
// wrap coordinates into simulation box and/or join molecules //{{{
void WrapJoinCoordinates(SYSTEM *System, bool wrap, bool join) {
  if (System->Box.Volume != -1 && (wrap || join)) {
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(System);
    if (wrap) { // wrap coordinates into the simulation box
      RestorePBC(System);
    }
    if (join) { // join molecules by removing periodic boundary conditions
      RemovePBCMolecules(System);
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FromFractionalCoor(System);
  }
} //}}}
// distance between two beads; in the range <-BoxLength/2,BoxLength/2) //{{{
void Distance(const double id1[3], const double id2[3],
              const double BoxLength[3], double out[3]) {
  // remove periodic boundary conditions in x-direction
  for (int dd = 0; dd < 3; dd++) {
    out[dd] = id1[dd] - id2[dd];
    while (out[dd] >= (BoxLength[dd] / 2)) {
      out[dd] -= BoxLength[dd];
    }
    while (out[dd] < (-BoxLength[dd] / 2)) {
      out[dd] += BoxLength[dd];
    }
  }
} //}}}
// calculate centre of mass for a list of beads //{{{
void CentreOfMass(int n, const int list[], SYSTEM System, double com[3]) {
  for (int dd = 0; dd < 3; dd++) {
    com[dd] = 0;
  }
  double mass = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    BEAD *b = &System.Bead[id];
    BEADTYPE *bt = &System.BeadType[b->Type];
    for (int dd = 0; dd < 3; dd++) {
      com[dd] += b->Position[dd] * bt->Mass;
    }
    mass += bt->Mass;
  }
  for (int dd = 0; dd < 3; dd++) {
    com[dd] /= mass;
  }
} //}}}
// calculate geometric centre for a list of beads //{{{
void GeomCentre(int n, const int *list, BEAD *Bead, double gc[3]) {
  InitDoubleArray(gc, 3, 0);
  int count = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    if (Bead[id].InTimestep) {
      for (int dd = 0; dd < 3; dd++) {
        gc[dd] += Bead[id].Position[dd];
      }
      count++;
    }
  }
  for (int dd = 0; dd < 3; dd++) {
    gc[dd] /= count;
  }
} //}}}

// identify input coordinate and structure files //{{{
bool InputCoorStruct(int argc, char *argv[], SYS_FILES *f) {
  // input structure file (-i option)
  if (FileOption(argc, argv, "-i", f->stru.name)) {
    f->stru.type = StructureFileType(f->stru.name);
  }
  f->coor.type = CoordinateFileType(f->coor.name);
  // set default structure file if -i option not used
  if (f->stru.name[0] == '\0') {
    if (f->coor.type == VCF_FILE) { // use vcf file with .vsf ending
      int last = -1;
      for (int i = 0; i < strnlen(f->coor.name, LINE); i++) {
        if (f->coor.name[i] == '.') {
          last = i;
        }
      }
      s_strcpy(f->stru.name, f->coor.name, LINE);
      f->stru.name[last+2] = 's';
      f->stru.type = VSF_FILE;
    } else if (f->coor.type == VTF_FILE ||   //
               f->coor.type == XYZ_FILE ||   // use both as a coordinate and
               f->coor.type == LDATA_FILE || // a structure files
               f->coor.type == LTRJ_FILE) {  //
      s_strcpy(f->stru.name, f->coor.name, LINE);
      f->stru.type = f->coor.type;
    } else {
      err_msg("missing structure file; should never happen!");
      PrintError();
      exit(1);
    }
  }
  return true;
} //}}}
int StructureFileType(char name[]) { //{{{
  // a) check for FIELD file
  if (strcasecmp(name, "FIELD") == 0) {
    return FIELD_FILE;
  }
  // b) check for extension
  // copy the name as it's destroyed by strrchr()
  char orig[LINE];
  snprintf(orig, LINE, "%s", name);
  // check for known extensions
  int ext = 6;
  char extension[ext][EXTENSION];
  s_strcpy(extension[0], ".vsf", EXTENSION);
  s_strcpy(extension[1], ".vtf", EXTENSION);
  s_strcpy(extension[2], ".xyz", EXTENSION);
  s_strcpy(extension[3], ".lammpstrj", EXTENSION);
  s_strcpy(extension[4], ".data", EXTENSION);
  s_strcpy(extension[5], ".field", EXTENSION);
  char *dot = strrchr(name, '.');
  for (int i = 0; i < ext; i++) {
    if (dot && strcasecmp(dot, extension[i]) == 0) {
      ext = i;
      break;
    }
  }
  switch (ext) {
    case 0:
      return VSF_FILE;
    case 1:
      return VTF_FILE;
    case 2:
      return XYZ_FILE;
    case 3:
      return LTRJ_FILE;
    case 4:
      return LDATA_FILE;
    case 5:
      return FIELD_FILE;
    default:
      err_msg("Unknown structure file type");
      PrintErrorFile(name, "\0", "\0");
      exit(1);
  }
} //}}}
int CoordinateFileType(char name[]) { //{{{
  // copy the name as it's destroyed by strrchr()
  char orig[LINE];
  snprintf(orig, LINE, "%s", name);
  // check for known extensions
  int ext = 5;
  char extension[ext][EXTENSION];
  s_strcpy(extension[0], ".vcf", EXTENSION);
  s_strcpy(extension[1], ".vtf", EXTENSION);
  s_strcpy(extension[2], ".xyz", EXTENSION);
  s_strcpy(extension[3], ".lammpstrj", EXTENSION);
  s_strcpy(extension[4], ".data", EXTENSION);
  char *dot = strrchr(name, '.');
  for (int i = 0; i < ext; i++) {
    if (dot && strcasecmp(dot, extension[i]) == 0) {
      ext = i;
      break;
    }
  }
  switch (ext) {
    case 0:
      return VCF_FILE;
    case 1:
      return VTF_FILE;
    case 2:
      return XYZ_FILE;
    case 3:
      return LTRJ_FILE;
    case 4:
      return LDATA_FILE;
    default:
      err_msg("Unknown coordinate file type");
      PrintErrorFile(name, "\0", "\0");
      exit(1);
  }
} //}}}
int FileType(char name[]) { //{{{
  // a) check for FIELD/CONFIG file
  if (strcasecmp(name, "FIELD") == 0) {
    return FIELD_FILE;
  } else if (strcasecmp(name, "CONFIG") == 0) {
    return CONFIG_FILE;
  }
  // b) check for extension
  // copy the name as it's destroyed by strrchr()
  char orig[LINE];
  snprintf(orig, LINE, "%s", name);
  // check for known extensions
  int ext = 8;
  char extension[ext][EXTENSION];
  s_strcpy(extension[0], ".vtf", EXTENSION);
  s_strcpy(extension[1], ".vsf", EXTENSION);
  s_strcpy(extension[2], ".vcf", EXTENSION);
  s_strcpy(extension[3], ".xyz", EXTENSION);
  s_strcpy(extension[4], ".data", EXTENSION);
  s_strcpy(extension[5], ".lammpstrj", EXTENSION);
  s_strcpy(extension[6], ".field", EXTENSION);
  s_strcpy(extension[7], ".config", EXTENSION);
  char *dot = strrchr(name, '.');
  for (int i = 0; i < ext; i++) {
    if (dot && strcasecmp(dot, extension[i]) == 0) {
      ext = i;
      break;
    }
  }
  switch (ext) {
    case 0:
      return VTF_FILE;
    case 1:
      return VSF_FILE;
    case 2:
      return VCF_FILE;
    case 3:
      return XYZ_FILE;
    case 4:
      return LDATA_FILE;
    case 5:
      return LTRJ_FILE;
    case 6:
      return FIELD_FILE;
    case 7:
      return CONFIG_FILE;
    default:
      err_msg("Unknown coordinate file type");
      PrintErrorFile(name, "\0", "\0");
      exit(1);
  }
} //}}}

// create a cell-linked list //{{{
void LinkedList(SYSTEM System, int **Head, int **Link, double rcut,
                int n_cells[3], int Dc[14][3]) {
  double (*box)[3] = &System.Box.Length;
  COUNT *Count = &System.Count;
  double rl[3];
  for (int dd = 0; dd < 3; dd++) {
    rl[dd] = (*box)[dd] / rcut;
    n_cells[dd] = (int)(rl[dd]);
  }
  if (n_cells[0] < 3 || n_cells[1] < 3 || n_cells[2] < 3) {
    err_msg("cell size too small for cut-off in linked list");
    PrintError();
    exit(1);
  }
  for (int dd = 0; dd < 3; dd++) {
    rl[dd] = (double)n_cells[dd] / (*box)[dd];
  }
  // allocate arrays
  *Head = malloc(sizeof **Head * n_cells[0] * n_cells[1] * n_cells[2]);
  *Link = malloc(sizeof **Link * Count->BeadCoor);
  for (int i = 0; i < (n_cells[0] * n_cells[1] * n_cells[2]); i++) {
    (*Head)[i] = -1;
  }
  // sort beads into cells
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    int cell = (int)(bead->Position[0] * rl[0]) +
               (int)(bead->Position[1] * rl[1]) * n_cells[0] +
               (int)(bead->Position[2] * rl[2]) * n_cells[0] * n_cells[1];
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  }
  // coordinates of adjoining cells
  int x[14] = {0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};
  int y[14] = {0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
  int z[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 14; i++) {
    Dc[i][0] = x[i];
    Dc[i][1] = y[i];
    Dc[i][2] = z[i];
  }
}
int SelectCell1(const int c1[3], const int n_cells[3]) {
  return c1[0] + c1[1] * n_cells[0] + c1[2] * n_cells[0] * n_cells[1];
}
int SelectCell2(const int c1[3], const int n_cells[3],
                const int Dc[14][3], int n) {
  int c2[3];
  for (int dd = 0; dd < 3; dd++) {
    c2[dd] = c1[dd] + Dc[n][dd];
  }
  // periodic boundary conditions for cells
  if (c2[0] >= n_cells[0])
    c2[0] -= n_cells[0];
  else if (c2[0] < 0)
    c2[0] += n_cells[0];

  if (c2[1] >= n_cells[1])
    c2[1] -= n_cells[1];
  else if (c2[1] < 0)
    c2[1] += n_cells[1];

  if (c2[2] >= n_cells[2])
    c2[2] -= n_cells[2];

  return c2[0] + c2[1] * n_cells[0] + c2[2] * n_cells[0] * n_cells[1];
} //}}}

// TODO: use Jacobi method
// TODO: use SYSTEM
// calculate gyration tensor and various shape descriptors //{{{
void Gyration(int n, int *list, COUNT Counts, BEADTYPE *BeadType,
              BEAD **Bead, double eigen[3]) {
  // gyration tensor (3x3 array)
  // use long double to ensure precision -- previous problem with truncation in
  // short chains
  // struct Tensor {
  //   LONGVECTOR x, y, z;
  // } GyrationTensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  long double GyrationTensor[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  double cog[3];
  GeomCentre(n, list, *Bead, cog);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    for (int dd = 0; dd < 3; dd++) {
      (*Bead)[list[i]].Position[dd] -= cog[dd];
    }
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    int id = list[i];
    GyrationTensor[0][0] += Bead[id]->Position[0] * Bead[id]->Position[0];
    GyrationTensor[0][1] += Bead[id]->Position[0] * Bead[id]->Position[1];
    GyrationTensor[0][2] += Bead[id]->Position[0] * Bead[id]->Position[2];
    GyrationTensor[1][1] += Bead[id]->Position[1] * Bead[id]->Position[1];
    GyrationTensor[1][2] += Bead[id]->Position[1] * Bead[id]->Position[2];
    GyrationTensor[2][2] += Bead[id]->Position[2] * Bead[id]->Position[2];
  }
  GyrationTensor[0][0] /= n;
  GyrationTensor[0][1] /= n;
  GyrationTensor[0][2] /= n;
  GyrationTensor[1][1] /= n;
  GyrationTensor[1][2] /= n;
  GyrationTensor[2][2] /= n; //}}}

  // characteristic polynomial:
  // a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0
  long double a_cube = -1;
  long double b_cube = GyrationTensor[0][0] +
                       GyrationTensor[1][1] +
                       GyrationTensor[2][2];
  long double c_cube = -GyrationTensor[0][0] * GyrationTensor[1][1] -
                       GyrationTensor[0][0] * GyrationTensor[2][2] -
                       GyrationTensor[1][1] * GyrationTensor[2][2] +
                       SQR(GyrationTensor[1][2]) +
                       SQR(GyrationTensor[0][1]) +
                       SQR(GyrationTensor[0][2]);
  long double d_cube =
      +GyrationTensor[0][0] * GyrationTensor[1][1] * GyrationTensor[2][2] +
      2 * GyrationTensor[0][1] * GyrationTensor[1][2] * GyrationTensor[0][2] -
      SQR(GyrationTensor[0][2]) * GyrationTensor[1][1] -
      SQR(GyrationTensor[0][1]) * GyrationTensor[2][2] -
      SQR(GyrationTensor[1][2]) * GyrationTensor[0][0];

  // first root: either 0 or Newton's iterative method to get it //{{{
  long double root0 = 0;
  if (fabsl(d_cube) > 0.0000000001L) {
    // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
    long double a_deriv = 3 * a_cube;
    long double b_deriv = 2 * b_cube;
    long double c_deriv = c_cube;

    long double root1 = 1;

    while (fabsl(root0 - root1) > 0.0000000001L) {
      long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) +
                             c_cube * root0 + d_cube);
      long double f_deriv_root0 =
          (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
      root1 = root0 - f_root0 / f_deriv_root0;

      // swap root0 and root1 for the next iteration
      long double tmp = root0;
      root0 = root1;
      root1 = tmp;
    }
  } //}}}

  // find parameters of a quadratic equation a_quad*x^2+b_quad*x+c_quad=0 //{{{
  /*
   * derived by division:
   * (x^3+(b_cube/a_cube)*x^2+(c_cube/a_cube)*x+(d_cube/a_cube)):(x-root0)
   */
  long double a_quad = 1;
  long double b_quad = b_cube / a_cube + root0;
  long double c_quad =
      SQR(root0) + b_cube / a_cube * root0 + c_cube / a_cube; //}}}
  // calculate & sort eigenvalues
  long double e[3];
  e[0] = root0; // found out by Newton's method
  // roots of the quadratic equation
  e[1] = (-b_quad + sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
  e[2] = (-b_quad - sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
  for (int dd = 0; dd < 3; dd++) {
    eigen[dd] = e[dd];
  }
  SortArrayDouble(eigen, 3, 0);
} //}}}

// evaluate contacts between molecules, creating aggregates //{{{
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      int contacts, int **contact) {
  COUNT *Count = &System->Count;
  // go over all pairs of molecules
  for (int i = 1; i < Count->Molecule; i++) {
    for (int j = 0; j < i; j++) {
      if (System->Molecule[i].InTimestep && System->Molecule[j].InTimestep) {
        int agg_i = System->Molecule[i].Aggregate,
            agg_j = System->Molecule[j].Aggregate;
        // if molecules 'i' and 'j' are in contact, put them into one aggregate
        if (contact[i][j] >= contacts) { //{{{
          // create new aggregate if molecule 'j' isn'it in any
          if (agg_j == -1) {
            agg_j = Count->Aggregate;
            System->Molecule[j].Aggregate = agg_j;
            Aggregate[agg_j].nMolecules = 1;
            Aggregate[agg_j].Molecule[0] = j;
            Count->Aggregate++;
          }
          /*
           * add molecule 'i' to aggregate 'j' (which contains molecule 'j')
           * if molecule 'i' isn't in any aggregate
           */
          if (agg_i == -1) {
            int mols = Aggregate[agg_j].nMolecules;
            Aggregate[agg_j].Molecule[mols] = i;
            Aggregate[agg_j].nMolecules++;
            System->Molecule[i].Aggregate = agg_j;
          }
          /*
           * if molecules 'i' and 'j' are in different aggregate,
           * unite those aggregates
           */
          if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {
            // add molecules from aggregate 'i' to aggregate 'j'
            int n_mol_old = Aggregate[agg_j].nMolecules;
            Aggregate[agg_j].nMolecules += Aggregate[agg_i].nMolecules;
            for (int k = n_mol_old; k < Aggregate[agg_j].nMolecules; k++) {
              int mol = Aggregate[agg_i].Molecule[k-n_mol_old];
              Aggregate[agg_j].Molecule[k] = mol;
              System->Molecule[mol].Aggregate = agg_j;
            }
            // move aggregates with id greater then agg_i to id-1
            for (int k = (agg_i + 1); k < Count->Aggregate; k++) {
              Aggregate[k-1].nMolecules = Aggregate[k].nMolecules;
              // move every molecule from aggregate 'k' to aggregate 'k-1'
              for (int l = 0; l < Aggregate[k].nMolecules; l++) {
                int mol = Aggregate[k].Molecule[l];
                Aggregate[k-1].Molecule[l] = mol;
                System->Molecule[mol].Aggregate = k - 1;
              }
            }
            // reduce number of aggregates since the two aggregates were merged
            Count->Aggregate--;
          } //}}}
        /*
         * or if molecules 'i' and 'j' aren't in contact, and molecule 'j' isn't
         * in any aggregate, create new aggregate for molecule 'j'
         */
        } else if (agg_j == -1) { //{{{
          agg_j = Count->Aggregate;
          System->Molecule[j].Aggregate = agg_j;
          Aggregate[agg_j].nMolecules = 1;
          Aggregate[agg_j].Molecule[0] = j;
          Count->Aggregate++;
        } //}}}
      }
    }
  }
  // check if the highest id residue is in an aggregate //{{{
  bool test = false;
  for (int i = 0; i < Count->Aggregate; i++) {
    for (int j = 1; j < Aggregate[i].nMolecules; j++) {
      if (Aggregate[i].Molecule[j] == (Count->Molecule - 1)) {
        test = 1;
      }
    }
  } //}}}
  // if the highest id residue isn't in an aggregate, create separate one //{{{
  if (!test) {
    int aggs = Count->Aggregate;
    Aggregate[aggs].nMolecules = 1;
    Aggregate[aggs].Molecule[0] = Count->Molecule - 1;

    Count->Aggregate++;
  } //}}}
} //}}}

// RemovePBCAggregates() //{{{
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate,
                         SYSTEM *System) {
  COUNT *Count = &System->Count;

  int **mol_eligible_beads = malloc(Count->MoleculeType * sizeof(int *));
  int *count_eligible_beads = malloc(Count->MoleculeType *
                                     sizeof *count_eligible_beads);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    mol_eligible_beads[i] = malloc(mt->nBeads * sizeof(int));
    count_eligible_beads[i] = 0;
    for (int j = 0; j < mt->nBeads; j++) {
      if (System->BeadType[mt->Bead[j]].Flag) {
        mol_eligible_beads[i][count_eligible_beads[i]] = j;
        count_eligible_beads[i]++;
      }
    }
  }

  double (*box)[3] = &System->Box.Length;
  // helper array indicating whether molecules already moved
  int *list_moved = calloc(Count->Molecule, sizeof *list_moved),
      *list_unmoved = calloc(Count->Molecule, sizeof *list_unmoved);
  // go through aggregates larger than As=1, knitting together //{{{
  for (int i = 0; i < Count->Aggregate; i++) {
    int count_moved = 0,
    count_unmoved = Aggregate[i].nMolecules - 1;
    // set first molecule as already moved
    list_moved[count_moved++] = 0;
    // set all other molecules as unmoved
    for (int j = 0; j < (Aggregate[i].nMolecules - 1); j++) {
      list_unmoved[j] = j + 1;
    }
    while (count_unmoved > 0) {
      // go through all molecule pairs
      for (int jj = 0; jj < count_moved; jj++) {
        int j = list_moved[jj];
        for (int kk = 0; kk < count_unmoved; kk++) {
          int k = list_unmoved[kk];
          bool moved = false;
          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          int mol1 = Aggregate[i].Molecule[j],
              mol2 = Aggregate[i].Molecule[k];
          int mtype1 = System->Molecule[mol1].Type,
              mtype2 = System->Molecule[mol2].Type;

          // go through all bead pairs in the two molecules
          for (int ll = 0; ll < count_eligible_beads[mtype1]; ll++) {
            int l = mol_eligible_beads[mtype1][ll];
            int bead1 = System->Molecule[mol1].Bead[l];
            BEAD *b1 = &System->Bead[bead1];
            for (int mm = 0; mm < count_eligible_beads[mtype2]; mm++) {
              int m = mol_eligible_beads[mtype2][mm];
              int bead2 = System->Molecule[mol2].Bead[m];
              BEAD *b2 = &System->Bead[bead2];
              // calculate distance between 'bead1' and 'bead2'
              double dist[3];
              Distance(b1->Position, b2->Position, *box, dist);
              dist[0] = VECTORLENGTH(dist);
              // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
              if (dist[0] <= distance) {
                // distance vector between 'bead1' and 'bead2'
                for (int dd = 0; dd < 3; dd++) {
                  dist[dd] = b1->Position[dd] - b2->Position[dd];
                }
                // if 'bead1' and 'bead2' are too far, move 'mol2' //{{{
                for (int dd = 0; dd < 3; dd++) {
                  while (dist[dd] > ((*box)[dd] / 2)) {
                    for (int n = 0; n < System->MoleculeType[mtype2].nBeads; n++) {
                      int id = System->Molecule[mol2].Bead[n];
                      System->Bead[id].Position[dd] += (*box)[dd];
                    }
                    dist[dd] = b1->Position[dd] - b2->Position[dd];
                  }
                  while (dist[dd] <= -((*box)[dd] / 2)) {
                    for (int n = 0; n < System->MoleculeType[mtype2].nBeads; n++) {
                      int id = System->Molecule[mol2].Bead[n];
                      System->Bead[id].Position[dd] -= (*box)[dd];
                    }
                    dist[dd] = b1->Position[dd] - b2->Position[dd];
                  }
                } //}}}
                moved = true;
                for (int x = kk; x < count_unmoved; x++) {
                  list_unmoved[x] = list_unmoved[x+1];
                }
                count_unmoved--;
                list_moved[count_moved++] = k;
                // skip remainder of 'mol2' (or 'k')
                break;
              }
            }
            if (moved) {
              break;
            }
          }
        }
      }
    }
  } //}}}
  // free(moved);
  free(list_moved);
  free(list_unmoved);
  free(count_eligible_beads);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(mol_eligible_beads[i]);
  }
  free(mol_eligible_beads);
  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Count->Aggregate; i++) {
    double com[3];
    CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, *System, com);
    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    int move[3];
    for (int dd = 0; dd < 3; dd++) {
      move[dd] = com[dd] / (*box)[dd];
      if (com[dd] < 0) {
        move[dd]--;
      }
    }
    // move all the beads
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[bead].Position[dd] -= move[dd] * (*box)[dd];
      }
    }
  } //}}}
} //}}}

// should given step be used for calculations? //{{{
bool UseStep(COMMON_OPT opt, int step) {
  if (step >= opt.start &&
      (step <= opt.end || opt.end == -1) &&
      ((step - opt.start) % opt.skip) == 0) {
    return true;
  } else {
    return false;
  }
} //}}}
