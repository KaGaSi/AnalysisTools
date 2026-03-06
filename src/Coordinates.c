#include "Coordinates.h"
#include "Errors.h"
#include "General.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>

// STATIC DECLARATIONS
// TODO: the whole fractional stuff - not sure if it works
static void FractionalCoor(SYSTEM *System, const int mode);
// comparison function for ascending order (used by qsort in Gyration)
static int Compare(const void *a, const void *b);

static int Compare(const void *a, const void *b) {
  return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
}

// STATIC IMPLEMENTATIONS
// transform to/from fractional coordinates //{{{
/* HOW TO CALCULATE DISTANCE IN TRICLINIC SYSTEM //{{{
//double dist[3];
//dist[0] = (*Bead)[0].Position[0] - (*Bead)[10].Position[0];
//dist[1] = (*Bead)[0].Position[1] - (*Bead)[10].Position[1];
//dist[2] = (*Bead)[0].Position[2] - (*Bead)[10].Position[2];
//printf("dist1 = (%lf, %lf, %lf) = %lf\n", dist[0], dist[1], dist[2],
sqrt(Square(dist[0])+Square(dist[1])+Square(dist[2])));

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
sqrt(Square(dist[0])+Square(dist[1])+Square(dist[2])));
*/ //}}}
// mode=0 ... to fractional; mode=1 ... from fractional
static void FractionalCoor(SYSTEM *System, const int mode) {
  if (mode != 0 && mode != 1) {
    err_msg("mode in FractionalCoor() must be 0 or 1");
    PrintError();
    exit(1);
  }
  BOX *box = &System->Box;
  if (fabs(box->alpha - 90) > 0.00001 ||
      fabs(box->beta - 90) > 0.00001 ||
      fabs(box->gamma - 90) > 0.00001) {
    for (int i = 0; i < System->Count.BeadCoor; i++) {
      int id = System->BeadCoor[i];
      BEAD *b = &System->Bead[id];
      double new[3] = {0, 0, 0};
      for (int dd = 0; dd < 3; dd++) {
        if (mode == 0) {
          new[dd] = box->inverse[dd][0] * b->Position.v[0] +
                    box->inverse[dd][1] * b->Position.v[1] +
                    box->inverse[dd][2] * b->Position.v[2];
        } else {
          vec3d fraction = { .v = { b->Position.v[0] / box->Length.v[dd],
                                    b->Position.v[1] / box->Length.v[dd],
                                    b->Position.v[2] / box->Length.v[dd]}};
          new[dd] = box->transform[dd][0] * fraction.v[0] +
                    box->transform[dd][1] * fraction.v[1] +
                    box->transform[dd][2] * fraction.v[2];
        }
      }
      for (int dd = 0; dd < 3; dd++) {
        b->Position.v[dd] = new[dd] * System->Box.Length.v[dd];
      }
    }
  }
} //}}}

// put given vector into range <0,BoxLength) //{{{
vec3d RestorePBC(const vec3d coor, const vec3d BoxLength) {
  vec3d out;
  for (int dd = 0; dd < 3; dd++) {
    // by how many boxlength should the distance be changed?
    int move = floor(coor.v[dd] / BoxLength.v[dd]);
    // transform it into <0,BoxLength) range
    out.v[dd] = coor.v[dd] - move * BoxLength.v[dd];
  }
  return out;
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
void RemovePBCMolecule(int mol_id, SYSTEM *System) {
  MOLECULE *mol = &System->Molecule[mol_id];
  MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
  // do nothing for molecule that isn't in the timestep or has no bonds
  if (!mol->InTimestep || mt->nBonds == 0) {
    return;
  }
  BOX *box = &System->Box;
  // arrays holding bonds already connected and yet unconnected
  int *connected = calloc(mt->nBonds, sizeof *connected);
  int *unconnected = calloc(mt->nBonds, sizeof *unconnected);
  int count_unconnected = 0;
  // 1)
  for (int i = 0; i < mt->nBonds; i++) {
    int id[2] = {mol->Bead[mt->Bond[i][0]], mol->Bead[mt->Bond[i][1]]};
    BEAD *b_1 = &System->Bead[id[0]];
    BEAD *b_2 = &System->Bead[id[1]];
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
    return;
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
        int bond[2]; // the connected and unconnected bonds
        int con[2]; // beads in the already connected bond (bond[0])
        int uncon[2]; // beads in the yet unconneced bond (bond[1])
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
      vec3d dist;
      if (!moved[id[0]] && moved[id[1]]) {
        dist = Distance(b_2->Position.v, b_1->Position.v, box->OrthoLength);
        for (int dd = 0; dd < 3; dd++) {
          b_1->Position.v[dd] = b_2->Position.v[dd] - dist.v[dd];
        }
        moved[id[0]] = true;
      } else if (moved[id[0]] && !moved[id[1]]) {
        dist = Distance(b_1->Position.v, b_2->Position.v, box->OrthoLength);
        for (int dd = 0; dd < 3; dd++) {
          b_2->Position.v[dd] = b_1->Position.v[dd] - dist.v[dd];
        }
        moved[id[1]] = true;
      }
    }
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
    move[dd] = cog[dd] / box->OrthoLength.v[dd];
    if (cog[dd] < 0) {
      move[dd]--;
    }
  }
  for (int j = 0; j < mt->nBeads; j++) {
    int bead = mol->Bead[j];
    for (int dd = 0; dd < 3; dd++) {
      System->Bead[bead].Position.v[dd] -= move[dd] * box->OrthoLength.v[dd];
    }
  } //}}}
} //}}}
// wrap coordinates into simulation box and/or join molecules //{{{
void WrapJoinCoordinates(SYSTEM *System, const bool wrap, const bool join) {
  if (System->Box.Volume != -1 && (wrap || join)) {
    // transform coordinates into fractional ones for non-orthogonal box
    FractionalCoor(System, 0); // TODO: does fracional work?
    if (wrap) { // wrap coordinates into the simulation box
      for (int i = 0; i < System->Count.BeadCoor; i++) {
        int id = System->BeadCoor[i];
        BEAD *bead = &System->Bead[id];
        bead->Position = RestorePBC(bead->Position, System->Box.OrthoLength);
      }
    }
    if (join) { // join molecules by removing periodic boundary conditions
      for (int i = 0; i < System->Count.Molecule; i++) {
        RemovePBCMolecule(i, System);
      }
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FractionalCoor(System, 1); // TODO: does fracional work?
  }
} //}}}
// distance between two beads; in the range <-BoxLength/2,BoxLength/2) //{{{
vec3d Distance(const double id1[3], const double id2[3],
               const vec3d BoxLength) {
  vec3d out;
  // calculate distance, transforming it into <0,BoxLength) range
  for (int dd = 0; dd < 3; dd++) {
    out.v[dd] = id1[dd] - id2[dd] + BoxLength.v[dd] / 2;
  }
  out = RestorePBC(out, BoxLength);
  // transform the distance back to <-BoxLength/2,BoxLength/2) range
  for (int dd = 0; dd < 3; dd++) {
    out.v[dd] -= BoxLength.v[dd] / 2;
  }
  return out;
} //}}}
// calculate centre of mass for a list of beads //{{{
void CentreOfMass(const int n, const int *list,
                  const SYSTEM System, double com[3]) {
  for (int dd = 0; dd < 3; dd++) {
    com[dd] = 0;
  }
  double mass = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    BEAD *b = &System.Bead[id];
    BEADTYPE *bt = &System.BeadType[b->Type];
    if (bt->Mass == MASS) {
      for (int dd = 0; dd < 3; dd++) {
        com[dd] = 0;
        if (snprintf(ERROR_MSG, LINE, "unspecified mass: bead %s%d%s (%s%s%s)",
                     ErrYellow(), id, ErrCyan(),
                     ErrYellow(), bt->Name, ErrCyan()) < 0) {
          ErrorSnprintf();
        }
        PrintWarning();
        return;
      }
    }
    for (int dd = 0; dd < 3; dd++) {
      com[dd] += b->Position.v[dd] * bt->Mass;
    }
    mass += bt->Mass;
  }
  for (int dd = 0; dd < 3; dd++) {
    com[dd] /= mass;
  }
} //}}}
// calculate geometric centre for a list of beads //{{{
void GeomCentre(const int n, const int *list, const BEAD *Bead, double gc[3]) {
  InitDoubleArray(gc, 3, 0);
  int count = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    if (Bead[id].InTimestep) {
      for (int dd = 0; dd < 3; dd++) {
        gc[dd] += Bead[id].Position.v[dd];
      }
      count++;
    }
  }
  for (int dd = 0; dd < 3; dd++) {
    gc[dd] /= count;
  }
} //}}}
// calculate gyration tensor and various shape descriptors //{{{
vec3d Gyration(const int n, const int *list, SYSTEM *System) {
  // gyration tensor (3x3 array)
  long double GyrationTensor[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  double cog[3];
  GeomCentre(n, list, System->Bead, cog);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    for (int dd = 0; dd < 3; dd++) {
      System->Bead[list[i]].Position.v[dd] -= cog[dd];
    }
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    int id = list[i];
    vec3d *pos = &System->Bead[id].Position;
    GyrationTensor[0][0] += pos->v[0] * pos->v[0];
    GyrationTensor[0][1] += pos->v[0] * pos->v[1];
    GyrationTensor[0][2] += pos->v[0] * pos->v[2];
    GyrationTensor[1][1] += pos->v[1] * pos->v[1];
    GyrationTensor[1][2] += pos->v[1] * pos->v[2];
    GyrationTensor[2][2] += pos->v[2] * pos->v[2];
  }
  GyrationTensor[0][0] /= n;
  GyrationTensor[0][1] /= n;
  GyrationTensor[0][2] /= n;
  GyrationTensor[1][1] /= n;
  GyrationTensor[1][2] /= n;
  GyrationTensor[2][2] /= n; //}}}
  // Define the symmetric matrix (example 3x3 matrix)
  gsl_matrix *A = gsl_matrix_alloc(3, 3);
  gsl_vector *eigenvalues = gsl_vector_alloc(3);
  gsl_matrix *eigenvectors = gsl_matrix_alloc(3, 3);
  for (int dd1 = 0; dd1 < 3; dd1++) {
    for (int dd2 = 0; dd2 < 3; dd2++) {
      gsl_matrix_set(A, dd1, dd2, GyrationTensor[dd1][dd2]);
    }
  }
  // Perform the eigendecomposition using the Jacobi method
  gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(3);
  gsl_eigen_symmv(A, eigenvalues, eigenvectors, workspace);

  vec3d eigen;
  for (int dd = 0; dd < 3; dd++) {
    eigen.v[dd] = gsl_vector_get(eigenvalues, dd);
    if (fabs(eigen.v[dd]) < 1e-5) {
      eigen.v[dd] = 0;
    }
  }
  qsort(eigen.v, 3, sizeof(eigen.v[0]), Compare);

  // Free allocated memory
  gsl_matrix_free(A);
  gsl_vector_free(eigenvalues);
  gsl_matrix_free(eigenvectors);
  gsl_eigen_symmv_free(workspace);

  // error for negative eigenvalues - shouldn't happen
  if (eigen.x < 0 || eigen.y < 0 || eigen.z < 0) {
    snprintf(ERROR_MSG, LINE, "negative eigenvalues (%s%lf%s, %s%lf%s, "
             "%s%lf%s)", ErrYellow(), eigen.x, ErrCyan(), ErrYellow(),
             eigen.y, ErrCyan(), ErrYellow(), eigen.z, ErrCyan());
    PrintWarning();
  }
  return eigen;
} //}}}
