#include "Pairs.h"
#include "Errors.h"
#include "General.h"
#include "MathUtils.h"

/*
 * norm_axis convention (shared by all functions below): -1 means full 3D
 * binning; 0/1/2 means the given axis is not binned (a single cell spans it),
 * so pairs are found by their in-plane separation only. Useful for 2D
 * (slit-like) systems where the distance along norm_axis is irrelevant.
 */
// cross product of two 3D vectors
static inline vec3d Cross(const vec3d a, const vec3d b) {
  return (vec3d){ .v = { a.y * b.z - a.z * b.y,
                         a.z * b.x - a.x * b.z,
                         a.x * b.y - a.y * b.x } };
}
// is the simulation box orthogonal? (same epsilon as DistancePBC())
static inline bool BoxIsOrthogonal(const BOX *box) {
  return fabs(box->alpha - 90) < 1e-5 &&
         fabs(box->beta  - 90) < 1e-5 &&
         fabs(box->gamma - 90) < 1e-5;
}
/*
 * Perpendicular widths of the (possibly triclinic) cell, i.e. the spacings
 * between opposite periodic faces. This is the largest cell extent along each
 * lattice direction that still lets the half-shell +/-1 neighbour search see
 * every pair within the cut-off. For orthogonal boxes these are the box
 * lengths; for triclinic boxes they are Volume / |edge_j x edge_k|.
 */
static vec3d CellWidths(const BOX *box);
// create a cell linked list
static vec3i LinkedList(const SYSTEM System, int **Head, int **Link,
                        const double cell_size, const int norm_axis);
// fill neighbour cell offsets (home + half-shell); returns their count
static int FillNeighbours(vec3i neighbour[14], const int norm_axis);
// cell selectors for bead pair linked list traversal
static inline int SelectCell1(const vec3i c1, const vec3i n_cells);
static inline int SelectCell2(const vec3i c1, const vec3i n_cells,
                              const vec3i neighbour[14], int n);
// traversal core shared by the 3D and 2D variants
static void TraverseLL(const SYSTEM System, const double cell_size,
                       const int norm_axis,
                       pair_cb_t pair_callback, void *pair_ud,
                       check_cb_t check_callback, void *check_ud);

// linked list traversal //{{{
static void TraverseLL(const SYSTEM System, const double cell_size,
                       const int norm_axis,
                       pair_cb_t pair_callback, void *pair_ud,
                       check_cb_t check_callback, void *check_ud) {
  int *Head, *Link;
  vec3i n_cells = LinkedList(System, &Head, &Link, cell_size, norm_axis);
  vec3i neighbour[14];
  int n_neigh = FillNeighbours(neighbour, norm_axis);
  vec3i c1;
  for (c1.z = 0; c1.z < n_cells.z; c1.z++) {
    for (c1.y = 0; c1.y < n_cells.y; c1.y++) {
      for (c1.x = 0; c1.x < n_cells.x; c1.x++) {
        int cell1 = SelectCell1(c1, n_cells);
        int i = Head[cell1];
        while (i != -1) {
          if (!check_callback(i, System, check_ud)) {
            i = Link[i];
            continue;
          }
          // loop over all neighbour offsets (home + half-shell)
          for (int k = 0; k < n_neigh; k++) {
            int cell2 = SelectCell2(c1, n_cells, neighbour, k);

            int j = Head[cell2];
            while (j != -1) {
              if (!check_callback(j, System, check_ud)) {
                j = Link[j];
                continue;
              }
              // avoid double-counting in home cell
              if (cell1 != cell2 || i < j) {
                pair_callback(i, j, System, pair_ud);
              }
              j = Link[j];
            }
          }
          i = Link[i];
        }
      }
    }
  }
  free(Head);
  free(Link);
}
void TraverseLinkedListPairs(const SYSTEM System, const double cell_size,
                             pair_cb_t pair_callback, void *pair_ud,
                             check_cb_t check_callback, void *check_ud) {
  TraverseLL(System, cell_size, -1, pair_callback, pair_ud,
             check_callback, check_ud);
}
void TraverseLinkedListPairs2D(const SYSTEM System, const double cell_size,
                               const int norm_axis,
                               pair_cb_t pair_callback, void *pair_ud,
                               check_cb_t check_callback, void *check_ud) {
  TraverseLL(System, cell_size, norm_axis, pair_callback, pair_ud,
             check_callback, check_ud);
} //}}}
// brute O(N^2) traversal //{{{
void TraverseBrutePairs(const SYSTEM System,
                        pair_cb_t pair_callback, void *pair_ud,
                        check_cb_t check_callback, void *check_ud) {
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    if (!check_callback(i, System, check_ud)) {
      continue;
    }
    for (int j = (i + 1); j < System.Count.BeadCoor; j++) {
      if (!check_callback(j, System, check_ud)) {
        continue;
      }
      pair_callback(i, j, System, pair_ud);
    }
  }
} //}}}
// brute traversal for high cell_size, linked list traversal otherwise //{{{
void TraversePairs(const SYSTEM System, const double cell_size,
                   pair_cb_t pair_callback, void *pair_ud,
                   check_cb_t check_callback, void *check_ud) {
  bool linked = true;
  // decide on the perpendicular cell widths (equal box lengths when orthogonal)
  vec3d width = CellWidths(&System.Box);
  if ((Min3(width.x, width.y, width.z) / 3) < cell_size) {
    linked = false;
  }
  if (linked) {
    TraverseLinkedListPairs(System, cell_size, pair_callback, pair_ud,
                            check_callback, check_ud);
  } else {
    TraverseBrutePairs(System, pair_callback, pair_ud,
                       check_callback, check_ud);
  }
}
void TraversePairs2D(const SYSTEM System, const double cell_size,
                     const int norm_axis,
                     pair_cb_t pair_callback, void *pair_ud,
                     check_cb_t check_callback, void *check_ud) {
  if (norm_axis < 0 || norm_axis > 2) {
    err_msg("TraversePairs2D(): norm_axis must be 0, 1, or 2");
    PrintError();
    exit(1);
  }
  // only the two in-plane axes need to fit at least 3 cells
  vec3d width = CellWidths(&System.Box);
  double min_plane = -1;
  for (int dd = 0; dd < 3; dd++) {
    if (dd != norm_axis &&
        (min_plane == -1 || width.v[dd] < min_plane)) {
      min_plane = width.v[dd];
    }
  }
  if ((min_plane / 3) < cell_size) {
    TraverseBrutePairs(System, pair_callback, pair_ud,
                       check_callback, check_ud);
  } else {
    TraverseLinkedListPairs2D(System, cell_size, norm_axis,
                              pair_callback, pair_ud,
                              check_callback, check_ud);
  }
} //}}}

static vec3d CellWidths(const BOX *box) {
  if (BoxIsOrthogonal(box)) {
    return box->Length;
  }
  // lattice edge vectors are the columns of the transform matrix
  vec3d a = { .v = {box->transform[0][0], box->transform[1][0],
                    box->transform[2][0]} };
  vec3d b = { .v = {box->transform[0][1], box->transform[1][1],
                    box->transform[2][1]} };
  vec3d c = { .v = {box->transform[0][2], box->transform[1][2],
                    box->transform[2][2]} };
  vec3d w;
  w.v[0] = box->Volume / VectLength(Cross(b, c));
  w.v[1] = box->Volume / VectLength(Cross(c, a));
  w.v[2] = box->Volume / VectLength(Cross(a, b));
  return w;
}
// create a cell linked list //{{{
/*
 * Cells are laid out in fractional (lattice) space so the same code handles
 * orthogonal and triclinic boxes: a bead's fractional coordinate is r/Length
 * per axis when orthogonal, or inverse*r (as in DistancePBC) when triclinic.
 * Fractional coordinates are wrapped into [0,1) before binning, so a bead that
 * sits outside [0,Length) - which is normal for the real Cartesian coordinates
 * of a tilted cell - can never index outside the cell array.
 */
static vec3i LinkedList(const SYSTEM System, int **Head, int **Link,
                        const double cell_size, const int norm_axis) {
  const BOX *box = &System.Box;
  const COUNT *Count = &System.Count;
  const bool ortho = BoxIsOrthogonal(box);
  const vec3d width = CellWidths(box);
  vec3i n_cells;
  // compute number of cells along each axis
  for (int dd = 0; dd < 3; dd++) {
    if (dd == norm_axis) { // single cell spans the non-binned axis
      n_cells.v[dd] = 1;
      continue;
    }
    n_cells.v[dd] = (int)(width.v[dd] / cell_size);
    if (n_cells.v[dd] < 3) {
      err_msg("cell size too small for cut-off in linked list");
      PrintError();
      exit(1);
    }
  }
  // allocate lists
  int cells = n_cells.x * n_cells.y * n_cells.z;
  *Head = malloc(sizeof **Head * cells);
  *Link = malloc(sizeof **Link * Count->BeadCoor);
  if (!*Head || !*Link) {
    ErrorAlloc("linked list Head/Link");
  }
  for (int i = 0; i < cells; i++) {
    (*Head)[i] = -1;
  }
  // insert beads
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    const vec3d pos = System.Bead[id].Position;
    // fractional coordinate: orthogonal is diagonal (r/Length), triclinic uses
    // the Cartesian->fractional matrix (inverse might be unset for orthogonal
    // boxes read without CalculateBoxData(), hence the split)
    vec3d s;
    if (ortho) {
      for (int dd = 0; dd < 3; dd++) {
        s.v[dd] = pos.v[dd] / box->Length.v[dd];
      }
    } else {
      for (int dd = 0; dd < 3; dd++) {
        s.v[dd] = box->inverse[dd][0] * pos.v[0] +
                  box->inverse[dd][1] * pos.v[1] +
                  box->inverse[dd][2] * pos.v[2];
      }
    }
    vec3i c;
    for (int dd = 0; dd < 3; dd++) {
      if (dd == norm_axis) { // non-binned axis maps to the single cell 0
        c.v[dd] = 0;
        continue;
      }
      double f = s.v[dd] - floor(s.v[dd]); // wrap fractional coord into [0,1)
      c.v[dd] = (int)(f * n_cells.v[dd]);
      if (c.v[dd] >= n_cells.v[dd]) { // guard FP boundary (f -> 1.0)
        c.v[dd] = n_cells.v[dd] - 1;
      } else if (c.v[dd] < 0) { // guard FP boundary (f -> tiny negative)
        c.v[dd] = 0;
      }
    }
    int cell = c.x + c.y * n_cells.x + c.z * n_cells.x * n_cells.y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  }
  return n_cells;
}
static int FillNeighbours(vec3i neighbour[14], const int norm_axis) {
  if (norm_axis == -1) { // 3D: home cell + 13 half-shell
    const vec3i neigh3d[14] = {
      { .v = { 0, 0, 0} },
      { .v = { 1, 0, 0} },
      { .v = { 1, 1, 0} },
      { .v = {-1, 1, 0} },
      { .v = { 0, 1, 0} },
      { .v = { 0, 0, 1} },
      { .v = {-1, 0, 1} },
      { .v = { 1, 0, 1} },
      { .v = {-1,-1, 1} },
      { .v = { 0,-1, 1} },
      { .v = { 1,-1, 1} },
      { .v = {-1, 1, 1} },
      { .v = { 0, 1, 1} },
      { .v = { 1, 1, 1} },
    };
    for (int i = 0; i < 14; i++) {
      neighbour[i] = neigh3d[i];
    }
    return 14;
  }
  // 2D: home cell + 4 half-shell offsets placed in the two in-plane axes
  const int neigh2d[5][2] = { {0, 0}, {1, 0}, {1, 1}, {-1, 1}, {0, 1} };
  int plane[2], n = 0;
  for (int dd = 0; dd < 3; dd++) {
    if (dd != norm_axis) {
      plane[n++] = dd;
    }
  }
  for (int i = 0; i < 5; i++) {
    neighbour[i] = (vec3i){ .v = {0, 0, 0} };
    neighbour[i].v[plane[0]] = neigh2d[i][0];
    neighbour[i].v[plane[1]] = neigh2d[i][1];
  }
  return 5;
}
static inline int SelectCell1(const vec3i c1, const vec3i n_cells) {
  return c1.x + c1.y * n_cells.x + c1.z * n_cells.x * n_cells.y;
}
static inline int SelectCell2(const vec3i c1, const vec3i n_cells,
                              const vec3i neighbour[14], int n) {
  vec3i c2;
  for (int dd = 0; dd < 3; dd++) {
    c2.v[dd] = (c1.v[dd] + neighbour[n].v[dd] + n_cells.v[dd]) % n_cells.v[dd];
  }
  return SelectCell1(c2, n_cells);
}
//}}}
