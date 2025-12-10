#include "AnalysisTools.h"
#include "Pairs.h"

// create a cell linked list
static vec3i LinkedList(const SYSTEM System, int **Head, int **Link,
                        const double cell_size);
// cell selectors for bead pair linked list traversal
static inline int SelectCell1(const vec3i c1, const vec3i n_cells);
static inline int SelectCell2(const vec3i c1, const vec3i n_cells,
                              const vec3i neighbour[13], int n);

// linked list traversal //{{{
void TraverseLinkedListPairs(const SYSTEM System, const double cell_size,
                             pair_cb_t pair_callback, void *pair_ud,
                             check_cb_t check_callback, void *check_ud) {
  int *Head, *Link;
  vec3i n_cells = LinkedList(System, &Head, &Link, cell_size);

  // neighbour offsets: home cell + 13 half-shell
  const vec3i neighbour[14] = {
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
          // loop over all 14 neighbour offsets (home + 13 neighbours)
          for (int k = 0; k < 14; k++) {
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
// brute traversal for high cell_size but linked list traversal otherwise //{{{
void TraversePairs(const SYSTEM System, const double cell_size,
                   pair_cb_t pair_callback, void *pair_ud,
                   check_cb_t check_callback, void *check_ud) {
  bool linked = true;
  if ((Min3(System.Box.Length[0],
            System.Box.Length[1],
            System.Box.Length[2]) / 3) < cell_size) {
    linked = false;
  }
  if (linked) {
    TraverseLinkedListPairs(System, cell_size, pair_callback, pair_ud,
                            check_callback, check_ud);
  } else {
    TraverseBrutePairs(System, pair_callback, pair_ud,
                       check_callback, check_ud);
  }
} //}}}

// create a cell linked list //{{{
static vec3i LinkedList(const SYSTEM System, int **Head, int **Link,
                        const double cell_size) {
  const double (*box)[3] = &System.Box.Length;
  const COUNT *Count = &System.Count;
  vec3d rl;
  vec3i n_cells;
  // compute number of cells along each axis
  for (int dd = 0; dd < 3; dd++) {
    rl.v[dd] = (*box)[dd] / cell_size;
    n_cells.v[dd] = (int)(rl.v[dd]);
    if (n_cells.v[dd] < 3) {
      err_msg("cell size too small for cut-off in linked list");
      PrintError();
      exit(1);
    }
    rl.v[dd] = (double)n_cells.v[dd] / (*box)[dd]; // inverse length
  }
  // allocate lists
  int cells = n_cells.x * n_cells.y * n_cells.z;
  *Head = malloc(sizeof **Head * cells);
  *Link = malloc(sizeof **Link * Count->BeadCoor);
  for (int i = 0; i < cells; i++) {
    (*Head)[i] = -1;
  }
  // insert beads
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    vec3i c;
    for (int dd = 0; dd < 3; dd++) {
      c.v[dd] = (int)(bead->Position.v[dd] * rl.v[dd]);
      if (c.v[dd] == n_cells.v[dd]) { // guard FP boundary
        c.v[dd] = n_cells.v[dd] - 1;
      }
    }
    int cell = c.x + c.y * n_cells.x + c.z * n_cells.x * n_cells.y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  }
  return n_cells;
}
static inline int SelectCell1(const vec3i c1, const vec3i n_cells) {
  return c1.x + c1.y * n_cells.x + c1.z * n_cells.x * n_cells.y;
}
static inline int SelectCell2(const vec3i c1, const vec3i n_cells,
                              const vec3i neighbour[13], int n) {
  vec3i c2;
  for (int dd = 0; dd < 3; dd++) {
    c2.v[dd] = (c1.v[dd] + neighbour[n].v[dd] + n_cells.v[dd]) % n_cells.v[dd];
  }
  return SelectCell1(c2, n_cells);
}
//}}}
