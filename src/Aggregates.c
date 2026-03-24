#include "Aggregates.h"
#include "Coordinates.h"
#include "Errors.h"
#include "General.h"
#include "MathUtils.h"
#include "Options.h"
#include "PairHash.h"

// TODO: AggPickerOptions should be in Options.c, no?

// evaluate bead contacts to assign molecules to aggregates using DBSCAN //{{{
/*
 * DBSCAN at the molecule level:
 *   - epsilon-neighbourhood of molecule m = all molecules with bead-bead
 *     contact count >= 'contacts'
 *   - core point: degree(m) >= min_pts
 *   - border point: degree(m) < neighbours but reachable from a core point
 *   - noise / isolated: not reachable from any core point -> singleton
 *
 * With neighbours = 1 the algorithm reduces to simple
 */
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      const int contacts, const int neighbours,
                      PairHash *contact) {
  COUNT *Count = &System->Count;

  // Step 1: build a CSR adjacency list from the contact hash.
  // Only pairs where both molecules are InTimestep and whose bead-contact
  // count meets the threshold are treated as edges.
  int *degree = calloc(Count->Molecule, sizeof *degree);
  for (khiter_t it = kh_begin(contact); it != kh_end(contact); ++it) {
    if (!kh_exist(contact, it)) {
      continue;
    }
    if (kh_val(contact, it) < (uint8_t)contacts) {
      continue;
    }
    uint64_t key = kh_key(contact, it);
    int i = PairHash_mol_i(key);
    int j = PairHash_mol_j(key);
    if (!System->Molecule[i].InTimestep ||
        !System->Molecule[j].InTimestep) {
      continue;
    }
    degree[i]++;
    degree[j]++;
  }
  int *offset = malloc((Count->Molecule + 1) * sizeof *offset);
  offset[0] = 0;
  for (int i = 0; i < Count->Molecule; i++) {
    offset[i + 1] = offset[i] + degree[i];
  }
  int total = offset[Count->Molecule];
  int *nbrs = malloc((total > 0 ? total : 1) * sizeof *nbrs);
  int *fill = calloc(Count->Molecule, sizeof *fill);
  for (khiter_t it = kh_begin(contact); it != kh_end(contact); ++it) {
    if (!kh_exist(contact, it)) {
      continue;
    }
    if (kh_val(contact, it) < (uint8_t)contacts) {
      continue;
    }
    uint64_t key = kh_key(contact, it);
    int i = PairHash_mol_i(key);
    int j = PairHash_mol_j(key);
    if (!System->Molecule[i].InTimestep ||
        !System->Molecule[j].InTimestep) continue;
    nbrs[offset[i] + fill[i]++] = j;
    nbrs[offset[j] + fill[j]++] = i;
  }
  free(fill);

  // Step 2: DBSCAN over the molecule graph.
  // label[i]: -1 = unvisited, -2 = noise/border candidate, >= 0 = cluster id
  int *label = malloc(Count->Molecule * sizeof *label);
  for (int i = 0; i < Count->Molecule; i++) {
    label[i] = -1;
  }
  int *queue = malloc(Count->Molecule * sizeof *queue);
  int cluster_id = 0;
  for (int i = 0; i < Count->Molecule; i++) {
    if (!System->Molecule[i].InTimestep || label[i] != -1) continue;
    if (degree[i] < neighbours) {
      label[i] = -2; // noise or potential border point
      continue;
    }
    // core point: seed a new cluster via BFS
    label[i] = cluster_id;
    int head = 0, tail = 0;
    queue[tail++] = i;
    while (head < tail) {
      int m = queue[head++];
      for (int k = offset[m]; k < offset[m + 1]; k++) {
        int n = nbrs[k];
        if (label[n] == -1) {
          // unvisited neighbour: assign to cluster
          label[n] = cluster_id;
          if (degree[n] >= neighbours) {
            queue[tail++] = n; // also a core point, expand further
          }
        } else if (label[n] == -2) {
          label[n] = cluster_id; // border point absorbed into cluster
        }
      }
    }
    cluster_id++;
  }
  free(queue);
  free(degree);
  free(offset);
  free(nbrs);

  // Step 3: assign cluster members to Aggregate structs.
  Count->Aggregate = 0;
  int *cluster_size = calloc(cluster_id, sizeof *cluster_size);
  int *cluster_fill = calloc(cluster_id, sizeof *cluster_fill);
  for (int i = 0; i < Count->Molecule; i++) {
    if (System->Molecule[i].InTimestep && label[i] >= 0) {
      cluster_size[label[i]]++;
    }
  }
  for (int c = 0; c < cluster_id; c++) {
    int agg = Count->Aggregate++;
    Aggregate[agg].nMolecules = cluster_size[c];
    Aggregate[agg].Molecule = s_realloc(Aggregate[agg].Molecule,
                                        cluster_size[c] *
                                        sizeof *Aggregate[agg].Molecule);
  }
  for (int i = 0; i < Count->Molecule; i++) {
    if (!System->Molecule[i].InTimestep || label[i] < 0) continue;
    int agg = label[i];
    Aggregate[agg].Molecule[cluster_fill[agg]++] = i;
    System->Molecule[i].Aggregate = agg;
  }
  free(cluster_size);
  free(cluster_fill);

  // singleton aggregates: noise molecules (label == -2) and any remaining
  // uncontacted InTimestep molecules (label == -1)
  for (int i = 0; i < Count->Molecule; i++) {
    if (System->Molecule[i].InTimestep && label[i] < 0) {
      int agg = Count->Aggregate++;
      System->Molecule[i].Aggregate = agg;
      Aggregate[agg].nMolecules = 1;
      Aggregate[agg].Molecule[0] = i;
    }
  }

  free(label);
} //}}}
// RemovePBCAggregates() //{{{
void RemovePBCAggregates(const double distance, const AGGREGATE *Aggregate,
                         SYSTEM *System, const bool *use_bt) {
  COUNT *Count = &System->Count;
  int **mol_eligible_beads = malloc(Count->MoleculeType * sizeof(int *));
  int *count_eligible_beads = malloc(Count->MoleculeType *
                                     sizeof *count_eligible_beads);
  WrapJoinCoordinates(System, false, true);
  bool eligible = false;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    mol_eligible_beads[i] = malloc(mt->nBeads * sizeof(int));
    count_eligible_beads[i] = 0;
    for (int j = 0; j < mt->nBeads; j++) {
      if (use_bt[mt->Bead[j]]) {
        mol_eligible_beads[i][count_eligible_beads[i]] = j;
        count_eligible_beads[i]++;
        eligible = true;
      }
    }
  }
  if (!eligible) {
    err_msg("RemovePBCAggregates(): no bead types for joining");
    PrintError();
    exit(1);
  }

  vec3d *box = &System->Box.Length;
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
      bool any_moved = false;
      // go through all molecule pairs
      for (int jj = 0; jj < count_moved; jj++) {
        int j = list_moved[jj];
        for (int kk = 0; kk < count_unmoved; kk++) {
          int k = list_unmoved[kk];
          bool moved = false;
          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          int mol1 = Aggregate[i].Molecule[j];
          int mol2 = Aggregate[i].Molecule[k];
          int mtype1 = System->Molecule[mol1].Type;
          int mtype2 = System->Molecule[mol2].Type;

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
              vec3d dist = Distance(b1->Position, b2->Position, *box);
              // move 'mol2' if 'bead1' and 'bead2' are in contact
              if (VectLength(dist) <= distance) {
                // multiples of box lengths that place b2 at b1 - dist
                vec3d shift;
                for (int dd = 0; dd < 3; dd++) {
                  shift.v[dd] = b1->Position.v[dd] - b2->Position.v[dd]
                                - dist.v[dd];
                }
                for (int n = 0; n < System->MoleculeType[mtype2].nBeads; n++) {
                  int id = System->Molecule[mol2].Bead[n];
                  for (int dd = 0; dd < 3; dd++) {
                    System->Bead[id].Position.v[dd] += shift.v[dd];
                  }
                }
                moved = true;
                for (int x = kk; x < count_unmoved; x++) {
                  list_unmoved[x] = list_unmoved[x+1];
                }
                count_unmoved--;
                list_moved[count_moved++] = k;
                any_moved = true;
                kk--;
                break;
              }
            }
            if (moved) {
              break;
            }
          }
        }
      }
      // guard against infinite loop if no bead pair is within range
      if (!any_moved) {
        snprintf(ERROR_MSG, LINE, "cannot join all molecules in aggregate "
                 "%s%d%s - no eligible bead pair within distance %s%g%s",
                 ErrYellow(), i, ErrCyan(), ErrYellow(), distance, ErrCyan());
        PrintWarning();
        break;
      }
    }
  } //}}}
  free(list_moved);
  free(list_unmoved);
  free(count_eligible_beads);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(mol_eligible_beads[i]);
  }
  free(mol_eligible_beads);
  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Count->Aggregate; i++) {
    vec3d com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, *System);
    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    int move[3];
    for (int dd = 0; dd < 3; dd++) {
      move[dd] = com.v[dd] / (*box).v[dd];
      if (com.v[dd] < 0) {
        move[dd]--;
      }
    }
    // move all the beads
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[bead].Position.v[dd] -= move[dd] * (*box).v[dd];
      }
    }
  } //}}}
} //}}}

// based on options, should an aggregate be used for calculations? //{{{
bool UseAggregate(SYSTEM System, AGGREGATE *Aggregate, int id,
                  AGG_PICKER agg, int *size, double *mass) {
  bool only = true;
  bool x = false;
  *size = 0;
  *mass = 0;
  for (int j = 0; j < Aggregate[id].nMolecules; j++) {
    MOLECULE *mol = &System.Molecule[Aggregate[id].Molecule[j]];
    int mtype = mol->Type;
    MOLECULETYPE *mt = &System.MoleculeType[mtype];
    if (agg.m[mtype]) {
      (*size)++;
      *mass += mt->Mass;
    }
    // if at least one unwanted molecule is present, don't use aggregate
    if (!agg.only[mtype]) {
      only = false;
      break;
    }
    // if at least one molecule isn't exluded, use aggregate
    if (!agg.x[mtype]) {
      x = true;
    }
  }
  if (*size == 0 || // no molecules remained in aggregate
      *size < agg.range[0] || *size > agg.range[1] || // -n: not in range
      !only || // -only: found molecule that weren't supposed to be in
      !x) { // -x: didn't find any un-excluded molecules
    return false;
  } else {
    return true;
  }
} //}}}
// detect -m, -x, -only, and -n options //{{{
void AggPickerOptions(const int argc, char **argv, AGG_PICKER *opt,
                      SYSTEM System) {
  COUNT *Count = &System.Count;
  // '-n' option
  opt->range[0] = 1;
  opt->range[1] = Count->Molecule;
  TwoNumbersOption(argc, argv, "-n", opt->range, 'i');
  if (opt->range[0] > opt->range[1]) {
    SwapInt(&opt->range[0], &opt->range[1]);
  }
  // '-m' option - define aggregate size as sum of those molecule types
  opt->m = calloc(Count->MoleculeType, sizeof *opt->m);
  opt->m_flag = true;
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->m, System)) {
    opt->m_flag = false;
    InitBoolArray(opt->m, Count->MoleculeType, true);
  }
  // '-only' option - use aggregates composed only of specified molecule types
  opt->only = calloc(Count->MoleculeType, sizeof *opt->only);
  opt->only_flag = true;
  if (!TypeOption(argc, argv, "-only", 'm', true, opt->only, System)) {
    opt->only_flag = false;
    InitBoolArray(opt->only, Count->MoleculeType, true);
  }
  // '-x' option - exclude aggregates composed only of specified molecule types
  opt->x = calloc(Count->MoleculeType, sizeof *opt->x);
  opt->x_flag = true;
  if (!TypeOption(argc, argv, "-x", 'm', true, opt->x, System)) {
    opt->x_flag = false;
  }
  // error checking //{{{
  // all molecule specified
  if (opt->x_flag) {
    bool overlap = true; // are all molecule types specified by -x?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (!opt->x[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("with all molecules listed, no aggregates would be detected");
      PrintErrorOption("-x");
      exit(1);
    }
  }
  // molecules specified by -m and -only do not overlap
  if (opt->m_flag && opt->only_flag) {
    bool overlap = false;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->m[i] && opt->only[i]) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      err_msg("for any aggregate to be used, at least one molecule "
              "must be specified in both options");
      PrintErrorOption("-m/-only");
      exit(1);
    }
  }
  // molecules specified by -only and -x must differ
  if (opt->only_flag && opt->x_flag) {
    bool overlap = true; // do the two array fully overlap?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->x[i] != opt->only[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("the lists of molecules must be different");
      PrintErrorOption("-x/-only");
      exit(1);
    }
  } //}}}
} //}}}

