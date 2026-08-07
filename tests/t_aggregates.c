/*
 * Aggregate detection (DBSCAN) tests.
 *
 * EvaluateContacts() is the scientific core of the Aggregates utility and had
 * no coverage. The main test here is a differential one, in the same spirit as
 * t_traversal.c: with neighbours == 1 the DBSCAN reduces to plain connected
 * components, so a random contact graph can be partitioned independently by a
 * union-find and the two partitions compared. Any molecule placed in the wrong
 * cluster, dropped, or counted twice shows up immediately.
 *
 * For neighbours > 1 there is no such simple reference, so the structural
 * invariants are asserted instead: the core/border split, cluster coverage,
 * and the fact that every multi-molecule cluster is seeded by a core point.
 *
 * The rest covers UseAggregate()'s option truth table, FillAggregateBeads(),
 * SortAggStruct(), and the PairHash the contact counting is built on.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>

static pcg32_random_t rng;

// ---- synthetic system -----------------------------------------------------

// build n_mol molecules of one type, each with beads_per_mol beads //{{{
/*
 * EvaluateContacts only reads Count.Molecule and Molecule[].InTimestep, but
 * the aggregate helpers also need molecule types, bead lists, and masses, so
 * build a system complete enough for all of them.
 */
static SYSTEM make_system(int n_mol, int beads_per_mol) {
  SYSTEM S;
  InitSystem(&S);
  COUNT *Count = &S.Count;
  Count->BeadType = 1;
  S.BeadType = s_realloc(S.BeadType, sizeof *S.BeadType);
  InitBeadType(&S.BeadType[0]);
  s_strcpy(S.BeadType[0].Name, "A", BEAD_NAME);
  S.BeadType[0].Mass = 1;
  S.BeadType[0].Charge = 0;
  S.BeadType[0].Number = n_mol * beads_per_mol;

  Count->MoleculeType = 1;
  S.MoleculeType = s_realloc(S.MoleculeType, sizeof *S.MoleculeType);
  MOLECULETYPE *mt = &S.MoleculeType[0];
  InitMoleculeType(mt);
  s_strcpy(mt->Name, "mol", MOL_NAME);
  mt->Number = n_mol;
  mt->nBeads = beads_per_mol;
  mt->Bead = malloc(beads_per_mol * sizeof *mt->Bead);
  for (int i = 0; i < beads_per_mol; i++) {
    mt->Bead[i] = 0; // every bead is of type 0
  }
  mt->Mass = beads_per_mol * S.BeadType[0].Mass;
  mt->Charge = 0;
  mt->Index = malloc(n_mol * sizeof *mt->Index);

  Count->Bead = n_mol * beads_per_mol;
  Count->Bonded = Count->Bead;
  Count->BeadCoor = Count->Bead;
  S.Bead = s_realloc(S.Bead, Count->Bead * sizeof *S.Bead);
  S.BeadCoor = s_realloc(S.BeadCoor, Count->Bead * sizeof *S.BeadCoor);
  for (int i = 0; i < Count->Bead; i++) {
    InitBead(&S.Bead[i]);
    S.Bead[i].Type = 0;
    S.Bead[i].Molecule = i / beads_per_mol;
    S.Bead[i].InTimestep = true;
    S.BeadCoor[i] = i;
  }

  Count->Molecule = n_mol;
  Count->MoleculeCoor = n_mol;
  Count->HighestResid = n_mol - 1;
  S.Molecule = s_realloc(S.Molecule, n_mol * sizeof *S.Molecule);
  for (int i = 0; i < n_mol; i++) {
    InitMolecule(&S.Molecule[i]);
    S.Molecule[i].Type = 0;
    S.Molecule[i].Index = i;
    S.Molecule[i].InTimestep = true;
    S.Molecule[i].Aggregate = -1;
    S.Molecule[i].Bead = malloc(beads_per_mol * sizeof *S.Molecule[i].Bead);
    for (int j = 0; j < beads_per_mol; j++) {
      S.Molecule[i].Bead[j] = i * beads_per_mol + j;
    }
    mt->Index[i] = i;
  }
  return S;
} //}}}

// ---- union-find reference -------------------------------------------------

static int uf_find(int *parent, int x) { //{{{
  while (parent[x] != x) {
    parent[x] = parent[parent[x]]; // path halving
    x = parent[x];
  }
  return x;
}
static void uf_union(int *parent, int a, int b) {
  int ra = uf_find(parent, a);
  int rb = uf_find(parent, b);
  if (ra != rb) {
    parent[ra] = rb;
  }
} //}}}

// ---- helpers over the produced aggregates ---------------------------------

// map molecule -> aggregate id, checking each appears exactly once //{{{
/*
 * Returns a malloc'd array with -1 for molecules in no aggregate. Also the
 * coverage check: a molecule listed in two aggregates, or listed twice in one,
 * is a failure regardless of what the partition looks like.
 */
static int *agg_membership(const SYSTEM S, const AGGREGATE *Agg,
                           const char *label) {
  int *of = malloc(S.Count.Molecule * sizeof *of);
  for (int i = 0; i < S.Count.Molecule; i++) {
    of[i] = -1;
  }
  for (int a = 0; a < S.Count.Aggregate; a++) {
    CHECK(Agg[a].nMolecules == Agg[a].nCore + Agg[a].nBorder);
    for (int j = 0; j < Agg[a].nMolecules; j++) {
      int mol = AggGetMol(&Agg[a], j);
      if (mol < 0 || mol >= S.Count.Molecule) {
        g_failures++;
        g_checks++;
        fprintf(stderr, "  FAIL %s: aggregate %d holds molecule id %d\n",
                label, a, mol);
        continue;
      }
      if (of[mol] != -1) {
        g_failures++;
        g_checks++;
        fprintf(stderr, "  FAIL %s: molecule %d in aggregates %d and %d\n",
                label, mol, of[mol], a);
      }
      of[mol] = a;
    }
  }
  // every molecule in the timestep must have landed somewhere
  for (int i = 0; i < S.Count.Molecule; i++) {
    if (S.Molecule[i].InTimestep && of[i] == -1) {
      g_failures++;
      g_checks++;
      fprintf(stderr, "  FAIL %s: molecule %d in no aggregate\n", label, i);
    }
  }
  return of;
} //}}}

// random contact graph; edges are also recorded for the reference partition //{{{
struct graph {
  PairHash *hash;
  int *edge;  // 2*n_edges molecule ids
  int n_edge;
  int *degree; // number of distinct neighbours per molecule
};

static struct graph random_graph(int n_mol, int n_edge_target, int contacts) {
  struct graph g = {0};
  g.hash = PairHashAlloc();
  g.edge = malloc(2 * n_edge_target * sizeof *g.edge);
  g.degree = calloc(n_mol, sizeof *g.degree);
  // track which pairs already exist so degrees stay exact
  bool *seen = calloc((size_t)n_mol * n_mol, sizeof *seen);
  for (int e = 0; e < n_edge_target; e++) {
    int i = pcg32Rand0Int(&rng, n_mol);
    int j = pcg32Rand0Int(&rng, n_mol);
    if (i == j) {
      continue;
    }
    SortPairAsc(&j, &i); // PairHash wants the larger index first
    if (seen[(size_t)i * n_mol + j]) {
      continue;
    }
    seen[(size_t)i * n_mol + j] = true;
    // raise the count to exactly the threshold so the pair counts as an edge
    for (int c = 0; c < contacts; c++) {
      PairHashIncrement(g.hash, i, j);
    }
    g.edge[2 * g.n_edge] = i;
    g.edge[2 * g.n_edge + 1] = j;
    g.n_edge++;
    g.degree[i]++;
    g.degree[j]++;
  }
  free(seen);
  return g;
}
static void free_graph(struct graph *g) {
  PairHashFree(g->hash);
  free(g->edge);
  free(g->degree);
} //}}}

// ---- the differential test ------------------------------------------------

// neighbours == 1: DBSCAN degenerates to connected components //{{{
static void test_dbscan_equals_connected_components(void) {
  const int n_mol = 60;
  for (int trial = 0; trial < 200; trial++) {
    SYSTEM S = make_system(n_mol, 2);
    // vary the density so trials span "mostly isolated" to "one big blob"
    int n_edge = pcg32Rand0Int(&rng, 120);
    struct graph g = random_graph(n_mol, n_edge, 1);

    AGGREGATE *Agg = nullptr;
    InitAggregate(S, &Agg);
    EvaluateContacts(Agg, &S, 1, 1, g.hash);

    int *of = agg_membership(S, Agg, "components");
    // reference partition over the same edge set
    int *parent = malloc(n_mol * sizeof *parent);
    for (int i = 0; i < n_mol; i++) {
      parent[i] = i;
    }
    for (int e = 0; e < g.n_edge; e++) {
      uf_union(parent, g.edge[2 * e], g.edge[2 * e + 1]);
    }
    // same aggregate  <=>  same connected component
    for (int i = 0; i < n_mol; i++) {
      for (int j = i + 1; j < n_mol; j++) {
        bool same_agg = (of[i] == of[j]);
        bool same_comp = (uf_find(parent, i) == uf_find(parent, j));
        if (same_agg != same_comp) {
          g_failures++;
          fprintf(stderr, "  FAIL components: molecules %d,%d agg %d,%d "
                  "but components %d,%d\n", i, j, of[i], of[j],
                  uf_find(parent, i), uf_find(parent, j));
        }
        g_checks++;
      }
    }
    // Molecule[].Aggregate must agree with the Aggregate array
    for (int i = 0; i < n_mol; i++) {
      CHECK(S.Molecule[i].Aggregate == of[i]);
    }
    free(parent);
    free(of);
    FreeAggregate(S.Count, Agg);
    free_graph(&g);
    FreeSystem(&S);
  }
} //}}}

// neighbours > 1: structural invariants of the core/border split //{{{
static void test_dbscan_core_border_invariants(void) {
  const int n_mol = 50;
  for (int trial = 0; trial < 120; trial++) {
    int neighbours = pcg32RandIntInt(&rng, 2, 4);
    SYSTEM S = make_system(n_mol, 1);
    struct graph g = random_graph(n_mol, pcg32Rand0Int(&rng, 200), 1);

    AGGREGATE *Agg = nullptr;
    InitAggregate(S, &Agg);
    EvaluateContacts(Agg, &S, 1, neighbours, g.hash);

    int *of = agg_membership(S, Agg, "core/border");
    for (int a = 0; a < S.Count.Aggregate; a++) {
      if (Agg[a].nMolecules == 1) {
        /*
         * Noise molecules become singleton aggregates and are filed under
         * Core even when their degree is below the threshold; that is the
         * documented shape of the output, pinned in its own test below.
         */
        continue;
      }
      // a cluster is always seeded by a core point
      CHECK(Agg[a].nCore >= 1);
      // Core holds exactly the molecules meeting the neighbour threshold
      for (int j = 0; j < Agg[a].nCore; j++) {
        CHECK(g.degree[Agg[a].Core[j]] >= neighbours);
      }
      for (int j = 0; j < Agg[a].nBorder; j++) {
        CHECK(g.degree[Agg[a].Border[j]] < neighbours);
      }
      // every border point is reachable from a core point of the same cluster
      for (int j = 0; j < Agg[a].nBorder; j++) {
        int b = Agg[a].Border[j];
        bool touches_core = false;
        for (int e = 0; e < g.n_edge && !touches_core; e++) {
          int u = g.edge[2 * e];
          int v = g.edge[2 * e + 1];
          int other = -1;
          if (u == b) {
            other = v;
          } else if (v == b) {
            other = u;
          }
          if (other >= 0 && of[other] == a &&
              g.degree[other] >= neighbours) {
            touches_core = true;
          }
        }
        CHECK(touches_core);
      }
    }
    free(of);
    FreeAggregate(S.Count, Agg);
    free_graph(&g);
    FreeSystem(&S);
  }
} //}}}

// molecules below the threshold and unconnected ones become singletons //{{{
static void test_dbscan_singletons(void) {
  SYSTEM S = make_system(4, 1);
  PairHash *h = PairHashAlloc();
  // a single contact between molecules 0 and 1; 2 and 3 stay isolated
  PairHashIncrement(h, 1, 0);
  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  // neighbours == 2 puts 0 and 1 below the core threshold as well
  EvaluateContacts(Agg, &S, 1, 2, h);
  CHECK(S.Count.Aggregate == 4); // every molecule alone
  for (int a = 0; a < S.Count.Aggregate; a++) {
    CHECK(Agg[a].nMolecules == 1);
    CHECK(Agg[a].nCore == 1);   // singletons are filed as core
    CHECK(Agg[a].nBorder == 0);
  }
  int *of = agg_membership(S, Agg, "singletons");
  free(of);
  FreeAggregate(S.Count, Agg);
  PairHashFree(h);
  FreeSystem(&S);
} //}}}

// a fully connected chain collapses into one aggregate //{{{
static void test_dbscan_single_cluster(void) {
  const int n_mol = 12;
  SYSTEM S = make_system(n_mol, 3);
  PairHash *h = PairHashAlloc();
  for (int i = 1; i < n_mol; i++) {
    PairHashIncrement(h, i, i - 1);
  }
  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  EvaluateContacts(Agg, &S, 1, 1, h);
  CHECK(S.Count.Aggregate == 1);
  CHECK(Agg[0].nMolecules == n_mol);
  int *of = agg_membership(S, Agg, "chain");
  for (int i = 0; i < n_mol; i++) {
    CHECK(of[i] == 0);
  }
  free(of);
  FreeAggregate(S.Count, Agg);
  PairHashFree(h);
  FreeSystem(&S);
} //}}}

// the contact-count threshold is inclusive //{{{
/*
 * A pair with exactly 'contacts' bead contacts is an edge; one contact fewer
 * is not. An off-by-one here silently changes every aggregate size the
 * Aggregates utility reports, so pin both sides of the boundary.
 */
static void test_dbscan_contact_threshold(void) {
  for (int contacts = 1; contacts <= 4; contacts++) {
    // exactly at the threshold -> one aggregate of two molecules
    SYSTEM S = make_system(2, 1);
    PairHash *h = PairHashAlloc();
    for (int c = 0; c < contacts; c++) {
      PairHashIncrement(h, 1, 0);
    }
    AGGREGATE *Agg = nullptr;
    InitAggregate(S, &Agg);
    EvaluateContacts(Agg, &S, contacts, 1, h);
    CHECK(S.Count.Aggregate == 1);
    CHECK(Agg[0].nMolecules == 2);
    FreeAggregate(S.Count, Agg);
    PairHashFree(h);
    FreeSystem(&S);

    // one contact short -> two singletons
    SYSTEM S2 = make_system(2, 1);
    PairHash *h2 = PairHashAlloc();
    for (int c = 0; c < (contacts - 1); c++) {
      PairHashIncrement(h2, 1, 0);
    }
    AGGREGATE *Agg2 = nullptr;
    InitAggregate(S2, &Agg2);
    EvaluateContacts(Agg2, &S2, contacts, 1, h2);
    CHECK(S2.Count.Aggregate == 2);
    FreeAggregate(S2.Count, Agg2);
    PairHashFree(h2);
    FreeSystem(&S2);
  }
} //}}}

// molecules outside the timestep never appear in an aggregate //{{{
static void test_dbscan_skips_absent_molecules(void) {
  const int n_mol = 10;
  SYSTEM S = make_system(n_mol, 1);
  PairHash *h = PairHashAlloc();
  // chain them all together, then drop the odd-numbered molecules
  for (int i = 1; i < n_mol; i++) {
    PairHashIncrement(h, i, i - 1);
  }
  for (int i = 1; i < n_mol; i += 2) {
    S.Molecule[i].InTimestep = false;
  }
  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  EvaluateContacts(Agg, &S, 1, 1, h);
  int *of = agg_membership(S, Agg, "absent");
  for (int i = 0; i < n_mol; i++) {
    if (!S.Molecule[i].InTimestep) {
      CHECK(of[i] == -1);
    } else {
      CHECK(of[i] != -1);
    }
  }
  // with every other molecule gone the chain falls apart into singletons
  CHECK(S.Count.Aggregate == 5);
  free(of);
  FreeAggregate(S.Count, Agg);
  PairHashFree(h);
  FreeSystem(&S);
} //}}}

// ---- aggregate bookkeeping ------------------------------------------------

// FillAggregateBeads: bead count and membership //{{{
static void test_fill_aggregate_beads(void) {
  const int n_mol = 6;
  const int per_mol = 4;
  SYSTEM S = make_system(n_mol, per_mol);
  PairHash *h = PairHashAlloc();
  // two clusters: {0,1,2} and {3,4,5}
  PairHashIncrement(h, 1, 0);
  PairHashIncrement(h, 2, 1);
  PairHashIncrement(h, 4, 3);
  PairHashIncrement(h, 5, 4);
  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  EvaluateContacts(Agg, &S, 1, 1, h);
  CHECK(S.Count.Aggregate == 2);
  FillAggregateBeads(Agg, S);
  for (int a = 0; a < S.Count.Aggregate; a++) {
    CHECK(Agg[a].nBeads == Agg[a].nMolecules * per_mol);
    // every listed bead must belong to a molecule of this aggregate
    for (int j = 0; j < Agg[a].nBeads; j++) {
      int bead = Agg[a].Bead[j];
      CHECK(bead >= 0 && bead < S.Count.Bead);
      if (bead >= 0 && bead < S.Count.Bead) {
        int mol = S.Bead[bead].Molecule;
        CHECK(S.Molecule[mol].Aggregate == a);
      }
    }
  }
  FreeAggregate(S.Count, Agg);
  PairHashFree(h);
  FreeSystem(&S);
} //}}}

// SortAggStruct: ordered by first member, including border-only aggregates //{{{
/*
 * The sort key falls back to Border[0] when an aggregate has no core
 * molecules, a branch nothing else exercises.
 */
static void test_sort_agg_struct(void) {
  SYSTEM S = make_system(6, 1);
  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  S.Count.Aggregate = 3;
  // deliberately out of order: keys 4, 0, 2
  // aggregate 0: core {4, 5}
  Agg[0].nCore = 2;
  Agg[0].nBorder = 0;
  Agg[0].nMolecules = 2;
  Agg[0].Core = s_realloc(Agg[0].Core, 2 * sizeof *Agg[0].Core);
  Agg[0].Core[0] = 4;
  Agg[0].Core[1] = 5;
  // aggregate 1: core {0, 1}
  Agg[1].nCore = 2;
  Agg[1].nBorder = 0;
  Agg[1].nMolecules = 2;
  Agg[1].Core = s_realloc(Agg[1].Core, 2 * sizeof *Agg[1].Core);
  Agg[1].Core[0] = 0;
  Agg[1].Core[1] = 1;
  // aggregate 2: border only, key comes from Border[0]
  Agg[2].nCore = 0;
  Agg[2].nBorder = 2;
  Agg[2].nMolecules = 2;
  Agg[2].Border = s_realloc(Agg[2].Border, 2 * sizeof *Agg[2].Border);
  Agg[2].Border[0] = 2;
  Agg[2].Border[1] = 3;

  SortAggStruct(Agg, S);
  CHECK(Agg[0].nCore == 2 && Agg[0].Core[0] == 0);
  CHECK(Agg[1].nCore == 0 && Agg[1].nBorder == 2);
  CHECK(Agg[1].nBorder == 2 && Agg[1].Border[0] == 2);
  CHECK(Agg[2].nCore == 2 && Agg[2].Core[0] == 4);
  // the counts must travel with the arrays
  for (int a = 0; a < 3; a++) {
    CHECK(Agg[a].nMolecules == Agg[a].nCore + Agg[a].nBorder);
  }
  FreeAggregate(S.Count, Agg);
  FreeSystem(&S);
} //}}}

// UseAggregate: the -m / -x / -only / -n truth table //{{{
/*
 * Two molecule types, so each option can include one and exclude the other.
 */
static void test_use_aggregate(void) {
  SYSTEM S = make_system(3, 2);
  // add a second molecule type and make molecule 2 use it
  S.Count.MoleculeType = 2;
  S.MoleculeType = s_realloc(S.MoleculeType,
                             2 * sizeof *S.MoleculeType);
  MOLECULETYPE *mt = &S.MoleculeType[1];
  InitMoleculeType(mt);
  s_strcpy(mt->Name, "other", MOL_NAME);
  mt->Number = 1;
  mt->nBeads = 2;
  mt->Bead = malloc(2 * sizeof *mt->Bead);
  mt->Bead[0] = 0;
  mt->Bead[1] = 0;
  mt->Mass = 5;
  mt->Index = malloc(sizeof *mt->Index);
  mt->Index[0] = 2;
  S.MoleculeType[0].Number = 2;
  S.Molecule[2].Type = 1;

  AGGREGATE *Agg = nullptr;
  InitAggregate(S, &Agg);
  S.Count.Aggregate = 1;
  // one aggregate holding all three molecules: types {0, 0, 1}
  Agg[0].nCore = 3;
  Agg[0].nBorder = 0;
  Agg[0].nMolecules = 3;
  Agg[0].Core = s_realloc(Agg[0].Core, 3 * sizeof *Agg[0].Core);
  for (int i = 0; i < 3; i++) {
    Agg[0].Core[i] = i;
  }

  bool m[2], only[2], x[2];
  AGG_PICKER p = {.m = m, .only = only, .x = x};
  int size = 0;
  double mass = 0;

  // everything allowed: size counts all three molecules
  m[0] = m[1] = true;
  only[0] = only[1] = true;
  x[0] = x[1] = false;
  p.range[0] = 1;
  p.range[1] = 3;
  CHECK(UseAggregate(S, Agg, 0, p, &size, &mass));
  CHECK(size == 3);
  CHECK_CLOSE(mass, 2 * S.MoleculeType[0].Mass + 5, 1e-9);

  // -m selecting only type 0: size counts two molecules, mass follows
  m[0] = true;
  m[1] = false;
  CHECK(UseAggregate(S, Agg, 0, p, &size, &mass));
  CHECK(size == 2);
  CHECK_CLOSE(mass, 2 * S.MoleculeType[0].Mass, 1e-9);

  // -n range excluding the resulting size
  p.range[0] = 3;
  p.range[1] = 3;
  CHECK(!UseAggregate(S, Agg, 0, p, &size, &mass));
  p.range[0] = 1;
  p.range[1] = 3;

  // -only listing just type 0: the type-1 molecule disqualifies the aggregate
  m[0] = m[1] = true;
  only[0] = true;
  only[1] = false;
  CHECK(!UseAggregate(S, Agg, 0, p, &size, &mass));
  only[0] = only[1] = true;

  // -x excluding every type present: nothing un-excluded remains
  x[0] = x[1] = true;
  CHECK(!UseAggregate(S, Agg, 0, p, &size, &mass));
  // excluding only one of the two types still leaves the aggregate usable
  x[1] = false;
  CHECK(UseAggregate(S, Agg, 0, p, &size, &mass));
  x[0] = x[1] = false;

  // -m selecting nothing: size stays 0, so the aggregate is unusable
  m[0] = m[1] = false;
  CHECK(!UseAggregate(S, Agg, 0, p, &size, &mass));
  CHECK(size == 0);

  FreeAggregate(S.Count, Agg);
  FreeSystem(&S);
} //}}}

// ---- PairHash -------------------------------------------------------------

// key encoding round trip and the 255 saturation cap //{{{
static void test_pair_hash(void) {
  // key round trip, including indices past the sign bit of an int32
  const int pairs[][2] = {
    {1, 0}, {5, 3}, {100000, 99999}, {2000000000, 1999999999},
  };
  for (size_t t = 0; t < sizeof pairs / sizeof *pairs; t++) {
    uint64_t key = PairHashKey(pairs[t][0], pairs[t][1]);
    CHECK(PairHash_mol_i(key) == pairs[t][0]);
    CHECK(PairHash_mol_j(key) == pairs[t][1]);
  }
  // distinct pairs must produce distinct keys
  CHECK(PairHashKey(3, 5) != PairHashKey(5, 3));

  PairHash *h = PairHashAlloc();
  // the count saturates at 255 instead of wrapping around the uint8_t
  for (int i = 0; i < 300; i++) {
    PairHashIncrement(h, 9, 2);
  }
  khiter_t it = kh_get(contacts, h, PairHashKey(9, 2));
  CHECK(it != kh_end(h));
  if (it != kh_end(h)) {
    CHECK(kh_val(h, it) == 255);
  }
  // a pair that was never incremented is absent
  CHECK(kh_get(contacts, h, PairHashKey(9, 3)) == kh_end(h));
  // counts are kept per pair, not shared
  for (int i = 0; i < 7; i++) {
    PairHashIncrement(h, 4, 1);
  }
  it = kh_get(contacts, h, PairHashKey(4, 1));
  CHECK(it != kh_end(h));
  if (it != kh_end(h)) {
    CHECK(kh_val(h, it) == 7);
  }
  PairHashFree(h);
} //}}}

int main(void) {
  pcg32Seed(&rng, 20260807u);
  RUN(test_pair_hash);
  RUN(test_dbscan_equals_connected_components);
  RUN(test_dbscan_core_border_invariants);
  RUN(test_dbscan_singletons);
  RUN(test_dbscan_single_cluster);
  RUN(test_dbscan_contact_threshold);
  RUN(test_dbscan_skips_absent_molecules);
  RUN(test_fill_aggregate_beads);
  RUN(test_sort_agg_struct);
  RUN(test_use_aggregate);
  return test_main_end();
}
