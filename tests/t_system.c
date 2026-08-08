/*
 * SYSTEM transformation tests.
 *
 * System.c is the largest file in the codebase and every one of these
 * functions rewrites index arrays in place - bead ids, bead type ids, molecule
 * ids and the bond/angle tables that point at them all have to stay in step. A
 * mistake shows up as silently wrong topology rather than a crash, which is
 * why the assertions below are mostly invariants ("every bond index is in
 * range", "the multiset of bonds is unchanged") rather than expected values.
 *
 * Base systems come from the committed fixtures instead of being hand-built,
 * so the input to each transformation is guaranteed to satisfy whatever
 * preconditions the readers establish. The exception is SortAll, which is pure
 * enough to drive with random tables.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>

#ifndef FIXTURES_DIR
#define FIXTURES_DIR "."
#endif
#define FIX(name) FIXTURES_DIR "/" name

static pcg32_random_t rng;

// read a fixture, guaranteeing a valid starting point //{{{
static SYSTEM read_fixture(const char *path, int type) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, path, LINE);
  in.stru.type = type;
  in.coor = in.stru;
  return ReadStructure(in, false);
} //}}}

// ---- SortAll --------------------------------------------------------------

// canonical key comparison over the first (num-1) entries of a tuple //{{{
static int tuple_cmp(const int *a, const int *b, int num) {
  for (int i = 0; i < (num - 1); i++) {
    if (a[i] != b[i]) {
      if (a[i] < b[i]) {
        return -1;
      } else {
        return 1;
      }
    }
  }
  return 0;
} //}}}
// multiset of bead-id sets, invariant under the in-tuple canonical swap //{{{
/*
 * Sorting is allowed to reorder ids inside a tuple and to reorder the tuples,
 * but it must not invent, drop, or alter a bond/angle/dihedral. Comparing the
 * per-tuple sorted id sets checks exactly that without re-implementing the
 * canonicalisation rule the code uses.
 */
static void id_set(const int *tuple, int n_ids, int *out) {
  for (int i = 0; i < n_ids; i++) {
    out[i] = tuple[i];
  }
  for (int i = 0; i < n_ids; i++) { // insertion sort, n_ids <= 4
    for (int j = i + 1; j < n_ids; j++) {
      if (out[j] < out[i]) {
        SwapInt(&out[i], &out[j]);
      }
    }
  }
}
static bool same_tuple_multiset(int (*before)[5], int (*after)[5],
                                int n, int n_ids) {
  bool *used = calloc(n, sizeof *used);
  bool ok = true;
  for (int i = 0; i < n; i++) {
    int want[4];
    id_set(before[i], n_ids, want);
    bool found = false;
    for (int j = 0; j < n && !found; j++) {
      if (used[j]) {
        continue;
      }
      int have[4];
      id_set(after[j], n_ids, have);
      bool equal = true;
      for (int k = 0; k < n_ids; k++) {
        if (want[k] != have[k]) {
          equal = false;
          break;
        }
      }
      // the trailing entry carries the bond/angle type and must travel along
      if (equal && before[i][n_ids] != after[j][n_ids]) {
        equal = false;
      }
      if (equal) {
        used[j] = true;
        found = true;
      }
    }
    if (!found) {
      ok = false;
      break;
    }
  }
  free(used);
  return ok;
} //}}}
// SortAll as a property test over random tables //{{{
static void test_sort_all_property(void) {
  for (int trial = 0; trial < 300; trial++) {
    MOLECULETYPE mt;
    InitMoleculeType(&mt);
    mt.nBeads = 8;
    mt.Bead = calloc(mt.nBeads, sizeof *mt.Bead);
    mt.nBonds = pcg32RandIntInt(&rng, 1, 10);
    mt.nAngles = pcg32RandIntInt(&rng, 1, 8);
    mt.nDihedrals = pcg32RandIntInt(&rng, 1, 6);
    mt.nImpropers = pcg32RandIntInt(&rng, 1, 6);
    mt.Bond = calloc(mt.nBonds, sizeof *mt.Bond);
    mt.Angle = calloc(mt.nAngles, sizeof *mt.Angle);
    mt.Dihedral = calloc(mt.nDihedrals, sizeof *mt.Dihedral);
    mt.Improper = calloc(mt.nImpropers, sizeof *mt.Improper);
    // random tuples, with a distinguishable type value in the trailing slot
    const struct { int (*arr)[5]; int n; int n_ids; } tables[] = {
      {mt.Bond, mt.nBonds, 2},
      {mt.Angle, mt.nAngles, 3},
      {mt.Dihedral, mt.nDihedrals, 4},
      {mt.Improper, mt.nImpropers, 4},
    };
    for (size_t t = 0; t < 4; t++) {
      for (int i = 0; i < tables[t].n; i++) {
        for (int k = 0; k < tables[t].n_ids; k++) {
          tables[t].arr[i][k] = pcg32Rand0Int(&rng, mt.nBeads);
        }
        tables[t].arr[i][tables[t].n_ids] = i; // type id, must follow the tuple
      }
    }
    // snapshot the input
    int (*before[4])[5];
    for (size_t t = 0; t < 4; t++) {
      before[t] = malloc(tables[t].n * sizeof *before[t]);
      memcpy(before[t], tables[t].arr, tables[t].n * sizeof *before[t]);
    }

    SortAll(&mt);

    for (size_t t = 0; t < 4; t++) {
      int (*arr)[5] = tables[t].arr;
      int n = tables[t].n;
      int n_ids = tables[t].n_ids;
      // 1) each tuple is canonical: first id <= last id
      for (int i = 0; i < n; i++) {
        CHECK(arr[i][0] <= arr[i][n_ids - 1]);
      }
      // 2) the table is sorted on the leading ids
      for (int i = 0; i < (n - 1); i++) {
        CHECK(tuple_cmp(arr[i], arr[i + 1], n_ids + 1) <= 0);
      }
      // 3) nothing invented, dropped, or retyped
      CHECK(same_tuple_multiset(before[t], arr, n, n_ids));
    }
    // 4) sorting an already sorted table changes nothing
    int (*again[4])[5];
    for (size_t t = 0; t < 4; t++) {
      again[t] = malloc(tables[t].n * sizeof *again[t]);
      memcpy(again[t], tables[t].arr, tables[t].n * sizeof *again[t]);
    }
    SortAll(&mt);
    for (size_t t = 0; t < 4; t++) {
      CHECK(memcmp(again[t], tables[t].arr,
                   tables[t].n * sizeof *again[t]) == 0);
      free(again[t]);
      free(before[t]);
    }
    FreeMoleculeTypeEssentials(&mt);
  }
} //}}}

// ---- type identity and lookup ---------------------------------------------

// SameBeadType, including the HIGHNUM wildcard //{{{
static void test_same_bead_type(void) {
  BEADTYPE a;
  InitBeadType(&a);
  s_strcpy(a.Name, "A", BEAD_NAME);
  a.Charge = 1.0;
  a.Mass = 2.0;
  a.Radius = 0.5;
  BEADTYPE b = a;
  CHECK(SameBeadType(a, b, true));
  CHECK(SameBeadType(a, b, false));

  // a different name only matters when name checking is on
  s_strcpy(b.Name, "B", BEAD_NAME);
  CHECK(!SameBeadType(a, b, true));
  CHECK(SameBeadType(a, b, false));

  // any differing numeric property makes them distinct
  b = a;
  b.Charge = -1.0;
  CHECK(!SameBeadType(a, b, true));
  b = a;
  b.Mass = 9.0;
  CHECK(!SameBeadType(a, b, true));
  b = a;
  b.Radius = 9.0;
  CHECK(!SameBeadType(a, b, true));

  /*
   * HIGHNUM in either operand acts as "unspecified" and matches anything.
   * This is what lets a type read from a file without a mass merge with one
   * that has it.
   */
  b = a;
  b.Mass = HIGHNUM;
  CHECK(SameBeadType(a, b, true));
  CHECK(SameBeadType(b, a, true));
  b = a;
  b.Charge = HIGHNUM;
  CHECK(SameBeadType(a, b, true));
  b = a;
  b.Radius = HIGHNUM;
  CHECK(SameBeadType(a, b, true));
} //}}}
// FindBeadType / FindMoleculeName //{{{
static void test_find_types(void) {
  SYSTEM S = read_fixture(FIX("struct.vtf"), VTF_FILE);
  // the fixture defines bead types A, B and the default C
  int a = FindBeadType("A", S);
  int b = FindBeadType("B", S);
  CHECK(a != -1);
  CHECK(b != -1);
  CHECK(a != b);
  CHECK(strcmp(S.BeadType[a].Name, "A") == 0);
  // lookups are exact, not prefix or case-insensitive
  CHECK(FindBeadType("Z", S) == -1);
  CHECK(FindBeadType("a", S) == -1);
  CHECK(FindBeadType("", S) == -1);

  int m = FindMoleculeName("mol", S);
  CHECK(m != -1);
  CHECK(strcmp(S.MoleculeType[m].Name, "mol") == 0);
  CHECK(FindMoleculeName("nope", S) == -1);
  FreeSystem(&S);
} //}}}
// FillMoleculeTypeBType collects the distinct bead types of a molecule //{{{
static void test_fill_molecule_type_btype(void) {
  MOLECULETYPE mt;
  InitMoleculeType(&mt);
  mt.nBeads = 6;
  mt.Bead = malloc(mt.nBeads * sizeof *mt.Bead);
  // types 3, 1, 3, 7, 1, 3 -> distinct set {3, 1, 7} in first-seen order
  int types[6] = {3, 1, 3, 7, 1, 3};
  for (int i = 0; i < 6; i++) {
    mt.Bead[i] = types[i];
  }
  FillMoleculeTypeBType(&mt);
  CHECK(mt.nBTypes == 3);
  CHECK(mt.BType[0] == 3);
  CHECK(mt.BType[1] == 1);
  CHECK(mt.BType[2] == 7);
  free(mt.BType);
  free(mt.Bead);

  // a single-bead molecule has exactly one bead type
  MOLECULETYPE one;
  InitMoleculeType(&one);
  one.nBeads = 1;
  one.Bead = malloc(sizeof *one.Bead);
  one.Bead[0] = 5;
  FillMoleculeTypeBType(&one);
  CHECK(one.nBTypes == 1);
  CHECK(one.BType[0] == 5);
  free(one.BType);
  free(one.Bead);
} //}}}
// RenameBeadTypes / RenameMoleculeTypes disambiguate duplicates //{{{
static void test_rename_types(void) {
  SYSTEM S;
  InitSystem(&S);
  S.Count.BeadType = 0;
  // three types sharing a name, plus one that does not
  NewBeadType(&S.BeadType, &S.Count.BeadType, "X", 0, 1, 1);
  NewBeadType(&S.BeadType, &S.Count.BeadType, "X", 0, 1, 1);
  NewBeadType(&S.BeadType, &S.Count.BeadType, "Y", 0, 1, 1);
  NewBeadType(&S.BeadType, &S.Count.BeadType, "X", 0, 1, 1);
  RenameBeadTypes(&S);
  // the first keeps the name; later duplicates are numbered in order
  CHECK(strcmp(S.BeadType[0].Name, "X") == 0);
  CHECK(strcmp(S.BeadType[1].Name, "X1") == 0);
  CHECK(strcmp(S.BeadType[2].Name, "Y") == 0);
  CHECK(strcmp(S.BeadType[3].Name, "X2") == 0);
  // every name is now unique, which is the point of the pass
  for (int i = 0; i < S.Count.BeadType; i++) {
    for (int j = i + 1; j < S.Count.BeadType; j++) {
      CHECK(strcmp(S.BeadType[i].Name, S.BeadType[j].Name) != 0);
    }
  }
  FreeSystem(&S);
} //}}}

// ---- copying --------------------------------------------------------------

// CopySystem must deep-copy: mutating the copy cannot touch the original //{{{
static void test_copy_system_independence(void) {
  SYSTEM a = read_fixture(FIX("system.data"), LDATA_FILE);
  SYSTEM b = CopySystem(a);
  // same content to start with
  CHECK(a.Count.Bead == b.Count.Bead);
  CHECK(a.Count.Molecule == b.Count.Molecule);
  CHECK(a.Count.Bond == b.Count.Bond);
  // the arrays must not be shared
  CHECK(a.Bead != b.Bead);
  CHECK(a.BeadType != b.BeadType);
  CHECK(a.Molecule != b.Molecule);
  CHECK(a.MoleculeType != b.MoleculeType);

  // snapshot a few values from the original
  double pos = a.Bead[0].Position.x;
  char name[BEAD_NAME];
  s_strcpy(name, a.BeadType[0].Name, BEAD_NAME);
  int mol_bead = a.Molecule[0].Bead[0];
  int mt_bead = a.MoleculeType[0].Bead[0];
  int bond0 = a.MoleculeType[0].Bond[0][0];

  // scribble all over the copy
  b.Bead[0].Position.x = pos + 12345;
  s_strcpy(b.BeadType[0].Name, "ZZZ", BEAD_NAME);
  b.Molecule[0].Bead[0] = 999;
  b.MoleculeType[0].Bead[0] = 888;
  b.MoleculeType[0].Bond[0][0] = 777;
  b.Count.Bead = -1;

  // the original is untouched
  CHECK(a.Bead[0].Position.x == pos);
  CHECK(strcmp(a.BeadType[0].Name, name) == 0);
  CHECK(a.Molecule[0].Bead[0] == mol_bead);
  CHECK(a.MoleculeType[0].Bead[0] == mt_bead);
  CHECK(a.MoleculeType[0].Bond[0][0] == bond0);
  CHECK(a.Count.Bead != -1);

  // both must be independently freeable (ASan catches a double free here)
  FreeSystem(&b);
  CHECK(a.Bead[0].Position.x == pos);
  FreeSystem(&a);
} //}}}
// CopyMoleculeType deep-copies the topology tables //{{{
static void test_copy_molecule_type(void) {
  SYSTEM S = read_fixture(FIX("system.data"), LDATA_FILE);
  CHECK(S.Count.MoleculeType > 0);
  MOLECULETYPE orig = S.MoleculeType[0];
  MOLECULETYPE copy = CopyMoleculeType(orig);
  CHECK(copy.nBeads == orig.nBeads);
  CHECK(copy.nBonds == orig.nBonds);
  CHECK(copy.nAngles == orig.nAngles);
  CHECK(copy.Bead != orig.Bead); // not aliased
  for (int i = 0; i < orig.nBeads; i++) {
    CHECK(copy.Bead[i] == orig.Bead[i]);
  }
  if (orig.nBonds > 0) {
    CHECK(copy.Bond != orig.Bond);
    for (int i = 0; i < orig.nBonds; i++) {
      for (int k = 0; k < 3; k++) {
        CHECK(copy.Bond[i][k] == orig.Bond[i][k]);
      }
    }
    // mutating the copy leaves the original alone
    int keep = orig.Bond[0][0];
    copy.Bond[0][0] = 4242;
    CHECK(orig.Bond[0][0] == keep);
  }
  FreeMoleculeType(&copy);
  FreeSystem(&S);
} //}}}

// ---- structural transformations -------------------------------------------

// every index in the system points somewhere valid //{{{
static void check_indices_in_range(const char *label, const SYSTEM S) {
  for (int i = 0; i < S.Count.Bead; i++) {
    if (S.Bead[i].Type < 0 || S.Bead[i].Type >= S.Count.BeadType) {
      g_failures++;
      g_checks++;
      fprintf(stderr, "  FAIL %s: Bead[%d].Type = %d (BeadType = %d)\n",
              label, i, S.Bead[i].Type, S.Count.BeadType);
      break;
    }
  }
  for (int i = 0; i < S.Count.Molecule; i++) {
    CHECK(S.Molecule[i].Type >= 0 &&
          S.Molecule[i].Type < S.Count.MoleculeType);
  }
  for (int i = 0; i < S.Count.BeadCoor; i++) {
    CHECK(S.BeadCoor[i] >= 0 && S.BeadCoor[i] < S.Count.Bead);
  }
  // bond ids are molecule-relative, so they must fit inside the molecule
  for (int i = 0; i < S.Count.MoleculeType; i++) {
    MOLECULETYPE *mt = &S.MoleculeType[i];
    for (int j = 0; j < mt->nBonds; j++) {
      for (int k = 0; k < 2; k++) {
        CHECK(mt->Bond[j][k] >= 0 && mt->Bond[j][k] < mt->nBeads);
      }
    }
    for (int j = 0; j < mt->nBeads; j++) {
      CHECK(mt->Bead[j] >= 0 && mt->Bead[j] < S.Count.BeadType);
    }
  }
  // bead counts have to add up the way CheckSystem expects
  CHECK(S.Count.Unbonded + S.Count.Bonded == S.Count.Bead);
  int sum = 0;
  for (int i = 0; i < S.Count.BeadType; i++) {
    sum += S.BeadType[i].Number;
  }
  CHECK(sum == S.Count.Bead);
} //}}}
// PruneSystem drops beads absent from the timestep and reindexes the rest //{{{
/*
 * Only unbonded beads are removed: they can go without disturbing any
 * molecule, which keeps the expected result something we can state exactly.
 */
static void test_prune_system(void) {
  SYSTEM S = read_fixture(FIX("system.data"), LDATA_FILE);
  int bead_before = S.Count.Bead;
  int mol_before = S.Count.Molecule;
  int bonded_before = S.Count.Bonded;
  CHECK(S.Count.BeadCoor == bead_before);

  // collect the unbonded beads and drop the first two from the timestep
  int dropped[2];
  int n_dropped = 0;
  for (int i = 0; i < S.Count.Bead && n_dropped < 2; i++) {
    if (S.Bead[i].Molecule == -1) {
      dropped[n_dropped++] = i;
    }
  }
  CHECK(n_dropped == 2);
  for (int d = 0; d < n_dropped; d++) {
    S.Bead[dropped[d]].InTimestep = false;
  }
  // rebuild BeadCoor without them
  int keep = 0;
  for (int i = 0; i < S.Count.Bead; i++) {
    if (S.Bead[i].InTimestep) {
      S.BeadCoor[keep++] = i;
    }
  }
  S.Count.BeadCoor = keep;

  int *map = malloc(bead_before * sizeof *map);
  PruneSystem(&S, map);

  CHECK(S.Count.Bead == bead_before - 2);
  CHECK(S.Count.Molecule == mol_before); // molecules untouched
  CHECK(S.Count.Bonded == bonded_before);
  check_indices_in_range("prune", S);

  // the mapping: dropped beads are -1, kept beads map injectively into range
  CHECK(map[dropped[0]] == -1);
  CHECK(map[dropped[1]] == -1);
  int mapped = 0;
  bool *seen = calloc(S.Count.Bead, sizeof *seen);
  for (int i = 0; i < bead_before; i++) {
    if (map[i] != -1) {
      CHECK(map[i] >= 0 && map[i] < S.Count.Bead);
      if (map[i] >= 0 && map[i] < S.Count.Bead) {
        CHECK(!seen[map[i]]); // injective
        seen[map[i]] = true;
      }
      mapped++;
    }
  }
  CHECK(mapped == S.Count.Bead);
  free(seen);
  free(map);
  FreeSystem(&S);
} //}}}
// pruning a system with nothing to drop is a no-op on the counts //{{{
static void test_prune_system_keeps_all(void) {
  SYSTEM S = read_fixture(FIX("system.data"), LDATA_FILE);
  int bead = S.Count.Bead;
  int mol = S.Count.Molecule;
  int bond = S.Count.Bond;
  int *map = malloc(bead * sizeof *map);
  PruneSystem(&S, map);
  CHECK(S.Count.Bead == bead);
  CHECK(S.Count.Molecule == mol);
  CHECK(S.Count.Bond == bond);
  check_indices_in_range("prune-all", S);
  // with nothing removed the mapping is the identity
  for (int i = 0; i < bead; i++) {
    CHECK(map[i] == i);
  }
  free(map);
  FreeSystem(&S);
} //}}}
// ConcatenateSystems appends a second system with every index offset //{{{
static void test_concatenate_systems(void) {
  SYSTEM out = read_fixture(FIX("system.data"), LDATA_FILE);
  SYSTEM add = read_fixture(FIX("system.data"), LDATA_FILE);
  int bead = out.Count.Bead;
  int mol = out.Count.Molecule;
  int bonded = out.Count.Bonded;
  int unbonded = out.Count.Unbonded;
  int mol_coor = out.Count.MoleculeCoor;
  int bead_coor = out.Count.BeadCoor;
  BOX box = out.Box;
  CHECK(mol_coor == mol); // the fixture has every molecule in the timestep

  ConcatenateSystems(&out, add, box, false);

  CHECK(out.Count.Bead == 2 * bead);
  CHECK(out.Count.Molecule == 2 * mol);
  CHECK(out.Count.Bonded == 2 * bonded);
  CHECK(out.Count.Unbonded == 2 * unbonded);
  CHECK(out.Count.BeadCoor == 2 * bead_coor);
  /*
   * Every *Coor count must track the entries actually appended to its array.
   * MoleculeCoor is the one that used to be left behind, which made the
   * concatenated system claim fewer in-timestep molecules than the array held.
   */
  CHECK(out.Count.MoleculeCoor == 2 * mol_coor);
  check_indices_in_range("concat", out);
  // every listed molecule is a real molecule of the combined system
  for (int i = 0; i < out.Count.MoleculeCoor; i++) {
    CHECK(out.MoleculeCoor[i] >= 0);
    CHECK(out.MoleculeCoor[i] < out.Count.Molecule);
  }
  // the appended half points at the molecules that came from the second system
  for (int i = mol_coor; i < out.Count.MoleculeCoor; i++) {
    CHECK(out.MoleculeCoor[i] >= mol);
  }
  /*
   * InTimestep is carried across unchanged rather than recomputed. It is
   * deliberately not tied to MoleculeCoor membership here: the LAMMPS data
   * reader fills MoleculeCoor with every molecule while leaving InTimestep
   * false, so asserting that coupling would be testing the reader, not this
   * function.
   */
  for (int i = 0; i < mol; i++) {
    CHECK(out.Molecule[mol + i].InTimestep == add.Molecule[i].InTimestep);
  }
  // no molecule is listed twice
  for (int i = 0; i < out.Count.MoleculeCoor; i++) {
    for (int j = i + 1; j < out.Count.MoleculeCoor; j++) {
      CHECK(out.MoleculeCoor[i] != out.MoleculeCoor[j]);
    }
  }
  // beads from the appended system point at appended molecules
  for (int i = bead; i < out.Count.Bead; i++) {
    if (out.Bead[i].Molecule != -1) {
      CHECK(out.Bead[i].Molecule >= mol);
      CHECK(out.Bead[i].Molecule < 2 * mol);
    }
  }
  // and appended molecules own appended beads
  for (int i = mol; i < out.Count.Molecule; i++) {
    MOLECULETYPE *mt = &out.MoleculeType[out.Molecule[i].Type];
    for (int j = 0; j < mt->nBeads; j++) {
      CHECK(out.Molecule[i].Bead[j] >= bead);
      CHECK(out.Molecule[i].Bead[j] < 2 * bead);
    }
  }
  // the source system is not consumed and stays independently freeable
  CHECK(add.Count.Bead == bead);
  FreeSystem(&add);
  CHECK(out.Count.Bead == 2 * bead);
  FreeSystem(&out);
} //}}}
// concatenating with pruning collapses the duplicate types //{{{
static void test_concatenate_with_prune(void) {
  SYSTEM out = read_fixture(FIX("system.data"), LDATA_FILE);
  SYSTEM add = read_fixture(FIX("system.data"), LDATA_FILE);
  int bead = out.Count.Bead;
  int btype = out.Count.BeadType;
  BOX box = out.Box;

  ConcatenateSystems(&out, add, box, true);

  CHECK(out.Count.Bead == 2 * bead);
  // the two halves use the same bead types, so pruning must not double them
  CHECK(out.Count.BeadType == btype);
  check_indices_in_range("concat-prune", out);
  FreeSystem(&add);
  FreeSystem(&out);
} //}}}
// ChangeBoxByLow shifts coordinates by Box.Low and back //{{{
static void test_change_box_by_low(void) {
  SYSTEM S = read_fixture(FIX("system.data"), LDATA_FILE);
  S.Box.Low = (vec3d){ .v = {1.5, -2.0, 0.25} };
  // snapshot every coordinate
  vec3d *saved = malloc(S.Count.Bead * sizeof *saved);
  for (int i = 0; i < S.Count.Bead; i++) {
    saved[i] = S.Bead[i].Position;
  }
  ChangeBoxByLow(&S, 1);
  for (int i = 0; i < S.Count.BeadCoor; i++) {
    int id = S.BeadCoor[i];
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(S.Bead[id].Position.v[dd],
                  saved[id].v[dd] + S.Box.Low.v[dd], 1e-12);
    }
  }
  // the inverse restores the originals exactly
  ChangeBoxByLow(&S, -1);
  for (int i = 0; i < S.Count.BeadCoor; i++) {
    int id = S.BeadCoor[i];
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(S.Bead[id].Position.v[dd], saved[id].v[dd], 1e-12);
    }
  }
  free(saved);
  FreeSystem(&S);
} //}}}
// NewBeadType / NewMolType append without disturbing what is already there //{{{
static void test_new_types(void) {
  SYSTEM S;
  InitSystem(&S);
  S.Count.BeadType = 0;
  NewBeadType(&S.BeadType, &S.Count.BeadType, "First", 1.0, 2.0, 3.0);
  NewBeadType(&S.BeadType, &S.Count.BeadType, "Second", -1.0, 4.0, 5.0);
  CHECK(S.Count.BeadType == 2);
  CHECK(strcmp(S.BeadType[0].Name, "First") == 0);
  CHECK(S.BeadType[0].Charge == 1.0);
  CHECK(S.BeadType[0].Mass == 2.0);
  CHECK(S.BeadType[0].Radius == 3.0);
  CHECK(S.BeadType[0].Number == 0); // a new type starts empty
  CHECK(strcmp(S.BeadType[1].Name, "Second") == 0);
  CHECK(S.BeadType[1].Charge == -1.0);

  S.Count.MoleculeType = 0;
  NewMolType(&S.MoleculeType, &S.Count.MoleculeType, "m1", 4, 3, 2, 1, 0);
  CHECK(S.Count.MoleculeType == 1);
  MOLECULETYPE *mt = &S.MoleculeType[0];
  CHECK(strcmp(mt->Name, "m1") == 0);
  CHECK(mt->nBeads == 4);
  CHECK(mt->nBonds == 3);
  CHECK(mt->nAngles == 2);
  CHECK(mt->nDihedrals == 1);
  CHECK(mt->nImpropers == 0);
  CHECK(mt->Number == 1);
  // the allocated tables are zeroed
  for (int i = 0; i < mt->nBeads; i++) {
    CHECK(mt->Bead[i] == 0);
  }
  // a second type does not disturb the first
  NewMolType(&S.MoleculeType, &S.Count.MoleculeType, "m2", 2, 1, 0, 0, 0);
  CHECK(S.Count.MoleculeType == 2);
  CHECK(strcmp(S.MoleculeType[0].Name, "m1") == 0);
  CHECK(S.MoleculeType[0].nBeads == 4);
  CHECK(strcmp(S.MoleculeType[1].Name, "m2") == 0);
  /*
   * NewMolType sets Number = 1 but leaves Index unallocated, and
   * FreeMoleculeType frees Index whenever Number > 0. FillMoleculeTypeIndex
   * is the step that allocates it, so a type built this way is not safe to
   * free until it has been called - which is why it appears here rather than
   * going straight to FreeSystem.
   */
  FillMoleculeTypeIndex(&S);
  FreeSystem(&S);
} //}}}

int main(void) {
  pcg32Seed(&rng, 20260807u);
  RUN(test_sort_all_property);
  RUN(test_same_bead_type);
  RUN(test_find_types);
  RUN(test_fill_molecule_type_btype);
  RUN(test_rename_types);
  RUN(test_copy_system_independence);
  RUN(test_copy_molecule_type);
  RUN(test_prune_system);
  RUN(test_prune_system_keeps_all);
  RUN(test_concatenate_systems);
  RUN(test_concatenate_with_prune);
  RUN(test_change_box_by_low);
  RUN(test_new_types);
  return test_main_end();
}
