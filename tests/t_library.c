/*
 * Molecule-library tests.
 *
 * ReadLibrary.c is the largest unit in the codebase with no coverage, and it
 * is all parsing plus index arithmetic - the combination that fails quietly.
 * AddToSystem builds whole systems out of it, so a mis-mapped bond type or a
 * counterion added the wrong number of times produces a plausible-looking file
 * rather than an error.
 *
 * Two fixture sets are used, deliberately:
 *   - tests/fixtures/library/ is a tiny purpose-built library. Exact counts
 *     and parameter values are asserted against it, so it must stay stable.
 *   - Examples/library/ is the real library shipped with the project. Only
 *     structural invariants are asserted there, so the examples can grow
 *     without breaking these tests while still being checked for parseability.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>

#ifndef FIXTURES_DIR
#define FIXTURES_DIR "."
#endif
#ifndef EXAMPLES_DIR
#define EXAMPLES_DIR "."
#endif
#define LIB_DIR FIXTURES_DIR "/library"
#define EXAMPLE_LIB EXAMPLES_DIR "/library"

// complete a library's system before freeing it //{{{
/*
 * ReadLibraryMolecule() raises
 * BeadType[].Number and MoleculeType[].Number but does not allocate the
 * matching Index arrays, while FreeSystem() frees Index whenever Number > 0.
 * FillSystemNonessentials() is the step that allocates them, and it is what
 * AddToSystem calls after building its molecules, so mirror that sequence
 * here. Freeing a library straight after adding a molecule is not safe.
 */
static void finish_and_free(LIBRARY *lib) {
  FillSystemNonessentials(&lib->System, true);
  FreeLibrary(lib);
} //}}}

// ---- list parsing ---------------------------------------------------------

// bead types, bond types, angle types and interactions from the list files //{{{
static void test_read_library(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;

  // bead types come from list_parameters.txt, in file order
  CHECK(S->Count.BeadType == 3);
  CHECK(strcmp(S->BeadType[0].Name, "SOL") == 0);
  CHECK(strcmp(S->BeadType[1].Name, "MID") == 0);
  CHECK(strcmp(S->BeadType[2].Name, "ION") == 0);
  CHECK_CLOSE(S->BeadType[0].Mass, 1.0, 1e-9);
  CHECK_CLOSE(S->BeadType[1].Mass, 2.0, 1e-9);
  CHECK_CLOSE(S->BeadType[1].Charge, 0.5, 1e-9);
  CHECK_CLOSE(S->BeadType[2].Charge, -1.0, 1e-9);
  // the bead radius is taken from the Rc column
  CHECK_CLOSE(S->BeadType[1].Radius, 0.9, 1e-9);

  // bond types: k is stored doubled, r0 as written
  CHECK(S->Count.BondType == 3);
  CHECK(lib.n_bond_ids == 3);
  CHECK(strcmp(lib.bond_id[0].id, "b01") == 0);
  CHECK(strcmp(lib.bond_id[1].id, "b02") == 0);
  CHECK(strcmp(lib.bond_id[2].id, "b03") == 0);
  CHECK(lib.bond_id[0].index == 0);
  CHECK(lib.bond_id[1].index == 1);
  CHECK(lib.bond_id[2].index == 2);
  CHECK_CLOSE(S->BondType[0].a, 2 * 10.0, 1e-9);
  CHECK_CLOSE(S->BondType[0].b, 0.25, 1e-9); // 0.25d0
  CHECK_CLOSE(S->BondType[1].a, 2 * 20.0, 1e-9);
  CHECK_CLOSE(S->BondType[1].b, 0.5, 1e-9);
  CHECK_CLOSE(S->BondType[2].b, 3.0, 1e-9);

  // angle types, likewise doubled
  CHECK(S->Count.AngleType == 2);
  CHECK(lib.n_angle_ids == 2);
  CHECK(strcmp(lib.angle_id[0].id, "a01") == 0);
  CHECK_CLOSE(S->AngleType[0].a, 2 * 1.5, 1e-9);
  CHECK_CLOSE(S->AngleType[0].b, 180.0, 1e-9);
  CHECK_CLOSE(S->AngleType[1].b, 90.0, 1e-9);

  /*
   * Interactions are the self-interactions (one per bead type, taken from the
   * A/Rc columns of list_parameters.txt) followed by the cross terms.
   */
  CHECK(lib.n_inter == 3 + 2);
  // self terms first, in bead-type order
  for (int i = 0; i < 3; i++) {
    CHECK(strcmp(lib.inter[i].name1, lib.inter[i].name2) == 0);
    CHECK(strcmp(lib.inter[i].name1, S->BeadType[i].Name) == 0);
  }
  CHECK_CLOSE(lib.inter[0].A, 25.0, 1e-9);
  CHECK_CLOSE(lib.inter[0].Rc, 1.0, 1e-9);
  CHECK_CLOSE(lib.inter[1].A, 20.0, 1e-9);
  // then the cross terms, in file order
  CHECK(strcmp(lib.inter[3].name1, "SOL") == 0);
  CHECK(strcmp(lib.inter[3].name2, "MID") == 0);
  CHECK_CLOSE(lib.inter[3].A, 15.0, 1e-9);
  CHECK_CLOSE(lib.inter[4].A, 18.0, 1e-9);
  // every entry gets the same fixed friction coefficient
  for (int i = 0; i < lib.n_inter; i++) {
    CHECK_CLOSE(lib.inter[i].gamma, 4.5, 1e-9);
  }
  FreeLibrary(&lib);
} //}}}
// LibraryMoleculeInfo reads n_beads and the counterion name //{{{
static void test_molecule_info(void) {
  // no counterion: the two bead counts agree and the cion list is empty
  LIB_MOL_INFO sol = LibraryMoleculeInfo(LIB_DIR, "sol");
  CHECK(sol.n_beads == 1);
  CHECK(sol.n_beads_total == 1);
  CHECK(sol.n_cion == 0);
  CHECK(!sol.bilayer);
  CHECK_CLOSE(sol.M_w, 1.0, 1e-9);

  LIB_MOL_INFO tri = LibraryMoleculeInfo(LIB_DIR, "tri");
  CHECK(tri.n_beads == 3);
  CHECK(tri.n_beads_total == 3);
  CHECK(tri.n_cion == 0);
  CHECK(tri.bilayer);

  /*
   * 'pair' has two beads of its own plus a one-bead counterion. n_beads is the
   * molecule's own, n_beads_total counts the counterion too; confusing the two
   * is the mistake this pins.
   */
  LIB_MOL_INFO pair = LibraryMoleculeInfo(LIB_DIR, "pair");
  CHECK(pair.n_beads == 2);
  CHECK(pair.n_beads_total == 3);
  CHECK(pair.n_cion == 1);
  CHECK(strcmp(pair.cion[0].name, "ion") == 0);
  CHECK(pair.cion[0].count == 1);

  // an unknown molecule is reported as -1 rather than 0
  LIB_MOL_INFO none = LibraryMoleculeInfo(LIB_DIR, "does_not_exist");
  CHECK(none.n_beads == -1);
  CHECK(none.n_cion == 0);
  // matching is exact
  CHECK(LibraryMoleculeInfo(LIB_DIR, "SOL").n_beads == -1);
  CHECK(LibraryMoleculeInfo(LIB_DIR, "so").n_beads == -1);
} //}}}
// a trailing slash on the directory is accepted //{{{
static void test_trailing_slash(void) {
  LIB_MOL_INFO a = LibraryMoleculeInfo(LIB_DIR, "tri");
  LIB_MOL_INFO b = LibraryMoleculeInfo(LIB_DIR "/", "tri");
  CHECK(a.n_beads == b.n_beads);
  CHECK(b.n_beads == 3);
} //}}}

// ---- molecule construction ------------------------------------------------

// a one-bead molecule becomes free beads, not a molecule type //{{{
/*
 * The n_beads == 1 early return in ReadLibraryMolecule is a separate code
 * path with its own bookkeeping, and it is the one used for solvent.
 */
static void test_single_bead_molecule(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  CHECK(S->Count.Bead == 0);

  ReadLibraryMolecule(LIB_DIR, "sol", 5, true, &lib);
  CHECK(S->Count.Bead == 5);
  CHECK(S->Count.Unbonded == 5);
  CHECK(S->Count.Bonded == 0);
  CHECK(S->Count.Molecule == 0);     // no molecule is created
  CHECK(S->Count.MoleculeType == 0);
  int sol = FindBeadType("SOL", *S);
  CHECK(sol != -1);
  CHECK(S->BeadType[sol].Number == 5);
  for (int i = 0; i < S->Count.Bead; i++) {
    CHECK(S->Bead[i].Type == sol);
    CHECK(S->Bead[i].Molecule == -1); // unbonded
    CHECK(S->Bead[i].InTimestep);
  }
  finish_and_free(&lib);
} //}}}
// a multi-bead molecule builds a molecule type with mapped topology //{{{
static void test_chain_molecule(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  const int n_mols = 4;
  ReadLibraryMolecule(LIB_DIR, "tri", n_mols, true, &lib);

  CHECK(S->Count.MoleculeType == 1);
  MOLECULETYPE *mt = &S->MoleculeType[0];
  CHECK(strcmp(mt->Name, "tri") == 0);
  CHECK(mt->Number == n_mols);
  CHECK(mt->nBeads == 3);
  CHECK(mt->nBonds == 2);
  CHECK(mt->nAngles == 1);
  CHECK(mt->Named);

  CHECK(S->Count.Molecule == n_mols);
  CHECK(S->Count.Bead == n_mols * 3);
  CHECK(S->Count.Bonded == n_mols * 3);
  CHECK(S->Count.Unbonded == 0);

  // bead sequence SOL-MID-SOL
  int sol = FindBeadType("SOL", *S);
  int mid = FindBeadType("MID", *S);
  CHECK(mt->Bead[0] == sol);
  CHECK(mt->Bead[1] == mid);
  CHECK(mt->Bead[2] == sol);
  CHECK(S->BeadType[sol].Number == 2 * n_mols);
  CHECK(S->BeadType[mid].Number == n_mols);

  // bond bead indices are converted from the file's 1-based to 0-based
  CHECK(mt->Bond[0][0] == 0 && mt->Bond[0][1] == 1);
  CHECK(mt->Bond[1][0] == 1 && mt->Bond[1][1] == 2);
  // and the bond ids are resolved to BondType indices
  CHECK(mt->Bond[0][2] == 0); // b01
  CHECK(mt->Bond[1][2] == 1); // b02
  // angle likewise
  CHECK(mt->Angle[0][0] == 0);
  CHECK(mt->Angle[0][1] == 1);
  CHECK(mt->Angle[0][2] == 2);
  CHECK(mt->Angle[0][3] == 0); // a01

  // each molecule owns a contiguous run of beads, in order
  for (int m = 0; m < n_mols; m++) {
    MOLECULE *mol = &S->Molecule[m];
    CHECK(mol->Type == 0);
    CHECK(mol->InTimestep);
    for (int b = 0; b < mt->nBeads; b++) {
      CHECK(mol->Bead[b] == m * 3 + b);
      CHECK(S->Bead[mol->Bead[b]].Molecule == m);
      CHECK(S->Bead[mol->Bead[b]].Type == mt->Bead[b]);
    }
  }
  finish_and_free(&lib);
} //}}}
// repeated calls accumulate rather than overwrite //{{{
static void test_multiple_molecule_types(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  ReadLibraryMolecule(LIB_DIR, "tri", 2, true, &lib);
  ReadLibraryMolecule(LIB_DIR, "sol", 3, true, &lib);
  ReadLibraryMolecule(LIB_DIR, "tri", 1, true, &lib);

  // the second 'tri' call adds a second molecule type, it does not merge
  CHECK(S->Count.MoleculeType == 2);
  CHECK(S->MoleculeType[0].Number == 2);
  CHECK(S->MoleculeType[1].Number == 1);
  CHECK(S->Count.Molecule == 3);
  CHECK(S->Count.Bead == 2 * 3 + 3 + 1 * 3);
  CHECK(S->Count.Unbonded == 3);
  CHECK(S->Count.Bonded == 2 * 3 + 1 * 3);

  // every molecule's beads still point back at it
  for (int m = 0; m < S->Count.Molecule; m++) {
    MOLECULE *mol = &S->Molecule[m];
    MOLECULETYPE *mt = &S->MoleculeType[mol->Type];
    for (int b = 0; b < mt->nBeads; b++) {
      CHECK(mol->Bead[b] >= 0 && mol->Bead[b] < S->Count.Bead);
      CHECK(S->Bead[mol->Bead[b]].Molecule == m);
    }
  }
  finish_and_free(&lib);
} //}}}
// a counterion is added as its own species, not appended to its parent //{{{
/*
 * The layout this pins is the one every "first bead"/"last bead" default in the
 * codebase depends on. While a counterion was appended to its parent, the last
 * bead of a molecule was the counterion rather than its tail, and a monoatomic
 * one sits at the parent's origin - so anything anchoring on it placed the
 * molecule upside down while looking entirely plausible.
 */
static void test_molecule_with_counterion(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  const int n_mols = 3;
  ReadLibraryMolecule(LIB_DIR, "pair", n_mols, true, &lib);

  // the parent keeps its own two beads and nothing else
  CHECK(S->Count.MoleculeType == 1);
  MOLECULETYPE *mt = &S->MoleculeType[0];
  CHECK(strcmp(mt->Name, "pair") == 0);
  CHECK(mt->nBeads == 2);
  CHECK(mt->Number == n_mols);
  CHECK(mt->nBonds == 1);

  int mid = FindBeadType("MID", *S);
  int ion = FindBeadType("ION", *S);
  CHECK(mt->Bead[0] == mid);
  CHECK(mt->Bead[1] == mid);
  // the last bead of the molecule is the molecule's own, not the counterion
  CHECK(mt->Bead[mt->nBeads - 1] != ion);

  // the monoatomic counterion became free beads, one per parent molecule
  CHECK(S->Count.Molecule == n_mols);
  CHECK(S->Count.Bead == n_mols * 3);
  CHECK(S->Count.Bonded == n_mols * 2);
  CHECK(S->Count.Unbonded == n_mols);
  CHECK(S->BeadType[mid].Number == 2 * n_mols);
  CHECK(S->BeadType[ion].Number == n_mols);

  CHECK(mt->Bond[0][0] == 0 && mt->Bond[0][1] == 1);
  CHECK(mt->Bond[0][2] == 1); // b02

  // every ION bead is free, every MID bead belongs to a molecule
  for (int i = 0; i < S->Count.Bead; i++) {
    if (S->Bead[i].Type == ion) {
      CHECK(S->Bead[i].Molecule == -1);
    } else {
      CHECK(S->Bead[i].Molecule >= 0);
      CHECK(S->Bead[i].Molecule < n_mols);
    }
  }
  finish_and_free(&lib);
} //}}}
// with_cion=false leaves the counterion for the caller to place //{{{
/*
 * The bilayer workflow needs this: a counterion added alongside its parent
 * would inherit the leaflet's placement constraint and start inside the
 * hydrophobic core, when it belongs in the water.
 */
static void test_molecule_without_counterion(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  ReadLibraryMolecule(LIB_DIR, "pair", 3, false, &lib);

  CHECK(S->Count.MoleculeType == 1);
  CHECK(S->MoleculeType[0].nBeads == 2);
  CHECK(S->Count.Bead == 3 * 2); // no counterion beads at all
  CHECK(S->Count.Unbonded == 0);
  CHECK(S->BeadType[FindBeadType("ION", *S)].Number == 0);
  finish_and_free(&lib);
} //}}}
// a counterion with a count above one is added that many times //{{{
static void test_counterion_count(void) {
  LIB_MOL_INFO info = LibraryMoleculeInfo(LIB_DIR, "pair");
  CHECK(info.n_cion == 1);
  // the fixture declares no count, which means one
  CHECK(info.cion[0].count == 1);
  CHECK(info.n_beads_total == info.n_beads + info.cion[0].count * 1);
} //}}}

// ---- potential table ------------------------------------------------------

// FillPotFromLibrary is symmetric and falls back for unlisted pairs //{{{
static void test_fill_pot(void) {
  LIBRARY lib = ReadLibrary(LIB_DIR);
  SYSTEM *S = &lib.System;
  int n = S->Count.BeadType;
  ArrNDd *pot = CreateArr3Dd(n, n, 3);
  FillPotFromLibrary(&lib, S, pot);

  int sol = FindBeadType("SOL", *S);
  int mid = FindBeadType("MID", *S);
  int ion = FindBeadType("ION", *S);

  // self terms come from the parameter file
  CHECK_CLOSE(GetArr3D(pot, sol, sol, 0), 25.0, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, mid, mid, 0), 20.0, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, ion, ion, 0), 30.0, 1e-9);
  // cross terms from the cross-interaction file
  CHECK_CLOSE(GetArr3D(pot, sol, mid, 0), 15.0, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, sol, mid, 1), 0.95, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, mid, ion, 0), 18.0, 1e-9);
  /*
   * SOL-ION is absent from the cross list, so the hard-coded defaults apply.
   * Pinned because a silently defaulted pair is easy to miss in a real run.
   */
  CHECK_CLOSE(GetArr3D(pot, sol, ion, 0), 25.0, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, sol, ion, 1), 1.0, 1e-9);
  CHECK_CLOSE(GetArr3D(pot, sol, ion, 2), 4.5, 1e-9);

  // the table must be symmetric in its two type indices
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < 3; k++) {
        CHECK_CLOSE(GetArr3D(pot, i, j, k), GetArr3D(pot, j, i, k), 0);
      }
    }
  }
  FreeArrND(pot);
  FreeLibrary(&lib);
} //}}}

// ---- the shipped example library ------------------------------------------

// Examples/library must stay parseable and self-consistent //{{{
/*
 * Only invariants are asserted, not exact counts, so the shipped library can
 * be extended without touching this test - but a malformed edit still fails.
 */
static void test_example_library_consistent(void) {
  LIBRARY lib = ReadLibrary(EXAMPLE_LIB);
  SYSTEM *S = &lib.System;
  CHECK(S->Count.BeadType > 0);
  CHECK(S->Count.BondType > 0);
  CHECK(S->Count.AngleType > 0);
  // every bond/angle id maps to a type that exists
  CHECK(lib.n_bond_ids == S->Count.BondType);
  for (int i = 0; i < lib.n_bond_ids; i++) {
    CHECK(lib.bond_id[i].index >= 0);
    CHECK(lib.bond_id[i].index < S->Count.BondType);
    CHECK(lib.bond_id[i].id[0] != '\0');
    // a bond length of zero would mean the 'd' exponent failed to parse
    CHECK(S->BondType[lib.bond_id[i].index].b > 0);
  }
  CHECK(lib.n_angle_ids == S->Count.AngleType);
  for (int i = 0; i < lib.n_angle_ids; i++) {
    CHECK(lib.angle_id[i].index >= 0);
    CHECK(lib.angle_id[i].index < S->Count.AngleType);
  }
  // there is one self-interaction per bead type, so at least that many entries
  CHECK(lib.n_inter >= S->Count.BeadType);
  for (int i = 0; i < lib.n_inter; i++) {
    CHECK(lib.inter[i].A > 0);
    CHECK(lib.inter[i].Rc > 0);
  }
  // bead types are uniquely named, which the name-based lookups rely on
  for (int i = 0; i < S->Count.BeadType; i++) {
    for (int j = i + 1; j < S->Count.BeadType; j++) {
      CHECK(strcmp(S->BeadType[i].Name, S->BeadType[j].Name) != 0);
    }
  }
  FreeLibrary(&lib);
} //}}}
// the example molecules build, and their declared bead counts hold up //{{{
static void test_example_molecules_build(void) {
  const char *mols[] = {"water", "FA_C16", "CTAC", "Cl"};
  for (size_t t = 0; t < sizeof mols / sizeof *mols; t++) {
    LIB_MOL_INFO info = LibraryMoleculeInfo(EXAMPLE_LIB, mols[t]);
    CHECK(info.n_beads > 0);
    if (info.n_beads <= 0) {
      continue;
    }
    LIBRARY lib = ReadLibrary(EXAMPLE_LIB);
    SYSTEM *S = &lib.System;
    ReadLibraryMolecule(EXAMPLE_LIB, mols[t], 2, true, &lib);
    // whatever the counterions are, the totals must add up
    CHECK(S->Count.Bead == 2 * info.n_beads_total);
    if (info.n_beads > 1) {
      CHECK(S->MoleculeType[0].nBeads == info.n_beads);
    }
    // whatever was built, every bead has a valid type
    for (int i = 0; i < S->Count.Bead; i++) {
      CHECK(S->Bead[i].Type >= 0 && S->Bead[i].Type < S->Count.BeadType);
    }
    // and every bond points inside its molecule and at a real bond type
    for (int i = 0; i < S->Count.MoleculeType; i++) {
      MOLECULETYPE *mt = &S->MoleculeType[i];
      for (int j = 0; j < mt->nBonds; j++) {
        CHECK(mt->Bond[j][0] >= 0 && mt->Bond[j][0] < mt->nBeads);
        CHECK(mt->Bond[j][1] >= 0 && mt->Bond[j][1] < mt->nBeads);
        CHECK(mt->Bond[j][2] >= 0 && mt->Bond[j][2] < S->Count.BondType);
      }
      for (int j = 0; j < mt->nAngles; j++) {
        for (int k = 0; k < 3; k++) {
          CHECK(mt->Angle[j][k] >= 0 && mt->Angle[j][k] < mt->nBeads);
        }
        CHECK(mt->Angle[j][3] >= 0 && mt->Angle[j][3] < S->Count.AngleType);
      }
    }
    finish_and_free(&lib);
  }
} //}}}

int main(void) {
  RUN(test_read_library);
  RUN(test_molecule_info);
  RUN(test_trailing_slash);
  RUN(test_single_bead_molecule);
  RUN(test_chain_molecule);
  RUN(test_multiple_molecule_types);
  RUN(test_molecule_with_counterion);
  RUN(test_molecule_without_counterion);
  RUN(test_counterion_count);
  RUN(test_fill_pot);
  RUN(test_example_library_consistent);
  RUN(test_example_molecules_build);
  return test_main_end();
}
