/*
 * Coordinate-manipulation tests.
 *
 * t_geometry.c covers the pure box maths; this covers what happens to actual
 * molecules: joining them across periodic boundaries, wrapping them back into
 * the cell, and the centre/shape descriptors built on top.
 *
 * RemovePBCMolecule is the interesting one. Its bond-walking loop has to cope
 * with bond graphs that are not simple chains - rings, molecules whose bonds
 * fall into two disconnected pieces, and beads with no bonds at all - and
 * none of those shapes appear in the committed fixtures. The invariant used
 * throughout is that a correctly joined molecule has, for every bond, a raw
 * coordinate difference equal to the minimum-image difference: that is what
 * "no longer wrapped" means.
 *
 * Systems here are built and freed by hand rather than via InitSystem /
 * FreeSystem, so a test owns exactly what it allocates.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>

static pcg32_random_t rng;
static double randr(double lo, double hi) {
  return lo + (hi - lo) * (pcg32Rand0Int(&rng, 1u << 24) / (double)(1u << 24));
}

// ---- system builders ------------------------------------------------------

static BOX ortho_box(double x, double y, double z) { //{{{
  BOX b = InitBox;
  b.Length = (vec3d){ .v = {x, y, z} };
  b.alpha = 90;
  b.beta = 90;
  b.gamma = 90;
  CHECK(CalculateBoxData(&b, 0));
  return b;
} //}}}
// a system of n loose beads, one bead type each so masses can differ //{{{
static SYSTEM make_beads(BOX box, int n, const double *mass,
                         const vec3d *pos) {
  SYSTEM s = (SYSTEM){0};
  s.Box = box;
  s.Count = InitCount;
  s.Count.Bead = n;
  s.Count.BeadCoor = n;
  s.Count.BeadType = n;
  s.Count.Unbonded = n;
  s.BeadType = malloc(n * sizeof *s.BeadType);
  s.Bead = malloc(n * sizeof *s.Bead);
  s.BeadCoor = malloc(n * sizeof *s.BeadCoor);
  CHECK(s.BeadType && s.Bead && s.BeadCoor);
  for (int i = 0; i < n; i++) {
    InitBeadType(&s.BeadType[i]);
    snprintf(s.BeadType[i].Name, BEAD_NAME, "B%d", i);
    s.BeadType[i].Mass = mass[i];
    s.BeadType[i].Charge = 0;
    s.BeadType[i].Number = 0; // no Index array is allocated, so keep this 0
    InitBead(&s.Bead[i]);
    s.Bead[i].Type = i;
    s.Bead[i].Molecule = -1;
    s.Bead[i].Position = pos[i];
    s.Bead[i].InTimestep = true;
    s.BeadCoor[i] = i;
  }
  return s;
}
static void free_beads(SYSTEM *s) {
  free(s->BeadType);
  free(s->Bead);
  free(s->BeadCoor);
} //}}}
// a system holding one molecule with the given bond list //{{{
/*
 * All beads share a single bead type of unit mass; only the topology and the
 * positions matter for the joining tests.
 */
static SYSTEM make_molecule(BOX box, int n_beads, int n_bonds,
                            const int (*bonds)[2], const vec3d *pos) {
  SYSTEM s = (SYSTEM){0};
  s.Box = box;
  s.Count = InitCount;
  s.Count.Bead = n_beads;
  s.Count.BeadCoor = n_beads;
  s.Count.BeadType = 1;
  s.Count.Bonded = n_beads;
  s.Count.Molecule = 1;
  s.Count.MoleculeType = 1;
  s.Count.HighestResid = 0;

  s.BeadType = malloc(sizeof *s.BeadType);
  InitBeadType(&s.BeadType[0]);
  s_strcpy(s.BeadType[0].Name, "A", BEAD_NAME);
  s.BeadType[0].Mass = 1;
  s.BeadType[0].Charge = 0;
  s.BeadType[0].Number = 0; // see make_beads

  s.MoleculeType = malloc(sizeof *s.MoleculeType);
  MOLECULETYPE *mt = &s.MoleculeType[0];
  InitMoleculeType(mt);
  s_strcpy(mt->Name, "mol", MOL_NAME);
  mt->Number = 1;
  mt->nBeads = n_beads;
  mt->Bead = malloc(n_beads * sizeof *mt->Bead);
  for (int i = 0; i < n_beads; i++) {
    mt->Bead[i] = 0;
  }
  mt->nBonds = n_bonds;
  if (n_bonds > 0) {
    mt->Bond = malloc(n_bonds * sizeof *mt->Bond);
    for (int i = 0; i < n_bonds; i++) {
      mt->Bond[i][0] = bonds[i][0];
      mt->Bond[i][1] = bonds[i][1];
      mt->Bond[i][2] = -1; // no bond type
    }
  }

  s.Molecule = malloc(sizeof *s.Molecule);
  InitMolecule(&s.Molecule[0]);
  s.Molecule[0].Type = 0;
  s.Molecule[0].Index = 0;
  s.Molecule[0].InTimestep = true;
  s.Molecule[0].Bead = malloc(n_beads * sizeof *s.Molecule[0].Bead);

  s.Bead = malloc(n_beads * sizeof *s.Bead);
  s.BeadCoor = malloc(n_beads * sizeof *s.BeadCoor);
  for (int i = 0; i < n_beads; i++) {
    InitBead(&s.Bead[i]);
    s.Bead[i].Type = 0;
    s.Bead[i].Molecule = 0;
    s.Bead[i].Position = pos[i];
    s.Bead[i].InTimestep = true;
    s.BeadCoor[i] = i;
    s.Molecule[0].Bead[i] = i;
  }
  return s;
}
static void free_molecule(SYSTEM *s) {
  free(s->MoleculeType[0].Bead);
  if (s->MoleculeType[0].nBonds > 0) {
    free(s->MoleculeType[0].Bond);
  }
  free(s->MoleculeType);
  free(s->Molecule[0].Bead);
  free(s->Molecule);
  free(s->BeadType);
  free(s->Bead);
  free(s->BeadCoor);
} //}}}

// every bond must span its true minimum-image distance once joined //{{{
/*
 * If a bond still crosses a periodic boundary, its raw coordinate difference
 * differs from the minimum-image one by a whole box length in some axis.
 */
static void check_bonds_joined(const char *label, const SYSTEM s) {
  const MOLECULETYPE *mt = &s.MoleculeType[0];
  for (int i = 0; i < mt->nBonds; i++) {
    int a = s.Molecule[0].Bead[mt->Bond[i][0]];
    int b = s.Molecule[0].Bead[mt->Bond[i][1]];
    vec3d raw;
    for (int dd = 0; dd < 3; dd++) {
      raw.v[dd] = s.Bead[a].Position.v[dd] - s.Bead[b].Position.v[dd];
    }
    vec3d mi = Distance(s.Bead[a].Position, s.Bead[b].Position,
                        s.Box.OrthoLength);
    for (int dd = 0; dd < 3; dd++) {
      if (fabs(raw.v[dd] - mi.v[dd]) > 1e-9) {
        g_failures++;
        g_checks++;
        fprintf(stderr, "  FAIL %s: bond %d (beads %d-%d) still wrapped on "
                "axis %d: raw %.6g vs min-image %.6g\n",
                label, i, a, b, dd, raw.v[dd], mi.v[dd]);
        return;
      }
      g_checks++;
    }
  }
} //}}}

// ---- RemovePBCMolecule ----------------------------------------------------

// a linear chain laid across a boundary is pulled back together //{{{
static void test_join_linear_chain(void) {
  BOX box = ortho_box(10, 10, 10);
  // four beads 1 apart along x, starting near the upper edge and wrapping
  vec3d pos[4] = {
    { .v = {9.0, 5, 5} }, { .v = {0.0, 5, 5} },
    { .v = {1.0, 5, 5} }, { .v = {2.0, 5, 5} },
  };
  const int bonds[3][2] = {{0, 1}, {1, 2}, {2, 3}};
  SYSTEM s = make_molecule(box, 4, 3, bonds, pos);
  RemovePBCMolecule(0, &s);
  check_bonds_joined("chain", s);
  // consecutive beads end up exactly 1 apart, as they were placed
  for (int i = 0; i < 3; i++) {
    double d = fabs(s.Bead[i + 1].Position.x - s.Bead[i].Position.x);
    CHECK_CLOSE(d, 1.0, 1e-9);
  }
  free_molecule(&s);
} //}}}
// a ring: the bond graph has a cycle, so a bond is revisited //{{{
/*
 * The connected-bond walk adds a bond as soon as it shares a bead with one
 * already connected, so a cycle makes the last bond arrive when both of its
 * beads have already been placed. That branch does nothing, and the ring must
 * still come out joined.
 */
static void test_join_ring(void) {
  BOX box = ortho_box(10, 10, 10);
  // unit square in the xy plane, straddling the x boundary
  vec3d pos[4] = {
    { .v = {9.5, 5.0, 5} }, { .v = {0.5, 5.0, 5} },
    { .v = {0.5, 6.0, 5} }, { .v = {9.5, 6.0, 5} },
  };
  const int bonds[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  SYSTEM s = make_molecule(box, 4, 4, bonds, pos);
  RemovePBCMolecule(0, &s);
  check_bonds_joined("ring", s);
  // the square keeps its shape: opposite corners stay 1 apart in each axis
  CHECK_CLOSE(fabs(s.Bead[1].Position.x - s.Bead[0].Position.x), 1.0, 1e-9);
  CHECK_CLOSE(fabs(s.Bead[2].Position.y - s.Bead[1].Position.y), 1.0, 1e-9);
  free_molecule(&s);
} //}}}
// two disconnected bond components inside one molecule //{{{
/*
 * The outer while-loop runs once per component, and each new component is
 * anchored near the centre of what has already been placed. Both halves must
 * come out internally joined.
 */
static void test_join_two_components(void) {
  BOX box = ortho_box(10, 10, 10);
  vec3d pos[6] = {
    // component A, wrapped across x
    { .v = {9.5, 2, 2} }, { .v = {0.5, 2, 2} }, { .v = {1.5, 2, 2} },
    // component B, wrapped across y
    { .v = {5, 9.5, 8} }, { .v = {5, 0.5, 8} }, { .v = {5, 1.5, 8} },
  };
  const int bonds[4][2] = {{0, 1}, {1, 2}, {3, 4}, {4, 5}};
  SYSTEM s = make_molecule(box, 6, 4, bonds, pos);
  RemovePBCMolecule(0, &s);
  check_bonds_joined("two components", s);
  // each component keeps its internal spacing
  CHECK_CLOSE(fabs(s.Bead[1].Position.x - s.Bead[0].Position.x), 1.0, 1e-9);
  CHECK_CLOSE(fabs(s.Bead[2].Position.x - s.Bead[1].Position.x), 1.0, 1e-9);
  CHECK_CLOSE(fabs(s.Bead[4].Position.y - s.Bead[3].Position.y), 1.0, 1e-9);
  CHECK_CLOSE(fabs(s.Bead[5].Position.y - s.Bead[4].Position.y), 1.0, 1e-9);
  free_molecule(&s);
} //}}}
// beads with no bonds are placed next to the bonded cluster //{{{
/*
 * The final block of RemovePBCMolecule handles beads the bond walk never
 * touched. They must end up at their minimum-image position relative to the
 * cluster, not left on the far side of the box.
 */
static void test_join_unbonded_beads(void) {
  BOX box = ortho_box(10, 10, 10);
  vec3d pos[5] = {
    { .v = {9.0, 5, 5} }, { .v = {0.0, 5, 5} }, { .v = {1.0, 5, 5} },
    // two beads in the molecule that no bond mentions, far across the box
    { .v = {9.6, 5, 5} }, { .v = {0.4, 5, 5} },
  };
  const int bonds[2][2] = {{0, 1}, {1, 2}};
  SYSTEM s = make_molecule(box, 5, 2, bonds, pos);
  RemovePBCMolecule(0, &s);
  check_bonds_joined("unbonded", s);
  // the bonded cluster's centre, and the unbonded beads measured against it
  vec3d ref = {0};
  for (int i = 0; i < 3; i++) {
    for (int dd = 0; dd < 3; dd++) {
      ref.v[dd] += s.Bead[i].Position.v[dd] / 3;
    }
  }
  for (int i = 3; i < 5; i++) {
    vec3d raw;
    for (int dd = 0; dd < 3; dd++) {
      raw.v[dd] = s.Bead[i].Position.v[dd] - ref.v[dd];
    }
    vec3d mi = Distance(s.Bead[i].Position, ref, box.OrthoLength);
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(raw.v[dd], mi.v[dd], 1e-9);
    }
  }
  free_molecule(&s);
} //}}}
// molecules with no bonds, or absent from the timestep, are left alone //{{{
static void test_join_early_returns(void) {
  BOX box = ortho_box(10, 10, 10);
  vec3d pos[2] = { { .v = {9.5, 1, 1} }, { .v = {0.5, 1, 1} } };
  const int bonds[1][2] = {{0, 1}};

  // no bonds at all: nothing to join, positions untouched
  SYSTEM a = make_molecule(box, 2, 0, nullptr, pos);
  RemovePBCMolecule(0, &a);
  CHECK_CLOSE(a.Bead[0].Position.x, 9.5, 1e-12);
  CHECK_CLOSE(a.Bead[1].Position.x, 0.5, 1e-12);
  free_molecule(&a);

  // molecule not in the timestep: also untouched
  SYSTEM b = make_molecule(box, 2, 1, bonds, pos);
  b.Molecule[0].InTimestep = false;
  RemovePBCMolecule(0, &b);
  CHECK_CLOSE(b.Bead[0].Position.x, 9.5, 1e-12);
  CHECK_CLOSE(b.Bead[1].Position.x, 0.5, 1e-12);
  free_molecule(&b);
} //}}}
// joining is idempotent: an already-joined molecule does not drift //{{{
static void test_join_idempotent(void) {
  BOX box = ortho_box(10, 10, 10);
  vec3d pos[4] = {
    { .v = {9.0, 5, 5} }, { .v = {0.0, 5, 5} },
    { .v = {1.0, 5, 5} }, { .v = {2.0, 5, 5} },
  };
  const int bonds[3][2] = {{0, 1}, {1, 2}, {2, 3}};
  SYSTEM s = make_molecule(box, 4, 3, bonds, pos);
  RemovePBCMolecule(0, &s);
  vec3d after[4];
  for (int i = 0; i < 4; i++) {
    after[i] = s.Bead[i].Position;
  }
  RemovePBCMolecule(0, &s);
  for (int i = 0; i < 4; i++) {
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(s.Bead[i].Position.v[dd], after[i].v[dd], 1e-9);
    }
  }
  free_molecule(&s);
} //}}}

// ---- WrapJoinCoordinates --------------------------------------------------

// wrapping puts every coordinate into [0, OrthoLength) //{{{
static void test_wrap_coordinates(void) {
  BOX box = ortho_box(10, 12, 8);
  const int n = 200;
  double *mass = malloc(n * sizeof *mass);
  vec3d *pos = malloc(n * sizeof *pos);
  for (int i = 0; i < n; i++) {
    mass[i] = 1;
    pos[i] = (vec3d){ .v = { randr(-30, 30), randr(-30, 30), randr(-30, 30) } };
  }
  SYSTEM s = make_beads(box, n, mass, pos);
  WrapJoinCoordinates(&s, true, false);
  for (int i = 0; i < n; i++) {
    for (int dd = 0; dd < 3; dd++) {
      CHECK(s.Bead[i].Position.v[dd] >= 0);
      CHECK(s.Bead[i].Position.v[dd] < box.OrthoLength.v[dd]);
    }
  }
  // wrapping again changes nothing
  vec3d first = s.Bead[0].Position;
  WrapJoinCoordinates(&s, true, false);
  for (int dd = 0; dd < 3; dd++) {
    CHECK_CLOSE(s.Bead[0].Position.v[dd], first.v[dd], 1e-12);
  }
  free_beads(&s);
  free(mass);
  free(pos);
} //}}}
// with no box, and with both flags false, nothing happens //{{{
static void test_wrap_no_box_or_no_flags(void) {
  BOX box = ortho_box(10, 10, 10);
  double mass[2] = {1, 1};
  vec3d pos[2] = { { .v = {15, -3, 25} }, { .v = {-1, 2, 3} } };

  // Volume == -1 marks "no box information", so coordinates must be left be
  SYSTEM a = make_beads(box, 2, mass, pos);
  a.Box.Volume = -1;
  WrapJoinCoordinates(&a, true, true);
  CHECK_CLOSE(a.Bead[0].Position.x, 15, 1e-12);
  CHECK_CLOSE(a.Bead[0].Position.z, 25, 1e-12);
  free_beads(&a);

  // neither wrapping nor joining requested
  SYSTEM b = make_beads(box, 2, mass, pos);
  WrapJoinCoordinates(&b, false, false);
  CHECK_CLOSE(b.Bead[0].Position.x, 15, 1e-12);
  CHECK_CLOSE(b.Bead[0].Position.y, -3, 1e-12);
  free_beads(&b);
} //}}}

// ---- centres --------------------------------------------------------------

// CentreOfMass: weighting, and translation equivariance //{{{
static void test_centre_of_mass(void) {
  BOX box = ortho_box(100, 100, 100);
  double mass[3] = {1, 1, 2};
  vec3d pos[3] = {
    { .v = {0, 0, 0} }, { .v = {2, 0, 0} }, { .v = {4, 0, 0} },
  };
  int list[3] = {0, 1, 2};
  SYSTEM s = make_beads(box, 3, mass, pos);
  // (1*0 + 1*2 + 2*4) / 4 = 2.5
  vec3d com = CentreOfMass(3, list, s);
  CHECK_CLOSE(com.x, 2.5, 1e-12);
  CHECK_CLOSE(com.y, 0, 1e-12);

  // shifting every bead shifts the centre by the same amount
  for (int i = 0; i < 3; i++) {
    s.Bead[i].Position.x += 7;
    s.Bead[i].Position.y -= 3;
  }
  vec3d moved = CentreOfMass(3, list, s);
  CHECK_CLOSE(moved.x, 2.5 + 7, 1e-12);
  CHECK_CLOSE(moved.y, -3, 1e-12);

  // a single bead sits at its own position
  int one[1] = {1};
  vec3d single = CentreOfMass(1, one, s);
  CHECK_CLOSE(single.x, s.Bead[1].Position.x, 1e-12);
  free_beads(&s);
} //}}}
// an unspecified mass makes CentreOfMass bail out //{{{
/*
 * MASS is the "not set" sentinel. The function warns and returns early rather
 * than dividing by a wrong total, so the caller gets a value that is not a
 * centre of mass at all - worth pinning, since nothing in the signature says
 * the result can be meaningless.
 */
static void test_centre_of_mass_unspecified(void) {
  BOX box = ortho_box(100, 100, 100);
  double mass[3] = {MASS, 1, 1};
  vec3d pos[3] = {
    { .v = {1, 1, 1} }, { .v = {2, 2, 2} }, { .v = {3, 3, 3} },
  };
  int list[3] = {0, 1, 2};
  SYSTEM s = make_beads(box, 3, mass, pos);
  // the very first bead is massless, so nothing is accumulated at all
  vec3d com = CentreOfMass(3, list, s);
  CHECK_CLOSE(com.x, 0, 1e-12);
  CHECK_CLOSE(com.y, 0, 1e-12);
  CHECK_CLOSE(com.z, 0, 1e-12);
  free_beads(&s);
} //}}}
// GeomCentre averages the beads that are in the timestep //{{{
static void test_geom_centre(void) {
  BOX box = ortho_box(100, 100, 100);
  double mass[4] = {1, 1, 1, 1};
  vec3d pos[4] = {
    { .v = {0, 0, 0} }, { .v = {2, 0, 0} },
    { .v = {0, 4, 0} }, { .v = {50, 50, 50} },
  };
  int list[4] = {0, 1, 2, 3};
  SYSTEM s = make_beads(box, 4, mass, pos);
  // the outlier is not in the timestep and must be skipped entirely
  s.Bead[3].InTimestep = false;
  vec3d gc = GeomCentre(4, list, s.Bead);
  CHECK_CLOSE(gc.x, 2.0 / 3, 1e-12);
  CHECK_CLOSE(gc.y, 4.0 / 3, 1e-12);
  CHECK_CLOSE(gc.z, 0, 1e-12);

  // with every bead present it is the plain average
  s.Bead[3].InTimestep = true;
  vec3d all = GeomCentre(4, list, s.Bead);
  CHECK_CLOSE(all.x, 52.0 / 4, 1e-12);

  /*
   * No bead in the timestep divides by zero. Pinned as NaN rather than left
   * to chance: callers get a value that will poison anything it touches.
   */
  for (int i = 0; i < 4; i++) {
    s.Bead[i].InTimestep = false;
  }
  vec3d none = GeomCentre(4, list, s.Bead);
  CHECK(isnan(none.x));
  CHECK(isnan(none.y));
  CHECK(isnan(none.z));
  free_beads(&s);
} //}}}

// ---- Gyration -------------------------------------------------------------

// Gyration recentres the beads it is given - an undocumented side effect //{{{
/*
 * The header gives no hint that the positions are modified, but Gyration
 * subtracts the geometric centre from every listed bead before building the
 * tensor. Anything calling it mid-analysis sees its coordinates move.
 */
static void test_gyration_mutates_positions(void) {
  BOX box = ortho_box(100, 100, 100);
  double mass[3] = {1, 1, 1};
  vec3d pos[3] = {
    { .v = {10, 20, 30} }, { .v = {11, 20, 30} }, { .v = {12, 20, 30} },
  };
  int list[3] = {0, 1, 2};
  SYSTEM s = make_beads(box, 3, mass, pos);
  Gyration(3, list, &s);
  // the beads are now centred on the origin
  vec3d gc = GeomCentre(3, list, s.Bead);
  CHECK_CLOSE(gc.x, 0, 1e-9);
  CHECK_CLOSE(gc.y, 0, 1e-9);
  CHECK_CLOSE(gc.z, 0, 1e-9);
  // and the original coordinates are gone
  CHECK(fabs(s.Bead[0].Position.x - 10) > 1e-6);
  free_beads(&s);
} //}}}
// analytically known shapes //{{{
static void test_gyration_shapes(void) {
  BOX box = ortho_box(100, 100, 100);
  {
    // a rod along x: one non-zero eigenvalue, returned last (ascending)
    double mass[3] = {1, 1, 1};
    vec3d pos[3] = {
      { .v = {-1, 0, 0} }, { .v = {0, 0, 0} }, { .v = {1, 0, 0} },
    };
    int list[3] = {0, 1, 2};
    SYSTEM s = make_beads(box, 3, mass, pos);
    vec3d e = Gyration(3, list, &s);
    CHECK_CLOSE(e.x, 0, 1e-9);
    CHECK_CLOSE(e.y, 0, 1e-9);
    CHECK_CLOSE(e.z, 2.0 / 3, 1e-9); // <x^2> = (1 + 0 + 1)/3
    // eigenvalues come back in ascending order
    CHECK(e.x <= e.y);
    CHECK(e.y <= e.z);
    free_beads(&s);
  }
  {
    // the eight corners of a cube: fully isotropic, all eigenvalues equal
    double mass[8];
    vec3d pos[8];
    int list[8];
    int k = 0;
    for (int i = -1; i <= 1; i += 2) {
      for (int j = -1; j <= 1; j += 2) {
        for (int l = -1; l <= 1; l += 2) {
          mass[k] = 1;
          pos[k] = (vec3d){ .v = {i, j, l} };
          list[k] = k;
          k++;
        }
      }
    }
    SYSTEM s = make_beads(box, 8, mass, pos);
    vec3d e = Gyration(8, list, &s);
    CHECK_CLOSE(e.x, 1.0, 1e-9);
    CHECK_CLOSE(e.y, 1.0, 1e-9);
    CHECK_CLOSE(e.z, 1.0, 1e-9);
    free_beads(&s);
  }
} //}}}
// a shape that is not axis-aligned, where the off-diagonal terms matter //{{{
/*
 * Two beads on the xy diagonal give the tensor [[1,1,0],[1,1,0],[0,0,0]],
 * whose eigenvalues are 0, 0 and 2. Reading the diagonal instead - which is
 * what happens if the tensor's lower triangle is left unfilled, since
 * gsl_eigen_symmv takes its input from there - would give 0, 1, 1 and a rod
 * would be mistaken for a disc. Every axis-aligned test case passes either
 * way, so this is the one that distinguishes them.
 */
static void test_gyration_off_diagonal(void) {
  BOX box = ortho_box(100, 100, 100);
  double mass[2] = {1, 1};
  vec3d pos[2] = { { .v = {1, 1, 0} }, { .v = {-1, -1, 0} } };
  int list[2] = {0, 1};
  SYSTEM s = make_beads(box, 2, mass, pos);
  vec3d e = Gyration(2, list, &s);
  CHECK_CLOSE(e.x, 0.0, 1e-9);
  CHECK_CLOSE(e.y, 0.0, 1e-9);
  CHECK_CLOSE(e.z, 2.0, 1e-9);
  free_beads(&s);

  // the same rod rotated onto the x axis must give the same eigenvalues
  double r = sqrt(2.0); // same length as the diagonal pair above
  vec3d aligned[2] = { { .v = {r, 0, 0} }, { .v = {-r, 0, 0} } };
  SYSTEM t = make_beads(box, 2, mass, aligned);
  vec3d ea = Gyration(2, list, &t);
  CHECK_CLOSE(ea.x, 0.0, 1e-9);
  CHECK_CLOSE(ea.y, 0.0, 1e-9);
  CHECK_CLOSE(ea.z, 2.0, 1e-9);
  free_beads(&t);
} //}}}
// the eigenvalues do not depend on how the shape is oriented //{{{
static void test_gyration_rotation_invariant(void) {
  BOX box = ortho_box(1000, 1000, 1000);
  const int n = 12;
  double mass[12];
  vec3d pos[12];
  int list[12];
  for (int i = 0; i < n; i++) {
    mass[i] = 1;
    pos[i] = (vec3d){ .v = { randr(-5, 5), randr(-5, 5), randr(-5, 5) } };
    list[i] = i;
  }
  SYSTEM a = make_beads(box, n, mass, pos);
  vec3d ref = Gyration(n, list, &a);
  free_beads(&a);

  // rotate by 40 degrees about z, then 25 about x, and translate far away
  double t1 = 40 * PI / 180, t2 = 25 * PI / 180;
  vec3d rot[12];
  for (int i = 0; i < n; i++) {
    double x = pos[i].x * cos(t1) - pos[i].y * sin(t1);
    double y = pos[i].x * sin(t1) + pos[i].y * cos(t1);
    double z = pos[i].z;
    double y2 = y * cos(t2) - z * sin(t2);
    double z2 = y * sin(t2) + z * cos(t2);
    rot[i] = (vec3d){ .v = {x + 100, y2 - 50, z2 + 200} };
  }
  SYSTEM b = make_beads(box, n, mass, rot);
  vec3d e = Gyration(n, list, &b);
  free_beads(&b);

  CHECK_CLOSE(e.x, ref.x, 1e-7);
  CHECK_CLOSE(e.y, ref.y, 1e-7);
  CHECK_CLOSE(e.z, ref.z, 1e-7);
} //}}}

int main(void) {
  pcg32Seed(&rng, 20260808u);
  RUN(test_join_linear_chain);
  RUN(test_join_ring);
  RUN(test_join_two_components);
  RUN(test_join_unbonded_beads);
  RUN(test_join_early_returns);
  RUN(test_join_idempotent);
  RUN(test_wrap_coordinates);
  RUN(test_wrap_no_box_or_no_flags);
  RUN(test_centre_of_mass);
  RUN(test_centre_of_mass_unspecified);
  RUN(test_geom_centre);
  RUN(test_gyration_mutates_positions);
  RUN(test_gyration_shapes);
  RUN(test_gyration_off_diagonal);
  RUN(test_gyration_rotation_invariant);
  return test_main_end();
}
