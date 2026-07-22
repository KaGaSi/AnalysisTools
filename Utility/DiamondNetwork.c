#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Generate a 3D tetrafunctional polymer network based on the diamond cubic "
  "lattice. Each junction ('J') has exactly 4 nearest neighbours arranged "
  "tetrahedrally (functionality 4). The conventional unit cell contains 8 "
  "junctions: 4 on FCC sites and 4 on tetrahedral interstitial sites. "
  "Adjacent junctions are connected by strands of 'n' intermediate beads "
  "('S'). The box size is set with -b; the number of unit cells in each "
  "direction is determined by rounding box/A_ideal to the nearest integer, "
  "where A_ideal = 4*(n+1)*l/sqrt(3). The lattice is then squished uniformly "
  "per axis to fill the box exactly. "
  "Output: one single-frame coordinate file; all beads form one molecule.",

  "Usage: DiamondNetwork <output> [options]",
  .args = 1,
  .all = 8,
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<output>", nullptr, "output coordinate file", OPT_ARG},
  {"-n", "<int>", "beads per strand between adjacent junctions (default: 3)",
    OPT_EXTRA},
  {"-b", "3*<float>", "box size Lx Ly Lz; unit cells = round(L/A_ideal), "
    "then squished to fit (default: 2*A_ideal in each direction)", OPT_EXTRA},
  {"-l", "<float>", "bond length l; ideal lattice param A = 4*(n+1)*l/sqrt(3) "
    "(default: 1.0)", OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int n_strand;   // -n
  vec3d box_size; // -b
  double bond_l;  // -l
}; //}}}

/*
 * Diamond cubic lattice geometry
 * ================================
 * Sublattice fractional positions within one conventional unit cell
 *   0: A0 = (0,   0,   0  )   FCC corner
 *   1: A1 = (1/2, 1/2, 0  )   FCC face xy
 *   2: A2 = (1/2, 0,   1/2)   FCC face xz
 *   3: A3 = (0,   1/2, 1/2)   FCC face yz
 *   4: B0 = (1/4, 1/4, 1/4)   tetrahedral site
 *   5: B1 = (3/4, 3/4, 1/4)
 *   6: B2 = (3/4, 1/4, 3/4)
 *   7: B3 = (1/4, 3/4, 3/4)
 *
 * Every bond goes from an A-type atom to a B-type atom. The four bond
 * directions are the same for every A atom:
 *   b=0: (+1/4, +1/4, +1/4)
 *   b=1: (+1/4, -1/4, -1/4)
 *   b=2: (-1/4, +1/4, -1/4)
 *   b=3: (-1/4, -1/4, +1/4)
 *
 * nbr[tA][b] = {B_sublattice, dcx, dcy, dcz} gives, for A-type atom tA and bond
 * direction b, the sublattice index of the target B atom and the cell offset
 * (with PBC).
 */

// sublattice fractional positions (units of A) //{{{
static const vec3d frac_pos[8] = {
  {.v = {0.00, 0.00, 0.00}}, // 0: A0
  {.v = {0.50, 0.50, 0.00}}, // 1: A1
  {.v = {0.50, 0.00, 0.50}}, // 2: A2
  {.v = {0.00, 0.50, 0.50}}, // 3: A3
  {.v = {0.25, 0.25, 0.25}}, // 4: B0
  {.v = {0.75, 0.75, 0.25}}, // 5: B1
  {.v = {0.75, 0.25, 0.75}}, // 6: B2
  {.v = {0.25, 0.75, 0.75}}, // 7: B3
}; //}}}

// A→B bond connectivity table //{{{
static const int nbr[4][4][4] = {
  {{4, 0, 0, 0}, {7, 0,-1,-1}, {6,-1, 0,-1}, {5,-1,-1, 0}}, // A0
  {{5, 0, 0, 0}, {6, 0, 0,-1}, {7, 0, 0,-1}, {4, 0, 0, 0}}, // A1
  {{6, 0, 0, 0}, {5, 0,-1, 0}, {4, 0, 0, 0}, {7, 0,-1, 0}}, // A2
  {{7, 0, 0, 0}, {4, 0, 0, 0}, {5,-1, 0, 0}, {6,-1, 0, 0}}, // A3
}; //}}}

// bond vectors from any A atom (units of A), indexed by b //{{{
static const vec3d bond_vec[4] = {
  {.v = {+0.25, +0.25, +0.25}}, // b=0
  {.v = {+0.25, -0.25, -0.25}}, // b=1
  {.v = {-0.25, +0.25, -0.25}}, // b=2
  {.v = {-0.25, -0.25, +0.25}}, // b=3
}; //}}}

int main(int argc, char *argv[]) {

  // command-line arguments //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  FILE_TYPE fw_coor;
  snprintf(fw_coor.name, LINE, "%s", argv[++count]);
  fw_coor.type = CoordinateFileType(fw_coor.name);

  SYS_FILES in = InitSysFiles;
  COMMON_OPT commons = CommonOptions(argc, argv, in);

  opt.n_strand = 3;
  OneNumberOption(argc, argv, "-n", &opt.n_strand, 'i');
  opt.bond_l = 1.0;
  OneNumberOption(argc, argv, "-l", &opt.bond_l, 'd');
  // default box = 2 ideal unit cells in each direction
  double A_default = 4.0 * (opt.n_strand + 1) * opt.bond_l / sqrt(3.0);
  for (int dd = 0; dd < 3; dd++) {
    opt.box_size.v[dd] = 2.0 * A_default;
  }
  ThreeNumbersOption(argc, argv, "-b", opt.box_size.v, 'd'); //}}}

  // validate //{{{
  if (opt.n_strand < 1) {
    snprintf(ERROR_MSG, LINE, "-n must be >= 1 (got %d)", opt.n_strand);
    PrintErrorOption("-n");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  if (opt.bond_l <= 0) {
    snprintf(ERROR_MSG, LINE, "-l must be > 0 (got %g)", opt.bond_l);
    PrintErrorOption("-l");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  for (int dd = 0; dd < 3; dd++) {
    if (opt.box_size.v[dd] <= 0) {
      snprintf(ERROR_MSG, LINE, "-b dimensions must be > 0 (got %g)",
               opt.box_size.v[dd]);
      PrintErrorOption("-b");
      Help(true, HelpDesc, opts);
      exit(1);
    }
  } //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // derived sizes //{{{
  const int n = opt.n_strand;
  const double l = opt.bond_l;

  // ideal lattice parameter for the given n and l
  const double A_ideal = 4.0 * (n + 1) * l / sqrt(3.0);

  // number of unit cells: how many ideal cells fit in the box (ceiling), min 1
  vec3d L;
  for (int dd = 0; dd < 3; dd++) {
    L.v[dd] = opt.box_size.v[dd];
  }
  vec3i N;
  for (int dd = 0; dd < 3; dd++) {
    N.v[dd] = (int)ceil(L.v[dd] / A_ideal);
    if (N.v[dd] < 1) N.v[dd] = 1;
  }

  // per-axis lattice parameters: squish to fill the box exactly
  vec3d A;
  for (int dd = 0; dd < 3; dd++) {
    A.v[dd] = L.v[dd] / N.v[dd];
  }

  const int Ncells = N.v[0] * N.v[1] * N.v[2];
  const int N_junc = 8 * Ncells;
  // bonds between junctions = 4 A-types * 4 bonds * Ncells
  // = 16*Ncells = 2*N_junc (each bond counted once, A→B)
  const int N_strands = 16 * Ncells;
  const int N_strand_beads = N_strands * n;
  const int N_total = N_junc + N_strand_beads;
  const int N_bonds = N_strands * (n + 1); //}}}

  // build SYSTEM //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  Count->Bead = N_total;
  Count->BeadCoor = N_total;
  Count->Bonded = N_total;
  Count->BondedCoor = N_total;
  Count->Unbonded = 0;
  Count->Molecule = 1;
  Count->HighestResid = 0;
  Count->BondType = 1;

  for (int dd = 0; dd < 3; dd++) {
    System.Box.Length.v[dd] = L.v[dd];
  }

  System.Bead = realloc(System.Bead,    N_total * sizeof *System.Bead);
  System.BeadCoor = realloc(System.BeadCoor,N_total * sizeof *System.BeadCoor);

  // bead types
  NewBeadType(&System.BeadType, &Count->BeadType, "J", 0, 1.0, 0.5);
  System.BeadType[0].Number = N_junc;
  NewBeadType(&System.BeadType, &Count->BeadType, "S", 0, 1.0, 0.5);
  System.BeadType[1].Number = N_strand_beads;

  System.BondType[0].b = l;

  NewMolType(&System.MoleculeType, &Count->MoleculeType, "NET",
             N_total, N_bonds, 0, 0, 0);
  MOLECULETYPE *mt = &System.MoleculeType[0];
  mt->Number = 1; //}}}

  // bead-type sequence in molecule //{{{
  int bt_J = FindBeadType("J", System);
  int bt_S = FindBeadType("S", System);
  for (int i = 0; i < N_junc; i++) {
    mt->Bead[i] = bt_J;
  }
  for (int i = N_junc; i < N_total; i++) {
    mt->Bead[i] = bt_S;
  } //}}}

  // bond connectivity //{{{
  /*
   * Junction global index:
   *   j(cx,cy,cz,t) = (cz*Ny*Nx + cy*Nx + cx) * 8 + t
   *
   * Strand index (for strand from A-atom tA, bond b, in cell cell_idx):
   *   strand_idx = cell_idx * 16 + tA * 4 + b
   *
   * First strand bead: N_junc + strand_idx * n
   *
   * Topology: j_start -- s[0] -- ... -- s[n-1] -- j_end  → (n+1) bonds
   */
  int bi = 0;
  for (int iz = 0; iz < N.v[2]; iz++) {
    for (int iy = 0; iy < N.v[1]; iy++) {
      for (int ix = 0; ix < N.v[0]; ix++) {
        int cell_idx = iz*N.v[1]*N.v[0] + iy*N.v[0] + ix;
        vec3i cell = {.v = {ix, iy, iz}};
        for (int tA = 0; tA < 4; tA++) {
          int j_start = cell_idx * 8 + tA;
          for (int b = 0; b < 4; b++) {
            int B_sub = nbr[tA][b][0];
            vec3i nc;
            for (int dd = 0; dd < 3; dd++) {
              int tmp = cell.v[dd] + nbr[tA][b][dd+1];
              nc.v[dd] = (tmp % N.v[dd] + N.v[dd]) % N.v[dd];
            }
            int tmp = nc.v[2]*N.v[1]*N.v[0] + nc.v[1]*N.v[0] + nc.v[0];
            int j_end = tmp * 8 + B_sub;

            int s0 = N_junc + (cell_idx * 16 + tA * 4 + b) * n;

            mt->Bond[bi][0] = j_start;
            mt->Bond[bi][1] = s0;
            mt->Bond[bi][2] = 0; bi++;
            for (int k = 0; k < n-1; k++) {
              mt->Bond[bi][0] = s0 + k;
              mt->Bond[bi][1] = s0 + k + 1;
              mt->Bond[bi][2] = 0;
              bi++;
            }
            mt->Bond[bi][0] = s0 + n - 1;
            mt->Bond[bi][1] = j_end;
            mt->Bond[bi][2] = 0;
            bi++;
          }
        }
      }
    }
  } //}}}

  // fill Molecule and Bead structs //{{{
  MOLECULE *mol = &System.Molecule[0];
  InitMolecule(mol);
  mol->Type  = 0;
  mol->Index = 0;
  mol->Bead  = calloc(N_total, sizeof *mol->Bead);
  for (int i = 0; i < N_total; i++) {
    BEAD *b = &System.Bead[i];
    InitBead(b);
    b->Type = bt_S;
    if (i < N_junc) {
      b->Type = bt_J;
    }
    b->Molecule = 0;
    mol->Bead[i]      = i;
    System.BeadCoor[i] = i;
  } //}}}

  FinishSystem(&System);

  // assign straight-line coordinates //{{{
  /*
   * Junction (cx,cy,cz,t): position = ((cx + frac[t][0])*Ax,
   *                                    (cy + frac[t][1])*Ay,
   *                                    (cz + frac[t][2])*Az)
   * Strand bead k along bond b from junction j_start:
   *   pos = j_start + (k+1)/(n+1) * (bond_vec[b][0]*Ax,
   *                                   bond_vec[b][1]*Ay,
   *                                   bond_vec[b][2]*Az)
   * Ax/Ay/Az may differ (squished box), so bond lengths are uniform only
   * when Ax=Ay=Az.
   */
  for (int iz = 0; iz < N.v[2]; iz++) {
    for (int iy = 0; iy < N.v[1]; iy++) {
      for (int ix = 0; ix < N.v[0]; ix++) {
        int cell_idx = iz * N.v[1] * N.v[0] + iy * N.v[0] + ix;
        vec3i cell = {.v = {ix, iy, iz}};
        // junctions
        for (int t = 0; t < 8; t++) {
          BEAD *b = &System.Bead[cell_idx*8+t];
          for (int dd = 0; dd < 3; dd++) {
            b->Position.v[dd] = (cell.v[dd] + frac_pos[t].v[dd]) * A.v[dd];
          }
        }
        // strand beads
        for (int tA = 0; tA < 4; tA++) {
          vec3d origin;
          for (int dd = 0; dd < 3; dd++) {
            origin.v[dd] = (cell.v[dd] + frac_pos[tA].v[dd]) * A.v[dd];
          }
          for (int b = 0; b < 4; b++) {
            int s0 = N_junc + (cell_idx * 16 + tA * 4 + b) * n;
            vec3d bv;
            for (int dd = 0; dd < 3; dd++) {
              bv.v[dd] = bond_vec[b].v[dd] * A.v[dd];
            }
            for (int k = 0; k < n; k++) {
              double f = (double)(k + 1) / (n + 1);
              for (int dd = 0; dd < 3; dd++) {
                System.Bead[s0+k].Position.v[dd] = origin.v[dd] + f * bv.v[dd];
              }
            }
          }
        }
      }
    }
  } //}}}

  if (commons.verbose) {
    VerboseOutput(System);
  }

  WrapJoinCoordinates(&System, true, false);
  InitOutputCoorFile(fw_coor, System, argc, argv);
  WriteTimestepAll(fw_coor, System, 1, argc, argv);

  FreeSystem(&System);

  return 0;
}
