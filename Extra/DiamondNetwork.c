#include "../src/AnalysisTools.h"
#include <stdbool.h>

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Generate a 3D tetrafunctional polymer network based on the diamond cubic "
  "lattice. Each junction ('J') has exactly 4 nearest neighbours arranged "
  "tetrahedrally (functionality 4). The conventional unit cell contains 8 "
  "junctions: 4 on FCC sites and 4 on tetrahedral interstitial sites. "
  "Adjacent junctions are connected by strands of 'n' intermediate beads "
  "('S'). The lattice parameter is A = 4*(n+1)*l/sqrt(3), so every strand is "
  "a straight chain of (n+1) bonds of length l. The full system contains "
  "Nx x Ny x Nz conventional unit cells with periodic boundaries. "
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
  {"<output>", NULL, "output coordinate file (.vtf, .vcf, .xyz, .data, etc.)", OPT_ARG},
  {"-n", "<int>", "beads per strand between adjacent junctions, >=1 (default: 3)", OPT_EXTRA},
  {"-box", "<int> <int> <int>", "unit cells in x, y, z (default: 2 2 2)", OPT_EXTRA},
  {"-l", "<float>", "bond length between consecutive beads (default: 1.0)", OPT_EXTRA},
  {NULL}
}; //}}}

// structure for options //{{{
struct OPT {
  int n_strand; // -n
  int box[3];   // -box
  double bond_l;// -l
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
static const double frac_pos[8][3] = {
  {0.00, 0.00, 0.00}, // 0: A0
  {0.50, 0.50, 0.00}, // 1: A1
  {0.50, 0.00, 0.50}, // 2: A2
  {0.00, 0.50, 0.50}, // 3: A3
  {0.25, 0.25, 0.25}, // 4: B0
  {0.75, 0.75, 0.25}, // 5: B1
  {0.75, 0.25, 0.75}, // 6: B2
  {0.25, 0.75, 0.75}, // 7: B3
}; //}}}

// A→B bond connectivity table //{{{
static const int nbr[4][4][4] = {
  {{4, 0, 0, 0}, {7, 0,-1,-1}, {6,-1, 0,-1}, {5,-1,-1, 0}}, // A0
  {{5, 0, 0, 0}, {6, 0, 0,-1}, {7, 0, 0,-1}, {4, 0, 0, 0}}, // A1
  {{6, 0, 0, 0}, {5, 0,-1, 0}, {4, 0, 0, 0}, {7, 0,-1, 0}}, // A2
  {{7, 0, 0, 0}, {4, 0, 0, 0}, {5,-1, 0, 0}, {6,-1, 0, 0}}, // A3
}; //}}}

// bond vectors from any A atom (units of A), indexed by b //{{{
static const double bond_vec[4][3] = {
  {+0.25, +0.25, +0.25}, // b=0
  {+0.25, -0.25, -0.25}, // b=1
  {-0.25, +0.25, -0.25}, // b=2
  {-0.25, -0.25, +0.25}, // b=3
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
  opt.box[0] = opt.box[1] = opt.box[2] = 2;
  ThreeNumbersOption(argc, argv, "-box", opt.box, 'i');
  opt.bond_l = 1.0;
  OneNumberOption(argc, argv, "-l", &opt.bond_l, 'd'); //}}}

  // validate //{{{
  if (opt.n_strand < 1) {
    snprintf(ERROR_MSG, LINE, "-n must be >= 1 (got %d)", opt.n_strand);
    PrintErrorOption("-n");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  for (int dd = 0; dd < 3; dd++) {
    if (opt.box[dd] < 1) {
      snprintf(ERROR_MSG, LINE, "-box must be >= 1 in all dimensions (got %d)",
               opt.box[dd]);
      PrintErrorOption("-box");
      Help(true, HelpDesc, opts);
      exit(1);
    }
  }
  if (opt.bond_l <= 0) {
    snprintf(ERROR_MSG, LINE, "-l must be > 0 (got %g)", opt.bond_l);
    PrintErrorOption("-l");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}

  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  // derived sizes //{{{
  const int Nx = opt.box[0], Ny = opt.box[1], Nz = opt.box[2];
  const int n = opt.n_strand;
  const double l = opt.bond_l;

  // lattice parameter: |bond_vec|*A = A*sqrt(3)/4 = (n+1)*l
  const double A = 4.0 * (n + 1) * l / sqrt(3.0);

  const int Ncells        = Nx * Ny * Nz;
  const int N_junc        = 8 * Ncells;
  // bonds between junctions = 4 A-types * 4 bonds * Ncells
  // = 16*Ncells = 2*N_junc (each bond counted once, A→B)
  const int N_strands     = 16 * Ncells;
  const int N_strand_beads= N_strands * n;
  const int N_total       = N_junc + N_strand_beads;
  const int N_bonds       = N_strands * (n + 1);

  const double Lx = Nx * A;
  const double Ly = Ny * A;
  const double Lz = Nz * A; //}}}

  // build SYSTEM //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;

  Count->Bead       = N_total;
  Count->BeadCoor   = N_total;
  Count->Bonded     = N_total;
  Count->BondedCoor = N_total;
  Count->Unbonded   = 0;
  Count->Molecule   = 1;
  Count->HighestResid = 0;
  Count->BondType   = 1;

  System.Box.Length.v[0] = Lx;
  System.Box.Length.v[1] = Ly;
  System.Box.Length.v[2] = Lz;

  System.Bead    = realloc(System.Bead,    N_total * sizeof *System.Bead);
  System.BeadCoor= realloc(System.BeadCoor,N_total * sizeof *System.BeadCoor);

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
  for (int iz = 0; iz < Nz; iz++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int ix = 0; ix < Nx; ix++) {
        int cell_idx = iz*Ny*Nx + iy*Nx + ix;
        for (int tA = 0; tA < 4; tA++) {
          int j_start = cell_idx * 8 + tA;
          for (int b = 0; b < 4; b++) {
            int B_sub  = nbr[tA][b][0];
            int nx_c = ((ix + nbr[tA][b][1]) % Nx + Nx) % Nx;
            int ny_c = ((iy + nbr[tA][b][2]) % Ny + Ny) % Ny;
            int nz_c = ((iz + nbr[tA][b][3]) % Nz + Nz) % Nz;
            int j_end = (nz_c*Ny*Nx + ny_c*Nx + nx_c) * 8 + B_sub;

            int s0 = N_junc + (cell_idx * 16 + tA * 4 + b) * n;

            mt->Bond[bi][0] = j_start; mt->Bond[bi][1] = s0;     mt->Bond[bi][2] = 0; bi++;
            for (int k = 0; k < n-1; k++) {
              mt->Bond[bi][0] = s0+k; mt->Bond[bi][1] = s0+k+1;  mt->Bond[bi][2] = 0; bi++;
            }
            mt->Bond[bi][0] = s0+n-1; mt->Bond[bi][1] = j_end;   mt->Bond[bi][2] = 0; bi++;
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
    b->Type     = (i < N_junc) ? bt_J : bt_S;
    b->Molecule = 0;
    mol->Bead[i]      = i;
    System.BeadCoor[i] = i;
  } //}}}

  FinishSystem(&System);

  // assign straight-line coordinates //{{{
  /*
   * Junction positions: (cx + frac_pos[t][0], ...) * A
   * Strand bead k (0..n-1) for strand (cell, tA, b):
   *   pos = j_start_pos + (k+1)/(n+1) * bond_vec[b] * A
   * This places beads evenly between the two junctions.
   * Bond length check: |bond_vec[b]| * A / (n+1) = (sqrt(3)/4)*A/(n+1) = l ✓
   */
  for (int iz = 0; iz < Nz; iz++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int ix = 0; ix < Nx; ix++) {
        int cell_idx = iz*Ny*Nx + iy*Nx + ix;
        // junctions
        for (int t = 0; t < 8; t++) {
          BEAD *b = &System.Bead[cell_idx * 8 + t];
          b->Position.v[0] = (ix + frac_pos[t][0]) * A;
          b->Position.v[1] = (iy + frac_pos[t][1]) * A;
          b->Position.v[2] = (iz + frac_pos[t][2]) * A;
        }
        // strand beads
        for (int tA = 0; tA < 4; tA++) {
          double ox = (ix + frac_pos[tA][0]) * A;
          double oy = (iy + frac_pos[tA][1]) * A;
          double oz = (iz + frac_pos[tA][2]) * A;
          for (int b = 0; b < 4; b++) {
            int s0 = N_junc + (cell_idx * 16 + tA * 4 + b) * n;
            double bvx = bond_vec[b][0] * A;
            double bvy = bond_vec[b][1] * A;
            double bvz = bond_vec[b][2] * A;
            for (int k = 0; k < n; k++) {
              double f = (double)(k + 1) / (n + 1);
              System.Bead[s0 + k].Position.v[0] = ox + f * bvx;
              System.Bead[s0 + k].Position.v[1] = oy + f * bvy;
              System.Bead[s0 + k].Position.v[2] = oz + f * bvz;
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
