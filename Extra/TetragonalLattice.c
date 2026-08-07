#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "Generate a 3D fully interconnected polymer network lattice with tetragonal "
  "symmetry (a=b!=c). Junction beads ('J') sit on a regular Nx x Ny x Nz grid. "
  "Adjacent junctions are connected by strands of intermediate beads ('S') "
  "along all three axes with periodic boundary conditions. In-plane (x,y) "
  "strands contain nxy beads and z-strands contain nz beads, giving lattice "
  "parameters a=b=(nxy+1)*l and c=(nz+1)*l; set nxy!=nz for tetragonal, "
  "nxy=nz for cubic. The output is a single-frame coordinate file with all "
  "beads in one molecule (the network).",

  "Usage: TetragonalLattice <output> [options]",
  .args = 1,
  .all = 8,
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_VERSION],
  {"<output>", nullptr, "output coordinate file", OPT_ARG},
  {"-n", "<nxy> <nz>", "strand beads in xy and z directions (default: 3 5)",
    OPT_EXTRA},
  {"-box", "<int> <int> <int>", "junction grid Nx Ny Nz, each >=2 "
    "(default: 4 4 4)", OPT_EXTRA},
  {"-l", "<float>", "bond length between consecutive beads (default: 1.0)",
    OPT_EXTRA},
  {nullptr}
}; //}}}

// structure for options //{{{
struct OPT {
  int n[2];       // -n nxy nz
  int box[3];     // -box Nx Ny Nz
  double bond_l;  // -l
}; //}}}

int main(int argc, char *argv[]) {

  // command-line arguments //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <output> - output coordinate file
  FILE_TYPE fw_coor;
  snprintf(fw_coor.name, LINE, "%s", argv[++count]);
  fw_coor.type = CoordinateFileType(fw_coor.name);

  SYS_FILES in = InitSysFiles;
  COMMON_OPT commons = CommonOptions(argc, argv, in);

  // -n: strand beads in xy plane and along z
  opt.n[0] = 3;
  opt.n[1] = 5;
  TwoNumbersOption(argc, argv, "-n", opt.n, 'i');
  // -box: junction grid dimensions
  opt.box[0] = opt.box[1] = opt.box[2] = 4;
  ThreeNumbersOption(argc, argv, "-box", opt.box, 'i');
  // -l: bond length
  opt.bond_l = 1.0;
  OneNumberOption(argc, argv, "-l", &opt.bond_l, 'd'); //}}}

  // validate options //{{{
  if (opt.n[0] < 1 || opt.n[1] < 1) {
    snprintf(ERROR_MSG, LINE, "both values of -n must be >= 1 (got %d %d)",
             opt.n[0], opt.n[1]);
    PrintErrorOption("-n");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  for (int dd = 0; dd < 3; dd++) {
    if (opt.box[dd] < 2) {
      snprintf(ERROR_MSG, LINE, "-box dimensions must be >= 2 (got %d)",
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
  const int nxy = opt.n[0]; // strand beads per x- or y-strand
  const int nz  = opt.n[1]; // strand beads per z-strand
  const double l = opt.bond_l;

  const int N_junc    = Nx * Ny * Nz;
  const int N_sxy     = 2 * N_junc * nxy; // total xy-strand beads
  const int N_sz      = N_junc * nz;       // total z-strand beads
  const int N_strand  = N_sxy + N_sz;
  const int N_total   = N_junc + N_strand;
  const int N_bonds   = N_junc * (2*(nxy+1) + (nz+1));

  const double a = (nxy + 1) * l; // in-plane lattice parameter
  const double c = (nz  + 1) * l; // out-of-plane lattice parameter
  const double Lx = Nx * a;
  const double Ly = Ny * a;
  const double Lz = Nz * c; //}}}

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

  // allocate bead arrays
  System.Bead = realloc(System.Bead, N_total * sizeof *System.Bead);
  System.BeadCoor = realloc(System.BeadCoor, N_total * sizeof *System.BeadCoor);

  // bead types: J (junction) and S (strand)
  NewBeadType(&System.BeadType, &Count->BeadType, "J", 0, 1.0, 0.5);
  System.BeadType[0].Number = N_junc;
  NewBeadType(&System.BeadType, &Count->BeadType, "S", 0, 1.0, 0.5);
  System.BeadType[1].Number = N_strand;

  // bond type
  System.BondType[0].b = l;

  // molecule type: one molecule "NET" containing all beads
  NewMolType(&System.MoleculeType, &Count->MoleculeType, "NET",
             N_total, N_bonds, 0, 0, 0);
  MOLECULETYPE *mt = &System.MoleculeType[0];
  mt->Number = 1; //}}}

  // bead-type sequence within the molecule //{{{
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
   * Bead index layout  (j0 = iz*Nx*Ny + iy*Nx + ix):
   *   junctions [0 .. N_junc-1]
   *   x-strands [N_junc+j0*nxy .. N_junc+j0*nxy+nxy-1]
   *   y-strands [N_junc+N_junc*nxy+j0*nxy .. N_junc+N_junc*nxy+j0*nxy+nxy-1]
   *   z-strands [N_junc+2*N_junc*nxy+j0*nz .. +j0*nz+nz-1]
   *
   * Each strand j0->j_next has (n+1) bonds:
   *   j0 -- s[0] -- s[1] -- ... -- s[n-1] -- j_next
   */
  int bi = 0;
  for (int iz = 0; iz < Nz; iz++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int ix = 0; ix < Nx; ix++) {
        int j0 = iz*Nx*Ny + iy*Nx + ix;

        // x-direction strand: j0 -> j((ix+1)%Nx, iy, iz)
        int jxp = iz*Nx*Ny + iy*Nx + (ix+1)%Nx;
        int sx0 = N_junc + j0*nxy;
        mt->Bond[bi][0] = j0;
        mt->Bond[bi][1] = sx0;
        mt->Bond[bi][2] = 0;
        bi++;
        for (int k = 0; k < nxy-1; k++) {
          mt->Bond[bi][0] = sx0+k;
          mt->Bond[bi][1] = sx0+k+1;
          mt->Bond[bi][2] = 0;
          bi++;
        }
        mt->Bond[bi][0] = sx0+nxy-1;
        mt->Bond[bi][1] = jxp;
        mt->Bond[bi][2] = 0;
        bi++;

        // y-direction strand: j0 -> j(ix, (iy+1)%Ny, iz)
        int jyp = iz*Nx*Ny + ((iy+1)%Ny)*Nx + ix;
        int sy0 = N_junc + N_junc*nxy + j0*nxy;
        mt->Bond[bi][0] = j0;
        mt->Bond[bi][1] = sy0;
        mt->Bond[bi][2] = 0;
        bi++;
        for (int k = 0; k < nxy-1; k++) {
          mt->Bond[bi][0] = sy0+k;
          mt->Bond[bi][1] = sy0+k+1;
          mt->Bond[bi][2] = 0;
          bi++;
        }
        mt->Bond[bi][0] = sy0+nxy-1;
        mt->Bond[bi][1] = jyp;
        mt->Bond[bi][2] = 0;
        bi++;

        // z-direction strand: j0 -> j(ix, iy, (iz+1)%Nz)
        int jzp = ((iz+1)%Nz)*Nx*Ny + iy*Nx + ix;
        int sz0 = N_junc + 2*N_junc*nxy + j0*nz;
        mt->Bond[bi][0] = j0;
        mt->Bond[bi][1] = sz0;
        mt->Bond[bi][2] = 0;
        bi++;
        for (int k = 0; k < nz-1; k++) {
          mt->Bond[bi][0] = sz0+k;
          mt->Bond[bi][1] = sz0+k+1;
          mt->Bond[bi][2] = 0;
          bi++;
        }
        mt->Bond[bi][0] = sz0+nz-1;
        mt->Bond[bi][1] = jzp;
        mt->Bond[bi][2] = 0;
        bi++;
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
    mol->Bead[i]     = i;
    System.BeadCoor[i] = i;
  } //}}}

  FinishSystem(&System);

  // assign straight-line coordinates //{{{
  /*
   * Strand beads lie on a straight line between the two junctions they connect.
   * xy-strands use spacing 'a', z-strands use spacing 'c'.
   */
  for (int iz = 0; iz < Nz; iz++) {
    for (int iy = 0; iy < Ny; iy++) {
      for (int ix = 0; ix < Nx; ix++) {
        int j0 = iz*Nx*Ny + iy*Nx + ix;
        double ox = ix * a;
        double oy = iy * a;
        double oz = iz * c;

        // junction position
        System.Bead[j0].Position.v[0] = ox;
        System.Bead[j0].Position.v[1] = oy;
        System.Bead[j0].Position.v[2] = oz;

        // x-strand beads: offset along x
        for (int k = 0; k < nxy; k++) {
          int id = N_junc + j0*nxy + k;
          System.Bead[id].Position.v[0] = ox + (k+1)*l;
          System.Bead[id].Position.v[1] = oy;
          System.Bead[id].Position.v[2] = oz;
        }
        // y-strand beads: offset along y
        for (int k = 0; k < nxy; k++) {
          int id = N_junc + N_junc*nxy + j0*nxy + k;
          System.Bead[id].Position.v[0] = ox;
          System.Bead[id].Position.v[1] = oy + (k+1)*l;
          System.Bead[id].Position.v[2] = oz;
        }
        // z-strand beads: offset along z
        for (int k = 0; k < nz; k++) {
          int id = N_junc + 2*N_junc*nxy + j0*nz + k;
          System.Bead[id].Position.v[0] = ox;
          System.Bead[id].Position.v[1] = oy;
          System.Bead[id].Position.v[2] = oz + (k+1)*l;
        }
      }
    }
  } //}}}

  // info output //{{{
  if (!commons.silent) {
    fprintf(stdout, "Tetragonal polymer network lattice:\n");
    fprintf(stdout, "  Junction grid:      %d x %d x %d = %d\n",
            Nx, Ny, Nz, N_junc);
    fprintf(stdout, "  Strand beads xy/z:  %d / %d  (bond length: %g)\n",
            nxy, nz, l);
    fprintf(stdout, "  Lattice param a=b:  %g  (= %d * %g)\n", a, nxy+1, l);
    fprintf(stdout, "  Lattice param c:    %g  (= %d * %g)\n", c, nz+1,  l);
    fprintf(stdout, "  Box size:           %g x %g x %g\n", Lx, Ly, Lz);
    fprintf(stdout, "  Total beads:        %d  (J: %d, S: %d)\n",
            N_total, N_junc, N_strand);
    fprintf(stdout, "  Total bonds:        %d\n", N_bonds);
  }
  if (commons.verbose) {
    VerboseOutput(System);
  } //}}}

  InitOutputCoorFile(fw_coor, System, argc, argv);
  WriteTimestepAll(fw_coor, System, 1, argc, argv);

  FreeSystem(&System);

  return 0;
}
