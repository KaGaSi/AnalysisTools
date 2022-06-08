#include "Structs.h"

void InitBeadType(BEADTYPE *bt) {
  (*bt).Name[0] = '\0';
  (*bt).Number = 0;
  (*bt).Charge = CHARGE;
  (*bt).Mass = MASS;
  (*bt).Radius = RADIUS;
//(*bt).Index = malloc(sizeof *(*bt).Index); // TODO only if Number > 0?
//(*bt).Index[0] = -1;
}
void InitBead(BEAD *b) {
  (*b).Type = -1;
  (*b).Molecule = -1;
  (*b).InTimestep = false;
}
void InitMoleculeType(MOLECULETYPE *mt) {
  (*mt).Name[0] = '\0';
  (*mt).Number = 0;
  (*mt).nBeads = 0;
  (*mt).nBonds = 0;
  (*mt).nAngles = 0;
  (*mt).nDihedrals = 0;
  (*mt).nBTypes = 0;
  (*mt).Mass = MASS;
  (*mt).Charge = CHARGE;
  (*mt).InVcf = false;
//(*mt).Index = malloc(sizeof *(*mt).Index); // TODO only if Number > 0?
}
void InitMolecule(MOLECULE *mol) {
  (*mol).Type = -1;
  (*mol).Index = -1;
//(*mol).Bead = malloc(sizeof *(*mol).Bead); // TODO only if its bt.Number > 0?
//(*mol).Bead[0] = -1;
}
// TODO add stuff?
// TODO InitStuff() - possibly initializing twice (second in Read.c)
// TODO why allocate it all? ...not in InitMolecule(), etc. Should it be?
void InitSystem(SYSTEM *System) {
  System->Box = InitBox;
  System->Count = InitCount;
  System->BeadType = calloc(1, sizeof *(*System).BeadType);
  System->Bead = calloc(1, sizeof *(*System).Bead);
  System->MoleculeType = calloc(1, sizeof *(*System).MoleculeType);
  System->Molecule = calloc(1, sizeof *(*System).Molecule);
  System->BondType = calloc(1, sizeof *(*System).BondType);
  System->AngleType = calloc(1, sizeof *(*System).AngleType);
  System->DihedralType = calloc(1, sizeof *(*System).DihedralType);
  System->Index_mol = calloc(1, sizeof *(*System).Index_mol);
  System->Bonded = calloc(1, sizeof *(*System).Bonded);
  System->BondedCoor = calloc(1, sizeof *(*System).BondedCoor);
  System->Unbonded = calloc(1, sizeof *(*System).Unbonded);
  System->UnbondedCoor = calloc(1, sizeof *(*System).UnbondedCoor);
  System->BeadCoor = calloc(1, sizeof *(*System).BeadCoor);
//InitBeadType(&(*System).BeadType[0]);
//InitBead(&(*System).Bead[0]);
//InitMoleculeType(&(*System).MoleculeType[0]);
//InitMolecule(&(*System).Molecule[0]);
} //}}}
