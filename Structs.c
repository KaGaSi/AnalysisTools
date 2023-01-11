#include "Structs.h"

void InitBeadType(BEADTYPE *bt) { //{{{
  bt->Name[0] = '\0';
  bt->Number = 0;
  bt->Charge = CHARGE;
  bt->Mass = MASS;
  bt->Radius = RADIUS;
  // bt->Index = calloc(1, sizeof *bt->Index);
} //}}}
void InitBead(BEAD *b) { //{{{
  b->Type = -1;
  b->Molecule = -1;
  b->InTimestep = false;
  b->Position.x = 0;
  b->Position.y = 0;
  b->Position.z = 0;
  b->Velocity.x = 0;
  b->Velocity.y = 0;
  b->Velocity.z = 0;
  b->Force.x = 0;
  b->Force.y = 0;
  b->Force.z = 0;
} //}}}
void InitMoleculeType(MOLECULETYPE *mt) { //{{{
  mt->Name[0] = '\0';
  mt->Number = 0;
  mt->nBeads = 0;
  mt->nBonds = 0;
  mt->nAngles = 0;
  mt->nDihedrals = 0;
  mt->nImpropers = 0;
  mt->nBTypes = 0;
  mt->Mass = MASS;
  mt->Charge = CHARGE;
  mt->InVcf = false;
  // mt->BType = calloc(1, sizeof *mt->BType);
} //}}}
void InitMolecule(MOLECULE *mol) { //{{{
  mol->Type = -1;
  mol->Index = -1;
} //}}}
void InitSystem(SYSTEM *System) { //{{{
  System->Box = InitBox;
  System->Count = InitCount;
  System->BeadType =     calloc(1, sizeof System->BeadType);
  System->Bead =         calloc(1, sizeof System->Bead);
  System->MoleculeType = calloc(1, sizeof System->MoleculeType);
  System->Molecule =     calloc(1, sizeof System->Molecule);
  System->BondType =     calloc(1, sizeof System->BondType);
  System->AngleType =    calloc(1, sizeof System->AngleType);
  System->DihedralType = calloc(1, sizeof System->DihedralType);
  System->ImproperType = calloc(1, sizeof System->ImproperType);
  System->Index_mol =    calloc(1, sizeof System->Index_mol);
  System->Bonded =       calloc(1, sizeof System->Bonded);
  System->BondedCoor =   calloc(1, sizeof System->BondedCoor);
  System->Unbonded =     calloc(1, sizeof System->Unbonded);
  System->UnbondedCoor = calloc(1, sizeof System->UnbondedCoor);
  System->BeadCoor =     calloc(1, sizeof System->BeadCoor);
} //}}}
