#include "Structs.h"

void InitBeadType(BEADTYPE *bt) { //{{{
  bt->Number = 0;
  bt->InCoor = 0;
  bt->Charge = CHARGE;
  bt->Mass = MASS;
  bt->Radius = RADIUS;
  bt->Flag = false;
} //}}}
void InitBead(BEAD *b) { //{{{
  b->Type = -1;
  b->Molecule = -1;
  b->Aggregate = -1;
  for (int dd = 0; dd < 3; dd++) {
    b->Position[dd] = 0;
    b->Velocity[dd] = 0;
    b->Force[dd] = 0;
  }
  b->InTimestep = false;
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
  mt->Flag = false;
} //}}}
void InitMolecule(MOLECULE *mol) { //{{{
  mol->Type = -1;
  mol->Index = -1;
  mol->Aggregate = -1;
  mol->InTimestep = false;
} //}}}
void InitSystem(SYSTEM *System) { //{{{
  System->Box = InitBox;
  System->Count = InitCount;
  System->BeadType =     calloc(1, sizeof(struct BeadType));
  System->Bead =         calloc(1, sizeof(struct Bead));
  System->MoleculeType = calloc(1, sizeof(struct MoleculeType));
  System->Molecule =     calloc(1, sizeof(struct Molecule));
  System->BondType =     calloc(1, sizeof(struct Params));
  System->AngleType =    calloc(1, sizeof(struct Params));
  System->DihedralType = calloc(1, sizeof(struct Params));
  System->ImproperType = calloc(1, sizeof(struct Params));
  System->MoleculeCoor = calloc(1, sizeof *System->MoleculeCoor);
  System->Bonded =       calloc(1, sizeof *System->Bonded);
  System->BondedCoor =   calloc(1, sizeof *System->BondedCoor);
  System->Unbonded =     calloc(1, sizeof *System->Unbonded);
  System->UnbondedCoor = calloc(1, sizeof *System->UnbondedCoor);
  System->BeadCoor =     calloc(1, sizeof *System->BeadCoor);
} //}}}
void InitAggregate(SYSTEM System, AGGREGATE **Aggregate) { //{{{
  COUNT *Count = &System.Count;
  *Aggregate = malloc(Count->Molecule * sizeof **Aggregate);
  for (int i = 0; i < Count->Molecule; i++) {
    (*Aggregate)[i].nMolecules = 0;
    (*Aggregate)[i].nBeads = 0;
    (*Aggregate)[i].Molecule = calloc(1, sizeof *Aggregate[i]->Molecule);
    (*Aggregate)[i].Bead = calloc(1, sizeof *Aggregate[i]->Bead);
  }
} //}}}
void ReInitAggregate(SYSTEM System, AGGREGATE *Aggregate) { //{{{
  COUNT *Count = &System.Count;
  for (int i = 0; i < Count->Molecule; i++) {
    Aggregate[i].nMolecules = 0;
    Aggregate[i].nBeads = 0;
    Aggregate[i].Molecule = s_realloc(Aggregate[i].Molecule,
                                      1 * sizeof *Aggregate[i].Molecule);
    Aggregate[i].Bead = s_realloc(Aggregate[i].Bead,
                                  1 * sizeof *Aggregate[i].Bead);
  }
} //}}}
