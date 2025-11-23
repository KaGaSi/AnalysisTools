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
    b->Position.v[dd] = 0;
    b->Velocity.v[dd] = 0;
    b->Force.v[dd] = 0;
  }
  for (int dd = 0; dd < 6; dd++) {
    b->Extra[dd] = HIGHNUM;
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
// free the agg_picker struct //{{{
void FreeAggPicker(AGG_PICKER *opt) {
  free(opt->m);
  free(opt->only);
  free(opt->x);
} //}}}
void FreeSystem(SYSTEM *System) { //{{{
  free(System->MoleculeCoor);
  free(System->BeadCoor);
  free(System->Bonded);
  free(System->BondedCoor);
  free(System->Unbonded);
  free(System->UnbondedCoor);
  free(System->Bead);
  for (int i = 0; i < System->Count.BeadType; i++) {
    if (System->BeadType[i].Number > 0) {
      free(System->BeadType[i].Index);
    }
  }
  free(System->BeadType);
  for (int i = 0; i < System->Count.Molecule; i++) {
    free(System->Molecule[i].Bead);
  }
  free(System->Molecule);
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    FreeMoleculeType(&System->MoleculeType[i]);
  }
  free(System->MoleculeType);
  free(System->BondType);
  free(System->AngleType);
  free(System->DihedralType);
  free(System->ImproperType);
}; //}}}
void FreeMoleculeType(MOLECULETYPE *MoleculeType) { //{{{
  FreeMoleculeTypeEssentials(MoleculeType);
  if (MoleculeType->nBTypes > 0) {
    free(MoleculeType->BType);
  }
  if (MoleculeType->Number > 0) {
    free(MoleculeType->Index);
  }
} //}}}
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType) { //{{{
  free(MoleculeType->Bead);
  if (MoleculeType->nBonds > 0) {
    free(MoleculeType->Bond);
  }
  if (MoleculeType->nAngles > 0) {
    free(MoleculeType->Angle);
  }
  if (MoleculeType->nDihedrals > 0) {
    free(MoleculeType->Dihedral);
  }
  if (MoleculeType->nImpropers > 0) {
    free(MoleculeType->Improper);
  }
} //}}}
void FreeAggregate(COUNT Count, AGGREGATE *Aggregate) { //{{{
  for (int i = 0; i < Count.Molecule; i++) {
    free(Aggregate[i].Molecule);
    free(Aggregate[i].Bead);
  }
  free(Aggregate);
} //}}}
// realloc some System.*{,Coor} arrays //{{{
void ReallocBead(SYSTEM *System) {
  COUNT *Count = &System->Count;
  System->Bead = s_realloc(System->Bead, sizeof *System->Bead * Count->Bead);
  System->BeadCoor = s_realloc(System->BeadCoor,
                               sizeof *System->BeadCoor * Count->Bead);
}
void ReallocBonded(SYSTEM *System) {
  COUNT *Count = &System->Count;
  if (Count->Bonded > 0) {
    System->Bonded = s_realloc(System->Bonded,
                               sizeof *System->Bonded * Count->Bonded);
    System->BondedCoor = s_realloc(System->BondedCoor,
                                   sizeof *System->BondedCoor * Count->Bonded);
  }
}
void ReallocUnbonded(SYSTEM *System) {
  COUNT *Count = &System->Count;
  if (Count->Unbonded > 0) {
    System->Unbonded = s_realloc(System->Unbonded,
                                 sizeof *System->Unbonded * Count->Unbonded);
    System->UnbondedCoor = s_realloc(System->UnbondedCoor,
                                     sizeof *System->UnbondedCoor *
                                     Count->Unbonded);
  }
}
void ReallocMolecule(SYSTEM *System) {
  COUNT *Count = &System->Count;
  if (Count->Molecule > 0) {
    System->Molecule = s_realloc(System->Molecule,
                                sizeof *System->Molecule * Count->Molecule);
    System->MoleculeCoor = s_realloc(System->MoleculeCoor, Count->Molecule *
                                    sizeof *System->MoleculeCoor);
  }
}
//}}}
