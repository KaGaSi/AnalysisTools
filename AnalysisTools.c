#include "AnalysisTools.h"

// STATIC DEFINITIONS
// functions to transform to/from fractional coordinates
static VECTOR ToFractional(VECTOR coor, BOX Box);
static void ToFractionalCoor(SYSTEM *System);
static VECTOR FromFractional(VECTOR coor, BOX Box);
static void FromFractionalCoor(SYSTEM *System);
// remove pbc for molecules by joining the molecules
static void RemovePBCMolecules(SYSTEM *System);
// restore pbc by wrapping all coordinates inside the simulation box
static void RestorePBC(SYSTEM *System);
// Copy System structure; assumes new unallocated SYSTEM
static SYSTEM CopySystem(SYSTEM S_in);
// calculate geometric centre for a list of beads
static VECTOR GeomCentre(int n, int *list, BEAD *Bead);
// calculate centre of mass for a list of beads (UNUSED)
// static VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType);

// STATIC IMPLEMENTATIONS
// transform to/from fractional coordinates
/* HOW TO CALCULATE DISTANCE IN TRICLINIC SYSTEM //{{{
//VECTOR dist;
//dist.x = (*Bead)[0].Position.x - (*Bead)[10].Position.x;
//dist.y = (*Bead)[0].Position.y - (*Bead)[10].Position.y;
//dist.z = (*Bead)[0].Position.z - (*Bead)[10].Position.z;
//printf("dist1 = (%lf, %lf, %lf) = %lf\n", dist.x, dist.y, dist.z,
sqrt(SQR(dist.x)+SQR(dist.y)+SQR(dist.z)));

//VECTOR new;
//new.x = Box.transform[0][0] * dist.x +
//        Box.transform[0][1] * dist.y +
//        Box.transform[0][2] * dist.z;
//new.y = Box.transform[1][0] * dist.x +
//        Box.transform[1][1] * dist.y +
//        Box.transform[1][2] * dist.z;
//new.z = Box.transform[2][0] * dist.x +
//        Box.transform[2][1] * dist.y +
//        Box.transform[2][2] * dist.z;
//dist.x = new.x / a;
//dist.y = new.y / b;
//dist.z = new.z / c;
//printf("dist2 = (%lf, %lf, %lf) = %lf\n", dist.x, dist.y, dist.z,
sqrt(SQR(dist.x)+SQR(dist.y)+SQR(dist.z)));
*/ //}}}
static VECTOR ToFractional(VECTOR coor, BOX Box) { //{{{
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    VECTOR new = {0, 0, 0};
    new.x = Box.inverse[0][0] * coor.x + Box.inverse[0][1] * coor.y +
            Box.inverse[0][2] * coor.z;
    new.y = Box.inverse[1][0] * coor.x + Box.inverse[1][1] * coor.y +
            Box.inverse[1][2] * coor.z;
    new.z = Box.inverse[2][0] * coor.x + Box.inverse[2][1] * coor.y +
            Box.inverse[2][2] * coor.z;
    coor.x = new.x *Box.Length.x;
    coor.y = new.y *Box.Length.y;
    coor.z = new.z *Box.Length.z;
  }
  return coor;
} //}}}
static void ToFractionalCoor(SYSTEM *System) { //{{{
  if (System->Box.alpha != 90 || System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.Bead; i++) {
      BEAD *bead = &System->Bead[i];
      if (bead->InTimestep) {
        bead->Position = ToFractional(bead->Position, System->Box);
      }
    }
  }
} //}}}
static VECTOR FromFractional(VECTOR coor, BOX Box) { //{{{
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    coor.x /= Box.Length.x;
    coor.y /= Box.Length.y;
    coor.z /= Box.Length.z;
    VECTOR new = {0, 0, 0};
    new.x = Box.transform[0][0] * coor.x + Box.transform[0][1] * coor.y +
            Box.transform[0][2] * coor.z;
    new.y = Box.transform[1][0] * coor.x + Box.transform[1][1] * coor.y +
            Box.transform[1][2] * coor.z;
    new.z = Box.transform[2][0] * coor.x + Box.transform[2][1] * coor.y +
            Box.transform[2][2] * coor.z;
    coor.x = new.x;
    coor.y = new.y;
    coor.z = new.z;
  }
  return coor;
} //}}}
static void FromFractionalCoor(SYSTEM *System) { //{{{
  if (System->Box.alpha != 90 || System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.Bead; i++) {
      BEAD *bead = &System->Bead[i];
      if (bead->InTimestep) {
        bead->Position = FromFractional(bead->Position, System->Box);
      }
    }
  }
} //}}}
bool ConnectedMolecule(MOLECULE mol, SYSTEM System) {
  MOLECULETYPE *mt = &System.MoleculeType[mol.Type];
  if (mt->nBonds == 0) {
    return false;
  }
  bool *connected = calloc(mt->nBeads, sizeof *connected);
  // find first bead in the first bond present in the timestep
  for (int i = 0; i < mt->nBonds; i++) {
    int id = mt->Bond[i][0];
    BEAD *b = &System.Bead[mol.Bead[id]];
    if (b->InTimestep) {
      connected[id] = true;
      break;
    }
  }
  for (int i = 0; i < mt->nBonds; i++) {
    int id1 = mt->Bond[i][0],
        id2 = mt->Bond[i][1];
    BEAD *b_1 = &System.Bead[mol.Bead[id1]],
         *b_2 = &System.Bead[mol.Bead[id2]];
    if (b_1->InTimestep && b_2->InTimestep && connected[id1]) {
      connected[id2] = true;
    }
  }
  for (int i = 0; i < mt->nBeads; i++) {
    BEAD *b = &System.Bead[mol.Bead[i]];
    if (b->InTimestep && !connected[i]) {
      free(connected);
      return false;
    }
  }
  free(connected);
  return true;
}
/*
 * TODO what about if the molecule has bonds, but is in more 'pieces'; should
 *      it really just try joining it 1000 times when we know it's impossible?
 *      There should be some check and then the separate pieces should be
 *      connected (with a warning issued).
 *      ...working on it - now, unconnected molecules aren't joined
 *      TODO: maybe split molecule to connected pieces, connect those and place
 *            their centres of mass nearest each other
 *      TODO: check the ConnectedMolecule() and possible make the joining
 *            algorithm along the same lines (go over bond after bond instead of
 *            that while loop)
 */
// remove pbc for molecules by joining the molecules //{{{
static void RemovePBCMolecules(SYSTEM *System) {
  BOX *box = &System->Box;
  // go through all molecules
  for (int i = 0; i < System->Count.Molecule; i++) {
    MOLECULE *mol_i = &System->Molecule[i];
    MOLECULETYPE *mt_i = &System->MoleculeType[mol_i->Type];
    if (!ConnectedMolecule(*mol_i, *System)) {
      snprintf(ERROR_MSG, LINE, "molecule %s%d%s (%s%s%s) does not have all "
              "beads connected", ErrYellow(), mol_i->Index, ErrCyan(),
              ErrYellow(), mt_i->Name, ErrCyan());
      PrintWarning();
    } else {
      // int type = System->Molecule[i].Type;
      // do nothing if the molecule has no bonds
      if (mt_i->nBonds == 0) {
        continue;
      }
      for (int j = 0; j < mt_i->nBeads; j++) {
        int id = mol_i->Bead[j];
        System->Bead[id].Flag = false; // no beads moved yet
      }
      // first bead in the first bond is considered moved
      for (int j = 0; j < mt_i->nBonds; j++) {
        int id1 = mt_i->Bond[j][0], // bead position in the molecule
            id2 = mt_i->Bond[j][1]; //
        id1 = mol_i->Bead[id1];     // bead index in the System
        id2 = mol_i->Bead[id2];     //
        if (System->Bead[id1].InTimestep && System->Bead[id2].InTimestep) {
          System->Bead[id1].Flag = true;
          break;
        }
      }
      bool done = false;
      int test = 0; // if too many loops, just leave the loop with error
      while (!done && test < 1000) {
        for (int j = 0; j < mt_i->nBonds; j++) {
          int id1 = mt_i->Bond[j][0], // bead position in the molecule
              id2 = mt_i->Bond[j][1]; //
          id1 = mol_i->Bead[id1];     // bead index in the System
          id2 = mol_i->Bead[id2];     //
          BEAD *b_1 = &System->Bead[id1], *b_2 = &System->Bead[id2];
          if (b_1->InTimestep && b_2->InTimestep) {
            // move id1, if id2 is moved already
            if (!b_1->Flag && b_2->Flag) {
              VECTOR dist = Distance(b_2->Position, b_1->Position, box->Length);
              b_1->Position.x = b_2->Position.x - dist.x;
              b_1->Position.y = b_2->Position.y - dist.y;
              b_1->Position.z = b_2->Position.z - dist.z;
              b_1->Flag = true;
              // move id2, if id1 was moved already
            } else if (b_1->Flag && !b_2->Flag) {
              VECTOR dist = Distance(b_1->Position, b_2->Position, box->Length);
              b_2->Position.x = b_1->Position.x - dist.x;
              b_2->Position.y = b_1->Position.y - dist.y;
              b_2->Position.z = b_1->Position.z - dist.z;
              b_2->Flag = true;
            }
          }
        }

        // break while loop if all beads have moved
        done = true;
        for (int j = 1; j < mt_i->nBeads; j++) {
          int id = mol_i->Bead[j];
          if (System->Bead[id].InTimestep && !System->Bead[id].Flag) {
            done = false;
            break;
          }
        }
        test++;
      }
      if (test == 1000) {
        snprintf(ERROR_MSG, LINE,
                 "unable to join molecule %s%s%s (resid "
                 "%s%d%s)\n     Either all beads are not connected or "
                 "some beads are missing from the timestep",
                 ErrYellow(), mt_i->Name, ErrCyan(), ErrYellow(), i + 1,
                 ErrCyan());
        PrintWarning();
      }

      // put molecule's geometric centre into the simulation box //{{{
      VECTOR cog = GeomCentre(mt_i->nBeads, mol_i->Bead, System->Bead);
      // by how many BoxLength's should cog be moved?
      // for distant molecules - it shouldn't happen, but better safe than sorry
      INTVECTOR move;
      move.x = cog.x / box->Length.x;
      move.y = cog.y / box->Length.y;
      move.z = cog.z / box->Length.z;
      if (cog.x < 0) {
        move.x--;
      }
      if (cog.y < 0) {
        move.y--;
      }
      if (cog.z < 0) {
        move.z--;
      }
      for (int j = 0; j < mt_i->nBeads; j++) {
        int bead = mol_i->Bead[j];
        System->Bead[bead].Position.x -= move.x * box->Length.x;
        System->Bead[bead].Position.y -= move.y * box->Length.y;
        System->Bead[bead].Position.z -= move.z * box->Length.z;
      } //}}}
    }
  }
} //}}}
// restore pbc by wrapping all coordinates inside the simulation box //{{{
static void RestorePBC(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    int id = System->BeadCoor[i];
    BEAD *bead = &System->Bead[id];
    BOX *box = &System->Box;
    // x direction
    while (bead->Position.x >= box->Length.x) {
      bead->Position.x -= box->Length.x;
    }
    while (bead->Position.x < 0) {
      bead->Position.x += box->Length.x;
    }
    // y direction
    while (bead->Position.y >= box->Length.y) {
      bead->Position.y -= box->Length.y;
    }
    while (bead->Position.y < 0) {
      bead->Position.y += box->Length.y;
    }
    // z direction
    while (bead->Position.z >= box->Length.z) {
      bead->Position.z -= box->Length.z;
    }
    while (bead->Position.z < 0) {
      bead->Position.z += box->Length.z;
    }
  }
} //}}}
// copy System structure; assumes new unallocated SYSTEM //{{{
static SYSTEM CopySystem(SYSTEM S_in) {
  SYSTEM S_out;
  InitSystem(&S_out);
  S_out.Box = S_in.Box;
  if (S_in.Count.Bead > 0) {
    S_out.Count = S_in.Count;
    // BeadType //{{{
    if (S_out.Count.BeadType > 0) {
      S_out.BeadType = realloc(S_out.BeadType,
                               sizeof *S_out.BeadType * S_out.Count.BeadType);
      for (int i = 0; i < S_out.Count.BeadType; i++) {
        S_out.BeadType[i] = S_in.BeadType[i];
        if (S_out.BeadType[i].Number > 0) {
          S_out.BeadType[i].Index = malloc(sizeof *S_out.BeadType[i].Index *
                                           S_out.BeadType[i].Number);
          for (int j = 0; j < S_out.BeadType[i].Number; j++) {
            S_out.BeadType[i].Index[j] = S_in.BeadType[i].Index[j];
          }
        }
      }
    } else {
      strcpy(ERROR_MSG, "no bead types to copy; should never happen!");
      PrintWarning();
      return S_out;
    } //}}}
    // Bead & BeadCoor //{{{
    S_out.Bead = realloc(S_out.Bead, sizeof(BEAD) * S_out.Count.Bead);
    S_out.BeadCoor =
        realloc(S_out.BeadCoor, sizeof *S_out.BeadCoor * S_out.Count.Bead);
    memcpy(S_out.Bead, S_in.Bead, sizeof(BEAD) * S_in.Count.Bead);
    memcpy(S_out.BeadCoor, S_in.BeadCoor,
           sizeof *S_in.BeadCoor * S_in.Count.Bead); //}}}
    // Bonded & BondedCoor //{{{
    if (S_out.Count.Bonded > 0) {
      S_out.Bonded =
          realloc(S_out.Bonded, sizeof *S_out.Bonded * S_out.Count.Bonded);
      S_out.BondedCoor = realloc(S_out.BondedCoor,
                                 sizeof *S_out.BondedCoor * S_out.Count.Bonded);
      memcpy(S_out.Bonded, S_in.Bonded,
             sizeof *S_in.Bonded * S_in.Count.Bonded);
      memcpy(S_out.BondedCoor, S_in.BondedCoor,
             sizeof *S_in.BondedCoor * S_in.Count.Bonded);
    } //}}}
    // Unbonded & UnbondedCoor //{{{
    if (S_out.Count.Unbonded > 0) {
      S_out.Unbonded = realloc(S_out.Unbonded,
                               sizeof *S_out.Unbonded * S_out.Count.Unbonded);
      S_out.UnbondedCoor =
          realloc(S_out.UnbondedCoor,
                  sizeof *S_out.UnbondedCoor * S_out.Count.Unbonded);
      memcpy(S_out.Unbonded, S_in.Unbonded,
             sizeof *S_in.Unbonded * S_in.Count.Unbonded);
      memcpy(S_out.UnbondedCoor, S_in.UnbondedCoor,
             sizeof *S_in.UnbondedCoor * S_in.Count.Unbonded);
    } //}}}
    // MoleculeType //{{{
    if (S_out.Count.MoleculeType > 0) {
      S_out.MoleculeType = realloc(
          S_out.MoleculeType, sizeof(MOLECULETYPE) * S_out.Count.MoleculeType);
      for (int i = 0; i < S_out.Count.MoleculeType; i++) {
        S_out.MoleculeType[i] = CopyMoleculeType(S_in.MoleculeType[i]);
      }
    } //}}}
    // Molecule & Index_mol //{{{
    if (S_out.Count.Molecule > 0) {
      S_out.Molecule =
          realloc(S_out.Molecule, sizeof(MOLECULE) * S_out.Count.Molecule);
      for (int i = 0; i < S_out.Count.Molecule; i++) {
        S_out.Molecule[i] = S_in.Molecule[i];
        // Molecule[].Bead array
        int type = S_out.Molecule[i].Type;
        if (S_out.MoleculeType[type].nBeads > 0) {
          S_out.Molecule[i].Bead = malloc(sizeof *S_out.Molecule[i].Bead *
                                          S_out.MoleculeType[type].nBeads);
          memcpy(S_out.Molecule[i].Bead, S_in.Molecule[i].Bead,
                 sizeof *S_in.Molecule[i].Bead *
                     S_in.MoleculeType[type].nBeads);
        }
      }
      // Index_mol
      S_out.Index_mol = realloc(S_out.Index_mol,
                                sizeof *S_out.Index_mol * S_out.Count.Molecule);
      memcpy(S_out.Index_mol, S_in.Index_mol,
             sizeof *S_in.Index_mol * S_in.Count.Molecule);
    } //}}}
    // BondType //{{{
    if (S_out.Count.BondType > 0) {
      S_out.BondType =
          realloc(S_out.BondType, sizeof(PARAMS) * S_out.Count.BondType);
      memcpy(S_out.BondType, S_in.BondType,
             sizeof *S_out.BondType * S_in.Count.BondType);
    } //}}}
    // AngleType //{{{
    if (S_out.Count.AngleType > 0) {
      S_out.AngleType = realloc(S_out.AngleType, sizeof *S_out.AngleType *
                                                     S_out.Count.AngleType);
      memcpy(S_out.AngleType, S_in.AngleType,
             sizeof *S_out.AngleType * S_in.Count.AngleType);
    } //}}}
    // DihedralType //{{{
    if (S_out.Count.DihedralType > 0) {
      S_out.DihedralType = realloc(S_out.DihedralType,
                                   sizeof(PARAMS) * S_out.Count.DihedralType);
      memcpy(S_out.DihedralType, S_in.DihedralType,
             sizeof *S_out.DihedralType * S_in.Count.DihedralType);
    } //}}}
    // ImproperType //{{{
    if (S_out.Count.ImproperType > 0) {
      S_out.ImproperType = realloc(S_out.ImproperType,
                                   sizeof(PARAMS) * S_out.Count.ImproperType);
      memcpy(S_out.ImproperType, S_in.ImproperType,
             sizeof *S_out.ImproperType * S_in.Count.ImproperType);
    } //}}}
  }
  return S_out;
} //}}}
// calculate geometric centre for a list of beads //{{{
static VECTOR GeomCentre(int n, int *list, BEAD *Bead) {
  VECTOR cog = {0, 0, 0};
  int count = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    if (Bead[id].InTimestep) {
      cog.x += Bead[id].Position.x;
      cog.y += Bead[id].Position.y;
      cog.z += Bead[id].Position.z;
      count++;
    }
  }
  cog.x /= count;
  cog.y /= count;
  cog.z /= count;
  return cog;
} //}}}
// calculate centre of mass for a list of beads (UNUSED) //{{{
// static VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType)
// {
//   VECTOR com = {0, 0, 0};
//   double mass = 0;
//   for (int i = 0; i < n; i++) {
//     int id = list[i];
//     int btype = Bead[id].Type;
//     com.x += Bead[id].Position.x * BeadType[btype].Mass;
//     com.y += Bead[id].Position.y * BeadType[btype].Mass;
//     com.z += Bead[id].Position.z * BeadType[btype].Mass;
//     mass += BeadType[btype].Mass;
//   }
//   com.x /= mass;
//   com.y /= mass;
//   com.z /= mass;
//   return com;
// } //}}}

// Helper functions for dealing with SYSTEM structure
// fill some System arrays and some such
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType) { //{{{
  MoleculeType->nBTypes = 0;
  MoleculeType->BType = malloc(sizeof *MoleculeType->BType * 1);
  for (int j = 0; j < MoleculeType->nBeads; j++) {
    bool new = true;
    for (int k = 0; k < MoleculeType->nBTypes; k++) {
      if (MoleculeType->Bead[j] == MoleculeType->BType[k]) {
        new = false;
        break;
      }
    }
    if (new) {
      int type = MoleculeType->nBTypes++;
      MoleculeType->BType =
          realloc(MoleculeType->BType,
                  sizeof *MoleculeType->BType * MoleculeType->nBTypes);
      MoleculeType->BType[type] = MoleculeType->Bead[j];
    }
  }
} //}}}
void FillMoleculeTypeChargeMass(MOLECULETYPE *mtype, BEADTYPE btype[]) { //{{{
  mtype->Mass = 0;
  mtype->Charge = 0;
  for (int j = 0; j < mtype->nBeads; j++) {
    int bt = mtype->Bead[j];
    // charge
    if (btype[bt].Charge != CHARGE) {
      mtype->Charge += btype[bt].Charge;
    } else {
      mtype->Charge = CHARGE;
    }
    // mass
    if (btype[bt].Mass != MASS) {
      mtype->Mass += btype[bt].Mass;
    } else {
      mtype->Mass = MASS;
    }
  }
} //}}}
void FillBeadTypeIndex(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  // allocate memory for Index arrays
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System->BeadType[i];
    if (bt_i->Number > 0) {
      bt_i->Index = malloc(sizeof *bt_i->Index * bt_i->Number);
    }
  }
  // fill the Index arrays
  int *count_id = calloc(Count->BeadType, sizeof *count_id);
  for (int i = 0; i < Count->Bead; i++) {
    int type = System->Bead[i].Type;
    if (System->BeadType[type].Number < count_id[type]) {
      fprintf(stderr, "...hmm; error count_id[%d]=%d (%d)\n", type,
              count_id[type], System->BeadType[type].Number);
    }
    System->BeadType[type].Index[count_id[type]] = i;
    count_id[type]++;
  }
  free(count_id);
} //}}}
void FillMoleculeTypeIndex(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    mt_i->Index = malloc(sizeof *mt_i->Index * mt_i->Number);
  }
  int *count_id = calloc(Count->MoleculeType, sizeof *count_id);
  for (int i = 0; i < Count->Molecule; i++) {
    int type = System->Molecule[i].Type;
    System->MoleculeType[type].Index[count_id[type]] = i;
    count_id[type]++;
  }
  free(count_id);
} //}}}
void CountBondAngleDihedralImproper(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  Count->Bond = 0;
  Count->Angle = 0;
  Count->Dihedral = 0;
  Count->Improper = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    Count->Bond += mt_i->nBonds * mt_i->Number;
    Count->Angle += mt_i->nAngles * mt_i->Number;
    Count->Dihedral += mt_i->nDihedrals * mt_i->Number;
    Count->Improper += mt_i->nImpropers * mt_i->Number;
  }
} //}}}
// Appends _# to bead/molecule types with the same name
void RenameBeadTypes(SYSTEM *System) { //{{{
  for (int i = 0; i < (System->Count.BeadType - 1); i++) {
    int count = 0;
    for (int j = (i + 1); j < System->Count.BeadType; j++) {
      if (strcmp(System->BeadType[i].Name, System->BeadType[j].Name) == 0) {
        count++;
        char name[BEAD_NAME];
        strncpy(name, System->BeadType[j].Name, BEAD_NAME);
        // shorten name if necessary
        if (count < 10) {
          name[BEAD_NAME - 3] = '\0';
        } else if (count < 100) {
          name[BEAD_NAME - 4] = '\0';
        } else if (count < 1000) {
          name[BEAD_NAME - 5] = '\0';
        }
        if (snprintf(System->BeadType[j].Name, BEAD_NAME, "%s_%d", name,
                     count) < 0) {
          strcpy(ERROR_MSG, "something wrong with snprintf()");
          PrintError();
          exit(1);
        }
      }
    }
  }
} //}}}
void RenameMoleculeTypes(SYSTEM *System) { //{{{
  for (int i = 0; i < (System->Count.MoleculeType - 1); i++) {
    int count = 0;
    for (int j = (i + 1); j < System->Count.MoleculeType; j++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i],
                   *mt_j = &System->MoleculeType[j];
      if (strcmp(mt_i->Name, mt_j->Name) == 0) {
        count++;
        char name[MOL_NAME];
        strncpy(name, mt_j->Name, MOL_NAME);
        // shorten name if necessary
        if (count < 10) {
          name[MOL_NAME - 3] = '\0';
        } else if (count < 100) {
          name[MOL_NAME - 4] = '\0';
        } else if (count < 1000) {
          name[MOL_NAME - 5] = '\0';
        }
        if (snprintf(mt_j->Name, MOL_NAME, "%s_%d", name,
                     count) < 0) {
          strcpy(ERROR_MSG, "something wrong with snprintf()");
          PrintError();
          exit(1);
        }
      }
    }
  }
} //}}}
// get bead indices from bonds/angles/dihedrals/impropers //{{{
/*
 * They all return an array of size 2*n where n is the number of bead indices
 * (2 for bonds, 3 for angles, 4 for dihedrals and impropers); elements 0..n-1
 * contain internal bead indices for the corresponding array in MoleculeType
 * (Bond, Angle, etc.), while elements n..2*n-1 contain bead indices for in
 * Bead.
 */
int *BondIndices(SYSTEM System, int mol, int bond) { //{{{
  static int index[4]; // first are 'real' ids, then the intramolecular ones
  int n = 2, type = System.Molecule[mol].Type;
  // intramolecular id
  for (int i = 0; i < n; i++) {
    index[i + n] = System.MoleculeType[type].Bond[bond][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i + n] < 0 || index[i + n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in a bond; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in bond %s%d%s:%s", ErrYellow(),
            bond, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i + n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i + n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in a bond; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in bond %s%d%s:%s", ErrYellow(), bond, ErrRed(),
            ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *AngleIndices(SYSTEM System, int mol, int angle) { //{{{
  static int index[6]; // first are 'real' ids, then the intramolecular ones
  int n = 3, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i + n] = System.MoleculeType[type].Angle[angle][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i + n] < 0 || index[i + n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an angle; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in angle %s%d%s:%s", ErrYellow(),
            angle, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i + n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i + n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an angle; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in angle %s%d%s:%s", ErrYellow(), angle, ErrRed(),
            ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *DihedralIndices(SYSTEM System, int mol, int dihed) { //{{{
  static int index[8]; // first are 'real' ids, then the intramolecular ones
  int n = 4, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i + n] = System.MoleculeType[type].Dihedral[dihed][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i + n] < 0 || index[i + n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an dihedral; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in dihedral %s%d%s:%s",
            ErrYellow(), dihed, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i + n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i + n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in a dihedral; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in dihedral %s%d%s:%s", ErrYellow(), dihed,
            ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int *ImproperIndices(SYSTEM System, int mol, int improper) { //{{{
  static int index[8]; // first are 'real' ids, then the intramolecular ones
  int n = 4, type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i + n] = System.MoleculeType[type].Improper[improper][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i + n] < 0 || index[i + n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an improper; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),", ErrRed(),
            ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in improper %s%d%s:%s",
            ErrYellow(), improper, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i + n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i + n]];
  }
  // error if wrong 'true' id //{{{
  err = false;
  for (int i = 0; i < n; i++) {
    if (index[i] < 0 || index[i] >= System.Count.Bead) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in a improper; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),", ErrRed(), ErrYellow(), mol,
            ErrRed(), ErrYellow(), System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in improper %s%d%s:%s", ErrYellow(), improper,
            ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
  //}}}
// ChangeMolecules() //{{{
/*
 * Function to add bonds, angles, dihedrals, and/or impropers as well as their
 * types (creating new ones) to a molecule type, possibly also changing the
 * bead types in MoleculeType[].Bead array for new ones.
 */
void ChangeMolecules(SYSTEM *Sys_orig, SYSTEM Sys_add, bool beads, bool name) {
  COUNT *Count_orig = &Sys_orig->Count, *Count_add = &Sys_add.Count,
        count_old = *Count_orig;
  // replace bead types from Sys_orig with those from Sys_add if required //{{{
  if (beads) {
    // append bead types from Sys_add to Sys_orig
    Count_orig->BeadType += Count_add->BeadType;
    Sys_orig->BeadType = realloc(
        Sys_orig->BeadType, sizeof *Sys_orig->BeadType * Count_orig->BeadType);
    for (int i = 0; i < Count_add->BeadType; i++) {
      int new = i + count_old.BeadType;
      Sys_orig->BeadType[new] = Sys_add.BeadType[i];
      Sys_orig->BeadType[new].Number = 0;
    }
    /*
     * For molecules with the same number of beads and for those beads, rewrite:
     * 1) MoleculeType[].Bead array in Sys_orig with bead types from Sys_add
     * 2) Bead[].Type in Sys_orig with the same new types
     * 3) adjust numbers of beads in BeadType[].Number in Sys_orig
     */
    for (int i = 0; i < Count_orig->MoleculeType; i++) {
      int mtype_add = FindMoleculeType(*Sys_orig, Sys_orig->MoleculeType[i],
                                       Sys_add, 1, name);
      if (mtype_add != -1) {
        MOLECULETYPE *mt_orig = &Sys_orig->MoleculeType[i],
                     *mt_add = &Sys_add.MoleculeType[mtype_add];
        for (int j = 0; j < mt_orig->nBeads; j++) {
          int bt_orig = mt_orig->Bead[j];
          mt_orig->Bead[j] = mt_add->Bead[j] + count_old.BeadType; // 1)
          // go through molecules of type i
          for (int k = 0; k < mt_orig->Number; k++) {
            int mol_id = mt_orig->Index[k];
            MOLECULE *mol = &Sys_orig->Molecule[mol_id];
            int bead = mol->Bead[j];
            Sys_orig->Bead[bead].Type = mt_orig->Bead[j];  // 2)
            Sys_orig->BeadType[mt_orig->Bead[j]].Number++; // 3)
            Sys_orig->BeadType[bt_orig].Number--;          //
          }
        }
      }
    }
    // remake BeadType[].Index arrays as number of bead types changed
    for (int i = 0; i < count_old.BeadType; i++) {
      free(Sys_orig->BeadType[i].Index);
    }
    FillBeadTypeIndex(Sys_orig);
  } //}}}
  // add bond/angle/dihedral/improper types from Sys_add to Sys_orig //{{{
  if (Count_add->BondType > 0) {
    Count_orig->BondType += Count_add->BondType;
    Sys_orig->BondType = realloc(
        Sys_orig->BondType, sizeof *Sys_orig->BondType * Count_orig->BondType);
    memcpy(Sys_orig->BondType + count_old.BondType, Sys_add.BondType,
           sizeof *Sys_orig->BondType * Count_add->BondType);
  }
  if (Count_add->AngleType > 0) {
    Count_orig->AngleType += Count_add->AngleType;
    Sys_orig->AngleType =
        realloc(Sys_orig->AngleType,
                sizeof *Sys_orig->AngleType * Count_orig->AngleType);
    memcpy(Sys_orig->AngleType + count_old.AngleType, Sys_add.AngleType,
           sizeof *Sys_orig->AngleType * Count_add->AngleType);
  }
  if (Count_add->DihedralType > 0) {
    Count_orig->DihedralType += Count_add->DihedralType;
    Sys_orig->DihedralType =
        realloc(Sys_orig->DihedralType,
                sizeof *Sys_orig->DihedralType * Count_orig->DihedralType);
    memcpy(Sys_orig->DihedralType + count_old.DihedralType,
           Sys_add.DihedralType,
           sizeof *Sys_orig->DihedralType * Count_add->DihedralType);
  }
  if (Count_add->ImproperType > 0) {
    Count_orig->ImproperType += Count_add->ImproperType;
    Sys_orig->ImproperType =
        realloc(Sys_orig->ImproperType,
                sizeof *Sys_orig->ImproperType * Count_orig->ImproperType);
    memcpy(Sys_orig->ImproperType + count_old.ImproperType,
           Sys_add.ImproperType,
           sizeof *Sys_orig->ImproperType * Count_add->ImproperType);
  }                                                    //}}}
  for (int i = 0; i < Count_orig->MoleculeType; i++) { //{{{
    MOLECULETYPE *mt_orig = &Sys_orig->MoleculeType[i];
    int type = FindMoleculeType(*Sys_orig, Sys_orig->MoleculeType[i], Sys_add,
                                2, name);
    if (type != -1) {
      MOLECULETYPE *mt_add = &Sys_add.MoleculeType[type];
      // add bonds, if there are none in the original molecule type... //{{{
      if (mt_add->nBonds > 0 && mt_orig->nBonds == 0) {
        mt_orig->nBonds = mt_add->nBonds;
        mt_orig->Bond = malloc(sizeof *mt_orig->Bond * mt_orig->nBonds);
        memcpy(mt_orig->Bond, mt_add->Bond,
               sizeof *mt_add->Bond * mt_add->nBonds); //}}}
        // ...or just add bond types where missing //{{{
      } else if (Count_add->BondType > 0) {
        for (int j = 0; j < mt_orig->nBonds; j++) {
          for (int k = 0; k < mt_add->nBonds; k++) {
            if (mt_orig->Bond[j][0] == mt_add->Bond[k][0] &&
                mt_orig->Bond[j][1] == mt_add->Bond[k][1] &&
                mt_orig->Bond[j][2] == -1 && mt_add->Bond[k][2] != -1) {
              mt_orig->Bond[j][2] = mt_add->Bond[k][2] + count_old.BondType;
              break;
            }
          }
        }
      } //}}}
      // add angles, if there are none in the original molecule type... //{{{
      if (mt_add->nAngles > 0 && mt_orig->nAngles == 0) {
        mt_orig->nAngles = mt_add->nAngles;
        mt_orig->Angle = malloc(sizeof *mt_orig->Angle * mt_orig->nAngles);
        memcpy(mt_orig->Angle, mt_add->Angle,
               sizeof *mt_add->Angle * mt_add->nAngles); //}}}
        // ...or just add angle types where missing //{{{
      } else if (Count_add->AngleType > 0) {
        for (int j = 0; j < mt_orig->nAngles; j++) {
          for (int k = 0; k < mt_orig->nAngles; k++) {
            if (mt_orig->Angle[j][0] == mt_add->Angle[k][0] &&
                mt_orig->Angle[j][1] == mt_add->Angle[k][1] &&
                mt_orig->Angle[j][2] == mt_add->Angle[k][2] &&
                mt_orig->Angle[j][3] == -1 && mt_add->Angle[k][3] != -1) {
              mt_orig->Angle[j][3] = mt_add->Angle[k][3] + count_old.AngleType;
              break;
            }
          }
        }
      } //}}}
      // add dihedrals, if there are none in the original molecule type... //{{{
      if (mt_add->nDihedrals > 0 && mt_orig->nDihedrals == 0) {
        mt_orig->nDihedrals = mt_add->nDihedrals;
        mt_orig->Dihedral =
            malloc(sizeof *mt_orig->Dihedral * mt_orig->nDihedrals);
        memcpy(mt_orig->Dihedral, mt_add->Dihedral,
               sizeof *mt_add->Dihedral * mt_add->nDihedrals); //}}}
        // ...or just add dihedral types where missing //{{{
      } else if (Count_add->DihedralType > 0) {
        for (int j = 0; j < mt_orig->nDihedrals; j++) {
          for (int k = 0; k < mt_orig->nDihedrals; k++) {
            if (mt_orig->Dihedral[j][0] == mt_add->Dihedral[k][0] &&
                mt_orig->Dihedral[j][1] == mt_add->Dihedral[k][1] &&
                mt_orig->Dihedral[j][2] == mt_add->Dihedral[k][2] &&
                mt_orig->Dihedral[j][3] == mt_add->Dihedral[k][3] &&
                mt_orig->Dihedral[j][4] == -1 && mt_add->Dihedral[j][4] != -1) {
              mt_orig->Dihedral[j][4] =
                  mt_add->Dihedral[j][4] + count_old.DihedralType;
              break;
            }
          }
        }
      } //}}}
      // add impropers, if there are none in the original molecule type... //{{{
      if (mt_add->nImpropers > 0 && mt_orig->nImpropers == 0) {
        mt_orig->nImpropers = mt_add->nImpropers;
        mt_orig->Improper =
            malloc(sizeof *mt_orig->Improper * mt_orig->nImpropers);
        memcpy(mt_orig->Improper, mt_add->Improper,
               sizeof *mt_add->Improper * mt_add->nImpropers); //}}}
        // ...or just add improper types where missing //{{{
      } else if (Count_add->ImproperType > 0) {
        for (int j = 0; j < mt_orig->nImpropers; j++) {
          for (int k = 0; k < mt_orig->nImpropers; k++) {
            if (mt_orig->Improper[j][0] == mt_add->Improper[k][0] &&
                mt_orig->Improper[j][1] == mt_add->Improper[k][1] &&
                mt_orig->Improper[j][2] == mt_add->Improper[k][2] &&
                mt_orig->Improper[j][3] == mt_add->Improper[k][3] &&
                mt_orig->Improper[j][4] == -1 && mt_add->Improper[j][4] != -1) {
              mt_orig->Improper[j][4] =
                  mt_add->Improper[j][4] + count_old.ImproperType;
              break;
            }
          }
        }
      } //}}}
    }
  } //}}}
  // make sure all stuff is properly counted and there's nothing extra
  CountBondAngleDihedralImproper(Sys_orig);
  PruneSystem(Sys_orig);
} //}}}
// test whether two bead types are identical //{{{
bool SameBeadType(BEADTYPE bt_1, BEADTYPE bt_2) {
  if (strcmp(bt_1.Name, bt_2.Name) == 0 && bt_1.Charge == bt_2.Charge &&
      bt_1.Mass == bt_2.Mass && bt_1.Radius == bt_2.Radius) {
    return true;
  } else {
    return false;
  }
} //}}}
// create new bead/molecule type, realloc'ing the appropriate array
// NewBeadType() //{{{
void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char name[],
                 double charge, double mass, double radius) {
  int btype = *number_of_types;
  (*number_of_types)++;
  *BeadType = realloc(*BeadType, sizeof **BeadType * (btype + 1));
  strcpy((*BeadType)[btype].Name, name);
  (*BeadType)[btype].Number = 0;
  (*BeadType)[btype].Charge = charge;
  (*BeadType)[btype].Mass = mass;
  (*BeadType)[btype].Radius = radius;
}; //}}}
// NewMolType() //{{{
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char *name,
                int n_beads, int n_bonds, int n_angles, int n_dihedrals,
                int n_impropers) {
  int mtype = (*n_types)++;
  *MoleculeType = realloc(*MoleculeType, sizeof **MoleculeType * (*n_types));
  // copy new name to MoleculeType[].Name
  strncpy((*MoleculeType)[mtype].Name, name, MOL_NAME);
  // initialize struct members
  (*MoleculeType)[mtype].Number = 1;
  (*MoleculeType)[mtype].nBeads = n_beads;
  (*MoleculeType)[mtype].Bead =
      calloc(n_beads, sizeof *(*MoleculeType)[mtype].Bead);
  (*MoleculeType)[mtype].nBonds = n_bonds;
  if (n_bonds > 0) {
    (*MoleculeType)[mtype].Bond =
        calloc(n_bonds, sizeof *(*MoleculeType)[mtype].Bond);
  }
  (*MoleculeType)[mtype].nAngles = n_angles;
  if (n_angles > 0) {
    (*MoleculeType)[mtype].Angle =
        calloc(n_angles, sizeof *(*MoleculeType)[mtype].Angle);
  }
  (*MoleculeType)[mtype].nDihedrals = n_dihedrals;
  if (n_dihedrals > 0) {
    (*MoleculeType)[mtype].Dihedral =
        calloc(n_dihedrals, sizeof *(*MoleculeType)[mtype].Dihedral);
  }
  (*MoleculeType)[mtype].nImpropers = n_impropers;
  if (n_impropers > 0) {
    (*MoleculeType)[mtype].Improper =
        calloc(n_impropers, sizeof *(*MoleculeType)[mtype].Improper);
  }
  (*MoleculeType)[mtype].nBTypes = 0;
}; //}}}
// identify bead type based on name //{{{
int FindBeadType(char name[], SYSTEM System) {
  int type;
  for (int i = 0; i < System.Count.BeadType; i++) {
    if (strcmp(name, System.BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in BeadType struct
  return -1;
} //}}}
// identify molecule type based on name //{{{
int FindMoleculeName(char name[], SYSTEM System) {
  int type = -1;
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    if (strcmp(name, System.MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return type;
} //}}}
// identify molecule type based on name only or on other parameters too //{{{
/*
 * If name==true, then find molecule type according to mode:
 *   0 - check nothing else
 *   1 - check number of beads
 *   2 - check number and types of beads
 *   3 - check everything
 * name = true/false for checking/ignoring molecule name
 */
int FindMoleculeType(SYSTEM Sys1, MOLECULETYPE mt, SYSTEM Sys2, int mode,
                     bool name) {
  // just to be sure the function's mode parameter is correct //{{{
  if (mode < 0 || mode > 3) {
    strcpy(ERROR_MSG, "FindMoleculeType() - mode parameter must be <0,3>\n");
    PrintError();
    exit(1);
  } //}}}
  MOLECULETYPE *mt_1 = &mt;
  // find the same name
  for (int i = 0; i < Sys2.Count.MoleculeType; i++) {
    MOLECULETYPE *mt_2 = &Sys2.MoleculeType[i];
    if (strcmp(mt_1->Name, mt_2->Name) == 0 || !name) {
      if (mode == 0) { // only name checked
        return i;
      }
      // check number of beads //{{{
      if (mt_1->nBeads != mt_2->nBeads) {
        goto end_loop; // not 'continue;' to avoid break/continue combo later
      } else if (mode == 1) {
        return i;
      } //}}}
      // check bead order //{{{
      for (int j = 0; j < mt_1->nBeads; j++) {
        int bt_1 = mt_1->Bead[j], bt_2 = mt_2->Bead[j];
        if (!SameBeadType(Sys1.BeadType[bt_1], Sys2.BeadType[bt_2])) {
          goto end_loop;
        }
      }
      if (mode == 2) {
        return i;
      } //}}}
      // check bonds //{{{
      if (mt_1->nBonds != mt_2->nBonds) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nBonds; j++) {
        if (!SameArray(mt_1->Bond[j], mt_2->Bond[j], 3)) {
          goto end_loop;
        }
      } //}}}
      // check angles //{{{
      if (mt_1->nAngles != mt_2->nAngles) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nAngles; j++) {
        if (!SameArray(mt_1->Angle[j], mt_2->Angle[j], 4)) {
          goto end_loop;
        }
      } //}}}
      // check dihedrals //{{{
      if (mt_1->nDihedrals != mt_2->nDihedrals) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nDihedrals; j++) {
        if (!SameArray(mt_1->Dihedral[j], mt_2->Dihedral[j], 5)) {
          goto end_loop;
        }
      } //}}}
      // check impropers //{{{
      if (mt_1->nImpropers != mt_2->nImpropers) {
        goto end_loop;
      }
      for (int j = 0; j < mt_1->nImpropers; j++) {
        if (!SameArray(mt_1->Improper[j], mt_2->Improper[j], 5)) {
          goto end_loop;
        }
      }         //}}}
      return i; // assumes mode=3, obviously
    }
  end_loop:;
  }
  return -1;
} //}}}
// cleanse System by removing molecule/bead types with .Number=0, etc. //{{{
void PruneSystem(SYSTEM *System) {
  SYSTEM S_old = CopySystem(*System);
  FreeSystem(System);
  InitSystem(System);
  COUNT *Count = &System->Count, *Count_old = &S_old.Count;
  System->Box = S_old.Box;
  *Count = *Count_old; // some counts will change later
  // allocate memory for bead count arrays
  System->Bead = realloc(System->Bead, sizeof(BEAD) * Count->Bead);
  System->BeadCoor = realloc(System->BeadCoor,
                             sizeof *System->BeadCoor * Count->Bead);
  System->Bonded = realloc(System->Bonded,
                           sizeof *System->Bonded * Count->Bonded);
  System->BondedCoor = realloc(System->BondedCoor,
                               sizeof *System->BondedCoor * Count->Bonded);
  System->Unbonded = realloc(System->Unbonded,
                             sizeof *System->Unbonded * Count->Unbonded);
  System->UnbondedCoor = realloc(System->UnbondedCoor,
                                 sizeof *System->UnbondedCoor *
                                 Count->Unbonded);
  // copy Bead/Unbonded/Bonded arrays & create new BeadType array //{{{
  int count_unbonded = 0, count_bonded = 0, count_all = 0;
  // arrays for mapping old bead ids/types to new ones
  int *b_id_old_to_new = calloc(Count_old->Bead, sizeof *b_id_old_to_new),
      *bt_old_to_new = calloc(Count_old->BeadType, sizeof *bt_old_to_new);
  Count->BeadType = 0;
  for (int i = 0; i < Count_old->Bead; i++) {
    if (S_old.Bead[i].InTimestep) {
      // create new bead type if it doesn't exist yet in the pruned system
      int old_type = S_old.Bead[i].Type, new_type = -1;
      BEADTYPE *bt_old = &S_old.BeadType[old_type];
      for (int j = 0; j < Count->BeadType; j++) {
        if (SameBeadType(*bt_old, System->BeadType[j])) {
          new_type = j;
          break;
        }
      }
      if (new_type == -1) {
        new_type = Count->BeadType;
        NewBeadType(&System->BeadType, &Count->BeadType, bt_old->Name,
                    bt_old->Charge, bt_old->Mass, bt_old->Radius);
      }

      System->Bead[count_all] = S_old.Bead[i];

      if (System->Bead[count_all].Molecule == -1) {
        System->Unbonded[count_unbonded] = count_all;
        System->UnbondedCoor[count_unbonded] = count_all;
        count_unbonded++;
      } else {
        System->Bonded[count_bonded] = count_all;
        System->BondedCoor[count_bonded] = count_all;
        count_bonded++;
      }
      System->BeadCoor[count_all] = count_all;

      System->Bead[count_all].Type = new_type;
      b_id_old_to_new[i] = count_all;

      System->BeadType[new_type].Number++;
      bt_old_to_new[old_type] = new_type;

      count_all++;
    }
  }
  RenameBeadTypes(System);
  Count->Bead = count_all;
  Count->BeadCoor = Count->Bead;
  Count->Bonded = count_bonded;
  Count->BondedCoor = Count->Bonded;
  Count->Unbonded = count_unbonded;
  Count->UnbondedCoor = Count->Unbonded; //}}}
  // copy Molecule array & create a new MoleculeType array //{{{
  Count->MoleculeType = 0;
  Count->Molecule = 0;
  for (int i = 0; i < Count_old->Molecule; i++) {
    MOLECULE *mol_old = &S_old.Molecule[i];
    MOLECULETYPE *mt_old = &S_old.MoleculeType[mol_old->Type];
    int c_bead = 0;
    for (int j = 0; j < mt_old->nBeads; j++) {
      int id = mol_old->Bead[j];
      if (S_old.Bead[id].InTimestep) {
        c_bead++;
      }
    }
    if (c_bead > 0) { // should the molecule be in the pruned system?
      // create new type for mt_old as some beads may be missing
      MOLECULETYPE mt_old_new;
      InitMoleculeType(&mt_old_new);
      strcpy(mt_old_new.Name, mt_old->Name);
      mt_old_new.Number = 1;
      mt_old_new.nBeads = c_bead;
      mt_old_new.Bead = malloc(sizeof *mt_old_new.Bead * mt_old_new.nBeads);
      c_bead = 0;
      // TODO: array to renumber beads ...huh? It's b_old_to_old_new, isn't it?
      int b_old_to_old_new[mt_old->nBeads];
      for (int j = 0; j < mt_old->nBeads; j++) {
        b_old_to_old_new[j] = -1;
        int id = mol_old->Bead[j];
        if (S_old.Bead[id].InTimestep) {
          mt_old_new.Bead[c_bead] = S_old.Bead[id].Type;
          b_old_to_old_new[j] = c_bead;
          c_bead++;
        }
      }
      // copy bonds to the mt_old_new molecule type
      int count = 0;
      if (mt_old->nBonds > 0) {
        mt_old_new.Bond = malloc(sizeof *mt_old_new.Bond * mt_old->nBonds);
        for (int j = 0; j < mt_old->nBonds; j++) {
          int id1 = mt_old->Bond[j][0], id2 = mt_old->Bond[j][1];
          id1 = b_old_to_old_new[id1];
          id2 = b_old_to_old_new[id2];
          if (id1 != -1 && id2 != -1) {
            mt_old_new.Bond[count][0] = id1;
            mt_old_new.Bond[count][1] = id2;
            mt_old_new.Bond[count][2] = mt_old->Bond[j][2];
            count++;
          }
        }
        if (count == 0) {
          free(mt_old_new.Bond);
        }
      }
      mt_old_new.nBonds = count;
      // copy angles to the mt_old_new molecule type
      count = 0;
      if (mt_old->nAngles > 0) {
        mt_old_new.Angle = malloc(sizeof *mt_old_new.Angle * mt_old->nAngles);
        for (int j = 0; j < mt_old->nAngles; j++) {
          int id1 = mt_old->Angle[j][0], id2 = mt_old->Angle[j][1],
              id3 = mt_old->Angle[j][2];
          id1 = b_old_to_old_new[id1];
          id2 = b_old_to_old_new[id2];
          id3 = b_old_to_old_new[id3];
          if (id1 != -1 && id2 != -1) {
            mt_old_new.Angle[count][0] = id1;
            mt_old_new.Angle[count][1] = id2;
            mt_old_new.Angle[count][2] = id3;
            mt_old_new.Angle[count][3] = mt_old->Angle[j][3];
            count++;
          }
        }
        if (count == 0) {
          free(mt_old_new.Angle);
        }
      }
      mt_old_new.nAngles = count;
      // copy dihedrals to the mt_old_new molecule type
      count = 0;
      if (mt_old->nDihedrals > 0) {
        mt_old_new.Dihedral =
            malloc(sizeof *mt_old_new.Dihedral * mt_old->nDihedrals);
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int id1 = mt_old->Dihedral[j][0], id2 = mt_old->Dihedral[j][1],
              id3 = mt_old->Dihedral[j][2], id4 = mt_old->Dihedral[j][3];
          id1 = b_old_to_old_new[id1];
          id2 = b_old_to_old_new[id2];
          id3 = b_old_to_old_new[id3];
          id4 = b_old_to_old_new[id4];
          if (id1 != -1 && id2 != -1) {
            mt_old_new.Dihedral[count][0] = id1;
            mt_old_new.Dihedral[count][1] = id2;
            mt_old_new.Dihedral[count][2] = id3;
            mt_old_new.Dihedral[count][3] = id4;
            mt_old_new.Dihedral[count][4] = mt_old->Dihedral[j][4];
            count++;
          }
        }
        if (count == 0) {
          free(mt_old_new.Dihedral);
        }
      }
      mt_old_new.nDihedrals = count;
      // copy impropers to the mt_old_new molecule type
      count = 0;
      if (mt_old->nImpropers > 0) {
        mt_old_new.Improper =
            malloc(sizeof *mt_old_new.Improper * mt_old->nImpropers);
        for (int j = 0; j < mt_old->nImpropers; j++) {
          int id1 = mt_old->Improper[j][0], id2 = mt_old->Improper[j][1],
              id3 = mt_old->Improper[j][2], id4 = mt_old->Improper[j][3];
          id1 = b_old_to_old_new[id1];
          id2 = b_old_to_old_new[id2];
          id3 = b_old_to_old_new[id3];
          id4 = b_old_to_old_new[id4];
          if (id1 != -1 && id2 != -1) {
            mt_old_new.Improper[count][0] = id1;
            mt_old_new.Improper[count][1] = id2;
            mt_old_new.Improper[count][2] = id3;
            mt_old_new.Improper[count][3] = id4;
            mt_old_new.Improper[count][4] = mt_old->Improper[j][4];
            count++;
          }
        }
        if (count == 0) {
          free(mt_old_new.Improper);
        }
      }
      mt_old_new.nImpropers = count;

      int new_id = Count->Molecule;
      Count->Molecule++;
      System->Molecule =
          realloc(System->Molecule, sizeof(MOLECULE) * Count->Molecule);
      MOLECULE *mol_new = &System->Molecule[new_id];
      *mol_new = *mol_old;
      mol_new->Bead = calloc(c_bead, sizeof *mol_new->Bead);
      c_bead = 0;
      for (int j = 0; j < mt_old->nBeads; j++) {
        int id = mol_old->Bead[j];
        if (S_old.Bead[id].InTimestep) {
          mol_new->Bead[c_bead] = b_id_old_to_new[id];
          System->Bead[b_id_old_to_new[id]].Molecule = new_id;
          c_bead++;
        }
      }

      /*
       * Is the molecule type already in the pruned system (check based on all
       * molecule type information)?
       */
      int new_type = FindMoleculeType(S_old, mt_old_new, *System, 3, true);
      FreeMoleculeTypeEssentials(&mt_old_new);
      if (new_type != -1) { // yes, the molecule type is in the pruned system
        mol_new->Type = new_type;
        System->MoleculeType[new_type].Number++;
      } else { // no, it isn't; create a new one
        int new_new_type = Count->MoleculeType, c_bond = 0, c_angle = 0,
            c_dihedral = 0, c_improper = 0;
        // count bonds in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nBonds; j++) {
          int *id = BondIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep) {
            c_bond++;
          }
        } //}}}
        // count angles in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nAngles; j++) {
          int *id = AngleIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            c_angle++;
          }
        } //}}}
        // count dihedrals in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep && S_old.Bead[id[3]].InTimestep) {
            c_dihedral++;
          }
        } //}}}
        // count impropers in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nImpropers; j++) {
          int *id = ImproperIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep && S_old.Bead[id[3]].InTimestep) {
            c_improper++;
          }
        } //}}}
        NewMolType(&System->MoleculeType, &Count->MoleculeType, mt_old->Name,
                   c_bead, c_bond, c_angle, c_dihedral, c_improper);
        System->Molecule[new_id].Type = new_new_type;
        MOLECULETYPE *mt_new = &System->MoleculeType[new_new_type];
        // copy beads to the new molecule type //{{{
        // map internal MoleculeType[].Bead ids to new ones (some may disappear)
        int *id_old_to_new = calloc(mt_old->nBeads, sizeof *id_old_to_new);
        for (int j = 0; j < mt_old->nBeads; j++) {
          id_old_to_new[j] = -1;
        }
        c_bead = 0;
        for (int j = 0; j < mt_old->nBeads; j++) {
          int id = mol_old->Bead[j], old_btype = mt_old->Bead[j];
          if (S_old.Bead[id].InTimestep) {
            mt_new->Bead[c_bead] = bt_old_to_new[old_btype];
            id_old_to_new[j] = c_bead;
            c_bead++;
          }
        } //}}}
        // copy bonds to the new molecule type //{{{
        c_bond = 0;
        for (int j = 0; j < mt_old->nBonds; j++) {
          int *id = BondIndices(S_old, i, j), last = mt_old->Bond[j][2];
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep) {
            mt_new->Bond[c_bond][0] = id_old_to_new[id[2]];
            mt_new->Bond[c_bond][1] = id_old_to_new[id[3]];
            mt_new->Bond[c_bond][2] = last;
            c_bond++;
          }
        } //}}}
        // copy angles to the new molecule type //{{{
        c_angle = 0;
        for (int j = 0; j < mt_old->nAngles; j++) {
          int *id = AngleIndices(S_old, i, j), last = mt_old->Angle[j][3];
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            mt_new->Angle[c_angle][0] = id_old_to_new[id[3]];
            mt_new->Angle[c_angle][1] = id_old_to_new[id[4]];
            mt_new->Angle[c_angle][2] = id_old_to_new[id[5]];
            mt_new->Angle[c_angle][3] = last;
            c_angle++;
          }
        } //}}}
        // copy dihedrals to the new molecule type //{{{
        c_dihedral = 0;
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j), last = mt_old->Dihedral[j][4];
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep && S_old.Bead[id[3]].InTimestep) {
            mt_new->Dihedral[c_dihedral][0] = id_old_to_new[id[4]];
            mt_new->Dihedral[c_dihedral][1] = id_old_to_new[id[5]];
            mt_new->Dihedral[c_dihedral][2] = id_old_to_new[id[6]];
            mt_new->Dihedral[c_dihedral][3] = id_old_to_new[id[7]];
            mt_new->Dihedral[c_dihedral][4] = last;
            c_dihedral++;
          }
        } //}}}
        // copy impropers to the new molecule type //{{{
        c_improper = 0;
        for (int j = 0; j < mt_old->nImpropers; j++) {
          int *id = ImproperIndices(S_old, i, j), last = mt_old->Improper[j][4];
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep && S_old.Bead[id[3]].InTimestep) {
            mt_new->Improper[c_improper][0] = id_old_to_new[id[4]];
            mt_new->Improper[c_improper][1] = id_old_to_new[id[5]];
            mt_new->Improper[c_improper][2] = id_old_to_new[id[6]];
            mt_new->Improper[c_improper][3] = id_old_to_new[id[7]];
            mt_new->Improper[c_improper][4] = last;
            c_improper++;
          }
        } //}}}
        free(id_old_to_new);
      }
    }
  } //}}}
  RenameBeadTypes(System);
  RenameMoleculeTypes(System);
  FillBeadTypeIndex(System);
  FillMoleculeTypeIndex(System);
  // rewrite Index_mol //{{{
  System->Index_mol = realloc(System->Index_mol, sizeof *System->Index_mol *
                                                     (Count->HighestResid + 1));
  for (int i = 0; i <= Count->HighestResid; i++) {
    System->Index_mol[i] = -1;
  }
  for (int i = 0; i < Count->Molecule; i++) {
    System->Index_mol[System->Molecule[i].Index] = i;
  } //}}}
  // copy bond/angle/dihedral/improper types //{{{
  // prune bond types //{{{
  if (Count_old->BondType > 0) {
    Count->BondType = 0;
    int *type_old_to_new = calloc(Count_old->BondType, sizeof *type_old_to_new);
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nBonds; j++) {
        int old_tbond = mt_i->Bond[j][2];
        if (old_tbond != -1) {
          bool new = true;
          PARAMS *tbond = &S_old.BondType[old_tbond];
          for (int k = 0; k < Count->BondType; k++) {
            if (tbond->a == System->BondType[k].a &&
                tbond->b == System->BondType[k].b &&
                tbond->c == System->BondType[k].c) {
              mt_i->Bond[j][2] = k;
              new = false;
              break;
            }
          }
          if (new) {
            int type = Count->BondType;
            Count->BondType++;
            System->BondType = realloc(
                System->BondType, Count->BondType * sizeof *System->BondType);
            System->BondType[type] = S_old.BondType[old_tbond];
            type_old_to_new[old_tbond] = type;
          }
        }
      }
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nBonds; j++) {
        int old_tbond = mt_i->Bond[j][2];
        if (old_tbond != -1) {
          mt_i->Bond[j][2] = type_old_to_new[old_tbond];
        } else {
          mt_i->Bond[j][2] = -1;
        }
      }
    }
    free(type_old_to_new);
  } //}}}
  // prune angle types //{{{
  if (Count->AngleType > 0) {
    Count->AngleType = 0;
    int *type_old_to_new =
        calloc(Count_old->AngleType, sizeof *type_old_to_new);
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nAngles; j++) {
        int old_tangle = mt_i->Angle[j][3];
        bool new = true;
        PARAMS *tangle = &S_old.AngleType[old_tangle];
        for (int k = 0; k < Count->AngleType; k++) {
          if (tangle->a == System->AngleType[k].a &&
              tangle->b == System->AngleType[k].b &&
              tangle->c == System->AngleType[k].c) {
            mt_i->Angle[j][3] = k;
            new = false;
            break;
          }
        }
        if (new) {
          int type = Count->AngleType;
          Count->AngleType++;
          System->AngleType = realloc(
              System->AngleType, Count->AngleType * sizeof *System->AngleType);
          System->AngleType[type] = S_old.AngleType[old_tangle];
          type_old_to_new[old_tangle] = type;
        }
      }
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nAngles; j++) {
        int old_tangle = mt_i->Angle[j][3];
        mt_i->Angle[j][3] = type_old_to_new[old_tangle];
      }
    }
    free(type_old_to_new);
  } //}}}
  // prune dihedral types //{{{
  if (Count->DihedralType > 0) {
    Count->DihedralType = 0;
    int *type_old_to_new =
        calloc(Count_old->DihedralType, sizeof *type_old_to_new);
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        int old_tdihed = mt_i->Dihedral[j][4];
        bool new = true;
        PARAMS *tdihed = &S_old.DihedralType[old_tdihed];
        for (int k = 0; k < Count->DihedralType; k++) {
          if (tdihed->a == System->DihedralType[k].a &&
              tdihed->b == System->DihedralType[k].b &&
              tdihed->c == System->DihedralType[k].c) {
            mt_i->Dihedral[j][4] = k;
            new = false;
            break;
          }
        }
        if (new) {
          int type = Count->DihedralType;
          Count->DihedralType++;
          System->DihedralType =
              realloc(System->DihedralType,
                      Count->DihedralType * sizeof *System->DihedralType);
          System->DihedralType[type] = S_old.DihedralType[old_tdihed];
          type_old_to_new[old_tdihed] = type;
        }
      }
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        int old_tdihed = mt_i->Dihedral[j][4];
        mt_i->Dihedral[j][4] = type_old_to_new[old_tdihed];
      }
    }
    free(type_old_to_new);
  } //}}}
  // prune improper types //{{{
  if (Count->ImproperType > 0) {
    Count->ImproperType = 0;
    int *type_old_to_new =
        calloc(Count_old->ImproperType, sizeof *type_old_to_new);
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nImpropers; j++) {
        int old_timpro = mt_i->Improper[j][4];
        bool new = true;
        PARAMS *timpro = &S_old.ImproperType[old_timpro];
        for (int k = 0; k < Count->ImproperType; k++) {
          if (timpro->a == System->ImproperType[k].a &&
              timpro->b == System->ImproperType[k].b &&
              timpro->c == System->ImproperType[k].c) {
            mt_i->Improper[j][4] = k;
            new = false;
            break;
          }
        }
        if (new) {
          int type = Count->ImproperType;
          Count->ImproperType++;
          System->ImproperType =
              realloc(System->ImproperType,
                      Count->ImproperType * sizeof *System->ImproperType);
          System->ImproperType[type] = S_old.ImproperType[old_timpro];
          type_old_to_new[old_timpro] = type;
        }
      }
    }
    for (int i = 0; i < Count->MoleculeType; i++) {
      MOLECULETYPE *mt_i = &System->MoleculeType[i];
      for (int j = 0; j < mt_i->nImpropers; j++) {
        int old_timpro = mt_i->Improper[j][4];
        mt_i->Improper[j][4] = type_old_to_new[old_timpro];
      }
    }
    free(type_old_to_new);
  } //}}}
  //}}}
  CountBondAngleDihedralImproper(System);
  for (int i = 0; i < Count->MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  FreeSystem(&S_old);
  free(b_id_old_to_new);
  free(bt_old_to_new);
} //}}}
// copy molecule type //{{{
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old) { //{{{
  MOLECULETYPE mt_new = CopyMoleculeTypeEssentials(mt_old);
  // MoleculeType[].Index array
  mt_new.Index = malloc(sizeof *mt_new.Index * mt_new.Number);
  memcpy(mt_new.Index, mt_old.Index, sizeof *mt_old.Index * mt_old.Number);
  // MoleculeType[].BType array
  if (mt_new.nBTypes > 0) {
    mt_new.BType = malloc(sizeof *mt_new.BType * mt_new.nBTypes);
    memcpy(mt_new.BType, mt_old.BType, sizeof *mt_old.BType * mt_old.nBTypes);
  }
  return mt_new;
} //}}}
MOLECULETYPE CopyMoleculeTypeEssentials(MOLECULETYPE mt_old) { //{{{
  MOLECULETYPE mt_new;
  mt_new = mt_old;
  // MoleculeType[].Bead array
  if (mt_new.nBeads > 0) {
    mt_new.Bead = malloc(sizeof *mt_new.Bead * mt_new.nBeads);
    memcpy(mt_new.Bead, mt_old.Bead, sizeof *mt_old.Bead * mt_old.nBeads);
  } else {
    // should never happen
    snprintf(ERROR_MSG, LINE,
             "molecule type without beads (%s%s%s); should "
             "never happen!",
             ErrYellow(), mt_new.Name, ErrYellow());
    PrintWarning();
  }
  // MoleculeType[].Bond array
  if (mt_new.nBonds > 0) {
    mt_new.Bond = malloc(sizeof *mt_new.Bond * mt_new.nBonds);
    memcpy(mt_new.Bond, mt_old.Bond, sizeof *mt_old.Bond * mt_old.nBonds);
  }
  // MoleculeType[].Angle array
  if (mt_new.nAngles > 0) {
    mt_new.Angle = malloc(sizeof *mt_new.Angle * mt_new.nAngles);
    memcpy(mt_new.Angle, mt_old.Angle, sizeof *mt_old.Angle * mt_old.nAngles);
  }
  // MoleculeType[].Dihedral array
  if (mt_new.nDihedrals > 0) {
    mt_new.Dihedral = malloc(sizeof *mt_new.Dihedral * mt_new.nDihedrals);
    memcpy(mt_new.Dihedral, mt_old.Dihedral,
           sizeof *mt_old.Dihedral * mt_old.nDihedrals);
  }
  // MoleculeType[].Improper array
  if (mt_new.nImpropers > 0) {
    mt_new.Improper = malloc(sizeof *mt_new.Improper * mt_new.nImpropers);
    memcpy(mt_new.Improper, mt_old.Improper,
           sizeof *mt_old.Improper * mt_old.nImpropers);
  }
  return mt_new;
} //}}}
  //}}}
// concatenateSystems() //{{{
/*
 * Assumes S_out needs reallocating memory to accommodate S_in.
 */
// TODO: some warning about S_in being empty
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box) {
  COUNT Count_old = S_out->Count; // copy the original COUNT
  COUNT *Count_out = &S_out->Count, *Count_in = &S_in.Count;
  S_out->Box = Box;
  // BeadType //{{{
  if (Count_in->BeadType > 0) {
    Count_out->BeadType += Count_in->BeadType;
    S_out->BeadType =
        realloc(S_out->BeadType, sizeof *S_out->BeadType * Count_out->BeadType);
    //  memcpy(S_out->BeadType + Count_old.BeadType, S_in.BeadType,
    //         sizeof (BEADTYPE) * Count_in->BeadType);
    for (int i = 0; i < Count_in->BeadType; i++) {
      int new = i + Count_old.BeadType;
      S_out->BeadType[new] = S_in.BeadType[i];
      S_out->BeadType[new].Index = malloc(sizeof *S_out->BeadType[new].Index *
                                          S_out->BeadType[new].Number);
      for (int j = 0; j < S_out->BeadType[new].Number; j++) {
        S_out->BeadType[new].Index[j] =
            S_in.BeadType[i].Index[j] + Count_old.Bead;
      }
    }
  } else {
    strcpy(ERROR_MSG, "no bead types to add to the system");
    PrintWarning();
    return;
  } //}}}
  // Bead //{{{
  if (Count_in->Bead > 0) {
    Count_out->Bead += Count_in->Bead;
    S_out->Bead = realloc(S_out->Bead, sizeof(BEAD) * Count_out->Bead);
    for (int i = 0; i < Count_in->Bead; i++) {
      int new = i + Count_old.Bead;
      S_out->Bead[new] = S_in.Bead[i];
      S_out->Bead[new].Type += Count_old.BeadType;
      if (S_out->Bead[new].Molecule != -1) {
        S_out->Bead[new].Molecule += Count_old.Molecule;
      }
    }
  } else {
    strcpy(ERROR_MSG, "no beads to add to the system");
    PrintWarning();
    return;
  } //}}}
  // Bonded //{{{
  if (Count_in->Bonded > 0) {
    Count_out->Bonded += Count_in->Bonded;
    S_out->Bonded =
        realloc(S_out->Bonded, sizeof *S_out->Bonded * Count_out->Bonded);
    for (int i = 0; i < Count_in->Bonded; i++) {
      int new = i + Count_old.Bonded;
      S_out->Bonded[new] = S_in.Bonded[i] + Count_old.Bead;
    }
  } //}}}
  // BondedCoor //{{{
  if (Count_in->BondedCoor > 0) {
    Count_out->BondedCoor += Count_in->BondedCoor;
    S_out->BondedCoor = realloc(S_out->BondedCoor,
                                sizeof *S_out->BondedCoor * Count_out->Bonded);
    for (int i = 0; i < Count_in->BondedCoor; i++) {
      int new = i + Count_old.BondedCoor;
      S_out->BondedCoor[new] = S_in.BondedCoor[i] + Count_old.Bead;
    }
  } //}}}
  // Unonded //{{{
  if (Count_in->Unbonded > 0) {
    Count_out->Unbonded += Count_in->Unbonded;
    S_out->Unbonded =
        realloc(S_out->Unbonded, sizeof *S_out->Unbonded * Count_out->Unbonded);
    for (int i = 0; i < Count_in->Unbonded; i++) {
      int new = i + Count_old.Unbonded;
      S_out->Unbonded[new] = S_in.Unbonded[i] + Count_old.Bead;
    }
  } //}}}
  // UnondedCoor //{{{
  if (Count_in->UnbondedCoor > 0) {
    Count_out->UnbondedCoor += Count_in->UnbondedCoor;
    S_out->UnbondedCoor = realloc(
        S_out->UnbondedCoor, sizeof *S_out->UnbondedCoor * Count_out->Unbonded);
    for (int i = 0; i < Count_in->UnbondedCoor; i++) {
      int new = i + Count_old.UnbondedCoor;
      S_out->UnbondedCoor[new] = S_in.UnbondedCoor[i] + Count_old.Bead;
    }
  } //}}}
  // BeadCoor //{{{
  if (Count_in->BeadCoor > 0) {
    Count_out->BeadCoor += Count_in->BeadCoor;
    S_out->BeadCoor =
        realloc(S_out->BeadCoor, sizeof *S_out->BeadCoor * Count_out->Bead);
    for (int i = 0; i < Count_in->BeadCoor; i++) {
      int new = i + Count_old.BeadCoor;
      S_out->BeadCoor[new] = S_in.BeadCoor[i] + Count_old.Bead;
    }
  } //}}}
  // MoleculeType //{{{
  if (Count_in->MoleculeType > 0) {
    Count_out->MoleculeType += Count_in->MoleculeType;
    S_out->MoleculeType = realloc(
        S_out->MoleculeType, sizeof(MOLECULETYPE) * Count_out->MoleculeType);
    for (int i = 0; i < Count_in->MoleculeType; i++) {
      int new = i + Count_old.MoleculeType;
      MOLECULETYPE *mt_out = &S_out->MoleculeType[new],
                   *mt_in = &S_in.MoleculeType[i];
      *mt_out = *mt_in;
      // MoleculeType[].Bead
      mt_out->Bead = calloc(mt_out->nBeads, sizeof *mt_out->Bead);
      for (int j = 0; j < mt_out->nBeads; j++) {
        mt_out->Bead[j] = mt_in->Bead[j] + Count_old.BeadType;
      }
      // MoleculeType[].Bond
      if (mt_out->nBonds > 0) {
        mt_out->Bond = calloc(mt_out->nBonds, sizeof *mt_out->Bond);
        memcpy(mt_out->Bond, mt_in->Bond, sizeof *mt_in->Bond * mt_in->nBonds);
        for (int j = 0; j < mt_out->nBonds; j++) {
          mt_out->Bond[j][2] += Count_old.BondType;
        }
      }
      // MoleculeType[].Angle
      if (mt_out->nAngles > 0) {
        mt_out->Angle = calloc(mt_out->nAngles, sizeof *mt_out->Angle);
        memcpy(mt_out->Angle, mt_in->Angle,
               sizeof *mt_in->Angle * mt_in->nAngles);
        for (int j = 0; j < mt_out->nAngles; j++) {
          mt_out->Angle[j][3] += Count_old.AngleType;
        }
      }
      // MoleculeType[].Dihedral
      if (mt_out->nDihedrals > 0) {
        mt_out->Dihedral = calloc(mt_out->nDihedrals, sizeof *mt_out->Dihedral);
        memcpy(mt_out->Dihedral, mt_in->Dihedral,
               sizeof *mt_in->Dihedral * mt_in->nDihedrals);
        for (int j = 0; j < mt_out->nDihedrals; j++) {
          mt_out->Dihedral[j][4] += Count_old.DihedralType;
        }
      }
      // MoleculeType[].Improper
      if (mt_out->nImpropers > 0) {
        mt_out->Improper = calloc(mt_out->nImpropers, sizeof *mt_out->Improper);
        memcpy(mt_out->Improper, mt_in->Improper,
               sizeof *mt_in->Improper * mt_in->nImpropers);
        for (int j = 0; j < mt_out->nImpropers; j++) {
          mt_out->Improper[j][4] += Count_old.ImproperType;
        }
      }
      // MoleculeType[].BType
      mt_out->BType = calloc(mt_out->nBTypes, sizeof *mt_out->BType);
      for (int j = 0; j < mt_out->nBTypes; j++) {
        mt_out->BType[j] = mt_in->BType[j] + Count_old.BeadType;
      }
      // MoleculeType[].Index
      mt_out->Index = malloc(sizeof *mt_out->Index * mt_out->Number);
      for (int j = 0; j < mt_out->Number; j++) {
        mt_out->Index[j] = mt_in->Index[j] + Count_old.Molecule;
      }
    }
  } //}}}
  // Molecule & Index_mol //{{{
  if (Count_in->Molecule > 0) {
    Count_out->Molecule += Count_in->Molecule;
    S_out->Molecule =
        realloc(S_out->Molecule, sizeof(MOLECULE) * Count_out->Molecule);
    for (int i = 0; i < Count_in->Molecule; i++) {
      MOLECULE *mol_out = &S_out->Molecule[i + Count_old.Molecule],
               *mol_in = &S_in.Molecule[i];
      int type = mol_in->Type + Count_old.MoleculeType;
      *mol_out = *mol_in;
      mol_out->Type = type;
      mol_out->Index += Count_old.HighestResid;
      mol_out->Bead =
          malloc(sizeof *mol_out->Bead * S_out->MoleculeType[type].nBeads);
      for (int j = 0; j < S_out->MoleculeType[type].nBeads; j++) {
        mol_out->Bead[j] = mol_in->Bead[j] + Count_old.Bead;
      }
    }
    Count_out->HighestResid += Count_in->HighestResid;
    S_out->Index_mol =
        realloc(S_out->Index_mol,
                sizeof *S_out->Index_mol * (Count_out->HighestResid + 1));
    for (int i = 0; i <= Count_out->HighestResid; i++) {
      S_out->Index_mol[i] = -1;
    }
    for (int i = 0; i < Count_out->Molecule; i++) {
      S_out->Index_mol[S_out->Molecule[i].Index] = i;
    }
  } //}}}
  // BondType //{{{
  if (Count_in->BondType > 0) {
    Count_out->BondType += Count_in->BondType;
    S_out->BondType =
        realloc(S_out->BondType, sizeof *S_out->BondType * Count_out->BondType);
    memcpy(S_out->BondType + Count_old.BondType, S_in.BondType,
           sizeof *S_out->BondType * Count_out->BondType);
  } //}}}
  // AngleType //{{{
  if (Count_in->AngleType > 0) {
    Count_out->AngleType += Count_in->AngleType;
    S_out->AngleType = realloc(S_out->AngleType,
                               sizeof *S_out->AngleType * Count_out->AngleType);
    memcpy(S_out->AngleType + Count_old.AngleType, S_in.AngleType,
           sizeof *S_out->AngleType * Count_out->AngleType);
  } //}}}
  // DihedralType //{{{
  if (Count_in->DihedralType > 0) {
    Count_out->DihedralType += Count_in->DihedralType;
    S_out->DihedralType =
        realloc(S_out->DihedralType,
                sizeof *S_out->DihedralType * Count_out->DihedralType);
    memcpy(S_out->DihedralType + Count_old.DihedralType, S_in.DihedralType,
           sizeof *S_out->DihedralType * Count_out->DihedralType);
  } //}}}
  // ImproperType //{{{
  if (Count_in->ImproperType > 0) {
    Count_out->ImproperType += Count_in->ImproperType;
    S_out->ImproperType =
        realloc(S_out->ImproperType,
                sizeof *S_out->ImproperType * Count_out->ImproperType);
    memcpy(S_out->ImproperType + Count_old.ImproperType, S_in.ImproperType,
           sizeof *S_out->ImproperType * Count_out->ImproperType);
  } //}}}
} //}}}
// TODO: split CheckSystem to CheckCount, CheckBeadType, etc.
// check that the System struct doesn't contain an error //{{{
void CheckSystem(SYSTEM System, char file[]) {
  COUNT *Count = &System.Count;
  if (Count->Molecule > 0 &&
      Count->Molecule > (Count->HighestResid + 1)) { //{{{
    strcpy(ERROR_MSG, "Count.HighestResid is lower than Count.Molecules");
    PrintErrorFile(file, "\0", "\0");
    fprintf(stderr, "%s, Count.HighestResid = %s%d", ErrRed(), ErrYellow(),
            Count->HighestResid);
    fprintf(stderr, "%s, Count.Molecule = %s%d%s\n", ErrRed(), ErrYellow(),
            Count->Molecule, ErrColourReset());
  } //}}}
  // total number of beads //{{{
  // i) just unbonded+bonded
  int count = Count->Unbonded + Count->Bonded;
  if (count != Count->Bead) {
    strcpy(ERROR_MSG, "unbonded and bonded beads do not add up properly!");
    PrintErrorFile(file, "\0", "\0");
    fprintf(stderr, "%s, unbonded: %s%d%s", ErrRed(), ErrYellow(),
            Count->Unbonded, ErrRed());
    fprintf(stderr, ", bonded: %s%d%s", ErrYellow(), Count->Bonded, ErrRed());
    fprintf(stderr, ", sum should be: %s%d%s\n", ErrYellow(), Count->Bead,
            ErrColourReset());
  }
  // ii) from BeadType
  count = 0;
  for (int i = 0; i < Count->BeadType; i++) {
    count += System.BeadType[i].Number;
  }
  if (count != Count->Bead) {
    strcpy(ERROR_MSG, "number of beads in bead types do not add up properly!");
    PrintErrorFile(file, "\0", "\0");
    fprintf(stderr, "%s, sum should be: %s%d%s", ErrRed(), ErrYellow(),
            Count->Bead, ErrRed());
    fprintf(stderr, " but is: %s%d%s\n", ErrYellow(), count, ErrRed());
    fprintf(stderr, "Number | Name%s\n", ErrYellow());
    for (int i = 0; i < Count->BeadType; i++) {
      fprintf(stderr, "%6d %s|%s %s\n", System.BeadType[i].Number, ErrRed(),
              ErrYellow(), System.BeadType[i].Name);
    }
    fputs(ErrColourReset(), stderr);
  } //}}}
  // BeadCoor array //{{{
  for (int i = 0; i < Count->BeadCoor; i++) {
    if (System.BeadCoor[i] < 0 || System.BeadCoor[i] >= Count->Bead) {
      strcpy(ERROR_MSG, "incorrect index in BeadCoor array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, BeadCoor[%s%d%s] = %s%d%s\n", ErrRed(), ErrYellow(),
              i, ErrRed(), ErrYellow(), System.BeadCoor[i], ErrColourReset());
      break;
    }
  } //}}}
  // Bonded array //{{{
  for (int i = 0; i < Count->Bonded; i++) {
    if (System.Bonded[i] < 0 || System.Bonded[i] >= Count->Bead) {
      strcpy(ERROR_MSG, "incorrect index in Bonded array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Bonded[%s%d%s] = %s%d%s\n", ErrRed(), ErrYellow(), i,
              ErrRed(), ErrYellow(), System.Bonded[i], ErrColourReset());
      break;
    }
  } //}}}
  // BondedCoor array //{{{
  for (int i = 0; i < Count->BondedCoor; i++) {
    if (System.BondedCoor[i] < 0 || System.BondedCoor[i] >= Count->Bead) {
      strcpy(ERROR_MSG, "incorrect index in BondedCoor array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, BondedCoor[%s%d%s] = %s%d%s\n", ErrRed(),
              ErrYellow(), i, ErrRed(), ErrYellow(), System.BondedCoor[i],
              ErrColourReset());
      break;
    }
  } //}}}
  // Unbonded array //{{{
  for (int i = 0; i < Count->Unbonded; i++) {
    if (System.Unbonded[i] < 0 || System.Unbonded[i] >= Count->Bead) {
      strcpy(ERROR_MSG, "incorrect index in Unbonded array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Unbonded[%s%d%s] = %s%d%s\n", ErrRed(), ErrYellow(),
              i, ErrRed(), ErrYellow(), System.Unbonded[i], ErrColourReset());
      break;
    }
  } //}}}
  // UnbondedCoor array //{{{
  for (int i = 0; i < Count->UnbondedCoor; i++) {
    if (System.UnbondedCoor[i] < 0 || System.UnbondedCoor[i] >= Count->Bead) {
      strcpy(ERROR_MSG, "incorrect index in UnbondedCoor array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, UnbondedCoor[%s%d%s] = %s%d%s\n", ErrRed(),
              ErrYellow(), i, ErrRed(), ErrYellow(), System.UnbondedCoor[i],
              ErrColourReset());
      break;
    }
  } //}}}
  // Index_mol array //{{{
  for (int i = 0; i < Count->HighestResid; i++) {
    if (System.Index_mol[i] < -1 || System.Index_mol[i] >= Count->Molecule) {
      strcpy(ERROR_MSG, "incorrect index in Index_mol array");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Index_mol[%s%d%s] = %s%d%s\n", ErrRed(), ErrYellow(),
              i, ErrRed(), ErrYellow(), System.Index_mol[i], ErrColourReset());
      break;
    }
  } //}}}
  // BeadType[].Index array //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System.BeadType[i];
    for (int j = 0; j < bt_i->Number; j++) {
      if (bt_i->Index[j] < 0 || bt_i->Index[j] >= Count->Bead) {
        strcpy(ERROR_MSG, "incorrect index in BeadType[].Index array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, BeadType[%s%d%s].Index[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), bt_i->Index[j], ErrColourReset());
        break;
      }
    }
  } //}}}
  // MoleculeType[]. arrays //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    // Index array //{{{
    for (int j = 0; j < mt_i->Number; j++) {
      if (mt_i->Index[j] < 0 || mt_i->Index[j] >= Count->Molecule) {
        strcpy(ERROR_MSG, "incorrect index in MoleculeType[].Index array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Index[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt_i->Index[j], ErrColourReset());
        break;
      }
    } //}}}
    // Bead array //{{{
    for (int j = 0; j < mt_i->nBeads; j++) {
      if (mt_i->Bead[j] < 0 || mt_i->Bead[j] >= Count->BeadType) {
        strcpy(ERROR_MSG, "incorrect index in Bead array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bead[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt_i->Bead[j], ErrColourReset());
        break;
      }
    } //}}}
    // Bond array //{{{
    for (int j = 0; j < mt_i->nBonds; j++) {
      if (mt_i->Bond[j][0] < 0 || mt_i->Bond[j][0] >= mt_i->nBeads ||
          mt_i->Bond[j][1] < 0 || mt_i->Bond[j][1] >= mt_i->nBeads ||
          mt_i->Bond[j][0] == mt_i->Bond[j][1]) {
        strcpy(ERROR_MSG, "incorrect index in Bond array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bond[%s%d%s][0..1]", ErrRed(),
                ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed());
        fprintf(stderr, " = %s%d %d%s\n", ErrYellow(), mt_i->Bond[j][0],
                mt_i->Bond[j][1], ErrColourReset());
        break;
      }
      if (mt_i->Bond[j][2] < -1 || mt_i->Bond[j][2] >= Count->BondType) {
        strcpy(ERROR_MSG, "incorrect bond type in Bond array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bond[2] = %s%d%s\n", ErrRed(),
                ErrYellow(), i, ErrRed(), ErrYellow(), mt_i->Bond[j][2],
                ErrColourReset());
        break;
      }
    } //}}}
    // Angle array //{{{
    for (int j = 0; j < mt_i->nAngles; j++) {
      if (mt_i->Angle[j][0] < 0 || mt_i->Angle[j][0] >= mt_i->nBeads ||
          mt_i->Angle[j][1] < 0 || mt_i->Angle[j][1] >= mt_i->nBeads ||
          mt_i->Angle[j][2] < 0 || mt_i->Angle[j][2] >= mt_i->nBeads ||
          mt_i->Angle[j][0] == mt_i->Angle[j][1] ||
          mt_i->Angle[j][0] == mt_i->Angle[j][2] ||
          mt_i->Angle[j][1] == mt_i->Angle[j][2]) {
        strcpy(ERROR_MSG, "incorrect index in Angle array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Angle[0..2] = %s%d %d %d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt_i->Angle[j][0], mt_i->Angle[j][1], mt_i->Angle[j][2],
                ErrColourReset());
        break;
      }
      if (mt_i->Angle[j][3] < -1 || mt_i->Angle[j][3] >= Count->AngleType) {
        strcpy(ERROR_MSG, "incorrect angle type in Angle array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Angle[3] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt_i->Angle[j][3], ErrColourReset());
        break;
      }
    } //}}}
    // Dihedral array //{{{
    for (int j = 0; j < mt_i->nDihedrals; j++) {
      if (mt_i->Dihedral[j][0] < 0 || mt_i->Dihedral[j][0] >= mt_i->nBeads ||
          mt_i->Dihedral[j][1] < 0 || mt_i->Dihedral[j][1] >= mt_i->nBeads ||
          mt_i->Dihedral[j][2] < 0 || mt_i->Dihedral[j][2] >= mt_i->nBeads ||
          mt_i->Dihedral[j][3] < 0 || mt_i->Dihedral[j][3] >= mt_i->nBeads ||
          mt_i->Dihedral[j][0] == mt_i->Dihedral[j][1] ||
          mt_i->Dihedral[j][0] == mt_i->Dihedral[j][2] ||
          mt_i->Dihedral[j][0] == mt_i->Dihedral[j][3] ||
          mt_i->Dihedral[j][1] == mt_i->Dihedral[j][2] ||
          mt_i->Dihedral[j][1] == mt_i->Dihedral[j][3] ||
          mt_i->Dihedral[j][2] == mt_i->Dihedral[j][3]) {
        strcpy(ERROR_MSG, "incorrect index in Dihedral array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Dihedral[0..3] = %s",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow());
        fprintf(stderr, "%d %d %d %d%s\n", mt_i->Dihedral[j][0],
                mt_i->Dihedral[j][1], mt_i->Dihedral[j][2],
                mt_i->Dihedral[j][3], ErrColourReset());
        break;
      }
      if (mt_i->Dihedral[j][4] < -1 ||
          mt_i->Dihedral[j][4] >= Count->DihedralType) {
        strcpy(ERROR_MSG, "incorrect dihedral type in Dihedral array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Dihedral[4] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt_i->Dihedral[j][4], ErrColourReset());
        break;
      }
    } //}}}
    // Improper array //{{{
    for (int j = 0; j < mt_i->nImpropers; j++) {
      if (mt_i->Improper[j][0] < 0 || mt_i->Improper[j][0] >= mt_i->nBeads ||
          mt_i->Improper[j][1] < 0 || mt_i->Improper[j][1] >= mt_i->nBeads ||
          mt_i->Improper[j][2] < 0 || mt_i->Improper[j][2] >= mt_i->nBeads ||
          mt_i->Improper[j][3] < 0 || mt_i->Improper[j][3] >= mt_i->nBeads ||
          mt_i->Improper[j][0] == mt_i->Improper[j][1] ||
          mt_i->Improper[j][0] == mt_i->Improper[j][2] ||
          mt_i->Improper[j][0] == mt_i->Improper[j][3] ||
          mt_i->Improper[j][1] == mt_i->Improper[j][2] ||
          mt_i->Improper[j][1] == mt_i->Improper[j][3] ||
          mt_i->Improper[j][2] == mt_i->Improper[j][3]) {
        strcpy(ERROR_MSG, "incorrect index in Improper array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Improper[0..3] = %s",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow());
        fprintf(stderr, "%d %d %d %d%s\n", mt_i->Improper[j][0],
                mt_i->Improper[j][1], mt_i->Improper[j][2],
                mt_i->Improper[j][3], ErrColourReset());
        break;
      }
      if (mt_i->Improper[j][4] < -1 ||
          mt_i->Improper[j][4] >= Count->ImproperType) {
        strcpy(ERROR_MSG, "incorrect improper type in Improper array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Improper[4] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt_i->Improper[j][4], ErrColourReset());
        break;
      }
    } //}}}
    // BType array //{{{
    for (int j = 0; j < mt_i->nBTypes; j++) {
      if (mt_i->BType[j] < 0 || mt_i->BType[j] >= Count->BeadType) {
        strcpy(ERROR_MSG, "incorrect index in BType array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].BType[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt_i->BType[j], ErrColourReset());
        break;
      }
    } //}}}
  }   //}}}
  // bead types & Bead[].Molecule //{{{
  for (int i = 0; i < Count->Bead; i++) {
    int type = System.Bead[i].Type;
    if (type < 0 || type >= Count->BeadType) {
      strcpy(ERROR_MSG, "incorrect bead type for a bead");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Bead[%s%d%s].Type = %s%d%s\n", ErrRed(), ErrYellow(),
              i, ErrRed(), ErrYellow(), type, ErrColourReset());
      break;
    }
    int mol = System.Bead[i].Molecule;
    if (mol < -1 || mol >= Count->Molecule) {
      strcpy(ERROR_MSG, "incorrect molecule index for a bead");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Bead[%s%d%s].Molecule = %s%d%s\n", ErrRed(),
              ErrYellow(), i, ErrRed(), ErrYellow(), mol, ErrColourReset());
      break;
    }
  } //}}}
  // molecule types & Molecule[].Bead & Index arrays //{{{
  int *test = calloc(Count->HighestResid + 1, sizeof *test);
  for (int i = 0; i <= Count->HighestResid; i++) {
    test[i] = -1;
  }
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol_i = &System.Molecule[i];
    int type = mol_i->Type;
    if (type < 0 || type >= Count->MoleculeType) {
      strcpy(ERROR_MSG, "incorrect molecule type for a molecule");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, Molecule[%s%d%s].Type = %s%d%s\n", ErrRed(),
              ErrYellow(), i, ErrRed(), ErrYellow(), type, ErrColourReset());
      break;
    }
    for (int j = 0; j < System.MoleculeType[type].nBeads; j++) {
      if (mol_i->Bead[j] < 0 || mol_i->Bead[j] >= Count->Bead) {
        strcpy(ERROR_MSG, "incorrect index in Bead array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, Molecule[%s%d%s].Bead[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mol_i->Bead[j], ErrColourReset());
        break;
      }
    }
    if (test[mol_i->Index] > -1) {
      strcpy(ERROR_MSG, "same molecule index with multiple molecules");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, index %s%d%s with molecuces %s%d %d%s\n", ErrRed(),
              ErrYellow(), mol_i->Index, ErrRed(), ErrYellow(),
              test[mol_i->Index], i, ErrColourReset());
      break;
    }
    test[mol_i->Index] = i;
  }
  free(test); //}}}
} //}}}

// Helper functions for manipulating coordinates
// wrap coordinates into simulation box and/or join molecules //{{{
void WrapJoinCoordinates(SYSTEM *System, bool wrap, bool join) {
  if (System->Box.Volume != -1 && (wrap || join)) {
    // transform coordinates into fractional ones for non-orthogonal box
    ToFractionalCoor(System);
    if (wrap) { // wrap coordinates into the simulation box
      RestorePBC(System);
    }
    if (join) { // join molecules by removing periodic boundary conditions
      RemovePBCMolecules(System);
    }
    // transform back to 'normal' coordinates for non-orthogonal box
    FromFractionalCoor(System);
  }
} //}}}
// distance between two beads; in the range <-BoxLength/2,BoxLength/2) //{{{
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength) {
  // distance vector
  VECTOR rij;
  rij.x = id1.x - id2.x;
  rij.y = id1.y - id2.y;
  rij.z = id1.z - id2.z;
  // remove periodic boundary conditions in x-direction
  while (rij.x >= (BoxLength.x / 2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x < -(BoxLength.x / 2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y >= (BoxLength.y / 2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y < -(BoxLength.y / 2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z >= (BoxLength.z / 2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z < -(BoxLength.z / 2))
    rij.z = rij.z + BoxLength.z;
  return rij;
} //}}}

// identify input coordinate and structure files //{{{
bool InputCoorStruct(int argc, char *argv[], char coor_file[], int *coor_type,
                     char struct_file[], int *struct_type) {
  int ext;
  char extension[6][EXTENSION];
  // input structure file (-i option) //{{{
  if (FileOption(argc, argv, "-i", struct_file, LINE)) {
    exit(1);
  }
  if (struct_file[0] != '\0') { // -i option is present
    if (strcasecmp(struct_file, "FIELD") == 0) {
      *struct_type = FIELD_FILE;
    } else {
      ext = 6;
      strcpy(extension[0], ".vsf");
      strcpy(extension[1], ".vtf");
      strcpy(extension[2], ".xyz");
      strcpy(extension[3], ".lammpstrj");
      strcpy(extension[4], ".data");
      strcpy(extension[5], ".field");
      ext = ErrorExtension(struct_file, ext, extension);
      switch (ext) {
      case 0:
        *struct_type = VSF_FILE;
        break;
      case 1:
        *struct_type = VSF_FILE;
        break;
      case 2:
        *struct_type = XYZ_FILE;
        break;
      case 3:
        *struct_type = LTRJ_FILE;
        break;
      case 4:
        *struct_type = LDATA_FILE;
        break;
      case 5:
        *struct_type = FIELD_FILE;
        break;
      default: // something wrong; should never happen
        return false;
      }
    }
  } //}}}
  // input coordinate file //{{{
  ext = 4;
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  strcpy(extension[2], ".xyz");
  strcpy(extension[3], ".lammpstrj");
  ext = ErrorExtension(coor_file, ext, extension);
  // define coordinate type and possibly vtf structure file
  switch (ext) {
  case 0:
    *coor_type = VCF_FILE;
    // copy to 'structure' with vsf ending, if 'structure' is empty
    if (struct_file[0] == '\0') {
      int last = -1;
      for (int i = 0; i < strlen(coor_file); i++) {
        if (coor_file[i] == '.') {
          last = i;
        }
      }
      strncpy(struct_file, coor_file, last);
      strcat(struct_file, ".vsf");
      *struct_type = VSF_FILE;
    }
    break;
  case 1: // vtf full file
    *coor_type = VCF_FILE;
    // use also as a 'structure' file, if 'structure' is empty
    if (struct_file[0] == '\0') {
      strcpy(struct_file, coor_file);
      *struct_type = VSF_FILE;
    }
    break;
  case 2: // xyz
    *coor_type = XYZ_FILE;
    // use also as a 'structure' file, if 'structure' is empty
    if (struct_file[0] == '\0') {
      strcpy(struct_file, coor_file);
      *struct_type = XYZ_FILE;
    }
    break;
  case 3: // lammpstrj
    *coor_type = LTRJ_FILE;
    // use also as a 'structure' file, if 'structure' is empty
    if (struct_file[0] == '\0') {
      strcpy(struct_file, coor_file);
      *struct_type = LTRJ_FILE;
    }
    break;
  default: // something wrong; should never happen
    return false;
  } //}}}
  return true;
} //}}}

// create a cell-linked list //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead, int **Head,
                int **Link, double cell_size, INTVECTOR *n_cells, int *Dcx,
                int *Dcy, int *Dcz) {
  n_cells->x = ceil(BoxLength.x / cell_size),
  n_cells->y = ceil(BoxLength.y / cell_size),
  n_cells->z = ceil(BoxLength.z / cell_size);
  // allocate arrays
  *Head = malloc(sizeof **Head * n_cells->x * n_cells->y * n_cells->z);
  *Link = malloc(sizeof **Link * Counts.BeadsCoor);
  for (int i = 0; i < (n_cells->x * n_cells->y * n_cells->z); i++) {
    (*Head)[i] = -1;
  }
  // sort beads into cells
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int cell = (int)(Bead[i].Position.x / cell_size) +
               (int)(Bead[i].Position.y / cell_size) * n_cells->x +
               (int)(Bead[i].Position.z / cell_size) * n_cells->x * n_cells->y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  }
  // coordinates of adjoining cells
  int x[14] = {0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};
  int y[14] = {0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
  int z[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 14; i++) {
    Dcx[i] = x[i];
    Dcy[i] = y[i];
    Dcz[i] = z[i];
  }
} //}}}

// verbose output (print various structures and some such)
void VerboseOutput(SYSTEM System) { //{{{
  PrintCount(System.Count);
  PrintBeadType(System);
  PrintBondType(System);
  PrintAngleType(System);
  PrintDihedralType(System);
  PrintImproperType(System);
  PrintMoleculeType(System);
  if (System.Box.Volume != -1) {
    putchar('\n');
    PrintBox(System.Box);
  }
  putchar('\n');
} //}}}
void PrintCount(COUNT Count) { //{{{
  bool coor = false;
  if (Count.Bead != Count.BeadCoor && Count.BeadCoor > 0) {
    coor = true;
  }
  fprintf(stdout, "\nCounts of\n");
  fprintf(stdout, "  Bead Types:     %d\n", Count.BeadType);
  fprintf(stdout, "  All Beads:      %d\n", Count.Bead);
  if (coor && Count.Bead > 0) {
    fprintf(stdout, "    In Coor File: %d\n", Count.BeadCoor);
  }
  fprintf(stdout, "  Bonded Beads:   %d\n", Count.Bonded);
  if (coor && Count.Bonded > 0) {
    fprintf(stdout, "    In Coor File: %d\n", Count.BondedCoor);
  }
  fprintf(stdout, "  Unbonded Beads: %d\n", Count.Unbonded);
  if (coor && Count.Unbonded > 0) {
    fprintf(stdout, "    In Coor File: %d\n", Count.UnbondedCoor);
  }
  fprintf(stdout, "  Molecule Types: %d\n", Count.MoleculeType);
  fprintf(stdout, "  Molecules:      %d", Count.Molecule);
  if (Count.BondType > 0) {
    fprintf(stdout, "\n  Bond Types:     %d", Count.BondType);
  }
  if (Count.Bonded > 0) {
    fprintf(stdout, "\n  Bonds:          %d", Count.Bond);
  }
  if (Count.AngleType > 0) {
    fprintf(stdout, "\n  Angle Types:    %d", Count.AngleType);
  }
  if (Count.Angle > 0) {
    fprintf(stdout, "\n  Angles:         %d", Count.Angle);
  }
  if (Count.DihedralType > 0) {
    fprintf(stdout, "\n  Dihedral Types: %d", Count.DihedralType);
  }
  if (Count.Dihedral > 0) {
    fprintf(stdout, "\n  Dihedrals:      %d", Count.Dihedral);
  }
  if (Count.ImproperType > 0) {
    fprintf(stdout, "\n  Improper Types: %d", Count.ImproperType);
  }
  if (Count.Improper > 0) {
    fprintf(stdout, "\n  Impropers:      %d", Count.Improper);
  }
  fprintf(stdout, "\n\n");
} //}}}
void PrintBeadType(SYSTEM System) { //{{{
  // some stuff to properly align the fields //{{{
  int precision = 3,     // number of decimal digits
      longest_name = 0,  // longest bead type name
      max_number = 0,    // maximum number of beads
      max_q = 0,         // maximum charge
      max_m = 0,         // maximum mass
      max_r = 0;         // maximum radius
  bool negative = false; // extra space for '-' if there's negative charge
  for (int i = 0; i < System.Count.BeadType; i++) {
    int length = strlen(System.BeadType[i].Name);
    if (length > longest_name) {
      longest_name = length;
    }
    if (System.BeadType[i].Number > max_number) {
      max_number = System.BeadType[i].Number;
    }
    if (System.BeadType[i].Charge < 0) {
      negative = true;
    }
    if (System.BeadType[i].Charge != CHARGE &&
        fabs(System.BeadType[i].Charge) > max_q) {
      max_q = floor(fabs(System.BeadType[i].Charge));
    }
    if (System.BeadType[i].Mass != MASS && System.BeadType[i].Mass > max_m) {
      max_m = floor(System.BeadType[i].Mass);
    }
    if (System.BeadType[i].Radius != RADIUS &&
        System.BeadType[i].Radius > max_r) {
      max_r = floor(System.BeadType[i].Radius);
    }
  }
  // number of digits of the highest_number
  if (max_number == 0) {
    max_number = 1;
  } else {
    max_number = floor(log10(max_number)) + 1;
  }
  // number of digits of the charge
  if (max_q == 0) {
    max_q = 1;
  } else {
    max_q = floor(log10(max_q)) + 1;
  }
  max_q += 1 + precision; // +1 for the decimal point
  if (negative) {
    max_q++; // extra space for minus sign
  }
  // number of digits of the mass
  if (max_m == 0) {
    max_m = 1;
  } else {
    max_m = floor(log10(max_m)) + 1 + precision + 1;
  }
  // number of digits of the radius
  if (max_r == 0) {
    max_r = 1;
  } else {
    max_r = floor(log10(max_m)) + 1 + precision + 1;
  }
  // number of digits of the number of types
  int types_digits = floor(log10(System.Count.BeadType)) + 1;
  //}}}
  // print the information
  for (int i = 0; i < System.Count.BeadType; i++) {
    fprintf(stdout, "BeadType[%*d] = {", types_digits, i);
    fprintf(stdout, ".Name = %*s, ", longest_name, System.BeadType[i].Name);
    fprintf(stdout, ".Number = %*d, ", max_number, System.BeadType[i].Number);
    // print charge
    fprintf(stdout, ".Charge = ");
    if (System.BeadType[i].Charge != CHARGE) {
      fprintf(stdout, "%*.*f, ", max_q, precision, System.BeadType[i].Charge);
    } else {
      for (int j = 0; j < (max_q - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    // print mass
    fprintf(stdout, ".Mass = ");
    if (System.BeadType[i].Mass != MASS) {
      fprintf(stdout, "%*.*f, ", max_m, precision, System.BeadType[i].Mass);
    } else {
      for (int j = 0; j < (max_m - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    fprintf(stdout, ".Radius = ");
    // print radius
    if (System.BeadType[i].Radius != RADIUS) {
      fprintf(stdout, "%*.*f", max_r, precision, System.BeadType[i].Radius);
    } else {
      for (int j = 0; j < (max_r - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a");
    }
    fprintf(stdout, "}\n");
  }
  putchar('\n');
} //}}}
void PrintMoleculeType(SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    fprintf(stdout, "MoleculeType[%d] = {\n", i);
    fprintf(stdout, "  .Name       = %s,\n", System.MoleculeType[i].Name);
    fprintf(stdout, "  .Number     = %d,\n", System.MoleculeType[i].Number);
    // print bead types (list all beads) //{{{
    fprintf(stdout, "  .nBeads     = %d,\n", System.MoleculeType[i].nBeads);
    fprintf(stdout, "  .Bead       = {");
    for (int j = 0; j < System.MoleculeType[i].nBeads; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      int type = System.MoleculeType[i].Bead[j];
      fprintf(stdout, "%s", System.BeadType[type].Name);
    }
    fprintf(stdout, "},\n"); //}}}
    // print bonds if there are any //{{{
    if (System.MoleculeType[i].nBonds > 0) {
      fprintf(stdout, "  .nBonds     = %d,\n", System.MoleculeType[i].nBonds);
      fprintf(stdout, "  .Bond       = {");
      for (int j = 0; j < System.MoleculeType[i].nBonds; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d", System.MoleculeType[i].Bond[j][0] + 1,
                System.MoleculeType[i].Bond[j][1] + 1);
        if (System.MoleculeType[i].Bond[j][2] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Bond[j][2] + 1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print angles if there are any //{{{
    if (System.MoleculeType[i].nAngles > 0) {
      fprintf(stdout, "  .nAngles    = %d,\n", System.MoleculeType[i].nAngles);
      fprintf(stdout, "  .Angle      = {");
      for (int j = 0; j < System.MoleculeType[i].nAngles; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d", System.MoleculeType[i].Angle[j][0] + 1,
                System.MoleculeType[i].Angle[j][1] + 1,
                System.MoleculeType[i].Angle[j][2] + 1);
        if (System.MoleculeType[i].Angle[j][3] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Angle[j][3] + 1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print dihedrals if there are any //{{{
    if (System.MoleculeType[i].nDihedrals > 0) {
      fprintf(stdout, "  .nDihedrals = %d,\n  .Dihedral   = {",
              System.MoleculeType[i].nDihedrals);
      for (int j = 0; j < System.MoleculeType[i].nDihedrals; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d-%d",
                System.MoleculeType[i].Dihedral[j][0] + 1,
                System.MoleculeType[i].Dihedral[j][1] + 1,
                System.MoleculeType[i].Dihedral[j][2] + 1,
                System.MoleculeType[i].Dihedral[j][3] + 1);
        if (System.MoleculeType[i].Dihedral[j][4] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Dihedral[j][4] + 1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print impropers if there are any //{{{
    if (System.MoleculeType[i].nImpropers > 0) {
      fprintf(stdout, "  .nImpropers = %d,\n  .Improper   = {",
              System.MoleculeType[i].nImpropers);
      for (int j = 0; j < System.MoleculeType[i].nImpropers; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d-%d",
                System.MoleculeType[i].Improper[j][0] + 1,
                System.MoleculeType[i].Improper[j][1] + 1,
                System.MoleculeType[i].Improper[j][2] + 1,
                System.MoleculeType[i].Improper[j][3] + 1);
        if (System.MoleculeType[i].Improper[j][4] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Improper[j][4] + 1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print bead types (just the which are present) //{{{
    fprintf(stdout, "  .nBTypes    = %d\n", System.MoleculeType[i].nBTypes);
    fprintf(stdout, "  .BType      = {");
    for (int j = 0; j < System.MoleculeType[i].nBTypes; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s",
              System.BeadType[System.MoleculeType[i].BType[j]].Name);
    } //}}}
    if (System.MoleculeType[i].Mass != MASS) {
      fprintf(stdout, "},\n  .Mass       = %.5f,\n",
              System.MoleculeType[i].Mass);
    } else {
      fprintf(stdout, "},\n  .Mass       = n/a,\n");
    }
    if (System.MoleculeType[i].Charge != CHARGE) {
      fprintf(stdout, "  .Charge     = %.5f\n}\n",
              System.MoleculeType[i].Charge);
    } else {
      fprintf(stdout, "  .Charge     = n/a\n}\n");
    }
  }
} //}}}
void PrintMolecule(SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.Molecule; i++) {
    int type = System.Molecule[i].Type;
    fprintf(stdout, "Molecule %3d (%d, %s):\n", i + 1, System.Molecule[i].Index,
            System.MoleculeType[type].Name);
    fprintf(stdout, " BEAD INDICES (%d): ", System.MoleculeType[type].nBeads);
    fputs("intramolecular; input file\n", stdout);
    for (int j = 0; j < System.MoleculeType[type].nBeads; j++) {
      int id = System.Molecule[i].Bead[j];
      fprintf(stdout, "   %3d; %5d\n", j + 1, id);
    }
  }
  fprintf(stdout, "\n");
} //}}}
void PrintBead(SYSTEM System) { //{{{
  fprintf(stdout, "Beads\n");
  fprintf(stdout, "<bead id>");
  fprintf(stdout, " (<bead type id>);");
  fprintf(stdout, " <molecule id>");
  fprintf(stdout, " (<molecule type id>);");
  fprintf(stdout, " <in coor>");
  putchar('\n');
  for (int i = 0; i < System.Count.Bead; i++) {
    BEAD *b = &System.Bead[i];
    fprintf(stdout, " %6d", i);
    fprintf(stdout, " (%3d);", b->Type);
    if (b->Molecule == -1) {
      fprintf(stdout, " %4s", "None");
      fprintf(stdout, "      ;");
    } else {
      fprintf(stdout, " %4d", System.Molecule[b->Molecule].Index);
      fprintf(stdout, " (%3d);", System.Molecule[b->Molecule].Type);
    }
    fprintf(stdout, " %s", b->InTimestep ? "yes" : " no");
    putchar('\n');
  }
} //}}}
void PrintBondType(SYSTEM System) { //{{{
  if (System.Count.BondType > 0) {
    fprintf(stdout, "Bond types\n");
    for (int i = 0; i < System.Count.BondType; i++) {
      fprintf(stdout, "   %lf %lf\n", System.BondType[i].a,
              System.BondType[i].b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintAngleType(SYSTEM System) { //{{{
  if (System.Count.AngleType > 0) {
    fprintf(stdout, "Angle types\n");
    for (int i = 0; i < System.Count.AngleType; i++) {
      fprintf(stdout, "   %lf %lf\n", System.AngleType[i].a,
              System.AngleType[i].b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintDihedralType(SYSTEM System) { //{{{
  if (System.Count.DihedralType > 0) {
    fprintf(stdout, "Dihedral types\n");
    for (int i = 0; i < System.Count.DihedralType; i++) {
      PARAMS *dihed = &System.DihedralType[i];
      fprintf(stdout, "   %lf %lf %lf\n", dihed->a, dihed->b, dihed->c);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintImproperType(SYSTEM System) { //{{{
  if (System.Count.ImproperType > 0) {
    fprintf(stdout, "Improper types\n");
    for (int i = 0; i < System.Count.ImproperType; i++) {
      PARAMS *imp = &System.ImproperType[i];
      fprintf(stdout, "   %lf %lf %lf\n", imp->a, imp->b, imp->c);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintBox(BOX Box) { //{{{
  fprintf(stdout, "Box = {\n");
  fprintf(stdout, "  .Length = (%lf, %lf, %lf)", Box.Length.x, Box.Length.y,
          Box.Length.z);
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    fprintf(stdout, ",  .alpha = %lf,", Box.alpha);
    fprintf(stdout, " .beta = %lf,", Box.beta);
    fprintf(stdout, " .gamma = %lf,\n", Box.gamma);
  }
  fprintf(stdout, "  .OrthoLength = (%lf, %lf, %lf),\n", Box.OrthoLength.x,
          Box.OrthoLength.y, Box.OrthoLength.z);
  fprintf(stdout, "  .alpha = %lf,\n", Box.alpha);
  fprintf(stdout, "  .beta  = %lf,\n", Box.beta);
  fprintf(stdout, "  .gamma = %lf,\n", Box.gamma);
  fprintf(stdout, "  .Volume = %lf,\n", Box.Volume);
  // print transform matrix //{{{
  for (int i = 0; i < 3; i++) {
    if (i == 0) {
      fprintf(stdout, "  .transform = (");
    } else {
      fprintf(stdout, "               (");
    }
    for (int j = 0; j < 3; j++) {
      if (Box.transform[i][j] >= 0) {
        putchar(' ');
      }
      fprintf(stdout, "%e", Box.transform[i][j]);
      if (j < 2) {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, ")\n");
  } //}}}
  // print inverse matrix //{{{
  for (int i = 0; i < 3; i++) {
    if (i == 0) {
      fprintf(stdout, "  .inverse = (");
    } else {
      fprintf(stdout, "             (");
    }
    for (int j = 0; j < 3; j++) {
      if (Box.inverse[i][j] >= 0) {
        putchar(' ');
      }
      fprintf(stdout, "%e", Box.inverse[i][j]);
      if (j < 2) {
        fprintf(stdout, ", ");
      }
    }
    fprintf(stdout, ")\n");
  } //}}}
  fprintf(stdout, "\n}\n");
} //}}}
void PrintByline(FILE *ptr, int argc, char *argv[]) { //{{{
  fprintf(ptr, "# Created by AnalysisTools v%s ", VERSION);
  fprintf(ptr, " (https://github.com/KaGaSi/AnalysisTools)\n");
  fprintf(ptr, "# command: ");
  PrintCommand(ptr, argc, argv);
} //}}}

// TODO: use Jacobi method
// calculate gyration tensor and various shape descriptors //{{{
VECTOR Gyration(int n, int *list, COUNTS Counts, BEADTYPE *BeadType,
                BEAD **Bead) {
  // gyration tensor (3x3 array)
  // use long double to ensure precision -- previous problem with truncation in
  // short chains
  struct Tensor {
    LONGVECTOR x, y, z;
  } GyrationTensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  VECTOR cog = GeomCentre(n, list, *Bead);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= cog.x;
    (*Bead)[list[i]].Position.y -= cog.y;
    (*Bead)[list[i]].Position.z -= cog.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    int id = list[i];
    GyrationTensor.x.x += Bead[id]->Position.x * Bead[id]->Position.x;
    GyrationTensor.x.y += Bead[id]->Position.x * Bead[id]->Position.y;
    GyrationTensor.x.z += Bead[id]->Position.x * Bead[id]->Position.z;
    GyrationTensor.y.y += Bead[id]->Position.y * Bead[id]->Position.y;
    GyrationTensor.y.z += Bead[id]->Position.y * Bead[id]->Position.z;
    GyrationTensor.z.z += Bead[id]->Position.z * Bead[id]->Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n; //}}}

  // characteristic polynomial:
  // a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0
  long double a_cube = -1;
  long double b_cube =
      GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
  long double c_cube = -GyrationTensor.x.x * GyrationTensor.y.y -
                       GyrationTensor.x.x * GyrationTensor.z.z -
                       GyrationTensor.y.y * GyrationTensor.z.z +
                       SQR(GyrationTensor.y.z) + SQR(GyrationTensor.x.y) +
                       SQR(GyrationTensor.x.z);
  long double d_cube =
      +GyrationTensor.x.x * GyrationTensor.y.y * GyrationTensor.z.z +
      2 * GyrationTensor.x.y * GyrationTensor.y.z * GyrationTensor.x.z -
      SQR(GyrationTensor.x.z) * GyrationTensor.y.y -
      SQR(GyrationTensor.x.y) * GyrationTensor.z.z -
      SQR(GyrationTensor.y.z) * GyrationTensor.x.x;

  // first root: either 0 or Newton's iterative method to get it //{{{
  long double root0 = 0;
  if (fabsl(d_cube) > 0.0000000001L) {
    // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
    long double a_deriv = 3 * a_cube;
    long double b_deriv = 2 * b_cube;
    long double c_deriv = c_cube;

    long double root1 = 1;

    while (fabsl(root0 - root1) > 0.0000000001L) {
      long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) +
                             c_cube * root0 + d_cube);
      long double f_deriv_root0 =
          (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
      root1 = root0 - f_root0 / f_deriv_root0;

      // swap root0 and root1 for the next iteration
      long double tmp = root0;
      root0 = root1;
      root1 = tmp;
    }
  } //}}}

  // find parameters of a quadratic equation a_quad*x^2+b_quad*x+c_quad=0 //{{{
  /*
   * derived by division:
   * (x^3+(b_cube/a_cube)*x^2+(c_cube/a_cube)*x+(d_cube/a_cube)):(x-root0)
   */
  long double a_quad = 1;
  long double b_quad = b_cube / a_cube + root0;
  long double c_quad =
      SQR(root0) + b_cube / a_cube * root0 + c_cube / a_cube; //}}}
  // calculate & sort eigenvalues //{{{
  LONGVECTOR eigen;
  eigen.x = root0; // found out by Newton's method
  // roots of the quadratic equation
  eigen.y = (-b_quad + sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);
  eigen.z = (-b_quad - sqrt(SQR(b_quad) - 4 * a_quad * c_quad)) / (2 * a_quad);

  VECTOR eigen2; // change LONGVECTOR to VECTOR
  eigen2.x = eigen.x;
  eigen2.y = eigen.y;
  eigen2.z = eigen.z;
  eigen2 = SortVector(eigen2); //}}}

  return eigen2;
} //}}}

// memory-freeing functions
void FreeSystem(SYSTEM *System) { //{{{
  free(System->Index_mol);
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
};                                                  //}}}
void FreeMoleculeType(MOLECULETYPE *MoleculeType) { //{{{
  FreeMoleculeTypeEssentials(MoleculeType);
  if (MoleculeType->nBTypes > 0) {
    free(MoleculeType->BType);
  }
  free(MoleculeType->Index);
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

#if 0  //{{{
// TODO redo
// EvaluateContacts() //{{{
/**
 * Function evaluating contacts for aggregate detection for Aggregate and
 * Aggregate-NotSameBeads utilities.
 */
void EvaluateContacts(COUNTS *Counts, AGGREGATE **Aggregate,
                      MOLECULE **Molecule,
                      int contacts, int **contact) {
  // first molecule
  for (int i = 1; i < (*Counts).Molecules; i++) {
    // second molecule
    for (int j = 0; j < i; j++) {
      int agg_i = (*Molecule)[i].Aggregate;
      int agg_j = (*Molecule)[j].Aggregate;
      // molecules 'i' and 'j' are in contact //{{{
      if (contact[i][j] >= contacts) {
        // create new aggregate if 'j' isn'it in any //{{{
        if (agg_j == -1) {
          agg_j = (*Counts).Aggregates;
          (*Molecule)[j].Aggregate = agg_j;

          (*Aggregate)[agg_j].nMolecules = 1;
          (*Aggregate)[agg_j].Molecule[0] = j;

          (*Counts).Aggregates++;
        } //}}}

        // add 'mol_i' to 'agg_j' aggregate (that contains 'mol_j' molecule) if 'i' isn't in an agg //{{{
        if (agg_i == -1) {
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].Molecule[mols] = i;
          (*Aggregate)[agg_j].nMolecules++;

          (*Molecule)[i].Aggregate = agg_j;
        } //}}}

        // 'mol_i' and 'mol_j' molecules are in different aggregate => unite aggregates
        if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {

          // add molecules from aggregate 'agg_i' to 'agg_j' //{{{
          int mols = (*Aggregate)[agg_j].nMolecules;
          (*Aggregate)[agg_j].nMolecules += (*Aggregate)[agg_i].nMolecules;

          // copy molecule ids from Aggregate[agg_i-1] to Aggregate[agg_j-1]
          int id1 = 0;
          for (int k = mols; k < (*Aggregate)[agg_j].nMolecules; k++) {
            int mol = (*Aggregate)[agg_i].Molecule[id1];
            (*Aggregate)[agg_j].Molecule[k] = mol;
            (*Molecule)[mol].Aggregate = agg_j;
            id1++;
          } //}}}

          // move aggregates with id greater then agg_i to id-1 //{{{
          for (int k = (agg_i+1); k < (*Counts).Aggregates; k++) {

            (*Aggregate)[k-1].nMolecules = (*Aggregate)[k].nMolecules;

            // move every molecule from aggregate 'k' to aggregate 'k-1'
            for (int l = 0; l < (*Aggregate)[k].nMolecules; l++) {
              int mol = (*Aggregate)[k].Molecule[l];
              (*Aggregate)[k-1].Molecule[l] = mol;
              (*Molecule)[mol].Aggregate = k - 1;
            }
          } //}}}

          // reduce number of aggregates (two aggregates were merged)
          (*Counts).Aggregates--;
        } //}}}
      } else if (agg_j == -1) { // or 'i' and 'j' aren't in contact and 'j' isn't in any aggregate =>  new aggregate for 'j' */ //{{{
        agg_j = (*Counts).Aggregates;
        (*Molecule)[j].Aggregate = agg_j;

        (*Aggregate)[agg_j].nMolecules = 1;
        (*Aggregate)[agg_j].Molecule[0] = j;

        (*Counts).Aggregates++;
      } //}}}
    }
  }

  // check if highest id residue is in aggregate //{{{
  bool test = false;
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    for (int j = 1; j < (*Aggregate)[i].nMolecules; j++) {
      if ((*Aggregate)[i].Molecule[j] == ((*Counts).Molecules-1)) {
        test = 1;
      }
    }
  } //}}}
  // if highest id residue isn't in any aggregate, create separate one //{{{
  if (!test) {
    int aggs = (*Counts).Aggregates;
    (*Aggregate)[aggs].nMolecules = 1;
    (*Aggregate)[aggs].Molecule[0] = (*Counts).Molecules - 1;

    (*Counts).Aggregates++;
  } //}}}
} //}}}
// SortAggStruct() //{{{
/**
 * Sort an Aggregate struct using the bubble sort algorithm. The resulting
 * struct is arranged so that aggregates with the first molecule's lower id
 * come first.
 */
void SortAggStruct(AGGREGATE **Aggregate, COUNTS Counts,
                   MOLECULE *Molecule, MOLECULETYPE *MoleculeType,
                   BEAD **Bead, BEADTYPE *BeadType) {
  for (int i = 0; i < (Counts.Aggregates-1); i++) {
    bool done = true;
    for (int j = 0; j < (Counts.Aggregates-i-1); j++) {
      if ((*Aggregate)[j].Molecule[0] > (*Aggregate)[j+1].Molecule[0]) {
        SwapInt(&(*Aggregate)[j].nMolecules, &(*Aggregate)[j+1].nMolecules);
        // switch the whole Aggregate[].Molecule array
        int mols; // number of molecules in the larger aggregate
        if ((*Aggregate)[j].nMolecules > (*Aggregate)[j+1].nMolecules) {
          mols = (*Aggregate)[j].nMolecules;
        } else {
          mols = (*Aggregate)[j+1].nMolecules;
        }
        for (int k = 0; k < mols; k++) {
          SwapInt(&(*Aggregate)[j].Molecule[k], &(*Aggregate)[j+1].Molecule[k]);
        }
        // switch bonded beads array
        SwapInt(&(*Aggregate)[j].nBeads, &(*Aggregate)[j+1].nBeads);
        int beads; // number of beads in the larger aggregate
        if ((*Aggregate)[j].nBeads > (*Aggregate)[j+1].nBeads) {
          beads = (*Aggregate)[j].nBeads;
        } else {
          beads = (*Aggregate)[j+1].nBeads;
        }
        for (int k = 0; k < beads; k++) {
          SwapInt(&(*Aggregate)[j].Bead[k], &(*Aggregate)[j+1].Bead[k]);
        }
        // switch monomer beads array
        SwapInt(&(*Aggregate)[j].nMonomers, &(*Aggregate)[j+1].nMonomers);
        int mons; // larger number of monomers of the two aggregates
        if ((*Aggregate)[j].nMonomers > (*Aggregate)[j+1].nMonomers) {
          mons = (*Aggregate)[j].nMonomers;
        } else {
          mons = (*Aggregate)[j+1].nMonomers;
        }
        for (int k = 0; k < mons; k++) {
          SwapInt(&(*Aggregate)[j].Monomer[k], &(*Aggregate)[j+1].Monomer[k]);
        }
        done = false;
      }
    }
    if (done)
      break;
  }

  // re-assign aggregate id to every bonded bead in the aggregate, correcting after sorting //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      int mtype = Molecule[mol].Type;
      for (int k = 0; k < MoleculeType[mtype].nBeads; k++) {
//      int id = Molecule[mol].Bead[k];
//      (*Bead)[id].nAggregates = 1;
//      (*Bead)[id].Aggregatexxx[0] = i;
      }
    }
  } //}}}
} //}}}
// RemovePBCAggregates() //{{{
/**
 * Function to remove periodic boundary conditions from all aggregates,
 * thus joining them.
 */
void RemovePBCAggregates(double distance, AGGREGATE *Aggregate, COUNTS Counts,
                         VECTOR BoxLength, BEADTYPE *BeadType, BEAD **Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {
  // helper array indicating whether molecules already moved
  bool *moved = malloc(sizeof *moved * Counts.Molecules);
  // go through all aggregates larger than unimers and put all molecules together //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {

    // negate moved array, but the first molecule is not to move //{{{
    for (int j = 1; j < Counts.Molecules; j++) {
      moved[j] = false;
    }
    moved[0] = true; //}}}

    bool done = false;
    int test = 0; // if too many loops, just exit with error
    while (!done && test < 1000) {

      // go through all molecule pairs
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        for (int k = 0; k < Aggregate[i].nMolecules; k++) {

          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          if (moved[j] && !moved[k]) { // automatically follows that j != k
            int mol1 = Aggregate[i].Molecule[j],
                mol2 = Aggregate[i].Molecule[k],
                mol1_type = Molecule[mol1].Type,
                mol2_type = Molecule[mol2].Type;

            // go through all bead pairs in the two molecules
            for (int l = 0; l < MoleculeType[mol1_type].nBeads; l++) {
              for (int m = 0; m < MoleculeType[mol2_type].nBeads; m++) {
                int bead1 = Molecule[mol1].Bead[l];
                int bead2 = Molecule[mol2].Bead[m];

                // use only bead types that were used to assign molecules to aggregates
                if (BeadType[(*Bead)[bead1].Type].Use &&
                    BeadType[(*Bead)[bead2].Type].Use) {

                  // calculate distance between 'bead1' and 'bead2'
                  VECTOR dist = Distance((*Bead)[bead1].Position, (*Bead)[bead2].Position, BoxLength);
                  dist.x = Length(dist);

                  // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
                  if (dist.x <= distance) {

                    // distance vector between 'bead1' and 'bead2' //{{{
                    dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z; //}}}
                    // if 'bead1' and 'bead2' are too far in x-direction, move 'mol2' in x-direction //{{{
                    while (dist.x > (BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x += BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    }
                    while (dist.x <= -(BoxLength.x/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.x -= BoxLength.x;
                      }
                      dist.x = (*Bead)[bead1].Position.x - (*Bead)[bead2].Position.x;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in y-direction, move 'mol2' in y-direction //{{{
                    while (dist.y > (BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y += BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    }
                    while (dist.y <= -(BoxLength.y/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.y -= BoxLength.y;
                      }
                      dist.y = (*Bead)[bead1].Position.y - (*Bead)[bead2].Position.y;
                    } //}}}
                    // if 'bead1' and 'bead2' are too far in z-direction, move 'mol2' in x-direction //{{{
                    while (dist.z > (BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z += BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    }
                    while (dist.z <= -(BoxLength.z/2)) {
                      for (int n = 0; n < MoleculeType[mol2_type].nBeads; n++) {
                        (*Bead)[Molecule[mol2].Bead[n]].Position.z -= BoxLength.z;
                      }
                      dist.z = (*Bead)[bead1].Position.z - (*Bead)[bead2].Position.z;
                    } //}}}

                    moved[k] = true;

                    // skip remainder of 'mol2' (or 'k')
                    break;
                  }
                }
              }
              // if molekule 'k' (or 'mol2') has been moved, skip also remainder of molecules 'mol1'
              if (moved[k]) {
                break;
              }
            }
          }
        }
      }

      // check if all molecules have moved //{{{
      done = true;
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        if (!moved[j]) {
          done = false;
          break;
        }
      } //}}}
      test++;
    }
    if (test == 1000) {
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: unable to 'join' aggregate with these ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%d", Aggregate[i].nMolecules);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " molecules:\n");
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        ColourChange(STDERR_FILENO, CYAN);
        fprintf(stderr, " %d", Aggregate[i].Molecule[j]);
      }
      fprintf(stderr, "\n");
      ColourReset(STDERR_FILENO);
    }
  }
  free(moved); //}}}
  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    VECTOR com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead,
                              *Bead, BeadType);

    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    INTVECTOR move;
    move.x = com.x / BoxLength.x;
    move.y = com.y / BoxLength.y;
    move.z = com.z / BoxLength.z;
    if (com.x < 0) {
      move.x--;
    }
    if (com.y < 0) {
      move.y--;
    }
    if (com.z < 0) {
      move.z--;
    }
    // move all the beads
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * BoxLength.x;
      (*Bead)[bead].Position.y -= move.y * BoxLength.y;
      (*Bead)[bead].Position.z -= move.z * BoxLength.z;
    }
  } //}}}
} //}}}
// PrintAggregate() //{{{
/**
 * Function printing Aggregate structure.
 */
void PrintAggregate(COUNTS Counts, int *Index,
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                    BEAD *Bead, BEADTYPE *BeadType, AGGREGATE *Aggregate) {
  fprintf(stdout, "Aggregates: %d\n", Counts.Aggregates);
  for (int i = 0; i < Counts.Aggregates; i++) {
    // print molecules
    fprintf(stdout, " %d mols:", Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      int type = Molecule[mol].Type;
      fprintf(stdout, " %d (%d)", mol, type);
      if (j != (Aggregate[i].nMolecules-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print bonded beads
    fprintf(stdout, " %d bonded beads:", Aggregate[i].nBeads);
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      fprintf(stdout, " %d", bead);
      if (j != (Aggregate[i].nBeads-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
    // print monomeric beads
    fprintf(stdout, " %d free beads:", Aggregate[i].nMonomers);
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int bead = Aggregate[i].Monomer[j];
      fprintf(stdout, " %d", bead);
      if (j != (Aggregate[i].nMonomers-1)) {
        putchar(',');
      } else {
        putchar('\n');
      }
    }
  }
} //}}}
// FreeAggregate() //{{{
/**
 * Free memory allocated for Aggregate struct array. This function makes it
 * easier other arrays to the Aggregate struct in the future
 */
void FreeAggregate(COUNTS Counts, AGGREGATE **Aggregate) {
  for (int i = 0; i < Counts.Molecules; i++) {
    free((*Aggregate)[i].Molecule);
    free((*Aggregate)[i].Bead);
    free((*Aggregate)[i].Monomer);
  }
  free(*Aggregate);
} //}}}
#endif //}}}
