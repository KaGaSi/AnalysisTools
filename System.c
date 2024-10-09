#include "System.h"

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
}
void ReFillMoleculeTypeBType(SYSTEM *System) {
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    free(System->MoleculeType[i].BType);
    FillMoleculeTypeBType(&System->MoleculeType[i]);
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
// void FillBeadTypeIndex(SYSTEM *System) { //{{{
//   COUNT *Count = &System->Count;
//   // allocate memory for Index arrays
//   for (int i = 0; i < Count->BeadType; i++) {
//     BEADTYPE *bt = &System->BeadType[i];
//     if (bt->Number > 0) {
//       bt->Index = malloc(sizeof *bt->Index * bt->Number);
//     }
//   }
//   // fill the Index arrays
//   int *count_id = calloc(Count->BeadType, sizeof *count_id);
//   for (int i = 0; i < Count->Bead; i++) {
//     int type = System->Bead[i].Type;
//     if (System->BeadType[type].Number < count_id[type]) {
//       fprintf(stderr, "...hmm; error count_id[%d]=%d (%d)\n", type,
//               count_id[type], System->BeadType[type].Number);
//     }
//     System->BeadType[type].Index[count_id[type]] = i;
//     count_id[type]++;
//   }
//   free(count_id);
// } //}}}
// TODO: bt->Index alloc must be changed if it's supposed to be run whenever
//       coordinates are read; two ways:
//       a) add a switch saying when to allocate (only when called through some
//          structure-reading function)
//       b) put the allocation directly into the structure-reading functions
//          (where other allocations are - well, should be, I think)
void FillBeadTypeIndex(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  // allocate memory for Index arrays
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt = &System->BeadType[i];
    bt->InCoor = 0;
  }
  // fill the Index arrays
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System->BeadCoor[i],
        type = System->Bead[id].Type;
    BEADTYPE *bt = &System->BeadType[type];
    if (bt->Number <= bt->InCoor) {
      fprintf(stderr, "...hmm; error count_id[%d]=%d (%d), bead %d\n", type,
              bt->InCoor, bt->Number, id);
    }
    bt->Index[bt->InCoor] = id;
    bt->InCoor++;
  }
}
void AllocFillBeadTypeIndex(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadType; i++) {
    BEADTYPE *bt = &System->BeadType[i];
    if (bt->Number > 0) {
      bt->Index = calloc(bt->Number, sizeof *bt->Index);
    }
  }
  FillBeadTypeIndex(System);
}
void RefillBeadTypeIndex(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadType; i++) {
    BEADTYPE *bt = &System->BeadType[i];
    if (bt->Number > 0) {
      free(bt->Index);
    }
  }
  AllocFillBeadTypeIndex(System);
} //}}}
void FillMoleculeTypeIndex(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    if (mt->Number > 0) {
      mt->Index = malloc(sizeof *mt->Index * mt->Number);
    }
  }
  int *count_id = calloc(Count->MoleculeType, sizeof *count_id);
  for (int i = 0; i < Count->Molecule; i++) {
    int type = System->Molecule[i].Type;
    MOLECULETYPE *mt = &System->MoleculeType[type];
    if (System->MoleculeType[type].Number < count_id[type]) {
      fprintf(stderr, "...hmm; %s error count_id[%d]=%d (%d)\n",
              System->MoleculeType[type].Name, type, count_id[type],
              System->MoleculeType[type].Number);
    }
    mt->Index[count_id[type]] = i;
    count_id[type]++;
  }
  free(count_id);
}
void ReFillMoleculeTypeIndex(SYSTEM *System) {
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    free(System->MoleculeType[i].Index);
  }
  FillMoleculeTypeIndex(System);
} //}}}
// void FillIndexMol(SYSTEM *System) { //{{{
//   COUNT *Count = &System->Count;
//   if (Count->Molecule > 0) {
//     System->Index_mol = realloc(System->Index_mol,
//                                 sizeof *System->Index_mol *
//                                 (Count->HighestResid + 1));
//     for (int i = 0; i <= Count->HighestResid; i++) {
//       System->Index_mol[i] = -1;
//     }
//     for (int i = 0; i < Count->Molecule; i++) {
//       System->Index_mol[System->Molecule[i].Index] = i;
//     }
//   }
// } //}}}
void FillBondedUnbonded(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  if (Count->Bonded > 0) {
    System->Bonded = realloc(System->Bonded,
                             sizeof *System->Bonded * Count->Bonded);
    System->BondedCoor = realloc(System->BondedCoor,
                                 sizeof *System->BondedCoor * Count->Bonded);
  }
  if (Count->Unbonded > 0) {
    System->Unbonded = realloc(System->Unbonded,
                               sizeof *System->Unbonded * Count->Unbonded);
    System->UnbondedCoor = realloc(System->UnbondedCoor,
                                   sizeof *System->UnbondedCoor *
                                   Count->Unbonded);
  }
  int c_bonded = 0, c_unbonded = 0;
  for (int i = 0; i < Count->Bead; i++) {
    if (System->Bead[i].Molecule > -1) {
      System->Bonded[c_bonded] = i;
      c_bonded++;
    } else {
      System->Unbonded[c_unbonded] = i;
      c_unbonded++;
    }
  }
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
void FillSystemNonessentials(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  AllocFillBeadTypeIndex(System);
  FillMoleculeTypeIndex(System);
  FillBondedUnbonded(System);
  CountBondAngleDihedralImproper(System);
  // sort bonds, angles, dihedrals, and impropers
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    SortBonds(mt_i->Bond, mt_i->nBonds);
    SortAngles(mt_i->Angle, mt_i->nAngles);
    SortDihImp(mt_i->Dihedral, mt_i->nDihedrals);
    SortDihImp(mt_i->Improper, mt_i->nImpropers);
  }
} //}}}
void FillInCoor(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->BeadType; i++) {
    System->BeadType[i].InCoor = 0;
  }
  Count->BondedCoor = 0;
  Count->UnbondedCoor = 0;
  Count->MoleculeCoor = 0;
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System->BeadCoor[i];
    BEAD *b = &System->Bead[id];
    b->InTimestep = true;
    if (b->Molecule != -1) {
      if (!System->Molecule[b->Molecule].InTimestep) {
        System->Molecule[b->Molecule].InTimestep = true;
        System->MoleculeCoor[Count->MoleculeCoor] = b->Molecule;
        Count->MoleculeCoor++;
      }
      System->BondedCoor[Count->BondedCoor] = id;
      Count->BondedCoor++;
    } else {
      System->UnbondedCoor[Count->UnbondedCoor] = id;
      Count->UnbondedCoor++;
    }
    BEADTYPE *bt = &System->BeadType[b->Type];
    bt->Index[bt->InCoor] = id;
    bt->InCoor++;
  }
} //}}}
// fill in BOX structure //{{{
/*
* mode: 0 ... knowing angles
*       1 ... knowing tilt vector
*/
bool CalculateBoxData(BOX *Box, int mode) {
  // calculate angles and tilt vectors or tilt vectors and OrthoLength //{{{
  switch (mode) {
    case 0: // angles & Length given //{{{
      for (int dd = 0; dd < 3; dd++) {
        Box->OrthoLength[dd] = Box->Length[dd];
      }
      Box->Volume = Box->Length[0] * Box->Length[1] *  Box->Length[2];
      if (Box->alpha != 90 || Box->beta != 90 || Box->gamma != 90) {
        double a = Box->Length[0], b = Box->Length[1], c = Box->Length[2];
        double c_a = cos(Box->alpha * PI / 180),
               c_b = cos(Box->beta * PI / 180),
               c_g = cos(Box->gamma * PI / 180),
               s_g = sin(Box->gamma * PI / 180);
        // cell volume
        double sqr = 1 - SQR(c_a) - SQR(c_b) - SQR(c_g) + 2 * c_a * c_b * c_g;
        if (sqr < 0) {
          strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
          return false;
        }
        Box->Volume *= sqrt(sqr);
        // transformation matrix fractional -> Cartesian coordinates
        double vol = Box->Volume;
        Box->transform[0][0] = a;
        Box->transform[0][1] = b * c_g;
        Box->transform[1][1] = b * s_g;
        Box->transform[0][2] = c * c_b;
        Box->transform[1][2] = c * (c_a - c_b * c_g) / s_g;
        Box->transform[2][2] = vol / (a * b * s_g);
        // transformation matrix Cartesian -> fractional coordinates
        Box->inverse[0][0] = 1 / a;
        Box->inverse[0][1] = -c_g / (a * s_g);
        Box->inverse[1][1] = 1 / (b * s_g);
        Box->inverse[0][2] =
            b * c * (c_g * (c_a - c_b * c_g) / (s_g * vol) - c_b * s_g / vol);
        Box->inverse[1][2] = -a * c * (c_a - c_b * c_g) / (vol * s_g);
        Box->inverse[2][2] = a * b * s_g / vol;
        // orthogonal box size
        // x direaction
        Box->OrthoLength[0] = a;
        // y direaction
        sqr = SQR(b) - SQR(Box->transform[0][1]);
        if (sqr < 0) {
          strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
          return false;
        }
        Box->OrthoLength[1] = sqrt(sqr);
        // z direaction
        sqr = SQR(c) - SQR(Box->transform[0][2]) - SQR(Box->transform[1][2]);
        if (sqr < 0) {
          strcpy(ERROR_MSG, "wrong simulation box dimensions");
          PrintError();
          exit(1);
        }
        Box->OrthoLength[2] = sqrt(sqr);
        // see https://docs.lammps.org/Howto_triclinic.html
        double xy = Box->transform[0][1], xz = Box->transform[0][2],
               yz = Box->transform[1][2],
               xyz = Box->transform[0][1] + Box->transform[0][2];
        Box->Bounding[0] = Box->OrthoLength[0] -
                           Max3(0, xy, Max3(0, xz, xyz)) +
                           Min3(0, xy, Min3(0, xz, xyz));
        Box->Bounding[1] = Box->OrthoLength[1] -
                           Min3(0, 0, yz) + Max3(0, 0, yz);
        Box->Bounding[2] = Box->OrthoLength[2];
      }
      break; //}}}
    case 1:  // tilt & OrthoLength given //{{{
      for (int dd = 0; dd < 3; dd++) {
        Box->Length[dd] = Box->OrthoLength[dd];
      }
      Box->Volume = Box->Length[0] * Box->Length[1] *  Box->Length[2];
      if (Box->transform[0][1] != 0 ||
          Box->transform[0][2] != 0 ||
          Box->transform[1][2] != 0) {
        double a = Box->OrthoLength[0],
               b = sqrt(SQR(Box->OrthoLength[1]) + SQR(Box->transform[0][1])),
               c = sqrt(SQR(Box->OrthoLength[2]) +
                        SQR(Box->transform[0][2]) +
                        SQR(Box->transform[1][2]));
        double c_a = (Box->transform[0][1] * Box->transform[0][2] +
                      Box->OrthoLength[1] * Box->transform[1][2]) /
                     (b * c),
               c_b = Box->transform[0][2] / c,
               c_g = Box->transform[0][1] / b,
               s_g = sin(Box->gamma * PI / 180);
        // cell length
        Box->Length[0] = a;
        Box->Length[1] = b;
        Box->Length[2] = c;
        // cell angles
        Box->alpha = acos(c_a) / PI * 180;
        Box->beta = acos(c_b) / PI * 180;
        Box->gamma = acos(c_g) / PI * 180;
        // cell volume
        double sqr = 1 - SQR(c_a) - SQR(c_b) - SQR(c_g) + 2 * c_a * c_b * c_g;
        if (sqr < 0) {
          strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
          return false;
        }
        Box->Volume *= sqrt(sqr);
        // finish transformation matrix fractional -> Cartesian coordinates
        double vol = Box->Volume;
        Box->transform[0][0] = a;
        Box->transform[1][1] = b * s_g;
        Box->transform[2][2] = vol / (a * b * s_g);
        // transformation matrix Cartesian -> fractional coordinates
        Box->inverse[0][0] = 1 / a;
        Box->inverse[0][1] = -c_g / (a * s_g);
        Box->inverse[1][1] = 1 / (b * s_g);
        Box->inverse[0][2] =
            b * c * (c_g * (c_a - c_b * c_g) / (s_g * vol) - c_b * s_g / vol);
        Box->inverse[1][2] = -a * c * (c_a - c_b * c_g) / (vol * s_g);
        Box->inverse[2][2] = a * b * s_g / vol;
      }
      break; //}}}
    default:
      strcpy(ERROR_MSG, "TriclinicCellData(): mode parameters must be 0 or 1");
      PrintError();
      exit(1);
  } //}}}
  // transformation matrices for orthogonal box //{{{
  if (Box->alpha == 90 && Box->beta == 90 && Box->gamma == 90) {
    Box->Volume = Box->Length[0] * Box->Length[1] * Box->Length[2];
    Box->transform[0][0] = Box->Length[0];
    Box->transform[1][1] = Box->Length[1];
    Box->transform[2][2] = Box->Length[2];
    Box->inverse[0][0] = 1 / Box->Length[0];
    Box->inverse[1][1] = 1 / Box->Length[1];
    Box->inverse[2][2] = 1 / Box->Length[2];
  } //}}}
  // maximum size of the the bounding box //{{{
  // see https://docs.lammps.org/Howto_triclinic.html
  double xy = Box->transform[0][1],
         xz = Box->transform[0][2],
         yz = Box->transform[1][2],
         xyz = Box->transform[0][1] + Box->transform[0][2];
  Box->Bounding[0] = Box->OrthoLength[0] -
                     Max3(0, xy, Max3(0, xz, xyz)) +
                     Min3(0, xy, Min3(0, xz, xyz));
  Box->Bounding[1] = Box->OrthoLength[1] -
                     Min3(0, 0, yz) + Max3(0, 0, yz);
  Box->Bounding[2] = Box->OrthoLength[2]; //}}}
  if (Box->Volume == 0) { //{{{
    strcpy(ERROR_MSG, "not all box dimensions are non-zero:");
    PrintError();
    PrintBox(*Box);
    exit(1);
  } //}}}
  return true;
} //}}}

// merge identical bead/molecule types
// MergeBeadTypes() //{{{
/* Merge bead types either using names only (detailed=false) or using names,
 * charge, mass, and radius (detailed=true). In the simpler case, the values
 * for charge mass, and radius are each taken from the first bead type that has
 * the corresponding value well defined; the more complicated case is described
 * below in some detail.
 * The function assumes BeadType[].Index arrays are unallocated (i.e., nothing
 * is freed for the remaining extra BeadType array elements).
 */
/* detailed=true description: //{{{
 * First, identify bead types based on name, charge, masse, and radius,
 * e.g., lines
 *   atom 0 n x q 1 m 1
 *   atom 1 n x q 2 m 1
 * will be of two different types. This can create an excess of bead types,
 * so some may have to be merged.
 *
 * What is to be merged:
 * i) If a keyword is missing in one line but present in another, that does
 * not count as a different type, e.g., lines
 *       atom 0 n x q 1 m 1
 *       atom 1 n x     m 1
 *    are of the same type (both with charge +1);
 * ii) however, there can be ambiguities, so e.g., lines
 *        atom 0 n x q 1 m 1
 *        atom 1 n x     m 1
 *        atom 2 n x q 0 m 1
 *     remain three distinct types (atom 1 has undefined charge);
 * iii) but only some lines may be ambiguous, e.g., lines
 *        atom 0 n x q 1 m 1
 *        atom 1 n x     m 1
 *        atom 2 n x q 0 m 1
 *        atom 3 n x q 0
 *      are still three different types (the last two should be considered
 *      the same because there is no ambiguity because all beads have the
 *      same mass)
 * iv) note that sometimes the charge/mass/radius can remain undefined
 *     even though there's only one well defined value; e.g., lines
 *       atom 0 n x q 1 m 1
 *       atom 1 n x     m 1
 *       atom 2 n x q 0 m 1
 *       atom 3 n x q 0
 *       atom 4 n x q 0 m 1 r 1
 *     will make radius well defined (with value 1) only for beads sharing
 *     the type with atom 4 (i.e., atoms 2, 3, and 4), while the first two
 *     atoms will still have undefined radius. What should the radius of
 *     atoms 0 and 1 be when the charge is different/unspecified to that of the
 *     last atom?
 *
 * Merging procedure:
 * 1) for each unique name, find values of charge/mass/radius, noting
 *    ambiguities (i.e., when more than one well defined value exists)
 * 2) create 2D boolean array of size <unique names>*<unique names> to
 *    see what should be merged based on points 1) and 2):
 *    i) pick two bead types sharing a name (or the same bead type twice
 *       if it does not share a name with any other), say 'i' and 'k'.
 *    ii) check every bead type (say 'j') against i and k; if i and
 *        j should be merged (i.e., share a name), check k's value of
 *        diff_q/m/r - if it is a proper value, merge i and j; if not,
 *        merge i and j only if they have the same diff_q/m/r value.
 * 3) merge the types and count the number of unique types
 *    i) create a new type when a diagonal element of the array is true
 *    ii) check the remaining types against and merge those that should
 *        be merge with that new type, making the diagonal element for
 *        that merged type false so that no new type is created when
 *        its time comes in i)
 * 4) reorder the types so that types sharing the name are next to each
 *    other
 */ //}}}
void MergeBeadTypes(SYSTEM *System, bool detailed) {

  COUNT *Count = &System->Count;
  int count_bt_old = Count->BeadType,
      *old_to_new = calloc(count_bt_old, sizeof *old_to_new);

  // merge bead types that are definitely the same //{{{
  int count_bt_new = 1;
  for (int i = 1; i < count_bt_old; i++) { // 1 as bt[0] remains bt[0]
    BEADTYPE *bt_old = &System->BeadType[i];
    bool new = true;
    for (int j = 0; j < count_bt_new; j++) {
      BEADTYPE *bt_new = &System->BeadType[j];
      if (SameBeadType(*bt_old, *bt_new)) {
        bt_new->Number += bt_old->Number;
        old_to_new[i] = j;
        new = false;
        break;
      }
    }
    if (new) {
      System->BeadType[count_bt_new] = *bt_old;
      old_to_new[i] = count_bt_new;
      count_bt_new++;
    }
  } //}}}

  // relabel bead types in arrays //{{{
  // Bead[].Type
  for (int i = 0; i < Count->Bead; i++) {
    int old_type = System->Bead[i].Type;
    System->Bead[i].Type = old_to_new[old_type];
  }
  // MoleculeType[].Bead[]
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < System->MoleculeType[i].nBeads; j++) {
      int old_type = System->MoleculeType[i].Bead[j];
      System->MoleculeType[i].Bead[j] = old_to_new[old_type];
    }
  } //}}}

  Count->BeadType = count_bt_new;
  count_bt_old = count_bt_new;
  old_to_new = realloc(old_to_new, sizeof *old_to_new * count_bt_old);

  // find the unique bead names //{{{
  int count_bnames = 0;
  char(*bname)[BEAD_NAME] = malloc(sizeof *bname);
  for (int i = 0; i < Count->BeadType; i++) {
    bool new = true;
    for (int j = 0; j < count_bnames; j++) {
      if (strcmp(System->BeadType[i].Name, bname[j]) == 0) {
        new = false;
        break;
      }
    }
    if (new) {
      int n = count_bnames;
      count_bnames++;
      bname = realloc(bname, sizeof *bname * count_bnames);
      strncpy(bname[n], System->BeadType[i].Name, BEAD_NAME);
    }
  }               //}}}
  if (detailed) { // use name as well as charge, mass, and radius
    // 1) find charge/mass/radius values for each unique name //{{{
    // arrays values of charge, mass, and radius for each bead type
    double diff_q[count_bnames], diff_m[count_bnames], diff_r[count_bnames];
    // initialize arrays: assign values from the last type with each name //{{{
    for (int i = 0; i < count_bnames; i++) {
      for (int j = 0; j < Count->BeadType; j++) {
        if (strcmp(bname[i], System->BeadType[j].Name) == 0) {
          diff_q[i] = System->BeadType[j].Charge;
          diff_m[i] = System->BeadType[j].Mass;
          diff_r[i] = System->BeadType[j].Radius;
          break;
        }
      }
    } //}}}
    // find the proper values for charge/mass/radius for each type //{{{
    /*
     * diff_q/m/r = high ... more than one value for beads with that name
     * diff_q/m/r = <value> ... exactly that one value;
     *                          if both proper and undefined values exist
     *                          (i.e., when there's really one value, but it's
     *                          not written in each atom line), the proper
     *                          value is assigned
     */
    // high, impossible number to indicate multiple values of charge/mass/radius
    int high = 1e6;
    // go through all bead type pairs (including self-pairs)
    for (int i = 0; i < count_bnames; i++) {
      for (int j = 0; j < Count->BeadType; j++) {
        BEADTYPE *bt_j = &System->BeadType[j];
        // only consider type pairs with the same name
        if (strcmp(bname[i], bt_j->Name) == 0) {
          // charge
          if (diff_q[i] != bt_j->Charge) {
            if (diff_q[i] != CHARGE && bt_j->Charge != CHARGE) {
              diff_q[i] = high;
            } else if (diff_q[i] == CHARGE) {
              diff_q[i] = bt_j->Charge;
            }
          }
          // mass
          if (diff_m[i] != bt_j->Mass) {
            if (diff_m[i] != MASS && bt_j->Mass != MASS) {
              diff_m[i] = high;
            } else if (diff_m[i] == MASS) {
              diff_m[i] = bt_j->Mass;
            }
          }
          // radius
          if (diff_r[i] != bt_j->Radius) {
            if (diff_r[i] != RADIUS && bt_j->Radius != RADIUS) {
              diff_r[i] = high;
            } else if (diff_r[i] == RADIUS) {
              diff_r[i] = bt_j->Radius;
            }
          }
        }
      }
    } //}}}
    //}}}
    // 2) create 2D merge array //{{{
    // initialize merge array by assuming nothing will be merged
    bool **merge = malloc(sizeof *merge * Count->BeadType);
    for (int i = 0; i < Count->BeadType; i++) {
      merge[i] = malloc(sizeof *merge[i] * Count->BeadType);
      for (int j = 0; j < Count->BeadType; j++) {
        merge[i][j] = false; // 'i' and 'j' aren't to be merged
        merge[i][i] = true;  // 'i' and 'i' are to be merged/copied
      }
    }
    // assume same-name bead types are to be merged
    for (int i = 0; i < (Count->BeadType - 1); i++) {
      for (int j = (i + 1); j < Count->BeadType; j++) {
        if (strcmp(System->BeadType[i].Name, System->BeadType[j].Name) == 0) {
          merge[i][j] = true;
        }
      }
    }
    // go through each bead type and compare it to two others
    for (int i = 0; i < (Count->BeadType - 1); i++) {
      BEADTYPE *bt_i = &System->BeadType[i];
      // i)
      int k = 0;
      for (; k < count_bnames; k++) {
        if (strcmp(bname[k], bt_i->Name) == 0) {
          break;
        }
      }
      // ii)
      for (int j = (i + 1); j < Count->BeadType; j++) {
        BEADTYPE *bt_j = &System->BeadType[j];
        // check charge //{{{
        if (merge[i][j]) {
          if (diff_q[k] == high) {
            if (bt_i->Charge == bt_j->Charge) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        } //}}}
        // check mass //{{{
        if (merge[i][j]) {
          if (diff_m[k] == high) {
            if (bt_i->Mass == bt_j->Mass) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        } //}}}
        // check radius //{{{
        if (merge[i][j]) {
          if (diff_r[k] == high) {
            if (bt_i->Radius == bt_j->Radius) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        } //}}}
      }
    } //}}}
    // 3) merge the types, counting the new types //{{{
    BEADTYPE *temp = calloc(Count->BeadType, sizeof *temp);
    int old_bt_count = Count->BeadType, count = 0,
        *bt_older_to_old = calloc(old_bt_count, sizeof *bt_older_to_old);
    for (int i = 0; i < Count->BeadType; i++) {
      if (merge[i][i]) { // i)
        temp[count] = System->BeadType[i];
        bt_older_to_old[i] = count;
        for (int j = (i + 1); j < Count->BeadType; j++) {
          BEADTYPE *bt_j = &System->BeadType[j];
          if (merge[i][j]) { // ii)
            temp[count].Number += bt_j->Number;
            bt_older_to_old[j] = count;
            if (temp[count].Charge == CHARGE) {
              temp[count].Charge = bt_j->Charge;
            }
            if (temp[count].Mass == MASS) {
              temp[count].Mass = bt_j->Mass;
            }
            if (temp[count].Radius == RADIUS) {
              temp[count].Radius = bt_j->Radius;
            }
            merge[j][j] = false;
          }
        }
        count++;
      }
    }
    for (int i = 0; i < Count->BeadType; i++) {
      free(merge[i]);
    }
    free(merge);
    Count->BeadType = count; //}}}
    // // 4) add 'bt' to unnamed types //{{{
    // TODO: this shouldn't have any effect now, right?
    // for (int i = 0; i < Count->BeadType; i++) {
    //   if (strcmp(System->BeadType[i].Name, NON) == 0) {
    //     strcpy(System->BeadType[i].Name, "bt");
    //   }
    // } //}}}
    // 5) reorder the types, placing same-named ones next to each other //{{{
    // copy all bead types temporarily to bt struct
    int *bt_old_to_new = malloc(sizeof *bt_old_to_new * Count->BeadType);
    bool *copied = calloc(Count->BeadType, sizeof *copied);
    for (int i = 0; i < Count->BeadType; i++) {
      System->BeadType[i] = temp[i];
    }
    // copy the bead types back to temp array in a proper order
    count = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt_i = &System->BeadType[i];
      if (!copied[i]) {
        temp[count] = *bt_i;
        bt_old_to_new[i] = count;
        count++;
        copied[i] = true;
        for (int j = (i + 1); j < Count->BeadType; j++) {
          BEADTYPE *bt_j = &System->BeadType[j];
          if (strcmp(bt_i->Name, bt_j->Name) == 0 && !copied[j]) {
            temp[count] = *bt_j;
            bt_old_to_new[j] = count;
            copied[j] = true;
            count++;
          }
        }
      }
    }
    free(copied);
    // finally, copy the types from the temporary array back to bt array
    for (int i = 0; i < Count->BeadType; i++) {
      System->BeadType[i] = temp[i];
    }
    free(temp);
    //}}}
    RenameBeadTypes(System);
    // fill array to relabel bead types in arrays //{{{
    for (int i = 0; i < count_bt_old; i++) {
      old_to_new[i] = bt_old_to_new[bt_older_to_old[i]];
    }
    free(bt_older_to_old);
    free(bt_old_to_new); //}}}
  } else { // use name only
    Count->BeadType = count_bnames;
    // go through old types, creating a new one for each unique name //{{{
    int count_bt_new = 0;
    for (int i = 0; i < count_bt_old; i++) {
      BEADTYPE *bt_i = &System->BeadType[i];
      bool new = true;
      int j = 0;
      for (; j < count_bt_new; j++) {
        BEADTYPE *bt_j = &System->BeadType[j];
        if (strcmp(bt_i->Name, bt_j->Name) == 0) {
          if (i != j) {
            bt_j->Number += bt_i->Number;
            // assign charge/mass/radius, if yet undefined
            if (bt_j->Charge == CHARGE) {
              bt_j->Charge = bt_i->Charge;
            }
            if (bt_j->Mass == MASS) {
              bt_j->Mass = bt_i->Mass;
            }
            if (bt_j->Radius == RADIUS) {
              bt_j->Radius = bt_i->Radius;
            }
          }
          old_to_new[i] = j;
          new = false;
          break;
        }
      }
      if (new) {      // create new type...
        if (i != j) { // ...unless its id is the same as the old one's
          strncpy(System->BeadType[count_bt_new].Name, bt_i->Name, BEAD_NAME);
          System->BeadType[count_bt_new] = *bt_i;
        }
        old_to_new[i] = count_bt_new;
        count_bt_new++;
      }
    } //}}}
  }
  free(bname);
  // relabel bead types in arrays //{{{
  // Bead[].Type
  for (int i = 0; i < Count->Bead; i++) {
    int old_type = System->Bead[i].Type;
    System->Bead[i].Type = old_to_new[old_type];
  }
  // MoleculeType[].Bead[]
  for (int i = 0; i < Count->MoleculeType; i++) {
    for (int j = 0; j < System->MoleculeType[i].nBeads; j++) {
      int old_type = System->MoleculeType[i].Bead[j];
      System->MoleculeType[i].Bead[j] = old_to_new[old_type];
    }
  }
  free(old_to_new); //}}}
  // warning - test count bead types; should never happen //{{{
  int *count_test = calloc(Count->BeadType, sizeof *count_test);
  for (int i = 0; i < Count->Bead; i++) {
    int type = System->Bead[i].Type;
    count_test[type]++;
  }
  for (int i = 0; i < Count->BeadType; i++) {
    if (count_test[i] != System->BeadType[i].Number) {
      strcpy(ERROR_MSG, "something went wrong with bead type differentiation; "
                        "this should never happen!");
      PrintWarning();
      fprintf(stderr, "%sBead count for %s%s%s type: %s%d%s and %s%d%s\n",
              ErrCyan(), ErrYellow(), System->BeadType[i].Name, ErrCyan(),
              ErrYellow(), System->BeadType[i].Number, ErrCyan(), ErrYellow(),
              count_test[i], ErrColourReset());
    }
  }
  free(count_test); //}}}
} //}}}
// MergeMoleculeTypes() //{{{
/*
 * Molecules of one type must share:
 * i) molecule name and numbers of beads, bonds, angles, dihedrals,
 *    and impropers
 * ii) order of bead types
 * iii) connectivity
 * iv) same angles, dihedrals & impropers
 */
void MergeMoleculeTypes(SYSTEM *System) {
  COUNT *Count = &System->Count;
  int count = 0, *old_to_new = calloc(Count->MoleculeType, sizeof *old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) { //{{{
    bool new = true;
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    SortBonds(mt_i->Bond, mt_i->nBonds);
    SortAngles(mt_i->Angle, mt_i->nAngles);
    SortDihImp(mt_i->Dihedral, mt_i->nDihedrals);
    SortDihImp(mt_i->Improper, mt_i->nImpropers);
    int j = 0;
    for (; j < count; j++) {
      MOLECULETYPE *mt_j = &System->MoleculeType[j];
      // i) check numbers of stuff
      if (strcmp(mt_i->Name, mt_j->Name) == 0 &&
          mt_i->nBeads == mt_j->nBeads &&
          mt_i->nAngles == mt_j->nAngles &&
          mt_i->nDihedrals == mt_j->nDihedrals &&
          mt_i->nImpropers == mt_j->nImpropers) {
        // ii) check bead order
        bool same_mol = true; // assume i and j are the same molecule
        for (int k = 0; k < mt_i->nBeads; k++) {
          if (mt_i->Bead[k] != mt_j->Bead[k]) {
            same_mol = false; // i and j aren't the same
          }
        }
        if (!same_mol) {
          continue;
        }
        // iii) check bonds, angles, etc.
        // bonds
        for (int k = 0; k < mt_j->nAngles; k++) {
          if (mt_i->Angle[k][0] != mt_j->Angle[k][0] ||
              mt_i->Angle[k][1] != mt_j->Angle[k][1]) {
            same_mol = false; // i and j aren't the same
            break;
          }
        }
        if (!same_mol) {
          continue;
        }
        // angles
        for (int k = 0; k < mt_j->nAngles; k++) {
          if (mt_i->Angle[k][0] != mt_j->Angle[k][0] ||
              mt_i->Angle[k][1] != mt_j->Angle[k][1] ||
              mt_i->Angle[k][2] != mt_j->Angle[k][2]) {
            same_mol = false; // i and j aren't the same
            break;
          }
        }
        if (!same_mol) {
          continue;
        }
        // dihedrals
        for (int k = 0; k < mt_j->nDihedrals; k++) {
          if (mt_i->Dihedral[k][0] != mt_j->Dihedral[k][0] ||
              mt_i->Dihedral[k][1] != mt_j->Dihedral[k][1] ||
              mt_i->Dihedral[k][2] != mt_j->Dihedral[k][2] ||
              mt_i->Dihedral[k][3] != mt_j->Dihedral[k][3]) {
            same_mol = false; // i and j aren't the same
            break;
          }
        }
        if (!same_mol) {
          continue;
        }
        // impropers
        for (int k = 0; k < mt_j->nImpropers; k++) {
          if (mt_i->Improper[k][0] != mt_j->Improper[k][0] ||
              mt_i->Improper[k][1] != mt_j->Improper[k][1] ||
              mt_i->Improper[k][2] != mt_j->Improper[k][2] ||
              mt_i->Improper[k][3] != mt_j->Improper[k][3]) {
            same_mol = false; // i and j aren't the same
            break;
          }
        }
        if (!same_mol) {
          continue;
        }
        // are molecule types i and j the same?
        if (same_mol) {
          if (i != j) {
            mt_j->Number += mt_i->Number;
            FreeMoleculeTypeEssentials(mt_i);
          }
          old_to_new[i] = j;
          new = false;
          break;
        }
      }
    }
    if (new) {      // create new type...
      if (i != j) { // ...unless its id is the same as the old one's
        strncpy(System->MoleculeType[count].Name, mt_i->Name, MOL_NAME);
        System->MoleculeType[count] = CopyMoleculeTypeEssentials(*mt_i);
        FreeMoleculeTypeEssentials(mt_i);
      }
      old_to_new[i] = count;
      count++;
    }
  }
  Count->MoleculeType = count; //}}}
  // relabel molecules with proper molecule types //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    int old_type = System->Molecule[i].Type;
    System->Molecule[i].Type = old_to_new[old_type];
  }
  free(old_to_new); //}}}
  // reorder the types, placing same-named ones next to each other //{{{
  // copy all molecule types to a temporary array & free the original array
  MOLECULETYPE *temp = malloc(sizeof(MOLECULETYPE) * Count->MoleculeType);
  for (int i = 0; i < Count->MoleculeType; i++) {
    temp[i] = CopyMoleculeTypeEssentials(System->MoleculeType[i]);
    temp[i].Flag = false;
    FreeMoleculeTypeEssentials(&System->MoleculeType[i]);
  }
  free(System->MoleculeType);
  System->MoleculeType = malloc(sizeof *System->MoleculeType *
                                Count->MoleculeType);
  // array to link old molecule type indices to new ones
  old_to_new = malloc(sizeof *old_to_new * Count->MoleculeType);
  // copy the molecule types back in a proper order
  count = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!temp[i].Flag) {
      System->MoleculeType[count] = CopyMoleculeTypeEssentials(temp[i]);
      old_to_new[i] = count;
      count++;
      temp[i].Flag = true;
      for (int j = (i + 1); j < Count->MoleculeType; j++) {
        if (strcmp(System->MoleculeType[i].Name, temp[j].Name) == 0 &&
            !temp[j].Flag) {
          System->MoleculeType[count] = CopyMoleculeTypeEssentials(temp[j]);
          old_to_new[j] = count;
          count++;
          temp[j].Flag = true;
        }
      }
    }
  }
  // free the temporary array
  for (int i = 0; i < Count->MoleculeType; i++) {
    FreeMoleculeTypeEssentials(&temp[i]);
  } //}}}
  free(temp);
  RenameMoleculeTypes(System);
  // relabel molecules with proper molecule types - again // TODO: why? //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    int old_type = System->Molecule[i].Type;
    System->Molecule[i].Type = old_to_new[old_type];
  }
  free(old_to_new); //}}}
} //}}}

// sort a bond/angle/dihedral/improper array in an ascending order
void SortBonds(int (*bond)[3], int n) { //{{{
  // first, check order in every bond
  for (int j = 0; j < n; j++) {
    if (bond[j][0] > bond[j][1]) {
      SwapInt(&bond[j][0], &bond[j][1]);
    }
  }
  // second, bubble sort bonds
  for (int j = 0; j < (n - 1); j++) {
    bool swap = false;
    for (int k = 0; k < (n - j - 1); k++) {
      /* swap bonds if
       * 1) first bead ids are in the wrong order or
       * 2) first bead ids are fine, but second ids are in the wrong order
       */
      if (bond[k][0] > bond[k + 1][0] ||   // 1)
          (bond[k][0] == bond[k + 1][0] && // 2)
           bond[k][1] > bond[k + 1][1])) { // 2)
        SwapInt(&bond[k][0], &bond[k + 1][0]);
        SwapInt(&bond[k][1], &bond[k + 1][1]);
        SwapInt(&bond[k][2], &bond[k + 1][2]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}
void SortAngles(int (*angle)[4], int n) { //{{{
  // first, check order of the 1st and 3rd id in every angle
  for (int j = 0; j < n; j++) {
    if (angle[j][0] > angle[j][2]) {
      SwapInt(&angle[j][0], &angle[j][2]);
    }
  }
  // second, bubble sort angles
  for (int j = 0; j < (n - 1); j++) {
    bool swap = false;
    for (int k = 0; k < (n - j - 1); k++) {
      /* swap angles if
       * 1) first bead ids are in the wrong order or
       * 2) first bead ids are fine, but second ids are in wrong order or
       * 3) first and second are fine, but third are in wrong order
       */
      if ((angle[k][0] > angle[k + 1][0]) || // 1)
          (angle[k][0] == angle[k + 1][0] && // 2)
           angle[k][1] > angle[k + 1][1]) || //
          (angle[k][0] == angle[k + 1][0] && // 3)
           angle[k][1] == angle[k + 1][1] && // 3)
           angle[k][2] > angle[k + 1][2])) { // 3)
        SwapInt(&angle[k][0], &angle[k + 1][0]);
        SwapInt(&angle[k][1], &angle[k + 1][1]);
        SwapInt(&angle[k][2], &angle[k + 1][2]);
        SwapInt(&angle[k][3], &angle[k + 1][3]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}
void SortDihImp(int (*dihimp)[5], int n) { //{{{
  // first, check order of the 1st and 4th id in every dihedral
  for (int j = 0; j < n; j++) {
    if (dihimp[j][0] > dihimp[j][3]) {
      SwapInt(&dihimp[j][0], &dihimp[j][3]);
      SwapInt(&dihimp[j][1], &dihimp[j][2]);
    }
  }
  // second, bubble sort dihedrals
  for (int j = 0; j < (n - 1); j++) {
    bool swap = false;
    for (int k = 0; k < (n - j - 1); k++) {
      /* swap angles if
       * 1) first bead ids are in the wrong order or
       * 2) first bead ids are fine, but second ids are in wrong order or
       * 3) first and second are fine, but third are in wrong order or
       * 4) only fourth ids are in wrong order
       */
      if ((dihimp[k][0] > dihimp[k + 1][0]) || // 1)
          (dihimp[k][0] == dihimp[k + 1][0] && // 2)
           dihimp[k][1] > dihimp[k + 1][1]) || //
          (dihimp[k][0] == dihimp[k + 1][0] && // 3)
           dihimp[k][1] == dihimp[k + 1][1] && //
           dihimp[k][2] > dihimp[k + 1][2]) || //
          (dihimp[k][0] == dihimp[k + 1][0] && // 4)
           dihimp[k][1] == dihimp[k + 1][1] && //
           dihimp[k][2] == dihimp[k + 1][2] && //
           dihimp[k][3] > dihimp[k + 1][3])) { //
        SwapInt(&dihimp[k][0], &dihimp[k + 1][0]);
        SwapInt(&dihimp[k][1], &dihimp[k + 1][1]);
        SwapInt(&dihimp[k][2], &dihimp[k + 1][2]);
        SwapInt(&dihimp[k][3], &dihimp[k + 1][3]);
        SwapInt(&dihimp[k][4], &dihimp[k + 1][4]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// Appends # to bead/molecule types with the same name
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
          name[BEAD_NAME-2] = '\0';
        } else if (count < 100) {
          name[BEAD_NAME-3] = '\0';
        } else if (count < 1000) {
          name[BEAD_NAME-4] = '\0';
        }
        if (snprintf(System->BeadType[j].Name, BEAD_NAME, "%s%d", name,
                     count) < 0) {
          ErrorSnprintf();
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
          name[MOL_NAME-2] = '\0';
        } else if (count < 100) {
          name[MOL_NAME-3] = '\0';
        } else if (count < 1000) {
          name[MOL_NAME-4] = '\0';
        }
        if (snprintf(mt_j->Name, MOL_NAME, "%s%d", name, count) < 0) {
          ErrorSnprintf();
        }
      }
    }
  }
} //}}}
bool SameBeadType(BEADTYPE bt_1, BEADTYPE bt_2) { //{{{
  if ((strcmp(bt_1.Name, bt_2.Name) == 0) &&
      (bt_1.Charge == bt_2.Charge ||
       bt_1.Charge == HIGHNUM || bt_2.Charge == HIGHNUM) &&
      (bt_1.Mass == bt_2.Mass ||
       bt_1.Mass == HIGHNUM || bt_2.Mass == HIGHNUM) &&
      (bt_1.Radius == bt_2.Radius ||
       bt_1.Radius == HIGHNUM || bt_2.Radius == HIGHNUM)) {
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
  snprintf((*MoleculeType)[mtype].Name, MOL_NAME, "%s", name);
  // initialize struct members
  (*MoleculeType)[mtype].Number = 1;
  (*MoleculeType)[mtype].nBeads = n_beads;
  (*MoleculeType)[mtype].Bead = calloc(n_beads,
                                       sizeof *(*MoleculeType)[mtype].Bead);
  (*MoleculeType)[mtype].nBonds = n_bonds;
  if (n_bonds > 0) {
    (*MoleculeType)[mtype].Bond = calloc(n_bonds,
                                         sizeof *(*MoleculeType)[mtype].Bond);
  }
  (*MoleculeType)[mtype].nAngles = n_angles;
  if (n_angles > 0) {
    (*MoleculeType)[mtype].Angle = calloc(n_angles,
                                          sizeof *(*MoleculeType)[mtype].Angle);
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

// copy molecule type //{{{
MOLECULETYPE CopyMoleculeType(MOLECULETYPE mt_old) { //{{{
  MOLECULETYPE mt_new = CopyMoleculeTypeEssentials(mt_old);
  // MoleculeType[].Index array
  if (mt_new.Number > 0) {
    mt_new.Index = malloc(sizeof *mt_new.Index * mt_new.Number);
    memcpy(mt_new.Index, mt_old.Index, sizeof *mt_old.Index * mt_old.Number);
  }
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
    snprintf(ERROR_MSG, LINE, "molecule type without beads (%s%s%s); should "
             "never happen!", ErrYellow(), mt_new.Name, ErrYellow());
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
// copy System structure; assumes new unallocated SYSTEM //{{{
SYSTEM CopySystem(SYSTEM S_in) {
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
        BEADTYPE *bt_out = &S_out.BeadType[i],
                 *bt_in = &S_in.BeadType[i];
        *bt_out = *bt_in;
        if (bt_out->Number > 0) {
          bt_out->Index = malloc(bt_out->Number * sizeof *bt_out->Index);
          for (int j = 0; j < bt_out->InCoor; j++) {
            bt_out->Index[j] = bt_in->Index[j];
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
      S_out.MoleculeType =
          realloc(S_out.MoleculeType,
                  sizeof *S_out.MoleculeType * S_out.Count.MoleculeType);
      for (int i = 0; i < S_out.Count.MoleculeType; i++) {
        S_out.MoleculeType[i] = CopyMoleculeType(S_in.MoleculeType[i]);
      }
    } //}}}
    // Molecule & MoleculeCoor //{{{
    if (S_out.Count.Molecule > 0) {
      S_out.Molecule = realloc(S_out.Molecule,
                               sizeof *S_out.Molecule * S_out.Count.Molecule);
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
      // MoleculeCoor
      S_out.MoleculeCoor = realloc(S_out.MoleculeCoor,
                                   sizeof *S_out.MoleculeCoor *
                                   S_out.Count.Molecule);
      memcpy(S_out.MoleculeCoor, S_in.MoleculeCoor,
             sizeof *S_in.MoleculeCoor * S_in.Count.Molecule);
    } //}}}
    // BondType //{{{
    if (S_out.Count.BondType > 0) {
      S_out.BondType = realloc(S_out.BondType,
                               sizeof *S_out.BondType * S_out.Count.BondType);
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

// cleanse System by removing molecule/bead types with .Number=0, etc.
void PruneBondTypes(SYSTEM S_old, SYSTEM *System) { //{{{
  COUNT *Count = &System->Count, *Count_old = &S_old.Count;
  Count->BondType = 0;
  int *type_old_to_new = calloc(Count_old->BondType, sizeof *type_old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    for (int j = 0; j < mt->nBonds; j++) {
      int old_tbond = mt->Bond[j][2];
      if (old_tbond != -1) {
        bool new = true;
        PARAMS *tbond = &S_old.BondType[old_tbond];
        for (int k = 0; k < Count->BondType; k++) {
          if (tbond->a == System->BondType[k].a &&
              tbond->b == System->BondType[k].b &&
              tbond->c == System->BondType[k].c) {
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
    MOLECULETYPE *mt = &System->MoleculeType[i];
    for (int j = 0; j < mt->nBonds; j++) {
      int old_tbond = mt->Bond[j][2];
      if (old_tbond != -1) {
        mt->Bond[j][2] = type_old_to_new[old_tbond];
      } else {
        mt->Bond[j][2] = -1;
      }
    }
  }
  free(type_old_to_new);
} //}}}
void PruneAngleTypes(SYSTEM S_old, SYSTEM *System) { //{{{
  COUNT *Count = &System->Count, *Count_old = &S_old.Count;
  Count->AngleType = 0;
  int *type_old_to_new = calloc(Count_old->AngleType, sizeof *type_old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nAngles; j++) {
      int old_tangle = mt_i->Angle[j][3];
      bool new = true;
      if (old_tangle != -1) {
        PARAMS *tangle = &S_old.AngleType[old_tangle];
        for (int k = 0; k < Count->AngleType; k++) {
          if (fabs(tangle->a - System->AngleType[k].a) < 1e-5 &&
              fabs(tangle->b - System->AngleType[k].b) < 1e-5 &&
              fabs(tangle->c - System->AngleType[k].c) < 1e-5) {
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
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nAngles; j++) {
      int old_tangle = mt_i->Angle[j][3];
      if (old_tangle != -1) {
        mt_i->Angle[j][3] = type_old_to_new[old_tangle];
      } else {
        mt_i->Angle[j][3] = -1;
      }
    }
  }
  free(type_old_to_new);
} //}}}
void PruneDihedralTypes(SYSTEM S_old, SYSTEM *System) { //{{{
  COUNT *Count = &System->Count, *Count_old = &S_old.Count;
  Count->DihedralType = 0;
  int *type_old_to_new =
      calloc(Count_old->DihedralType, sizeof *type_old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nDihedrals; j++) {
      int old_tdihed = mt_i->Dihedral[j][4];
      bool new = true;
      if (old_tdihed != -1) {
        PARAMS *tdihed = &S_old.DihedralType[old_tdihed];
        for (int k = 0; k < Count->DihedralType; k++) {
          if (tdihed->a == System->DihedralType[k].a &&
              tdihed->b == System->DihedralType[k].b &&
              tdihed->c == System->DihedralType[k].c) {
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
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nDihedrals; j++) {
      int old_tdihed = mt_i->Dihedral[j][4];
      if (old_tdihed != -1) {
        mt_i->Dihedral[j][4] = type_old_to_new[old_tdihed];
      } else {
        mt_i->Dihedral[j][4] = -1;
      }
    }
  }
  free(type_old_to_new);
} //}}}
void PruneImproperTypes(SYSTEM S_old, SYSTEM *System) { //{{{
  COUNT *Count = &System->Count, *Count_old = &S_old.Count;
  Count->ImproperType = 0;
  int *type_old_to_new =
      calloc(Count_old->ImproperType, sizeof *type_old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nImpropers; j++) {
      int old_timpro = mt_i->Improper[j][4];
      bool new = true;
      if (old_timpro != -1) {
        PARAMS *timpro = &S_old.ImproperType[old_timpro];
        for (int k = 0; k < Count->ImproperType; k++) {
          if (timpro->a == System->ImproperType[k].a &&
              timpro->b == System->ImproperType[k].b &&
              timpro->c == System->ImproperType[k].c) {
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
  }
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    for (int j = 0; j < mt_i->nImpropers; j++) {
      int old_timpro = mt_i->Improper[j][4];
      if (old_timpro != -1) {
        mt_i->Improper[j][4] = type_old_to_new[old_timpro];
      } else {
        mt_i->Improper[j][4] = -1;
      }
    }
  }
  free(type_old_to_new);
} //}}}
void PruneSystem(SYSTEM *System) { //{{{
  SYSTEM S_old = CopySystem(*System);
  FreeSystem(System);
  InitSystem(System);
  COUNT *Count = &System->Count,
        *Count_old = &S_old.Count;
  System->Box = S_old.Box;
  *Count = *Count_old; // some counts will change later
  // allocate memory for bead count arrays
  System->Bead = realloc(System->Bead, sizeof *System->Bead * Count->Bead);
  System->BeadCoor = realloc(System->BeadCoor,
                             sizeof *System->BeadCoor * Count->Bead);
  if (Count->Bonded > 0) {
    System->Bonded = realloc(System->Bonded,
                             sizeof *System->Bonded * Count->Bonded);
    System->BondedCoor = realloc(System->BondedCoor,
                                 sizeof *System->BondedCoor * Count->Bonded);
  }
  if (Count->Unbonded > 0) {
    System->Unbonded = realloc(System->Unbonded,
                               sizeof *System->Unbonded * Count->Unbonded);
    System->UnbondedCoor = realloc(System->UnbondedCoor,
                                   sizeof *System->UnbondedCoor *
                                   Count->Unbonded);
  }
  // copy bond/angle/dihedral/improper types //{{{
  if (Count->BondType > 0) {
    System->BondType = realloc(System->BondType, Count->BondType *
                               sizeof *System->BondType);
    for (int i = 0; i < Count->BondType; i++) {
      System->BondType[i] = S_old.BondType[i];
    }
  }
  if (Count->AngleType > 0) {
    System->AngleType = realloc(System->AngleType, Count->AngleType *
                                sizeof *System->AngleType);
    for (int i = 0; i < Count->AngleType; i++) {
      System->AngleType[i] = S_old.AngleType[i];
    }
  }
  if (Count->DihedralType > 0) {
    System->DihedralType = realloc(System->DihedralType, Count->DihedralType *
                                   sizeof *System->DihedralType);
    for (int i = 0; i < Count->DihedralType; i++) {
      System->DihedralType[i] = S_old.DihedralType[i];
    }
  }
  if (Count->ImproperType > 0) {
    System->ImproperType = realloc(System->ImproperType, Count->ImproperType *
                                   sizeof *System->ImproperType);
    for (int i = 0; i < Count->ImproperType; i++) {
      System->ImproperType[i] = S_old.ImproperType[i];
    }
  } //}}}
  // copy Bead/Unbonded/Bonded arrays & create new BeadType array //{{{
  int count_unbonded = 0, count_bonded = 0, count_all = 0;
  // arrays for mapping old bead ids/types to new ones
  int *b_id_old_to_new = calloc(Count_old->Bead, sizeof *b_id_old_to_new),
      *bt_old_to_new = calloc(Count_old->BeadType, sizeof *bt_old_to_new);
  Count->BeadType = 0;
  // TODO: why use InTimestep instead of BeadCoor?
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
      // copy bonds to the mt_old_new molecule type //{{{
      int count = 0;
      if (mt_old->nBonds > 0) {
        mt_old_new.Bond = malloc(sizeof *mt_old_new.Bond * mt_old->nBonds);
        for (int j = 0; j < mt_old->nBonds; j++) {
          int id[2] = {b_old_to_old_new[mt_old->Bond[j][0]],
                       b_old_to_old_new[mt_old->Bond[j][1]]};
          if (id[0] != -1 && id[1] != -1) {
            mt_old_new.Bond[count][0] = id[0];
            mt_old_new.Bond[count][1] = id[1];
            mt_old_new.Bond[count][2] = mt_old->Bond[j][2];
            count++;
          }
        }
        if (count == 0) {
          free(mt_old_new.Bond);
        }
      }
      mt_old_new.nBonds = count; //}}}
      // copy angles to the mt_old_new molecule type //{{{
      count = 0;
      if (mt_old->nAngles > 0) {
        mt_old_new.Angle = malloc(sizeof *mt_old_new.Angle * mt_old->nAngles);
        for (int j = 0; j < mt_old->nAngles; j++) {
          int id1 = mt_old->Angle[j][0],
              id2 = mt_old->Angle[j][1],
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
      mt_old_new.nAngles = count; //}}}
      // copy dihedrals to the mt_old_new molecule type //{{{
      count = 0;
      if (mt_old->nDihedrals > 0) {
        mt_old_new.Dihedral =
            malloc(sizeof *mt_old_new.Dihedral * mt_old->nDihedrals);
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int id1 = mt_old->Dihedral[j][0],
              id2 = mt_old->Dihedral[j][1],
              id3 = mt_old->Dihedral[j][2],
              id4 = mt_old->Dihedral[j][3];
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
      mt_old_new.nDihedrals = count; //}}}
      // copy impropers to the mt_old_new molecule type //{{{
      count = 0;
      if (mt_old->nImpropers > 0) {
        mt_old_new.Improper =
            malloc(sizeof *mt_old_new.Improper * mt_old->nImpropers);
        for (int j = 0; j < mt_old->nImpropers; j++) {
          int id1 = mt_old->Improper[j][0],
              id2 = mt_old->Improper[j][1],
              id3 = mt_old->Improper[j][2],
              id4 = mt_old->Improper[j][3];
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
      mt_old_new.nImpropers = count; //}}}

      int new_id = Count->Molecule;
      Count->Molecule++;
      System->Molecule = realloc(System->Molecule,
                                 sizeof *System->Molecule * Count->Molecule);
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
      int new_type = FindMoleculeType(S_old, mt_old_new, *System, 3);
      // printf("AAA %d %s %d\n", mt_old_new.Bond[0][2], mt_old_new.Name, new_type);
      // if (new_type > -1) {
      //   printf("%d %d\n", mt_old_new.Bond[0][2],
      //          System->MoleculeType[new_type].Bond[0][2]);
      // }
      FreeMoleculeTypeEssentials(&mt_old_new);
      if (new_type != -1) { // yes, the molecule type is in the pruned system
        mol_new->Type = new_type;
        System->MoleculeType[new_type].Number++;
      } else { // no, it isn't; create a new one
        int new_new_type = Count->MoleculeType,
            c_bond = 0, c_angle = 0, c_dihedral = 0, c_improper = 0;
        int mtype = S_old.Molecule[i].Type, id[10];
        // count bonds in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nBonds; j++) {
          for (int aa = 0; aa < 2; aa++) {
            id[aa] = S_old.MoleculeType[mtype].Bond[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa]];
          }
          if (S_old.Bead[id[0]].InTimestep && S_old.Bead[id[1]].InTimestep) {
            c_bond++;
          }
        } //}}}
        // count angles in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nAngles; j++) {
          for (int aa = 0; aa < 3; aa++) {
            id[aa] = S_old.MoleculeType[mtype].Angle[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa]];
          }
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            c_angle++;
          }
        } //}}}
        // count dihedrals in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          for (int aa = 0; aa < 4; aa++) {
            id[aa] = S_old.MoleculeType[mtype].Dihedral[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa]];
          }
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            c_dihedral++;
          }
        } //}}}
        // count impropers in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nImpropers; j++) {
          for (int aa = 0; aa < 4; aa++) {
            id[aa] = S_old.MoleculeType[mtype].Improper[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa]];
          }
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            c_improper++;
          }
        } //}}}
        NewMolType(&System->MoleculeType, &Count->MoleculeType, mt_old->Name,
                   c_bead, c_bond, c_angle, c_dihedral, c_improper);
        System->Molecule[new_id].Type = new_new_type;
        System->Molecule[new_id].Aggregate = mol_old->Aggregate;
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
          for (int aa = 0; aa < 2; aa++) {
            id[aa+2] = S_old.MoleculeType[mtype].Bond[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa+2]];
          }
          int last = mt_old->Bond[j][2];
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
          int last = mt_old->Angle[j][3];
          for (int aa = 0; aa < 3; aa++) {
            id[aa+3] = S_old.MoleculeType[mtype].Angle[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa+3]];
          }
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
          int last = mt_old->Dihedral[j][4];
          for (int aa = 0; aa < 4; aa++) {
            id[aa+4] = S_old.MoleculeType[mtype].Dihedral[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa+4]];
          }
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
          int last = mt_old->Improper[j][4];
          for (int aa = 0; aa < 4; aa++) {
            id[aa+4] = S_old.MoleculeType[mtype].Improper[j][aa];
            id[aa] = S_old.Molecule[i].Bead[id[aa+4]];
          }
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
  MergeMoleculeTypes(System);
  AllocFillBeadTypeIndex(System);
  FillMoleculeTypeIndex(System);
  // prune bond/angle/dihedral/improper types //{{{
  if (Count_old->BondType > 0) {
    PruneBondTypes(S_old, System);
  }
  if (Count->AngleType > 0) {
    PruneAngleTypes(S_old, System);
  }
  if (Count->DihedralType > 0) {
    PruneDihedralTypes(S_old, System);
  }
  if (Count->ImproperType > 0) {
    PruneImproperTypes(S_old, System);
  } //}}}
  CountBondAngleDihedralImproper(System);
  for (int i = 0; i < Count->MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  System->MoleculeCoor = realloc(System->MoleculeCoor, Count->Molecule *
                                 sizeof *System->MoleculeCoor);
  Count->MoleculeCoor = 0;
  for (int i = 0; i < Count->Molecule; i++) {
    if (System->Molecule[i].InTimestep) {
      System->MoleculeCoor[Count->MoleculeCoor] = i;
      Count->MoleculeCoor++;
    }
  }
  FreeSystem(&S_old);
  free(b_id_old_to_new);
  free(bt_old_to_new);
} //}}}
// ConcatenateSystems() //{{{
// ...assumes S_out needs reallocating memory to accommodate S_in.
// TODO: some warning about S_in being empty
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box, bool prune) {
  COUNT Count_old = S_out->Count; // copy the original COUNT
  COUNT *Count_out = &S_out->Count,
        *Count_in = &S_in.Count;
  S_out->Box = Box;
  // BeadType //{{{
  if (Count_in->BeadType > 0) {
    Count_out->BeadType += Count_in->BeadType;
    S_out->BeadType = realloc(S_out->BeadType,
                              sizeof *S_out->BeadType * Count_out->BeadType);
    for (int i = 0; i < Count_in->BeadType; i++) {
      int new = i + Count_old.BeadType;
      BEADTYPE *bt_new = &S_out->BeadType[new];
      *bt_new = S_in.BeadType[i];
      bt_new->Index = malloc(sizeof *bt_new->Index * bt_new->Number);
      for (int j = 0; j < bt_new->Number; j++) {
        bt_new->Index[j] = S_in.BeadType[i].Index[j] + Count_old.Bead;
      }
    }
  } else {
    strcpy(ERROR_MSG, "no bead types to add to the system");
    PrintWarning();
    return;
  } //}}}
  // Bead & BeadCoor //{{{
  if (Count_in->Bead > 0) {
    Count_out->Bead += Count_in->Bead;
    S_out->Bead = realloc(S_out->Bead, sizeof *S_out->Bead * Count_out->Bead);
    for (int i = 0; i < Count_in->Bead; i++) {
      int new = i + Count_old.Bead;
      S_out->Bead[new] = S_in.Bead[i];
      S_out->Bead[new].Type += Count_old.BeadType;
      if (S_out->Bead[new].Molecule != -1) {
        S_out->Bead[new].Molecule += Count_old.Molecule;
      }
    }
    Count_out->BeadCoor += Count_in->BeadCoor;
    S_out->BeadCoor = realloc(S_out->BeadCoor,
                              sizeof *S_out->BeadCoor * Count_out->Bead);
    for (int i = 0; i < Count_in->BeadCoor; i++) {
      int new = i + Count_old.BeadCoor;
      S_out->BeadCoor[new] = S_in.BeadCoor[i] + Count_old.Bead;
    }
  } else {
    strcpy(ERROR_MSG, "no beads to add to the system");
    PrintWarning();
    return;
  } //}}}
  // Bonded & BondedCoor //{{{
  if (Count_in->Bonded > 0) {
    Count_out->Bonded += Count_in->Bonded;
    S_out->Bonded =
        realloc(S_out->Bonded, sizeof *S_out->Bonded * Count_out->Bonded);
    for (int i = 0; i < Count_in->Bonded; i++) {
      int new = i + Count_old.Bonded;
      S_out->Bonded[new] = S_in.Bonded[i] + Count_old.Bead;
    }
    Count_out->BondedCoor += Count_in->BondedCoor;
    S_out->BondedCoor = realloc(S_out->BondedCoor,
                                sizeof *S_out->BondedCoor * Count_out->Bonded);
    for (int i = 0; i < Count_in->BondedCoor; i++) {
      int new = i + Count_old.BondedCoor;
      S_out->BondedCoor[new] = S_in.BondedCoor[i] + Count_old.Bead;
    }
  } //}}}
  // Unbonded & UnbondedCoor //{{{
  if (Count_in->Unbonded > 0) {
    Count_out->Unbonded += Count_in->Unbonded;
    S_out->Unbonded =
        realloc(S_out->Unbonded, sizeof *S_out->Unbonded * Count_out->Unbonded);
    for (int i = 0; i < Count_in->Unbonded; i++) {
      int new = i + Count_old.Unbonded;
      S_out->Unbonded[new] = S_in.Unbonded[i] + Count_old.Bead;
    }
    Count_out->UnbondedCoor += Count_in->UnbondedCoor;
    S_out->UnbondedCoor = realloc(S_out->UnbondedCoor, Count_out->UnbondedCoor *
                                  sizeof *S_out->UnbondedCoor);
    for (int i = 0; i < Count_in->UnbondedCoor; i++) {
      int new = i + Count_old.UnbondedCoor;
      S_out->UnbondedCoor[new] = S_in.UnbondedCoor[i] + Count_old.Bead;
    }
  } //}}}
  // MoleculeType //{{{
  if (Count_in->MoleculeType > 0) {
    Count_out->MoleculeType += Count_in->MoleculeType;
    S_out->MoleculeType = realloc(S_out->MoleculeType,
                                  sizeof *S_out->MoleculeType *
                                  Count_out->MoleculeType);
    for (int i = 0; i < Count_in->MoleculeType; i++) {
      int new = i + Count_old.MoleculeType;
      MOLECULETYPE *mt_out = &S_out->MoleculeType[new],
                   *mt_in = &S_in.MoleculeType[i];
      *mt_out = CopyMoleculeType(*mt_in);
      for (int j = 0; j < mt_out->nBeads; j++) {
        mt_out->Bead[j] = mt_in->Bead[j] + Count_old.BeadType;
      }
      if (mt_out->nBonds > 0) {
        for (int j = 0; j < mt_out->nBonds; j++) {
          mt_out->Bond[j][2] += Count_old.BondType;
        }
      }
      if (mt_out->nAngles > 0) {
        for (int j = 0; j < mt_out->nAngles; j++) {
          mt_out->Angle[j][3] += Count_old.AngleType;
        }
      }
      if (mt_out->nDihedrals > 0) {
        for (int j = 0; j < mt_out->nDihedrals; j++) {
          mt_out->Dihedral[j][4] += Count_old.DihedralType;
        }
      }
      if (mt_out->nImpropers > 0) {
        for (int j = 0; j < mt_out->nImpropers; j++) {
          mt_out->Improper[j][4] += Count_old.ImproperType;
        }
      }
      for (int j = 0; j < mt_out->nBTypes; j++) {
        mt_out->BType[j] = mt_in->BType[j] + Count_old.BeadType;
      }
      // MoleculeType[].Index
      for (int j = 0; j < mt_out->Number; j++) {
        mt_out->Index[j] = mt_in->Index[j] + Count_old.Molecule;
      }
    }
  } //}}}
  // Molecule & MoleculeCoor //{{{
  if (Count_in->Molecule > 0) {
    Count_out->Molecule += Count_in->Molecule;
    S_out->Molecule = realloc(S_out->Molecule,
                              sizeof *S_out->Molecule * Count_out->Molecule);
    for (int i = 0; i < Count_in->Molecule; i++) {
      MOLECULE *mol_out = &S_out->Molecule[i+Count_old.Molecule],
               *mol_in = &S_in.Molecule[i];
      int type = mol_in->Type + Count_old.MoleculeType;
      mol_out->Type = type;
      // destroys info about S_in molecules' resids, but who cares
      mol_out->Index = Count_old.HighestResid + i + 1;
      mol_out->Bead = malloc(S_out->MoleculeType[type].nBeads *
                             sizeof *mol_out->Bead);
      for (int j = 0; j < S_out->MoleculeType[type].nBeads; j++) {
        mol_out->Bead[j] = mol_in->Bead[j] + Count_old.Bead;
      }
    }
    Count_out->HighestResid += Count_in->Molecule;
    S_out->MoleculeCoor = realloc(S_out->MoleculeCoor, Count_out->Molecule *
                                  sizeof *S_out->MoleculeCoor);
    for (int i = 0; i <= Count_in->MoleculeCoor; i++) {
      int new = i + Count_old.MoleculeCoor;
      S_out->MoleculeCoor[new] = S_in.MoleculeCoor[i] + Count_old.Molecule;
    }
  } //}}}
  // BondType //{{{
  if (Count_in->BondType > 0) {
    Count_out->BondType += Count_in->BondType;
    S_out->BondType = realloc(S_out->BondType,
                              Count_out->BondType * sizeof *S_out->BondType);
    for (int i = Count_old.BondType; i < Count_out->BondType; i++) {
      S_out->BondType[i] = S_in.BondType[i - Count_old.BondType];
    }
  } //}}}
  // AngleType //{{{
  if (Count_in->AngleType > 0) {
    Count_out->AngleType += Count_in->AngleType;
    S_out->AngleType = realloc(S_out->AngleType,
                               sizeof *S_out->AngleType * Count_out->AngleType);
    for (int i = Count_old.AngleType; i < Count_out->AngleType; i++) {
      S_out->AngleType[i] = S_in.AngleType[i - Count_old.AngleType];
    }
  } //}}}
  // DihedralType //{{{
  if (Count_in->DihedralType > 0) {
    Count_out->DihedralType += Count_in->DihedralType;
    S_out->DihedralType = realloc(S_out->DihedralType, Count_out->DihedralType *
                                  sizeof *S_out->DihedralType);
    for (int i = Count_old.DihedralType; i < Count_out->DihedralType; i++) {
      S_out->DihedralType[i] = S_in.DihedralType[i - Count_old.DihedralType];
    }
  } //}}}
  // ImproperType //{{{
  if (Count_in->ImproperType > 0) {
    Count_out->ImproperType += Count_in->ImproperType;
    S_out->ImproperType = realloc(S_out->ImproperType, Count_out->ImproperType *
                                  sizeof *S_out->ImproperType);
    for (int i = Count_old.ImproperType; i < Count_out->ImproperType; i++) {
      S_out->ImproperType[i] = S_in.ImproperType[i-Count_old.ImproperType];
    }
  } //}}}
  if (prune) {
    PruneSystem(S_out);
  }
} //}}}

// TODO: split CheckSystem to CheckCount, CheckBeadType, etc.
// check that the System struct doesn't contain an error //{{{
void CheckSystem(SYSTEM System, char file[]) {
  COUNT *Count = &System.Count;
  if (Count->Molecule > 0 &&
      Count->Molecule > (Count->HighestResid + 1)) { //{{{
    strcpy(ERROR_MSG, "Count.HighestResid is lower than Count.Molecule");
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
  // for (int i = 0; i < Count->HighestResid; i++) {
  //   if (System.Index_mol[i] < -1 || System.Index_mol[i] >= Count->Molecule) {
  //     strcpy(ERROR_MSG, "incorrect index in Index_mol array");
  //     PrintErrorFile(file, "\0", "\0");
  //     fprintf(stderr, "%s, Index_mol[%s%d%s] = %s%d%s\n", ErrRed(), ErrYellow(),
  //             i, ErrRed(), ErrYellow(), System.Index_mol[i], ErrColourReset());
  //     break;
  //   }
  // } //}}}
  // BeadType[].Index array //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System.BeadType[i];
    for (int j = 0; j < bt_i->InCoor; j++) {
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
    MOLECULETYPE *mt = &System.MoleculeType[i];
    // Index array //{{{
    for (int j = 0; j < mt->Number; j++) {
      if (mt->Index[j] < 0 || mt->Index[j] >= Count->Molecule) {
        strcpy(ERROR_MSG, "incorrect index in MoleculeType[].Index array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Index[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt->Index[j], ErrColourReset());
        break;
      }
    } //}}}
    // Bead array //{{{
    for (int j = 0; j < mt->nBeads; j++) {
      if (mt->Bead[j] < 0 || mt->Bead[j] >= Count->BeadType) {
        strcpy(ERROR_MSG, "incorrect index in Bead array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bead[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt->Bead[j], ErrColourReset());
        break;
      }
    } //}}}
    // Bond array //{{{
    for (int j = 0; j < mt->nBonds; j++) {
      if (mt->Bond[j][0] < 0 || mt->Bond[j][0] >= mt->nBeads ||
          mt->Bond[j][1] < 0 || mt->Bond[j][1] >= mt->nBeads ||
          mt->Bond[j][0] == mt->Bond[j][1]) {
        strcpy(ERROR_MSG, "incorrect index in Bond array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bond[%s%d%s][0..1]", ErrRed(),
                ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed());
        fprintf(stderr, " = %s%d %d%s\n", ErrYellow(), mt->Bond[j][0],
                mt->Bond[j][1], ErrColourReset());
        break;
      }
      if (mt->Bond[j][2] < -1 || mt->Bond[j][2] >= Count->BondType) {
        strcpy(ERROR_MSG, "incorrect bond type in Bond array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Bond[2] = %s%d%s\n", ErrRed(),
                ErrYellow(), i, ErrRed(), ErrYellow(), mt->Bond[j][2],
                ErrColourReset());
        break;
      }
    } //}}}
    // Angle array //{{{
    for (int j = 0; j < mt->nAngles; j++) {
      if (mt->Angle[j][0] < 0 || mt->Angle[j][0] >= mt->nBeads ||
          mt->Angle[j][1] < 0 || mt->Angle[j][1] >= mt->nBeads ||
          mt->Angle[j][2] < 0 || mt->Angle[j][2] >= mt->nBeads ||
          mt->Angle[j][0] == mt->Angle[j][1] ||
          mt->Angle[j][0] == mt->Angle[j][2] ||
          mt->Angle[j][1] == mt->Angle[j][2]) {
        strcpy(ERROR_MSG, "incorrect index in Angle array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Angle[0..2] = %s%d %d %d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt->Angle[j][0], mt->Angle[j][1], mt->Angle[j][2],
                ErrColourReset());
        break;
      }
      if (mt->Angle[j][3] < -1 || mt->Angle[j][3] >= Count->AngleType) {
        strcpy(ERROR_MSG, "incorrect angle type in Angle array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Angle[3] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt->Angle[j][3], ErrColourReset());
        break;
      }
    } //}}}
    // Dihedral array //{{{
    for (int j = 0; j < mt->nDihedrals; j++) {
      if (mt->Dihedral[j][0] < 0 || mt->Dihedral[j][0] >= mt->nBeads ||
          mt->Dihedral[j][1] < 0 || mt->Dihedral[j][1] >= mt->nBeads ||
          mt->Dihedral[j][2] < 0 || mt->Dihedral[j][2] >= mt->nBeads ||
          mt->Dihedral[j][3] < 0 || mt->Dihedral[j][3] >= mt->nBeads ||
          mt->Dihedral[j][0] == mt->Dihedral[j][1] ||
          mt->Dihedral[j][0] == mt->Dihedral[j][2] ||
          mt->Dihedral[j][0] == mt->Dihedral[j][3] ||
          mt->Dihedral[j][1] == mt->Dihedral[j][2] ||
          mt->Dihedral[j][1] == mt->Dihedral[j][3] ||
          mt->Dihedral[j][2] == mt->Dihedral[j][3]) {
        strcpy(ERROR_MSG, "incorrect index in Dihedral array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Dihedral[0..3] = %s",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow());
        fprintf(stderr, "%d %d %d %d%s\n", mt->Dihedral[j][0],
                mt->Dihedral[j][1], mt->Dihedral[j][2],
                mt->Dihedral[j][3], ErrColourReset());
        break;
      }
      if (mt->Dihedral[j][4] < -1 ||
          mt->Dihedral[j][4] >= Count->DihedralType) {
        strcpy(ERROR_MSG, "incorrect dihedral type in Dihedral array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Dihedral[4] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt->Dihedral[j][4], ErrColourReset());
        break;
      }
    } //}}}
    // Improper array //{{{
    for (int j = 0; j < mt->nImpropers; j++) {
      if (mt->Improper[j][0] < 0 || mt->Improper[j][0] >= mt->nBeads ||
          mt->Improper[j][1] < 0 || mt->Improper[j][1] >= mt->nBeads ||
          mt->Improper[j][2] < 0 || mt->Improper[j][2] >= mt->nBeads ||
          mt->Improper[j][3] < 0 || mt->Improper[j][3] >= mt->nBeads ||
          mt->Improper[j][0] == mt->Improper[j][1] ||
          mt->Improper[j][0] == mt->Improper[j][2] ||
          mt->Improper[j][0] == mt->Improper[j][3] ||
          mt->Improper[j][1] == mt->Improper[j][2] ||
          mt->Improper[j][1] == mt->Improper[j][3] ||
          mt->Improper[j][2] == mt->Improper[j][3]) {
        strcpy(ERROR_MSG, "incorrect index in Improper array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Improper[0..3] = %s",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow());
        fprintf(stderr, "%d %d %d %d%s\n", mt->Improper[j][0],
                mt->Improper[j][1], mt->Improper[j][2],
                mt->Improper[j][3], ErrColourReset());
        break;
      }
      if (mt->Improper[j][4] < -1 ||
          mt->Improper[j][4] >= Count->ImproperType) {
        strcpy(ERROR_MSG, "incorrect improper type in Improper array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].Improper[4] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(),
                mt->Improper[j][4], ErrColourReset());
        break;
      }
    } //}}}
    // BType array //{{{
    for (int j = 0; j < mt->nBTypes; j++) {
      if (mt->BType[j] < 0 || mt->BType[j] >= Count->BeadType) {
        strcpy(ERROR_MSG, "incorrect index in BType array");
        PrintErrorFile(file, "\0", "\0");
        fprintf(stderr, "%s, MoleculeType[%s%d%s].BType[%s%d%s] = %s%d%s\n",
                ErrRed(), ErrYellow(), i, ErrRed(), ErrYellow(), j, ErrRed(),
                ErrYellow(), mt->BType[j], ErrColourReset());
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
    // TODO: use snprintf to create ERROR_MSG
    if (test[mol_i->Index] > -1) {
      strcpy(ERROR_MSG, "same molecule index with multiple molecules");
      PrintErrorFile(file, "\0", "\0");
      fprintf(stderr, "%s, index %s%d%s with molecules %s%d %d%s\n", ErrRed(),
              ErrYellow(), mol_i->Index, ErrRed(), ErrYellow(),
              test[mol_i->Index], i, ErrColourReset());
      break;
    }
    test[mol_i->Index] = i;
  }
  free(test); //}}}
} //}}}

// simplify system for vtf output - remove stuff vtf does not support //{{{
void VtfSystem(SYSTEM *System) {
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    // remove angles
    if (System->MoleculeType[i].nAngles > 0) {
      System->MoleculeType[i].nAngles = 0;
      free(System->MoleculeType[i].Angle);
    }
    // remove dihedrals
    if (System->MoleculeType[i].nDihedrals > 0) {
      System->MoleculeType[i].nDihedrals = 0;
      free(System->MoleculeType[i].Dihedral);
    }
    // remove impropers
    if (System->MoleculeType[i].nImpropers > 0) {
      System->MoleculeType[i].nImpropers = 0;
      free(System->MoleculeType[i].Improper);
    }
    // remove bond types
    for (int j = 0; j < System->MoleculeType[i].nBonds; j++) {
      System->MoleculeType[i].Bond[j][2] = -1;
    }
    // remove angle types
    for (int j = 0; j < System->MoleculeType[i].nAngles; j++) {
      System->MoleculeType[i].Angle[j][3] = -1;
    }
    // remove dihedral types
    for (int j = 0; j < System->MoleculeType[i].nDihedrals; j++) {
      System->MoleculeType[i].Dihedral[j][4] = -1;
    }
    // remove improper types
    for (int j = 0; j < System->MoleculeType[i].nImpropers; j++) {
      System->MoleculeType[i].Improper[j][4] = -1;
    }
  }
} //}}}

// FinishSystem() //{{{
void FinishSystem(SYSTEM *System) {
  FillMoleculeTypeIndex(System);
  AllocFillBeadTypeIndex(System);
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  FillBondedUnbonded(System);
  CountBondAngleDihedralImproper(System);
  FillInCoor(System);
  if (System->Box.transform[0][1] != 0 ||
      System->Box.transform[0][2] != 0 ||
      System->Box.transform[1][2] != 0) {
    CalculateBoxData(&System->Box, 1);
  } else {
    CalculateBoxData(&System->Box, 0);
  }
  CheckSystem(*System, "\0");
} //}}}

// ChangeMolecules() //{{{
/*
 * Function to add bonds, angles, dihedrals, and/or impropers as well as their
 * types (creating new ones) to a molecule type, possibly also changing the
 * bead types in MoleculeType[].Bead array for new ones.
 */
void ChangeMolecules(SYSTEM *S_orig, SYSTEM S_add, bool name) {
  COUNT *C_orig = &S_orig->Count, *C_add = &S_add.Count,
        count_old = *C_orig;
  // add bond/angle/dihedral/improper types from Sys_add to Sys_orig //{{{
  if (C_add->BondType > 0) {
    C_orig->BondType += C_add->BondType;
    S_orig->BondType = realloc(S_orig->BondType,
                               sizeof *S_orig->BondType * C_orig->BondType);
    memcpy(S_orig->BondType + count_old.BondType, S_add.BondType,
           sizeof *S_orig->BondType * C_add->BondType);
  }
  if (C_add->AngleType > 0) {
    C_orig->AngleType += C_add->AngleType;
    S_orig->AngleType = realloc(S_orig->AngleType,
                                sizeof *S_orig->AngleType * C_orig->AngleType);
    memcpy(S_orig->AngleType + count_old.AngleType, S_add.AngleType,
           sizeof *S_orig->AngleType * C_add->AngleType);
  }
  if (C_add->DihedralType > 0) {
    C_orig->DihedralType += C_add->DihedralType;
    S_orig->DihedralType = realloc(S_orig->DihedralType,
                                   sizeof *S_orig->DihedralType *
                                   C_orig->DihedralType);
    memcpy(S_orig->DihedralType + count_old.DihedralType, S_add.DihedralType,
           sizeof *S_orig->DihedralType * C_add->DihedralType);
  }
  if (C_add->ImproperType > 0) {
    C_orig->ImproperType += C_add->ImproperType;
    S_orig->ImproperType = realloc(S_orig->ImproperType,
                                   sizeof *S_orig->ImproperType *
                                   C_orig->ImproperType);
    memcpy(S_orig->ImproperType + count_old.ImproperType, S_add.ImproperType,
           sizeof *S_orig->ImproperType * C_add->ImproperType);
  }                                                    //}}}
  for (int i = 0; i < C_orig->MoleculeType; i++) { //{{{
    MOLECULETYPE *mt_orig = &S_orig->MoleculeType[i];
    int type = FindMoleculeType(*S_orig, S_orig->MoleculeType[i], S_add, 2);
    if (type != -1) {
      MOLECULETYPE *mt_add = &S_add.MoleculeType[type];
      // // add name should the original molecule be unnamed
      // TODO: this is useless now, right?
      // if (strcmp(mt_orig->Name, NON) == 0) {
      //   strcpy(mt_orig->Name, mt_add->Name);
      // }
      // add bonds, if there are none in the original molecule type... //{{{
      if (mt_add->nBonds > 0 && mt_orig->nBonds == 0) {
        mt_orig->nBonds = mt_add->nBonds;
        mt_orig->Bond = malloc(sizeof *mt_orig->Bond * mt_orig->nBonds);
        memcpy(mt_orig->Bond, mt_add->Bond,
               sizeof *mt_add->Bond * mt_add->nBonds); //}}}
        // ...or just add bond types where missing //{{{
      } else if (C_add->BondType > 0) {
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
      } else if (C_add->AngleType > 0) {
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
      } else if (C_add->DihedralType > 0) {
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
      } else if (C_add->ImproperType > 0) {
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
  CountBondAngleDihedralImproper(S_orig);
  PruneSystem(S_orig);
} //}}}

// Add/subtract Box.Low to coordinates //{{{
void AddLow(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    int id = System->BeadCoor[i];
    for (int dd = 0; dd < 3; dd++) {
      System->Bead[id].Position[dd] += System->Box.Low[dd];
    }
  }
}
void SubtractLow(SYSTEM *System) {
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    int id = System->BeadCoor[i];
    for (int dd = 0; dd < 3; dd++) {
      System->Bead[id].Position[dd] -= System->Box.Low[dd];
    }
  }
} //}}}

void SortAggStruct(AGGREGATE *Aggregate, SYSTEM System) { //{{{
  COUNT *Count = &System.Count;
  for (int i = 0; i < (Count->Aggregate - 1); i++) {
    bool done = true;
    for (int j = 0; j < (Count->Aggregate - i - 1); j++) {
      if (Aggregate[j].Molecule[0] > Aggregate[j+1].Molecule[0]) {
        SwapInt(&Aggregate[j].nMolecules, &Aggregate[j+1].nMolecules);
        // switch the whole Aggregate[].Molecule array
        int mols; // number of molecules in the larger aggregate
        if (Aggregate[j].nMolecules > Aggregate[j+1].nMolecules) {
          mols = Aggregate[j].nMolecules;
        } else {
          mols = Aggregate[j+1].nMolecules;
        }
        for (int k = 0; k < mols; k++) {
          SwapInt(&Aggregate[j].Molecule[k], &Aggregate[j+1].Molecule[k]);
        }
        // switch bonded beads array
        SwapInt(&Aggregate[j].nBeads, &Aggregate[j+1].nBeads);
        int beads; // number of beads in the larger aggregate
        if (Aggregate[j].nBeads > Aggregate[j+1].nBeads) {
          beads = Aggregate[j].nBeads;
        } else {
          beads = Aggregate[j+1].nBeads;
        }
        for (int k = 0; k < beads; k++) {
          SwapInt(&Aggregate[j].Bead[k], &Aggregate[j+1].Bead[k]);
        }
        done = false;
      }
    }
    if (done)
      break;
  }
} //}}}

// memory-freeing functions
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
