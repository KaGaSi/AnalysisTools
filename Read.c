#include "Read.h"
// TODO test output of snprintf() to get rid of a warning; see
//      https://stackoverflow.com/questions/51534284/how-to-circumvent-format-truncation-warning-in-gcc

// RemoveExtraTypes()  //{{{
/*
 * Remove bead and molecule types with Number=0. It assumes the allocated
 * memory for BeadType and MoleculeType arrays of structures correspond to the
 * number of beads and molecules, respectively (i.e., not to the number of
 * types).
 */
void RemoveExtraTypes(SYSTEM *System) {
  COUNT *Count = &System->Count;
  // BeadType & Bead[].Type
  int count = 0;
  int *bt_old_to_new = malloc(sizeof *bt_old_to_new * Count->BeadType);
  for (int i = 0; i < Count->BeadType; i++) {
    if (System->BeadType[i].Number != 0) {
      int bt_id = count;
      count++;
      if (bt_id != i) {
        System->BeadType[bt_id] = System->BeadType[i];
      }
      bt_old_to_new[i] = bt_id;
    }
  }
  Count->BeadType = count;
  for (int i = 0; i < Count->Bead; i++) {
    int old_type = System->Bead[i].Type;
    System->Bead[i].Type = bt_old_to_new[old_type];
  }
  // sync MoleculeType[].Bead (i.e., bead types) with the new BeadType
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (System->MoleculeType[i].Number != 0) {
      for (int j = 0; j < System->MoleculeType[i].nBeads; j++) {
        int id = System->MoleculeType[i].Bead[j];
        System->MoleculeType[i].Bead[j] = bt_old_to_new[id];
      }
    }
  }
  free(bt_old_to_new);
  // MoleculeType & Molecule // TODO angles etc.
  count = 0;
  Count->Molecule = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    if (mt_i->Number != 0) {
      Count->Molecule += mt_i->Number;
      int mt_id = count;
      count++;
      if (mt_id != i) {
        // MoleculeType struct
        MOLECULETYPE *mt_new = &System->MoleculeType[mt_id];
        *mt_new = *mt_i;
        mt_new->Bead = malloc(sizeof *mt_new->Bead * mt_new->nBeads);
        memcpy(mt_new->Bead, mt_i->Bead, sizeof *mt_i->Bead * mt_i->nBeads);
        free(mt_i->Bead);
        if (mt_new->nBonds > 0) {
          mt_new->Bond = malloc(sizeof *mt_new->Bond * mt_new->nBonds);
          memcpy(mt_new->Bond, mt_i->Bond, sizeof *mt_i->Bond * mt_i->nBonds);
          free(mt_i->Bond);
        }
        if (mt_new->nAngles > 0) {
          mt_new->Angle = malloc(sizeof *mt_new->Angle * mt_new->nAngles);
          memcpy(mt_new->Angle, mt_i->Angle,
                 sizeof *mt_i->Angle * mt_i->nAngles);
          free(mt_i->Angle);
        }
        if (mt_new->nDihedrals > 0) {
          mt_new->Dihedral = malloc(sizeof *mt_new->Dihedral *
                                    mt_new->nDihedrals);
          memcpy(mt_new->Dihedral, mt_i->Dihedral,
                 sizeof *mt_i->Dihedral * mt_i->nDihedrals);
          free(mt_i->Dihedral);
        }
        if (mt_new->nImpropers > 0) {
          mt_new->Improper = malloc(sizeof *mt_new->Improper *
                                    mt_new->nImpropers);
          memcpy(mt_new->Improper, mt_i->Improper,
                 sizeof *mt_i->Improper * mt_i->nImpropers);
          free(mt_i->Improper);
        }
        // Molecule struct
        MOLECULE *mol_new = &System->Molecule[mt_id];
        mol_new->Type = mt_id;
        mol_new->Index = i;
        mol_new->Bead = malloc(sizeof *mol_new->Bead * mt_new->nBeads);
        for (int j = 0; j < mt_new->nBeads; j++) {
          int id = System->Molecule[i].Bead[j];
          System->Molecule[mt_id].Bead[j] = id;
          System->Bead[id].Molecule = mt_id;
        }
        free(System->Molecule[i].Bead);
      }
    }
  }
  Count->MoleculeType = count;
} //}}}
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
 *     atoms 0 and 1 be when the mass/charge are different to that of the
 *     last atom?
 *     TODO is that right? shouldn't all those beads below have r=1?
 *     Similarly in a simpler case: e.g., lines
 *       atom 0 n x q 1 m 1
 *       atom 1 n x     m 1 r 1
 *       atom 2 n x     m 1
 *     will make radius well defined for (with value 1) only for beads 1
 *     and 2 as the charge is different for bead 0.
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
      *old_to_new = malloc(sizeof *old_to_new * count_bt_old);
  // find the unique bead names //{{{
  int count_bnames = 0;
  char (*bname)[BEAD_NAME] = malloc(sizeof *bname);
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
  } //}}}
  if (detailed) { // use name as well as charge, mass, and radius
    // 1) find charge/mass/radius values for each unique name //{{{
    // arrays values of charge, mass, and radius for each bead type
    double diff_q[count_bnames],
           diff_m[count_bnames],
           diff_r[count_bnames];
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
    int high = 1000000;
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
        merge[i][i] = true; // 'i' and 'i' is to be merged/copied
      }
    }
    // assume same-name bead types are to be merged
    for (int i = 0; i < (Count->BeadType-1); i++) {
      for (int j = (i+1); j < Count->BeadType; j++) {
        if (strcmp(System->BeadType[i].Name, System->BeadType[j].Name) == 0) {
          merge[i][j] = true;
        }
      }
    }
    // go through each bead type and compare it to two others
    for (int i = 0; i < (Count->BeadType-1); i++) {
      BEADTYPE *bt_i = &System->BeadType[i];
      // i)
      int k = 0;
      for (; k < count_bnames; k++) {
        if (strcmp(bname[k], bt_i->Name) == 0) {
          break;
        }
      }
      // ii)
      for (int j = (i+1); j < Count->BeadType; j++) {
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
    }
    //}}}
    // 3) merge the types, counting the new types //{{{
    BEADTYPE *temp = calloc(Count->BeadType, sizeof (BEADTYPE));
    int old_bt_count = Count->BeadType, count = 0,
        *bt_older_to_old = calloc(old_bt_count, sizeof *bt_older_to_old);
    for (int i = 0; i < Count->BeadType; i++) {
      if (merge[i][i]) { // i)
        temp[count] = System->BeadType[i];
        bt_older_to_old[i] = count;
        for (int j = (i+1); j < Count->BeadType; j++) {
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
    // 4) reorder the types, placing same-named ones next to each other //{{{
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
        for (int j = (i+1); j < Count->BeadType; j++) {
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
    // 5) rename same-named bead types //{{{
    for (int i = 0; i < (Count->BeadType-1); i++) {
      count = 0;
      for (int j = (i+1); j < Count->BeadType; j++) {
        if (strcmp(System->BeadType[i].Name, System->BeadType[j].Name) == 0) {
          count++;
          char name[BEAD_NAME];
          strncpy(name, System->BeadType[j].Name, BEAD_NAME);
          // shorten name if necessary
          if (count < 10) {
            name[BEAD_NAME-3] = '\0';
          } else if (count < 100) {
            name[BEAD_NAME-4] = '\0';
          } else if (count < 1000) {
            name[BEAD_NAME-5] = '\0';
          }
          snprintf(System->BeadType[j].Name, BEAD_NAME, "%s_%d", name, count);
        }
      }
    } //}}}
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
      if (new) { // create new type...
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
      strcpy(ERROR_MSG, "something went wrong with bead type differentiation; \
this should never happen!");
      PrintWarning();
      fprintf(stderr, "%sBead count for %s%s%s type: %s%d%s and %s%d%s\n",
              ErrCyan(), ErrYellow(), System->BeadType[i].Name, ErrCyan(),
              ErrYellow(), System->BeadType[i].Number, ErrCyan(),
              ErrYellow(), count_test[i], ErrColourReset());
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
   * iv) TODO same angles, dihedrals & impropers
   */
void MergeMoleculeTypes(SYSTEM *System) {
  COUNT *Count = &System->Count;
  int count = 0,
      *old_to_new = calloc(Count->MoleculeType, sizeof *old_to_new);
  for (int i = 0; i < Count->MoleculeType; i++) { //{{{
    bool new = true;
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    int j = 0;
    for (; j < count; j++) {
      MOLECULETYPE *mt_j = &System->MoleculeType[j];
      // i) check numbers of stuff
      if (strcmp(mt_i->Name, mt_j->Name) == 0 &&
          mt_i->nBeads == mt_j->nBeads &&
          mt_i->nBonds == mt_j->nBonds &&
          mt_i->nAngles == mt_j->nAngles &&
          mt_i->nDihedrals == mt_j->nDihedrals &&
          mt_i->nImpropers == mt_j->nImpropers) {
        // ii) check bead order
        bool same_beads = true; // assume the bead order is the same
        for (int k = 0; k < mt_i->nBeads; k++) {
          if (mt_i->Bead[k] != mt_j->Bead[k]) {
            same_beads = false;
          }
        }
        // iii) check bonds
        bool same_bonds = true; // assume molecule i has j type's connectivity
        for (int k = 0; k < mt_j->nBonds; k++) {
          if (mt_i->Bond[k][0] != mt_j->Bond[k][0] ||
              mt_i->Bond[k][1] != mt_j->Bond[k][1]) {
            same_bonds = false; // nope, it doesn't; i is not type j
            break;
          }
        }
        // are molecule types i and j the same?
        if (same_beads && same_bonds) {
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
    if (new) { // create new type...
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
  MOLECULETYPE *temp = malloc(sizeof (MOLECULETYPE) * Count->MoleculeType);
  for (int i = 0; i < Count->MoleculeType; i++) {
    temp[i] = CopyMoleculeTypeEssentials(System->MoleculeType[i]);
    temp[i].Use = false;
    FreeMoleculeTypeEssentials(&System->MoleculeType[i]);
  }
  free(System->MoleculeType);
  System->MoleculeType = malloc(sizeof (MOLECULETYPE) * Count->MoleculeType);
  // array to link old molecule type indices to new ones
  old_to_new = malloc(sizeof *old_to_new * Count->MoleculeType);
  // copy the molecule types back in a proper order
  count = 0;
  for (int i = 0; i < Count->MoleculeType; i++) {
    if (!temp[i].Use) {
      System->MoleculeType[count] = CopyMoleculeTypeEssentials(temp[i]);
      old_to_new[i] = count;
      count++;
      temp[i].Use = true;
      for (int j = (i+1); j < Count->MoleculeType; j++) {
        if(strcmp(System->MoleculeType[i].Name, temp[j].Name) == 0 &&
           !temp[j].Use) {
          System->MoleculeType[count] = CopyMoleculeTypeEssentials(temp[j]);
          old_to_new[j] = count;
          count++;
          temp[j].Use = true;
        }
      }
    }
  } //}}}
  // rename same-named molecule types & free the temporary array //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    count = 0;
    for (int j = (i+1); j < Count->MoleculeType; j++) {
      if (strcmp(System->MoleculeType[i].Name,
                 System->MoleculeType[j].Name) == 0) {
        count++;
        char name[MOL_NAME];
        strncpy(name, System->MoleculeType[j].Name, MOL_NAME);
        // shorten name if necessary
        if (count < 10) {
          name[MOL_NAME-3] = '\0';
        } else if (count < 100) {
          name[MOL_NAME-4] = '\0';
        } else if (count < 1000) {
          name[MOL_NAME-5] = '\0';
        }
        snprintf(System->MoleculeType[j].Name, MOL_NAME, "%s_%d", name, count);
      }
    }
    FreeMoleculeTypeEssentials(&temp[i]);
  }
  free(temp); //}}}
  // relabel molecules with proper molecule types; yep - again //{{{
  for (int i = 0; i < Count->Molecule; i++) {
    int old_type = System->Molecule[i].Type;
    System->Molecule[i].Type = old_to_new[old_type];
  }
  free(old_to_new); //}}}
} //}}}

// Read vtf files //{{{
// VtfReadStruct() //{{{
/*
 * Read system information from vsf/vtf structure file. It can recognize bead
 * and molecule types based either on name only or on all information (name,
 * mass, charge, and radius for bead types; bead order, bonds, angles, and
 * dihedrals for molecule types).
 */
// TODO check that there are any molecules before filling the structures
SYSTEM VtfReadStruct(char struct_file[], bool detailed) {
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  FILE *vsf = OpenFile(struct_file, "r");
  // define variables and structures and arrays //{{{
  int file_line_count = 0, // total number of lines
      count_atoms = 0, // number of 'atom <id>'
      default_atom = 0, // line number of the first 'atom default' line
      count_bonds = 0; // number of bonds
  bool warned = false; // has 'a[tom] default' line warning already been issued?
  BEADTYPE bt_def;
  InitBeadType(&bt_def);
//struct bond {
//  int index1, index2; // indices from vsf file
//} *bond = calloc(1, sizeof *bond); //}}}
  int (*bond)[2] = calloc(1, sizeof *bond);
  // read struct_file line by line, saving all atom and bond lines //{{{
  /* Do something based on to line type
   *   a) atom line: save bead information
   *      i) atom default line (the first one): save into separate bt_def struct
   *      ii) atom <id> line: save as a separate Sys.BeadType & Sys.Bead
   *      iii) resid keyword: create new Sys.MoleculeType & Sys.Molecule if
   *           this resid wasn't used yet, otherwise add to existing ones
   *   b) bond line: save bonded bead indices into a bond struct
   *   c) timestep line: break the loop (end of vtf file's structure section)
   *   d) coordinate line: exit program as coordinates cannot be inside a
   *      structure section
   *   e) anything else besides pbc, blank, or comment line: exit program as
   *      unrecognised line was encountered
   */
  char line[LINE], *split[SPL_STR];
  int words;
  while (ReadAndSplitLine(vsf, LINE, line, &words, split, SPL_STR, " \t\n")) {
    file_line_count++;
    // read line
    int ltype = VtfCheckLineType(words, split, struct_file, file_line_count);
    if (ltype == ATOM_LINE) { // a)
      if (strcmp(split[1], "default") == 0) { // 'a[tom] default' line
        if (default_atom != 0 && !warned) { // warn of multiple defaults //{{{
          warned = true; // warn only once
          strcpy(ERROR_MSG, "multiple 'a[tom] default' lines");
          PrintWarning();
          WarnPrintFile(struct_file, "\0");
          fprintf(stderr, "%s, using line %s%d%s as the default line%s\n",
                  ErrCyan(), ErrYellow(), default_atom, ErrCyan(),
                  ErrColourReset()); //}}}
        } else { // save line number of the first 'atom default' line //{{{
          default_atom = file_line_count;
          // save values for the default bead type
          int *values = VtfAtomLineValues(words, split);
          strncpy(bt_def.Name, split[values[0]], BEAD_NAME);
          if (values[1] != -1) {
            bt_def.Mass = atof(split[values[1]]);
          }
          if (values[2] != -1) {
            bt_def.Charge = atof(split[values[2]]);
          }
          if (values[3] != -1) {
            bt_def.Radius = atof(split[values[3]]);
          }
        } //}}}
      } else { // 'a[tom] <id>' line //{{{
        int id = atoi(split[1]);
        // warning - repeated atom line //{{{
        if (id < Count->Bead && Sys.BeadType[id].Number != 0) {
          strcpy(ERROR_MSG, "atom defined multiple times; \
disregarding the following line");
          PrintWarningFileLine(struct_file, "\0",
                               file_line_count, split, words);
          continue; // go to reading the next line
        } //}}}
        count_atoms++;
        // highest bead index? (corresponds to the number of beads in vsf)
        if (id >= Count->Bead) {
          Sys.BeadType = realloc(Sys.BeadType, sizeof (BEADTYPE) * (id + 1));
          Sys.Bead = realloc(Sys.Bead, sizeof (BEAD) * (id + 1));
          for (int i = Count->Bead; i <= id; i++) {
            InitBeadType(&Sys.BeadType[i]);
            InitBead(&Sys.Bead[i]);
          }
          Count->Bead = id + 1; // +1 as bead ids start from 0 in vsf
        }
        // save values from the 'a[tom] <id>' line
        int *values = VtfAtomLineValues(words, split);
        strncpy(Sys.BeadType[id].Name, split[values[0]], BEAD_NAME);
        Sys.BeadType[id].Name[BEAD_NAME-1] = '\0'; // ensure null-termination
        if (values[1] != -1) {
          Sys.BeadType[id].Mass = atof(split[values[1]]);
        }
        if (values[2] != -1) {
          Sys.BeadType[id].Charge = atof(split[values[2]]);
        }
        if (values[3] != -1) {
          Sys.BeadType[id].Radius = atof(split[values[3]]);
        }
        Sys.BeadType[id].Number = 1;
        Sys.Bead[id].Type = id;
        // is the bead in a molecule?
        if (values[5] > -1 ) {
          int resid = atoi(split[values[5]]);
          // highest molecule id?
          if (resid > Count->HighestResid) {
            Sys.MoleculeType = realloc(Sys.MoleculeType,
                                       sizeof (MOLECULETYPE) * (resid + 1));
            Sys.Molecule = realloc(Sys.Molecule, sizeof (MOLECULE) * (resid + 1));
            for (int i = (Count->HighestResid+1); i <= resid; i++) {
              InitMoleculeType(&Sys.MoleculeType[i]);
              InitMolecule(&Sys.Molecule[i]);
            }
            Count->HighestResid = resid; // goes from 0
            Count->Molecule = resid + 1;
          }
          MOLECULETYPE *mt_resid = &Sys.MoleculeType[resid];
          if (mt_resid->Number == 0) { // new molecule type
            strncpy(mt_resid->Name, split[values[4]], MOL_NAME);
            mt_resid->Name[MOL_NAME-1] = '\0'; // null-terminate!
            mt_resid->Number = 1;
            mt_resid->nBeads = 1;
            mt_resid->Bead = malloc(sizeof *mt_resid->Bead);
            mt_resid->Bead[0] = id; // bead type = bead index
          } else { // not new moleclue type
            int bead = mt_resid->nBeads;
            mt_resid->nBeads++;
            mt_resid->Bead = realloc(mt_resid->Bead,
                                     sizeof *mt_resid->Bead * mt_resid->nBeads);
            mt_resid->Bead[bead] = id; // bead type = bead index
          }
          Sys.Bead[id].Molecule = resid;
        }
      } //}}}
    } else if (ltype == BOND_LINE) { // b) //{{{
      bond = realloc(bond, sizeof *bond * (count_bonds + 1));
      long val;
      if (words == 2) { // case 'bond <id>:<id>'
        char *index[SPL_STR];
        SplitLine(SPL_STR, index, split[1], ":");
        IsInteger(index[0], &val);
        bond[count_bonds][0] = val;
        IsInteger(index[1], &val);
        bond[count_bonds][1] = val;
      } else { // case 'bond <id>: <id>'
        IsInteger(split[1], &val);
        bond[count_bonds][0] = val;
        IsInteger(split[2], &val);
        bond[count_bonds][1] = val;
      }
      // assure index1<index2 (may be unnecessary, but definitely won't hurt)
      if (bond[count_bonds][0] > bond[count_bonds][1]) {
        SwapInt(&bond[count_bonds][0], &bond[count_bonds][1]);
      }
      count_bonds++; //}}}
    } else if (ltype == TIME_LINE_I || ltype == TIME_LINE_O) { // c)
      break;
    } else if (ltype == COOR_LINE_I || ltype == COOR_LINE_O) { // d)
      strcpy(ERROR_MSG, "encountered a coordinate-like line \
inside the structure block ");
      PrintErrorFileLine(struct_file, "\0", file_line_count, split, words);
      exit(1);
    } else if (ltype != BLANK_LINE && ltype != COMMENT_LINE &&
               ltype != PBC_LINE && ltype != PBC_LINE_ANGLES) { // e)
      // proper error message already established in VtfCheckLineType()
      PrintErrorFileLine(struct_file, "\0", file_line_count, split, words);
      exit(1);
    }
  }
  fclose(vsf); //}}}
  Count->BeadType = Count->Bead;
  Count->MoleculeType = Count->Molecule;
  // error - no default line and too few atom lines //{{{
  if (default_atom == 0 && count_atoms != Count->Bead) {
    strcpy(ERROR_MSG, "not all beads defined ('atom default' line is omitted)");
    PrintError();
    ErrorPrintFile(struct_file, "\0");
    int undefined = Count->Bead - count_atoms;
    fprintf(stderr, "%s, %s%d%s bead(s) undefined%s\n", ErrRed(), ErrYellow(),
            undefined, ErrRed(), ErrColourReset());
    exit(1);
  } //}}}
  // assign atom default to default beads & count bonded/unbonded beads //{{{
  // find first unused bead type and make it the default
  int def, count_def = 0;
  for (int i = 0; i < Count->Bead; i++) {
    if (Sys.BeadType[i].Number == 0) {
     def = i;
     Sys.BeadType[def] = bt_def;
     Sys.BeadType[def].Number = Count->Bead - count_atoms;
     Sys.Bead[def].Type = def;
     count_def = 1;
     break;
    }
  }
  // add the default beads to their proper type & count Bonded/Unbonded
  for (int i = 0; i < Count->Bead; i++) {
    if (Sys.BeadType[i].Number == 0) { // default bead?
      Sys.Bead[i].Type = def;
      count_def++;
    }
    if (Sys.Bead[i].Molecule == -1) { // unbonded bead?
      Count->Unbonded++;
      Sys.Unbonded = realloc(Sys.Unbonded,
                             sizeof *Sys.Unbonded * Count->Unbonded);
      Sys.Unbonded[Count->Unbonded-1] = i;
    } else {
      Count->Bonded++; // bonded bead?
      Sys.Bonded = realloc(Sys.Bonded, sizeof *Sys.Bonded * Count->Bonded);
      Sys.Bonded[Count->Bonded-1] = i;
    }
  }
  // just check that it counts the beads correctly
  if (Count->Bead != (Count->Unbonded + Count->Bonded)) {
    strcpy(ERROR_MSG, "something went wrong with bead counting; \
contact developper\n");
    PrintError();
    exit(1);
  } //}}}
  FillMoleculeBeads(&Sys);
  FillMoleculeTypeBonds(&Sys, bond, count_bonds);
  RemoveExtraTypes(&Sys);
  MergeBeadTypes(&Sys, detailed);
  MergeMoleculeTypes(&Sys);
  FillSystemNonessentials(&Sys);
  CheckSystem(Sys, struct_file);
  WarnChargedSystem(Sys, struct_file, "\0");
  return Sys;
} //}}}
// VtfReadPBC() //{{{
/*
 * Get the first pbc line from a vcf/vtf coordinate file. If a coordinate line
 * is encountered before the pbc one, the function exits with an error.
 */
void VtfReadPBC(char input_vcf[], char input_vsf[], BOX *Box) {
  // open the coordinate file
  FILE *coor = OpenFile(input_vcf, "r");
  int file_line_count = 0;
  // read input_vcf line by line
  while (true) {
    file_line_count++;
    // read line & split it via whitespace //{{{
    char *split[SPL_STR], line[LINE];
    int words;
    if (!ReadAndSplitLine(coor, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "missing pbc line (and any coordinates)");
      PrintError();
      ErrorPrintFile(input_vcf, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    int ltype = VtfCheckLineType(words, split, input_vcf, file_line_count);
    // pbc line //{{{
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) {
      (*Box).Length.x = atof(split[1]);
      (*Box).Length.y = atof(split[2]);
      (*Box).Length.z = atof(split[3]);
      // angles in pbc line; possibly triclinic cell
      if (ltype == PBC_LINE_ANGLES) {
        (*Box).alpha = atof(split[4]);
        (*Box).beta = atof(split[5]);
        (*Box).gamma = atof(split[6]);
      } else { // definitely orthogonal box
        (*Box).alpha = 90;
        (*Box).beta = 90;
        (*Box).gamma = 90;
      }
      break; //}}}
    // error - coordinate line //{{{
    } else if (ltype == COOR_LINE_I || ltype == COOR_LINE_O) {
      strcpy(ERROR_MSG, "encountered coordinate line before pbc line");
      PrintErrorFileLine(input_vcf, "\0", file_line_count, split, words);
      exit(1); //}}}
    // warning - unrecognised line //{{{
    } else if (ltype == ERROR_LINE) {
      strcpy(ERROR_MSG, "ignoring unrecognised line while \
searching for a pbc line");
      PrintWarningFileLine(input_vcf, input_vsf, file_line_count, split, words);
    } //}}}
  };
  fclose(coor);
} //}}}
// VtfReadTimestep() //{{{
/*
 * Read a single timestep from a vcf/vtf coordinate line - first the preamble
 * (getting timestep type and possibly pbc, ignoring any unrecognised line),
 * then the coordinates (if a wrong line is in the coordinate section, the
 * function skips the remaining coordinate lines and returns to the beginning
 * to read the next timestep). Returns false only if a line couldn't be read
 * from the preamble (e.g., on eof); on eof within the coordinate block, true
 * is retuerned as some coordinates were read (i.e., a valid timestep).
 */
bool VtfReadTimestep(FILE *vcf, char vcf_file[], char vsf_file[],
                     SYSTEM *System, int *file_line_count,
                     int step_count, char stuff[]) {
  start_function: ; // return here when a bad line is encountered
  stuff[0] = '\0';
  char *cur = stuff, * const end = stuff + LINE; // to properly snprintf stuff[]
  // read timestep preamble //{{{
  bool coor_found = false;
  int ltype, timestep = ERROR_LINE;
  fpos_t position; // to save file position
  while (true) {
    // save file pointer position for when it's the first coordinate line
    fgetpos(vcf, &position);
    (*file_line_count)++;
    // read line & split it via whitespace //{{{
    char *split[SPL_STR], line[LINE];
    int words;
    if (!ReadAndSplitLine(vcf, LINE, line, &words, split, SPL_STR, " \t\n")) {
      return false;
    } //}}}
    ltype = VtfCheckLineType(words, split, vcf_file, *file_line_count);
    // do something based on what line it is
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) { //{{{
      BOX *Box = &System->Box;
      Box->Length.x = atof(split[1]);
      Box->Length.y = atof(split[2]);
      Box->Length.z = atof(split[3]);
      if (ltype == PBC_LINE_ANGLES) {
        Box->alpha = atof(split[4]);
        Box->beta = atof(split[5]);
        Box->gamma = atof(split[6]);
        if (!TriclinicCellData(&(*System).Box)) {
          ErrorPrintFull2(vcf_file, *file_line_count, split, words);
          exit(1);
        }
      } else {
        Box->alpha = 90;
        Box->beta = 90;
        Box->gamma = 90;
      }
     //}}}
    } else if (ltype == TIME_LINE_I || ltype == TIME_LINE_O) { //{{{
      // warn: multiple timestep lines
      if (timestep != ERROR_LINE) {
        strcpy(ERROR_MSG, "extra timestep line \
(or timestep without any coordinates)");
        PrintWarningFileLine(vcf_file, vsf_file, *file_line_count,
                             split, words);
      }
      timestep = ltype;
     //}}}
    } else if (ltype == COOR_LINE_I || ltype == COOR_LINE_O) { //{{{
      // warn: missing timestep line - read next timestep //{{{
      if (timestep == ERROR_LINE) {
        strcpy(ERROR_MSG, "found a coordinate line before a timestep line; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, vsf_file,
                             *file_line_count, split, words);
        // skip the remaining coordinate lines (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      coor_found = true;
      break; //}}}
    } else if (ltype == COMMENT_LINE) { //{{{
      // add the comment to the stuff array
      if ((end-cur) > 0) {
        cur += snprintf(cur, end-cur, "%s", split[0]);
      }
        for (int i = 1; i < (words-1) && cur < end; i++) {
          if ((end-cur) > 0) {
            cur += snprintf(cur, end-cur, " %s", split[i]);
          } else {
            break;
          }
        }
        if ((end-cur) > 0) {
          cur += snprintf(cur, end-cur, " %s\n", split[words-1]);
        }
    //}}}
    } else if (ltype == ERROR_LINE) { //{{{
      strcpy(ERROR_MSG, "ignoring unrecognised line in a timestep preamble");
      PrintWarningFileLine(vcf_file, vsf_file, *file_line_count, split, words);
    } //}}}
  } //}}}
  // warning - coordinate line encountered; should never trigger //{{{
  if (!coor_found) {
    strcpy(ERROR_MSG, "no coordinates; this warning should never trigger!");
    PrintWarning();
  } //}}}
  // set 'not in timestep' to all beads //{{{
  for (int i = 0; i < (*System).Count.Bead; i++) {
    System->Bead[i].InTimestep = false;
  } //}}}
  // read coordinates  //{{{
  COUNT *Count = &System->Count;
  Count->BeadCoor = 0;
  Count->UnbondedCoor = 0;
  Count->BondedCoor = 0;
  // restore file pointer to before the first coordinate line
  fsetpos(vcf, &position);
  (*file_line_count)--; // the first coordinate line will be re-read
  while (true) {
    // save file pointer position for when it's the first coordinate line
    fgetpos(vcf, &position);
    (*file_line_count)++;
    // read line & split it via whitespace //{{{
    char line[LINE], *split[SPL_STR];
    int words;
    if (!ReadAndSplitLine(vcf, LINE, line, &words, split, SPL_STR, " \t\n")) {
      return true;
    } //}}}
    ltype = VtfCheckLineType(words, split, vcf_file, *file_line_count);
    if (ltype != COOR_LINE_O && ltype != COOR_LINE_I) {
      // warn: unrecognised line - read next timestep //{{{
      if (ltype == ERROR_LINE) {
        strcpy(ERROR_MSG, "unrecognised line in a timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, vsf_file,
                             *file_line_count, split, words);
        // skip the remaining coordinate lines (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      break;
    }
    int id, indexed = 0;
    if (timestep == TIME_LINE_I) { // 'timestep indexed' coordinate line
      if (ltype == COOR_LINE_I) {
        id = atoi(split[0]);
        indexed = 1;
        // warn: bead index is too high - read next timestep //{{{
        if (id >= Count->Bead) {
          strcpy(ERROR_MSG, "bead index too high; \
using next timestep instead of this one");
          PrintWarningFileLine(vcf_file, vsf_file,
                               *file_line_count, split, words);
          // skip the remaing coordinate line (if any)
          do {
            fgetpos(vcf, &position);
            (*file_line_count)++;
          } while (VtfSkipCoorOrderedLine(vcf));
          fsetpos(vcf, &position);
          (*file_line_count)--; // the first non-coordinate line will be re-read
          goto start_function;
        } //}}}
        // warn: repeated bead index - read next timestep //{{{
        if (System->Bead[id].InTimestep) {
          strcpy(ERROR_MSG, "multiple bead entries with the same index; \
using next timestep instead of this one");
          PrintWarningFileLine(vcf_file, vsf_file,
                               *file_line_count, split, words);
          // skip the remaing coordinate line (if any)
          do {
            fgetpos(vcf, &position);
            (*file_line_count)++;
          } while (VtfSkipCoorOrderedLine(vcf));
          fsetpos(vcf, &position);
          (*file_line_count)--; // the first non-coordinate line will be re-read
          goto start_function;
        } //}}}
      } else {
        // warn: ordered line in indexed timestep - read next timestep //{{{
        strcpy(ERROR_MSG, "ordered coordinate line in indexed timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, vsf_file,
                             *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function; //}}}
      }
    } else { // 'timestep ordered' coordinate line
      // warn: extra bead in ordered timestep - read next timestep //{{{
      if (Count->BeadCoor == Count->Bead) {
        strcpy(ERROR_MSG, "too many beads in an ordered timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, vsf_file,
                             *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      id = (*System).Count.BeadCoor;
    }
    BEAD *bead_id = &System->Bead[id];
    bead_id->Position.x = atof(split[0+indexed]);
    bead_id->Position.y = atof(split[1+indexed]);
    bead_id->Position.z = atof(split[2+indexed]);
    bead_id->InTimestep = true;
    VECTOR vel;
    if (words >= 7 && IsReal(split[3+indexed], &vel.x) &&
                      IsReal(split[4+indexed], &vel.y) &&
                      IsReal(split[5+indexed], &vel.z)) {
      bead_id->Velocity.x = vel.x;
      bead_id->Velocity.y = vel.y;
      bead_id->Velocity.z = vel.z;
    }
    if (bead_id->Molecule == -1) {
      System->UnbondedCoor[Count->UnbondedCoor] = id;
      Count->UnbondedCoor++;
    } else {
      System->BondedCoor[Count->BondedCoor] = id;
      Count->BondedCoor++;
    }
    System->BeadCoor[Count->BeadCoor] = id;
    Count->BeadCoor++;
  }
  // TODO in file writing, also write vsf file!
  // restore file pointer to before the first non-coordinate line
  fsetpos(vcf, &position); //}}}
  // warn: too few beads in an ordered timestep - read next timestep //{{{
  if (timestep == TIME_LINE_O && Count->BeadCoor < Count->Bead) {
    strcpy(ERROR_MSG, "insufficient number of beads for ordered timestep; \
using next timestep instead of this one");
    PrintWarning();
    WarnPrintFile(vcf_file, vsf_file);
    fprintf(stderr, "%s, last line of the timestep: %s%d%s\n",
            ErrCyan(), ErrYellow(), *file_line_count, ErrColourReset());
    goto start_function;
  } //}}}
  (*file_line_count)--; // the last line will be re-read next time
  // error - test count of beads incorrect; should never trigger //{{{
  int count = Count->UnbondedCoor + Count->BondedCoor;
  if (count != Count->BeadCoor) {
    strcpy(ERROR_MSG, "unbonded and bonded beads in a coordinate file do not \
add up properly!");
    PrintError();
    ErrorPrintFile(vcf_file, vsf_file);
    fprintf(stderr, "%s, unbonded: %s%d%s",
            ErrRed(), ErrYellow(), Count->UnbondedCoor, ErrRed());
    fprintf(stderr, ", bonded: %s%d%s",
            ErrYellow(), Count->BondedCoor, ErrRed());
    fprintf(stderr, ", sum should be: %s%d%s\n",
            ErrYellow(), Count->BeadCoor, ErrColourReset());
  } //}}}
  return true; // coordinates read properly
} //}}}
// VtfSkipTimestep() //{{{
/*
 * Discard a single timestep from a vcf/vtf coordinate file. It assumes the
 * coordinates are ordered lines (i.e., three numbers); if a timestep line is
 * missing, it prints only a warning (and skips the ensuing coordinate lines).
 * Returns false only if a line cannot be read (e.g., on eof).
 */
bool VtfSkipTimestep(FILE *vcf, char vcf_file[], char vsf_file[],
                     int *file_line_count, int step_count) {
  // skip preamble - i.e., read until the first coordinate line
  fpos_t position;
  int ltype;
  bool timestep = false; // is timestep line present?
  do {
    char line[LINE];
    fgetpos(vcf, &position);
    (*file_line_count)++;
    if (!ReadLine(vcf, LINE, line)) {
      return false;
    }
    char *split[SPL_STR];
    int words = SplitLine(SPL_STR, split, line, " \t\n");
    ltype = VtfCheckLineType(words, split, vcf_file, *file_line_count);
    if (ltype == TIME_LINE_I || ltype == TIME_LINE_O) {
      timestep = true;
    }
  } while (ltype != COOR_LINE_O && ltype != COOR_LINE_I);
  // when timestep line is missing
  if (!timestep) {
    strcpy(ERROR_MSG, "found no timestep line in a timestp supposed \
to be skipped; skipping only this invalid timestep");
    PrintWarning();
    WarnPrintFile(vcf_file, vsf_file);
    fprintf(stderr, "%s, line %s%d%s\n", ErrCyan(), ErrYellow(),
                                         *file_line_count, ErrColourReset());
  }
  // skip coordinate lines - i.e., read until the first non-coordinate line
  do {
    fgetpos(vcf, &position);
    (*file_line_count)++;
  } while (VtfSkipCoorOrderedLine(vcf));
  (*file_line_count)--;
  // return to before the first non-coordinate line
  fsetpos(vcf, &position);
  return true;
} //}}}
// VtfSkipCoorOrderedLine() //{{{
/*
 * Helper funciton to discard an ordered coordinate line (e.i., a line starting
 * with three real numbers).
 */
bool VtfSkipCoorOrderedLine(FILE *fr) {
  char line[LINE];
  if (!ReadLine(fr, LINE, line)) {
    return false; // error/EOF
  }
  int strings = 3;
  char *split[strings];
  int words = SplitLine(strings, split, line, " \t\n");
  if (VtfCheckCoorOrderedLine(words, split) == COOR_LINE_O) {
    return true;
  } else {
    return false;
  }
} //}}}
// VtfCheckLineType() //{{{
/*
 * Function returning line type according to the codes in Read.h.
 */
int VtfCheckLineType(int words, char *split[], char file[], int line) {
  ERROR_MSG[0] = '\0'; // clear error message array
  // blank line
  if (words == 0) {
    return BLANK_LINE;
  }
  // comment line
  if (split[0][0] == '#') {
    return COMMENT_LINE;
  }
  // coordinate line
  int test = VtfCheckCoordinateLine(words, split);
  if (test != ERROR_LINE) {
    return test;
  }
  // timestep line
  test = VtfCheckTimestepLine(words, split);
  if (test != ERROR_LINE) {
    return test;
  }
  // pbc line
  test = VtfCheckPbcLine(words, split);
  if (test != ERROR_LINE) {
    return test;
  }
  // atom line (vsf)
  if (VtfCheckAtomLine(words, split)) {
    return ATOM_LINE;
  }
  // bond line (vsf)
  if (VtfCheckBondLine(words, split)) {
    return BOND_LINE;
  }
  if (ERROR_MSG[0] == '\0') {
    strcpy(ERROR_MSG, "unrecognised line");
  }
  return ERROR_LINE;
} //}}}
void FillMoleculeBeads(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    MOLECULE *mol_i = &System->Molecule[i];
    if (mt_i->Number > 0) {
      mol_i->Type = i;
      mol_i->Index = i;
      mol_i->Bead = malloc(sizeof *mol_i->Bead * mt_i->nBeads);
      memcpy(mol_i->Bead, mt_i->Bead, sizeof *mt_i->Bead * mt_i->nBeads);
    }
  }
} //}}}
void FillMoleculeTypeBonds(SYSTEM *System, int (*bond)[2], int nbonds) { //{{{
  COUNT *Count = &System->Count;
  // fill MoleculeType[].Bond array with bead indices
  for (int i = 0; i < nbonds; i++) {
    int id1 = bond[i][0],
        id2 = bond[i][1],
        mol = System->Bead[id1].Molecule;
    // warning - bonded beads in different molecules (skip the bond)  //{{{
    if (mol != System->Bead[id2].Molecule || mol == -1) {
      strcpy(ERROR_MSG, "bonded beads in different molecules (or in none); \
discarding this bond");
      PrintWarning();
      fprintf(stderr, "%sBead (molecule):", ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s);", ErrYellow(), id1, ErrCyan(),
              ErrYellow(), System->Bead[id1].Molecule, ErrCyan());
      fprintf(stderr, " %s%d%s (%s%d%s)\n", ErrYellow(), id2, ErrCyan(),
              ErrYellow(), System->Bead[id2].Molecule, ErrCyan());
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int bond = mt_mol->nBonds;
    mt_mol->nBonds++;
    if (bond == 0) {
      mt_mol->Bond = malloc(sizeof *mt_mol->Bond);
    } else {
      mt_mol->Bond = realloc(mt_mol->Bond,
                             sizeof *mt_mol->Bond * mt_mol->nBonds);
    }
    mt_mol->Bond[bond][0] = id1;
    mt_mol->Bond[bond][1] = id2;
    mt_mol->Bond[bond][2] = -1;
  }
  // make the MoleculeType[].Bond bead indices go from 0 to nBeads
  for (int i = 0; i < Count->MoleculeType; i++) {
    int lowest = 1e7;
    for (int j = 0; j < System->MoleculeType[i].nBonds; j++) {
      if (System->MoleculeType[i].Bond[j][0] < lowest) {
        lowest = System->MoleculeType[i].Bond[j][0];
      }
    }
    for (int j = 0; j < System->MoleculeType[i].nBonds; j++) {
      System->MoleculeType[i].Bond[j][0] -= lowest;
      System->MoleculeType[i].Bond[j][1] -= lowest;
      // warning - too high an intramolecular bead index; shouldn't happen //{{{
      if (System->MoleculeType[i].Bond[j][0] > System->MoleculeType[i].nBeads ||
          System->MoleculeType[i].Bond[j][0] > System->MoleculeType[i].nBeads) {
        strcpy(ERROR_MSG, "something went wrong in bead indices in bond; \
should never happen!");
        PrintWarning();
      } //}}}
    }
  }
  free(bond);
} //}}}
// Helper functions to check whether provided line is of a given type
int VtfCheckCoorOrderedLine(int words, char *split[]) { //{{{
  double val_d;
  if (words > 2 && IsReal(split[0], &val_d) &&
                   IsReal(split[1], &val_d) &&
                   IsReal(split[2], &val_d)) {
    return COOR_LINE_O;
  }
  return ERROR_LINE;
} //}}}
int VtfCheckCoorIndexedLine(int words, char *split[]) { //{{{
  long val_i;
  double val_d;
  // indexed line (may also be ordered)
  if (words > 3 && IsNatural(split[0], &val_i) &&
                   IsReal(split[1], &val_d) &&
                   IsReal(split[2], &val_d) &&
                   IsReal(split[3], &val_d)) {
    return COOR_LINE_I;
  }
  return ERROR_LINE;
} //}}}
int VtfCheckCoordinateLine(int words, char *split[]) { //{{{
  // indexed line (may also be ordered)
  if (VtfCheckCoorIndexedLine(words, split) == COOR_LINE_I) {
    return COOR_LINE_I;
  }
  // definitely ordered line
  if (VtfCheckCoorOrderedLine(words, split) == COOR_LINE_O) {
    return COOR_LINE_O;
  }
  return ERROR_LINE;
} //}}}
int VtfCheckTimestepLine(int words, char *split[]) { //{{{
  // there are several possibilities how the timestep line can look
  /* ordered timestep:
   *   1) 't[imestep]'
   *   2) 't[imestep] o[rdered] ...'
   *   3) 'o[rdered] ...'
   */
  if ((words == 1 && split[0][0] == 't') || // 1)
      (words > 1 && split[0][0] == 't' && split[1][0] == 'o') || // 2)
       split[0][0] == 'o') { // 3)
    return TIME_LINE_O;
  }
  /* indexed timestep:
   *   1) 't[imestep] i[ndexed] ...'
   *   2) 'i[ndexed] ...'
   */
  if ((words > 1 && split[0][0] == 't' && split[1][0] == 'i') || // 1)
       split[0][0] == 'i') { // 2)
    return TIME_LINE_I;
  }
  return ERROR_LINE; // not a timestep line
} //}}}
int VtfCheckPbcLine(int words, char *split[]) { //{{{
  // valid line: pbc <x> <y> <z> [<alpha> <beta> <gamm>]
  // unrecognised line
  double val;
  if (words < 4 || strcmp(split[0], "pbc") != 0 ||
      !IsPosReal(split[1], &val) ||
      !IsPosReal(split[2], &val) ||
      !IsPosReal(split[3], &val)) {
    return ERROR_LINE;
  } else if (words > 6 && IsPosReal(split[4], &val) &&
                          IsPosReal(split[5], &val) &&
                          IsPosReal(split[6], &val)) {
    return PBC_LINE_ANGLES;
  } else {
    return PBC_LINE;
  }
} //}}}
bool VtfCheckAtomLine(int words, char *split[]) { //{{{
  long val_i;
  double val_d;
  // error - line not starting with a[tom] default/<id> //{{{
  if (split[0][0] != 'a' ||
      (strcmp(split[1], "default") != 0 && !IsInteger(split[1], &val_i))) {
    return false;
  } //}}}
  // error - odd number of strings //{{{
  if ((words%2) != 0) {
    strcpy(ERROR_MSG, "atom line with odd number of strings ");
    return false;
  } //}}}
  // check <keyword> <value> pairs
  bool name = false, resid = false, resname = false;
  for (int i = 2; i < words; i+=2) {
    // is n[ame] keyword present?
    if (split[i][0] == 'n') {
      name = true;
    }
    int r_id = strcmp(split[i], "resid"); // resid cannot be shortened
    int r_name = strncmp(split[i], "res", 3); // res[name] can be shortened
    // is resid keyword present? //{{{
    if (r_id == 0) {
      resid = true;
      // resid must be followed by non-negative integer
      if (!IsNatural(split[i+1], &val_i)) {
        strcpy(ERROR_MSG, "atom line: 'resid' not followed by natural number");
        return false;
      }
    } //}}}
    // is res[name] keyword present? //{{{
    if (r_id != 0 && r_name == 0) {
      resname = true;
    } //}}}
    // error - charge|q //{{{
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') &&
        !IsReal(split[i+1], &val_d)) {
      strcpy(ERROR_MSG, "atom line: 'charge|q' not followed by real number ");
      return false; //}}}
    // error - r[adius] not followed by positive number //{{{
    } else if (split[i][0] == 'r' && r_name != 0 &&
               !IsPosReal(split[i+1], &val_d)) {
      strcpy(ERROR_MSG, "atom line: 'r[adius]]]' not followed by \
positive real number ");
      return false; //}}}
    // error - m[ass] not followed by positive number //{{{
    } else if (split[i][0] == 'm' && !IsPosReal(split[i+1], &val_d)) {
      strcpy(ERROR_MSG, "atom line: 'm[ass]' not followed by \
positive real number ");
      return false;
    } //}}}
  }
  // error - missing the mandatory n[ame] keyword or //{{{
  if (!name) {
    strcpy(ERROR_MSG, "atom line: missing 'n[ame]' keyword ");
    return false;
  } //}}}
  // error - if res[name] is present, there must be resid as well //{{{
  if ((!resid && resname) || (resid && !resname)) {
    strcpy(ERROR_MSG, "atom line: if 'res[name]' is present, \
'resid' must be too ");
//  ErrorPrintFull(file, file_line_count, split, words);
    return false;
  } //}}}
  // valid atom line
  return true;
} //}}}
bool VtfCheckBondLine(int words, char *split[]) { //{{{
  long val_i;
  // valid line b[ond] '<id>:[  ]<id> anything'
  // error - only one string or missing 'b[ond]' keyword
  if (words < 2 || split[0][0] != 'b') {
    return false;
  }
  // two strings - assume '<int>:<int>' and test it
  if (words == 2) {
    // ':'-split the second string
    char index[SPL_STR][SPL_LEN], string[SPL_LEN];
    strcpy(string, split[1]);
    int strings = SplitLine_old(index, string, ":");
    if (strings != 2 || !IsInteger(index[0], &val_i) || val_i < 0 ||
                        !IsInteger(index[1], &val_i) || val_i < 0) {
      strcpy(ERROR_MSG, "bond line: only 'b[ond] <int>:<int>' or \
'b[ond] <int>: <int>' is valid (for now)\n");
      return false;
    }
  }
  // more than two strings - assume '<int>: <int>' and test it
  if (words > 2) {
    split[1][strlen(split[1])-1] = '\0';
    if (!IsInteger(split[1], &val_i) || val_i < 0 ||
        !IsInteger(split[2], &val_i) || val_i < 0) {
      strcpy(ERROR_MSG, "bond line: only 'b[ond] <int>:<int>' or \
'b[ond] <int>: <int>' is valid (for now)\n");
      return false;
    }
  }
  return true;
} //}}}
// save data from a vsf line
// VtfAtomLineValues() //{{{
/*
 * Go through a confirmed atom line, picking out positions of values
 * corresponding to certain keywords. Returns 'array' where subscripts
 * correspond to n-th split[] for: 0..name, 1..mass, 2..charge, 3..radius,
 * 4..resame, 5..resid If not present, the corresponding element has - 1.
 */
int * VtfAtomLineValues(int words, char *split[]) {
  static int value[6];
  for (int i = 0; i < 6; i++) {
    value[i] = -1;
  }
  for (int i = 2; i < words; i += 2) {
    if (split[i][0] == 'n') {
      value[0] = i + 1;
    } else if (split[i][0] == 'm') {
      value[1] = i + 1;
    } else if (split[i][0] == 'q' || strcmp(split[i], "charge") == 0) {
      value[2] = i + 1;
    } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0) {
      value[3] = i + 1;
    } else if (strncmp(split[i], "resname", 3) == 0 &&
               strcmp(split[i], "resid") != 0) {
      value[4] = i + 1;
    } else if (strcmp(split[i], "resid") == 0) {
      value[5] = i + 1;
    }
  }
  return value;
} //}}}
 //}}}

SYSTEM FieldReadFull(char field_file[]) { //{{{
  SYSTEM System;
  InitSystem(&System);
  COUNT *Count = &System.Count;
  // read species
  FieldReadSpecies(field_file, &System);
  // fill System.Bead & System.Unbonded //{{{
  if (System.Unbonded > 0) {
    System.Bead = realloc(System.Bead, sizeof (BEAD) * Count->Bead);
    System.Unbonded = realloc(System.Unbonded,
                              sizeof *System.Unbonded * Count->Unbonded);
    int count = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      for (int j = 0; j < System.BeadType[i].Number; j++) {
        System.Bead[count].Type = i;
        System.Bead[count].Molecule = -1;
        System.Unbonded[count] = count;
        count++;
      }
    }
    if (count != Count->Unbonded) {
      strcpy(ERROR_MSG, "wrong number of unbonded beads; should never happen!");
      PrintError();
      exit(1);
    }
  } //}}}
  FieldReadMolecules(field_file, &System);
  RemoveExtraTypes(&System);
  MergeBeadTypes(&System, true);
  MergeMoleculeTypes(&System);
  FillSystemNonessentials(&System);
  CheckSystem(System, field_file);
  PrintCount(*Count);
  PrintBeadType(System);
  PrintMoleculeType(System);
  VtfWriteStruct("field.vsf", System);
  return System;
} //}}}
void FieldReadSpecies(char field_file[], SYSTEM *System) { //{{{
  int file_line_count = 0, words;
  char line[LINE], *split[SPL_STR];
  FILE *fr = OpenFile(field_file, "r");
  // skip till line starting with 'Species' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, LINE, line, &words,
                                  split, SPL_STR, " \t\n"))) {
    file_line_count++;
    if (words > 0 && strncasecmp(split[0], "species", 6) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    strcpy(ERROR_MSG, "premature end of file");
    PrintError();
    ErrorPrintFile(field_file, "\0");
    putc('\n', stderr);
    exit(1);
  } //}}}
  COUNT *Count = &System->Count;
  // read number of bead types //{{{
  long types;
  if (words < 2 || !IsNatural(split[1], &types)) {
    strcpy(ERROR_MSG, "incorrect 'Species' keyword line");
    PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
    exit(1);
  } //}}}
  // read bead types //{{{
  for (int i = 0; i < types; i++) {
    file_line_count++;
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (not enough 'species' lines)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    }
    double mass, charge;
    long number;
    // error - illegal 'species' line //{{{
    if (words < 4 || !IsPosReal(split[1], &mass) ||
        !IsReal(split[2], &charge) || !IsNatural(split[3], &number)) {
      strcpy(ERROR_MSG, "incorrect 'Species' line");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
    NewBeadType(&System->BeadType, &Count->BeadType, split[0], charge, mass,
                RADIUS);
    System->BeadType[Count->BeadType-1].Number = number;
    Count->Bead += number;
    Count->Unbonded += number;
  } //}}}
  fclose(fr);
} //}}}
void FieldReadMolecules(char field_file[], SYSTEM *System) { //{{{
  int file_line_count = 0, words;
  char line[LINE], *split[SPL_STR];
  FILE *fr = OpenFile(field_file, "r");
  // skip till line starting with 'Molecule' keyword //{{{
  bool test;
  while ((test = ReadAndSplitLine(fr, LINE, line, &words,
                                  split, SPL_STR, " \t\n"))) {
    file_line_count++;
    if (words > 0 && strncasecmp(split[0], "molecules", 7) == 0) {
      break;
    }
  }
  // error - missing 'Species' line
  if (!test) {
    strcpy(ERROR_MSG, "premature end of file");
    PrintError();
    ErrorPrintFile(field_file, "\0");
    putc('\n', stderr);
    exit(1);
  } //}}}
  // read number of types //{{{
  long val;
  if (words < 2 || !IsNatural(split[1], &val)) {
    strcpy(ERROR_MSG, "incorrect 'Molecules' keyword line");
    PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
    exit(1);
  } //}}}
  COUNT *Count = &System->Count;
  Count->MoleculeType = val;
  System->MoleculeType = realloc(System->MoleculeType,
                                 sizeof (MOLECULETYPE) * Count->MoleculeType);
  // read molecule types
  /*
   * Order of entries:
   *   1) <molecule name>
   *   2) nummol[s] <int>
   *   3) bead[s] <int>
   *      <int> lines: <bead name> <x> <y> <z>
   *   4) bond[s] <int> (optional)
   *      <int> lines: harm <id1> <id2> <k> <r_0>
   *   5) angle[s] <int> (optional)
   *      <int> lines: harm <id1> <id2> <id3> <k> <theta_0>
   *   6) dihed[rals] <int> (optional)
   *      <int> lines: harm <id1> <id2> <id3> <id4> <k> <theta_0>
   *   7) finish keyword
   */
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    InitMoleculeType(mt_i);
    // 1) name //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (missing a 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words == 0) {
      strcpy(ERROR_MSG, "missing molecule name");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      fprintf(stderr, "%s, line %d%s\n", ErrRed(),
              file_line_count, ErrColourReset());
      exit(1);
    }
    strncpy(mt_i->Name, split[0], MOL_NAME);
    mt_i->Name[MOL_NAME-1] = '\0'; // null-terminate! //}}}
    // 2) number of molecules //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words < 2 || strcasecmp(split[0], "nummols") != 0 ||
        !IsPosInteger(split[1], &val)) {
      strcpy(ERROR_MSG, "incorrect 'nummols' line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    }
    mt_i->Number = val; //}}}
    // 3) beads in the molecule //{{{
    // a) number of beads //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words < 2 || strncasecmp(split[0], "beads", 4) != 0 ||
        !IsPosInteger(split[1], &val)) {
      strcpy(ERROR_MSG, "incorrect 'beads' line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    }
    mt_i->nBeads = val; //}}}
    mt_i->Bead = malloc(sizeof *mt_i->Bead * mt_i->nBeads);
    VECTOR *coor = malloc(sizeof (VECTOR) * mt_i->nBeads);
    // b) beads themselves //{{{
    for (int j = 0; j < mt_i->nBeads; j++) {
      file_line_count++;
      // read a line //{{{
      if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
        strcpy(ERROR_MSG, "premature end of file \
(incomplete beads section for molecule entry)");
        PrintError();
        ErrorPrintFile(field_file, "\0");
        putc('\n', stderr);
        exit(1);
      } //}}}
      // error - incorrect line //{{{
      if (words < 4 || !IsReal(split[1], &coor[i].x) ||
                       !IsReal(split[2], &coor[i].y) ||
                       !IsReal(split[3], &coor[i].z)) {
        strcpy(ERROR_MSG, "incorrect bead coordinate line in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        exit(1);
      } //}}}
      int btype = FindBeadType(split[0], *System);
      // error - unknown type //{{{
      if (btype == -1) {
        strcpy(ERROR_MSG, "unknown bead in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        ErrorBeadType(*System);
        exit(1);
      } //}}}
      mt_i->Bead[j] = btype;
    } //}}}
    int count = Count->Bead,         // for filling Molecule[] & Bead[]
        mol_count = Count->Molecule; // for filling Molecule[] & Bead[]
    Count->Molecule += mt_i->Number;
    Count->Bead += mt_i->Number * mt_i->nBeads ;
    Count->Bonded += mt_i->Number * mt_i->nBeads ;
    System->Molecule = realloc(System->Molecule,
                               sizeof (MOLECULE) * Count->Molecule);
    System->Bead = realloc(System->Bead, sizeof (BEAD) * Count->Bead);
    System->Bonded = realloc(System->Bonded,
                             sizeof System->Bonded * Count->Bonded);
    // c) fill Molecule[] & Bead[] //{{{
    for (int j = 0; j < mt_i->Number; j++) {
      MOLECULE *mol = &System->Molecule[mol_count];
      mol->Type = i;
      mol->Index = mol_count;
      mol->Bead = malloc(sizeof *mol->Bead * mt_i->nBeads);
      for (int k = 0; k < mt_i->nBeads; k++) {
        BEAD *bead = &System->Bead[count];
        bead->Type = mt_i->Bead[k];
        bead->Molecule = mol_count;
        bead->Position.x = coor[k].x;
        bead->Position.y = coor[k].y;
        bead->Position.z = coor[k].z;
        System->BeadType[bead->Type].Number++;
        System->Bonded[count-Count->Unbonded] = count;
        mol->Bead[k] = count;
        count++;
      }
      mol_count++;
    }
    // pro forma checks //{{{
    if (count != Count->Bead) {
      strcpy(ERROR_MSG, "wrong number of beads; should never happen!");
      PrintError();
      exit(1);
    }
    if (mol_count != Count->Molecule) {
      strcpy(ERROR_MSG, "wrong number of molecules; should never happen!");
      PrintError();
      exit(1);
    } //}}}
    //}}}
    free(coor); //}}}
    // 4) bonds in the molecule (if present) //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words > 1 && strncasecmp(split[0], "bonds", 4) == 0) {
      // a) number of bonds //{{{
      if (!IsPosInteger(split[1], &val)) {
        strcpy(ERROR_MSG, "incorrect 'bonds' line in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        exit(1);
      }
      mt_i->nBonds = val; //}}}
      mt_i->Bond = malloc(sizeof *mt_i->Bond * mt_i->nBonds);
      // b) bonds themselves & bond types //{{{
      // TODO: for now, only harmonic bonds are considered
      for (int j = 0; j < mt_i->nBonds; j++) {
        file_line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, LINE, line, &words,
                              split, SPL_STR, " \t\n")) {
          strcpy(ERROR_MSG, "premature end of file \
(incomplete bonds section for molecule entry)");
          PrintError();
          ErrorPrintFile(field_file, "\0");
          putc('\n', stderr);
          exit(1);
        } //}}}
        long beads[2];
        PARAMS values;
        // error - incorrect line //{{{
        if (words < 5 || !IsPosInteger(split[1], &beads[0]) ||
                         !IsPosInteger(split[2], &beads[1]) ||
                         !IsReal(split[3], &values.a) || values.a < 0 ||
                         !IsReal(split[4], &values.b) || values.b < 0) {
          strcpy(ERROR_MSG, "incorrect bond line in a molecule entry");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        // error - bead index is too high //{{{
        if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads) {
          strcpy(ERROR_MSG, "bead index in a bond is too high");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        mt_i->Bond[j][0] = beads[0] - 1; // in FIELD, bead indices start from 1
        mt_i->Bond[j][1] = beads[1] - 1; //
        // find bond type //{{{
        int bond_type = -1;
        // find if this bond type already exists
        for (int k = 0; k < Count->BondType; k++) {
          if (System->BondType[k].a == values.a &&
              System->BondType[k].b == values.b) {
            bond_type = k;
            break;
          }
        }
        // create a new bond type if necessary
        if (bond_type == -1) {
          bond_type = Count->BondType;
          Count->BondType++;
          System->BondType = realloc(System->BondType,
                                     sizeof (PARAMS) * Count->BondType);
          System->BondType[bond_type].a = values.a;
          System->BondType[bond_type].b = values.b;
        } //}}}
        mt_i->Bond[j][2] = bond_type;
      } //}}}
    } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
      continue;
    } else {
      strcpy(ERROR_MSG, "unrecognised line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
    // 5) angles in the molecule (if present) //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words > 1 && strncasecmp(split[0], "angles", 4) == 0) {
      // a) number of angles //{{{
      if (!IsPosInteger(split[1], &val)) {
        strcpy(ERROR_MSG, "incorrect 'angles' line in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        exit(1);
      }
      mt_i->nAngles = val; //}}}
      mt_i->Angle = malloc(sizeof *mt_i->Angle * mt_i->nAngles);
      // b) angles themselves & angle types //{{{
      // TODO: for now, only harmonic angles are considered
      for (int j = 0; j < mt_i->nAngles; j++) {
        file_line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, LINE, line, &words,
                              split, SPL_STR, " \t\n")) {
          strcpy(ERROR_MSG, "premature end of file \
(incomplete angles section for molecule entry)");
          PrintError();
          ErrorPrintFile(field_file, "\0");
          putc('\n', stderr);
          exit(1);
        } //}}}
        long beads[3];
        PARAMS values;
        // error - incorrect line //{{{
        if (words < 5 || !IsPosInteger(split[1], &beads[0]) ||
                         !IsPosInteger(split[2], &beads[1]) ||
                         !IsPosInteger(split[3], &beads[2]) ||
                         !IsReal(split[4], &values.a) || values.a < 0 ||
                         !IsReal(split[5], &values.b) || values.b < 0) {
          strcpy(ERROR_MSG, "incorrect angle line in a molecule entry");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        // error - bead index is too high //{{{
        if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads) {
          strcpy(ERROR_MSG, "bead index in a angle is too high");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        mt_i->Angle[j][0] = beads[0] - 1; // in FIELD, bead indices start from 1
        mt_i->Angle[j][1] = beads[1] - 1; //
        mt_i->Angle[j][2] = beads[2] - 1; //
        // find angle type //{{{
        int angle_type = -1;
        // find if this angle type already exists
        for (int k = 0; k < Count->AngleType; k++) {
          if (System->AngleType[k].a == values.a &&
              System->AngleType[k].b == values.b) {
            angle_type = k;
            break;
          }
        }
        // create a new angle type if necessary
        if (angle_type == -1) {
          angle_type = Count->AngleType;
          Count->AngleType++;
          System->AngleType = realloc(System->AngleType,
                                      sizeof (PARAMS) * Count->AngleType);
          System->AngleType[angle_type].a = values.a;
          System->AngleType[angle_type].b = values.b;
        } //}}}
        mt_i->Angle[j][3] = angle_type;
      } //}}}
    } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
      continue;
    } else {
      strcpy(ERROR_MSG, "unrecognised line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
    // 6) dihedrals in the molecule (if present) //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words > 1 && strncasecmp(split[0], "dihedrals", 5) == 0) {
      // a) number of dihedrals //{{{
      if (!IsPosInteger(split[1], &val)) {
        strcpy(ERROR_MSG, "incorrect 'dihedrals' line in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        exit(1);
      }
      mt_i->nDihedrals = val; //}}}
      mt_i->Dihedral = malloc(sizeof *mt_i->Dihedral * mt_i->nDihedrals);
      // b) dihedrals themselves & dihedral types //{{{
      // TODO: for now, only harmonic dihedrals are considered
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        file_line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, LINE, line, &words,
                              split, SPL_STR, " \t\n")) {
          strcpy(ERROR_MSG, "premature end of file \
(incomplete dihedrals section for molecule entry)");
          PrintError();
          ErrorPrintFile(field_file, "\0");
          putc('\n', stderr);
          exit(1);
        } //}}}
        long beads[4];
        PARAMS values;
        // error - incorrect line //{{{
        if (words < 5 || !IsPosInteger(split[1], &beads[0]) ||
                         !IsPosInteger(split[2], &beads[1]) ||
                         !IsPosInteger(split[3], &beads[2]) ||
                         !IsPosInteger(split[4], &beads[3]) ||
                         !IsReal(split[5], &values.a) || values.a < 0 ||
                         !IsReal(split[6], &values.b) || values.b < 0) {
          strcpy(ERROR_MSG, "incorrect dihedral line in a molecule entry");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        // error - bead index is too high //{{{
        if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
            beads[2] > mt_i->nBeads || beads[3] > mt_i->nBeads) {
          strcpy(ERROR_MSG, "bead index in a dihedral is too high");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        mt_i->Dihedral[j][0] = beads[0] - 1; // in FIELD, indices start from 1
        mt_i->Dihedral[j][1] = beads[1] - 1; //
        mt_i->Dihedral[j][2] = beads[2] - 1; //
        mt_i->Dihedral[j][3] = beads[3] - 1; //
        // find dihedral type //{{{
        int dihedral_type = -1;
        // find if this dihedral type already exists
        for (int k = 0; k < Count->DihedralType; k++) {
          if (System->DihedralType[k].a == values.a &&
              System->DihedralType[k].b == values.b) {
            dihedral_type = k;
            break;
          }
        }
        // create a new dihedral type if necessary
        if (dihedral_type == -1) {
          dihedral_type = Count->DihedralType;
          Count->DihedralType++;
          System->DihedralType = realloc(System->DihedralType,
                                         sizeof (PARAMS) * Count->DihedralType);
          System->DihedralType[dihedral_type].a = values.a;
          System->DihedralType[dihedral_type].b = values.b;
        } //}}}
        mt_i->Dihedral[j][4] = dihedral_type;
      } //}}}
    } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
      continue;
    } else {
      strcpy(ERROR_MSG, "unrecognised line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
    // 6) impropers in the molecule (if present) //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words > 1 && strncasecmp(split[0], "impropers", 6) == 0) {
      // a) number of impropers //{{{
      if (!IsPosInteger(split[1], &val)) {
        strcpy(ERROR_MSG, "incorrect 'impropers' line in a molecule entry");
        PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
        exit(1);
      }
      mt_i->nImpropers = val; //}}}
      mt_i->Improper = malloc(sizeof *mt_i->Improper * mt_i->nImpropers);
      // b) impropers themselves & improper types //{{{
      // TODO: for now, only harmonic impropers are considered
      for (int j = 0; j < mt_i->nImpropers; j++) {
        file_line_count++;
        // read a line //{{{
        if (!ReadAndSplitLine(fr, LINE, line, &words,
                              split, SPL_STR, " \t\n")) {
          strcpy(ERROR_MSG, "premature end of file \
(incomplete impropers section for molecule entry)");
          PrintError();
          ErrorPrintFile(field_file, "\0");
          putc('\n', stderr);
          exit(1);
        } //}}}
        long beads[4];
        PARAMS values;
        // error - incorrect line //{{{
        if (words < 5 || !IsPosInteger(split[1], &beads[0]) ||
                         !IsPosInteger(split[2], &beads[1]) ||
                         !IsPosInteger(split[3], &beads[2]) ||
                         !IsPosInteger(split[4], &beads[3]) ||
                         !IsReal(split[5], &values.a) || values.a < 0 ||
                         !IsReal(split[6], &values.b) || values.b < 0) {
          strcpy(ERROR_MSG, "incorrect improper line in a molecule entry");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        // error - bead index is too high //{{{
        if (beads[0] > mt_i->nBeads || beads[1] > mt_i->nBeads ||
            beads[2] > mt_i->nBeads || beads[3] > mt_i->nBeads) {
          strcpy(ERROR_MSG, "bead index in a improper is too high");
          PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
          exit(1);
        } //}}}
        mt_i->Improper[j][0] = beads[0] - 1; // in FIELD, indices start from 1
        mt_i->Improper[j][1] = beads[1] - 1; //
        mt_i->Improper[j][2] = beads[2] - 1; //
        mt_i->Improper[j][3] = beads[3] - 1; //
        // find improper type //{{{
        int improper_type = -1;
        // find if this improper type already exists
        for (int k = 0; k < Count->ImproperType; k++) {
          if (System->ImproperType[k].a == values.a &&
              System->ImproperType[k].b == values.b) {
            improper_type = k;
            break;
          }
        }
        // create a new improper type if necessary
        if (improper_type == -1) {
          improper_type = Count->ImproperType;
          Count->ImproperType++;
          System->ImproperType = realloc(System->ImproperType,
                                         sizeof (PARAMS) * Count->ImproperType);
          System->ImproperType[improper_type].a = values.a;
          System->ImproperType[improper_type].b = values.b;
        } //}}}
        mt_i->Improper[j][4] = improper_type;
      } //}}}
    } else if (words > 0 && strcasecmp(split[0], "finish") == 0) {
      continue;
    } else {
      strcpy(ERROR_MSG, "unrecognised line in a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
    // finish keyword //{{{
    file_line_count++;
    // read a line //{{{
    if (!ReadAndSplitLine(fr, LINE, line, &words, split, SPL_STR, " \t\n")) {
      strcpy(ERROR_MSG, "premature end of file (incomplete 'molecule' entry)");
      PrintError();
      ErrorPrintFile(field_file, "\0");
      putc('\n', stderr);
      exit(1);
    } //}}}
    if (words == 0 || strcasecmp(split[0], "finish") != 0) {
      strcpy(ERROR_MSG, "missing 'finish' at the end of a molecule entry");
      PrintErrorFileLine(field_file, "\0", file_line_count, split, words);
      exit(1);
    } //}}}
  }
  Count->HighestResid = Count->Molecule;
  fclose(fr);
} //}}}

#if 0
// TODO will be changed - FIELD file
// ReadFieldPbc() //{{{
/*
 * Function reading box size from the first line of the FIELD-like file.
 */
bool ReadFieldPbc(char *field, VECTOR *BoxLength) {
  bool pbc = false; // assume the first line doesn't contain box size
  // open FIELD-like file
  FILE *fr = OpenFile(field, "r");
  // read first line
  char line[LINE], split[SPL_STR][SPL_LEN];
  fgets(line, sizeof line, fr);
  int words = SplitLine_old(split, line, " \t");
  /*
   * box size must be 'pbc <double> <double> <double>'; test
   * 1) number of strings
   * 2) 'pbc' string
   * 3) positive doubles
   */
  if (words >= 4 && // 1)
      strcasecmp(split[0], "pbc") == 0 && // 2)
      IsPosReal_old(split[1]) && //
      IsPosReal_old(split[2]) && // 3)
      IsPosReal_old(split[3])) { //
    (*BoxLength).x = atof(split[1]);
    (*BoxLength).y = atof(split[2]);
    (*BoxLength).z = atof(split[3]);
    pbc = true;
  }
  fclose(fr);
  return pbc;
} //}}}
// ReadFieldBeadType() //{{{
/*
 * Function reading the species section of the FIELD-like file.
 */
void ReadFieldBeadType(char *field, COUNTS *Counts,
                       BEADTYPE *BeadType[], BEAD *Bead[]) {
  // open FIELD-like file
  FILE *fr = OpenFile(field, "r");
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read number of bead types //{{{
  bool missing = true; // assume species keyword is missing
  while(fgets(line, sizeof line, fr)) {
    int words = SplitLine_old(split, line, " \t");
    if (strcasecmp(split[0], "species") == 0) {
      missing = false; // species keyword is present
      // check if the next string is a number
      if (words < 2 || !IsInteger_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing number of species\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", field);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'species' keyword\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof (BEADTYPE));
  // read info about bead types //{{{
  (*Counts).Unbonded = 0;
  (*Counts).BeadsCoor = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    fgets(line, sizeof line, fr);
    int words = SplitLine_old(split, line, " \t");
    // Error on the line //{{{
    /*
     * species lines must be <name> <mass> <charge> <number> and cannot contain
     * blamk lines, so error when:
     *   1) empty line
     *   2) fewer than four strings
     *   3) second string isn't a positive double (mass)
     *   4) third string isn't a double (charge)
     *   5) fifth string isn't an integer (unbonded beads)
     */
    if (words == 1 && split[0][0] == '\0') { // 1)
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - no blank lines permitted in the species section\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    } else if (words < 4 ||                  // 2)
               !IsPosReal_old(split[1]) ||     // 3)
               !IsReal_old(split[2]) ||        // 4)
               !IsInteger_old(split[3])) {       // 5)
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - wrong species line\n");
      ColourReset(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
//  P_IGNORE(-Wformat-truncation);
    snprintf((*BeadType)[i].Name,BEAD_NAME, "%s", split[0]);
//  P_POP;
    (*BeadType)[i].Mass = atof(split[1]);
    (*BeadType)[i].Charge = atof(split[2]);
    (*BeadType)[i].Number = atoi(split[3]);
    (*Counts).Unbonded += (*BeadType)[i].Number;
  } //}}}
  (*Counts).BeadsCoor = (*Counts).Unbonded;
  // allocate & fill Bead array //{{{
  *Bead = calloc((*Counts).Unbonded, sizeof (BEAD));
  int count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    for (int j = 0; j < (*BeadType)[i].Number; j++) {
      (*Bead)[count].Type = i;
      (*Bead)[count].Molecule = -1;
      count++;
    }
  } //}}}
  fclose(fr);
} //}}}
// ReadFieldMolecules() //{{{
/*
 * Function reading the molecules section of the FIELD-like file.
 */
void ReadFieldMolecules(char *field, COUNTS *Counts,
                        BEADTYPE *BeadType[], BEAD *Bead[],
                        MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[],
                        PARAMS *bond_type[], PARAMS *angle_type[],
                        PARAMS *dihedral_type[]) {
  FILE *fr = OpenFile(field, "r");
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read number of molecule types //{{{
  bool missing = true; // assume molecule keyword is missing
  while(fgets(line, sizeof line, fr)) {
    int words = SplitLine_old(split, line, " \t");
    if (strncasecmp(split[0], "molecule", 8) == 0) {
      missing = false; // molecule keyword is present
      // error - 'molecule' not followed by an integer
      if (!IsInteger_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - missing number of molecule types\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  // error - no molecule keyword
  if (missing) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", field);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'molecule' line\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}
  // allocate molecule type struct
  *MoleculeType = calloc((*Counts).TypesOfMolecules, sizeof (MOLECULETYPE));
  // test proper number of 'finish' keywords //{{{
  fpos_t position;
  fgetpos(fr, &position); // save pointer position in field
  int count = 0;
  // count 'finish' keywords
  while(fgets(line, sizeof line, fr)) {
    SplitLine_old(split, line, " \t");
    if (strcmp(split[0], "finish") == 0) {
      count++;
    }
  }
  // error - fewer 'finish'es than molecule types
  if (count < (*Counts).TypesOfMolecules) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", field);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'finish' keyword\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  // restore the file pointer position
  fsetpos(fr, &position); //}}}

  // read info about molecule types //{{{
  /*
   * Order of entries:
   *   1) <molecule name>
   *   2) nummol[s] <int>
   *   3) bead[s] <int>
   *      <int> lines: <bead name> <x> <y> <z>
   *   4) bond[s] <int> (optional)
   *      <int> lines: harm <id1> <id2> <k> <r_0>
   *   5) angle[s] <int> (optional)
   *      <int> lines: harm <id1> <id2> <id3> <k> <theta_0>
   *   6) dihed[rals] <int> (optional)
   *      <int> lines: harm <id1> <id2> <id3> <id4> <k> <theta_0>
   *   7) finish keyword
   */
  fgetpos(fr, &position); // save pointer position
  (*Counts).TypesOfBonds = 0;
  (*Counts).TypesOfAngles = 0;
  (*Counts).TypesOfDihedrals = 0;
  *bond_type = malloc(sizeof (PARAMS) * 1);
  *angle_type = malloc(sizeof (PARAMS) * 1);
  *dihedral_type = malloc(sizeof (PARAMS) * 1);
  // stored bond & angle types - temporary //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // 1) //{{{
    fgets(line, sizeof line, fr);
    SplitLine_old(split, line, " \t");
    if (split[0][0] == '\0') {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - expecting molecule name\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    }
    // copy name to MOLECULETYPE
//  P_IGNORE(-Wformat-truncation);
    snprintf((*MoleculeType)[i].Name, MOL_NAME, "%s", split[0]);
//  P_POP; //}}}
    // 2) //{{{
    fgets(line, sizeof line, fr);
    int words = SplitLine_old(split, line, " \t");
    if (strncasecmp(split[0], "nummols", 6) != 0 ||
        words < 2 || !IsInteger_old(split[1])) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - expecting 'nummol[s] <number of molecules>'\n");
      ColourReset(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].Number = atoi(split[1]); //}}}
    // 3) //{{{
    // read number of beads
    fgets(line, sizeof line, fr);
    words = SplitLine_old(split, line, " \t");
    // error - wrong keyword line //{{{
    if (strncasecmp(split[0], "beads", 4) != 0 ||
        words < 2 || !IsInteger_old(split[1])) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - expecting 'bead[s] <number of beads>'\n");
      ColourReset(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
    (*MoleculeType)[i].nBeads = atoi(split[1]);
    (*MoleculeType)[i].Bead = malloc(sizeof *(*MoleculeType)[i].Bead *
                                     (*MoleculeType)[i].nBeads);
    memset((*MoleculeType)[i].Bead, 0,
           sizeof *(*MoleculeType)[i].Bead *(*MoleculeType)[i].nBeads);
    // read bead info
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof line, fr);
      words = SplitLine_old(split, line, " \t");
      // error - bead line must be '<name> <double> <double> <double>' //{{{
      if (words < 4 ||
          !IsReal_old(split[1]) || !IsReal_old(split[2]) || !IsReal_old(split[3])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - wrong bead line in molecule ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      int type = FindBeadType_old(split[0], *Counts, *BeadType);
      // error - unknown bead type //{{{
      if (type == -1) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - non-existent bead name ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", split[0]);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " in molecule ");
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ColourReset(STDERR_FILENO);
        ErrorBeadType_old(*Counts, *BeadType);
        exit(1);
      } //}}}
      (*MoleculeType)[i].Bead[j] = type;
      // TODO: Read coordinates?
    } //}}}
    // 4) //{{{
    // read number of bonds in the molecule
    fgets(line, sizeof line, fr);
    words = SplitLine_old(split, line, " \t");
    // are bonds present?
    if (strncasecmp(split[0], "bonds", 4) == 0) {
      // error - wrong keyword line //{{{
      if (strncasecmp(split[0], "bonds", 4) != 0 ||
          words < 2 || !IsInteger_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - expecting 'bond[s] <number of bonds>'\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      (*MoleculeType)[i].nBonds = atoi(split[1]);
      (*MoleculeType)[i].Bond = malloc(sizeof *(*MoleculeType)[i].Bond *
                                       (*MoleculeType)[i].nBonds);
      memset((*MoleculeType)[i].Bond, 0,
             sizeof *(*MoleculeType)[i].Bond * (*MoleculeType)[i].nBonds);
      // allocate memory for bonds
      for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
        (*MoleculeType)[i].Bond[j][2] = -1; // no bond type assigned
      }
      // read bond info
      for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error - bead line must be '<name> <id1> <id2> <double> <double>' //{{{
        if (words < 5 || !IsInteger_old(split[1]) || !IsInteger_old(split[2]) ||
            !IsPosReal_old(split[3]) || !IsPosReal_old(split[4])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", field);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - wrong bond line in molecule ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        // error - wrong bead index //{{{
        if (atoi(split[1]) > (*MoleculeType)[i].nBeads || atoi(split[1]) < 1 ||
            atoi(split[2]) > (*MoleculeType)[i].nBeads || atoi(split[2]) < 1) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", field);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - wrong bead index in a bond in molecule ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*MoleculeType)[i].Bond[j][0] = atoi(split[1]) - 1;
        (*MoleculeType)[i].Bond[j][1] = atoi(split[2]) - 1;
        // find if the bond type is new
        bool known = false; // assume it's new
        for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
          if ((*bond_type)[k].a == atof(split[3]) &&
              (*bond_type)[k].b == atof(split[4])) {
            known = true; // it's not new
            break;
          }
        }
        if (!known) { // add new bond type if necessary
          (*Counts).TypesOfBonds++;
          *bond_type = realloc(*bond_type, sizeof (PARAMS) *
                               (*Counts).TypesOfBonds);
          (*bond_type)[(*Counts).TypesOfBonds-1].a = atof(split[3]);
          (*bond_type)[(*Counts).TypesOfBonds-1].b = atof(split[4]);
        }
      }
    } //}}}
    // 5) //{{{
    fpos_t position2;
    fgetpos(fr, &position2); // save file pointer
    fgets(line, sizeof line, fr);
    words = SplitLine_old(split, line, " \t");
    // are angles present?
    if (strncasecmp(split[0], "angles", 5) == 0) {
      // error - missing number of angles //{{{
      if (words < 2 || !IsInteger_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - expecting 'angle[s] <number of angles>'\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      (*MoleculeType)[i].nAngles = atoi(split[1]);
      (*MoleculeType)[i].Angle = malloc(sizeof *(*MoleculeType)[i].Angle *
                                        (*MoleculeType)[i].nAngles);
      memset((*MoleculeType)[i].Angle, 0,
             sizeof *(*MoleculeType)[i].Angle * (*MoleculeType)[i].nAngles);
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        (*MoleculeType)[i].Angle[j][3] = -1; // no angle type assigned
      }
      // read info about angles
      for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error - wrong angle line //{{{
        // '<type> 3x<id> <double> <double>'
        if (words < 6 || !IsInteger_old(split[1]) ||
            !IsInteger_old(split[2]) || !IsInteger_old(split[3]) ||
            !IsReal_old(split[4]) || !IsReal_old(split[5])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", field);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - wrong angle line in molecule ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*MoleculeType)[i].Angle[j][0] = atoi(split[1]) - 1;
        (*MoleculeType)[i].Angle[j][1] = atoi(split[2]) - 1;
        (*MoleculeType)[i].Angle[j][2] = atoi(split[3]) - 1;
        // find whether the angle type is known
        bool known = false; // assume it's a new angle type
        for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
          if ((*angle_type)[k].a == atof(split[4]) &&
              (*angle_type)[k].b == atof(split[5])) {
            (*MoleculeType)[i].Angle[j][3] = k;
            known = true; // it's not a new angle type
            break;
          }
        }
        if (!known) { // add a new angle type if necessary
          (*Counts).TypesOfAngles++;
          *angle_type = realloc(*angle_type, sizeof (PARAMS) *
                                (*Counts).TypesOfAngles);
          (*angle_type)[(*Counts).TypesOfAngles-1].a = atof(split[4]);
          (*angle_type)[(*Counts).TypesOfAngles-1].b = atof(split[5]);
          (*MoleculeType)[i].Angle[j][3] = (*Counts).TypesOfDihedrals - 1;
        }
      }
    // error - extra bonds (from previous section)
    } else if (strncasecmp(split[0], "harm", 4) == 0) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - extra bond line in molecule ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
      ErrorPrintLine(split, words);
      exit(1);
    } else {
      // reset file pointer as the angle section isn't present
      fsetpos(fr, &position2);
    } //}}}
    // 6) //{{{
    fgetpos(fr, &position2);
    fgets(line, sizeof line, fr); // save file pointer
    words = SplitLine_old(split, line, " \t");
    // are dihedrals present?
    if (strncasecmp(split[0], "dihedrals", 5) == 0) {
      // get number of dihedrals
      // error - wrong number of dihedrals //{{{
      if (words < 2 || !IsInteger_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", field);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - expecting 'dihed[rals] <number of angles>'\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // allocate memory for dihedrals
      (*MoleculeType)[i].nDihedrals = atoi(split[1]);
      (*MoleculeType)[i].Dihedral = malloc(sizeof *(*MoleculeType)[i].Dihedral *
                                           (*MoleculeType)[i].nDihedrals);
      memset((*MoleculeType)[i].Dihedral, 0,
             sizeof *(*MoleculeType)[i].Dihedral *
             (*MoleculeType)[i].nDihedrals);
      for (int j = 0; j < (*MoleculeType)[i].nDihedrals; j++) {
        (*MoleculeType)[i].Dihedral[j][4] = -1; // no dihedral type assigned
      }
      // get info about dihedrals
      for (int j = 0; j < (*MoleculeType)[i].nDihedrals; j++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error - wrong dihedral line //{{{
        // '<type> 4x<id> <double> <double>'
        // TODO: do something about having 4x<id> 0 0 (IsPosReal() checks > 0)
        if (words < 7 || !IsInteger_old(split[1]) || !IsInteger_old(split[2]) ||
            !IsInteger_old(split[3]) || !IsInteger_old(split[4]) ||
            !IsReal_old(split[5]) || !IsReal_old(split[6])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", field);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - wrong dihedral line in molecule ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*MoleculeType)[i].Dihedral[j][0] = atoi(split[1]) - 1;
        (*MoleculeType)[i].Dihedral[j][1] = atoi(split[2]) - 1;
        (*MoleculeType)[i].Dihedral[j][2] = atoi(split[3]) - 1;
        (*MoleculeType)[i].Dihedral[j][3] = atoi(split[4]) - 1;
        // find if this dihedral type is known
        bool known = false; // assume it's a new type
        for (int k = 0; k < (*Counts).TypesOfDihedrals; k++) {
          if ((*dihedral_type)[k].a == atof(split[5]) &&
              (*dihedral_type)[k].b == atof(split[6])) {
            (*MoleculeType)[i].Dihedral[j][4] = k;
            known = true; // it's not a new type
            break;
          }
        }
        if (!known) { // create a new dihedral type if necessary
          (*Counts).TypesOfDihedrals++;
          *dihedral_type = realloc(*dihedral_type, sizeof (PARAMS) *
                                   (*Counts).TypesOfDihedrals);
          (*dihedral_type)[(*Counts).TypesOfDihedrals-1].a = atof(split[5]);
          (*dihedral_type)[(*Counts).TypesOfDihedrals-1].b = atof(split[6]);
          (*MoleculeType)[i].Dihedral[j][4] = (*Counts).TypesOfDihedrals - 1;
        }
      }
    // error - extra bonds or angles (from previous section)
    } else if (strncasecmp(split[0], "harm", 4) == 0) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - extra bond or angle line in molecule ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
      ErrorPrintLine(split, words);
      exit(1);
    } else {
      // reset file pointer as the dihedral section isn't present
      fsetpos(fr, &position2);
    } //}}}
    (*Counts).Bonded += (*MoleculeType)[i].Number * (*MoleculeType)[i].nBeads;
    (*Counts).Molecules += (*MoleculeType)[i].Number;
    // skip till 'finish' //{{{
    fsetpos(fr, &position2);
    while(fgets(line, sizeof line, fr)) {
      SplitLine_old(split, line, " \t");
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}
  //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  // calculate molecule masses //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*MoleculeType)[i].Mass += (*BeadType)[btype].Mass;
    }
  } //}}}

  // update the number beads in the system //{{{
  (*Counts).BeadsCoor += (*Counts).Bonded;
  (*Counts).BeadsTotal = (*Counts).BeadsCoor;
  *Bead = realloc(*Bead, sizeof (BEAD) * (*Counts).BeadsCoor);
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      (*BeadType)[btype].Number += (*MoleculeType)[i].Number;
    }
  } //}}}

  // read coordinates of bonded beads & assign bond and angle types //{{{
  // get the coordinates from FIELD to each molecule
  count = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    // skip name, nummols, and beads lines //{{{
    fgets(line, sizeof line, fr);
    fgets(line, sizeof line, fr);
    fgets(line, sizeof line, fr); //}}}
    // read bead lines //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      fgets(line, sizeof line, fr);
      SplitLine_old(split, line, " \t");
      for (int k = 0; k < (*MoleculeType)[i].Number; k++) {
        int id = count + k * (*MoleculeType)[i].nBeads;
        (*Bead)[id].Position.x = atof(split[1]);
        (*Bead)[id].Position.y = atof(split[2]);
        (*Bead)[id].Position.z = atof(split[3]);
      }
      count++;
    } //}}}
    // set count to the first bead of the next molecule type
    count += ((*MoleculeType)[i].Number - 1) * (*MoleculeType)[i].nBeads;
    // skip bonds line
    fgets(line, sizeof line, fr);
    // read bond types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nBonds; j++) {
      fgets(line, sizeof line, fr);
      SplitLine_old(split, line, " \t");
      for (int k = 0; k < (*Counts).TypesOfBonds; k++) {
        if ((*bond_type)[k].a == atof(split[3]) && (*bond_type)[k].b == atof(split[4])) {
          (*MoleculeType)[i].Bond[j][2] = k;
          break;
        }
      }
    } //}}}
    // read angle types //{{{
    for (int j = 0; j < (*MoleculeType)[i].nAngles; j++) {
      // skip 'angles' line at the beginning of the angle section //{{{
      if (j == 0) {
        fgets(line, sizeof line, fr);
      } //}}}
      fgets(line, sizeof line, fr);
      SplitLine_old(split, line, " \t");
      for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
        if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
          (*MoleculeType)[i].Angle[j][3] = k;
          break;
        }
      }
    } //}}}
    // skip till 'finish' //{{{
    while(fgets(line, sizeof line, fr)) {
      SplitLine_old(split, line, " \t");
      if (strcasecmp(split[0], "finish") == 0) {
        break;
      }
    } //}}}
  } //}}}

  // return the file pointer to the beginning of the molecules section
  fsetpos(fr, &position);

  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    SortBonds((*MoleculeType)[i].Bond, (*MoleculeType)[i].nBonds);
    SortAngles((*MoleculeType)[i].Angle, (*MoleculeType)[i].nAngles);
  }

  fclose(fr);
} //}}}
// ReadField() //{{{
/*
 * Function reading the FIELD-like file; it completely fills provided structs,
 * overriding any possible data in there. If the FIELD-like file is a source of
 * some additional structure information, new structs must be used, and the
 * data copied from there afterwards.
 */
void ReadField(char *field, VECTOR *BoxLength, COUNTS *Counts,
               BEADTYPE *BeadType[], BEAD *Bead[], int *Index[],
               MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[],
               PARAMS *bond_type[], PARAMS *angle_type[],
               PARAMS *dihedral_type[]) {

  // read pbc if required //{{{
  if (BoxLength != NULL && !ReadFieldPbc(field, BoxLength)) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", field);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - first line must start with box size, ");
    fprintf(stderr, "i.e., 'pbc <double> <double> <double>'\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}
  ReadFieldBeadType(field, Counts, BeadType, Bead);
  // FIELD doesn't contain bead radius, so fill it with impossible values
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Radius = RADIUS;
  }
  ReadFieldMolecules(field, Counts, BeadType, Bead, MoleculeType, Molecule,
                     bond_type, angle_type, dihedral_type);
//// allocate Bead[].Aggregate array - needed only to free() //{{{
//for (int i = 0; i < (*Counts).BeadsCoor; i++) {
//  (*Bead)[i].Aggregatexxx = calloc(10, sizeof *(*Bead)[i].Aggregatexxx);
//} //}}}
  FillMolBTypes_old((*Counts).TypesOfMolecules, MoleculeType);
  // fill Molecule & Bead structs //{{{
  *Molecule = calloc((*Counts).Molecules, sizeof (MOLECULE));
  int count_mol = 0, count_bead = (*Counts).Unbonded;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].Number; j++) {
      (*Molecule)[count_mol].Type = i;
      (*Molecule)[count_mol].Bead = malloc(sizeof *(*Molecule)[count_mol].Bead *
                                           (*MoleculeType)[i].nBeads);
      for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
        int btype = (*MoleculeType)[i].Bead[k];
        (*Molecule)[count_mol].Bead[k] = count_bead;
        (*Bead)[count_bead].Molecule = count_mol;
        (*Bead)[count_bead].Type = btype;
        count_bead++;
      }
      count_mol++;
    }
  } //}}}
  // fill Index - I don't thinks it's used right now //{{{
  *Index = malloc(sizeof **Index * (*Counts).BeadsCoor);
  for (int i = 0; i < (*Counts).BeadsCoor; i++) {
    (*Index)[i] = i;
  } //}}}
  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, field);
} //}}}
// TODO will be changed - lammps data file
// ReadLmpData() //{{{
void ReadLmpData(char *data_file, int *bonds, PARAMS *bond_type[],
                 int *angles, PARAMS *angle_type[],
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE *BeadType[], BEAD *Bead[], int *Index[],
                 MOLECULETYPE *MoleculeType[], MOLECULE *Molecule[]) {
  FILE *fr = OpenFile(data_file, "r");

  // ignore first line (comment) //{{{
  char line[LINE];
  fgets(line, sizeof line, fr);
  // if the line is too long, skip the rest of it
  if (strcspn(line, "\n") == (LINE-1)) {
    while (getc(fr) != '\n')
      ;
  } //}}}

  // read data file header //{{{
  // data file header lines must start with a number (or '#' for comment),
  // therefore read until something else is encountered
  char split[SPL_STR][SPL_LEN];
  int words;
  do {
    fgets(line, sizeof line, fr);
    words = SplitLine_old(split, line, " \t");
    // number of atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'atoms' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).BeadsTotal = atoi(split[0]);
      (*Counts).BeadsCoor = atoi(split[0]);
    } //}}}
    // number of bonds //{{{
    if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'bonds' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *bonds = atoi(split[0]);
    } //}}}
    // number of angles //{{{
    if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'angles' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *angles = atoi(split[0]);
    } //}}}
    // number of bead types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'atom types' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[0]);
    } //}}}
    // number of bond types //{{{
    if (words > 2 && strcmp(split[1], "bond") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'bond types' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBonds = atoi(split[0]);
    } //}}}
    // number of angle types //{{{
    if (words > 2 && strcmp(split[1], "angle") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger_old(split[0])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'angle types' keyword must be preceded by integer\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfAngles = atoi(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "xlo") == 0 && strcmp(split[3], "xhi") == 0) {
      if (!IsReal_old(split[0]) || !IsReal_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'xlo xhi' keyword must be preceded by two floats\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).x = atof(split[1]) - atof(split[0]);
      (*box_lo).x = atof(split[0]);
    } //}}}
    // box length in y //{{{
    if (words > 3 && strcmp(split[2], "ylo") == 0 && strcmp(split[3], "yhi") == 0) {
      if (!IsReal_old(split[0]) || !IsReal_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'ylo yhi' keyword must be preceded by two floats\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).y = atof(split[1]) - atof(split[0]);
      (*box_lo).y = atof(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "zlo") == 0 && strcmp(split[3], "zhi") == 0) {
      if (!IsReal_old(split[0]) || !IsReal_old(split[1])) {
        ErrorPrintError_old();
        ColourChange(STDERR_FILENO, YELLOW);
        fprintf(stderr, "%s", data_file);
        ColourChange(STDERR_FILENO, RED);
        fprintf(stderr, " - 'zlo zhi' keyword must be preceded by two floats\n");
        ColourReset(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).z = atof(split[1]) - atof(split[0]);
      (*box_lo).z = atof(split[0]);
    } //}}}
  } while (words == 0 ||
           split[0][0] == '#' ||
           IsReal_old(split[0]) ||
           IsInteger_old(split[0])); //}}}

  // some error checking //{{{
  if ((*Counts).TypesOfBeads == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", data_file);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'atom types' line (or is 0)\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  if ((*Counts).BeadsTotal == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", data_file);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'atoms' line (or is 0)\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).x == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", data_file);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'xlo xhi' line (or is 0)\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).y == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", data_file);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'ylo yhi' line (or is 0)\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).z == 0) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "%s", data_file);
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - missing 'zlo zhi' line (or is 0)\n\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof (BEADTYPE));
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    sprintf((*BeadType)[i].Name, "bead%d", i+1);
//  (*BeadType)[i].Use = true;
//  (*BeadType)[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  *Bead = calloc((*Counts).BeadsCoor, sizeof (BEAD));
//for (int i = 0; i < (*Counts).BeadsCoor; i++) {
//  (*Bead)[i].Aggregatexxx = malloc(sizeof *(*Bead)[i].Aggregatexxx * 1);
//} //}}}

  *Index = calloc((*Counts).BeadsCoor, sizeof **Index);
  *bond_type = calloc((*Counts).TypesOfBonds, sizeof (PARAMS));
  *angle_type = calloc((*Counts).TypesOfAngles, sizeof (PARAMS));

  // read body of data file //{{{
  int test,
      *mols = NULL, // number of beads in each molecule
      monomer = 0; // monomer beads designated by mol_ID = 0 in data file
  while ((test = getc(fr)) != EOF) {
    ungetc(test, fr);
    // atom masses //{{{
    if (words > 0 && strcmp(split[0], "Masses") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 2 || !IsInteger_old(split[0]) || !IsPosReal_old(split[1])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each line in 'Masses' section must start with '<int> <float>'\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*BeadType)[i].Mass = atof(split[1]);
        // if there's a comment at the end of the line, consider it bead name
        if (words > 2 && split[2][0] == '#') {
          if (strlen(split[2]) == 1 && words > 3) { // comment form '# name'
//          P_IGNORE(-Wformat-truncation);
            // BEAD_NAME is max string length, i.e., array is longer
            snprintf((*BeadType)[i].Name, BEAD_NAME+1, "%s", split[3]);
//          P_POP;
          } else if (strlen(split[2]) > 1) { // comment form '#name'
            for (int j = 0; j < strlen(split[2]); j++) {
              split[2][j] = split[2][j+1];
            }
            strncpy((*BeadType)[i].Name, split[2], BEAD_NAME);
          }
        }
      }
    } //}}}
    // bond coefficients //{{{
    if (words > 1 && strcmp(split[0], "Bond") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfBonds; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger_old(split[0]) || !IsPosReal_old(split[1]) || !IsPosReal_old(split[2])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each line in 'Bond Coeffs' section must start with '<int> <float> <float>'\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*bond_type)[i].a = atof(split[1]);
        (*bond_type)[i].b = atof(split[2]);
      }
    } //}}}
    // angle coefficients //{{{
    if (words > 1 && strcmp(split[0], "Angle") == 0 && strcmp(split[1], "Coeffs") == 0) {
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // get mass of every bead
      for (int i = 0; i < (*Counts).TypesOfAngles; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger_old(split[0]) || !IsPosReal_old(split[1]) || !IsPosReal_old(split[2])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each line in 'Angle Coeffs' section must start with '<int> <float> <float>'\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        (*angle_type)[i].a = atof(split[1]);
        (*angle_type)[i].b = atof(split[2]);
      }
    } //}}}
    // atoms section //{{{
    if (words > 0 && strcmp(split[0], "Atoms") == 0) {
      // array for number of beads in each molecule
      mols = calloc((*Counts).BeadsCoor, sizeof *mols);
      // array for list of beads in each molecule
      // bead_mols[i] ... molecule's id; bead_mols[][i] ... molecule's beads
      int **bead_mols = calloc((*Counts).BeadsCoor, sizeof **bead_mols);
      for (int i = 0; i < (*Counts).BeadsCoor; i++) {
        bead_mols[i] = malloc(sizeof *bead_mols[i] * 1);
      }
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < (*Counts).BeadsCoor; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // format of each line: <id> <mol_id> <btype> <charge> <x> <y> <z>
        // Error - incorrect format //{{{
        if (words < 7 ||
            !IsInteger_old(split[0]) || !IsInteger_old(split[1]) || !IsInteger_old(split[2]) ||
            !IsReal_old(split[3]) || !IsReal_old(split[4]) || !IsReal_old(split[5]) || !IsReal_old(split[6])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each 'Atoms' line must be <id> <mol_id> <bead type> <charge> <x> <y> <z>\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int id = atoi(split[0]) - 1, // in lammps, these start at 1
            mol_id = atoi(split[1]) - 1, // in lammps, molecules start with 1; unbonded atoms can be 0
            btype = atoi(split[2]) - 1; // in lammps, these start at 1

        if (mol_id == -1) { // corresponds to 0 in the data file
          monomer++;
          (*Bead)[id].Molecule = -1;
        } else { // possibly in a molecule (if more beads share its mol_id)
          mols[mol_id]++;
          bead_mols[mol_id] = realloc(bead_mols[mol_id],
                                      sizeof *bead_mols[mol_id] * mols[mol_id]);
          bead_mols[mol_id][mols[mol_id]-1] = id;
          (*Bead)[id].Molecule = mol_id;
        }
        (*BeadType)[btype].Charge = atof(split[3]);
        (*BeadType)[btype].Number++;
        (*Bead)[id].Position.x = atof(split[4]) - (*box_lo).x;
        (*Bead)[id].Position.y = atof(split[5]) - (*box_lo).y;
        (*Bead)[id].Position.z = atof(split[6]) - (*box_lo).z;
        (*Bead)[id].Type = btype;
        (*Index)[id] = id;
      } //}}}
      (*Counts).Unbonded = monomer;
      // go through possible molecules and remove 1-bead molecules //{{{
      int count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).BeadsCoor; i++) {
        if (mols[i] == 1) {
          int bead = bead_mols[i][0];
          (*Bead)[bead].Molecule = -1;
          (*Counts).Unbonded++;
        } else if (mols[i] > 1){
          (*Counts).Molecules++;
          (*Counts).Bonded += mols[i];
          for (int j = 0; j < mols[i]; j++) {
            int bead = bead_mols[i][j];
            (*Bead)[bead].Molecule = count;
          }
          count++;
        }
      } //}}}
      // TODO: is the 'remove 1-bead...' necessary? Join with 'allocate Molecule struct...'
      // remove single-bead molecules from mols array/{{{
      count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).BeadsCoor; i++) {
        if (mols[i] > 1) {
          mols[count] = mols[i];
          bead_mols[count] = realloc(bead_mols[count],
                                     sizeof *bead_mols[count] * mols[count]);
          for (int j = 0; j < mols[i]; j++) {
            bead_mols[count][j] = bead_mols[i][j];
          }
          // sort molecules in bead_mols[count][] according to ascending id
          SortArray(bead_mols[count], mols[i], 0);
          count++;
        }
      } //}}}
      // zeroize unused part of mols array - just to be on the save side //{{{
      for (int i = (*Counts).Molecules; i < (*Counts).BeadsCoor; i++) {
        mols[i] = 0;
      } //}}}
      // allocate Molecule struct and fill Molecule[].Bead array //{{{
      *Molecule = calloc((*Counts).Molecules, sizeof (MOLECULE));
      for (int i = 0; i < (*Counts).Molecules; i++) {
        (*Molecule)[i].Bead = malloc(sizeof *(*Molecule)[i].Bead * mols[i]);
        for (int j = 0; j < mols[i]; j++) {
          (*Molecule)[i].Bead[j] = bead_mols[i][j];
        }
      } //}}}
      // free helper array  //{{{
      for (int i = 0; i < (*Counts).BeadsCoor; i++) {
        free(bead_mols[i]);
      }
      free(bead_mols); //}}}
      free(mols);
    } //}}}
    // bonds section //{{{
    if (words > 0 && strcmp(split[0], "Bonds") == 0) {
      // allocate helper arrays to hold bond info //{{{
      // number of bonds in each molecule
      int *bonds_per_mol = calloc((*Counts).Molecules, sizeof *bonds_per_mol);
      // bond list for each molecule
      // [i][j][] ... bond id in the molecule 'i'
      // [i][][0] & [i][][1] ... ids of connected beads in molecule 'i'
      // [i][][2] ... bond type
      // TODO: this ain't right - some struct with (*array)[3] akin to connectivity
      int (**mol_bonds)[3] = calloc((*Counts).Molecules, sizeof (**mol_bonds)[3]);
      for (int i = 0; i < (*Counts).Molecules; i++) {
         mol_bonds[i] = calloc(1, sizeof(int *));
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // read all bonds //{{{
      for (int i = 0; i < *bonds; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // format of each line: <bond id> <bond type> <bead1> <bead2>
        // Error - incorrect format //{{{
        if (words < 4 ||
            !IsInteger_old(split[0]) || !IsInteger_old(split[1]) ||
            !IsInteger_old(split[2]) || !IsInteger_old(split[3])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each 'Bonds' line must be <bond id> <bond type> <bead1d> <bead2>\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - beads in ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "bond %d", atoi(split[0]));
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, "are in different molecules\n\n");
          ColourReset(STDERR_FILENO);
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        bonds_per_mol[mol]++;
        int bond = bonds_per_mol[mol];
        // add bonded beads to the molecule they belong to
        mol_bonds[mol] = realloc(mol_bonds[mol],
                                 sizeof *mol_bonds[mol] * bond);
        mol_bonds[mol][bond-1][0] = bead1;
        mol_bonds[mol][bond-1][1] = bead2;
        mol_bonds[mol][bond-1][2] = type;
      } //}}}
      // sort bonds according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortBonds(mol_bonds[i], bonds_per_mol[i]);
      } //}}}
      // minimize mol_bonds based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).BeadsCoor; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mol_bonds[i][j][0] -= lowest;
          mol_bonds[i][j][1] -= lowest;
        }
      } //}}}
      // identify molecule type based on bead order and connectivity //{{{
      *MoleculeType = malloc(sizeof *MoleculeType * 1);
      // number of molecule types differing in connectivity but not in beads - naming purposes
      int *diff_conn = malloc(sizeof *diff_conn * 1);
      // number of molecule types differing in bead order - naming purposes
      int mtype_bead_order = 0;
      for (int i = 0; i < (*Counts).Molecules; i++) {
        // is molecule 'i' of known type?
        bool exists = false;
        // do 'i' and given molecule type share connectivity and bead order?
        bool same_conn = true, same_bead = true;
        // molecule type with which 'i' shares bead order - naming purposes
        int type_bead_order = -1;
        for (int j = 0; j < (*Counts).TypesOfMolecules; j++) { //{{{
          if ((*MoleculeType)[j].nBeads == mols[i] && // same number of molecules?
              (*MoleculeType)[j].nBonds == bonds_per_mol[i]) { // same number of bonds?
            // same connectivity?
            same_conn = true;
            for (int k = 0; k < (*MoleculeType)[j].nBonds; k++) {
              if ((*MoleculeType)[j].Bond[k][0] != mol_bonds[i][k][0] ||
                  (*MoleculeType)[j].Bond[k][1] != mol_bonds[i][k][1] ||
                  (*MoleculeType)[j].Bond[k][2] != mol_bonds[i][k][2]) {
                same_conn = false;
                break;
              }
            }
            // same bead types?
            same_bead = true;
            for (int k = 0; k < (*MoleculeType)[j].nBeads; k++) {
              int btype = (*Bead)[(*Molecule)[i].Bead[k]].Type;
              if ((*MoleculeType)[j].Bead[k] != btype) {
                same_bead = false;
                break;
              }
            }
            if (same_bead && type_bead_order == -1) {
              type_bead_order = j;
            }
            // if the molecule has the same connectivity and bead order, it's of known type
            if (same_conn && same_bead) {
              exists = true;
              (*MoleculeType)[j].Number++;
              (*Molecule)[i].Type = j;
              break;
            }
          }
        } //}}}
        // add new type? //{{{
        if (!exists) {
          int mtype = (*Counts).TypesOfMolecules;
          *MoleculeType = realloc(*MoleculeType,
                                  sizeof (MOLECULETYPE) * (mtype + 1));
          diff_conn = realloc(diff_conn, sizeof *diff_conn * (mtype + 1));
          diff_conn[mtype] = 0;
          // molecule name
          if (same_bead && !same_conn) { // same beads - mol<type_bead_order>-b#
            diff_conn[type_bead_order]++;
            // shorten name if necessary to append '-b<int>'
            char name[MOL_NAME+1];
            strcpy(name, (*MoleculeType)[type_bead_order].Name);
            if (diff_conn[type_bead_order] < 10) {
              name[MOL_NAME-3] = '\0';
            } else if (diff_conn[type_bead_order] < 100) {
              name[MOL_NAME-4] = '\0';
            }
//          P_IGNORE(-Wformat-truncation);
            // MOL_NAME is max string length, i.e., array is longer
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "%s-b%d", name, diff_conn[type_bead_order]);
//          P_POP;
          } else { // same connectivity or both different - mol#
//          P_IGNORE(-Wformat-truncation);
            mtype_bead_order++;
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "mol%d", mtype_bead_order);
//          P_POP;
          }
          (*Molecule)[i].Type = mtype;
          (*MoleculeType)[mtype].Number = 1;
          // copy bead sequence and determine BTypes stuff
          (*MoleculeType)[mtype].nBeads = mols[i];
          (*MoleculeType)[mtype].Bead =
              malloc(sizeof *(*MoleculeType)[mtype].Bead *
                     (*MoleculeType)[mtype].nBeads);
          (*MoleculeType)[mtype].nBTypes = 0;
          (*MoleculeType)[mtype].BType =
              malloc(sizeof *(*MoleculeType)[mtype].BType * 1);
          for (int j = 0; j < (*MoleculeType)[mtype].nBeads; j++) {
            int btype = (*Bead)[(*Molecule)[i].Bead[j]].Type;
            (*MoleculeType)[mtype].Bead[j] = btype;
            exists = false; // recycling the bool to check if btype is already in BType[]
            for (int k = 0; k < (*MoleculeType)[mtype].nBTypes; k++) {
              if ((*MoleculeType)[mtype].BType[k] == btype) {
                exists = true;
              }
            }
            if (!exists) { // recycled exists
              int types = (*MoleculeType)[mtype].nBTypes;
              (*MoleculeType)[mtype].nBTypes++;
              (*MoleculeType)[mtype].BType =
                  realloc((*MoleculeType)[mtype].BType,
                          sizeof *(*MoleculeType)[mtype].BType * (types + 1));
              (*MoleculeType)[mtype].BType[types] = btype;
            }
          }
          // copy bonds
          (*MoleculeType)[mtype].nBonds = bonds_per_mol[i];
          (*MoleculeType)[mtype].Bond =
              malloc(sizeof *(*MoleculeType)[mtype].Bond *
                     (*MoleculeType)[mtype].nBonds);
          for (int j = 0; j < (*MoleculeType)[mtype].nBonds; j++) {
            (*MoleculeType)[mtype].Bond[j][0] = mol_bonds[i][j][0];
            (*MoleculeType)[mtype].Bond[j][1] = mol_bonds[i][j][1];
            (*MoleculeType)[mtype].Bond[j][2] = mol_bonds[i][j][2];
          }
          (*MoleculeType)[mtype].nAngles = 0;
          (*MoleculeType)[mtype].Angle =
              malloc(sizeof *(*MoleculeType)[mtype].Angle * 1);
          (*MoleculeType)[mtype].Write = true;
          (*Counts).TypesOfMolecules++;
        } //}}}
      } //}}}
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          free(mol_bonds[i][j]);
        }
        free(mol_bonds[i]);
      }
      free(mol_bonds);
      free(bonds_per_mol);
      free(diff_conn); //}}}
    } //}}}
    // angles section //{{{
    if (words > 0 && strcmp(split[0], "Angles") == 0) {
      // allocate helper arrays to hold angle info //{{{
      // number of angles in each molecule
      int *angles_per_mol = calloc((*Counts).Molecules, sizeof *angles_per_mol);
      // bond list for each molecule
      // temp[i].angle[j][] ... bond id in the molecule 'i'
      // temp[i].angle[][0] & [i][][1] & [i][][2] ... ids of connected beads in molecule 'i'
      // temp[i].angle[][3] ... angle type
      struct temp {
        int (*angle)[4];
      } *molec = calloc((*Counts).Molecules, sizeof *molec);
      for (int i = 0; i < (*Counts).Molecules; i++) {
         molec[i].angle = calloc(1, sizeof *molec[i].angle);
      } //}}}
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // read all angles //{{{
      for (int i = 0; i < *angles; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine_old(split, line, " \t");
        // format of each line: <angle id> <angle type> <bead1> <bead2> <bead3>
        // Error - incorrect format //{{{
        if (words < 5 ||
            !IsInteger_old(split[0]) || !IsInteger_old(split[1]) ||
            !IsInteger_old(split[2]) || !IsInteger_old(split[3]) || !IsInteger_old(split[4])) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - each 'Angles' line must be <angle id> <angle type> <bead1d> <bead2> <bead3>\n");
          ColourReset(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        int type = atoi(split[1]) - 1; // in lammps, bond types start at 1
        int bead1 = atoi(split[2]) - 1; // in lammps, atom ids start at 1
        int bead2 = atoi(split[3]) - 1;
        int bead3 = atoi(split[4]) - 1;
        // assign molecule to the bond
        int mol = (*Bead)[bead1].Molecule;
        // error when the second bead is in different molecule //{{{
        if (mol != (*Bead)[bead2].Molecule || mol != (*Bead)[bead3].Molecule) {
          ErrorPrintError_old();
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "%s", data_file);
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " - atoms in ");
          ColourChange(STDERR_FILENO, YELLOW);
          fprintf(stderr, "angle %d", atoi(split[0]));
          ColourChange(STDERR_FILENO, RED);
          fprintf(stderr, " are in different molecules\n\n");
          ColourReset(STDERR_FILENO);
          exit(1);
        } //}}}
        // increment number of bonds in the molecule
        angles_per_mol[mol]++;
        int num = angles_per_mol[mol];
        // add angle beads to the molecule they belong to
        molec[mol].angle = realloc(molec[mol].angle,
                                   sizeof *molec[mol].angle * num);
        molec[mol].angle[num-1][0] = bead1;
        molec[mol].angle[num-1][1] = bead2;
        molec[mol].angle[num-1][2] = bead3;
        molec[mol].angle[num-1][3] = type;
      } //}}}
      // sort angles according to the id of the first bead in a bond //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        SortAngles(molec[i].angle, angles_per_mol[i]);
      } //}}}
      // minimize mol_angles based on lowest id in each molecule //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int lowest = (*Counts).BeadsCoor; // just some high number
        for (int j = 0; j < mols[i]; j++) {
          if ((*Molecule)[i].Bead[j] < lowest) {
            lowest = (*Molecule)[i].Bead[j];
          }
        }
        for (int j = 0; j < angles_per_mol[i]; j++) {
          molec[i].angle[j][0] -= lowest;
          molec[i].angle[j][1] -= lowest;
          molec[i].angle[j][2] -= lowest;
        }
      } //}}}
      /*
       * add angles to existing molecule types (or create a new one differing only in angles)
       * ...'basic' molecule types already generated in bonds section
       */
      // number of molecule types differing in angles
      int *extra = calloc((*Counts).TypesOfMolecules, sizeof *extra);
      for (int i = 0; i < (*Counts).Molecules; i++) {
        int mtype = (*Molecule)[i].Type;
        if ((*MoleculeType)[mtype].nAngles == 0) { // add angles if there are no angles in the molecule type //{{{
          (*MoleculeType)[mtype].nAngles = angles_per_mol[i];
          (*MoleculeType)[mtype].Angle =
              realloc((*MoleculeType)[mtype].Angle,
                      sizeof *(*MoleculeType)[mtype].Angle *
                      (*MoleculeType)[mtype].nAngles);
          for (int j = 0; j < (*MoleculeType)[mtype].nAngles; j++) {
            (*MoleculeType)[mtype].Angle[j][0] = molec[i].angle[j][0];
            (*MoleculeType)[mtype].Angle[j][1] = molec[i].angle[j][1];
            (*MoleculeType)[mtype].Angle[j][2] = molec[i].angle[j][2];
            (*MoleculeType)[mtype].Angle[j][3] = molec[i].angle[j][3];
          } //}}}
        } else { // if angles already present, check whether they're the same //{{{
          bool exists = true;
          // same number of angles? //{{{
          if ((*MoleculeType)[mtype].nAngles != angles_per_mol[i]) {
            exists = false;
          } //}}}
          // same angles as in mtype? //{{{
          for (int j = 0; j < angles_per_mol[i] && j < (*MoleculeType)[mtype].nAngles; j++) {
            if ((*MoleculeType)[mtype].Angle[j][0] != molec[i].angle[j][0] ||
                (*MoleculeType)[mtype].Angle[j][1] != molec[i].angle[j][1] ||
                (*MoleculeType)[mtype].Angle[j][2] != molec[i].angle[j][2] ||
                (*MoleculeType)[mtype].Angle[j][3] != molec[i].angle[j][3] ) {
              exists = false;
            }
          } //}}}
          // check against angle-generated molecule types //{{{
          // If its angles aren't the same as in mtype, check against other
          // types (i.e., against newly generated thanks to different angles)
          if (!exists) {
            for (int j = (mtype+1); j < (*Counts).TypesOfMolecules; j++) {
              if ((*MoleculeType)[j].nAngles == angles_per_mol[i]) {
                int count = 0;
                for (int k = 0; k < (*MoleculeType)[j].nAngles; k++) {
                  if ((*MoleculeType)[j].Angle[k][0] == molec[i].angle[k][0] &&
                      (*MoleculeType)[j].Angle[k][1] == molec[i].angle[k][1] &&
                      (*MoleculeType)[j].Angle[k][2] == molec[i].angle[k][2] &&
                      (*MoleculeType)[j].Angle[k][3] == molec[i].angle[k][3] ) {
                    count++;
                  }
                }
                if (count == angles_per_mol[i]) {
                  exists = true;
                  (*MoleculeType)[j].Number++;
                  (*MoleculeType)[mtype].Number--;
                  (*Molecule)[i].Type = j;
                  break;
                }
              }
            }
          } //}}}
          // create a new molecule type //{{{
          if (!exists) {
            int new = (*Counts).TypesOfMolecules;
            extra[mtype]++;
            (*Molecule)[i].Type = new;
            (*Counts).TypesOfMolecules++;
            *MoleculeType = realloc(*MoleculeType,
                                    sizeof (MOLECULETYPE) * (new + 1));
            // shorten name if necessary to append '-a<int>'
            char name[MOL_NAME+1];
            strcpy(name, (*MoleculeType)[mtype].Name);
            if (extra[mtype] < 10) {
              name[MOL_NAME-3] = '\0';
            } else if (extra[mtype] < 100) {
              name[MOL_NAME-4] = '\0';
            }
//          P_IGNORE(-Wformat-truncation);
            // MOL_NAME is max string length, i.e., array is longer
            snprintf((*MoleculeType)[mtype].Name, MOL_NAME+1, "%s-a%d",
                      name, extra[mtype]);
//          P_POP;
            (*MoleculeType)[new].Number = 1;
            (*MoleculeType)[mtype].Number--;
            (*MoleculeType)[new].nBeads = (*MoleculeType)[mtype].nBeads;
            (*MoleculeType)[new].Bead =
                malloc(sizeof *(*MoleculeType)[new].Bead *
                       (*MoleculeType)[new].nBeads);
            for (int j = 0; j < (*MoleculeType)[new].nBeads; j++) {
              (*MoleculeType)[new].Bead[j] = (*MoleculeType)[mtype].Bead[j];
            }
            (*MoleculeType)[new].nBonds = (*MoleculeType)[mtype].nBonds;
            (*MoleculeType)[new].Bond =
                malloc(sizeof *(*MoleculeType)[new].Bond *
                       (*MoleculeType)[new].nBonds);
            for (int j = 0; j < (*MoleculeType)[new].nBonds; j++) {
              (*MoleculeType)[new].Bond[j][0] = (*MoleculeType)[mtype].Bond[j][0];
              (*MoleculeType)[new].Bond[j][1] = (*MoleculeType)[mtype].Bond[j][1];
              (*MoleculeType)[new].Bond[j][2] = (*MoleculeType)[mtype].Bond[j][2];
            }
            (*MoleculeType)[new].nAngles = (*MoleculeType)[mtype].nAngles;
            (*MoleculeType)[new].Angle =
                malloc(sizeof *(*MoleculeType)[new].Angle *
                       (*MoleculeType)[new].nAngles);
            for (int j = 0; j < (*MoleculeType)[new].nAngles; j++) {
              (*MoleculeType)[new].Angle[j][0] = molec[i].angle[j][0];
              (*MoleculeType)[new].Angle[j][1] = molec[i].angle[j][1];
              (*MoleculeType)[new].Angle[j][2] = molec[i].angle[j][2];
              (*MoleculeType)[new].Angle[j][3] = molec[i].angle[j][3];
            }
            (*MoleculeType)[new].nBTypes = (*MoleculeType)[mtype].nBTypes;
            (*MoleculeType)[new].BType =
                malloc(sizeof *(*MoleculeType)[new].BType *
                       (*MoleculeType)[new].nBTypes);
            for (int j = 0; j < (*MoleculeType)[new].nBTypes; j++) {
              (*MoleculeType)[new].BType[j] = (*MoleculeType)[mtype].BType[j];
            }
            (*MoleculeType)[new].Mass = (*MoleculeType)[mtype].Mass;
            (*MoleculeType)[new].Write = true;
          } //}}}
        } //}}}
      }
      // free helper arrays //{{{
      for (int i = 0; i < (*Counts).Molecules; i++) {
        free(molec[i].angle);
      }
      free(molec);
      free(angles_per_mol);
      free(extra); //}}}
    } //}}}
    // read and split next line
    fgets(line, sizeof line, fr);
    words = SplitLine_old(split, line, " \t");
  } //}}}
  free(mols);

  // calculate molecular mass //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      (*MoleculeType)[i].Mass += (*BeadType)[(*MoleculeType)[i].Bead[j]].Mass;
    }
  } //}}}

  WarnElNeutrality(*Counts, *BeadType, data_file);

  fclose(fr);
} //}}}
// TODO will be changed - agg files
// TODO restructure aggregates - don't include in BEAD
// ReadAggregates() //{{{
/*
 * Function reading information about aggregates from agg file generated by
 * Aggregates utility.
 */
void ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts,
                    AGGREGATE *Aggregate[], BEADTYPE *BeadType, BEAD *Bead[],
                    MOLECULETYPE *MoleculeType, MOLECULE *Molecule[],
                    int *Index) {
  // TODO counting lines
  int file_line_count = -1;
  char line[LINE], *split[SPL_STR];
  // read 'Step|Last Step' line
  fgets(line, sizeof line, fr);
  int words = SplitLine(SPL_STR, split, line, " \t\n");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  if (split[0][0] == 'L') {
    strcpy(ERROR_MSG, "premature end of file");
    PrintError();
    ErrorPrintFile(agg_file, "\0");
    putc('\n', stderr);
    exit(1);
  } else if (words < 2 || strcmp(split[0], "Step:") != 0 ||
             !IsInteger_old(split[1])) {
//  ErrorPrintError_old();
//  ColourChange(STDERR_FILENO, YELLOW);
//  fprintf(stderr, "%s", agg_file);
//  ColourChange(STDERR_FILENO, RED);
//  fprintf(stderr, " - wrong 'Step' line\n");
//  ColourReset(STDERR_FILENO);
//  ErrorPrintLine(split, words);

    strcpy(ERROR_MSG, "missing pbc line (and any coordinates)");
    PrintError();
    PrintErrorFileLine(agg_file, "\0", file_line_count, split, words);
    putc('\n', stderr);
    exit(1);
  } //}}}
  // initialize array of number of aggregates per bead //{{{
  for (int i = 0; i < (*Counts).BeadsCoor; i++) {
    (*Bead)[i].nAggregates = 0;
  } //}}}
  // get number of aggregates //{{{
  fgets(line, sizeof line, fr);
  words = SplitLine(SPL_STR, split, line, " \t\n");
  // error - the number of aggregates must be <int>
  if (words == 0 || !IsInteger_old(split[0])) {
//  ErrorPrintError_old();
//  ColourChange(STDERR_FILENO, YELLOW);
//  fprintf(stderr, "%s", agg_file);
//  ColourChange(STDERR_FILENO, RED);
//  fprintf(stderr, " - number of aggregates must be a whole number\n");
//  ColourReset(STDERR_FILENO);
//  ErrorPrintLine(split, words);
    exit(1);
  }
  (*Counts).Aggregates = atoi(split[0]); //}}}
  // skip blank line - error if not a blank line //{{{
  fgets(line, sizeof line, fr);
  words = SplitLine(SPL_STR, split, line, " \t\n");
  if (words > 0) {
//  ErrorPrintError_old();
//  ColourChange(STDERR_FILENO, YELLOW);
//  fprintf(stderr, "%s", agg_file);
//  ColourChange(STDERR_FILENO, RED);
//  fprintf(stderr, " - missing a blank line after number of aggregates\n");
//  ColourReset(STDERR_FILENO);
//  ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // go through all aggregates
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    // TODO: error catching - use fgets() & SafeStrcat() & SplitLine()
    // read molecules in Aggregate 'i' //{{{
    fscanf(fr, "%d :", &(*Aggregate)[i].nMolecules);
    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol;
      fscanf(fr, "%d", &mol);
      mol--; // in agg file, the numbers correspond to vmd
      (*Aggregate)[i].Molecule[j] = mol;
      (*Molecule)[mol].Aggregate = i;
    }
    while (getc(fr) != '\n')
     ; //}}}
    // read monomeric beads in Aggregate 'i' //{{{
    int count;
    fscanf(fr, "%d :", &count);
    // their number will be counted according to which beads are in vcf
    (*Aggregate)[i].nMonomers = 0;
    for (int j = 0; j < count; j++) {
      int vsf_id;
      fscanf(fr, "%d", &vsf_id); // monomer index in vsf file
      int id = Index[vsf_id]; // monomer index in Bead structure
      if (id > -1) {
        int beads = (*Aggregate)[i].nMonomers++;
        (*Aggregate)[i].Monomer[beads] = id;
//      (*Bead)[id].nAggregates++;
//      (*Bead)[id].Aggregatexxx = realloc((*Bead)[id].Aggregatexxx,
//                                      sizeof *(*Bead)[id].Aggregatexxx *
//                                      (*Bead)[id].nAggregates);
//      (*Bead)[id].Aggregatexxx[(*Bead)[id].nAggregates-1] = i;
      }
    }
    while (getc(fr) != '\n')
     ; //}}}
  }
  // skip blank line at the end of every entry
  while (getc(fr) != '\n')
    ;
  // fill Aggregate[].Bead, Bead[].Aggregate, and Molecule[].Aggregate arrays //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    (*Aggregate)[i].nBeads = 0;

    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int mol = (*Aggregate)[i].Molecule[j];
      (*Molecule)[mol].Aggregate = i;
      for (int k = 0; k < MoleculeType[(*Molecule)[mol].Type].nBeads; k++) {
        int bead = (*Molecule)[mol].Bead[k];
        (*Aggregate)[i].Bead[(*Aggregate)[i].nBeads+k] = bead;
//      (*Bead)[bead].Aggregatexxx[(*Bead)[bead].nAggregates] = i;
        (*Bead)[bead].nAggregates++;
      }

      // increment number of beads in aggregate 'i'
      (*Aggregate)[i].nBeads += MoleculeType[(*Molecule)[mol].Type].nBeads;
    }
  }

  // fill Bead[].Aggregate for monomeric Beads
//  for (int i = 0, i < (*Counts).Bead)
  //}}}
  // calculate aggregates' masses //{{{
  for (int i = 0; i < (*Counts).Aggregates; i++) {
    (*Aggregate)[i].Mass = 0;

    for (int j = 0; j < (*Aggregate)[i].nMolecules; j++) {
      int type = (*Molecule)[(*Aggregate)[i].Molecule[j]].Type;

      (*Aggregate)[i].Mass += MoleculeType[type].Mass;
    }
  } //}}}
} //}}}
// ReadAggCommand() //{{{
/*
 * Function to read the Aggregate command from an agg file. The command must be
 * on the first line.
 */
void ReadAggCommand(BEADTYPE *BeadType, COUNTS Counts,
                    char *input_coor, char *input_agg,
                    double *distance, int *contacts) {
  // open input aggregate file
  FILE *agg = OpenFile(input_agg, "r");
  // read first line (Aggregate command)
  char line[LINE];
  fgets(line, sizeof line, agg);
  if (!ReadLine(agg, LINE, line)) {
    strcpy(ERROR_MSG, "cannot read the first line from an aggregate file");
    PrintError();
    ErrorPrintFile(input_agg, "\0");
    putc('\n', stderr);
    exit(1);
  }
  fclose(agg);
  char *split[SPL_STR];
  int words = SplitLine(SPL_STR, split, line, " \t\n");
  // error - not enough strings for a proper Aggregate command //{{{
  if (words < 6) {
    strcpy(ERROR_MSG, "first line must contain a valid Aggregates command");
    PrintError();
    ErrorPrintFile(input_agg, "\0");
    putc('\n', stderr);
    exit(1);
  } //}}}
  // read <distance> argument from Aggregates command //{{{
  if (!IsPosReal(split[2], distance)) {
    strcpy(ERROR_MSG, "Aggregate command: <distance> must be \
a non-negative number");
    PrintError();
    ErrorPrintFile(input_agg, "\0");
    ErrorPrintLine2(split, words);
    exit(1);
  } //}}}
  // read <contacts> argument from Aggregates command //{{{
  long val;
  if (!IsPosInteger(split[3], &val)) {
    strcpy(ERROR_MSG, "Aggregate command: <contacts> must be \
a positive integer");
    PrintError();
    ErrorPrintFile(input_agg, "\0");
    ErrorPrintLine2(split, words);
    exit(1);
  }
  *contacts = val; //}}}
  // warning - differently named vcf file than the one in agg file //{{{
  if (strcmp(split[1], input_coor) != 0) {
    strcpy(ERROR_MSG, "different coordinate file in the Aggregate command to \
the one used here can lead to undefined behaviour in case of mismatch between \
beads present in the two files");
    PrintWarning();
  } //}}}
  // read <type names> from Aggregates command //{{{
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType_old(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      strcpy(ERROR_MSG, "non-existent bead name in the Aggregaates command");
      PrintError();
      ErrorPrintFile(input_agg, "\0");
      fprintf(stderr, "%s, bead %s%s%s\n", Red(), Yellow(), split[i],
                                           ColourReset());
      ErrorBeadType_old(Counts, BeadType);
      exit(1);
    }
  } //}}}
} //}}}
// SkipAgg() //{{{
/*
 * Function to skip one timestep in aggregate file.
 */
/*
 * TODO: all agg file stuff - consider what & how to read it; do I want to
 * exit(1) if there's fewer steps in agg than in coor file? That seems to be
 * what's happening now...
 */
void SkipAgg(FILE *agg, char *agg_file) {
  int file_line_count = -1; // TODO line counting
  char line[LINE];
  fgets(line, sizeof line, agg);
  if (!ReadLine(agg, LINE, line)) {
    strcpy(ERROR_MSG, "cannot read a line from the aggregate file");
    PrintError();
    ErrorPrintFile(agg_file, "\0");
    putc('\n', stderr);
    exit(1);
  }
  char *split[SPL_STR];
  int words = SplitLine(SPL_STR, split, line, " \t\n");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  long val;
  if (split[0][0] == 'L') {
    strcpy(ERROR_MSG, "premature end of file");
    ErrorPrintFile(agg_file, "\0");
    putc('\n', stderr);
    exit(1);
  } else if (words < 2 || !IsNatural(split[1], &val)) {
    strcpy(ERROR_MSG, "invalid 'Step' line");
    PrintErrorFileLine(agg_file, "\0", file_line_count, split, words);
    exit(1);
  } //}}}
  // get number of aggregates
  fgets(line, sizeof line, agg);
  if (!ReadLine(agg, LINE, line)) {
    strcpy(ERROR_MSG, "cannot read a line from the aggregate file");
    PrintError();
    ErrorPrintFile(agg_file, "\0");
    putc('\n', stderr);
    exit(1);
  }
  words = SplitLine(SPL_STR, split, line, " \t\n");
  // Error - number of aggregates must be <int> //{{{
  if (words != 0 && !IsPosInteger(split[0], &val)) {
    strcpy(ERROR_MSG, "number of aggregates must be a natural number\n");
    PrintErrorFileLine(agg_file, "\0", file_line_count, split, words);
    exit(1);
  } //}}}
  int aggs = atoi(split[0]);
  // skip the blank line and the aggregate lines
  for (int i = 0; i < (1+2*aggs); i++) {
    int test;
    while ((test = getc(agg)) != '\n') {
      if (test == EOF) {
        strcpy(ERROR_MSG, "premature end of file");
        ErrorPrintFile(agg_file, "\0");
        fprintf(stderr, "%s, line %s%d%s\n", ErrRed(), ErrYellow(),
                                             file_line_count, ErrColourReset());
        exit(1);
      }
    }
  }
  // skip empty line at the end
  fgets(line, sizeof line, agg);
  if (!ReadLine(agg, LINE, line)) {
    strcpy(ERROR_MSG, "cannot read a line from the aggregate file");
    PrintError();
    ErrorPrintFile(agg_file, "\0");
    putc('\n', stderr);
    exit(1);
  }
  words = SplitLine(SPL_STR, split, line, " \t\n");
  if (feof(agg) == EOF) {
    strcpy(ERROR_MSG, "premature end of file");
    ErrorPrintFile(agg_file, "\0");
    exit(1);
  }
} //}}}
// TODO remove
// VtfReadTimestep_old() //{{{
/*
 * Read a single timestep from a vcf/vtf coordinate line - first the preamble
 * (getting timestep type and possibly pbc, ignoring any unrecognised line),
 * then the coordinates (if a wrong line is in the coordinate section, the
 * function skips the remaining coordinate lines and returns to the beginning
 * to read the next timestep). Returns false only if a line couldn't be read
 * (e.g., on eof).
 */
bool VtfReadTimestep_old(FILE *vcf, char vcf_file[], BOX *Box, COUNTS *Counts,
                     BEADTYPE BeadType[], BEAD *Bead[], int Index[],
                     MOLECULETYPE MoleculeType[], MOLECULE Molecule[],
                     int *InFile[], int *file_line_count,
                     int step_count, char stuff[]) {
  start_function: ; // return here when a bad line is encountered
  char *cur = stuff, * const end = stuff + LINE; // to properly snprintf stuff[]
  // set 'not in timestep' to all beads //{{{
  for (int i = 0; i < (*Counts).BeadsTotal; i++) {
    (*Bead)[i].InTimestep = false;
  } //}}}
  // read timestep preamble //{{{
  (*Counts).BeadsCoor = -1; // no coordinate line found yet
  int ltype, timestep = ERROR_LINE;
  fpos_t position; // to save file position
  while (true) {
    // save file pointer position for when it's the first coordinate line
    fgetpos(vcf, &position);
    (*file_line_count)++;
    // read line & split it via whitespace //{{{
    char *split[SPL_STR], line[LINE];
    int words;
    if (!ReadAndSplitLine(vcf, LINE, line, &words, split, SPL_STR, " \t\n")) {
      return false;
    } //}}}
    ltype = VtfCheckLineType(words, split, vcf_file, *file_line_count);
    // do something based on what line it is
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) { //{{{
      (*Box).Length.x = atof(split[1]);
      (*Box).Length.y = atof(split[2]);
      (*Box).Length.z = atof(split[3]);
      if (ltype == PBC_LINE_ANGLES) {
        (*Box).alpha = atof(split[4]);
        (*Box).beta = atof(split[5]);
        (*Box).gamma = atof(split[6]);
        if (!TriclinicCellData(Box)) {
          ErrorPrintFull2(vcf_file, *file_line_count, split, words);
          exit(1);
        }
      } else {
        (*Box).alpha = 90;
        (*Box).beta = 90;
        (*Box).gamma = 90;
      }
     //}}}
    } else if (ltype == TIME_LINE_I || ltype == TIME_LINE_O) { //{{{
      // warn: multiple timestep lines
      if (timestep != ERROR_LINE) {
        strcpy(ERROR_MSG, "extra timestep line \
(or timestep without any coordinates)");
        PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
      }
      timestep = ltype;
     //}}}
    } else if (ltype == COOR_LINE_I || ltype == COOR_LINE_O) { //{{{
      // warn: missing timestep line - read next timeste //{{{
      if (timestep == ERROR_LINE) {
        strcpy(ERROR_MSG, "found a coordinate line before a timestep line; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      (*Counts).BeadsCoor = 0;
      break; //}}}
    } else if (ltype == COMMENT_LINE) {
      // TODO clean this; well, do something about it, anyway
//    printf(">|%s|<\n", stuff);
//    printf("%ld %ld\n", sizeof stuff, strlen(stuff));
//    if ((end-cur) < LINE) {
      if ((end-cur) > 0) {
        cur += snprintf(cur, end-cur, "%s", split[0]);
      }
        for (int i = 1; i < (words-1) && cur < end; i++) {
          if ((end-cur) > 0) {
            cur += snprintf(cur, end-cur, " %s", split[i]);
          } else {
            break;
          }
        }
//      printf("%ld %ld %d |%s|\n", end-cur, strlen(stuff), LINE, stuff);
        if ((end-cur) > 0) {
          cur += snprintf(cur, end-cur, " %s\n", split[words-1]);
        }
//      printf("x %ld %ld %d |%s|\n", end-cur, strlen(stuff), LINE, stuff);
//    } else {
//      printf("%ld %d |%s|\n", end-cur, LINE, stuff);
//    }
//    printf("|%s| %ld\n", stuff, strlen(stuff));
    } else if (ltype == ERROR_LINE) {
      strcpy(ERROR_MSG, "ignoring unrecognised line in a timestep preamble");
      PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
    }
  } //}}}
  // return 'false' if no coordinate line encountered //{{{
  if ((*Counts).BeadsCoor == -1) {
    strcpy(ERROR_MSG, "no coordinates; this warning should never trigger!");
    PrintWarning();
  } //}}}
  // read coordinates  //{{{
  (*Counts).UnbondedCoor = 0;
  (*Counts).BondedCoor = 0;
  // restore file pointer to before the first coordinate line
  fsetpos(vcf, &position);
  (*file_line_count)--; // the first coordinate line will be re-read
  while (true) {
    // save file pointer position for when it's the first coordinate line
    fgetpos(vcf, &position);
    (*file_line_count)++;
    // read line & split it via whitespace //{{{
    char line[LINE], *split[SPL_STR];
    int words;
    if (!ReadAndSplitLine(vcf, LINE, line, &words, split, SPL_STR, " \t\n")) {
      return true;
    } //}}}
    ltype = VtfCheckLineType(words, split, vcf_file, *file_line_count);
    if (ltype != COOR_LINE_O && ltype != COOR_LINE_I) {
      // warn: unrecognised line - read next timestep //{{{
      if (ltype == ERROR_LINE) {
        strcpy(ERROR_MSG, "unrecognised line in a timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      break;
    }
    int id;
    if (timestep == TIME_LINE_I) { // 'timestep indexed' coordinate line
      if (ltype == COOR_LINE_I) {
        id = atoi(split[0]);
        // warn: bead index is too high - read next timestep //{{{
        if (id >= (*Counts).BeadsTotal) {
          strcpy(ERROR_MSG, "bead index too high; \
using next timestep instead of this one");
          PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
          // skip the remaing coordinate line (if any)
          do {
            fgetpos(vcf, &position);
            (*file_line_count)++;
          } while (VtfSkipCoorOrderedLine(vcf));
          fsetpos(vcf, &position);
          (*file_line_count)--; // the first non-coordinate line will be re-read
          goto start_function;
        } //}}}
        id = Index[id];
        // warn: repeated bead index - read next timestep //{{{
        if ((*Bead)[id].InTimestep) {
          strcpy(ERROR_MSG, "multiple bead entries with the same index; \
using next timestep instead of this one");
          PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
          // skip the remaing coordinate line (if any)
          do {
            fgetpos(vcf, &position);
            (*file_line_count)++;
          } while (VtfSkipCoorOrderedLine(vcf));
          fsetpos(vcf, &position);
          (*file_line_count)--; // the first non-coordinate line will be re-read
          goto start_function;
        } //}}}
        (*Bead)[id].Position.x = atof(split[1]);
        (*Bead)[id].Position.y = atof(split[2]);
        (*Bead)[id].Position.z = atof(split[3]);
        (*Bead)[id].InTimestep = true;
        VECTOR vel;
        if (words >= 7 && IsReal(split[4], &vel.x) &&
                          IsReal(split[5], &vel.y) &&
                          IsReal(split[6], &vel.z)) { // bead velocities, if there
          (*Bead)[id].Velocity.x = vel.x;
          (*Bead)[id].Velocity.y = vel.y;
          (*Bead)[id].Velocity.z = vel.z;
        }
        (*InFile)[(*Counts).BeadsCoor] = id;
      } else {
        // warn: ordered line in indexed timestep - read next timestep //{{{
        strcpy(ERROR_MSG, "ordered coordinate line in indexed timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function; //}}}
      }
    } else { // 'timestep ordered' coordinate line
      // warn: extra bead in ordered timestep - read next timestep //{{{
      if ((*Counts).BeadsCoor == (*Counts).BeadsTotal) {
        strcpy(ERROR_MSG, "too many beads in an ordered timestep; \
using next timestep instead of this one");
        PrintWarningFileLine(vcf_file, "\0", *file_line_count, split, words);
        // skip the remaing coordinate line (if any)
        do {
          fgetpos(vcf, &position);
          (*file_line_count)++;
        } while (VtfSkipCoorOrderedLine(vcf));
        fsetpos(vcf, &position);
        (*file_line_count)--; // the first non-coordinate line will be re-read
        goto start_function;
      } //}}}
      id = Index[(*Counts).BeadsCoor];
      (*Bead)[id].Position.x = atof(split[0]);
      (*Bead)[id].Position.y = atof(split[1]);
      (*Bead)[id].Position.z = atof(split[2]);
      (*Bead)[id].InTimestep = true;
      VECTOR val;
      if (words >= 6 && IsReal(split[3], &val.x) &&
          IsReal(split[4], &val.y) && IsReal(split[5], &val.z)) { // bead velocities, if present
        (*Bead)[id].Velocity.x = val.x;
        (*Bead)[id].Velocity.y = val.y;
        (*Bead)[id].Velocity.z = val.z;
      }
      (*InFile)[(*Counts).BeadsCoor] = id;
    }
    if ((*Bead)[id].Molecule == -1) {
      (*Counts).UnbondedCoor++;
    } else {
      (*Counts).BondedCoor++;
    }
    (*Counts).BeadsCoor++;
  }
  // restore file pointer to before the first non-coordinate line
  fsetpos(vcf, &position); //}}}
  // warn: too few beads in an ordered timestep - read next timestep //{{{
  if (timestep == TIME_LINE_O && (*Counts).BeadsCoor < (*Counts).BeadsTotal) {
    strcpy(ERROR_MSG, "insufficient number of beads for ordered timestep; \
using next timestep instead of this one");
    PrintWarning();
    WarnPrintFile(vcf_file, "\0");
    fprintf(stderr, "%s, last line of the timestep: %s%d%s\n",
            ErrCyan(), ErrYellow(), *file_line_count, ErrColourReset());
    goto start_function;
  } //}}}
  (*file_line_count)--; // the last line will be re-read next time
  return true; // coordinates read properly
} //}}}
// VtfReadStruct_old() //{{{
/*
 * Read system information from vsf/vtf structure file. It can recognize bead
 * and molecule types based either on name only or on all information (name,
 * mass, charge, and radius for bead types; bead order, bonds, angles, and
 * dihedrals for molecule types).
 */
// TODO check that there are any molecules before filling the structures
SYSTEM VtfReadStruct_old(char struct_file[], bool detailed) {
  SYSTEM Sys;
  InitSystem(&Sys);
  FILE *vsf = OpenFile(struct_file, "r");
  // define variables and structures and arrays //{{{
  int file_line_count = 0, // total number of lines
      count_atoms = 0, // number of 'atom <id>'
      default_atom = 0, // line number of the first 'atom default' line
      atom_names = 0, res_names = 0, // number of unieque bead & molecule names
      count_bonds = 0; // number of bonds
  bool warned = false; // has 'a[tom] default' line warning already been issued?
  int *mol_id = calloc(1, sizeof *mol_id); // to know which resids are used
  mol_id[0] = -1; // initialize by assuming 'resid 0' isn't in struct_file
  char (*atom_name)[BEAD_NAME+1] = calloc(1, sizeof *atom_name),
       (*res_name)[MOL_NAME+1] = calloc(1, sizeof *res_name);
  struct atom {
    int name, resid, resname;
    double charge, mass, radius;
  } atom_def, *atom = calloc(1, sizeof *atom);
  struct bond {
    int index1, index2; // indices from vsf file
  } *bond = calloc(1, sizeof *bond); //}}}
  // 1) read struct_file line by line, saving all atom and bond lines //{{{
  /* Do something based on to line type
   *   a) atom line: save bead information
   *   b) bond line: save bonded bead indices
   *   c) timestep line: break the loop (end of structure part of vsf file)
   *   d) coordinate line: exit program as coordinates cannot be inside
   *      structure file
   *   e) anything else besides pbc, blank, or comment line: exit program as
   *      unrecognised line was encountered
   */
  char line[LINE], *split[SPL_STR];
  int words;
  while (ReadAndSplitLine(vsf, LINE, line, &words, split, SPL_STR, " \t\n")) {
    file_line_count++;
    // read line
    file_line_count++;
    int ltype = VtfCheckLineType(words, split, struct_file, file_line_count);
    if (ltype == ATOM_LINE) { // a)
      if (strcmp(split[1], "default") == 0) { // 'a[tom] default' line
        // warning - multiple 'atom default' lines (warn only once)
        if (default_atom != 0 && !warned) { //{{{
          warned = true;
          strcpy(ERROR_MSG, "multiple 'a[tom] default' lines");
          PrintWarning();
          WarnPrintFile(struct_file, "\0");
          fprintf(stderr, "%s, using line %s%d%s as the default line%s\n",
                  ErrCyan(), ErrYellow(), default_atom, ErrCyan(),
                  ErrColourReset());
          //}}}
        } else { // save line number of the first 'atom default' line //{{{
          default_atom = file_line_count;
          // save values for the default bead type
          int *values = VtfAtomLineValues(words, split);
          // save name if unique //{{{
          int type = -1; // assume the name is yet unknown
          for (int i = 0; i < atom_names; i++) {
            if (strncmp(split[values[0]], atom_name[i], BEAD_NAME) == 0) {
              type = i;
              break;
            }
          }
          if (type == -1) { // unknown (i.e., new) bead name
            atom_name = realloc(atom_name, sizeof *atom_name *
                                           (atom_names + 1));
            strncpy(atom_name[atom_names], split[values[0]], BEAD_NAME);
            type = atom_names;
            atom_names++;
          }
          atom_def.name = type; //}}}
          // save mass //{{{
          if (values[1] != -1) {
            atom_def.mass = atof(split[values[1]]);
          } else {
            atom_def.mass = MASS;
          } //}}}
          // save charge //{{{
          if (values[2] != -1) {
            atom_def.charge = atof(split[values[2]]);
          } else {
            atom_def.charge = CHARGE;
          } //}}}
          // save radius //{{{
          if (values[3] != -1) {
            atom_def.radius = atof(split[values[3]]);
          } else {
            atom_def.radius = RADIUS;
          } //}}}
          atom_def.resname = -1;
          atom_def.resid = -1;
        } //}}}
      } else { // 'a[tom] <id>' line
        count_atoms++;
        int id = atoi(split[1]);
        // highest bead index? (corresponds to the number of beads in vsf) //{{{
        if (id >= Sys.Count.Bead) {
          atom = realloc(atom, sizeof *atom * (id + 1));
          // assign 'possible default bead' to all ids between old and new
          for (int i = Sys.Count.Bead; i < id; i++) {
            atom[i].name = -1;
          }
          Sys.Count.Bead = id + 1; // +1 as bead ids start from 0 in vsf
        } //}}}
        // save values from the 'a[tom] <id>' line
        int *values = VtfAtomLineValues(words, split);
        // save bead name if unique //{{{
        int type = -1; // assume the name is yet unknown
        for (int i = 0; i < atom_names; i++) {
          if (strncmp(split[values[0]], atom_name[i], BEAD_NAME) == 0) {
            type = i;
            break;
          }
        }
        if (type == -1) { // unknown (i.e., new) bead name
          atom_name = realloc(atom_name, sizeof *atom_name * (atom_names + 1));
          strncpy(atom_name[atom_names], split[values[0]], BEAD_NAME);
          type = atom_names;
          atom_names++;
        }
        atom[id].name = type; //}}}
        // save mass //{{{
        if (values[1] != -1) {
          atom[id].mass = atof(split[values[1]]);
        } else {
          atom[id].mass = MASS;
        } //}}}
        // save charge //{{{
        if (values[2] != -1) {
          atom[id].charge = atof(split[values[2]]);
        } else {
          atom[id].charge = CHARGE;
        } //}}}
        // save radius //{{{
        if (values[3] != -1) {
          atom[id].radius = atof(split[values[3]]);
        } else {
          atom[id].radius = RADIUS;
        } //}}}
        // save resid (-1 if the bead isn't in a molecule) //{{{
        if (values[5] > -1 ) {
          atom[id].resid = atoi(split[values[5]]);
          // save molecule name if unique
          if (values[4] != -1) {
            type = -1; // assume the name is yet unknown
            for (int i = 0; i < res_names; i++) {
              if (strncmp(split[values[4]], res_name[i], MOL_NAME) == 0) {
                type = i;
                break;
              }
            }
            if (type == -1) { // unknown (i.e., new) bead name
              res_name = realloc(res_name, sizeof *res_name * (res_names + 2));
              strncpy(res_name[res_names], split[values[4]], MOL_NAME);
              type = res_names;
              res_names++;
            }
            atom[id].resname = type;
          }
          // highest molecule id?
          if (atom[id].resid > Sys.Count.HighestResid) {
            mol_id = realloc(mol_id, sizeof *mol_id * (atom[id].resid + 1));
            for (int i = (Sys.Count.HighestResid+1); i <= atom[id].resid; i++) {
              mol_id[i] = -1;
            }
            Sys.Count.HighestResid = atom[id].resid;
          }
          if (mol_id[atom[id].resid] == -1) {
            mol_id[atom[id].resid] = Sys.Count.Molecule;
            Sys.Count.Molecule++;
          }
        } else {
          atom[id].resid = -1;
        } //}}}
      }
    } else if (ltype == BOND_LINE) { // b)
      bond = realloc(bond, sizeof *atom * (count_bonds + 1));
      long val;
      if (words == 2) { // case 'bond <id>:<id>'
        char *index[SPL_STR];
        SplitLine(SPL_STR, index, split[1], ":");
        IsInteger(index[0], &val);
        bond[count_bonds].index1 = val;
        IsInteger(index[1], &val);
        bond[count_bonds].index2 = val;
      } else { // case 'bond <id>: <id>'
        IsInteger(split[1], &val);
        bond[count_bonds].index1 = val;
        IsInteger(split[2], &val);
        bond[count_bonds].index2 = val;
      }
      count_bonds++;
    } else if (ltype == TIME_LINE_I || ltype == TIME_LINE_O) { // c)
      break;
    } else if (ltype == COOR_LINE_I || ltype == COOR_LINE_O) { // d)
      strcpy(ERROR_MSG, "encountered a coordinate-like line \
inside the structure block ");
      PrintErrorFileLine(struct_file, "\0", file_line_count, split, words);
      exit(1);
    } else if (ltype != BLANK_LINE && ltype != COMMENT_LINE &&
               ltype != PBC_LINE && ltype != PBC_LINE_ANGLES) { // e)
      // proper error message already established in VtfCheckLineType()
      PrintErrorFileLine(struct_file, "\0", file_line_count, split, words);
      exit(1);
    }
  }
  fclose(vsf); //}}}
  // error - no default line and too few atom lines //{{{
  if (default_atom == 0 && count_atoms != Sys.Count.Bead) {
    strcpy(ERROR_MSG, "not all beads defined ('atom default' line is omitted)");
    PrintError();
    ErrorPrintFile(struct_file, "\0");
    int undefined = Sys.Count.Bead-count_atoms;
    fprintf(stderr, "%s, %s%d%s bead(s) undefined%s\n", ErrRed(), ErrYellow(),
                                                        undefined, ErrRed(),
                                                        ErrColourReset());
    exit(1);
  } //}}}
  // create array linking internal and vsf resids //{{{
  /*
   * mol_id[vsf resid] = internal resid
   * mol_id_internal[internal resid] = vsf resid
   * therefore, mol_id_internal[mol_id[vsf resid]] = vsf resid
   */
  int *mol_id_internal = calloc(Sys.Count.Molecule, sizeof *mol_id_internal);
  int count = 0;
  for (int i = 0; i <= Sys.Count.HighestResid; i++) {
    int id = mol_id[i];
    if (id != -1) {
      count++;
      mol_id_internal[id] = i;
    }
  }
  // check the number of molecules; it shouldn't ever be wrong
  if (Sys.Count.Molecule != count) {
    strcpy(ERROR_MSG, "something went wrong with per bead type bead counting; \
contact developper\n");
    ErrorPrintError();
    exit(1);
  } //}}}
  // 2) assign atom default to default beads & count bonded/unbonded beads //{{{
  for (int i = 0; i < Sys.Count.Bead; i++) {
    if (atom[i].name == -1) { // default bead?
      atom[i] = atom_def;
    }
    if (atom[i].resid == -1) { // unbonded bead?
      Sys.Count.Unbonded++;
    } else {
      Sys.Count.Bonded++; // bonded bead?
    }
  }
  // just check that it counts the beads correctly
  if (Sys.Count.Bead != (Sys.Count.Unbonded + Sys.Count.Bonded)) {
    strcpy(ERROR_MSG, "something went wrong with bead counting; \
contact developper\n");
    ErrorPrintError();
    exit(1);
  } //}}}
  // 3) identify bead types //{{{
  if (detailed) { // check other stuff besides name //{{{
    /*
     * First, identify bead types based on name, charge, masse, and radius,
     * e.g., lines
     *   atom 0 n x q 1 m 1
     *   atom 1 n x q 2 m 1
     * will be of two different types. This can create an excess of bead
     * types, so some may have to be merged.
     *
     * What is to be merged:
     * i) If a keyword is missing in one line but present in another, that
     *    does not count as a different type, e.g., lines
     *       atom 0 n x q 1 m 1
     *       atom 1 n x     m 1
     *    are of the same type (both with charge +1);
     * ii) however, there can be ambiguities, so e.g., lines
     *        atom 0 n x q 1 m 1
     *        atom 1 n x     m 1
     *        atom 2 n x q 0 m 1
     *     remain three distinct types (atom 1 has undefined charge);
     * iii) but only some lines can be ambiguous, e.g., lines
     *        atom 0 n x q 1 m 1
     *        atom 1 n x     m 1
     *        atom 2 n x q 0 m 1
     *        atom 3 n x q 0
     *      are still three different types (the last two should be
     *      considered the same because there is no ambiguity because all
     *      beads have the same mass)
     * iv) note that sometimes the charge/mass/radius can remain undefined
     *     even though there's only one well defined value; e.g., lines
     *       atom 0 n x q 1 m 1
     *       atom 1 n x     m 1
     *       atom 2 n x q 0 m 1
     *       atom 3 n x q 0
     *       atom 4 n x q 0 m 1 r 1
     *     will make radius well defined (with value 1) only for beads
     *     sharing the type with atom 4 (i.e., atoms 2, 3, and 4), while
     *     the first two atoms will still have undefined radius. What
     *     should the radius of atoms 0 and 1 be when the mass/charge are
     *     different to that of the last atom?
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
     */
    // create (possibly too many) bead types according to bead properties //{{{
    if (default_atom != 0) { // create 'atom default' bead type first
      NewBeadType(&Sys.BeadType, &Sys.Count.BeadType, atom_name[atom_def.name],
                  atom_def.charge, atom_def.mass, atom_def.radius);
    }
    for (int i = 0; i < Sys.Count.Bead; i++) {
      int btype = -1;
      for (int j = 0; j < Sys.Count.BeadType; j++) {
        if (strcmp(Sys.BeadType[j].Name, atom_name[atom[i].name]) == 0 &&
            Sys.BeadType[j].Charge == atom[i].charge &&
            Sys.BeadType[j].Mass == atom[i].mass &&
            Sys.BeadType[j].Radius == atom[i].radius) {
          btype = j;
        }
      }
      if (btype == -1) { // new bead type?
        NewBeadType(&Sys.BeadType, &Sys.Count.BeadType, atom_name[atom[i].name],
                    atom[i].charge, atom[i].mass, atom[i].radius);
        btype = Sys.Count.BeadType - 1;
      }
      Sys.BeadType[btype].Number++;
    }
    // count number of beads of default type if 'atom default' is present
    if (default_atom != 0) {
      Sys.BeadType[0].Number = Sys.Count.Bead;
      for (int i = 1; i < Sys.Count.BeadType; i++) {
        Sys.BeadType[0].Number -= Sys.BeadType[i].Number;
      }
    } //}}}
//PrintBeadType2((*Counts).TypesOfBeads, Sys.BeadType);
    // Merging //{{{
    // 1) //{{{
    // arrays values of charge, mass, and radius for each bead type //{{{
    double diff_q[atom_names],
           diff_m[atom_names],
           diff_r[atom_names]; //}}}
    // initialize arrays: assign values from the last type with each name //{{{
    for (int i = 0; i < atom_names; i++) {
      for (int j = 0; j < Sys.Count.BeadType; j++) {
        if (strcmp(atom_name[i], Sys.BeadType[j].Name) == 0) {
          diff_q[i] = Sys.BeadType[j].Charge;
          diff_m[i] = Sys.BeadType[j].Mass;
          diff_r[i] = Sys.BeadType[j].Radius;
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
    int high = 1000000;
    // go through all bead type pairs (including self-pairs)
    for (int i = 0; i < atom_names; i++) {
      for (int j = 0; j < Sys.Count.BeadType; j++) {
        // only consider type pairs with the same name
        if (strcmp(atom_name[i], Sys.BeadType[j].Name) == 0) {
          // charge
          if (diff_q[i] != Sys.BeadType[j].Charge) {
            if (diff_q[i] != CHARGE && Sys.BeadType[j].Charge != CHARGE) {
              diff_q[i] = high;
            } else if (diff_q[i] == CHARGE) {
              diff_q[i] = Sys.BeadType[j].Charge;
            }
          }
          // mass
          if (diff_m[i] != Sys.BeadType[j].Mass) {
            if (diff_m[i] != MASS && Sys.BeadType[j].Mass != MASS) {
              diff_m[i] = high;
            } else if (diff_m[i] == MASS) {
              diff_m[i] = Sys.BeadType[j].Mass;
            }
          }
          // radius
          if (diff_r[i] != Sys.BeadType[j].Radius) {
            if (diff_r[i] != RADIUS && Sys.BeadType[j].Radius != RADIUS) {
              diff_r[i] = high;
            } else if (diff_r[i] == RADIUS) {
              diff_r[i] = Sys.BeadType[j].Radius;
            }
          }
        }
      }
    } //}}}
    //}}}
    /* test print diff_q/m/r //{{{
    for (int i = 0; i < atom_names; i++) {
      printf("%s: q=%10.1f; m=%10.1f; r=%10.1f\n", atom_name[i], diff_q[i], diff_m[i], diff_r[i]);
    }
    */ //}}}
    // 2) //{{{
    // initialize merge array by assuming nothing will be merged //{{{
    bool merge[Sys.Count.BeadType][Sys.Count.BeadType];
    for (int i = 0; i < Sys.Count.BeadType; i++) {
      for (int j = 0; j < Sys.Count.BeadType; j++) {
        merge[i][j] = false; // 'i' and 'j' aren't to be merged
        merge[i][i] = true; // 'i' and 'i' is to be merged/copied
      }
    } //}}}
    // assume same-name bead types are to be merged //{{{
    for (int i = 0; i < (Sys.Count.BeadType-1); i++) {
      for (int j = (i+1); j < Sys.Count.BeadType; j++) {
        if (strcmp(Sys.BeadType[i].Name, Sys.BeadType[j].Name) == 0) {
          merge[i][j] = true;
        }
      }
    } //}}}
    // go through each bead type and compare it to two others //{{{
    for (int i = 0; i < (Sys.Count.BeadType-1); i++) {
      // i)
      int k = 0;
      for (; k < Sys.Count.BeadType; k++) {
        if (strcmp(atom_name[k], Sys.BeadType[i].Name) == 0) {
          break;
        }
      }
      // ii)
      for (int j = (i+1); j < Sys.Count.BeadType; j++) {
        // check charge
        if (merge[i][j]) {
          if (diff_q[k] == high) {
            if (Sys.BeadType[i].Charge == Sys.BeadType[j].Charge) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        }
        // check mass
        if (merge[i][j]) {
          if (diff_m[k] == high) {
            if (Sys.BeadType[i].Mass == Sys.BeadType[j].Mass) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        }
        // check radius
        if (merge[i][j]) {
          if (diff_r[k] == high) {
            if (Sys.BeadType[i].Radius == Sys.BeadType[j].Radius) {
              merge[i][j] = true;
            } else {
              merge[i][j] = false;
            }
          } else {
            merge[i][j] = true;
          }
        }
      }
    } //}}}
    //}}}
    /* test print merge matrix //{{{
    printf("   ");
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      printf("%s ", Sys.BeadType[i].Name);
    }
    putchar('\n');
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      printf("%s ", Sys.BeadType[i].Name);
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        printf(" %d", merge[i][j]);
      }
      putchar('\n');
    }
    */ //}}}
    // 3) //{{{
    BEADTYPE *temp = calloc(Sys.Count.BeadType, sizeof (BEADTYPE));
    int count = 0;
    for (int i = 0; i < Sys.Count.BeadType; i++) {
      if (merge[i][i]) { // i)
        temp[count] = Sys.BeadType[i];
        for (int j = (i+1); j < Sys.Count.BeadType; j++) {
          if (merge[i][j]) { // ii)
            temp[count].Number += Sys.BeadType[j].Number;
            if (temp[count].Charge == CHARGE) {
              temp[count].Charge = Sys.BeadType[j].Charge;
            }
            if (temp[count].Mass == MASS) {
              temp[count].Mass = Sys.BeadType[j].Mass;
            }
            if (temp[count].Radius == RADIUS) {
              temp[count].Radius = Sys.BeadType[j].Radius;
            }
            merge[j][j] = false;
          }
        }
        count++;
      }
    }
    Sys.Count.BeadType = count; //}}}
    //}}}
    // 4) //{{{
    // copy all bead types temporarily to bt struct
    for (int i = 0; i < Sys.Count.BeadType; i++) {
      Sys.BeadType[i] = temp[i];
//    Sys.BeadType[i].Use = false; // will indicate that it wasn't copied yet
    }
    // copy the bead types back to temp array in a proper order
    count = 0;
    for (int i = 0; i < Sys.Count.BeadType; i++) {
      if (!Sys.BeadType[i].Use) {
        temp[count++] = Sys.BeadType[i];
        Sys.BeadType[i].Use = true;
        for (int j = (i+1); j < Sys.Count.BeadType; j++) {
          if (strcmp(Sys.BeadType[i].Name, Sys.BeadType[j].Name) == 0 &&
              !Sys.BeadType[j].Use) {
            temp[count++] = Sys.BeadType[j];
            Sys.BeadType[j].Use = true;
          }
        }
      }
    }
    // finally, copy the types from the temporary array back to bt array
    for (int i = 0; i < Sys.Count.BeadType; i++) {
      Sys.BeadType[i] = temp[i];
    }
    free(temp); //}}}
    //}}}
  } else { // based on name only  //{{{
    // create new bead type from each unique atom_name
    for (int i = 0; i < atom_names; i++) {
      // if this is a 'atom default' name, use its charge, mass & radius
      if (default_atom != 0 && atom_def.name == i) {
        NewBeadType(&Sys.BeadType, &Sys.Count.BeadType, atom_name[i],
                    atom_def.charge, atom_def.mass, atom_def.radius);
      } else { // otherwise, use 'undefined' values
        NewBeadType(&Sys.BeadType, &Sys.Count.BeadType, atom_name[i],
                    CHARGE, MASS, RADIUS);
      }
    }
    // go through all beads to count them & assign stuff to bead types
    for (int i = 0; i < Sys.Count.Bead; i++) {
      int name = atom[i].name;
      int btype = FindBeadType2(atom_name[name],
                                Sys.Count.BeadType, Sys.BeadType);
      // take the first defined values of charge, mass, and radius
      if (btype != -1) {
        Sys.BeadType[btype].Number++;
        if (Sys.BeadType[btype].Mass == MASS && atom[i].mass != MASS) {
          Sys.BeadType[btype].Mass = atom[i].mass;
        }
        if (Sys.BeadType[btype].Charge == CHARGE && atom[i].charge != CHARGE) {
          Sys.BeadType[btype].Charge = atom[i].charge;
        }
        if (Sys.BeadType[btype].Radius == RADIUS && atom[i].radius != RADIUS) {
          Sys.BeadType[btype].Radius = atom[i].radius;
        }
      } else { // weird error that should never happen //{{{
        strcpy(ERROR_MSG, "something went with recognizing bead types; \
    contact developper\n");
        ErrorPrintError();
        exit(1);
      } //}}}
    }
  } //}}}
  // recheck the total number of beads; it shouldn't ever be wrong //{{{
  count = 0;
  for (int i = 0; i < Sys.Count.BeadType; i++) {
    count += Sys.BeadType[i].Number;
  }
  if (Sys.Count.Bead != count) {
    strcpy(ERROR_MSG, "something went wrong with per bead type bead counting; \
contact developper\n");
    ErrorPrintError();
    exit(1);
  } //}}}
  //}}}
  // 4) fill BEAD struct, putting unbonded beads first //{{{
  Sys.Bead = realloc(Sys.Bead, sizeof (BEAD) * Sys.Count.Bead);
  Sys.Bonded = realloc(Sys.Bonded, sizeof *Sys.Bonded * Sys.Count.Bonded);
  Sys.Unbonded = realloc(Sys.Unbonded, sizeof *Sys.Unbonded * Sys.Count.Unbonded);
  int c_bonded = 0, c_unbonded = 0;
  for (int i = 0; i < Sys.Count.Bead; i++) {
    if (atom[i].resid == -1) {
      Sys.Bead[i].Molecule = -1;
      Sys.Unbonded[c_unbonded] = i;
      c_unbonded++;
    } else {
      Sys.Bead[i].Molecule = mol_id[atom[i].resid];
      Sys.Bonded[c_bonded] = i;
      c_bonded++;
    }
    // assign type based on the 'detailed' parameter
    if (detailed) { //{{{
      Sys.Bead[i].Type = -1;
      // 1) if a bead shares all values with a bead type, it is of that type
      for (int j = 0; j < Sys.Count.BeadType; j++) {
        if (strcmp(Sys.BeadType[j].Name, atom_name[atom[i].name]) == 0 &&
            Sys.BeadType[j].Charge == atom[i].charge &&
            Sys.BeadType[j].Mass == atom[i].mass &&
            Sys.BeadType[j].Radius == atom[i].radius) {
          Sys.Bead[i].Type = j;
          break;
        }
      }
      /*
       * 2) if no type was assigned, check for undefined values of the bead's
       *    charge/mass/radius
       */
      if (Sys.Bead[i].Type == -1) {
        for (int j = 0; j < Sys.Count.BeadType; j++) {
          // only check if the bead type and the bead share name
          if (strcmp(Sys.BeadType[j].Name, atom_name[atom[i].name]) == 0) {
            // check charge
            if (Sys.BeadType[j].Charge != atom[i].charge &&
                atom[i].charge != CHARGE) {
              continue;
            }
            // check mass
            if (Sys.BeadType[j].Mass != atom[i].mass &&
                atom[i].mass != MASS) {
              continue;
            }
            // check radius
            if (Sys.BeadType[j].Radius != atom[i].radius &&
                atom[i].radius != RADIUS) {
              continue;
            }
            // assign bead type if all checks passed
            Sys.Bead[i].Type = j;
            break;
          }
        }
      } //}}}
    } else {
      Sys.Bead[i].Type = atom[i].name;
    }
    Sys.Bead[i].nAggregates = 1; // TODO will be removed
  } //}}}
  // 5) if detailed, rename the bead types with the same name //{{{
  if (detailed) {
    for (int i = 0; i < (Sys.Count.BeadType-1); i++) {
      count = 0;
      for (int j = (i+1); j < Sys.Count.BeadType; j++) {
        if (strcmp(Sys.BeadType[i].Name, Sys.BeadType[j].Name) == 0) {
          count++;
          char name[BEAD_NAME+1];
          // shorten name if necessary
          if (count < 10) {
            strncpy(name, Sys.BeadType[j].Name, BEAD_NAME-2);
          } else if (count < 100) {
            strncpy(name, Sys.BeadType[j].Name, BEAD_NAME-3);
          } else if (count < 1000) {
            strncpy(name, Sys.BeadType[j].Name, BEAD_NAME-4);
          }
          // BEAD_NAME is max string length, i.e., array is longer
          snprintf(Sys.BeadType[j].Name, BEAD_NAME+1, "%s_%d", name, count);
        }
      }
    }
  } //}}}
  // 6) identify molecule types based on all data //{{{
  /*
   * Molecules of one type must share:
   * i) molecule name
   * ii) number of beads and bonds
   * iii) order of bead types
   * iv) connectivity
   */
  // count beads in each molecule for ii) //{{{
  int *atoms_per_mol = calloc(Sys.Count.Molecule, sizeof *atoms_per_mol);
  for (int i = 0; i < Sys.Count.Bonded; i++) {
    int id = Sys.Bonded[i];
    atoms_per_mol[Sys.Bead[id].Molecule]++;
  } //}}}
  // allocate MOLECULE struct & and fill with indices for iv) //{{{
  Sys.Molecule = realloc(Sys.Molecule,
                         sizeof (MOLECULE) * Sys.Count.Molecule);
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    Sys.Molecule[i].Bead = calloc(atoms_per_mol[i], sizeof *Sys.Molecule[i].Bead);
    Sys.Molecule[i].Aggregate = -1; // in no aggregate // TODO: will probably disappear
    Sys.Molecule[i].Index = mol_id_internal[i];
  }
  // fill to Molecule[].Bead array with bead indices
  int *count_in_mol = calloc(Sys.Count.Molecule, sizeof *count_in_mol);
  // go through bonded beads
  for (int i = 0; i < Sys.Count.Bonded; i++) {
    int b_id = Sys.Bonded[i],
        m_id = Sys.Bead[b_id].Molecule;
    Sys.Molecule[m_id].Bead[count_in_mol[m_id]] = b_id;
    count_in_mol[m_id]++; // count number of beads placed in each molecule
  }
  // recheck bead numbers in molecules; this should never trigger
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    if (count_in_mol[i] != atoms_per_mol[i]) {
      strcpy(ERROR_MSG, "something went wrong with per molecule bead numbers; \
contact developper\n");
      ErrorPrintError();
    }
  } //}}}
  // count bonds in all molecules for ii) //{{{
  int *bonds_per_mol = calloc(Sys.Count.Molecule, sizeof *bonds_per_mol);
  for (int i = 0; i < count_bonds; i++) {
    int id1 = bond[i].index1, // vsf bead id
        id2 = bond[i].index2; // vsf bead id
    int m_id1 = Sys.Bead[id1].Molecule, // internal molecule id
        m_id2 = Sys.Bead[id2].Molecule; // internal molecule id
    // error - beads from one bond are in different molecules //{{{
    if (m_id1 != m_id2) {
      strcpy(ERROR_MSG, "bonded beads not in the same molecule");
      PrintError();
      ErrorPrintFile(struct_file, "\0");
      fprintf(stderr, "%s, atom id (resid): %s%d %s(%s", ErrRed(), ErrYellow(),
                                                         bond[i].index1,
                                                         ErrRed(), ErrYellow());
      if (m_id1 != -1) {
        fprintf(stderr, "%d", mol_id_internal[m_id1]);
      } else {
        fprintf(stderr, "none");
      }
      fprintf(stderr, "%s); %s%d %s(%s", ErrRed(), ErrYellow(), bond[i].index2,
                                         ErrRed(), ErrYellow());
      if (m_id2 != -1) {
        fprintf(stderr, "%d", mol_id_internal[m_id2]);
      } else {
        fprintf(stderr, "none");
      }
      fprintf(stderr, "%s)%s\n", ErrRed(), ErrColourReset());
      exit(1);
    } //}}}
    bonds_per_mol[m_id1]++;
  }
  //}}}
  // save connectivity for each molecule for iv) //{{{
  // allocate connectivity array //{{{
  /*
   * in molec[i].connect[j][k]:
   *   i ... molecule id
   *   j ... bond id
   *   k ... 0 and 1: connected bead indices; 2: bond type (not used now)
   */
  struct connectivity {
    int (*connect)[3];
  } *molec = malloc(sizeof *molec * Sys.Count.Molecule);
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    count_in_mol[i] = 0; // to count bonds already placed in all molecules
    molec[i].connect = malloc(sizeof *molec[i].connect * bonds_per_mol[i]);
    for (int j = 0; j < bonds_per_mol[i]; j++) {
      molec[i].connect[j][0] = -1;
      molec[i].connect[j][1] = -1;
    }
  } //}}}
  // save ids of bonded beads
  for (int i = 0; i < count_bonds; i++) { // go through all bonds
    int id1 = bond[i].index1, // vsf bead id
        id2 = bond[i].index2; // vsf bead id
    int m_id = Sys.Bead[id1].Molecule; // internal molecule id
    bool done[2] = {false}; // so as not to go through the whole molecule
    // go through beads in m_id to identify molecularly internal bead ids
    for (int j = 0; j < atoms_per_mol[m_id]; j++) {
      if (Sys.Molecule[m_id].Bead[j] == id1) {
        molec[m_id].connect[count_in_mol[m_id]][0] = j;
        done[0] = true;
      }
      if (Sys.Molecule[m_id].Bead[j] == id2) {
        molec[m_id].connect[count_in_mol[m_id]][1] = j;
        done[1] = true;
      }
      if (done[0] && done[1]) {
        break;
      }
    }
    count_in_mol[m_id]++;
  }
  // sort bonds from lowest id to the highest
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    SortBonds(molec[i].connect, bonds_per_mol[i]);
  }
  //}}}
  // count molecule types based on i) through iv) //{{{
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    int id = Sys.Molecule[i].Bead[0]; // vsf id of the first bead of 'i'
    int name = atom[id].resname; // index in res_name array
    bool add = false; // assume no type exists that 'i' can be added to
    // go through all molecule types to find a match for molecule 'i' //{{{
    for (int j = 0; j < Sys.Count.MoleculeType; j++) {
      /* check that molecule 'i' shares with molecule type 'j':
       *   a) resname
       *   b) number of beads
       *   c) number of bonds
       *   d) bead order
       *   e) connectivity
       */
      if (strcmp(res_name[name], Sys.MoleculeType[j].Name) == 0 && // a)
          atoms_per_mol[i] == Sys.MoleculeType[j].nBeads && // b)
          bonds_per_mol[i] == Sys.MoleculeType[j].nBonds) { // c)
        // d)
        bool same_beads = true; // assume molecule 'i' has 'j' type's bead order
        for (int k = 0; k < Sys.MoleculeType[j].nBeads; k++) {
          if (Sys.Bead[Sys.Molecule[i].Bead[k]].Type != Sys.MoleculeType[j].Bead[k]) {
            same_beads = false; // nope, it doesn't; is not type j
            break;
          }
        }
        // e)
        bool same_bonds = true; // assume molecule i has j type's connectivity
        for (int k = 0; k < Sys.MoleculeType[j].nBonds; k++) {
          if (molec[i].connect[k][0] != Sys.MoleculeType[j].Bond[k][0] ||
              molec[i].connect[k][1] != Sys.MoleculeType[j].Bond[k][1]) {
            same_bonds = false; // nope, it doesn't; i is not type j
            break;
          }
        }
        if (same_beads && same_bonds) {
          add = true; // add to existing molecule type
        }
      }
      if (add) { // found matching molecule type
        Sys.MoleculeType[j].Number++;
        Sys.Molecule[i].Type = j;
        break;
      }
    } //}}}
    // create new type if no match found for 'i' or no type exists yet //{{{
    if (!add || Sys.Count.MoleculeType == 0) {
      NewMolType(&Sys.MoleculeType, &Sys.Count.MoleculeType, res_name[name],
                 atoms_per_mol[i], bonds_per_mol[i], 0, 0);
      int type = Sys.Count.MoleculeType - 1;
      for (int j = 0; j < Sys.MoleculeType[type].nBeads; j++) {
        int id = Sys.Molecule[i].Bead[j];
        Sys.MoleculeType[type].Bead[j] = Sys.Bead[id].Type;
      }
      if (bonds_per_mol[i] > 0) {
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          Sys.MoleculeType[type].Bond[j][0] = molec[i].connect[j][0];
          Sys.MoleculeType[type].Bond[j][1] = molec[i].connect[j][1];
          Sys.MoleculeType[type].Bond[j][2] = -1; // no bond types
        }
      }
      // TODO angles & dihedrals

      Sys.Molecule[i].Type = type;
    } //}}}
  } //}}}
  // rename molecule types with the same name //{{{
  for (int i = 0; i < (Sys.Count.MoleculeType-1); i++) {
    count = 0; // number of types sharing the name with type i
    for (int j = (i+1); j < Sys.Count.MoleculeType; j++) {
      if (strcmp(Sys.MoleculeType[i].Name, Sys.MoleculeType[j].Name) == 0) {
        count++;
        // shorten name if necessary to append '_<int>'
        char name[MOL_NAME+1];
        strcpy(name, Sys.MoleculeType[j].Name);
        if (count < 10) {
          name[MOL_NAME-2] = '\0';
        } else if (count < 100) {
          name[MOL_NAME-3] = '\0';
        } else if (count < 1000) {
          name[MOL_NAME-4] = '\0';
        }
        snprintf(Sys.MoleculeType[j].Name, MOL_NAME+1, "%s_%d", name, count);
      }
    }
  } //}}}
  FillMolType(Sys.Count.MoleculeType, Sys.BeadType, &Sys.MoleculeType);
  //}}}
  // 7) copy everything back to their 'proper' arrays and structures //{{{
  Sys.Index_mol = realloc(Sys.Index_mol, sizeof *Sys.Index_mol *
                          (Sys.Count.HighestResid+1));
  for (int i = 0; i <= Sys.Count.HighestResid; i++) {
    Sys.Index_mol[i] = -1;
  }
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    Sys.Index_mol[Sys.Molecule[i].Index] = i;
  } //}}}
  // free memory //{{{
  free(atom);
  free(atom_name);
  free(res_name);
  free(mol_id);
  free(mol_id_internal);
  free(bond);
  free(atoms_per_mol);
  free(bonds_per_mol);
  free(count_in_mol);
  for (int i = 0; i < Sys.Count.Molecule; i++) {
    free(molec[i].connect);
  }
  free(molec);
  //}}}
  // assume empty/none-existent coordinate file
  Sys.Count.BeadCoor = 0;
  Sys.BeadCoor = realloc(Sys.BeadCoor,
                          sizeof *Sys.BeadCoor * Sys.Count.Bead);
  Sys.Count.BondedCoor = 0;
  if (Sys.Count.Bonded > 0) {
    Sys.BondedCoor = realloc(Sys.BondedCoor,
                             sizeof *Sys.BondedCoor * Sys.Count.Bonded);
  }
  Sys.Count.UnbondedCoor = 0;
  if (Sys.Count.Unbonded > 0) {
    Sys.UnbondedCoor = realloc(Sys.UnbondedCoor,
                               sizeof *Sys.UnbondedCoor * Sys.Count.Unbonded);
  }
  WarnChargedSystem(Sys, struct_file, "\0");
  return Sys;
} //}}}
// TODO remove, I think
// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined.
 */
void FillMolMass(int number_of_types,
                 BEADTYPE *BeadType, MOLECULETYPE *MoleculeType[]) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      // undefined mass for a bead type => undefined molecule mass
      if (BeadType[type].Mass == MASS) {
        (*MoleculeType)[i].Mass = MASS;
        break;
      } else {
        (*MoleculeType)[i].Mass += BeadType[type].Mass;
      }
    }
    if ((*MoleculeType)[i].Mass == MASS) {
      strcpy(ERROR_MSG, "molecule type with undefined mass");
      PrintWarning();
      fprintf(stderr, "%sMolecule %s%s%s\n", ErrCyan(), ErrYellow(),
                                             (*MoleculeType)[i].Name,
                                             ErrColourReset());
    }
  }
} //}}}
// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined.
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType,
                   MOLECULETYPE *MoleculeType[]) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].Charge = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int type = (*MoleculeType)[i].Bead[j];
      // undefined charge for a bead type => undefined molecule charge
      if (BeadType[type].Charge == CHARGE) {
        (*MoleculeType)[i].Charge = CHARGE;
        break;
      } else {
        (*MoleculeType)[i].Charge += BeadType[type].Charge;
      }
    }
    if ((*MoleculeType)[i].Charge == CHARGE) {
      strcpy(ERROR_MSG, "molecule type with undefined charge");
      PrintWarning();
      fprintf(stderr, "%sMolecule %s%s%s\n", ErrCyan(), ErrYellow(),
                                             (*MoleculeType)[i].Name,
                                             ErrColourReset());
    }
  }
} //}}}
// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type.
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType,
                 MOLECULETYPE *MoleculeType[]) {
  FillMolBTypes_old(number_of_types, MoleculeType);
  FillMolMass(number_of_types, BeadType, MoleculeType);
  FillMolCharge(number_of_types, BeadType, MoleculeType);
} //}}}
#endif
