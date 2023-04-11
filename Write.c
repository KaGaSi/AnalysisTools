#include "Write.h"

// STATIC DEFINITIONS
static void VtfWriteCoorIndexed(FILE *fw, bool write[], SYSTEM System);
static void XyzWriteCoor(FILE *fw, bool write[], SYSTEM System);
static void LtrjWriteCoor(FILE *fw, int step, bool write[], SYSTEM System);
static void WriteConfig(SYSTEM System, char file[]);
static void VtfWriteStruct(char file[], SYSTEM System, int type_def);
static void WriteLmpData(SYSTEM System, char file[], bool mass);
static void WriteField(SYSTEM System, char file_field[]);

// STATIC IMPLEMENTATIONS
static void VtfWriteCoorIndexed(FILE *fw, bool write[], SYSTEM System) { //{{{
  // print box size if present //{{{
  BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "pbc %lf %lf %lf", box->Length.x, box->Length.y, box->Length.z);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, "    %lf %lf %lf", box->alpha, box->beta, box->gamma);
    }
    putc('\n', fw);
  } //}}}
  fprintf(fw, "indexed\n");
  bool none = true;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    if (bead->InTimestep && write[id]) {
      none = false;
      fprintf(fw, "%8d %8.4f %8.4f %8.4f\n", id, bead->Position.x,
                                                 bead->Position.y,
                                                 bead->Position.z);
    }
  }
  if (none) {
    strcpy(ERROR_MSG, "no beads to save");
    PrintWarning();
  }
} //}}}
static void XyzWriteCoor(FILE *fw, bool write[], SYSTEM System) { //{{{
  // find out number of beads to save
  int count = 0;
  bool none = true; // to make sure there are beads to save
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (System.Bead[id].InTimestep && write[id]) {
      none = false;
      count++;
    }
  }
  if (none) {
    strcpy(ERROR_MSG, "no beads to save");
    PrintWarning();
  } else {
    // TODO: write pbc on the second line?
    fprintf(fw, "%d\n", count);
    BOX *box = &System.Box;
    if (box->Volume != -1) {
      fprintf(fw, "%lf %lf %lf", box->Length.x, box->Length.y, box->Length.z);
      if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
        fprintf(fw, " %lf %lf %lf", box->alpha, box->beta, box->gamma);
      }
    }
    putc('\n', fw);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *bead = &System.Bead[id];
      if (bead->InTimestep && write[id]) {
        int type = bead->Type;
        fprintf(fw, "%8s %8.4f %8.4f %8.4f\n", System.BeadType[type].Name,
                bead->Position.x, bead->Position.y, bead->Position.z);
      }
    }
  }
} //}}}
static void LtrjWriteCoor(FILE *fw, int step, bool write[], SYSTEM System) { //{{{
  // find out number of beads to save and if velocity/force should be saved
  int count = 0;
  bool vel = false, force = false;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *b = &System.Bead[id];
    if (b->InTimestep && write[id]) {
      count++;
      if (b->Velocity.x != 0 || b->Velocity.y != 0 || b->Velocity.z != 0) {
        vel = true;
      }
      if (b->Force.x != 0 || b->Force.y != 0 || b->Force.z != 0) {
        force = true;
      }
    }
  }
  // print the step
  if (count > 0) {
    BOX *box = &System.Box;
    fprintf(fw, "ITEM: TIMESTEP\n%d\n", step);
    fprintf(fw, "ITEM: NUMBER OF ATOMS\n%d\n", count);
    if (box->Volume == -1) {
      strcpy(ERROR_MSG, "unspecified box dimensions");
      PrintWarning();
    }
    // orthogonal box
    if (box->alpha == 90 && box->beta == 90 && box->gamma == 90) {
      fprintf(fw, "ITEM: BOX BOUNDS pp pp pp\n");
      fprintf(fw, "0.0 %lf\n", box->Length.x);
      fprintf(fw, "0.0 %lf\n", box->Length.y);
      fprintf(fw, "0.0 %lf\n", box->Length.z);
    } else {
      fprintf(fw, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding.x, box->transform[0][1]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding.y, box->transform[0][2]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding.z, box->transform[1][2]);
    }
    fprintf(fw, "ITEM: ATOMS id element x y z");
    if (vel) {
      fprintf(fw, " vx vy vz");
    }
    if (force) {
      fprintf(fw, " fx fy fz");
    }
    putc('\n', fw);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *b = &System.Bead[id];
      if (b->InTimestep && write[id]) {
        int type = b->Type;
        fprintf(fw, "%8d %8s %8.4f %8.4f %8.4f", id + 1,
                System.BeadType[type].Name,
                b->Position.x, b->Position.y, b->Position.z);
        if (vel) {
          fprintf(fw, " %8.4f %8.4f %8.4f",
                  b->Velocity.x, b->Velocity.y, b->Velocity.z);
        }
        if (force) {
          fprintf(fw, " %8.4f %8.4f %8.4f", b->Force.x, b->Force.y, b->Force.z);
        }
        putc('\n', fw);
      }
    }
  } else {
    strcpy(ERROR_MSG, "no beads to save");
    WarnPrintWarning();
  }
} //}}}
static void WriteConfig(SYSTEM System, char file[]) { //{{{
  FILE *out = OpenFile(file, "w");
  // TODO: check triclinic box in dl_meso
  // print CONFIG file initial stuff
  fprintf(out, "NAME\n 0 1\n"); // not sure what 0 1 is...
  fprintf(out, "%lf 0.000000 0.000000\n", System.Box.Length.x);
  fprintf(out, "0.000000 %lf 0.000000\n", System.Box.Length.y);
  fprintf(out, "0.000000 0.000000 %lf\n", System.Box.Length.z);

  // bead coordinates
  // unbonded beads must be first (dl_meso requirement)
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i], btype = System.Bead[id].Type;
    fprintf(out, "%s %d\n", System.BeadType[btype].Name, id + 1);
    fprintf(out, "%lf %lf %lf\n", System.Bead[id].Position.x,
            System.Bead[id].Position.y, System.Bead[id].Position.z);
  }
  fclose(out);
} //}}}
static void VtfWriteStruct(char file[], SYSTEM System, int type_def) { //{{{
  FILE *fw = OpenFile(file, "w");
  COUNT *Count = &System.Count;
  // default bead type //{{{
  if (type_def == -1) {
    // find most common type of bead and make it default
    int *count = calloc(Count->BeadType, sizeof *count);
    for (int i = 0; i < Count->Bead; i++) {
      if (System.Bead[i].Molecule == -1) {
        int type = System.Bead[i].Type;
        count[type]++;
      }
    }
    int max = 0;
    for (int i = 0; i < Count->BeadType; i++) {
      if (count[i] > max) {
        max = count[i];
        type_def = i;
      }
    }
    free(count);
  } //}}}
  // print default bead type //{{{
  if (type_def != -1) {
    BEADTYPE *bt = &System.BeadType[type_def];
    fprintf(fw, "atom default name %8s", bt->Name);
    if (bt->Mass != MASS) {
      fprintf(fw, " mass %lf", bt->Mass);
    }
    if (bt->Charge != CHARGE) {
      fprintf(fw, " charge %lf", bt->Charge);
    }
    if (bt->Radius != RADIUS) {
      fprintf(fw, " radius %lf", bt->Radius);
    }
    putc('\n', fw);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    int btype = System.Bead[i].Type, mol = System.Bead[i].Molecule;
    // print beads that are non-default, in a molecule, or have the highest id
    bool print = false;
    if (btype != type_def || mol != -1 || i == (Count->Bead - 1)) {
      print = true;
    }
    BEADTYPE *bt = &System.BeadType[btype];
    if (print) {
      fprintf(fw, "atom %7d name %8s", i, bt->Name);
      if (bt->Mass != MASS) {
        fprintf(fw, " mass %lf ", bt->Mass);
      }
      if (bt->Charge != CHARGE) {
        fprintf(fw, " charge %lf", bt->Charge);
      }
      if (bt->Radius != RADIUS) {
        fprintf(fw, " radius %lf", bt->Radius);
      }
      if (mol != -1) {
        int mtype = System.Molecule[mol].Type,
            id = System.Molecule[mol].Index;
        fprintf(fw, " resname %10s ", System.MoleculeType[mtype].Name);
        fprintf(fw, "resid %5d", id);
      }
      putc('\n', fw);
    }
  } //}}}
  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Count->Molecule; i++) {
    fprintf(fw, "# resid %d\n", i + 1); // in VMD resid start with 1
    int mol_type = System.Molecule[i].Type;
    for (int j = 0; j < System.MoleculeType[mol_type].nBonds; j++) {
      int *bead = BondIndices(System, i, j);
      fprintf(fw, "bond %6d: %6d\n", bead[0], bead[1]);
    }
  } //}}}
  // close structure file
  fclose(fw);
} //}}}
static void WriteLmpData(SYSTEM System, char file[], bool mass) { //{{{
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "Created via AnalysisTools v%s "
          "(https://github.com/KaGaSi/AnalysisTools)\n\n", VERSION);
  COUNT *Count = &System.Count;
  // PrintCount(*Count);
  // count bonds according to Bead[].InTimestep
  // create new SYSTEM structure to figure out bead types if mass == true //{{{
  int mass_types = 0;
  int *bt_masstype_to_old = calloc(Count->BeadType, sizeof *bt_masstype_to_old);
  int *bt_old_to_masstype = calloc(Count->BeadType, sizeof *bt_old_to_masstype);
  if (mass) {
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt_i = &System.BeadType[i];
      bool found = false;
      for (int j = 0; j < mass_types; j++) {
        BEADTYPE *bt_j = &System.BeadType[bt_masstype_to_old[j]];
        if (bt_i->Mass == bt_j->Mass) {
          // bt_masstype_to_old[j] = i;
          bt_old_to_masstype[i] = j;
          found = true;
          break;
        }
      }
      if (!found) {
        bt_masstype_to_old[mass_types] = i;
        bt_old_to_masstype[i] = mass_types;
        mass_types++;
      }
    }
  } else {
    mass_types = Count->BeadType;
    for (int i = 0; i < Count->BeadType; i++) {
      bt_old_to_masstype[i] = i;
      bt_masstype_to_old[i] = i;
    }
  } //}}}
  // print counts //{{{
  fprintf(fw, "%10d atoms\n", Count->Bead);
  fprintf(fw, "%10d bonds\n", Count->Bond);
  fprintf(fw, "%10d angles\n", Count->Angle);
  fprintf(fw, "%10d dihedrals\n", Count->Dihedral);
  fprintf(fw, "%10d impropers\n", Count->Improper);
  putc('\n', fw);
  fprintf(fw, "%10d atom types\n", Count->BeadType);
  if (Count->Bond > 0 && Count->BondType == 0) {
    fprintf(fw, "       ???");
  } else {
    fprintf(fw, "%10d", Count->BondType);
  }
  fprintf(fw, " bond types\n");
  if (Count->Angle > 0 && Count->AngleType == 0) {
    fprintf(fw, "       ???");
  } else {
    fprintf(fw, "%10d", Count->AngleType);
  }
  fprintf(fw, " angle types\n");
  if (Count->Dihedral > 0 && Count->DihedralType == 0) {
    fprintf(fw, "       ???");
  } else {
    fprintf(fw, "%10d", Count->DihedralType);
  }
  fprintf(fw, " dihedral types\n");
  if (Count->Improper > 0 && Count->ImproperType == 0) {
    fprintf(fw, "       ???");
  } else {
    fprintf(fw, "%10d", Count->ImproperType);
  }
  fprintf(fw, " improper types\n");
  putc('\n', fw); //}}}
  // print box size //{{{
  fprintf(fw, "0.0 %lf xlo xhi\n", System.Box.OrthoLength.x);
  fprintf(fw, "0.0 %lf ylo yhi\n", System.Box.OrthoLength.y);
  fprintf(fw, "0.0 %lf zlo zhi\n", System.Box.OrthoLength.z);
  if (System.Box.alpha != 90 ||
      System.Box.beta != 90 ||
      System.Box.gamma != 90) {
    fprintf(fw, "%lf %lf %lf xy xz yz\n", System.Box.transform[0][1],
                                          System.Box.transform[0][2],
                                          System.Box.transform[1][2]);
  }
  putc('\n', fw); //}}}
  // print bead type masses //{{{
  fprintf(fw, "Masses\n\n");
  for (int i = 0; i < mass_types; i++) {
    BEADTYPE *bt = &System.BeadType[bt_masstype_to_old[i]];
    if (bt->Mass == MASS) {
      fprintf(fw, "%5d ??? # %s\n", i + 1, bt->Name);
    } else {
      fprintf(fw, "%5d %lf # %s\n", i + 1, bt->Mass, bt->Name);
    }
  } //}}}
  // print various coeffs //{{{
  fprintf(fw, "\nPair Coeffs\n\n");
  for (int i = 0; i < Count->BeadType; i++) {
    fprintf(fw, "%5d ???\n", i + 1);
  }
  if (Count->Bond > 0) {
    fprintf(fw, "\nBond Coeffs\n\n");
    if (Count->BondType > 0) {
      for (int i = 0; i < Count->BondType; i++) {
        fprintf(fw, "%5d %lf %lf\n", i + 1, System.BondType[i].a / 2,
                System.BondType[i].b);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count->Angle > 0) {
    fprintf(fw, "\nAngle Coeffs\n\n");
    if (Count->AngleType > 0) {
      for (int i = 0; i < Count->AngleType; i++) {
        fprintf(fw, "%5d %lf %lf\n", i + 1, System.AngleType[i].a / 2,
                System.AngleType[i].b);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count->Dihedral > 0) {
    fprintf(fw, "\nDihedral Coeffs\n\n");
    if (Count->DihedralType > 0) {
      for (int i = 0; i < Count->DihedralType; i++) {
        fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.DihedralType[i].a / 2,
                System.DihedralType[i].b, System.DihedralType[i].c);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count->Improper > 0) {
    fprintf(fw, "\nImproper Coeffs\n\n");
    if (Count->ImproperType > 0) {
      for (int i = 0; i < Count->ImproperType; i++) {
        fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.ImproperType[i].a / 2,
                System.ImproperType[i].b, System.ImproperType[i].c);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  } //}}}
  // print atoms //{{{
  fprintf(fw, "\nAtoms\n\n");
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    // <bead id>
    fprintf(fw, "%7d", id+1);
    // <molecule id (-1 for no molecule)>
    int mol = bead->Molecule;
    if (mol != -1) {
      fprintf(fw, "%5d", System.Molecule[mol].Index);
    } else {
      fprintf(fw, "%5d", -1);
    }
    // <bead type id>
    fprintf(fw, "%5d", bt_old_to_masstype[bead->Type]+1);
    // <charge> from original System as the charge can differ in lmp data file
    int type = bead->Type;
    double q = System.BeadType[type].Charge;
    if (q == CHARGE) {
      fprintf(fw, "%15s", "???");
    } else {
      fprintf(fw, "%15f", q);
    }
    // coordinates
    fprintf(fw, "%15f %15f %15f", bead->Position.x,
                                  bead->Position.y,
                                  bead->Position.z);
    putc('\n', fw);
  } //}}}
  // print velocities (if at least one non-zero) //{{{
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    VECTOR *vel = &System.Bead[id].Velocity;
    if (fabs(vel->x) > 1e-5 || fabs(vel->y) > 1e-5 || fabs(vel->z) > 1e-5) {
      fprintf(fw, "\nVelocities\n\n");
      for (int j = 0; j < Count->BeadCoor; j++) {
        id = System.BeadCoor[j];
        vel = &System.Bead[id].Velocity;
        fprintf(fw, "%7d", id + 1);
        fprintf(fw, "%15f %15f %15f", vel->x, vel->y, vel->z);
        putc('\n', fw);
      }
      break;
    }
  } //}}}
  // print bonds //{{{
  if (Count->Bond > 0) {
    fprintf(fw, "\nBonds\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = System.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nBonds; j++) {
        count++;
        int *val = BondIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep && System.Bead[val[1]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->BondType > 0) {
            fprintf(fw, "%6d", mt_i->Bond[j][2] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d\n", val[0] + 1, val[1] + 1);
        }
      }
    }
  } //}}}
  // print angles //{{{
  if (Count->Angle > 0) {
    fprintf(fw, "\nAngles\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = System.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nAngles; j++) {
        count++;
        int *val = AngleIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep && System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->AngleType > 0) {
            fprintf(fw, "%6d", mt_i->Angle[j][3] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d\n", val[0] + 1, val[1] + 1, val[2] + 1);
        }
      }
    }
  }
  //}}}
  // print dihedrals //{{{
  if (Count->Dihedral > 0) {
    fprintf(fw, "\nDihedrals\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = System.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        count++;
        int *val = DihedralIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep && System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep && System.Bead[val[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->DihedralType > 0) {
            fprintf(fw, "%6d", mt_i->Dihedral[j][4] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n", val[0] + 1, val[1] + 1, val[2] + 1,
                  val[3] + 1);
        }
      }
    }
  } //}}}
  // print impropers //{{{
  if (Count->Improper > 0) {
    fprintf(fw, "\nImpropers\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = System.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nImpropers; j++) {
        count++;
        int *val = ImproperIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep &&
            System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep &&
            System.Bead[val[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->DihedralType > 0) {
            fprintf(fw, "%6d", mt_i->Dihedral[j][4] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n", val[0] + 1,
                                            val[1] + 1,
                                            val[2] + 1,
                                            val[3] + 1);
        }
      }
    }
  } //}}}
  free(bt_masstype_to_old);
  free(bt_old_to_masstype);
  fclose(fw);
} //}}}
static void WriteField(SYSTEM System, char file_field[]) { //{{{
  FILE *fw = OpenFile(file_field, "w");
  fprintf(fw, "Created via AnalysisTools v%s"
          "(https://github.com/KaGaSi/AnalysisTools)\n\n", VERSION);
  COUNT *Count = &System.Count;
  // print species section //{{{
  fprintf(fw, "species %d <name> <m> <q> <# of unbonded>\n", Count->BeadType);
  // count unbonded beads of each type
  int *unbonded = calloc(Count->BeadType, sizeof *unbonded);
  for (int i = 0; i < Count->Unbonded; i++) {
    int id = System.Unbonded[i];
    int btype = System.Bead[id].Type;
    unbonded[btype]++;
  }
  // print the lines
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System.BeadType[i];
    fprintf(fw, "%16s", bt_i->Name);
    if (bt_i->Mass == MASS) {
      fprintf(fw, " %8s", "???");
    } else {
      fprintf(fw, " %8.5f", bt_i->Mass);
    }
    if (bt_i->Charge == CHARGE) {
      fprintf(fw, " %8s", "???");
    } else {
      fprintf(fw, " %8.5f", bt_i->Charge);
    }
    fprintf(fw, " %5d\n", unbonded[i]);
  }
  free(unbonded); //}}}
  // print molecules section //{{{
  fprintf(fw, "molecule %d\n", Count->MoleculeType);
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System.MoleculeType[i];
    fprintf(fw, "%s\n", mt_i->Name);
    fprintf(fw, "nummols %d\n", mt_i->Number);
    fprintf(fw, "beads %d\n", mt_i->nBeads);
    int mol = mt_i->Index[0];
    // beads
    for (int j = 0; j < mt_i->nBeads; j++) {
      int id = System.Molecule[mol].Bead[j];
      int bt = mt_i->Bead[j];
      VECTOR *pos = &System.Bead[id].Position;
      fprintf(fw, "%16s %8.5f %8.5f %8.5f\n", System.BeadType[bt].Name, pos->x,
              pos->y, pos->z);
    }
    // bonds (if present)
    if (mt_i->nBonds > 0) {
      fprintf(fw, "bonds %d\n", mt_i->nBonds);
      for (int j = 0; j < mt_i->nBonds; j++) {
        // TODO harm only for now
        fprintf(fw, "harm %5d %5d", mt_i->Bond[j][0] + 1, mt_i->Bond[j][1] + 1);
        int type = mt_i->Bond[j][2];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n", System.BondType[type].a,
                  System.BondType[type].b);
        } else {
          fprintf(fw, " ???   ???\n");
        }
      }
    }
    // angles (if present)
    if (mt_i->nAngles > 0) {
      fprintf(fw, "angles %d\n", mt_i->nAngles);
      for (int j = 0; j < mt_i->nAngles; j++) {
        // TODO harm only for now
        fprintf(fw, "harm %5d %5d %5d", mt_i->Angle[j][0] + 1,
                mt_i->Angle[j][1] + 1, mt_i->Angle[j][2] + 1);
        int type = mt_i->Angle[j][3];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n", System.AngleType[type].a,
                  System.AngleType[type].b);
        } else {
          fprintf(fw, " ???   ???\n");
        }
      }
    }
    // dihedrals (if present)
    if (mt_i->nDihedrals > 0) {
      fprintf(fw, "dihedrals %d ", mt_i->nDihedrals);
      fprintf(fw, "# lammps' harmonic style\n");
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        // TODO harm (lammps) only for now
        fprintf(fw, "harm %5d %5d %5d %5d", mt_i->Dihedral[j][0] + 1,
                mt_i->Dihedral[j][1] + 1, mt_i->Dihedral[j][2] + 1,
                mt_i->Dihedral[j][3] + 1);
        int type = mt_i->Dihedral[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.DihedralType[type].a,
                  System.DihedralType[type].b, System.DihedralType[type].c);
        } else {
          fprintf(fw, " ???\n");
        }
      }
    }
    // impropers (if present)
    if (mt_i->nImpropers > 0) {
      fprintf(fw, "impropers %d ", mt_i->nImpropers);
      fprintf(fw, "# lammps' cvff style\n");
      for (int j = 0; j < mt_i->nImpropers; j++) {
        // TODO harm only for now
        fprintf(fw, "cvff %5d %5d %5d %5d", mt_i->Improper[j][0],
                mt_i->Improper[j][1], mt_i->Improper[j][2],
                mt_i->Improper[j][3]);
        int type = mt_i->Improper[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.ImproperType[type].a,
                  System.ImproperType[type].b, System.ImproperType[type].c);
        } else {
          fprintf(fw, " ???\n");
        }
      }
    }
    fprintf(fw, "finish\n");
  } //}}}
  fclose(fw);
} //}}}

// Write a single timestep to output file based on the file type //{{{
void WriteTimestep(int coor_type, char file[], SYSTEM System,
                   int count_step, bool write[]) {
  FILE *fw = OpenFile(file, "a");
  switch (coor_type) {
    case VCF_FILE:
      VtfWriteCoorIndexed(fw, write, System);
      break;
    case XYZ_FILE:
      XyzWriteCoor(fw, write, System);
      break;
    case LTRJ_FILE:
      LtrjWriteCoor(fw, count_step, write, System);
      break;
    default:
      strcpy(ERROR_MSG, "Inexistent output coor_type; should never happen!");
      PrintError();
      exit(1);
  }
  fclose(fw);
} //}}}
// Create a structure file based on the file type (including dl_meso CONFIG) //{{{
void WriteStructure(int struct_type, char file[], SYSTEM System,
                    int vsf_def_type, bool lmp_mass) {
  FILE *fw = OpenFile(file, "w");
  switch (struct_type) {
    case VSF_FILE:
      VtfWriteStruct(file, System, vsf_def_type);
      break;
    case LDATA_FILE:
      WriteLmpData(System, file, lmp_mass);
      break;
    case CONFIG_FILE:
      WriteConfig(System, file);
      break;
    case FIELD_FILE:
      WriteField(System, file);
      break;
    default:
      strcpy(ERROR_MSG, "Inexistent output struct_type; should never happen!");
      PrintError();
      exit(1);
  }
  fclose(fw);
} //}}}

#if 0  //{{{
// TODO will change
// WriteAggregates() //{{{
/*
 * Append aggregate information from a timestep to an .agg file
 */
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                     MOLECULETYPE *MoleculeType, BEAD *Bead,
                     AGGREGATE *Aggregate) {
  // get number of aggregates to write to agg_file //{{{
  int number_of_aggs = 0;
  for (int i = 0; i < Counts.Aggregates; i++) {
    if (Aggregate[i].Use) {
      number_of_aggs++;
    }
  } //}}}

  FILE *fw = OpenFile(agg_file, "a");

  // print number of aggregates to agg file //{{{
  fprintf(fw, "\nStep: %d\n%d\n\n", step_count, number_of_aggs);
  // go through all aggregates
  for (int i = 0; i < Counts.Aggregates; i++) {
    // write only those that aren't excluded
    if (Aggregate[i].Use) {
      // go through all molecules in aggregate 'i'
      fprintf(fw, "%d :", Aggregate[i].nMolecules);
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        fprintf(fw, " %d", Aggregate[i].Molecule[j] + 1);
      }
      putc('\n', fw);

      // go through all monomeric beads in aggregate 'i'
      fprintf(fw, "   %d :", Aggregate[i].nMonomers);
      for (int j = 0; j < Aggregate[i].nMonomers; j++) {
        fprintf(fw, " %d", Aggregate[i].Monomer[j]);
      }
      putc('\n', fw);
    }
  } //}}}

  fclose(fw);
} //}}}
#endif //}}}
