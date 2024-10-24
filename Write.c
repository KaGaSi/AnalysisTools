#include "Write.h"
#include "AnalysisTools.h"
#include "General.h"
#include <stdio.h>
#include <unistd.h>

// STATIC DEFINITIONS
static void VtfWriteCoorIndexed(FILE *fw, const bool write[], SYSTEM System);
static void XyzWriteCoor(FILE *fw, const bool write[], SYSTEM System);
static void LtrjWriteCoor(FILE *fw, int step, const bool write[],
                          SYSTEM System);
static void WriteConfig(SYSTEM System, char file[]);
static void VtfWriteStruct(char file[], SYSTEM System, int type_def,
                           int argc, char *argv[]);
static void WriteLmpData(SYSTEM System, char file[], bool mass,
                         int argc, char *argv[]);
static void WriteField(SYSTEM System, char file_field[],
                       int argc, char *argv[]);
static void SimplifyResid(SYSTEM *System);

// TODO: renumber molecules so that the lowest id is 1 while creating new output
//       structure files

// STATIC IMPLEMENTATIONS
// VtfWriteCoorIndexed() //{{{
static void VtfWriteCoorIndexed(FILE *fw, const bool write[], SYSTEM System) {
  fprintf(fw, "indexed\n");
  // print box size if present //{{{
  BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "pbc %lf %lf %lf", box->Length[0],
                                   box->Length[1],
                                   box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, "    %.3f %.3f %.3f", box->alpha, box->beta, box->gamma);
    }
    putc('\n', fw);
  } //}}}
  bool none = true;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    if (write[id]) {
      none = false;
      fprintf(fw, "%8d %8.4f %8.4f %8.4f\n", id, bead->Position[0],
                                                 bead->Position[1],
                                                 bead->Position[2]);
    }
  }
  if (none) {
    err_msg("no beads to save");
    PrintWarning();
  }
} //}}}
static void XyzWriteCoor(FILE *fw, const bool write[], SYSTEM System) { //{{{
  // find out number of beads to save
  int count = 0;
  bool none = true; // to make sure there are beads to save
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (write[id]) {
      none = false;
      count++;
    }
  }
  if (none) {
    err_msg("no beads to save");
    PrintWarning();
    return;
  }
  // write pbc on the second line
  fprintf(fw, "%d\n", count);
  BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "%.3f %.3f %.3f", box->Length[0],
                                  box->Length[1],
                                  box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, " %lf %lf %lf", box->alpha, box->beta, box->gamma);
    }
  }
  putc('\n', fw);
  // write the coodinates
  // for (int i = 0; i < System.Count.BeadCoor; i++) {
  //   int id = System.BeadCoor[i];
  for (int id = 0; id < System.Count.Bead; id++) {
    BEAD *bead = &System.Bead[id];
    if (write[id] && bead->InTimestep) {
      int type = bead->Type;
      fprintf(fw, "%8s %8.4f %8.4f %8.4f\n", System.BeadType[type].Name,
              bead->Position[0], bead->Position[1], bead->Position[2]);
    }
  }
} //}}}
// TODO: possibly an option to have bead ids to go from 1 to number of beads
// LtrjWriteCoor() //{{{
static void LtrjWriteCoor(FILE *fw, int step, const bool write[],
                          SYSTEM System) {
  // find out number of beads to save and if velocity/force should be saved
  int count_write = 0;
  bool vel = false;
  bool force = false;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *b = &System.Bead[id];
    if (write[id]) {
      count_write++;
      if (b->Velocity[0] != 0 ||
          b->Velocity[1] != 0 ||
          b->Velocity[2] != 0) {
        vel = true;
      }
      if (b->Force[0] != 0 ||
          b->Force[1] != 0 ||
          b->Force[2] != 0) {
        force = true;
      }
    }
  }
  // print the step
  if (count_write > 0) {
    BOX *box = &System.Box;
    fprintf(fw, "ITEM: TIMESTEP\n%d\n", step);
    fprintf(fw, "ITEM: NUMBER OF ATOMS\n%d\n", count_write);
    if (box->Volume == -1) {
      err_msg("unspecified box dimensions");
      PrintWarning();
    }
    // orthogonal box
    if (box->alpha == 90 && box->beta == 90 && box->gamma == 90) {
      fprintf(fw, "ITEM: BOX BOUNDS pp pp pp\n");
      fprintf(fw, "%lf %lf\n", box->Low[0], box->Length[0]+box->Low[0]);
      fprintf(fw, "%lf %lf\n", box->Low[1], box->Length[1]+box->Low[1]);
      fprintf(fw, "%lf %lf\n", box->Low[2], box->Length[2]+box->Low[2]);
    } else {
      fprintf(fw, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[0], box->transform[0][1]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[1], box->transform[0][2]);
      fprintf(fw, "0.0 %lf %lf\n", box->Bounding[2], box->transform[1][2]);
    }
    fprintf(fw, "ITEM: ATOMS id element x y z");
    if (vel) {
      fprintf(fw, " vx vy vz");
    }
    if (force) {
      fprintf(fw, " fx fy fz");
    }
    // fprintf(fw, " mol");
    putc('\n', fw);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *b = &System.Bead[id];
      if (write[id]) {
        int type = b->Type;
        fprintf(fw, "%8d %8s %8.4f %8.4f %8.4f", id + 1,
                System.BeadType[type].Name,
                b->Position[0]+box->Low[0],
                b->Position[1]+box->Low[1],
                b->Position[2]+box->Low[2]);
        if (vel) {
          for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %8.4f", b->Velocity[dd]);
          }
        }
        if (force) {
          for (int dd = 0; dd < 3; dd++) {
            fprintf(fw, " %8.4f", b->Force[dd]);
          }
        }
        // fprintf(fw, " %5d", b->Molecule);
        putc('\n', fw);
      }
    }
  } else {
    err_msg("no beads to save");
    PrintWarning();
  }
} //}}}
static void WriteConfig(SYSTEM System, char file[]) { //{{{
  FILE *out = OpenFile(file, "w");
  // TODO: check triclinic box in dl_meso
  // print CONFIG file initial stuff
  fprintf(out, "NAME\n 0 1\n"); // not sure what 0 1 is...
  fprintf(out, "%.3f 0.000 0.000\n", System.Box.Length[0]);
  fprintf(out, "0.000 %.3f 0.000\n", System.Box.Length[1]);
  fprintf(out, "0.000 0.000 %.3f\n", System.Box.Length[2]);

  // bead coordinates
  // unbonded beads must be first (dl_meso requirement)
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    int btype = System.Bead[id].Type;
    fprintf(out, "%s %d\n", System.BeadType[btype].Name, id + 1);
    fprintf(out, "%lf %lf %lf\n", System.Bead[id].Position[0],
                                  System.Bead[id].Position[1],
                                  System.Bead[id].Position[2]);
  }
  fclose(out);
} //}}}
// VtfWriteStruct() //{{{
static void VtfWriteStruct(char file[], SYSTEM System, int type_def,
                           int argc, char *argv[]) {
  SimplifyResid(&System);
  PrintByline(file, argc, argv);
  FILE *fw = OpenFile(file, "a");
  COUNT *Count = &System.Count;
  BOX *box = &System.Box;
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
    fprintf(fw, "atom default");
    fprintf(fw, " name %8s", bt->Name);
    if (bt->Mass != MASS && bt->Mass != HIGHNUM) {
      fprintf(fw, " mass %12f", bt->Mass);
    }
    if (bt->Charge != CHARGE && bt->Charge != HIGHNUM) {
      fprintf(fw, " charge %12f", bt->Charge);
    }
    if (bt->Radius != RADIUS && bt->Radius != HIGHNUM) {
      fprintf(fw, " radius %12f", bt->Radius);
    }
    putc('\n', fw);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    int btype = System.Bead[i].Type;
    int mol = System.Bead[i].Molecule;
    // print beads that are non-default, in a molecule, or have the highest id
    bool print = false;
    if (btype != type_def || mol != -1 || i == (Count->Bead - 1)) {
      print = true;
    }
    BEADTYPE *bt = &System.BeadType[btype];
    if (print) {
      fprintf(fw, "atom %7d", i);
      fprintf(fw, " name %8s", bt->Name);
      if (bt->Mass != MASS && bt->Mass != HIGHNUM) {
        fprintf(fw, " mass %15f ", bt->Mass);
      }
      if (bt->Charge != CHARGE && bt->Charge != HIGHNUM) {
        fprintf(fw, " charge %12f", bt->Charge);
      }
      if (bt->Radius != RADIUS && bt->Radius != HIGHNUM) {
        fprintf(fw, " radius %12f", bt->Radius);
      }
      if (mol != -1) {
        int mtype = System.Molecule[mol].Type;
        int id = System.Molecule[mol].Index;
        char name[8];
        if (snprintf(name, 8, "%s", System.MoleculeType[mtype].Name) < 0) {
          ErrorSnprintf();
        }
        fprintf(fw, " resname %10s", name);
        fprintf(fw, " resid %5d", id);
      }
      putc('\n', fw);
    }
  } //}}}
  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULE *mol = &System.Molecule[i];
    MOLECULETYPE *mt = &System.MoleculeType[mol->Type];
    if (mt->nBonds > 0) {
      fprintf(fw, "# resid %d\n", i + 1); // in VMD resid start with 1
      for (int j = 0; j < mt->nBonds; j++) {
        int id[2];
        for (int aa = 0; aa < 2; aa++) {
          id[aa] = mt->Bond[j][aa];
          id[aa] = mol->Bead[id[aa]];
        }
        fprintf(fw, "bond %6d: %6d\n", id[0], id[1]);
      }
    }
  } //}}}
  // print box size, if present //{{{
  if (box->Volume != -1 && StructureFileType(file) == VSF_FILE) {
    fprintf(fw, "pbc %.3f %.3f %.3f",
            box->Length[0], box->Length[1], box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 90) {
      fprintf(fw, " %lf %lf %lf", box->alpha, box->beta, box->gamma);
    }
    putc('\n', fw);
  } //}}}
  // close structure file
  fclose(fw);
} //}}}
// WriteLmpData() //{{{
static void WriteLmpData(SYSTEM System, char file[], bool mass,
                         int argc, char *argv[]) {
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "Created via AnalysisTools v%s ", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools); Command: ");
  PrintCommand(fw, argc, argv);
  putc('\n', fw);
  COUNT *Count = &System.Count;
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
  // // add one atom type for extra (possibly srp)
  // fprintf(fw, "%10d atom types\n", Count->BeadType + 1);
  fprintf(fw, "%10d atom types\n", mass_types);
  if (Count->Bond > 0 && Count->BondType > 0) {
    fprintf(fw, "%10d bond types\n", Count->BondType);
  }
  if (Count->Angle > 0 && Count->AngleType > 0) {
    fprintf(fw, "%10d angle types\n", Count->AngleType);
  }
  if (Count->Dihedral > 0 && Count->DihedralType > 0) {
    fprintf(fw, "%10d dihedral types\n", Count->DihedralType);
  }
  if (Count->Improper > 0 && Count->ImproperType > 0) {
    fprintf(fw, "%10d improper types\n", Count->ImproperType);
  }
  putc('\n', fw); //}}}
  // print box size //{{{
  fprintf(fw, "%.3f", System.Box.Low[0]);
  fprintf(fw, " %.3f xlo xhi\n", System.Box.Low[0] + System.Box.OrthoLength[0]);
  fprintf(fw, "%.3f", System.Box.Low[1]);
  fprintf(fw, " %.3f ylo yhi\n", System.Box.Low[1] + System.Box.OrthoLength[1]);
  fprintf(fw, "%.3f", System.Box.Low[2]);
  fprintf(fw, " %.3f zlo zhi\n", System.Box.Low[2] + System.Box.OrthoLength[2]);
  if (System.Box.alpha != 90 ||
      System.Box.beta != 90 ||
      System.Box.gamma != 90) {
    fprintf(fw, "%.3f %.3f %.3f xy xz yz\n", System.Box.transform[0][1],
                                             System.Box.transform[0][2],
                                             System.Box.transform[1][2]);
  }
  putc('\n', fw); //}}}
  // print bead type masses //{{{
  fprintf(fw, "Masses\n\n");
  for (int i = 0; i < mass_types; i++) {
    BEADTYPE *bt = &System.BeadType[bt_masstype_to_old[i]];
    fprintf(fw, "%5d", i + 1);
    if (bt->Mass == MASS) {
      fprintf(fw, " ???");
    } else {
      fprintf(fw, " %lf", bt->Mass);
    }
    fprintf(fw, " # %s\n", bt->Name);
  } //}}}
  // print various coeffs //{{{
  if (Count->BondType > 0) {
    fprintf(fw, "\nBond Coeffs\n\n");
    for (int i = 0; i < Count->BondType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.BondType[i].a / 2,
                                          System.BondType[i].b);
    }
  }
  if (Count->AngleType > 0) {
    fprintf(fw, "\nAngle Coeffs\n\n");
    for (int i = 0; i < Count->AngleType; i++) {
      fprintf(fw, "%5d %lf %lf\n", i + 1, System.AngleType[i].a / 2,
                                          System.AngleType[i].b);
    }
  }
  if (Count->DihedralType > 0) {
    fprintf(fw, "\nDihedral Coeffs\n\n");
    for (int i = 0; i < Count->DihedralType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.DihedralType[i].a / 2,
                                              System.DihedralType[i].b,
                                              System.DihedralType[i].c);
    }
  }
  if (Count->ImproperType > 0) {
    fprintf(fw, "\nImproper Coeffs\n\n");
    for (int i = 0; i < Count->ImproperType; i++) {
      fprintf(fw, "%5d %lf %lf %lf\n", i + 1, System.ImproperType[i].a / 2,
                                              System.ImproperType[i].b,
                                              System.ImproperType[i].c);
    }
  } //}}}
  // print atoms //{{{
  // if there is 0 molecule index, saved indices will get +1
  // TODO: why would I need to go from 1?
  bool zero = false;
  // for (int i = 0; i < Count->Molecule; i++) {
  //   if (System.Molecule[i].Index == 0) {
  //     zero = true;
  //     break;
  //   }
  // }
  fprintf(fw, "\nAtoms # full\n\n");
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    // <bead id>
    fprintf(fw, "%7d", id + 1);
    // <molecule id (-1 for no molecule)>
    int mol = bead->Molecule;
    if (mol != -1) {
      int id = System.Molecule[mol].Index;
      if (zero) {
        id++;
      }
      fprintf(fw, " %5d", id);
    } else {
      fprintf(fw, " %5d", -1);
    }
    // <bead type id>
    fprintf(fw, " %5d", bt_old_to_masstype[bead->Type] + 1);
    // <charge> from original System as the charge can differ in lmp data file
    int type = bead->Type;
    double q = System.BeadType[type].Charge;
    if (q == CHARGE) {
      fprintf(fw, " %15s", "???");
    } else {
      fprintf(fw, " %15f", q);
    }
    // coordinates
    for (int dd = 0; dd < 3; dd++) {
      fprintf(fw, " %15f", bead->Position[dd] + System.Box.Low[dd]);
    }
    // molecule name
    if (mol != -1) {
      int type = System.Molecule[mol].Type;
      fprintf(fw, " # %s", System.MoleculeType[type].Name);
    }
    putc('\n', fw);
  } //}}}
  // print velocities (if at least one non-zero) //{{{
  for (int i = 0; i < Count->BeadCoor; i++) {
    int id = System.BeadCoor[i];
    double (*vel)[3] = &System.Bead[id].Velocity;
    if (fabs((*vel)[0]) > 1e-5 ||
        fabs((*vel)[1]) > 1e-5 ||
        fabs((*vel)[2]) > 1e-5) {
      fprintf(fw, "\nVelocities\n\n");
      for (int j = 0; j < Count->BeadCoor; j++) {
        id = System.BeadCoor[j];
        vel = &System.Bead[id].Velocity;
        fprintf(fw, "%7d", id + 1);
        for (int dd = 0; dd < 3; dd++) {
          fprintf(fw, " %15f", (*vel)[dd]);
        }
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
        int id[2];
        for (int aa = 0; aa < 2; aa++) {
          id[aa] = System.MoleculeType[mtype].Bond[j][aa];
          id[aa] = System.Molecule[i].Bead[id[aa]];
        }
        if (System.Bead[id[0]].InTimestep &&
            System.Bead[id[1]].InTimestep) {
          fprintf(fw, "%7d", count);
          // if (Count->BondType > 0) {
          if (mt_i->Bond[j][2] != -1 ) {
            fprintf(fw, " %6d", mt_i->Bond[j][2] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d\n", id[0] + 1, id[1] + 1);
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
        int id[3];
        for (int aa = 0; aa < 3; aa++) {
          id[aa] = System.MoleculeType[mtype].Angle[j][aa];
          id[aa] = System.Molecule[i].Bead[id[aa]];
        }
        if (System.Bead[id[0]].InTimestep &&
            System.Bead[id[1]].InTimestep &&
            System.Bead[id[2]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->AngleType > 0) {
            fprintf(fw, " %6d", mt_i->Angle[j][3] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d\n", id[0] + 1, id[1] + 1, id[2] + 1);
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
        int id[4];
        for (int aa = 0; aa < 4; aa++) {
          id[aa] = System.MoleculeType[mtype].Dihedral[j][aa];
          id[aa] = System.Molecule[i].Bead[id[aa]];
        }
        if (System.Bead[id[0]].InTimestep &&
            System.Bead[id[1]].InTimestep &&
            System.Bead[id[2]].InTimestep &&
            System.Bead[id[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->DihedralType > 0) {
            fprintf(fw, " %6d", mt_i->Dihedral[j][4] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n",
                  id[0] + 1, id[1] + 1, id[2] + 1, id[3] + 1);
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
        int id[4];
        for (int aa = 0; aa < 4; aa++) {
          id[aa] = System.MoleculeType[mtype].Improper[j][aa];
          id[aa] = System.Molecule[i].Bead[id[aa]];
        }
        if (System.Bead[id[0]].InTimestep &&
            System.Bead[id[1]].InTimestep &&
            System.Bead[id[2]].InTimestep &&
            System.Bead[id[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count->DihedralType > 0) {
            fprintf(fw, " %6d", mt_i->Improper[j][4] + 1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n",
                  id[0] + 1, id[1] + 1, id[2] + 1, id[3] + 1);
        }
      }
    }
  } //}}}
  free(bt_masstype_to_old);
  free(bt_old_to_masstype);
  fclose(fw);
} //}}}
// WriteField() //{{{
static void WriteField(SYSTEM System, char file_field[],
                       int argc, char *argv[]) {
  FILE *fw = OpenFile(file_field, "w");
  BOX *box = &System.Box;
  if (box->Volume != -1) {
    fprintf(fw, "%.3f %.3f %.3f ",
            box->Length[0], box->Length[1], box->Length[2]);
    if (box->alpha != 90 || box->beta != 90 || box->gamma != 0) {
      fprintf(fw, "%lf %lf %lf ", box->alpha, box->beta, box->gamma);
    }
  }
  fprintf(fw, "Created via AnalysisTools v%s", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools); Command: ");
  PrintCommand(fw, argc, argv);
  putc('\n', fw);
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
      double (*pos)[3] = &System.Bead[id].Position;
      fprintf(fw, "%16s %8.5f %8.5f %8.5f\n", System.BeadType[bt].Name,
              (*pos)[0], (*pos)[1], (*pos)[2]);
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
                mt_i->Angle[j][1] + 1,
                mt_i->Angle[j][2] + 1);
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
                mt_i->Dihedral[j][1] + 1,
                mt_i->Dihedral[j][2] + 1,
                mt_i->Dihedral[j][3] + 1);
        int type = mt_i->Dihedral[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.DihedralType[type].a,
                  System.DihedralType[type].b,
                  System.DihedralType[type].c);
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
        fprintf(fw, "cvff %5d %5d %5d %5d", mt_i->Improper[j][0] + 1,
                mt_i->Improper[j][1] + 1,
                mt_i->Improper[j][2] + 1,
                mt_i->Improper[j][3] + 1);
        int type = mt_i->Improper[j][4];
        if (type != -1) {
          fprintf(fw, " %lf %lf %lf\n", System.ImproperType[type].a,
                  System.ImproperType[type].b,
                  System.ImproperType[type].c);
        } else {
          fprintf(fw, " ???\n");
        }
      }
    }
    fprintf(fw, "finish\n");
  } //}}}
  fclose(fw);
} //}}}
static void SimplifyResid(SYSTEM *System) { //{{{
  int lowest = 1e6;
  for (int i = 0; i < System->Count.Molecule; i++) {
    if (System->Molecule[i].Index < lowest) {
      lowest = System->Molecule[i].Index;
    }
  }
  for (int i = 0; i < System->Count.Molecule; i++) {
    System->Molecule[i].Index += -lowest + 1; // start from 1
  }
} //}}}

// InitCoorFile() //{{{
void InitCoorFile(FILE_TYPE file, SYSTEM System, int argc, char *argv[]) {
  if (file.type == VCF_FILE) {
    PrintByline(file.name, argc, argv);
  } else if (file.type == VTF_FILE) {
    WriteStructure(file, System, -1, false, argc, argv);
  } else {
    FILE *out = OpenFile(file.name, "w");
    fclose(out);
  }
} //}}}

// write structure and/or coordinates to a new file (can be any format) //{{{
void WriteOutput(SYSTEM System, const bool write[], FILE_TYPE fw,
                 bool lmp_mass, int vsf_def, int argc, char *argv[]) {
  if (fw.type == VCF_FILE) { // create vsf file if output file is vcf format
    PrintByline(fw.name, argc, argv); // byline to vcf file
    fw.name[strnlen(fw.name, LINE)-2] = 's';
    fw.type = VSF_FILE;
    WriteStructure(fw, System, vsf_def, lmp_mass, argc, argv);
    fw.name[strnlen(fw.name, LINE)-2] = 'c';
    fw.type = VCF_FILE;
  } else if (fw.type == VTF_FILE ||
    fw.type == VSF_FILE ||
    fw.type == FIELD_FILE ||
    fw.type == CONFIG_FILE ||
    fw.type == LDATA_FILE) {
    WriteStructure(fw, System, vsf_def, lmp_mass, argc, argv);
  }
  // write coordinates if the file is of coordinate type
  if (fw.type == VTF_FILE ||
    fw.type == VCF_FILE ||
    fw.type == XYZ_FILE ||
    fw.type == LTRJ_FILE) {
    // ensure it's a new file if the coordinate file is usually appended
    if (fw.type != VTF_FILE) {
      FILE *out = OpenFile(fw.name, "w");
      fclose(out);
    }
    WriteTimestep(fw, System, 1, write, argc, argv);
  }
}
void WriteOutputAll(SYSTEM System, FILE_TYPE fw, bool lmp_mass,
                    int vsf_def, int argc, char *argv[]) {
  bool *write = malloc(System.Count.Bead * sizeof *write);
  InitBoolArray(write, System.Count.Bead, true);
  WriteOutput(System, write, fw, lmp_mass, vsf_def, argc, argv);
  free(write);
} //}}}
// Write a single timestep to output file based on the file type //{{{
void WriteTimestep(FILE_TYPE f, SYSTEM System, int count_step,
                   const bool write[], int argc, char *argv[]) {
  FILE *fw = OpenFile(f.name, "a");
  switch (f.type) {
    case VCF_FILE:
    case VTF_FILE:
      VtfWriteCoorIndexed(fw, write, System);
      break;
    case XYZ_FILE:
      XyzWriteCoor(fw, write, System);
      break;
    case LTRJ_FILE:
      LtrjWriteCoor(fw, count_step, write, System);
      break;
    case LDATA_FILE:
      WriteLmpData(System, f.name, false, argc, argv);
      break;
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for output coor_type %s%d",
               ErrYellow(), f.type);
      PrintError();
      exit(1);
  }
  fclose(fw);
}
void WriteTimestepAll(FILE_TYPE f, SYSTEM System, int count_step,
                      int argc, char *argv[]) {
  bool *write = malloc(System.Count.Bead * sizeof *write);
  InitBoolArray(write, System.Count.Bead, true);
  WriteTimestep(f, System, count_step, write, argc, argv);
  free(write);
}
//}}}
// Create a structure file based on the file type (including dl_meso CONFIG) //{{{
void WriteStructure(FILE_TYPE f, SYSTEM System, int vsf_def_type,
                    bool lmp_mass, int argc, char *argv[]) {
  // ensure the output file is new
  switch (f.type) {
    case VSF_FILE:
    case VTF_FILE:
      VtfWriteStruct(f.name, System, vsf_def_type, argc, argv);
      break;
    case LDATA_FILE:
      WriteLmpData(System, f.name, lmp_mass, argc, argv);
      break;
    case CONFIG_FILE:
      WriteConfig(System, f.name);
      break;
    case FIELD_FILE:
      WriteField(System, f.name, argc, argv);
      break;
    case LTRJ_FILE:
      if (System.Count.BeadCoor == 0) {
        err_msg("no data to save into lammpstrj file (no coordinates loaded)");
        PrintError();
        exit(1);
      }
      bool *write = calloc(System.Count.Bead, sizeof *write);
      InitBoolArray(write, System.Count.Bead, true);
      FILE *fw = OpenFile(f.name, "w");
      LtrjWriteCoor(fw, 0, write, System);
      fclose(fw);
      free(write);
      break;
    default:
      err_msg("Inexistent output struct_type; should never happen!");
      PrintError();
      exit(1);
  }
} //}}}

// WriteAggregates() //{{{
void WriteAggregates(int step_count, char *agg_file, SYSTEM System,
                     AGGREGATE *Aggregate) {
  // get number of aggregates to write to agg_file
  int number_of_aggs = 0;
  for (int i = 0; i < System.Count.Aggregate; i++) {
    if (Aggregate[i].Flag) {
      number_of_aggs++;
    }
  }
  FILE *fw = OpenFile(agg_file, "a");
  // print number of aggregates to agg file
  fprintf(fw, "Step: %d\n%d\n", step_count, number_of_aggs);
  // go through all aggregates
  for (int i = 0; i < System.Count.Aggregate; i++) {
    // write only those that aren't excluded
    if (Aggregate[i].Flag) {
      // go through all molecules in aggregate 'i'
      fprintf(fw, "%d :", Aggregate[i].nMolecules);
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol = Aggregate[i].Molecule[j];
        fprintf(fw, " %d", System.Molecule[mol].Index);
      }
      putc('\n', fw);
    }
  }
  fclose(fw);
} //}}}

// verbose output (print various structures and some such)
void VerboseOutput(SYSTEM System) { //{{{
  PrintCount(System.Count);
  PrintBeadType(System);
  PrintAllMolTypes(System);
  PrintBondType(System);
  PrintAngleType(System);
  PrintDihedralType(System);
  PrintImproperType(System);
  if (System.Box.Volume != -1) {
    PrintBox(System.Box);
  }
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
  // if (Count.Molecule > 0) {
  //   fprintf(stdout, "  HighestResid:   %d", Count.HighestResid);
  // }
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
  int precision = 3;     // number of decimal digits
  int longest_name = 0;  // longest bead type name
  int max_number = 0;    // maximum number of beads
  int max_q = 0;         // maximum charge
  int max_m = 0;         // maximum mass
  int max_r = 0;         // maximum radius
  bool negative = false; // extra space for '-' if there's negative charge
  // determine length of values to have a nice-looking output
  for (int i = 0; i < System.Count.BeadType; i++) {
    BEADTYPE *bt = &System.BeadType[i];
    int length = strnlen(bt->Name, BEAD_NAME);
    if (length > longest_name) {
      longest_name = length;
    }
    if (bt->Number > max_number) {
      max_number = bt->Number;
    }
    if (bt->Charge < 0) {
      negative = true;
    }
    if (bt->Charge != CHARGE && bt->Charge != HIGHNUM && fabs(bt->Charge) > max_q) {
      max_q = floor(fabs(bt->Charge));
    }
    if (bt->Mass != MASS && bt->Mass != HIGHNUM && bt->Mass > max_m) {
      max_m = floor(bt->Mass);
    }
    if (bt->Radius != RADIUS && bt->Radius != HIGHNUM && bt->Radius > max_r) {
      max_r = floor(bt->Radius);
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
    BEADTYPE *bt = &System.BeadType[i];
    fprintf(stdout, "BeadType[%*d] = {", types_digits, i);
    fprintf(stdout, ".Name = %*s ", longest_name, bt->Name);
    fprintf(stdout, ".Number = %*d ", max_number, bt->Number);
    fprintf(stdout, ".Charge = ");
    if (bt->Charge != CHARGE && bt->Charge != HIGHNUM) {
      fprintf(stdout, "%*.*f ", max_q, precision, bt->Charge);
    } else {
      for (int j = 0; j < (max_q - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a ");
    }
    fprintf(stdout, ".Mass = ");
    if (bt->Mass != MASS && bt->Mass != HIGHNUM) {
      fprintf(stdout, "%*.*f ", max_m, precision, bt->Mass);
    } else {
      for (int j = 0; j < (max_m - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a ");
    }
    fprintf(stdout, ".Radius = ");
    if (bt->Radius != RADIUS && bt->Radius != HIGHNUM) {
      fprintf(stdout, "%*.*f", max_r, precision, bt->Radius);
    } else {
      for (int j = 0; j < (max_r - 3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a");
    }
    fprintf(stdout, " }\n");
  }
  putchar('\n');
} //}}}
void PrintOneMolType(SYSTEM System, int n) { //{{{
  int line = 80; // maximum printed line length
  MOLECULETYPE *mt = &System.MoleculeType[n];
  fprintf(stdout, "MoleculeType[%d] = {\n", n);
  fprintf(stdout, "  .Name       = %s\n", mt->Name);
  fprintf(stdout, "  .Number     = %d\n", mt->Number);
  // print bead types (list all beads) //{{{
  fprintf(stdout, "  .nBeads     = %d\n", mt->nBeads);
  int count = fprintf(stdout, "  .Bead       = {");
  for (int j = 0; j < mt->nBeads; j++) {
    count += fprintf(stdout, " %d", mt->Bead[j]);
    if (count >= line) {
      count = fprintf(stdout, "\n                 ") - 1;
    }
  }
  fprintf(stdout, " }\n"); //}}}
  // print bonds if there are any //{{{
  if (mt->nBonds > 0) {
    fprintf(stdout, "  .nBonds     = %d\n", mt->nBonds);
    count = fprintf(stdout, "  .Bond       = {");
    for (int j = 0; j < mt->nBonds; j++) {
      count += fprintf(stdout, " %d-%d", mt->Bond[j][0] + 1,
                       mt->Bond[j][1] + 1);
      if (mt->Bond[j][2] != -1) {
        count += fprintf(stdout, " (%d)", mt->Bond[j][2] + 1);
        if (j != (mt->nBonds - 1)) {
          putchar(',');
        }
      }
      if (count >= line) {
        count = fprintf(stdout, "\n                 ") - 1;
      }
    }
    fprintf(stdout, " }\n");
  } //}}}
  // print angles if there are any //{{{
  if (mt->nAngles > 0) {
    fprintf(stdout, "  .nAngles    = %d\n", mt->nAngles);
    count = fprintf(stdout, "  .Angle      = {");
    for (int j = 0; j < mt->nAngles; j++) {
      count += fprintf(stdout, " %d-%d-%d", mt->Angle[j][0] + 1,
                       mt->Angle[j][1] + 1,
                       mt->Angle[j][2] + 1);
      if (mt->Angle[j][3] != -1) {
        count += fprintf(stdout, " (%d)", mt->Angle[j][3] + 1);
        if (j != (mt->nAngles - 1)) {
          putchar(',');
        }
      }
      if (count >= 80) {
        count = fprintf(stdout, "\n                 ") - 1;
      }
    }
    fprintf(stdout, " }\n");
  } //}}}
  // print dihedrals if there are any //{{{
  if (mt->nDihedrals > 0) {
    fprintf(stdout, "  .nDihedrals = %d\n", mt->nDihedrals);
    count = fprintf(stdout, "  .Dihedral   = {");
    for (int j = 0; j < mt->nDihedrals; j++) {
      count += fprintf(stdout, " %d-%d-%d-%d", mt->Dihedral[j][0] + 1,
                       mt->Dihedral[j][1] + 1,
                       mt->Dihedral[j][2] + 1,
                       mt->Dihedral[j][3] + 1);
      if (mt->Dihedral[j][4] != -1) {
        count += fprintf(stdout, " (%d)", mt->Dihedral[j][4] + 1);
        if (j != (mt->nDihedrals - 1)) {
          putchar(',');
        }
      }
      if (count >= line) {
        count =fprintf(stdout, "\n                 ") - 1;
      }
    }
    fprintf(stdout, " }\n");
  } //}}}
  // print impropers if there are any //{{{
  if (mt->nImpropers > 0) {
    fprintf(stdout, "  .nImpropers = %d\n", mt->nImpropers);
    count = fprintf(stdout, "  .Improper   = {");
    for (int j = 0; j < mt->nImpropers; j++) {
      count += fprintf(stdout, " %d-%d-%d-%d", mt->Improper[j][0] + 1,
                       mt->Improper[j][1] + 1,
                       mt->Improper[j][2] + 1,
                       mt->Improper[j][3] + 1);
      if (mt->Improper[j][4] != -1) {
        count += fprintf(stdout, " (%d)", mt->Improper[j][4] + 1);
        if (j != (mt->nImpropers - 1)) {
          putchar(',');
        }
      }
      if (count >= line) {
        count = fprintf(stdout, "\n                 ") - 1;
      }
    }
    fprintf(stdout, " }\n");
  } //}}}
  // print bead types (just the which are present) //{{{
  fprintf(stdout, "  .nBTypes    = %d\n", mt->nBTypes);
  count = fprintf(stdout, "  .BType      = {");
  for (int j = 0; j < mt->nBTypes; j++) {
    count += fprintf(stdout, " %d", mt->BType[j]);
    if (count >= 80) {
      count = fprintf(stdout, "\n                 ") - 1;
    }
  }
  fprintf(stdout, " }\n"); //}}}
  if (mt->Mass != MASS) {
    fprintf(stdout, "  .Mass       = %.5f\n", mt->Mass);
  } else {
    fprintf(stdout, "  .Mass       = n/a\n");
  }
  if (mt->Charge != CHARGE) {
    fprintf(stdout, "  .Charge     = %.5f\n}\n", mt->Charge);
  } else {
    fprintf(stdout, "  .Charge     = n/a\n}\n");
  }
} //}}}
void PrintAllMolTypes(SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    PrintOneMolType(System, i);
  }
  if (System.Count.MoleculeType > 0) {
    putchar('\n');
  }
} //}}}
void Print1Molecule(SYSTEM System, int n) { //{{{
  MOLECULE *mol = &System.Molecule[n];
  MOLECULETYPE *mtype = &System.MoleculeType[mol->Type];
  fprintf(stdout, "Molecule %3d (%d, %s):\n", n + 1, mol->Index, mtype->Name);
  fprintf(stdout, " BEAD INDICES (%d): ", mtype->nBeads);
  fputs("intramolecular; input file\n", stdout);
  for (int j = 0; j < mtype->nBeads; j++) {
    fprintf(stdout, "   %3d; %5d\n", j + 1, mol->Bead[j]);
  }
} //}}}
void PrintMolecules(SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.Molecule; i++) {
    Print1Molecule(System, i);
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
    // TODO: eventually, there should be more than cvff
    // find highest value to get proper width //{{{
    PARAMS high = InitParams;
    for (int i = 0; i < System.Count.BondType; i++) {
      PARAMS *b = &System.BondType[i];
      if (b->a > high.a) {
        high.a = b->a;
      }
      if (b->b > high.b) {
        high.b = b->b;
      }
      if (b->c > high.c) {
        high.c = b->c;
      }
      if (b->d > high.d) {
        high.d = b->d;
      }
    }
    int wid[4];
    wid[0] = snprintf(NULL, 0, "%.5f", high.a);
    wid[1] = snprintf(NULL, 0, "%.0f", high.b);
    wid[2] = snprintf(NULL, 0, "%.0f", high.c);
    wid[3] = snprintf(NULL, 0, "%.0f", high.d);
    //}}}
    fprintf(stdout, "Bond types");
    fprintf(stdout, " (lammps style 'harm')\n");
    for (int i = 0; i < System.Count.BondType; i++) {
      PARAMS *b = &System.BondType[i];
      fprintf(stdout, "  %*.5f %*.5f\n", wid[0], b->a, wid[1], b->b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintAngleType(SYSTEM System) { //{{{
  if (System.Count.AngleType > 0) {
    // TODO: eventually, there should be more than cvff
    // find highest value to get proper width //{{{
    PARAMS high = InitParams;
    for (int i = 0; i < System.Count.AngleType; i++) {
      PARAMS *ang = &System.AngleType[i];
      if (ang->a > high.a) {
        high.a = ang->a;
      }
      if (ang->b > high.b) {
        high.b = ang->b;
      }
      if (ang->c > high.c) {
        high.c = ang->c;
      }
      if (ang->d > high.d) {
        high.d = ang->d;
      }
    }
    int wid[4];
    wid[0] = snprintf(NULL, 0, "%.5f", high.a);
    wid[1] = snprintf(NULL, 0, "%.0f", high.b);
    wid[2] = snprintf(NULL, 0, "%.0f", high.c);
    wid[3] = snprintf(NULL, 0, "%.0f", high.d);
    //}}}
    fprintf(stdout, "Angle types");
    fprintf(stdout, " (lammps style 'harm')\n");
    for (int i = 0; i < System.Count.AngleType; i++) {
      PARAMS *ang = &System.AngleType[i];
      fprintf(stdout, "  %*.5f %*.5f\n", wid[0], ang->a, wid[1], ang->b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintDihedralType(SYSTEM System) { //{{{
  if (System.Count.DihedralType > 0) {
    // TODO: eventually, there should be more than harm
    // find highest value to get proper width //{{{
    PARAMS high = InitParams;
    for (int i = 0; i < System.Count.DihedralType; i++) {
      PARAMS *dih = &System.DihedralType[i];
      if (dih->a > high.a) {
        high.a = dih->a;
      }
      if (dih->b > high.b) {
        high.b = dih->b;
      }
      if (dih->c > high.c) {
        high.c = dih->c;
      }
      if (dih->d > high.d) {
        high.d = dih->d;
      }
    }
    int wid[4];
    wid[0] = snprintf(NULL, 0, "%.5f", high.a);
    wid[1] = snprintf(NULL, 0, "%.0f", high.b);
    wid[2] = snprintf(NULL, 0, "%.0f", high.c);
    wid[3] = snprintf(NULL, 0, "%.0f", high.d);
    //}}}
    fprintf(stdout, "Dihedral types");
    // TODO: eventually, there should be more
    fprintf(stdout, " (lammps style 'harm')\n");
    for (int i = 0; i < System.Count.DihedralType; i++) {
      PARAMS *dih = &System.DihedralType[i];
      fprintf(stdout, "  %*.5f %*.0f %*.0f\n",
              wid[0], dih->a, wid[1], dih->b, wid[2], dih->c);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintImproperType(SYSTEM System) { //{{{
  if (System.Count.ImproperType > 0) {
    // TODO: eventually, there should be more than cvff
    // find highest value to get proper width //{{{
    PARAMS high = InitParams;
    for (int i = 0; i < System.Count.ImproperType; i++) {
      PARAMS *imp = &System.ImproperType[i];
      if (imp->a > high.a) {
        high.a = imp->a;
      }
      if (imp->b > high.b) {
        high.b = imp->b;
      }
      if (imp->c > high.c) {
        high.c = imp->c;
      }
      if (imp->d > high.d) {
        high.d = imp->d;
      }
    }
    int wid[4];
    wid[0] = snprintf(NULL, 0, "%.5f", high.a);
    wid[1] = snprintf(NULL, 0, "%.0f", high.b);
    wid[2] = snprintf(NULL, 0, "%.0f", high.c);
    wid[3] = snprintf(NULL, 0, "%.0f", high.d);
    //}}}
    fprintf(stdout, "Improper types");
    fprintf(stdout, " (lammps style 'cvff')\n");
    for (int i = 0; i < System.Count.ImproperType; i++) {
      PARAMS *imp = &System.ImproperType[i];
      fprintf(stdout, "  %*.5f %*.0f %*.0f\n",
              wid[0], imp->a, wid[1], imp->b, wid[2], imp->c);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintBox(BOX Box) { //{{{
  fprintf(stdout, "Box = {\n");
  if (Box.Low[0] != 0 || Box.Low[1] != 0 || Box.Low[2] != 0) {
    fprintf(stdout, "  .Low = ( %lf %lf %lf )\n",
            Box.Low[0], Box.Low[1], Box.Low[2]);
  }
  fprintf(stdout, "  .Length = ( %lf %lf %lf )\n",
          Box.Length[0], Box.Length[1], Box.Length[2]);
  if (Box.alpha != 0 || Box.beta != 90 || Box.gamma != 90) {
    fprintf(stdout, "  .alpha = %lf\n", Box.alpha);
    fprintf(stdout, "  .beta  = %lf\n", Box.beta);
    fprintf(stdout, "  .gamma = %lf\n", Box.gamma);
    fprintf(stdout, "  .OrthoLength = ( %lf %lf %lf )\n", Box.OrthoLength[0],
                                                          Box.OrthoLength[1],
                                                          Box.OrthoLength[2]);
    fprintf(stdout, "  .Bounding = ( %lf %lf %lf )\n", Box.Bounding[0],
                                                       Box.Bounding[1],
                                                       Box.Bounding[2]);
    fprintf(stdout, "  .transform = ( %lf %lf %lf)\n", Box.transform[0][0],
                                                       Box.transform[0][1],
                                                       Box.transform[0][2]);
    fprintf(stdout, "               ( %lf %lf %lf)\n", Box.transform[1][0],
                                                       Box.transform[1][1],
                                                       Box.transform[1][2]);
    fprintf(stdout, "               ( %lf %lf %lf)\n", Box.transform[2][0],
                                                       Box.transform[2][1],
                                                       Box.transform[2][2]);
  }
  fprintf(stdout, "  .Volume = %lf\n", Box.Volume);
  fprintf(stdout, "}\n");
} //}}}
void PrintByline(char file[], int argc, char *argv[]) { //{{{
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "# Created by AnalysisTools v%s ", VERSION);
  fprintf(fw, " (https://github.com/KaGaSi/AnalysisTools)\n");
  fprintf(fw, "# command: ");
  PrintCommand(fw, argc, argv);
  fclose(fw);
} //}}}
void PrintStep(int *count_coor, int start, bool silent) { //{{{
  (*count_coor)++;
  if (!silent && isatty(STDOUT_FILENO)) {
    if (*count_coor < start) {
      fprintf(stdout, "\rDiscarding step: %d", *count_coor);
    } else {
      if (*count_coor == start) {
        fprintf(stdout, "\rStarting step: %d    \n", start);
      }
      fprintf(stdout, "\rStep: %d", *count_coor);
    }
    fflush(stdout);
  }
} //}}}
void PrintAggregate(SYSTEM System, AGGREGATE Aggregate[]) { //{{{
  COUNT *Count = &System.Count;
  fprintf(stdout, "Aggregates: %d\n", Count->Aggregate);
  for (int i = 0; i < Count->Aggregate; i++) {
    // print molecules
    fprintf(stdout, " %d mols:", Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = Aggregate[i].Molecule[j];
      int type = System.Molecule[mol].Type;
      fprintf(stdout, " %d (%d)", mol, type);
      if (j != (Aggregate[i].nMolecules - 1)) {
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
  }
} //}}}
