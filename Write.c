#include "Write.h"

//VtfWriteCoorIndexed() //{{{
void VtfWriteCoorIndexed(FILE *vcf, char stuff[],
                         bool write[], SYSTEM System) {
  // print comment at the beginning of a timestep if present in initial vcf file
  if (stuff[0] != '\0') {
    fprintf(vcf, "%s\n", stuff);
  }
  // print box size (or 1 1 1 if unspecified box dimensions) //{{{
  BOX *box = &System.Box;
  if (box->Volume == -1) {
    strcpy(ERROR_MSG, "undefined simulation box size");
    PrintWarning();
    fprintf(stderr, "%sdimensions:%s %lf %lf %lf%s;", ErrCyan(), ErrYellow(),
            box->Length.x, box->Length.y, box->Length.z, ErrCyan());
    fprintf(stderr, " using 1 1 1 instead%s\n", ErrColourReset());
    fprintf(vcf, "pbc 1 1 1");
  } else {
    fprintf(vcf, "pbc %lf %lf %lf",
            box->Length.x, box->Length.y, box->Length.z);
  }
  fprintf(vcf, "    %lf %lf %lf\n", box->alpha, box->beta, box->gamma);
  //}}}
  fprintf(vcf, "indexed\n");
  bool none = true;
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    BEAD *bead = &System.Bead[id];
    if (bead->InTimestep && write[id]) {
      none = false;
      fprintf(vcf, "%8d %8.4f %8.4f %8.4f\n", id,
              bead->Position.x, bead->Position.y, bead->Position.z);
    }
  }
  if (none) {
    strcpy(ERROR_MSG, "no beads to save");
    PrintWarning();
  }
} //}}}
void XyzWriteCoor(FILE *xyz, bool write[], char *stuff, SYSTEM System) { //{{{
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
    WarnPrintWarning();
  } else {
    fprintf(xyz, "%d\n%s\n", count, stuff);
    for (int i = 0; i < System.Count.BeadCoor; i++) {
      int id = System.BeadCoor[i];
      BEAD *bead = &System.Bead[id];
      if (bead->InTimestep && write[id]) {
        int type = bead->Type;
        fprintf(xyz, "%8s %8.4f %8.4f %8.4f\n", System.BeadType[type].Name,
                bead->Position.x, bead->Position.y, bead->Position.z);
      }
    }
  }
} //}}}
void VtfWriteStruct(char file[], SYSTEM System, int type_def) { //{{{
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
    putc('\n', fw);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Count->Bead; i++) {
    int btype = System.Bead[i].Type,
        mol = System.Bead[i].Molecule;
    // print beads that are non-default, in a molecule, or have the highest id
    bool print = false;
    if (btype != type_def || mol != -1 || i == (Count->Bead-1)) {
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
    fprintf(fw, "# resid %d\n", i+1); // in VMD resid start with 1
    int mol_type = System.Molecule[i].Type;
    for (int j = 0; j < System.MoleculeType[mol_type].nBonds; j++) {
      int *bead = BondIndices(System, i, j);
      fprintf(fw, "bond %6d: %6d\n", bead[0], bead[1]);
    }
  } //}}}
  // close structure file
  fclose(fw);
} //}}}
void WriteLmpData(SYSTEM System, char file_lmp[], bool srp, bool mass) { //{{{
  FILE *fw = OpenFile(file_lmp, "w");
  fprintf(fw, "Created via AnalysisTools v%s \
  (https://github.com/KaGaSi/AnalysisTools)\n\n", VERSION);
  COUNT *Count = &System.Count;
  // create new SYSTEM structure to figure out bead types if mass==true //{{{
  SYSTEM Sys_print = CopySystem(System);
  COUNT *Count_print = &Sys_print.Count;
  if (mass) {
    Count_print->BeadType = 0;
    int *bt_old_to_new = calloc(Count->BeadType, sizeof *bt_old_to_new);
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt_i = &System.BeadType[i];
      bool found = false;
      for (int j = 0; j < Count_print->BeadType; j++) {
        BEADTYPE *bt_print_j = &Sys_print.BeadType[j];
        if (bt_i->Mass == bt_print_j->Mass) {
          bt_old_to_new[i] = j;
          found = true;
          break;
        }
      }
      if (!found) {
        BEADTYPE *bt_new = &Sys_print.BeadType[Count_print->BeadType];
        bt_new->Charge = bt_i->Charge;
        bt_new->Mass = bt_i->Mass;
        bt_new->Radius = bt_i->Radius;
        strcpy(bt_new->Name, bt_i->Name);
        bt_old_to_new[i] = Count_print->BeadType;
        Count_print->BeadType++;
      }
    }
    // change bead types to new ones
    for (int i = 0; i < Count_print->Bead; i++) {
      int old_type = Sys_print.Bead[i].Type;
      Sys_print.Bead[i].Type = bt_old_to_new[old_type];
    }
    for (int i = 0; i < Count_print->MoleculeType; i++) {
      for (int j = 0; j < Sys_print.MoleculeType[i].nBeads; j++) {
        int old_type = Sys_print.MoleculeType[i].Bead[j];
        Sys_print.MoleculeType[i].Bead[j] = bt_old_to_new[old_type];
      }
      free(Sys_print.MoleculeType[i].BType);
      FillMoleculeTypeBType(&Sys_print.MoleculeType[i]);
    }
    free(bt_old_to_new);
    // recount numbers of bead type - not really needed
    for (int i = 0; i < Count->BeadType; i++) {
      Sys_print.BeadType[i].Number = 0;
    }
    for (int i = 0; i < Count_print->Bead; i++) {
      int type = Sys_print.Bead[i].Type;
      Sys_print.BeadType[type].Number++;
    }
    // free extraneous Index arrays
    for (int i = 0; i < Count->BeadType; i++) {
      if (Sys_print.BeadType[i].Number == 0) {
        free(Sys_print.BeadType[i].Index);
      }
    }
  } //}}}
  // print counts //{{{
  fprintf(fw, "%10d atoms\n", Count_print->Bead);
  fprintf(fw, "%10d bonds\n", Count_print->Bond);
  fprintf(fw, "%10d angles\n", Count_print->Angle);
  fprintf(fw, "%10d dihedrals\n", Count_print->Dihedral);
  fprintf(fw, "%10d impropers\n", Count_print->Improper);
  putc('\n', fw);
  fprintf(fw, "%10d atom types\n", Count_print->BeadType);
  if (Count_print->Bond > 0 && Count_print->BondType == 0) {
    fprintf(fw,  "       ???");
  } else {
    fprintf(fw, "%10d", Count_print->BondType);
  }
  fprintf(fw, " bond types\n");
  if (Count_print->Angle > 0 && Count_print->AngleType == 0) {
    fprintf(fw,  "       ???");
  } else {
    fprintf(fw, "%10d", Count_print->AngleType);
  }
  fprintf(fw, " angle types\n");
  if (Count_print->Dihedral > 0 && Count_print->DihedralType == 0) {
    fprintf(fw,  "       ???");
  } else {
    fprintf(fw, "%10d", Count_print->DihedralType);
  }
  fprintf(fw, " dihedral types\n");
  if (Count_print->Improper > 0 && Count_print->ImproperType == 0) {
    fprintf(fw,  "       ???");
  } else {
    fprintf(fw, "%10d", Count_print->ImproperType);
  }
  fprintf(fw, " improper types\n");
  putc('\n', fw); //}}}
  // print box size //{{{
  fprintf(fw, "0 %lf xlo xhi\n", Sys_print.Box.OrthoLength.x);
  fprintf(fw, "0 %lf ylo yhi\n", Sys_print.Box.OrthoLength.y);
  fprintf(fw, "0 %lf zlo zhi\n", Sys_print.Box.OrthoLength.z);
  if (System.Box.alpha != 90 ||
      System.Box.beta != 90 ||
      System.Box.gamma != 90) {
    fprintf(fw, "%lf %lf %lf xy xz yz\n", Sys_print.Box.Tilt[0],
            Sys_print.Box.Tilt[1], Sys_print.Box.Tilt[2]);
  }
  putc('\n', fw); //}}}
  // print bead type masses //{{{
  fprintf(fw, "Masses\n\n");
  for (int i = 0; i < Count_print->BeadType; i++) {
    BEADTYPE *bt_i = &Sys_print.BeadType[i];
    fprintf(fw, "%5d %lf # %s\n", i+1, bt_i->Mass, bt_i->Name);
  }
  //}}}
  // print various coeffs //{{{
  fprintf(fw, "\nPair Coeffs\n\n");
  for (int i = 0; i < Count_print->BeadType; i++) {
    fprintf(fw, "%5d ???\n", i+1);
  }
  if (Count_print->Bond > 0) {
    fprintf(fw, "\nBond Coeffs\n\n");
    if (Count_print->BondType > 0) {
      for (int i = 0; i < Count_print->BondType; i++) {
        fprintf(fw, "%5d %lf %lf\n", i+1,
                Sys_print.BondType[i].a/2, Sys_print.BondType[i].b);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count_print->Angle > 0) {
    fprintf(fw, "\nAngle Coeffs\n\n");
    if (Count_print->AngleType > 0) {
      for (int i = 0; i < Count_print->AngleType; i++) {
        fprintf(fw, "%5d %lf %lf\n", i+1,
                Sys_print.AngleType[i].a/2, Sys_print.AngleType[i].b);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count_print->Dihedral > 0) {
    fprintf(fw, "\nDihedral Coeffs\n\n");
    if (Count_print->DihedralType > 0) {
      for (int i = 0; i < Count_print->DihedralType; i++) {
        fprintf(fw, "%5d %lf %lf %lf\n", i+1, Sys_print.DihedralType[i].a/2,
                Sys_print.DihedralType[i].b, Sys_print.DihedralType[i].c);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  }
  if (Count_print->Improper > 0) {
    fprintf(fw, "\nImproper Coeffs\n\n");
    if (Count_print->ImproperType > 0) {
      for (int i = 0; i < Count_print->ImproperType; i++) {
        fprintf(fw, "%5d %lf %lf %lf\n", i+1, Sys_print.ImproperType[i].a/2,
                Sys_print.ImproperType[i].b, Sys_print.ImproperType[i].c);
      }
    } else {
      fprintf(fw, "  ???\n");
    }
  } //}}}
  // print atoms //{{{
  fprintf(fw, "\nAtoms\n\n");
  for (int i = 0; i < Count_print->BeadCoor; i++) {
    int id = Sys_print.BeadCoor[i];
    BEAD *bead = &Sys_print.Bead[id];
    // <bead id>
    fprintf(fw, "%7d", id+1);
    // <molecule id (-1 for no molecule)>
    int mol = bead->Molecule;
    if (mol != -1) {
      fprintf(fw, "%5d", Sys_print.Molecule[mol].Index);
    } else {
      fprintf(fw, "%5d", -1);
    }
    // <bead type id>
    fprintf(fw, "%5d", bead->Type+1);
    // <charge> from original System as the charge can differ in lmp data file
    int type = System.Bead[id].Type;
    fprintf(fw, "%15f", System.BeadType[type].Charge);
    // coordinates
    fprintf(fw, "%15f %15f %15f", bead->Position.x,
            bead->Position.y, bead->Position.z);
    putc('\n', fw);
  } //}}}
  // print velocities (if at least one non-zero) //{{{
  for (int i = 0; i < Count_print->BeadCoor; i++) {
    int id = Sys_print.BeadCoor[i];
    VECTOR *vel = &Sys_print.Bead[id].Velocity;
    if (fabs(vel->x) > 1e-5 || fabs(vel->y) > 1e-5 || fabs(vel->z) > 1e-5) {
      fprintf(fw, "\nVelocities\n\n");
      for (int j = 0; j < Count_print->BeadCoor; j++) {
        id = Sys_print.BeadCoor[j];
        vel = &Sys_print.Bead[id].Velocity;
        fprintf(fw, "%7d", id+1);
        fprintf(fw, "%15f %15f %15f", vel->x, vel->y, vel->z);
        putc('\n', fw);
      }
      break;
    }
  } //}}}
  // print bonds //{{{
  if (Count_print->Bond > 0) {
    fprintf(fw, "\nBonds\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = Sys_print.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nBonds; j++) {
        count++;
        int *val = BondIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep &&
            System.Bead[val[1]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count_print->BondType > 0) {
            fprintf(fw, "%6d", mt_i->Bond[j][2]+1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d\n", val[0]+1, val[1]+1);
        }
      }
    }
  } //}}}
  // print angles //{{{
  if (Count_print->Angle > 0) {
    fprintf(fw, "\nAngles\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = Sys_print.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nAngles; j++) {
        count++;
        int *val = AngleIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep &&
            System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count_print->AngleType > 0) {
            fprintf(fw, "%6d", mt_i->Angle[j][3]+1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d\n", val[0]+1, val[1]+1, val[2]+1);
        }
      }
    }
  }
  //}}}
  // print dihedrals //{{{
  if (Count_print->Dihedral > 0) {
    fprintf(fw, "\nDihedrals\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = Sys_print.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nDihedrals; j++) {
        count++;
        int *val = DihedralIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep &&
            System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep &&
            System.Bead[val[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count_print->DihedralType > 0) {
            fprintf(fw, "%6d", mt_i->Dihedral[j][4]+1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n",
                  val[0]+1, val[1]+1, val[2]+1, val[3]+1);
        }
      }
    }
  } //}}}
  // print impropers //{{{
  if (Count_print->Improper > 0) {
    fprintf(fw, "\nImpropers\n\n");
    int count = 0;
    for (int i = 0; i < Count->Molecule; i++) {
      int mtype = Sys_print.Molecule[i].Type;
      MOLECULETYPE *mt_i = &System.MoleculeType[mtype];
      for (int j = 0; j < mt_i->nImpropers; j++) {
        count++;
        int *val = ImproperIndices(System, i, j);
        if (System.Bead[val[0]].InTimestep &&
            System.Bead[val[1]].InTimestep &&
            System.Bead[val[2]].InTimestep &&
            System.Bead[val[3]].InTimestep) {
          fprintf(fw, "%7d", count);
          if (Count_print->DihedralType > 0) {
            fprintf(fw, "%6d", mt_i->Dihedral[j][4]+1);
          } else {
            fprintf(fw, "   ???");
          }
          fprintf(fw, " %5d %5d %5d %5d\n",
                  val[0]+1, val[1]+1, val[2]+1, val[3]+1);
        }
      }
    }
  } //}}}
  FreeSystem(&Sys_print);
  fclose(fw);
} //}}}
void WriteField(SYSTEM System, char file_field[]) { //{{{
  FILE *fw = OpenFile(file_field, "w");
  fprintf(fw, "Created via AnalysisTools v%s \
  (https://github.com/KaGaSi/AnalysisTools)\n\n", VERSION);
  COUNT *Count = &System.Count;
  // print species section //{{{
  fprintf(fw, "species %d\n", Count->BeadType);
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
    fprintf(fw, "%16s %8.5f %8.5f %5d\n",
            bt_i->Name, bt_i->Mass, bt_i->Charge, unbonded[i]);
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
      fprintf(fw, "%16s %8.5f %8.5f %8.5f\n",
              System.BeadType[bt].Name, pos->x, pos->y, pos->z);
    }
    // bonds (if present)
    if (mt_i->nBonds > 0) {
      fprintf(fw, "bonds %d\n", mt_i->nBonds);
      for (int j = 0; j < mt_i->nBonds; j++) {
        // TODO harm only for now
        fprintf(fw, "harm %5d %5d", mt_i->Bond[j][0], mt_i->Bond[j][1]);
        int type = mt_i->Bond[j][2];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n",
                  System.BondType[type].a, System.BondType[type].b);
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
        fprintf(fw, "harm %5d %5d %5d",
                mt_i->Angle[j][0], mt_i->Angle[j][1], mt_i->Angle[j][2]);
        int type = mt_i->Angle[j][3];
        if (type != -1) {
          fprintf(fw, " %lf %lf\n",
                  System.AngleType[type].a, System.AngleType[type].b);
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
        fprintf(fw, "harm %5d %5d %5d %5d", mt_i->Dihedral[j][0],
                                            mt_i->Dihedral[j][1],
                                            mt_i->Dihedral[j][2],
                                            mt_i->Dihedral[j][3]);
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
        fprintf(fw, "cvff %5d %5d %5d %5d", mt_i->Improper[j][0],
                                            mt_i->Improper[j][1],
                                            mt_i->Improper[j][2],
                                            mt_i->Improper[j][3]);
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

// TODO will change
// WriteAggregates() //{{{
/*
* Append aggregate information from a timestep to an .agg file
*/
void WriteAggregates(int step_count, char *agg_file, COUNTS Counts,
                   MOLECULETYPE *MoleculeType, BEAD *Bead, AGGREGATE *Aggregate) {

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
        fprintf(fw, " %d", Aggregate[i].Molecule[j]+1);
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

#if 0 //{{{
// TODO will change
// WriteField() //{{{
/*
 * Create a new dl_meso FIELD file
 */
void WriteField(char *field, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                PARAMS *bond_type, PARAMS *angle_type, PARAMS *dihedral_type) {
  // count unbonded beads for each bead type
  int unbonded[Counts.TypesOfBeads];
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    unbonded[i] = 0;
  }
  for (int i = 0; i < Counts.Unbonded; i++) {
    if (Bead[i].Molecule == -1) {
      int btype = Bead[i].Type;
      unbonded[btype]++;
    }
  }

  FILE *fw = OpenFile(field, "w");
  fprintf(fw, "Created by AnalysisTools version %s ", VERSION);
  fprintf(fw, "(https://github.com/KaGaSi/AnalysisTools/releases)\n");
  fprintf(fw, "\nspecies %d\n", Counts.TypesOfBeads);
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(fw, "%s %lf %lf %d\n", BeadType[i].Name,
                                   BeadType[i].Mass,
                                   BeadType[i].Charge,
                                   unbonded[i]);
  }
  fprintf(fw, "\nmolecule %d\n", Counts.TypesOfMolecules);
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(fw, "%s\n", MoleculeType[i].Name);
    fprintf(fw, "nummols %d\n", MoleculeType[i].Number);
    // error - no beads //{{{
    if (MoleculeType[i].nBeads < 1) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - molecule ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", MoleculeType[i].Name);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " contains no beads\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    } //}}}
    // find first molecule of the given type
    int mol = -1;
    for (int j = 0; j < Counts.Molecules; j++) {
      if (Molecule[j].Type == i) {
        mol = j;
        break;
      }
    }
    // error - no molecule of given type //{{{
    if (MoleculeType[i].Number < 1 || mol == -1) {
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " - molecule ");
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", MoleculeType[i].Name);
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " does not exist\n\n");
      ColourReset(STDERR_FILENO);
      exit(1);
    } //}}}
    /*
     * print beads with coordinates of the first molecule, changing them so
     * that the first bead is (0,0,0)
     */
    fprintf(fw, "beads %d\n", MoleculeType[i].nBeads);
    int first = Molecule[mol].Bead[0];
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      int btype = MoleculeType[i].Bead[j];
      int bead = Molecule[mol].Bead[j];
      VECTOR pos;
      pos.x = Bead[bead].Position.x - Bead[first].Position.x;
      pos.y = Bead[bead].Position.y - Bead[first].Position.y;
      pos.z = Bead[bead].Position.z - Bead[first].Position.z;
      fprintf(fw, "%s %lf %lf %lf\n", BeadType[btype].Name,
                                      pos.x, pos.y, pos.z);
    }
    if (MoleculeType[i].nBonds > 0) {
      fprintf(fw, "bonds %d\n", MoleculeType[i].nBonds);
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        fprintf(fw, "harm %d %d ", MoleculeType[i].Bond[j][0]+1,
                                   MoleculeType[i].Bond[j][1]+1);
        int bond = MoleculeType[i].Bond[j][2];
        if (bond != -1) {
          fprintf(fw, "%lf %lf\n", bond_type[bond].a, bond_type[bond].b);
        } else {
          fprintf(fw, "0 0\n");
        }
      }
    }
    if (MoleculeType[i].nAngles > 0) {
      fprintf(fw, "angles %d\n", MoleculeType[i].nAngles);
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        fprintf(fw, "harm %d %d %d ", MoleculeType[i].Angle[j][0],
                                      MoleculeType[i].Angle[j][1],
                                      MoleculeType[i].Angle[j][2]);
        int angle = MoleculeType[i].Angle[j][3];
        if (angle != -1) {
          fprintf(fw, "%lf %lf\n", angle_type[angle].a, angle_type[angle].b);
        } else {
          fprintf(fw, "0 0\n");
        }
      }
    }
    if (MoleculeType[i].nDihedrals > 0) {
      fprintf(fw, "dihedrals %d\n", MoleculeType[i].nDihedrals);
      for (int j = 0; j < MoleculeType[i].nDihedrals; j++) {
        fprintf(fw, "harm %d %d %d %d ", MoleculeType[i].Dihedral[j][0],
                                         MoleculeType[i].Dihedral[j][1],
                                         MoleculeType[i].Dihedral[j][2],
                                         MoleculeType[i].Dihedral[j][3]);
        int dihedral = MoleculeType[i].Dihedral[j][4];
        if (dihedral != -1) {
          fprintf(fw, "%lf %lf\n", dihedral_type[dihedral].a,
                                   dihedral_type[dihedral].b);
        } else {
          fprintf(fw, "0 0\n");
        }
      }
    }
    fprintf(fw, "finish\n");
  }
  fclose(fw);
} //}}}
// TODO remove
// WriteCoorIndexed() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Write flag
 * in BeadType and MoleculeType structures only certain bead types will be
 * saved into the indexed timestep in .vcf file.
 */
void WriteCoorIndexed(FILE *vcf_file, COUNTS Counts,
                         BEADTYPE *BeadType, BEAD *Bead,
                         MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                         char *stuff, BOX Box) {
  // print comment at the beginning of a timestep if present in initial vcf file
  fprintf(vcf_file, "%s\n", stuff);
  // print box size
  fprintf(vcf_file, "pbc %lf %lf %lf  ", Box.Length.x,
                                         Box.Length.y,
                                         Box.Length.z);
  fprintf(vcf_file, "    %lf %lf %lf\n", Box.alpha, Box.beta, Box.gamma);
  // print 'indexed' on the next
  fprintf(vcf_file, "indexed\n");

  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int btype = Bead[i].Type;
    if (BeadType[btype].Write) {
      if (Bead[i].Molecule != -1) { // bead in a molecule
        int mtype = Molecule[Bead[i].Molecule].Type;
        if (MoleculeType[mtype].Write) {
          fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", i,
                                                       Bead[i].Position.x,
                                                       Bead[i].Position.y,
                                                       Bead[i].Position.z);
        }
      } else { // monomer bead
        fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", i,
                                                     Bead[i].Position.x,
                                                     Bead[i].Position.y,
                                                     Bead[i].Position.z);
      }
    }
  }
} //}}}
// WriteVsf_old //{{{
void WriteVsf_old(char *input_vsf, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
              MOLECULETYPE *MoleculeType, MOLECULE *Molecule, bool change) {

  FILE *fw = OpenFile(input_vsf, "w");
  // TODO integrate InFile & BeadsCoor
  // find most common type of bead and make it default //{{{
  int type_def = -1, count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    bool use = true;
    for (int j = 0; j < Counts.TypesOfMolecules; j++) {
      for (int k = 0; k < MoleculeType[j].nBTypes; k++) {
        if (MoleculeType[j].BType[k] == i) {
          use = false;
        }
      }
    }
    if (use && BeadType[i].Number >= count) {
      count = BeadType[i].Number;
      type_def = i;
    }
  } //}}}
  // print default bead type //{{{
  if (type_def != -1) {
    fprintf(fw, "atom default name %8s ", BeadType[type_def].Name);
    fprintf(fw, "mass %lf ", BeadType[type_def].Mass);
    fprintf(fw, "charge %lf\n", BeadType[type_def].Charge);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int btype = Bead[i].Type;
    int mol = Bead[i].Molecule;
    // don't print beads with type 'type_def'
    if (btype != type_def || mol != -1) {
      fprintf(fw, "atom %7d ", i);
      if (mol != -1 && change) {
        int mtype = Molecule[mol].Type;
        int n = -1;
        for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
          if (i == Molecule[mol].Bead[j]) {
            n = j;
            break;
          }
        }
        btype = MoleculeType[mtype].Bead[n];
        fprintf(fw, "name %8s ", BeadType[btype].Name);
        fprintf(fw, "mass %lf ", BeadType[btype].Mass);
        fprintf(fw, "charge %lf", BeadType[btype].Charge);
      } else {
        fprintf(fw, "name %8s ", BeadType[btype].Name);
        fprintf(fw, "mass %lf ", BeadType[btype].Mass);
        fprintf(fw, "charge %lf", BeadType[btype].Charge);
      }
      if (mol != -1) {
        int mtype = Molecule[mol].Type;
        fprintf(fw, " resname %10s ", MoleculeType[mtype].Name);
        fprintf(fw, "resid %5d", mol+1);
      }
      putc('\n', fw);
    // print highest bead id even if it's default type
    } else if (i == (Counts.BeadsTotal-1)) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[btype].Name);
      fprintf(fw, "mass %lf ", BeadType[btype].Mass);
      fprintf(fw, "charge %lf", BeadType[btype].Charge);
      putc('\n', fw);
    }
  } //}}}
  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Counts.Molecules; i++) {
    fprintf(fw, "# resid %d\n", i+1); // in VMD resid start with 1
    int mol_type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mol_type].nBonds; j++) {
      int bead1 = MoleculeType[mol_type].Bond[j][0];
      int bead2 = MoleculeType[mol_type].Bond[j][1];
      bead1 = Molecule[i].Bead[bead1];
      bead2 = Molecule[i].Bead[bead2];
      fprintf(fw, "bond %6d: %6d\n", bead1, bead1);
    }
  } //}}}
  // close structure file
  fclose(fw);
} //}}}
// Append an indexed timestep to a vcf/vtf coordinate file //{{{
void VtfWriteCoorIndexed_old(FILE *vcf, char stuff[], int InFile[],
                         COUNTS Counts, BEAD Bead[], BOX Box) {
  // print comment at the beginning of a timestep if present in initial vcf file
  if (stuff[0] != '\0') {
    fprintf(vcf, "%s\n", stuff);
  }
  // print box size
  fprintf(vcf, "pbc %lf %lf %lf  ", Box.Length.x, Box.Length.y, Box.Length.z);
  fprintf(vcf, "    %lf %lf %lf\n", Box.alpha, Box.beta, Box.gamma);
  // print 'indexed' on the next
  fprintf(vcf, "indexed\n");

  bool none = true;
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int id = InFile[i];
    if (Bead[id].InTimestep) {
      none = false;
      fprintf(vcf, "%8d %8.4f %8.4f %8.4f\n", id,
                                              Bead[id].Position.x,
                                              Bead[id].Position.y,
                                              Bead[id].Position.z);
    }
  }
  if (none) {
    strcpy(ERROR_MSG, "no beads to save");
    PrintWarning();
  }
} //}}}
// Append a timestep to an xyz file //{{{
void XyzWriteCoor_old(FILE *xyz, COUNTS Counts, int InFile[],
                  BEADTYPE *BeadType, BEAD *Bead) {
  // find out number of beads to save
  int count = 0;
  bool none = true; // to make sure there are beads to save
  PrintCounts_old(Counts);
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int id = InFile[i];
    if (Bead[id].InTimestep) {
      none = false;
      count++;
    }
  }
  if (none) {
    strcpy(ERROR_MSG, "no beads to save");
    WarnPrintWarning();
  } else {
    fprintf(xyz, "%d\n\n", count);
    for (int i = 0; i < Counts.BeadsCoor; i++) {
      int id = InFile[i];
      if (Bead[id].InTimestep) {
        int type = Bead[id].Type;
        fprintf(xyz, "%8s %7.3f %7.3f %7.3f\n", BeadType[type].Name,
                                                Bead[id].Position.x,
                                                Bead[id].Position.y,
                                                Bead[id].Position.z);
      }
    }
  }
} //}}}
// Create a new vsf/vtf structure file. //{{{
void VtfWriteStruct_old(char file[], COUNTS Counts,
                    BEADTYPE BeadType[], BEAD Bead[],
                    MOLECULETYPE MoleculeType[], MOLECULE Molecule[]) {

  FILE *fw = OpenFile(file, "w");
  // TODO integrate InFile & BeadsCoor
  // find most common type of bead and make it default //{{{
  int *count = calloc(Counts.TypesOfBeads, sizeof *count);
  for (int i = 0; i < Counts.BeadsTotal; i++) {
    if (Bead[i].Molecule == -1) {
      int type = Bead[i].Type;
      count[type]++;
    }
  }
  int type_def = -1, max = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (count[i] > max) {
      max = count[i];
      type_def = i;
    }
  }
  PrintBeadType2(Counts.TypesOfBeads, BeadType);
  free(count); //}}}
  // print default bead type //{{{
  if (type_def != -1) {
    fprintf(fw, "atom default name %8s ", BeadType[type_def].Name);
    fprintf(fw, "mass %lf ", BeadType[type_def].Mass);
    fprintf(fw, "charge %lf\n", BeadType[type_def].Charge);
  } //}}}
  // print beads //{{{
  for (int i = 0; i < Counts.BeadsTotal; i++) {
    int btype = Bead[i].Type;
    int mol = Bead[i].Molecule;
    // don't print beads with type 'type_def'
    if (btype != type_def || mol != -1) {
      fprintf(fw, "atom %7d ", i);
      if (mol != -1) {
        int mtype = Molecule[mol].Type;
        int n = -1;
        for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
          if (i == Molecule[mol].Bead[j]) {
            n = j;
            break;
          }
        }
        btype = MoleculeType[mtype].Bead[n];
        fprintf(fw, "name %8s ", BeadType[btype].Name);
        fprintf(fw, "mass %lf ", BeadType[btype].Mass);
        fprintf(fw, "charge %lf", BeadType[btype].Charge);
      } else {
        fprintf(fw, "name %8s ", BeadType[btype].Name);
        fprintf(fw, "mass %lf ", BeadType[btype].Mass);
        fprintf(fw, "charge %lf", BeadType[btype].Charge);
      }
      if (mol != -1) {
        int mtype = Molecule[mol].Type;
        fprintf(fw, " resname %10s ", MoleculeType[mtype].Name);
        fprintf(fw, "resid %5d", mol+1);
      }
      putc('\n', fw);
    // print highest bead id even if it's default type
    } else if (i == (Counts.BeadsTotal-1)) {
      fprintf(fw, "atom %7d ", i);
      fprintf(fw, "name %8s ", BeadType[btype].Name);
      fprintf(fw, "mass %lf ", BeadType[btype].Mass);
      fprintf(fw, "charge %lf", BeadType[btype].Charge);
      putc('\n', fw);
    }
  } //}}}
  // print bonds //{{{
  putc('\n', fw);
  for (int i = 0; i < Counts.Molecules; i++) {
    fprintf(fw, "# resid %d\n", i+1); // in VMD resid start with 1
    int mol_type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mol_type].nBonds; j++) {
      int bead1 = MoleculeType[mol_type].Bond[j][0];
      int bead2 = MoleculeType[mol_type].Bond[j][1];
      bead1 = Molecule[i].Bead[bead1];
      bead2 = Molecule[i].Bead[bead2];
      fprintf(fw, "bond %6d: %6d\n", bead1, bead2);
    }
  } //}}}
  // close structure file
  fclose(fw);
} //}}}
#endif //}}}
