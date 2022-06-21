#include "AnalysisTools.h"
/* TODO: PrintBeadType() needs to find the longest name and the most number
 * of beads and insert white space accordingly
*/
// TODO add detailed switch FindBeadType/FindMoleculeType (check not just name)

// TODO make CheckSystem & CheckSystemCoor?
void CheckSystem(SYSTEM System, char file[]) {
  COUNT *Count = &System.Count;
  // count total number of beads //{{{
  // i) just unbonded+bonded
  int count = Count->Unbonded + Count->Bonded;
  if (count != Count->Bead) {
    strcpy(ERROR_MSG, "unbonded and bonded beads do not add up properly!");
    PrintError();
    ErrorPrintFile(file, "\0");
    fprintf(stderr, "%s, unbonded: %s%d%s",
            ErrRed(), ErrYellow(), Count->Unbonded, ErrRed());
    fprintf(stderr, ", bonded: %s%d%s",
            ErrYellow(), Count->Bonded, ErrRed());
    fprintf(stderr, ", sum should be: %s%d%s\n",
            ErrYellow(), Count->Bead, ErrColourReset());
  }
  // ii) from BeadType
  count = 0;
  for (int i = 0; i < Count->BeadType; i++) {
    count += System.BeadType[i].Number;
  }
  if (count != Count->Bead) {
    strcpy(ERROR_MSG, "number of beads in bead types do not add up properly!");
    PrintError();
    ErrorPrintFile(file, "\0");
    fprintf(stderr, "%s, sum should be: %s%d%s",
            ErrRed(), ErrYellow(), Count->Bead, ErrRed());
    fprintf(stderr, " but is: %s%d%s\n",
            ErrYellow(), count, ErrRed());
    fprintf(stderr, "Number | Name%s\n", ErrYellow());
    for (int i = 0; i < Count->BeadType; i++) {
      fprintf(stderr, "%6d %s|%s %s\n", System.BeadType[i].Number, ErrRed(),
              ErrYellow(), System.BeadType[i].Name);
    }
  } //}}}
}
int * BondIndices(SYSTEM System, int mol, int bond) { //{{{
  static int index[4]; // first are 'real' ids, then the intramolecular ones
  int n = 2,
      type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Bond[bond][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in a bond; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),",
            ErrRed(), ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in bond %s%d%s:%s",
            ErrYellow(), bond, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
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
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s),",
            ErrRed(), ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed());
    fprintf(stderr, " indices in bond %s%d%s:%s",
            ErrYellow(), bond, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int * AngleIndices(SYSTEM System, int mol, int angle) { //{{{
  static int index[6]; // first are 'real' ids, then the intramolecular ones
  int n = 3,
      type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Angle[angle][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an angle; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),",
            ErrRed(), ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in angle %s%d%s:%s",
            ErrYellow(), angle, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
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
    fprintf(stderr, " indices in angle %s%d%s:%s",
            ErrYellow(), angle, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}
int * DihedralIndices(SYSTEM System, int mol, int dihed) { //{{{
  static int index[8]; // first are 'real' ids, then the intramolecular ones
  int n = 4,
      type = System.Molecule[mol].Type;
  // intramolecular is
  for (int i = 0; i < n; i++) {
    index[i+n] = System.MoleculeType[type].Dihedral[dihed][i];
  }
  // error if wrong intramolecular id //{{{
  bool err = false;
  for (int i = 0; i < n; i++) {
    if (index[i+n] < 0 || index[i+n] >= System.MoleculeType[type].nBeads) {
      err = true;
      break;
    }
  }
  if (err) {
    strcpy(ERROR_MSG, "wrong index in an dihedral; should never happen!");
    PrintError();
    fprintf(stderr, "%sMolecule %s%d%s (%s%s%s with %s%d%s beads),",
            ErrRed(), ErrYellow(), mol, ErrRed(), ErrYellow(),
            System.MoleculeType[type].Name, ErrRed(), ErrYellow(),
            System.MoleculeType[type].nBeads, ErrRed());
    fprintf(stderr, " intramolecular indices in dihedral %s%d%s:%s",
            ErrYellow(), dihed, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i+n]);
    }
    fprintf(stderr, "%s\n", ErrColourReset());
  } //}}}
  for (int i = 0; i < n; i++) {
    index[i] = System.Molecule[mol].Bead[index[i+n]];
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
    fprintf(stderr, " indices in dihedral %s%d%s:%s",
            ErrYellow(), dihed, ErrRed(), ErrYellow());
    for (int i = 0; i < n; i++) {
      fprintf(stderr, " %d", index[i]);
    }
    fprintf(stderr, " %s(%s%d%s beads in the system)%s\n", ErrRed(),
            ErrYellow(), System.Count.Bead, ErrRed(), ErrColourReset());
  } //}}}
  return index;
} //}}}

/* HOW TO CALCULATE DISTANCE IN TRICLINIC SYSTEM //{{{
//VECTOR dist;
//dist.x = (*Bead)[0].Position.x - (*Bead)[10].Position.x;
//dist.y = (*Bead)[0].Position.y - (*Bead)[10].Position.y;
//dist.z = (*Bead)[0].Position.z - (*Bead)[10].Position.z;
//printf("dist1 = (%lf, %lf, %lf) = %lf\n", dist.x, dist.y, dist.z, sqrt(SQR(dist.x)+SQR(dist.y)+SQR(dist.z)));

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
//printf("dist2 = (%lf, %lf, %lf) = %lf\n", dist.x, dist.y, dist.z, sqrt(SQR(dist.x)+SQR(dist.y)+SQR(dist.z)));
*/ //}}}

// TriclinicCellData() //{{{
bool TriclinicCellData(BOX *Box) {
  // triclinic box //{{{
  if (Box->alpha != 90 || Box->beta != 90 || Box->gamma != 90 ) {
    double a = Box->Length.x,
           b = Box->Length.y,
           c = Box->Length.z;
    double c_a = cos(Box->alpha * PI / 180),
           c_b = cos(Box->beta * PI / 180),
           c_g = cos(Box->gamma * PI / 180),
           s_g = sin(Box->gamma * PI / 180);
    double sqr = 1 - SQR(c_a) - SQR(c_b) - SQR(c_g) + 2 * c_a * c_b * c_g;
    if (sqr < 0) {
      strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
      return false;
    }
    double vol = a * b * c * sqrt(sqr);
    Box->Volume = vol;

    Box->transform[0][0] = a;
    Box->transform[1][0] = 0;
    Box->transform[2][0] = 0;
    Box->transform[0][1] = b * c_g;
    Box->transform[1][1] = b * s_g;
    Box->transform[2][1] = 0;
    Box->transform[0][2] = c * c_b;
    Box->transform[1][2] = c * (c_a - c_b * c_g) / s_g;
    Box->transform[2][2] = vol / (a * b * s_g);

    Box->inverse[0][0] = 1 / a;
    Box->inverse[1][0] = 0;
    Box->inverse[2][0] = 0;
    Box->inverse[0][1] = -c_g / (a * s_g);
    Box->inverse[1][1] = 1 / (b * s_g);
    Box->inverse[2][1] = 0;
    Box->inverse[0][2] = b * c * (c_g * (c_a - c_b * c_g) / (s_g * vol) -
                           c_b * s_g / vol);
    Box->inverse[1][2] = -a * c * (c_a - c_b * c_g) / (vol * s_g);
    Box->inverse[2][2] = a * b * s_g / vol;
    // tilt for triclinic axes (xy, xz, and yz) & lx, ly, and lz (lammps labels)
    Box->TriLength.x = a;
    Box->TriTilt[0] = b * c_g; // xy
    Box->TriTilt[1] = c * c_b; // xz
    sqr = SQR(b) - SQR(Box->TriTilt[0]);
    if (sqr < 0) {
      strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
      return false;
    }
    Box->TriLength.y = sqrt(sqr);
    Box->TriTilt[2] = (b * c_a - Box->TriTilt[0] * Box->TriTilt[1]) /
                        Box->TriLength.y;
    sqr = SQR(c) - SQR(Box->TriTilt[1]) - SQR(Box->TriTilt[2]);
    if (sqr < 0) {
      ErrorPrintError_old();
      //TODO coloured output
      fprintf(stderr, "Error - wrong dimensions");
      exit(1);
    }
    Box->TriLength.z = sqrt(sqr);
    // make tilt component zero if they're close to zero
    for (int i = 0; i < 3; i++) {
      if (fabs(Box->TriTilt[i]) < 0.001) {
        Box->TriTilt[i] = 0;
      }
    }
//  printf("\nLength:    %lf %lf %lf\n", Box->Length.x, Box->Length.y, Box->Length.z);
//  printf("TirLength: %lf %lf %lf\n", Box->TriLength.x, Box->TriLength.y, Box->TriLength.z);
//  printf("TriTilt:   %lf %lf %lf\n", Box->TriTilt[0], Box->TriTilt[1], Box->TriTilt[2]);
  //}}}
  } else { // orthogonal box //{{{
    Box->Volume = Box->Length.x * Box->Length.y * Box->Length.z;

    Box->transform[0][0] = Box->Length.x;
    Box->transform[1][0] = 0;
    Box->transform[2][0] = 0;
    Box->transform[0][1] = 0;
    Box->transform[1][1] = Box->Length.y;
    Box->transform[2][1] = 0;
    Box->transform[0][2] = 0;
    Box->transform[1][2] = 0;
    Box->transform[2][2] = Box->Length.z;

    Box->inverse[0][0] = 1 / Box->Length.x;
    Box->inverse[1][0] = 0;
    Box->inverse[2][0] = 0;
    Box->inverse[0][1] = 0;
    Box->inverse[1][1] = 1 / Box->Length.y;
    Box->inverse[2][1] = 0;
    Box->inverse[0][2] = 0;
    Box->inverse[1][2] = 0;
    Box->inverse[2][2] = 1 / Box->Length.z;
    // lx, ly, and lz for lammps
    Box->TriLength = Box->Length;
    // tilt for triclinic axes (xy, xz, and yz)
    Box->TriTilt[0] = 0;
    Box->TriTilt[1] = 0;
    Box->TriTilt[2] = 0;
  } //}}}
  // test print the matrices //{{{
//printf("Transformation matrix:\n");
//printf("   %lf %lf %lf\n", Box->transform[0][0], Box->transform[0][1], Box->transform[0][2]);
//printf("   %lf %lf %lf\n", Box->transform[1][0], Box->transform[1][1], Box->transform[1][2]);
//printf("   %lf %lf %lf\n", Box->transform[2][0], Box->transform[2][1], Box->transform[2][2]);
//printf("Inverse of transformation matrix:\n");
//printf("   %lf %lf %lf\n", Box->inverse[0][0], Box->inverse[0][1], Box->inverse[0][2]);
//printf("   %lf %lf %lf\n", Box->inverse[1][0], Box->inverse[1][1], Box->inverse[1][2]);
//printf("   %lf %lf %lf\n", Box->inverse[2][0], Box->inverse[2][1], Box->inverse[2][2]); //}}}
  return true;
} //}}}

// ToFractional() //{{{
void ToFractional(VECTOR *coor, BOX Box) {
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    VECTOR new = {0, 0, 0};
    new.x = Box.inverse[0][0] * coor->x +
            Box.inverse[0][1] * coor->y +
            Box.inverse[0][2] * coor->z;
    new.y = Box.inverse[1][0] * coor->x +
            Box.inverse[1][1] * coor->y +
            Box.inverse[1][2] * coor->z;
    new.z = Box.inverse[2][0] * coor->x +
            Box.inverse[2][1] * coor->y +
            Box.inverse[2][2] * coor->z;
    coor->x = new.x * Box.Length.x;
    coor->y = new.y * Box.Length.y;
    coor->z = new.z * Box.Length.z;
  }
} //}}}

// ToFractionalCoor() //{{{
void ToFractionalCoor(SYSTEM *System) {
  if (System->Box.alpha != 90 ||
      System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.Bead; i++) {
      if (System->Bead[i].InTimestep) {
        ToFractional(&System->Bead[i].Position, System->Box);
      }
    }
  }
} //}}}

// FromFractional() //{{{
VECTOR FromFractional(VECTOR coor, BOX Box) {
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    coor.x /= Box.Length.x;
    coor.y /= Box.Length.y;
    coor.z /= Box.Length.z;
    VECTOR new = {0, 0, 0};
    new.x = Box.transform[0][0] * coor.x +
            Box.transform[0][1] * coor.y +
            Box.transform[0][2] * coor.z;
    new.y = Box.transform[1][0] * coor.x +
            Box.transform[1][1] * coor.y +
            Box.transform[1][2] * coor.z;
    new.z = Box.transform[2][0] * coor.x +
            Box.transform[2][1] * coor.y +
            Box.transform[2][2] * coor.z;
    coor.x = new.x;
    coor.y = new.y;
    coor.z = new.z;
  }
  return coor;
} //}}}

// FromFractionalCoor() //{{{
void FromFractionalCoor(SYSTEM *System) {
  if (System->Box.alpha != 90 ||
      System->Box.beta != 90 ||
      System->Box.gamma != 90) {
    for (int i = 0; i < System->Count.Bead; i++) {
      BEAD *bead = &System->Bead[i];
      bead->Position = FromFractional(bead->Position, System->Box);
    }
  }
} //}}}

// InputCoor() //{{{
/**
 * Function to test whether input coordinate file is vtf or vcf and assign
 * default structure file name as either the vtf or traject.vsf
 */
bool InputCoor(bool *vtf, char *file_coor, char *file_struct) {
  int ext = 2;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  strcpy(extension[1], ".vtf");
  ext = ErrorExtension(file_coor, ext, extension);
  *vtf = false; // file_coor is vcf by default
  if (ext == -1) {
    return false; // wrong extension to file_coor
  } else if (ext == 1) {
    *vtf = true; // file_coor is vtf
  }
  // if vtf, copy to input_vsf
  if (*vtf) {
    strcpy(file_struct, file_coor);
  } else {
    strcpy(file_struct, "traject.vsf");
  }
  return true;
} //}}}

// PrintBox()  //{{{
/**
 * Function printing Counts structure.
 */
void PrintBox(BOX Box) {
  fprintf(stdout, "Box = {\n");
  fprintf(stdout, "  .Length = (%lf, %lf, %lf),\n", Box.Length.x,
                                                    Box.Length.y,
                                                    Box.Length.z);
  fprintf(stdout, "  .TriLength = (%lf, %lf, %lf),\n", Box.TriLength.x,
                                                       Box.TriLength.y,
                                                       Box.TriLength.z);
  fprintf(stdout, "  .TriTilt = (%lf, %lf, %lf),\n", Box.TriTilt[0],
                                                     Box.TriTilt[1],
                                                     Box.TriTilt[2]);
  fprintf(stdout, "  .alpha = %lf,\n", Box.alpha);
  fprintf(stdout, "  .beta  = %lf,\n", Box.beta);
  fprintf(stdout, "  .gamma = %lf,\n", Box.gamma);
  fprintf(stdout, "  .transform = (%9.5f, %9.5f, %9.5f)\n",
          Box.transform[0][0], Box.transform[0][1], Box.transform[0][2]);
  fprintf(stdout, "               (%9.5f, %9.5f, %9.5f)\n",
          Box.transform[1][0], Box.transform[1][1], Box.transform[1][2]);
  fprintf(stdout, "               (%9.5f, %9.5f, %9.5f)\n",
          Box.transform[2][0], Box.transform[2][1], Box.transform[2][2]);
  fprintf(stdout, "  .inverse = (%9.5f, %9.5f, %9.5f)\n",
          Box.inverse[0][0], Box.inverse[0][1], Box.inverse[0][2]);
  fprintf(stdout, "             (%9.5f, %9.5f, %9.5f)\n",
          Box.inverse[1][0], Box.inverse[1][1], Box.inverse[1][2]);
  fprintf(stdout, "             (%9.5f, %9.5f, %9.5f)\n",
          Box.inverse[2][0], Box.inverse[2][1], Box.inverse[2][2]);
  fprintf(stdout, "  .Volume = %lf,\n", Box.Volume);
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

// PrintBondTypes() //{{{
void PrintBondTypes(COUNTS Counts, PARAMS *bond_type) {
  for (int i = 0; i < Counts.TypesOfBonds; i++) {
    fprintf(stdout, "bond %2d: k = %lf, r_0 = %lf\n", i+1, bond_type[i].a, bond_type[i].b);
  }
  putc('\n', stdout);
} //}}}

// PrintBondTypes2() //{{{
void PrintBondTypes2(int number_of_bonds, PARAMS *bond_type) {
  for (int i = 0; i < number_of_bonds; i++) {
    fprintf(stdout, "BondType[%d] = {", i);
    fprintf(stdout, ".k = %9.5f, ", bond_type[i].a);
    fprintf(stdout, ".r_0 = %9.5f", bond_type[i].b);
    fprintf(stdout, "}\n");
  }
} //}}}

// PrintAngleTypes() //{{{
void PrintAngleTypes(COUNTS Counts, PARAMS *angle_type) {
  for (int i = 0; i < Counts.TypesOfAngles; i++) {
    fprintf(stdout, "angle %2d: k = %lf, r_0 = %lf\n", i+1, angle_type[i].a, angle_type[i].b);
  }
  putc('\n', stdout);
} //}}}

// PrintAngleTypes2() //{{{
void PrintAngleTypes2(int number_of_angles, PARAMS *angle_type) {
  for (int i = 0; i < number_of_angles; i++) {
    fprintf(stdout, "AngleType[%d] = {", i);
    fprintf(stdout, ".k = %9.5f, ", angle_type[i].a);
    fprintf(stdout, ".theta_0 = %9.5f", angle_type[i].b);
    fprintf(stdout, "}\n");
  }
} //}}}

// PrintDihedralTypes2() //{{{
void PrintDihedralTypes2(int number_of_dihedrals, PARAMS *dihedral_type) {
  for (int i = 0; i < number_of_dihedrals; i++) {
    fprintf(stdout, "DihedralType[%d] = {", i);
    fprintf(stdout, ".k = %9.5f, ", dihedral_type[i].a);
    fprintf(stdout, ".theta_0 = %9.5f", dihedral_type[i].b);
    fprintf(stdout, "}\n");
  }
} //}}}

int FindBeadType(char name[], SYSTEM System) { //{{{
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
int FindMoleculeType(char name[], SYSTEM System) { //{{{
  int type;
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    if (strcmp(name, System.MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return(-1);
} //}}}

// fill the some System arrays //{{{
void FillSystemNonessentials(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  FillBeadTypeIndex(System);
  FillMoleculeTypeIndex(System);
  System->Index_mol = realloc(System->Index_mol, sizeof *System->Index_mol *
                              (Count->HighestResid + 1));
  for (int i = 0; i <= Count->HighestResid; i++) {
    System->Index_mol[i] = -1;
  }
  for (int i = 0; i < Count->Molecule; i++) {
    System->Index_mol[System->Molecule[i].Index] = i;
  }
  System->BeadCoor = realloc(System->BeadCoor,
                             sizeof *System->BeadCoor * Count->Bead);
  if (Count->Bonded > 0) {
    System->BondedCoor = realloc(System->BondedCoor,
                                 sizeof *System->BondedCoor * Count->Bonded);
  }
  if (Count->Unbonded > 0) {
    System->UnbondedCoor = realloc(System->UnbondedCoor,
                                   sizeof *System->UnbondedCoor *
                                   Count->Unbonded);
  }
} //}}}
// molecule type's nBTypes and BType array //{{{
void FillMoleculeTypeBType(MOLECULETYPE *MoleculeType) {
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
      MoleculeType->BType = realloc(MoleculeType->BType,
                                    sizeof *MoleculeType->BType *
                                    MoleculeType->nBTypes);
      MoleculeType->BType[type] = MoleculeType->Bead[j];
    }
  }
} //}}}
// molecule type's Charge and Mass //{{{
void FillMoleculeTypeChargeMass(MOLECULETYPE *MoleculeType,
                                BEADTYPE BeadType[]) {
  MoleculeType->Mass = 0;
  MoleculeType->Charge = 0;
  for (int j = 0; j < MoleculeType->nBeads; j++) {
    int btype = MoleculeType->Bead[j];
    // charge
    if (BeadType[btype].Charge != CHARGE) {
      MoleculeType->Charge += BeadType[btype].Charge;
    } else {
      MoleculeType->Charge = CHARGE;
      break;
    }
    // mass
    if (BeadType[btype].Mass != MASS) {
      MoleculeType->Mass += BeadType[btype].Mass;
    } else {
      MoleculeType->Mass = MASS;
      break;
    }
  }
} //}}}
// all bead types' Index array //{{{
void FillBeadTypeIndex(SYSTEM *System) {
  COUNT *Count = &System->Count;
  // allocate memory for Index arrays
  for (int i = 0; i < Count->BeadType; i++) {
    BEADTYPE *bt_i = &System->BeadType[i];
    bt_i->Index = malloc(sizeof *bt_i->Index * bt_i->Number);
  }
  // fill the Index arrays
  int *count_id = calloc(Count->BeadType, sizeof *count_id);
  for (int i = 0; i < Count->Bead; i++) {
    int type = System->Bead[i].Type;
    System->BeadType[type].Index[count_id[type]] = i;
    count_id[type]++;
  }
  free(count_id);
} //}}}
// all molecule types' Index array //{{{
void FillMoleculeTypeIndex(SYSTEM *System) {
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
 //}}}

// Distance() //{{{
/**
 * Function calculating distance vector between two beads. It removes
 * periodic boundary conditions and returns x, y, and z distances in the
 * range <0, BoxLength/2).
 */
VECTOR Distance(VECTOR id1, VECTOR id2, VECTOR BoxLength) {
  // distance vector
  VECTOR rij;
  rij.x = id1.x - id2.x;
  rij.y = id1.y - id2.y;
  rij.z = id1.z - id2.z;
  // remove periodic boundary conditions in x-direction
  while (rij.x >= (BoxLength.x/2))
    rij.x = rij.x - BoxLength.x;
  while (rij.x < -(BoxLength.x/2))
    rij.x = rij.x + BoxLength.x;
  // in y-direction
  while (rij.y >= (BoxLength.y/2))
    rij.y = rij.y - BoxLength.y;
  while (rij.y < -(BoxLength.y/2))
    rij.y = rij.y + BoxLength.y;
  // in z-direction
  while (rij.z >= (BoxLength.z/2))
    rij.z = rij.z - BoxLength.z;
  while (rij.z < -(BoxLength.z/2))
    rij.z = rij.z + BoxLength.z;
  return rij;
} //}}}

/* TODO what about if the molecule has bonds, but is in more 'pieces'; should
 *      it really just try joining it 1000 times when we know it's impossible?
 *      There should be some check and then the separate pieces should be
 *      connected (with a warning issued).
 */
// RemovePBCMolecules() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them. The function requires orthogonal box, i.e.,
 * for triclinic box, the supplied coordinates must first be transformed.
 */
void RemovePBCMolecules(SYSTEM *System) {
  // go through all molecules
  for (int i = 0; i < (*System).Count.Molecule; i++) {
    int type = (*System).Molecule[i].Type;
    // do nothing if the molecule has no bonds
    if ((*System).MoleculeType[type].nBonds == 0) {
      continue;
    }
    for (int j = 0; j < (*System).MoleculeType[type].nBeads; j++) {
      int id = (*System).Molecule[i].Bead[j];
      (*System).Bead[id].Flag = false; // no beads moved yet
    }
    // first bead in the first bond is considered moved
    for (int j = 0; j < (*System).MoleculeType[type].nBonds; j++) {
      int id1 = (*System).MoleculeType[type].Bond[j][0];
      id1 = (*System).Molecule[i].Bead[id1];
      int id2 = (*System).MoleculeType[type].Bond[j][1];
          id2 = (*System).Molecule[i].Bead[id2];
      if ((*System).Bead[id1].InTimestep && (*System).Bead[id2].InTimestep) {
        (*System).Bead[id1].Flag = true;
        break;
      }
    }
    bool done = false;
    int test = 0; // if too many loops, just leave the loop with error
    while (!done && test < 1000) {
      for (int j = 0; j < (*System).MoleculeType[type].nBonds; j++) {
        int id1 = (*System).MoleculeType[type].Bond[j][0];
        id1 = (*System).Molecule[i].Bead[id1];
        int id2 = (*System).MoleculeType[type].Bond[j][1];
        id2 = (*System).Molecule[i].Bead[id2];
        if ((*System).Bead[id1].InTimestep && (*System).Bead[id2].InTimestep) {
          // move id1, if id2 is moved already
          if (!(*System).Bead[id1].Flag && (*System).Bead[id2].Flag) {
            VECTOR dist = Distance((*System).Bead[id2].Position,
                                   (*System).Bead[id1].Position,
                                   (*System).Box.Length);
            (*System).Bead[id1].Position.x = (*System).Bead[id2].Position.x -
                                             dist.x;
            (*System).Bead[id1].Position.y = (*System).Bead[id2].Position.y -
                                             dist.y;
            (*System).Bead[id1].Position.z = (*System).Bead[id2].Position.z -
                                             dist.z;
            (*System).Bead[id1].Flag = true;
          // move id2, if id1 was moved already
          } else if ((*System).Bead[id1].Flag && !(*System).Bead[id2].Flag) {
            VECTOR dist = Distance((*System).Bead[id1].Position, (*System).Bead[id2].Position,
                                   (*System).Box.Length);
            (*System).Bead[id2].Position.x = (*System).Bead[id1].Position.x -
                                             dist.x;
            (*System).Bead[id2].Position.y = (*System).Bead[id1].Position.y -
                                             dist.y;
            (*System).Bead[id2].Position.z = (*System).Bead[id1].Position.z -
                                             dist.z;
            (*System).Bead[id2].Flag = true;
          }
        }
      }

      // break while loop if all beads have moved
      done = true;
      for (int j = 1; j < (*System).MoleculeType[type].nBeads; j++) {
        int id = (*System).Molecule[i].Bead[j];
        if ((*System).Bead[id].InTimestep && !(*System).Bead[id].Flag) {
          done = false;
          break;
        }
      }
      test++;
    }
    if (test == 1000) {
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: unable to 'join' molecule");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%s", (*System).MoleculeType[type].Name);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " (resid ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%d", i+1);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " )\n");
      fprintf(stderr, "           Maybe not all beads are connected?\n");
      ColourReset(STDERR_FILENO);
    }

    // put molecule's geometric centre into the simulation box //{{{
    VECTOR cog = GeomCentre((*System).MoleculeType[type].nBeads,
                            (*System).Molecule[i].Bead,
                            (*System).Bead);
    // by how many BoxLength's should cog be moved?
    // for distant molecules - it shouldn't happen, but better safe than sorry
    INTVECTOR move;
    move.x = cog.x / (*System).Box.Length.x;
    move.y = cog.y / (*System).Box.Length.y;
    move.z = cog.z / (*System).Box.Length.z;
    if (cog.x < 0) {
      move.x--;
    }
    if (cog.y < 0) {
      move.y--;
    }
    if (cog.z < 0) {
      move.z--;
    }
    for (int j = 0; j < (*System).MoleculeType[type].nBeads; j++) {
      int bead = (*System).Molecule[i].Bead[j];
      (*System).Bead[bead].Position.x -= move.x * (*System).Box.Length.x;
      (*System).Bead[bead].Position.y -= move.y * (*System).Box.Length.y;
      (*System).Bead[bead].Position.z -= move.z * (*System).Box.Length.z;
    } //}}}
  }
} //}}}

// RestorePBC() //{{{
/**
 * Function to restore removed periodic boundary conditions. Used also in case
 * of cell linked lists, because they need coordinates <0, BoxLength>.
 */
void RestorePBC(int number_of_beads, BOX Box, BEAD **Bead) {
  for (int i = 0; i < number_of_beads; i++) {
    // TODO: whiles - really? There's remainder() in math.h, I think, (or
    //       fmod() or some such)
    // x direction
    while ((*Bead)[i].Position.x >= Box.Length.x) {
      (*Bead)[i].Position.x -= Box.Length.x;
    }
    while ((*Bead)[i].Position.x < 0) {
      (*Bead)[i].Position.x += Box.Length.x;
    }
    // y direction
    while ((*Bead)[i].Position.y >= Box.Length.y) {
      (*Bead)[i].Position.y -= Box.Length.y;
    }
    while ((*Bead)[i].Position.y < 0) {
      (*Bead)[i].Position.y += Box.Length.y;
    }
    // z direction
    while ((*Bead)[i].Position.z >= Box.Length.z) {
      (*Bead)[i].Position.z -= Box.Length.z;
    }
    while ((*Bead)[i].Position.z < 0) {
      (*Bead)[i].Position.z += Box.Length.z;
    }
  }
} //}}}

// CentreOfMass() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
// TODO: what if bead's mass is undefined?
VECTOR CentreOfMass(int n, int *list, BEAD *Bead, BEADTYPE *BeadType) {
  VECTOR com = {0, 0, 0};
  double mass = 0;
  for (int i = 0; i < n; i++) {
    int id = list[i];
    int btype = Bead[id].Type;
    com.x += Bead[id].Position.x * BeadType[btype].Mass;
    com.y += Bead[id].Position.y * BeadType[btype].Mass;
    com.z += Bead[id].Position.z * BeadType[btype].Mass;
    mass += BeadType[btype].Mass;
  }
  com.x /= mass;
  com.y /= mass;
  com.z /= mass;
  return com;
} //}}}

// GeomCentre() //{{{
/**
 * Function to calculate centre of mass for a given list of beads.
 */
VECTOR GeomCentre(int n, int *list, BEAD *Bead) {
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

// Gyration() //{{{
/**
 * Function to calculate the principle moments of the gyration tensor.
 */
VECTOR Gyration(int n, int *list, COUNTS Counts,
                BEADTYPE *BeadType, BEAD **Bead) {
  // gyration tensor (3x3 array)
  // use long double to ensure precision -- previous problem with truncation in short chains
  struct Tensor {
    LONGVECTOR x, y, z;
  } GyrationTensor = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

  // TODO: shouldn't there be centre of mass?!?
  VECTOR com = GeomCentre(n, list, *Bead);

  // move centre of mass to [0,0,0] //{{{
  for (int i = 0; i < n; i++) {
    (*Bead)[list[i]].Position.x -= com.x;
    (*Bead)[list[i]].Position.y -= com.y;
    (*Bead)[list[i]].Position.z -= com.z;
  } //}}}

  // calculate gyration tensor //{{{
  for (int i = 0; i < n; i++) {
    int id = list[i];
    GyrationTensor.x.x += (*Bead)[id].Position.x * (*Bead)[id].Position.x;
    GyrationTensor.x.y += (*Bead)[id].Position.x * (*Bead)[id].Position.y;
    GyrationTensor.x.z += (*Bead)[id].Position.x * (*Bead)[id].Position.z;
    GyrationTensor.y.y += (*Bead)[id].Position.y * (*Bead)[id].Position.y;
    GyrationTensor.y.z += (*Bead)[id].Position.y * (*Bead)[id].Position.z;
    GyrationTensor.z.z += (*Bead)[id].Position.z * (*Bead)[id].Position.z;
  }
  GyrationTensor.x.x /= n;
  GyrationTensor.x.y /= n;
  GyrationTensor.x.z /= n;
  GyrationTensor.y.y /= n;
  GyrationTensor.y.z /= n;
  GyrationTensor.z.z /= n; //}}}

  // char polynomial: a_cube * x^3 + b_cube * x^2 + c_cube * x + d_cube = 0 //{{{
  long double a_cube = -1;
  long double b_cube = GyrationTensor.x.x + GyrationTensor.y.y + GyrationTensor.z.z;
  long double c_cube = - GyrationTensor.x.x * GyrationTensor.y.y
                  - GyrationTensor.x.x * GyrationTensor.z.z
                  - GyrationTensor.y.y * GyrationTensor.z.z
                  + SQR(GyrationTensor.y.z)
                  + SQR(GyrationTensor.x.y)
                  + SQR(GyrationTensor.x.z);
  long double d_cube = + GyrationTensor.x.x * GyrationTensor.y.y * GyrationTensor.z.z
                  + 2 * GyrationTensor.x.y * GyrationTensor.y.z * GyrationTensor.x.z
                  - SQR(GyrationTensor.x.z) * GyrationTensor.y.y
                  - SQR(GyrationTensor.x.y) * GyrationTensor.z.z
                  - SQR(GyrationTensor.y.z) * GyrationTensor.x.x; //}}}

  // first root: either 0 or Newton's iterative method to get it //{{{
  long double root0 = 0;
  if (fabsl(d_cube) > 0.0000000001L) {
    // derivative of char. polynomial: a_deriv * x^2 + b_deriv * x + c_deriv
    long double a_deriv = 3 * a_cube;
    long double b_deriv = 2 * b_cube;
    long double c_deriv = c_cube;

    long double root1 = 1;

    while (fabsl(root0-root1) > 0.0000000001L) {
      long double f_root0 = (a_cube * CUBE(root0) + b_cube * SQR(root0) + c_cube * root0 + d_cube);
      long double f_deriv_root0 = (a_deriv * SQR(root0) + b_deriv * root0 + c_deriv);
      root1 = root0 - f_root0 / f_deriv_root0;

      // swap root0 and root1 for the next iteration
      long double tmp = root0;
      root0 = root1;
      root1 = tmp;
    }
  } //}}}

  // determine parameters of quadratic equation a_quad * x^2 + b_quad * x + c_quad = 0 //{{{
  // derived by division: (x^3 + (b_cube/a_cube) * x^2 + (c_cube/a_cube) * x + (d_cube/a_cube)):(x - root0)
  long double a_quad = 1;
  long double b_quad = b_cube / a_cube + root0;
  long double c_quad = SQR(root0) + b_cube / a_cube * root0 + c_cube/a_cube; //}}}
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

// LinkedList() //{{{
void LinkedList(VECTOR BoxLength, COUNTS Counts, BEAD *Bead,
                int **Head, int **Link, double cell_size, INTVECTOR *n_cells,
                int *Dcx, int *Dcy, int *Dcz) {

  (*n_cells).x = ceil(BoxLength.x/cell_size),
  (*n_cells).y = ceil(BoxLength.y/cell_size),
  (*n_cells).z = ceil(BoxLength.z/cell_size);

  // allocate arrays
  *Head = malloc(sizeof **Head * (*n_cells).x * (*n_cells).y * (*n_cells).z);
  *Link = malloc(sizeof **Link * Counts.BeadsCoor);
  for (int i = 0; i < ((*n_cells).x*(*n_cells).y*(*n_cells).z); i++) {
    (*Head)[i] = -1;
  }

  // sort beads into cells //{{{
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int cell = (int)(Bead[i].Position.x / cell_size)
             + (int)(Bead[i].Position.y / cell_size) * (*n_cells).x
             + (int)(Bead[i].Position.z / cell_size) * (*n_cells).x * (*n_cells).y;
    (*Link)[i] = (*Head)[cell];
    (*Head)[cell] = i;
  } //}}}

  // coordinates of adjoining cells //{{{
  int x[14] = {0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
  int y[14] = {0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
  int z[14] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  for (int i = 0; i < 14; i++) {
    Dcx[i] = x[i];
    Dcy[i] = y[i];
    Dcz[i] = z[i];
  } //}}}
} //}}}

// SortBonds() //{{{
/**
 * Function to sort a bond array. As each bond contains a 2-member array, first
 * check that first bead id is lower than the second. Then sort the bonds
 * according to the index of the first bead in each bond.
 */
void SortBonds(int (*bond)[3], int number_of_bonds) {
  // first, check order in every bond
  for (int j = 0; j < number_of_bonds; j++) {
    if (bond[j][0] > bond[j][1]) {
      SwapInt(&bond[j][0], &bond[j][1]);
    }
  }
  // second, bubble sort bonds
  for (int j = 0; j < (number_of_bonds-1); j++) {
    bool swap = false;
    for (int k = 0; k < (number_of_bonds-j-1); k++) {
      if ((bond[k][0] > bond[k+1][0]) || // swap if first beads are in wrong order
          (bond[k][0] == bond[k+1][0] && // or if they're the same, but second ones are in wrong order
          bond[k][1] > bond[k+1][1])) {
        SwapInt(&bond[k][0], &bond[k+1][0]);
        SwapInt(&bond[k][1], &bond[k+1][1]);
        SwapInt(&bond[k][2], &bond[k+1][2]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// SortAngles() //{{{
/**
 * Function to sort an angle array. As each angle contains a 3-member array
 * with the middle number being the 'centre' of the angle (i.e., it must remain
 * in the middle). Sort it so that the first index is lower than the third one
 * and then ascendingly according to the first indices.
 */
void SortAngles(int (*angle)[4], int length) {
  // first, check order of the 1st and 3rd id in every angle
  for (int j = 0; j < length; j++) {
    if (angle[j][0] > angle[j][2]) {
      SwapInt(&angle[j][0], &angle[j][2]);
    }
  }
  // second, bubble sort angles
  for (int j = 0; j < (length-1); j++) {
    bool swap = false;
    for (int k = 0; k < (length-j-1); k++) {
      if ((angle[k][0] > angle[k+1][0]) || // swap if first beads are in wrong order
          (angle[k][0] == angle[k+1][0] &&
           angle[k][1] > angle[k+1][1]) || // or if they're the same, but 2nd ones are in wrong order
          (angle[k][0] == angle[k+1][0] &&
           angle[k][1] == angle[k+1][1] &&
           angle[k][2] > angle[k+1][2])) { // same for 3rd
        SwapInt(&angle[k][0], &angle[k+1][0]);
        SwapInt(&angle[k][1], &angle[k+1][1]);
        SwapInt(&angle[k][2], &angle[k+1][2]);
        SwapInt(&angle[k][3], &angle[k+1][3]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// SortDihedrals() //{{{
/**
 * Function to sort an angle array. As each angle contains a 3-member array
 * with the middle number being the 'centre' of the angle (i.e., it must remain
 * in the middle). Sort it so that the first index is lower than the third one
 * and then ascendingly according to the first indices.
 */
void SortDihedrals(int (*dihedral)[5], int length) {
  // first, check order of the 1st and 4th id in every dihedral
  for (int j = 0; j < length; j++) {
    if (dihedral[j][0] > dihedral[j][3]) {
      SwapInt(&dihedral[j][0], &dihedral[j][3]);
      SwapInt(&dihedral[j][1], &dihedral[j][2]);
    }
  }
  // second, bubble sort dihedrals
  for (int j = 0; j < (length-1); j++) {
    bool swap = false;
    for (int k = 0; k < (length-j-1); k++) {
      if ((dihedral[k][0] > dihedral[k+1][0]) || // swap if first beads are in wrong order
          (dihedral[k][0] == dihedral[k+1][0] &&
           dihedral[k][1] > dihedral[k+1][1]) || // or if they're the same, but 2nd ones are in wrong order
          (dihedral[k][0] == dihedral[k+1][0] &&
           dihedral[k][1] == dihedral[k+1][1] &&
           dihedral[k][2] > dihedral[k+1][2]) || // same for 3rd...
          (dihedral[k][0] == dihedral[k+1][0] &&
           dihedral[k][1] == dihedral[k+1][1] &&
           dihedral[k][2] == dihedral[k+1][2] &&
           dihedral[k][3] > dihedral[k+1][3])) { // ...and for 4th.
        SwapInt(&dihedral[k][0], &dihedral[k+1][0]);
        SwapInt(&dihedral[k][1], &dihedral[k+1][1]);
        SwapInt(&dihedral[k][2], &dihedral[k+1][2]);
        SwapInt(&dihedral[k][3], &dihedral[k+1][3]);
        SwapInt(&dihedral[k][4], &dihedral[k+1][4]);
        swap = true;
      }
    }
    // if no swap was made, the array is sorted
    if (!swap) {
      break;
    }
  }
} //}}}

// CopyBeadType() //{{{
/**
 * Function to copy BEADTYPE structure into a new. Memory for the new one will be
 * reallocated in this function. Memory management for the output structure:
 * sufficient memory is already allocated (mode=0), the structure needs freeing
 * and allocating (mode=1), no memory wasn't yet allocated at all (mode=2), or
 * it needs reallocating (mode=3).
 */
void CopyBeadType(int number_of_types, BEADTYPE **bt_out,
                  BEADTYPE *bt_in, int mode) {
  // bt_out memory management //{{{
  switch (mode) {
    case 0:
      break;
    case 1:
      free(*bt_out);
      *bt_out = malloc(sizeof (BEADTYPE) * number_of_types);
      break;
    case 2:
      *bt_out = malloc(sizeof (BEADTYPE) * number_of_types);
      break;
    case 3:
      *bt_out = realloc(*bt_out, sizeof (BEADTYPE) * number_of_types);
      break;
    default:
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "CopyBeadType()");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " function requires mode=0, 1, 2, or 3\n");
      ColourReset(STDERR_FILENO);
      exit(1);
  } //}}}
  for (int i = 0; i < number_of_types; i++) {
    (*bt_out)[i] = bt_in[i];
  }
} //}}}

// CopySystem() //{{{
/*
 * Assumes unallocated SYSTEM.
 */
SYSTEM CopySystem(SYSTEM S_in) {
  SYSTEM S_out;
  InitSystem(&S_out);
  S_out.Box = S_in.Box;
  if (S_in.Count.Bead > 0) {
    S_out.Count = S_in.Count;
    // BeadType //{{{
    if (S_out.Count.BeadType > 0) {
      S_out.BeadType = realloc(S_out.BeadType,
                               sizeof (BEADTYPE) * S_out.Count.BeadType);
      for (int i = 0; i < S_out.Count.BeadType; i++) {
        S_out.BeadType[i] = S_in.BeadType[i];
        S_out.BeadType[i].Index = malloc(sizeof *S_out.BeadType[i].Index *
                                         S_out.BeadType[i].Number);
        memcpy(S_out.BeadType[i].Index, S_in.BeadType[i].Index,
               sizeof *S_in.BeadType[i].Index * S_in.BeadType[i].Number);
      }
    } else {
      strcpy(ERROR_MSG, "no bead types to copy; should never happen!");
      PrintWarning();
      return S_out;
    } //}}}
    // Bead & BeadCoor //{{{
    S_out.Bead = realloc(S_out.Bead,
                         sizeof (BEAD) * S_out.Count.Bead);
    S_out.BeadCoor = realloc(S_out.BeadCoor,
                              sizeof *S_out.BeadCoor * S_out.Count.Bead);
    memcpy(S_out.Bead, S_in.Bead, sizeof (BEAD) * S_in.Count.Bead);
    memcpy(S_out.BeadCoor, S_in.BeadCoor,
           sizeof *S_in.BeadCoor * S_in.Count.Bead); //}}}
    // Bonded & BondedCoor //{{{
    if (S_out.Count.Bonded > 0) {
      S_out.Bonded = realloc(S_out.Bonded,
                             sizeof *S_out.Bonded * S_out.Count.Bonded);
      S_out.BondedCoor = realloc(S_out.BondedCoor,
                                 sizeof *S_out.BondedCoor * S_out.Count.Bonded);
      memcpy(S_out.Bonded, S_in.Bonded,
             sizeof *S_in.Bonded * S_in.Count.Bonded);
      memcpy(S_out.BondedCoor, S_in.BondedCoor,
             sizeof *S_in.BondedCoor * S_in.Count.Bonded);
    } //}}}
    // Unbonded & UnbondedCoor //{{{
    if (S_out.Count.Unbonded > 0) {
      S_out.Unbonded = realloc(S_out.Unbonded, sizeof *S_out.Unbonded *
                               S_out.Count.Unbonded);
      S_out.UnbondedCoor = realloc(S_out.UnbondedCoor,
                                   sizeof *S_out.UnbondedCoor *
                                   S_out.Count.Unbonded);
      memcpy(S_out.Unbonded, S_in.Unbonded,
             sizeof *S_in.Unbonded * S_in.Count.Unbonded);
      memcpy(S_out.UnbondedCoor, S_in.UnbondedCoor,
             sizeof *S_in.UnbondedCoor * S_in.Count.Unbonded);
    } //}}}
    // MoleculeType //{{{
    if (S_out.Count.MoleculeType > 0) {
      S_out.MoleculeType = realloc(S_out.MoleculeType,
                                   sizeof (MOLECULETYPE) * S_out.Count.MoleculeType);
      for (int i = 0; i < S_out.Count.MoleculeType; i++) {
        S_out.MoleculeType[i] = CopyMoleculeType(S_in.MoleculeType[i]);
      }
    } //}}}
    // Molecule & Index_mol //{{{
    if (S_out.Count.Molecule > 0) {
      S_out.Molecule = realloc(S_out.Molecule,
                               sizeof (MOLECULE) * S_out.Count.Molecule);
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
      S_out.BondType = realloc(S_out.BondType,
                               sizeof (PARAMS) * S_out.Count.BondType);
      memcpy(S_out.BondType, S_in.BondType,
             sizeof (PARAMS) * S_in.Count.BondType);
    } //}}}
    // AngleType //{{{
    if (S_out.Count.AngleType > 0) {
      S_out.AngleType = realloc(S_out.AngleType,
                                sizeof (PARAMS) * S_out.Count.AngleType);
      memcpy(S_out.AngleType, S_in.AngleType,
             sizeof (PARAMS) * S_in.Count.AngleType);
    } //}}}
    // DihedralType //{{{
    if (S_out.Count.DihedralType > 0) {
      S_out.DihedralType = realloc(S_out.DihedralType,
                                   sizeof (PARAMS) * S_out.Count.DihedralType);
      memcpy(S_out.DihedralType, S_in.DihedralType,
             sizeof (PARAMS) * S_in.Count.DihedralType);
    } //}}}
  }
  return S_out;
} //}}}
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
    strcpy(ERROR_MSG, "molecule type without beads");
    PrintWarning();
    fprintf(stderr, "%sMolecule %s%s%s\n\n", ErrCyan(),
            ErrYellow(), mt_new.Name, ErrColourReset());
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
  return mt_new;
} //}}}
void PruneSystem(SYSTEM *System) { //{{{
  SYSTEM S_old = CopySystem(*System);
  FreeSystem(System);
  InitSystem(System);
  COUNT *Count = &System->Count,
        *Count_old = &S_old.Count;
  // copy bead counts
  System->Box = S_old.Box;
  Count->Bead         = Count_old->BeadCoor;
  Count->BeadCoor     = Count_old->BeadCoor;
  Count->Bonded       = Count_old->BondedCoor;
  Count->BondedCoor   = Count_old->BondedCoor;
  Count->Unbonded     = Count_old->UnbondedCoor;
  Count->UnbondedCoor = Count_old->UnbondedCoor;
  // allocate memory for bead count arrays
  System->Bead = realloc(System->Bead, sizeof (BEAD) * Count->Bead);
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
  int count_unbonded = 0, count_bonded = 0, count_all = 0,
      *connect = calloc(Count_old->Bead, sizeof *connect);
  for (int i = 0; i < Count_old->Bead; i++) {
    if (S_old.Bead[i].InTimestep) {
      System->Bead[count_all] = S_old.Bead[i];
      System->BeadCoor[count_all] = count_all;
      connect[i] = count_all;
      if (System->Bead[count_all].Molecule == -1) {
        System->Unbonded[count_unbonded] = count_all;
        System->UnbondedCoor[count_unbonded] = count_all;
        count_unbonded++;
      } else {
        System->Bonded[count_bonded] = count_all;
        System->BondedCoor[count_bonded] = count_all;
        count_bonded++;
      }
      // create new bead type if it doesn't exist yet in the pruned system
      bool new = true;
      int old_type = S_old.Bead[count_all].Type;
      for (int j = 0; j < Count->BeadType; j++) {
        int new_type = FindBeadType(S_old.BeadType[old_type].Name, *System);
        if (new_type != -1) {
          System->Bead[count_all].Type = new_type;
          System->BeadType[new_type].Number++;
          new = false;
          break;
        }
      }
      if (new) {
        int type = Count->BeadType;
        NewBeadType(&System->BeadType, &Count->BeadType,
                    S_old.BeadType[old_type].Name,
                    S_old.BeadType[old_type].Charge,
                    S_old.BeadType[old_type].Mass,
                    S_old.BeadType[old_type].Radius);
        System->BeadType[type].Number = 1;
        System->Bead[count_all].Type = type;
      }
      count_all++;
    }
  } //}}}
  // copy Molecule array & create a new MoleculeType array //{{{
  // TODO: Molecule[].Index must be adapted
  for (int i = 0; i < Count_old->Molecule; i++) {
    MOLECULE *mol_old = &S_old.Molecule[i];
    int old_type = mol_old->Type;
    MOLECULETYPE *mt_old = &S_old.MoleculeType[old_type];
    // find if the molecule is to be in the pruned system
    int c_bead = 0;
    for (int j = 0; j < mt_old->nBeads; j++) {
      int id = mol_old->Bead[j];
      if (S_old.Bead[id].InTimestep) {
        c_bead++;
      }
    }
    if (c_bead > 0) { // should the molecule be in the pruned system?
      int new_id = Count->Molecule;
      Count->Molecule++;
      System->Molecule = realloc(System->Molecule,
                                 sizeof (MOLECULE) * Count->Molecule);
      MOLECULE *mol_new = &System->Molecule[new_id];
      *mol_new = *mol_old;
      mol_new->Bead = calloc(c_bead, sizeof *mol_new->Bead);
      c_bead = 0;
      for (int j = 0; j < mt_old->nBeads; j++) {
        int id = mol_old->Bead[j];
        if (S_old.Bead[id].InTimestep) {
          mol_new->Bead[c_bead] = connect[id];
          System->Bead[connect[id]].Molecule = new_id;
          c_bead++;
        }
      }
      // find if the molecule already exists in the pruned system
      bool new = true;
      for (int j = 0; j < Count->BeadType; j++) {
        int new_type = FindMoleculeType(mt_old->Name,
                                        (*System));
        if (new_type != -1) {
          mol_new->Type = new_type;
          System->MoleculeType[new_type].Number++;
          new = false;
          break;
        }
      }
      if (new) { // create new molecule type if it doesn't exist yet
        int new_type = Count->MoleculeType,
            c_bond = 0, c_angle = 0, c_dihedral = 0;
        // count bonds in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nBonds; j++) {
          int *id = BondIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep) {
            c_bond++;
          }
        } //}}}
        // count angles in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nAngles; j++) {
          int *id = AngleIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            c_angle++;
          }
        } //}}}
        // count dihedrals in the pruned molecule type //{{{
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            c_dihedral++;
          }
        } //}}}
        NewMolType(&System->MoleculeType, &Count->MoleculeType,
                   mt_old->Name, c_bead, c_bond, c_angle, c_dihedral);
        System->Molecule[new_id].Type = new_type;
        MOLECULETYPE *mt_new = &System->MoleculeType[new_type];
        // copy beads to the new molecule type //{{{
        int *id_old_to_new = calloc(mt_old->nBeads, sizeof *id_old_to_new);
        for (int j = 0; j < mt_old->nBeads; j++) {
          id_old_to_new[j] = -1;
        }
        c_bead = 0;
        for (int j = 0; j < mt_old->nBeads; j++) {
          int id = mol_old->Bead[j],
              btype = mt_old->Bead[j];
          if (S_old.Bead[id].InTimestep) {
            mt_new->Bead[c_bead] = FindBeadType(S_old.BeadType[btype].Name,
                                                *System);
            id_old_to_new[j] = c_bead;
            c_bead++;
          }
        } //}}}
        // copy bonds to the new molecule type //{{{
        c_bond = 0;
        for (int j = 0; j < mt_old->nBonds; j++) {
          int *id = BondIndices(S_old, i, j),
              last = mt_old->Bond[j][2];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep) {
            mt_new->Bond[c_bond][0] = id_old_to_new[id[2]];
            mt_new->Bond[c_bond][1] = id_old_to_new[id[3]];
            mt_new->Bond[c_bond][2] = last;
            c_bond++;
          }
        } //}}}
        // copy angles to the new molecule type //{{{
        c_angle = 0;
        for (int j = 0; j < mt_old->nAngles; j++) {
          int *id = AngleIndices(S_old, i, j),
              last = mt_old->Angle[j][3];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            mt_new->Angle[c_angle][0] = id_old_to_new[id[3]];
            mt_new->Angle[c_angle][1] = id_old_to_new[id[4]];
            mt_new->Angle[c_angle][2] = id_old_to_new[id[5]];
            mt_new->Angle[c_angle][3] = last;
            c_angle++;
          }
        } //}}}
        // copy bonds to the new molecule type //{{{
        c_bond = 0;
        for (int j = 0; j < mt_old->nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j),
              last = mt_old->Dihedral[j][4];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            mt_new->Dihedral[c_dihedral][0] = id_old_to_new[id[4]];
            mt_new->Dihedral[c_dihedral][1] = id_old_to_new[id[5]];
            mt_new->Dihedral[c_dihedral][2] = id_old_to_new[id[6]];
            mt_new->Dihedral[c_dihedral][3] = id_old_to_new[id[7]];
            mt_new->Dihedral[c_dihedral][4] = last;
            c_dihedral++;
          }
        } //}}}
        free(id_old_to_new);
      }
    }
  }
  // TODO: FillSystemNonessentials()?
  for (int i = 0; i < Count->MoleculeType; i++) {
    FillMoleculeTypeBType(&System->MoleculeType[i]);
    FillMoleculeTypeChargeMass(&System->MoleculeType[i], System->BeadType);
  }
  //}}}
  // fill BeadType[].Index //{{{
  for (int i = 0; i < Count->BeadType; i++) {
    System->BeadType[i].Index = malloc(sizeof *System->BeadType[i].Index *
                                       System->BeadType[i].Number);
  }
  int *count = calloc(Count->BeadType, sizeof *count);
  for (int i = 0; i < Count->Bead; i++) {
    int type = System->Bead[i].Type;
    System->BeadType[type].Index[count[type]] = i;
    count[type]++;
  }
  free(count); //}}}
  // fill MoleculeType[].Index //{{{
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    mt_i->Index = malloc(sizeof *mt_i->Index * mt_i->Number);
  }
  count = calloc(Count->MoleculeType, sizeof *count);
  for (int i = 0; i < Count->Molecule; i++) {
    int type = System->Molecule[i].Type;
    System->MoleculeType[type].Index[count[type]] = i;
    count[type]++;
  }
  free(count); //}}}
  FreeSystem(&S_old);
  free(connect);
} //}}}
// ConcatenateSystems() //{{{
/*
 * Assumes S_out needs reallocating memory to accommodate S_in.
 */
// TODO: some warning about S_in being empty
// TODO: bond/angle/dihedral types - in Bead[].Bond/Angle/Dihedral, the type
//       must change if S_out and S_in contain types (akin to molecule types)
void ConcatenateSystems(SYSTEM *S_out, SYSTEM S_in, BOX Box) {
  COUNT Count_old = S_out->Count; // copy the original COUNT
  COUNT *Count_out = &S_out->Count,
        *Count_in = &S_in.Count;
  S_out->Box = Box;
  // BeadType //{{{
  if (Count_in->BeadType > 0) {
    Count_out->BeadType += Count_in->BeadType;
    S_out->BeadType = realloc(S_out->BeadType, sizeof (BEADTYPE) *
                                Count_out->BeadType);
//  memcpy(S_out->BeadType + Count_old.BeadType, S_in.BeadType,
//         sizeof (BEADTYPE) * Count_in->BeadType);
    for (int i = 0; i < Count_in->BeadType; i++) {
      int new = i + Count_old.BeadType;
      S_out->BeadType[new] = S_in.BeadType[i];
      S_out->BeadType[new].Index =
        malloc(sizeof *S_out->BeadType[new].Index *
               S_out->BeadType[new].Number);
      for (int j = 0; j < S_out->BeadType[new].Number; j++) {
        S_out->BeadType[new].Index[j] = S_in.BeadType[i].Index[j] +
                                        Count_old.Bead;
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
    S_out->Bead = realloc(S_out->Bead, sizeof (BEAD) * Count_out->Bead);
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
    S_out->Bonded = realloc(S_out->Bonded,
                            sizeof *S_out->Bonded * Count_out->Bonded);
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
    S_out->Unbonded = realloc(S_out->Unbonded,
                              sizeof *S_out->Unbonded * Count_out->Unbonded);
    for (int i = 0; i < Count_in->Unbonded; i++) {
      int new = i + Count_old.Unbonded;
      S_out->Unbonded[new] = S_in.Unbonded[i] + Count_old.Bead;
    }
  } //}}}
  // UnondedCoor //{{{
  if (Count_in->UnbondedCoor > 0) {
    Count_out->UnbondedCoor += Count_in->UnbondedCoor;
    S_out->UnbondedCoor = realloc(S_out->UnbondedCoor,
                                  sizeof *S_out->UnbondedCoor *
                                  Count_out->Unbonded);
    for (int i = 0; i < Count_in->UnbondedCoor; i++) {
      int new = i + Count_old.UnbondedCoor;
      S_out->UnbondedCoor[new] = S_in.UnbondedCoor[i] + Count_old.Bead;
    }
  } //}}}
  // BeadCoor //{{{
  if (Count_in->BeadCoor > 0) {
    Count_out->BeadCoor += Count_in->BeadCoor;
    S_out->BeadCoor = realloc(S_out->BeadCoor,
                                 sizeof *S_out->BeadCoor * Count_out->Bead);
    for (int i = 0; i < Count_in->BeadCoor; i++) {
      int new = i + Count_old.BeadCoor;
      S_out->BeadCoor[new] = S_in.BeadCoor[i] + Count_old.Bead;
    }
  } //}}}
  // MoleculeType //{{{
  if (Count_in->MoleculeType > 0) {
    Count_out->MoleculeType += Count_in->MoleculeType;
    S_out->MoleculeType = realloc(S_out->MoleculeType, sizeof (MOLECULETYPE) *
                                  Count_out->MoleculeType);
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
      // MoleculeType[].Bond // TODO + Count_old.BondType (like BType)
      if (mt_out->nBonds > 0) {
        mt_out->Bond = calloc(mt_out->nBonds, sizeof *mt_out->Bond);
        memcpy(mt_out->Bond, mt_in->Bond, sizeof *mt_in->Bond * mt_in->nBonds);
      }
      // MoleculeType[].Angle // TODO + Count_old.AngleType (like BType)
      if (mt_out->nAngles > 0) {
        mt_out->Angle = calloc(mt_out->nAngles, sizeof *mt_out->Angle);
        memcpy(mt_out->Angle, mt_in->Angle,
               sizeof *mt_in->Angle * mt_in->nAngles);
      }
      // MoleculeType[].Dihedral // TODO + Count_old.DihedralType (like BType)
      if (mt_out->nDihedrals > 0) {
        mt_out->Dihedral = calloc(mt_out->nDihedrals, sizeof *mt_out->Dihedral);
        memcpy(mt_out->Dihedral, mt_in->Dihedral,
               sizeof *mt_in->Dihedral * mt_in->nDihedrals);
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
    S_out->Molecule = realloc(S_out->Molecule, sizeof (MOLECULE) *
                                Count_out->Molecule);
    for (int i = 0; i < Count_in->Molecule; i++) {
      MOLECULE *mol_out = &S_out->Molecule[i+Count_old.Molecule],
               *mol_in = &S_in.Molecule[i];
      int type = mol_in->Type + Count_old.MoleculeType;
      *mol_out = *mol_in;
      mol_out->Type = type;
      mol_out->Index += Count_old.HighestResid;
      mol_out->Bead = malloc(sizeof *mol_out->Bead *
                             S_out->MoleculeType[type].nBeads);
      for (int j = 0; j < S_out->MoleculeType[type].nBeads; j++) {
        mol_out->Bead[j] = mol_in->Bead[j] + Count_old.Bead;
      }
    }
    Count_out->HighestResid += Count_in->HighestResid;
    S_out->Index_mol = realloc(S_out->Index_mol, sizeof *S_out->Index_mol *
                                 (Count_out->HighestResid+1));
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
      S_out->BondType = realloc(S_out->BondType, sizeof *S_out->BondType *
                                  Count_out->BondType);
      memcpy(S_out->BondType + Count_old.BondType, S_in.BondType,
             sizeof *S_out->BondType * Count_out->BondType);
  } //}}}
  // AngleType //{{{
  if (Count_in->AngleType > 0) {
      Count_out->AngleType += Count_in->AngleType;
      S_out->AngleType = realloc(S_out->AngleType,
                                   sizeof *S_out->AngleType *
                                   Count_out->AngleType);
      memcpy(S_out->AngleType + Count_old.AngleType, S_in.AngleType,
             sizeof *S_out->AngleType * Count_out->AngleType);
  } //}}}
  // DihedralType //{{{
  if (Count_in->DihedralType > 0) {
      Count_out->DihedralType += Count_in->DihedralType;
      S_out->DihedralType = realloc(S_out->DihedralType,
                                      sizeof *S_out->DihedralType *
                                      Count_out->DihedralType);
      memcpy(S_out->DihedralType + Count_old.DihedralType, S_in.DihedralType,
             sizeof *S_out->DihedralType * Count_out->DihedralType);
  } //}}}
} //}}}

// NewBeadType() //{{{
/*
 * Function to add a new bead type to a BEADTYPE struct
 * (and increment the number of bead types).
 */
void NewBeadType(BEADTYPE *BeadType[], int *number_of_types, char *name,
                 double charge, double mass, double radius) {
  int btype = *number_of_types;
  (*number_of_types)++;
  *BeadType = realloc(*BeadType, sizeof (BEADTYPE) * (btype + 1));
  strcpy((*BeadType)[btype].Name, name);
  (*BeadType)[btype].Number = 0;
  (*BeadType)[btype].Charge = charge;
  (*BeadType)[btype].Mass = mass;
  (*BeadType)[btype].Radius = radius;
}; //}}}
// NewMolType() //{{{
/*
 * Function to create a new molecule type in a MOLECULETYPE struct.
 */
void NewMolType(MOLECULETYPE *MoleculeType[], int *n_types, char *name,
                int n_beads, int n_bonds, int n_angles, int n_dihedrals) {
  int mtype = (*n_types)++;
  *MoleculeType = realloc(*MoleculeType,
                          sizeof (MOLECULETYPE) * (*n_types));
  // copy new name to MoleculeType[].Name
  strncpy((*MoleculeType)[mtype].Name, name, MOL_NAME);
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
    (*MoleculeType)[mtype].Dihedral = calloc(n_dihedrals,
                                       sizeof *(*MoleculeType)[mtype].Dihedral);
  }
  (*MoleculeType)[mtype].nBTypes = 0;
}; //}}}

void VerboseOutput(SYSTEM System) { //{{{
  PrintCount(System.Count);
  PrintBeadType(System);
  PrintMoleculeType(System);
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
    fprintf(stdout, ",\n  Bond Types:     %d", Count.BondType);
  }
  if (Count.AngleType > 0) {
    fprintf(stdout, ",\n  Angle Types:    %d", Count.AngleType);
  }
  if (Count.DihedralType > 0) {
    fprintf(stdout, ",\n  Dihedral Types: %d", Count.DihedralType);
  }
  fprintf(stdout, "\n\n");
} //}}}
void PrintBeadType(SYSTEM System) { //{{{
  // some stuff to properly align the fields //{{{
  int precision = 3, // number of decimal digits
      longest_name = 0, // longest bead type name
      most_beads = 0, // maximum number of beads
      max_q = 0, // maximum charge
      max_mass = 0, // maximum mass
      max_r = 0; // maximum radius
  bool negative = false; // extra space for '-' if there's negative charge
  for (int i = 0; i < System.Count.BeadType; i++) {
    int length = strlen(System.BeadType[i].Name);
    if (length > longest_name) {
      longest_name = length;
    }
    if (System.BeadType[i].Number > most_beads) {
      most_beads = System.BeadType[i].Number;
    }
    if (System.BeadType[i].Charge < 0) {
      negative = true;
    }
    if (System.BeadType[i].Charge != CHARGE && fabs(System.BeadType[i].Charge) > max_q) {
      max_q = floor(fabs(System.BeadType[i].Charge));
    }
    if (System.BeadType[i].Mass != MASS && System.BeadType[i].Mass > max_mass) {
      max_mass = floor(System.BeadType[i].Mass);
    }
    if (System.BeadType[i].Radius != RADIUS && System.BeadType[i].Radius > max_r) {
      max_r = floor(System.BeadType[i].Radius);
    }
  }
  // number of digits of the highest_number
  most_beads = floor(log10(most_beads)) + 1;
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
  if (max_mass == 0) {
    max_mass = 1;
  } else {
    max_mass = floor(log10(max_mass)) + 1 + precision + 1;
  }
  // number of digits of the radius
  if (max_r == 0) {
    max_r = 1;
  } else {
    max_r = floor(log10(max_mass)) + 1 + precision + 1;
  }
  // number of digits of the number of types
  int types_digits = floor(log10(System.Count.BeadType)) + 1;
  //}}}
  // print the information
  for (int i = 0; i < System.Count.BeadType; i++) {
    fprintf(stdout, "BeadType[%*d] = {", types_digits, i);
    fprintf(stdout, ".Name = %*s, ", longest_name, System.BeadType[i].Name);
    fprintf(stdout, ".Number = %*d, ", most_beads, System.BeadType[i].Number);
    // print charge
    fprintf(stdout, ".Charge = ");
    if (System.BeadType[i].Charge != CHARGE) {
      fprintf(stdout, "%*.*f, ", max_q, precision, System.BeadType[i].Charge);
    } else {
      for (int j = 0; j < (max_q-3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    // print mass
    fprintf(stdout, ".Mass = ");
    if (System.BeadType[i].Mass != MASS) {
      fprintf(stdout, "%*.*f, ", max_mass, precision, System.BeadType[i].Mass);
    } else {
      for (int j = 0; j < (max_mass-3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    fprintf(stdout, ".Radius = ");
    // print radius
    if (System.BeadType[i].Radius != RADIUS) {
      fprintf(stdout, "%*.*f", max_r, precision, System.BeadType[i].Radius);
    } else {
      for (int j = 0; j < (max_r-3); j++) {
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
        fprintf(stdout, "%d-%d", System.MoleculeType[i].Bond[j][0]+1,
                                 System.MoleculeType[i].Bond[j][1]+1);
        if (System.MoleculeType[i].Bond[j][2] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Bond[j][2]+1);
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
        fprintf(stdout, "%d-%d-%d", System.MoleculeType[i].Angle[j][0]+1,
                                    System.MoleculeType[i].Angle[j][1]+1,
                                    System.MoleculeType[i].Angle[j][2]+1);
        if (System.MoleculeType[i].Angle[j][3] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Angle[j][3]+1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print dihedrals if there are any //{{{
    if (System.MoleculeType[i].nDihedrals > 0) {
      fprintf(stdout, "  .nDihedrals = %d,\n  .Dihedral   = {", System.MoleculeType[i].nDihedrals);
      for (int j = 0; j < System.MoleculeType[i].nDihedrals; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d-%d", System.MoleculeType[i].Dihedral[j][0]+1,
                                       System.MoleculeType[i].Dihedral[j][1]+1,
                                       System.MoleculeType[i].Dihedral[j][2]+1,
                                       System.MoleculeType[i].Dihedral[j][3]+1);
        if (System.MoleculeType[i].Dihedral[j][4] != -1) {
          fprintf(stdout, "(%d)", System.MoleculeType[i].Dihedral[j][4]+1);
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
      fprintf(stdout, "%s", System.BeadType[System.MoleculeType[i].BType[j]].Name);
    } //}}}
    if (System.MoleculeType[i].Mass != MASS) {
      fprintf(stdout, "},\n  .Mass       = %.5f,\n", System.MoleculeType[i].Mass);
    } else {
      fprintf(stdout, "},\n  .Mass       = n/a,\n");
    }
    if (System.MoleculeType[i].Charge != CHARGE) {
      fprintf(stdout, "  .Charge     = %.5f\n}\n", System.MoleculeType[i].Charge);
    } else {
      fprintf(stdout, "  .Charge     = n/a\n}\n");
    }
  }
} //}}}
void PrintMolecule(SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.Molecule; i++) {
    int type = System.Molecule[i].Type;
    fprintf(stdout, "Molecule %3d (%d, %s):\n", i+1, System.Molecule[i].Index,
                                                System.MoleculeType[type].Name);
    fprintf(stdout, " BEAD INDICES (%d): ", System.MoleculeType[type].nBeads);
    fputs("intramolecular; input file\n", stdout);
    for (int j = 0; j < System.MoleculeType[type].nBeads; j++) {
      int id = System.Molecule[i].Bead[j];
      fprintf(stdout, "   %3d; %5d\n", j+1, id);
    }
  }
  fprintf(stdout, "\n");
} //}}}
void PrintBead(SYSTEM System) { //{{{
  fprintf(stdout, "Beads\n<input file id>, <bead type id> ");
  fprintf(stdout, "(<name>, <charge>, <mass>, <radius>), <molecule id>\n");
  for (int i = 0; i < System.Count.Bead; i++) {
    int type = System.Bead[i].Type;
    fprintf(stdout, "   %6d, type %3d (%8s, ", i, type,
                                               System.BeadType[type].Name);
    if (System.BeadType[type].Charge != CHARGE) {
      fprintf(stdout, "q=%5.2f, ", System.BeadType[type].Charge);
    } else {
      fprintf(stdout, "q=  n/a, ");
    }
    if (System.BeadType[type].Mass != MASS) {
      fprintf(stdout, "m=%5.2f, ", System.BeadType[type].Mass);
    } else {
      fprintf(stdout, "m=  n/a, ");
    }
    if (System.BeadType[type].Radius != RADIUS) {
      fprintf(stdout, "r=%5.2f)", System.BeadType[type].Radius);
    } else {
      fprintf(stdout, "r=  n/a)");
    }
    fprintf(stdout, ", molecule: ");
    if (System.Bead[i].Molecule == -1) {
      fprintf(stdout, "None\n");
    } else {
      int id = System.Bead[i].Molecule;
      fprintf(stdout, "%6d (%d)\n", System.Molecule[id].Index, id);
    }
  }
} //}}}

void FreeSystem(SYSTEM *System) { //{{{
  free(System->Index_mol);
  free(System->BeadCoor);
  free(System->Bonded);
  free(System->BondedCoor);
  free(System->Unbonded);
  free(System->UnbondedCoor);
  free(System->Bead);
  for (int i = 0; i < System->Count.BeadType; i++) {
    free(System->BeadType[i].Index);
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
};
void FreeMoleculeType(MOLECULETYPE *MoleculeType) {
  FreeMoleculeTypeEssentials(MoleculeType);
  if (MoleculeType->nBTypes > 0) {
    free(MoleculeType->BType);
  }
  free(MoleculeType->Index);
}
void FreeMoleculeTypeEssentials(MOLECULETYPE *MoleculeType) {
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
} //}}}

#if 0
// TODO redo
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
// TODO remove, I think
// CopyMolecule() //{{{
/*
 * Function to copy a MOLECULE struct into a new one. Memory management for
 * the output structure: sufficient memory is already allocated (mode=0), the
 * structure needs freeing and allocating (mode=1), no memory was yet allocated
 * at all (mode=2), or it just needs reallocating - only for cases when memory
 * only for the structure itself was allocated, not for the arrays within the
 * structure (mode=3).
 */
void CopyMolecule(int number_of_molecules, MOLECULETYPE *mt,
                  MOLECULE **m_out, MOLECULE *m_in, int mode) {
  // mt_out memory management //{{{
  switch (mode) {
    case 1:
      FreeMolecule(number_of_molecules, m_out);
      *m_out = malloc(sizeof (MOLECULE) * number_of_molecules);
      break;
    case 2:
      *m_out = malloc(sizeof (MOLECULE) * number_of_molecules);
      break;
    case 3:
      *m_out = realloc(*m_out, sizeof (MOLECULE) * number_of_molecules);
      break;
    default:
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "CopyMolecule()");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " function requires mode=0, 1, 2, or 3\n");
      ColourReset(STDERR_FILENO);
      exit(1);
  } //}}}
  for (int i = 0; i < number_of_molecules; i++) {
    (*m_out)[i] = m_in[i];
    int mtype = m_in[i].Type;
    (*m_out)[i].Bead = malloc(sizeof *(*m_out)[i].Bead * mt[mtype].nBeads);
    for (int j = 0; j < mt[mtype].nBeads; j++) {
      (*m_out)[i].Bead[j] = m_in[i].Bead[j];
    }
  }
} //}}}
// TODO remove
// RestorePBC_old() //{{{
/**
 * Function to restore removed periodic boundary conditions. Used also in case
 * of cell linked lists, because they need coordinates <0, BoxLength>. The
 * function requires orthogonal box, i.e., for triclinic box, the supplied
 * coordinates must first be transformed.
 */
void RestorePBC_old(COUNTS Counts, VECTOR BoxLength, BEAD **Bead) {

  for (int i = 0; i < Counts.BeadsCoor; i++) {
    // x direction
    while ((*Bead)[i].Position.x >= BoxLength.x) {
      (*Bead)[i].Position.x -= BoxLength.x;
    }
    while ((*Bead)[i].Position.x < 0) {
      (*Bead)[i].Position.x += BoxLength.x;
    }
    // y direction
    while ((*Bead)[i].Position.y >= BoxLength.y) {
      (*Bead)[i].Position.y -= BoxLength.y;
    }
    while ((*Bead)[i].Position.y < 0) {
      (*Bead)[i].Position.y += BoxLength.y;
    }
    // z direction
    while ((*Bead)[i].Position.z >= BoxLength.z) {
      (*Bead)[i].Position.z -= BoxLength.z;
    }
    while ((*Bead)[i].Position.z < 0) {
      (*Bead)[i].Position.z += BoxLength.z;
    }
  }
} //}}}
// VerboseOutput_old() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput_old(char *input_vcf, COUNTS Counts, VECTOR BoxLength,
                   BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {

  putchar('\n');
  if (BoxLength.x != -1) {
    fprintf(stdout, "Box size: %lf x %lf x %lf\n\n", BoxLength.x, BoxLength.y, BoxLength.z);
  } else {
    fprintf(stdout, "Unknown box size because no coordinate file is provided.\n\n");
  }
  PrintCounts_old(Counts);
  PrintBeadType2(Counts.TypesOfBeads, BeadType);
  PrintMoleculeType2(Counts.TypesOfMolecules, BeadType, MoleculeType);
  putchar('\n');
} //}}}
// VerboseOutput_oldish() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput_oldish(COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {
  PrintCounts_old(Counts);
  PrintBeadType2(Counts.TypesOfBeads, BeadType);
  PrintMoleculeType2(Counts.TypesOfMolecules, BeadType, MoleculeType);
} //}}}
// PrintCounts_old()  //{{{
/**
 * Function printing Counts structure.
 */
void PrintCounts_old(COUNTS Counts) {
  fprintf(stdout, "Counts = {\n");
  fprintf(stdout, "  .TypesOfBeads     = %d,\n", Counts.TypesOfBeads);
  fprintf(stdout, "  .Bonded           = %d,\n", Counts.Bonded);
  if (Counts.Bonded != Counts.BondedCoor) {
    fprintf(stdout, "  .BondedCoor       = %d,\n", Counts.BondedCoor);
  }
  fprintf(stdout, "  .Unbonded         = %d,\n", Counts.Unbonded);
  if (Counts.Unbonded != Counts.UnbondedCoor) {
    fprintf(stdout, "  .UnbondedCoor     = %d,\n", Counts.UnbondedCoor);
  }
  fprintf(stdout, "  .BeadsTotal       = %d,\n", Counts.BeadsTotal);
  if (Counts.BeadsTotal != Counts.BeadsCoor) {
    fprintf(stdout, "  .BeadsCoor        = %d,\n", Counts.BeadsCoor);
  }
  fprintf(stdout, "  .TypesOfMolecules = %d,\n", Counts.TypesOfMolecules);
  fprintf(stdout, "  .Molecules        = %d", Counts.Molecules);
  if (Counts.TypesOfBonds != -1) {
    fprintf(stdout, ",\n  .TypesOfBonds     = %d", Counts.TypesOfBonds);
  }
  if (Counts.TypesOfAngles != -1) {
    fprintf(stdout, ",\n  .TypesOfAngles    = %d", Counts.TypesOfAngles);
  }
  if (Counts.TypesOfDihedrals != -1) {
    fprintf(stdout, ",\n  .TypesOfDihedrals = %d", Counts.TypesOfDihedrals);
  }
  fprintf(stdout, "\n}\n\n");
} //}}}
// PrintBeadType_old()  //{{{
/**
 * Function printing BeadType structure.
 */
void PrintBeadType_old(COUNTS Counts, BEADTYPE *BeadType) {
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    fprintf(stdout, "BeadType[%2d] = {", i);
    fprintf(stdout, ".Name =%10s, ", BeadType[i].Name);
    fprintf(stdout, ".Number =%7d, ", BeadType[i].Number);
    fprintf(stdout, ".Charge = %lf, ", BeadType[i].Charge);
    fprintf(stdout, ".Mass = %lf}\n", BeadType[i].Mass);
//  fprintf(stdout, "Use = %3s, ", BeadType[i].Use? "Yes":"No");
//  fprintf(stdout, "Write = %3s}\n", BeadType[i].Write? "Yes":"No");
  }
  putchar('\n');
} //}}}
// PrintBeadType2()  //{{{
/**
 * Function printing BeadType structure.
 */
void PrintBeadType2(int number, BEADTYPE *BeadType) {
  // some stuff to properly align the fields //{{{
  int precision = 3, // number of decimal digits
      longest_name = 0, // longest bead type name
      most_beads = 0, // maximum number of beads
      max_q = 0, // maximum charge
      max_mass = 0, // maximum mass
      max_r = 0; // maximum radius
  bool negative = false; // extra space for '-' if there's negative charge
  for (int i = 0; i < number; i++) {
    int length = strlen(BeadType[i].Name);
    if (length > longest_name) {
      longest_name = length;
    }
    if (BeadType[i].Number > most_beads) {
      most_beads = BeadType[i].Number;
    }
    if (BeadType[i].Charge < 0) {
      negative = true;
    }
    if (BeadType[i].Charge != CHARGE && fabs(BeadType[i].Charge) > max_q) {
      max_q = floor(fabs(BeadType[i].Charge));
    }
    if (BeadType[i].Mass != MASS && BeadType[i].Mass > max_mass) {
      max_mass = floor(BeadType[i].Mass);
    }
    if (BeadType[i].Radius != RADIUS && BeadType[i].Radius > max_r) {
      max_r = floor(BeadType[i].Radius);
    }
  }
  // number of digits of the highest_number
  most_beads = floor(log10(most_beads)) + 1;
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
  if (max_mass == 0) {
    max_mass = 1;
  } else {
    max_mass = floor(log10(max_mass)) + 1 + precision + 1;
  }
  // number of digits of the radius
  if (max_r == 0) {
    max_r = 1;
  } else {
    max_r = floor(log10(max_mass)) + 1 + precision + 1;
  }
  // number of digits of the number of types
  int types_digits = floor(log10(number)) + 1;
  //}}}
  // print the information
  for (int i = 0; i < number; i++) {
    fprintf(stdout, "BeadType[%*d] = {", types_digits, i);
    fprintf(stdout, ".Name = %*s, ", longest_name, BeadType[i].Name);
    fprintf(stdout, ".Number = %*d, ", most_beads, BeadType[i].Number);
    // print charge
    fprintf(stdout, ".Charge = ");
    if (BeadType[i].Charge != CHARGE) {
      fprintf(stdout, "%*.*f, ", max_q, precision, BeadType[i].Charge);
    } else {
      for (int j = 0; j < (max_q-3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    // print mass
    fprintf(stdout, ".Mass = ");
    if (BeadType[i].Mass != MASS) {
      fprintf(stdout, "%*.*f, ", max_mass, precision, BeadType[i].Mass);
    } else {
      for (int j = 0; j < (max_mass-3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a, ");
    }
    fprintf(stdout, ".Radius = ");
    // print radius
    if (BeadType[i].Radius != RADIUS) {
      fprintf(stdout, "%*.*f", max_r, precision, BeadType[i].Radius);
    } else {
      for (int j = 0; j < (max_r-3); j++) {
        putchar(' ');
      }
      fprintf(stdout, "n/a");
    }
    fprintf(stdout, "}\n");
  }
  putchar('\n');
} //}}}
// PrintMoleculeType_old()  //{{{
/**
 * Function printing MoleculeType structure.
 */
void PrintMoleculeType_old(COUNTS Counts, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType) {
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    fprintf(stdout, "MoleculeType[%2d] = {\n", i);
    fprintf(stdout, "  .Name    = %s,\n", MoleculeType[i].Name);
    fprintf(stdout, "  .Number  = %d,\n", MoleculeType[i].Number);
    // print bead types (list all beads) //{{{
    fprintf(stdout, "  .nBeads  = %d,\n", MoleculeType[i].nBeads);
    fprintf(stdout, "  .Bead    = {");
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].Bead[j]].Name);
    }
    fprintf(stdout, "},\n"); //}}}
    // print bonds if there are any //{{{
    if (MoleculeType[i].nBonds > 0) {
      fprintf(stdout, "  .nBonds  = %d,\n", MoleculeType[i].nBonds);
      fprintf(stdout, "  .Bond    = {");
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d", MoleculeType[i].Bond[j][0]+1, MoleculeType[i].Bond[j][1]+1);
        if (MoleculeType[i].Bond[j][2] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Bond[j][2]+1);
        }
      }
    } //}}}
    // print angles if there are any //{{{
    if (MoleculeType[i].nAngles > 0) {
      fprintf(stdout, "},\n  .nAngles = %d,\n  .Angle   = {", MoleculeType[i].nAngles);
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d", MoleculeType[i].Angle[j][0]+1,
                                    MoleculeType[i].Angle[j][1]+1,
                                    MoleculeType[i].Angle[j][2]+1);
        if (MoleculeType[i].Angle[j][3] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Angle[j][3]+1);
        }
      }
    } //}}}
    // print bead types (just the which are present) //{{{
    fprintf(stdout, "},\n  .nBTypes = %d,\n  .BType   = {", MoleculeType[i].nBTypes);
    for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].BType[j]].Name);
    } //}}}
    if (MoleculeType[i].Mass != MASS) {
      fprintf(stdout, "},\n  .Mass    = %.5f,\n", MoleculeType[i].Mass);
    } else {
      fprintf(stdout, "},\n  .Mass    = n/a,\n");
    }
    if (MoleculeType[i].Charge != CHARGE) {
      fprintf(stdout, "  .Charge  = %.5f\n}\n", MoleculeType[i].Charge);
    } else {
      fprintf(stdout, "  .Charge  = n/a\n}\n");
    }
  }
} //}}}
// PrintMolecule_old() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule_old(int number_of_molecules,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                   BEADTYPE *BeadType, BEAD *Bead) {
  for (int i = 0; i < number_of_molecules; i++) {
    int type = Molecule[i].Type;
    fprintf(stdout, "Molecule %3d (%d, %s):\n", i+1, Molecule[i].Index,
                                                MoleculeType[type].Name);
    fprintf(stdout, " BEAD INDICES (%d): ", MoleculeType[type].nBeads);
    fputs("intramolecular; internal; input file\n", stdout);
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      fprintf(stdout, "   %3d; %5d\n", j+1, id);
    }
  }
  fprintf(stdout, "\n");
} //}}}
// PrintBead_old() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead_old(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead) {
  fprintf(stdout, "Beads - <i> (<Bead[i].Index>; <Index[i]>)\n");
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int type = Bead[i].Type;
    fprintf(stdout, "   %6d (%6d; %6d) %8s molecule: ", i, i, Index[i], BeadType[type].Name);
    if (Bead[i].Molecule == -1) {
      fprintf(stdout, "None\n");
    } else {
      fprintf(stdout, "%6d\n", Bead[i].Molecule+1);
    }
  }
} //}}}
// PrintBead2() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead2(int number_of_beads, int *Index,
                BEADTYPE *BeadType, BEAD *Bead) {
  fprintf(stdout, "Beads\n<input file id> (<internal id>), <bead type id> ");
  fprintf(stdout, "(<name>, <charge>, <mass>, <radius>), <molecule id>\n");
  for (int i = 0; i < number_of_beads; i++) {
    int type = Bead[i].Type;
    fprintf(stdout, "   %6d (%6d), type %3d (%8s, ", i, i, type,
                                                     BeadType[type].Name);
    if (BeadType[type].Charge != CHARGE) {
      fprintf(stdout, "q=%5.2f, ", BeadType[type].Charge);
    } else {
      fprintf(stdout, "q=  n/a, ");
    }
    if (BeadType[type].Mass != MASS) {
      fprintf(stdout, "m=%5.2f, ", BeadType[type].Mass);
    } else {
      fprintf(stdout, "m=  n/a, ");
    }
    if (BeadType[type].Radius != RADIUS) {
      fprintf(stdout, "r=%5.2f)", BeadType[type].Radius);
    } else {
      fprintf(stdout, "r=  n/a)");
    }
    fprintf(stdout, ", molecule: ");
    if (Bead[i].Molecule == -1) {
      fprintf(stdout, "None\n");
    } else {
      fprintf(stdout, "%6d\n", Bead[i].Molecule);
    }
  }
} //}}}
// FreeSystemInfo() //{{{
/**
 * Free memory for all standard arrays and structures of arrays.
 */
void FreeSystemInfo(COUNTS Counts, MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
                    BEADTYPE **BeadType, BEAD **Bead, int **Index) {
  free(*Index);
  FreeBead(Counts.BeadsTotal, Bead);
  free(*BeadType);
  FreeMolecule(Counts.Molecules, Molecule);
  FreeMoleculeType_old(Counts.TypesOfMolecules, MoleculeType);
} //}}}
// FreeSystemInfo2() //{{{
/**
 * Free memory for all standard arrays and structures of arrays.
 */
void FreeSystemInfo2(COUNTS Counts, MOLECULETYPE **MoleculeType,
                     MOLECULE **Molecule, int **Index_mol,
                     BEADTYPE **BeadType, BEAD **Bead, int **Index) {
  free(*Index);
  free(*Index_mol);
  FreeBead(Counts.BeadsTotal, Bead);
  free(*BeadType);
  FreeMolecule(Counts.Molecules, Molecule);
  FreeMoleculeType_old(Counts.TypesOfMolecules, MoleculeType);
} //}}}
// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(int number_of_beads, BEAD **Bead) {
  for (int i = 0; i < number_of_beads; i++) {
//  free((*Bead)[i].Aggregatexxx);
  }
  free(*Bead);
} //}}}
// FreeMoleculeType() //{{{
/**
 * Free memory allocated for MoleculeType struct array. This function makes
 * it easier other arrays to the MoleculeType struct in the future
 */
void FreeMoleculeType_old(int number_of_types, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    free((*MoleculeType)[i].Bead);
    if ((*MoleculeType)[i].nBonds > 0) {
      free(*(*MoleculeType)[i].Bond);
    }
    if ((*MoleculeType)[i].nAngles > 0) {
      free(*(*MoleculeType)[i].Angle);
    }
    if ((*MoleculeType)[i].nDihedrals > 0) {
      free(*(*MoleculeType)[i].Dihedral);
    }
    free((*MoleculeType)[i].BType);
  }
  free(*MoleculeType);
} //}}}
// FindBeadType2() //{{{
/* Function to identify type of bead from its name; returns -1 on non-existent
 * bead name.
 */
int FindBeadType2(char *name, int types_of_beads, BEADTYPE *BeadType) {
  int type;
  for (int i = 0; i < types_of_beads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in BeadType struct
  return -1;
} //}}}
// FindMoleculeType2() //{{{
int FindMoleculeType2(char *name, int number_of_types, MOLECULETYPE *MoleculeType) {
  int type;
  for (int i = 0; i < number_of_types; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return(-1);
} //}}}
// RemovePBCMolecules_old() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them.
 */
void RemovePBCMolecules_old(COUNTS Counts, VECTOR BoxLength,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {

  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      (*Bead)[Molecule[i].Bead[j]].Flag = false; // no beads moved yet
    }
    // first bead in the first bond is considered moved
    (*Bead)[Molecule[i].Bead[MoleculeType[type].Bond[0][0]]].Flag = true;
    bool done = false;
    int test = 0; // if too many loops, just leave the loop with error
    while (!done && test < 1000) {
      for (int j = 0; j < MoleculeType[type].nBonds; j++) {
        int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
        // move id1, if id2 is moved already
        if (!(*Bead)[id1].Flag && (*Bead)[id2].Flag) {
          VECTOR dist = Distance((*Bead)[id2].Position, (*Bead)[id1].Position, BoxLength);
          (*Bead)[id1].Position.x = (*Bead)[id2].Position.x - dist.x;
          (*Bead)[id1].Position.y = (*Bead)[id2].Position.y - dist.y;
          (*Bead)[id1].Position.z = (*Bead)[id2].Position.z - dist.z;
          (*Bead)[id1].Flag = true;
        // move id2, if id1 was moved already
        } else if ((*Bead)[id1].Flag && !(*Bead)[id2].Flag) {
          VECTOR dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position, BoxLength);

          (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
          (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
          (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
          (*Bead)[id2].Flag = true;
        }
      }

      // break while loop if all beads have moved
      done = true;
      for (int j = 1; j < MoleculeType[type].nBeads; j++) {
        if (!(*Bead)[Molecule[i].Bead[j]].Flag) {
          done = false;
          break;
        }
      }
      test++;
    }
    if (test == 1000) {
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: unable to 'join' molecule ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%s", MoleculeType[type].Name);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " (resid ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%d", i+1);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " )\n");
      ColourReset(STDERR_FILENO);
    }

    // put molecule's centre of mass into the simulation box //{{{
    VECTOR com = GeomCentre(MoleculeType[type].nBeads, Molecule[i].Bead, *Bead);
    // TODO: use math.h remainder()
    // by how many BoxLength's should com be moved?
    // for distant molecules - it shouldn't happen, but better safe than sorry
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
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int bead = Molecule[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * BoxLength.x;
      (*Bead)[bead].Position.y -= move.y * BoxLength.y;
      (*Bead)[bead].Position.z -= move.z * BoxLength.z;
    } //}}}
  }
} //}}}
// RemovePBCMolecules2() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them. The function requires orthogonal box, i.e.,
 * for triclinic box, the supplied coordinates must first be transformed.
 */
void RemovePBCMolecules2(COUNTS Counts, BOX Box,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {
  // TODO useless warning? //{{{
//for (int i = 0; i < Counts.TypesOfMolecules; i++) {
//  if (MoleculeType[i].nBonds == 0) {
//    ColourText(STDERR_FILENO, YELLOW);
//    fprintf(stderr, "\nWarning: molecule type ");
//    ColourText(STDERR_FILENO, CYAN);
//    fprintf(stderr, "%s", MoleculeType[i].Name);
//    ColourText(STDERR_FILENO, YELLOW);
//    fprintf(stderr, " has no bonds, so it cannot be 'joined'\n");
//    ColourReset(STDERR_FILENO);
//  }
//} //}}}
  // go through all molecules
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    // do nothing if the molecule has no bonds
    if (MoleculeType[type].nBonds == 0) {
      continue;
    }
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      (*Bead)[Molecule[i].Bead[j]].Flag = false; // no beads moved yet
    }
    // first bead in the first bond is considered moved
    (*Bead)[Molecule[i].Bead[MoleculeType[type].Bond[0][0]]].Flag = true;
    bool done = false;
    int test = 0; // if too many loops, just leave the loop with error
    while (!done && test < 1000) {
      for (int j = 0; j < MoleculeType[type].nBonds; j++) {
        int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
        // move id1, if id2 is moved already
        if (!(*Bead)[id1].Flag && (*Bead)[id2].Flag) {
          VECTOR dist = Distance((*Bead)[id2].Position, (*Bead)[id1].Position,
                                 Box.Length);
          (*Bead)[id1].Position.x = (*Bead)[id2].Position.x - dist.x;
          (*Bead)[id1].Position.y = (*Bead)[id2].Position.y - dist.y;
          (*Bead)[id1].Position.z = (*Bead)[id2].Position.z - dist.z;
          (*Bead)[id1].Flag = true;
        // move id2, if id1 was moved already
        } else if ((*Bead)[id1].Flag && !(*Bead)[id2].Flag) {
          VECTOR dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position,
                                 Box.Length);
          (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
          (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
          (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
          (*Bead)[id2].Flag = true;
        }
      }

      // break while loop if all beads have moved
      done = true;
      for (int j = 1; j < MoleculeType[type].nBeads; j++) {
        if (!(*Bead)[Molecule[i].Bead[j]].Flag) {
          done = false;
          break;
        }
      }
      test++;
    }
    if (test == 1000) {
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: unable to 'join' molecule");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%s", MoleculeType[type].Name);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " (resid ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%d", i+1);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " )\n");
      fprintf(stderr, "           Maybe not all beads are connected?\n");
      ColourReset(STDERR_FILENO);
    }

    // put molecule's centre of mass into the simulation box //{{{
    VECTOR com = GeomCentre(MoleculeType[type].nBeads, Molecule[i].Bead, *Bead);
    // by how many BoxLength's should com be moved?
    // for distant molecules - it shouldn't happen, but better safe than sorry
    INTVECTOR move;
    move.x = com.x / Box.Length.x;
    move.y = com.y / Box.Length.y;
    move.z = com.z / Box.Length.z;
    if (com.x < 0) {
      move.x--;
    }
    if (com.y < 0) {
      move.y--;
    }
    if (com.z < 0) {
      move.z--;
    }
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int bead = Molecule[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * Box.Length.x;
      (*Bead)[bead].Position.y -= move.y * Box.Length.y;
      (*Bead)[bead].Position.z -= move.z * Box.Length.z;
    } //}}}
  }
} //}}}
// RemovePBCMolecules_new() //{{{
/**
 * Function to remove periodic boundary conditions from all individual
 * molecules, thus joining them. The function requires orthogonal box, i.e.,
 * for triclinic box, the supplied coordinates must first be transformed.
 */
void RemovePBCMolecules_new(COUNTS Counts, BOX Box,
                        BEADTYPE *BeadType, BEAD **Bead,
                        MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {
  // go through all molecules
  for (int i = 0; i < Counts.Molecules; i++) {
    int type = Molecule[i].Type;
    // do nothing if the molecule has no bonds
    if (MoleculeType[type].nBonds == 0) {
      continue;
    }
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      (*Bead)[Molecule[i].Bead[j]].Flag = false; // no beads moved yet
    }
    // first bead in the first bond is considered moved
    for (int j = 0; j < MoleculeType[type].nBonds; j++) {
      int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]],
          id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
      if ((*Bead)[id1].InTimestep && (*Bead)[id2].InTimestep) {
        (*Bead)[id1].Flag = true;
        break;
      }
    }
    bool done = false;
    int test = 0; // if too many loops, just leave the loop with error
    while (!done && test < 1000) {
      for (int j = 0; j < MoleculeType[type].nBonds; j++) {
        int id1 = Molecule[i].Bead[MoleculeType[type].Bond[j][0]];
        int id2 = Molecule[i].Bead[MoleculeType[type].Bond[j][1]];
        if ((*Bead)[id1].InTimestep && (*Bead)[id2].InTimestep) {
          // move id1, if id2 is moved already
          if (!(*Bead)[id1].Flag && (*Bead)[id2].Flag) {
            VECTOR dist = Distance((*Bead)[id2].Position, (*Bead)[id1].Position,
                                   Box.Length);
            (*Bead)[id1].Position.x = (*Bead)[id2].Position.x - dist.x;
            (*Bead)[id1].Position.y = (*Bead)[id2].Position.y - dist.y;
            (*Bead)[id1].Position.z = (*Bead)[id2].Position.z - dist.z;
            (*Bead)[id1].Flag = true;
          // move id2, if id1 was moved already
          } else if ((*Bead)[id1].Flag && !(*Bead)[id2].Flag) {
            VECTOR dist = Distance((*Bead)[id1].Position, (*Bead)[id2].Position,
                                   Box.Length);
            (*Bead)[id2].Position.x = (*Bead)[id1].Position.x - dist.x;
            (*Bead)[id2].Position.y = (*Bead)[id1].Position.y - dist.y;
            (*Bead)[id2].Position.z = (*Bead)[id1].Position.z - dist.z;
            (*Bead)[id2].Flag = true;
          }
        }
      }

      // break while loop if all beads have moved
      done = true;
      for (int j = 1; j < MoleculeType[type].nBeads; j++) {
        int id = Molecule[i].Bead[j];
        if ((*Bead)[id].InTimestep && !(*Bead)[id].Flag) {
          done = false;
          break;
        }
      }
      test++;
    }
    if (test == 1000) {
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "\nWarning: unable to 'join' molecule");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%s", MoleculeType[type].Name);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " (resid ");
      ColourChange(STDERR_FILENO, CYAN);
      fprintf(stderr, "%d", i+1);
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, " )\n");
      fprintf(stderr, "           Maybe not all beads are connected?\n");
      ColourReset(STDERR_FILENO);
    }

    // put molecule's geometric centre into the simulation box //{{{
    VECTOR cog = GeomCentre(MoleculeType[type].nBeads, Molecule[i].Bead, *Bead);
    // by how many BoxLength's should cog be moved?
    // for distant molecules - it shouldn't happen, but better safe than sorry
    INTVECTOR move;
    move.x = cog.x / Box.Length.x;
    move.y = cog.y / Box.Length.y;
    move.z = cog.z / Box.Length.z;
    if (cog.x < 0) {
      move.x--;
    }
    if (cog.y < 0) {
      move.y--;
    }
    if (cog.z < 0) {
      move.z--;
    }
    for (int j = 0; j < MoleculeType[type].nBeads; j++) {
      int bead = Molecule[i].Bead[j];
      (*Bead)[bead].Position.x -= move.x * Box.Length.x;
      (*Bead)[bead].Position.y -= move.y * Box.Length.y;
      (*Bead)[bead].Position.z -= move.z * Box.Length.z;
    } //}}}
  }
} //}}}
// CopyBead_old() //{{{
/**
 * Function to copy BEAD structure into a new. Memory for the new one will be
 * reallocated in this function. Memory management for the output structure:
 * sufficient memory is already allocated (mode=0), the structure needs freeing
 * and allocating (mode=1), no memory wasn't yet allocated at all (mode=2), or
 * it needs reallocating - only for the case when the BEAD array is allocated,
 * but not the internal arrays (mode=3).
 */
void CopyBead_old(int number_of_beads, BEAD **b_out, BEAD *b_in, int mode) {
  // bt_out memory management //{{{
  switch (mode) {
    case 1:
      free(*b_out);
      *b_out = malloc(sizeof (BEAD) * number_of_beads);
      break;
    case 2:
      *b_out = malloc(sizeof (BEAD) * number_of_beads);
      break;
    case 3:
      *b_out = realloc(*b_out, sizeof (BEAD) * number_of_beads);
      break;
    default:
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "CopyBead()");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " function requires mode=0, 1, 2, or 3\n");
      ColourReset(STDERR_FILENO);
      exit(1);
  } //}}}
  for (int i = 0; i < number_of_beads; i++) {
    (*b_out)[i] = b_in[i];
//  (*b_out)[i].Aggregatexxx = calloc(1, sizeof *(*b_out)[i].Aggregatexxx);
  }
} //}}}
void FreeSystem_old(SYSTEM *System) { //{{{
  free((*System).Index_mol);
  free((*System).BeadCoor);
  // System.Bead
//for (int i = 0; i < (*System).nBeadsTotal; i++) {
//  free((*System).Bead[i].Aggregatexxx);
//}
  free((*System).Bead);
  // System.BeadType
  free((*System).BeadType);
  // System.Molecule
  for (int i = 0; i < (*System).Count.Molecule; i++) {
    free((*System).Molecule[i].Bead);
  }
  free((*System).Molecule);
  // System.MoleculeType
  for (int i = 0; i < (*System).Count.MoleculeType; i++) {
    free((*System).MoleculeType[i].Bead);
    if ((*System).MoleculeType[i].nBonds > 0) {
      free(*(*System).MoleculeType[i].Bond);
    }
    if ((*System).MoleculeType[i].nAngles > 0) {
      free(*(*System).MoleculeType[i].Angle);
    }
    if ((*System).MoleculeType[i].nDihedrals > 0) {
      free(*(*System).MoleculeType[i].Dihedral);
    }
    free((*System).MoleculeType[i].BType);
  }
  free((*System).MoleculeType);
  // bond & angle types
  if ((*System).Count.BondType > 0) {
    free((*System).BondType);
  }
  if ((*System).Count.AngleType > 0) {
    free((*System).AngleType);
  }
  if ((*System).Count.DihedralType > 0) {
    free((*System).DihedralType);
  }
}; //}}}
// CopySystem_old() //{{{
/*
 * Assumes unallocated SYSTEM.
 */
SYSTEM CopySystem_old(SYSTEM S_in) {
  SYSTEM S_out;
  InitSystem(&S_out);
  S_out.Box = S_in.Box;
  if (S_in.Count.Bead > 0) {
    S_out.Count = S_in.Count;
    // BeadType //{{{
    if (S_out.Count.BeadType > 0) {
      S_out.BeadType = realloc(S_out.BeadType,
                               sizeof (BEADTYPE) * S_out.Count.BeadType);
      memcpy(S_out.BeadType, S_in.BeadType,
             sizeof (BEADTYPE) * S_out.Count.BeadType);
    } else {
      strcpy(ERROR_MSG, "no bead types to copy; should never happen!");
      PrintWarning();
      return S_out;
    } //}}}
    // Bead //{{{
    S_out.Bead = realloc(S_out.Bead,
                         sizeof (BEAD) * S_out.Count.Bead);
    memcpy(S_out.Bead, S_in.Bead, sizeof (BEAD) * S_out.Count.Bead);
    //}}}
    // BeadsCoor //{{{
    if (S_out.Count.BeadCoor > 0) {
      S_out.BeadCoor = realloc(S_out.BeadCoor,
                                sizeof *S_out.BeadCoor * S_out.Count.Bead);
      memcpy(S_out.BeadCoor, S_in.BeadCoor,
             sizeof *S_out.BeadCoor * S_out.Count.Bead);
    } //}}}
    // Bonded //{{{
    if (S_out.Count.Bonded > 0) {
      S_out.Bonded = realloc(S_out.Bonded,
                             sizeof *S_out.Bonded * S_out.Count.Bonded);
      memcpy(S_out.Bonded, S_in.Bonded, sizeof *S_out.Bonded * S_out.Count.Bonded);
    } //}}}
    // BondedCoor //{{{
    if (S_out.Count.BondedCoor > 0) {
      S_out.BondedCoor = realloc(S_out.BondedCoor,
                                 sizeof *S_out.BondedCoor * S_out.Count.BondedCoor);
      memcpy(S_out.BondedCoor, S_in.BondedCoor,
             sizeof *S_out.BondedCoor * S_out.Count.BondedCoor);
    } //}}}
    // Unbonded //{{{
    if (S_out.Count.Unbonded > 0) {
      S_out.Unbonded = realloc(S_out.Unbonded, sizeof *S_out.Unbonded *
                               S_out.Count.Unbonded);
      memcpy(S_out.Unbonded, S_in.Unbonded,
             sizeof *S_out.Unbonded * S_out.Count.Unbonded);
    } //}}}
    // UnbondedCoor //{{{
    if (S_out.Count.UnbondedCoor > 0) {
      S_out.UnbondedCoor = realloc(S_out.UnbondedCoor,
                                   sizeof *S_out.UnbondedCoor *
                                   S_out.Count.UnbondedCoor);
      memcpy(S_out.UnbondedCoor, S_in.UnbondedCoor,
             sizeof *S_out.UnbondedCoor * S_out.Count.UnbondedCoor);
    } //}}}
    // MoleculeType //{{{
    if (S_out.Count.MoleculeType > 0) {
      S_out.MoleculeType = realloc(S_out.MoleculeType,
                                   sizeof (MOLECULETYPE) * S_out.Count.MoleculeType);
      for (int i = 0; i < S_out.Count.MoleculeType; i++) {
        S_out.MoleculeType[i] = S_in.MoleculeType[i];
        // MoleculeType[].Bead array
        if (S_out.MoleculeType[i].nBeads > 0) {
          S_out.MoleculeType[i].Bead = malloc(sizeof *S_out.MoleculeType[i].Bead *
                                              S_out.MoleculeType[i].nBeads);
          memcpy(S_out.MoleculeType[i].Bead, S_in.MoleculeType[i].Bead,
                 sizeof *S_out.MoleculeType[i].Bead *
                 S_out.MoleculeType[i].nBeads);
        } else {
          // should never happen
          strcpy(ERROR_MSG, "molecule type without beads");
          PrintWarning();
          fprintf(stderr, "%sMolecule %s%s%s\n\n", ErrCyan(), ErrYellow(),
                                                   S_out.MoleculeType[i].Name,
                                                   ErrColourReset());
        }
        // MoleculeType[].Bond array
        if (S_out.MoleculeType[i].nBonds > 0) {
          S_out.MoleculeType[i].Bond = malloc(sizeof *S_out.MoleculeType[i].Bond *
                                              S_out.MoleculeType[i].nBonds);
          memcpy(S_out.MoleculeType[i].Bond, S_in.MoleculeType[i].Bond,
                 sizeof *S_out.MoleculeType[i].Bond *
                 S_out.MoleculeType[i].nBonds);
        }
        // MoleculeType[].Angle array
        if (S_out.MoleculeType[i].nAngles > 0) {
          S_out.MoleculeType[i].Angle =
            malloc(sizeof *S_out.MoleculeType[i].Angle *
                   S_out.MoleculeType[i].nAngles);
          memcpy(S_out.MoleculeType[i].Angle, S_in.MoleculeType[i].Angle,
                 sizeof *S_out.MoleculeType[i].Angle *
                 S_out.MoleculeType[i].nAngles);
        }
        // MoleculeType[].Dihedral array
        if (S_out.MoleculeType[i].nDihedrals > 0) {
          S_out.MoleculeType[i].Dihedral =
            malloc(sizeof *S_out.MoleculeType[i].Dihedral *
                   S_out.MoleculeType[i].nDihedrals);
          memcpy(S_out.MoleculeType[i].Dihedral, S_in.MoleculeType[i].Dihedral,
                 sizeof *S_out.MoleculeType[i].Dihedral *
                 S_out.MoleculeType[i].nDihedrals);
        }
        // MoleculeType[].BType array
        if (S_out.MoleculeType[i].nBTypes > 0) {
          S_out.MoleculeType[i].BType =
            malloc(sizeof *S_out.MoleculeType[i].BType *
                   S_out.MoleculeType[i].nBTypes);
          memcpy(S_out.MoleculeType[i].BType, S_in.MoleculeType[i].BType,
                 sizeof *S_out.MoleculeType[i].BType *
                 S_out.MoleculeType[i].nBTypes);
        }
      }
    } //}}}
    // Molecule & Index_mol //{{{
    if (S_out.Count.Molecule > 0) {
      S_out.Molecule = realloc(S_out.Molecule,
                               sizeof (MOLECULE) * S_out.Count.Molecule);
      for (int i = 0; i < S_out.Count.Molecule; i++) {
        S_out.Molecule[i] = S_in.Molecule[i];
        // Molecule[].Bead array
        int type = S_out.Molecule[i].Type;
        if (S_out.MoleculeType[type].nBeads > 0) {
          S_out.Molecule[i].Bead = malloc(sizeof *S_out.Molecule[i].Bead *
                                          S_out.MoleculeType[type].nBeads);
          memcpy(S_out.Molecule[i].Bead, S_in.Molecule[i].Bead,
                 sizeof *S_out.Molecule[i].Bead *
                 S_out.MoleculeType[type].nBeads);
        }
      }
      // Index_mol
      S_out.Index_mol = realloc(S_out.Index_mol,
                                sizeof *S_out.Index_mol * S_out.Count.Molecule);
      memcpy(S_out.Index_mol, S_in.Index_mol,
             sizeof *S_out.Index_mol * S_out.Count.Molecule);
    } //}}}
    // BondType //{{{
    if (S_out.Count.BondType > 0) {
      S_out.BondType = realloc(S_out.BondType, sizeof (PARAMS) * S_out.Count.BondType);
      memcpy(S_out.BondType, S_in.BondType,
             sizeof (PARAMS) * S_out.Count.BondType);
    } //}}}
    // AngleType //{{{
    if (S_out.Count.AngleType > 0) {
      S_out.AngleType = realloc(S_out.AngleType,
                                sizeof (PARAMS) * S_out.Count.AngleType);
      memcpy(S_out.AngleType, S_in.AngleType,
             sizeof (PARAMS) * S_out.Count.AngleType);
    } //}}}
    // DihedralType //{{{
    if (S_out.Count.DihedralType > 0) {
      S_out.DihedralType = realloc(S_out.DihedralType,
                                   sizeof (PARAMS) * S_out.Count.DihedralType);
      memcpy(S_out.DihedralType, S_in.DihedralType,
             sizeof (PARAMS) * S_out.Count.DihedralType);
    } //}}}
  }
  return S_out;
} //}}}
void PruneSystem_old(SYSTEM *System) { //{{{
  SYSTEM S_old = CopySystem_old(*System);
  // re-initialize the system
  FreeSystem(System);
  InitSystem(System);
  // copy bead counts
  (*System).Box = S_old.Box;
  (*System).Count.Bead      = S_old.Count.BeadCoor;
  (*System).Count.BeadCoor = S_old.Count.BeadCoor;
  (*System).Count.Bonded     = S_old.Count.BondedCoor;
  (*System).Count.BondedCoor = S_old.Count.BondedCoor;
  (*System).Count.Unbonded     = S_old.Count.UnbondedCoor;
  (*System).Count.UnbondedCoor = S_old.Count.UnbondedCoor;
  // allocate memory for bead count arrays
  (*System).Bead = realloc((*System).Bead, sizeof (BEAD) *
                           (*System).Count.Bead);
  (*System).BeadCoor = realloc((*System).BeadCoor,
                                sizeof *(*System).BeadCoor *
                                (*System).Count.Bead);
  (*System).Bonded = realloc((*System).Bonded, sizeof *(*System).Bonded *
                             (*System).Count.Bonded);
  (*System).BondedCoor = realloc((*System).BondedCoor,
                                 sizeof *(*System).BondedCoor *
                                 (*System).Count.Bonded);
  (*System).Unbonded = realloc((*System).Unbonded, sizeof *(*System).Unbonded *
                               (*System).Count.Unbonded);
  (*System).UnbondedCoor = realloc((*System).UnbondedCoor,
                                   sizeof *(*System).UnbondedCoor *
                                   (*System).Count.Unbonded);
  // copy Bead/Unbonded/Bonded arrays & create new BeadType array //{{{
  int count_unbonded = 0, count_bonded = 0,
      *connect = calloc(S_old.Count.Bead, sizeof *connect);
  for (int i = 0; i < S_old.Count.BeadCoor; i++) {
    int old_id = S_old.BeadCoor[i];
    (*System).Bead[i] = S_old.Bead[old_id];
    (*System).BeadCoor[i] = i;
    connect[old_id] = i;
    if ((*System).Bead[i].Molecule == -1) {
      (*System).Unbonded[count_unbonded] = i;
      (*System).UnbondedCoor[count_unbonded] = i;
      count_unbonded++;
    } else {
      (*System).Bonded[count_bonded] = i;
      (*System).BondedCoor[count_bonded] = i;
      count_bonded++;
    }
    // create new bead type if it doesn't exist yet in the pruned system
    bool new = true;
    int old_type = S_old.Bead[old_id].Type;
    for (int j = 0; j < (*System).Count.BeadType; j++) {
      int new_type = FindBeadType(S_old.BeadType[old_type].Name, (*System));
      if (new_type != -1) {
        (*System).Bead[i].Type = new_type;
        (*System).BeadType[new_type].Number++;
        new = false;
        break;
      }
    }
    if (new) {
      int type = (*System).Count.BeadType;
      NewBeadType(&(*System).BeadType, &(*System).Count.BeadType,
                  S_old.BeadType[old_type].Name,
                  S_old.BeadType[old_type].Charge,
                  S_old.BeadType[old_type].Mass,
                  S_old.BeadType[old_type].Radius);
      (*System).BeadType[type].Number = 1;
      (*System).Bead[i].Type = type;
    }
  } //}}}
  // copy Molecule array & create a new MoleculeType array //{{{
  for (int i = 0; i < S_old.Count.Molecule; i++) {
    int old_type = S_old.Molecule[i].Type;
    // find if the molecule is to be in the pruned system
    int c_bead = 0;
    for (int j = 0; j < S_old.MoleculeType[old_type].nBeads; j++) {
      int id = S_old.Molecule[i].Bead[j];
      if (S_old.Bead[id].InTimestep) {
        c_bead++;
      }
    }
    if (c_bead > 0) { // should the molecule be in the pruned system?
      int new_id = (*System).Count.Molecule,
          old_type = S_old.Molecule[i].Type;
      (*System).Count.Molecule++;
      (*System).Molecule = realloc((*System).Molecule, sizeof (MOLECULE) *
                                   (*System).Count.Molecule);
      (*System).Molecule[new_id] = S_old.Molecule[i];
      (*System).Molecule[new_id].Bead =
        calloc(c_bead, sizeof *(*System).Molecule[new_id].Bead);
      c_bead = 0;
      for (int j = 0; j < S_old.MoleculeType[old_type].nBeads; j++) {
        int id = S_old.Molecule[i].Bead[j];
        if (S_old.Bead[id].InTimestep) {
          (*System).Molecule[new_id].Bead[c_bead] = connect[id];
          (*System).Bead[connect[id]].Molecule = new_id;
          c_bead++;
        }
      }
      // find if the molecule already exists in the pruned system
      bool new = true;
      for (int j = 0; j < (*System).Count.BeadType; j++) {
        int new_type = FindMoleculeType(S_old.MoleculeType[old_type].Name,
                                        (*System));
        if (new_type != -1) {
          (*System).Molecule[new_id].Type = new_type;
          (*System).MoleculeType[new_type].Number++;
          new = false;
          break;
        }
      }
      if (new) { // create new molecule type if it doesn't exist yet
        int type = (*System).Count.MoleculeType,
            c_bond = 0, c_angle = 0, c_dihedral = 0;
        // count bonds in the pruned molecule type //{{{
        for (int j = 0; j < S_old.MoleculeType[old_type].nBonds; j++) {
          int *id = BondIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep) {
            c_bond++;
          }
        } //}}}
//      printf("XXX c_bond=%d %d %s\n", c_bond, i, S_old.MoleculeType[old_type].Name);
        // count angles in the pruned molecule type //{{{
        for (int j = 0; j < S_old.MoleculeType[old_type].nAngles; j++) {
          int *id = AngleIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            c_angle++;
          }
        } //}}}
        // count dihedrals in the pruned molecule type //{{{
        for (int j = 0; j < S_old.MoleculeType[old_type].nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j);
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            c_dihedral++;
          }
        } //}}}
        NewMolType(&(*System).MoleculeType, &(*System).Count.MoleculeType,
                   S_old.MoleculeType[old_type].Name, c_bead, c_bond,
                   c_angle, c_dihedral);
        (*System).Molecule[new_id].Type = type;
        // copy beads to the new molecule type //{{{
        int *id_old_to_new = calloc(S_old.MoleculeType[old_type].nBeads,
                                    sizeof *id_old_to_new);
        for (int j = 0; j < S_old.MoleculeType[old_type].nBeads; j++) {
          id_old_to_new[j] = -1;
        }
        c_bead = 0;
        for (int j = 0; j < S_old.MoleculeType[old_type].nBeads; j++) {
          int id = S_old.Molecule[i].Bead[j],
              btype = S_old.MoleculeType[old_type].Bead[j];
          if (S_old.Bead[id].InTimestep) {
            (*System).MoleculeType[type].Bead[c_bead] =
              FindBeadType(S_old.BeadType[btype].Name, *System);
            id_old_to_new[j] = c_bead;
            c_bead++;
          }
        } //}}}
        // copy bonds to the new molecule type //{{{
        c_bond = 0;
        for (int j = 0; j < S_old.MoleculeType[old_type].nBonds; j++) {
          int *id = BondIndices(S_old, i, j),
              last = S_old.MoleculeType[old_type].Bond[j][2];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep) {
            // TODO the internal beads should have changed! - create array somewhere up there linking S_old.MoleculeType bead ids with (*System).MoleculeType bead ids...
            (*System).MoleculeType[type].Bond[c_bond][0] = id_old_to_new[id[2]];
            (*System).MoleculeType[type].Bond[c_bond][1] = id_old_to_new[id[3]];
            (*System).MoleculeType[type].Bond[c_bond][2] = last;
            c_bond++;
          }
        } //}}}
        // copy angles to the new molecule type //{{{
        c_angle = 0;
        for (int j = 0; j < S_old.MoleculeType[old_type].nAngles; j++) {
          int *id = AngleIndices(S_old, i, j),
              last = S_old.MoleculeType[old_type].Angle[j][3];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep) {
            (*System).MoleculeType[type].Angle[c_angle][0] =
              id_old_to_new[id[3]];
            (*System).MoleculeType[type].Angle[c_angle][1] =
              id_old_to_new[id[4]];
            (*System).MoleculeType[type].Angle[c_angle][2] =
              id_old_to_new[id[5]];
            (*System).MoleculeType[type].Angle[c_angle][3] = last;
            c_angle++;
          }
        } //}}}
        // copy bonds to the new molecule type //{{{
        c_bond = 0;
        for (int j = 0; j < S_old.MoleculeType[old_type].nDihedrals; j++) {
          int *id = DihedralIndices(S_old, i, j),
              last = S_old.MoleculeType[old_type].Dihedral[j][4];
          if (S_old.Bead[id[0]].InTimestep &&
              S_old.Bead[id[1]].InTimestep &&
              S_old.Bead[id[2]].InTimestep &&
              S_old.Bead[id[3]].InTimestep) {
            (*System).MoleculeType[type].Dihedral[c_dihedral][0] =
              id_old_to_new[id[4]];
            (*System).MoleculeType[type].Dihedral[c_dihedral][1] =
              id_old_to_new[id[5]];
            (*System).MoleculeType[type].Dihedral[c_dihedral][2] =
              id_old_to_new[id[6]];
            (*System).MoleculeType[type].Dihedral[c_dihedral][3] =
              id_old_to_new[id[7]];
            (*System).MoleculeType[type].Dihedral[c_dihedral][4] = last;
            c_dihedral++;
          }
        } //}}}
        free(id_old_to_new);
      }
    }
  }
  FillMolBTypes_old((*System).Count.MoleculeType, &(*System).MoleculeType);
  FillMolMassCharge_old((*System).Count.MoleculeType, &(*System).MoleculeType,
                    (*System).BeadType);
  //}}}
//PrintMoleculeType(*System);
//PrintMoleculeType(S_old);
//PrintMolecule(*System);
  FreeSystem(&S_old);
  free(connect);
} //}}}
// CopyMoleculeType_old() //{{{
/**
 * Function to copy MOLECULETYPE structure into a new one. Memory management
 * for the output structure: sufficient memory is already allocated (mode=0),
 * the structure needs freeing and allocating (mode=1), no memory was yet
 * allocated at all (mode=2), or it just needs reallocating - only for cases
 * when memory only for the structure itself was allocated, not for the arrays
 * within the structure (mode=3).
 */
// TODO if mode=0, then there shouldn't be any allocations at all; probably
//      remove this mode - when do I need straight out copy the arrays?
void CopyMoleculeType_old(int number_of_types, MOLECULETYPE **mt_out,
                      MOLECULETYPE *mt_in, int mode) {
  // mt_out memory management //{{{
  switch (mode) {
    case 1:
      FreeMoleculeType_old(number_of_types, mt_out);
      *mt_out = malloc(sizeof (MOLECULETYPE) * number_of_types);
      break;
    case 2:
      *mt_out = malloc(sizeof (MOLECULETYPE) * number_of_types);
      break;
    case 3:
      *mt_out = realloc(*mt_out, sizeof (MOLECULETYPE) * number_of_types);
      break;
    default:
      ErrorPrintError_old();
      ColourChange(STDERR_FILENO, YELLOW);
      fprintf(stderr, "CopyMoleculeType()");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " function requires mode=0, 1, 2, or 3\n");
      ColourReset(STDERR_FILENO);
      exit(1);
  } //}}}
  for (int i = 0; i < number_of_types; i++) {
    (*mt_out)[i] = mt_in[i]; // copy simple variables
    // allocate & copy Bead array
    (*mt_out)[i].Bead = malloc(sizeof *(*mt_out)[i].Bead * (*mt_out)[i].nBeads);
    for (int j = 0; j < (*mt_out)[i].nBeads; j++) {
      (*mt_out)[i].Bead[j] = mt_in[i].Bead[j];
    }
    // allocate & copy Bond array (if bonds are present)
    if ((*mt_out)[i].nBonds > 0) {
      (*mt_out)[i].nBonds = mt_in[i].nBonds;
      (*mt_out)[i].Bond = malloc(sizeof *(*mt_out)[i].Bond *
                                 (*mt_out)[i].nBonds);
      for (int j = 0; j < (*mt_out)[i].nBonds; j++) {
        for (int k = 0; k < 3; k++) {
          (*mt_out)[i].Bond[j][k] = mt_in[i].Bond[j][k];
        }
      }
    }
    // allocate & copy Angle array (if angles are present)
    if ((*mt_out)[i].nAngles > 0) {
      (*mt_out)[i].Angle = malloc(sizeof *(*mt_out)[i].Angle *
                                  (*mt_out)[i].nAngles);
      for (int j = 0; j < (*mt_out)[i].nAngles; j++) {
        for (int k = 0; k < 4; k++) {
          (*mt_out)[i].Angle[j][k] = mt_in[i].Angle[j][k];
        }
      }
    }
    // allocate & copy Dihedral array (if dihedrals are present)
    if ((*mt_out)[i].nDihedrals > 0) {
      (*mt_out)[i].Dihedral = malloc(sizeof *(*mt_out)[i].Dihedral *
                                     (*mt_out)[i].nDihedrals);
      for (int j = 0; j < (*mt_out)[i].nDihedrals; j++) {
        for (int k = 0; k < 5; k++) {
          (*mt_out)[i].Dihedral[j][k] = mt_in[i].Dihedral[j][k];
        }
      }
    }
  }
  FillMolBTypes_old(number_of_types, mt_out);
} //}}}
// FillMolBTypes_old() //{{{
/*
 * Function to fill MoleculeType[].BType array based on MoleculeType[].Bead array.
 */
void FillMolBTypes_old(int number_of_types, MOLECULETYPE **MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    (*MoleculeType)[i].nBTypes = 0;
    (*MoleculeType)[i].BType = malloc(sizeof *(*MoleculeType)[i].BType * 1);
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      bool new = true;
      for (int k = 0; k < (*MoleculeType)[i].nBTypes; k++) {
        if ((*MoleculeType)[i].Bead[j] == (*MoleculeType)[i].BType[k]) {
          new = false;
          break;
        }
      }
      if (new) {
        int type = (*MoleculeType)[i].nBTypes++;
        (*MoleculeType)[i].BType = realloc((*MoleculeType)[i].BType,
                                           sizeof *(*MoleculeType)[i].BType *
                                           (*MoleculeType)[i].nBTypes);
        (*MoleculeType)[i].BType[type] = (*MoleculeType)[i].Bead[j];
      }
    }
  }
} //}}}
// FindMoleculeType_old() //{{{
int FindMoleculeType_old(char *name, COUNTS Counts, MOLECULETYPE *MoleculeType) {
  int type;
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    if (strcmp(name, MoleculeType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in MoleculeType struct
  return(-1);
} //}}}
// FindBeadType_old() //{{{
/* Function to identify type of bead from its name; returns -1 on non-existent
 * bead name.
 */
int FindBeadType_old(char *name, COUNTS Counts, BEADTYPE *BeadType) {
  int type;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (strcmp(name, BeadType[i].Name) == 0) {
      type = i;
      return type;
    }
  }
  // name isn't in BeadType struct
  return -1;
} //}}}
// FreeMolecule() //{{{
/**
 * Free memory allocated for Molecule struct array. This function makes it
 * easier other arrays to the Molecule struct in the future
 */
void FreeMolecule(int number_of_molecules, MOLECULE **Molecule) {
  for (int i = 0; i < number_of_molecules; i++) {
    free((*Molecule)[i].Bead);
  }
  free(*Molecule);
} //}}}
// FillMolMassCharge_old() //{{{
/*
 * Function to calculate total mass and charge of molecules.
 */
void FillMolMassCharge_old(int number_of_types, MOLECULETYPE **MoleculeType,
                       BEADTYPE *BeadType) {
  for (int i = 0; i < number_of_types; i++) {
    // mass
    (*MoleculeType)[i].Mass = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      if (BeadType[btype].Mass != MASS) {
        (*MoleculeType)[i].Mass += BeadType[btype].Mass;
      } else {
        (*MoleculeType)[i].Mass = MASS;
        break;
      }
    }
    // charge
    (*MoleculeType)[i].Charge = 0;
    for (int j = 0; j < (*MoleculeType)[i].nBeads; j++) {
      int btype = (*MoleculeType)[i].Bead[j];
      if (BeadType[btype].Charge != CHARGE) {
        (*MoleculeType)[i].Charge += BeadType[btype].Charge;
      } else {
        (*MoleculeType)[i].Charge = CHARGE;
        break;
      }
    }
  }
} //}}}
// PrintMoleculeType2()  //{{{
/**
 * Function printing MoleculeType structure.
 */
void PrintMoleculeType2(int number_of_types, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType) {
  for (int i = 0; i < number_of_types; i++) {
    fprintf(stdout, "MoleculeType[%d] = {\n", i);
    fprintf(stdout, "  .Name       = %s,\n", MoleculeType[i].Name);
    fprintf(stdout, "  .Number     = %d,\n", MoleculeType[i].Number);
    // print bead types (list all beads) //{{{
    fprintf(stdout, "  .nBeads     = %d,\n", MoleculeType[i].nBeads);
    fprintf(stdout, "  .Bead       = {");
    for (int j = 0; j < MoleculeType[i].nBeads; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].Bead[j]].Name);
    }
    fprintf(stdout, "},\n"); //}}}
    // print bonds if there are any //{{{
    if (MoleculeType[i].nBonds > 0) {
      fprintf(stdout, "  .nBonds     = %d,\n", MoleculeType[i].nBonds);
      fprintf(stdout, "  .Bond       = {");
      for (int j = 0; j < MoleculeType[i].nBonds; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d", MoleculeType[i].Bond[j][0]+1,
                                 MoleculeType[i].Bond[j][1]+1);
        if (MoleculeType[i].Bond[j][2] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Bond[j][2]+1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print angles if there are any //{{{
    if (MoleculeType[i].nAngles > 0) {
      fprintf(stdout, "  .nAngles    = %d,\n", MoleculeType[i].nAngles);
      fprintf(stdout, "  .Angle      = {");
      for (int j = 0; j < MoleculeType[i].nAngles; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d", MoleculeType[i].Angle[j][0]+1,
                                    MoleculeType[i].Angle[j][1]+1,
                                    MoleculeType[i].Angle[j][2]+1);
        if (MoleculeType[i].Angle[j][3] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Angle[j][3]+1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print dihedrals if there are any //{{{
    if (MoleculeType[i].nDihedrals > 0) {
      fprintf(stdout, "  .nDihedrals = %d,\n  .Dihedral   = {", MoleculeType[i].nDihedrals);
      for (int j = 0; j < MoleculeType[i].nDihedrals; j++) {
        if (j != 0) {
          fprintf(stdout, ", ");
        }
        fprintf(stdout, "%d-%d-%d-%d", MoleculeType[i].Dihedral[j][0]+1,
                                       MoleculeType[i].Dihedral[j][1]+1,
                                       MoleculeType[i].Dihedral[j][2]+1,
                                       MoleculeType[i].Dihedral[j][3]+1);
        if (MoleculeType[i].Dihedral[j][4] != -1) {
          fprintf(stdout, "(%d)", MoleculeType[i].Dihedral[j][4]+1);
        }
      }
      fprintf(stdout, "},\n");
    } //}}}
    // print bead types (just the which are present) //{{{
    fprintf(stdout, "  .nBTypes    = %d\n", MoleculeType[i].nBTypes);
    fprintf(stdout, "  .BType      = {");
    for (int j = 0; j < MoleculeType[i].nBTypes; j++) {
      if (j != 0) {
        fprintf(stdout, ", ");
      }
      fprintf(stdout, "%s", BeadType[MoleculeType[i].BType[j]].Name);
    } //}}}
    if (MoleculeType[i].Mass != MASS) {
      fprintf(stdout, "},\n  .Mass       = %.5f,\n", MoleculeType[i].Mass);
    } else {
      fprintf(stdout, "},\n  .Mass       = n/a,\n");
    }
    if (MoleculeType[i].Charge != CHARGE) {
      fprintf(stdout, "  .Charge     = %.5f\n}\n", MoleculeType[i].Charge);
    } else {
      fprintf(stdout, "  .Charge     = n/a\n}\n");
    }
  }
} //}}}
#endif
