#include "AnalysisTools.h"
/* TODO: PrintBeadType() needs to find the longest name and the most number
 * of beads and insert white space accordingly
*/

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
  if ((*Box).alpha != 90 || (*Box).beta != 90 || (*Box).gamma != 90 ) {
    double a = (*Box).Length.x,
           b = (*Box).Length.y,
           c = (*Box).Length.z;
    double c_a = cos((*Box).alpha * PI / 180),
           c_b = cos((*Box).beta * PI / 180),
           c_g = cos((*Box).gamma * PI / 180),
           s_g = sin((*Box).gamma * PI / 180);
    double sqr = 1 - SQR(c_a) - SQR(c_b) - SQR(c_g) + 2 * c_a * c_b * c_g;
    if (sqr < 0) {
      strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
      return false;
    }
    double vol = a * b * c * sqrt(sqr);
    (*Box).Volume = vol;

    (*Box).transform[0][0] = a;
    (*Box).transform[1][0] = 0;
    (*Box).transform[2][0] = 0;
    (*Box).transform[0][1] = b * c_g;
    (*Box).transform[1][1] = b * s_g;
    (*Box).transform[2][1] = 0;
    (*Box).transform[0][2] = c * c_b;
    (*Box).transform[1][2] = c * (c_a - c_b * c_g) / s_g;
    (*Box).transform[2][2] = vol / (a * b * s_g);

    (*Box).inverse[0][0] = 1 / a;
    (*Box).inverse[1][0] = 0;
    (*Box).inverse[2][0] = 0;
    (*Box).inverse[0][1] = -c_g / (a * s_g);
    (*Box).inverse[1][1] = 1 / (b * s_g);
    (*Box).inverse[2][1] = 0;
    (*Box).inverse[0][2] = b * c * (c_g * (c_a - c_b * c_g) / (s_g * vol) -
                           c_b * s_g / vol);
    (*Box).inverse[1][2] = -a * c * (c_a - c_b * c_g) / (vol * s_g);
    (*Box).inverse[2][2] = a * b * s_g / vol;
    // tilt for triclinic axes (xy, xz, and yz) & lx, ly, and lz (lammps labels)
    (*Box).TriLength.x = a;
    (*Box).TriTilt[0] = b * c_g; // xy
    (*Box).TriTilt[1] = c * c_b; // xz
    sqr = SQR(b) - SQR((*Box).TriTilt[0]);
    if (sqr < 0) {
      strcpy(ERROR_MSG, "wrong dimensions for triclinic cell");
      return false;
    }
    (*Box).TriLength.y = sqrt(sqr);
    (*Box).TriTilt[2] = (b * c_a - (*Box).TriTilt[0] * (*Box).TriTilt[1]) /
                        (*Box).TriLength.y;
    sqr = SQR(c) - SQR((*Box).TriTilt[1]) - SQR((*Box).TriTilt[2]);
    if (sqr < 0) {
      ErrorPrintError_old();
      //TODO coloured output
      fprintf(stderr, "Error - wrong dimensions");
      exit(1);
    }
    (*Box).TriLength.z = sqrt(sqr);
    // make tilt component zero if they're close to zero
    for (int i = 0; i < 3; i++) {
      if (fabs((*Box).TriTilt[i]) < 0.001) {
        (*Box).TriTilt[i] = 0;
      }
    }
//  printf("\nLength:    %lf %lf %lf\n", (*Box).Length.x, (*Box).Length.y, (*Box).Length.z);
//  printf("TirLength: %lf %lf %lf\n", (*Box).TriLength.x, (*Box).TriLength.y, (*Box).TriLength.z);
//  printf("TriTilt:   %lf %lf %lf\n", (*Box).TriTilt[0], (*Box).TriTilt[1], (*Box).TriTilt[2]);
  //}}}
  } else { // orthogonal box //{{{
    (*Box).Volume = (*Box).Length.x * (*Box).Length.y * (*Box).Length.z;

    (*Box).transform[0][0] = (*Box).Length.x;
    (*Box).transform[1][0] = 0;
    (*Box).transform[2][0] = 0;
    (*Box).transform[0][1] = 0;
    (*Box).transform[1][1] = (*Box).Length.y;
    (*Box).transform[2][1] = 0;
    (*Box).transform[0][2] = 0;
    (*Box).transform[1][2] = 0;
    (*Box).transform[2][2] = (*Box).Length.z;

    (*Box).inverse[0][0] = 1 / (*Box).Length.x;
    (*Box).inverse[1][0] = 0;
    (*Box).inverse[2][0] = 0;
    (*Box).inverse[0][1] = 0;
    (*Box).inverse[1][1] = 1 / (*Box).Length.y;
    (*Box).inverse[2][1] = 0;
    (*Box).inverse[0][2] = 0;
    (*Box).inverse[1][2] = 0;
    (*Box).inverse[2][2] = 1 / (*Box).Length.z;
    // lx, ly, and lz for lammps
    (*Box).TriLength = (*Box).Length;
    // tilt for triclinic axes (xy, xz, and yz)
    (*Box).TriTilt[0] = 0;
    (*Box).TriTilt[1] = 0;
    (*Box).TriTilt[2] = 0;
  } //}}}
  // test print the matrices //{{{
//printf("Transformation matrix:\n");
//printf("   %lf %lf %lf\n", (*Box).transform[0][0], (*Box).transform[0][1], (*Box).transform[0][2]);
//printf("   %lf %lf %lf\n", (*Box).transform[1][0], (*Box).transform[1][1], (*Box).transform[1][2]);
//printf("   %lf %lf %lf\n", (*Box).transform[2][0], (*Box).transform[2][1], (*Box).transform[2][2]);
//printf("Inverse of transformation matrix:\n");
//printf("   %lf %lf %lf\n", (*Box).inverse[0][0], (*Box).inverse[0][1], (*Box).inverse[0][2]);
//printf("   %lf %lf %lf\n", (*Box).inverse[1][0], (*Box).inverse[1][1], (*Box).inverse[1][2]);
//printf("   %lf %lf %lf\n", (*Box).inverse[2][0], (*Box).inverse[2][1], (*Box).inverse[2][2]); //}}}
  return true;
} //}}}

// ToFractional() //{{{
void ToFractional(VECTOR *coor, BOX Box) {
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    VECTOR new = {0, 0, 0};
    new.x = Box.inverse[0][0] * (*coor).x +
            Box.inverse[0][1] * (*coor).y +
            Box.inverse[0][2] * (*coor).z;
    new.y = Box.inverse[1][0] * (*coor).x +
            Box.inverse[1][1] * (*coor).y +
            Box.inverse[1][2] * (*coor).z;
    new.z = Box.inverse[2][0] * (*coor).x +
            Box.inverse[2][1] * (*coor).y +
            Box.inverse[2][2] * (*coor).z;
    (*coor).x = new.x * Box.Length.x;
    (*coor).y = new.y * Box.Length.y;
    (*coor).z = new.z * Box.Length.z;
  }
} //}}}

// ToFractionalCoor() //{{{
void ToFractionalCoor(int number_of_beads, BEAD **Bead, BOX Box) {
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    for (int i = 0; i < number_of_beads; i++) {
      ToFractional(&(*Bead)[i].Position, Box);
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
void FromFractionalCoor(int number_of_beads, BEAD **Bead, BOX Box) {
  if (Box.alpha != 90 || Box.beta != 90 || Box.gamma != 90) {
    for (int i = 0; i < number_of_beads; i++) {
      (*Bead)[i].Position = FromFractional((*Bead)[i].Position, Box);
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
  PrintCounts(Counts);
  PrintBeadType2(Counts.TypesOfBeads, BeadType);
  PrintMoleculeType2(Counts.TypesOfMolecules, BeadType, MoleculeType);
  putchar('\n');
} //}}}
// TODO: print stuff from Counts only if > 0
// VerboseOutput() //{{{
/**
 * Function providing standard verbose output (for cases when verbose
 * option is used). It prints most of the information about used system.
 */
void VerboseOutput(COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
                   MOLECULETYPE *MoleculeType, MOLECULE *Molecule) {
  PrintCounts(Counts);
  PrintBeadType2(Counts.TypesOfBeads, BeadType);
  PrintMoleculeType2(Counts.TypesOfMolecules, BeadType, MoleculeType);
} //}}}

// PrintCounts()  //{{{
/**
 * Function printing Counts structure.
 */
void PrintCounts(COUNTS Counts) {
  fprintf(stdout, "Counts = {\n");
  fprintf(stdout, "  .TypesOfBeads     = %d,\n", Counts.TypesOfBeads);
  fprintf(stdout, "  .Bonded           = %d,\n", Counts.Bonded);
  fprintf(stdout, "  .Unbonded         = %d,\n", Counts.Unbonded);
  fprintf(stdout, "  .BeadsCoor        = %d,\n", Counts.BeadsCoor);
  if (Counts.BeadsTotal != Counts.BeadsCoor) {
    fprintf(stdout, "  .BeadsTotal       = %d,\n", Counts.BeadsTotal);
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

// PrintBeadType()  //{{{
/**
 * Function printing BeadType structure.
 */
void PrintBeadType(COUNTS Counts, BEADTYPE *BeadType) {
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

// PrintMoleculeType()  //{{{
/**
 * Function printing MoleculeType structure.
 */
void PrintMoleculeType(COUNTS Counts, BEADTYPE *BeadType, MOLECULETYPE *MoleculeType) {
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

// PrintBead() //{{{
/**
 * Function printing Bead structure.
 */
void PrintBead(COUNTS Counts, int *Index, BEADTYPE *BeadType, BEAD *Bead) {
  fprintf(stdout, "Beads - <i> (<Bead[i].Index>; <Index[i]>)\n");
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int type = Bead[i].Type;
    fprintf(stdout, "   %6d (%6d; %6d) %8s molecule: ", i, Bead[i].Index, Index[i], BeadType[type].Name);
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
    fprintf(stdout, "   %6d (%6d), type %3d (%8s, ", Bead[i].Index, i, type,
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

// PrintMolecule() //{{{
/**
 * Function printing Molecule structure.
 */
void PrintMolecule(int number_of_molecules,
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
      fprintf(stdout, "   %3d; %5d; %5d\n", j+1, id, Bead[id].Index);
    }
  }
  fprintf(stdout, "\n");
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
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
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
      fprintf(stdout, " %d (%d)", bead, Bead[bead].Index);
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

// FindBeadType() //{{{
/* Function to identify type of bead from its name; returns -1 on non-existent
 * bead name.
 */
int FindBeadType(char *name, COUNTS Counts, BEADTYPE *BeadType) {
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

// FindMoleculeType() //{{{
int FindMoleculeType(char *name, COUNTS Counts, MOLECULETYPE *MoleculeType) {
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

// FillMolBTypes //{{{
/*
 * Function to fill MoleculeType[].BType array based on MoleculeType[].Bead array.
 */
void FillMolBTypes(int number_of_types, MOLECULETYPE **MoleculeType) {
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

// FillMolMassCharge() //{{{
/*
 * Function to calculate total mass and charge of molecules.
 */
void FillMolMassCharge(int number_of_types, MOLECULETYPE **MoleculeType,
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
void RemovePBCMolecules(COUNTS Counts, BOX Box,
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
        int id = Molecule[mol].Bead[k];
        (*Bead)[id].nAggregates = 1;
        (*Bead)[id].Aggregate[0] = i;
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

// CopyBead() //{{{
/**
 * Function to copy BEAD structure into a new. Memory for the new one will be
 * reallocated in this function. Memory management for the output structure:
 * sufficient memory is already allocated (mode=0), the structure needs freeing
 * and allocating (mode=1), no memory wasn't yet allocated at all (mode=2), or
 * it needs reallocating - only for the case when the BEAD array is allocated,
 * but not the internal arrays (mode=3).
 */
void CopyBead(int number_of_beads, BEAD **b_out, BEAD *b_in, int mode) {
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
      fprintf(stderr, "CopyBeadType()");
      ColourChange(STDERR_FILENO, RED);
      fprintf(stderr, " function requires mode=0, 1, 2, or 3\n");
      ColourReset(STDERR_FILENO);
      exit(1);
  } //}}}
  for (int i = 0; i < number_of_beads; i++) {
    (*b_out)[i] = b_in[i];
    (*b_out)[i].Aggregate = calloc(1, sizeof *(*b_out)[i].Aggregate);
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

// CopyMoleculeType() //{{{
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
void CopyMoleculeType(int number_of_types, MOLECULETYPE **mt_out,
                      MOLECULETYPE *mt_in, int mode) {
  // mt_out memory management //{{{
  switch (mode) {
    case 1:
      FreeMoleculeType(number_of_types, mt_out);
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
  FillMolBTypes(number_of_types, mt_out);
} //}}}

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

// CopySystem() //{{{
/*
 * Function to copy the whole system - COUNTS, BEADTYPE, BEAD, MOLECULETYPE,
 * and MOLECULE structures and Index array. Memory management for the output:
 * sufficient memory is already allocated (mode=0), the structure needs freeing
 * and allocating (mode=1), no memory was yet allocated at all (mode=2), or it
 * just needs reallocating - only for cases when memory only for the structure
 * itself was allocated, not for the arrays within the structure (mode=3).
 */
void CopySystem(COUNTS *Counts_out, COUNTS Counts_in,
                BEADTYPE **bt_out, BEADTYPE *bt_in,
                BEAD **bead_out, BEAD *bead_in, int **index_out, int *index_in,
                MOLECULETYPE **mt_out, MOLECULETYPE *mt_in,
                MOLECULE **mol_out, MOLECULE *mol_in, int mode) {
  *Counts_out = Counts_in;
  // BEADTYPE
  CopyBeadType((*Counts_out).TypesOfBeads, bt_out, bt_in, mode);
  // BEAD & Index
  *bead_out = malloc(sizeof (BEAD) * (*Counts_out).BeadsCoor);
  *index_out = malloc(sizeof *index_out * (*Counts_out).BeadsCoor);
  for (int i = 0; i < (*Counts_out).BeadsCoor; i++) {
    (*bead_out)[i] = bead_in[i];
    (*bead_out)[i].Aggregate = malloc(sizeof *(*bead_out)[i].Aggregate *
                                   1); // just to free later
    (*index_out)[i] = index_in[i];
  }
  // MOLECULETYPE
  CopyMoleculeType((*Counts_out).TypesOfMolecules, mt_out, mt_in, mode);
  // MOLECULE
  CopyMolecule((*Counts_out).Molecules, mt_in, mol_out, mol_in, mode);
} //}}}

// FreeBead() //{{{
/**
 * Free memory allocated for Bead struct array. This function makes it
 * easier to add other arrays to the Bead struct in the future
 */
void FreeBead(int number_of_beads, BEAD **Bead) {
  for (int i = 0; i < number_of_beads; i++) {
    free((*Bead)[i].Aggregate);
  }
  free(*Bead);
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

// FreeMoleculeType() //{{{
/**
 * Free memory allocated for MoleculeType struct array. This function makes
 * it easier other arrays to the MoleculeType struct in the future
 */
void FreeMoleculeType(int number_of_types, MOLECULETYPE **MoleculeType) {
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
  FreeMoleculeType(Counts.TypesOfMolecules, MoleculeType);
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
  FreeMoleculeType(Counts.TypesOfMolecules, MoleculeType);
} //}}}
