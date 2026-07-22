#include "Debug.h"
#include "Aggregates.h"
#include "Errors.h"
#include "General.h"
#include <errno.h>

// STATIC DECLARATIONS
static void PrintBeadHeader();
static void PrintOneBead(const SYSTEM System, const int id);
// find highest values from {Bond,Angle,Dihedral,Improper}Type params
static int * HighestParam(const SYSTEM System, const int type);

static void PrintBeadHeader() {
  fprintf(stdout, "Beads\n");
  fprintf(stdout, "<bead id>");
  fprintf(stdout, " (<bead type id>);");
  fprintf(stdout, " <molecule id>");
  fprintf(stdout, " (<molecule type id>);");
  fprintf(stdout, " <in coor>");
  putchar('\n');
}
static void PrintOneBead(const SYSTEM System, const int id) {
  BEAD *b = &System.Bead[id];
  fprintf(stdout, " %6d", id);
  fprintf(stdout, " (%3d);", b->Type);
  if (b->Molecule == -1) {
    fprintf(stdout, " %4s", "None");
    fprintf(stdout, "      ;");
  } else {
    fprintf(stdout, " %4d", System.Molecule[b->Molecule].Index);
    fprintf(stdout, " (%3d);", System.Molecule[b->Molecule].Type);
  }
  const char *in_ts = " no";
  if (b->InTimestep) {
    in_ts = "yes";
  }
  fprintf(stdout, " %s", in_ts);
  putchar('\n');
}
// find highest values from {Bond,Angle,Dihedral,Improper}Type params //{{{
static int * HighestParam(const SYSTEM System, const int type) {
  int num = 0;
  PARAMS *param;
  if (type == 0) {
    num = System.Count.BondType;
    param = System.BondType;
  } else if (type == 1) {
    num = System.Count.AngleType;
    param = System.AngleType;
  } else if (type == 2) {
    num = System.Count.DihedralType;
    param = System.DihedralType;
  } else if (type == 3) {
    num = System.Count.ImproperType;
    param = System.ImproperType;
  } else {
    err_msg("Highest(): 'type' must be 0 to 3");
  }
  PARAMS high = InitParams;
  for (int i = 0; i < num; i++) {
    if (param->a > high.a) {
      high.a = param->a;
    }
    if (param->b > high.b) {
      high.b = param->b;
    }
    if (param->c > high.c) {
      high.c = param->c;
    }
    if (param->d > high.d) {
      high.d = param->d;
    }
  }
  static int width[4];
  width[0] = snprintf(nullptr, 0, "%.5f", high.a);
  width[1] = snprintf(nullptr, 0, "%.0f", high.b);
  width[2] = snprintf(nullptr, 0, "%.0f", high.c);
  width[3] = snprintf(nullptr, 0, "%.0f", high.d);
  return width;
} //}}}

void VerboseOutput(const SYSTEM System) { //{{{
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
void PrintCount(const COUNT Count) { //{{{
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
  if (Count.Bond > 0) {
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
void PrintBeadType(const SYSTEM System) { //{{{
  // some stuff to properly align the fields //{{{
  int precision = 4;     // number of decimal digits
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
  int types_digits = floor(log10(System.Count.BeadType));
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
void PrintOneMolType(const SYSTEM System, const int n) { //{{{
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
void PrintAllMolTypes(const SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.MoleculeType; i++) {
    PrintOneMolType(System, i);
  }
  if (System.Count.MoleculeType > 0) {
    putchar('\n');
  }
} //}}}
void Print1Molecule(const SYSTEM System, const int n) { //{{{
  MOLECULE *mol = &System.Molecule[n];
  MOLECULETYPE *mtype = &System.MoleculeType[mol->Type];
  fprintf(stdout, "Molecule %3d (%d, %s):\n", n + 1, mol->Index, mtype->Name);
  fprintf(stdout, " BEAD INDICES (%d): ", mtype->nBeads);
  fputs("intramolecular; input file\n", stdout);
  for (int j = 0; j < mtype->nBeads; j++) {
    fprintf(stdout, "   %3d; %5d\n", j + 1, mol->Bead[j]);
  }
} //}}}
void PrintMolecules(const SYSTEM System) { //{{{
  for (int i = 0; i < System.Count.Molecule; i++) {
    Print1Molecule(System, i);
  }
  fprintf(stdout, "\n");
} //}}}
void PrintBead(const SYSTEM System) {
  PrintBeadHeader();
  for (int i = 0; i < System.Count.Bead; i++) {
    PrintOneBead(System, i);
  }
}
void PrintBeadCoor(const SYSTEM System) {
  PrintBeadHeader();
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    PrintOneBead(System, System.BeadCoor[i]);
  }
}
void PrintBondType(const SYSTEM System) { //{{{
  if (System.Count.BondType > 0) {
    // TODO: eventually, there should be more than harm
    int *wid = HighestParam(System, 0);
    fprintf(stdout, "Bond types");
    fprintf(stdout, " (lammps style 'harm')\n");
    for (int i = 0; i < System.Count.BondType; i++) {
      PARAMS *b = &System.BondType[i];
      fprintf(stdout, "  %*.5f %*.5f\n", wid[0], b->a, wid[1], b->b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintAngleType(const SYSTEM System) { //{{{
  if (System.Count.AngleType > 0) {
    // TODO: eventually, there should be more than cvff
    int *wid = HighestParam(System, 1);
    fprintf(stdout, "Angle types");
    fprintf(stdout, " (lammps style 'harm')\n");
    for (int i = 0; i < System.Count.AngleType; i++) {
      PARAMS *ang = &System.AngleType[i];
      fprintf(stdout, "  %*.5f %*.5f\n", wid[0], ang->a, wid[1], ang->b);
    }
    fprintf(stdout, "\n");
  }
} //}}}
void PrintDihedralType(const SYSTEM System) { //{{{
  if (System.Count.DihedralType > 0) {
    // TODO: eventually, there should be more than harm
    int *wid = HighestParam(System, 2);
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
void PrintImproperType(const SYSTEM System) { //{{{
  if (System.Count.ImproperType > 0) {
    // TODO: eventually, there should be more than cvff
    int *wid = HighestParam(System, 2);
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
void PrintBox(const BOX Box) { //{{{
  fprintf(stdout, "Box = {\n");
  if (Box.Low.x != 0 || Box.Low.y != 0 || Box.Low.z != 0) {
    fprintf(stdout, "  .Low = ( %lf %lf %lf )\n",
            Box.Low.x, Box.Low.y, Box.Low.z);
  }
  fprintf(stdout, "  .Length = ( %lf %lf %lf )\n",
          Box.Length.x, Box.Length.y, Box.Length.z);
  if (Box.alpha != 0 || Box.beta != 90 || Box.gamma != 90) {
    fprintf(stdout, "  .alpha = %lf\n", Box.alpha);
    fprintf(stdout, "  .beta  = %lf\n", Box.beta);
    fprintf(stdout, "  .gamma = %lf\n", Box.gamma);
    fprintf(stdout, "  .OrthoLength = ( %lf %lf %lf )\n", Box.OrthoLength.x,
                                                          Box.OrthoLength.y,
                                                          Box.OrthoLength.z);
    fprintf(stdout, "  .Bounding = ( %lf %lf %lf )\n", Box.Bounding.x,
                                                       Box.Bounding.y,
                                                       Box.Bounding.z);
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
void PrintStep(int *count_coor, const int start, const bool silent) { //{{{
  (*count_coor)++;
  if (!silent) {
    int saved_errno = errno;
    if (isatty(STDOUT_FILENO)) {
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
    errno = saved_errno;
  }
} //}}}
// TODO: coor & used -> STEP struct
void PrintLastStep(const int coor, const int used, const bool silent) { //{{{
  if (!silent) {
    int saved_errno = errno;
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    errno = saved_errno;
    fprintf(stdout, "Last Step: %d (used %d)\n", coor, used);
  }
} //}}}
void PrintAggregate(const SYSTEM System, const AGGREGATE *Aggregate) { //{{{
  const COUNT *Count = &System.Count;
  fprintf(stdout, "Aggregates: %d\n", Count->Aggregate);
  for (int i = 0; i < Count->Aggregate; i++) {
    // print molecules
    fprintf(stdout, " %d mols:", Aggregate[i].nMolecules);
    for (int j = 0; j < Aggregate[i].nMolecules; j++) {
      int mol = AggGetMol(&Aggregate[i], j);
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
