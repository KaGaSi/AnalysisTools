#include "Read.h"
// TODO: since VtfCountStructLines is obsolete, we must somehow check for
//       invalid atom/bond lines in purely vcf file; or do we? It may not
//       actually matter to me; it will just be skipped - just write a note
//       into the manual
// TODO: perror() for file openings
// TODO test output of snprintf() to get rid of a warning; see https://stackoverflow.com/questions/51534284/how-to-circumvent-format-truncation-warning-in-gcc

// TODO: what's gonna happen with Counts.Beads & Counts.BeadsInVsf; I guess it
//       should stay and be used strictly - InVsf means in vsf, the other in
//       vcf
// TODO impropers vs dihedrals - add improper section to FIELD
/* TODO: PrintBeadType() needs to find the longest name and the most number
 * of beads and insert white space accordingly
*/

// VtfReadPBC() //{{{
/*
 * Function to get box dimensions from the provided coordinate file, that
 * is, search for a 'pbc <x> <y> <z> [<alpha> <beta> <gamma>]' line. The
 * provided file can either be a vcf coordinate only file or a full vtf one.
 */
void VtfReadPBC(char *input_vcf, BOX *Box) {
  // open the coordinate file
  FILE *coor;
  if ((coor = fopen(input_vcf, "r")) == NULL) {
    ErrorFileOpen(input_vcf, 'r');
    exit(1);
  }
  int file_line_count = 0;
  int words;
  char split[SPL_STR][SPL_LEN];
  // read input_vcf line by line
  while (ReadAndSplitLine(coor, &words, split)) {
    file_line_count++;
    int pbc = CheckVtLineType(words, split, false, input_vcf, file_line_count);
    // pbc line //{{{
    if (pbc == PBC_LINE || pbc == PBC_LINE_ANGLES) {
      (*Box).Length.x = atof(split[1]);
      (*Box).Length.y = atof(split[2]);
      (*Box).Length.z = atof(split[3]);
      // assume orthogonal box
      (*Box).alpha = 90;
      (*Box).beta = 90;
      (*Box).gamma = 90;
      // angles in pbc line; possibly triclinic cell
      if (pbc == PBC_LINE_ANGLES) {
        (*Box).alpha = atof(split[4]);
        (*Box).beta = atof(split[5]);
        (*Box).gamma = atof(split[6]);
      }
      break; //}}}
    // error - coordinate line //{{{
    } else if (pbc == COOR_LINE_I || pbc == COOR_LINE_O) {
      ErrorPrintError();
      ErrorPrintFileLine(input_vcf, file_line_count);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - encountered coordinate line before pbc line\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1); //}}}
    // error - unrecognised line //{{{
    } else if (pbc == ERROR_LINE) {
      ErrorPrintError();
      ErrorPrintFileLine(input_vcf, file_line_count);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - invalide line\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
  };
  fclose(coor);
} //}}}

// TODO not used anymore
// CountVtfStructLines() //{{{
/*
 * Function to count lines in the structure part of a vtf file (i.e., it
 * returns the number of the last line containing atom/bond keyword).
 */
int VtfCountStructLines(bool vtf, char *input) {
  if (vtf) {
    // open input file
    FILE *vtf;
    if ((vtf = fopen(input, "r")) == NULL) {
      ErrorFileOpen(input, 'r');
      exit(1);
    }
    int file_line_count = 0,
        last_line = -1;
    char split[SPL_STR][SPL_LEN];
    int words;
    while (ReadAndSplitLine(vtf, &words, split)) {
      int ltype = CheckVtLineType(words, split, false, input, file_line_count);
      if (ltype == TIME_LINE) {
        break;
      }
      file_line_count++;
      if (ltype == ATOM_LINE || ltype == BOND_LINE) {
        last_line = file_line_count;
      } else if (ltype == ERROR_LINE) {
        ErrorPrintError();
        ErrorPrintFileLine(input, file_line_count);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - invalide line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
    }
    fclose(vtf);
    return last_line;
  } else {
    return -1; // if not vtf, the number doesn't matter
  }
} //}}}

// SkipVtfStructure() //{{{
/*
 * Function to skip a structure part of a vtf file (i.e., combined vsf/vcf
 * file).
 */
void SkipVtfStructure(FILE *vcf, int struct_lines) {
  for (int i = 0; i < struct_lines; i++) {
    while (getc(vcf) != '\n')
      ;
  }
} //}}}

// NewBeadType() //{{{
/*
 * Function to add a new bead type to a BEADTYPE struct
 * (and increment the number of bead types).
 */
void NewBeadType(BEADTYPE **BeadType, int *number_of_types, char *name,
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
 * TODO: not used anywhere
 */
void NewMolType(char *name, int *number_of_types,
                MOLECULETYPE **MoleculeType, char *vsf_file) {
  int mtype = (*number_of_types)++;
  *MoleculeType = realloc(*MoleculeType,
                          sizeof (MOLECULETYPE) * (*number_of_types));
  // copy new name to MoleculeType[].Name
  strncpy((*MoleculeType)[mtype].Name, name, MOL_NAME);
  // initialize struct members
  (*MoleculeType)[mtype].Number = 0;
  (*MoleculeType)[mtype].nBonds = 0;
  (*MoleculeType)[mtype].nBeads = 0;
  (*MoleculeType)[mtype].nBTypes = 0;
}; //}}}

// FillMolMass //{{{
/*
 * Function to calculate mass of all molecules. If at least one bead has
 * undefined mass, the mass of the molecule is also undefined.
 */
void FillMolMass(int number_of_types,
                 BEADTYPE *BeadType, MOLECULETYPE **MoleculeType) {
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
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: molecule type ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%s", (*MoleculeType)[i].Name);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " has undefined mass\n");
      ResetColour(STDERR_FILENO);
    }
  }
} //}}}

// FillMolCharge //{{{
/*
 * Function to calculate charge of all molecules. If at least one bead has
 * undefined charge, the charge of the molecule is also undefined.
 */
void FillMolCharge(int number_of_types, BEADTYPE *BeadType,
                   MOLECULETYPE **MoleculeType) {
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
      YellowText(STDERR_FILENO);
      fprintf(stderr, "\nWarning: molecule type ");
      CyanText(STDERR_FILENO);
      fprintf(stderr, "%s", (*MoleculeType)[i].Name);
      YellowText(STDERR_FILENO);
      fprintf(stderr, " has undefined charge\n");
      ResetColour(STDERR_FILENO);
    }
  }
} //}}}

// FillMolType //{{{
/*
 * Function to fill BType array and mass and charge for each molecule type.
 */
void FillMolType(int number_of_types, BEADTYPE *BeadType,
                 MOLECULETYPE **MoleculeType) {
  FillMolBTypes(number_of_types, MoleculeType);
  FillMolMass(number_of_types, BeadType, MoleculeType);
  FillMolCharge(number_of_types, BeadType, MoleculeType);
} //}}}

// CheckVtfTimestep() //{{{
/*
 * Function to find what beads and molecules are in a timestep. For now, the
 * function is used only in FullVtfRead to determine these things only once for
 * all timesteps. Presumably, later something similar will be used to deremine
 * system composition of each step.
 */
bool CheckVtfTimestep(FILE *vcf, char *vcf_file, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

//PrintMoleculeType2((*Counts).TypesOfMolecules, *BeadType, *MoleculeType);
  bool indexed; // is the timestep indexed? ...returned by this function
  // skip timestep preamble and determine if a it's ordered/indexed //{{{
  char *stuff = calloc(LINE, sizeof *stuff);
  // TODO: just added boxlength for variable box size
  BOX Box;
  int lines = ReadVtfTimestepPreamble(&indexed, vcf_file, vcf, &stuff, &Box, true);
  free(stuff);
  for (int i = 0; i < lines; i++) {
    while (getc(vcf) != '\n')
      ;
  } //}}}
  // assume no beads are present //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    (*Bead)[i].Flag = false;
  } //}}}
  // count beads & save their ids //{{{
  lines = 0;
  while (true) {
    char split[SPL_STR][SPL_LEN], line[LINE];
    fgets(line, sizeof line, vcf);
    // if the line is too long, skip the rest of it
    if (strcspn(line, "\n") == (LINE-1)) {
      while (getc(vcf) != '\n')
        ;
    }
    // break loop on the end of the coordinate file
    if (feof(vcf)) {
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // break loop if the line isn't a coordinate line
    if (!CheckVtfCoordinateLine_old(words, split, line, indexed)) {
      // only blank, comment, or pbc line can follow the last coordinate line
      if (words == 0 || split[0][0] == '#' ||
          strcasecmp(split[0], "pbc") == 0) {
        break;
      } else {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", vcf_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - %s\n", line);
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
    }
    int id;
    if (indexed) { // indexed timestep => 1st is bead id
      id = (*Index)[atoi(split[0])];
    } else { // ordered timestep => bead id is based on the line count
      id = lines;
    }
    (*Bead)[id].Flag = true;
    lines++;
  } //}}}
  // error - not all beads present for ordered timestep //{{{
  if (!indexed && lines != (*Counts).Beads) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", vcf_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - ordered timestep requires all beads in the vcf file\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // if all beads are in the coordinate file, nothing more to do
  if (lines == (*Counts).Beads) {
    return indexed;
  }
  // otherwise, identify what's present
  COUNTS c_new = InitCounts;
  c_new.Beads = lines;
  c_new.BeadsInVsf = (*Counts).BeadsInVsf;
  // copy beads present in the timestep to a new BEAD struct & Index array //{{{
  int count = 0;
  BEAD *b_new = calloc(c_new.Beads, sizeof (BEAD));
  int *index_new = malloc(sizeof *index_new * (*Counts).BeadsInVsf);
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if ((*Bead)[i].Flag) {
      b_new[count] = (*Bead)[i];
      index_new[b_new[count].Index] = count;
      // count bonded & unbonded beads
      if (b_new[count].Molecule == -1) {
        c_new.Unbonded++;
      } else {
        c_new.Bonded++;
      }
      count++;
    }
  } //}}}
  // count numbers of beads of different types //{{{
  int *numbers = calloc((*Counts).TypesOfBeads, sizeof *numbers);
  for (int i = 0; i < c_new.Beads; i++) {
    int btype = b_new[i].Type;
    numbers[btype]++;
  } //}}}
  // count present bead types //{{{
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    // error - not all beads of a given type are in the timestep
    if (numbers[i] != 0 && numbers[i] != (*BeadType)[i].Number) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", vcf_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - not all ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", (*BeadType)[i].Name);
      RedText(STDERR_FILENO);
      fprintf(stderr, " beads are in the coordinate file\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    if (numbers[i] > 0) {
      c_new.TypesOfBeads++;
    }
  } //}}}
  // copy present bead types to a new BEADTYPE struct //{{{
  BEADTYPE *bt_new = calloc(c_new.TypesOfBeads, sizeof (BEADTYPE));
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    if (numbers[i] != 0) {
      bt_new[count] = (*BeadType)[i];
      count++;
    }
  }
  c_new.TypesOfBeads = count; //}}}
  free(numbers);
  // update bead type ids in the new BEAD according to the new BEADTYPE //{{{
  for (int i = 0; i < c_new.Beads; i++) {
    int btype = b_new[i].Type;
    btype = FindBeadType2((*BeadType)[btype].Name, c_new.TypesOfBeads, bt_new);
    b_new[i].Type = btype;
  } //}}}
  // find molecule types (and numbers of molecules) in the timestep //{{{
  /*
   * A molecule type is present in the 'new' system if it contains at least one
   * bead type present in the timestep.
   * 1) count number of molecule types and of molecules
   * 2) copy the molecule types present in the coordinate file to a new
   *    MOLECULETYPE struct
   * 3) prune the molecules to contain only beads present in the timestep
   *    3a) use only present beads
   *    3b) use only bonds between present beads
   *    3c) use only angles between present beads
   *        TODO: test that it works (once extra info is read from other files)
   *    3d) use only dihedrals between present beads
   *        TODO: test that it works (once extra info is read from other files)
   * 4) determine BTypes, charge, and mass of the 'new' molecule types (taking
   *    into account only present bead types)
   */
  // 1) //{{{
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int old_type = (*MoleculeType)[i].BType[j];
      if (FindBeadType2((*BeadType)[old_type].Name,
                        c_new.TypesOfBeads, bt_new) != -1) {
        c_new.TypesOfMolecules++;
        c_new.Molecules += (*MoleculeType)[i].Number;
        break;
      }
    }
  } //}}}
  // 2) & 3) //{{{
  MOLECULETYPE *mt_new = calloc(c_new.TypesOfMolecules, sizeof (MOLECULETYPE));
  count = 0;
  for (int i = 0; i < (*Counts).TypesOfMolecules; i++) {
    for (int j = 0; j < (*MoleculeType)[i].nBTypes; j++) {
      int old_type = (*MoleculeType)[i].BType[j];
      if (FindBeadType2((*BeadType)[old_type].Name,
                        c_new.TypesOfBeads, bt_new) != -1) {
        strcpy(mt_new[count].Name, (*MoleculeType)[i].Name);
        mt_new[count].Number = (*MoleculeType)[i].Number;
        // 3a) //{{{
        mt_new[count].nBeads = 0;
        mt_new[count].Bead = malloc(sizeof *mt_new[count].Bead * 1);
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          int btype = (*MoleculeType)[i].Bead[k];
          btype = FindBeadType2((*BeadType)[btype].Name,
                                c_new.TypesOfBeads, bt_new);
          if (btype != -1) {
            int id = mt_new[count].nBeads++;
            mt_new[count].Bead = realloc(mt_new[count].Bead,
                                         sizeof *mt_new[count].Bead *
                                         mt_new[count].nBeads);
            mt_new[count].Bead[id] = btype;
          }
        } //}}}
        // 3b) //{{{
        if ((*MoleculeType)[i].nBonds > 0) {
          mt_new[count].nBonds = 0;
          mt_new[count].Bond = malloc(sizeof *mt_new[count].Bond * 1);
          // go through all original types, testing each pair of beads
          for (int k = 0; k < (*MoleculeType)[i].nBonds; k++) {
            int id1 = (*MoleculeType)[i].Bond[k][0],
                id2 = (*MoleculeType)[i].Bond[k][1];
            int old_type1 = (*MoleculeType)[i].Bead[id1],
                old_type2 = (*MoleculeType)[i].Bead[id2];
            int btype1 = FindBeadType2((*BeadType)[old_type1].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype2 = FindBeadType2((*BeadType)[old_type2].Name,
                                       c_new.TypesOfBeads, bt_new);
            if (btype1 != -1 && btype2 != -1) {
              int id = mt_new[count].nBonds++;
              mt_new[count].Bond = realloc(mt_new[count].Bond,
                                           mt_new[count].nBonds *
                                           sizeof *mt_new[count].Bond);
              /* these indices correspond to the old MOLECULETYPE that contains
               * all beads - will be adjusted later */
              mt_new[count].Bond[id][0] = (*MoleculeType)[i].Bond[k][0];
              mt_new[count].Bond[id][1] = (*MoleculeType)[i].Bond[k][1];
              mt_new[count].Bond[id][2] = (*MoleculeType)[i].Bond[k][2];
            }
          }
          // if there are no bonds, free the array
          if (mt_new[count].nBonds == 0) {
            free(mt_new[count].Bond);
          }
        } //}}}
        // 3c) //{{{
        if ((*MoleculeType)[i].nAngles > 0) {
          mt_new[count].nAngles = 0;
          mt_new[count].Angle = malloc(sizeof *mt_new[count].Angle * 1);
          // go through all original types, testing each pair of beads
          for (int k = 0; k < (*MoleculeType)[i].nAngles; k++) {
            int id1 = (*MoleculeType)[i].Angle[k][0],
                id2 = (*MoleculeType)[i].Angle[k][1],
                id3 = (*MoleculeType)[i].Angle[k][2];
            int old_type1 = (*MoleculeType)[i].Bead[id1],
                old_type2 = (*MoleculeType)[i].Bead[id2],
                old_type3 = (*MoleculeType)[i].Bead[id3];
            int btype1 = FindBeadType2((*BeadType)[old_type1].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype2 = FindBeadType2((*BeadType)[old_type2].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype3 = FindBeadType2((*BeadType)[old_type3].Name,
                                       c_new.TypesOfBeads, bt_new);
            if (btype1 != -1 && btype2 != -1 && btype3 != -1) {
              int id = mt_new[count].nBonds++;
              mt_new[count].Angle = realloc(mt_new[count].Angle,
                                            sizeof *mt_new[count].Angle *
                                            mt_new[count].nAngles);
              /* these indices correspond to the old MOLECULETYPE that contains
               * all beads - will be adjusted later */
              mt_new[count].Angle[id][0] = (*MoleculeType)[i].Angle[k][0];
              mt_new[count].Angle[id][1] = (*MoleculeType)[i].Angle[k][1];
              mt_new[count].Angle[id][2] = (*MoleculeType)[i].Angle[k][2];
              mt_new[count].Angle[id][3] = (*MoleculeType)[i].Angle[k][3];
            }
          }
          // if there are no angles, free the array
          if (mt_new[count].nAngles == 0) {
            free(mt_new[count].Angle);
          }
        } //}}}
        // 3d) //{{{
        if ((*MoleculeType)[i].nDihedrals > 0) {
          mt_new[count].nDihedrals = 0;
          mt_new[count].Dihedral = malloc(sizeof *mt_new[count].Dihedral * 1);
          // go through all original types, testing each pair of beads
          for (int k = 0; k < (*MoleculeType)[i].nDihedrals; k++) {
            int id1 = (*MoleculeType)[i].Dihedral[k][0],
                id2 = (*MoleculeType)[i].Dihedral[k][1],
                id3 = (*MoleculeType)[i].Dihedral[k][2],
                id4 = (*MoleculeType)[i].Dihedral[k][3];
            int old_type1 = (*MoleculeType)[i].Bead[id1],
                old_type2 = (*MoleculeType)[i].Bead[id2],
                old_type3 = (*MoleculeType)[i].Bead[id3],
                old_type4 = (*MoleculeType)[i].Bead[id4];
            int btype1 = FindBeadType2((*BeadType)[old_type1].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype2 = FindBeadType2((*BeadType)[old_type2].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype3 = FindBeadType2((*BeadType)[old_type3].Name,
                                       c_new.TypesOfBeads, bt_new),
                btype4 = FindBeadType2((*BeadType)[old_type4].Name,
                                       c_new.TypesOfBeads, bt_new);
            if (btype1 != -1 && btype2 != -1 && btype3 != -1 && btype4 != -1) {
              int id = mt_new[count].nDihedrals++;
              mt_new[count].Dihedral = realloc(mt_new[count].Dihedral,
                                               sizeof *mt_new[count].Dihedral *
                                               mt_new[count].nDihedrals);
              /* these indices correspond to the old MOLECULETYPE that contains
               * all beads - will be adjusted later */
              mt_new[count].Dihedral[id][0] = (*MoleculeType)[i].Dihedral[k][0];
              mt_new[count].Dihedral[id][1] = (*MoleculeType)[i].Dihedral[k][1];
              mt_new[count].Dihedral[id][2] = (*MoleculeType)[i].Dihedral[k][2];
              mt_new[count].Dihedral[id][3] = (*MoleculeType)[i].Dihedral[k][3];
              mt_new[count].Dihedral[id][4] = (*MoleculeType)[i].Dihedral[k][4];
            }
          }
          // if there are no angles, free the array
          if (mt_new[count].nDihedrals == 0) {
            free(mt_new[count].Dihedral);
          }
        } //}}}
        // rebase bonded beads' indices according to the new MOLECULETYPE //{{{
        /*
         * All indices must go from 0 to (n-1) with n being the number of beads
         * in the molecule type
         * i) count how much to subtract from old indices in the new
         *    MOLECULETYPE's Bond array by going through the old MOLECULETYPE's
         *    index of bead types and testing if they are in the timestep
         * ii) subtract the required amount from each index in every bond
         * iii) subtract the required amount from each index in every angle
         * iv) subtract the required amount from each index in every dihedral
         */
        int subtract[(*MoleculeType)[i].nBeads];
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          subtract[k] = 0;
        }
        // i)
        for (int k = 0; k < (*MoleculeType)[i].nBeads; k++) {
          int old_type = (*MoleculeType)[i].Bead[k];
          int btype = FindBeadType2((*BeadType)[old_type].Name,
                                    c_new.TypesOfBeads, bt_new);
          if (btype == -1) {
            for (int l = k; l < (*MoleculeType)[i].nBeads; l++) {
              subtract[l]++;
            }
          }
        }
        // ii)
        for (int k = 0; k < mt_new[count].nBonds; k++) {
          int id1 = mt_new[count].Bond[k][0],
              id2 = mt_new[count].Bond[k][1];
          mt_new[count].Bond[k][0] -= subtract[id1];
          mt_new[count].Bond[k][1] -= subtract[id2];
        }
        // iii)
        for (int k = 0; k < mt_new[count].nAngles; k++) {
          int id1 = mt_new[count].Angle[k][0],
              id2 = mt_new[count].Angle[k][1],
              id3 = mt_new[count].Angle[k][2];
          mt_new[count].Angle[k][0] -= subtract[id1];
          mt_new[count].Angle[k][1] -= subtract[id2];
          mt_new[count].Angle[k][2] -= subtract[id3];
        }
        // iv)
        for (int k = 0; k < mt_new[count].nDihedrals; k++) {
          int id1 = mt_new[count].Dihedral[k][0],
              id2 = mt_new[count].Dihedral[k][1],
              id3 = mt_new[count].Dihedral[k][2],
              id4 = mt_new[count].Dihedral[k][3];
          mt_new[count].Dihedral[k][0] -= subtract[id1];
          mt_new[count].Dihedral[k][1] -= subtract[id2];
          mt_new[count].Dihedral[k][2] -= subtract[id3];
          mt_new[count].Dihedral[k][3] -= subtract[id4];
        } //}}}
        count++; // increment number of molecule types copied
        break;
      }
    }
  } //}}}
  // 4)
  FillMolBTypes(c_new.TypesOfMolecules, &mt_new);
  FillMolMass(c_new.TypesOfMolecules, bt_new, &mt_new);
  FillMolCharge(c_new.TypesOfMolecules, bt_new, &mt_new);
  //}}}
  // copy molecules present in the timestep to a new MOLECULE struct //{{{
  MOLECULE *m_new = calloc(c_new.Molecules, sizeof (MOLECULE));
  count = 0; // count molecules
  for (int i = 0; i < (*Counts).Molecules; i++) {
    int old_mtype = (*Molecule)[i].Type;
    int new_mtype = FindMoleculeType2((*MoleculeType)[old_mtype].Name,
                                      c_new.TypesOfMolecules, mt_new);
    if (new_mtype != -1) {
      m_new[count].Type = new_mtype;
      m_new[count].Bead = calloc(mt_new[new_mtype].nBeads,
                                 sizeof *m_new[count].Bead);
      int count2 = 0; // counts present beads in count-th molecule
      for (int j = 0; j < (*MoleculeType)[old_mtype].nBeads; j++) {
        int old_btype = (*MoleculeType)[old_mtype].Bead[j];
        if (FindBeadType2((*BeadType)[old_btype].Name,
                          c_new.TypesOfBeads, bt_new) != -1) {
          int id = (*Molecule)[i].Bead[j];
          int new_id = index_new[(*Bead)[id].Index];
          m_new[count].Bead[count2] = new_id;
          b_new[new_id].Molecule = count;
          count2++;
        }
      }
      count++;
    }
  } //}}}
  // copy data back to original structs & realloc the structs //{{{
  // BEAD struct - realloc and copy; redo Index array
  *Bead = realloc(*Bead, sizeof (BEAD) * c_new.Beads);
  for (int i = 0; i < c_new.Beads; i++) {
    (*Bead)[i] = b_new[i];
    (*Index)[(*Bead)[i].Index] = i;
  }
  // BEADTYPE struct - bt_new into BeadType
  CopyBeadType(c_new.TypesOfBeads, BeadType, bt_new, 3);
  // MOLECULETYPE - mt_new into MoleculeType
  CopyMoleculeType(c_new.TypesOfMolecules, MoleculeType, mt_new, 1);
  // MOLECULE - m_new into Molecule
  CopyMolecule(c_new.Molecules, *MoleculeType, Molecule, m_new, 1);
  // copy the new Counts struct back to the original one
  *Counts = c_new; //}}}
  // free memory
  FreeSystemInfo(*Counts, &mt_new, &m_new, &bt_new, &b_new, &index_new);
  return indexed;
} //}}}

// ReadVtfStructure() //{{{
/*
 * Function to read vtf structure file. It can recognize bead and molecule
 * types either according to name only (taking everything with the same name
 * for the first bead/molecule with that name) or according to all information
 * (name, mass, charge, and radius for bead types; bead order, bonds, angles,
 * and dihedrals for molecule types).
 */
/*
 * TODO: split into more functions and incorporate possible extra information
 * from lammps data/dl_meso FIELD files.
 */
void ReadVtfStructure(char *struct_file, bool detailed, COUNTS *Counts,
                      BEADTYPE **BeadType, BEAD **Bead, int **Index,
                      MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {

  (*Counts) = InitCounts; // zeroize
  // test if struct_file is a full vtf or just vsf //{{{
  bool struct_vtf = false;
  char *dot = strrchr(struct_file, '.');
  if (strcmp(dot, ".vtf") == 0) {
    struct_vtf = true;
  } //}}}
  // open struct_file //{{{
  FILE *vsf;
  if ((vsf = fopen(struct_file, "r")) == NULL) {
    ErrorFileOpen(struct_file, 'r');
    exit(1);
  } //}}}
  fpos_t pos;
  fgetpos(vsf, &pos); // save file pointer
  // 1) read through the structure part to find basic info: //{{{
  /*
   * i) check for errors
   * ii) find number of beads and molecules
   * iii) save unique bead and molecule names
   */
  int count_atom_lines = 0, // number of atom lines
      default_atom_line = -1, // line number of the first 'atom default' line
      count_bond_lines = 0, // number of bond lines
      count_comment_lines = 0, // number of comments (#) or blank lines
      atom_names = 0, res_names = 0, // number of unieque bead & molecule names
      last_struct = 0; // number of the last structure line (i.e., bond or atom)
  char (*atom_name)[BEAD_NAME+1] = malloc(sizeof *atom_name * 1);
  char (*res_name)[MOL_NAME+1] = malloc(sizeof *res_name * 1);
  while(true) {
    char error[LINE] = {'\0'}, // error message; no error for strlen(error) == 0
         split[SPL_STR][SPL_LEN], line[LINE];
    fgets(line, sizeof line, vsf);
    if (strcspn(line, "\n") == (LINE-1)) {
      while (getc(vsf) != '\n')
        ;
    }
//buffer[strcspn(buffer, "\r\n")] = '\0';
    if (feof(vsf)) { // break while loop on end of file
      break;
    }
    int words = SplitLine(split, line, "\t ");
    // in case of a vtf structure file, stop reading at the 'timestep' line
    if (struct_vtf && CheckVtfTimestepLine(words, split)) {
      break;
    }
    // skip blank, comment, and pbc lines
    if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
      count_comment_lines++;
      continue;
    }
    // check first character of the line for a[tom]/b[ond] line
    switch(split[0][0]) {
      case 'a': // atom line //{{{
        if (CheckVtfAtomLine_old(words, split, error)) {
          count_atom_lines++;
          last_struct = count_atom_lines +
                        count_bond_lines +
                        count_comment_lines;
        } else { // error - not a proper atom line //{{{
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", struct_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - %s\n", error);
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        // check if default line, otherwise count the number of beads //{{{
        if (strcmp(split[1], "default") == 0) {
          // warning - more than one 'default line'
          if (default_atom_line != -1) {
            YellowText(STDERR_FILENO);
            fprintf(stderr, "\nWarning: ");
            CyanText(STDERR_FILENO);
            fprintf(stderr, "%s", struct_file);
            YellowText(STDERR_FILENO);
            fprintf(stderr, " - multiple 'atom default lines'\n");
            fprintf(stderr, "         Using line number ");
            CyanText(STDERR_FILENO);
            fprintf(stderr, "%d\n", default_atom_line+1);
            ResetColour(STDERR_FILENO);
          } else { // save line number of the default line
            default_atom_line = count_atom_lines +
                                count_comment_lines +
                                count_bond_lines - 1; // 1st line should be 0th
          }
        } else {
          // check for highest index (i.e., the number of beads in vsf)
          if (atoi(split[1]) >= (*Counts).BeadsInVsf) {
            // vtf atom indices start at 0, so add 1
            (*Counts).BeadsInVsf = atoi(split[1]) + 1;
          }
        } //}}}
        // save unique bead and molecule names & count molecules //{{{
        for (int i = 0; i < words; i+= 2) {
          if (split[i][0] == 'n') { // atom name
            bool new = true; // assume the name is yet unknown
            for (int j = 0; j < atom_names; j++) {
              if (strncmp(split[i+1], atom_name[j], BEAD_NAME) == 0) {
                new = false;
                break;
              }
            }
            if (new) { // unknown (i.e., new) name
              atom_name = realloc(atom_name, sizeof *atom_name *
                                  (atom_names + 1));
              strncpy(atom_name[atom_names], split[i+1], BEAD_NAME);
              atom_names++;
            }
          } else if (strcmp(split[i], "resid") == 0) { // molecule index
            // check for highest index (i.e., the number of molecules)
            if (atoi(split[i+1]) >= (*Counts).Molecules) {
              // vtf resid indices start at 1, so don't add 1
              (*Counts).Molecules = atoi(split[i+1]);
            }
          } else if (strncmp(split[i], "res", 3) == 0) { // molecule name
            bool new = true; // assume the name is yet unknown
            for (int j = 0; j < res_names; j++) {
              if (strncmp(split[i+1], res_name[j], MOL_NAME) == 0) {
                new = false;
                break;
              }
            }
            if (new) { // unknown (i.e., new) name
              res_name = realloc(res_name, sizeof *res_name * (res_names + 1));
              strncpy(res_name[res_names], split[i+1], MOL_NAME);
              res_names++;
            }
          }
        } //}}}
        break; //}}}
      case 'b': // bond line //{{{
        if (CheckVtfBondLine_old(words, split, error)) {
          count_bond_lines++;
          last_struct = count_atom_lines +
                        count_bond_lines +
                        count_comment_lines;
        } else { // error - not a proper bond line
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", struct_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - %s\n", error);
          ErrorPrintLine(split, words);
          exit(1);
        }
        break; //}}}
      default: // unrecognised structure line //{{{
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", struct_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unrecognised line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1); //}}}
    }
  } //}}}
  // error - no default line and too few atom lines //{{{
  if (default_atom_line == -1 && count_atom_lines != (*Counts).BeadsInVsf) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", struct_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing atom line(s)\n");
    fprintf(stderr, "       (i.e., when 'atom default' is omitted, \
there must be a line for each atom)\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // total number of lines in the structure section //{{{
  /*
   * File type dependent: for a full vtf file, number of structure lines
   * corresponds to the last atom/bond lines; for a separate vsf file, the
   * number corresponds to the total number of lines in the file.
   */
  int total_lines = last_struct;
  if (!struct_vtf) {
    total_lines = count_atom_lines + count_bond_lines + count_comment_lines;
  } //}}}
  // structures and arrays to hold all vsf info //{{{
  // info from atom lines
  struct atom {
    int index, name, resid, resname;
    double charge, mass, radius;
  } *atom = calloc(count_atom_lines, sizeof *atom);
  for (int i = 0; i < count_atom_lines; i++) {
    atom[i].resid = -1; // i.e., not in a molecule
    atom[i].charge = CHARGE; // i.e., missing charge in the atom line
    atom[i].mass = MASS; // i.e., missing mass in the atom line
    atom[i].radius = RADIUS; // i.e., missing radius in the atom line
  }
  // array to connect atom ids in vsf with 'struct atom' numbering
  int *atom_id = malloc(sizeof *atom_id * (*Counts).BeadsInVsf);
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    atom_id[i] = -1; // if it stays -1, it's 'atom default'
  }
  // info from bond lines
  struct bond {
    int index1, index2; // indices from vsf file
  } *bond = calloc(count_bond_lines, sizeof *bond); //}}}
  // 2) save all vsf lines //{{{
  fsetpos(vsf, &pos); // restore file pointer to the beginning of struct_file
  int count_atoms = 0, count_bonds = 0;
  for (int line_no = 0; line_no < total_lines; line_no++) {
    char split[SPL_STR][SPL_LEN], line[LINE];
    fgets(line, sizeof line, vsf);
    // if the line is too long, skip the rest of it
    if (strcspn(line, "\n") == (LINE-1)) {
      while (getc(vsf) != '\n')
        ;
    }
    int words = SplitLine(split, line, "\t ");
    // skip blank, comment, and pbc lines
    if (words == 0 || split[0][0] == '#' || strcasecmp(split[0], "pbc") == 0) {
      continue;
    }
    switch (split[0][0]) {
      case 'a': // save atom line into atom struct //{{{
        // ignore multiple 'atom default' lines (use only the first one)
        if (strcmp(split[1], "default") != 0 || line_no == default_atom_line) {
          int id; // linetype-based bead index
          if (line_no == default_atom_line) { // 'atom default' line
            id = count_atom_lines - 1;
            atom[id].index = -1;
          } else { // 'atom <id>' line
            id = count_atoms++;
            atom_id[atoi(split[1])] = id;
            atom[id].index = atoi(split[1]);
          }
          // go through the line (first two strings are 'a[tom] <id>/default')
          for (int i = 2; i < words; i+=2) {
            if (strncmp(split[i], "resid", 3) == 0 &&
                strcmp(split[i], "resname") != 0) { // resid index
              atom[id].resid = atoi(split[i+1]) - 1; // vsf resid starts with 1
            } else if (split[i][0] == 'n') { // bead name
              for (int j = 0; j < atom_names; j++) {
                if (strncmp(split[i+1], atom_name[j], 16) == 0) {
                  atom[id].name = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "resname") == 0) { // molecule name
              for (int j = 0; j < res_names; j++) {
                if (strncmp(split[i+1], res_name[j], 8) == 0) {
                  atom[id].resname = j;
                  break;
                }
              }
            } else if (strcmp(split[i], "charge") == 0 ||
                       strcmp(split[i], "q") == 0) { // bead charge
              atom[id].charge = atof(split[i+1]);
            } else if (split[i][0] == 'm') { // bead mass
              atom[id].mass = atof(split[i+1]);
            } else if (split[i][0] == 'r' &&
                       strncmp(split[i], "res", 3) != 0) { // bead radius
              atom[id].radius = atof(split[i+1]);
            }
          }
        }
        break; //}}}
      case 'b': // save bond line into bond struct //{{{
        if (words == 2) {
          char index[SPL_STR][SPL_LEN];
          SplitLine(index, split[1], ":");
          bond[count_bonds].index1 = atoi(index[0]);
          bond[count_bonds].index2 = atoi(index[1]);
        } else {
          split[1][strlen(split[1])-1] = '\0';
          bond[count_bonds].index1 = atoi(split[1]);
          bond[count_bonds].index2 = atoi(split[2]);
        }
        count_bonds++;
        break; //}}}
    }
  } //}}}
  fclose(vsf); // everything is in memory now
  // error - bond with bead id higher than the number of beads //{{{
  for (int i = 0; i < count_bond_lines; i++) {
    if (bond[i].index1 >= (*Counts).BeadsInVsf ||
        bond[i].index2 >= (*Counts).BeadsInVsf) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", struct_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - bead index too high in ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "bond %d:%d\n", bond[i].index1, bond[i].index2);
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  // fill atom_id for 'atom default' beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (atom_id[i] == -1) {
      atom_id[i] = count_atom_lines - 1;
    }
  } //}}}
  // 3) identify bead types & fill BEADTYPE struct //{{{
  BEADTYPE *bt = malloc(sizeof (BEADTYPE) * 1);
  // save default type - if there is one
  if (atom[count_atom_lines-1].index == -1) {
    (*Counts).TypesOfBeads = 1;
    strcpy(bt[0].Name, atom_name[atom[count_atom_lines-1].name]);
    bt[0].Number = 0;
    bt[0].Charge = atom[count_atom_lines-1].charge;
    bt[0].Mass = atom[count_atom_lines-1].mass;
    bt[0].Radius = atom[count_atom_lines-1].radius;
  }
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
     * 1) identify unique names that will never merge
     * 2) for each unique name, find values of charge/mass/radius, noting
     *    ambiguities (i.e., when more than one well defined value exists)
     * 3) create 2D boolean array of size <unique names>*<unique names> to
     *    see what should be merged based on points 1) and 2):
     *    i) pick two bead types sharing a name (or the same bead type twice
     *       if it does not share a name with any other), say 'i' and 'k'.
     *    ii) check every bead type (say 'j') against i and k; if i and
     *        j should be merged (i.e., share a name), check k's value of
     *        diff_q/m/r - if it is a proper value, merge i and j; if not,
     *        merge i and j only if they have the same diff_q/m/r value.
     * 4) merge the types and count the number of unique types
     *    i) create a new type when a diagonal element of the array is true
     *    ii) check the remaining types against and merge those that should
     *        be merge with that new type, making the diagonal element for
     *        that merged type false so that no new type is created when
     *        its time comes in i)
     * 5) reorder the types so that types sharing the name are next to each
     *    other
     */
    // create (possibly too many) bead types according to bead properties //{{{
    for (int i = 0; i < count_atoms; i++) {
      int btype = -1;
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[j].Name, atom_name[atom[i].name]) == 0 &&
            bt[j].Charge == atom[i].charge &&
            bt[j].Mass == atom[i].mass &&
            bt[j].Radius == atom[i].radius) {
          btype = j;
        }
      }
      if (btype == -1) { // new bead type?
        NewBeadType(&bt, &(*Counts).TypesOfBeads, atom_name[atom[i].name],
                    atom[i].charge, atom[i].mass, atom[i].radius);
        btype = (*Counts).TypesOfBeads - 1;
      }
      bt[btype].Number++;
    } //}}}
//PrintBeadType2((*Counts).TypesOfBeads, bt);
    // count number of beads of default type if 'atom default' is present //{{{
    if (default_atom_line != -1) {
      bt[0].Number = (*Counts).BeadsInVsf;
      for (int i = 1; i < (*Counts).TypesOfBeads; i++) {
        bt[0].Number -= bt[i].Number;
      }
    } //}}}
//PrintBeadType2((*Counts).TypesOfBeads, bt);
    // Merging //{{{
    // 1) count and save unigue names //{{{
    int number_of_names = 0;
    char name[(*Counts).TypesOfBeads][BEAD_NAME+1];
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      bool exists = false;
      for (int j = 0; j < i; j++) {
        if (strcmp(bt[i].Name, bt[j].Name) == 0) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        strcpy(name[number_of_names], bt[i].Name);
        number_of_names++;
      }
    } //}}}
    // 2) //{{{
    // arrays for holding the value of charge/mass/radius for each bead type //{{{
    double diff_q[number_of_names],
           diff_m[number_of_names],
           diff_r[number_of_names]; //}}}
    // initialize arrays: assign values from the last type with each name //{{{
    for (int i = 0; i < number_of_names; i++) {
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(name[i], bt[j].Name) == 0) {
          diff_q[i] = bt[j].Charge;
          diff_m[i] = bt[j].Mass;
          diff_r[i] = bt[j].Radius;
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
    for (int i = 0; i < number_of_names; i++) {
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        // only consider type pairs with the same name
        if (strcmp(name[i], bt[j].Name) == 0) {
          // charge
          if (diff_q[i] != bt[j].Charge) {
            if (diff_q[i] != CHARGE && bt[j].Charge != CHARGE) {
              diff_q[i] = high;
            } else if (diff_q[i] == CHARGE) {
              diff_q[i] = bt[j].Charge;
            }
          }
          // mass
          if (diff_m[i] != bt[j].Mass) {
            if (diff_m[i] != MASS && bt[j].Mass != MASS) {
              diff_m[i] = high;
            } else if (diff_m[i] == MASS) {
              diff_m[i] = bt[j].Mass;
            }
          }
          // radius
          if (diff_r[i] != bt[j].Radius) {
            if (diff_r[i] != RADIUS && bt[j].Radius != RADIUS) {
              diff_r[i] = high;
            } else if (diff_r[i] == RADIUS) {
              diff_r[i] = bt[j].Radius;
            }
          }
        }
      }
    } //}}}
    //}}}
    /* test print diff_q/m/r //{{{
    for (int i = 0; i < number_of_names; i++) {
      printf("%s: q=%10.1f; m=%10.1f; r=%10.1f\n", name[i], diff_q[i], diff_m[i], diff_r[i]);
    }
    */ //}}}
    // 3) //{{{
    // initialize merge array by assuming nothing will be merged //{{{
    bool merge[(*Counts).TypesOfBeads][(*Counts).TypesOfBeads];
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        merge[i][j] = false; // 'i' and 'j' aren't to be merged
        merge[i][i] = true; // 'i' and 'i' is to be merged/copied
      }
    } //}}}
    // assume same-name bead types are to be merged //{{{
    for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[i].Name, bt[j].Name) == 0) {
          merge[i][j] = true;
        }
      }
    } //}}}
    // go through each bead type and compare it to two others //{{{
    for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
      // i)
      int k = 0;
      for (; k < (*Counts).TypesOfBeads; k++) {
        if (strcmp(name[k], bt[i].Name) == 0) {
          break;
        }
      }
      // ii)
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        // check charge
        if (merge[i][j]) {
          if (diff_q[k] == high) {
            if (bt[i].Charge == bt[j].Charge) {
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
            if (bt[i].Mass == bt[j].Mass) {
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
            if (bt[i].Radius == bt[j].Radius) {
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
    for (int i = 0; i < 5; i++) {
      printf("%s ", bt[i].Name);
    }
    putchar('\n');
    for (int i = 0; i < 5; i++) {
      printf("%s ", bt[i].Name);
      for (int j = 0; j < 5; j++) {
        printf(" %d", merge[i][j]);
      }
      putchar('\n');
    }
    */ //}}}
    // 4) //{{{
    BEADTYPE *temp = calloc((*Counts).TypesOfBeads, sizeof (BEADTYPE));
    int count = 0;
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      if (merge[i][i]) { // i)
        temp[count] = bt[i];
        for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
          if (merge[i][j]) { // ii)
            temp[count].Number += bt[j].Number;
            if (temp[count].Charge == CHARGE) {
              temp[count].Charge = bt[j].Charge;
            }
            if (temp[count].Mass == MASS) {
              temp[count].Mass = bt[j].Mass;
            }
            if (temp[count].Radius == RADIUS) {
              temp[count].Radius = bt[j].Radius;
            }
            merge[j][j] = false;
          }
        }
        count++;
      }
    }
    (*Counts).TypesOfBeads = count; //}}}
    //}}}
    // 5) //{{{
    // copy all bead types temporarily to bt struct
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      bt[i] = temp[i];
      bt[i].Use = false; // will indicate that it wasn't copied yet
    }
    // copy the bead types back to temp array in a proper order
    count = 0;
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      if (!bt[i].Use) {
        temp[count++] = bt[i];
        bt[i].Use = true;
        for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(bt[i].Name, bt[j].Name) == 0 && !bt[j].Use) {
            temp[count++] = bt[j];
            bt[j].Use = true;
          }
        }
      }
    }
    // finally, copy the types from the temporary array back to bt array
    for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
      bt[i] = temp[i];
    }
    free(temp); //}}}
    //}}}
  } else { // check only bead type name //{{{
    /*
     * charge/mass/radius is taken from the first bead with the given name
     * (even if undefined)
     */
    // go through all non-default atoms
    for (int i = 0; i < count_atoms; i++) {
      int btype = -1; // assume bead i is of new type
      for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[j].Name, atom_name[atom[i].name]) == 0) {
          btype = j; // nope, not a new type
        }
      }
      if (btype == -1) { // new bead type
        NewBeadType(&bt, &(*Counts).TypesOfBeads, atom_name[atom[i].name],
                    atom[i].charge, atom[i].mass, atom[i].radius);
        btype = (*Counts).TypesOfBeads - 1;
      }
      bt[btype].Number++;
    }
    // count number of beads of default type if atom default present
    if (atom[count_atom_lines-1].index == -1) {
      bt[0].Number = (*Counts).BeadsInVsf;
      for (int i = 1; i < (*Counts).TypesOfBeads; i++) {
        bt[0].Number -= bt[i].Number;
      }
    }
  } //}}}
  //}}}
  // 4) fill BEAD struct, putting unbonded beads first //{{{
  BEAD *bead_all = calloc((*Counts).BeadsInVsf, sizeof (BEAD));
  // connect to internal indexing (atom_id)
  int *index_all = malloc(sizeof *index_all * (*Counts).BeadsInVsf);
  // unbonded beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int internal_id = atom_id[i];
    if (atom[internal_id].resid == -1) {
      int id = (*Counts).Unbonded++;
      bead_all[id].Molecule = -1;
      index_all[id] = internal_id;
      if (detailed) { // find bead type based on all information //{{{
        bead_all[id].Type = -1;
        // 1) if a bead shares all values with a bead type, it is of that type
        for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0 &&
              bt[j].Charge == atom[internal_id].charge &&
              bt[j].Mass == atom[internal_id].mass &&
              bt[j].Radius == atom[internal_id].radius) {
            bead_all[id].Type = j;
            break;
          }
        }
        /*
         * 2) if no type was assigned, check for undefined values of the bead's
         *    charge/mass/radius
         */
        if (bead_all[id].Type == -1) {
          for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
            // only check if the bead type and bead share name
            if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0) {
              // check charge
              if (bt[j].Charge != atom[internal_id].charge &&
                  atom[internal_id].charge != CHARGE) {
                continue;
              }
              // check mass
              if (bt[j].Mass != atom[internal_id].mass &&
                  atom[internal_id].mass != MASS) {
                continue;
              }
              // check radius
              if (bt[j].Radius != atom[internal_id].radius &&
                  atom[internal_id].radius != RADIUS) {
                continue;
              }
              // assign bead type if all checks passed
              bead_all[id].Type = j;
              break;
            }
          }
        } //}}}
      } else { // find bead based only on name //{{{
        bead_all[id].Type = FindBeadType2(atom_name[atom[internal_id].name],
                                          (*Counts).TypesOfBeads, bt);
      } //}}}
      if (internal_id == (count_atom_lines-1)) { // default bead
        bead_all[id].Index = -1; // to be filled later
      } else {
        bead_all[id].Index = atom[internal_id].index;
      }
    }
  } //}}}
  // bonded beads //{{{
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    int internal_id = atom_id[i];
    if (atom[internal_id].resid != -1) {
      int id = (*Counts).Unbonded + (*Counts).Bonded; // place after unbonded
      (*Counts).Bonded++;
      bead_all[id].Molecule = atom[internal_id].resid;
      if (detailed) { // find bead type based on all information //{{{
        bead_all[id].Type = -1;
        // 1) if a bead shares all values with a bead type, it is of that type
        for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
          if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0 &&
              bt[j].Charge == atom[internal_id].charge &&
              bt[j].Mass == atom[internal_id].mass &&
              bt[j].Radius == atom[internal_id].radius) {
            bead_all[id].Type = j;
            break;
          }
        }
        /*
         * 2) if no type was assigned, check for undefined values of the bead's
         *    charge/mass/radius
         */
        if (bead_all[id].Type == -1) {
          for (int j = 0; j < (*Counts).TypesOfBeads; j++) {
            // only check if the bead type and bead share name
            if (strcmp(bt[j].Name, atom_name[atom[internal_id].name]) == 0) {
              // check charge
              if (bt[j].Charge != atom[internal_id].charge &&
                  atom[internal_id].charge != CHARGE) {
                continue;
              }
              // check mass
              if (bt[j].Mass != atom[internal_id].mass &&
                  atom[internal_id].mass != MASS) {
                continue;
              }
              // check radius
              if (bt[j].Radius != atom[internal_id].radius &&
                  atom[internal_id].radius != RADIUS) {
                continue;
              }
              // assign bead type if all checks passed
              bead_all[id].Type = j;
              break;
            }
          }
        } //}}}
      } else { // find bead based only on name //{{{
        bead_all[id].Type = FindBeadType2(atom_name[atom[internal_id].name],
                                          (*Counts).TypesOfBeads, bt);
      } //}}}
      bead_all[id].Index = atom[internal_id].index;
      index_all[id] = internal_id;
    }
  } //}}}
  // fill BEAD[].Index array for default beads //{{{
  // i) identify used indices (i.e., those explicitly written in the vsf file)
  bool *filled = calloc((*Counts).BeadsInVsf, sizeof *filled);
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (bead_all[i].Index != -1) {
      filled[bead_all[i].Index] = true;
    }
  }
  // ii) assign unused indices to beads of default type
  int count = 0;
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    if (bead_all[i].Index == -1) {
      for (; count < (*Counts).BeadsInVsf; count++) {
        if (!filled[count]) {
          bead_all[i].Index = count;
          count++;
          break;
        }
      }
    }
  }
  free(filled); //}}}
  //}}}
  // construct 'proper' Index array //{{{
  /*
   * until now, index_all was connected to atom struct, but from now on it will
   * be connected to the BEAD bead_all struct, because atom struct is no longer
   * used anywhere
   */
  for (int i = 0; i < (*Counts).BeadsInVsf; i++) {
    index_all[bead_all[i].Index] = i;
  } //}}}
  // if detailed, rename the bead types with the same name //{{{
  if (detailed) {
    for (int i = 0; i < ((*Counts).TypesOfBeads-1); i++) {
      count = 0;
      for (int j = (i+1); j < (*Counts).TypesOfBeads; j++) {
        if (strcmp(bt[i].Name, bt[j].Name) == 0) {
          count++;
          char name[BEAD_NAME+1];
          // shorten name if necessary
          if (count < 10) {
            strncpy(name, bt[j].Name, BEAD_NAME-2);
          } else if (count < 100) {
            strncpy(name, bt[j].Name, BEAD_NAME-3);
          }
//        P_IGNORE(-Wformat-truncation);
          // BEAD_NAME is max string length, i.e., array is longer
          snprintf(bt[j].Name, BEAD_NAME+1, "%s_%d", name, count);
//        P_POP;
        }
      }
    }
  } //}}}
  // 5) identify molecule types based on all data //{{{
  /*
   * Molecules of one type must share:
   * i) molecule name
   * ii) number of beads and bonds
   * iii) order of bead types
   * iv) connectivity
   */
  // count beads in each molecule for ii) //{{{
  int *atoms_per_mol = calloc((*Counts).Molecules, sizeof *atoms_per_mol);
  for (int i = (*Counts).Unbonded; i < (*Counts).BeadsInVsf; i++) {
    atoms_per_mol[bead_all[i].Molecule]++;
  } //}}}
  // error - resid numbering in vsf must be continuous //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    if (atoms_per_mol[i] == 0) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", struct_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - 'resid' numbering is not continuous\n");
      fprintf(stderr, "       Molecule ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", i);
      RedText(STDERR_FILENO);
      fprintf(stderr, " is missing\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  /*
  // test print number of beads in each molecule //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    printf("%d: %d\n", i, atoms_per_mol[i]);
  }
  putchar('\n'); //}}}
  */
  // allocate MOLECULE struct & and fill with indices for iv) //{{{
  MOLECULE *mol = calloc((*Counts).Molecules, sizeof (MOLECULE));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    mol[i].Bead = calloc(atoms_per_mol[i], sizeof *mol[i].Bead);
    mol[i].Aggregate = -1; // in no aggregate
  }
  // fill to Molecule[].Bead array with bead indices
  int *count_in_mol = calloc((*Counts).Molecules, sizeof *count_in_mol);
  for (int i = (*Counts).Unbonded; i < (*Counts).BeadsInVsf; i++) {
    int m_id = bead_all[i].Molecule;
    mol[m_id].Bead[count_in_mol[m_id]] = i;
    count_in_mol[m_id]++;
  }
  // pro-forma: test that the two bead counts in individual molecules agree
  // unless things are somehow weirdly bad, the error is never triggered
  for (int i = 0; i < (*Counts).Molecules; i++) {
    if (count_in_mol[i] != atoms_per_mol[i]) {
      ErrorPrintError();
      RedText(STDERR_FILENO);
      fprintf(stderr, "something wrong with bead count in molecule");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", i);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d!=%d\n", count_in_mol[i], atoms_per_mol[i]);
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  } //}}}
  /*
  // test print beads in molecules //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    printf("%d:", i);
    for (int j = 0; j < atoms_per_mol[i]; j++) {
      printf(" %d", mol[i].Bead[j]);
    }
    putchar('\n');
  }
  putchar('\n'); //}}}
  */
  // count bonds in all molecules for ii) //{{{
//P_IGNORE(-Walloc-size-larger-than=);
  int *bonds_per_mol = calloc((*Counts).Molecules, sizeof *bonds_per_mol);
//P_POP
  for (int i = 0; i < count_bond_lines; i++) {
    int m_id1 = bead_all[index_all[bond[i].index1]].Molecule;
    int m_id2 = bead_all[index_all[bond[i].index2]].Molecule;
    // error - beads from one bond are in different molecules
    if (m_id1 != m_id2) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", struct_file);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - bonded beads in different molecules (bead ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", bond[i].index1);
      RedText(STDERR_FILENO);
      fprintf(stderr, " in molecule ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", m_id1+1); // resid start at 1 in vsf
      RedText(STDERR_FILENO);
      fprintf(stderr, " and bead ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", bond[i].index2);
      RedText(STDERR_FILENO);
      fprintf(stderr, " in molecule ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%d", m_id2+1); // resid start at 1 in vsf
      RedText(STDERR_FILENO);
      fprintf(stderr, ")\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    bonds_per_mol[m_id1]++;
  } //}}}
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
  } *molec = malloc(sizeof *molec * (*Counts).Molecules);
  for (int i = 0; i < (*Counts).Molecules; i++) {
    count_in_mol[i] = 0;
    molec[i].connect = malloc(sizeof *molec[i].connect * bonds_per_mol[i]);
    for (int j = 0; j < bonds_per_mol[i]; j++) {
      molec[i].connect[j][0] = -1;
      molec[i].connect[j][1] = -1;
    }
  } //}}}
  // save ids of bonded beads
  for (int i = 0; i < count_bond_lines; i++) {
    int id1 = index_all[bond[i].index1],
        id2 = index_all[bond[i].index2];
    int m_id = bead_all[id1].Molecule;
    bool done[2] = {false}; // so as not to go through the whole molecule
    for (int j = 0; j < atoms_per_mol[m_id]; j++) {
      if (mol[m_id].Bead[j] == id1) {
        molec[m_id].connect[count_in_mol[m_id]][0] = j;
        done[0] = true;
      }
      if (mol[m_id].Bead[j] == id2) {
        molec[m_id].connect[count_in_mol[m_id]][1] = j;
        done[1] = true;
      }
      if (done[0] && done[1]) {
        break;
      }
    }
    count_in_mol[m_id]++;
  }
  // sort bonds -- pro forma as bead ids in atom struct are sorted
  for (int i = 0; i < (*Counts).Molecules; i++) {
    SortBonds(molec[i].connect, bonds_per_mol[i]);
  }
  //}}}
  /*
  // test print bonds in molecules //{{{
  for (int i = 0; i < (*Counts).Molecules; i++) {
    printf("%d:", i);
    for (int j = 0; j < bonds_per_mol[i]; j++) {
      printf(" %d-%d", connectivity[i][j][0]+1, connectivity[i][j][1]+1);
    }
    putchar('\n');
  }
  putchar('\n'); //}}}
  */
  // count molecule types based //{{{
  /*
   * a) if detailed, check i) through iv)
   * b) if !detailed, check only i) and exit if a molecule contains too few or
   *    too many beads
   */
  MOLECULETYPE *mt = calloc(1, sizeof (MOLECULETYPE));
  for (int i = 0; i < (*Counts).Molecules; i++) {
    // atom struct is necessary because the name is not recorded anywhere else
    int id = bead_all[mol[i].Bead[0]].Index; // id of the molecule's first bead
    int name = atom[atom_id[id]].resname; // name index in res_name array
    bool add = false; // assume no type exists that i can be added to
    // go through all molecule types to find a match for molecule i //{{{
    for (int j = 0; j < (*Counts).TypesOfMolecules; j++) {
      if (detailed) { // a)
        if (strcmp(res_name[name], mt[j].Name) == 0 && // check name,
            atoms_per_mol[i] == mt[j].nBeads && // number of beads, and
            bonds_per_mol[i] == mt[j].nBonds) { // number of bonds
          bool same_beads = true; // assume molecule i has j type's bead order
          for (int k = 0; k < mt[j].nBeads; k++) {
            if (bead_all[mol[i].Bead[k]].Type != mt[j].Bead[k]) {
              same_beads = false; // nope, it doesn't; is not type j
              break;
            }
          }
          bool same_bonds = true; // assume molecule i has j type's connectivity
          for (int k = 0; k < mt[j].nBonds; k++) {
            if (molec[i].connect[k][0] != mt[j].Bond[k][0] ||
                molec[i].connect[k][1] != mt[j].Bond[k][1]) {
              same_bonds = false; // nope, it doesn't; i is not type j
              break;
            }
          }
          if (same_beads && same_bonds) {
            add = true; // add to existing molecule type
          }
        }
      } else { // b)
        if (strcmp(res_name[name], mt[j].Name) == 0) {
          add = true;
          // error - too few/too many beads in a molecule
          if (atoms_per_mol[i] != mt[j].nBeads) {
            ErrorPrintError();
            YellowText(STDERR_FILENO);
            fprintf(stderr, "%s", struct_file);
            RedText(STDERR_FILENO);
            fprintf(stderr, " - molecule ");
            YellowText(STDERR_FILENO);
            fprintf(stderr, "resid %d", i);
            RedText(STDERR_FILENO);
            fprintf(stderr, " contains too ");
            if (atoms_per_mol[i] > mt[j].nBeads) {
              fprintf(stderr, "many ");
            } else {
              fprintf(stderr, "few ");
            }
            fprintf(stderr, "beads (%d instead of %d)\n", atoms_per_mol[i], mt[j].nBeads);
            ResetColour(STDERR_FILENO);
            exit(1);
          }
        }
      }
      if (add) { // found matching molecule type
        mt[j].Number++;
        mol[i].Type = j;
        break;
      }
    } //}}}
    // new molecule type if no match found for i or no type exists yet //{{{
    if (!add || (*Counts).TypesOfMolecules == 0) {
      int type = (*Counts).TypesOfMolecules++;
      mol[i].Type = type;
      mt = realloc(mt, sizeof (MOLECULETYPE) * (*Counts).TypesOfMolecules);
      strcpy(mt[type].Name, res_name[name]);
      mt[type].Number = 1;
      mt[type].nBeads = atoms_per_mol[i];
      mt[type].Bead = malloc(sizeof *mt[type].Bead * mt[type].nBeads);
      for (int j = 0; j < mt[type].nBeads; j++) {
        int bead = mol[i].Bead[j];
        int id = bead;
        mt[type].Bead[j] = bead_all[id].Type;
      }
      mt[type].nBonds = bonds_per_mol[i];
      if (bonds_per_mol[i] > 0) {
        mt[type].Bond = malloc(sizeof *mt[type].Bond * bonds_per_mol[i]);
        for (int j = 0; j < bonds_per_mol[i]; j++) {
          mt[type].Bond[j][0] = molec[i].connect[j][0];
          mt[type].Bond[j][1] = molec[i].connect[j][1];
          mt[type].Bond[j][2] = -1; // no bond types
        }
      }
      mt[type].nAngles = 0;
      mt[type].nDihedrals = 0;
    } //}}}
  } //}}}
  // if detailed, rename the molecule types with the same name //{{{
  if (detailed) {
    for (int i = 0; i < ((*Counts).TypesOfMolecules-1); i++) {
      count = 0; // number of types sharing the name with type i
      for (int j = (i+1); j < (*Counts).TypesOfMolecules; j++) {
        if (strcmp(mt[i].Name, mt[j].Name) == 0) {
          count++;
          // shorten name if necessary to append '_<int>'
          char name[MOL_NAME+1];
          strcpy(name, mt[j].Name);
          if (count < 10) {
            name[MOL_NAME-2] = '\0';
          } else if (count < 100) {
            name[MOL_NAME-3] = '\0';
          } else if (count < 1000) {
            name[MOL_NAME-4] = '\0';
          }
//        P_IGNORE(-Wformat-truncation);
          snprintf(mt[j].Name, MOL_NAME+1, "%s_%d", name, count);
//        P_POP;
        }
      }
    }
  } //}}}
  // calculate molecules' mass and charge and fill their BType array
  FillMolType((*Counts).TypesOfMolecules, bt, &mt); //}}}
  // 6) copy everything back to their 'proper' arrays and structures //{{{
  (*Counts).Beads = (*Counts).BeadsInVsf;
  *Index = malloc(sizeof **Index * (*Counts).Beads);
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Index)[i] = index_all[i];
  }
  *BeadType = malloc(sizeof (BEADTYPE) * (*Counts).TypesOfBeads);
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i] = bt[i];
  }
  *Bead = malloc(sizeof (BEAD) * (*Counts).Beads);
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i] = bead_all[i];
  }
  // molecule types
  CopyMoleculeType((*Counts).TypesOfMolecules, MoleculeType, mt, 2);
  CopyMolecule((*Counts).Molecules, *MoleculeType, Molecule, mol, 2); //}}}
  // free memory //{{{
  free(bt);
  free(index_all);
  FreeBead((*Counts).BeadsInVsf, &bead_all);
  FreeMolecule((*Counts).Molecules, &mol);
  FreeMoleculeType((*Counts).TypesOfMolecules, &mt);
  for (int i = 0; i < (*Counts).Molecules; i++) {
//  for (int j = 0; j < bonds_per_mol[i]; j++) {
//    free(connectivity[i][j]);
//  }
    free(molec[i].connect);
  }
  free(molec);
//free(connectivity);
//for (int i = 0; i < atom_names; i++) {
//  free(atom_name[i]);
//}
  free(atom_name);
  free(atom_id);
//for (int i = 0; i < res_names; i++) {
//  free(res_name[i]);
//}
  free(res_name);
  free(atom);
  free(bond);
  free(atoms_per_mol);
  free(bonds_per_mol);
  free(count_in_mol); //}}}
} //}}}

// TODO: struct_lines no longer relevant
// FullVtfRead() //{{{
/*
 * Function to read vtf structure and (optionally) detect what beads are
 * present in a coordinate file. First, read structure file; second, check
 * coordinate file if present; third TODO: check FIELD/lammps file for extra
 * info (bond/angle stuff). The function warns if the system isn't electrically
 * neutral.
 */
void FullVtfRead(char *struct_file, char *vcf_file, bool detailed, bool vtf,
                 bool *indexed, int *struct_lines, BOX *Box, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // read the whole structure section
  ReadVtfStructure(struct_file, detailed, Counts, BeadType, Bead, Index,
                   MoleculeType, Molecule);
  // check coordinate file if provided //{{{
  if (vcf_file[0] != '\0') {
    VtfReadPBC(vcf_file, Box);
    // number of structure lines (or -1 if coordinate file is not vtf)
//  *struct_lines = VtfCountStructLines(vtf, vcf_file);
    // get timestep type & contained beads from the first timestep
    FILE *vcf;
    if ((vcf = fopen(vcf_file, "r")) == NULL) {
      ErrorFileOpen(vcf_file, 'r');
      exit(1);
    }
//  SkipVtfStructure(vcf, *struct_lines);
    *indexed = CheckVtfTimestep(vcf, vcf_file, Counts, BeadType, Bead, Index,
                                MoleculeType, Molecule);
    fclose(vcf);
  } //}}}
  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, struct_file);
  // allocate memory for aggregates - TODO: that'll be scrapped at some time
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(1, sizeof *(*Bead)[i].Aggregate);
  }
  FillMolMassCharge((*Counts).TypesOfMolecules, MoleculeType, *BeadType);
} //}}}

// ReadVtfTimestepPreamble_old() //{{{
/*
 * Function to read timestep preamble, i.e., comments, pbc, and timestep lines.
 */
int ReadVtfTimestepPreamble_old(bool *indexed, char *input_coor, FILE *vcf_file,
                            char **stuff, VECTOR *BoxLength, bool quit) {
  return 0;
//// save pointer position in vcf_file
//fpos_t position;
//fgetpos(vcf_file, &position); // get file pointer position
//(*stuff)[0] = '\0'; // start with empty string
//int words, count_lines = 0;
//char split[SPL_STR][SPL_LEN];
//bool timestep = false;
//// read lines until float-started line is encountered
//do {
//  char line[LINE] = {0}, line2[LINE];
//  fgets(line, sizeof line, vcf_file);
//  // if the line is too long, skip the rest of it
//  if (strcspn(line, "\n") == (LINE-1)) {
//    while (getc(vcf_file) != '\n')
//      ;
//  }
//  if (feof(vcf_file)) {
//    return -1;
//  }
//  strcpy(line2, line); // save 'unsplit' line
//  words = SplitLine(split, line, " \t");
//  // comment line - copy to stuff array //{{{
//  if (split[0][0] == '#') {
//    strncat(*stuff, line2, LINE-strlen(*stuff)-1);
//  //}}}
//  // error (as long as we care) - incorrect line //{{{
//  /*
//   * error if not
//   * 1) blank line
//   * 2) timestep line starting with with 't[imestep]|i[ndexed]|o[rdered]',
//   * 3) pbc line starting with 'pbc',
//   * 4) float-initiated line, i.e., first coordinate line
//   */
//  } else if (quit && // should the funciton exit on error?
//             words != 0 && // 1)
//             split[0][0] != 't' && //
//             split[0][0] != 'i' && // 2)
//             split[0][0] != 'o' && //
//             strcasecmp(split[0], "pbc") != 0 && // 3)
//             !IsReal(split[0])) { // 4)
//    RedText(STDERR_FILENO);
//    fprintf(stderr, "\nError: ");
//    YellowText(STDERR_FILENO);
//    fprintf(stderr, "%s", input_coor);
//    RedText(STDERR_FILENO);
//    fprintf(stderr, " - unrecognised line in the timestep preamble\n");
//    ResetColour(STDERR_FILENO);
//    ErrorPrintLine(split, words);
//    exit(1);
//  } //}}}
//  // change BoxLength, if correct pbc line is present
//  if (strcasecmp(split[0], "pbc") == 0 &&
//      IsReal(split[1]) && IsReal(split[2]) && IsReal(split[3])) {
//    (*BoxLength).x = atof(split[1]);
//    (*BoxLength).y = atof(split[2]);
//    (*BoxLength).z = atof(split[3]);
//  }
//  // if the line is longer than LINE, copy the next part to stuff array //{{{
//  while (strlen(line2) == (LINE-1) && line2[LINE-1] == '\n') {
//    fgets(line2, sizeof line2, vcf_file);
//    // if the line is too long, skip the rest of it
//    if (strcspn(line, "\n") == (LINE-1)) {
//      while (getc(vcf_file) != '\n')
//        ;
//    }
//    if (feof(vcf_file)) {
//      return -1;
//    }
//    strncat(*stuff, line2, LINE-strlen(*stuff)-1);
//    count_lines++;
//  } //}}}
//  // test for a t(imestep) i(ndexed)/o(rdered) line //{{{
//  if ((words > 1 && split[0][0] == 't' && split[1][0] == 'i') ||
//      split[0][0] == 'i') {
//    timestep = true;
//    *indexed = true; // indexed timestep present
//  } else if ((words > 1 && split[0][0] == 't' && split[1][0] == 'o') ||
//             split[0][0] == 'o' || (words == 1 && split[0][0] == 't')) {
//    timestep = true;
//    *indexed = false; // ordered timestep present
//  } //}}}
//  count_lines++;
//} while (words == 0 || !IsReal(split[0]));
//count_lines--; // the last counted line contained the first coordinate line
//// error (as long as we care) - missing timestep line //{{{
//if (quit && !timestep) {
//  RedText(STDERR_FILENO);
//  fprintf(stderr, "\nError: ");
//  YellowText(STDERR_FILENO);
//  fprintf(stderr, "%s", input_coor);
//  RedText(STDERR_FILENO);
//  fprintf(stderr, " - missing t(imestep) o(rdered)/i(indexed) line\n");
//  if (indexed) {
//    fprintf(stderr, "       or more indexed coordinate lines");
//    fprintf(stderr, " than in the first timestep\n");
//  }
//  ResetColour(STDERR_FILENO);
//  exit(1);
//} //}}}
//fsetpos(vcf_file, &position); // restore pointer position
//return count_lines;
} //}}}

// ReadVtfTimestepPreamble() //{{{
/*
 * Function to read timestep preamble, i.e., comments, pbc, and timestep lines.
 */
int ReadVtfTimestepPreamble(bool *indexed, char *input_coor, FILE *vcf_file,
                            char **stuff, BOX *Box, bool quit) {
  // save pointer position in vcf_file
  fpos_t position;
  fgetpos(vcf_file, &position); // get file pointer position
  (*stuff)[0] = '\0'; // start with empty string
  int words, count_lines = 0;
  char split[SPL_STR][SPL_LEN];
  bool timestep = false;
  // read lines until float-started line is encountered
  do {
    char line[LINE] = {0}, line2[LINE];
    fgets(line, sizeof line, vcf_file);
    // if the line is too long, skip the rest of it
    if (strcspn(line, "\n") == (LINE-1)) {
      while (getc(vcf_file) != '\n')
        ;
    }
    if (feof(vcf_file)) {
      return -1;
    }
    strcpy(line2, line); // save 'unsplit' line
    words = SplitLine(split, line, " \t");
    int ltype = CheckVtLineType(words, split, false, input_coor, count_lines);
    // pbc line - get box dimensions //{{{
    if (ltype == PBC_LINE || ltype == PBC_LINE_ANGLES) {
      (*Box).Length.x = atof(split[1]);
      (*Box).Length.y = atof(split[2]);
      (*Box).Length.z = atof(split[3]);
      (*Box).alpha = 90;
      (*Box).beta = 90;
      (*Box).gamma = 90;
      if (ltype == PBC_LINE_ANGLES) {
        (*Box).alpha = atof(split[4]);
        (*Box).beta = atof(split[5]);
        (*Box).gamma = atof(split[6]);
      } //}}}
    // timestep line - save the timestep type
    } else if (ltype == TIME_LINE) {
      timestep = true;
    //}}}
    // comment line - copy to stuff[] //{{{
    } else if (ltype == COMMENT_LINE) {
      // TODO: thoroughly test & and add some stuff to suppress warnings
      strncat(*stuff, line2, LINE-strlen(*stuff)-1); //}}}
    // error - only if we care \TODO why would we care? //{{{
    } else if (quit && ltype == ERROR_LINE) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - unrecognised line in a timestep preamble\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
    // test for a t(imestep) i(ndexed)/o(rdered) line //{{{
    if ((words > 1 && split[0][0] == 't' && split[1][0] == 'i') ||
        split[0][0] == 'i') {
      timestep = true;
      *indexed = true; // indexed timestep present
    } else if ((words > 1 && split[0][0] == 't' && split[1][0] == 'o') ||
               split[0][0] == 'o' || (words == 1 && split[0][0] == 't')) {
      timestep = true;
      *indexed = false; // ordered timestep present
    } //}}}
    count_lines++;
  } while (words == 0 || !IsReal(split[0]));
  count_lines--; // the last counted line contained the first coordinate line
  // error (as long as we care) - missing timestep line //{{{
  if (quit && !timestep) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing t(imestep) o(rdered)/i(indexed) line\n");
    if (indexed) {
      fprintf(stderr, "       or more indexed coordinate lines");
      fprintf(stderr, " than in the first timestep\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  fsetpos(vcf_file, &position); // restore pointer position
  return count_lines;
} //}}}

// LastStep() //{{{
/*
 * Function to check that there are no more data in vcf/agg files.
 */
bool LastStep(FILE *vcf_file, FILE *agg_file) {
  fpos_t position;
  char split[SPL_STR][SPL_LEN];
  int words = 0;
  // check coordinate file //{{{
  /*
   * the first float-started line is assumed to be the first coordinate line of
   * the new timestep.
   */
  fgetpos(vcf_file, &position); // get coor file pointer position
  do {
    char line[LINE] = {0};
    fgets(line, sizeof line, vcf_file);
    // if the line is too long, skip the rest of it
    if (strcspn(line, "\n") == (LINE-1)) {
      while (getc(vcf_file) != '\n')
        ;
    }
    words = SplitLine(split, line, " \t");
    if (feof(vcf_file)) {
      return true;
    }
  } while (words == 0 || !IsReal(split[0]));
  fsetpos(vcf_file, &position); // restore coor file pointer position //}}}
  // if aggregate file is provided, check that too //{{{
  if (agg_file != NULL) {
    fgetpos(agg_file, &position); // get agg file pointer position
    // check for aggregate file ending
    if (feof(agg_file)) {
      return true;
    }
    fsetpos(agg_file, &position); // restore agg file pointer position
  } //}}}
  return false;
} //}}}

// ReadCoordinates_old() //{{{
/*
 * Function reading coordinates from vtf file with indexed timesteps (\ref
 * IndexedCoorFile). TODO: deprecated
 */
void ReadCoordinates_old(bool indexed, char *input_coor, FILE *vcf_file, COUNTS Counts, int *Index, BEAD **Bead, char **stuff) {
  bool test_indexed; // to test if the present timestep type is the same as detected by ReadStructure()
  // TODO: added box for variable box size
  BOX Box;
  int count_lines = ReadVtfTimestepPreamble(&test_indexed, input_coor, vcf_file,
                                            stuff, &Box, true);
  // error - wrong type of step (indexed vs. ordered) //{{{
  if (test_indexed != indexed) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    if (test_indexed) {
      fprintf(stderr, " - indexed timestep instead of an ordered one\n");
    } else {
      fprintf(stderr, " - ordered timestep instead of an indexed one\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // skip the preamble lines //{{{
  for (int i = 0; i < count_lines; i++) {
    while (getc(vcf_file) != '\n')
      ;
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      // error - end of file
      if (feof(vcf_file)) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unexpected end of file\n");
        fprintf(stderr, "       possibly fewer beads in a timestep than in the first timestep\n");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double>
      if (words < 4 || !IsInteger(split[0]) || !IsReal(split[1]) ||
          !IsReal(split[2]) || !IsReal(split[3])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read a coordinate line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      int index = atoi(split[0]);
      // bead coordinates
      (*Bead)[Index[index]].Position.x = atof(split[1]);
      (*Bead)[Index[index]].Position.y = atof(split[2]);
      (*Bead)[Index[index]].Position.z = atof(split[3]);
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double>
      if (words < 3 || !IsReal(split[0]) ||
          !IsReal(split[1]) || !IsReal(split[2])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read coordinates\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
    }
  } //}}}
} //}}}

// ReadVcfCoordinates_old() //{{{
/*
 * Function reading coordinates from vtf coordinate file.
 */
void ReadVcfCoordinates_old(bool indexed, char *input_coor, FILE *vcf_file,
                        VECTOR *BoxLength, COUNTS Counts,
                        int *Index, BEAD **Bead, char **stuff) {
  // count preamble lines (exiting on an unrecognised line)
  bool test_indexed; // is current timestep ordered or indexed?
  BOX Box;
  int preamble_lines = ReadVtfTimestepPreamble(&test_indexed, input_coor,
                                               vcf_file, stuff, &Box,
                                               true);
  // error - wrong type of step (indexed vs. ordered) //{{{
  if (test_indexed != indexed) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    if (test_indexed) {
      fprintf(stderr, " - indexed timestep instead of an ordered one\n");
    } else {
      fprintf(stderr, " - ordered timestep instead of an indexed one\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // skip the preamble lines //{{{
  for (int i = 0; i < preamble_lines; i++) {
    while (getc(vcf_file) != '\n')
      ;
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      // error - end of file //{{{
      if (feof(vcf_file)) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unexpected end of file\n");
        fprintf(stderr, "       maybe fewer beads in a timestep");
        fprintf(stderr, " than in the first timestep\n");
        ResetColour(STDERR_FILENO);
        exit(1);
      } //}}}
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double> //{{{
      if (words < 4 || !IsInteger(split[0]) ||
          !IsReal(split[1]) || !IsReal(split[2]) || !IsReal(split[3])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read a coordinate line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      int index = atoi(split[0]);
      // bead coordinates
      int id = Index[index];
      (*Bead)[id].Position.x = atof(split[1]);
      (*Bead)[id].Position.y = atof(split[2]);
      (*Bead)[id].Position.z = atof(split[3]);
      if (words >= 7) { // bead velocities, if present
        (*Bead)[id].Velocity.x = atof(split[4]);
        (*Bead)[id].Velocity.y = atof(split[5]);
        (*Bead)[id].Velocity.z = atof(split[6]);
      }
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double> //{{{
      if (words < 3 ||
          !IsReal(split[0]) || !IsReal(split[1]) || !IsReal(split[2])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read coordinates\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
      if (words >= 6) { // bead velocities, if present
        (*Bead)[i].Velocity.x = atof(split[3]);
        (*Bead)[i].Velocity.y = atof(split[4]);
        (*Bead)[i].Velocity.z = atof(split[5]);
      }
    }
  } //}}}
} //}}}

// ReadVcfCoordinates() //{{{
/*
 * Function reading coordinates from vtf coordinate file.
 */
void ReadVcfCoordinates(bool indexed, char *input_coor, FILE *vcf_file,
                        BOX *Box, COUNTS Counts,
                        int *Index, BEAD **Bead, char **stuff) {
  // count preamble lines (exiting on an unrecognised line)
  bool test_indexed; // is current timestep ordered or indexed?
  int preamble_lines = ReadVtfTimestepPreamble(&test_indexed, input_coor,
                                               vcf_file, stuff, Box, true);
  TriclinicCellData(Box);
  // error - wrong type of step (indexed vs. ordered) //{{{
  if (test_indexed != indexed) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    if (test_indexed) {
      fprintf(stderr, " - indexed timestep instead of an ordered one\n");
    } else {
      fprintf(stderr, " - ordered timestep instead of an indexed one\n");
    }
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  // skip the preamble lines //{{{
  for (int i = 0; i < preamble_lines; i++) {
    while (getc(vcf_file) != '\n')
      ;
  } //}}}
  if (indexed) { // indexed timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      // error - end of file //{{{
      if (feof(vcf_file)) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - unexpected end of file\n");
        fprintf(stderr, "       maybe fewer beads in a timestep");
        fprintf(stderr, " than in the first timestep\n");
        ResetColour(STDERR_FILENO);
        exit(1);
      } //}}}
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <int> <double> <double> <double> //{{{
      if (words < 4 || !IsInteger(split[0]) ||
          !IsReal(split[1]) || !IsReal(split[2]) || !IsReal(split[3])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read a coordinate line\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      int index = atoi(split[0]);
      // bead coordinates
      int id = Index[index];
      (*Bead)[id].Position.x = atof(split[1]);
      (*Bead)[id].Position.y = atof(split[2]);
      (*Bead)[id].Position.z = atof(split[3]);
      if (words >= 7) { // bead velocities, if present
        (*Bead)[id].Velocity.x = atof(split[4]);
        (*Bead)[id].Velocity.y = atof(split[5]);
        (*Bead)[id].Velocity.z = atof(split[6]);
      }
    } //}}}
  } else { // ordered timestep //{{{
    for (int i = 0; i < Counts.Beads; i++) {
      char line[LINE], split[SPL_STR][SPL_LEN];
      fgets(line, sizeof line, vcf_file);
      // if the line is too long, skip the rest of it
      if (strcspn(line, "\n") == (LINE-1)) {
        while (getc(vcf_file) != '\n')
          ;
      }
      int words = SplitLine(split, line, " \t");
      // error - coordinate line must be <double> <double> <double> //{{{
      if (words < 3 ||
          !IsReal(split[0]) || !IsReal(split[1]) || !IsReal(split[2])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", input_coor);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - cannot read coordinates\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      // bead coordinates
      (*Bead)[i].Position.x = atof(split[0]);
      (*Bead)[i].Position.y = atof(split[1]);
      (*Bead)[i].Position.z = atof(split[2]);
      if (words >= 6) { // bead velocities, if present
        (*Bead)[i].Velocity.x = atof(split[3]);
        (*Bead)[i].Velocity.y = atof(split[4]);
        (*Bead)[i].Velocity.z = atof(split[5]);
      }
    }
  } //}}}
} //}}}

// SkipVcfCoor() //{{{
/*
 * Function to skip one timestep in a vtf coordinate file.
 */
void SkipVcfCoor(FILE *vcf_file, char *input_coor,
                 COUNTS Counts, char **stuff) {
  bool rubbish; // testing timestep type - not used here
  // skip timestep preamble (don't exit if it includes something wrong)
  BOX dust; // not used
  int preamble_lines = ReadVtfTimestepPreamble(&rubbish, input_coor, vcf_file,
                                               stuff, &dust, false);
  for (int i = 0; i < preamble_lines; i++) {
    while (getc(vcf_file) != '\n')
      ;
  }
  // skip coordinate lines
  for (int i = 0; i < Counts.Beads; i++) {
    while (getc(vcf_file) != '\n')
      ;
    // premature end of file?
    if (feof(vcf_file) == EOF) {
      ErrorPrintError();
      RedText(STDERR_FILENO);
      fprintf(stderr, "premature end of ");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_coor);
      RedText(STDERR_FILENO);
      fprintf(stderr, " file\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
  }
} //}}}

// TODO: rewrite to use fgets & SplitLine instead of fscanf
// TODO: restructure aggregates - don't include in BEAD
// ReadAggregates() //{{{
/*
 * Function reading information about aggregates from agg file generated by
 * Aggregates utility.
 */
void ReadAggregates(FILE *fr, char *agg_file, COUNTS *Counts,
                    AGGREGATE **Aggregate, BEADTYPE *BeadType, BEAD **Bead,
                    MOLECULETYPE *MoleculeType, MOLECULE **Molecule,
                    int *Index) {
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read 'Step|Last Step' line
  fgets(line, sizeof line, fr);
  int words = SplitLine(split, line, " \t");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  if (split[0][0] == 'L') {
    ErrorPrintError();
    RedText(STDERR_FILENO);
    fprintf(stderr, "premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
  } else if (words < 2 || strcmp(split[0], "Step:") != 0 ||
             !IsInteger(split[1])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - wrong 'Step' line\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // initialize array of number of aggregates per bead //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].nAggregates = 0;
  } //}}}
  // get number of aggregates //{{{
  fgets(line, sizeof line, fr);
  words = SplitLine(split, line, " \t");
  // error - the number of aggregates must be <int>
  if (words == 0 || !IsInteger(split[0])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - number of aggregates must be a whole number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  (*Counts).Aggregates = atoi(split[0]); //}}}
  // skip blank line - error if not a blank line //{{{
  fgets(line, sizeof line, fr);
  words = SplitLine(split, line, " \t");
  if (words > 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing a blank line after number of aggregates\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
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
        (*Bead)[id].nAggregates++;
        (*Bead)[id].Aggregate = realloc((*Bead)[id].Aggregate,
                                        sizeof *(*Bead)[id].Aggregate *
                                        (*Bead)[id].nAggregates);
        (*Bead)[id].Aggregate[(*Bead)[id].nAggregates-1] = i;
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
        (*Bead)[bead].Aggregate[(*Bead)[bead].nAggregates] = i;
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
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read (Last) Step line
  fgets(line, sizeof line, agg);
  int words = SplitLine(split, line, " \t");
  // error if the first line is 'L[ast Step]' or isn't 'Step: <int>'//{{{
  if (split[0][0] == 'L') {
    ErrorPrintError();
    RedText(STDERR_FILENO);
    fprintf(stderr, "premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
  } else if (words < 2 || !IsInteger(split[1])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - wrong 'Step' line\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // get number of aggregates
  fgets(line, sizeof line, agg);
  words = SplitLine(split, line, " \t");
  // Error - number of aggregates must be <int> //{{{
  if (words != 0 && !IsInteger(split[0])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - number of aggregates must be a whole number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  int aggs = atoi(split[0]);
  // skip the blank line and the aggregate lines
  for (int i = 0; i < (1+2*aggs); i++) {
    int test;
    while ((test = getc(agg)) != '\n') {
      if (test == EOF) {
        ErrorPrintError();
        RedText(STDERR_FILENO);
        fprintf(stderr, "premature end of ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", agg_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " file");
        ResetColour(STDERR_FILENO);
        exit(1);
      }
    }
  }
  // skip empty line at the end
  fgets(line, sizeof line, agg);
  words = SplitLine(split, line, " \t");
  if (feof(agg) == EOF) {
    ErrorPrintError();
    RedText(STDERR_FILENO);
    fprintf(stderr, "premature end of ");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", agg_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " file");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
} //}}}

// ReadFieldPbc() //{{{
/*
 * Function reading box size from the first line of the FIELD-like file.
 */
bool ReadFieldPbc(char *field, VECTOR *BoxLength) {
  bool pbc = false; // assume the first line doesn't contain box size
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}
  // read first line
  char line[LINE], split[SPL_STR][SPL_LEN];
  fgets(line, sizeof line, fr);
  int words = SplitLine(split, line, " \t");
  /*
   * box size must be 'pbc <double> <double> <double>'; test
   * 1) number of strings
   * 2) 'pbc' string
   * 3) positive doubles
   */
  if (words >= 4 && // 1)
      strcasecmp(split[0], "pbc") == 0 && // 2)
      IsPosReal(split[1]) && //
      IsPosReal(split[2]) && // 3)
      IsPosReal(split[3])) { //
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
                       BEADTYPE **BeadType, BEAD **Bead) {
  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read number of bead types //{{{
  bool missing = true; // assume species keyword is missing
  while(fgets(line, sizeof line, fr)) {
    int words = SplitLine(split, line, " \t");
    if (strcasecmp(split[0], "species") == 0) {
      missing = false; // species keyword is present
      // check if the next string is a number
      if (words < 2 || !IsInteger(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - missing number of species\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[1]);
      break;
    }
  }
  if (missing) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'species' keyword\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof (BEADTYPE));
  // read info about bead types //{{{
  (*Counts).Unbonded = 0;
  (*Counts).Beads = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    fgets(line, sizeof line, fr);
    int words = SplitLine(split, line, " \t");
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
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - no blank lines permitted in the species section\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    } else if (words < 4 ||                  // 2)
               !IsPosReal(split[1]) ||     // 3)
               !IsReal(split[2]) ||        // 4)
               !IsInteger(split[3])) {       // 5)
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - wrong species line\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    } //}}}
//  P_IGNORE(-Wformat-truncation);
    snprintf((*BeadType)[i].Name, BEAD_NAME+1, "%s", split[0]);
//  P_POP;
    (*BeadType)[i].Mass = atof(split[1]);
    (*BeadType)[i].Charge = atof(split[2]);
    (*BeadType)[i].Number = atoi(split[3]);
    (*Counts).Unbonded += (*BeadType)[i].Number;
  } //}}}
  (*Counts).Beads = (*Counts).Unbonded;
  // allocate & fill Bead array //{{{
  *Bead = calloc((*Counts).Unbonded, sizeof (BEAD));
  int count = 0;
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    for (int j = 0; j < (*BeadType)[i].Number; j++) {
      (*Bead)[count].Type = i;
      (*Bead)[count].Molecule = -1;
      (*Bead)[count].Index = count;
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
                        BEADTYPE **BeadType, BEAD **Bead,
                        MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
                        PARAMS **bond_type, PARAMS **angle_type,
                        PARAMS **dihedral_type) {

  // open FIELD-like file //{{{
  FILE *fr;
  if ((fr = fopen(field, "r")) == NULL) {
    ErrorFileOpen(field, 'r');
    exit(1);
  } //}}}
  char line[LINE], split[SPL_STR][SPL_LEN];
  // read number of molecule types //{{{
  bool missing = true; // assume molecule keyword is missing
  while(fgets(line, sizeof line, fr)) {
    int words = SplitLine(split, line, " \t");
    if (strncasecmp(split[0], "molecule", 8) == 0) {
      missing = false; // molecule keyword is present
      // error - 'molecule' not followed by an integer
      if (!IsInteger(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - missing number of molecule types\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfMolecules = atoi(split[1]);
      break;
    }
  }
  // error - no molecule keyword
  if (missing) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'molecule' line\n\n");
    ResetColour(STDERR_FILENO);
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
    SplitLine(split, line, " \t");
    if (strcmp(split[0], "finish") == 0) {
      count++;
    }
  }
  // error - fewer 'finish'es than molecule types
  if (count < (*Counts).TypesOfMolecules) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'finish' keyword\n\n");
    ResetColour(STDERR_FILENO);
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
    SplitLine(split, line, " \t");
    if (split[0][0] == '\0') {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - expecting molecule name\n\n");
      ResetColour(STDERR_FILENO);
      exit(1);
    }
    // copy name to MOLECULETYPE
//  P_IGNORE(-Wformat-truncation);
    // MOL_NAME is max string length, i.e., array is longer
    snprintf((*MoleculeType)[i].Name, MOL_NAME+1, "%s", split[0]);
//  P_POP; //}}}
    // 2) //{{{
    fgets(line, sizeof line, fr);
    int words = SplitLine(split, line, " \t");
    if (strncasecmp(split[0], "nummols", 6) != 0 ||
        words < 2 || !IsInteger(split[1])) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - expecting 'nummol[s] <number of molecules>'\n");
      ResetColour(STDERR_FILENO);
      ErrorPrintLine(split, words);
      exit(1);
    }
    (*MoleculeType)[i].Number = atoi(split[1]); //}}}
    // 3) //{{{
    // read number of beads
    fgets(line, sizeof line, fr);
    words = SplitLine(split, line, " \t");
    // error - wrong keyword line //{{{
    if (strncasecmp(split[0], "beads", 4) != 0 ||
        words < 2 || !IsInteger(split[1])) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - expecting 'bead[s] <number of beads>'\n");
      ResetColour(STDERR_FILENO);
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
      words = SplitLine(split, line, " \t");
      // error - bead line must be '<name> <double> <double> <double>' //{{{
      if (words < 4 ||
          !IsReal(split[1]) || !IsReal(split[2]) || !IsReal(split[3])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - wrong bead line in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      } //}}}
      int type = FindBeadType(split[0], *Counts, *BeadType);
      // error - unknown bead type //{{{
      if (type == -1) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - non-existent bead name ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", split[0]);
        RedText(STDERR_FILENO);
        fprintf(stderr, " in molecule ");
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
        ResetColour(STDERR_FILENO);
        ErrorBeadType(*Counts, *BeadType);
        exit(1);
      } //}}}
      (*MoleculeType)[i].Bead[j] = type;
      // TODO: Read coordinates?
    } //}}}
    // 4) //{{{
    // read number of bonds in the molecule
    fgets(line, sizeof line, fr);
    words = SplitLine(split, line, " \t");
    // are bonds present?
    if (strncasecmp(split[0], "bonds", 4) == 0) {
      // error - wrong keyword line //{{{
      if (strncasecmp(split[0], "bonds", 4) != 0 ||
          words < 2 || !IsInteger(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - expecting 'bond[s] <number of bonds>'\n");
        ResetColour(STDERR_FILENO);
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
        words = SplitLine(split, line, " \t");
        // error - bead line must be '<name> <id1> <id2> <double> <double>' //{{{
        if (words < 5 || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsPosReal(split[3]) || !IsPosReal(split[4])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", field);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - wrong bond line in molecule ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ResetColour(STDERR_FILENO);
          ErrorPrintLine(split, words);
          exit(1);
        } //}}}
        // error - wrong bead index //{{{
        if (atoi(split[1]) > (*MoleculeType)[i].nBeads || atoi(split[1]) < 1 ||
            atoi(split[2]) > (*MoleculeType)[i].nBeads || atoi(split[2]) < 1) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", field);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - wrong bead index in a bond in molecule ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s\n", (*MoleculeType)[i].Name);
          ResetColour(STDERR_FILENO);
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
    words = SplitLine(split, line, " \t");
    // are angles present?
    if (strncasecmp(split[0], "angles", 5) == 0) {
      // error - missing number of angles //{{{
      if (words < 2 || !IsInteger(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - expecting 'angle[s] <number of angles>'\n");
        ResetColour(STDERR_FILENO);
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
        words = SplitLine(split, line, " \t");
        // error - wrong angle line //{{{
        // '<type> 3x<id> <double> <double>'
        if (words < 6 || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3]) ||
            !IsPosReal(split[4]) || !IsPosReal(split[5])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", field);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - wrong angle line in molecule ");
          YellowText(STDERR_FILENO);
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
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - extra bond line in molecule ");
      YellowText(STDERR_FILENO);
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
    words = SplitLine(split, line, " \t");
    // are dihedrals present?
    if (strncasecmp(split[0], "dihedrals", 5) == 0) {
      // get number of dihedrals
      // error - wrong number of dihedrals //{{{
      if (words < 2 || !IsInteger(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", field);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - expecting 'dihed[rals] <number of angles>'\n");
        ResetColour(STDERR_FILENO);
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
        words = SplitLine(split, line, " \t");
        // error - wrong dihedral line //{{{
        // '<type> 4x<id> <double> <double>'
        if (words < 7 || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsInteger(split[3]) || !IsInteger(split[4]) ||
            !IsPosReal(split[5]) || !IsPosReal(split[6])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", field);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - wrong dihedral line in molecule ");
          YellowText(STDERR_FILENO);
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
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", field);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - extra bond or angle line in molecule ");
      YellowText(STDERR_FILENO);
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
      SplitLine(split, line, " \t");
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
  (*Counts).Beads += (*Counts).Bonded;
  (*Counts).BeadsInVsf = (*Counts).Beads;
  *Bead = realloc(*Bead, sizeof (BEAD) * (*Counts).Beads);
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
      SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
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
      SplitLine(split, line, " \t");
      for (int k = 0; k < (*Counts).TypesOfAngles; k++) {
        if ((*angle_type)[k].a == atof(split[4]) && (*angle_type)[k].b == atof(split[5])) {
          (*MoleculeType)[i].Angle[j][3] = k;
          break;
        }
      }
    } //}}}
    // skip till 'finish' //{{{
    while(fgets(line, sizeof line, fr)) {
      SplitLine(split, line, " \t");
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
               BEADTYPE **BeadType, BEAD **Bead, int **Index,
               MOLECULETYPE **MoleculeType, MOLECULE **Molecule,
               PARAMS **bond_type, PARAMS **angle_type,
               PARAMS **dihedral_type) {

  // read pbc if required //{{{
  if (BoxLength != NULL && !ReadFieldPbc(field, BoxLength)) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", field);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - first line must start with box size, ");
    fprintf(stderr, "i.e., 'pbc <double> <double> <double>'\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  ReadFieldBeadType(field, Counts, BeadType, Bead);
  // FIELD doesn't contain bead radius, so fill it with impossible values
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    (*BeadType)[i].Radius = RADIUS;
  }
  ReadFieldMolecules(field, Counts, BeadType, Bead, MoleculeType, Molecule,
                     bond_type, angle_type, dihedral_type);
  // allocate Bead[].Aggregate array - needed only to free() //{{{
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = calloc(10, sizeof *(*Bead)[i].Aggregate);
  } //}}}
  FillMolBTypes((*Counts).TypesOfMolecules, MoleculeType);
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
        (*Bead)[count_bead].Index = count_bead;
        count_bead++;
      }
      count_mol++;
    }
  } //}}}
  // fill Index - I don't thinks it's used right now //{{{
  *Index = malloc(sizeof **Index * (*Counts).Beads);
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Index)[i] = i;
  } //}}}
  // check electroneutrality
  WarnElNeutrality(*Counts, *BeadType, field);
} //}}}

// ReadLmpData() //{{{
void ReadLmpData(char *data_file, int *bonds, PARAMS **bond_type,
                 int *angles, PARAMS **angle_type,
                 VECTOR *BoxLength, VECTOR *box_lo, COUNTS *Counts,
                 BEADTYPE **BeadType, BEAD **Bead, int **Index,
                 MOLECULETYPE **MoleculeType, MOLECULE **Molecule) {
  // open data file //{{{
  FILE *fr;
  if ((fr = fopen(data_file, "r")) == NULL) {
    ErrorFileOpen(data_file, 'r');
    exit(1);
  } //}}}

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
    words = SplitLine(split, line, " \t");
    // number of atoms //{{{
    if (words > 1 && strcmp(split[1], "atoms") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'atoms' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).BeadsInVsf = atoi(split[0]);
      (*Counts).Beads = atoi(split[0]);
    } //}}}
    // number of bonds //{{{
    if (words > 1 && strcmp(split[1], "bonds") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'bonds' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *bonds = atoi(split[0]);
    } //}}}
    // number of angles //{{{
    if (words > 1 && strcmp(split[1], "angles") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'angles' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      *angles = atoi(split[0]);
    } //}}}
    // number of bead types //{{{
    if (words > 2 && strcmp(split[1], "atom") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'atom types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBeads = atoi(split[0]);
    } //}}}
    // number of bond types //{{{
    if (words > 2 && strcmp(split[1], "bond") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'bond types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfBonds = atoi(split[0]);
    } //}}}
    // number of angle types //{{{
    if (words > 2 && strcmp(split[1], "angle") == 0 && strcmp(split[2], "types") == 0) {
      if (!IsInteger(split[0])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'angle types' keyword must be preceded by integer\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*Counts).TypesOfAngles = atoi(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "xlo") == 0 && strcmp(split[3], "xhi") == 0) {
      if (!IsReal(split[0]) || !IsReal(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'xlo xhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).x = atof(split[1]) - atof(split[0]);
      (*box_lo).x = atof(split[0]);
    } //}}}
    // box length in y //{{{
    if (words > 3 && strcmp(split[2], "ylo") == 0 && strcmp(split[3], "yhi") == 0) {
      if (!IsReal(split[0]) || !IsReal(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'ylo yhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).y = atof(split[1]) - atof(split[0]);
      (*box_lo).y = atof(split[0]);
    } //}}}
    // box length in x //{{{
    if (words > 3 && strcmp(split[2], "zlo") == 0 && strcmp(split[3], "zhi") == 0) {
      if (!IsReal(split[0]) || !IsReal(split[1])) {
        ErrorPrintError();
        YellowText(STDERR_FILENO);
        fprintf(stderr, "%s", data_file);
        RedText(STDERR_FILENO);
        fprintf(stderr, " - 'zlo zhi' keyword must be preceded by two floats\n");
        ResetColour(STDERR_FILENO);
        ErrorPrintLine(split, words);
        exit(1);
      }
      (*BoxLength).z = atof(split[1]) - atof(split[0]);
      (*box_lo).z = atof(split[0]);
    } //}}}
  } while (words == 0 ||
           split[0][0] == '#' ||
           IsReal(split[0]) ||
           IsInteger(split[0])); //}}}

  // some error checking //{{{
  if ((*Counts).TypesOfBeads == 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'atom types' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*Counts).BeadsInVsf == 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'atoms' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).x == 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'xlo xhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).y == 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'ylo yhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  }
  if ((*BoxLength).z == 0) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", data_file);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - missing 'zlo zhi' line (or is 0)\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}

  // fill something in BeadType struct //{{{
  *BeadType = calloc((*Counts).TypesOfBeads, sizeof (BEADTYPE));
  for (int i = 0; i < (*Counts).TypesOfBeads; i++) {
    sprintf((*BeadType)[i].Name, "bead%d", i+1);
    (*BeadType)[i].Use = true;
    (*BeadType)[i].Write = true;
  } //}}}

  // bead struct memory allocation //{{{
  *Bead = calloc((*Counts).Beads, sizeof (BEAD));
  for (int i = 0; i < (*Counts).Beads; i++) {
    (*Bead)[i].Aggregate = malloc(sizeof *(*Bead)[i].Aggregate * 1);
  } //}}}

  *Index = calloc((*Counts).Beads, sizeof **Index);
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
        words = SplitLine(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 2 || !IsInteger(split[0]) || !IsPosReal(split[1])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Masses' section must start with '<int> <float>'\n");
          ResetColour(STDERR_FILENO);
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
        words = SplitLine(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosReal(split[1]) || !IsPosReal(split[2])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Bond Coeffs' section must start with '<int> <float> <float>'\n");
          ResetColour(STDERR_FILENO);
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
        words = SplitLine(split, line, " \t");
        // error if incorrect line //{{{
        if (words < 3 || !IsInteger(split[0]) || !IsPosReal(split[1]) || !IsPosReal(split[2])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each line in 'Angle Coeffs' section must start with '<int> <float> <float>'\n");
          ResetColour(STDERR_FILENO);
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
      mols = calloc((*Counts).Beads, sizeof *mols);
      // array for list of beads in each molecule
      // bead_mols[i] ... molecule's id; bead_mols[][i] ... molecule's beads
      int **bead_mols = calloc((*Counts).Beads, sizeof **bead_mols);
      for (int i = 0; i < (*Counts).Beads; i++) {
        bead_mols[i] = malloc(sizeof *bead_mols[i] * 1);
      }
      // skip one line (mandatory in lammps data format)
      fgets(line, sizeof line, fr);
      // go through atom section to get basic info //{{{
      fpos_t pos; // set file counter
      fgetpos(fr, &pos); // save file pointer
      for (int i = 0; i < (*Counts).Beads; i++) {
        fgets(line, sizeof line, fr);
        words = SplitLine(split, line, " \t");
        // format of each line: <id> <mol_id> <btype> <charge> <x> <y> <z>
        // Error - incorrect format //{{{
        if (words < 7 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) || !IsInteger(split[2]) ||
            !IsReal(split[3]) || !IsReal(split[4]) || !IsReal(split[5]) || !IsReal(split[6])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Atoms' line must be <id> <mol_id> <bead type> <charge> <x> <y> <z>\n");
          ResetColour(STDERR_FILENO);
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
        (*Bead)[id].Index = id;
        (*Index)[id] = id;
      } //}}}
      (*Counts).Unbonded = monomer;
      // go through possible molecules and remove 1-bead molecules //{{{
      int count = 0; // count real molecules (i.e., those with more than 1 bead)
      for (int i = 0; i < (*Counts).Beads; i++) {
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
      for (int i = 0; i < (*Counts).Beads; i++) {
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
      for (int i = (*Counts).Molecules; i < (*Counts).Beads; i++) {
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
      for (int i = 0; i < (*Counts).Beads; i++) {
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
        words = SplitLine(split, line, " \t");
        // format of each line: <bond id> <bond type> <bead1> <bead2>
        // Error - incorrect format //{{{
        if (words < 4 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Bonds' line must be <bond id> <bond type> <bead1d> <bead2>\n");
          ResetColour(STDERR_FILENO);
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
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - beads in ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "bond %d", atoi(split[0]));
          RedText(STDERR_FILENO);
          fprintf(stderr, "are in different molecules\n\n");
          ResetColour(STDERR_FILENO);
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
        int lowest = (*Counts).Beads; // just some high number
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
        words = SplitLine(split, line, " \t");
        // format of each line: <angle id> <angle type> <bead1> <bead2> <bead3>
        // Error - incorrect format //{{{
        if (words < 5 ||
            !IsInteger(split[0]) || !IsInteger(split[1]) ||
            !IsInteger(split[2]) || !IsInteger(split[3]) || !IsInteger(split[4])) {
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - each 'Angles' line must be <angle id> <angle type> <bead1d> <bead2> <bead3>\n");
          ResetColour(STDERR_FILENO);
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
          ErrorPrintError();
          YellowText(STDERR_FILENO);
          fprintf(stderr, "%s", data_file);
          RedText(STDERR_FILENO);
          fprintf(stderr, " - atoms in ");
          YellowText(STDERR_FILENO);
          fprintf(stderr, "angle %d", atoi(split[0]));
          RedText(STDERR_FILENO);
          fprintf(stderr, " are in different molecules\n\n");
          ResetColour(STDERR_FILENO);
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
        int lowest = (*Counts).Beads; // just some high number
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
    words = SplitLine(split, line, " \t");
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

// SkipCoorSteps() { //{{{
/*
 * Function to skip timesteps (coordinate file only) from the beginning (-st
 * option).
 */
int SkipCoorSteps(FILE *vcf, char *input_coor, COUNTS Counts, int start, bool silent) {
  int test;
  int count = 0;
  char *stuff = malloc(sizeof *stuff * LINE); // just for SkipVcfCoor
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}
    SkipVcfCoor(vcf, input_coor, Counts, &stuff);
  }
  // print starting step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Starting step: %d\n", start);
  } //}}}
  free(stuff);
  // error if last step //{{{
  if (LastStep(vcf, NULL)) {
    fflush(stdout);
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - starting timestep (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", start);
    RedText(STDERR_FILENO);
    fprintf(stderr, ") is higher than the total number of steps (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", count);
    RedText(STDERR_FILENO);
    fprintf(stderr, ")\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  return count;
} //}}}

// SkipCoorAggSteps() { //{{{
/*
 * Function to skip timesteps (coordinate and aggregate files) from the
 * beginning (-st option).
 */
int SkipCoorAggSteps(FILE *vcf, char *input_coor, FILE *agg, char *input_agg,
                     COUNTS Counts, int start, bool silent) {
  int test, count = 0;
  char *stuff = malloc(sizeof *stuff * LINE); // just for SkipVcfCoor
  // TODO: use feof()
  for (int i = 1; i < start && (test = getc(vcf)) != EOF; i++) {
    ungetc(test, vcf);
    count++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rDiscarding step: %d", count);
    } //}}}
    SkipAgg(agg, input_agg);
    SkipVcfCoor(vcf, input_coor, Counts, &stuff);
  }
  free(stuff);
  // print starting step? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Starting step: %d\n", start);
  } //}}}
  // error if last step //{{{
  if (LastStep(vcf, NULL)) {
    fflush(stdout);
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - starting timestep (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", start);
    RedText(STDERR_FILENO);
    fprintf(stderr, ") is higher than the total number of steps (");
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%d", count);
    RedText(STDERR_FILENO);
    fprintf(stderr, ")\n\n");
    ResetColour(STDERR_FILENO);
    exit(1);
  } //}}}
  return count;
} //}}}

// CheckVtfTimestepLine() //{{{
/*
 * Function to check if the provided line is a timestep line.
 *
 */
bool CheckVtfTimestepLine(int words, char split[SPL_STR][SPL_LEN]) {
  // there are several possibilities how the timestep line can look
  // 1) 't[imestep]' - ordered
  if (words == 1 && strncasecmp(split[0], "timestep", 1) == 0) {
    return true;
  }
  // 2) 't[imestep] i[indexed]/o[rdered] any extra stuff'
  if (words > 1 && split[0][0] == 't' &&
      (split[1][0] == 'i' ||
       split[1][0] == 'o')) {
    return true;
  }
  // 3) 'i[indexed]/o[rdered] any extra stuff'
  if (words > 0 &&
      (split[0][0] == 'i' ||
       split[0][0] == 'o')) {
    return true;
  }
  return false; // not a timestep line
} //}}}

// CheckVtfPbcLine() //{{{
/*
 * Function to check if the provided line is a pbc line. It returns -1 if not
 * and 0 or 1 if pbc line without or with angles, respectively.
 */
int CheckVtfPbcLine(int words, char split[SPL_STR][SPL_LEN], char *file,
                    int line_count) {
  // valid line: pbc <x> <y> <z> [<alpha> <beta> <gamm>]
  // invalid line
  if (words < 4 || strcmp(split[0], "pbc") != 0 || !IsPosReal(split[1]) ||
                                                   !IsPosReal(split[2]) ||
                                                   !IsPosReal(split[3])) {
    return ERROR_LINE;
  }
  if (words > 6 && IsPosReal(split[4]) &&
                   IsPosReal(split[5]) &&
                   IsPosReal(split[6])) {
    return PBC_LINE_ANGLES;
  } else {
    return PBC_LINE;
  }
} //}}}

// CheckVtfAtomLine_old() //{{{
/*
 * Function to check if the provided line is a proper vtf structure atom line.
 */
bool CheckVtfAtomLine_old(int words, char split[SPL_STR][SPL_LEN], char *error) {
  // is there even number of strings?
  if ((words%2) != 0) {
    strcpy(error, "odd number of strings in an atom line");
    return false;
  }
  // is the line starting with a[tom] <id>/default?
  if (strcmp(split[1], "default") != 0 && !IsInteger(split[1])) {
    strcpy(error, "'atom' must be followed by <int> or 'default'");
    return false;
  }
  // is the mandatory n[ame] keyword present?
  bool exists = false;
  for (int i = 2; i < words; i+=2) {
    if (split[i][0] == 'n') {
      exists = true;
      break;
    }
  }
  if (!exists) {
    strcpy(error, "missing 'name' in atom line");
    return false;
  }
  // if res[name] is present, there must be resid as well
  for (int i = 2; i < words; i+=2) {
    if (strncmp(split[i], "res", 3) == 0) {
      if (strcmp(split[i], "resid") == 0) { // search for resname
        exists = false;
        for (int j = 2; j < words; j+=2) {
          if (strncmp(split[j], "resname", 3) == 0 &&
              strcmp(split[j], "resid") != 0) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          strcpy(error, "if 'resid' is present, 'resname' must be too");
          return false;
        }
        break;
      } else { // search for resid
        exists = false;
        for (int j = 2; j < words; j+=2) {
          if (strcmp(split[j], "resid") == 0) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          strcpy(error, "if 'resname' is present, 'resid' must be too");
          return false;
        }
        break;
      }
    }
  }
  // if mass, charge, or radius are present, they must be followed by <float>
  // if resid is present, it must be followed by an <int>
  for (int i = 2; i < words; i+=2) {
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') &&
        !IsReal(split[i+1])) {
      strcpy(error, "'charge|q' must be followed by <float>");
      return false;
    } else if (split[i][0] == 'r' && strncmp(split[i], "res", 3) != 0 &&
               (!IsPosReal(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "'radius' must be followed by positive <float>");
      return false;
    } else if (split[i][0] == 'm' &&
               (!IsPosReal(split[i+1]) || atof(split[i+1]) <= 0)) {
      strcpy(error, "'mass' must be followed by positive <float>");
      return false;
    } else if (strcmp(split[i], "resid") == 0 && !IsInteger(split[i+1])) {
      strcpy(error, "'resid' must be followed by a whole number");
      return false;
    }
  }
  return true;
} //}}}

// CheckVtfAtomLine() //{{{
/*
 * Function to check if the provided line is a proper vtf structure atom line.
 */
bool CheckVtfAtomLine(int words, char split[SPL_STR][SPL_LEN]) {
  // error - odd number of strings
  if ((words%2) != 0) {
    return false;
  }
  // error - line not starting with a[tom] default/<id>
  if (split[0][0] != 'a' ||
      (strcmp(split[1], "default") != 0 && !IsInteger(split[1]))) {
    return false;
  }
  // check <keyword> <value> pairs
  bool name = false, resid = false, resname = false;
  for (int i = 2; i < words; i+=2) {
    // is n[ame] keyword present?
    if (split[i][0] == 'n') {
      name = true;
    }
    int r_id = strcmp(split[i], "resid"); // resid cannot be shortened
    int r_name = strncmp(split[i], "res", 3); // res[name] can be shortened
    // is resid keyword present?
    if (r_id == 0) {
      resid = true;
      // resid must be followed by positive integer
      // TODO check when IsPosInteger/IsNatural exists
      if (!IsInteger(split[i+1]) || atoi(split[i+1]) > 0) {
        return false;
      }
    }
    // is res[name] keyword present?
    if (r_id != 0 && r_name == 0) {
      resname = true;
    }
    // error - charge|q
    if ((strcmp(split[i], "charge") == 0 || split[i][0] == 'q') &&
        !IsReal(split[i+1])) {
      return false;
    // error - r[adius] not followed by positive number
    } else if (split[i][0] == 'r' && r_name != 0 &&
               !IsPosReal(split[i+1])) {
      return false;
    // error - m[ass] not followed by positive number
    } else if (split[i][0] == 'm' && !IsPosReal(split[i+1])) {
      return false;
    }
  }
  // error - missing the mandatory n[ame] keyword or
  //         if res[name] is present, there must be resid as well
  if (!name || (!resid && resname) || (resid && !resname)) {
    return false;
  }
  // valid atom line
  return true;
} //}}}

// CheckVtfBondLine_old() //{{{
/*
 * Function to check if the provided line is a proper vtf structure bond line.
 */
bool CheckVtfBondLine_old(int words, char split[SPL_STR][SPL_LEN], char *error) {
  // the line must be: b[ond] <id>:<id> (with possible blanks after ':')
  switch(words) {
    case 1: // there must be at least 1 more string after 'bond'
      strcpy(error, "missing a string after 'bond'");
      return false;
    case 2: // if there's only 1 string after 'bond', it must be '<int>:<int>'
      for (int i = 1; i < (strlen(split[1])-2); i++) {
        if ((split[1][i] == ':' && split[1][i+1] == ':') ||
            split[1][i] == ',' || split[1][i+1] == ',') {
          strcpy(error, "for now, only bond lines in the format <int>:<int> \
are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
          return false;
        }
      }
      char index[SPL_STR][SPL_LEN], string[SPL_LEN];
      strcpy(string, split[1]);
      // split '<int>:<int>' into two strings and test those
      int strings = SplitLine(index, string, ":");
      if (strings < 2 || !IsInteger(index[0]) || !IsInteger(index[1]) ||
          split[1][0] == ':' || split[1][strlen(split[1])-1] == ':') {
        strcpy(error, "missing two colon separated integers in a bond line");
        return false;
      }
      return true;
    case 3: // if there are 2 strings after 'bond', they must be '<int>: <int>'
      if (split[1][strlen(split[1])-1] == ':' && IsInteger(split[2])) {
        if (split[1][strlen(split[1])-2] == ':') { // two colons aren't allowed
          strcpy(error, "for now, only bond lines in the format <int>:<int> \
are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
          return false;
        }
        // test that <int> preceeds the colon
        strcpy(string, split[1]);
        string[strlen(string)-1] = '\0';
        if (!IsInteger(string)) {
          strcpy(error, "missing an integer before colon in a bond line");
          return false;
        }
      } else {
        strcpy(error, "unrecognised bond line (note that whitespace \
is not allowed before the colon; i.e., '<int> : <int>' is not allowed)");
        return false;
      }
      return true;
    default: // more then 2 strings after 'bond' isn't allowed (for now)
      strcpy(error, "unrecognised bond line (note that for now, only \
the format <int>:<int> is allowed; i.e., no comma separated list \
or 'two-colon' string of bonds)");
      return false;
  }
} //}}}

// CheckVtfBondLine() //{{{
/*
 * Function to check if the provided line is a proper vtf structure bond line.
 */
bool CheckVtfBondLine(int words, char split[SPL_STR][SPL_LEN]) {
  // valid line b[ond] '<id>:[  ]<id> anything'
  // error - only one string
  if (words < 2 || split[0][0] != 'b') {
    return false;
  }
  // two strings - assume '<int>:<int>' and test it
  if (words == 2) {
    // ':'-split the second string
    char index[SPL_STR][SPL_LEN], string[SPL_LEN];
    strcpy(string, split[1]);
    int strings = SplitLine(index, string, ":");
    if (strings != 2 || !IsInteger(index[0]) || !IsInteger(index[1])) {
      return false;
    }
  }
  return true;
//    char index[SPL_STR][SPL_LEN], string[SPL_LEN];
//    strcpy(string, split[1]);
//    // split '<int>:<int>' into two strings and test those
//    int strings = SplitLine(index, string, ":");
//    if (strings < 2 || !IsInteger(index[0]) || !IsInteger(index[1]) ||
//        split[1][0] == ':' || split[1][strlen(split[1])-1] == ':') {
//      strcpy(error, "missing two colon separated integers in a bond line");
//      return false;
//    }
//    return true;
//  case 3: // if there are 2 strings after 'bond', they must be '<int>: <int>'
//    if (split[1][strlen(split[1])-1] == ':' && IsInteger(split[2])) {
//      if (split[1][strlen(split[1])-2] == ':') { // two colons aren't allowed
//        strcpy(error, "for now, only bond lines in the format <int>:<int>
//are allowed (i.e., no comma separated list or 'two-colon' string of bonds)");
//        return false;
//      }
//      // test that <int> preceeds the colon
//      strcpy(string, split[1]);
//      string[strlen(string)-1] = '\0';
//      if (!IsInteger(string)) {
//        strcpy(error, "missing an integer before colon in a bond line");
//        return false;
//      }
//    } else {
//      strcpy(error, "unrecognised bond line (note that whitespace
//is not allowed before the colon; i.e., '<int> : <int>' is not allowed)");
//      return false;
//    }
//    return true;
//  default: // more then 2 strings after 'bond' isn't allowed (for now)
//    strcpy(error, "unrecognised bond line (note that for now, only
//the format <int>:<int> is allowed; i.e., no comma separated list
//or 'two-colon' string of bonds)");
//    return false;
//}
} //}}}

// CheckVtfCoordinateLine_old() //{{{
/*
 * Function to check if the provided line is a proper vtf coordinate line.
 */
bool CheckVtfCoordinateLine_old(int words, char split[SPL_STR][SPL_LEN],
                                char *error, bool indexed) {
// the line must be: <x> <y> <z>, preceded by <int> in case of indexed timestep
  // not enough strings
  if ((words < 3 && !indexed) ||
      (words < 4 && indexed)) {
    strcpy(error, "wrong coordinate line");
    return false;
  }
  // incorrect string types
  if (indexed && (!IsInteger(split[0]) ||
                  !IsReal(split[1]) ||
                  !IsReal(split[2]) ||
                  !IsReal(split[3]))) {
    strcpy(error, "wrong coordinate line for an indexed timestep");
    return false;
  }
  if (!indexed && (!IsReal(split[0]) ||
                   !IsReal(split[1]) ||
                   !IsReal(split[2]))) {
    strcpy(error, "wrong coordinate line for an ordered timestep");
    return false;
  }
  // correct line
  return true;
} //}}}

// CheckVtfCoordinateLine() //{{{
/*
 * Function to check if the provided line is a proper vtf coordinate line. It
 * returns -1 if not and 0 or 1 if line for indexed or ordered timestep.
 */
int CheckVtfCoordinateLine(int words, char split[SPL_STR][SPL_LEN]) {
  // valid line: [<id>] <x> <y> <z> ... with <id>, it is indexed timestep
  // indexed timestep
  if (words >= 4 && IsInteger(split[0]) &&
      IsReal(split[1]) && IsReal(split[2]) && IsReal(split[3])) {
    return COOR_LINE_I;
  }
  // ordered timestep
  if (words >= 3 && IsReal(split[0]) && IsReal(split[1]) && IsReal(split[2])) {
    return COOR_LINE_O;
  }
  // non-coordinate line
  return -1;
} //}}}

// CheckVtLine() //{{{
/*
 * Function to check what line from a vtf file is passed to it.
 */
int CheckVtLineType(int words, char split[SPL_STR][SPL_LEN], bool indexed,
                    char *file, int line_count) {
  // blank line
  if (words == 0) {
    return BLANK_LINE;
  }
  // comment line
  if (split[0][0] == '#') {
    return COMMENT_LINE;
  }
  // pbc line
  int test = CheckVtfPbcLine(words, split, file, line_count);
  if (test != ERROR_LINE) {
    return test;
  }
  // coordinate line
  test = CheckVtfCoordinateLine(words, split);
  if (test != ERROR_LINE) {
    return test;
  }
  // atom line (vsf)
  if (CheckVtfAtomLine(words, split)) {
    return ATOM_LINE;
  }
  // bond line (vsf)
  if (CheckVtfBondLine(words, split)) {
    return BOND_LINE;
  }
  // timestep line (vcf)
  if (CheckVtfTimestepLine(words, split)) {
    return TIME_LINE;
  }
  return ERROR_LINE;
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
  FILE *agg;
  if ((agg = fopen(input_agg, "r")) == NULL) {
    ErrorFileOpen(input_agg, 'r');
    exit(1);
  }
  // read first line (Aggregate command)
  char line[LINE], split[SPL_STR][SPL_LEN];
  fgets(line, sizeof line, agg);
  // if the line is too long (which it should never be), skip the rest of it
  if (strcspn(line, "\n") == (LINE-1)) {
    while (getc(agg) != '\n')
      ;
  }
  int words = SplitLine(split, line, " \t");
  // error - not enough strings for a proper Aggregate command //{{{
  if (words < 6) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - first line must contain a valid Aggregates command\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  } //}}}
  // read <distance> argument from Aggregates command //{{{
  if (!IsPosReal(split[2])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <distance> in the Aggregate command \
must be non-negative real number\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *distance = atof(split[2]); //}}}
  // read <contacts> argument from Aggregates command //{{{
  if (!IsInteger(split[3])) {
    ErrorPrintError();
    YellowText(STDERR_FILENO);
    fprintf(stderr, "%s", input_agg);
    RedText(STDERR_FILENO);
    fprintf(stderr, " - <contacts> in the Aggregate command \
must be a non-negative integer\n");
    ResetColour(STDERR_FILENO);
    ErrorPrintLine(split, words);
    exit(1);
  }
  *contacts = atof(split[3]); //}}}
  // warning - differently named vcf file than the one in agg file //{{{
  if (strcmp(split[1], input_coor) != 0) {
    YellowText(STDERR_FILENO);
    fprintf(stderr, "\nWarning: coordinate file ");
    CyanText(STDERR_FILENO);
    fprintf(stderr, "%s", input_coor);
    YellowText(STDERR_FILENO);
    fprintf(stderr, " is different to the one in the aggregate file (");
    CyanText(STDERR_FILENO);
    fprintf(stderr, "%s", split[1]);
    YellowText(STDERR_FILENO);
    fprintf(stderr, ")\n");
    fprintf(stderr, "         Mismatch between beads present in both files \
can lead to undefined behaviour.\n");
    ResetColour(STDERR_FILENO);
  } //}}}
  // read <type names> from Aggregates command //{{{
  for (int i = 5; i < words && split[i][0] != '-'; i++) {
    int type = FindBeadType(split[i], Counts, BeadType);
    // Error - specified bead type name not in vcf input file
    if (type == -1) {
      ErrorPrintError();
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", input_agg);
      RedText(STDERR_FILENO);
      fprintf(stderr, " - non-existent bead name (");
      YellowText(STDERR_FILENO);
      fprintf(stderr, "%s", split[i]);
      RedText(STDERR_FILENO);
      fprintf(stderr, ") in Aggregate command\n");
      ResetColour(STDERR_FILENO);
      ErrorBeadType(Counts, BeadType);
      exit(1);
    }
    BeadType[type].Use = true;
  } //}}}
  fclose(agg);
} //}}}
