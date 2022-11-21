#include "../AnalysisTools.h"
int *InFile;

// TODO: add that it's only for cuboid simulation box
// TODO add --silent option

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(stdout, "\
GenLayers reads information from a FIELD-like file and generates two \
mirror monolayers in xy planes of the simulation box specified distance \
from edges of the box (in z direction). Only the first molecule type in \
the FIELD-like file is used. The molecules are \
generated on a square grid with distance between anchoring points specified \
either explicitly or implicitly (by specifying the number of molecules per \
layer).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <out.vsf> <out.vcf> [options]\n\n", cmd);
  fprintf(ptr, "   <out.vsf>   output structure file (vsf format)\n");
  fprintf(ptr, "   <out.vcf>   output coordinate file (vcf format)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -s <float> <float>   spacing in x and y directions \
(default: 1 1)\n");
  fprintf(ptr, "      -nm <int>            total number molecules per side \
(rewrites spacing to fit)\n");
  fprintf(ptr, "      -g <float>           gap between walls and the layers\n");
  fprintf(ptr, "      -f <name>            FIELD-like file (default: FIELD)\n");
  fprintf(ptr, "      -v                   verbose output\n");
  fprintf(ptr, "      -h                   print this help and exit\n");
  fprintf(ptr, "      --version            print version number and exit\n");
} //}}}

int main(int argc, char *argv[]) {

  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 2; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count+1) < argc && argv[count+1][0] != '-') {
    count++;
  }

  if (count < req_args) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' &&
        strcmp(argv[i], "-s") != 0 &&
        strcmp(argv[i], "-nm") != 0 &&
        strcmp(argv[i], "-g") != 0 &&
        strcmp(argv[i], "-f") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-h") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count arguments

  // <out.vsf> - output structure file //{{{
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]);
  // test if <output.vsf> ends with '.vsf'
  int ext = 1;
  char extension[2][EXTENSION];
  strcpy(extension[0], ".vsf");
  if (ErrorExtension(output, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <out.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vsf'
  ext = 1;
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // print command to stdout
  PrintCommand(stdout, argc, argv);

  // options before reading system data //{{{
  // '-s' option - spacing in x and y directions //{{{
  int test = 2;
  double spacing[100]; // 100 in case of too many arguments provided
  spacing[0] = 1;
  spacing[1] = 1;
  if (MultiDoubleOption(argc, argv, "-s", &test, spacing)) {
    exit(1);
  }
  if (test != 2) {
    ErrorPrintError_old();
    ColourChange(STDERR_FILENO, YELLOW);
    fprintf(stderr, "-s");
    ColourChange(STDERR_FILENO, RED);
    fprintf(stderr, " - two numeric arguments required\n");
    ColourReset(STDERR_FILENO);
    exit(1);
  } //}}}

  // '-nm' option - custom total number of chains //{{{
  int number_of_mols = 0;
  if (IntegerOption(argc, argv, "-nm", &number_of_mols)) {
    exit(1);
  } //}}}

  // '-g' option - gap from both sides of the box in z direction //{{{
  double gap = 0;
  if (DoubleOption(argc, argv, "-g", &gap)) {
    exit(1);
  } //}}}

  // FIELD-like file //{{{
  char input[LINE] = "";
  if (FileOption(argc, argv, "-f", input, LINE)) {
    exit(1);
  }
  if (input[0] == '\0') {
    strcpy(input, "FIELD");
  } //}}}
  bool verbose = BoolOption(argc, argv, "-v");
  //}}}

  // read info from the FIELD file //{{{
  BEADTYPE *BeadType_field;
  MOLECULETYPE *MoleculeType_field;
  MOLECULE *Molecule_field;
  BEAD *Bead_field;
  COUNTS Counts_field = InitCounts;
  BOX Box = InitBox; // triclinic box dimensions and angles
  int *Index_field;
  PARAMS *bond_type_field;
  PARAMS *angle_type_field;
  PARAMS *dihedral_type_field;

  FieldRead(input, &Box.Length, &Counts_field, &BeadType_field, &Bead_field,
            &Index_field, &MoleculeType_field, &Molecule_field,
            &bond_type_field, &angle_type_field, &dihedral_type_field); //}}}

  // set all bead types to not use to later determine which are used //{{{
  for (int i = 0; i < Counts_field.TypesOfBeads; i++) {
    BeadType_field[i].Use = false;
  } //}}}

  // structures for info for final system //{{{
  BEADTYPE *BeadType;
  MOLECULETYPE *MoleculeType;
  MOLECULE *Molecule;
  BEAD *Bead;
  COUNTS Counts = InitCounts;
  PARAMS *bond_type;
  PARAMS *angle_type;
  PARAMS *dihedral_type;
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i) //}}}

  // find bead types to be used & copy them //{{{
  /*
   * 1) signal which are to be used via Use=true & count them
   * 2) copy the bead types to a new array
   */
  // 1)
  Counts.TypesOfBeads = 0;
  for (int j = 0; j < MoleculeType_field[0].nBeads; j++) {
    int btype = MoleculeType_field[0].Bead[j];
    if (!BeadType_field[btype].Use) {
      Counts.TypesOfBeads++;
    }
    BeadType_field[btype].Use = true;
  }
  // 2)
  BeadType = calloc(Counts.TypesOfBeads, sizeof(BEADTYPE));
  count = 0;
  for (int i = 0; i < Counts_field.TypesOfBeads; i++) {
    if (BeadType_field[i].Use) {
      strcpy(BeadType[count].Name, BeadType_field[i].Name);
      BeadType[count].Charge = BeadType_field[i].Charge;
      BeadType[count].Mass = BeadType_field[i].Mass;
      BeadType[count].Number = 0; // will be counted later
      count++;
    }
  } //}}}

  // only one molecule type is used - copy it to a new struct //{{{
  Counts.TypesOfMolecules = 1;
  MoleculeType = calloc(Counts.TypesOfMolecules, sizeof(MOLECULETYPE));
  strcpy(MoleculeType[0].Name, MoleculeType_field[0].Name);
  // copy beads (using the new bead type indices)
  MoleculeType[0].nBeads = MoleculeType_field[0].nBeads;
  MoleculeType[0].Bead = calloc(MoleculeType[0].nBeads, sizeof(int));
  for (int j = 0; j < MoleculeType[0].nBeads; j++) {
    int bt_field = MoleculeType_field[0].Bead[j];
    MoleculeType[0].Bead[j] = FindBeadType2(BeadType_field[bt_field].Name,
                                            Counts.TypesOfBeads, BeadType);
  }
  // copy bonds & count bond types
  Counts.TypesOfBonds = 0;
  MoleculeType[0].nBonds = MoleculeType_field[0].nBonds;
  MoleculeType[0].Bond = malloc(sizeof *MoleculeType[0].Bond *
                                MoleculeType[0].nBonds);
  for (int j = 0; j < MoleculeType[0].nBonds; j++) {
    MoleculeType[0].Bond[j][0] = MoleculeType_field[0].Bond[j][0];
    MoleculeType[0].Bond[j][1] = MoleculeType_field[0].Bond[j][1];
    MoleculeType[0].Bond[j][2] = MoleculeType_field[0].Bond[j][2];
    // number of bond types is the highest type id in the molecule
    if (MoleculeType[0].Bond[j][2] >= Counts.TypesOfBonds) {
      Counts.TypesOfBonds = MoleculeType[0].Bond[j][2] + 1;
    }
  }
  // copy angles & count angles types
  Counts.TypesOfAngles = 0;
  MoleculeType[0].nAngles = MoleculeType_field[0].nAngles;
  MoleculeType[0].Angle = malloc(sizeof *MoleculeType[0].Angle *
                                 MoleculeType[0].nAngles);
  for (int j = 0; j < MoleculeType[0].nAngles; j++) {
    MoleculeType[0].Angle[j][0] = MoleculeType_field[0].Angle[j][0];
    MoleculeType[0].Angle[j][1] = MoleculeType_field[0].Angle[j][1];
    MoleculeType[0].Angle[j][2] = MoleculeType_field[0].Angle[j][2];
    MoleculeType[0].Angle[j][3] = MoleculeType_field[0].Angle[j][3];
    // number of angle types is the highest angle type id in the molecule
    if (MoleculeType[0].Angle[j][3] >= Counts.TypesOfAngles) {
      Counts.TypesOfAngles = MoleculeType[0].Angle[j][3] + 1;
    }
  }
  // copy dihedrals & count dihedrals types
  Counts.TypesOfDihedrals = 0;
  MoleculeType[0].nDihedrals = MoleculeType_field[0].nDihedrals;
  MoleculeType[0].Dihedral = malloc(sizeof *MoleculeType[0].Dihedral *
                                    MoleculeType[0].nDihedrals);
  for (int j = 0; j < MoleculeType[0].nDihedrals; j++) {
    MoleculeType[0].Dihedral[j][0] = MoleculeType_field[0].Dihedral[j][0];
    MoleculeType[0].Dihedral[j][1] = MoleculeType_field[0].Dihedral[j][1];
    MoleculeType[0].Dihedral[j][2] = MoleculeType_field[0].Dihedral[j][2];
    MoleculeType[0].Dihedral[j][3] = MoleculeType_field[0].Dihedral[j][3];
    MoleculeType[0].Dihedral[j][4] = MoleculeType_field[0].Dihedral[j][4];
    // number of dihedral types is the highest dihedral type id in the molecule
    if (MoleculeType[0].Dihedral[j][4] >= Counts.TypesOfDihedrals) {
      Counts.TypesOfDihedrals = MoleculeType[0].Dihedral[j][4] + 1;
    }
  }
  // copy mass & charge & fill the BType array
  MoleculeType[0].Mass = MoleculeType_field[0].Mass;
  MoleculeType[0].Charge = MoleculeType_field[0].Charge;
  // fill BTypes struct using the new bead type indices
  FillMolBTypes_old(Counts.TypesOfMolecules, &MoleculeType); //}}}

  // copy bond & angle & dihedral types to be used //{{{
  bond_type = calloc(Counts.TypesOfBonds, sizeof(PARAMS));
  for (int i = 0; i < Counts.TypesOfBonds; i++) {
    bond_type[i].a = bond_type_field[i].a;
    bond_type[i].b = bond_type_field[i].b;
  }
  angle_type = calloc(Counts.TypesOfAngles, sizeof(PARAMS));
  for (int i = 0; i < Counts.TypesOfAngles; i++) {
    angle_type[i].a = angle_type_field[i].a;
    angle_type[i].b = angle_type_field[i].b;
  }
  dihedral_type = calloc(Counts.TypesOfDihedrals, sizeof(PARAMS));
  for (int i = 0; i < Counts.TypesOfDihedrals; i++) {
    dihedral_type[i].a = dihedral_type_field[i].a;
    dihedral_type[i].b = dihedral_type_field[i].b;
  } //}}}

  // number of molecules in x and y directions //{{{
  int mols[2]; // maximum number of molecules in x and y directions per wall
  if (number_of_mols == 0) { // if -nm option is not used
    mols[0] = round(Box.Length.x / spacing[0]);
    mols[1] = round(Box.Length.y / spacing[1]);
    Counts.Molecules = mols[0] * mols[1] * 2;
  } else { // if -nm option is used
    double mols_sqrt = sqrt(number_of_mols);
    mols[0] = (int)(mols_sqrt);
    mols[1] = mols[0] + (number_of_mols - SQR(mols[0])) / mols[0];
    if ((mols[0]*mols[1]-number_of_mols) != 0) {
      mols[0]++;
    }
    spacing[0] = Box.Length.x / mols[0];
    spacing[1] = Box.Length.y / mols[1];
    Counts.Molecules = 2 * number_of_mols;
  }
  MoleculeType[0].Number = Counts.Molecules; //}}}

  // determine number of beads of each type //{{{
  for (int i = 0; i < MoleculeType[0].nBeads; i++) {
    int btype = MoleculeType[0].Bead[i];
    BeadType[btype].Number += MoleculeType[0].Number;
  } //}}}

  // total number of beads //{{{
  Counts.Bonded = MoleculeType[0].Number * MoleculeType[0].nBeads;
  Counts.BeadsTotal = Counts.Bonded;
  Counts.BeadsCoor = Counts.Bonded;
  Counts.Unbonded = 0; //}}}

  // allocate memory
  Molecule = calloc(Counts.Molecules, sizeof(MOLECULE));
  Bead = calloc(Counts.BeadsCoor, sizeof(BEAD));
  Index = calloc(Counts.BeadsCoor, sizeof(int));

  // fill MOLECULE struct //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    Molecule[i].Type = 0;
    Molecule[i].Bead = malloc(MoleculeType[0].nBeads*sizeof(int));
    for (int j = 0; j < MoleculeType[0].nBeads; j++) {
      Molecule[i].Bead[j] = 0;
    }
    for (int j = 0; j < MoleculeType[0].nBeads; j++) {
      count = i * MoleculeType[0].nBeads + j;
      Molecule[i].Bead[j] = count;
      Bead[count].Type = MoleculeType[0].Bead[j];
      Bead[count].Molecule = i;
      Bead[count].Indexxx = count;
      Bead[count].Aggregatexxx = calloc(1, sizeof(int)); // just to free later
      Index[Bead[count].Indexxx] = count;
    }
  } //}}}

  // set all molecules & beads to write //{{{
  MoleculeType[0].Write = true;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    BeadType[i].Write = true;
  } //}}}

  // move the molecule, so its first bead is [0,0,0] //{{{
  for (int i = 1; i < MoleculeType[0].nBeads; i++) {
    int id_field = Molecule_field[0].Bead[i];
    Bead_field[id_field].Position.x -= Bead_field[id_field-i].Position.x;
  }
  // the first bead's coordinates aren't used, so this is just to be sure
  int id_field = Molecule_field[0].Bead[0];
  Bead_field[id_field].Position.x = 0;
  Bead_field[id_field].Position.y = 0;
  Bead_field[id_field].Position.z = 0; //}}}

  // generate brush on the grid & fill remaining Bead[] array //{{{
  count = 0;
  for (int layer = 0; layer < 2; layer++) { // 0 for z=0; 1 for z=BoxLength.z
    int count_mols = 0; // per layer
    for (int i = 0; i < mols[0]; i++) {
      for (int j = 0; j < mols[1]; j++) {
        // first bead of a molecule
        Bead[count].Position.x = spacing[0] / 2 + spacing[0] * i;
        Bead[count].Position.y = spacing[1] / 2 + spacing[1] * j;
        if (layer == 0) {
          Bead[count].Position.z = 0;
        } else {
          Bead[count].Position.z = Box.Length.z;
        }
        Bead[count].Indexxx = count;
        Bead[count].Type = MoleculeType[0].Bead[0];
        // remaining beads of the molecule
        for (int k = 1; k < MoleculeType[0].nBeads; k++) {
          int id_field = Molecule_field[0].Bead[k];
          Bead[count+k].Position.x = Bead[count].Position.x +
                                     Bead_field[id_field].Position.x;
          Bead[count+k].Position.y = Bead[count].Position.y +
                                     Bead_field[id_field].Position.y;
          if(layer == 0) {
            Bead[count+k].Position.z = Bead[count].Position.z +
                                       Bead_field[id_field].Position.z;
          } else {
            Bead[count+k].Position.z = Bead[count].Position.z -
                                       Bead_field[id_field].Position.z;
          }
          Bead[count+k].Indexxx = count + k;
          Bead[count+k].Type = MoleculeType[0].Bead[k];
        }
        count += MoleculeType[0].nBeads;
        count_mols++;
        // stop placing molecules (-nm option)
        if (count_mols == number_of_mols) {
          break;
        }
      }
      if (count_mols == number_of_mols) {
        break;
      }
    }
  } //}}}

  // adjust coordinates for gap from box side (-g option) //{{{
  for (int i = 0; i < Counts.Molecules; i++) {
    int mtype = Molecule[i].Type;
    for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
      int id = Molecule[i].Bead[j];
      if (i < (Counts.Molecules/2)) {
        Bead[id].Position.z += gap;
      } else {
        Bead[id].Position.z -= gap;
      }
    }
  } //}}}

  // write output vsf file
  WriteVsf_old(output, Counts, BeadType, Bead, MoleculeType, Molecule, false);

  // open output .vcf file for writing
  FILE *out = OpenFile(output_vcf, "w");
  PrintByline(out, argc, argv);
  WriteCoorIndexed(out, Counts, BeadType, Bead,
                   MoleculeType, Molecule, "\0", Box);
  fclose(out);

  // print information - verbose option //{{{
  if (verbose) {
    fprintf(stdout, "\nNote that only the first molecule in the ");
    fprintf(stdout, "%s file (and the associated beads, ", input);
    fprintf(stdout, "bonds, angles, and dihedrals) are used.\n");
    fprintf(stdout, "\nGrid of %d x %d ", mols[0], mols[1]);
    if (number_of_mols != 0) {
      fprintf(stdout, "(%d molecules) ", number_of_mols);
    }
    fprintf(stdout, "molecules on each wall\n");
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
    PrintBondTypes2(Counts.TypesOfBonds, bond_type);
    PrintAngleTypes2(Counts.TypesOfAngles, angle_type);
    PrintDihedralTypes2(Counts.TypesOfDihedrals, dihedral_type);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts_field, &MoleculeType_field, &Molecule_field,
                 &BeadType_field, &Bead_field, &Index_field);
  FreeSystemInfo(Counts, &MoleculeType, &Molecule,
                 &BeadType, &Bead, &Index);
  free(bond_type_field);
  free(angle_type_field);
  free(dihedral_type_field);
  free(bond_type);
  free(angle_type);
  free(dihedral_type); //}}}

  return 0;
}
