#include "../AnalysisTools.h"

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
Info analyzes the provided input structure file, printing system composition \
to standard output and, optionally, producing an output structure file \
of specified format (-o option). If some information required in the output \
file is missing, '???\' is printed instead. The system from the input file \
can be modified using a second structure file (-i option) and/or \
a coordinate file (-c option); see Examples/Info folder for details.\
\n\n");
  }

  fprintf(ptr, "Usage: %s <input> [options]\n\n", cmd);
  fprintf(ptr, "<input>             input structure file\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  -i <file>         secondary structure file\n");
  fprintf(ptr, "  -c <file>         input coordinate file\n");
  fprintf(ptr, "  --detailed        use name as well as charge, mass, and "
          "radius to identfy bead types (input vtf structure files only)\n");
  fprintf(ptr, "  -o <file>         output structure file\n");
  fprintf(ptr, "  --unique          make all bead/molecule names unique\n");
  fprintf(ptr, "  -def <bead name>  default bead type "
          "(output vtf structure file only)\n");
  fprintf(ptr, "  --mol             make unbonded beads into molecules\n");
  fprintf(ptr, "  --mass            define lammps atom types by mass, but "
          "print per-atom charges in Atoms section "
          "(output lammps data file only)\n");
  fprintf(ptr, "  -ebt <int>        number of extra bead types "
          "(output lammps data file only)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  int vsf_def, b_mol, ebt; // -def --mol -ebt
  bool lmp_mass, detailed; // --mass --detailed
  FILE_TYPE fout;          // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  // define options & check their validity
  int common = 5, all = common + 9, count = 0, req_arg = 1;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, true, option,
               "-st", "--verbose", "--silent", "--help",
               "--version", "-i", "-c", "-o", "--unique",
               "--detailed", "-def", "--mol", "--mass", "-ebt");

  count = 0; // count arguments
  OPT *opt = opt_create();
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, argv[++count], LINE);
  in.stru.type = StructureFileType(in.stru.name);

  PrintCommand(stdout, argc, argv);

  // -i option //{{{
  SYS_FILES extra = InitSysFiles;
  if (FileOption(argc, argv, "-i", extra.stru.name)) {
    extra.stru.type = StructureFileType(extra.stru.name);
  } //}}}
  // input coordinate file (-c option) //{{{
  FileOption(argc, argv, "-c", in.coor.name);
  if (in.coor.name[0] != '\0') {
    in.coor.type = CoordinateFileType(in.coor.name);
    extra.coor = in.coor;
  } //}}}
  // output file (-o option) //{{{
  opt->fout = InitFile;
  if (FileOption(argc, argv, "-o", opt->fout.name)) {
    opt->fout.type = FileType(opt->fout.name);
  } //}}}
  opt->c = CommonOptions(argc, argv, in);
  // extra bead types for data output (-ebt option)
  opt->ebt = 0;
  OneNumberOption(argc, argv, "-ebt", &opt->ebt, 'i');
  // use mass only for atom type definition for data output file
  opt->lmp_mass = BoolOption(argc, argv, "--mass");
  // make unbonded beads into molecules (vtf output only)
  opt->b_mol = BoolOption(argc, argv, "--mol");
  // base bead types on name, charge, mass, and radius (vtf input file)
  opt->detailed = BoolOption(argc, argv, "--detailed");

  // read information from input file(s) //{{{
  SYSTEM System = ReadStructure(in, opt->detailed);
  if (in.stru.type == GROM_FILE) {
    FreeSystem(&System);
    free(opt);
    return 0;
  }
  COUNT *Count = &System.Count;
  // use coordinate from a separate file (-c option)
  if (in.coor.type != -1) {
    int line_count = 0;
    FILE *fr = OpenFile(in.coor.name, "r");
    for (int i = 1; i < opt->c.start; i++) { // from 1 as timestep=1 is the first
      SkipTimestep(in, fr, &line_count);
    }
    ReadTimestep(in, fr, &System, &line_count);
    fclose(fr);
  } else {
    // all beads are in the timestep
    for (int i = 0; i < Count->Bead; i++) {
      System.Bead[i].InTimestep = true;
    }
  }
  // print initial system information only if extra file(s) are present
  if (opt->c.verbose && (extra.stru.type != -1 || in.coor.type != -1)) {
    fprintf(stdout, "\n==================================================");
    printf("\nSystem in %s", in.stru.name);
    if (in.coor.type != -1) {
      printf(" (coordinates: %s)", in.coor.name);
    }
    fprintf(stdout, "\n==================================================\n");
    VerboseOutput(System);
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecules(System);
  }
  SYSTEM Sys_extra;
  if (extra.stru.name[0] != '\0') {
    Sys_extra = ReadStructure(extra, opt->detailed);
    if (opt->c.verbose) {
      fprintf(stdout, "\n==================================================");
      printf("\nSystem in extra file (%s)", extra.stru.name);
      fprintf(stdout, "\n==================================================\n");
      VerboseOutput(Sys_extra);
      fprintf(stdout, "Information about every bead:\n");
      PrintBead(Sys_extra);
      fprintf(stdout, "\nInformation about every molecule:\n");
      PrintMolecules(Sys_extra);
    }
    // add charge, mass, and radius to bead types if possible
    for (int i = 0; i < Count->BeadType; i++) {
      BEADTYPE *bt = &System.BeadType[i];
      int type_extra = FindBeadType(bt->Name, Sys_extra);
      if (type_extra != -1) {
        BEADTYPE *bt_extra = &Sys_extra.BeadType[type_extra];
        if (bt->Charge == CHARGE) {
          bt->Charge = bt_extra->Charge;
        }
        if (bt->Mass == MASS) {
          bt->Mass = bt_extra->Mass;
        }
        if (bt->Radius == RADIUS) {
          bt->Radius = bt_extra->Radius;
        }
      }
    }
    ChangeMolecules(&System, Sys_extra, false);
    CheckSystem(System, extra.stru.name);
  } //}}}

  // -def option (for vsf output file) //{{{
  bool *def_type = calloc(Count->BeadType, sizeof *def_type);
  TypeOption(argc, argv, "-def", 'b', false, def_type, System);
  opt->vsf_def = -1;
  for (int i = 0; i < Count->BeadType; i++) {
    if (def_type[i]) {
      opt->vsf_def = i;
      break;
    }
  }
  free(def_type); //}}}

  if (Count->Bead > 0) {
    PruneSystem(&System);
  }

  // make unbonded beads into molecules //{{{
  if (opt->b_mol) {
    // first new molid is the one higher than the last one
    int mol_id = Count->Molecule;
    // first new resid is one higher than the old one, so +1
    int resid = Count->HighestResid + 1;
    Count->HighestResid += Count->Unbonded;
    Count->Molecule += Count->Unbonded;
    System.Molecule = realloc(System.Molecule,
                              Count->Molecule * sizeof *System.Molecule);
    for (int i = 0; i < Count->BeadType; i++) {
      count = System.BeadType[i].Number;
      bool first = false;
      for (int j = 0; j < count; j++) {
        int id = System.BeadType[i].Index[j];
        BEAD *b = &System.Bead[id];
        BEADTYPE *bt = &System.BeadType[b->Type];
        if (b->Molecule == -1) {
          int n_mt = Count->MoleculeType;
          MOLECULE *mol = &System.Molecule[mol_id];
          mol->Bead = calloc(1, sizeof *mol->Bead);
          if (!first) {
            NewMolType(&System.MoleculeType, &Count->MoleculeType,
                       bt->Name, 1, 0, 0, 0, 0);
            MOLECULETYPE *mt = &System.MoleculeType[n_mt];
            mt->Number = 1;
            mt->Mass = bt->Mass;
            mt->Charge = bt->Charge;
            mt->Bead[0] = b->Type;
            mt->nBTypes = 1;
            mt->BType = malloc(sizeof *mt->BType);
            mt->BType[0] = b->Type;
            mol->Type = n_mt;
            mt->Index = malloc(sizeof *mt->Index);
            first = true;
          } else {
            MOLECULETYPE *mt = &System.MoleculeType[n_mt-1];
            mt->Number++;
            mol->Type = n_mt - 1;
          }
          b->Molecule = mol_id;
          mol->Index = resid;
          mol->Bead[0] = id;
          mol->Aggregate = -1;
          mol->InTimestep = true;
          resid++;
          mol_id++;
        }
      }
    }
    Count->Bonded += Count->Unbonded;
    Count->Unbonded = 0;
    FillBondedUnbonded(&System);
    // fill all molecule type indices
    ReFillMoleculeTypeIndex(&System);
  } //}}}

  // print the system information //{{{
  fprintf(stdout, "\n==================================================");
  if (extra.stru.type != -1 || in.coor.type != -1) {
    printf("\nFinal system composition");
  } else {
    printf("\nSystem composition");
  }
  fprintf(stdout, "\n==================================================\n");
  VerboseOutput(System);
  if (opt->c.verbose) { // -v option
    fprintf(stdout, "Information about every bead:\n");
    PrintBead(System);
    fprintf(stdout, "\nInformation about every molecule:\n");
    PrintMolecules(System);
  } //}}}

  // write the output file if required (-o option) //{{{
  if (opt->fout.name[0] != '\0') {
    if (opt->fout.type == LDATA_FILE && opt->ebt > 0) {
      NewBeadType(&System.BeadType, &Count->BeadType, "extra", 0, 1, 1);
    }
    bool *write = malloc(sizeof *write * Count->Bead);
    InitBoolArray(write, Count->Bead, true);
    WriteOutput(System, write, opt->fout, opt->lmp_mass, opt->vsf_def,
                argc, argv);
    free(write);
  } //}}}

  FreeSystem(&System);
  if (extra.stru.name[0] != '\0') {
    FreeSystem(&Sys_extra);
  }
  free(opt);

  return 0;
}
