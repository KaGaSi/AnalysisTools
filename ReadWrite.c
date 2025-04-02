#include "ReadWrite.h"
#include "General.h"
#include "ReadWriteVtf.h"
#include "ReadWriteXyz.h"
#include "ReadWriteLtrj.h"
#include "ReadWriteLdata.h"
#include "ReadWriteField.h"
#include "ReadWriteConfig.h"
#include "ReadWriteGromacs.h"
#include "AnalysisTools.h"
#include "Errors.h"
#include "Structs.h"
#include "System.h"

static void CopyAndFreeStuff(const int n, int (**old)[5], int (**new)[5]);
static void CopyAndFreeAllStuff(MOLECULETYPE *mt_old, MOLECULETYPE *mt_new);
static void MinimizeOneMtypeStuffIds(const int num, int (**arr)[5],
                                     const int n, const int *link_bead_ids);
static void FillMTypeStuff(SYSTEM *System, const int type, const int size,
                           const int (*all)[5], const int n);
static void PrintBeadHeader();
static void PrintOneBead(const SYSTEM System, const int id);

// General helper functions
// print initial stuff to output coordinate file //{{{
void InitOutputCoorFile(const FILE_TYPE file, const SYSTEM System,
                        const int argc, char **argv) {
  if (file.type == VCF_FILE) {
    PrintByline(file.name, argc, argv);
  } else if (file.type == VTF_FILE) {
    WriteStructure(file, System, -1, false, argc, argv);
  } else {
    FILE *out = OpenFile(file.name, "w");
    fclose(out);
  }
} //}}}
void CopyMoleculeTypeBeadsToMoleculeBeads(SYSTEM *System) { //{{{
  COUNT *Count = &System->Count;
  for (int i = 0; i < Count->Molecule; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    MOLECULE *mol_i = &System->Molecule[i];
    if (mt_i->nBeads == 1 && !mt_i->Flag) { // remove 'fake' molecules
      mt_i->Number = 0;
      mt_i->nBeads = 0;
      System->Bead[mt_i->Bead[0]].Molecule = -1;
      free(mt_i->Bead);
      Count->Bonded--;
      Count->Unbonded++;
    }
    if (mt_i->Number == 1) { // copy beads only if molecule with id 'i' exists
      mol_i->Type = i;
      mol_i->Index = i;
      mol_i->Bead = malloc(sizeof *mol_i->Bead * mt_i->nBeads);
      for (int j = 0; j < mt_i->nBeads; j++) {
        mol_i->Bead[j] = mt_i->Bead[j];
      }
    }
  }
} //}}}
// make the MoleculeType[].Stuff bead indices be between 0 and nBeads //{{{
// Assumes that there's one Molecule per MoleculeType
void MinimizeMTypeStuffIds(SYSTEM *System) {
  /*
   * Allocate memory for linking bead indices in the molecule (i.e., bead
   * indices in the input structure file) to bead indices in the molecule
   * type (i.e., indices 0 to nBeads).
   */
  int *link_ids = calloc(System->Count.Bead, sizeof *link_ids);
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    MOLECULETYPE *mt_i = &System->MoleculeType[i];
    int count_bead = 0;
    for (int j = 0; j < mt_i->nBeads; j++) {
      int id = System->Molecule[i].Bead[j];
      link_ids[id] = count_bead;
      count_bead++;
    }
    MinimizeOneMtypeStuffIds(2, &mt_i->Bond, mt_i->nBonds, link_ids);
    MinimizeOneMtypeStuffIds(3, &mt_i->Angle, mt_i->nAngles, link_ids);
    MinimizeOneMtypeStuffIds(4, &mt_i->Dihedral, mt_i->nDihedrals, link_ids);
    MinimizeOneMtypeStuffIds(4, &mt_i->Improper, mt_i->nImpropers, link_ids);
  }
  free(link_ids);
}
static void MinimizeOneMtypeStuffIds(const int num, int (**arr)[5],
                                     const int n, const int *link_bead_ids) {
  if (n > 0) { // only for the molecule types with that stuff
    // use the linking array to assign proper bead ids
    for (int j = 0; j < n; j++) {
      for (int aa = 0; aa < num; aa++) {
        (*arr)[j][aa] = link_bead_ids[(*arr)[j][aa]];
      }
    }
  }
} //}}}
// fill MoleculeType Bond/Angle/Dihedral/Improper arrays //{{{
void FillAllMTypeStuff(SYSTEM *System, const int (*bond)[5], const int
                       (*angle)[5], const int (*dih)[5], const int (*imp)[5]) {
  FillMTypeStuff(System, 0, 3, bond, System->Count.Bond);
  FillMTypeStuff(System, 1, 4, angle, System->Count.Angle);
  FillMTypeStuff(System, 2, 5, dih, System->Count.Dihedral);
  FillMTypeStuff(System, 3, 5, imp, System->Count.Improper);
  MinimizeMTypeStuffIds(System);
}
static void FillMTypeStuff(SYSTEM *System, const int type, const int size,
                           const int (*all)[5], const int n) {
  // fill MoleculeType[].Stuff array with bead indices
  for (int i = 0; i < n; i++) {
    int id[size];
    for (int aa = 0; aa < size; aa++) {
      id[aa] = all[i][aa];
    }
    int mol = System->Bead[id[0]].Molecule;
    // warning - beads in different molecules (skip it)  //{{{
    bool err = false;
    for (int aa = 1; aa < (size - 1); aa++) {
      if (mol != System->Bead[id[aa]].Molecule || mol == -1) {
        err_msg("Discarding bond/angle/dihedral/improper with beads that"
                " do not share a molecule");
        PrintWarning();
        fprintf(stderr, "%sBead (molecule):", ErrCyan());
        for (int bb = 1; bb < (size - 1); bb++) {
          fprintf(stderr, " %s%d%s (%s%d%s)", ErrYellow(), id[bb], ErrCyan(),
                  ErrYellow(), System->Bead[id[bb]].Molecule, ErrCyan());
        }
        putc('\n', stderr);
        err = true;
        break;
      }
    }
    if (err) {
      continue;
    } //}}}
    MOLECULETYPE *mt_mol = &System->MoleculeType[mol];
    int (**arr)[5] = NULL;
    int n_stuff = 0;
    int *count = NULL;
    if (type == 0) {
      arr = &mt_mol->Bond;
      n_stuff = mt_mol->nBonds;
      count = &mt_mol->nBonds;
    } else if (type == 1) {
      arr = &mt_mol->Angle;
      n_stuff = mt_mol->nAngles;
      count = &mt_mol->nAngles;
    } else if (type == 2) {
      arr = &mt_mol->Dihedral;
      n_stuff = mt_mol->nDihedrals;
      count = &mt_mol->nDihedrals;
    } else if (type == 3) {
      arr = &mt_mol->Improper;
      n_stuff = mt_mol->nImpropers;
      count = &mt_mol->nImpropers;
    }
    (*count)++;
    if (n_stuff == 0) {
      *arr = malloc(sizeof **arr);
    } else {
      *arr = s_realloc(*arr, sizeof **arr * *count);
    }
    for (int aa = 0; aa < size; aa++) {
      (*arr)[n_stuff][aa] = id[aa];
    }
  }
} //}}}
// copy bonds/angles/dihedrals/impropers to a new molecule type //{{{
// Also frees the Stuff array from the old molecule type
static void CopyAndFreeStuff(const int n, int (**old)[5], int (**new)[5]) {
  if (n > 0) {
    *new = malloc(sizeof **new * n);
    for (int j = 0; j < n; j++) {
      for (int aa = 0; aa < 5; aa++) {
        (*new)[j][aa] = (*old)[j][aa];
      }
    }
    free(*old);
  }
}
// Also frees the Bond/Angle/etc/ arrays
static void CopyAndFreeAllStuff(MOLECULETYPE *mt_old, MOLECULETYPE *mt_new) {
  CopyAndFreeStuff(mt_new->nBonds, &mt_old->Bond, &mt_new->Bond);
  CopyAndFreeStuff(mt_new->nAngles, &mt_old->Angle, &mt_new->Angle);
  CopyAndFreeStuff(mt_new->nDihedrals, &mt_old->Dihedral, &mt_new->Dihedral);
  CopyAndFreeStuff(mt_new->nImpropers, &mt_old->Improper, &mt_new->Improper);
} //}}}
// RemoveExtraTypes() { //{{{
/*
 * Remove bead and molecule types with .Number=0. It assumes the allocated
 * memory for BeadType and MoleculeType arrays of structures correspond to the
 * number of beads and molecules, respectively (i.e., not to the number of
 * types).
 */
void RemoveExtraTypes(SYSTEM *System) {
  COUNT *Count = &System->Count;
  if (Count->Bead > 0) {
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
    // MoleculeType & Molecule
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
          for (int j = 0; j < mt_new->nBeads; j++) {
            mt_new->Bead[j] = mt_i->Bead[j];
          }
          free(mt_i->Bead);
          CopyAndFreeAllStuff(mt_i, mt_new);
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
  }
} //}}}
void WriteBoxLengthAngles(FILE *fw, const BOX box) { //{{{
  if (box.Volume != -1) {
    fprintf(fw, "pbc %.3f %.3f %.3f", box.Length[0],
                                      box.Length[1],
                                      box.Length[2]);
    if (fabs(box.alpha - 90) > 0.00001 ||
        fabs(box.beta - 90) > 0.00001 ||
        fabs(box.gamma - 90) > 0.00001) {
      fprintf(fw, " %lf %lf %lf", box.alpha, box.beta, box.gamma);
    }
  }
} //}}}

// Read a structure file and extra info from coordinate file if necessary //{{{
SYSTEM ReadStructure(const SYS_FILES f, const bool detailed) {
  SYSTEM System;
  switch (f.stru.type) {
    case VTF_FILE:
    case VSF_FILE:
      System = VtfReadStruct(f.stru.name, detailed);
      break;
    case XYZ_FILE:
      System = XyzReadStruct(f.stru.name);
      break;
    case LTRJ_FILE:
      System = LtrjReadStruct(f.stru.name);
      break;
    case LDATA_FILE:
      System = LmpDataReadStruct(f.stru.name);
      break;
    case FIELD_FILE:
      System = FieldRead(f.stru.name);
      break;
    case ITP_FILE:
      System = ItpReadStruct(f.stru.name);
      break;
    case PDB_FILE:
      System = PdbReadStruct(f.stru.name);
      break;
    default:
      err_msg("unspecified structure file; should never happen!");
      PrintError();
      exit(1);
  }
  // read extra stuff from coordinate file if necessary
  if (f.coor.type == LTRJ_FILE && System.Box.Volume == -1) {
    System.Box = LtrjReadPBC(f.coor.name);
  }
  if (f.coor.type == VCF_FILE) {
    System.Count.BeadCoor = VtfReadNumberOfBeads(f.coor.name);
    if (System.Count.BeadCoor < 0) {
      exit(1);
    }
    if (System.Box.Volume == -1) {
      System.Box = VtfReadPBC(f.coor.name);
    }
  }
  WarnChargedSystem(System, f.stru.name, "\0", "\0");
  // warn if missing box dimensions (unless it's a pbc-less file type)
  if (System.Box.Volume == -1 &&
      f.stru.type != VSF_FILE &&
      f.stru.type != FIELD_FILE &&
      f.stru.type != XYZ_FILE &&
      f.stru.type != ITP_FILE) {
    err_msg("unspecified box dimensions in structure definition");
    PrintWarnFile(f.stru.name, "\0", "\0");
  }
  return System;
} //}}} //}}}
// Read a single timestep from the provided coordinate file //{{{
bool ReadTimestep(const SYS_FILES f, FILE *fr,
                  SYSTEM *System, int *line_count) {
  switch (f.coor.type) {
    case VTF_FILE:
    case VCF_FILE:
      if (VtfReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case XYZ_FILE:
      if (XyzReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case LTRJ_FILE:;
      int line = *line_count;
      if (LtrjReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      // skip this step if it's a first one that contain only zeroes
      // ...huh? Why would this be a thing?
      if (line == 0 && System->Count.BeadCoor == System->Count.Bead) {
        bool zeroes = true;
        for (int i = 0; i < System->Count.Bead; i++) {
          if (System->Bead[i].Position[0] != 0 ||
              System->Bead[i].Position[1] != 0 ||
              System->Bead[i].Position[2] != 0) {
            zeroes = false;
            break;
          }
        }
        if (zeroes &&
            LtrjReadTimestep(fr, f.coor.name, System, line_count) < 0) {
          return false;
        }
      }
      break;
    case LDATA_FILE:
      if (LmpDataReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for coor_type %s%d",
               ErrYellow(), f.coor.type);
      PrintError();
      exit(1);
  }
  return true;
} //}}}
// Skip a single timestep from the provided coordinate file //{{{
bool SkipTimestep(const SYS_FILES f, FILE *fr, int *line_count) {
  switch (f.coor.type) {
    case VTF_FILE:
    case VCF_FILE:
      if (VtfSkipTimestep(fr, f.coor.name, f.stru.name, line_count) < 0) {
        return false;
      }
      break;
    case XYZ_FILE:
      if (!XyzSkipTimestep(fr, f.coor.name, line_count)) {
        return false;
      }
      break;
    case LTRJ_FILE:
      if (LtrjSkipTimestep(fr, f.coor.name, line_count) < 0) {
        return false;
      }
      break;
    case LDATA_FILE:
      err_msg("lammps data file contains only one step; should never trigger!");
      PrintWarnFile(f.coor.name, "\0", "\0");
      return false;
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for coor_type %s%d",
               ErrYellow(), f.coor.type);
      PrintError();
      exit(1);
  }
  return true;
} //}}}
// Read aggregates from a single timestep //{{{
int ReadAggregates(FILE *fr, const char *file, SYSTEM *System,
                   AGGREGATE *Aggregate, int *line_count) {
  COUNT *Count = &System->Count;
  // read Step:/Last Step: line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE,
             "premature end of %s%s%s file; "
             "stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return -2;
  }
  if (words < 1) {
    snprintf(ERROR_MSG, LINE,
             "blank line instead of Step:/Last Step: line; "
             "stopping reading %s%s",
             ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return -2;
  }
  //}}}
  if (strcasecmp(split[0], "Last") == 0) {
    return -1; // no error - just the end of the timesteps
  }
  // read number of aggregates //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE, "premature end of %s%s%s file; stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return -2;
  }
  long int val;
  if (words < 1 || !IsWholeNumber(split[0], &val)) {
    snprintf(ERROR_MSG, LINE,
             "incorrect line with number of aggregates stopping reading %s%s",
             ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return -2;
  } //}}}
  Count->Aggregate = val;
  // go through all aggregates, filling in the molecules
  for (int i = 0; i < Count->Aggregate; i++) {
    (*line_count)++;
    AGGREGATE *Agg = &Aggregate[i];
    fscanf(fr, "%d :", &Agg->nMolecules);
    Agg->Molecule = s_realloc(Agg->Molecule,
                              Agg->nMolecules * sizeof *Agg->Molecule);
    for (int j = 0; j < Agg->nMolecules; j++) {
      int mol;
      fscanf(fr, "%d", &mol);
      mol--; // in agg file, the numbers correspond to vmd
      Agg->Molecule[j] = mol;
      System->Molecule[mol].Aggregate = i;
    }
    int ch;
    while ((ch = getc(fr)) != '\n' && ch != EOF)
      ;
  }
  // fill the rest of aggregate info
  for (int i = 0; i < Count->Aggregate; i++) {
    int count = 0;
    double mass = 0;
    AGGREGATE *Agg = &Aggregate[i];
    for (int j = 0; j < Agg->nMolecules; j++) {
      MOLECULE *mol = &System->Molecule[Agg->Molecule[j]];
      MOLECULETYPE *mt = &System->MoleculeType[mol->Type];
      if (mt->Mass == MASS || mass == -1) {
        mass = -1;
      } else {
        mass += mt->Mass;
      }
      mol->Aggregate = i;
      Agg->Bead = s_realloc(Agg->Bead,
                            (count + mt->nBeads) * sizeof *Agg->Bead);
      for (int k = 0; k < mt->nBeads; k++) {
        Agg->Bead[count] = mol->Bead[k];
        count++;
      }
    }
    Agg->nBeads = count;
    if (mass == -1) {
      Agg->Mass = MASS; // unspecified mass
    } else {
      Agg->Mass = mass; // valid mass
    }
  }
  return 1;
} //}}}
bool SkipAggregates(FILE *fr, const char *file, int *line_count) { //{{{
  // read Step:/Last Step: line //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE, "premature end of %s%s%s file; stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return false;
  }
  if (words < 1) {
    snprintf(ERROR_MSG, LINE, "blank line instead of Step:/Last Step: line; "
             "stopping reading %s%s", ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return false;
  }
  //}}}
  if (strcasecmp(split[0], "Last") == 0) {
    return false; // no error - just the end of the timesteps
  }
  // read number of aggregates //{{{
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    snprintf(ERROR_MSG, LINE, "premature end of %s%s%s file; stopping reading",
             ErrYellow(), file, ErrCyan());
    PrintWarnFile(file, "\0", "\0");
    return false;
  }
  long int num_agg;
  if (words < 1 || !IsWholeNumber(split[0], &num_agg)) {
    snprintf(ERROR_MSG, LINE, "incorrect line with number of aggregates; "
             "stopping reading %s%s", ErrYellow(), file);
    PrintWarnFileLine(file, *line_count);
    return false;
  } //}}}
  // go through all aggregates, filling in the molecules
  for (int i = 0; i < num_agg; i++) {
    (*line_count)++;
    int ch;
    while ((ch = getc(fr)) != '\n' && ch != EOF)
      ;
  }
  return true;
} //}}}

// write structure and/or coordinates to a new file (can be any format) //{{{
void WriteOutput(const SYSTEM System, const bool *write, FILE_TYPE fw,
                 const bool lmp_mass, const int vsf_def,
                 const int argc, char **argv) {
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
void WriteOutputAll(const SYSTEM System, FILE_TYPE fw, const bool lmp_mass,
                    const int vsf_def, const int argc, char **argv) {
  bool *write = malloc(System.Count.Bead * sizeof *write);
  InitBoolArray(write, System.Count.Bead, true);
  WriteOutput(System, write, fw, lmp_mass, vsf_def, argc, argv);
  free(write);
} //}}}
// Write a single timestep to output file based on the file type //{{{
void WriteTimestep(const FILE_TYPE f, const SYSTEM System, const int count_step,
                   const bool *write, const int argc, char **argv) {
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
    case CONFIG_FILE:
      WriteConfig(System, f.name);
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
                      int argc, char **argv) {
  bool *write = malloc(System.Count.Bead * sizeof *write);
  InitBoolArray(write, System.Count.Bead, true);
  WriteTimestep(f, System, count_step, write, argc, argv);
  free(write);
}
//}}}
// Create a structure file based on the file type (including dl_meso CONFIG) //{{{
void WriteStructure(FILE_TYPE f, SYSTEM System, const int vsf_def_type,
                    const bool lmp_mass, const int argc, char **argv) {
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
void WriteAggregates(const int step_count, const char *agg_file,
                     const SYSTEM System, const AGGREGATE *Aggregate) {
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
  // if (Count.Molecule > 0) {
  //   fprintf(stdout, "  HighestResid:   %d", Count.HighestResid);
  // }
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
// print bead info //{{{
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
  fprintf(stdout, " %s", b->InTimestep ? "yes" : " no");
  putchar('\n');
}
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
} //}}}
// find highest values from {Bond,Angle,Dihedral,Improper}Type params //{{{
int * HighestParam(const SYSTEM System, const int type) {
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
  width[0] = snprintf(NULL, 0, "%.5f", high.a);
  width[1] = snprintf(NULL, 0, "%.0f", high.b);
  width[2] = snprintf(NULL, 0, "%.0f", high.c);
  width[3] = snprintf(NULL, 0, "%.0f", high.d);
  return width;
} //}}}
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
void PrintByline(const char *file, const int argc, char **argv) { //{{{
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "# Created by AnalysisTools v%s ", VERSION);
  fprintf(fw, " (https://github.com/KaGaSi/AnalysisTools)\n");
  fprintf(fw, "# command: ");
  PrintCommand(fw, argc, argv);
  fclose(fw);
} //}}}
void PrintStep(int *count_coor, const int start, const bool silent) { //{{{
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
void PrintLastStep(const int coor, const int used, const bool silent) { //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
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

// WriteFormatedDataLine //{{{
void WriteFormatedDataLine(FILE *fw, const int columns, const double *data,
                           const int (*digits)[2]) {
  for (int col = 0; col < columns; col++) {
    Fprintf1(fw, data[col], digits[col]);
  }
  putc('\n', fw);
} //}}}
