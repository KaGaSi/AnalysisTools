#include "ReadWrite.h"
#include "Aggregates.h"
#include "Errors.h"
#include "General.h"
#include "Options.h"
#include "ReadWriteVtf.h"
#include "ReadWriteXyz.h"
#include "ReadWriteLtrj.h"
#include "ReadWriteLdata.h"
#include "ReadWriteField.h"
#include "ReadWriteConfig.h"
#include "ReadWriteGromacs.h"

// Format registry: one entry per supported file format.
// When adding a new format:
//   1. Add its enum value to Globals.h
//   2. Add a row here (that's all for detection/routing)
//   3. Add a case to ReadStructure() if it's a structure format
//   4. Add a case to ReadTimestep() / SkipTimestep() if it's a coordinate format
typedef struct {
  int         type;
  const char *extension;   // e.g. ".xyz"; NULL if matched only by basename
  const char *basename;    // e.g. "FIELD"; NULL if matched only by extension
  const char *string;      // primary name for FileTypeFromString / -ft flag
  const char *string_alt;  // secondary alias (e.g. "lammpstrj"); NULL if none
  bool        is_structure;
  bool        is_coordinate;
  bool        self_contained; // file serves as both struct + coor (no -i needed)
} FORMAT_INFO;
static const FORMAT_INFO FORMATS[] = {
  { VTF_FILE,    ".vtf",       NULL,     "vtf",    NULL,        true,  true,  true  },
  { VSF_FILE,    ".vsf",       NULL,     "vsf",    NULL,        true,  false, false },
  { VCF_FILE,    ".vcf",       NULL,     "vcf",    NULL,        false, true,  false },
  { XYZ_FILE,    ".xyz",       NULL,     "xyz",    NULL,        true,  true,  true  },
  { LDATA_FILE,  ".data",      NULL,     "data",   NULL,        true,  true,  true  },
  { LTRJ_FILE,   ".lammpstrj", NULL,     "ltrj",   "lammpstrj", true,  true,  true  },
  { FIELD_FILE,  ".field",     "FIELD",  "field",  NULL,        true,  false, false },
  { CONFIG_FILE, ".config",    "CONFIG", "config", NULL,        false, true,  false },
  { ITP_FILE,    ".itp",       NULL,     "itp",    NULL,        true,  false, false },
  { PDB_FILE,    ".pdb",       NULL,     "pdb",    NULL,        true,  false, false },
};
static const int N_FORMATS = (int)(sizeof FORMATS / sizeof *FORMATS);

static void CopyAndFreeStuff(const int n, int (**old)[5], int (**new)[5]);
static void CopyAndFreeAllStuff(MOLECULETYPE *mt_old, MOLECULETYPE *mt_new);
static void MinimizeOneMtypeStuffIds(const int num, int (**arr)[5],
                                     const int n, const int *link_bead_ids);
static void FillMTypeStuff(SYSTEM *System, const int type, const int size,
                           const int (*all)[5], const int n);
static int FindFileType(const char *name);
static bool ReadAggMolLine(FILE *fr, const char *file, const int line_count,
                           const int n_mol, const int *resid_to_mol,
                           const int max_index, int *n, int **arr);

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
    if (mt_i->nBeads == 1 && !mt_i->Named) { // remove 'fake' molecules
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
    if (!bt_old_to_new) {
      ErrorAlloc("bt_old_to_new");
    }
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
    fprintf(fw, "pbc %.3f %.3f %.3f", box.Length.x,
                                      box.Length.y,
                                      box.Length.z);
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
      err_msg("vcf file without data");
      PrintErrorFile(f.coor.name, "\0", "\0");
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
    case LTRJ_FILE:
      if (LtrjReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case LDATA_FILE:
      if (LmpDataReadTimestep(fr, f.coor.name, System, line_count) < 0) {
        return false;
      }
      break;
    case CONFIG_FILE:
      err_msg("CONFIG format reading not yet implemented");
      PrintError();
      exit(1);
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
      if (XyzSkipTimestep(fr, f.coor.name, line_count) < 0) {
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
    case CONFIG_FILE:
      err_msg("CONFIG format reading not yet implemented");
      PrintError();
      exit(1);
    default:
      snprintf(ERROR_MSG, LINE, "no action specified for coor_type %s%d",
               ErrYellow(), f.coor.type);
      PrintError();
      exit(1);
  }
  return true;
} //}}}
// Read aggregates from a single timestep //{{{
/*
 * Returns 1 on success, -1 on end of data (Last Step line, or file ending at
 * a step boundary, e.g., when the Aggregates run was interrupted), and -2 on
 * malformed data (mid-record end of file, invalid counts or molecule ids).
 */
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
    // end of file at a step boundary: truncated but consistent data
    return -1;
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
  // build resid -> compact molecule index lookup table
  int max_index = 0;
  for (int i = 0; i < Count->Molecule; i++) {
    if (System->Molecule[i].Index > max_index) {
      max_index = System->Molecule[i].Index;
    }
  }
  int *resid_to_mol = malloc((max_index + 1) * sizeof *resid_to_mol);
  if (!resid_to_mol) {
    ErrorAlloc("resid_to_mol");
  }
  for (int i = 0; i <= max_index; i++) {
    resid_to_mol[i] = -1;
  }
  for (int i = 0; i < Count->Molecule; i++) {
    int idx = System->Molecule[i].Index;
    if (idx >= 0 && idx <= max_index) {
      resid_to_mol[idx] = i;
    }
  }
  // go through all aggregates, reading core and border molecule lines
  for (int i = 0; i < Count->Aggregate; i++) {
    AGGREGATE *Agg = &Aggregate[i];
    // line 1: core molecules
    (*line_count)++;
    if (!ReadAggMolLine(fr, file, *line_count, Count->Molecule, resid_to_mol,
                        max_index, &Agg->nCore, &Agg->Core)) {
      free(resid_to_mol);
      return -2;
    }
    // line 2: border molecules
    (*line_count)++;
    if (!ReadAggMolLine(fr, file, *line_count, Count->Molecule, resid_to_mol,
                        max_index, &Agg->nBorder, &Agg->Border)) {
      free(resid_to_mol);
      return -2;
    }
    Agg->nMolecules = Agg->nCore + Agg->nBorder;
  }
  free(resid_to_mol);
  // fill the rest of aggregate info
  for (int i = 0; i < Count->Aggregate; i++) {
    int count = 0;
    double mass = 0;
    AGGREGATE *Agg = &Aggregate[i];
    for (int j = 0; j < Agg->nMolecules; j++) {
      MOLECULE *mol = &System->Molecule[AggGetMol(Agg, j)];
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
// Read and validate single aggregate line //{{{
static bool ReadAggMolLine(FILE *fr, const char *file, const int line_count,
                           const int n_mol, const int *resid_to_mol,
                           const int max_index, int *n, int **arr) {
  if (fscanf(fr, "%d :", n) != 1 || *n < 0 || *n > n_mol) {
    snprintf(ERROR_MSG, LINE, "invalid number of aggregate molecules; "
             "stopping reading %s%s%s", ErrYellow(), file, ErrCyan());
    PrintWarnFileLine(file, line_count);
    return false;
  }
  int all = *n;
  if (all <= 0) {
    all = 1;
  }
  *arr = s_realloc(*arr, all * sizeof **arr);
  for (int j = 0; j < *n; j++) {
    int resid;
    if (fscanf(fr, "%d", &resid) != 1) {
      snprintf(ERROR_MSG, LINE, "incomplete aggregate molecule line; "
               "stopping reading %s%s%s", ErrYellow(), file, ErrCyan());
      PrintWarnFileLine(file, line_count);
      return false;
    }
    if (resid < 0 || resid > max_index || resid_to_mol[resid] == -1) {
      snprintf(ERROR_MSG, LINE, "molecule id %s%d%s not present in the "
               "system; stopping reading %s%s%s", ErrYellow(), resid,
               ErrCyan(), ErrYellow(), file, ErrCyan());
      PrintWarnFileLine(file, line_count);
      return false;
    }
    (*arr)[j] = resid_to_mol[resid];
  }
  // skip the rest of the line
  int ch;
  while ((ch = getc(fr)) != '\n' && ch != EOF)
    ;
  return true;
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
  // skip two lines per aggregate (core line + border line)
  for (int i = 0; i < num_agg; i++) {
    for (int l = 0; l < 2; l++) {
      (*line_count)++;
      int ch;
      while ((ch = getc(fr)) != '\n' && ch != EOF)
        ;
    }
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
                     const SYSTEM System, const AGGREGATE *Aggregate,
                     const bool *use_agg) {
  // get number of aggregates to write to agg_file
  int number_of_aggs = 0;
  for (int i = 0; i < System.Count.Aggregate; i++) {
    if (use_agg[i]) {
      number_of_aggs++;
    }
  }
  FILE *fw = OpenFile(agg_file, "a");
  // print number of aggregates to agg file
  fprintf(fw, "Step: %d\n%d\n", step_count, number_of_aggs);
  // go through all aggregates
  for (int i = 0; i < System.Count.Aggregate; i++) {
    // write only those that aren't excluded
    if (use_agg[i]) {
      // line 1: core molecules
      fprintf(fw, "%d :", Aggregate[i].nCore);
      for (int j = 0; j < Aggregate[i].nCore; j++) {
        fprintf(fw, " %d", System.Molecule[Aggregate[i].Core[j]].Index);
      }
      putc('\n', fw);
      // line 2: border molecules
      fprintf(fw, "%d :", Aggregate[i].nBorder);
      for (int j = 0; j < Aggregate[i].nBorder; j++) {
        fprintf(fw, " %d", System.Molecule[Aggregate[i].Border[j]].Index);
      }
      putc('\n', fw);
    }
  }
  fclose(fw);
} //}}}

void PrintByline(const char *file, const int argc, char **argv) { //{{{
  FILE *fw = OpenFile(file, "w");
  fprintf(fw, "# Created by AnalysisTools v%s ", VERSION);
  fprintf(fw, " (https://github.com/KaGaSi/AnalysisTools)\n");
  fprintf(fw, "# Command: ");
  PrintCommand(fw, argc, argv);
  fclose(fw);
} //}}}
FILE * PrintBylineOpenFile(const char *f, const int argc, char **argv) { //{{{
  PrintByline(f, argc, argv);
  FILE *ptr = OpenFile(f, "a");
  return ptr;
} //}}}
// file type detection
static int FindFileType(const char *name) { //{{{
  // a) basename match (e.g. "FIELD", "CONFIG")
  for (int i = 0; i < N_FORMATS; i++) {
    if (FORMATS[i].basename && strcasecmp(name, FORMATS[i].basename) == 0) {
      return FORMATS[i].type;
    }
  }
  // b) extension match
  const char *dot = strrchr(name, '.');
  if (dot) {
    for (int i = 0; i < N_FORMATS; i++) {
      if (FORMATS[i].extension && strcasecmp(dot, FORMATS[i].extension) == 0) {
        return FORMATS[i].type;
      }
    }
  }
  return -1;
} //}}}
int FileTypeFromString(const char *str) { //{{{
  for (int i = 0; i < N_FORMATS; i++) {
    if (strcasecmp(str, FORMATS[i].string) == 0) {
      return FORMATS[i].type;
    }
    if (FORMATS[i].string_alt && strcasecmp(str, FORMATS[i].string_alt) == 0) {
      return FORMATS[i].type;
    }
  }
  return -1;
} //}}}
// identify input coordinate and structure files //{{{
// TODO: no return false - so why bool?
bool InputCoorStruct(const int argc, char **argv, SYS_FILES *f) {
  // -ft option: override file type without extension-based detection
  int ft_override = -1;
  char ft_str[LINE];
  if (FileOption(argc, argv, COMMON_OPTS[C_FT].opt, ft_str)) {
    ft_override = FileTypeFromString(ft_str);
    if (ft_override == -1 || ft_override == VSF_FILE) {
      if (snprintf(ERROR_MSG, LINE, "unknown coordinate file type '%s'",
                   ft_str) < 0) {
        ErrorSnprintf();
      }
      PrintErrorOption(COMMON_OPTS[C_FT].opt);
      exit(1);
    }
  }
  // input structure file (-i option) with optional type string as 2nd argument
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      if ((i+1) >= argc || argv[i+1][0] == '-') {
        s_strcpy(ERROR_MSG,
                 "missing file name (or file name begins with a dash)", LINE);
        PrintErrorOption("-i");
        exit(1);
      }
      s_strcpy(f->stru.name, argv[i+1], LINE);
      // optional type string: consume only if it's a recognized structure type
      if ((i+2) < argc && argv[i+2][0] != '-') {
        int ft = FileTypeFromString(argv[i+2]);
        bool is_struct = false;
        for (int j = 0; j < N_FORMATS; j++) {
          if (FORMATS[j].type == ft && FORMATS[j].is_structure) {
            is_struct = true;
            break;
          }
        }
        if (is_struct) {
          f->stru.type = ft;
        } else if (ft != -1) {
          if (snprintf(ERROR_MSG, LINE, "not a structure file type: '%s%s%s'",
                       ErrYellow(), argv[i+2], ErrRed()) < 0) {
            ErrorSnprintf();
          }
          PrintErrorOption("-i");
          exit(1);
        } else {
          // not a format string or a type hint, detect from filename
          f->stru.type = StructureFileType(f->stru.name);
        }
      } else {
        f->stru.type = StructureFileType(f->stru.name);
      }
      break;
    }
  }
  if (ft_override != -1) {
    f->coor.type = ft_override;
  } else {
    f->coor.type = CoordinateFileType(f->coor.name);
  }
  // set default structure file if -i option not used
  if (f->stru.name[0] == '\0') {
    if (f->coor.type == VCF_FILE) { // use vcf file with .vsf ending
      int last = -1;
      for (int i = 0; i < strnlen(f->coor.name, LINE); i++) {
        if (f->coor.name[i] == '.') {
          last = i;
        }
      }
      s_strcpy(f->stru.name, f->coor.name, LINE);
      f->stru.name[last+2] = 's';
      f->stru.type = VSF_FILE;
    } else {
      bool self_contained = false;
      for (int j = 0; j < N_FORMATS; j++) {
        if (FORMATS[j].type == f->coor.type && FORMATS[j].self_contained) {
          self_contained = true;
          break;
        }
      }
      if (self_contained) {
        s_strcpy(f->stru.name, f->coor.name, LINE);
        f->stru.type = f->coor.type;
      } else {
        err_msg("missing structure file; should never happen!");
        PrintError();
        exit(1);
      }
    }
  }
  return true;
} //}}}
int StructureFileType(const char *path) { //{{{
  const char *name = StripPath(path);
  int ft = FindFileType(name);
  for (int i = 0; i < N_FORMATS; i++) {
    if (FORMATS[i].type == ft && FORMATS[i].is_structure) {
      return ft;
    }
  }
  err_msg("Not a structure file");
  PrintErrorFile(path, "\0", "\0");
  exit(1);
} //}}}
int CoordinateFileType(const char *path) { //{{{
  const char *name = StripPath(path);
  int ft = FindFileType(name);
  for (int i = 0; i < N_FORMATS; i++) {
    if (FORMATS[i].type == ft && FORMATS[i].is_coordinate) {
      return ft;
    }
  }
  err_msg("Not a coordinate file");
  PrintErrorFile(path, "\0", "\0");
  exit(1);
} //}}}
int FileType(const char *name) { //{{{
  int ft = FindFileType(name);
  if (ft != -1) {
    return ft;
  } else {
    err_msg("Unknown file type");
    PrintErrorFile(name, "\0", "\0");
    exit(1);
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
