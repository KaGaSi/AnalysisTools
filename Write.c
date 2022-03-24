#include "Write.h"

// WriteCoorIndexed_old() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Write flag
 * in BeadType and MoleculeType structures only certain bead types will be
 * saved into the indexed timestep in .vcf file.
 */
void WriteCoorIndexed_old(FILE *vcf_file, COUNTS Counts,
                      BEADTYPE *BeadType, BEAD *Bead,
                      MOLECULETYPE *MoleculeType, MOLECULE *Molecule,
                      char *stuff, VECTOR BoxLength) {
  // print comment at the beginning of a timestep if present in initial vcf file
  fprintf(vcf_file, "%s\n", stuff);
  // print box size
  fprintf(vcf_file, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  // print 'indexed' on the next
  fprintf(vcf_file, "indexed\n");

  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int btype = Bead[i].Type;
    if (BeadType[btype].Write) {
      if (Bead[i].Molecule != -1) { // bead in a molecule
        int mtype = Molecule[Bead[i].Molecule].Type;
        if (MoleculeType[mtype].Write) {
          fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", Bead[i].Index,
                                                       Bead[i].Position.x,
                                                       Bead[i].Position.y,
                                                       Bead[i].Position.z);
        }
      } else { // monomer bead
        fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", Bead[i].Index,
                                                     Bead[i].Position.x,
                                                     Bead[i].Position.y,
                                                     Bead[i].Position.z);
      }
    }
  }
} //}}}

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
          fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", Bead[i].Index,
                                                       Bead[i].Position.x,
                                                       Bead[i].Position.y,
                                                       Bead[i].Position.z);
        }
      } else { // monomer bead
        fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", Bead[i].Index,
                                                     Bead[i].Position.x,
                                                     Bead[i].Position.y,
                                                     Bead[i].Position.z);
      }
    }
  }
} //}}}

// VtfWriteCoorIndexed() //{{{
/**
 * Function writing coordinates to a `.vcf` file. According to the Write flag
 * in BeadType and MoleculeType structures only certain bead types will be
 * saved into the indexed timestep in .vcf file.
 */
void VtfWriteCoorIndexed(FILE *vcf_file, char *stuff,
                         COUNTS Counts, BEAD *Bead, BOX Box) {
  // print comment at the beginning of a timestep if present in initial vcf file
  fprintf(vcf_file, "%s\n", stuff);
  // print box size
  fprintf(vcf_file, "pbc %lf %lf %lf  ", Box.Length.x,
                                         Box.Length.y,
                                         Box.Length.z);
  fprintf(vcf_file, "    %lf %lf %lf\n", Box.alpha, Box.beta, Box.gamma);
  // print 'indexed' on the next
  fprintf(vcf_file, "indexed\n");

//for (int i = 0; i < Counts.BeadsTotal; i++) {
//  int id = i;
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int id = InFile[i];
    if (Bead[id].InTimestep && Bead[id].Use) {
      fprintf(vcf_file, "%8d %8.4f %8.4f %8.4f\n", Bead[id].Index,
                                                   Bead[id].Position.x,
                                                   Bead[id].Position.y,
                                                   Bead[id].Position.z);
    }
  }
} //}}}

// WriteCoorXYZ() //{{{
void WriteCoorXYZ(FILE *xyz_file, COUNTS Counts,
                  BEADTYPE *BeadType, BEAD *Bead) {

  // count beads to write
  int count = 0;
  for (int i = 0; i < Counts.TypesOfBeads; i++) {
    if (BeadType[i].Write) {
      count += BeadType[i].Number;
    }
  }

  // print number of beads to file
  fprintf(xyz_file, "%d\n\n", count);

  // print coordinates
  for (int i = 0; i < Counts.BeadsCoor; i++) {
    int type = Bead[i].Type;
    if (BeadType[type].Write) {
      fprintf(xyz_file, "%8s %7.3f %7.3f %7.3f\n", BeadType[type].Name, Bead[i].Position.x, Bead[i].Position.y, Bead[i].Position.z);
    }
  }
} //}}}

// TODO: if q<0, are the lines nicely lined up?
// WriteVsf() //{{{
/*
 * Function creating `.vsf` structure file for use in conjunction with
 * `.vcf` coordinate file for better visualisation via VMD program. Note that
 * if not all bead types from the original system are not present in the
 * structures, the system may significantly differ from the original one.
 */
void WriteVsf(char *input_vsf, COUNTS Counts, BEADTYPE *BeadType, BEAD *Bead,
              MOLECULETYPE *MoleculeType, MOLECULE *Molecule, bool change) {

  // opten structure file //{{{
  FILE *fw;
  if ((fw = fopen(input_vsf, "w")) == NULL) {
    ErrorFileOpen(input_vsf, 'w');
    exit(1);
  } //}}}
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
      fprintf(fw, "atom %7d ", Bead[i].Index);
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
      fprintf(fw, "bond %6d: %6d\n", Bead[bead1].Index,
                                     Bead[bead2].Index);
    }
  } //}}}

  // close structure file
  fclose(fw);
} //}}}

// WriteAggregates() //{{{
/**
 * Function writiing (appending to .agg file) information about aggregates from
 * given timestep.
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

  // open .agg file for appending //{{{
  FILE *fw;
  if ((fw = fopen(agg_file, "a")) == NULL) {
    ErrorFileOpen(agg_file, 'a');
    exit(1);
  } //}}}

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
        fprintf(fw, " %d", Bead[Aggregate[i].Monomer[j]].Index);
      }
      putc('\n', fw);
    }
  } //}}}

  fclose(fw);
} //}}}

// WriteField() //{{{
/*
 * Function writing a dl_meso FIELD file.
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

  // opten FIELD file //{{{
  FILE *fw;
  if ((fw = fopen(field, "w")) == NULL) {
    ErrorFileOpen(field, 'w');
    exit(1);
  } //}}}
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
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourText(STDERR_FILENO, RED);
      fprintf(stderr, " - molecule ");
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", MoleculeType[i].Name);
      ColourText(STDERR_FILENO, RED);
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
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", field);
      ColourText(STDERR_FILENO, RED);
      fprintf(stderr, " - molecule ");
      ColourText(STDERR_FILENO, YELLOW);
      fprintf(stderr, "%s", MoleculeType[i].Name);
      ColourText(STDERR_FILENO, RED);
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

// PrintByline() //{{{
void PrintByline(FILE *ptr, int argc, char *argv[]) {
  fprintf(ptr, "# Created by AnalysisTools v%s ", VERSION);
  fprintf(ptr, " (https://github.com/KaGaSi/AnalysisTools)\n");
  fprintf(ptr, "# command: ");
  PrintCommand(ptr, argc, argv);
} //}}}
