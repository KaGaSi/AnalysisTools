#ifndef STRUCTS_H
#define STRUCTS_H

#define _POSIX_C_SOURCE 200809L

#include "General.h"

#define VERSION "4.0"
#define DATE "TBD"
#define OPT_LENGTH 15

// structures for options //{{{
typedef struct OPT OPT;
OPT * opt_create(void);
typedef struct common_opt {
  bool verbose, silent;
  int start, end, skip;
} COMMON_OPT; //}}}
// structures for file names and types //{{{
typedef struct file_type {
  char name[LINE];
  int type;
} FILE_TYPE;
static const FILE_TYPE InitFile = {
  .name = "",
  .type = -1,
};
typedef struct sys_files {
  FILE_TYPE coor, stru;
} SYS_FILES;
static const SYS_FILES InitSysFiles = {
  .coor = InitFile,
  .stru = InitFile,
}; //}}}
typedef struct Box { //{{{
  double Length[3],
         OrthoLength[3],
         Bounding[3],
         Low[3],
         alpha, beta, gamma, // angles - all 90 for orthogonal box
         transform[3][3], // transformation matrix
         inverse[3][3], // inverse of the transformation matrix
         Volume;
} BOX;
// Initialize Box
static const BOX InitBox = {
  .Length = {-1, -1, -1},
  .OrthoLength = {-1, -1, -1},
  .Bounding = {-1, -1, -1},
  .Low = {0, 0, 0},
  .alpha = 90,
  .beta = 90,
  .gamma = 90,
  .transform = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
  .inverse = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
  .Volume = -1,
}; //}}}
typedef struct Count { //{{{
  int BeadType, MoleculeType,
      BondType, AngleType, DihedralType, ImproperType,
      Bead, // total number of beads in the system (e.g., in vsf file)
      BeadCoor, // number of beads in the coordinate file (e.g., in vcf file)
      Bonded, // total number of beads in all molecules
      BondedCoor, // total number of bonded beads in in the coordinate file
      Unbonded, // total number of monomeric beads
      UnbondedCoor, // total number of monomeric beads in the coordinate file
      Molecule, // total number of molecules
      MoleculeCoor, // total number of molecules in the coordinate file
      Aggregate, // number of aggregates
      HighestResid, // highest id in a file (discontinuous molecule counting)
      Bond, Angle, Dihedral, Improper;
} COUNT;
// Initialize Count
static const COUNT InitCount = {
  .BeadType = 0,
  .MoleculeType = 0,
  .BondType = 0,
  .AngleType = 0,
  .DihedralType = 0,
  .ImproperType = 0,
  .Bead = 0,
  .BeadCoor = 0,
  .Bonded = 0,
  .BondedCoor = 0,
  .Unbonded = 0,
  .UnbondedCoor = 0,
  .Molecule = 0,
  .MoleculeCoor = 0,
  .Aggregate = 0,
  .HighestResid = -1,
  .Bond = 0,
  .Angle = 0,
  .Dihedral = 0,
  .Improper = 0,
}; //}}}
typedef struct Params { //{{{
  double a, b, c, d;
} PARAMS;
static const PARAMS InitParams = {
  .a = 0,
  .b = 0,
  .c = 0,
  .d = 0,
}; //}}}
typedef struct BeadType { //{{{
  char Name[BEAD_NAME]; // name of given bead type
  int Number, // number of beads of given type
      InCoor, // nubmer of beads in the coordinate file
      *Index; // array of Bead[] indices
  double Charge, // charge of every bead of given type
         Mass, // mass of every bead of given type
         Radius; // radius of every bead of the given type
  bool Flag; // general-purpose flag
} BEADTYPE;
void InitBeadType(BEADTYPE *bt); //}}}
typedef struct Bead { //{{{
  int Type, // type of bead corresponding to index in BeadType struct
      Molecule, // id corresponding to Molecule struct (-1 for monomeric bead)
      Aggregate; // aggregate id the molecule is in (-1 for none)
  double Position[3], Velocity[3], Force[3];
  // VECTOR Position; // cartesian coordinates of the bead
         // Velocity; // velocity of the bead
         // Force; // force acting on the bead
  bool InTimestep; // is the bead in the present timestep?
  bool Flag; // general-purpose flag; TODO: remove?
} BEAD;
void InitBead(BEAD *b); //}}}
typedef struct MoleculeType { //{{{
  char Name[MOL_NAME]; // name of given molecule type

  int Number, // number of molecules of given type
      *Index, // array of Molecule[] indices
      nBeads, // number of beads in every molecule of given type
      *Bead, // ids of bead types of every molecule bead
      nBonds, // number of bonds in every molecule of given type
      (*Bond)[5], // pair of ids for every bond (with relative bead numbers from 0 to nBeads-1)
                  // has to be sorted; size: [MoleculeType[i].Bond[3]
      nAngles, // number of angles in every molecule of given type
      (*Angle)[5], // trio of ids for every angle (with relative bead numbers from 0 to nBeads-1)
                  // has to be sorted; size: [MoleculeType[i].Angle[4]
      nDihedrals, // number of dihedrals in every molecule of given type
      (*Dihedral)[5], // fourtet of ids for every dihedral (with relative bead numbers from 0 to nBeads-1)
      nImpropers, // number of improper dihedrals in every molecule of given type
      (*Improper)[5], // fourtet of ids for every improper dihedral (with relative bead numbers from 0 to nBeads-1)
                      // has to be sorted; size: [MoleculeType[i].Improper[5]
      nBTypes, // number of bead types in every molecule of given type
      *BType; // ids of bead types in every molecule of given type (corresponds to indices in BeadType struct)

  double Mass, // total mass of every molecule of given type
         Charge; // total charge of every molecule of given type

  bool InVcf, // is molecule type in vcf file? TODO: useless?
       Flag; // general-purpose flag
} MOLECULETYPE;
void InitMoleculeType(MOLECULETYPE *mt); //}}}
typedef struct Molecule { //{{{
  int Type, // type of molecule corresponding to index in MoleculeType struct
      *Bead, // ids of beads in the molecule
      Index, // resid according to input file
      Aggregate; // aggregate id the molecule is in (-1 for none)
  bool InTimestep; // is the molecule in a timestep?
} MOLECULE;
void InitMolecule(MOLECULE *mol); //}}}
typedef struct System { //{{{
  BOX Box;
  COUNT Count;
  BEADTYPE *BeadType;
  BEAD *Bead;
  MOLECULETYPE *MoleculeType;
  MOLECULE *Molecule;
  PARAMS *BondType;
  PARAMS *AngleType;
  PARAMS *DihedralType;
  PARAMS *ImproperType;
  int *MoleculeCoor, // array of molecule ids with InTimestep=true
      *Bonded, // array of Bead[] ids of in-molecule beads
      *BondedCoor, // array of Bead[] ids of in-molecule beads in a timestep
      *Unbonded, // array of Bead[] ids of in-molecule beads
      *UnbondedCoor, // array of Bead[] ids of in-molecule beads in a timestep
      *BeadCoor; // array of internal ids for beads with InTimestep=true
} SYSTEM;
void InitSystem(SYSTEM *System); //}}}
typedef struct Aggregate { //{{{
  int nMolecules, // number of molecules in aggregate
      *Molecule, // ids of molecules in aggregate
      nBeads, // number of bonded beads in aggregate
      *Bead; // ids of bonded beads in aggregate
  double Mass; // total mass of the aggregate
  bool Flag; // should aggregate be used for calculation?
} AGGREGATE;
void InitAggregate(SYSTEM System, AGGREGATE **Aggregate);
void ReInitAggregate(SYSTEM System, AGGREGATE *Aggregate); //}}}
#endif
