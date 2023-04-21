#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#include "General.h"

#define CHARGE 10000.0 // 'impossible' charge
#define MASS 0.0 // 'impossible' mass
#define RADIUS 0.0 // 'impossible' radius
#define MOL_NAME 9 // maximum molecule name length (with null terminator)
#define BEAD_NAME 17 // maximum bead name length (with null terminator)

typedef struct Box { //{{{
  VECTOR Length, // side lengths (a, b, c for triclinic cell)
         OrthoLength, // orthogonal length (lx, ly, lz in lammps speak)
         Bounding; // maxium orthogonal length (x_bound, etc. in lammps speak)
  double alpha, beta, gamma, // angles - all 90 for orthogonal box
         transform[3][3], // transformation matrix
         inverse[3][3], // inverse of the transformation matrix
         Volume;
} BOX;
// Initialize Box
static const BOX InitBox = {
  .Length = {-1, -1, -1},
  .OrthoLength = {-1, -1, -1},
  .Bounding = {-1, -1, -1},
  .alpha = 90,
  .beta = 90,
  .gamma = 90,
  .transform = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
  .inverse = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
  .Volume = -1,
}; //}}}
typedef struct Count { //{{{
  int BeadType, // number of bead types
      MoleculeType, // number of molecule types
      BondType, // number of bond types; -1 if not read from anywhere
      AngleType, // number of bond types; -1 if not read from anywhere
      DihedralType, // number of dihedral types; -1 if not read from anywhere
      ImproperType, // number of dihedral types; -1 if not read from anywhere
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
      Bond,
      Angle,
      Dihedral,
      Improper;
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
  double a, b, c;
} PARAMS;
static const PARAMS InitParams = {
  .a = 0,
  .b = 0,
  .c = 0,
}; //}}}
typedef struct BeadType { //{{{
  char Name[BEAD_NAME]; // name of given bead type
  int Number, // number of beads of given type
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

  VECTOR Position, // cartesian coordinates of the bead
         Velocity, // velocity of the bead
         Force; // force acting on the bead

  bool InTimestep, // is the bead in the present timestep?
       Flag; // general-purpose flag
} BEAD;
void InitBead(BEAD *b); //}}}
typedef struct MoleculeType { //{{{
  char Name[MOL_NAME]; // name of given molecule type

  int Number, // number of molecules of given type
      *Index, // array of Molecule[] indices
      nBeads, // number of beads in every molecule of given type
      *Bead, // ids of bead types of every molecule bead
      nBonds, // number of bonds in every molecule of given type
      (*Bond)[3], // pair of ids for every bond (with relative bead numbers from 0 to nBeads-1)
                  // has to be sorted; size: [MoleculeType[i].Bond[3]
      nAngles, // number of angles in every molecule of given type
      (*Angle)[4], // trio of ids for every angle (with relative bead numbers from 0 to nBeads-1)
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
  int *Index_mol, // link between indices (i.e., Index_mol[Molecule[i].Index]=i)
      *Bonded, // array of Bead[] ids of in-molecule beads
      *BondedCoor, // array of Bead[] ids of in-molecule beads in a timestep
      *Unbonded, // array of Bead[] ids of in-molecule beads
      *UnbondedCoor, // array of Bead[] ids of in-molecule beads in a timestep
      *BeadCoor; // array of internal ids for beads with InTimestep=true
} SYSTEM;
void InitSystem(SYSTEM *System); //}}}

// TODO remake/remove?
typedef struct Aggregate { //{{{
  int nMolecules, // number of molecules in aggregate
      *Molecule, // ids of molecules in aggregate
      nBeads, // number of bonded beads in aggregate
      *Bead, // ids of bonded beads in aggregate
      nMonomers, // that's probably gonna go
      *Monomer; // that's probably gonna go

  double Mass; // total mass of the aggregate

  bool Use; // should aggregate be used for calculation?
} AGGREGATE; //}}}
#if 0
// TODO remove
typedef struct Counts { //{{{
  int TypesOfBeads, // number of bead types
      TypesOfMolecules, // number of molecule types
      TypesOfBonds, // number of bond types; -1 if not read from anywhere
      TypesOfAngles, // number of bond types; -1 if not read from anywhere
      TypesOfDihedrals, // number of dihedral types; -1 if not read from anywhere
      Bonded, // total number of beads in all molecules
      BondedCoor, // total number of bonded beads in in the coordinate file
      Unbonded, // total number of monomeric beads
      UnbondedCoor, // total number of monomeric beads in the coordinate file
      BeadsTotal, // total number of beads in the system (e.g., in vsf file)
      BeadsCoor, // number of beads in the coordinate file (e.g., in vcf file)
      Molecules, // total number of molecules
      HighestResid, // highest id in a file (discontinuous molecule counting)
                    // TODO is that useful outside reading vsf?
      Aggregates; // total number of aggregates
} COUNTS;

// Initialize Counts
static const COUNTS InitCounts = {
  .TypesOfBeads = 0,
  .TypesOfMolecules = 0,
  .TypesOfBonds = -1,
  .TypesOfAngles = -1,
  .TypesOfDihedrals = -1,
  .Bonded = 0,
  .BondedCoor = 0,
  .Unbonded = 0,
  .UnbondedCoor = 0,
  .BeadsTotal = 0,
  .BeadsCoor = 0,
  .Molecules = 0,
  .HighestResid = -1,
  .Aggregates = 0,
}; //}}}
#endif
#endif
