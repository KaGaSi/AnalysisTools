/**
 * \file
 * \brief Structures for utilities
 */

#ifndef _STRUCTS_H_
#define _STRUCTS_H_

#include "General.h"

#define CHARGE 10000.0 // 'impossible' charge to define a given bead type has charge specified in an input file
#define MASS 0.0 // 'impossible' mass to define a given bead type has charge specified in an input file
#define RADIUS 0.0 // 'impossible' radius to define a given bead type has charge specified in an input file
#define MOL_NAME 8 // maximum molecule name length (array size is +1)
#define BEAD_NAME 16 // maximum bead name length (array size is +1)

extern int *InFile;

// struct Counts //{{{
/**
 * \brief Total numbers of various things.
 */
typedef struct Counts {
  int TypesOfBeads, // number of bead types
      TypesOfMolecules, // number of molecule types
      TypesOfBonds, // number of bond types; -1 if not read from anywhere
      TypesOfAngles, // number of bond types; -1 if not read from anywhere
      TypesOfDihedrals, // number of dihedral types; -1 if not read from anywhere
      Bonded, // total number of beads in all molecules
      Unbonded, // total number of monomeric beads
      BeadsCoor, // number of beads in the coordinate file (e.g., in vcf file)
      BeadsTotal, // total number of beads in the system (e.g., in vsf file)
      Molecules, // total number of molecules
      HighestResid, // highest id in a file (discontinuous molecule counting)
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
  .Unbonded = 0,
  .BeadsCoor = 0,
  .BeadsTotal = 0,
  .Molecules = 0,
  .HighestResid = -1,
  .Aggregates = 0,
}; //}}}

// struct Box //{{{
/**
 * \brief Total numbers of various things.
 */
typedef struct Box {
  VECTOR Length; // side lengths (a, b, c for triclinic cell)
  VECTOR TriLength; // lx, ly, lz for lammps (triclinic cell)
  double TriTilt[3]; // xy, xz, and yz tilts for triclinic box
  double alpha, beta, gamma; // angles - all 90 for orthogonal box
  double transform[3][3], // transformation matrix
         inverse[3][3]; // inverse of the transformation matrix
  double Volume;
} BOX;

// TODO: is that used at all?
// Initialize Box
static const BOX InitBox = {
  .Length.x = -1,
  .Length.y = -1,
  .Length.z = -1,
  .TriLength.x = -1,
  .TriLength.y = -1,
  .TriLength.z = -1,
  .TriTilt = {0},
  .alpha = 90,
  .beta = 90,
  .gamma = 90,
  .Volume = -1,
}; //}}}

// struct Params //{{{
/**
 * \brief Parameters for various things (bonds and angles for now)
 */
typedef struct Params {
  double a, b;
} PARAMS;

// Initialize Params
static const PARAMS InitParams = {
  .a = -1,
  .b = -1,
}; //}}}

// struct BeadType //{{{
/**
 * \brief Information about bead types.
 */
typedef struct BeadType {
  char Name[BEAD_NAME+1]; // name of given bead type

  int Number; // number of beads of given type

  bool Use, // should bead type in .vcf file be used for calculation?
       Write; // should bead type in .vcf file be written to output .vcf?

  double Charge, // charge of every bead of given type
         Mass, // mass of every bead of given type
         Radius; // radius of every bead of the given type
} BEADTYPE; //}}}

// struct Bead //{{{
/**
 * \brief Information about every bead.
 */
typedef struct Bead {
  int Type, // type of bead corresponding to index in BeadType struct
      Molecule, // index number of molecule corresponding to Molecule struct (-1 for monomeric bead)
      nAggregates, // number of aggregates the bead is in (only monomeric beads can be in more aggregates - allocated memory for 10)
      *Aggregate, // index numbers of aggregates corresponding to Aggregate struct (-1 for bead in no aggregate)
      Index; // index of the bead according to .vsf file (needed for indexed timesteps)

  VECTOR Position; // cartesian coordinates of the bead
  VECTOR Velocity; // velocity of the bead

  bool InTimestep, // is the bead in the present timestep?
       Use; // general-purpose flag
  bool Flag; // general-purpose flag
} BEAD; //}}}

// struct MoleculeType //{{{
/**
 * \brief Information about molecule types.
 */
typedef struct MoleculeType {
  char Name[MOL_NAME+1]; // name of given molecule type

  int Number, // number of molecules of given type
      nBeads, // number of beads in every molecule of given type
      *Bead, // ids of bead types of every molecule bead
      nBonds, // number of bonds in every molecule of given type
      (*Bond)[3], // pair of ids for every bond (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Bond[3]
      nAngles, // number of angles in every molecule of given type
      (*Angle)[4], // trio of ids for every angle (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Angle[4]
      nDihedrals, // number of dihedrals in every molecule of given type
      (*Dihedral)[5], // fourtet of ids for every dihedral (with relative bead numbers from 0 to nBeads)
               // has to be sorted; size: [MoleculeType[i].Dihedral[5]
      nBTypes, // number of bead types in every molecule of given type
      *BType; // ids of bead types in every molecule of given type (corresponds to indices in BeadType struct)

  double Mass, // total mass of every molecule of given type
         Charge; // total charge of every molecule of given type

  bool InVcf, // is molecule type in vcf file?
       Use, // should molecule type be used for calculation?
       Write; // should molecule type be used for calculation?
} MOLECULETYPE; //}}}

// struct Molecule //{{{
/**
 * \brief Information about every molecule.
 */
typedef struct Molecule {
  int Type, // type of molecule corresponding to index in MoleculeType struct
      *Bead, // ids of beads in the molecule
      Aggregate, // id of aggregate molecule is in (corresponding to index in Aggregate struct)
      Index; // resid according to .vsf file
} MOLECULE; //}}}

// struct Aggregate //{{{
/**
 * \brief Information about every aggregate.
 */
typedef struct Aggregate {
  int nMolecules, // number of molecules in aggregate
      *Molecule, // ids of molecules in aggregate
      nBeads, // number of bonded beads in aggregate
      *Bead, // ids of bonded beads in aggregate
      nMonomers, // number of monomeric beads in aggregate
      *Monomer; // ids of monomeric beads in aggregate

  double Mass; // total mass of the aggregate

  bool Use; // should aggregate be used for calculation?
} AGGREGATE; //}}}
#endif
