#include "Aggregates.h"
#include "Coordinates.h"
#include "Errors.h"
#include "General.h"
#include "MathUtils.h"
#include "Options.h"
#include "PairHash.h"

// TODO: AggPickerOptions should be in Options.c, no?

// STATIC DECLARATIONS
static int NewAgg(AGGREGATE *Aggregate, SYSTEM *System,
                  const int i, const int j);

// evaluate bead contacts to assign molecules to aggregates //{{{
void EvaluateContacts(AGGREGATE *Aggregate, SYSTEM *System,
                      const int contacts, PairHash *contact) {
  COUNT *Count = &System->Count;
  /*
   * iterate only over the entries stored in the hash table, i.e., pairs that
   * have at least one bead-bead contact detected.
   *
   * khash iteration pattern:
   *   kh_begin(h) / kh_end(h) ... first / one-past-last bucket index
   *   kh_exist(h, it) ... true if bucket 'it' holds a live entry
   *   kh_key(h, it) ... the uint64_t key at 'it'
   *   kh_val(h, it) ... the uint8_t contact count at 'it'
   *
   * Each key encodes a molecule pair (i, j) with i > j; use the
   * PairHash_mol_i / PairHash_mol_j helpers to decode them.
   */
  for (khiter_t it = kh_begin(contact); it != kh_end(contact); ++it) {
    // skip empty buckets (open-addressing tables have gaps)
    if (!kh_exist(contact, it)) {
      continue;
    }
    // decode the pair
    uint64_t key = kh_key(contact, it);
    int i = PairHash_mol_i(key); // larger mol index
    int j = PairHash_mol_j(key); // smaller mol index
    if (System->Molecule[i].InTimestep && System->Molecule[j].InTimestep) {
      int agg_i = System->Molecule[i].Aggregate,
          agg_j = System->Molecule[j].Aggregate;
      // if molecules 'i' and 'j' are in contact, put them into one aggregate
      if (kh_val(contact, it) >= (uint8_t)contacts) { //{{{
          // create new aggregate if molecule 'j' isn'it in any
          if (agg_j == -1) {
            agg_j = NewAgg(Aggregate, System, i, j);
          }
          /*
           * add molecule 'i' to aggregate 'j' ()
           * if molecule 'i' isn't in any aggregate
           */
          if (agg_i == -1) {
            int mols = Aggregate[agg_j].nMolecules;
            AGGREGATE *Agg = &Aggregate[agg_j];
            Agg->nMolecules++;
            Agg->Molecule = s_realloc(Agg->Molecule,
                                      Agg->nMolecules * sizeof *Agg->Molecule);
            Agg->Molecule[mols] = i;
            System->Molecule[i].Aggregate = agg_j;
          }
          /*
           * if molecules 'i' and 'j' are in different aggregates,
           * unite those aggregates
           */
          if (agg_i != -1 && agg_j != -1 && agg_i != agg_j) {
            // add molecules from aggregate 'i' to aggregate 'j'
            int n_mol_old = Aggregate[agg_j].nMolecules;
            AGGREGATE *Agg_i = &Aggregate[agg_i];
            AGGREGATE *Agg_j = &Aggregate[agg_j];
            Agg_j->nMolecules += Agg_i->nMolecules;
            Agg_j->Molecule = s_realloc(Agg_j->Molecule,
                                      Agg_j->nMolecules * sizeof *Agg_j->Molecule);
            for (int k = n_mol_old; k < Agg_j->nMolecules; k++) {
              int mol = Agg_i->Molecule[k-n_mol_old];
              Agg_j->Molecule[k] = mol;
              System->Molecule[mol].Aggregate = agg_j;
            }
            // move aggregates with id greater then agg_i to id-1
            for (int k = (agg_i + 1); k < Count->Aggregate; k++) {
              AGGREGATE *Agg_k = &Aggregate[k];
              AGGREGATE *Agg_k1 = &Aggregate[k-1];
              Agg_k1->nMolecules = Agg_k->nMolecules;
              Agg_k1->Molecule = s_realloc(Agg_k1->Molecule,
                                           Agg_k1->nMolecules *
                                           sizeof *Agg_k1->Molecule);
              // move every molecule from aggregate 'k' to aggregate 'k-1'
              for (int l = 0; l < Agg_k->nMolecules; l++) {
                int mol = Agg_k->Molecule[l];
                Agg_k1->Molecule[l] = mol;
                System->Molecule[mol].Aggregate = k - 1;
              }
            }
            // reduce number of aggregates since the two aggregates were merged
            Count->Aggregate--;
          } //}}}
        /*
         * or if molecules 'i' and 'j' aren't in contact, and molecule 'j' isn't
         * in any aggregate, create new aggregate for molecule 'j'
         */
        } else if (agg_j == -1) { //{{{
          NewAgg(Aggregate, System, i, j);
        } //}}}
    }
  }
  // single-molecule aggregates (hash table ignores 0-contact pairs)
  for (int i = 0; i < Count->Molecule; i++) {
    if (System->Molecule[i].InTimestep &&
        System->Molecule[i].Aggregate == -1) {
      int agg = Count->Aggregate;
      System->Molecule[i].Aggregate = agg;
      Aggregate[agg].nMolecules = 1;
      Aggregate[agg].Molecule[0] = i;
      Count->Aggregate++;
    }
  }
} //}}}
// RemovePBCAggregates() //{{{
void RemovePBCAggregates(const double distance, const AGGREGATE *Aggregate,
                         SYSTEM *System, const bool *use_bt) {
  COUNT *Count = &System->Count;
  int **mol_eligible_beads = malloc(Count->MoleculeType * sizeof(int *));
  int *count_eligible_beads = malloc(Count->MoleculeType *
                                     sizeof *count_eligible_beads);
  WrapJoinCoordinates(System, false, true);
  bool eligible = false;
  for (int i = 0; i < Count->MoleculeType; i++) {
    MOLECULETYPE *mt = &System->MoleculeType[i];
    mol_eligible_beads[i] = malloc(mt->nBeads * sizeof(int));
    count_eligible_beads[i] = 0;
    for (int j = 0; j < mt->nBeads; j++) {
      if (use_bt[mt->Bead[j]]) {
        mol_eligible_beads[i][count_eligible_beads[i]] = j;
        count_eligible_beads[i]++;
        eligible = true;
      }
    }
  }
  if (!eligible) {
    err_msg("RemovePBCAggregates(): no bead types for joining");
    PrintError();
    exit(1);
  }

  vec3d *box = &System->Box.Length;
  // helper array indicating whether molecules already moved
  int *list_moved = calloc(Count->Molecule, sizeof *list_moved),
      *list_unmoved = calloc(Count->Molecule, sizeof *list_unmoved);
  // go through aggregates larger than As=1, knitting together //{{{
  for (int i = 0; i < Count->Aggregate; i++) {
    int count_moved = 0,
    count_unmoved = Aggregate[i].nMolecules - 1;
    // set first molecule as already moved
    list_moved[count_moved++] = 0;
    // set all other molecules as unmoved
    for (int j = 0; j < (Aggregate[i].nMolecules - 1); j++) {
      list_unmoved[j] = j + 1;
    }
    while (count_unmoved > 0) {
      // go through all molecule pairs
      for (int jj = 0; jj < count_moved; jj++) {
        int j = list_moved[jj];
        for (int kk = 0; kk < count_unmoved; kk++) {
          int k = list_unmoved[kk];
          bool moved = false;
          // use only moved molecule 'mol1' and unmoved molecule 'mol2'
          int mol1 = Aggregate[i].Molecule[j];
          int mol2 = Aggregate[i].Molecule[k];
          int mtype1 = System->Molecule[mol1].Type;
          int mtype2 = System->Molecule[mol2].Type;

          // go through all bead pairs in the two molecules
          for (int ll = 0; ll < count_eligible_beads[mtype1]; ll++) {
            int l = mol_eligible_beads[mtype1][ll];
            int bead1 = System->Molecule[mol1].Bead[l];
            BEAD *b1 = &System->Bead[bead1];
            for (int mm = 0; mm < count_eligible_beads[mtype2]; mm++) {
              int m = mol_eligible_beads[mtype2][mm];
              int bead2 = System->Molecule[mol2].Bead[m];
              BEAD *b2 = &System->Bead[bead2];
              // calculate distance between 'bead1' and 'bead2'
              vec3d dist = Distance(b1->Position, b2->Position, *box);
              dist.v[0] = VectLength(dist);
              // move 'mol2' (or 'k') if 'bead1' and 'bead2' are in contact
              if (dist.v[0] <= distance) {
                // distance vector between 'bead1' and 'bead2'
                for (int dd = 0; dd < 3; dd++) {
                  dist.v[dd] = b1->Position.v[dd] - b2->Position.v[dd];
                }
                // if 'bead1' and 'bead2' are too far, move 'mol2' //{{{
                for (int dd = 0; dd < 3; dd++) {
                  while (dist.v[dd] > ((*box).v[dd] / 2)) {
                    for (int n = 0; n < System->MoleculeType[mtype2].nBeads; n++) {
                      int id = System->Molecule[mol2].Bead[n];
                      System->Bead[id].Position.v[dd] += (*box).v[dd];
                    }
                    dist.v[dd] = b1->Position.v[dd] - b2->Position.v[dd];
                  }
                  while (dist.v[dd] <= -((*box).v[dd] / 2)) {
                    for (int n = 0; n < System->MoleculeType[mtype2].nBeads; n++) {
                      int id = System->Molecule[mol2].Bead[n];
                      System->Bead[id].Position.v[dd] -= (*box).v[dd];
                    }
                    dist.v[dd] = b1->Position.v[dd] - b2->Position.v[dd];
                  }
                } //}}}
                moved = true;
                for (int x = kk; x < count_unmoved; x++) {
                  list_unmoved[x] = list_unmoved[x+1];
                }
                count_unmoved--;
                list_moved[count_moved++] = k;
                // skip remainder of 'mol2' (or 'k')
                break;
              }
            }
            if (moved) {
              break;
            }
          }
        }
      }
    }
  } //}}}
  free(list_moved);
  free(list_unmoved);
  free(count_eligible_beads);
  for (int i = 0; i < Count->MoleculeType; i++) {
    free(mol_eligible_beads[i]);
  }
  free(mol_eligible_beads);
  // put aggregates' centre of mass into the simulation box //{{{
  for (int i = 0; i < Count->Aggregate; i++) {
    vec3d com = CentreOfMass(Aggregate[i].nBeads, Aggregate[i].Bead, *System);
    // by how many BoxLength's should com by moved?
    // for distant aggregates - it shouldn't happen, but better safe than sorry
    int move[3];
    for (int dd = 0; dd < 3; dd++) {
      move[dd] = com.v[dd] / (*box).v[dd];
      if (com.v[dd] < 0) {
        move[dd]--;
      }
    }
    // move all the beads
    for (int j = 0; j < Aggregate[i].nBeads; j++) {
      int bead = Aggregate[i].Bead[j];
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[bead].Position.v[dd] -= move[dd] * (*box).v[dd];
      }
    }
  } //}}}
} //}}}

// based on options, should an aggregate be used for calculations? //{{{
bool UseAggregate(SYSTEM System, AGGREGATE *Aggregate, int id,
                  AGG_PICKER agg, int *size, double *mass) {
  bool only = true;
  bool x = false;
  *size = 0;
  *mass = 0;
  for (int j = 0; j < Aggregate[id].nMolecules; j++) {
    MOLECULE *mol = &System.Molecule[Aggregate[id].Molecule[j]];
    int mtype = mol->Type;
    MOLECULETYPE *mt = &System.MoleculeType[mtype];
    if (agg.m[mtype]) {
      (*size)++;
      *mass += mt->Mass;
    }
    // if at least one unwanted molecule is present, don't use aggregate
    if (!agg.only[mtype]) {
      only = false;
      break;
    }
    // if at least one molecule isn't exluded, use aggregate
    if (!agg.x[mtype]) {
      x = true;
    }
  }
  if (*size == 0 || // no molecules remained in aggregate
      *size < agg.range[0] || *size > agg.range[1] || // -n: not in range
      !only || // -only: found molecule that weren't supposed to be in
      !x) { // -x: didn't find any un-excluded molecules
    return false;
  } else {
    return true;
  }
} //}}}
// detect -m, -x, -only, and -n options //{{{
void AggPickerOptions(const int argc, char **argv, AGG_PICKER *opt,
                      SYSTEM System) {
  COUNT *Count = &System.Count;
  // '-n' option
  opt->range[0] = 1;
  opt->range[1] = Count->Molecule;
  TwoNumbersOption(argc, argv, "-n", opt->range, 'i');
  if (opt->range[0] > opt->range[1]) {
    SwapInt(&opt->range[0], &opt->range[1]);
  }
  // '-m' option - define aggregate size as sum of those molecule types
  opt->m = calloc(Count->MoleculeType, sizeof *opt->m);
  opt->m_flag = true;
  if (!TypeOption(argc, argv, "-m", 'm', true, opt->m, System)) {
    opt->m_flag = false;
    InitBoolArray(opt->m, Count->MoleculeType, true);
  }
  // '-only' option - use aggregates composed only of specified molecule types
  opt->only = calloc(Count->MoleculeType, sizeof *opt->only);
  opt->only_flag = true;
  if (!TypeOption(argc, argv, "-only", 'm', true, opt->only, System)) {
    opt->only_flag = false;
    InitBoolArray(opt->only, Count->MoleculeType, true);
  }
  // '-x' option - exclude aggregates composed only of specified molecule types
  opt->x = calloc(Count->MoleculeType, sizeof *opt->x);
  opt->x_flag = true;
  if (!TypeOption(argc, argv, "-x", 'm', true, opt->x, System)) {
    opt->x_flag = false;
  }
  // error checking //{{{
  // all molecule specified
  if (opt->x_flag) {
    bool overlap = true; // are all molecule types specified by -x?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (!opt->x[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("with all molecules listed, no aggregates would be detected");
      PrintErrorOption("-x");
      exit(1);
    }
  }
  // molecules specified by -m and -only do not overlap
  if (opt->m_flag && opt->only_flag) {
    bool overlap = false;
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->m[i] && opt->only[i]) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      err_msg("for any aggregate to be used, at least one molecule "
              "must be specified in both options");
      PrintErrorOption("-m/-only");
      exit(1);
    }
  }
  // molecules specified by -only and -x must differ
  if (opt->only_flag && opt->x_flag) {
    bool overlap = true; // do the two array fully overlap?
    for (int i = 0; i < Count->MoleculeType; i++) {
      if (opt->x[i] != opt->only[i]) {
        overlap = false;
        break;
      }
    }
    if (overlap) {
      err_msg("the lists of molecules must be different");
      PrintErrorOption("-x/-only");
      exit(1);
    }
  } //}}}
} //}}}

// create a new aggregate //{{{
static int NewAgg(AGGREGATE *Aggregate, SYSTEM *System,
                  const int i, const int j) {
  int agg_j = System->Count.Aggregate;
  System->Molecule[j].Aggregate = agg_j;
  Aggregate[agg_j].nMolecules = 1;
  Aggregate[agg_j].Molecule[0] = j;
  System->Count.Aggregate++;
  return agg_j;
} //}}}
