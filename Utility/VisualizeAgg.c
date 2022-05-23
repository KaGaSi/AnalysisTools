#include "../AnalysisTools.h"
int *InFile;
// TODO: reading agg file

/*
 * TODO: move monomeric beads towards the aggregate
  // put monomeric beads in contact with their aggregates //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    // go through monomeric beads in the aggregate
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int id1 = Aggregate[i].Monomer[j]; // id_move_to = -1;
      // find smallest distance between the monomeric bead and bonded beads //{{{
      double min_dist = 1000000;
      for (int k = 0; k < Aggregate[i].nBeads; k++) {
        int id2 = Aggregate[i].Bead[k];
        if (BeadType[(*Bead)[id2].Type].Use) {
          VECTOR dist;
          dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
          dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
          dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
          double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
          if (d < min_dist) {
            min_dist = d;
          }
          // stop if the monomeric bead is confirmed close to the aggregate
          if (min_dist <= distance) {
            break;
          }
        }
      } //}}}

      // move monomer if it's too far from aggregate //{{{
      while (min_dist > distance) {
        double d, min_dist_2;
        if (id1 == 10) {
          printf("before -x %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.x //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.x -= BoxLength.x;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            } //}}}
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.x += BoxLength.x;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +x %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.x //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.x += BoxLength.x;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.x -= BoxLength.x;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before -y %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.y //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.y -= BoxLength.y;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.y += BoxLength.y;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +y %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.y //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.y += BoxLength.y;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.y -= BoxLength.y;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before -z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by -BoxLength.z //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.z -= BoxLength.z;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.z += BoxLength.z;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("before +z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }
        // test moving by +BoxLength.z //{{{
        if (min_dist > distance) {
          (*Bead)[id1].Position.z += BoxLength.z;
          // find smallest distance between the monomeric bead and bonded beads //{{{
          min_dist_2 = 1000000;
          for (int k = 0; k < Aggregate[i].nBeads; k++) {
            int id2 = Aggregate[i].Bead[k];
            VECTOR dist;
            dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
            dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
            dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
            d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
            if (d < min_dist_2) {
              min_dist_2 = d;
            }
            // stop if the monomeric bead is confirmed close to the aggregate
            if (min_dist_2 <= distance) {
              break;
            }
          } //}}}
          if (min_dist_2 > min_dist) {
            (*Bead)[id1].Position.z -= BoxLength.z;
          } else {
            min_dist = min_dist_2;
          }
        } //}}}
        if (id1 == 10) {
          printf("after +z %lf, (%lf, %lf, %lf)\n", min_dist,
                                                    (*Bead)[10].Position.x,
                                                    (*Bead)[10].Position.y,
                                                    (*Bead)[10].Position.z);
        }

        // find smallest distance between the monomeric bead and bonded beads
        min_dist = 1000000;
        for (int k = 0; k < Aggregate[i].nBeads; k++) {
          int id2 = Aggregate[i].Bead[k];
          VECTOR dist;
          dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
          dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
          dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
          double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
          if (d < min_dist) {
            min_dist = d;
          }
        }
        if (id1 == 10) {
          printf("%3d %4d %lf (%lf, %lf, %lf)\n", i, id1, min_dist,
                                                  (*Bead)[id1].Position.x,
                                                  (*Bead)[id1].Position.y,
                                                  (*Bead)[id1].Position.z);
        }
        if (id1 == 10) {
          printf("%3d %4d %lf (%lf, %lf, %lf)\n", i, id1, min_dist,
                                                  (*Bead)[id1].Position.x,
                                                  (*Bead)[id1].Position.y,
                                                  (*Bead)[id1].Position.z);
        }
      }
    }
  } //}}}

  // test monomeric beads //{{{
  for (int i = 0; i < Counts.Aggregates; i++) {
    for (int j = 0; j < Aggregate[i].nMonomers; j++) {
      int id1 = Aggregate[i].Monomer[j];
      int near;
      double min_dist = 1000000;
      VECTOR near_d;
      for (int k = 0; k < Aggregate[i].nBeads; k++) {
        int id2 = Aggregate[i].Bead[k];
        VECTOR dist;
        dist.x = (*Bead)[id1].Position.x - (*Bead)[id2].Position.x;
        dist.y = (*Bead)[id1].Position.y - (*Bead)[id2].Position.y;
        dist.z = (*Bead)[id1].Position.z - (*Bead)[id2].Position.z;
        double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
        if (d < min_dist) {
          min_dist = d;
          near = id2;
          near_d.x = dist.x;
          near_d.y = dist.y;
          near_d.z = dist.z;
        }
      }
      if (min_dist > distance) {
        printf("%3d %4d %4d %lf (%lf, %lf, %lf); (%lf, %lf, %lf)\n", i, id1, near, min_dist, near_d.x, near_d.y, near_d.z,
                                                                     (*Bead)[id1].Position.x,
                                                                     (*Bead)[id1].Position.y,
                                                                     (*Bead)[id1].Position.z);
      }
    }
  } //}}}
*/

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
VisualizeAgg writes specified aggregates into a new vcf file. This file does \
not have to contain all beads of any type, so it cannot be used for further \
analysis using AnalysisTools utilities.\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <input.agg> <output> <agg size(s)> \
[options]\n\n", cmd);

  fprintf(ptr, "   <input>         input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <input.agg>     input agg file\n");
  fprintf(ptr, "   <output>        output vcf file(s) - one per aggregate size \
(automatic ending <size>.vcf)\n");
  fprintf(ptr, "   <agg size(s)>   aggregate size(s) to save\n");
  fprintf(ptr, "   <options>\n");
  fprintf(ptr, "      --joined       specify that <input> contains joined \
coordinates\n");
  fprintf(ptr, "      -m <name(s)>   agg size means number of <name(s)> \
molecules in an aggregate\n");
  fprintf(ptr, "      -x <name(s)>   exclude aggregates containing only \
specified molecule(s)\n");
  fprintf(ptr, "      -st <int>      starting timestep for calculation\n");
  fprintf(ptr, "      -e <end>       ending timestep for calculation\n");
  CommonHelp(error);
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
  int req_args = 4; //}}}

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
        strcmp(argv[i], "-i") != 0 &&
        strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 &&
        strcmp(argv[i], "-h") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-x") != 0 &&
        strcmp(argv[i], "-st") != 0 &&
        strcmp(argv[i], "-e") != 0) {

      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", input_vsf[LINE] = "";
  snprintf(input_coor, LINE, "%s", argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <input.agg> - input agg file //{{{
  char input_agg[LINE] = "";
  snprintf(input_agg, LINE, "%s", argv[++count]);
  // test if <input.agg> ends with '.agg'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".agg");
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <output> - vcf file(s) for aggregates
  char output[LINE] = "";
  snprintf(output, LINE, "%s", argv[++count]);

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are provided coordinates joined?
  bool joined = BoolOption(argc, argv, "--joined");
  int start, end;
  StartEndTime(argc, argv, &start, &end); //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  bool indexed; // indexed timestep?
  int struct_lines; // number of structure lines (relevant for vtf)
  FullVtfRead(input_vsf, input_coor, false, vtf, &indexed, &struct_lines,
              &Box, &Counts, &BeadType, &Bead, &Index,
              &MoleculeType, &Molecule); //}}}

  // '-m' option //{{{
  int specific_moltype_for_size[Counts.TypesOfMolecules];
  // set all to be used when '-m' option is missing
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    specific_moltype_for_size[i] = 1;
  }
  if (MoleculeTypeOption2(argc, argv, "-m", specific_moltype_for_size,
                          Counts, &MoleculeType)) {
    exit(1);
  } //}}}

  // '-x' option //{{{
  if (ExcludeOption_old(argc, argv, Counts, &MoleculeType)) {
    exit(1);
  }
  // copy Use flag to Write (for '-x' option)
  for (int i = 0; i < Counts.TypesOfMolecules; i++) {
    MoleculeType[i].Write = MoleculeType[i].Use;
  } //}}}

  // <agg sizes> - aggregate sizes to write //{{{
  int *agg_sizes = calloc(Counts.Molecules, sizeof *agg_sizes);
  int aggs = 0;
  while (++count < argc && argv[count][0] != '-') {

    // Error - non-numeric argument //{{{
    if (!IsInteger_old(argv[count]) || atoi(argv[count]) == 0) {
      ErrorNaN("<agg sizes>");
      exit(1);
    } //}}}
    agg_sizes[aggs] = atoi(argv[count]);
    // ensure output string isn't too long for attaching <size>.vcf
    if (agg_sizes[aggs] < 10) {
      output[LINE-1-4] = '\0';
    } else if (agg_sizes[aggs] < 100) {
      output[LINE-2-4] = '\0';
    } else if (agg_sizes[aggs] < 1000) {
      output[LINE-3-4] = '\0';
    } else if (agg_sizes[aggs] < 10000) {
      output[LINE-4-4] = '\0';
    } else {
      output[LINE-100] = '\0';
    }
    aggs++; // number of aggregate sizes
  } //}}}

  double distance; // <distance> parameter from Aggregate command
  int contacts; // <contacts> parameter from Aggregate command - not used here
  ReadAggCommand(BeadType, Counts, input_coor, input_agg, &distance, &contacts);

  // TODO: will change when the agg format changes (at least when Byline is
  //       added to Aggregates*)
  // open input aggregate file and skip the first two lines
  FILE *agg = OpenFile(input_agg, "r");
  char line[LINE];
  // TODO go for while(gets()); approach
  fgets(line, sizeof line, agg);
  fgets(line, sizeof line, agg); //}}}

  // write initial stuff to output density file //{{{
  for (int i = 0; i < aggs; i++) {
    char str[LINE];
    snprintf(str, LINE, "%s%d.vcf", output, agg_sizes[i]);
    FILE *out = OpenFile(str, "w");
    PrintByline(out, argc, argv);
    // print agg size and periodic boundary conditions
    fprintf(out, "# aggregate size: %d\n", agg_sizes[i]);
    fclose(out);
  } //}}}

  // allocate Aggregate struct //{{{
  AGGREGATE *Aggregate = calloc(Counts.Molecules, sizeof (AGGREGATE));
  for (int i = 0; i < Counts.Molecules; i++) {
    // assumes all monomeric beads can be near one aggregate;
    // memory-heavy, but reliable
    Aggregate[i].Monomer = calloc(Counts.Unbonded, sizeof *Aggregate[i].Monomer);
    // assumes all bonded beads can be in one aggregate;
    // memory-heavy, but reliable
    Aggregate[i].Bead = calloc(Counts.Bonded, sizeof *Aggregate[i].Bead);
    // maximum of all molecules can be in one aggregate
    Aggregate[i].Molecule = calloc(Counts.Molecules,
                                   sizeof *Aggregate[i].Molecule);
  } //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput_oldish(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  SkipVtfStructure(vcf, struct_lines);

  count = SkipCoorAggSteps(vcf, input_coor, agg,
                           input_agg, Counts, start, silent);

  // main loop //{{{
  count = 0; // count timesteps in the main loop
  int count_vcf = start - 1; // count timesteps from the beginning
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  while (true) {
    count++;
    count_vcf++;
    // print step? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\rStep: %d", count_vcf);
    } //}}}

    // TODO: will change (probably)
    ReadAggregates(agg, input_agg, &Counts, &Aggregate, BeadType, &Bead,
                   MoleculeType, &Molecule, Index);
    ReadVcfCoordinates(indexed, input_coor, vcf, &Box,
                       Counts, Index, &Bead, &stuff);
    if (!joined) {
      // transform coordinates into fractional ones for non-orthogonal box
      ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      RemovePBCMolecules2(Counts, Box, BeadType, &Bead,
                         MoleculeType, Molecule);
      RemovePBCAggregates(distance, Aggregate, Counts, Box.Length,
                          BeadType, &Bead, MoleculeType, Molecule);
      // transform back to 'normal' coordinates for non-orthogonal box
      FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
    }

    // find correct aggregates and save them //{{{
    for (int i = 0; i < Counts.Aggregates; i++) {
      // test if aggregate 'i' should be used //{{{
      int size = 0;
      // agg size = number of molecules of type 'specific_moltype_for_size'
      for (int j = 0; j < Aggregate[i].nMolecules; j++) {
        int mol_type = Molecule[Aggregate[i].Molecule[j]].Type;
        if (specific_moltype_for_size[mol_type]) {
          size++;
        }
      }
      // is 'size' in provided aggregate sizes?
      int correct_size = -1;
      for (int j = 0; j < aggs; j++) {
        if (agg_sizes[j] == size) {
          correct_size = j;
        }
      } //}}}
      // if '-x' is used, dismiss aggs with only specified molecule(s) //{{{
      bool test = false;
      for (int j = 0; j < size; j++) {
        int moltype = Molecule[Aggregate[i].Molecule[j]].Type;
        if (MoleculeType[moltype].Write) {
          test = true; // a molecule that shouldn't be in agg 'i' is there
          break;
        }
      }
      if (!test) { // should the rest of the for loop agg i be skipped?
        continue;
      } //}}}
      // save the aggregate it should be saved
      if (correct_size != -1) {
        char str[LINE];
        snprintf(str, LINE, "%s%d.vcf", output, agg_sizes[correct_size]);
        FILE *out = OpenFile(str, "r");
        fprintf(out, "%s\n", stuff);
        fprintf(out, "pbc %lf %lf %lf", Box.Length.x,
                                        Box.Length.y,
                                        Box.Length.z);
        fprintf(out, "    %lf %lf %lf\n", Box.alpha, Box.beta, Box.gamma);
        fprintf(out, "indexed\n");
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          fprintf(out, "%6d %lf %lf %lf\n", Bead[id].Indexxx,
                                            Bead[id].Position.x,
                                            Bead[id].Position.y,
                                            Bead[id].Position.z);
        }
        fclose(out);
      }
    } //}}}

    if (LastStep(vcf, NULL) || end == count_vcf) {
      break;
    }
  }
  fclose(vcf);
  fclose(agg);
  // print last step count?
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(stuff); //}}}

  return 0;
}
