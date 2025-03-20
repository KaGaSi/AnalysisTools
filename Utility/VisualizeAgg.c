#include "../AnalysisTools.h"
// TODO: create arrays mapping Bead[id] to BeadCoor[i]=id (and for other xCoor?)
// TODO: implement -x, -m, and -only options
// TODO: something with monomeric beads? //{{{
/*
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
*/ //}}}

// Help() //{{{
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]) {
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
ExtractAgg writes specified aggregates into a coordinate file, placing each \
aggregate into its own timestep. This is, therefore, useful only for \
visualization or further analysis by utilities that do not distinguish \
per-step aggregates (e.g., Surface).\n\n");
  }

  fprintf(ptr, "Usage: %s <in.coor> <in.agg> <output> "
          "<agg size(s)> [options]\n\n", cmd);

  fprintf(ptr, "<in.coor>           input coordinate file\n");
  fprintf(ptr, "<in.agg>            input aggregate file\n");
  fprintf(ptr, "<output>            output coordinate file\n");
  fprintf(ptr, "<agg size(s)>       aggregate size(s) to save\n");
  fprintf(ptr, "[options]\n");
  fprintf(ptr, "  --join            join aggregates (remove pbc)\n");
  // fprintf(ptr, "  -m <name(s)>      agg size means number of <name(s)> "
  //         "molecules in an aggregate\n");
  // fprintf(ptr, "  -x <name(s)>      exclude aggregates containing only "
  //         "specified molecule(s)\n");
  // fprintf(ptr, "  -only <name(s)>   use only aggregates composed of "
  //         "specified molecule type(s)\n");
  fprintf(ptr, "  --range           the first two aggregate sizes specify "
          "a range of aggregate sizes (any following numbers are ignored)\n");
  CommonHelp(error, n, opt);
} //}}}

// structure for options //{{{
struct OPT {
  bool join, range;   // --join --range
  FILE_TYPE fout;     // -o
  COMMON_OPT c;
};
OPT * opt_create(void) {
  return malloc(sizeof(OPT));
} //}}}

int main(int argc, char *argv[]) {

  int common = 8, all = common + 2, count = 0,
      req_arg = 4;
  char option[all][OPT_LENGTH];
  OptionCheck(argc, argv, req_arg, common, all, false, option,
               "-st", "-e", "-sk", "-i", "--verbose", "--silent",
               "--help", "--version", "--join", "--range");

  count = 0; // count mandatory arguments
  OPT *opt = opt_create();
  // arguments & options before reading system data //{{{
  // <in.coor> - input coordinate (and structure) file
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  }
  // <in.agg> - input agg file
  char in_agg[LINE] = "";
  s_strcpy(in_agg, argv[++count], LINE);
  // <output> - output coordinate file
  FILE_TYPE fout;
  s_strcpy(fout.name, argv[++count], LINE);
  fout.type = CoordinateFileType(fout.name);
  if (fout.type == LDATA_FILE) {
    err_msg("lammps data file not allowed as output coordinate file");
    exit(1);
  }

  opt->c = CommonOptions(argc, argv, in);
  // are provided coordinates joined?
  opt->join = BoolOption(argc, argv, "--join");
  opt->range = BoolOption(argc, argv, "--range"); //}}}

  if (!opt->c.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  // <agg sizes> - aggregate sizes to write //{{{
  int *agg_sizes = calloc(Count->Molecule, sizeof *agg_sizes);
  int aggs = 0;
  // use range of aggreate numbers, if --range switch specified
  if (opt->range) {
    long val[2];
    if ((count+2) > argc ||
        !IsNaturalNumber(argv[count+1], &val[0]) ||
        !IsNaturalNumber(argv[count+2], &val[1]) ||
        val[0] == val[1]) {
      err_msg("two different positive numbers needed for size range");
      PrintError();
      Help(argv[0], true, common, option);
      exit(1);
    }
    agg_sizes[0] = val[0];
    agg_sizes[1] = val[1];
    if (agg_sizes[0] > agg_sizes[1]) {
      SwapInt(&agg_sizes[0], &agg_sizes[1]);
    }
    aggs = agg_sizes[1] - agg_sizes[0] + 1;
  } else {
    while (++count < argc && argv[count][0] != '-') {
      // Error - non-numeric argument //{{{
      long val;
      if (!IsNaturalNumber(argv[count], &val)) {
        ErrorNaN("<agg size(s)>");
        Help(argv[0], true, common, option);
        exit(1);
      } //}}}
      agg_sizes[aggs] = atoi(argv[count]);
      aggs++; // number of aggregate sizes
    }
  } //}}}

  if (opt->c.verbose) {
    VerboseOutput(System);
  }

  AGGREGATE *Aggregate;
  InitAggregate(System, &Aggregate);

  double distance = 1;
  // open <in.agg> and skip the first two lines //{{{
  FILE *agg = OpenFile(in_agg, "r");
  while (getc(agg) != '\n')
    ;
  // read Aggregates command if --join is used or...
  if (opt->join) {
    ReadAndSplitLine(agg, SPL_STR, " \t\n");
    // find & flag bead types
    for (count = 5; count < words && split[count][0] != '-'; count++) {
      int btype = FindBeadType(split[count], System);
      if (btype == -1) {
        snprintf(ERROR_MSG, LINE, "bead type %s%s%s from Aggregate command "
                 "does not exist in the system",
                 ErrYellow(), split[count], ErrRed());
        PrintErrorFile(in.stru.name, in_agg, "\0");
        exit(1);
      }
      System.BeadType[btype].Flag = true;
    }
    for (; count < words; count++) {
      if (strcmp(split[count], "-d") == 0) {
        if ((count+1) >= words || !IsRealNumber(split[count+1], &distance)) {
          err_msg("wrong distance in Aggregate command (-d option); "
                  "using 1 for joining the aggregates");
          PrintWarnFile(in_agg, "\0", "\0");
          distance = 1;
        }
      }
    }
  } else { //...skip the command if --join is not used
    while (getc(agg) != '\n')
      ;
  } //}}}

  // array for holding which beads to save
  bool *write = calloc(Count->Bead, sizeof *write);

  // print initial stuff to output coordinate file //{{{
  if (fout.type == VCF_FILE) {
    PrintByline(fout.name, argc, argv);
  } else if (fout.type == VTF_FILE) {
    WriteStructure(fout, System, -1, false, argc, argv);
  } else {
    FILE *out = OpenFile(fout.name, "w");
    fclose(out);
  }
  // write empty lammpstrj timestep containing all beads (vmd needs it)
  if (fout.type == LTRJ_FILE) {
    InitBoolArray(write, Count->Bead, true);
    Count->BeadCoor = Count->Bead;
    for (int i = 0; i < Count->Bead; i++) {
      System.BeadCoor[i] = i;
    }
    WriteTimestep(fout, System, 0, write, argc, argv);
  }
  //}}}

  // open input coordinate file
  FILE *coor = OpenFile(in.coor.name, "r");
  // main loop //{{{
  int count_step = 0;  // count timesteps from the beginning
  int count_saved = 0; // count timesteps (i.e., aggregates) in output file
  int coor_line_count = 0; // count lines in the coor file
  int count_agg_lines = 0; // count lines in the agg file
  while (true) {
    PrintStep(&count_step, opt->c.start, opt->c.silent);
    if (ReadAggregates(agg, in_agg, &System, Aggregate, &count_agg_lines) < 0) {
      count_step--;
      break;
    }
    // decide whether to use this timestep (based on -st/-sk/-e) //{{{
    bool use = false;
    if (UseStep(opt->c, count_step)) {
      use = true;
    } //}}}
    if (use) { //{{{
      if (!ReadTimestep(in, coor, &System, &coor_line_count)) {
        count_step--;
        break;
      }
      if (opt->join) {
        WrapJoinCoordinates(&System, false, true);
        // TODO: distance=1 for now; read agg command (check for -d opt)
        RemovePBCAggregates(distance, Aggregate, &System);
      }
      for (int i = 0; i < Count->Aggregate; i++) {
        // use the aggregate?
        use = false;
        if (opt->range &&
            Aggregate[i].nMolecules >= agg_sizes[0] &&
            Aggregate[i].nMolecules <= agg_sizes[1]) {
          use = true;
        } else {
          for (int j = 0; j < aggs; j++) {
            if (agg_sizes[j] == Aggregate[i].nMolecules) {
              use = true;
              break;
            }
          }
        }
        if (use) {
          InitBoolArray(write, Count->Bead, false);
          for (int j = 0; j < Count->BondedCoor; j++) {
            int id = System.BondedCoor[j];
            for (int k = 0; k < Aggregate[i].nBeads; k++) {
              if (Aggregate[i].Bead[k] == id) {
                if (System.Bead[id].InTimestep) {
                  write[id] = true;
                }
                break;
              }
            }
          }
          WriteTimestep(fout, System, count_step, write, argc, argv);
          count_saved++;
        }
      } //}}}
    } else { //{{{
      if (!SkipTimestep(in, coor, &coor_line_count)) {
        count_step--;
        break;
      }
    } //}}}
    // exit the main loop if reached user-specied end timestep
    if (count_step == opt->c.end) {
      break;
    }
  }
  fclose(coor);
  fclose(agg);
  // print last step count?
  if (!opt->c.silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d ", count_step);
    fprintf(stdout, "(%d aggregates saved)\n", count_saved);
  } //}}}

  // free memory - to make valgrind happy //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  free(write);
  free(agg_sizes);
  free(opt);
  //}}}

  return 0;
}
