#include "../src/AnalysisTools.h"

// Help message //{{{
const struct HelpHelp HelpDesc = {
  "AggCentres calculates centres of aggregates, saving them to an output "
  "coordinate file, and determines distribution of distances to the nearest N "
  "other aggregate centres for each aggregate centre.",

  "Usage: AggCentres <input> <in.agg> <neighbours> <width> <out.vtf> <out> [options]",
  .args = 6,
  .all = 20,
};
static const struct OptSpec opts[] = {
  COMMON_OPTS[C_I],
  COMMON_OPTS[C_FT],
  COMMON_OPTS[C_ST],
  COMMON_OPTS[C_E],
  COMMON_OPTS[C_SK],
  COMMON_OPTS[C_VERBOSE],
  COMMON_OPTS[C_SILENT],
  COMMON_OPTS[C_HELP],
  COMMON_OPTS[C_VERSION],
  {"<input>", nullptr, "input coordinate file", OPT_ARG},
  {"<in.agg>", nullptr, "input agg file", OPT_ARG},
  {"<neighbours>", nullptr, "number of nearest neighbours", OPT_ARG},
  {"<width>", nullptr, "bin width for the distribution of distances", OPT_ARG},
  {"<out.vtf>", nullptr, "vtf output file with aggregate centres", OPT_ARG},
  {"<out>", nullptr, "output file with distribution of minimal distances",
    OPT_ARG},
  {"--joined", nullptr, "<input> contains joined coordinates", OPT_EXTRA},
  {"-bt", nullptr, "bead types used for calculation (default: all)", OPT_EXTRA},
  {"-x", "<name(s)>", "exclude aggregates containing only specified molecules",
    OPT_EXTRA},
  {"-only", "<name(s)>", "use only aggregates composed of specified molecules",
    OPT_EXTRA},
  {"-n", "<int> <int>", "calculate for aggregate sizes in given range",
    OPT_EXTRA},
  {nullptr},
}; //}}}

struct OPT {
  AGG_PICKER agg;
  bool join,
       *bt;
};

int main(int argc, char *argv[]) {

  // command line arguments before reading the structure //{{{
  OptionCheck(argc, argv, true, HelpDesc, opts);
  OPT opt;
  int count = 0;
  // <input> - input coordinate (and structure) file //{{{
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.coor.name, argv[++count], LINE);
  if (!InputCoorStruct(argc, argv, &in)) {
    exit(1);
  } //}}}
  // <in.agg> - input aggregate file //{{{
  char input_agg[LINE] = "";
  s_strcpy(input_agg, argv[++count], LINE);
  int ext = 1;
  char extension[2][EXTENSION];
  s_strcpy(extension[0], ".agg", EXTENSION);
  if (ErrorExtension(input_agg, ext, extension) == -1) {
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <neighbours> - number of nearest neighbours //{{{
  long tmp;
  if (!IsWholeNumber(argv[++count], &tmp)) {
    ErrorNaN("<neighbours>");
    Help(true, HelpDesc, opts);
    exit(1);
  }
  int neighbours = (int)tmp; //}}}
  // <width> - bin width for distance distribution //{{{
  double width;
  if (!IsPosRealNumber(argv[++count], &width)) {
    ErrorNaN("<width>");
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <out.vtf> - centres of aggregates //{{{
  char output_vtf[LINE];
  s_strcpy(output_vtf, argv[++count], LINE);
  s_strcpy(extension[0], ".vtf", EXTENSION);
  if (ErrorExtension(output_vtf, ext, extension) == -1) {
    Help(true, HelpDesc, opts);
    exit(1);
  } //}}}
  // <out> - distribution of min distances //{{{
  char output[LINE];
  s_strcpy(output, argv[++count], LINE); //}}}

  // options before reading system data
  COMMON_OPT commons = CommonOptions(argc, argv, in);
  // --joined option //{{{
  if (BoolOption(argc, argv, "--joined")) {
    opt.join = false; // joined coordinates supplied, no need to join
  } else {
    opt.join = true;
  } //}}}
  //}}}

  // print command to stdout
  if (!commons.silent) {
    PrintCommand(stdout, argc, argv);
  }

  SYSTEM System = ReadStructure(in, false);
  COUNT *Count = &System.Count;

  AggPickerOptions(argc, argv, &opt.agg, System);

  // -bt option //{{{
  if (!(opt.bt = calloc(Count->BeadType, sizeof *opt.bt))) {
    ErrorAlloc("opt.bt");
  }
  if (!TypeOption(argc, argv, "-bt", 'b', true, opt.bt, System)) {
    InitBoolArray(opt.bt, Count->BeadType, true);
  } //}}}

  // open agg file, read first line, parse second line for distance //{{{
  FILE *agg = OpenFile(input_agg, "r");
  char line_buf[LINE];
  fgets(line_buf, sizeof line_buf, agg);
  if (!ReadAndSplitLine(agg, SPL_STR, " \t\n")) {
    if (snprintf(ERROR_MSG, LINE, "empty %s%s%s line",
                 ErrYellow(), input_agg, ErrRed()) < 0) {
      ErrorSnprintf();
    }
    PrintError();
    exit(1);
  }
  double distance = 1;
  for (int i = 5; i < words; i++) {
    if (strcmp(split[i], "-d") == 0 && (i + 1) < words) {
      if (!IsPosRealNumber(split[i+1], &distance)) {
        distance = 1;
      }
      break;
    }
  } //}}}

  AGGREGATE *Aggregate = nullptr;
  InitAggregate(System, &Aggregate);

  if (commons.verbose) {
    VerboseOutput(System);
  }

  vec3d BoxLength = System.Box.OrthoLength;
  int bins = (int)(Max3(BoxLength.x, BoxLength.y, BoxLength.z) / width);

  // open vtf output file and write header //{{{
  FILE *fw = OpenFile(output_vtf, "w");
  fprintf(fw, "atom 0:%d name C\n", Count->Molecule);
  fprintf(fw, "pbc %lf %lf %lf\n", BoxLength.x, BoxLength.y, BoxLength.z);
  //}}}

  vec3d *centres = calloc(Count->Molecule, sizeof *centres);
  double **distr_min_dist = malloc(bins * sizeof *distr_min_dist);
  for (int i = 0; i < bins; i++) {
    distr_min_dist[i] = calloc(neighbours, sizeof(double));
  }
  int total_count_agg = 0;

  // main loop //{{{
  FILE *coor = OpenFile(in.coor.name, "r");
  int count_step = 0,
      count_used = 0,
      line_count = 0,
      line_count_agg = 2; // first two lines already skipped
  while (true) {
    PrintStep(&count_step, commons.start, commons.silent);
    bool use = UseStep(commons, count_step);
    if (use) { //{{{
      if (!ReadTimestep(in, coor, &System, &line_count) ||
          ReadAggregates(agg, input_agg, &System, Aggregate, &line_count_agg) < 0) {
        count_step--;
        break;
      }
      count_used++;

      // TODO: enable when RemovePBCAggregates() is reliable
      // if (opt.join) {
      //   RemovePBCAggregates(distance, Aggregate, &System, opt.bt);
      // }

      fprintf(fw, "indexed\n");

      // calculate aggregate centres //{{{
      int count_agg = 0;
      for (int i = 0; i < Count->Aggregate; i++) {
        int agg_size;
        double rubbish;
        if (!UseAggregate(System, Aggregate, i, opt.agg, &agg_size, &rubbish)) {
          continue;
        }

        // build list of bead ids matching -bt selection //{{{
        int *list = malloc(Aggregate[i].nBeads * sizeof *list);
        if (!list) {
          ErrorAlloc("list");
        }
        int n = 0;
        for (int j = 0; j < Aggregate[i].nBeads; j++) {
          int id = Aggregate[i].Bead[j];
          if (opt.bt[System.Bead[id].Type]) {
            list[n++] = id;
          }
        } //}}}

        vec3d cog = GeomCentre(n, list, System.Bead);
        free(list);
        cog = RestorePBC(cog, BoxLength);
        centres[count_agg++] = cog;
      } //}}}

      // nearest-neighbour distance distribution //{{{
      for (int i = 0; i < count_agg; i++) {
        double min_dist[neighbours];
        for (int j = 0; j < neighbours; j++) {
          min_dist[j] = 1e7;
        }
        for (int j = 0; j < count_agg; j++) {
          if (i == j) {
            continue;
          }
          double d = VectLength(Distance(centres[i], centres[j], BoxLength));
          for (int k = 0; k < neighbours; k++) {
            if (d < min_dist[k]) {
              for (int l = (neighbours - 1); l > k; l--) {
                min_dist[l] = min_dist[l-1];
              }
              min_dist[k] = d;
              break;
            }
          }
        }
        for (int j = 0; j < neighbours; j++) {
          int bin = (int)(min_dist[j] / width);
          if (bin < bins) {
            distr_min_dist[bin][j]++;
          }
        }
        total_count_agg++;
      } //}}}

      // write centres to vtf output file //{{{
      for (int i = 0; i < count_agg; i++) {
        fprintf(fw, "%5d %lf %lf %lf\n", i, centres[i].x, centres[i].y, centres[i].z);
      }
      for (int j = count_agg; j < Count->Molecule; j++) {
        fprintf(fw, "%5d -1 -1 -1\n", j);
      } //}}}

      ReInitAggregate(System, Aggregate);
    //}}}
    } else {
      if (!SkipTimestep(in, coor, &line_count) ||
          !SkipAggregates(agg, input_agg, &line_count_agg)) {
        count_step--;
        break;
      }
    }

    if (count_step == commons.end) {
      break;
    }
  }
  fclose(fw);
  fclose(coor);
  fclose(agg);
  PrintLastStep(count_step, count_used, commons.silent); //}}}

  // write distance distribution to output file //{{{
  fw = PrintBylineOpenFile(output, argc, argv);
  count = 1;
  fprintf(fw, "# column: (%d) middle of the bin's distance", count++);
  fprintf(fw, ", (%d-%d) nearest", count, count + neighbours - 1);
  putc('\n', fw);
  for (int i = 0; i < bins; i++) {
    fprintf(fw, " %lf", width * (2 * i + 1) / 2);
    for (int j = 0; j < neighbours; j++) {
      fprintf(fw, " %lf", distr_min_dist[i][j] / total_count_agg);
    }
    putc('\n', fw);
  }
  fclose(fw); //}}}

  // free memory //{{{
  FreeAggregate(*Count, Aggregate);
  FreeSystem(&System);
  FreeAggPicker(&opt.agg);
  free(opt.bt);
  free(centres);
  for (int i = 0; i < bins; i++) {
    free(distr_min_dist[i]);
  }
  free(distr_min_dist);
  //}}}

  return 0;
}
