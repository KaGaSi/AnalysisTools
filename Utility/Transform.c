#include "../AnalysisTools.h"
int nanosleep(const struct timespec *req, struct timespec *rem);
int *InFile;

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
<mode> == ref 3(or more)x<int> => reflect molecules; \
first 2x<int>: beads (in structurre file ids) specifying normal to plane; \
remaining <int>(s): beads specifying point on the plane (either the one bead \
or arithmetic mean of provided beads' positions)\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <output.vcf> <mode> [options]\n\n", cmd);

  fprintf(ptr, "   <input>        input coordinate file (vcf or vtf format)\n");
  fprintf(ptr, "   <output.vcf>   output coordinate file (vcf format)\n");
  fprintf(ptr, "   <mode>         how to transform molecules (see the help \
text or manual for details)");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      --joined       are provided coordinates joined?\n");
  fprintf(ptr, "      -n <int(s)>    use specified molecules (vsf indices)\n");
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
  int req_args = 3; //}}}

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
        strcmp(argv[i], "--version") != 0 &&
        strcmp(argv[i], "--joined") != 0 &&
        strcmp(argv[i], "-n") != 0) {
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

  // <output.vcf> - output vcf file //{{{
  char output_vcf[LINE] = "";
  snprintf(output_vcf, LINE, "%s", argv[++count]);
  // test if <output.vcf> ends with '.vcf'
  int ext = 1;
  char extension[2][5];
  strcpy(extension[0], ".vcf");
  if (ErrorExtension(output_vcf, ext, extension) == -1) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <mode> - how to transform molecules //{{{
  count++;
  int max_ids = 100,
      transform_ids = 0; // number of values in index array
  long index[max_ids];
  for (int i = 0; i < max_ids; i++) {
    index[i] = -1;
  }
  // <mode> 'ref[lect]' //{{{
  /* the reflecting plane is specified by a normal vector constructed from two
   * bead indices and a point taken as either the third provided bead index or
   * average position from indices 3+ (when more than three are provided)
   * TODO: what about normal perpendicular to a molecule or some such?
   */
  if (strncmp(argv[count],"reflect", 3) == 0) {
    for (int i = 0; ++count < argc && argv[count][0] != '-' && i < 100; i++) {
      // error - not a non-negative number
      if (!IsInteger2(argv[count], &index[i]) || index[i] < 0) {
        strcpy(ERROR_MSG, "<mode> 'reflect': at least three \
non-negative whole numbers are required as arguments");
        ErrorPrintError();
        exit(1);
      }
    }
    // error - less than three arguments
    if (index[2] == -1) {
      strcpy(ERROR_MSG, "<mode> 'reflect': at least three \
non-negative whole numbers are required as arguments");
      ErrorPrintError();
      exit(1);
    }
    // error - the first the numbers are identical
    if (index[0] == index[1]) {
      strcpy(ERROR_MSG, "<mode> 'reflect': the first two numbers \
must be different because they define a normal vector to the reflecting plane");
      ErrorPrintError();
      exit(1);
    }
    // count number of ids
    for (int i = 0; i < max_ids; i++) {
      if (index[i] == -1) {
        break;
      }
      transform_ids++;
    }
  } //}}}
  //}}}

  // options before reading system data //{{{
  bool silent;
  bool verbose;
  CommonOptions(argc, argv, input_vsf, &verbose, &silent, LINE);
  // are joined coordinates provided?
  bool joined = BoolOption(argc, argv, "--joined");
  // use the last step?
  bool last = BoolOption(argc, argv, "--last");
  //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // read information from vtf file(s) //{{{
  BEADTYPE *BeadType; // structure with info about all bead types
  MOLECULETYPE *MoleculeType; // structure with info about all molecule types
  BEAD *Bead; // structure with info about every bead
  int *Index; // link between indices (i.e., Index[Bead[i].Index]=i)
  int *Index_mol; // same as Index, but for molecules
  MOLECULE *Molecule; // structure with info about every molecule
  COUNTS Counts = InitCounts; // structure with number of beads, molecules, etc.
  BOX Box = InitBox; // triclinic box dimensions and angles
  VtfReadStruct(input_vsf, false, &Counts, &BeadType, &Bead, &Index,
                &MoleculeType, &Molecule, &Index_mol);
  InFile = calloc(Counts.BeadsTotal, sizeof *InFile); //}}}

  // '-n' option - specify molecule ids //{{{
  int n_opt_id[100] = {0}, n_opt_number = -1;
  if (MultiIntegerOption(argc, argv, "-n", &n_opt_number, n_opt_id)) {
    exit(1);
  }
  SortArray(n_opt_id, n_opt_number, 0);
  // test the provided indices are good
  for (int i = 0; i < n_opt_number; i++) {
    int id = n_opt_id[i];
    if (id > Counts.HighestResid || Index_mol[id] == -1) {
      strcpy(ERROR_MSG, "option '-n'");
      ErrorPrintError();
      fprintf(stderr, "%sMolecule %s%d%s does not exist%s\n", RED2, YELLOW2, id,
                                                              RED2, C_RESET);
      exit(1);
    }
  } //}}}

  // save all beads
  for (int i = 0; i < Counts.BeadsTotal; i++) {
    Bead[i].Use = true;
  }

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(Counts, BeadType, Bead, MoleculeType, Molecule);
  } //}}}

  // print initial stuff to output vcf file
  FILE *out = OpenFile(output_vcf, "w");
  PrintByline(out, argc, argv);
  fclose(out);

  // main loop //{{{
  int count_vcf = 0, // count steps in the vcf file
      file_line_count = 0; // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  // open input coordinate file
  FILE *vcf = OpenFile(input_coor, "r");
  while (true) {
    count_vcf++;
    // print step info? //{{{
    if (!silent && isatty(STDOUT_FILENO)) {
      if (last) {
        fprintf(stdout, "\rDiscarding step: %d", count_vcf);
      } else {
        fprintf(stdout, "\rStep: %d", count_vcf);
        fflush(stdout);
      }
    } //}}}
    // decide whether this timestep is to be saved
    bool use = true; // TODO maybe later there'll be something
    // read and write the timestep, if it should be saved //{{{
    if (use) {
      if (!VtfReadTimestep(vcf, input_coor, &Box, &Counts, BeadType, &Bead,
                           Index, MoleculeType, Molecule,
                           &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
      // remove pbc if provided molecules aren't joined
      if (!joined) {
        ToFractionalCoor(Counts.BeadsCoor, &Bead, Box);
        RemovePBCMolecules_new(Counts, Box, BeadType, &Bead,
                           MoleculeType, Molecule);
        FromFractionalCoor(Counts.BeadsCoor, &Bead, Box);
      }

      // mirror the molecule
      for (int i = 0; i < Counts.Molecules; i++) {
        use = false;
        for (int j = 0; j < n_opt_number; j++) {
          if (Molecule[i].Index == n_opt_id[j]) {
            use = true;
            break;
          }
        }
        if (use) {
          int mtype = Molecule[i].Type;
          // normal to mirroring plane
//        int id1 = Molecule[i].Bead[0],
//            id2 = Molecule[i].Bead[1];
          // TODO error when index[] is higher than i's number of beads
          int id1 = Molecule[i].Bead[index[0]],
              id2 = Molecule[i].Bead[index[1]];
          VECTOR normal;
          normal.x = Bead[id1].Position.x - Bead[id2].Position.x;
          normal.y = Bead[id1].Position.y - Bead[id2].Position.y;
          normal.z = Bead[id1].Position.z - Bead[id2].Position.z;
          // point on the plane
          VECTOR point = {0, 0, 0};
          for (int j = 2; j < transform_ids; j++) {
            int id = Molecule[i].Bead[index[j]];
            point.x += Bead[id].Position.x;
            point.y += Bead[id].Position.y;
            point.z += Bead[id].Position.z;
          }
          point.x /= transform_ids - 2;
          point.y /= transform_ids - 2;
          point.z /= transform_ids - 2;
          // d from ax+by+cz+d=0 - equation of reflecting plane
          double d_plane = -(normal.x * point.x +
                             normal.y * point.y +
                             normal.z * point.z);
  //      printf("plane: %lfx + %lfy + %lfz + %lf = 0\n", normal.x, normal.y,
  //                                                      normal.z, d_plane);
  //      printf("point on plane: (%lf, %lf, %lf)\n", point.x, point.y, point.z);
          for (int j = 0; j < MoleculeType[mtype].nBeads; j++) {
            int id = Molecule[i].Bead[j];
            // intersection of the line containing bead id with the plane
            double t = -(d_plane + normal.x * Bead[id].Position.x +
                                   normal.y * Bead[id].Position.y +
                                   normal.z * Bead[id].Position.z);
            t = t / (SQR(normal.x) + SQR(normal.y) + SQR(normal.z));
            VECTOR intersect;
            intersect.x = Bead[id].Position.x + normal.x * t;
            intersect.y = Bead[id].Position.y + normal.y * t;
            intersect.z = Bead[id].Position.z + normal.z * t;
  //        printf("line: x = %lf + %lft\n", Bead[id].Position.x, normal.x);
  //        printf("      y = %lf + %lft\n", Bead[id].Position.y, normal.y);
  //        printf("      z = %lf + %lft\n", Bead[id].Position.z, normal.z);
  //        printf("intersection: (%lf, %lf, %lf)\n", intersect.x,
  //                                                  intersect.y,
  //                                                  intersect.z);
            Bead[id].Position.x = 2 * intersect.x - Bead[id].Position.x;
            Bead[id].Position.y = 2 * intersect.y - Bead[id].Position.y;
            Bead[id].Position.z = 2 * intersect.z - Bead[id].Position.z;
          }
        }
      }

      // write to output .vcf file
      out = OpenFile(output_vcf, "a");
      VtfWriteCoorIndexed(out, stuff, Counts, Bead, Box);
      fclose(out);
      //}}}
    // skip the timestep, if it shouldn't be saved //{{{
    } else {
      if (!VtfSkipTimestep(vcf, input_coor, &file_line_count, count_vcf)) {
        count_vcf--;
        break;
      }
    } //}}}
    // definitely break the loop - later, there'll be option for specific step
    break;
  }
  fclose(vcf); //}}}

  // print last step count? //{{{
  if (!silent) {
    if (isatty(STDOUT_FILENO)) {
      fflush(stdout);
      fprintf(stdout, "\r                          \r");
    }
    fprintf(stdout, "Last Step: %d\n", count_vcf);
  } //}}}

  // free memory
  FreeSystemInfo(Counts, &MoleculeType, &Molecule, &BeadType, &Bead, &Index);
  free(Index_mol);
  free(stuff);
  free(InFile);

  return 0;
}
