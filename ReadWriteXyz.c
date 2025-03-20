#include "ReadWriteXyz.h"
#include "General.h"
#include "System.h"

static bool XyzCheckCoorLine(double coor[3]);
static long ReadFirstLine(const char *f, FILE *fr, int *line_count);

SYSTEM XyzReadStruct(const char *file) { //{{{
  SYSTEM Sys;
  InitSystem(&Sys);
  COUNT *Count = &Sys.Count;
  int line_count = 0;
  FILE *fr = OpenFile(file, "r");
  long val = ReadFirstLine(file, fr, &line_count);
  Count->Bead = val;
  Count->BeadCoor = val;
  Count->Unbonded = val;
  Count->UnbondedCoor = val;
  ReallocBead(&Sys);
  // read next line //{{{
  line_count++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    exit(1);
  } //}}}
  // read atoms //{{{
  for (int i = 0; i < Count->Bead; i++) {
    line_count++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      ErrorEOF(file, "incomplete timestep");
      exit(1);
    }
    double coor[3];
    if (!XyzCheckCoorLine(coor)) {
      err_msg("wrong coordinate line");
      PrintErrorFileLine(file, line_count);
      exit(1);
    }
    BEAD *b = &Sys.Bead[i];
    InitBead(b);
    for (int dd = 0; dd < 3; dd++) {
      b->Position[dd] = coor[dd];
    }
    b->InTimestep = true;
    Sys.BeadCoor[i] = i;
    // assign bead type or create a new one
    bool new = true;
    for (int j = 0; j < Count->BeadType; j++) {
      BEADTYPE *bt = &Sys.BeadType[j];
      if (strcmp(split[0], bt->Name) == 0) {
        bt->Number++;
        b->Type = j;
        new = false;
        break;
      }
    }
    if (new) {
      int type = Count->BeadType;
      NewBeadType(&Sys.BeadType, &Count->BeadType, split[0],
                  HIGHNUM, HIGHNUM, HIGHNUM);
      BEADTYPE *bt_new = &Sys.BeadType[type];
      bt_new->Number = 1;
      b->Type = type;
    }
  } //}}}
  fclose(fr);
  FillSystemNonessentials(&Sys, false);
  CheckSystem(Sys, file);
  return Sys;
} //}}}
// XYZReadTimestep() //{{{
int XyzReadTimestep(FILE *fr, const char *file,
                    SYSTEM *System, int *line_count) {
  System->Count.MoleculeCoor = 0;
  long val = ReadFirstLine(file, fr, line_count);
  if (val == -2) {
    return -2; // proper eof - before the first line of the timestep
  }
  if (val > System->Count.Bead) {
    snprintf(ERROR_MSG, LINE,
             "too many beads in the timestep (maximum number is %s%d%s)",
             ErrYellow(), System->Count.Bead, ErrRed());
    PrintErrorFileLine(file, *line_count);
    return -1;
  }
  System->Count.BeadCoor = val;
  // ignore next line //{{{
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    ErrorEOF(file, "missing comment line");
    return -2;
  }
  (*line_count)++; //}}}
  // read coordinates //{{{
  for (int i = 0; i < System->Count.BeadCoor; i++) {
    (*line_count)++;
    if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
      snprintf(ERROR_MSG, LINE, "premature end of file "
               "(%s%d%s coordinate lines instead of %s%d%s)", ErrYellow(),
               i, ErrRed(), ErrYellow(), System->Count.BeadCoor, ErrRed());
      PrintErrorFileLine(file, *line_count);
      return -2;
    }
    double coor[3];
    if (XyzCheckCoorLine(coor)) {
      for (int dd = 0; dd < 3; dd++) {
        System->Bead[i].Position[dd] = coor[dd];
      }
      System->BeadCoor[i] = i;
      double vel[3];
      if (words > 6 &&
          IsRealNumber(split[4], &vel[0]) &&
          IsRealNumber(split[5], &vel[1]) &&
          IsRealNumber(split[6], &vel[2])) {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = vel[dd];
        }
      } else {
        for (int dd = 0; dd < 3; dd++) {
          System->Bead[i].Velocity[dd] = 0;
        }
      }
    } else {
      err_msg("unrecognized line (not enough valid coordinate lines)");
      PrintErrorFileLine(file, *line_count);
      return -1;
    }
  } //}}}
  FillInCoor(System);
  return 1;
} //}}}
bool XyzSkipTimestep(FILE *fr, const char *file, int *line_count) { //{{{
  long val = ReadFirstLine(file, fr, line_count);
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2;
  }
  for (int i = 0; i < val; i++) {
    (*line_count)++;
    int a;
    while ((a = getc(fr)) != '\n' && a != EOF)
      ;
  }
  return true;
} //}}}
static bool XyzCheckCoorLine(double coor[3]) { //{{{
  if (words > 3 &&
      IsRealNumber(split[1], &coor[0]) &&
      IsRealNumber(split[2], &coor[1]) &&
      IsRealNumber(split[3], &coor[2])) {
    return true;
  } else {
    return false;
  }
} //}}}
static long ReadFirstLine(const char *f, FILE *fr, int *line_count) { //{{{
  // read a line
  (*line_count)++;
  if (!ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    return -2; // proper eof - before the first line of a timestep
  }
  // get number of beads
  long val;
  if (words == 0 || !IsNaturalNumber(split[0], &val)) {
    err_msg("wrong first line of an xyz file");
    PrintErrorFileLine(f, *line_count);
    exit(1);
  }
  return val;
} //}}}

void XyzWriteCoor(FILE *fw, const bool *write, const SYSTEM System) { //{{{
  // find out number of beads to save
  int count = 0;
  bool none = true; // to make sure there are beads to save
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    if (write[id]) {
      none = false;
      count++;
    }
  }
  if (none) {
    err_msg("no beads to save");
    PrintWarning();
    return;
  }
  // write pbc on the second line
  fprintf(fw, "%d\n", count);
  WriteBoxLengthAngles(fw, System.Box);
  putc('\n', fw);
  // write the coodinates
  // for (int i = 0; i < System.Count.BeadCoor; i++) {
  //   int id = System.BeadCoor[i];
  for (int id = 0; id < System.Count.Bead; id++) {
    BEAD *bead = &System.Bead[id];
    if (write[id] && bead->InTimestep) {
      int type = bead->Type;
      fprintf(fw, "%8s %8.4f %8.4f %8.4f\n", System.BeadType[type].Name,
              bead->Position[0], bead->Position[1], bead->Position[2]);
    }
  }
} //}}}
