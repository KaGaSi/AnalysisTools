#include "ReadWriteSysInfo.h"

void WriteSysInfo(const char *filename, const SYSTEM *System) { //{{{
  FILE *fw = OpenFile(filename, "w");
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    fprintf(fw, "%-20s %7d\n", System->MoleculeType[i].Name,
                               System->MoleculeType[i].Number);
  }
  fclose(fw);
} //}}}
void ReadSysInfo(const char *filename, SYSTEM *System) { //{{{
  FILE *fr = OpenFile(filename, "r");
  int i = 0;
  while (ReadAndSplitLine(fr, SPL_STR, " \t\n")) {
    if (words == 0 || split[0][0] == '#' || words < 2) {
      continue;
    }
    long n;
    // skip lines without a number
    if (!IsWholeNumber(split[1], &n)) {
      continue;
    }
    if (i >= System->Count.MoleculeType) {
      if (snprintf(ERROR_MSG, LINE,
                   "%s has more entries than molecule types in system",
                   filename) < 0) {
        ErrorSnprintf();
      }
      PrintWarning();
      break;
    }
    s_strcpy(System->MoleculeType[i].Name, split[0], MOL_NAME);
    i++;
  }
  fclose(fr);
  if (i < System->Count.MoleculeType) {
    if (snprintf(ERROR_MSG, LINE,
                 "%s has %d entries but system has %d molecule types",
                 filename, i, System->Count.MoleculeType) < 0) {
      ErrorSnprintf();
    }
    PrintWarning();
  }
} //}}}
