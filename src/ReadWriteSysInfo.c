#include "ReadWriteSysInfo.h"
#include "General.h"

void WriteSysInfo(const char *filename, const SYSTEM *System) { //{{{
  FILE *fw = OpenFile(filename, "w");
  for (int i = 0; i < System->Count.MoleculeType; i++) {
    fprintf(fw, "%-20s %7d\n", System->MoleculeType[i].Name,
                               System->MoleculeType[i].Number);
  }
  fclose(fw);
} //}}}
// sequentially (re)name molecule types based on -sys file's entries //{{{
void ReadSysInfo(const char *filename, SYSTEM *System) {
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
      if (snprintf(ERROR_MSG, LINE, "%s%s%s: more entries than molecule types "
                   "in the system (only the first %s%d%s entries used)",
                   ErrCyan(), filename, ErrYellow(),
                   ErrCyan(), System->Count.MoleculeType, ErrYellow()) < 0) {
        ErrorSnprintf();
      }
      PrintWarning();
      break;
    }
    s_strcpy(System->MoleculeType[i].Name, split[0], MOL_NAME);
    // an explicit name, same as one from the input file
    System->MoleculeType[i].Named = true;
    i++;
  }
  fclose(fr);
  if (i < System->Count.MoleculeType) {
    if (snprintf(ERROR_MSG, LINE, "%s%s%s: %s%d%s entries for "
                 "a system with %s%d%s molecule types",
                 ErrCyan(), filename, ErrYellow(), ErrCyan(), i, ErrYellow(),
                 ErrCyan(), System->Count.MoleculeType, ErrYellow()) < 0) {
      ErrorSnprintf();
    }
    PrintWarning();
  }
} //}}}
