#include "AnalysisTools.h"
#include "MainLoop.h"

// TODO: add callback function for step printing?
// TODO: add callback function for while loop breaking?
// TODO: add callback function for use timestep?

// main loop //{{{
int MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
             callback callback_func, void *ud) {
  FILE *fr = OpenFile(in.coor.name, "r");
  int count_coor = 0, // count steps in the vcf file
      count_used = 0, // count steps in output file
      line_count = 0; // count lines in the vcf file
  while (true) {
    PrintStep(&count_coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, count_coor)) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(in, fr, System, &line_count)) {
        count_coor--;
        break;
      }
      count_used++;
      callback_func(System, ud);
    } else {
      if (!SkipTimestep(in, fr, &line_count)) {
        count_coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (count_coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(count_coor, count_used, commons.silent);
  return count_used;
}; //}}}
