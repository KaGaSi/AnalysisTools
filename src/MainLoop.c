#include "MainLoop.h"
#include "Debug.h"
#include "General.h"
#include "ReadWrite.h"

// TODO: add callback function for step printing?
// TODO: add callback function for while loop breaking?
// TODO: add callback function for use timestep?

// should the given step be used for calculations? //{{{
bool UseStep(const COMMON_OPT opt, const int step) {
  if (step >= opt.start &&
      (step <= opt.end || opt.end == -1) &&
      ((step - opt.start) % opt.skip) == 0) {
    return true;
  } else {
    return false;
  }
} //}}}

// main loop //{{{
void MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
                  STEP *step, callback callback_func, void *ud) {
  FILE *fr = OpenFile(in.coor.name, "r");
  while (true) {
    PrintStep(&step->coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    bool use = false;
    if (UseStep(commons, step->coor)) {
      use = true;
    }
    if (use) {
      if (!ReadTimestep(in, fr, System, &step->line_count)) {
        step->coor--;
        break;
      }
      step->used++;
      callback_func(System, step, ud);
    } else {
      if (!SkipTimestep(in, fr, &step->line_count)) {
        step->coor--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (step->coor == commons.end) {
      break;
    }
  }
  fclose(fr);
  PrintLastStep(step->coor, step->used, commons.silent);
}; //}}}
