#include "MainLoop.h"
#include "Debug.h"
#include "Errors.h"
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

// main loop reading coordinate and aggregate files in lockstep //{{{
void MainLoopCoorAgg(SYSTEM *System, SYS_FILES in, const char *agg_file,
                     COMMON_OPT commons, STEP *step, AGGREGATE *Aggregate,
                     callback callback_func, void *ud) {
  FILE *coor = OpenFile(in.coor.name, "r");
  FILE *agg = OpenFile(agg_file, "r");
  // skip the agg file header (byline + Aggregates command line)
  for (int i = 0; i < 2; i++) {
    step->line_count_agg++;
    int ch;
    while ((ch = getc(agg)) != '\n' && ch != EOF)
      ;
  }
  while (true) {
    PrintStep(&step->coor, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    if (UseStep(commons, step->coor)) {
      if (!ReadTimestep(in, coor, System, &step->line_count)) {
        step->coor--;
        break;
      }
      int agg_ret = ReadAggregates(agg, agg_file, System,
                                   Aggregate, &step->line_count_agg);
      if (agg_ret == -1) { // end of aggregate data - no more full steps
        step->coor--;
        break;
      } else if (agg_ret < 0) { // malformed data - only exit is safe, as
                                // truncating would silently skew results
        err_msg("malformed aggregate data");
        PrintErrorFile(agg_file, "\0", "\0");
        exit(1);
      }
      step->agg++;
      step->used++;
      callback_func(System, step, ud);
    } else {
      if (!SkipTimestep(in, coor, &step->line_count) ||
          !SkipAggregates(agg, agg_file, &step->line_count_agg)) {
        step->coor--;
        break;
      }
      step->agg++;
    }
    // exit the main loop if reached user-specied end timestep
    if (step->coor == commons.end) {
      break;
    }
  }
  fclose(coor);
  fclose(agg);
  PrintLastStep(step->coor, step->used, commons.silent);
}; //}}}

// main loop reading only an aggregate file //{{{
void MainLoopAgg(SYSTEM *System, const char *agg_file, COMMON_OPT commons,
                 STEP *step, AGGREGATE *Aggregate,
                 callback callback_func, void *ud) {
  FILE *agg = OpenFile(agg_file, "r");
  // skip the agg file header (byline + Aggregates command line)
  for (int i = 0; i < 2; i++) {
    step->line_count_agg++;
    int ch;
    while ((ch = getc(agg)) != '\n' && ch != EOF)
      ;
  }
  while (true) {
    PrintStep(&step->agg, commons.start, commons.silent);
    // use every skip-th timestep between start and end
    if (UseStep(commons, step->agg)) {
      int agg_ret = ReadAggregates(agg, agg_file, System,
                                   Aggregate, &step->line_count_agg);
      if (agg_ret == -1) { // end of aggregate data
        step->agg--;
        break;
      } else if (agg_ret < 0) { // malformed data - only exit is safe, as
                                // truncating would silently skew results
        err_msg("malformed aggregate data");
        PrintErrorFile(agg_file, "\0", "\0");
        exit(1);
      }
      step->used++;
      callback_func(System, step, ud);
    } else {
      if (!SkipAggregates(agg, agg_file, &step->line_count_agg)) {
        step->agg--;
        break;
      }
    }
    // exit the main loop if reached user-specied end timestep
    if (step->agg == commons.end) {
      break;
    }
  }
  fclose(agg);
  PrintLastStep(step->agg, step->used, commons.silent);
}; //}}}
