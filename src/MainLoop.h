#ifndef MAINLOOP_H
#define MAINLOOP_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

typedef struct {
  int used, // number of used steps (coor or agg file)
      coor, // number of steps read from coor file
      agg, // number of steps read from agg file (should equal coor if used)
      line_count, // line from the file's beginning - for error purposes
      line_count_agg; // ditto for the agg file
} STEP;
static const STEP InitStep = {
  .used = 0,
  .coor = 0,
  .agg = 0,
  .line_count = 0,
  .line_count_agg = 0,
};

typedef void (*callback)(SYSTEM *System, STEP *step, void *userdata);
void MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
                  STEP *step, callback callback_func, void *ud);
/*
 * As MainLoopCoor(), but reads a coordinate and an aggregate timestep in
 * lockstep (so the two files cannot get out of sync). Opens agg_file itself
 * and skips its two header lines. Stops at the end of either file (or at
 * the agg file's Last Step line); exits on malformed aggregate data, as
 * truncating the analysis there would silently skew the results.
 */
void MainLoopCoorAgg(SYSTEM *System, SYS_FILES in, const char *agg_file,
                     COMMON_OPT commons, STEP *step, AGGREGATE *Aggregate,
                     callback callback_func, void *ud);
// As MainLoopCoorAgg(), but for utilities that need no coordinate file
void MainLoopAgg(SYSTEM *System, const char *agg_file, COMMON_OPT commons,
                 STEP *step, AGGREGATE *Aggregate,
                 callback callback_func, void *ud);
// should the given step be used for calculations?
bool UseStep(const COMMON_OPT opt, int step);

#endif
