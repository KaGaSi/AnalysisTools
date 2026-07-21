#ifndef MAINLOOP_H
#define MAINLOOP_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

typedef struct {
  int used, // number of used steps (coor or agg file)
      coor, // number of steps read from coor file
      agg, // number of steps read from agg file (should equal coor if used)
      line_count, // line from the file's beginning - for error purposes
      line_count_agg; // same for the agg file
} STEP;
static const STEP InitStep = {
  .used = 0,
  .coor = 0,
  .agg = 0,
  .line_count = 0,
  .line_count_agg = 0,
};

typedef void (*callback)(SYSTEM *System, STEP *step, void *userdata);
// main loop for utilities reading coordinate files only
void MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
                  STEP *step, callback callback_func, void *ud);
// main loop for utilities reading both coordinate and aggregate files
void MainLoopCoorAgg(SYSTEM *System, SYS_FILES in, const char *agg_file,
                     COMMON_OPT commons, STEP *step, AGGREGATE *Aggregate,
                     callback callback_func, void *ud);
// main loop for utilities reading aggregate files only
void MainLoopAgg(SYSTEM *System, const char *agg_file, COMMON_OPT commons,
                 STEP *step, AGGREGATE *Aggregate,
                 callback callback_func, void *ud);
// should the given step be used for calculations?
bool UseStep(const COMMON_OPT opt, int step);

#endif
