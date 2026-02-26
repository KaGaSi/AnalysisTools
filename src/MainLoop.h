#ifndef MAINLOOP_H
#define MAINLOOP_H

#define _POSIX_C_SOURCE 200809L

#include "Structs.h"

typedef struct {
  int used, // number of used steps (coor or agg file)
      coor, // number of steps read from coor file
      agg, // number of steps read from agg file (should equal coor if used)
      line_count; // line from the file's beginning - for error purposes
} STEP;
static const STEP InitStep = {
  .used = 0,
  .coor = 0,
  .agg = 0,
  .line_count = 0,
};

typedef void (*callback)(SYSTEM *System, STEP *step, void *userdata);
void MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
                  STEP *step, callback callback_func, void *ud);

#endif
