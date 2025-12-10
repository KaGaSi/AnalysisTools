#ifndef MAINLOOP_H
#define MAINLOOP_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

typedef void (*callback)(SYSTEM *System, void *userdata);
int MainLoopCoor(SYSTEM *System, SYS_FILES in, COMMON_OPT commons,
             callback callback_func, void *ud);

#endif
