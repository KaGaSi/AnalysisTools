#include "Globals.h"

const double PI = 3.14159265359;

char line[LINE];
char *split[SPL_STR];
int words;

char ERROR_MSG[LINE];

const char *BLACK =   "\033[1;30m";
const char *RED =     "\033[1;31m";
const char *GREEN =   "\033[1;32m";
const char *YELLOW =  "\033[1;33m";
const char *BLUE =    "\033[1;34m";
const char *MAGENTA = "\033[1;35m";
const char *CYAN =    "\033[1;36m";
const char *WHITE =   "\033[1;37m";
const char *C_RESET = "\033[0m";

// 'impossible' values
const double CHARGE = 10000.0;
const double MASS = 0.0;
const double RADIUS = 0.0;
const double HIGHNUM = 1000000;
