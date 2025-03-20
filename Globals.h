#ifndef CONSTANTS_H
#define CONSTANTS_H

extern const double PI;

// maximum splits and string lengths
enum { SPL_STR = 32, SPL_LEN = 64, LINE = 1024 };
// file type signatures
enum { VTF_FILE = 0,
       VSF_FILE = 1,
       VCF_FILE = 2,
       XYZ_FILE = 3,
       LDATA_FILE = 4,
       LTRJ_FILE = 5,
       FIELD_FILE = 6,
       CONFIG_FILE = 7,
       GROM_FILE = 8 };
// maxumum bead/molecule name and file extension lengths
enum { MOL_NAME = 21, BEAD_NAME = 21, EXTENSION = 16 };
// 'impossible' values
extern const double CHARGE;
extern const double MASS;
extern const double RADIUS;
extern const double HIGHNUM; // just some high number

extern char line[LINE];
extern char *split[SPL_STR];
extern int words;

extern char ERROR_MSG[LINE];

extern const char *BLACK;
extern const char *RED;
extern const char *GREEN;
extern const char *YELLOW;
extern const char *BLUE;
extern const char *MAGENTA;
extern const char *CYAN;
extern const char *WHITE;
extern const char *C_RESET;
#endif
