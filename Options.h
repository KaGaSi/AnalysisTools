#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"

// print help - function body in each utility
void Help(char cmd[50], bool error, int n, char opt[n][OPT_LENGTH]);
// version/help printing and initial check of provided options
int OptionCheck(int argc, char *argv[], int auto_c, int req, int common,
                int all, char opt[all][OPT_LENGTH], bool check_extra);
int OptionCheck2(int argc, char *argv[], int req, int common, int all,
                 bool check_extra, char opt[all][OPT_LENGTH], ...);
// print version/help and exit
void HelpVersionOption(int argc, char *argv[]);
// print help for common options
void CommonHelp(bool error, int n, char option[n][OPT_LENGTH]);
// detect options common for most utilities
// void CommonOptions(int argc, char *argv[], int length, bool *verbose,
//                    bool *silent, bool *detailed,
//                    int *start, int *end, int *skip);
COMMON_OPT CommonOptions(int argc, char *argv[], int length, SYS_FILES f);

// exclude specified molecule names (-x <mol name(s)>)
bool ExcludeOption(int argc, char *argv[], SYSTEM *System);
// join aggregates, saving the coordinates (-j <filename>)
bool JoinCoorOption(int argc, char *argv[], int *coor_type, char file[]);
// tag which bead types to use (if not present, set to specified value)
bool BeadTypeOption(int argc, char *argv[], char opt[],
                    bool use, bool flag[], SYSTEM System);
// tag which molecule types to use (if not present, set to specified value)
bool MoleculeTypeOption(int argc, char *argv[], char opt[],
                        bool use, bool flag[], SYSTEM System);

// general boolean option
bool BoolOption(int argc, char *argv[], char opt[]);
// general option with multiple integer arguments (up to 'max')
bool IntegerOption(int argc, char *argv[], int max,
                   char opt[], int *count, int *values);
bool IntegerOption1(int argc, char *argv[], char opt[], int *value);
bool IntegerOption2(int argc, char *argv[], char opt[], int value[2]);
// general option with multiple double arguments (up to 'max')
bool DoubleOption(int argc, char *argv[], int max,
                  char opt[], int *count, double values[max]);
bool DoubleOption1(int argc, char *argv[], char opt[], double *value);
bool DoubleOption2(int argc, char *argv[], char opt[], double value[2]);
bool DoubleOption3(int argc, char *argv[], char opt[], double value[3]);
// general option with filename and integer(s) arguments
bool FileIntegerOption(int argc, char *argv[], int min, int max, char opt[],
                       int *values, int *count, char file[]);
bool FileOption(int argc, char *argv[], char opt[], char file[]);
// general option with filename and double(s) arguments
bool FileDoubleOption(int argc, char *argv[], int max, char opt[],
                       double *values, int *count, char file[]);

#if 0 //{{{
// TODO redo
bool MoleculeTypeOption(int argc, char *argv[], char opt[], int *moltype,
                        COUNTS counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeOption2(int argc, char *argv[], char opt[], int *moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeIntOption(int argc, int i, char *argv[], char opt[],
                           int *moltype, int *value, COUNTS Counts,
                           MOLECULETYPE *MoleculeType);
#endif //}}}
#endif
