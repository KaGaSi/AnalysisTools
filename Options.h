/**
 * \file
 * \brief Options usable in utilities
 */

#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "AnalysisTools.h"

#define VERSION "3.3"
#define DATE "August 18, 2020"
#define OPT_LENGTH 16

// print help - function body in each utility
void Help(char cmd[50], bool error);
// print help for common options
void CommonHelp(bool error, int n, char option[n][OPT_LENGTH]);
// detect options common for most utilities
void CommonOptions(int argc, char *argv[], int length, bool *verbose,
                   bool *silent, bool *detailed, bool *vtf_var_coor,
                   int *pbc_xyz, int *ltrj_start_id, int *start,
                   int *end, int *skip);

// print AnalysisTools version number (--version)
bool VersionOption(int argc, char *argv[]);
// exclude specified molecule names (-x <mol name(s)>)
bool ExcludeOption(int argc, char *argv[], SYSTEM *System);
// join aggregates, saving the coordinates (-j <filename>)
bool JoinCoorOption(int argc, char *argv[], int *coor_type, char file[]);
// tag which bead types to use (if not present, set to specified value)
bool BeadTypeOption(int argc, char *argv[], char *opt,
                    bool use, bool flag[], SYSTEM *System);
// tag which molecule types to use (if not present, set to specified value)
bool MoleculeTypeOption(int argc, char *argv[], char *opt,
                        bool use, bool flag[], SYSTEM System);

// general boolean option
bool BoolOption(int argc, char *argv[], char *opt);
// general option with multiple integer arguments (up to 'max')
bool IntegerOption(int argc, char *argv[], int max,
                   char *opt, int *count, int *values);
// general option with multiple double arguments (up to 'max')
bool DoubleOption(int argc, char *argv[], int max,
                  char *opt, int *count, double *values);
// general option with filename and integer(s) arguments
bool FileIntegerOption(int argc, char *argv[], int max, char *opt,
                       int *values, int *count, char *file);

#if 0 //{{{
// TODO redo
bool MoleculeTypeOption(int argc, char *argv[], char *opt, int *moltype,
                        COUNTS counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeOption2(int argc, char *argv[], char *opt, int *moltype,
                         COUNTS Counts, MOLECULETYPE **MoleculeType);
bool MoleculeTypeIntOption(int argc, int i, char *argv[], char *opt,
                           int *moltype, int *value, COUNTS Counts,
                           MOLECULETYPE *MoleculeType);
#endif //}}}
#endif
