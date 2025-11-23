#ifndef OPTIONS_H
#define OPTIONS_H

#define _POSIX_C_SOURCE 200809L

#include "AnalysisTools.h"

// version/help printing and initial check of provided options
int OptionCheck(const int argc, char **argv, const int req, const int common,
                const int all, const bool check_extra,
                char opt[all][OPT_LENGTH], ...);
// print help for common options
void CommonHelp(const bool error, const int n,
                const char option[n][OPT_LENGTH]);
// detect options common for most utilities
COMMON_OPT CommonOptions(const int argc, char **argv, const SYS_FILES f);
// exclude specified molecule names (-x <mol name(s)>)
bool ExcludeOption(const int argc, char **argv, SYSTEM *System);
// tag bead/molecule types to use
bool TypeOption(const int argc, char **argv, const char opt[], const int mode,
                const bool use, bool *flag, const SYSTEM System);
// general boolean option
bool BoolOption(const int argc, char **argv, const char *opt);
// general option with multiple integer/double arguments
bool NumbersOption(const int argc, char **argv, const int max, const char *opt,
                   int *count, void *values, const char type);
bool OneNumberOption(const int argc, char **argv,
                      const char *opt, void *value, const char type);
bool TwoNumbersOption(const int argc, char **argv,
                      const char *opt, void *value, const char type);
bool ThreeNumbersOption(const int argc, char **argv,
                        const char *opt, void *value, const char type);
// general option with filename and integer(s)/double(s) arguments
bool FileNumbersOption(const int argc, char **argv, const int min,
                       const int max, const char *opt, void *values,
                       int *count, char *file, const char type);
// general option with filename argument
bool FileOption(const int argc, char **argv, const char *opt, char *file);
// print help - function body in each utility
void Help(const char cmd[50], const bool error,
          const int n, const char opt[n][OPT_LENGTH]);

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
