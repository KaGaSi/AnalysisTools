/**
 * \file
 * \brief Error prints
 */

#ifndef ERRORS_H
#define ERRORS_H

#define _POSIX_C_SOURCE 200809L

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "General.h"
#include "Structs.h"
#include "AnalysisTools.h"

// print 'warning - <ERROR_MSG>\n' in cyan
void PrintWarning();
// print 'error - <ERROR_MSG>\n' in red
void PrintError();
// print 'error: <option> - <ERROR_MSG>' in red and yellow
void PrintErrorOption(const char *opt);
// print 'warning: <option> - <ERROR_MSG>' in cyan and yellow
void PrintWarnOption(const char *opt);
// print 'error: - <ERROR_MSG>\nFile <file(s)>'
void PrintErrorFile(const char *file1, const char *file2, const char *file3);
// print 'warning: - <ERROR_MSG>\nFile <file(s)>'
void PrintWarnFile(const char *file1, const char *file2, const char *file3);
// print 'error: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintErrorFileLine(const char *file, const int count);
// print 'warning: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintWarnFileLine(const char *file, const int count);
// print 'file <name(s)>' in given colour (Warn: cyan, Error: red)
void WarnPrintFile(const char *file1, const char *file2, const char *file3);
void ErrorPrintFile(const char *file1, const char *file2, const char *file3);
void WarnPrintLine();
void ErrorPrintLine();
void ErrorEOF(const char *file, char *msg);
void ErrorSnprintf();
void ErrorArgNumber(const int count, const int need);
int ErrorExtension(const char *file, const int number,
                   const char extension[][EXTENSION]);
void ErrorOption(const char *option);
void ErrorNaN(const char *option);
void ErrorBeadType(const char *name, const SYSTEM System);
void ErrorMoleculeType(const char *name, const SYSTEM System);
void WarnChargedSystem(const SYSTEM System, const char *file1,
                       const char *file2, const char *file3);
void ErrorStartEnd(const int start, const int end);

void err_msg(const char *str);
#endif
