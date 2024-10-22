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

extern char ERROR_MSG[LINE];

// print 'warning - <ERROR_MSG>\n' in cyan
void PrintWarning();
// print 'error - <ERROR_MSG>\n' in red
void PrintError();
// print 'error: <option> - <ERROR_MSG>' in red and yellow
void PrintErrorOption(char *opt);
// print 'warning: <option> - <ERROR_MSG>' in cyan and yellow
void PrintWarnOption(char *opt);
// print 'error: - <ERROR_MSG>\nFile <file(s)>'
void PrintErrorFile(char file1[], char file2[], char file3[]);
// print 'warning: - <ERROR_MSG>\nFile <file(s)>'
void PrintWarnFile(char file1[], char file2[], char file3[]);
// print 'error: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintErrorFileLine(char file[], int count);
// print 'warning: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintWarnFileLine(char file[], int count);
// print 'file <name(s)>' in given colour
void WarnPrintFile(char file1[], char file2[], char file3[]); // in cyan
void ErrorPrintFile(char file1[], char file2[], char file3[]); // in red
void WarnPrintLine(); // in cyan
void ErrorPrintLine();
void ErrorEOF(char file[], char msg[]);
void ErrorSnprintf();
void ErrorArgNumber(int count, int need);
int ErrorExtension(char *file, int number, char extension[][EXTENSION]);
void ErrorOption(char *option);
void ErrorNaN(char *option);
void ErrorBeadType(char name[], SYSTEM System);
void ErrorMoleculeType(char name[], SYSTEM System);
void WarnChargedSystem(SYSTEM System, char file1[], char file2[], char file3[]);
void ErrorStartEnd(int start, int end);

void err_msg(char *str);
#endif

// // FilePrintFile() //{{{
// void FilePrintFile(char *file, char *colour); //}}}
//
// // PrintFileLine() //{{{
// void PrintFileLine(char *file, int line,
//                    char split[SPL_STR][SPL_LEN], int words); //}}}
//
// // ErrorPrintFull() //{{{
// void ErrorPrintFull(char *file, int line,
//                     char split[SPL_STR][SPL_LEN], int words); //}}}
// // ErrorPrintFull2() //{{{
// void ErrorPrintFull2(char *file, int line,
//                     char *split[SPL_STR], int words); //}}}
//
// // WarnStopReading() //{{{
// /*
//  * Warning when stopping file reading due to some error in the file
//  */
// void WarnStopReading(char *vcf_file, int line_count, int step_count,
//                      char split[SPL_STR][SPL_LEN], int words); //}}}
// // WarnStopReading2() //{{{
// /*
//  * Warning when stopping file reading due to some error in the file
//  */
// void WarnStopReading2(char *vcf_file, int line_count, int step_count,
//                      char *split[SPL_STR], int words); //}}}
//
// // ErrorPrintError_old() //{{{
// void ErrorPrintError_old(); //}}}
// // ErrorPrintError() //{{{
// void ErrorPrintError(); //}}}
// // ErrorDiscard() //{{{
// /**
//  * \brief Starting timestep is higher than the number of steps
//  *
//  * \param [in] start  starting timestep
//  * \param [in] step   number of steps read
//  * \param [in] file   coordinate filename
//  * \param [in] coor   pointer to the coordinate file
//  * \return 'true' if the starting step is too high, 'false' otherwise
//  */
// bool ErrorDiscard(int start, int step, char *file, FILE *coor); //}}}
// // WarnPrintWarning() //{{{
// void WarnPrintWarning(); //}}}
// // WarnElNeutrality() //{{{
// /**
//  * \brief Function warning about charged system
//  *
//  * \param [in] Counts    numbers of beads, molecules, etc.
//  * \param [in] BeadType  informationn about bead types
//  * \param [in] file      file name containing the system data
//  */
// void WarnElNeutrality(COUNT Counts, BEADTYPE *BeadType, char *file); //}}}
// // ErrorBeadType_old() //{{{
// /**
//  * Error when non-existent bead is used.
//  *
//  * \param [in] Counts      numbers of beads, molecules, etc.
//  * \param [in] BeadType    information about bead types
//  */
// void ErrorBeadType_old(COUNT Counts, BEADTYPE *BeadType); //}}}
// // ErrorMoleculeType_old() //{{{
// /**
//  * Error when non-existent bead is used.
//  *
//  * \param [in] Counts        numbers of beads, molecules, etc.
//  * \param [in] MoleculeType  information about molecule types
//  */
// void ErrorMoleculeType_old(COUNT Counts, MOLECULETYPE *MoleculeType); //}}}
