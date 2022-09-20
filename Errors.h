/**
 * \file
 * \brief Error prints
 */

#ifndef _ERRORS_H_
#define _ERRORS_H_

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

// simple messages //{{{
// print 'WARNING - <ERROR_MSG>' in cyan
void PrintWarning();
// print 'ERROR - <ERROR_MSG>' in red
void PrintError();
// print 'WARNING: <option> - <ERROR_MSG>' in cyan and yellow
void PrintWarningOption(char *opt);
// print 'ERROR: <option> - <ERROR_MSG>' in red and yellow
void PrintErrorOption(char *opt);
// print 'ERROR: - <ERROR_MSG>\nFile <file(s)>'
void PrintErrorFile(char file1[], char file2[], char file3[]);
// print 'ERROR: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintErrorFileLine(char file1[], int count,
                        char *split[SPL_STR], int words);
// print 'WARNING: - <ERROR_MSG>\nFile <file(s)>, line <count>:\n<line>'
void PrintWarningFileLine(char file1[], char file2[], char file3[], int count,
                          char *split[SPL_STR], int words);
// print 'FILE <name(s)>' in given colour
void PrintFile(FILE *f, char file[], char colour[]);
void WarnPrintFile(char file1[], char file2[], char file3[]); // in cyan
void ErrorPrintFile(char file1[], char file2[], char file3[]); // in red
// print 'Line: <line>|(blank)' in given colours
void PrintLine2(FILE *f, char *split[SPL_STR], int words,
                char *colour1, char *colour2);
void ErrorPrintLine2(char *split[SPL_STR], int words); // in red
void WarnPrintLine(char *split[SPL_STR], int words); // in cyan
 //}}}

void PrintLine(char split[SPL_STR][SPL_LEN], int words,
               int col_line, int col_blank);

// ErrorCoorRead() //{{{
/**
 * \brief Incorrect reading of vcf file
 *
 * \param [in] input_vcf  .vcf coordinate file
 * \param [in] bead       bead's line in its timestep in .vcf file where error occurred
 * \param [in] step       timestep when error occurred
 * \param [in] stuff      comment line of a timestep when error occurred
 */
void ErrorCoorRead(char *input_vcf, int bead, int step, char *stuff); //}}}

// ErrorArgNumber() //{{{
/**
 * \brief Insufficient number of arguments
 *
 * \param [in] count  number of supplied arguments
 * \param [in] need   minimum number of required arguments
 */
void ErrorArgNumber(int count, int need); //}}}

// ErrorExtension() //{{{
/**
 * \brief Wrong file extension
 *
 * \param [in] file       filename
 * \param [in] number     number of correct extension(s)
 * \param [in] extension  correct extension(s)
 * \return <int> of the extension (according to the input array) or -1 if missing a correct one
 */
int ErrorExtension(char *file, int number, char extension[][5]); //}}}

// ErrorNaN() //{{{
/**
 * \brief Non-numeric argument
 *
 * \param [in] option   the option with wrong argument
 */
void ErrorNaN(char *option); //}}}

// ErrorOption() //{{{
/**
 * \brief Unknown option
 *
 * \param [in] option   the unknown option
 */
void ErrorOption(char *option); //}}}

void ErrorBeadType(SYSTEM System);
void ErrorMoleculeType(SYSTEM System);

// ErrorPrintLine() //{{{
/**
 * \brief Print provided strings to error output.
 *
 * \param [in] split     array of strings to prints
 * \param [in] words     number of strings in the split array
 */
void ErrorPrintLine(char split[SPL_STR][SPL_LEN], int words); //}}}

void ErrorEOF(char file1[]);

void WarnChargedSystem(SYSTEM System, char file1[], char file2[], char file3[]);

// ErrorStartEnd() //{{{
/**
 * \brief Error when ending timestep is higher than the starting one
 *
 * \param [in] start   starting timestep
 * \param [in] end     ending timestep
 */
void ErrorStartEnd(int start, int end); //}}}

// FilePrintFile() //{{{
void FilePrintFile(char *file, char *colour); //}}}

// PrintFileLine() //{{{
void PrintFileLine(char *file, int line,
                   char split[SPL_STR][SPL_LEN], int words); //}}}

// ErrorPrintFull() //{{{
void ErrorPrintFull(char *file, int line,
                    char split[SPL_STR][SPL_LEN], int words); //}}}
// ErrorPrintFull2() //{{{
void ErrorPrintFull2(char *file, int line,
                    char *split[SPL_STR], int words); //}}}

// WarnStopReading() //{{{
/*
 * Warning when stopping file reading due to some error in the file
 */
void WarnStopReading(char *vcf_file, int line_count, int step_count,
                     char split[SPL_STR][SPL_LEN], int words); //}}}
// WarnStopReading2() //{{{
/*
 * Warning when stopping file reading due to some error in the file
 */
void WarnStopReading2(char *vcf_file, int line_count, int step_count,
                     char *split[SPL_STR], int words); //}}}

// ErrorPrintError_old() //{{{
void ErrorPrintError_old(); //}}}
// ErrorPrintError() //{{{
void ErrorPrintError(); //}}}
// ErrorDiscard() //{{{
/**
 * \brief Starting timestep is higher than the number of steps
 *
 * \param [in] start  starting timestep
 * \param [in] step   number of steps read
 * \param [in] file   coordinate filename
 * \param [in] coor   pointer to the coordinate file
 * \return 'true' if the starting step is too high, 'false' otherwise
 */
bool ErrorDiscard(int start, int step, char *file, FILE *coor); //}}}
// WarnPrintWarning() //{{{
void WarnPrintWarning(); //}}}
// WarnElNeutrality() //{{{
/**
 * \brief Function warning about charged system
 *
 * \param [in] Counts    numbers of beads, molecules, etc.
 * \param [in] BeadType  informationn about bead types
 * \param [in] file      file name containing the system data
 */
void WarnElNeutrality(COUNTS Counts, BEADTYPE *BeadType, char *file); //}}}
// ErrorBeadType_old() //{{{
/**
 * Error when non-existent bead is used.
 *
 * \param [in] Counts      numbers of beads, molecules, etc.
 * \param [in] BeadType    information about bead types
 */
void ErrorBeadType_old(COUNTS Counts, BEADTYPE *BeadType); //}}}
// ErrorMoleculeType_old() //{{{
/**
 * Error when non-existent bead is used.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] MoleculeType  information about molecule types
 */
void ErrorMoleculeType_old(COUNTS Counts, MOLECULETYPE *MoleculeType); //}}}
#endif
