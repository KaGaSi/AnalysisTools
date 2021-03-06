/**
 * \file
 * \brief Error prints
 */

#ifndef _ERRORS_H_
#define _ERRORS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "General.h"
#include "Structs.h"

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

// ErrorExtension() //{{{
/**
 * \brief Wrong file extension
 *
 * \param [in] file       filename
 * \param [in] number     number of correct extension(s)
 * \param [in] extension  correct extension(s)
 * \return 'true' if wrong extension, 'false' otherwise
 */
bool ErrorExtension(char *file, int number, char extension[][5]); //}}}

// ErrorFileOpen() //{{{
/**
 * \brief Cannot open file
 *
 * \param [in] file  filename
 * \param [in] mode  open mode - r(ead), w(rite), a(ppend)
 */
void ErrorFileOpen(char *file, char mode); //}}}

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

// ErrorBeadType() //{{{
/**
 * Error when non-existent bead is used.
 *
 * \param [in] Counts      numbers of beads, molecules, etc.
 * \param [in] BeadType    information about bead types
 */
void ErrorBeadType(COUNTS Counts, BEADTYPE *BeadType); //}}}

// ErrorMoleculeType() //{{{
/**
 * Error when non-existent bead is used.
 *
 * \param [in] Counts        numbers of beads, molecules, etc.
 * \param [in] MoleculeType  information about molecule types
 */
void ErrorMoleculeType(COUNTS Counts, MOLECULETYPE *MoleculeType); //}}}

// ErrorPrintLine() //{{{
/**
 * \brief Print provided strings to error output.
 *
 * \param [in] split     array of strings to prints
 * \param [in] words     number of strings in the split array
 */
void ErrorPrintLine(char split[30][100], int words); //}}}

// WarnElNeutrality() //{{{
/**
 * \brief Function warning about charged system
 *
 * \param [in] Counts    numbers of beads, molecules, etc.
 * \param [in] BeadType  informationn about bead types
 * \param [in] file      file name containing the system data
 */
void WarnElNeutrality(COUNTS Counts, BEADTYPE *BeadType, char *file); //}}}
#endif
