/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : gmap.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Part of a graph static mapper.          **/
/**                These lines are the data declaration    **/
/**                for the main routine.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 05 jan 1993     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     30 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     08 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     09 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     15 aug 1995     **/
/**                # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     10 nov 1995     **/
/**                # Version 3.2  : from : 04 oct 1996     **/
/**                                 to     18 jul 1997     **/
/**                # Version 3.3  : from : 07 oct 1998     **/
/**                                 to   : 31 may 1999     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to   : 16 jan 2004     **/
/**                # Version 5.0  : from : 12 jun 2008     **/
/**                                 to   : 18 jun 2008     **/
/**                # Version 5.1  : from : 28 aug 2010     **/
/**                                 to   : 18 jul 2011     **/
/**                # Version 6.0  : from : 29 may 2010     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   7             /* Number of files in list */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* Source graph input file name        */
#define C_filenametgtinp            fileBlockName (C_fileTab, 1) /* Target architecture input file name */
#define C_filenamemapout            fileBlockName (C_fileTab, 2) /* Mapping result output file name     */
#define C_filenamelogout            fileBlockName (C_fileTab, 3) /* Log file name                       */
#define C_filenamevfxinp            fileBlockName (C_fileTab, 4) /* Fixed vertex file                   */
#define C_filenamemaoinp            fileBlockName (C_fileTab, 5) /* Old mapping file                    */
#define C_filenamevmlinp            fileBlockName (C_fileTab, 6) /* Vertex migration load file          */

#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* Source graph input file        */
#define C_filepntrtgtinp            fileBlockFile (C_fileTab, 1) /* Target architecture input file */
#define C_filepntrmapout            fileBlockFile (C_fileTab, 2) /* Mapping result output file     */
#define C_filepntrlogout            fileBlockFile (C_fileTab, 3) /* Log file                       */
#define C_filepntrvfxinp            fileBlockFile (C_fileTab, 4) /* Fixed vertex file                   */
#define C_filepntrmaoinp            fileBlockFile (C_fileTab, 5) /* Old mapping file                    */
#define C_filepntrvmlinp            fileBlockFile (C_fileTab, 6) /* Vertex migration load file          */

/*+ Process flags. +*/

#define C_FLAGNONE                  0x0000        /* No flags                   */
#define C_FLAGPART                  0x0001        /* Partitioning               */
#define C_FLAGPARTOVL               0x0002        /* Partitioning with overlap  */
#define C_FLAGVERBSTR               0x0004        /* Verbose flags              */
#define C_FLAGVERBTIM               0x0008
#define C_FLAGVERBMAP               0x0010
#define C_FLAGKBALVAL               0x0020        /* Imbalance tolerance        */
#define C_FLAGCLUSTER               0x0040        /* Clustering                 */
#define C_FLAGFIXED                 0x0080        /* Fixed vertices input file  */
#define C_FLAGRMAPOLD               0x0100        /* Old mapping file           */
#define C_FLAGRMAPRAT               0x0200        /* Edge migration ratio       */
#define C_FLAGRMAPCST               0x0400        /* Vertex migration cost file */

/*
**  The type and structure definitions.
*/

/*+ This structure stores part lists. +*/

typedef struct C_PartList_ {
  SCOTCH_Num                vertnum;              /*+ Number of vertex of which part is neighbor +*/
  SCOTCH_Num                nextidx;              /*+ Pointer to index of next recorded neighbor +*/
} C_PartList;

/*
**  The function prototypes.
*/

void                        C_partSave          (SCOTCH_Graph * restrict const, SCOTCH_Num * restrict const, FILE * const);
void                        C_partViewOvl       (SCOTCH_Graph * restrict const, SCOTCH_Num * restrict const, FILE * const);
