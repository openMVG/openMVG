/* Copyright 2008,2010,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgmap.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel static mapper.       **/
/**                These lines are the data declaration    **/
/**                for the main routine.                   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 12 jun 2008     **/
/**                                 to   : 18 jul 2011     **/
/**                # Version 6.0  : from : 10 nov 2014     **/
/**                                 to   : 10 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   4             /* Number of files in list */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* Source graph input file name        */
#define C_filenametgtinp            fileBlockName (C_fileTab, 1) /* Target architecture input file name */
#define C_filenamemapout            fileBlockName (C_fileTab, 2) /* Mapping result output file name     */
#define C_filenamelogout            fileBlockName (C_fileTab, 3) /* Log file name                       */

#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* Source graph input file        */
#define C_filepntrtgtinp            fileBlockFile (C_fileTab, 1) /* Target architecture input file */
#define C_filepntrmapout            fileBlockFile (C_fileTab, 2) /* Mapping result output file     */
#define C_filepntrlogout            fileBlockFile (C_fileTab, 3) /* Log file                       */

/*+ Process flags. +*/

#define C_FLAGNONE                  0x0000        /* No flags            */
#define C_FLAGPART                  0x0001        /* Partitioning        */
#define C_FLAGVERBSTR               0x0002        /* Verbose flags       */
#define C_FLAGVERBTIM               0x0004
#define C_FLAGVERBMAP               0x0008
#define C_FLAGVERBMEM               0x0010
#define C_FLAGDEBUG                 0x0020        /* Debugging           */
#define C_FLAGKBALVAL               0x0040        /* Imbalance tolerance */
#define C_FLAGCLUSTER               0x0080        /* Clustering          */

/*
**  The function prototypes.
*/

void                        dgmapStatReduceOp   (double *, double *, int *, MPI_Datatype *);
