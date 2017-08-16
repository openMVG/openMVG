/* Copyright 2004,2007,2008,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gord.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a sparse matrix graph orderer.  **/
/**                This module contains the data declara-  **/
/**                tions for the main routine.             **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 aug 1996     **/
/**                                 to   : 20 oct 1997     **/
/**                # Version 3.3  : from : 07 oct 1998     **/
/**                                 to   : 31 may 1999     **/
/**                # Version 4.0  : from : 11 dec 2002     **/
/**                                 to   : 27 dec 2004     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   5             /* Number of files in list                */
#define C_FILEARGNBR                3             /* Number of files which can be arguments */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* Source graph input file name */
#define C_filenameordout            fileBlockName (C_fileTab, 1) /* Ordering output file name    */
#define C_filenamelogout            fileBlockName (C_fileTab, 2) /* Log file name                */
#define C_filenamemapout            fileBlockName (C_fileTab, 3) /* Separator mapping file name  */
#define C_filenametreout            fileBlockName (C_fileTab, 4) /* Separator tree file name     */
				                                
#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* Source graph input file */
#define C_filepntrordout            fileBlockFile (C_fileTab, 1) /* Ordering output file    */
#define C_filepntrlogout            fileBlockFile (C_fileTab, 2) /* Log file                */
#define C_filepntrmapout            fileBlockFile (C_fileTab, 3) /* Separator mapping file  */
#define C_filepntrtreout            fileBlockFile (C_fileTab, 4) /* Separator tree file     */

/*+ Process flags. +*/

#define C_FLAGNONE                  0x0000        /* No flags                   */
#define C_FLAGMAPOUT                0x0001        /* Output mapping data        */
#define C_FLAGTREOUT                0x0002        /* Output separator tree data */
#define C_FLAGVERBSTR               0x0004        /* Output strategy string     */
#define C_FLAGVERBTIM               0x0008        /* Output timing information  */
