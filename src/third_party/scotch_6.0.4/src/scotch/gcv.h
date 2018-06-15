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
/**   NAME       : gcv.h                                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a graph file converter.         **/
/**                This module contains the data declara-  **/
/**                tions for the main module.              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 apr 1993     **/
/**                                 to     02 apr 1993     **/
/**                # Version 2.0  : from : 28 oct 1994     **/
/**                                 to     16 nov 1994     **/
/**                # Version 3.0  : from : 08 sep 1995     **/
/**                                 to     17 sep 1995     **/
/**                # Version 3.1  : from : 22 mar 1996     **/
/**                                 to     22 mar 1996     **/
/**                # Version 3.2  : from : 04 oct 1996     **/
/**                                 to     04 mar 1997     **/
/**                # Version 3.3  : from : 06 oct 1998     **/
/**                                 to   : 06 oct 1998     **/
/**                # Version 3.4  : from : 13 oct 1999     **/
/**                                 to   : 14 oct 1999     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to   : 29 nov 2003     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines
*/

/*+ File name aliases. +*/

#define C_FILENBR                   3            /* Number of files in list                */
#define C_FILEARGNBR                3            /* Number of files which can be arguments */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* External graph input file name  */
#define C_filenamesrcout            fileBlockName (C_fileTab, 1) /* Source graph output file name   */
#define C_filenamegeoout            fileBlockName (C_fileTab, 2) /* Source graph geometry file name */

#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* External graph input file  */
#define C_filepntrsrcout            fileBlockFile (C_fileTab, 1) /* Source graph output file   */
#define C_filepntrgeoout            fileBlockFile (C_fileTab, 2) /* Source graph geometry file */

/*
**  The type and structure definitions.
*/

/*+ This structure defines the method array element. +*/

typedef struct C_Format_ {
  char                      code;                /* Format type code */
  int                    (* func) ();            /* Function to call */
} C_Format;
