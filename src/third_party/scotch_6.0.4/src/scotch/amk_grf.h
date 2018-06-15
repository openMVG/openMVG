/* Copyright 2004,2007,2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : amk_grf.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Decomposition architecture builder.     **/
/**                These lines are the data declarations   **/
/**                for the program routines.               **/
/**                                                        **/
/**   DATES      : # Version 3.0  : from : 06 jul 1995     **/
/**                                 to   : 02 oct 1995     **/
/**                # Version 3.1  : from : 26 mar 1996     **/
/**                                 to   : 26 mar 1996     **/
/**                # Version 3.2  : from : 23 apr 1997     **/
/**                                 to   : 02 jun 1997     **/
/**                # Version 3.3  : from : 15 may 1999     **/
/**                                 to   : 15 may 1999     **/
/**                # Version 5.1  : from : 17 jul 2011     **/
/**                                 to   : 17 jul 2011     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                 3              /* Number of files in list                */
#define C_FILEARGNBR              2              /* Number of files which can be arguments */

#define C_filenamegrfinp          fileBlockName (C_fileTab, 0) /* Source graph input file name  */
#define C_filenametgtout          fileBlockName (C_fileTab, 1) /* Architecture output file name */
#define C_filenamevrtinp          fileBlockName (C_fileTab, 2) /* Vertex list input file name   */

#define C_filepntrgrfinp          fileBlockFile (C_fileTab, 0) /* Source graph input file  */
#define C_filepntrtgtout          fileBlockFile (C_fileTab, 1) /* Architecture output file */
#define C_filepntrvrtinp          fileBlockFile (C_fileTab, 2) /* Vertex list input file   */

/*+ Process flags. +*/

#define C_FLAGVRTINP              0x0001         /* Input vertex list */

#define C_FLAGNONE                0x0000          /* Default flags */

/*
**  The type and structure definitions.
*/

/*+ The sort structure, used to sort graph vertices by label. +*/

typedef struct C_VertSort_ {
  SCOTCH_Num                vlblnum;              /*+ Vertex label  +*/
  SCOTCH_Num                vertnum;              /*+ Vertex number +*/
} C_VertSort;
