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
/**   NAME       : amk_hy.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the target architecture file    **/
/**                for hypercube graphs.                   **/
/**                Here are the data declaration for the   **/
/**                target machine architecture functions.  **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 14 nov 1994     **/
/**                                 to   : 14 nov 1994     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to   : 07 jul 1995     **/
/**                # Version 3.1  : from : 30 may 1996     **/
/**                                 to   : 30 may 1996     **/
/**                # Version 3.2  : from : 31 may 1997     **/
/**                                 to   : 31 may 1997     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 02 oct 1998     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/** File name aliases. **/

#define C_FILENBR                   1             /* Number of files in list                */
#define C_FILEARGNBR                1             /* Number of files which can be arguments */

#define C_filenametgtout            fileBlockName (C_fileTab, 0) /* Architecture output file name */

#define C_filepntrtgtout            fileBlockFile (C_fileTab, 0) /* Architecture output file */

/*
**  The macro definitions.
*/

#ifndef abs
#define abs(a)                      (((a) >= 0) ? (a) : -(a))
#endif /* abs */
