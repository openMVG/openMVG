/* Copyright 2008,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dggath.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the distributed graph file    **/
/**                gathering program.                      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 26 oct 2008     **/
/**                                 to   : 26 oct 2008     **/
/**                # Version 6.0  : from : 10 nov 2014     **/
/**                                 to   : 10 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   2             /* Number of files in list                */
#define C_FILEARGNBR                2             /* Number of files which can be arguments */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* Source graph input file name  */
#define C_filenamesrcout            fileBlockName (C_fileTab, 1) /* Source graph output file name */

#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* Source graph input file  */
#define C_filepntrsrcout            fileBlockFile (C_fileTab, 1) /* Source graph output file */

/*+ Process flags. +*/

#define C_FLAGNONE                  0x0000        /* No flags                        */
#define C_FLAGCHECK                 0x0001        /* Check distributed source graph  */
#define C_FLAGDEBUG                 0x0002        /* Enable easy debugger attachment */
