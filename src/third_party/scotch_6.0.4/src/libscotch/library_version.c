/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_version.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains miscellaneous      **/
/**                routines for handling API- and          **/
/**                version-related routines in the         **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 nov 2010     **/
/**                                 to     16 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* miscellaneous routines.          */
/*                                  */
/************************************/

/*+ This routine returns the version,
*** release and patchlevel numbers of
*** the library being used (useful for
*** dynamic libraries).
*** It returns:
*** - void  : in all cases.
+*/

void
SCOTCH_version (
int * const                 versptr,
int * const                 relaptr,
int * const                 patcptr)
{
  *versptr = SCOTCH_VERSION;
  *relaptr = SCOTCH_RELEASE;
  *patcptr = SCOTCH_PATCHLEVEL;
}
