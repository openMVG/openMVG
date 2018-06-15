/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_arch_build_f.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the target architecture building        **/
/**                routine of the libSCOTCH library.       **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 17 mar 2005     **/
/**                                 to     17 mar 2005     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/****************************************/
/*                                      */
/* These routines are the Fortran API   */
/* for the target architecture building */
/* routines.                            */
/*                                      */
/****************************************/

/*
**
*/

FORTRAN (                                           \
SCOTCHFSTRATGRAPHBIPART, scotchfstratgraphbipart, ( \
SCOTCH_Strat * const        stratptr,               \
const char * const          string,                 \
int * const                 revaptr,                \
const int                   strnbr),                \
(stratptr, string, revaptr, strnbr))
{
  char * restrict     strtab;                     /* Pointer to null-terminated string */

  if ((strtab = (char *) memAlloc (strnbr + 1)) == NULL) { /* Allocate temporary space */
    errorPrint ("SCOTCHFSTRATGRAPHBIPART: out of memory (1)");
    *revaptr = 1;
  }
  memCpy (strtab, string, strnbr);                /* Copy string contents */
  strtab[strnbr] = '\0';                          /* Terminate string     */

  *revaptr = SCOTCH_stratGraphBipart (stratptr, strtab); /* Call original routine */

  memFree (strtab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFARCHBUILD, scotchfarchbuild, (           \
SCOTCH_Arch * const         archptr,            \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    listnbr,            \
const SCOTCH_Num * const    listptr,            \
SCOTCH_Strat * const        stratptr,           \
int * const                 revaptr),           \
(archptr, grafptr, listnbr, listptr, stratptr, revaptr))
{
  *revaptr = SCOTCH_archBuild (archptr, grafptr, *listnbr, listptr, stratptr);
}
