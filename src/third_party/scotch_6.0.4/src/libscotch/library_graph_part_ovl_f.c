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
/**   NAME       : library_graph_part_ovl_f.c              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                graph partitioning routines with        **/
/**                overlap of the libSCOTCH library.       **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 29 may 2010     **/
/**                                 to     17 oct 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the partitioning routines.     */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHPARTOVL, scotchfgraphpartovl, (     \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    partptr,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          parttab,            \
int * const                 revaptr),           \
(grafptr, partptr, straptr, parttab, revaptr))
{
  *revaptr = SCOTCH_graphPartOvl (grafptr, *partptr, straptr, parttab);
}

/* String lengths are passed at the very
** end of the argument list.
*/

FORTRAN (                                             \
SCOTCHFSTRATGRAPHPARTOVL, scotchfstratgraphpartovl, ( \
SCOTCH_Strat * const        straptr,                  \
const char * const          string,                   \
int * const                 revaptr,                  \
const int                   strnbr),                  \
(straptr, string, revaptr, strnbr))
{
  char * restrict     strtab;                     /* Pointer to null-terminated string */

  if ((strtab = (char *) memAlloc (strnbr + 1)) == NULL) { /* Allocate temporary space */
    errorPrint ("SCOTCHFSTRATGRAPHPARTOVL: out of memory (1)");
    *revaptr = 1;
  }
  memCpy (strtab, string, strnbr);                /* Copy string contents */
  strtab[strnbr] = '\0';                          /* Terminate string     */

  *revaptr = SCOTCH_stratGraphPartOvl (straptr, strtab); /* Call original routine */

  memFree (strtab);
}

/*
**
*/

FORTRAN (                                                       \
SCOTCHFSTRATGRAPHPARTOVLBUILD, scotchfstratgraphpartovlbuild, ( \
SCOTCH_Strat * const        straptr,                            \
const SCOTCH_Num * const    flagptr,                            \
const SCOTCH_Num * const    partptr,                            \
const double * const        balrptr,                            \
int * const                 revaptr),                           \
(straptr, flagptr, partptr, balrptr, revaptr))
{
  *revaptr = SCOTCH_stratGraphPartOvlBuild (straptr, *flagptr, *partptr, *balrptr); /* Call original routine */
}
