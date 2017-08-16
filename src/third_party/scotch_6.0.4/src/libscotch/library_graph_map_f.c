/* Copyright 2004,2007,2010-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_map_f.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                mapping routines of the libSCOTCH       **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 02 dec 1999     **/
/**                                 to     15 nov 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     12 dec 2005     **/
/**                # Version 5.1  : from : 27 mar 2010     **/
/**                                 to     31 aug 2011     **/
/**                # Version 6.0  : from : 17 apr 2011     **/
/**                                 to     23 nov 2012     **/
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
/* for the mapping routines.          */
/*                                    */
/**************************************/

/*
**
*/
/* TODO update fortran interface... */
FORTRAN (                                       \
SCOTCHFGRAPHMAPINIT, scotchfgraphmapinit, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mappptr,            \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Num * const          mapptab,            \
int * const                 revaptr),           \
(grafptr, mappptr, archptr, mapptab, revaptr))
{
  *revaptr = SCOTCH_graphMapInit (grafptr, mappptr, archptr, mapptab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPEXIT, scotchfgraphmapexit, (     \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Mapping * const      mappptr),           \
(grafptr, mappptr))
{
  SCOTCH_graphMapExit (grafptr, mappptr);
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFGRAPHMAPCOMPUTE, scotchfgraphmapcompute, ( \
SCOTCH_Graph * const        grafptr,              \
SCOTCH_Mapping * const      mappptr,              \
SCOTCH_Strat * const        straptr,              \
int * const                 revaptr),             \
(grafptr, mappptr, straptr, revaptr))
{
  *revaptr = SCOTCH_graphMapCompute (grafptr, mappptr, straptr);
}

/*
**
*/

FORTRAN (                                                   \
SCOTCHFGRAPHMAPFIXEDCOMPUTE, scotchfgraphmapfixedcompute, ( \
SCOTCH_Graph * const        grafptr,                        \
SCOTCH_Mapping * const      mappptr,                        \
SCOTCH_Strat * const        straptr,                        \
int * const                 revaptr),                       \
(grafptr, mappptr, straptr, revaptr))
{
  *revaptr = SCOTCH_graphMapFixedCompute (grafptr, mappptr, straptr);
}

/*
**
*/

FORTRAN (                                             \
SCOTCHFGRAPHREMAPCOMPUTE, scotchfgraphremapcompute, ( \
SCOTCH_Graph * const        grafptr,                  \
SCOTCH_Mapping * const      mappptr,                  \
SCOTCH_Mapping * const      mapoptr,                  \
const double * const        emraptr,                  \
const SCOTCH_Num * const    vmlotab,                  \
SCOTCH_Strat * const        straptr,                  \
int * const                 revaptr),                 \
(grafptr, mappptr, mapoptr, emraptr, vmlotab, straptr, revaptr))
{
  *revaptr = SCOTCH_graphRemapCompute (grafptr, mappptr, mapoptr, *emraptr, vmlotab, straptr);
}

/*
**
*/

FORTRAN (                                                       \
SCOTCHFGRAPHREMAPFIXEDCOMPUTE, scotchfgraphremapfixedcompute, ( \
SCOTCH_Graph * const        grafptr,                            \
SCOTCH_Mapping * const      mappptr,                            \
SCOTCH_Mapping * const      mapoptr,                            \
const double * const        emraptr,                            \
const SCOTCH_Num * const    vmlotab,                            \
SCOTCH_Strat * const        straptr,                            \
int * const                 revaptr),                           \
(grafptr, mappptr, mapoptr, emraptr, vmlotab, straptr, revaptr))
{
  *revaptr = SCOTCH_graphRemapFixedCompute (grafptr, mappptr, mapoptr, *emraptr, vmlotab, straptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAP, scotchfgraphmap, (             \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          parttab,            \
int * const                 revaptr),           \
(grafptr, archptr, straptr, parttab, revaptr))
{
  *revaptr = SCOTCH_graphMap (grafptr, archptr, straptr, parttab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHMAPFIXED, scotchfgraphmapfixed, (   \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          parttab,            \
int * const                 revaptr),           \
(grafptr, archptr, straptr, parttab, revaptr))
{
  *revaptr = SCOTCH_graphMapFixed (grafptr, archptr, straptr, parttab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHREMAP, scotchfgraphremap, (         \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Arch * const   archptr,            \
SCOTCH_Num * const          parotab,            \
const double * const        emraptr,            \
const SCOTCH_Num * const    vmlotab,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          parttab,            \
int * const                 revaptr),           \
(grafptr, archptr, parotab, emraptr, vmlotab, straptr, parttab, revaptr))
{
  *revaptr = SCOTCH_graphRemap (grafptr, archptr, parotab, *emraptr, vmlotab, straptr, parttab);
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFGRAPHREMAPFIXED, scotchfgraphremapfixed, ( \
SCOTCH_Graph * const        grafptr,              \
const SCOTCH_Arch * const   archptr,              \
SCOTCH_Num * const          parotab,              \
const double * const        emraptr,              \
const SCOTCH_Num * const    vmlotab,              \
SCOTCH_Strat * const        straptr,              \
SCOTCH_Num * const          parttab,              \
int * const                 revaptr),             \
(grafptr, archptr, parotab, emraptr, vmlotab, straptr, parttab, revaptr))
{
  *revaptr = SCOTCH_graphRemapFixed (grafptr, archptr, parotab, *emraptr, vmlotab, straptr, parttab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHPART, scotchfgraphpart, (           \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    partptr,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          mapptab,            \
int * const                 revaptr),           \
(grafptr, partptr, straptr, mapptab, revaptr))
{
  *revaptr = SCOTCH_graphPart (grafptr, *partptr, straptr, mapptab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHPARTFIXED, scotchfgraphpartfixed, ( \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    partptr,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          mapptab,            \
int * const                 revaptr),           \
(grafptr, partptr, straptr, mapptab, revaptr))
{
  *revaptr = SCOTCH_graphPartFixed (grafptr, *partptr, straptr, mapptab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHREPART, scotchfgraphrepart, (       \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    partptr,            \
SCOTCH_Num * const          parotab,            \
const double * const        emraptr,            \
const SCOTCH_Num * const    vmlotab,            \
SCOTCH_Strat * const        straptr,            \
SCOTCH_Num * const          mapptab,            \
int * const                 revaptr),           \
(grafptr, partptr, parotab, emraptr, vmlotab, straptr, mapptab, revaptr))
{
  *revaptr = SCOTCH_graphRepart (grafptr, *partptr, parotab, *emraptr, vmlotab, straptr, mapptab);
}

/*
**
*/

FORTRAN (                                           \
SCOTCHFGRAPHREPARTFIXED, scotchfgraphrepartfixed, ( \
SCOTCH_Graph * const        grafptr,                \
const SCOTCH_Num * const    partptr,                \
SCOTCH_Num * const          parotab,                \
const double * const        emraptr,                \
const SCOTCH_Num * const    vmlotab,                \
SCOTCH_Strat * const        straptr,                \
SCOTCH_Num * const          mapptab,                \
int * const                 revaptr),               \
(grafptr, partptr, parotab, emraptr, vmlotab, straptr, mapptab, revaptr))
{
  *revaptr = SCOTCH_graphRepartFixed (grafptr, *partptr, parotab, *emraptr, vmlotab, straptr, mapptab);
}

/* String lengths are passed at the very
** end of the argument list.
*/

FORTRAN (                                       \
SCOTCHFSTRATGRAPHMAP, scotchfstratgraphmap, (   \
SCOTCH_Strat * const        stratptr,           \
const char * const          string,             \
int * const                 revaptr,            \
const int                   strnbr),            \
(stratptr, string, revaptr, strnbr))
{
  char * restrict     strtab;                     /* Pointer to null-terminated string */

  if ((strtab = (char *) memAlloc (strnbr + 1)) == NULL) { /* Allocate temporary space */
    errorPrint ("SCOTCHFSTRATGRAPHMAP: out of memory (1)");
    *revaptr = 1;
  }
  memCpy (strtab, string, strnbr);                /* Copy string contents */
  strtab[strnbr] = '\0';                          /* Terminate string     */

  *revaptr = SCOTCH_stratGraphMap (stratptr, strtab); /* Call original routine */

  memFree (strtab);
}

/*
**
*/

FORTRAN (                                               \
SCOTCHFSTRATGRAPHMAPBUILD, scotchfstratgraphmapbuild, ( \
SCOTCH_Strat * const        stratptr,                   \
const SCOTCH_Num * const    flagval,                    \
const SCOTCH_Num * const    partnbr,                    \
const double * const        kbalptr,                    \
int * const                 revaptr),                   \
(stratptr, flagval, partnbr, kbalptr, revaptr))
{
  *revaptr = SCOTCH_stratGraphMapBuild (stratptr, *flagval, *partnbr, *kbalptr);
}

/*
**
*/

FORTRAN (                                                       \
SCOTCHFSTRATGRAPHCLUSTERBUILD, scotchfstratgraphclusterbuild, ( \
SCOTCH_Strat * const        stratptr,                           \
const SCOTCH_Num * const    flagval,                            \
const SCOTCH_Num * const    pwgtval,                            \
const double * const        densval,                            \
const double * const        bbalptr,                            \
int * const                 revaptr),                           \
(stratptr, flagval, pwgtval, densval, bbalptr, revaptr))
{
  *revaptr = SCOTCH_stratGraphClusterBuild (stratptr, *flagval, *pwgtval, *densval, *bbalptr);
}
