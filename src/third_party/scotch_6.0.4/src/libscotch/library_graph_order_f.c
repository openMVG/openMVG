/* Copyright 2004,2007,2008,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_order_f.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                graph ordering routines of the          **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 02 feb 2000     **/
/**                                 to     15 nov 2001     **/
/**                # Version 4.0  : from : 02 feb 2002     **/
/**                                 to     13 dec 2005     **/
/**                # Version 5.0  : from : 04 aug 2007     **/
/**                                 to     31 may 2008     **/
/**                # Version 5.1  : from : 27 mar 2010     **/
/**                                 to     25 jul 2010     **/
/**                # Version 6.0  : from : 08 jan 2012     **/
/**                                 to     08 jan 2012     **/
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
/* for the ordering routines.         */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHORDERINIT, scotchfgraphorderinit, ( \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Ordering * const     ordeptr,            \
SCOTCH_Num * const          permtab,            \
SCOTCH_Num * const          peritab,            \
SCOTCH_Num * const          cblkptr,            \
SCOTCH_Num * const          rangtab,            \
SCOTCH_Num * const          treetab,            \
int * const                 revaptr),           \
(grafptr, ordeptr, permtab, peritab,            \
 cblkptr, rangtab, treetab, revaptr))
{
  *revaptr = SCOTCH_graphOrderInit (grafptr, ordeptr, permtab, peritab, cblkptr, rangtab, treetab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHORDEREXIT, scotchfgraphorderexit, ( \
const SCOTCH_Graph * const  grafptr,            \
SCOTCH_Ordering * const     ordeptr),           \
(grafptr, ordeptr))
{
  SCOTCH_graphOrderExit (grafptr, ordeptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHORDERSAVE, scotchfgraphordersave, ( \
const SCOTCH_Graph * const    grafptr,          \
const SCOTCH_Ordering * const ordeptr,          \
int * const                   fileptr,          \
int * const                   revaptr),         \
(grafptr, ordeptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFGRAPHORDERSAVE: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "w")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFGRAPHORDERSAVE: cannot open output stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_graphOrderSave (grafptr, ordeptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                             \
SCOTCHFGRAPHORDERSAVEMAP, scotchfgraphordersavemap, ( \
const SCOTCH_Graph * const    grafptr,                \
const SCOTCH_Ordering * const ordeptr,                \
int * const                   fileptr,                \
int * const                   revaptr),               \
(grafptr, ordeptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFGRAPHORDERSAVEMAP: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "w")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFGRAPHORDERSAVEMAP: cannot open output stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_graphOrderSaveMap (grafptr, ordeptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                               \
SCOTCHFGRAPHORDERSAVETREE, scotchfgraphordersavetree, ( \
const SCOTCH_Graph * const    grafptr,                  \
const SCOTCH_Ordering * const ordeptr,                  \
int * const                   fileptr,                  \
int * const                   revaptr),                 \
(grafptr, ordeptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFGRAPHORDERSAVETREE: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "w")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFGRAPHORDERSAVETREE: cannot open output stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_graphOrderSaveTree (grafptr, ordeptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                             \
SCOTCHFGRAPHORDERCOMPUTE, scotchfgraphordercompute, ( \
SCOTCH_Graph * const        grafptr,                  \
SCOTCH_Ordering * const     ordeptr,                  \
SCOTCH_Strat * const        stratptr,                 \
int * const                 revaptr),                 \
(grafptr, ordeptr, stratptr, revaptr))
{
  *revaptr = SCOTCH_graphOrderCompute (grafptr, ordeptr, stratptr);
}

/*
**
*/

FORTRAN (                                                     \
SCOTCHFGRAPHORDERCOMPUTELIST, scotchfgraphordercomputelist, ( \
SCOTCH_Graph * const        grafptr,                          \
SCOTCH_Ordering * const     ordeptr,                          \
const SCOTCH_Num *          listptr,                          \
const SCOTCH_Num * const    listtab,                          \
SCOTCH_Strat * const        stratptr,                         \
int * const                 revaptr),                         \
(grafptr, ordeptr, listptr, listtab, stratptr, revaptr))
{
  *revaptr = SCOTCH_graphOrderComputeList (grafptr, ordeptr, *listptr, listtab, stratptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHORDER, scotchfgraphorder, (         \
SCOTCH_Graph * const        grafptr,            \
SCOTCH_Strat * const        stratptr,           \
SCOTCH_Num * const          permtab,            \
SCOTCH_Num * const          peritab,            \
SCOTCH_Num * const          cblkptr,            \
SCOTCH_Num * const          rangtab,            \
SCOTCH_Num * const          treetab,            \
int * const                 revaptr),           \
(grafptr, stratptr, permtab, peritab,           \
 cblkptr, rangtab, treetab, revaptr))
{
  *revaptr = SCOTCH_graphOrder (grafptr, stratptr, permtab, peritab, cblkptr, rangtab, treetab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFGRAPHORDERLIST, scotchfgraphorderlist, ( \
SCOTCH_Graph * const        grafptr,            \
const SCOTCH_Num * const    listptr,            \
const SCOTCH_Num * const    listtab,            \
SCOTCH_Strat * const        stratptr,           \
SCOTCH_Num * const          permtab,            \
SCOTCH_Num * const          peritab,            \
SCOTCH_Num * const          cblkptr,            \
SCOTCH_Num * const          rangtab,            \
SCOTCH_Num * const          treetab,            \
int * const                 revaptr),           \
(grafptr, listptr, listtab, stratptr,           \
 permtab, peritab, cblkptr, rangtab, treetab, revaptr))
{
  *revaptr = SCOTCH_graphOrderList (grafptr, *listptr, listtab, stratptr, permtab, peritab, cblkptr, rangtab, treetab);
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFGRAPHORDERCHECK, scotchfgraphordercheck, ( \
const SCOTCH_Graph * const    grafptr,            \
const SCOTCH_Ordering * const ordeptr,            \
int * const                   revaptr),           \
(grafptr, ordeptr, revaptr))
{
  *revaptr = SCOTCH_graphOrderCheck (grafptr, ordeptr);
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFSTRATGRAPHORDER, scotchfstratgraphorder, ( \
SCOTCH_Strat * const        stratptr,             \
const char * const          string,               \
int * const                 revaptr,              \
const int                   strnbr),              \
(stratptr, string, revaptr, strnbr))
{
  char * restrict     strtab;                     /* Pointer to null-terminated string */

  if ((strtab = (char *) memAlloc (strnbr + 1)) == NULL) { /* Allocate temporary space */
    errorPrint ("SCOTCHFSTRATGRAPHORDER: out of memory (1)");
    *revaptr = 1;
  }
  memCpy (strtab, string, strnbr);                /* Copy string contents */
  strtab[strnbr] = '\0';                          /* Terminate string     */

  *revaptr = SCOTCH_stratGraphOrder (stratptr, strtab); /* Call original routine */

  memFree (strtab);                               /* Prevent compiler warnings */
}

/*
**
*/

FORTRAN (                                                   \
SCOTCHFSTRATGRAPHORDERBUILD, scotchfstratgraphorderbuild, ( \
SCOTCH_Strat * const        stratptr,                       \
const SCOTCH_Num * const    flagval,                        \
const SCOTCH_Num * const    levlnbr,                        \
const double * const        balrat,                         \
int * const                 revaptr),                       \
(stratptr, flagval, levlnbr, balrat, revaptr))
{
  *revaptr = SCOTCH_stratGraphOrderBuild (stratptr, *levlnbr, *flagval, *balrat);
}
