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
/**   NAME       : bgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the bipartition    **/
/**                graph data structure handling           **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     09 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     01 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     15 aug 1995     **/
/**                # Version 3.1  : from : 15 nov 1995     **/
/**                                 to     16 nov 1995     **/
/**                # Version 3.2  : from : 24 aug 1996     **/
/**                                 to   : 14 oct 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     19 oct 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     31 aug 2004     **/
/**                # Version 5.0  : from : 17 dec 2006     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 08 oct 2008     **/
/**                                 to     18 mar 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     25 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* bipartition graphs.   */
/*                       */
/*************************/

/* This routine builds the active graph
** corresponding to the given bipartitioning
** job parameters.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
bgraphInit (
Bgraph * restrict const         actgrafptr,       /*+ Active graph                  +*/
const Graph * restrict const    srcgrafptr,       /*+ Source graph                  +*/
const Arch * restrict const     archptr,          /*+ Target architecture           +*/
const ArchDom * restrict const  domnsubtab,       /*+ Array of the two subdomains   +*/
const Gnum * restrict const     vflowgttab)       /*+ Array of vertex weight biases +*/
{
  Anum                domndist;                   /* Distance between both subdomains   */
  Anum                domnwght0;                  /* Processor workforce in each domain */
  Anum                domnwght1;

  domndist  = archDomDist (archptr, &domnsubtab[0], &domnsubtab[1]); /* Get distance between subdomains */
  domnwght0 = archDomWght (archptr, &domnsubtab[0]); /* Get weights of subdomains                       */
  domnwght1 = archDomWght (archptr, &domnsubtab[1]);

  actgrafptr->s         = *srcgrafptr;            /* Get source graph data */
  actgrafptr->s.flagval = ((srcgrafptr->flagval & GRAPHBITSUSED) & ~GRAPHFREETABS) | BGRAPHFREEFRON | BGRAPHFREEPART; /* Graph is a clone with own grouped bipartitioning arrays */
  actgrafptr->s.vlbltax = NULL;                   /* Remove vertex labels    */
  actgrafptr->veextax   = NULL;                   /* No external gains (yet) */

  if (((actgrafptr->parttax = memAlloc (actgrafptr->s.vertnbr * sizeof (GraphPart))) == NULL) ||
      ((actgrafptr->frontab = memAlloc (actgrafptr->s.vertnbr * sizeof (Gnum)))      == NULL)) {
    errorPrint ("bgraphInit: out of memory");
    if (actgrafptr->parttax != NULL)
      memFree (actgrafptr->parttax);
    return (1);
  }
  actgrafptr->parttax -= actgrafptr->s.baseval;

  bgraphInit2 (actgrafptr, domndist, domnwght0, domnwght1, vflowgttab[0], vflowgttab[1]);

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (actgrafptr) != 0) {
    errorPrint ("bgraphInit: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}

/* This routine fills the internal fields of
** the given active graph.
** It returns:
** - void  : in all cases.
*/

void
bgraphInit2 (
Bgraph * restrict const         grafptr,          /*+ Active graph                       +*/
const Anum                      domndist,         /*+ Distance between both subdomains   +*/
const Anum                      domnwght0,        /*+ Processor workforce in each domain +*/
const Anum                      domnwght1,
const Gnum                      vfixload0,        /*+ Fixed load in each subdomain       +*/
const Gnum                      vfixload1)
{
  grafptr->fronnbr       = 0;                     /* No frontier vertices                         */
  grafptr->compload0min  = 0;                     /* No external constraints on bipartition (yet) */
  grafptr->compload0max  = grafptr->s.velosum;
  grafptr->compload0avg  = (Gnum) (((double) (grafptr->s.velosum + vfixload0 + vfixload1) * (double) domnwght0) / (double) (domnwght0 + domnwght1)) - vfixload0;
  grafptr->compload0dlt  = grafptr->s.velosum - grafptr->compload0avg;
  grafptr->compload0     = grafptr->s.velosum;
  grafptr->compsize0     = grafptr->s.vertnbr;    /* Size does not account for fixed vertices */
  grafptr->commload      = 0;
  grafptr->commloadextn0 = 0;
  grafptr->commgainextn  = 0;
  grafptr->commgainextn0 = 0;
  grafptr->domndist      = domndist;
  grafptr->domnwght[0]   = domnwght0;
  grafptr->domnwght[1]   = domnwght1;
  grafptr->vfixload[0]   = vfixload0;
  grafptr->vfixload[1]   = vfixload1;
  grafptr->bbalval       = (double) grafptr->compload0dlt / (double) grafptr->compload0avg;
  grafptr->levlnum       = 0;

  memSet (grafptr->parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */
}

/* This routine frees the contents
** of the given active graph.
** It returns:
** - void  : in all cases.
*/

void
bgraphExit (
Bgraph * restrict const     grafptr)
{
  if ((grafptr->veextax != NULL) &&               /* External gain array is private */
      ((grafptr->s.flagval & BGRAPHFREEVEEX) != 0))
    memFree (grafptr->veextax + grafptr->s.baseval);
  if ((grafptr->frontab != NULL) &&
      ((grafptr->s.flagval & BGRAPHFREEFRON) != 0))
    memFree (grafptr->frontab);
  if ((grafptr->parttax != NULL) &&
      ((grafptr->s.flagval & BGRAPHFREEPART) != 0))
    memFree (grafptr->parttax + grafptr->s.baseval);

  graphExit (&grafptr->s);                        /* Free re-allocated arrays of cloned source graph, if any */

#ifdef SCOTCH_DEBUG_BGRAPH2
  memSet (grafptr, ~0, sizeof (Bgraph));
#endif /* SCOTCH_DEBUG_BGRAPH2 */
}

/* This routine swaps all of the graph
** vertices from one part to another, and
** recomputes the resulting gains.
** It returns:
** - void  : in all cases.
*/

void
bgraphSwal (
Bgraph * restrict const     grafptr)
{
  Gnum                vertnnd;
  Gnum                vertnum;
  Gnum                comploadsum;

  GraphPart * restrict const  parttax = grafptr->parttax;

  for (vertnum = grafptr->s.baseval, vertnnd = grafptr->s.vertnnd; vertnum < vertnnd; vertnum ++)
    parttax[vertnum] ^= 1;

  comploadsum = grafptr->s.velosum + grafptr->vfixload[0] + grafptr->vfixload[1]; /* Overall load sum, including fixed vertex loads */
  grafptr->compload0    =   comploadsum - grafptr->compload0;
  grafptr->compload0dlt =   comploadsum - grafptr->compload0dlt - 2 * grafptr->compload0avg;
  grafptr->compsize0    =   grafptr->s.vertnbr - grafptr->compsize0; /* Size does not account for fixed vertices */
  grafptr->commload    +=   grafptr->commgainextn;
  grafptr->commgainextn = - grafptr->commgainextn;
}

/* This routine moves all of the graph
** vertices to the first part, and
** computes the resulting gains.
** It returns:
** - void  : in all cases.
*/

void
bgraphZero (
Bgraph * restrict const     grafptr)
{
  Gnum                compload0;

  compload0 = grafptr->s.velosum + grafptr->vfixload[0];
  grafptr->fronnbr      = 0;                      /* No frontier vertices */
  grafptr->compload0dlt = compload0 - grafptr->compload0avg;
  grafptr->compload0    = compload0;
  grafptr->compsize0    = grafptr->s.vertnbr;
  grafptr->commload     = grafptr->commloadextn0; /* Initialize communication load */
  grafptr->commgainextn = grafptr->commgainextn0;
  grafptr->bbalval      = (double) grafptr->compload0dlt / (double) grafptr->compload0avg;

  memSet (grafptr->parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */
}
