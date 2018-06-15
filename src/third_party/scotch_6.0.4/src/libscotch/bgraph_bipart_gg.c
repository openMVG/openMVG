/* Copyright 2004,2007,2009,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_gg.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Luca SCARANO (v3.1)                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of   **/
/**                a bipartition graph by multiple runs of **/
/**                the greedy graph growing algorithm.     **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 jan 1996     **/
/**                                 to     07 jun 1996     **/
/**                # Version 3.2  : from : 20 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     01 jun 2001     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     01 sep 2004     **/
/**                # Version 5.0  : from : 02 jan 2007     **/
/**                                 to     04 feb 2007     **/
/**                # Version 5.1  : from : 21 nov 2007     **/
/**                                 to     22 feb 2011     **/
/**                # Version 6.0  : from : 23 fev 2011     **/
/**                                 to     08 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH_BIPART_GG

#define SCOTCH_TABLE_GAIN

#include "module.h"
#include "common.h"
#include "gain.h"
#include "fibo.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_gg.h"

/*
**  The static variables.
*/

static const Gnum           bgraphbipartggloadone  = 1;
static const Gnum           bgraphbipartggloadzero = 0;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

#ifndef SCOTCH_TABLE_GAIN

/* bgraphBipartFmCmpFunc(a,b) must return a negative
** number if a is "better" than b. The smaller, the
** better.
*/

static
int
bgraphBipartGgCmpFunc (
const FiboNode * const      data0ptr,             /* TRICK: BgraphBipartFmLink is FIRST in BgraphBipartFmVertex */
const FiboNode * const      data1ptr)
{
  const BgraphBipartGgVertex * const  node0ptr = (BgraphBipartGgVertex *) data0ptr;
  const BgraphBipartGgVertex * const  node1ptr = (BgraphBipartGgVertex *) data1ptr;

  if (node0ptr->commgain < node1ptr->commgain)
    return (-1);
  if (node0ptr->commgain > node1ptr->commgain)
    return (1);
  return (0);
}
 
#endif /* SCOTCH_TABLE_GAIN */

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bgraphBipartGg (
Bgraph * restrict const           grafptr,        /*+ Active graph      +*/
const BgraphBipartGgParam * const paraptr)        /*+ Method parameters +*/
{
  BgraphBipartGgTabl      tabldat;                /* Gain table                                   */
  BgraphBipartGgVertex *  vexxtax;                /* Extended vertex array [norestrict]           */
  BgraphBipartGgVertex *  vexxptr;                /* Pointer to current vertex to swap            */
  const Gnum * restrict   veexptr;                /* Pointer to external gain of current vertex   */
  Gnum * restrict         permtab;                /* Permutation table for finding new roots      */
  Gnum                    permnum;                /* Current permutation index                    */
  const Gnum * restrict   velobax;                /* Data for handling of optional arrays         */
  Gnum                    velomsk;
  const byte * restrict   veexbab;                /* Un-based array for external gains            */
  int                     veexsiz;
  const Gnum *            edlobax;                /* Pointer to array or dummy value [norestrict] */
  Gnum                    edlomsk;
  byte * restrict         flagtax;
  Gnum                    vertnum;
  Gnum                    fronnum;
  Gnum                    compsize1;
  Gnum                    commgainextn;
  unsigned int            passnum;
  Anum                    domndist;
  Anum                    domndist2;              /* Two times domndist */

  const Gnum * restrict const verttax = grafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * const          edlotax = grafptr->s.edlotax; /* [norestrict] */
  const Gnum * restrict const veextax = grafptr->veextax;

  if ((bgraphBipartGgTablInit (&tabldat) != 0) ||
      ((vexxtax = (BgraphBipartGgVertex *) memAlloc (grafptr->s.vertnbr * sizeof (BgraphBipartGgVertex))) == NULL)) {
    errorPrint ("bgraphBipartGg: out of memory (1)");
    bgraphBipartGgTablExit (&tabldat);
    return (1);
  }
  vexxtax -= grafptr->s.baseval;                  /* Base access to vexxtax                */
  permtab  = NULL;                                /* Do not allocate permutation array yet */

  domndist  = grafptr->domndist;
  domndist2 = grafptr->domndist * 2;

  if (edlotax == NULL) {                          /* If graph has no edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                commload;

      commload = (vendtax[vertnum] - verttax[vertnum]) * domndist;
      vexxtax[vertnum].commgain0 = (veextax == NULL) ? commload : commload + veextax[vertnum];
    }

    edlobax = &bgraphbipartggloadone;
    edlomsk = 0;
  }
  else {                                          /* Graph has edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                commload;
      Gnum                edgenum;

      for (edgenum = verttax[vertnum], commload = 0;
           edgenum < vendtax[vertnum]; edgenum ++)
        commload += edlotax[edgenum];
      commload *= domndist;

      vexxtax[vertnum].commgain0 = (veextax == NULL) ? commload : commload + veextax[vertnum];
    }

    edlobax = edlotax;
    edlomsk = ~((Gnum) 0);                        /* TRICK: assume that ~0 is -1 */
  }
  if (grafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &bgraphbipartggloadone;             /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }
  if (veextax == NULL) {
    veexbab = (byte *) &bgraphbipartggloadzero;
    veexsiz = 0;
  }
  else {
    veexbab = (byte *) (veextax + grafptr->s.baseval);
    veexsiz = sizeof (Gnum);
  }

  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) { /* For all passes */
    Gnum                vertnum;
    Gnum                commload;
    Gnum                compload0dlt;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) { /* Reset extended vertex array */
      bgraphBipartGgSetFree (&vexxtax[vertnum]);
      vexxtax[vertnum].commgain = vexxtax[vertnum].commgain0;
    }
    bgraphBipartGgTablFree (&tabldat);            /* Reset gain table                          */
    permnum      = 0;                             /* No permutation built yet                  */
    compload0dlt = grafptr->s.velosum - grafptr->compload0avg; /* Reset bipartition parameters */
    commload     = grafptr->commloadextn0;

    vexxptr = vexxtax + (grafptr->s.baseval + intRandVal (grafptr->s.vertnbr)); /* Randomly select first root vertex */

    do {                                          /* For all root vertices, till balance  */
#ifdef SCOTCH_TABLE_GAIN
      vexxptr->gainlink.next =                    /* TRICK: allow deletion of root vertex */
      vexxptr->gainlink.prev = (GainLink *) vexxptr;
#ifdef SCOTCH_DEBUG_BGRAPH2
      vexxptr->gainlink.tabl = NULL;
#endif /* SCOTCH_DEBUG_BGRAPH2 */
#endif /* SCOTCH_TABLE_GAIN    */

      do {                                        /* As long as vertices can be retrieved */
        const Gnum * restrict       edgeptr;      /* Pointer to current end vertex index  */
        const Gnum * restrict       edgetnd;      /* Pointer to end of edge array         */
        const Gnum * restrict       edloptr;      /* Pointer to current edge load         */
        Gnum                        vertnum;      /* Number of current vertex             */
        Gnum                        veloval;      /* Load of selected vertex              */

#ifndef SCOTCH_TABLE_GAIN
        if (bgraphBipartGgIsTabl (vexxptr))
#endif /* SCOTCH_TABLE_GAIN */
        bgraphBipartGgTablDel (&tabldat, vexxptr); /* Remove vertex from table */

        vertnum = vexxptr - vexxtax;              /* Get number of selected vertex */
        veloval = velobax[vertnum & velomsk];

        if (abs (compload0dlt - veloval) >= abs (compload0dlt)) { /* If swapping would cause imbalance */
#ifndef SCOTCH_TABLE_GAIN
          bgraphBipartGgNext (vexxptr) = BGRAPHBIPARTGGSTATELINK; /* Vertex belongs to frontier of part 0 */
#endif /* SCOTCH_TABLE_GAIN */
          permnum = grafptr->s.vertnbr;           /* Terminate swapping process */
          vexxptr = NULL;
          break;
        }

        bgraphBipartGgSetUsed (vexxptr);          /* Mark it as swapped          */
        compload0dlt -= veloval;                  /* Update partition parameters */
        commload     += vexxptr->commgain;
        for (edgeptr = edgetax + verttax[vertnum], /* (Re-)link neighbors */
             edgetnd = edgetax + vendtax[vertnum],
             edloptr = edlobax + (verttax[vertnum] & edlomsk);
             edgeptr < edgetnd; edgeptr ++, edloptr -= edlomsk) { /* TRICK: assume that ~0 is -1 */
          BgraphBipartGgVertex *  vexxend;        /* Pointer to end vertex of current edge       */

          vexxend = vexxtax + *edgeptr;           /* Point to end vertex            */
          if (! bgraphBipartGgIsUsed (vexxend)) { /* If vertex needs to be updated  */
            vexxend->commgain -= *edloptr * domndist2; /* Adjust gain value         */
            if (bgraphBipartGgIsTabl (vexxend))   /* If vertex is linked            */
              bgraphBipartGgTablDel (&tabldat, vexxend); /* Remove it from table    */
            bgraphBipartGgTablAdd (&tabldat, vexxend); /* (Re-)link vertex in table */
          }
        }
      } while ((vexxptr = (BgraphBipartGgVertex *) bgraphBipartGgTablFrst (&tabldat)) != NULL);

      if (permnum == 0) {                         /* If permutation has not been built yet  */
        if (permtab == NULL) {                    /* If permutation array not allocated yet */
          if ((permtab = (Gnum *) memAlloc (grafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
            errorPrint ("bgraphBipartGg: out of memory (3)");
            memFree    (vexxtax + grafptr->s.baseval);
            bgraphBipartGgTablExit (&tabldat);
            return (1);
          }
          intAscn (permtab, grafptr->s.vertnbr, grafptr->s.baseval); /* Initialize based permutation array */
        }
        intPerm (permtab, grafptr->s.vertnbr);    /* Build random permutation */
      }
      for ( ; permnum < grafptr->s.vertnbr; permnum ++) { /* Find next root vertex */
        if (bgraphBipartGgIsFree (&vexxtax[permtab[permnum]])) {
          vexxptr = vexxtax + permtab[permnum ++];
          break;
        }
      }
    } while (vexxptr != NULL);

    if ((passnum == 0) ||                         /* If first try                  */
        ( (grafptr->commload >  commload) ||      /* Or if better solution reached */
         ((grafptr->commload == commload) &&
          (abs (grafptr->compload0dlt) > abs (compload0dlt))))) {
      Gnum                vertnum;

      grafptr->compload0dlt = compload0dlt;       /* Set graph parameters */
      grafptr->commload     = commload;

      for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) /* Copy bipartition state with flag 2 for tabled vertices */
        grafptr->parttax[vertnum] = (bgraphBipartGgIsTabl (&vexxtax[vertnum])) ? 2 : (GraphPart) ((intptr_t) bgraphBipartGgNext (&vexxtax[vertnum]));
    }
  }

  flagtax = (byte *) (vexxtax + grafptr->s.baseval) - grafptr->s.baseval; /* Re-use extended vertex array for flag array */
  memSet (flagtax + grafptr->s.baseval, ~0, grafptr->s.vertnbr * sizeof (byte));
  for (vertnum = grafptr->s.baseval, veexptr = (Gnum *) veexbab, fronnum = 0, compsize1 = 0, commgainextn = grafptr->commgainextn0;
       vertnum < grafptr->s.vertnnd; vertnum ++, veexptr = (Gnum *) ((byte *) veexptr + veexsiz)) {
    int                 partval;

    partval = grafptr->parttax[vertnum];
    if (partval > 1) {                            /* If vertex belongs to frontier of part 0 */
      Gnum                edgenum;
      Gnum                frontmp;                /* Temporary count value for frontier */

      grafptr->frontab[fronnum ++] = vertnum;     /* Then it belongs to the frontier */
      grafptr->parttax[vertnum]    = 0;           /* And it belongs to part 0        */
      for (edgenum = verttax[vertnum], frontmp = 1;
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = edgetax[edgenum];
        if (grafptr->parttax[vertend] == 1) {     /* If vertex belongs to other part       */
          frontmp = 0;                            /* Then first frontier vertex was useful */
          if (flagtax[vertend] != 0) {            /* If vertex has not yet been flagged    */
            grafptr->frontab[fronnum ++] = vertend; /* Then add it to the frontier         */
            flagtax[vertend] = 0;                 /* Flag it                               */
          }
        }
      }
      fronnum -= frontmp;                         /* Remove vertex from frontier if it was useless */
    }
    partval      &= 1;
    compsize1    += partval;
    commgainextn -= partval * 2 * *veexptr;
  }
  grafptr->fronnbr      = fronnum;
  grafptr->compload0    = grafptr->compload0avg + grafptr->compload0dlt;
  grafptr->compsize0    = grafptr->s.vertnbr - compsize1;
  grafptr->commgainextn = commgainextn;
  grafptr->bbalval      = (double) ((grafptr->compload0dlt < 0) ? (- grafptr->compload0dlt) : grafptr->compload0dlt) / (double) grafptr->compload0avg;

  if (permtab != NULL)                            /* Free work arrays */
    memFree (permtab);
  memFree (vexxtax + grafptr->s.baseval);
  bgraphBipartGgTablExit (&tabldat);

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (grafptr) != 0) {
    errorPrint ("bgraphBipartGg: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
