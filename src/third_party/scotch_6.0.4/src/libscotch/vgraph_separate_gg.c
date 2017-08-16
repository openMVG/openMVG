/* Copyright 2004,2007,2008,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_gg.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates an active         **/
/**                graph using a vertex-oriented version   **/
/**                of the Greedy Graph Growing algorithm.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 10 nov 1997     **/
/**                                 to     15 jul 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     22 jan 2004     **/
/**                # Version 5.0  : from : 02 jan 2007     **/
/**                                 to     24 mar 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                # Version 6.0  : from : 04 feb 2012     **/
/**                                 to     04 feb 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_GG

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_gg.h"

/*
**  The static variables.
*/

static const Gnum           vgraphseparateggloadone = 1;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
vgraphSeparateGg (
Vgraph * restrict const             grafptr,      /*+ Separation graph  +*/
const VgraphSeparateGgParam * const paraptr)      /*+ Method parameters +*/
{
  GainTabl * restrict               tablptr;      /* Pointer to gain table             */
  VgraphSeparateGgVertex * restrict vexxtax;      /* Complementary vertex array        */
  Gnum                              vertnum;      /* Index of current vertex           */
  Gnum * restrict                   permtab;      /* Table for finding new roots       */
  Gnum                              permnum;      /* Current permutation index         */
  Gnum                              fronnum;
  INT                               passnum;
  const Gnum * restrict             velobax;      /* Data for handling optional arrays */
  Gnum                              velomsk;      /* Mask for handling optional arrays */
  Gnum                              comploaddlt;
  Gnum                              compload2;
  Gnum                              compsize1;
  Gnum                              compsize2;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * restrict const edlotax = grafptr->s.edlotax;
  GraphPart * restrict const  parttax = grafptr->parttax;
  Gnum * restrict const       frontab = grafptr->frontab;

  if (((tablptr = gainTablInit (GAIN_LINMAX, VGRAPHSEPAGGSUBBITS)) == NULL) || /* Use logarithmic array only */
      ((vexxtax = (VgraphSeparateGgVertex *) memAlloc (grafptr->s.vertnbr * sizeof (VgraphSeparateGgVertex))) == NULL)) {
    errorPrint ("vgraphSeparateGg: out of memory (1)");
    if (tablptr != NULL)
      gainTablExit (tablptr);
    return (1);
  }
  vexxtax -= grafptr->s.baseval;                  /* Base access to vexxtax                */
  permtab  = NULL;                                /* Do not allocate permutation array yet */

  if (grafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &vgraphseparateggloadone;           /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) { /* For all passes        */
    VgraphSeparateGgVertex *  vexxptr;            /* Pointer to current vertex to swap */

    memSet (vexxtax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (VgraphSeparateGgVertex)); /* All vertices to part 0 */
    gainTablFree (tablptr);                       /* Reset gain table            */
    permnum     = 0;                              /* No permutation built yet    */
    comploaddlt = grafptr->s.velosum;             /* Reset separation parameters */
    compload2   = 0;

    vexxptr = vexxtax + (grafptr->s.baseval + intRandVal (grafptr->s.vertnbr)); /* Randomly select first root vertex */

    do {                                          /* Loop on root vertices    */
      Gnum                        vertnum;        /* Number of current vertex */
      Gnum                        veloval;        /* Load of selected vertex  */
      Gnum                        compgain2;

      vexxptr->gainlink.next =                    /* TRICK: allow deletion of root vertex */
      vexxptr->gainlink.prev = (GainLink *) vexxptr;
#ifdef SCOTCH_DEBUG_GAIN2
      vexxptr->gainlink.tabl = NULL;
#endif /* SCOTCH_DEBUG_GAIN2 */

      vertnum = vexxptr - vexxtax;                /* Get root vertex based number */
      if (velomsk == 0) {                         /* If vertices are not weighted */
        veloval   = 1;
        compgain2 = vendtax[vertnum] - verttax[vertnum] - 1;
      }
      else {                                      /* Graph vertices are weighted */
        Gnum                        edgenum;

        veloval   = velobax[vertnum];
        compgain2 = - veloval;
        for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++)
          compgain2 += velobax[edgetax[edgenum]];
      }
      vexxptr->compgain2 = compgain2;             /* Set root gain (root not in separator) */
      comploaddlt -= veloval;                     /* Move vertex from part 0 to separator  */
      compload2   += veloval;

      do {                                        /* While vertices can be retrieved */
        VgraphSeparateGgVertex *    sepaptr;      /* List of vertices in separator   */
        Gnum                        veloval;      /* Load of selected vertex         */
        Gnum                        edgenum;

        vertnum = vexxptr - vexxtax;              /* Get number of selected vertex */
        veloval = velobax[vertnum & velomsk];

        if (comploaddlt < abs (comploaddlt - veloval)) { /* If swapping would cause imbalance */
          permnum = grafptr->s.vertnbr;           /* Terminate swapping process               */
          vexxptr = NULL;
          break;
        }
        gainTablDel (tablptr, (GainLink *) vexxptr); /* Remove vertex from table */
        vexxptr->gainlink.next = VGRAPHSEPAGGSTATEPART1; /* Put vertex in part 1 */
        compload2   += vexxptr->compgain2;        /* Update partition parameters */
        comploaddlt -= vexxptr->compgain2 + 2 * veloval;           

        sepaptr = NULL;                           /* No separator vertices to relink yet */
        for (edgenum = verttax[vertnum];          /* (Re-)link neighbor vertices         */
             edgenum < vendtax[vertnum]; edgenum ++) {
          Gnum                        vertend;
          VgraphSeparateGgVertex *    vexxend;

          vertend = edgetax[edgenum];             /* Point to end vertex */
          vexxend = vexxtax + vertend;
          if (vexxend->gainlink.next == VGRAPHSEPAGGSTATEPART0) { /* If end in part 0 */
            Gnum                veloend;
            Gnum                edgtnum;
            Gnum                compgain2;

            vexxend->gainlink.next = VGRAPHSEPAGGSTATEPART2; /* Move vertex to separator */
            vexxend->gainlink.prev = (GainLink *) sepaptr; /* Chain vertex               */
            sepaptr                = vexxend;

            veloend   = velobax[vertend & velomsk];
            compgain2 = - veloend;
            for (edgtnum = verttax[vertend];
                 edgtnum < vendtax[vertend]; edgtnum ++) {
              Gnum                        vertent;
              VgraphSeparateGgVertex *    vexxent;

              vertent = edgetax[edgtnum];         /* Point to end vertex */
              vexxent = vexxtax + vertent;
              if (vexxent->gainlink.next == VGRAPHSEPAGGSTATEPART0)
                compgain2 += velobax[vertent & velomsk];
              else if (vexxent->gainlink.next >= VGRAPHSEPAGGSTATEPART2) {
                vexxent->compgain2 -= veloend;
                if (vexxent->gainlink.next >= VGRAPHSEPAGGSTATELINK) {
                  gainTablDel (tablptr, (GainLink *) vexxent); /* Unlink vertex    */
                  vexxent->gainlink.next = VGRAPHSEPAGGSTATEPART2; /* Chain vertex */
                  vexxent->gainlink.prev = (GainLink *) sepaptr;
                  sepaptr                = vexxent;
                }
              }
            }
            vexxend->compgain2 = compgain2;
          }
        }
        while (sepaptr != NULL) {                 /* For all vertices in chain list */
          vexxptr = sepaptr;                      /* Unlink vertex from list        */
          sepaptr = (VgraphSeparateGgVertex *) vexxptr->gainlink.prev;
          gainTablAdd (tablptr, (GainLink *) vexxptr, vexxptr->compgain2); /* Relink it */
        }
      } while ((vexxptr = (VgraphSeparateGgVertex *) gainTablFrst (tablptr)) != NULL);

      if (permnum == 0) {                         /* If permutation has not been built yet  */
        if (permtab == NULL) {                    /* If permutation array not allocated yet */
          if ((permtab = (Gnum *) memAlloc (grafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
            errorPrint   ("vgraphSeparateGg: out of memory (2)");
            memFree      (vexxtax + grafptr->s.baseval);
            gainTablExit (tablptr);
            return       (1);
          }
          intAscn (permtab, grafptr->s.vertnbr, grafptr->s.baseval); /* Initialize based permutation array */
        }
        intPerm (permtab, grafptr->s.vertnbr);    /* Build random permutation */
      }
      for ( ; permnum < grafptr->s.vertnbr; permnum ++) { /* Find next root vertex */
        if (vexxtax[permtab[permnum]].gainlink.next == VGRAPHSEPAGGSTATEPART0) {
          vexxptr = vexxtax + permtab[permnum ++];
          break;
        }
      }
    } while (vexxptr != NULL);

    if ((passnum == 0) ||                         /* If first try                  */
        ( (grafptr->compload[2] >  compload2) ||  /* Or if better solution reached */
         ((grafptr->compload[2] == compload2) &&
          (abs (grafptr->comploaddlt) > abs (comploaddlt))))) {
      Gnum                vertnum;

      grafptr->comploaddlt = comploaddlt;         /* Set graph parameters */
      grafptr->compload[2] = compload2;

      for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) /* Copy bipartition state */
        parttax[vertnum] = (vexxtax[vertnum].gainlink.next <= VGRAPHSEPAGGSTATEPART2) ? (GraphPart) (intptr_t) vexxtax[vertnum].gainlink.next : (GraphPart) 2;
    }
  }

  if (permtab != NULL)                            /* Free work arrays */
    memFree (permtab);
  memFree      (vexxtax + grafptr->s.baseval);
  gainTablExit (tablptr);

  grafptr->compload[0] = (grafptr->s.velosum + grafptr->comploaddlt - grafptr->compload[2]) / 2;
  grafptr->compload[1] = grafptr->s.velosum - grafptr->compload[2] - grafptr->compload[0];
  compsize1 =
  compsize2 = 0;
  for (vertnum = grafptr->s.baseval, fronnum = 0;
       vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                partval;

    partval    = (Gnum) parttax[vertnum];
    compsize1 += (partval & 1);                   /* Superscalar update */
    compsize2 += (partval >> 1);
    if (partval == 2)                             /* If vertex belongs to frontier */
      frontab[fronnum ++] = vertnum;              /* Record it in frontier array   */
  }
  grafptr->compsize[0] = grafptr->s.vertnbr - compsize1 - compsize2;
  grafptr->compsize[1] = compsize1;
  grafptr->fronnbr     = compsize2;

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (grafptr) != 0) {
    errorPrint ("vgraphSeparateGg: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  return (0);
}
