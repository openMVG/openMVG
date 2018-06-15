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
/**   NAME       : vgraph_separate_gp.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates an active         **/
/**                graph using a vertex-oriented version   **/
/**                of the Gibbs-Poole-Stockmeyer           **/
/**                algorithm.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 15 may 2004     **/
/**                                 to     17 may 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     12 sep 2007     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                # Version 6.0  : from : 10 feb 2011     **/
/**                                 to     10 feb 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_GP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_gp.h"

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
vgraphSeparateGp (
Vgraph * restrict const             grafptr,      /*+ Separation graph  +*/
const VgraphSeparateGpParam * const paraptr)      /*+ Method parameters +*/
{
  VgraphSeparateGpQueue             queudat;      /* Vertex queue               */
  VgraphSeparateGpVertex * restrict vexxtax;      /* Complementary vertex array */
  Gnum                              rootnum;
  Gnum                              vertnum;
  Gnum                              fronnum;
  Gnum                              compsize1;
  Gnum                              compsize2;
  Gnum                              compload2;
  Gnum                              comploaddlt;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  GraphPart * restrict const  parttax = grafptr->parttax;
  Gnum * restrict const       frontab = grafptr->frontab;

  if (grafptr->compload[0] != grafptr->s.velosum) /* If not all vertices already in part 0 */
    vgraphZero (grafptr);                         /* Move all graph vertices to part 0     */

  if (memAllocGroup ((void **) (void *)
                     &queudat.queutab, (size_t) (grafptr->s.vertnbr * sizeof (Gnum)),
                     &vexxtax,         (size_t) (grafptr->s.vertnbr * sizeof (VgraphSeparateGpVertex)), NULL) == NULL) {
    errorPrint ("vgraphSeparateGp: out of memory");
    return     (1);
  }
  memSet (vexxtax, 0, grafptr->s.vertnbr * sizeof (VgraphSeparateGpVertex)); /* Initialize pass numbers */
  vexxtax -= grafptr->s.baseval;

  compload2   = 0;                                /* All vertices to part 0 */
  comploaddlt = grafptr->s.velosum;
  for (rootnum = grafptr->s.baseval;              /* Loop on connected components */
       (rootnum < grafptr->s.vertnnd) && (comploaddlt > 0); rootnum ++) {
    Gnum                passnum;                  /* Pass number                                        */
    Gnum                diamnum;                  /* Number of current diameter vertex                  */
    Gnum                diamval;                  /* Current diameter value                             */
    Gnum                diamdeg;                  /* Degree of current diameter vertex                  */
    int                 diamflag;                 /* Flag set if improvement in diameter between passes */
    Gnum                veloval;

    while (vexxtax[rootnum].passnum != 0)         /* Find first unallocated vertex */
      rootnum ++;

    for (diamnum = rootnum, diamval = diamdeg = 0, diamflag = 1, passnum = 1; /* Start from root   */
         (passnum < paraptr->passnbr) && (diamflag -- != 0); passnum ++) { /* Loop if improvements */
      vgraphSeparateGpQueueFlush (&queudat);      /* Flush vertex queue                            */
      vgraphSeparateGpQueuePut   (&queudat, diamnum); /* Start from diameter vertex                */
      vexxtax[diamnum].passnum = passnum;         /* It has been enqueued                          */
      vexxtax[diamnum].distval = 0;

      do {                                        /* Loop on vertices in queue */
        Gnum                vertnum;
        Gnum                distval;
        Gnum                edgenum;

        vertnum = vgraphSeparateGpQueueGet (&queudat); /* Get vertex from queue */
        distval = vexxtax[vertnum].distval;       /* Get vertex distance        */

        if ((distval > diamval) ||                /* If vertex increases diameter         */
            ((distval == diamval) &&              /* Or is at diameter distance           */
             ((vendtax[vertnum] - verttax[vertnum]) < diamdeg))) { /* With smaller degree */
          diamnum  = vertnum;                     /* Set it as new diameter vertex        */
          diamval  = distval;
          diamdeg  = vendtax[vertnum] - verttax[vertnum];
          diamflag = 1;
        }

        distval ++;                               /* Set neighbor distance */
        for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
          Gnum                vertend;            /* End vertex number */

          vertend = edgetax[edgenum];
          if (vexxtax[vertend].passnum < passnum) { /* If vertex not yet queued      */
            vgraphSeparateGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex */
            vexxtax[vertend].passnum = passnum;
            vexxtax[vertend].distval = distval;
          }
        }
      } while (! vgraphSeparateGpQueueEmpty (&queudat)); /* As long as queue is not empty */
    }

    vgraphSeparateGpQueueFlush (&queudat);        /* Flush vertex queue           */
    vgraphSeparateGpQueuePut   (&queudat, diamnum); /* Start from diameter vertex */
    vexxtax[diamnum].passnum = passnum;           /* It has been enqueued         */
    vexxtax[diamnum].distval = 0;
    veloval = (velotax != NULL) ? velotax[diamnum] : 1;
    parttax[diamnum] = 2;                         /* Move diameter vertex to separator */
    comploaddlt -= veloval;
    compload2   += veloval;

    do {                                          /* Loop on vertices in queue */
      Gnum                vertnum;
      Gnum                veloval;
      Gnum                distval;
      Gnum                edgenum;

      vertnum = vgraphSeparateGpQueueGet (&queudat); /* Get vertex from queue */
      veloval = (velotax != NULL) ? velotax[vertnum] : 1;
      distval = vexxtax[vertnum].distval + 1;
      parttax[vertnum] = 1;                       /* Move selected vertex from separator to part 1 */
      comploaddlt -= veloval;
      compload2   -= veloval;

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;              /* End vertex number */
        Gnum                veloval;

        vertend = edgetax[edgenum];
        veloval = (velotax != NULL) ? velotax[vertend] : 1;
        if (vexxtax[vertend].passnum < passnum) { /* If vertex not yet queued      */
          vgraphSeparateGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex */
          vexxtax[vertend].passnum = passnum;
          vexxtax[vertend].distval = distval;
          parttax[vertend] = 2;                   /* Move neighbor vertex to separator */
          comploaddlt -= veloval;
          compload2   += veloval;
        }
      }
    } while ((comploaddlt > 0) && (! vgraphSeparateGpQueueEmpty (&queudat))); /* As long as balance not achieved and queue is not empty */
  }
  grafptr->compload[0] = (grafptr->s.velosum + comploaddlt - compload2) / 2;
  grafptr->compload[1] = grafptr->s.velosum - compload2 - grafptr->compload[0];
  grafptr->compload[2] = compload2;
  grafptr->comploaddlt = comploaddlt;

  memFree (queudat.queutab);                      /* Free group leader */

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
    errorPrint ("vgraphSeparateGp: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  return (0);
}
