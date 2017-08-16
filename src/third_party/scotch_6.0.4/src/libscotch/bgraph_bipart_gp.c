/* Copyright 2004,2007,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bipart_gp.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module bipartitions an active      **/
/**                graph using the Gibbs, Poole, and       **/
/**                Stockmeyer algorithm.                   **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 02 jun 1994     **/
/**                                 to     05 oct 1994     **/
/**                # Version 3.1  : from : 02 may 1996     **/
/**                                 to     02 may 1996     **/
/**                # Version 3.2  : from : 21 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     01 jun 2001     **/
/**                # Version 4.0  : from : 04 nov 2003     **/
/**                                 to     27 nov 2006     **/
/**                # Version 5.0  : from : 10 sep 2007     **/
/**                                 to     22 feb 2011     **/
/**                # Version 6.0  : from : 08 aug 2014     **/
/**                                 to     08 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH_BIPART_GP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_gp.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bgraphBipartGp (
Bgraph * restrict const           grafptr,
const BgraphBipartGpParam * const paraptr)        /*+ Method parameters +*/
{
  BgraphBipartGpQueue             queudat;        /* Neighbor queue               */
  BgraphBipartGpVertex * restrict vexxtax;        /* Complementary vertex array   */
  Gnum                            compload0dlt;
  Gnum                            compsize0;
  Gnum                            commloadintn;
  Gnum                            commloadextn;
  Gnum                            commgainextn;
  Gnum                            rootnum;        /* Index of potential next root */

  const Gnum * restrict const verttax = grafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * restrict const edlotax = grafptr->s.edlotax;
  const Gnum * restrict const veextax = grafptr->veextax;

  if (grafptr->compload0 != grafptr->s.velosum)   /* If not all vertices already in part 0 */
    bgraphZero (grafptr);                         /* Move all graph vertices to part 0     */

  if (memAllocGroup ((void **) (void *)
                     &queudat.queutab, (size_t) (grafptr->s.vertnbr * sizeof (Gnum)),
                     &vexxtax,         (size_t) (grafptr->s.vertnbr * sizeof (BgraphBipartGpVertex)), NULL) == NULL) {
    errorPrint ("bgraphBipartGp: out of memory");
    return     (1);
  }

  memSet (vexxtax, 0, grafptr->s.vertnbr * sizeof (BgraphBipartGpVertex)); /* Initialize pass numbers */
  vexxtax -= grafptr->s.baseval;

  compsize0    = grafptr->s.vertnbr;              /* All vertices in part zero */
  compload0dlt = grafptr->s.velosum - grafptr->compload0avg;
  commloadintn = 0;
  commloadextn = grafptr->commloadextn0;
  commgainextn = grafptr->commgainextn0;
  for (rootnum = grafptr->s.baseval;              /* Loop on connected components */
       (rootnum < grafptr->s.vertnnd) && (compload0dlt > 0); rootnum ++) {
    Gnum                passnum;                  /* Pass number                                        */
    Gnum                diamnum;                  /* Number of current diameter vertex                  */
    Gnum                diamval;                  /* Current diameter value                             */
    Gnum                diamdeg;                  /* Degree of current diameter vertex                  */
    int                 diamflag;                 /* Flag set if improvement in diameter between passes */

    while (vexxtax[rootnum].passnum != 0)         /* Find first unallocated vertex */
      rootnum ++;

    for (diamnum = rootnum, diamval = diamdeg = 0, diamflag = 1, passnum = 1; /* Start from root   */
         (passnum < paraptr->passnbr) && (diamflag -- != 0); passnum ++) { /* Loop if improvements */
      bgraphBipartGpQueueFlush (&queudat);        /* Flush vertex queue                            */
      bgraphBipartGpQueuePut   (&queudat, diamnum); /* Start from diameter vertex                  */
      vexxtax[diamnum].passnum = passnum;         /* It has been enqueued                          */
      vexxtax[diamnum].distval = 0;

      do {                                        /* Loop on vertices in queue */
        Gnum                vertnum;
        Gnum                distval;
        Gnum                edgenum;

        vertnum = bgraphBipartGpQueueGet (&queudat); /* Get vertex from queue */
        distval = vexxtax[vertnum].distval;       /* Get vertex distance      */

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
            bgraphBipartGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex */
            vexxtax[vertend].passnum = passnum;
            vexxtax[vertend].distval = distval;
          }
        }
      } while (! bgraphBipartGpQueueEmpty (&queudat)); /* As long as queue is not empty */
    }

    bgraphBipartGpQueueFlush (&queudat);          /* Flush vertex queue         */
    bgraphBipartGpQueuePut   (&queudat, diamnum); /* Start from diameter vertex */
    vexxtax[diamnum].passnum = passnum;           /* It has been enqueued       */
    vexxtax[diamnum].distval = 0;

    do {                                          /* Loop on vertices in queue */
      Gnum                vertnum;
      Gnum                veloval;
      Gnum                veexval;
      Gnum                distval;
      Gnum                edgenum;

      vertnum = bgraphBipartGpQueueGet (&queudat); /* Get vertex from queue */
      veloval = (velotax != NULL) ? velotax[vertnum] : 1;
      veexval = (veextax != NULL) ? veextax[vertnum] : 0;
      grafptr->parttax[vertnum] = 1;              /* Move selected vertex to part 1 */
      compsize0    --;
      compload0dlt -= veloval;
      commloadextn += veexval;
      commgainextn -= veexval * 2;

      distval = vexxtax[vertnum].distval + 1;
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;              /* End vertex number */

        vertend = edgetax[edgenum];
        if (vexxtax[vertend].passnum < passnum) { /* If vertex not yet queued    */
          bgraphBipartGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex */
          vexxtax[vertend].passnum = passnum;
          vexxtax[vertend].distval = distval;
        }
      }
    } while ((compload0dlt > 0) && (! bgraphBipartGpQueueEmpty (&queudat))); /* As long as balance not achieved and queue is not empty */

    if (! bgraphBipartGpQueueEmpty (&queudat)) {  /* If frontier non empty */
      Gnum                edloval;
      Gnum                fronnbr;

      fronnbr = 0;                                /* No frontier yet      */
      edloval = 1;                                /* Assume no edge loads */
      do {
        Gnum                vertnum;
        Gnum                edgenum;

        vertnum = bgraphBipartGpQueueGet (&queudat); /* Get vertex from queue */
        grafptr->frontab[fronnbr ++] = vertnum;
#ifdef SCOTCH_DEBUG_BGRAPH2
        if (grafptr->parttax[vertnum] != 0) {
          errorPrint ("bgraphBipartGp: internal error");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

        for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
          Gnum                vertend;            /* End vertex number */

          vertend = edgetax[edgenum];
          if (grafptr->parttax[vertend] == 1) {   /* If vertex belongs to other part */
            if (edlotax != NULL)
              edloval = edlotax[edgenum];
            commloadintn += edloval;
            if (vexxtax[vertend].distval != ~0) { /* If neighbor vertex not already put in frontier */
              grafptr->frontab[fronnbr ++] = vertend; /* Record it in frontier                      */
              vexxtax[vertend].distval = ~0;      /* Set it as recorded                             */
            }
          }
        }
      } while (! bgraphBipartGpQueueEmpty (&queudat));
      grafptr->fronnbr = fronnbr;
      break;                                      /* No need to process rest of graph */
    }                                             /* Else grafptr->fronnbr = 0 anyway */
  }

  grafptr->compload0    = grafptr->compload0avg + compload0dlt;
  grafptr->compload0dlt = compload0dlt;
  grafptr->compsize0    = compsize0;
  grafptr->commload     = commloadintn * grafptr->domndist + commloadextn;
  grafptr->commgainextn = commgainextn;
  grafptr->bbalval      = (double) ((grafptr->compload0dlt < 0) ? (- grafptr->compload0dlt) : grafptr->compload0dlt) / (double) grafptr->compload0avg;

  memFree (queudat.queutab);                      /* Free group leader */

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (grafptr) != 0) {
    errorPrint ("bgraphBipartGp: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
