/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hgraph_order_gp.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders a subgraph (most     **/
/**                likely a separator) using the Gibbs,    **/
/**                Poole, and Stockmeyer algorithm.        **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 31 oct 1996     **/
/**                                 to     27 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 02 oct 1998     **/
/**                # Version 4.0  : from : 28 jun 2002     **/
/**                                 to   : 01 dec 2003     **/
/**                # Version 4.0  : from : 10 sep 2007     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 01 oct 2009     **/
/**                                 to   : 01 oct 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_GP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_gp.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the ordering.
** It returns:
** - 0   : if the ordering could be computed.
** - !0  : on error.
*/

int
hgraphOrderGp (
const Hgraph * restrict const             grafptr,
Order * restrict const                    ordeptr,
const Gnum                                ordenum,
OrderCblk * restrict const                cblkptr, /*+ Single column-block +*/
const HgraphOrderGpParam * restrict const paraptr)
{
  HgraphOrgerGpQueue              queudat;        /* Neighbor queue                 */
  HgraphOrderGpVertex * restrict  vexxtax;        /* Based access to vertex array   */
  Gnum                            passnum;        /* Pass number                    */
  Gnum                            rootnum;        /* Number of root vertex          */
  Gnum                            diamnum;        /* Vertex which achieves diameter */
  int                             diamflag;       /* Flag set if diameter changed   */
  Gnum                            diamdist;       /* Maximum diameter value found   */
  Gnum                            vertdist;       /* DIstance of current vertex     */
  Gnum                            vertnum;        /* Number of current vertex       */
  Gnum                            edgenum;        /* Number of current edge         */
  Gnum                            ordeval;        /* Current ordering value         */
  Gnum                            ordevnd;        /* End value of ordering          */

  const Gnum * restrict const     verttax = grafptr->s.verttax;
  const Gnum * restrict const     vnumtax = grafptr->s.vnumtax;
  const Gnum * restrict const     vnhdtax = grafptr->vnhdtax;
  const Gnum * restrict const     edgetax = grafptr->s.edgetax;

  if (memAllocGroup ((void **) (void *)
                     &queudat.qtab, (size_t) (grafptr->vnohnbr * sizeof (Gnum)),
                     &vexxtax,      (size_t) (grafptr->vnohnbr * sizeof (HgraphOrderGpVertex)), NULL) == NULL) {
    errorPrint ("hgraphOrderGp: out of memory");
    return     (1);
  }
  memSet (vexxtax, 0, grafptr->vnohnbr * sizeof (HgraphOrderGpVertex)); /* Initialize pass numbers */
  vexxtax -= grafptr->s.baseval;

#ifdef SCOTCH_DEBUG_ORDER2
  memSet (ordeptr->peritab + ordenum, ~0, grafptr->vnohnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_ORDER2 */

  for (ordeval = ordenum, rootnum = grafptr->s.baseval, /* For all connected components */
       ordevnd = ordeval + grafptr->vnohnbr;
       ordeval < ordevnd; ) {
    while (vexxtax[rootnum].passnum != 0) {       /* Find first unallocated root */
      rootnum ++;
#ifdef SCOTCH_DEBUG_ORDER2
      if (rootnum >= grafptr->vnohnnd) {
        errorPrint ("hgraphOrderGp: internal error (1)");
        memFree    (queudat.qtab);                /* Free group leader */
        return     (1);
      }
#endif /* SCOTCH_DEBUG_ORDER2 */
    }

    diamnum  = rootnum;                           /* Start from found root */
    diamdist = 0;
    for (diamflag = 0, passnum = 1;               /* Loop if modifications */
         (diamflag ++ == 0) && (passnum <= paraptr->passnbr); passnum ++) {
      Gnum                  diamdegr;             /* Degree of current pseudo-peripherial vertex */

      hgraphOrderGpQueueFlush (&queudat);         /* Flush vertex queue          */
      hgraphOrderGpQueuePut   (&queudat, diamnum); /* Start from diameter vertex */
      vexxtax[diamnum].passnum  = passnum;        /* It has been enqueued        */
      vexxtax[diamnum].vertdist = 0;              /* It is at distance zero      */
      diamdegr = vnhdtax[diamnum] - verttax[diamnum];

      do {                                        /* Loop on vertices in queue */
        vertnum = hgraphOrderGpQueueGet (&queudat); /* Get vertex from queue   */
#ifdef SCOTCH_DEBUG_ORDER2
        if ((vertnum < grafptr->s.baseval) || (vertnum >= grafptr->vnohnnd)) {
          errorPrint ("hgraphOrderGp: internal error (2)");
          memFree    (queudat.qtab);              /* Free group leader */
          return     (1);
        }
#endif /* SCOTCH_DEBUG_ORDER2 */
        vertdist = vexxtax[vertnum].vertdist;     /* Get vertex distance */

        if ((vertdist > diamdist) ||              /* If vertex increases diameter          */
            ((vertdist == diamdist) &&            /* Or is at diameter distance            */
             ((vnhdtax[vertnum] - verttax[vertnum]) < diamdegr))) { /* With smaller degree */
          diamnum  = vertnum;                     /* Set it as new diameter vertex         */
          diamdist = vertdist;
          diamdegr = vnhdtax[vertnum] - verttax[vertnum];
          diamflag = 0;
        }

        vertdist ++;                              /* Set neighbor distance */
        for (edgenum = verttax[vertnum]; edgenum < vnhdtax[vertnum]; edgenum ++) {
          Gnum                  vertend;

          vertend = edgetax[edgenum];
#ifdef SCOTCH_DEBUG_ORDER2
          if ((vertend < grafptr->s.baseval) || (vertend >= grafptr->vnohnnd)) {
            errorPrint ("hgraphOrderGp: internal error (3)");
            memFree    (queudat.qtab);              /* Free group leader */
            return     (1);
          }
#endif /* SCOTCH_DEBUG_ORDER2 */

          if (vexxtax[vertend].passnum < passnum) { /* If vertex not queued yet   */
            hgraphOrderGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex */
            vexxtax[vertend].passnum  = passnum;
            vexxtax[vertend].vertdist = vertdist;
          }
        }
      } while (! hgraphOrderGpQueueEmpty (&queudat)); /* As long as queue is not empty */
    }

    hgraphOrderGpQueueFlush (&queudat);           /* Flush vertex queue         */
    hgraphOrderGpQueuePut   (&queudat, diamnum);  /* Start from diameter vertex */
    vexxtax[diamnum].passnum = passnum;           /* Vertex has been enqueued   */

    do {                                          /* Loop on vertices in queue */
      vertnum = hgraphOrderGpQueueGet (&queudat); /* Get vertex from queue     */

      if (vexxtax[vertnum].passnum > passnum)     /* If vertex already ordered (by-level ordering) */
        continue;                                 /* Skip to next vertex in queue                  */

      vertdist = vexxtax[vertnum].vertdist;       /* Get vertex distance       */
      do {                                        /* Loop on vertices in layer */
        Gnum                  edgennd;            /* End of edge sub-array     */

        ordeptr->peritab[ordeval ++] = (vnumtax == NULL) ? vertnum : vnumtax[vertnum];
        vexxtax[vertnum].passnum = passnum + 1;   /* Set vertex as ordered */

        for (edgenum = verttax[vertnum], edgennd = vnhdtax[vertnum], vertnum = ~0;
             edgenum < edgennd; edgenum ++) {     /* Need edgennd because vertnum is overwritten */
          Gnum                  vertend;

          vertend = edgetax[edgenum];

          if ((vexxtax[vertend].vertdist == vertdist) && /* If neighbor vertex in same layer        */
              (vexxtax[vertend].passnum <= passnum)) { /* And not yet ordered                       */
            vertnum = vertend;                    /* Set neighbor as next vertex                    */
            edgenum ++;                           /* Process next neighbors, not this one again     */
            break;                                /* Process next neighbors without further testing */
          }
          if (vexxtax[vertend].passnum < passnum) { /* Else if vertex not yet enqueued */
            hgraphOrderGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex      */
            vexxtax[vertend].passnum = passnum; /* Set it as enqueued                  */
          }
        }
        for ( ; edgenum < edgennd; edgenum ++) {  /* Enqueue remaining neighbors */
          Gnum                  vertend;

          vertend = edgetax[edgenum];

          if (vexxtax[vertend].passnum < passnum) { /* If neighbor not yet enqueued */
            hgraphOrderGpQueuePut (&queudat, vertend); /* Enqueue neighbor vertex   */
            vexxtax[vertend].passnum = passnum; /* Set it as enqueued               */
          }
        }
      } while (vertnum != ~0);
    } while (! hgraphOrderGpQueueEmpty (&queudat)); /* As long as queue is not empty */
  }

#ifdef SCOTCH_DEBUG_ORDER2
  for (ordeval = ordenum; ordeval < ordenum + grafptr->vnohnbr; ordeval ++) {
    if (ordeptr->peritab[ordeval] == ~0) {
      errorPrint ("hgraphOrderGp: internal error (4)");
      memFree    (queudat.qtab);                  /* Free group leader */
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  memFree (queudat.qtab);                         /* Group freeing */

  return (0);
}
