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
/**   NAME       : hmesh_order_gp.c                        **/
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
/**                # Version 4.0  : from : 05 nov 2002     **/
/**                                 to   : 27 jan 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to   : 12 sep 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_ORDER_GP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"
#include "mesh.h"
#include "hmesh.h"
#include "hmesh_order_gp.h"

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
hmeshOrderGp (
const Hmesh * restrict const              meshptr,
Order * restrict const                    ordeptr,
const Gnum                                ordenum,
OrderCblk * restrict const                cblkptr, /*+ Single column-block +*/
const HmeshOrderGpParam * restrict const  paraptr)
{
  HmeshOrderGpQueue             queue;            /* Neighbor queue                 */
  HmeshOrderGpVertex * restrict vexxtax;          /* Based access to vertex array   */
  HmeshOrderGpVertex *          rootptr;          /* Pointer to root vertex         */
  Gnum                          passnum;          /* Pass number                    */
  int                           passflag;         /* Flag set if diameter changed   */
  Gnum                          vdianum;          /* Vertex which achieves diameter */
  Gnum                          vdiadist;         /* Maximum diameter value found   */
  Gnum                          vnodnbr;          /* Number of vertices found yet   */
  Gnum                          ordeval;          /* Current ordering value         */

  if (memAllocGroup ((void **) (void *)
        &queue.qtab, (size_t) ((meshptr->vnohnnd - meshptr->m.baseval)   * sizeof (Gnum)),
        &vexxtax,    (size_t) ((meshptr->m.velmnbr + meshptr->m.vnodnbr) * sizeof (HmeshOrderGpVertex)), NULL) == NULL) {
    errorPrint ("hmeshOrderGp: out of memory");
    return     (1);
  }
  vexxtax -= meshptr->m.baseval;                /* Base vexxtab array */
  memSet (vexxtax + meshptr->m.velmbas, 0, meshptr->m.velmnbr                      * sizeof (HmeshOrderGpVertex)); /* Initialize pass numbers for */
  memSet (vexxtax + meshptr->m.vnodbas, 0, (meshptr->vnohnnd - meshptr->m.vnodbas) * sizeof (HmeshOrderGpVertex)); /* All but halo node vertices  */

  for (vnodnbr = 0, ordeval = ordenum, rootptr = vexxtax + meshptr->m.vnodbas, passnum = 1; /* For all connected components */
       vnodnbr < meshptr->vnohnbr; passnum ++) {
    while (rootptr->passnum != 0)                 /* Find first unallocated root */
      rootptr ++;

    vdianum  = rootptr - vexxtax;                 /* Start from found root */
    vdiadist = 0;
    for (passflag = 1; (passflag -- != 0) && (passnum <= paraptr->passnbr); passnum ++) { /* Loop if modifications */
      hmeshOrderGpQueueFlush (&queue);            /* Flush vertex queue         */
      hmeshOrderGpQueuePut   (&queue, vdianum);   /* Start from diameter vertex */
      vexxtax[vdianum].passnum  = passnum;        /* It has been enqueued       */
      vexxtax[vdianum].vertdist = 0;              /* It is at distance zero     */

      do {                                        /* Loop on vertices in queue  */
        Gnum                  vnodnum;            /* Number of current vertex   */
        Gnum                  vnoddist;           /* Distance of current vertex */
        Gnum                  enodnum;

        vnodnum  = hmeshOrderGpQueueGet (&queue); /* Get vertex from queue */
        vnoddist = vexxtax[vnodnum].vertdist;     /* Get vertex distance   */

        if ((vnoddist > vdiadist) ||              /* If vertex increases diameter                  */
            ((vnoddist == vdiadist) &&            /* Or is at diameter distance                    */
             ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]) < /* With smaller degree */
              (meshptr->m.vendtax[vdianum] - meshptr->m.verttax[vdianum])))) {
          vdianum  = vnodnum;                     /* Set it as new diameter vertex */
          vdiadist = vnoddist;
          passflag = 1;
        }

        vnoddist ++;                              /* Set neighbor distance */

        for (enodnum = meshptr->m.verttax[vnodnum]; enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
          Gnum                  velmnum;
          Gnum                  eelmnum;

          velmnum = meshptr->m.edgetax[enodnum];  /* Get neighboring element */

          if (vexxtax[velmnum].passnum >= passnum) /* If element already scanned */
            continue;                             /* Skip to next element        */

          vexxtax[velmnum].passnum = passnum;     /* Set element as scanned */

          for (eelmnum = meshptr->m.verttax[velmnum]; /* For all neighboring non-halo nodes */
               eelmnum < meshptr->vehdtax[velmnum]; eelmnum ++) {
            Gnum                  vnodend;        /* Neighboring node */

            vnodend = meshptr->m.edgetax[eelmnum]; /* Get neighboring node */

            if (vexxtax[vnodend].passnum < passnum) { /* If node vertex not yet enqueued */
              hmeshOrderGpQueuePut (&queue, vnodend); /* Enqueue neighbor vertex         */
              vexxtax[vnodend].passnum  = passnum;
              vexxtax[vnodend].vertdist = vnoddist;
            }
          }
        }
      } while (! hmeshOrderGpQueueEmpty (&queue)); /* As long as queue is not empty */
    }

    hmeshOrderGpQueueFlush (&queue);              /* Flush vertex queue         */
    hmeshOrderGpQueuePut   (&queue, vdianum);     /* Start from diameter vertex */
    vexxtax[vdianum].passnum = passnum;           /* It has been enqueued       */

    do {                                          /* Loop on vertices in queue  */
      Gnum                  vnodnum;              /* Number of current vertex   */
      Gnum                  vnoddist;             /* Distance of current vertex */

      vnodnum = hmeshOrderGpQueueGet (&queue);    /* Get vertex from queue        */
      if (vexxtax[vnodnum].passnum > passnum)     /* If vertex already ordered    */
        continue;                                 /* Skip to next vertex in queue */

      vnoddist = vexxtax[vnodnum].vertdist;       /* Get vertex distance       */
      do {                                        /* Loop on vertices in layer */
        Gnum                  enodnum;
        Gnum                  enodnnd;

        ordeptr->peritab[ordeval] = (meshptr->m.vnumtax == NULL) /* Order node vertex */
                                    ? vnodnum - (meshptr->m.vnodbas - meshptr->m.baseval)
                                    : meshptr->m.vnumtax[vnodnum];
#ifdef SCOTCH_DEBUG_ORDER2
        if ((ordeptr->peritab[ordeval] <   ordeptr->baseval) ||
            (ordeptr->peritab[ordeval] >= (ordeptr->baseval + ordeptr->vnodnbr))) {
          errorPrint ("hmeshOrderGp: invalid permutation index");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_ORDER2 */
        ordeval ++;

        vexxtax[vnodnum].passnum = (passnum + 1); /* Set vertex as ordered */
        vnodnbr ++;                               /* Count it              */

        for (enodnum = meshptr->m.verttax[vnodnum], enodnnd = meshptr->m.vendtax[vnodnum], vnodnum = ~0;
             enodnum < enodnnd; enodnum ++) {     /* Order node vertices with high locality */
          Gnum                  velmnum;          /* Neighboring element                    */
          Gnum                  eelmnum;

          velmnum = meshptr->m.edgetax[enodnum];  /* Get neighboring element */

          if (vexxtax[velmnum].passnum >= passnum) /* If element already scanned */
            continue;                             /* Skip to next element        */

          vexxtax[velmnum].passnum = passnum;     /* Set element as scanned */

          for (eelmnum = meshptr->m.verttax[velmnum]; /* For all neighboring non-halo nodes */
               eelmnum < meshptr->vehdtax[velmnum]; eelmnum ++) {
            Gnum                  vnodend;        /* Neighboring node */

            vnodend = meshptr->m.edgetax[eelmnum]; /* Get neighboring node */

            if (vexxtax[vnodend].passnum <= passnum) { /* If vertex not ordered yet */
              if ((vnodnum == ~0) &&              /* If no next vertex set yet      */
                  (vexxtax[vnodend].vertdist == vnoddist)) /* And in same layer     */
                vnodnum = vnodend;                /* Set neighbor as next vertex    */
              else if (vexxtax[vnodend].passnum < passnum) { /* If not enqueued yet */
                hmeshOrderGpQueuePut (&queue, vnodend); /* Enqueue neighbor vertex  */
                vexxtax[vnodend].passnum = passnum; /* Set it as enqueued           */
              }
            }
          }
        }
      } while (vnodnum != ~0);
    } while (! hmeshOrderGpQueueEmpty (&queue));  /* As long as queue is not empty */
  }

#ifdef SCOTCH_DEBUG_ORDER2
  for (ordeval = ordenum; ordeval < (ordenum + meshptr->vnohnbr); ordeval ++) {
    if (ordeptr->peritab[ordeval] == ~0) {
      errorPrint ("hmeshOrderGp: internal error");
      memFree    (queue.qtab);                    /* Free group leader */
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  memFree (queue.qtab);                           /* Free group leader */

  return (0);
}
