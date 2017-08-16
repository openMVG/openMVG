/* Copyright 2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_order_sq.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders distributed graphs   **/
/**                by centralizing them on a single        **/
/**                process and running a sequential        **/
/**                ordering strategy.                      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 nov 2008     **/
/**                                 to     11 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HDGRAPH_ORDER_SQ

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_st.h"
#include "dgraph.h"
#include "dorder.h"
#include "hdgraph.h"
#include "hdgraph_order_sq.h"

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
hdgraphOrderSq (
Hdgraph * restrict const                    grafptr,
DorderCblk * restrict const                 cblkptr,
const HdgraphOrderSqParam * restrict const  paraptr)
{
  Hgraph                    cgrfdat;              /* Centralized halo graph data       */
  Hgraph *                  cgrfptr;              /* Pointer to centralized halo graph */
  int                       o;

  cgrfptr = (grafptr->s.proclocnum == 0) ? &cgrfdat : NULL; /* Set root process */
  if (hdgraphGather (grafptr, cgrfptr) != 0) { /* Gather centralized subgraph   */
    errorPrint ("hdgraphOrderSq: cannot create centralized graph");
    return     (1);
  }

  o = 0;
  if (cgrfptr != NULL) {
    o = hdgraphOrderSq2 (&cgrfdat, cblkptr, paraptr->ordstratseq);
    hgraphFree (&cgrfdat);
  }

  return (o);
}

int
hdgraphOrderSq2 (
Hgraph * restrict const         grafptr,
DorderCblk * restrict const     cblkptr,
const Strat * restrict const    stratptr)
{
  Order                     corddat;              /* Centralized ordering structure */
  Gnum * restrict           vnumtax;
  Gnum * restrict           peritab;
  int                       o;

  if (orderInit (&corddat, grafptr->s.baseval, cblkptr->vnodglbnbr, NULL) != 0) {
    errorPrint ("hdgraphOrderSq2: cannot initialize centralized ordering");
    return     (1);
  }

  vnumtax = grafptr->s.vnumtax;                   /* Save number array of subgraph to order               */
  grafptr->s.vnumtax = NULL;                      /* Assume graph does not have one (fake original graph) */

  if (hgraphOrderSt (grafptr, &corddat, 0, &corddat.cblktre, stratptr) != 0) { /* Compute sequential ordering */
    orderExit (&corddat);
    return    (1);
  }
#ifdef SCOTCH_DEBUG_HDGRAPH2
  if (orderCheck (&corddat) != 0) {
    errorPrint ("hdgraphOrderSq2: invalid centralized ordering");
    orderExit  (&corddat);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  peritab = corddat.peritab;                      /* Get permutation array */

  if (vnumtax != NULL) {
    Gnum                      perinum;

    grafptr->s.vnumtax = vnumtax;                 /* Restore vertex number array            */
    for (perinum = 0; perinum < grafptr->vnohnbr; perinum ++) /* Adjust inverse permutation */
      peritab[perinum] = vnumtax[peritab[perinum]];
  }

  cblkptr->typeval = DORDERCBLKLEAF;              /* Fill node as leaf */
  cblkptr->data.leaf.ordelocval = cblkptr->ordeglbval;
  cblkptr->data.leaf.vnodlocnbr = cblkptr->vnodglbnbr;
  cblkptr->data.leaf.periloctab = peritab;
  cblkptr->data.leaf.nodelocnbr = corddat.treenbr - 1; /* Get number of tree nodes, save for root */
  o = 0;
  if (corddat.treenbr > 1) {
    cblkptr->data.leaf.cblklocnum = dorderNewSequIndex (cblkptr, corddat.treenbr - 1); /* Reserve local indices for local nodes */
    if ((cblkptr->data.leaf.nodeloctab = hdgraphOrderSqTree (&corddat)) == NULL) {
      errorPrint ("hdgraphOrderSq2: cannot import centralized separation tree");
      o = 1;
    }
    if (corddat.cblktre.typeval == ORDERCBLKNEDI) /* If root of centralized tree is a nested dissection node */
      cblkptr->typeval |= DORDERCBLKNEDI;         /* Distributed leaf is also a nested dissection node       */
  }
  else
    cblkptr->data.leaf.nodeloctab = NULL;

  corddat.flagval = ORDERNONE;                    /* Do not free permutation array */
  orderExit (&corddat);                           /* Free permutation tree         */

  return (o);
}

/* This routine builds the distributed part of
** a distributed halo graph. This is a distinct
** routine to allow for multi-threading.
*/

static
DorderNode *
hdgraphOrderSqTree (
const Order * const             cordptr)
{
  DorderNode *        nodetab;
  Gnum                nodenum;
  Gnum                cblknum;

  if ((nodetab = memAlloc ((cordptr->treenbr - 1) * sizeof (DorderNode))) == NULL) { /* "- 1" as root of tree will not be copied */
    errorPrint ("hdgraphOrderSqTree: out of memory");
    return     (NULL);
  }

  nodenum = 0;                                    /* Start labeling nodes from 0 */
  for (cblknum = 0; cblknum < cordptr->cblktre.cblknbr; cblknum ++)
    hdgraphOrderSqTree2 (nodetab, &nodenum, &cordptr->cblktre.cblktab[cblknum], -1, cblknum); /* Root of tree is labeled "-1" */

#ifdef SCOTCH_DEBUG_HDGRAPH2
  if (nodenum != (cordptr->treenbr - 1)) {
    errorPrint ("hdgraphOrderSqTree: internal error");
    return     (NULL);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  return (nodetab);
}

static
void
hdgraphOrderSqTree2 (
DorderNode * const              nodetab,
Gnum * const                    nodeptr,
const OrderCblk * const         cblkptr,
const Gnum                      fathnum,
const Gnum                      fcbknum)
{
  Gnum                nodenum;
  DorderNode *        nodetmp;
  Gnum                cblknum;

  nodenum = (*nodeptr) ++;
  nodetmp = &nodetab[nodenum];
  nodetmp->fathnum = fathnum;
  nodetmp->typeval = (Gnum) cblkptr->typeval;
  nodetmp->vnodnbr = cblkptr->vnodnbr;
  nodetmp->cblknum = fcbknum;

  for (cblknum = 0; cblknum < cblkptr->cblknbr; cblknum ++)
    hdgraphOrderSqTree2 (nodetab, nodeptr, &cblkptr->cblktab[cblknum], nodenum, cblknum);
}
