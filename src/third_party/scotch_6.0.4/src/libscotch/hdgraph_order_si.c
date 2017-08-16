/* Copyright 2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_order_si.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders halo distributed     **/
/**                graph vertices using a simple method.   **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 15 apr 2006     **/
/**                                 to     25 jul 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HDGRAPH_ORDER_SI

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dorder.h"
#include "hdgraph.h"
#include "hdgraph_order_si.h"

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
hdgraphOrderSi (
const Hdgraph * restrict const  grafptr,
DorderCblk * restrict const     cblkptr)          /*+ Single column-block +*/
{
  Gnum * restrict     periloctab;
  Gnum * restrict     periloctax;
  Gnum                vnohlocnbr;
  Gnum                vertlocnum;
#ifdef SCOTCH_DEBUG_HDGRAPH1
  int                 cheklocval;
  int                 chekglbval;

  cheklocval = 0;
#endif /* SCOTCH_DEBUG_HDGRAPH1 */
  vnohlocnbr = grafptr->s.vertlocnbr;             /* Get number of local non-halo vertices */
  if ((periloctab = (Gnum *) memAlloc (vnohlocnbr * sizeof (Gnum))) == NULL) { /* Allocate local fragment */
    errorPrint ("hdgraphOrderSi: out of memory");
#ifndef SCOTCH_DEBUG_HDGRAPH1
    return (1);
  }
#else /* SCOTCH_DEBUG_HDGRAPH1 */
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) { /* Communication cannot be merged with a useful one */
    errorPrint ("hdgraphOrderSi: communication error (1)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (periloctab != NULL)
      memFree (periloctab);
    return (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH1 */

  cblkptr->typeval = DORDERCBLKLEAF;              /* Fill node as leaf */
  cblkptr->data.leaf.ordelocval = cblkptr->ordeglbval + grafptr->s.procdsptab[grafptr->s.proclocnum] - grafptr->s.baseval;;
  cblkptr->data.leaf.vnodlocnbr = vnohlocnbr;
  cblkptr->data.leaf.periloctab = periloctab;
  cblkptr->data.leaf.nodelocnbr = 0;
  cblkptr->data.leaf.nodeloctab = NULL;
#ifdef SCOTCH_DEBUG_HDGRAPH2
  cblkptr->data.leaf.cblklocnum = -1;
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  periloctax = periloctab - grafptr->s.baseval;
  if (grafptr->s.vnumloctax == NULL) {            /* If graph is original graph */
    Gnum                vertglbadj;

    vertglbadj = grafptr->s.procdsptab[grafptr->s.proclocnum] - grafptr->s.baseval; /* Set adjustement for global ordering */
    for (vertlocnum = grafptr->s.baseval; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++)
      periloctax[vertlocnum] = vertlocnum + vertglbadj;
  }
  else {                                          /* Graph is not original graph */
    for (vertlocnum = grafptr->s.baseval; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++)
      periloctax[vertlocnum] = grafptr->s.vnumloctax[vertlocnum];
  }

  return (0);
}
