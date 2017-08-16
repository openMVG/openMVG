/* Copyright 2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_kp.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a block ordering   **/
/**                from a k-way edge partition.            **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 17 oct 2012     **/
/**                                 to   : 17 oct 2012     **/
/**                # Version 6.0  : from : 23 aug 2014     **/
/**                                 to   : 23 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_KP

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_kp.h"
#include "hgraph_order_si.h"
#include "kgraph.h"
#include "kgraph_map_st.h"
#include "scotch.h"

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
hgraphOrderKp (
const Hgraph * restrict const             grafptr,
Order * restrict const                    ordeptr,
const Gnum                                ordenum, /*+ Zero-based ordering number +*/
OrderCblk * restrict const                cblkptr, /*+ Single column-block        +*/
const HgraphOrderKpParam * restrict const paraptr)
{
  Kgraph              actgrafdat;
  Gnum * restrict     ordetab;
  Gnum                ordeadj;
  Anum * restrict     parttax;
  Gnum                partnbr;
  Gnum                partnum;
  Gnum * restrict     peritab;
  Gnum                vertnnd;
  Gnum                vertnum;
  Gnum                cblknbr;

  if ((paraptr->partsiz < 1) ||                   /* If nothing to do, order consecutively */
      ((partnbr = grafptr->vnohnbr / paraptr->partsiz) <= 1))
    return (hgraphOrderSi (grafptr, ordeptr, ordenum, cblkptr));

  if ((cblkptr->cblktab = (OrderCblk *) memAlloc (partnbr * sizeof (OrderCblk))) == NULL) { /* Allocate first as it will remain */
    errorPrint ("hgraphOrderKp: out of memory (1)");
    return     (1);
  }

  hgraphUnhalo (grafptr, &actgrafdat.s);          /* Extract non-halo part of given graph       */
  actgrafdat.s.vnumtax = NULL;                    /* Do not keep numbers from nested dissection */

  SCOTCH_archCmplt ((SCOTCH_Arch *) &actgrafdat.a, (SCOTCH_Num) partnbr); /* Build complete graph architecture */

  if ((kgraphInit (&actgrafdat, &actgrafdat.s, &actgrafdat.a, NULL, 0, NULL, NULL, 1, 1, NULL) != 0) ||
      (kgraphMapSt (&actgrafdat, paraptr->strat) != 0)) {
    errorPrint ("hgraphOrderKp: cannot compute partition");
    memFree    (cblkptr->cblktab);
    cblkptr->cblktab = NULL;
    return (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &ordetab, (size_t) (partnbr          * sizeof (Gnum)),
                     &parttax, (size_t) (grafptr->vnohnbr * sizeof (Anum)), NULL) == NULL) {
    errorPrint ("hgraphOrderKp: out of memory (2)");
    memFree    (cblkptr->cblktab);
    cblkptr->cblktab = NULL;
    return (1);
  }
  parttax -= actgrafdat.s.baseval;

  mapTerm (&actgrafdat.m, parttax);               /* Get result of partitioning as terminal part array */

  memSet (ordetab, 0, partnbr * sizeof (Gnum));   /* Reset part count array */
  for (vertnum = actgrafdat.s.baseval, vertnnd = actgrafdat.s.vertnnd;
       vertnum < vertnnd; vertnum ++) {
    ordetab[parttax[vertnum]] ++;                 /* Count number of vertices in each part */
  }

  for (partnum = 0, cblknbr = 0, ordeadj = ordenum; /* For all potential column blocks */
       partnum < partnbr; partnum ++) {
    Gnum                ordetmp;

    ordetmp = ordetab[partnum];
    ordetab[partnum] = ordeadj;
    ordeadj += ordetmp;
    if (ordetmp != 0) {                           /* If part is not empty, one more column block */
      cblkptr->cblktab[cblknbr].typeval = ORDERCBLKOTHR;
      cblkptr->cblktab[cblknbr].vnodnbr = ordetmp;
      cblkptr->cblktab[cblknbr].cblknbr = 0;
      cblkptr->cblktab[cblknbr].cblktab = NULL;
      cblknbr ++;
    }
  }

  ordeptr->treenbr += cblknbr;                    /* These more number of tree nodes    */
  ordeptr->cblknbr += cblknbr - 1;                /* These more number of column blocks */
  cblkptr->cblknbr  = cblknbr;

  peritab = ordeptr->peritab;
  if (grafptr->s.vnumtax == NULL) {               /* If graph is original graph */
    for (vertnum = actgrafdat.s.baseval; vertnum < vertnnd; vertnum ++)
      peritab[ordetab[parttax[vertnum]] ++] = vertnum;
  }
  else {                                          /* Graph is not original graph */
    const Gnum * restrict vnumtax;

    vnumtax = grafptr->s.vnumtax;
    for (vertnum = actgrafdat.s.baseval; vertnum < vertnnd; vertnum ++)
      peritab[ordetab[parttax[vertnum]] ++] = vnumtax[vertnum];
  }

  memFree    (ordetab);                           /* Free group leader */
  kgraphExit (&actgrafdat);

  return (0);
}
