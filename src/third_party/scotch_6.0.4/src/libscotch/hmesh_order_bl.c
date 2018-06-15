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
/**   NAME       : hmesh_order_bl.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module resizes block data using    **/
/**                the block splitting post-processing     **/
/**                algorithm.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 28 sep 2002     **/
/**                                 to     09 feb 2005     **/
/**                # Version 5.0  : from : 25 jul 2007     **/
/**                                 to   : 25 jul 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_ORDER_BL

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "mesh.h"
#include "hmesh.h"
#include "hmesh_order_bl.h"
#include "hmesh_order_st.h"

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
hmeshOrderBl (
const Hmesh * restrict const              meshptr,
Order * restrict const                    ordeptr,
const Gnum                                ordenum, /*+ Zero-based ordering number +*/
OrderCblk * restrict const                cblkptr, /*+ Single column-block        +*/
const HmeshOrderBlParam * restrict const  paraptr)
{
  Gnum                cblknbr;                    /* Number of old column blocks before splitting */
  Gnum                cblknum;                    /* Number of current column block               */

  if (paraptr->cblkmin <= 0) {
    errorPrint ("hmeshOrderBl: invalid minimum block size");
    return     (1);
  }

  if (hmeshOrderSt (meshptr, ordeptr, ordenum, cblkptr, paraptr->strat) != 0) /* Perform ordering strategy */
    return (1);

  if (cblkptr->cblktab == NULL) {                 /* If single column block    */
    if (cblkptr->vnodnbr < (2 * paraptr->cblkmin)) /* If block cannot be split */
	return (0);

    cblknbr = cblkptr->vnodnbr / paraptr->cblkmin; /* Get new number of blocks */

    if ((cblkptr->cblktab = (OrderCblk *) memAlloc (cblknbr * sizeof (OrderCblk))) == NULL) {
      errorPrint ("hgraphOrderBl: out of memory");
      return     (1);
    }
    ordeptr->treenbr += cblknbr;                  /* These more number of tree nodes    */
    ordeptr->cblknbr += cblknbr - 1;              /* These more number of column blocks */
    cblkptr->cblknbr  = cblknbr;

    for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
      cblkptr->cblktab[cblknum].typeval = ORDERCBLKOTHR;
      cblkptr->cblktab[cblknum].vnodnbr = ((cblkptr->vnodnbr + cblknbr - 1) - cblknum) / cblknbr;
      cblkptr->cblktab[cblknum].cblknbr = 0;
      cblkptr->cblktab[cblknum].cblktab = NULL;
    }
  }
  else {                                          /* Block already partitioned */
    for (cblknum = 0; cblknum < cblkptr->cblknbr; cblknum ++) {
      if (hmeshOrderBl (meshptr, ordeptr, ordenum, cblkptr->cblktab + cblknum, paraptr) != 0)
        return (1);
    }
  }

  return (0);
}
