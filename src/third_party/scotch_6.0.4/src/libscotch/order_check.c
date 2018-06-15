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
/**   NAME       : order_check.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module checks the consistency of   **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     20 nov 2003     **/
/**                # Version 5.0  : from : 26 jul 2007     **/
/**                                 to     26 jul 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER_CHECK

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine checks the consistency
** of the given ordering.
** It returns:
** - 0   : if ordering data are consistent.
** - !0  : on error.
*/

static
int
orderCheck2 (
const OrderCblk * restrict const  cblkptr,
Gnum * const                      cblknbr,
Gnum * const                      treenbr)
{
  if (cblkptr->vnodnbr < 1) {
    errorPrint ("orderCheck2: invalid number of vertex nodes (1)");
    return     (1);
  }

  if (cblkptr->cblktab != NULL) {                 /* If node has sons */
    Gnum                vnodnbr;
    Gnum                cblknum;

    if (cblkptr->cblknbr <= 0) {
      errorPrint ("orderCheck2: invalid number of column blocks (1)");
      return     (1);
    }
    *cblknbr += cblkptr->cblknbr - 1;
    *treenbr += cblkptr->cblknbr;
    for (cblknum = vnodnbr = 0; cblknum < cblkptr->cblknbr; cblknum ++) {
      vnodnbr += cblkptr->cblktab[cblknum].vnodnbr;
      if (orderCheck2 (&cblkptr->cblktab[cblknum], cblknbr, treenbr) != 0)
        return (1);
    }
    if (vnodnbr != cblkptr->vnodnbr) {
      errorPrint ("orderCheck2: invalid number of vertex nodes (2)");
      return     (1);
    }
  }
  else if (cblkptr->cblknbr != 0) {
    errorPrint ("orderCheck2: invalid number of column blocks (2)");
    return     (1);
  }

  return (0);
}

int
orderCheck (
const Order * restrict const  ordeptr)
{
  Gnum * restrict     permtab;
  Gnum * restrict     permtax;
  Gnum                treenbr;
  Gnum                cblknbr;
  Gnum                vertnnd;
  Gnum                vertnum;

  if (ordeptr->vnodnbr != ordeptr->cblktre.vnodnbr) {
    errorPrint ("orderCheck: invalid vertex count");
    return     (1);
  }
  if ((ordeptr->cblknbr < 0) || (ordeptr->cblknbr > ordeptr->treenbr)) {
    errorPrint ("orderCheck: invalid column block count (1)");
    return     (1);
  }

  if ((permtab = (Gnum *) memAlloc (ordeptr->vnodnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("orderCheck: out of memory");
    return     (1);
  }
  memSet (permtab, ~0, ordeptr->cblktre.vnodnbr * sizeof (Gnum));
  permtax = permtab - ordeptr->baseval;
  vertnnd = ordeptr->baseval + ordeptr->vnodnbr;

  for (vertnum = 0; vertnum < ordeptr->vnodnbr; vertnum ++) {
    if ((ordeptr->peritab[vertnum] <  ordeptr->baseval) || /* If index not in range */
        (ordeptr->peritab[vertnum] >= vertnnd)) {
      errorPrint ("orderCheck: invalid index");
      memFree    (permtab);
      return     (1);
    } 
    if (permtax[ordeptr->peritab[vertnum]] != ~0) { /* If index already used */
      errorPrint ("orderCheck: duplicate index");
      memFree    (permtab);
      return     (1);
    }
    permtax[ordeptr->peritab[vertnum]] = vertnum; /* Set who updated index */
  }
  for (vertnum = 0; vertnum < ordeptr->vnodnbr; vertnum ++) {
    if (permtab[vertnum] == ~0) {                 /* If index not used */
      errorPrint ("orderCheck: missing index");
      memFree    (permtab);
      return     (1);
    }
  }

  memFree (permtab);

  treenbr =                                       /* Assume there is just a root node */
  cblknbr = 1;
  if (orderCheck2 (&ordeptr->cblktre, &cblknbr, &treenbr) != 0)
    return (1);
  if (cblknbr != ordeptr->cblknbr) {
    errorPrint ("orderCheck: invalid number of column blocks");
    return     (1);
  }
  if (treenbr != ordeptr->treenbr) {
    errorPrint ("orderCheck: invalid number of tree nodes");
    return     (1);
  }

  return  (0);
}
