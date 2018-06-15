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
/**   NAME       : order.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles generic orderings.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 oct 1996     **/
/**                                 to     27 aug 1998     **/
/**                # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     26 dec 2004     **/
/**                # Version 5.0  : from : 25 jul 2007     **/
/**                                 to     25 jul 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine initializes an ordering
** with respect to a given source graph.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
orderInit (
Order * restrict const      ordeptr,
const Gnum                  baseval,
const Gnum                  vnodnbr,
Gnum * restrict const       peritab)
{
  ordeptr->flagval         = ORDERNONE;
  ordeptr->baseval         = baseval;
  ordeptr->vnodnbr         = vnodnbr;
  ordeptr->treenbr         =                      /* Initialize a simple blocking */
  ordeptr->cblknbr         = 1;
  ordeptr->cblktre.typeval = ORDERCBLKOTHR;
  ordeptr->cblktre.vnodnbr = vnodnbr;
  ordeptr->cblktre.cblknbr = 0;
  ordeptr->cblktre.cblktab = NULL;
  ordeptr->peritab         = peritab;

  if (ordeptr->peritab == NULL) {                 /* Inverse permutation must be allocated */
    ordeptr->flagval |= ORDERFREEPERI;            /* Flag it so it will be freed           */
    if ((ordeptr->peritab = (Gnum *) memAlloc (vnodnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("orderInit: out of memory");
      return     (1);
    }
  }

#ifdef SCOTCH_DEBUG_ORDER2
  memSet (ordeptr->peritab, ~0, vnodnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_ORDER2 */

  return (0);
}

/* This routine frees the contents
** of the given ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

void
orderExit (
Order * restrict const      ordeptr)
{
  if (ordeptr->cblktre.cblktab != NULL)           /* Free column block tree */
    orderExit2 (ordeptr->cblktre.cblktab, ordeptr->cblktre.cblknbr);

  if ((ordeptr->peritab != NULL) && ((ordeptr->flagval & ORDERFREEPERI) != 0)) /* If peritab is group leader */
    memFree (ordeptr->peritab);                   /* Free group leader */

#ifdef SCOTCH_DEBUG_ORDER2
  memSet (ordeptr, ~0, sizeof (Order));
#endif /* SCOTCH_DEBUG_ORDER2 */
}

static
void
orderExit2 (
OrderCblk * restrict const  cblktab,
const Gnum                  cblknbr)
{
  Gnum                cblknum;

  for (cblknum = 0; cblknum < cblknbr; cblknum ++) {
    if (cblktab[cblknum].cblktab != NULL)
      orderExit2 (cblktab[cblknum].cblktab, cblktab[cblknum].cblknbr);
  }
  memFree (cblktab);
}

/* This routine computes the inverse permutation
** of the given permutation, according to the
** direct and inverse base values.
** It returns:
** - VOID  : in all cases.
*/

void
orderPeri (
const Gnum * restrict const   permtab,            /* Permutation to invert             */
const Gnum                    permbas,            /* Permutation base value            */
const Gnum                    permnbr,            /* Number of permutation indices     */
Gnum * restrict const         peritab,            /* Array of inverse permutation      */
const Gnum                    peribas)            /* Base value of inverse permutation */
{
  Gnum                permnum;

  for (permnum = 0; permnum < permnbr; permnum ++)
    peritab[permtab[permnum] - permbas] = permnum + peribas;
}

/* This routine computes the column block
** range array of the given ordering.
** It returns:
** - VOID  : in all cases.
*/

void
orderRang (
const Order * restrict const  ordeptr,            /* Ordering                      */
Gnum * restrict const         rangtab)            /* Column block range array [+1] */
{
  Gnum *              rangptr;
  Gnum                ordenum;

  rangptr = rangtab;                              /* Set beginning of range array */
  ordenum = ordeptr->baseval;                     /* Set initial number           */
  orderRang2 (&rangptr, &ordenum, &ordeptr->cblktre);
  *rangptr = ordenum;                             /* Set end of range array */
}

static
void
orderRang2 (
Gnum ** restrict const            rangppt,
Gnum * restrict const             ordeppt,
const OrderCblk * restrict const  cblkptr)
{
  Gnum                cblknum;
#ifdef SCOTCH_DEBUG_ORDER2
  Gnum * restrict     rangtmp;

  if (cblkptr->vnodnbr < 1)
    errorPrint ("orderRang2: internal error (1)");
#endif /* SCOTCH_DEBUG_ORDER2 */

  if (cblkptr->cblktab == NULL) {                 /* If leaf of column block tree  */
    *(*rangppt) ++ = *ordeppt;                    /* Set beginning of column block */
    *ordeppt      += cblkptr->vnodnbr;            /* Advance by column block size  */
  }
  else {
#ifdef SCOTCH_DEBUG_ORDER2
    rangtmp = *rangppt;
#endif /* SCOTCH_DEBUG_ORDER2 */
    for (cblknum = 0; cblknum < cblkptr->cblknbr; cblknum ++)
      orderRang2 (rangppt, ordeppt, &cblkptr->cblktab[cblknum]);
#ifdef SCOTCH_DEBUG_ORDER2
    if ((*ordeppt - *rangtmp) != cblkptr->vnodnbr)
      errorPrint ("orderRang2: internal error (2)");
#endif /* SCOTCH_DEBUG_ORDER2 */
  }
}

/* This routine computes the separator tree
** array of the given ordering.
** It returns:
** - VOID  : in all cases.
*/

void
orderTree (
const Order * restrict const  ordeptr,            /* Ordering                          */
Gnum * restrict const         treetab)            /* Column block separator tree array */
{
  Gnum                cblanum;

  cblanum = ordeptr->cblknbr + ordeptr->baseval - 1; /* Set number of last column block */
  orderTree2 (treetab - ordeptr->baseval, &cblanum, &ordeptr->cblktre, -1);

#ifdef SCOTCH_DEBUG_ORDER2
  if (cblanum != ordeptr->baseval - 1)
    errorPrint ("orderTree: internal error");
#endif /* SCOTCH_DEBUG_ORDER2 */
}

static
void
orderTree2 (
Gnum * restrict const             treetax,        /* Based access to tree table                                          */
Gnum * restrict const             cblaptr,        /* Pointer to current number of last column block, in descending order */
const OrderCblk * restrict const  cblkptr,        /* Current column block tree node                                      */
Gnum                              cbfanum)        /* Current number of ancestor separator column block                   */
{
#ifdef SCOTCH_DEBUG_ORDER2
  if (cblkptr->vnodnbr < 1)
    errorPrint ("orderTree2: internal error (1)");
#endif /* SCOTCH_DEBUG_ORDER2 */

  if (cblkptr->cblktab == NULL)                   /* If leaf of column block tree */
    treetax[(*cblaptr) --] = cbfanum;             /* Set its ancestor             */
  else {                                          /* Node has sub-nodes           */
    Gnum                cblknum;

    cblknum = cblkptr->cblknbr - 1;               /* Assume all column blocks will be scanned */
    if ((cblkptr->cblknbr == 3) &&                /* If node is a nested dissection node      */
        (cblkptr->typeval == ORDERCBLKNEDI)) {    /* With a non-empty separator               */
      Gnum                cblanum;

      cblanum = *cblaptr;                         /* Save number of last column block of separator   */
      orderTree2 (treetax, cblaptr, &cblkptr->cblktab[cblknum], cbfanum); /* Scan separator apart    */
      cbfanum = cblanum;                          /* Separator becomes most recent ancestor of parts */
      cblknum = 1;                                /* Only scan the two parts, not the separator      */
    }
      
    for ( ; cblknum >= 0; cblknum --) {
      orderTree2 (treetax, cblaptr, &cblkptr->cblktab[cblknum], cbfanum);
#ifdef SCOTCH_DEBUG_ORDER2
      if (*cblaptr < -1)
        errorPrint ("orderTree2: internal error (2)");
#endif /* SCOTCH_DEBUG_ORDER2 */
    }
  }
}
