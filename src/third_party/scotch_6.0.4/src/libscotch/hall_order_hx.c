/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hall_order_hx.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains service routines   **/
/**                for interfacing the halo ordering       **/
/**                routines provided by Patrick Amestoy    **/
/**                with the ones of libScotch.             **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 13 jan 2003     **/
/**                                 to   : 28 dec 2004     **/
/**                # Version 5.0  : from : 25 jul 2007     **/
/**                                 to   : 29 may 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HALL_ORDER_HX

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hall_order_hx.h"

/*********************/
/*                   */
/* Service routines. */
/*                   */
/*********************/

/* This routine post-processes the elimination tree
** produced by Amestoy's halo ordering routines.
** On input, the elimination tree is described by
** nvtab and petab.
** On output, we build an un-based elimination tree
** (that is, based to 0), structured as follows:
** - cblknbr : number of column blocks (that is,
**             principal variables).
** - sizetab : sizetab[i] holds the total number of
**   variables (principal and secondary) if i is
**   a principal variable, or ~0 if i is a secondary
**   variable or in halo. (sizetab = lentab)
** - petab : petab[i] holds the father of principal
**   variable i, or ~0 if i is a root.
** - frsttab : when i is a principal variable,
**   frsttab[i] holds the number of the first principal
**   son variable of principal variable i, or ~0 if none.
**   If i is a secondary variable or in halo, frsttab[i]
**   holds ~0 too.
** - nexttab : linked list of principal son variables.
**   When i is a principal variable, j=frsttab[i] is
**   its first son (if any), and k=nexttab[j] is the
**   second son, and l=nexttab[k] is its third son, and
**   so on, or ~0 if none.
** - secntab : linked list of secondary variables.
**   When i is a principal variable, secntab[i] holds
**   the number of the first secondary variable of
**   principal variable i, or ~0 if none. When i is
**   a secondary variable, secntab[i] holds the next
**   secondary variable in linked list, or ~0 if none.
** It returns:
** - 0   : if data structures have been initialized.
** - !0  : on error.
*/

int
hallOrderHxBuild (
const Gnum                  baseval,              /*+ Base value of graph and permutation             +*/
const Gnum                  vertnbr,              /*+ Number of vertices in considered graph          +*/
const Gnum                  vnohnbr,              /*+ Number of non-halo vertices in considered graph +*/
const Gnum * restrict const vnumtax,              /*+ Vertex number array of subgraph, if subgraph    +*/
Order * restrict const      ordeptr,              /*+ Ordering to update                              +*/
OrderCblk * restrict const  cblkptr,              /*+ Multiple column-block of ordering               +*/
Gnum * restrict const       nvartax,
Gnum * restrict const       sizetax,
Gnum * restrict const       fathtax,              /*+ Was petab +*/
Gnum * restrict const       frsttax,
Gnum * restrict const       nexttax,
Gnum * restrict const       secntax,
Gnum * restrict const       desctax,              /*+ Was iwtab                            +*/
Gnum * restrict const       permtax,              /*+ Based direct permutation array       +*/
Gnum * restrict const       peritab,              /*+ Un-based inverse permutation array   +*/
Gnum * restrict const       leaftab,              /*+ Un-based array for storage of leaves +*/
const Gnum                  colmin,
const Gnum                  colmax,
const float                 fillrat)
{
  Gnum                vnohnnd;
  Gnum                cblknbr;
  Gnum                cblknum;
  Gnum                leafnbr;
  Gnum                leafnum;
  Gnum                rootnum;
  Gnum                ordetmp;
  Gnum                i, j, k;

  memSet (desctax + baseval,  0, vertnbr * sizeof (Gnum));
  memSet (sizetax + baseval,  0, vertnbr * sizeof (Gnum));
  memSet (frsttax + baseval, ~0, vertnbr * sizeof (Gnum));
  memSet (secntax + baseval, ~0, vertnbr * sizeof (Gnum));

  vnohnnd = vnohnbr + baseval;

#ifdef SCOTCH_DEBUG_ORDER2
  for (i = baseval; i < vnohnnd; i ++) {
    if ((fathtax[i] > 0) || (fathtax[i] < - vertnbr)) {
      errorPrint ("hallOrderHxBuild: elimination tree out of bounds");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_ORDER2 */
  
  for (i = baseval, cblknbr = 0, rootnum = ~0;    /* Assume no root found yet */
       i < vnohnnd; i ++) {
    if (nvartax[i] != 0) {                        /* If principal variable         */
      cblknbr ++;                                 /* One more column block         */
      sizetax[i] ++;                              /* One more column               */
      if ((fathtax[i] < 0) &&                     /* If not root of tree           */
          (fathtax[i] > - (vnohnbr + 1))) {       /* And father not in halo        */
        fathtax[i]          = baseval - (fathtax[i] + 1); /* Re-base father number */
        nexttax[i]          = frsttax[fathtax[i]]; /* Link vertex to tree          */
        frsttax[fathtax[i]] = i;                  /* Variable is first son         */
        desctax[fathtax[i]] ++;                   /* Father has one more son       */
      }
      else {
        fathtax[i] = ~0;                          /* Father is (pseudo-)root */
        rootnum    = i;                           /* Record (last) root      */
      }
    }
    else {                                        /* If secondary variable */
      fathtax[i] = baseval - (fathtax[i] + 1);    /* Re-base father number */
      if (fathtax[i] >= vnohnnd) {                /* If father in halo     */
        if (frsttax[fathtax[i]] == ~0) {          /* If first such vertex  */
          cblknbr ++;                             /* One more column block */
          sizetax[i] = 1;                         /* One more column       */
          nvartax[i] = 1;                         /* Make it principal     */
          frsttax[fathtax[i]] = i;                /* Record it as root     */
          fathtax[i] = ~0;                        /* Make it (pseudo-)root */
          rootnum    = i;                         /* Record (last) root    */
          continue;                               /* Skip to next vertex   */
        }
        else {
          fathtax[i] = frsttax[fathtax[i]];       /* Get first such vertex as root   */
          nvartax[fathtax[i]] ++;                 /* Record us as secondary variable */
        }
      }
      sizetax[fathtax[i]] ++;                     /* One more column         */
      secntax[i] = secntax[fathtax[i]];           /* Link secondary variable */
      secntax[fathtax[i]] = i;
    }
  }

  for (i = baseval, leafnbr = 0;                  /* Build leaf list for amalgamation */
       i < vnohnnd; i ++) {
    if ((fathtax[i] != ~0) &&                     /* If node has a father        */
        (nvartax[i] !=  0) &&                     /* And is a principal variable */
        (frsttax[i] == ~0))                       /* And is a leaf               */
      leaftab[leafnbr ++] = i;                    /* Add it to leaf list         */
  }

  for (leafnum = 0; leafnum < leafnbr; leafnum ++) { /* As long as candidate leaves exist */
    i = leaftab[leafnum];
    j = fathtax[i];

    if ((sizetax[i] + sizetax[j]) <= colmax) {    /* If will not be too large  */
      if ((sizetax[i] < colmin) ||                /* If column block too small */
          (((float) (2 * sizetax[i]) * (float) (nvartax[j] - nvartax[i] + sizetax[i])) <
           (float) nvartax[j] * (float) nvartax[j] * fillrat)) {
        nvartax[j] += sizetax[i];
        sizetax[j] += sizetax[i];
        nvartax[i]  = 0;
        if (secntax[i] == ~0)                     /* If node had no secondary variables   */
          secntax[i] = secntax[j];                /* Make it take the ones of its father  */
        else if (secntax[j] != ~0) {              /* Else if there is something to append */
          for (k = secntax[i]; secntax[k] != ~0; k = secntax[k]) ; /* Find last node      */
          secntax[k] = secntax[j];                /* Append father list to node list      */
        }
        secntax[j] = i;                           /* Now he is a secondary variable of it */
        if (frsttax[j] == i) {                    /* If node is first son of father       */
          if (frsttax[i] == ~0)                   /* If node has no sons                  */
            frsttax[j] = nexttax[i];              /* First son is now next node           */
          else {
            frsttax[j] = frsttax[i];
            for (k = frsttax[i]; nexttax[k] != ~0; k = nexttax[k])
              fathtax[k] = j;
            fathtax[k] = j;
            nexttax[k] = nexttax[i];
          }
        }
        else {                                    /* Else unlink node from son chain */
          for (k = frsttax[j]; nexttax[k] != i; k = nexttax[k]) ;
          if (frsttax[i] == ~0)                   /* If node has no sons */
            nexttax[k] = nexttax[i];
          else {
            nexttax[k] = frsttax[i];
            for (k = frsttax[i]; nexttax[k] != ~0; k = nexttax[k])
              fathtax[k] = j;
            fathtax[k] = j;
            nexttax[k] = nexttax[i];
          }
        }
        cblknbr --;                               /* One less column block */
      }
    }
    if (((-- desctax[j]) <= 0) &&                 /* If all descendents processed */
        (fathtax[j] != ~0))                       /* And node has a father        */
      leaftab[leafnbr ++] = j;                    /* Enqueue father               */
  }

#ifdef SCOTCH_DEBUG_ORDER2
  memSet (peritab, ~0, vnohnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_ORDER2 */
  ordetmp = hallOrderHxTree (frsttax, nexttax, secntax, peritab, 0, rootnum);
  if (ordetmp < vnohnbr) {                        /* If not all nodes ordered          */
    for (i = baseval; i < rootnum; i ++) {        /* For all potential remaining roots */
      if (fathtax[i] == ~0)                       /* If node is a pseudo-root          */
        ordetmp = hallOrderHxTree (frsttax, nexttax, secntax, peritab, ordetmp, i);
    }
  }
#ifdef SCOTCH_DEBUG_ORDER2
  if (ordetmp != vnohnbr) {
    errorPrint ("hallOrderHxBuild: incomplete elimination tree");
    return     (1);
  }

  memSet (permtax + baseval, ~0, vnohnbr * sizeof (Gnum));

  for (i = 0; i < vnohnbr; i ++) {
    if ((peritab[i] < baseval) || (peritab[i] >= vnohnnd)) {
      errorPrint ("hallOrderHxBuild: permutation out of bounds");
      return     (1);
    } 
    if (permtax[peritab[i]] != ~0) {
      errorPrint ("hallOrderHxBuild: duplicate permutation index");
      return     (1);
    }
    permtax[peritab[i]] = i;
  }
  for (i = baseval; i < vnohnnd; i ++) {
    if (permtax[i] == ~0) {
      errorPrint ("hallOrderHxBuild: unused permutation index");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  if (cblknbr != 1) {                             /* If more than one column block in the end, create subtree */
    if ((cblkptr->cblktab = (OrderCblk *) memAlloc (cblknbr * sizeof (OrderCblk))) == NULL) {
      errorPrint ("hallOrderHxBuild: out of memory");
      return     (1);
    }
    cblkptr->cblknbr  = cblknbr;
    ordeptr->cblknbr += cblknbr - 1;              /* These more column blocks created */
    ordeptr->treenbr += cblknbr;                  /* These more tree nodes created    */

    for (i = 0, cblknum = 0; i < vnohnbr; i ++) {
      if (nvartax[peritab[i]] == 0)               /* If secondary variable      */
        continue;                                 /* Skip to next vertex        */
      cblkptr->cblktab[cblknum].typeval = ORDERCBLKOTHR; /* Build column blocks */
      cblkptr->cblktab[cblknum].vnodnbr = sizetax[peritab[i]];
      cblkptr->cblktab[cblknum].cblknbr = 0;
      cblkptr->cblktab[cblknum].cblktab = NULL;
      cblknum ++;                                 /* One more column block created */
    }
  }

  if (vnumtax != NULL) {                          /* If graph is not original graph */
    for (i = 0; i < vnohnbr; i ++)                /* Finalize inverse permutation   */
      peritab[i] = vnumtax[peritab[i]];
  }

  return (0);
}

/*+ This routine computes the inverse
*** permutation according to the
*** elimination tree.
*** It returns:
*** - >0  : next index to be used to order, in all cases.
+*/

Gnum
hallOrderHxTree (
const Gnum * restrict const frsttax,
const Gnum * restrict const nexttax,
const Gnum * restrict const secntax,
Gnum * restrict const       peritab,
const Gnum                  ordenum,
const Gnum                  nodenum)
{
  Gnum                ordetmp;
  Gnum                nodetmp;

  ordetmp = ordenum;
  for (nodetmp = frsttax[nodenum]; nodetmp != ~0; nodetmp = nexttax[nodetmp])
    ordetmp = hallOrderHxTree (frsttax, nexttax, secntax, peritab, ordetmp, nodetmp);

  peritab[ordetmp ++] = nodenum;                  /* Order principal variable */
  for (nodetmp = secntax[nodenum]; nodetmp != ~0; nodetmp = secntax[nodetmp]) {
    peritab[ordetmp ++] = nodetmp;                /* Order secondary variables */
  }

  return (ordetmp);
}
