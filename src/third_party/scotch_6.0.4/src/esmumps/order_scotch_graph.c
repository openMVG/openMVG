/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : order_scotch_graph.c                    **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This is the interface module with the   **/
/**                libSCOTCH matrix ordering library.      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 20 aug 1998     **/
/**                                 to     18 may 1999     **/
/**                # Version 1.0  : from : 18 mar 2003     **/
/**                                 to     21 jan 2004     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     04 jan 2005     **/
/**                # Version 2.1  : from : 21 jun 2007     **/
/**                                 to     21 jun 2007     **/
/**                # Version 5.0  : from : 08 feb 2008     **/
/**                                 to     01 jun 2008     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     02 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER_GRAPH

#include "common.h"
#ifdef SCOTCH_PTSCOTCH
#include "ptscotch.h"
#else /* SCOTCH_PTSCOTCH */
#include "scotch.h"
#endif /* SCOTCH_PTSCOTCH */
#include "graph.h"
#include "order.h"

/****************************/
/*                          */
/* Graph ordering routines. */
/*                          */
/****************************/

/*+ This routine orders the given
*** graph using the Emilio default
*** ordering strategy.
*** It returns:
*** - 0   : if ordering succeeded.
*** - !0  : on error.
+*/

int
orderGraph (
Order * restrict const      ordeptr,              /*+ Ordering to compute   +*/
Graph * restrict const      grafptr)              /*+ Graph matrix to order +*/
{
  INT                 vertnbr;

  SCOTCH_graphSize (grafptr, &vertnbr, NULL);
  return (orderGraphList (ordeptr, grafptr, vertnbr, NULL));
}

/*+ This routine orders the subgraph of
*** the given graph induced by the given
*** vertex list, using the Emilio default
*** ordering strategy.
*** It returns:
*** - 0   : if ordering succeeded.
*** - !0  : on error.
+*/

int
orderGraphList (
Order * restrict const      ordeptr,              /*+ Ordering to compute        +*/
Graph * restrict const      grafptr,              /*+ Graph matrix to order      +*/
const INT                   listnbr,              /*+ Number of vertices in list +*/
const INT * restrict const  listtab)              /*+ Vertex list array          +*/
{
  return (orderGraphListStrat (ordeptr, grafptr, listnbr, listtab,
                               "c{rat=0.7,cpr=n{sep=/(vert>120)?m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=(|h{pass=10})f{bal=0.2}}}|m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=(|h{pass=10})f{bal=0.2}}};,ole=f{cmin=0,cmax=100000,frat=0.0},ose=g},unc=n{sep=/(vert>120)?m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=(|h{pass=10})f{bal=0.2}}}|m{type=h,rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=0.2},org=(|h{pass=10})f{bal=0.2}}};,ole=f{cmin=15,cmax=100000,frat=0.0},ose=g}}"));
}

/*+ This routine orders the given
*** graph using the given ordering
*** strategy.
*** It returns:
*** - 0   : if ordering succeeded.
*** - !0  : on error.
+*/

int
orderGraphStrat (
Order * restrict const      ordeptr,              /*+ Ordering to compute   +*/
Graph * restrict const      grafptr,              /*+ Graph matrix to order +*/
const char * restrict const stratptr)             /*+ Ordering strategy     +*/
{
  INT                 vertnbr;

  SCOTCH_graphSize (grafptr, &vertnbr, NULL);
  return (orderGraphListStrat (ordeptr, grafptr, vertnbr, NULL, stratptr));
}

/*+ This routine orders the subgraph of
*** the given graph induced by the given
*** vertex list, using the given ordering
*** strategy.
*** It returns:
*** - 0   : if ordering succeeded.
*** - !0  : on error.
+*/

int
orderGraphListStrat (
Order * restrict const      ordeptr,              /*+ Ordering to compute        +*/
Graph * restrict const      grafptr,              /*+ Graph matrix to order      +*/
const INT                   listnbr,              /*+ Number of vertices in list +*/
const INT * restrict const  listtab,              /*+ Vertex list array          +*/
const char * restrict const stratptr)             /*+ Ordering strategy          +*/
{
  SCOTCH_Strat        scotstrat;                  /* Scotch ordering strategy */
  INT                 baseval;
  INT                 vertnbr;
  INT                 edgenbr;
  int                 o;

  if (sizeof (INT) != sizeof (SCOTCH_Num)) {      /* Check integer consistency */
    errorPrint ("orderGraphListStrat: inconsistent integer types");
    return     (1);
  }

  SCOTCH_graphData (grafptr, &baseval, &vertnbr, NULL, NULL, NULL, NULL, &edgenbr, NULL, NULL);

  if (((ordeptr->permtab = (INT *) memAlloc ( vertnbr      * sizeof (INT))) == NULL) ||
      ((ordeptr->peritab = (INT *) memAlloc ( vertnbr      * sizeof (INT))) == NULL) ||
      ((ordeptr->rangtab = (INT *) memAlloc ((vertnbr + 1) * sizeof (INT))) == NULL)) {
    errorPrint ("orderGraphListStrat: out of memory");
    orderExit  (ordeptr);
    orderInit  (ordeptr);
    return     (1);
  }

  SCOTCH_stratInit (&scotstrat);                  /* Initialize default ordering strategy */

  o = SCOTCH_stratGraphOrder (&scotstrat, stratptr);
  if (o == 0)
    o = SCOTCH_graphOrderList (grafptr,           /* Compute graph ordering */
                               (SCOTCH_Num) listnbr, (SCOTCH_Num *) listtab, &scotstrat,
                               (SCOTCH_Num *) ordeptr->permtab,  (SCOTCH_Num *) ordeptr->peritab,
                               (SCOTCH_Num *) &ordeptr->cblknbr, (SCOTCH_Num *) ordeptr->rangtab, NULL);

  SCOTCH_stratExit (&scotstrat);

  if (o != 0) {                                   /* If something failed in Scotch */
    orderExit (ordeptr);                          /* Free ordering arrays          */
    orderInit (ordeptr);
    return    (1);
  }
#ifdef ORDER_DEBUG
  if ((ordeptr->rangtab[0]                != baseval)           ||
      (ordeptr->rangtab[ordeptr->cblknbr] != baseval + vertnbr) ||
      (orderCheck (ordeptr) != 0)) {
    errorPrint ("orderGraphListStrat: invalid ordering");
  }
#endif /* ORDER_DEBUG */

  ordeptr->rangtab = (INT *) memRealloc (ordeptr->rangtab, (ordeptr->cblknbr + 1) * sizeof (INT));

  return (0);
}
