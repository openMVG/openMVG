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
/**   NAME       : hgraph_order_hx.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains service routines   **/
/**                for the hgraphOrderH{d|f} ordering      **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 23 jan 2004     **/
/**                                 to   : 28 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_HX

#include "module.h"
#include "common.h"
#include "graph.h"
#include "hgraph.h"
#include "hgraph_order_hx.h"

/***********************************/
/*                                 */
/* These are the service routines. */
/*                                 */
/***********************************/

/* This routine fills the input arrays for
** the graph ordering routines.
** It returns:
** - void  : in all cases.
*/

void
hgraphOrderHxFill (
const Hgraph * restrict const             grafptr,
Gnum * restrict const                     petab,
Gnum * restrict const                     lentab,
Gnum * restrict const                     iwtab,
Gnum * restrict const                     elentab,
Gnum * restrict const                     pfreptr)
{
  Gnum * restrict             petax;
  Gnum * restrict             iwtax;
  Gnum * restrict             lentax;
  Gnum * restrict             elentax;
  Gnum                        vertadj;            /* Index adjustment for vertices */
  Gnum                        vertnum;
  Gnum                        vertnew;
  Gnum                        edgenew;

  petax   = petab - 1;                            /* Base HAMF arrays at base 1 */
  iwtax   = iwtab - 1;
  lentax  = lentab - 1;
  elentax = elentab - 1;

  vertadj = 1 - grafptr->s.baseval;
  for (vertnum = grafptr->s.baseval, vertnew = edgenew = 1; /* Process non-halo vertices */
       vertnum < grafptr->vnohnnd; vertnum ++, vertnew ++) {
    Gnum                      degrval;
    Gnum                      edgenum;

    degrval = grafptr->s.vendtax[vertnum] - grafptr->s.verttax[vertnum];
    petax[vertnew]   = edgenew;
    lentax[vertnew]  = degrval;
    elentax[vertnew] = degrval;

    for (edgenum = grafptr->s.verttax[vertnum];
         edgenum < grafptr->s.vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = grafptr->s.edgetax[edgenum] + vertadj;
  }
  for ( ; vertnum < grafptr->s.vertnnd; vertnum ++, vertnew ++) { /* Process halo vertices */
    Gnum                      degrval;
    Gnum                      edgenum;

    degrval = grafptr->s.verttax[vertnum] - grafptr->s.vendtax[vertnum];
    petax[vertnew]   = edgenew;
    lentax[vertnew]  = (degrval != 0) ? degrval : (-1 - grafptr->s.vertnbr);
    elentax[vertnew] = 0;

    for (edgenum = grafptr->s.verttax[vertnum];
         edgenum < grafptr->s.vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = grafptr->s.edgetax[edgenum] + vertadj;
  }

  *pfreptr = edgenew;                             /* Set index to first free area */
}
