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
/**   NAME       : hgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 17 jan 2002     **/
/**                                 to     01 dec 2003     **/
/**                # Version 5.0  : from : 19 dec 2006     **/
/**                                 to     19 dec 2006     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "hgraph.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

/* This routine checks the consistency
** of the given halo graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
hgraphCheck (
const Hgraph * restrict const grafptr)
{
  Gnum                vertnum;                    /* Number of current vertex */
  Gnum                edgenum;                    /* Number of current edge   */
  Gnum                enohsum;

  if (graphCheck (&grafptr->s) != 0) {
    errorPrint ("hgraphCheck: invalid graph structure in halo graph");
    return     (1);
  }

  if ((grafptr->vnohnbr < 0)                                        ||
      (grafptr->vnohnbr > grafptr->s.vertnbr)                       ||
      (grafptr->vnohnnd != (grafptr->vnohnbr + grafptr->s.baseval)) ||
      (grafptr->vnlosum > grafptr->s.velosum)                       ||
      (grafptr->enohnbr > grafptr->s.edgenbr)                       ||
      (grafptr->enohsum < grafptr->enohnbr)) {
    errorPrint ("hgraphCheck: invalid halo graph parameters");
    return     (1);
  }

  enohsum = (grafptr->s.edlotax == NULL) ? grafptr->enohnbr : 0;
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->vnohnnd; vertnum ++) { /* For all non-halo vertices */
    if ((grafptr->vnhdtax[vertnum] < grafptr->s.verttax[vertnum]) ||
        (grafptr->vnhdtax[vertnum] > grafptr->s.vendtax[vertnum])) {
      errorPrint ("hgraphCheck: invalid non-halo end vertex array");
      return     (1);
    }

    if (grafptr->s.edlotax != NULL) {
      Gnum                edgenum;

      for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->vnhdtax[vertnum]; edgenum ++)
        enohsum += grafptr->s.edlotax[edgenum];
    }
  }

  if (grafptr->enohsum != enohsum) {
    errorPrint ("hgraphCheck: invalid non-halo edge load sum");
    return     (1);
  }

  for ( ; vertnum < grafptr->s.vertnnd; vertnum ++) { /* For all halo vertices */
    for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
      if (grafptr->s.edgetax[edgenum] >= grafptr->vnohnnd) { /* If two halo vertices connected together */
        errorPrint ("hgraphCheck: halo vertices should not be connected together");
        return     (1);
      }
    }
  }

  return (0);
}
