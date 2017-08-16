/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : symbol_fax_graph.c                      **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This is the block symbolic factoriza-   **/
/**                tion routine for graphs.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     29 sep 1998     **/
/**                # Version 0.2  : from : 08 may 2000     **/
/**                                 to     09 may 2000     **/
/**                # Version 1.0  : from : 01 jun 2002     **/
/**                                 to     03 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 2.0  : from : 21 mar 2003     **/
/**                                 to     21 mar 2003     **/
/**                # Version 3.0  : from : 02 mar 2004     **/
/**                                 to     02 mar 2004     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                                                        **/
/**   NOTES      : # symbolFaxGraph() could have called    **/
/**                  symbolFax() in the regular way, as    **/
/**                  do all of the grid-like factorization **/
/**                  routines. However, for efficiency     **/
/**                  reasons, we have decided to inline    **/
/**                  symbolFax(), to avoid a function call **/
/**                  for every arc.                        **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_FAX
#define SYMBOL_FAX_GRAPH

#include "common.h"
#ifdef SCOTCH_PTSCOTCH
#include "ptscotch.h"
#else /* SCOTCH_PTSCOTCH */
#include "scotch.h"
#endif /* SCOTCH_PTSCOTCH */
#include "graph.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "symbol_fax.h"

/***********************************/
/*                                 */
/* Symbolic factorization routine. */
/*                                 */
/***********************************/

/*+ This routine computes the block symbolic
*** factorization of the given matrix graph
*** according to the given vertex ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolFaxGraph (
SymbolMatrix * const        symbptr,              /*+ Symbolic block matrix [based]        +*/
const Graph * const         grafptr,              /*+ Matrix adjacency structure [based]   +*/
const Order * const         ordeptr)              /*+ Matrix ordering                      +*/
{
  INT                   baseval;
  INT                   vertnbr;
  INT *                 verttab;
  const INT * restrict  verttax;
  INT                   edgenbr;
  INT                   edgenum;
  INT *                 edgetab;
  const INT * restrict  edgetax;

  SCOTCH_graphData (grafptr, &baseval, &vertnbr, &verttab, NULL, NULL, NULL, &edgenbr, &edgetab, NULL);
  verttax = verttab - baseval;
  edgetax = edgetab - baseval;

#define SYMBOL_FAX_ITERATOR(ngbdptr, vertnum, vertend) \
                                    for (edgenum = verttax[vertnum];     \
                                         edgenum < verttax[vertnum + 1]; \
                                         edgenum ++) {                   \
                                      vertend = edgetax[edgenum];

#define SYMBOL_FAX_VERTEX_DEGREE(ngbdptr, vertnum) \
                                    (verttax[(vertnum) + 1] - verttax[(vertnum)])

  {
#define SYMBOL_FAX_INCLUDED
#include "symbol_fax.c"
  }
}
