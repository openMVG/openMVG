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
/**   NAME       : bgraph_bipart_ex.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tries to balance the        **/
/**                subgraphs of the partition as best as   **/
/**                it can.                                 **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 25 oct 1994     **/
/**                                 to     03 nov 1994     **/
/**                # Version 3.0  : from : 18 nov 1995     **/
/**                                 to     18 nov 1995     **/
/**                # Version 3.1  : from : 20 nov 1995     **/
/**                                 to     29 nov 1995     **/
/**                # Version 3.2  : from : 15 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     01 jun 2001     **/
/**                # Version 4.0  : from : 11 dec 2003     **/
/**                                 to     11 dec 2003     **/
/**                # Version 5.1  : from : 30 nov 2007     **/
/**                                 to     30 nov 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes
*/

#define BGRAPH_BIPART_EX

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_ex.h"
#include "bgraph_bipart_fm.h"
#include "bgraph_bipart_gg.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bgraphBipartEx (
Bgraph * restrict const     grafptr)
{
  BgraphBipartFmParam parafmdat;                  /* Parameter area for the Fiduccia-Mattheyses algorithm */

  if (grafptr->compload0dlt == 0)                 /* Return if nothing to do */
    return (0);

  parafmdat.movenbr = grafptr->s.vertnbr;
  parafmdat.passnbr = ~0;
  parafmdat.deltval = 0.0L;                       /* Exact balance required */
  if (bgraphBipartFm (grafptr, &parafmdat) != 0)  /* Return if error        */
    return (1);

  if ((grafptr->s.vertnbr > 1) &&                 /* If graph has several vertices but is completely imbalanced */
      ((grafptr->compload0 == 0) || (grafptr->compload0 == grafptr->s.velosum))) {
    BgraphBipartGgParam paraggdat;                /* Parameter area for the Greedy Graph Growing algorithm */

    paraggdat.passnbr = 4;
    if (bgraphBipartGg (grafptr, &paraggdat) != 0) /* Return if error */
      return (1);
  }

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (grafptr) != 0) {
    errorPrint ("bgraphBipartEx: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
