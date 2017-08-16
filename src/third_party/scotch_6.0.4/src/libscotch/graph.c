/* Copyright 2004,2007,2011,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     31 oct 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     08 jun 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     15 sep 1998     **/
/**                # Version 3.3  : from : 22 sep 1998     **/
/**                                 to     31 dec 1998     **/
/**                # Version 4.0  : from : 24 nov 2001     **/
/**                                 to     22 apr 2004     **/
/**                # Version 5.1  : from : 08 mar 2011     **/
/**                                 to     08 mar 2011     **/
/**                # Version 6.0  : from : 09 sep 2012     **/
/**                                 to     09 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH

#include "module.h"
#include "common.h"
#include "graph.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

/* This routine initializes a source graph
** structure.
** It returns:
** - 0  : in all cases.
*/

int
graphInit (
Graph * const               grafptr)
{
  memSet (grafptr, 0, sizeof (Graph));            /* Initialize graph fields     */
  grafptr->flagval = GRAPHFREETABS;               /* By default, free all arrays */

  return (0);
}

/* This routine frees a source graph structure.
** It returns:
** - VOID  : in all cases.
*/

void
graphExit (
Graph * const               grafptr)
{
  graphFree (grafptr);                            /* Free graph data */

#ifdef SCOTCH_DEBUG_GRAPH2
  memSet (grafptr, ~0, sizeof (Graph));           /* Purge graph fields */
#endif /* SCOTCH_DEBUG_GRAPH2 */
}

/* This routine frees the graph data.
** It returns:
** - VOID  : in all cases.
*/

void
graphFree (
Graph * const               grafptr)
{
  if (((grafptr->flagval & GRAPHFREEEDGE) != 0) && /* If edgetab must be freed */
      (grafptr->edgetax != NULL))                 /* And if it exists          */
    memFree (grafptr->edgetax + grafptr->baseval); /* Free it                  */

  if ((grafptr->flagval & GRAPHFREEVERT) != 0) {  /* If verttab/vendtab must be freed                            */
    if ((grafptr->vendtax != NULL) &&             /* If vendtax is distinct from verttab                         */
        (grafptr->vendtax != grafptr->verttax + 1) && /* (if vertex arrays grouped, vendtab not distinct anyway) */
        ((grafptr->flagval & GRAPHVERTGROUP) == 0))
      memFree (grafptr->vendtax + grafptr->baseval); /* Then free vendtax                                 */
    if (grafptr->verttax != NULL)                 /* Free verttab anyway, as it is the array group leader */
      memFree (grafptr->verttax + grafptr->baseval);
  }
  if ((grafptr->flagval & GRAPHFREEVNUM) != 0) {  /* If vnumtab must be freed         */
    if ((grafptr->vnumtax != NULL) &&             /* And is not in vertex array group */
        ((grafptr->flagval & GRAPHVERTGROUP) == 0))
      memFree (grafptr->vnumtax + grafptr->baseval);
  }
  if ((grafptr->flagval & GRAPHFREEOTHR) != 0) {  /* If other arrays must be freed */
    if ((grafptr->velotax != NULL) &&             /* Free graph tables             */
        ((grafptr->flagval & GRAPHVERTGROUP) == 0))
      memFree (grafptr->velotax + grafptr->baseval);
    if ((grafptr->vlbltax != NULL) &&
        ((grafptr->flagval & GRAPHVERTGROUP) == 0))
      memFree (grafptr->vlbltax + grafptr->baseval);
    if ((grafptr->edlotax != NULL) &&
        ((grafptr->flagval & GRAPHEDGEGROUP) == 0))
      memFree (grafptr->edlotax + grafptr->baseval);
  }

#ifdef SCOTCH_DEBUG_GRAPH2
  memSet (grafptr, ~0, sizeof (Graph));           /* Purge graph fields */
#endif /* SCOTCH_DEBUG_GRAPH2 */
  grafptr->flagval = GRAPHNONE;                   /* Allow to double-call graphFree or call graphExit */
}
