/* Copyright 2004,2007,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**                                 to     30 may 2008     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to     04 aug 2014     **/
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

/* This routine initializes a source graph
** structure.
** It returns:
** - 0  : in all cases.
*/

int
hgraphInit (
Hgraph * restrict const     grafptr)
{
  memSet (grafptr, 0, sizeof (Hgraph));           /* Initialize graph fields     */
  grafptr->s.flagval = GRAPHFREETABS;             /* By default, free all arrays */

  return (0);
}

/* This routine frees a source graph structure.
** It returns:
** - VOID  : in all cases.
*/

void
hgraphExit (
Hgraph * restrict const     grafptr)
{
  hgraphFree (grafptr);
}

/* This routine frees a source graph structure.
** It returns:
** - VOID  : in all cases.
*/

void
hgraphFree (
Hgraph * restrict const     grafptr)
{
  if ((grafptr->vnhdtax != NULL) &&               /* Free end vertex array for non-halo vertices */
      ((grafptr->s.flagval & HGRAPHFREEVNHD) != 0))
    memFree (grafptr->vnhdtax + grafptr->s.baseval);

  graphFree (&grafptr->s);                        /* Free graph data */

#ifdef SCOTCH_DEBUG_HGRAPH2
  memSet (grafptr, ~0, sizeof (Hgraph));          /* Purge graph fields */
#endif /* SCOTCH_DEBUG_HGRAPH2 */
}

/* This routine creates a non-halo graph from a
** halo graph.
** It returns:
** - VOID  : in all cases.
*/

void
hgraphUnhalo (
const Hgraph * restrict const grafptr,
Graph * restrict const        ugrfptr)
{
  ugrfptr->flagval = grafptr->s.flagval & (GRAPHBITSUSED & ~GRAPHFREETABS); /* Remove extended graph class flags and do not allow freeing */
  ugrfptr->baseval = grafptr->s.baseval;
  ugrfptr->vertnbr = grafptr->vnohnbr;
  ugrfptr->vertnnd = grafptr->vnohnnd;
  ugrfptr->verttax = grafptr->s.verttax;
  ugrfptr->vendtax = grafptr->vnhdtax;
  ugrfptr->velotax = grafptr->s.velotax;
  ugrfptr->velosum = grafptr->vnlosum;
  ugrfptr->vnumtax = grafptr->s.vnumtax;
  ugrfptr->vlbltax = NULL;
  ugrfptr->edgenbr = grafptr->enohnbr;
  ugrfptr->edgetax = grafptr->s.edgetax;
  ugrfptr->edlotax = grafptr->s.edlotax;
  ugrfptr->edlosum = grafptr->enohsum;
  ugrfptr->degrmax = grafptr->s.degrmax;          /* Upper bound */
}
