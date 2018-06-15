/* Copyright 2007-2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : wgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                Charles-Edmond BICHOT (5.1b)            **/
/**                                                        **/
/**   FUNCTION   : This module contains the data structure **/
/**                handling routines for the vertex overl- **/
/**                apped graph partitioning                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 28 may 2010     **/
/**                                 to   : 29 may 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "wgraph.h"

/***********************************/
/*                                 */
/* Active graph handling routines. */
/*                                 */
/***********************************/

/* This routine initialize the active graph
** corresponding to the source graph
** It returns:
** - VOID  : in all cases.
*/

void
wgraphInit (
Wgraph * restrict const         actgrafptr,     /* Active graph */
const Graph * restrict const    srcgrafptr,     /* Source graph */
const Anum                      partnbr)
{
  actgrafptr->s          = *srcgrafptr;         /* Clone source graph       */
  actgrafptr->s.flagval &= ~GRAPHFREETABS;      /* Do not free its contents */
  actgrafptr->compload   = NULL;                /* No group leader yet      */
  actgrafptr->partnbr    = partnbr;
  actgrafptr->parttax    = NULL;                /* No part array yet */
}

/* This routine frees the contents
** of the given active graph
** It returns:
** - VOID  : in all cases.
*/

void
wgraphExit (
Wgraph * const              grafptr)              /* Active graph */
{
  if (grafptr->compload != NULL)
    memFree (grafptr->compload);                  /* Free group leader */

  graphExit (&grafptr->s);                        /* Free source graph */

#ifdef SCOTCH_DEBUG_WGRAPH2
  memSet (grafptr, ~0, sizeof (Wgraph));
#endif /* SCOTCH_DEBUG_WGRAPH2 */
}

/* This routine builds the active graph
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
wgraphAlloc (
Wgraph * const              grafptr)              /* Active graph */
{
  Gnum                partsiz;
  Gnum * restrict     parttab;

  partsiz = (grafptr->parttax == NULL) ? grafptr->s.vertnbr : 0; /* Allocate part array only if not already */

  if (memAllocGroup ((void **) (void *)
                     &grafptr->compload, (size_t) (grafptr->partnbr   * sizeof (Gnum)),
                     &grafptr->compsize, (size_t) (grafptr->partnbr   * sizeof (Gnum)),
		     &parttab,           (size_t) (partsiz            * sizeof (Anum)),
		     &grafptr->frontab,  (size_t) (grafptr->s.vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("wgraphAlloc: out of memory (1)");
    return     (1);
  }
  if (grafptr->parttax == NULL)                   /* Part array does not need flag as will be group freed */
    grafptr->parttax = parttab - grafptr->s.baseval;

  return (0);
}

/* This routine moves all of the graph
** vertices to the first part.
** It returns:
** - VOID  : in all cases.
*/

void
wgraphZero (
Wgraph * const              grafptr)
{
  memSet (grafptr->compload, 0, grafptr->partnbr * sizeof (Gnum));
  memSet (grafptr->compsize, 0, grafptr->partnbr * sizeof (Gnum));

  grafptr->compload[0] = grafptr->s.velosum;
  grafptr->compsize[0] = grafptr->s.vertnbr;
  grafptr->fronload    = 0;
  grafptr->fronnbr     = 0;

  memSet (grafptr->parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (Anum));  /* Set all vertices to part 0 */
}
