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
/**   NAME       : vgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the separator      **/
/**                handling routines.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 aug 1996     **/
/**                                 to     03 nov 1997     **/
/**                # Version 4.0  : from : 12 dec 2001     **/
/**                                 to     08 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* separator graphs.     */
/*                       */
/*************************/

/* This routine frees the contents
** of the given active graph.
** It returns:
** - VOID  : in all cases.
*/

void
vgraphExit (
Vgraph * const              grafptr)
{
  if (grafptr->parttax != NULL)
    memFree (grafptr->parttax + grafptr->s.baseval);
  if (grafptr->frontab != NULL)
    memFree (grafptr->frontab);

  graphFree (&grafptr->s);                        /* Free source graph */

#ifdef SCOTCH_DEBUG_VGRAPH2
  memSet (grafptr, ~0, sizeof (Vgraph));
#endif /* SCOTCH_DEBUG_VGRAPH2 */
}

/* This routine moves all of the graph
** vertices to the first part.
** It returns:
** - VOID  : in all cases.
*/

void
vgraphZero (
Vgraph * const              grafptr)
{
  memSet (grafptr->parttax + grafptr->s.baseval, 0, grafptr->s.vertnbr * sizeof (GraphPart)); /* Set all vertices to part 0 */

  grafptr->compload[0] = grafptr->s.velosum;      /* No frontier vertices */
  grafptr->compload[1] =
  grafptr->compload[2] = 0;
  grafptr->comploaddlt = grafptr->s.velosum;
  grafptr->compsize[0] = grafptr->s.vertnbr;
  grafptr->compsize[1] =
  grafptr->fronnbr     = 0;
}
