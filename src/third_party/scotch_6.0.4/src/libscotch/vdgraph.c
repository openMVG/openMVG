/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                separator handling routines.            **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 07 feb 2006     **/
/**                                 to     13 mar 2006     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "vdgraph.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* separator graphs.                 */
/*                                   */
/*************************************/

/* This routine initializes a distributed
** separator graph structure. As for the Dgraph
** structure, in order to avoid collective
** communication whenever possible, the allocation
** of send and receive index arrays is not performed
** in the routine itself, but rather delegated to
** subsequent routines such as dgraphBuild.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
vdgraphInit (
Vdgraph * restrict const    grafptr,              /* Distributed separator graph structure          */
MPI_Comm                    proccomm)             /* Communicator to be used for all communications */
{
  memSet (grafptr, 0, sizeof (Vdgraph));          /* Clear public and private graph fields */

  grafptr->s.proccomm = proccomm;                 /* Set private fields      */
  MPI_Comm_size (proccomm, &grafptr->s.procglbnbr); /* Get communicator data */
  MPI_Comm_rank (proccomm, &grafptr->s.proclocnum);

  return (0);
}

/* This routine frees the contents
** of the given distributed active graph.
** It returns:
** - VOID  : in all cases.
*/

void
vdgraphExit (
Vdgraph * const             grafptr)
{
  if (grafptr->partgsttax != NULL)
    memFree (grafptr->partgsttax + grafptr->s.baseval);
  if (grafptr->fronloctab != NULL)
    memFree (grafptr->fronloctab);

  dgraphExit (&grafptr->s);                       /* Free distributed source graph and its private data (flagval may be corrupted afterwards) */

#ifdef SCOTCH_DEBUG_VDGRAPH2
  memSet (grafptr, ~0, sizeof (Vdgraph));
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
}

/* This routine moves all of the graph
** vertices to the first part.
** It returns:
** - VOID  : in all cases.
*/

void
vdgraphZero (
Vdgraph * const             grafptr)
{
  memSet (grafptr->partgsttax + grafptr->s.baseval, 0, grafptr->s.vertgstnbr * sizeof (GraphPart)); /* Set all local and ghost vertices to part 0 */

  grafptr->compglbloaddlt = grafptr->s.veloglbsum;
  grafptr->compglbload[0] = grafptr->s.veloglbsum; /* No frontier vertices */
  grafptr->compglbload[1] =
  grafptr->compglbload[2] = 0;
  grafptr->compglbsize[0] = grafptr->s.vertglbnbr;
  grafptr->compglbsize[1] =
  grafptr->compglbsize[2] = 0;
  grafptr->complocload[0] = grafptr->s.velolocsum;
  grafptr->complocload[1] =
  grafptr->complocload[2] = 0;
  grafptr->complocsize[0] = grafptr->s.vertlocnbr;
  grafptr->complocsize[1] =
  grafptr->complocsize[2] = 0;
}
