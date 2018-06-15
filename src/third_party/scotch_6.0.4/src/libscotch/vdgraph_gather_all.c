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
/**   NAME       : vdgraph_gather_all.c                    **/
/**                                                        **/
/**   AUTHORS    : Cedric CHEVALIER                        **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized Vgraph on all      **/
/**                processes by gathering the pieces of a  **/
/**                distributed Vdgraph.                    **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 29 apr 2006     **/
/**                                 to     01 mar 2008     **/
/**                # Version 5.1  : from : 18 apr 2009     **/
/**                                 to     30 jul 2010     **/
/**                                                        **/
/**   NOTES      : # The definitions of MPI_Gather and     **/
/**                  MPI_Gatherv indicate that elements in **/
/**                  the receive array should not be       **/
/**                  written more than once. Great care    **/
/**                  should be taken to enforce this rule, **/
/**                  especially when the number of         **/
/**                  vertices in the centralized graph is  **/
/**                  smaller than the number of            **/
/**                  processes.                            **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define VDGRAPH

#include "module.h"
#include "common.h"
#include "comm.h"
#include "graph.h"
#include "vgraph.h"
#include "dgraph.h"
#include "vdgraph.h"

/* This function gathers on all processes
** the pieces of a distributed Vdgraph to
** build a centralized Vgraph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
vdgraphGatherAll (
const Vdgraph * restrict const dgrfptr,           /* Distributed graph */
Vgraph * restrict              cgrfptr)           /* Centralized graph */
{
  int * restrict     froncnttab;                  /* Count array for gather operations        */
  int * restrict     frondsptab;                  /* Displacement array for gather operations */
  int                fronlocnbr;                  /* Also int to enforce MPI standard         */
  int                cheklocval;
#ifdef SCOTCH_DEBUG_VDGRAPH1
  int                chekglbval;
#endif /* SCOTCH_DEBUG_VDGRAPH1 */
  int                procnum;

  cheklocval = 0;
#ifdef SCOTCH_DEBUG_VDGRAPH1
  if (cgrfptr == NULL)                            /* Centralized graphs should be provided by all */
    cheklocval = 1;
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (1)");
    return     (1);
  }
  if (chekglbval != 0) {
    errorPrint ("vdgraphGatherAll: centralized graphs should be provided on every process");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH1 */

  if (dgraphGatherAll (&dgrfptr->s, &cgrfptr->s) != 0) {
    errorPrint ("vdgraphGatherAll: cannot build centralized graph");
    return     (1);
  }

  cgrfptr->parttax = NULL;                        /* In case of error */
  cgrfptr->frontab = NULL;
  if (((cgrfptr->parttax = (GraphPart *) memAlloc (cgrfptr->s.vertnbr * sizeof (GraphPart))) == NULL) ||
      ((cgrfptr->parttax -= cgrfptr->s.baseval,
        cgrfptr->frontab = (Gnum *) memAlloc (cgrfptr->s.vertnbr * sizeof (Gnum))) == NULL)) {
    errorPrint ("vdgraphGatherAll: out of memory (1)");
#ifndef SCOTCH_DEBUG_VDGRAPH1
    vgraphExit (cgrfptr);
    return     (1);
  }
#else /* SCOTCH_DEBUG_VDGRAPH1 */
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (2)");
    return     (1);
  }
  if (chekglbval != 0) {
    vgraphExit (cgrfptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH1 */

  cgrfptr->levlnum = dgrfptr->levlnum;            /* Set level of separation graph as level of halo graph */

  if (dgrfptr->partgsttax == NULL) {              /* If distributed graph does not have a part array yet */
    vgraphZero (cgrfptr);
    return     (0);
  }

  if (memAllocGroup ((void **) (void *)           /* Allocate tempory arrays to gather separator vertices */
                     &froncnttab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)),
                     &frondsptab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)), NULL) == NULL) {
    errorPrint ("vdgraphGatherAll: out of memory (2)");
#ifndef SCOTCH_DEBUG_VDGRAPH1
    vgraphExit (cgrfptr);
    return     (1);
  }
#else /* SCOTCH_DEBUG_VDGRAPH1 */
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (3)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (froncnttab != NULL)
      memFree (froncnttab);
    vgraphExit (cgrfptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH1 */

  if (commAllgatherv (dgrfptr->partgsttax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GRAPHPART_MPI, /* Get parttax of distributed graph */
                      cgrfptr->parttax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GRAPHPART_MPI, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (4)");
    return     (1);
  }

  fronlocnbr = (int) dgrfptr->complocsize[2];
  if (MPI_Allgather (&fronlocnbr, 1, MPI_INT,     /* Compute how separator vertices are distributed */
                     froncnttab, 1, MPI_INT, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (5)");
    return     (1);
  }
  frondsptab[0] = 0;                              /* Offset 0 for first process                                                    */
  for (procnum = 1; procnum < dgrfptr->s.procglbnbr; procnum ++) /* Adjust index sub-arrays for all processes except the first one */
    frondsptab[procnum] = frondsptab[procnum - 1] + froncnttab[procnum - 1];

  if (MPI_Allgatherv (dgrfptr->fronloctab, fronlocnbr, GNUM_MPI, /* Gather separator vertices */
                      cgrfptr->frontab, froncnttab, frondsptab, GNUM_MPI, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphGatherAll: communication error (6)");
    return     (1);
  }

  for (procnum = 1; procnum < dgrfptr->s.procglbnbr; procnum ++) { /* Adjust index sub-arrays for all processes except the first one */
    Gnum               vertnum;
    Gnum               vertnnd;

    for (vertnum = (Gnum) frondsptab[procnum], vertnnd = vertnum + (Gnum) froncnttab[procnum];
         vertnum < vertnnd; vertnum ++)
      cgrfptr->frontab[vertnum] += (Gnum) dgrfptr->s.procdsptab[procnum] - dgrfptr->s.baseval;
  }

  memFree (froncnttab);                           /* Free group leader */

  for (procnum = 0; procnum < dgrfptr->s.proclocnum; procnum ++) /* Desynchronize random generators across processes */
    cheklocval = intRandVal (2);
  intPerm (cgrfptr->frontab, dgrfptr->compglbsize[2]); /* Compute permutation of frontier array to have different solutions on every process */

  cgrfptr->compload[0] = dgrfptr->compglbload[0]; /* Update other fields */
  cgrfptr->compload[1] = dgrfptr->compglbload[1];
  cgrfptr->compload[2] = dgrfptr->compglbload[2];
  cgrfptr->comploaddlt = dgrfptr->compglbloaddlt;
  cgrfptr->compsize[0] = dgrfptr->compglbsize[0];
  cgrfptr->compsize[1] = dgrfptr->compglbsize[1];
  cgrfptr->fronnbr     = dgrfptr->compglbsize[2];

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vgraphCheck (cgrfptr) != 0) {
    errorPrint ("vdgraphGatherAll: internal error");
    vgraphExit (cgrfptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  return (0);
}
