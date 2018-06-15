/* Copyright 2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_gather_all.c                    **/
/**                                                        **/
/**   AUTHORS    : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized Bgraph on all      **/
/**                processors by gathering the pieces of   **/
/**                a distributed Bdgraph.                  **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 21 dec 2007     **/
/**                                 to     14 apr 2011     **/
/**                # Version 6.0  : from : 29 aug 2014     **/
/**                                 to     31 aug 2014     **/
/**                                                        **/
/**   NOTES      : # The definitions of MPI_Gather and     **/
/**                  MPI_Gatherv indicate that elements in **/
/**                  the receive array should not be       **/
/**                  written more than once. Great care    **/
/**                  should be taken to enforce this rule, **/
/**                  especially when the number of         **/
/**                  vertices in the centralized graph is  **/
/**                  smaller than the number of            **/
/**                  processors.                           **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "comm.h"
#include "arch.h"
#include "graph.h"
#include "bgraph.h"
#include "dgraph.h"
#include "bdgraph.h"

/* This function gathers on all processors
** the pieces of a distributed Bdgraph to
** build a centralized Bgraph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
bdgraphGatherAll (
const Bdgraph * restrict const dgrfptr,            /* Distributed graph */
Bgraph * restrict              cgrfptr)            /* Centralized graph */
{
  int * restrict     froncnttab;                   /* Count array for gather operations        */
  int * restrict     fronvrttab;                   /* Displacement array for gather operations */
  int                fronlocnbr;                   /* Also int to enforce MPI standard         */
  int                cheklocval;
#ifdef SCOTCH_DEBUG_BDGRAPH1
  int                chekglbval;
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  int                procnum;

  cheklocval = 0;
#ifdef SCOTCH_DEBUG_BDGRAPH1
  if (cgrfptr == NULL)                            /* Centralized graphs should be provided by all */
    cheklocval = 1;
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphGatherAll: communication error (1)");
    return     (1);
  }
  if (chekglbval != 0) {
    errorPrint ("bdgraphGatherAll: centralized graphs should be provided on every process");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH1 */

  if (dgraphGatherAll (&dgrfptr->s, &cgrfptr->s) != 0) {
    errorPrint ("bdgraphGatherAll: cannot build centralized graph");
    return     (1);
  }

  cgrfptr->s.flagval |= BGRAPHFREEFRON | BGRAPHFREEPART | BGRAPHFREEVEEX;
  cgrfptr->veextax = NULL;                        /* In case of error */
  cgrfptr->parttax = NULL;
  cgrfptr->frontab = NULL;
  if ((cgrfptr->frontab = (Gnum *) memAlloc (cgrfptr->s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("bdgraphGatherAll: out of memory (1)");
#ifndef SCOTCH_DEBUG_BDGRAPH1
    bgraphExit (cgrfptr);
    return     (1);
#else /* SCOTCH_DEBUG_BDGRAPH1 */
    cheklocval = 1;
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  }
  else if ((cgrfptr->parttax = (GraphPart *) memAlloc (cgrfptr->s.vertnbr * sizeof (GraphPart))) == NULL) {
    errorPrint ("bdgraphGatherAll: out of memory (2)");
#ifndef SCOTCH_DEBUG_BDGRAPH1
    bgraphExit (cgrfptr);
    return     (1);
#else /* SCOTCH_DEBUG_BDGRAPH1 */
    cheklocval = 1;
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  }
  else {
    cgrfptr->parttax -= cgrfptr->s.baseval;
  
    if (dgrfptr->veexloctax != NULL) {
      if ((cgrfptr->veextax = (Gnum *) memAlloc (cgrfptr->s.vertnbr * sizeof (Gnum))) == NULL) {
        errorPrint ("bdgraphGatherAll: out of memory (3)");
#ifndef SCOTCH_DEBUG_BDGRAPH1
        bgraphExit (cgrfptr);
        return     (1);
#else /* SCOTCH_DEBUG_BDGRAPH1 */
        cheklocval = 1;
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
      }
      else
        cgrfptr->veextax -= cgrfptr->s.baseval;
    }
  }

#ifdef SCOTCH_DEBUG_BDGRAPH1
  if (cheklocval == 0) {
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
    if (memAllocGroup ((void **) (void *)           /* Allocate tempory arrays to gather frontiers */
                       &froncnttab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)),
                       &fronvrttab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("bdgraphGatherAll: out of memory (4)");
#ifndef SCOTCH_DEBUG_BDGRAPH1
      bgraphExit (cgrfptr);
      return     (1);
    }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
      cheklocval = 1;
    }
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphGatherAll: communication error (2)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (froncnttab != NULL)
      memFree (froncnttab);                       /* Free group leader */
    bgraphExit (cgrfptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH1 */

  cgrfptr->compload0min  = dgrfptr->compglbload0min; /* Set constant fields of the centralized graph as those of the distibuted graph */
  cgrfptr->compload0max  = dgrfptr->compglbload0max;
  cgrfptr->compload0avg  = dgrfptr->compglbload0avg;
  cgrfptr->commloadextn0 = dgrfptr->commglbloadextn0; 
  cgrfptr->commgainextn0 = dgrfptr->commglbgainextn0;
  cgrfptr->domndist      = dgrfptr->domndist; 
  cgrfptr->domnwght[0]   = dgrfptr->domnwght[0]; 
  cgrfptr->domnwght[1]   = dgrfptr->domnwght[1]; 
  cgrfptr->vfixload[0]   =                        /* Fixed vertices will soon be available in PT-Scotch */
  cgrfptr->vfixload[1]   = 0;
  cgrfptr->levlnum       = dgrfptr->levlnum;           

  if (dgrfptr->partgsttax == NULL) {              /* If distributed graph does not have a part array yet */
    bgraphZero (cgrfptr);
    memFree    (froncnttab);                      /* Free group leader */
    return     (0);
  }

  if (commAllgatherv (dgrfptr->partgsttax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GRAPHPART_MPI, /* Get parttax of distributed graph */
                      cgrfptr->parttax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GRAPHPART_MPI, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphGatherAll: communication error (4)");
    return     (1);
  }
  
  if (dgrfptr->veexloctax != NULL) {
    if (commAllgatherv (dgrfptr->veexloctax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GNUM_MPI, /* Get veextax of distributed graph */
                        cgrfptr->veextax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GNUM_MPI, dgrfptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("bdgraphGatherAll: communication error (5)");
      return     (1);
    }
  }

  fronlocnbr = (int) dgrfptr->fronlocnbr;
  if (MPI_Allgather (&fronlocnbr, 1, MPI_INT,     /* Compute how frontiers are distributed */
                     froncnttab, 1, MPI_INT, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphGatherAll: communication error (6)");
    return     (1);
  }
  fronvrttab[0] = 0;                              /* Offset 0 for first process                                                     */
  for (procnum = 1; procnum < dgrfptr->s.procglbnbr; procnum ++) /* Adjust index sub-arrays for all processors except the first one */
    fronvrttab[procnum] = fronvrttab[procnum - 1] + froncnttab[procnum - 1];

  if (MPI_Allgatherv (dgrfptr->fronloctab, (int) dgrfptr->fronlocnbr, GNUM_MPI, /* Gather frontiers */
                      cgrfptr->frontab, froncnttab, fronvrttab, GNUM_MPI, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphGatherAll: communication error (7)");
    return     (1);
  }

  for (procnum = 1; procnum < dgrfptr->s.procglbnbr; procnum ++) { /* Adjust index sub-arrays for all processors except the first one */
    Gnum               vertnum;
    Gnum               vertnnd;

    for (vertnum = (Gnum) fronvrttab[procnum], vertnnd = (Gnum) fronvrttab[procnum] + (Gnum) froncnttab[procnum];
         vertnum < vertnnd; vertnum ++)
      cgrfptr->frontab[vertnum] += (Gnum) dgrfptr->s.procdsptab[procnum] - dgrfptr->s.baseval;
  }

  memFree (froncnttab);                           /* Free group leader */

  for (procnum = 0; procnum < dgrfptr->s.proclocnum; procnum ++) /* Desynchronize random generators across processes */
    cheklocval = intRandVal (2);
  intPerm (cgrfptr->frontab, dgrfptr->fronglbnbr); /* Compute permutation of frontier array to have different solutions on every process */

  cgrfptr->compload0     = dgrfptr->compglbload0; /* Update other fields */
  cgrfptr->compload0dlt  = dgrfptr->compglbload0dlt;
  cgrfptr->compsize0     = dgrfptr->compglbsize0;
  cgrfptr->commload      = dgrfptr->commglbload;
  cgrfptr->commgainextn  = dgrfptr->commglbgainextn;
  cgrfptr->commgainextn0 = dgrfptr->commglbgainextn0; 
  cgrfptr->fronnbr       = dgrfptr->fronglbnbr;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bgraphCheck (cgrfptr) != 0) {
    errorPrint ("bdgraphGatherAll: internal error");
    bgraphExit (cgrfptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}
