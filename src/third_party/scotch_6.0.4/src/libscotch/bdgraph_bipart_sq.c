/* Copyright 2007,2008,2010,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_sq.c                     **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of a **/
/**                given bipartitioned distributed graph   **/
/**                by moving all (interesting) vertices of **/
/**                the given graph to every processor,     **/
/**                executing a sequential bipartition      **/
/**                routine, and projecting back the best   **/
/**                result obtained.                        **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 27 dec 2007     **/
/**                                 to     14 apr 2011     **/
/**                # Version 6.0  : from : 27 dec 2007     **/
/**                                 to     31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_BIPART_SQ

#include "module.h"
#include "common.h"
#include "comm.h"
#include "arch.h"
#include "parser.h"
#include "graph.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "dgraph.h"
#include "bdgraph.h"
#include "bdgraph_bipart_sq.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine is the reduction-loc operator which
** returns in inout[2] the rank of the process which
** holds the best partition.
** It returns:
** - void  : in all cases.
*/

static
void
bdgraphBipartSqOpBest (
const Gnum * const          in,                   /* First operand                              */
Gnum * const                inout,                /* Second and output operand                  */
const int * const           len,                  /* Number of instances; should be 1, not used */
const MPI_Datatype * const  typedat)              /* MPI datatype; not used                     */
{
  inout[5] |= in[5];                              /* Propagate errors */

  inout[4] += in[4];                              /* Count cases for which a bipartitioning error occured */
  if (inout[3] == 1) {                            /* Handle cases when at least one of them is erroneous  */
    if (in[3] == 1)
      return;

    inout[0] = in[0];
    inout[1] = in[1];
    inout[2] = in[2];
    inout[3] = in[3];
    return;
  }
  else if (in[3] == 1)
    return;

  if ((in[0] < inout[0]) ||                       /* Select best partition */
      ((in[0] == inout[0]) &&
       ((in[1] < inout[1]) ||
        ((in[1] == inout[1]) && (in[2] < inout[2]))))) {
    inout[0] = in[0];
    inout[1] = in[1];
    inout[2] = in[2];
  }
}

/* This routine computes a partition of the
** given distributed graph by gathering as many
** copies of the graph as there are processes
** sharing the distributed graph, running a
** sequential algorithm on them, and collecting
** the best solution found.
** It returns:
** - 0   : if the bipartition could be computed.
** - !0  : on error.
*/

int
bdgraphBipartSq (
Bdgraph * const                     dgrfptr,      /*+ Distributed graph +*/
const BdgraphBipartSqParam * const  paraptr)      /*+ Method parameters +*/
{
  Bgraph            cgrfdat;                      /* Centralized bipartitioned graph structure           */
  Gnum              reduloctab[6];                /* Local array for best bipartition data (7 for Bcast) */
  Gnum              reduglbtab[6];                /* Global array for best bipartition data              */
  MPI_Datatype      besttypedat;                  /* Data type for finding best bipartition              */
  MPI_Op            bestoperdat;                  /* Handle of MPI operator for finding best bipartition */
  int               bestprocnum;                  /* Rank of process holding best partition              */
  Gnum * restrict   vnumloctax;
  Gnum              vertlocnum;
  Gnum              complocsize1;
  Gnum              complocload1;
  Gnum              fronlocnbr;
  int               o;

  if ((MPI_Type_contiguous (6, GNUM_MPI, &besttypedat)                              != MPI_SUCCESS) ||
      (MPI_Type_commit (&besttypedat)                                               != MPI_SUCCESS) ||
      (MPI_Op_create ((MPI_User_function *) bdgraphBipartSqOpBest, 1, &bestoperdat) != MPI_SUCCESS)) {
    errorPrint ("bdgraphBipartSq: communication error (1)");
    return     (1);
  }

  reduloctab[0] =                                 /* In case of error, maximum communication load */
  reduloctab[1] = GNUMMAX;                        /* And maximum load imbalance                   */
  reduloctab[2] = dgrfptr->s.proclocnum;
  reduloctab[3] =                                 /* Assume sequential bipartioning went fine */
  reduloctab[4] = 0;
  reduloctab[5] = 0;                              /* Assume no errors */

  vnumloctax = dgrfptr->s.vnumloctax;             /* No need for vertex number array when centralizing graph */
  dgrfptr->s.vnumloctax = NULL;
  o = bdgraphGatherAll (dgrfptr, &cgrfdat);
  dgrfptr->s.vnumloctax = vnumloctax;             /* Restore vertex number array */
  if (o != 0) {
    errorPrint ("bdgraphBipartSq: cannot build centralized graph");
    return     (1);
  }

  if (bgraphBipartSt (&cgrfdat, paraptr->strat) != 0) { /* Bipartition centralized graph */
    errorPrint ("bdgraphBipartSq: cannot bipartition centralized graph");
    reduloctab[3] =
    reduloctab[4] = 1;
  }
  else {                                          /* Fill local array with local bipartition data */
    reduloctab[0] = ((cgrfdat.fronnbr != 0) || ((cgrfdat.compsize0 != 0) && ((cgrfdat.s.vertnbr - cgrfdat.compsize0) != 0)))
                    ? cgrfdat.commload
                    : GNUMMAX; /* Partitions with empty bipartitions unwanted if they are completely unbalanced */
    reduloctab[1] = cgrfdat.compload0dlt;
  }

  if (dgrfptr->partgsttax == NULL) {
    if (dgraphGhst (&dgrfptr->s) != 0) {          /* Compute ghost edge array if not already present, before copying graph fields */
      errorPrint ("bdgraphBipartSq: cannot compute ghost edge array");
      reduloctab[5] = 1;
    }
    else {
      if ((dgrfptr->partgsttax = (GraphPart *) memAlloc (dgrfptr->s.vertgstnbr * sizeof (GraphPart))) == NULL) {
        errorPrint ("bdgraphBipartSq: out of memory (1)");
        reduloctab[5] = 1;                        /* Allocated data will be freed along with graph structure */
      }
      dgrfptr->partgsttax -= dgrfptr->s.baseval;
    }
    if ((dgrfptr->fronloctab = (Gnum *) memAlloc (dgrfptr->s.vertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("bdgraphBipartSq: out of memory (2)");
      reduloctab[5] = 1;
    }
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 1, besttypedat, bestoperdat, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartSq: communication error (2)");
    return     (1);
  }
  if ((reduloctab[4] != 0) && (reduloctab[4] != dgrfptr->s.procglbnbr)) {
    errorPrint ("bdgraphBipartSq: internal error");
    return     (1);
  }

  if ((MPI_Op_free   (&bestoperdat) != MPI_SUCCESS) ||
      (MPI_Type_free (&besttypedat) != MPI_SUCCESS)) {
    errorPrint ("bdgraphBipartSq: communication error (3)");
    return     (1);
  }

  if (reduglbtab[3] != 0) {                       /* If none of the sequential methods succeeded */
    bgraphExit (&cgrfdat);
    return     (1);
  }

  bestprocnum = (int) reduglbtab[2];
  if (dgrfptr->s.proclocnum == bestprocnum) {     /* If process holds best partition */
    reduloctab[0] = cgrfdat.compload0;            /* Global values to share          */
    reduloctab[1] = cgrfdat.compsize0;
    reduloctab[2] = cgrfdat.commload;
    reduloctab[3] = cgrfdat.commgainextn;
    reduloctab[4] = cgrfdat.fronnbr;
  }
  if (MPI_Bcast (reduloctab, 5, GNUM_MPI, bestprocnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartSq: communication error (4)");
    return     (1);
  }
  dgrfptr->compglbload0    = reduloctab[0];
  dgrfptr->compglbload0dlt = reduloctab[0] - dgrfptr->compglbload0avg;
  dgrfptr->compglbsize0    = reduloctab[1];
  dgrfptr->commglbload     = reduloctab[2];
  dgrfptr->commglbgainextn = reduloctab[3];
  dgrfptr->fronglbnbr      = reduloctab[4];

  if (commScatterv (cgrfdat.parttax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GRAPHPART_MPI, /* No base for sending as procdsptab holds based values */
                    dgrfptr->partgsttax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GRAPHPART_MPI,
                    bestprocnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartSq: communication error (5)");
    return     (1);
  }

  if (dgraphHaloSync (&dgrfptr->s, (byte *) (dgrfptr->partgsttax + dgrfptr->s.baseval), GRAPHPART_MPI) != 0) {
    errorPrint ("bdgraphBipartSq: cannot perform halo exchange");
    return     (1);
  }

  complocsize1 = 
  complocload1 = 0;
  for (vertlocnum = dgrfptr->s.baseval, fronlocnbr = 0;
       vertlocnum < dgrfptr->s.vertlocnnd; vertlocnum ++) {
    int               partval;
    Gnum              partval1;
    Gnum              commcut;
    Gnum              edgelocnum;

    partval  = dgrfptr->partgsttax[vertlocnum];
    partval1 = partval & 1;
    complocsize1 += partval1;                     /* Superscalar update */
    if (dgrfptr->s.veloloctax != NULL) {
      Gnum              veloval;

      veloval       = dgrfptr->s.veloloctax[vertlocnum];
      complocload1 += (-partval1) & veloval;      /* Superscalar update */
    }
    for (edgelocnum = dgrfptr->s.vertloctax[vertlocnum], commcut = 0;
       	 edgelocnum < dgrfptr->s.vendloctax[vertlocnum]; edgelocnum ++) { /* Build local frontier */
      int                 partend;
      int                 partdlt;

      partend  = dgrfptr->partgsttax[dgrfptr->s.edgegsttax[edgelocnum]];
      partdlt  = partval ^ partend;
      commcut |= partdlt;
    }
    if (commcut != 0)
      dgrfptr->fronloctab[fronlocnbr ++] = vertlocnum;
  }
  dgrfptr->fronlocnbr   = fronlocnbr;
  dgrfptr->complocsize0 = dgrfptr->s.vertlocnbr - complocsize1;
  dgrfptr->complocload0 = (dgrfptr->s.veloloctax != NULL) ? (dgrfptr->s.velolocsum - complocload1) : dgrfptr->complocsize0;
  
  bgraphExit (&cgrfdat);

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (dgrfptr) != 0) {
    errorPrint ("bdgraphBipartSq: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}
