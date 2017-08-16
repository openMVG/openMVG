/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph_separate_sq.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a separator of the **/
/**                given distributed separator graph by    **/
/**                moving all (interesting) vertices of    **/
/**                the given graph to every processor,     **/
/**                running a sequential vertex separation  **/
/**                computing, and projecting back the      **/
/**                best result obtained.                   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 15 feb 2006     **/
/**                                 to     30 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_SEPARATE_SQ

#include "module.h"
#include "common.h"
#include "comm.h"
#include "parser.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "dgraph.h"
#include "vdgraph.h"
#include "vdgraph_separate_sq.h"

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
vdgraphSeparateSqOpBest (
const Gnum * const          in,                   /* First operand                              */
Gnum * const                inout,                /* Second and output operand                  */
const int * const           len,                  /* Number of instances; should be 1, not used */
const MPI_Datatype * const  typedat)              /* MPI datatype; not used                     */
{
  if (inout[3] == 1) {                            /* Handle cases when at least one of them is erroneous */
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
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
vdgraphSeparateSq (
Vdgraph * const                       dgrfptr,    /*+ Distributed graph +*/
const VdgraphSeparateSqParam * const  paraptr)    /*+ Method parameters +*/
{
  Vgraph            cgrfdat;                      /* Centralized vertex separator graph structure      */
  Gnum              reduloctab[7];                /* Local array for best separator data (7 for Bcast) */
  Gnum              reduglbtab[4];                /* Global array for best separator data              */
  MPI_Datatype      besttypedat;                  /* Data type for finding best separator              */
  MPI_Op            bestoperdat;                  /* Handle of MPI operator for finding best separator */
  int               bestprocnum;                  /* Rank of process holding best partition            */
  Gnum * restrict   vnumloctax;
  Gnum              vertlocnum;
  Gnum              complocsize1;
  Gnum              complocload1;
  Gnum              complocload2;
  Gnum              fronlocnbr;
  int               o;

  if ((MPI_Type_contiguous (4, GNUM_MPI, &besttypedat)                                != MPI_SUCCESS) ||
      (MPI_Type_commit (&besttypedat)                                                 != MPI_SUCCESS) ||
      (MPI_Op_create ((MPI_User_function *) vdgraphSeparateSqOpBest, 1, &bestoperdat) != MPI_SUCCESS)) {
    errorPrint ("vdgraphSeparateSq: communication error (1)");
    return     (1);
  }

  reduloctab[0] =                                 /* In case of error, maximum frontier size */
  reduloctab[1] = GNUMMAX;                        /* And maximum load imbalance              */
  reduloctab[2] = dgrfptr->s.proclocnum;
  reduloctab[3] = 0;                              /* Assume sequential separation went fine */

  vnumloctax = dgrfptr->s.vnumloctax;             /* No need for vertex number array when centralizing graph */
  dgrfptr->s.vnumloctax = NULL;
  o = vdgraphGatherAll (dgrfptr, &cgrfdat);
  dgrfptr->s.vnumloctax = vnumloctax;             /* Restore vertex number array */
  if (o != 0) {
    errorPrint ("vdgraphSeparateSq: cannot build centralized graph");
    return     (1);
  }

  if (vgraphSeparateSt (&cgrfdat, paraptr->strat) != 0) { /* Separate centralized graph */
    errorPrint ("vdgraphSeparateSq: cannot separate centralized graph");
    reduloctab[3] = 1;
  }
  else {                                          /* Fill local array with local separator data */
    reduloctab[0] = ((cgrfdat.fronnbr != 0) || ((cgrfdat.compload[0] != 0) && (cgrfdat.compload[1] != 0)))
                    ? cgrfdat.fronnbr
                    : (cgrfdat.fronnbr + cgrfdat.s.vertnbr); /* Partitions with empty separators unwanted if they are completely unbalanced */
    reduloctab[1] = cgrfdat.comploaddlt;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 1, besttypedat, bestoperdat, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateSq: communication error (2)");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (MPI_Allreduce (&reduglbtab[3], &reduloctab[3], 1, GNUM_MPI, MPI_SUM, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateSq: communication error (3)");
    return     (1);
  }
  if ((reduloctab[3] != 0) && (reduloctab[3] != dgrfptr->s.procglbnbr)) {
    errorPrint ("vdgraphSeparateSq: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  if ((MPI_Op_free   (&bestoperdat) != MPI_SUCCESS) ||
      (MPI_Type_free (&besttypedat) != MPI_SUCCESS)) {
    errorPrint ("vdgraphSeparateSq: communication error (4)");
    return     (1);
  }

  if (reduglbtab[3] != 0) {                       /* If none of the sequential methods succeeded */
    vgraphExit (&cgrfdat);
    return     (1);
  }

  bestprocnum = (int) reduglbtab[2];
  if (dgrfptr->s.proclocnum == bestprocnum) {     /* If process holds best partition */
    reduloctab[0] = cgrfdat.compload[0];          /* Global values to share          */
    reduloctab[1] = cgrfdat.compload[1];
    reduloctab[2] = cgrfdat.compload[2];
    reduloctab[3] = cgrfdat.comploaddlt;
    reduloctab[4] = cgrfdat.compsize[0];
    reduloctab[5] = cgrfdat.compsize[1];
    reduloctab[6] = cgrfdat.fronnbr;
  }
  if (MPI_Bcast (reduloctab, 7, GNUM_MPI, bestprocnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateSq: communication error (5)");
    return     (1);
  }
  dgrfptr->compglbload[0] = reduloctab[0];
  dgrfptr->compglbload[1] = reduloctab[1];
  dgrfptr->compglbload[2] = reduloctab[2];
  dgrfptr->compglbloaddlt = reduloctab[3];
  dgrfptr->compglbsize[0] = reduloctab[4];
  dgrfptr->compglbsize[1] = reduloctab[5];
  dgrfptr->compglbsize[2] = reduloctab[6];

  if (commScatterv (cgrfdat.parttax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GRAPHPART_MPI, /* No base for sending as procdsptab holds based values */
                    dgrfptr->partgsttax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GRAPHPART_MPI,
                    bestprocnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateSq: communication error (6)");
    return     (1);
  }

  complocsize1 = 
  complocload1 = 
  complocload2 = 0;
  for (vertlocnum = dgrfptr->s.baseval, fronlocnbr = 0;
       vertlocnum < dgrfptr->s.vertlocnnd; vertlocnum ++) {
    int               partval;
    Gnum              partval1;

    partval  = dgrfptr->partgsttax[vertlocnum];
    partval1 = partval & 1;
    complocsize1 += partval1;                     /* Superscalar update   */
    if (partval == 2)                             /* Build local frontier */
      dgrfptr->fronloctab[fronlocnbr ++] = vertlocnum;
    if (dgrfptr->s.veloloctax != NULL) {
      Gnum              partval2;
      Gnum              veloval;

      veloval       = dgrfptr->s.veloloctax[vertlocnum];
      partval2      = (partval >> 1) & 1;
      complocload1 += (-partval1) & veloval;      /* Superscalar update */
      complocload2 += (-partval2) & veloval;      /* Superscalar update */
    }
  }
  dgrfptr->complocsize[0] = dgrfptr->s.vertlocnbr - fronlocnbr - complocsize1;
  dgrfptr->complocsize[1] = complocsize1;
  dgrfptr->complocsize[2] = fronlocnbr;
  if (dgrfptr->s.veloloctax != NULL) {
    dgrfptr->complocload[0] = dgrfptr->s.velolocsum - complocload1 - complocload2;
    dgrfptr->complocload[1] = complocload1;
    dgrfptr->complocload[2] = complocload2;
  }
  else {
    dgrfptr->complocload[0] = dgrfptr->complocsize[0];
    dgrfptr->complocload[1] = dgrfptr->complocsize[1];
    dgrfptr->complocload[2] = fronlocnbr;
  }

  vgraphExit (&cgrfdat);

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (dgrfptr) != 0) {
    errorPrint ("vdgraphSeparateSq: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  return (0);
}
