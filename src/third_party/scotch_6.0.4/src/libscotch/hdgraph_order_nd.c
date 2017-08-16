/* Copyright 2007,2008,2012 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_order_nd.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders distributed graphs   **/
/**                using the nested dissection algorithm.  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 16 apr 2006     **/
/**                                 to     01 mar 2008     **/
/**                # Version 5.1  : from : 27 sep 2008     **/
/**                                 to     11 nov 2008     **/
/**                # Version 6.0  : from : 12 sep 2012     **/
/**                                 to     12 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HDGRAPH_ORDER_ND

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_st.h"
#include "dgraph.h"
#include "dorder.h"
#include "hdgraph.h"
#include "hdgraph_order_nd.h"
#include "hdgraph_order_sq.h"
#include "hdgraph_order_st.h"
#include "vdgraph.h"
#include "vdgraph_separate_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine builds either a centralized or a
** distributed subgraph, according to the number
** of processes in the given part. The calling
** conventions of this routine have been designed
** so as to allow for multi-threading.
*/

static
void *
hdgraphOrderNdFold2 (
void * const                    dataptr)          /* Pointer to thread data */
{
  HdgraphOrderNdData *            fldthrdptr;     /* Thread input parameters        */
  HdgraphOrderNdGraph * restrict  fldgrafptr;     /* Pointer to folded graph area   */
  Hdgraph                         indgrafdat;     /* Induced distributed halo graph */
  void *                          o;

  fldthrdptr = (HdgraphOrderNdData *) dataptr;
  fldgrafptr = fldthrdptr->fldgrafptr;

  if (hdgraphInduceList (fldthrdptr->orggrafptr, fldthrdptr->indlistnbr, /* Compute unfinished induced subgraph on all processes */
                         fldthrdptr->indlisttab, &indgrafdat) != 0)
    return ((void *) 1);

  o = ((void *) 0);
  if (fldthrdptr->fldprocnbr > 1) {               /* If subpart has several processes, fold a distributed graph  */
    if (hdgraphFold2 (&indgrafdat, fldthrdptr->fldpartval, /* Fold temporary induced subgraph from all processes */
                      &fldgrafptr->data.dgrfdat, fldthrdptr->fldproccomm) != 0)
      o = ((void *) 1);
  }
  else {                                          /* Create a centralized graph */
    Hgraph * restrict         fldcgrfptr;

    fldcgrfptr = (fldthrdptr->fldprocnum == 0) ? &fldgrafptr->data.cgrfdat : NULL; /* See if we are the receiver */
    if (hdgraphGather (&indgrafdat, fldcgrfptr) != 0) /* Gather centralized subgraph from all other processes    */
      o = ((void *) 1);
  }
  hdgraphExit (&indgrafdat);                      /* Free temporary induced graph */

  return (o);
}

static
int
hdgraphOrderNdFold (
Hdgraph * restrict const              orggrafptr,
const Gnum                            indlistnbr0, /* Number of vertices in subgraph 0 */
const Gnum * restrict const           indlisttab0, /* List of vertices in subgraph 0   */
const Gnum                            indlistnbr1, /* Number of vertices in subgraph 1 */
const Gnum * restrict const           indlisttab1, /* List of vertices in subgraph 1   */
HdgraphOrderNdGraph * restrict const  fldgrafptr)
{
  HdgraphOrderNdData        fldthrdtab[2];
  MPI_Comm                  fldproccomm;
  int                       fldprocnbr;
  int                       fldprocnum;
  int                       fldproccol;
  int                       fldpartval;
#ifdef SCOTCH_PTHREAD
  Hdgraph                   orggrafdat;           /* Structure for copying graph fields except communicator */
  pthread_t                 thrdval;              /* Data of second thread                                  */
#endif /* SCOTCH_PTHREAD */
  int                       o;

  if (dgraphGhst (&orggrafptr->s) != 0) {         /* Compute ghost edge array if not already present, before copying graph fields */
    errorPrint ("hdgraphOrderNdFold: cannot compute ghost edge array");
    return     (1);
  }

  fldprocnbr = (orggrafptr->s.procglbnbr + 1) / 2; /* Median cut on number of processors */
  fldthrdtab[0].fldprocnbr = fldprocnbr;
  fldthrdtab[1].fldprocnbr = orggrafptr->s.procglbnbr - fldprocnbr;
  if (orggrafptr->s.proclocnum < fldprocnbr) {    /* Compute color and rank in two subparts */
    fldpartval = 0;
    fldprocnum = orggrafptr->s.proclocnum;
    fldthrdtab[0].fldprocnum  = fldprocnum;
    fldthrdtab[1].fldprocnum  = -1;
    fldthrdtab[1].fldproccomm = MPI_COMM_NULL;
  }
  else {
    fldpartval = 1;
    fldprocnum = orggrafptr->s.proclocnum - fldprocnbr;
    fldprocnbr = orggrafptr->s.procglbnbr - fldprocnbr;
    fldthrdtab[0].fldproccomm = MPI_COMM_NULL;
    fldthrdtab[0].fldprocnum  = -1;
    fldthrdtab[1].fldprocnum  = fldprocnum;
  }
  fldgrafptr->typeval = HDGRAPHORDERNDTYPEDIST;   /* Assume we belong to a distributed subpart */
  fldproccol = fldpartval;                        /* Split color is the part value             */
  if (fldprocnbr <= 1) {                          /* If our part will have only one processor  */
    fldproccol = MPI_UNDEFINED;                   /* Do not create any sub-communicator for it */
    fldgrafptr->typeval = HDGRAPHORDERNDTYPECENT; /* We will host a centralized subgraph       */
  }

  if (MPI_Comm_split (orggrafptr->s.proccomm, fldproccol, fldprocnum, &fldproccomm) != MPI_SUCCESS) {
    errorPrint  ("hdgraphOrderNdFold: communication error");
    return      (1);
  }
  fldthrdtab[fldpartval].fldproccomm = fldproccomm; /* Assign folded communicator to proper part */

  fldthrdtab[0].orggrafptr  = orggrafptr;         /* Load data to pass to the subgraph building routines */
  fldthrdtab[0].indlistnbr  = indlistnbr0;
  fldthrdtab[0].indlisttab  = indlisttab0;
  fldthrdtab[0].fldgrafptr  = fldgrafptr;
  fldthrdtab[0].fldpartval  = 0;
  fldthrdtab[1].indlistnbr  = indlistnbr1;
  fldthrdtab[1].indlisttab  = indlisttab1;
  fldthrdtab[1].fldgrafptr  = fldgrafptr;
  fldthrdtab[1].fldpartval  = 1;

#ifdef SCOTCH_PTHREAD
  orggrafdat = *orggrafptr;                       /* Create a separate graph structure to change its communicator */
  fldthrdtab[1].orggrafptr = &orggrafdat;
  MPI_Comm_dup (orggrafptr->s.proccomm, &orggrafdat.s.proccomm); /* Duplicate communicator to avoid interferences in communications */

  if (pthread_create (&thrdval, NULL, hdgraphOrderNdFold2, (void *) &fldthrdtab[1]) != 0) /* If could not create thread */
    o = ((int) (intptr_t) hdgraphOrderNdFold2 ((void *) &fldthrdtab[0])) || /* Perform inductions in sequence           */
        ((int) (intptr_t) hdgraphOrderNdFold2 ((void *) &fldthrdtab[1]));
  else {                                          /* Newly created thread is processing subgraph 1, so let's process subgraph 0 */
    void *                    o2;

    o = (int) (intptr_t) hdgraphOrderNdFold2 ((void *) &fldthrdtab[0]); /* Work on copy with private communicator */

    pthread_join (thrdval, &o2);
    o |= (int) (intptr_t) o2;
  }
  MPI_Comm_free (&orggrafdat.s.proccomm);
#else /* SCOTCH_PTHREAD */
  fldthrdtab[1].orggrafptr = orggrafptr;

  o = ((int) (intptr_t) hdgraphOrderNdFold2 ((void *) &fldthrdtab[0])) || /* Perform inductions in sequence */
      ((int) (intptr_t) hdgraphOrderNdFold2 ((void *) &fldthrdtab[1]));
#endif /* SCOTCH_PTHREAD */

  return (o);
}

/* This routine performs the ordering.
** It returns:
** - 0   : if the ordering could be computed.
** - !0  : on error.
*/

int
hdgraphOrderNd (
Hdgraph * restrict const                    grafptr,
DorderCblk * restrict const                 cblkptr,
const HdgraphOrderNdParam * restrict const  paraptr)
{
  Vdgraph                   vspgrafdat;           /* Vertex separation graph data       */
  Gnum                      vspvertlocnum;        /* Current vertex in separation graph */
  Gnum * restrict           vspvnumtab[2];        /* Lists of separated parts           */
  Gnum * restrict           vspvnumptr0;
  Gnum * restrict           vspvnumptr1;
  Gnum                      ordeglbval;
  Gnum                      vnodglbnbr;
  Gnum                      cblkfthnum;
  HdgraphOrderNdGraph       indgrafdat01;         /* Induced folded graph area          */
  DorderCblk *              cblkptr01;
  int                       partmax;              /* Induced part having most vertices  */
  int                       procnbr0;             /* Number of processes in first part  */
  int                       cheklocval;
  int                       chekglbval;
  int                       o;

#ifdef SCOTCH_DEBUG_HDGRAPH2
  if (cblkptr->vnodglbnbr != grafptr->s.vertglbnbr) {
    errorPrint ("hdgraphOrderNd: inconsistent parameters");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  if (grafptr->s.procglbnbr == 1) {               /* If we are running on a single process */
    HdgraphOrderSqParam   paradat;

    paradat.ordstratseq = paraptr->ordstratseq;
    return (hdgraphOrderSq (grafptr, cblkptr, &paradat)); /* Run sequentially */
  }

  if (dgraphGhst (&grafptr->s) != 0) {            /* Compute ghost edge array if not already present, to have vertgstnbr (and procsidtab) */
    errorPrint ("hdgraphOrderNd: cannot compute ghost edge array");
    return     (1);
  }

  vspgrafdat.s            = grafptr->s;           /* Get non-halo part of halo distributed graph  */
  vspgrafdat.s.flagval   &= ~DGRAPHFREEALL;       /* Do not free contents of separation graph     */
  vspgrafdat.s.vlblloctax = NULL;                 /* Never mind about vertex labels in the future */
  cheklocval = 0;
  if ((vspgrafdat.fronloctab = (Gnum *) memAlloc (vspgrafdat.s.vertlocnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("hdgraphOrderNd: out of memory (1)");
    vspgrafdat.partgsttax = NULL;
    cheklocval            = 1;
  }
  else if ((vspgrafdat.partgsttax = (GraphPart *) memAlloc (vspgrafdat.s.vertgstnbr * sizeof (GraphPart))) == NULL) {
    errorPrint ("hdgraphOrderNd: out of memory (2)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_HDGRAPH1                      /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphOrderNd: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_HDGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_HDGRAPH1 */
  if (chekglbval != 0) {
    if (vspgrafdat.fronloctab != NULL) {
      if (vspgrafdat.partgsttax != NULL)
        memFree (vspgrafdat.partgsttax);
      memFree (vspgrafdat.fronloctab);
    }
    return (1);
  }

  vspgrafdat.partgsttax -= vspgrafdat.s.baseval;
  vspgrafdat.levlnum     = grafptr->levlnum;      /* Set level of separation graph as level of halo graph */
  vdgraphZero (&vspgrafdat);                      /* Set all local vertices to part 0                     */

  if (vdgraphSeparateSt (&vspgrafdat, paraptr->sepstrat) != 0) { /* Separate vertex-separation graph */
    memFree (vspgrafdat.partgsttax + vspgrafdat.s.baseval);
    memFree (vspgrafdat.fronloctab);
    return  (1);
  }

  if ((vspgrafdat.compglbsize[0] == 0) ||         /* If could not separate more */
      (vspgrafdat.compglbsize[1] == 0)) {
    memFree (vspgrafdat.partgsttax + vspgrafdat.s.baseval); /* Free useless space */
    memFree (vspgrafdat.fronloctab);

    hdgraphOrderSt (grafptr, cblkptr, paraptr->ordstratlea); /* Order this leaf */
    return (0);                                   /* Leaf has been processed    */
  }

  vspvnumtab[0] = vspgrafdat.fronloctab + vspgrafdat.complocsize[2]; /* Build vertex lists within frontier array */
  vspvnumtab[1] = vspvnumtab[0] + vspgrafdat.complocsize[0];
  vspvnumptr0   = vspvnumtab[0];
  vspvnumptr1   = vspvnumtab[1];
  for (vspvertlocnum = vspgrafdat.s.baseval; vspvertlocnum < vspgrafdat.s.vertlocnnd; vspvertlocnum ++) { /* Fill lists */
    GraphPart                 partval;

    partval = vspgrafdat.partgsttax[vspvertlocnum];
    if (partval == 0)
      *vspvnumptr0 ++ = vspvertlocnum;
    else if (partval == 1)
      *vspvnumptr1 ++ = vspvertlocnum;
  }
  memFree (vspgrafdat.partgsttax + vspgrafdat.s.baseval); /* Free useless space */
#ifdef SCOTCH_DEBUG_HDGRAPH2
  vspgrafdat.partgsttax = NULL;                   /* Will cause bug if re-read */

  if ((vspvnumptr0 != vspvnumtab[0] + vspgrafdat.complocsize[0]) ||
      (vspvnumptr1 != vspvnumtab[1] + vspgrafdat.complocsize[1])) {
    errorPrint  ("hdgraphOrderNd: internal error (1)");
    memFree     (vspgrafdat.fronloctab);          /* Free useless space */
    return      (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  cblkptr->typeval = DORDERCBLKNEDI;              /* Node becomes a nested dissection node */

  o = 0;
  if (vspgrafdat.compglbsize[2] != 0) {           /* If separator not empty */
    DorderCblk *        cblkptr2;
    Hdgraph             indgrafdat2;

    cblkptr2 = dorderNew (cblkptr, grafptr->s.proccomm); /* Create separator node */
    cblkptr2->ordeglbval = cblkptr->ordeglbval + grafptr->s.vertglbnbr - vspgrafdat.compglbsize[2];
    cblkptr2->vnodglbnbr = vspgrafdat.compglbsize[2];
    cblkptr2->cblkfthnum = 2;

    cblkptr->data.nedi.cblkglbnbr = 3;            /* It is a three-cell node */

    dgraphInit (&indgrafdat2.s, grafptr->s.proccomm); /* Re-use original graph communicator                                                   */
    if (dgraphInduceList (&grafptr->s, vspgrafdat.complocsize[2], /* Perform non-halo induction for separator, as it will get highest numbers */
                          vspgrafdat.fronloctab, &indgrafdat2.s) != 0) {
      errorPrint ("hdgraphOrderNd: cannot build induced subgraph (1)");
      memFree    (vspgrafdat.fronloctab);         /* Free remaining space */
      return     (1);
    }
    indgrafdat2.vhallocnbr = 0;                   /* No halo on graph */
    indgrafdat2.vhndloctax = indgrafdat2.s.vendloctax;
    indgrafdat2.ehallocnbr = 0;
    indgrafdat2.levlnum    = 0;                   /* Separator graph is at level zero not to be suppressed as an intermediate graph */

    o = hdgraphOrderSt (&indgrafdat2, cblkptr2, paraptr->ordstratsep);
    hdgraphExit   (&indgrafdat2);
    dorderDispose (cblkptr2);                     /* Dispose of separator column block (may be kept as leaf) */
    if (o != 0) {
      memFree (vspgrafdat.fronloctab);            /* Free remaining space */
      return  (1);
    }
  }
  else                                            /* Separator is empty         */
    cblkptr->data.nedi.cblkglbnbr = 2;            /* It is a two-cell tree node */

  partmax  = (vspgrafdat.compglbsize[0] >= vspgrafdat.compglbsize[1]) ? 0 : 1; /* Get part of largest subgraph      */
  procnbr0 = (grafptr->s.procglbnbr + 1) / 2;     /* Get number of processes in part 0 (always more than in part 1) */

  if (grafptr->s.proclocnum < procnbr0) {         /* If process will handle part 0 (bigger part if odd number of processes) */
    ordeglbval = cblkptr->ordeglbval;
    vnodglbnbr = vspgrafdat.compglbsize[partmax];
    cblkfthnum = 0;
  }
  else {                                          /* If process will handle part 1 (smaller part if odd number of processes) */
    ordeglbval = cblkptr->ordeglbval + vspgrafdat.compglbsize[partmax];
    vnodglbnbr = vspgrafdat.compglbsize[partmax ^ 1];
    cblkfthnum = 1;
  }

  o = hdgraphOrderNdFold (grafptr, vspgrafdat.complocsize[partmax], vspvnumtab[partmax],
                          vspgrafdat.complocsize[partmax ^ 1], vspvnumtab[partmax ^ 1], &indgrafdat01);

  if (o == 0) {
    switch (indgrafdat01.typeval) {
      case HDGRAPHORDERNDTYPECENT :
        if ((cblkptr01 = dorderNewSequ (cblkptr)) == NULL) {
          o = 1;
          break;
        }
        if (grafptr->levlnum > 0) {              /* If intermediate level nested dissection graph */
          hdgraphExit   (grafptr);               /* Free graph before going to next level         */
          dorderDispose (cblkptr);               /* Dispose of column block node too              */
        }
        cblkptr01->ordeglbval = ordeglbval;
        cblkptr01->vnodglbnbr = vnodglbnbr;
        cblkptr01->cblkfthnum = cblkfthnum;
        o = hdgraphOrderSq2 (&indgrafdat01.data.cgrfdat, cblkptr01, paraptr->ordstratseq);
        hgraphExit (&indgrafdat01.data.cgrfdat);  /* Free centralized graph here as it is last level                              */
        break;                                    /* No need to dispose of final column block as locally created by dorderNewSequ */
#ifdef SCOTCH_DEBUG_HDGRAPH2
      case HDGRAPHORDERNDTYPEDIST :
#else /* SCOTCH_DEBUG_HDGRAPH2 */
      default :
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
        if ((cblkptr01 = dorderNew (cblkptr, indgrafdat01.data.dgrfdat.s.proccomm)) == NULL) {
          o = 1;
          break;
        }
        if (grafptr->levlnum > 0) {              /* If intermediate level nested dissection graph */
          hdgraphExit   (grafptr);               /* Free graph before going to next level         */
          dorderDispose (cblkptr);               /* Dispose of column block node too              */
        }
        cblkptr01->ordeglbval = ordeglbval;
        cblkptr01->vnodglbnbr = vnodglbnbr;
        cblkptr01->cblkfthnum = cblkfthnum;
        o = hdgraphOrderNd (&indgrafdat01.data.dgrfdat, cblkptr01, paraptr);
        break;
#ifdef SCOTCH_DEBUG_HDGRAPH2
      default :
        errorPrint ("hdgraphOrderNd: internal error (2)");
        o = 1;
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
    }
  }

  memFree (vspgrafdat.fronloctab);                /* Free remaining space */

  return (o);
}
