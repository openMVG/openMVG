/* Copyright 2004,2007,2009,2011,2012,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_match.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source graph   **/
/**                matching functions, generated from the  **/
/**                generic pattern.                        **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 oct 2012     **/
/**                                 to     26 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_MATCH

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "graph_match.h"

/*
**  The static variables.
*/

static void              (* graphmatchfuncseqtab[]) (GraphCoarsenThread *) = { /* Array of sequential matching routines */
                              GRAPHMATCHFUNCBLOCK (Seq)
                            };

#ifdef GRAPHMATCHTHREAD

static void              (* graphmatchfuncthrbegtab[]) (GraphCoarsenThread *) = { /* Array of threaded matching start routines */
                              GRAPHMATCHFUNCBLOCK (ThrBeg)
                            };

static void              (* graphmatchfuncthrmidtab[]) (GraphCoarsenThread *) = { /* Array of threaded matching intermediate routines */
                              GRAPHMATCHFUNCBLOCK (ThrMid)
                            };

static void              (* graphmatchfuncthrendtab[]) (GraphCoarsenThread *) = { /* Array of threaded matching end routines */
                              GRAPHMATCHFUNCBLOCK (ThrEnd)
                            };

#endif /* GRAPHMATCHTHREAD */

/***************************/
/*                         */
/* The sequential matching */
/* subroutines.            */
/*                         */
/***************************/

#define GRAPHMATCHSCANP1INPERT                    /* Perturbation scan for first pass  */
#define GRAPHMATCHSCANP2INPERT                    /* Perturbation scan for second pass */

#define GRAPHMATCHSCANNAME          graphMatchSeqNfNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchSeqNfNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchSeqNfVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchSeqNfVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB

#define GRAPHMATCHSCANPFIXTAB
#define GRAPHMATCHSCANNAME          graphMatchSeqFxNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchSeqFxNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchSeqFxVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchSeqFxVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB
#undef GRAPHMATCHSCANPFIXTAB

#undef GRAPHMATCHSCANP1INPERT
#undef GRAPHMATCHSCANP2INPERT

/*************************/
/*                       */
/* The threaded matching */
/* start subroutines.    */
/*                       */
/*************************/

#ifdef GRAPHMATCHTHREAD

/* Start subroutines
*/

#define GRAPHMATCHSCANP1INPERT                    /* Perturbation scan for first pass  */
#define GRAPHMATCHSCANP2INPERT                    /* Perturbation scan for second pass */
#define GRAPHMATCHSCANP2OUTQUEUE                  /* Queue storage for second pass     */

#define GRAPHMATCHSCANNAME          graphMatchThrBegNfNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrBegNfNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrBegNfVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrBegNfVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB

#define GRAPHMATCHSCANPFIXTAB
#define GRAPHMATCHSCANNAME          graphMatchThrBegFxNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrBegFxNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrBegFxVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrBegFxVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB
#undef GRAPHMATCHSCANPFIXTAB

#undef GRAPHMATCHSCANP1INPERT
#undef GRAPHMATCHSCANP2INPERT
#undef GRAPHMATCHSCANP2OUTQUEUE

/* Intermediate subroutines
*/

#define GRAPHMATCHSCANP2INQUEUE                   /* Read queue for second (only) pass */
#define GRAPHMATCHSCANP2OUTQUEUE                  /* Queue storage for second pass     */

#define GRAPHMATCHSCANNAME          graphMatchThrMidNfNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrMidNfNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrMidNfVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrMidNfVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB

#define GRAPHMATCHSCANPFIXTAB
#define GRAPHMATCHSCANNAME          graphMatchThrMidFxNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrMidFxNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrMidFxVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrMidFxVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB
#undef GRAPHMATCHSCANPFIXTAB

#undef GRAPHMATCHSCANP2INQUEUE
#undef GRAPHMATCHSCANP2OUTQUEUE

/* End subroutines
*/

#define GRAPHMATCHSCANP2INQUEUE                   /* Read queue for second (only) pass */

#define GRAPHMATCHSCANNAME          graphMatchThrEndNfNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrEndNfNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrEndNfVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrEndNfVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB

#define GRAPHMATCHSCANPFIXTAB
#define GRAPHMATCHSCANNAME          graphMatchThrEndFxNvNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANEDLOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrEndFxNvEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME

#define GRAPHMATCHSCANVELOTAB
#define GRAPHMATCHSCANNAME          graphMatchThrEndFxVlEl
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANEDLOTAB

#define GRAPHMATCHSCANNAME          graphMatchThrEndFxVlNe
#include "graph_match_scan.c"
#undef GRAPHMATCHSCANNAME
#undef GRAPHMATCHSCANVELOTAB
#undef GRAPHMATCHSCANPFIXTAB

#undef GRAPHMATCHSCANP2INQUEUE

#endif /* GRAPHMATCHTHREAD */

/*************************/
/*                       */
/* The matching routine. */
/*                       */
/*************************/

/* This routine performs the sequential
** initialization of the global mating
** data structures, so as to indicate
** that no mating will be performed.
** It returns:
** - 0  : in all cases.
*/

void
graphMatchNone (
GraphCoarsenData * restrict coarptr)
{
#ifdef SCOTCH_PTHREAD
  coarptr->finelocktax = NULL;
  coarptr->finequeutab = NULL;
  coarptr->fendptr     = (void (*) (void *)) NULL;
  coarptr->fmidptr     = (void (*) (void *)) NULL;
#endif /* SCOTCH_PTHREAD */
  coarptr->fbegptr     = (void (*) (void *)) NULL;
}

/* This routine performs the sequential
** initialization of the global mating
** data structures, before the threads
** are launched.
** It returns:
** - 0  : if initialization could be performed.
** - 1  : on error.
*/

int
graphMatchInit (
GraphCoarsenData * restrict coarptr)
{
  int                 flagval;

  const Graph * restrict const  finegrafptr = coarptr->finegrafptr;
#if ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC))
  const Gnum                    finevertnbr = finegrafptr->vertnbr;
  const Gnum                    baseval     = finegrafptr->baseval;
  const int                     thrdnbr     = coarptr->thrddat.thrdnbr;
#endif /* ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC)) */

  flagval = (finegrafptr->edlotax != NULL) ? 1 : 0;
  if (finegrafptr->velotax != NULL)
    flagval |= 2;
  if ((coarptr->finevfixnbr > 0) || (coarptr->fineparotax != NULL))
    flagval |= 4;

#if ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC))
  if (thrdnbr > 1) {
    if (memAllocGroup ((void **) (void *)
                       &coarptr->finequeutab, (size_t) (finevertnbr * sizeof (Gnum)),
                       &coarptr->finelocktax, (size_t) (finevertnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("graphMatchInit: out of memory");
      return     (1);
    }
    coarptr->finelocktax -= baseval;

    coarptr->fbegptr = (void (*) (void *)) graphmatchfuncthrbegtab[flagval];
    coarptr->fmidptr = (void (*) (void *)) graphmatchfuncthrmidtab[flagval];
    coarptr->fendptr = (void (*) (void *)) graphmatchfuncthrendtab[flagval];
  }
  else  {
    coarptr->finequeutab = NULL;
    coarptr->finelocktax = NULL;                  /* If deterministic behavior wanted, no threaded mating */
    coarptr->fbegptr = (void (*) (void *)) graphmatchfuncseqtab[flagval];
#ifdef SCOTCH_DEBUG_GRAPH2
    coarptr->fmidptr = (void (*) (void *)) NULL;
    coarptr->fendptr = (void (*) (void *)) NULL;
#endif /* SCOTCH_DEBUG_GRAPH2 */
  }
#else /* ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC)) */
  coarptr->fbegptr = (void (*) (void *)) graphmatchfuncseqtab[flagval];
#endif /* ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC)) */

  return (0);
}

/* This routine merges the results of two mating
** threads and re-launches a mating operations
** if necessary.
*/

#ifdef GRAPHMATCHTHREAD
static
void
graphMatchReduce (
GraphCoarsenThread * restrict const tlocptr,      /* Pointer to local thread */
void * restrict const               vlocptr,      /* Pointer to local value  */
void * restrict const               vremptr)      /* Pointer to remote value */
{
  GraphCoarsenData * restrict const   coarptr = (GraphCoarsenData *) (tlocptr->thrddat.grouptr);
  GraphCoarsenThread * restrict const tremptr = (GraphCoarsenThread *) vremptr;
  const int                           thrdnbr = coarptr->thrddat.thrdnbr;
  const int                           thrdnum = tlocptr->thrddat.thrdnum;
  Gnum                                qremnbr;

  qremnbr = tremptr->finequeunnd - tremptr->finequeubas; /* Number of enqueued fine vertices in second thread */

  memMov (coarptr->finequeutab + tlocptr->finequeunnd, /* Merge queues */
          coarptr->finequeutab + tremptr->finequeubas,
          qremnbr * sizeof (Gnum));
  tlocptr->finequeunnd += qremnbr;
  tlocptr->coarvertnbr += tremptr->coarvertnbr;

  if ((thrdnum == 0) && (((tremptr - tlocptr) << 1) >= thrdnbr)) /* If last join */
    coarptr->fendptr (tlocptr);                   /* Call end match routine      */
  else
    coarptr->fmidptr (tlocptr);                   /* Call intermediate match routine */
}
#endif /* GRAPHMATCHTHREAD */

/* This routine matches the vertices of the given
** graph, according to various constraints. The
** matching can be either single-threaded or
** multi-threaded.
** It returns:
** - 0  : if matching could be performed.
** - 1  : on error.
*/

void
graphMatch (
GraphCoarsenThread * restrict thrdptr)            /*+ Pointer to incomplete match data array +*/
{
  Gnum                finevertsiz;
#ifdef SCOTCH_DEBUG_GRAPH2
  Gnum                finevertnum;
#endif /* SCOTCH_DEBUG_GRAPH2 */

  GraphCoarsenData * restrict const coarptr     = (GraphCoarsenData *) (thrdptr->thrddat.grouptr);
  const Graph * restrict const      finegrafptr = coarptr->finegrafptr;
  const Gnum                        finevertbas = thrdptr->finevertbas; /* Get fine vertex range */
  const Gnum                        finevertnnd = thrdptr->finevertnnd;
  Gnum * restrict const             finematetax = coarptr->finematetax;
  const Gnum                        baseval     = finegrafptr->baseval;

  if (coarptr->fbegptr == NULL)                   /* If user-provided mating, nothing to do */
    return;

  thrdptr->finequeubas = finevertbas;             /* Assume matching range is fine vertex processing range */
  thrdptr->finequeunnd = finevertnnd;
  thrdptr->coarvertnbr = 0;                       /* No coarse vertices created yet */

  finevertsiz = finevertnnd - finevertbas;        /* Compute fine vertex range */
  memSet (finematetax + finevertbas, ~0, finevertsiz * sizeof (Gnum));
#if ((defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC))
  if (coarptr->thrddat.thrdnbr > 1) {
    memSet (coarptr->finelocktax + finevertbas, 0, finevertsiz * sizeof (int)); /* Initialize local part of lock array for concurrent accesses */
    threadBarrier (thrdptr);                      /* finematetax and finelocktax must have been globally initialized before we can go on       */

    coarptr->fbegptr (thrdptr);                   /* Perform bulk on local part                                   */
    threadReduce (thrdptr, thrdptr, (ThreadReduceFunc) graphMatchReduce, 0); /* Reduce work on remaining vertices */

    if (thrdptr->thrddat.thrdnum == 0) {
      coarptr->coarvertnbr = thrdptr->coarvertnbr;  /* Global number of coarse vertices is reduced number */
      memFree (coarptr->finequeutab);             /* Free group leader of matching data                   */
    }

    threadBarrier (thrdptr);                      /* coarptr->coarvertnbr must be known to all */
  }
  else
#else /* (defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC) */
#ifdef GRAPHCOARSENTHREAD                         /* If matching was called from a threaded environment */
  if (coarptr->thrddat.thrdnbr > 1) {
    threadBarrier (thrdptr);                      /* finematetax must have been fully initialized before we go on */

    thrdptr->finequeubas = finegrafptr->baseval;  /* Thread 0 will handle all of fine graph vertices */
    thrdptr->finequeunnd = finegrafptr->vertnnd;
    if (thrdptr->thrddat.thrdnum == 0) {          /* Only thread 0 will do the job sequentially                                  */
      coarptr->fbegptr (thrdptr);                 /* Call sequential mating routine                                              */
      coarptr->coarvertnbr = thrdptr->coarvertnbr; /* Global number of coarse vertices is that computed by (sequential) thread 0 */
    }

    threadBarrier (thrdptr);                      /* coarptr->coarvertnbr must be known to all */
  }
  else
#endif /* GRAPHCOARSENTHREAD */
#endif /* (defined GRAPHMATCHTHREAD) && ! (defined SCOTCH_DETERMINISTIC) */
  {
    coarptr->fbegptr (thrdptr);                   /* Call sequential mating routine                                             */
    coarptr->coarvertnbr = thrdptr->coarvertnbr;  /* Global number of coarse vertices is that computed by (sequential) thread 0 */
  }

#ifdef SCOTCH_DEBUG_GRAPH2
  for (finevertnum = finevertbas; finevertnum < finevertnnd; finevertnum ++) {
    if (finematetax[finevertnum] == ~0) {         /* If matching not aborted, this should not happen */
      errorPrint ("graphMatch: internal error");
      coarptr->coarvertnbr = coarptr->coarvertmax;
    }
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */
}
