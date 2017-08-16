/* Copyright 2004,2007,2009,2011-2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_coarsen.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source graph   **/
/**                coarsening functions.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     31 oct 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     08 jun 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     17 sep 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to     31 aug 2005     **/
/**                # Version 5.0  : from : 13 dec 2006     **/
/**                                 to     24 mar 2008     **/
/**                # Version 5.1  : from : 30 oct 2009     **/
/**                                 to     30 oct 2009     **/
/**                # Version 6.0  : from : 09 mar 2011     **/
/**                                 to     28 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_COARSEN

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "graph_match.h"

/***************************/
/*                         */
/* The coarsening routine. */
/*                         */
/***************************/

#ifdef GRAPHCOARSENTHREAD

/* This routine aggregates a sum and max
** reduction of partial coarse graph
** parameters computed by multiple
** threads.
*/

static
void
graphCoarsenReduce (
GraphCoarsenThread * restrict const tlocptr,      /* Pointer to local thread */
void * restrict const               vlocptr,      /* Pointer to local value  */
void * restrict const               vremptr)      /* Pointer to remote value */
{
  GraphCoarsenThread * restrict const tremptr = (GraphCoarsenThread *) vremptr;

  tlocptr->coaredloadj += tremptr->coaredloadj;   /* Sum edge load sum adjustments */
  if (tremptr->coardegrmax > tlocptr->coardegrmax) /* Take maximum of degrees      */
    tlocptr->coardegrmax = tremptr->coardegrmax;
}

/* This routine performs a perfix scan
** sum operation on a single Gnum value.
*/

static
void
graphCoarsenScan (
GraphCoarsenThread * restrict const tlocptr,      /* Pointer to local thread */
Gnum * restrict const               vlocptr,      /* Pointer to local value  */
Gnum * restrict const               vremptr,      /* Pointer to remote value */
const int                           phasval)      /* Phase index             */
{
  vlocptr[1 - phasval] = vlocptr[phasval] + ((vremptr == NULL) ? 0 : vremptr[phasval]);
}

#endif /* GRAPHCOARSENTHREAD */

/* This routine is the threaded core of the building
** of the coarse graph from the fine graph.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

static
int
graphCoarsen2 (
void *                      dataptr)
{
  Gnum                          baseval;
  Gnum                          finevertbas;
  Gnum                          finevertnnd;
  Gnum                          finevertnbr;
  Gnum                          finevertnum;
  Gnum                          coarvertnbr;
  Gnum                          coarvertnum;
  Gnum                          coarhashnbr;      /* Size of neighbor vertex hash table */
  GraphCoarsenMulti * restrict  coarmulttax;
  Gnum                          coarmultsiz;      /* Size of embedded multinode array   */
#ifdef GRAPHCOARSENTHREAD
  int                           thrdnbr;
  int                           thrdnum;
#endif /* GRAPHCOARSENTHREAD */

  GraphCoarsenThread * restrict const         thrdptr = (GraphCoarsenThread *) dataptr;
  volatile GraphCoarsenData * restrict const  coarptr = (GraphCoarsenData *) (thrdptr->thrddat.grouptr);
  const Graph * restrict const                finegrafptr = coarptr->finegrafptr;
  volatile Gnum * restrict const              finecoartax = coarptr->finematetax;
  volatile Graph * restrict const             coargrafptr = coarptr->coargrafptr;

  graphMatch (thrdptr);                           /* Perform (threaded) matching */

  coarvertnbr = coarptr->coarvertnbr;             /* Get number of vertices actually created            */
  if (coarvertnbr >= coarptr->coarvertmax)        /* If matching failed or if coarsened graph too large */
    return  (1);                                  /* Do not proceed any further                         */

  coarmultsiz = (coarptr->coarmulttab == NULL) ? coarvertnbr : 0; /* Determine whether coarmulttab is iser-provided */

  baseval = finegrafptr->baseval;
#ifdef GRAPHCOARSENTHREAD
  thrdnbr = coarptr->thrddat.thrdnbr;
  thrdnum = thrdptr->thrddat.thrdnum;

  if (thrdnum == 0)                               /* Thread 0 populates the graph data structure */
#endif /* GRAPHCOARSENTHREAD */
  {
    memSet (coargrafptr, 0, sizeof (Graph));      /* Initialize coarse graph on thread 0 */
    coargrafptr->flagval = GRAPHFREEVERT | GRAPHVERTGROUP | GRAPHEDGEGROUP;
    coargrafptr->baseval = baseval;
    coargrafptr->vertnbr = coarvertnbr;
    coargrafptr->vertnnd = coarvertnbr + baseval;
    coargrafptr->velosum = finegrafptr->velosum;    /* Keep load of finer graph */
    if (memAllocGroup ((void **) (void *)
                       &coargrafptr->verttax, (size_t) ((coarvertnbr + 1)    * sizeof (Gnum)),
                       &coargrafptr->velotax, (size_t) (coarvertnbr          * sizeof (Gnum)),
                       &coarmulttax,          (size_t) (coarmultsiz          * sizeof (GraphCoarsenMulti)),
                       &coargrafptr->edgetax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)), /* Pre-allocate space for edge arrays */ 
                       &coargrafptr->edlotax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("graphCoarsen2: out of memory (1)"); /* Allocate coarser graph structure */
      return     (2);
    }
    if (coarmultsiz > 0)                          /* If array created internally, record its location */
      coarptr->coarmulttab = coarmulttax;         /* Record un-based array location                   */
    coarmulttax = coarptr->coarmulttab - baseval; /* Only thread 0 knows coarptr->coarmulttab         */

    coargrafptr->verttax -= baseval;              /* Base coarse graph arrays */
    coargrafptr->velotax -= baseval;
    coargrafptr->edgetax -= baseval;
    coargrafptr->edlotax -= baseval;
  }

  finevertnnd = thrdptr->finevertnnd;             /* Will be used by many loops */
#ifdef GRAPHCOARSENTHREAD
  if (thrdnbr > 1) {                              /* If more than one thread */
    Gnum                coarvertnnd;

    for (finevertnum = thrdptr->finevertbas, coarvertnum = 0;
         finevertnum < finevertnnd; finevertnum ++) {
      Gnum                finematenum;            /* Number of current mate vertex */

      finematenum = finecoartax[finevertnum];     /* Get mate number                  */
      if (finematenum >= finevertnum)             /* If mate has larger number        */
        coarvertnum ++;                           /* One more local multinode created */
    }

    thrdptr->coarvertbas = (thrdnum == 0) ? (coarvertnum + baseval) : coarvertnum;
    threadScan (thrdptr, &thrdptr->coarvertbas, (ThreadScanFunc) graphCoarsenScan); /* Compute start indices for multinodes; barrier for coarptr->coarmulttab */
    coarvertnum = thrdptr->coarvertbas - coarvertnum;
    coarmulttax = coarptr->coarmulttab - baseval; /* All threads know coarptr->coarmulttab */

    for (finevertnum = thrdptr->finevertbas;
         finevertnum < finevertnnd; finevertnum ++) {
      Gnum                finematenum;            /* Number of current mate vertex */

      finematenum = finecoartax[finevertnum];     /* Get mate number                      */
      if (finematenum >= finevertnum) {           /* If mate has larger number            */
        coarmulttax[coarvertnum].vertnum[0] = finevertnum; /* Build new multinode         */
        coarmulttax[coarvertnum].vertnum[1] = finematenum; /* Second index always biggest */
        coarvertnum ++;                           /* One more local multinode created     */
      }
    }

    thrdptr->coarvertbas = DATASCAN (coargrafptr->vertnbr, thrdnbr, thrdnum) + baseval; /* Set bounds for coarse vertex processing */
    thrdptr->coarvertnnd = DATASIZE (coargrafptr->vertnbr, thrdnbr, thrdnum) + thrdptr->coarvertbas;

    threadBarrier (thrdptr);                      /* Ensure all of coarmulttax has been written */

    for (coarvertnum = thrdptr->coarvertbas, coarvertnnd = thrdptr->coarvertnnd;
         coarvertnum < coarvertnnd; coarvertnum ++) {
      finecoartax[coarmulttax[coarvertnum].vertnum[0]] = /* Build fine-to-coarse array */
      finecoartax[coarmulttax[coarvertnum].vertnum[1]] = coarvertnum;
    }
  }
  else
#endif /* GRAPHCOARSENTHREAD */
  {
    for (finevertnum = thrdptr->finevertbas, coarvertnum = baseval; /* Finalize finecoartab array */
         finevertnum < finevertnnd; finevertnum ++) {
      Gnum                finematenum;            /* Number of current mate vertex */

      finematenum = finecoartax[finevertnum];     /* Get mate number                               */
      if (finematenum >= finevertnum) {           /* If mate has larger number                     */
        coarmulttax[coarvertnum].vertnum[0] = finevertnum; /* Build new multinode                  */
        coarmulttax[coarvertnum].vertnum[1] = finematenum; /* Second index always biggest          */
        finecoartax[finematenum] =                /* Point to coarse vertex                        */
        finecoartax[finevertnum] = coarvertnum;   /* Always valid since coarvertnum <= finevertnum */
        coarvertnum ++;                           /* One more multinode created                    */
      }
    }

    thrdptr->coarvertbas = baseval;               /* Set bounds for coarse vertex processing */
    thrdptr->coarvertnnd = coarvertnbr + baseval;
  }

  coarhashnbr = coarptr->coarhashmsk + 1;
  if ((thrdptr->coarhashtab = memAlloc (coarhashnbr * sizeof (GraphCoarsenHash))) == NULL) { /* Allocate local thread memory */
    errorPrint ("graphCoarsen2: out of memory (2)");
    return     (2);
  }
  memSet (thrdptr->coarhashtab, ~0, coarhashnbr * sizeof (GraphCoarsenHash)); /* Initialize (local) hash table */

#ifdef GRAPHCOARSENTHREAD
  if (thrdnbr > 1) {                              /* If more than one thread */
    Gnum                coaredgenbr;

    threadBarrier (thrdptr);                      /* Ensure all of finecoartax has been written */

    thrdptr->coaredgebas = 0;                     /* No coarse edges accounted for yet            */
    graphCoarsenEdgeCt (thrdptr);                 /* Count number of coarse edges for each thread */

    coaredgenbr = thrdptr->coaredgebas;           /* Save number of local coarse edges       */
    if (thrdnum == 0)                             /* Prepare start index for edge index scan */
      thrdptr->coaredgebas = coaredgenbr + baseval;
    threadScan (thrdptr, &thrdptr->coaredgebas, (ThreadScanFunc) graphCoarsenScan); /* Compute scan on coarse edge indices */
    thrdptr->coaredgebas -= coaredgenbr;          /* Adjust value to have real edge start index                            */

    memSet (thrdptr->coarhashtab, ~0, coarhashnbr * sizeof (GraphCoarsenHash)); /* Re-initialize (local) hash table */
  }
  else
#endif /* GRAPHCOARSENTHREAD */
    thrdptr->coaredgebas = baseval;               /* We start from the beginning */

  ((finegrafptr->edlotax != NULL) ? graphCoarsenEdgeLl : graphCoarsenEdgeLu) (thrdptr); /* Build coarse graph edge array */

  memFree (thrdptr->coarhashtab);                 /* Free used (local) hash table */

#ifdef GRAPHCOARSENTHREAD
  if (thrdnbr > 1)
    threadReduce (thrdptr, thrdptr, (ThreadReduceFunc) graphCoarsenReduce, 0); /* Sum edloadj and get maximum of degrmax */
  if (thrdnum == 0)
#endif /* GRAPHCOARSENTHREAD */
  {
    coargrafptr->edlosum = thrdptr->coaredloadj + finegrafptr->edlosum;
    coargrafptr->degrmax = thrdptr->coardegrmax;
  }
#ifdef GRAPHCOARSENTHREAD
  if (thrdnum == (thrdnbr - 1))
#endif /* GRAPHCOARSENTHREAD */
    coargrafptr->verttax[coargrafptr->vertnnd] = thrdptr->coaredgebas; /* Mark end of edge array */

  return (0);                                     /* Joining all treads will serve as synchronization for coarse graph data */
}

/* This routine is the sequential core of the
** matching and coarse graph building process.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

static
int
graphCoarsen3 (
GraphCoarsenData * restrict const     coarptr,
GraphCoarsenMulti * restrict * const  coarmultptr) /*+ Pointer to un-based multinode table to build +*/
{
  Gnum                          coarvertnbr;      /* Number of coarse vertices                */
  Gnum                          coarvertnum;      /* Number of current multinode vertex       */
  Gnum                          coarvfixnbr;      /* Coarse number of fixed vertices          */
  GraphCoarsenMulti * restrict  coarmulttax;      /* Multinode array                          */
  Gnum                          coarmultsiz;      /* Size of embedded multinode array         */
  Gnum *                        finecoartab;      /* Fine vertex mating / indexing array      */
  Gnum                          finevertnum;      /* Number of currently selected fine vertex */
  Gnum                          coarhashmsk;      /* Mask for access to hash table            */
  size_t                        coarmultoftval;
  size_t                        coarvelooftval;
  size_t                        coaredgeoftval;
  size_t                        coaredlooftval;
#ifdef GRAPHCOARSENTHREAD
  int                           thrdnbr;
#endif /* GRAPHCOARSENTHREAD */
  int                           o;

  Graph * restrict const        coargrafptr = coarptr->coargrafptr;
  const Graph * restrict const  finegrafptr = coarptr->finegrafptr;
  const Gnum                    finevertnbr = finegrafptr->vertnbr;
  const Gnum                    baseval     = finegrafptr->baseval;

  for (coarhashmsk = 31; coarhashmsk < finegrafptr->degrmax; coarhashmsk = coarhashmsk * 2 + 1) ; /* Compute size of hash table */
  coarptr->coarhashmsk = coarhashmsk * 4 + 3;     /* Record it for (local) hash table allocation */

#ifdef GRAPHCOARSENTHREAD
  thrdnbr                  =
  coarptr->thrddat.thrdnbr = SCOTCH_PTHREAD_NUMBER; /* Required by graphMatchInit */
#else /* GRAPHCOARSENTHREAD */
  coarptr->thrddat.thrdnbr = 1;
#endif /* GRAPHCOARSENTHREAD */

  if (coarptr->finematetax == NULL) {             /* If no user-provided mating array           */
    if (graphMatchInit (coarptr) != 0)            /* Initialize global data needed for matching */
      return (1);

    if ((finecoartab = (Gnum *) memAlloc (finevertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("graphCoarsen3: out of memory (1)"); /* Allocate coarse graph mating and indexing array */
      return     (2);
    }
    coarptr->finematetax = finecoartab - baseval;   /* Set based access to finematetab / finecoartab */
  }
  else {
    graphMatchNone (coarptr);                     /* Initialize global data to avoid matching */
    finecoartab = NULL;                           /* No need to allocate local mating array   */
  }

  coarptr->coarmulttab = *coarmultptr;

#ifdef GRAPHCOARSENTHREAD
  if (thrdnbr > 1) {
    GraphCoarsenThread * restrict thrdtab;
    int                           thrdnum;
    Gnum                          finevertbas;

    if ((thrdtab = memAlloc (thrdnbr * sizeof (GraphCoarsenThread))) == NULL) {
      errorPrint ("graphCoarsen3: out of memory (2)");
      memFree    (finecoartab);
      return     (2);
    }

    for (thrdnum = 0, finevertbas = baseval;
         thrdnum < thrdnbr; thrdnum ++) {
      thrdtab[thrdnum].randval = intRandVal (INT_MAX);
      thrdtab[thrdnum].finevertbas = finevertbas;
      thrdtab[thrdnum].finevertnnd = finevertbas += DATASIZE (finevertnbr, thrdnbr, thrdnum);
      thrdtab[thrdnum].coarvertnbr = 0;           /* No coarse vertices yet */
    }

    o = threadLaunch (coarptr, thrdtab, sizeof (GraphCoarsenThread), (ThreadLaunchStartFunc) graphCoarsen2, (ThreadLaunchJoinFunc) NULL,
                      thrdnbr, THREADCANBARRIER | THREADCANREDUCE);

    memFree (thrdtab);                            /* Free group leader */
  }
  else
#endif /* GRAPHCOARSENTHREAD */
  {
    GraphCoarsenThread  thrddat;

#ifdef GRAPHCOARSENTHREAD
    thrddat.thrddat.thrdnum = 0;                  /* Thread 0 of 1 */
#endif /* GRAPHCOARSENTHREAD */
    thrddat.thrddat.grouptr = (void *) coarptr;
    thrddat.randval         = intRandVal (INT_MAX);
    thrddat.finevertbas     = baseval;
    thrddat.finevertnnd     = baseval + finevertnbr;

    o = graphCoarsen2 (&thrddat);
  }

  if (finecoartab != NULL)
    memFree (finecoartab);

  if (o != 0)                                     /* If coarsened graph is too small, abort here */
    return (1);

  coargrafptr->edgenbr = coargrafptr->verttax[coargrafptr->vertnnd] - baseval; /* Set exact number of edges */

  coarvertnbr = coargrafptr->vertnbr;
  coarmultsiz = (*coarmultptr == NULL) ? coarvertnbr : 0; /* Tell whether we have to resize multloctab within coarse graph vertex group */
  coarvelooftval = coargrafptr->velotax - coargrafptr->verttax;
  coaredgeoftval = coargrafptr->edgetax - coargrafptr->verttax;
  coaredlooftval = coargrafptr->edlotax - coargrafptr->verttax;
  coarmultoftval = (Gnum *) coarptr->coarmulttab - coargrafptr->verttax;
  if (memReallocGroup ((void *) (coargrafptr->verttax + baseval), /* Re-allocate data, wiping temporary arrays */
                       &coargrafptr->verttax, (size_t) ((coarvertnbr + 1)    * sizeof (Gnum)),
                       &coargrafptr->velotax, (size_t) (coarvertnbr          * sizeof (Gnum)),
                       &coarptr->coarmulttab, (size_t) (coarmultsiz          * sizeof (GraphCoarsenMulti)),
                       &coargrafptr->edgetax, (size_t) (finegrafptr->edgenbr * sizeof (Gnum)),
                       &coargrafptr->edlotax, (size_t) (coargrafptr->edgenbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("graphCoarsen3: cannot reallocate memory"); /* Allocate coarser graph structure */
    return (2);
  }
  coargrafptr->verttax -= baseval;
  coargrafptr->vendtax  = coargrafptr->verttax + 1; /* Use compact representation of arrays */
  coargrafptr->velotax  = coargrafptr->verttax + coarvelooftval;
  coargrafptr->edgetax  = coargrafptr->verttax + coaredgeoftval;
  coargrafptr->edlotax  = coargrafptr->verttax + coaredlooftval;
  if (coarmultsiz > 0)                            /* If multinode array not user-provided                                         */
    *coarmultptr = ((GraphCoarsenMulti *) (coargrafptr->verttax + coarmultoftval)); /* Return pointer to un-based multinode array */
  if (coarptr->coarvfixptr != NULL)
    *coarptr->coarvfixptr = coarvfixnbr = coarptr->finevfixnbr; /* TODO: compute real number ! */

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (coargrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("graphCoarsen3: inconsistent graph data");
    graphFree  (coargrafptr);
    return     (2);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}

/* This routine coarsens the given "finegraph" into
** "coargraph", as long as the coarsening ratio remains
** below some threshold value and the coarsened graph
** is not too small.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

int
graphCoarsen (
const Graph * restrict const          finegrafptr, /*+ Graph to coarsen                             +*/
Graph * restrict const                coargrafptr, /*+ Coarse graph to build                        +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to un-based multinode table to build +*/
const Gnum                            coarvertnbr, /*+ Minimum number of coarse vertices            +*/
const double                          coarval,     /*+ Maximum contraction ratio                    +*/
const Anum * restrict const           fineparotax,
const Anum * restrict const           finepfixtax,
const Gnum                            finevfixnbr,
Gnum * restrict const                 coarvfixptr)
{
  GraphCoarsenData              coardat;          /* Graph coarsening global data */

#ifdef SCOTCH_DEBUG_GRAPH1
  if (coarval < 0.5L)                             /* If impossible coarsening ratio wanted */
    return (1);                                   /* We will never succeed                 */
#endif /* SCOTCH_DEBUG_GRAPH1 */

  coardat.coarvertmax = (Gnum) ((double) (finegrafptr->vertnbr - finevfixnbr) * coarval) + finevfixnbr; /* Maximum number of coarse vertices */
  if (coardat.coarvertmax < coarvertnbr)          /* If there will be too few vertices in graph */
    return (1);                                   /* It is useless to go any further            */

  coardat.finegrafptr = finegrafptr;              /* Fill caller part of matching data structure */
  coardat.fineparotax = fineparotax;
  coardat.finepfixtax = finepfixtax;
  coardat.finevfixnbr = finevfixnbr;
  coardat.finematetax = NULL;                     /* No user-provided mating array */
  coardat.coargrafptr = coargrafptr;
  coardat.coarvfixptr = coarvfixptr;

  return (graphCoarsen3 (&coardat, coarmultptr));
}

/* This routine coarsens the given "finegraph" into
** "coargraph", as long as the coarsening ratio remains
** below some threshold value and the coarsened graph
** is not too small.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

int
graphCoarsenBuild (
const Graph * restrict const          finegrafptr, /*+ Graph to coarsen                               +*/
Graph * restrict const                coargrafptr, /*+ Coarse graph to build                          +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to un-based multinode table to build   +*/
const Gnum                            coarvertnbr, /*+ User-provided number of coarse vertices        +*/
Gnum * restrict const                 finematetab) /*+ Pointer to un-based user-provided mating array +*/
{
  GraphCoarsenData              coardat;          /* Graph coarsening global data */

  coardat.finegrafptr = finegrafptr;              /* Fill caller part of matching data structure */
  coardat.fineparotax = NULL;                     /* TODO: arrays not handled yet                */
  coardat.finepfixtax = NULL;
  coardat.finevfixnbr = 0;
  coardat.finematetax = finematetab - finegrafptr->baseval; /* User-provided mating array */
  coardat.coargrafptr = coargrafptr;
  coardat.coarvertmax = finegrafptr->vertnbr + 1; /* Value that always succeeds    */
  coardat.coarvertnbr = coarvertnbr;              /* Set number of coarse vertices */
  coardat.coarvfixptr = NULL;                     /* TODO: arrays not handled yet  */

  return (graphCoarsen3 (&coardat, coarmultptr));
}

/****************************************/
/*                                      */
/* The edge array building subroutines. */
/*                                      */
/****************************************/

#define GRAPHCOARSENEDGENAME        graphCoarsenEdgeLl
#define GRAPHCOARSENEDGEINIT        const Gnum * restrict const fineedlotax = finegrafptr->edlotax
#define GRAPHCOARSENEDGEEDLOINIT    coaredlotax[coaredgenum] = fineedlotax[fineedgenum]
#define GRAPHCOARSENEDGEEDLOADD     coaredlotax[coarhashtab[h].edgenum] += fineedlotax[fineedgenum]
#define GRAPHCOARSENEDGEEDLOSUB     coaredloadj -= fineedlotax[fineedgenum]
#include "graph_coarsen_edge.c"
#undef GRAPHCOARSENEDGENAME
#undef GRAPHCOARSENEDGEINIT
#undef GRAPHCOARSENEDGEEDLOINIT
#undef GRAPHCOARSENEDGEEDLOADD
#undef GRAPHCOARSENEDGEEDLOSUB

#define GRAPHCOARSENEDGENAME        graphCoarsenEdgeLu
#define GRAPHCOARSENEDGEINIT
#define GRAPHCOARSENEDGEEDLOINIT    coaredlotax[coaredgenum] = 1
#define GRAPHCOARSENEDGEEDLOADD     coaredlotax[coarhashtab[h].edgenum] ++
#define GRAPHCOARSENEDGEEDLOSUB     coaredloadj --
#include "graph_coarsen_edge.c"
#undef GRAPHCOARSENEDGENAME
#undef GRAPHCOARSENEDGEINIT
#undef GRAPHCOARSENEDGEEDLOINIT
#undef GRAPHCOARSENEDGEEDLOADD
#undef GRAPHCOARSENEDGEEDLOSUB

#ifdef GRAPHCOARSENTHREAD
#define GRAPHCOARSENEDGECOUNT                     /* Local coarse edge count routine */
#define GRAPHCOARSENEDGENAME        graphCoarsenEdgeCt
#define GRAPHCOARSENEDGEINIT
#define GRAPHCOARSENEDGEEDLOINIT
#define GRAPHCOARSENEDGEEDLOADD
#define GRAPHCOARSENEDGEEDLOSUB
#include "graph_coarsen_edge.c"
#undef GRAPHCOARSENEDGENAME
#undef GRAPHCOARSENEDGEINIT
#undef GRAPHCOARSENEDGEEDLOINIT
#undef GRAPHCOARSENEDGEEDLOADD
#undef GRAPHCOARSENEDGEEDLOSUB
#undef GRAPHCOARSENEDGECOUNT
#endif /* GRAPHCOARSENTHREAD */
