/* Copyright 2004,2007-2009,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb_map.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm.       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     08 dec 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     07 nov 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     06 mar 2005     **/
/**                # Version 5.1  : from : 22 nov 2007     **/
/**                                 to     04 feb 2009     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     29 aug 2014     **/
/**                                                        **/
/**   NOTES      : # This code is a complete rewrite of    **/
/**                  the original code of kgraphMapRb(),   **/
/**                  hence the kept history.               **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB_MAP

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_rb.h"
#include "kgraph_map_rb_map.h"

/*
**  The static variables.
*/

static KgraphMapRbMapPoolLink kgraphmaprbmappooldummy; /* Dummy links for pool routines; TRICK */

/************************************/
/*                                  */
/* These routines handle job pools. */
/*                                  */
/************************************/

/* This routine initializes the job pool
** structures.
** It returns:
** - 0   : in case of success.
** - !0  : on error.
*/

static
int
kgraphMapRbMapPoolInit (
KgraphMapRbMapPoolData * restrict const poolptr,
const KgraphMapRbData * restrict const  dataptr)
{
  int                 flagval;

  Mapping * restrict const  mappptr = dataptr->mappptr;

  flagval = 0;
  if (archVar (mappptr->archptr) != 0)
    flagval |= KGRAPHMAPRBMAPARCHVAR;
  if (archPart (mappptr->archptr) != 0) {
    flagval |= KGRAPHMAPRBMAPARCHCMPLT;
    poolptr->polival = KGRAPHMAPRBPOLILEVEL;      /* A simple policy will do            */
    poolptr->grafptr = NULL;                      /* We don't need top level graph data */
  }
  else {
    poolptr->polival = dataptr->paraptr->polival; /* Enforce original policy           */
    poolptr->grafptr = dataptr->grafptr;          /* We will need top-level graph data */
  }
  poolptr->pfixtax = dataptr->pfixtax;

  poolptr->linktab[0].prev =                      /* Initialize doubly linked list as empty, pointing to the dummy element */
  poolptr->linktab[0].next =
  poolptr->linktab[1].prev =
  poolptr->linktab[1].next = &kgraphmaprbmappooldummy;
  poolptr->pooltab[0] = &poolptr->linktab[0];
  poolptr->pooltab[1] = (dataptr->paraptr->flagjobtie != 0) ? &poolptr->linktab[0] : &poolptr->linktab[1];

  if ((poolptr->jobtab = (KgraphMapRbMapJob *) memAlloc (mappptr->domnmax * sizeof (KgraphMapRbMapJob))) == NULL) {
    errorPrint ("kgraphMapRbMapPoolInit: out of memory (2)");
    return     (1);
  }
  poolptr->jobtab[0].poolflag = 0;                /* In case kgraphMapRbPoolExit() is called just afterwards on single-domain mapping */

  poolptr->mappptr = mappptr;

  poolptr->domntab[0] = mappptr->domntab;         /* Use original domain array                   */
  if (dataptr->paraptr->flagmaptie != 0) {        /* If mappings are tied, use same domain array */
    poolptr->domntab[1] = mappptr->domntab;
    flagval |= KGRAPHMAPRBMAPPARTHALF;            /* Updates will only involve half of the vertices */
  }
  else {
    if ((poolptr->domntab[1] = (ArchDom *) memAlloc (mappptr->domnmax * sizeof (ArchDom))) == NULL) {
      errorPrint ("kgraphMapRbMapPoolInit: out of memory (3)");
      memFree    (poolptr->jobtab);
      return     (1);
    }
  }

  poolptr->flagval = flagval;

  return (0);
}

/* This routine frees all of the internal arrays
** involved in the DRB algorithms. Great care
** should be taken that this routine always
** succeeds, whatever part of the algorithm it
** is called from.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolExit (
KgraphMapRbMapPoolData * restrict const poolptr)
{
  Anum                domnnbr;
  Anum                jobnum;

  Mapping * restrict const  mappptr = poolptr->mappptr;

  domnnbr = mappptr->domnnbr;

  for (jobnum = 0; jobnum < domnnbr; jobnum ++) { /* For all potential jobs in both pools             */
    if (poolptr->jobtab[jobnum].poolflag != 0) {  /* If job slot is active                            */
      graphExit (&poolptr->jobtab[jobnum].grafdat); /* Free job graph, if not clone of original graph */
    }
  }

  if (mappptr->domntab != poolptr->domntab[1]) {  /* If current mapping domain array is not original domain array */
    if ((mappptr->flagval & MAPPINGFREEDOMN) != 0) /* If mapping domain array was privately owned, free it        */
      memFree (mappptr->domntab);
    mappptr->flagval |= MAPPINGFREEDOMN;          /* Keep current domain array as private mapping domain array */
    mappptr->domntab  = poolptr->domntab[1];      
  }

  memFree (poolptr->jobtab);
}

/* This routine swaps the internal arrays
** involved in the DRB algorithms.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolSwap (
KgraphMapRbMapPoolData * restrict const poolptr)
{
  KgraphMapRbMapPoolLink *  linktmp;
  ArchDom *                 domntmp;

  linktmp             = poolptr->pooltab[0];
  poolptr->pooltab[0] = poolptr->pooltab[1];
  poolptr->pooltab[1] = linktmp;

  domntmp             = poolptr->domntab[0];
  poolptr->domntab[0] = poolptr->domntab[1];
  poolptr->domntab[1] = domntmp;
}

/* This routine doubles the size all of the arrays
** involved in handling the target architecture,
** to make room for new domains of variable-sized
** architectures.
** It returns:
** - 0   : if resize succeeded.
** - !0  : if out of memory.
*/

static
int
kgraphMapRbMapPoolResize (
KgraphMapRbMapPoolData * restrict const poolptr)
{
  KgraphMapRbMapJob * restrict  jobtab;           /* Pointer to (new) job array      */
  Anum                          domnnbr;          /* Current (max) number of domains */
  Anum                          domnmax;          /* New maximum number of domains   */
  int                           i;

  domnnbr = poolptr->mappptr->domnmax;            /* Current max size     */
  domnmax = domnnbr + (domnnbr >> 2) + 8;         /* Increase size by 25% */

  if ((jobtab = (KgraphMapRbMapJob *) memRealloc (poolptr->jobtab, domnmax * sizeof (KgraphMapRbMapJob))) == NULL) {
    errorPrint ("kgraphMapRbMapPoolResize: out of memory (1)");
    return     (1);
  }
  if (jobtab != poolptr->jobtab) {                /* If job array moved                */
    KgraphMapRbMapJob *       joboldtab;          /* Pointer to start of old job array */
    KgraphMapRbMapJob *       joboldtnd;          /* Pointer to end of old job array   */
    Anum                      jobnum;             /* Temporary job index               */
    intptr_t                  jobdlt;             /* Address delta value               */

    joboldtab = poolptr->jobtab;
    joboldtnd = joboldtab + domnnbr;
    jobdlt = (byte *) jobtab - (byte *) joboldtab; /* Compute delta between new and old addresses */

    for (jobnum = 0; jobnum < domnnbr; jobnum ++) {
      if ((jobtab[jobnum].poollink.prev >= (KgraphMapRbMapPoolLink *) joboldtab) && /* If old pointers within bounds of old array, adjust them */
          (jobtab[jobnum].poollink.prev <  (KgraphMapRbMapPoolLink *) joboldtnd))
        jobtab[jobnum].poollink.prev = (KgraphMapRbMapPoolLink *) ((byte *) jobtab[jobnum].poollink.prev + jobdlt);
      if ((jobtab[jobnum].poollink.next >= (KgraphMapRbMapPoolLink *) joboldtab) &&
          (jobtab[jobnum].poollink.next <  (KgraphMapRbMapPoolLink *) joboldtnd))
        jobtab[jobnum].poollink.next = (KgraphMapRbMapPoolLink *) ((byte *) jobtab[jobnum].poollink.next + jobdlt);
    }
    if (poolptr->linktab[0].next != &kgraphmaprbmappooldummy) /* Update first pool pointer */
      poolptr->linktab[0].next = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->linktab[0].next + jobdlt);
    if (poolptr->pooltab[0] != poolptr->pooltab[1]) { /* If job pools not tied                */
      if (poolptr->linktab[1].next != &kgraphmaprbmappooldummy) /* Update second pool pointer */
        poolptr->linktab[1].next = (KgraphMapRbMapPoolLink *) ((byte *) poolptr->linktab[1].next + jobdlt);
    }

    poolptr->jobtab = jobtab;                     /* Set new memory location of job array */
  }

  i = (poolptr->domntab[1] == poolptr->mappptr->domntab) ? 1 : 0; /* Find which domain array is that of the mapping */
  if (mapResize (poolptr->mappptr, domnmax) != 0) {
    errorPrint ("kgraphMapRbMapPoolResize: out of memory (2)");
    return     (1);
  }
  if (poolptr->domntab[1] != poolptr->domntab[0]) { /* If two domain arrays present */
    ArchDom *           domntab;

    if ((domntab = (ArchDom *) memRealloc (poolptr->domntab[i ^ 1], domnmax * sizeof (ArchDom))) == NULL) { /* Reallocate other domain array */
      errorPrint ("kgraphMapRbMapPoolResize: out of memory (3)");
      return (1);
    }
    poolptr->domntab[i ^ 1] = domntab;            /* Set (possibly new) memory location of other domain array */
  }
  else
    poolptr->domntab[i ^ 1] = poolptr->mappptr->domntab; /* Both domain arrays point to the same (new) location  */
  poolptr->domntab[i] = poolptr->mappptr->domntab; /* Set (possibly new) memory location of mapping domain array */

  return (0);
}

/**********************************************/
/*                                            */
/* These routines handle bipartitioning jobs. */
/*                                            */
/**********************************************/


/* This routine adds a job to pool 1 of the
** given pool data structure.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolAdd (
KgraphMapRbMapPoolLink * restrict const linkptr,
KgraphMapRbMapJob * const               jobptr)
{
  jobptr->poollink.prev = linkptr;                /* Link job in pool: TRICK */
  jobptr->poollink.next = linkptr->next;
  jobptr->poolflag      = 1;                      /* Job is in pool    */
  jobptr->poolptr       = linkptr;                /* Point to the pool */
  linkptr->next->prev   = &jobptr->poollink;
  linkptr->next         = &jobptr->poollink;
}

/* This routine gets the best job available from
** the given pool, according to the given policy.
** It returns:
** - !NULL  : pointer to the job.
** - NULL   : if the pool is empty.
*/

static
KgraphMapRbMapJob *
kgraphMapRbMapPoolGet (
KgraphMapRbMapPoolData * const  poolptr)
{
  KgraphMapRbMapJob * jobbest;                    /* Best job found */
  KgraphMapRbMapJob * jobptr;

  jobbest = (KgraphMapRbMapJob *) poolptr->pooltab[0]->next;  /* Get first job in pool */
  for (jobptr  = jobbest;                         /* For all jobs in pool              */
       jobptr != (KgraphMapRbMapJob *) (void *) &kgraphmaprbmappooldummy;
       jobptr  = (KgraphMapRbMapJob *) jobptr->poollink.next) {
    if (jobptr->priolvl > jobbest->priolvl)       /* If the current job has stronger priority */
      jobbest = jobptr;                           /* Select it as the best job                */
  }

  if (jobbest != (KgraphMapRbMapJob *) (void *) &kgraphmaprbmappooldummy) { /* If job found */
    jobbest->poollink.next->prev = jobbest->poollink.prev; /* Remove it from pool           */
    jobbest->poollink.prev->next = jobbest->poollink.next; /* But do not mark it unused     */
  }
  else                                            /* Dummy job means no job found */
    jobbest = NULL;

  return (jobbest);
}

/* This routine adds a job to the given pool
** as the first bipartitioning job.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolFrst (
KgraphMapRbMapPoolData * const  poolptr,
KgraphMapRbMapJob * const       jobptr)           /* Job to be added */
{
  switch (poolptr->polival) {                     /* Set job priority value */
    case KGRAPHMAPRBPOLIRANDOM :
      jobptr->prioval =
      jobptr->priolvl = intRandVal (INTVALMAX);
      break;
    case KGRAPHMAPRBPOLILEVEL   :
    case KGRAPHMAPRBPOLINGLEVEL :
      jobptr->prioval = jobptr->grafdat.vertnbr;
      jobptr->priolvl = 0;
      break;
    case KGRAPHMAPRBPOLISIZE   :
    case KGRAPHMAPRBPOLINGSIZE :
      jobptr->prioval =
      jobptr->priolvl = jobptr->grafdat.vertnbr;
      break;
#ifdef SCOTCH_DEBUG_KGRAPH2
    default :
      errorPrint ("kgraphMapRbMapPoolFrst: unknown job selection policy");
      jobptr->prioval = 0;
      jobptr->priolvl = 0;
      return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  }

  kgraphMapRbMapPoolAdd (poolptr->pooltab[0], jobptr); /* Add job to pool */
}

/* This routine updates the given job
** table with both of the given subjob
** data.
** This routine can be called only if
** the parent jobs of the vertices to
** be updated still exist.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolUpdt1 (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr,        /* Job to be removed      */
const GraphPart * const         parttax,
KgraphMapRbMapJob * const       jobnewptr,        /* Its only active subjob */
const GraphPart                 partval)
{
  Gnum                prioval;
  Gnum                priolvl;

  priolvl = 0;                                    /* Prepare for neighbor updating methods */

  switch (poolptr->polival) {                     /* Set job priority value */
    case KGRAPHMAPRBPOLIRANDOM :
      prioval =
      priolvl = intRandVal (INTVALMAX);
      break;
    case KGRAPHMAPRBPOLILEVEL :
      priolvl = joboldptr->priolvl + 1;
    case KGRAPHMAPRBPOLINGLEVEL :
      prioval = joboldptr->prioval - 1;
      break;
    case KGRAPHMAPRBPOLISIZE :
      priolvl = jobnewptr->grafdat.vertnbr;
    case KGRAPHMAPRBPOLINGSIZE :
      prioval = jobnewptr->grafdat.vertnbr;
      break;
#ifdef SCOTCH_DEBUG_KGRAPH2
    default :
      errorPrint ("kgraphMapRbMapPoolUpdt1: unknown job selection policy");
      jobnewptr->prioval = 0;
      jobnewptr->priolvl = 0;
      return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  }

  jobnewptr->prioval = prioval;

  if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be updated */
    Gnum                prioold;

    KgraphMapRbMapJob * restrict const  jobtab     = poolptr->jobtab;
    const Anum * restrict const         mapparttax = poolptr->mappptr->parttax; /* Based pointer to mapping part array */
    const Anum * restrict const         toppfixtax = poolptr->pfixtax;
    const Gnum * restrict const         topverttax = poolptr->grafptr->verttax; /* Point to top-level graph arrays */
    const Gnum * restrict const         topvendtax = poolptr->grafptr->vendtax;
    const Gnum * restrict const         topedgetax = poolptr->grafptr->edgetax;

    prioold = joboldptr->prioval;

    if (joboldptr->grafdat.vertnbr < poolptr->grafptr->vertnbr) { /* If subgraph is not top graph, change priority of neighboring jobs of old job */
      Gnum                jobvertnnd;
      Gnum                jobvertnum;

      const Gnum * restrict const jobverttax = joboldptr->grafdat.verttax;
      const Gnum * restrict const jobvendtax = joboldptr->grafdat.vendtax;
      const Gnum * restrict const jobvnumtax = joboldptr->grafdat.vnumtax;

      jobnewptr->poolflag = 0;                    /* TRICK: avoid new job being considered for update */

      for (jobvertnum = joboldptr->grafdat.baseval, jobvertnnd = joboldptr->grafdat.vertnnd;
           jobvertnum < jobvertnnd; jobvertnum ++) {
        Gnum                topvertnum;
        Gnum                topedgenum;

        if (parttax[jobvertnum] == partval)       /* If vertex belongs to part which is still alive */
          continue;                               /* Do not consider update part as removed         */

        topvertnum = jobvnumtax[jobvertnum];      /* If graph is smaller than top graph, then vnumtax must exist */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * restrict  jobnghbptr; /* (Old ?) job of neighbor vertex */
          Gnum                          topvertend;

          topvertend = topedgetax[topedgenum];
          if ((toppfixtax != NULL) && (toppfixtax[topvertend] >= 0)) {
#ifdef SCOTCH_DEBUG_KGRAPH2
            if (mapparttax[topvertend] != ~0) {
              errorPrint ("kgraphMapRbMapPoolUpdt1: internal error (1)");
              return;
            }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
            continue;
          }

          jobnghbptr = &jobtab[mapparttax[topvertend]]; /* Get pointer to neighboring job */

          if ((jobnghbptr->poolflag != 0) &&      /* If neighbor is active                   */
              (jobnghbptr->prioval <= prioold))   /* And had not already a stronger priority */
            jobnghbptr->priolvl ++;               /* Update neighbor priority                */
        }
      }

      jobnewptr->poolflag = 1;                    /* TRICK: new job is active again */
    }

    if (jobnewptr->grafdat.vertnbr < poolptr->grafptr->vertnbr) { /* If subgraph is not top graph, update priority of neighbors of new job only */
      Gnum                jobvertnnd;
      Gnum                jobvertnum;

      const Gnum * restrict const jobverttax = jobnewptr->grafdat.verttax;
      const Gnum * restrict const jobvendtax = jobnewptr->grafdat.vendtax;
      const Gnum * restrict const jobvnumtax = jobnewptr->grafdat.vnumtax;

      for (jobvertnum = jobnewptr->grafdat.baseval, jobvertnnd = jobnewptr->grafdat.vertnnd;
           jobvertnum < jobvertnnd; jobvertnum ++) {
        Gnum                          topvertnum;
        Gnum                          topedgenum;

        topvertnum = jobvnumtax[jobvertnum];      /* For subjobs jobvnumtax always exists */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * restrict  jobnghbptr; /* (Old ?) job of neighbor vertex */
          Gnum                          topvertend;

          topvertend = topedgetax[topedgenum];
          if ((toppfixtax != NULL) && (toppfixtax[topvertend] >= 0)) {
#ifdef SCOTCH_DEBUG_KGRAPH2
            if (mapparttax[topvertend] != ~0) {
              errorPrint ("kgraphMapRbMapPoolUpdt1: internal error (2)");
              return;
            }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
            continue;
          }

          jobnghbptr = &jobtab[mapparttax[topvertend]]; /* Get pointer to neighboring job               */
          if (jobnghbptr == jobnewptr)            /* If it is the current job, do not consider the edge */
            continue;

          if ((jobnghbptr->poolflag == 0) ||      /* If neighbor is not active            */
              (prioval > jobnghbptr->prioval))    /* Or if we have higher priority        */
            priolvl ++;                           /* Increase our priority                */
          else if ((prioval <  jobnghbptr->prioval) && /* Else if neighbor has higher one */
                   (prioold >= jobnghbptr->prioval)) /* Which it did not already have     */
            jobnghbptr->priolvl ++;               /* Update neighbor priority             */
        }
      }
    }
  }

  jobnewptr->priolvl = priolvl;

  kgraphMapRbMapPoolAdd (poolptr->pooltab[1], jobnewptr); /* Add job to pool */
}

static
void
kgraphMapRbMapPoolUpdt2 (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr,        /* Job to be removed */
const GraphPart * const         parttax,
KgraphMapRbMapJob * const       jobnewptr0,       /* Its two subjobs   */
KgraphMapRbMapJob * const       jobnewptr1)
{
  KgraphMapRbMapJob * restrict  jobnewtab[2];
  int                           i;

  jobnewtab[0] = jobnewptr0;
  jobnewtab[1] = jobnewptr1;

  for (i = 1; i >= 0; i --) {
    KgraphMapRbMapJob * jobnewptr;
    Gnum                prioval;
    Gnum                priolvl;

    jobnewptr = jobnewtab[i];                     /* Get concerned subjob */

    priolvl = 0;                                  /* Prepare for neighbor updating methods */

    switch (poolptr->polival) {                   /* Set job priority value */
      case KGRAPHMAPRBPOLIRANDOM :
        prioval =
        priolvl = intRandVal (INTVALMAX);
        break;
      case KGRAPHMAPRBPOLILEVEL :
        priolvl = joboldptr->priolvl + 1;
      case KGRAPHMAPRBPOLINGLEVEL :
        prioval = joboldptr->prioval - 1;
        break;
      case KGRAPHMAPRBPOLISIZE :
        priolvl = jobnewptr->grafdat.vertnbr;
      case KGRAPHMAPRBPOLINGSIZE :
        prioval = jobnewptr->grafdat.vertnbr;
        break;
#ifdef SCOTCH_DEBUG_KGRAPH2
      default :
        errorPrint ("kgraphMapRbMapPoolUpdt2: unknown job selection policy");
        jobnewptr->prioval = 0;
        jobnewptr->priolvl = 0;
        return;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }

    jobnewptr0->prioval = prioval + 1;            /* TRICK: when processing subdomain 1, subdomain 0 has higher priority value */
    jobnewptr->prioval  = prioval;                /* Then in its turn subdomain 0 will have its proper priority value          */

    if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be updated */
      Gnum                jobvertnnd;
      Gnum                jobvertnum;
      Gnum                prioold;

      KgraphMapRbMapJob * restrict const  jobtab     = poolptr->jobtab;
      const Anum * restrict const         mapparttax = poolptr->mappptr->parttax; /* Based pointer to mapping part array */
      const Anum * restrict const         toppfixtax = poolptr->pfixtax;
      const Gnum * restrict const         topverttax = poolptr->grafptr->verttax; /* Point to top-level graph arrays */
      const Gnum * restrict const         topvendtax = poolptr->grafptr->vendtax;
      const Gnum * restrict const         topedgetax = poolptr->grafptr->edgetax;
      const Gnum * restrict const         jobverttax = jobnewptr->grafdat.verttax;
      const Gnum * restrict const         jobvendtax = jobnewptr->grafdat.vendtax;
      const Gnum * restrict const         jobvnumtax = jobnewptr->grafdat.vnumtax;

      prioold = joboldptr->prioval;

      for (jobvertnum = jobnewptr->grafdat.baseval, jobvertnnd = jobnewptr->grafdat.vertnnd;
           jobvertnum < jobvertnnd; jobvertnum ++) {
        Gnum                          topvertnum;
        Gnum                          topedgenum;

        topvertnum = jobvnumtax[jobvertnum];      /* For subjobs jobvnumtax always exists */

        if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
            (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
          continue;

        for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
          KgraphMapRbMapJob * jobnghbptr;         /* (Old ?) job of neighbor vertex */
          Gnum                topvertend;

          topvertend = topedgetax[topedgenum];
          if ((toppfixtax != NULL) && (toppfixtax[topvertend] >= 0)) {
#ifdef SCOTCH_DEBUG_KGRAPH2
            if (mapparttax[topvertend] != ~0) {
              errorPrint ("kgraphMapRbMapPoolUpdt2: internal error");
              return;
            }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
            continue;
          }

          jobnghbptr = &jobtab[mapparttax[topvertend]]; /* Get pointer to neighboring job */

          if ((jobnghbptr->poolflag != 0)      && /* If neighbor is in active job  */
              (jobnghbptr->prioval >  prioval) && /* Which gained priority over us */
              (jobnghbptr->prioval <= prioold)) {
            jobnghbptr->priolvl ++;               /* Update neighbor priority */
          }
          if ((jobnghbptr->poolflag == 0) ||      /* If neighbor is fully known    */
              (jobnghbptr->prioval < prioval))    /* Or has smaller priority value */
            priolvl ++;                           /* Then we should be processed   */
        }
      }
    }

    jobnewptr->priolvl = priolvl;                 /* Set new priority          */
    kgraphMapRbMapPoolAdd (poolptr->pooltab[1], jobnewptr); /* Add job to pool */
  }
}

/*
** This routine removes the influence of the
** given job from its neighbor jobs.
** It returns:
** - VOID  : in all cases.
*/

static
void
kgraphMapRbMapPoolRemv (
KgraphMapRbMapPoolData * const  poolptr,
const KgraphMapRbMapJob * const joboldptr)        /* Job to be removed */
{
  KgraphMapRbMapJob * restrict  jobtab;
  const Anum * restrict         mapparttax;       /* Based pointer to mapping part array */
  const Gnum * restrict         jobvnumtax;
  const Gnum * restrict         jobverttax;
  const Gnum * restrict         jobvendtax;
  const Gnum * restrict         topverttax;
  const Gnum * restrict         topvendtax;
  const Gnum * restrict         topedgetax;

  if (poolptr->polival >= KGRAPHMAPRBPOLINEIGHBOR) { /* If neighbors have to be modified */
    Gnum                jobvertnnd;
    Gnum                jobvertnum;
    Gnum                prioold;

    KgraphMapRbMapJob * restrict const  jobtab     = poolptr->jobtab;
    const Anum * restrict const         mapparttax = poolptr->mappptr->parttax; /* Based pointer to mapping part array */
    const Anum * restrict const         toppfixtax = poolptr->pfixtax;
    const Gnum * restrict const         topverttax = poolptr->grafptr->verttax; /* Point to top-level graph arrays */
    const Gnum * restrict const         topvendtax = poolptr->grafptr->vendtax;
    const Gnum * restrict const         topedgetax = poolptr->grafptr->edgetax;
    const Gnum * restrict const         jobverttax = joboldptr->grafdat.verttax;
    const Gnum * restrict const         jobvendtax = joboldptr->grafdat.vendtax;
    const Gnum * restrict const         jobvnumtax = joboldptr->grafdat.vnumtax;

    prioold = joboldptr->prioval;

    for (jobvertnum = joboldptr->grafdat.baseval, jobvertnnd = joboldptr->grafdat.vertnnd;
         jobvertnum < jobvertnnd; jobvertnum ++) {
      Gnum                topvertnum;             /* Source graph vertex number */
      Gnum                topedgenum;             /* Source graph edge number   */

      topvertnum = (jobvnumtax == NULL) ? jobvertnum : jobvnumtax[jobvertnum];

      if ((topvendtax[topvertnum] - topverttax[topvertnum]) == /* If vertex is internal, skip it */
          (jobvendtax[jobvertnum] - jobverttax[jobvertnum]))
        continue;

      for (topedgenum = topverttax[topvertnum]; topedgenum < topvendtax[topvertnum]; topedgenum ++) {
        KgraphMapRbMapJob * jobnghbptr;           /* (Old ?) job of neighbor vertex */
        Gnum                topvertend;

        topvertend = topedgetax[topedgenum];
        if ((toppfixtax != NULL) && (toppfixtax[topvertend] >= 0)) {
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (mapparttax[topvertend] != ~0) {
            errorPrint ("kgraphMapRbMapPoolRemv: internal error (1)");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          continue;
        }

        jobnghbptr = &jobtab[mapparttax[topvertend]]; /* Get pointer to neighboring job */

        if ((jobnghbptr->poolflag != 0) &&        /* If neighbor job is active                       */
            (jobnghbptr->prioval <= prioold))     /* And had not already a stronger priority         */
          jobnghbptr->priolvl ++;                 /* Increase its priority since we are now inactive */
      }
    }
  }
}

/**********************************************/
/*                                            */
/* These routines handle the pool part array. */
/*                                            */
/**********************************************/

static
void
kgraphMapRbMapPartBoth (
KgraphMapRbMapPoolData * restrict const poolptr,
const Bgraph * restrict const           actgrafptr,
const Anum * restrict const             jobsubnum)
{
  Gnum                       actvertnum;
  const GraphPart * restrict actparttax;
  Anum * restrict            mapparttax;
  Anum                       mappartval1;
  Anum                       mappartdlt;

  actparttax  = actgrafptr->parttax;
  mapparttax  = poolptr->mappptr->parttax;
  mappartval1 = jobsubnum[1];
  mappartdlt  = jobsubnum[0] - jobsubnum[1];

  if (actgrafptr->s.vnumtax != NULL) {
    const Gnum * restrict      actvnumtax;

    actvnumtax = actgrafptr->s.vnumtax;
    for (actvertnum = actgrafptr->s.baseval; actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (mapparttax[actvnumtax[actvertnum]] == ~0) {
        errorPrint ("kgraphMapRbMapPartBoth: internal error (1)");
        return;
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      mapparttax[actvnumtax[actvertnum]] = mappartval1 + ((((Anum) actparttax[actvertnum]) - 1) & mappartdlt);
    }
  }
  else {
    for (actvertnum = actgrafptr->s.baseval; actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (mapparttax[actvertnum] == ~0) {
        errorPrint ("kgraphMapRbMapPartBoth: internal error (2)");
        return;
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      mapparttax[actvertnum] = mappartval1 + ((((Anum) actparttax[actvertnum]) - 1) & mappartdlt);
    }
  }
}

static
void
kgraphMapRbMapPartOne (
KgraphMapRbMapPoolData * restrict const poolptr,
const Bgraph * restrict const           actgrafptr,
const Anum                              jobsubnum1)
{
  Gnum                       actvertnum;
  const GraphPart * restrict actparttax;
  Anum * restrict            mapparttax;

  actparttax  = actgrafptr->parttax;
  mapparttax  = poolptr->mappptr->parttax;

  if (actgrafptr->s.vnumtax != NULL) {
    const Gnum * restrict      actvnumtax;

    actvnumtax = actgrafptr->s.vnumtax;
    for (actvertnum = actgrafptr->s.baseval; actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
      if (actparttax[actvertnum] == 1) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (mapparttax[actvnumtax[actvertnum]] == ~0) {
          errorPrint ("kgraphMapRbMapPartOne: internal error (1)");
          return;
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        mapparttax[actvnumtax[actvertnum]] = jobsubnum1;
      }
    }
  }
  else {
    for (actvertnum = actgrafptr->s.baseval; actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
      if (actparttax[actvertnum] == 1) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (mapparttax[actvertnum] == ~0) {
          errorPrint ("kgraphMapRbMapPartOne: internal error (2)");
          return;
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        mapparttax[actvertnum] = jobsubnum1;
      }
    }
  }
}

/********************************************/
/*                                          */
/* This is the entry point for the Dual     */
/* Recursive Bipartitioning mapping method. */
/*                                          */
/********************************************/

/* This routine runs the Dual Recursive
** Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbMap (
const KgraphMapRbData * restrict const  dataptr,  /*+ Global mapping data                  +*/
const Graph * restrict const            grafptr,  /*+ Graph to map, without fixed vertices +*/
const Anum                              vflonbr,  /*+ Number of fixed vertex load slots    +*/
KgraphMapRbVflo * restrict const        vflotab)  /*+ Array of fixed vertex load slots     +*/
{
  KgraphMapRbMapPoolData  pooldat;                /* Data for handling jobs and job pools */
  ArchDom                 domnsubtab[2];          /* Subdomains of current job domain     */
  KgraphMapRbMapJob       joborgdat;              /* Aera to save original job data       */
  Anum                    jobsubnum[2];           /* Number of subjob slots in job array  */
  Gnum                    jobsubsiz[2];           /* Sizes of subjobs                     */
  Bgraph                  actgrafdat;             /* Bipartition graph                    */
  double                  comploadmin;            /* Minimum vertex load per target load  */
  double                  comploadmax;            /* Maximum vertex load per target load  */
  int                     i;

  Mapping * restrict const  mappptr = dataptr->mappptr;

  mapFrst (mappptr);                              /* Initialize mapping */
#ifdef SCOTCH_DEBUG_KGRAPH2
  if (dataptr->pfixtax != NULL) {                 /* In debug mode, fixed vertex parts are set to ~0 */
    Gnum                vertnnd;
    Gnum                vertnum;

    Anum * restrict const       parttax = mappptr->parttax;
    const Anum * restrict const pfixtax = dataptr->pfixtax;

    for (vertnum = dataptr->grafptr->baseval, vertnnd = dataptr->grafptr->vertnnd; vertnum < vertnnd; vertnum ++) {
      if (pfixtax[vertnum] >= 0)
        parttax[vertnum] = ~0;
    }
  }
  mappptr->domnmax = 1;                           /* Force resizing of job arrays, for debugging */
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  if (kgraphMapRbMapPoolInit (&pooldat, dataptr) != 0) /* Initialize pool data; done first for kgraphMapRbMapPoolExit() to succeed afterwards */
    return (1);

  if ((((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) && (archDomSize (mappptr->archptr, &mappptr->domnorg) <= 1)) || /* If single-vertex domain   */
      (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) && (grafptr->vertnbr <= 1))) { /* Or if variable-sized architecture with single vertex graph */
    kgraphMapRbMapPoolExit (&pooldat);
    return (0);                                   /* Job already done */
  }

  pooldat.jobtab[0].domnorg = mappptr->domnorg;   /* Build first job                         */
  pooldat.jobtab[0].grafdat = *grafptr;           /* Clone induced graph as first job graph  */
  pooldat.jobtab[0].grafdat.flagval &= ~GRAPHFREETABS; /* Do not free its arrays on exit     */
  pooldat.jobtab[0].vflonbr = vflonbr;            /* Record initial list of fixed load slots */
  pooldat.jobtab[0].vflotab = vflotab;
  kgraphMapRbMapPoolFrst (&pooldat, &pooldat.jobtab[0]); /* Add initial job */

  comploadmin = (1.0 - dataptr->paraptr->kbalval) * dataptr->comploadrat; /* Ratio can have been tilted when working on subgraph */
  comploadmax = (1.0 + dataptr->paraptr->kbalval) * dataptr->comploadrat;

  while (! kgraphMapRbMapPoolEmpty (&pooldat)) {  /* For all non-empty pools */
    KgraphMapRbMapJob * joborgptr;                /* Pointer to current job  */

    while ((joborgptr = kgraphMapRbMapPoolGet (&pooldat)) != NULL) { /* For all jobs in pool */
      Gnum                vflonbrtab[2];
      Gnum                vflowgttab[2];
      int                 partval;

      jobsubnum[0] = joborgptr - pooldat.jobtab;  /* Get current (and first son) job slot number before possible move of pointers */
      joborgdat = *joborgptr;                     /* Save current job data (clone graph)                                          */

      if (archDomBipart (mappptr->archptr, &joborgdat.domnorg, &domnsubtab[0], &domnsubtab[1]) != 0) {
        errorPrint ("kgraphMapRbMap: cannot bipartition domain");
        kgraphMapRbMapPoolExit (&pooldat);        /* Copied graph will be freed as not yet removed */
        return (1);
      }

      kgraphMapRbVfloSplit (mappptr->archptr, domnsubtab, /* Split fixed vertex load slots, if any */
                            joborgdat.vflonbr, joborgdat.vflotab, vflonbrtab, vflowgttab);
      if (kgraphMapRbBgraph (dataptr, &actgrafdat, &joborgdat.grafdat, pooldat.mappptr, domnsubtab, vflowgttab) != 0) { /* Create bipartition graph */
        errorPrint ("kgraphMapRbMap: cannot create bipartition graph");
        kgraphMapRbMapPoolExit (&pooldat);        /* Copied graph will be freed as not yet removed */
        return (1);
      }

      actgrafdat.s.flagval |= (joborgdat.grafdat.flagval & GRAPHFREETABS); /* Bipartition graph is responsible for freeing the cloned graph data fields */
      joborgptr->poolflag = 0;                    /* Original slot is now considered unused so that cloned graph data will not be freed twice           */

      if ((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) { /* If not variable-sized, impose constraints on bipartition */
        double              comploadavg;

        comploadavg = (double) actgrafdat.s.velosum / (double) archDomWght (mappptr->archptr, &joborgdat.domnorg);
        actgrafdat.compload0min = actgrafdat.compload0avg -
                                  (Gnum) MIN ((comploadmax - comploadavg) * (double) actgrafdat.domnwght[0],
                                              (comploadavg - comploadmin) * (double) actgrafdat.domnwght[1]);
        actgrafdat.compload0max = actgrafdat.compload0avg +
                                  (Gnum) MIN ((comploadavg - comploadmin) * (double) actgrafdat.domnwght[0],
                                              (comploadmax - comploadavg) * (double) actgrafdat.domnwght[1]);
      }

      if (bgraphBipartSt (&actgrafdat, dataptr->paraptr->strat) != 0) { /* Perform bipartitioning */
        errorPrint             ("kgraphMapRbMap: cannot bipartition job");
        bgraphExit             (&actgrafdat);
        kgraphMapRbMapPoolExit (&pooldat);
        return                 (1);
      }

      if ((partval = 1, actgrafdat.compsize0 == 0) || /* If no bipartition found */
          (partval = 0, actgrafdat.compsize0 == actgrafdat.s.vertnbr)) {
        if (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) || /* If architecture is variable-sized     */
            (archDomSize (mappptr->archptr, &domnsubtab[partval]) <= 1)) { /* Or if domain is terminal   */
          pooldat.domntab[0][jobsubnum[0]] = joborgdat.domnorg; /* Update domain in next pool            */
          kgraphMapRbMapPoolRemv (&pooldat, &joborgdat); /* Remove job from pool as long as graph exists */
          bgraphExit (&actgrafdat);               /* Free bipartitioning data as well as current graph   */
          continue;                               /* Process next job in current pool                    */
        }
        else {                                    /* Re-use job slot and graph for further bipartitioning */
          pooldat.domntab[0][jobsubnum[0]] =      /* Update domain in next pool                           */
          joborgptr->domnorg = domnsubtab[partval]; /* New job takes same graph and non-empty subdomain   */
          joborgptr->vflonbr = vflonbrtab[partval];
          joborgptr->vflotab = joborgdat.vflotab + (partval * vflonbrtab[0]); /* Point to proper sub-array           */
          kgraphMapRbMapPoolUpdt1 (&pooldat, &joborgdat, actgrafdat.parttax, joborgptr, partval); /* Add job to pool */
          actgrafdat.s.flagval &= ~GRAPHFREETABS; /* Since graph will be re-used, never free its internal arrays     */
          bgraphExit (&actgrafdat);               /* Free bipartitioning data                                        */
          continue;                               /* Process next job in current pool                                */
        }
      }

      if ((pooldat.mappptr->domnnbr == pooldat.mappptr->domnmax) && /* If all job slots busy and if cannot resize */
          (kgraphMapRbMapPoolResize (&pooldat) != 0)) {
        errorPrint             ("kgraphMapRbMap: cannot resize structures");
        kgraphMapRbMapPoolExit (&pooldat);
        return                 (1);
      }

      jobsubnum[1] = pooldat.mappptr->domnnbr ++; /* Get slot number of new subdomain */
      jobsubsiz[1] = actgrafdat.s.vertnbr - actgrafdat.compsize0;
      jobsubsiz[0] = actgrafdat.compsize0;

      pooldat.jobtab[jobsubnum[1]].poolflag = 0;  /* Assume that new job is inactive in case of premature freeing                           */
      pooldat.domntab[1][jobsubnum[1]] = joborgdat.domnorg; /* Copy original domain to new subdomain as old mapping shares parttax with new */
      pooldat.domntab[0][jobsubnum[0]] = domnsubtab[0]; /* Set subdomains of second mapping before relinking subjobs in pool                */
      pooldat.domntab[0][jobsubnum[1]] = domnsubtab[1];

      if ((pooldat.flagval & KGRAPHMAPRBMAPPARTHALF) != 0) /* If can only update second half */
        kgraphMapRbMapPartOne (&pooldat, &actgrafdat, jobsubnum[1]);
      else
        kgraphMapRbMapPartBoth (&pooldat, &actgrafdat, jobsubnum);

      for (i = 1; i >= 0; i --) {                 /* For both subdomains */
        KgraphMapRbMapJob * jobsubptr;

        jobsubptr = &pooldat.jobtab[jobsubnum[i]]; /* Point to subdomain job slot                                */
        jobsubptr->poollink.prev =                /* Prevent Valgrind from yelling in kgraphMapRbMapPoolResize() */
        jobsubptr->poollink.next = NULL;
        jobsubptr->prioval =                      /* Prevent Valgrind from yelling in kgraphMapRbMapPoolRemv()/Updt1()/Updt2() */
        jobsubptr->priolvl = 0;

        if ((((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) == 0) && (archDomSize (mappptr->archptr, &domnsubtab[i]) <= 1)) || /* If single-vertex domain  */
            (((pooldat.flagval & KGRAPHMAPRBMAPARCHVAR) != 0) && (jobsubsiz[i] <= 1))) { /* Or if variable-sized architecture with single vertex graph */
          jobsubsiz[i] = 0;                       /* Cancel subjob */
          continue;
        }

        partval = i;                              /* At least this subjob works */

        if (graphInducePart (&actgrafdat.s, actgrafdat.parttax, jobsubsiz[i], (GraphPart) i, &jobsubptr->grafdat) != 0) {
          errorPrint             ("kgraphMapRbMap: cannot create induced subgraph");
          bgraphExit             (&actgrafdat);
          kgraphMapRbMapPoolExit (&pooldat);
          return                 (1);
        }
        jobsubptr->poolflag = 1;                  /* So that graph is freed in case of error on other part */
        jobsubptr->domnorg  = domnsubtab[i];
        jobsubptr->vflonbr  = vflonbrtab[i];
        jobsubptr->vflotab  = joborgdat.vflotab + (i * vflonbrtab[0]); /* Point to proper sub-array */
      }

      if ((jobsubsiz[0] | jobsubsiz[1]) == 0)     /* If both subjobs do not need further processing */
        kgraphMapRbMapPoolRemv (&pooldat, &joborgdat);
      else if (jobsubsiz[1 - partval] == 0)       /* If one of the subjobs only needs further processing */
        kgraphMapRbMapPoolUpdt1 (&pooldat, &joborgdat, actgrafdat.parttax, &pooldat.jobtab[jobsubnum[partval]], (GraphPart) partval);
      else
        kgraphMapRbMapPoolUpdt2 (&pooldat, &joborgdat, actgrafdat.parttax, &pooldat.jobtab[jobsubnum[0]], &pooldat.jobtab[jobsubnum[1]]);

      bgraphExit (&actgrafdat);                   /* Free bipartition graph data */
    }

    kgraphMapRbMapPoolSwap (&pooldat);            /* Swap current and next levels */
  }

  kgraphMapRbMapPoolExit (&pooldat);              /* Free internal structures and propagate back new partition */

  return (0);
}
