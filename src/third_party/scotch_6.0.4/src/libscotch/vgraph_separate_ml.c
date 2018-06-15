/* Copyright 2004,2007,2009,2011,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_ml.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module separates a separator       **/
/**                graph using a multi-level scheme.       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 28 oct 1997     **/
/**                                 to     05 nov 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to     20 mar 2005     **/
/**                # Version 5.1  : from : 11 nov 2009     **/
/**                                 to     11 nov 2009     **/
/**                # Version 6.0  : from : 09 mar 2011     **/
/**                                 to     27 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "graph_coarsen.h"
#include "vgraph.h"
#include "vgraph_separate_ml.h"
#include "vgraph_separate_st.h"

/*********************************************/
/*                                           */
/* The coarsening and uncoarsening routines. */
/*                                           */
/*********************************************/

/* This routine builds a coarser graph from the
** graph that is given on input. The coarser
** graphs differ at this stage from classical
** active graphs as their internal gains are not
** yet computed.
** It returns:
** - 0  : if the coarse graph has been built.
** - 1  : if threshold achieved or on error.
*/

static
int
vgraphSeparateMlCoarsen (
const Vgraph * restrict const         finegrafptr, /*+ Finer graph                                  +*/
Vgraph * restrict const               coargrafptr, /*+ Coarser graph to build                       +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to un-based multinode table to build +*/
const VgraphSeparateMlParam * const   paraptr)    /*+ Method parameters                             +*/
{
  *coarmultptr = NULL;                            /* Allocate coarmulttab along with coarse graph */
  if (graphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr, paraptr->coarnbr, paraptr->coarval, NULL, NULL, 0, NULL) != 0)
    return (1);                                   /* Return if coarsening failed */

  coargrafptr->parttax = NULL;                    /* Do not allocate partition data yet      */
  coargrafptr->frontab = finegrafptr->frontab;    /* Re-use frontier array for coarser graph */
  coargrafptr->levlnum = finegrafptr->levlnum + 1; /* Graph level is coarsening level        */

  return (0);
}

/* This routine propagates the separation of the
** coarser graph back to the finer graph, according
** to the multinode table of collapsed vertices.
** After the separation is propagated, it finishes
** to compute the parameters of the finer graph that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse graph data has been propagated to fine graph.
** - !0  : on error.
*/

static
int
vgraphSeparateMlUncoarsen (
Vgraph * restrict const                   finegrafptr, /*+ Finer graph              +*/
const Vgraph * restrict const             coargrafptr, /*+ Coarser graph            +*/
const GraphCoarsenMulti * restrict const  coarmulttab) /*+ Un-based multinode array +*/
{
  Gnum                coarvertnbr;
  Gnum                coarvertnum;                /* Number of current coarse vertex           */
  Gnum                finefronnbr;                /* Number of frontier vertices in fine graph */

  if (finegrafptr->parttax == NULL) {             /* If partition array not yet allocated */
    if ((finegrafptr->parttax = (GraphPart *) memAlloc (finegrafptr->s.vertnbr * sizeof (GraphPart))) == NULL) {
      errorPrint ("vgraphSeparateMlUncoarsen: out of memory");
      return     (1);                             /* Allocated data will be freed along with graph structure */
    }
    finegrafptr->parttax -= finegrafptr->s.baseval;
  }

  if (coargrafptr != NULL) {                      /* If coarser graph provided */
    GraphPart * restrict  fineparttax;
    Gnum                  finesize1;              /* Number of vertices in fine part 1 */

    const GraphPart * restrict const  coarparttab = coargrafptr->parttax + coargrafptr->s.baseval;
    Gnum * restrict const             finefrontab = finegrafptr->frontab;

    finesize1   = coargrafptr->compsize[1];       /* Pre-allocate size */
    fineparttax = finegrafptr->parttax;
    for (coarvertnum = finefronnbr = 0, coarvertnbr = coargrafptr->s.vertnbr;
         coarvertnum < coarvertnbr; coarvertnum ++) {
      Gnum                finevertnum0;           /* First multinode vertex          */
      Gnum                finevertnum1;           /* Second multinode vertex         */
      GraphPart           coarpartval;            /* Value of current multinode part */

      finevertnum0 = coarmulttab[coarvertnum].vertnum[0];
      finevertnum1 = coarmulttab[coarvertnum].vertnum[1];
      coarpartval  = coarparttab[coarvertnum];

      fineparttax[finevertnum0] = coarpartval;
      if (coarpartval != 2) {                     /* If vertex is not in separator */
        if (finevertnum0 != finevertnum1) {
          fineparttax[finevertnum1] = coarpartval;
          finesize1 += (Gnum) coarpartval;        /* One extra vertex accounted for in part 1 if (coarpartval == 1) */
        }
      }
      else {                                      /* Vertex is in separator */
        finefrontab[finefronnbr ++] = finevertnum0;
        if (finevertnum0 != finevertnum1) {
          fineparttax[finevertnum1] = coarpartval;
          finefrontab[finefronnbr ++] = finevertnum1; /* One extra vertex in separator */
        }
      }
    }

    finegrafptr->fronnbr     = finefronnbr;
    finegrafptr->compload[0] = coargrafptr->compload[0];
    finegrafptr->compload[1] = coargrafptr->compload[1];
    finegrafptr->compload[2] = coargrafptr->compload[2];
    finegrafptr->comploaddlt = coargrafptr->comploaddlt;
    finegrafptr->compsize[0] = finegrafptr->s.vertnbr - finefronnbr - finesize1;
    finegrafptr->compsize[1] = finesize1;
  }
  else                                            /* No coarse graph provided      */
    vgraphZero (finegrafptr);                     /* Assign all vertices to part 0 */

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (finegrafptr) != 0) {
    errorPrint ("vgraphSeparateMlUncoarsen: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  return (0);
}

/* This routine recursively performs the
** separation recursion.
** It returns:
** - 0   : if separator could be computed.
** - !0  : on error.
*/

static
int
vgraphSeparateMl2 (
Vgraph * restrict const             grafptr,      /* Vertex-separation graph */
const VgraphSeparateMlParam * const paraptr)      /* Method parameters       */
{
  Vgraph                        coargrafdat;
  GraphCoarsenMulti * restrict  coarmulttab;
  int                           o;

  if (vgraphSeparateMlCoarsen (grafptr, &coargrafdat, &coarmulttab, paraptr) == 0) {
    if (((o = vgraphSeparateMl2         (&coargrafdat, paraptr))              == 0) &&
        ((o = vgraphSeparateMlUncoarsen (grafptr, &coargrafdat, coarmulttab)) == 0) &&
        ((o = vgraphSeparateSt          (grafptr, paraptr->stratasc))         != 0)) /* Apply ascending strategy */
      errorPrint ("vgraphSeparateMl2: cannot apply ascending strategy");
    coargrafdat.frontab = NULL;                   /* Prevent frontab of fine graph from being freed */
    vgraphExit (&coargrafdat);
  }
  else {                                          /* Cannot coarsen due to lack of memory or error */
    if (((o = vgraphSeparateMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph          */
        ((o = vgraphSeparateSt          (grafptr, paraptr->stratlow)) != 0)) /* Apply low strategy */
      errorPrint ("vgraphSeparateMl2: cannot apply low strategy");
  }

  return (o);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the muti-level separation.
** It returns:
** - 0 : if separator could be computed.
** - 1 : on error.
*/

int
vgraphSeparateMl (
Vgraph * const                      grafptr,      /*+ Vertex-separation graph +*/
const VgraphSeparateMlParam * const paraptr)      /*+ Method parameters       +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = grafptr->levlnum;                     /* Save graph level               */
  grafptr->levlnum = 0;                           /* Initialize coarsening level    */
  o = vgraphSeparateMl2 (grafptr, paraptr);       /* Perform multi-level separation */
  grafptr->levlnum = levlnum;                     /* Restore graph level            */

  return (o);
}
