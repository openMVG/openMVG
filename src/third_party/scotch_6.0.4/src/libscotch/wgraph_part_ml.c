/* Copyright 2007-2011,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_ml.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module conducts multilevel framew- **/
/**                ork for the vertex overlapped graph pa- **/
/**                rtitioning.                             **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to     27 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "graph_coarsen.h"
#include "wgraph.h"
#include "wgraph_part_ml.h"
#include "wgraph_part_st.h" 

/*
**  The static variables.
*/

static const Gnum           wgraphpartmlloadone = 1;

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
wgraphPartMlCoarsen (
const Wgraph * restrict const         finegrafptr, /*+ Finer graph                                  +*/
Wgraph * restrict const               coargrafptr, /*+ Coarser graph to build                       +*/
GraphCoarsenMulti * restrict * const  coarmultptr, /*+ Pointer to un-based multinode table to build +*/
const WgraphPartMlParam * const       paraptr)     /*+ Method parameters                            +*/
{
  *coarmultptr = NULL;                            /* Allocate coarmulttab along with coarse graph */
  if (graphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr,
                    (paraptr->coarnbr * finegrafptr->partnbr),
                    paraptr->coarval, NULL, NULL, 0, NULL) != 0)
    return (1);                                   /* Return if coarsening failed */

  coargrafptr->parttax  = NULL;                   /* Do not allocate partition data yet */
  coargrafptr->compload = NULL;
  coargrafptr->partnbr  = finegrafptr->partnbr;
  coargrafptr->levlnum  = finegrafptr->levlnum + 1; /* Graph level is coarsening level */

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
wgraphPartMlUncoarsen (
Wgraph * restrict const                   finegrafptr, /*+ Finer graph              +*/
const Wgraph * restrict const             coargrafptr, /*+ Coarser graph            +*/
const GraphCoarsenMulti * restrict const  coarmulttab) /*+ Un-based multinode array +*/
{
  Gnum                          coarvertnbr;
  Gnum                          coarvertnum;      /* Number of current coarse vertex           */
  const Anum * restrict         coarparttab;
  Anum * restrict               fineparttax;
  Gnum * restrict               finefrontab;
  Gnum                          finefronnbr;      /* Number of frontier vertices in fine graph */
  Gnum                          finevertnum;
  WgraphPartList * restrict     finelisttab;
  Gnum                          finevelomsk;
  const Gnum * restrict         finevelobax;      /* Data for handling of optional arrays      */
  Gnum * restrict               finecompload;
  Gnum * restrict               finecompsize;

  const Gnum * restrict const fineverttax = finegrafptr->s.verttax;
  const Gnum * restrict const finevendtax = finegrafptr->s.vendtax;
  const Gnum * restrict const fineedgetax = finegrafptr->s.edgetax;

  if ((finegrafptr->levlnum > 0) &&               /* If partition data not yet allocated */
      (wgraphAlloc (finegrafptr) != 0)) {         /* Allocate tables before processing   */
    errorPrint ("wgraphPartMlUncoarsen: out of memory (1)");
    return     (1);
  }

  if (coargrafptr == NULL) {                      /* If no coarse graph provided */
    wgraphZero (finegrafptr);
    return     (0);
  }

  finecompload = finegrafptr->compload;
  finecompsize = finegrafptr->compsize;

  if ((finelisttab = (WgraphPartList *) memAlloc ((finegrafptr->partnbr + 1) * sizeof (WgraphPartList))) == NULL) { /* TRICK: "+1" to create slot for a "-1" index */
    errorPrint ("wgraphPartMlUncoarsen: out of memory (2)");
    return     (1);
  }
  finelisttab ++;                                 /* TRICK: Trim array so that finelisttab[-1] is valid */
  memSet (finelisttab, ~0, finegrafptr->partnbr * sizeof (WgraphPartList)); /* Set vertex indices to ~0 */

  memSet (finecompload, 0, finegrafptr->partnbr * sizeof (Gnum)); /* Reset load arrays to 0 */
  memSet (finecompsize, 0, finegrafptr->partnbr * sizeof (Gnum));

  if (finegrafptr->s.velotax == NULL) {           /* Set accesses to optional arrays */
    finevelobax = &wgraphpartmlloadone;           /* In case vertices not weighted   */
    finevelomsk = 0;
  }
  else {
    finevelobax = finegrafptr->s.velotax;
    finevelomsk = ~((Gnum) 0);
  }

  finefronnbr = 0;
  finefrontab = finegrafptr->frontab;
  fineparttax = finegrafptr->parttax;
  coarparttab = coargrafptr->parttax + coargrafptr->s.baseval;
  for (coarvertnum = 0, coarvertnbr = coargrafptr->s.vertnbr;
       coarvertnum < coarvertnbr; coarvertnum ++) {
    Anum           coarpartval;                   /* Value of current multinode part */
    Gnum           finevertnum0;
    Gnum           finevertnum1;

    coarpartval  = coarparttab[coarvertnum];
    finevertnum0 = coarmulttab[coarvertnum].vertnum[0];
    finevertnum1 = coarmulttab[coarvertnum].vertnum[1];

    fineparttax[finevertnum0] = coarpartval;
    if (coarpartval >= 0) {                       /* If vertex is not in separator */
      if (finevertnum0 != finevertnum1)
        fineparttax[finevertnum1] = coarpartval;
    }
    else {                                        /* Vertex is in separator */
      finefrontab[finefronnbr ++] = finevertnum0;
      if (finevertnum0 != finevertnum1) {
        fineparttax[finevertnum1] = coarpartval;
        finefrontab[finefronnbr ++] = finevertnum1; /* One extra vertex in separator */
      }
    }
  }
  finegrafptr->fronnbr  = finefronnbr;
  finegrafptr->fronload = coargrafptr->fronload;

  for (finevertnum = finegrafptr->s.baseval; finevertnum < finegrafptr->s.vertnnd; finevertnum ++) {
    Anum                finepartval;

    finepartval = fineparttax[finevertnum];
    if (finepartval >= 0) {
      finecompload[finepartval] += finevelobax[finevertnum & finevelomsk];
      finecompsize[finepartval] ++;
    }
    else {                                        /* Fine vertex is in separator  */
      Gnum                finelistidx;            /* Index of first neighbor part */
      Gnum                fineedgenum;
      Gnum                fineveloval;

      finelistidx = -1;                           /* No neighboring parts recorded yet          */
      finelisttab[-1].vertnum = finevertnum;      /* Separator neighbors will not be considered */
      for (fineedgenum = fineverttax[finevertnum];
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) { /* Compute gain */
        Gnum                finevertend;
        Anum                finepartend;

        finevertend = fineedgetax[fineedgenum];
        finepartend = fineparttax[finevertend];
        if (finelisttab[finepartend].vertnum != finevertnum) { /* If part not yet considered */
          finelisttab[finepartend].vertnum = finevertnum; /* Link it in list of neighbors    */
          finelisttab[finepartend].nextidx = finelistidx;
          finelistidx = finepartend;
        }
      }

      fineveloval = finevelobax[finevertnum & finevelomsk];

      while (finelistidx != -1) {                 /* For all neighboring parts found      */
        finecompload[finelistidx] += fineveloval; /* Add load of separator vertex to part */
        finecompsize[finelistidx] ++;
        finelistidx = finelisttab[finelistidx].nextidx;
      }
    }
  }

  memFree (finelisttab - 1);                      /* TRICK: free array using its real beginning */

#ifdef SCOTCH_DEBUG_WGRAPH2
  if (wgraphCheck (finegrafptr) != 0) {
    errorPrint ("wgraphPartMlUncoarsen: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_WGRAPH2 */

  return (0);
}

/* This routine recursively performs the partitioning.
** It returns:
** - 0   : if separator could be computed.
** - !0  : on error.
*/

static
int
wgraphPartMl2 (
Wgraph * restrict const         grafptr,
const WgraphPartMlParam * const paraptr)
{
  Wgraph                        coargrafdat;
  GraphCoarsenMulti * restrict  coarmulttab;
  int                           o;

  if (wgraphPartMlCoarsen (grafptr, &coargrafdat, &coarmulttab, paraptr) == 0) {
    if (((o = wgraphPartMl2         (&coargrafdat, paraptr)) == 0)              &&
        ((o = wgraphPartMlUncoarsen (grafptr, &coargrafdat, coarmulttab)) == 0) &&
	((o = wgraphPartSt          (grafptr, paraptr->stratasc)) != 0)) /* Apply ascending strategy */
      errorPrint ("wgraphPartMl2: cannot apply ascending strategy");
    wgraphExit (&coargrafdat);
  }
  else {                                          /* Cannot coarsen due to lack of memory or error */
    if (((o = wgraphPartMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph              */
        ((o = wgraphPartSt          (grafptr, paraptr->stratlow)) != 0)) /* Apply low strategy     */
      errorPrint ("wgraphPartMl2: cannot apply low strategy");
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
wgraphPartMl (
Wgraph * const                  wgrafptr,         /*+ Vertex-separation graph +*/
const WgraphPartMlParam * const paraptr)          /*+ Method parameters       +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = wgrafptr->levlnum;                    /* Save graph level               */
  wgrafptr->levlnum = 0;                          /* Initialize coarsening level    */
  o = wgraphPartMl2 (wgrafptr, paraptr);          /* Perform multi-level separation */
  wgrafptr->levlnum = levlnum;                    /* Restore graph level            */

  return (o);
}
