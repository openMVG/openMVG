/* Copyright 2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb_part.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm for    **/
/**                (eventually weighted) complete graph    **/
/**                target architectures.                   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 sep 2008     **/
/**                                 to     31 aug 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     16 sep 2014     **/
/**                                                        **/
/**   NOTES      : # This is a rewrite of kgraphMapRb()    **/
/**                  for complete-graph target topologies. **/
/**                  Its advantage over kgraphMapRbMap()   **/
/**                  is that no job arrays are allocated,  **/
/**                  which can save space for instance for **/
/**                  using the variable-sized complete     **/
/**                  graph architecture.                   **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB_PART

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
#include "kgraph_map_rb_part.h"

/********************************************/
/*                                          */
/* This is the entry point for the Dual     */
/* Recursive Bipartitioning mapping method. */
/*                                          */
/********************************************/

/* This routine updates partial mappings
** according to the result of the recursive
** bipartitioning process.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
void
kgraphMapRbPart3 (
const Graph * restrict const      srcgrafptr,     /* Graph to induce and bipartition         */
const GraphPart * restrict const  srcparttax,     /* Part array of original graph            */
const GraphPart                   indpartval,     /* Part of graph to consider               */
const int                         domnnum,        /* Index of domain onto which map the part */
Mapping * restrict const          mappptr)        /* Final mapping                           */
{
  Gnum               vertnum;

  const Gnum * restrict const srcvnumtax = srcgrafptr->vnumtax;
  Anum * restrict const       mapparttax = mappptr->parttax;

  if (srcparttax == NULL) {                       /* If graph is full graph */
#ifdef SCOTCH_DEBUG_KGRAPH2
    if (domnnum != 0) {
      errorPrint ("kgraphMapRbPart3: internal error (1)");
      return;
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    if (srcvnumtax == NULL)                       /* If full graph doesn't have fixed vertices */
      memSet (mapparttax + srcgrafptr->baseval, 0, srcgrafptr->vertnbr * sizeof (Anum));
    else {
      Gnum                vertnnd;

      for (vertnum = srcgrafptr->baseval, vertnnd = srcgrafptr->vertnnd;
           vertnum < vertnnd; vertnum ++) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (mapparttax[srcvnumtax[vertnum]] == ~0) {
          errorPrint ("kgraphMapRbPart3: internal error (2)");
          return;
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        mapparttax[srcvnumtax[vertnum]] = domnnum;
      }
    }
  }
  else {                                          /* Graph to consider is a subgraph of the original graph */
    if (srcvnumtax == NULL) {                     /* If original graph is not itself a subgraph            */
      Gnum                vertnnd;

      for (vertnum = srcgrafptr->baseval, vertnnd = srcgrafptr->vertnnd;
           vertnum < vertnnd; vertnum ++) {
        if (srcparttax[vertnum] == indpartval) {  /* If vertex belongs to the right part */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (mapparttax[vertnum] == ~0) {
            errorPrint ("kgraphMapRbPart3: internal error (3)");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          mapparttax[vertnum] = domnnum;
        }
      }
    }
    else {
      Gnum                vertnnd;

      for (vertnum = srcgrafptr->baseval, vertnnd = srcgrafptr->vertnnd;
           vertnum < vertnnd; vertnum ++) {
        if (srcparttax[vertnum] == indpartval) {  /* If vertex belongs to the right part */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (mapparttax[srcvnumtax[vertnum]] == ~0) {
            errorPrint ("kgraphMapRbPart3: internal error (4)");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          mapparttax[srcvnumtax[vertnum]] = domnnum;
        }
      }
    }
  }
}

/* This routine is the core of the degenerated,
** graph partitioning version of the Dual Recursive
** Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
kgraphMapRbPart2 (
const KgraphMapRbData * restrict const  dataptr,  /*+ Global mapping data                        +*/
const Graph * restrict const            srcgrafptr, /*+ Graph to induce and bipartition          +*/
const GraphPart * restrict const        srcparttax, /*+ Part array of original graph to consider +*/
const GraphPart                         indpartval, /*+ Part of graph to consider                +*/
const Gnum                              indvertnbr, /*+ Number of vertices in part or in graph   +*/
const Anum                              domnnum,  /*+ Index of domain onto which to map the part +*/
const Anum                              vflonbr,  /*+ Number of fixed vertex load slots          +*/
KgraphMapRbVflo * restrict const        vflotab)  /*+ Array of fixed vertex load slots           +*/
{
  Graph               indgrafdat;
  const Graph *       indgrafptr;
  Bgraph              actgrafdat;
  Anum                domnsubidx;
  Anum                domnsubdlt;
  ArchDom             domnsubtab[2];              /* Target subdomains                           */
  Anum                domnidxtab[2];              /* Index of subdomains in mapping              */
  Gnum                vertnbrtab[2];              /* Number of vertices in subgraphs             */
  Anum                vflonbrtab[2];              /* Number of fixed vertex slots in subdomains  */
  Gnum                vflowgttab[2];              /* Weights of fixed vertex slots in subdomains */
  Mapping * restrict  mappptr;
  int                 avarval;                    /* Flag set if variable-sized                  */
  int                 i;
  int                 o;

  mappptr = dataptr->mappptr;
  avarval = archVar (mappptr->archptr);
  o = (avarval &&                                 /* If architecture is variable-sized   */
       (indvertnbr <= 1))                         /* And source subgraph of minimal size */
      ? 1                                         /* Then do not bipartition target more */
      : archDomBipart (mappptr->archptr, &mappptr->domntab[domnnum], &domnsubtab[0], &domnsubtab[1]);

  switch (o) {
    case 1 :                                      /* If target domain is terminal                           */
      kgraphMapRbPart3 (srcgrafptr, srcparttax, indpartval, domnnum, mappptr); /* Update mapping and return */
      return (0);
    case 2 :                                      /* On error */
      errorPrint ("kgraphMapRbPart2: cannot bipartition domain");
      return     (1);
  }

  indgrafptr = srcgrafptr;                        /* Assume we will work on the original graph */
  if ((srcparttax != NULL) &&                     /* If not the case, build induced subgraph   */
      (indvertnbr < srcgrafptr->vertnbr)) {
    indgrafptr = &indgrafdat;
    if (graphInducePart (srcgrafptr, srcparttax, indvertnbr, indpartval, &indgrafdat) != 0) {
      errorPrint ("kgraphMapRbPart2: cannot induce graph");
      return     (1);
    }
  }

  kgraphMapRbVfloSplit (mappptr->archptr, domnsubtab, vflonbr, vflotab, vflonbrtab, vflowgttab);

  if (kgraphMapRbBgraph (dataptr, &actgrafdat, indgrafptr, mappptr, domnsubtab, vflowgttab) != 0) { /* Create active graph */
    errorPrint ("kgraphMapRbPart2: cannot create bipartition graph");
    return     (1);
  }

  if (! avarval) {                                /* If not variable-sized, impose constraints on bipartition */
    double              comploadavg;

    comploadavg = (double) (actgrafdat.s.velosum + vflowgttab[0] + vflowgttab[1]) /
                  (double) archDomWght (mappptr->archptr, &mappptr->domntab[domnnum]);
    actgrafdat.compload0min = actgrafdat.compload0avg -
                              (Gnum) MIN ((dataptr->comploadmax - comploadavg) * (double) actgrafdat.domnwght[0],
                                          (comploadavg - dataptr->comploadmin) * (double) actgrafdat.domnwght[1]);
    actgrafdat.compload0max = actgrafdat.compload0avg +
                              (Gnum) MIN ((comploadavg - dataptr->comploadmin) * (double) actgrafdat.domnwght[0],
                                          (dataptr->comploadmax - comploadavg) * (double) actgrafdat.domnwght[1]);
  }

  if (bgraphBipartSt (&actgrafdat, dataptr->paraptr->strat) != 0) { /* Perform bipartitioning */
    errorPrint ("kgraphMapRbPart2: cannot bipartition graph");
    bgraphExit (&actgrafdat);
    return     (1);
  }
  memFree (actgrafdat.frontab);                   /* Frontier array of bipartitioning graph is no longer necessary */
  actgrafdat.s.flagval &= ~BGRAPHFREEFRON;

  if (archVar (mappptr->archptr)) {               /* If architecture is variable-sized */
    if ((actgrafdat.compsize0 == 0) ||            /* If bipartition failed             */
        (actgrafdat.compsize0 == actgrafdat.s.vertnbr)) {
      bgraphExit (&actgrafdat);                   /* Free bipartition graph (that is, parttax)                        */
      if (indgrafptr == &indgrafdat)              /* If an induced subgraph had been created                          */
        graphExit (&indgrafdat);                  /* Free it                                                          */
      kgraphMapRbPart3 (srcgrafptr, srcparttax, indpartval, domnnum, mappptr); /* Update mapping with original domain */
      return (0);
    }
  }

  domnsubdlt = mappptr->domnnbr - domnnum;        /* Increment in domain number      */
  domnsubidx = domnnum - domnsubdlt;              /* Place where to insert subdomain */
  mappptr->domnnbr --;                            /* One less subdomain as for now   */
  vertnbrtab[0] = actgrafdat.compsize0;
  vertnbrtab[1] = actgrafdat.s.vertnbr - actgrafdat.compsize0;

  o = 0;
  for (i = 1; i >= 0; i --) {                     /* For all subparts             */
    if (vertnbrtab[i] <= 0)                       /* If subpart is empty, skip it */
      continue;

    mappptr->domnnbr ++;                          /* One more subdomain to account for */
    if (mappptr->domnnbr > mappptr->domnmax) {
      if ((o = mapResize (mappptr, mappptr->domnmax + (mappptr->domnmax >> 2) + 8)) != 0) { /* Increase size by 25% */
        errorPrint ("kgraphMapRbPart: cannot resize structures");
        break;
      }
    }
    domnsubidx   += domnsubdlt;                   /* Compute location of subdomain */
    domnidxtab[i] = domnsubidx;                   /* Record it before recursion    */
    mappptr->domntab[domnsubidx] = domnsubtab[i]; /* Write it at this (new) place  */
  }

  if (o == 0) {
    for (i = 1; i >= 0; i --) {                   /* For all subparts             */
      if (vertnbrtab[i] <= 0)                     /* If subpart is empty, skip it */
        continue;

      if ((o = kgraphMapRbPart2 (dataptr, indgrafptr, actgrafdat.parttax, (GraphPart) i, vertnbrtab[i],
                                 domnidxtab[i], vflonbrtab[i], vflotab + (i * vflonbrtab[0]))) != 0)
        return (1);                               /* If problem in recursion, stop */
    }
  }

  bgraphExit (&actgrafdat);                       /* Free bipartition graph (that is, parttax) */
  if (indgrafptr == &indgrafdat)                  /* If an induced subgraph had been created   */
    graphExit (&indgrafdat);                      /* Free it                                   */

  return (o);
}

/* This routine is the entry point for
** the degenerated, graph partitioning
** version, of the Dual Recursive
** Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbPart (
const KgraphMapRbData * restrict const  dataptr,  /*+ Global mapping data                  +*/
const Graph * restrict const            grafptr,  /*+ Graph to map, without fixed vertices +*/
const Anum                              vflonbr,  /*+ Number of fixed vertex load slots    +*/
KgraphMapRbVflo * restrict const        vflotab)  /*+ Array of fixed vertex load slots     +*/
{
  Mapping * restrict const  mappptr = dataptr->mappptr;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (dataptr->pfixtax != NULL) {                 /* In debug mode, fixed vertex parts are set to ~0 */
    Gnum                vertnum;

    for (vertnum = dataptr->grafptr->baseval; vertnum < dataptr->grafptr->vertnnd; vertnum ++)
      mappptr->parttax[vertnum] = (dataptr->pfixtax[vertnum] >= 0) ? ~0 : 0;
  }
  else
    memSet (mappptr->parttax + dataptr->grafptr->baseval, 0, dataptr->grafptr->vertnbr * sizeof (Anum));
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  mappptr->domntab[0] = mappptr->domnorg;         /* Initialize mapping */
  mappptr->domnnbr    = 1;
  return (kgraphMapRbPart2 (dataptr, grafptr, NULL, 0, grafptr->vertnbr, 0, vflonbr, vflotab));
}
