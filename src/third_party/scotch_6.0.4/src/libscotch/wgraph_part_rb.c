/* Copyright 2010,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_rb.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs the vertex overla- **/
/**                pped graph partitioning based on recur- **/
/**                sive bipartitioning approach.           **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 16 mar 2010     **/
/**                                 to     12 aug 2014     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  kgraph_map_rb_part.c for the vertex   **/
/**                  overlapped graph partitioning.        **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_RB

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_cmplt.h"
#include "mapping.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "vgraph_separate_zr.h"
#include "wgraph.h"
#include "wgraph_part_rb.h"
#include "scotch.h"

/*
**  The static variables.
*/

static const Gnum           wgraphpartrbloadone = 1;

/********************************************/
/*                                          */
/* This is the entry point for the vertex   */
/* overlapped graph partitioning based on   */
/* recursive bipartitioning approach.       */
/*                                          */
/********************************************/

/* This routine runs recursive 
** bipartitioning approach.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
wgraphPartRb3 (
const Graph * restrict const     orggrafptr,    /* Graph to induce and bipartition         */
const GraphPart * restrict const orgparttax,    /* Part array of original graph            */
const GraphPart                  indpartval,    /* Part of graph to consider               */
const int                        domnnum,       /* Index of domain onto which map the part */
Mapping * restrict const         mappptr)       /* Final mapping                           */
{
  Gnum               vertnum;

  if (orgparttax == NULL) {                       /* If graph is full graph */
#ifdef SCOTCH_DEBUG_WGRAPH2
    if ((orggrafptr->vnumtax != NULL) || (domnnum != 0)) {
      errorPrint ("wgraphPartRb3: internal error");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_WGRAPH2 */
    memSet (mappptr->parttax + mappptr->grafptr->baseval, 0, orggrafptr->vertnbr * sizeof (ArchDomNum));
  }
  else {                                          /* Graph to consider is a subgraph of the original graph       */
    if (orggrafptr->vnumtax == NULL) {            /* If original graph is not itself a subgraph                  */
      for (vertnum = orggrafptr->baseval; vertnum < orggrafptr->vertnnd; vertnum ++) { /* For all graph vertices */
        if (orgparttax[vertnum] == indpartval)    /* If vertex belongs to the right part                         */
          mappptr->parttax[vertnum] = domnnum;
      }
    }
    else {
      for (vertnum = orggrafptr->baseval; vertnum < orggrafptr->vertnnd; vertnum ++) { /* For all graph vertices */
        if (orgparttax[vertnum] == indpartval)    /* If vertex belongs to the right part                         */
          mappptr->parttax[orggrafptr->vnumtax[vertnum]] = domnnum;
      }
    }
  }

  return (0);
}

static
int
wgraphPartRb2 (
WgraphPartRbData * restrict const dataptr,        /* Top-level graph and partition data       */
Graph * restrict const            orggrafptr,     /* Graph to induce and bipartition          */
const GraphPart * restrict const  orgparttax,     /* Part array of original graph to consider */
const GraphPart                   indpartval,     /* Part of graph to consider                */
const int                         indvertnbr,     /* Number of vertices in part or in graph   */
const int                         domnnum)        /* Index of domain onto which map the part  */
{
  Graph                  indgrafdat;
  Graph *                indgrafptr;
  Vgraph                 actgrafdat;
  Anum                   domnsubidx;
  Anum                   domnsubdlt;
  ArchDom                domnsubtab[2];           /* Target subdomains              */
  Anum                   domnsubnum[2];           /* Index of subdomains in mapping */
  Gnum                   grafsubsiz[2];
  Gnum                   vertnum;
  Mapping * restrict     mappptr;
  int                    i;
  int                    o;

  mappptr = &dataptr->mappdat;
  o = archDomBipart (mappptr->archptr, &mappptr->domntab[domnnum], &domnsubtab[0], &domnsubtab[1]);

  switch (o) {
    case 1 :                                      /* If target domain is terminal */
      return (wgraphPartRb3 (orggrafptr, orgparttax, indpartval, domnnum, mappptr)); /* Update mapping and return */
    case 2 :                                      /* On error */
      errorPrint ("wgraphPartRb2: cannot bipartition domain");
      return     (1);
  }

  indgrafptr = orggrafptr;                        /* Assume we will work on the original graph */
  if (orgparttax != NULL) {                       /* If not the case, build induced subgraph   */
    indgrafptr = &indgrafdat;
    if (graphInducePart (orggrafptr, orgparttax, indvertnbr, indpartval, &indgrafdat) != 0) {
      errorPrint ("wgraphPartRb2: cannot induce graph");
      return     (1);
    }
  }

  actgrafdat.s = *indgrafptr;
  actgrafdat.s.vlbltax = NULL;
  if ((actgrafdat.frontab = (Gnum *) memAlloc (actgrafdat.s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("wgraphPartRb2: out of memory (1)");
    return     (1);
  }
  if ((actgrafdat.parttax = (GraphPart *) memAlloc (actgrafdat.s.vertnbr * sizeof (GraphPart))) == NULL) {
    errorPrint ("wgraphPartRb2: out of memory (2)");
    memFree    (actgrafdat.frontab);
    return     (1);
  }
  actgrafdat.parttax -= actgrafdat.s.baseval;
  vgraphZero (&actgrafdat);                       /* Create active graph */
  if (vgraphSeparateSt (&actgrafdat, dataptr->stratptr) != 0) { /* Perform bipartitioning */
    errorPrint ("wgraphPartRb2: cannot bipartition graph");
    vgraphExit (&actgrafdat);
    return     (1);
  }

  if (actgrafdat.s.vnumtax == NULL) {             /* If the active graph is not itself a subgraph                  */
    for (vertnum = actgrafdat.s.baseval; vertnum < actgrafdat.s.vertnnd; vertnum ++) { /* For all graph vertices */
      if (actgrafdat.parttax[vertnum] == 2) {     /* If vertex belongs to frontier */
	mappptr->parttax[vertnum]   = -1;
	actgrafdat.parttax[vertnum] = 3;
      }
    }
  }
  else {
    for (vertnum = actgrafdat.s.baseval; vertnum < actgrafdat.s.vertnnd; vertnum ++) { /* For all graph vertices */
      if (actgrafdat.parttax[vertnum] == 2) {     /* If vertex belongs to frontier */
	mappptr->parttax[actgrafdat.s.vnumtax[vertnum]]= -1;
	actgrafdat.parttax[vertnum] = 3;
      }
    }
  }


  domnsubdlt = mappptr->domnnbr - domnnum;        /* Increment in domain number */
  domnsubidx = domnnum - domnsubdlt;              /* Place where to insert subdomain */
  mappptr->domnnbr --;                            /* One less subdomain as for now   */
  grafsubsiz[0] = actgrafdat.compsize[0];
  grafsubsiz[1] = actgrafdat.compsize[1];

  o = 0;
  for (i = 1; i >= 0; i --) {                     /* For all subparts             */
    if (grafsubsiz[i] <= 0)                       /* If subpart is empty, skip it */
      continue;
    mappptr->domnnbr ++;                          /* One more subdomain to account for */
    domnsubidx   += domnsubdlt;                   /* Compute location of subdomain */
    domnsubnum[i] = domnsubidx;                   /* Record it before recursion    */
    mappptr->domntab[domnsubidx] = domnsubtab[i]; /* Write it at this place        */
  }

  if (o == 0) {
    for (i = 1; i >= 0; i --) {                   /* For all subparts             */
      if (grafsubsiz[i] <= 0)                     /* If subpart is empty, skip it */
        continue;

      if ((o = wgraphPartRb2 (dataptr, indgrafptr, actgrafdat.parttax, (GraphPart) i, grafsubsiz[i], domnsubnum[i])) != 0)
        return (1);                               /* If problem in recursion, stop */
    }
  }

  memFree (actgrafdat.frontab);                   /* Frontier array of bipartitioning graph is no longer necessary      */
  memFree (actgrafdat.parttax + actgrafdat.s.baseval); /* Frontier array of bipartitioning graph is no longer necessary */
  if (indgrafptr == &indgrafdat)                  /* If an induced subgraph had been created                            */
    graphExit (indgrafptr);                       /* Free it                                                            */

  return (o);
}

int
wgraphPartRb (
Wgraph * restrict const                   grafptr,
const WgraphPartRbParam * restrict const  paraptr)
{
  const Anum * restrict         parttax;
  Anum				partval;
  Gnum                          vertnum;
  Gnum                          velomsk;
  const Gnum * restrict         velobax;              /* Data for handling of optional arrays  */
  Gnum * restrict               frontab;
  Gnum                          fronnbr;
  Gnum                          fronload;
  Gnum * restrict               compload;
  Gnum * restrict               compsize;
  WgraphPartRbData              datadat;
  Arch                          archdat;
  WgraphPartList * restrict     listtab;

  const Gnum * restrict const   verttax = grafptr->s.verttax;
  const Gnum * restrict const   vendtax = grafptr->s.vendtax;
  const Gnum * restrict const   edgetax = grafptr->s.edgetax;

  if ((listtab = (WgraphPartList *) memAlloc ((grafptr->partnbr + 1) * sizeof (WgraphPartList))) == NULL) { /* TRICK: "+1" to create slot for a "-1" index */
    errorPrint ("wgraphPartRb: out of memory (1)");
    return     (1);
  }
  listtab ++;                                     /* TRICK: Trim array so that listtab[-1] is valid */
  memSet (listtab, ~0, grafptr->partnbr * sizeof (WgraphPartList)); /* Set vertex indices to ~0     */
  
  datadat.grafptr  = &grafptr->s;
  datadat.frontab  = grafptr->frontab;            /* Re-use frontier array */
  datadat.fronnbr  = 0;
  datadat.stratptr = paraptr->stratptr;
  datadat.mappdat.grafptr = &grafptr->s;
  datadat.mappdat.parttax = grafptr->parttax;     /* Re-use part array */
  datadat.mappdat.domnmax = grafptr->partnbr + 1;
  datadat.mappdat.domnnbr = 1;

  SCOTCH_archCmplt ((SCOTCH_Arch *) &archdat, grafptr->partnbr); /* Create a complete graph architecture */
  datadat.mappdat.archptr = &archdat;

  archDomFrst (datadat.mappdat.archptr, &datadat.mappdat.domnorg); /* Get first domain of architecture */
  if ((datadat.mappdat.domntab = (ArchDom *) memAlloc ((grafptr->partnbr + 2) * sizeof (ArchDom))) == NULL) {
    errorPrint ("wgraphPartRb: out of memory (2)");
    memFree    (listtab - 1);                     /* TRICK: free array using its real beginning */
    return     (1);
  }
  datadat.mappdat.domntab[0] = datadat.mappdat.domnorg; /* Set first domain */

  if (wgraphPartRb2 (&datadat, &grafptr->s, NULL, 0, grafptr->s.vertnbr, 0) != 0) {
    errorPrint ("wgraphPartRb: internal error (1)");
    return     (1);
  }

  if (grafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &wgraphpartrbloadone;               /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  compload = grafptr->compload;
  compsize = grafptr->compsize;
  memSet (compload, 0, grafptr->partnbr * sizeof (Gnum));
  memSet (compsize, 0, grafptr->partnbr * sizeof (Gnum));

  parttax  = grafptr->parttax;
  frontab  = grafptr->frontab;
  fronnbr  =
  fronload = 0;
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    Gnum                partval;

    partval = parttax[vertnum];
    if (partval >= 0) {
      compload[partval] += velobax[vertnum & velomsk];
      compsize[partval] ++;
    }
    else {                                        /* Vertex is in separator       */
      Gnum                listidx;                /* Index of first neighbor part */
      Gnum                edgenum;
      Gnum                veloval;

      frontab[fronnbr ++] = vertnum;              /* Add vertex to frontier */
      fronload           += velobax[vertnum & velomsk];

      listidx = -1;                               /* No neighboring parts recorded yet          */
      listtab[-1].vertnum = vertnum;              /* Separator neighbors will not be considered */
      for (edgenum = verttax[vertnum];
           edgenum < vendtax[vertnum]; edgenum ++) { /* Compute gain */
        Gnum                vertend;
        Gnum                partend;

        vertend = edgetax[edgenum];
        partend = parttax[vertend];
        if (listtab[partend].vertnum != vertnum) { /* If part not yet considered  */
          listtab[partend].vertnum = vertnum;     /* Link it in list of neighbors */
          listtab[partend].nextidx = listidx;
          listidx = partend;
        }
      }

      veloval = velobax[vertnum & velomsk];

      while (listidx != -1) {                     /* For all neighboring parts found      */
        compload[listidx] += veloval;             /* Add load of separator vertex to part */
        compsize[listidx] ++;
        listidx = listtab[listidx].nextidx;
      }
    }
  }
  grafptr->fronnbr  = fronnbr;
  grafptr->fronload = fronload;

#if 0 /* TODO REMOVE */
  for (partval = 0; partval < grafptr->partnbr; partval ++)
    printf("\033[0;33mcompload[%d] %d %d\033[0m\n", partval, grafptr->compload[partval], grafptr->compsize[partval]);
#endif

  memFree (datadat.mappdat.domntab);              /* Free only newly allocated array of mapping */
  memFree (listtab - 1);                          /* TRICK: free array using its real beginning */

#ifdef SCOTCH_DEBUG_WGRAPH2
  if (wgraphCheck (grafptr) != 0) {
    errorPrint ("wgraphPartRb: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_WGRAPH2 */

  return (0);
}
