/* Copyright 2004,2007,2008,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_es.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a matrix ordering software.     **/
/**                This module computes the node separator **/
/**                of a graph based on the edge-separation **/
/**                module of "bgraph_bipart_st.c".         **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 17 oct 1996     **/
/**                                 to   : 07 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 18 aug 2004     **/
/**                                 to     20 aug 2004     **/
/**                # Version 5.0  : from : 24 jan 2007     **/
/**                                 to     12 sep 2007     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                                                        **/
/**   NOTES      : # This algorithm comes from:            **/
/**                  "Computing the Block Triangular form  **/
/**                   of a Sparse Matrix", A. Pothen and   **/
/**                  C.-J. Fan, ACM Trans. on Mathematical **/
/**                  Software, 16 (4), pp 303-324, 1990.   **/
/**                  and from:                             **/
/**                  "Implementations of $O(n^{1/2}\tau)$  **/
/**                   assignment algorithms", I. Duff and  **/
/**                  T. Wieberg, ACM Trans. on Math.       **/
/**                  Software, 4, pp 267-287, 1988.        **/
/**                                                        **/
/**                # The choice of the separator to take,  **/
/**                  either HR u SC u VC or HR u SR u VC,  **/
/**                  is made regarding the size of the     **/
/**                  separator only, irrespective of its   **/
/**                  balance. This choice is made because  **/
/**                  else an imbalance ratio should be     **/
/**                  provided for this method, and because **/
/**                  it is assumed that the edge biparti-  **/
/**                  tioning method is assumed to have     **/
/**                  reached suitable balance. When they   **/
/**                  are equal, the choice is biased       **/
/**                  towards SR, because the xC block is   **/
/**                  the one which has less vertices so    **/
/**                  removing more separator vertices from **/
/**                  it would mean increasing imbalance.   **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_ES

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "vgraph.h"
#include "vgraph_separate_es.h"

/*********************************************/
/*                                           */
/* These routines compute a vertex separator */
/* from an edge separator represented as a   */
/* bipartite graph, by minimum covering.     */
/*                                           */
/*********************************************/

/* This routine computes a vertex separator
** from an edge separator represented as a
** bipartite graph, by minimum covering.
** It returns:
** - 0   : if a separator could be computed.
** - !0  : on error.
*/

static
int
vgraphSeparateEsCover (
const Graph * restrict const  grafptr,            /* Bipartite graph to cover         */
const Gnum                    partnbr,            /* Number of vertices in first part */
Gnum * const                  sepatab,            /* Array of covering vertices       */
Gnum * const                  sepaptr)            /* Pointer to size of the array     */
{
  Gnum * restrict                 levltax;        /* Array of vertex level values      */
  Gnum                            levlmax;        /* Maximum level searched            */
  Gnum * restrict                 listtab;        /* List of reachable augmenting rows */
  Gnum                            listnbr;        /* Number of items in list           */
  Gnum * restrict                 matetax;        /* Matching array                    */
  Gnum *                          queutab;        /* Queue of (free) column nodes      */
  Gnum * restrict                 queuhead;       /* Head of queue                     */
  Gnum * restrict                 queutail;       /* Tail of queue                     */
  VgraphSeparateEsTrav * restrict travtax;        /* Array of traversal flag values    */
  VgraphSeparateEsType * restrict typetax;        /* Vertex type in the graph          */
  Gnum                            loadcval;       /* Load of subset (HR u SC u VC)     */
  Gnum                            loadrval;       /* Load of subset (HR u SR u VC)     */
  Gnum                            sizecval;       /* Load of subset (HR u SC u VC)     */
  Gnum                            sizerval;       /* Load of subset (HR u SR u VC)     */
  Gnum                            vertnum;

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (sizeof (VgraphSeparateEsType) > sizeof (VgraphSeparateEsTrav)) {  /* Assert next trick will work */
    errorPrint ("vgraphSeparateEsCover: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
  if (memAllocGroup ((void **) (void *)
                     &travtax, (size_t) (grafptr->vertnbr * sizeof (VgraphSeparateEsTrav)), /* TRICK: VgraphSeparateEsType should also fit */
                     &matetax, (size_t) (grafptr->vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("vgraphSeparateEsCover: out of memory (1)");
    return     (1);
  }
  if (memAllocGroup ((void **) (void *)
                     &queutab, (size_t) (partnbr          * sizeof (Gnum)),
                     &levltax, (size_t) (grafptr->vertnbr * sizeof (Gnum)),
                     &listtab, (size_t) (grafptr->vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("vgraphSeparateEsCover: out of memory (2)");
    memFree    (travtax);                         /* Not based yet */
    return     (1);
  }
  travtax -= grafptr->baseval;
  matetax -= grafptr->baseval;
  levltax -= grafptr->baseval;

  memSet (matetax + (partnbr + grafptr->baseval), ~0, (grafptr->vertnbr - partnbr) * sizeof (Gnum));
  for (vertnum = grafptr->baseval;                /* Compute a cheap matching */
       vertnum < (partnbr + grafptr->baseval); vertnum ++) {
    Gnum              edgenum;
    Gnum              matenum;

    for (edgenum = grafptr->verttax[vertnum], matenum = ~0; /* Search a matching end vertex */
         edgenum < grafptr->vendtax[vertnum]; edgenum ++) {
      Gnum              vertend;

      vertend = grafptr->edgetax[edgenum];
      if (matetax[vertend] == ~0) {               /* If an unmatched end vertex is found */
        matenum          = vertend;
        matetax[vertend] = vertnum;
        break;
      }
    }
    matetax[vertnum] = matenum;
  }

  do {                                            /* Matching augmentation loop */
    queuhead =                                    /* Flush the data structures  */
    queutail = queutab;
    listnbr  = 0;
    memSet (levltax + grafptr->baseval, 0, grafptr->vertnbr * sizeof (Gnum));
    memSet (travtax + grafptr->baseval, 0, grafptr->vertnbr * sizeof (VgraphSeparateEsTrav));
    levlmax = ~0;

    for (vertnum = grafptr->baseval;              /* Enqueue unmatched column nodes */
         vertnum < (partnbr + grafptr->baseval); vertnum ++) {
      if (matetax[vertnum] == ~0) {
        *queuhead ++ = vertnum;
        levltax[vertnum] = 1;
      }
    }

    while (queuhead > queutail) {                 /* As long as there are free columns */
      Gnum              vertcol;

      vertcol = *queutail ++;                     /* Get the free column vertex */
      if (levltax[vertcol] < levlmax) {
        Gnum              edgenum;

        travtax[vertcol] = VGRAPHSEPAESTRAVUSED;  /* Column has been reached */

        for (edgenum = grafptr->verttax[vertcol]; /* For all neighboring rows */
             edgenum < grafptr->vendtax[vertcol]; edgenum ++) {
          Gnum              vertrow;

          vertrow = grafptr->edgetax[edgenum];
          if (travtax[vertrow] == VGRAPHSEPAESTRAVFREE) { /* If row not yet reached yet */
            travtax[vertrow] = VGRAPHSEPAESTRAVUSED; /* Now it is                       */
            if (matetax[vertrow] == ~0) {         /* If row is unmatched                */
              listtab[listnbr ++] = vertrow;      /* Put it in list                     */
              levlmax = levltax[vertcol];         /* Do not go any further              */
            }
            else {                                /* Row is matched              */
              *queuhead ++ = matetax[vertrow];    /* Enqueue its matching column */
              levltax[matetax[vertrow]] = levltax[vertcol] + 1;
            }
          }
        }
      }
    }

    if (listnbr <= 0)                             /* If no free rows could be reached */
      break;                                      /* Then the matching is maximal     */

    while (-- listnbr >= 0)                       /* For all rows in list, try to augment the matching */
      vgraphSeparateEsCoverAugment (levltax, levlmax, matetax, travtax, grafptr->verttax, grafptr->vendtax, grafptr->edgetax, listtab[listnbr]);
  } while (1);

  memFree (queutab);                              /* Free group leader of arrays no longer in use */
  typetax = (VgraphSeparateEsType *) travtax;     /* TRICK: re-use traversal table as type table  */

  for (vertnum = grafptr->baseval; vertnum < (partnbr + grafptr->baseval); vertnum ++) /* Pre-set vertex types */
    typetax[vertnum] = VGRAPHSEPAESTYPESC;
  for ( ; vertnum < grafptr->vertnnd; vertnum ++)
    typetax[vertnum] = VGRAPHSEPAESTYPESR;
  for (vertnum = grafptr->baseval; vertnum < (partnbr + grafptr->baseval); vertnum ++) /* For all column vertices */
    if (matetax[vertnum] == ~0)                   /* If vertex is unmatched */
      vgraphSeparateEsCoverCol (matetax, typetax, grafptr->verttax, grafptr->vendtax, grafptr->edgetax, vertnum); /* Find HC and HR */
  for ( ; vertnum < grafptr->vertnnd; vertnum ++) /* For all row vertices   */
    if (matetax[vertnum] == ~0)                   /* If vertex is unmatched */
      vgraphSeparateEsCoverRow (matetax, typetax, grafptr->verttax, grafptr->vendtax, grafptr->edgetax, vertnum); /* Find VC and VR */

  sizecval =                                      /* Reset sizes */
  sizerval = 0;
  if (grafptr->velotax != NULL) {                 /* If graph vertices are weighted */
    Gnum              vertnum;

    loadcval =                                    /* Reset loads */
    loadrval = 0;
    for (vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++) { /* Accumulate loads */
      VgraphSeparateEsType  typeval;
      Gnum                  veloval;
      Gnum                  bitcval;
      Gnum                  bitrval;

      typeval = typetax[vertnum];
      veloval = grafptr->velotax[vertnum];
      bitcval = (typeval >> VGRAPHSEPAESTYPEBITC) & 1;
      bitrval =  typeval >> VGRAPHSEPAESTYPEBITR; /* TRICK: highest bit so does not need mask  */

      loadcval += bitcval * veloval;              /* Superscalar update */
      loadrval += bitrval * veloval;
      sizecval += bitcval;
      sizerval += bitrval;
    }
  }
  else {                                          /* Graph vertices are not weighted */
    Gnum              vertnum;

    for (vertnum = grafptr->baseval; vertnum < grafptr->vertnnd; vertnum ++) { /* Accumulate vertex sizes */
      sizecval += (typetax[vertnum] >> VGRAPHSEPAESTYPEBITC) & 1; /* Superscalar update                   */
      sizerval +=  typetax[vertnum] >> VGRAPHSEPAESTYPEBITR; /* TRICK: highest bit so does not need mask  */
    }
    loadcval = sizecval;                          /* Loads equal sizes */
    loadrval = sizerval;
  }

  if (loadcval < loadrval) {                      /* If separator with SC is smaller */
    Gnum              vertnum;
    Gnum              sepanum;

    *sepaptr = sizecval;

    for (vertnum = grafptr->baseval, sepanum = 0;
         vertnum < grafptr->vertnnd; vertnum ++) {
      if ((typetax[vertnum] & VGRAPHSEPAESTYPEHRSCVC) != 0) {
#ifdef SCOTCH_DEBUG_VGRAPH2
        if ((sepanum >= sizecval) ||
            ((typetax[vertnum] != VGRAPHSEPAESTYPEHR) &&
             (typetax[vertnum] != VGRAPHSEPAESTYPESC) &&
             (typetax[vertnum] != VGRAPHSEPAESTYPEVC))) {
          errorPrint ("vgraphSeparateEsCover: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
        sepatab[sepanum ++] = vertnum;
      }
    }
#ifdef SCOTCH_DEBUG_VGRAPH2
    if (sepanum != sizecval) {
      errorPrint ("vgraphSeparateEsCover: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
  }
  else {                                          /* If separator with SR is smaller */
    Gnum              vertnum;
    Gnum              sepanum;

    *sepaptr = sizerval;

    for (vertnum = grafptr->baseval, sepanum = 0;
         vertnum < grafptr->vertnnd; vertnum ++) {
      if ((typetax[vertnum] & VGRAPHSEPAESTYPEHRSRVC) != 0) {
#ifdef SCOTCH_DEBUG_VGRAPH2
        if ((sepanum >= sizerval) ||
            ((typetax[vertnum] != VGRAPHSEPAESTYPEHR) &&
             (typetax[vertnum] != VGRAPHSEPAESTYPESR) &&
             (typetax[vertnum] != VGRAPHSEPAESTYPEVC))) {
          errorPrint ("vgraphSeparateEsCover: internal error (4)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
        sepatab[sepanum ++] = vertnum;
      }
    }
#ifdef SCOTCH_DEBUG_VGRAPH2
    if (sepanum != sizerval) {
      errorPrint ("vgraphSeparateEsCover: internal error (5)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
  }

  memFree (travtax + grafptr->baseval);           /* Free group leader of remaining arrays */

  return (0);
}

/* This routine augments the current matching
** by performing a backtracking depth-first
** search from a free row vertex to a free
** column vertex, guided by the level values.
** It returns:
** - 0   : backtracking succeeded.
** - !0  : could not find a valid return path.
*/

static
int
vgraphSeparateEsCoverAugment (
const Gnum * restrict const           levltax,
const Gnum                            levlcur,    /* Current backtracking level */
Gnum * restrict const                 matetax,
VgraphSeparateEsTrav * restrict const travtax,
const Gnum * restrict const           verttax,
const Gnum * restrict const           vendtax,
const Gnum * restrict const           edgetax,
const Gnum                            vertrow)    /* Row vertex to backtrack from */
{
  Gnum                edgenum;

  travtax[vertrow] = VGRAPHSEPAESTRAVDRTY;        /* Never re-use this row */

  for (edgenum = verttax[vertrow]; edgenum < vendtax[vertrow]; edgenum ++) {
    Gnum                vertcol;

    vertcol = edgetax[edgenum];                   /* Get column vertex                                     */
    if ((travtax[vertcol] == VGRAPHSEPAESTRAVUSED) && /* If this column may be a backtracking path         */
        (levltax[vertcol] == levlcur)) {          /* At the proper distance from a free column             */
      travtax[vertcol] = VGRAPHSEPAESTRAVDRTY;    /* Never re-use this column                              */
      if ((levlcur == 1) ||                       /* If we have (recursively) reached a free column vertex */
          (vgraphSeparateEsCoverAugment (levltax, levlcur - 1, matetax, travtax, verttax, vendtax, edgetax, matetax[vertcol]) == 0)) {
        matetax[vertcol] = vertrow;               /* Switch the edges of the augmenting path */
        matetax[vertrow] = vertcol;
        return (0);                               /* Backtracking process is under way */
      }
    }
  }

  return (1);                                     /* No improvement could be done */
}

/* Starting from unmatched column and row vertices,
** these routines perform depth-first traversals of
** the bipartite graph, following alternating paths.
** It is assumed that the matchings are sufficently
** large, so that the depth of the trees is small
** and the stack will not overflow.
** They return:
** - VOID  : in all cases.
*/

static
void
vgraphSeparateEsCoverCol (
const Gnum * restrict const           matetax,
VgraphSeparateEsType * restrict const typetax,
const Gnum * restrict const           verttax,
const Gnum * restrict const           vendtax,
const Gnum * restrict const           edgetax,
const Gnum                            vertcol)    /* Column vertex index */
{
  Gnum              edgenum;

  if (typetax[vertcol] == VGRAPHSEPAESTYPEHC)     /* If vertex already traversed */
    return;

  typetax[vertcol] = VGRAPHSEPAESTYPEHC;

  for (edgenum = verttax[vertcol]; edgenum < vendtax[vertcol]; edgenum ++) {
    Gnum              vertrow;

    vertrow = edgetax[edgenum];
    if (typetax[vertrow] == VGRAPHSEPAESTYPEHR)   /* If end vertex already traversed */
      continue;                                   /* Skip to next vertex             */
    typetax[vertrow] = VGRAPHSEPAESTYPEHR;
    if (matetax[vertrow] != ~0)                   /* If end vertex matched */
      vgraphSeparateEsCoverCol (matetax, typetax, verttax, vendtax, edgetax, matetax[vertrow]);
  }
}

static
void
vgraphSeparateEsCoverRow (
const Gnum * restrict const           matetax,
VgraphSeparateEsType * restrict const typetax,
const Gnum * restrict const           verttax,
const Gnum * restrict const           vendtax,
const Gnum * restrict const           edgetax,
const Gnum                            vertrow)    /* Row vertex index */
{
  Gnum              edgenum;

  if (typetax[vertrow] == VGRAPHSEPAESTYPEVR)     /* If vertex already traversed */
    return;

  typetax[vertrow] = VGRAPHSEPAESTYPEVR;

  for (edgenum = verttax[vertrow]; edgenum < vendtax[vertrow]; edgenum ++) {
    Gnum              vertcol;

    vertcol = edgetax[edgenum];
    if (typetax[vertcol] == VGRAPHSEPAESTYPEVC)   /* If end vertex already traversed */
      continue;                                   /* Skip to next vertex             */
    typetax[vertcol] = VGRAPHSEPAESTYPEVC;
    if (matetax[vertcol] != ~0)                   /* If end vertex matched */
      vgraphSeparateEsCoverRow (matetax, typetax, verttax, vendtax, edgetax, matetax[vertcol]);
  }
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine separates the given graph by first
** computing an edge separator, according to the
** given bipartitioning strategy, and then turning
** it into a vertex separator.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
vgraphSeparateEs (
Vgraph * restrict const             grafptr,      /*+ Active graph      +*/
const VgraphSeparateEsParam * const paraptr)      /*+ Method parameters +*/
{
  Bgraph            actgrafdat;                   /* Active graph  structure   */
  Graph             bipgrafdat;                   /* Bipartite graph structure */

  actgrafdat.s         = grafptr->s;              /* Initialize active graph */
  actgrafdat.s.flagval = grafptr->s.flagval & ~(GRAPHFREETABS | BGRAPHFREEPART | BGRAPHFREEFRON);
  actgrafdat.s.vnumtax = NULL;
  actgrafdat.s.vlbltax = NULL;
  actgrafdat.veextax   = NULL;                    /* No external gains                           */
  actgrafdat.parttax   = grafptr->parttax;        /* Inherit arrays from vertex separation graph */
  actgrafdat.frontab   = grafptr->frontab;
  bgraphInit2 (&actgrafdat, 1, 1, 1, 0, 0);       /* Complete initialization and set all vertices to part 0 */

  if (bgraphBipartSt (&actgrafdat, paraptr->strat) != 0) { /* Bipartition active subgraph */
    errorPrint ("vgraphSeparateEs: cannot bipartition active graph");
    return     (1);
  }

  grafptr->compload[0] = actgrafdat.compload0;    /* Reset vertex counts */
  grafptr->compload[1] = actgrafdat.s.velosum - actgrafdat.compload0;
  grafptr->compsize[0] = actgrafdat.compsize0;
  grafptr->compsize[1] = actgrafdat.s.vertnbr - actgrafdat.compsize0;

  if (actgrafdat.fronnbr > 0) {                   /* If edge separator is not empty     */
    if (paraptr->widtval == VGRAPHSEPAESWIDTHTHIN) { /* If thin vertex separator wanted */
      Gnum * restrict   actvnumtax;
      Gnum              actfronnum;
      Gnum              bipvelosum;
      Gnum              bipedgenbr;               /* Number of edges in bipartite graph (i.e. arcs)   */
      Gnum              bipedgenbr0;              /* Number of edges adjacent to part 0               */
      Gnum              bipedgenbr1;              /* Number of edges adjacent to part 1               */
      Gnum              bipvertnbr0;              /* Number of vertices in part 0                     */
      Gnum              bipvertnbr1;              /* Number of vertices in part 1                     */
      Gnum              bipvertnbrp;              /* Number of vertices in part put in first place    */
      Gnum              bippartval;               /* Part of bipartite graph to be put in first place */
      Gnum              compsizep;                /* Number of vertices to be removed from part p     */
      Gnum              compload01;               /* Load of vertices to be removed from both parts   */
      Gnum              comploadp;                /* Load of vertices to be removed from part p       */

      if ((actvnumtax = (Gnum *) memAlloc (actgrafdat.s.vertnbr * sizeof (Gnum))) == NULL) {
        errorPrint ("vgraphSeparateEs: out of memory (1)");
        return     (1);
      }
#ifdef SCOTCH_DEBUG_VGRAPH2
      memSet (actvnumtax, ~0, actgrafdat.s.vertnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_VGRAPH2 */
      actvnumtax -= actgrafdat.s.baseval;

      bipedgenbr  = 0;                            /* Initialize bipartite graph counts  */
      bipvertnbr0 =
      bipvertnbr1 = 0;
      for (actfronnum = 0; actfronnum < actgrafdat.fronnbr; actfronnum ++) { /* For all frontier vertices */
        Gnum              actvertnum;
        int               actpartval;
        Gnum              actedgenum;

        actvertnum = grafptr->frontab[actfronnum];
        actpartval = grafptr->parttax[actvertnum];

        if (actpartval == 0) {                    /* Count separator edges only for nodes of one side and multply by 2 */
          for (actedgenum = actgrafdat.s.verttax[actvertnum];
               actedgenum < actgrafdat.s.vendtax[actvertnum]; actedgenum ++)
            bipedgenbr += (actpartval ^ grafptr->parttax[actgrafdat.s.edgetax[actedgenum]]);
        }

        actvnumtax[actvertnum] = actpartval * (bipvertnbr1 - bipvertnbr0) + bipvertnbr0; /* Count and number separator vertices on each side */
        bipvertnbr0 += actpartval ^ 1;            /* Superscalar update */
        bipvertnbr1 += actpartval;
      }
      bipedgenbr *= 2;                            /* Count both sides of arcs */

      bipgrafdat.flagval = GRAPHFREEVERT | GRAPHVERTGROUP; /* Initialize bipartite graph structure */
      bipgrafdat.baseval = 0;                     /* Base bipartite graph from 0                   */
      bipgrafdat.vertnbr =
      bipgrafdat.vertnnd = bipvertnbr0 + bipvertnbr1;
      if (memAllocGroup ((void **) (void *)
                         &bipgrafdat.verttax, (size_t) ((bipgrafdat.vertnbr + 1) * sizeof (Gnum)),
                         &bipgrafdat.velotax, (size_t) ((actgrafdat.s.velotax != NULL) ? (bipgrafdat.vertnbr * sizeof (Gnum)) : 0),
                         &bipgrafdat.vnumtax, (size_t) (bipgrafdat.vertnbr * sizeof (Gnum)),
                         &bipgrafdat.edgetax, (size_t) (bipedgenbr * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("vgraphSeparateEs: out of memory (2)");
        memFree    (actvnumtax + actgrafdat.s.baseval);
        return     (1);
      }
      bipgrafdat.vendtax = bipgrafdat.verttax + 1;
      if (actgrafdat.s.velotax == NULL)
        bipgrafdat.velotax = NULL;
      bipgrafdat.vlbltax = NULL;
      bipgrafdat.edgenbr = bipedgenbr;
      bipgrafdat.edlotax = NULL;
      bipgrafdat.edlosum = bipedgenbr;
      bipgrafdat.degrmax = grafptr->s.degrmax;

      bippartval = (bipvertnbr0 <= bipvertnbr1) ? 0 : 1; /* Select smallest part to be placed first */
      if (bippartval == 0) {
        bipvertnbrp = bipvertnbr0;
        bipedgenbr0 = 0;
        bipedgenbr1 = bipedgenbr / 2;
      }
      else {
        bipvertnbrp = bipvertnbr1;
        bipedgenbr0 = bipedgenbr / 2;
        bipedgenbr1 = 0;
      }

      bipvelosum = 0;
      for (actfronnum = 0; actfronnum < actgrafdat.fronnbr; actfronnum ++) { /* For all frontier vertices */
        Gnum              actvertnum;
        int               actpartval;
        Gnum              bipvertnum;
        Gnum              actedgenum;

        actvertnum = grafptr->frontab[actfronnum];
        actpartval = grafptr->parttax[actvertnum];

        bipvertnum  = (actpartval ^ bippartval) * bipvertnbrp + actvnumtax[actvertnum];

        if (bipgrafdat.velotax != NULL) {
          Gnum              actveloval;

          actveloval  = actgrafdat.s.velotax[actvertnum];
          bipvelosum += actveloval;
          bipgrafdat.velotax[bipvertnum] = actveloval;
        }
        bipgrafdat.vnumtax[bipvertnum] = actvertnum;
        bipgrafdat.verttax[bipvertnum] = actpartval * (bipedgenbr1 - bipedgenbr0) + bipedgenbr0;

        for (actedgenum = actgrafdat.s.verttax[actvertnum]; /* Count separator edges */
             actedgenum < actgrafdat.s.vendtax[actvertnum]; actedgenum ++) {
          Gnum              actvertend;
          int               actpartend;

          actvertend = actgrafdat.s.edgetax[actedgenum];
          actpartend = grafptr->parttax[actvertend];
          if (actpartend != actpartval) {
            Gnum              bipedgenum;
            
#ifdef SCOTCH_DEBUG_VGRAPH2
            if (actvnumtax[actvertend] == ~0) {
              errorPrint ("vgraphSeparateEs: internal error (1)");
              graphExit  (&bipgrafdat);
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
            bipedgenum   = actpartval * (bipedgenbr1 - bipedgenbr0) + bipedgenbr0;
            bipedgenbr0 += actpartval ^ 1;        /* Superscalar update */
            bipedgenbr1 += actpartval;
            bipgrafdat.edgetax[bipedgenum] = actvnumtax[actvertend] + (actpartend ^ bippartval) * bipvertnbrp;
          }
        }
      }
      bipgrafdat.verttax[bipgrafdat.vertnbr] = bipgrafdat.edgenbr;
      bipgrafdat.velosum = (bipgrafdat.velotax != NULL) ? bipvelosum : bipgrafdat.vertnbr;

      memFree (actvnumtax + actgrafdat.s.baseval);

#ifdef SCOTCH_DEBUG_VGRAPH2
      if (((bipedgenbr0 - bipedgenbr1) * bippartval + bipedgenbr1) != bipgrafdat.edgenbr) {
        errorPrint ("vgraphSeparateEs: internal error (2)");
        graphExit  (&bipgrafdat);
        return     (1);
      }
      if (graphCheck (&bipgrafdat) != 0) {
        errorPrint ("vgraphSeparateEs: internal error (3)");
        graphExit  (&bipgrafdat);
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

      if (vgraphSeparateEsCover (&bipgrafdat, bipvertnbrp, grafptr->frontab, &grafptr->fronnbr) != 0) {
        errorPrint ("vgraphSeparateEs: cannot compute cover");
        graphExit  (&bipgrafdat);
        return     (1);
      }

      compsizep = 0;
      if (actgrafdat.s.velotax != NULL) {         /* If vertices are weighted */
        Gnum              fronnum;

        compload01 =
        comploadp  = 0;
        for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
          Gnum              bipvertnum;
          Gnum              actvertnum;
          Gnum              actveloval;

          bipvertnum = grafptr->frontab[fronnum];
          actvertnum = bipgrafdat.vnumtax[bipvertnum];
          actveloval = actgrafdat.s.velotax[actvertnum];
          grafptr->frontab[fronnum] = actvertnum; /* Express vertices with respect to original graph */

          grafptr->parttax[actvertnum] = 2;       /* Write separator part for global renumbering */
          compload01 += actveloval;
          if (bipvertnum < bipvertnbrp) {         /* Update separator vertices */
            compsizep ++;                         /* Superscalar update        */
            comploadp += actveloval;
          }
        }
      }
      else {                                      /* Vertices are not weighted */
        Gnum              fronnum;

        for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
          Gnum              bipvertnum;
          Gnum              actvertnum;

          bipvertnum = grafptr->frontab[fronnum];
          actvertnum = bipgrafdat.vnumtax[bipvertnum];
          grafptr->frontab[fronnum] = actvertnum; /* Express vertices with respect to original graph */

          grafptr->parttax[actvertnum] = 2;       /* Write separator part for global renumbering */
          if (bipvertnum < bipvertnbrp)           /* Update separator vertices                   */
            compsizep ++;                         /* Superscalar update                          */
        }
        compload01 = grafptr->fronnbr;            /* Loads are equivalent to sizes */
        comploadp  = compsizep;
      }
      grafptr->compsize[bippartval]     -= compsizep;
      grafptr->compsize[bippartval ^ 1] -= grafptr->fronnbr - compsizep;
      grafptr->compload[bippartval]     -= comploadp;
      grafptr->compload[bippartval ^ 1] -= compload01 - comploadp;

      graphExit (&bipgrafdat);
    }
    else {                                        /* Fat separator wanted                           */
      Gnum              compsize1;                /* Number of vertices to be removed from part 1   */
      Gnum              compload01;               /* Load of vertices to be removed from both parts */
      Gnum              compload1;                /* Load of vertices to be removed from part 1     */

      compsize1 = 0;
      grafptr->fronnbr = actgrafdat.fronnbr;      /* Keep separator as is */

      if (actgrafdat.s.velotax != NULL) {         /* If vertices are weighted */
        Gnum              fronnum;

        compload01 =
        compload1  = 0;
        for (fronnum = 0; fronnum < actgrafdat.fronnbr; fronnum ++) {
          Gnum              vertnum;
          Gnum              veloval;
          int               partval;

          vertnum = grafptr->frontab[fronnum];
          partval = grafptr->parttax[vertnum];
          veloval = grafptr->s.velotax[vertnum];

          compsize1  += partval;                  /* Superscalar update */
          compload01 += veloval;
          compload1  += partval * veloval;
          grafptr->parttax[vertnum] = 2;          /* Write separator part for global renumbering */
        }
      }
      else {                                      /* Vertices are not weighted */
        Gnum              fronnum;

        for (fronnum = 0; fronnum < actgrafdat.fronnbr; fronnum ++) {
          Gnum              vertnum;
          int               partval;

          vertnum = grafptr->frontab[fronnum];
          partval = grafptr->parttax[vertnum];

          compsize1 += partval;
          grafptr->parttax[vertnum] = 2;          /* Write separator part for global renumbering */
        }

        compload01 = actgrafdat.fronnbr;          /* Loads are equivalent to sizes */
        compload1  = compsize1;
      }

      grafptr->compsize[0] -= actgrafdat.fronnbr - compsize1; /* Update graph properties */
      grafptr->compsize[1] -= compsize1;
      grafptr->compload[0] -= compload01 - compload1;
      grafptr->compload[1] -= compload1;
    }
  }

  grafptr->comploaddlt = grafptr->compload[0] - grafptr->compload[1];
  grafptr->compload[2] = grafptr->s.velosum - grafptr->compload[0] - grafptr->compload[1];
  grafptr->fronnbr     = grafptr->s.vertnbr - grafptr->compsize[0] - grafptr->compsize[1];

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (grafptr) != 0) {
    errorPrint ("vgraphSeparateEs: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  return (0);
}
