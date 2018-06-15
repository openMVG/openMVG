/* Copyright 2004,2007-2009,2011,2013-2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_induce.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                subgraph-making functions.              **/
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
/**                # Version 4.0  : from : 28 nov 2001     **/
/**                                 to     17 apr 2006     **/
/**                # Version 5.0  : from : 14 dec 2006     **/
/**                                 to     11 jun 2008     **/
/**                # Version 5.1  : from : 01 jan 2009     **/
/**                                 to     01 jan 2009     **/
/**                # Version 6.0  : from : 29 mar 2011     **/
/**                                 to     28 feb 2015     **/
/**                                                        **/
/**   NOTES      : # Several algorithms, such as the       **/
/**                  active graph building routine of      **/
/**                  bgraphInit2, assume that, for every   **/
/**                  vertex, remaining edges will be kept  **/
/**                  in the same order as in the original  **/
/**                  graph. This must be enforced.         **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH
#define GRAPH_INDUCE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_induce.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

/* This routine builds the graph induced
** by the original graph and the list of
** selected vertices.
** The induced vnumtab array is the list
** array if the original graph does not have
** a vnumtab, or the proper subset of the
** original vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphInduceList (
const Graph * restrict const    orggrafptr,
const VertList * restrict const indlistptr,
Graph * restrict const          indgrafptr)
{
  Gnum * restrict       orgindxtax;               /* Based access to vertex translation array       */
  Gnum                  indvertnbr;               /* Number of vertices in induced graph            */
  Gnum                  indvertnnd;
  Gnum                  indvertnum;               /* Number of current vertex in induced graph      */
  const Gnum * restrict indvnumtax;
  Gnum * restrict       indedgetab;               /* Pointer to pre-allocated edge array            */
  Gnum                  indedgenbr;               /* (Approximate) number of edges in induced graph */

  const Gnum * restrict const orgverttax = orggrafptr->verttax;
  const Gnum * restrict const orgvendtax = orggrafptr->vendtax;

  indvertnbr = indlistptr->vnumnbr;

  memSet (indgrafptr, 0, sizeof (Graph));         /* Initialize graph fields */
  indgrafptr->flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  indgrafptr->baseval = orggrafptr->baseval;

  if (orggrafptr->velotax != NULL) {
    if (memAllocGroup ((void **) (void *)
                       &indgrafptr->verttax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)),
                       &indgrafptr->vnumtax, (size_t) ( indvertnbr      * sizeof (Gnum)),
                       &indgrafptr->velotax, (size_t) ( indvertnbr      * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("graphInduceList: out of memory (1)");
      return     (1);                             /* Nothing to free because group allocation failed */
    }
    indgrafptr->velotax -= indgrafptr->baseval;
  }
  else {
    if (memAllocGroup ((void **) (void *)
                       &indgrafptr->verttax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)),
                       &indgrafptr->vnumtax, (size_t) ( indvertnbr      * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("graphInduceList: out of memory (2)");
      return     (1);
    }
  }
  indgrafptr->verttax -= indgrafptr->baseval;     /* Adjust base of arrays */
  indgrafptr->vnumtax -= indgrafptr->baseval;
  indgrafptr->vertnbr  = indvertnbr;
  indgrafptr->vertnnd  = indvertnbr + indgrafptr->baseval;

  indedgenbr = orggrafptr->edgenbr;               /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->degrmax > 0) && (indvertnbr < (indedgenbr / orggrafptr->degrmax)))
    indedgenbr = indvertnbr * orggrafptr->degrmax;
  if (orggrafptr->edlotax != NULL)                /* If graph has edge weights */
    indedgenbr *= 2;                              /* Account for edge weights  */

  if (memAllocGroup ((void *)
                     &indedgetab, (size_t) (indedgenbr          * sizeof (Gnum)), /* Pre-allocate space for edgetab (and edlotab)          */
                     &orgindxtax, (size_t) (orggrafptr->vertnbr * sizeof (Gnum)), NULL) == NULL) { /* orgindxtab is at the end of the heap */
    errorPrint ("graphInduceList: out of memory (3)");
    graphExit  (indgrafptr);
    return     (1);
  }
  orgindxtax -= orggrafptr->baseval;

  memCpy (indgrafptr->vnumtax + indgrafptr->baseval, /* Copy vertex number array from list */
          indlistptr->vnumtab, indvertnbr * sizeof (Gnum));

  memSet (orgindxtax + orggrafptr->baseval, ~0, orggrafptr->vertnbr * sizeof (Gnum)); /* Preset index array */

  indvnumtax = indgrafptr->vnumtax;
  for (indvertnum = indgrafptr->baseval, indvertnnd = indvertnum + indvertnbr, indedgenbr = 0; /* Fill index array */
       indvertnum < indvertnnd; indvertnum ++) {
    Gnum                orgvertnum;

    orgvertnum = indvnumtax[indvertnum];

    orgindxtax[orgvertnum] = indvertnum;          /* Mark selected vertices */
    indedgenbr += orgvendtax[orgvertnum] - orgverttax[orgvertnum];
  }

  return (graphInduce2 (orggrafptr, indgrafptr, indvertnbr, indedgenbr, indedgetab, orgindxtax));
}

/* This routine builds the graph induced
** by the original graph and the vector
** of selected vertices.
** The induced vnumtab array is the list of
** selected vertices if the original graph
** does not have a vnumtab, or the proper
** subset of the original vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphInducePart (
const Graph * restrict const  orggrafptr,         /* Pointer to original graph             */
const GraphPart * const       orgparttax,         /* Based array of vertex partition flags */
const Gnum                    indvertnbr,         /* Number of vertices in selected part   */
const GraphPart               indpartval,         /* Partition value of vertices to keep   */
Graph * restrict const        indgrafptr)         /* Pointer to induced subgraph           */
{
  Gnum * restrict     orgindxtax;                 /* Based access to vertex translation array       */
  Gnum                indvertnum;                 /* Number of current vertex in induced graph      */
  Gnum * restrict     indvnumtax;
  Gnum * restrict     indedgetab;                 /* Pointer to pre-allocated edge array            */
  Gnum                indedgenbr;                 /* (Approximate) number of edges in induced graph */
  Gnum                orgvertnum;

  const Gnum * restrict const orgverttax = orggrafptr->verttax;
  const Gnum * restrict const orgvendtax = orggrafptr->vendtax;

  memSet (indgrafptr, 0, sizeof (Graph));         /* Initialize graph fields */
  indgrafptr->flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  indgrafptr->baseval = orggrafptr->baseval;

  if (orggrafptr->velotax != NULL) {
    if (memAllocGroup ((void **) (void *)
                       &indgrafptr->verttax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)),
                       &indgrafptr->vnumtax, (size_t) ( indvertnbr      * sizeof (Gnum)),
                       &indgrafptr->velotax, (size_t) ( indvertnbr      * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("graphInducePart: out of memory (1)");
      return     (1);                             /* Nothing to free because group allocation failed */
    }
    indgrafptr->velotax -= indgrafptr->baseval;
  }
  else {
    if (memAllocGroup ((void **) (void *)
                       &indgrafptr->verttax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)),
                       &indgrafptr->vnumtax, (size_t) ( indvertnbr      * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("graphInducePart: out of memory (2)");
      return     (1);
    }
  }
  indgrafptr->verttax -= indgrafptr->baseval;     /* Adjust base of arrays */
  indgrafptr->vnumtax -= indgrafptr->baseval;
  indgrafptr->vertnbr  = indvertnbr;
  indgrafptr->vertnnd  = indvertnbr + indgrafptr->baseval;

  indedgenbr = orggrafptr->edgenbr;               /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->degrmax > 0) && (indvertnbr < (indedgenbr / orggrafptr->degrmax)))
    indedgenbr = indvertnbr * orggrafptr->degrmax;
  if (orggrafptr->edlotax != NULL)                /* If graph has edge weights */
    indedgenbr *= 2;                              /* Account for edge weights  */

  if (memAllocGroup ((void *)
                     &indedgetab, (size_t) (indedgenbr          * sizeof (Gnum)), /* Pre-allocate space for edgetab (and edlotab)          */
                     &orgindxtax, (size_t) (orggrafptr->vertnbr * sizeof (Gnum)), NULL) == NULL) { /* orgindxtab is at the end of the heap */
    errorPrint ("graphInducePart: out of memory (3)");
    graphExit  (indgrafptr);
    return     (1);
  }
  orgindxtax -= orggrafptr->baseval;

  indvnumtax = indgrafptr->vnumtax;
  for (orgvertnum = indvertnum = orggrafptr->baseval, indedgenbr = 0; /* Fill index array */
       orgvertnum < orggrafptr->vertnnd; orgvertnum ++) {
    if (orgparttax[orgvertnum] == indpartval) {   /* If vertex should be kept */
      orgindxtax[orgvertnum] = indvertnum;        /* Mark selected vertex     */
      indvnumtax[indvertnum] = orgvertnum;
      indedgenbr += orgvendtax[orgvertnum] - orgverttax[orgvertnum];
      indvertnum ++;                              /* One more induced vertex created */
    }
    else
      orgindxtax[orgvertnum] = ~0;
  }
#ifdef SCOTCH_DEBUG_GRAPH2
  if ((indvertnum - indgrafptr->baseval) != indvertnbr) {
    errorPrint ("graphInducePart: inconsistent data");
    memFree    (indedgetab);
    graphExit  (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (graphInduce2 (orggrafptr, indgrafptr, indvertnbr, indedgenbr, indedgetab, orgindxtax));
}

/* This routine finalizes the building
** of the induced subgraph.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
graphInduce2 (
const Graph * restrict const  orggrafptr,         /* Pointer to original graph                          */
Graph * restrict const        indgrafptr,         /* Pointer to induced graph                           */
const Gnum                    indvertnbr,         /* Number of vertices in induced graph                */
const Gnum                    indedgenbr,         /* (Upper bound of) number of edges in induced graph  */
Gnum * const                  indedgetab,         /* Pointer to pre-allocated edge and edge load arrays */
const Gnum * restrict const   orgindxtax)         /* Array of numbers of selected vertices              */
{
  Gnum                indvertnum;                 /* Current induced vertex number              */
  Gnum                indvelosum;                 /* Overall induced vertex load                */
  Gnum                indedlosum;                 /* Overall induced edge load                  */
  Gnum                indedgenum;                 /* Number of current induced edge             */
  Gnum                orgvertnum;                 /* Number of current vertex in original graph */
  Gnum                orgedgenum;                 /* Number of current edge in original graph   */

  const Gnum * restrict const orgverttax = orggrafptr->verttax;
  const Gnum * restrict const orgvendtax = orggrafptr->vendtax;
  const Gnum * restrict const orgvelotax = orggrafptr->velotax;
  const Gnum * restrict const orgvnumtax = orggrafptr->vnumtax;
  const Gnum * restrict const orgedgetax = orggrafptr->edgetax;
  const Gnum * restrict const orgedlotax = orggrafptr->edlotax;
  Gnum * restrict const       indverttax = indgrafptr->verttax;
  Gnum * restrict const       indvelotax = indgrafptr->velotax;
  Gnum * restrict const       indvnumtax = indgrafptr->vnumtax;
  Gnum * restrict             indedgetax;
  Gnum * restrict             indedlotax;

  if (orgedlotax != NULL) {
    memOffset ((void *) indedgetab,
               &indedgetax, (size_t) (indedgenbr * sizeof (Gnum)),
               &indedlotax, (size_t) (indedgenbr * sizeof (Gnum)), NULL);
    indedgetax -= indgrafptr->baseval;
    indedlotax -= indgrafptr->baseval;
  }
  else {
    indedgetax = indedgetab - indgrafptr->baseval;
    indedlotax = NULL;
  }

  indvelosum = (indvelotax == NULL) ? indgrafptr->vertnbr : 0;
  indedlosum = 0;
  for (indvertnum = indedgenum = indgrafptr->baseval;
       indvertnum < indgrafptr->vertnnd; indvertnum ++) {
    orgvertnum = indvnumtax[indvertnum];
    indverttax[indvertnum] = indedgenum;
    if (indvelotax != NULL) {                     /* If graph has vertex weights */
      indvelosum +=                               /* Accumulate vertex loads     */
      indvelotax[indvertnum] = orgvelotax[orgvertnum];
    }

    if (indedlotax != NULL) {                     /* If graph has edge weights */
      for (orgedgenum = orgverttax[orgvertnum];
           orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
        if (orgindxtax[orgedgetax[orgedgenum]] != ~0) { /* If edge should be kept */
          indedlosum                     +=
          indedlotax[indedgenum] = orgedlotax[orgedgenum];
          indedgetax[indedgenum] = orgindxtax[orgedgetax[orgedgenum]];
          indedgenum ++;
        }
      }
    }
    else {
      for (orgedgenum = orgverttax[orgvertnum];
           orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
        if (orgindxtax[orgedgetax[orgedgenum]] != ~0) { /* If edge should be kept */
          indedgetax[indedgenum] = orgindxtax[orgedgetax[orgedgenum]];
          indedgenum ++;
        }
      }
    }
  }
  indverttax[indvertnum] = indedgenum;            /* Mark end of edge array */

  indgrafptr->vendtax = indgrafptr->verttax + 1;  /* Use compact representation of vertex arrays */
  indgrafptr->vertnbr = indvertnum - indgrafptr->baseval;
  indgrafptr->vertnnd = indvertnum;
  indgrafptr->velosum = indvelosum;
  indgrafptr->edgenbr = indedgenum - indgrafptr->baseval; /* Set actual number of edges */
  indgrafptr->edlosum = (indedlotax != NULL) ? indedlosum : indgrafptr->edgenbr;
  indgrafptr->degrmax = orggrafptr->degrmax;      /* Induced maximum degree is likely to be the one of the original graph */

  if (orggrafptr->vnumtax != NULL) {              /* Adjust vnumtax */
    for (indvertnum = indgrafptr->baseval; indvertnum < indgrafptr->vertnnd; indvertnum ++)
      indvnumtax[indvertnum] = orgvnumtax[indvnumtax[indvertnum]];
  }

  if (indedlotax != NULL) {                       /* Re-allocate arrays and delete orgindxtab             */
    size_t              indedlooftval;            /* Offset of edge load array with respect to edge array */

    indedlooftval = indedlotax - indedgetax;
    indgrafptr->edgetax = (Gnum *) memRealloc (indedgetab, (indedlooftval + indgrafptr->edgenbr) * sizeof (Gnum)) - indgrafptr->baseval;
    indgrafptr->edlotax = indgrafptr->edgetax + indedlooftval; /* Use old index into old array as new index */
  }
  else
    indgrafptr->edgetax = (Gnum *) memRealloc (indedgetab, indgrafptr->edgenbr * sizeof (Gnum)) - indgrafptr->baseval;

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (indgrafptr) != 0) {             /* Check graph consistency */
    errorPrint ("graphInduce2: inconsistent graph data");
    graphExit  (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}
