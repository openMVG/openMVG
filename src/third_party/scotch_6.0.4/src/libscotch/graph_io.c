/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_io.c                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                input/output functions.                 **/
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
/**                                 to     15 mar 1999     **/
/**                # Version 4.0  : from : 25 nov 2001     **/
/**                                 to     21 jan 2004     **/
/**                # Version 5.0  : from : 13 dec 2006     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_IO

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_io.h"

/*******************************************/
/*                                         */
/* These routines handle source graph I/O. */
/*                                         */
/*******************************************/

/* This routine loads a source graph from
** the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphLoad (
Graph * restrict const      grafptr,              /* Graph structure to fill              */
FILE * const                stream,               /* Stream from which to read graph data */
const Gnum                  baseval,              /* Base value (-1 means keep file base) */
const GraphFlag             flagval)              /* Graph loading flags                  */
{
  Gnum                edgenum;                    /* Number of edges really allocated */
  Gnum                edgennd;
  Gnum                vlblnbr;                    /* = vertnbr if vertex labels       */
  Gnum                vlblmax;                    /* Maximum vertex label number      */
  Gnum                velonbr;                    /* = vertnbr if vertex loads wanted */
  Gnum                velosum;                    /* Sum of vertex loads              */
  Gnum                edlonbr;                    /* = edgenbr if edge loads wanted   */
  Gnum                edlosum;                    /* Sum of edge loads                */
  Gnum                edgeval;                    /* Value where to read edge end     */
  Gnum                baseadj;
  Gnum                versval;
  Gnum                degrmax;
  Gnum                propval;
  char                proptab[4];
  Gnum                vertnum;

  memSet (grafptr, 0, sizeof (Graph));

  if (intLoad (stream, &versval) != 1) {          /* Read version number */
    errorPrint ("graphLoad: bad input (1)");
    return     (1);
  }
  if (versval != 0) {                             /* If version not zero */
    errorPrint ("graphLoad: old-style graph format no longer supported");
    return     (1);
  }

  if ((intLoad (stream, &grafptr->vertnbr) != 1) || /* Read rest of header */
      (intLoad (stream, &grafptr->edgenbr) != 1) ||
      (intLoad (stream, &baseadj)          != 1) ||
      (intLoad (stream, &propval)          != 1) ||
      (propval < 0)                              ||
      (propval > 111)) {
    errorPrint ("graphLoad: bad input (2)");
    return     (1);
  }
  sprintf (proptab, "%3.3d", (int) propval);      /* Compute file properties */
  proptab[0] -= '0';                              /* Vertex labels flag      */
  proptab[1] -= '0';                              /* Edge weights flag       */
  proptab[2] -= '0';                              /* Vertex loads flag       */

  grafptr->flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  if (baseval == -1) {                            /* If keep file graph base     */
    grafptr->baseval = baseadj;                   /* Set graph base as file base */
    baseadj          = 0;                         /* No base adjustment needed   */
  }
  else {                                          /* If set graph base     */
    grafptr->baseval = baseval;                   /* Set wanted graph base */
    baseadj          = baseval - baseadj;         /* Update base adjust    */
  }
  if (proptab[0] != 0)                            /* If vertex labels, no base adjust */
    baseadj = 0;

  velonbr = ((proptab[2] != 0) && ((flagval & GRAPHIONOLOADVERT) == 0)) ? grafptr->vertnbr : 0;
  vlblnbr = (proptab[0] != 0) ? grafptr->vertnbr : 0;
  edlonbr = ((proptab[1] != 0) && ((flagval & GRAPHIONOLOADEDGE) == 0)) ? grafptr->edgenbr : 0;

  if ((memAllocGroup ((void **) (void *)
                      &grafptr->verttax, (size_t) ((grafptr->vertnbr + 1) * sizeof (Gnum)),
                      &grafptr->velotax, (size_t) (velonbr                * sizeof (Gnum)),
                      &grafptr->vlbltax, (size_t) (vlblnbr                * sizeof (Gnum)), NULL) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &grafptr->edgetax, (size_t) (grafptr->edgenbr       * sizeof (Gnum)),
                      &grafptr->edlotax, (size_t) (edlonbr                * sizeof (Gnum)), NULL) == NULL)) {
    if (grafptr->verttax != NULL)
      memFree (grafptr->verttax);
    errorPrint ("graphLoad: out of memory");
    graphFree  (grafptr);
    return     (1);
  }
  grafptr->vertnnd  = grafptr->vertnbr + grafptr->baseval;
  grafptr->verttax -= grafptr->baseval;
  grafptr->vendtax  = grafptr->verttax + 1;       /* Use compact vertex array */
  grafptr->velotax  = (velonbr != 0) ? (grafptr->velotax - grafptr->baseval) : NULL;
  grafptr->vlbltax  = (vlblnbr != 0) ? (grafptr->vlbltax - grafptr->baseval) : NULL;
  grafptr->edgetax -= grafptr->baseval;
  grafptr->edlotax  = (edlonbr != 0) ? (grafptr->edlotax - grafptr->baseval) : NULL;

  vlblmax = grafptr->vertnnd - 1;                 /* No vertex labels known */
  velosum = (grafptr->velotax == NULL) ? grafptr->vertnbr : 0;
  edlosum = (grafptr->edlotax == NULL) ? grafptr->edgenbr : 0;
  edgennd = grafptr->edgenbr + grafptr->baseval;
  degrmax = 0;                                    /* No maximum degree yet */

  for (vertnum = edgenum = grafptr->baseval; vertnum < grafptr->vertnnd; vertnum ++) {
    Gnum                degrval;

    if (grafptr->vlbltax != NULL) {               /* If must read label               */
      Gnum                vlblval;                /* Value where to read vertex label */

      if (intLoad (stream, &vlblval) != 1) {      /* Read label data */
        errorPrint ("graphLoad: bad input (3)");
        graphFree  (grafptr);
        return     (1);
      }
      grafptr->vlbltax[vertnum] = vlblval;
      if (grafptr->vlbltax[vertnum] > vlblmax)    /* Get maximum vertex label */
        vlblmax = grafptr->vlbltax[vertnum];
    }
    if (proptab[2] != 0) {                        /* If must read vertex load        */
      Gnum                veloval;                /* Value where to read vertex load */

      if (intLoad (stream, &veloval) != 1) {      /* Read vertex load data    */
        errorPrint ("graphLoad: bad input (4)");
        graphFree  (grafptr);
        return     (1);
      }
      if (grafptr->velotax != NULL)
        velosum                  +=
        grafptr->velotax[vertnum] = veloval;
    }
    if (intLoad (stream, &degrval) != 1) {        /* Read vertex degree */
      errorPrint ("graphLoad: bad input (5)");
      graphFree  (grafptr);
      return     (1);
    }
    if (degrmax < degrval)                        /* Set maximum degree */
      degrmax = degrval;

    grafptr->verttax[vertnum] = edgenum;          /* Set index in edge array */
    degrval += edgenum;
    if (degrval > edgennd) {                      /* Check if edge array overflows */
      errorPrint ("graphLoad: invalid arc count (1)");
      graphFree  (grafptr);
      return     (1);
    }

    for ( ; edgenum < degrval; edgenum ++) {
      if (proptab[1] != 0) {                      /* If must read edge load        */
        Gnum                edloval;              /* Value where to read edge load */

        if (intLoad (stream, &edloval) != 1) {    /* Read edge load data    */
          errorPrint ("graphLoad: bad input (6)");
          graphFree  (grafptr);
          return     (1);
        }
        if (grafptr->edlotax != NULL)
          edlosum                  +=
          grafptr->edlotax[edgenum] = (Gnum) edloval;
      }
      if (intLoad (stream, &edgeval) != 1) {      /* Read edge data */
        errorPrint ("graphLoad: bad input (7)");
        graphFree  (grafptr);
        return     (1);
      }
      grafptr->edgetax[edgenum] = edgeval + baseadj;
    }
  }
  grafptr->verttax[vertnum] = edgenum;            /* Set end of edge array */
  if (edgenum != edgennd) {                       /* Check if number of edges is valid */
    errorPrint ("graphLoad: invalid arc count (2)");
    graphFree  (grafptr);
    return     (1);
  }
  grafptr->velosum = velosum;
  grafptr->edlosum = edlosum;
  grafptr->degrmax = degrmax;

  if (grafptr->vlbltax != NULL) {                 /* If vertex label renaming necessary       */
    if (graphLoad2 (grafptr->baseval, grafptr->vertnnd, grafptr->verttax, /* Rename edge ends */
                    grafptr->vendtax, grafptr->edgetax, vlblmax, grafptr->vlbltax) != 0) {
      errorPrint ("graphLoad: cannot relabel vertices");
      graphFree  (grafptr);
      return     (1);
    }
  }

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (grafptr) != 0) {                /* Check graph consistency */
    errorPrint ("graphLoad: inconsistent graph data");
    graphFree  (grafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}

int
graphLoad2 (
const Gnum                  baseval,
const Gnum                  vertnnd,
const Gnum * const          verttax,
const Gnum * const          vendtax,
Gnum * restrict const       edgetax,
const Gnum                  vlblmax,
const Gnum * const          vlbltax)
{
  Gnum                vertnum;                    /* Number of current vertex        */
  Gnum * restrict     indxtab;                    /* Vertex label/number index table */

  if ((indxtab = (Gnum *) memAlloc ((vlblmax + 1) * sizeof (Gnum))) == NULL) {
    errorPrint  ("graphLoad2: out of memory");
    return      (1);
  }

  memSet (indxtab, ~0, (vlblmax + 1) * sizeof (Gnum)); /* Assume labels not used */
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    if (indxtab[vlbltax[vertnum]] != ~0) {        /* If vertex label already used */
      errorPrint  ("graphLoad2: duplicate vertex label");
      memFree     (indxtab);
      return      (1);
    }
    indxtab[vlbltax[vertnum]] = vertnum;          /* Set vertex number index */
  }
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    Gnum                edgenum;                  /* Number of current edge */

    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
      if (edgetax[edgenum] > vlblmax) {           /* If invalid edge end number */
        errorPrint ("graphLoad2: invalid arc end number (1)");
        memFree    (indxtab);
        return     (1);
      }
      if (indxtab[edgetax[edgenum]] == ~0) {      /* If unused edge end number */
        errorPrint ("graphLoad2: invalid arc end number (2)");
        memFree    (indxtab);
        return     (1);
      }
      edgetax[edgenum] = indxtab[edgetax[edgenum]]; /* Replace label by number */
    }
  }

  memFree (indxtab);                              /* Free index array */

  return (0);
}

/* This routine saves a source graph to
** the given stream, in the new-style
** graph format.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphSave (
const Graph * const         grafptr,
FILE * const                stream)
{
  Gnum                vertnum;
  char                propstr[4];                 /* Property string */
  int                 o;

  propstr[0] = (grafptr->vlbltax != NULL) ? '1' : '0'; /* Set property string */
  propstr[1] = (grafptr->edlotax != NULL) ? '1' : '0';
  propstr[2] = (grafptr->velotax != NULL) ? '1' : '0';
  propstr[3] = '\0';

  if (fprintf (stream, "0\n" GNUMSTRING "\t" GNUMSTRING "\n" GNUMSTRING "\t%3s\n", /* Write file header */
               (Gnum) grafptr->vertnbr,
               (Gnum) grafptr->edgenbr,
               (Gnum) grafptr->baseval,
               propstr) == EOF) {
    errorPrint ("graphSave: bad output (1)");
    return     (1);
  }

  for (vertnum = grafptr->baseval, o = 0;
       (vertnum < grafptr->vertnnd) && (o == 0); vertnum ++) {
    Gnum                edgenum;

    if (grafptr->vlbltax != NULL)                 /* Write vertex label if necessary */
      o  = (fprintf (stream, GNUMSTRING "\t", (Gnum) grafptr->vlbltax[vertnum]) == EOF);
    if (grafptr->velotax != NULL)                 /* Write vertex load if necessary */
      o |= (fprintf (stream, GNUMSTRING "\t", (Gnum) grafptr->velotax[vertnum]) == EOF);

    o |= (fprintf (stream, GNUMSTRING, (Gnum) (grafptr->vendtax[vertnum] - grafptr->verttax[vertnum])) == EOF); /* Write vertex degree */

    for (edgenum = grafptr->verttax[vertnum];
         (edgenum < grafptr->vendtax[vertnum]) && (o == 0); edgenum ++) {
      Gnum                vertend;

      o |= (putc ('\t', stream) == EOF);
      if (grafptr->edlotax != NULL)               /* Write edge load if necessary */
        o |= (fprintf (stream, GNUMSTRING "\t", (Gnum) grafptr->edlotax[edgenum]) == EOF);
      vertend = grafptr->edgetax[edgenum];
      o |= (fprintf (stream, GNUMSTRING, (Gnum) ((grafptr->vlbltax != NULL) ? grafptr->vlbltax[vertend] : vertend)) == EOF); /* Write edge end */
    }
    o |= (putc ('\n', stream) == EOF);
  }

  if (o != 0)
    errorPrint ("graphSave: bad output (2)");

  return (o);
}
