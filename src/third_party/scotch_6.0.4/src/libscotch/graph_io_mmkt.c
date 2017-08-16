/* Copyright 2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_io_mmkt.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : This module contains the input/output   **/
/**                routines for handling the Matrix Market **/
/**                format.                                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 17 jan 2008     **/
/**                                 to   : 21 mar 2008     **/
/**                # Version 5.1  : from : 27 apr 2010     **/
/**                                 to   : 11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_IO_MMKT

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "graph_io_mmkt.h"

/********************************************************/
/*                                                      */
/* These routines handle source Matrix Market matrices. */
/*                                                      */
/********************************************************/

/* This routine loads a source graph from
** the given stream, corresponding to a MatrixMarket file.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphGeomLoadMmkt (
Graph * restrict const      grafptr,              /* Graph to load    */
Geom * restrict const       geomptr,              /* Geometry to load */
FILE * const                filesrcptr,           /* Topological data */
FILE * const                filegeoptr,           /* No use           */
const char * const          dataptr)              /* Fake base value  */
{
  Gnum                baseval;
  Gnum                mrownbr;
  Gnum                mcolnbr;
  Gnum                linenbr;
  Gnum                linenum;
  Gnum                vertnum;
  Gnum                verttmp;
  GraphGeomMmktEdge * sorttab;
  Gnum                sortnbr;
  Gnum                sortnum;
  Gnum *              edgetax;
  Gnum                edgenum;
  Gnum                edgetmp;
  Gnum                degrmax;
  char                linetab[1025];
  char *              lineptr;
  char                c;

  baseval = 1;                                    /* Regular MatrixMarket indices start from 1 */

  if ((dataptr != NULL)                        && /* If base value provided */
      (dataptr[0] != '\0')                     &&
      ((baseval = (Gnum) atol (dataptr)) == 0) && /* Get base value */
      (dataptr[0] != '0')) {
    errorPrint ("graphGeomLoadMmkt: invalid parameter");
    return     (1);
  }

  if (fgets (linetab, 1025, filesrcptr) == NULL) { /* Read header lines */
    errorPrint ("graphGeomLoadMmkt: bad input (1)");
    return     (1);
  }
  if (strncmp (linetab, "%%MatrixMarket", 14) != 0) {
    errorPrint ("graphGeomLoadMmkt: invalid header");
    return     (1);
  }

  for (lineptr = linetab + 14; *lineptr != '\0'; lineptr ++)
    *lineptr = tolower (*lineptr);

  if (strstr (linetab + 14, "matrix") == NULL) {
    errorPrint ("graphGeomLoadMmkt: only matrix types supported");
    return (1);
  }

  while ((c = fgetc (filesrcptr)) == '%') {       /* Skip additional comment lines */
    if (fgets (linetab, 1025, filesrcptr) == NULL) {
      errorPrint ("graphGeomLoadMmkt: bad input (2)");
      return (1);
    }
  }
  ungetc (c, filesrcptr);

  if ((intLoad (filesrcptr, &mrownbr) != 1) ||    /* Read number of rows    */
      (intLoad (filesrcptr, &mcolnbr) != 1) ||    /* Read number of columns */
      (intLoad (filesrcptr, &linenbr) != 1)) {    /* Read number of lines   */
    errorPrint ("graphGeomLoadMmkt: bad input (3)");
    return     (1);
  }
  if (mrownbr != mcolnbr) {                       /* If not a square matrix */
    errorPrint ("graphGeomLoadMmkt: not a square matrix");
    return     (1);
  }

  memSet (grafptr, 0, sizeof (Graph));

  grafptr->flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  grafptr->baseval = baseval;
  grafptr->vertnbr = mrownbr;
  grafptr->vertnnd = grafptr->vertnbr + baseval;

  if ((grafptr->verttax = memAlloc ((grafptr->vertnbr + 1) * sizeof (Gnum))) == NULL) {
    errorPrint ("graphGeomLoadMmkt: out of memory (1)");
    graphExit  (grafptr);
    return     (1);
  }
  grafptr->verttax -= baseval;
  grafptr->vendtax  = grafptr->verttax + 1;
  grafptr->velosum  = grafptr->vertnbr;

  if ((sorttab = (GraphGeomMmktEdge *) memAlloc (2 * linenbr * sizeof (GraphGeomMmktEdge))) == NULL) { /* Twice the space for symmetric edges */
    errorPrint ("graphGeomLoadMmkt: out of memory (2)");
    graphExit  (grafptr);
    return     (1);
  }
  grafptr->edgetax = ((Gnum *) sorttab) - baseval; /* TRICK: will be freed if graph is freed */

  for (linenum = sortnbr = 0; linenum < linenbr; linenum ++) {
    if ((intLoad (filesrcptr, &sorttab[sortnbr].vertnum[0]) != 1) || /* Read edge ends */
        (intLoad (filesrcptr, &sorttab[sortnbr].vertnum[1]) != 1) ||
        (fgets (linetab, 1025, filesrcptr) == NULL)) { /* Skip end of line */
      errorPrint ("graphGeomLoadMmkt: bad input (4)");
      graphExit  (grafptr);
      return     (1);
    }

    if ((sorttab[sortnbr].vertnum[0] < baseval) || (sorttab[sortnbr].vertnum[0] >= (mrownbr + baseval)) ||
        (sorttab[sortnbr].vertnum[1] < baseval) || (sorttab[sortnbr].vertnum[1] >= (mrownbr + baseval))) {
      errorPrint ("graphGeomLoadMmkt: bad input (5)");
      graphExit  (grafptr);
      return     (1);
    }

    if (sorttab[sortnbr].vertnum[0] != sorttab[sortnbr].vertnum[1]) { /* If not loop edge  */
      sorttab[sortnbr + 1].vertnum[0] = sorttab[sortnbr].vertnum[1]; /* Add symmetric edge */
      sorttab[sortnbr + 1].vertnum[1] = sorttab[sortnbr].vertnum[0];
      sortnbr += 2;
    }
  }

  intSort2asc2 (sorttab, sortnbr);                /* Sort edges by increasing indices */

  edgetax = grafptr->edgetax;                     /* TRICK: point to beginning of sorted edge array for re-use */
  for (sortnum = degrmax = 0, vertnum = baseval - 1, edgetmp = edgenum = baseval;
       sortnum < sortnbr; sortnum ++) {
    Gnum                vertend;

    if (vertnum < sorttab[sortnum].vertnum[0]) {  /* If change of vertex index, that is, first edge end */
      edgetmp = edgenum - edgetmp;                /* Compute degree and see if it is maximum degree     */
      if (edgetmp > degrmax)
        degrmax = edgetmp;
      edgetmp = edgenum;

      grafptr->verttax[++ vertnum] = edgenum;     /* Set beginning of new edge sub-array */

      while (vertnum < sorttab[sortnum].vertnum[0]) /* Fill gaps with isolated vertices */
        grafptr->verttax[++ vertnum] = edgenum;

      verttmp = baseval - 1;                      /* Make sure next edge will be considered as never seen before */
    }

    vertend = sorttab[sortnum].vertnum[1];        /* Get end of current edge                */
    if (vertend != verttmp)                       /* If edge differs from previous one      */
      edgetax[edgenum ++] = verttmp = vertend;    /* Add it to array and prevent duplicates */
  }
  edgetmp = edgenum - edgetmp;                    /* Compute degree and see if it is maximum degree */
  if (edgetmp > degrmax)
    degrmax = edgetmp;
  while (vertnum < mrownbr)                       /* Fill gaps with isolated vertices and mark beginning of new one */
    grafptr->verttax[++ vertnum] = edgenum;
  grafptr->verttax[++ vertnum] = edgenum;         /* Mark end of array */

#ifdef SCOTCH_DEBUG_GRAPH2
  if (vertnum != grafptr->vertnnd) {
    errorPrint ("graphGeomLoadMmkt: internal error (1)");
    graphExit  (grafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  grafptr->edgenbr = edgenum - baseval;
  grafptr->edgetax = ((Gnum *) memRealloc (edgetax + baseval, grafptr->edgenbr * sizeof (Gnum))) - baseval; /* TRICK: keep only useful space in re-used array */
  grafptr->edlotax = NULL;
  grafptr->edlosum = grafptr->edgenbr;
  grafptr->degrmax = degrmax;

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (grafptr) != 0) {                /* Check graph consistency */
    errorPrint ("graphGeomLoadMmkt: internal error (2)");
    graphExit  (grafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}

/* This routine saves the geometrical graph
** in the Matrix Market symmetric graph format.
** It returns:
** - 0   : on succes
** - !0  : on error.
*/

int
graphGeomSaveMmkt (
const Graph * restrict const  grafptr,            /* Graph to save    */
const Geom * restrict const   geomptr,            /* Geometry to save */
FILE * const                  filesrcptr,         /* Topological data */
FILE * const                  filegeoptr,         /* No use           */
const char * const            dataptr)            /* No use           */
{
  Gnum              baseadj;                      /* Base adjustment  */
  Gnum              vertnum;                      /* Current vertex   */
  Gnum              edgenum;                      /* Current edge     */
  int               o;

  baseadj = 1 - grafptr->baseval;                 /* Output base is always 1 */

  o = (fprintf (filesrcptr, "%%%%MatrixMarket matrix coordinate pattern symmetric\n%% Produced by Scotch graphGeomSaveMmkt\n" GNUMSTRING " " GNUMSTRING " " GNUMSTRING "\n", /* Write graph header */
                (Gnum) grafptr->vertnbr,
                (Gnum) grafptr->vertnbr,
                (Gnum) ((grafptr->edgenbr / 2) + grafptr->vertnbr)) == EOF);

  for (vertnum = grafptr->baseval; (o == 0) && (vertnum < grafptr->vertnnd); vertnum ++) {
    Gnum              vlblnum;                    /* Vertex label to output */

    vlblnum = ((grafptr->vlbltax != NULL) ? grafptr->vlbltax[vertnum] : vertnum) + baseadj;

    if (fprintf (filesrcptr, GNUMSTRING " " GNUMSTRING "\n", /* Write diagonal term */
                 (Gnum) vlblnum,
                 (Gnum) vlblnum) < 0) {
      o = 1;
      break;
    }

    for (edgenum = grafptr->verttax[vertnum]; edgenum < grafptr->vendtax[vertnum]; edgenum ++) {
      Gnum              vlblend;                  /* End vertex label to output */

      vlblend = grafptr->edgetax[edgenum];
      if (grafptr->vlbltax != NULL)
        vlblend = grafptr->vlbltax[vlblend];
      vlblend += baseadj;

      if (vlblend < vlblnum) {
        if (fprintf (filesrcptr, GNUMSTRING " " GNUMSTRING "\n",
                     (Gnum) vlblnum,
                     (Gnum) vlblend) < 0) {
          o = 1;
          break;
        }
      }
    }
  }
  if (o != 0)
    errorPrint ("graphGeomSaveMmkt: bad output");

  return (o);
}
