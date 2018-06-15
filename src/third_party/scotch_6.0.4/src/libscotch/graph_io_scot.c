/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_io_scot.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the I/O routines   **/
/**                for handling the Scotch graph format.   **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 06 nov 1997     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 13 dec 1998     **/
/**                                 to     21 dec 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     22 dec 2005     **/
/**                # Version 5.0  : from : 13 sep 2006     **/
/**                                 to     27 feb 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_IO_SCOT

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "graph_io_scot.h"

/* This routine loads the geometrical graph
** in the Scotch graph format, and allocates
** the proper structures.
** - 0   : on success.
** - !0  : on error.
*/

int
graphGeomLoadScot (
Graph * restrict const      grafptr,              /* Graph to load    */
Geom * restrict const       geomptr,              /* Geometry to load */
FILE * const                filesrcptr,           /* Topological data */
FILE * const                filegeoptr,           /* No use           */
const char * const          dataptr)              /* No use           */
{
  void *                        coorfileptr;      /* Temporary pointer to comply with C99 rules */
  double * restrict             coorfiletab;      /* Pointer to geometric data read from file   */
  GraphGeomScotSort * restrict  coorsorttab;      /* Pointer to geometric data sorting array    */
  int                           coorsortflag;     /* Flag set if geometric data sorted by label */
  Gnum                          coornbr;          /* Number of geometric coordinates in file    */
  Gnum                          coornum;          /* Number of current coordinate               */
  GraphGeomScotSort * restrict  vertsorttab;      /* Pointer to graph sorting array             */
  int                           vertsortflag;     /* Flag set if graph data sorted by label     */
  Gnum                          vertnum;          /* Current graph vertex                       */
  Gnum                          dimnnbr;          /* Dimension of geometry file                 */
  int                           o;

  if (filesrcptr != NULL) {
    if (graphLoad (grafptr, filesrcptr, -1, 0) != 0)
      return (1);
  }

  if (filegeoptr == NULL)
    return (0);

  if ((intLoad (filegeoptr, &dimnnbr) != 1) ||    /* Read type and number of geometry items */
      (intLoad (filegeoptr, &coornbr) != 1) ||
      (dimnnbr < 1)                         ||
      (dimnnbr > 3)) {
    errorPrint ("graphGeomLoadScot: bad input (1)");
    return     (1);
  }
  if ((filesrcptr != NULL) && (grafptr->vertnbr != coornbr)) {
    errorPrint ("graphGeomLoadScot: inconsistent number of vertices");
    return     (1);
  }

  if (grafptr->vertnbr == 0)
    return (0);

  if ((geomptr->geomtab == NULL) &&               /* Allocate geometry if necessary */
      ((geomptr->geomtab = (double *) memAlloc (grafptr->vertnbr * dimnnbr * sizeof (double))) == NULL)) {
    errorPrint ("graphGeomLoadScot: out of memory (1)");
    return     (1);
  }

  if (memAllocGroup ((void **)
                     &coorfileptr, (size_t) (coornbr * dimnnbr * sizeof (double)),
                     &coorsorttab, (size_t) (coornbr           * sizeof (GraphGeomScotSort)),
                     &vertsorttab, (size_t) (grafptr->vertnbr  * sizeof (GraphGeomScotSort)), NULL) == NULL) {
    errorPrint ("graphGeomLoadScot: out of memory (2)");
    return     (1);
  }
  coorfiletab = coorfileptr;

  o = 0;
  coorsortflag = 1;                               /* Assume geometry data sorted */
  for (coornum = 0; (o == 0) && (coornum < coornbr); coornum ++) {
    Gnum                vlblnum;

    o = 1 - intLoad (filegeoptr, &vlblnum);
    coorsorttab[coornum].labl = vlblnum;
    coorsorttab[coornum].num  = coornum;
    if ((coornum > 0) &&                          /* Check if geometry data sorted */
        (coorsorttab[coornum].labl < coorsorttab[coornum - 1].labl))
      coorsortflag = 0;                           /* Geometry data not sorted */

    o |= 1 - fscanf (filegeoptr, "%lf",           /* Read X coordinate */
                     &coorfiletab[coornum * dimnnbr]);
    if (dimnnbr > 1) {
      o |= 1 - fscanf (filegeoptr, "%lf",         /* Read Y coordinate */
                       &coorfiletab[(coornum * dimnnbr) + 1]);
      if (dimnnbr > 2)
        o |= 1 - fscanf (filegeoptr, "%lf",       /* Read Z coordinate */
                         &coorfiletab[(coornum * dimnnbr) + 2]);
    }
  }
  if (o != 0) {
    errorPrint ("graphGeomLoadScot: bad input (2)");
    memFree    (coorfiletab);                     /* Free group leader */
    return     (1);
  }

  if (coorsortflag != 1)                          /* If geometry data not sorted        */
    intSort2asc1 (coorsorttab, coornbr);          /* Sort sort area by ascending labels */

  for (coornum = 1; coornum < coornbr; coornum ++) { /* Check geometric data integrity */
    if (coorsorttab[coornum].labl == coorsorttab[coornum - 1].labl) {
      errorPrint ("graphGeomLoadScot: duplicate vertex label");
      memFree    (coorfiletab);                   /* Free group leader */
      return     (1);
    }
  }

  if (grafptr->vlbltax != NULL) {                 /* If graph has vertex labels */
    vertsortflag = 1;                             /* Assume graph data sorted   */
    for (vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++) {
      vertsorttab[vertnum].labl = grafptr->vlbltax[vertnum + grafptr->baseval];
      vertsorttab[vertnum].num  = vertnum;
      if ((vertnum > 0) &&                        /* Check if graph data sorted */
          (vertsorttab[vertnum].labl < vertsorttab[vertnum - 1].labl))
        vertsortflag = 0;                         /* Graph data not sorted */
    }
    if (vertsortflag != 1)                        /* If graph data not sorted             */
      intSort2asc1 (vertsorttab, grafptr->vertnbr); /* Sort sort area by ascending labels */
  }
  else {                                          /* Graph does not have vertex labels */
    for (vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++)
      vertsorttab[vertnum].labl =
      vertsorttab[vertnum].num  = vertnum;
  }

  for (coornum = vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++) { /* For all vertices in graph */
    while ((coornum < coornbr) && (coorsorttab[coornum].labl < vertsorttab[vertnum].labl))
      coornum ++;                                 /* Search geometry vertex with same label                           */
    if ((coornum >= coornbr) || (coorsorttab[coornum].labl > vertsorttab[vertnum].labl)) { /* If label does not exist */
      errorPrint ("graphGeomLoadScot: vertex geometry data not found (%d)",
                  vertsorttab[vertnum].labl);
      memFree    (coorfiletab);                   /* Free group leader */
      return     (1);
    }
    memCpy (&geomptr->geomtab[vertsorttab[vertnum].num * dimnnbr], &coorfiletab[coorsorttab[coornum ++].num * dimnnbr], dimnnbr * sizeof (double));
  }

  memFree (coorfiletab);                          /* Free group leader */

  return (0);
}

/* This routine saves the source process graph
** in the Scotch source and geometry formats.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphGeomSaveScot (
const Graph * restrict const  grafptr,            /* Graph to save    */
const Geom * restrict const   geomptr,            /* Geometry to save */
FILE * const                  filesrcptr,         /* Topological data */
FILE * const                  filegeoptr,         /* No use           */
const char * const            dataptr)            /* No use           */
{
  Gnum              vertnum;
  int               dimnnbr;
  int               o;

  if (filesrcptr != NULL) {
    if (graphSave (grafptr, filesrcptr) != 0)     /* Save graph structural data */
      return (1);
  }

  dimnnbr = geomptr->dimnnbr;

  o = 0;
  if (geomptr->geomtab != NULL) {                 /* If geometrical data present     */
    o = (fprintf (filegeoptr, GNUMSTRING "\n" GNUMSTRING "\n", /* Output file header */
                  (Gnum) geomptr->dimnnbr,
                  (Gnum) grafptr->vertnbr) == EOF);

    switch (dimnnbr) {                            /* Output geometry data */
      case 1 :
        for (vertnum = grafptr->baseval; (o == 0) && (vertnum < grafptr->vertnnd); vertnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\n",
                         (Gnum) ((grafptr->vlbltax != NULL) ? grafptr->vlbltax[vertnum] : vertnum),
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr]) == EOF);
        break;
      case 2 :
        for (vertnum = grafptr->baseval; (o == 0) && (vertnum < grafptr->vertnnd); vertnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\t%lf\n",
                         (Gnum) ((grafptr->vlbltax != NULL) ? grafptr->vlbltax[vertnum] : vertnum),
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr],
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr + 1]) == EOF);
        break;
      case 3 :
        for (vertnum = grafptr->baseval; (o == 0) && (vertnum < grafptr->vertnnd); vertnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\t%lf\t%lf\n",
                         (Gnum) ((grafptr->vlbltax != NULL) ? grafptr->vlbltax[vertnum] : vertnum),
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr],
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr + 1],
                         (double) geomptr->geomtab[(vertnum - grafptr->baseval) * dimnnbr + 2]) == EOF);
        break;
#ifdef SCOTCH_DEBUG_GRAPH2
      default :
        errorPrint ("graphGeomSaveScot: invalid geometry type");
        return     (1);
#endif /* SCOTCH_DEBUG_GRAPH2 */
    }

    if (o != 0) {
      errorPrint ("graphGeomSaveScot: bad output");
    }
  }

  return (o);
}
