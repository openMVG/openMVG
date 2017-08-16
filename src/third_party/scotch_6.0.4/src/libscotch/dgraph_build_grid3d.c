/* Copyright 2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_build_grid3d.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (5.0)                  **/
/**                                                        **/
/**   FUNCTION   : These lines are the distributed source  **/
/**                graph building routines for 3D grid     **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 21 jul 2005     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 05 jun 2010     **/
/**                                 to   : 06 jun 2010     **/
/**                                                        **/
/************************************************************/

#define DGRAPH_BUILD_GRID3D

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_build_grid3d.h"

/**************************************/
/*                                    */
/* Vertex neighbor handling routines. */
/*                                    */
/**************************************/

/*
**
*/

static
Gnum
dgraphBuildGrid3Dvertex26M (
const DgraphBuildGrid3DData * restrict const  dataptr,
const Gnum                                    vertglbnum,
Gnum                                          edgelocnum,
const Gnum                                    posxval,
const Gnum                                    posyval,
const Gnum                                    poszval)
{
  Gnum                ngbxmin;
  Gnum                ngbxmax;
  Gnum                ngbxval;
  Gnum                ngbymin;
  Gnum                ngbymax;
  Gnum                ngbyval;
  Gnum                ngbzmin;
  Gnum                ngbzmax;
  Gnum                ngbzval;

  ngbxmin = (posxval > 0) ? -1 : 0;
  ngbymin = (posyval > 0) ? -1 : 0;
  ngbzmin = (poszval > 0) ? -1 : 0;
  ngbxmax = (posxval < (dataptr->dimxval - 1)) ? 1 : 0;
  ngbymax = (posyval < (dataptr->dimyval - 1)) ? 1 : 0;
  ngbzmax = (poszval < (dataptr->dimzval - 1)) ? 1 : 0;

  for (ngbzval = ngbzmin; ngbzval <= ngbzmax; ngbzval ++) {
    for (ngbyval = ngbymin; ngbyval <= ngbymax; ngbyval ++) {
      for (ngbxval = ngbxmin; ngbxval <= ngbxmax; ngbxval ++) {
        if ((ngbxval | ngbyval | ngbzval) != 0)   /* If not loop edge */
          DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++,
                                (posxval + dataptr->dimxval + ngbxval) % dataptr->dimxval,
                                (posyval + dataptr->dimyval + ngbyval) % dataptr->dimyval,
                                (poszval + dataptr->dimzval + ngbzval) % dataptr->dimzval);
      }
    }
  }

  return (edgelocnum);
}

/*
**
*/

static
Gnum
dgraphBuildGrid3Dvertex26T (
const DgraphBuildGrid3DData * restrict const  dataptr,
const Gnum                                    vertglbnum,
Gnum                                          edgelocnum,
const Gnum                                    posxval,
const Gnum                                    posyval,
const Gnum                                    poszval)
{
  Gnum                ngbxmin;
  Gnum                ngbxmax;
  Gnum                ngbxval;
  Gnum                ngbymin;
  Gnum                ngbymax;
  Gnum                ngbyval;
  Gnum                ngbzmin;
  Gnum                ngbzmax;
  Gnum                ngbzval;

  ngbxmin = dataptr->t26.ngbxmin;
  ngbxmax = dataptr->t26.ngbxmax;
  ngbymin = dataptr->t26.ngbymin;
  ngbymax = dataptr->t26.ngbymax;
  ngbzmin = dataptr->t26.ngbzmin;
  ngbzmax = dataptr->t26.ngbzmax;

  for (ngbzval = ngbzmin; ngbzval <= ngbzmax; ngbzval ++) {
    for (ngbyval = ngbymin; ngbyval <= ngbymax; ngbyval ++) {
      for (ngbxval = ngbxmin; ngbxval <= ngbxmax; ngbxval ++) {
        Gnum                vertglbend;

        vertglbend = (((poszval + ngbzval) % dataptr->dimzval)  * dataptr->dimyval +
                      ((posyval + ngbyval) % dataptr->dimyval)) * dataptr->dimxval +
                      ((posxval + ngbxval) % dataptr->dimxval)  + dataptr->baseval;
        if (vertglbend != vertglbnum) {           /* If not loop edge */
          if (dataptr->edloloctax != NULL)
            dataptr->edloloctax[edgelocnum] = ((vertglbend + vertglbnum) % 16) + 1;
          dataptr->edgeloctax[edgelocnum ++] = vertglbend;
        }
      }
    }
  }

  return (edgelocnum);
}

/*
**
*/

static
Gnum
dgraphBuildGrid3Dvertex6M (
const DgraphBuildGrid3DData * restrict const  dataptr,
const Gnum                                    vertglbnum,
Gnum                                          edgelocnum,
const Gnum                                    posxval,
const Gnum                                    posyval,
const Gnum                                    poszval)
{
  Gnum                ngbxval;
  Gnum                ngbyval;
  Gnum                ngbzval;

  ngbxval = posxval - 1;
  if (ngbxval >= 0)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, ngbxval, posyval, poszval);
  ngbxval = posxval + 1;
  if (ngbxval < dataptr->dimxval)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, ngbxval, posyval, poszval);

  ngbyval = posyval - 1;
  if (ngbyval >= 0)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, ngbyval, poszval);
  ngbyval = posyval + 1;
  if (ngbyval < dataptr->dimyval)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, ngbyval, poszval);

  ngbzval = poszval - 1;
  if (ngbzval >= 0)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, posyval, ngbzval);
  ngbzval = poszval + 1;
  if (ngbzval < dataptr->dimzval)
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, posyval, ngbzval);

  return (edgelocnum);
}

/*
**
*/

static
Gnum
dgraphBuildGrid3Dvertex6T (
const DgraphBuildGrid3DData * restrict const  dataptr,
const Gnum                                    vertglbnum,
Gnum                                          edgelocnum,
const Gnum                                    posxval,
const Gnum                                    posyval,
const Gnum                                    poszval)
{
  Gnum                ngbxval;
  Gnum                ngbyval;
  Gnum                ngbzval;

  if (dataptr->dimxval > 1) {
    ngbxval = (posxval + 1) % dataptr->dimxval;
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, ngbxval, posyval, poszval);
    if (dataptr->dimxval > 2) {
      ngbxval = (posxval + dataptr->dimxval - 1) % dataptr->dimxval;
      DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, ngbxval, posyval, poszval);
    }
  }

  if (dataptr->dimyval > 1) {
    ngbyval = (posyval + 1) % dataptr->dimyval;
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, ngbyval, poszval);
    if (dataptr->dimyval > 2) {
      ngbyval = (posyval + dataptr->dimyval - 1) % dataptr->dimyval;
      DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, ngbyval, poszval);
    }
  }

  if (dataptr->dimzval > 1) {
    ngbzval = (poszval + 1) % dataptr->dimzval;
    DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, posyval, ngbzval);
    if (dataptr->dimzval > 2) {
      ngbzval = (poszval + dataptr->dimzval - 1) % dataptr->dimzval;
      DGRAPHBUILDGRID3DNGB (dataptr, vertglbnum, edgelocnum ++, posxval, posyval, ngbzval);
    }
  }

  return (edgelocnum);
}

/*******************************************/
/*                                         */
/* The distributed graph building routine. */
/*                                         */
/*******************************************/

/* This routine builds a distrbuted grid graph
** of the given dimensions.
** hashval is the increment between two vertex
** indices (1 for sliced meshes).
** flagval is a combilation of:
** - 1  : 26-neighbor mesh (default: 6-neighbor mesh).
** - 2  : torus (default: mesh)
** - 4  : weighted vertices (default: no weights).
** - 8  : weighted edges (default: no weights).
** It returns:
** - 0   : graph created.
** - !0  : on error.
*/

int
dgraphBuildGrid3D (
Dgraph * restrict const     grafptr,              /* Graph            */
const Gnum                  baseval,              /* Base value       */
const Gnum                  dimxval,              /* First dimension  */
const Gnum                  dimyval,              /* Second dimension */
const Gnum                  dimzval,              /* Third dimension  */
const Gnum                  incrval,              /* Increment step   */
const int                   flagval)              /* Grid type        */
{
  DgraphBuildGrid3DData datadat;                  /* Data structure for creating vertices   */
  Gnum                  proclocadj;               /* Number of processes with most vertices */
  Gnum                  vertglbmin;               /* Minimum global index of local vertices */
  Gnum                  vertglbnbr;
  Gnum                  vertlocnbr;
  Gnum                  vertlocnnd;
  Gnum                  vertlocnum;
  Gnum *                vertloctax;
  Gnum                  velolocsiz;
  Gnum                  velolocsum;
  Gnum *                veloloctax;
  Gnum *                vlblloctax;
  Gnum                  vlbllocsiz;
  Gnum                  edgelocsiz;
  Gnum                  edgelocnum;
  Gnum *                edgeloctab;
  Gnum                  edlolocsiz;
  Gnum *                edloloctab;
  Gnum                  degrglbmax;

#ifdef SCOTCH_DEBUG_DGRAPH1
  if ((dimxval < 1) || (dimyval < 1) || (dimzval < 1)) { /* At least one vertex */
    errorPrint ("dgraphBuildGrid3D: invalid parameters (1)");
    return     (1);
  }
  if (incrval < 1) {
    errorPrint ("dgraphBuildGrid3D: invalid parameters (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  vertglbnbr = dimxval * dimyval * dimzval;
  vertlocnbr = DATASIZE (vertglbnbr, grafptr->procglbnbr, grafptr->proclocnum);

  if ((flagval & 1) != 0) {                       /* If 26-neighbor mesh */
    degrglbmax = 26;
    if ((flagval & 2) != 0) {                     /* If torus graph */
      datadat.t26.ngbxmin = (dimxval > 1) ? (dimxval - 1) : dimxval; /* Avoid loop edges */
      datadat.t26.ngbxmax = (dimxval > 2) ? (dimxval + 1) : dimxval;
      datadat.t26.ngbymin = (dimyval > 1) ? (dimyval - 1) : dimyval;
      datadat.t26.ngbymax = (dimyval > 2) ? (dimyval + 1) : dimyval;
      datadat.t26.ngbzmin = (dimzval > 1) ? (dimzval - 1) : dimzval;
      datadat.t26.ngbzmax = (dimzval > 2) ? (dimzval + 1) : dimzval;

      datadat.funcvrtptr = dgraphBuildGrid3Dvertex26T;
    }
    else
      datadat.funcvrtptr = dgraphBuildGrid3Dvertex26M;
  }
  else {                                          /* If 6-neighbor mesh */
    degrglbmax = 6;
    datadat.funcvrtptr = ((flagval & 2) != 0) ? dgraphBuildGrid3Dvertex6T : dgraphBuildGrid3Dvertex6M;
  }
  edgelocsiz  = vertlocnbr * degrglbmax;          /* (Possibly upper bound on) number of edges */
  vlbllocsiz  = (incrval != 1) ? vertlocnbr : 0;  /* If no hashing, no need for vertex labels  */
  velolocsiz  = ((flagval & 4) != 0) ? vertlocnbr : 0;
  edlolocsiz  = ((flagval & 8) != 0) ? edgelocsiz : 0;

  if (memAllocGroup ((void **) (void *)
                     &vertloctax, (size_t) ((vertlocnbr + 1) * sizeof (Gnum)), /* +1 to indicate end of array */
                     &veloloctax, (size_t) (velolocsiz       * sizeof (Gnum)),
                     &vlblloctax, (size_t) (vlbllocsiz       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphBuildGrid3D: out of memory (1)");
    return     (1);
  }
  if (memAllocGroup ((void **) (void *)
                     &edgeloctab, (size_t) (edgelocsiz * sizeof (Gnum)),
                     &edloloctab, (size_t) (edlolocsiz * sizeof (Gnum)), NULL) == NULL) {
    memFree    (vertloctax);
    errorPrint ("dgraphBuildGrid3D: out of memory (2)");
    return     (1);
  }

  datadat.baseval = baseval;
  datadat.dimxval = dimxval;
  datadat.dimyval = dimyval;
  datadat.dimzval = dimzval;
  datadat.edgeloctax = edgeloctab - baseval;
  datadat.edloloctax = ((flagval & 8) != 0) ? (edloloctab - baseval) : NULL;
  vertloctax = vertloctax - baseval;
  veloloctax = ((flagval & 4) != 0) ? (veloloctax - baseval) : NULL;
  vlblloctax = (incrval != 1) ? (vlblloctax - baseval) : NULL;

  proclocadj = vertglbnbr % grafptr->procglbnbr;  /* Number of processes with +1 number of vertices */
  vertglbmin = (vertglbnbr / grafptr->procglbnbr) * grafptr->proclocnum + MIN (grafptr->proclocnum, proclocadj);

  edgelocnum =
  vertlocnum = baseval;
  vertlocnnd = baseval + vertlocnbr;
  velolocsum = (veloloctax == NULL) ? vertlocnbr : 0;
  if (incrval != 1) {                             /* If strided or pseudo-randomly distributed mesh   */
    Gnum                vertglbidx;               /* Un-based global index of current vertex          */
    Gnum                rondlocnbr;               /* Number of already completed rounds of increments */
    Gnum                a;
    Gnum                b;

    a = (vertglbnbr > incrval) ? vertglbnbr : incrval; /* Get biggest of the two */
    b = (vertglbnbr + incrval) - a;               /* Get smallest of the two     */
    do {
      Gnum                t;

      t = a % b;
      if (t == 0)
        break;
      a = b;
      b = t;
    } while (b > 1);                              /* Compute GCD of vertglbnbr and incrval in b */

    rondlocnbr = (vertglbmin * b) / vertglbnbr;
    vertglbidx = (vertglbmin * incrval + rondlocnbr) % vertglbnbr; /* Compute skewed index, with rounds */

    for ( ; vertlocnum < vertlocnnd; vertlocnum ++) {
      Gnum                vertglbnum;
      Gnum                positmp;
      Gnum                posxval;
      Gnum                posyval;
      Gnum                poszval;

      poszval = vertglbidx / (dimxval * dimyval);
      positmp = vertglbidx % (dimxval * dimyval);
      posyval = positmp / dimxval;
      posxval = positmp % dimxval;

      vertglbnum = vertglbidx + baseval;
      vertloctax[vertlocnum] = edgelocnum;
      vlblloctax[vertlocnum] = vertglbnum;
      if (veloloctax != NULL) {
        velolocsum +=
        veloloctax[vertlocnum] = (vertglbnum % 16) + 1;
      }
      edgelocnum = datadat.funcvrtptr (&datadat, vertglbnum, edgelocnum, posxval, posyval, poszval);
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (edgelocnum > (edgelocsiz + baseval)) {
        errorPrint ("dgraphBuildGrid3D: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      vertglbidx = (vertglbidx + incrval) % vertglbnbr; /* Add increment to global index        */
      if (vertglbidx == rondlocnbr) {             /* If we looped back to the current beginning */
        rondlocnbr ++;                            /* Start a new round of increments            */
        vertglbidx = rondlocnbr;
      }
    }
  }
  else {                                          /* Regularly sliced mesh      */
    Gnum                vertglbnum;               /* Based global vertex number */
    Gnum                positmp;
    Gnum                posxval;
    Gnum                posyval;
    Gnum                poszval;

    poszval = vertglbmin / (dimxval * dimyval);
    positmp = vertglbmin % (dimxval * dimyval);
    posyval = positmp / dimxval;
    posxval = positmp % dimxval;

    for (vertglbnum = vertglbmin + baseval; vertlocnum < vertlocnnd; vertlocnum ++, vertglbnum ++) {
      vertloctax[vertlocnum] = edgelocnum;
      if (veloloctax != NULL) {
        velolocsum +=
        veloloctax[vertlocnum] = (vertglbnum % 16) + 1;
      }

      edgelocnum = datadat.funcvrtptr (&datadat, vertglbnum, edgelocnum, posxval, posyval, poszval);
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (edgelocnum > (edgelocsiz + baseval)) {
        errorPrint ("dgraphBuildGrid3D: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      if (++ posxval >= dimxval) {
        posxval = 0;
        if (++ posyval >= dimyval) {
          posyval = 0;
          poszval ++;
#ifdef SCOTCH_DEBUG_DGRAPH2
          if ((poszval >= dimzval) &&
              (vertglbnum < (vertglbnbr + baseval - 1))){
            errorPrint ("dgraphBuildGrid3D: internal error (X)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        }
      }
    }
  }
  vertloctax[vertlocnum] = edgelocnum;            /* Mark end of local vertex array */

  grafptr->flagval = (DGRAPHFREETABS | DGRAPHVERTGROUP | DGRAPHEDGEGROUP); /* All arrays will be freed on exit */

  if (dgraphBuild2 (grafptr, baseval,             /* Build distributed graph */
                    vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, veloloctax, velolocsum, NULL, vlblloctax,
                    edgelocnum - baseval, edgelocsiz, datadat.edgeloctax, NULL, datadat.edloloctax, degrglbmax) != 0) {
    memFree (datadat.edgeloctax + baseval);       /* Free memory group leaders */
    memFree (vertloctax + baseval);
    return  (1);
  }

  return (0);
}
