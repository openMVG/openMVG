/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_build_grid3d.h                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the data             **/
/**                declarations for the distributed 3D     **/
/**                grid graph building routine.            **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 06 jun 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
** The type and structure definitions.
*/

/*+ The multinode table element, which contains
    pairs of based indices of collapsed vertices.
    Both values are equal for uncollapsed vertices. +*/

typedef struct DgraphBuildGrid3DData_ {
  Gnum                      baseval;
  Gnum                      dimxval;
  Gnum                      dimyval;
  Gnum                      dimzval;
  Gnum *                    edgeloctax;
  Gnum *                    edloloctax;
  Gnum                   (* funcvrtptr) (const struct DgraphBuildGrid3DData_ * restrict const, const Gnum, Gnum, const Gnum, const Gnum, const Gnum);
  struct {                                        /* Pre-computed data for 26-neighbor torus */
    Gnum                    ngbxmin;
    Gnum                    ngbxmax;
    Gnum                    ngbymin;
    Gnum                    ngbymax;
    Gnum                    ngbzmin;
    Gnum                    ngbzmax;
  } t26;
} DgraphBuildGrid3DData;

/*
** The function prototypes.
*/

#ifndef DGRAPH_BUILD_GRID3D
#define static
#endif

static Gnum                 dgraphBuildGrid3Dvertex26M (const DgraphBuildGrid3DData * restrict const, const Gnum, Gnum, const Gnum, const Gnum, const Gnum);
static Gnum                 dgraphBuildGrid3Dvertex26T (const DgraphBuildGrid3DData * restrict const, const Gnum, Gnum, const Gnum, const Gnum, const Gnum);
static Gnum                 dgraphBuildGrid3Dvertex6M (const DgraphBuildGrid3DData * restrict const, const Gnum, Gnum, const Gnum, const Gnum, const Gnum);
static Gnum                 dgraphBuildGrid3Dvertex6T (const DgraphBuildGrid3DData * restrict const, const Gnum, Gnum, const Gnum, const Gnum, const Gnum);

#undef static

/*
** The macro definitions.
*/

#define DGRAPHBUILDGRID3DNGB(d,v,e,x,y,z) {                                                                        \
                                      Gnum                edgeloctmp;                                              \
                                      Gnum                vertglbend;                                              \
                                      edgeloctmp = (e);                                                            \
                                      vertglbend = ((z) * (d)->dimyval + (y)) * (d)->dimxval + (x) + (d)->baseval; \
                                      (d)->edgeloctax[edgeloctmp] = vertglbend;                                    \
                                      if ((d)->edloloctax != NULL)                                                 \
                                        (d)->edloloctax[edgeloctmp] = ((vertglbend + (v)) % 16) + 1;               \
                                    }

