/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : mesh_coarsen.h                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source mesh coarsening          **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 oct 2003     **/
/**                                 to     04 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/** Prime number for cache-friendly perturbations. **/

#define MESHCOARSENPERTPRIME        31            /* Prime number */

/** Prime number for hashing vertex numbers. **/

#define MESHCOARSENHASHPRIME        17            /* Prime number */

/*
**  The type and structure definitions.
*/

/*+ Here are the edge matching function types for coarsening. +*/

typedef enum MeshCoarsenType_ {
  MESHCOARSENNGB,                                 /*+ Most neighbors matching   +*/
  MESHCOARSENNBR                                  /*+ Number of matching types  +*/
} MeshCoarsenType;

/*+ A table made of such elements is used during
    coarsening to build the edge array of the new
    mesh, after the labeling of the vertices.     +*/

typedef struct MeshCoarsenMult_ {
  Gnum                      finevelmnum[2];
} MeshCoarsenMult;

/*+ A table made of such cells is used during
    coarsening to build the edge array of the
    elements of the new mesh.                 +*/

typedef struct MeshCoarsenHngb_ {
  Gnum                      coarvelmnum;          /*+ Coarse origin element vertex (i.e. pass) number +*/
  Gnum                      coarvnodnum;          /*+ Neighbor fine node vertex number                +*/
} MeshCoarsenHngb;

/*+ A table made of such cells is used during
    coarsening to detect and merge bridge nodes,
    that is, nodes that connect only two coarse
    nodes together.                              +*/

typedef struct MeshCoarsenHbdg_ {
  Gnum                      coarvelmnum;          /*+ Coarse origin element vertex (i.e. pass) number +*/
  Gnum                      coarvelmend;          /*+ Coarse end element vertex number                +*/
  Gnum                      coarvnodnum;          /*+ Number of connecting coarse node                +*/
} MeshCoarsenHbdg;

/*+ A table made of such elements is used during
    coarsening to build the edge array of the new
    mesh, after the labeling of the vertices.     +*/

typedef struct MeshCoarsenNgHash_ {
  Gnum                      velmnum;              /*+ Origin element vertex (i.e. pass) number   +*/
  Gnum                      velmend;              /*+ End element vertex number in fine mesh     +*/
  Gnum                      vnngnbr;              /*+ Number of shared neighboring node vertices +*/
  Gnum                      vnbgnbr;              /*+ Number of bridge neighboring node vertices +*/
} MeshCoarsenNgHash;

/*
**  The function prototypes.
*/

#ifndef MESH_COARSEN
#define static
#endif

int                         meshCoarsen         (const Mesh * restrict const, Mesh * restrict const, Gnum * restrict * const, const Gnum, const double, const MeshCoarsenType);

static void                 meshCoarsenMatchNg  (const Mesh * restrict const, MeshCoarsenMult * restrict const, Gnum * restrict const, Gnum * restrict const, Gnum * restrict const, Gnum * restrict const);

#undef static
