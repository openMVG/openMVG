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
/**   NAME       : hdgraph.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source halo distributed graph   **/
/**                structure.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 15 apr 2006     **/
/**                                 to     16 jun 2007     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to     04 nov 2010     **/
/**                                                        **/
/************************************************************/

#define HDGRAPH_H

/*
**  The defines.
*/

/*+ Graph option flags. +*/

#define HDGRAPHFREEVHND             0x0400        /* Free vnhdtab array */
#define HDGRAPHFREETABS             (DGRAPHFREETABS | HGRAPHFREEVHND)

/*
**  The type and structure definitions.
*/

/*+ Halo distributed graph structure. In order to
    keep efficiency, distributed halo graphs are not
    considered as regular graphs as sequential halo
    graphs were. Halo distributed graphs have a compact
    vertex array, with halo edges added at the end of
    each vertex sub-array. They are not visible when
    considering the vertlocnbr, vertloctax (which is in
    fact most often of size vhallocnbr + 1 when the graph
    is compact, as in this case we have vnhdloctax =
    vertloctax + 1) and vendloctax (which is of size
    vertlocnbr) of the embedded distributed graph, but
    can be accessed through vendloctax and vnhdloctax.
    Halo vertex ends are stored only in edgeloctax, not
    in edgegsttax, except when graph has only an edgegsttax
    and no edgeloctax. Since halo vertices have no real
    existence in distributed graphs, they are simply
    numbered from baseval. They are converted into real
    vertices when a distributed halo graph is turned into
    a sequential halo graph.                                */

typedef struct Hdgraph_ {
  Dgraph                    s;                    /*+ Source distributed graph                       +*/
  Gnum                      vhallocnbr;           /*+ Local number of halo end vertices              +*/
  Gnum *                    vhndloctax;           /*+ End vertex array including halo vertex indices +*/
  Gnum                      ehallocnbr;           /*+ Local number of halo edges                     +*/
  Gnum                      levlnum;              /*+ Nested dissection level                        +*/
} Hdgraph;

/*
**  The function prototypes.
*/

#ifndef HDGRAPH
#define static
#endif

int                         hdgraphInit         (Hdgraph * const);
void                        hdgraphExit         (Hdgraph * const);
void                        hdgraphFree         (Hdgraph * const);
int                         hdgraphFold         (const Hdgraph *, const int, Hdgraph * const);
int                         hdgraphFold2        (const Hdgraph *, const int, Hdgraph * const, MPI_Comm);
int                         hdgraphCheck        (const Hdgraph *);
#ifdef HGRAPH_H
int                         hdgraphGather       (Hdgraph *, Hgraph *);
#endif /* HGRAPH_H */
int                         hdgraphInduceList   (Hdgraph * restrict const, const Gnum, const Gnum * restrict const, Hdgraph * restrict const);

#undef static
