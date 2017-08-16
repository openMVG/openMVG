/* Copyright 2008,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kdgraph_map_rb_part.h                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the Parallel Dual Recursive         **/
/**                Bipartitioning mapping algorithm.       **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 23 jun 2008     **/
/**                                 to     31 aug 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds folded graph data, whether centralized or distributed +*/

typedef struct KdgraphMapRbPartGraph_ {
  ArchDom                   domnorg;              /*+ Domain to bipartition at this stage +*/
  int                       procnbr;              /*+ Number of processes holding graph   +*/
  INT                       levlnum;              /*+ Level number                        +*/
  union {
    Graph                   cgrfdat;              /*+ Centralized graph +*/
    Dgraph                  dgrfdat;              /*+ Distributed graph +*/
  } data;
} KdgraphMapRbPartGraph;

/*+ This structure holds the data passed to the subgraph building threads. +*/

typedef struct KdgraphMapRbPartThread_ {
  Dmapping *                mappptr;              /*+ Pointer to mapping structure                   +*/
  Dgraph *                  orggrafptr;           /*+ Pointer to original graph                      +*/
  const ArchDom *           inddomnptr;           /*+ Pointer to subjob domain                       +*/
  Gnum                      indvertnbr;           /*+ Local number of vertices in subgraph           +*/
  GraphPart                 indpartval;           /*+ Graph part from which to extract subgraph      +*/
  GraphPart *               indparttax;           /*+ Based local vertex partition flags in subgraph +*/
  KdgraphMapRbPartGraph *   fldgrafptr;           /*+ Pointer to folded graph union area             +*/
  int                       fldpartval;           /*+ Part of processor array to which to fold to    +*/
  int                       fldprocnbr;           /*+ Number of processes in folded communicator     +*/
  int                       fldprocnum;           /*+ Rank of process in folded communicator, or -1  +*/
  MPI_Comm                  fldproccomm;          /*+ Communicator for the folded graph, if any      +*/
} KdgraphMapRbPartThread;

/*+ This structure holds the data passed to each bipartitioning job. +*/

typedef struct KdgraphMapRbPartData_ {
  Dmapping *                mappptr;
  const KdgraphMapRbParam * paraptr;
  double                    comploadrat;          /*+ Ideal vertex load per target load   +*/
  double                    comploadmin;          /*+ Minimum vertex load per target load +*/
  double                    comploadmax;          /*+ Maximum vertex load per target load +*/
} KdgraphMapRbPartData;

/*
**  The function prototypes.
*/

#ifndef KDGRAPH_MAP_RB
#define static
#endif

int                         kdgraphMapRbPart    (Kdgraph * const, Kdmapping * const, const KdgraphMapRbParam * const);

#undef static
