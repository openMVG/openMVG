/* Copyright 2008,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kdgraph.h                               **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a static mapper.                **/
/**                These lines are the data declarations   **/
/**                for the parallel k-way graph partiton-  **/
/**                ing structures and routines.            **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 31 mar 2008     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 29 aug 2012     **/
/**                                 to     29 aug 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The dynamic mapping structure. +*/

typedef struct Kdmapping_ {
  Dmapping *                mappptr;              /*+ Resulting mapping +*/
  ArchDom                   domnorg;              /*+ Initial domain    +*/
} Kdmapping;

/*+ The graph structure. +*/

typedef struct Kdgraph_ {
  Dgraph                    s;                    /*+ Source graph +*/
  Kdmapping                 m;                    /*+ Mapping      +*/
  INT                       levlnum;
} Kdgraph;

/*
**  The function prototypes.
*/

#ifndef KDGRAPH
#define static
#endif

int                         kdgraphInit         (Kdgraph * const, const Dgraph * restrict const, Dmapping * restrict const);
void                        kdgraphExit         (Kdgraph * const);
int                         kdgraphFold         (const Kdgraph *, const int, Kdgraph * const);
int                         kdgraphFold2        (const Kdgraph *, const int, Kdgraph * const, MPI_Comm);
#ifdef KGRAPH_H
int                         kdgraphGather       (Kdgraph *, Kgraph *);
#endif /* KGRAPH_H */

#undef static
