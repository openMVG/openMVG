/* Copyright 2004,2007,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source halo graph structure.    **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 02 jan 2002     **/
/**                                 to     30 apr 2004     **/
/**                # Version 5.0  : from : 19 dec 2006     **/
/**                                 to     19 dec 2006     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to     17 oct 2012     **/
/**                                                        **/
/************************************************************/

#define HGRAPH_H

/*
**  The defines.
*/

/*+ Graph option flags. +*/

#define HGRAPHFREEVNHD              0x0400        /* Free vnhdtab array */
#define HGRAPHFREETABS              (GRAPHFREETABS | HGRAPHFREEVNHD)

/*
**  The type and structure definitions.
*/

/*+ Halo graph structure. +*/

typedef struct Hgraph_ {
  Graph                     s;                    /*+ Source graph                                                   +*/
  Gnum                      vnohnbr;              /*+ Number of non-halo vertices                                    +*/
  Gnum                      vnohnnd;              /*+ Based number of first halo vertex in graph (s.vertnnd if none) +*/
  Gnum *                    vnhdtax;              /*+ End vertex array for non-halo vertices [vnohnbr, based]        +*/
  Gnum                      vnlosum;              /*+ Sum of vertex loads for non-halo vertices only (<= s.velosum)  +*/
  Gnum                      enohnbr;              /*+ Number of non-halo edges                                       +*/
  Gnum                      enohsum;              /*+ Sum of non-halo edge loads                                     +*/
  Gnum                      levlnum;              /*+ Nested dissection level                                        +*/
} Hgraph;

/*
**  The function prototypes.
*/

#ifndef HGRAPH
#define static
#endif

int                         hgraphInit          (Hgraph * const);
void                        hgraphExit          (Hgraph * const);
void                        hgraphFree          (Hgraph * const);
Gnum                        hgraphBase          (Hgraph * const, const Gnum);
int                         hgraphCheck         (const Hgraph *);
int                         hgraphInduceList    (const Hgraph * const, const VertList * const, const Gnum, Hgraph * const);
void                        hgraphUnhalo        (const Hgraph * const, Graph * const);

#undef static
