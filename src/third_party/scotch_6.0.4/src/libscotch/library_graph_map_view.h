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
/**   NAME       : library_graph_map_view.h                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the declarations for    **/
/**                the mapping viewing routines.           **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 04 feb 2007     **/
/**                                 to     04 feb 2007     **/
/**                                                        **/
/************************************************************/

#define LIBRARY_GRAPH_MAP_VIEW_H

/*
**  The type definitions.
*/

/*+ Complementary vertex structure. +*/

typedef struct GraphMapViewVertex_ {
  Gnum                      passnum;              /*+ Number of pass when vertex selected   +*/
  Gnum                      vertdist;             /*+ Current distance from diameter vertex +*/
} GraphMapViewVertex;

/*+ Neighbor queue. +*/

typedef struct GraphMapViewQueue_ {
  Gnum *                    head;                 /*+ Head of distance queue  +*/
  Gnum *                    tail;                 /*+ Tail of distance queue  +*/
  Gnum *                    qtab;                 /*+ Array of queue elements +*/
} GraphMapViewQueue;

/*
**  The function prototypes.
*/

#ifndef LIBRARY_GRAPH_MAP_VIEW
#define static
#endif

static Gnum                 graphMapView3       (const Graph * const, const Anum * const, const Anum);

#undef static

/*
**  The macro definitions.
*/

#define graphMapViewQueueFlush(queue) ((queue)->head = (queue)->tail = (queue)->qtab)
#define graphMapViewQueueEmpty(queue) ((queue)->head <= (queue)->tail)
#define graphMapViewQueuePut(queue,vnum) (* ((queue)->head ++) = (vnum))
#define graphMapViewQueueGet(queue) (* ((queue)->tail ++))
