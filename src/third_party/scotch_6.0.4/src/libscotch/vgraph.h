/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for vertex separation routines.   **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 aug 1996     **/
/**                                 to   : 17 oct 1997     **/
/**                # Version 3.3  : from : 13 mar 1999     **/
/**                                 to   : 13 mar 1999     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to   : 07 jan 2002     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

#define VGRAPH_H

/*
**  The type and structure definitions.
*/

/*+ Active graph structure. +*/

typedef struct Vgraph_ {
  Graph                     s;                    /*+ Source graph                                       +*/
  GraphPart *               parttax;              /*+ Based part array: 0,1: part; 2: separator          +*/
  Gnum                      compload[3];          /*+ Size of both parts and separator                   +*/
  Gnum                      comploaddlt;          /*+ Load difference between both parts                 +*/
  Gnum                      compsize[2];          /*+ Number of vertices in parts (separator is fronnbr) +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices; TRICK: compsize[2]    +*/
  Gnum *                    frontab;              /*+ Array of frontier vertex numbers                   +*/
  Gnum                      levlnum;              /*+ Nested dissection or coarsening level              +*/
} Vgraph;

/*+ The graph separator storing structure. +*/

typedef struct VgraphStore_ {
  Gnum                      fronnbr;              /*+ Number of frontier nodes     +*/
  Gnum                      comploaddlt;          /*+ Difference from the average  +*/
  Gnum                      compload[2];          /*+ Load in both parts           +*/
  Gnum                      compsize0;            /*+ Number of vertices in part 0 +*/
  byte *                    datatab;              /*+ Variable-sized data array    +*/
} VgraphStore;

/*
**  The function prototypes.
*/

#ifndef VGRAPH
#define static
#endif

void                        vgraphExit          (Vgraph * const);
void                        vgraphZero          (Vgraph * const);
int                         vgraphCheck         (const Vgraph * const);

int                         vgraphStoreInit     (const Vgraph * const, VgraphStore * const);
void                        vgraphStoreExit     (VgraphStore * const);
void                        vgraphStoreSave     (const Vgraph * const , VgraphStore * const);
void                        vgraphStoreUpdt     (Vgraph * const, const VgraphStore * const);

#undef static
