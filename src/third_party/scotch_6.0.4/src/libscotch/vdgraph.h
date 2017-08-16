/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for distributed vertex separation **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 06 feb 2006     **/
/**                                 to   : 29 apr 2006     **/
/**                # Version 5.1  : from : 07 nov 2007     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Active graph structure. +*/

typedef struct Vdgraph_ {
  Dgraph                    s;                    /*+ Source distributed graph                                                 +*/
  GraphPart *               partgsttax;           /*+ Based local part array: 0,1: part; 2: separator                          +*/
  Gnum                      compglbloaddlt;       /*+ Load difference between both parts                                       +*/
  Gnum                      compglbload[3];       /*+ Global loads of both parts and of separator; TRICK: before compglbsize[] +*/
  Gnum                      compglbsize[3];       /*+ Number of vertices in parts; compglbsize[2] is fronglbnbr, the separator +*/
  Gnum                      complocload[3];       /*+ Local loads of both parts and of separator; TRICK: before complocsize[]  +*/
  Gnum                      complocsize[3];       /*+ Number of vertices in parts; complocsize[2] is fronlocnbr, the separator +*/
  Gnum *                    fronloctab;           /*+ Array of local frontier vertex numbers                                   +*/
  Gnum                      levlnum;              /*+ Nested dissection or coarsening level                                    +*/
} Vdgraph;

/*+ The graph separator storing structure. +*/

typedef struct VdgraphStore_ {
  Gnum                      fronglbnbr;           /*+ Number of frontier nodes          +*/
  Gnum                      compglbloaddlt;       /*+ Difference from the average       +*/
  Gnum                      compglbload[2];       /*+ Load in both parts                +*/
  Gnum                      compglbsize0;         /*+ Number of vertices in part 0      +*/
  Gnum                      complocsize0;         /*+ Number of vertices in parts       +*/
  Gnum                      fronlocnbr;           /*+ Number of local frontier vertices +*/
  byte *                    datatab;              /*+ Variable-sized data array         +*/
} VdgraphStore;

/*
**  The function prototypes.
*/

#ifndef VDGRAPH
#define static
#endif

int                         vdgraphInit         (Vdgraph * restrict const, MPI_Comm);
void                        vdgraphExit         (Vdgraph * const);
void                        vdgraphZero         (Vdgraph * const);
int                         vdgraphCheck        (const Vdgraph * const);
#ifdef VGRAPH_H
int                         vdgraphGatherAll    (const Vdgraph * restrict const, Vgraph * restrict);
#endif /* VGRAPH_H */

int                         vdgraphStoreInit    (const Vdgraph * const, VdgraphStore * const);
void                        vdgraphStoreExit    (VdgraphStore * const);
void                        vdgraphStoreSave    (const Vdgraph * const , VdgraphStore * const);
void                        vdgraphStoreUpdt    (Vdgraph * const, const VdgraphStore * const);

#undef static
