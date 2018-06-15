/* Copyright 2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph.h                               **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for distributed edge bipartition- **/
/**                ing routines.                           **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 10 sep 2007     **/
/**                                 to   : 14 apr 2011     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to   : 31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Graph option flags. +*/

#define BDGRAPHFREEFRON             (DGRAPHBITSNOTUSED) /* Free part array               */
#define BDGRAPHFREEPART             (DGRAPHBITSNOTUSED << 1) /* Free frontier array      */
#define BDGRAPHFREEVEEX             (DGRAPHBITSNOTUSED << 2) /* Free external gain array */

/*+ Active graph structure. +*/

typedef struct Bdgraph_ {
  Dgraph                    s;                    /*+ Distributed source graph                           +*/
  Gnum *                    veexloctax;           /*+ Local vertex external gain array if moved to 1     +*/
  Gnum                      veexglbsum;           /*+ Global sum of veexloctax array cells               +*/
  GraphPart *               partgsttax;           /*+ Based local part array: 0,1: part                  +*/
  Gnum *                    fronloctab;           /*+ Array of local frontier vertex numbers             +*/
  Gnum                      fronlocnbr;           /*+ Number of local frontier vertices                  +*/
  Gnum                      fronglbnbr;           /*+ Number of global frontier vertices                 +*/
  Gnum                      complocload0;         /*+ Local load of part 0                               +*/
  Gnum                      compglbload0;         /*+ Global load of part 0                              +*/
  Gnum                      compglbload0min;      /*+ Minimum allowed load in part 0 (strategy variable) +*/
  Gnum                      compglbload0max;      /*+ Maximum allowed load in part 0 (strategy variable) +*/
  Gnum                      compglbload0avg;      /*+ Global average load of part 0                      +*/
  Gnum                      compglbload0dlt;      /*+ Load difference from the average                   +*/
  Gnum                      complocsize0;         /*+ Number of local vertices in part 0                 +*/
  Gnum                      compglbsize0;         /*+ Number of global vertices in part 0                +*/
  Gnum                      commglbload;          /*+ Global communication load                          +*/
  Gnum                      commglbgainextn;      /*+ Global external gain if all swapped                +*/
  Gnum                      commglbloadextn0;     /*+ Global communication load if all moved to part 0   +*/
  Gnum                      commglbgainextn0;     /*+ Global external gain if all swapped from part 0    +*/
  double                    bbalglbval;           /*+ Bipartitioning imbalance ratio (strategy variable) +*/
  Anum                      domndist;             /*+ Distance between subdomains                        +*/
  Anum                      domnwght[2];          /*+ Weights for each subdomain                         +*/
  INT                       levlnum;              /*+ Graph coarsening level                             +*/
} Bdgraph;

/*+ The distributed save graph structure. +*/

typedef struct BdgraphStore_ {
  Gnum                      fronlocnbr;           /*+ Number of local frontier vertices   +*/
  Gnum                      fronglbnbr;           /*+ Number of frontier vertices         +*/
  Gnum                      complocload0;         /*+ Local load in part 0                +*/
  Gnum                      compglbload0;         /*+ Load in part 0                      +*/ 
  Gnum                      compglbload0dlt;      /*+ Difference from the average         +*/
  Gnum                      complocsize0;         /*+ Number of local vertices in part 0  +*/
  Gnum                      compglbsize0;         /*+ Number of global vertices in part 0 +*/
  Gnum                      commglbload;
  Gnum                      commglbgainextn;
  byte *                    datatab;              /*+ Variable-sized data array           +*/
} BdgraphStore;

/*
**  The function prototypes.
*/

#ifndef BDGRAPH
#define static
#endif

#ifdef DMAPPING_H
int                         bdgraphInit         (Bdgraph * const, const Dgraph * const, const Dgraph * const, const Arch * const, const ArchDom[]);
#endif /* DMAPPING_H */
void                        bdgraphInit2        (Bdgraph * const, const Anum, const Anum, const Anum);
#ifdef DMAPPING_H
int                         bdgraphInit3        (Bdgraph * const, const Dgraph * const, const Dmapping * const, const ArchDom[]);
#endif /* DMAPPING_H */
void                        bdgraphExit         (Bdgraph * restrict const);
void                        bdgraphFree         (Bdgraph * restrict const);
void                        bdgraphZero         (Bdgraph * restrict const);
int                         bdgraphCheck        (const Bdgraph * restrict const);
#ifdef BGRAPH_H
int                         bdgraphGatherAll    (const Bdgraph * restrict const, Bgraph * restrict);
#endif /* BGRAPH_H */

int                         bdgraphStoreInit    (const Bdgraph * const, BdgraphStore * const);
void                        bdgraphStoreExit    (BdgraphStore * const);
void                        bdgraphStoreSave    (const Bdgraph * const , BdgraphStore * const);
void                        bdgraphStoreUpdt    (Bdgraph * const, const BdgraphStore * const);

#undef static
