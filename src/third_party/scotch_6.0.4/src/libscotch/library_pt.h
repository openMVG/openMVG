/* Copyright 2004,2007-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_pt.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Declaration file for the LibPtscotch    **/
/**                parallel static mapping and sparse      **/
/**                matrix block ordering library.          **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     22 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     31 may 1999     **/
/**                # Version 3.4  : from : 10 oct 1999     **/
/**                                 to     15 nov 2001     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     20 dec 2005     **/
/**                # Version 5.0  : from : 26 apr 2006     **/
/**                                 to   : 20 feb 2008     **/
/**                # Version 5.1  : from : 30 nov 2007     **/
/**                                 to   : 07 aug 2011     **/
/**                # Version 6.0  : from : 12 sep 2008     **/
/**                                 to     01 dec 2012     **/
/**                                                        **/
/************************************************************/

#ifndef PTSCOTCH_H
#define PTSCOTCH_H

/*
**  The defines and includes.
*/

#ifndef SCOTCH_H
#include "scotch.h"
#endif /* SCOTCH_H */

/*
**  The type and structure definitions.
*/

/*+ Parallel processing flag. +*/

#ifndef SCOTCH_PTSCOTCH
#define SCOTCH_DUMMYPTFLAG
#endif /* SCOTCH_PTSCOTCH */

/*+ Opaque objects. The dummy sizes of these
objects, computed at compile-time by program
"dummysizes", are given as double values for
proper padding                               +*/

typedef struct {
  double                    dummy[DUMMYSIZEDGRAPH];
} SCOTCH_Dgraph;

typedef struct {
  double                    dummy[DUMMYSIZEDGRAPHHALOREQ];
} SCOTCH_DgraphHaloReq;

typedef struct {
  double                    dummy[DUMMYSIZEDMAP];
} SCOTCH_Dmapping;

typedef struct {
  double                    dummy[DUMMYSIZEDORDER];
} SCOTCH_Dordering;

/*
**  The function prototypes.
*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

SCOTCH_Dgraph *             SCOTCH_dgraphAlloc  (void);
int                         SCOTCH_dgraphInit   (SCOTCH_Dgraph * const, MPI_Comm);
void                        SCOTCH_dgraphExit   (SCOTCH_Dgraph * const);
void                        SCOTCH_dgraphFree   (SCOTCH_Dgraph * const);
int                         SCOTCH_dgraphLoad   (SCOTCH_Dgraph * const, FILE * const, const SCOTCH_Num, const SCOTCH_Num);
int                         SCOTCH_dgraphSave   (SCOTCH_Dgraph * const, FILE * const);
int                         SCOTCH_dgraphCheck  (const SCOTCH_Dgraph * const);
int                         SCOTCH_dgraphBand   (SCOTCH_Dgraph * const, const SCOTCH_Num, SCOTCH_Num * const, const SCOTCH_Num, SCOTCH_Dgraph * const);
int                         SCOTCH_dgraphBuild  (SCOTCH_Dgraph * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
int                         SCOTCH_dgraphBuildGrid3D (SCOTCH_Dgraph * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const int);
int                         SCOTCH_dgraphCoarsen (SCOTCH_Dgraph * const, const SCOTCH_Num, const double, const SCOTCH_Num, SCOTCH_Dgraph * const, SCOTCH_Num * const);
int                         SCOTCH_dgraphGather (const SCOTCH_Dgraph * const, SCOTCH_Graph * const);
int                         SCOTCH_dgraphGrow   (SCOTCH_Dgraph * const, const SCOTCH_Num, SCOTCH_Num * const, const SCOTCH_Num, SCOTCH_Num * const);
int                         SCOTCH_dgraphInducePart (SCOTCH_Dgraph * const, const SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num, SCOTCH_Dgraph * const);
int                         SCOTCH_dgraphScatter (SCOTCH_Dgraph * const, const SCOTCH_Graph * const);
int                         SCOTCH_dgraphRedist (SCOTCH_Dgraph * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num, SCOTCH_Dgraph * const);
void                        SCOTCH_dgraphSize   (const SCOTCH_Dgraph * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
void                        SCOTCH_dgraphData   (const SCOTCH_Dgraph * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num ** const, SCOTCH_Num ** const, SCOTCH_Num ** const, MPI_Comm * const);
int                         SCOTCH_dgraphStat   (const SCOTCH_Dgraph * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, double * const, double * const);
int                         SCOTCH_dgraphGhst   (SCOTCH_Dgraph * const);
int                         SCOTCH_dgraphHalo   (SCOTCH_Dgraph * const, void * const, const MPI_Datatype);
int                         SCOTCH_dgraphHaloAsync (SCOTCH_Dgraph * const, void * const, const MPI_Datatype, SCOTCH_DgraphHaloReq * const);
SCOTCH_DgraphHaloReq *      SCOTCH_dgraphHaloReqAlloc (void);
int                         SCOTCH_dgraphHaloWait (SCOTCH_DgraphHaloReq * const);
int                         SCOTCH_dgraphMapInit (const SCOTCH_Dgraph * const, SCOTCH_Dmapping * const, const SCOTCH_Arch * const, SCOTCH_Num * const);
void                        SCOTCH_dgraphMapExit (const SCOTCH_Dgraph * const, SCOTCH_Dmapping * const);
int                         SCOTCH_dgraphMapSave (const SCOTCH_Dgraph * const, const SCOTCH_Dmapping * const, FILE * const);
int                         SCOTCH_dgraphMapView (SCOTCH_Dgraph * const, const SCOTCH_Dmapping * const, FILE * const);
int                         SCOTCH_dgraphMapCompute (SCOTCH_Dgraph * const, SCOTCH_Dmapping * const, SCOTCH_Strat * const);
int                         SCOTCH_dgraphMap     (SCOTCH_Dgraph * const, const SCOTCH_Arch * const, SCOTCH_Strat * const, SCOTCH_Num * const);
int                         SCOTCH_dgraphPart    (SCOTCH_Dgraph * const, const SCOTCH_Num, SCOTCH_Strat * const, SCOTCH_Num * const);
int                         SCOTCH_dgraphCorderInit (const SCOTCH_Dgraph * const, SCOTCH_Ordering * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const);
void                        SCOTCH_dgraphCorderExit (const SCOTCH_Dgraph * const, SCOTCH_Ordering * const);

int                         SCOTCH_dgraphOrderInit (const SCOTCH_Dgraph * const, SCOTCH_Dordering * const);
void                        SCOTCH_dgraphOrderExit (const SCOTCH_Dgraph * const, SCOTCH_Dordering * const);
int                         SCOTCH_dgraphOrderSave (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, FILE * const);
int                         SCOTCH_dgraphOrderSaveBlock (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, FILE * const);
int                         SCOTCH_dgraphOrderSaveMap (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, FILE * const);
int                         SCOTCH_dgraphOrderSaveTree (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, FILE * const);
int                         SCOTCH_dgraphOrderPerm (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, SCOTCH_Num * const);
SCOTCH_Num                  SCOTCH_dgraphOrderCblkDist (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const);
int                         SCOTCH_dgraphOrderTreeDist (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, SCOTCH_Num * const, SCOTCH_Num * const);
int                         SCOTCH_dgraphOrderCompute (SCOTCH_Dgraph * const, SCOTCH_Dordering * const, SCOTCH_Strat * const);
int                         SCOTCH_dgraphOrderComputeList (SCOTCH_Dgraph * const, SCOTCH_Dordering * const, const SCOTCH_Num, const SCOTCH_Num * const, SCOTCH_Strat * const);
int                         SCOTCH_dgraphOrderGather (const SCOTCH_Dgraph * const, const SCOTCH_Dordering * const, SCOTCH_Ordering * const);

SCOTCH_Dmapping *           SCOTCH_dmapAlloc    (void);

SCOTCH_Dordering *          SCOTCH_dorderAlloc  (void);

int                         SCOTCH_stratDgraphMap (SCOTCH_Strat * const, const char * const);
int                         SCOTCH_stratDgraphMapBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const double);
int                         SCOTCH_stratDgraphClusterBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const double, const double);
int                         SCOTCH_stratDgraphOrder (SCOTCH_Strat * const, const char * const);
int                         SCOTCH_stratDgraphOrderBuild (SCOTCH_Strat * const, const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num, const double);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* PTSCOTCH_H */
