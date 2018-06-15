/* Copyright 2004,2007,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Part of a static mapper.                **/
/**                These lines are the data declarations   **/
/**                for the k-way graph mapping structures  **/
/**                and routines.                           **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 12 sep 1997     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     12 mar 1999     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     16 feb 2005     **/
/**                # Version 5.0  : from : 17 jun 2008     **/
/**                                 to     17 jun 2008     **/
/**                # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     31 aug 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     23 aug 2014     **/
/**                                                        **/
/**   NOTES      : # The comploadavg and comploaddlt       **/
/**                  should always be allocated together,  **/
/**                  with comploaddlt just after           **/
/**                  comploadavg.                          **/
/**                                                        **/
/**                # To date, the only k-way graph mapping **/
/**                  method is recursive bipartitioning,   **/
/**                  which does not use the comploadavg,   **/
/**                  comploaddlt and frontab parameters.   **/
/**                  Consequently, although allocated,     **/
/**                  these parameters are not set by this  **/
/**                  method. However, when coupled with    **/
/**                  the kdgraphMapSq() parallel static    **/
/**                  mapping routine, these parameters     **/
/**                  have to be computed. This is why the  **/
/**                  kgraphCost() routine is called within **/
/**                  kdgraphMapSq(). When other sequential **/
/**                  mapping routines are proposed, it may **/
/**                  be more interesting to move           **/
/**                  kgraphCost() to kgraphMapRb().        **/
/**                                                        **/
/**                # When (pfixtax != NULL), we are        **/
/**                  handling fixed vertices.              **/
/**                                                        **/
/**                # When (r.m.parttax != NULL), we are    **/
/**                  doing repartitioning.                 **/
/**                                                        **/
/************************************************************/

#define KGRAPH_H

/*
**  The defines.
*/

/*+ Graph option flags. +*/

#define KGRAPHFREEFRON              (GRAPHBITSNOTUSED)      /*+ Free frontier array              +*/
#define KGRAPHFREECOMP              (GRAPHBITSNOTUSED << 1) /*+ Free computational loads array   +*/
#define KGRAPHFREEPFIX              (GRAPHBITSNOTUSED << 2) /*+ Free fixed vertex array          +*/
#define KGRAPHFREEVMLO              (GRAPHBITSNOTUSED << 3) /*+ Free vertex migration cost array +*/
#define KGRAPHHASANCHORS            (GRAPHBITSNOTUSED << 4) /*+ The graph is a band graph        +*/

/*
**  The type and structure definitions.
*/

/*+ The graph structure. +*/

typedef struct Kgraph_ {
  Graph                     s;                    /*+ Current graph                                     +*/
  Arch                      a;                    /*+ Current architecture                              +*/
  Mapping                   m;                    /*+ Current mapping of graph vertices                 +*/
  struct {                                        /*+ Remapping structure                               +*/
    Mapping                 m;                    /*+ Old mapping                                       +*/
    Gnum                    crloval;              /*+ Coefficient load for regular edges                +*/
    Gnum                    cmloval;              /*+ Coefficient load for migration edges; may be zero +*/
    const Gnum *            vmlotax;              /*+ Vertex migration cost array                       +*/
  }                         r;
  Gnum                      vfixnbr;              /*+ Number of fixed vertices                          +*/
  const Anum *              pfixtax;              /*+ Fixed terminal part array                         +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices                       +*/
  Gnum *                    frontab;              /*+ Array of frontier vertex numbers                  +*/
  Gnum *                    comploadavg;          /*+ Array of target average loads                     +*/
  Gnum *                    comploaddlt;          /*+ Array of target imbalances                        +*/
  double                    comploadrat;          /*+ Ideal load balance per weight unit                +*/
  double                    kbalval;              /*+ Last k-way imbalance ratio                        +*/
  Gnum                      commload;             /*+ Communication load                                +*/
  Gnum                      levlnum;              /*+ Graph coarsening level                            +*/
} Kgraph;

/*+ The save graph structure. +*/

typedef struct KgraphStore_ {
  Gnum                      partnbr;              /*+ Number of parts                  +*/
  int                       mflaval;              /*+ Mapping properties               +*/
  Anum *                    parttab;              /*+ Mapping array [vertnbr]          +*/
  ArchDom *                 domntab;              /*+ Array of domains [termmax]       +*/
  Anum                      domnnbr;              /*+ Current number of domains        +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices      +*/
  Gnum *                    frontab;              /*+ Array of frontier vertex numbers +*/
  Gnum *                    comploadavg;          /*+ Array of target average loads    +*/
  Gnum *                    comploaddlt;          /*+ Array of target imbalances       +*/
  double                    kbalval;              /*+ Last k-way imbalance ratio       +*/
  Gnum                      commload;             /*+ Communication load               +*/
} KgraphStore;

/*
**  The function prototypes.
*/

#ifndef KGRAPH
#define static
#endif

int                         kgraphInit          (Kgraph * restrict const, const Graph * restrict const, const Arch * restrict const, const ArchDom * restrict const, const Gnum, const Anum * restrict const, const Anum * restrict const, const Gnum, const Gnum, const Gnum * restrict const);
void                        kgraphExit          (Kgraph * const);
void                        kgraphFrst          (Kgraph * const);
int                         kgraphCheck         (const Kgraph * const);
void                        kgraphCost          (Kgraph * const);
void                        kgraphFron          (Kgraph * const);
int                         kgraphBand          (Kgraph * restrict const, const Gnum, Kgraph * restrict const, Gnum * const, Gnum * restrict * restrict const);

int                         kgraphStoreInit     (const Kgraph * const, KgraphStore * const);
void                        kgraphStoreExit     (KgraphStore * const);
void                        kgraphStoreSave     (const Kgraph * const, KgraphStore * const);
void                        kgraphStoreUpdt     (Kgraph * const, const KgraphStore * const);

#undef static
