/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb_map.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the Dual Recursive Bipartitioning   **/
/**                mapping algorithm.                      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 23 mar 1993     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     09 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     04 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 15 nov 1995     **/
/**                                 to     15 nov 1995     **/
/**                # Version 3.2  : from : 01 oct 1996     **/
/**                                 to     10 jun 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     17 may 1999     **/
/**                # Version 3.4  : from : 12 sep 2001     **/
/**                                 to     06 nov 2001     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     05 may 2006     **/
/**                # Version 5.1  : from : 30 sep 2008     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     28 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Dual recursive bipartitioning option flags. +*/

#define KGRAPHMAPRBMAPARCHVAR       0x0001        /* Variable-sized architecture                         */
#define KGRAPHMAPRBMAPARCHCMPLT     0x0002        /* Complete-graph architecture                         */
#define KGRAPHMAPRBMAPPARTHALF      0x0004        /* Only update half of part array as mappings are tied */

/*
**  The type and structure definitions.
*/

/*+ Job pool structures. +*/

typedef struct KgraphMapRbMapPoolLink_ {
  struct KgraphMapRbMapPoolLink_ *  prev;         /*+ Pointer to previous link +*/
  struct KgraphMapRbMapPoolLink_ *  next;         /*+ Pointer to next link     +*/
} KgraphMapRbMapPoolLink;

/*+ This structure defines a job to be
    performed with respect to a partial
    mapping of a source graph.          +*/

typedef struct KgraphMapRbMapJob_ {
  KgraphMapRbMapPoolLink    poollink;             /*+ Link to job pool; TRICK: FIRST           +*/
  KgraphMapRbMapPoolLink *  poolptr;              /*+ Pointer to last/current job pool         +*/
  int                       poolflag;             /*+ Flag set if job in pool                  +*/
  Gnum                      prioval;              /*+ Job priority value by policy             +*/
  Gnum                      priolvl;              /*+ Priority level computed for this job     +*/
  ArchDom                   domnorg;              /*+ Domain to which the vertices belong      +*/
  Graph                     grafdat;              /*+ Job graph data (may be clone of another) +*/
  Anum                      vflonbr;              /*+ Number of fixed vertex load slots        +*/
  KgraphMapRbVflo *         vflotab;              /*+ Partial array of fixed vertex load slots +*/
} KgraphMapRbMapJob;

/*+ This structure defines the working data,
    for easier parameter passing.            +*/

typedef struct KgraphMapRbMapPoolData_ {
  int                       flagval;              /*+ Pool flag value                            +*/
  KgraphMapRbPolicy         polival;              /*+ Job selection policy                       +*/
  const Graph *             grafptr;              /*+ Pointer to top graph                       +*/
  const Anum *              pfixtax;              /*+ Pointer to fixed part array                +*/
  KgraphMapRbMapPoolLink    linktab[2];           /*+ Lists of jobs in pools                     +*/
  KgraphMapRbMapPoolLink *  pooltab[2];           /*+ Pointer to pools (same if tied)            +*/
  ArchDom *                 domntab[2];           /*+ Pointer to domain arrays (same if tied)    +*/
  KgraphMapRbMapJob *       jobtab;               /*+ Job table                                  +*/
  Mapping *                 mappptr;              /*+ Pointer to original mapping: current state +*/
} KgraphMapRbMapPoolData;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_RB_MAP
#define static
#endif

static int                  kgraphMapRbMapPoolInit (KgraphMapRbMapPoolData * restrict const, const KgraphMapRbData * restrict const);
static void                 kgraphMapRbMapPoolExit (KgraphMapRbMapPoolData * restrict const poolptr);
static void                 kgraphMapRbMapPoolAdd (KgraphMapRbMapPoolLink * restrict const, KgraphMapRbMapJob * const);
static KgraphMapRbMapJob *  kgraphMapRbMapPoolGet (KgraphMapRbMapPoolData * restrict const);
static void                 kgraphMapRbMapPoolFrst (KgraphMapRbMapPoolData * const, KgraphMapRbMapJob * const);
static void                 kgraphMapRbMapPoolUpdt1 (KgraphMapRbMapPoolData * const, const KgraphMapRbMapJob * const, const GraphPart * const, KgraphMapRbMapJob * const, const GraphPart);
static void                 kgraphMapRbMapPoolUpdt2 (KgraphMapRbMapPoolData * const, const KgraphMapRbMapJob * const, const GraphPart * const, KgraphMapRbMapJob * const, KgraphMapRbMapJob * const);

int                         kgraphMapRbMap      (const KgraphMapRbData * restrict const, const Graph * restrict const, const Anum, KgraphMapRbVflo * restrict const);

static int                  kgraphMapRbMapPoolResize (KgraphMapRbMapPoolData * restrict const);

#undef static

/*
**  The macro definitions.
*/

#define kgraphMapRbMapPoolEmpty(poolptr) ((poolptr)->pooltab[0]->next == &kgraphmaprbmappooldummy)
