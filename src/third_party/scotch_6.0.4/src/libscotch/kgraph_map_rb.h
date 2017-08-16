/* Copyright 2004,2007,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
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
/**                # Version 5.1  : from : 07 oct 2008     **/
/**                                 to     28 mar 2011     **/
/**                # Version 6.0  : from : 07 aug 2014     **/
/**                                 to     24 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Prime number for hashing terminal domain numbers. +*/

#define KGRAPHMAPRBVFLOHASHPRIME    17            /*+ Prime number for hashing +*/

/*+ Kinds of external edge processing. +*/

#define KGRAPHMAPRBVEEXNONE         0x0000        /* No options set */

#define KGRAPHMAPRBVEEXMAPP         0x0001        /* Graph mapping  */
#define KGRAPHMAPRBVEEXVFIX         0x0002        /* Fixed vertices */
#define KGRAPHMAPRBVEEXREMA         0x0004        /* Remapping      */

#define KGRAPHMAPRBVEEXEDGE         (KGRAPHMAPRBVEEXMAPP | KGRAPHMAPRBVEEXVFIX)
#define KGRAPHMAPRBVEEXVERT         (KGRAPHMAPRBVEEXREMA)

/*
**  The type and structure definitions.
*/

/*+ Job selection policy types. +*/

typedef enum KgraphMapRbPolicy_ {
  KGRAPHMAPRBPOLIRANDOM = 0,                      /*+ Random job selection policy                       +*/
  KGRAPHMAPRBPOLILEVEL,                           /*+ Select job with highest level                     +*/
  KGRAPHMAPRBPOLISIZE,                            /*+ Select job with largest size                      +*/
  KGRAPHMAPRBPOLINEIGHBOR,                        /*+ Priority level computed with respect to neighbors +*/
  KGRAPHMAPRBPOLINGLEVEL,                         /*+ Select job with most neighbors of higher level    +*/
  KGRAPHMAPRBPOLINGSIZE,                          /*+ Select job with most neighbors of smaller size    +*/
  KGRAPHMAPRBPOLIOLD                              /*+ Select job in old style (version 2.x)             +*/
} KgraphMapRbPolicy;

/*+ Method parameters. +*/

typedef struct KgraphMapRbParam_ {
  int                       flagjobtie;           /*+ Flag set of job pools are tied +*/
  int                       flagmaptie;           /*+ Flag set if mappings are tied  +*/
  KgraphMapRbPolicy         polival;              /*+ Job selection policy           +*/
  Strat *                   strat;                /*+ Bipartitioning strategy used   +*/
  double                    kbalval;              /*+ K-way imbalance ratio          +*/
} KgraphMapRbParam;

/*+ This structure holds the data passed to each bipartitioning job. +*/

typedef struct KgraphMapRbData_ {
  const Graph *             grafptr;              /*+ Pointer to top-level graph, possibly with fixed vertices +*/
  Mapping *                 mappptr;              /*+ Mapping to compute                                       +*/
  struct {                                        /*+ Remapping structure                                      +*/
    const Mapping *         mappptr;              /*+ Old mapping (for remapping only)                         +*/  
    const Gnum *            vmlotax;              /*+ Array of vertex migration costs                          +*/
    Gnum                    cmloval;              /*+ Migration edge load for remapping                        +*/
    Gnum                    crloval;              /*+ Regular edge load for mapping                            +*/
  } r;
  const Anum *              pfixtax;              /*+ Fixed vertex partition array                             +*/
  const KgraphMapRbParam *  paraptr;              /*+ Pointer to mapping parameter structure                   +*/
  double                    comploadrat;          /*+ Ideal load balance per weight unit                       +*/
  double                    comploadmin;          /*+ Minimum vertex load per target load                      +*/
  double                    comploadmax;          /*+ Maximum vertex load per target load                      +*/
} KgraphMapRbData;

/*+ Fixed vertex load type. An array of such
    cells stores the loads of strictly positive
    fixed vertices (zero ones are discarded) that
    must be assigned to some subdomain of the
    current domain to be bipartitioned.           +*/

typedef struct KgraphMapRbVflo_ {
  Anum                      termnum;              /*+ Terminal domain number +*/
  Gnum                      veloval;              /*+ Vertex load            +*/
} KgraphMapRbVflo;

/*+ Hash structure for merging fixed vertex
    domains with non-fixed vertex domains.  +*/

typedef struct KgraphMapRbVfloHash_ {
  Anum                      termnum;              /*+ Terminal domain number        +*/
  Anum                      domnnum;              /*+ Domain number in domain array +*/
} KgraphMapRbVfloHash;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_RB
#define static
#endif

int                         kgraphMapRb         (Kgraph * const, const KgraphMapRbParam * const);

int                         kgraphMapRbVfloBuild (const Arch * restrict const, const Graph * restrict const, const Gnum, const Anum * restrict const, Graph * restrict const, Anum * restrict const, KgraphMapRbVflo * restrict * restrict const);
void                        kgraphMapRbVfloSplit (const Arch * restrict const, const ArchDom * restrict const, const Anum, KgraphMapRbVflo * restrict const, Anum * restrict const, Gnum * restrict const);
int                         kgraphMapRbVfloMerge (Mapping * restrict const, const Gnum, const Anum * restrict const, const Anum);

int                         kgraphMapRbBgraph   (const KgraphMapRbData * restrict const, Bgraph * restrict const, const Graph * restrict const, const Mapping * restrict const, const ArchDom * restrict const, const Gnum * restrict const);

#undef static
