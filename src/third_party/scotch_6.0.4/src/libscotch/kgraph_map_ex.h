/* Copyright 2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_map_ex.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the load balancing graph mapping    **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 08 jun 2011     **/
/**                                 to     08 jun 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct KgraphMapExParam_ {
  double                    kbalval;              /*+ Imbalance ratio +*/
} KgraphMapExParam;

/*+ This structure holds extended domain
    information. sonstab[0] == -1 if node
    is terminal. Else, sonstab[1] == -1
    means only one branch is considered
    at this level.                        +*/

typedef struct KgraphMapExDom_ {
  Anum                      treenum;              /*+ Tree node index for this domain +*/
  Anum                      domnwght;             /*+ Domain weight                   +*/
  Gnum                      compload;             /*+ Current load in domain          +*/
  Gnum                      comploadmax;          /*+ Maximum load allowed in domain  +*/
} KgraphMapExDom;

/*+ This structure records the best
    candidate domain found to date. +*/

typedef struct KgraphMapExFind_ {
  Gnum                      comploaddlt;          /*+ Best imbalance +*/
  Anum                      domnnum;              /*+ Domain number  +*/
} KgraphMapExFind;

/*+ This structure allows one to
    sort vertices by vertex weight. +*/

typedef struct KgraphMapExSort_ {
  Gnum                      veloval;              /*+ Vertex load   +*/
  Gnum                      vertnum;              /*+ Vertex number +*/
} KgraphMapExSort;

/*+ This structure allows one to sort
    existing terminal domains by terminal
    number.                               +*/

typedef struct KgraphMapExTerm_ {
  Anum                      termnum;              /*+ Domain terminal number  +*/
  Anum                      domnnum;              /*+ Domain index in domntab +*/
} KgraphMapExTerm;

/*+ This structure holds a recursive bi-mapping tree node. +*/

typedef struct KgraphMapExTree_ {
  Anum                      fathnum;              /*+ Index of father node; -1 if root      +*/
  Anum                      sonstab[2];           /*+ Index of sons; [0] == -1 for terminal +*/
  ArchDom                   domndat;              /*+ Subdomain data                        +*/
} KgraphMapExTree;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_EX
#define static
#endif

int                         kgraphMapEx         (Kgraph * restrict const, const KgraphMapExParam * const);

static Anum                 kgraphMapExTree     (const Arch * restrict const, const KgraphMapExTerm * restrict const, const Anum, KgraphMapExDom * restrict const, KgraphMapExTree * restrict const, Anum * restrict const, const ArchDom * restrict const);
static Anum                 kgraphMapExFind     (const Arch * restrict const, const KgraphMapExTree * restrict const, const KgraphMapExDom * restrict const, const Anum, const Gnum);
static int                  kgraphMapExFind2    (const Arch * restrict const, const KgraphMapExTree * restrict const, const KgraphMapExDom * restrict const, KgraphMapExFind * restrict const, const Anum, const Anum, const Gnum);

#undef static
