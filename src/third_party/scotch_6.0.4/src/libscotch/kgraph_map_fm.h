/* Copyright 2004,2010-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_fm.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the k-way Fiduccia-Mattheyses       **/
/**                mapping algorithm.                      **/
/**                                                        **/
/**   DATES      : # Version 3.3  : from : 10 may 1999     **/
/**                                 to     20 jun 1999     **/
/**                # Version 3.4  : from : 27 sep 1999     **/
/**                                 to     13 nov 1999     **/
/**                # Version 5.0  : from : 11 oct 2006     **/
/**                                 to     12 oct 2006     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     19 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Prime number for hashing vertex numbers. +*/

#define KGRAPHMAPFMHASHPRIME        17            /*+ Prime number for hashing +*/

/*+ Ratio of unused slots before compacting extended edge array. +*/

#define KGRAPHMAPFMEDXXCOMP         5             /*+ Compact if edxunbr > (edxxnbr / KGRAPHMAPFMEDXXCOMP) */ 

/*+ Save type identifier +*/

#define KGRAPHMAPPFMSAVEVEXX        0
#define KGRAPHMAPPFMSAVEEDXX        1
#define KGRAPHMAPPFMSAVELINK        2             /* Bit value for KGRAPHMAPPFMSAVELINKADD and KGRAPHMAPPFMSAVELINKDEL */
#define KGRAPHMAPPFMSAVELINKDEL     2
#define KGRAPHMAPPFMSAVELINKADD     3

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct KgraphMapFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltval;              /*+ Maximum weight imbalance ratio                       +*/
} KgraphMapFmParam;

/*+ The extended edge structure. In fact, this
    structure does not represent edges, but
    possible moves in another domain.
    For trick reasons, the gain table data
    structure must be the first field of the
    structure.                                 +*/

#ifdef SCOTCH_TABLE_GAIN

typedef GainTabl * KgraphMapFmTabl;
typedef GainLink KgraphMapFmLink;

#else /* SCOTCH_TABLE_GAIN */

typedef FiboTree KgraphMapFmTabl;
typedef FiboNode KgraphMapFmLink;

#endif /* SCOTCH_TABLE_GAIN */

typedef struct KgraphMapFmEdge_ {
  KgraphMapFmLink           gainlink;             /*+ Gain link; TRICK: FIRST                 +*/
  Gnum                      commgain;             /*+ Communication gain                      +*/
  Gnum                      cmiggain;             /*+ Migration communication gain            +*/
  Gnum                      cmigmask;             /*+ Migration communication mask            +*/
  Gnum                      edlosum;              /*+ Sum of edge loads linking to the domain +*/
  Gnum                      edgenbr;              /*+ Number of edges linking to the domain   +*/
  Anum                      domnnum;              /*+ Destination domain index                +*/
  Anum                      distval;              /*+ Distance between the two domains        +*/
  Gnum                      vexxidx;              /*+ Index of owner vertex in vertex array   +*/
  Gnum                      edxxidx;              /*+ Index of next edge in edge array        +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded +*/
} KgraphMapFmEdge;

/*+ The hash vertex structure. +*/

typedef struct KgraphMapFmVertex_ {
  struct KgraphMapFmVertex_ * lockptr;            /*+ Pointer to next vertex in lock list or NULL; FIRST +*/
  Gnum                        vertnum;            /*+ Number of vertex                                   +*/
  Gnum                        cmigload;           /*+ Migration communication load                       +*/
  Gnum                        edlosum;            /*+ Sum of edge loads linking to self domain           +*/
  Gnum                        edgenbr;            /*+ Number of edges linking to self domain             +*/
  Anum                        domnnum;            /*+ Domain number                                      +*/
  ArchDom *                   domoptr;            /*+ Domain in old mapping (for repartitioning)         +*/
  Gnum                        veloval;            /*+ Vertex load; negative when locked                  +*/
  Gnum                        edxxidx;            /*+ Index of first element in edge array               +*/
  Gnum                        mswpnum;            /*+ Number of move sweep when data recorded            +*/
} KgraphMapFmVertex;

/*+ The move recording structures. +*/

typedef struct KgraphMapFmSave_ {
  Gnum                      type;
  union {
    struct {
      Gnum                  vexxidx;              /*+ Index of vertex slot in hash table (vexxhab) +*/
      Gnum                  veloval;              /*+ Vertex load                                  +*/
      Anum                  domnnum;              /*+ Original vertex domain                       +*/
      Gnum                  commload;             /*+ Communication load for current domain        +*/
      Gnum                  cmigload;             /*+ Migration communication load                 +*/
      Gnum                  edlosum;              /*+ Sum of edge loads linking to self domain     +*/
      Gnum                  edgenbr;              /*+ Number of edges linking to self domain       +*/
    } vexxdat;
    struct {
      Gnum                  edxxidx;              /*+ Index of edge in edge array                  +*/
      Anum                  domnnum;              /*+ Destination domain index                     +*/
      Anum                  distval;              /*+ Distance between the two domains             +*/
      Gnum                  commgain;             /*+ Communication gain                           +*/
      Gnum                  cmiggain;             /*+ Migration communication gain                 +*/
      Gnum                  edlosum;              /*+ Sum of edge loads linking to the domain      +*/
      Gnum                  edgenbr;              /*+ Number of edges linking to the domain        +*/
    } edxxdat;
    struct {
      Gnum                  edxxidx;              /*+ Index of extended vertex or edge             +*/
      Gnum                  vexxidx;              /*+ Index of vertex slot in hash table (vexxhab) +*/
    } linkdat;
  } u;
} KgraphMapFmSave;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_FM
#define static
#endif

int                         kgraphMapFm         (Kgraph * restrict const, const KgraphMapFmParam * const);

#undef static

/*
**  The macro definitions.
*/

#ifdef SCOTCH_TABLE_GAIN

/*+ Gain table subbits. +*/

#define KGRAPHMAPFMSUBBITS       4

/*+ Service routines. +*/

#define kgraphMapFmTablInit(t)      (((*(t)) = gainTablInit (GAINMAX, KGRAPHMAPFMSUBBITS)) == NULL)
#define kgraphMapFmTablFree(t)      gainTablFree (*(t))
#define kgraphMapFmTablExit(t)      do {                     \
                                      if (*(t) != NULL)      \
                                        gainTablExit (*(t)); \
                                    } while (0)
#define kgraphMapFmTablAdd(t,e)     gainTablAdd ((*(t)), &(e)->gainlink, ((e)->commgain + (((e)->cmiggain) & ((e)->cmigmask))) * (e)->distval)
#define kgraphMapFmTablDel(t,e)     gainTablDel ((*(t)), &(e)->gainlink)
#else /* SCOTCH_TABLE_GAIN */

/*+ Service routines. +*/

#define kgraphMapFmTablInit(t)      (fiboTreeInit ((t), kgraphMapFmCmpFunc))
#define kgraphMapFmTablFree(t)      fiboTreeFree (t)
#define kgraphMapFmTablExit(t)      fiboTreeExit (t)
#define kgraphMapFmTablAdd(t,e)     fiboTreeAdd ((t), &(e)->gainlink)
#define kgraphMapFmTablDel(t,e)     fiboTreeDel ((t), &(e)->gainlink)

#endif /* SCOTCH_TABLE_GAIN */

#define kgraphMapFmLock(l,v)        do {                                        \
                                      (v)->lockptr = (KgraphMapFmVertex *) (l); \
                                      (l) = (v);                                \
                                    } while (0)
#define kgraphMapFmLockNext(v)      ((KgraphMapFmVertex *) (v)->lockptr)
