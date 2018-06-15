/* Copyright 2004,2007,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_fm.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for our Improved Fiduccia-Mattheyses    **/
/**                bipartitioning algorithm.               **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 30 sep 1993     **/
/**                                 to     09 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     13 apr 1994     **/
/**                # Version 2.0  : from : 04 jul 1994     **/
/**                                 to     25 nov 1994     **/
/**                # Version 3.0  : from : 06 jul 1995     **/
/**                                 to     06 jul 1995     **/
/**                # Version 3.1  : from : 06 nov 1995     **/
/**                                 to     07 jun 1996     **/
/**                # Version 3.2  : from : 21 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     12 mar 1999     **/
/**                # Version 4.0  : from : 27 aug 2004     **/
/**                                 to     27 aug 2004     **/
/**                # Version 6.0  : from : 23 fev 2011     **/
/**                                 to     05 sep 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Prime number for hashing vertex numbers. +*/

#define BGRAPHBIPARTFMHASHPRIME     17            /*+ Prime number for hashing +*/

/*
**  The type and structure definitions.
*/

/*+ Move type. +*/

typedef enum BgraphBipartFmType_ {
  BGRAPHBIPARTFMTYPEALL,                          /*+ All vertices           +*/
  BGRAPHBIPARTFMTYPEBOUNDARY                      /*+ Boundary vertices only +*/
} BgraphBipartFmType;

/*+ This structure holds the method parameters. +*/

typedef struct BgraphBipartFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltval;              /*+ Maximum weight imbalance ratio                       +*/
  BgraphBipartFmType        typeval;              /*+ Whether considered vertices are boundary or all      +*/
} BgraphBipartFmParam;

#ifdef BGRAPH_BIPART_FM                           /* Private part of the module */

/*+ The hash vertex structure. For trick
    reasons, the gain table data structure
    must be the first field of the structure. +*/

#ifdef SCOTCH_TABLE_GAIN

typedef GainTabl * BgraphBipartFmTabl;
typedef GainLink BgraphBipartFmLink;

#else /* SCOTCH_TABLE_GAIN */

typedef FiboTree BgraphBipartFmTabl;
typedef FiboNode BgraphBipartFmLink;

#endif /* SCOTCH_TABLE_GAIN */

typedef struct BgraphBipartFmVertex_ {
  BgraphBipartFmLink        gainlink;             /*+ Gain link: FIRST                        +*/
  Gnum                      vertnum;              /*+ Number of vertex                        +*/
  int                       partval;              /*+ Vertex part                             +*/
  Gnum                      compgain;             /*+ Computation gain                        +*/
  Gnum                      commgain;             /*+ Communication gain                      +*/
  Gnum                      commcut;              /*+ Cut edges                               +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded +*/
} BgraphBipartFmVertex;

/*+ The move recording structure. +*/

typedef struct BgraphBipartFmSave_ {
  Gnum                      hashnum;              /*+ Number of vertex slot +*/
  int                       partval;              /*+ Vertex part           +*/
  Gnum                      compgain;             /*+ Computation gain      +*/
  Gnum                      commgain;             /*+ Communication gain    +*/
  Gnum                      commcut;              /*+ Cut edges             +*/
} BgraphBipartFmSave;

/*
**  The function prototypes.
*/

static BgraphBipartFmVertex * bgraphBipartFmTablGet (BgraphBipartFmTabl * restrict const, const Gnum, const Gnum, const Gnum);

static int                  bgraphBipartFmResize (BgraphBipartFmVertex * restrict *, Gnum * restrict const, Gnum * const, BgraphBipartFmSave * restrict *, const Gnum, BgraphBipartFmTabl * const, BgraphBipartFmVertex ** const);
#ifdef SCOTCH_DEBUG_BGRAPH3
static int                  bgraphBipartFmCheck (const Bgraph * restrict const, const BgraphBipartFmVertex * restrict const, const Gnum, const int, const Gnum, const Gnum, const Gnum);
#endif /* SCOTCH_DEBUG_BGRAPH3 */

#endif /* BGRAPH_BIPART_FM */

int                         bgraphBipartFm      (Bgraph * restrict const, const BgraphBipartFmParam * const);

/*
**  The macro definitions.
*/

#ifdef SCOTCH_TABLE_GAIN

/*+ Gain table subbits. +*/

#define BGRAPHBIPARTFMSUBBITS       4

/** Gain table vertex status. **/

#define BGRAPHBIPARTFMSTATEFREE     ((GainLink *) 0) /*+ Vertex in initial state           +*/
#define BGRAPHBIPARTFMSTATEUSED     ((GainLink *) 1) /*+ Swapped vertex                    +*/
#define BGRAPHBIPARTFMSTATELINK     ((GainLink *) 2) /*+ Currently in gain table if higher +*/

/*+ Service routines. +*/

#define bgraphBipartFmTablInit(t)   (((*(t)) = gainTablInit (GAINMAX, BGRAPHBIPARTFMSUBBITS)) == NULL)
#define bgraphBipartFmTablFree(t)   gainTablFree (*(t))
#define bgraphBipartFmTablExit(t)   do {                     \
                                      if (*(t) != NULL)      \
                                        gainTablExit (*(t)); \
                                    } while (0)
#define bgraphBipartFmTablAdd(t,v)  gainTablAdd ((*(t)), &(v)->gainlink, (v)->commgain)
#define bgraphBipartFmTablDel(t,v)  gainTablDel ((*(t)), &(v)->gainlink)
#define bgraphBipartFmIsFree(v)     ((v)->gainlink.next == BGRAPHBIPARTFMSTATEFREE)
#define bgraphBipartFmIsTabl(v)     ((v)->gainlink.next >= BGRAPHBIPARTFMSTATELINK)
#define bgraphBipartFmIsUsed(v)     ((v)->gainlink.next == BGRAPHBIPARTFMSTATEUSED)
#define bgraphBipartFmSetFree(v)    do {                                            \
                                      (v)->gainlink.next = BGRAPHBIPARTFMSTATEFREE; \
                                    } while (0)
#define bgraphBipartFmSetUsed(v)    do {                                            \
                                      (v)->gainlink.next = BGRAPHBIPARTFMSTATEUSED; \
                                    } while (0)
#define bgraphBipartFmChain(l,v)    do {                                      \
                                      (v)->gainlink.prev = (GainLink *) *(l); \
                                      *(l) = (v);                             \
                                    } while (0)
#define bgraphBipartFmChainNext(v)  ((BgraphBipartFmVertex *) (v)->gainlink.prev)

#else /* SCOTCH_TABLE_GAIN */

/** Gain table vertex status. **/

#define BGRAPHBIPARTFMSTATEFREE     ((FiboNode *) 0) /*+ Vertex in initial state           +*/
#define BGRAPHBIPARTFMSTATEUSED     ((FiboNode *) 1) /*+ Swapped vertex                    +*/
#define BGRAPHBIPARTFMSTATELINK     ((FiboNode *) 2) /*+ Currently in gain table if higher +*/

/*+ Service routines. +*/

#define bgraphBipartFmTablInit(t)   (fiboTreeInit ((t), bgraphBipartFmCmpFunc))
#define bgraphBipartFmTablFree(t)   fiboTreeFree (t)
#define bgraphBipartFmTablExit(t)   fiboTreeExit (t)
#define bgraphBipartFmTablAdd(t,v)  fiboTreeAdd ((t), &(v)->gainlink)
#define bgraphBipartFmTablDel(t,v)  fiboTreeDel ((t), &(v)->gainlink)
#define bgraphBipartFmIsFree(v)     ((v)->gainlink.linkdat.nextptr == BGRAPHBIPARTFMSTATEFREE)
#define bgraphBipartFmIsTabl(v)     ((v)->gainlink.linkdat.nextptr >= BGRAPHBIPARTFMSTATELINK)
#define bgraphBipartFmIsUsed(v)     ((v)->gainlink.linkdat.nextptr == BGRAPHBIPARTFMSTATEUSED)
#define bgraphBipartFmSetFree(v)    do {                                                       \
                                      (v)->gainlink.linkdat.nextptr = BGRAPHBIPARTFMSTATEFREE; \
                                    } while (0)
#define bgraphBipartFmSetUsed(v)    do {                                                       \
                                      (v)->gainlink.linkdat.nextptr = BGRAPHBIPARTFMSTATEUSED; \
                                    } while (0)
#define bgraphBipartFmChain(l,v)    do {                                                 \
                                      (v)->gainlink.linkdat.prevptr = (FiboNode *) *(l); \
                                      *(l) = (v);                                        \
                                    } while (0)
#define bgraphBipartFmChainNext(v)  ((BgraphBipartFmVertex *) (v)->gainlink.linkdat.prevptr)

#endif /* SCOTCH_TABLE_GAIN */
