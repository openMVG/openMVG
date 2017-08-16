/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vgraph_separate_fm.h                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the improved Fiduccia-Mattheyses    **/
/**                graph separation routine.               **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 02 nov 1997     **/
/**                                 to     20 nov 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     28 dec 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to     18 aug 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Gain table subbits. +*/

#define VGRAPHSEPAFMGAINBITS        4

/*+ Prime number for hashing vertex numbers. +*/

#define VGRAPHSEPAFMHASHPRIME       17            /*+ Prime number for hashing +*/

/*+ Gain table vertex status. +*/

#define VGRAPHSEPAFMSTATEFREE       ((GainLink *) 0) /*+ Vertex is free or separator-chained  +*/
#define VGRAPHSEPAFMSTATESUCH       ((GainLink *) 1) /*+ Separator vertex is used and chained +*/
#define VGRAPHSEPAFMSTATEUSED       ((GainLink *) 2) /*+ Vertex already swapped once          +*/
#define VGRAPHSEPAFMSTATELINK       ((GainLink *) 3) /*+ Currently in gain table if higher    +*/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct VgraphSeparateFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltrat;              /*+ Maximum weight imbalance ratio                       +*/
} VgraphSeparateFmParam;

/*+ The hash vertex structure. For trick reasons,
    one of the gain table data structures is followed
    by a negative integer, and the other by a positive
    one. Thus, one can deduce the value of the pointer
    to the structure from a pointer to any of the gain
    table data structures.
    Moreover, some fields have special meaning:
    - gainlink0.next: state of vertex (see
      VGRAPHSEPAFMSTATEXXXX).
    - gainlink0.prev: simple chaining for separator
      vertices, if vertex is in chained state
      (((vertpart == 2) &&
       (gainlink0.next == VGRAPHSEPAFMSTATEFREE)) ||
      (gainlink0.next == VGRAPHSEPAFMSTATESUCH)).
    - gainlink1: double chained list of locked vertices,
      if ((gainlink0.next == VGRAPHSEPAFMSTATESUCH) ||
      (gainlink0.next == VGRAPHSEPAFMSTATEUSED)).        +*/

typedef struct VgraphSeparateFmVertex_ {
  GainLink                  gainlink0;            /*+ Gain link if moved to part 0; FIRST                     +*/
  Gnum                      veloval;              /*+ TRICK: opposite of vertex load                          +*/
  GainLink                  gainlink1;            /*+ Gain link if moved to part 1; TRICK: before vertpart    +*/
  Gnum                      partval;              /*+ Vertex part TRICK: same type as vertload                +*/
  Gnum                      compgain[2];          /*+ Separator gain if moved to given part; TRICK: not first +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded                 +*/
  Gnum                      vertnum;              /*+ Number of vertex in hash table                          +*/
} VgraphSeparateFmVertex;

/*+ The move recording structure. +*/

typedef struct VgraphSeparateFmSave_ {
  Gnum                      hashnum;              /*+ Number of hash slot for saved vertex +*/
  int                       partval;              /*+ Saved vertex part value              +*/
  Gnum                      compgain[2];          /*+ Saved vertex gain                    +*/
} VgraphSeparateFmSave;

/*
**  The function prototypes.
*/

#ifndef VGRAPH_SEPARATE_FM
#define static
#endif

int                         vgraphSeparateFm    (Vgraph * const, const VgraphSeparateFmParam * const);

static int                  vgraphSeparateFmResize (VgraphSeparateFmVertex * restrict * hashtabptr, Gnum * const, Gnum * const, VgraphSeparateFmSave * restrict *, const Gnum, GainTabl * const, GainLink * const);
#ifdef SCOTCH_DEBUG_VGRAPH3
static int                  vgraphSeparateFmCheck (const Vgraph * const, const VgraphSeparateFmVertex * restrict const, const Gnum, const Gnum, const Gnum);
#endif /* SCOTCH_DEBUG_VGRAPH3 */
static GainLink *           vgraphSeparateFmTablGet (GainTabl * const, const Gnum, const Gnum, const int);

#undef static
