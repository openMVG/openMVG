/* Copyright 2007-2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_fm.h                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the improved Fiduccia-Mattheyses    **/
/**                refinement routine for the vertex over- **/
/**                lapped graph partitioning.              **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 30 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Gain table subbits. +*/

#define WGRAPHSEPAFMGAINBITS        4

/*+ Prime number for hashing vertex numbers. +*/

#define WGRAPHSEPAFMHASHPRIME       17            /*+ Prime number for hashing +*/

/*+ +*/

#define WGRAPHSEPAFMMAXPARTLOAD     ((Gnum) (((Gunum) ~0) >> 1))

/*+ Gain table vertex status. +*/

#define WGRAPHSEPAFMSTATEFREE       ((GainLink *) 0) /*+ Vertex is free or separator-chained  +*/
#define WGRAPHSEPAFMSTATESUCH       ((GainLink *) 1) /*+ Separator vertex is used and chained +*/
#define WGRAPHSEPAFMSTATEUSED       ((GainLink *) 2) /*+ Vertex already swapped once          +*/
#define WGRAPHSEPAFMSTATELINK       ((GainLink *) 3) /*+ Currently in gain table if higher    +*/

/*+ Save type identifier +*/

#define WGRAPHSEPAFMSAVEMOVE        0
#define WGRAPHSEPAFMSAVELINKDEL     1
#define WGRAPHSEPAFMSAVELINKADD     2
#define WGRAPHSEPAFMSAVELOAD        3

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct WgraphPartFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltrat;              /*+ Maximum weight imbalance ratio                       +*/
} WgraphPartFmParam;


typedef struct WgraphPartFmPartList_ {
  struct WgraphPartFmPartList_ *   prev;
  Gnum                             gain;
  struct WgraphPartFmPartList_ *   loadprev;
  Gnum                             loadgain;
  Gnum                             sizegain;
  Gnum                             isinloadlist;
} WgraphPartFmPartList;

typedef struct WgraphPartFmVertex_ {
  Gnum                          partval;          /*+ Vertex part TRICK: same type as vertload        +*/
  Gnum                          partval2;         /*+ Vertex part into which the current vertex moves +*/
  Gnum                          vertnum;          /*+ Number of vertex in hash table                  +*/
  struct WgraphPartFmVertex_ *  prev;
  struct WgraphPartFmLink_ *    linklist;
  Gnum                          linked;
  struct WgraphPartFmVertex_ *  lockprev;
} WgraphPartFmVertex;

/*+ The move recording structure. +*/

typedef struct WgraphPartFmLink_ {
  GainLink                    gainlink;           /*+ Gain link: FIRST +*/
  Gnum                        partval;
  Gnum                        hashnum;
  Gnum                        gain;
  struct WgraphPartFmLink_ *  next;
  Gnum                        minloadpartval;
  Gnum                        minloadpartload;  
} WgraphPartFmLink;

typedef struct WgraphPartFmSave_ {
  Gnum type;
  union {
      struct {
          Gnum              hashnum;              /*+ Number of hash slot for saved vertex +*/
          Gnum              partval;              /*+ Saved vertex part value              +*/
      } movedata;
      struct {
          Gnum              linknum;
          Gnum              gain;
      } linkdata;
      struct {
          Gnum              loaddiff;
          Gnum              sizediff;
          Gnum              partval;
      } loaddata;
  } u;
} WgraphPartFmSave;

/*
**  The function prototypes.
*/

#ifndef WGRAPH_PART_FM
#define static
#endif

int                         wgraphPartFm        (Wgraph * const, const WgraphPartFmParam * const);

static int                  wgraphPartFmResize  ();
static int                  wgraphPartFmCheck   (const Wgraph * restrict const, const WgraphPartFmVertex * restrict const, Gnum, const Gnum, const Gnum);

#undef static
