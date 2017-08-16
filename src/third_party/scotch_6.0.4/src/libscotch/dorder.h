/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dorder.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data           **/
/**                declarations for the generic            **/
/**                distributed ordering structure.         **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 15 apr 2006     **/
/**                                 to     14 oct 2007     **/
/**                # Version 5.1  : from : 28 nov 2007     **/
/**                                 to     04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Tag for MPI communications. +*/

#define DORDERTAGPERI               0

/*+ Column block separation tree cell flags.
    The DORDERCBLKNEDI value must correspond
    to a single bit and be equal to the
    ORDERCBLKNEDI value.                     +*/

#define DORDERCBLKNONE              0x0000        /*+ Not yet assigned                 +*/
#define DORDERCBLKNEDI              0x0001        /*+ Nested dissection separator node +*/
#define DORDERCBLKLEAF              0x0002        /*+ Distributed leaf                 +*/

/*
**  The type and structure definitions.
*/

/*+ Distributed index of column block
    tree cell.                        +*/

typedef struct DorderIndex_ {
  int                       proclocnum;           /*+ Number of this process      +*/
  Gnum                      cblklocnum;           /*+ Local index of column block +*/
} DorderIndex;

/*+ Link structure to double-chain all column blocks
    into the distributed ordering structure. Nodes
    are inserted at end of list such that a simple
    traversal gives nodes in ascending creation
    order, which is essential for locally-rooted
    nodes when gathering them on a centralized
    ordering.                                        +*/

typedef struct DorderLink_ {
  struct DorderLink_ *              nextptr;      /*+ Pointer to previous column block +*/
  struct DorderLink_ *              prevptr;      /*+ Pointer to next column block     +*/
} DorderLink;

/*+ Centralized column block node. +*/

typedef struct DorderNode_ {
  Gnum                              fathnum;      /*+ Number of father in centralized node array        +*/
  int                               typeval;      /*+ Centralized type of tree node                     +*/
  Gnum                              vnodnbr;      /*+ Number of nodes in this column block              +*/
  Gnum                              cblknum;      /*+ Rank of column block in father column block array +*/
} DorderNode;

/*+ Distributed column-block tree cell. Each
    cell defines a distributed column block,
    which is either a nested dissection node,
    with its two subgraphs and its separator,
    or a leaf. Leaves which are located on a
    single process can be nested dissection
    sequential nodes, with the sequential tree
    folded as a node array.
    Column blocks are given in ascending order
    within all sub-arrays, for proper infix
    traversal.                                 +*/

typedef struct DorderCblk_ {
  DorderLink                        linkdat;      /*+ Link to other blocks. TRICK: FIRST             +*/
  struct Dorder_ *                  ordelocptr;   /*+ Pointer to local distributed ordering          +*/
  int                               typeval;      /*+ Distributed type of tree node                  +*/
  DorderIndex                       fathnum;      /*+ Master index of parent column block            +*/
  DorderIndex                       cblknum;      /*+ Master index of this column block              +*/
  Gnum                              ordeglbval;   /*+ Un-based starting index of inverse permutation +*/
  Gnum                              vnodglbnbr;   /*+ Number of node vertices in subtree             +*/
  Gnum                              cblkfthnum;   /*+ Rank of node in father column block array      +*/
  union {
    struct {                                      /*+ Fragment of inverse permutation                +*/
      Gnum                          ordelocval;   /*+ Starting index of inverse permutation          +*/
      Gnum                          vnodlocnbr;   /*+ Number of node vertices in fragment            +*/
      Gnum *                        periloctab;   /*+ Pointer to inverse permutation fragment        +*/
      Gnum                          nodelocnbr;   /*+ Number of local column blocks                  +*/
      DorderNode *                  nodeloctab;   /*+ Array of local column blocks                   +*/
      Gnum                          cblklocnum;   /*+ Local number of first local column block       +*/
    } leaf;
    struct {                                      /*+ Fragment of inverse permutation                +*/
      Gnum                          cblkglbnbr;   /*+ Number of descendent nodes (2 or 3)            +*/
    } nedi;
  } data;
} DorderCblk;

/*+ Distributed ordering structure. A distributed
    block ordering is defined by fragments of its
    inverse permutation, distributed across all
    of the participating processes.
    For the sake of consistency between orderings
    that have been produced either from graphs or
    meshes, whether centralized or distributed, all
    ordering values are based from baseval.         +*/

typedef struct Dorder_ {
  Gnum                      baseval;              /*+ Base value for structures                                                      +*/
  Gnum                      vnodglbnbr;           /*+ Global number of node vertices                                                 +*/
  Gnum                      cblklocnbr;           /*+ Local number of unique locally-rooted distributed and sequential column blocks +*/
  DorderLink                linkdat;              /*+ Link to column blocks                                                          +*/
  MPI_Comm                  proccomm;             /*+ Ordering global communicator                                                   +*/
  int                       proclocnum;           /*+ Rank of this process in the communicator                                       +*/
#ifdef SCOTCH_PTHREAD
  pthread_mutex_t           mutelocdat;           /*+ Local mutex for counter and link updates                                       +*/
#endif /* SCOTCH_PTHREAD */
} Dorder;

/*
**  The function prototypes.
*/

#ifndef DORDER
#define static
#endif

int                         dorderInit          (Dorder * const, const Gnum, const Gnum, MPI_Comm);
void                        dorderExit          (Dorder * const);
void                        dorderFree          (Dorder * const);
#ifdef DGRAPH_H
int                         dorderPerm          (const Dorder * const, const Dgraph * const, Gnum * const);
int                         dorderSave          (const Dorder * const, const Dgraph * const, FILE * const);
int                         dorderSaveBlock     (const Dorder * const, const Dgraph * const, FILE * const);
int                         dorderSaveMap       (const Dorder * const, const Dgraph * const, FILE * const);
int                         dorderSaveTree      (const Dorder * const, const Dgraph * const, FILE * const);
#ifdef ORDER_H
int                         dorderSaveTree2     (const Dorder * restrict const, const Dgraph * restrict const, FILE * restrict const, int (*) (const Order * const, const Gnum * const, FILE * const));
#endif /* ORDER_H */
#endif /* DGRAPH_H */
Gnum                        dorderCblkDist      (const Dorder * restrict const);
int                         dorderTreeDist      (const Dorder * restrict const, const Dgraph * restrict const, Gnum * restrict const, Gnum * restrict const);
#ifdef ORDER_H
int                         dorderGather        (const Dorder * const, Order * const);
int                         dorderGatherTree    (const Dorder * const, Order * const, const int);
#endif /* ORDER_H */

DorderCblk *                dorderFrst          (Dorder * const);
DorderCblk *                dorderNew           (DorderCblk * const, MPI_Comm);
DorderCblk *                dorderNewSequ       (DorderCblk * const);
Gnum                        dorderNewSequIndex  (DorderCblk * const, const Gnum);
void                        dorderDispose       (DorderCblk * const);

#undef static
