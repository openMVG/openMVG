/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : order.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data           **/
/**                declarations for the generic ordering   **/
/**                structure.                              **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 oct 1996     **/
/**                                 to   : 21 aug 1998     **/
/**                # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     28 dec 2004     **/
/**                # Version 5.0  : from : 25 jul 2007     **/
/**                                 to     25 jul 2007     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to     04 nov 2010     **/
/**                                                        **/
/************************************************************/

#define ORDER_H

/*
**  The defines.
*/

/*+ Ordering option flags. +*/

#define ORDERNONE                   0x0000        /* No options set                 */
#define ORDERFREEPERI               0x0001        /* Free inverse permutation array */

/*+ Column block separation tree cell flags.
    The ORDERCBLKNEDI value must correspond
    to a single bit and be equal to the
    DORDERCBLKNEDI value.                    +*/

#define ORDERCBLKOTHR               0x0000        /*+ Other ordering node              +*/
#define ORDERCBLKNEDI               0x0001        /*+ Nested dissection separator node +*/

/*
**  The type and structure definitions.
*/

/*+ Column-block tree node. Each node
    defines a column block, which is either
    a separator or a leaf built by nested
    dissection, or a super-variable built
    by minimum-degree algorithms. Column
    blocks are given in ascending order
    within sub-arrays, for proper infix
    traversal.                              +*/

typedef struct OrderCblk_ {
  int                       typeval;              /*+ Type of tree node                  +*/
  Gnum                      vnodnbr;              /*+ Number of node vertices in subtree +*/
  Gnum                      cblknbr;              /*+ Number of descendent column blocks +*/
  struct OrderCblk_ *       cblktab;              /*+ Sub-array of column-blocks         +*/
} OrderCblk;

/*+ Ordering structure. A block ordering is
    defined by its inverse permutation peritab
    and by the tree of permuted ordering blocks,
    which, once flattened, defines the blocks
    of the ordering. For the sake of consistency
    between orderings that have been produced
    either from graphs or meshes, all ordering
    values are based from baseval.               +*/

typedef struct Order_ {
  int                       flagval;              /*+ Flag value                          +*/
  Gnum                      baseval;              /*+ Base value for structures           +*/
  Gnum                      vnodnbr;              /*+ Number of node vertices             +*/
  Gnum                      treenbr;              /*+ Number of column block tree nodes   +*/
  Gnum                      cblknbr;              /*+ Number of column blocks             +*/
  OrderCblk                 cblktre;              /*+ Root of column block tree           +*/
  Gnum *                    peritab;              /*+ Inverse permutation array [vnodnbr] +*/
} Order;

/*
**  The function prototypes.
*/

#ifndef ORDER
#define static
#endif

int                         orderInit           (Order * const, const Gnum, const Gnum, Gnum * const);
void                        orderExit           (Order * const);
static void                 orderExit2          (OrderCblk * const, const Gnum);
int                         orderLoad           (Order * restrict const, const Gnum * restrict const, FILE * restrict const);
int                         orderSave           (const Order * restrict const, const Gnum * restrict const, FILE * restrict const);
int                         orderSaveMap        (const Order * const, const Gnum * restrict const, FILE * restrict const);
int                         orderSaveTree       (const Order * const, const Gnum * restrict const, FILE * restrict const);
void                        orderPeri           (const Gnum * const, const Gnum, const Gnum, Gnum * const, const Gnum);
void                        orderRang           (const Order * const, Gnum * const);
static void                 orderRang2          (Gnum ** const, Gnum * const, const OrderCblk * const);
void                        orderTree           (const Order * restrict const, Gnum * restrict const);
static void                 orderTree2          (Gnum * restrict const, Gnum * restrict const, const OrderCblk * restrict const, Gnum);
int                         orderCheck          (const Order * const);

#undef static
