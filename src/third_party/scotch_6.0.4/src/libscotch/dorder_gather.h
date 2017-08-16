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
/**   NAME       : dorder_gather.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the distributed ordering gathering  **/
/**                function.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 20 jul 2007     **/
/**                                 to   : 27 jul 2007     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Tag for MPI communications. +*/

#define DORDERGATHERNODESIZE        (sizeof (DorderGatherNode) / sizeof (Gnum))

/*
**  The type and structure definitions.
*/

/*+ This structure holds the tree definition. +*/

typedef struct DorderGatherLeaf_ {
  Gnum                      ordelocval;           /*+ Starting index of inverse permutation +*/
  Gnum                      vnodlocnbr;           /*+ Number of node vertices in fragment   +*/
} DorderGatherLeaf;

/*+ This structure holds the separator tree structure.
    Because arrays of this structure is to be sent as
    a single contiguous array, all its fields must be
    of the same type.                                  +*/

typedef struct DorderGatherNode_ {
  Gnum                      fathnum;              /*+ Global number of father node              +*/
  Gnum                      typeval;              /*+ Node type                                 +*/
  Gnum                      vnodnbr;              /*+ Number of node vertices                   +*/
  Gnum                      cblknum;              /*+ Rank of node in father column block array +*/
} DorderGatherNode;

/*+ This structure holds the separator tree structure. +*/

typedef struct DorderGatherCblk_ {
  Gnum                      cblknbr;              /*+ Number of sons       +*/
  OrderCblk *               cblktab;              /*+ Pointer to sub-array +*/
} DorderGatherCblk;
