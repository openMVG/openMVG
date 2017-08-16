/* Copyright 2007,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_fold_comm.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed source graph        **/
/**                building routines.                      **/
/**                                                        **/
/**                # Version 5.0  : from : 23 may 2006     **/
/**                                 to     19 aug 2006     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to     03 jan 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Starting maximum number of communications per process +*/

#define DGRAPHFOLDCOMMNBR           4             /* Starting maximum number of communications per process */

/*
**  The type and structure definitions.
*/

/*+ Process communication type. Receivers that have
    more local vertices than needed can also be
    senders, hence the " = 1" to allow for or-ing
    both values.                                    +*/

typedef enum DgraphFoldCommType_ {
  DGRAPHFOLDCOMMRECV = 1,                         /*+ Process is a receiver +*/
  DGRAPHFOLDCOMMSEND                              /*+ Process is a sender   +*/
} DgraphFoldCommType;

/* Sort structure for processors. */

typedef struct DgraphFoldCommData_ {
  Gnum             vertnbr;                       /* Number of vertices; TRICK: FIRST   */
  Gnum             procnum;                       /* Processor index (Gnum for sorting) */
} DgraphFoldCommData;

/*
** The function prototypes.
*/

int                         dgraphFoldComm      (const Dgraph * restrict const, const int, int * restrict const, int * restrict const, DgraphFoldCommData * restrict * restrict const, Gnum * restrict * restrict const, Gnum * restrict const, int * restrict const, Gnum * restrict * restrict const, Gnum * restrict * restrict const);
