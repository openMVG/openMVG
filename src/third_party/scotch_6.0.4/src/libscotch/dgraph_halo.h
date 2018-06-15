/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_halo.h                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the data declara-    **/
/**                tions for the asynchronous halo         **/
/**                exchange routine.                       **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 28 dec 2007     **/
/**                                 to   : 05 feb 2008     **/
/**                # Version 5.1  : from : 28 aug 2008     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
** The defines.
*/

/* procsidtab-related values. */

#define DGRAPHGHSTSIDMAX            ((int) ((unsigned int) (1 << (sizeof (int) * 8 - 1)) - 2U)) /* Maximum leap value for procsidtab entries */

/*
** The type and structure definitions.
*/

/* Sort structure for ghost edges. */

typedef struct DgraphHaloRequest_ {
  int                       flagval;
#ifdef SCOTCH_PTHREAD
  Dgraph *                  grafptr;              /* Pointer to graph data       */
  void *                    attrgsttab;           /* Attribute array to share    */
  MPI_Datatype              attrglbtype;          /* Attribute datatype          */
  pthread_t                 thrdval;              /* Data of asynchronous thread */
#else /* SCOTCH_PTHREAD */
#ifdef SCOTCH_MPI_ASYNC_COLL
  byte *                    attrsndtab;           /* Group leader for memory freeing        */
  MPI_Request               requval;              /* MPI asynchronous communication request */
#endif /* SCOTCH_MPI_ASYNC_COLL */
#endif /* SCOTCH_PTHREAD */
} DgraphHaloRequest;

/*
** The function prototypes.
*/

#ifndef DGRAPH_HALO
#define static
#endif

#ifdef SCOTCH_PTHREAD
static void *               dgraphHaloAsync2    (DgraphHaloRequest * restrict);
#endif /* SCOTCH_PTHREAD */
void                        dgraphHaloAsync     (Dgraph * restrict const, void * restrict const, const MPI_Datatype, DgraphHaloRequest * restrict);
int                         dgraphHaloWait      (DgraphHaloRequest * restrict);

int                         dgraphHaloCheck     (const Dgraph * restrict const);

#undef static
