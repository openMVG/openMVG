/* Copyright 2007,2009,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_halo.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source graph handling routines of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 17 jul 2007     **/
/**                                 to     02 aug 2007     **/
/**                # Version 5.1  : from : 02 jul 2008     **/
/**                                 to     17 nov 2010     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph.h"
#include "dgraph_halo.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/*+ This routine requests the computation
*** of the ghost edge array.
*** It returns:
*** - 0   : if the computation succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphGhst (
SCOTCH_Dgraph * const       grafptr)
{
  return (dgraphGhst ((Dgraph *) grafptr));
}

/*+ This routine requests the computation of the
*** ghost edge array in replacement of the global
*** edge array.
*** It returns:
*** - 0   : if the computation succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphGhstReplace (
SCOTCH_Dgraph * const       grafptr)
{
  Dgraph * restrict   srcgrafptr;                 /* Pointer to scotch graph */
  DgraphFlag          srcflagval;                 /* Graph properties        */
  int                 o;

  srcgrafptr = (Dgraph *) grafptr;
  srcflagval = srcgrafptr->flagval;
  srcgrafptr->flagval |= DGRAPHFREETABS;          /* If edge array was not allocated internally, assume it was */

  o = dgraphGhstReplace (srcgrafptr);

  srcgrafptr->flagval = (srcgrafptr->flagval & ~DGRAPHFREETABS) | srcflagval; /* Restore original allocation flag */

  return (o);
}

/*+ This routine spreads local information
*** borne by local vertices across the ghost
*** vertices of the neighboring processes.
*** It returns:
*** - 0   : if the exchange succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphHalo (
SCOTCH_Dgraph * const       grafptr,
void * const                datatab,
const MPI_Datatype          typeval)
{
  return (dgraphHaloSync ((Dgraph *) grafptr, (byte *) datatab, typeval));
}

/*+ This routine spreads local information
*** borne by local vertices across the ghost
*** vertices of the neighboring processes, in
*** an asynchronous way.
*** It returns:
*** - 0   : if the exchange succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphHaloAsync (
SCOTCH_Dgraph * const         grafptr,
void * const                  datatab,
const MPI_Datatype            typeval,
SCOTCH_DgraphHaloReq * const  requptr)
{
  dgraphHaloAsync ((Dgraph *) grafptr, (byte *) datatab, typeval, (DgraphHaloRequest *) requptr);
  return (0);
}

/*+ This routine waits for the termination of
*** an asynchronous halo request.
*** It returns:
*** - 0   : if the exchange succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphHaloWait (
SCOTCH_DgraphHaloReq * const  requptr)
{
  return (dgraphHaloWait ((DgraphHaloRequest *) requptr));
}

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** halo request structure.
*** It returns:
*** - !NULL  : if the initialization succeeded.
*** - NULL   : on error.
+*/

SCOTCH_DgraphHaloReq *
SCOTCH_dgraphHaloReqAlloc ()
{
  return ((SCOTCH_DgraphHaloReq *) memAlloc (sizeof (SCOTCH_DgraphHaloReq)));
}
