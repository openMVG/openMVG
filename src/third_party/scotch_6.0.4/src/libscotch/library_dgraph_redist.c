/* Copyright 2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_redist.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the data declara-    **/
/**                tions for the graph redistribution      **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 28 mar 2012     **/
/**                                 to   : 29 nov 2012     **/
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
#include "dgraph_redist.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/*+ This routine computes a distributed graph
*** that matches the provided partition.
*** It returns:
*** - 0   : if redistributed graph created.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphRedist (
SCOTCH_Dgraph * const       orggrafptr,
const SCOTCH_Num * const    partloctab,           /* Array of part numbers for each local vertex */
const SCOTCH_Num * const    permgsttab,           /* Redistribution permutation array            */
const SCOTCH_Num            vertlocdlt,           /* Extra size of local vertex array            */
const SCOTCH_Num            edgelocdlt,           /* Extra size of local edge array              */
SCOTCH_Dgraph * const       redgrafptr)
{
  SCOTCH_Num          baseval;
#ifdef SCOTCH_DEBUG_LIBRARY1
  int                 o;

  MPI_Comm_compare (((Dgraph * restrict const) orggrafptr)->proccomm,
                    ((Dgraph * restrict const) redgrafptr)->proccomm, &o);
  if ((o != MPI_IDENT) && (o != MPI_CONGRUENT)) {
    errorPrint ("SCOTCH_dgraphRedist: communicators are not congruent");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  baseval = ((Dgraph *) orggrafptr)->baseval;

  return (dgraphRedist ((Dgraph *) orggrafptr,
                        ((partloctab != NULL) && (partloctab != (SCOTCH_Num *) orggrafptr)) ? (const Gnum * restrict const) (partloctab - baseval) : NULL,
                        ((permgsttab != NULL) && (permgsttab != (SCOTCH_Num *) orggrafptr)) ? (const Gnum * restrict const) (permgsttab - baseval) : NULL,
                        (vertlocdlt < 0) ? 0 : vertlocdlt, (edgelocdlt < 0) ? 0 : edgelocdlt, (Dgraph *) redgrafptr));
}
