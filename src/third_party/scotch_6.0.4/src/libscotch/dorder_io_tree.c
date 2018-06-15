/* Copyright 2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dorder_io_tree.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 26 jul 2007     **/
/**                                 to     26 jul 2007     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to     30 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "comm.h"
#include "dgraph.h"
#include "order.h"
#include "dorder.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine saves a distributed ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dorderSaveTree2 (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
FILE * restrict const         stream,
int                        (* funcptr) (const Order * const, const Gnum * const, FILE * const))
{
  Order                 corddat;                  /* Centralized ordering for tree structure */
  Gnum * restrict       vlbltab;
  int                   procglbnbr;
  int                   protnum;
  int                   reduloctab[3];
  int                   reduglbtab[3];
  int                   cheklocval;
  int                   chekglbval;

  if (stream != NULL) {                           /* If file provided         */
    reduloctab[0] = 1;                            /* This process is the root */
    reduloctab[1] = ordeptr->proclocnum;          /* Get its rank             */
  }
  else {
    reduloctab[0] =                               /* This process is not the root */
    reduloctab[1] = 0;
  }
  reduloctab[2] = (grafptr->vlblloctax != NULL) ? 1 : 0; /* See if vertex labels provided */
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, MPI_INT, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderSaveTree2: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("dorderSaveTree2: should have only one root");
    return     (1);
  }
  MPI_Comm_size (ordeptr->proccomm, &procglbnbr);
  if ((reduglbtab[2] != 0) && (reduglbtab[2] != procglbnbr)) {
    errorPrint ("dorderSaveTree2: inconsistent parameters");
    return     (1);
  }
  protnum = (int) reduglbtab[1];                  /* Get rank of root process */

  cheklocval = 0;
  vlbltab = NULL;
  if (reduglbtab[2] != 0) {
    if (protnum == ordeptr->proclocnum)
      if ((vlbltab = memAlloc (ordeptr->vnodglbnbr * sizeof (Gnum))) == NULL) {
        errorPrint ("dorderSaveTree2: out of memory");
        cheklocval = 1;
      }
#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
    if (MPI_Bcast (&cheklocval, 1, MPI_INT, protnum, ordeptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dorderSaveTree2: communication error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER1 */
    if (cheklocval != 0)
      return (1);
    if (commGatherv (grafptr->vlblloctax + grafptr->baseval, grafptr->vertlocnbr, GNUM_MPI,
                     vlbltab, grafptr->proccnttab, grafptr->procdsptab, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dorderSaveTree2: communication error (3)");
      return     (1);
    }
  }

  if (protnum == ordeptr->proclocnum)
    cheklocval = orderInit (&corddat, ordeptr->baseval, ordeptr->vnodglbnbr, NULL);
#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Bcast (&cheklocval, 1, MPI_INT, protnum, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderSaveTree2: communication error (4)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DORDER1 */
  if (cheklocval != 0)
    return (1);

  if (protnum == ordeptr->proclocnum) {
    cheklocval = dorderGather (ordeptr, &corddat); /* Need inverse permutation too */
    if (cheklocval == 0)
      cheklocval = funcptr (&corddat, vlbltab, stream);
    orderExit (&corddat);
  }
  else
    cheklocval = dorderGather (ordeptr, NULL);

  if (vlbltab != NULL)
    memFree (vlbltab);

#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderSaveTree2: communication error (3)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DORDER1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DORDER1 */

  return  (chekglbval);
}

/* This routine saves the separator tree
** data of the given distributed ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dorderSaveTree (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
FILE * restrict const         stream)
{
  return (dorderSaveTree2 (ordeptr, grafptr, stream, orderSaveTree));
}

/* This routine saves the column block
** mapping data of the given distributed
** ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dorderSaveMap (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
FILE * restrict const         stream)
{
  return (dorderSaveTree2 (ordeptr, grafptr, stream, orderSaveMap));
}
