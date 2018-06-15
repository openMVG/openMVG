/* Copyright 2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_allreduce.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the distributed source  **/
/**                graph building routines.                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 28 aug 2006     **/
/**                                 to   : 28 aug 2006     **/
/**                                                        **/
/************************************************************/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"

/* This routine creates an allreduce operator from
** the function pointer that is passed to it, and
** use it to perform an allreduce operation on the
** given data.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphAllreduceMaxSum2 (
Gnum *                      reduloctab,           /* Pointer to array of local Gnum data   */
Gnum *                      reduglbtab,           /* Pointer to array of reduced Gnum data */
int                         redumaxsumnbr,        /* Number of max + sum Gnum operations   */
MPI_User_function *         redufuncptr,          /* Pointer to operator function          */
MPI_Comm                    proccomm)             /* Communicator to be used for reduction */
{
  MPI_Datatype      redutypedat;                  /* Data type for finding best separator              */
  MPI_Op            reduoperdat;                  /* Handle of MPI operator for finding best separator */

  if ((MPI_Type_contiguous (redumaxsumnbr, GNUM_MPI, &redutypedat) != MPI_SUCCESS) ||
      (MPI_Type_commit (&redutypedat)                              != MPI_SUCCESS) ||
      (MPI_Op_create (redufuncptr, 1, &reduoperdat)                != MPI_SUCCESS)) {
    errorPrint ("dgraphAllreduceMaxSum: communication error (1)");
    return     (1);
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 1, redutypedat, reduoperdat, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphAllreduceMaxSum: communication error (2)");
    return     (1);
  }

  if ((MPI_Op_free   (&reduoperdat) != MPI_SUCCESS) ||
      (MPI_Type_free (&redutypedat) != MPI_SUCCESS)) {
    errorPrint ("dgraphAllreduceMaxSum: communication error (3)");
    return     (1);
  }

  return (0);
}
