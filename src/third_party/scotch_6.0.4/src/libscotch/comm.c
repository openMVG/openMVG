/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : comm.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the large size     **/
/**                communication handling routines.        **/
/**                                                        **/
/**    DATES     : # Version 5.1  : from : 30 jul 2010     **/
/**                                 to   : 30 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define COMM

#include "module.h"
#include "common.h"
#include "comm.h"

/************************************/
/*                                  */
/* These routines handle large size */
/* communications                   */
/*                                  */
/************************************/

/*
**
*/

int
commAllgatherv (
void * const                senddattab,
const Gnum                  sendcntnbr,
MPI_Datatype                sendtypval,
void * const                recvdattab,
const Gnum * const          recvcnttab,
const Gnum * const          recvdsptab,
MPI_Datatype                recvtypval,
MPI_Comm                    comm)
{
  int * restrict      ircvcnttab;
  int * restrict      ircvdsptab;
  int                 procglbnbr;
  int                 procnum;
  int                 o;

  MPI_Comm_size (comm, &procglbnbr);
  if (memAllocGroup ((void **) (void *)
                     &ircvcnttab, (size_t) (procglbnbr * sizeof (int)),
                     &ircvdsptab, (size_t) (procglbnbr * sizeof (int)), NULL) == NULL) {
    errorPrint ("commAllgatherv: out of memory");
    return     (MPI_ERR_OTHER);
  }

  for (procnum = 0; procnum < procglbnbr; procnum ++) {
    ircvcnttab[procnum] = (int) recvcnttab[procnum];
    ircvdsptab[procnum] = (int) recvdsptab[procnum];
    if (((Gnum) ircvcnttab[procnum] != recvcnttab[procnum]) ||
        ((Gnum) ircvdsptab[procnum] != recvdsptab[procnum])) {
      errorPrint ("commAllgatherv: communication indices out of range");
      memFree    (ircvcnttab);
      return     (MPI_ERR_ARG);
    }
  }

  o = MPI_Allgatherv (senddattab, sendcntnbr, sendtypval,
                      recvdattab, ircvcnttab, ircvdsptab, recvtypval, comm);

  memFree (ircvcnttab);

  return (o);
}

/*
**
*/

int
commGatherv (
void * const                senddattab,
const Gnum                  sendcntnbr,
MPI_Datatype                sendtypval,
void * const                recvdattab,
const Gnum * const          recvcnttab,
const Gnum * const          recvdsptab,
MPI_Datatype                recvtypval,
const int                   rootnum,
MPI_Comm                    comm)
{
  int * restrict      ircvcnttab;
  int * restrict      ircvdsptab;
  int                 proclocnum;
  int                 o;

  MPI_Comm_rank (comm, &proclocnum);

  ircvcnttab = NULL;

  if (rootnum == proclocnum) {
    int                 procglbnbr;
    int                 procnum;

    MPI_Comm_size (comm, &procglbnbr);
    if (memAllocGroup ((void **) (void *)
                       &ircvcnttab, (size_t) (procglbnbr * sizeof (int)),
                       &ircvdsptab, (size_t) (procglbnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("commGatherv: out of memory");
      return     (MPI_ERR_OTHER);
    }

    for (procnum = 0; procnum < procglbnbr; procnum ++) {
      ircvcnttab[procnum] = (int) recvcnttab[procnum];
      ircvdsptab[procnum] = (int) recvdsptab[procnum];
      if (((Gnum) ircvcnttab[procnum] != recvcnttab[procnum]) ||
          ((Gnum) ircvdsptab[procnum] != recvdsptab[procnum])) {
        errorPrint ("commGatherv: communication indices out of range");
        memFree    (ircvcnttab);
        return     (MPI_ERR_ARG);
      }
    }
  }

  o = MPI_Gatherv (senddattab, sendcntnbr, sendtypval,
                   recvdattab, ircvcnttab, ircvdsptab, recvtypval, rootnum, comm);

  if (ircvcnttab != NULL)
    memFree (ircvcnttab);

  return (o);
}

/*
**
*/

int
commScatterv (
void * const                senddattab,
const Gnum * const          sendcnttab,
const Gnum * const          senddsptab,
MPI_Datatype                sendtypval,
void * const                recvdattab,
const Gnum                  recvcntnbr,
MPI_Datatype                recvtypval,
const int                   rootnum,
MPI_Comm                    comm)
{
  int * restrict      isndcnttab;
  int * restrict      isnddsptab;
  int                 proclocnum;
  int                 o;

  MPI_Comm_rank (comm, &proclocnum);

  isndcnttab = NULL;

  if (rootnum == proclocnum) {
    int                 procglbnbr;
    int                 procnum;

    MPI_Comm_size (comm, &procglbnbr);
    if (memAllocGroup ((void **) (void *)
                       &isndcnttab, (size_t) (procglbnbr * sizeof (int)),
                       &isnddsptab, (size_t) (procglbnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("commScatterv: out of memory");
      return     (MPI_ERR_OTHER);
    }

    for (procnum = 0; procnum < procglbnbr; procnum ++) {
      isndcnttab[procnum] = (int) sendcnttab[procnum];
      isnddsptab[procnum] = (int) senddsptab[procnum];
      if (((Gnum) isndcnttab[procnum] != sendcnttab[procnum]) ||
          ((Gnum) isnddsptab[procnum] != senddsptab[procnum])) {
        errorPrint ("commScatterv: communication indices out of range");
        memFree    (isndcnttab);
        return     (MPI_ERR_ARG);
      }
    }
  }

  o = MPI_Scatterv (senddattab, isndcnttab, isnddsptab, sendtypval,
                    recvdattab, (int) recvcntnbr, recvtypval, rootnum, comm);

  if (isndcnttab != NULL)
    memFree (isndcnttab);

  return (o);
}
