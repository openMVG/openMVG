/* Copyright 2007-2009 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_fold_dup.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module folds a distributed graph   **/
/**                into two distinct copies, which may     **/
/**                be different when the number of         **/
/**                processes is odd.                       **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 aug 2006     **/
/**                                 to   : 20 jun 2007     **/
/**                # Version 5.1  : from : 14 nov 2008     **/
/**                                 to   : 28 oct 2009     **/
/**                # Version 6.0  : from : 28 sep 2014     **/
/**                                 to   : 28 sep 2014     **/
/**                                                        **/
/************************************************************/

#define DGRAPH
#define DGRAPH_FOLD_DUP

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_fold_dup.h"

/******************************/
/*                            */
/* This routine handles       */
/* distributed source graphs. */
/*                            */
/******************************/

/* This routine builds two folded graphs of
** a given graph on each of the two halves
** of the processes. The number of processes
** does not need to be even. There is a
** multi-threaded version, as well as a
** sequential one.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

#ifdef SCOTCH_PTHREAD
static
void *
dgraphFoldDup2 (
void * const                    dataptr)          /* Pointer to thread data */
{
  DgraphFoldDupData *             fldthrdptr;

  fldthrdptr = (DgraphFoldDupData *) dataptr;

  return ((void *) (intptr_t) dgraphFold2 (fldthrdptr->orggrafptr, fldthrdptr->partval, fldthrdptr->fldgrafptr,
                                           fldthrdptr->fldproccomm, fldthrdptr->orgdataptr,
                                           fldthrdptr->flddataptr, fldthrdptr->datatype));
}
#endif /* SCOTCH_PTHREAD */

int
dgraphFoldDup (
const Dgraph * restrict const orggrafptr,
Dgraph * restrict const       fldgrafptr,
void * restrict const         orgdataptr,         /* Un-based array of data which must be folded, e.g. coarmulttab */
void ** restrict const        flddataptr,         /* Un-based array of data which must be folded, e.g. coarmulttab */
MPI_Datatype                  datatype)
{
  int                       fldprocnbr;
  int                       fldprocnum;
  int                       fldproccol;
  MPI_Comm                  fldproccommtab[2];
#ifdef SCOTCH_PTHREAD
  Dgraph                    orggrafdat;
  DgraphFoldDupData         fldthrdtab[2];
  pthread_t                 thrdval;              /* Data of second thread */
#endif /* SCOTCH_PTHREAD */
  int                       o;

  fldprocnbr = (orggrafptr->procglbnbr + 1) / 2;  /* Median cut on number of processors     */
  if (orggrafptr->proclocnum < fldprocnbr) {      /* Compute color and rank in two subparts */
    fldproccol = 0;
    fldprocnum = orggrafptr->proclocnum;
    fldproccommtab[1] = MPI_COMM_NULL;
  }
  else {
    fldproccol = 1;
    fldprocnum = orggrafptr->proclocnum - fldprocnbr;
    fldproccommtab[0] = MPI_COMM_NULL;
  }
  if (MPI_Comm_split (orggrafptr->proccomm, fldproccol, fldprocnum, &fldproccommtab[fldproccol]) != MPI_SUCCESS) {
    errorPrint  ("dgraphFoldDup: communication error (1)");
    return      (1);
  }

#ifdef SCOTCH_PTHREAD
  orggrafdat = *orggrafptr;                       /* Create a separate graph structure to change its communicator */

  fldthrdtab[0].orggrafptr  = orggrafptr;
  fldthrdtab[0].fldgrafptr  = fldgrafptr;
  fldthrdtab[0].fldproccomm = fldproccommtab[0];
  fldthrdtab[0].partval     = 0;
  fldthrdtab[0].orgdataptr  = orgdataptr;
  fldthrdtab[0].flddataptr  = flddataptr;
  fldthrdtab[0].datatype    = datatype;
  fldthrdtab[1].orggrafptr  = &orggrafdat;
  fldthrdtab[1].fldgrafptr  = fldgrafptr;
  fldthrdtab[1].fldproccomm = fldproccommtab[1];
  fldthrdtab[1].partval     = 1;
  fldthrdtab[1].orgdataptr  = orgdataptr;
  fldthrdtab[1].flddataptr  = flddataptr;
  fldthrdtab[1].datatype    = datatype;

  if (MPI_Comm_dup (orggrafptr->proccomm, &orggrafdat.proccomm) != MPI_SUCCESS) { /* Duplicate communicator to avoid interferences in communications */
    errorPrint ("dgraphFoldDup: communication error (2)");
    return     (1);
  }

  if (pthread_create (&thrdval, NULL, dgraphFoldDup2, (void *) &fldthrdtab[1]) != 0) /* If could not create thread */
    o = (int) (intptr_t) dgraphFold2 (orggrafptr, 0, fldgrafptr, fldproccommtab[0], orgdataptr, flddataptr, datatype) || /* Call routines in sequence */
        (int) (intptr_t) dgraphFold2 (orggrafptr, 1, fldgrafptr, fldproccommtab[1], orgdataptr, flddataptr, datatype);
  else {                                          /* Newly created thread is processing subgraph 1, so let's process subgraph 0 */
    void *                    o2;

    o = (int) (intptr_t) dgraphFoldDup2 ((void *) &fldthrdtab[0]); /* Work on copy with private communicator */

    pthread_join (thrdval, &o2);
    o |= (int) (intptr_t) o2;
  }
  MPI_Comm_free (&orggrafdat.proccomm);

#else /* SCOTCH_PTHREAD */
  o = (dgraphFold2 (orggrafptr, 0, fldgrafptr, fldproccommtab[0], orgdataptr, flddataptr, datatype) || /* Call routines in sequence */
       dgraphFold2 (orggrafptr, 1, fldgrafptr, fldproccommtab[1], orgdataptr, flddataptr, datatype));
#endif /* SCOTCH_PTHREAD */

  fldgrafptr->prockeyval = fldproccol;            /* Discriminate between folded communicators at same level */

  return (o);
}
