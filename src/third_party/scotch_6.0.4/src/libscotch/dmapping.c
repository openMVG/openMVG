/* Copyright 2008,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dmapping.c                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                                                        **/
/**   FUNCTION   : This module handles (partial) mappings. **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 31 mar 2008     **/
/**                                 to     09 nov 2008     **/
/**                # Version 6.0  : from : 03 sep 2013     **/
/**                                 to     03 sep 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DMAPPING

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "arch.h"
#include "dmapping.h"

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine builds a mapping.
** It returns:
** - 0   : if mapping successfully initialized.
** - !0  : on error.
*/

int
dmapInit (
Dmapping * restrict const     dmapptr,
const Arch * restrict const   archptr)
{
  dmapptr->fragptr    = NULL;
  dmapptr->fragnbr    =
  dmapptr->vertlocmax =
  dmapptr->vertlocnbr = 0;
  dmapptr->archdat    = *archptr;

#ifdef SCOTCH_PTHREAD
  pthread_mutex_init (&dmapptr->mutelocdat, NULL); /* Initialize local mutex */
#endif /* SCOTCH_PTHREAD */

  return (0);
}

/* This routine frees the contents of the given
** mapping. The architecture data is never freed
** as it is usually a copy of an existing Arch
** structure.
** It returns:
** - VOID  : in all cases.
*/

void
dmapExit (
Dmapping * const             dmapptr)
{
  DmappingFrag *      fragptr;
  DmappingFrag *      fragtmp;

  for (fragptr = dmapptr->fragptr; fragptr != NULL; fragptr = fragtmp) {
    memFree (fragptr->vnumtab);
    memFree (fragptr->parttab);
    memFree (fragptr->domntab);
    fragtmp = fragptr->nextptr;
    memFree (fragptr);
  }

#ifdef SCOTCH_PTHREAD
  pthread_mutex_destroy (&dmapptr->mutelocdat);   /* Destroy local mutex */
#endif /* SCOTCH_PTHREAD */

#ifdef SCOTCH_DEBUG_DMAP2
  memSet (dmapptr, ~0, sizeof (Dmapping));
#endif /* SCOTCH_DEBUG_DMAP2 */
}

/* This routine adds a fragment to the given
** distributed mapping.
** It returns:
** - void  : in all cases.
*/

void
dmapAdd (
Dmapping * restrict const     dmapptr,
DmappingFrag * restrict const fragptr)
{
#ifdef SCOTCH_PTHREAD
  pthread_mutex_lock (&dmapptr->mutelocdat);      /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */

  if (dmapptr->vertlocmax < fragptr->vertnbr)
    dmapptr->vertlocmax = fragptr->vertnbr;
  dmapptr->vertlocnbr += fragptr->vertnbr;

  dmapptr->fragnbr ++;
  fragptr->nextptr = dmapptr->fragptr;            /* Link fragment to mapping */
  dmapptr->fragptr = fragptr;

#ifdef SCOTCH_PTHREAD
  pthread_mutex_unlock (&dmapptr->mutelocdat);    /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */
}

/* This routine propagates back distributed mapping
** information to a part array associated with a
** distributed graph structure.
** It returns:
** - 0   : if partition data successfully obtained.
** - !0  : on error.
*/

int
dmapTerm (
const Dmapping * restrict const dmapptr,
const Dgraph * restrict const   grafptr,
Gnum * restrict const           termloctab)
{
  Gnum * restrict             termloctax;
  int * restrict              sendcnttab;
  int * restrict              senddsptab;
  int * restrict              recvcnttab;
  int * restrict              recvdsptab;
  DmappingTermSort * restrict sortsndtab;
  DmappingTermSort * restrict sortrcvtab;
  Gnum                        vertlocnum;
  int                         vertrcvnbr;
  int                         vertsndnbr;
  int                         procnum;
  DmappingFrag * restrict     fragptr;
  Gnum                        reduloctab[2];
  Gnum                        reduglbtab[2];

  reduloctab[0] = dmapptr->vertlocnbr;
  reduloctab[1] = 0;
  if (memAllocGroup ((void **) (void *)
                     &senddsptab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &sendcnttab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &recvdsptab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &recvcnttab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &sortsndtab, (size_t) ((dmapptr->vertlocnbr + 1) * sizeof (DmappingTermSort)), /* "+1" for end marker */
                     &sortrcvtab, (size_t) (grafptr->vertlocnbr       * sizeof (DmappingTermSort)), NULL) == NULL) {
    errorPrint ("dmapTerm: out of memory");
    reduloctab[1] = 1;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 2, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmapTerm: communication error (1)");
    reduglbtab[1] = 1;
  }
  if (reduglbtab[1] != 0) {
    if (senddsptab != NULL)
      memFree (senddsptab);                       /* Free group leader */
    return (1);
  }

  if (reduglbtab[0] == 0) {                       /* If mapping structure is empty, create an empty mapping */
    memSet  (termloctab, 0, grafptr->vertlocnbr * sizeof (Gnum));
    memFree (senddsptab);                         /* Free group leader */
    return (0);
  }
  if (reduglbtab[0] != grafptr->vertglbnbr) {
    errorPrint ("dmapTerm: invalid mapping (1)");
    memFree    (senddsptab);                      /* Free group leader */
    return     (1);
  }

  for (fragptr = dmapptr->fragptr, vertlocnum = 0; fragptr != NULL; fragptr = fragptr->nextptr) {
    Gnum                      fraglocnum;

    for (fraglocnum = 0; fraglocnum < fragptr->vertnbr; fraglocnum ++, vertlocnum ++) {
#ifdef SCOTCH_DEBUG_DMAP2
      if ((vertlocnum >= dmapptr->vertlocnbr) || (fragptr->parttab[fraglocnum] < 0) || (fragptr->parttab[fraglocnum] >= fragptr->domnnbr)) {
        errorPrint ("dmapTerm: invalid mapping (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DMAP2 */
      sortsndtab[vertlocnum].vertnum = fragptr->vnumtab[fraglocnum];
      sortsndtab[vertlocnum].termnum = (Gnum) archDomNum (&dmapptr->archdat, &fragptr->domntab[fragptr->parttab[fraglocnum]]);
    }
  }
#ifdef SCOTCH_DEBUG_DMAP2
  if (vertlocnum != dmapptr->vertlocnbr) {
    errorPrint ("dmapTerm: invalid mapping (3)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DMAP2 */

  sortsndtab[vertlocnum].vertnum =                /* Set end marker */
  sortsndtab[vertlocnum].termnum = GNUMMAX;
  intSort2asc1 (sortsndtab, dmapptr->vertlocnbr); /* Sort mapping array by original vertex numbers, without marker */

  for (vertlocnum = 0, procnum = 0; procnum < grafptr->procglbnbr; ) {
    Gnum                  vertsndnbr;
    Gnum                  procvrtval;

    vertsndnbr = 0;
    procvrtval = grafptr->procvrttab[procnum + 1];
    while (sortsndtab[vertlocnum].vertnum < procvrtval) {
      vertsndnbr ++;
      vertlocnum ++;
#ifdef SCOTCH_DEBUG_DMAP2
      if (vertlocnum > dmapptr->vertlocnbr) {     /* If beyond regular indices plus end marker */
        errorPrint ("dmapTerm: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DMAP2 */
    }
    sendcnttab[procnum ++] = (int) (vertsndnbr * 2); /* Communication array for MPI, so (int), and "*2" because a Sort is 2 Gnums */
  }
#ifdef SCOTCH_DEBUG_DMAP2
  if (vertlocnum != dmapptr->vertlocnbr) {
    errorPrint ("dmapTerm: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DMAP2 */

  if (MPI_Alltoall (sendcnttab, 1, MPI_INT, recvcnttab, 1, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmapTerm: communication error (2)");
    return     (1);
  }

  for (procnum = 0, vertrcvnbr = vertsndnbr = 0; procnum < grafptr->procglbnbr; procnum ++) { /* Accumulate send and receive indices */
    recvdsptab[procnum] = vertrcvnbr;
    vertrcvnbr += recvcnttab[procnum];            /* Accumulate "*2" values as counts */
    senddsptab[procnum] = vertsndnbr;
    vertsndnbr += sendcnttab[procnum];
  }

  if (MPI_Alltoallv (sortsndtab, sendcnttab, senddsptab, GNUM_MPI, sortrcvtab, recvcnttab, recvdsptab, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmapTerm: communication error (3)");
    return     (1);
  }

  memSet (termloctab, ~0, grafptr->vertlocnbr * sizeof (Gnum));

  termloctax = termloctab - grafptr->procvrttab[grafptr->proclocnum]; /* Base local array through global indices */
  for (vertlocnum = 0; vertlocnum < grafptr->vertlocnbr; vertlocnum ++) {
#ifdef SCOTCH_DEBUG_DMAP2
    if (termloctax[sortrcvtab[vertlocnum].vertnum] != ~0) {
      errorPrint ("dmapTerm: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DMAP2 */
    termloctax[sortrcvtab[vertlocnum].vertnum] = sortrcvtab[vertlocnum].termnum;
  }
#ifdef SCOTCH_DEBUG_DMAP2
  for (vertlocnum = 0; vertlocnum < grafptr->vertlocnbr; vertlocnum ++) {
    if (termloctab[vertlocnum] == ~0) {
      errorPrint ("dmapTerm: internal error (4)");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DMAP2 */

  memFree (senddsptab);                           /* Free group leader */

  return (0);
}
