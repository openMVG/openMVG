/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dorder_gather.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 13 oct 2007     **/
/**                                 to     21 oct 2007     **/
/**                # Version 5.1  : from : 26 sep 2008     **/
/**                                 to     26 sep 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dorder.h"
#include "dorder_perm.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This function builds a distributed direct
** permutation from the information stored
** in the distributed ordering structure.
** It returns:
** - 0   : if the distributed permutation could be computed.
** - !0  : on error.
*/

int
dorderPerm (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
Gnum * restrict const         permloctab)
{
  Gnum * restrict             permloctax;
  int * restrict              sendcnttab;
  int * restrict              senddsptab;
  int * restrict              recvcnttab;
  int * restrict              recvdsptab;
  DorderPermSort * restrict   sortsndtab;
  DorderPermSort * restrict   sortrcvtab;
  const DorderLink * restrict linklocptr;
  Gnum                        vnodlocnbr;
  Gnum                        vnodlocnum;
  int                         vnodrcvnbr;
  int                         vnodsndnbr;
  int                         procnum;
  Gnum                        reduloctab[2];
  Gnum                        reduglbtab[2];

  for (linklocptr = ordeptr->linkdat.nextptr, vnodlocnbr = 0; /* For all nodes in local ordering structure */
       linklocptr != &ordeptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST               */
    if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0) /* If node is leaf         */
      vnodlocnbr += cblklocptr->data.leaf.vnodlocnbr; /* And more node vertices */
#ifdef SCOTCH_DEBUG_DORDER2
    else if (cblklocptr->typeval != DORDERCBLKNEDI) {
      errorPrint ("dorderPerm: invalid parameters (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
  }

  reduloctab[0] = vnodlocnbr;
  reduloctab[1] = 0;
  if (memAllocGroup ((void **) (void *)
                     &senddsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &sendcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &recvdsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &recvcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &sortsndtab, (size_t) ((vnodlocnbr + 1)    * sizeof (DorderPermSort)), /* "+1" for end marker */
                     &sortrcvtab, (size_t) (grafptr->vertlocnbr * sizeof (DorderPermSort)), NULL) == NULL) {
    errorPrint ("dorderPerm: out of memory");
    reduloctab[1] = 1;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 2, GNUM_MPI, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderPerm: communication error (1)");
    reduglbtab[1] = 1;
  }
  if (reduglbtab[1] != 0) {
    if (senddsptab != NULL)
      memFree (senddsptab);                       /* Free group leader */
    return (1);
  }

  if (reduglbtab[0] == 0) {                       /* If ordering structure is empty */
    Gnum                  ordelocval;             /* Based permutation start index  */

    memFree (senddsptab);                         /* Free group leader */

    for (vnodlocnum = 0, ordelocval = grafptr->procvrttab[grafptr->proclocnum]; /* Build identity permutation */
         vnodlocnum < grafptr->vertlocnbr; vnodlocnum ++)
      permloctab[vnodlocnum] = ordelocval ++;

    return (0);
  }
  if (reduglbtab[0] != grafptr->vertglbnbr) {
    errorPrint ("dorderPerm: invalid parameters (2)");
    memFree    (senddsptab);                      /* Free group leader */
    return     (1);
  }

  for (linklocptr = ordeptr->linkdat.nextptr, vnodlocnum = 0; /* For all nodes in local ordering structure */
       linklocptr != &ordeptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST         */
    if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0) { /* If node is leaf */
      Gnum                  leaflocnbr;
      Gnum                  leaflocnum;
      Gnum                  ordelocval;           /* Based permutation start index */

      for (leaflocnum = 0, leaflocnbr = cblklocptr->data.leaf.vnodlocnbr, ordelocval = cblklocptr->data.leaf.ordelocval + ordeptr->baseval;
           leaflocnum < leaflocnbr; leaflocnum ++, vnodlocnum ++) {
        sortsndtab[vnodlocnum].vertnum = cblklocptr->data.leaf.periloctab[leaflocnum];
        sortsndtab[vnodlocnum].permnum = ordelocval + leaflocnum;
#ifdef SCOTCH_DEBUG_DORDER2
        if ((sortsndtab[vnodlocnum].vertnum <  ordeptr->baseval) ||
            (sortsndtab[vnodlocnum].vertnum > (ordeptr->baseval + ordeptr->vnodglbnbr)) ||
            (sortsndtab[vnodlocnum].permnum <  ordeptr->baseval) ||
            (sortsndtab[vnodlocnum].permnum > (ordeptr->baseval + ordeptr->vnodglbnbr))) {
          errorPrint ("dorderPerm: internal error (1)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DORDER2 */
      }
    }
  }
  sortsndtab[vnodlocnbr].vertnum =                /* Set end marker */
  sortsndtab[vnodlocnbr].permnum = GNUMMAX;
  intSort2asc1 (sortsndtab, vnodlocnbr);          /* Sort permutation array by original vertex numbers, without marker */

  for (vnodlocnum = 0, procnum = 0; procnum < grafptr->procglbnbr; ) {
    Gnum                  vnodsndnbr;
    Gnum                  procdspval;

    vnodsndnbr = 0;
    procdspval = grafptr->procdsptab[procnum + 1];
    while (sortsndtab[vnodlocnum].vertnum < procdspval) {
      vnodsndnbr ++;
      vnodlocnum ++;
#ifdef SCOTCH_DEBUG_DORDER2
      if (vnodlocnum > vnodlocnbr) {              /* If beyond regular indices plus end marker */
        errorPrint ("dorderPerm: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DORDER2 */
    }
    sendcnttab[procnum ++] = (int) (vnodsndnbr * 2); /* Communication array for MPI, so (int), and "*2" because a Sort is 2 Gnums */
  }
#ifdef SCOTCH_DEBUG_DORDER2
  if (vnodlocnum != vnodlocnbr) {
    errorPrint ("dorderPerm: internal error (3)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DORDER2 */

  if (MPI_Alltoall (sendcnttab, 1, MPI_INT, recvcnttab, 1, MPI_INT, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderPerm: communication error (2)");
    return     (1);
  }

  for (procnum = 0, vnodrcvnbr = vnodsndnbr = 0; procnum < grafptr->procglbnbr; procnum ++) { /* Accumulate send and receive indices */
    recvdsptab[procnum] = vnodrcvnbr;
    vnodrcvnbr += recvcnttab[procnum];            /* Accumulate "*2" values as counts */
    senddsptab[procnum] = vnodsndnbr;
    vnodsndnbr += sendcnttab[procnum];
  }

  if (MPI_Alltoallv (sortsndtab, sendcnttab, senddsptab, GNUM_MPI, sortrcvtab, recvcnttab, recvdsptab, GNUM_MPI, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderPerm: communication error (3)");
    return     (1);
  }

#ifdef SCOTCH_DEBUG_DORDER2
  memSet (permloctab, ~0, grafptr->vertlocnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_DORDER2 */

  permloctax = permloctab - grafptr->procdsptab[grafptr->proclocnum]; /* Base local array through global indices */
  for (vnodlocnum = 0; vnodlocnum < grafptr->vertlocnbr; vnodlocnum ++) {
#ifdef SCOTCH_DEBUG_DORDER2
    if (permloctax[sortrcvtab[vnodlocnum].vertnum] != ~0) {
      errorPrint ("dorderPerm: internal error (4)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
    permloctax[sortrcvtab[vnodlocnum].vertnum] = sortrcvtab[vnodlocnum].permnum;
  }
#ifdef SCOTCH_DEBUG_DORDER2
  for (vnodlocnum = 0; vnodlocnum < grafptr->vertlocnbr; vnodlocnum ++) {
    if (permloctab[vnodlocnum] == ~0) {
      errorPrint ("dorderPerm: internal error (5)");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DORDER2 */

  memFree (senddsptab);                           /* Free group leader */

  return (0);
}
