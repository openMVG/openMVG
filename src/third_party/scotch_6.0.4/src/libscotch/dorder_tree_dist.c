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
/**   NAME       : dorder_tree_dist.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 28 nov 2007     **/
/**                                 to     09 may 2008     **/
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

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This function returns to all processes the
** number of distributed leaf column blocks
** possessed by the ordering.
** It returns:
** - >=0 : number of distributed column blocks.
** - <0  : on error.
*/

Gnum
dorderCblkDist (
const Dorder * restrict const ordeptr)
{
  const DorderLink * restrict linklocptr;
  Gnum                        dblklocnbr;         /* Local number of locally-rooted distributed column blocks */
  Gnum                        dblkglbnbr;

  for (linklocptr = ordeptr->linkdat.nextptr, dblklocnbr = 0; /* For all nodes in local ordering structure */
       linklocptr != &ordeptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST */
    if (cblklocptr->cblknum.proclocnum == ordeptr->proclocnum)
      dblklocnbr ++;
  }

  if (MPI_Allreduce (&dblklocnbr, &dblkglbnbr, 1, GNUM_MPI, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderCblkDist: communication error");
    return     ((Gnum) -1);
  }

  return (dblkglbnbr);
}

/* This function returns on all of the procesors the
** distributed part of the distributed structure of
** the given distributed ordering. The two array
** pointers which must be passed should both point to
** arrays of size dorderCblkDist().
** It returns:
** - 0   : if the distributed tree structure could be computed.
** - !0  : on error.
*/

int
dorderTreeDist (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
Gnum * restrict const         treeglbtab,
Gnum * restrict const         sizeglbtab)
{
  const DorderLink * restrict linklocptr;
  Gnum * restrict             dataloctab;
  Gnum * restrict             dataglbtab;
  Gnum                        dblklocnum;
  Gnum                        dblklocnbr;         /* Local number of distributed column blocks  */
  Gnum                        dblkglbnbr;         /* Global number of distributed column blocks */
  Gnum                        dblkglbnum;
  Gnum                        dblkglbtmp;
  int * restrict              dblkcnttab;
  int * restrict              dblkdsptab;
  int * restrict              cblkdsptab;
  Gnum                        cblkglbtmp;
  Gnum * restrict             srt1glbtab;
  Gnum * restrict             srt2glbtab;
  int                         procglbnbr;
  int                         procnum;
  Gnum                        reduloctab[3];
  Gnum                        reduglbtab[3];

  for (linklocptr = ordeptr->linkdat.nextptr, dblklocnbr = 0; /* For all nodes in local ordering structure */
       linklocptr != &ordeptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST */
    if (cblklocptr->cblknum.proclocnum == ordeptr->proclocnum) {
#ifdef SCOTCH_DEBUG_DORDER2
      Gnum                        cblklocnum;

      cblklocnum = cblklocptr->cblknum.cblklocnum;
      if ((cblklocnum < 0) || (cblklocnum >= ordeptr->cblklocnbr)) {
        errorPrint ("dorderTreeDist: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DORDER2 */
      dblklocnbr ++;
    }
  }
  if (MPI_Allreduce (&dblklocnbr, &dblkglbnbr, 1, GNUM_MPI, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) { /* Get overall number of distributed blocks */
    errorPrint ("dorderTreeDist: communication error (1)");
    return     (1);
  }

  MPI_Comm_size (ordeptr->proccomm, &procglbnbr);

  reduloctab[0] =
  reduloctab[1] =
  reduloctab[2] = 0;
  if (memAllocGroup ((void **) (void *)
                     &dblkcnttab, (size_t) ( procglbnbr      * sizeof (int)),
                     &dblkdsptab, (size_t) ( procglbnbr      * sizeof (int)), /* TRICK: cblkdsptab used as secondary array after cblkcnttab */
                     &cblkdsptab, (size_t) ((procglbnbr + 1) * sizeof (int)), /* TRICK: have an array at least of size 2                    */
                     &dataloctab, (size_t) ( dblklocnbr * 4  * sizeof (Gnum)),
                     &dataglbtab, (size_t) ( dblkglbnbr * 4  * sizeof (Gnum)),
                     &srt1glbtab, (size_t) ( dblkglbnbr * 2  * sizeof (Gnum)), /* TRICK: one more slot for root node                  */
                     &srt2glbtab, (size_t) ( dblkglbnbr * 2  * sizeof (Gnum)), NULL) == NULL) { /* TRICK: one more slot for root node */
    errorPrint ("dorderTreeDist: out of memory");
    reduloctab[0] = 1;                            /* Memory error */
  }
  else {
    if (treeglbtab != NULL)
      reduloctab[1] = 1;                          /* Compute the "or" of any array being non-null */
    if (sizeglbtab != NULL) {
      reduloctab[2] = reduloctab[1];              /* Compute the "and" of any array being non-null */
      reduloctab[1] = 1;
    }
  }
#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderTreeDist: communication error (1)");
    reduglbtab[0] =                               /* Post-process error below      */
    reduglbtab[1] =                               /* Prevent Valgrind from yelling */
    reduglbtab[2] = 1;
  }
#else /* SCOTCH_DEBUG_DORDER1 */
  reduglbtab[0] = reduloctab[0];
  reduglbtab[1] = procglbnbr - 1 + reduloctab[1];
  reduglbtab[2] = procglbnbr - 1 + reduloctab[2];
#endif /* SCOTCH_DEBUG_DORDER1 */

  if (reduglbtab[1] != reduglbtab[2]) {           /* If not both arrays provided on each of the candidate processors */
    if (reduloctab[1] != reduloctab[2])
      errorPrint ("dorderTreeDist: invalid parameters (1)");
    reduglbtab[0] = 1;
  }
  if (reduglbtab[2] != procglbnbr) {
    errorPrint ("dorderTreeDist: invalid parameters (2)");
    reduglbtab[0] = 1;
  }
  if (reduglbtab[0] != 0) {
    if (dblkcnttab != NULL)
      memFree (dblkcnttab);                       /* Free group leader */
    return (1);
  }

  cblkdsptab[0] = (int) dblklocnbr;               /* MPI only supports int as count type     */
  cblkdsptab[1] = (int) ordeptr->cblklocnbr;      /* TRICK: cblkdsptab is at least of size 2 */
  if (MPI_Allgather (cblkdsptab, 2, MPI_INT, dblkcnttab, 2, MPI_INT, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderTreeDist: communication error (2)");
    return     (1);
  }
  for (procnum = cblkglbtmp = 0; procnum < procglbnbr; procnum ++) { /* Accumulate un-based global start indices for all column blocks */
    cblkdsptab[procnum] = cblkglbtmp;
    dblkcnttab[procnum] = dblkcnttab[2 * procnum] * 4; /* Four times for dataloctab */
    cblkglbtmp         += dblkcnttab[2 * procnum + 1];
  }
  for (procnum = dblkglbtmp = 0; procnum < procglbnbr; procnum ++) { /* Accumulate un-based global start indices for distributed column blocks */
    dblkdsptab[procnum] = dblkglbtmp;
    dblkglbtmp         += dblkcnttab[procnum];
  }

  for (linklocptr = ordeptr->linkdat.nextptr, dblklocnum = 0; /* For all nodes in local ordering structure */
       linklocptr != &ordeptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST                    */
    if (cblklocptr->cblknum.proclocnum == ordeptr->proclocnum) { /* If node is local */
      dataloctab[4 * dblklocnum]     = cblkdsptab[ordeptr->proclocnum] + cblklocptr->cblknum.cblklocnum;
      dataloctab[4 * dblklocnum + 1] = cblklocptr->ordeglbval;
      dataloctab[4 * dblklocnum + 2] = cblkdsptab[cblklocptr->fathnum.proclocnum] + cblklocptr->fathnum.cblklocnum;
      dataloctab[4 * dblklocnum + 3] = cblklocptr->vnodglbnbr;
      dblklocnum ++;
    }
  }
  if (MPI_Allgatherv (dataloctab, 4 * dblklocnbr, GNUM_MPI, dataglbtab, dblkcnttab, dblkdsptab, GNUM_MPI, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderTreeDist: communication error (3)");
    return     (1);
  }

  for (dblkglbnum = 0; dblkglbnum < dblkglbnbr; dblkglbnum ++) {
    srt1glbtab[2 * dblkglbnum]     = dataglbtab[4 * dblkglbnum + 1];
    srt1glbtab[2 * dblkglbnum + 1] = dataglbtab[4 * dblkglbnum];
  }
  intSort2asc2 (srt1glbtab, dblkglbnbr);          /* Sort nodes by ascending inverse start index to get permutation of column block indices */
  for (dblkglbnum = 0; dblkglbnum < dblkglbnbr; dblkglbnum ++) {
    srt1glbtab[2 * dblkglbnum]     = srt1glbtab[2 * dblkglbnum + 1];
    srt1glbtab[2 * dblkglbnum + 1] = dblkglbnum;
  }
  intSort2asc2 (srt1glbtab, dblkglbnbr);          /* Sort nodes by ascending column block index to match with the ones of dataglbtab */

  for (dblkglbnum = 0; dblkglbnum < dblkglbnbr; dblkglbnum ++) {
    srt2glbtab[2 * dblkglbnum]     = dataglbtab[4 * dblkglbnum + 2];
    srt2glbtab[2 * dblkglbnum + 1] = dblkglbnum;
  }
  intSort2asc2 (srt2glbtab, dblkglbnbr);          /* Sort father indices by ascending column block indices */
#ifdef SCOTCH_DEBUG_DORDER2
  if (srt2glbtab[0] != -1) {                      /* If tree has no root */
    errorPrint ("dorderTreeDist: internal error (2)");
    memFree    (dblkcnttab);                      /* Free group leader */
    return     (1);
  }
  if ((dblkglbnbr > 1) && (srt2glbtab[2] == -1)) { /* If tree has multiple roots */
    errorPrint ("dorderTreeDist: internal error (3)");
    memFree    (dblkcnttab);                      /* Free group leader */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DORDER2 */
  for (dblkglbnum = 1, dblkglbtmp = 0; dblkglbnum < dblkglbnbr; ) { /* Replace in block data the father column block indices by the new permuted indices */
    if (srt2glbtab[2 * dblkglbnum] == srt1glbtab[2 * dblkglbtmp])
      dataglbtab[4 * srt2glbtab[2 * (dblkglbnum ++) + 1] + 2] = srt1glbtab[2 * dblkglbtmp + 1];
    else {
#ifdef SCOTCH_DEBUG_DORDER2
      if ((srt2glbtab[2 * dblkglbnum] < srt1glbtab[2 * dblkglbtmp]) || /* If column block index not found in table */
          (dblkglbtmp >= (dblkglbnbr - 1))) {
        errorPrint ("dorderTreeDist: internal error (4)");
        memFree    (dblkcnttab);                  /* Free group leader */
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DORDER2 */
      dblkglbtmp ++;
    }
  }

  for (dblkglbnum = 0; dblkglbnum < dblkglbnbr; dblkglbnum ++) {
    srt2glbtab[2 * dblkglbnum]     = dataglbtab[4 * dblkglbnum];
    srt2glbtab[2 * dblkglbnum + 1] = dblkglbnum;
  }
  intSort2asc2 (srt2glbtab, dblkglbnbr);          /* Sort father indices by ascending column block indices */
  for (dblkglbnum = 0; dblkglbnum < dblkglbnbr; dblkglbnum ++) {
#ifdef SCOTCH_DEBUG_DORDER2
    if (srt1glbtab[2 * dblkglbnum] != srt2glbtab[2 * dblkglbnum]) {
      errorPrint ("dorderTreeDist: internal error (5)");
      memFree    (dblkcnttab);                    /* Free group leader */
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
    treeglbtab[srt1glbtab[2 * dblkglbnum + 1]] = dataglbtab[4 * srt2glbtab[2 * dblkglbnum + 1] + 2];
    sizeglbtab[srt1glbtab[2 * dblkglbnum + 1]] = dataglbtab[4 * srt2glbtab[2 * dblkglbnum + 1] + 3];
  }

  memFree (dblkcnttab);                           /* Free group leader */

  return (0);
}
