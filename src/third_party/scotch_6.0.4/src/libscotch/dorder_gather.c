/* Copyright 2007,2008,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   DATES      : # Version 5.0  : from : 19 jul 2007     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 28 sep 2008     **/
/**                                 to     28 sep 2008     **/
/**                # Version 6.0  : from : 10 oct 2013     **/
/**                                 to     10 oct 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_allreduce.h"
#include "order.h"
#include "dorder.h"
#include "dorder_gather.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This function gathers the pieces of
** a distributed ordering to build a
** centralized ordering.
** It returns:
** - 0   : if ordering data are consistent.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (1, 1)

int
dorderGather (
const Dorder * restrict const dordptr,
Order * restrict const        cordptr)
{
  Gnum                        leaflocnbr;
  int                         leafrcvnbr;
  DorderGatherLeaf * restrict leafrcvtab;
  int                         leafsndnbr;         /* "int" since used as count in MPI_Gatherv */
  DorderGatherLeaf * restrict leafsndtab;
  Gnum * restrict             perircvtab;
  int                         perisndnbr;         /* "int" since used as count in MPI_Gatherv */
  Gnum * restrict             perisndtab;
  int * restrict              recvcnttab;
  int * restrict              recvdsptab;
  const DorderLink * restrict linklocptr;
  Gnum                        vnodlocnbr;
  int                         procglbnbr;
  int                         protnum;
  Gnum                        reduloctab[2];
  Gnum                        reduglbtab[2];
  int                         cheklocval;
  int                         chekglbval;

#ifdef SCOTCH_DEBUG_DORDER2
  if ((DORDERCBLKNEDI == 0) || (DORDERCBLKNEDI != ORDERCBLKNEDI)) {
    errorPrint ("dorderGather: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DORDER2 */

  for (linklocptr = dordptr->linkdat.nextptr, leaflocnbr = vnodlocnbr = 0; /* For all nodes in local ordering structure */
       linklocptr != &dordptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST               */
    if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0) { /* If node is leaf       */
      leaflocnbr ++;                              /* One more leaf fragment     */
      vnodlocnbr += cblklocptr->data.leaf.vnodlocnbr; /* And more node vertices */
    }
#ifdef SCOTCH_DEBUG_DORDER2
    else if (cblklocptr->typeval != DORDERCBLKNEDI) {
      errorPrint ("dorderGather: invalid parameters");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
  }

  MPI_Comm_size (dordptr->proccomm, &procglbnbr);

  if (cordptr != NULL) {
    Gnum                  vnodglbnbr;

    reduloctab[0] = (Gnum) dordptr->proclocnum;
    reduloctab[1] = 1;

    vnodglbnbr = 2 * procglbnbr;                  /* TRICK: use perircvtab array as gather array        */
    if (vnodglbnbr < (dordptr->vnodglbnbr - vnodlocnbr)) /* But should receive permutation indices too! */
      vnodglbnbr = dordptr->vnodglbnbr - vnodlocnbr; /* TRICK: root will not receive from itself        */

    if (memAllocGroup ((void **) (void *)
                       &recvcnttab, (size_t) (procglbnbr * sizeof (int)),
                       &recvdsptab, (size_t) (procglbnbr * sizeof (int)),
                       &perircvtab, (size_t) (vnodglbnbr * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("dorderGather: out of memory (1)");
      reduloctab[0] = (Gnum) procglbnbr;          /* Indicate memory error */
    }
  }
  else {
    recvcnttab    = NULL;                         /* Prepare possible freeing on error */
    reduloctab[0] =
    reduloctab[1] = 0;
  }
  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 1, 1, dordptr->proccomm) != 0) {
    errorPrint ("dorderGather: communication error (1)");
    return     (1);
  }
  if (reduglbtab[1] != 1) {
    errorPrint ("dorderGather: should have only one root");
    reduglbtab[0] = (Gnum) procglbnbr;
  }
  if (reduglbtab[0] >= (Gnum) procglbnbr) {
    if (recvcnttab != NULL)
      memFree (recvcnttab);
    return (1);
  }
  protnum = (int) reduglbtab[0];

  reduloctab[0] = leaflocnbr;
  reduloctab[1] = vnodlocnbr;
  if (MPI_Gather (reduloctab, 2, GNUM_MPI, perircvtab, 2, GNUM_MPI, protnum, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGather: communication error (2)");
    return     (1);
  }

  if (dordptr->proclocnum == protnum) {
    int                   procnum;

    perircvtab[2 * protnum] = 0;                  /* TRICK: root will not send to nor receive from itself to avoid unnecessary memory copy */
    for (procnum = 0, leafrcvnbr = 0; procnum < procglbnbr; procnum ++) {
      recvdsptab[procnum] = leafrcvnbr;
      recvcnttab[procnum] = (int) (perircvtab[2 * procnum] * 2); /* TRICK: DorderGatherLeaf structures are made of 2 GNUM_MPI fields */
      leafrcvnbr         += recvcnttab[procnum];
    }
    leafrcvnbr /= 2;                              /* TRICK: restore real number of leaf structures to be received */
    leafsndnbr  = 0;
    perisndnbr  = 0;
  }
  else {
    leafrcvnbr = 0;
    leafsndnbr = (int) leaflocnbr;
    perisndnbr = (int) vnodlocnbr;
  }

  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)
                     &leafrcvtab, (size_t) (leafrcvnbr * sizeof (DorderGatherLeaf)),
                     &leafsndtab, (size_t) (leafsndnbr * sizeof (DorderGatherLeaf)),
                     &perisndtab, (size_t) (perisndnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dorderGather: out of memory (2)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGather: communication error (3)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DORDER1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DORDER1 */
  if (chekglbval != 0) {
    if (recvcnttab != NULL)
      memFree (recvcnttab);
    return (1);
  }

  if (dordptr->proclocnum == protnum) {           /* If root process */
#ifdef SCOTCH_DEBUG_DORDER2
    memSet (cordptr->peritab, ~0, dordptr->vnodglbnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_DORDER2 */

    for (linklocptr = dordptr->linkdat.nextptr; linklocptr != &dordptr->linkdat; linklocptr = linklocptr->nextptr) { /* For all nodes */
      const DorderCblk * restrict cblklocptr;

      cblklocptr = (DorderCblk *) linklocptr;     /* TRICK: FIRST                              */
      if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0)  /* If tree node is leaf, copy fragment */
        memCpy (cordptr->peritab + cblklocptr->data.leaf.ordelocval, cblklocptr->data.leaf.periloctab, cblklocptr->data.leaf.vnodlocnbr * sizeof (Gnum));
    }
  }
  else {
    Gnum                  leaflocnum;
    Gnum                  vnodlocnum;

    for (linklocptr = dordptr->linkdat.nextptr, leaflocnum = vnodlocnum = 0;
         linklocptr != &dordptr->linkdat; linklocptr = linklocptr->nextptr) { /* For all nodes */
      const DorderCblk * restrict cblklocptr;

      cblklocptr = (DorderCblk *) linklocptr;     /* TRICK: FIRST           */
      if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0) { /* If node is leaf */
        leafsndtab[leaflocnum].ordelocval = cblklocptr->data.leaf.ordelocval; /* Fill send structures with permutation data */
        leafsndtab[leaflocnum].vnodlocnbr = cblklocptr->data.leaf.vnodlocnbr;
        memCpy (perisndtab + vnodlocnum, cblklocptr->data.leaf.periloctab, cblklocptr->data.leaf.vnodlocnbr * sizeof (Gnum));
        vnodlocnum += cblklocptr->data.leaf.vnodlocnbr;
        leaflocnum ++;
      }
    }
    leafsndnbr *= 2;                              /* TRICK: DorderGatherLeaf structures are made of 2 GNUM_MPI fields */
  }

  if (MPI_Gatherv (leafsndtab, leafsndnbr, GNUM_MPI, leafrcvtab, recvcnttab, recvdsptab, GNUM_MPI, protnum, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGather: communication error (4)");
    return     (1);
  }

  if (dordptr->proclocnum == protnum) {
    int                   vnodglbnbr;
    int                   procnum;

    perircvtab[2 * protnum + 1] = 0;              /* TRICK: root will not send to nor receive from itself to avoid unnecessary memory copy */
    for (procnum = 0, vnodglbnbr = 0; procnum < procglbnbr; procnum ++) {
      recvdsptab[procnum] = vnodglbnbr;
      recvcnttab[procnum] = (int) perircvtab[2 * procnum + 1];
      vnodglbnbr         += recvcnttab[procnum];
    }
#ifdef SCOTCH_DEBUG_DORDER2
    if (((Gnum) vnodglbnbr + vnodlocnbr) != dordptr->vnodglbnbr) {
      errorPrint ("dorderGather: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
  }

  if (MPI_Gatherv (perisndtab, perisndnbr, GNUM_MPI, perircvtab, recvcnttab, recvdsptab, GNUM_MPI, protnum, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGather: communication error (5)");
    return     (1);
  }

  if (dordptr->proclocnum == protnum) {           /* If root process */
    int                   leafglbnum;
    int                   vnodglbnum;

    for (leafglbnum = vnodglbnum = 0; leafglbnum < leafrcvnbr; leafglbnum ++) {
      memCpy (cordptr->peritab + leafrcvtab[leafglbnum].ordelocval, perircvtab + vnodglbnum, leafrcvtab[leafglbnum].vnodlocnbr * sizeof (Gnum));
      vnodglbnum += leafrcvtab[leafglbnum].vnodlocnbr;
    }

    memFree (recvcnttab);                         /* Free group leader */
  }
  memFree (leafrcvtab);                           /* Free group leader */

  if (dorderGatherTree (dordptr, cordptr, protnum) != 0) /* Gather ordering tree */
    return (1);

#ifdef SCOTCH_DEBUG_DORDER2
  if (dordptr->proclocnum == protnum) {
    if (orderCheck (cordptr) != 0) {
      errorPrint ("dorderGather: invalid centralized ordering");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DORDER2 */

  return (0);
}

/* This function gathers the pieces of
** a distributed ordering tree to build a
** centralized ordering tree.
** It returns:
** - 0   : if ordering data are consistent.
** - !0  : on error.
*/

int
dorderGatherTree (
const Dorder * restrict const dordptr,
Order * restrict const        cordptr,
const int                     protnum)
{
  int                         treelocnbr;         /* "int" since used as way to fill count array in MPI_Allgather */
  Gnum                        treeglbnbr;
  DorderGatherNode * restrict treercvtab;
  int                         treesndnbr;         /* "int" since used as count in MPI_Gatherv */
  DorderGatherNode *          treesndtab;
  DorderGatherNode * restrict treesndptr;
  int * restrict              treecnttab;
  int * restrict              treedsptab;
  DorderGatherCblk * restrict cblkglbtab;
  const DorderLink * restrict linklocptr;
  int                         procglbnbr;
  int                         procnum;
  int                         cheklocval;
  int                         chekglbval;

  for (linklocptr = dordptr->linkdat.nextptr, treelocnbr = 0; /* Count only purely local nodes */
       linklocptr != &dordptr->linkdat; linklocptr = linklocptr->nextptr) {
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST */
#ifdef SCOTCH_DEBUG_DORDER2
    if ((cblklocptr->cblknum.proclocnum != dordptr->proclocnum) && /* Local sub-nodes of non-locally rooted node not implemented */
        ((cblklocptr->typeval & DORDERCBLKLEAF) != 0)           &&
        (cblklocptr->data.leaf.nodelocnbr != 0)) {
      errorPrint ("dorderGatherTree: not implemented");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */

    if (cblklocptr->cblknum.proclocnum == dordptr->proclocnum) {
      treelocnbr ++;
      if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0)
        treelocnbr += (int) cblklocptr->data.leaf.nodelocnbr;
    }
  }

  MPI_Comm_size (dordptr->proccomm, &procglbnbr);

  treesndnbr = (dordptr->proclocnum == protnum) ? 0 : treelocnbr; /* TRICK: root will not send nor receive */
  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)
                     &treecnttab, (size_t) (procglbnbr * sizeof (int)),
                     &treedsptab, (size_t) (procglbnbr * sizeof (int)),
                     &treesndtab, (size_t) (treesndnbr * sizeof (DorderGatherNode)), NULL) == NULL) {
    errorPrint ("dorderGatherTree: out of memory (1)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGatherTree: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DORDER1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DORDER1 */
  if (chekglbval != 0) {
    if (treecnttab != NULL)
      memFree (treecnttab);
    return (1);
  }

  if (MPI_Allgather (&treelocnbr, 1, MPI_INT, treecnttab, 1, MPI_INT, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGatherTree: communication error (2)");
    return     (1);
  }

  for (procnum = 0, treeglbnbr = 0; procnum < procglbnbr; procnum ++) { /* Compute prefix sum of local numbers for global numbering */
    treedsptab[procnum] = treeglbnbr;
    treeglbnbr         += treecnttab[procnum];
  }
  if (dordptr->proclocnum == protnum) {
    treecnttab[protnum] = 0;                      /* TRICK: root will not send to nor receive from itself to avoid unnecessary memory copy */

    cordptr->treenbr = treeglbnbr;

    if (memAllocGroup ((void **) (void *)
                       &treercvtab, (size_t) (treeglbnbr * sizeof (DorderGatherNode)),
                       &cblkglbtab, (size_t) (treeglbnbr * sizeof (DorderGatherCblk)), NULL) == NULL) {
      errorPrint ("dorderGatherTree: out of memory (2)");
      cheklocval = 1;
    }
    treesndptr = treercvtab + treedsptab[protnum]; /* TRICK: root process will build its column blocks in place as if received */
  }
  else
    treesndptr = treesndtab;

#ifdef SCOTCH_DEBUG_DORDER1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGather: communication error (3)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DORDER1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DORDER1 */
  if (chekglbval != 0) {
    memFree (treecnttab);
    return  (1);
  }

  for (linklocptr = dordptr->linkdat.nextptr; linklocptr != &dordptr->linkdat; linklocptr = linklocptr->nextptr) { /* For all nodes */
    const DorderCblk * restrict cblklocptr;

    cblklocptr = (DorderCblk *) linklocptr;       /* TRICK: FIRST                      */
    if (cblklocptr->cblknum.proclocnum != dordptr->proclocnum) /* Skip non-local nodes */
      continue;

    treesndptr->fathnum = treedsptab[cblklocptr->fathnum.proclocnum] + cblklocptr->fathnum.cblklocnum; /* If node is (part of) the root node */
    treesndptr->typeval = (Gnum) (((cblklocptr->typeval & DORDERCBLKNEDI) != 0) ? ORDERCBLKNEDI : ORDERCBLKOTHR);
    treesndptr->vnodnbr = cblklocptr->vnodglbnbr;
    treesndptr->cblknum = cblklocptr->cblkfthnum;
    treesndptr ++;

    if ((cblklocptr->typeval & DORDERCBLKLEAF) != 0) {  /* If node is a distributed leaf */
      Gnum                        cblkglbnum;
      Gnum                        cblkglbadj;
      const DorderNode * restrict nodelocptr;
      const DorderNode * restrict nodeloctnd;

      cblkglbnum = treedsptab[cblklocptr->cblknum.proclocnum] + cblklocptr->cblknum.cblklocnum;
      cblkglbadj = treedsptab[cblklocptr->cblknum.proclocnum] + cblklocptr->data.leaf.cblklocnum;

      for (nodelocptr = cblklocptr->data.leaf.nodeloctab, nodeloctnd = nodelocptr + cblklocptr->data.leaf.nodelocnbr;
           nodelocptr < nodeloctnd; nodelocptr ++) { /* Build nodes for all local nodes */
        treesndptr->fathnum = (nodelocptr->fathnum == -1) ? cblkglbnum : (nodelocptr->fathnum + cblkglbadj);
        treesndptr->typeval = (Gnum) nodelocptr->typeval;
        treesndptr->vnodnbr = nodelocptr->vnodnbr;
        treesndptr->cblknum = nodelocptr->cblknum;
        treesndptr ++;
      }
    }
#ifdef SCOTCH_DEBUG_DORDER2
    else if (cblklocptr->typeval != DORDERCBLKNEDI) {
      errorPrint ("dorderGatherTree: invalid column block type");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DORDER2 */
  }
#ifdef SCOTCH_DEBUG_DORDER2
  if (treesndptr != ((dordptr->proclocnum == protnum) ? (treercvtab + treedsptab[protnum]) : treesndtab) + treelocnbr) {
    errorPrint ("dorderGatherTree: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DORDER2 */

  if (dordptr->proclocnum == protnum) {           /* If node is root, adjust displacements in terms of Gnum and not DorderGatherNode */
    for (procnum = 0; procnum < procglbnbr; procnum ++) {
      treecnttab[procnum] *= DORDERGATHERNODESIZE;
      treedsptab[procnum] *= DORDERGATHERNODESIZE;
    }
  }
  if (MPI_Gatherv (treesndtab, treesndnbr * DORDERGATHERNODESIZE, GNUM_MPI,
                   treercvtab, treecnttab, treedsptab, GNUM_MPI, protnum, dordptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderGatherTree: communication error (4)");
    return     (1);
  }

  if (dordptr->proclocnum == protnum) {
    Gnum                  treeglbnum;
    Gnum                  cblkglbnbr;

    memSet (cblkglbtab, 0, treeglbnbr * sizeof (DorderGatherCblk)); /* Set all son counters to zero and all array pointers to NULL */

    for (treeglbnum = 1; treeglbnum < treeglbnbr; treeglbnum ++) { /* For all local and received tree nodes except root node */
      Gnum                  cblkfthnum;

      cblkfthnum = treercvtab[treeglbnum].fathnum;
#ifdef SCOTCH_DEBUG_DORDER2
      if ((cblkfthnum < 0) ||                     /* Father of non-root node cannot be -1                 */
          (cblkfthnum >= treeglbnum)) {           /* Father should always have smaller global node number */
        errorPrint ("dorderGatherTree: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DORDER2 */
      cblkglbtab[cblkfthnum].cblknbr ++;          /* Add a son to its father */
    }

    for (treeglbnum = 0, cblkglbnbr = treeglbnbr; treeglbnum < treeglbnbr; treeglbnum ++) { /* For all local and received tree nodes */
      if (cblkglbtab[treeglbnum].cblknbr > 0) {
#ifdef SCOTCH_DEBUG_DORDER2
        if (cblkglbtab[treeglbnum].cblknbr < 2) { /* Descendent nodes should comprise at least two column block slots */
          errorPrint ("dorderGatherTree: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DORDER2 */
        cblkglbnbr --;                            /* One new subblock means one more shared frontier, so one less column block than nodes */
        if ((cblkglbtab[treeglbnum].cblktab = memAlloc (cblkglbtab[treeglbnum].cblknbr * sizeof (OrderCblk))) == NULL) {
          errorPrint ("dorderGather: out of memory (3)");
          while (-- treeglbnum >= 0) {
            if (cblkglbtab[treeglbnum].cblktab != NULL)
              memFree (cblkglbtab[treeglbnum].cblktab);
          }
          memFree (treercvtab);
          memFree (treecnttab);
          return (1);
        }
      }
    }
    cordptr->cblknbr = cblkglbnbr;

    cordptr->cblktre.typeval = (int) treercvtab[0].typeval; /* Process root node of separator tree */
    cordptr->cblktre.vnodnbr = treercvtab[0].vnodnbr;
    cordptr->cblktre.cblknbr = cblkglbtab[0].cblknbr;
    cordptr->cblktre.cblktab = cblkglbtab[0].cblktab; /* Link its sons array  */

    for (treeglbnum = 1; treeglbnum < treeglbnbr; treeglbnum ++) { /* For all nodes except the root */
      Gnum                  cblkfthnum;
      OrderCblk * restrict  cblksonptr;

      cblkfthnum = treercvtab[treeglbnum].cblknum;
#ifdef SCOTCH_DEBUG_DORDER2
      if ((cblkfthnum < 0) || (cblkfthnum >= cblkglbtab[treercvtab[treeglbnum].fathnum].cblknbr)) {
        errorPrint ("dorderGatherTree: internal error (4)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DORDER2 */
      cblksonptr = &cblkglbtab[treercvtab[treeglbnum].fathnum].cblktab[cblkfthnum]; /* Point to son's slot in father array */
      cblksonptr->typeval = (int) treercvtab[treeglbnum].typeval;
      cblksonptr->vnodnbr = treercvtab[treeglbnum].vnodnbr;
      cblksonptr->cblknbr = cblkglbtab[treeglbnum].cblknbr; /* Link son column block array to column block structure */
      cblksonptr->cblktab = cblkglbtab[treeglbnum].cblktab;
    }

    memFree (treercvtab);
  }
  memFree (treecnttab);

  return (0);
}
