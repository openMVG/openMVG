/* Copyright 2007-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_fold.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the halo distribu-  **/
/**                ted graph folding function.             **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 27 jun 2008     **/
/**                                 to   : 04 jan 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HDGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_fold_comm.h"
#include "hdgraph.h"
#include "hdgraph_fold.h"

/******************************/
/*                            */
/* These routines handle halo */
/* distributed source graphs. */
/*                            */
/******************************/

/* This routine builds a folded graph by
** merging graph data to the processes of
** the first half or to the second half
** of the communicator.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hdgraphFold (
const Hdgraph * restrict const  orggrafptr,
const int                       partval,          /* 0 for first half, 1 for second half */
Hdgraph * restrict const        fldgrafptr)
{
  int               fldprocglbnbr;
  int               fldproclocnum;                /* Index of local process in folded communicator   */
  int               fldproccol;                   /* Color of receiver or not wanted in communicator */
  MPI_Comm          fldproccomm;                  /* Communicator of folded part                     */

  fldprocglbnbr = (orggrafptr->s.procglbnbr + 1) / 2;
  if (partval == 1) {
    fldproclocnum = orggrafptr->s.proclocnum - fldprocglbnbr;
    fldprocglbnbr = orggrafptr->s.procglbnbr - fldprocglbnbr; 
  }
  else
    fldproclocnum = orggrafptr->s.proclocnum;

  fldproccol = ((fldproclocnum >= 0) && (fldproclocnum < fldprocglbnbr)) ? 0 : MPI_UNDEFINED;

  if (MPI_Comm_split (orggrafptr->s.proccomm, fldproccol, fldproclocnum, &fldproccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphFold: communication error");
    return     (1);
  }

  return (hdgraphFold2 (orggrafptr, partval, fldgrafptr, fldproccomm));
}

int
hdgraphFold2 (
const Hdgraph * restrict const  orggrafptr,
const int                       partval,          /* 0 for first half, 1 for second half */
Hdgraph * restrict const        fldgrafptr,
MPI_Comm                        fldproccomm)      /* Pre-computed communicator */
{
  int                           fldcommtypval;    /* Type of communication for this process              */
  DgraphFoldCommData * restrict fldcommdattab;    /* Array of two communication data                     */
  Gnum * restrict               fldcommvrttab;    /* Starting global send indices of communications      */
  Gnum * restrict               fldvertidxtab;    /* Start indices of vertex arrays                      */
  Gnum * restrict               fldvendidxtab;    /* Adjustment value for end vertex arrays              */
  Gnum * restrict               fldedgeidxtab;    /* Start indices of edge arrays                        */
  Gnum * restrict               fldedgecnttab;    /* Number of edges exchanged during each communication */
  Gnum                          fldvertlocnbr;    /* Number of vertices in local folded part             */
  Gnum                          fldedgelocsiz;    /* (Upper bound of) number of edges in folded graph    */
  int                           fldprocglbnbr;
  int                           fldproclocnum;    /* Index of local process in folded communicator       */
  int                           fldvertadjnbr;
  Gnum * restrict               fldvertadjtab;    /* Array of index adjustments for original vertices    */
  Gnum * restrict               fldvertdlttab;    /* Array of index adjustments for original vertices    */
  Gnum * restrict               fldvhalloctax;    /* Index array for remote halo vertex renumbering      */
  int                           cheklocval;
  int                           chekglbval;
  int                           commmax;
  int                           commnbr;
  int                           requnbr;
  MPI_Request * restrict        requtab;

#ifdef SCOTCH_DEBUG_HDGRAPH2
  if (orggrafptr->vhndloctax != (orggrafptr->s.vertloctax + 1)) {
    errorPrint ("hdgraphFold2: halo graph must be compact");
    return     (1);
  }
  if (orggrafptr->s.vendloctax < (orggrafptr->s.vertloctax + orggrafptr->s.vertlocnbr)) { /* MPI_Isend calls should not overlap */
    errorPrint ("hdgraphFold2: halo graph must have distinct arrays");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  fldprocglbnbr = (orggrafptr->s.procglbnbr + 1) / 2;
  if (partval == 1) {
    fldproclocnum = orggrafptr->s.proclocnum - fldprocglbnbr;
    fldprocglbnbr = orggrafptr->s.procglbnbr - fldprocglbnbr; 
  }
  else
    fldproclocnum = orggrafptr->s.proclocnum;

  fldcommtypval = ((fldproclocnum >= 0) && (fldproclocnum < fldprocglbnbr)) ? DGRAPHFOLDCOMMRECV : DGRAPHFOLDCOMMSEND;

  cheklocval    = 0;
  fldvertidxtab = NULL;
  fldcommdattab = NULL;
  if (fldcommtypval == DGRAPHFOLDCOMMRECV) {      /* If we are going to receive */
#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (fldgrafptr == NULL) {
      errorPrint ("hdgraphFold2: invalid parameters (1)");
      return     (1);
    }
    if (fldproccomm == MPI_COMM_NULL) {
      errorPrint ("hdgraphFold2: invalid parameters (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

    memSet (fldgrafptr, 0, sizeof (Hdgraph));     /* Pre-initialize graph fields */

    fldgrafptr->s.proccomm   = fldproccomm;
    fldgrafptr->s.procglbnbr = fldprocglbnbr;
    fldgrafptr->s.proclocnum = fldproclocnum;
    fldgrafptr->s.flagval    = DGRAPHFREEALL | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* For premature freeing on error; do not free vhndloctab as it is grouped with vertloctab */

    if (memAllocGroup ((void **) (void *)         /* Allocate distributed graph private data */
                       &fldgrafptr->s.procdsptab, (size_t) ((fldprocglbnbr + 1) * sizeof (Gnum)),
                       &fldgrafptr->s.proccnttab, (size_t) (fldprocglbnbr       * sizeof (Gnum)),
                       &fldgrafptr->s.procngbtab, (size_t) (fldprocglbnbr       * sizeof (int)),
                       &fldgrafptr->s.procrcvtab, (size_t) (fldprocglbnbr       * sizeof (int)),
                       &fldgrafptr->s.procsndtab, (size_t) (fldprocglbnbr       * sizeof (int)), NULL) == NULL) {
      errorPrint ("hdgraphFold2: out of memory (1)");
      cheklocval = 1;
    }
    else if (dgraphFoldComm (&orggrafptr->s, partval, &commmax, &fldcommtypval, &fldcommdattab, &fldcommvrttab, /* Process can become a sender receiver */
                             fldgrafptr->s.proccnttab, &fldvertadjnbr, &fldvertadjtab, &fldvertdlttab) != 0) {
      errorPrint ("hdgraphFold2: cannot compute folding communications (1)");
      cheklocval = 1;
    }
    else {
      Gnum              fldvelolocnbr;

      if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver */
        int               i;

        for (i = 0, fldvertlocnbr = orggrafptr->s.vertlocnbr; (i < commmax) && (fldcommdattab[i].procnum != -1); i ++)
          fldvertlocnbr += fldcommdattab[i].vertnbr;
        commnbr = i;

        fldedgelocsiz = orggrafptr->s.edgelocsiz + orggrafptr->s.edgeglbsmx * i; /* Upper bound on local edges (degree useless since only for non-halo vertices) */
      }
      else {                                      /* Process is a sender receiver */
        fldvertlocnbr = fldcommvrttab[0] - orggrafptr->s.procvrttab[orggrafptr->s.proclocnum]; /* Communications will remove vertices     */
        fldedgelocsiz = orggrafptr->s.vertloctax[fldvertlocnbr + orggrafptr->s.baseval] - orggrafptr->s.baseval; /* Exact number of edges */

        fldgrafptr->s.edgelocsiz = fldedgelocsiz;
      }
      fldvelolocnbr = (orggrafptr->s.veloloctax != NULL) ? fldvertlocnbr : 0;

      if (memAllocGroup ((void **) (void *)       /* Allocate distributed graph public data */
                         &fldgrafptr->s.vertloctax, (size_t) ((fldvertlocnbr + 1) * sizeof (Gnum)),
                         &fldgrafptr->s.vendloctax, (size_t) ( fldvertlocnbr      * sizeof (Gnum)), /* Vertex end array for non-halo vertices */
                         &fldgrafptr->s.vnumloctax, (size_t) ( fldvertlocnbr      * sizeof (Gnum)),
                         &fldgrafptr->s.veloloctax, (size_t) ( fldvelolocnbr      * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("hdgraphFold2: out of memory (2)");
        cheklocval = 1;
      }
      else if (fldgrafptr->s.vertloctax -= orggrafptr->s.baseval,
               fldgrafptr->s.vendloctax -= orggrafptr->s.baseval,
               fldgrafptr->s.vnumloctax -= orggrafptr->s.baseval,
               fldgrafptr->s.veloloctax = ((orggrafptr->s.veloloctax != NULL) ? fldgrafptr->s.veloloctax - orggrafptr->s.baseval : NULL),
               memAllocGroup ((void **) (void *)
                              &fldgrafptr->s.edgeloctax, (size_t) (fldedgelocsiz            * sizeof (Gnum)),
                              &fldvhalloctax,            (size_t) (orggrafptr->s.edgeglbsmx * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("hdgraphFold2: out of memory (3)");
        cheklocval = 1;
      }
      else {
        fldgrafptr->s.edgeloctax -= orggrafptr->s.baseval;
        fldvhalloctax            -= orggrafptr->s.baseval;
      }
    }
  }
  else {                                          /* Process is a sender */
#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (fldproccomm != MPI_COMM_NULL) {
      errorPrint ("hdgraphFold2: invalid parameters (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

    if (dgraphFoldComm (&orggrafptr->s, partval, &commmax, &fldcommtypval, &fldcommdattab, &fldcommvrttab, NULL, NULL, NULL, NULL) != 0) {
      errorPrint ("hdgraphFold2: cannot compute folding communications (2)");
      cheklocval = 1;
    }
  }

  if ((cheklocval == 0) &&
      (memAllocGroup ((void **) (void *)         /* Allocate folding data */
                      &fldvertidxtab, (size_t) (commmax * sizeof (Gnum)),
                      &fldvendidxtab, (size_t) (commmax * sizeof (Gnum)),
                      &fldedgeidxtab, (size_t) (commmax * sizeof (Gnum)),
                      &fldedgecnttab, (size_t) (commmax * sizeof (Gnum)),
                      &requtab,       (size_t) (commmax * HDGRAPHFOLDTAGNBR * sizeof (MPI_Request)), NULL) == NULL)) {
    errorPrint ("hdgraphFold2: out of memory (4)");
    cheklocval = 1;
  }

#ifdef SCOTCH_DEBUG_HDGRAPH1                      /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, orggrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphFold2: communication error (1)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_HDGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_HDGRAPH1 */
  if (chekglbval != 0) {
    if ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) {
      hdgraphExit (fldgrafptr);
      if (fldvertidxtab != NULL)
        memFree (fldvertidxtab);                  /* Free group leaders */
      if (fldcommdattab != NULL)
        memFree (fldcommdattab);
    }
    return (1);
  }

  requnbr = 0;                                    /* Communications without further processing are placed at beginning of array */

  if ((fldcommtypval & DGRAPHFOLDCOMMSEND) != 0) { /* If process is (also) a sender */
    Gnum              vertsndbas;
    Gnum              vertsndnbr;
    int               i;

    vertsndnbr = ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) ? (fldcommvrttab[0] - orggrafptr->s.procvrttab[orggrafptr->s.proclocnum]) : 0; /* If process is also a receiver, start sending after kept vertices */

    for (i = 0, requnbr = 0, vertsndbas = orggrafptr->s.baseval; /* For all send communications to perform */
         (i < commmax) && (fldcommdattab[i].procnum != -1); i ++) {
      Gnum              edgelocsiz;

      vertsndbas += vertsndnbr;
      vertsndnbr  = fldcommdattab[i].vertnbr;
      edgelocsiz  = orggrafptr->s.vertloctax[vertsndbas + vertsndnbr] - orggrafptr->s.vertloctax[vertsndbas]; /* Graph is compact */

      fldvertidxtab[i] = vertsndbas;
      fldedgeidxtab[i] = orggrafptr->s.vertloctax[vertsndbas];
      fldedgecnttab[i] = edgelocsiz;
      if (MPI_Isend (&edgelocsiz, 1, GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGVLBLLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (2)");
        cheklocval = 1;
      }
    }
    commnbr = i;

    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      if (MPI_Isend (orggrafptr->s.vertloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGVERTLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (3)");
        cheklocval = 1;
      }
    }
    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      if (MPI_Isend (orggrafptr->s.vendloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGVENDLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (4)");
        cheklocval = 1;
      }
    }
    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      if (MPI_Isend (orggrafptr->s.edgeloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGEDGELOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (5)");
        cheklocval = 1;
      }
    }
    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      Gnum              vertsndbas;
      int               vertsndnbr;
      int               procsndnum;               /* Rank of process to send to */

      vertsndbas = fldvertidxtab[i];
      vertsndnbr = (int) fldcommdattab[i].vertnbr;
      procsndnum = (int) fldcommdattab[i].procnum;
      if ((orggrafptr->s.veloloctax != NULL) &&
          (MPI_Isend (orggrafptr->s.veloloctax + vertsndbas, vertsndnbr, GNUM_MPI, procsndnum,
                      TAGFOLD + TAGVELOLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
        errorPrint ("hdgraphFold2: communication error (6)");
        cheklocval = 1;
      }
      else if ((orggrafptr->s.vnumloctax != NULL) &&
               (MPI_Isend (orggrafptr->s.vnumloctax + vertsndbas, vertsndnbr, GNUM_MPI, procsndnum,
                           TAGFOLD + TAGVNUMLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
        errorPrint ("hdgraphFold2: communication error (7)");
        cheklocval = 1;
      }
    }                                             /* Communications of sender-receivers will be completed in the receiving phase */
  }

  if ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) { /* If process is (also) a receiver */
    Gnum                orgvertlocnbr;
    Gnum                orgvertlocnnd;
    Gnum                fldvertlocadj;
    Gnum                fldvelolocsum;
    Gnum                fldedgelocnum;
    Gnum                fldvhallocnum;
    Gnum                fldehallocnbr;
    int                 fldprocnum;
    int                 i;

    const Gnum * restrict const orgvertloctax = orggrafptr->s.vertloctax;
    const Gnum * restrict const orgvendloctax = orggrafptr->s.vendloctax;
    const Gnum * restrict const orgedgeloctax = orggrafptr->s.edgeloctax;

    fldgrafptr->s.procvrttab = fldgrafptr->s.procdsptab; /* Graph does not have holes                           */
    fldgrafptr->s.procdsptab[0] = orggrafptr->s.baseval; /* Build private data of folded graph and array        */
    for (fldprocnum = 0; fldprocnum < fldprocglbnbr; fldprocnum ++) /* New subdomain indices start from baseval */
      fldgrafptr->s.procdsptab[fldprocnum + 1] = fldgrafptr->s.procdsptab[fldprocnum] + fldgrafptr->s.proccnttab[fldprocnum];

    if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver */
      Gnum                orgvertlocmin;
      Gnum                orgvertlocmax;
      Gnum                fldvertlocnum;
      Gnum                fldedgelocbas;
      Gnum                fldvertrcvbas;
      Gnum                fldvertrcvnbr;
      int                 procngbmin;
      int                 procngbmax;

      Gnum * restrict const fldedgeloctax = fldgrafptr->s.edgeloctax;

      for (i = 0, fldvertrcvbas = orggrafptr->s.vertlocnnd, fldvertrcvnbr = 0; /* For all receive communications to perform */
           (i < commnbr) && (cheklocval == 0); i ++) {
        fldvertrcvbas += fldvertrcvnbr;
        fldvertrcvnbr  = fldcommdattab[i].vertnbr;

        fldvertidxtab[i] = fldvertrcvbas;
        if (MPI_Irecv (&fldedgecnttab[i], 1, GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGVLBLLOCTAB, orggrafptr->s.proccomm, &requtab[HDGRAPHFOLDTAGENBR * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("hdgraphFold2: communication error (8)");
          cheklocval = 1;
        }
      }

      for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) { /* Let these communications progress while we process the edge size messages */
        if (MPI_Irecv (fldgrafptr->s.vertloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGVERTLOCTAB, orggrafptr->s.proccomm, &requtab[HDGRAPHFOLDTAGVERT * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("hdgraphFold2: communication error (9)");
          cheklocval = 1;
        }
      }
      for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
        if (MPI_Irecv (fldgrafptr->s.vendloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGVENDLOCTAB, orggrafptr->s.proccomm, &requtab[HDGRAPHFOLDTAGVEND * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("hdgraphFold2: communication error (10)");
          cheklocval = 1;
        }
      }

      MPI_Waitall (commnbr, &requtab[HDGRAPHFOLDTAGENBR * commmax], MPI_STATUSES_IGNORE);

      for (i = 0, fldedgelocbas = orggrafptr->s.vertloctax[orggrafptr->s.vertlocnnd]; (i < commnbr) && (cheklocval == 0); i ++) {
        fldedgeidxtab[i] = fldedgelocbas;
        fldedgelocbas += fldedgecnttab[i];

        if (MPI_Irecv (fldedgeloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGEDGELOCTAB, orggrafptr->s.proccomm, &requtab[HDGRAPHFOLDTAGEDGE * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("hdgraphFold2: communication error (11)");
          cheklocval = 1;
        }
      }
      fldgrafptr->s.edgelocsiz = fldedgelocbas - orggrafptr->s.baseval; /* Get number of local and halo edges */

      if (orggrafptr->s.veloloctax != NULL) {
        for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
          if (MPI_Irecv (fldgrafptr->s.veloloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                         TAGFOLD + TAGVELOLOCTAB, orggrafptr->s.proccomm, &requtab[HDGRAPHFOLDTAGVELO * commmax + i]) != MPI_SUCCESS) {
            errorPrint ("hdgraphFold2: communication error (12)");
            cheklocval = 1;
          }
        }
      }
      if (orggrafptr->s.vnumloctax != NULL) {
        for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
          if (MPI_Irecv (fldgrafptr->s.vnumloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                         TAGFOLD + TAGVNUMLOCTAB, orggrafptr->s.proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
            errorPrint ("hdgraphFold2: communication error (13)");
            cheklocval = 1;
          }
        }
      }

      orgvertlocnbr = orggrafptr->s.vertlocnbr;   /* Process all local vertices */
      orgvertlocnnd = orggrafptr->s.vertlocnnd;

      if (orggrafptr->s.vnumloctax == NULL) {     /* If original graph does not have vertex numbers, create remote parts of vertex number array */
        Gnum              fldvertlocnum;
        Gnum              fldvertlocadj;
        int               i;

        Gnum * restrict const fldvnumloctax = fldgrafptr->s.vnumloctax;

        for (i = 0, fldvertlocnum = orgvertlocnnd; i < commnbr; i ++) {
          Gnum              fldvertlocnnd;

          for (fldvertlocnnd = fldvertlocnum + fldcommdattab[i].vertnbr, fldvertlocadj = fldcommvrttab[i];
               fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
            fldvnumloctax[fldvertlocnum] = fldvertlocadj ++;
        }
      }

      for (procngbmin = 0, procngbmax = fldvertadjnbr; /* Initialize search accelerator */
           procngbmax - procngbmin > 1; ) {
        int               procngbmed;

        procngbmed = (procngbmax + procngbmin) / 2;
        if (fldvertadjtab[procngbmed] <= orggrafptr->s.procvrttab[orggrafptr->s.proclocnum])
          procngbmin = procngbmed;
        else
          procngbmax = procngbmed;
      }
      orgvertlocmin = fldvertadjtab[procngbmin];
      orgvertlocmax = fldvertadjtab[procngbmax];
      fldvertlocadj = fldvertdlttab[procngbmin];
      for (fldvertlocnum = fldedgelocnum = orggrafptr->s.baseval; /* Adjust local part of edge array */
           fldvertlocnum < orgvertlocnnd; ) {
        for ( ; fldedgelocnum < orgvendloctax[fldvertlocnum]; fldedgelocnum ++) { /* Reorder end vertices */
          Gnum              orgvertlocend;

#ifdef SCOTCH_DEBUG_HDGRAPH2
          if (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval)) {
            errorPrint  ("hdgraphFold2: internal error (1)");
            return      (1);
          }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

          orgvertlocend = orgedgeloctax[fldedgelocnum];

          if ((orgvertlocend >= orgvertlocmin) && /* If end vertex is local */
              (orgvertlocend <  orgvertlocmax))
            fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertlocadj;
          else {                                  /* End vertex is not local */
            int               procngbmin;
            int               procngbmax;

            for (procngbmin = 0, procngbmax = fldvertadjnbr;
                 procngbmax - procngbmin > 1; ) {
              int               procngbnum;

              procngbnum = (procngbmax + procngbmin) / 2;
              if (fldvertadjtab[procngbnum] <= orgvertlocend)
                procngbmin = procngbnum;
              else
                procngbmax = procngbnum;
            }
            fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertdlttab[procngbmin];
          }
        }
        fldvertlocnum ++;
        for ( ; fldedgelocnum < orgvertloctax[fldvertlocnum]; fldedgelocnum ++) { /* Copy halo part as is */
#ifdef SCOTCH_DEBUG_HDGRAPH2
          if ((orgedgeloctax[fldedgelocnum] < orggrafptr->s.baseval) ||
              (orgedgeloctax[fldedgelocnum] >= (orggrafptr->vhallocnbr + orggrafptr->s.baseval))) {
            errorPrint  ("hdgraphFold2: internal error (2)");
            return      (1);
          }
          if (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval)) {
            errorPrint  ("hdgraphFold2: internal error (3)");
            return      (1);
          }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
          fldedgeloctax[fldedgelocnum] = orgedgeloctax[fldedgelocnum];
        }
      }

      fldvelolocsum = orggrafptr->s.velolocsum;   /* In case there are vertex loads, we keep all of existing load    */
      fldehallocnbr = orggrafptr->ehallocnbr;     /* Normal receivers have at least all of their local halo vertices */
      fldvhallocnum = orggrafptr->vhallocnbr + orggrafptr->s.baseval; /* Index of next halo vertex number to assign  */
    }
    else {                                        /* Receiver process is also a sender */
      Gnum              orgvertlocmin;
      Gnum              orgvertlocmax;
      Gnum              fldvertlocnum;
      Gnum              fldvertlocadj;
      Gnum              fldvhallocmax;            /* Maximum current size of halo vertex array */
      int               procngbmin;
      int               procngbmax;

      Gnum * restrict const fldedgeloctax = fldgrafptr->s.edgeloctax;

      orgvertlocnbr = fldvertlocnbr;              /* Process only remaining local vertices */
      orgvertlocnnd = fldvertlocnbr + orggrafptr->s.baseval;

      for (procngbmin = 0, procngbmax = fldvertadjnbr; /* Initialize search accelerator */
           procngbmax - procngbmin > 1; ) {
        int               procngbmed;

        procngbmed = (procngbmax + procngbmin) / 2;
        if (fldvertadjtab[procngbmed] <= orggrafptr->s.procvrttab[orggrafptr->s.proclocnum])
          procngbmin = procngbmed;
        else
          procngbmax = procngbmed;
      }
      orgvertlocmin = fldvertadjtab[procngbmin];
      orgvertlocmax = fldvertadjtab[procngbmax];
      fldvertlocadj = fldvertdlttab[procngbmin];
      fldvhallocmax = orggrafptr->s.baseval - 1;  /* Reset halo vertex array for local part as halo vertices may have disappeared */
      fldehallocnbr = 0;                          /* Recount all remaining halo vertices and edges                                */
      fldvhallocnum = orggrafptr->s.baseval;
      for (fldvertlocnum = fldedgelocnum = orggrafptr->s.baseval; /* Copy remaining local part of edge array */
           fldvertlocnum < orgvertlocnnd; ) {
        for ( ; fldedgelocnum < orgvendloctax[fldvertlocnum]; fldedgelocnum ++) { /* Reorder end vertices */
          Gnum              orgvertlocend;

#ifdef SCOTCH_DEBUG_HDGRAPH2
          if (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval)) {
            errorPrint  ("hdgraphFold2: internal error (4)");
            return      (1);
          }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

          orgvertlocend = orgedgeloctax[fldedgelocnum];

          if ((orgvertlocend >= orgvertlocmin) && /* If end vertex is local */
              (orgvertlocend <  orgvertlocmax))
            fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertlocadj;
          else {                                  /* End vertex is not local */
            int               procngbnum;
            int               procngbmax;

            for (procngbnum = 0, procngbmax = fldvertadjnbr;
                 procngbmax - procngbnum > 1; ) {
              int               procngbmed;

              procngbmed = (procngbmax + procngbnum) / 2;
              if (fldvertadjtab[procngbmed] <= orgvertlocend)
                procngbnum = procngbmed;
              else
                procngbmax = procngbmed;
            }
            fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertdlttab[procngbnum];
          }
        }
        fldvertlocnum ++;
        fldehallocnbr += orgvertloctax[fldvertlocnum] - fldedgelocnum;
        for ( ; fldedgelocnum < orgvertloctax[fldvertlocnum]; fldedgelocnum ++) { /* Renumber halo part */
          Gnum              orgverthalend;
          Gnum              fldvhallocend;

          orgverthalend = orgedgeloctax[fldedgelocnum];
#ifdef SCOTCH_DEBUG_HDGRAPH2
          if ((orgverthalend < orggrafptr->s.baseval)                               ||
              (orgverthalend >= (orggrafptr->vhallocnbr   + orggrafptr->s.baseval)) ||
              (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval))) {
            errorPrint ("hdgraphFold2: internal error (5)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

          while (fldvhallocmax < orgverthalend)   /* Expand halo vertex index array whenever necessary */
            fldvhalloctax[++ fldvhallocmax] = ~0;
          fldvhallocend = fldvhalloctax[orgverthalend]; /* Get renumbered halo vertex */
          if (fldvhallocend < 0) {                /* If new halo vertex not yet given */
            fldvhallocend                =        /* Allocate it                      */
            fldvhalloctax[orgverthalend] = fldvhallocnum ++;
          }
          fldedgeloctax[fldedgelocnum] = fldvhallocend;
        }
      }

      if (orggrafptr->s.veloloctax != NULL) {     /* If original graph has vertex loads */
        Gnum                fldvertlocnum;

        for (fldvertlocnum = orggrafptr->s.baseval, fldvelolocsum = 0; /* Accumulate load sum of remaining part */
             fldvertlocnum < orgvertlocnnd; fldvertlocnum ++)
          fldvelolocsum += orggrafptr->s.veloloctax[fldvertlocnum];
      }

      commnbr = 0;                                /* Turn sender-receiver into normal receiver without any communications to perform */
    }

    if (orggrafptr->s.veloloctax != NULL)         /* If original graph has vertex loads                 */
      memCpy (fldgrafptr->s.veloloctax + orggrafptr->s.baseval, /* Copy local part of vertex load array */
              orggrafptr->s.veloloctax + orggrafptr->s.baseval, orgvertlocnbr * sizeof (Gnum));

    if (orggrafptr->s.vnumloctax != NULL)         /* If original graph has vertex numbers                 */
      memCpy (fldgrafptr->s.vnumloctax + orggrafptr->s.baseval, /* Copy local part of vertex number array */
              orggrafptr->s.vnumloctax + orggrafptr->s.baseval, orgvertlocnbr * sizeof (Gnum));
    else {                                        /* Build local part of vertex number array */
      Gnum              fldvertlocnum;
      Gnum              fldvertlocadj;

      for (fldvertlocnum = orggrafptr->s.baseval,
           fldvertlocadj = orggrafptr->s.procvrttab[orggrafptr->s.proclocnum];
           fldvertlocnum < orgvertlocnnd; fldvertlocnum ++)
        fldgrafptr->s.vnumloctax[fldvertlocnum] = fldvertlocadj ++;
    }

    memCpy (fldgrafptr->s.vertloctax + orggrafptr->s.baseval, /* Copy local part of vertex arrays, since they are compact */
            orggrafptr->s.vertloctax + orggrafptr->s.baseval, orgvertlocnbr * sizeof (Gnum)); /* Last value not copied    */
    fldgrafptr->s.vertloctax[fldvertlocnbr + orggrafptr->s.baseval] = fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval;
    memCpy (fldgrafptr->s.vendloctax + orggrafptr->s.baseval,
            orggrafptr->s.vendloctax + orggrafptr->s.baseval, orgvertlocnbr * sizeof (Gnum));

    for (i = 0; i < commnbr; i ++) {
      int               j;

      if (MPI_Waitany (commnbr, &requtab[HDGRAPHFOLDTAGVERT * commmax], &j, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (14)");
        cheklocval = 1;
      }
      else {                                      /* Adjust first remote part of vertex array */
        Gnum              fldvertlocnum;
        Gnum              fldvertlocnnd;
        Gnum              fldvertlocadj;

        Gnum * restrict const fldvertloctax = fldgrafptr->s.vertloctax;

        fldvertlocnum = fldvertidxtab[j];
        fldvertlocadj = fldedgeidxtab[j] - fldgrafptr->s.vertloctax[fldvertlocnum];
        fldvendidxtab[j] = fldvertlocadj;         /* Record updated adjust value for vendloctab pass */

        for (fldvertlocnnd = fldvertlocnum + fldcommdattab[j].vertnbr; fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
          fldvertloctax[fldvertlocnum] += fldvertlocadj;
      }
    }

    for (i = 0; i < commnbr; i ++) {
      int               j;

      if (MPI_Waitany (commnbr, &requtab[HDGRAPHFOLDTAGVEND * commmax], &j, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (15)");
        cheklocval = 1;
      }
      else {                                      /* Adjust first remote part of vertex array */
        Gnum              fldvendlocnum;
        Gnum              fldvendlocnnd;
        Gnum              fldvendlocadj;

        Gnum * restrict const fldvendloctax = fldgrafptr->s.vendloctax;

        fldvendlocnum = fldvertidxtab[j];
        fldvendlocadj = fldvendidxtab[j];         /* Get updated adjust from above vertloctab pass */

        for (fldvendlocnnd = fldvendlocnum + fldcommdattab[j].vertnbr; fldvendlocnum < fldvendlocnnd; fldvendlocnum ++)
          fldvendloctax[fldvendlocnum] += fldvendlocadj;
      }
    }

    for (i = 0; i < commnbr; i ++) {
      MPI_Status        statdat;
      int               j;

      if (MPI_Waitany (commnbr, &requtab[HDGRAPHFOLDTAGEDGE * commmax], &j, &statdat) != MPI_SUCCESS) {
        errorPrint ("hdgraphFold2: communication error (16)");
        cheklocval = 1;
      }
      else if (cheklocval == 0) {                 /* Adjust remote part(s) of edge array */
        Gnum              orgvertlocmin;
        Gnum              orgvertlocmax;
        Gnum              fldvertlocnum;
        Gnum              fldvertlocnnd;
        Gnum              fldvertlocadj;
        Gnum              fldvhallocmax;          /* Maximum current size of halo vertex array */
        int               procngbmin;
        int               procngbmax;

        Gnum * restrict const fldvertloctax = fldgrafptr->s.vertloctax;
        Gnum * restrict const fldvendloctax = fldgrafptr->s.vendloctax;
        Gnum * restrict const fldedgeloctax = fldgrafptr->s.edgeloctax;

#ifdef SCOTCH_DEBUG_HDGRAPH2
        int               fldedgercvnbr;

        MPI_Get_count (&statdat, GNUM_MPI, &fldedgercvnbr);
        if (fldedgercvnbr != fldedgecnttab[j]) {
          errorPrint  ("hdgraphFold2: internal error (6)");
          return      (1);
        }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

        for (procngbmin = 0, procngbmax = fldvertadjnbr; /* Initialize search accelerator */
             procngbmax - procngbmin > 1; ) {
          int               procngbmed;

          procngbmed = (procngbmax + procngbmin) / 2;
          if (fldvertadjtab[procngbmed] <= fldcommvrttab[j])
            procngbmin = procngbmed;
          else
            procngbmax = procngbmed;
        }
        orgvertlocmin = fldvertadjtab[procngbmin];
        orgvertlocmax = fldvertadjtab[procngbmax];
        fldvertlocadj = fldvertdlttab[procngbmin];
        fldvhallocmax = orggrafptr->s.baseval - 1; /* Reset halo vertex array for each remote part                     */
        for (fldvertlocnum = fldvertidxtab[j], fldedgelocnum = fldedgeidxtab[j], /* Update received part of edge array */
             fldvertlocnnd = fldvertlocnum + fldcommdattab[j].vertnbr;
             fldvertlocnum < fldvertlocnnd; ) {
          for ( ; fldedgelocnum < fldvendloctax[fldvertlocnum]; fldedgelocnum ++) { /* Reorder end vertices */
            Gnum              orgvertlocend;

#ifdef SCOTCH_DEBUG_HDGRAPH2
            if (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval)) {
              errorPrint  ("hdgraphFold2: internal error (7)");
              return      (1);
            }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

            orgvertlocend = fldedgeloctax[fldedgelocnum];

            if ((orgvertlocend >= orgvertlocmin) && /* If end vertex is local */
                (orgvertlocend <  orgvertlocmax))
              fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertlocadj;
            else {
              int               procngbnum;
              int               procngbmax;

              for (procngbnum = 0, procngbmax = fldvertadjnbr;
                   procngbmax - procngbnum > 1; ) {
                int               procngbmed;

                procngbmed = (procngbmax + procngbnum) / 2;
                if (fldvertadjtab[procngbmed] <= orgvertlocend)
                  procngbnum = procngbmed;
                else
                  procngbmax = procngbmed;
              }
              fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertdlttab[procngbnum];
            }
          }
          fldvertlocnum ++;
          fldehallocnbr += fldvertloctax[fldvertlocnum] - fldedgelocnum;
          for ( ; fldedgelocnum < fldvertloctax[fldvertlocnum]; fldedgelocnum ++) { /* Renumber halo part */
            Gnum              orgverthalend;
            Gnum              fldvhallocend;

            orgverthalend = fldedgeloctax[fldedgelocnum];
#ifdef SCOTCH_DEBUG_HDGRAPH2
            if ((orgverthalend < orggrafptr->s.baseval) ||
                (orgverthalend >= (orggrafptr->s.edgeglbsmx + orggrafptr->s.baseval)) ||
                (fldedgelocnum >= (fldgrafptr->s.edgelocsiz + orggrafptr->s.baseval))) {
              errorPrint ("hdgraphFold2: internal error (8)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

            while (fldvhallocmax < orgverthalend) /* Expand halo vertex index array whenever necessary */
              fldvhalloctax[++ fldvhallocmax] = ~0;
            fldvhallocend = fldvhalloctax[orgverthalend]; /* Get renumbered halo vertex */
            if (fldvhallocend < 0) {              /* If new halo vertex not yet given   */
              fldvhallocend                =      /* Allocate it                        */
              fldvhalloctax[orgverthalend] = fldvhallocnum ++;
            }
            fldedgeloctax[fldedgelocnum] = fldvhallocend;
          }
        }
      }
    }

    if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver, edge arrays may have been oversized */
      fldgrafptr->s.edgeloctax  = memRealloc (fldgrafptr->s.edgeloctax + orggrafptr->s.baseval, fldgrafptr->s.edgelocsiz * sizeof (Gnum));
      fldgrafptr->s.edgeloctax -= orggrafptr->s.baseval;
    }

    fldgrafptr->vhallocnbr = fldvhallocnum - orggrafptr->s.baseval;
    fldgrafptr->vhndloctax = fldgrafptr->s.vertloctax + 1; /* Compact edge array with halo vertices */
    fldgrafptr->ehallocnbr = fldehallocnbr;
    fldgrafptr->levlnum    = orggrafptr->levlnum; /* Folded graph is of same level */

    if (orggrafptr->s.veloloctax == NULL)         /* If no vertex loads, reset graph vertex load to number of vertices */
      fldvelolocsum = fldvertlocnbr;
    else {                                        /* Graph has vertex loads and load of local part has already been computed */
      for (i = 0; i < commnbr; i ++) {
        int               j;

        if (MPI_Waitany (commnbr, &requtab[HDGRAPHFOLDTAGVELO * commmax], &j, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
          errorPrint ("hdgraphFold2: communication error (17)");
          cheklocval = 1;
        }
        else if (cheklocval == 0) {               /* Accumulate vertex loads for received vertex load array */
          Gnum              fldvertlocnum;
          Gnum              fldvertlocnnd;

          for (fldvertlocnum = fldvertidxtab[j], fldvertlocnnd = fldvertlocnum + fldcommdattab[j].vertnbr;
               fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
            fldvelolocsum += fldgrafptr->s.veloloctax[fldvertlocnum];
        }
      }
    }

    fldgrafptr->s.baseval    = orggrafptr->s.baseval;
    fldgrafptr->s.vertlocnbr = fldvertlocnbr;
    fldgrafptr->s.vertlocnnd = fldvertlocnbr + orggrafptr->s.baseval;
    fldgrafptr->s.velolocsum = fldvelolocsum;
    fldgrafptr->s.edgelocnbr = fldgrafptr->s.edgelocsiz - fldehallocnbr;
    fldgrafptr->s.degrglbmax = orggrafptr->s.degrglbmax;
    if (dgraphBuild4 (&fldgrafptr->s) != 0) {
      errorPrint  ("hdgraphFold2: cannot build folded graph");
      hdgraphExit (fldgrafptr);
      return      (1);
    }

#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (hdgraphCheck (fldgrafptr) != 0) {         /* Check graph consistency; vnumloctab is not checked so no need to wait for it */
      errorPrint  ("hdgraphFold2: internal error (9)");
      hdgraphExit (fldgrafptr);
      return      (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
  }

  memFree (fldcommdattab);                        /* Free group leader */

  if (MPI_Waitall (requnbr, requtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait for all graph data to arrive because graph could be freed afterwards */
    errorPrint ("hdgraphFold2: communication error (18)");
    cheklocval = 1;
  }

  memFree (fldvertidxtab);                        /* Free group leader including request array */

#ifdef SCOTCH_DEBUG_HDGRAPH1                      /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, orggrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphFold2: communication error (19)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_HDGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_HDGRAPH1 */

  return (chekglbval);
}
