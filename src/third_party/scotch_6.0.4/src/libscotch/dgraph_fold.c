/* Copyright 2007-2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_fold.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the distributed     **/
/**                source graph folding function.          **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 aug 2006     **/
/**                                 to   : 27 jun 2008     **/
/**                # Version 5.1  : from : 12 nov 2008     **/
/**                                 to   : 04 jan 2011     **/
/**                # Version 6.0  : from : 28 sep 2014     **/
/**                                 to   : 28 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_fold.h"
#include "dgraph_fold_comm.h"

/******************************/
/*                            */
/* This routine handles       */
/* distributed source graphs. */
/*                            */
/******************************/

/* This routine builds a folded graph by merging graph
** data to the processes of the first half or to the
** second half of the communicator.
** The key value of the folded communicator is not
** changed as it is not relevant.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphFold (
const Dgraph * restrict const orggrafptr,
const int                     partval,            /* 0 for first half, 1 for second half */
Dgraph * restrict const       fldgrafptr,
const void * restrict const   orgdataptr,         /* Un-based array of data which must be folded, e.g. coarmulttab */
void ** restrict const        flddataptr,         /* Un-based array of data which must be folded, e.g. coarmulttab */
MPI_Datatype                  datatype)
{
  int               fldprocnbr;
  int               fldprocnum;                   /* Index of local process in folded communicator   */
  int               fldproccol;                   /* Color of receiver or not wanted in communicator */
  MPI_Comm          fldproccomm;                  /* Communicator of folded part                     */
  int               o;

  fldprocnbr = (orggrafptr->procglbnbr + 1) / 2;
  fldprocnum = orggrafptr->proclocnum;
  if (partval == 1) {
    fldprocnum = fldprocnum - fldprocnbr;
    fldprocnbr = orggrafptr->procglbnbr - fldprocnbr;
  }
  fldproccol = ((fldprocnum >= 0) && (fldprocnum < fldprocnbr)) ? 0 : MPI_UNDEFINED;

  if (MPI_Comm_split (orggrafptr->proccomm, fldproccol, fldprocnum, &fldproccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFold: communication error");
    return     (1);
  }

  o = dgraphFold2 (orggrafptr, partval, fldgrafptr, fldproccomm, orgdataptr, flddataptr, datatype);
  fldgrafptr->prockeyval = fldproccol;            /* Key of folded communicator is always zero if no duplication occurs */

  return (o);
}

int
dgraphFold2 (
const Dgraph * restrict const orggrafptr,
const int                     partval,            /* 0 for first half, 1 for second half */
Dgraph * restrict const       fldgrafptr,
MPI_Comm                      fldproccomm,
const void * restrict const   orgdataptr,         /* Un-based array of data which must be kept, e.g. coarmulttab */
void ** restrict const        flddataptr,         /* Un-based array of data which must be kept, e.g. coarmulttab */
MPI_Datatype                  datatype)
{
  int                           fldcommtypval;    /* Type of communication for this process                 */
  DgraphFoldCommData * restrict fldcommdattab;    /* Array of two communication data                        */
  Gnum * restrict               fldcommvrttab;    /* Starting global send indices of communications         */
  Gnum * restrict               fldvertidxtab;    /* Start indices of vertex arrays                         */
  Gnum * restrict               fldedgeidxtab;    /* Start indices of edge arrays                           */
  Gnum * restrict               fldedgecnttab;    /* Number of edges exchanged during each communication    */
  Gnum * restrict               fldedgecnptab;    /* Temporary save for fldedgecnttab for MPI standard      */
  Gnum                          fldvertlocnbr;    /* Number of vertices in local folded part                */
  Gnum                          fldedgelocsiz;    /* (Upper bound of) number of edges in folded graph       */
  Gnum                          fldedlolocsiz;    /* (Upper bound of) number of edge loads in folded graph  */
  int                           fldprocglbnbr;
  int                           fldproclocnum;    /* Index of local process in folded communicator          */
  int                           fldvertadjnbr;
  Gnum * restrict               fldvertadjtab;    /* Array of global start indices for adjustment slots     */
  Gnum * restrict               fldvertdlttab;    /* Array of index adjustments for original global indices */
  int                           cheklocval;
  int                           chekglbval;
  int                           commmax;
  int                           commnbr;
  int                           requnbr;
  MPI_Request * restrict        requtab;
  int                           infosiz;          /* Size of one information                                */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (orggrafptr->vendloctax != (orggrafptr->vertloctax + 1)) {
    errorPrint ("dgraphFold2: graph must be compact");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  fldprocglbnbr = (orggrafptr->procglbnbr + 1) / 2;
  if (partval == 1) {
    fldproclocnum = orggrafptr->proclocnum - fldprocglbnbr;
    fldprocglbnbr = orggrafptr->procglbnbr - fldprocglbnbr;
  }
  else
    fldproclocnum = orggrafptr->proclocnum;

  fldcommtypval = ((fldproclocnum >= 0) && (fldproclocnum < fldprocglbnbr)) ? DGRAPHFOLDCOMMRECV : DGRAPHFOLDCOMMSEND;
  if (orgdataptr != NULL)
    MPI_Type_size (datatype, &infosiz);

  cheklocval    = 0;
  fldcommdattab = NULL;
  fldvertidxtab = NULL;
  if (fldcommtypval == DGRAPHFOLDCOMMRECV) {      /* If we are going to receive */
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (fldgrafptr == NULL) {
      errorPrint ("dgraphFold2: invalid parameters (1)");
      return     (1);
    }
    if (fldproccomm == MPI_COMM_NULL) {
      errorPrint ("dgraphFold2: invalid parameters (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    memSet (fldgrafptr, 0, sizeof (Dgraph));      /* Pre-initialize graph fields */

    fldgrafptr->proccomm   = fldproccomm;
    fldgrafptr->procglbnbr = fldprocglbnbr;
    fldgrafptr->proclocnum = fldproclocnum;
    fldgrafptr->flagval    = DGRAPHFREEALL | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* For premature freeing on error */

    if (memAllocGroup ((void **) (void *)         /* Allocate distributed graph private data */
                       &fldgrafptr->procdsptab, (size_t) ((fldprocglbnbr + 1) * sizeof (Gnum)),
                       &fldgrafptr->proccnttab, (size_t) (fldprocglbnbr       * sizeof (Gnum)),
                       &fldgrafptr->procngbtab, (size_t) (fldprocglbnbr       * sizeof (int)),
                       &fldgrafptr->procrcvtab, (size_t) (fldprocglbnbr       * sizeof (int)),
                       &fldgrafptr->procsndtab, (size_t) (fldprocglbnbr       * sizeof (int)), NULL) == NULL) {
      errorPrint ("dgraphFold2: out of memory (1)");
      cheklocval = 1;
    }
    else if (dgraphFoldComm (orggrafptr, partval, &commmax, &fldcommtypval, &fldcommdattab, &fldcommvrttab, /* Process can become a sender receiver */
                             fldgrafptr->proccnttab, &fldvertadjnbr, &fldvertadjtab, &fldvertdlttab) != 0) {
      errorPrint ("dgraphFold2: cannot compute folding communications (1)");
      cheklocval = 1;
    }
    else {
      Gnum              fldvelolocnbr;

      if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver */
        int               i;

        for (i = 0, fldvertlocnbr = 0; (i < commmax) && (fldcommdattab[i].procnum != -1); i ++)
          fldvertlocnbr += fldcommdattab[i].vertnbr;
        commnbr = i;

        fldedgelocsiz = orggrafptr->edgeglbsmx * i; /* Upper bound on received edges */
        if ((orggrafptr->degrglbmax > 0) && (fldvertlocnbr < (fldedgelocsiz / orggrafptr->degrglbmax))) /* Choose best upper bound on number of edges (avoid multiply overflow) */
          fldedgelocsiz = fldvertlocnbr * orggrafptr->degrglbmax;

        fldedgelocsiz += orggrafptr->edgelocnbr;  /* Add local edges and vertices */
        fldvertlocnbr += orggrafptr->vertlocnbr;
      }
      else {                                      /* Process is a sender receiver */
        fldvertlocnbr = fldcommvrttab[0] - orggrafptr->procvrttab[orggrafptr->proclocnum]; /* Communications will remove vertices   */
        fldedgelocsiz = orggrafptr->vertloctax[fldvertlocnbr + orggrafptr->baseval] - orggrafptr->baseval; /* Exact number of edges */

        fldgrafptr->edgelocnbr =
        fldgrafptr->edgelocsiz = fldedgelocsiz;
      }
      fldvelolocnbr = (orggrafptr->veloloctax != NULL) ? fldvertlocnbr : 0;

      if (memAllocGroup ((void **) (void *)       /* Allocate distributed graph public data */
                         &fldgrafptr->vertloctax, (size_t) ((fldvertlocnbr + 1) * sizeof (Gnum)),
                         &fldgrafptr->vnumloctax, (size_t) ( fldvertlocnbr      * sizeof (Gnum)),
                         &fldgrafptr->veloloctax, (size_t) ( fldvelolocnbr      * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("dgraphFold2: out of memory (2)");
        cheklocval = 1;
      }
      else if (fldgrafptr->vertloctax -= orggrafptr->baseval,
               fldgrafptr->vnumloctax -= orggrafptr->baseval,
               fldgrafptr->vendloctax  = fldgrafptr->vertloctax + 1, /* Folded graph is compact */
               fldgrafptr->veloloctax  = ((orggrafptr->veloloctax != NULL) ? (fldgrafptr->veloloctax - orggrafptr->baseval) : NULL),
               fldedlolocsiz = ((orggrafptr->edloloctax != NULL) ? fldedgelocsiz : 0),
               (fldgrafptr->edgeloctax = memAlloc ((fldedgelocsiz + fldedlolocsiz) * sizeof (Gnum))) == NULL) { /* Allocate single array for both edge arrays */
        errorPrint ("dgraphFold2: out of memory (3)");
        cheklocval = 1;
      }
      else {
        if (orgdataptr != NULL) {
          if ((*flddataptr = (byte *) memAlloc (fldvertlocnbr * infosiz)) == NULL) {
            errorPrint ("dgraphFold2: out of memory (4)");
            cheklocval = 1;
          }
        }
        fldgrafptr->edgeloctax -= orggrafptr->baseval; /* Do not care about the validity of edloloctax at this stage */
      }
    }
  }
  else {                                          /* Process is a sender */
#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (fldproccomm != MPI_COMM_NULL) {
      errorPrint ("dgraphFold2: invalid parameters (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

    if (dgraphFoldComm (orggrafptr, partval, &commmax, &fldcommtypval, &fldcommdattab, &fldcommvrttab, NULL, NULL, NULL, NULL) != 0) {
      errorPrint ("dgraphFold2: cannot compute folding communications (2)");
      cheklocval = 1;
    }
  }

  if ((cheklocval == 0) &&
      (memAllocGroup ((void **) (void *)          /* Allocate folding data */
                      &fldvertidxtab, (size_t) (commmax * sizeof (Gnum)),
                      &fldedgeidxtab, (size_t) (commmax * sizeof (Gnum)),
                      &fldedgecnttab, (size_t) (commmax * sizeof (Gnum)),
                      &fldedgecnptab, (size_t) (commmax * sizeof (Gnum)),
                      &requtab,       (size_t) (commmax * DGRAPHFOLDTAGNBR * sizeof (MPI_Request)), NULL) == NULL)) {
    errorPrint ("dgraphFold2: out of memory (5)");
    cheklocval = 1;
  }

#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, orggrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFold2: communication error (1)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0) {
    if ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) {
      if (fldvertidxtab != NULL)
        memFree (fldvertidxtab);                  /* Free group leader */
      if (fldcommdattab != NULL)
        memFree (fldcommdattab);
      dgraphExit (fldgrafptr);
    }
    return (1);
  }

  requnbr = 0;                                    /* Communications without further processing are placed at beginning of array */

  if ((fldcommtypval & DGRAPHFOLDCOMMSEND) != 0) { /* If process is (also) a sender */
    Gnum              vertsndbas;
    Gnum              vertsndnbr;
    int               i;

    vertsndnbr = ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) ? (fldcommvrttab[0] - orggrafptr->procvrttab[orggrafptr->proclocnum]) : 0; /* If process is also a receiver, start sending after kept vertices */

    for (i = 0, vertsndbas = orggrafptr->baseval; /* For all send communications to perform */
         (i < commmax) && (fldcommdattab[i].procnum != -1) && (cheklocval == 0); i ++) {
      vertsndbas += vertsndnbr;
      vertsndnbr  = fldcommdattab[i].vertnbr;

      fldvertidxtab[i] = vertsndbas;
      fldedgeidxtab[i] = orggrafptr->vertloctax[vertsndbas];
      fldedgecnptab[i] =                          /* Save fldedgecnttab in temporary array to read it while MPI communication in progress */
      fldedgecnttab[i] = orggrafptr->vertloctax[vertsndbas + vertsndnbr] - orggrafptr->vertloctax[vertsndbas]; /* Graph is compact        */
      if (MPI_Isend (&fldedgecnptab[i], 1, GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGVLBLLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("dgraphFold2: communication error (2)");
        cheklocval = 1;
      }
    }
    commnbr = i;

    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      if (MPI_Isend (orggrafptr->vertloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGVERTLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("dgraphFold2: communication error (3)");
        cheklocval = 1;
      }
    }
    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      if (MPI_Isend (orggrafptr->edgeloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, fldcommdattab[i].procnum,
                     TAGFOLD + TAGEDGELOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("dgraphFold2: communication error (4)");
        cheklocval = 1;
      }
    }
    if (orggrafptr->veloloctax != NULL) {
      for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
        if  (MPI_Isend (orggrafptr->veloloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                        TAGFOLD + TAGVELOLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
          errorPrint ("dgraphFold2: communication error (5)");
          cheklocval = 1;
        }
      }
    }
    for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
      int               procsndnum;               /* Rank of process to send to */

      procsndnum = fldcommdattab[i].procnum;
      if ((orggrafptr->edloloctax != NULL) &&
          (MPI_Isend (orggrafptr->edloloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, procsndnum,
                      TAGFOLD + TAGEDLOLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
        errorPrint ("dgraphFold2: communication error (6)");
        cheklocval = 1;
      }
      else if ((orggrafptr->vnumloctax != NULL) &&
               (MPI_Isend (orggrafptr->vnumloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, procsndnum,
                           TAGFOLD + TAGVNUMLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
        errorPrint ("dgraphFold2: communication error (7)");
        cheklocval = 1;
      }
      else if ((orgdataptr != NULL)  &&
               (MPI_Isend ((byte *) orgdataptr + ((fldvertidxtab[i] - orggrafptr->baseval) * infosiz), fldcommdattab[i].vertnbr, datatype, procsndnum,
                           TAGFOLD + TAGDATALOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
        errorPrint ("dgraphFold2: communication error (8)");
        cheklocval = 1;
      }
    }
  }                                               /* Communications of sender-receivers will be completed in the receiving phase */

  if ((fldcommtypval & DGRAPHFOLDCOMMRECV) != 0) { /* If process is (also) a receiver */
    Gnum                orgvertlocnbr;
    Gnum                orgvertlocnnd;
    Gnum                orgvertlocmin;
    Gnum                orgvertlocmax;
    Gnum                fldvertlocadj;
    Gnum                fldvelolocsum;
    Gnum                fldedgelocnum;
    Gnum                fldedgelocnnd;
    int                 fldprocnum;
    int                 procngbmin;
    int                 procngbmax;
    int                 i;

    const Gnum * restrict const orgedgeloctax = orggrafptr->edgeloctax;
    Gnum * restrict const       fldedgeloctax = fldgrafptr->edgeloctax;

    fldgrafptr->procvrttab = fldgrafptr->procdsptab; /* Graph does not have holes                               */
    fldgrafptr->procdsptab[0] = orggrafptr->baseval; /* Build private data of folded graph and array            */
    for (fldprocnum = 0; fldprocnum < fldprocglbnbr; fldprocnum ++) /* New subdomain indices start from baseval */
      fldgrafptr->procdsptab[fldprocnum + 1] = fldgrafptr->procdsptab[fldprocnum] + fldgrafptr->proccnttab[fldprocnum];

    if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver */
      Gnum                fldedgelocbas;
      Gnum                fldvertrcvbas;
      Gnum                fldvertrcvnbr;

      for (i = 0, fldvertrcvbas = orggrafptr->vertlocnnd, fldvertrcvnbr = 0; /* For all receive communications to perform */
           (i < commnbr) && (cheklocval == 0); i ++) {
        fldvertrcvbas += fldvertrcvnbr;
        fldvertrcvnbr  = fldcommdattab[i].vertnbr;

        fldvertidxtab[i] = fldvertrcvbas;
        if (MPI_Irecv (&fldedgecnttab[i], 1, GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGVLBLLOCTAB, orggrafptr->proccomm, &requtab[DGRAPHFOLDTAGENBR * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("dgraphFold2: communication error (9)");
          cheklocval = 1;
        }
      }

      for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) { /* Let these communications progress while we process the edge size messages */
        if (MPI_Irecv (fldgrafptr->vertloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGVERTLOCTAB, orggrafptr->proccomm, &requtab[DGRAPHFOLDTAGVERT * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("dgraphFold2: communication error (10)");
          cheklocval = 1;
        }
      }

      MPI_Waitall (commnbr, &requtab[DGRAPHFOLDTAGENBR * commmax], MPI_STATUSES_IGNORE);

      for (i = 0, fldedgelocbas = orggrafptr->vertloctax[orggrafptr->vertlocnnd]; (i < commnbr) && (cheklocval == 0); i ++) {
        fldedgeidxtab[i] = fldedgelocbas;
        fldedgelocbas += fldedgecnttab[i];

        if (MPI_Irecv (fldgrafptr->edgeloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, fldcommdattab[i].procnum,
                       TAGFOLD + TAGEDGELOCTAB, orggrafptr->proccomm, &requtab[DGRAPHFOLDTAGEDGE * commmax + i]) != MPI_SUCCESS) {
          errorPrint ("dgraphFold2: communication error (11)");
          cheklocval = 1;
        }
      }
      fldgrafptr->edgelocnbr =                    /* Get number of local edges */
      fldgrafptr->edgelocsiz = fldedgelocbas - orggrafptr->baseval;

      if (orggrafptr->veloloctax != NULL) {
        for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
          if (MPI_Irecv (fldgrafptr->veloloctax + fldvertidxtab[i], fldcommdattab[i].vertnbr, GNUM_MPI, fldcommdattab[i].procnum,
                         TAGFOLD + TAGVELOLOCTAB, orggrafptr->proccomm, &requtab[DGRAPHFOLDTAGVELO * commmax + i]) != MPI_SUCCESS) {
            errorPrint ("dgraphFold2: communication error (12)");
            cheklocval = 1;
          }
        }
      }
      if (orggrafptr->edloloctax != NULL) {
        fldgrafptr->edloloctax = fldgrafptr->edgeloctax + fldgrafptr->edgelocnbr; /* Set start index of edge load array */

        for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
          if (MPI_Irecv (fldgrafptr->edloloctax + fldedgeidxtab[i], fldedgecnttab[i], GNUM_MPI, fldcommdattab[i].procnum,
                         TAGFOLD + TAGEDLOLOCTAB, orggrafptr->proccomm, &requtab[DGRAPHFOLDTAGEDLO * commmax + i]) != MPI_SUCCESS) {
            errorPrint ("dgraphFold2: communication error (13)");
            cheklocval = 1;
          }
        }
      }
      for (i = 0; (i < commnbr) && (cheklocval == 0); i ++) {
        int               procrcvnum;             /* Rank of process to receive from */
        Gnum              vertrcvnbr;

        procrcvnum = fldcommdattab[i].procnum;
        vertrcvnbr = fldcommdattab[i].vertnbr;
        if ((orggrafptr->vnumloctax != NULL) &&
            (MPI_Irecv (fldgrafptr->vnumloctax + fldvertidxtab[i], vertrcvnbr, GNUM_MPI, procrcvnum,
                        TAGFOLD + TAGVNUMLOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
          errorPrint ("dgraphFold2: communication error (14)");
          cheklocval = 1;
        }
        else if ((orgdataptr != NULL) &&
                 (MPI_Irecv ((byte *) (*flddataptr) + ((fldvertidxtab[i] - orggrafptr->baseval) * infosiz), vertrcvnbr, datatype, procrcvnum,
                             TAGFOLD + TAGDATALOCTAB, orggrafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS)) {
          errorPrint ("dgraphFold2: communication error (15)");
          cheklocval = 1;
        }
      }

      orgvertlocnbr = orggrafptr->vertlocnbr;     /* Process all local vertices */
      orgvertlocnnd = orggrafptr->vertlocnnd;

      if (orggrafptr->vnumloctax == NULL) {       /* If original graph does not have vertex numbers, create remote parts of vertex number array */
        Gnum              fldvertlocnum;
        Gnum              fldvertlocadj;
        int               i;

        Gnum * restrict const fldvnumloctax = fldgrafptr->vnumloctax;

        for (i = 0, fldvertlocnum = orgvertlocnnd; i < commnbr; i ++) {
          Gnum              fldvertlocnnd;

          for (fldvertlocnnd = fldvertlocnum + fldcommdattab[i].vertnbr, fldvertlocadj = fldcommvrttab[i];
               fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
            fldvnumloctax[fldvertlocnum] = fldvertlocadj ++;
        }
      }

      fldedgelocnnd = orggrafptr->vertloctax[orggrafptr->vertlocnnd];
      fldvelolocsum = orggrafptr->velolocsum;     /* In case there are vertex loads, we keep all of existing load */
    }
    else {                                        /* Receiver process is also a sender     */
      orgvertlocnbr = fldvertlocnbr;              /* Process only remaining local vertices */
      orgvertlocnnd = fldvertlocnbr + orggrafptr->baseval;

      if (orggrafptr->veloloctax != NULL) {       /* If original graph has vertex loads */
        Gnum                fldvertlocnum;

        for (fldvertlocnum = orggrafptr->baseval, fldvelolocsum = 0; /* Accumulate load sum of remaining part */
             fldvertlocnum < orgvertlocnnd; fldvertlocnum ++)
          fldvelolocsum += orggrafptr->veloloctax[fldvertlocnum];
      }
      fldedgelocnnd = orggrafptr->vertloctax[orgvertlocnnd]; /* Reorder remaining local part of edge array */
      if (orggrafptr->edloloctax != NULL)
        fldgrafptr->edloloctax = fldgrafptr->edgeloctax + fldgrafptr->edgelocnbr; /* Set start index of edge load array */

      commnbr = 0;                                /* Turn sender-receiver into normal receiver without any communications to perform */
    }

    for (procngbmin = 0, procngbmax = fldvertadjnbr; /* Initialize search accelerator */
         procngbmax - procngbmin > 1; ) {
      int               procngbmed;

      procngbmed = (procngbmax + procngbmin) / 2;
      if (fldvertadjtab[procngbmed] <= orggrafptr->procvrttab[orggrafptr->proclocnum])
        procngbmin = procngbmed;
      else
        procngbmax = procngbmed;
    }
    orgvertlocmin = fldvertadjtab[procngbmin];
    orgvertlocmax = fldvertadjtab[procngbmax];
    fldvertlocadj = fldvertdlttab[procngbmin];
    for (fldedgelocnum = orggrafptr->baseval; fldedgelocnum < fldedgelocnnd; fldedgelocnum ++) {
      Gnum              orgvertlocend;

      orgvertlocend = orgedgeloctax[fldedgelocnum];

      if ((orgvertlocend >= orgvertlocmin) &&     /* If end vertex is local */
          (orgvertlocend <  orgvertlocmax))
        fldedgeloctax[fldedgelocnum] = orgvertlocend + fldvertlocadj;
      else {                                      /* End vertex is not local */
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

    if (orggrafptr->veloloctax != NULL)           /* If original graph has vertex loads             */
      memCpy (fldgrafptr->veloloctax + orggrafptr->baseval, /* Copy local part of vertex load array */
              orggrafptr->veloloctax + orggrafptr->baseval, orgvertlocnbr * sizeof (Gnum));

    if (orggrafptr->edloloctax != NULL)           /* If original graph has edge loads             */
      memCpy (fldgrafptr->edloloctax + orggrafptr->baseval, /* Copy local part of edge load array */
              orggrafptr->edloloctax + orggrafptr->baseval,
              (orggrafptr->vertloctax[orgvertlocnnd] - orggrafptr->baseval) * sizeof (Gnum));

    if (orggrafptr->vnumloctax != NULL)           /* If original graph has vertex numbers             */
      memCpy (fldgrafptr->vnumloctax + orggrafptr->baseval, /* Copy local part of vertex number array */
              orggrafptr->vnumloctax + orggrafptr->baseval, orgvertlocnbr * sizeof (Gnum));
    else {                                        /* Build local part of vertex number array */
      Gnum              fldvertlocnum;
      Gnum              fldvertlocadj;

      for (fldvertlocnum = orggrafptr->baseval, fldvertlocadj = orggrafptr->procvrttab[orggrafptr->proclocnum];
           fldvertlocnum < orgvertlocnnd; fldvertlocnum ++)
        fldgrafptr->vnumloctax[fldvertlocnum] = fldvertlocadj ++;
    }

    memCpy (fldgrafptr->vertloctax + orggrafptr->baseval, /* Copy local part of vertex array, since it is compact     */
            orggrafptr->vertloctax + orggrafptr->baseval, orgvertlocnbr * sizeof (Gnum)); /* Last value is not copied */
    fldgrafptr->vertloctax[fldvertlocnbr + orggrafptr->baseval] = fldgrafptr->edgelocnbr + orggrafptr->baseval;

    if (orgdataptr != NULL)                       /* If additional data present */
      memCpy ((byte *) (*flddataptr), (byte *) orgdataptr, orgvertlocnbr * infosiz); /* Copy local part */

    for (i = 0; i < commnbr; i ++) {
      int               j;

      if (MPI_Waitany (commnbr, &requtab[DGRAPHFOLDTAGVERT * commmax], &j, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        errorPrint ("dgraphFold2: communication error (16)");
        cheklocval = 1;
      }
      else {                                      /* Adjust first remote part of vertex array */
        Gnum              fldvertlocnum;
        Gnum              fldvertlocnnd;
        Gnum              fldvertlocadj;

        Gnum * restrict const fldvertloctax = fldgrafptr->vertloctax;

        fldvertlocnum = fldvertidxtab[j];
        fldvertlocadj = fldedgeidxtab[j] - fldgrafptr->vertloctax[fldvertlocnum];

        for (fldvertlocnnd = fldvertlocnum + fldcommdattab[j].vertnbr; fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
          fldvertloctax[fldvertlocnum] += fldvertlocadj;
      }
    }

    for (i = 0; i < commnbr; i ++) {
      MPI_Status        statdat;
      int               j;

      if (MPI_Waitany (commnbr, &requtab[DGRAPHFOLDTAGEDGE * commmax], &j, &statdat) != MPI_SUCCESS) {
        errorPrint ("dgraphFold2: communication error (17)");
        cheklocval = 1;
      }
      else if (cheklocval == 0) {                 /* Adjust remote part(s) of edge array */
        Gnum              orgvertlocmin;
        Gnum              orgvertlocmax;
        Gnum              fldvertlocadj;
        int               procngbnum;
        int               procngbmax;

        Gnum * restrict const fldedgeloctax = fldgrafptr->edgeloctax;

#ifdef SCOTCH_DEBUG_DGRAPH2
        int               fldedgercvnbr;

        MPI_Get_count (&statdat, GNUM_MPI, &fldedgercvnbr);
        if (fldedgercvnbr != fldedgecnttab[j]) {
          errorPrint  ("dgraphFold2: internal error (1)");
          return      (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        
        for (procngbnum = 0, procngbmax = fldvertadjnbr; /* Initialize search accelerator */
             procngbmax - procngbnum > 1; ) {
          int               procngbmed;

          procngbmed = (procngbmax + procngbnum) / 2;
          if (fldvertadjtab[procngbmed] <= fldcommvrttab[j])
            procngbnum = procngbmed;
          else
            procngbmax = procngbmed;
        }
        orgvertlocmin = fldvertadjtab[procngbnum];
        orgvertlocmax = fldvertadjtab[procngbmax];
        fldvertlocadj = fldvertdlttab[procngbnum];
        for (fldedgelocnum = fldedgeidxtab[j], fldedgelocnnd = fldedgelocnum + fldedgecnttab[j];
             fldedgelocnum < fldedgelocnnd; fldedgelocnum ++) { /* Reorder end vertices */
          Gnum              orgvertlocend;

#ifdef SCOTCH_DEBUG_DGRAPH2
          if (fldedgelocnum >= (fldgrafptr->edgelocnbr + orggrafptr->baseval)) {
            errorPrint  ("dgraphFold2: internal error (2)");
            return      (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

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
      }
    }

    if (orggrafptr->veloloctax == NULL)           /* If no vertex loads, reset graph vertex load to number of vertices */
      fldvelolocsum = fldvertlocnbr;
    else {                                        /* Graph has vertex loads and load of local part has already been computed */
      for (i = 0; i < commnbr; i ++) {
        int               j;

        if (MPI_Waitany (commnbr, &requtab[DGRAPHFOLDTAGVELO * commmax], &j, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
          errorPrint ("dgraphFold2: communication error (18)");
          cheklocval = 1;
        }
        else if (cheklocval == 0) {               /* Accumulate vertex loads for received vertex load array */
          Gnum              fldvertlocnum;
          Gnum              fldvertlocnnd;

          for (fldvertlocnum = fldvertidxtab[j], fldvertlocnnd = fldvertlocnum + fldcommdattab[j].vertnbr;
               fldvertlocnum < fldvertlocnnd; fldvertlocnum ++)
            fldvelolocsum += fldgrafptr->veloloctax[fldvertlocnum];
        }
      }
    }

    if ((fldcommtypval & DGRAPHFOLDCOMMSEND) == 0) { /* If process is a normal receiver, edge arrays may have been oversized */
      Gnum                fldedgeloctmp;

      fldedgeloctmp = fldgrafptr->edgelocnbr;
      if (orggrafptr->edloloctax != NULL) {
        fldedgeloctmp *= 2;
        if (MPI_Waitall (commnbr, &requtab[DGRAPHFOLDTAGEDLO * commmax], MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait for edge load sub-arrays */
          errorPrint ("dgraphFold2: communication error (19)");
          cheklocval = 1;
        }
      }

      fldgrafptr->edgeloctax  = memRealloc (fldgrafptr->edgeloctax + orggrafptr->baseval, fldedgeloctmp * sizeof (Gnum));
      fldgrafptr->edgeloctax -= orggrafptr->baseval;
      if (orggrafptr->edloloctax != NULL)
        fldgrafptr->edloloctax = fldgrafptr->edgeloctax + fldgrafptr->edgelocnbr;
    }

    fldgrafptr->baseval    = orggrafptr->baseval;
    fldgrafptr->vertlocnbr = fldvertlocnbr;
    fldgrafptr->vertlocnnd = fldvertlocnbr + orggrafptr->baseval;
    fldgrafptr->velolocsum = fldvelolocsum;
    fldgrafptr->degrglbmax = orggrafptr->degrglbmax;
    if (dgraphBuild4 (fldgrafptr) != 0) {
      errorPrint ("dgraphFold2: cannot build folded graph");
      dgraphExit (fldgrafptr);
      return     (1);
    }

#ifdef SCOTCH_DEBUG_DGRAPH2
    if (dgraphCheck (fldgrafptr) != 0) {          /* Check graph consistency; vnumloctab is not checked so no need to wait for it */
      errorPrint ("dgraphFold2: internal error (3)");
      dgraphExit (fldgrafptr);
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  }

  memFree (fldcommdattab);                        /* Free group leader */

  if (MPI_Waitall (requnbr, requtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait for all graph data to arrive because graph could be freed afterwards */
    errorPrint ("dgraphFold2: communication error (20)");
    cheklocval = 1;
  }

  memFree (fldvertidxtab);                        /* Free group leader including request array */

#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, orggrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFold2: communication error (21)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  return (chekglbval);
}
