/* Copyright 2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_match_sync_coll.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This routine synchronizes the fragments **/
/**                of a distributed matching by means of   **/
/**                collective communications.              **/
/**                                                        **/
/**    DATES     : # Version 5.1  : from : 06 feb 2009     **/
/**                                 to   : 22 apr 2009     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_MATCH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "dgraph_match.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* source graphs.                    */
/*                                   */
/*************************************/

/* This routine performs a round of communication
** to synchronize enqueued matching requests across
** processors.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphMatchSyncColl (
DgraphMatchData * restrict const  mateptr)
{
  Gnum                queulocnbr;
  Gnum                queulocnum;
  Gnum                matelocnbr;
  Gnum                multlocnbr;
  Gnum                vertlocadj;
  Gnum                edgekptnbr;
  int                 procngbnbr;
  int                 procngbidx;
  int                 procngbnum;
  int * restrict      vsndcnttab;
  int * restrict      vrcvcnttab;
  int * restrict      vsnddsptab;
  int * restrict      vrcvdsptab;

  Dgraph * restrict const             grafptr    = mateptr->c.finegrafptr;
  const int * restrict const          procngbtab = grafptr->procngbtab;
  int * restrict const                procgsttax = mateptr->c.procgsttax;
  const Gnum * restrict const         procvgbtab = mateptr->procvgbtab;
  const Gnum * restrict const         vertloctax = grafptr->vertloctax;
  const Gnum * restrict const         vendloctax = grafptr->vendloctax;
  const Gnum * restrict const         edgeloctax = grafptr->edgeloctax;
  const Gnum * restrict const         edgegsttax = grafptr->edgegsttax;
  Gnum * restrict const               queuloctab = mateptr->queuloctab;
  Gnum * restrict const               mategsttax = mateptr->mategsttax;
  DgraphCoarsenMulti * restrict const multloctab = mateptr->c.multloctab;
  int * restrict const                nsndidxtab = mateptr->c.nsndidxtab;
  DgraphCoarsenVert * const           vsnddattab = mateptr->c.vsnddattab; /* [norestrict:async] */

  procngbnbr = grafptr->procngbnbr;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (edgeloctax == NULL) {
    errorPrint ("dgraphMatchSyncColl: not implemented");
    return     (1);
  }
  if (MPI_Barrier (grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchSyncColl: communication error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if ((vsnddsptab = memAlloc (4 * grafptr->procglbnbr * sizeof (int))) == NULL) {
    errorPrint ("dgraphMatchSyncColl: out of memory");
    return     (1);
  }
  vsndcnttab = vsnddsptab + grafptr->procglbnbr;  /* TRICK: put vsnddsptab, vsndcnttab, vrcvdsptab in order for memSet() */
  vrcvdsptab = vsndcnttab + grafptr->procglbnbr;
  vrcvcnttab = vrcvdsptab + grafptr->procglbnbr;

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) /* Reset indices for sending messages */
    nsndidxtab[procngbnum] = mateptr->c.vsnddsptab[procngbtab[procngbnum]];

  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;
  for (queulocnum = 0, queulocnbr = mateptr->queulocnbr; queulocnum < queulocnbr; queulocnum ++) {
    Gnum                vertlocnum;
    Gnum                vertgstnum;
    Gnum                edgelocnum;
    Gnum                mategstnum;
    Gnum                mateglbnum;
    int                 procngbnum;
    int                 vsndidxnum;

    vertlocnum = queuloctab[queulocnum];          /* Get local vertex index  */
    mategstnum = mategsttax[vertlocnum];          /* Get mate (edge ?) index */

    if (mategstnum >= -1)                         /* If vertex not willing to mate or matched locally after being considered during matching phase */
      continue;

    edgelocnum = -2 - mategstnum;                 /* Get local edge to mate ghost vertex */
#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((edgelocnum < grafptr->baseval) ||
        (edgelocnum >= (grafptr->edgelocsiz + grafptr->baseval)) ||
        (mategsttax[edgegsttax[edgelocnum]] != -1)) {
      errorPrint ("dgraphMatchSyncColl: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    mateglbnum = edgeloctax[edgelocnum];

    vertgstnum = edgegsttax[edgelocnum];
    procngbnum = procgsttax[vertgstnum];          /* Find neighbor owner process  */
    if (procngbnum < 0) {                         /* If neighbor not yet computed */
      int                 procngbmax;

      procngbnum = 0;
      procngbmax = procngbnbr;
      while ((procngbmax - procngbnum) > 1) {     /* Find owner process by dichotomy on procvgbtab */
        int                 procngbmed;

        procngbmed = (procngbmax + procngbnum) / 2;
        if (procvgbtab[procngbmed] > mateglbnum)
          procngbmax = procngbmed;
        else
          procngbnum = procngbmed;
      }
      procgsttax[vertgstnum] = procngbnum;
    }
#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((grafptr->procvrttab[procngbtab[procngbnum]]     >  mateglbnum) ||
        (grafptr->procvrttab[procngbtab[procngbnum] + 1] <= mateglbnum)) {
      errorPrint ("dgraphMatchSyncColl: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    vsndidxnum = nsndidxtab[procngbnum] ++;       /* Get position of message in send array */
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (vsndidxnum >= mateptr->c.vsnddsptab[procngbtab[procngbnum] + 1]) {
      errorPrint ("dgraphMatchSyncColl: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    vsnddattab[vsndidxnum].datatab[0] = vertlocnum + vertlocadj;
    vsnddattab[vsndidxnum].datatab[1] = mateglbnum;
  }

  memSet (vsnddsptab, 0, 3 * grafptr->procglbnbr * sizeof (int)); /* TRICK: resets vsnddsptab, vsndcnttab, vrcvdsptab */
  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = procngbtab[procngbnum];
    vrcvdsptab[procglbnum] = 2 * mateptr->c.vrcvdsptab[procglbnum]; /* Times 2 because a "DgraphCoarsenVert" is two "Gnum"s */
    vsnddsptab[procglbnum] = 2 * mateptr->c.vsnddsptab[procglbnum];
    vsndcnttab[procglbnum] = 2 * (nsndidxtab[procngbnum] - mateptr->c.vsnddsptab[procglbnum]);
  }

  if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, vrcvcnttab, 1, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchSyncColl: communication error (2)");
    return     (1);
  }
  if (MPI_Alltoallv (vsnddattab,            vsndcnttab, vsnddsptab, GNUM_MPI,
                     mateptr->c.vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchSyncColl: communication error (3)");
    return     (1);
  }

  matelocnbr = mateptr->matelocnbr;
  multlocnbr = mateptr->c.multlocnbr;
  edgekptnbr = mateptr->c.edgekptnbr;

  for (procngbidx = 0; procngbidx < procngbnbr; procngbidx ++) {
    int                 procngbnum;
    int                 procglbnum;
    int                 vrcvidxnnd;
    int                 requrcvnum;
    int                 requnxtnum;               /* Index of location where to pack requests to process when all messages arrive */

    procngbnum = (procngbidx + mateptr->c.procngbnxt) % procngbnbr;
    procglbnum = procngbtab[procngbnum];
    vrcvidxnnd = mateptr->c.vrcvdsptab[procglbnum];
    if (vrcvcnttab[procglbnum] > 0) {             /* If query message is not empty                        */
      Gnum                vertsndnbr;             /* Number of vertices to be sent to requesting neighbor */
      Gnum                edgesndnbr;             /* Number of edges to be sent to requesting neighbor    */

      DgraphCoarsenVert * restrict const  vrcvdattab = mateptr->c.vrcvdattab; /* Local restrict pointer only when data available */

      vertsndnbr =
      edgesndnbr = 0;
      for (requrcvnum = requnxtnum = vrcvidxnnd, vrcvidxnnd += (vrcvcnttab[procglbnum] / 2); /* TRICK: each message item costs 2 Gnum's */
           requrcvnum < vrcvidxnnd; requrcvnum ++) {
        Gnum                vertglbnum;           /* Our global number (the one seen as mate by sender)   */
        Gnum                vertlocnum;           /* Our local number (the one seen as mate by sender)    */
        Gnum                vmatglbnum;           /* Global number of requesting mate (sender of message) */
        Gnum                mategstnum;           /* The mate we wanted to ask for                        */

        vmatglbnum = vrcvdattab[requrcvnum].datatab[0]; /* Names are opposite because receiving side */
        vertglbnum = vrcvdattab[requrcvnum].datatab[1];
        vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocnum <  grafptr->baseval) ||   /* If matching request is not directed towards our process */
            (vertlocnum >= grafptr->vertlocnnd)) {
          errorPrint ("dgraphMatchSyncColl: internal error (5)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

        mategstnum = mategsttax[vertlocnum];      /* Get our local mating decision data        */
        if (mategstnum == -1) {                   /* If local vertex wanted for mating is free */
          Gnum                edgelocnum;

          for (edgelocnum = vertloctax[vertlocnum]; edgeloctax[edgelocnum] != vmatglbnum; edgelocnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH2
            if (edgelocnum >= vendloctax[vertlocnum]) {
              errorPrint ("dgraphMatchSyncColl: internal error (6)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          }
          mategsttax[edgegsttax[edgelocnum]] = vertglbnum; /* We are no longer free        */
          mategsttax[vertlocnum] = vmatglbnum;    /* Leave message as is to acknowledge it */
          matelocnbr ++;
          vertsndnbr ++;
          edgesndnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum];
        }
        else if (mategstnum < -1) {               /* If local vertex is also asking for mating */
          Gnum                edgelocnum;
          Gnum                mateglbnum;

          edgelocnum = -2 - mategstnum;
          mateglbnum = edgeloctax[edgelocnum];    /* Get global number of our remote mate       */
          if (mateglbnum == vmatglbnum) {         /* If it is with the sender                   */
            Gnum                flagval;          /* Flag for choosing side to create multinode */

            mategsttax[vertlocnum] = mateglbnum;  /* Say we are mated to inform future requesting processes in same pass */
            mategsttax[edgegsttax[edgelocnum]] = vertglbnum;
            flagval = (mateglbnum > vertglbnum) ? 1 : 0; /* Compute pseudo-random flag always opposite for both ends */
            flagval = ((mateglbnum + (mateglbnum - vertglbnum) * flagval) & 1) ^ flagval;
            if (flagval == 0) {                   /* If flag is even, create multinode */
              multloctab[multlocnbr].vertglbnum[0] = vertglbnum;
              multloctab[multlocnbr].vertglbnum[1] = mategstnum; /* Remote mate: negative value */
              multlocnbr ++;                      /* One more coarse vertex created             */
            edgekptnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum];
            }
            else {                                /* If flag is odd, prepare to send vertex data at build time */
              vertsndnbr ++;
              edgesndnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum];
            }                                     /* Go on by destroying message in all cases since both ends know what it is about */

            vrcvdattab[requrcvnum --] = vrcvdattab[-- vrcvidxnnd]; /* Replace current message with another one and process it */
            matelocnbr ++;                        /* One more local vertex mated on each side; no messages will tell it       */
          }
          else {                                  /* If willing to mate but not with the sender, tell later with whom */
            DgraphCoarsenVert   vertdat;          /* Temporary storage data for swapping vertices                     */

            vertdat = vrcvdattab[requnxtnum];     /* Pack requests to process later at beginning of message */
            vrcvdattab[requnxtnum].datatab[0] = vmatglbnum;
            vrcvdattab[requnxtnum].datatab[1] = -2 - vertlocnum; /* Build appropriate answer to mating request later, when all messages arrived */
            if (requnxtnum ++ != requrcvnum)
              vrcvdattab[requrcvnum] = vertdat;   /* Swap vertices if not already at the right place */
          }
        }
        else                                      /* If already matched, inform sender */
          vrcvdattab[requrcvnum].datatab[1] = mategstnum;
      }
      mateptr->c.dcntloctab[procglbnum].vertsndnbr += vertsndnbr;
      mateptr->c.dcntloctab[procglbnum].edgesndnbr += edgesndnbr;
    }
    mateptr->c.nrcvidxtab[procngbnum] = vrcvidxnnd;
  }

  for (procngbidx = 0; procngbidx < procngbnbr; procngbidx ++) {
    int                 procngbnum;
    int                 procglbnum;
    int                 vsndidxnnd;
    int                 vsndidxnum;

    DgraphCoarsenVert * restrict const  vrcvdattab = mateptr->c.vrcvdattab; /* Local restrict pointer only once data received */

    procngbnum = (procngbidx + mateptr->c.procngbnxt) % procngbnbr;
    procglbnum = procngbtab[procngbnum];
    vsndidxnnd = mateptr->c.nrcvidxtab[procngbnum]; /* Re-send the messages we have received to acknowledge */

    for (vsndidxnum = mateptr->c.vrcvdsptab[procglbnum]; /* Finalize unfinished messages */
         vsndidxnum < vsndidxnnd; vsndidxnum ++) {
      Gnum                vertlocnum;
      Gnum                mateglbnum;

      vertlocnum = vrcvdattab[vsndidxnum].datatab[1];
      if (vertlocnum >= 0)                        /* If no more unfinished messages to process, quit scanning */
        break;
      vertlocnum = -2 - vertlocnum;
      mateglbnum = mategsttax[vertlocnum];
      if (mateglbnum >= 0)                        /* If vertex we wanted to mate with has been mated in this round */
        vrcvdattab[vsndidxnum].datatab[1] = mateglbnum; /* Propagate this information back to the requester        */
      else {                                      /* Vertex mating data not yet available (maybe in answer)        */
        vrcvdattab[vsndidxnum] = vrcvdattab[-- vsndidxnnd]; /* Remove message as no reply means not willing        */
        if (vrcvdattab[vsndidxnum].datatab[1] < 0) /* If replacing message is also to be processed                 */
          vsndidxnum --;                          /* Do not skip replaced message in next iteration                */
      }
    }

    if (vsndidxnnd < (mateptr->c.vrcvdsptab[procglbnum] + (vrcvcnttab[procglbnum] / 2))) { /* If space created */
      vrcvdattab[vsndidxnnd].datatab[0] =         /* Indicate end of useful data */
      vrcvdattab[vsndidxnnd].datatab[1] = -1;
    }
  }

  if (MPI_Alltoallv (mateptr->c.vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI,
                     vsnddattab,            vsndcnttab, vsnddsptab, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchSyncColl: communication error (3)");
    return     (1);
  }

  for (procngbidx = 0; procngbidx < procngbnbr; procngbidx ++) {
    int                 procngbnum;
    int                 procglbnum;
    int                 vrcvidxnnd;
    int                 vrcvidxnum;

    procngbnum = (procngbidx + mateptr->c.procngbnxt) % procngbnbr;
    procglbnum = procngbtab[procngbnum];
    for (vrcvidxnum = mateptr->c.vsnddsptab[procglbnum],
         vrcvidxnnd = vrcvidxnum + (vsndcnttab[procglbnum] / 2); /* TRICK: each message item costs 2 Gnum's */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                edgelocnum;
      Gnum                vertglbnum;             /* Our global number (the one seen as mate by sender) */
      Gnum                vertlocnum;             /* Our local number (the one seen as mate by sender)  */
      Gnum                vmatglbnum;             /* Global number of vertex to which the mate is mated */
      Gnum                mategstnum;             /* The mate we wanted to ask for                      */

      vertglbnum = vsnddattab[vrcvidxnum].datatab[0];
      if (vertglbnum == -1)                       /* If end of useful space reached */
        break;

      vmatglbnum = vsnddattab[vrcvidxnum].datatab[1];
      vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((vertlocnum <  grafptr->baseval) ||     /* If matching reply is not directed towards our process */
          (vertlocnum >= grafptr->vertlocnnd)) {
        errorPrint ("dgraphMatchSyncColl: internal error (8)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      mategstnum = mategsttax[vertlocnum];        /* Get our local mating decision data */
      edgelocnum = -2 - mategstnum;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((mategstnum >= -1) ||                   /* If we did not ask anything or if we were already matched, no reply message should come to us */
          ((mategsttax[edgegsttax[edgelocnum]] >= 0) && /* Also, if our prospective mate was itself already set as matched by a previous reply    */
           (mategsttax[edgegsttax[edgelocnum]] != vertglbnum) && /* And this message is not the positive reply which acknowledges this mating     */
           (mategsttax[edgegsttax[edgelocnum]] != vmatglbnum))) { /* Or an informative negative reply which gives again the mate of the ghost     */
        errorPrint ("dgraphMatchSyncColl: internal error (9)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      if (edgeloctax[edgelocnum] == vmatglbnum) { /* If positive answer from the mate we wanted */
        mategsttax[vertlocnum] = vmatglbnum;      /* Set local vertex as matched with the mate  */
        mategsttax[edgegsttax[edgelocnum]] = vertglbnum; /* Update state of ghost mate          */
        multloctab[multlocnbr].vertglbnum[0] = vertglbnum;
        multloctab[multlocnbr].vertglbnum[1] = mategstnum; /* Remote mate: negative value */
        multlocnbr ++;                            /* One more coarse vertex created       */
        matelocnbr ++;
        edgekptnbr += vendloctax[vertlocnum] - vertloctax[vertlocnum];
      }
      else {                                      /* If negative answer from the mate we wanted  */
        mategsttax[vertlocnum] = -1;              /* Reset local vertex as free for mating       */
        mategsttax[edgegsttax[edgelocnum]] = vmatglbnum; /* Update state of unwilling ghost mate */
      }
    }
  }

  mateptr->matelocnbr   = matelocnbr;
  mateptr->c.multlocnbr = multlocnbr;
  mateptr->c.edgekptnbr = edgekptnbr;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Barrier (grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchSyncColl: communication error (11)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  memFree (vsnddsptab);

  return (0);
}
