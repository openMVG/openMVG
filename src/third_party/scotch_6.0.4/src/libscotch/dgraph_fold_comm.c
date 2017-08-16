/* Copyright 2007,2009-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_fold_comm.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes the communication  **/
/**                pattern of the distributed graph        **/
/**                folding process.                        **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 may 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 18 jan 2009     **/
/**                                 to   : 10 sep 2011     **/
/**                                                        **/
/************************************************************/

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_fold_comm.h"

/* This routine computes an optimized communication
** scheme for folding the data of a distributed graph.
** It is currently based on a maximum fixed number of
** communications per process. If this maximum is reached,
** the algorithm will fail.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphFoldComm (
const Dgraph * restrict const                 grafptr,
const int                                     partval, /* 0 for first half, 1 for second half                        */
int * restrict const                          commptr, /* Pointer to maximum number of communications per process    */
int * restrict const                          commtypval, /* Process will be sender or receiver                      */
DgraphFoldCommData *restrict * restrict const commdatptr, /* Slots for communication                                 */
Gnum *restrict * restrict const               commvrtptr, /* Slots of starting global vertex send indices            */
Gnum * restrict const                         proccnttab, /* Receive count array, for receivers                      */
int * restrict const                          vertadjnbrptr, /* Number of adjustment ranges, for receivers           */
Gnum * restrict * restrict const              vertadjptr, /* Pointer to global index adjustment array, for receivers */
Gnum * restrict * restrict const              vertdltptr) /* Pointer to global delta adjustment array, for receivers */
{
  int                           commmax;          /* Maximum number of communications per process        */
  int                           procnum;
  Gnum * restrict               commvrttab;
  DgraphFoldCommData * restrict commdattab;
  DgraphFoldCommData * restrict procsrttab;       /* Sort array                                          */
  int                           procsndbas;
  int                           procsndmnd;
  int                           procsndidx;
  int                           procrcvbas;
  int                           procrcvnnd;
  int                           procrcvidx;
  int                           fldprocnbr;
  Gnum * restrict               vertadjtab;
  Gnum * restrict               vertdlttab;
  int * restrict                vertprmtab;       /* Permutation array for adjustment range computations */
  int                           i;
#ifdef SCOTCH_DEBUG_DGRAPH2
  DgraphFoldCommData * restrict procchktab;
  int                           chekloctab[2];
  int                           chekglbtab[2];

  if (grafptr->procglbnbr < 2) {
    errorPrint ("dgraphFoldComm: invalid parameters");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if ((procsrttab = memAlloc (grafptr->procglbnbr * sizeof (DgraphFoldCommData))) == NULL) {
    errorPrint ("dgraphFoldComm: out of memory (1)");
    return     (1);
  }

  for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
    procsrttab[procnum].procnum = procnum;
    procsrttab[procnum].vertnbr = grafptr->proccnttab[procnum];
  }
  fldprocnbr = (grafptr->procglbnbr + 1) / 2;     /* Get number of processes in part 0 (always more than in part 1) */

  intSort2asc1 (procsrttab, fldprocnbr);          /* Sort both parts of processor array */
  intSort2asc1 (procsrttab + fldprocnbr, grafptr->procglbnbr - fldprocnbr);

  if (partval == 0) {                             /* If part 0 will receive the data                 */
    procrcvbas = 0;                               /* Receive by ascending weight order in first part */
    procrcvnnd = fldprocnbr;
    procsndbas = grafptr->procglbnbr;             /* Send by descending weight order in other part */
    procsndmnd = fldprocnbr;
  }
  else {                                          /* Part 1 will receive the data                    */
    procrcvbas = fldprocnbr;                      /* Receive by ascending weight order in first part */
    procrcvnnd = grafptr->procglbnbr;
    procsndbas = fldprocnbr;                      /* Send by descending weight order in other part */
    procsndmnd = 0;
    fldprocnbr = grafptr->procglbnbr - fldprocnbr;
  }
  *commtypval = ((grafptr->proclocnum >= procrcvbas) && (grafptr->proclocnum < procrcvnnd)) ? DGRAPHFOLDCOMMRECV : DGRAPHFOLDCOMMSEND;
#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((*commtypval == DGRAPHFOLDCOMMRECV) &&
      ((proccnttab == NULL) || (vertadjptr == NULL) || (vertdltptr == NULL))) {
    errorPrint ("dgraphFoldComm: internal error (1)");
    memFree    (procsrttab);                      /* Free group leader */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (commmax = DGRAPHFOLDCOMMNBR, commdattab = NULL; ; commmax ++) { /* Start with small number of communications per process */
    int                   procrcvnum;
    int                   procsndnum;
    int                   procsndnxt;             /* Index of next sender to process after current sender                */
    Gnum                  vertglbsum;             /* Global number of vertices already sent from all previous senders    */
    Gnum                  vertglbavg;             /* Global average value of sent vertices to reach during current round */
    Gnum                  vertsndrmn;             /* Remaining number of vertices to send from current sending process   */
    int                   flagrcvval;             /* Number of messages already sent by receiver and sender              */
    int                   flagsndval;
    int                   commtmp;

    if (commdattab != NULL)
      memFree (commdattab);

    commtmp = (proccnttab == NULL) ? 0 : commmax;
    if (memAllocGroup ((void **) (void *)
                       &commdattab, (size_t) (commmax * sizeof (DgraphFoldCommData)),
                       &commvrttab, (size_t) (commmax * sizeof (Gnum)),
                       &vertadjtab, (size_t) (commtmp * grafptr->procglbnbr * sizeof (Gnum)),
                       &vertdlttab, (size_t) (commtmp * grafptr->procglbnbr * sizeof (Gnum)),
#ifdef SCOTCH_DEBUG_DGRAPH2
                       &procchktab, (size_t) (commmax * grafptr->procglbnbr * sizeof (DgraphFoldCommData)),
#endif /* SCOTCH_DEBUG_DGRAPH2 */
                       &vertprmtab, (size_t) ((commmax + 1) * fldprocnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("dgraphFoldComm: out of memory (2)");
      memFree    (procsrttab);
      return     (1);
    }

    for (i = 0; i < commmax; i ++) {              /* Assume we will perform no communication at all */
      commdattab[i].procnum = -1;
      commdattab[i].vertnbr = 0;
    }
    if (proccnttab != NULL) {                     /* If we are a receiver process, start filling count and adjustment arrays */
      int                   procrcvtmp;

      memSet (vertprmtab, ~0, (commmax + 1) * fldprocnbr * sizeof (int)); /* Reset adjustment index arrays */
      memSet (vertadjtab, ~0, commmax * grafptr->procglbnbr * sizeof (Gnum));

      for (procrcvtmp = procrcvbas; procrcvtmp < procrcvnnd; procrcvtmp ++) { /* Fill count and adjustment arrays for receiver processes slots */
        vertadjtab[procrcvtmp * commmax]    = grafptr->procvrttab[procrcvtmp];
        vertdlttab[procrcvtmp * commmax]    =
        proccnttab[procrcvtmp - procrcvbas] = grafptr->proccnttab[procrcvtmp];
        vertprmtab[(procrcvtmp - procrcvbas) * (commmax + 1)] = procrcvtmp * commmax;
      }
    }

    for (procrcvidx = procrcvbas; procrcvidx < (procrcvnnd - 1); procrcvidx ++, procrcvnnd --) { /* Overloaded receiver vertices will re-send some of their load */
      Gnum                  vertsndnbr;           /* Potential overload of receiver process */

      procsndidx = procrcvnnd - 1;
      vertsndnbr = procsrttab[procsndidx].vertnbr - DATASIZE (grafptr->vertglbnbr, fldprocnbr, procsndidx - procrcvbas);
      if (vertsndnbr <= 0)                        /* If no further overload to improve, exit loop */
        break;

      procrcvnum = procsrttab[procrcvidx].procnum;
      procsndnum = procsrttab[procsndidx].procnum;

      procsrttab[procrcvidx].vertnbr = - (procsrttab[procrcvidx].vertnbr + vertsndnbr); /* Flag as having had a communication */
      if (proccnttab != NULL) {                   /* If we are a receiver process, fill count and adjustment arrays           */
        proccnttab[procrcvnum - procrcvbas] += vertsndnbr;
        proccnttab[procsndnum - procrcvbas] -= vertsndnbr; /* Vertices are transferred between receiver processes */
        vertadjtab[procsndnum * commmax + 1] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndnbr;
        vertdlttab[procsndnum * commmax + 1] = vertsndnbr;
        vertdlttab[procsndnum * commmax]    -= vertsndnbr;
        vertprmtab[(procrcvnum - procrcvbas) * (commmax + 1) + 1] = procsndnum * commmax + 1;

        if (procsndnum == grafptr->proclocnum) {  /* If we are the sender receiver process */
          *commtypval = DGRAPHFOLDCOMMSEND | DGRAPHFOLDCOMMRECV; /* Indicate it            */
          commdattab[0].procnum = procrcvnum;     /* Record communication                  */
          commdattab[0].vertnbr = vertsndnbr;
          commvrttab[0] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndnbr; /* Set starting global index of vertices to be sent */
        }
        else if (procrcvnum == grafptr->proclocnum) { /* If we are the receiver receiver process */
          commdattab[0].procnum = procsndnum;     /* Record communication                        */
          commdattab[0].vertnbr = vertsndnbr;
          commvrttab[0] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndnbr; /* Set starting global index of vertices to be sent */
        }
      }
    }

    while ((procsndmnd < procsndbas) && (procsrttab[procsndmnd].vertnbr == 0)) /* Do not account for empty sender processes */
      procsndmnd ++;

    if (procsndmnd >= procsndbas)                 /* If there are no more sender processes       */
      break;                                      /* We have completed communication computation */

    procsndidx = procsndbas - 1;
    procsndnxt = procsndidx - 1;
    vertsndrmn = procsrttab[procsndidx].vertnbr;
    flagsndval = 0;
    procrcvidx = procrcvbas;
    procrcvnum = procsrttab[procrcvidx].procnum;
    flagrcvval = 0;
    vertglbavg = DATASIZE (grafptr->vertglbnbr, fldprocnbr, 0);
    vertglbsum = procsrttab[procrcvidx].vertnbr;  /* Account for vertices already present in receiver process */
    if (vertglbsum < 0) {                         /* If first receiver has already received a communication   */
      flagrcvval = 1;
      vertglbsum = - vertglbsum;
      procsrttab[procrcvidx].vertnbr = vertglbsum; /* Un-flag */
    }

    while (1) {
      Gnum                vertrcvrmn;             /* Remaining number of vertices to receive           */
      Gnum                vertsndnbr;             /* Number of vertices actually sent and received     */
      Gnum                vertsndnxt;             /* Adjustment to add to vertex count of next sender  */
      int                 mesgrcvrmn;             /* Number of message slots to receive small messages */
      int                 flagsndnxt;

#ifdef SCOTCH_DEBUG_DGRAPH2
      if (flagrcvval + flagsndval >= (commmax * 2)) {
        errorPrint ("dgraphFoldComm: internal error (2)");
        memFree    (commdattab);                  /* Free group leader */
        memFree    (procsrttab);
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      vertrcvrmn = (vertglbavg > vertglbsum) ? (vertglbavg - vertglbsum) : 0; /* Remaining space on receiver */
      mesgrcvrmn = (vertsndrmn >= vertrcvrmn) ? (commmax - 1) : (commmax - 2);
      if ((procsndidx > procsndmnd) &&            /* If there remains small messages to be considered          */
          (((flagsndval == 0) &&                  /* If sender has not sent anything yet or just started       */
            (flagrcvval < mesgrcvrmn)) ||         /* And receiver has space for larger messages to come        */
           ((flagrcvval == 0) &&                  /* Or if receiver has not received anything yet              */
            (flagsndval < (commmax - 2)))) &&     /* And sender has enough slots to send                       */
          (vertrcvrmn >= procsrttab[procsndmnd].vertnbr)) { /* And if receiver can hold small message entirely */
        procsndnxt = procsndidx;                  /* Current large message will be processed next time         */
        procsndidx = procsndmnd;                  /* Process smallest message available to date                */
        flagsndnxt = flagsndval;                  /* Record location of next message to send                   */
        flagsndval = 0;
        vertsndnxt = procsrttab[procsndnxt].vertnbr - vertsndrmn; /* Record vertices already sent  */
        vertsndnbr =                              /* All of its contents will be sent in one piece */
        vertsndrmn = procsrttab[procsndmnd].vertnbr;
        procsndmnd ++;                            /* Small message has been processed */
      }
      else {
        flagsndnxt = 0;
        vertsndnxt = 0;                           /* Next sender will not have been interrupted by a small message */
        vertsndnbr = ((flagsndval >= (commmax - 1)) || /* If last chance to send for this process                  */
                      (((procrcvnnd - procrcvidx) * (commmax - 1) - flagrcvval) <= (procsndidx - procsndmnd))) /* Or if too few communications remain with sender receivers accounted for */
                     ? vertsndrmn                 /* Send all of the vertices to be sent       */
                     : MIN (vertsndrmn, vertrcvrmn); /* Else just send what the receiver needs */
      }
      procsndnum = procsrttab[procsndidx].procnum;

      if (vertsndnbr > 0) {                       /* If useful communication can take place                         */
        if (proccnttab != NULL) {                 /* If we are a receiver process, fill count and adjustment arrays */
          proccnttab[procrcvnum - procrcvbas] += vertsndnbr;
          vertadjtab[procsndnum * commmax + flagsndval] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndrmn;
          vertdlttab[procsndnum * commmax + flagsndval] = vertsndnbr;
          vertprmtab[(procrcvnum - procrcvbas) * (commmax + 1) + 1 + flagrcvval] = procsndnum * commmax + flagsndval;

          if (procrcvnum == grafptr->proclocnum) { /* If we are the receiver process */
            commdattab[flagrcvval].procnum = procsndnum; /* Record communication     */
            commdattab[flagrcvval].vertnbr = vertsndnbr;
            commvrttab[flagrcvval] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndrmn; /* Set starting global index of vertices to be sent */
          }
        }
        else if (procsndnum == grafptr->proclocnum) { /* If we are the sending process */
          commdattab[flagsndval].procnum = procrcvnum; /* Record communication         */
          commdattab[flagsndval].vertnbr = vertsndnbr;
          commvrttab[flagsndval] = grafptr->procvrttab[procsndnum] + grafptr->proccnttab[procsndnum] - vertsndrmn; /* Set starting global index of vertices to be sent */
        }
        vertglbsum += vertsndnbr;                 /* Account for vertices sent and received */
        vertsndrmn -= vertsndnbr;
        flagsndval ++;
        flagrcvval ++;
      }
      if (vertsndrmn <= 0) {                      /* If sending process has sent everything */
        procsndidx = procsndnxt;                  /* Process next vertex to send            */
        if (procsndidx < procsndmnd)              /* If was last sending process, end loop  */
          break;
        procsndnxt = procsndidx - 1;              /* Prepare next sender process  */
        flagsndval = flagsndnxt;                  /* Skip to next sending process */
        vertsndrmn = procsrttab[procsndidx].vertnbr - vertsndnxt;
      }
      if ((flagrcvval >= commmax) ||              /* If receiver cannot receive more                                 */
          ((vertglbsum >= vertglbavg) &&          /* Or has received what it needed and is not forced to accept more */
           (((procrcvnnd - procrcvidx) * (commmax - 1) - flagrcvval) > (procsndidx - procsndmnd)))) {
        if (++ procrcvidx >= procrcvnnd)          /* If was last receiver, exit loop and go finalizing communication arrays */
          break;

        procrcvnum  = procsrttab[procrcvidx].procnum; /* Skip to next receiver process */
        vertglbavg += DATASIZE (grafptr->vertglbnbr, fldprocnbr, procrcvidx - procrcvbas);
        if (procsrttab[procrcvidx].vertnbr >= 0) { /* If receiver did not receive from a sender receiver */
          flagrcvval  = 0;
          vertglbsum += procsrttab[procrcvidx].vertnbr; /* Account for vertices already present in receiver process */
        }
        else {                                    /* Receiver already received from a sender receiver */
          flagrcvval  = 1;                        /* Already a communication performed                */
          vertglbsum -= procsrttab[procrcvidx].vertnbr; /* Use negative value                         */
        }
      }
    }

    if ((procsndidx <= procsndmnd) && (vertsndrmn <= 0)) /* If no sender vertex which has something to send remains */
      break;                                      /* Exit the loop on increasing number of communications           */
  }

#ifdef SCOTCH_DEBUG_DGRAPH2
  chekloctab[0] = - commmax;
  chekloctab[1] =   commmax;
  if (MPI_Allreduce (chekloctab, chekglbtab, 2, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFoldComm: communication error");
    memFree    (commdattab);                      /* Free group leader */
    return     (1);
  }
  if ((chekglbtab[0] != chekloctab[0]) ||
      (chekglbtab[1] != chekloctab[1])) {
    errorPrint ("dgraphFoldComm: internal error (3)");
    memFree    (commdattab);                      /* Free group leader */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (proccnttab != NULL) {                       /* If we are a receiver process                  */
    int                   vertadjnbr;             /* Number of adjustment slots                    */
    Gnum                  vertadjsum;             /* Current new starting position of current slot */

    for (procnum = 0, vertadjsum = grafptr->baseval; procnum < ((commmax + 1) * fldprocnbr); procnum ++) {
      int                   vertprmnum;
      Gnum                  vertadjtmp;

      vertprmnum = vertprmtab[procnum];
      if (vertprmnum == ~0)                       /* Skip empty slots */
        continue;

      vertadjtmp = vertdlttab[vertprmnum];        /* Accumulate new slot indices and compute adjustments */
      vertdlttab[vertprmnum] = vertadjsum - vertadjtab[vertprmnum];
      vertadjsum += vertadjtmp;
    }

    for (procnum = vertadjnbr = 0; procnum < (commmax * grafptr->procglbnbr); procnum ++) { /* Compact vertex adjustment arrays */
      if (vertadjtab[procnum] != -1) {
        vertadjtab[vertadjnbr] = vertadjtab[procnum];
        vertdlttab[vertadjnbr] = vertdlttab[procnum];
        vertadjnbr ++;
      }
    }
    vertadjtab[vertadjnbr] = grafptr->procvrttab[grafptr->procglbnbr]; /* Set upper bound on global vertex indices */

    *vertadjnbrptr = vertadjnbr;
    *vertadjptr    = vertadjtab;
    *vertdltptr    = vertdlttab;
  }
  *commdatptr = commdattab;                       /* Set group leader */
  *commvrtptr = commvrttab;
  *commptr    = commmax;

  memFree (procsrttab);

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (proccnttab == NULL) {                       /* If we are a sender process       */
    Gnum                  vertsndnbr;             /* Number of vertices actually sent */
    int                   i;

    for (i = 0, vertsndnbr = 0; (i < commmax) && (commdattab[i].procnum >= 0); i ++)
      vertsndnbr += commdattab[i].vertnbr;
    if (vertsndnbr != grafptr->vertlocnbr) {
      errorPrint ("dgraphFoldComm: internal error (4)");
      memFree    (commdattab);                    /* Free group leader */
      return     (1);
    }
  }
  if (MPI_Allgather (commdattab, 2 * commmax, GNUM_MPI,
                     procchktab, 2 * commmax, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFoldComm: communication error");
    memFree    (commdattab);                      /* Free group leader */
    return     (1);
  }

  for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
    int                   commnum;

    for (commnum = 0; (commnum < commmax) && (procchktab[commmax * procnum + commnum].procnum != -1); commnum ++) {
      Gnum                  procend;
      Gnum                  vertnbr;
      int                   commend;

      procend = procchktab[commmax * procnum + commnum].procnum;
      vertnbr = procchktab[commmax * procnum + commnum].vertnbr;

      for (commend = 0; commend < commmax; commend ++) {
        if ((procchktab[commmax * procend + commend].procnum == procnum) &&
            (procchktab[commmax * procend + commend].vertnbr == vertnbr))
          break;
      }
      if (commend >= commmax) {
        errorPrint ("dgraphFoldComm: internal error (5)");
        memFree    (commdattab);                  /* Free group leader */
        return     (1);
      }
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}
