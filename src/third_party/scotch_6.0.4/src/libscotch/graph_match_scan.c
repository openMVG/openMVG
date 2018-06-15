/* Copyright 2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_match_scan.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the stub of the    **/
/**                threaded and un-threaded centralized    **/
/**                graph matching functions.               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 01 oct 2012     **/
/**                                 to     04 jun 2014     **/
/**                                                        **/
/**   NOTES      : # This code partly derives from the     **/
/**                  code of graph_match.c partly updated  **/
/**                  in the early stages of release 6.0    **/
/**                  so as to account for fixed vertices   **/
/**                  and an old partition.                 **/
/**                                                        **/
/************************************************************/

/***************************/
/*                         */
/* The matching subroutine */
/* pattern.                */
/*                         */
/***************************/

/* This routine matches the vertices of the given
** centralized graph according to various constraints.
** It returns:
** - void  : in all cases
*/

static
void
GRAPHMATCHSCANNAME (
GraphCoarsenThread * restrict thrdptr)            /* Thread-dependent data */
{
  GraphCoarsenData * restrict const coarptr = (GraphCoarsenData *) (thrdptr->thrddat.grouptr);
  const Graph * restrict const      finegrafptr = coarptr->finegrafptr;
  Gnum                              finevertnum;  /* Current vertex index        */
  Gnum                              finevertnnd;  /* Current end of vertex array */
#ifdef GRAPHMATCHSCANVELOTAB
  Gnum                    finevelomin = (1 * finegrafptr->velosum) / (4 * finegrafptr->vertnbr); /* Minimum load of neighbor */
  Gnum                    coarvelomax = (4 * finegrafptr->velosum) / (1 * (coarptr->coarvertmax - coarptr->finevfixnbr)) + 1;
#endif /* GRAPHMATCHSCANVELOTAB */
  Gnum                    coarvertnbr = thrdptr->coarvertnbr; /* Current number of multinode vertices */
#if ((defined GRAPHMATCHSCANP1INPERT) || defined (GRAPHMATCHSCANP2INPERT))
  const Gnum              degrmax = finegrafptr->degrmax;
  Gnum                    pertbas;
  Gnum                    pertnnd;
  Gnum                    pertnbr;
  Gunum                   randval = thrdptr->randval; /* Unsigned to avoid negative random numbers */
#endif /* ((defined GRAPHMATCHSCANP1INPERT) || defined (GRAPHMATCHSCANP2INPERT)) */
#if ((defined GRAPHMATCHSCANP1INQUEUE) || defined (GRAPHMATCHSCANP2INQUEUE) || defined (GRAPHMATCHSCANP2OUTQUEUE))
  Gnum * const            queutab = coarptr->finequeutab;
  volatile int * const    locktax = coarptr->finelocktax;
#endif /* ((defined GRAPHMATCHSCANP1INQUEUE) || defined (GRAPHMATCHSCANP2INQUEUE) || defined (GRAPHMATCHSCANP2OUTQUEUE)) */
#if ((defined GRAPHMATCHSCANP1INQUEUE) || defined (GRAPHMATCHSCANP2INQUEUE))
  Gnum                    queunnd;
  Gnum                    queunum;
#endif /* ((defined GRAPHMATCHSCANP1INQUEUE) || defined (GRAPHMATCHSCANP2INQUEUE)) */
#if ((defined GRAPHMATCHSCANP1OUTQUEUE) || defined (GRAPHMATCHSCANP2OUTQUEUE))
  Gnum                    queunew;
#endif /* ((defined GRAPHMATCHSCANP1OUTQUEUE) || defined (GRAPHMATCHSCANP2OUTQUEUE)) */

  const Gnum * restrict const     fineverttax = finegrafptr->verttax;
  const Gnum * restrict const     finevendtax = finegrafptr->vendtax;
#ifdef GRAPHMATCHSCANVELOTAB
  const Gnum * restrict const     finevelotax = finegrafptr->velotax;
#endif /* GRAPHMATCHSCANVELOTAB */
  const Gnum * restrict const     fineedgetax = finegrafptr->edgetax;
#ifdef GRAPHMATCHSCANEDLOTAB
  const Gnum * restrict const     fineedlotax = finegrafptr->edlotax;
#endif /* GRAPHMATCHSCANEDLOTAB */
  volatile Gnum * restrict const  finematetax = coarptr->finematetax;
#ifdef GRAPHMATCHSCANPFIXTAB
  const Gnum * restrict const     fineparotax = coarptr->fineparotax;
  const Gnum * restrict const     finepfixtax = coarptr->finepfixtax;
#endif /* GRAPHMATCHSCANPFIXTAB */

#ifdef GRAPHMATCHSCANP1INPERT                     /* First pass is only for start or sequential routines */
#ifdef GRAPHMATCHSCANVELOTAB
#ifdef GRAPHMATCHSCANP1OUTQUEUE
  queunew = thrdptr->queubas;
#endif /* GRAPHMATCHSCANP1OUTQUEUE */
#ifdef GRAPHMATCHSCANP1INQUEUE
  for (queunum = thrdptr->queubas, queunnd = thrdptr->queunnd;
       queunum < queunnd; queunum ++) {
    finevertnum = queutab[queunum];
#endif /* GRAPHMATCHSCANP1INQUEUE */
#ifdef GRAPHMATCHSCANP1INPERT
  for (pertbas = thrdptr->finequeubas, pertnnd = thrdptr->finequeunnd;
       pertbas < pertnnd; pertbas += pertnbr) {   /* Run cache-friendly perturbation    */
    Gnum                pertval;                  /* Current index in perturbation area */

    pertnbr = degrmax * 2 + randval % (degrmax + 1) + 1; /* Compute perturbation area size (avoid DIV0 in random) */
    if (pertnbr >= GRAPHMATCHSCANPERTPRIME)
      pertnbr = 32 + randval % (GRAPHMATCHSCANPERTPRIME - 34);
    if (pertbas + pertnbr > pertnnd)
      pertnbr = pertnnd - pertbas;

    pertval = 0;                                  /* Start from first perturbation vertex */
    do {                                          /* Loop on perturbation vertices        */
      finevertnum = pertbas + pertval;            /* Compute corresponding vertex number  */
#endif /* GRAPHMATCHSCANP1INPERT */
  {
    Gnum                fineedgenum;
    Gnum                finevertbst;
#ifdef GRAPHMATCHSCANEDLOTAB
    Gnum                fineedlobst = -1;         /* Edge load of current best neighbor */
#endif /* GRAPHMATCHSCANEDLOTAB */

    if (finematetax[finevertnum] >= 0)            /* If vertex already mated, skip it without remembering it */
      goto loop1;

    finevertbst = finevertnum;                    /* Assume we match with ourselves */

    if ((finevelotax[finevertnum] >= finevelomin) || /* If vertex is not too light, skip it and record it  */
        (fineverttax[finevertnum] == finevendtax[finevertnum])) /* Same for isolated vertex: process later */
#ifdef GRAPHMATCHSCANP1OUTQUEUE
      goto record1;
#else /* GRAPHMATCHSCANP1OUTQUEUE */
      goto loop1;
#endif /* GRAPHMATCHSCANP1OUTQUEUE */

    for (fineedgenum = fineverttax[finevertnum];
         fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
      Gnum                finevertend;

      finevertend = fineedgetax[fineedgenum];

      if ((finematetax[finevertend] < 0)          /* If unmatched vertex */
#ifdef GRAPHMATCHSCANPFIXTAB
          && ((finepfixtax == NULL) || (finepfixtax[finevertend] == finepfixtax[finevertnum])) /* We can only mate if potential mate has same value */
          && ((fineparotax == NULL) || (fineparotax[finevertend] == fineparotax[finevertnum])) /* And is in the same old part                       */
#endif /* GRAPHMATCHSCANPFIXTAB */
#ifdef GRAPHMATCHSCANEDLOTAB
          && (fineedlotax[fineedgenum] > fineedlobst) /* And is better candidate */
#endif /* GRAPHMATCHSCANEDLOTAB */
      ) {
        finevertbst = finevertend;
#ifdef GRAPHMATCHSCANEDLOTAB
        fineedlobst = fineedlotax[fineedgenum];
#else /* GRAPHMATCHSCANEDLOTAB */
        break;
#endif /* GRAPHMATCHSCANEDLOTAB */
      }
    }

#ifdef GRAPHMATCHSCANP2OUTQUEUE
    if (__sync_lock_test_and_set (&locktax[finevertnum], 1)) /* If could not acquire local vertex                      */
      goto loop1;                                 /* Do not remember it as some other vertex has already acquired both */

    if (finevertbst != finevertnum) {             /* If we mated with another vertex                 */
      if (__sync_lock_test_and_set (&locktax[finevertbst], 1)) { /* If could not acquire mate vertex */
        __sync_lock_release (&locktax[finevertnum]); /* Release lock on local vertex                 */
#ifdef GRAPHMATCHSCANP1OUTQUEUE
record1:
        queutab[queunew ++] = finevertnum;        /* Postpone processing to next pass */
#endif /* GRAPHMATCHSCANP1OUTQUEUE */
        goto loop1;
      }
      finematetax[finevertbst] = finevertnum;
    }
#else /* GRAPHMATCHSCANP2OUTQUEUE */
    finematetax[finevertbst] = finevertnum;       /* If last (sequential) pass, always record */
#endif /* GRAPHMATCHSCANP2OUTQUEUE */
    finematetax[finevertnum] = finevertbst;       /* At this point or if last pass, record */
    coarvertnbr ++;                               /* One more coarse vertex created        */
loop1: ;
  }
#ifdef GRAPHMATCHSCANP1INPERT
      pertval = (pertval + GRAPHMATCHSCANPERTPRIME) % pertnbr; /* Compute next perturbation index */
    } while (pertval != 0);
    randval += finevertnum;                       /* Perturbation for thread-dependent pseudo-random number */
  }
#endif /* GRAPHMATCHSCANP1INPERT */
#ifdef GRAPHMATCHSCANP1INQUEUE
  }
#endif /* GRAPHMATCHSCANP1INQUEUE */
#ifdef GRAPHMATCHSCANP1OUTQUEUE
  thrdptr->queunnd = queunew;                     /* Record queue index for next pass */
#endif /* GRAPHMATCHSCANP1OUTQUEUE */
#endif /* GRAPHMATCHSCANVELOTAB    */
#endif /* GRAPHMATCHSCANP1INPERT   */

#ifdef GRAPHMATCHSCANP2OUTQUEUE
  queunew = thrdptr->finequeubas;
#endif /* GRAPHMATCHSCANP2OUTQUEUE */
#ifdef GRAPHMATCHSCANP2INQUEUE
  for (queunum = thrdptr->finequeubas, queunnd = thrdptr->finequeunnd;
       queunum < queunnd; queunum ++) {
    finevertnum = queutab[queunum];
#endif /* GRAPHMATCHSCANP2INQUEUE */
#ifdef GRAPHMATCHSCANP2INPERT
  for (pertbas = thrdptr->finequeubas, pertnnd = thrdptr->finequeunnd;
       pertbas < pertnnd; pertbas += pertnbr) {   /* Run cache-friendly perturbation    */
    Gnum                pertval;                  /* Current index in perturbation area */

    pertnbr = degrmax * 2 + randval % (degrmax + 1) + 1; /* Compute perturbation area size (avoid DIV0 in random) */
    if (pertnbr >= GRAPHMATCHSCANPERTPRIME)
      pertnbr = 32 + randval % (GRAPHMATCHSCANPERTPRIME - 34);
    if (pertbas + pertnbr > pertnnd)
      pertnbr = pertnnd - pertbas;

    pertval = 0;                                  /* Start from first perturbation vertex */
    do {                                          /* Loop on perturbation vertices        */
      finevertnum = pertbas + pertval;            /* Compute corresponding vertex number  */
#endif /* GRAPHMATCHSCANP2INPERT */
  {
    Gnum                fineedgenum;
    Gnum                finevertbst;

    if (finematetax[finevertnum] >= 0)            /* If vertex already mated, skip it without remembering it */
      goto loop2;

    if (fineverttax[finevertnum] == finevendtax[finevertnum]) { /* If isolated vertex */
#ifdef GRAPHMATCHSCANP2INPERT
      Gnum                perttmp = pertnnd;
      do
        finevertbst = -- perttmp;
#endif /* GRAPHMATCHSCANP2INPERT */
#ifdef GRAPHMATCHSCANP2INQUEUE
      Gnum                queutmp = queunnd;
      do
        finevertbst = queutab[-- queutmp];
#endif /* GRAPHMATCHSCANP2INQUEUE */
      while ((finematetax[finevertbst] >= 0)      /* No test for overflow; we will always mate ourselves as we are isolated */
#ifdef GRAPHMATCHSCANPFIXTAB
             || ((finepfixtax != NULL) && (finepfixtax[finevertbst] != fineparotax[finevertnum]))
             || ((fineparotax != NULL) && (fineparotax[finevertbst] != fineparotax[finevertnum])));
#else /* GRAPHMATCHSCANPFIXTAB */
      );
#ifdef GRAPHMATCHSCANP2INPERT
      pertnnd = perttmp;                          /* If no extra conditions on mating, no longer consider traversed array */
#endif /* GRAPHMATCHSCANP2INPERT */
#ifdef GRAPHMATCHSCANP2INQUEUE
      queunnd = queutmp;
#endif /* GRAPHMATCHSCANP2INQUEUE */
#endif /* GRAPHMATCHSCANPFIXTAB   */
    }
    else {
      Gnum                fineedgenum;
#ifdef GRAPHMATCHSCANVELOTAB
      Gnum                finevelodlt = coarvelomax - finevelotax[finevertnum];
#endif /* GRAPHMATCHSCANVELOTAB */
#ifdef GRAPHMATCHSCANEDLOTAB
      Gnum                fineedlobst = -1;
#endif /* GRAPHMATCHSCANEDLOTAB */

      finevertbst = finevertnum;                  /* No matching neighbor found yet */
      for (fineedgenum = fineverttax[finevertnum]; /* For all adjacent vertices     */
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        Gnum                finevertend;

        finevertend = fineedgetax[fineedgenum];
        if ((finematetax[finevertend] < 0)        /* If unmatched vertex */
#ifdef GRAPHMATCHSCANPFIXTAB
            && ((finepfixtax == NULL) || (finepfixtax[finevertend] == finepfixtax[finevertnum])) /* And is in the same part */
            && ((fineparotax == NULL) || (fineparotax[finevertend] == fineparotax[finevertnum]))
#endif /* GRAPHMATCHSCANPFIXTAB */
#ifdef GRAPHMATCHSCANVELOTAB
            && (finevelodlt >= finevelotax[finevertend]) /* And does not create overloads */
#endif /* GRAPHMATCHSCANVELOTAB */
#ifdef GRAPHMATCHSCANEDLOTAB
            && (fineedlotax[fineedgenum] > fineedlobst)
#endif /* GRAPHMATCHSCANEDLOTAB */
           ) { /* And is better candidate */
          finevertbst = finevertend;
#ifdef GRAPHMATCHSCANEDLOTAB
          fineedlobst = fineedlotax[fineedgenum];
#else /* GRAPHMATCHSCANEDLOTAB */
          break;                                  /* First match will do if no edge loads */
#endif /* GRAPHMATCHSCANEDLOTAB */
        }
      }
    }

#ifdef GRAPHMATCHSCANP2OUTQUEUE
    if (__sync_lock_test_and_set (&locktax[finevertnum], 1)) /* If could not acquire local vertex                      */
      goto loop2;                                 /* Do not remember it as some other vertex has already acquired both */

    if (finevertbst != finevertnum) {             /* If we mated with another vertex                 */
      if (__sync_lock_test_and_set (&locktax[finevertbst], 1)) { /* If could not acquire mate vertex */
        __sync_lock_release (&locktax[finevertnum]); /* Release lock on local vertex                 */
        queutab[queunew ++] = finevertnum;        /* Postpone processing to next pass                */
        goto loop2;
      }
      finematetax[finevertbst] = finevertnum;
    }
#else /* GRAPHMATCHSCANP2OUTQUEUE */
    finematetax[finevertbst] = finevertnum;       /* If last (sequential) pass, always record */
#endif /* GRAPHMATCHSCANP2OUTQUEUE */
    finematetax[finevertnum] = finevertbst;       /* At this point or if last pass, record */
    coarvertnbr ++;                               /* One more coarse vertex created        */
loop2: ;                                          /* We may do something before looping    */
  }
#ifdef GRAPHMATCHSCANP2INPERT
      pertval = (pertval + GRAPHMATCHSCANPERTPRIME) % pertnbr; /* Compute next perturbation index */
    } while (pertval != 0);
    randval += finevertnum;                       /* Perturbation for thread-dependent pseudo-random number */
  }
#endif /* GRAPHMATCHSCANP2INPERT */
#ifdef GRAPHMATCHSCANP2INQUEUE
  }
#endif /* GRAPHMATCHSCANP2INQUEUE */
#ifdef GRAPHMATCHSCANP2OUTQUEUE
  thrdptr->finequeunnd = queunew;                 /* Record queue index for next pass */
#endif /* GRAPHMATCHSCANP2OUTQUEUE */
  thrdptr->coarvertnbr = coarvertnbr;             /* Record subsequent number of multinode vertices */
}
