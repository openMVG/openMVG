/* Copyright 2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_ex.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module bipartitions a distributed  **/
/**                active graph using a parallel gradient  **/
/**                method to balance the two parts as      **/
/**                evenly as possible.                     **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 jul 2010     **/
/**                                 to   : 15 apr 2011     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to   : 31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_BIPART_EX

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "arch.h"
#include "bdgraph.h"
#include "bdgraph_bipart_ex.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bdgraphBipartEx (
Bdgraph * restrict const                    grafptr, /*+ Active graph      +*/
const BdgraphBipartExParam * restrict const paraptr) /*+ Method parameters +*/
{
  int                   sbbtnbr;                  /* Number of subbits               */
  Gnum                  sbbtmsk;                  /* Subbit mask                     */
  int                   gainsiz;                  /* Size of gain array              */
  int                   gainidx;
  Gnum * restrict       gainloctab;
  Gnum * restrict       gainglbtab;
  Gnum * restrict       frstloctab;
  Gnum * restrict       nextloctab;
  Gnum * restrict       loadglbtab;
  const Gnum * restrict edgegsttax;
  Gnum                  edlolocval;
  Gnum                  fronlocnum;
  Gnum                  vertgstnum;
  Gnum                  vertlocnnd;
  Gnum                  complocsizedlt;           /* Count of vertices moved locally */
  Gnum                  complocloaddlt;           /* Load of vertices moved locally  */
  Gnum                  compglbloaddlt;
  Gnum                  compglbloaddltmax;
  Gnum                  compglbloaddltmat;
  Gnum                  commlocgain;
  Gnum                  commlocgainextn;
  Gnum                  fronlocnbr;
  int * restrict        movegsttab;
  int * restrict        flagloctab;
  int                   cheklocval;
  int                   chekglbval;
  BdgraphBipartExSort * sorttab;
  size_t                sortsiz;
  Gnum                  reduloctab[5];
  Gnum                  reduglbtab[5];
  Gnum                  domndist;
  Gnum                  partval;

  const Gnum * restrict const vertloctax = grafptr->s.vertloctax; /* Fast accesses */
  const Gnum * restrict const vendloctax = grafptr->s.vendloctax;
  const Gnum * restrict const veloloctax = grafptr->s.veloloctax;
  const Gnum * restrict const veexloctax = grafptr->veexloctax;
  const Gnum * restrict const edloloctax = grafptr->s.edloloctax;
  Gnum * restrict const       fronloctab = grafptr->fronloctab;
  GraphPart * restrict const  partgsttax = grafptr->partgsttax;

  partval = (grafptr->compglbload0dlt > 0) ? 1 : 0; /* Get number of underloaded part to receive vertices */

  compglbloaddltmax = (Gnum) ((double) grafptr->compglbload0avg * paraptr->deltval);
  compglbloaddltmat = (partval == 0)
                      ? (grafptr->compglbload0avg - grafptr->compglbload0min)
                      : (grafptr->compglbload0max - grafptr->compglbload0avg);
  if (compglbloaddltmax > compglbloaddltmat)
    compglbloaddltmax = compglbloaddltmat;

  if ((abs (grafptr->compglbload0dlt) < compglbloaddltmax) || /* If nothing to do */
      (grafptr->fronglbnbr == 0))                 /* Or if no current frontier    */
    return (0);                                   /* This algorithm is useless    */

  if (dgraphGhst (&grafptr->s) != 0) {            /* Compute ghost edge array if not already present */
    errorPrint ("bdgraphBipartEx: cannot compute ghost edge array");
    return     (1);
  }

  sbbtnbr = (int) paraptr->sbbtnbr;
  sbbtmsk = (1 << sbbtnbr) - 1;
  gainsiz = ((sizeof (Gnum) << 3) - sbbtnbr) << (sbbtnbr + 1); /* Compute gain array size  */
  sortsiz = MAX ((grafptr->s.procglbnbr * sizeof (BdgraphBipartExSort)), /* TRICK: recycle */
                 (grafptr->fronlocnbr   * sizeof (BdgraphBipartExMove)));

  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)
                     &gainglbtab, (size_t) (gainsiz * sizeof (Gnum)),
                     &gainloctab, (size_t) (gainsiz * sizeof (Gnum)),
                     &frstloctab, (size_t) (gainsiz * sizeof (Gnum)),
                     &nextloctab, (size_t) (grafptr->fronlocnbr   * sizeof (Gnum)),
                     &loadglbtab, (size_t) (grafptr->s.procglbnbr * sizeof (Gnum)),
                     &movegsttab, (size_t) (flagSize (grafptr->s.vertgstnbr + grafptr->s.baseval) * sizeof (int)), /* TRICK: use based vertices as flag array indices */
                     &flagloctab, (size_t) (flagSize (grafptr->s.vertlocnbr + grafptr->s.baseval) * sizeof (int)),
                     &sorttab,    (size_t) (sortsiz), NULL) == NULL) {
    errorPrint ("bdgraphBipartEx: out of memory");
    cheklocval = 1;
  }
  else {
    memSet (gainloctab,  0, gainsiz * sizeof (Gnum)); /* Initialize gain array  */
    memSet (frstloctab, ~0, gainsiz * sizeof (Gnum)); /* Initialize linked list */
    memSet (nextloctab, ~0, grafptr->fronlocnbr * sizeof (Gnum));
    memSet (movegsttab,  0, flagSize (grafptr->s.vertgstnbr + grafptr->s.baseval) * sizeof (int)); /* TRICK: based sizes */
    memSet (flagloctab,  0, flagSize (grafptr->s.vertlocnbr + grafptr->s.baseval) * sizeof (int));
  }

#ifdef SCOTCH_DEBUG_BDGRAPH1
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartEx: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  if (chekglbval != 0) {
    if (gainglbtab != NULL)
      memFree (gainglbtab);                       /* Free group leader */
    return (1);
  }

  domndist = (Gnum) grafptr->domndist;

  edgegsttax = grafptr->s.edgegsttax;
  edlolocval = 1;                                 /* Assume no edge loads */
  for (fronlocnum = 0; fronlocnum < grafptr->fronlocnbr; fronlocnum ++) {
    Gnum                vertlocnum;
    Gnum                velolocval;
    Gnum                edgelocnum;
    Gnum                commgain;
    int                 gainidx;
    int                 i;
    Gnum                j;

    vertlocnum = fronloctab[fronlocnum];

    if ((Gnum) partgsttax[vertlocnum] == partval) /* If vertex belongs to lighter part, skip it */
      continue;

    for (edgelocnum = vertloctax[vertlocnum], commgain = 0;
         edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
      Gnum                vertlocend;
      Gnum                partend;
      Gnum                partdlt;

      vertlocend = edgegsttax[edgelocnum];
      partend = (Gnum) partgsttax[vertlocend];
      if (edloloctax != NULL)
        edlolocval = edloloctax[edgelocnum];

      partdlt   = partval ^ partend;              /* Inverse of partdlt, because "partval" is the opposite */
      commgain += (2 * partdlt - 1) * edlolocval; /* Since partdlt has reversed meaning, reverse gain too  */
    }
    commgain *= domndist;                         /* Adjust internal gains with respect to external gains */
    if (veexloctax != NULL)
      commgain += (2 * partval - 1) * veexloctax[vertlocnum]; /* Partval has reversed meaning */

    velolocval = (veloloctax != NULL) ? veloloctax[vertlocnum] : 1;

    if (commgain >= 0) {                          /* Compute table entry for gain */
      for (i = 0, j = commgain; j > sbbtmsk; i ++, j >>= 1) ;
      i = (i << sbbtnbr) + (int) j;
    }
    else {
      for (i = 0, j = - (commgain + 1); j > sbbtmsk; i ++, j >>= 1) ;
      i = - ((i << sbbtnbr) + (int) j + 1);
    }

    gainidx = (gainsiz >> 1) + i;
#ifdef SCOTCH_DEBUG_BDGRAPH2
    if ((gainidx < 0) || (gainidx >= gainsiz)) {
      errorPrint ("bdgraphBipartEx: internal error");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    gainloctab[gainidx]   += velolocval;          /* Accumulate gain in proper cell     */
    nextloctab[fronlocnum] = frstloctab[gainidx]; /* Chain frontier vertex in gain list */
    frstloctab[gainidx]    = fronlocnum;
  }

  if (MPI_Allreduce (gainloctab, gainglbtab, gainsiz, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartEx: communication error (2)");
    return     (1);
  }

  complocloaddlt = 0;                             /* No moved vertices yet                              */
  compglbloaddlt = abs (grafptr->compglbload0dlt); /* We want to reduce the absolute value of imbalance */
  for (gainidx = 0; gainidx < gainsiz; gainidx ++) {
    Gnum                compglbloadtmp;
    Gnum                gainglbval;
    Gnum                fronlocnum;

    gainglbval = gainglbtab[gainidx];
    if (gainglbval <= 0)
      continue;

    compglbloadtmp = compglbloaddlt - gainglbval;
    if (compglbloadtmp < compglbloaddltmax)
      break;

    compglbloaddlt  = compglbloadtmp;
    complocloaddlt += gainloctab[gainidx];        /* We moved that much load locally */

    for (fronlocnum = frstloctab[gainidx]; fronlocnum != ~0; /* For all vertices in swapped gain slot */
         fronlocnum = nextloctab[fronlocnum])
      partgsttax[fronloctab[fronlocnum]] = (GraphPart) (partval | 2); /* Swap vertex part and flag vertex */
  }

  if ((gainidx < gainsiz) &&                      /* If we can make further adjustments */
      (compglbloaddlt > compglbloaddltmax)) {     /* And if there are some to make      */
    Gnum                loadglbmax;
    int                 procglbnbr;
    int                 proclocnum;
    int                 procnum;
    int                 procmax;
    int                 sortnbr;
    int                 sortnum;

    if (MPI_Allgather (gainloctab + gainidx, 1, GNUM_MPI,
                       loadglbtab, 1, GNUM_MPI, grafptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("bdgraphBipartEx: communication error (3)");
      return     (1);
    }

    for (procnum = procmax = sortnbr = 0, loadglbmax = 0;   /* For all potential contributions */
         procnum < grafptr->s.procglbnbr; procnum ++) {
      if (loadglbtab[procnum] > 0) {
        if (loadglbtab[procnum] > loadglbmax) {   /* Find maximum contribution index as starting point */
          loadglbmax = loadglbtab[procnum];
          procmax    = procnum;
        }
        sorttab[sortnbr].veloval = loadglbtab[procnum];
        sorttab[sortnbr].procnum = (Gnum) procnum;
        sortnbr ++;
      }
    }

    procglbnbr = grafptr->s.procglbnbr;
    for (sortnum = 0; sortnum < sortnbr; sortnum ++) /* Increase priority value from found maximum */
      sorttab[sortnum].prioval = (sorttab[sortnum].procnum + procglbnbr - procmax) % procglbnbr;

    intSort3asc2 (sorttab, sortnbr);              /* Sort contributions array unambiguously */

    proclocnum = grafptr->s.proclocnum;
    for (sortnum = sortnbr - 1; sortnum >= 0; sortnum --) { /* Process contributions by descending load */
      Gnum                compglbloadtmp;

      compglbloadtmp = compglbloaddlt - sorttab[sortnum].veloval;
      if (compglbloadtmp < compglbloaddltmax) {   /* If entire move would cause imbalance */
        Gnum                  fronlocnum;
        BdgraphBipartExMove * movetab;
        Gnum                  movenbr;
        Gnum                  movenum;

        if (sorttab[sortnum].procnum != proclocnum) /* If this is not our job to handle it */
          break;                                  /* Nothing more to do for us             */

        movetab = (BdgraphBipartExMove *) sorttab; /* TRICK: recycle sorttab array as move array */

        for (fronlocnum = frstloctab[gainidx], movenbr = 0; /* For all vertices in swapped gain slot */
             fronlocnum != ~0; fronlocnum = nextloctab[fronlocnum], movenbr ++) {
          Gnum                vertlocnum;

          vertlocnum = fronloctab[fronlocnum];
          movetab[movenbr].veloval = (veloloctax != NULL) ? veloloctax[vertlocnum] : 1;
          movetab[movenbr].vertnum = vertlocnum;
        }

        intSort2asc1 (movetab, movenbr);          /* Sort local moves by ascending load order */

        for (movenum = movenbr - 1; movenum >= 0; movenum --) { /* For all potential moves by descending weight */
          Gnum                compglbloadtmp;

          compglbloadtmp = compglbloaddlt - movetab[movenum].veloval;
          if (compglbloadtmp >= compglbloaddltmax) { /* If move reduces imbalance                                */
            partgsttax[movetab[movenum].vertnum] = (GraphPart) (partval | 2); /* Swap vertex part and flag vertex */
            compglbloaddlt  = compglbloadtmp;
            complocloaddlt += movetab[movenum].veloval; /* We moved that much load locally */
            if (compglbloaddlt <= compglbloaddltmax) /* If nothing more to do, exit loop   */
              break;
          }
        }

        break;                                    /* Nothing more to do */
      }

      compglbloaddlt = compglbloadtmp;            /* Accept move entirely */

      if (sorttab[sortnum].procnum == proclocnum) { /* If we are the process to do it */
        Gnum                fronlocnum;

        complocloaddlt += sorttab[sortnum].veloval; /* We moved all of our local load for this gain */

        for (fronlocnum = frstloctab[gainidx]; fronlocnum != ~0; /* For all vertices in swapped gain slot */
             fronlocnum = nextloctab[fronlocnum])
          partgsttax[fronloctab[fronlocnum]] = (GraphPart) (partval | 2); /* Swap vertex part and flag vertex */

        break;                                    /* We did our job; don't care about the rest */
      }
    }
  }
  grafptr->complocload0 -= (2 * partval - 1) * complocloaddlt; /* Update according to load of moved vertices */

  if (dgraphHaloSync (&grafptr->s, partgsttax + grafptr->s.baseval, GRAPHPART_MPI) != 0) {
    errorPrint ("bdgraphBipartEx: communication error (4)");
    return     (1);
  }

  for (vertgstnum = grafptr->s.vertlocnnd;        /* For all received ghosts */
       vertgstnum < grafptr->s.vertgstnnd; vertgstnum ++) {
    if ((partgsttax[vertgstnum] & 2) != 0) {      /* If ghost vertex changed part */
      partgsttax[vertgstnum] &= 1;                /* Put it back in normal state  */
      flagSet (movegsttab, vertgstnum);           /* While recording state change */
    }
  }

  complocsizedlt = 0;                             /* No difference to number of vertices yet */
  for (fronlocnum = 0, fronlocnbr = grafptr->fronlocnbr;
       fronlocnum < fronlocnbr; fronlocnum ++) {
    Gnum                vertlocnum;

    vertlocnum = fronloctab[fronlocnum];

    flagSet (flagloctab, vertlocnum);             /* Record vertex as already seen   */
    if ((partgsttax[vertlocnum] & 2) != 0) {      /* If frontier vertex changed part */
      partgsttax[vertlocnum] &= 1;                /* Put it back in normal state     */
      flagSet (movegsttab, vertlocnum);           /* While recording state change    */
      complocsizedlt ++;                          /* One more vertex changed of part */
    }
  }
  grafptr->complocsize0 -= (2 * partval - 1) * complocsizedlt; /* Update according to number of moved vertices */

  if (grafptr->s.procsidnbr != 0) {               /* Add potential new frontiers to frontier array */
    Gnum                vertlocnum;
    Gnum                procsidnbr;
    Gnum                procsidnum;
    int                 procsidval;

    const int * restrict const  procsidtab = grafptr->s.procsidtab;

    vertlocnum = grafptr->s.baseval;
    procsidnbr = grafptr->s.procsidnbr;
    procsidnum = 0;
    procsidval = procsidtab[procsidnum ++];  

    while (1) {                                   /* Scan all vertices which have foreign neighbors */
      while (procsidval < 0) {
        vertlocnum -= (Gnum) procsidval;
        procsidval  = procsidtab[procsidnum ++];  
      }

      if (flagVal (flagloctab, vertlocnum) == 0) { /* If vertex not already processed */
        flagSet (flagloctab, vertlocnum);
        fronloctab[fronlocnbr ++] = vertlocnum;   /* Record candidate frontier vertex */
      }

      do {
        if (procsidnum >= procsidnbr)
          goto loop_exit;
      } while ((procsidval = procsidtab[procsidnum ++]) >= 0);
    }
loop_exit : ;
  }

  edlolocval      = 1;                            /* Assume no edge loads */
  commlocgain     =
  commlocgainextn = 0;
  for (fronlocnum = 0, vertlocnnd = grafptr->s.vertlocnnd; /* For all potential frontier vertices */
       fronlocnum < fronlocnbr; ) {
    Gnum                vertlocnum;
    Gnum                edgelocnum;
    Gnum                commcut;
    Gnum                flagval;

    vertlocnum = fronloctab[fronlocnum];
    partval    = partgsttax[vertlocnum];
    flagval    = flagVal (movegsttab, vertlocnum);

    for (edgelocnum = vertloctax[vertlocnum], commcut = 0;
         edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
      Gnum                vertlocend;
      Gnum                flagend;
      Gnum                partend;
      Gnum                partdlt;

      vertlocend = edgegsttax[edgelocnum];
      partend = (Gnum) partgsttax[vertlocend];
      flagend = (Gnum) flagVal (movegsttab, vertlocend);
      if (edloloctax != NULL)
        edlolocval = edloloctax[edgelocnum];

      partdlt  = partval ^ partend;               /* Compute difference between new parts */
      commcut |= partdlt;                         /* Accumulate difference                */

      if ((partdlt != 0) &&                       /* If vertices belong to different parts */
          (vertlocend < vertlocnnd) &&            /* And end vertex is local               */
          (flagVal (flagloctab, vertlocend) == 0)) { /* And end vertex not already queued  */
        fronloctab[fronlocnbr ++] = vertlocend;   /* Add end vertex to frontier queue      */
        flagSet (flagloctab, vertlocend);         /* Flag it to avoid multiple insertions  */
      }

      commlocgain += (((- partdlt) & edlolocval) - /* Compute difference between new and old communication loads */
                      ((- (partdlt ^ flagval ^ flagend)) & edlolocval));
    }

    if (veexloctax != NULL)
      commlocgainextn += (partval - (partval ^ flagval)) * veexloctax[vertlocnum]; /* Compute half of difference in external load */

    if (commcut == 0) {                           /* If vertex no longer belongs to frontier */
      fronloctab[fronlocnum] = fronloctab[-- fronlocnbr]; /* Replace it by another one       */
      continue;                                   /* Process replacement vertex              */
    }

    fronlocnum ++;                                /* Process next vertex */
  }
  grafptr->fronlocnbr = fronlocnbr;               /* Set new number of frontier vertices */

  memFree (gainglbtab);                           /* Free group leader */

  reduloctab[0] = grafptr->fronlocnbr;
  reduloctab[1] = grafptr->complocload0;
  reduloctab[2] = grafptr->complocsize0;
  reduloctab[3] = commlocgain * domndist;         /* Send internal gain */
  reduloctab[4] = commlocgainextn * 2;            /* Send external gain */
  if (MPI_Allreduce (reduloctab, reduglbtab, 5, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartEx: communication error (5)");
    return     (1);
  }
  grafptr->fronglbnbr       = reduglbtab[0];
  grafptr->compglbload0     = reduglbtab[1];
  grafptr->compglbload0dlt  = reduglbtab[1] - grafptr->compglbload0avg;
  grafptr->compglbsize0     = reduglbtab[2];
  grafptr->commglbload     += reduglbtab[3] / 2 + reduglbtab[4]; /* Add modifications, counted twice for internal gain */
  grafptr->commglbgainextn -= reduglbtab[4];      /* Account for modifications in external gain                        */
  grafptr->bbalglbval       = (double) ((grafptr->compglbload0dlt < 0) ? (- grafptr->compglbload0dlt) : grafptr->compglbload0dlt) / (double) grafptr->compglbload0avg;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (grafptr) != 0) {
    errorPrint ("bdgraphBipartEx: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}
