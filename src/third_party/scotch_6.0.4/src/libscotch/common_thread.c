/* Copyright 2012-2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : common_thread.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module provides routines to ease   **/
/**                the use of Posix threads.               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 04 jul 2012     **/
/**                                 to     02 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifdef SCOTCH_PTHREAD_AFFINITY_LINUX
#define _GNU_SOURCE
#include <sched.h>
#endif /* SCOTCH_PTHREAD_AFFINITY_LINUX */

#define COMMON_THREAD

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"

/*****************************/
/*                           */
/* Thread handling routines. */
/*                           */
/*****************************/

#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))

/* These routines implement the classical barrier
** operations for systems that do not provide them.
*/

#ifdef COMMON_PTHREAD_BARRIER

int
threadBarrierInit (
ThreadBarrier *             barrptr,
void *                      attrptr,              /* Not used */
int                         thrdnbr)
{
  barrptr->thrdnbr = thrdnbr;
  barrptr->thrdcur = 0;
  barrptr->instnum = 0;
  pthread_mutex_init (&barrptr->mutedat, NULL);
  pthread_cond_init  (&barrptr->conddat, NULL);

  return (0);
}

/*
**
*/

int
threadBarrierDestroy (
ThreadBarrier *             barrptr)
{
  pthread_cond_destroy  (&barrptr->conddat);
  pthread_mutex_destroy (&barrptr->mutedat);

  return (0);
}

/*
**
*/

int
threadBarrierWait (
ThreadBarrier *             barrptr)
{
  int                 instnum;
  int                 thrdcur;
  int                 o;

  pthread_mutex_lock (&barrptr->mutedat);

  thrdcur = barrptr->thrdcur + 1;
  instnum = barrptr->instnum;

  o = 0;                                          /* Assume thread will not be the last one */

  if (thrdcur == barrptr->thrdnbr) {
    barrptr->thrdcur = 0;
    barrptr->instnum = instnum + 1;
    pthread_cond_broadcast (&barrptr->conddat);
    o = PTHREAD_BARRIER_SERIAL_THREAD;            /* Last thread returns special value */
  }
  else {
    barrptr->thrdcur = thrdcur;
    do
      pthread_cond_wait (&barrptr->conddat, &barrptr->mutedat);
    while (barrptr->instnum == instnum);
  }

  pthread_mutex_unlock (&barrptr->mutedat);

  return (o);
}

#endif /* COMMON_PTHREAD_BARRIER */

/* This routine performs a synchronous
** reduction operation on the given block
**  of threads. The routine is called only
** if the two threads in the reduction binary
** tree exist.
** A final, global barrier may be necessary for
** all threads to benefit from the result of the
** reduction operation.
** It returns:
** - void  : in all cases.
*/

void
threadReduce (
void * const                dataptr,              /* Per-thread data block        */
void * const                contptr,              /* Pointer to thread contents   */
ThreadReduceFunc const      redfptr,              /* Pointer to reduction routine */
int                         rootnum)              /* Root of reduction            */
{
  ThreadHeader * restrict const       thrdptr = (ThreadHeader *) dataptr;
  ThreadGroupHeader * restrict const  grouptr = thrdptr->grouptr;
  const size_t                        datasiz = grouptr->datasiz;
  const int                           thrdnbr = grouptr->thrdnbr;
  const int                           thrdnum = thrdptr->thrdnum;
  int                                 thrdnsk;    /* Rank of thread in skewed reduction tree */
  int                                 thrdmsk;

  thrdnsk = (thrdnum + thrdnbr - rootnum) % thrdnbr;
  for (thrdmsk = 1; thrdmsk < thrdnbr; thrdmsk <<= 1) {
    int                 thrdesk;                  /* Skewed rank of end thread */

    threadBarrierWait (&grouptr->barrdat);

    thrdesk = thrdnsk ^ thrdmsk;                  /* Get skewed rank of end thread */

    if (thrdesk < thrdnbr) {                      /* If end thread exists            */
      if (thrdesk > thrdnsk) {                    /* If we are on the receiving side */
        int                 thrdend;
        int                 thrddlt;

        thrdend = (thrdesk + rootnum) % thrdnbr;
        thrddlt = thrdend - thrdnum;
        redfptr (dataptr, contptr, (void *) ((byte *) contptr + thrddlt * datasiz)); /* Call reduction routine */
      }
      else                                        /* We are on the sending side       */
        thrdnsk += thrdnbr;                       /* Make sure we will no longer work */
    }
  }
}

/* This routine performs a synchronous
** scan operation on the given block of
** threads. It requires a dummy area for
** storing every other result, hence the
** phase number that is passed to the
** auxiliary routine.
** It returns:
** - void  : in all cases.
*/

void
threadScan (
void * const                dataptr,              /* Per-thread data block      */
void * const                contptr,              /* Pointer to thread contents */
ThreadScanFunc const        scafptr)              /* Scan function              */
{
  ThreadHeader * restrict const       thrdptr = (ThreadHeader *) dataptr;
  ThreadGroupHeader * restrict const  grouptr = thrdptr->grouptr;
  const size_t                        datasiz = grouptr->datasiz;
  const int                           thrdnbr = grouptr->thrdnbr;
  const int                           thrdnum = thrdptr->thrdnum;
  int                                 thrdmsk;
  int                                 i;

  for (thrdmsk = 1, i = 0; thrdmsk < thrdnbr; thrdmsk <<= 1, i ^= 1) ; /* Determine number of steps to go         */
  if (i != 0)                                     /* If number of steps is odd                                    */
    scafptr (dataptr, contptr, NULL, 0);          /* Pre-copy to swap area so that it will end at the right place */

  for (thrdmsk = 1; thrdmsk < thrdnbr; thrdmsk <<= 1, i ^= 1) {
    int                 thrdend;

    threadBarrierWait (&grouptr->barrdat);        /* Barrier on all threads, even those which do not participate */

    thrdend = thrdnum - thrdmsk;                  /* Get rank of end thread */
    scafptr (dataptr, contptr, (thrdend >= 0) ? (void *) ((byte *) contptr - thrdmsk * datasiz) : NULL, i); /* If end thread exists, perform scan, else just copy */
  }
}

/* This routine actually launches the
** initial thread routine, and performs
** the necessary joins. It is its task
** to set its affinity mask, as in some
** thread implementations threads can
** only set affinity for themselves.
*/

static
void *
threadLaunch2 (
void *                      dataptr)              /* Per-thread data block */
{
  ThreadHeader * restrict const       thrdptr = (ThreadHeader *) dataptr;
  ThreadGroupHeader * restrict const  grouptr = thrdptr->grouptr;
  const size_t                        datasiz = grouptr->datasiz;
  const int                           thrdnbr = grouptr->thrdnbr;
  const int                           thrdnum = thrdptr->thrdnum;
  int                                 thrdmsk;
  int                                 o;
#ifdef SCOTCH_PTHREAD_AFFINITY_LINUX
  cpu_set_t                           cpuset;
#endif /* SCOTCH_PTHREAD_AFFINITY_LINUX */

#ifdef SCOTCH_PTHREAD_AFFINITY_LINUX
  CPU_ZERO (&cpuset);
  CPU_SET  (thrdnum, &cpuset);                    /* Thread sets its own affinity */
  pthread_setaffinity_np (thrdptr->thidval, sizeof (cpu_set_t), &cpuset);
#endif /* SCOTCH_PTHREAD_AFFINITY_LINUX */

  o = grouptr->stafptr (dataptr);                 /* Call start routine */

  for (thrdmsk = 1; thrdmsk < thrdnbr; thrdmsk <<= 1) {
    volatile ThreadHeader * restrict  thrdtmp;    /* Pointer to thread header of other thread */
    int                               thrdend;

    thrdend = thrdnum ^ thrdmsk;                  /* Get rank of end thread       */
    if (thrdend >= thrdnbr)                       /* If end thread does not exist */
      continue;

    thrdtmp = (ThreadHeader *) ((byte *) dataptr + grouptr->datasiz * (thrdend - thrdnum));
    while (thrdtmp->thrdnum == -1) ;              /* Spin-lock until end thread created */

    if (thrdnum > thrdend) {                      /* If we are on the sending side        */
      if (thrdtmp->thrdnum < 0) {                 /* If end thread could not be created   */
        pthread_detach (thrdptr->thidval);        /* Detach since nobody will join for us */
        o = 1;                                    /* Set (useless) error return value     */
      }

      pthread_exit ((void *) (intptr_t) o);       /* Exit anyway */
    }
    else {
      if (thrdtmp->thrdnum < 0)                   /* If end thread could not be created */
        o = 1;                                    /* Set error return value             */
      else {
        void *              o2;

        pthread_join (thrdtmp->thidval, &o2);     /* Get return value from end thread */
        o |= (int) (intptr_t) o2;                 /* Amalgamate return status         */

        if ((grouptr->joifptr != NULL) &&         /* If we have something to do      */
            (o == 0))                             /* And if no error in both threads */
          o |= grouptr->joifptr (dataptr, (void *) ((byte *) dataptr + thrdmsk * datasiz)); /* Call join routine */
      }
    }
  }

  return ((void *) (intptr_t) o);                 /* Thread of rank 0 returns global status */
}

/* This routine launches the given
** number of threads that will run the
** given start routine, and will end up
** calling the given join routine in a
** binary tree fashion.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
threadLaunch (
void * const                gdatptr,              /* Pointer to thread group data block           */
void * const                tdattab,              /* Array of thread data                         */
const size_t                datasiz,              /* Size of individual data array cell           */
ThreadLaunchStartFunc       stafptr,              /* Pointer to start routine                     */
ThreadLaunchJoinFunc        joifptr,              /* Pointer to join routine                      */
const int                   thrdnbr,              /* Number of threads to run (including current) */
const int                   flagval)              /* Flag for thread operation data structures    */
{
  ThreadGroupHeader * const grouptr = (ThreadGroupHeader *) gdatptr;
  ThreadHeader *            thrdptr;
  int                       thrdnum;
  byte *                    dataptr;
  void *                    o;

  grouptr->flagval = flagval;
  grouptr->datasiz = datasiz;
  grouptr->thrdnbr = thrdnbr;
  grouptr->stafptr = stafptr;
  grouptr->joifptr = joifptr;

  if ((flagval & THREADHASBARRIER) != 0) {
    if (threadBarrierInit (&grouptr->barrdat, NULL, thrdnbr) != 0) {
      errorPrint ("threadLaunch: cannot initialize barrier (1)");
      return     (1);
    }
  }

  for (thrdnum = 0, dataptr = (byte *) tdattab;   /* Prepare threads for launching */
       thrdnum < thrdnbr; thrdnum ++, dataptr += datasiz) {
    ThreadHeader *      thrdptr;

    thrdptr = (ThreadHeader *) dataptr;
    thrdptr->thrdnum = -1;                        /* Set threads as not yet launched */
  }

  __sync_synchronize ();                          /* Full memory barrier */

  for (thrdnum = 1, dataptr = (byte *) tdattab + datasiz; /* Launch threads from 1 to (thrdnbr - 1) */
       thrdnum < thrdnbr; thrdnum ++, dataptr += datasiz) {
    ThreadHeader *      thrdptr;

    thrdptr = (ThreadHeader *) dataptr;
    thrdptr->grouptr = gdatptr;
    thrdptr->thrdnum = thrdnum;

    if (pthread_create (&thrdptr->thidval, NULL, threadLaunch2, (void *) dataptr) != 0) {
      errorPrint ("threadLaunch: cannot launch thread (%d)", thrdnum);
      return     (1);
    }
  }

  thrdptr = (ThreadHeader *) tdattab;             /* Run thread 0 */
  thrdptr->grouptr = gdatptr;
  thrdptr->thidval = pthread_self ();
  thrdptr->thrdnum = 0;

  o = threadLaunch2 (tdattab);

  if ((flagval & THREADHASBARRIER) != 0)          /* Free allocated resources */
    threadBarrierDestroy (&grouptr->barrdat);

  return ((int) (intptr_t) o);
}

#else /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

/**********************************/
/*                                */
/* Thread handling routine stubs. */
/*                                */
/**********************************/

void
threadReduce (
void * const                dataptr,              /* Per-thread data block        */
void * const                contptr,              /* Pointer to thread contents   */
ThreadReduceFunc const      redfptr,              /* Pointer to reduction routine */
int                         rootnum)              /* Root of reduction            */
{
  errorPrint ("threadReduce: Scotch not compiled with either COMMON_PTHREAD or SCOTCH_PTHREAD");
}

void
threadScan (
void * const                dataptr,              /* Per-thread data block      */
void * const                contptr,              /* Pointer to thread contents */
ThreadScanFunc const        scafptr)              /* Scan function              */
{
  errorPrint ("threadScan: Scotch not compiled with either COMMON_PTHREAD or SCOTCH_PTHREAD");
}

int
threadLaunch (
void * const                gdatptr,              /* Pointer to thread group data block           */
void * const                tdattab,              /* Array of thread data                         */
const size_t                datasiz,              /* Size of individual data array cell           */
ThreadLaunchStartFunc       stafptr,              /* Pointer to start routine                     */
ThreadLaunchJoinFunc        joifptr,              /* Pointer to join routine                      */
const int                   thrdnbr,              /* Number of threads to run (including current) */
const int                   flagval)              /* Flag for thread operation data structures    */
{
  errorPrint ("threadLaunch: Scotch not compiled with either COMMON_PTHREAD or SCOTCH_PTHREAD");
  return     (1);
}

#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */
