/* Copyright 2012,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_common_thread.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the sequential        **/
/**                strategy building routines.             **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 04 nov 2012     **/
/**                                 to     01 mar 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE               600
#endif /* _XOPEN_SOURCE */
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K                             /* For POSIX pthread_barrier_t */
#endif /* __USE_XOPEN2K */

#include <stdio.h>
#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
#include <pthread.h>
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

#include "../libscotch/module.h"
#include "../libscotch/common.h"

#define COMPVAL(n)                  (((n) * ((n) + 1)) / 2)

/*
**  The type and structure definitions.
*/

/*+ The block data structure +*/

typedef struct TestThreadGroup_ {
  ThreadGroupHeader         thrddat;              /*+ Thread handling data          +*/
  int                       redusum;              /*+ Value to compare reduction to +*/
} TestThreadGroup;

/*+ The thread-specific data block. +*/

typedef struct TestThread_ {
  ThreadHeader              thrddat;              /*+ Thread management data +*/
  int                       reduval;              /*+ Value to reduce        +*/
  int                       scanval;              /*+ Value for scan         +*/
  int                       dummval;              /*+ Dummy value for scan   +*/
} TestThread;

/*************************/
/*                       */
/* The threaded routine. */
/*                       */
/*************************/

#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))

static
void
testReduce (
TestThread * restrict const tlocptr,              /* Pointer to local thread */
void * restrict const       vlocptr,              /* Pointer to local value  */
void * restrict const       vremptr)              /* Pointer to remote value */
{
  TestThread * restrict const tremptr = (TestThread *) vremptr;

  tlocptr->reduval += tremptr->reduval;
}

static
void
testScan (
TestThread * restrict const tlocptr,              /* Pointer to local thread */
int * restrict const        vlocptr,              /* Pointer to local value  */
int * restrict const        vremptr,              /* Pointer to remote value */
const int                   phasval)              /* Phase index             */
{
  vlocptr[1 - phasval] = vlocptr[phasval] + ((vremptr == NULL) ? 0 : vremptr[phasval]);
}

static
int
testThreads (
TestThread * restrict   thrdptr)
{
  TestThreadGroup * restrict const  grouptr = (TestThreadGroup *) (thrdptr->thrddat.grouptr);
  const int                         thrdnbr = grouptr->thrddat.thrdnbr;
  const int                         thrdnum = thrdptr->thrddat.thrdnum;
  int                               o;

  printf ("%d: running\n", thrdnum);

  o = 0;

  if (thrdnum == 0)
    printf ("Performing reduction\n");

  threadBarrier (thrdptr);

  thrdptr->reduval = 1 + thrdptr->thrddat.thrdnum;
  threadReduce (thrdptr, thrdptr, (ThreadReduceFunc) testReduce, 0);

  if ((thrdnum == 0) &&                           /* Test reduction result on thread 0 */
      (thrdptr->reduval != grouptr->redusum)) {
    printf ("0: invalid reduction operator\n");
    o = 1;
  }

  if (thrdnum == 0)
    printf ("Performing scan\n");

  threadBarrier (thrdptr);

  thrdptr->scanval = 1 + thrdptr->thrddat.thrdnum;
  threadScan (thrdptr, &thrdptr->scanval, (ThreadScanFunc) testScan);

  if (thrdptr->scanval != COMPVAL (thrdnum + 1)) {
    printf ("%d: invalid scan operator\n", thrdnum);
    o = 1;
  }

  threadBarrier (thrdptr);

  return (o);
}

#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                 argc,
char *              argv[])
{
  TestThreadGroup       groudat;
#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
  TestThread * restrict thrdtab;
  int                   thrdnbr;
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

  SCOTCH_errorProg (argv[0]);

#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
  thrdnbr = SCOTCH_PTHREAD_NUMBER;

  groudat.redusum = COMPVAL (thrdnbr);

  if ((thrdtab = malloc (thrdnbr * sizeof (TestThread))) == NULL) {
    errorPrint ("main: out of memory");
    return     (1);
  }

  if (threadLaunch (&groudat, thrdtab, sizeof (TestThread), (ThreadLaunchStartFunc) testThreads, (ThreadLaunchJoinFunc) NULL,
                    thrdnbr, THREADCANBARRIER | THREADCANREDUCE | THREADCANSCAN) != 0) {
    errorPrint ("main: cannot launch or run threads");
    return     (1);
  }

  free (thrdtab);
#else /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */
  printf ("Scotch not compiled with either COMMON_PTHREAD or SCOTCH_PTHREAD\n");
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

  return (0);
}
