/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : common.c                                **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Pierre RAMET                            **/
/**                Yves SECRETAN (v5.1)                    **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are common routines used    **/
/**                by all modules.                         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 08 may 1998     **/
/**                                 to     14 sep 1998     **/
/**                # Version 2.0  : from : 27 sep 2004     **/
/**                                 to     27 sep 2004     **/
/**                # Version 5.1  : from : 27 jun 2010     **/
/**                                 to     23 nov 2010     **/
/**                # Version 6.0  : from : 21 sep 2013     **/
/**                                 to     21 sep 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"

/*******************/
/*                 */
/* Timing routine. */
/*                 */
/*******************/

double
clockGet (void)
{
#ifdef MPI_INT
  return (MPI_Wtime ());
#else /* MPI_INT */
#if defined COMMON_WINDOWS
  double              res = 0.0;
  LARGE_INTEGER       fq;
  if (QueryPerformanceFrequency (&fq) == 0) {
    FILETIME            ft;
    ULARGE_INTEGER      t;

    GetSystemTimeAsFileTime (&ft);
    t.LowPart  = ft.dwLowDateTime;
    t.HighPart = ft.dwHighDateTime;
    res = (double) t.QuadPart / 10000000.0;
  }
  else {
    LARGE_INTEGER       pc;

    QueryPerformanceCounter (&pc);
    res = (double) pc.QuadPart / (double) fq.QuadPart;
  }
  return (res);
#elif defined COMMON_TIMING_OLD                   /* Old Unix timing routine */
  struct rusage       data;

  getrusage (RUSAGE_SELF, &data);

  return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
          ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) * 1.0e-6L);
#else /* COMMON_TIMING_OLD */
#if defined (_POSIX_TIMERS) && (_POSIX_TIMERS >= 200112L)
  struct timespec     tp;

  clock_gettime (CLOCK_REALTIME, &tp);            /* Elapsed time */

  return ((double) tp.tv_sec + (double) tp.tv_nsec * 1.0e-9L);
#else /* defined (_POSIX_TIMERS) && (_POSIX_TIMERS >= 200112L) */
  struct timeval      tv;

  gettimeofday (&tv, NULL);

 return ((double) tv.tv_sec + (double) tv.tv_usec * 1.0e-6L);
#endif /* defined (_POSIX_TIMERS) && (_POSIX_TIMERS >= 200112L) */
#endif /* COMMON_TIMING_OLD */
#endif /* MPI_INT */
}

/***************************/
/*                         */
/* Usage printing routine. */
/*                         */
/***************************/

void
usagePrint (
FILE * const                stream,
const char ** const         data)
{
  const char **       cptr;

  fprintf (stream, "Usage is:\n");
  for (cptr = data; *cptr != NULL; cptr ++)
    fprintf (stream, "  %s\n", *cptr);
}
