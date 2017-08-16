/* Copyright 2004,2007-2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : common.h                                **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Pierre RAMET                            **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the common data         **/
/**                declarations for all modules.           **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 08 may 1998     **/
/**                                 to   : 08 jan 2001     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to   : 06 jun 2002     **/
/**                # Version 2.0  : from : 13 jun 2005     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 23 nov 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     01 mar 2015     **/
/**                                                        **/
/************************************************************/

#define COMMON_H

/*
** The includes.
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE               600
#endif /* _XOPEN_SOURCE */
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K                             /* For POSIX pthread_barrier_t */
#endif /* __USE_XOPEN2K */

#include            <ctype.h>
#include            <fcntl.h>                     /* Fow Windows _pipe () call */
#include            <math.h>
#include            <memory.h>
#include            <stdio.h>
#include            <stdarg.h>
#include            <stdlib.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include            <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#ifdef HAVE_MALLOC_H
#include            <malloc.h>                    /* Deprecated, but required on some old systems */
#endif /* HAVE_MALLOC_H */
#include            <string.h>
#include            <strings.h>
#include            <time.h>                      /* For the effective calls to clock () */
#include            <limits.h>
#include            <float.h>
#include            <sys/types.h>
#if ((defined COMMON_TIMING_OLD) || (defined HAVE_SYS_TIME_H))
#include            <sys/time.h>
#endif /* ((defined COMMON_TIMING_OLD) || (defined HAVE_SYS_TIME_H)) */
#if ((defined COMMON_TIMING_OLD) || (defined HAVE_SYS_RESOURCE_H))
#include            <sys/resource.h>
#endif /* ((defined COMMON_TIMING_OLD) || (defined HAVE_SYS_RESOURCE_H)) */
#if ((defined COMMON_WINDOWS) || (defined HAVE_WINDOWS_H))
#include            <io.h>                        /* For _pipe () */
#include            <stddef.h>                    /* For intptr_t */
#include            <windows.h>
#endif /* ((defined COMMON_WINDOWS) || (defined HAVE_WINDOWS_H)) */
#if ((! defined COMMON_WINDOWS) && (! defined HAVE_NOT_UNISTD_H))
#include            <unistd.h>
#endif /* ((! defined COMMON_WINDOWS) && (! defined HAVE_NOT_UNISTD_H)) */

#ifdef SCOTCH_PTSCOTCH
#include            <mpi.h>
#endif /* SCOTCH_PTSCOTCH */

#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
#include            <pthread.h>
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */

/*
**  Working definitions.
*/

#ifdef COMMON_MEMORY_TRACE
#define memAlloc(size)              memAllocRecord ((size) | 8)
#define memRealloc(ptr,size)        memReallocRecord ((ptr), ((size) | 8))
#define memFree(ptr)                (memFreeRecord ((void *) (ptr)), 0)
#else /* COMMON_MEMORY_TRACE */
#define memAlloc(size)              malloc ((size) | 8) /* For platforms which return NULL for malloc(0) */
#define memRealloc(ptr,size)        realloc ((ptr),((size) | 8))
#define memFree(ptr)                (free ((char *) (ptr)), 0)
#endif /* COMMON_MEMORY_TRACE */

#define memSet(ptr,val,siz)         memset ((void *) (ptr), (val), (siz))
#define memCpy(dst,src,siz)         memcpy ((void *) (dst), (void *) (src), (siz))
#define memMov(dst,src,siz)         memmove ((void *) (dst), (void *) (src), (siz))

#define MIN(x,y)                    (((x) < (y)) ? (x) : (y))
#define MAX(x,y)                    (((x) < (y)) ? (y) : (x))
#define ABS(x)                      MAX ((x), -(x))
#define SIGN(x)                     (((x) < 0) ? -1 : 1)

/*
**  Handling of generic types.
*/

#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T))
#define UINT32                      uint32_t
#else /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T)) */
#define UINT32                      u_int32_t
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T)) */

#ifndef INT                                       /* If type not externally overriden */
#ifdef INTSIZE32
#define INT                         int32_t
#define UINT                        UINT32
#define COMM_INT                    MPI_INTEGER4
#define INTSTRING                   "%d"
#else /* INTSIZE32 */
#ifdef INTSIZE64
#define INT                         int64_t
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T))
#define UINT                        uint64_t
#else /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T)) */
#define UINT                        u_int64_t
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_UINT_T)) */
#define COMM_INT                    MPI_LONG_LONG
#define INTSTRING                   "%lld"
#else /* INTSIZE64 */
#ifdef LONG                                       /* Better not use it */
#define INT                         long          /* Long integer type */
#define UINT                        unsigned long
#define COMM_INT                    MPI_LONG
#define INTSTRING                   "%ld"
#else /* LONG */
#define INT                         int           /* Default integer type */
#define UINT                        unsigned int
#define COMM_INT                    MPI_INT       /* Generic MPI integer type */
#define INTSTRING                   "%d"
#endif /* LONG      */
#endif /* INTSIZE64 */
#endif /* INTSIZE32 */
#endif /* INT       */

#ifndef IDX                                       /* If type not externally overriden */
#ifdef IDXSIZE32
#define IDX                         int32_t
#else /* IDXSIZE32 */
#ifdef IDXSIZE64
#define IDX                         int64_t
#else /* IDXSIZE64 */
#define IDX                         INT
#endif /* IDXSIZE64 */
#endif /* IDXSIZE32 */
#endif /* IDX       */

#ifndef INTSIZEBITS
#define INTSIZEBITS                 (sizeof (INT) << 3)
#endif /* INTSIZEBITS */

#define INTVALMAX                   ((INT) (((UINT) 1 << (INTSIZEBITS - 1)) - 1))

#define byte unsigned char                        /* Byte type */
#ifndef BYTE
#define BYTE                        byte
#endif /* BYTE */
#ifndef COMM_BYTE
#define COMM_BYTE                   MPI_BYTE
#endif /* COMM_BYTE */
#define COMM_PART                   COMM_BYTE

/*
**  Handling of pseudo-random numbers.
*/

/* The pseudo-random state structure. It is
   based on a Mersenne twister generator, also
   referred to as MT19937. */

typedef struct IntRandState_ {
  UINT32                    randtab[624];         /* State vector */
  int                       randnum;              /* Index value  */
} IntRandState;

/*
**  Handling of flag arrays.
*/

#define flagSize(n)                 (((n) + (sizeof (int) << 3) - 1) / (sizeof (int) << 3))
#define flagVal(a,n)                (((a)[(n) / (sizeof (int) << 3)] >> ((n) & ((sizeof (int) << 3) - 1))) & 1)
#define flagSet(a,n)                (a)[(n) / (sizeof (int) << 3)] |= (1 << ((n) & ((sizeof (int) << 3) - 1)))

/*
**  Handling of timers.
*/

/** The clock type. **/

typedef struct Clock_ {
  double                    time[2];              /*+ The start and accumulated times +*/
} Clock;

/*
**  Handling of Windows constructs.
*/

#if defined COMMON_WINDOWS
#define pipe(fd)                    _pipe (fd, 32768, O_BINARY)
#endif /* COMMON_WINDOWS */

/*
**  Handling of threads.
*/

/** The thread creation flags **/

#define THREADNONE                  0x0000        /* Thread capabilities */

#define THREADHASBARRIER            0x0001

#define THREADCANBARRIER            THREADHASBARRIER
#define THREADCANSCAN               THREADHASBARRIER
#define THREADCANREDUCE             THREADHASBARRIER

/** The thread barrier structure and routines **/

#ifdef COMMON_PTHREAD_BARRIER

#ifndef PTHREAD_BARRIER_SERIAL_THREAD
#define PTHREAD_BARRIER_SERIAL_THREAD -1
#endif /* PTHREAD_BARRIER_SERIAL_THREAD */

typedef struct ThreadBarrier_ {
  int                       thrdnbr;              /*+ Number of threads to wait for       +*/
  volatile int              thrdcur;              /*+ Number of threads currently blocked +*/
  volatile int              instnum;              /*+ Number of barrier instance          +*/
  pthread_mutex_t           mutedat;
  pthread_cond_t            conddat;
} ThreadBarrier;

int                         threadBarrierDestroy (ThreadBarrier *);
int                         threadBarrierInit   (ThreadBarrier *, void *, int); /* Thread attribute not used */
int                         threadBarrierWait   (ThreadBarrier *);

#else /* COMMON_PTHREAD_BARRIER */

#define ThreadBarrier               pthread_barrier_t

#define threadBarrierDestroy        pthread_barrier_destroy
#define threadBarrierInit           pthread_barrier_init
#define threadBarrierWait           pthread_barrier_wait

#endif /* COMMON_PTHREAD_BARRIER */

#define threadBarrier(t)            threadBarrierWait (&(((ThreadGroupHeader *) (((ThreadHeader *) (void *) (t))->grouptr))->barrdat))

/** The thread service routines auxiliary function types **/

typedef int (* ThreadLaunchJoinFunc) (void * const, void * const);
typedef int (* ThreadLaunchStartFunc) (void * const);
typedef void (* ThreadReduceFunc) (void * const, void * const, void * const);
typedef void (* ThreadScanFunc)   (void * const, void * const, void * const, const int);

/** The thread group header block. **/

typedef struct ThreadGroupHeader_ {
#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
  int                       flagval;              /*+ Thread block flags       +*/
  size_t                    datasiz;              /*+ Size of data array cell  +*/
  int                       thrdnbr;              /*+ Number of threads        +*/
  ThreadLaunchStartFunc     stafptr;              /*+ Pointer to start routine +*/
  ThreadLaunchJoinFunc      joifptr;              /*+ Pointer to join routine  +*/
  ThreadBarrier             barrdat;              /*+ Barrier data structure   +*/
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */
} ThreadGroupHeader;

/** The thread header block. **/

typedef struct ThreadHeader_ {
  void *                    grouptr;              /*+ Pointer to thread group +*/
#if ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD))
  pthread_t                 thidval;              /*+ Thread ID               +*/
  int                       thrdnum;              /*+ Thread instance number  +*/
#endif /* ((defined COMMON_PTHREAD) || (defined SCOTCH_PTHREAD)) */
} ThreadHeader;

/** The number of threads **/

#ifdef SCOTCH_PTHREAD

#ifndef SCOTCH_PTHREAD_NUMBER
#define SCOTCH_PTHREAD_NUMBER       1
#endif /* SCOTCH_PTHREAD_NUMBER */

#else /* SCOTCH_PTHREAD */

#ifdef SCOTCH_PTHREAD_NUMBER
#undef SCOTCH_PTHREAD_NUMBER
#endif /* SCOTCH_PTHREAD_NUMBER */
#define SCOTCH_PTHREAD_NUMBER       1

#endif /* SCOTCH_PTHREAD */

/*
**  Handling of files.
*/

/** The file structure. **/

typedef struct File_ {
  char *                    modeptr;              /*+ Opening mode  +*/
  char *                    nameptr;              /*+ File name     +*/
  FILE *                    fileptr;              /*+ File pointer  +*/
  char *                    dataptr;              /*+ Array to free +*/
} File;

/*
**  Function prototypes.
*/

void *                      memAllocGroup       (void **, ...);
void *                      memReallocGroup     (void *, ...);
void *                      memOffset           (void *, ...);
#ifdef COMMON_MEMORY_TRACE
void *                      memAllocRecord      (size_t);
void *                      memReallocRecord    (void * const, size_t);
void                        memFreeRecord       (void * const);
IDX                         memCur              (); /* What is internally an intptr_t has to be turned into an interface type */
IDX                         memMax              ();
#endif /* COMMON_MEMORY_TRACE */

void                        usagePrint          (FILE * const, const char (* []));

void                        fileBlockInit       (File * const, const int);
int                         fileBlockOpen       (File * const, const int);
int                         fileBlockOpenDist   (File * const, const int, const int, const int, const int);
void                        fileBlockClose      (File * const, const int);
FILE *                      fileCompress        (FILE * const, const int);
int                         fileCompressType    (const char * const);
FILE *                      fileUncompress      (FILE * const, const int);
int                         fileUncompressType  (const char * const);
int                         fileNameDistExpand  (char ** const, const int, const int, const int);

void                        errorProg           (const char * const);
void                        errorPrint          (const char * const, ...);
void                        errorPrintW         (const char * const, ...);

int                         intLoad             (FILE * const, INT * const);
int                         intSave             (FILE * const, const INT);
void                        intAscn             (INT * const, const INT, const INT);
void                        intPerm             (INT * const, const INT);
void                        intRandInit         (void);
void                        intRandProc         (int);
void                        intRandReset        (void);
void                        intRandSeed         (INT);
#ifndef COMMON_RANDOM_SYSTEM
INT                         intRandVal          (INT);
#endif /* COMMON_RANDOM_SYSTEM */
void                        intSort1asc1        (void * const, const INT);
void                        intSort2asc1        (void * const, const INT);
void                        intSort2asc2        (void * const, const INT);
void                        intSort3asc1        (void * const, const INT);
void                        intSort3asc2        (void * const, const INT);
INT                         intSearchDicho      (const INT * const, const INT, const INT, const INT);
INT                         intGcd              (INT, INT);

void                        clockInit           (Clock * const);
void                        clockStart          (Clock * const);
void                        clockStop           (Clock * const);
double                      clockVal            (Clock * const);
double                      clockGet            (void);

void                        stringSubst         (char * const, const char * const, const char * const);

#ifdef COMMON_PTHREAD
int                         threadLaunch        (void * const, void * const, const size_t, int (*) (void *), int (*) (void *, void *), const int, const int);
void                        threadReduce        (void * const, void * const, ThreadReduceFunc const, const int);
void                        threadScan          (void * const, void * const, ThreadScanFunc const);
#endif /* COMMON_PTHREAD */

/*
**  Macro definitions.
*/

#define clockInit(clk)              ((clk)->time[0]  = (clk)->time[1] = 0)
#define clockStart(clk)             ((clk)->time[0]  = clockGet ())
#define clockStop(clk)              ((clk)->time[1] += (clockGet () - (clk)->time[0]))
#define clockVal(clk)               ((clk)->time[1])

#define fileBlockFile(b,i)          ((b)[i].fileptr)
#define fileBlockMode(b,i)          ((b)[i].modeptr)
#define fileBlockName(b,i)          ((b)[i].nameptr)

#ifdef COMMON_RANDOM_SYSTEM
#ifdef COMMON_RANDOM_RAND
#define intRandVal(ival)            ((INT) (((UINT) rand ()) % ((UINT) (ival))))
#else /* COMMON_RANDOM_RAND */
#define intRandVal(ival)            ((INT) (((UINT) random ()) % ((UINT) (ival))))
#endif /* COMMON_RANDOM_RAND */
#endif /* COMMON_RANDOM_SYSTEM */

#define DATASIZE(n,p,i)             ((INT) (((n) + ((p) - 1 - (i))) / (p)))
#define DATASCAN(n,p,i)             ((i) * ((INT) (n) / (INT) (p)) + (((i) > ((n) % (p))) ? ((n) % (p)) : (i)))

#define FORTRAN(nu,nl,pl,pc)        FORTRAN2(REPLACE(nu),REPLACE(nl),pl,pc)
#define FORTRAN2(nu,nl,pl,pc)                    \
void nu pl;                                      \
void nl pl                                       \
{ nu pc; }                                       \
void GLUE(nl,_) pl                               \
{ nu pc; }                                       \
void GLUE(nl,__) pl                              \
{ nu pc; }                                       \
void nu pl

#define REPLACE(s)                  s
#define GLUE(p,s)                   p##s

#define STRINGIFY2(n)               #n
#define STRINGIFY(n)                STRINGIFY2(n)
