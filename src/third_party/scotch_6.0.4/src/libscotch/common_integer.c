/* Copyright 2004,2007-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : common_integer.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the generic integer **/
/**                type.                                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 sep 1998     **/
/**                                 to     22 sep 1998     **/
/**                # Version 0.1  : from : 07 jan 2002     **/
/**                                 to     17 jan 2003     **/
/**                # Version 1.0  : from : 23 aug 2005     **/
/**                                 to   : 19 dec 2006     **/
/**                # Version 2.0  : from : 26 feb 2008     **/
/**                                 to   : 26 feb 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 16 jul 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     13 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON_INTEGER

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"

/********************************/
/*                              */
/* Basic routines for fast I/O. */
/*                              */
/********************************/

/* Fast read for INT values.
** It returns:
** - 1  : on success.
** - 0  : on error.
*/

int
intLoad (
FILE * const                stream,               /*+ Stream to read from     +*/
INT * const                 valptr)               /*+ Area where to put value +*/
{
  int                 sign;                       /* Sign flag      */
  int                 car;                        /* Character read */
  INT                 val;                        /* Value          */

  sign = 0;                                       /* Assume positive constant     */
  for ( ; ; ) {                                   /* Consume whitespaces and sign */
    car = getc (stream);
    if (isspace (car))
      continue;
    if ((car >= '0') && (car <= '9'))
      break;
    if (car == '-') {
      sign = 1;
      car  = getc (stream);
      break;
    }
    if (car == '+') {
      car = getc (stream);
      break;
    }
    return (0);
  }
  if ((car < '0') || (car > '9'))                 /* If first char is non numeric */
    return (0);                                   /* Then it is an error          */
  val = car - '0';                                /* Get first digit              */
  for ( ; ; ) {
    car = getc (stream);
    if ((car < '0') || (car > '9')) {
      ungetc (car, stream);
      break;
    }
    val = val * 10 + (car - '0');                 /* Accumulate digits */
  }
  *valptr = (sign != 0) ? (- val) : val;          /* Set result */

  return (1);
}

/* Write routine for INT values.
** It returns:
** - 1  : on success.
** - 0  : on error.
*/

int
intSave (
FILE * const                stream,               /*+ Stream to write to +*/
const INT                   val)                  /*+ Value to write     +*/
{
  return ((fprintf (stream, INTSTRING, (INT) val) == EOF) ? 0 : 1);
}

/**********************************/
/*                                */
/* Permutation building routines. */
/*                                */
/**********************************/

/* This routine fills an array with
** consecutive INT values, in
** ascending order.
** It returns:
** - VOID  : in all cases.
*/

void
intAscn (
INT * const                 permtab,              /*+ Permutation array to build +*/
const INT                   permnbr,              /*+ Number of entries in array +*/
const INT                   baseval)              /*+ Base value                 +*/
{
  INT *               permtax;
  INT                 permnum;
  INT                 permnnd;

  for (permnum = baseval, permnnd = baseval + permnbr, permtax = permtab - baseval;
       permnum < permnnd; permnum ++)
    permtax[permnum] = permnum;
}

/* This routine computes a random permutation
** of an array of INT values.
** It returns:
** - VOID  : in all cases.
*/

void
intPerm (
INT * const                 permtab,              /*+ Permutation array to build +*/
const INT                   permnbr)              /*+ Number of entries in array +*/
{
  INT *               permptr;
  INT                 permrmn;

  for (permptr = permtab, permrmn = permnbr;      /* Perform random permutation */
       permrmn > 0; permptr ++, permrmn --) {
    INT                 permnum;
    INT                 permtmp;

    permnum          = intRandVal (permrmn);      /* Select index to swap       */
    permtmp          = permptr[0];                /* Swap it with current index */
    permptr[0]       = permptr[permnum];
    permptr[permnum] = permtmp;
  }
}

/*************************************/
/*                                   */
/* Pseudo-random generator routines. */
/*                                   */
/*************************************/

static volatile int         intrandflag = 0;      /*+ Flag set if generator already initialized +*/
static UINT32               intrandproc = 0;      /*+ Process number                            +*/
static UINT32               intrandseed = 1;      /*+ Pseudo-random seed                        +*/

/* This routine sets the process number that is
** used to generate a different seed across all
** processes. In order for this number to be
** taken into account, it must be followed by
** a subsequent call to intRandInit(),
** intRandReset() or intRandSeed().
** It returns:
** - VOID  : in all cases.
*/

void
intRandProc (
int                         procnum)
{
  intrandproc = (UINT32) procnum;                 /* Set process number */
}

/* This routine initializes the seed used by Scotch
** with the provided value. Hence, all subsequent
** calls to intRandInit() will start from this seed.
** It returns:
** - VOID  : in all cases.
*/

#ifndef COMMON_RANDOM_SYSTEM
static IntRandState         intrandstat;          /*+ Pseudo-random state value +*/

static
void
intRandSeed3 (
IntRandState * restrict     randptr,
UINT32                      randval)
{
  UINT32              randtmp;
  UINT32              i;

  UINT32 * restrict const randtab = randptr->randtab; /* Fast access */

  randtmp    = (UINT32) randval;
  randtab[0] = randtmp;                           /* Reset array contents */
  for (i = 1; i < 623; i ++) {
    randtmp = 0x6c078965 * randtmp ^ (randtmp >> 30) + i;
    randtab[i] = randtmp;
  }
  randptr->randnum = 0;                           /* Reset array index */
}
#endif /* COMMON_RANDOM_SYSTEM */

static
void
intRandSeed2 (
UINT32                      seedval)
{
  UINT32              randtmp;

  randtmp = seedval * (intrandproc + 1);          /* Account for process index */

#ifdef COMMON_RANDOM_SYSTEM
#ifdef COMMON_RANDOM_RAND
  srand ((unsigned int) randtmp);
#else /* COMMON_RANDOM_RAND */
  srandom ((unsigned int) randtmp);
#endif /* COMMON_RANDOM_RAND */
#else /* COMMON_RANDOM_SYSTEM */
  intRandSeed3 (&intrandstat, randtmp);           /* Initialize state vector from random seed */
#endif /* COMMON_RANDOM_SYSTEM */
}

void
intRandSeed (
INT                         seedval)
{
  intrandflag = 1;                                /* Generator has been initialized */
  intrandseed = (UINT32) seedval;                 /* Save new seed                  */

  intRandSeed2 (intrandseed);                     /* Initialize pseudo-random seed */
}

/* This routine initializes the pseudo-random
** generator if necessary. In order for multi-sequential
** programs to have exactly the same behavior on any
** process, the random seed does not depend on process
** rank. This routine is not really thread-safe, so it
** should not be called concurrently when it has never
** been initialized before.
** It returns:
** - VOID  : in all cases.
*/

void
intRandInit (void)
{
  if (intrandflag == 0) {                         /* Non thread-safe check          */
    intrandflag = 1;                              /* Generator has been initialized */

#if ! ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED) || (defined SCOTCH_DETERMINISTIC))
    intrandseed = (UINT32) time (NULL);           /* Set random seed if needed */
#endif /* ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED) || (defined SCOTCH_DETERMINISTIC)) */
    intRandSeed2 (intrandseed);                   /* Initialize state vector from seed */
  }
}

/* This routine reinitializes the pseudo-random
** generator to its initial value. This routine
** is not thread-safe.
** It returns:
** - VOID  : in all cases.
*/

void
intRandReset (void)
{
  if (intrandflag == 0)                           /* Keep seed computed during first initialization */
    intRandInit ();

  intRandSeed2 (intrandseed);
}

/* This routine computes a new pseudo-random
** 32bit value from the state that is passed
** to it.
** For speed and reproducibility reasons,
** this routine is not thread-safe. Providing
** a thread-safe routine would mean determinism
** could not be achieved in caller routines.
** It is the responsibility of application
** routines to call intRandVal() in a way that
** avoids concurrent execution and potentially
** enforces reproducibility.
** It returns:
** - x  : pseudo-random value.
*/

#ifndef COMMON_RANDOM_SYSTEM
static
UINT32
intRandVal2 (
IntRandState * restrict     randptr)
{
  int                 randnum;
  UINT32              randval;

  UINT32 * restrict const randtab = randptr->randtab; /* Fast access */

#ifdef COMMON_DEBUG
  if (intrandflag == 0) {
    errorPrint ("intRandVal2: random generator not initialized");
    return     (~0);
  }
#endif /* COMMON_DEBUG */

  randnum = randptr->randnum;
  if (randnum == 0) {
    int                 i;

    for (i = 0; i < 624; i ++) {
      UINT32              randtmp;

      randtmp = (randtab[i] & 0x80000000) + (randtab[(i + 1) % 624] & 0x7FFFFFFF);
      randtmp = randtab[(i + 397) % 624] ^ (randtmp >> 1);
      if ((randtmp & 1) != 0)
        randtmp ^= 0x9908B0DF;

      randtab[i] = randtmp;
    }
  }

  randval  = randtab[randnum];
  randval ^= (randval >> 11);
  randval ^= (randval >> 7) & 0x9D2C5680;
  randval ^= (randval >> 15) & 0xEFC60000;
  randval ^= (randval >> 18);
  randptr->randnum = (randnum + 1) % 624;

  return (randval);
}
#endif /* COMMON_RANDOM_SYSTEM */

/* This routine returns a pseudo-random integer
** value in the range [0..randmax[. This routine
** is not thread-safe as it uses a global state
** variable.
** It returns:
** - x  : pseudo-random value.
*/

#ifndef COMMON_RANDOM_SYSTEM
INT
intRandVal (
INT                         randmax)
{
  return (((UINT) intRandVal2 (&intrandstat)) % randmax);
}
#endif /* COMMON_RANDOM_SYSTEM */

/*********************/
/*                   */
/* Sorting routines. */
/*                   */
/*********************/

/* This routine sorts an array of
** INT values in ascending order
** by their first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort1asc1
#define INTSORTSIZE                 (sizeof (INT))
#define INTSORTSWAP(p,q)            do { INT t; t = *((INT *) (p)); *((INT *) (p)) = *((INT *) (q)); *((INT *) (q)) = t; } while (0)
#define INTSORTCMP(p,q)             (*((INT *) (p)) < *((INT *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of pairs of
** INT values in ascending order by their
** first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort2asc1
#define INTSORTSIZE                 (2 * sizeof (INT))
#define INTSORTSWAP(p,q)            do { INT t, u; t = *((INT *) (p)); u = *((INT *) (p) + 1); *((INT *) (p)) = *((INT *) (q)); *((INT *) (p) + 1) = *((INT *) (q) + 1); *((INT *) (q)) = t; *((INT *) (q) + 1) = u; } while (0)
#define INTSORTCMP(p,q)             (*((INT *) (p)) < *((INT *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of pairs of
** INT values in ascending order by both
** of their values, used as primary and
** secondary keys.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort2asc2
#define INTSORTSIZE                 (2 * sizeof (INT))
#define INTSORTSWAP(p,q)            do { INT t, u; t = *((INT *) (p)); u = *((INT *) (p) + 1); *((INT *) (p)) = *((INT *) (q)); *((INT *) (p) + 1) = *((INT *) (q) + 1); *((INT *) (q)) = t; *((INT *) (q) + 1) = u; } while (0)
#define INTSORTCMP(p,q)             ((*((INT *) (p)) < *((INT *) (q))) || ((*((INT *) (p)) == *((INT *) (q))) && (*((INT *) (p) + 1) < *((INT *) (q) + 1))))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of 3-uples of
** INT values in ascending order by their
** first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort3asc1
#define INTSORTSIZE                 (3 * sizeof (INT))
#define INTSORTSWAP(p,q)            do { INT t, u, v; t = *((INT *) (p)); u = *((INT *) (p) + 1); v = *((INT *) (p) + 2); *((INT *) (p)) = *((INT *) (q)); *((INT *) (p) + 1) = *((INT *) (q) + 1); *((INT *) (p) + 2) = *((INT *) (q) + 2); *((INT *) (q)) = t; *((INT *) (q) + 1) = u; *((INT *) (q) + 2) = v; } while (0)
#define INTSORTCMP(p,q)             (*((INT *) (p)) < *((INT *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of 3-uples of
** INT values in ascending order by their
** first and second values, used as primary
** and secondary keys.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort3asc2
#define INTSORTSIZE                 (3 * sizeof (INT))
#define INTSORTSWAP(p,q)            do { INT t, u, v; t = *((INT *) (p)); u = *((INT *) (p) + 1); v = *((INT *) (p) + 2); *((INT *) (p)) = *((INT *) (q)); *((INT *) (p) + 1) = *((INT *) (q) + 1); *((INT *) (p) + 2) = *((INT *) (q) + 2); *((INT *) (q)) = t; *((INT *) (q) + 1) = u; *((INT *) (q) + 2) = v; } while (0)
#define INTSORTCMP(p,q)             ((*((INT *) (p)) < *((INT *) (q))) || ((*((INT *) (p)) == *((INT *) (q))) && (*((INT *) (p) + 1) < *((INT *) (q) + 1))))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine computes the greatest common
** divisor of two non-negative integers u and v.
** It returns:
** - x  : the GCD of u and v.
*/

INT
intGcd (
INT                         u,
INT                         v)
{
  INT                 t;

  if (v < u) {                                    /* u should always be the biggest */
    t = u;
    u = v;
    v = t;
  }

  while (v != 0) {
    t = v;
    v = u % v;
    u = t;
  }

  return (u);
}
