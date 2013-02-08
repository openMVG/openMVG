/* ========================================================================== */
/* === SuiteSparse_config =================================================== */
/* ========================================================================== */

/* Configuration file for SuiteSparse: a Suite of Sparse matrix packages
 * (AMD, COLAMD, CCOLAMD, CAMD, CHOLMOD, UMFPACK, CXSparse, and others).
 *
 * SuiteSparse_config.h provides the definition of the long integer.  On most
 * systems, a C program can be compiled in LP64 mode, in which long's and
 * pointers are both 64-bits, and int's are 32-bits.  Windows 64, however, uses
 * the LLP64 model, in which int's and long's are 32-bits, and long long's and
 * pointers are 64-bits.
 *
 * SuiteSparse packages that include long integer versions are
 * intended for the LP64 mode.  However, as a workaround for Windows 64
 * (and perhaps other systems), the long integer can be redefined.
 *
 * If _WIN64 is defined, then the __int64 type is used instead of long.
 *
 * The long integer can also be defined at compile time.  For example, this
 * could be added to SuiteSparse_config.mk:
 *
 * CFLAGS = -O -D'SuiteSparse_long=long long' \
 *  -D'SuiteSparse_long_max=9223372036854775801' -D'SuiteSparse_long_idd="lld"'
 *
 * This file defines SuiteSparse_long as either long (on all but _WIN64) or
 * __int64 on Windows 64.  The intent is that a SuiteSparse_long is always a
 * 64-bit integer in a 64-bit code.  ptrdiff_t might be a better choice than
 * long; it is always the same size as a pointer.
 *
 * This file also defines the SUITESPARSE_VERSION and related definitions.
 *
 * Copyright (c) 2012, Timothy A. Davis.  No licensing restrictions apply
 * to this file or to the SuiteSparse_config directory.
 * Author: Timothy A. Davis.
 */

#ifndef _SUITESPARSECONFIG_H
#define _SUITESPARSECONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <limits.h>
#include <stdlib.h>

/* ========================================================================== */
/* === SuiteSparse_long ===================================================== */
/* ========================================================================== */

#ifndef SuiteSparse_long

#ifdef _WIN64

#define SuiteSparse_long __int64
#define SuiteSparse_long_max _I64_MAX
#define SuiteSparse_long_idd "I64d"

#else

#define SuiteSparse_long long
#define SuiteSparse_long_max LONG_MAX
#define SuiteSparse_long_idd "ld"

#endif
#define SuiteSparse_long_id "%" SuiteSparse_long_idd
#endif

/* For backward compatibility with prior versions of SuiteSparse.  The UF_*
 * macros are deprecated and will be removed in a future version. */
#ifndef UF_long
#define UF_long     SuiteSparse_long
#define UF_long_max SuiteSparse_long_max
#define UF_long_idd SuiteSparse_long_idd
#define UF_long_id  SuiteSparse_long_id
#endif

/* ========================================================================== */
/* === SuiteSparse_config parameters and functions ========================== */
/* ========================================================================== */

/* SuiteSparse-wide parameters will be placed in this struct. */

typedef struct SuiteSparse_config_struct
{
    void *(*malloc_memory) (size_t) ;           /* pointer to malloc */
    void *(*realloc_memory) (void *, size_t) ;  /* pointer to realloc */
    void (*free_memory) (void *) ;              /* pointer to free */
    void *(*calloc_memory) (size_t, size_t) ;   /* pointer to calloc */

} SuiteSparse_config ;

void *SuiteSparse_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc (>=1 is enforced) */
    size_t size_of_item,    /* sizeof each item */
    int *ok,                /* TRUE if successful, FALSE otherwise */
    SuiteSparse_config *config        /* SuiteSparse-wide configuration */
) ;

void *SuiteSparse_free      /* always returns NULL */
(
    void *p,                /* block to free */
    SuiteSparse_config *config        /* SuiteSparse-wide configuration */
) ;

void SuiteSparse_tic    /* start the timer */
(
    double tic [2]      /* output, contents undefined on input */
) ;

double SuiteSparse_toc  /* return time in seconds since last tic */
(
    double tic [2]      /* input: from last call to SuiteSparse_tic */
) ;

double SuiteSparse_time  /* returns current wall clock time in seconds */
(
    void
) ;

/* determine which timer to use, if any */
#ifndef NTIMER
#ifdef _POSIX_C_SOURCE
#if    _POSIX_C_SOURCE >= 199309L
#define SUITESPARSE_TIMER_ENABLED
#endif
#endif
#endif

/* ========================================================================== */
/* === SuiteSparse version ================================================== */
/* ========================================================================== */

/* SuiteSparse is not a package itself, but a collection of packages, some of
 * which must be used together (UMFPACK requires AMD, CHOLMOD requires AMD,
 * COLAMD, CAMD, and CCOLAMD, etc).  A version number is provided here for the
 * collection itself.  The versions of packages within each version of
 * SuiteSparse are meant to work together.  Combining one packge from one
 * version of SuiteSparse, with another package from another version of
 * SuiteSparse, may or may not work.
 *
 * SuiteSparse contains the following packages:
 *
 *  SuiteSparse_config version 4.0.2 (version always the same as SuiteSparse)
 *  AMD             version 2.3.1
 *  BTF             version 1.2.0
 *  CAMD            version 2.3.1
 *  CCOLAMD         version 2.8.0
 *  CHOLMOD         version 2.0.1
 *  COLAMD          version 2.8.0
 *  CSparse         version 3.1.1
 *  CXSparse        version 3.1.1
 *  KLU             version 1.2.1
 *  LDL             version 2.1.0
 *  RBio            version 2.1.1
 *  SPQR            version 1.3.1 (full name is SuiteSparseQR)
 *  UMFPACK         version 5.6.1
 *  MATLAB_Tools    various packages & M-files
 *
 * Other package dependencies:
 *  BLAS            required by CHOLMOD and UMFPACK
 *  LAPACK          required by CHOLMOD
 *  METIS 4.0.1     required by CHOLMOD (optional) and KLU (optional)
 */

#define SUITESPARSE_DATE "July 17, 2012"
#define SUITESPARSE_VER_CODE(main,sub) ((main) * 1000 + (sub))
#define SUITESPARSE_MAIN_VERSION 4
#define SUITESPARSE_SUB_VERSION 0
#define SUITESPARSE_SUBSUB_VERSION 2
#define SUITESPARSE_VERSION \
    SUITESPARSE_VER_CODE(SUITESPARSE_MAIN_VERSION,SUITESPARSE_SUB_VERSION)

#ifdef __cplusplus
}
#endif
#endif
