/** @file generic.h
 ** @brief Generic (@ref generic)
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_GENERIC_H
#define VL_GENERIC_H

#include "host.h"
#include "random.h"

#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <assert.h>

#if defined(VL_OS_WIN)
#include <Windows.h>
#endif

#if ! defined(VL_DISABLE_THREADS) && defined(VL_THREADS_POSIX)
#include <pthread.h>
#endif

/** @brief Library version string */
#define VL_VERSION_STRING "0.9.16"

/** @brief Maximum length (in characters) of an error message */
#define VL_ERR_MSG_LEN 1024

/** @name Type identidifers for atomic data types
 ** @{ */

#define VL_TYPE_FLOAT   1     /**< @c float type */
#define VL_TYPE_DOUBLE  2     /**< @c double type */
#define VL_TYPE_INT8    3     /**< @c ::vl_int8 type */
#define VL_TYPE_UINT8   4     /**< @c ::vl_uint8 type */
#define VL_TYPE_INT16   5     /**< @c ::vl_int16 type */
#define VL_TYPE_UINT16  6     /**< @c ::vl_uint16 type */
#define VL_TYPE_INT32   7     /**< @c ::vl_int32 type */
#define VL_TYPE_UINT32  8     /**< @c ::vl_uint32 type */
#define VL_TYPE_INT64   9     /**< @c ::vl_int64 type */
#define VL_TYPE_UINT64  10    /**< @c ::vl_uint64 type */

typedef vl_uint32 vl_type ;

/** @brief Get the name of a data type.
 ** @param type data type.
 ** @return data name of the data type.
 **
 ** @c type is one of ::VL_TYPE_FLOAT, ::VL_TYPE_DOUBLE,
 ** ::VL_TYPE_INT8, ::VL_TYPE_INT16, ::VL_TYPE_INT32, ::VL_TYPE_INT64,
 ** ::VL_TYPE_UINT8, ::VL_TYPE_UINT16, ::VL_TYPE_UINT32, ::VL_TYPE_UINT64.
 **/

VL_INLINE char const *
vl_get_type_name (vl_type type)
{
  switch (type) {
    case VL_TYPE_FLOAT   : return "float"  ;
    case VL_TYPE_DOUBLE  : return "double" ;
    case VL_TYPE_INT8    : return "int8"   ;
    case VL_TYPE_INT16   : return "int16"  ;
    case VL_TYPE_INT32   : return "int32"  ;
    case VL_TYPE_INT64   : return "int64"  ;
    case VL_TYPE_UINT8   : return "int8"   ;
    case VL_TYPE_UINT16  : return "int16"  ;
    case VL_TYPE_UINT32  : return "int32"  ;
    case VL_TYPE_UINT64  : return "int64"  ;
    default: return NULL ;
  }
}

/** @brief Get data type size.
 ** @param type data type.
 ** @return size (in byte)
 **
 ** @c type is one of ::VL_TYPE_FLOAT, ::VL_TYPE_DOUBLE,
 ** ::VL_TYPE_INT8, ::VL_TYPE_INT16, ::VL_TYPE_INT32, ::VL_TYPE_INT64,
 ** ::VL_TYPE_UINT8, ::VL_TYPE_UINT16, ::VL_TYPE_UINT32, ::VL_TYPE_UINT64.
 **/

VL_INLINE vl_size
vl_get_type_size (vl_type type)
{
  vl_size dataSize = 0 ;
  switch (type) {
    case VL_TYPE_DOUBLE : dataSize = sizeof(double) ; break ;
    case VL_TYPE_FLOAT  : dataSize = sizeof(float) ; break ;
    case VL_TYPE_INT64  : case VL_TYPE_UINT64 : dataSize = sizeof(vl_int64) ; break ;
    case VL_TYPE_INT32  : case VL_TYPE_UINT32 : dataSize = sizeof(vl_int32) ; break ;
    case VL_TYPE_INT16  : case VL_TYPE_UINT16 : dataSize = sizeof(vl_int16) ; break ;
    case VL_TYPE_INT8   : case VL_TYPE_UINT8  : dataSize = sizeof(vl_int8)  ; break ;
    default:
      abort() ;
  }
  return dataSize ;
}

/** @} */

/** ------------------------------------------------------------------
 ** @name Library state and configuration
 ** @{ */

/** @internal @brief VLFeat thread state */
typedef struct _VlThreadSpecificState
{
  /* errors */
  int lastError ;
  char lastErrorMessage [VL_ERR_MSG_LEN] ;

  /* random number generator */
  VlRand rand ;

  /* time */
#if defined(VL_OS_WIN)
  LARGE_INTEGER ticFreq ;
  LARGE_INTEGER ticMark ;
#else
  clock_t ticMark ;
#endif
} VlThreadSpecificState ;

/** @internal @brief VLFeat global state */
typedef struct _VlState
{

#if ! defined(VL_DISABLE_THREADS)
#if   defined(VL_THREADS_POSIX)
  pthread_key_t threadKey ;
  pthread_mutex_t mutex ;
  pthread_t mutexOwner ;
  pthread_cond_t mutexCondition ;
  size_t mutexCount ;
#elif defined(VL_THREADS_WIN)
  DWORD tlsIndex ;
  CRITICAL_SECTION mutex ;
#endif
#else
  VlThreadSpecificState * threadState ;
#endif

  int   (*printf_func)  (char const * format, ...) ;
  void *(*malloc_func)  (size_t) ;
  void *(*realloc_func) (void*,size_t) ;
  void *(*calloc_func)  (size_t, size_t) ;
  void  (*free_func)    (void*) ;

#if defined(VL_ARCH_IX86) || defined(VL_ARCH_X64) || defined(VL_ARCH_IA64)
  VlX86CpuInfo cpuInfo ;
#endif
  vl_size numCPUs ;

  vl_bool simdEnabled ;
  vl_size maxNumThreads ;

} VlState ;

/** @internal @brief VLFeat global state */
VL_EXPORT VlState _vl_state ;

VL_INLINE VlState * vl_get_state () ;
VL_INLINE VlThreadSpecificState * vl_get_thread_specific_state () ;
VL_EXPORT void vl_lock_state () ;
VL_EXPORT void vl_unlock_state () ;
VL_EXPORT VlThreadSpecificState * vl_thread_specific_state_new () ;
VL_EXPORT void vl_thread_specific_state_delete (VlThreadSpecificState * self) ;
VL_EXPORT char const * vl_get_version_string () ;
VL_EXPORT char * vl_configuration_to_string_copy () ;
VL_INLINE void vl_set_simd_enabled (vl_bool x) ;
VL_INLINE vl_bool vl_get_simd_enabled () ;
VL_INLINE vl_bool vl_cpu_has_sse3 () ;
VL_INLINE vl_bool vl_cpu_has_sse2 () ;
VL_INLINE vl_size vl_get_num_cpus () ;
VL_EXPORT VlRand * vl_get_rand () ;

/** @} */

/** ------------------------------------------------------------------
 ** @name Error handling
 ** @{ */

#define VL_ERR_OK       0  /**< No error */
#define VL_ERR_OVERFLOW 1  /**< Buffer overflow error */
#define VL_ERR_ALLOC    2  /**< Resource allocation error */
#define VL_ERR_BAD_ARG  3  /**< Bad argument or illegal data error */
#define VL_ERR_IO       4  /**< Input/output error */
#define VL_ERR_EOF      5  /**< End-of-file or end-of-sequence error */
#define VL_ERR_NO_MORE  5  /**< End-of-sequence @deprecated */

VL_INLINE int vl_get_last_error () ;
VL_INLINE char const *  vl_get_last_error_message () ;
VL_EXPORT int vl_set_last_error (int error, char const * errorMessage, ...) ;

/** @} */

/** ------------------------------------------------------------------
 ** @name Memory allocation
 ** @{ */

VL_EXPORT void
vl_set_alloc_func (void *(*malloc_func)  (size_t),
                   void *(*realloc_func) (void*,size_t),
                   void *(*calloc_func)  (size_t, size_t),
                   void  (*free_func)    (void*)) ;
VL_INLINE void *vl_malloc  (size_t n) ;
VL_INLINE void *vl_realloc (void *ptr, size_t n) ;
VL_INLINE void *vl_calloc  (size_t n, size_t size) ;
VL_INLINE void  vl_free    (void* ptr) ;

/** @} */

/** ------------------------------------------------------------------
 ** @name Logging
 ** @{ */

/** ------------------------------------------------------------------
 ** @brief Customizable printf function pointer type */
typedef int(*printf_func_t) (char const *format, ...) ;

/** @brief Set printf function
 ** @param printf_func  pointer to @c printf.
 ** Let @c print_func be NULL to disable printf.
 **/
VL_EXPORT void vl_set_printf_func (printf_func_t printf_func) ;

/** @def VL_PRINTF
 ** @brief Call user-customizable @c printf function
 **
 ** The function calls the user customizable @c printf.
 **/
#define VL_PRINTF (*vl_get_state()->printf_func)

/** @def VL_PRINT
 ** @brief Same as ::VL_PRINTF (legacy code)
 **/
#define VL_PRINT (*vl_get_state()->printf_func)

/** @} */

/** ------------------------------------------------------------------
 ** @name Common operations
 ** @{ */

/** @brief Min operation
 ** @param x value
 ** @param y value
 ** @return the minimum of @a x and @a y.
 **/
#define VL_MIN(x,y) (((x)<(y))?(x):(y))

/** @brief Max operation
 ** @param x value.
 ** @param y value.
 ** @return the maximum of @a x and @a y.
 **/
#define VL_MAX(x,y) (((x)>(y))?(x):(y))

/** @brief Signed left shift operation
 **
 ** The macro is equivalent to the builtin @c << operator, but it
 ** supports negative shifts too.
 **
 ** @param x value.
 ** @param n number of shift positions.
 ** @return @c x << n .
 **/
#define VL_SHIFT_LEFT(x,n) (((n)>=0)?((x)<<(n)):((x)>>-(n)))
/* @} */

/** ------------------------------------------------------------------
 ** @name Measuring time
 ** @{
 **/
VL_EXPORT void vl_tic () ;
VL_EXPORT double vl_toc () ;
VL_EXPORT double vl_get_cpu_time () ;
/** @} */

VL_EXPORT void vl_constructor () ;
VL_EXPORT void vl_destructor () ;


/* -------------------------------------------------------------------
 *                                                    Inline functions
 * ---------------------------------------------------------------- */



VL_INLINE VlState *
vl_get_state ()
{
  return &_vl_state ;
}

VL_INLINE VlThreadSpecificState *
vl_get_thread_specific_state ()
{
#ifdef VL_DISABLE_THREADS
  return vl_get_state()->threadState ;
#else
  VlState * state ;
  VlThreadSpecificState * threadState ;

  vl_lock_state() ;
  state = vl_get_state() ;

#if defined(VL_THREADS_POSIX)
  threadState = (VlThreadSpecificState *) pthread_getspecific(state->threadKey) ;
#elif defined(VL_THREADS_WIN)
  threadState = (VlThreadSpecificState *) TlsGetValue(state->tlsIndex) ;
#endif

  if (! threadState) {
    threadState = vl_thread_specific_state_new () ;
  }

#if defined(VL_THREADS_POSIX)
  pthread_setspecific(state->threadKey, threadState) ;
#elif defined(VL_THREADS_WIN)
  TlsSetValue(state->tlsIndex, threadState) ;
#endif

  vl_unlock_state() ;
  return threadState ;
#endif
}

VL_INLINE void
vl_set_simd_enabled (vl_bool x)
{
  vl_get_state()->simdEnabled = x ;
}

VL_INLINE vl_bool
vl_get_simd_enabled ()
{
  return vl_get_state()->simdEnabled ;
}

VL_INLINE vl_bool
vl_cpu_has_sse3 ()
{
#if defined(VL_ARCH_IX86) || defined(VL_ARCH_X64) || defined(VL_ARCH_IA64)
  return vl_get_state()->cpuInfo.hasSSE3 ;
#else
  return 0 ;
#endif
}

VL_INLINE vl_bool
vl_cpu_has_sse2 ()
{
#if defined(VL_ARCH_IX86) || defined(VL_ARCH_X64) || defined(VL_ARCH_IA64)
  return vl_get_state()->cpuInfo.hasSSE2 ;
#else
  return 0 ;
#endif
}

VL_INLINE vl_size
vl_get_num_cpus ()
{
  return vl_get_state()->numCPUs ;
}

VL_INLINE int
vl_get_last_error () {
  return vl_get_thread_specific_state()->lastError ;
}

VL_INLINE char const *
vl_get_last_error_message ()
{
  return vl_get_thread_specific_state()->lastErrorMessage ;
}

VL_INLINE void*
vl_malloc (size_t n)
{
  return (vl_get_state()->malloc_func)(n) ;
}

VL_INLINE void*
vl_realloc (void* ptr, size_t n)
{
  return (vl_get_state()->realloc_func)(ptr, n) ;
}

VL_INLINE void*
vl_calloc (size_t n, size_t size)
{
  return (vl_get_state()->calloc_func)(n, size) ;
}

VL_INLINE void
vl_free (void *ptr)
{
  (vl_get_state()->free_func)(ptr) ;
}

/* VL_GENERIC_H */
#endif
