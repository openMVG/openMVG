/** @file generic.c
 ** @brief Generic - Definition
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

/**
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@mainpage VLFeat -- Vision Lab Features Library
@version __VLFEAT_VERSION__
@author The VLFeat Team
@par Copyright &copy; 2007-12 Andrea Vedaldi and Brian Fulkerson
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<em>VLFeat C library contains implementations of common computer
vision algorithms, with a special focus on visual features for
matching image regions. Applications include structure from motion and
object and category detection and recognition.

We strive to make the library free of clutter, portable (VLFeat is
largely C-89 compatible), and self- documented. Different parts of the
library are weakly interdependent, simplifying understanding and
extraction of code.</em>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section main-contents Contents
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

- <b>Algorithms</b>
  - @ref sift
  - @ref dsift
  - @ref mser
  - @ref kmeans
  - @ref ikmeans.h  "Integer K-means (IKM)"
  - @ref hikmeans.h "Hierarchical Integer K-means (HIKM)"
  - @ref aib
  - @ref kdtree
  - @ref homkermap
  - @ref pegasos
  - @ref slic
  - @ref quickshift
  - @ref hog
  - @ref covdet
- <b>Support functionalities</b>
  - @ref host.h      "Platform abstraction"
  - @ref generic
  - @ref random.h    "Random number generator"
  - @ref mathop.h    "Math operations"
  - @ref heap-def.h  "Generic heap object (priority queue)"
  - @ref stringop.h  "String operations"
  - @ref imopv.h     "Image operations"
  - @ref pgm.h       "PGM reading and writing"
  - @ref rodrigues.h "Rodrigues formula"
  - @ref mexutils.h  "MATLAB MEX helper functions"
  - @ref getopt_long.h "Drop-in @c getopt_long replacement"
- @ref design
  - @ref design-objects     "Objects"
  - @ref design-resources   "Memory and resource management"
  - @ref design-threads     "Multi-threading"
  - @ref design-portability "Portability"
- @ref dev
- @ref main-glossary

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@page design VLFeat design concepts
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat is designed to be portable and simple to integrate with high
level languages such as MATLAB. This section illustrates and
motivates the aspects of the design that are relevant to the users of
the library.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section design-resources Memory and resource handling
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Some VLFeat functions return pointers to memory blocks or
objects. Only ::vl_malloc, ::vl_calloc, ::vl_realloc, or functions
whose name contains either the keywords @c new or @c copy,
transfer the ownership of the memory block or object to the
caller. The caller must dispose explicitly of all the resources he
owns (by calling ::vl_free for a memory block, or the appropriate
destructor for an object).

The memory allocation functions can be customized by
::vl_set_alloc_func (which sets the implementations of ::vl_malloc,
::vl_realloc, ::vl_calloc and ::vl_free). Remapping the memory
allocation functions can be done only if there are no currently
allocated VLFeat memory blocks or objects. The memory allocation
functions are common to all threads.

VLFeat uses three rules that simplify handling exceptions:

- The library allocates local memory only through the reprogrammable
  ::vl_malloc, ::vl_calloc, and ::vl_realloc functions.

- The only resource referenced by VLFeat objects is memory (for
  instance, it is illegal for an object to reference an open file).
  Other resources such as files or threads may be allocated within a
  VLFeat function call, but they are all released before the function
  ends, or their ownership is directly transferred to the caller.

- The global library state is an exception. It cannot reference any
  local object created by the caller and uses the standard C memory
  allocation functions.

In this way, the VLFeat local state can be reset at any time simply
by disposing of all the memory allocated so far. The latter can be
done easily by mapping the memory allocation functions to
implementations that track the memory blocks allocated, and simply
disposing of all such blocks. Since the global state does not
reference any local object nor uses the remapped memory functions, it
is unaffected by such an operation; conversely, since no VLFeat
object references anything but memory, this guarantees that all
allocated resources are properly disposed (no leaks). This is used
extensively in the design of MATLAB MEX files (see @ref
design-matlab).

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section design-objects Objects
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Many VLFeat algorithms are available in the form of
&ldquo;objects&rdquo;. Notice that the C language, used by VLFeat,
does not support objects explicitly. Here an object indicates an
opaque data structure along with a number of functions (methods)
operationg on it.

Object names are capitalised and start with the <code>Vl</code>
prefix (for example ::VlSiftFilt). Object methods are lowercase and
start with the <code>vl_<object_name>_</code> suffix
(e.g. ::vl_sift_new). Object methods typically include a constructor
(e.g. ::vl_sift_new), a destructor (::vl_sift_delete), some getter
methods (::vl_sift_get_octave_index), and some setter methods
(::vl_sift_set_magnif).

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section design-threads Multi-threading
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat can be used with multiple threads by treating appropriately
different kinds of operations:

- <b>Local operations (reentrancy).</b> Most VLFeat operations are
  reentrant, in the sense that they operate on local data. It is safe
  to execute such operations concurrently from multiple threads as
  long as each thread operates on an independent sets of objects or
  memory blocks.  By contrast, operating on the same object or memory
  block from multiple threads requires proper synchronization by the
  user.

- <b>Task-specific operations.</b> A few operations are intrinsically
  non-reentrant but thread-specific. These include: retrieving the
  last error by ::vl_get_last_error and obtaining the thread-specific
  random number generator by ::vl_get_rand. VLFeat makes such
  operations thread-safe by operating on task-specific data.

- <b>Global operations.</b> A small number of operations are
  non-reentrant <em>and</em> affect all threads simultaneously. These
  are restricted to changing certain global configuration parameters,
  such as the memory allocation functions by
  ::vl_set_alloc_func. These operations are <em>not</em> thread safe
  and should be executed before multiple threads are started.

Some VLFeat algorithms are randomised. Each thread has his own random
number generator (an instance of ::VlRand) accessed by
::vl_get_rand. To make calculations reproducible the random number
generator must be seeded appropriately in each thread. Notice also
that using the same VLFeat object from multiple threads (with
appropriate synchronization) will cause multiple random number
generators to be intermixed.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section design-portability Portability features
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Platform dependent details are isolated in the @ref host.h
library module. These include:

- Atomic types (e.g. ::vl_int32).
- Special syntaxes for the declaration of symbols exported by the library
  and inline functions (e.g. ::VL_EXPORT).
- Host-dependent conversion of data endianess
  (e.g. ::vl_swap_host_big_endianness_8()).

VLFeat uses processor specific features (e.g. Intel SSE) if those are
available at run time.

<!-- @see http://www.macresearch.org/how_to_properly_use_sse3_and_ssse3_and_future_intel_vector_extensions_0  -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section design-matlab MATLAB integration issues
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The VLFeat C library is designed to integrate seamlessly with MATLAB.
Binary compatibility is simplified by the use of the C language
(rather than C++). In addition, the library design follows certain
restrictions that make it compatible with the MATLAB MEX interface.

The main issue in calling a library function from a MATLAB MEX
function is that MATLAB can abort the execution of the MEX function
at any point, either due to an error, or directly upon a user request
(Ctrl-C) (empirically, however, a MEX function seems to be
interruptible only during the invocation of certain functions of the
MEX API such as @c mexErrMsgTxt).

When a MEX function is interrupted, resources (memory blocks or
objects) whose ownership was transferred from VLFeat to the MEX
function may be leaked. Notice that interrupting a MEX function would
similarly leak any memory block allocated within the MEX function. To
solve this issue, MATLAB provides his own memory manager (@c
mxMalloc, @c mxRealloc, ...). When a MEX file is interrupted or ends,
all memory blocks allocated by using one of such functions are
released, preventing leakage.

In order to integrate VLFeat with this model in the most seamless
way, VLFeat memory allocation functions (::vl_malloc, ::vl_realloc,
::vl_calloc) are mapped to the corresponding MEX memory allocation
functions. Such functions automatically dispose of all the memory
allocated by a MEX function when the function ends (even because of
an exception). Because of the restrictions of the library design
illustrated in @ref design-resources, this operation is safe and
correctly dispose of VLFeat local state. As a consequence, it is
possible to call @c mexErrMsgTxt at any point in the MEX function
without worrying about leaking resources.

This however comes at the price of some limitations. Beyond the
restrictions illustrated in @ref design-resources, here we note that no
VLFeat local resource (memory blocks or objects) can persist across
MEX file invocations. This implies that any result produced by a
VLFeat MEX function must be converted back to a MATLAB object such as
a vector or a structure. In particular, there is no direct way of
creating an object within a MEX file, returning it to MATLAB, and
passing it again to another MEX file.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section main-metaprogramming Preprocessor metaprogramming
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Part of VLFeat code uses a simple form of perprocessor metaprogramming.
This technique is used, similarly to C++ templates, to instantiate
multiple version of a given algorithm for different data types
(e.g. @c float and @c double).

In most cases preprocessor metaprogramming is invisible to the library
user, as it is used only internally.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@page main-glossary Glossary
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

 - <b>Column-major.</b> A <em>M x N </em> matrix <em>A</em> is
 stacked with column-major order as the sequence \f$(A_{11}, A_{21},
 \dots, A_{12}, \dots)\f$. More in general, when stacking a multi
 dimensional array this indicates that the first index is the one
 varying most quickly, with the other followed in the natural order.
 - <b>Opaque structure.</b> A structure is opaque if the user is not supposed
 to access its member directly, but through appropriate interface functions.
 Opaque structures are commonly used to define objects.
 - <b>Row-major.</b> A <em>M x N </em> matrix <em>A</em> is
 stacked with row-major order as the sequence \f$(A_{11}, A_{12},
 \dots, A_{21}, \dots)\f$. More in general, when stacking a multi
 dimensional array this indicates that the last index is the one
 varying most quickly, with the other followed in reverse order.
 - <b>Feature frame.</b> A <em>feature frame</em> is the geometrical
 description of a visual features. For instance, the frame of
 a @ref sift.h "SIFT feature" is oriented disk and the frame of
 @ref mser.h "MSER feature" is either a compact and connected set or
 a disk.
 - <b>Feature descriptor.</b> A <em>feature descriptor</em> is a quantity
 (usually a vector) which describes compactly the appearance of an
 image region (usually corresponding to a feature frame).
**/

/**

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@page dev Developing the library
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section dev-doc Coding style
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<ul>

<li><b>Develop a sensibility for good-looking code.</b> The coding
style should be as uniform as possible throughout the library. The
style is specified by this set of rules. However, the quality of the
code can only be guaranteed by a reasonable application of the rules
and combined with one's eye for good code.</li>

<li><b>No whitespaces at the end of lines.</b> Whitespaces introduce
invisible changes in the code that are however picked up by control
version systems such as Git.</li>

<li><b>Descriptive variable names.</b> Most variable names start with
a lower case letter and are capitalized, e.g., @c numElements. Only
the following abbreviations are considered acceptable: @c num. The @c
dimension of a vector is the number of elements it contains (for other
objects that could be a @c size, a @c length, or a @c
numElements). For multi-dimensional arrays, @c dimensions could
indicate the array with each of the @c numDimensions dimensions.</li>

<li><b>Short variable names.</b> For indexes in short for loops it is
fine to use short index names such as @c i, @c j, and @c k. For example:

@code
for (i = 0 ; i < numEntries ; ++i) values[i] ++ ;
@endcode

is considered acceptable.</li>

</ul>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section dev-doc Documenting the code
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

The VLFeat C library code contains its own in documentation <a
href='http://www.stack.nl/~dimitri/doxygen/'>Doxygen</a> format. The
documentation consists in generic pages, such as the @ref index
"index" and the page you are reading, and documentations for each
specified library module, usually corresponding to a certain header
file.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection devl-doc-modules Documenting the library modules
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

A library module groups a number of data types and functions that
implement a certain functionality of VLFeat. Consider a module called
<em>Example Module</em>. Then one would typically have:

<ul>
<li>A header or declaration file @c example-module.h. Such a file has an
heading of the type:

@verbinclude example-module-doc.h

This comment block contains a file directive, causing the file to be
included in the documentation, a brief directive, specifying a short
description of what the file is, and a list of authors. A
(non-Doxygen) comment block with a short the copyright notice follows.
The brief directive should include a <code>@@ref</code> directive to point
to the main documentation page describing the module, if there is one.
</li>

<li> An implementation or definition file @c example-module.c. This file
has an heading of the type:

@verbinclude example-module-doc.c

This is similar to the declaration file, except for the content of the
brief comment.
</li>
</ul>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection devl-doc-functions Documenting the library functions
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@subsection devl-doc-bib Bibliographic references
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Since version 0.9.14, the VLFeat C library documentation makes use of
a proper bibliographic reference in BibTeX format (see the file @c
docsrc/vlfeat.bib). Doxygen uses this file when it sees instances of
the <code>@@cite{xyz}</code> command.  Here @c xyz is a BibTeX
key. For example, @c vlfeat.bib file contains the entry:

<pre>
@@inproceedings{martin97the-det-curve,
	Author = {A. Martin and G. Doddington and T. Kamm and M. Ordowski and M. Przybocki},
	Booktitle = {Proc. Conf. on Speech Communication and Technology},
	Title = {The {DET} curve in assessment of detection task performance},
	Year = {1997}}
</pre>

For example, the Doxygen directive
<code>@@cite{martin97the-det-curve}</code> generates the output
@cite{martin97the-det-curve}, which is a link to the corresponding
entry in the bibliography.

**/

/**

@file generic.h

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@page generic Preprocessor, library state, etc.
@author Andrea Vedaldi
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

@ref generic.h provides the following functionalities:

- @ref generic-preproc
- @ref generic-state
- @ref generic-error
- @ref generic-memory
- @ref generic-logging
- @ref generic-time

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-preproc C preprocessor helpers
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat provides a few C preprocessor macros of general
utility. These include stringification (::VL_STRINGIFY,
::VL_XSTRINGIFY) and concatenation (::VL_CAT, ::VL_XCAT) of symbols.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-state VLFeat state and configuration parameters
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat has some global configuration parameters that can
changed. Changing the configuration is thread unsave
(@ref design-threads). Use ::vl_set_simd_enabled to toggle the use of
a SIMD unit (Intel SSE code), ::vl_set_alloc_func to change
the memory allocation functions, and ::vl_set_printf_func
to change the logging function.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-error Error handling
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

Some VLFeat functions signal errors in a way similar to the
standard C library. In case of error, a VLFeat function
may return an error code directly,
or an invalid result (for instance a negative file descriptor or a null
pointer). Then ::vl_get_last_error and ::vl_get_last_error_message can be used
to retrieve further details about the error (these functions should be
used right after an error has occurred, before any other VLFeat call).

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-memory Memory allocation
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat uses the ::vl_malloc, ::vl_realloc, ::vl_calloc and ::vl_free
functions to allocate memory. Normally these functions are mapped to
the underlying standard C library implementations. However
::vl_set_alloc_func can be used to map them to other
implementations.  For instance, in MATALB MEX files these functions
are mapped to the MATLAB equivalent which has a garbage collection
mechanism to cope with interruptions during execution.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-logging Logging
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat uses the macros ::VL_PRINT and ::VL_PRINTF to print progress
or debug informations. These functions are normally mapped to the @c
printf function of the underlying standard C library. However
::vl_set_printf_func can be used to map it to a different
implementation. For instance, in MATLAB MEX files this function is
mapped to @c mexPrintf. Setting the function to @c NULL disables
logging.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->
@section generic-time Measuring time
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  -->

VLFeat provides ::vl_tic and ::vl_toc as an easy way of measuring
elapsed time.

**/

#include "generic.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#if defined(VL_OS_WIN)
#include <Windows.h>
#elif defined(_POSIX_THREADS)
#include <pthread.h>
#endif

#if defined(VL_OS_MACOSX) || defined(VL_OS_LINUX)
#include <unistd.h>
#endif

/** ------------------------------------------------------------------
 ** @brief Get version string
 ** @return library version string
 **/

VL_EXPORT char const *
vl_get_version_string ()
{
  return VL_VERSION_STRING ;
}

/** ------------------------------------------------------------------
 ** @brief Human readable library configuration
 ** @return a new string with the library configuration.
 **
 ** The function returns a new string with a human readable
 ** rendition of the library configuration.
 **/

VL_EXPORT char *
vl_configuration_to_string_copy ()
{
  char * string = 0 ;
  int length = 0 ;
  char * staticString = vl_static_configuration_to_string_copy() ;
  char * cpuString =
#if defined(VL_ARCH_IX86) || defined(VL_ARCH_X64) || defined(VL_ARCH_IA64)
  _vl_x86cpu_info_to_string_copy(&vl_get_state()->cpuInfo) ;
#else
  "Generic CPU" ;
#endif
#if defined(DEBUG)
  int const debug = 1 ;
#else
  int const debug = 0 ;
#endif

  while (string == 0) {
    if (length > 0) {
      string = vl_malloc(sizeof(char) * length) ;
      if (string == NULL) break ;
    }
    length = snprintf(string, length,
                      "VLFeat version %s\n"
                      "    Static config: %s\n"
                      "    %" VL_FMT_SIZE " CPU(s): %s\n"
                      "    Debug: %s\n",
                      vl_get_version_string (),
                      staticString,
                      vl_get_num_cpus(), cpuString,
                      VL_YESNO(debug)) ;
    length += 1 ;
  }

  if (staticString) vl_free(staticString) ;
  if (cpuString) vl_free(cpuString) ;
  return string ;
}

/** @internal @brief A printf that does not do anything */
static int
do_nothing_printf (char const* format VL_UNUSED, ...)
{
  return 0 ;
}

/** --------------------------------------------------------------- */

VlState _vl_state ;

/** ------------------------------------------------------------------
 ** @internal @brief Lock VLFeat state
 **
 ** The function locks VLFeat global state mutex.
 **
 ** The mutex is recursive: locking multiple times from the same thread
 ** is a valid operations, but requires an equivalent number
 ** of calls to ::vl_unlock_state.
 **
 ** @sa ::vl_unlock_state
 **/

VL_EXPORT void
vl_lock_state ()
{
#if ! defined(VL_DISABLE_THREADS)
#if   defined(VL_THREADS_POSIX)
  VlState * state = vl_get_state () ;
  pthread_t thisThread = pthread_self () ;
  pthread_mutex_lock (&state->mutex) ;
  if (state->mutexCount >= 1 &&
      pthread_equal (state->mutexOwner, thisThread)) {
    state->mutexCount ++ ;
  } else {
    while (state->mutexCount >= 1) {
      pthread_cond_wait (&state->mutexCondition, &state->mutex) ;
    }
    state->mutexOwner = thisThread ;
    state->mutexCount = 1 ;
  }
  pthread_mutex_unlock (&state->mutex) ;
#elif defined(VL_THREADS_WIN)
  EnterCriticalSection (&vl_get_state()->mutex) ;
#endif
#endif
}

/** ------------------------------------------------------------------
 ** @internal @brief Unlock VLFeat state
 **
 ** The function unlocks VLFeat global state mutex.
 **
 ** @sa ::vl_lock_state
 **/

VL_EXPORT void
vl_unlock_state ()
{
#if ! defined(VL_DISABLE_THREADS)
#if   defined(VL_THREADS_POSIX)
  VlState * state = vl_get_state () ;
  pthread_mutex_lock (&state->mutex) ;
  -- state->mutexCount ;
  if (state->mutexCount == 0) {
    pthread_cond_signal (&state->mutexCondition) ;
  }
  pthread_mutex_unlock (&state->mutex) ;
#elif defined(VL_THREADS_WIN)
  LeaveCriticalSection (&vl_get_state()->mutex) ;
#endif
#endif
}

/** @internal @fn ::vl_get_state()
 ** @brief Return VLFeat global state
 **
 ** The function returns a pointer to VLFeat global state.
 **
 ** @return pointer to the global state structure.
 **/

/** @internal @fn ::vl_get_thread_specific_state()
 ** @brief Get VLFeat thread state
 **
 ** The function returns a pointer to VLFeat thread state.
 **
 ** @return pointer to the thread state structure.
 **/

/** @fn ::vl_malloc(size_t)
 ** @brief Call customizable @c malloc function
 ** @param n number of bytes to allocate.
 **
 ** The function calls the user customizable @c malloc.
 **
 ** @return result of @c malloc
 **/

/** @fn ::vl_realloc(void*,size_t)
 ** @brief Call customizable @c resize function
 **
 ** @param ptr buffer to reallocate.
 ** @param n   number of bytes to allocate.
 **
 ** The function calls the user-customizable @c realloc.
 **
 ** @return result of the user-customizable @c realloc.
 **/

/** @fn ::vl_calloc(size_t,size_t)
 ** @brief Call customizable @c calloc function
 **
 ** @param n    size of each element in byte.
 ** @param size size of the array to allocate (number of elements).
 **
 ** The function calls the user-customizable @c calloc.
 **
 ** @return result of the user-customizable @c calloc.
 **/

/** @fn ::vl_free(void*)
 ** @brief Call customizable @c free function
 **
 ** @param ptr buffer to free.
 **
 ** The function calls the user customizable @c free.
 **/

/** @fn ::vl_get_last_error()
 ** @brief Get VLFeat last error code
 **
 ** The function returns the code of the last error generated
 ** by VLFeat.
 **
 ** @return laste error code.
 **/

/** @internal @fn ::vl_get_last_error_message()
 ** @brief Get VLFeat last error message
 **
 ** The function returns the message of the last
 ** error generated by VLFeat.
 **
 ** @return last error message.
 **/

/** @fn ::vl_set_simd_enabled(vl_bool)
 ** @brief Toggle usage of SIMD instructions
 ** @param x @c true if SIMD instructions are used.
 **
 ** Notice that SIMD instructions are used only if the CPU model
 ** supports them. Note also that data alignment may restrict the use
 ** of such instructions.
 **
 ** @see ::vl_cpu_has_sse2(), ::vl_cpu_has_sse3(), etc.
 **/

/** @fn ::vl_get_simd_enabled()
 ** @brief Are SIMD instructons enabled?
 ** @return @c true is SIMD instructions are enabled.
 **/

/** @fn ::vl_cpu_has_sse3()
 ** @brief Check for SSE3 instruction set
 ** @return @c true if SSE3 is present.
 **/

/** @fn ::vl_cpu_has_sse2()
 ** @brief Check for SSE2 instruction set
 ** @return @c true if SSE2 is present.
 **/

/** ------------------------------------------------------------------
 ** @internal @brief Set last VLFeat error
 **
 ** The function sets the code and optionally the error message
 ** of the last encountered error. @a errorMessage is the message
 ** format. It uses the @c printf convention and is followed by
 ** the format arguments. The maximum lenght of the error message is
 ** given by ::VL_ERR_MSG_LEN (longer messages are truncated).
 **
 ** Passing @c NULL as @a errorMessage
 ** sets the error message to the empty string.
 **
 ** @param error error code.
 ** @param errorMessage error message format string.
 ** @param ... format string arguments.
 ** @return error code.
 **/

VL_EXPORT int
vl_set_last_error (int error, char const * errorMessage, ...)
{
  VlThreadSpecificState * state = vl_get_thread_specific_state() ;
  va_list args;
  va_start(args, errorMessage) ;
  if (errorMessage) {
#ifdef VL_COMPILER_LCC
    vsprintf(state->lastErrorMessage, errorMessage, args) ;
#else
    vsnprintf(state->lastErrorMessage,
              sizeof(state->lastErrorMessage)/sizeof(char),
              errorMessage, args) ;
#endif
  } else {
    state->lastErrorMessage[0] = 0 ;
  }
  state->lastError = error ;
  va_end(args) ;
  return error ;
}

/** ------------------------------------------------------------------
 ** @brief Set memory allocation functions
 ** @param malloc_func  pointer to @c malloc.
 ** @param realloc_func pointer to @c realloc.
 ** @param calloc_func  pointer to @c calloc.
 ** @param free_func    pointer to @c free.
 **/

VL_EXPORT void
vl_set_alloc_func (void *(*malloc_func)  (size_t),
                   void *(*realloc_func) (void*, size_t),
                   void *(*calloc_func)  (size_t, size_t),
                   void  (*free_func)    (void*))
{
  VlState * state ;
  vl_lock_state () ;
  state = vl_get_state() ;
  state->malloc_func  = malloc_func ;
  state->realloc_func = realloc_func ;
  state->calloc_func  = calloc_func ;
  state->free_func    = free_func ;
  vl_unlock_state () ;
}

VL_EXPORT void
vl_set_printf_func (printf_func_t printf_func)
{
  vl_get_state()->printf_func = printf_func ? printf_func : do_nothing_printf ;
}


/** ------------------------------------------------------------------
 ** @brief Get processor time
 ** @return processor time.
 ** @sa ::vl_tic, ::vl_toc
 **/

VL_EXPORT double
vl_get_cpu_time ()
{
  #ifdef VL_OS_WIN
  VlThreadSpecificState * threadState = vl_get_thread_specific_state() ;
  LARGE_INTEGER mark ;
  QueryPerformanceCounter (&mark) ;
  return (double)mark.QuadPart / (double)threadState->ticFreq.QuadPart ;
#else
  return (double)clock() / (double)CLOCKS_PER_SEC ;
#endif
}

/** ------------------------------------------------------------------
 ** @brief Reset processor time reference
 ** The function resets VLFeat TIC/TOC time reference.
 ** @sa ::vl_get_cpu_time, ::vl_toc.
 **/

VL_EXPORT void
vl_tic ()
{
  VlThreadSpecificState * threadState = vl_get_thread_specific_state() ;
#ifdef VL_OS_WIN
  QueryPerformanceCounter (&threadState->ticMark) ;
#else
  threadState->ticMark = clock() ;
#endif
}

/** ------------------------------------------------------------------
 ** @brief Get elapsed time since tic
 **
 ** The function
 ** returns the processor time elapsed since ::vl_tic was called last.
 **
 ** @remark In multi-threaded applications, there is an independent
 ** timer for each execution thread.
 **
 ** @remark On UNIX, this function uses the @c clock() system call.
 ** On Windows, it uses the @c QueryPerformanceCounter() system call,
 ** which is more accurate than @c clock() on this platform.
 **
 ** @return elapsed time in seconds.
 **/

VL_EXPORT double
vl_toc ()
{
  VlThreadSpecificState * threadState = vl_get_thread_specific_state() ;
#ifdef VL_OS_WIN
  LARGE_INTEGER tocMark ;
  QueryPerformanceCounter(&tocMark) ;
  return (double) (tocMark.QuadPart - threadState->ticMark.QuadPart) /
    threadState->ticFreq.QuadPart ;
#else
  return (double) (clock() - threadState->ticMark) / CLOCKS_PER_SEC ;
#endif
}

/** ------------------------------------------------------------------
 ** @brief Get the random number generator for this thread
 ** @return random number generator.
 **
 ** The function returns a pointer to the random number genrator
 ** for this thread.
 **/

VL_EXPORT VlRand *
vl_get_rand ()
{
  return &vl_get_thread_specific_state()->rand ;
}

/* -------------------------------------------------------------------
 *                       Library construction and destruction routines
 *  --------------------------------------------------------------- */

VL_EXPORT VlThreadSpecificState *
vl_thread_specific_state_new ()
{
  VlThreadSpecificState * self ;
#if defined(DEBUG)
  printf("VLFeat thread constructor called\n") ;
#endif
  self = malloc(sizeof(VlThreadSpecificState)) ;
  self->lastError = 0 ;
  self->lastErrorMessage[0] = 0 ;
#if defined(VL_OS_WIN)
  QueryPerformanceFrequency (&self->ticFreq) ;
  self->ticMark.QuadPart = 0 ;
#else
  self->ticMark = 0 ;
#endif
  vl_rand_init (&self->rand) ;

  return self ;
}

VL_EXPORT void
vl_thread_specific_state_delete (VlThreadSpecificState * self)
{
#if defined(DEBUG)
  printf("VLFeat thread destructor called\n") ;
#endif
  free (self) ;
}

#if (defined(VL_OS_LINUX) || defined(VL_OS_MACOSX)) && defined(VL_COMPILER_GNUC)

//static void vl_constructor () __attribute__ ((constructor)) ;
//static void vl_destructor () __attribute__ ((destructor))  ;

#elif defined(VL_OS_WIN)

//static void vl_constructor () ;
//static void vl_destructor () ;

BOOL WINAPI DllMain(
    HINSTANCE hinstDLL,  // handle to DLL module
    DWORD fdwReason,     // reason for calling function
    LPVOID lpReserved )  // reserved
{
  VlState * state ;
  VlThreadSpecificState * threadState ;
  switch (fdwReason) {
    case DLL_PROCESS_ATTACH:
      /* Initialize once for each new process */
      vl_constructor () ;
      break ;

    case DLL_THREAD_ATTACH:
      /* Do thread-specific initialization */
      break ;

    case DLL_THREAD_DETACH:
      /* Do thread-specific cleanup */
#if ! defined(VL_DISABLE_THREADS) && defined(VL_THREADS_WIN)
      state = vl_get_state() ;
      threadState = (VlThreadSpecificState*) TlsGetValue(state->tlsIndex) ;
      if (threadState) {
        vl_thread_specific_state_delete (threadState) ;
      }
#endif
      break;

    case DLL_PROCESS_DETACH:
      /* Perform any necessary cleanup */
      vl_destructor () ;
      break;
    }
    return TRUE ; /* Successful DLL_PROCESS_ATTACH */
}

#endif

/** @internal @brief Initialize VLFeat */
void
vl_constructor ()
{
  VlState * state ;
#if defined(DEBUG)
  printf("VLFeat DEBUG: constructor begins.\n") ;
#endif

  state = vl_get_state() ;

#if ! defined(VL_DISABLE_THREADS)
#if defined(DEBUG)
  printf("VLFeat DEBUG: constructing thread specific state.\n") ;
#endif
#if   defined(VL_THREADS_POSIX)
  {
    typedef void (*destructorType)(void * );
    pthread_key_create (&state->threadKey,
                        (destructorType)
                          vl_thread_specific_state_delete) ;
    pthread_mutex_init (&state->mutex, NULL) ;
    pthread_cond_init (&state->mutexCondition, NULL) ;
  }
#elif defined(VL_THREADS_WIN)
  InitializeCriticalSection (&state->mutex) ;
  state->tlsIndex = TlsAlloc () ;
#endif
#else
  /* threading support disabled */
#if defined(DEBUG)
  printf("VLFeat DEBUG: constructing the generic thread state instance (threading support disabled).\n") ;
#endif
  vl_get_state()->threadState = vl_thread_specific_state_new() ;
#endif

  state->malloc_func  = malloc ;
  state->realloc_func = realloc ;
  state->calloc_func  = calloc ;
  state->free_func    = free ;
  state->printf_func  = printf ;

#if defined(VL_ARCH_IX86) || defined(VL_ARCH_X64) || defined(VL_ARCH_IA64)
  _vl_x86cpu_info_init (&state->cpuInfo) ;
#endif

#if defined(VL_OS_WIN)
  {
    SYSTEM_INFO info;
    GetSystemInfo (&info) ;
    state->numCPUs = info.dwNumberOfProcessors ;
  }
#elif defined(VL_OS_MACOSX) || defined(VL_OS_LINUX)
  state->numCPUs = sysconf(_SC_NPROCESSORS_ONLN) ;
#else
  state->numCPUs = 1 ;
#endif
  state->simdEnabled = VL_TRUE ;
  state->maxNumThreads = 1 ;
#if defined(DEBUG)
  printf("VLFeat DEBUG: constructor ends.\n") ;
#endif
}

/** @internal @brief Destruct VLFeat */
void
vl_destructor ()
{
#if defined(DEBUG)
  printf("VLFeat DEBUG: destructor begins.\n") ;
#endif

#if ! defined(VL_DISABLE_THREADS)
  VlState * state = vl_get_state() ;
#if defined(DEBUG)
  printf("VLFeat DEBUG: destroying a thread specific state instance.\n") ;
#endif
#if   defined(VL_THREADS_POSIX)
  {
    /* Delete the thread state of this thread as the
       destructor is not called by pthread_key_delete or after
       the key is deleted. When the library
       is unloaded, this thread should also be the last one
       using the library, so this is fine.
     */
    VlThreadSpecificState * threadState =
       pthread_getspecific(state->threadKey) ;
    if (threadState) {
      vl_thread_specific_state_delete (threadState) ;
      pthread_setspecific(state->threadKey, NULL) ;
    }
  }
  pthread_cond_destroy (&state->mutexCondition) ;
  pthread_mutex_destroy (&state->mutex) ;
  pthread_key_delete (state->threadKey) ;
#elif defined(VL_THREADS_WIN)
 {
    /* Delete the thread state of this thread as the
       destructor is not called by pthread_key_delete or after
       the key is deleted. When the library
       is unloaded, this thread should also be the last one
       using the library, so this is fine.
     */
    VlThreadSpecificState * threadState =
       TlsGetValue(state->tlsIndex) ;
    if (threadState) {
      vl_thread_specific_state_delete (threadState) ;
      TlsSetValue(state->tlsIndex, NULL) ;
    }
  }
  TlsFree (state->tlsIndex) ;
  DeleteCriticalSection (&state->mutex) ;
#endif
#else
#if defined(DEBUG)
  printf("VLFeat DEBUG: destroying the generic thread state instance (threading support disabled).\n") ;
#endif
  vl_thread_specific_state_delete(vl_get_state()->threadState) ;
#endif

#if defined(DEBUG)
  printf("VLFeat DEBUG: destructor ends.\n") ;
#endif
}
