/*
 * Copyright (c) 2004-2012 Samuel Thibault <samuel.thibault@ens-lyon.org>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Samuel Thibault ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* This is a minimal pthread implementation based on windows functions.
 * It is *not* intended to be complete - just complete enough to get
 * simple applications running.
 */

#ifndef _WIN_PTHREAD_PTHREAD_H
#define _WIN_PTHREAD_PTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <windows.h>
#include <stdio.h>
#include <errno.h>

#ifdef __CYGWIN32__
#include <sys/cygwin.h>
#define unixErrno() cygwin_internal(CW_GET_ERRNO_FROM_WINERROR, (GetLastError())
#else
#define unixErrno() win_toErrno(GetLastError())
#endif

#define setSystemErrno() errno = unixErrno()
#define winPthreadAssertWindows(expr) do { if (!(expr)) { setSystemErrno(); return errno; } } while (0)
#define winPthreadAssertPthread(expr) do { int ret = (expr); if (ret) return ret; } while (0)
#define winPthreadAssert(expr) do { if (!(expr)) return EIO; } while (0)

/***********
 * threads *
 ***********/

typedef DWORD pthread_attr_t;
typedef HANDLE pthread_t;

static inline pthread_t pthread_self(void) {
  return GetCurrentThread();
}

static inline int pthread_equal(pthread_t t1, pthread_t t2) {
  return t1 == t2;
}

static inline int pthread_attr_init (pthread_attr_t *attr) {
  *attr = 0;
  return 0;
}

#define PTHREAD_CREATE_DETACHED 1
static inline int pthread_attr_setdetachstate (pthread_attr_t *attr, int yes) {
  (void)attr;
  (void)yes;
  /* not supported, ignore */
  return 0;
}

static inline int pthread_attr_setstacksize (pthread_attr_t *attr, size_t stacksize) {
  (void)attr;
  (void)stacksize;
  /* not supported, ignore */
  return 0;
}

static inline int pthread_attr_destroy (pthread_attr_t *attr) {
  (void)attr;
  return 0;
}

/* "real" cleanup handling not yet implemented */
typedef struct {
  void (*routine) (void *);
  void *arg;
} __pthread_cleanup_handler;

void pthread_cleanup_push (void (*routine) (void *), void *arg);
#define pthread_cleanup_push(routine, arg) do { \
  __pthread_cleanup_handler __cleanup_handler = {routine, arg};

void pthread_cleanup_pop (int execute);
#define pthread_cleanup_pop(execute) \
  if (execute) __cleanup_handler.routine(__cleanup_handler.arg); \
} while (0);

static inline int pthread_create (
  pthread_t *thread, const pthread_attr_t *attr,
  void * (*fun) (void *), void *arg
) {
  if (attr && *attr)
    return EINVAL;
  winPthreadAssertWindows(*thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) fun, arg, 0, NULL));
  return 0;
}

static inline int pthread_setcancelstate (int state, int *oldstate) {
  (void)state;
  (void)oldstate;
  /* not yet implemented :( */
  return 0;
}

static inline int pthread_cancel (pthread_t thread) {
  /* This is quite harsh :( */
  winPthreadAssertWindows(TerminateThread(thread, 0));
  return 0;
}

static inline void pthread_exit (void *res) {
  ExitThread((DWORD) (DWORD_PTR) res);
}

static inline int pthread_join (pthread_t thread, void **res) {
again:
  switch (WaitForSingleObject(thread, INFINITE)) {
    default:
    case WAIT_FAILED:
      return unixErrno();
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      break;
    case WAIT_TIMEOUT:
      goto again;
    }
  if (res) {
    DWORD _res;
    if (GetExitCodeThread(thread, &_res))
      *res = (void *)(DWORD_PTR)_res;
  }
  return 0;
}

/***********
 * mutexes *
 ***********/

#define PTHREAD_MUTEX_INITIALIZER NULL
typedef HANDLE pthread_mutex_t;
#define PTHREAD_MUTEX_RECURSIVE 1
#define PTHREAD_MUTEX_ERRORCHECK 2
typedef int pthread_mutexattr_t;

static inline int pthread_mutexattr_init(pthread_mutexattr_t *attr) {
  *attr = PTHREAD_MUTEX_RECURSIVE;
  return 0;
}

static inline int pthread_mutexattr_settype(pthread_mutexattr_t *attr, int type) {
  if (type != PTHREAD_MUTEX_RECURSIVE && type != PTHREAD_MUTEX_ERRORCHECK)
    return EINVAL;
  *attr = type;
  return 0;
}

static inline int pthread_mutex_init (pthread_mutex_t *mutex, pthread_mutexattr_t *attr) {
  if (attr && *attr!=PTHREAD_MUTEX_RECURSIVE)
    return EINVAL;
  winPthreadAssertWindows(*mutex = CreateMutex(NULL, FALSE, NULL));
  return 0;
}

static inline int pthread_mutex_unlock (pthread_mutex_t *mutex) {
  winPthreadAssertWindows(ReleaseMutex(*mutex));
  return 0;
}

static inline int pthread_mutex_lock (pthread_mutex_t *mutex);
static inline int __pthread_mutex_alloc_concurrently (pthread_mutex_t *mutex) {
    HANDLE mutex_init_mutex;
  /* Get access to one global named mutex to serialize mutex initialization */
  winPthreadAssertWindows((mutex_init_mutex = CreateMutex(NULL, FALSE, "StarPU mutex init")));
  winPthreadAssertPthread(pthread_mutex_lock(&mutex_init_mutex));
  /* Now we are the one that can initialize it */
  if (!*mutex)
    winPthreadAssertPthread(pthread_mutex_init(mutex,NULL));
  winPthreadAssertPthread(pthread_mutex_unlock(&mutex_init_mutex));
    winPthreadAssertWindows(CloseHandle(mutex_init_mutex));
  return 0;
}

static inline int pthread_mutex_lock (pthread_mutex_t *mutex) {
  if (!*mutex)
    __pthread_mutex_alloc_concurrently (mutex);
again:
  switch (WaitForSingleObject(*mutex, INFINITE)) {
    default:
    case WAIT_FAILED:
      return unixErrno();;
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      return 0;
    case WAIT_TIMEOUT:
      goto again;
  }
}

static inline int pthread_mutex_trylock (pthread_mutex_t *mutex) {
  if (!*mutex)
    __pthread_mutex_alloc_concurrently (mutex);
  switch (WaitForSingleObject(*mutex, 0)) {
    default:
    case WAIT_FAILED:
      return unixErrno();
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      return 0;
    case WAIT_TIMEOUT:
      return EBUSY;
  }
}

static inline int pthread_mutex_destroy (pthread_mutex_t *mutex) {
  winPthreadAssertWindows(CloseHandle(*mutex));
  *mutex = INVALID_HANDLE_VALUE;
  return 0;
}

/********************************************
 * rwlock                                   *
 * VERY LAZY, don't even look at it please! *
 * Should be fine unoptimized for now.      *
 * TODO: FIXME, using conds for instance?   *
 ********************************************/

#define PTHREAD_RWLOCK_INITIALIZER NULL
typedef pthread_mutex_t pthread_rwlock_t;
#define pthread_rwlock_init(lock, attr) pthread_mutex_init(lock, NULL)
#define pthread_rwlock_wrlock(lock) pthread_mutex_lock(lock)
#define pthread_rwlock_rdlock(lock) pthread_mutex_lock(lock)
#define pthread_rwlock_unlock(lock) pthread_mutex_unlock(lock)
#define pthread_rwlock_destroy(lock) pthread_mutex_destroy(lock)

/**************
 * conditions *
 **************/

typedef struct {
  HANDLE sem;
  volatile unsigned nbwait;
} pthread_cond_t;
#define PTHREAD_COND_INITIALIZER { NULL, 0}

#ifndef _TIMESPEC_DEFINED
struct timespec {
  time_t  tv_sec;  /* Seconds */
  long    tv_nsec; /* Nanoseconds */
};
#endif

typedef unsigned pthread_condattr_t;

static inline int pthread_cond_init (pthread_cond_t *cond, const pthread_condattr_t *attr) {
  if (attr)
    return EINVAL;
  winPthreadAssertWindows(cond->sem = CreateSemaphore(NULL, 0, MAXLONG, NULL));
  cond->nbwait = 0;
  return 0;
}

static inline int pthread_cond_timedwait (pthread_cond_t *cond, pthread_mutex_t *mutex, const struct timespec *time) {
  if (!cond->sem)
    winPthreadAssertPthread(pthread_cond_init(cond,NULL));
  cond->nbwait++;
  winPthreadAssertPthread(pthread_mutex_unlock(mutex));
again:
  switch (WaitForSingleObject(cond->sem, time->tv_sec*1000+time->tv_nsec/1000)) {
    default:
    case WAIT_FAILED:
    {
      int error = unixErrno();
      winPthreadAssertPthread(pthread_mutex_lock(mutex));
      return error;
    }
    case WAIT_TIMEOUT:
      goto again;
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      break;
  }
  winPthreadAssertPthread(pthread_mutex_lock(mutex));
  cond->nbwait--;
  return 0;
}

static inline int pthread_cond_wait (pthread_cond_t *cond, pthread_mutex_t *mutex) {
  if (!cond->sem)
    winPthreadAssertPthread(pthread_cond_init(cond,NULL));
  cond->nbwait++;
  winPthreadAssertPthread(pthread_mutex_unlock(mutex));
again:
  switch (WaitForSingleObject(cond->sem, INFINITE)) {
    case WAIT_FAILED:
    {
      int error;
      error = unixErrno();
      winPthreadAssertPthread(pthread_mutex_lock(mutex));
      return error;
    }
    case WAIT_TIMEOUT:
      goto again;
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      break;
  }
  winPthreadAssertPthread(pthread_mutex_lock(mutex));
  cond->nbwait--;
  return 0;
}

static inline int pthread_cond_signal (pthread_cond_t *cond) {
  if (!cond->sem)
    winPthreadAssertPthread(pthread_cond_init(cond,NULL));
  if (cond->nbwait)
    ReleaseSemaphore(cond->sem, 1, NULL);
  return 0;
}

static inline int pthread_cond_broadcast (pthread_cond_t *cond) {
  if (!cond->sem)
    winPthreadAssertPthread(pthread_cond_init(cond,NULL));
  ReleaseSemaphore(cond->sem, cond->nbwait, NULL);
  return 0;
}

static inline int pthread_cond_destroy (pthread_cond_t *cond) {
  if (cond->sem) {
  winPthreadAssertWindows(CloseHandle(cond->sem));
    cond->sem = NULL;
  }
  return 0;
}

/*******
 * TLS *
 *******/

typedef DWORD pthread_key_t;
#define PTHREAD_ONCE_INIT {PTHREAD_MUTEX_INITIALIZER, 0}
typedef struct {
  pthread_mutex_t mutex;
  unsigned done;
} pthread_once_t;

static inline int pthread_once (pthread_once_t *once, void (*oncefun)(void)) {
  winPthreadAssertPthread(pthread_mutex_lock(&once->mutex));
  if (!once->done) {
    oncefun();
    once->done = 1;
  }
  winPthreadAssertPthread(pthread_mutex_unlock(&once->mutex));
  return 0;
}

static inline int pthread_key_create (pthread_key_t *key, void (*freefun)(void *)) {
  (void)freefun;
  DWORD res;
  winPthreadAssertWindows((res = TlsAlloc()) != 0xFFFFFFFF);
  *key = res;
  return 0;
}

static inline int pthread_key_delete (pthread_key_t key) {
  winPthreadAssertWindows(TlsFree(key));
  return 0;
}

static inline void *pthread_getspecific (pthread_key_t key) {
  return TlsGetValue(key);
}

static inline int pthread_setspecific (pthread_key_t key, const void *data) {
  winPthreadAssertWindows(TlsSetValue(key, (LPVOID) data));
  return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _WIN_PTHREAD_PTHREAD_H */
