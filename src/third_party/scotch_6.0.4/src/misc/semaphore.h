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

#ifndef _WIN_PTHREAD_SEMAPHORE_H
#define _WIN_PTHREAD_SEMAPHORE_H

#include "pthread.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**************
 * semaphores *
 **************/

typedef HANDLE sem_t;

static inline int sem_init(sem_t *sem, int pshared, unsigned int value) {
  (void)pshared;
  winPthreadAssertWindows(*sem = CreateSemaphore(NULL, value, MAXLONG, NULL));
  return 0;
}

static inline int do_sem_wait(sem_t *sem, DWORD timeout) {
  switch (WaitForSingleObject(*sem, timeout)) {
    default:
    case WAIT_FAILED:
      setSystemErrno();
      return -1;
    case WAIT_TIMEOUT:
      errno = EAGAIN;
      return -1;
    case WAIT_ABANDONED:
    case WAIT_OBJECT_0:
      return 0;
  }
}

#define sem_wait(sem) do_sem_wait(sem, INFINITE)
#define sem_trywait(sem) do_sem_wait(sem, 0)

static inline int sem_post(sem_t *sem) {
  winPthreadAssertWindows(ReleaseSemaphore(*sem, 1, NULL));
  return 0;
}

static inline int sem_destroy(sem_t *sem) {
  winPthreadAssertWindows(CloseHandle(*sem));
  return 0;
}

#endif /* _WIN_PTHREAD_SEMAPHORE_H */
