/* Copyright 2007-2009,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_halo.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                Nicolas GICQUEL (P0.1)                  **/
/**                Jerome LACOSTE (P0.1)                   **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel static mapper.       **/
/**                This module contains the halo update    **/
/**                routines.                               **/
/**                                                        **/
/**                # Version P0.0 : from : 01 apr 1997     **/
/**                                 to     20 jun 1997     **/
/**                # Version P0.1 : from : 14 apr 1998     **/
/**                                 to     20 jun 1998     **/
/**                # Version 5.0  : from : 28 feb 2006     **/
/**                                 to     05 feb 2008     **/
/**                # Version 5.1  : from : 28 aug 2008     **/
/**                                 to     22 feb 2011     **/
/**                # Version 6.0  : from : 29 oct 2014     **/
/**                                 to     29 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_HALO

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_halo.h"

/* These routines fill the send arrays used by
** all of the halo routines.
** They return:
** - void  : in all cases.
*/

#define DGRAPHHALOFILLNAME          dgraphHaloFillGeneric
#define DGRAPHHALOFILLSIZE          attrglbsiz
#define DGRAPHHALOFILLCOPY(d,s,n)   memCpy ((d), (s), (n))
#include "dgraph_halo_fill.c"
#undef DGRAPHHALOFILLNAME
#undef DGRAPHHALOFILLSIZE
#undef DGRAPHHALOFILLCOPY

#define DGRAPHHALOFILLNAME          dgraphHaloFillGnum
#define DGRAPHHALOFILLSIZE          sizeof (Gnum)
#define DGRAPHHALOFILLCOPY(d,s,n)   *((Gnum *) (d)) = *((Gnum *) (s))
#include "dgraph_halo_fill.c"
#undef DGRAPHHALOFILLNAME
#undef DGRAPHHALOFILLSIZE
#undef DGRAPHHALOFILLCOPY

#define DGRAPHHALOFILLNAME          dgraphHaloFillGraphPart
#define DGRAPHHALOFILLSIZE          sizeof (GraphPart)
#define DGRAPHHALOFILLCOPY(d,s,n)   *((GraphPart *) (d)) = *((GraphPart *) (s))
#include "dgraph_halo_fill.c"
#undef DGRAPHHALOFILLNAME
#undef DGRAPHHALOFILLSIZE
#undef DGRAPHHALOFILLCOPY

#define DGRAPHHALOFILLNAME          dgraphHaloFillInt /* In case Gnum is not int */
#define DGRAPHHALOFILLSIZE          sizeof (int)
#define DGRAPHHALOFILLCOPY(d,s,n)   *((int *) (d)) = *((int *) (s))
#include "dgraph_halo_fill.c"
#undef DGRAPHHALOFILLNAME
#undef DGRAPHHALOFILLSIZE
#undef DGRAPHHALOFILLCOPY

static
void
dgraphHaloFill (
const Dgraph * restrict const grafptr,
const void * restrict const   attrgsttab,         /* Attribute array to diffuse     */
int                           attrglbsiz,         /* Type extent of attribute       */
byte * restrict const         attrsndtab,         /* Array for packing data to send */
int * const                   senddsptab,         /* Temporary displacement array   */
const int * restrict const    sendcnttab)         /* Count array                    */
{
  int                 procnum;
  byte **             attrdsptab;

  attrdsptab = (byte **) senddsptab;              /* TRICK: use senddsptab (int *) as attrdsptab (byte **) */
  attrdsptab[0] = attrsndtab;                     /* Pre-set send arrays for send buffer filling routines */
  for (procnum = 1; procnum < grafptr->procglbnbr; procnum ++)
    attrdsptab[procnum] = attrdsptab[procnum - 1] + sendcnttab[procnum - 1] * attrglbsiz;

  if (attrglbsiz == sizeof (Gnum))
    dgraphHaloFillGnum (grafptr, attrgsttab, attrglbsiz, attrdsptab);
  else if (attrglbsiz == sizeof (GraphPart))
    dgraphHaloFillGraphPart (grafptr, attrgsttab, attrglbsiz, attrdsptab);
  else if (attrglbsiz == sizeof (int))            /* In case Gnum is not int (suitable for float's too) */
    dgraphHaloFillInt (grafptr, attrgsttab, attrglbsiz, attrdsptab);
  else                                            /* Generic but slower fallback routine */
    dgraphHaloFillGeneric (grafptr, attrgsttab, attrglbsiz, attrdsptab);

  senddsptab[0] = 0;                              /* Pre-set send arrays for data sending routines */
  for (procnum = 1; procnum < grafptr->procglbnbr; procnum ++)
    senddsptab[procnum] = senddsptab[procnum - 1] + sendcnttab[procnum - 1];
}

/* This function checks that the data of proc{snd,rcv}tab
** are consistent.
*/

#ifdef SCOTCH_DEBUG_DGRAPH2
int
dgraphHaloCheck (
const Dgraph * restrict const grafptr)
{
  int *               proctab;                    /* Array to collect data */
  int                 procnum;
  int                 o;

  if ((proctab = memAlloc (grafptr->procglbnbr * sizeof (int))) == NULL) {
    errorPrint ("dgraphHaloCheck: out of memory");
    return     (1);
  }

  if (MPI_Alltoall (grafptr->procsndtab, 1, MPI_INT, proctab, 1, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphHaloCheck: communication error");
    memFree    (proctab);                         /* Free group leader */
    return     (1);
  }

  o = 0;
  for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
    if (proctab[procnum] != grafptr->procrcvtab[procnum]) {
      errorPrint ("dgraphHaloCheck: data error");
      o = 1;
      break;
    }
  }

  memFree (proctab);                              /* Free group leader */

  return  (o);
}
#endif /* SCOTCH_DEBUG_DGRAPH2 */

/* These functions perform a synchronous collective
** halo diffusion operation on the ghost array given
** on input.
** It returns:
** - 0   : if the halo has been successfully propagated.
** - !0  : on error.
*/

static
int
dgraphHaloSync2 (
Dgraph * restrict const       grafptr,
void * restrict const         attrgsttab,         /* Attribute array to share                          */
const MPI_Datatype            attrglbtype,        /* Attribute datatype                                */
byte ** const                 attrsndptr,         /* Pointer to array for packing data to send         */
int ** const                  senddspptr,         /* Pointers to communication displacement arrays     */
int ** const                  recvdspptr,
MPI_Request ** const          requptr)            /* Pointer to local request array for point-to-point */
{
#if ! ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
  MPI_Aint              attrglbtmp;               /* Lower bound of attribute datatype (not used)            */
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  MPI_Aint              attrglbsiz;               /* Extent of attribute datatype                            */
  int                   procngbsiz;               /* Size of request array for point-to-point communications */
  int                   procngbnum;
  int * restrict        recvdsptab;
  const int * restrict  procrcvtab;

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint ("dgraphHaloSync2: cannot compute ghost edge array");
    return     (1);
  }

  procngbsiz = ((grafptr->flagval & DGRAPHCOMMPTOP) != 0) ? grafptr->procngbnbr : 0;

#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
  MPI_Type_extent (attrglbtype, &attrglbsiz);     /* Get type extent */
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  MPI_Type_get_extent (attrglbtype, &attrglbtmp, &attrglbsiz); /* Get type extent */
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  if (memAllocGroup ((void **) (void *)
                     attrsndptr, (size_t) (grafptr->procsndnbr * attrglbsiz),
                     senddspptr, (size_t) (grafptr->procglbnbr * MAX (sizeof (int), sizeof (byte *))), /* TRICK: use senddsptab (int *) as attrdsptab (byte **) */
                     recvdspptr, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     requptr,    (size_t) (procngbsiz * 2      * sizeof (MPI_Request)), NULL) == NULL) {
    errorPrint ("dgraphHaloSync2: out of memory");
    return     (1);
  }

  dgraphHaloFill (grafptr, attrgsttab, attrglbsiz, *attrsndptr, *senddspptr, grafptr->procsndtab); /* Fill data arrays */

  recvdsptab = *recvdspptr;
  procrcvtab = grafptr->procrcvtab;
  recvdsptab[0] = grafptr->vertlocnbr;
  for (procngbnum = 1; procngbnum < grafptr->procglbnbr; procngbnum ++)
    recvdsptab[procngbnum] = recvdsptab[procngbnum - 1] + procrcvtab[procngbnum - 1];

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphHaloCheck (grafptr) != 0) {
    errorPrint ("dgraphHaloSync2: invalid communication data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}

int
dgraphHaloSync (
Dgraph * restrict const       grafptr,
void * restrict const         attrgsttab,         /* Attribute array to share */
const MPI_Datatype            attrglbtype)        /* Attribute datatype       */
{
  byte *              attrsndtab;                 /* Array for packing data to send */
  int *               senddsptab;
  int *               recvdsptab;
  MPI_Request *       requtab;
  int                 o;

  if (dgraphHaloSync2 (grafptr, attrgsttab, attrglbtype, &attrsndtab, &senddsptab, &recvdsptab, &requtab) != 0) /* Prepare communication arrays */
    return (1);

  o = 0;                                          /* Assume success                               */
  if ((grafptr->flagval & DGRAPHCOMMPTOP) != 0) { /* If point-to-point exchange                   */
#if ! ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
    MPI_Aint              attrglbtmp;             /* Lower bound of attribute datatype (not used) */
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
    MPI_Aint              attrglbsiz;             /* Extent of attribute datatype                 */
    const int * restrict  procrcvtab;
    const int * restrict  procsndtab;
    const int * restrict  procngbtab;
    int                   procngbnbr;
    int                   procngbnum;
    MPI_Comm              proccomm;
    int                   requnbr;

    proccomm   = grafptr->proccomm;
    procngbtab = grafptr->procngbtab;
    procngbnbr = grafptr->procngbnbr;
    procrcvtab = grafptr->procrcvtab;
#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
    MPI_Type_extent (attrglbtype, &attrglbsiz);   /* Get type extent */
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
    MPI_Type_get_extent (attrglbtype, &attrglbtmp, &attrglbsiz); /* Get type extent */
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
    for (procngbnum = procngbnbr - 1, requnbr = 0; procngbnum >= 0; procngbnum --, requnbr ++) { /* Post receives first */
      int                 procglbnum;

      procglbnum = procngbtab[procngbnum];
      if (MPI_Irecv ((byte *) attrgsttab + recvdsptab[procglbnum] * attrglbsiz, procrcvtab[procglbnum],
                     attrglbtype, procglbnum, TAGHALO, proccomm, requtab + requnbr) != MPI_SUCCESS) {
        errorPrint ("dgraphHaloSync: communication error (1)");
        o = 1;
        break;
      }
    }

    procsndtab = grafptr->procsndtab;
    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++, requnbr ++) { /* Post sends afterwards */
      int                 procglbnum;

      procglbnum = procngbtab[procngbnum];
      if (MPI_Isend (attrsndtab + senddsptab[procglbnum] * attrglbsiz, procsndtab[procglbnum],
                     attrglbtype, procglbnum, TAGHALO, proccomm, requtab + requnbr) != MPI_SUCCESS) {
        errorPrint ("dgraphHaloSync: communication error (2)");
        o = 1;
        break;
      }
    }
    if (MPI_Waitall (requnbr, requtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
      errorPrint ("dgraphHaloSync: communication error (3)");
      o = 1;
    }
  }
  else {                                          /* Collective communication */
    if (MPI_Alltoallv (attrsndtab, grafptr->procsndtab, senddsptab, attrglbtype, /* Perform diffusion */
                       attrgsttab, grafptr->procrcvtab, recvdsptab, attrglbtype,
                       grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphHaloSync: communication error (4)");
      o = 1;
    }
  }

  memFree (attrsndtab);                           /* Free group leader */

  return (o);
}

/* This function performs an asynchronous collective
** halo diffusion operation on the ghost array given
** on input. It fills the given request structure with
** the relevant data.
** It returns:
** - 0   : if the halo has been successfully propagated.
** - !0  : on error.
*/

#ifdef SCOTCH_PTHREAD
static
void *
dgraphHaloAsync2 (
DgraphHaloRequest * restrict  requptr)
{
  return ((void *) (intptr_t) dgraphHaloSync (requptr->grafptr, requptr->attrgsttab, requptr->attrglbtype));
}
#endif /* SCOTCH_PTHREAD */

void
dgraphHaloAsync (
Dgraph * restrict const       grafptr,
void * restrict const         attrgsttab,         /* Attribute array to share */
const MPI_Datatype            attrglbtype,        /* Attribute datatype       */
DgraphHaloRequest * restrict  requptr)
{
#ifndef SCOTCH_PTHREAD
#ifdef SCOTCH_MPI_ASYNC_COLL
  int *               senddsptab;
  int *               recvdsptab;
#endif /* SCOTCH_MPI_ASYNC_COLL */
#endif /* SCOTCH_PTHREAD */

#ifdef SCOTCH_PTHREAD
  requptr->flagval     = -1;                      /* Assume thread will be successfully launched */
  requptr->grafptr     = grafptr;
  requptr->attrgsttab  = attrgsttab;
  requptr->attrglbtype = attrglbtype;

  if (pthread_create (&requptr->thrdval, NULL, (void * (*) (void *)) dgraphHaloAsync2, (void *) requptr) != 0) /* If could not create thread */
    requptr->flagval = (int) (intptr_t) dgraphHaloAsync2 (requptr); /* Call function synchronously */
#else /* SCOTCH_PTHREAD */
#ifdef SCOTCH_MPI_ASYNC_COLL
  requptr->flagval    = 1;                        /* Assume error */
  requptr->attrsndtab = NULL;                     /* No memory    */

  if (dgraphHaloSync2 (grafptr, attrgsttab, attrglbtype, &requptr->attrsndtab, &senddsptab, &recvdsptab) != 0) /* Prepare communication arrays */
    return;

  if (MPE_Ialltoallv (requptr->attrsndtab, grafptr->procsndtab, senddsptab, attrglbtype, /* Perform asynchronous collective communication */
                      attrgsttab, grafptr->procrcvtab, recvdsptab, attrglbtype,
                      grafptr->proccomm, &requptr->requval) != MPI_SUCCESS) {
    errorPrint ("dgraphHaloAsync: communication error"); /* Group leader will be freed on wait routine */
    return;
  }
  requptr->flagval = -1;                          /* Communication successfully launched */
#else /* SCOTCH_MPI_ASYNC_COLL */
  requptr->flagval = dgraphHaloSync (grafptr, attrgsttab, attrglbtype); /* Last resort is synchronous communication */
#endif /* SCOTCH_MPI_ASYNC_COLL */
#endif /* SCOTCH_PTHREAD */
}

/* This function performs an asynchronous collective
** halo diffusion operation on the ghost array given
** on input. It fills the given request structure with
** the relevant data.
** It returns:
** - 0   : if the halo has been successfully propagated.
** - !0  : on error.
*/

int
dgraphHaloWait (
DgraphHaloRequest * restrict  requptr)
{
#ifdef SCOTCH_PTHREAD
  void *                    o;

  if (requptr->flagval == -1) {                   /* If thread launched              */
    pthread_join (requptr->thrdval, &o);          /* Wait for its completion         */
    requptr->flagval = (int) (intptr_t) o;        /* Get thread return value         */
  }                                               /* Else return value already known */
#else /* SCOTCH_PTHREAD */
#ifdef SCOTCH_MPI_ASYNC_COLL
  MPI_Status                statval;

  if (requptr->flagval == -1) {                   /* If communication launched                                    */
    MPI_Wait (&requptr->requval, &statval);       /* Wait for completion of asynchronous collective communication */
    requptr->flagval = (statval.MPI_ERROR == MPI_SUCCESS) ? 0 : 1; /* Set return value accordingly                */
  }
  if (requptr->attrsndtab != NULL)                /* Free group leader if it was successfully allocated before */
    memFree (requptr->attrsndtab);
#endif /* SCOTCH_MPI_ASYNC_COLL */
#endif /* SCOTCH_PTHREAD */

  return (requptr->flagval);                      /* Return asynchronous or synchronous error code */
}
