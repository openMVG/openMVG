/* Copyright 2010-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_df.c                         **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a k-way partition  **/
/**                of the given mapping graph by applying  **/
/**                a diffusion method to what is assumed   **/
/**                to be a band graph.                     **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 jan 2010     **/
/**                                 to   : 04 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_DF

#ifdef SCOTCH_PTHREAD
#define KGRAPHMAPDFTHREAD
#endif /* SCOTCH_PTHREAD */

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "mapping.h"
#include "kgraph.h"
#include "kgraph_map_df.h"

/************************/
/*                      */
/* The sorting routine. */
/*                      */
/************************/

/* This routine sorts an array of KgraphMapDfVertex
** values in descending order by their amount of liquid.
** By nature of the sorting algorithm, data are left in
** place in case of equality. Therefore, the original
** part of the vertex, which is put first in the sort
** array during the diffusion process, is always preserved
** when all liquid amounts are equal.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTQUAL                 static
#define INTSORTNAME                 kgraphMapDfSort
#define INTSORTSIZE                 (sizeof (KgraphMapDfSort))
#define INTSORTSWAP(p,q)            do {                                                       \
                                      KgraphMapDfSort t;                                       \
                                      t = *((KgraphMapDfSort *) (p));                          \
                                      *((KgraphMapDfSort *) (p)) = *((KgraphMapDfSort *) (q)); \
                                      *((KgraphMapDfSort *) (q)) = t;                          \
                                    } while (0)
#define INTSORTCMP(p,q)             (((KgraphMapDfSort *) (p))->diffval > ((KgraphMapDfSort *) (q))->diffval)
#include "common_sort.c"
#undef INTSORTQUAL
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/********************************/
/*                              */
/* The sequential loop routine. */
/*                              */
/********************************/

#define KGRAPHMAPDFLOOPNAME         kgraphMapDfSeq
#include "kgraph_map_df_loop.c"
#undef KGRAPHMAPDFLOOPNAME

/******************************/
/*                            */
/* The threaded loop routine. */
/*                            */
/******************************/

#ifdef KGRAPHMAPDFTHREAD

#define KGRAPHMAPDFLOOPTHREAD
#define KGRAPHMAPDFLOOPNAME         kgraphMapDfThr
#include "kgraph_map_df_loop.c"
#undef KGRAPHMAPDFLOOPNAME
#undef KGRAPHMAPDFLOOPTHREAD

#endif /* KGRAPHMAPDFTHREAD */

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a k-way partition
** by diffusion across what is assumed
** to be a k-way band graph.
** It returns:
** - 0   : if the k-partition could be computed.
** - !0  : on error.
*/

int
kgraphMapDf (
Kgraph * restrict const        grafptr,           /*+ Active graph      +*/
const KgraphMapDfParam * const paraptr)           /*+ Method parameters +*/
{
  KgraphMapDfData                  loopdat;       /* Diffusion loop data */
  Gnum                             vertnbr;
  Gnum                             vancnbr;
  Gnum                             domnnbr;
#ifdef KGRAPHMAPDFTHREAD                          /* Threads can be accepted even when SCOTCH_DETERMINISTIC set */
  int                              thrdnbr;
#endif /* KGRAPHMAPDFTHREAD */

  domnnbr = grafptr->m.domnnbr;
  vertnbr = grafptr->s.vertnbr;
  vancnbr = vertnbr - domnnbr;
  if (memAllocGroup ((void **) (void *)
                     &loopdat.vanctab, (size_t) (domnnbr * sizeof (float)),
                     &loopdat.valotab, (size_t) (domnnbr * sizeof (Gnum)),
                     &loopdat.velstax, (size_t) (vertnbr * sizeof (Gnum)),
                     &loopdat.difntax, (size_t) (vertnbr * sizeof (KgraphMapDfVertex)),
                     &loopdat.difotax, (size_t) (vertnbr * sizeof (KgraphMapDfVertex)), NULL) == NULL) {
    errorPrint ("kgraphMapDf: out of memory (1)");
    return     (1);
  }
  loopdat.grafptr  = grafptr;
  loopdat.velstax -= grafptr->s.baseval;
  loopdat.difntax -= grafptr->s.baseval;
  loopdat.difotax -= grafptr->s.baseval;
  loopdat.passnbr  = paraptr->passnbr;

#ifdef KGRAPHMAPDFTHREAD                          /* Threads can be accepted even when SCOTCH_DETERMINISTIC set */
  thrdnbr = SCOTCH_PTHREAD_NUMBER;

  loopdat.abrtval = 0;                            /* No one wants to abort yet */

  if (thrdnbr > 1) {
    KgraphMapDfThread * restrict  thrdtab;
    int                           thrdnum;
    Gnum                          vertbas;
    Anum                          domnbas;

    if ((thrdtab = memAlloc (thrdnbr * sizeof (KgraphMapDfThread))) == NULL) {
      errorPrint ("kgraphMapDf: out of memory (2)");
      memFree    (loopdat.vanctab);
      return     (1);
    }

    for (thrdnum = 0, vertbas = grafptr->s.baseval, domnbas = 0;
         thrdnum < thrdnbr; thrdnum ++) {
      thrdtab[thrdnum].vertbas = vertbas;
      thrdtab[thrdnum].vertnnd = vertbas += DATASIZE (vancnbr, thrdnbr, thrdnum);
      thrdtab[thrdnum].domnbas = domnbas;
      thrdtab[thrdnum].domnnnd = domnbas += DATASIZE (domnnbr, thrdnbr, thrdnum);
    }

    threadLaunch (&loopdat, thrdtab, sizeof (KgraphMapDfThread),
                  (ThreadLaunchStartFunc) kgraphMapDfThr,
                  (ThreadLaunchJoinFunc)  NULL, thrdnbr, THREADCANBARRIER);

    memFree (thrdtab);                            /* Free group leader */
  }
  else
#endif /* KGRAPHMAPDFTHREAD */
  {
    KgraphMapDfThread   thrddat;

    thrddat.thrddat.grouptr = &loopdat;
    thrddat.vertbas = grafptr->s.baseval;
    thrddat.vertnnd = vertnbr + grafptr->s.baseval - domnnbr;
    thrddat.domnbas = 0;
    thrddat.domnnnd = domnnbr;
#ifdef KGRAPHMAPDFTHREAD
    thrddat.thrddat.thrdnum = 0;                  /* Thread is thread 0 of 1 */
#endif /* KGRAPHMAPDFTHREAD */

    kgraphMapDfSeq (&thrddat);
  }

  memFree (loopdat.vanctab);                      /* Free group leader */

  kgraphFron (grafptr);
  kgraphCost (grafptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapDf: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  return (0);
}
