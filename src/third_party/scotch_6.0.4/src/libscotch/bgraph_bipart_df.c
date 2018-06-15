/* Copyright 2004,2007,2008,2011-2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_df.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of   **/
/**                a bipartition graph by using a          **/
/**                diffusion scheme.                       **/
/**                                                        **/
/**   NOTES      : # This algorithm has been designed to   **/
/**                  work on band graphs only, for which   **/
/**                  the two anchor vertices are the two   **/
/**                  last vertices, the before-last as     **/
/**                  anchor of part 0, and the last as     **/
/**                  anchor of part 1.                     **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 09 jan 2007     **/
/**                                 to     10 sep 2007     **/
/**                # Version 5.1  : from : 29 oct 2007     **/
/**                                 to     27 mar 2011     **/
/**                # Version 6.0  : from : 07 nov 2011     **/
/**                                 to     08 aug 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH_BIPART_DF

#ifdef SCOTCH_PTHREAD
#define BGRAPHBIPARTDFTHREAD
#endif /* SCOTCH_PTHREAD */

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_df.h"

/************************************/
/*                                  */
/* The threaded reduction routines. */
/*                                  */
/************************************/

#ifdef BGRAPHBIPARTDFTHREAD

static
void
bgraphBipartDfReduceVanc (
BgraphBipartDfThread * restrict const tlocptr,    /* Pointer to local thread */
void * restrict const                 vlocptr,    /* Pointer to local value  */
void * restrict const                 vremptr)    /* Pointer to remote value */
{
  BgraphBipartDfThread * restrict const tremptr = (BgraphBipartDfThread *) vremptr;

  tlocptr->vanctab[0] += tremptr->vanctab[0];     /* Accumulate external gains */
  tlocptr->vanctab[1] += tremptr->vanctab[1];
}

static
void
bgraphBipartDfReduceVeex (
BgraphBipartDfThread * restrict const tlocptr,    /* Pointer to local thread */
void * restrict const                 vlocptr,    /* Pointer to local value  */
void * restrict const                 vremptr)    /* Pointer to remote value */
{
  BgraphBipartDfThread * restrict const tremptr = (BgraphBipartDfThread *) vremptr;

  tlocptr->veexsum  += tremptr->veexsum;          /* Accumulate external gains */
  tlocptr->veexsum1 += tremptr->veexsum1;
}

#endif /* BGRAPHBIPARTDFTHREAD */

/********************************/
/*                              */
/* The sequential loop routine. */
/*                              */
/********************************/

#define BGRAPHBIPARTDFLOOPNAME      bgraphBipartDfSeq
#include "bgraph_bipart_df_loop.c"
#undef BGRAPHBIPARTDFLOOPNAME

/******************************/
/*                            */
/* The threaded loop routine. */
/*                            */
/******************************/

#ifdef BGRAPHBIPARTDFTHREAD

#define BGRAPHBIPARTDFLOOPTHREAD
#define BGRAPHBIPARTDFLOOPNAME      bgraphBipartDfThr
#include "bgraph_bipart_df_loop.c"
#undef BGRAPHBIPARTDFLOOPNAME
#undef BGRAPHMAPARTDFLOOPTHREAD

#endif /* BGRAPHBIPARTDFTHREAD */

/******************************/
/*                            */
/* The threaded join routine. */
/*                            */
/******************************/

#ifdef BGRAPHBIPARTDFTHREAD

int
bgraphBipartDfJoin (
BgraphBipartDfThread *  tlocptr,
BgraphBipartDfThread *  tremptr)
{
  BgraphBipartDfData * restrict const loopptr = (BgraphBipartDfData *) tlocptr->thrddat.grouptr;
  Bgraph * restrict const             grafptr = loopptr->grafptr;
  Gnum * restrict const               frontab = grafptr->frontab;
  Gnum                                fronnbr;

  fronnbr = tremptr->fronnnd - (tremptr->vertbas - grafptr->s.baseval);
  memMov (frontab + tlocptr->fronnnd,             /* Aggregate frontier array; TODO: we should do differently for large number of threads */
          frontab + (tremptr->vertbas - grafptr->s.baseval), fronnbr * sizeof (Gnum));
  tlocptr->fronnnd      += fronnbr;               /* Accumulate frontier vertices */
  tlocptr->compload1    += tremptr->compload1;    /* Accumulate graph properties  */
  tlocptr->compsize1    += tremptr->compsize1;
  tlocptr->commloadextn += tremptr->commloadextn;
  tlocptr->commloadintn += tremptr->commloadintn;
  tlocptr->commgainextn += tremptr->commgainextn;

  return (0);
}

#endif /* BGRAPHBIPARTDFTHREAD */

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bgraphBipartDf (
Bgraph * restrict const           grafptr,        /*+ Active graph      +*/
const BgraphBipartDfParam * const paraptr)        /*+ Method parameters +*/
{
  BgraphBipartDfData  loopdat;
  Gnum                compload0;
  Gnum                compload1;
  Gnum                compsize1;
  Gnum                commloadintn;
  Gnum                commloadextn;
  Gnum                commgainextn;
#ifdef BGRAPHBIPARTDFTHREAD                       /* Threads can be accepted even when SCOTCH_DETERMINISTIC set */
  int                 thrdnbr;
#endif /* BGRAPHBIPARTDFTHREAD */
  Gnum                fronnbr;

#ifdef SCOTCH_DEBUG_BGRAPH1
  if ((grafptr->s.flagval & BGRAPHHASANCHORS) == 0) { /* Method valid only if graph has anchors */
    errorPrint ("bgraphBipartDf: graph does not have anchors");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH1 */

  if (memAllocGroup ((void **) (void *)
                     &loopdat.difotax, (size_t) (grafptr->s.vertnbr * sizeof (float)),
                     &loopdat.difntax, (size_t) (grafptr->s.vertnbr * sizeof (float)), NULL) == NULL) {
    errorPrint ("bgraphBipartDf: out of memory (1)");
    return     (1);
  }

  loopdat.grafptr  = grafptr;
  loopdat.difotax -= grafptr->s.baseval;
  loopdat.difntax -= grafptr->s.baseval;
  loopdat.passnbr  = paraptr->passnbr;

  compload0 = (paraptr->typeval == BGRAPHBIPARTDFTYPEBAL) /* If balanced parts wanted */
              ? grafptr->compload0avg             /* Target is average                */
              : ( (grafptr->compload0 < grafptr->compload0min) ? grafptr->compload0min : /* Else keep load if not off balance */
                 ((grafptr->compload0 > grafptr->compload0max) ? grafptr->compload0max : grafptr->compload0));
  loopdat.vanctab[0] = (float) - compload0;       /* Values to be injected to anchor vertices at every iteration                */
  loopdat.vanctab[1] = (float) (grafptr->s.velosum - compload0)- BGRAPHBIPARTDFEPSILON; /* Slightly tilt value to add to part 1 */

#ifdef BGRAPHBIPARTDFTHREAD                       /* Threads can be accepted even when SCOTCH_DETERMINISTIC set */
  thrdnbr = SCOTCH_PTHREAD_NUMBER;

  loopdat.abrtval = 0;                            /* No one wants to abort yet */

  if (thrdnbr > 1) {
    BgraphBipartDfThread * restrict thrdtab;
    int                             thrdnum;
    Gnum                            vertbas;

    if ((thrdtab = memAlloc (thrdnbr * sizeof (BgraphBipartDfThread))) == NULL) {
      errorPrint ("bgraphBipartDf: out of memory (2)");
      memFree    (loopdat.difotax + grafptr->s.baseval);
      return     (1);
    }

    for (thrdnum = 0, vertbas = grafptr->s.baseval; /* For all threads except the last one */
         thrdnum < (thrdnbr - 1); thrdnum ++) {
      thrdtab[thrdnum].vertbas = vertbas;
      thrdtab[thrdnum].vertnnd = vertbas += DATASIZE ((grafptr->s.vertnbr - 2), thrdnbr, thrdnum); /* Do not count anchors in distribution */
    }
    thrdtab[thrdnum].vertbas = vertbas;
    thrdtab[thrdnum].vertnnd = grafptr->s.vertnnd; /* Both anchors will always be on the same thread */

    threadLaunch (&loopdat, thrdtab, sizeof (BgraphBipartDfThread),
                  (ThreadLaunchStartFunc) bgraphBipartDfThr,
                  (ThreadLaunchJoinFunc)  bgraphBipartDfJoin, thrdnbr, THREADCANBARRIER | THREADCANREDUCE);

    fronnbr      = thrdtab[0].fronnnd;
    compload1    = thrdtab[0].compload1;
    compsize1    = thrdtab[0].compsize1;
    commloadextn = thrdtab[0].commloadextn;
    commloadintn = thrdtab[0].commloadintn;
    commgainextn = thrdtab[0].commgainextn;

    memFree (thrdtab);                            /* Free group leader */
  }
  else
#endif /* BGRAPHBIPARTDFTHREAD */
  {
    BgraphBipartDfThread  thrddat;

    thrddat.thrddat.grouptr = &loopdat;
    thrddat.vertbas = grafptr->s.baseval;         /* Process all vertices, including both anchors */
    thrddat.vertnnd = grafptr->s.vertnnd;
#ifdef BGRAPHBIPARTDFTHREAD
    loopdat.thrddat.thrdnbr = 1;                  /* Thread is thread 0 of 1 */
    thrddat.thrddat.thrdnum = 0;
#endif /* BGRAPHBIPARTDFTHREAD */

    bgraphBipartDfSeq (&thrddat);

    fronnbr      = thrddat.fronnnd;
    compload1    = thrddat.compload1;
    compsize1    = thrddat.compsize1;
    commloadextn = thrddat.commloadextn;
    commloadintn = thrddat.commloadintn;
    commgainextn = thrddat.commgainextn;
  }

  memFree (loopdat.difotax + grafptr->s.baseval); /* Free group leader */

  grafptr->fronnbr      = fronnbr;
  grafptr->compload0    = grafptr->s.velosum - compload1;
  grafptr->compload0dlt = grafptr->compload0 - grafptr->compload0avg;
  grafptr->compsize0    = grafptr->s.vertnbr - compsize1;
  grafptr->commload     = commloadextn + (commloadintn / 2) * grafptr->domndist;
  grafptr->commgainextn = commgainextn;
  grafptr->bbalval      = (double) ((grafptr->compload0dlt < 0) ? (- grafptr->compload0dlt) : grafptr->compload0dlt) / (double) grafptr->compload0avg;

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (grafptr) != 0) {
    errorPrint ("bgraphBipartDf: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
