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
/**   NAME       : bgraph_bipart_df_loop.c                 **/
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
/**                                 to   : 08 aug 2014     **/
/**                                                        **/
/************************************************************/

/****************************/
/*                          */
/* The diffusion subroutine */
/* pattern.                 */
/*                          */
/****************************/

/* This routine computes the diffusion of two
** liquids on the given part of the bipartition
** graph.
** It returns:
** - void  : in all cases
*/

static
int
BGRAPHBIPARTDFLOOPNAME (
BgraphBipartDfThread * restrict thrdptr)          /* Thread-dependent data */
{
  float * restrict      ielstax;                  /* Inverse of edge load sum array   */
  float * restrict      difotax;                  /* Old diffusion value array        */
  float * restrict      difntax;                  /* New diffusion value array        */
  Gnum                  vertnum;
  Gnum                  fronnum;
  Gnum                  compload0;
  Gnum                  compload1;
  Gnum                  compsize1;
  Gnum                  commloadintn;
  Gnum                  commloadextn;
  Gnum                  commgainextn;
  Gnum                  veexnbr;
  Gnum                  veexval;
  Gnum                  veexval1;                 /* Negative external gain to part 1 */
  Gnum                  veexsum;
  Gnum                  veexsum1;                 /* Sum of negative external gains   */
  Gnum                  veloval;
  float                 velfval;
  INT                   passnum;
  Anum                  distval;

  BgraphBipartDfData * restrict const loopptr = (BgraphBipartDfData *) thrdptr->thrddat.grouptr;
  Bgraph * restrict const             grafptr = loopptr->grafptr;
  const Gnum                          vertbas = thrdptr->vertbas;
  const Gnum                          vertnnd = thrdptr->vertnnd;
  const Gnum                          vancnnd = MIN (vertnnd, grafptr->s.vertnnd - 2); /* End of regular vertices */
  const Gnum * restrict const         verttax = grafptr->s.verttax;
  const Gnum * restrict const         vendtax = grafptr->s.vendtax;
  const Gnum * restrict const         velotax = grafptr->s.velotax;
  const Gnum * restrict const         edgetax = grafptr->s.edgetax;
  const Gnum * restrict const         edlotax = grafptr->s.edlotax;
  const Gnum * restrict const         veextax = grafptr->veextax;
  GraphPart * restrict const          parttax = grafptr->parttax;
#ifdef BGRAPHBIPARTDFLOOPTHREAD
  const int                           thrdlst = loopptr->thrddat.thrdnbr - 1;
#endif /* BGRAPHBIPARTDFLOOPTHREAD */


  difotax = loopptr->difotax;
  difntax = loopptr->difntax;

  if ((ielstax = memAlloc ((vertnnd - vertbas) * sizeof (float))) == NULL) { /* Allocate here for memory affinity as it is a private array */
    errorPrint (STRINGIFY (BGRAPHBIPARTDFLOOPNAME) ": out of memory");
#ifdef BGRAPHBIPARTDFLOOPTHREAD
    loopptr->abrtval = 1;
#else /* BGRAPHBIPARTDFLOOPTHREAD */
    return (1);
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
  }
  else {
    ielstax -= vertbas;                           /* Base access to local part of edge load sum array */

    distval  = grafptr->domndist;
    veexval  =                                    /* Assume no external gains */
    veexval1 = 0;
    veexsum  =
    veexsum1 = 0;
    for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) { /* Process all vertices, including anchors */
      Gnum                edlosum;

      if (edlotax == NULL)                        /* If graph doesn't have edge weights */
        edlosum = vendtax[vertnum] - verttax[vertnum];
      else {
        Gnum                edgenum;
        Gnum                edgennd;

        for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum], edlosum = 0;
             edgenum < edgennd; edgenum ++)
          edlosum += edlotax[edgenum];
      }
      edlosum *= distval;

      if (veextax != NULL) {
        veexval  = veextax[vertnum];
        veexval1 = veexval & BGRAPHBIPARTDFGNUMSGNMSK (veexval); /* Get negative external gain only, by superscalar update */
#ifdef SCOTCH_DEBUG_BGRAPH2
        if (((veexval >= 0) && (veexval1 != 0)) ||
            ((veexval <  0) && (veexval1 != veexval))) {
          errorPrint (STRINGIFY (BGRAPHBIPARTDFLOOPNAME) ": internal error");
#ifdef BGRAPHBIPARTDFLOOPTHREAD
          loopptr->abrtval = 1;
#else /* BGRAPHBIPARTDFLOOPTHREAD */
          return (1);
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
        }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

        veexsum  += veexval;                      /* Sum all external gains, positive and negative */
        veexsum1 += veexval1;                     /* Sum all negative gains                        */
      }

      difotax[vertnum] = 0.0;
      ielstax[vertnum] = 1.0F / (float) (edlosum + veexval - 2 * veexval1); /* Add absolute value of veexval */
    }
    if (veextax != NULL) {
#ifdef BGRAPHBIPARTDFLOOPTHREAD
      thrdptr->veexsum  = veexsum;
      thrdptr->veexsum1 = veexsum1;
      threadReduce (thrdptr, thrdptr, (ThreadReduceFunc) bgraphBipartDfReduceVeex, thrdlst);
      veexsum  = thrdptr->veexsum;                /* Will be useful for thread (thrdlst) only */
      veexsum1 = thrdptr->veexsum1;

      if (thrdptr->thrddat.thrdnum == thrdlst)    /* Last thread will handle anchors as root of reduction */
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
      {
        ielstax[vertnnd - 2] = 1.0F / (1.0F / ielstax[vertnnd - 2] + (float) (veexsum - veexsum1));
        ielstax[vertnnd - 1] = 1.0F / (1.0F / ielstax[vertnnd - 1] - (float) veexsum1);
      }
    }
#ifdef BGRAPHBIPARTDFLOOPTHREAD
    if (thrdptr->thrddat.thrdnum == thrdlst)      /* Last thread will handle anchors as root of reduction */
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
    {
      difotax[vertnnd - 2] = loopptr->vanctab[0] * ielstax[vertnnd - 2]; /* Load anchor vertices for first pass */
      difotax[vertnnd - 1] = loopptr->vanctab[1] * ielstax[vertnnd - 1];
    }
  }

#ifdef BGRAPHBIPARTDFLOOPTHREAD
  threadBarrier (thrdptr);                        /* Wait until all array values have been computed */

  if (loopptr->abrtval == 1) {                    /* If process alone or some decided to quit */
    if (ielstax != NULL)                          /* Free local array if necessary            */
      memFree (ielstax + vertbas);
    return  (1);
  }
#endif /* BGRAPHBIPARTDFLOOPTHREAD */

  velfval = 1.0F;                                 /* Assume no vertex loads     */
  for (passnum = loopptr->passnbr; passnum > 0; passnum --) { /* For all passes */
    Gnum                vertnum;
    Gnum                vancnnt;
    float               vancval;                  /* Value to load vertex with if anchor   */
    float *             difttax;                  /* Temporary swap value                  */
    float               vancold0 = difotax[grafptr->s.vertnnd - 2]; /* Get for all threads */
    float               vancold1 = difotax[grafptr->s.vertnnd - 1];
    float               vancval0;                 /* External gain contributions from regular vertices to anchors */
    float               vancval1;

    vancval0 =
    vancval1 = 0.0F;
    vancval  = 0.0F;                              /* At first vertices are not anchors           */
    vertnum  = vertbas;                           /* Start processing regular vertices, then see */
    vancnnt  = vancnnd;                           /* Loop until end of (regular) vertex block    */
    while (1) {
      for ( ; vertnum < vancnnt; vertnum ++) {
        Gnum                edgenum;
        Gnum                edgennd;
        float               diffval;

        edgenum = verttax[vertnum];
        edgennd = vendtax[vertnum];
        diffval = 0.0F;
        if (edlotax != NULL)
          for ( ; edgenum < edgennd; edgenum ++)
            diffval += difotax[edgetax[edgenum]] * (float) edlotax[edgenum];
        else
          for ( ; edgenum < edgennd; edgenum ++)
            diffval += difotax[edgetax[edgenum]];

        diffval *= (float) distval;

        if (veextax != NULL) {
          Gnum                veexval;

          veexval = veextax[vertnum];
          if (veexval != 0) {
            float               veextmp;
            float               vanctmp;

            veextmp = (float) veexval;
            vanctmp = veextmp * difotax[vertnum];

            if (veexval > 0) {                    /* If external gain links to part 0  */
              diffval  += veextmp * vancold0;     /* Spread contribution from anchor 0 */
              vancval0 += vanctmp;
            }
            else {                                /* If external gain links to part 1 */
              diffval  -= veextmp * vancold1;     /* Take opposite of negative value  */
              vancval1 -= vanctmp;
            }
          }
        }

        diffval += vancval;                       /* Add anchor contribution if anchor vertex */

        if (velotax != NULL)
          velfval = (float) velotax[vertnum];
        if (diffval >= 0.0F) {
          diffval -= velfval;
          if (diffval <= 0.0F)
            diffval = +BGRAPHBIPARTDFEPSILON;
        }
        else {
          diffval += velfval;
          if (diffval >= 0.0F)
            diffval = -BGRAPHBIPARTDFEPSILON;
        }
        if (isnan (diffval)) {                    /* If overflow occured (because of avalanche process) */
#ifdef SCOTCH_DEBUG_BGRAPH2
          errorPrintW (STRINGIFY (BGRAPHBIPARTDFLOOPNAME) ": overflow");
#endif /* SCOTCH_DEBUG_BGRAPH2 */
#ifdef BGRAPHBIPARTDFLOOPTHREAD
          loopptr->abrtval = 1;                   /* Threads need to halt                      */
          vertnum = vancnnt;                      /* Skip regular computations but synchronize */
#else /* BGRAPHBIPARTDFLOOPTHREAD */
          goto abort;                             /* Exit this loop without swapping arrays */
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
        }

        difntax[vertnum] = diffval * ielstax[vertnum];
      }
      if (vertnum == vancnnd) {                   /* If first time we reach the end of regular vertices */
        thrdptr->vanctab[0] = vancval0;
        thrdptr->vanctab[1] = vancval1;
#ifdef BGRAPHBIPARTDFLOOPTHREAD
        if (veextax != NULL)
          threadReduce (thrdptr, thrdptr, (ThreadReduceFunc) bgraphBipartDfReduceVanc, thrdlst);
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
      }

      if (vertnum >= vertnnd)                     /* If all vertices processed in range array, exit intermediate infinite loop */
        break;

      vancnnt ++;                                 /* Prepare to go only for one more run, to be done twice                    */
      vancval = loopptr->vanctab[vertnum - vancnnd] + thrdptr->vanctab[vertnum - vancnnd]; /* Load variable with anchor value */
    }

    difttax = (float *) difntax;                  /* Swap old and new diffusion arrays          */
    difntax = (float *) difotax;                  /* Casts to prevent IBM compiler from yelling */
    difotax = (float *) difttax;
#ifdef BGRAPHBIPARTDFLOOPTHREAD
    threadBarrier (thrdptr);

    if (loopptr->abrtval == 1) {                  /* If all threads need to abort       */
      difotax = (float *) difntax;                /* Roll-back to keep last valid array */
      break;
    }
#endif /* BGRAPHBIPARTDFLOOPTHREAD */
  }
abort : ;

  for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) /* Update part according to diffusion state */
    parttax[vertnum] = (difotax[vertnum] <= 0.0F) ? 0 : 1;

#ifdef BGRAPHBIPARTDFLOOPTHREAD
  threadBarrier (thrdptr);
#endif /* BGRAPHBIPARTDFLOOPTHREAD */

  veloval = 1;
  veexval = 0;
  for (vertnum = vertbas, fronnum = vertbas - grafptr->s.baseval,
       commloadextn = commgainextn = commloadintn = compload1 = compsize1 = 0;
       vertnum < vertnnd; vertnum ++) {
    Gnum                edgenum;
    Gnum                partval;
    Gnum                commload;                 /* Vertex internal communication load */

    partval = (Gnum) parttax[vertnum];
    if (velotax != NULL)
      veloval = velotax[vertnum];
    if (veextax != NULL)
      veexval = veextax[vertnum];
    compsize1    += partval;
    compload1    += partval * veloval;
    commloadextn += partval * veexval;
    commgainextn += (1 - 2 * partval) * veexval;
    commload      = 0;
    if (edlotax != NULL) {
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                partend;

        partend   = (Gnum) parttax[edgetax[edgenum]];
        commload += (partval ^ partend) * edlotax[edgenum];
      }
    }
    else {
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++)
        commload += partval ^ (Gnum) parttax[edgetax[edgenum]];
    }
    commloadintn += commload;                     /* Internal loads will be added twice */
    if (commload != 0)                            /* If end vertex is in the other part */
      grafptr->frontab[fronnum ++] = vertnum;     /* Then it belongs to the frontier    */
  }
  thrdptr->fronnnd      = fronnum;                /* Save state */
  thrdptr->compload1    = compload1;
  thrdptr->compsize1    = compsize1;
  thrdptr->commloadextn = commloadextn;
  thrdptr->commloadintn = commloadintn;
  thrdptr->commgainextn = commgainextn;

  memFree (ielstax + vertbas);                    /* Free local part of (local part of) edge load sum array */

  return (0);
}
