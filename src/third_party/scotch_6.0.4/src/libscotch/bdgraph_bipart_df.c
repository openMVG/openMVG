/* Copyright 2007,2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_df.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of   **/
/**                the given distributed separator graph   **/
/**                by applying a diffusion method to what  **/
/**                is assumed to be a distributed band     **/
/**                graph.                                  **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 nov 2007     **/
/**                                 to   : 19 jul 2011     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to   : 31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_BIPART_DF

#include "module.h"
#include "common.h"
#include "arch.h"
#include "dgraph.h"
#include "bdgraph.h"
#include "bdgraph_bipart_df.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed bipartition
** by diffusion across what is assumed to be a
** distributed band graph.
** It returns:
** - 0   : if the bipartition could be computed.
** - !0  : on error.
*/

int
bdgraphBipartDf (
Bdgraph * const                     grafptr,      /*+ Distributed graph +*/
const BdgraphBipartDfParam * const  paraptr)      /*+ Method parameters +*/
{
  float * restrict        ielsloctax;             /* Inverse of degree array   */
  float * restrict        veexloctax;             /* Veexval over domndist     */
  float * restrict        difogsttax;             /* Old diffusion value array */
  float * restrict        difngsttax;             /* New diffusion value array */
  const Gnum * restrict   edgegsttax;
  Gnum                    fronlocnum;
  Gnum                    veexlocnbr;
  float                   vanclocval[2];
  float                   valolocval[2];          /* Fraction of load to remove from anchor vertices at each step */
  Gnum                    vanclocnnd;
  Gnum                    vertlocnum;
  Gnum                    complocload1;
  Gnum                    complocsize1;
  Gnum                    commlocloadintn;
  Gnum                    commlocloadextn;
  Gnum                    commlocgainextn;
  Gnum                    reduloctab[6];
  Gnum                    reduglbtab[6];
  Gnum                    passnum;
  float                   cdifval;
  float                   cremval;
  int                     ovflval;                /* Overflow flag value */

  const Gnum * restrict const veloloctax = grafptr->s.veloloctax;
  const Gnum * restrict const edloloctax = grafptr->s.edloloctax;

  if (dgraphGhst (&grafptr->s) != 0) {            /* Compute ghost edge array if not already present */
    errorPrint ("bdgraphBipartDf: cannot compute ghost edge array");
    return     (1);
  }

  reduloctab[0] = grafptr->s.vendloctax[grafptr->s.vertlocnnd - 2] - grafptr->s.vertloctax[grafptr->s.vertlocnnd - 2] - (grafptr->s.procglbnbr - 1); /* Local degree of both anchor vertices, minus edges to other anchors */
  reduloctab[1] = grafptr->s.vendloctax[grafptr->s.vertlocnnd - 1] - grafptr->s.vertloctax[grafptr->s.vertlocnnd - 1] - (grafptr->s.procglbnbr - 1); /* Anchor edges have load 1 even for weighted graphs                  */
  if (grafptr->s.veloloctax == NULL)
    reduloctab[2] =                               /* Weights of anchors */
    reduloctab[3] = 1;
  else {
    reduloctab[2] = grafptr->s.veloloctax[grafptr->s.vertlocnnd - 2];
    reduloctab[3] = grafptr->s.veloloctax[grafptr->s.vertlocnnd - 1];
  }

  veexlocnbr = (grafptr->veexloctax != NULL) ? grafptr->s.vertlocnbr : 0;
  if (memAllocGroup ((void **) (void *)
                     &ielsloctax, (size_t) (grafptr->s.vertlocnbr * sizeof (float)),
                     &veexloctax, (size_t) (veexlocnbr            * sizeof (float)),
                     &difogsttax, (size_t) (grafptr->s.vertgstnbr * sizeof (float)),
                     &difngsttax, (size_t) (grafptr->s.vertgstnbr * sizeof (float)), NULL) == NULL) {
    errorPrint ("bdgraphBipartDf: out of memory");
    reduloctab[0] = -1;
  }
  ielsloctax -= grafptr->s.baseval;
  difogsttax -= grafptr->s.baseval;
  difngsttax -= grafptr->s.baseval;
  veexloctax  = (grafptr->veexloctax != NULL) ? (veexloctax - grafptr->s.baseval) : NULL;

  if (MPI_Allreduce (reduloctab, reduglbtab, 4, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartDf: communication error (1)");
    return     (1);
  }

  if (reduglbtab[0] < 0) {
    if (ielsloctax != NULL)
      memFree (ielsloctax + grafptr->s.baseval);  /* Free group leader */
    return (1);
  }
  if ((reduglbtab[0] == 0) ||                     /* If graph is too small to have any usable anchors, leave partition as is */
      (reduglbtab[1] == 0)) {
    memFree (ielsloctax + grafptr->s.baseval);

    if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (grafptr->partgsttax + grafptr->s.baseval), GRAPHPART_MPI) != 0) {
      errorPrint ("bdgraphBipartDf: cannot propagate part data (1)");
      return     (1);
    }

    return  (0);
  }

  vanclocval[0] = (float) ((paraptr->typeval == BDGRAPHBIPARTDFTYPEBAL) /* If balanced parts wanted */
                           ? grafptr->compglbload0avg /* Target is average                          */
                           : ( (grafptr->compglbload0 < grafptr->compglbload0min) ? grafptr->compglbload0min : /* Else keep load if not off balance */
                              ((grafptr->compglbload0 > grafptr->compglbload0max) ? grafptr->compglbload0max : grafptr->compglbload0)));
  vanclocval[1] = (float) grafptr->s.veloglbsum - vanclocval[0];
  vanclocval[0] = - vanclocval[0];                /* Part 0 holds negative values                         */
  valolocval[0] = (float) reduglbtab[2];          /* Compute values to remove from anchor vertices        */
  valolocval[1] = (float) reduglbtab[3] - BDGRAPHBIPARTDFEPSILON; /* Slightly tilt value to add to part 1 */

  vanclocnnd = grafptr->s.vertlocnnd - 2;         /* Do not account for anchor vertices in diffusion computations */
  if (grafptr->s.edloloctax != NULL) {
    for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) {
      Gnum                edgelocnum;
      Gnum                edgelocnnd;
      Gnum                edlolocsum;

#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((grafptr->s.vendloctax[vertlocnum] - grafptr->s.vertloctax[vertlocnum]) == 0) {
        errorPrint ("bdgraphBipartDf: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      difogsttax[vertlocnum] = 0.0F;
      for (edgelocnum = grafptr->s.vertloctax[vertlocnum], edgelocnnd = grafptr->s.vendloctax[vertlocnum], edlolocsum = 0;
           edgelocnum < edgelocnnd; edgelocnum ++)
        edlolocsum += grafptr->s.edloloctax[edgelocnum];

      ielsloctax[vertlocnum] = 1.0F / (float) edlolocsum;
    }
  }
  else {                                          /* Graph has no edge loads */
    for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) {
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((grafptr->s.vendloctax[vertlocnum] - grafptr->s.vertloctax[vertlocnum]) == 0) {
        errorPrint ("bdgraphBipartDf: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      ielsloctax[vertlocnum] = 1.0F / (float) (grafptr->s.vendloctax[vertlocnum] - grafptr->s.vertloctax[vertlocnum]);
      difogsttax[vertlocnum] = 0.0F;
    }
  }
  ielsloctax[vanclocnnd]     = 1.0F / (float) reduglbtab[0];
  ielsloctax[vanclocnnd + 1] = 1.0F / (float) reduglbtab[1];
  difogsttax[vanclocnnd]     = vanclocval[0] * ielsloctax[vanclocnnd]; /* Load anchor vertices for first pass */
  difogsttax[vanclocnnd + 1] = vanclocval[1] * ielsloctax[vanclocnnd + 1];
  difngsttax[vanclocnnd]     =                    /* In case of isolated anchors, do not risk overflow because of NaN */
  difngsttax[vanclocnnd + 1] = 0.0F;

  if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (difogsttax + grafptr->s.baseval), MPI_FLOAT) != 0) { /* Perform initial diffusion (and build communication structures) */
    errorPrint ("bdgraphBipartDf: cannot propagate diffusion data (1)");
    memFree    (ielsloctax + grafptr->s.baseval); /* Free group leader */
    return     (1);
  }

  ovflval    = 0;
  cdifval    = paraptr->cdifval;
  cremval    = paraptr->cremval;
  edgegsttax = grafptr->s.edgegsttax;
  for (passnum = 0; ; ) {                         /* For all passes         */
    if (ovflval == 0) {                           /* If no overflow occured */
      float *             diftgsttax;             /* Temporary swap value   */
      Gnum                vertlocnum;
      float               veloval;

      veloval = 1.0F;                             /* Assume no vertex loads */
      for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) {
        Gnum                edgelocnum;
        Gnum                edgelocnnd;
        float               diffval;

        diffval    = 0.0F;
        edgelocnum = grafptr->s.vertloctax[vertlocnum];
        edgelocnnd = grafptr->s.vendloctax[vertlocnum];
        if (grafptr->s.edloloctax != NULL)
          for ( ; edgelocnum < edgelocnnd; edgelocnum ++)
            diffval += difogsttax[edgegsttax[edgelocnum]] * (float) grafptr->s.edloloctax[edgelocnum];
        else
          for ( ; edgelocnum < edgelocnnd; edgelocnum ++)
            diffval += difogsttax[edgegsttax[edgelocnum]];

        diffval *= cdifval;
        diffval += (difogsttax[vertlocnum] * cremval) / ielsloctax[vertlocnum];

        if (grafptr->s.veloloctax != NULL)
          veloval = (float) grafptr->s.veloloctax[vertlocnum];
        if (diffval >= 0.0F) {
          diffval = (diffval - veloval) * ielsloctax[vertlocnum];
          if (diffval <= 0.0F)
            diffval = +BDGRAPHBIPARTDFEPSILON;
        }
        else {
          diffval = (diffval + veloval) * ielsloctax[vertlocnum];
          if (diffval >= 0.0F)
            diffval = -BDGRAPHBIPARTDFEPSILON;
        }
        if (isnan (diffval)) {                    /* If overflow occured                    */
          ovflval = 1;                            /* We are in state of overflow            */
          goto abort;                             /* Exit this loop without swapping arrays */
        }
        difngsttax[vertlocnum] = diffval;
      }
      for ( ; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) { /* For the two local anchor vertices */
        Gnum                edgelocnum;
        Gnum                edgelocnnd;
        float               diffval;

        diffval    = 0.0F;
        edgelocnum = grafptr->s.vertloctax[vertlocnum] + grafptr->s.procglbnbr - 1; /* Skip links to other anchors */
        edgelocnnd = grafptr->s.vendloctax[vertlocnum];
        if (edgelocnum == edgelocnnd)             /* If isolated anchor */
          continue;                               /* Barrel is empty    */

        for ( ; edgelocnum < edgelocnnd; edgelocnum ++) /* Anchor edges have load 1 even for weighted graphs */
          diffval += difogsttax[edgegsttax[edgelocnum]];

        diffval *= cdifval;
        diffval += vanclocval[vertlocnum - vanclocnnd] + (difogsttax[vertlocnum] * cremval) / ielsloctax[vertlocnum];
        if (diffval >= 0.0F) {
          diffval = (diffval - valolocval[vertlocnum - vanclocnnd]) * ielsloctax[vertlocnum];
          if (diffval <= 0.0F)
            diffval = +BDGRAPHBIPARTDFEPSILON;
        }
        else {
          diffval = (diffval + valolocval[vertlocnum - vanclocnnd]) * ielsloctax[vertlocnum];
          if (diffval >= 0.0F)
            diffval = -BDGRAPHBIPARTDFEPSILON;
        }
        if (isnan (diffval)) {                    /* If overflow occured                    */
          ovflval = 1;                            /* We are in state of overflow            */
          goto abort;                             /* Exit this loop without swapping arrays */
        }
        difngsttax[vertlocnum] = diffval;
      }

      diftgsttax = (float *) difngsttax;          /* Swap old and new diffusion arrays          */
      difngsttax = (float *) difogsttax;          /* Casts to prevent IBM compiler from yelling */
      difogsttax = (float *) diftgsttax;
    }
abort :                                           /* If overflow occured, resume here    */
    if (++ passnum >= paraptr->passnbr)           /* If maximum number of passes reached */
      break;                                      /* Exit main loop                      */

    if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (difogsttax + grafptr->s.baseval), MPI_FLOAT) != 0) {
      errorPrint ("bdgraphBipartDf: cannot propagate diffusion data (2)");
      memFree    (ielsloctax + grafptr->s.baseval); /* Free group leader */
      return     (1);
    }
  }

  for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) /* Set new part distribution */
    grafptr->partgsttax[vertlocnum] = (difogsttax[vertlocnum] <= 0.0F) ? 0 : 1;
  grafptr->partgsttax[vanclocnnd]     = 0;        /* Set up parts in case anchors are isolated */
  grafptr->partgsttax[vanclocnnd + 1] = 1;

  memFree (ielsloctax + grafptr->s.baseval);      /* Free group leader */

  if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (grafptr->partgsttax + grafptr->s.baseval), GRAPHPART_MPI) != 0) {
    errorPrint ("bdgraphBipartDf: cannot propagate part data (2)");
    return     (1);
  }

  commlocloadintn =
  commlocloadextn =
  commlocgainextn = 0;
  for (vertlocnum = grafptr->s.baseval, fronlocnum = complocsize1 = complocload1 = 0;
       vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
    Gnum                edgelocnum;
    Gnum                edgelocnnd;
    Gnum                veloval;
    Gnum                partval;
    Gnum                flagval;

#ifdef SCOTCH_DEBUG_BDGRAPH2
    if (grafptr->partgsttax[vertlocnum] > 1) {
      errorPrint ("bdgraphBipartDf: internal error (3)");
      break;                                      /* Do not break upcoming collective communications */
    }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    partval = (Gnum) grafptr->partgsttax[vertlocnum];
    veloval = (veloloctax != NULL) ? veloloctax[vertlocnum] : 1;
    if (grafptr->veexloctax != NULL) {
      commlocloadextn += grafptr->veexloctax[vertlocnum] * partval;
      commlocgainextn += grafptr->veexloctax[vertlocnum] * (1 - partval * 2);
    }
    complocsize1 += partval;
    complocload1 += partval * veloval;

    flagval = 0;
    for (edgelocnum = grafptr->s.vertloctax[vertlocnum], edgelocnnd = grafptr->s.vendloctax[vertlocnum];
         edgelocnum < edgelocnnd; edgelocnum ++) {
      Gnum                edloval;
      Gnum                partend;

      partend = (Gnum) grafptr->partgsttax[edgegsttax[edgelocnum]];
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if (partend > 1) {
        errorPrint ("bdgraphBipartDf: internal error (4)");
        vertlocnum = grafptr->s.vertlocnnd;
        break;                                    /* Do not break upcoming collective communications */
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      edloval  = (edloloctax != NULL) ? edloloctax[edgelocnum] : 1;
      flagval |= partval ^ partend;
      commlocloadintn += (partval ^ partend) * edloval; /* Internal load is accounted for twice */
    }
    if (flagval != 0)                             /* If vertex has neighbors in other part */
      grafptr->fronloctab[fronlocnum ++] = vertlocnum; /* Record it as member of separator */
  }
  grafptr->complocload0 = grafptr->s.velolocsum - complocload1;
  grafptr->complocsize0 = grafptr->s.vertlocnbr - complocsize1;
  grafptr->fronlocnbr   = fronlocnum;

  reduloctab[0] = fronlocnum;
  reduloctab[1] = grafptr->complocload0;
  reduloctab[2] = grafptr->complocsize0;
  reduloctab[3] = commlocloadintn;                /* Twice the internal load; sum globally before dividing by two */
  reduloctab[4] = commlocloadextn;
  reduloctab[5] = commlocgainextn;
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 6, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartDf: communication error (2)");
    return     (1);
  }
  grafptr->fronglbnbr      = reduglbtab[0];
  grafptr->compglbload0    = reduglbtab[1];
  grafptr->compglbload0dlt = grafptr->compglbload0 - grafptr->compglbload0avg;
  grafptr->compglbsize0    = reduglbtab[2];
  grafptr->commglbload     = (reduglbtab[3] / 2) * grafptr->domndist + reduglbtab[4];
  grafptr->commglbgainextn = reduglbtab[5];
  grafptr->bbalglbval      = (double) ((grafptr->compglbload0dlt < 0) ? (- grafptr->compglbload0dlt) : grafptr->compglbload0dlt) / (double) grafptr->compglbload0avg;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (grafptr) != 0) {
    errorPrint ("bdgraphBipartDf: internal error (5)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}
