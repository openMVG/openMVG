/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph_separate_df.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a separator of the **/
/**                given distributed separator graph by    **/
/**                applying a diffusion method to what is  **/
/**                assumed to be a distributed band graph. **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 05 nov 2007     **/
/**                                 to   : 09 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_SEPARATE_DF

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "vdgraph.h"
#include "vdgraph_separate_df.h"

/*
**  The static variables.
*/

static const Gnum           vdgraphseparatedfloadone = 1;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed separator
** by diffusion across what is assumed to be a
** distributed band graph.
** It returns:
** - 0   : if the separator could be computed.
** - !0  : on error.
*/

int
vdgraphSeparateDf (
Vdgraph * const                       grafptr,    /*+ Distributed graph +*/
const VdgraphSeparateDfParam * const  paraptr)    /*+ Method parameters +*/
{
  float * restrict                ielsloctax;     /* Inverse of degree array   */
  float * restrict                difogsttax;     /* Old diffusion value array */
  float * restrict                difngsttax;     /* New diffusion value array */
  const Gnum * restrict           edgegsttax;
  Gnum                            fronlocnum;
  float                           compglbavg;
  Gnum                            complocload1;
  Gnum                            complocload2;
  Gnum                            complocsize1;
  const Gnum * restrict           velolocbax;
  Gnum                            velolocmsk;
  float                           vanclocval[2];
  float                           valolocval[2];  /* Fraction of load to remove from anchor vertices at each step */
  Gnum                            vanclocnnd;
  Gnum                            vertlocnum;
  Gnum                            reduloctab[4];  /* Local degree of both anchor vertices, minus edges to other anchors, and their loads */
  Gnum                            reduglbtab[4];
  Gnum                            passnum;
  float                           cdifval;
  float                           cremval;
  Gnum                            psepval;        /* Separator part      */
  int                             ovflval;        /* Overflow flag value */

  if (dgraphGhst (&grafptr->s) != 0) {            /* Compute ghost edge array if not already present */
    errorPrint ("vdgraphSeparateDf: cannot compute ghost edge array");
    return     (1);
  }

  reduloctab[0] = grafptr->s.vendloctax[grafptr->s.vertlocnnd - 2] - grafptr->s.vertloctax[grafptr->s.vertlocnnd - 2] - (grafptr->s.procglbnbr - 1);
  reduloctab[1] = grafptr->s.vendloctax[grafptr->s.vertlocnnd - 1] - grafptr->s.vertloctax[grafptr->s.vertlocnnd - 1] - (grafptr->s.procglbnbr - 1);
  if (grafptr->s.veloloctax == NULL)
    reduloctab[2] = 
    reduloctab[3] = 1;
  else {
    reduloctab[2] = grafptr->s.veloloctax[grafptr->s.vertlocnnd - 2];
    reduloctab[3] = grafptr->s.veloloctax[grafptr->s.vertlocnnd - 1];
  }

  if (memAllocGroup ((void **) (void *)
                     &ielsloctax, (size_t) (grafptr->s.vertlocnbr * sizeof (float)),
                     &difogsttax, (size_t) (grafptr->s.vertgstnbr * sizeof (float)),
                     &difngsttax, (size_t) (grafptr->s.vertgstnbr * sizeof (float)), NULL) == NULL) {
    errorPrint ("vdgraphSeparateDf: out of memory");
    reduloctab[0] = -1;
  }
  else {
    ielsloctax -= grafptr->s.baseval;
    difogsttax -= grafptr->s.baseval;
    difngsttax -= grafptr->s.baseval;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 4, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateDf: communication error (1)");
    return     (1);
  }

  if (reduglbtab[0] < 0) {                        /* If memory error */
    if (ielsloctax != NULL)
      memFree (ielsloctax + grafptr->s.baseval);  /* Free group leader */
  }
  if ((reduglbtab[0] == 0) ||                     /* If graph is too small to have any usable anchors, leave partition as is */
      (reduglbtab[1] == 0)) {
    memFree (ielsloctax + grafptr->s.baseval);    /* Free group leader */
    return  (0);
  }

  psepval       = paraptr->partval & 1;           /* Coerce part in the {0,1} range */
  compglbavg    = (float) (grafptr->compglbload[0] + grafptr->compglbload[1]) * 0.5F;
  vanclocval[0] = (float) grafptr->compglbload[0];
  if (vanclocval[0] < (compglbavg * (1.0F - (float) paraptr->deltval))) /* Enforce balance constraint */
    vanclocval[0] = compglbavg * (1.0F - (float) paraptr->deltval);
  else if (vanclocval[0] > (compglbavg * (1.0F + (float) paraptr->deltval)))
    vanclocval[0] = compglbavg * (1.0F + (float) paraptr->deltval);
  vanclocval[1] = (float) (grafptr->compglbload[0] + grafptr->compglbload[1]) - vanclocval[0];
  valolocval[0] = (float) reduglbtab[2];          /* Compute values to remove from anchor vertices */
  valolocval[1] = (float) reduglbtab[3];
  if (vanclocval[0] < valolocval[0])              /* If anchor in part 0 too large to reduce imbalance        */
    psepval = 1;                                  /* Separator must be taken from part 1 to stick to anchor 0 */
  else if (vanclocval[1] < valolocval[1])         /* Else if anchor in part 1 too large to reduce imbalance   */
    psepval = 0;                                  /* It is from part 0 that separator must be extracted       */
  vanclocval[psepval] += (float) grafptr->compglbload[2]; /* Aggregate separator to proper part               */
  vanclocval[0]  = - vanclocval[0];               /* Part 0 holds negative values                             */
  vanclocval[1] -= VDGRAPHSEPARATEDFEPSILON;      /* Slightly tilt value to add to part 1                     */

  for (vertlocnum = grafptr->s.baseval, vanclocnnd = grafptr->s.vertlocnnd - 2; /* Do not account for anchor vertices in diffusion computations */
       vertlocnum < vanclocnnd; vertlocnum ++) {
#ifdef SCOTCH_DEBUG_VDGRAPH2
    if ((grafptr->s.vendloctax[vertlocnum] - grafptr->s.vertloctax[vertlocnum]) == 0) {
      errorPrint ("vdgraphSeparateDf: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
    ielsloctax[vertlocnum] = 1.0F / (float) (grafptr->s.vendloctax[vertlocnum] - grafptr->s.vertloctax[vertlocnum]);
    difogsttax[vertlocnum] = 0.0F;
  }
  ielsloctax[vanclocnnd]     = 1.0F / (float) reduglbtab[0];
  ielsloctax[vanclocnnd + 1] = 1.0F / (float) reduglbtab[1];
  difogsttax[vanclocnnd]     = vanclocval[0] * ielsloctax[vanclocnnd]; /* Load anchor vertices for first pass */
  difogsttax[vanclocnnd + 1] = vanclocval[1] * ielsloctax[vanclocnnd + 1];
  difngsttax[vanclocnnd]     =                    /* In case of isolated anchors, do not risk overflow because of NaN */
  difngsttax[vanclocnnd + 1] = 0.0F;

  if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (difogsttax + grafptr->s.baseval), MPI_FLOAT) != 0) { /* Perform initial diffusion (and build communication structures) */
    errorPrint ("vdgraphSeparateDf: cannot propagate diffusion data (1)");
    memFree    (ielsloctax + grafptr->s.baseval); /* Free group leader */
    return     (1);
  }

  ovflval    = 0;
  cdifval    = paraptr->cdifval;
  cremval    = paraptr->cremval;
  edgegsttax = grafptr->s.edgegsttax;
  for (passnum = 0; ; ) {                         /* For all passes         */
    if (ovflval == 0) {                           /* If no overflow occured */
      Gnum                vertlocnum;
      float *             diftgsttax;             /* Temporary swap value */
      float               veloval;

      veloval = 1.0F;                             /* Assume no vertex loads */
      for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) {
        Gnum                edgelocnum;
        Gnum                edgelocnnd;
        float               diffval;

        diffval = 0.0F;
        for (edgelocnum = grafptr->s.vertloctax[vertlocnum], edgelocnnd = grafptr->s.vendloctax[vertlocnum];
             edgelocnum < edgelocnnd; edgelocnum ++)
          diffval += difogsttax[edgegsttax[edgelocnum]];

        diffval *= cdifval;
        diffval += (difogsttax[vertlocnum] * cremval) / ielsloctax[vertlocnum];

        if (grafptr->s.veloloctax != NULL)
          veloval = (float) grafptr->s.veloloctax[vertlocnum];
        if (diffval >= 0.0F) {
          diffval = (diffval - veloval) * ielsloctax[vertlocnum];
          if (diffval <= 0.0F)
            diffval = +VDGRAPHSEPARATEDFEPSILON;
        }
        else {
          diffval = (diffval + veloval) * ielsloctax[vertlocnum];
          if (diffval >= 0.0F)
            diffval = -VDGRAPHSEPARATEDFEPSILON;
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

        for ( ; edgelocnum < edgelocnnd; edgelocnum ++)
          diffval += difogsttax[edgegsttax[edgelocnum]];

        diffval *= cdifval;
        diffval += vanclocval[vertlocnum - vanclocnnd] + (difogsttax[vertlocnum] * cremval) / ielsloctax[vertlocnum];
        if (diffval >= 0.0F) {
          diffval = (diffval - valolocval[vertlocnum - vanclocnnd]) * ielsloctax[vertlocnum];
          if (diffval <= 0.0F)
            diffval = +VDGRAPHSEPARATEDFEPSILON;
        }
        else {
          diffval = (diffval + valolocval[vertlocnum - vanclocnnd]) * ielsloctax[vertlocnum];
          if (diffval >= 0.0F)
            diffval = -VDGRAPHSEPARATEDFEPSILON;
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
      errorPrint ("vdgraphSeparateDf: cannot propagate diffusion data (2)");
      memFree    (ielsloctax + grafptr->s.baseval); /* Free group leader */
      return     (1);
    }
  }

  for (vertlocnum = grafptr->s.baseval; vertlocnum < vanclocnnd; vertlocnum ++) /* Pre-set parts without separator */
    grafptr->partgsttax[vertlocnum] = (difogsttax[vertlocnum] <= 0.0F) ? 0 : 1;
  grafptr->partgsttax[vanclocnnd]     = 0;        /* Set up parts in case anchors are isolated */
  grafptr->partgsttax[vanclocnnd + 1] = 1;

  if (grafptr->s.veloloctax != NULL) {
    velolocbax = grafptr->s.veloloctax;
    velolocmsk = ~((Gnum) 0);
  }
  else {
    velolocbax = &vdgraphseparatedfloadone;
    velolocmsk = 0;
  }

  memFree (ielsloctax + grafptr->s.baseval);      /* Free group leader */

  if (dgraphHaloSync (&grafptr->s, (byte *) (void *) (grafptr->partgsttax + grafptr->s.baseval), GRAPHPART_MPI) != 0) {
    errorPrint ("vdgraphSeparateDf: cannot propagate part data");
    return     (1);
  }

  for (vertlocnum = grafptr->s.baseval, fronlocnum = complocsize1 = complocload1 = complocload2 = 0;
       vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
    Gnum                partval;
    GraphPart           partend;
    Gnum                veloval;

#ifdef SCOTCH_DEBUG_VDGRAPH2
    if (grafptr->partgsttax[vertlocnum] > 1) {
      errorPrint ("vdgraphSeparateDf: internal error (2)");
      break;                                      /* Do not break upcoming collective communications */
    }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
    partend =        grafptr->partgsttax[vertlocnum] ^ 1;
    partval = (Gnum) grafptr->partgsttax[vertlocnum];
    veloval = velolocbax[vertlocnum & velolocmsk];
    complocsize1 += partval;                      /* Here, part is 0 or 1 only */
    complocload1 += partval * veloval;
    if (partval == psepval) {                     /* Only vertices of aggregated part can be in separator */
      Gnum                edgelocnum;

      for (edgelocnum = grafptr->s.vertloctax[vertlocnum];
           edgelocnum < grafptr->s.vendloctax[vertlocnum]; edgelocnum ++) {
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (grafptr->partgsttax[edgegsttax[edgelocnum]] > 2) {
          errorPrint ("vdgraphSeparateDf: internal error (3)");
          vertlocnum = grafptr->s.vertlocnnd;
          break;                                  /* Do not break upcoming collective communications */
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        if (grafptr->partgsttax[edgegsttax[edgelocnum]] == partend) { /* If end vertex is in other part (and not in separator) */
          grafptr->fronloctab[fronlocnum ++] = vertlocnum; /* Record it as member of the separator                             */
          grafptr->partgsttax[vertlocnum]    = 2;
          complocload2 += veloval;
          break;                                  /* No need to go further */
        }
      }
    }
  }
  grafptr->complocload[0] = grafptr->s.velolocsum - complocload1;
  grafptr->complocload[1] = complocload1;
  grafptr->complocload[2] = complocload2;
  grafptr->complocload[psepval] -= complocload2;
  grafptr->complocsize[0] = grafptr->s.vertlocnbr - complocsize1;
  grafptr->complocsize[1] = complocsize1;
  grafptr->complocsize[psepval] -= fronlocnum;
  grafptr->complocsize[2] = fronlocnum;

  if (MPI_Allreduce (&grafptr->complocload[0], &grafptr->compglbload[0], 6, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) { /* TRICK: all arrays */
    errorPrint ("vdgraphSeparateDf: communication error (2)");
    return     (1);
  }
  grafptr->compglbloaddlt = grafptr->compglbload[0] - grafptr->compglbload[1];

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (grafptr) != 0) {
    errorPrint ("vdgraphSeparateDf: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

if (grafptr->s.proclocnum == 0)
fprintf (stderr, "BROL " GNUMSTRING "," GNUMSTRING "," GNUMSTRING "(" GNUMSTRING ")\n",
         (Gnum) grafptr->compglbload[0],
         (Gnum) grafptr->compglbload[1],
         (Gnum) grafptr->compglbload[2],
         (Gnum) grafptr->compglbloaddlt);

  return (0);
}
