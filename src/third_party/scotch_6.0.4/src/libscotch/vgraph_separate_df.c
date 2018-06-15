/* Copyright 2007,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_df.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a separator of     **/
/**                a separation graph by using a diffusion **/
/**                scheme.                                 **/
/**                                                        **/
/**   NOTES      : # This algorithm has been designed to   **/
/**                  work on band graphs only, for which   **/
/**                  the two anchor vertices are the two   **/
/**                  last vertices, the before-last as     **/
/**                  anchor of part 0, and the last as     **/
/**                  anchor of part 1.                     **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 29 oct 2007     **/
/**                                 to     24 may 2008     **/
/**                # Version 6.0  : from : 24 dec 2013     **/
/**                                 to     24 dec 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_DF

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_df.h"

/*
**  The static variables.
*/

static const Gnum           vgraphseparatedfloadone = 1;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the separation.
** It returns:
** - 0 : if the separator could be computed.
** - 1 : on error.
*/

int
vgraphSeparateDf (
Vgraph * restrict const             grafptr,      /*+ Active graph      +*/
const VgraphSeparateDfParam * const paraptr)      /*+ Method parameters +*/
{
  float * restrict      edlstax;                  /* Degree array              */
  float * restrict      difotax;                  /* Old diffusion value array */
  float * restrict      difntax;                  /* New diffusion value array */
  float * restrict      difttax;                  /* Temporary swap value      */
  float                 cdifval;
  float                 cremval;
  Gnum                  fronnum;
  Gnum                  compload0avg;
  Gnum                  compload2;
  Gnum                  compsize1;
  float                 veloval;
  float                 vanctab[2];
  INT                   movenum;
  INT                   passnum;

  const Gnum * restrict const verttax = grafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  Gnum *       restrict const frontab = grafptr->frontab;
  GraphPart *  restrict const parttax = grafptr->parttax;

  if (memAllocGroup ((void **) (void *)
                     &edlstax, (size_t) (grafptr->s.vertnbr * sizeof (float)),
                     &difotax, (size_t) (grafptr->s.vertnbr * sizeof (float)),
                     &difntax, (size_t) (grafptr->s.vertnbr * sizeof (float)), NULL) == NULL) {
    errorPrint ("vgraphSeparateDf: out of memory");
    return     (1);
  }
  edlstax -= grafptr->s.baseval;                  /* Base access to edlstax and diffusion arrays */
  difotax -= grafptr->s.baseval;
  difntax -= grafptr->s.baseval;

  if (grafptr->s.edlotax == NULL) {               /* If graph has no edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++)
      edlstax[vertnum] = (float) (vendtax[vertnum] - verttax[vertnum]);
  }
  else {                                          /* If graph has edge weights */
    Gnum                vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                edgenum;
      Gnum                edlosum;

      for (edgenum = verttax[vertnum], edlosum = 0;
           edgenum < vendtax[vertnum]; edgenum ++)
        edlosum += grafptr->s.edlotax[edgenum];

      edlstax[vertnum] = (float) edlosum;
    }
  }

  compload0avg = grafptr->compload[0] + grafptr->compload[2] / 2;

  passnum = 0;
  do {
    const Gnum * restrict velobax;
    Gnum                  velomsk;
    Gnum                  vertnum;
    Gnum                  compload0;
    Gnum                  compload1;
    int                   rootval;                /* Root part for separator vertices */
    
    compload0  = compload0avg - grafptr->compload[2] / 2;
    compload1  = grafptr->s.velosum - compload0avg - (grafptr->compload[2] + 1) / 2;
    vanctab[0] = (float) (- compload0);           /* Values to be injected to anchor vertices at every iteration  */
    vanctab[1] = (float)    compload1 - VGRAPHSEPARATEDFEPSILON; /* Slightly tilt value to add to part 1          */
    rootval    = (paraptr->partval + passnum) & 1; /* Compute index of part which will receive separator vertices */
    if (rootval == 0)                             /* If separator must be aggregated to part 0                    */
      vanctab[0] -= (float) grafptr->compload[2];
    else
      vanctab[1] += (float) grafptr->compload[2];

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd - 2; vertnum ++)
      difotax[vertnum] = 0.0F;
    difotax[grafptr->s.vertnnd - 2] = vanctab[0] / edlstax[grafptr->s.vertnnd - 2]; /* Load anchor vertices for first move */
    difotax[grafptr->s.vertnnd - 1] = vanctab[1] / edlstax[grafptr->s.vertnnd - 1];

    veloval = 1.0F;                               /* Assume no vertex loads */
    cdifval = paraptr->cdifval;
    cremval = paraptr->cremval;
    for (movenum = 0; movenum < paraptr->movenbr; movenum ++) { /* For all moves */
      Gnum                vertnum;
      Gnum                vertnnd;
      float               vancval;                /* Value to load vertex with if anchor */

      vancval = 0.0F;                             /* At first vertices are not anchors */
      vertnum = grafptr->s.baseval;
      vertnnd = grafptr->s.vertnnd - 2;
      while (1) {
        for ( ; vertnum < vertnnd; vertnum ++) {
          Gnum                edgenum;
          Gnum                edgennd;
          float               diffval;

          edgenum = verttax[vertnum];
          edgennd = vendtax[vertnum];
          diffval = 0.0F;
          if (grafptr->s.edlotax != NULL)
            for ( ; edgenum < edgennd; edgenum ++)
              diffval += difotax[edgetax[edgenum]] * (float) grafptr->s.edlotax[edgenum];
          else
            for ( ; edgenum < edgennd; edgenum ++)
              diffval += difotax[edgetax[edgenum]];

          if (grafptr->s.velotax != NULL)
            veloval = (float) grafptr->s.velotax[vertnum];
          diffval *= cdifval;
          diffval += difotax[vertnum] * cremval * edlstax[vertnum];
          if (diffval >= 0.0F) {
            diffval -= veloval;
            if (diffval <= 0.0F)
              diffval = +VGRAPHSEPARATEDFEPSILON;
          }
          else {
            diffval += veloval;
            if (diffval >= 0.0F)
              diffval = -VGRAPHSEPARATEDFEPSILON;
          }
          if (isnan (diffval))                    /* If overflow occured                                                       */
            goto abort;                           /* Exit main loop without swapping arrays so as to keep last valid iteration */

          difntax[vertnum] = diffval / edlstax[vertnum]; /* Prepare vertex for diffusion */
        }
        if (vertnum == grafptr->s.vertnnd)        /* If all vertices processed, exit intermediate infinite loop */
          break;

        vertnnd ++;                               /* Prepare to go only for one more run        */
        vancval = vanctab[vertnum - grafptr->s.vertnnd + 2]; /* Load variable with anchor value */
      }

      difttax = difntax;                          /* Swap old and new diffusion arrays */
      difntax = difotax;
      difotax = difttax;
    }
abort :                                           /* If overflow occured, resume here */

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) /* Pre-set parts without separator */
      parttax[vertnum] = (difotax[vertnum] <= 0.0F) ? 0 : 1;

    if (grafptr->s.velotax != NULL) {
      velobax = grafptr->s.velotax;
      velomsk = ~((Gnum) 0);
    }
    else {
      velobax = &vgraphseparatedfloadone;
      velomsk = 0;
    }

    veloval = 1.0F;                               /* Assume no vertex loads */
    for (vertnum = grafptr->s.baseval, fronnum = compsize1 = compload1 = compload2 = 0;
         vertnum < grafptr->s.vertnnd; vertnum ++) {
      Gnum                partval;
      GraphPart           partend;
      Gnum                veloval;

      partend = parttax[vertnum] ^ 1;
      partval = (Gnum) parttax[vertnum];
      veloval = velobax[vertnum & velomsk];
      compsize1 += partval;                       /* Here, part is 0 or 1 only */
      compload1 += partval * veloval;
      if (partval == (Gnum) rootval) {            /* Only vertices of aggregated part can be in separator */
        Gnum                edgenum;

        for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
          if (parttax[edgetax[edgenum]] == partend) { /* If end vertex is in other part (and not in separator) */
            frontab[fronnum ++] = vertnum;        /* Record it                                                 */
            parttax[vertnum]    = 2;
            compload2 += veloval;
            break;                                /* No need to go further */
          }
        }
      }
    }
    grafptr->compload[0] = grafptr->s.velosum - compload1;
    grafptr->compload[1] = compload1;
    grafptr->compload[2] = compload2;
    grafptr->compload[rootval] -= compload2;
    grafptr->comploaddlt = grafptr->compload[0] - grafptr->compload[1];
    grafptr->compsize[0] = grafptr->s.vertnbr - compsize1;
    grafptr->compsize[1] = compsize1;
    grafptr->compsize[rootval] -= fronnum;
    grafptr->fronnbr     = fronnum;
  } while (++ passnum < paraptr->passnbr);        /* As long as not all passes performed */

  memFree (edlstax + grafptr->s.baseval);         /* Free group leader */

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (grafptr) != 0) {
    errorPrint ("vgraphSeparateDf: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  return (0);
}
