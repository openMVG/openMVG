/* Copyright 2007-2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : wgraph_part_gg.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : This module is the vertex overlapped    **/
/**                graph partitioning rountine based on    **/
/**                a vertex-oriented version of the Greedy **/
/**                Graph Growing algorithm.                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 14 mar 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_GG

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "wgraph.h"
#include "wgraph_part_gg.h"

/*
**  The static variables.
*/

static const Gnum           wgraphpartggloadone = 1;

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
wgraphPartGg (
Wgraph * restrict const         wgrafptr,   /*+ Separation graph  +*/
const WgraphPartGgParam * const paraptr)    /*+ Method parameters +*/
{
  Gnum                              i;
  Gnum                              gain;
  Gnum                              fronnbr;     /* Number of frontier vertices */
  Gnum                              fronload;    /* Load of frontier vertices   */
  Gnum                              frlobst;     /* Best frontier load found    */
  Gnum                              fronnum;
  Gnum                              partval;
  Gnum                              vertnum;
  Gnum                              vertnum2;
  Gnum                              vertnum3;
  Gnum                              edgenum;
  Gnum                              edgenum2;
  Gnum                              palooth;     /* Load of vertices in the remained unassigned part            */ 
  Gnum                              paloexc;     /* Load of vertices in the current part excepting the frontier */
  Gnum                              frloprt;     /* Load of vertices in the frontier of the current part        */
  Gnum                              passnum;
  Gnum                              velomsk;
  const Gnum * restrict             velobax;     /* Data for handling of optional arrays  */
  Gnum * restrict                   permtab;     /* permutation table */
  Gnum * restrict                   parttax;
  Gnum * restrict                   compload;    /* Array of part load            */
  Gnum * restrict                   compsize;    /* Array of part vertices number */
  GainTabl * restrict               tabl;        /* Pointer to gain table                 */
  WgraphPartGgVertex *              vexxtax;     /* Complementary vertex array            */
  GainLink * restrict               gainlinkptr;
  WgraphPartGgVertex * restrict     vertlist;    /* List of vertices                                               */

  printf ("GG (" GNUMSTRING ")\n", wgrafptr->s.vertnbr);

  if (((tabl = gainTablInit (GAIN_LINMAX, WGRAPHSEPAGGSUBBITS)) == NULL) || /* Use logarithmic array only */
      memAllocGroup((void **) (void *)
                    &vexxtax,  (size_t) (wgrafptr->s.vertnbr * sizeof (WgraphPartGgVertex)),
                    &compload, (size_t) (wgrafptr->partnbr * sizeof (Gnum)),
                    &compsize, (size_t) (wgrafptr->partnbr * sizeof (Gnum)),
                    &parttax,  (size_t) (wgrafptr->s.vertnbr * sizeof (Gnum)),
		    &permtab,  (size_t) (wgrafptr->s.vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("wgraphPartGg: out of memory (1)");
    gainTablExit (tabl);
    return (1);
  }

  vexxtax -= wgrafptr->s.baseval;                  /* Base access to vexxtax */
  parttax -= wgrafptr->s.baseval;                  /* Base access to parttax */


  for (vertnum = 0; vertnum < wgrafptr->s.vertnbr; vertnum ++) { /* Initialization of the permutation table */
    i = intRandVal (vertnum + 1);
    permtab[vertnum] = permtab[i];
    permtab[i] = vertnum;
  }

  if (wgrafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &wgraphpartggloadone;                /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = wgrafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  frlobst = wgrafptr->s.velosum + 1;
  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) {

    fronload     =
    fronnbr      = 0;
    vertnum      = -1;
    palooth      = wgrafptr->s.velosum;

    memSet (compload, 0, wgrafptr->partnbr * sizeof (Gnum));
    memSet (compsize, 0, wgrafptr->partnbr * sizeof (Gnum));
    memSet (parttax + wgrafptr->s.baseval, 0, wgrafptr->s.vertnbr * sizeof (Gnum));
    memSet (vexxtax + wgrafptr->s.baseval, 0, wgrafptr->s.vertnbr * sizeof (WgraphPartGgVertex));

    gainTablFree (tabl);

    for (partval = 0; partval < wgrafptr->partnbr; partval ++) {

      paloexc =
      frloprt = 0;

      gainlinkptr = gainTablFrst (tabl);             /* Try to take a vertex from the frontier of the last part */
      gainTablFree (tabl);
      if (gainlinkptr != NULL) {                     /* If the table was not empty */
        vertnum = (WgraphPartGgVertex *) gainlinkptr - vexxtax;
	
        for (edgenum = wgrafptr->s.verttax[vertnum]; /* search a neighbor vertex that is in the part 0 */
             edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
          vertnum2 = wgrafptr->s.edgetax[edgenum];
          if (parttax[vertnum2] == 0) {
            gainTablAdd(tabl, &(vexxtax[vertnum2].gainlink), 0); /* Add it in the table in order to be selected */
                                                                 /* (the table will contain only this vertex)   */
	    
            fronnbr  ++;
            parttax[vertnum2]  = -1;                  /* Move it in the seperator */
            fronload          += velobax[vertnum2 & velomsk];
            compload[partval] += velobax[vertnum2 & velomsk];
            compsize[partval] ++;
            break;
          }
        }
      }

      do {                                          /* While the part is not big enought */

        gainlinkptr = gainTablFrst(tabl);
        if (gainlinkptr != NULL) {
          vertnum = (WgraphPartGgVertex *) gainlinkptr - vexxtax;
          gainTablDel(tabl, gainlinkptr);
        }
        else {                                      /* If the table was empty */
          if ((2 * (paloexc + frloprt / 2)) * (wgrafptr->partnbr - partval + 1) <= palooth) { /* If the part load is not big enought */
	    for (i = 0; i < wgrafptr->s.vertnbr; i ++) { /* select a random vertex in the part 0 */
	      Gnum pos = i + intRandVal (wgrafptr->s.vertnbr - i);

	      vertnum = permtab[pos];
	      permtab[pos] = permtab[i];
	      permtab[i] = vertnum;
              vertnum += wgrafptr->s.baseval;
              if (parttax[vertnum] == 0) {
                for (edgenum = wgrafptr->s.verttax[vertnum]; edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
                  vertnum2 = wgrafptr->s.edgetax[edgenum];
                  if (parttax[vertnum2] > 0) {
                    break;
                  }
                }
		break;
	      }
	    }
            fronload          += velobax[vertnum & velomsk];
            frloprt           += velobax[vertnum & velomsk];
            compload[partval] += velobax[vertnum & velomsk];
            parttax[vertnum]   = -1;
            compsize[partval] ++;
            fronnbr           ++;
          }
          else
            break;
        }

        fronnbr --;
        frloprt  -= velobax[vertnum & velomsk];
        fronload -= velobax[vertnum & velomsk];
        paloexc  += velobax[vertnum & velomsk];
        parttax[vertnum] = partval;
        vertlist         = NULL;

        for (edgenum = wgrafptr->s.verttax[vertnum]; edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
          vertnum2 = wgrafptr->s.edgetax[edgenum];
          if ((parttax[vertnum2] != -1) && (parttax[vertnum2] != partval)) { /* If the vertex is in a different part */
            fronnbr ++;
            parttax[vertnum2]  = -1;                                     /* Move the vertex in the separator */
            fronload          += velobax[vertnum2 & velomsk];
            frloprt           += velobax[vertnum2 & velomsk];
            compload[partval] += velobax[vertnum2 & velomsk];
            compsize[partval] ++;

            vexxtax[vertnum2].partlvl = partval;                         /* Label the vertex for this part */
            vexxtax[vertnum2].prev    = vertlist;                        /* Add the vertex in the list */
            vertlist                  = vexxtax + vertnum2;

          } else if (parttax[vertnum2] == -1) {
            if (vexxtax[vertnum2].partlvl == partval) {                  /* If the vertex is labeled for the current part */
              if (vexxtax[vertnum2].gainlink.next != (GainLink * )~0) {    /* If the vertex is in the table */
                gainTablDel(tabl, &(vexxtax[vertnum2].gainlink));        /* Remove it from table */

                vexxtax[vertnum2].prev = vertlist;                       /* Add the vertex in the list */
                vertlist = vexxtax + vertnum2;
              }
            }
            else {
              frloprt           += velobax[vertnum2 & velomsk];
              compload[partval] += velobax[vertnum2 & velomsk];
              compsize[partval] ++;
              vexxtax[vertnum2].partlvl       = partval;                 /* Label it for the part */
              vexxtax[vertnum2].gainlink.next = (GainLink * )~0;
            }
          }
        }

        while (vertlist != NULL) {                                       /* For eack linked vertices */

          vertnum2 = vertlist - vexxtax;
          gain     = - velobax[vertnum2 & velomsk];                          /* Compute gain */
          for (edgenum2 = wgrafptr->s.verttax[vertnum2];
               edgenum2 < wgrafptr->s.vendtax[vertnum2]; edgenum2 ++) {
            vertnum3 = wgrafptr->s.edgetax[edgenum2];
            if ((parttax[vertnum3] != -1) && (parttax[vertnum3] != partval)) {
              gain += velobax[vertnum3 & velomsk];
            }
          }

          gainTablAdd(tabl, &(vertlist->gainlink), gain);                /* Add the vertex in the table */
          vertlist = vertlist->prev;
        }
      } while ((paloexc + frloprt / 2) * (wgrafptr->partnbr - partval + 1) <= palooth); /* While the part is not big enought */
      
      palooth -= (paloexc + frloprt / 2);
    }

    compload[0] =
    compsize[0] = 0;
    for (vertnum = wgrafptr->s.baseval; vertnum < wgrafptr->s.vertnnd; vertnum++) { /* Recompute load and size of part 0 */
      if (parttax[vertnum] == 0) {
        compload[0] += velobax[vertnum & velomsk];
        compsize[0] ++;
      } 
      else if (parttax[vertnum] == -1) {
        for (edgenum = wgrafptr->s.verttax[vertnum];
             edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
          vertnum2 = wgrafptr->s.edgetax[edgenum];
          if (parttax[vertnum2] == 0) {
            compload[0] += velobax[vertnum & velomsk];
            compsize[0] ++;
            break;
          }
        }
      }
    }

    if (frlobst > fronload) {                    /* If the pass frontier load is better than the better one */
      
      wgrafptr->fronnbr  = fronnbr;
      wgrafptr->fronload = fronload;
      memCpy (wgrafptr->compload, compload, sizeof (Gnum) * wgrafptr->partnbr);
      memCpy (wgrafptr->compsize, compsize, sizeof (Gnum) * wgrafptr->partnbr);
      memCpy (wgrafptr->parttax + wgrafptr->s.baseval, parttax + wgrafptr->s.baseval, sizeof (Gnum) * wgrafptr->s.vertnbr);
      
      
      for (vertnum = wgrafptr->s.baseval, fronnum = 0; vertnum < wgrafptr->s.vertnnd; vertnum ++) { /* Recompute frontab */
        if (parttax[vertnum] == -1)
          wgrafptr->frontab[fronnum ++] = vertnum;
      }

      frlobst = fronload;                       /* This frontier load will the best found */
    }
  }

  for (partval = 0; partval < wgrafptr->partnbr; partval ++)
    printf("\033[0;33mcompload[" GNUMSTRING "] " GNUMSTRING " " GNUMSTRING "\033[0m\n", partval, wgrafptr->compload[partval], wgrafptr->compsize[partval]);
  memFree(vexxtax + wgrafptr->s.baseval);       /* Free work arrays */
  gainTablExit (tabl);

  return (0);
}
