/* Copyright 2007-2010,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_gp.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : This module is the vertex overlapped    **/
/**                graph partitioning rountine based on    **/
/**                a vertex-oriented version of the Gibbs- **/
/**                Poole-Stockmeyer algorithm.             **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 29 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_GP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "wgraph.h"
#include "wgraph_part_gp.h"

/*
**  The static variables.
*/

static const Gnum           wgraphpartgploadone = 1;

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
wgraphPartGp (
Wgraph * restrict const         wgrafptr,     /*+ Separation graph  +*/
const WgraphPartGpParam * const paraptr)      /*+ Method parameters +*/
{
  Gnum                              i;
  Gnum                              head;
  Gnum                              tail;
  Gnum                              fronnum;
  Gnum                              fronnbr;  /* number of frontier vertices */
  Gnum                              fronload; /* load of frontier vertices   */
  Gnum                              frlobst;  /* best frontier load found    */
  Gnum                              partval;
  Gnum                              vertnum;
  Gnum                              vertnum2;
  Gnum                              edgenum;
  Gnum                              palooth;  /* load of vertices in the remained unassigned part            */ 
  Gnum                              paloexc;  /* load of vertices in the current part excepting the frontier */
  Gnum                              frloprt;  /* load of vertices in the frontier of the current part        */
  Gnum                              passnum;
  Gnum                              velomsk;
  const Gnum * restrict             velobax;      /* Data for handling of optional arrays  */
  Gnum * restrict                   compload;     /* array of part load            */
  Gnum * restrict                   compsize;     /* array of part vertices number */
  Gnum * restrict                   permtab;      /* permutation table */
  Gnum * restrict                   stack;        /* Pointer to stack */
  Gnum * restrict                   parttax;
  WgraphPartGpVertex *              vexxtax;      /* Complementary vertex array            */
  WgraphPartGpVertex * restrict     vertlist;     /* list of vertices                                               */

  if (memAllocGroup((void **) (void *)
                    &vexxtax,  (size_t) (wgrafptr->s.vertnbr * sizeof (WgraphPartGpVertex)),
                    &compload, (size_t) (wgrafptr->partnbr * sizeof (Gnum)), 
                    &compsize, (size_t) (wgrafptr->partnbr * sizeof (Gnum)), 
                    &stack ,   (size_t) (wgrafptr->s.vertnbr * sizeof (Gnum)),
                    &parttax,  (size_t) (wgrafptr->s.vertnbr * sizeof (Gnum)),
		    &permtab,  (size_t) (wgrafptr->s.vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("wgraphPartGp: out of memory (1)");
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
    velobax = &wgraphpartgploadone;  /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = wgrafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  frlobst = wgrafptr->s.velosum + 1;           /* The first solution that will be found will be better than frlobst */

  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) {
    fronnbr  =                                 /* Begin with empty separators */
    fronload = 0;
    vertnum  = -1;
    palooth  = wgrafptr->s.velosum;

    memSet (compload, 0, wgrafptr->partnbr * sizeof (Gnum));
    memSet (compsize, 0, wgrafptr->partnbr * sizeof (Gnum));
    memSet (parttax + wgrafptr->s.baseval, 0, wgrafptr->s.vertnbr * sizeof (Gnum));
    memSet (vexxtax + wgrafptr->s.baseval, 0, wgrafptr->s.vertnbr * sizeof (WgraphPartGpVertex));

    head =
    tail = 0;
    for (partval = 0; partval < wgrafptr->partnbr; partval ++) {
      frloprt =
      paloexc = 0;
      if (tail > head)                            /* try to take a vertex from the frontier of the last part */
        vertnum = stack[tail % wgrafptr->s.vertnbr];
      else
        vertnum = -1;

      head =
      tail = 0;
      if (vertnum != -1) {                  /* if the stack was not empty */
        for (edgenum = wgrafptr->s.verttax[vertnum]; /* search a neighbor vertex that is in the part 0 */
             edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
          vertnum2 = wgrafptr->s.edgetax[edgenum];
          if (parttax[vertnum2] == 0) {
            fronnbr ++;
            compsize[partval] ++;
            parttax[vertnum2]  = -1;               /* move it in the seperator */
            fronload          += velobax[vertnum2 & velomsk]; /* update load and size of parts */
            compload[partval] += velobax[vertnum2 & velomsk];
            vexxtax[vertnum2].isinstack = 0;
            stack[(head ++) % (wgrafptr->s.vertnbr)] = vertnum2;
            break;
          }
        }
      }

      do {                                        /* while the part is not big enought */
        if (head > tail)                            /* select a vertex */
          vertnum = stack[tail % wgrafptr->s.vertnbr];
        else
          vertnum = -1;

        if (vertnum != -1) {
          tail ++;
        }
        else                                      /* if the stack was empty */
        {
          if ((2 * (paloexc + frloprt / 2)) * (wgrafptr->partnbr - partval + 1) <= palooth) { /* if the part load is not big enought */
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
            fronnbr ++;
            parttax[vertnum]   = -1;
            frloprt           += velobax[vertnum & velomsk]; 
            fronload          += velobax[vertnum & velomsk];
            compload[partval] += velobax[vertnum & velomsk];
            compsize[partval] ++;
          }
          else
            break;
        }
        fronnbr --;
        vertlist         = NULL;
        parttax[vertnum] = partval;                /* Move selected vertex in the current part */
        paloexc         += velobax[vertnum & velomsk];
        frloprt         -= velobax[vertnum & velomsk];
        fronload        -= velobax[vertnum & velomsk];

        for (edgenum = wgrafptr->s.verttax[vertnum]; edgenum < wgrafptr->s.vendtax[vertnum]; edgenum ++) {
          vertnum2 = wgrafptr->s.edgetax[edgenum]; /* for each neighbor vertex */

          if (parttax[vertnum2] != -1 && parttax[vertnum2] != partval) {
            parttax[vertnum2]  = -1;                /* Move the part in the separator    */
            fronnbr ++;
            fronload          += velobax[vertnum2 & velomsk];
            frloprt           += velobax[vertnum2 & velomsk];
            compload[partval] += velobax[vertnum2 & velomsk];
            compsize[partval] ++;
            vexxtax[vertnum2].partlvl = partval;    /* Label the vertex for this part    */
            vexxtax[vertnum2].prev    = vertlist;   /* Add the vertex in the list        */
            vertlist                  = vexxtax + vertnum2;

          } else if (parttax[vertnum2] == -1) {     /* If the vertex is in the separator */
            if (vexxtax[vertnum2].partlvl == partval) {
              if (vexxtax[vertnum2].isinstack != 1) {   /* If the vertex is in the stack */
              }
            }
            else {                  /* If the vertex is not labeled for the current part */
              frloprt           += velobax[vertnum2 & velomsk];
              compload[partval] += velobax[vertnum2 & velomsk];
              compsize[partval] ++;
              vexxtax[vertnum2].partlvl   = partval;            /* Label it for the part */
              vexxtax[vertnum2].isinstack = 1;
            }
          }
        }

        while (vertlist != NULL) {                           /* For eack linked vertices */
          vertnum2 = vertlist - vexxtax;
          stack[head ++ % wgrafptr->s.vertnbr] = vertnum2;
          vexxtax[vertnum2].isinstack = 0;
          vertlist = vertlist->prev;
        }

      } while ((paloexc + frloprt / 2) * (wgrafptr->partnbr - partval + 1) <= palooth); /* While the part is not big enought */
      palooth -= (paloexc + frloprt / 2);
    }

    compload[0] =
    compsize[0] = 0;
    for (vertnum = wgrafptr->s.baseval; vertnum < wgrafptr->s.vertnnd; vertnum ++) { /* Recompute load and size of part 0 */
      if (parttax[vertnum] == 0) {
        compload[0] += velobax[vertnum & velomsk];
        compsize[0] ++;
      } else if (parttax[vertnum] == -1) {
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

    if (frlobst > fronload) {                 /* If the pass frontier load is better than the better one */
      wgrafptr->fronnbr  = fronnbr;
      wgrafptr->fronload = fronload;
      memCpy (wgrafptr->compload, compload, sizeof (Gnum) * wgrafptr->partnbr);
      memCpy (wgrafptr->compsize, compsize, sizeof (Gnum) * wgrafptr->partnbr);
      memCpy (wgrafptr->parttax + wgrafptr->s.baseval, parttax + wgrafptr->s.baseval, sizeof (Gnum) * wgrafptr->s.vertnbr);
      for (vertnum = wgrafptr->s.baseval, fronnum = 0; vertnum < wgrafptr->s.vertnnd; vertnum ++) { /* Recompute frontab */
        if (parttax[vertnum] == -1)
          wgrafptr->frontab[fronnum ++] = vertnum;
      }
      frlobst = fronload;
    }
  }
  for (partval = 0; partval < wgrafptr->partnbr; partval ++) /* for each part */
    printf("\033[0;33mcompload[" GNUMSTRING "] " GNUMSTRING " " GNUMSTRING "\033[0m\n", partval, wgrafptr->compload[partval], wgrafptr->compsize[partval]);
  printf("\033[0;33mfronload " GNUMSTRING " " GNUMSTRING "\033[0m\n", wgrafptr->fronload, wgrafptr->fronnbr);
  memFree(vexxtax + wgrafptr->s.baseval);          /* Free work arrays */

  return (0);
}
