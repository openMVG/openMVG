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
/**   NAME       : wgraph_check.c                          **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : This module check the graph consistency **/
/**                for the vertex overlapped graph partit- **/
/**                ioning.                                 **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 10 mar 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_CHECK

#include "module.h"
#include "common.h"
#include "graph.h"
#include "wgraph.h"

/*
**  The static variables.
*/

static const Gnum           wgraphcheckloadone = 1;

/*************************/
/*                       */
/* These routines handle */
/* separator graphs.     */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given separator graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
wgraphCheck (
const Wgraph * const        grafptr)
{
  Gnum                      vertnbr;              /* number of vertex                     */
  Gnum                      vertnum;              /* Number of current vertex             */
  const Gnum * restrict     velobax;              /* Data for handling of optional arrays */
  Gnum                      velomsk;
  Gnum                      edgenum;
  Gnum                      partnum;
  Gnum                      fronnum;              /* Current number of frontier vertex    */
  Gnum                      fronload;
  Gnum * restrict           compload;
  Gnum * restrict           compsize;
  Gnum * restrict           flagtab;

  vertnbr = grafptr->s.vertnbr;

  if (memAllocGroup ((void **) (void *)
                     &flagtab,  (size_t) (grafptr->partnbr * sizeof (Gnum)),
                     &compload, (size_t) (grafptr->partnbr * sizeof (Gnum)),
                     &compsize, (size_t) (grafptr->partnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("wgraphCheck: out of memory (1)");
    return     (1);
  }

  if (grafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &wgraphcheckloadone;                /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  fronnum  =
  fronload = 0;
  memSet (compload, 0, grafptr->partnbr * sizeof (Gnum)); /* Reset loads */
  memSet (compsize, 0, grafptr->partnbr * sizeof (Gnum));
  memSet (flagtab, ~0, grafptr->partnbr * sizeof (Gnum)); /* Reset flag array */

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    if ((grafptr->parttax[vertnum] >= grafptr->partnbr) ||
       	(grafptr->parttax[vertnum] < -1)) {
      errorPrint ("wgraphCheck: invalid part array");
      memFree    (flagtab);                       /* Free group leader */
      return     (1);
    }

    if (grafptr->parttax[vertnum] == -1) {
      fronnum ++;
      fronload += velobax[vertnum & velomsk];

      for (edgenum = grafptr->s.verttax[vertnum];
           edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
        Gnum               vertend;

        vertend = grafptr->s.edgetax[edgenum];
        if ((grafptr->parttax[vertend] != -1) &&
            (flagtab[grafptr->parttax[vertend]] != vertnum)) {
          compload[grafptr->parttax[vertend]] += velobax[vertnum & velomsk];
          compsize[grafptr->parttax[vertend]] ++;
          flagtab[grafptr->parttax[vertend]] = vertnum;
        }
      }
    }
    else {
      compload[grafptr->parttax[vertnum]] += velobax[vertnum & velomsk];
      compsize[grafptr->parttax[vertnum]] ++;
    }
  }

  for (partnum = 0; partnum < grafptr->partnbr; partnum ++) {
    if (grafptr->compsize[partnum] != compsize[partnum]) {
      errorPrint ("wgraphCheck: invalid part size %d %d %d", grafptr->compsize[partnum], compsize[partnum], partnum);
      memFree    (flagtab);
      return     (1);
    }
    if (grafptr->compload[partnum] != compload[partnum]) {
      errorPrintW ("wgraphCheck: invalid part load %d %d %d", grafptr->compload[partnum], compload[partnum], partnum);
      memFree     (flagtab);
      return      (1);
    }
  }

  if (grafptr->fronload != fronload) {
    errorPrint ("wgraphCheck: invalid frontier load %d %d", grafptr->fronload, fronload);
    memFree    (flagtab);
    return     (1);
  }
  if (grafptr->fronnbr != fronnum) {
    errorPrint ("wgraphCheck: invalid frontier size %d %d", grafptr->fronnbr, fronnum);
    memFree    (flagtab);
    return     (1);
  }

  for(fronnum = 0; fronnum < grafptr->fronnbr; fronnum++) {
    vertnum = grafptr->frontab[fronnum];
    if (grafptr->parttax[vertnum] != -1) {
      errorPrint ("wgraphCheck: invalid frontab");
      memFree    (flagtab);
      return     (1);
    }
  } 
  
  memFree (flagtab);

  return (0);
}
