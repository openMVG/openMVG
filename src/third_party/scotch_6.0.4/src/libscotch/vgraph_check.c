/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vgraph_check.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the separator      **/
/**                graph consistency checking routine.     **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 aug 1996     **/
/**                                 to     03 nov 1997     **/
/**                # Version 4.0  : from : 12 dec 2001     **/
/**                                 to     08 jan 2004     **/
/**                # Version 5.0  : from : 16 sep 2006     **/
/**                                 to   : 16 sep 2006     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 09 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"

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
vgraphCheck (
const Vgraph * const        grafptr)
{
  Gnum                vertnum;                    /* Number of current vertex  */
  Gnum                fronnum;                    /* Number of frontier vertex */
  Gnum                compload[3];
  Gnum                compsize[3];
  Gnum                commcut[3];

  if (grafptr->comploaddlt != (grafptr->compload[0] - grafptr->compload[1])) {
    errorPrint ("vgraphCheck: invalid balance");
    return     (1);
  }

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    if (grafptr->parttax[vertnum] > 2) {
      errorPrint ("vgraphCheck: invalid part array");
      return     (1);
    }
  }

  if ((grafptr->fronnbr < 0) ||
      (grafptr->fronnbr > grafptr->s.vertnbr)) {
    errorPrint ("vgraphCheck: invalid number of frontier vertices");
    return     (1);
  }
  for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
    Gnum                vertnum;

    vertnum = grafptr->frontab[fronnum];
    if ((vertnum < grafptr->s.baseval) || (vertnum >= grafptr->s.vertnnd)) {
      errorPrint ("vgraphCheck: invalid vertex index in frontier array");
      return     (1);
    }
    if (grafptr->parttax[vertnum] != 2) {
      errorPrint ("vgraphCheck: invalid vertex in frontier array");
      return     (1);
    }
  }

  compload[0] =
  compload[1] =
  compload[2] = 0;
  compsize[0] =
  compsize[1] =
  compsize[2] = 0;
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    int                 partnum;                  /* Part of current vertex */
    Gnum                edgenum;                  /* Number of current edge */

    partnum = (int) grafptr->parttax[vertnum];

    compload[partnum] += (grafptr->s.velotax == NULL) ? 1 : grafptr->s.velotax[vertnum];
    compsize[partnum] ++;

    commcut[0] =
    commcut[1] =
    commcut[2] = 0;
    if ((grafptr->s.verttax[vertnum] < grafptr->s.baseval) ||
        (grafptr->s.verttax[vertnum] > grafptr->s.vendtax[vertnum])) {
      errorPrint ("vgraphCheck: invalid graph structure (1)");
      return     (1);
    }
    for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
      Gnum                vertend;

      vertend = grafptr->s.edgetax[edgenum];
      if ((vertend <  grafptr->s.baseval) ||
          (vertend >= grafptr->s.vertnnd)) {
        errorPrint ("vgraphCheck: invalid graph structure (2)");
        return     (1);
      }
      commcut[grafptr->parttax[vertend]] ++;
    }

#ifdef SCOTCH_DEBUG_VGRAPH3
    if (partnum == 2) {
      if ((commcut[0] == 0) ||
          (commcut[1] == 0))
        errorPrintW ("vgraphCheck: no-use separator vertex%s (%ld)", /* Warning only */
                     ((grafptr->levlnum == 0) ? " at level 0" : ""),
                     (long) vertnum);
    }
    else {
#else
    if (partnum != 2) {
#endif /* SCOTCH_DEBUG_VGRAPH3 */
      if (commcut[1 - partnum] != 0) {
        errorPrint ("vgraphCheck: vertex should be in separator (%ld)", (long) vertnum);
        return     (1);
      }
    }
  }

  if ((grafptr->compload[0] != compload[0]) ||
      (grafptr->compload[1] != compload[1]) ||
      (grafptr->compload[2] != compload[2])) {
    errorPrint ("vgraphCheck: invalid part loads");
    return     (1);
  }
  if (grafptr->comploaddlt != (grafptr->compload[0] - grafptr->compload[1])) {
    errorPrint ("vgraphCheck: invalid balance");
    return     (1);
  }
  if ((grafptr->compsize[0] != compsize[0]) ||
      (grafptr->compsize[1] != compsize[1]) ||
      (grafptr->fronnbr     != compsize[2])) {
    errorPrint ("vgraphCheck: invalid part sizes");
    return     (1);
  }

  return (0);
}
