/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : separate_th.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module thins a vertex separator.   **/
/**                                                        **/
/**   DATES      : # Version 3.3  : from : 17 oct 1998     **/
/**                                 to     17 oct 1998     **/
/**                # Version 4.0  : from : 12 dec 2001     **/
/**                                 to     06 jan 2002     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_TH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_th.h"

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
vgraphSeparateTh (
Vgraph * const              grafptr)
{
  Gnum                fronnbr;                    /* Current number of frontier vertices */
  Gnum                fronnum;                    /* Number of current frontier vertex   */
  Gnum                commcut[3];                 /* Cut count array ([3] for halo)      */

  fronnbr = grafptr->fronnbr;                     /* Get current number of frontier vertices */
  for (fronnum = 0; fronnum < fronnbr; ) {
    Gnum                vertnum;
    Gnum                edgenum;

    vertnum    = grafptr->frontab[fronnum];       /* Get vertex number */
    commcut[0] =
    commcut[1] =
    commcut[2] = 0;
    for (edgenum = grafptr->s.verttax[vertnum]; edgenum < grafptr->s.vendtax[vertnum]; edgenum ++)
      commcut[grafptr->parttax[grafptr->s.edgetax[edgenum]]] ++;

    if (commcut[0] == 0) {
      grafptr->parttax[vertnum] = 1;
      grafptr->compload[1] += (grafptr->s.velotax == NULL) ? 1 : grafptr->s.velotax[vertnum];
      grafptr->compsize[1] ++;
      grafptr->frontab[fronnum] = grafptr->frontab[-- fronnbr]; /* Replace frontier vertex by another */
    }
    else if (commcut[1] == 0) {
      grafptr->parttax[vertnum] = 0;
      grafptr->compload[0] += (grafptr->s.velotax == NULL) ? 1 : grafptr->s.velotax[vertnum];
      grafptr->compsize[0] ++;
      grafptr->frontab[fronnum] = grafptr->frontab[-- fronnbr]; /* Replace frontier vertex by another */
    }
    else
      fronnum ++;                                 /* Keep vertex in separator */
  }
  grafptr->fronnbr     = fronnbr;                 /* Set new frontier parameters */
  grafptr->compload[2] = grafptr->s.velosum - (grafptr->compload[0] + grafptr->compload[1]);
  grafptr->comploaddlt = grafptr->compload[0] - grafptr->compload[1];

  return (0);
}
