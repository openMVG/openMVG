/* Copyright 2007-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_part.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS partitioning      **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 sep 2006     **/
/**                                 to     07 jun 2007     **/
/**                # Version 5.1  : from : 06 jun 2009     **/
/**                                 to     30 jun 2010     **/
/**                # Version 6.0  : from : 23 dec 2011     **/
/**                                 to     13 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */

/************************************/
/*                                  */
/* These routines are the C API for */
/* MeTiS graph ordering routines.   */
/*                                  */
/************************************/

/* This routine is the interface between MeTiS
** and Scotch. It computes the partition of a
** weighted or unweighted graph.
** It returns:
** - 0   : if the partition could be computed.
** - !0  : on error.
*/

static
int
_SCOTCH_METIS_PartGraph2 (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
SCOTCH_Num * const          part,
SCOTCH_Num                  flagval,
double                      kbalval)
{
  SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat        stradat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertnbr;
  int                 o;

  SCOTCH_graphInit (&grafdat);

  baseval = *numflag;
  vertnbr = *n;

  o = 1;                                          /* Assume something will go wrong */
  if (SCOTCH_graphBuild (&grafdat, baseval, vertnbr, xadj, xadj + 1, vwgt, NULL,
                         xadj[vertnbr] - baseval, adjncy, adjwgt) == 0) {
    SCOTCH_stratInit          (&stradat);
    SCOTCH_stratGraphMapBuild (&stradat, flagval, *nparts, kbalval);
#ifdef SCOTCH_DEBUG_ALL
    if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
    o = SCOTCH_graphPart (&grafdat, *nparts, &stradat, part);
    SCOTCH_stratExit (&stradat);
  }
  SCOTCH_graphExit (&grafdat);

  if (o != 0)
    return (1);

  if (baseval != 0) {                             /* MeTiS part array is based, Scotch is not */
    SCOTCH_Num          vertnum;

    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      part[vertnum] += baseval;
  }

  return (0);
}

/*
**
*/

static
void
_SCOTCH_METIS_PartGraph (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part,
SCOTCH_Num                  flagval,
double                      kbalval)
{
  const SCOTCH_Num *          vwgt2;
  const SCOTCH_Num *          adjwgt2;
  const SCOTCH_Num * restrict parttax;
  const SCOTCH_Num * restrict verttax;
  const SCOTCH_Num * restrict edgetax;
  SCOTCH_Num                  vertnnd;
  SCOTCH_Num                  vertnum;
  SCOTCH_Num                  edgenum;
  SCOTCH_Num                  commcut;

  vwgt2   = (((*wgtflag & 2) != 0) ? vwgt   : NULL);
  adjwgt2 = (((*wgtflag & 1) != 0) ? adjwgt : NULL);

  if (_SCOTCH_METIS_PartGraph2 (n, xadj, adjncy, vwgt2, adjwgt2, numflag, nparts, part, flagval, kbalval) != 0) {
    *edgecut = -1;                                /* Indicate error */
    return;
  }

  parttax = part   - *numflag;
  verttax = xadj   - *numflag;
  edgetax = adjncy - *numflag;
  edgenum = *numflag;
  vertnum = *numflag;
  vertnnd = *n + vertnum;
  commcut = 0;

  if (adjwgt2 == NULL) {                          /* If graph does not have edge weights */
    for ( ; vertnum < vertnnd; vertnum ++) {
      SCOTCH_Num          edgennd;
      SCOTCH_Num          partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        if (parttax[edgetax[edgenum]] != partval)
          commcut ++;
      }
    }
  }
  else {                                          /* Graph has edge weights */
    const SCOTCH_Num * restrict edlotax;

    edlotax = adjwgt2 - *numflag;
    for ( ; vertnum < vertnnd; vertnum ++) {
      SCOTCH_Num          edgennd;
      SCOTCH_Num          partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        SCOTCH_Num          vertend;

        vertend = edgetax[edgenum];
        if (parttax[vertend] != partval)
          commcut += edlotax[edgenum];
      }
    }
  }
  *edgecut = commcut / 2;
}

/*
**
*/

void
METISNAMEU (METIS_PartGraphKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  _SCOTCH_METIS_PartGraph (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part,
                           SCOTCH_STRATDEFAULT, 0.01);
}

void
METISNAMEU (METIS_PartGraphRecursive) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  _SCOTCH_METIS_PartGraph (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part,
                           SCOTCH_STRATRECURSIVE, 0.01);
}

/* Scotch does not directly consider communication volume.
** Instead, wertex communication loads are added to the edge
** loads so as to emulate this behavior: heavily weighted
** edges, connected to heavily communicating vertices, will
** be less likely to be cut.
*/

void
METISNAMEU (METIS_PartGraphVKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          volume,
SCOTCH_Num * const          part)
{
  SCOTCH_Num                  baseval;
  const SCOTCH_Num *          vwgt2;
  const SCOTCH_Num *          vsize2;
  SCOTCH_Num                  vsizval;            /* Communication volume of current vertex */
  SCOTCH_Num                  vertnbr;
  SCOTCH_Num                  vertnum;
  SCOTCH_Num                  edgenum;
  const SCOTCH_Num * restrict edgetax;
  const SCOTCH_Num * restrict parttax;
  SCOTCH_Num * restrict       nghbtab;
  SCOTCH_Num                  commvol;

  vsize2  = ((*wgtflag & 1) != 0) ? vsize : NULL;
  vwgt2   = ((*wgtflag & 2) != 0) ? vwgt  : NULL;
  baseval = *numflag;
  vertnbr = *n;
  edgetax = adjncy - baseval;

  if (vsize2 == NULL) {                           /* If no communication load data provided */
    if (_SCOTCH_METIS_PartGraph2 (n, xadj, adjncy, vwgt2, NULL, numflag, nparts, part, SCOTCH_STRATDEFAULT, 0.01) != 0)
      return;
  }
  else {                                          /* Will have to turn communication volumes into edge loads */
    const SCOTCH_Num * restrict vsiztax;
    SCOTCH_Num                  edgenbr;
    SCOTCH_Num * restrict       edlotax;
    int                         o;

    edgenbr = xadj[vertnbr] - baseval;
    if ((edlotax = memAlloc (edgenbr * sizeof (SCOTCH_Num))) == NULL)
      return;
    edlotax -= baseval;                           /* Base access to edlotax */
    vsiztax  = vsize2 - baseval;

    for (vertnum = 0, edgenum = baseval;          /* Un-based scan of vertex array xadj */
         vertnum < vertnbr; vertnum ++) {
      SCOTCH_Num          vsizval;                /* Communication size of current vertex */
      SCOTCH_Num          edgennd;

      vsizval = vsize2[vertnum];
      for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
        SCOTCH_Num          vertend;              /* Based end vertex number                                     */

        vertend = edgetax[edgenum];
        edlotax[edgenum] = vsizval + vsiztax[vertend];
      }
    }

    o = _SCOTCH_METIS_PartGraph2 (n, xadj, adjncy, vwgt2, edlotax + baseval, numflag, nparts, part,
                                  SCOTCH_STRATDEFAULT, 0.01);

    memFree (edlotax + baseval);

    if (o != 0)
      return;
  }

  if ((nghbtab = memAlloc (*nparts * sizeof (SCOTCH_Num))) == NULL)
    return;
  memSet (nghbtab, ~0, *nparts * sizeof (SCOTCH_Num));

  parttax = part - baseval;
  vsizval = 1;                                      /* Assume no vertex communication sizes */
  for (vertnum = 0, edgenum = baseval, commvol = 0; /* Un-based scan of vertex array xadj   */
       vertnum < vertnbr; vertnum ++) {
    SCOTCH_Num          partval;
    SCOTCH_Num          edgennd;

    partval = part[vertnum];
    nghbtab[partval] = vertnum;                   /* Do not count local neighbors in communication volume */
    if (vsize2 != NULL)
      vsizval = vsize2[vertnum];

    for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
      SCOTCH_Num          vertend;                /* Based end vertex number                                   */
      SCOTCH_Num          partend;

      vertend = edgetax[edgenum];
      partend = parttax[vertend];
      if (nghbtab[partend] != vertnum) {          /* If first neighbor in this part */
        nghbtab[partend] = vertnum;               /* Set part as accounted for      */
        commvol += vsizval;
      }
    }
  }
  *volume = commvol;

  memFree (nghbtab);
}
