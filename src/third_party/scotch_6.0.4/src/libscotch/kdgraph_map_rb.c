/* Copyright 2008,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kdgraph_map_rb.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm        **/ 
/**                in parallel.                            **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 apr 2008     **/
/**                                 to     01 jul 2008     **/
/**                # Version 6.0  : from : 03 oct 2012     **/
/**                                 to     10 oct 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KDGRAPH_MAP_RB

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "dgraph.h"
#include "dmapping.h"
#include "kdgraph.h"
#include "kdgraph_map_rb.h"
#include "kdgraph_map_rb_map.h"
#include "kdgraph_map_rb_part.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* These routines add mapping fragments to the distributed
** mapping structure.
*/

DmappingFrag *
kdgraphMapRbAdd2 (
const Gnum                vertnbr,
const Anum                domnnbr)
{
  DmappingFrag * restrict   fragptr;

  if ((fragptr = memAlloc (sizeof (DmappingFrag))) == NULL) {
    errorPrint  ("kdgraphMapRbAdd2: out of memory (1)");
    return      (NULL);
  }
  if (((fragptr->vnumtab = memAlloc (vertnbr * sizeof (Gnum)))    == NULL) ||
      ((fragptr->parttab = memAlloc (vertnbr * sizeof (Anum)))    == NULL) ||
      ((fragptr->domntab = memAlloc (domnnbr * sizeof (ArchDom))) == NULL)) {
    errorPrint ("kdgraphMapRbAdd2: out of memory (2)");
    if (fragptr->vnumtab != NULL) {
      if (fragptr->parttab != NULL)
        memFree (fragptr->parttab);
      memFree (fragptr->vnumtab);
    }
    return (NULL);
  }
  fragptr->vertnbr = vertnbr;
  fragptr->domnnbr = domnnbr;

  return (fragptr);
}

int
kdgraphMapRbAddBoth (
const Dgraph * restrict const     grafptr,
Dmapping * restrict const         mappptr,
const ArchDom * restrict const    domnptr,        /*+ Pointer to both subdomains   +*/
const GraphPart * restrict const  parttab)        /*+ Bipartition graph part array +*/
{
  DmappingFrag * restrict fragptr;
  Anum * restrict         fragparttab;
  Gnum                    vertlocnum;

  if ((fragptr = kdgraphMapRbAdd2 (grafptr->vertlocnbr, 2)) == NULL) /* Two domains */
    return (1);

  fragptr->domntab[0] = domnptr[0];
  fragptr->domntab[1] = domnptr[1];
  if (parttab == NULL)                            /* If bipartition part array not set */
    memSet (fragptr->parttab, 0, grafptr->vertlocnbr * sizeof (Anum));
  else {
    fragparttab = fragptr->parttab;
    for (vertlocnum = 0; vertlocnum < grafptr->vertlocnbr; vertlocnum ++)
      fragparttab[vertlocnum] = (Anum) parttab[vertlocnum];
  }

  if (grafptr->vnumloctax != NULL)
    memCpy (fragptr->vnumtab, grafptr->vnumloctax + grafptr->baseval, fragptr->vertnbr * sizeof (Gnum));
  else {
    Gnum * restrict     fragvnumtab;
    Gnum                vertlocadj;
    Gnum                vertlocnum;

    fragvnumtab = fragptr->vnumtab;
    for (vertlocnum = 0, vertlocadj = grafptr->procvrttab[grafptr->proclocnum]; vertlocnum < grafptr->vertlocnbr; vertlocnum ++)
      fragvnumtab[vertlocnum] = vertlocadj + vertlocnum;
  }

  dmapAdd (mappptr, fragptr);

  return (0);
}

int
kdgraphMapRbAddOne (
const Dgraph * restrict const   grafptr,
Dmapping * restrict const       mappptr,
const ArchDom * restrict const  domnptr)
{
  DmappingFrag * restrict   fragptr;

  if ((fragptr = kdgraphMapRbAdd2 (grafptr->vertlocnbr, 1)) == NULL) /* Only one domain */
    return (1);

  fragptr->domntab[0] = *domnptr;                 /* Only one domain for this mapping fragment */
  memSet (fragptr->parttab, 0, fragptr->vertnbr * sizeof (Anum)); /* All vertices mapped to it */
  if (grafptr->vnumloctax != NULL)
    memCpy (fragptr->vnumtab, grafptr->vnumloctax + grafptr->baseval, fragptr->vertnbr * sizeof (Gnum));
  else {
    Gnum * restrict     fragvnumtab;
    Gnum                vertlocadj;
    Gnum                vertlocnum;

    fragvnumtab = fragptr->vnumtab;
    for (vertlocnum = 0, vertlocadj = grafptr->procvrttab[grafptr->proclocnum]; vertlocnum < grafptr->vertlocnbr; vertlocnum ++)
      fragvnumtab[vertlocnum] = vertlocadj + vertlocnum;
  }

  dmapAdd (mappptr, fragptr);

  return (0);
}

int
kdgraphMapRbAddPart (
const Dgraph * restrict const   grafptr,
Dmapping * restrict const       mappptr,
const ArchDom * restrict const  domnptr,          /*+ Pointer to one subdomain +*/
const Gnum                      vertnbr,
const GraphPart * const         parttab,
const GraphPart                 partval)
{
  DmappingFrag * restrict fragptr;
  Gnum * restrict         fragvnumtab;
  Gnum                    vertlocnum;
  Gnum                    partlocnum;

  if ((fragptr = kdgraphMapRbAdd2 (vertnbr, 1)) == NULL) /* Only one domain and a limited number of vertices */
    return (1);

  fragptr->domntab[0] = *domnptr;                 /* Only one domain for this mapping fragment */
  memSet (fragptr->parttab, 0, fragptr->vertnbr * sizeof (Anum)); /* All vertices mapped to it */

  fragvnumtab = fragptr->vnumtab;

  if (grafptr->vnumloctax != NULL) {
    const Gnum * restrict vnumtab;

    for (vertlocnum = partlocnum = 0, vnumtab = grafptr->vnumloctax + grafptr->baseval; vertlocnum < grafptr->vertlocnbr; vertlocnum ++) {
      if (parttab[vertlocnum] == partval) {
#ifdef SCOTCH_DEBUG_KDMAP2
        if (partlocnum >= vertnbr) {
          errorPrint ("kdgraphMapRbAddPart: invalid parameters (1)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KDMAP2 */
        fragvnumtab[partlocnum ++] = vnumtab[vertlocnum];
      }
    }
  }
  else {
    Gnum              vertlocadj;

    for (vertlocnum = partlocnum = 0, vertlocadj = grafptr->procvrttab[grafptr->proclocnum];
         vertlocnum < grafptr->vertlocnbr; vertlocnum ++) {
      if (parttab[vertlocnum] == partval) {
#ifdef SCOTCH_DEBUG_KDMAP2
        if (partlocnum >= vertnbr) {
          errorPrint ("kdgraphMapRbAddPart: invalid parameters (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KDMAP2 */
        fragvnumtab[partlocnum ++] = vertlocadj + vertlocnum;
      }
    }
  }
#ifdef SCOTCH_DEBUG_KDMAP2
  if (partlocnum != vertnbr) {
    errorPrint ("kdgraphMapRbAddPart: invalid parameters (3)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KDMAP2 */

  dmapAdd (mappptr, fragptr);

  return (0);
}

/*
** This routine runs the parallel Dual
** Recursive Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kdgraphMapRb (
Kdgraph * restrict const                 grafptr,
Kdmapping * restrict const               mappptr,
const KdgraphMapRbParam * restrict const paraptr)
{
  if (grafptr->s.vertglbnbr == 0)                 /* If nothing to do, return immediately */
    return (0);

  return (archPart (&mappptr->mappptr->archdat)   /* If target architecture is some flavor of complete graph */
          ? kdgraphMapRbPart (grafptr, mappptr, paraptr)
          : kdgraphMapRbMap  (grafptr, mappptr, paraptr));
}
