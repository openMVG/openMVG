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
/**   NAME       : bdgraph_store.c                         **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the save data      **/
/**                structure handling routines for dis-    **/
/**                tributed bipartition graphs.            **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 10 sep 2007     **/
/**                                 to     22 oct 2008     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to     11 sep 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_STORE

#include "module.h"
#include "common.h"
#include "arch.h"
#include "dgraph.h"
#include "bdgraph.h"

/**********************************/
/*                                */
/* Store graph handling routines. */
/*                                */
/**********************************/

/* This routine builds a save structure
** for the given active graph.
** It returns:
** - 0   : if allocation succeeded.
** - !0  : on error.
*/

int
bdgraphStoreInit (
const Bdgraph * restrict const  grafptr,
BdgraphStore * restrict const   storptr)
{
  Gnum                savsize;

  savsize = grafptr->s.vertlocnbr * (sizeof (GraphPart) + sizeof (Gnum)); /* Compute size for frontier and part arrays */

  if ((storptr->datatab = (byte *) memAlloc (savsize)) == NULL) { /* Allocate save structure */
    errorPrint ("bdgraphStoreInit: out of memory");
    return     (1);
  }

  return (0);
}

/* This routine frees a save structure.
** It returns:
** - VOID  : in all cases.
*/

void
bdgraphStoreExit (
BdgraphStore * const        storptr)
{
  memFree (storptr->datatab);
#ifdef SCOTCH_DEBUG_BDGRAPH2
  storptr->datatab = NULL;
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
}

/* This routine saves partition data from the
** given active graph to the given save structure.
** It returns:
** - VOID  : in all cases.
*/

void
bdgraphStoreSave (
const Bdgraph * const       grafptr,
BdgraphStore * const        storptr)
{
  byte *              partloctab;                 /* Pointer to part data save area     */
  byte *              fronloctab;                 /* Pointer to frontier data save area */

  storptr->fronlocnbr      = grafptr->fronlocnbr; /* Save partition parameters */
  storptr->fronglbnbr      = grafptr->fronglbnbr;
  storptr->complocload0    = grafptr->complocload0;
  storptr->compglbload0    = grafptr->compglbload0;
  storptr->compglbload0dlt = grafptr->compglbload0dlt;
  storptr->complocsize0    = grafptr->complocsize0;
  storptr->compglbsize0    = grafptr->compglbsize0;
  storptr->commglbload     = grafptr->commglbload;
  storptr->commglbgainextn = grafptr->commglbgainextn;

  fronloctab = storptr->datatab;                  /* Compute data offsets within save structure */
  partloctab = fronloctab + grafptr->fronlocnbr * sizeof (Gnum);

  if (grafptr->fronloctab != NULL)                /* If frontier array allocated */
    memCpy (fronloctab, grafptr->fronloctab, grafptr->fronlocnbr * sizeof (Gnum));
#ifdef SCOTCH_DEBUG_BDGRAPH2
  else if (grafptr->fronglbnbr != 0)
    errorPrint ("bdgraphStoreSave: inconsistent graph data (1)");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
  if (grafptr->partgsttax != NULL)
    memCpy (partloctab, grafptr->partgsttax + grafptr->s.baseval, grafptr->s.vertlocnbr * sizeof (GraphPart));
  else {
#ifdef SCOTCH_DEBUG_BDGRAPH2
    if (grafptr->compglbload0 != grafptr->s.veloglbsum)
      errorPrint ("bdgraphStoreSave: inconsistent graph data (2)");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    memSet (partloctab, 0, grafptr->s.vertlocnbr * sizeof (GraphPart)); /* In case part array is allocated before update */
  }
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
bdgraphStoreUpdt (
Bdgraph * const             grafptr,
const BdgraphStore * const  storptr)
{
  byte *              fronloctab;                 /* Pointer to frontier data save area */
  byte *              partloctab;                 /* Pointer to part data save area     */

  grafptr->fronlocnbr      = storptr->fronlocnbr; /* Save partition parameters */
  grafptr->fronglbnbr      = storptr->fronglbnbr;
  grafptr->complocload0    = storptr->complocload0;
  grafptr->compglbload0    = storptr->compglbload0;
  grafptr->compglbload0dlt = storptr->compglbload0dlt;
  grafptr->complocsize0    = storptr->complocsize0;
  grafptr->compglbsize0    = storptr->compglbsize0;
  grafptr->commglbload     = storptr->commglbload;
  grafptr->commglbgainextn = storptr->commglbgainextn;
  grafptr->bbalglbval      = (double) ((grafptr->compglbload0dlt < 0) ? (- grafptr->compglbload0dlt) : grafptr->compglbload0dlt) / (double) grafptr->compglbload0avg;

  fronloctab = storptr->datatab;                  /* Compute data offsets within save structure */
  partloctab = fronloctab + grafptr->fronlocnbr * sizeof (Gnum);

  if (grafptr->fronloctab != NULL)
    memCpy (grafptr->fronloctab, fronloctab, grafptr->fronlocnbr * sizeof (Gnum));
#ifdef SCOTCH_DEBUG_BDGRAPH2
  else if (grafptr->fronglbnbr != 0)
    errorPrint ("bdgraphStoreUpdt: inconsistent graph data (1)");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  if (grafptr->partgsttax != NULL)
    memCpy (grafptr->partgsttax + grafptr->s.baseval, partloctab, grafptr->s.vertlocnbr * sizeof (GraphPart));
#ifdef SCOTCH_DEBUG_BDGRAPH2
  else if (grafptr->compglbload0 != grafptr->s.veloglbsum)
    errorPrint ("bdgraphStoreUpdt: inconsistent graph data (2)");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (grafptr) != 0)
    errorPrint ("bdgraphStoreUpdt: inconsistent graph data (3)");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
}
