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
/**   NAME       : vdgraph_store.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the save data      **/
/**                structure handling routines for         **/
/**                distributed separation graphs.          **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 08 mar 2006     **/
/**                                 to   : 01 mar 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_STORE

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "vdgraph.h"

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
vdgraphStoreInit (
const Vdgraph * restrict const  grafptr,
VdgraphStore * restrict const   storptr)
{
  Gnum                savsize;

  savsize = grafptr->s.vertlocnbr * (sizeof (GraphPart) + sizeof (Gnum)); /* Compute size for frontier and part arrays */

  if ((storptr->datatab = (byte *) memAlloc (savsize)) == NULL) { /* Allocate save structure */
    errorPrint ("vdgraphStoreInit: out of memory");
    return     (1);
  }

  return (0);
}

/* This routine frees a save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vdgraphStoreExit (
VdgraphStore * const        storptr)
{
  memFree (storptr->datatab);
#ifdef SCOTCH_DEBUG_VDGRAPH2
  storptr->datatab = NULL;
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
}

/* This routine saves partition data from the
** given active graph to the given save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vdgraphStoreSave (
const Vdgraph * const       grafptr,
VdgraphStore * const        storptr)
{
  byte *              partloctab;                 /* Pointer to part data save area     */
  byte *              fronloctab;                 /* Pointer to frontier data save area */

  storptr->fronglbnbr     = grafptr->compglbsize[2]; /* Save partition parameters */
  storptr->compglbloaddlt = grafptr->compglbloaddlt;
  storptr->compglbload[0] = grafptr->compglbload[0];
  storptr->compglbload[1] = grafptr->compglbload[1];
  storptr->compglbsize0   = grafptr->compglbsize[0];
  storptr->complocsize0   = grafptr->complocsize[0];
  storptr->fronlocnbr     = grafptr->complocsize[2];

  fronloctab = storptr->datatab;                  /* Compute data offsets within save structure */
  partloctab = fronloctab + grafptr->complocsize[2] * sizeof (Gnum);

  memCpy (fronloctab, grafptr->fronloctab, grafptr->complocsize[2] * sizeof (Gnum));
  memCpy (partloctab, grafptr->partgsttax + grafptr->s.baseval, grafptr->s.vertlocnbr * sizeof (GraphPart));
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
vdgraphStoreUpdt (
Vdgraph * const             grafptr,
const VdgraphStore * const  storptr)
{
  byte *              fronloctab;                 /* Pointer to frontier data save area */
  byte *              partloctab;                 /* Pointer to part data save area     */

  grafptr->compglbload[0] = storptr->compglbload[0]; /* Load partition parameters */
  grafptr->compglbload[1] = storptr->compglbload[1];
  grafptr->compglbload[2] = grafptr->s.veloglbsum - (storptr->compglbload[0] + storptr->compglbload[1]);
  grafptr->compglbloaddlt = storptr->compglbloaddlt;
  grafptr->compglbsize[0] = storptr->compglbsize0;
  grafptr->compglbsize[1] = grafptr->s.vertglbnbr - (storptr->compglbsize0 + storptr->fronglbnbr);
  grafptr->compglbsize[2] = storptr->fronglbnbr;
  grafptr->complocsize[0] = storptr->complocsize0;
  grafptr->complocsize[1] = grafptr->s.vertlocnbr - (storptr->complocsize0 + storptr->fronlocnbr);
  grafptr->complocsize[2] = storptr->fronlocnbr;

  fronloctab = storptr->datatab;                  /* Compute data offsets within save structure */
  partloctab = fronloctab + grafptr->complocsize[2] * sizeof (Gnum);

  memCpy (grafptr->fronloctab, fronloctab, grafptr->complocsize[2] * sizeof (Gnum));
  memCpy (grafptr->partgsttax + grafptr->s.baseval, partloctab, grafptr->s.vertlocnbr * sizeof (GraphPart));

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (grafptr) != 0)
    errorPrint ("vdgraphStoreUpdt: inconsistent graph data");
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
}
