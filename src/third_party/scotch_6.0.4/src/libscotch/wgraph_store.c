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
/**   NAME       : wgraph_store.c                          **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : This module contains the data store re- **/
/**                lated rountines for the vertex overlap- **/
/**                ped graph partitioning.                 **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2010     **/
/**                                 to   : 30 may 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_STORE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "wgraph.h"

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
wgraphStoreInit (
const Wgraph * const        grafptr,
WgraphStore * const         storptr)
{
  Gnum                savsize;
  
  savsize = 2 * grafptr->partnbr * sizeof (Gnum) + /* Compute size for frontier, part arrays, communication load and size */
            grafptr->s.vertnbr * (sizeof (Gnum) + sizeof (Anum));

  if ((storptr->datatab = (byte *) memAlloc (savsize)) == NULL) { /* Allocate save structure */
    errorPrint ("wgraphStoreInit: out of memory");
    return     (1);
  }

  return (0);
}

/* This routine frees a save structure.
** It returns:
** - VOID  : in all cases.
*/

void
wgraphStoreExit (
WgraphStore * const         storptr)
{
  memFree (storptr->datatab);
#ifdef SCOTCH_DEBUG_WGRAPH2
  storptr->datatab = NULL;
#endif /* SCOTCH_DEBUG_WGRAPH2 */
}

/* This routine saves partition data from the
** given active graph to the given save structure.
** It returns:
** - VOID  : in all cases.
*/

void
wgraphStoreSave (
const Wgraph * const        grafptr,
WgraphStore * const         storptr)
{
  byte *              compload;                   /* Pointer to part load data save area */
  byte *              compsize;                   /* Pointer to part size data save area */
  byte *              frontab;                    /* Pointer to frontier data save area  */
  byte *              parttab;                    /* Pointer to partition data save area */

  storptr->partnbr  = grafptr->partnbr;           /* Save partition parameters */
  storptr->fronnbr  = grafptr->fronnbr;
  storptr->fronload = grafptr->fronload;

  compload = storptr->datatab;                    /* Compute data offsets within save structure */
  compsize = compload + grafptr->partnbr * sizeof (Gnum);
  frontab  = compsize + grafptr->partnbr * sizeof (Gnum);
  parttab  = frontab  + grafptr->fronnbr * sizeof (Gnum);

  memCpy (compload, grafptr->compload, grafptr->partnbr * sizeof (Gnum));
  memCpy (compsize, grafptr->compsize, grafptr->partnbr * sizeof (Gnum));
  memCpy (frontab,  grafptr->frontab,  grafptr->fronnbr * sizeof (Gnum));
  memCpy (parttab,  grafptr->parttax + grafptr->s.baseval, grafptr->s.vertnbr * sizeof (Anum));
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
wgraphStoreUpdt (
Wgraph * const              grafptr,
const WgraphStore * const   storptr)
{
  byte *              compload;                   /* Pointer to part load data save area */
  byte *              compsize;                   /* Pointer to part size data save area */
  byte *              frontab;                    /* Pointer to frontier data save area  */
  byte *              parttab;                    /* Pointer to partition data save area */

  grafptr->partnbr  = storptr->partnbr;           /* Load partition parameters */
  grafptr->fronnbr  = storptr->fronnbr;
  grafptr->fronload = storptr->fronload;

  compload = storptr->datatab;                    /* Compute data offsets within save structure */
  compsize = compload + grafptr->partnbr * sizeof (Gnum);
  frontab  = compsize + grafptr->partnbr * sizeof (Gnum);
  parttab  = frontab  + grafptr->fronnbr * sizeof (Gnum);

  memCpy (grafptr->compload, compload, grafptr->partnbr * sizeof (Gnum));
  memCpy (grafptr->compsize, compsize, grafptr->partnbr * sizeof (Gnum));
  memCpy (grafptr->frontab,  frontab,  grafptr->fronnbr * sizeof (Gnum));
  memCpy (grafptr->parttax + grafptr->s.baseval, parttab, grafptr->s.vertnbr * sizeof (Anum));

#ifdef SCOTCH_DEBUG_WGRAPH2
  if (wgraphCheck (grafptr) != 0)
    errorPrint ("wgraphStoreUpdt: inconsistent graph data");
#endif /* SCOTCH_DEBUG_WGRAPH2 */
}
