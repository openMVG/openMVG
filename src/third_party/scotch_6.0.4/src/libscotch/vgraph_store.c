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
/**   NAME       : vgraph_store.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the save data      **/
/**                structure handling routines for separa- **/
/**                tion graphs.                            **/
/**                                                        **/
/**   DATES      : # Version 3.3  : from : 17 oct 1998     **/
/**                                 to     17 oct 1998     **/
/**                # Version 3.4  : from : 11 dec 2001     **/
/**                                 to   : 11 dec 2001     **/
/**                # Version 4.0  : from : 01 jan 2002     **/
/**                                 to   : 06 jan 2002     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_STORE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "vgraph.h"

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
vgraphStoreInit (
const Vgraph * restrict const grafptr,
VgraphStore * restrict const  storptr)
{
  Gnum                savsize;

  savsize = grafptr->s.vertnbr * (sizeof (GraphPart) + sizeof (Gnum)); /* Compute size for frontier and part arrays */

  if ((storptr->datatab = (byte *) memAlloc (savsize)) == NULL) { /* Allocate save structure */
    errorPrint ("vgraphStoreInit: out of memory");
    return     (1);
  }

  return (0);
}

/* This routine frees a save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vgraphStoreExit (
VgraphStore * const         storptr)
{
  memFree (storptr->datatab);
#ifdef SCOTCH_DEBUG_VGRAPH2
  storptr->datatab = NULL;
#endif /* SCOTCH_DEBUG_VGRAPH2 */
}

/* This routine saves partition data from the
** given active graph to the given save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vgraphStoreSave (
const Vgraph * const        grafptr,
VgraphStore * const         storptr)
{
  byte *              parttab;                    /* Pointer to part data save area     */
  byte *              frontab;                    /* Pointer to frontier data save area */

  storptr->fronnbr     = grafptr->fronnbr;        /* Save partition parameters */
  storptr->comploaddlt = grafptr->comploaddlt;
  storptr->compload[0] = grafptr->compload[0];
  storptr->compload[1] = grafptr->compload[1];
  storptr->compsize0   = grafptr->compsize[0];

  frontab = storptr->datatab;                     /* Compute data offsets within save structure */
  parttab = frontab + grafptr->fronnbr * sizeof (Gnum);

  memCpy (frontab, grafptr->frontab, grafptr->fronnbr * sizeof (Gnum));
  memCpy (parttab, grafptr->parttax + grafptr->s.baseval, grafptr->s.vertnbr * sizeof (GraphPart));
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
vgraphStoreUpdt (
Vgraph * const              grafptr,
const VgraphStore * const   storptr)
{
  byte *              frontab;                    /* Pointer to frontier data save area */
  byte *              parttab;                    /* Pointer to part data save area     */

  grafptr->compload[0] = storptr->compload[0];    /* Load partition parameters */
  grafptr->compload[1] = storptr->compload[1];
  grafptr->compload[2] = grafptr->s.velosum - (storptr->compload[0] + storptr->compload[1]);
  grafptr->comploaddlt = storptr->comploaddlt;
  grafptr->compsize[0] = storptr->compsize0;
  grafptr->compsize[1] = grafptr->s.vertnbr - (storptr->compsize0 + storptr->fronnbr);
  grafptr->fronnbr     = storptr->fronnbr;

  frontab = storptr->datatab;                     /* Compute data offsets within save structure */
  parttab = frontab + grafptr->fronnbr * sizeof (Gnum);

  memCpy (grafptr->frontab, frontab, grafptr->fronnbr * sizeof (Gnum));
  memCpy (grafptr->parttax + grafptr->s.baseval, parttab, grafptr->s.vertnbr * sizeof (GraphPart));

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (grafptr) != 0)
    errorPrint ("vgraphStoreUpdt: inconsistent graph data");
#endif /* SCOTCH_DEBUG_VGRAPH2 */
}
