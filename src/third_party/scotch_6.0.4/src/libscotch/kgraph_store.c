/* Copyright 2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kgraph_store.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the save data      **/
/**                structure handling routines for k-par-  **/
/**                tition graphs.                          **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 16 aug 2011     **/
/**                                 to     16 aug 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_STORE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "kgraph.h"

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
kgraphStoreInit (
const Kgraph * const        grafptr,
KgraphStore * const         storptr)
{
  ArchDom             domnfrst;                   /* Largest domain in the architecture */

  archDomFrst (&grafptr->a, &domnfrst);           /* Get architecture domain               */
  storptr->partnbr = (Gnum) archDomSize (&grafptr->a, &domnfrst); /* Get architecture size */

  if (memAllocGroup ((void **) (void *)           /* Allocate save structure */
                    &storptr->parttab, (size_t) (grafptr->s.vertnbr * sizeof (Anum)),
                    &storptr->domntab, (size_t) (grafptr->m.domnmax * sizeof (ArchDom)),
                    &storptr->frontab, (size_t) (grafptr->s.vertnbr * sizeof (Gnum)),
                    &storptr->comploadavg, (size_t) (storptr->partnbr * sizeof (Gnum)),
                    &storptr->comploaddlt, (size_t) (storptr->partnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint   ("kgraphStoreInit out of memory (1)");
    return       (1);
  }

  return (0);
}

/* This routine frees a save structure.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphStoreExit (
KgraphStore * const         storptr)
{
  memFree (storptr->parttab);                     /* Free group leader */
#ifdef SCOTCH_DEBUG_KGRAPH2
  storptr->parttab     = NULL;
  storptr->domntab     = NULL;
  storptr->frontab     = NULL;
  storptr->comploadavg = NULL;
  storptr->comploaddlt = NULL;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
}

/* This routine saves partition data from the
** given active graph to the given save structure.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphStoreSave (
const Kgraph * const        grafptr,
KgraphStore * const         storptr)
{
  storptr->mflaval  = grafptr->commload;
  storptr->domnnbr  = grafptr->m.domnnbr;
  storptr->fronnbr  = grafptr->fronnbr;
  storptr->kbalval  = grafptr->kbalval;
  storptr->commload = grafptr->commload;

  memCpy (storptr->parttab,     grafptr->m.parttax + grafptr->s.baseval, grafptr->s.vertnbr * sizeof (Anum));
  memCpy (storptr->domntab,     grafptr->m.domntab,                      grafptr->m.domnnbr * sizeof (ArchDom));
  memCpy (storptr->frontab,     grafptr->frontab,                        storptr->fronnbr * sizeof (Gnum));
  memCpy (storptr->comploadavg, grafptr->comploadavg,                    storptr->partnbr * sizeof (Gnum));
  memCpy (storptr->comploaddlt, grafptr->comploaddlt,                    storptr->partnbr * sizeof (Gnum));
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
kgraphStoreUpdt (
Kgraph * const              grafptr,
const KgraphStore * const   storptr)
{
  grafptr->commload  = storptr->mflaval;
  grafptr->m.domnnbr = storptr->domnnbr;
  grafptr->fronnbr   = storptr->fronnbr;
  grafptr->kbalval   = storptr->kbalval;
  grafptr->commload  = storptr->commload;

  memCpy (grafptr->m.parttax + grafptr->s.baseval, storptr->parttab,     grafptr->s.vertnbr * sizeof (Anum));
  memCpy (grafptr->m.domntab,                      storptr->domntab,     grafptr->m.domnnbr * sizeof (ArchDom));
  memCpy (grafptr->frontab,                        storptr->frontab,     storptr->fronnbr * sizeof (Gnum));
  memCpy (grafptr->comploadavg,                    storptr->comploadavg, storptr->partnbr * sizeof (Gnum));
  memCpy (grafptr->comploaddlt,                    storptr->comploaddlt, storptr->partnbr * sizeof (Gnum));

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0)
    errorPrint ("kgraphStoreUpdt: inconsistent graph data");
#endif /* SCOTCH_DEBUG_KGRAPH2 */
}
