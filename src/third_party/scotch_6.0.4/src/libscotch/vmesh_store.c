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
/**   NAME       : vmesh_store.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the save data      **/
/**                structure handling routines for node    **/
/**                separation meshes.                      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 10 sep 2002     **/
/**                                 to   : 10 sep 2002     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_STORE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "vmesh.h"

/*********************************/
/*                               */
/* Store mesh handling routines. */
/*                               */
/*********************************/

/* This routine builds a save structure
** for the given node separation mesh.
** It returns:
** - 0   : if allocation succeeded.
** - !0  : on error.
*/

int
vmeshStoreInit (
const Vmesh * const         meshptr,
VmeshStore * const          storptr)
{
  Gnum                savsize;

  savsize = (meshptr->m.velmnbr + meshptr->m.vnodnbr) * (sizeof (GraphPart) + sizeof (Gnum)); /* Compute size for frontier and part arrays */

  if ((storptr->datatab = (byte *) memAlloc (savsize)) == NULL) { /* Allocate save structure */
    errorPrint ("vmeshStoreInit: out of memory");
    return     (1);
  }

  return (0);
}

/* This routine frees a mesh node
** separation save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vmeshStoreExit (
VmeshStore * const          storptr)
{
  memFree (storptr->datatab);
#ifdef SCOTCH_DEBUG_VMESH2
  storptr->datatab = NULL;
#endif /* SCOTCH_DEBUG_VMESH2 */
}

/* This routine saves partition data from the
** given node separation mesh to the given
** save structure.
** It returns:
** - VOID  : in all cases.
*/

void
vmeshStoreSave (
const Vmesh * const         meshptr,
VmeshStore * const          storptr)
{
  byte *              parttab;                    /* Pointer to part data save area     */
  byte *              frontab;                    /* Pointer to frontier data save area */

  storptr->ecmpsize[0] = meshptr->ecmpsize[0];    /* Save partition parameters */
  storptr->ecmpsize[1] = meshptr->ecmpsize[1];
  storptr->ncmpload[0] = meshptr->ncmpload[0];
  storptr->ncmpload[1] = meshptr->ncmpload[1];
  storptr->ncmpload[2] = meshptr->ncmpload[2];
  storptr->ncmploaddlt = meshptr->ncmploaddlt;
  storptr->ncmpsize[0] = meshptr->ncmpsize[0];
  storptr->ncmpsize[1] = meshptr->ncmpsize[1];
  storptr->fronnbr     = meshptr->fronnbr;

  frontab = storptr->datatab;                     /* Compute data offsets within save structure */
  parttab = frontab + meshptr->fronnbr * sizeof (Gnum);

  memCpy (frontab, meshptr->frontab, meshptr->fronnbr * sizeof (Gnum));
  memCpy (parttab, meshptr->parttax + meshptr->m.baseval, (meshptr->m.velmnbr + meshptr->m.vnodnbr) * sizeof (GraphPart));
}

/* This routine updates partition data of the
** given active graph, using the given save graph.
** It returns:
** - VOID  : in all cases.
*/

void
vmeshStoreUpdt (
Vmesh * const               meshptr,
const VmeshStore * const    storptr)
{
  byte *              frontab;                    /* Pointer to frontier data save area */
  byte *              parttab;                    /* Pointer to part data save area     */

  meshptr->ecmpsize[0] = storptr->ecmpsize[0];    /* Load partition parameters */
  meshptr->ecmpsize[1] = storptr->ecmpsize[1];
  meshptr->ncmpload[0] = storptr->ncmpload[0];
  meshptr->ncmpload[1] = storptr->ncmpload[1];
  meshptr->ncmpload[2] = storptr->ncmpload[2];
  meshptr->ncmploaddlt = storptr->ncmploaddlt;
  meshptr->ncmpsize[0] = storptr->ncmpsize[0];
  meshptr->ncmpsize[1] = storptr->ncmpsize[1];
  meshptr->fronnbr     = storptr->fronnbr;

  frontab = storptr->datatab;                     /* Compute data offsets within save structure */
  parttab = frontab + storptr->fronnbr * sizeof (Gnum);

  memCpy (meshptr->frontab, frontab, storptr->fronnbr * sizeof (Gnum));
  memCpy (meshptr->parttax + meshptr->m.baseval, parttab, (meshptr->m.velmnbr + meshptr->m.vnodnbr) * sizeof (GraphPart));
}
