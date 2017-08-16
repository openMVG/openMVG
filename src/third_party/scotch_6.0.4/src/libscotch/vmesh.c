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
/**   NAME       : vmesh.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the separator      **/
/**                handling routines.                      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 06 feb 2003     **/
/**                                 to     05 mar 2003     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "vmesh.h"

/*************************/
/*                       */
/* These routines handle */
/* separator meshes.     */
/*                       */
/*************************/

/* This routine frees the contents
** of the given active mesh.
** It returns:
** - VOID  : in all cases.
*/

void
vmeshExit (
Vmesh * const               meshptr)
{
  if (meshptr->parttax != NULL)                   /* Free leader of group (parttab + frontab) */
    memFree (meshptr->parttax + meshptr->m.baseval);

  meshFree (&meshptr->m);                         /* Free source mesh */

#ifdef SCOTCH_DEBUG_VMESH2
  memSet (meshptr, ~0, sizeof (Vmesh));
#endif /* SCOTCH_DEBUG_VMESH2 */
}

/* This routine moves all of the mesh
** elements to the first part.
** It returns:
** - VOID  : in all cases.
*/

void
vmeshZero (
Vmesh * const               meshptr)
{
  memSet (meshptr->parttax + meshptr->m.baseval, 0, (meshptr->m.velmnbr + meshptr->m.vnodnbr) * sizeof (GraphPart));

  meshptr->ecmpsize[0] = meshptr->m.velmnbr;
  meshptr->ecmpsize[1] = 0;
  meshptr->ncmpload[0] = meshptr->m.vnlosum;
  meshptr->ncmpload[1] =
  meshptr->ncmpload[2] = 0;
  meshptr->ncmploaddlt = meshptr->m.vnlosum;
  meshptr->ncmpsize[0] = meshptr->m.vnodnbr;
  meshptr->ncmpsize[1] =
  meshptr->fronnbr     = 0;
}
