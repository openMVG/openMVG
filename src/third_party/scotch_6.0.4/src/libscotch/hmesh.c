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
/**   NAME       : hmesh.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the halo source     **/
/**                mesh functions.                         **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 12 sep 2002     **/
/**                                 to     10 feb 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "hmesh.h"

/****************************************/
/*                                      */
/* These routines handle source meshes. */
/*                                      */
/****************************************/

/* This routine frees a source halo mesh structure.
** It returns:
** - VOID  : in all cases.
*/

void
hmeshExit (
Hmesh * const               meshptr)
{
  if ((meshptr->vehdtax != NULL) &&               /* Exit halo mesh data */
      (meshptr->vehdtax != (meshptr->m.vendtax + (meshptr->m.baseval - meshptr->m.velmbas))) &&
      ((meshptr->m.flagval & MESHVERTGROUP) == 0))
    memFree (meshptr->vehdtax + meshptr->m.velmbas);
  meshExit (&meshptr->m);                         /* Exit mesh data */

#ifdef SCOTCH_DEBUG_HMESH2
  memSet (meshptr, ~0, sizeof (Hmesh));           /* Purge halo mesh fields */
#endif /* SCOTCH_DEBUG_HMESH2 */
}

/* This routine sets the base of the given
** halo mesh to the given base value, and
** returns the old base value.
** It returns:
** - old base value : in all cases.
*/

Gnum
hmeshBase (
Hmesh * const               meshptr,
const Gnum                  baseval)
{
  Gnum                baseold;                    /* Old base value  */
  Gnum                baseadj;                    /* Base adjustment */
  Gnum                velmnum;

  if (meshptr->m.baseval == baseval)              /* If nothing to do */
    return (baseval);

  baseold = meshptr->m.baseval;                   /* Record old base value */
  baseadj = baseval - baseold;                    /* Compute adjustment    */

  meshBase (&meshptr->m, baseval);                /* Change base of mesh */

  for (velmnum = meshptr->m.velmbas; velmnum < meshptr->m.velmnnd; velmnum ++)
    meshptr->vehdtax[velmnum] += baseadj;         /* Change base of array */

  meshptr->vnohnnd += baseadj;
  meshptr->vehdtax -= baseadj;

  return (baseold);                               /* Return old base value */
}
