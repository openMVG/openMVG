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
/**   NAME       : mesh.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source mesh     **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 29 dec 2001     **/
/**                                 to     05 may 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"

/****************************************/
/*                                      */
/* These routines handle source meshes. */
/*                                      */
/****************************************/

/* This routine initializes a source mesh
** structure.
** It returns:
** - 0  : in all cases.
*/

int
meshInit (
Mesh * const                meshptr)
{
  memSet (meshptr, 0, sizeof (Mesh));             /* Initialize mesh fields      */
  meshptr->flagval = MESHFREETABS;                /* By default, free all arrays */

  return (0);
}

/* This routine frees a source mesh structure.
** It returns:
** - VOID  : in all cases.
*/

void
meshExit (
Mesh * const                meshptr)
{
  meshFree (meshptr);                             /* Exit mesh data */

#ifdef SCOTCH_DEBUG_MESH2
  memSet (meshptr, ~0, sizeof (Mesh));            /* Purge mesh fields */
#endif /* SCOTCH_DEBUG_MESH2 */
}

/* This routine frees the mesh data. Because
** vertex load arrays are either passed by
** the user, or grouped with other arrays,
** they are not considered for explicit
** freeing. This is also much simpler, as
** load arrays can be grouped or not.
** It returns:
** - VOID  : in all cases.
*/

void
meshFree (
Mesh * const                meshptr)
{
  if (((meshptr->flagval & MESHFREEEDGE) != 0) && /* If edgetab must be freed */
      (meshptr->edgetax != NULL))                 /* And if it exists          */
    memFree (meshptr->edgetax + meshptr->baseval); /* Free it                  */

  if ((meshptr->flagval & MESHFREEVEND) != 0) {   /* If vendtab must be freed                                    */
    if ((meshptr->vendtax != NULL) &&             /* If vendtab is distinct from verttab                         */
        (meshptr->vendtax != meshptr->verttax + 1) && /* (if vertex arrays grouped, vendtab not distinct anyway) */
        ((meshptr->flagval & MESHVERTGROUP) == 0))
      memFree (meshptr->vendtax + meshptr->baseval); /* Then free vendtab */
  }
  if ((meshptr->flagval & MESHFREEVERT) != 0) {   /* If verttab must be freed                             */
    if (meshptr->verttax != NULL)                 /* Free verttab anyway, as it is the array group leader */
      memFree (meshptr->verttax + meshptr->baseval);
  }
#ifdef SCOTCH_DEBUG_MESH2
  if ((meshptr->flagval & MESHFREEVNUM) != 0) {   /* If vnumtab must be freed         */
    if ((meshptr->vnumtax != NULL) &&             /* And is not in vertex array group */
        ((meshptr->flagval & MESHVERTGROUP) == 0))
      errorPrint ("meshFree: vnumtab should never be freed as its base may vary according to creation routines");
  }
#endif /* SCOTCH_DEBUG_MESH2 */
  if ((meshptr->flagval & MESHFREEOTHR) != 0) {   /* If other arrays must be freed */
    if (meshptr->vlbltax != NULL)
      memFree (meshptr->vlbltax + meshptr->baseval);
  }

#ifdef SCOTCH_DEBUG_MESH2
  memSet (meshptr, ~0, sizeof (Mesh));            /* Purge mesh fields */
#endif /* SCOTCH_DEBUG_MESH2 */
}

/* This routine sets the base of the given
** mesh to the given base value, and returns
** the old base value.
** It returns:
** - old base value : in all cases.
*/

Gnum
meshBase (
Mesh * const                meshptr,
const Gnum                  baseval)
{
  Gnum                baseold;                    /* Old base value  */
  Gnum                baseadj;                    /* Base adjustment */
  Gnum                vertnum;
  Gnum                edgenum;

  if (meshptr->baseval == baseval)                /* If nothing to do */
    return (baseval);

  baseold = meshptr->baseval;                     /* Record old base value */
  baseadj = baseval - baseold;                    /* Compute adjustment    */

  for (vertnum = meshptr->baseval; vertnum < (meshptr->velmnbr + meshptr->vnodnbr + meshptr->baseval); vertnum ++) {
    for (edgenum = meshptr->verttax[vertnum]; edgenum < meshptr->vendtax[vertnum]; edgenum ++)
      meshptr->edgetax[edgenum] += baseadj;
    meshptr->verttax[vertnum] += baseadj;
  }
  if (meshptr->vendtax != meshptr->verttax + 1) { /* If distinct vertex end array */
    for (vertnum = meshptr->baseval; vertnum < (meshptr->velmnbr + meshptr->vnodnbr + meshptr->baseval); vertnum ++)
      meshptr->vendtax[vertnum] += baseadj;
  }
  else                                            /* If same vertex end array (of size +1)                               */
    meshptr->verttax[meshptr->velmnbr + meshptr->vnodnbr + meshptr->baseval] += baseadj; /* Adjust last entry of verttab */

  meshptr->verttax -= baseadj;                    /* Adjust array accesses */
  meshptr->vendtax -= baseadj;
  meshptr->edgetax -= baseadj;

  if (meshptr->vnumtax != NULL)
    meshptr->vnumtax -= baseadj;
  if (meshptr->vlbltax != NULL)
    meshptr->vlbltax -= baseadj;

  meshptr->baseval  = baseval;                    /* Set new base value    */
  meshptr->velmbas += baseadj;                    /* Adjust mesh parameter */
  meshptr->velmnnd += baseadj;
  meshptr->vnodbas += baseadj;
  meshptr->vnodnnd += baseadj;

  return (baseold);                               /* Return old base value */
}
