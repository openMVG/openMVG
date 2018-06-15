/* Copyright 2007,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_build_grid3d.c           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : These lines are the distributed source  **/
/**                graph building routines for 3D grid     **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 21 jul 2005     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 05 jun 2010     **/
/**                                 to   : 06 jun 2010     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/*+ This routine builds a distributed
*** 3D grid or torus graph structure.
*** It returns:
*** - 0   : if the creation succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphBuildGrid3D (
SCOTCH_Dgraph * const       grafptr,
const SCOTCH_Num            baseval,              /* Base value       */
const SCOTCH_Num            dimx,                 /* First dimension  */
const SCOTCH_Num            dimy,                 /* Second dimension */
const SCOTCH_Num            dimz,                 /* Third dimension  */
const SCOTCH_Num            incrval,              /* Increment value  */
const int                   flagval)              /* Flag value       */
{
  return (dgraphBuildGrid3D ((Dgraph *) grafptr, baseval, dimx, dimy, dimz, incrval, flagval));
}
