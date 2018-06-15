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
/**   NAME       : library_mesh_io_habo.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the Harwell- **/
/**                Boeing geometry and mesh handling       **/
/**                routines of the libSCOTCH library.      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 jan 2004     **/
/**                                 to     19 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "mesh.h"
#include "scotch.h"

/*************************************/
/*                                   */
/* These routines are the C API for  */
/* the Harwell-Boeing mesh and       */
/* geometry handling routines.       */
/*                                   */
/*************************************/

/*+ This routine loads the given opaque geom
*** structure with the data of the given stream.
*** - 0   : if loading succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_meshGeomLoadHabo (
SCOTCH_Mesh * restrict const  meshptr,
SCOTCH_Geom * restrict const  geomptr,
FILE * const                  filesrcptr,
FILE * const                  filegeoptr,
const char * const            dataptr)            /* No use */
{
  return (meshGeomLoadHabo ((Mesh *) meshptr, (Geom *) geomptr, filesrcptr, filegeoptr, dataptr));
}
