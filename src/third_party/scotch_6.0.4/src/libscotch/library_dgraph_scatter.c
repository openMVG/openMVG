/* Copyright 2007,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_scatter.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source graph handling routines of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 28 apr 2006     **/
/**                                 to     12 jul 2007     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/*+ This routine scatters the data of a
*** centralized graph on a distributed graph.
*** It returns:
*** - 0   : if the scattering succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphScatter (
SCOTCH_Dgraph * const       dgrfptr,
const SCOTCH_Graph * const  cgrfptr)
{
  const Graph *       srccgrfptr;

  srccgrfptr = (((void *) cgrfptr) == ((void *) dgrfptr)) ? NULL : (const Graph *) cgrfptr; /* Consider same pointers as flag for non-root process */

  return (dgraphScatter ((Dgraph *) dgrfptr, srccgrfptr));
}
