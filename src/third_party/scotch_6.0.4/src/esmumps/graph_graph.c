/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_graph.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These module holds the array graph      **/
/**                building routine.                       **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 14 oct 2003     **/
/**                                 to     22 jan 2004     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     06 dec 2004     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH

#include "common.h"
#ifdef SCOTCH_PTSCOTCH
#include "ptscotch.h"
#else /* SCOTCH_PTSCOTCH */
#include "scotch.h"
#endif /* SCOTCH_PTSCOTCH */
#include "graph.h"

/********************************/
/*                              */
/* The graph handling routines. */
/*                              */
/********************************/

/* This routine builds a graph
** structure from the given
** arrays.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphBuildGraph (
Graph * const               grafptr,              /*+ Graph to build                             +*/
const INT                   baseval,              /*+ Base value                                 +*/
const INT                   vertnbr,              /*+ Number of vertices                         +*/
const INT                   edgenbr,              /*+ Number of arcs                             +*/
INT * restrict              verttab,              /*+ Vertex array                               +*/
INT * restrict              velotab,              /*+ Array of vertex weights (DOFs) if not NULL +*/
INT * restrict              edgetab)              /*+ Edge array                                 +*/
{
  if (sizeof (INT) != sizeof (SCOTCH_Num)) {      /* Check integer consistency */
    errorPrint ("graphBuildGraph: inconsistent integer types");
    return     (1);
  }

  SCOTCH_graphBuild (grafptr, baseval, vertnbr, verttab, NULL, velotab,
                     NULL, edgenbr, edgetab, NULL);

#ifdef GRAPH_DEBUG
  if (graphCheck (grafptr) != 0) {                /* Check graph consistency */
    errorPrint ("graphBuildGraph: inconsistent graph data");
    return     (1);
  }
#endif /* GRAPH_DEBUG */

  return (0);
}

/* This routine builds a graph
** structure from the given
** arrays, with full libScotch
** features.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
graphBuildGraph2 (
Graph * const               grafptr,              /*+ Graph to build                             +*/
const INT                   baseval,              /*+ Base value                                 +*/
const INT                   vertnbr,              /*+ Number of vertices                         +*/
const INT                   edgenbr,              /*+ Number of arcs                             +*/
INT * restrict              verttab,              /*+ Vertex array                               +*/
INT * restrict              vendtab,              /*+ Vertex end array                           +*/
INT * restrict              velotab,              /*+ Array of vertex weights (DOFs) if not NULL +*/
INT * restrict              vlbltab,              /*+ Array of vertex labels if not NULL         +*/
INT * restrict              edgetab,              /*+ Edge array                                 +*/
INT * restrict              edlotab)              /*+ Edge load array                            +*/
{
  if (sizeof (INT) != sizeof (SCOTCH_Num)) {      /* Check integer consistency */
    errorPrint ("graphBuildGraph2: inconsistent integer types");
    return     (1);
  }

  SCOTCH_graphBuild (grafptr, baseval, vertnbr, verttab, vendtab, velotab,
                     vlbltab, edgenbr, edgetab, edlotab);

#ifdef GRAPH_DEBUG
  if (graphCheck (grafptr) != 0) {                /* Check graph consistency */
    errorPrint ("graphBuildGraph2: inconsistent graph data");
    return     (1);
  }
#endif /* GRAPH_DEBUG */

  return (0);
}
