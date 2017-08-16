/* Copyright 2004,2007,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_base.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the graph base     **/
/**                changing routine.                       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     31 oct 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     08 jun 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     15 sep 1998     **/
/**                # Version 3.3  : from : 22 sep 1998     **/
/**                                 to     31 dec 1998     **/
/**                # Version 4.0  : from : 24 nov 2001     **/
/**                                 to     22 apr 2004     **/
/**                # Version 6.0  : from : 05 aug 2014     **/
/**                                 to     05 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH

#include "module.h"
#include "common.h"
#include "graph.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

/* This routine sets the base of the given
** graph to the given base value, and returns
** the old base value.
** It returns:
** - old base value : in all cases.
*/

Gnum
graphBase (
Graph * const               grafptr,
const Gnum                  baseval)
{
  Gnum                baseold;                    /* Old base value  */
  Gnum                baseadj;                    /* Base adjustment */
  Gnum                vertnum;
  Gnum                edgenum;

  if (grafptr->baseval == baseval)                /* If nothing to do */
    return (baseval);

  baseold = grafptr->baseval;                     /* Record old base value */
  baseadj = baseval - baseold;                    /* Compute adjustment    */

  for (vertnum = grafptr->baseval; vertnum < grafptr->vertnnd; vertnum ++) {
    for (edgenum = grafptr->verttax[vertnum]; edgenum < grafptr->vendtax[vertnum]; edgenum ++)
      grafptr->edgetax[edgenum] += baseadj;
    grafptr->verttax[vertnum] += baseadj;
  }
  if (grafptr->vendtax != grafptr->verttax + 1) { /* If distinct vertex end array */
    for (vertnum = grafptr->baseval; vertnum < grafptr->vertnnd; vertnum ++)
      grafptr->vendtax[vertnum] += baseadj;
  }
  else                                            /* If same vertex end array (of size +1) */
    grafptr->verttax[grafptr->vertnnd] += baseadj; /* Adjust last entry of verttax         */

  grafptr->verttax -= baseadj;                    /* Adjust array accesses */
  grafptr->vendtax -= baseadj;
  grafptr->edgetax -= baseadj;

  if (grafptr->velotax != NULL)
    grafptr->velotax -= baseadj;
  if (grafptr->vnumtax != NULL)
    grafptr->vnumtax -= baseadj;
  if (grafptr->vlbltax != NULL)
    grafptr->vlbltax -= baseadj;
  if (grafptr->edlotax != NULL)
    grafptr->edlotax -= baseadj;

  grafptr->baseval  = baseval;                    /* Set new base value */
  grafptr->vertnnd += baseadj;

  return (baseold);                               /* Return old base value */
}
