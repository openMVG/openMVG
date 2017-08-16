/* Copyright 2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_color.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the graph    **/
/**                coloring routine of the libSCOTCH       **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 02 jan 2012     **/
/**                                 to     19 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "graph_coarsen.h"
#include "scotch.h"

/*********************************/
/*                               */
/* This routine is the C API for */
/* the graph coloring routine.   */
/*                               */
/*********************************/

/*+ This routine creates a color array for the
*** given graph.
*** It returns:
*** - 0  : if the graph has been coarsened.
*** - 1  : if the graph could not be coarsened.
*** - 2  : on error.
+*/

int
SCOTCH_graphColor (
const SCOTCH_Graph * restrict const grafptr,      /* Graph to color              */
SCOTCH_Num * restrict const         colotab,      /* Pointer to color array      */
SCOTCH_Num * restrict const         coloptr,      /* Pointer to number of colors */
const SCOTCH_Num                    flagval)      /* Flag value (not used)       */
{
  Gnum                baseval;
  Gnum                vertnum;
  Gnum                vertnbr;
  Gnum                vertnnd;
  Gnum                queunnd;
  Gnum * restrict     queutax;
  Gnum * restrict     randtax;
  Gnum                colonum;
  Gnum * restrict     colotax;

  const Gnum * restrict const verttax = ((Graph *) grafptr)->verttax;
  const Gnum * restrict const vendtax = ((Graph *) grafptr)->vendtax;
  const Gnum * restrict const edgetax = ((Graph *) grafptr)->edgetax;

  baseval = ((Graph *) grafptr)->baseval;
  vertnbr = ((Graph *) grafptr)->vertnbr;
  vertnnd = vertnbr + baseval;

  memSet (colotab, ~0, vertnbr * sizeof (Gnum));
  colotax = ((Gnum *) colotab) - baseval;

  if (memAllocGroup ((void **) (void *)
                     &queutax, (size_t) (vertnbr * sizeof (Gnum)),
                     &randtax, (size_t) (vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("SCOTCH_graphColor: out of memory");
    return     (1);
  }
  queutax -= baseval;
  randtax -= baseval;

  intRandInit ();                                 /* Check that random number generator is initialized */

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++)
    randtax[vertnum] = intRandVal (32768);

  queunnd = vertnnd;
  for (colonum = 0; queunnd > baseval; colonum ++) { /* Color numbers are not based */
    Gnum                queuold;
    Gnum                queunew;

    for (queunew = queuold = baseval; queuold < queunnd; queuold ++) {
      Gnum                vertnum;
      Gnum                edgenum;
      Gnum                edgennd;
      Gnum                randval;

      vertnum = (queunnd == vertnnd) ? queuold : queutax[queuold];
      randval = randtax[vertnum];
      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
        Gnum                vertend;
        Gnum                randend;

        vertend = edgetax[edgenum];

        if (colotax[vertend] >= 0)
          continue;

        randend = randtax[vertend];
        if ((randend > randval) ||
            ((randend == randval) && (vertend > vertnum))) /* Tie breaking when same random value */
          break;
      }
      if (edgenum >= edgennd)
        colotax[vertnum] = colonum;
      else
        queutax[queunew ++] = vertnum;
    }
    queunnd = queunew;
  }

  *coloptr = colonum;                             /* Set number of colors found */

  memFree (queutax + baseval);

  return (0);
}
