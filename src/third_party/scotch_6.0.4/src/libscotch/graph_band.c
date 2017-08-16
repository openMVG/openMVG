/* Copyright 2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_band.c                            **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a band graph       **/
/**                from the given frontier array.          **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 jan 2010     **/
/**                                 to   : 22 sep 2011     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  dgraph_band.c in version 5.1.         **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_BAND

#include "module.h"
#include "common.h"
#include "graph.h"

/**********************************/
/*                                */
/* Distance computation routines. */
/*                                */
/**********************************/

/* This routine computes an index array
** of given width around the current separator.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
graphBand (
const Graph * restrict const      grafptr,        /*+ Graph                                                    +*/
const Gnum                        queunbr,        /*+ Number of frontier vertices, start size for vertex queue +*/
Gnum * restrict const             queutab,        /*+ Array of frontier vertices, re-used as queue array       +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices                 +*/
Gnum * restrict * restrict const  vnumptr,        /*+ Pointer to vnumtax                                       +*/
Gnum * restrict const             bandvertlvlptr, /*+ Pointer to based start index of last level               +*/
Gnum * restrict const             bandvertptr,    /*+ Pointer to bandvertnbr                                   +*/
Gnum * restrict const             bandedgeptr,    /*+ Pointer to bandedgenbr                                   +*/
const Gnum * restrict const       pfixtax,        /*+ Fixed partition array                                    +*/
Gnum * restrict const             bandvfixptr)    /*+ Pointer to bandvfixnbr                                   +*/
{
  Gnum                    queunum;
  Gnum * restrict         vnumtax;                /* Index array for vertices kept in band graph */
  Gnum                    queuheadidx;            /* Index of head of queue                      */
  Gnum                    queutailidx;            /* Index of tail of queue                      */
  Gnum                    bandvertlvlnum;
  Gnum                    bandvertnum;
  Gnum                    bandedgenbr;
  Gnum                    distval;
  Gnum                    bandvfixnbr;            /* Number of band fixed vertices */

  const Gnum * restrict const verttax = grafptr->verttax;
  const Gnum * restrict const vendtax = grafptr->vendtax;
  const Gnum * restrict const edgetax = grafptr->edgetax;

  if ((vnumtax = memAlloc (grafptr->vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("graphBand: out of memory (1)");
    return     (1);
  }

  bandvertlvlnum =                                /* Start index of last level is start index */
  bandvertnum    = grafptr->baseval;              /* Reset number of band vertices            */
  bandedgenbr    =
  bandvfixnbr    = 0;
  memSet (vnumtax, ~0, grafptr->vertnbr * sizeof (Gnum)); /* Reset part array */
  vnumtax -= grafptr->baseval;

  for (queunum = 0; queunum < queunbr; queunum ++) { /* All frontier vertices will be first vertices of band graph */
    Gnum              vertnum;

    vertnum = queutab[queunum];
    if ((pfixtax != NULL) && (pfixtax[vertnum] != -1)) { /* It is a fixed vertex */
      vnumtax[vertnum] = -2;                      /* Set vertex as fixed         */
      bandvfixnbr ++;
    }
    else
      vnumtax[vertnum] = bandvertnum ++;          /* Keep frontier vertex in band */
    bandedgenbr += vendtax[vertnum] - verttax[vertnum]; /* Account for its edges  */
  }
  queuheadidx = 0;                                /* No queued vertex read yet                  */
  queutailidx = queunbr;                          /* All frontier vertices are already in queue */

  for (distval = 0; ++ distval <= distmax; ) {
    Gnum              queunextidx;                /* Tail index for enqueuing vertices of next band */

    bandvertlvlnum = bandvertnum;
    *bandvertlvlptr = bandvertlvlnum;             /* Save start index of current level, based */

    for (queunextidx = queutailidx; queuheadidx < queutailidx; ) { /* For all vertices in queue */
      Gnum              vertnum;
      Gnum              edgenum;

      vertnum = queutab[queuheadidx ++];          /* Dequeue vertex */
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum              vertend;

        vertend = edgetax[edgenum];
        if (vnumtax[vertend] != ~0)               /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        if ((pfixtax != NULL) && (pfixtax[vertend] != -1)) { /* If fixed vertex */
          vnumtax[vertend] = -2;                  /* Set vertex as fixed        */
          bandvfixnbr ++;
        }
        else
          vnumtax[vertend] = bandvertnum ++;      /* Enqueue vertex label */

        bandedgenbr += vendtax[vertend] - verttax[vertend]; /* Account for its edges */
        queutab[queunextidx ++] = vertend;        /* Enqueue vertex for next pass    */
      }
    }

    queutailidx = queunextidx;                    /* Prepare queue for next sweep */
  }

  *vnumptr     = vnumtax;
  *bandvfixptr = bandvfixnbr;
  *bandvertptr = bandvertnum - grafptr->baseval;
  *bandedgeptr = bandedgenbr;

  return (0);
}
