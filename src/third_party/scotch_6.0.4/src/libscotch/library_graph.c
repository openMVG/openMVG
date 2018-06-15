/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_graph.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the source   **/
/**                graph handling routines of the          **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 18 aug 1998     **/
/**                                 to     18 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     09 dec 2005     **/
/**                # Version 5.0  : from : 10 sep 2006     **/
/**                                 to     03 apr 2008     **/
/**                # Version 5.1  : from : 17 nov 2010     **/
/**                                 to     17 nov 2010     **/
/**                # Version 6.0  : from : 04 dec 2012     **/
/**                                 to     04 dec 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_io.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** centralized graph structure.
*** It returns:
*** - !NULL  : if the initialization succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Graph *
SCOTCH_graphAlloc ()
{
  return ((SCOTCH_Graph *) memAlloc (sizeof (SCOTCH_Graph)));
}

/*+ This routine initializes the opaque
*** graph structure used to handle graphs
*** in the Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_graphInit (
SCOTCH_Graph * const        grafptr)
{
  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint ("SCOTCH_graphInit: internal error (1)");
    return     (1);
  }
  if (sizeof (SCOTCH_Graph) < sizeof (Graph)) {
    errorPrint ("SCOTCH_graphInit: internal error (2)");
    return     (1);
  }

  return (graphInit ((Graph *) grafptr));
}

/*+ This routine frees the contents of the
*** given opaque graph structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphExit (
SCOTCH_Graph * const        grafptr)
{
  graphExit ((Graph *) grafptr);
}

/*+ This routine frees the contents of the
*** given opaque graph structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphFree (
SCOTCH_Graph * const        grafptr)
{
  graphFree ((Graph *) grafptr);
}

/*+ This routine loads the given opaque graph
*** structure with the data of the given stream.
*** The base value allows the user to set the
*** graph base to 0 or 1, or to the base value
*** of the stream if the base value is equal
*** to -1. On input, vertex loads are discarded if
*** flagval is 1, edge loads are discarded if flagval
*** is 2, and both if flagval is set to 3.
*** It returns:
*** - 0   : if the loading succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_graphLoad (
SCOTCH_Graph * const        grafptr,
FILE * const                stream,
const SCOTCH_Num            baseval,
const SCOTCH_Num            flagval)
{
  GraphFlag           srcgrafflag;                /* Graph flags */

  if ((baseval < -1) || (baseval > 1)) {
    errorPrint ("SCOTCH_graphLoad: invalid base parameter");
    return     (1);
  }
  if ((flagval < 0) || (flagval > 3)) {
    errorPrint ("SCOTCH_graphLoad: invalid flag parameter");
    return     (1);
  }

  srcgrafflag = (((flagval & 1) != 0) ? GRAPHIONOLOADVERT : 0) +
                (((flagval & 2) != 0) ? GRAPHIONOLOADEDGE : 0);

  return (graphLoad ((Graph * const) grafptr, stream, (Gnum) baseval, srcgrafflag));
}

/*+ This routine saves the contents of the given
*** opaque graph structure to the given stream.
*** It returns:
*** - 0   : if the saving succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_graphSave (
const SCOTCH_Graph * const  grafptr,
FILE * const                stream)
{
  return (graphSave ((const Graph * const) grafptr, stream));
}

/*+ This routine fills the contents of the given
*** opaque graph structure with the data provided
*** by the user. The base value allows the user to
*** set the graph base to 0 or 1.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphBuild (
SCOTCH_Graph * const        grafptr,              /* Graph structure to fill             */
const SCOTCH_Num            baseval,              /* Base value                          */
const SCOTCH_Num            vertnbr,              /* Number of vertices                  */
const SCOTCH_Num * const    verttab,              /* Vertex array [vertnbr or vertnbr+1] */
const SCOTCH_Num * const    vendtab,              /* Vertex end array [vertnbr]          */
const SCOTCH_Num * const    velotab,              /* Vertex load array                   */
const SCOTCH_Num * const    vlbltab,              /* Vertex label array                  */
const SCOTCH_Num            edgenbr,              /* Number of edges (arcs)              */
const SCOTCH_Num * const    edgetab,              /* Edge array [edgenbr]                */
const SCOTCH_Num * const    edlotab)              /* Edge load array                     */
{
  Graph *             srcgrafptr;                 /* Pointer to source graph structure */
  Gnum                vertnum;                    /* Current vertex number             */
  Gnum                degrmax;                    /* Maximum degree                    */

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Graph) < sizeof (Graph)) {
    errorPrint ("SCOTCH_graphBuild: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
  if ((baseval < 0) || (baseval > 1)) {
    errorPrint ("SCOTCH_graphBuild: invalid base parameter");
    return     (1);
  }

  srcgrafptr = (Graph *) grafptr;                 /* Use structure as source graph */
  srcgrafptr->flagval = GRAPHNONE;
  srcgrafptr->baseval = baseval;
  srcgrafptr->vertnbr = vertnbr;
  srcgrafptr->vertnnd = vertnbr + baseval;
  srcgrafptr->verttax = (Gnum *) verttab - baseval;
  srcgrafptr->vendtax = ((vendtab == NULL) || (vendtab == verttab)) ? srcgrafptr->verttax + 1 : (Gnum *) vendtab - baseval;
  srcgrafptr->velotax = ((velotab == NULL) || (velotab == verttab)) ? NULL : (Gnum *) velotab - baseval;
  srcgrafptr->vnumtax = NULL;
  srcgrafptr->vlbltax = ((vlbltab == NULL) || (vlbltab == verttab)) ? NULL : (Gnum *) vlbltab - baseval;
  srcgrafptr->edgenbr = edgenbr;
  srcgrafptr->edgetax = (Gnum *) edgetab - baseval;
  srcgrafptr->edlotax = ((edlotab == NULL) || (edlotab == edgetab)) ? NULL : (Gnum *) edlotab - baseval;

  if (srcgrafptr->velotax == NULL)                /* Compute vertex load sum */
    srcgrafptr->velosum = vertnbr;
  else {
    Gnum                velosum;                  /* Sum of vertex loads */

    for (vertnum = srcgrafptr->baseval, velosum = 0; vertnum < srcgrafptr->vertnnd; vertnum ++)
      velosum += srcgrafptr->velotax[vertnum];
    srcgrafptr->velosum = velosum;
  }

  if (srcgrafptr->edlotax == NULL) {              /* If no edge loads       */
    srcgrafptr->edlosum = srcgrafptr->edgenbr;    /* Edge load sum is known */

    for (vertnum = srcgrafptr->baseval, degrmax = 0; /* Compute maximum degree only */
         vertnum < srcgrafptr->vertnnd; vertnum ++) {
      Gnum                degrval;                /* Degree of current vertex */

      degrval = srcgrafptr->vendtax[vertnum] - srcgrafptr->verttax[vertnum];
      if (degrval > degrmax)
        degrmax = degrval;
    }
  }
  else {                                          /* Graph has edge loads, compute edge load sum */
    Gnum                edlosum;

    for (vertnum = srcgrafptr->baseval, edlosum = degrmax = 0;
         vertnum < srcgrafptr->vertnnd; vertnum ++) {
      Gnum                edgenum;
      Gnum                degrval;                /* Degree of current vertex */

      degrval = srcgrafptr->vendtax[vertnum] - srcgrafptr->verttax[vertnum];
      if (degrval > degrmax)
        degrmax = degrval;
      for (edgenum = srcgrafptr->verttax[vertnum]; edgenum < srcgrafptr->vendtax[vertnum]; edgenum ++)
        edlosum += srcgrafptr->edlotax[edgenum];
    }
    srcgrafptr->edlosum = edlosum;
  }
  srcgrafptr->degrmax = degrmax;

  return (0);
}

/*+ This routine accesses graph size data.
*** NULL pointers on input indicate unwanted
*** data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphSize (
const SCOTCH_Graph * const  grafptr,
SCOTCH_Num * const          vertnbr,
SCOTCH_Num * const          edgenbr)
{
  const Graph *       srcgrafptr;

  srcgrafptr = (Graph *) grafptr;

  if (vertnbr != NULL)
    *vertnbr = (SCOTCH_Num) (srcgrafptr->vertnbr);
  if (edgenbr != NULL)
    *edgenbr = (SCOTCH_Num) srcgrafptr->edgenbr;
}

/*+ This routine accesses all of the graph data.
*** NULL pointers on input indicate unwanted
*** data. NULL pointers on output indicate
*** unexisting arrays.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphData (
const SCOTCH_Graph * const  grafptr,              /* Graph structure to read  */
SCOTCH_Num * const          baseptr,              /* Base value               */
SCOTCH_Num * const          vertptr,              /* Number of vertices       */
SCOTCH_Num ** const         verttab,              /* Vertex array [vertnbr+1] */
SCOTCH_Num ** const         vendtab,              /* Vertex array [vertnbr]   */
SCOTCH_Num ** const         velotab,              /* Vertex load array        */
SCOTCH_Num ** const         vlbltab,              /* Vertex label array       */
SCOTCH_Num * const          edgeptr,              /* Number of edges (arcs)   */
SCOTCH_Num ** const         edgetab,              /* Edge array [edgenbr]     */
SCOTCH_Num ** const         edlotab)              /* Edge load array          */
{
  const Graph *       srcgrafptr;                 /* Pointer to source graph structure */

  srcgrafptr = (const Graph *) grafptr;

  if (baseptr != NULL)
    *baseptr = srcgrafptr->baseval;
  if (vertptr != NULL)
    *vertptr = srcgrafptr->vertnbr;
  if (verttab != NULL)
    *verttab = srcgrafptr->verttax + srcgrafptr->baseval;
  if (vendtab != NULL)
    *vendtab = srcgrafptr->vendtax + srcgrafptr->baseval;
  if (velotab != NULL)
    *velotab = (srcgrafptr->velotax != NULL) ? srcgrafptr->velotax + srcgrafptr->baseval : NULL;
  if (vlbltab != NULL)
    *vlbltab = (srcgrafptr->vlbltax != NULL) ? srcgrafptr->vlbltax + srcgrafptr->baseval : NULL;
  if (edgeptr != NULL)
    *edgeptr = srcgrafptr->edgenbr;
  if (edgetab != NULL)
    *edgetab = srcgrafptr->edgetax + srcgrafptr->baseval;
  if (edlotab != NULL)
    *edlotab = (srcgrafptr->edlotax != NULL) ? srcgrafptr->edlotax + srcgrafptr->baseval : NULL;
}

/*+ This routine computes statistics
*** on the given graph.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphStat (
const SCOTCH_Graph * const  grafptr,
SCOTCH_Num * const          velominptr,
SCOTCH_Num * const          velomaxptr,
SCOTCH_Num * const          velosumptr,
double *                    veloavgptr,
double *                    velodltptr,
SCOTCH_Num * const          degrminptr,
SCOTCH_Num * const          degrmaxptr,
double *                    degravgptr,
double *                    degrdltptr,
SCOTCH_Num * const          edlominptr,
SCOTCH_Num * const          edlomaxptr,
SCOTCH_Num * const          edlosumptr,
double *                    edloavgptr,
double *                    edlodltptr)
{
  const Graph *       srcgrafptr;
  Gnum                vertnum;
  Gnum                vertnbr;
  Gnum                velomin;
  Gnum                velomax;
  double              veloavg;
  double              velodlt;
  Gnum                degrval;
  Gnum                degrmin;
  Gnum                degrmax;
  double              degravg;
  double              degrdlt;
  Gnum                edgenum;
  Gnum                edlomin;
  Gnum                edlomax;
  Gnum                edlosum;
  double              edloavg;
  double              edlodlt;

  srcgrafptr = (Graph *) grafptr;

  vertnbr = srcgrafptr->vertnnd - srcgrafptr->baseval;

  velodlt = 0.0L;
  if (vertnbr > 0) {
    if (srcgrafptr->velotax != NULL) {            /* If graph has vertex loads */
      velomin = GNUMMAX;
      velomax = 0;
      veloavg = (double) srcgrafptr->velosum / (double) vertnbr;

      for (vertnum = srcgrafptr->baseval; vertnum < srcgrafptr->vertnnd; vertnum ++) {
        if (srcgrafptr->velotax[vertnum] < velomin) /* Account for vertex loads */
          velomin = srcgrafptr->velotax[vertnum];
        if (srcgrafptr->velotax[vertnum] > velomax)
          velomax = srcgrafptr->velotax[vertnum];
        velodlt += fabs ((double) srcgrafptr->velotax[vertnum] - veloavg);
      }
      velodlt /= (double) vertnbr;
    }
    else {
      velomin =
      velomax = 1;
      veloavg = 1.0L;
    }
  }
  else {
    velomin =
    velomax = 0;
    veloavg = 0.0L;
  }

  if (velominptr != NULL)
    *velominptr = (SCOTCH_Num) velomin;
  if (velomaxptr != NULL)
    *velomaxptr = (SCOTCH_Num) velomax;
  if (velosumptr != NULL)
    *velosumptr = (SCOTCH_Num) srcgrafptr->velosum;
  if (veloavgptr != NULL)
    *veloavgptr = (double) veloavg;
  if (velodltptr != NULL)
    *velodltptr = (double) velodlt;

  degrmax  = 0;
  degrdlt  = 0.0L;
  if (vertnbr > 0) {
    degrmin = GNUMMAX;
    degravg = (double) srcgrafptr->edgenbr / (double) vertnbr;
    for (vertnum = srcgrafptr->baseval; vertnum < srcgrafptr->vertnnd; vertnum ++) {
      degrval = srcgrafptr->vendtax[vertnum] - srcgrafptr->verttax[vertnum]; /* Get vertex degree */
      if (degrval < degrmin)
        degrmin = degrval;
      if (degrval > degrmax)
        degrmax = degrval;
      degrdlt += fabs ((double) degrval - degravg);
    }
    degrdlt /= (double) vertnbr;
  }
  else {
    degrmin = 0;
    degravg = 0.0L;
  }

  if (degrminptr != NULL)
    *degrminptr = (SCOTCH_Num) degrmin;
  if (degrmaxptr != NULL)
    *degrmaxptr = (SCOTCH_Num) degrmax;
  if (degravgptr != NULL)
    *degravgptr = (double) degravg;
  if (degrdltptr != NULL)
    *degrdltptr = (double) degrdlt;

  edlodlt = 0.0L;
  if (srcgrafptr->edgenbr > 0) {
    if (srcgrafptr->edlotax != NULL) {            /* If graph has edge loads */
      edlomin = GNUMMAX;
      edlomax = 0;
      edlosum = 0;

      for (vertnum = srcgrafptr->baseval; vertnum < srcgrafptr->vertnnd; vertnum ++) {
        for (edgenum = srcgrafptr->verttax[vertnum]; edgenum < srcgrafptr->vendtax[vertnum]; edgenum ++) { /* For all edges */
          if (srcgrafptr->edlotax[edgenum] < edlomin) /* Account for edge load */
            edlomin = srcgrafptr->edlotax[edgenum];
          if (srcgrafptr->edlotax[edgenum] > edlomax)
            edlomax = srcgrafptr->edlotax[edgenum];
          edlosum += srcgrafptr->edlotax[edgenum];
        }
      }
      edloavg = (double) edlosum /
                (double) (2 * srcgrafptr->edgenbr);

      for (vertnum = srcgrafptr->baseval; vertnum < srcgrafptr->vertnnd; vertnum ++) {
        for (edgenum = srcgrafptr->verttax[vertnum]; edgenum < srcgrafptr->vendtax[vertnum]; edgenum ++) /* For all edges */
          edlodlt += fabs ((double) srcgrafptr->edlotax[edgenum] - edloavg);
      }
      edlodlt /= (double) srcgrafptr->edgenbr;
    }
    else {
      edlomin =
      edlomax = 1;
      edlosum = srcgrafptr->edgenbr / 2;
      edloavg = 1.0L;
    }
  }
  else {
    edlomin =
    edlomax = 0;
    edlosum = 0;
    edloavg = 0.0L;
  }

  if (edlominptr != NULL)
    *edlominptr = (SCOTCH_Num) edlomin;
  if (edlomaxptr != NULL)
    *edlomaxptr = (SCOTCH_Num) edlomax;
  if (edlosumptr != NULL)
    *edlosumptr = (SCOTCH_Num) edlosum;
  if (edloavgptr != NULL)
    *edloavgptr = (double) edloavg;
  if (edlodltptr != NULL)
    *edlodltptr = (double) edlodlt;
}
