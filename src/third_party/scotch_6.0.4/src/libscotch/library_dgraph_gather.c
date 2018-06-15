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
/**   NAME       : library_dgraph_gather.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source graph handling routines of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 12 jul 2007     **/
/**                                 to     17 jul 2007     **/
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

/*+ This routine gathers the data of a
*** distributed graph on a centralized graph.
*** It returns:
*** - 0   : if the centralization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphGather (
const SCOTCH_Dgraph * const dgrfptr,
SCOTCH_Graph * const        cgrfptr)
{
  Dgraph * restrict   srcdgrfptr;
  Gnum                reduloctab[3];
  Gnum                reduglbtab[3];

  srcdgrfptr = (Dgraph *) dgrfptr;

  if ((cgrfptr != NULL) && (((void *) cgrfptr) != ((void *) dgrfptr))) { /* If centralized graph provided */
    reduloctab[0] = 1;                            /* Process is a potential root                          */
    reduloctab[1] = (Gnum) srcdgrfptr->proclocnum;
  }
  else {                                          /* Process is not a root */
    reduloctab[0] = 0;
    reduloctab[1] = 0;
  }

  if (srcdgrfptr->edloloctax == NULL)             /* Compute sum of edge loads for access to low-level routines */
    reduloctab[2] = srcdgrfptr->edgelocnbr;
  else {
    Gnum                vertlocnum;
    Gnum                edlolocsum;

    for (vertlocnum = srcdgrfptr->baseval, edlolocsum = 0;
         vertlocnum < srcdgrfptr->vertlocnnd; vertlocnum ++) {
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      for (edgelocnum = srcdgrfptr->vertloctax[vertlocnum],
           edgelocnnd = srcdgrfptr->vendloctax[vertlocnum];
           edgelocnum < edgelocnnd; edgelocnum ++)
        edlolocsum += srcdgrfptr->edloloctax[edgelocnum];
    }
    reduloctab[2] = edlolocsum;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_SUM, srcdgrfptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("SCOTCH_dgraphGather: communication error");
    return     (1);
  }
  if (reduglbtab[0] == 1)                         /* If only one single root */
    return (dgraphGatherAll2 (srcdgrfptr, (Graph *) cgrfptr, reduglbtab[2], (int) reduglbtab[1]));
  else if (reduglbtab[0] == srcdgrfptr->procglbnbr) /* If all processes are roots */
    return (dgraphGatherAll2 (srcdgrfptr, (Graph *) cgrfptr, reduglbtab[2], -1));

  errorPrint ("SCOTCH_dgraphGather: invalid number of roots");
  return     (1);
}
