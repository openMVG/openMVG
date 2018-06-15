/* Copyright 2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_gather.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Francois CHATENET (0.0)                 **/
/**                Sebastien FOUCAULT (0.0)                **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized graph by gathering **/
/**                the pieces of a distributed graph.      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 apr 1997     **/
/**                                 to   : 20 jun 1997     **/
/**                # Version 0.1  : from : 14 apr 1998     **/
/**                                 to   : 20 jun 1998     **/
/**                # Version 0.2  : from : 20 may 1999     **/
/**                                 to   : 20 may 1999     **/
/**                # Version 5.0  : from : 07 feb 2006     **/
/**                                 to   : 16 jul 2007     **/
/**                                                        **/
/**   NOTES      : # The definitions of MPI_Gather and     **/
/**                  MPI_Gatherv indicate that elements in **/
/**                  the receive array should not be       **/
/**                  written more than once. Great care    **/
/**                  should be taken to enforce this rule, **/
/**                  especially when the number of         **/
/**                  vertices in the centralized graph is  **/
/**                  smaller than the number of            **/
/**                  processors.                           **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_GATHER

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph.h"

/* This function gathers the pieces of
** a distributed graph to build a
** centralized graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
dgraphGather (
const Dgraph * restrict const dgrfptr,            /* Distributed graph */
Graph * restrict              cgrfptr)            /* Centralized graph */
{
  Gnum                reduloctab[3];
  Gnum                reduglbtab[3];

  if (dgrfptr->edloloctax == NULL)                /* Compute sum of edge loads */
    reduloctab[2] = dgrfptr->edgelocnbr;
  else {
    Gnum                vertlocnum;
    Gnum                edlolocsum;

    for (vertlocnum = dgrfptr->baseval, edlolocsum = 0;
         vertlocnum < dgrfptr->vertlocnnd; vertlocnum ++) {
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      for (edgelocnum = dgrfptr->vertloctax[vertlocnum],
           edgelocnnd = dgrfptr->vendloctax[vertlocnum];
           edgelocnum < edgelocnnd; edgelocnum ++)
        edlolocsum += dgrfptr->edloloctax[edgelocnum];
    }
    reduloctab[2] = edlolocsum;
  }

  if (cgrfptr != NULL) {                          /* If centralized graph provided */
    reduloctab[0] = 1;                            /* This process is the root      */
    reduloctab[1] = (Gnum) dgrfptr->proclocnum;   /* Get its number                */
  }
  else {
    reduloctab[0] =                               /* This process is not the root */
    reduloctab[1] = 0;
  }
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_SUM, dgrfptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphGather: communication error");
    return     (1);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("dgraphGather: should have only one root");
    return     (1);
  }

  return (dgraphGatherAll2 (dgrfptr, cgrfptr, reduglbtab[2], (int) reduglbtab[1]));
}
