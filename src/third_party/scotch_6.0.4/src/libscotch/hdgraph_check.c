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
/**   NAME       : hdgraph_check.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel static mapper.       **/
/**                This module contains the distributed    **/
/**                graph consistency checking routine.     **/
/**                                                        **/
/**                # Version 5.0  : from : 21 apr 2006     **/
/**                                 to   : 29 apr 2006     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define HDGRAPH_CHECK

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "hdgraph.h"

/******************************/
/*                            */
/* These routines handle halo */
/* distributed source graphs. */
/*                            */
/******************************/

/* This function checks the consistency
** of the given halo distributed graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
hdgraphCheck (
const Hdgraph * restrict const  grafptr)
{
  Gnum                vertlocnum;
  int * restrict      vhalloctax;                 /* Flag array for halo vertices */
  Gnum                vhallocnnd;
  Gnum                vhallocnum;
  Gnum                ehallocnbr;
  int                 cheklocval;                 /* Local consistency flag       */
  int                 chekglbval;                 /* Global consistency flag      */

  cheklocval = 0;
  for (vertlocnum = grafptr->s.baseval, ehallocnbr = 0; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
    if ((grafptr->vhndloctax[vertlocnum] < grafptr->s.vendloctax[vertlocnum]) ||
        (grafptr->vhndloctax[vertlocnum] > (grafptr->s.edgelocsiz + grafptr->s.baseval))) {
      errorPrint ("hdgraphCheck: inconsistent local vertex arrays");
      cheklocval = 1;
    }
    ehallocnbr += grafptr->vhndloctax[vertlocnum] - grafptr->s.vendloctax[vertlocnum];
  }
  if (ehallocnbr != grafptr->ehallocnbr) {
    errorPrint ("hdgraphCheck: invalid local number of halo edges");
    cheklocval = 1;
  }

  if ((grafptr->vhallocnbr < 0) || (grafptr->vhallocnbr > grafptr->s.edgelocsiz)) {
    errorPrint ("hdgraphCheck: invalid local number of halo vertices");
    cheklocval = 1;
  }

  vhalloctax = NULL;
  if ((cheklocval == 0) &&
      ((vhalloctax = (int *) memAlloc (grafptr->vhallocnbr * sizeof (int))) == NULL)) {
    errorPrint ("hdgraphCheck: out of memory");
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphCheck: communication error (1)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (vhalloctax != NULL)
      memFree (vhalloctax);
    return (1);
  }

  memSet (vhalloctax, ~0, grafptr->vhallocnbr * sizeof (int));
  vhalloctax -= grafptr->s.baseval;
  vhallocnnd  = grafptr->vhallocnbr + grafptr->s.baseval;
  for (vertlocnum = grafptr->s.baseval; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
    Gnum                edgelocnum;

    for (edgelocnum = grafptr->s.vendloctax[vertlocnum];
         edgelocnum < grafptr->vhndloctax[vertlocnum]; edgelocnum ++) {
      Gnum                vhallocend;

      vhallocend = grafptr->s.edgeloctax[edgelocnum];
      if ((vhallocend < grafptr->s.baseval) || (vhallocend >= vhallocnnd)) {
        errorPrint ("hdgraphCheck: invalid halo vertex number");
        vertlocnum = grafptr->s.vertlocnnd;       /* Avoid unwanted cascaded error messages */
        cheklocval = 1;
        break;
      }
      vhalloctax[vhallocend] = 0;                 /* Flag halo vertex as used */
    }
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphCheck: communication error (2)");
    return     (1);
  }
  if (chekglbval != 0) {
    memFree (vhalloctax + grafptr->s.baseval);
    return (1);
  }

  for (vhallocnum = grafptr->s.baseval; vhallocnum < vhallocnnd; vhallocnum ++) {
    if (vhalloctax[vhallocnum] != 0) {            /* If halo vertex index not used in graph */
      errorPrint ("hdgraphCheck: unused halo vertex number");
      cheklocval = 1;
      break;
    }
  }
  memFree (vhalloctax + grafptr->s.baseval);
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphCheck: communication error (3)");
    return     (1);
  }
  if (chekglbval != 0)
    return (1);

  return (dgraphCheck (&grafptr->s));
}
