/* Copyright 2007,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hdgraph_gather.c                        **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized halo graph by      **/
/**                gathering the pieces of a distributed   **/
/**                halo graph.                             **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 19 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to   : 30 jul 2010     **/
/**                # Version 6.0  : from : 27 nov 2012     **/
/**                                 to     27 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define HDGRAPH_GATHER

#include "module.h"
#include "common.h"
#include "comm.h"
#include "graph.h"
#include "hgraph.h"
#include "dgraph.h"
#include "hdgraph.h"

/******************************/
/*                            */
/* These routines handle halo */
/* distributed source graphs. */
/*                            */
/******************************/

/* This function gathers the pieces of
** a distributed halo graph to build a
** centralized halo graph.
** There is no gathered vnumtab array if
** the original graph did not have one, as
** vertices are gathered in global order, or
** else the original vnumloctab is gathered.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
hdgraphGather (
Hdgraph * restrict const    dgrfptr,              /* Distributed halo graph */
Hgraph * restrict const     cgrfptr)              /* Centralized halo graph */
{
  Gnum               vertlocnum;
  Gnum               vertlocadj;                  /* Local vertex array adjust                  */
  Gnum               vhallocnnd;
  int                vhallocnbr;                  /* Local copy for sending as a MPI_INT        */
  Gnum * restrict    verthaltax;
  Gnum * restrict    edgehaltax;
  Gnum               edgehalnum;
  int                ehallocnbr;                  /* Local copy for sending as a MPI_INT        */
  int                rootnum;                     /* Index of root process                      */
  Gnum               reduloctab[4];               /* Arrays for reductions                      */
  Gnum               reduglbtab[4];
  int * restrict     recvcnttab;                  /* Arrays for parametrizing gather operations */
  int * restrict     recvdsptab;
  int                cheklocval;
  int                chekglbval;
  Gnum               degrmax;

  if (cgrfptr != NULL) {                          /* If centralized graph provided */
    reduloctab[0] = 1;                            /* This process is the root      */
    reduloctab[1] = (Gnum) dgrfptr->s.proclocnum; /* Get its rank                  */
  }
  else {
    reduloctab[0] =                               /* This process is not the root */
    reduloctab[1] = 0;
  }
  reduloctab[2] = dgrfptr->vhallocnbr;
  reduloctab[3] = dgrfptr->ehallocnbr;
  if (MPI_Allreduce (reduloctab, reduglbtab, 4, GNUM_MPI, MPI_SUM, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphGather: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("hdgraphGather: should have only one root");
    return     (1);
  }
  rootnum = (int) reduglbtab[1];                  /* Get rank of root process                           */
  degrmax = dgrfptr->s.degrglbmax;                /* Distributed degree does not account for halo edges */

  cheklocval = 0;
  if (cgrfptr != NULL) {                          /* If process is root */
    Gnum               vnohnbr;
    Gnum               vertnbr;
    Gnum               velonbr;
    Gnum               vnumnbr;
    Gnum * restrict    velotax;
    Gnum * restrict    vnumtax;
    Gnum               edgenbr;

    vnohnbr = dgrfptr->s.vertglbnbr;
    vertnbr = vnohnbr + reduglbtab[2];
    velonbr = (dgrfptr->s.veloloctax != NULL) ? vertnbr : 0;
    vnumnbr = (dgrfptr->s.vnumloctax != NULL) ? vnohnbr : 0; /* Vertex numbers only serve for non-halo vertices */
    edgenbr = dgrfptr->s.edgeglbnbr + 2 * reduglbtab[3]; /* Twice since halo vertices will be created for real  */

    cgrfptr->s.flagval = GRAPHFREEEDGE | GRAPHEDGEGROUP | GRAPHFREEVERT | GRAPHVERTGROUP; /* In case of premature freeing on error */
    recvcnttab = NULL;
    if (memAllocGroup ((void **) (void *)
                       &cgrfptr->s.verttax, (size_t) ((vertnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                       &velotax,            (size_t) (velonbr       * sizeof (Gnum)),
                       &vnumtax,            (size_t) (vnumnbr       * sizeof (Gnum)),
                       &cgrfptr->vnhdtax,   (size_t) (vnohnbr       * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("hdgraphGather: out of memory (1)");
      cheklocval = 1;
    }
    else if (cgrfptr->s.verttax -= dgrfptr->s.baseval,
             cgrfptr->s.velotax  = (dgrfptr->s.veloloctax != NULL) ? velotax - dgrfptr->s.baseval : NULL,
             cgrfptr->s.vnumtax  = (dgrfptr->s.vnumloctax != NULL) ? vnumtax - dgrfptr->s.baseval : NULL,
             cgrfptr->vnhdtax   -= dgrfptr->s.baseval,
             ((cgrfptr->s.edgetax = (Gnum *) memAlloc (edgenbr * sizeof (Gnum))) == NULL)) {
      errorPrint ("hdgraphGather: out of memory (2)");
      cheklocval = 1;
    }
    else if (cgrfptr->s.edgetax -= dgrfptr->s.baseval,
             memAllocGroup ((void **) (void *)
                            &recvcnttab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)),
                            &recvdsptab, (size_t) (dgrfptr->s.procglbnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("hdgraphGather: out of memory (3)");
      cheklocval = 1;
    }
    else {
      cgrfptr->s.baseval = dgrfptr->s.baseval;
      cgrfptr->s.vertnbr = vertnbr;
      cgrfptr->s.vertnnd = vertnbr + dgrfptr->s.baseval;
      cgrfptr->s.vendtax = cgrfptr->s.verttax + 1; /* Compact edge array                                   */
      cgrfptr->s.velosum = dgrfptr->s.veloglbsum + reduglbtab[2]; /* Halo vertices have unity vertex loads */
      cgrfptr->s.vlbltax = NULL;
      cgrfptr->s.edgenbr = edgenbr;
      cgrfptr->s.edlotax = NULL;
      cgrfptr->s.edlosum = edgenbr;
      cgrfptr->s.procptr = NULL;                  /* Not a multi-sequential gather: no communication possible */
      cgrfptr->vnohnbr   = vnohnbr;
      cgrfptr->vnohnnd   = vnohnbr + dgrfptr->s.baseval;
      cgrfptr->vnlosum   = dgrfptr->s.veloglbsum;
      cgrfptr->enohnbr   =
      cgrfptr->enohsum   = dgrfptr->s.edgeglbnbr;
      cgrfptr->levlnum   = dgrfptr->levlnum;
    }
  }
  if ((cheklocval == 0) &&
      (memAllocGroup ((void **) (void *)
                      &verthaltax, (size_t) (dgrfptr->vhallocnbr * sizeof (Gnum)),
                      &edgehaltax, (size_t) (dgrfptr->ehallocnbr * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("hdgraphGather: out of memory (4)");
    cheklocval = 1;
  }
  else {
    verthaltax -= dgrfptr->s.baseval;
    edgehaltax -= dgrfptr->s.baseval;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphGather: communication error (2)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (verthaltax != NULL)
      memFree (verthaltax + dgrfptr->s.baseval);
    if (cgrfptr != NULL) {                        /* If data were previously allocated */
      if (recvcnttab != NULL)
        memFree (recvcnttab);
      hgraphExit (cgrfptr);
    }
    return (1);
  }

  if (dgrfptr->vhndloctax == dgrfptr->s.vertloctax + 1) { /* If distributed halo graph is compact */
    Gnum               procglbnum;
    Gnum               edgenum;

    if (cgrfptr != NULL) {
      Gnum               vertnum;

      cgrfptr->s.verttax[dgrfptr->s.baseval] = dgrfptr->s.baseval;
      if (commGatherv (dgrfptr->s.vertloctax + 1 + dgrfptr->s.baseval, /* Do not send first index, it is always equal to baseval */
                       dgrfptr->s.vertlocnbr, GNUM_MPI,
                       cgrfptr->s.verttax + 1,      /* First index will always be equal to baseval too, and procdsptab holds based values */
                       dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (3)");
        return     (1);
      }
      if (commGatherv (dgrfptr->s.vendloctax + dgrfptr->s.baseval,
                       dgrfptr->s.vertlocnbr, GNUM_MPI, cgrfptr->vnhdtax, /* procdsptab holds based values */
                       dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (4)");
        return     (1);
      }

      for (procglbnum = 1, vertnum = dgrfptr->s.procdsptab[1] + 1; /* Adjust index sub-arrays for all processors except the first one */
           procglbnum < dgrfptr->s.procglbnbr; procglbnum ++) {
        Gnum               vertnnd;
        Gnum               edgeadj;

        for (vertnnd = dgrfptr->s.procdsptab[procglbnum + 1] + 1,
             edgeadj = cgrfptr->s.verttax[vertnum - 1] - dgrfptr->s.baseval;
             vertnum < vertnnd; vertnum ++) {
          cgrfptr->s.verttax[vertnum]   += edgeadj;
          cgrfptr->vnhdtax[vertnum - 1] += edgeadj;
        }
      }

      for (procglbnum = 0, edgenum = dgrfptr->s.baseval; /* Build arrays for MPI_Gatherv on edge arrays */
           procglbnum < dgrfptr->s.procglbnbr; procglbnum ++) {
        recvcnttab[procglbnum] = cgrfptr->s.verttax[dgrfptr->s.procdsptab[procglbnum + 1]] -
                                 cgrfptr->s.verttax[dgrfptr->s.procdsptab[procglbnum]]; /* verttax used twice since centralized graph is compact */
        recvdsptab[procglbnum] = edgenum;
        edgenum += recvcnttab[procglbnum];
      }
      if (MPI_Gatherv (dgrfptr->s.edgeloctax + dgrfptr->s.baseval, /* Gather edge arrays with global vertex indices */
                       (int) (dgrfptr->s.edgelocnbr + dgrfptr->ehallocnbr), GNUM_MPI, cgrfptr->s.edgetax,
                       recvcnttab, recvdsptab, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (5)");
        return     (1);
      }
    }
    else {
      if (MPI_Gatherv (dgrfptr->s.vertloctax + 1 + dgrfptr->s.baseval, /* Do not send first index, it is always equal to baseval */
                       (int) dgrfptr->s.vertlocnbr, GNUM_MPI, NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (6)");
        return     (1);
      }
      if (MPI_Gatherv (dgrfptr->s.vendloctax + dgrfptr->s.baseval,
                       (int) dgrfptr->s.vertlocnbr, GNUM_MPI, NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (7)");
        return     (1);
      }
      if (MPI_Gatherv (dgrfptr->s.edgeloctax + dgrfptr->s.baseval, /* Gather edge arrays with global vertex indices */
                       (int) (dgrfptr->s.edgelocnbr + dgrfptr->ehallocnbr), GNUM_MPI,
                       NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (8)");
        return     (1);
      }
    }
  }
  else {
    errorPrint ("hdgraphGather: Not implemented"); /* Not really necessary as all Hdgraph structures created by Scotch itself are compact */
    return     (1);
  }

  memSet (verthaltax + dgrfptr->s.baseval, 0, dgrfptr->vhallocnbr * sizeof (Gnum)); /* Initialize halo end vertex count array */
  for (vertlocnum = dgrfptr->s.baseval; vertlocnum < dgrfptr->s.vertlocnnd; vertlocnum ++) {
    Gnum               edgelocnum;

    for (edgelocnum = dgrfptr->s.vendloctax[vertlocnum]; edgelocnum < dgrfptr->vhndloctax[vertlocnum]; edgelocnum ++)
      verthaltax[dgrfptr->s.edgeloctax[edgelocnum]] ++; /* One more edge to this halo vertex */
  }
  vhallocnbr = (int) dgrfptr->vhallocnbr;
  if (MPI_Gather (&vhallocnbr, 1, MPI_INT, recvcnttab, 1, MPI_INT, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphGather: communication error (9)");
    return     (1);
  }
  if (cgrfptr != NULL) {                          /* Build gather parameter array to receive halo edge counts */
    Gnum               procglbnum;
    Gnum               vertnum;

    for (procglbnum = 0, vertnum = 0; procglbnum < dgrfptr->s.procglbnbr; procglbnum ++) { /* Displacements start from zero because adjusted later */
      recvdsptab[procglbnum] = vertnum;
      vertnum += recvcnttab[procglbnum];
    }
    if (MPI_Gatherv (verthaltax + dgrfptr->s.baseval, (int) dgrfptr->vhallocnbr, GNUM_MPI, /* Gather count arrays of halo vertices */
                     cgrfptr->s.verttax + cgrfptr->vnohnnd + 1, recvcnttab, recvdsptab, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("hdgraphGather: communication error (10)");
      return     (1);
    }

    for (procglbnum = 0, vertnum = dgrfptr->s.baseval; /* Adjust end vertex indices for halo edges */
         procglbnum < dgrfptr->s.procglbnbr; procglbnum ++) {
      Gnum               vertnnd;
      Gnum               vertadj;

      for (vertnnd = dgrfptr->s.procdsptab[procglbnum + 1], vertadj = cgrfptr->vnohnbr + recvdsptab[procglbnum];
           vertnum < vertnnd; vertnum ++) {
        Gnum               edgenum;

        if (degrmax < (cgrfptr->s.vendtax[vertnum] - cgrfptr->s.verttax[vertnum])) /* Account for halo edges in maximum degree */
          degrmax = (cgrfptr->s.vendtax[vertnum] - cgrfptr->s.verttax[vertnum]);
        for (edgenum = cgrfptr->vnhdtax[vertnum]; edgenum < cgrfptr->s.vendtax[vertnum]; edgenum ++)
          cgrfptr->s.edgetax[edgenum] += vertadj;
      }
    }
  }
  else {
    if (MPI_Gatherv (verthaltax + dgrfptr->s.baseval, (int) dgrfptr->vhallocnbr, GNUM_MPI, /* Gather count arrays of halo vertices */
                     NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("hdgraphGather: communication error (11)");
      return     (1);
    }
  }
  for (vertlocnum = edgehalnum = dgrfptr->s.baseval, vhallocnnd = dgrfptr->vhallocnbr + dgrfptr->s.baseval;
       vertlocnum < vhallocnnd; vertlocnum ++) { /* Prepare index array for edge collection */
    Gnum               degrlocval;

    degrlocval = verthaltax[vertlocnum];
    verthaltax[vertlocnum] = edgehalnum;
    edgehalnum += degrlocval;
  }
  vertlocadj = dgrfptr->s.procdsptab[dgrfptr->s.proclocnum] - dgrfptr->s.baseval;
  for (vertlocnum = dgrfptr->s.baseval; vertlocnum < dgrfptr->s.vertlocnnd; vertlocnum ++) { /* Collect halo edge ends */
    Gnum               edgelocnum;

    for (edgelocnum = dgrfptr->s.vendloctax[vertlocnum]; edgelocnum < dgrfptr->vhndloctax[vertlocnum]; edgelocnum ++)
      edgehaltax[verthaltax[dgrfptr->s.edgeloctax[edgelocnum]] ++] = vertlocnum + vertlocadj;
  }
  ehallocnbr = (int) dgrfptr->ehallocnbr;
  if (MPI_Gather (&ehallocnbr, 1, MPI_INT, recvcnttab, 1, MPI_INT, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) { /* Gather halo edge counts */
    errorPrint ("hdgraphGather: communication error (12)");
    return     (1);
  }
  if (cgrfptr != NULL) {                          /* Compute receive arrays for edge sub-arrays of halo vertices */
    Gnum               procglbnum;
    Gnum               edgeadj;

    for (procglbnum = 0, edgeadj = 0; procglbnum < dgrfptr->s.procglbnbr; procglbnum ++) {
      recvdsptab[procglbnum] = edgeadj;
      edgeadj += recvcnttab[procglbnum];
    }
    if (MPI_Gatherv (edgehaltax + dgrfptr->s.baseval, (int) dgrfptr->ehallocnbr, GNUM_MPI, /* Gather edge arrays of halo vertices */
                     cgrfptr->s.edgetax + cgrfptr->enohnbr + reduglbtab[3] + dgrfptr->s.baseval, recvcnttab, recvdsptab, GNUM_MPI,
                     rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("hdgraphGather: communication error (13)");
      return     (1);
    }
  }
  else {
    if (MPI_Gatherv (edgehaltax + dgrfptr->s.baseval, (int) dgrfptr->ehallocnbr, GNUM_MPI, /* Gather edge arrays of halo vertices */
                     NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("hdgraphGather: communication error (14)");
      return     (1);
    }
  }

  memFree (verthaltax + dgrfptr->s.baseval);      /* Free group leader */

  if (cgrfptr != NULL) {                          /* Finalize vertex and edge arrays of centralized graph */
    Gnum               vertnum;
    Gnum               edgeadj;

    if (dgrfptr->s.veloloctax != NULL) {            /* Get vertex loads if any */
      if (commGatherv (dgrfptr->s.veloloctax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GNUM_MPI,
                       cgrfptr->s.velotax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GNUM_MPI,
                       rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (15)");
        return     (1);
      }

      for (vertnum = cgrfptr->vnohnnd; vertnum < cgrfptr->s.vertnnd; vertnum ++) /* complete filling of vertex load array */
        cgrfptr->s.velotax[vertnum] = 1;
    }
    if (dgrfptr->s.vnumloctax != NULL) {            /* Get vertex numbers if any */
      if (commGatherv (dgrfptr->s.vnumloctax + dgrfptr->s.baseval, dgrfptr->s.vertlocnbr, GNUM_MPI,
                       cgrfptr->s.vnumtax, dgrfptr->s.proccnttab, dgrfptr->s.procdsptab, GNUM_MPI,
                       rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (16)");
        return     (1);
      }
    }

    memFree (recvcnttab);                         /* Free group leader */

    for (vertnum = cgrfptr->vnohnnd + 1, edgeadj = cgrfptr->s.verttax[cgrfptr->vnohnnd]; /* Adjust vertex array for halo vertices */
         vertnum <= cgrfptr->s.vertnnd; vertnum ++) {
      Gnum               degrval;

      degrval = cgrfptr->s.verttax[vertnum];
      if (degrmax < degrval)                      /* Account for halo edges in maximum degree */
        degrmax = degrval;
      edgeadj += degrval;
      cgrfptr->s.verttax[vertnum] = edgeadj;
    }
    cgrfptr->s.degrmax = degrmax;
  }
  else {
    if (dgrfptr->s.veloloctax != NULL) {          /* Get vertex loads if any */
      if (MPI_Gatherv (dgrfptr->s.veloloctax + dgrfptr->s.baseval, (int) dgrfptr->s.vertlocnbr, GNUM_MPI,
                       NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (17)");
        return     (1);
      }
    }
    if (dgrfptr->s.vnumloctax != NULL) {          /* Get vertex numbers if any */
      if (MPI_Gatherv (dgrfptr->s.vnumloctax + dgrfptr->s.baseval, (int) dgrfptr->s.vertlocnbr, GNUM_MPI,
                       NULL, NULL, NULL, GNUM_MPI, rootnum, dgrfptr->s.proccomm) != MPI_SUCCESS) {
        errorPrint ("hdgraphGather: communication error (18)");
        return     (1);
      }
    }
  }

#ifdef SCOTCH_DEBUG_HDGRAPH2
  cheklocval = (cgrfptr != NULL) ? hgraphCheck (cgrfptr) : 0;
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("hdgraphGather: communication error (19)");
    return     (1);
  }
  if (chekglbval != 0) {
    errorPrint ("hdgraphGather: internal error");
    if (cgrfptr != NULL)
      hgraphExit (cgrfptr);
    return (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  return (0);
}
