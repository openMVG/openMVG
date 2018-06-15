/* Copyright 2007,2008,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_gather_all.c                     **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized graph on all       **/
/**                processes by gathering the pieces of    **/
/**                a distributed graph.                    **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 07 feb 2006     **/
/**                                 to     17 jun 2008     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to     30 jul 2010     **/
/**                # Version 6.0  : from : 27 nov 2012     **/
/**                                 to     27 nov 2012     **/
/**                                                        **/
/**   NOTES      : # The definitions of MPI_Gather and     **/
/**                  MPI_Gatherv indicate that elements in **/
/**                  the receive array should not be       **/
/**                  written more than once. Great care    **/
/**                  should be taken to enforce this rule, **/
/**                  especially when the number of         **/
/**                  vertices in the centralized graph is  **/
/**                  smaller than the number of            **/
/**                  processes.                            **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "comm.h"
#include "graph.h"
#include "dgraph.h"

/* This function gathers on all processes
** the pieces of a distributed graph to
** build a centralized graph. This function
** does not compute edlosum on the centralized
** graphs when it is already given in the passed
** value, as a non-negative number.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

static
int
dgraphGatherAll3 (
Gnum * const                senddattab,
const Gnum                  sendcntnbr,
Gnum * const                recvdattab,
Gnum * const                recvcnttab,
Gnum * const                recvdsptab,
const int                   rootnum,
MPI_Comm                    comm)
{
  if (rootnum == -1)                              /* If collective communication wanted */
    return (commAllgatherv (senddattab, sendcntnbr, GNUM_MPI, recvdattab, recvcnttab, recvdsptab, GNUM_MPI, comm));
  else
    return (commGatherv (senddattab, sendcntnbr, GNUM_MPI, recvdattab, recvcnttab, recvdsptab, GNUM_MPI, rootnum, comm));
}

int
dgraphGatherAll2 (
const Dgraph * restrict const dgrfptr,            /* Distributed graph  */
Graph * restrict              cgrfptr,            /* Centralized graph  */
const Gnum                    edlosum,            /* -1 means recompute */
const int                     protnum)            /* -1 means allgather */
{
  Gnum                baseval;
  Gnum * restrict     verttax;                    /* Target vertex array for root, dummy for non-roots        */
  Gnum * restrict     velotax;                    /* Target vertex load array for root, dummy for non-roots   */
  Gnum * restrict     vnumtax;                    /* Target vertex index array for root, dummy for non-roots  */
  Gnum * restrict     vlbltax;                    /* Target vertex label array for root, dummy for non-roots  */
  Gnum * restrict     edgetax;                    /* Target edge array for root, dummy for non-roots          */
  Gnum * restrict     edlotax;                    /* Target edge load array for root, dummy for non-roots     */
  Gnum                vertlocnbr;                 /* Size of temporary distributed vertex array               */
  Gnum * restrict     vertloctax;                 /* Temporary vertex array if graph is not compact           */
  Gnum                edgelocnbr;                 /* Size of temporary distributed edge array                 */
  Gnum * restrict     edgeloctab;                 /* Temporary edge array if distributed graph is not compact */
  Gnum * restrict     recvcnttab;                 /* Count array for gather operations                        */
  Gnum * restrict     recvdsptab;                 /* Displacement array for gather operations                 */
  int                 cheklocval;
  int                 chekglbval;

  const Gnum * restrict const edgeloctax = dgrfptr->edgeloctax;

#ifdef SCOTCH_DEBUG_DGRAPH1
  cheklocval = 0;
  if (cgrfptr != NULL)                            /* Centralized graphs should be provided by all */
    cheklocval = 1;
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, dgrfptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphGatherAll2: communication error (1)");
    return     (1);
  }
  if (protnum == -1) {                            /* If collective gathering wanted */
    if (chekglbval != dgrfptr->procglbnbr) {
      errorPrint ("dgraphGatherAll2: centralized graphs should be provided on every process");
      return     (1);
    }
  }
  else {                                          /* Single gathering wanted */
    if (chekglbval != 1) {
      errorPrint ("dgraphGatherAll2: should have only one root");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  baseval = dgrfptr->baseval;

  cheklocval = 0;
  if (cgrfptr != NULL) {                          /* If root process */
    Gnum                velonbr;
    Gnum                vnumnbr;
    Gnum                vlblnbr;
    Gnum                edlonbr;

    velonbr = (dgrfptr->veloloctax != NULL) ? dgrfptr->vertglbnbr : 0;
    vnumnbr = (dgrfptr->vnumloctax != NULL) ? dgrfptr->vertglbnbr : 0;
    vlblnbr = (dgrfptr->vlblloctax != NULL) ? dgrfptr->vertglbnbr : 0;
    edlonbr = (dgrfptr->edloloctax != NULL) ? dgrfptr->edgeglbnbr : 0;

    if (memAllocGroup ((void **) (void *)
                       &cgrfptr->verttax, (size_t) ((dgrfptr->vertglbnbr + 1) * sizeof (Gnum)),
                       &cgrfptr->velotax, (size_t) (velonbr                   * sizeof (Gnum)),
                       &cgrfptr->vnumtax, (size_t) (vnumnbr                   * sizeof (Gnum)),
                       &cgrfptr->vlbltax, (size_t) (vlblnbr                   * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("dgraphGatherAll2: out of memory (1)");
      cheklocval = 1;
    }
    else if (memAllocGroup ((void **) (void *)
                            &cgrfptr->edgetax, (size_t) (dgrfptr->edgeglbnbr * sizeof (Gnum)),
                            &cgrfptr->edlotax, (size_t) (edlonbr             * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("dgraphGatherAll2: out of memory (2)");
      cheklocval = 1;
    }
  }

  if (dgrfptr->vendloctax == (dgrfptr->vertloctax + 1)) { /* If distributed graph is compact */
    vertlocnbr =                                  /* No need to recompact arrays             */
    edgelocnbr = 0;
  }
  else {                                          /* Need extra space to compact vertex and edge arrays before sending */
    vertlocnbr = dgrfptr->vertlocnbr;
    edgelocnbr = dgrfptr->edgelocnbr;
  }

  if (cheklocval == 0) {
    if (memAllocGroup ((void **) (void *)
                       &recvcnttab, (size_t) (dgrfptr->procglbnbr * sizeof (Gnum)), /* Allocated for non-roots too but don't care as these are very small */
                       &recvdsptab, (size_t) (dgrfptr->procglbnbr * sizeof (Gnum)),
                       &vertloctax, (size_t) (vertlocnbr          * sizeof (Gnum)),
                       &edgeloctab, (size_t) (edgelocnbr          * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("dgraphGatherAll2: out of memory (3)");
      cheklocval = 1;
    }
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be merged with a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphGatherAll2: communication error (2)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0) {
    if (recvcnttab != NULL)
      memFree (recvcnttab);
    if (cgrfptr->verttax != NULL) {
      if (cgrfptr->edgetax != NULL)
        memFree (cgrfptr->edgetax);               /* Arrays are not based yet */
      memFree (cgrfptr->verttax);
    }
    return (1);
  }

  if (cgrfptr != NULL) {
    verttax = cgrfptr->verttax - baseval;
    velotax = (dgrfptr->veloloctax != NULL) ? (cgrfptr->velotax - baseval) : NULL;
    vnumtax = (dgrfptr->vnumloctax != NULL) ? (cgrfptr->vnumtax - baseval) : NULL;
    vlbltax = (dgrfptr->vlblloctax != NULL) ? (cgrfptr->vlbltax - baseval) : NULL;
    edgetax = cgrfptr->edgetax - baseval;
    edlotax = (dgrfptr->edloloctax != NULL) ? (cgrfptr->edlotax - baseval) : NULL;

    cgrfptr->flagval = GRAPHFREEVERT | GRAPHVERTGROUP | GRAPHFREEEDGE | GRAPHEDGEGROUP; /* Other arrays are grouped, too */
    cgrfptr->baseval = baseval;
    cgrfptr->vertnbr = dgrfptr->vertglbnbr;
    cgrfptr->vertnnd = dgrfptr->vertglbnbr + baseval;
    cgrfptr->verttax = verttax;
    cgrfptr->vendtax = verttax + 1;               /* Compact edge array */
    cgrfptr->velotax = velotax;
    cgrfptr->velosum = dgrfptr->veloglbsum;
    cgrfptr->vnumtax = vnumtax;
    cgrfptr->vlbltax = vlbltax;
    cgrfptr->edgenbr = dgrfptr->edgeglbnbr;
    cgrfptr->edgetax = edgetax;
    cgrfptr->edlotax = edlotax;
    cgrfptr->edlosum = edlosum;
    cgrfptr->degrmax = dgrfptr->degrglbmax;
    cgrfptr->procptr = NULL;                      /* This field exists only when compiled with SCOTCH_PTSCOTCH */
  }
#ifdef SCOTCH_DEBUG_DGRAPH2                       /* Prevent Valgrind from yelling */
  else {                                          /* Process is not root           */
    verttax =
    velotax =
    vlbltax =
    edgetax =
    edlotax = NULL;
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (dgrfptr->vendloctax == (dgrfptr->vertloctax + 1)) { /* If distributed graph is compact                                                */
    if (dgraphGatherAll3 (dgrfptr->vertloctax + baseval + 1, dgrfptr->vertlocnbr, /* Do not send first index, it is always equal to baseval */
                          verttax + 1,            /* First index will always be equal to baseval too, and procdsptab holds based values     */
                          dgrfptr->proccnttab, dgrfptr->procdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (3)");
      return     (1);
    }

    if (cgrfptr != NULL) {
      Gnum                procnum;

      verttax[baseval] = baseval;
      for (procnum = 1; procnum < dgrfptr->procglbnbr; procnum ++) { /* Adjust index sub-arrays for all processes except the first one */
        Gnum                vertnum;
        Gnum                vertnnd;
        Gnum                edgedlt;

        for (vertnum = dgrfptr->procdsptab[procnum] + 1,
             vertnnd = dgrfptr->proccnttab[procnum] + vertnum,
             edgedlt = verttax[vertnum - 1] - baseval;
             vertnum < vertnnd; vertnum ++)
          verttax[vertnum] += edgedlt;
      }
    }
  }
  else {                                          /* Distributed graph is not compact */
    Gnum                vertlocnum;
    Gnum * restrict     edgelocptr;

    vertloctax -= baseval;                        /* Base temporary vertex array   */
    for (vertlocnum = baseval, edgelocptr = edgeloctab; /* Build vertex send array */
         vertlocnum < dgrfptr->vertlocnnd; vertlocnum ++) {
      Gnum                edgelocnum;

      vertloctax[vertlocnum] = dgrfptr->vendloctax[vertlocnum] - dgrfptr->vertloctax[vertlocnum]; /* Get edge counts */

      for (edgelocnum = dgrfptr->vertloctax[vertlocnum]; edgelocnum < dgrfptr->vendloctax[vertlocnum]; edgelocnum ++)
        *edgelocptr ++ = edgeloctax[edgelocnum];
    }

    if (dgraphGatherAll3 (vertloctax + baseval, dgrfptr->vertlocnbr,
                          verttax + 1,            /* First index will always be equal to baseval, and procdsptab holds based values */
                          dgrfptr->proccnttab, dgrfptr->procdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (4)");
      return     (1);
    }

    if (cgrfptr != NULL) {
      Gnum                vertnum;
      Gnum                edgenum;

      verttax[baseval] = baseval;
      for (vertnum = baseval + 1, edgenum = baseval; /* Create compact centralized vertex array */
           vertnum <= cgrfptr->vertnnd; vertnum ++) {
        edgenum += verttax[vertnum];
        verttax[vertnum] = edgenum;
      }
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (verttax[cgrfptr->vertnnd] != (cgrfptr->edgenbr + baseval)) {
        errorPrint ("dgraphGatherAll2: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    }
  }

  if (dgrfptr->veloloctax != NULL) {
    if (dgraphGatherAll3 (dgrfptr->veloloctax + baseval, dgrfptr->vertlocnbr,
                          velotax,                /* Based array since procdsptab holds based values */
                          dgrfptr->proccnttab, dgrfptr->procdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (6)");
      return     (1);
    }
  }
  if (dgrfptr->vnumloctax != NULL) {
    if (dgraphGatherAll3 (dgrfptr->vnumloctax + baseval, dgrfptr->vertlocnbr,
                          vnumtax,                /* Based array since procdsptab holds based values */
                          dgrfptr->proccnttab, dgrfptr->procdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (5)");
      return     (1);
    }
  }
  if (dgrfptr->vlblloctax != NULL) {
    if (dgraphGatherAll3 (dgrfptr->vlblloctax + baseval, dgrfptr->vertlocnbr,
                          vlbltax,                /* Based array since procdsptab holds based values */
                          dgrfptr->proccnttab, dgrfptr->procdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (7)");
      return     (1);
    }
  }

  if (cgrfptr != NULL) {
    Gnum                procnum;
    Gnum                edgenum;

    for (procnum = 0, edgenum = baseval;          /* Build arrays for MPI_Gatherv on edge arrays */
         procnum < dgrfptr->procglbnbr; procnum ++) {
      recvcnttab[procnum] = verttax[dgrfptr->procdsptab[procnum] + dgrfptr->proccnttab[procnum]] -
                            verttax[dgrfptr->procdsptab[procnum]]; /* verttax used twice since centralized graph is compact */
      recvdsptab[procnum] = edgenum;
      edgenum += recvcnttab[procnum];
    }
#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((recvdsptab[dgrfptr->procglbnbr - 1] + recvcnttab[dgrfptr->procglbnbr - 1]) != (cgrfptr->edgenbr + baseval)) {
      errorPrint ("dgraphGatherAll2: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  }

  if (dgrfptr->vendloctax == (dgrfptr->vertloctax + 1)) { /* If distributed graph is compact         */
    if (dgraphGatherAll3 (dgrfptr->edgeloctax + baseval, dgrfptr->edgelocnbr, /* Send global indices */
                          edgetax,                /* Based array as recvdsptab holds based values    */
                          recvcnttab, recvdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (8)");
      return     (1);
    }

    if (dgrfptr->edloloctax != NULL) {
      if (dgraphGatherAll3 (dgrfptr->edloloctax + baseval, dgrfptr->edgelocnbr,
                            edlotax,              /* Based array as recvdsptab holds based values */
                            recvcnttab, recvdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphGatherAll2: communication error (9)");
        return     (1);
      }
    }
  }
  else {                                          /* Distributed graph is not compact */
    if (dgraphGatherAll3 (edgeloctab, dgrfptr->edgelocnbr,
                          edgetax,                /* Based array as recvdsptab holds based values */
                          recvcnttab, recvdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphGatherAll2: communication error (10)");
      return     (1);
    }

    if (dgrfptr->edloloctax != NULL) {
      Gnum                vertlocnum;
      Gnum * restrict     edlolocptr;

      for (vertlocnum = baseval, edlolocptr = edgeloctab; /* Recycle edge send array to build edge load send array */
           vertlocnum < dgrfptr->vertlocnnd; vertlocnum ++) {
        Gnum                edgelocnum;

        for (edgelocnum = dgrfptr->vertloctax[vertlocnum]; edgelocnum < dgrfptr->vendloctax[vertlocnum]; edgelocnum ++)
          *edlolocptr ++ = dgrfptr->edloloctax[edgelocnum];
      }

      if (dgraphGatherAll3 (edgeloctab, dgrfptr->edgelocnbr, /* Send compacted edge load array    */
                            edlotax,              /* Based array as recvdsptab holds based values */
                            recvcnttab, recvdsptab, protnum, dgrfptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphGatherAll2: communication error (11)");
        return     (1);
      }
    }
  }

  if (cgrfptr != NULL) {
    if ((dgrfptr->procdsptab[dgrfptr->procglbnbr] != /* If graph has holes, relabel end vertices */
         dgrfptr->procvrttab[dgrfptr->procglbnbr])) {
      Gnum                procnum;

      for (procnum = 0; procnum < dgrfptr->procglbnbr; procnum ++) { /* Accelerate search per sender process */
        Gnum                vertlocmin;
        Gnum                vertlocmax;
        Gnum                vertlocadj;
        Gnum                edgelocnum;
        Gnum                edgelocnnd;

        vertlocmin = dgrfptr->procvrttab[procnum];  /* Initialize search accelerator */
        vertlocmax = dgrfptr->procvrttab[procnum + 1];
        vertlocadj = dgrfptr->procdsptab[procnum] - vertlocmin;

        for (edgelocnum = recvdsptab[procnum], edgelocnnd = edgelocnum + recvcnttab[procnum];
             edgelocnum < edgelocnnd; edgelocnum ++) {
          Gnum                vertlocend;

          vertlocend = cgrfptr->edgetax[edgelocnum];

          if ((vertlocend >= vertlocmin) &&       /* If end vertex is local with respect to current process */
              (vertlocend <  vertlocmax))
            cgrfptr->edgetax[edgelocnum] = vertlocend + vertlocadj;
          else {                                  /* End vertex is not local */
            int                 procngbmin;
            int                 procngbmax;

            for (procngbmin = 0, procngbmax = dgrfptr->procglbnbr;
                 procngbmax - procngbmin > 1; ) {
              int               procngbnum;

              procngbnum = (procngbmax + procngbmin) / 2;
              if (dgrfptr->procvrttab[procngbnum] <= vertlocend)
                procngbmin = procngbnum;
              else
                procngbmax = procngbnum;
            }
            cgrfptr->edgetax[edgelocnum] = vertlocend + dgrfptr->procdsptab[procngbmin] - dgrfptr->procvrttab[procngbmin];
          }
        }
      }
    }

    if (cgrfptr->edlotax == NULL)                 /* If no edge loads         */
      cgrfptr->edlosum = cgrfptr->edgenbr;        /* Edge load sum is trivial */
    else {
      if (edlosum >= 0)                           /* If edge load sum already computed by library call */
        cgrfptr->edlosum = edlosum;
      else {                                      /* Compute it from scratch on every root process (small graph assumed) */
        Gnum                edgenum;
        Gnum                edgennd;
        Gnum                edlotmp;

        for (edgenum = cgrfptr->baseval, edgennd = edgenum + cgrfptr->edgenbr, edlotmp = 0; /* Edge load array is always compact */
             edgenum < edgennd; edgenum ++)
          edlotmp += cgrfptr->edlotax[edgenum];

        cgrfptr->edlosum = edlotmp;
      }
    }
  }

  memFree (recvcnttab);

#ifdef SCOTCH_DEBUG_DGRAPH2
  cheklocval = (cgrfptr != NULL) ? graphCheck (cgrfptr) : 0;
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, dgrfptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphGatherAll2: communication error (12)");
    return     (1);
  }
  if (chekglbval != 0) {
    errorPrint ("dgraphGatherAll2: inconsistent centralized graph data");
    if (cgrfptr != NULL)
      graphFree  (cgrfptr);
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}

/* This function gathers on all processes
** the pieces of a distributed graph to
** build a centralized graph.
** Since the resulting centralized graphs are
** supposed to be small in the general case,
** edlosum is computed without communication
** on each of the processors.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
dgraphGatherAll (
const Dgraph * restrict const dgrfptr,            /* Distributed graph */
Graph * restrict              cgrfptr)            /* Centralized graph */
{
  return (dgraphGatherAll2 (dgrfptr, cgrfptr, -1, -1));
}
