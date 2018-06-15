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
/**   NAME       : dgraph_scatter.c                        **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine that   **/
/**                builds a distributed graph by evenly    **/
/**                distributing the pieces of a central-   **/
/**                ized graph across processors.           **/
/**                                                        **/
/**                # Version P0.0 : from : 01 apr 1997     **/
/**                                 to     20 jun 1997     **/
/**                # Version P0.1 : from : 14 apr 1998     **/
/**                                 to     20 jun 1998     **/
/**                # Version P0.2 : from : 19 may 1999     **/
/**                                 to     19 may 1999     **/
/**                # Version 5.0  : from : 27 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                                                        **/
/**   NOTES      : # The definitions of MPI_Scatter and    **/
/**                  MPI_Scatterv indicate that elements   **/
/**                  in the send array should not be read  **/
/**                  more than once. Great care should be  **/
/**                  taken to enforce this rule, especial- **/
/**                  ly when the number of vertices in the **/
/**                  centralized graph is smaller than the **/
/**                  number of processors.                 **/
/**                                                        **/
/**                # When the source graph is not compact, **/
/**                  compacted arrays are created prior to **/
/**                  sending parts of them. In a really    **/
/**                  efficient implementation, these       **/
/**                  should be created by pieces and sent  **/
/**                  in a one-to-one way so as to save as  **/
/**                  much memory as possible. This is yet  **/
/**                  to be done.                           **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_SCATTER

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph.h"

/* Service function which creates compact
** arrays from non-compact ones.
*/

static
void
dgraphScatter2 (
const Graph * restrict const  cgrfptr,
Gnum * restrict               verttax,
Gnum * restrict               edgetax,
Gnum * restrict               edlotax)
{
  Gnum                vertnum;                    /* Current vertex number in compacted arrays */
  Gnum                edgenum;                    /* Current edge number in compacted arrays   */

  for (vertnum = edgenum = cgrfptr->baseval;
       vertnum < cgrfptr->vertnnd; vertnum ++) {
    Gnum                edgetmp;
    Gnum                edgetnd;

    verttax[vertnum] = edgenum;
    edgetmp = cgrfptr->verttax[vertnum];
    edgetnd = cgrfptr->vendtax[vertnum];

    for ( ; edgetmp < edgetnd; edgetmp ++, edgenum ++)
      edgetax[edgenum] = cgrfptr->edgetax[edgetmp];

    if (edlotax != NULL) {
      for (edgetmp = cgrfptr->verttax[vertnum], edgenum = verttax[vertnum];
           edgetmp < edgetnd; edgetmp ++, edgenum ++)
        edlotax[edgenum] = cgrfptr->edlotax[edgetmp];
    }
  }
  verttax[vertnum] = edgenum;
}

/* This function evenly distributes the pieces
** of a centralized graph across processors.
** It returns:
** - 0   : if scattering has succeeded.
** - !0  : on error.
*/

int
dgraphScatter (
Dgraph * restrict const       grafptr,            /* Distributed graph            */
const Graph * restrict const  cgrfptr)            /* Centralized graph to scatter */
{
  Gnum                baseval;                    /* Base value                                           */
  Gnum * restrict     verttax;                    /* Array of vertices when edge array is not compact     */
  Gnum * restrict     edgetax;                    /* Compact array of edges when edge aray is not compact */
  Gnum * restrict     edlotax;                    /* Compact array of edges weights                       */
  Gnum                vertlocnum;                 /* Current local vertex number                          */
  Gnum                vertlocnbr;                 /* Number of local vertices                             */
  Gnum                vertlocnnd;
  Gnum * restrict     vertloctax;                 /* Array of local vertices                              */
  Gnum * restrict     veloloctax;                 /* Array of local vertex weights                        */
  Gnum                velolocnbr;
  Gnum                vlbllocnbr;
  Gnum * restrict     vlblloctax;                 /* Array of local vertex labels                         */
  Gnum                edgelocnbr;                 /* Number of local edges                                */
  Gnum                edlolocnbr;
  Gnum * restrict     edgeloctax;                 /* Array of local edges                                 */
  Gnum * restrict     edloloctax;                 /* Array of local edge weights                          */
  int * restrict      attrdsptab;                 /* Displacement array for scatter operations            */
  int * restrict      attrcnttab;                 /* Count array for scatter operations                   */
  Gnum * restrict     attrdattab;                 /* Temporary array to avoid multiple scatter reads      */
  Gnum                reduloctab[9];              /* Arrays for reductions                                */
  Gnum                reduglbtab[9];
  Gnum                vertlocadj;                 /* Local vertex array adjust                            */
  int                 protnum;                    /* Root process                                         */

  if (cgrfptr != NULL) {                          /* If centralized graph provided    */
    if (cgrfptr->vendtax != (cgrfptr->verttax + 1)) { /* If edge array is not compact */
      Gnum                edlonbr;

      edlonbr = (cgrfptr->edlotax != NULL) ? cgrfptr->edgenbr : 0;
      if (memAllocGroup ((void **) (void *)
                         &verttax, (size_t) ((cgrfptr->vertnbr + 1) * sizeof (Gnum)),
                         &edgetax, (size_t) (cgrfptr->edgenbr       * sizeof (Gnum)),
                         &edlotax, (size_t) (edlonbr                * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("dgraphScatter: out of memory (1)");
        return     (1);
      }
      verttax -= cgrfptr->baseval;
      edgetax -= cgrfptr->baseval;
      edlotax  = (cgrfptr->edlotax != NULL) ? (edlotax - cgrfptr->baseval) : NULL;
      dgraphScatter2 (cgrfptr, verttax, edgetax, edlotax);
    }
    else {
      verttax = cgrfptr->verttax;
      edgetax = cgrfptr->edgetax;
      edlotax = cgrfptr->edlotax;
    }
    reduloctab[0] = 1;                            /* This process is the root */
    reduloctab[1] = (Gnum) grafptr->proclocnum;   /* Get its number           */
    reduloctab[2] = cgrfptr->baseval;
    reduloctab[3] = cgrfptr->vertnbr;
    reduloctab[4] = cgrfptr->edgenbr;
    reduloctab[5] = cgrfptr->velosum;
    reduloctab[6] = (cgrfptr->velotax != NULL) ? 1 : 0;
    reduloctab[7] = (cgrfptr->vlbltax != NULL) ? 1 : 0;
    reduloctab[8] = (cgrfptr->edlotax != NULL) ? 1 : 0;
  }
  else {
    reduloctab[0] =                               /* This process is not the root */
    reduloctab[1] =
    reduloctab[2] =
    reduloctab[3] =
    reduloctab[4] =
    reduloctab[5] =
    reduloctab[6] =
    reduloctab[7] =
    reduloctab[8] = 0;
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 9, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphScatter: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("dgraphScatter: should have only one root");
    return     (1);
  }

  baseval    = reduglbtab[2];
  vertlocnbr = DATASIZE (reduglbtab[3], grafptr->procglbnbr, grafptr->proclocnum);
  velolocnbr = (reduglbtab[6] != 0) ? vertlocnbr : 0;
  vlbllocnbr = (reduglbtab[7] != 0) ? vertlocnbr : 0;
  if (memAllocGroup ((void **) (void *)
                     &vertloctax, (size_t) ((vertlocnbr + 1) * sizeof (Gnum)),
                     &veloloctax, (size_t) (velolocnbr       * sizeof (Gnum)),
                     &vlblloctax, (size_t) (vlbllocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphScatter: out of memory (2)");
    if ((cgrfptr != NULL) && (cgrfptr->verttax != verttax))
      memFree (verttax + baseval);                /* Free group leader */
    return (1);
  }
  vertloctax -= baseval;
  veloloctax  = (reduglbtab[6] != 0) ? (veloloctax - baseval) : NULL;
  vlblloctax  = (reduglbtab[7] != 0) ? (vlblloctax - baseval) : NULL;

  protnum = (int) reduglbtab[1];
  if (cgrfptr != NULL) {                          /* If root process */
    Gnum                procnum;

    if (memAllocGroup ((void **) (void *)
                       &attrdattab, (size_t) (grafptr->procglbnbr * sizeof (Gnum)),
                       &attrdsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                       &attrcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("dgraphScatter: out of memory (3)");
      memFree (vertloctax + baseval);
      if (cgrfptr->verttax != verttax)
        memFree (verttax + baseval);              /* Free group leader */
      return (1);
    }

    attrdsptab[0] = 0;                            /* Build arrays for MPI_Scatterv */
    attrcnttab[0] = DATASIZE (reduglbtab[3], grafptr->procglbnbr, 0);
    attrdattab[0] = verttax[attrdsptab[0] + attrcnttab[0] + baseval];
    for (procnum = 1; procnum < grafptr->procglbnbr; procnum ++) {
      attrdsptab[procnum] = attrdsptab[procnum - 1] + attrcnttab[procnum - 1];
      attrcnttab[procnum] = DATASIZE (reduglbtab[3], grafptr->procglbnbr, procnum);
      attrdattab[procnum] = verttax[attrdsptab[procnum] + attrcnttab[procnum] + baseval];
    }

    if (MPI_Scatterv (verttax + baseval, attrcnttab, attrdsptab, GNUM_MPI, /* Perform two scatters since cannot avoid multiple reads with only one scatter */
                      vertloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (2)");
      return     (1);
    }
    if (MPI_Scatter  (attrdattab, 1, GNUM_MPI,
                      vertloctax + baseval + vertlocnbr, 1, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (3)");
      return     (1);
    }
    if (reduglbtab[6] != 0) {                     /* Scatter vertex loads */
      if (MPI_Scatterv (cgrfptr->velotax + baseval, attrcnttab, attrdsptab, GNUM_MPI,
                        veloloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (4)");
        return     (1);
      }
    }
    if (reduglbtab[7] != 0) {                     /* Scatter labels */
      if (MPI_Scatterv (cgrfptr->vlbltax + baseval, attrcnttab, attrdsptab, GNUM_MPI,
                        vlblloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (5)");
        return     (1);
      }
    }
  }
  else {                                          /* Process is not root */
    if (MPI_Scatterv (NULL, NULL, NULL, GNUM_MPI,
                      vertloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (6)");
      return     (1);
    }
    if (MPI_Scatter  (NULL, 1, GNUM_MPI, vertloctax + baseval + vertlocnbr, 1, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (7)");
      return     (1);
    }
    if (reduglbtab[6] != 0) {                     /* Scatter vertex loads */
      if (MPI_Scatterv (NULL, NULL, NULL, GNUM_MPI,
                        veloloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (8)");
        return     (1);
      }
    }
    if (reduglbtab[7] != 0) {                     /* Scatter labels */
      if (MPI_Scatterv (NULL, NULL, NULL, GNUM_MPI,
                        vlblloctax + baseval, vertlocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (9)");
        return     (1);
      }
    }
  }

  vertlocadj = vertloctax[baseval] - baseval;     /* Compute local indices */
  for (vertlocnum = baseval, vertlocnnd = vertlocnbr + baseval;
       vertlocnum <= vertlocnnd; vertlocnum ++)
     vertloctax[vertlocnum] -= vertlocadj;

  edgelocnbr = vertloctax[vertlocnnd] - vertloctax[baseval];
  edlolocnbr = (reduglbtab[8] != 0) ? edgelocnbr : 0;
  if (memAllocGroup ((void **) (void *)
                     &edgeloctax, (size_t) (edgelocnbr * sizeof (Gnum)),
                     &edloloctax, (size_t) (edlolocnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphScatter: out of memory (4)");
    if (cgrfptr != NULL) {
      memFree (attrdattab);                       /* Free group leader */
      if (cgrfptr->verttax != verttax)
        memFree (verttax + baseval);              /* Free group leader */
    }
    memFree (vertloctax + baseval);
    return  (1);
  }
  edgeloctax -= baseval;
  edloloctax  = (reduglbtab[8] != 0) ? edloloctax - baseval : NULL;

  if (cgrfptr != NULL) {                          /* If root process */
    Gnum                procnum;

    for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) { /* Build arrays for MPI_Scatterv */
      attrcnttab[procnum] = verttax[attrdsptab[procnum] + attrcnttab[procnum]+baseval] -
                            verttax[attrdsptab[procnum] + baseval];
      attrdsptab[procnum] = verttax[attrdsptab[procnum] + baseval] - baseval;
    }

    if (MPI_Scatterv (edgetax + baseval, attrcnttab, attrdsptab, GNUM_MPI,
                      edgeloctax + baseval, edgelocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (10)");
      return     (1);
    }
    if (reduglbtab[8] != 0) {
      if (MPI_Scatterv (edlotax + baseval, attrcnttab, attrdsptab, GNUM_MPI,
                        edloloctax + baseval, edgelocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (11)");
        return     (1);
      }
    }
     memFree (attrdattab);                        /* Free group leader */
     if (cgrfptr->verttax != verttax)
       memFree (verttax + baseval);
  }
  else {                                          /* Process is not root */
    if (MPI_Scatterv (NULL, NULL, NULL, GNUM_MPI,
                      edgeloctax + baseval , edgelocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphScatter: communication error (12)");
      return     (1);
    }
    if (reduglbtab[8] != 0) {
      if (MPI_Scatterv (NULL, NULL, NULL, GNUM_MPI,
                        edloloctax + baseval , edgelocnbr, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dgraphScatter: communication error (13)");
        return     (1);
      }
    }
  }

  if (dgraphBuild (grafptr, baseval,
                   vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, veloloctax, NULL, NULL,
                   edgelocnbr, edgelocnbr, edgeloctax, NULL, edloloctax) != 0) {
    memFree (edgeloctax + baseval);
    memFree (vertloctax + baseval);
    return  (1);
  }

  grafptr->flagval   |= DGRAPHFREETABS | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* Give ownership of arrays to graph */
  grafptr->vlblloctax = vlblloctax;               /* Add labels afterwards, since relabeling already done        */

  return (0);
}
