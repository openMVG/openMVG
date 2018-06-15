/* Copyright 2007,2008,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_check.c                         **/
/**                                                        **/
/**   AUTHORS    : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                bipartition graph consistency checking  **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 10 sep 2007     **/
/**                                 to     22 jul 2008     **/
/**                # Version 6.0  : from : 03 sep 2011     **/
/**                                 to     31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "arch.h"
#include "dgraph.h"
#include "bdgraph.h"

int
bdgraphCheck (
const Bdgraph * restrict const grafptr)
{
  Dgraph                    grafdat;              /* Dummy graph for ghost edge array */
  MPI_Comm                  proccomm;             /* Graph communicator               */
  int * restrict            flagloctax;           /* Frontier flag array              */
  GraphPart * restrict      partgsttax;
  Gnum                      fronlocnum;
  Gnum                      vertlocnum;           /* Number of current vertex         */
  Gnum                      complocload[2];
  Gnum                      complocsize[2];
  Gnum                      commcut[2];
  Gnum                      commlocloadintn;
  Gnum                      commlocloadextn;
  Gnum                      commlocgainextn;
  Gnum                      edlolocval;
  Gnum                      reduloctab[21];       /* Arrays for reductions            */
  Gnum                      reduglbtab[21];
  int                       chekglbval;           /* Global consistency flag          */
  int                       cheklocval;           /* Local consistency flag           */

  proccomm = grafptr->s.proccomm;
  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize */
    errorPrint ("bdgraphCheck: communication error (1)");
    return     (1);
  }

  cheklocval = 0;                                 /* Assume everything is all right */
  if ((grafptr->compglbload0min < 0) ||
      (grafptr->compglbload0max > grafptr->s.veloglbsum)) {
    errorPrint ("bdgraphCheck: invalid extrema loads");
    cheklocval = 1;
  }

  if (grafptr->compglbload0 != (grafptr->compglbload0avg + grafptr->compglbload0dlt)) {
    errorPrint ("bdgraphCheck: invalid global balance");
    cheklocval |= 2;
  }

  if ((grafptr->fronlocnbr < 0) ||
      (grafptr->fronlocnbr > grafptr->s.vertlocnbr)) {
    errorPrint ("bdgraphCheck: invalid number of local frontier vertices");
    cheklocval |= 4;
  }

  if (grafptr->partgsttax != NULL) {
    for (vertlocnum = grafptr->s.baseval; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
      if ((grafptr->partgsttax[vertlocnum] | 1) != 1) { /* If part is neither 0 nor 1 */
        errorPrint ("bdgraphCheck: invalid local part array");
        cheklocval |= 8;
        break;
      }
    }
  }

  grafdat = grafptr->s;                           /* Copy minimal distributed graph data      */
  if (dgraphGhst (&grafdat) != 0) {               /* Create ghost edge array if did not exist */
    errorPrint ("bdgraphCheck: cannot compute ghost edge array");
    cheklocval |= 16;
  }

  if (memAllocGroup ((void **) (void *)
                     &partgsttax, (size_t) (grafdat.vertgstnbr    * sizeof (GraphPart)),
                     &flagloctax, (size_t) (grafptr->s.vertlocnbr * sizeof (int)), NULL) == NULL) {
    errorPrint ("bdgraphCheck: out of memory");
    cheklocval |= 32;
  }
  else {
    memSet (flagloctax, ~0, grafptr->s.vertlocnbr * sizeof (int));
    flagloctax -= grafptr->s.baseval;

    for (fronlocnum = 0; fronlocnum < grafptr->fronlocnbr; fronlocnum ++) {
      Gnum                vertlocnum;

      vertlocnum = grafptr->fronloctab[fronlocnum];
      if ((vertlocnum < grafptr->s.baseval) || (vertlocnum >= grafptr->s.vertlocnnd)) {
        errorPrint ("bdgraphCheck: invalid vertex index in frontier array");
        cheklocval |= 64;
        break;
      }
      if (flagloctax[vertlocnum] != ~0) {
        errorPrint ("bdgraphCheck: duplicate vertex in frontier array");
        cheklocval |= 128;
        break;
      }
      flagloctax[vertlocnum] = 0;
    }
  }

  reduloctab[0]   =   grafptr->commglbload;
  reduloctab[1]   = - grafptr->commglbload;
  reduloctab[2]   =   grafptr->compglbload0;
  reduloctab[3]   = - grafptr->compglbload0;
  reduloctab[4]   =   grafptr->s.veloglbsum - grafptr->compglbload0;
  reduloctab[5]   = - grafptr->s.veloglbsum + grafptr->compglbload0;
  reduloctab[6]   =   grafptr->compglbsize0;
  reduloctab[7]   = - grafptr->compglbsize0;
  reduloctab[8]   =   grafptr->s.vertglbnbr - grafptr->compglbsize0;
  reduloctab[9]   = - grafptr->s.vertglbnbr + grafptr->compglbsize0;
  reduloctab[10]  =   grafptr->commglbgainextn;
  reduloctab[11]  = - grafptr->commglbgainextn;
  reduloctab[12]  =   grafptr->commglbgainextn0;
  reduloctab[13]  = - grafptr->commglbgainextn0;
  reduloctab[14]  =   grafptr->commglbloadextn0;
  reduloctab[15]  = - grafptr->commglbloadextn0;
  reduloctab[16]  =   grafptr->fronglbnbr;
  reduloctab[17]  = - grafptr->fronglbnbr;
  reduloctab[18]  =   grafptr->levlnum;
  reduloctab[19]  = - grafptr->levlnum;
  reduloctab[20]  =   cheklocval;

  if (MPI_Allreduce (reduloctab, reduglbtab, 21, GNUM_MPI, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphCheck: communication error (2)");
    return     (1);
  }

  if (reduglbtab[20] != 0) {                      /* Exit and Return information of previous errors */
    if (partgsttax != NULL)
      memFree (partgsttax);                       /* Free yet unbased group leader */
    return ((int) reduglbtab[20]);
  }

  if ((reduglbtab[1]  != - reduglbtab[0])  ||
      (reduglbtab[3]  != - reduglbtab[2])  ||
      (reduglbtab[5]  != - reduglbtab[4])  ||
      (reduglbtab[7]  != - reduglbtab[6])  ||
      (reduglbtab[9]  != - reduglbtab[8])  ||
      (reduglbtab[11] != - reduglbtab[10]) ||
      (reduglbtab[13] != - reduglbtab[12]) ||
      (reduglbtab[15] != - reduglbtab[14]) ||
      (reduglbtab[17] != - reduglbtab[16]) ||
      (reduglbtab[19] != - reduglbtab[18])) {
    errorPrint ("bdgraphCheck: inconsistent global graph data");
    return     (1);
  }

  if (grafptr->partgsttax != NULL)
    memCpy (partgsttax, grafptr->partgsttax + grafptr->s.baseval, grafptr->s.vertlocnbr * sizeof (GraphPart)); /* Copy local part data */
  else
    memSet (partgsttax, 0, grafptr->s.vertlocnbr * sizeof (GraphPart));
  dgraphHaloSync (&grafdat, partgsttax, GRAPHPART_MPI); /* Spread yet unbased halo part data across neighboring processes */
  partgsttax -= grafptr->s.baseval;
  cheklocval  = 0;

  for (fronlocnum = 0; fronlocnum < grafptr->fronlocnbr; fronlocnum ++) {
    Gnum                vertlocnum;
    Gnum                edgelocnum;
    GraphPart           commcut;
    GraphPart           partval;

    vertlocnum = grafptr->fronloctab[fronlocnum];
    partval    = partgsttax[vertlocnum];

    for (edgelocnum = grafptr->s.vertloctax[vertlocnum], commcut = 0;
         edgelocnum < grafptr->s.vendloctax[vertlocnum]; edgelocnum ++) {
      GraphPart           partdlt;

      partdlt  = partgsttax[grafdat.edgegsttax[edgelocnum]] ^ partval;
      commcut |= partdlt;
    }
    if (commcut == 0) {
      errorPrint ("bdgraphCheck: invalid vertex in frontier array");
      cheklocval |= 1;
      break;
    }
  }

  complocload[0]   =
  complocload[1]   = 0;
  complocsize[0]   =
  complocsize[1]   = 0;
  commlocloadintn  = 0;
  commlocloadextn  = 0;
  commlocgainextn  = 0;
  edlolocval       = 1;                            /* Assume edges are not weighted */
  for (vertlocnum = grafptr->s.baseval; vertlocnum < grafptr->s.vertlocnnd; vertlocnum ++) {
    Gnum                partval;                  /* Part of current vertex */
    Gnum                edgelocnum;               /* Number of current edge */

    partval = (Gnum) partgsttax[vertlocnum];
    if (grafptr->veexloctax != NULL) {
      Gnum                veexval;

      veexval = grafptr->veexloctax[vertlocnum];
      commlocloadextn += veexval * partval;
      commlocgainextn += veexval * (1 - 2 * partval);
    }

    complocload[partval] += (grafptr->s.veloloctax == NULL) ? 1 : grafptr->s.veloloctax[vertlocnum];
    complocsize[partval] ++;

    commcut[0] =
    commcut[1] = 0;
    for (edgelocnum = grafptr->s.vertloctax[vertlocnum]; edgelocnum < grafptr->s.vendloctax[vertlocnum]; edgelocnum ++) {
      int                 partend;
      int                 partdlt;

      if (grafptr->s.edloloctax != NULL)
        edlolocval = grafptr->s.edloloctax[edgelocnum];
      partend = partgsttax[grafdat.edgegsttax[edgelocnum]];
      partdlt = partval ^ partend;
      commcut[partend] ++;
      commlocloadintn += partdlt * edlolocval;    /* Internal load is accounted for twice */
    }

    if ((commcut[0] != 0) && (commcut[1] != 0) && /* If vertex should be in separator */
        (flagloctax[vertlocnum] != 0)) {
      errorPrint ("bdgraphCheck: vertex should be in separator");
      cheklocval |= 2;
    }
  }

  if (grafptr->s.edgegsttax != grafdat.edgegsttax) /* If ghost edge array was allocated here, free it manually */
    memFree (grafdat.edgegsttax + grafptr->s.baseval);
  if (grafptr->s.procsidtab != grafdat.procsidtab) /* The same for procsidtab */
    memFree (grafdat.procsidtab);
  memFree (partgsttax + grafptr->s.baseval);      /* Free group leader */

  if ((cheklocval == 0) &&
      ((complocsize[0] != grafptr->complocsize0) ||
       (complocsize[1] != (grafptr->s.vertlocnbr - grafptr->complocsize0)))) {
    errorPrint ("bdgraphCheck: invalid local part size");
    cheklocval |= 4;
  }

  reduloctab[0] = complocload[0];
  reduloctab[1] = complocsize[0];
  reduloctab[2] = commlocloadintn;                /* Twice the internal load; sum globally before dividing by two */
  reduloctab[3] = commlocloadextn;
  reduloctab[4] = commlocgainextn;
  reduloctab[5] = cheklocval;

  if (MPI_Allreduce (reduloctab, reduglbtab, 6, GNUM_MPI, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphCheck: communication error (3)");
    return     (1);
  }

  if (reduglbtab[5] != 0)                         /* Return from previous errors */
    return (1);

  if (grafptr->compglbload0 != reduglbtab[0]) {
    errorPrint ("bdgraphCheck: invalid global part loads");
    cheklocval |= 8;
  }

  if (grafptr->compglbsize0 != reduglbtab[1]) {
    errorPrint ("bdgraphCheck: invalid global part sizes");
    cheklocval |= 16;
  }

  if (grafptr->commglbload != ((reduglbtab[2] / 2) * grafptr->domndist + reduglbtab[3] + grafptr->commglbloadextn0)) {
    errorPrint ("bdgraphCheck: invalid global communication loads");
    cheklocval |= 32;
  }
 
  if (grafptr->commglbgainextn != reduglbtab[4]) {
    errorPrint ("bdgraphCheck: invalid global communication gains");
    cheklocval |= 64;
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphCheck: communication error (4)");
    return     (1);
  }

  return (chekglbval);
}
