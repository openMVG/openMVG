/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph_check.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                separator graph consistency checking    **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 07 feb 2006     **/
/**                                 to     01 mar 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "vdgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* separator graphs.     */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given distributed separator graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
vdgraphCheck (
const Vdgraph * const       grafptr)
{
  Dgraph                grafdat;                  /* Dummy graph for ghost edge array */
  MPI_Comm              proccomm;                 /* Graph communicator               */
  Gnum                  vertnum;                  /* Number of current vertex         */
  Gnum                  fronnum;                  /* Number of frontier vertex        */
  Gnum                  complocload[3];
  Gnum                  complocsize[3];
  Gnum                  commcut[3];
  GraphPart * restrict  partgsttax;
  Gnum                  reduloctab[11];           /* Arrays for reductions            */
  Gnum                  reduglbtab[11];
  int                   cheklocval;               /* Local consistency flag           */
  int                   chekglbval;               /* Global consistency flag          */

  proccomm = grafptr->s.proccomm;
  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize */
    errorPrint ("vdgraphCheck: communication error (1)");
    return     (1);
  }

  cheklocval = 0;                                 /* Assume everything is all right */

  if ((grafptr->compglbload[0]  + grafptr->compglbload[1] + grafptr->compglbload[2]) != grafptr->s.veloglbsum) {
    errorPrint ("vdgraphCheck: invalid global load sum");
    cheklocval = 1;
  }

  if (grafptr->compglbloaddlt != (grafptr->compglbload[0] - grafptr->compglbload[1])) {
    errorPrint ("vdgraphCheck: invalid global balance");
    cheklocval |= 2;
  }

  if ((grafptr->compglbsize[0] + grafptr->compglbsize[1] + grafptr->compglbsize[2]) != grafptr->s.vertglbnbr) {
    errorPrint ("vdgraphCheck: invalid global size sum");
    cheklocval |= 4;
  }

  if ((grafptr->complocsize[0] + grafptr->complocsize[1] + grafptr->complocsize[2]) != grafptr->s.vertlocnbr) {
    errorPrint ("vdgraphCheck: invalid local size sum");
    cheklocval |= 8;
  }

  if ((grafptr->complocsize[2] < 0) ||
      (grafptr->complocsize[2] > grafptr->s.vertlocnbr)) {
    errorPrint ("vdgraphCheck: invalid number of local frontier vertices");
    cheklocval |= 16;
  }

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertlocnnd; vertnum ++) {
    if (grafptr->partgsttax[vertnum] > 2) {
      errorPrint ("vdgraphCheck: invalid local part array");
      cheklocval |= 32;
      break;
    }
  }

  for (fronnum = 0; fronnum < grafptr->complocsize[2]; fronnum ++) {
    Gnum                vertnum;

    vertnum = grafptr->fronloctab[fronnum];
    if ((vertnum < grafptr->s.baseval) || (vertnum >= grafptr->s.vertlocnnd)) {
      errorPrint ("vdgraphCheck: invalid vertex index in frontier array");
      cheklocval |= 64;
      break;
    }
    if (grafptr->partgsttax[vertnum] != 2) {
      errorPrint ("vdgraphCheck: invalid vertex in frontier array");
      cheklocval |= 64;
      break;
    }
  }

  grafdat = grafptr->s;                           /* Copy minimal distributed graph data      */
  if (dgraphGhst (&grafdat) != 0) {               /* Create ghost edge array if did not exist */
    errorPrint ("vdgraphCheck: cannot compute ghost edge array");
    cheklocval |= 128;
  }

  if ((partgsttax = memAlloc (grafdat.vertgstnbr * sizeof (byte))) == NULL) {
    errorPrint ("vdgraphCheck: out of memory");
    cheklocval |= 256;
  }

  reduloctab[0]  =   grafptr->compglbload[0];
  reduloctab[1]  = - grafptr->compglbload[0];
  reduloctab[2]  =   grafptr->compglbload[1];
  reduloctab[3]  = - grafptr->compglbload[1];
  reduloctab[4]  =   grafptr->compglbload[2];
  reduloctab[5]  = - grafptr->compglbload[2];
  reduloctab[6]  =   grafptr->compglbsize[2];
  reduloctab[7]  = - grafptr->compglbsize[2];
  reduloctab[8]  =   grafptr->levlnum;
  reduloctab[9]  = - grafptr->levlnum;
  reduloctab[10] =   cheklocval;

  if (MPI_Allreduce (reduloctab, reduglbtab, 11, GNUM_MPI, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphCheck: communication error (2)");
    return     (1);
  }

  if (reduglbtab[10] != 0) {                      /* Return from previous errors */
    if (partgsttax != NULL)
      memFree (partgsttax);
    return (1);
  }

  if ((reduglbtab[1] != - reduglbtab[0]) ||
      (reduglbtab[3] != - reduglbtab[2]) ||
      (reduglbtab[5] != - reduglbtab[4]) ||
      (reduglbtab[7] != - reduglbtab[6]) ||
      (reduglbtab[9] != - reduglbtab[8])) {
    errorPrint ("vdgraphCheck: inconsistent global graph data");
    return     (1);
  }

  memCpy (partgsttax, grafptr->partgsttax + grafptr->s.baseval, grafptr->s.vertlocnbr); /* Copy local part data           */
  dgraphHaloSync (&grafdat, partgsttax, GRAPHPART_MPI); /* Spread yet unbased halo part data across neighboring processes */
  partgsttax -= grafptr->s.baseval;

  complocload[0] =
  complocload[1] =
  complocload[2] = 0;
  complocsize[0] =
  complocsize[1] =
  complocsize[2] = 0;
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertlocnnd; vertnum ++) {
    int                 partnum;                  /* Part of current vertex */
    Gnum                edgenum;                  /* Number of current edge */

    partnum = (int) partgsttax[vertnum];

    complocload[partnum] += (grafptr->s.veloloctax == NULL) ? 1 : grafptr->s.veloloctax[vertnum];
    complocsize[partnum] ++;

    commcut[0] =
    commcut[1] =
    commcut[2] = 0;
    for (edgenum = grafptr->s.vertloctax[vertnum]; edgenum < grafptr->s.vendloctax[vertnum]; edgenum ++) {
      if (grafdat.edgegsttax[edgenum] < grafptr->s.vertlocnnd) /* Check only for local ends, as ghost part might be inaccurate */
        commcut[partgsttax[grafdat.edgegsttax[edgenum]]] ++;
    }
    if (partnum != 2) {
      if (commcut[1 - partnum] != 0) {
        errorPrint ("vdgraphCheck: vertex should be in separator (%ld)", (long) vertnum);
        cheklocval = 1;
        break;
      }
    }
  }
  if (grafptr->s.edgegsttax != grafdat.edgegsttax) /* If ghost edge array was allocated here, free it manually */
    memFree (grafdat.edgegsttax + grafptr->s.baseval);
  if (grafptr->s.procsidtab != grafdat.procsidtab) /* The same for procsidtab */
    memFree (grafdat.procsidtab);
  memFree (partgsttax + grafptr->s.baseval);

  if ((cheklocval == 0) &&
      ((complocsize[0] != grafptr->complocsize[0]) ||
       (complocsize[1] != grafptr->complocsize[1]) ||
       (complocsize[2] != grafptr->complocsize[2]))) {
    errorPrint ("vgraphCheck: invalid local part sizes");
    cheklocval = 1;
  }

  reduloctab[0] = complocload[0];
  reduloctab[1] = complocload[1];
  reduloctab[2] = complocload[2];
  reduloctab[3] = complocsize[0];
  reduloctab[4] = complocsize[1];
  reduloctab[5] = complocsize[2];
  reduloctab[6] = cheklocval;

  if (MPI_Allreduce (reduloctab, reduglbtab, 7, GNUM_MPI, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphCheck: communication error (3)");
    return     (1);
  }
  if (reduglbtab[6] != 0)                         /* Return from previous errors */
    return (1);

  if ((grafptr->compglbload[0] != reduglbtab[0]) ||
      (grafptr->compglbload[1] != reduglbtab[1]) ||
      (grafptr->compglbload[2] != reduglbtab[2])) {
    errorPrint ("vdgraphCheck: invalid global part loads");
    cheklocval = 1;
  }
  if ((grafptr->compglbsize[0] != reduglbtab[3]) ||
      (grafptr->compglbsize[1] != reduglbtab[4]) ||
      (grafptr->compglbsize[2] != reduglbtab[5])) {
    errorPrint ("vgraphCheck: invalid global part sizes");
    cheklocval = 1;
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphCheck: communication error (4)");
    return     (1);
  }

  return (chekglbval);
}
