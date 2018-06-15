/* Copyright 2008,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_match_check.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed graph matching      **/
/**                routines.                               **/
/**                                                        **/
/**    DATES     : # Version 5.1  : from : 25 dec 2008     **/
/**                                 to   : 08 apr 2009     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_MATCH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "dgraph_match.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* matchings.                        */
/*                                   */
/*************************************/

/* This routine checks the consistency of a
** given complete matching.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphMatchCheck (
DgraphMatchData * restrict const    mateptr)
{
  Gnum                baseval;
  Gnum * restrict     flaggsttax;
  int                 procngbnum;
  Gnum                multlocnbr;
  Gnum                multlocnum;
  Gnum                vertglbnnd;
  Gnum                vertlocnbr;
  Gnum                vertlocnnd;
  Gnum                vertlocnum;
  Gnum                vertlocadj;
  int                 cheklocval;
  int                 chekglbval;

  Dgraph * restrict const                   grafptr    = mateptr->c.finegrafptr;
  const int * restrict const                procngbtab = grafptr->procngbtab;
  const Gnum * restrict const               mategsttax = mateptr->mategsttax;
  DgraphCoarsenVert * restrict const        vsnddattab = mateptr->c.vsnddattab;
  const DgraphCoarsenMulti * restrict const multloctab = mateptr->c.multloctab;
  const int * restrict const                procgsttax = mateptr->c.procgsttax;
  const Gnum * restrict const               edgeloctax = grafptr->edgeloctax;
  const Gnum * restrict const               edgegsttax = grafptr->edgegsttax;
  const Gnum * restrict const               vertloctax = grafptr->vertloctax;
  const Gnum * restrict const               vendloctax = grafptr->vendloctax;
  int * restrict const                      nsndidxtab = mateptr->c.nsndidxtab;

  baseval = grafptr->baseval;

  cheklocval = 0;

  multlocnbr = mateptr->c.multlocnbr;
  if ((multlocnbr < 0) || (multlocnbr > grafptr->vertlocnbr)) {
    errorPrint ("dgraphMatchCheck: invalid number of multinodes");
    cheklocval = 1;
  }

  vertlocnbr = grafptr->vertlocnbr;
  for (vertlocnum = baseval; vertlocnum < vertlocnbr; vertlocnum ++) {
    if (mategsttax[vertlocnum] < 0) {
      errorPrint ("dgraphMatchCheck: unmatched local vertex");
      cheklocval = 1;
      break;
    }
  }

  if ((flaggsttax = memAlloc (grafptr->vertgstnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphMatchCheck: out of memory");
    cheklocval = 1;
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, mateptr->c.finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchCheck: communication error (1)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    if (flaggsttax != NULL)
      memFree (flaggsttax);
    return (1);
  }

  for (procngbnum = 0; procngbnum < grafptr->procngbnbr; procngbnum ++) /* Reset indices for sending messages */
    nsndidxtab[procngbnum] = mateptr->c.vsnddsptab[procngbtab[procngbnum]];

  memSet (flaggsttax, ~0, grafptr->vertgstnbr * sizeof (Gnum));
  flaggsttax -= baseval;

  vertglbnnd = grafptr->vertglbnbr + baseval;
  vertlocnnd = grafptr->vertlocnnd;
  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - baseval;
  for (multlocnum = 0; multlocnum < multlocnbr; multlocnum ++) {
    Gnum                vertglbnum;
    Gnum                vertlocnum;
    Gnum                vertglbend;

    vertglbnum = multloctab[multlocnum].vertglbnum[0];
    vertlocnum = vertglbnum - vertlocadj;         /* First vertex is always local */
    if ((vertlocnum < baseval) || (vertlocnum >= vertlocnnd)) {
      errorPrint ("dgraphMatchCheck: invalid multinode vertex (1)");
      goto abort;
    }
    if (flaggsttax[vertlocnum] != -1) {
      errorPrint ("dgraphMatchCheck: duplicate multinode vertex (1)");
      goto abort;
    }
    flaggsttax[vertlocnum] = multlocnum + vertlocadj;

    vertglbend = multloctab[multlocnum].vertglbnum[1];
    if (vertglbend < 0) {                         /* If end vertex is remote */
      Gnum                edgelocnum;
      Gnum                vertgstend;
      int                 vsndidxnum;
      int                 procngbnum;

      edgelocnum = -2 - vertglbend;
      if ((edgelocnum < grafptr->baseval) ||
          (edgelocnum >= (grafptr->edgelocsiz + grafptr->baseval))) {
        errorPrint ("dgraphMatchCheck: invalid multinode vertex (2)");
        goto abort;
      }

      vertglbend = edgeloctax[edgelocnum];

      if (mategsttax[vertlocnum] != vertglbend) {
        errorPrint ("dgraphMatchCheck: invalid mate array (1)");
        goto abort;
      }

      vertgstend = edgegsttax[edgelocnum];

      if (flaggsttax[vertgstend] != -1) {
        errorPrint ("dgraphMatchCheck: duplicate multinode vertex (2)");
        goto abort;
      }
      flaggsttax[vertgstend] = multlocnum + vertlocadj;

      if (mategsttax[vertgstend] != vertglbnum) {
        errorPrint ("dgraphMatchCheck: invalid mate array (2)");
        goto abort;
      }

      procngbnum = procgsttax[vertgstend];        /* Find neighbor owner process                                      */
      if ((procngbnum < 0) || (procngbnum >= grafptr->procngbnbr)) { /* If neighbor had not been computed or is wrong */
        errorPrint ("dgraphMatchCheck: internal error (1)");
        goto abort;
      }
      if ((grafptr->procvrttab[procngbtab[procngbnum]]     >  vertglbend) ||
          (grafptr->procvrttab[procngbtab[procngbnum] + 1] <= vertglbend)) {
        errorPrint ("dgraphMatchCheck: internal error (2)");
        goto abort;
      }

      vsndidxnum = nsndidxtab[procngbnum] ++;     /* Get position of message in send array */
      if (vsndidxnum >= mateptr->c.vsnddsptab[procngbtab[procngbnum] + 1]) {
        errorPrint ("dgraphMatchCheck: internal error (3)");
        goto abort;
      }
      vsnddattab[vsndidxnum].datatab[0] = vertglbnum;
      vsnddattab[vsndidxnum].datatab[1] = vertglbend;
    }
    else {                                        /* End vertex is local */
      Gnum                vertlocend;
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      if (mategsttax[vertlocnum] != vertglbend) {
        errorPrint ("dgraphMatchCheck: invalid mate array (3)");
        goto abort;
      }

      if (vertglbend == vertglbnum)               /* If single multinode */
        continue;

      vertlocend = vertglbend - vertlocadj;
      if ((vertlocend < baseval) || (vertlocend >= vertlocnnd)) {
        errorPrint ("dgraphMatchCheck: invalid multinode vertex (3)");
        goto abort;
      }

      edgelocnum = vertloctax[vertlocnum];
      edgelocnnd = vendloctax[vertlocnum];
      if (edgelocnum != edgelocnnd) {             /* If first multinode vertex is not an isolated vertex */
        for ( ; ; edgelocnum ++) {                /* Loop on edges of first multinode vertex             */
          if (edgelocnum >= edgelocnnd) {         /* If not a valid neighbor                             */
            errorPrint ("dgraphMatchCheck: invalid multinode vertex (4)");
            goto abort;
          }
          if (edgeloctax[edgelocnum] == vertglbend) /* If edge to end vertex found */
            break;
        }
      }

      if (flaggsttax[vertlocend] != -1) {
        errorPrint ("dgraphMatchCheck: duplicate multinode vertex (3)");
        goto abort;
      }
      flaggsttax[vertlocend] = multlocnum + vertlocadj;

      if (mategsttax[vertlocend] != vertglbnum) {
        errorPrint ("dgraphMatchCheck: invalid mate array (4)");
        goto abort;
      }
    }
  }
  cheklocval = -1;
abort:
  cheklocval ++;

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, mateptr->c.finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphMatchCheck: communication error (2)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    memFree (flaggsttax + baseval);
    return (1);
  }

/* TODO: Send messages and check consistency of matching on the receiving side */

  memFree (flaggsttax + baseval);

  return (0);
}
