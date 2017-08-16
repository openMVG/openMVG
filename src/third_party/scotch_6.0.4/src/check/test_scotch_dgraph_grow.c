/* Copyright 2012,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_dgraph_grow.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_dgraphGrow() routine.        **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 26 sep 2012     **/
/**                                 to     02 mar 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <mpi.h>
#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <pthread.h>
#include <unistd.h>

#include "ptscotch.h"

#define errorProg                   SCOTCH_errorProg
#define errorPrint                  SCOTCH_errorPrint

void                        _SCOTCHintRandInit  (void);
SCOTCH_Num                  _SCOTCHintRandVal   (SCOTCH_Num);

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (

int                 argc,
char *              argv[])
{
  MPI_Status            statdat;
  MPI_Comm              proccomm;
  int                   procglbnbr;               /* Number of processes sharing graph data */
  int                   proclocnum;               /* Number of this process                 */
  long                  vertlocadj;
  SCOTCH_Num            vertglbnbr;
  SCOTCH_Num            vertlocnbr;
  SCOTCH_Num            vertgstnbr;
  SCOTCH_Num *          seedloctab;
  SCOTCH_Num *          partgsttab;
  SCOTCH_Num            baseval;
  SCOTCH_Dgraph         grafdat;
  FILE *                file;
  int                   procnum;
#ifdef SCOTCH_PTHREAD
  int                 thrdlvlreqval;
  int                 thrdlvlproval;
#endif /* SCOTCH_PTHREAD */

  errorProg (argv[0]);

#ifdef SCOTCH_PTHREAD
  thrdlvlreqval = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdlvlreqval, &thrdlvlproval) != MPI_SUCCESS)
    errorPrint ("main: Cannot initialize (1)");
  if (thrdlvlreqval > thrdlvlproval)
    errorPrint ("main: MPI implementation is not thread-safe: recompile without SCOTCH_PTHREAD");
#else /* SCOTCH_PTHREAD */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    errorPrint ("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD */

  if (argc != 2) {
    errorPrint ("main: invalid number of parameters");
    exit       (1);
  }

  proccomm = MPI_COMM_WORLD;
  MPI_Comm_size (proccomm, &procglbnbr);          /* Get communicator data */
  MPI_Comm_rank (proccomm, &proclocnum);

  fprintf (stderr, "Proc %2d of %2d, pid %d\n", proclocnum, procglbnbr, getpid ());

#ifdef SCOTCH_CHECK_NOAUTO
  if (proclocnum == 0) {                          /* Synchronize on keybord input */
    char           c;

    printf ("Waiting for key press...\n");
    scanf ("%c", &c);
  }
#endif /* SCOTCH_CHECK_NOAUTO */

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    errorPrint ("main: cannot communicate");
    return     (1);
  }

  if (SCOTCH_dgraphInit (&grafdat, proccomm) != 0) { /* Initialize source graph */
    errorPrint ("main: cannot initialize graph");
    return     (1);
  }

  file = NULL;
  if ((proclocnum == 0) &&
      ((file = fopen (argv[1], "r")) == NULL)) {
    errorPrint ("main: cannot open graph file");
    return     (1);
  }

  if (SCOTCH_dgraphLoad (&grafdat, file, 0, 0) != 0) {
    errorPrint ("main: cannot load graph");
    return     (1);
  }

  if (file != NULL)
    fclose (file);

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    errorPrint ("main: cannot communicate");
    return     (1);
  }

  if (SCOTCH_dgraphGhst (&grafdat) != 0) {
    errorPrint ("main: cannot compute ghost edge array");
    return     (1);
  }

  SCOTCH_dgraphData (&grafdat, &baseval, &vertglbnbr, &vertlocnbr, NULL, &vertgstnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  if ((seedloctab = malloc (vertlocnbr * sizeof (SCOTCH_Num))) == NULL) {
    errorPrint ("main: cannot allocate seed array");
    return     (1);
  }

  if ((partgsttab = malloc (vertgstnbr * sizeof (SCOTCH_Num))) == NULL) {
    errorPrint ("main: cannot allocate part array");
    return     (1);
  }

  memset (partgsttab, ~0, vertgstnbr * sizeof (SCOTCH_Num));

  _SCOTCHintRandInit ();
  seedloctab[0] = _SCOTCHintRandVal (vertlocnbr);
  seedloctab[1] = _SCOTCHintRandVal (vertlocnbr);
  seedloctab[2] = _SCOTCHintRandVal (vertlocnbr);

  partgsttab[seedloctab[0] - baseval] = 0;
  partgsttab[seedloctab[1] - baseval] = 1;
  partgsttab[seedloctab[2] - baseval] = 2;

  if (SCOTCH_dgraphGrow (&grafdat, 3, seedloctab, 4, partgsttab) != 0) {
    errorPrint ("main: cannot compute grown regions");
    return     (1);
  }

  free (seedloctab);

  for (procnum = 0; procnum < procglbnbr; procnum ++) {
    SCOTCH_Num          vertlocnum;

    MPI_Barrier (proccomm);

    if (procnum == proclocnum) {
      if ((file = fopen ("/tmp/test_scotch_dgraph_grow.map", (procnum == 0) ? "w" : "a+")) == NULL) {
        errorPrint ("main: cannot open mapping file");
        return     (1);
      }

      if (procnum == 0) {
        fprintf (file, "%ld\n", (long) vertglbnbr);
        vertlocadj = (long) baseval;
      }
      else
        MPI_Recv (&vertlocadj, 1, MPI_LONG, procnum - 1, 0, MPI_COMM_WORLD, &statdat);

      for (vertlocnum = 0; vertlocnum < vertlocnbr; vertlocnum ++)
        fprintf (file, "%ld\t%ld\n", vertlocadj + (long) vertlocnum, (long) partgsttab[vertlocnum]);

      fclose (file);

      if (procnum < (procglbnbr - 1)) {
        vertlocadj += (long) vertlocnbr;
        MPI_Send (&vertlocadj, 1, MPI_LONG, procnum + 1, 0, MPI_COMM_WORLD);
      }
    }
  }

  SCOTCH_dgraphExit (&grafdat);

  MPI_Finalize ();
  exit         (0);
}
