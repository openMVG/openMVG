/* Copyright 2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_graph_coarsen.c             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphColor() routine.        **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 15 oct 2014     **/
/**                                 to     15 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>

#include "scotch.h"

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
  SCOTCH_Num          baseval;
  FILE *              fileptr;
  SCOTCH_Num          finevertnbr;
  SCOTCH_Graph        finegrafdat;
  double              coarrat;
  SCOTCH_Num          coarvertmax;
  SCOTCH_Graph        coargrafdat;
  SCOTCH_Num *        coarmulttab;
  int                 o;

  SCOTCH_errorProg (argv[0]);

  if (SCOTCH_graphInit (&finegrafdat) != 0) {     /* Initialize source, fine graph */
    SCOTCH_errorPrint ("main: cannot initialize fine graph");
    return            (1);
  }

  if (SCOTCH_graphInit (&coargrafdat) != 0) {     /* Initialize coarse graph */
    SCOTCH_errorPrint ("main: cannot initialize coarse graph");
    return            (1);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    return            (1);
  }

  if (SCOTCH_graphLoad (&finegrafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    return            (1);
  }

  fclose (fileptr);

  SCOTCH_graphData (&finegrafdat, &baseval, &finevertnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  coarrat     = 0.8;
  coarvertmax = (SCOTCH_Num) (coarrat * (double) finevertnbr);

  if ((coarmulttab = malloc (coarvertmax * 2 * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (1)");
    return            (1);
  }

  if ((o = SCOTCH_graphCoarsen (&finegrafdat, &coargrafdat, coarmulttab, 1, coarrat)) >= 2) {
    SCOTCH_errorPrint ("main: cannot coarsen graph");
    return            (1);
  }

  if (o == 0) {
    SCOTCH_Num          finevertnnd;
    SCOTCH_Num          coarvertnbr;
    SCOTCH_Num          coarvertnum;
    SCOTCH_Num          coarverttmp;

    SCOTCH_graphSize (&coargrafdat, &coarvertnbr, NULL);

    printf ("Graph coarsened with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

    for (coarvertnum = 0, coarverttmp = finevertnbr, finevertnnd = finevertnbr + baseval;
         coarvertnum < coarvertnbr; coarvertnum ++) {
      SCOTCH_Num          finevertnum0;
      SCOTCH_Num          finevertnum1;

      finevertnum0 = coarmulttab[2 * coarvertnum];
      finevertnum1 = coarmulttab[2 * coarvertnum + 1];
      if ((finevertnum0 <  baseval)     ||
          (finevertnum0 >= finevertnnd) ||
          (finevertnum1 <  baseval)     ||
          (finevertnum1 >= finevertnnd)) {
        SCOTCH_errorPrint ("main: invalid multinode array (1)");
        return            (1);
      }
      if (finevertnum0 != finevertnum1)
        coarverttmp --;
    }
    if (coarverttmp != coarvertnbr) {
      SCOTCH_errorPrint ("main: invalid multinode array (2)");
      return            (1);
    }
  }
  else
    printf ("Graph could not be coarsened with a ratio of %lg\n", coarrat);

  free (coarmulttab);
  SCOTCH_graphExit (&coargrafdat);
  SCOTCH_graphExit (&finegrafdat);

  return (0);
}
