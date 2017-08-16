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
/**   NAME       : test_scotch_graph_order.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphOrderCompute*()         **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 aug 2014     **/
/**                                 to     29 aug 2014     **/
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
  FILE *              fileptr;
  SCOTCH_Graph        grafdat;
  SCOTCH_Ordering     ordedat;
  SCOTCH_Strat        stradat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num          vertnum;
  SCOTCH_Num          listnbr;
  SCOTCH_Num          listnum;
  SCOTCH_Num *        listtab;

  SCOTCH_errorProg (argv[0]);

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    return            (1);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file (1)");
    return            (1);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    return            (1);
  }

  fclose (fileptr);

  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  listnbr = (vertnbr + 1) / 2;                    /* Only keep half of the vertices in induced graph */
  if ((listtab = malloc (listnbr * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (1)");
    return            (1);
  }
  for (listnum = 0, vertnum = baseval + (listnbr / 4); /* Keep only middle half of the vertices */
       listnum < listnbr; listnum ++, vertnum ++)
    listtab[listnum] = vertnum;

  if ((fileptr = tmpfile ()) == NULL) {           /* Open temporary file for resulting output */
    SCOTCH_errorPrint ("main: cannot open file (2)");
    return            (1);
  }

  if (SCOTCH_stratInit (&stradat) != 0) {         /* Initialize ordering strategy */
    SCOTCH_errorPrint ("main: cannot initialize strategy");
    return            (1);
  }

  if (SCOTCH_graphOrderInit (&grafdat, &ordedat, NULL, NULL, NULL, NULL, NULL) != 0) { /* Initialize ordering */
    SCOTCH_errorPrint ("main: cannot initialize ordering (1)");
    return            (1);
  }

  if (SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat) != 0) {
    SCOTCH_errorPrint ("main: cannot order graph");
    return            (1);
  }

  if (SCOTCH_graphOrderCheck (&grafdat, &ordedat) != 0) {
    SCOTCH_errorPrint ("main: invalid ordering (1)");
    return            (1);
  }

  SCOTCH_graphOrderSave     (&grafdat, &ordedat, fileptr); /* Test ordering data output routines */
  SCOTCH_graphOrderSaveMap  (&grafdat, &ordedat, fileptr);
  SCOTCH_graphOrderSaveTree (&grafdat, &ordedat, fileptr);

  SCOTCH_graphOrderExit (&grafdat, &ordedat);     /* Free computed ordering */

  if (SCOTCH_graphOrderInit (&grafdat, &ordedat, NULL, NULL, NULL, NULL, NULL) != 0) { /* Initialize ordering again */
    SCOTCH_errorPrint ("main: cannot initialize ordering (2)");
    return            (1);
  }

  if (SCOTCH_graphOrderComputeList (&grafdat, &ordedat, listnbr, listtab, &stradat) != 0) {
    SCOTCH_errorPrint ("main: cannot order induced graph");
    return            (1);
  }

  if (SCOTCH_graphOrderCheck (&grafdat, &ordedat) != 0) {
    SCOTCH_errorPrint ("main: invalid ordering (2)");
    return            (1);
  }

  SCOTCH_graphOrderSave     (&grafdat, &ordedat, fileptr); /* Test ordering data output routines */
  SCOTCH_graphOrderSaveMap  (&grafdat, &ordedat, fileptr);
  SCOTCH_graphOrderSaveTree (&grafdat, &ordedat, fileptr);

  free (listtab);
  SCOTCH_stratExit      (&stradat);
  SCOTCH_graphOrderExit (&grafdat, &ordedat);
  SCOTCH_graphExit      (&grafdat);

  return (0);
}
