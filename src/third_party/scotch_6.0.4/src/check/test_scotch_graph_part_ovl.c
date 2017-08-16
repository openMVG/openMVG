/* Copyright 2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_graph_part_ovl.c            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the sequential        **/
/**                graph partitioning with overlap         **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 20 sep 2014     **/
/**                                 to     20 sep 2014     **/
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
  SCOTCH_Graph          grafdat;
  SCOTCH_Strat          stradat;
  SCOTCH_Num            baseval;
  SCOTCH_Num            partnbr;
  SCOTCH_Num            partnum;
  SCOTCH_Num * restrict parttax;
  SCOTCH_Num            vertnbr;
  SCOTCH_Num            vertnum;
  SCOTCH_Num *          verttab;
  SCOTCH_Num *          vendtab;
  SCOTCH_Num *          velotab;
  SCOTCH_Num *          vlbltab;
  SCOTCH_Num *          edgetax;
  SCOTCH_Num * restrict flagtab;
  SCOTCH_Num * restrict loadtab;
  SCOTCH_Num            loadmin;
  SCOTCH_Num            loadmax;
  SCOTCH_Num            loadsum;
  double                loadavg;
  FILE *                fileptr;

  SCOTCH_errorProg (argv[0]);

  if ((argc < 4) || (argc > 5)) {
    SCOTCH_errorPrint ("main: usage is \"%s <nparts> <input source graph file> <output mapping file> [<strategy>]\"\n", argv[0]);
    exit              (1);
  }

  if ((partnbr = (SCOTCH_Num) atoi (argv[1])) < 1) {
    SCOTCH_errorPrint ("main: invalid number of parts (\"%s\")", argv[1]);
    return            (1);
  }

  if (SCOTCH_stratInit (&stradat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize strategy");
    return            (1);
  }
  if (argc == 5) {
    if (SCOTCH_stratGraphPartOvl (&stradat, argv[4]) != 0) {
      SCOTCH_errorPrint ("main: invalid user-provided strategy");
      return            (1);
    }
  }

  if (SCOTCH_graphInit (&grafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize graph");
    return            (1);
  }

  if ((fileptr = fopen (argv[2], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file (1)");
    return            (1);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) {
    SCOTCH_errorPrint ("main: cannot load graph");
    return            (1);
  }

  fclose (fileptr);

  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, &verttab, &vendtab, &velotab, &vlbltab, NULL, &edgetax, NULL);

  if (((parttax = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((flagtab = malloc (partnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((loadtab = malloc (partnbr * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    return            (1);
  }

  if (SCOTCH_graphPartOvl (&grafdat, partnbr, &stradat, parttax) != 0) { /* Parttax is not based yet */
    SCOTCH_errorPrint ("main: cannot compute mapping");
    return            (1);
  }

  edgetax -= baseval;
  parttax -= baseval;

  memset (loadtab,  0, partnbr * sizeof (SCOTCH_Num)); /* Part loads set to 0                */
  memset (flagtab, ~0, partnbr * sizeof (SCOTCH_Num)); /* Flags set to invalid vertex number */

  for (vertnum = 0; vertnum < vertnbr; vertnum ++) {
    SCOTCH_Num          veloval;
    SCOTCH_Num          partval;

    veloval = (velotab == NULL) ? 1 : velotab[vertnum];
    partval = parttax[vertnum + baseval];         /* vertnum is not based               */
    if (partval >= 0)                             /* If vertex belongs to one part only */
      loadtab[partval] += veloval;                /* Add vertex load to this part       */
    else {                                        /* Vertex belongs to several parts    */
      SCOTCH_Num          edgenum;

      for (edgenum = verttab[vertnum]; edgenum < vendtab[vertnum]; edgenum ++) {
        SCOTCH_Num          vertend;
        SCOTCH_Num          partend;

        vertend = edgetax[edgenum];
        partend = parttax[vertend];               /* vertend is based                     */
        if (partend < 0)                          /* If neighbor has no identifiable part */
          continue;
        if (flagtab[partend] == vertnum)          /* If neighbor part already accounted for, skip it */
          continue;

        loadtab[partend] += veloval;              /* Vertex load contributes to this part */
        flagtab[partend]  = vertnum;              /* Record a contribution has been made  */
      }
    }
  }

  loadsum =
  loadmax = 0;
  loadmin = SCOTCH_NUMMAX;
  for (partnum = 0; partnum < partnbr; partnum ++) {
    loadsum += loadtab[partnum];
    if (loadtab[partnum] > loadmax)
      loadmax = loadtab[partnum];
    if (loadtab[partnum] < loadmin)
      loadmin = loadtab[partnum];

    printf ("M\tCompload[%02ld]\t%ld\n",
            (long) partnum,
            (long) loadtab[partnum]);
  }

  loadavg = (double) loadsum / (double) partnbr;
  printf ("M\tCompLoadAvg\t%g\n",
	  (double) loadavg);
  printf ("M\tCompLoadMax/Avg\t%g\n",
	  (double) loadmax / loadavg);

  if ((fileptr = fopen (argv[3], "w")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file (2)");
    return            (1);
  }

  if (fprintf (fileptr, "%ld\n", (long) vertnbr) == EOF) {
    SCOTCH_errorPrint ("main: bad output (1)");
    return            (1);
  }

  for (vertnum = 0; vertnum < vertnbr; vertnum ++) {
    if (fprintf (fileptr, "%ld\t%ld\n",
                 (long) ((vlbltab == NULL) ? vertnum : vlbltab[vertnum]),
                 (long) parttax[vertnum + baseval]) == EOF) {
      SCOTCH_errorPrint ("main: bad output (2)");
      return            (1);
    }
  }

  fclose (fileptr);

  free (loadtab);
  free (flagtab);
  free (parttax + baseval);

  SCOTCH_stratExit (&stradat);
  SCOTCH_graphExit (&grafdat);

  exit (0);
}
