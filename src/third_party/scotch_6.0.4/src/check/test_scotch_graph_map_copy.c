/* Copyright 2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_graph_map_copy.c            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphMap*() routines in the  **/
/**                specific case of the "copy" method.     **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 15 oct 2014     **/
/**                                 to     28 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <math.h>
#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>

#include "scotch.h"

#define STRANBR                     3

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
  SCOTCH_Mapping          mappdat;                /* Mapping to compute */
  SCOTCH_Mapping          mapodat;                /* Old mapping        */
  FILE *                  fileptr;
  SCOTCH_Graph            grafdat;
  SCOTCH_Num              xdimsiz;
  int                     archnum;
  SCOTCH_Arch             archdat;
  SCOTCH_Strat            stratab[STRANBR];
  int                     stranum;
  int                     typenum;
  SCOTCH_Num              baseval;
  SCOTCH_Num              vertnbr;
  SCOTCH_Num              vertnum;
  SCOTCH_Num *            parttab;
  SCOTCH_Num *            parotab;

  SCOTCH_errorProg (argv[0]);

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    return            (1);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) { /* Read the givel graph */
    SCOTCH_errorPrint ("main: cannot open file (1)");
    return            (1);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    return            (1);
  }

  fclose (fileptr);

  SCOTCH_graphSize (&grafdat, &vertnbr, NULL);

  if (((parttab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((parotab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    return            (1);
  }

  for (stranum = 0; stranum < STRANBR; stranum ++) { /* Initialize mapping strategies */
    if (SCOTCH_stratInit (&stratab[stranum]) != 0) {
      SCOTCH_errorPrint ("main: cannot initialize strategy");
      return            (1);
    }
  }
  SCOTCH_stratGraphMap (&stratab[0], "cf{move=10000,pass=-1,bal=0.05}");
  SCOTCH_stratGraphMap (&stratab[1], "m{vert=120,low=cf{move=10000,pass=-1,bal=0.05},asc=b{bnd=f{move=10000,pass=-1,bal=0.05},org=f{move=10000,pass=-1,bal=0.05}}}");

  if (SCOTCH_archInit (&archdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize architecture");
    return            (1);
  }
  SCOTCH_archCmplt (&archdat, 5);

  for (stranum = 0; stranum < (STRANBR - 1); stranum ++) {
    for (typenum = 0; typenum < 2; typenum ++) {
      int                 i;
      int                 o;

      printf ("Strat %d, type %d\n", stranum, typenum);

      switch (typenum) {
        case 0 :                                  /* Plain mapping */
          if (SCOTCH_graphMapInit (&grafdat, &mappdat, &archdat, parttab) != 0) { /* Initialize new mapping */
            SCOTCH_errorPrint ("main: cannot initialize mapping (1)");
            return            (1);
          }

          o = SCOTCH_graphMapCompute (&grafdat, &mappdat, &stratab[STRANBR - 1]); /* Last strategy is plain mapping strategy */
          memcpy (parotab, parttab, vertnbr * sizeof (SCOTCH_Num)); /* Use plain mapping as old mapping in the following     */
          break;
        case 1 :                                  /* Remapping with copy of the old partition array         */
          if (SCOTCH_graphMapInit (&grafdat, &mapodat, &archdat, parotab) != 0) { /* Initialize old mapping */
            SCOTCH_errorPrint ("main: cannot initialize mapping (2)");
            return            (1);
          }

          o = SCOTCH_graphRemapCompute (&grafdat, &mappdat, &mapodat, 0, NULL, &stratab[stranum]);
          break;
      }

      if (o != 0) {
        SCOTCH_errorPrint ("main: cannot compute mapping");
        return (1);
      }
    }

    SCOTCH_graphMapExit (&grafdat, &mapodat);
    SCOTCH_graphMapExit (&grafdat, &mappdat);
  }


  SCOTCH_archExit (&archdat);

  for (stranum = 0; stranum < STRANBR; stranum ++)
    SCOTCH_stratExit (&stratab[stranum]);

  free             (parotab);
  free             (parttab);
  SCOTCH_graphExit (&grafdat);

  return (0);
}
