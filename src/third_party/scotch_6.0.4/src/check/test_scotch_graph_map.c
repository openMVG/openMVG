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
/**   NAME       : test_scotch_graph_map.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphMap*() routines.        **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 12 aug 2014     **/
/**                                 to     20 sep 2014     **/
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

#define ARCHNBR                     4
#define STRANBR                     2

#define COORD(x,y)                  ((y) * xdimsiz + (x))

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
  SCOTCH_Arch             archtab[ARCHNBR];
  SCOTCH_Strat            stratab[STRANBR];
  int                     stranum;
  int                     typenum;
  SCOTCH_Num              baseval;
  SCOTCH_Num              vertnbr;
  SCOTCH_Num              vertnum;
  SCOTCH_Num *            parttab;
  SCOTCH_Num *            parotab;
  SCOTCH_Num *            vmlotab;
  SCOTCH_Num *            vmloptr;                /* vmlotab or NULL */

  SCOTCH_errorProg (argv[0]);

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    return            (1);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) { /* Read a square 2D grid graph */
    SCOTCH_errorPrint ("main: cannot open file (1)");
    return            (1);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    return            (1);
  }

  fclose (fileptr);

  SCOTCH_graphSize (&grafdat, &vertnbr, NULL);
  xdimsiz = (SCOTCH_Num) sqrt ((double) vertnbr);
  if (vertnbr != (xdimsiz * xdimsiz)) {
    SCOTCH_errorPrint ("main: graph is not a square grid");
    return            (1);
  }

  if (((parttab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((parotab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((vmlotab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    return            (1);
  }

  for (vertnum = 0; vertnum < vertnbr; vertnum ++) /* Fill vertex migration load array */
    vmlotab[vertnum] = vertnum % 3;

  for (stranum = 0; stranum < STRANBR; stranum ++) { /* Initialize mapping strategies */
    if (SCOTCH_stratInit (&stratab[stranum]) != 0) {
      SCOTCH_errorPrint ("main: cannot initialize strategy");
      return            (1);
    }
  }
  SCOTCH_stratGraphMapBuild (&stratab[0], SCOTCH_STRATRECURSIVE, 4, 0.05);
  SCOTCH_stratGraphMapBuild (&stratab[1], SCOTCH_STRATDEFAULT,   4, 0.05);

  for (archnum = 0; archnum < ARCHNBR; archnum ++) { /* Initialize architectures */
    if (SCOTCH_archInit (&archtab[archnum]) != 0) {
      SCOTCH_errorPrint ("main: cannot initialize architecture");
      return            (1);
    }
  }
  SCOTCH_archCmplt (&archtab[0], 5);
  SCOTCH_archMesh2 (&archtab[1], 2, 2);
  SCOTCH_archMesh2 (&archtab[2], xdimsiz * 2, xdimsiz * 2); /* Oversized architecture */
  SCOTCH_archVhcub (&archtab[3]);

  if ((fileptr = tmpfile ()) == NULL) {           /* Open temporary file for resulting output */
    SCOTCH_errorPrint ("main: cannot open file (2)");
    return            (1);
  }

  for (stranum = 0; stranum < STRANBR; stranum ++) {
    for (archnum = 0; archnum < ARCHNBR; archnum ++) {
      SCOTCH_Num              archsiz;

      if (SCOTCH_graphMapInit (&grafdat, &mappdat, &archtab[archnum], parttab) != 0) { /* Initialize new mapping */
        SCOTCH_errorPrint ("main: cannot initialize mapping (1)");
        return            (1);
      }
      if (SCOTCH_graphMapInit (&grafdat, &mapodat, &archtab[archnum], parotab) != 0) { /* Initialize old mapping */
        SCOTCH_errorPrint ("main: cannot initialize mapping (2)");
        return            (1);
      }

      archsiz = SCOTCH_archSize (&archtab[archnum]);

      for (typenum = 0; typenum < 6; typenum ++) {
        int                 i;
        int                 o;

        memset (parttab, ~0, xdimsiz * xdimsiz * sizeof (SCOTCH_Num)); /* Assume all vertices are not fixed */
        if (archnum < 2) {                        /* For fixed-size architectures                           */
          for (i = 0; i < (xdimsiz - 1); i ++) {    /* Place fixed vertices at all four sides               */
            parttab[COORD (0, i)] = 0;
            parttab[COORD (i + 1, 0)] = 1;
            parttab[COORD (xdimsiz - 1, i + 1)] = archsiz - 2;
            parttab[COORD (i, xdimsiz - 1)] = archsiz - 1;
          }
        }
        else {                                    /* For variable-sized architectures       */
          for (i = 0; i < (xdimsiz - 1); i ++) {  /* Place fixed vertices at all four sides */
            parttab[COORD (0, i)] = vertnbr - 2;
            parttab[COORD (i + 1, 0)] = vertnbr - 1;
            parttab[COORD (xdimsiz - 1, i + 1)] = vertnbr;
            parttab[COORD (i, xdimsiz - 1)] = vertnbr + 1;
          }
        }

        printf ("Strat %d, arch %d, type %d\n", stranum, archnum, typenum);

        vmloptr = vmlotab;
        switch (typenum) {
          case 0 :                                /* Plain mapping */
            o = SCOTCH_graphMapCompute (&grafdat, &mappdat, &stratab[stranum]);
            memcpy (parotab, parttab, vertnbr * sizeof (SCOTCH_Num)); /* Use plain mapping as old mapping in the following */
            break;
          case 1 :                                /* Plain mapping with fixed vertices */
            o = SCOTCH_graphMapFixedCompute (&grafdat, &mappdat, &stratab[stranum]);
            break;
          case 2 :                                /* Remapping without vertex migration load array */
            vmloptr = NULL;
          case 3 :                                /* Remapping with vertex migration load array */
            o = SCOTCH_graphRemapCompute (&grafdat, &mappdat, &mapodat, 0.2, vmloptr, &stratab[stranum]);
            break;
          case 4 :                                /* Remapping with fixed vertices and without vertex migration load array */
            vmloptr = NULL;
          case 5 :                                /* Remapping with fixed vertices and with vertex migration load array */
            o = SCOTCH_graphRemapFixedCompute (&grafdat, &mappdat, &mapodat, 0.2, vmloptr, &stratab[stranum]);
            break;
        }

        if (o != 0) {
          SCOTCH_errorPrint ("main: cannot compute mapping");
          return (1);
        }
      }

      SCOTCH_graphMapSave (&grafdat, &mappdat, fileptr);

      SCOTCH_graphMapExit (&grafdat, &mapodat);
      SCOTCH_graphMapExit (&grafdat, &mappdat);
    }
  }


  for (archnum = 0; archnum < ARCHNBR; archnum ++)
    SCOTCH_archExit (&archtab[archnum]);

  for (stranum = 0; stranum < STRANBR; stranum ++)
    SCOTCH_stratExit (&stratab[stranum]);

  free             (vmlotab);
  free             (parotab);
  free             (parttab);
  SCOTCH_graphExit (&grafdat);

  fclose (fileptr);

  return (0);
}
