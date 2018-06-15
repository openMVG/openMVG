/* Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : mmk_m3.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the source meshes for           **/
/**                tridimensional mesh source graphs.      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 26 sep 2002     **/
/**                                 to   : 17 feb 2004     **/
/**                # Version 5.0  : from : 13 dec 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/**   NOTES      : # The nodes and elements of the         **/
/**                  (dX,dY,dZ) mesh are numbered so that  **/
/**                  t(0,0,0) = 0, t(1,0,0) = 1,           **/
/**                  t(dX - 1, 0, 0) = dX - 1, t(0,1,0) =  **/
/**                  dX, t (0, 0, 1) = dX * dY - 1,        **/
/**                  and t(x,y,z) = (z * dY + y) * dX + x. **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MMK_M3

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "mmk_m3.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "w" },
                              { "w" } };

static const int            C_nghbTab[4] = { 8, 4, 2, 1 };

static const char *         C_usageList[] = {
  "mmk_m3 <dimX> [<dimY> [<dimZ> [<output mesh file>]]] <options>",
  "  -g<file>  : Output mesh geometry to <file>",
  "  -h        : Display this help",
  "  -V        : Print program version and copyright",
  NULL };

/****************************************/
/*                                      */
/* The main routine, which computes the */
/* source graph description.            */
/*                                      */
/****************************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Num          e[3] = { 1, 1, 1 };         /* Mesh element dimensions */
  SCOTCH_Num          n[3];                       /* Mesh node dimensions    */
  SCOTCH_Num          c[3];                       /* Vertex coordinates      */
  SCOTCH_Num          velmnbr;                    /* First node number       */
  int                 flagval;                    /* Process flags           */
  int                 i;

  errorProg ("mmk_m3");

  flagval = C_FLAGDEFAULT;                        /* Set default flags */

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 3) {                        /* If number of parameters not reached              */
        if ((e[C_paraNum ++] = atoi (argv[i])) < 1) { /* Get the dimension                            */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
          return     (1);
        }
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else {
        errorPrint ("main: too many file names given");
        return     (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'G' :                                /* Output mesh geometry */
        case 'g' :
          flagval |= C_FLAGGEOOUT;
          if (argv[i][2] != '\0')
            C_filenamegeoout = &argv[i][2];
          break;
        case 'H' :                                /* Give program usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "mmk_m3, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          return     (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  n[0] = e[0] + 1;
  n[1] = e[1] + 1;
  n[2] = e[2] + 1;
  velmnbr = e[0] * e[1] * e[2];

  fprintf (C_filepntrmshout, "1\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n0\t" SCOTCH_NUMSTRING "\t000\n", /* Print mesh file header */
             (SCOTCH_Num) velmnbr,
             (SCOTCH_Num) (n[0] * n[1] * n[2]),
             (SCOTCH_Num) ((velmnbr + (n[0] * n[1] * n[2]) -
                                      (n[0] * n[1] + n[0] * n[2] + n[1] * n[2]) +
                                       n[0] + n[1] + n[2] - 1) * 8),
             (SCOTCH_Num) velmnbr);

  for (c[2] = 0; c[2] < e[2]; c[2] ++) {          /* Output element neighbor list */
    for (c[1] = 0; c[1] < e[1]; c[1] ++) {
      for (c[0] = 0; c[0] < e[0]; c[0] ++)
        fprintf (C_filepntrmshout, "8\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n", /* Output neighbors of element */
                  (SCOTCH_Num) ((c[2] * n[1] + c[1]) * n[0] + c[0]),
                  (SCOTCH_Num) ((c[2] * n[1] + c[1]) * n[0] + (c[0] + 1)),
                  (SCOTCH_Num) ((c[2] * n[1] + (c[1] + 1)) * n[0] + c[0]),
                  (SCOTCH_Num) ((c[2] * n[1] + (c[1] + 1)) * n[0] + (c[0] + 1)),
                  (SCOTCH_Num) (((c[2] + 1) * n[1] + c[1]) * n[0] + c[0]),
                  (SCOTCH_Num) (((c[2] + 1) * n[1] + c[1]) * n[0] + (c[0] + 1)),
                  (SCOTCH_Num) (((c[2] + 1) * n[1] + (c[1] + 1)) * n[0] + c[0]),
                  (SCOTCH_Num) (((c[2] + 1) * n[1] + (c[1] + 1)) * n[0] + (c[0] + 1)));
    }
  }
  for (c[2] = 0; c[2] < n[2]; c[2] ++) {          /* Output node neighbor list */
    for (c[1] = 0; c[1] < n[1]; c[1] ++) {
      for (c[0] = 0; c[0] < n[0]; c[0] ++) {
        fprintf (C_filepntrmshout, "%d",          /* Output number of neighboring elements */
                  C_nghbTab[(((c[0] != 0) && (c[0] != e[0])) ? 0 : 1) +
                            (((c[1] != 0) && (c[1] != e[1])) ? 0 : 1) +
                            (((c[2] != 0) && (c[2] != e[2])) ? 0 : 1)]);
        if (c[2] != 0) {                          /* Output neighbors of nodes */
          if (c[1] != 0) {
            if (c[0] != 0)
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) (((c[2] - 1) * e[1] + (c[1] - 1)) * e[0] + (c[0] - 1)));
            if (c[0] != e[0])
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) (((c[2] - 1) * e[1] + (c[1] - 1)) * e[0] + c[0]));
          }
          if (c[1] != e[1]) {
            if (c[0] != 0)
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) (((c[2] - 1) * e[1] + c[1]) * e[0] + (c[0] - 1)));
            if (c[0] != e[0])
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) (((c[2] - 1) * e[1] + c[1]) * e[0] + c[0]));
          }
        }
        if (c[2] != e[2]) {
          if (c[1] != 0) {
            if (c[0] != 0)
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) ((c[2] * e[1] + (c[1] - 1)) * e[0] + (c[0] - 1)));
            if (c[0] != e[0])
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) ((c[2] * e[1] + (c[1] - 1)) * e[0] + c[0]));
          }
          if (c[1] != e[1]) {
            if (c[0] != 0)
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) ((c[2] * e[1] + c[1]) * e[0] + (c[0] - 1)));
            if (c[0] != e[0])
              fprintf (C_filepntrmshout, "\t" SCOTCH_NUMSTRING,
                        (SCOTCH_Num) ((c[2] * e[1] + c[1]) * e[0] + c[0]));
          }
        }
        fprintf (C_filepntrmshout, "\n");
      }
    }
  }

  if (flagval & C_FLAGGEOOUT) {                   /* If geometry is wanted                 */
    fprintf (C_filepntrgeoout, "3\n" SCOTCH_NUMSTRING "\n", /* Output geometry file header */
              (SCOTCH_Num) (velmnbr + n[0] * n[1] * n[2]));

    for (c[2] = 0; c[2] < e[2]; c[2] ++) {        /* Output element coordinates */
      for (c[1] = 0; c[1] < e[1]; c[1] ++) {
        for (c[0] = 0; c[0] < e[0]; c[0] ++)
          fprintf (C_filepntrgeoout, SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING ".5\t" SCOTCH_NUMSTRING ".5\t" SCOTCH_NUMSTRING ".5\n",
                    (SCOTCH_Num) (((c[2] * e[1]) + c[1]) * e[0] + c[0]),
                    (SCOTCH_Num) c[0],
                    (SCOTCH_Num) (e[1] - 1 - c[1]),
                    (SCOTCH_Num) c[2]);
      }
    }
    for (c[2] = 0; c[2] <= e[2]; c[2] ++) {       /* Output node coordinates */
      for (c[1] = 0; c[1] <= e[1]; c[1] ++) {
        for (c[0] = 0; c[0] <= e[0]; c[0] ++)
          fprintf (C_filepntrgeoout, SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n",
                    (SCOTCH_Num) (velmnbr + ((c[2] * n[1]) + c[1]) * n[0] + c[0]),
                    (SCOTCH_Num) c[0],
                    (SCOTCH_Num) (e[1] - c[1]),
                    (SCOTCH_Num) c[2]);
      }
    }
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
