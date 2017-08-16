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
/**   NAME       : amk_fft2.c                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the distance map for FFT        **/
/**                graphs, to be used to build the archi-  **/
/**                tecture description files for these     **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 19 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to   : 19 sep 1995     **/
/**                # Version 3.2  : from : 07 may 1997     **/
/**                                 to   : 07 may 1997     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 02 oct 1998     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define AMK_FFT2

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "amk_fft2.h"

/*
**  The static and global definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "w" } };

static C_VertDist *         C_distaTab;           /* Pointer to distance map table */
static C_Queue              C_distaQueue;         /* Distance queue                */

static const char *         C_usageList[] = {     /* Usage list */
  "amk_fft2 [<dim> [<output target file>]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/*************************************************/
/*                                               */
/* The main routine, which computes the distance */
/* triangular table.                             */
/*                                               */
/*************************************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Num          fdim;                       /* FFT dimension               */
  SCOTCH_Num          fnbr;                       /* Number of FFT vertices      */
  SCOTCH_Num          fmax;                       /* Maximum terminal number     */
  SCOTCH_Num          fmsk;                       /* Position bit mask           */
  C_Vertex            v, w, x;                    /* A FFT vertex (lvl, pos)     */
  SCOTCH_Num          b, d;                       /* Mask and bit variables      */
  SCOTCH_Num          i;                          /* Loop counter                */
  SCOTCH_Num          t;                          /* Vertex terminal value       */

  errorProg ("amk_fft2");

  fdim = 2;

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 1) {                        /* If number of parameters not reached              */
        if ((fdim = atoi (argv[i])) < 1) {        /* Get the dimension                                */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
          return     (1);
        }
        C_paraNum ++;
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILENBR)                  /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else {
        errorPrint ("main: too many file names given");
        return     (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "amk_fft2, version " SCOTCH_VERSION_STRING "\n");
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

  fnbr = (fdim + 1) * (1 << fdim);                /* Compute number of vertices      */
  fmax = (1 << (fdim * 2 - 1)) - 1;               /* Compute maximum terminal number */
  fmsk = (1 << fdim) - 1;                         /* Get maximum position number     */

  fprintf (C_filepntrarcout, "deco\n0\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n", /* Print file header */
           (SCOTCH_Num) fnbr,                     /* Print number of terminal domains */
           (SCOTCH_Num) fmax);                    /* Print the biggest terminal value */

  for (v.lvl = 0; v.lvl <= fdim; v.lvl ++) {      /* For all vertices */
    for (v.pos = 0; v.pos <= fmsk; v.pos ++) {
      for (i = v.lvl, b = 1 << (fdim - 1), t = 1; /* Recurse through the vertical + horizontal cuts */
           i <= fdim;
           i ++, b >>= 1) {
        t <<= 1;                                  /* Vertical cut: tell if vertex is in left or...   */
        t |= (v.pos & b) ? 1 : 0;                 /* right part from the position heaviest bits      */
        t <<= 1;                                  /* Vertex is still in upper part of horizontal cut */
      }
      if (v.lvl == 0)                             /* If vertex is in the first level...              */
        t >>= 2;                                  /* We have gone one step too far                   */
      else {                                      /* Else                                            */
        t |= 1;                                   /* This time vertex is in the lower part           */
        t <<= (v.lvl - 1);                        /* Make space for the chain bipartition            */
        t |= v.pos & ((1 << (v.lvl - 1)) - 1);    /* Bipartition the chain following the lowest bits */
      }

      printf (((v.lvl == fdim) && (v.pos == fmsk)) /* Print terminal domain number */
               ? SCOTCH_NUMSTRING "\n\n" : SCOTCH_NUMSTRING " ", t);
    }
  }

  if ((C_queueInit (&C_distaQueue, fnbr) != 0) || /* Allocate distance array */
      ((C_distaTab = (C_VertDist *) memAlloc (fnbr * sizeof (C_VertDist))) == NULL)) {
    errorPrint ("main: out of memory");
    return     (1);
  }

  for (v.lvl = 0; v.lvl <= fdim; v.lvl ++) {      /* For all vertices */
    for (v.pos = 0; v.pos <= fmsk; v.pos ++) {
      for (i = 0; i < fnbr; i ++)                 /* Initialize the vertex table */
        C_distaTab[i].queued = 0;                 /* Vertex not queued yet       */

      C_distaRoot (&v);                           /* Set the queue with root v */

      while (C_distaGet (&w, &d)) {               /* As long as the queue is not empty */
        C_distaTab[C_vertLabl (&w)].dist = d;     /* Keep the distance information     */

        d ++;                                     /* Search for neighbors at next level */
        if (w.lvl > 0) {                          /* Add new neighbors to the queue     */
          x.lvl = w.lvl - 1;
          x.pos = w.pos;
          C_distaPut (&x, d);
          x.pos = w.pos ^ (1 << (w.lvl - 1));
          C_distaPut (&x, d);
        }
        if (w.lvl < fdim) {
          x.lvl = w.lvl + 1;
          x.pos = w.pos;
          C_distaPut (&x, d);
          x.pos = w.pos ^ (1 << w.lvl);
          C_distaPut (&x, d);
        }
      }

      if (v.lvl + v.pos > 0) {                    /* Print the distance triangular map line */
        fprintf (C_filepntrarcout, SCOTCH_NUMSTRING,
                 (SCOTCH_Num) C_distaTab[0].dist);
        for (i = 1; i < (v.lvl << fdim) + v.pos; i ++)
          fprintf (C_filepntrarcout, " " SCOTCH_NUMSTRING,
                   (SCOTCH_Num) C_distaTab[i].dist);
        fprintf (C_filepntrarcout, "\n");
      }
    }
  }

  C_queueExit (&C_distaQueue);
  memFree     (C_distaTab);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
