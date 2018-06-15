/* Copyright 2004,2007,2008,2010-2012,2014 Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gmk_hy.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the source graph for hypercube  **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 03 nov 1994     **/
/**                                 to     03 nov 1994     **/
/**                # Version 3.0  : from : 11 jul 1995     **/
/**                                 to     11 jul 1995     **/
/**                # Version 3.2  : from : 03 jun 1997     **/
/**                                 to   : 03 jun 1997     **/
/**                # Version 3.3  : from : 06 oct 1998     **/
/**                                 to   : 06 oct 1998     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 16 map 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GMK_HY

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmk_hy.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "w" } };

static const char *         C_usageList[] = {
  "gmk_hy <dim> [<output source file>] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
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
  SCOTCH_Num          hdim = 1;                   /* Graph dimension      */
  SCOTCH_Num          hnbr;                       /* Number of vertices   */
  SCOTCH_Num          hbit;                       /* Most significant bit */
  SCOTCH_Num          hvrt;
  SCOTCH_Num          hngb;
  int                 i;

  errorProg ("gmk_hy");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    exit       (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 1) {                        /* If number of parameters not reached              */
        if ((hdim = atoi (argv[i])) < 1) {        /* Get the dimension                                */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
          exit       (1);
        }
        C_paraNum ++;
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else {
        errorPrint ("main: too many file names given");
        exit       (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          exit       (0);
        case 'V' :
          fprintf (stderr, "gmk_hy, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          exit       (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  hnbr = 1 <<  hdim;                              /* Compute number of vertices */
  hbit = 1 << (hdim - 1);                         /* Compute highest bit value  */

  fprintf (C_filepntrsrcout, "0\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n0\t000\n",
           (SCOTCH_Num) hnbr,                     /* Print number of vertices     */
           (SCOTCH_Num) (hdim * hnbr));           /* Print number of edges (arcs) */

  for (hvrt = 0; hvrt < hnbr; hvrt ++) {          /* For all vertices */
    fprintf (C_filepntrsrcout, "" SCOTCH_NUMSTRING "",
             (SCOTCH_Num) hdim);                  /* Output number of neighbors         */
    for (hngb = hbit; hngb > 0; hngb >>= 1)       /* For all vertex bits                */
      fprintf (C_filepntrsrcout, "\t" SCOTCH_NUMSTRING, /* Write corresponding neighbor */
               (SCOTCH_Num) (hvrt ^ hngb));
    fprintf (C_filepntrsrcout, "\n");
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
