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
/**   NAME       : gmk_ub2.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the source graph for undirected **/
/**                de Bruijn graphs, to be used to build   **/
/**                the architecture description files for  **/
/**                these graphs.                           **/
/**                                                        **/
/**   DATES      : # Version 1.2  : from : 11 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 2.0  : from : 05 nov 1994     **/
/**                                 to     05 nov 1994     **/
/**                # Version 3.0  : from : 11 jul 1995     **/
/**                                 to     12 jul 1995     **/
/**                # Version 3.2  : from : 03 jun 1997     **/
/**                                 to   : 03 jun 1997     **/
/**                # Version 3.3  : from : 07 jun 1999     **/
/**                                 to   : 07 jun 1999     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 5.0  : from : 22 jan 2008     **/
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

#define GMK_UB2

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmk_ub2.h"

#define ngbadd(v)                 if ((v) != vertnum) {           \
                                    int                   k;      \
                                    for (k = 0; k < ngbnbr; k ++) \
                                      if ((v) == ngbtab[k])       \
                                        break;                    \
                                    if (k == ngbnbr)              \
                                      ngbtab[ngbnbr ++] = (v);    \
                                  }

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "w" } };

static const char *         C_usageList[] = {
  "gmk_ub2 <dim> [<output source file>] <options>",
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
  SCOTCH_Num          ubdim = 1;                  /* Graph dimension             */
  SCOTCH_Num          ubnbr;                      /* Number of vertices          */
  SCOTCH_Num          ubbit;                      /* Most significant bit        */
  SCOTCH_Num          ngbtab[4];                  /* Array of neighbors          */
  int                 ngbnbr;                     /* Current number of neighbors */
  SCOTCH_Num          vertnum;                    /* Current vertex number       */
  int                 i, j;

  errorProg ("gmk_ub2");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 1) {                        /* If number of parameters not reached              */
        if ((ubdim = (SCOTCH_Num) atol (argv[i])) < 1) { /* Get dimension                             */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
          return     (1);
        }
        C_paraNum ++;
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
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "gmk_ub2, version " SCOTCH_VERSION_STRING "\n");
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

  ubnbr = 1 <<  ubdim;                            /* Compute number of vertices */
  ubbit = 1 << (ubdim - 1);                       /* Bit to add on the left     */

  fprintf (C_filepntrsrcout, "0\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n0\t000\n",
           (SCOTCH_Num) ubnbr,                    /* Print number of vertices     */
           (SCOTCH_Num) (4 * ubnbr - 6));         /* Print number of edges (arcs) */

  for (vertnum = 0; vertnum < ubnbr; vertnum ++) { /* For all vertices         */
    ngbnbr = 0;                                   /* No neighbors defined yet  */
    ngbadd  ((vertnum << 1) & (ubnbr - 1));       /* Register vertex neighbors */
    ngbadd (((vertnum << 1) & (ubnbr - 1)) | 1);
    ngbadd  ((vertnum >> 1) & (ubnbr - 1));
    ngbadd (((vertnum >> 1) & (ubnbr - 1)) | ubbit);

    fprintf (C_filepntrsrcout, "%d", ngbnbr);     /* Output number of neighbors */
    for (j = 0; j < ngbnbr; j ++)
      fprintf (C_filepntrsrcout, "\t" SCOTCH_NUMSTRING,
               (SCOTCH_Num) ngbtab[j]);
    fprintf (C_filepntrsrcout, "\n");
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
