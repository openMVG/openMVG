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
/**   NAME       : gtst.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program gives statistics on source **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 31 oct 1994     **/
/**                                 to     03 nov 1994     **/
/**                # Version 3.0  : from : 15 sep 1995     **/
/**                                 to     19 sep 1995     **/
/**                # Version 3.2  : from : 03 jun 1997     **/
/**                                 to   : 25 may 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to   : 19 oct 1998     **/
/**                # Version 3.4  : from : 10 oct 1999     **/
/**                                 to     12 oct 1999     **/
/**                # Version 4.0  : from : 10 sep 2003     **/
/**                                 to   : 10 sep 2003     **/
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

#define GTST

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gtst.h"

/*
**  The static definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {
  "gtst [<input graph file> [<output data file>]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Graph        grafdat;                    /* Source graph */
  SCOTCH_Num          vertnbr;
  SCOTCH_Num          velomin;
  SCOTCH_Num          velomax;
  SCOTCH_Num          velosum;
  double              veloavg;
  double              velodlt;
  SCOTCH_Num          degrmin;
  SCOTCH_Num          degrmax;
  double              degravg;
  double              degrdlt;
  SCOTCH_Num          edgenbr;
  SCOTCH_Num          edlomin;
  SCOTCH_Num          edlomax;
  SCOTCH_Num          edlosum;
  double              edloavg;
  double              edlodlt;
  int                 i;

  errorProg ("gtst");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    exit       (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
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
          fprintf (stderr, "gtst, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          exit       (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit  (&grafdat);
  SCOTCH_graphLoad  (&grafdat, C_filepntrsrcinp, -1, 0);
  SCOTCH_graphCheck (&grafdat);

  SCOTCH_graphSize  (&grafdat, &vertnbr, &edgenbr);
  SCOTCH_graphStat  (&grafdat, &velomin, &velomax, &velosum, &veloavg, &velodlt,
                     &degrmin, &degrmax, &degravg, &degrdlt,
                     &edlomin, &edlomax, &edlosum, &edloavg, &edlodlt);

  if (C_filepntrdatout != NULL) {
    fprintf (C_filepntrdatout, "S\tVertex\tnbr=" SCOTCH_NUMSTRING "\n",
             (SCOTCH_Num) vertnbr);
    fprintf (C_filepntrdatout, "S\tVertex load\tmin=" SCOTCH_NUMSTRING "\tmax=" SCOTCH_NUMSTRING "\tsum=" SCOTCH_NUMSTRING "\tavg=%g\tdlt=%g\n",
             (SCOTCH_Num) velomin, (SCOTCH_Num) velomax, (SCOTCH_Num) velosum, veloavg, velodlt);
    fprintf (C_filepntrdatout, "S\tVertex degree\tmin=" SCOTCH_NUMSTRING "\tmax=" SCOTCH_NUMSTRING "\tsum=" SCOTCH_NUMSTRING "\tavg=%g\tdlt=%g\n",
             (SCOTCH_Num) degrmin, (SCOTCH_Num) degrmax, (SCOTCH_Num) edgenbr, degravg, degrdlt);
    fprintf (C_filepntrdatout, "S\tEdge\tnbr=" SCOTCH_NUMSTRING "\n",
             (SCOTCH_Num) (edgenbr / 2));
    fprintf (C_filepntrdatout, "S\tEdge load\tmin=" SCOTCH_NUMSTRING "\tmax=" SCOTCH_NUMSTRING "\tsum=" SCOTCH_NUMSTRING "\tavg=%g\tdlt=%g\n",
             (SCOTCH_Num) edlomin, (SCOTCH_Num) edlomax, (SCOTCH_Num) edlosum, edloavg, edlodlt);
  }

  SCOTCH_graphExit (&grafdat);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
