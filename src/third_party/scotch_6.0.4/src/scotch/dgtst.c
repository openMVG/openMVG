/* Copyright 2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgtst.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program gives statistics on        **/
/**                distributed source graphs.              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 jun 2007     **/
/**                                 to   : 16 jun 2008     **/
/**                # Version 5.1  : from : 26 oct 2008     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGTST
#define SCOTCH_PTSCOTCH

#include "module.h"
#include "common.h"
#include "ptscotch.h"
#include "dgtst.h"

/*
**  The static and global definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {
  "dgtst [<input graph file> [<output data file>]] <options>",
  "  -h       : Display this help",
  "  -r<num>  : Set root process for centralized files (default is 0)",
  "  -V       : Print program version and copyright",
  NULL };

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
  SCOTCH_Dgraph       grafdat;
  int                 procglbnbr;
  int                 proclocnum;
  int                 protglbnum;                 /* Root process */
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
  int                 flagval;
  int                 i;
#ifdef SCOTCH_PTHREAD
  int                 thrdlvlreqval;
  int                 thrdlvlproval;
#endif /* SCOTCH_PTHREAD */

  errorProg ("dgtst");

#ifdef SCOTCH_PTHREAD
  thrdlvlreqval = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdlvlreqval, &thrdlvlproval) != MPI_SUCCESS)
    errorPrint ("main: Cannot initialize (1)");
  if (thrdlvlreqval > thrdlvlproval)
    errorPrint ("main: MPI implementation is not thread-safe: recompile without SCOTCH_PTHREAD");
#else /* SCOTCH_PTHREAD */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    errorPrint ("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD */

  MPI_Comm_size (MPI_COMM_WORLD, &procglbnbr);    /* Get communicator data */
  MPI_Comm_rank (MPI_COMM_WORLD, &proclocnum);
  protglbnum = 0;                                 /* Assume root process is process 0 */

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  flagval = C_FLAGNONE;

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
#ifdef SCOTCH_DEBUG_ALL
        case 'D' :
        case 'd' :
          flagval |= C_FLAGDEBUG;
          break;
#endif /* SCOTCH_DEBUG_ALL */
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'R' :                                /* Root process (if necessary) */
        case 'r' :
          protglbnum = atoi (&argv[i][2]);
          if ((protglbnum < 0)           ||
              (protglbnum >= procglbnbr) ||
              ((protglbnum == 0) && (argv[i][2] != '0')))
            errorPrint ("main: invalid root process number");
          break;
        case 'V' :
        case 'v' :
          fprintf (stderr, "dgtst, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

#ifdef SCOTCH_DEBUG_ALL
  if ((flagval & C_FLAGDEBUG) != 0) {
    fprintf (stderr, "Proc %4d of %d, pid %d\n", proclocnum, procglbnbr, getpid ());
    if (proclocnum == protglbnum) {               /* Synchronize on keybord input */
      char           c;

      printf ("Waiting for key press...\n");
      scanf ("%c", &c);
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
#endif /* SCOTCH_DEBUG_ALL */

  fileBlockOpenDist (C_fileTab, C_FILENBR, procglbnbr, proclocnum, protglbnum); /* Open all files */

  SCOTCH_dgraphInit  (&grafdat, MPI_COMM_WORLD);
  SCOTCH_dgraphLoad  (&grafdat, C_filepntrsrcinp, -1, 0);
  SCOTCH_dgraphCheck (&grafdat);

  SCOTCH_dgraphSize (&grafdat, &vertnbr, NULL, &edgenbr, NULL);
  SCOTCH_dgraphStat (&grafdat, &velomin, &velomax, &velosum, &veloavg, &velodlt,
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

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_dgraphExit (&grafdat);

  MPI_Finalize ();
#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
