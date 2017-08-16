/* Copyright 2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dggath.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program gathers the fragments of a **/
/**                distributed graph and saves it as a     **/
/**                centralized source graph.               **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 26 oct 2008     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGGATH
#define SCOTCH_PTSCOTCH

#include "module.h"
#include "common.h"
#include "ptscotch.h"
#include "dggath.h"

/*
**  The static and global definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {
  "dggath [<input (distributed) source file> [<output centralized source file>]] <options>",
  "  -c       : Check the input graph after loading",
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
  SCOTCH_Graph *      cgrfptr;
  SCOTCH_Graph        cgrfdat;
  SCOTCH_Dgraph       dgrfdat;
  int                 procglbnbr;
  int                 proclocnum;
  int                 protglbnum;                 /* Root process */
  int                 flagval;
  int                 i;
  int                 reduloctab[2];
  int                 reduglbtab[2];
#ifdef SCOTCH_PTHREAD
  int                 thrdlvlreqval;
  int                 thrdlvlproval;
#endif /* SCOTCH_PTHREAD */

  errorProg ("dggath");

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

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'C' :
        case 'c' :
          flagval |= C_FLAGCHECK;
          break;
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
              ((protglbnum == 0) && (argv[i][2] != '0'))) {
            errorPrint ("main: invalid root process number");
          }
          break;
        case 'V' :
        case 'v' :
          fprintf (stderr, "dggath, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
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
      scanf  ("%c", &c);
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
#endif /* SCOTCH_DEBUG_ALL */

  fileBlockOpenDist (C_fileTab, C_FILENBR, procglbnbr, proclocnum, protglbnum); /* Open all files */

  if (C_filepntrsrcout == NULL) {
    cgrfptr = NULL;
    reduloctab[0] =
    reduloctab[1] = 0;
  }
  else {
    cgrfptr = &cgrfdat;
    reduloctab[0] = 1;
    reduloctab[1] = proclocnum;
  }
  if (MPI_Allreduce (reduloctab, reduglbtab, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    errorPrint ("main: communication error");

  if (reduglbtab[0] != 1)
    errorPrint ("main: should have only one root");
  if (reduglbtab[1] != protglbnum)
    errorPrint ("main: root process mismatch");

  SCOTCH_dgraphInit (&dgrfdat, MPI_COMM_WORLD);
  SCOTCH_dgraphLoad (&dgrfdat, C_filepntrsrcinp, -1, 0);

  if ((flagval & C_FLAGCHECK) != 0)
    SCOTCH_dgraphCheck (&dgrfdat);

  SCOTCH_graphInit    (&cgrfdat);
  SCOTCH_dgraphGather (&dgrfdat, cgrfptr);
  if (cgrfptr != NULL)
    SCOTCH_graphSave (cgrfptr, C_filepntrsrcout);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_graphExit  (&cgrfdat);
  SCOTCH_dgraphExit (&dgrfdat);

  MPI_Finalize ();
#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
