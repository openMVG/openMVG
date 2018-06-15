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
/**   NAME       : mord.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a mesh reordering software.     **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 15 nov 2002     **/
/**                                 to   : 06 jan 2006     **/
/**                # Version 5.0  : from : 22 jan 2008     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 08 sep 2008     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MORD

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "mord.h"

/*
**  The static and global definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" },
                              { "w" },
                              { "w" } };

static const char *         C_usageList[] = {
  "mord [<input source mesh file> [<output ordering file> [<output log file>]]] <options>",
  "  -c<opt>    : Choose default ordering strategy according to one or several of <opt>:",
  "                 b  : enforce load balance as much as possible",
  "                 q  : privilege quality over speed (default)",
  "                 s  : privilege speed over quality",
  "                 t  : enforce safety",
  "  -h         : Display this help",
  "  -m<file>   : Save column block mapping data to <file>",
  "  -o<strat>  : Use mesh ordering strategy <strat> (see user's manual)",
  "               (see default strategy with option '-vs')",
  "  -t<file>   : Save partitioning tree data to <file>",
  "  -V         : Print program version and copyright",
  "  -v<verb>   : Set verbose mode to <verb> :",
  "                 s  : strategy information",
  "                 t  : timing information",
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
  SCOTCH_Num          vnodnbr;                    /* Number of nodes   */
  SCOTCH_Mesh         meshdat;                    /* Source graph      */
  SCOTCH_Ordering     ordedat;                    /* Graph ordering    */
  SCOTCH_Num *        permtab;                    /* Permutation array */
  SCOTCH_Strat        stradat;                    /* Ordering strategy */
  SCOTCH_Num          straval;
  char *              straptr;
  int                 flagval;
  Clock               runtime[2];                 /* Timing variables  */
  int                 i, j;

  errorProg ("mord");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  flagval = C_FLAGNONE;                           /* Default behavior  */
  straval = 0;                                    /* No strategy flags */
  straptr = NULL;
  SCOTCH_stratInit (&stradat);

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else {
        errorPrint ("main: too many file names given");
        return     (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'C' :
        case 'c' :                                /* Strategy selection parameters */
          for (j = 2; argv[i][j] != '\0'; j ++) {
            switch (argv[i][j]) {
              case 'B' :
              case 'b' :
                straval |= SCOTCH_STRATBALANCE;
                break;
              case 'Q' :
              case 'q' :
                straval |= SCOTCH_STRATQUALITY;
                break;
              case 'S' :
              case 's' :
                straval |= SCOTCH_STRATSPEED;
                break;
              case 'T' :
              case 't' :
                straval |= SCOTCH_STRATSAFETY;
                break;
              default :
                errorPrint ("main: invalid strategy selection option '%c'", argv[i][j]);
            }
          }
          break;
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'M' :                                /* Output separator mapping */
        case 'm' :
          flagval |= C_FLAGMAPOUT;
          if (argv[i][2] != '\0')
            C_filenamemapout = &argv[i][2];
          break;
        case 'O' :                                /* Ordering strategy */
        case 'o' :
          straptr = &argv[i][2];
          SCOTCH_stratExit (&stradat);
          SCOTCH_stratInit (&stradat);
          if ((SCOTCH_stratMeshOrder (&stradat, straptr)) != 0) {
            errorPrint ("main: invalid ordering strategy");
            return     (1);
          }
          break;
        case 'V' :
          fprintf (stderr, "mord, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
          break;
        case 'v' :                               /* Output control info */
          for (j = 2; argv[i][j] != '\0'; j ++) {
            switch (argv[i][j]) {
              case 'S' :
              case 's' :
                flagval |= C_FLAGVERBSTR;
                break;
              case 'T' :
              case 't' :
                flagval |= C_FLAGVERBTIM;
                break;
              default :
                errorPrint ("main: unprocessed parameter '%c' in '%s'",
                            argv[i][j], argv[i]);
                return     (1);
            }
          }
          break;
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          return     (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  clockInit  (&runtime[0]);
  clockStart (&runtime[0]);

  SCOTCH_meshInit (&meshdat);                     /* Create mesh structure */
  SCOTCH_meshLoad (&meshdat, C_filepntrsrcinp, -1); /* Read source mesh    */
  SCOTCH_meshSize (&meshdat, NULL, &vnodnbr, NULL); /* Get number of nodes */

  if (straval != 0) {
    if (straptr != NULL)
      errorPrint ("main: options '-c' and '-o' are exclusive");

    SCOTCH_stratMeshOrderBuild (&stradat, straval, 0.1);
  }

  clockStop  (&runtime[0]);                       /* Get input time */
  clockInit  (&runtime[1]);
  clockStart (&runtime[1]);

  if ((permtab = (SCOTCH_Num *) memAlloc (vnodnbr * sizeof (SCOTCH_Num))) == NULL) {
    errorPrint ("main: out of memory");
    return     (1);
  }
  SCOTCH_meshOrderInit    (&meshdat, &ordedat, permtab, NULL, NULL, NULL, NULL); /* Create ordering */
  SCOTCH_meshOrderCompute (&meshdat, &ordedat, &stradat); /* Perform ordering */

  clockStop (&runtime[1]);                        /* Get ordering time */

#ifdef SCOTCH_DEBUG_ALL
  if (SCOTCH_meshOrderCheck (&meshdat, &ordedat) != 0)
    return (1);
#endif /* SCOTCH_DEBUG_ALL */

  clockStart (&runtime[0]);

  SCOTCH_meshOrderSave (&meshdat, &ordedat, C_filepntrordout); /* Write ordering     */
  if (flagval & C_FLAGMAPOUT)                     /* If mapping wanted               */
    SCOTCH_meshOrderSaveMap (&meshdat, &ordedat, C_filepntrmapout); /* Write mapping */

  clockStop  (&runtime[0]);                       /* Get output time */

  if (flagval & C_FLAGVERBSTR) {
    fprintf (C_filepntrlogout, "S\tStrat=");
    SCOTCH_stratSave (&stradat, C_filepntrlogout);
    putc ('\n', C_filepntrlogout);
  }
  if (flagval & C_FLAGVERBTIM) {
    fprintf (C_filepntrlogout, "T\tOrder\t\t%g\nT\tI/O\t\t%g\nT\tTotal\t\t%g\n",
             (double) clockVal (&runtime[1]),
             (double) clockVal (&runtime[0]),
             (double) clockVal (&runtime[0]) +
             (double) clockVal (&runtime[1]));
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_meshOrderExit (&meshdat, &ordedat);
  SCOTCH_stratExit     (&stradat);
  SCOTCH_meshExit      (&meshdat);
  memFree              (permtab);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}

/*******************************************/
/*                                         */
/* Stubs to avoid including target module. */
/*                                         */
/*******************************************/

SCOTCH_Num
T_domWght (
const void * const          arch,
const void * const          dom)
{
  errorPrint ("T_domWghtMord: internal error");
  return     (1);
}

SCOTCH_Num
_SCOTCHTdomWght (
const void * const          arch,
const void * const          dom)
{
  errorPrint ("T_domWghtMord: internal error");
  return     (1);
}

SCOTCH_Num
T_domDist (
const void * const          arch,
const void * const          dom0,
const void * const          dom1)
{
  errorPrint ("T_domDistMord: internal error");
  return     (1);
}

SCOTCH_Num
_SCOTCHTdomDist (
const void * const          arch,
const void * const          dom0,
const void * const          dom1)
{
  errorPrint ("T_domDistMord: internal error");
  return     (1);
}
