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
/**   NAME       : gmk_msh.c                               **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a mesh-to-graph converter.      **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 21 jan 2004     **/
/**                                 to   : 21 jan 2004     **/
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

#define GMK_MSH

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmk_msh.h"

/*
**  The static and global variables.
*/

static int                  C_fileNum    = 0;     /* Number of file in arg list  */
static File                 C_fileTab[2] = {      /* File array                  */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {
  "gmk_msh [<input mesh file> [<output graph file>]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/*****************************/
/*                           */
/* This is the main function */
/*                           */
/*****************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Mesh         meshdat;
  SCOTCH_Graph        grafdat;
  int                 i;

  errorProg ("gmk_msh");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

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
    else {                                       /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                               /* Give help */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "gmk_msh, version " SCOTCH_VERSION_STRING "\n");
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

  SCOTCH_meshInit  (&meshdat);
  SCOTCH_graphInit (&grafdat);

  SCOTCH_meshLoad  (&meshdat, C_filepntrmshinp, -1);
  SCOTCH_meshCheck (&meshdat);
  SCOTCH_meshGraph (&meshdat, &grafdat);
#ifdef SCOTCH_DEBUG_ALL
  if (SCOTCH_graphCheck (&grafdat) != 0) {
    errorPrint ("main: bad graph structure");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ALL */
  SCOTCH_graphSave (&grafdat, C_filepntrgrfout);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_graphExit (&grafdat);
  SCOTCH_meshExit  (&meshdat);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
