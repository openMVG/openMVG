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
/**   NAME       : gmtst.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program computes statistics on     **/
/**                graph mappings.                         **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 17 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 02 jun 1997     **/
/**                                 to   : 16 jul 1997     **/
/**                # Version 3.3  : from : 07 jun 1999     **/
/**                                 to   : 07 jun 1999     **/
/**                # Version 4.0  : from : 12 feb 2004     **/
/**                                 to     16 nov 2005     **/
/**                # Version 5.0  : from : 22 jan 2008     **/
/**                                 to     16 mar 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GMTST

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmtst.h"

/*
**  The static variables.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "r" },
                              { "r" },
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {     /* Usage */
  "gmtst [<input source file> [<input target file> [<input mapping file> [<output data file>]]]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  "",
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
  SCOTCH_Graph        grafdat;                    /* Source graph                      */
  SCOTCH_Num          vertnbr;                    /* Source graph size                 */
  SCOTCH_Num *        vlbltab;                    /* Source graph vertex label array   */
  SCOTCH_Arch         archdat;                    /* Target architecture               */
  SCOTCH_Num          archnbr;                    /* Size of the target architecture   */
  SCOTCH_Mapping      mappdat;                    /* Mapping data                      */
  int                 i;

  errorProg ("gmtst");

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
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'M' :                                /* No mapping flag */
        case 'm' :
          C_filenamemapinp = "-";                 /* Default name to avoid opening   */
          C_filepntrmapinp = NULL;                /* NULL file pointer means no file */
          break;
        case 'V' :
          fprintf (stderr, "gmtst, version " SCOTCH_VERSION_STRING "\n");
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

  SCOTCH_graphInit (&grafdat);                    /* Create graph structure    */
  SCOTCH_graphLoad (&grafdat, C_filepntrsrcinp, -1, 0); /* Read source graph   */
  SCOTCH_graphData (&grafdat, NULL,               /* Get graph characteristics */
                    &vertnbr, NULL, NULL, NULL, &vlbltab,
                    NULL, NULL, NULL);

  SCOTCH_archInit (&archdat);                     /* Create architecture structure                    */
  SCOTCH_archLoad (&archdat, C_filepntrtgtinp);   /* Read target architecture                         */
  if (strcmp (SCOTCH_archName (&archdat), "term") == 0) { /* If target architecture is variable-sized */
    errorPrint ("main: variable-sized architectures cannot be mapped");
    return     (1);
  }

  archnbr = SCOTCH_archSize (&archdat);           /* Get architecture size */

  SCOTCH_graphMapInit (&grafdat, &mappdat, &archdat, NULL); /* Create mapping structure */
  if (SCOTCH_graphMapLoad (&grafdat, &mappdat, C_filepntrmapinp) != 0)
    errorPrint ("main: bad input (1)");

  SCOTCH_graphMapView (&grafdat, &mappdat, C_filepntrdatout); /* Display mapping statistics */

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_graphMapExit (&grafdat, &mappdat);
  SCOTCH_archExit     (&archdat);
  SCOTCH_graphExit    (&grafdat);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
