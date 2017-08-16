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
/**   NAME       : gcv.c                                   **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Bruno MARCUSSEAU (v3.1)                 **/
/**                                                        **/
/**   FUNCTION   : Part of a graph file converter.         **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 apr 1993     **/
/**                                 to     02 apr 1993     **/
/**                # Version 2.0  : from : 28 oct 1994     **/
/**                                 to     16 nov 1994     **/
/**                # Version 3.0  : from : 08 sep 1995     **/
/**                                 to     20 sep 1995     **/
/**                # Version 3.1  : from : 22 mar 1996     **/
/**                                 to     22 mar 1996     **/
/**                # Version 3.2  : from : 04 oct 1996     **/
/**                                 to     26 may 1997     **/
/**                # Version 3.3  : from : 06 oct 1998     **/
/**                                 to   : 21 dec 1998     **/
/**                # Version 3.4  : from : 05 oct 1999     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to   : 19 jan 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 11 jun 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GCV

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gcv.h"

/*
**  The static and global variables.
*/

static int                  C_inpFormatType  = 0; /* Input graph format           */
static char *               C_inpFormatData  = "\0"; /* Pointer to auxiliary data */
static const C_Format       C_inpFormatTab[] = {  /* Table of input formats       */
                              { 'B',  SCOTCH_graphGeomLoadHabo },
                              { 'b',  SCOTCH_graphGeomLoadHabo },
                              { 'C',  SCOTCH_graphGeomLoadChac },
                              { 'c',  SCOTCH_graphGeomLoadChac },
                              { 'M',  SCOTCH_graphGeomLoadMmkt },
                              { 'm',  SCOTCH_graphGeomLoadMmkt },
                              { 'S',  SCOTCH_graphGeomLoadScot },
                              { 's',  SCOTCH_graphGeomLoadScot },
                              { '\0', NULL } };

static int                  C_outFormatType  = 4; /* Output graph format          */
static char *               C_outFormatData  = "\0"; /* Pointer to auxiliary data */
static C_Format             C_outFormatTab[] = {  /* Table of output formats      */
                              { 'C',  SCOTCH_graphGeomSaveChac },
                              { 'c',  SCOTCH_graphGeomSaveChac },
                              { 'M',  SCOTCH_graphGeomSaveMmkt },
                              { 'm',  SCOTCH_graphGeomSaveMmkt },
                              { 'S',  SCOTCH_graphGeomSaveScot },
                              { 's',  SCOTCH_graphGeomSaveScot },
                              { '\0', NULL } };

static int                  C_fileNum    = 0;     /* Number of file in arg list  */
static File                 C_fileTab[3] = {      /* File array                  */
                              { "r" },
                              { "w" },
                              { "w" } };

static const char *         C_usageList[] = {
  "gcv [<input graph file> [<output graph file> [<output geometry file>]]] <options>",
  "  -h          : Display this help",
  "  -i<format>  : Select input file format",
  "                  b  : Boeing-Harwell format (matrices)",
  "                  c  : Chaco v2.0 format (adjacency)",
  "                  m  : Matrix Market format (edges, symmetrized)",
  "                  s  : Scotch v3.0 format (adjacency)",
  "  -o<format>  : Select output file format",
  "                  c  : Chaco v2.0 format (adjacency)",
  "                  m  : Matrix Market symmetric pattern format (edges)",
  "                  s  : Scotch v3.0 format (adjacency)",
  "  -V          : Print program version and copyright",
  "",
  "Default option set is : '-Ib -Os'",
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
  SCOTCH_Graph        grafdat;
  SCOTCH_Geom         geomdat;
  int                 i, j;

  errorProg ("gcv");

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
        case 'I' :                               /* Select input file type */
        case 'i' :
          for (j = 0; C_inpFormatTab[j].code != '\0'; j ++) { /* Find proper format code */
            if (C_inpFormatTab[j].code == argv[i][2]) {
              C_inpFormatType = j;
              C_inpFormatData = &argv[i][3];
              break;
            }
          }
          if (C_inpFormatTab[j].code == '\0') {
            errorPrint ("main: unprocessed option '%s'", argv[i]);
            return     (1);
          }
          break;
        case 'O' :                               /* Select input file type */
        case 'o' :
          for (j = 0; C_outFormatTab[j].code != '\0'; j ++) { /* Find proper format code */
            if (C_outFormatTab[j].code == argv[i][2]) {
              C_outFormatType = j;
              C_outFormatData = &argv[i][3];
              break;
            }
          }
          if (C_inpFormatTab[j].code == '\0') {
            errorPrint ("main: unprocessed option '%s'", argv[i]);
            return     (1);
          }
          break;
        case 'V' :
          fprintf (stderr, "gcv, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 Universite de Bordeaux, INRIA & CNRS\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          return     (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit (&grafdat);
  SCOTCH_geomInit  (&geomdat);
  C_inpFormatTab[C_inpFormatType].func (&grafdat, &geomdat, C_filepntrsrcinp, NULL, C_inpFormatData);
#ifdef SCOTCH_DEBUG_ALL
  if (SCOTCH_graphCheck (&grafdat) != 0) {
    errorPrint ("main: bad graph structure");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ALL */
  C_outFormatTab[C_outFormatType].func (&grafdat, &geomdat, C_filepntrsrcout, C_filepntrgeoout, C_outFormatData);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_geomExit  (&geomdat);
  SCOTCH_graphExit (&grafdat);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
