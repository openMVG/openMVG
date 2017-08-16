/* Copyright 2009-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gscat.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program writes a centralized       **/
/**                source graph file in the form of a      **/
/**                distributed source graph.               **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 26 apr 2009     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GSCAT

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gscat.h"

/*
**  The static variables.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {     /* Usage */
  "gscat <nparts> <input source file> <output target file pattern> <options>",
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
  SCOTCH_Num          p[1] = { 1 };               /* Number of parts */
  int                 i;

  errorProg ("gscat");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 1) {                        /* If number of parameters not reached              */
        if ((p[C_paraNum ++] = atoi (argv[i])) < 1) { /* Get number of parts                          */
          errorPrint ("main: invalid number of parts '%s'", argv[i]);
          return     (1);
        }
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
          fprintf (stderr, "gscat, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2009-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
          return     (1);
      }
    }
  }

  fileBlockOpen (C_fileTab, 1);                   /* Open input graph file */

  C_graphScat (C_filepntrsrcinp, p[0], C_filenamesrcout);

  fileBlockClose (C_fileTab, 1);                  /* Always close explicitely to end eventual (un)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}

static
int
C_graphScat (
FILE * const                stream,
SCOTCH_Num                  procnbr,
char * const                nameptr)
{
  SCOTCH_Num          versval;
  SCOTCH_Num          propval;
  char                proptab[4];
  int                 flagtab[3];
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertglbnbr;
  SCOTCH_Num          edgeglbnbr;
  SCOTCH_Num          procnum;

  if (intLoad (stream, &versval) != 1) {          /* Read version number */
    errorPrint ("C_graphScat: bad input (1)");
    return     (1);
  }
  if (versval != 0) {                             /* If version not zero */
    errorPrint ("C_graphScat: only centralized graphs supported");
    return     (1);
  }

  if ((intLoad (stream, &vertglbnbr) != 1) ||     /* Read rest of header */
      (intLoad (stream, &edgeglbnbr) != 1) ||
      (intLoad (stream, &baseval)    != 1) ||
      (intLoad (stream, &propval)    != 1) ||
      (propval < 0)                        ||
      (propval > 111)) {
    errorPrint ("C_graphScat: bad input (2)");
    return     (1);
  }
  sprintf (proptab, "%3.3d", (int) propval);      /* Compute file properties */
  flagtab[0] = proptab[0] - '0';                  /* Vertex labels flag      */
  flagtab[1] = proptab[1] - '0';                  /* Edge weights flag       */
  flagtab[2] = proptab[2] - '0';                  /* Vertex loads flag       */

  for (procnum = 0; procnum < procnbr; procnum ++) {
    char *              nametmp;
    FILE *              ostream;
    SCOTCH_Num          vertlocnbr;
    SCOTCH_Num          vertlocnum;
    SCOTCH_Num          edgelocnbr;

    nametmp = nameptr;
    ostream = NULL;
    if (fileNameDistExpand (&nametmp, procnbr, procnum, -1) == 0) {
      ostream = fopen (nametmp, "w+");
      memFree (nametmp);                            /* Expanded name no longer needed anyway */
    }
    if (ostream == NULL) {
      errorPrint ("C_graphScat: cannot open file");
      return     (1);
    }

    vertlocnbr = DATASIZE (vertglbnbr, procnbr, procnum);

    if (fprintf (ostream, "2\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t%015d\n" SCOTCH_NUMSTRING "\t%3s\n", /* Write file header */
                 (SCOTCH_Num) procnbr,
                 (SCOTCH_Num) procnum,
                 (SCOTCH_Num) vertglbnbr,
                 (SCOTCH_Num) edgeglbnbr,
                 (SCOTCH_Num) vertlocnbr,
                 0,                               /* Number of edges not yet known */
                 (SCOTCH_Num) baseval,
                 proptab) == EOF) {
      errorPrint ("C_graphScat: bad output (1)");
      return     (1);
    }

    for (vertlocnum = edgelocnbr = 0; vertlocnum < vertlocnbr; vertlocnum ++) {
      SCOTCH_Num          degrval;

      if (flagtab[0] != 0) {                      /* If must read label               */
        SCOTCH_Num          vlblval;              /* Value where to read vertex label */

        if (intLoad (stream, &vlblval) != 1) {    /* Read label data */
          errorPrint ("C_graphScat: bad input (3)");
          return     (1);
        }
        intSave (ostream, vlblval);
        putc ('\t', ostream);
      }
      if (flagtab[2] != 0) {                      /* If must read vertex load        */
        SCOTCH_Num          veloval;              /* Value where to read vertex load */

        if (intLoad (stream, &veloval) != 1) {    /* Read vertex load data    */
          errorPrint ("C_graphScat: bad input (4)");
          return     (1);
        }
        intSave (ostream, veloval);
        putc ('\t', ostream);
      }
      if (intLoad (stream, &degrval) != 1) {      /* Read vertex degree */
        errorPrint ("C_graphScat: bad input (5)");
        return     (1);
      }
      intSave (ostream, degrval);

      edgelocnbr += degrval;

      for ( ; degrval > 0; degrval --) {
        SCOTCH_Num          edgeval;              /* Value where to read edge end */

        if (flagtab[1] != 0) {                    /* If must read edge load        */
          SCOTCH_Num          edloval;            /* Value where to read edge load */

          if (intLoad (stream, &edloval) != 1) {  /* Read edge load data    */
            errorPrint ("C_graphScat: bad input (6)");
            return     (1);
          }
          putc ('\t', ostream);
          intSave (ostream, edloval);
        }

        if (intLoad (stream, &edgeval) != 1) {    /* Read edge data */
          errorPrint ("C_graphScat: bad input (7)");
          return     (1);
        }
        putc ('\t', ostream);
        intSave (ostream, edgeval);
      }
      putc ('\n', ostream);
    }

    rewind (ostream);

    if (fprintf (ostream, "2\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t%015lld\n" SCOTCH_NUMSTRING "\t%3s\n", /* Write file header */
                 (SCOTCH_Num) procnbr,
                 (SCOTCH_Num) procnum,
                 (SCOTCH_Num) vertglbnbr,
                 (SCOTCH_Num) edgeglbnbr,
                 (SCOTCH_Num) vertlocnbr,
                 (long long)  edgelocnbr,         /* Now we know the exact number of edges */
                 (SCOTCH_Num) baseval,
                 proptab) == EOF) {
      errorPrint ("C_graphScat: bad output (2)");
      return     (1);
    }

    fclose (ostream);
  }

  return (0);
}
