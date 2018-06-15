/* Copyright 2004,2007-2010,2012,2014 Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dummysizes.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of the libScotch compilation job.  **/
/**                This small program processes files that **/
/**                are in fact pattern header files for    **/
/**                the libScotch library, and replaces     **/
/**                symbolic sizes of the opaque libScotch  **/
/**                by the proper integer values according  **/
/**                to the machine on which it is run.      **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 22 oct 2001     **/
/**                                 to   : 22 nov 2001     **/
/**                # Version 4.0  : from : 25 nov 2001     **/
/**                                 to   : 06 jan 2006     **/
/**                # Version 5.0  : from : 26 apr 2006     **/
/**                                 to   : 03 apr 2008     **/
/**                # Version 5.1  : from : 16 jun 2008     **/
/**                                 to   : 15 aug 2010     **/
/**                # Version 6.0  : from : 01 dec 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DUMMYSIZES

#define CHARMAX                     2048          /* Maximum line size */

#define SUBSMAX                     48            /* Maximum number of substitutions */

#define C_FILENBR                   2             /* Number of files in list                */
#define C_FILEARGNBR                2             /* Number of files which can be arguments */

#define C_filenamehedinp            C_fileTab[0].nameptr /* Source graph input file name */
#define C_filenamehedout            C_fileTab[1].nameptr /* Statistics output file name  */

#define C_filepntrhedinp            C_fileTab[0].fileptr /* Source graph input file */
#define C_filepntrhedout            C_fileTab[1].fileptr /* Statistics output file  */

#define EXPAND(s)                   EXPANDTWO(s)
#define EXPANDTWO(s)                #s

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "geom.h"
#include "mesh.h"
#include "arch.h"
#include "mapping.h"
#include "order.h"
#ifdef SCOTCH_PTSCOTCH
#include "dgraph.h"
#include "dgraph_halo.h"
#include "dmapping.h"
#include "dorder.h"
#include "library_dmapping.h"
#endif /* SCOTCH_PTSCOTCH */
#include "library_mapping.h"
#include "library_order.h"

/*
**  The static definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { "r" },
                              { "w" } };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

void
subsFill (
char *                      substab[2],
char *                      origptr,
int                         subsval)
{
  char *              subsptr;

  subsptr = malloc (32 * sizeof (char));
  sprintf (subsptr, "%d", (int) ((subsval + sizeof (double) - 1) / sizeof (double)));

  substab[0] = origptr;
  substab[1] = subsptr;
}

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
  char                chartab[CHARMAX];
  char                chartmp[CHARMAX];
  char *              substab[SUBSMAX][2];        /* Substitution array */
  int                 subsnbr;
  int                 i;

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    printf ("Usage is:\ndummysizes [<input pattern header file> [<output header file>]]\n");
    return (((argv[1][0] == '?') && argv[1][1] == '\0') ? 0 : 1);
  }

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].fileptr = (C_fileTab[i].modeptr[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        C_fileTab[C_fileNum ++].nameptr = argv[i];
      else {
        fprintf (stderr, "dummysizes: ERROR: main: too many file names given");
        exit    (1);
      }
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          printf ("Usage is:\ndummysizes [<input pattern header file> [<output header file>]]\n");
          exit       (0);
        case 'V' :
          fprintf (stderr, "dummysizes, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007-2010 ENSEIRB, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          fprintf (stderr, "dummysizes: ERROR: main: unprocessed option (\"%s\")", argv[i]);
          exit    (1);
      }
    }
  }

  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].nameptr[0] != '-') ||       /* If not standard stream */
        (C_fileTab[i].nameptr[1] != '\0')) {
      if ((C_fileTab[i].fileptr = fopen (C_fileTab[i].nameptr, C_fileTab[i].modeptr)) == NULL) { /* Open the file */
          fprintf (stderr, "dummysizes: ERROR: main: cannot open file (%d)", i);
          exit    (1);
      }
    }
  }

#ifdef SCOTCH_PTSCOTCH
  substab[0][0] = "library_pt.h";
  substab[0][1] = "ptscotch.h  ";
  substab[1][0] = "library_pt_f.h";
  substab[1][1] = "ptscotchf.h   ";
#else /* SCOTCH_PTSCOTCH */
  substab[0][0] = "library.h";
  substab[0][1] = "scotch.h ";
  substab[1][0] = "library_f.h";
  substab[1][1] = "scotchf.h  ";
#endif /* SCOTCH_PTSCOTCH */
  substab[2][0] = "DUMMYIDX";
  substab[2][1] = EXPAND (IDX);
  substab[3][0] = "DUMMYINT";
  substab[3][1] = EXPAND (INT);
  substab[4][0] = "DUMMYMAXINT";
  substab[4][1] = EXPAND (INTVALMAX);
  substab[5][0] = "DUMMYNUMSTRING";
  substab[5][1] = "\"" GNUMSTRING "\"";
  substab[6][0] = "DUMMYVERSION";
  substab[6][1] = EXPAND (SCOTCH_VERSION);
  substab[7][0] = "DUMMYRELEASE";
  substab[7][1] = EXPAND (SCOTCH_RELEASE);
  substab[8][0] = "DUMMYPATCHLEVEL";
  substab[8][1] = EXPAND (SCOTCH_PATCHLEVEL);
  subsnbr = 9;
  subsFill (substab[subsnbr ++], "DUMMYSIZEARCH",          sizeof (Arch));
  subsFill (substab[subsnbr ++], "DUMMYSIZEGEOM",          sizeof (Geom));
  subsFill (substab[subsnbr ++], "DUMMYSIZEGRAPH",         sizeof (Graph));
  subsFill (substab[subsnbr ++], "DUMMYSIZEMESH",          sizeof (Mesh));
  subsFill (substab[subsnbr ++], "DUMMYSIZEMAP",           sizeof (LibMapping));
  subsFill (substab[subsnbr ++], "DUMMYSIZEORDER",         sizeof (LibOrder));
  subsFill (substab[subsnbr ++], "DUMMYSIZESTRAT",         sizeof (Strat *));
#ifdef SCOTCH_PTSCOTCH
  subsFill (substab[subsnbr ++], "DUMMYSIZEDGRAPHHALOREQ", sizeof (DgraphHaloRequest)); /* TRICK: before DUMMYSIZEDGRAPH */
  subsFill (substab[subsnbr ++], "DUMMYSIZEDGRAPH",        sizeof (Dgraph));
  subsFill (substab[subsnbr ++], "DUMMYSIZEDMAP",          sizeof (LibDmapping));
  subsFill (substab[subsnbr ++], "DUMMYSIZEDORDER",        sizeof (Dorder));
#else /* SCOTCH_PTSCOTCH */
  subsFill (substab[subsnbr ++], "DUMMYSIZEDGRAPHHALOREQ", 1); /* TRICK: before DUMMYSIZEDGRAPH */
  subsFill (substab[subsnbr ++], "DUMMYSIZEDGRAPH",        1);
  subsFill (substab[subsnbr ++], "DUMMYSIZEDMAP",          1);
  subsFill (substab[subsnbr ++], "DUMMYSIZEDORDER",        1);
#endif /* SCOTCH_PTSCOTCH */

  while (fgets (chartab, CHARMAX, C_filepntrhedinp) != NULL) { /* Infinite loop on file lines */
    int                 charnbr;
    int                 subsnum;

    if (((charnbr = strlen (chartab)) >= (CHARMAX - 1)) && /* If line read is at least as long as maximum size     */
        (chartab[CHARMAX - 1] != '\n')) {         /* And last character is not a newline, that is, some is missing */
      fprintf (stderr, "dummysizes: ERROR: line too long\n");
      exit    (1);
    }

    for (subsnum = 0; subsnum < subsnbr; subsnum ++) { /* Perform substitutions */
      char *              charptr;                /* Place where token found    */

      while ((charptr = strstr (chartab, substab[subsnum][0])) != NULL) { /* As long as substitution can be performed */
        int                 charnbr;
        int                 charnum;

        charnum = charptr - chartab;              /* Position where token found */
        charnbr = strlen (substab[subsnum][0]);   /* Length of token            */

        strcpy (chartmp, charptr + charnbr);      /* Save end of line */

        sprintf (charptr, "%s%s", substab[subsnum][1], chartmp); /* Replace end of line with substituted token */
      }
    }

    fputs (chartab, C_filepntrhedout);            /* Output possibly updated line */
  }

#ifdef SCOTCH_DEBUG_MAIN1
  for (i = 0; i < C_FILENBR; i ++) {              /* For all file names     */
    if ((C_fileTab[i].name[0] != '-') ||          /* If not standard stream */
        (C_fileTab[i].name[1] != '\0')) {
      fclose (C_fileTab[i].pntr);                 /* Close the stream */
    }
  }
#endif /* SCOTCH_DEBUG_MAIN1 */

  exit (0);
}
