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
/**   NAME       : amk_grf.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the architecture description    **/
/**                file for any source graph.              **/
/**                                                        **/
/**   DATES      : # Version 3.0  : from : 06 jul 1995     **/
/**                                 to   : 02 oct 1995     **/
/**                # Version 3.1  : from : 26 mar 1996     **/
/**                                 to   : 26 mar 1996     **/
/**                # Version 3.2  : from : 23 apr 1997     **/
/**                                 to   : 03 jun 1998     **/
/**                # Version 3.3  : from : 15 may 1999     **/
/**                                 to   : 15 may 1999     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to   : 17 mar 2005     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 11 dec 2008     **/
/**                                 to   : 17 jul 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define AMK_GRF

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "amk_grf.h"

/*
**  The static variables.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" },
                              { "r" } };

static const char *         C_usageList[] = {     /* Usage */
  "amk_grf [<input source file> [<output target file>]] <options>",
  "  -b<strat>  : Apply bipartitioning strategy <strat>",
  "  -h         : Display this help",
  "  -l<file>   : Load vertex list from <file>",
  "  -V         : Print program version and copyright",
  "",
  "Default option set is : '-Bhf{move=1000}/((load0=load)|(load0=0))?x;'",
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
  SCOTCH_Strat        bipastrat;                  /* Bipartitioning strategy                   */
  SCOTCH_Arch         archdat;                    /* Target (terminal) architecture            */
  SCOTCH_Graph        grafdat;                    /* Source graph to turn into architecture    */
  SCOTCH_Num          vertnbr;                    /* Number of vertices in graph               */
  SCOTCH_Num *        vlbltab;                    /* Pointer to vertex label array, if present */
  SCOTCH_Num          listnbr;                    /* Size of list array                        */
  SCOTCH_Num *        listtab;                    /* Pointer to list array                     */
  C_VertSort *        sorttab;                    /* Vertex label sort area                    */
  SCOTCH_Num          baseval;
  int                 flagval;                    /* Process flags                             */
  int                 i;

  errorProg ("amk_grf");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  flagval = C_FLAGNONE;
  SCOTCH_stratInit (&bipastrat);

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
        case 'B' :                                /* Bipartitioning strategy */
        case 'b' :
          SCOTCH_stratExit (&bipastrat);
          SCOTCH_stratInit (&bipastrat);
          if ((SCOTCH_stratGraphBipart (&bipastrat, &argv[i][2])) != 0) {
            errorPrint ("main: invalid bipartitioning strategy");
            return     (1);
          }
          break;
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'L' :                                /* Input vertex list */
        case 'l' :
          flagval |= C_FLAGVRTINP;
          if (argv[i][2] != '\0')
            C_filenamevrtinp = &argv[i][2];
          break;
        case 'V' :
          fprintf (stderr, "amk_grf, version " SCOTCH_VERSION_STRING "\n");
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

  SCOTCH_graphInit (&grafdat);                    /* Create graph structure           */
  SCOTCH_graphLoad (&grafdat, C_filepntrgrfinp, -1, 0); /* Load source graph          */
  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, NULL, NULL, NULL, /* Get graph data */
                    &vlbltab, NULL, NULL, NULL);

  listnbr = 0;                                    /* Initialize vertex list */
  listtab = NULL;
  if (flagval & C_FLAGVRTINP) {                   /* If list of vertices provided */
    SCOTCH_Num          listnum;

    if ((intLoad (C_filepntrvrtinp, &listnbr) != 1) || /* Read list size */
        (listnbr < 0)                               ||
        (listnbr > vertnbr)) {
      errorPrint ("main: bad list input (1)");
      return     (1);
    }
    if ((listtab = (SCOTCH_Num *) memAlloc (listnbr * sizeof (SCOTCH_Num) + 1)) == NULL) {
      errorPrint ("main: out of memory (1)");
      return     (1);
    }
    for (listnum = 0; listnum < listnbr; listnum ++) { /* Read list data */
      if (intLoad (C_filepntrvrtinp, &listtab[listnum]) != 1) {
        errorPrint ("main: bad list input (2)");
        return     (1);
      }
    }
    intSort1asc1 (listtab, listnbr);
    for (listnum = 0; listnum < listnbr - 1; listnum ++) { /* Search for duplicates */
      if (listtab[listnum] == listtab[listnum + 1]) {
        errorPrint ("main: duplicate list labels");
        memFree    (listtab);
        return     (1);
      }
    }

    if (vlbltab != NULL) {                        /* If graph has vertex labels */
      SCOTCH_Num          vertnum;

      if ((sorttab = (C_VertSort *) memAlloc (vertnbr * sizeof (C_VertSort))) == NULL) {
        errorPrint ("main: out of memory (2)");
        memFree    (listtab);
        return     (1);
      }
      for (vertnum = 0; vertnum < vertnbr; vertnum ++) { /* Initialize sort area */
        sorttab[vertnum].vlblnum = vlbltab[vertnum];
        sorttab[vertnum].vertnum = vertnum;
      }
      intSort2asc1 (sorttab, vertnbr);            /* Sort by ascending labels */

      for (listnum = 0, vertnum = 0; listnum < listnbr; listnum ++) {  /* For all labels in list */
        while ((vertnum < vertnbr) && (sorttab[vertnum].vlblnum < listtab[listnum]))
          vertnum ++;                             /* Search vertex graph with corresponding label */
        if ((vertnum >= vertnbr) ||               /* If label not found                           */
            (sorttab[vertnum].vlblnum > listtab[listnum])) {
          errorPrint ("main: list label '" SCOTCH_NUMSTRING "' not in graph", (SCOTCH_Num) listtab[listnum]);
          memFree    (sorttab);
          memFree    (listtab);
          return     (1);
        }
        listtab[listnum] = sorttab[vertnum ++].vertnum; /* Replace label by number */
      }
      memFree (sorttab);                          /* Free sort area */
    }
  }

  SCOTCH_archInit  (&archdat);                    /* Initialize target architecture            */
  SCOTCH_archBuild (&archdat, &grafdat, listnbr, listtab, &bipastrat); /* Compute architecture */
  SCOTCH_archSave  (&archdat, C_filepntrtgtout);  /* Write target architecture                 */

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end potential (un)compression tasks */

  SCOTCH_graphExit (&grafdat);                    /* Free target graph        */
  SCOTCH_archExit  (&archdat);                    /* Free target architecture */
  SCOTCH_stratExit (&bipastrat);                  /* Free strategy string     */
  if (listtab != NULL)                            /* If vertex list provided  */
    memFree (listtab);                            /* Free it                  */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
