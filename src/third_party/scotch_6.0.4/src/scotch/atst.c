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
/**   NAME       : atst.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Target architecture graph analyzer.     **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 17 may 1994     **/
/**                                 to   : 17 may 1994     **/
/**                # Version 2.0  : from : 11 nov 1994     **/
/**                                 to   : 11 nov 1994     **/
/**                # Version 3.0  : from : 05 jul 1995     **/
/**                                 to   : 19 aug 1995     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to   : 12 may 1998     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 09 feb 2004     **/
/**                                 to   : 23 nov 2005     **/
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

#define ATST

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "arch.h"
#include "arch_deco.h"
#include "arch_mesh.h"
#include "atst.h"

/*
**  The static variables.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "r" },
                              { "w" } };

static const char *         C_usageList[] = {
  "atst [<input target file> [<output data file>]] <options>",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (argc, argv)
int                 argc;
char *              argv[];
{
  Arch                archdat;                    /* The architecture read */
  ArchDeco *          deco;                       /* Its decomposition     */
  ArchDecoDom         dom0;
  ArchDecoDom         dom1;
  Anum                dstval;
  Anum                dstmin;
  Anum                dstmax;
  Anum                dstsum;
  double              dstavg;
  double              dstdlt;
  int                 i;

  errorProg ("atst");

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
        case 'V' :
          fprintf (stderr, "atst, version " SCOTCH_VERSION_STRING "\n");
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

  archInit (&archdat);                            /* Initialize architecture structure */
  archLoad (&archdat, C_filepntrtgtinp);          /* Load architecture                 */
  if (strcmp (archName (&archdat), "deco") != 0) { /* If it is not a decomposition     */
    errorPrint ("main: architecture is not decomposition-defined");
    return     (1);
  }
  deco = (ArchDeco *) (void *) &archdat.data;     /* Point to the decomposition */

  dstmin = (Anum) (((unsigned long) ((Anum) -1)) >> 1); /* Set to maximum number in Anum */
  dstmax = 0;
  dstsum = 0;

  for (dom0.num = 1; dom0.num <= deco->domvertnbr; dom0.num ++) { /* For all pairs of vertices */
    if (archDecoDomSize (deco, &dom0) == 1) {     /* If vertex is a terminal                   */
      for (dom1.num = dom0.num + 1; dom1.num <= deco->domvertnbr; dom1.num ++) {
        if (archDecoDomSize (deco, &dom1) == 1) { /* If vertex is a terminal               */
          dstval = archDecoDomDist (deco, &dom0, &dom1); /* Compute distance between pairs */
          if (dstmin > dstval)
            dstmin = dstval;
          if (dstmax < dstval)
            dstmax = dstval;
          dstsum += dstval;                       /* Compute distance between pairs */
        }
      }
    }
  }
  dstavg = (deco->domtermnbr > 1)
           ? (double) dstsum / (double) (deco->domtermnbr * (deco->domtermnbr - 1) / 2)
           : 0.0L;
  dstdlt = 0.0L;
  for (dom0.num = 1; dom0.num <= deco->domvertnbr; dom0.num ++) { /* For all pairs of vertices */
    if (archDecoDomSize (deco, &dom0) == 1) {     /* If vertex is a terminal                   */
      for (dom1.num = dom0.num + 1; dom1.num <= deco->domvertnbr; dom1.num ++) {
        if (archDecoDomSize (deco, &dom1) == 1)   /* If vertex is a terminal */
          dstdlt += fabs (archDecoDomDist (deco, &dom0, &dom1) - dstavg);
      }
    }
  }
  if (deco->domtermnbr > 1)
    dstdlt /= (double) (deco->domtermnbr * (deco->domtermnbr - 1) / 2);

  fprintf (C_filepntrlogout, "A\tTerminals\tnbr=" SCOTCH_NUMSTRING "\n",
           (SCOTCH_Num) deco->domtermnbr);
  fprintf (C_filepntrlogout, "A\tDistance\tmin=" SCOTCH_NUMSTRING "\tmax=" SCOTCH_NUMSTRING "\tavg=%g\tdlt=%g\n",
           (SCOTCH_Num) dstmin, (SCOTCH_Num) dstmax, dstavg, dstdlt);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  archExit (&archdat);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}
