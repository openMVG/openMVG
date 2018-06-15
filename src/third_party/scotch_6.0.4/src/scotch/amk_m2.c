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
/**   NAME       : amk_m2.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the distance map for            **/
/**                bidimensional mesh graphs, to be used   **/
/**                to build the architecture description   **/
/**                files for these graphs.                 **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 21 apr 1994     **/
/**                                 to   : 21 apr 1994     **/
/**                # Version 2.0  : from : 12 jul 1994     **/
/**                                 to   : 12 nov 1994     **/
/**                # Version 3.0  : from : 17 jul 1995     **/
/**                                 to   : 19 sep 1995     **/
/**                # Version 3.2  : from : 31 may 1997     **/
/**                                 to   : 02 jun 1997     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 09 feb 2004     **/
/**                                 to   : 09 feb 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 21 apr 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/**   NOTES      : # The vertices of the (dX,dY) mesh are  **/
/**                  numbered as terminals so that         **/
/**                  t (0,0)=0, t (1,0) = 1,               **/
/**                  t (dX - 1, 0) = dX - 1,               **/
/**                  t (0,1) = dX, and                     **/
/**                  t(x,y) = (y * dX) + x.                **/
/**                                                        **/
/**                # The nested dissection method should   **/
/**                  behave like the architecture built in **/
/**                  the mapper.                           **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define AMK_M2

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "arch.h"
#include "arch_mesh.h"
#include "amk_m2.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { "w" } };

static const char *         C_usageList[] = {
  "amk_m2 <dimX> [<dimY> [<output target file>]] <options>",
  "  -h          : Display this help",
  "  -m<method>  : Decomposition method",
  "                  n  : Nested dissection (cut biggest dimension)",
  "                  o  : One-way dissection (y, then x)",
  "  -V          : Print program version and copyright",
  "",
  "Default option set is : '-Mn'",
  NULL };

/*************************************************/
/*                                               */
/* The main routine, which computes the distance */
/* triangular table.                             */
/*                                               */
/*************************************************/

int
main (
int                         argc,
char *                      argv[])
{
  ArchMesh2        arch;                          /* Mesh dimensions            */
  ArchMesh2Dom     dom;                           /* Initial domain             */
  C_MethType       methtype;                      /* Bipartitioning method      */
  unsigned int     termnbr;                       /* Number of terminal domains */
  unsigned int     termmax;                       /* Maximum terminal number    */
  unsigned int *   termtab;                       /* Terminal numbers table     */
  unsigned int     x0, y0, x1, y1;
  int              i;

  errorProg ("amk_m2");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  methtype  = C_METHNESTED;
  arch.c[0] =                                     /* Preset mesh dimensions */
  arch.c[1] = 1;

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 2) {                        /* If number of parameters not reached              */
        if ((arch.c[C_paraNum ++] = atoi (argv[i])) < 1) { /* Get the dimension                       */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
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
        case 'M' :                                /* Use a built-in method */
        case 'm' :
          switch (argv[i][2]) {
            case 'N' :                            /* Nested dissection */
            case 'n' :
              methtype = C_METHNESTED;
              break;
            case 'O' :                            /* One-way dissection */
            case 'o' :
              methtype = C_METHONEWAY;
              break;
            default :
              errorPrint ("main: unprocessed option '%s'", argv[i]);
              return     (1);
          }
          break;
        case 'H' :                               /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'V' :
          fprintf (stderr, "amk_m2, version " SCOTCH_VERSION_STRING "\n");
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

  dom.c[0][0] = 0;                                /* Set the initial domain */
  dom.c[0][1] = arch.c[0] - 1;
  dom.c[1][0] = 0;
  dom.c[1][1] = arch.c[1] - 1;

  termnbr = arch.c[0] * arch.c[1];                /* Compute number of terminals                                    */
  termmax = 0;                                    /* Maximum terminal value not known yet                           */
  if ((termtab = (unsigned int *) memAlloc (termnbr * sizeof (unsigned int))) == NULL) { /* Allocate terminal array */
    errorPrint ("main: out of memory");
    return     (1);
  }
  memset (termtab, -1, termnbr * sizeof (unsigned int)); /* Initilize mapping table */

  C_termBipart (&arch, &dom, 1, termtab, &termmax, /* Compute terminal numbers */
                (methtype == C_METHNESTED) ? archMesh2DomBipart : C_methBipartOne);

  fprintf (C_filepntrarcout, "deco\n0\n%u\t%u\n", /* Print file header                */
           termnbr,                               /* Print number of terminal domains */
           termmax);                              /* Print biggest terminal value     */
  for (i = 0; i < termnbr; i ++)                  /* For all terminals                */
    fprintf (C_filepntrarcout, "%u\t1\t%u\n",     /* Print terminal data              */
             i, termtab[i]);

  for (y0 = 0; y0 < arch.c[1]; y0 ++) {           /* For all vertices */
    for (x0 = 0; x0 < arch.c[0]; x0 ++) {
      for (y1 = 0; y1 <= y0; y1 ++) {             /* Compute distance to smaller vertices */
        for (x1 = 0; (x1 < arch.c[0]) && ((y1 < y0) || (x1 < x0)); x1 ++)
          fprintf (C_filepntrarcout,
                   ((x1 == 0) && (y1 == 0)) ? "%u" : " %u",
                   C_termDist (x0, y0, x1, y1));
      }
      fprintf (C_filepntrarcout, "\n");
    }
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  memFree (termtab);                              /* Free terminal number array */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}

/* This routine recursively determines the values
** of all the terminal vertices of the mesh domain,
** and puts them in table.
*/

void
C_termBipart (
ArchMesh2 *                 archptr,
ArchMesh2Dom *              domptr,
unsigned int                num,
unsigned int *              termtab,
unsigned int *              termmax,
int                     (*  methfunc) ())
{
  ArchMesh2Dom        dom0;
  ArchMesh2Dom        dom1;

  if (methfunc (archptr, domptr, &dom0, &dom1) == 0) { /* If we can bipartition                          */
    C_termBipart (archptr, &dom0, num + num,     termtab, termmax, methfunc); /* Bipartition recursively */
    C_termBipart (archptr, &dom1, num + num + 1, termtab, termmax, methfunc);
  }
  else {                                          /* If we have reached the end */
    termtab[domptr->c[1][0] * archptr->c[0] +     /* Set the terminal number    */
            domptr->c[0][0]] = num;
    if (*termmax < num)                           /* If we have reached a new maximum */
      *termmax = num;                             /* Record it                        */
  }
}

/* This function tries to split a rectangular domain into
** two subdomains using a one-way dissection method, by
** cutting first in the Y dimension, then in the X one.
** It returns:
** - 1  : if the domain has been bipartitioned.
** - 0  : else.
*/

int
C_methBipartOne (
const ArchMesh2 * const       archptr,
const ArchMesh2Dom * const    domptr,
ArchMesh2Dom * restrict const dom0ptr,
ArchMesh2Dom * restrict const dom1ptr)
{
  if ((domptr->c[0][0] == domptr->c[0][1]) &&     /* Return if cannot split more */
      (domptr->c[1][0] == domptr->c[1][1]))
    return (0);

  if (domptr->c[1][1] == domptr->c[1][0]) {       /* If the Y dimension cannot be cut */
    dom0ptr->c[0][0] = domptr->c[0][0];           /* Cut in the X dimension           */
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];
  }
  else {                                          /* If the Y dimension can be cut, cut it */
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];
  }

  return (1);
}
