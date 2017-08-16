/* Copyright 2004,2007,2010,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_vw.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module outputs the state of the    **/
/**                current partition on the form of a      **/
/**                Scotch mapping file.                    **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 18 may 2004     **/
/**                                 to     18 may 2004     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 10 oct 2013     **/
/**                                 to     10 oct 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_VW

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_vw.h"

/*
**  The static variables.
*/

static int                  vgraphseparatevwfilenum = 0; /* Number of file to output */

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine outputs the mapping file.
** It returns:
** - 0   : if the file could be produced.
** - !0  : on error.
*/

int
vgraphSeparateVw (
Vgraph * restrict const             grafptr)      /*+ Separation graph +*/
{
  char                nametab[64];                /* File name */
  FILE * restrict     fileptr;
  Gnum                vertnum;                    /* Vertex number */

  sprintf (nametab, "vgraphseparatevw_output_%08d.map", vgraphseparatevwfilenum ++);
  if ((fileptr = fopen (nametab, "w+")) == NULL) {
    errorPrint ("vgraphSeparateVw: cannot open partition file");
    return     (1);
  }

  fprintf (fileptr, GNUMSTRING "\n",              /* Output size of mapping; test if failure later, in main loop */
           (Gnum) grafptr->s.vertnbr);

  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
    if (fprintf (fileptr, GNUMSTRING "\t%d\n",
                 (Gnum) ((grafptr->s.vnumtax != NULL) ? grafptr->s.vnumtax[vertnum] : vertnum),
                 (int) grafptr->parttax[vertnum]) <= 0) {
      errorPrint ("vgraphSeparateVw: bad output");
      fclose     (fileptr);
      return     (1);
    }
  }

  fclose (fileptr);
  return (0);
}
