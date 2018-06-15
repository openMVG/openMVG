/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : mapping_io.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles (partial) mappings. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     30 mar 1999     **/
/**                # Version 3.4  : from : 11 sep 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to     14 nov 2005     **/
/**                # Version 5.0  : from : 13 sep 2006     **/
/**                                 to     27 feb 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     22 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MAPPING_IO

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "mapping_io.h"

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine reads the contents of the given mapping
** file to the given mapping, reordering vertices
** according to the given vertex label table if necessary.
** It returns:
** - 0   : if mapping successfully written.
** - 1   : on error.
** - 2   : variable-sized architectures cannot be loaded.
*/

/* TODO remove it */
int /* TODO copy the one from gout? */
mapLoad (
Mapping * restrict const        mappptr,
const Gnum * restrict const     vlbltab,
FILE * restrict const           stream)
{
  Gnum                  vertnum;
  Gnum                  mappnum;
  Gnum                  mappnbr;
  MappingLoadMap *      mapptab;                  /* Mapping array                     */
  MappingLoadPerm *     permtab;                  /* Array of sorted label/index pairs */
  Anum                  archnbr;                  /* Size of the target architecture   */
  ArchDom               fdomdat;                  /* First domain of architecture      */

  if (strcmp (archName (mappptr->archptr), "term") == 0) /* If target architecture is variable-sized */
    return (2);

  archDomFrst (mappptr->archptr, &fdomdat);      /* Get architecture size */
  archnbr = archDomSize (mappptr->archptr, &fdomdat);
  if (mappptr->domnmax < (archnbr + 1)) {         /* If mapping array too small to store mapping data */
    ArchDom * restrict    domntab;

    if ((domntab = (ArchDom *) memRealloc (mappptr->domntab, (archnbr + 1) * sizeof (ArchDom))) == NULL) { /* If cannot resize domain array */
      errorPrint ("mapLoad: out of memory (1)");
      return     (1);
    }

    mappptr->domnmax = archnbr + 1;               /* Point to new array */
    mappptr->domntab = domntab;
  }
  mappptr->domnnbr = archnbr + 1;                 /* One more for first domain, for unmapped vertices                 */
  archDomFrst (mappptr->archptr, &mappptr->domntab[0]); /* Set first domain with root domain data                    */
  for (mappnum = 0; mappnum < archnbr; mappnum ++) /* For all terminal domain numbers                                 */
    archDomTerm (mappptr->archptr, &mappptr->domntab[mappnum + 1], mappnum); /* Set domain with terminal domain data */

  if ((intLoad (stream, &mappnbr) != 1) ||        /* Read number of mapping entries */
      (mappnbr < 1)) {
    errorPrint ("mapLoad: bad input (1)");
    return     (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &mapptab, (size_t) (mappnbr          * sizeof (MappingLoadMap)),
                     &permtab, (size_t) (mappptr->grafptr->vertnbr * sizeof (MappingLoadPerm)), NULL) == NULL) {
    errorPrint ("mapLoad: out of memory (2)");
    return     (1);
  }

  for (mappnum = 0; mappnum < mappnbr; mappnum ++) { /* Load mapping array */
    if ((intLoad (stream, &mapptab[mappnum].slblnum) != 1) ||
        (intLoad (stream, &mapptab[mappnum].tlblnum) != 1)) {
      errorPrint ("mapLoad: bad input (2)");
      return     (1);
    }
  }
  intSort2asc1 (mapptab, mappnbr);                /* Sort mapping array by increasing source labels */

  if (vlbltab != NULL) {                          /* If graph has vertex labels */
    Gnum                vertnum;

    for (vertnum = 0; vertnum < mappptr->grafptr->vertnbr; vertnum ++) { /* Build inverse permutation */
      permtab[vertnum].vertnum = vertnum + mappptr->grafptr->baseval;
      permtab[vertnum].vlblnum = vlbltab[vertnum];
    }
    intSort2asc1 (permtab, mappptr->grafptr->vertnbr); /* Sort vertex array by increasing labels */
  }
  else {
    Gnum                vertnum;

    for (vertnum = 0; vertnum < mappptr->grafptr->vertnbr; vertnum ++) { /* Build identity permutation */
      permtab[vertnum].vertnum = vertnum + mappptr->grafptr->baseval;
      permtab[vertnum].vlblnum = vertnum + mappptr->grafptr->baseval;
    }
  }

  for (vertnum = 0, mappnum = 0;                  /* For all graph vertices */
       vertnum < mappptr->grafptr->vertnbr; vertnum ++) {
    while ((mappnum < mappnbr) &&                 /* Skip useless mapping data (if graph is subgraph of originally mapped graph) */
           (permtab[vertnum].vlblnum > mapptab[mappnum].slblnum))
      mappnum ++;
    if (mappnum >= mappnbr)                       /* If all mapping data exhausted */
      break;                                      /* Exit the matching loop        */

    if (permtab[vertnum].vlblnum == mapptab[mappnum].slblnum) { /* If matching mapping data found */
      if ((mapptab[mappnum].tlblnum >= 0) &&      /* If mapping valid                             */
          (mapptab[mappnum].tlblnum < archnbr))
        mappptr->parttax[permtab[vertnum].vertnum] = mapptab[mappnum].tlblnum + 1; /* Set mapping to terminal domain */
      mappnum ++;                                 /* Mapping pair has been used                               */
    }
  }

  memFree (mapptab);                              /* Free group leader */

  return (0);
}

/* This routine writes the contents of the
** given mapping to the given string.
** It returns:
** - 0   : if mapping successfully written.
** - !0  : on error.
*/

int
mapSave (
const Mapping * restrict const  mappptr,
FILE * restrict const           stream)
{

  Gnum                  vertnnd;
  Gnum                  vertnum;

  const Arch * restrict const     archptr = mappptr->archptr;
  const ArchDom * restrict const  domntab = mappptr->domntab;
  const Anum * restrict const     parttax = mappptr->parttax;
  const Gnum * restrict const     vlbltax = mappptr->grafptr->vlbltax;

  vertnum = mappptr->grafptr->baseval;
  vertnnd = mappptr->grafptr->vertnbr;            /* Un-based number at first */

  if (fprintf (stream, GNUMSTRING "\n",
               (Gnum) vertnnd) == EOF) {
    errorPrint ("mapSave: bad output (1)");
    return     (1);
  }

  for (vertnnd += vertnum; vertnum < vertnnd; vertnum ++) {
    if (fprintf (stream, GNUMSTRING "\t" ANUMSTRING "\n",
                 (Gnum) ((vlbltax != NULL) ? vlbltax[vertnum] : vertnum),
                 (Anum) (parttax != NULL) ? archDomNum (archptr, &domntab[parttax[vertnum]]) : -1) == EOF) {
      errorPrint ("mapSave: bad output (2)");
      return     (1);
    }
  }

  return (0);
}
