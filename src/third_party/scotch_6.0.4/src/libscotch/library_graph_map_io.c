/* Copyright 2011,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_map_io.c                  **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the API mappings    **/
/**                of the libSCOTCH library.               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 16 apr 2011     **/
/**                                 to     23 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "library_mapping.h"
#include "library_graph_map_io.h"
#include "scotch.h"

/*************************************/
/*                                   */
/* These routines are the C API for  */
/* the API mapping handling          */
/* routines.                         */
/*                                   */
/*************************************/

/*+ This routine loads the contents of the
*** given mapping array from the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphTabLoad (
const SCOTCH_Graph * const    actgrafptr,         /*+ Graph to map  +*/
SCOTCH_Num * const            parttab,            /*+ Array to load +*/
FILE * const                  stream)             /*+ Input stream  +*/
{
  VertSort *            vertsorttab;              /* Pointer to graph sorting array           */
  int                   vertsortflag;             /* Flag set if graph data sorted by label   */
  VertSort *            mappsorttab;              /* Pointer to mapping data sorting array    */
  int                   mappsortflag;             /* Flag set if mapping data sorted by label */
  Gnum                  mappsortval;
  Gnum                  mappfileval;
  Gnum                  mappfilenbr;              /* Number of mapping pairs in file          */
  Gnum                  mappfilenum;              /* Counter of mapping pairs in file         */
  Gnum *                mappfiletab;              /* Pointer to mapping data read from file   */
  Graph *               grafptr;
  Gnum                  vertnum;
  Gnum                  i, j;

  grafptr = (Graph *) actgrafptr;
  memSet (parttab, ~0, grafptr->vertnbr * sizeof (Anum)); /* Pre-initialize the partition array */

  if ((fscanf (stream, GNUMSTRING,&mappfileval) != 1) || /* Read number of mapping pairs */
      (mappfileval < 1)) {
    errorPrint ("SCOTCH_graphTabLoad: bad input (1)");
    return     (1);
  }
  mappfilenbr = (Gnum) mappfileval;

  if (memAllocGroup ((void **) (void *)           /* Allocate temporary data */
                     &mappfiletab, (size_t) (mappfilenbr      * sizeof (Gnum)),
                     &mappsorttab, (size_t) (mappfilenbr      * sizeof (VertSort)),
                     &vertsorttab, (size_t) (grafptr->vertnbr * sizeof (VertSort)), NULL) == NULL) {
    errorPrint ("SCOTCH_graphTabLoad: out of memory (1)");
    return (1);
  }

  mappsortflag = 1;                               /* Assume mapping data sorted */
  for (mappfilenum = 0; mappfilenum < mappfilenbr; mappfilenum ++) {
    if (fscanf (stream, GNUMSTRING GNUMSTRING,
                &mappsortval,
                &mappfileval) != 2) {
      errorPrint ("SCOTCH_graphTabLoad: bad input (2)");
      memFree    (mappfiletab);                   /* Free group leader */
      return     (1);
    }
    mappsorttab[mappfilenum].labl = mappsortval;
    mappsorttab[mappfilenum].num  = mappfilenum;
    mappfiletab[mappfilenum]      = mappfileval;

    if ((mappfilenum > 0) &&                      /* Check if mapping data sorted */
        (mappsorttab[mappfilenum].labl < mappsorttab[mappfilenum - 1].labl))
      mappsortflag = 0;                           /* Mapping data not sorted */
  }
  if (mappsortflag != 1)                          /* If mapping data not sorted                      */
    intSort2asc1 (mappsorttab, mappfilenbr);      /* Sort area by ascending labels                   */
  for (mappfilenum = 1; mappfilenum < mappfilenbr; mappfilenum ++) { /* Check mapping data integrity */
    if (mappsorttab[mappfilenum].labl == mappsorttab[mappfilenum - 1].labl) {
      errorPrint ("SCOTCH_graphTabLoad: duplicate vertex label");
      memFree    (mappfiletab);                   /* Free group leader */
      return     (1);
    }
  }

  if (grafptr->vlbltax != NULL) {                 /* If graph has vertex labels */
    vertsortflag = 1;                             /* Assume graph data sorted   */
    for (vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++) {
      vertsorttab[vertnum].labl = grafptr->vlbltax[vertnum];
      vertsorttab[vertnum].num  = vertnum;
      if ((vertnum > 0) &&                        /* Check if graph data sorted */
          (vertsorttab[vertnum].labl < vertsorttab[vertnum - 1].labl))
        vertsortflag = 0;                         /* Graph data not sorted */
    }
    if (vertsortflag != 1)                        /* If graph data not sorted             */
      intSort2asc1 (vertsorttab, grafptr->vertnbr); /* Sort sort area by ascending labels */
  }
  else {                                          /* Graph does not have vertex labels */
    for (vertnum = 0; vertnum < grafptr->vertnbr; vertnum ++) {
      vertsorttab[vertnum].labl = vertnum + mappsorttab[0].labl; /* Use first index as base value */
      vertsorttab[vertnum].num  = vertnum;
    }
  }

  for (vertnum = 0, mappfilenum = 0; vertnum < grafptr->vertnbr; vertnum ++) { /* For all vertices in graph */
    while ((mappfilenum < mappfilenbr) && (mappsorttab[mappfilenum].labl < vertsorttab[vertnum].labl))
      mappfilenum ++;                             /* Search mapping vertex with same label */
    if ((mappfilenum >= mappfilenbr) || (mappsorttab[mappfilenum].labl > vertsorttab[vertnum].labl)) /* If label does not exist */
      continue;                                   /* This vertex has no related mapping data */
    ((Anum *) parttab)[vertsorttab[vertnum].num] = mappfiletab[mappsorttab[mappfilenum ++].num];
  }

  memFree (mappfiletab);                          /* Free group leader */

  return (0);
}

/*+ This routine loads the contents of
*** the given user mapping from the 
*** given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapLoad (
const SCOTCH_Graph * const    actgrafptr,         /*+ Graph to map    +*/
const SCOTCH_Mapping * const  mappptr,            /*+ Mapping to save +*/
FILE * const                  stream)             /*+ Output stream   +*/
{
  Graph *               grafptr;
  LibMapping * restrict lmapptr;

  grafptr = (Graph *) actgrafptr;
  lmapptr = (LibMapping *) mappptr;
#ifdef SCOTCH_DEBUG_GRAPH2
  if (grafptr != lmapptr->grafptr) {
    errorPrint ("SCOTCH_graphMapLoad: mapping structure must derive from graph");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  if (lmapptr->parttab == NULL) {                 /* Allocate array if necessary */
    if ((lmapptr->parttab = (Gnum *) memAlloc (grafptr->vertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("SCOTCH_graphMapLoad: out of memory");
      return     (1);
    }
    lmapptr->flagval |= LIBMAPPINGFREEPART;       /* As the array was not allocated by the user, it will be freed */
  }

  return (SCOTCH_graphTabLoad (actgrafptr, (SCOTCH_Num *) lmapptr->parttab, stream));
}

/*+ This routine saves the contents of
*** the given user mapping to the given
*** stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapSave (
const SCOTCH_Graph * const    actgrafptr,         /*+ Graph to map    +*/
const SCOTCH_Mapping * const  mappptr,            /*+ Mapping to save +*/
FILE * const                  stream)             /*+ Output stream   +*/
{
  Graph *               grafptr;
  LibMapping * restrict lmapptr;
  Gnum                  vertnum;
  Gnum                  baseval;

  grafptr = (Graph *) actgrafptr;
  lmapptr = (LibMapping *) mappptr;
#ifdef SCOTCH_DEBUG_GRAPH2
  if (grafptr != lmapptr->grafptr) {
    errorPrint ("SCOTCH_graphMapSave: mapping structure must derive from graph");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  const Gnum * restrict vlbltax = grafptr->vlbltax;
  const Gnum * restrict parttab = lmapptr->parttab;

  if (fprintf (stream, GNUMSTRING "\n",
               (Gnum) grafptr->vertnbr) == EOF) {
    errorPrint ("SCOTCH_graphMapSave: bad output (1)");
    return     (1);
  }

  baseval = grafptr->baseval;
  for (vertnum = baseval; vertnum < grafptr->vertnnd; vertnum ++) {
    if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                 (Gnum) ((vlbltax != NULL) ? vlbltax[vertnum] : vertnum),
                 (Gnum) parttab[vertnum - baseval]) == EOF) {
      errorPrint ("SCOTCH_graphMapSave: bad output (2)");
      return     (1);
    }
  }

  return (0);
}
