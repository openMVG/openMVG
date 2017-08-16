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
/**   NAME       : gout_c.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a result viewer.                **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 06 oct 1994     **/
/**                                 to     23 dec 1994     **/
/**                # Version 3.0  : from : 14 jul 1995     **/
/**                                 to     11 oct 1995     **/
/**                # Version 3.1  : from : 27 mar 1996     **/
/**                                 to     03 apr 1996     **/
/**                # Version 3.2  : from : 02 dec 1996     **/
/**                                 to     05 jun 1998     **/
/**                # Version 3.3  : from : 29 may 1999     **/
/**                                 to   : 03 jun 1999     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     08 feb 2004     **/
/**                # Version 5.0  : from : 25 may 2007     **/
/**                                 to     25 may 2007     **/
/**                # Version 5.1  : from : 25 oct 2007     **/
/**                                 to     14 feb 2011     **/
/**                # Version 6.0  : from : 16 oct 2010     **/
/**                                 to     12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes
*/

#define GOUT

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gout_c.h"
#include "gout_o.h"

/*
**  The static and global variables
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
File                        C_fileTab[C_FILENBR] = { /* The file array; public  */
                              { "r" },
                              { "r" },
                              { "r" },
                              { "w" } };

static unsigned int         C_geoFlag = C_GEOFLAGDEFAULT; /* Geometry flag */

static const char *         C_usageList[] = {     /* Usage */
  "gout [<input source file> [<input geometry file> [<input mapping file> [<output picture file>]]]] <options>",
  "  -g<arguments>       : Geometry parameters :",
  "                          n  : do not read geometry data (matrix display)",
  "                          p  : permute Y and Z geometry dimensions",
  "                          r  : rotate geometry by 90 degrees",
  "  -h                  : Display this help",
  "  -mn                 : Do not read mapping data",
  "  -Oi[{<arguments>}]  : Open Inventor mesh file :",
  "                          c        : color output",
  "                          g        : gray level output",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "  -Om[{<arguments>}]  : PostScript matrix file :",
  "                          e        : EPSF-type output",
  "                          f        : full-page output",
  "  -Op[{<arguments>}]  : PostScript mesh file :",
  "                          c        : color output",
  "                          g        : gray level output",
  "                          e        : EPSF-type output",
  "                          f        : full-page output",
  "                          s        : short clipping (disks excluded)",
  "                          l        : large clipping (disks included)",
  "                          a        : avoid displaying disks",
  "                          d        : display disks",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "                          X=<val>  : maximum x clipping ratio (in [0.0;1.0])",
  "                          x=<val>  : minimum x clipping ratio",
  "                          Y=<val>  : maximum y clipping ratio",
  "                          y=<val>  : minimum y clipping ratio",
  "  -Ot[{<arguments>}]  : Tulip graph file :",
  "                          b        : b/w output",
  "                          c        : color output",
  "                          a        : avoid displaying disks",
  "                          d        : display disks",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "  -V                  : Print program version and copyright",
  "",
  "Default option set is : -Oi{c,v}",
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
  C_Graph            grafdat;                     /* Source graph   */
  C_Geometry         geo;                         /* Graph geometry */
  C_Mapping          map;                         /* Result mapping */
  int                i, j;

  errorProg ("gout");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'G' :                                /* Geometry parameters */
        case 'g' :
          if ((j = C_geoParse (&argv[i][2])) != 0)
            errorPrint ("main: error in geometry option string '%d'", j);
          break;
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'M' :                                /* No-mapping flag */
        case 'm' :
          if (((argv[i][2] != 'N') && (argv[i][2] != 'n')) || (argv[i][3] != '\0'))
            errorPrint ("main: error in mapping option string '%s'", &argv[i][2]);
          C_filenamemapinp = "-";                 /* Default name to avoid opening   */
          C_filepntrmapinp = NULL;                /* NULL file pointer means no file */
          break;
        case 'O' :                                /* Output parameters */
        case 'o' :
          if ((j = outDrawParse (&argv[i][2])) != 0)
            errorPrint ("main: error in output option string (%d)", j);
          break;
        case 'V' :
          fprintf (stderr, "gout, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, "Copyright 2004,2007,2008,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        default :
          errorPrint ("main: Unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit (&grafdat.grafdat);            /* Create graph structure         */
  SCOTCH_graphLoad (&grafdat.grafdat, C_filepntrsrcinp, 0, 3); /* Read source graph */
  SCOTCH_graphData (&grafdat.grafdat, &grafdat.baseval,
                    &grafdat.vertnbr, &grafdat.verttab, &grafdat.vendtab, NULL, &grafdat.vlbltab,
                    &grafdat.edgenbr, &grafdat.edgetab, NULL);

  C_geoInit (&geo, &grafdat);                     /* Create geometry structure */
  if (C_geoFlag & C_GEOFLAGUSE)                   /* If geometry is wanted     */
    C_geoLoad (&geo, C_filepntrgeoinp);           /* Read graph geometry       */

  C_mapInit (&map, &grafdat);                     /* Create mapping structure */
  C_mapLoad (&map, C_filepntrmapinp);             /* Read result mapping      */

  outDraw (&grafdat, &geo, &map, C_filepntrdatout); /* Build and write the output */

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  C_mapExit        (&map);                        /* Free data structures */
  C_geoExit        (&geo);
  SCOTCH_graphExit (&grafdat.grafdat);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);
}

/***********************************/
/*                                 */
/* These routines handle geometry. */
/*                                 */
/***********************************/

/* This routine parses the source graph
** option string.
** It returns:
** - 0  : if string successfully scanned.
** - 1  : if invalid options
** - 2  : if invalid option arguments.
** - 3  : if syntax error in string.
*/

int
C_geoParse (
const char * const          string)
{
  const char *        cptr;

  for (cptr = string; ; cptr ++) {
    switch (*cptr) {
      case 'N' :                                  /* Do not read geometry data */
      case 'n' :
        C_geoFlag &= ~C_GEOFLAGUSE;
        break;
      case 'P' :                                  /* Permute Y and Z */
      case 'p' :
        C_geoFlag |= C_GEOFLAGPERMUT;
        break;
      case 'R' :                                  /* If want to rotate */
      case 'r' :
        C_geoFlag |= C_GEOFLAGROTATE;
        break;
      case '\0' :
        return (0);
      default   :
        return (1);
    }
  }
}

/* This routine creates a geometry with
** respect to a given source graph.
** It returns:
** - VOID  : in all cases.
*/

void
C_geoInit (
C_Geometry * const          geomptr,
const C_Graph * const       grafptr)
{
  geomptr->grafptr = grafptr;
  geomptr->verttab = NULL;
}

/* This routine deletes a geometry.
** It returns:
** - VOID  : in all cases.
*/

void
C_geoExit (
C_Geometry * const          geomptr)
{
  if (geomptr->verttab != NULL)                   /* If there is a geometry array */
    memFree (geomptr->verttab);                   /* Free it                      */
}

/* This routine loads a mapping.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

/* This is the comparison function used by the
   quicksort algorithm, to sort by increasing
   labels.                                     */

static
int
C_geoLoad2 (
const C_VertSort * const    vert0,
const C_VertSort * const    vert1)
{
  return ((vert0->labl > vert1->labl) ?  1 : -1);
}

/** This is the loading routine. **/

int
C_geoLoad (
C_Geometry * const          geomptr,
FILE * const                stream)
{
  C_VertSort *        vertsorttab;                /* Pointer to graph sorting array             */
  int                 vertsortflag;               /* Flag set if graph data sorted by label     */
  C_VertSort *        geomsorttab;                /* Pointer to geometric data sorting array    */
  int                 geomsortflag;               /* Flag set if geometric data sorted by label */
  int                 geomfiletype;               /* Type of geometry file                      */
  SCOTCH_Num          geomfilenbr;                /* Number of geometric coordinates in file    */
  SCOTCH_Num          geomfileval;                /* Value of maximum size for compatibility    */
  C_GeoVert *         geomfiletab;                /* Pointer to geometric data read from file   */
  SCOTCH_Num          vertlablval;                /* Value of maximum size for compatibility    */
  SCOTCH_Num          i, j;
  int                 o;

  if ((geomptr->verttab == NULL) &&               /* Allocate geometry if necessary */
      ((geomptr->verttab = (C_GeoVert *) memAlloc (geomptr->grafptr->vertnbr * sizeof (C_GeoVert))) == NULL)) {
    errorPrint ("C_geoLoad: out of memory (1)");
    return     (1);
  }

  if ((fscanf (stream, "%d" SCOTCH_NUMSTRING,     /* Read type and number of geometry items */
               &geomfiletype,
               &geomfileval) != 2) ||
      (geomfiletype < 1)           ||
      (geomfiletype > 3)           ||
      (geomfileval  < 1)) {
    errorPrint ("C_geoLoad: bad input (1)");
    return     (1);
  }
  geomfilenbr = (SCOTCH_Num) geomfileval;

  if (memAllocGroup ((void **) (void *)
                     &geomfiletab, (size_t) (geomfilenbr               * sizeof (C_GeoVert)),
                     &geomsorttab, (size_t) (geomfilenbr               * sizeof (C_VertSort)),
                     &vertsorttab, (size_t) (geomptr->grafptr->vertnbr * sizeof (C_VertSort)), NULL) == NULL) {
    errorPrint ("C_geoLoad: out of memory (2)");
    return     (1);
  }

  o = 0;
  geomsortflag = 1;                               /* Assume geometry data sorted */
  switch (geomfiletype) {
    case 1 :                                      /* Load 2D coordinates array */
      for (i = 0; (i < geomfilenbr) && (o == 0); i ++) {
        if (fscanf (stream, SCOTCH_NUMSTRING "%lf",
                    &vertlablval,
                    &geomfiletab[i].x) != 2)
          o = 1;
        geomsorttab[i].labl = (SCOTCH_Num) vertlablval;
        geomfiletab[i].y    =                     /* No Y and Z coordinates */
        geomfiletab[i].z    = 0.0;
        geomsorttab[i].num  = i;

        if (C_geoFlag & C_GEOFLAGROTATE) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                = geomfiletab[i].y;
          geomfiletab[i].y = geomfiletab[i].x;
          geomfiletab[i].x = - t;
        }

        if ((i > 0) &&                            /* Check if geometry data sorted */
            (geomsorttab[i].labl < geomsorttab[i - 1].labl))
          geomsortflag = 0;                       /* Geometry data not sorted */
      }
      break;
    case 2 :                                      /* Load 2D coordinates array */
      for (i = 0; (i < geomfilenbr) && (o == 0); i ++) {
        if (fscanf (stream, SCOTCH_NUMSTRING "%lf%lf",
                    &vertlablval,
                    &geomfiletab[i].x,
                    &geomfiletab[i].y) != 3)
          o = 1;
        geomsorttab[i].labl = (SCOTCH_Num) vertlablval;
        geomfiletab[i].z    = 0.0;                /* No Z coordinate */
        geomsorttab[i].num  = i;

        if (C_geoFlag & C_GEOFLAGROTATE) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                = geomfiletab[i].y;
          geomfiletab[i].y = geomfiletab[i].x;
          geomfiletab[i].x = - t;
        }

        if ((i > 0) &&                            /* Check if geometry data sorted */
            (geomsorttab[i].labl < geomsorttab[i - 1].labl))
          geomsortflag = 0;                       /* Geometry data are not sorted */
      }
      break;
    case 3 :                                      /* Load 3D coordinates array */
      for (i = 0; (i < geomfilenbr) && (o == 0); i ++) {
        if (fscanf (stream, SCOTCH_NUMSTRING "%lf%lf%lf",
                    &vertlablval,
                    &geomfiletab[i].x,
                    &geomfiletab[i].y,
                    &geomfiletab[i].z) != 4)
          o = 1;
        geomsorttab[i].labl = (SCOTCH_Num) vertlablval;
        geomsorttab[i].num  = i;

        if (C_geoFlag & C_GEOFLAGPERMUT) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                = geomfiletab[i].z;
          geomfiletab[i].z = geomfiletab[i].y;
          geomfiletab[i].y = t;
        }
        if ((i > 0) &&                            /* Check if geometry data sorted */
            (geomsorttab[i].labl < geomsorttab[i - 1].labl))
          geomsortflag = 0;                       /* Geometry data not sorted */
      }
      break;
    default :
      errorPrint ("C_geoLoad: invalid geometry type (%d)", geomfiletype);
      memFree    (geomfiletab);                   /* Free group leader */
      return     (1);
  }
  if (o != 0) {
    errorPrint ("C_geoLoad: bad input (2)");
    memFree    (geomfiletab);                     /* Free group leader */
    return     (1);
  }

  if (geomsortflag != 1)                          /* If geometry data not sorted        */
    qsort ((char *) geomsorttab, geomfilenbr,     /* Sort sort area by ascending labels */
           sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  for (i = 1; i < geomfilenbr; i ++) {            /* Check geometric data integrity */
    if (geomsorttab[i].labl == geomsorttab[i - 1].labl) {
      errorPrint ("C_geoLoad: duplicate vertex label");
      memFree    (geomfiletab);                   /* Free group leader */
      return     (1);
    }
  }

  if (geomptr->grafptr->vlbltab != NULL) {        /* If graph has vertex labels */
    vertsortflag = 1;                             /* Assume graph data sorted   */
    for (i = 0; i < geomptr->grafptr->vertnbr; i ++) {
      vertsorttab[i].labl = geomptr->grafptr->vlbltab[i];
      vertsorttab[i].num  = i;
      if ((i > 0) &&                              /* Check if graph data sorted */
          (vertsorttab[i].labl < vertsorttab[i - 1].labl))
        vertsortflag = 0;                         /* Graph data not sorted */
    }
    if (vertsortflag != 1)                        /* If graph data not sorted                       */
      qsort ((char *) vertsorttab, geomptr->grafptr->vertnbr, /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  }
  else {                                          /* Graph does not have vertex labels */
    for (i = 0; i < geomptr->grafptr->vertnbr; i ++) {
      vertsorttab[i].labl = i + geomsorttab[0].labl; /* Use first index as base value */
      vertsorttab[i].num  = i;
    }
  }

  for (i = 0, j = 0; i < geomptr->grafptr->vertnbr; i ++) { /* For all vertices in graph */
    while ((j < geomfilenbr) && (geomsorttab[j].labl < vertsorttab[i].labl))
      j ++;                                       /* Search geometry vertex with same label             */
    if ((j >= geomfilenbr) || (geomsorttab[j].labl > vertsorttab[i].labl)) { /* If label does not exist */
      errorPrint ("C_geoLoad: vertex geometry data not found for label '" SCOTCH_NUMSTRING "'",
                  vertsorttab[i].labl);
      memFree    (geomfiletab);                   /* Free group leader */
      return     (1);
    }
    geomptr->verttab[vertsorttab[i].num] = geomfiletab[geomsorttab[j ++].num];
  }

  memFree (geomfiletab);                          /* Free group leader */

  return (0);
}

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine creates a mapping with
** respect to a given source graph.
** It returns:
** - VOID  : in all cases.
*/

void
C_mapInit (
C_Mapping * const           mapptr,
const C_Graph * const       grafptr)
{
  mapptr->grafptr = grafptr;
  mapptr->labltab = NULL;
}

/* This routine deletes a mapping.
** It returns:
** - VOID  : in all cases.
*/

void
C_mapExit (
C_Mapping * const           mapptr)
{
  if (mapptr->labltab != NULL)                    /* If there is a domain array */
    memFree (mapptr->labltab);                    /* Free it                    */
}

/* This routine loads a mapping.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
C_mapLoad (
C_Mapping * const           mapptr,
FILE * const                stream)
{
  C_VertSort *        vertsorttab;                /* Pointer to graph sorting array           */
  int                 vertsortflag;               /* Flag set if graph data sorted by label   */
  C_VertSort *        mapsorttab;                 /* Pointer to mapping data sorting array    */
  int                 mapsortflag;                /* Flag set if mapping data sorted by label */
  SCOTCH_Num          mapsortval;                 /* Value of maximum size for compatibility  */
  SCOTCH_Num          mapfileval;                 /* Value of maximum size for compatibility  */
  SCOTCH_Num          mapfilenbr;                 /* Number of mapping pairs in file          */
  SCOTCH_Num *        mapfiletab;                 /* Pointer to mapping data read from file   */
  SCOTCH_Num          i, j;

  if ((mapptr->labltab == NULL) &&                /* Allocate array if necessary */
      ((mapptr->labltab = (SCOTCH_Num *) memAlloc (mapptr->grafptr->vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    errorPrint ("C_mapLoad: out of memory (1)");
    return     (1);
  }

  memset (mapptr->labltab, ~0, mapptr->grafptr->vertnbr * sizeof (SCOTCH_Num)); /* Pre-initialize mapping */

  if (stream == NULL)                             /* If stream is invalid */
    return (0);

  if ((fscanf (stream, SCOTCH_NUMSTRING,          /* Read number of mapping pairs */
               &mapfileval) != 1) ||
      (mapfileval < 1)) {
    errorPrint ("C_mapLoad: bad input (1)");
    return     (1);
  }
  mapfilenbr = (SCOTCH_Num) mapfileval;

  if (memAllocGroup ((void **) (void *)
                     &mapfiletab,  (size_t) (mapfilenbr               * sizeof (SCOTCH_Num)),
                     &mapsorttab,  (size_t) (mapfilenbr               * sizeof (C_VertSort)),
                     &vertsorttab, (size_t) (mapptr->grafptr->vertnbr * sizeof (C_VertSort)), NULL) == NULL) {
    errorPrint ("C_mapLoad: out of memory (2)");
    return     (1);
  }

  mapsortflag = 1;                                /* Assume mapping data sorted */
  for (i = 0; i < mapfilenbr; i ++) {
    if (fscanf (stream, SCOTCH_NUMSTRING SCOTCH_NUMSTRING,
                &mapsortval,
                &mapfileval) != 2) {
      errorPrint ("C_mapLoad: bad input (2)");
      memFree    (mapfiletab);                    /* Free group leader */
      return     (1);
    }
    mapsorttab[i].labl = mapsortval;
    mapsorttab[i].num  = i;
    mapfiletab[i]      = mapfileval;

    if ((i > 0) &&                                /* Check if mapping data sorted */
        (mapsorttab[i].labl < mapsorttab[i - 1].labl))
      mapsortflag = 0;                            /* Mapping data not sorted */
  }
  if (mapsortflag != 1)                           /* If mapping data not sorted         */
      qsort ((char *) mapsorttab, mapfilenbr,     /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  for (i = 1; i < mapfilenbr; i ++) {             /* Check mapping data integrity */
    if (mapsorttab[i].labl == mapsorttab[i - 1].labl) {
      errorPrint ("C_mapLoad: duplicate vertex label");
      memFree    (mapfiletab);                    /* Free group leader */
      return     (1);
    }
  }

  if (mapptr->grafptr->vlbltab != NULL) {         /* If graph has vertex labels */
    vertsortflag = 1;                             /* Assume graph data sorted   */
    for (i = 0; i < mapptr->grafptr->vertnbr; i ++) {
      vertsorttab[i].labl = mapptr->grafptr->vlbltab[i];
      vertsorttab[i].num  = i;
      if ((i > 0) &&                              /* Check if graph data sorted */
          (vertsorttab[i].labl < vertsorttab[i - 1].labl))
        vertsortflag = 0;                         /* Graph data not sorted */
    }
    if (vertsortflag != 1)                        /* If graph data not sorted                      */
      qsort ((char *) vertsorttab, mapptr->grafptr->vertnbr, /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  }
  else {                                          /* Graph does not have vertex labels */
    for (i = 0; i < mapptr->grafptr->vertnbr; i ++) {
      vertsorttab[i].labl = i + mapptr->grafptr->baseval;
      vertsorttab[i].num  = i;
    }
  }

  for (i = 0, j = 0; i < mapptr->grafptr->vertnbr; i ++) { /* For all vertices in graph */
    while ((j < mapfilenbr) && (mapsorttab[j].labl < vertsorttab[i].labl))
      j ++;                                       /* Search mapping vertex with same label          */
    if ((j >= mapfilenbr) || (mapsorttab[j].labl > vertsorttab[i].labl)) /* If label does not exist */
      continue;                                   /* This vertex has no related mapping data        */
    mapptr->labltab[vertsorttab[i].num] = mapfiletab[mapsorttab[j ++].num];
  }

  memFree (mapfiletab);                           /* Free group leader */

  return (0);
}

/**************************************/
/*                                    */
/* The option string parsing routine. */
/*                                    */
/**************************************/

/* This routine parses an option string.
** It returns:
** - 0 : if string successfully scanned.
** - 1 : if invalid code name.
** - 2 : if invalid arguments for the code.
** - 3 : if syntax error in string.
*/

int
C_parse (
const C_ParseCode * const   codeptr,             /* Pointer to the code array          */
const C_ParseArg * const    argptr,              /* Pointer to the code argument array */
int * const                 codeval,             /* Pointer to the code value to set   */
char * const                string)              /* Pointer to the string to parse     */
{
  int                 code;                      /* Code found                       */
  int                 codelen;                   /* Code name length                 */
  char                argbuf[128];               /* Buffer for argument scanning     */
  int                 arglen;                    /* Length of the current argument   */
  char *              argbeg;                    /* Pointer to beginning of argument */
  char *              argend;                    /* Pointer to end of argument       */
  char *              argequ;                    /* Position of the '=' character    */
  int                 i, j;

  code    =
  codelen = 0;                                   /* No code recognized yet              */
  for (i = 0; codeptr[i].name != NULL; i ++) {   /* For all the codes                   */
    if ((strncasecmp (string,                    /* Find the longest matching code name */
                      codeptr[i].name,
                      j = strlen (codeptr[i].name)) == 0) &&
        (j > codelen)) {
      code    = codeptr[i].code;
      codelen = j;
    }
  }
  if (codelen == 0)                              /* If no code recognized  */
    return (1);                                  /* Return the error value */
  *codeval = code;                               /* Set the code value     */

  argbeg = string + codelen;                     /* Point to the end of the code name */
  if (*argbeg == '{') {                          /* If there are arguments            */
    argbeg ++;                                   /* Point to argument beginning       */
    do {                                         /* For all arguments                 */
      argend = strpbrk (argbeg, ",}\0");         /* Search for the argument end       */
      if (*argend == '\0')                       /* If there is no end delimiter      */
        return (3);                              /* Return the syntax error value     */

      arglen = ((argend - argbeg) < 127)         /* Get argument bounded length */
               ? (argend - argbeg)
               : 127;
      strncpy (argbuf, argbeg, arglen);          /* Copy the argument to the buffer */
      argbuf[arglen] = '\0';                     /* Mark the end of the argument    */
      argequ = strpbrk (argbuf, "=");            /* Search for the '=' character    */
      if (argequ != NULL)                        /* If it exists                    */
        *argequ++ = '\0';                        /* Turn it into a separating null  */

      for (i = 0, j = -1; argptr[i].name != NULL; i ++) { /* Scan all the possible arguments */
        if ((argptr[i].code == code) &&          /* If the proper name is found              */
            (strcmp (argbuf, argptr[i].name) == 0)) {
          j = i;                                 /* Record the position */
          break;                                 /* Exit the loop       */
        }
      }
      if (j == -1)                               /* If invalid argument     */
        return (2);                              /* Return the proper value */

      if (argptr[j].format != NULL) {            /* If there is a value to read    */
        if (argequ == NULL)                      /* If none has been given however */
          return (2);                            /* Return the error value         */
        if (sscanf (argequ,                      /* Try to read the argument value */
                    argptr[j].format,
                    argptr[j].ptr) != 1)
          return (2);                            /* Return if error                */
        if (argptr[j].func != NULL)              /* If there is a control function */
          if (argptr[j].func (argptr[j].ptr) != 0) /* If the function fails        */
            return (2);                          /* Return the error value         */
      }
      else {                                     /* If no value needed           */
        if (argequ != NULL)                      /* If there is one however      */
          return (2);                            /* Return the error code        */
        *((char *) argptr[j].ptr) = argbuf[0];   /* Copy the first argument char */
      }

      argbeg = argend + 1;                       /* Skip the processed argument         */
    } while (*argend != '}');                    /* Loop as long as there are arguments */
  }

  return ((*argbeg == '\0') ? 0 : 3);            /* Check if no extraneous characters */
}
