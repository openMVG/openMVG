/* Copyright 2004,2007,2008,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gout_c.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a result viewer.                **/
/**                This module contains the data declara-  **/
/**                tions for the main module.              **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 06 oct 1994     **/
/**                                 to     01 nov 1994     **/
/**                # Version 3.0  : from : 14 jul 1995     **/
/**                                 to     02 oct 1995     **/
/**                # Version 3.2  : from : 02 dec 1996     **/
/**                                 to     05 jun 1998     **/
/**                # Version 3.3  : from : 01 jun 1999     **/
/**                                 to     01 jun 1999     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     11 dec 2001     **/
/**                # Version 5.0  : from : 13 dec 2007     **/
/**                                 to     15 mar 2008     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   4             /* Number of files in list                */
#define C_FILEARGNBR                4             /* Number of files which can be arguments */

#define C_filenamesrcinp            fileBlockName (C_fileTab, 0) /* Source graph file name          */
#define C_filenamegeoinp            fileBlockName (C_fileTab, 1) /* Source graph geometry file name */
#define C_filenamemapinp            fileBlockName (C_fileTab, 2) /* Mapping result file name        */
#define C_filenamedatout            fileBlockName (C_fileTab, 3) /* Output data file name           */

#define C_filepntrsrcinp            fileBlockFile (C_fileTab, 0) /* Source graph input file    */
#define C_filepntrgeoinp            fileBlockFile (C_fileTab, 1) /* Source graph geometry file */
#define C_filepntrmapinp            fileBlockFile (C_fileTab, 2) /* Mapping result input file  */
#define C_filepntrdatout            fileBlockFile (C_fileTab, 3) /* Data output file           */

#define C_filemodemapinp            fileBlockMode (C_fileTab, 2) /* Mapping result mode */

/*+ Dimension definitions. +*/

#define x                           c[0]
#define y                           c[1]
#define z                           c[2]

/*+ Geometry flags. +*/

#define C_GEOFLAGDEFAULT            0x0001        /* Default geometry flag            */
#define C_GEOFLAGUSE                0x0001        /* Use geometry                     */
#define C_GEOFLAGROTATE             0x0002        /* Rotate the picture by 90 degrees */
#define C_GEOFLAGPERMUT             0x0004        /* Permute Y and Z dimensions       */

/*
**  The type and structure definitions.
*/

/*+ This structure defines a source graph. +*/

typedef struct C_Graph_ {
  SCOTCH_Graph            grafdat;                /*+ Source graph data  +*/
  SCOTCH_Num              baseval;                /*+ Base value         +*/
  SCOTCH_Num              vertnbr;                /*+ Number of vertices +*/
  SCOTCH_Num *            verttab;                /*+ Vertex array       +*/
  SCOTCH_Num *            vendtab;                /*+ Vertex end array   +*/
  SCOTCH_Num *            vlbltab;                /*+ Vertex label array +*/
  SCOTCH_Num              edgenbr;                /*+ Number of edges    +*/
  SCOTCH_Num *            edgetab;                /*+ Edge array         +*/
} C_Graph;

/*+ This structure defines a geometrical vertex. +*/

typedef struct C_GeoVert_ {
  double                    c[3];                 /*+ Vertex coordinates (x,y,z) +*/
} C_GeoVert;

/*+ This structure defines a geometrical
    mapping which contains the positions
    of the graph vertices.               +*/

typedef struct C_Geometry_ {
  const C_Graph *         grafptr;                /*+ Pointer to source graph      +*/
  C_GeoVert *             verttab;                /*+ Pointer to coordinates array +*/
} C_Geometry;

/*+ This structure defines a domain label
    mapping, which contains the reference
    to the mapping source graph.          +*/

typedef struct C_Mapping_ {
  const C_Graph *         grafptr;                /*+ Pointer to source graph +*/
  SCOTCH_Num *            labltab;                /*+ Pointer to label array  +*/
} C_Mapping;

/*+ The sort structure, used to sort graph vertices by label. +*/

typedef struct C_VertSort_ {
  SCOTCH_Num                labl;                 /*+ Vertex label  +*/
  SCOTCH_Num                num;                  /*+ Vertex number +*/
} C_VertSort;

/*+ This structure is the code
    name array entries.        +*/

typedef struct C_ParseCode_ {
  int                       code;                 /*+ Code value +*/
  char *                    name;                 /*+ Code name  +*/
} C_ParseCode;

/* This structure defines the
   code argument array entries. */

typedef struct C_ParseArg_ {
  const char *              name;                 /*+ Name of the argument                         +*/
  int                       code;                 /*+ Code value                                   +*/
  const char *              format;               /*+ scanf-like format; NULL means char, no value +*/
  const void *              ptr;                  /*+ Pointer to the argument location             +*/
  int                    (* func) ();             /*+ Pointer to the argument test function        +*/
} C_ParseArg;

/*
**  The global data declarations.
*/

extern File                 C_fileTab[C_FILENBR]; /*+ File array +*/

/*
**  The function prototypes.
*/

int                         C_geoParse          (const char * const);
void                        C_geoInit           (C_Geometry * const, const C_Graph * const);
void                        C_geoExit           (C_Geometry * const);
int                         C_geoLoad           (C_Geometry * const, FILE * const);

void                        C_mapInit           (C_Mapping * const, const C_Graph * const);
void                        C_mapExit           (C_Mapping * const);
int                         C_mapLoad           (C_Mapping * const, FILE * const);

int                         C_parse             (const C_ParseCode * const, const C_ParseArg * const, int * const, char * const);
