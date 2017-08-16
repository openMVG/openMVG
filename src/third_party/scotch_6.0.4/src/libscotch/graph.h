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
/**   NAME       : graph.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source graph functions.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     18 aug 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     28 nov 1995     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     15 sep 1998     **/
/**                # Version 3.3  : from : 28 sep 1998     **/
/**                                 to     23 mar 1999     **/
/**                # Version 3.4  : from : 20 mar 2000     **/
/**                                 to     20 mar 2000     **/
/**                # Version 4.0  : from : 24 nov 2001     **/
/**                                 to     03 mar 2006     **/
/**                # Version 5.0  : from : 03 mar 2006     **/
/**                                 to     01 jun 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     09 aug 2014     **/
/**                                                        **/
/************************************************************/

#define GRAPH_H

/*
**  The defines.
*/

/*+ Graph option flags. +*/

#define GRAPHNONE                   0x0000        /* No options set */

#define GRAPHFREEEDGE               0x0001        /* Free edgetab array        */
#define GRAPHFREEVERT               0x0002        /* Free verttab array        */
#define GRAPHFREEVNUM               0x0004        /* Free vnumtab array        */
#define GRAPHFREEOTHR               0x0008        /* Free all other arrays     */
#define GRAPHFREETABS               0x000F        /* Free all graph arrays     */
#define GRAPHVERTGROUP              0x0010        /* All vertex arrays grouped */
#define GRAPHEDGEGROUP              0x0020        /* All edge arrays grouped   */

#define GRAPHBITSUSED               0x003F        /* Significant bits for plain graph routines               */
#define GRAPHBITSNOTUSED            0x0040        /* Value above which bits not used by plain graph routines */

#define GRAPHIONOLOADVERT           1             /* Remove vertex loads on loading */
#define GRAPHIONOLOADEDGE           2             /* Remove edge loads on loading   */

/*
**  The type and structure definitions.
*/

#ifndef GNUMMAX                                   /* If dgraph.h not included    */
typedef INT                   Gnum;               /* Vertex and edge numbers     */
typedef UINT                  Gunum;              /* Unsigned type of same width */
#define GNUMMAX                     (INTVALMAX)   /* Maximum signed Gnum value   */
#define GNUMSTRING                  INTSTRING     /* String to printf a Gnum     */
#endif /* GNUMMAX */

/*+ The vertex part type, in compressed form. +*/

typedef byte GraphPart;

/*+ The vertex list structure. Since a vertex list
    always refers to a given graph, vertex indices
    contained in the vertex list array are based with
    respect to the base value of the associated graph.
    However, the array itself is not based.            +*/

typedef struct VertList_ {
  Gnum                      vnumnbr;              /*+ Number of vertices in list +*/
  Gnum *                    vnumtab;              /*+ Pointer to vertex array    +*/
} VertList;

/*+ The graph flag type. +*/

typedef int GraphFlag;                            /*+ Graph property flags +*/

/*+ The graph parallel context structure. +*/

typedef struct GraphProc_ {
#ifdef SCOTCH_PTSCOTCH
  MPI_Comm                  proccomm;             /*+ Communicator used for parallel algorithm +*/
  int                       procglbnbr;           /*+ Number of processes in communicator      +*/
  int                       proclocnum;           /*+ Rank of process in current communicator  +*/
#endif /* SCOTCH_PTSCOTCH */
} GraphProc;

/*+ The graph structure. +*/

typedef struct Graph_ {
  GraphFlag                 flagval;              /*+ Graph properties                          +*/
  Gnum                      baseval;              /*+ Base index for edge/vertex arrays         +*/
  Gnum                      vertnbr;              /*+ Number of vertices in graph               +*/
  Gnum                      vertnnd;              /*+ Number of vertices in graph, plus baseval +*/
  Gnum *                    verttax;              /*+ Vertex array [based]                      +*/
  Gnum *                    vendtax;              /*+ End vertex array [based]                  +*/
  Gnum *                    velotax;              /*+ Vertex load array (if present)            +*/
  Gnum                      velosum;              /*+ Overall graph vertex load                 +*/
  Gnum *                    vnumtax;              /*+ Vertex number in ancestor graph           +*/
  Gnum *                    vlbltax;              /*+ Vertex label (from file)                  +*/
  Gnum                      edgenbr;              /*+ Number of edges (arcs) in graph           +*/
  Gnum *                    edgetax;              /*+ Edge array [based]                        +*/
  Gnum *                    edlotax;              /*+ Edge load array (if present)              +*/
  Gnum                      edlosum;              /*+ Sum of edge (in fact arc) loads           +*/
  Gnum                      degrmax;              /*+ Maximum degree                            +*/
  GraphProc *               procptr;              /*+ Pointer to parallel context (if any)      +*/
} Graph;

/*
**  The function prototypes.
*/

#ifndef GRAPH
#define static
#endif

int                         listInit            (VertList *);
void                        listExit            (VertList *);
int                         listAlloc           (VertList *, Gnum);
int                         listFree            (VertList *);
int                         listLoad            (VertList *, FILE *);
int                         listSave            (VertList *, FILE *);
void                        listSort            (VertList *);
int                         listCopy            (VertList *, VertList *);

int                         graphInit           (Graph * const);
void                        graphExit           (Graph * const);
void                        graphFree           (Graph * const);
Gnum                        graphBase           (Graph * const, const Gnum);
int                         graphBand           (const Graph * restrict const, const Gnum, Gnum * restrict const, const Gnum, Gnum * restrict * restrict const, Gnum * restrict const, Gnum * restrict const, Gnum * restrict const, const Gnum * restrict const, Gnum * restrict const);
int                         graphCheck          (const Graph *);
int                         graphInduceList     (const Graph * const, const VertList * const, Graph * const);
int                         graphInducePart     (const Graph * const, const GraphPart *, const Gnum, const GraphPart, Graph * const);
int                         graphLoad           (Graph * const, FILE * const, const Gnum, const GraphFlag);
int                         graphLoad2          (const Gnum, const Gnum, const Gnum * const, const Gnum * const, Gnum * restrict const, const Gnum, const Gnum * const);
int                         graphSave           (const Graph * const, FILE * const);

#ifdef GEOM_H
int                         graphGeomLoadChac   (Graph * restrict const, Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomSaveChac   (const Graph * restrict const, const Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomLoadHabo   (Graph * restrict const, Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomLoadMmkt   (Graph * restrict const, Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomSaveMmkt   (const Graph * restrict const, const Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomLoadScot   (Graph * restrict const, Geom * restrict const, FILE * const, FILE * const, const char * const);
int                         graphGeomSaveScot   (const Graph * restrict const, const Geom * restrict const, FILE * const, FILE * const, const char * const);
#endif /* GEOM_H */

#undef static
