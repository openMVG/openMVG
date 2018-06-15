/* Copyright 2007-2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                Nicolas GICQUEL (P0.1)                  **/
/**                Jerome LACOSTE (P0.1)                   **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed source graph        **/
/**                structure.                              **/
/**                                                        **/
/**   DATES      : # Version P0.0 : from : 01 apr 1997     **/
/**                                 to   : 20 jun 1997     **/
/**                # Version P0.1 : from : 07 apr 1998     **/
/**                                 to   : 20 jun 1998     **/
/**                # Version P0.2 : from : 11 may 1999     **/
/**                                 to   : 02 feb 2000     **/
/**                # Version P0.3 : from : 16 jun 2005     **/
/**                                 to   : 16 jun 2005     **/
/**                # Version 5.0  : from : 22 jul 2005     **/
/**                                 to   : 03 aug 2007     **/
/**                # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 20 feb 2011     **/
/**                # Version 6.0  : from : 30 aug 2012     **/
/**                                 to   : 26 sep 2012     **/
/**                                                        **/
/************************************************************/

#define DGRAPH_H

#define PTSCOTCH_FOLD_DUP                         /* Activate folding on coarsening */

#ifndef SCOTCH_COMM_PTOP_RAT
#define SCOTCH_COMM_PTOP_RAT        0.25          /* Percentage under which point-to-point is allowed */
#endif /* SCOTCH_COMM_PTOP_RAT */

/*
** The defines.
*/

/* Graph flags. */

#define DGRAPHNONE                  0x0000        /* No options set */

#define DGRAPHFREEPRIV              0x0001        /* Set if private arrays freed on exit */
#define DGRAPHFREECOMM              0x0002        /* MPI communicator has to be freed    */
#define DGRAPHFREETABS              0x0004        /* Set if local arrays freed on exit   */
#define DGRAPHFREEPSID              0x0008        /* Set if procsidtab freed on exit     */
#define DGRAPHFREEEDGEGST           0x0010        /* Set if edgegsttab freed on exit     */
#define DGRAPHHASEDGEGST            0x0020        /* Edge ghost array computed           */
#define DGRAPHVERTGROUP             0x0040        /* All vertex arrays grouped           */
#define DGRAPHEDGEGROUP             0x0080        /* All edge arrays grouped             */
#define DGRAPHFREEALL               (DGRAPHFREEPRIV | DGRAPHFREECOMM | DGRAPHFREETABS | DGRAPHFREEPSID  | DGRAPHFREEEDGEGST)
#define DGRAPHCOMMPTOP              0x0100        /* Use point-to-point collective communication */

#define DGRAPHBITSUSED              0x01FF        /* Significant bits for plain distributed graph routines               */
#define DGRAPHBITSNOTUSED           0x0200        /* Value above which bits not used by plain distributed graph routines */

/* Used in algorithms */

#define COARPERTPRIME               31            /* Prime number */
#define COARHASHPRIME               179           /* Prime number */

/*
** The type and structure definitions.
*/

/* The graph basic types, which must be signed. */

#ifndef GNUMMAX                                   /* If graph.h not included */
typedef INT                 Gnum;                 /* Vertex or edge number   */
#define GNUMMAX                     (INTVALMAX)   /* Maximum Gnum value      */
#define GNUMSTRING                  INTSTRING     /* String to printf a Gnum */
#endif /* GNUMMAX */

#define GNUM_MPI                    COMM_INT      /* MPI type for Gnum is MPI type for INT */

#define GRAPHPART_MPI               COMM_BYTE     /* Raw byte type for graph parts */

/* Tags used for point-to-point communications. */

typedef enum DgraphTag_ {
  TAGPROCVRTTAB = 0,                              /*+ procvrttab message       +*/
  TAGVERTLOCTAB,                                  /*+ vertloctab message       +*/
  TAGVENDLOCTAB,                                  /*+ vendloctab message       +*/
  TAGVELOLOCTAB,                                  /*+ veloloctab message       +*/
  TAGVNUMLOCTAB,                                  /*+ vnumloctab message       +*/
  TAGVLBLLOCTAB,                                  /*+ vlblloctab message       +*/
  TAGEDGELOCTAB,                                  /*+ edgeloctab message       +*/
  TAGEDLOLOCTAB,                                  /*+ edloloctab message       +*/
  TAGDATALOCTAB,                                  /*+ Generic data message     +*/
  TAGOK,                                          /*+ Positive answer          +*/
  TAGBAD,                                         /*+ Negative answer          +*/
  TAGHALO    = 100,                               /*+ Tag class for halo       +*/
  TAGCOARSEN = 200,                               /*+ Tag class for coarsening +*/
  TAGMATCH   = 300,                               /*+ Tag class for matching   +*/
  TAGFOLD    = 400,                               /*+ Tag class for folding    +*/
  TAGBAND    = 500                                /*+ Tag class for band graph +*/
} DgraphTag;

/*+ The graph flag type. +*/

typedef int DgraphFlag;                           /*+ Graph property flags +*/

/*+ The vertex part type, in compressed form. From graph.h +*/

#ifndef GRAPH_H
typedef byte GraphPart;
#endif /* GRAPH_H */

/* The distributed graph structure. */

typedef struct Dgraph_ {
  DgraphFlag                flagval;              /*+ Graph properties                                          +*/
  Gnum                      baseval;              /*+ Base index for edge/vertex arrays                         +*/
  Gnum                      vertglbnbr;           /*+ Global number of vertices                                 +*/
  Gnum                      vertglbmax;           /*+ Maximum number of local vertices over all processes       +*/
  Gnum                      vertgstnbr;           /*+ Number of local + ghost vertices                          +*/
  Gnum                      vertgstnnd;           /*+ vertgstnbr + baseval                                      +*/
  Gnum                      vertlocnbr;           /*+ Local number of vertices                                  +*/
  Gnum                      vertlocnnd;           /*+ Local number of vertices + baseval                        +*/
  Gnum *                    vertloctax;           /*+ Local vertex beginning index array [based]                +*/
  Gnum *                    vendloctax;           /*+ Local vertex end index array [based]                      +*/
  Gnum *                    veloloctax;           /*+ Local vertex load array if present                        +*/
  Gnum                      velolocsum;           /*+ Local sum of all vertex loads                             +*/
  Gnum                      veloglbsum;           /*+ Global sum of all vertex loads                            +*/
  Gnum *                    vnumloctax;           /*+ Arrays of global vertex numbers in original graph         +*/
  Gnum *                    vlblloctax;           /*+ Arrays of vertex labels (when read from file)             +*/
  Gnum                      edgeglbnbr;           /*+ Global number of arcs                                     +*/
  Gnum                      edgeglbmax;           /*+ Maximum number of local edges over all processes          +*/
  Gnum                      edgelocnbr;           /*+ Number of local edges                                     +*/
  Gnum                      edgelocsiz;           /*+ Size of local edge array (= edgelocnbr when compact)      +*/
  Gnum                      edgeglbsmx;           /*+ Maximum size of local edge arrays over all processes      +*/
  Gnum *                    edgegsttax;           /*+ Edge array holding local indices of neighbors [based]     +*/
  Gnum *                    edgeloctax;           /*+ Edge array holding global neighbor numbers [based]        +*/
  Gnum *                    edloloctax;           /*+ Edge load array                                           +*/
  Gnum                      degrglbmax;           /*+ Maximum degree over all processes                         +*/
  MPI_Comm                  proccomm;             /*+ Graph communicator                                        +*/
  int                       prockeyval;           /*+ Communicator key value: folded communicators are distinct +*/
  int                       procglbnbr;           /*+ Number of processes sharing graph data                    +*/
  int                       proclocnum;           /*+ Number of this process                                    +*/
  Gnum *                    procvrttab;           /*+ Global array of vertex number ranges [+1,based]           +*/
  Gnum *                    proccnttab;           /*+ Count array for local number of vertices                  +*/
  Gnum *                    procdsptab;           /*+ Displacement array with respect to proccnttab [+1,based]  +*/
  int                       procngbnbr;           /*+ Number of neighboring processes                           +*/
  int                       procngbmax;           /*+ Maximum number of neighboring processes                   +*/
  int *                     procngbtab;           /*+ Array of neighbor process numbers [sorted]                +*/
  int *                     procrcvtab;           /*+ Number of vertices to receive in ghost vertex sub-arrays  +*/
  int                       procsndnbr;           /*+ Overall size of local send array                          +*/
  int *                     procsndtab;           /*+ Number of vertices to send in ghost vertex sub-arrays     +*/
  int *                     procsidtab;           /*+ Array of indices to build communication vectors (send)    +*/
  int                       procsidnbr;           /*+ Size of the send index array                              +*/
} Dgraph;

/*
** The function prototypes.
*/

int                         dgraphInit          (Dgraph * const, MPI_Comm);
void                        dgraphExit          (Dgraph * const);
void                        dgraphFree          (Dgraph * const);
int                         dgraphLoad          (Dgraph * const, FILE * const, const Gnum, const DgraphFlag);
int                         dgraphSave          (Dgraph * const, FILE * const);
int                         dgraphBuild         (Dgraph * const, const Gnum, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const);
int                         dgraphBuild2        (Dgraph * const, const Gnum, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const, const Gnum, Gnum * const, Gnum * const, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const, const Gnum);
int                         dgraphBuild3        (Dgraph * const, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const, const Gnum, Gnum * const, Gnum * const, const Gnum, const Gnum, Gnum * const, Gnum * const, Gnum * const, const Gnum);
int                         dgraphBuild4        (Dgraph * const);
int                         dgraphBuildHcub     (Dgraph * const, const Gnum, const Gnum, const Gnum);
int                         dgraphBuildGrid3D   (Dgraph * const, const Gnum, const Gnum, const Gnum, const Gnum, const Gnum, const int);
int                         dgraphCheck         (const Dgraph * const);
int                         dgraphView          (const Dgraph * const, FILE * const);
int                         dgraphGhst2         (Dgraph * const, const int);
int                         dgraphBand          (Dgraph * restrict const, const Gnum, Gnum * restrict const, const GraphPart * restrict const, const Gnum, const Gnum, Gnum, Dgraph * restrict const, Gnum * restrict * const, GraphPart * restrict * const, Gnum * const, Gnum * const, Gnum * const);

int                         dgraphFold          (const Dgraph * restrict const, const int, Dgraph * restrict const, const void * restrict const, void ** restrict const, MPI_Datatype);
int                         dgraphFold2         (const Dgraph * restrict const, const int, Dgraph * restrict const, MPI_Comm, const void * restrict const, void ** restrict const, MPI_Datatype);
int                         dgraphFoldDup       (const Dgraph * restrict const, Dgraph * restrict const, void * restrict const, void ** restrict const, MPI_Datatype);
int                         dgraphInduce2       (Dgraph * restrict const, Gnum (*) (Dgraph * restrict const, Dgraph * restrict const, const void * restrict const, Gnum * restrict const), const void * const, const Gnum, Gnum *, Dgraph * restrict const);

int                         dgraphInduceList    (Dgraph * const, const Gnum, const Gnum * const, Dgraph * const);
int                         dgraphInducePart    (Dgraph * const, const GraphPart * restrict const, const Gnum, const GraphPart, Dgraph * const);
#ifdef GRAPH_H
int                         dgraphGather        (const Dgraph * restrict const, Graph * restrict);
int                         dgraphGather2       (const Dgraph * restrict const, Graph * restrict, const int, const Gnum);
int                         dgraphGatherAll     (const Dgraph * restrict const, Graph * restrict);
int                         dgraphGatherAll2    (const Dgraph * restrict const, Graph * restrict, const Gnum, const int);
int                         dgraphScatter       (Dgraph * const, const Graph * const);
#endif /* GRAPH_H */

int                         dgraphHaloSync      (Dgraph * const, void * const, MPI_Datatype);

/*
**  The macro definitions.
*/

#define dgraphGhst(grafptr)         dgraphGhst2 (grafptr, 0) /* Build ghost edge array in addition to local edge array */
#define dgraphGhstReplace(grafptr)  dgraphGhst2 (grafptr, 1) /* Replace local edge array by ghost edge array           */
#define dgraphHasGhst(grafptr)      (((grafptr)->flagval & DGRAPHHASEDGEGST) != 0) /* If graph has a ghost edge array  */
