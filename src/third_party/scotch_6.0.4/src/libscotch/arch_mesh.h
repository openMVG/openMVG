/* Copyright 2004,2007,2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_mesh.h                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the mesh graph target architecture  **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to   : 24 mar 1993     **/
/**                # Version 1.2  : from : 04 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 1.3  : from : 20 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 12 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to   : 30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     17 aug 1995     **/
/**                # Version 3.1  : from : 22 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 16 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     09 jan 2004     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     21 jan 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_MESH_H_STRUCT
#define ARCH_MESH_H_STRUCT

/*+ The 2D-mesh definitions. +*/

typedef struct ArchMesh2_ {
  Anum                      c[2];                 /*+ Mesh dimensions +*/
} ArchMesh2;

typedef struct ArchMesh2Dom_ {
  Anum                      c[2][2];              /*+ Inclusive X and Y coordinates +*/
} ArchMesh2Dom;

/*+ The 3D-mesh definitions. +*/

typedef struct ArchMesh3_ {
  Anum                      c[3];                 /*+ Mesh dimensions +*/
} ArchMesh3;

typedef struct ArchMesh3Dom_ {
  Anum                      c[3][2];              /*+ Inclusive X, Y, and Z coordinates +*/
} ArchMesh3Dom;

#endif /* ARCH_MESH_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_MESH_H_PROTO
#define ARCH_MESH_H_PROTO

#ifndef ARCH_MESH
#define static
#endif

int                         archMesh2ArchLoad   (ArchMesh2 * restrict const, FILE * restrict const);
int                         archMesh2ArchSave   (const ArchMesh2 * const, FILE * restrict const);
#define archMesh2ArchFree           NULL
ArchDomNum                  archMesh2DomNum     (const ArchMesh2 * const, const ArchMesh2Dom * const);
int                         archMesh2DomTerm    (const ArchMesh2 * const, ArchMesh2Dom * restrict const, const ArchDomNum);
Anum                        archMesh2DomSize    (const ArchMesh2 * const, const ArchMesh2Dom * const);
#define archMesh2DomWght            archMesh2DomSize
Anum                        archMesh2DomDist    (const ArchMesh2 * const, const ArchMesh2Dom * const, const ArchMesh2Dom * const);
int                         archMesh2DomFrst    (const ArchMesh2 * const, ArchMesh2Dom * const);
int                         archMesh2DomLoad    (const ArchMesh2 * const, ArchMesh2Dom * const, FILE * restrict const);
int                         archMesh2DomSave    (const ArchMesh2 * const, const ArchMesh2Dom * const, FILE * restrict const);
int                         archMesh2DomBipart  (const ArchMesh2 * const, const ArchMesh2Dom * const, ArchMesh2Dom * restrict const, ArchMesh2Dom * restrict const);
int                         archMesh2DomBipartO (const ArchMesh2 * const, const ArchMesh2Dom * const, ArchMesh2Dom * restrict const, ArchMesh2Dom * restrict const);
int                         archMesh2DomBipartU (const ArchMesh2 * const, const ArchMesh2Dom * const, ArchMesh2Dom * restrict const, ArchMesh2Dom * restrict const);
int                         archMesh2DomIncl    (const ArchMesh2 * const, const ArchMesh2Dom * const, const ArchMesh2Dom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archMesh2DomMpiType (const ArchMesh2 * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

int                         archMesh3ArchLoad   (ArchMesh3 * restrict const, FILE * restrict const);
int                         archMesh3ArchSave   (const ArchMesh3 * const, FILE * restrict const);
#define archMesh3ArchFree           NULL
ArchDomNum                  archMesh3DomNum     (const ArchMesh3 * const, const ArchMesh3Dom * const);
int                         archMesh3DomTerm    (const ArchMesh3 * const, ArchMesh3Dom * restrict const, const ArchDomNum);
Anum                        archMesh3DomSize    (const ArchMesh3 * const, const ArchMesh3Dom * const);
#define archMesh3DomWght            archMesh3DomSize
Anum                        archMesh3DomDist    (const ArchMesh3 * const, const ArchMesh3Dom * const, const ArchMesh3Dom * const);
int                         archMesh3DomFrst    (const ArchMesh3 * const, ArchMesh3Dom * const);
int                         archMesh3DomLoad    (const ArchMesh3 * const, ArchMesh3Dom * const, FILE * restrict const);
int                         archMesh3DomSave    (const ArchMesh3 * const, const ArchMesh3Dom * const, FILE * restrict const);
int                         archMesh3DomBipart  (const ArchMesh3 * const, const ArchMesh3Dom * const, ArchMesh3Dom * restrict const, ArchMesh3Dom * restrict const);
int                         archMesh3DomIncl    (const ArchMesh3 * const, const ArchMesh3Dom * const, const ArchMesh3Dom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archMesh3DomMpiType (const ArchMesh3 * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_MESH_H_PROTO */
#endif /* ARCH_NOPROTO      */
