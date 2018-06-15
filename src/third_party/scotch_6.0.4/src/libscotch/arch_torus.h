/* Copyright 2004,2007,2008,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_torus.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the tori graph target architecture  **/
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
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     05 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     21 jan 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

#ifndef ARCH_TORUS_H_STRUCT
#define ARCH_TORUS_H_STRUCT

/*+ Maximum dimension. +*/

#define ARCHTORUSDIMMAX             5             /* Maximum dimension (at least 3) */

/*+ Data structure equivalence for target architecture array. +*/

#ifdef ARCH
#define ArchTorus2Dom               ArchTorusXDom
#define ArchTorus3Dom               ArchTorusXDom
#endif /* ARCH */

/*
**  The type and structure definitions.
*/

/*+ The xD-torus definitions. +*/

typedef struct ArchTorusX_ {
  Anum                      dimmax;               /*+ Number of torus dimensions +*/
  Anum                      c[ARCHTORUSDIMMAX];   /*+ Mesh dimensions            +*/
} ArchTorusX;

typedef struct ArchTorusXDom_ {
  Anum                      c[ARCHTORUSDIMMAX][2]; /*+ Inclusive X and Y coordinates +*/
} ArchTorusXDom;

#endif /* ARCH_TORUS_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_TORUS_H_PROTO
#define ARCH_TORUS_H_PROTO

#ifndef ARCH_TORUS
#define static
#endif

int                         archTorus2ArchLoad  (ArchTorusX * restrict const, FILE * restrict const);
int                         archTorus2ArchSave  (const ArchTorusX * const, FILE * restrict const);
#define archTorus2ArchFree          NULL
ArchDomNum                  archTorus2DomNum    (const ArchTorusX * const, const ArchTorusXDom * const);
int                         archTorus2DomTerm   (const ArchTorusX * const, ArchTorusXDom * restrict const, const ArchDomNum);
Anum                        archTorus2DomSize   (const ArchTorusX * const, const ArchTorusXDom * const);
#define archTorus2DomWght           archTorus2DomSize
Anum                        archTorus2DomDist   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
#define archTorus2DomFrst           archTorusXDomFrst
#define archTorus2DomLoad           archTorusXDomLoad
#define archTorus2DomSave           archTorusXDomSave
int                         archTorus2DomBipart (const ArchTorusX * const, const ArchTorusXDom * const, ArchTorusXDom * restrict const, ArchTorusXDom * restrict const);
int                         archTorus2DomBipartO (const ArchTorusX * const, const ArchTorusXDom * const, ArchTorusXDom * restrict const, ArchTorusXDom * restrict const);
int                         archTorus2DomBipartU (const ArchTorusX * const, const ArchTorusXDom * const, ArchTorusXDom * restrict const, ArchTorusXDom * restrict const);
int                         archTorus2DomIncl   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
#ifdef SCOTCH_PTSCOTCH
#define archTorus2DomMpiType        archTorusXDomMpiType
#endif /* SCOTCH_PTSCOTCH */

int                         archTorus3ArchLoad  (ArchTorusX * restrict const, FILE * restrict const);
int                         archTorus3ArchSave  (const ArchTorusX * const, FILE * restrict const);
#define archTorus3ArchFree          NULL
ArchDomNum                  archTorus3DomNum    (const ArchTorusX * const, const ArchTorusXDom * const);
int                         archTorus3DomTerm   (const ArchTorusX * const, ArchTorusXDom * restrict const, const ArchDomNum);
Anum                        archTorus3DomSize   (const ArchTorusX * const, const ArchTorusXDom * const);
#define archTorus3DomWght           archTorus3DomSize
Anum                        archTorus3DomDist   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
#define archTorus3DomFrst           archTorusXDomFrst
#define archTorus3DomLoad           archTorusXDomLoad
#define archTorus3DomSave           archTorusXDomSave
int                         archTorus3DomBipart (const ArchTorusX * const, const ArchTorusXDom * const, ArchTorusXDom * restrict const, ArchTorusXDom * restrict const);
int                         archTorus3DomIncl   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
#ifdef SCOTCH_PTSCOTCH
#define archTorus3DomMpiType        archTorusXDomMpiType
#endif /* SCOTCH_PTSCOTCH */

int                         archTorusXArchLoad  (ArchTorusX * restrict const, FILE * restrict const);
int                         archTorusXArchSave  (const ArchTorusX * const, FILE * restrict const);
#define archTorusXArchFree          NULL
ArchDomNum                  archTorusXDomNum    (const ArchTorusX * const, const ArchTorusXDom * const);
int                         archTorusXDomTerm   (const ArchTorusX * const, ArchTorusXDom * restrict const, const ArchDomNum);
Anum                        archTorusXDomSize   (const ArchTorusX * const, const ArchTorusXDom * const);
#define archTorusXDomWght           archTorusXDomSize
Anum                        archTorusXDomDist   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
int                         archTorusXDomFrst   (const ArchTorusX * const, ArchTorusXDom * const);
int                         archTorusXDomLoad   (const ArchTorusX * const, ArchTorusXDom * const, FILE * restrict const);
int                         archTorusXDomSave   (const ArchTorusX * const, const ArchTorusXDom * const, FILE * restrict const);
int                         archTorusXDomBipart (const ArchTorusX * const, const ArchTorusXDom * const, ArchTorusXDom * restrict const, ArchTorusXDom * restrict const);
int                         archTorusXDomIncl   (const ArchTorusX * const, const ArchTorusXDom * const, const ArchTorusXDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archTorusXDomMpiType (const ArchTorusX * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_TORUS_H_PROTO */
#endif /* ARCH_NOPROTO       */

