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
/**   NAME       : arch_cmplt.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the complete graph target           **/
/**                architecture functions.                 **/
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
/**                                 to     24 jul 1995     **/
/**                # Version 3.1  : from : 11 jun 1996     **/
/**                                 to     11 jun 1996     **/
/**                # Version 3.2  : from : 20 sep 1996     **/
/**                                 to     13 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     09 jan 2004     **/
/**                # Version 5.1  : from : 19 jan 2008     **/
/**                                 to     19 jan 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_CMPLT_H_STRUCT
#define ARCH_CMPLT_H_STRUCT

/*+ The complete graph definitions. +*/

typedef struct ArchCmplt_ {
  Anum                      numnbr;               /*+ Number of vertices +*/
} ArchCmplt;

typedef struct ArchCmpltDom_ {
  Anum                      nummin;               /*+ Minimum vertex number +*/
  Anum                      numnbr;               /*+ Number of vertices    +*/
} ArchCmpltDom;

#endif /* ARCH_CMPLT_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_CMPLT_H_PROTO
#define ARCH_CMPLT_H_PROTO

#ifndef ARCH_CMPLT
#define static
#endif

int                         archCmpltArchLoad   (ArchCmplt * restrict const, FILE * restrict const);
int                         archCmpltArchSave   (const ArchCmplt * const, FILE * restrict const);
#define archCmpltArchFree           NULL
ArchDomNum                  archCmpltDomNum     (const ArchCmplt * const, const ArchCmpltDom * const);
int                         archCmpltDomTerm    (const ArchCmplt * const, ArchCmpltDom * restrict const, const ArchDomNum);
Anum                        archCmpltDomSize    (const ArchCmplt * const, const ArchCmpltDom * const);
#define archCmpltDomWght            archCmpltDomSize
Anum                        archCmpltDomDist    (const ArchCmplt * const, const ArchCmpltDom * const, const ArchCmpltDom * const);
int                         archCmpltDomFrst    (const ArchCmplt * const, ArchCmpltDom * const);
int                         archCmpltDomLoad    (const ArchCmplt * const, ArchCmpltDom * const, FILE * const);
int                         archCmpltDomSave    (const ArchCmplt * const, const ArchCmpltDom * const, FILE * const);
int                         archCmpltDomBipart  (const ArchCmplt * const, const ArchCmpltDom * const, ArchCmpltDom * restrict const, ArchCmpltDom * restrict const);
int                         archCmpltDomIncl    (const ArchCmplt * const, const ArchCmpltDom * const, const ArchCmpltDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archCmpltDomMpiType (const ArchCmplt * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_CMPLT_H_PROTO */
#endif /* ARCH_NOPROTO       */
