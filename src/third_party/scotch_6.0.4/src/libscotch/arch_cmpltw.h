/* Copyright 2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_cmpltw.h                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the weighted complete graph target  **/
/**                architecture functions.                 **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 dec 2007     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     23 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_CMPLTW_H_STRUCT
#define ARCH_CMPLTW_H_STRUCT

/*+ The weighted target vertex. Since Anum's
    are INT's, they can be sorted, by means
    of the intSort2asc1 routine.             +*/

typedef struct ArchCmpltwLoad_ {
  Anum                      veloval;              /*+ Vertex load  +*/
  Anum                      vertnum;              /*+ Vertex index +*/
} ArchCmpltwLoad;

/*+ The weighted complete graph definitions. +*/

typedef struct ArchCmpltw_ {
  Anum                      vertnbr;              /*+ Number of vertices +*/
  ArchCmpltwLoad *          velotab;              /*+ Vertex index array +*/
  Anum                      velosum;              /*+ Sum of all weights +*/
} ArchCmpltw;

/*+ The weighted domain structure. +*/

typedef struct ArchCmpltwDom_ {
  Anum                      vertmin;              /*+ Minimum vertex number +*/
  Anum                      vertnbr;              /*+ Number of vertices    +*/
  Anum                      veloval;              /*+ Weight of subdomain   +*/
} ArchCmpltwDom;

#endif /* ARCH_CMPLTW_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_CMPLTW_H_PROTO
#define ARCH_CMPLTW_H_PROTO

#ifndef ARCH_CMPLTW
#define static
#endif

int                         archCmpltwArchBuild (ArchCmpltw * restrict const archptr, const Anum, const Anum * restrict const);
int                         archCmpltwArchLoad  (ArchCmpltw * restrict const, FILE * restrict const);
int                         archCmpltwArchSave  (const ArchCmpltw * const, FILE * restrict const);
int                         archCmpltwArchFree  (ArchCmpltw * restrict const);
ArchDomNum                  archCmpltwDomNum    (const ArchCmpltw * const, const ArchCmpltwDom * const);
int                         archCmpltwDomTerm   (const ArchCmpltw * const, ArchCmpltwDom * restrict const, const ArchDomNum);
Anum                        archCmpltwDomSize   (const ArchCmpltw * const, const ArchCmpltwDom * const);
Anum                        archCmpltwDomWght   (const ArchCmpltw * const, const ArchCmpltwDom * const);
Anum                        archCmpltwDomDist   (const ArchCmpltw * const, const ArchCmpltwDom * const, const ArchCmpltwDom * const);
int                         archCmpltwDomFrst   (const ArchCmpltw * const, ArchCmpltwDom * const);
int                         archCmpltwDomLoad   (const ArchCmpltw * const, ArchCmpltwDom * const, FILE * const);
int                         archCmpltwDomSave   (const ArchCmpltw * const, const ArchCmpltwDom * const, FILE * const);
int                         archCmpltwDomBipart (const ArchCmpltw * const, const ArchCmpltwDom * const, ArchCmpltwDom * restrict const, ArchCmpltwDom * restrict const);
int                         archCmpltwDomIncl   (const ArchCmpltw * const, const ArchCmpltwDom * const, const ArchCmpltwDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archCmpltwDomMpiType (const ArchCmpltw * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_CMPLTW_H_PROTO */
#endif /* ARCH_NOPROTO        */
