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
/**   NAME       : arch_vhcub.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the variable-sized hypercube        **/
/**                target architecture functions.          **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 08 nov 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 04 nov 2003     **/
/**                                 to     04 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     21 jan 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_VHCUB_H_STRUCT
#define ARCH_VHCUB_H_STRUCT

/*+ The variable-sized hypercube bipartitioning definitions. +*/

typedef struct ArchVhcub_ {
  int                       padding;              /*+ No data needed +*/
} ArchVhcub;

typedef struct ArchVhcubDom_ {
  Anum                      termlvl;              /*+ Terminal depth  +*/
  Anum                      termnum;              /*+ Terminal number +*/
} ArchVhcubDom;

#endif /* ARCH_VHCUB_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_VHCUB_H_PROTO
#define ARCH_VHCUB_H_PROTO

#ifndef ARCH_VHCUB
#define static
#endif

#define archVhcubArchLoad           NULL
#define archVhcubArchSave           NULL
#define archVhcubArchFree           NULL
ArchDomNum                  archVhcubDomNum     (const ArchVhcub * const, const ArchVhcubDom * const);
int                         archVhcubDomTerm    (const ArchVhcub * const, ArchVhcubDom * restrict const, const ArchDomNum);
Anum                        archVhcubDomSize    (const ArchVhcub * const, const ArchVhcubDom * const);
#define archVhcubDomWght            archVhcubDomSize
Anum                        archVhcubDomDist    (const ArchVhcub * const, const ArchVhcubDom * const, const ArchVhcubDom * const);
int                         archVhcubDomFrst    (const ArchVhcub * const, ArchVhcubDom * const);
int                         archVhcubDomLoad    (const ArchVhcub * const, ArchVhcubDom * const, FILE * const);
int                         archVhcubDomSave    (const ArchVhcub * const, const ArchVhcubDom * const, FILE * const);
int                         archVhcubDomBipart  (const ArchVhcub * const, const ArchVhcubDom * const, ArchVhcubDom * restrict const, ArchVhcubDom * restrict const);
int                         archVhcubDomIncl    (const ArchVhcub * const, const ArchVhcubDom * const, const ArchVhcubDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archVhcubDomMpiType (const ArchVhcub * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_VHCUB_H_PROTO */
#endif /* ARCH_NOPROTO       */
