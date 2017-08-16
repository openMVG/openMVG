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
/**   NAME       : arch_vcmplt.h                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the variable-sized complete graph   **/
/**                target architecture functions.          **/
/**                                                        **/
/**   DATES      : # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     09 aug 1995     **/
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     20 jul 1996     **/
/**                # Version 3.2  : from : 15 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 08 nov 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     05 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     21 jan 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     26 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_VCMPLT_H_STRUCT
#define ARCH_VCMPLT_H_STRUCT

/*+ The variable-sized complete graph bipartitioning definitions. +*/

typedef struct ArchVcmplt_ {
  int                       padding;              /*+ No data needed +*/
} ArchVcmplt;

typedef struct ArchVcmpltDom_ {
  Anum                      termlvl;              /*+ Terminal depth  +*/
  Anum                      termnum;              /*+ Terminal number +*/
} ArchVcmpltDom;

#endif /* ARCH_VCMPLT_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_VCMPLT_H_PROTO
#define ARCH_VCMPLT_H_PROTO

#ifndef ARCH_VCMPLT
#define static
#endif

#define archVcmpltArchLoad          NULL
#define archVcmpltArchSave          NULL
#define archVcmpltArchFree          NULL
ArchDomNum                  archVcmpltDomNum    (const ArchVcmplt * const, const ArchVcmpltDom * const);
int                         archVcmpltDomTerm   (const ArchVcmplt * const, ArchVcmpltDom * restrict const, const ArchDomNum);
Anum                        archVcmpltDomSize   (const ArchVcmplt * const, const ArchVcmpltDom * const);
#define archVcmpltDomWght           archVcmpltDomSize
Anum                        archVcmpltDomDist   (const ArchVcmplt * const, const ArchVcmpltDom * const, const ArchVcmpltDom * const);
int                         archVcmpltDomFrst   (const ArchVcmplt * const, ArchVcmpltDom * const);
int                         archVcmpltDomLoad   (const ArchVcmplt * const, ArchVcmpltDom * const, FILE * const);
int                         archVcmpltDomSave   (const ArchVcmplt * const, const ArchVcmpltDom * const, FILE * const);
int                         archVcmpltDomBipart (const ArchVcmplt * const, const ArchVcmpltDom * const, ArchVcmpltDom * restrict const, ArchVcmpltDom * restrict const);
int                         archVcmpltDomIncl   (const ArchVcmplt * const, const ArchVcmpltDom * const, const ArchVcmpltDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archVcmpltDomMpiType (const ArchVcmplt * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

#endif /* ARCH_VCMPLT_H_PROTO */
#endif /* ARCH_NOPROTO        */
