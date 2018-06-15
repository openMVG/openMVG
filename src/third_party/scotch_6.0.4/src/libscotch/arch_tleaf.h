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
/**   NAME       : arch_tleaf.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the tree-leaf pseudo-graph target   **/
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
/**                                 to     16 aug 1995     **/
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 10 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 10 dec 2003     **/
/**                                 to     10 dec 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     24 jun 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     03 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

#ifndef ARCH_TLEAF_H_STRUCT
#define ARCH_TLEAF_H_STRUCT

/** The Tree-Leaf graph definitions. **/

typedef struct ArchTleaf_ {
  Anum                      termnbr;              /*+ Number of terminal domains in architecture   +*/
  Anum                      levlnbr;              /*+ Number of levels                             +*/
  Anum *                    sizetab;              /*+ Array of cluster sizes, per descending level +*/
  Anum *                    linktab;              /*+ Value of extra-cluster link costs            +*/
  Anum                      permnbr;              /*+ Number of label permutation indices          +*/
  Anum *                    permtab;              /*+ Label permutation array                      +*/
  Anum *                    peritab;              /*+ Invertse label permutation array             +*/
} ArchTleaf;

typedef struct ArchTleafDom_ {
  Anum                      levlnum;              /*+ Current block level         +*/
  Anum                      indxmin;              /*+ Minimum index in level      +*/
  Anum                      indxnbr;              /*+ Number of indices in domain +*/
} ArchTleafDom;

#endif /* ARCH_TLEAF_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_TLEAF_H_PROTO
#define ARCH_TLEAF_H_PROTO

#ifndef ARCH_TLEAF
#define static
#endif

int                         archTleafArchLoad   (ArchTleaf * restrict const, FILE * restrict const);
int                         archTleafArchFree   (ArchTleaf * restrict const);
int                         archTleafArchSave   (const ArchTleaf * const, FILE * restrict const);
ArchDomNum                  archTleafDomNum     (const ArchTleaf * const, const ArchTleafDom * const);
int                         archTleafDomTerm    (const ArchTleaf * const, ArchTleafDom * restrict const, const ArchDomNum);
Anum                        archTleafDomSize    (const ArchTleaf * const, const ArchTleafDom * const);
#define archTleafDomWght            archTleafDomSize
Anum                        archTleafDomDist    (const ArchTleaf * const, const ArchTleafDom * const, const ArchTleafDom * const);
int                         archTleafDomFrst    (const ArchTleaf * const, ArchTleafDom * restrict const);
int                         archTleafDomLoad    (const ArchTleaf * const, ArchTleafDom * restrict const, FILE * restrict const);
int                         archTleafDomSave    (const ArchTleaf * const, const ArchTleafDom * const, FILE * restrict const);
int                         archTleafDomBipart  (const ArchTleaf * const, const ArchTleafDom * const, ArchTleafDom * restrict const, ArchTleafDom * restrict const);
int                         archTleafDomIncl    (const ArchTleaf * const, const ArchTleafDom * const, const ArchTleafDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archTleafDomMpiType (const ArchTleaf * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

int                         archLtleafArchLoad  (ArchTleaf * restrict const, FILE * restrict const);
int                         archLtleafArchSave  (const ArchTleaf * const, FILE * restrict const);
ArchDomNum                  archLtleafDomNum    (const ArchTleaf * const, const ArchTleafDom * const);
int                         archLtleafDomTerm   (const ArchTleaf * const, ArchTleafDom * restrict const, const ArchDomNum);
#define archLtleafDomWght           archLtleafDomSize

#undef static

/*
**  The macro definitions.
*/

#define ArchLtleaf                  ArchTleaf
#define ArchLtleafDom               ArchTleafDom

#define archLtleafArchFree          archTleafArchFree
#define archLtleafDomSize           archTleafDomSize
#define archLtleafDomDist           archTleafDomDist
#define archLtleafDomFrst           archTleafDomFrst
#define archLtleafDomLoad           archTleafDomLoad
#define archLtleafDomSave           archTleafDomSave
#define archLtleafDomBipart         archTleafDomBipart
#define archLtleafDomIncl           archTleafDomIncl
#define archLtleafDomMpiType        archTleafDomMpiType

#endif /* ARCH_TLEAF_H_PROTO */
#endif /* ARCH_NOPROTO       */
