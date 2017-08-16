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
/**   NAME       : arch_deco.h                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the decomposition-defined target    **/
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
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     20 jul 1996     **/
/**                # Version 3.2  : from : 11 sep 1996     **/
/**                                 to     28 sep 1998     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     14 jun 2004     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     27 sep 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

#ifndef ARCH_DECO_H_STRUCT
#define ARCH_DECO_H_STRUCT

/*+ Decomposition architecture flags. +*/

#define ARCHDECONONE                0x0000        /*+ No options set +*/

#define ARCHDECOFREE                0x0001        /*+ Free arrays    +*/

/*
**  The type and structure definitions.
*/

/*+ The decomposition-described terminal vertex definition. +*/

typedef struct ArchDecoTermVert_ {
  ArchDomNum                labl;                 /*+ Number for terminals, or ARCHDOMNOTTERM +*/
  Anum                      wght;                 /*+ Weight of the domain (processor load)   +*/
  Anum                      num;                  /*+ Number of the terminal                  +*/
} ArchDecoTermVert;

/*+ The decomposition-described architecture definitions. +*/

typedef struct ArchDecoVert_ {
  ArchDomNum                labl;                 /*+ Smallest number of included terminal  +*/
  Anum                      size;                 /*+ Number of processors in the domain    +*/
  Anum                      wght;                 /*+ Weight of the domain (processor load) +*/
} ArchDecoVert;

typedef struct ArchDeco_ {
  int                       flagval;              /*+ Flag value                 +*/
  Anum                      domtermnbr;           /*+ Number of terminal domains +*/
  Anum                      domvertnbr;           /*+ Number of domains          +*/
  ArchDecoVert *            domverttab;           /*+ Table of domain "vertices" +*/
  Anum *                    domdisttab;           /*+ Table of domain distances  +*/
} ArchDeco;

typedef struct ArchDecoDom_ {
  Anum                      num;                  /*+ Domain number in the decomposition +*/
} ArchDecoDom;

#endif /* ARCH_DECO_H_STRUCT */

/*
**  The function prototypes.
*/

#ifndef ARCH_NOPROTO
#ifndef ARCH_DECO_H_PROTO
#define ARCH_DECO_H_PROTO

#ifndef ARCH_DECO
#define static
#endif

int                         archDecoArchBuild   (ArchDeco * const, const Anum, const Anum, const ArchDecoTermVert * const, const Anum  * const);
int                         archDecoArchLoad    (ArchDeco * const, FILE * restrict const);
int                         archDecoArchSave    (const ArchDeco * const, FILE * restrict const);
int                         archDecoArchFree    (ArchDeco * const);
Anum                        archDecoArchSize    (ArchDeco * const, const Anum);
Anum                        archDecoArchDist    (ArchDeco * const, const Anum, const Anum);
Anum                        archDecoArchDistE   (ArchDeco * const, const Anum, const Anum);
ArchDomNum                  archDecoDomNum      (const ArchDeco * const, const ArchDecoDom * const);
int                         archDecoDomTerm     (const ArchDeco * const, ArchDecoDom * restrict const, const ArchDomNum);
Anum                        archDecoDomSize     (const ArchDeco * const, const ArchDecoDom * const);
Anum                        archDecoDomWght     (const ArchDeco * const, const ArchDecoDom * const);
Anum                        archDecoDomDist     (const ArchDeco * const, const ArchDecoDom * const, const ArchDecoDom * const);
int                         archDecoDomFrst     (const ArchDeco * const, ArchDecoDom * restrict const);
int                         archDecoDomLoad     (const ArchDeco * const, ArchDecoDom * restrict const, FILE * restrict const);
int                         archDecoDomSave     (const ArchDeco * const, const ArchDecoDom * const, FILE * restrict const);
int                         archDecoDomBipart   (const ArchDeco * const, const ArchDecoDom * const, ArchDecoDom * restrict const, ArchDecoDom * restrict const);
int                         archDecoDomIncl     (const ArchDeco * const, const ArchDecoDom * const, const ArchDecoDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archDecoDomMpiType  (const ArchDeco * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

/*
**  The macro definitions.
*/

#define archDecoArchSize(d,i)       ((d)->domverttab[(i) - 1].size)
#define archDecoArchDist(d,i,j)     ((d)->domdisttab[((i) >= (j)) ? (((i) - 1) * ((i) - 2)) / 2 + (j) - 1 \
                                                                  : (((j) - 1) * ((j) - 2)) / 2 + (i) - 1])
#define archDecoArchDistE(d,i,j)    (((i) == (j)) ? 0 : archDecoArchDist ((d), (i), (j)))

#endif /* ARCH_DECO_H_PROTO */
#endif /* ARCH_NOPROTO      */
