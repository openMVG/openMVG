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
/**   NAME       : mapping.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the declarations for    **/
/**                the mapping handling routines.          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15 dec 1992     **/
/**                                 to     01 apr 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     02 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     04 jul 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     06 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     30 mar 1999     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to     13 nov 2005     **/
/**                # Version 5.1  : from : 25 jun 2008     **/
/**                                 to     04 nov 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     14 sep 2014     **/
/**                                                        **/
/**   NOTES      : # While Anum and Gnum are different     **/
/**                  types, because architectures are      **/
/**                  most often much smaller than          **/
/**                  than graphs and require smaller       **/
/**                  integer ranges, for interface         **/
/**                  consistency reasons as well as for    **/
/**                  variable-sized architecture handling, **/
/**                  they will always amount to the same   **/
/**                  type.                                 **/
/**                                                        **/
/************************************************************/

#define MAPPING_H

/*
**  The defines.
*/

/*+ Prime number for hashing terminal domain numbers. +*/

#define MAPPINGHASHPRIME            17            /*+ Prime number for hashing +*/

/*+ Graph option flags. +*/

#define MAPPINGNONE                 0x0000        /* No options set */

#define MAPPINGFREEPART             0x0001        /* Free partition array */
#define MAPPINGFREEDOMN             0x0002        /* Free domain array    */

/*
**  The type definitions.
*/

/*+ This structure defines an (eventually
    partial) mapping of a source graph to
    a target architecture.                +*/

typedef struct Mapping_ {
  int                       flagval;              /*+ Mapping properties         +*/
  const Graph *             grafptr;              /*+ Graph data                 +*/
  const Arch *              archptr;              /*+ Architecture data          +*/
  Anum *                    parttax;              /*+ Mapping array [vertnbr]    +*/
  ArchDom *                 domntab;              /*+ Array of domains [domnmax] +*/
  Anum                      domnnbr;              /*+ Current number of domains  +*/
  Anum                      domnmax;              /*+ Maximum number of domains  +*/
  ArchDom                   domnorg;              /*+ Initial (sub)domain        +*/
} Mapping;

/*+ The target architecture sort structure, used
    to sort vertices by increasing label value.  +*/

typedef struct MappingHash_ {
  Anum                      termnum;              /*+ Terminal vertex number +*/
  Anum                      domnnum;              /*+ Domain number          +*/
} MappingHash;

/*+ The target architecture sort structure, used
    to sort vertices by increasing label value.  +*/

typedef struct MappingSort_ {
  Anum                      labl;                 /*+ Target architecture vertex label +*/
  Anum                      peri;                 /*+ Inverse permutation              +*/
} MappingSort;

/*
**  The function prototypes.
*/
void                        mapExit             (Mapping * const);

#ifndef MAPPING
#define static
#endif

void                        mapInit             (Mapping * restrict const, const Graph * restrict const, const Arch * restrict const, const ArchDom * restrict const);
void                        mapInit2            (Mapping * restrict const, const Graph * restrict const, const Arch * restrict const, const ArchDom * restrict const, const Anum, const Anum);
void                        mapExit             (Mapping * const);
int                         mapAlloc            (Mapping * const);
void                        mapFree             (Mapping * const);
int                         mapResize           (Mapping * restrict const, const Anum);
int                         mapResize2          (Mapping * restrict const, const Anum);
int                         mapCopy             (Mapping * const, const Mapping * const);
void                        mapFrst             (Mapping * const);
int                         mapBuild            (Mapping * restrict const, const Anum * restrict const);
int                         mapMerge            (Mapping * restrict const, const Anum * restrict const);
void                        mapTerm             (const Mapping * restrict const, Anum * restrict const);
int                         mapLoad             (Mapping * restrict const, const Gnum * restrict const, FILE * restrict const);
int                         mapSave             (const Mapping * restrict const, FILE * restrict const);
int                         mapView             (const Mapping * restrict const, const Graph * restrict const, FILE * const);

#undef static

/*
**  The macro definitions.
*/

#define mapDomain(map,idx)          (&((map)->domntab[(map)->parttax[(idx)]]))
