/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : mapping_io.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
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
/**                # Version 5.0  : from : 13 sep 2006     **/
/**                                 to     13 sep 2006     **/
/**                                                        **/
/************************************************************/

#define MAPPING_H

/*
**  The defines.
*/

/*+ Ordering option flags. +*/

#define MAPNONE                     0x0000        /* No options set       */
#define MAPFREEPART                 0x0001        /* Free partition array */

/*
**  The type definitions.
*/

/*+ This structure defines a source label
    to target label mapping element.      +*/

typedef struct MappingLoadMap_ {
  Gnum                      slblnum;              /*+ Source graph vertex label: FIRST +*/
  Gnum                      tlblnum;              /*+ Target architecture vertex label +*/
} MappingLoadMap;

/*+ The source graph sort structure, used to
    sort vertices by increasing label value. +*/

typedef struct MappingLoadPerm {
  Gnum                      vlblnum;              /*+ Vertex label: FIRST +*/
  Gnum                      vertnum;              /*+ Vertex number       +*/
} MappingLoadPerm;
