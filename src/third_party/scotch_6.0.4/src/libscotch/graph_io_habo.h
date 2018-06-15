/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_io_habo.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the Harwell-Boeing matrix     **/
/**                format I/O module.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 06 nov 1997     **/
/**                                 to     06 nov 1997     **/
/**                # Version 3.3  : from : 13 dec 1998     **/
/**                                 to     15 dec 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     19 jan 2004     **/
/**                # Version 5.0  : from : 06 jun 2007     **/
/**                                 to     06 jun 2007     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Prime number for hashing vertex numbers. +*/

#define GRAPHGEOMHABOHASHPRIME      7             /*+ Prime number +*/

/*
**  The type and structure definitions.
*/

/*+ This structure holds neighbor vertex hashing data. +*/

typedef struct GraphGeomHaboHash_ {
  Gnum                      vertnum;              /*+ Origin vertex (i.e. pass) number +*/
  Gnum                      vertend;              /*+ Adjacent end vertex number       +*/
} GraphGeomHaboHash;

/*+ This structure holds line formats for reading input data +*/

typedef struct GraphGeomHaboLine_ {
  int                       strtnbr;              /*+ Number of starting blank characters +*/
  int                       datanbr;              /*+ Number of integers par line         +*/
  int                       datalen;              /*+ Number of characters per integer    +*/
} GraphGeomHaboLine;

/*
**  The function prototypes.
*/

#ifndef GRAPH_IO_HABO
#define static
#endif

static int                  graphGeomLoadHaboFormat (GraphGeomHaboLine * restrict const, const char * const);

#undef static
