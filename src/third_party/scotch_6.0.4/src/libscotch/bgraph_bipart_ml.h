/* Copyright 2004,2007,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_ml.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Luca SCARANO (v3.1)                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the multi-level graph bipartitio-   **/
/**                ning routines.                          **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 24 oct 1995     **/
/**                                 to     03 jul 1996     **/
/**                # Version 3.2  : from : 20 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 12 dec 2003     **/
/**                                 to     20 mar 2005     **/
/**                # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     13 jul 2010     **/
/**                # Version 6.0  : from : 16 apr 2011     **/
/**                                 to     03 sep 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct BgraphBipartMlParam_ {
  INT                       coarnbr;              /*+ Minimum number of vertices   +*/
  double                    coarrat;              /*+ Coarsening ratio             +*/
  Strat *                   stratlow;             /*+ Strategy at lowest level     +*/
  Strat *                   stratasc;             /*+ Strategy at ascending levels +*/
} BgraphBipartMlParam;

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_ML
#define static
#endif

static int                  bgraphBipartMlCoarsen (const Bgraph * const, Bgraph * restrict const, GraphCoarsenMulti * restrict * const, const BgraphBipartMlParam * const);
static int                  bgraphBipartMlUncoarsen (Bgraph * restrict const, const Bgraph * restrict const, const GraphCoarsenMulti * const);

int                         bgraphBipartMl      (Bgraph * restrict const, const BgraphBipartMlParam * const);
static int                  bgraphBipartMl2     (Bgraph * restrict const, const BgraphBipartMlParam * const);

#undef static
