/* Copyright 2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_ml.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the multi-level graph mapping       **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 10 jul 2010     **/
/**                                 to     10 jul 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     01 jun 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct KgraphMapMlParam_ {
  INT                       coarnbr;              /*+ Minimum number of vertices   +*/
  double                    coarval;              /*+ Coarsening ratio             +*/
  Strat *                   stratlow;             /*+ Strategy at lowest level     +*/
  Strat *                   stratasc;             /*+ Strategy at ascending levels +*/
  int                       typeval;              /*+ Not used                     +*/
} KgraphMapMlParam;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_ML
#define static
#endif

static int                  kgraphMapMlCoarsen  (Kgraph * const, Kgraph * restrict const, GraphCoarsenMulti * restrict * const, const KgraphMapMlParam * const);
static int                  kgraphMapMlUncoarsen (Kgraph * restrict const, Kgraph * const, const GraphCoarsenMulti * const);

int                         kgraphMapMl         (Kgraph * restrict const, const KgraphMapMlParam * const);
static int                  kgraphMapMl2        (Kgraph * restrict const, const KgraphMapMlParam * const);

#undef static
