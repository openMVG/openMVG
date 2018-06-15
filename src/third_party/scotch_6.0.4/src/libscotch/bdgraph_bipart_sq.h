/* Copyright 2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_sq.h                     **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the sequential bipartitioning       **/
/**                routine for distributed graphs.         **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 19 nov 2007     **/
/**                                 to   : 20 nov 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct BdgraphBipartSqParam_ {
  Strat *                   strat;                /*+ Sequential bipartitioning strategy used +*/
} BdgraphBipartSqParam;

/*
**  The function prototypes.
*/

#ifndef BDGRAPH_BIPART_SQ
#define static
#endif

int                         bdgraphBipartSq   (Bdgraph * const, const BdgraphBipartSqParam * const);

static void                 bdgraphBipartSqOpBest (const Gnum * const, Gnum * const, const int * const, const MPI_Datatype * const);

#undef static
