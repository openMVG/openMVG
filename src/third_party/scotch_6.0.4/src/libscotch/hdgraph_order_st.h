/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_order_st.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data           **/
/**                declarations for the main distributed   **/
/**                graph ordering routine.                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 14 apr 2006     **/
/**                                 to   : 14 apr 2006     **/
/**                # Version 5.1  : from : 11 nov 2008     **/
/**                                 to   : 11 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The type definitions.
*/

/*+ Method types. +*/

typedef enum HdgraphOrderStMethodType_ {
  HDGRAPHORDERSTMETHND = 0,                       /*+ Nested Dissection +*/
  HDGRAPHORDERSTMETHSI,                           /*+ Simple            +*/
  HDGRAPHORDERSTMETHSQ,                           /*+ Sequential method +*/
  HDGRAPHORDERSTMETHNBR                           /*+ Number of methods +*/
} HdgraphOrderStMethodType;

/*
**  The external declarations.
*/

extern StratTab             hdgraphorderststratab;

/*
**  The function prototypes.
*/

#ifndef HDGRAPH_ORDER_ST
#define static
#endif

int                         hdgraphOrderSt      (const Hdgraph * restrict const, DorderCblk * restrict const, const Strat * const);

#undef static
