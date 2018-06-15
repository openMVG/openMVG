/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_st.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the strategy and method       **/
/**                tables and the generic entry point for  **/
/**                the graph bipartitioning methods.       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 08 oct 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 4.0  : from : 15 jan 2002     **/
/**                                 to     15 jan 2002     **/
/**                # Version 5.0  : from : 27 nov 2006     **/
/**                                 to     13 jan 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The type definitions.
*/

/** Method types. **/

typedef enum BgraphBipartStMethodType_ {
  BGRAPHBIPARTSTMETHBD = 0,                       /*+ Band                   +*/
  BGRAPHBIPARTSTMETHDF,                           /*+ Diffusion              +*/
  BGRAPHBIPARTSTMETHEX,                           /*+ Exactifying            +*/
  BGRAPHBIPARTSTMETHFM,                           /*+ Fiduccia-Mattheyses    +*/
  BGRAPHBIPARTSTMETHGG,                           /*+ Greedy Graph Growing   +*/
  BGRAPHBIPARTSTMETHGP,                           /*+ Gibbs-Poole-Stockmeyer +*/
  BGRAPHBIPARTSTMETHML,                           /*+ Multi-level (strategy) +*/
  BGRAPHBIPARTSTMETHZR,                           /*+ Move all to part zero  +*/
  BGRAPHBIPARTSTMETHNBR                           /*+ Number of methods      +*/
} BgraphBipartStMethodType;

/*
**  The external declarations.
*/

extern StratTab             bgraphbipartststratab;

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_ST
#define static
#endif

int                         bgraphBipartSt      (Bgraph * restrict const, const Strat * restrict const);

#undef static
