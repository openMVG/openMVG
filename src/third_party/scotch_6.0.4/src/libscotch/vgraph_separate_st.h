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
/**   NAME       : vgraph_separate_st.h                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the global separation         **/
/**                strategy and method tables.             **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 oct 1996     **/
/**                                 to     14 nov 1997     **/
/**                # Version 3.3  : from : 31 may 1999     **/
/**                                 to     31 may 1999     **/
/**                # Version 4.0  : from : 06 jan 2002     **/
/**                                 to     19 aug 2004     **/
/**                # Version 5.0  : from : 12 sep 2006     **/
/**                                 to   : 17 feb 2007     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to   : 30 oct 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The type definitions.
*/

/*+ Method types. +*/

typedef enum VgraphSeparateStMethodType_ {
  VGRAPHSEPASTMETHBD = 0,                         /*+ Banding strategy         +*/
  VGRAPHSEPASTMETHES,                             /*+ Edge separation strategy +*/
  VGRAPHSEPASTMETHFM,                             /*+ Fiduccia-Mattheyses      +*/
  VGRAPHSEPASTMETHGG,                             /*+ Greedy Graph Growing     +*/
  VGRAPHSEPASTMETHGP,                             /*+ Gibbs-Poole-Stockmeyer   +*/
  VGRAPHSEPASTMETHML,                             /*+ Multi-level separation   +*/
  VGRAPHSEPASTMETHVW,                             /*+ Partition viewer         +*/
  VGRAPHSEPASTMETHZR,                             /*+ Zero method              +*/
  VGRAPHSEPASTMETHNBR                             /*+ Number of methods        +*/
} VgraphSeparateStMethodType;

/*
**  The external declarations.
*/

extern StratTab             vgraphseparateststratab;

/*
**  The function prototypes.
*/

#ifndef VGRAPH_SEPARATE_ST
#define static
#endif

int                         vgraphSeparateSt    (Vgraph * const, const Strat * const);

#undef static
