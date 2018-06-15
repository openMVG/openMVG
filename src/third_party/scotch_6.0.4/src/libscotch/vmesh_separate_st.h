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
/**   NAME       : vmesh_separate_st.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the global mesh separation    **/
/**                strategy and method tables.             **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 20 sep 2002     **/
/**                                 to     31 oct 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The type definitions.
*/

/*+ Method types. +*/

typedef enum VmeshSeparateStMethodType_ {
  VMESHSEPASTMETHFM = 0,                          /*+ Fiduccia-Mattheyses    +*/
  VMESHSEPASTMETHGG,                              /*+ Greedy Mesh Growing    +*/
#ifdef SCOTCH_DEBUG_VMESH2
  VMESHSEPASTMETHGR,                              /*+ Graph partitioning     +*/
#endif /* SCOTCH_DEBUG_VMESH2 */
  VMESHSEPASTMETHML,                              /*+ Multi-level separation +*/
  VMESHSEPASTMETHZR,                              /*+ Zero method            +*/
  VMESHSEPASTMETHNBR                              /*+ Number of methods      +*/
} VmeshSeparateStMethodType;

/*
**  The external declarations.
*/

extern StratTab             vmeshseparateststratab;

/*
**  The function prototypes.
*/

#ifndef VMESH_SEPARATE_ST
#define static
#endif

int                         vmeshSeparateSt     (Vmesh * restrict const, const Strat * restrict const);

#undef static
