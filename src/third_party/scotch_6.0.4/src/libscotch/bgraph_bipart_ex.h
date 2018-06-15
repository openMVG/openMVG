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
/**   NAME       : bgraph_bipart_ex.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the exact-balance post-processing   **/
/**                module.                                 **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 25 oct 1994     **/
/**                                 to     26 oct 1994     **/
/**                # Version 3.0  : from : 18 nov 1995     **/
/**                                 to     20 nov 1995     **/
/**                # Version 3.1  : from : 20 nov 1995     **/
/**                                 to     20 nov 1995     **/
/**                # Version 3.2  : from : 15 sep 1996     **/
/**                                 to     13 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 11 dec 2003     **/
/**                                 to     11 dec 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ System-defined constants. +*/

#define BGRAPHBIPARTEXGAINTABLSUBBITS 1

#define BGRAPHBIPARTEXSTATEFREE     ((GainLink *) 0) /*+ Vertex in initial state (TRICK: must be 0) +*/
#define BGRAPHBIPARTEXSTATEUSED     ((GainLink *) 1) /*+ Swapped vertex                             +*/
#define BGRAPHBIPARTEXSTATELINK     ((GainLink *) 2) /*+ Currently in gain table if higher          +*/

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_EX
#define static
#endif

int                         bgraphBipartEx      (Bgraph * restrict const);

#undef static
