/* Copyright 2012,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_match.h                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the centralized source graph        **/
/**                matching routines.                      **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 02 oct 2012     **/
/**                                 to     25 feb 2015     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

#if (defined SCOTCH_PTHREAD) && (defined GRAPHCOARSENTHREAD)
#define GRAPHMATCHTHREAD
#endif /* (defined SCOTCH_PTHREAD) && (defined GRAPHCOARSENTHREAD) */

/** Prime number for cache-friendly perturbations. **/

#define GRAPHMATCHSCANPERTPRIME     179           /* Prime number */

/** Function block building macro. **/

#define GRAPHMATCHFUNCBLOCK(t)      graphMatch##t##NfNvNe, \
                                    graphMatch##t##NfNvEl, \
                                    graphMatch##t##NfVlNe, \
                                    graphMatch##t##NfVlEl, \
                                    graphMatch##t##FxNvNe, \
                                    graphMatch##t##FxNvEl, \
                                    graphMatch##t##FxVlNe, \
                                    graphMatch##t##FxVlEl

#define GRAPHMATCHFUNCDECL(t)       static void graphMatch##t##NfNvNe (GraphCoarsenThread *); \
                                    static void graphMatch##t##NfNvEl (GraphCoarsenThread *); \
                                    static void graphMatch##t##NfVlNe (GraphCoarsenThread *); \
                                    static void graphMatch##t##NfVlEl (GraphCoarsenThread *); \
                                    static void graphMatch##t##FxNvNe (GraphCoarsenThread *); \
                                    static void graphMatch##t##FxNvEl (GraphCoarsenThread *); \
                                    static void graphMatch##t##FxVlNe (GraphCoarsenThread *); \
                                    static void graphMatch##t##FxVlEl (GraphCoarsenThread *);

/*
**  The function prototypes.
*/

void                        graphMatchNone      (GraphCoarsenData *);
int                         graphMatchInit      (GraphCoarsenData *);
void                        graphMatch          (GraphCoarsenThread * restrict const);

#ifndef GRAPH_MATCH
#define static
#endif

GRAPHMATCHFUNCDECL (Seq);

GRAPHMATCHFUNCDECL (ThrBeg);

GRAPHMATCHFUNCDECL (ThrMid);

GRAPHMATCHFUNCDECL (ThrEnd);

#undef static
