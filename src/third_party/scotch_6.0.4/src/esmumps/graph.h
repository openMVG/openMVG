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
/**   NAME       : graph.h                                 **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the graph structure.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 jul 1998     **/
/**                                 to     24 jan 2004     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     23 apr 2004     **/
/**                                                        **/
/************************************************************/

#define GRAPH_H

/*
**  The defines.
*/

#define Graph                       SCOTCH_Graph

#define graphInit                   SCOTCH_graphInit
#define graphExit                   SCOTCH_graphExit
#define graphLoad                   SCOTCH_graphLoad
#define graphSave                   SCOTCH_graphSave
#define graphBase                   SCOTCH_graphBase
#define graphData                   SCOTCH_graphData
#define graphCheck                  SCOTCH_graphCheck

/*
**  The function prototypes.
*/

#ifndef GRAPH
#define static
#endif

int                         graphBuild          (Graph * const grafptr, const INT baseval, const INT vertnbr, const INT edgenbr, void * const ngbdptr, INT nghbfrstfunc (void * const, const INT), INT nghbnextfunc (void * const));
int                         graphBuildGraph     (Graph * const grafptr, const INT baseval, const INT vertnbr, const INT edgenbr, INT * restrict verttab, INT * restrict velotab, INT * restrict edgetab);
int                         graphBuildGraph2    (Graph * const grafptr, const INT baseval, const INT vertnbr, const INT edgenbr, INT * restrict verttab, INT * restrict vendtab, INT * restrict velotab, INT * restrict vlbltab, INT * restrict edgetab, INT * restrict edlotab);

#undef static
