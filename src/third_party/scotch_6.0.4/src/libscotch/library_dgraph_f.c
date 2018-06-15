/* Copyright 2007-2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_f.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the distributed source graph handling   **/
/**                routines of the libSCOTCH library.      **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 04 sep 2006     **/
/**                                 to     05 aug 2007     **/
/**                # Version 5.1  : from : 27 jul 2008     **/
/**                                 to     15 apr 2010     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "ptscotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the graph handling routines.   */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHINIT, scotchfdgraphinit, (         \
SCOTCH_Dgraph * const       grafptr,            \
const MPI_Fint * const      commptr,            \
int * const                 revaptr),           \
(grafptr, commptr, revaptr))
{
  MPI_Comm            commdat;

  commdat = MPI_Comm_f2c (*commptr);
  *revaptr = SCOTCH_dgraphInit (grafptr, commdat);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHEXIT, scotchfdgraphexit, (         \
SCOTCH_Dgraph * const       grafptr),           \
(grafptr))
{
  SCOTCH_dgraphExit (grafptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHSIZE, scotchfdgraphsize, (         \
const SCOTCH_Dgraph * const grafptr,            \
SCOTCH_Num * const          vertglbptr,         \
SCOTCH_Num * const          vertlocptr,         \
SCOTCH_Num * const          edgeglbptr,         \
SCOTCH_Num * const          edgelocptr),        \
(grafptr, vertglbptr, vertlocptr, edgeglbptr, edgelocptr))
{
  SCOTCH_dgraphSize (grafptr, vertglbptr, vertlocptr, edgeglbptr, edgelocptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHDATA, scotchfdgraphdata, (         \
const SCOTCH_Dgraph * const grafptr,            \
const SCOTCH_Num * const    indxptr,            \
SCOTCH_Num * const          baseptr,            \
SCOTCH_Num * const          vertglbptr,         \
SCOTCH_Num * const          vertlocptr,         \
SCOTCH_Num * const          vertlocptz,         \
SCOTCH_Num * const          vertgstptr,         \
SCOTCH_Idx * const          vertlocidx,         \
SCOTCH_Idx * const          vendlocidx,         \
SCOTCH_Idx * const          velolocidx,         \
SCOTCH_Idx * const          vlbllocidx,         \
SCOTCH_Num * const          edgeglbptr,         \
SCOTCH_Num * const          edgelocptr,         \
SCOTCH_Num * const          edgelocptz,         \
SCOTCH_Idx * const          edgelocidx,         \
SCOTCH_Idx * const          edgegstidx,         \
SCOTCH_Idx * const          edlolocidx,         \
MPI_Fint * const            commptr),           \
(grafptr, indxptr, baseptr,                     \
 vertglbptr, vertlocptr, vertlocptz,            \
 vertgstptr, vertlocidx, vendlocidx,            \
 velolocidx, vlbllocidx, edgeglbptr,            \
 edgelocptr, edgelocptz, edgelocidx,            \
 edgegstidx, edlolocidx, commptr))
{
  SCOTCH_Num *        vertloctab;                 /* Pointer to graph arrays */
  SCOTCH_Num *        vendloctab;
  SCOTCH_Num *        veloloctab;
  SCOTCH_Num *        vlblloctab;
  SCOTCH_Num *        edgeloctab;
  SCOTCH_Num *        edgegsttab;
  SCOTCH_Num *        edloloctab;
  MPI_Comm            commdat;

  SCOTCH_dgraphData (grafptr, baseptr, vertglbptr, vertlocptr, vertlocptz, vertgstptr,
                     &vertloctab, &vendloctab, &veloloctab, &vlblloctab,
                     edgeglbptr, edgelocptr, edgelocptz,
                     &edgeloctab, &edgegsttab, &edloloctab, &commdat);
  *vertlocidx = (vertloctab - indxptr) + 1;       /* Add 1 since Fortran indices start at 1 */
  *vendlocidx = (vendloctab - indxptr) + 1;
  *velolocidx = (veloloctab != NULL) ? (veloloctab - indxptr) + 1 : *vertlocidx;
  *vlbllocidx = (vlblloctab != NULL) ? (vlblloctab - indxptr) + 1 : *vertlocidx;
  *edgelocidx = (edgeloctab - indxptr) + 1;
  *edgegstidx = (edgegsttab != NULL) ? (edgegsttab - indxptr) + 1 : *vertlocidx;
  *edlolocidx = (edloloctab != NULL) ? (edloloctab - indxptr) + 1 : *vertlocidx;
  *commptr = MPI_Comm_c2f (commdat);
}
