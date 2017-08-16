/* Copyright 2007,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_order_gather_f.c         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                distributed ordering gathering routines **/
/**                of the libSCOTCH library.               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 21 jul 2007     **/
/**                                 to     22 jul 2007     **/
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
/* for the ordering routines.         */
/*                                    */
/**************************************/

FORTRAN (                                           \
SCOTCHFDGRAPHCORDERINIT, scotchfdgraphcorderinit, ( \
const SCOTCH_Dgraph * const grafptr,                \
SCOTCH_Ordering * const     ordeptr,                \
SCOTCH_Num * const          permtab,                \
SCOTCH_Num * const          peritab,                \
SCOTCH_Num * const          cblkptr,                \
SCOTCH_Num * const          rangtab,                \
SCOTCH_Num * const          treetab,                \
int * const                 revaptr),               \
(grafptr, ordeptr, permtab, peritab,                \
 cblkptr, rangtab, treetab, revaptr))
{
  *revaptr = SCOTCH_dgraphCorderInit (grafptr, ordeptr, permtab, peritab, cblkptr, rangtab, treetab);
}

/*
**
*/

FORTRAN (                                           \
SCOTCHFDGRAPHCORDEREXIT, scotchfdgraphcorderexit, ( \
const SCOTCH_Dgraph * const grafptr,                \
SCOTCH_Ordering * const     ordeptr),               \
(grafptr, ordeptr))
{
  SCOTCH_dgraphCorderExit (grafptr, ordeptr);
}

/*
**
*/

FORTRAN (                                             \
SCOTCHFDGRAPHORDERGATHER, scotchfdgraphordergather, ( \
const SCOTCH_Dgraph * const     grafptr,              \
const SCOTCH_Dordering * const  dordptr,              \
SCOTCH_Ordering * const         cordptr,              \
int * const                     revaptr),             \
(grafptr, dordptr, cordptr, revaptr))
{
  *revaptr = SCOTCH_dgraphOrderGather (grafptr, dordptr, cordptr);
}
