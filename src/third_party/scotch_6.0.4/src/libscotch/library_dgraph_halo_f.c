/* Copyright 2007,2009,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_halo_f.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the distributed source graph handling   **/
/**                routines of the libSCOTCH library.      **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 17 jul 2007     **/
/**                                 to     02 aug 2007     **/
/**                # Version 5.1  : from : 09 may 2009     **/
/**                                 to     10 may 2009     **/
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
/* for the distributed graph handling */
/* routines.                          */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHGHST, scotchfdgraphghst, (         \
SCOTCH_Dgraph * const       grafptr,            \
int * const                 revaptr),           \
(grafptr, revaptr))
{
  *revaptr = SCOTCH_dgraphGhst (grafptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHHALO, scotchfdgraphhalo, (         \
SCOTCH_Dgraph * const       grafptr,            \
void * const                datatab,            \
MPI_Fint * const            typeptr,            \
int * const                 revaptr),           \
(grafptr, datatab, typeptr, revaptr))
{
  MPI_Datatype        typeval;

  typeval = MPI_Type_f2c (*typeptr);
  *revaptr = SCOTCH_dgraphHalo (grafptr, datatab, typeval);
}

/*
**
*/

FORTRAN (                                         \
SCOTCHFDGRAPHHALOASYNC, scotchfdgraphhaloasync, ( \
SCOTCH_Dgraph * const         grafptr,            \
void * const                  datatab,            \
MPI_Fint * const              typeptr,            \
SCOTCH_DgraphHaloReq * const  requptr,            \
int * const                   revaptr),           \
(grafptr, datatab, typeptr, requptr, revaptr))
{
  MPI_Datatype        typeval;

  typeval = MPI_Type_f2c (*typeptr);
  *revaptr = SCOTCH_dgraphHaloAsync (grafptr, datatab, typeval, requptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFDGRAPHHALOWAIT, scotchfdgraphhalowait, ( \
SCOTCH_DgraphHaloReq * const  requptr,          \
int * const                   revaptr),         \
(requptr, revaptr))
{
  *revaptr = SCOTCH_dgraphHaloWait (requptr);
}
