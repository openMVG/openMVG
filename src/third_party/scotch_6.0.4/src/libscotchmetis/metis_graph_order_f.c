/* Copyright 2007,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_order_f.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API of   **/
/**                the compatibility library for the       **/
/**                MeTiS ordering routines.                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 sep 2006     **/
/**                                 to     10 sep 2006     **/
/**                # Version 5.1  : from : 30 jun 2010     **/
/**                                 to     30 jun 2010     **/
/**                # Version 6.0  : from : 13 sep 2012     **/
/**                                 to     13 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "common.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the graph ordering routines.   */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                             \
METISNAMEU(METIS_EDGEND), METISNAMEL(metis_edgend), ( \
const SCOTCH_Num * const    n,                        \
const SCOTCH_Num * const    xadj,                     \
const SCOTCH_Num * const    adjncy,                   \
const SCOTCH_Num * const    numflag,                  \
const SCOTCH_Num * const    options,                  \
SCOTCH_Num * const          perm,                     \
SCOTCH_Num * const          iperm),                   \
(n, xadj, adjncy, numflag, options, perm, iperm))
{
  METISNAMEU(METIS_EdgeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                                             \
METISNAMEU(METIS_NODEND), METISNAMEL(metis_nodend), ( \
const SCOTCH_Num * const    n,                        \
const SCOTCH_Num * const    xadj,                     \
const SCOTCH_Num * const    adjncy,                   \
const SCOTCH_Num * const    numflag,                  \
const SCOTCH_Num * const    options,                  \
SCOTCH_Num * const          perm,                     \
SCOTCH_Num * const          iperm),                   \
(n, xadj, adjncy, numflag, options, perm, iperm))
{
  METISNAMEU(METIS_NodeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/* When an input stream is built from the given
** file handle, it is set as unbuffered, so as to
** allow for multiple stream reads from the same
** file handle. If it were buffered, too many
** input characters would be read on the first
** block read.
*/

FORTRAN (                                               \
METISNAMEU(METIS_NODEWND), METISNAMEL(metis_nodewnd), ( \
const SCOTCH_Num * const    n,                          \
const SCOTCH_Num * const    xadj,                       \
const SCOTCH_Num * const    adjncy,                     \
const SCOTCH_Num * const    vwgt,                       \
const SCOTCH_Num * const    numflag,                    \
const SCOTCH_Num * const    options,                    \
SCOTCH_Num * const          perm,                       \
SCOTCH_Num * const          iperm),                     \
(n, xadj, adjncy, vwgt, numflag, options, perm, iperm))
{
  METISNAMEU(METIS_NodeWND) (n, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
