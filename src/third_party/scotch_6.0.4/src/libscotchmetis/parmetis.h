/*********************************************************
**                                                      **
**  WARNING: THIS IS NOT THE ORIGINAL INCLUDE FILE OF   **
**  THE ParMeTiS SOFTWARE PACKAGE.                      **
**  This file is a compatibility include file provided  **
**  as part of the Scotch software distribution.        **
**  Preferably use the original ParMeTiS include file   **
**  to keep definitions of routines not overloaded by   **
**  the libPTScotchMeTiS library.                       **
**                                                      **
*********************************************************/
/* Copyright 2007,2008,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : parmetis.h                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Compatibility declaration file for the  **/
/**                MeTiS interface routines provided by    **/
/**                the Scotch project.                     **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 17 oct 2007     **/
/**                                 to     18 oct 2007     **/
/**                # Version 5.1  : from : 19 jun 2008     **/
/**                                 to     30 jun 2010     **/
/**                # Version 6.0  : from : 13 sep 2012     **/
/**                                 to     13 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifndef __parmetis_h__
#define __parmetis_h__

#include <mpi.h>                                  /* Since ParMeTiS does it, do it too */

#endif /* __parmetis_h__ */

#ifdef SCOTCH_METIS_PREFIX
#define SCOTCH_METIS_PREFIXL        scotch_
#define SCOTCH_METIS_PREFIXU        SCOTCH_
#endif /* SCOTCH_METIS_PREFIX */

#ifndef SCOTCH_METIS_PREFIXL
#define SCOTCH_METIS_PREFIXL
#endif /* SCOTCH_METIS_PREFIXL */

#ifndef SCOTCH_METIS_PREFIXU
#define SCOTCH_METIS_PREFIXU
#endif /* SCOTCH_METIS_PREFIXU */

#ifndef METISNAMEL
#define METISNAMEL(s)               METISNAME2(METISNAME3(SCOTCH_METIS_PREFIXL),s)
#define METISNAMEU(s)               METISNAME2(METISNAME3(SCOTCH_METIS_PREFIXU),s)
#define METISNAME2(p,s)             METISNAME4(p,s)
#define METISNAME3(s)               s
#define METISNAME4(p,s)             p##s
#endif /* METISNAMEL */

/*
**  The function prototypes.
*/

void                        METISNAMEU(ParMETIS_V3_NodeND) (const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, MPI_Comm * const);
void                        METISNAMEU(ParMETIS_V3_PartGeomKway) (const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const float * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const float * const, const float * const, const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, MPI_Comm * const);
void                        METISNAMEU(ParMETIS_V3_PartKway) (const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const float * const, const float * const, const SCOTCH_Num * const, SCOTCH_Num * const, SCOTCH_Num * const, MPI_Comm * const);
