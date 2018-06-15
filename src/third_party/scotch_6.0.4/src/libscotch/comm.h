/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : comm.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for communication functions.            **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 30 jul 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

#define COMM_H

/*
**  The type and structure definitions.
*/

#ifndef GNUMMAX                                   /* If dgraph.h not included    */
typedef INT                   Gnum;               /* Vertex and edge numbers     */
typedef UINT                  Gunum;              /* Unsigned type of same width */
#define GNUMMAX                     (INTVALMAX)   /* Maximum signed Gnum value   */
#define GNUMSTRING                  INTSTRING     /* String to printf a Gnum     */
#endif /* GNUMMAX */

/*
**  The function prototypes.
*/

#ifndef COMM
#define static
#endif

int                         commAllgatherv      (void * const, const Gnum, MPI_Datatype, void * const, const Gnum * const, const Gnum * const, MPI_Datatype, MPI_Comm);
int                         commGatherv         (void * const, const Gnum, MPI_Datatype, void * const, const Gnum * const, const Gnum * const, MPI_Datatype, const int, MPI_Comm);
int                         commScatterv        (void * const, const Gnum * const, const Gnum * const, MPI_Datatype, void * const, const Gnum, MPI_Datatype, const int, MPI_Comm);

#undef static

/*
**  The macro definitions.
*/

#ifndef COMM
#ifndef INTSIZE64
#define commAllgatherv              MPI_Allgatherv
#define commGatherv                 MPI_Gatherv
#define commScatterv                MPI_Scatterv
#endif /* INTSIZE64 */
#endif /* COMM      */
