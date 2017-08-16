/* Copyright 2007 ENSEIRB, INRIA & CNRS
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
/**********************************************************/
/*                                                        */
/*   NAME       : dgraph_allreduce.h                      */
/*                                                        */
/*   AUTHOR     : Francois PELLEGRINI                     */
/*                                                        */
/*   FUNCTION   : These lines are the data declarations   */
/*                for the communication routines          */
/*                                                        */
/*                # Version 5.0  : from : 28 aug 2006     */
/*                                 to     29 aug 2006     */
/*                                                        */
/**********************************************************/

/*
**  The defines.
*/

/*+ Combined maximum-sum reduction operator +*/

#define DGRAPHALLREDUCEMAXSUMOP(m,s)                                                                \
static                                                                                              \
void                                                                                                \
dgraphAllreduceMaxSumOp##m##_##s (                                                                  \
const Gnum * const          in,                   /* First operand                               */ \
Gnum * const                inout,                /* Second and output operand                   */ \
const int * const           len,                  /* Number of instances ; should be 1, not used */ \
const MPI_Datatype * const  typedat)              /* MPI datatype ; not used                     */ \
{                                                                                                   \
  int               i;                                                                              \
                                                                                                    \
  for (i = 0; i < (m); i ++)                      /* Perform maximum on first part of data array */ \
    if (in[i] > inout[i])                                                                           \
      inout[i] = in[i];                                                                             \
                                                                                                    \
  for ( ; i < ((m) + (s)); i ++)                  /* Perform sum on second part of data array */    \
    inout[i] += in[i];                                                                              \
}

#define dgraphAllreduceMaxSum(rlt,rgt,m,s,comm) dgraphAllreduceMaxSum2 ((rlt), (rgt), (m) + (s), (MPI_User_function *) (dgraphAllreduceMaxSumOp##m##_##s), (comm))

/*
** The function prototypes.
*/

int                         dgraphAllreduceMaxSum2 (Gnum *, Gnum *, int, MPI_User_function *, MPI_Comm);
