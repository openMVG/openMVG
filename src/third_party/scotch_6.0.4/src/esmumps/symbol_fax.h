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
/**   NAME       : symbol_fax.h                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Jean ROMAN (v0.0)                       **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic factorization routine. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     24 sep 1998     **/
/**                # Version 0.1  : from : 04 apr 1999     **/
/**                                 to     01 may 1999     **/
/**                # Version 1.0  : from : 01 jun 2002     **/
/**                                 to     05 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 3.0  : from : 03 mar 2004     **/
/**                                 to     03 mar 2004     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_FAX_H

/*
**  The defines.
*/

/* Prime number for hashing vertex numbers. */

#define SYMBOL_FAX_HASHPRIME        17            /*+ Prime number for hashing +*/

/*
**  The type and structure definitions.
*/

/*+ The chained column block structure. These
    blocks are chained in a single linked list
    for block merge with blocks of left columns. +*/

typedef struct SymbolFaxTlok_ {
  INT                       frownum;              /*+ First row index            +*/
  INT                       lrownum;              /*+ Last row index (inclusive) +*/
  INT                       cblknum;              /*+ Facing column block        +*/
  INT                       nextnum;              /*+ Index of next block        +*/
} SymbolFaxTlok;
