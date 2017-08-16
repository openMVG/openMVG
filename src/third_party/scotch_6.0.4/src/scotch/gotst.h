/* Copyright 2004,2007,2008,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gotst.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a graph scalar factorizer.      **/
/**                This module contains the data declara-  **/
/**                tions for the main module.              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 27 jan 2004     **/
/**                                 to   : 27 jan 2004     **/
/**                # Version 5.0  : from : 25 jun 2007     **/
/**                                 to   : 25 jul 2007     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines
*/

/*+ File name aliases. +*/

#define C_FILENBR                   3            /* Number of files in list                */
#define C_FILEARGNBR                3            /* Number of files which can be arguments */

#define C_filenamegrfinp            fileBlockName (C_fileTab, 0) /* Graph input file name    */
#define C_filenamesrcout            fileBlockName (C_fileTab, 1) /* Ordering input file name */
#define C_filenamedatout            fileBlockName (C_fileTab, 2) /* Output data file name    */

#define C_filepntrgrfinp            fileBlockFile (C_fileTab, 0) /* Graph input file     */
#define C_filepntrordinp            fileBlockFile (C_fileTab, 1) /* Ordering output file */
#define C_filepntrdatout            fileBlockFile (C_fileTab, 2) /* Output data file     */

/*
**  The type and structure definitions.
*/

/* Factorization node */

typedef struct C_FactorNode_ {
  struct C_FactorNode_ *    linkdad;              /*+ Father of node    +*/
  struct C_FactorNode_ *    linkson;              /*+ First son of node +*/
  struct C_FactorNode_ *    linkbro;              /*+ Brother of node   +*/
} C_FactorNode;

/* Data structure for computing factored matrix statistics. */

typedef struct FactorStat_ {
  const SCOTCH_Num *      ldadtax;
  const SCOTCH_Num *      lsontax;
  const SCOTCH_Num *      lbrotax;
  SCOTCH_Num              heigmin;
  SCOTCH_Num              heigmax;
  SCOTCH_Num              heignbr;
  double                  heigavg;
  double                  heigdlt;
  const SCOTCH_Num *      fnnztax;
  double                  fnnzsum;
} FactorStat;

/*
**  The function prototypes.
*/

static int                  factorView          (const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, FILE * restrict const);
static int                  factorView2         (const SCOTCH_Num, const SCOTCH_Num, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, const SCOTCH_Num * const, SCOTCH_Num * restrict, SCOTCH_Num * restrict, SCOTCH_Num * restrict, SCOTCH_Num * restrict);
static void                 factorView3         (FactorStat * restrict const, SCOTCH_Num, SCOTCH_Num, double * restrict const);
static void                 factorView4         (FactorStat * restrict const, SCOTCH_Num, SCOTCH_Num, double * restrict const);
