/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : symbol.h                                **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                David GOUDIN (v0.0)                     **/
/**                Pascal HENON (v0.0)                     **/
/**                Pierre RAMET (v0.0)                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 21 mar 2002     **/
/**                                 to     21 mar 2002     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 apr 2003     **/
/**                                 to     10 jun 2003     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     03 mar 2005     **/
/**                # Version 5.1  : from : 05 nov 2010     **/
/**                                 to     05 nov 2010     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_H
#define SYMBOL_VERSION              1

/*
**  The type and structure definitions.
*/

/*+ The column block structure. +*/

typedef struct SymbolCblk_ {
  INT                       fcolnum;              /*+ First column index               +*/
  INT                       lcolnum;              /*+ Last column index (inclusive)    +*/
  INT                       bloknum;              /*+ First block in column (diagonal) +*/
} SymbolCblk;

/*+ The column block structure. +*/

typedef struct SymbolBlok_ {
  INT                       frownum;              /*+ First row index            +*/
  INT                       lrownum;              /*+ Last row index (inclusive) +*/
  INT                       cblknum;              /*+ Facing column block        +*/
  INT                       levfval;              /*+ Level-of-fill value        +*/
} SymbolBlok;

/*+ The symbolic block matrix. +*/

typedef struct SymbolMatrix_ {
  INT                       baseval;              /*+ Base value for numberings         +*/
  INT                       cblknbr;              /*+ Number of column blocks           +*/
  INT                       bloknbr;              /*+ Number of blocks                  +*/
  SymbolCblk *              cblktab;              /*+ Array of column blocks [+1,based] +*/
  SymbolBlok *              bloktab;              /*+ Array of blocks [based]           +*/
  INT                       nodenbr;              /*+ Number of nodes in matrix         +*/
} SymbolMatrix;

/*+ The type of cost computations. +*/

typedef enum SymbolCostType_ {
  SYMBOLCOSTLDLT                                  /*+ Crout (i.e. LDLt) cost function +*/
} SymbolCostType;

/* Structure for keeping track of selected
   blocks in the matrix pattern. The values
   of the tables are the remaining values
   for the yet unselected blocks.           */

typedef struct SymbolKeepBlok_ {
  INT                       levfval;              /*+ Values for incomplete factorisation +*/
  INT                       nupdval;
  INT                       ctrival;
  INT                       ctroval;
  INT                       hghtval;
} SymbolKeepBlok;

typedef struct SymbolKeep_ {
  INT                       levfmax;              /*+ Maximum values for incomplete fax +*/
  INT                       nupdmax;
  INT                       ctrimax;
  INT                       ctromax;
  INT                       hghtmax;
  byte *                    keeptab;              /*+ Flag array for kept blocks      +*/
  SymbolKeepBlok *          kblktab;              /*+ Block parameter array           +*/
  double *                  levftab;              /*+ Area arrays for selected blocks +*/
  double *                  nupdtab;
  double *                  ctritab;
  double *                  ctrotab;
  double *                  hghttab;
} SymbolKeep;

/*
**  The function prototypes.
*/

#ifndef SYMBOL
#define static
#endif

int                         symbolInit          (SymbolMatrix * const symbptr);
void                        symbolExit          (SymbolMatrix * const symbptr);
void                        symbolRealloc       (SymbolMatrix * const symbptr);
int                         symbolLoad          (SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolSave          (const SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolCheck         (const SymbolMatrix * const symbptr);
int                         symbolDraw          (const SymbolMatrix * const symbptr, FILE * const stream);
int                         symbolDrawFunc      (const SymbolMatrix * const symbptr, int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const), int (*) (const SymbolMatrix * const, const SymbolBlok * const, void * const, float * const), void * const, FILE * const stream);
void                        symbolDrawColor     (const INT labl, float * const coloptr);
#ifdef DOF_H
int                         symbolCost          (const SymbolMatrix * const symbptr, const Dof * const deofptr, const SymbolCostType typeval, double * const nnzptr, double * const opcptr);
int                         symbolCosti         (const SymbolMatrix * const symbptr, const Dof * const deofptr, const SymbolCostType typeval, const INT levfval, double * const nnzptr, double * const opcptr);
int                         symbolLevf          (const SymbolMatrix * const symbptr, INT * const levfmax, INT ** const levftab);
int                         symbolTree          (const SymbolMatrix * const symbptr, const Dof * const deofptr, INT * const leafnbr, INT * const heigmin, INT * const heigmax, double * const heigavg, double * const heigdlt);
int                         symbolNonzeros      (const SymbolMatrix * const symbptr, FILE * const stream);
#endif /* DOF_H */

int                         symbolKeepInit      (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr);
void                        symbolKeepExit      (SymbolKeep * restrict const keepptr);
void                        symbolKeepAdd       (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr, int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
void                        symbolKeepDel       (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr, int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
int                         symbolKeepCompute   (SymbolKeep * restrict const keepptr, const SymbolMatrix * const symbptr);
int                         symbolKeepHisto     (SymbolKeep * const keepptr, const SymbolMatrix * const, int (* funcptr) (const SymbolKeepBlok * const, void * const), void * dataptr);
int                         symbolKeepPurge     (SymbolKeep * restrict const keepptr, SymbolMatrix * restrict const symbptr);
int                         symbolKeepView      (const SymbolKeep * const keepptr, const double nnzlmax, const char * const nameptr);

#undef static
