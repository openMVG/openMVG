/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : esmumps.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains a MUMPS interface  **/
/**                for the ordering routines of the        **/
/**                libSCOTCH + Emilio libfax libraries.    **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 16 may 2001     **/
/**                                 to     04 jun 2001     **/
/**                # Version 0.1  : from : 13 feb 2002     **/
/**                                 to     13 feb 2002     **/
/**                # Version 1.0  : from : 06 dec 2004     **/
/**                                 to     06 dec 2004     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ESMUMPS

#include "common.h"
#ifdef SCOTCH_PTSCOTCH
#include "ptscotch.h"
#else /* SCOTCH_PTSCOTCH */
#include "scotch.h"
#endif /* SCOTCH_PTSCOTCH */
#include "graph.h"
#include "dof.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "esmumps.h"

/**************************************/
/*                                    */
/* This routine acts as an interface  */
/* between ordering software such as  */
/* MUMPS and Scotch+Emilio.           */
/*                                    */
/**************************************/

/* Meaning of the parameters :
** - n : order of the system (that is, number of columns).
** - iwlen : not used. Here for compatibility.
** - pe : on input, position in array iw of the extra-diagonal
**   terms for the considered column.
**   on output, -pe(i) is the father of node i in the elimination
**   tree if i is a principal variable, or it is the index of
**   the principal variable if i is a secondary variable.
** - pfree : number of extra-diagonal terms for the considered
**   node (that is, the number of arcs dans le graph for this
**   vertex).
** - len : array holding the number of extra-diagonal terms for
**   each column.
** - iw : array of extra-diagonal terms (preserved).
** - nv : on output, nv(i) = 0 if variable i is a secondary
**   variable, else nv(i) is the number of columns that are
**   merged into principal variable i.
** - elen : on output, direct permutation (for MUMPS; the
**   meaning of the "direct" and "inverse" permutations is
**   just the opposite for Scotch) :
**   k=elen(i) <==> column i is the k-th pivot.
** - last : on output, inverse permutation (for MUMPS) :
**   i=last(k) <==> column i est le k-th pivot
*/

int
esmumps (
const INT                   n,
const INT                   iwlen,                /* Not used, just here for consistency */
INT * restrict const        petab,
const INT                   pfree,
INT * restrict const        lentab,
INT * restrict const        iwtab,
INT * restrict const        nvtab,
INT * restrict const        elentab,              /* Permutations computed for debugging only */
INT * restrict const        lasttab)              /* Permutations computed for debugging only */
{
  INT                 baseval;                    /* Base value            */
  INT * restrict      vendtab;                    /* Vertex end array      */
  Graph               grafdat;                    /* Graph                 */
  Order               ordedat;                    /* Graph ordering        */
  SymbolMatrix        symbdat;                    /* Block factored matrix */
  Dof                 deofdat;                    /* Matrix DOF structure  */
  INT                 vertnum;
  INT                 cblknum;
  INT                 colnum;

  if ((vendtab = memAlloc (n * sizeof (INT))) == NULL) {
    errorPrint ("esmumps: out of memory");
    return     (1);
  }
  for (vertnum = 0; vertnum < n; vertnum ++)
    vendtab[vertnum] = petab[vertnum] + lentab[vertnum];

  baseval = 1;                                    /* Assume Fortran-based indexing */
  graphInit        (&grafdat);
  graphBuildGraph2 (&grafdat, baseval, n, pfree - 1, petab, vendtab, NULL, NULL, iwtab, NULL);

  dofInit     (&deofdat);
  dofConstant (&deofdat, 1, n, 1);                /* One DOF per node, Fortran-based indexing */

  orderInit  (&ordedat);
  orderGraph (&ordedat, &grafdat);                /* Compute ordering with Scotch */

#ifdef ESMUMPS_DEBUG                              /* Permutations are output for debugging only */
  memCpy (elentab, ordedat.permtab, n * sizeof (INT)); /* Copy permutations                     */
  memCpy (lasttab, ordedat.peritab, n * sizeof (INT));
#endif /* ESMUMPS_DEBUG */
 
  symbolInit     (&symbdat);
  symbolFaxGraph (&symbdat, &grafdat, &ordedat);  /* Compute block symbolic factorizaion */

  for (cblknum = 0; cblknum < symbdat.cblknbr; cblknum ++) { /* For all column blocks */
    INT                 degnbr;                   /* True degree of column block      */
    INT                 bloknum;

    for (bloknum = symbdat.cblktab[cblknum].bloknum, degnbr = 0;
         bloknum < symbdat.cblktab[cblknum + 1].bloknum; bloknum ++)
      degnbr += symbdat.bloktab[bloknum - baseval].lrownum -
                symbdat.bloktab[bloknum - baseval].frownum + 1;
    nvtab[ordedat.peritab[symbdat.cblktab[cblknum].fcolnum - baseval] - baseval] = degnbr; /* Set true block degree */

    for (colnum  = symbdat.cblktab[cblknum].fcolnum + 1; /* For all secondary variables */
         colnum <= symbdat.cblktab[cblknum].lcolnum; colnum ++) {
      nvtab[ordedat.peritab[colnum - baseval] - baseval] = 0; /* Set nv = 0 and pe = - principal variable */
      petab[ordedat.peritab[colnum - baseval] - baseval] =
        - ordedat.peritab[symbdat.cblktab[cblknum].fcolnum - baseval];
    }

    if (symbdat.cblktab[cblknum].bloknum ==       /* If column block has no extra-diagonals */
        symbdat.cblktab[cblknum + 1].bloknum - 1) /* Then mark block as root of subtree     */
      petab[ordedat.peritab[symbdat.cblktab[cblknum].fcolnum - baseval] - baseval] = 0;
    else
      petab[ordedat.peritab[symbdat.cblktab[cblknum].fcolnum - baseval] - baseval] =
        - ordedat.peritab[symbdat.cblktab[symbdat.bloktab[symbdat.cblktab[cblknum].bloknum + 1 - baseval].cblknum - baseval].fcolnum - baseval];
  }

  symbolExit (&symbdat);
  orderExit  (&ordedat);
  dofExit    (&deofdat);
  graphExit  (&grafdat);

  memFree (vendtab);

  return (0);
}
