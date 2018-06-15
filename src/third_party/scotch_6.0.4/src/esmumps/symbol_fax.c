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
/**   NAME       : symbol_fax.c                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Jean ROMAN (v0.0)                       **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This is the generic block symbolic      **/
/**                factorization routine.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     29 sep 1998     **/
/**                # Version 0.1  : from : 04 apr 1999     **/
/**                                 to     21 apr 1999     **/
/**                # Version 0.2  : from : 08 may 2000     **/
/**                                 to     09 may 2000     **/
/**                # Version 1.0  : from : 13 mar 2002     **/
/**                                 to     08 jun 2002     **/
/**                # Version 1.2  : from : 23 aug 2002     **/
/**                                 to     23 aug 2002     **/
/**                # Version 2.0  : from : 21 mar 2003     **/
/**                                 to     21 mar 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifndef SYMBOL_FAX_INCLUDED                       /* If included from other file */
#define SYMBOL_FAX

#include "common.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "symbol_fax.h"
#endif /* SYMBOL_FAX_INCLUDED */

/***********************************/
/*                                 */
/* Symbolic factorization routine. */
/*                                 */
/***********************************/

/*+ This routine computes the block symbolic
*** factorization of the given matrix
*** according to the given vertex ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
*** Algorithm:
*** The algorithm is implemented in a
*** cache-friendly manner, by using a single
*** dynamic array which grows along with the
*** number of computed blocks. The array is
*** decomposed in the following manner:
*** - In a first phase, a hash table and a
***   sort area are reserved at the end of
***   the space of already computed blocks.
***   The sort area is created far enough from
***   the end of the array of already computed
***   blocks such that if there are no contributing
***   blocks all new blocks can be created without
***   colliding with the sort area.
*** - Then, in a second phase, if the current
***   column block does have contributing column
***   blocks, an area for simply-linked temporary
***   blocks is reserved at least after the sort area,
***   leaving enough space to create all of the
***   corresponding potential new blocks
***   just after all the blocks of the previous
***   column block (right picture).
***      ___________
***     |ccccccccccc| <- bloktab (bloktax)
***     |ccccccccccc|
***     |ccccccccccc|                                :ccccccccccc:
***     |ccccccccccc| >- Computed blocks ----------< |ccccccccccc|
***     |ccccccccccc|                                |ccccccccccc|
***     |-----------|                                |:::::::::::|
***     |hhhhhhhhhhh| <- hashtab = bloknum --------> |bcbcbcbcbcb|
***     |hhhhhhhhhhh|                 |              |cbcbcbcbcbc|
***     |hhhhhhhhhhh|                 |              |bcbcbcbcbcb|
***     |hhhhhhhhhhh|                 |              |cbcbcbcbcbc|
***     |-----------|                 |              |bcbcbcbcbcb|
***     |           |                 |              |-----------|
***     |-----------| <- sorttab...... ------------> |           |
***     |sssssssssss|                                |           |
***     |sssssssssss|                                |           |
***     |-----------| <- ............................|           |
***     |           |                     tloktab -> |-----------|
***     |           |                                |ttttttttttt|
***     |           |                                |ttttttttttt|
***     :           :                                |-----------|
***     :___________:                                :___________:
***                   <- bloktab + blokmax
+*/

#ifndef SYMBOL_FAX_INCLUDED
#define SYMBOL_FAX_ITERATOR(ngbdptr, vertnum, vertend) \
                                    for (vertend  = ngbfrst ((ngbdptr), (vertnum)); \
                                         vertend >= baseval;                        \
                                         vertend  = ngbnext (ngbdptr)) {

#define SYMBOL_FAX_VERTEX_DEGREE(ngbdptr, vertnum) \
                                    (ngbdegr ((ngbdptr), (vertnum)))

int
symbolFax (
SymbolMatrix * const        symbptr,              /*+ Symbolic block matrix [based]                      +*/
const INT                   vertnbr,              /*+ Number of vertices                                 +*/
const INT                   edgenbr,              /*+ Number of edges                                    +*/
const INT                   baseval,              /*+ Base value                                         +*/
void * const                ngbdptr,              /*+ Neighbor bookkeeping area                          +*/
INT                         ngbfrst (void * const, const INT), /*+ First neighbor function               +*/
INT                         ngbnext (void * const), /*+ Next neighbor function                           +*/
INT                         ngbdegr (void * const, const INT),  /*+ Vertex degree function (upper bound) +*/
const Order * const         ordeptr)              /*+ Matrix ordering                                    +*/
#endif /* SYMBOL_FAX_INCLUDED */
{
  INT                       vertnum;              /* Vertex number of current column                   */
  INT                       vertend;              /* Current end vertex number                         */
  const INT * restrict      permtax;              /* Based access to direct permutation array          */
  const INT * restrict      peritax;              /* Based access to inverse permutation array         */
  const INT * restrict      rangtax;              /* Based access to column block range array          */
  INT * restrict            ctrbtax;              /* Based access to array of contribution chains      */
  SymbolCblk * restrict     cblktax;              /* Based access to column block array                */
  INT                       cblknum;              /* Based number of current column block              */
  INT                       cblkctr;              /* Based number of current contributing column block */
  SymbolBlok * restrict     bloktax;              /* Based access to block array                       */
  INT                       bloknum;              /* Based number of current first free block slot     */
  INT                       blokmax;              /* Maximum number of blocks in array                 */
  SymbolFaxTlok * restrict  tloktab;              /* Beginning of array of temporary blocks            */
  INT                       ctrbsum;              /* Number of contributing blocks for column block    */
  INT * restrict            sorttab;              /* Beginning of sort area                            */
  INT                       sortnbr;              /* Number of vertices in sort area and hash table    */
  INT * restrict            hashtab;              /* Hash vertex table                                 */
  INT                       hashmsk;              /* Mask for access to hash table                     */
  INT                       colend;               /* Column number of vertex neighbor                  */

  permtax = ordeptr->permtab - baseval;           /* Compute array bases */
  peritax = ordeptr->peritab - baseval;
  rangtax = ordeptr->rangtab - baseval;

  blokmax  = ordeptr->cblknbr * (2 + edgenbr / vertnbr) + 2; /* Estimate size of initial block array */

  {                                               /* Allocate arrays for factoring   */
    INT *               ctrbtab;                  /* Array for contribution chaining */
    SymbolCblk *        cblktab;                  /* Column block array              */
    SymbolBlok *        bloktab;                  /* Block array                     */

    if (((ctrbtab = (INT *)        memAlloc (ordeptr->cblknbr       * sizeof (INT)))        == NULL) ||
        ((cblktab = (SymbolCblk *) memAlloc ((ordeptr->cblknbr + 1) * sizeof (SymbolCblk))) == NULL) ||
        ((bloktab = (SymbolBlok *) memAlloc (blokmax                * sizeof (SymbolBlok))) == NULL)) {
      errorPrint ("symbolFax: out of memory (1)");
      if (ctrbtab != NULL) {
        if (cblktab != NULL)
          memFree (cblktab);
        memFree (ctrbtab);
      }
      return (1);
    }
    cblktax = cblktab - baseval;                  /* Set based accesses */
    bloktax = bloktab - baseval;
    ctrbtax = ctrbtab - baseval;

    memset (ctrbtab, ~0, ordeptr->cblknbr * sizeof (INT)); /* Initialize column block contributions link array */
  }

  bloknum = baseval;
  for (cblknum = baseval; cblknum < baseval + ordeptr->cblknbr; cblknum ++) { /* For all column blocks */
    INT                 colnum;                   /* Number of current column [based]                  */
    INT                 colmax;                   /* Maximum column index for current column block     */

    {                                             /* Compute offsets and check for array size */
      INT                 degrsum;
      INT                 hashsiz;
      INT                 hashmax;
      INT                 ctrbtmp;
      INT                 sortoft;                /* Offset of sort array                   */
      INT                 tlokoft;                /* Offset of temporary block array        */
      INT                 tlndoft;                /* Offset of end of temporary block array */
      INT                 tlokmax;

      colnum = rangtax[cblknum];
      colmax = rangtax[cblknum + 1];              /* Get maximum column value */

      cblktax[cblknum].fcolnum = colnum;          /* Set column block data */
      cblktax[cblknum].lcolnum = colmax - 1;
      cblktax[cblknum].bloknum = bloknum;

      degrsum = 0;
      for ( ; colnum < colmax; colnum ++) /* For all columns                                  */
        degrsum += SYMBOL_FAX_VERTEX_DEGREE (ngbdptr, peritax[colnum]); /* Add column degrees */

      for (hashmax = 256; hashmax < degrsum; hashmax *= 2) ; /* Get upper bound on hash table size */
      hashsiz = hashmax << 2;                     /* Fill hash table at 1/4 of capacity            */
      hashmsk = hashsiz - 1;

      for (ctrbsum = 0, ctrbtmp = ctrbtax[cblknum]; /* Follow chain of contributing column blocks */
           ctrbtmp != ~0; ctrbtmp = ctrbtax[ctrbtmp])
        ctrbsum += cblktax[ctrbtmp + 1].bloknum - cblktax[ctrbtmp].bloknum - 2; /* Sum contributing column blocks */

      tlokmax = degrsum + ctrbsum;
      sortoft = tlokmax * sizeof (SymbolBlok);
      if ((hashsiz * sizeof (INT)) > sortoft)     /* Compute offset of sort area */
        sortoft = (hashsiz * sizeof (INT));
      tlokoft = sortoft + degrsum * sizeof (INT); /* Compute offset of temporary block area */
      tlndoft = tlokoft + tlokmax * sizeof (SymbolFaxTlok); /* Compute end of area          */

      if (((byte *) (bloktax + bloknum) + tlndoft) > /* If not enough room */
          ((byte *) (bloktax + blokmax))) {
        SymbolBlok *        bloktmp;              /* Temporary pointer for array resizing */

        do
          blokmax = blokmax + (blokmax >> 2) + 4; /* Increase block array size by 25% as long as it does not fit */
        while (((byte *) (bloktax + bloknum) + tlndoft) > ((byte *) (bloktax + blokmax)));

        if ((bloktmp = (SymbolBlok *) memRealloc (bloktax + baseval, (blokmax * sizeof (SymbolBlok)))) == NULL) {
          errorPrint ("symbolFax: out of memory (2)");
          memFree    (bloktax + baseval);
          memFree    (cblktax + baseval);
          memFree    (ctrbtax + baseval);
          return     (1);
        }
        bloktax = bloktmp - baseval;
      }

      hashtab = (INT *)           (bloktax + bloknum);
      sorttab = (INT *)           ((byte *) hashtab + sortoft);
      tloktab = (SymbolFaxTlok *) ((byte *) hashtab + tlokoft);

      memset (hashtab, ~0, hashsiz * sizeof (INT)); /* Initialize hash table */
    }

    sortnbr = 0;                                  /* No vertices yet                 */
    for (colnum = rangtax[cblknum]; colnum < colmax; colnum ++) { /* For all columns */
      INT                 hashnum;

      vertnum = peritax[colnum];                  /* Get associated vertex      */
      SYMBOL_FAX_ITERATOR (ngbdptr, vertnum, vertend) /* For all adjacent edges */
        colend = permtax[vertend];                /* Get end column number      */

        if (colend < colmax)                      /* If end vertex number in left columns */
          continue;                               /* Skip to next neighbor                */

        for (hashnum = (colend * SYMBOL_FAX_HASHPRIME) & hashmsk; ; /* Search end column in hash table */
             hashnum = (hashnum + 1) & hashmsk) {
          INT *               hashptr;

          hashptr = hashtab + hashnum;            /* Point to hash slot           */
          if (*hashptr == colend)                 /* If end column in hash table  */
            break;                                /* Skip to next end column      */
          if (*hashptr == ~0) {                   /* If slot is empty             */
            *hashptr = colend;                    /* Set column in hash table     */
            sorttab[sortnbr ++] = colend;         /* Add end column to sort array */
            break;
          }
        }
      }                                           /* End of loop on neighbors */
    }                                             /* End of loop on columns   */

    intSort1asc1 (sorttab, sortnbr);              /* Sort neighbor array */

    cblkctr = cblknum;
    if (ctrbtax[cblknum] == ~0) {                 /* If column is not to be updated */
      INT                 sortnum;

      bloktax[bloknum].frownum = cblktax[cblknum].fcolnum; /* Build diagonal block */
      bloktax[bloknum].lrownum = cblktax[cblknum].lcolnum;
      bloktax[bloknum].cblknum = cblknum;
      bloktax[bloknum].levfval = 0;
      bloknum ++;

      for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */
        INT                 colend;               /* Column number of current entry  */

        colend = sorttab[sortnum];
        if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
          INT                 cblktmm;            /* Median value                       */
          INT                 cblktmx;            /* Maximum value                      */

          for (cblkctr ++,                        /* Find new column block by dichotomy */
               cblktmx = ordeptr->cblknbr + baseval;
               cblktmx - cblkctr > 1; ) {
            cblktmm = (cblktmx + cblkctr) >> 1;
            if (rangtax[cblktmm] <= colend)
              cblkctr = cblktmm;
            else
              cblktmx = cblktmm;
          }
        }

        bloktax[bloknum].frownum = colend;        /* Set beginning of new block */
        while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
               (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
               (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
        bloktax[bloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
        bloktax[bloknum].cblknum = cblkctr;
        bloktax[bloknum].levfval = 0;
        bloknum ++;                               /* One more block */
      }
    }
    else {                                        /* Column will be updated           */
      INT                 sortnum;                /* Current index in sort array      */
      INT                 tloknum;                /* Current index on temporary block */
      INT                 tlokfre;                /* Index of first free block        */

      tloktab->frownum = cblktax[cblknum].fcolnum; /* Build diagonal chained block */
      tloktab->lrownum = cblktax[cblknum].lcolnum;
      tloktab->cblknum = cblknum;
      tloktab->nextnum = 1;
      tloknum = 1;

      for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */
        INT                 colend;               /* Column number of current entry  */

        colend = sorttab[sortnum];
        if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
          INT                 cblktmm;            /* Median value                       */
          INT                 cblktmx;            /* Maximum value                      */

          for (cblkctr ++,                        /* Find new column block by dichotomy */
               cblktmx = ordeptr->cblknbr + baseval;
               cblktmx - cblkctr > 1; ) {
            cblktmm = (cblktmx + cblkctr) >> 1;
            if (rangtax[cblktmm] <= colend)
              cblkctr = cblktmm;
            else
              cblktmx = cblktmm;
          }
        }
        tloktab[tloknum].frownum = colend;        /* Set beginning of new block */
        while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
               (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
               (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
        tloktab[tloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
        tloktab[tloknum].cblknum = cblkctr;
        tloktab[tloknum].nextnum = tloknum + 1;   /* Chain block */
        tloknum = tloknum + 1;
      }
      tloktab[tloknum].frownum =                  /* Build trailing block */
      tloktab[tloknum].lrownum = vertnbr + baseval;
      tloktab[tloknum].cblknum = ordeptr->cblknbr + baseval;
      tloktab[tloknum].nextnum = 0;               /* Set end of chain (never chain to diagonal block) */

      tlokfre = ++ tloknum;                       /* Build free chain for possible contributing blocks */
      for ( ; tloknum < tlokfre + ctrbsum; tloknum = tloknum + 1)
        tloktab[tloknum].nextnum = tloknum + 1;
      tloktab[tloknum].nextnum = ~0;              /* Set end of free chain */

      for (cblkctr = ctrbtax[cblknum]; cblkctr != ~0; cblkctr = ctrbtax[cblkctr]) { /* Follow chain */
        INT                 blokctr;              /* Current index of contributing column block     */
        INT                 tloklst;              /* Index of previous temporary block              */

        tloklst = 0;                              /* Previous is diagonal block */
        tloknum = 0;                              /* Current is diagonal block  */

        for (blokctr = cblktax[cblkctr].bloknum + 2; /* For all blocks in contributing column block */
             blokctr < cblktax[cblkctr + 1].bloknum; blokctr ++) {
          while ((tloktab[tloknum].cblknum < bloktax[blokctr].cblknum) || /* Skip unmatched chained blocks */
                 (tloktab[tloknum].lrownum < bloktax[blokctr].frownum - 1)) {
            tloklst = tloknum;
            tloknum = tloktab[tloknum].nextnum;
          }

          if ((bloktax[blokctr].cblknum < tloktab[tloknum].cblknum) || /* If contributing block has no mate */
              (bloktax[blokctr].lrownum < tloktab[tloknum].frownum - 1)) {
            INT                 tloktmp;

#ifdef FAX_DEBUG
            if (tlokfre == ~0) {
              errorPrint ("symbolFax: internal error (1)");
              memFree    (bloktax + baseval);
              memFree    (cblktax + baseval);
              memFree    (ctrbtax + baseval);
              return     (1);
            }
#endif /* FAX_DEBUG */
            tloktmp                  =
            tloktab[tloklst].nextnum = tlokfre;   /* Chain new block                */
            tloktab[tlokfre].frownum = bloktax[blokctr].frownum; /* Copy block data */
            tloktab[tlokfre].lrownum = bloktax[blokctr].lrownum;
            tloktab[tlokfre].cblknum = bloktax[blokctr].cblknum;
            tlokfre                  = tloktab[tlokfre].nextnum;
            tloktab[tloktmp].nextnum = tloknum;   /* Complete chainimg                    */
            tloknum                  = tloktab[tloklst].nextnum; /* Resume from new block */
            continue;                             /* Process next block                   */
          }

          if ((bloktax[blokctr].lrownum >= tloktab[tloknum].frownum - 1) && /* Update chained block lower bound */
              (bloktax[blokctr].frownum <  tloktab[tloknum].frownum))
            tloktab[tloknum].frownum = bloktax[blokctr].frownum;

          if ((bloktax[blokctr].frownum <= tloktab[tloknum].lrownum + 1) && /* Update chained block upper bound */
              (bloktax[blokctr].lrownum >  tloktab[tloknum].lrownum)) {
            INT                 tloktmp;

            tloktab[tloknum].lrownum = bloktax[blokctr].lrownum;

            for (tloktmp = tloktab[tloknum].nextnum; /* Aggregate following chained blocks */
                 (tloktab[tloktmp].cblknum == tloktab[tloknum].cblknum) &&
                 (tloktab[tloktmp].frownum <= tloktab[tloknum].lrownum + 1);
                 tloktmp = tloktab[tloknum].nextnum) {
              if (tloktab[tloktmp].lrownum > tloktab[tloknum].lrownum) /* Merge aggregated block */
                tloktab[tloknum].lrownum = tloktab[tloktmp].lrownum;
              tloktab[tloknum].nextnum = tloktab[tloktmp].nextnum; /* Unlink aggregated block */
              tloktab[tloktmp].nextnum = tlokfre;
              tlokfre                  = tloktmp;
            }
          }
        }
      }

      for (tloknum = 0;                           /* For all chained blocks                    */
           tloktab[tloknum].nextnum != 0;         /* Until trailer block is reached            */
           tloknum = tloktab[tloknum].nextnum, bloknum ++) { /* Copy block data to block array */
        bloktax[bloknum].frownum = tloktab[tloknum].frownum;
        bloktax[bloknum].lrownum = tloktab[tloknum].lrownum;
        bloktax[bloknum].cblknum = tloktab[tloknum].cblknum;
        bloktax[bloknum].levfval = 0;
      }
    }
    if ((bloknum - cblktax[cblknum].bloknum) > 2) { /* If more than one extra-diagonal blocks exist                 */
      ctrbtax[cblknum] = ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].cblknum]; /* Link contributing column blocks */
      ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].cblknum] = cblknum;
    }
  }
  cblktax[cblknum].fcolnum =                      /* Set last column block data */
  cblktax[cblknum].lcolnum = vertnbr + baseval;
  cblktax[cblknum].bloknum = bloknum;

  memFree (ctrbtax + baseval);                    /* Free contribution link array */

  symbptr->baseval = baseval;                     /* Fill in matrix fields */
  symbptr->cblknbr = ordeptr->cblknbr;
  symbptr->bloknbr = bloknum - baseval;
  symbptr->cblktab = cblktax + baseval;
  symbptr->bloktab = memRealloc (bloktax + baseval, (bloknum - baseval) * sizeof (SymbolBlok)); /* Set array to its exact size */
  symbptr->nodenbr = vertnbr;

#ifdef FAX_DEBUG
  if (symbolCheck (symbptr) != 0) {
    errorPrint ("symbolFax: internal error (2)");
    symbolExit (symbptr);
    return     (1);
  }
#endif /* FAX_DEBUG */

  return (0);
}
