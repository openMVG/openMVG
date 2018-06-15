/* Copyright 2004,2007,2009,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_cp.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders vertices by compres- **/
/**                sing vertices with identical adjacency  **/
/**                structure.                              **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 29 aug 1998     **/
/**                                 to     12 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     03 jan 1999     **/
/**                # Version 4.0  : from : 01 jan 2003     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.0  : from : 29 dec 2006     **/
/**                                 to     22 may 2008     **/
/**                # Version 5.1  : from : 01 oct 2009     **/
/**                                 to   : 01 oct 2009     **/
/**                # Version 6.0  : from : 04 aug 2014     **/
/**                                 to   : 24 feb 2015     **/
/**                                                        **/
/**   NOTES      : # Pre-hashing proves itself extremely   **/
/**                  efficient, since for graphs that      **/
/**                  will be compressed very few writes    **/
/**                  will be performed in the pre-hashing  **/
/**                  array, and for others, for which pre- **/
/**                  hashing costs much more, it will save **/
/**                  time in the end.                      **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_CP

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_cp.h"
#include "hgraph_order_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the ordering.
** It returns:
** - 0   : if the ordering could be computed.
** - !0  : on error.
*/

int
hgraphOrderCp (
const Hgraph * restrict const     finegrafptr,
Order * restrict const            fineordeptr,
const Gnum                        ordenum,        /*+ Zero-based ordering number +*/
OrderCblk * restrict const        cblkptr,        /*+ Single column-block        +*/
const HgraphOrderCpParam * const  paraptr)
{
  Hgraph                        coargrafdat;      /* Compressed halo subgraph                                              */
  Order                         coarordedat;      /* Ordering of compressed halo subgraph                                  */
  Gnum *                        coarperitab;      /* Coarse permutation array                                              */
  const Gnum * restrict         coarperitax;      /* Temporary based access to coarperitab                                 */
  Gnum                          coarvertnbr;      /* Number of compressed vertices                                         */
  Gnum                          coarvertnum;      /* Number of current compressed vertex                                   */
  Gnum                          coarvnhdsiz;      /* Size of non-halo end vertex array; zero if graph has no halo          */
  Gnum                          coarvsizsiz;      /* Size of coarse vertex sizes array; zero if no fine vertex loads       */
  Gnum * restrict               coarvsiztax;      /* Array of coarse vertex sizes (as number of merged fine vertices)      */
  Gnum                          coaredgenbr;      /* Number of compressed edges                                            */
  Gnum                          coaredgenum;      /* Number of current compressed edge                                     */
  Gnum                          coarenohnnd;      /* Position in edge array of first edge of first halo vertex             */
  Gnum * restrict               coarvpostax;      /* Position in fine permutation of fine vertices merged into same vertex */
  Gnum * restrict               finecoartax;      /* Original to compressed vertex number array                            */
  HgraphOrderCpMate * restrict  finematetab;      /* Array of fine vertices that may be compressed with current vertex     */
  HgraphOrderCpHash * restrict  finehashtab;      /* Neighbor hash table                                                   */
  Gnum                          finehashmsk;      /* Mask for access to hash table                                         */
  int * restrict                finehasptab;      /* Pre-hashing table                                                     */
  Gnum                          finehaspmsk;      /* Mask for access to pre-hashing table                                  */
  Gnum * restrict               finehsumtax;      /* Array of hash values for each original vertex                         */
  Gnum                          finevertnbr;      /* Number of fine vertices in compressed elimination tree                */
  Gnum                          finevertnum;      /* Number of current original vertex                                     */
  Gnum                          finevsizsum;      /* Sum of compressed vertex sizes to build fine inverse permutation      */
  void *                        dataptr;          /* Flag of memory allocation success                                     */

  Gnum * restrict const         fineperitab = fineordeptr->peritab;
  const Gnum * restrict const   fineverttax = finegrafptr->s.verttax;
  const Gnum * restrict const   finevelotax = finegrafptr->s.velotax;
  const Gnum * restrict const   finevendtax = finegrafptr->s.vendtax;
  const Gnum * restrict const   finevnumtax = finegrafptr->s.vnumtax;
  const Gnum * restrict const   finevnhdtax = finegrafptr->vnhdtax;
  const Gnum * restrict const   fineedgetax = finegrafptr->s.edgetax;

  for (finehashmsk = 15;                          /* Set neighbor hash table sizes */
       finehashmsk < finegrafptr->s.degrmax;
       finehashmsk = finehashmsk * 2 + 1) ;
  finehashmsk = finehashmsk * 4 + 3;              /* Fill hash table at 1/4 of capacity */

  if (((finecoartax = (Gnum *) memAlloc (finegrafptr->s.vertnbr * sizeof (Gnum))) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &finehashtab, (size_t) ((finehashmsk + 1)      * sizeof (HgraphOrderCpHash)),
                      &finematetab, (size_t) (finegrafptr->s.degrmax * sizeof (HgraphOrderCpMate)), NULL) == NULL) ||
      ((finehsumtax = (Gnum *) memAlloc (finegrafptr->vnohnbr * sizeof (Gnum))) == NULL)) {
    errorPrint ("hgraphOrderCp: out of memory (1)");
    if (finecoartax != NULL) {
      if (finehashtab != NULL)
        memFree (finehashtab);
      memFree (finecoartax);
    }
    return (1);
  }
  finehsumtax -= finegrafptr->s.baseval;          /* TRICK: do not base finecoartax yet (see later) */

  finehasptab = (int *) finecoartax;              /* Use finecoartab as temporary pre-hash table */
  for (finehaspmsk = 1;                           /* Get pre-hash mask that fits in finecoartab  */
       finehaspmsk < finegrafptr->s.vertnbr;      /* Smallest (2^i)-1 value >= vertnbr           */
       finehaspmsk = finehaspmsk * 2 + 1) ;
  finehaspmsk >>= 1;                              /* Ensure masked data will always fit into finecoartab array */
  finehaspmsk = (finehaspmsk * (sizeof (Gnum) / sizeof (int))) + ((sizeof (Gnum) / sizeof (int)) - 1);
  if (finehaspmsk >= ((sizeof (int) << (3 + 1)) - 1)) /* Only use 1/8 of array for pre-hashing, for increased cache locality */
    finehaspmsk >>= 3;
  memSet (finehasptab, 0, (finehaspmsk + 1) * sizeof (int)); /* Initialize pre-hash table */

  for (finevertnum = finegrafptr->s.baseval, coarvertnbr = finegrafptr->vnohnbr; /* For all non-halo vertices */
       finevertnum < finegrafptr->vnohnnd; finevertnum ++) {
    Gnum                fineedgenum;              /* Current edge number */
    Gnum                finehsumval;              /* Hash sum value      */
    Gnum                finehsumbit;

    for (fineedgenum = fineverttax[finevertnum], finehsumval = finevertnum; /* For all edges, including halo edges */
         fineedgenum < finevendtax[finevertnum]; fineedgenum ++)
      finehsumval += fineedgetax[fineedgenum];

    finehsumtax[finevertnum] = finehsumval;

    finehsumbit = finehsumval & ((sizeof (int) << 3) - 1); /* Get bit mask and byte position (division should be optimized into a shift) */
    finehsumval /= (sizeof (int) << 3);
    finehsumval &= finehaspmsk;                   /* Make hash sum value fit into finehasptab                                                   */
    coarvertnbr -= (finehasptab[finehsumval] >> finehsumbit) & 1;  /* If hash value already in pre-hash table, maybe one more vertex compressed */
    finehasptab[finehsumval] |= (1 << finehsumbit); /* Put value into pre-hash table anyway                                                     */
  }

  if ((double) coarvertnbr > ((double) finegrafptr->vnohnbr * paraptr->comprat)) { /* If graph needs not be compressed */
    memFree (finehsumtax + finegrafptr->s.baseval);
    memFree (finehashtab);
    memFree (finecoartax);                        /* Not yet based */
    return (hgraphOrderSt (finegrafptr, fineordeptr, ordenum, cblkptr, paraptr->stratunc));
  }

  finecoartax -= finegrafptr->s.baseval;          /* Base finecoartab array */

  memSet (finehashtab, ~0, (finehashmsk + 1) * sizeof (HgraphOrderCpHash));

  hgraphInit (&coargrafdat);                      /* Initialize compressed halo graph structure                               */
  coargrafdat.s.baseval = 1;                      /* Base coarse graph to 1 because hgraphOrderHb and hgraphOrderHf prefer it */

  for (finevertnum = finegrafptr->s.baseval, coarvertnbr = coargrafdat.s.baseval, coaredgenbr = finegrafptr->s.edgenbr; /* For all non-halo vertices */
       finevertnum < finegrafptr->vnohnnd; finevertnum ++) {
    Gnum                finedegrval;              /* Degree of current fine vertex     */
    Gnum                finehsumval;              /* Current hash sum value            */
    Gnum                finematenbr;              /* Number of mates of current vertex */
    Gnum                fineedgenum;              /* Current edge number               */

    finedegrval = finevendtax[finevertnum] - fineverttax[finevertnum];
    finehsumval = finehsumtax[finevertnum];
    finematenbr = 0;                              /* Reset potential mate array */

    for (fineedgenum = fineverttax[finevertnum];  /* For all edges, including halo edges */
         fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
      Gnum                finevertend;

      finevertend = fineedgetax[fineedgenum];

      if ((finevertend < finevertnum) &&          /* If neighbor has same characteristics */
          (finehsumval == finehsumtax[finevertend]) &&
          (finedegrval == (finevendtax[finevertend] - fineverttax[finevertend]))) {
        Gnum                finematenum;
        Gnum                coarvertend;

        for (finematenum = 0, coarvertend = finecoartax[finevertend]; /* Search if end vertex has already been compressed with some mate */
             (finematenum < finematenbr) && (finematetab[finematenum].coarvertend != coarvertend); finematenum ++) ;

        if (finematenum == finematenbr) {         /* If new slot needed   */
          finematetab[finematenum].coarvertend = coarvertend; /* Build it */
          finematetab[finematenum].finevertend = finevertend;
          finematenbr ++;
        }
      }
    }

    finecoartax[finevertnum] = coarvertnbr ++;    /* Assume no mate found */

    if (finematenbr > 0) {                        /* If potential mates exist */
      Gnum                fineedgenum;            /* Current edge number      */
      Gnum                finehashnum;

      for (fineedgenum = fineverttax[finevertnum]; /* For all edges, including halo edges */
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        Gnum                finevertend;

        finevertend = fineedgetax[fineedgenum];   /* Add end vertex to hash table */

        for (finehashnum = (finevertend * HGRAPHORDERCPHASHPRIME) & finehashmsk; /* Search for empty slot in hash table */
             finehashtab[finehashnum].vertnum == finevertnum; finehashnum = (finehashnum + 1) & finehashmsk) ;
        finehashtab[finehashnum].vertnum = finevertnum;
        finehashtab[finehashnum].vertend = finevertend;
      }
      for (finehashnum = (finevertnum * HGRAPHORDERCPHASHPRIME) & finehashmsk;  /* Add current vertex to hash table */
           finehashtab[finehashnum].vertnum == finevertnum; finehashnum = (finehashnum + 1) & finehashmsk) ;
      finehashtab[finehashnum].vertnum = finevertnum;
      finehashtab[finehashnum].vertend = finevertnum;

      finematenbr --;                             /* Point to first potential mate */
      do {                                        /* For all potential mates       */
        Gnum                fineedgenum;          /* Current edge number           */
        Gnum                fineedgennd;

        for (fineedgenum = fineverttax[finematetab[finematenbr].finevertend], /* For all edges, including halo edges */
             fineedgennd = finevendtax[finematetab[finematenbr].finevertend];
             fineedgenum < fineedgennd; fineedgenum ++) {
          Gnum                finevertend;

          finevertend = fineedgetax[fineedgenum];

          for (finehashnum = (finevertend * HGRAPHORDERCPHASHPRIME) & finehashmsk; ;
               finehashnum = (finehashnum + 1) & finehashmsk) {
            if (finehashtab[finehashnum].vertnum != finevertnum) /* If mate neighbor not found in hash table  */
              goto loop_failed;                   /* Vertex cannot be merged to mate, so skip to next mate    */
            if (finehashtab[finehashnum].vertend == finevertend) /* Else if mate neighbor found in hash table */
              break;                              /* Skip to next mate neighbor to find                       */
          }
        }
        coarvertnbr --;                           /* Same adjacency structure          */
        finecoartax[finevertnum] = finematetab[finematenbr].coarvertend; /* Get number */
        coaredgenbr -= finedegrval + 1;           /* Remove exceeding edges            */
        break;
loop_failed: ;
      } while (finematenbr -- > 0);
    }
  }

  coargrafdat.vnohnnd = coarvertnbr;              /* Save number of non-halo vertices */

  memFree (finehsumtax + finegrafptr->s.baseval);

  if ((double) (coarvertnbr - coargrafdat.s.baseval) > ((double) finegrafptr->vnohnbr * paraptr->comprat)) { /* If graph needs not be compressed */
    memFree (finehashtab);
    memFree (finecoartax + finegrafptr->s.baseval);
    return (hgraphOrderSt (finegrafptr, fineordeptr, ordenum, cblkptr, paraptr->stratunc));
  }

  for ( ; finevertnum < finegrafptr->s.vertnnd; finevertnum ++) /* For all halo vertices */
    finecoartax[finevertnum] = coarvertnbr ++;    /* Halo vertices are never compressed  */

  coargrafdat.s.flagval = GRAPHFREETABS | GRAPHVERTGROUP; /* Do not set HGRAPHFREEVNHD since vnhdtax allocated in group */
  coargrafdat.s.vertnbr = coarvertnbr - coargrafdat.s.baseval;
  coargrafdat.s.vertnnd = coarvertnbr;
  coargrafdat.s.velosum = finegrafptr->s.velosum;
  coargrafdat.s.degrmax = finegrafptr->s.degrmax;
  coargrafdat.vnohnbr   = coargrafdat.vnohnnd - coargrafdat.s.baseval;
  coargrafdat.vnlosum   = finegrafptr->vnlosum;

  coarvnhdsiz = (finegrafptr->s.vertnbr == finegrafptr->vnohnbr) ? 0 : coargrafdat.vnohnbr; /* If no halo, no need for vnhdtax; will use vendtax */
  coarvsizsiz = (finevelotax == NULL) ? 0 : coarvertnbr; /* If no fine vertex loads, use coarse velotax as coarvsiztax                           */
  if ((dataptr = memAllocGroup ((void **) (void *)
                                &coargrafdat.s.verttax, (size_t) ((coarvertnbr + 1) * sizeof (Gnum)),
                                &coargrafdat.vnhdtax,   (size_t) (coarvnhdsiz       * sizeof (Gnum)),
                                &coargrafdat.s.velotax, (size_t) (coarvertnbr       * sizeof (Gnum)),
                                &coarvsiztax,           (size_t) (coarvsizsiz       * sizeof (Gnum)), NULL)) != NULL) {
    dataptr               =
    coargrafdat.s.edgetax = (Gnum *) memAlloc (coaredgenbr * sizeof (Gnum));
  }
  if (dataptr == NULL) {
    errorPrint ("hgraphOrderCp: out of memory (2)");
    hgraphExit (&coargrafdat);
    memFree    (finehashtab);
    memFree    (finecoartax + finegrafptr->s.baseval);
    return     (1);
  }
  coargrafdat.s.verttax -= coargrafdat.s.baseval;
  coargrafdat.s.vendtax  = coargrafdat.s.verttax + 1; /* Use compact representation of arrays */
  coargrafdat.s.velotax -= coargrafdat.s.baseval;
  coargrafdat.s.edgetax -= coargrafdat.s.baseval;
  coargrafdat.vnhdtax    = (finegrafptr->s.vertnbr == finegrafptr->vnohnbr) ? coargrafdat.s.vendtax : coargrafdat.vnhdtax - coargrafdat.s.baseval;
  coarvsiztax            = (finevelotax == NULL) ? coargrafdat.s.velotax : coarvsiztax - coargrafdat.s.baseval;

  memSet (finehashtab, ~0, (finehashmsk + 1) * sizeof (HgraphOrderCpHash));

  for (finevertnum = finegrafptr->s.baseval, coarvertnum = coaredgenum = coargrafdat.s.baseval; /* For all non-halo vertices */
       finevertnum < finegrafptr->vnohnnd; finevertnum ++) {
    Gnum                fineedgenum;              /* Current edge number */

    if (finecoartax[finevertnum] != coarvertnum)  /* Skip fine vertices until we find one that is part of current coarse vertex */
      continue;

    coargrafdat.s.verttax[coarvertnum] = coaredgenum;
    coarvsiztax[coarvertnum] = 1;                 /* Fill coargrafdat.s.velotax if finegrafptr has no vertex loads */

    for (fineedgenum = fineverttax[finevertnum];  /* For all non-halo edges of vertex */
         fineedgenum < finevnhdtax[finevertnum]; fineedgenum ++) {
      Gnum                finevertend;
      Gnum                finehashnum;

      finevertend = fineedgetax[fineedgenum];
      if (finecoartax[finevertend] == coarvertnum) { /* If neighbor is merged into us, merge load but do not write edge */
        coarvsiztax[coarvertnum] ++;              /* Fill coargrafdat.s.velotax if finegrafptr has no vertex loads      */
        continue;
      }
      for (finehashnum = (finecoartax[finevertend] * HGRAPHORDERCPHASHPRIME) & finehashmsk; ; /* Search for end vertex in hash table */
           finehashnum = (finehashnum + 1) & finehashmsk) {
        if (finehashtab[finehashnum].vertnum != coarvertnum) {
          finehashtab[finehashnum].vertnum = coarvertnum;
          finehashtab[finehashnum].vertend =
          coargrafdat.s.edgetax[coaredgenum ++] = finecoartax[finevertend];
          break;
        }
        if (finehashtab[finehashnum].vertend == finecoartax[finevertend])
          break;                                  /* If edge already exists */
      }
    }
    coargrafdat.vnhdtax[coarvertnum] = coaredgenum; /* Set end of non-halo edge sub-array */

    for ( ; fineedgenum < finegrafptr->s.vendtax[finevertnum]; fineedgenum ++) { /* For edges linking to halo vertices */
      Gnum                finevertend;

      finevertend = fineedgetax[fineedgenum];
      coargrafdat.s.edgetax[coaredgenum ++] = finecoartax[finevertend]; /* Halo vertices are always defined and unique */
    }
    coarvertnum ++;
  }
  for (coarenohnnd = coaredgenum; finevertnum < finegrafptr->s.vertnnd; finevertnum ++) { /* For all halo vertices */
    Gnum                fineedgenum;              /* Current edge number */

#ifdef SCOTCH_DEBUG_ORDER2
    if (finecoartax[finevertnum] != coarvertnum) {
      errorPrint ("hgraphOrderCp: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_ORDER2 */

    coargrafdat.s.verttax[coarvertnum] = coaredgenum;
    coarvsiztax[coarvertnum] = 1;                 /* Fill coargrafdat.s.velotax if finegrafptr has no vertex loads */

    for (fineedgenum = fineverttax[finevertnum];  /* For all edges of halo vertex */
         fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
      Gnum                finevertend;
      Gnum                finehashnum;

      finevertend = fineedgetax[fineedgenum];
#ifdef SCOTCH_DEBUG_ORDER2
      if (finecoartax[finevertend] == coarvertnum) { /* No neighbor can be merged into us since halo vertices are unique */
        errorPrint ("hgraphOrderCp: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_ORDER2 */
      for (finehashnum = (finecoartax[finevertend] * HGRAPHORDERCPHASHPRIME) & finehashmsk; ; /* Search for end vertex in hash table */
           finehashnum = (finehashnum + 1) & finehashmsk) {
        if (finehashtab[finehashnum].vertnum != coarvertnum) {
          finehashtab[finehashnum].vertnum = coarvertnum;
          finehashtab[finehashnum].vertend =
          coargrafdat.s.edgetax[coaredgenum ++] = finecoartax[finevertend];
          break;
        }
        if (finehashtab[finehashnum].vertend == finecoartax[finevertend])
          break;                                  /* If edge already exists */
      }
    }
    coarvertnum ++;
  }
  coargrafdat.s.verttax[coarvertnum] = coaredgenum; /* Set end of compact vertex array */
  coargrafdat.s.edlosum =
  coargrafdat.s.edgenbr = coaredgenum - coargrafdat.s.baseval;
  coargrafdat.enohsum   =
  coargrafdat.enohnbr   = coargrafdat.s.edgenbr - 2 * (coaredgenum - coarenohnnd);

  if (finevelotax != NULL) {                      /* If fine graph has vertex loads */
    memSet (coargrafdat.s.velotax + coargrafdat.s.baseval, 0, coargrafdat.s.vertnbr * sizeof (Gnum));

    for (finevertnum = finegrafptr->s.baseval; finevertnum < finegrafptr->s.vertnnd; finevertnum ++) /* Compute vertex loads for compressed graph */
      coargrafdat.s.velotax[finecoartax[finevertnum]] += finevelotax[finevertnum];
  }

  memFree (finehashtab);

  coargrafdat.s.edgetax = (Gnum *) memRealloc (coargrafdat.s.edgetax + coargrafdat.s.baseval, coargrafdat.s.edgenbr * sizeof (Gnum)) - coargrafdat.s.baseval;

#ifdef SCOTCH_DEBUG_ORDER2
  if (hgraphCheck (&coargrafdat) != 0) {
    errorPrint ("hgraphOrderCp: internal error (3)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  if ((coarperitab = memAlloc (coargrafdat.vnohnbr * sizeof (Gnum))) == NULL) { /* Coarse permutation only for non-halo vertices */
    errorPrint ("hgraphOrderCp: out of memory (3)");
    hgraphExit (&coargrafdat);
    memFree    (finecoartax + finegrafptr->s.baseval);
    return     (1);
  }
  orderInit (&coarordedat, coargrafdat.s.baseval, coargrafdat.vnohnbr, coarperitab); /* Build ordering of compressed subgraph */
  if (hgraphOrderSt (&coargrafdat, &coarordedat, 0, &coarordedat.cblktre, paraptr->stratcpr) != 0) {
    memFree    (coarperitab);
    hgraphExit (&coargrafdat);
    memFree    (finecoartax + finegrafptr->s.baseval);
    return     (1);
  }

  *cblkptr = coarordedat.cblktre;                 /* Link sub-tree to ordering         */
  coarordedat.cblktre.cblktab = NULL;             /* Unlink sub-tree from sub-ordering */
  finevertnbr = hgraphOrderCpTree (coarordedat.peritab, /* Expand sub-tree             */
                                   coarvsiztax, cblkptr, 0);
#ifdef SCOTCH_DEBUG_ORDER2
  if (finevertnbr != finegrafptr->vnohnbr) {
    errorPrint ("hgraphOrderCp: internal error (4)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ORDER2 */
  fineordeptr->treenbr += coarordedat.treenbr - 1; /* Adjust number of tree nodes    */
  fineordeptr->cblknbr += coarordedat.cblknbr - 1; /* Adjust number of column blocks */

  coarvpostax = coargrafdat.s.verttax;            /* Re-cycle verttab (not velotab as may be merged with coarvsiztab) */
  coarperitax = coarperitab - coargrafdat.s.baseval;

  for (coarvertnum = coargrafdat.s.baseval, finevsizsum = 0; /* Compute initial indices for inverse permutation expansion */
       coarvertnum < coargrafdat.vnohnnd; coarvertnum ++) {
    coarvpostax[coarperitax[coarvertnum]] = finevsizsum;
    finevsizsum += coarvsiztax[coarperitax[coarvertnum]];
  }
  if (finevnumtax == NULL) {                      /* If fine graph is original graph */
    for (finevertnum = finegrafptr->s.baseval; finevertnum < finegrafptr->vnohnnd; finevertnum ++) /* Compute fine permutation */
      fineperitab[coarvpostax[finecoartax[finevertnum]] ++] = finevertnum;
  }
  else {                                          /* Graph is not original graph */
    for (finevertnum = finegrafptr->s.baseval; finevertnum < finegrafptr->vnohnnd; finevertnum ++) /* Compute fine permutation */
      fineperitab[coarvpostax[finecoartax[finevertnum]] ++] = finevnumtax[finevertnum];
  }

  memFree    (coarperitab);
  memFree    (finecoartax + finegrafptr->s.baseval);
  orderExit  (&coarordedat);
  hgraphExit (&coargrafdat);                      /* Free coarvsiztab as part of vertex group */

  return (0);
}

/* This routine turns the coarse elimination
** tree produced by the ordering of the coarse
** graph into a fine elimination tree, according
** to the cardinality of the coarse vertices.
** It returns:
** - !0  : overall number of fine vertices, in all cases.
*/

static
Gnum
hgraphOrderCpTree (
const Gnum * restrict const coarperitab,          /* Coarse inverse permutation              */
const Gnum * restrict const coarvsiztax,          /* Array of fine sizes of coarse vertices  */
OrderCblk * restrict const  coficblkptr,          /* Current coarse/fine column block cell   */
const Gnum                  coarordenum)          /* Compressed vertex to start expansion at */
{
  Gnum                finevertnbr;                /* Number of fine vertices in subtree */

  finevertnbr = 0;                                /* No fine vertices yet */

  if (coficblkptr->cblktab == NULL) {             /* If leaf of column block tree */
    Gnum                coarvnumnum;

    for (coarvnumnum = coarordenum;
         coarvnumnum < coarordenum + coficblkptr->vnodnbr; coarvnumnum ++)
      finevertnbr += coarvsiztax[coarperitab[coarvnumnum]];   /* Sum-up fine vertices */
  }
  else {
    Gnum                coarvertnbr;              /* Number of coarse vertices in cell    */
    Gnum                coarvertsum;              /* Number of coarse vertices in subtree */
    Gnum                coficblknum;              /* Index in column block array          */

    for (coficblknum = 0, coarvertsum = coarordenum; /* Start at current coarse index */
         coficblknum < coficblkptr->cblknbr; coficblknum ++) {
      coarvertnbr  = coficblkptr->cblktab[coficblknum].vnodnbr; /* Save number of coarse vertices */
      finevertnbr += hgraphOrderCpTree (coarperitab, coarvsiztax, &coficblkptr->cblktab[coficblknum], coarvertsum);
      coarvertsum += coarvertnbr;                 /* Sum-up coarse vertices */
    }
  }
  coficblkptr->vnodnbr = finevertnbr;             /* Set number of fine vertices */

  return (finevertnbr);                           /* Return accumulated number */
}
