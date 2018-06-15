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
/**   NAME       : hmesh_order_cp.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module performs compession of mesh **/
/**                nodes, by merging nodes having the same **/
/**                adjacency structure.                    **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 08 feb 2004     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.0  : from : 25 jul 2007     **/
/**                                 to   : 12 sep 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_ORDER_CP

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "order.h"
#include "mesh.h"
#include "hmesh.h"
#include "hmesh_order_cp.h"
#include "hmesh_order_st.h"

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
hmeshOrderCp (
const Hmesh * restrict const              finemeshptr,
Order * restrict const                    fineordeptr,
const Gnum                                ordenum, /*+ Zero-based ordering number +*/
OrderCblk * restrict const                cblkptr, /*+ Single column-block        +*/
const HmeshOrderCpParam * restrict const  paraptr)
{
  Hmesh                         coarmeshdat;      /* Compressed halo submesh                                               */
  Order                         coarordedat;      /* Ordering of compressed halo submesh                                   */
  Gnum * restrict               coarperitab;      /* Coarse permutation array                                              */
  const Gnum * restrict         coarperitax;      /* Temporary based access to coarperitab                                 */
  Gnum                          coarvertnum;      /* Number of current compressed vertex                                   */
  Gnum * restrict               coarvsiztax;      /* Array of coarse vertex sizes (as number of merged fine vertices)      */
  Gnum                          coaredgenum;      /* Number of current compressed edge                                     */
  Gnum * restrict               coarvpostax;      /* Position in fine permutation of fine vertices merged into same vertex */
  Gnum                          coarvnodmax;
  Gnum                          coarvertmax;
  Gnum * restrict               finecoartax;      /* Original to compressed vertex number array                            */
  HmeshOrderCpHash * restrict   finehashtab;      /* Neighbor hash table                                                   */
  Gnum                          finehashmsk;      /* Mask for access to hash table                                         */
  int * restrict                finehasptab;      /* Pre-hashing table                                                     */
  Gnum                          finehaspmsk;      /* Mask for access to pre-hashing table                                  */
  Gnum * restrict               finehsumtax;      /* Array of hash values for each original vertex                         */
  Gnum                          finevertnbr;      /* Number of fine vertices in compressed elimination tree                */
  Gnum                          finevsizsum;      /* Sum of compressed vertex sizes to build fine inverse permutation      */
  Gnum                          coarvsizsiz;      /* Size of array of sizes of coarse nodes                                */
  Gnum                          coarvelmnbr;      /* Number of coarse element vertices                                     */
  Gnum * restrict               coarverttax;
  Gnum * restrict               coarvnlotax;
  Gnum * restrict               coaredgetax;
  Gnum * restrict               coarfinetax;
  Gnum                          finevnodnum;
  Gnum                          coarvnodnbr;
  Gnum                          finevelmnum;
  Gnum                          coarvnodnnd;
  Gnum                          coardegrmax;
  Gnum                          coarvnodnum;

  if (finemeshptr->vnohnbr != finemeshptr->m.vnodnbr) {
    errorPrint ("hmeshOrderCp: halo meshes not supported yet");
    return     (1);
  }

  coarvelmnbr = finemeshptr->m.velmnbr;           /* To date, keep isolated elements */
  coarvnodmax = (Gnum) ((double) finemeshptr->vnohnbr * paraptr->comprat) +
                (finemeshptr->m.vnodnbr - finemeshptr->vnohnbr);
  coarvertmax = coarvnodmax + coarvelmnbr;
  coarvsizsiz = (finemeshptr->m.vnlotax == NULL) ? 0 : coarvnodmax;

  for (finehashmsk = 15;                          /* Set neighbor hash table size */
       finehashmsk < finemeshptr->m.degrmax;
       finehashmsk = finehashmsk * 2 + 1) ;
  finehashmsk = finehashmsk * 4 + 3;              /* Fill hash table at 1/4 of capacity */

  if (memAllocGroup ((void **) (void *)
                     &coarverttax, (size_t) ((coarvertmax + 1)      * sizeof (Gnum)),
                     &coarvsiztax, (size_t) (coarvsizsiz            * sizeof (Gnum)), /* TRICK: if no vertex loads, coarvsiztax points to coarvnodtax */
                     &coarvnlotax, (size_t) (coarvnodmax            * sizeof (Gnum)), /* Only change node weights                                     */
                     &finecoartax, (size_t) (finemeshptr->m.vnodnbr * sizeof (Gnum)),
                     &coaredgetax, (size_t) (finemeshptr->m.edgenbr * sizeof (Gnum)),
                     &finehsumtax, (size_t) (finemeshptr->m.vnodnbr * sizeof (Gnum)),
                     &finehashtab, (size_t) ((finehashmsk + 1)      * sizeof (HmeshOrderCpHash)),
                     &coarperitab, (size_t) (coarvnodmax            * sizeof (Gnum)), /* TODO: move after resize */
                     &coarfinetax, (size_t) (coarvnodmax            * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("hmeshOrderCp: out of memory (1)");
    return     (1);
  }
/* TODO : resize after success */

  finehsumtax -= finemeshptr->m.vnodbas;          /* TRICK: do not base finecoartax yet (see later) */
  finehasptab  = (int *) finecoartax;             /* Use finecoartab as temporary pre-hash table    */
  for (finehaspmsk = 1;                           /* Get pre-hash mask that fits in finecoartab     */
       finehaspmsk <= finemeshptr->m.vnodnbr;     /* Smallest (2^i)-1 value > vertnbr               */
       finehaspmsk = finehaspmsk * 2 + 1) ;
  finehaspmsk >>= 1;                              /* Ensure masked data will always fit into finecoartab array */
  finehaspmsk = (finehaspmsk * (sizeof (Gnum) / sizeof (int))) + ((sizeof (Gnum) / sizeof (int)) - 1);
  if (finehaspmsk >= ((sizeof (int) << (3 + 1)) - 1)) /* Only use 1/8 of array for pre-hashing, for increased cache locality */
    finehaspmsk >>= 3;
  memSet (finehasptab, 0, (finehaspmsk + 1) * sizeof (int)); /* Initialize pre-hash table */

  for (finevnodnum = finemeshptr->m.vnodbas, coarvnodnbr = finemeshptr->m.vnodnbr; /* For all non-halo node vertices */
       finevnodnum < finemeshptr->vnohnnd; finevnodnum ++) {
    Gnum                fineenodnum;              /* Current edge number */
    Gnum                finehsumval;              /* Hash sum value      */
    Gnum                finehsumbit;

    for (fineenodnum = finemeshptr->m.verttax[finevnodnum], finehsumval = 0;
         fineenodnum < finemeshptr->m.vendtax[finevnodnum]; fineenodnum ++)
      finehsumval += finemeshptr->m.edgetax[fineenodnum];

    finehsumtax[finevnodnum] = finehsumval;

    finehsumbit = finehsumval & ((sizeof (int) << 3) - 1); /* Get bit mask and byte position (division should be optimized into a shift) */
    finehsumval /= (sizeof (int) << 3);
    finehsumval &= finehaspmsk;                   /* Make hash sum value fit into finehasptab                                                   */
    coarvnodnbr -= (finehasptab[finehsumval] >> finehsumbit) & 1;  /* If hash value already in pre-hash table, maybe one more vertex compressed */
    finehasptab[finehsumval] |= (1 << finehsumbit); /* Put value into pre-hash table anyway                                                     */
  }

  if (coarvnodnbr > coarvnodmax) {                /* If mesh needs not be compressed */
    memFree (coarverttax);                        /* Group leader not yet based      */
    return  (hmeshOrderSt (finemeshptr, fineordeptr, ordenum, cblkptr, paraptr->stratunc));
  }

  memSet (finecoartax, ~0, finemeshptr->m.vnodnbr * sizeof (Gnum));

  memSet (&coarmeshdat, 0, sizeof (Hmesh));       /* Initialize compressed halo mesh structure   */
  coarmeshdat.m.flagval = MESHFREEVERT | MESHVERTGROUP; /* Free only coarverttab as group leader */
  coarmeshdat.m.baseval = finemeshptr->m.baseval;
  coarmeshdat.m.velmbas = coarmeshdat.m.baseval;
  coarmeshdat.m.velmnbr = coarvelmnbr;            /* Mesh compression does not touch elements, apart from isolated elements */
  coarmeshdat.m.velmnnd =
  coarmeshdat.m.vnodbas = coarvelmnbr + coarmeshdat.m.baseval;
  coarmeshdat.m.veisnbr = finemeshptr->m.veisnbr; /* To date, keep all isolated element vertices, if any */
  coarmeshdat.m.velosum = finemeshptr->m.velosum;
  coarmeshdat.m.vnlosum = finemeshptr->m.vnlosum;

  coarverttax -= coarmeshdat.m.baseval;
  coarvsiztax -= coarmeshdat.m.vnodbas;           /* TRICK: if no vertex loads, coarvsiztax points to coarvnodtax */
  coarvnlotax -= coarmeshdat.m.vnodbas;
  coaredgetax -= coarmeshdat.m.baseval;
  finecoartax -= finemeshptr->m.vnodbas;

  coarmeshdat.m.verttax = coarverttax;
  coarmeshdat.m.vendtax = coarverttax + 1;        /* Use compact representation of arrays     */
  coarmeshdat.m.velotax = finemeshptr->m.velotax; /* Re-use element vertex load array, if any */
  coarmeshdat.m.vnlotax = coarvnlotax;
  coarmeshdat.m.edgetax = coaredgetax;

  coarfinetax -= coarmeshdat.m.vnodbas;

  memSet (finehashtab, ~0, (finehashmsk + 1) * sizeof (HmeshOrderCpHash));

  for (finevelmnum = finemeshptr->m.velmbas, coarvertnum = coaredgenum = coarmeshdat.m.baseval, /* Build element arrays */
       coarvnodnnd = coarmeshdat.m.baseval + finemeshptr->m.velmnbr, coardegrmax = 0;
       finevelmnum < finemeshptr->m.velmnnd; finevelmnum ++) {
    Gnum                fineeelmnum;
    Gnum                coardegrval;

    fineeelmnum = finemeshptr->m.verttax[finevelmnum];
#ifdef DEAD_CODE
    if (fineeelmnum == finemeshptr->m.vendtax[finevelmnum]) /* Skip isolated elements */
      continue;
#endif
#ifdef SCOTCH_DEBUG_ORDER2
    if (coarvertnum >= coarmeshdat.m.velmnnd) {   /* If too many elements declared   */
      errorPrint ("hmeshOrderCp: internal error (1)"); /* Maybe problem with veisnbr */
      return     (1);
    }
#endif /* SCOTCH_DEBUG_ORDER2 */

    coarverttax[coarvertnum] = coaredgenum;

    for ( ; fineeelmnum < finemeshptr->m.vendtax[finevelmnum]; fineeelmnum ++) {
      Gnum                finevnodnum;
      Gnum                coarvnodnum;

      finevnodnum = finemeshptr->m.edgetax[fineeelmnum];
      coarvnodnum = finecoartax[finevnodnum];
      if (coarvnodnum != ~0) {                    /* If fine node already considered   */
        if (coarvnodnum >= 0)                     /* If node is leader of cluster      */
          coaredgetax[coaredgenum ++] = coarvnodnum; /* Add it to coarse mesh          */
      }
      else {                                      /* Fine node not yet considered */
        Gnum                finehsumval;
        Gnum                finedegrval;
        Gnum                fineeelmngb;
        Gnum                coarvsizval;          /* Number of fine node vertices in coarse node vertex  */

        if (coarvnodnnd >= (coarvnodmax + coarmeshdat.m.vnodbas)) { /* If mesh needs not be compressed */
          memFree (coarverttax + coarmeshdat.m.baseval);
          return  (hmeshOrderSt (finemeshptr, fineordeptr, ordenum, cblkptr, paraptr->stratunc));
        }

        coarvsizval = 1;                          /* Cluster leader it at least alone     */
        coarvnodnum = coarvnodnnd ++;             /* Node is leader of future cluster     */
        finecoartax[finevnodnum]    = coarvnodnum; /* Record node as cluster leader       */
        coaredgetax[coaredgenum ++] = coarvnodnum; /* Add leader to coarse mesh           */
        coarfinetax[coarvnodnum]    = finevnodnum; /* Record node to build edge sub-array */
        finehsumval = finehsumtax[finevnodnum];   /* Get hash sum of cluster leader       */
        finedegrval = finemeshptr->m.vendtax[finevnodnum] - finemeshptr->m.verttax[finevnodnum];

        for (fineeelmngb = fineeelmnum + 1;       /* For all remaining edges of current element */
             fineeelmngb < finemeshptr->m.vendtax[finevelmnum]; fineeelmngb ++) {
          Gnum                finevnodngb;
          Gnum                fineenodngb;

          finevnodngb = finemeshptr->m.edgetax[fineeelmngb]; /* Get index of neighboring node */

          if ((finecoartax[finevnodngb] != ~0)          || /* If node has already been processed or */
              (finehsumval != finehsumtax[finevnodngb]) || /* If hash sum does not match, skip node */
              (finedegrval != (finemeshptr->m.vendtax[finevnodngb] - finemeshptr->m.verttax[finevnodngb])))
            continue;

          if (finehashtab[(finevelmnum * HMESHORDERCPHASHPRIME) & finehashmsk].vnodnum != finevnodnum) { /* If hash table not yet filled */
            Gnum                fineenodnum;

            for (fineenodnum = finemeshptr->m.verttax[finevnodnum];
                 fineenodnum < finemeshptr->m.vendtax[finevnodnum]; fineenodnum ++) {
              Gnum                finevelmend;
              Gnum                finehelmend;

              finevelmend = finemeshptr->m.edgetax[fineenodnum];
              for (finehelmend = (finevelmend * HMESHORDERCPHASHPRIME) & finehashmsk;
                   finehashtab[finehelmend].vnodnum == finevnodnum;
                   finehelmend = (finehelmend + 1) & finehashmsk) ;
              finehashtab[finehelmend].vnodnum = finevnodnum; /* Fill hash table with node adjacency */
              finehashtab[finehelmend].velmnum = finevelmend;
            }
          }

          for (fineenodngb = finemeshptr->m.verttax[finevnodngb];
               fineenodngb < finemeshptr->m.vendtax[finevnodngb]; fineenodngb ++) {
            Gnum                finevelmngb;
            Gnum                finehelmngb;

            finevelmngb = finemeshptr->m.edgetax[fineenodngb];
            for (finehelmngb = (finevelmngb * HMESHORDERCPHASHPRIME) & finehashmsk; ; finehelmngb = (finehelmngb + 1) & finehashmsk) {
              if (finehashtab[finehelmngb].vnodnum != finevnodnum) /* If adjacencies differ, break */
                goto loop_failed;
              if (finehashtab[finehelmngb].velmnum == finevelmngb) /* If neighbor found, process next neighbor */
                break;
            }
          }
          finecoartax[finevnodngb] = -2 - coarvnodnum; /* Set index of cluster non-leader */
          coarvsizval ++;                         /* One more non-leader in cluster       */
loop_failed: ;
        }
        coarvsiztax[coarvnodnum] = coarvsizval;
      }
    }

    coardegrval = coaredgenum - coarverttax[coarvertnum];
    if (coardegrval > coardegrmax)
      coardegrmax = coardegrval;

    coarvertnum ++;                               /* One more coarse element created */
  }
#ifdef SCOTCH_DEBUG_ORDER2
  if (coarvertnum != coarmeshdat.m.velmnnd) {     /* If too many elements declared */
    errorPrint ("hmeshOrderCp: internal error (2)"); /* Maybe problem with veisnbr */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ORDER2 */
  coarmeshdat.m.vnodnnd = coarvnodnnd;
  coarmeshdat.m.vnodnbr = coarvnodnnd - coarmeshdat.m.vnodbas;

#ifdef SCOTCH_DEBUG_ORDER2
  for (finevnodnum = finemeshptr->m.vnodbas;
       finevnodnum < finemeshptr->m.vnodnnd; finevnodnum ++) {
    if (finecoartax[finevnodnum] == ~0) {
      errorPrint ("hmeshOrderCp: internal error (3)");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  for ( ; coarvertnum < coarmeshdat.m.vnodnnd; coarvertnum ++) { /* Build node arrays */
    Gnum                finevnodnum;
    Gnum                fineenodnum;
    Gnum                coardegrval;

    coarverttax[coarvertnum] = coaredgenum;

    finevnodnum = coarfinetax[coarvertnum];
    coardegrval = finemeshptr->m.vendtax[finevnodnum] - finemeshptr->m.verttax[finevnodnum];
    for (fineenodnum = finemeshptr->m.verttax[finevnodnum];
         fineenodnum < finemeshptr->m.vendtax[finevnodnum]; fineenodnum ++)
      coaredgetax[coaredgenum ++] = finemeshptr->m.edgetax[fineenodnum] - (finemeshptr->m.velmbas - coarmeshdat.m.velmbas);

    if (coardegrval > coardegrmax)
      coardegrmax = coardegrval;
  }
  coarverttax[coarvertnum] = coaredgenum;         /* Set end of vertex array */

  coarmeshdat.m.edgenbr = coaredgenum - coarmeshdat.m.baseval;
  coarmeshdat.m.degrmax = coardegrmax;
  coarmeshdat.vnohnbr   = coarmeshdat.m.vnodnbr;  /* Halo meshes not yet supported */
  coarmeshdat.vnohnnd   = coarmeshdat.m.vnodnnd;
  coarmeshdat.vehdtax   = coarmeshdat.m.vendtax;  /* Only element part of vendtab will be accessed through vehdtab */
  coarmeshdat.vnhlsum   = coarmeshdat.m.vnlosum;
  coarmeshdat.enohnbr   = coarmeshdat.m.edgenbr;

  if (finemeshptr->m.vnlotax != NULL) {           /* If fine mesh has node vertex loads */
    memSet (coarmeshdat.m.vnlotax + coarmeshdat.m.vnodbas, 0, coarmeshdat.m.vnodnbr * sizeof (Gnum));

    for (finevnodnum = finemeshptr->m.vnodbas; finevnodnum < finemeshptr->m.vnodnnd; finevnodnum ++) { /* Compute vertex loads for compressed mesh */
      coarmeshdat.m.vnlotax[finecoartax[finevnodnum]] += finemeshptr->m.vnlotax[finevnodnum];
    }
  }

#ifdef SCOTCH_DEBUG_ORDER2
  if (hmeshCheck (&coarmeshdat) != 0) {
    errorPrint ("hmeshOrderCp: internal error (4)");
    hmeshExit  (&coarmeshdat);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ORDER2 */

  orderInit (&coarordedat, coarmeshdat.m.baseval, coarmeshdat.m.vnodnbr, coarperitab); /* Build ordering of compressed submesh */
  if (hmeshOrderSt (&coarmeshdat, &coarordedat, 0, &coarordedat.cblktre, paraptr->stratcpr) != 0) {
    hmeshExit (&coarmeshdat);
    return    (1);
  }

  coarvsiztax += (coarmeshdat.m.vnodbas - coarmeshdat.m.baseval); /* Adjust array to match permutation bounds */

  *cblkptr = coarordedat.cblktre;                 /* Link sub-tree to ordering         */
  coarordedat.cblktre.cblktab = NULL;             /* Unlink sub-tree from sub-ordering */
  finevertnbr = hmeshOrderCpTree (coarordedat.peritab, /* Expand sub-tree              */
                                  coarvsiztax, cblkptr, 0);
#ifdef SCOTCH_DEBUG_ORDER2
  if (finevertnbr != finemeshptr->m.vnodnbr) {
    errorPrint ("hmeshOrderCp: internal error (5)");
    hmeshExit  (&coarmeshdat);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ORDER2 */
  fineordeptr->treenbr += coarordedat.cblknbr;    /* Adjust number of tree nodes     */
  fineordeptr->cblknbr += coarordedat.cblknbr - 1; /* Adjust number of column blocks */

  coarvpostax = coarmeshdat.m.verttax;            /* Recycle verttab (not velotab as may be merged with coarvsiztab) */
  coarperitax = coarperitab - coarmeshdat.m.vnodbas;

  for (coarvnodnum = coarmeshdat.m.vnodbas, finevsizsum = 0; /* Compute initial indices for inverse permutation expansion */
       coarvnodnum < coarmeshdat.m.vnodnnd; coarvnodnum ++) {
    coarvpostax[coarperitax[coarvnodnum]] = finevsizsum;
    finevsizsum += coarvsiztax[coarperitax[coarvnodnum]];
  }
  coarvpostax = coarmeshdat.m.verttax + (coarmeshdat.m.baseval - coarmeshdat.m.vnodbas);
  for (finevnodnum = finemeshptr->m.vnodbas; finevnodnum < finemeshptr->m.vnodnnd; finevnodnum ++) { /* Compute fine permutation */
    Gnum                coarvnodnum;

    coarvnodnum = finecoartax[finevnodnum];       /* Get index of corresponding coarse node */
    if (coarvnodnum < 0)                          /* If node is not cluster leader          */
      coarvnodnum = -2 - coarvnodnum;             /* Get index of cluster leader            */
    fineordeptr->peritab[coarvpostax[coarvnodnum] ++] = finevnodnum + (finemeshptr->m.baseval - finemeshptr->m.vnodbas);
  }

  orderExit (&coarordedat);
  hmeshExit (&coarmeshdat);

  return (0);
}

/* This routine turns the coarse elimination
** tree produced by the ordering of the coarse
** mesh into a fine elimination tree, according
** to the cardinality of the coarse vertices.
** It returns:
** - !0  : overall number of fine vertices, in all cases.
*/

static
Gnum
hmeshOrderCpTree (
const Gnum * restrict const coarperitab,          /* Coarse inverse permutation              */
const Gnum * restrict const coarvsiztax,          /* Array of fine sizes of coarse vertices  */
OrderCblk * restrict const  coficblkptr,          /* Current coarse/fine column block cell   */
Gnum                        coarordenum)          /* Compressed vertex to start expansion at */
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
      finevertnbr += hmeshOrderCpTree (coarperitab, coarvsiztax, &coficblkptr->cblktab[coficblknum], coarvertsum);
      coarvertsum += coarvertnbr;                 /* Sum-up coarse vertices */
    }
  }
  coficblkptr->vnodnbr = finevertnbr;             /* Set number of fine vertices */

  return (finevertnbr);                           /* Return accumulated number */
}
