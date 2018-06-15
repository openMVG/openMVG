/* Copyright 2007-2010,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vdgraph_separate_ml.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : This module contains the multi-level    **/
/**                separation strategy.                    **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 07 mar 2006     **/
/**                                 to   : 01 mar 2008     **/
/**                # Version 5.1  : from : 14 dec 2008     **/
/**                                 to   : 26 aug 2010     **/
/**                # Version 6.0  : from : 11 sep 2012     **/
/**                                 to   : 28 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_SEPARATE_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "vdgraph.h"
#include "vdgraph_separate_ml.h"
#include "vdgraph_separate_st.h"

/*********************************************/
/*                                           */
/* The coarsening and uncoarsening routines. */
/*                                           */
/*********************************************/

/* This routine builds a coarser graph from the
** Dgraph that is given on input. The coarser
** Dgraphs differ at this stage from classical
** active Dgraphs as their internal gains are not
** yet computed.
** It returns:
** - 0  : if the coarse Dgraph has been built.
** - 1  : if threshold achieved or on error.
*/

static
int
vdgraphSeparateMlCoarsen (
Vdgraph * restrict const              finegrafptr, /*+ Finer graph                               +*/
Vdgraph * restrict const              coargrafptr, /*+ Coarser graph to build                    +*/
DgraphCoarsenMulti * restrict * const coarmultptr, /*+ Pointer to based multinode table to build +*/
const VdgraphSeparateMlParam * const  paraptr)     /*+ Method parameters                         +*/
{
  int                 foldval;

  switch (paraptr->foldval) {
    case 0 :
      foldval = DGRAPHCOARSENNONE;
      break;
    case 1 :
      foldval = DGRAPHCOARSENFOLD;
      break;
    case 2 :
      foldval = DGRAPHCOARSENFOLDDUP;
      break;
#ifdef SCOTCH_DEBUG_VDGRAPH2
    default :
      errorPrint ("vdgraphSeparateMlCoarsen: invalid parameter");
      return     (1);
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
  }
  if ((finegrafptr->s.vertglbnbr / finegrafptr->s.procglbnbr) > paraptr->foldmax) /* If no need to fold */
    foldval = DGRAPHCOARSENNONE;

  *coarmultptr = NULL;                            /* Let the routine create the multinode array */
  dgraphInit (&coargrafptr->s, finegrafptr->s.proccomm); /* Re-use fine graph communicator      */
  if (dgraphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr, paraptr->passnbr,
                     paraptr->coarnbr, paraptr->coarrat, foldval) != 0)
    return (1);                                   /* Return if coarsening failed */

  coargrafptr->fronloctab = NULL;
  coargrafptr->partgsttax = NULL;                 /* Do not allocate partition data yet */

  if (coargrafptr->s.procglbnbr == 0) {           /* Not a owner graph                        */
    coargrafptr->s.vertlocnbr = 0;                /* Set it to zero for vrcvdattab allocation */
    return (0);
  }

  coargrafptr->levlnum = finegrafptr->levlnum + 1; /* Graph level is coarsening level                 */
  if (coargrafptr->s.vertlocnbr <= finegrafptr->s.vertlocnbr) /* If (folded) coarser graph is smaller */
    coargrafptr->fronloctab = finegrafptr->fronloctab; /* Re-use frontier array for coarser graph     */
  else {                                          /* Else allocate new private frontier array         */
    if ((coargrafptr->fronloctab = memAlloc (coargrafptr->s.vertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("vdgraphSeparateMlCoarsen: out of memory");
      dgraphExit (&coargrafptr->s);               /* Only free Dgraph since fronloctab not allocated */
      memFree    (*coarmultptr);                  /* Free un-based array                             */
      return     (1);
    }
  }

  *coarmultptr -= coargrafptr->s.baseval;         /* Base the multinode array */

  return (0);
}

/* This routine is the reduction-loc operator which
** returns in inout[2] the rank of the process which
** holds the best partition.
** It returns:
** - void  : in all cases.
*/

static
void
vdgraphSeparateMlOpBest (
const Gnum * const          in,                   /* First operand                               */
Gnum * const                inout,                /* Second and output operand                   */
const int * const           len,                  /* Number of instances ; should be 1, not used */
const MPI_Datatype * const  typedat)              /* MPI datatype ; not used                     */
{
  inout[5] |= in[5];                              /* Memory error flag */

  if (inout[0] == 1) {                            /* Handle cases when at least one of them is erroneous */
    if (in[0] == 1) {
      if (inout[1] > in[1]) {                     /* To enforce commutativity, always keep smallest process number */
        inout[1] = in[1];
        inout[2] = in[2];
      }
      return;
    }

    inout[0] = in[0];                             /* Validity flag      */
    inout[1] = in[1];                             /* Lead process rank  */
    inout[2] = in[2];                             /* Lead process color */
    inout[3] = in[3];                             /* Separator size     */
    inout[4] = in[4];                             /* Parts imbalance    */
    return;
  }
  else if (in[0] == 1)
    return;

  if ((in[3] < inout[3]) ||                       /* Select best partition */
      ((in[3] == inout[3]) &&
       ((in[4] < inout[4]) ||
	((in[4] == inout[4]) && (in[1] < inout[1]))))) {
    inout[1] = in[1];
    inout[2] = in[2];
    inout[3] = in[3];
    inout[4] = in[4];
  }
}

/* This routine packs the neighbor data to be sent
** to one of the neighbors by part number.
** It returns:
** - void  : in all cases.
*/

static
void
vdgraphSeparateMlPack (
Gnum * restrict const       dataloctab,
const Gnum                  datalocidx,
Gnum * restrict             ssndcnttab)
{
  Gnum                finevertsndnbr0;
  Gnum                finevertsndnbr1;
  Gnum                finevertsndnbr2;
  Gnum                datalocnbr;

  finevertsndnbr0 = ssndcnttab[0];
  finevertsndnbr1 = ssndcnttab[1];
  finevertsndnbr2 = ssndcnttab[2];
  datalocnbr = finevertsndnbr0 + finevertsndnbr1 + finevertsndnbr2;

  if (datalocnbr <= datalocidx) {                 /* If arrays do not overlap */
    Gnum * restrict             datadsttab = dataloctab + datalocidx;
    const Gnum * restrict const datasrctab = dataloctab + datalocidx * 2;
    Gnum                        datasrcnum;
    Gnum                        partidxtab[3];

    partidxtab[0] = 0;
    partidxtab[1] = finevertsndnbr0;
    partidxtab[2] = finevertsndnbr0 + finevertsndnbr1;
    for (datasrcnum = 0, datalocnbr <<= 1; datasrcnum < datalocnbr; ) { /* Work on pairs of Gnum's */
      Gnum                finevertglbnum;
      Gnum                finepartval;

      finevertglbnum = datasrctab[datasrcnum ++];
      finepartval    = datasrctab[datasrcnum ++];
      datadsttab[partidxtab[finepartval] ++] = finevertglbnum;
    }
  }
  else {                                          /* Arrays do overlap */
    Gnum                datadstnum;
    Gnum                datasrcnum;
    Gnum                datasrcnnd;
    Gnum                datasrcnxt;

    datadstnum = datalocidx;
    for (datasrcnum = datalocidx << 1, datasrcnnd = datasrcnum + (finevertsndnbr0 << 1), datasrcnxt = datasrcnnd; /* Work on pairs of Gnum's */
         datasrcnum < datasrcnnd; ) {
      Gnum                finevertglbnum;
      Gnum                finepartval;

      finevertglbnum = dataloctab[datasrcnum ++];
      finepartval    = dataloctab[datasrcnum ++];
      if (finepartval != 0) {
        Gnum                finevertglbtmp;

#ifdef SCOTCH_DEBUG_VDGRAPH2
        if ((finepartval < 1) || (finepartval > 2)) {
          errorPrint ("vdgraphSeparateMlPack: internal error (1)");
          return;
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        while (dataloctab[datasrcnxt + 1] != 0) { /* Find first vertex of part zero in next block */
          datasrcnxt += 2;
#ifdef SCOTCH_DEBUG_VDGRAPH2
          if (datasrcnxt >= ((datalocidx + datalocnbr) << 1)) {
            errorPrint ("vdgraphSeparateMlPack: internal error (2)");
            return;
          }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        }

        finevertglbtmp = dataloctab[datasrcnxt];
        dataloctab[datasrcnxt ++] = finevertglbnum;
        dataloctab[datasrcnxt ++] = finepartval;
        finevertglbnum = finevertglbtmp;
      }
      dataloctab[datadstnum ++] = finevertglbnum;
    }

    for (datasrcnnd += finevertsndnbr1 << 1, datasrcnxt = datasrcnnd; /* Work on pairs of Gnum's */
         datasrcnum < datasrcnnd; ) {
      Gnum                finevertglbnum;
      Gnum                finepartval;

      finevertglbnum = dataloctab[datasrcnum ++];
      finepartval    = dataloctab[datasrcnum ++];
      if (finepartval != 1) {
        Gnum                finevertglbtmp;

#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (finepartval != 2) {
          errorPrint ("vdgraphSeparateMlPack: internal error (3)");
          return;
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        while (dataloctab[datasrcnxt + 1] != 1) { /* Find first vertex of part one in next block */
          datasrcnxt += 2;
#ifdef SCOTCH_DEBUG_VDGRAPH2
          if (datasrcnxt >= ((datalocidx + datalocnbr) << 1)) {
            errorPrint ("vdgraphSeparateMlPack: internal error (4)");
            return;
          }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        }

        finevertglbtmp = dataloctab[datasrcnxt];
        dataloctab[datasrcnxt ++] = finevertglbnum;
        dataloctab[datasrcnxt ++] = finepartval;
        finevertglbnum = finevertglbtmp;
      }
      dataloctab[datadstnum ++] = finevertglbnum;
    }

    for (datasrcnnd += finevertsndnbr2 << 1; datasrcnum < datasrcnnd; ) { /* Work on pairs of Gnum's */
      Gnum                finevertglbnum;

      finevertglbnum = dataloctab[datasrcnum];
#ifdef SCOTCH_DEBUG_VDGRAPH2
      if (dataloctab[datasrcnum + 1] != 2) {
        errorPrint ("vdgraphSeparateMlPack: internal error (5)");
        return;
      }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
      datasrcnum += 2;
      dataloctab[datadstnum ++] = finevertglbnum;
    }
  }
}

/* This routine propagates the separation of the
** coarser graph back to the finer graph, according
** to the multinode table of collapsed vertices.
** After the separation is propagated, it finishes
** to compute the parameters of the finer graph that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse graph data has been propagated to the fine graph.
** - !0  : on error.
*/

static
int
vdgraphSeparateMlUncoarsen (
Vdgraph * restrict                        finegrafptr, /*+ Finer graph           +*/
const Vdgraph * restrict const            coargrafptr, /*+ Coarser graph         +*/
const DgraphCoarsenMulti * restrict const coarmulttax) /*+ Based multinode array +*/
{
  Gnum                    coarvertnum;
  Gnum                    finevertlocadj;
  Gnum                    finecomplocload0;
  Gnum                    finecomplocload2;
  Gnum                    finecomplocsize1;
  Gnum                    finecomplocsize2;
  Gnum * restrict         srcvdattab;
  Gnum *                  ssnddattab;             /* TRICK: holds vrcvcnttab, vsnddsptab, vrcvdsptab                   */
  Gnum *                  vrcvdattab;             /* TRICK: overlaps with vsnddattab before packing [norestrict:async] */
  Gnum *                  vsnddattab;             /* [norestrict:async]                                                */
  int *                   vrcvcnttab;
  int *                   vsndcnttab;
  int *                   vrcvdsptab;
  int *                   vsnddsptab;
  int                     vrcvdspnum;
  int                     vsnddspnum;
  Gnum                    vrcvdatnum;
  MPI_Datatype            besttypedat;            /* Data type for finding best bipartition              */
  MPI_Op                  bestoperdat;            /* Handle of MPI operator for finding best bipartition */
  Gnum                    reduloctab[6];
  Gnum                    reduglbtab[6];
  int                     procnum;
  const Gnum * restrict   fineveloglbtax;
  GraphPart * restrict    finepartglbtax;

  Gnum * restrict const finefronloctab = finegrafptr->fronloctab;

  reduloctab[5] = 0;                              /* Assume everything is fine                      */
  if (finegrafptr->partgsttax == NULL) {          /* If partition array not yet allocated           */
    if (dgraphGhst (&finegrafptr->s) != 0) {      /* Create ghost edge array and compute vertgstnbr */
      errorPrint ("vdgraphSeparateMlUncoarsen: cannot compute ghost edge array");
      reduloctab[5] = 1;                          /* Allocated data will be freed along with graph structure */
    }
    else if ((finegrafptr->partgsttax = (GraphPart *) memAlloc (finegrafptr->s.vertgstnbr * sizeof (GraphPart))) == NULL) {
      errorPrint ("vdgraphSeparateMlUncoarsen: out of memory (1)");
      reduloctab[5] = 1;                          /* Allocated data will be freed along with graph structure */
    }
    else
      finegrafptr->partgsttax -= finegrafptr->s.baseval;
  }

  if (coargrafptr == NULL) {                      /* If coarser graph not provided                      */
#ifdef SCOTCH_DEBUG_BDGRAPH1                      /* Communication cannot be overlapped by a useful one */
    if (MPI_Allreduce (&reduloctab[5], &reduglbtab[5], 1, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("vdgraphSeparateMlUncoarsen: communication error (1)");
      return     (1);
    }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
    reduglbtab[5] = reduloctab[5];
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
    if (reduglbtab[5] != 0)
      return (1);

    vdgraphZero (finegrafptr);                    /* Assign all vertices to part 0 */

    return (0);
  }

  if (memAllocGroup ((void **) (void *)
                     &vsndcnttab, (size_t) (finegrafptr->s.procglbnbr *     sizeof (int)), /* TRICK: srcvdattab after ssnddattab, after vsndcnttab     */
                     &ssnddattab, (size_t) (finegrafptr->s.procglbnbr * 3 * sizeof (Gnum)), /* TRICK: ssnddattab is vrcvcnttab, vsnddsptab, vrcvdsptab */
                     &srcvdattab, (size_t) (finegrafptr->s.procglbnbr * 3 * sizeof (Gnum)),
                     &vsnddattab, (size_t) (coargrafptr->s.vertlocnbr * 2 * sizeof (Gnum)), /* TRICK: vsnddattab overlaps with vrcvdattab */
                     &vrcvdattab, (size_t) (MAX ((coargrafptr->s.vertlocnbr * 2), finegrafptr->s.vertlocnbr) * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("vdgraphSeparateMlUncoarsen: out of memory (2)");
    reduloctab[5] = 1;
  }

  if (coargrafptr->s.procglbnbr <= 0) {           /* If unused folded coargrafptr   */
    reduloctab[0] = 1;                            /* Set it as invalid              */
    reduloctab[1] = 0;                            /* Useless rank                   */
    reduloctab[2] = 1;                            /* Color is not the one of folded */
    reduloctab[3] =                               /* Prevent Valgrind from yelling  */
    reduloctab[4] = 0;
  }
  else {
    reduloctab[0] = ((coargrafptr->compglbsize[0] == 0) || /* Empty separated parts are deemed invalid */
                     (coargrafptr->compglbsize[1] == 0)) ? 1 : 0;
    reduloctab[1] = finegrafptr->s.proclocnum;    /* Set rank and color key according to coarse graph (sub)communicator */
    reduloctab[2] = finegrafptr->s.prockeyval;
    reduloctab[3] = coargrafptr->compglbsize[2];
    reduloctab[4] = coargrafptr->compglbloaddlt;
  }

  if ((MPI_Type_contiguous (6, GNUM_MPI, &besttypedat)                                != MPI_SUCCESS) ||
      (MPI_Type_commit (&besttypedat)                                                 != MPI_SUCCESS) ||
      (MPI_Op_create ((MPI_User_function *) vdgraphSeparateMlOpBest, 1, &bestoperdat) != MPI_SUCCESS)) {
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (2)");
    return     (1);
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 1, besttypedat, bestoperdat, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (3)");
    return     (1);
  }

  if ((MPI_Op_free   (&bestoperdat) != MPI_SUCCESS) ||
      (MPI_Type_free (&besttypedat) != MPI_SUCCESS)) {
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (4)");
    return     (1);
  }

  if (reduglbtab[5] != 0) {                       /* If memory error, return                     */
    if (vsndcnttab != NULL)                       /* Partgsttax will be freed at the above level */
      memFree (vsndcnttab);
    return (1);
  }

  if (reduglbtab[0] == 1) {                       /* If all possible partitions are invalid */
#ifdef SCOTCH_DEBUG_BDGRAPH2
    errorPrintW ("vdgraphSeparateMlUncoarsen: no valid partition");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    return (1);                                   /* All invalid partitions will lead to low method be applied at upper level */
  }

  finevertlocadj = finegrafptr->s.procvrttab[finegrafptr->s.proclocnum] - finegrafptr->s.baseval;
  fineveloglbtax = (finegrafptr->s.veloloctax != NULL) ? (finegrafptr->s.veloloctax - finevertlocadj) : NULL; /* Array can be indexed with global vertex numbers */
  finepartglbtax = finegrafptr->partgsttax - finevertlocadj;

  finegrafptr->complocload[0] =
  finegrafptr->complocload[1] =
  finegrafptr->complocload[2] =
  finegrafptr->complocsize[0] =
  finegrafptr->complocsize[1] =
  finegrafptr->complocsize[2] = 0;
  
#ifdef SCOTCH_DEBUG_VDGRAPH2
  memSet (finegrafptr->partgsttax + finegrafptr->s.baseval, 3, finegrafptr->s.vertgstnbr * sizeof (GraphPart)); /* Mark all vertices as unvisited */
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  memSet (vsndcnttab, 0, ((byte *) srcvdattab) - ((byte *) vsndcnttab)); /* TRICK: Assume process has nothing to send in vsndcnttab and ssnddattab */

  if (reduglbtab[2] == (Gnum) coargrafptr->s.prockeyval) { /* If we belong to the group of the lead process, we must browse and send local data */
    Gnum                fineveloval;
    Gnum                finevertsndnbr1;
    Gnum                finevertsndnbr2;
    Gnum                finevertglbmin;
    Gnum                finevertglbmax;
    Gnum                finevertglbnnd;
    Gnum                vsnddatnbr;
    Gnum                vsnddatnum;
    Gnum                vsnddattmp;

    const GraphPart * restrict const  coarpartgsttax = coargrafptr->partgsttax;

    fineveloval = 1;                              /* Assume no vertex loads */
    finevertglbmin   = finegrafptr->s.procvrttab[finegrafptr->s.proclocnum];
    finevertglbmax   = finegrafptr->s.procvrttab[finegrafptr->s.proclocnum + 1];
    finecomplocload0 =
    finecomplocload2 =
    finecomplocsize1 =
    finecomplocsize2 = 0;
    for (coarvertnum = coargrafptr->s.baseval, vsnddatnbr = 0;
         coarvertnum < coargrafptr->s.vertlocnnd; coarvertnum ++) {
      Gnum                finevertglbnum;
      GraphPart           coarpartval;

      coarpartval = coarpartgsttax[coarvertnum];
#ifdef SCOTCH_DEBUG_VDGRAPH2
      if ((coarpartval < 0) || (coarpartval > 2)) {
        errorPrint ("vdgraphSeparateMlUncoarsen: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

      finevertglbnum = coarmulttax[coarvertnum].vertglbnum[0];
      while (1) {                                 /* Loop on both fine vertices of multinode */
        Gnum                finepartval;

        finepartval = (Gnum) coarpartval;
        if ((finevertglbnum >= finevertglbmin) && (finevertglbnum < finevertglbmax)) { /* Vertex is a local one */
#ifdef SCOTCH_DEBUG_VDGRAPH2
          if (finepartglbtax[finevertglbnum] != 3) {
            errorPrint ("vdgraphSeparateMlUncoarsen: internal error (3)");
            return (1);
          }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
          finepartglbtax[finevertglbnum] = coarpartval;

          finecomplocsize1 += finepartval & 1;    /* One extra vertex created in part 1 if (coarpartval == 1) */
          if (fineveloglbtax != NULL)
            fineveloval = fineveloglbtax[finevertglbnum];
          if (coarpartval == 2) {
            finecomplocload2 += fineveloval;
            finefronloctab[finecomplocsize2 ++] = finevertglbnum - finevertlocadj;
          }
          else
            finecomplocload0 += fineveloval & (finepartval - 1);
        }
        else {                                    /* Non local vertex         */
          vsnddattab[vsnddatnbr ++] = finevertglbnum; /* Store index and part */
          vsnddattab[vsnddatnbr ++] = finepartval;
        }

        if (finevertglbnum == coarmulttax[coarvertnum].vertglbnum[1]) /* If single-vertex multinode or both vertices processed */
          break;
        finevertglbnum = coarmulttax[coarvertnum].vertglbnum[1]; /* Process second multinode */
      }
    }

    finegrafptr->complocload[0] = finecomplocload0; /* Account for local vertices already processed */
    finegrafptr->complocload[2] = finecomplocload2;
    finegrafptr->complocsize[1] = finecomplocsize1;
    finegrafptr->complocsize[2] = finecomplocsize2;

    intSort2asc1 ((void *) vsnddattab, vsnddatnbr >> 1); /* Sort vertices to send by ascending global numbers */

    finevertsndnbr1 =
    finevertsndnbr2 = 0;
    for (vsnddatnum = vsnddattmp = 0, procnum = 0, finevertglbnnd = finegrafptr->s.procvrttab[1];
         vsnddatnum < vsnddatnbr; ) {
      Gnum                finevertglbnum;
      Gnum                finepartval;

      finevertglbnum = vsnddattab[vsnddatnum];
      finepartval    = vsnddattab[vsnddatnum + 1];
      if (finevertglbnum >= finevertglbnnd) {
        Gnum                finevertsndnbr;

        finevertsndnbr = (vsnddatnum - vsnddattmp) >> 1;
        finevertsndnbr2 >>= 1;
        vsndcnttab[procnum] = (int) finevertsndnbr;
        ssnddattab[3 * procnum]     = finevertsndnbr - finevertsndnbr1 - finevertsndnbr2;
        ssnddattab[3 * procnum + 1] = finevertsndnbr1;
        ssnddattab[3 * procnum + 2] = finevertsndnbr2;
        vdgraphSeparateMlPack (vsnddattab, vsnddattmp >> 1, ssnddattab + (3 * procnum));

        do
          finevertglbnnd = finegrafptr->s.procvrttab[(++ procnum) + 1];
        while (finevertglbnum >= finevertglbnnd);

        vsnddattmp      = vsnddatnum;             /* Set startpoint for new neighbor */
        finevertsndnbr1 =
        finevertsndnbr2 = 0;
      }

      vsnddatnum += 2;
      finevertsndnbr1 += finepartval & 1;         /* Count number of vertices in part 1       */
      finevertsndnbr2 += finepartval & 2;         /* Count twice number of vertices in part 2 */
    }
    finevertsndnbr2 >>= 1;                        /* Complete data for last receiver process */
    vsndcnttab[procnum] = (int) ((vsnddatnum - vsnddattmp) >> 1);
    ssnddattab[3 * procnum]     = ((vsnddatnum - vsnddattmp) >> 1) - finevertsndnbr1 - finevertsndnbr2;
    ssnddattab[3 * procnum + 1] = finevertsndnbr1;
    ssnddattab[3 * procnum + 2] = finevertsndnbr2;
    vdgraphSeparateMlPack (vsnddattab, (Gnum) vsnddattmp >> 1, ssnddattab + (3 * procnum));

#ifdef SCOTCH_DEBUG_VDGRAPH2
    if ((ssnddattab[3 * finegrafptr->s.proclocnum]     != 0) || /* One should never send something to itself */
        (ssnddattab[3 * finegrafptr->s.proclocnum + 1] != 0) ||
        (ssnddattab[3 * finegrafptr->s.proclocnum + 2] != 0)) {
      errorPrint ("vdgraphSeparateMlUncoarsen: internal error (4)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
  }

  if (MPI_Alltoall (ssnddattab, 3, GNUM_MPI, srcvdattab, 3, GNUM_MPI, finegrafptr->s.proccomm) != MPI_SUCCESS) { /* Exchange sizes */
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (2)");
    return     (1);
  }

  vrcvcnttab = (int *) ssnddattab;                /* TRICK: re-use ssnddattab */
  vsnddsptab = vrcvcnttab + finegrafptr->s.procglbnbr;
  vrcvdsptab = vrcvcnttab + finegrafptr->s.procglbnbr * 2;
  for (procnum = 0, vsnddspnum = vrcvdspnum = 0; procnum < finegrafptr->s.procglbnbr; procnum ++) { /* Compute size of data to exchange */
    vrcvcnttab[procnum] = (int) (srcvdattab[3 * procnum] + srcvdattab[3 * procnum + 1] + srcvdattab[3 * procnum + 2]);
    vrcvdsptab[procnum] = vrcvdspnum;
    vsnddsptab[procnum] = vsnddspnum;
    vrcvdspnum += vrcvcnttab[procnum];
    vsnddspnum += vsndcnttab[procnum];
  }

  if (MPI_Alltoallv (vsnddattab, vsndcnttab, vsnddsptab, GNUM_MPI, /* Exchange data */
                     vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (3)");
    return     (1);
  }

  finecomplocload0 = finegrafptr->complocload[0];
  finecomplocload2 = finegrafptr->complocload[2];
  finecomplocsize1 = finegrafptr->complocsize[1];
  finecomplocsize2 = finegrafptr->complocsize[2];
  for (procnum = 0, vrcvdatnum = 0;               /* Process partition data per process number */
       procnum < finegrafptr->s.procglbnbr; procnum ++) {
    Gnum                vrcvdatnnd;

    vrcvdatnnd = vrcvdatnum + srcvdattab[3 * procnum];
    if (fineveloglbtax != NULL) {
      for ( ; vrcvdatnum < vrcvdatnnd; vrcvdatnum ++) {
        Gnum                finevertglbnum;

        finevertglbnum = vrcvdattab[vrcvdatnum];
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (finepartglbtax[finevertglbnum] != 3) {
          errorPrint ("vdgraphSeparateMlUncoarsen: internal error (5)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        finepartglbtax[finevertglbnum] = 0;
        finecomplocload0 += fineveloglbtax[finevertglbnum];
      }
    }
    else {
      finecomplocload0 += srcvdattab[3 * procnum];
      for ( ; vrcvdatnum < vrcvdatnnd; vrcvdatnum ++) {
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (finepartglbtax[vrcvdattab[vrcvdatnum]] != 3) {
          errorPrint ("vdgraphSeparateMlUncoarsen: internal error (6)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        finepartglbtax[vrcvdattab[vrcvdatnum]] = 0;
      }
    }

    finecomplocsize1 += srcvdattab[3 * procnum + 1];
    vrcvdatnnd = vrcvdatnum + srcvdattab[3 * procnum + 1];
    for ( ; vrcvdatnum < vrcvdatnnd; vrcvdatnum ++) {
#ifdef SCOTCH_DEBUG_VDGRAPH2
      if (finepartglbtax[vrcvdattab[vrcvdatnum]] != 3) {
        errorPrint ("vdgraphSeparateMlUncoarsen: internal error (7)");
        return (1);
      }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
      finepartglbtax[vrcvdattab[vrcvdatnum]] = 1;
    }

    vrcvdatnnd = vrcvdatnum + srcvdattab[3 * procnum + 2];
    if (fineveloglbtax != NULL) {
      for ( ; vrcvdatnum < vrcvdatnnd; vrcvdatnum ++) {
        Gnum                finevertglbnum;

        finevertglbnum = vrcvdattab[vrcvdatnum];
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (finepartglbtax[finevertglbnum] != 3) {
          errorPrint ("vdgraphSeparateMlUncoarsen: internal error (8)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        finefronloctab[finecomplocsize2 ++] = finevertglbnum - finevertlocadj;
        finepartglbtax[finevertglbnum] = 2;
        finecomplocload2 += fineveloglbtax[finevertglbnum];
      }
    }
    else {
      finecomplocload2 += srcvdattab[3 * procnum + 2];
      for ( ; vrcvdatnum < vrcvdatnnd; vrcvdatnum ++) {
        Gnum                finevertglbnum;

        finevertglbnum = vrcvdattab[vrcvdatnum];
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if (finepartglbtax[finevertglbnum] != 3) {
          errorPrint ("vdgraphSeparateMlUncoarsen: internal error (9)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        finefronloctab[finecomplocsize2 ++] = finevertglbnum - finevertlocadj;
        finepartglbtax[finevertglbnum] = 2;
      }
    }
  }

  finegrafptr->complocload[0] = finecomplocload0;
  finegrafptr->complocload[1] = finegrafptr->s.velolocsum - finecomplocload0 - finecomplocload2;
  finegrafptr->complocload[2] = finecomplocload2;
  finegrafptr->complocsize[0] = finegrafptr->s.vertlocnbr - finecomplocsize1 - finecomplocsize2;
  finegrafptr->complocsize[1] = finecomplocsize1;
  finegrafptr->complocsize[2] = finecomplocsize2;

  memFree (vsndcnttab);                           /* Free group leader */

  reduloctab[0] = finegrafptr->complocload[0];
  reduloctab[1] = finegrafptr->complocload[1];
  reduloctab[2] = finegrafptr->complocload[2];
  reduloctab[3] = finegrafptr->complocsize[0];
  reduloctab[4] = finegrafptr->complocsize[1];
  reduloctab[5] = finegrafptr->complocsize[2];

  if (MPI_Allreduce (reduloctab, reduglbtab, 6, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateMlUncoarsen: communication error (4)");
    return     (1);
  }

  finegrafptr->compglbload[0] = reduglbtab[0];
  finegrafptr->compglbload[1] = reduglbtab[1];
  finegrafptr->compglbload[2] = reduglbtab[2];
  finegrafptr->compglbsize[0] = reduglbtab[3];
  finegrafptr->compglbsize[1] = reduglbtab[4];
  finegrafptr->compglbsize[2] = reduglbtab[5];
  finegrafptr->compglbloaddlt = reduglbtab[0] - reduglbtab[1];

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (finegrafptr) != 0) {
    errorPrint ("vdgraphSeparateMlUncoarsen: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  return (0);
}

/* This routine recursively performs the
** separation recursion.
** It returns:
** - 0   : if separator could be computed.
** - !0  : on error.
*/

static
int
vdgraphSeparateMl2 (
Vdgraph * restrict const             grafptr,     /* Vertex-separation graph */
const VdgraphSeparateMlParam * const paraptr)     /* Method parameters       */
{
  Vdgraph                       coargrafdat;
  DgraphCoarsenMulti * restrict coarmulttax;
  int                           o;

  if (grafptr->s.procglbnbr <= 1) {               /* No need to stay parallel              */
    if (((o = vdgraphSeparateMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph */
        ((o = vdgraphSeparateSt (grafptr, paraptr->stratseq)) != 0)) {
#ifdef SCOTCH_DEBUG_VDGRAPH2
      errorPrintW ("vdgraphSeparateMl2: cannot apply sequential strategy");
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
    }
    return (o);
  }

  coarmulttax = NULL;                             /* Assume multinode array is not allocated */
  if (vdgraphSeparateMlCoarsen (grafptr, &coargrafdat, &coarmulttax, paraptr) == 0) {
    o = (coargrafdat.s.procglbnbr == 0) ? 0 : vdgraphSeparateMl2 (&coargrafdat, paraptr); /* Apply recursion on coarsened graph if it exists */

    if ((o == 0) &&
        ((o = vdgraphSeparateMlUncoarsen (grafptr, &coargrafdat, coarmulttax)) == 0) &&
        ((o = vdgraphSeparateSt          (grafptr, paraptr->stratasc))         != 0)) { /* Apply ascending strategy if uncoarsening worked */
#ifdef SCOTCH_DEBUG_VDGRAPH2
      errorPrintW ("vdgraphSeparateMl2: cannot apply ascending strategy");
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
    }

    if (coargrafdat.fronloctab == grafptr->fronloctab) /* If coarse graph shares fronloctab with fine graph */
      coargrafdat.fronloctab = NULL;              /* Prevent fronloctab of fine graph from being freed      */
    vdgraphExit (&coargrafdat);

    if (coarmulttax != NULL)                      /* If multinode array has been allocated */
      memFree (coarmulttax + grafptr->s.baseval); /* Free array                            */

    if (o == 0)                                   /* If multi-level failed, apply low strategy as fallback */
      return (o);
  }

  if (((o = vdgraphSeparateMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph            */
      ((o = vdgraphSeparateSt          (grafptr, paraptr->stratlow)) != 0)) { /* Apply low strategy */
#ifdef SCOTCH_DEBUG_VDGRAPH2
    errorPrintW ("vdgraphSeparateMl2: cannot apply low strategy");
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
  }

  return (o);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the muti-level separation.
** It returns:
** - 0 : if separator could be computed.
** - 1 : on error.
*/

int
vdgraphSeparateMl (
Vdgraph * const                       grafptr,    /*+ Vertex-separation graph +*/
const VdgraphSeparateMlParam * const  paraptr)    /*+ Method parameters       +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = grafptr->levlnum;                     /* Save graph level               */
  grafptr->levlnum = 0;                           /* Initialize coarsening level    */
  o = vdgraphSeparateMl2 (grafptr, paraptr);      /* Perform multi-level separation */
  grafptr->levlnum = levlnum;                     /* Restore graph level            */

  return (o);
}
