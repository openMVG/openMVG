/* Copyright 2007-2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_ml.c                     **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module bipartitions a distributed  **/
/**                graph using a multi-level scheme.       **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 30 oct 2007     **/
/**                                 to   : 14 apr 2011     **/
/**              : # Version 6.0  : from : 11 sep 2011     **/
/**                                 to   : 28 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_BIPART_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "arch.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "bdgraph.h"
#include "bdgraph_bipart_ml.h"
#include "bdgraph_bipart_st.h"

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
** - 1  : if threshold reached or on error.
*/

static
int
bdgraphBipartMlCoarsen (
Bdgraph * restrict const              finegrafptr, /*+ Finer graph                         +*/
Bdgraph * restrict const              coargrafptr, /*+ Coarser graph to build              +*/
DgraphCoarsenMulti * restrict * const coarmultptr, /*+ Pointer to multinode table to build +*/
const BdgraphBipartMlParam * const    paraptr)    /*+ Method parameters                    +*/
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
#ifdef SCOTCH_DEBUG_BGRAPH2
    default :
      errorPrint ("bdgraphBipartMlCoarsen: invalid parameter");
      return     (1);
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
  }
  if ((finegrafptr->s.vertglbnbr / finegrafptr->s.procglbnbr) > paraptr->foldmax) /* If no need to fold */
    foldval = DGRAPHCOARSENNONE;

  *coarmultptr = NULL;                            /* Let the routine create the multinode array */
  dgraphInit (&coargrafptr->s, finegrafptr->s.proccomm); /* Re-use fine graph communicator      */
  if (dgraphCoarsen (&finegrafptr->s, &coargrafptr->s, coarmultptr, paraptr->passnbr,
                     paraptr->coarnbr, paraptr->coarrat, foldval) != 0)
    return (1);                                   /* Return if coarsening failed */

  *coarmultptr -= coargrafptr->s.baseval;         /* Base pointer to multinode array    */
  coargrafptr->partgsttax = NULL;                 /* Do not allocate partition data yet */
  coargrafptr->fronloctab = NULL;
  coargrafptr->fronglbnbr = 0;

  if (coargrafptr->s.procglbnbr == 0) {           /* Not owner of graph */
    coargrafptr->veexloctax = NULL;
    return (0);
  }

  if (finegrafptr->veexloctax != NULL) {          /* Merge external gains for coarsened vertices */
    DgraphCoarsenMulti * restrict coarmulttax;
    Gnum * restrict               coarveexloctax;
    Gnum                          coarvertlocnum;

    if ((coarveexloctax = (Gnum *) memAlloc (coargrafptr->s.vertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("bdgraphBipartMlCoarsen: out of memory");
      dgraphExit (&coargrafptr->s);               /* Only free Dgraph since veexloctax not allocated */
      memFree    (*coarmultptr + coargrafptr->s.baseval);
      return     (1);
    }
    coarveexloctax -= coargrafptr->s.baseval;
    coargrafptr->veexloctax = coarveexloctax;
    coarmulttax = *coarmultptr;

    for (coarvertlocnum = coargrafptr->s.baseval; coarvertlocnum < coargrafptr->s.vertlocnnd; coarvertlocnum++) {
      Gnum		finevertnum0;             /* First multinode vertex  */
      Gnum              finevertnum1;             /* Second multinode vertex */

      finevertnum0 = coarmulttax[coarvertlocnum].vertglbnum[0];
      finevertnum1 = coarmulttax[coarvertlocnum].vertglbnum[1];
      coarveexloctax[coarvertlocnum] = (finevertnum0 != finevertnum1)
                                       ? finegrafptr->veexloctax[finevertnum0] + finegrafptr->veexloctax[finevertnum1] 
                                       : finegrafptr->veexloctax[finevertnum0];
    }
  }
  else                                            /* If fine graph does not have external gains */
    coargrafptr->veexloctax = NULL;               /* Coarse graph does not have external gains  */

  coargrafptr->veexglbsum       = finegrafptr->veexglbsum;
  coargrafptr->compglbload0min  = finegrafptr->compglbload0min; /* Only set constant partition parameters as others will be set on uncoarsening */
  coargrafptr->compglbload0max  = finegrafptr->compglbload0max;
  coargrafptr->compglbload0avg  = finegrafptr->compglbload0avg;
  coargrafptr->commglbloadextn0 = finegrafptr->commglbloadextn0;
  coargrafptr->commglbgainextn0 = finegrafptr->commglbgainextn0;
  coargrafptr->domndist         = finegrafptr->domndist;
  coargrafptr->domnwght[0]      = finegrafptr->domnwght[0];
  coargrafptr->domnwght[1]      = finegrafptr->domnwght[1];
  coargrafptr->levlnum          = finegrafptr->levlnum + 1;

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
bdgraphBipartMlOpBest (
const Gnum * const          in,                   /* First operand                               */
Gnum * const                inout,                /* Second and output operand                   */
const int * const           len,                  /* Number of instances ; should be 1, not used */
const MPI_Datatype * const  typedat)              /* MPI datatype ; not used                     */
{
  inout[5] |= in[5];                              /* Memory error flag */

  if (inout[0] == 1) {                            /* Handle cases when at least one of them is erroneous */
    if (in[0] == 1) {
      if (inout[1] > in[1])                       /* To enforce commutativity, always keep smallest process number */
        inout[1] = in[1];
        inout[2] = in[2];
      return;
    }

    inout[0] = in[0];                             /* Validity flag        */
    inout[1] = in[1];                             /* Lead process rank    */
    inout[2] = in[2];                             /* Lead process color   */
    inout[3] = in[3];                             /* Communication load   */
    inout[4] = in[4];                             /* Load imbalance       */
    return;
  }
  else if (in[0] == 1)
    return;

  if ((in[3] < inout[3]) ||                       /* Select best partition */
      ((in[3] == inout[3]) && ((in[4] < inout[4]) ||
                               ((in[4] == inout[4]) && (in[1] < inout[1]))))) {
    inout[1] = in[1];
    inout[2] = in[2];
    inout[3] = in[3];
    inout[4] = in[4];
  }
}

/* This routine propagates the bipartitioning of the
** coarser graph back to the finer graph, according
** to the multinode table of collapsed vertices.
** After the bipartitioning is propagated, it finishes
** to compute the parameters of the finer graph that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse graph data has been propagated to fine graph.
** - !0  : on error.
*/

static
int
bdgraphBipartMlUncoarsen (
Bdgraph * restrict                        finegrafptr, /*+ Finer graph     +*/
const Bdgraph * restrict const            coargrafptr, /*+ Coarser graph   +*/
const DgraphCoarsenMulti * restrict const coarmulttax) /*+ Multinode array +*/
{
  Gnum                            baseval;
  Gnum                            finefronlocnbr;
  Gnum                            finefronlocnum;
  Gnum                            fineedlolocval;
  Gnum                            finevertlocadj; /* Global vertex adjustment                            */
  Gnum                            finevertlocnum;
  Gnum                            finevertlocnnd; /* Index for frontier array fronloctab                 */
  Gnum                            finecomplocsize1;
  Gnum                            finecomplocload1;
  Gnum                            finecommlocloadintn;
  Gnum                            finecommlocloadextn;
  Gnum                            finecommlocgainextn;
  int                             vrcvdspnbr;
  int                             vsnddspnbr;
  int * restrict                  vrcvcnttab;
  int * restrict                  vsndcnttab;
  int * restrict                  vrcvdsptab;
  int * restrict                  vsnddsptab;
  Gnum * restrict                 vrcvdattab;
  Gnum * restrict                 vsnddattab;
  Gnum * restrict                 vsndidxtab;
  BdgraphBipartMlSort * restrict  sortloctab;     /* Array of vertices to send to their owner            */
  Gnum                            sortlocnbr;
  Gnum                            sortlocnum;
  int                             procnum;
  MPI_Datatype                    besttypedat;    /* Data type for finding best bipartition              */
  MPI_Op                          bestoperdat;    /* Handle of MPI operator for finding best bipartition */
  Gnum                            reduloctab[6];  /* "6": both for selecting best and propagating data   */
  Gnum                            reduglbtab[6];
  const Gnum * restrict           coarfronloctab;
  GraphPart * restrict            coarpartgsttax;
  GraphPart * restrict            finepartgsttax;
  Gnum * restrict                 finefronloctab;

  const int                   fineprocglbnbr = finegrafptr->s.procglbnbr;
  const Gnum * restrict const fineprocvrttab = finegrafptr->s.procvrttab;
  const Gnum * restrict const fineedgegsttax = finegrafptr->s.edgegsttax;
  const Gnum * restrict const finevertloctax = finegrafptr->s.vertloctax;
  const Gnum * restrict const finevendloctax = finegrafptr->s.vendloctax;
  const Gnum * restrict const fineveloloctax = finegrafptr->s.veloloctax;
  const Gnum * restrict const fineveexloctax = finegrafptr->veexloctax;
  const Gnum * restrict const fineedloloctax = finegrafptr->s.edloloctax;

  reduloctab[5] = 0;                              /* Assume everything is fine                      */
  if (finegrafptr->partgsttax == NULL) {          /* If partition array not yet allocated           */
    if (dgraphGhst (&finegrafptr->s) != 0) {      /* Create ghost edge array and compute vertgstnbr */
      errorPrint ("bdgraphBipartMlUncoarsen: cannot compute ghost edge array");
      reduloctab[5] = 1;                          /* Allocated data will be freed along with graph structure */
    }
    else if ((finegrafptr->partgsttax = (GraphPart *) memAlloc (finegrafptr->s.vertgstnbr * sizeof (GraphPart))) == NULL) {
      errorPrint ("bdgraphBipartMlUncoarsen: out of memory (1)");
      reduloctab[5] = 1;                          /* Allocated data will be freed along with graph structure */
    }
    else if (finegrafptr->partgsttax -= finegrafptr->s.baseval,
             (finegrafptr->fronloctab = (Gnum *) memAlloc (finegrafptr->s.vertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("bdgraphBipartMlUncoarsen: out of memory (2)");
      reduloctab[5] = 1;                          /* Allocated data will be freed along with graph structure */
    }
  }

  if (coargrafptr == NULL) {                      /* If coarser graph not provided                      */
#ifdef SCOTCH_DEBUG_BDGRAPH1                      /* Communication cannot be overlapped by a useful one */
    if (MPI_Allreduce (&reduloctab[5], &reduglbtab[5], 1, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("bdgraphBipartMlUncoarsen: communication error (1)");
      return     (1);
    }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
    reduglbtab[5] = reduloctab[5];
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
    if (reduglbtab[5] != 0)
      return (1);

    bdgraphZero (finegrafptr);                    /* Assign all vertices to part 0 */

#ifdef SCOTCH_DEBUG_BDGRAPH2
    if (bdgraphCheck (finegrafptr) != 0) {
      errorPrint ("bdgraphBipartMlUncoarsen: inconsistent graph data (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

    return (0);
  }

  if (coargrafptr->s.procglbnbr <= 0) {           /* If unused folded coargrafptr   */
    reduloctab[0] = 1;                            /* Set it as invalid              */
    reduloctab[1] = 0;                            /* Useless rank                   */
    reduloctab[2] = 1;                            /* Color is not the one of folded */
    reduloctab[3] =                               /* Prevent Valgrind from yelling  */
    reduloctab[4] = 0;
  }
  else {
    reduloctab[0] = ((coargrafptr->compglbload0 == 0) || /* Empty subdomains are deemed invalid */
                     (coargrafptr->compglbload0 == coargrafptr->s.veloglbsum)) ? 1 : 0;
    reduloctab[1] = finegrafptr->s.proclocnum;    /* Set rank and color key according to coarse graph (sub)communicator */
    reduloctab[2] = finegrafptr->s.prockeyval;
    reduloctab[3] = coargrafptr->commglbload;
    reduloctab[4] = coargrafptr->compglbload0dlt;
  }

  if ((MPI_Type_contiguous (6, GNUM_MPI, &besttypedat)                              != MPI_SUCCESS) ||
      (MPI_Type_commit (&besttypedat)                                               != MPI_SUCCESS) ||
      (MPI_Op_create ((MPI_User_function *) bdgraphBipartMlOpBest, 1, &bestoperdat) != MPI_SUCCESS)) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (2)");
    return     (1);
  }

  if (MPI_Allreduce (reduloctab, reduglbtab, 1, besttypedat, bestoperdat, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (3)");
    return     (1);
  }

  if ((MPI_Op_free   (&bestoperdat) != MPI_SUCCESS) ||
      (MPI_Type_free (&besttypedat) != MPI_SUCCESS)) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (4)");
    return     (1);
  }

  if (reduglbtab[5] != 0)                         /* If memory error, return */
    return (1);

  if (reduglbtab[0] == 1) {                       /* If all possible partitions are invalid */
#ifdef SCOTCH_DEBUG_BDGRAPH2
    errorPrintW ("bdgraphBipartMlUncoarsen: no valid partition");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    return (1);                                   /* All invalid partitions will lead to low method be applied at upper level */
  }

  if (memAllocGroup ((void **) (void *)
                     &vrcvcnttab, (size_t) (fineprocglbnbr * sizeof (int)),
                     &vrcvdsptab, (size_t) (fineprocglbnbr * sizeof (int)),
                     &vsnddsptab, (size_t) (fineprocglbnbr * sizeof (int)),
                     &vsndcnttab, (size_t) (fineprocglbnbr * sizeof (int)),
                     &vsndidxtab, (size_t) (fineprocglbnbr * sizeof (Gnum) * 4), /* TRICK: sortloctab after vsndidxtab after vsndcnttab */
                     &sortloctab, (size_t) (2 * coargrafptr->s.vertlocnbr * sizeof (BdgraphBipartMlSort)), NULL) == NULL) {
    errorPrint ("bdgraphBipartMlUncoarsen: out of memory (3)");
    reduloctab[5] = 1;
  }
#ifdef SCOTCH_DEBUG_BDGRAPH1                      /* Communication cannot be overlapped by a useful one */
  if (MPI_Allreduce (&reduloctab[5], &reduglbtab[5], 1, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (5)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
  reduglbtab[5] = reduloctab[5];
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  if (reduglbtab[5] != 0) {
    if (vrcvcnttab != NULL)
      memFree (vrcvcnttab);
    return (1);
  }

  memSet (vsndcnttab, 0, ((byte *) sortloctab) - ((byte *) vsndcnttab)); /* TRICK: sortloctab after vsndidxtab after vsndcnttab */

  baseval        = finegrafptr->s.baseval;
  coarfronloctab = coargrafptr->fronloctab;
  coarpartgsttax = coargrafptr->partgsttax;
  finepartgsttax = finegrafptr->partgsttax;
  finevertlocnnd = finegrafptr->s.vertlocnnd;
  finevertlocadj = finegrafptr->s.procvrttab[finegrafptr->s.proclocnum] - baseval;
  finefronloctab = finegrafptr->fronloctab;

  finecomplocsize1    = 0;
  finecomplocload1    = 0;
  finecommlocloadextn = 0;
  finecommlocgainextn = 0;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  memSet (finepartgsttax + baseval, ~0, finegrafptr->s.vertgstnbr * sizeof (GraphPart)); /* All vertices are unvisited */
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  finefronlocnbr = 0;
  sortlocnbr     = 0;
  if (reduglbtab[2] == (Gnum) coargrafptr->s.prockeyval) { /* If we belong to the group of the lead process, we must browse and send local data */
    Gnum                coarfronlocnum;
    Gnum                coarvertlocnum;

    for (coarfronlocnum = 0; coarfronlocnum < coargrafptr->fronlocnbr; coarfronlocnum ++)
      coarpartgsttax[coarfronloctab[coarfronlocnum]] |= 2; /* Flag vertex as belonging to frontier */

    for (coarvertlocnum = baseval; coarvertlocnum < coargrafptr->s.vertlocnnd; coarvertlocnum ++) {
      GraphPart           coarpartval;
      Gnum                coarpartmsk;
      Gnum                finevertglbnum;
      Gnum                finevertlocnum;
      int                 i;

      coarpartval = coarpartgsttax[coarvertlocnum];
      coarpartmsk = (Gnum) (coarpartval & 1);

      i = 0;
      do {
        finevertglbnum = coarmulttax[coarvertlocnum].vertglbnum[i];
        finevertlocnum = finevertglbnum - finevertlocadj;

        if ((finevertlocnum >= baseval) &&        /* If vertex is local */
            (finevertlocnum <  finevertlocnnd)) {
#ifdef SCOTCH_DEBUG_BDGRAPH2
          if (finepartgsttax[finevertlocnum] != ((GraphPart) ~0)) {
            errorPrint ("bdgraphBipartMlUncoarsen: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
          finepartgsttax[finevertlocnum] = (coarpartval & 1);
          finecomplocsize1 += coarpartmsk;        /* One extra vertex created in part 1 if (coarpartval == 1) */

          if ((coarpartval & 2) != 0)             /* If coarse frontier vertex, add fine vertex to fine frontier */
            finefronloctab[finefronlocnbr ++] = finevertlocnum;

          if (fineveloloctax != NULL) {
            Gnum                veloval;

            veloval = fineveloloctax[finevertlocnum]; 
            finecomplocload1 += veloval & (- coarpartmsk);
          }
          if (fineveexloctax != NULL) {
            Gnum                veexval;

            veexval = fineveexloctax[finevertlocnum];
            finecommlocloadextn += veexval * coarpartmsk;
            finecommlocgainextn += veexval * (1 - 2 * coarpartmsk);
          }
        }
        else {
          int               procngbnum;
          int               procngbmax;

          procngbnum = 0;
          procngbmax = fineprocglbnbr;
          while ((procngbmax - procngbnum) > 1) { /* Find owner process by dichotomy on procvgbtab */
            int                 procngbmed;

            procngbmed = (procngbmax + procngbnum) / 2;
            if (fineprocvrttab[procngbmed] > finevertglbnum)
              procngbmax = procngbmed;
            else
              procngbnum = procngbmed;
          }

          vsndidxtab[4 * procngbnum + coarpartval] ++; /* One of four counters per process number will be incremented */
          sortloctab[sortlocnbr].vertnum = finevertglbnum;
          sortloctab[sortlocnbr].procnum = ((procngbnum + (fineprocglbnbr * coarpartmsk)) ^ (- (Gnum) (coarpartval >> 1))); /* Embed part and frontier information */
          sortlocnbr ++;
        }

        i ++;                                     /* Process next multinode vertex                 */
      } while (finevertglbnum != coarmulttax[coarvertlocnum].vertglbnum[1]); /* If not single node */
    }

    for (procnum = 0; procnum < fineprocglbnbr; procnum ++) { /* Aggregate data to be sent */
      vsndcnttab[procnum] = vsndidxtab[4 * procnum]     + vsndidxtab[4 * procnum + 1] +
                            vsndidxtab[4 * procnum + 2] + vsndidxtab[4 * procnum + 3];

      if (vsndcnttab[procnum] != 0)               /* If we will send data to neighbor */
        vsndcnttab[procnum] += 3;                 /* Add control data to message size */
    }
  }

  if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, vrcvcnttab, 1, MPI_INT, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (3)");
    return     (1);
  }

  for (procnum = 0, vrcvdspnbr = vsnddspnbr = 0; /* Build communication index arrays */
       procnum < fineprocglbnbr; procnum ++) {
    vrcvdsptab[procnum] = vrcvdspnbr;
    vsnddsptab[procnum] = vsnddspnbr;
    vrcvdspnbr += vrcvcnttab[procnum];
    vsnddspnbr += vsndcnttab[procnum];
  }

  if (memAllocGroup ((void **) (void *)
                     &vrcvdattab, (size_t) (vrcvdspnbr * sizeof (Gnum)),
                     &vsnddattab, (size_t) (vsnddspnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("bdgraphBipartMlUncoarsen: out of memory (4)");
    reduloctab[5] = 1;
  }
#ifdef SCOTCH_DEBUG_BDGRAPH1                      /* Communication cannot be overlapped by a useful one */
  if (MPI_Allreduce (&reduloctab[5], &reduglbtab[5], 1, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm) != MPI_SUCCESS)  {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (4)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_BDGRAPH1 */
  reduglbtab[5] = reduloctab[5];
#endif /* SCOTCH_DEBUG_BDGRAPH1 */
  if (reduglbtab[5] != 0) {
    if (vrcvdattab != NULL)
      memFree (vrcvdattab);
    if (vrcvcnttab != NULL)
      memFree (vrcvcnttab);
    return (1);
  }

  for (procnum = 0; procnum < fineprocglbnbr; procnum ++) {
    Gnum              vsnddspval;

    vsnddspval = vsnddsptab[procnum];
    if (vsndcnttab[procnum] != 0) {
      Gnum              vsndidxnum;

      vsnddattab[vsnddspval]     = vsndidxtab[4 * procnum];
      vsnddattab[vsnddspval + 1] = vsndidxtab[4 * procnum + 1];
      vsnddattab[vsnddspval + 2] = vsndidxtab[4 * procnum + 2];

      vsnddspval += 3;                            /* Compute sub-array indices to pack vertices to be sent */
      vsndidxnum = vsndidxtab[4 * procnum];
      vsndidxtab[4 * procnum] = vsnddspval;
      vsnddspval += vsndidxnum;
      vsndidxnum = vsndidxtab[4 * procnum + 1];
      vsndidxtab[4 * procnum + 1] = vsnddspval;
      vsnddspval += vsndidxnum;
      vsndidxnum = vsndidxtab[4 * procnum + 2];
      vsndidxtab[4 * procnum + 2] = vsnddspval;
      vsnddspval += vsndidxnum;
      vsndidxtab[4 * procnum + 3] = vsnddspval;
    }
  }

  for (sortlocnum = 0; sortlocnum < sortlocnbr; sortlocnum ++) { /* For all vertices to send */
    Gnum              vertglbend;
    Gnum              procngbnum;
    int               partval;

    vertglbend = sortloctab[sortlocnum].vertnum;
    procngbnum = sortloctab[sortlocnum].procnum;

    partval = 0;                                  /* Extract frontier and part data from process number */
    if (procngbnum < 0) {
      partval = 2;
      procngbnum ^= (Gnum) -1;
    }
    if (procngbnum >= fineprocglbnbr) {
      partval |= 1;
      procngbnum -= fineprocglbnbr;
    }

#ifdef SCOTCH_DEBUG_BDGRAPH2
    if (((partval < 3) && (vsndidxtab[4 * procngbnum + partval] >= vsndidxtab[4 * procngbnum + partval + 1])) ||
        (vsndidxtab[4 * procngbnum + partval] >= (vsnddsptab[procngbnum] + vsndcnttab[procngbnum]))) {
      errorPrint ("bdgraphBipartMlUncoarsen: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    vsnddattab[vsndidxtab[4 * procngbnum + partval] ++] = vertglbend; /* Pack vertex in proper sub-array */
  }

  if (MPI_Alltoallv (vsnddattab, vsndcnttab, vsnddsptab, GNUM_MPI, 
                     vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI, finegrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartMlUncoarsen: communication error (5)");
    return     (1);
  }
    
  for (procnum = 0; procnum < fineprocglbnbr; ++ procnum) { /* Update local ones from the buffer for receiving data */
    Gnum                vrcvidxnum;
    Gnum                vrcvidxnnd;

    if (vrcvcnttab[procnum] == 0)                 /* If nothing received from this process, skip it */
      continue;

    finecomplocsize1 += (vrcvcnttab[procnum] - 3) - vrcvdattab[vrcvdsptab[procnum]] - vrcvdattab[vrcvdsptab[procnum] + 2];

    for (vrcvidxnum = vrcvdsptab[procnum] + 3, vrcvidxnnd = vrcvidxnum + vrcvdattab[vrcvdsptab[procnum]]; /* Vertices in sub-array 0 */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                finevertlocnum;

      finevertlocnum = vrcvdattab[vrcvidxnum] - finevertlocadj;
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((finevertlocnum < baseval) || (finevertlocnum >= finevertlocnnd)) {
        errorPrint ("bdgraphBipartMlUncoarsen: internal error (4)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      finepartgsttax[finevertlocnum] = 0;
    }

    for (vrcvidxnnd = vrcvidxnum + vrcvdattab[vrcvdsptab[procnum] + 1]; /* Vertices in sub-array 1 */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                finevertlocnum;

      finevertlocnum = vrcvdattab[vrcvidxnum] - finevertlocadj;
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((finevertlocnum < baseval) || (finevertlocnum >= finevertlocnnd)) {
        errorPrint ("bdgraphBipartMlUncoarsen: internal error (5)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      finepartgsttax[finevertlocnum] = 1;
      if (fineveloloctax != NULL)
        finecomplocload1 += fineveloloctax[finevertlocnum];
      if (fineveexloctax != NULL) {
        Gnum                veexval;

        veexval = fineveexloctax[finevertlocnum];
        finecommlocloadextn += veexval;
        finecommlocgainextn -= veexval;
      }
    }

    for (vrcvidxnnd = vrcvidxnum + vrcvdattab[vrcvdsptab[procnum] + 2]; /* Vertices in sub-array 2 */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                finevertlocnum;

      finevertlocnum = vrcvdattab[vrcvidxnum] - finevertlocadj;
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((finevertlocnum < baseval) || (finevertlocnum >= finevertlocnnd)) {
        errorPrint ("bdgraphBipartMlUncoarsen: internal error (6)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      finepartgsttax[finevertlocnum] = 0;
      finefronloctab[finefronlocnbr ++] = finevertlocnum;
    }

    for (vrcvidxnnd = vrcvdsptab[procnum] + vrcvcnttab[procnum]; /* Vertices in sub-array 3 */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                finevertlocnum;

      finevertlocnum = vrcvdattab[vrcvidxnum] - finevertlocadj;
#ifdef SCOTCH_DEBUG_BDGRAPH2
      if ((finevertlocnum < baseval) || (finevertlocnum >= finevertlocnnd)) {
        errorPrint ("bdgraphBipartMlUncoarsen: internal error (7)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
      finepartgsttax[finevertlocnum] = 1;
      finefronloctab[finefronlocnbr ++] = finevertlocnum;
      if (fineveloloctax != NULL)
        finecomplocload1 += fineveloloctax[finevertlocnum];

      if (fineveexloctax != NULL) {
        Gnum                veexval;

        veexval = fineveexloctax[finevertlocnum];
        finecommlocloadextn += veexval; 
        finecommlocgainextn -= veexval;
      }
    }
  }
#ifdef SCOTCH_DEBUG_BDGRAPH2
  for (finevertlocnum = baseval; finevertlocnum < finevertlocnnd; finevertlocnum ++) {
    if (finepartgsttax[finevertlocnum] == ((GraphPart) ~0)) {
      errorPrint ("bdgraphBipartMlUncoarsen: internal error (8)");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  if (dgraphHaloSync (&finegrafptr->s, (byte *) (finepartgsttax + baseval), GRAPHPART_MPI) != 0) {
    errorPrint ("bdgraphBipartMlUncoarsen: cannot perform halo exchange");
    return     (1);
  }
    
  finecommlocloadintn = 0;
  fineedlolocval      = 1;                        /* Assume edges are not weighted */
  for (finefronlocnum = 0; finefronlocnum < finefronlocnbr; finefronlocnum ++) {
    Gnum                finevertlocnum;
    Gnum                fineedgelocnum;
    Gnum                partval;
    Gnum                commcut;

    finevertlocnum = finefronloctab[finefronlocnum];
    partval = finepartgsttax[finevertlocnum];
    for (fineedgelocnum = finevertloctax[finevertlocnum], commcut = 0; 
         fineedgelocnum < finevendloctax[finevertlocnum]; fineedgelocnum ++) {
      Gnum                partdlt;

      partdlt  = partval ^ finepartgsttax[fineedgegsttax[fineedgelocnum]];
      commcut += partdlt;
      if (fineedloloctax != NULL)
        fineedlolocval = fineedloloctax[fineedgelocnum];
      finecommlocloadintn += partdlt * fineedlolocval; /* Counted in both part, should be divided by 2 in summing up phase */
    }
    if (commcut == 0)                             /* If vertex does not really belong to frontier       */
      finefronloctab[finefronlocnum --] = finefronloctab[-- finefronlocnbr]; /* Replace vertex and redo */
  }

  memFree (vrcvdattab);                         /* Free group leaders */
  memFree (vrcvcnttab);

  finegrafptr->fronlocnbr   = finefronlocnbr;
  finegrafptr->complocsize0 = finegrafptr->s.vertlocnbr - finecomplocsize1; 
  finegrafptr->complocload0 = (fineveloloctax == NULL) ? finegrafptr->complocsize0 : (finegrafptr->s.velolocsum - finecomplocload1);

  reduloctab[0] = finegrafptr->complocload0;
  reduloctab[1] = finegrafptr->complocsize0;
  reduloctab[2] = finegrafptr->fronlocnbr;
  reduloctab[3] = finecommlocloadintn;
  reduloctab[4] = finecommlocloadextn;
  reduloctab[5] = finecommlocgainextn;
  MPI_Allreduce (reduloctab, reduglbtab, 6, GNUM_MPI, MPI_SUM, finegrafptr->s.proccomm);

  finegrafptr->compglbload0    = reduglbtab[0];
  finegrafptr->compglbload0dlt = finegrafptr->compglbload0 - finegrafptr->compglbload0avg;
  finegrafptr->compglbsize0    = reduglbtab[1];
  finegrafptr->fronglbnbr      = reduglbtab[2];
  finegrafptr->commglbload     = ((reduglbtab[3] / 2) * finegrafptr->domndist + reduglbtab[4] + finegrafptr->commglbloadextn0);
  finegrafptr->commglbgainextn = reduglbtab[5];
  finegrafptr->bbalglbval      = coargrafptr->bbalglbval;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (finegrafptr) != 0) {
    errorPrint ("bdgraphBipartMlUncoarsen: inconsistent graph data (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}

/* This routine performs the
** bipartitioning recursion.
** It returns:
** - 0   : if bipartitioning could be computed.
** - !0  : on error.
*/

static
int
bdgraphBipartMl2 (
Bdgraph * restrict const            grafptr,      /* Active graph      */
const BdgraphBipartMlParam * const  paraptr)      /* Method parameters */
{
  Bdgraph                       coargrafdat;
  DgraphCoarsenMulti * restrict coarmulttax;
  int                           o;

  if (grafptr->s.procglbnbr <= 1) {               /* Enter into sequential mode          */
    if (((o = bdgraphBipartMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph */
        ((o = bdgraphBipartSt (grafptr, paraptr->stratseq)) != 0)) {
#ifdef SCOTCH_DEBUG_BDGRAPH2
      errorPrintW ("bdgraphBipartMl2: cannot apply sequential strategy");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    }
    return (o);
  }

  coarmulttax = NULL;                             /* Assume multinode array is not allocated */
  if (bdgraphBipartMlCoarsen (grafptr, &coargrafdat, &coarmulttax, paraptr) == 0) {
    o = (coargrafdat.s.procglbnbr == 0) ? 0 : bdgraphBipartMl2 (&coargrafdat, paraptr); /* Apply recursion on coarsened graph if it exists */

    if ((o == 0) &&
        ((o = bdgraphBipartMlUncoarsen (grafptr, &coargrafdat, coarmulttax)) == 0) &&
        ((o = bdgraphBipartSt          (grafptr, paraptr->stratasc)) != 0)) { /* Apply ascending strategy if uncoarsening worked */
#ifdef SCOTCH_DEBUG_BDGRAPH2
      errorPrintW ("bdgraphBipartMl2: cannot apply ascending strategy");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
    }

    bdgraphExit (&coargrafdat);
    if (coarmulttax != NULL)                      /* If multinode array has been allocated */
      memFree (coarmulttax + grafptr->s.baseval); /* Free array                            */

    if (o == 0)                                   /* If multi-level failed, apply low strategy as fallback */
      return (o);
  }

  if (((o = bdgraphBipartMlUncoarsen (grafptr, NULL, NULL)) == 0) && /* Finalize graph            */
      ((o = bdgraphBipartSt          (grafptr, paraptr->stratlow)) != 0)) { /* Apply low strategy */
#ifdef SCOTCH_DEBUG_BDGRAPH2
    errorPrintW ("bdgraphBipartMl2: cannot apply low strategy");
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
  }

  return (o);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the muti-level bipartitioning.
** It returns:
** - 0 : if bipartitioning could be computed.
** - 1 : on error.
*/

int
bdgraphBipartMl (
Bdgraph * const                       grafptr,    /*+ Active graph      +*/
const BdgraphBipartMlParam * const  paraptr)      /*+ Method parameters +*/
{
  Gnum                levlnum;                    /* Save value for graph level */
  int                 o;

  levlnum = grafptr->levlnum;                     /* Save graph level                   */
  grafptr->levlnum = 0;                           /* Initialize coarsening level        */
  o = bdgraphBipartMl2 (grafptr, paraptr);        /* Perform multi-level bipartitioning */
  grafptr->levlnum = levlnum;                     /* Restore graph level                */

  return (o);
}
