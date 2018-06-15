/* Copyright 2007-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_coarsen.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : This file implements the coarsening     **/
/**                phase of the multi-level method.        **/
/**                The implementation uses several         **/
/**                processes, which could have several     **/
/**                threads each (3 at this time).          **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 27 jul 2005     **/
/**                                 to   : 15 may 2008     **/
/**                # Version 5.1  : from : 23 jun 2008     **/
/**                                 to   : 20 feb 2011     **/
/**                # Version 6.0  : from : 11 sep 2012     **/
/**                                 to   : 28 sep 2014     **/
/**                                                        **/
/************************************************************/

#define DGRAPH_COARSEN

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_allreduce.h"
#include "dgraph_coarsen.h"
#include "dgraph_match.h"

/******************************/
/*                            */
/* Graph coarsening routines. */
/*                            */
/******************************/

static
int
dgraphCoarsenInit (
DgraphCoarsenData * restrict const  coarptr,      /*+ Coarsening data structure +*/
Dgraph * restrict const             finegrafptr,  /*+ Graph to coarsen          +*/
Dgraph * restrict const             coargrafptr)  /*+ Coarse graph to build     +*/
{
  int                 procglbnbr;
  int                 procglbnum;
  int                 procngbnbr;
  int                 procngbnum;
  int                 procngbnxt;
  int                 vertrcvnbr;
  int                 vertsndnbr;
  Gnum                vertlocnbr;
  Gnum                vertgstnbr;
  int                 vdsprcvnum;
  int                 vdspsndnum;
  byte *              bufftab;
  size_t              buffsiz;

  const int * restrict const  fineprocngbtab = finegrafptr->procngbtab;
  const int * restrict const  fineprocrcvtab = finegrafptr->procrcvtab;
  const int * restrict const  fineprocsndtab = finegrafptr->procsndtab;

  vertlocnbr = finegrafptr->vertlocnbr;
  vertgstnbr = finegrafptr->vertgstnbr;
  procglbnbr = finegrafptr->procglbnbr;
  procngbnbr = finegrafptr->procngbnbr;
  vertrcvnbr = vertgstnbr - vertlocnbr;
  vertsndnbr = finegrafptr->procsndnbr;

  coarptr->coarprvptr = NULL;                     /* Assume nothing to free on error */
  coarptr->multloctmp = NULL;
  coarptr->nsndidxtab = NULL;
  coarptr->nrcvidxtab = NULL;

  if ((coarptr->coarprvptr = memAllocGroup ((void **) (void *) /* Allocate distributed coarse graph private data */
                                            &coargrafptr->procdsptab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                                            &coargrafptr->proccnttab, (size_t) (procglbnbr       * sizeof (Gnum)),
                                            &coargrafptr->procngbtab, (size_t) (procglbnbr       * sizeof (int)),
                                            &coargrafptr->procrcvtab, (size_t) (procglbnbr       * sizeof (int)),
                                            &coargrafptr->procsndtab, (size_t) (procglbnbr       * sizeof (int)), NULL)) == NULL) {
    errorPrint ("dgraphCoarsenInit: out of memory (1)");
    return     (1);
  }
  coargrafptr->procvrttab = coargrafptr->procdsptab; /* Coarse graph has no holes */

  if (coarptr->multloctab == NULL) {              /* If no multinode array provided */
    if ((coarptr->multloctab = memAlloc (vertlocnbr * sizeof (DgraphCoarsenMulti))) == NULL) {
      errorPrint        ("dgraphCoarsenInit: out of memory (2)");
      dgraphCoarsenExit (coarptr);
      return            (1);
    }
    coarptr->multloctmp = coarptr->multloctab;    /* Array will have to be freed on error */
  }

  if (memAllocGroup ((void **) (void *)           /* Data used up to edge exchange phase at coarse graph build time */
                     &coarptr->nrcvidxtab, (size_t) (procngbnbr * sizeof (int)),
                     &coarptr->vrcvdsptab, (size_t) ((procglbnbr + 1) * sizeof (int)), /* TRICK: "+1" for size count */
                     &coarptr->coargsttax, (size_t) (vertgstnbr * sizeof (Gnum)),
                     &coarptr->procgsttax, (size_t) (vertrcvnbr * sizeof (int)), /* TRICK: Only purely ghost part of array will be used */
                     &coarptr->vrcvdattab, (size_t) (vertrcvnbr * sizeof (DgraphCoarsenVert)), NULL) == NULL) {
    errorPrint        ("dgraphCoarsenInit: out of memory (3)");
    dgraphCoarsenExit (coarptr);
    return            (1);
  }

  buffsiz = 2 * MAX ((procngbnbr * sizeof (MPI_Request)), (procglbnbr * sizeof (int)));
  if (memAllocGroup ((void **) (void *)           /* Data released after coarse vertex index exchange phase */
                     &coarptr->nsndidxtab, (size_t) (procngbnbr * sizeof (int)),
                     &coarptr->vsnddsptab, (size_t) ((procglbnbr + 1) * sizeof (int)), /* TRICK: "+1" for size count check */
                     &bufftab,             (size_t) buffsiz,
                     &coarptr->dcntloctab, (size_t) (procglbnbr * sizeof (DgraphCoarsenCount)),
                     &coarptr->dcntglbtab, (size_t) (procglbnbr * sizeof (DgraphCoarsenCount)),
                     &coarptr->vsnddattab, (size_t) (vertsndnbr * sizeof (DgraphCoarsenVert)), NULL) == NULL) {
    errorPrint        ("dgraphCoarsenInit: out of memory (4)");
    dgraphCoarsenExit (coarptr);
    return            (1);
  }
  coarptr->nrcvreqtab = (MPI_Request *) (void *) bufftab; /* TRICK: point-to-point requests and collective arrays share same space */
  coarptr->nsndreqtab = coarptr->nrcvreqtab + procngbnbr;
  coarptr->vrcvcnttab = (int *) (void *) bufftab;
  coarptr->vsndcnttab = coarptr->vrcvcnttab + procglbnbr;

  for (procglbnum = 0, vdsprcvnum = vdspsndnum = 0; /* Build communication index arrays */
       procglbnum < procglbnbr; procglbnum ++) {
    coarptr->vrcvdsptab[procglbnum] = vdsprcvnum;
    coarptr->vsnddsptab[procglbnum] = vdspsndnum;
    vdsprcvnum += fineprocrcvtab[procglbnum];
    vdspsndnum += fineprocsndtab[procglbnum];
  }
  coarptr->vrcvdsptab[procglbnum] = vdsprcvnum;   /* Mark end of communication index arrays */
  coarptr->vsnddsptab[procglbnum] = vdspsndnum;

  for (procngbnum = procngbnxt = 0; procngbnum < procngbnbr; procngbnum ++) {
    if ((procngbnxt == 0) && (fineprocngbtab[procngbnum] > finegrafptr->proclocnum)) { /* Find index of first neighbor of higher rank */
      procngbnxt = procngbnum;
      break;
    }
  }
  coarptr->procngbnxt = procngbnxt;

  coarptr->coargsttax -= finegrafptr->baseval;
  coarptr->finegrafptr = finegrafptr;
  coarptr->coargrafptr = coargrafptr;

  memSet (coarptr->dcntloctab, 0, procglbnbr * sizeof (DgraphCoarsenCount));

  memSet (coarptr->procgsttax, ~0, vertrcvnbr * sizeof (int)); /* Values have not yet been computed                       */
  coarptr->procgsttax -= vertlocnbr + finegrafptr->baseval; /* TRICK: base array such that only purely ghost part is used */

  coarptr->edgekptnbr = 0;

  return (0);
}

static
void
dgraphCoarsenExit (
DgraphCoarsenData * restrict const    coarptr)    /*+ Coarsening data structure +*/
{
  if (coarptr->nsndidxtab != NULL)                /* Auxiliary array is released after first phase of coarse graph building */
    memFree (coarptr->nsndidxtab);
  if (coarptr->nrcvidxtab != NULL)
    memFree (coarptr->nrcvidxtab);
  if (coarptr->multloctmp != NULL)                /* If multinode array not provided nor passed back to calling routine */
    memFree (coarptr->multloctmp);
  if (coarptr->coarprvptr != NULL)                /* If ownership of coarse graph private data not yet transferred to it */
    memFree (coarptr->coarprvptr);
}

static
int
dgraphCoarsenBuildColl (
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          vertlocadj;
  int                           procngbnbr;
  int                           procngbnum;

  MPI_Comm                      proccomm   = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const       grafptr    = coarptr->finegrafptr;
  const int * restrict const    procngbtab = grafptr->procngbtab;
  Gnum * restrict const         coargsttax = coarptr->coargsttax;
  int * restrict const          vsndcnttab = coarptr->vsndcnttab;
  int * restrict const          vrcvdsptab = coarptr->coargrafptr->procrcvtab; /* TRICK: use coarse graph procrcvtab and procsndtab */
  int * restrict const          vsnddsptab = coarptr->coargrafptr->procsndtab;
  int * restrict const          nrcvidxtab = coarptr->nrcvidxtab;
  int * restrict const          nsndidxtab = coarptr->nsndidxtab;

  procngbnbr = grafptr->procngbnbr;
  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;

  memSet (vsndcnttab, 0, grafptr->procglbnbr * sizeof (int));
  memSet (vrcvdsptab, 0, grafptr->procglbnbr * sizeof (int));
  memSet (vsnddsptab, 0, grafptr->procglbnbr * sizeof (int));
  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = procngbtab[procngbnum];
    vsndcnttab[procglbnum] = 2 * (nsndidxtab[procngbnum] - coarptr->vsnddsptab[procglbnum]);
    vrcvdsptab[procglbnum] = 2 * coarptr->vrcvdsptab[procglbnum];
    vsnddsptab[procglbnum] = 2 * coarptr->vsnddsptab[procglbnum];
  }

  if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, coarptr->vrcvcnttab, 1, MPI_INT, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuildColl: communication error (1)");
    return     (1);
  }
  if (MPI_Alltoallv (coarptr->vsnddattab, vsndcnttab,          vsnddsptab, GNUM_MPI,
                     coarptr->vrcvdattab, coarptr->vrcvcnttab, vrcvdsptab, GNUM_MPI, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuildColl: communication error (2)");
    return     (1);
  }

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* For all received data chunks */
    int                 vrcvidxnnd;
    int                 vrcvidxnum;
    int                 procglbnum;
    int                 statsiz;

    const DgraphCoarsenVert * restrict const  vrcvdattab = coarptr->vrcvdattab; /* After data is received */

    procglbnum = procngbtab[procngbnum];
    statsiz = coarptr->vrcvcnttab[procglbnum];
    for (vrcvidxnum = coarptr->vrcvdsptab[procglbnum], vrcvidxnnd = vrcvidxnum + (statsiz / 2); /* TRICK: each message item costs 2 Gnum's */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                vertglbnum;             /* Our global number (the one seen as mate by sender) */
      Gnum                vertlocnum;             /* Our local number (the one seen as mate by sender)  */
      Gnum                multglbnum;             /* Global number of coarse vertex                     */

      vertglbnum = vrcvdattab[vrcvidxnum].datatab[0];
      multglbnum = vrcvdattab[vrcvidxnum].datatab[1];
      vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((vertlocnum <  grafptr->baseval) ||     /* If matching request is not directed towards our process */
          (vertlocnum >= grafptr->vertlocnnd)) {
        errorPrint ("dgraphCoarsenBuildColl: internal error");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      coargsttax[vertlocnum] = multglbnum;
    }
    nrcvidxtab[procngbnum] = vrcvidxnnd;          /* Keep receive end index for preparing edge arrays */
  }

  return (0);
}

static
int
dgraphCoarsenBuildPtop (
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          vertlocadj;
  int                           procngbnbr;
  int                           procngbnum;
  int                           vrcvreqnbr;

  MPI_Comm                      proccomm   = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const       grafptr    = coarptr->finegrafptr;
  const int * restrict const    procngbtab = grafptr->procngbtab;
  Gnum * restrict const         coargsttax = coarptr->coargsttax;
  const int * restrict const    vrcvdsptab = coarptr->vrcvdsptab;
  const int * restrict const    vsnddsptab = coarptr->vsnddsptab;
  int * restrict const          nrcvidxtab = coarptr->nrcvidxtab;
  int * restrict const          nsndidxtab = coarptr->nsndidxtab;

  procngbnbr = grafptr->procngbnbr;
  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;

  if (procngbnbr > 0) {                           /* No communication else             */
    procngbnum = coarptr->procngbnxt;             /* Post receives in descending order */
    do {
      int                 procglbnum;

      procngbnum = (procngbnum + (procngbnbr - 1)) % procngbnbr; /* Pre-decrement neighbor rank */
      procglbnum = procngbtab[procngbnum];
      if (MPI_Irecv (coarptr->vrcvdattab + vrcvdsptab[procglbnum], 2 * (vrcvdsptab[procglbnum + 1] - vrcvdsptab[procglbnum]), GNUM_MPI,
                     procglbnum, TAGCOARSEN, proccomm, &coarptr->nrcvreqtab[procngbnum]) != MPI_SUCCESS) {
        errorPrint ("dgraphCoarsenBuildPtop: communication error (1)");
        return     (1);
      }
    } while (procngbnum != coarptr->procngbnxt);

    procngbnum = coarptr->procngbnxt;             /* Post sends in ascending order */
    do {
      int                 procglbnum;

      procglbnum = procngbtab[procngbnum];
      if (MPI_Isend (coarptr->vsnddattab + vsnddsptab[procglbnum], 2 * (nsndidxtab[procngbnum] - vsnddsptab[procglbnum]), GNUM_MPI,
                     procglbnum, TAGCOARSEN, proccomm, &coarptr->nsndreqtab[procngbnum]) != MPI_SUCCESS) {
        errorPrint ("dgraphCoarsenBuildPtop: communication error (2)");
        return     (1);
      }
      procngbnum = (procngbnum + 1) % procngbnbr; /* Post-increment neighbor rank */
    } while (procngbnum != coarptr->procngbnxt);
  }

  for (vrcvreqnbr = procngbnbr; vrcvreqnbr > 0; vrcvreqnbr --) { /* For all pending receive requests */
    int                 vrcvidxnnd;
    int                 vrcvidxnum;
    int                 procngbnum;
    MPI_Status          statdat;
    int                 statsiz;
    int                 o;

#ifdef SCOTCH_DETERMINISTIC
    procngbnum = vrcvreqnbr - 1;
    o = MPI_Wait (&coarptr->nrcvreqtab[procngbnum], &statdat);
#else /* SCOTCH_DETERMINISTIC */
    o = MPI_Waitany (procngbnbr, coarptr->nrcvreqtab, &procngbnum, &statdat);
#endif /* SCOTCH_DETERMINISTIC */
    if ((o != MPI_SUCCESS) ||
        (MPI_Get_count (&statdat, GNUM_MPI, &statsiz) != MPI_SUCCESS)) {
      errorPrint ("dgraphCoarsenBuildPtop: communication error (3)");
      return     (1);
    }
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (statdat.MPI_SOURCE != procngbtab[procngbnum]) {
      errorPrint ("dgraphCoarsenBuildPtop: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    {
      const DgraphCoarsenVert * restrict const  vrcvdattab = coarptr->vrcvdattab; /* After data is received */

      for (vrcvidxnum = vrcvdsptab[procngbtab[procngbnum]], vrcvidxnnd = vrcvidxnum + (statsiz / 2); /* TRICK: each message item costs 2 Gnum's */
           vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
        Gnum                vertglbnum;           /* Our global number (the one seen as mate by sender) */
        Gnum                vertlocnum;           /* Our local number (the one seen as mate by sender)  */
        Gnum                multglbnum;           /* Global number of coarse vertex                     */

        vertglbnum = vrcvdattab[vrcvidxnum].datatab[0];
        multglbnum = vrcvdattab[vrcvidxnum].datatab[1];
        vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocnum <  grafptr->baseval) ||   /* If matching request is not directed towards our process */
            (vertlocnum >= grafptr->vertlocnnd)) {
          errorPrint ("dgraphCoarsenBuildPtop: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        coargsttax[vertlocnum] = multglbnum;
      }
      nrcvidxtab[procngbnum] = vrcvidxnnd;        /* Keep receive end index for preparing edge arrays */
    }
  }

  if (MPI_Waitall (procngbnbr, coarptr->nsndreqtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait for send requests to complete */
    errorPrint ("dgraphCoarsenBuildPtop: communication error (4)");
    return     (1);
  }

  return (0);
}

/* This routine performs the coarsening of edges
** with respect to the coarmulttax array computed
** by dgraphMatch. All data must be available when
** running (all receptions done). This function is
** inspired by libscotch/src/graph_coarsen_edge.c.
*/

DGRAPHALLREDUCEMAXSUMOP (3, 1)

static
int
dgraphCoarsenBuild (
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          vertlocnum;
  Gnum                          vertlocadj;
  Gnum                          edgelocnbr;
  Gnum                          edlolocval;
  Gnum * restrict               ercvdattab;
  Gnum * restrict               esnddattab;
  int * restrict                ercvcnttab;
  int * restrict                esndcnttab;
  int * restrict                ercvdsptab;
  int * restrict                esnddsptab;
  int                           ercvdspval;
  int                           esnddspval;
  int                           ercvdatsiz;
  int                           esnddatsiz;
  DgraphCoarsenMulti * restrict multloctax;
  Gnum                          multlocnum;
  Gnum                          multlocadj;
  int                           procngbnbr;
  int                           procngbnum;
  int                           procnum;
  Gnum                          coarvertglbnum;
  Gnum                          coarvertlocnum;
  Gnum                          coarvertlocnnd;
  Gnum * restrict               coarvertloctax;
  Gnum * restrict               coarveloloctax;
  Gnum                          coarvelolocsum;
  Gnum                          coardegrlocmax;
  Gnum                          coaredgelocnum;
  Gnum * restrict               coaredgeloctax;
  Gnum * restrict               coaredloloctax;
  Gnum                          coarhashnbr;      /* Size of hash table                 */
  Gnum                          coarhashmsk;      /* Mask for access hash table         */
  DgraphCoarsenHash * restrict  coarhashtab;      /* Table of edges to other multinodes */
  Gnum                          reduloctab[4];
  Gnum                          reduglbtab[4];
  int                           cheklocval;
  int                           chekglbval;
#ifdef SCOTCH_DEBUG_DGRAPH2
  int * restrict                ercvdbgtab;
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  MPI_Comm                                  proccomm    = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const                   grafptr     = coarptr->finegrafptr;
  Dgraph * restrict const                   coargrafptr = coarptr->coargrafptr;
  Gnum * restrict const                     coargsttax  = coarptr->coargsttax;
  const int * restrict const                procngbtab  = grafptr->procngbtab;
  const int * restrict const                procgsttax  = coarptr->procgsttax;
  const Gnum * restrict const               vertloctax  = grafptr->vertloctax;
  const Gnum * restrict const               vendloctax  = grafptr->vendloctax;
  const Gnum * restrict const               veloloctax  = grafptr->veloloctax;
  const Gnum * restrict const               edgeloctax  = grafptr->edgeloctax;
  const Gnum * restrict const               edgegsttax  = grafptr->edgegsttax;
  const Gnum * restrict const               edloloctax  = grafptr->edloloctax;
  const DgraphCoarsenMulti * restrict const multloctab  = coarptr->multloctab;
  DgraphCoarsenVert * const                 vrcvdattab  = coarptr->vrcvdattab; /* [norestrict:async] */
  DgraphCoarsenVert * restrict const        vsnddattab  = coarptr->vsnddattab;
  int * restrict const                      nsndidxtab  = coarptr->nsndidxtab;

#ifdef SCOTCH_DEBUG_DGRAPH2
  memSet (coargsttax + grafptr->baseval, ~0, grafptr->vertgstnbr * sizeof (Gnum));
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  procngbnbr = grafptr->procngbnbr;

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) /* Reset indices for sending messages */
    nsndidxtab[procngbnum] = coarptr->vsnddsptab[procngbtab[procngbnum]];

  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;
  multlocadj = coarptr->coargrafptr->procdsptab[grafptr->proclocnum];
  for (multlocnum = 0; multlocnum < coarptr->multlocnbr; multlocnum ++) {
    Gnum                vertlocnum0;
    Gnum                vertlocnum1;

    vertlocnum0 = multloctab[multlocnum].vertglbnum[0] - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((vertlocnum0 <  grafptr->baseval) ||
        (vertlocnum0 >= grafptr->vertlocnnd)) {
      errorPrint ("dgraphCoarsenBuild: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    coargsttax[vertlocnum0] = multlocnum + multlocadj;

    vertlocnum1 = multloctab[multlocnum].vertglbnum[1];
    if (vertlocnum1 >= 0) {                       /* If second vertex is local */
      vertlocnum1 -= vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((vertlocnum1 <  grafptr->baseval) ||
          (vertlocnum1 >= grafptr->vertlocnnd)) {
        errorPrint ("dgraphCoarsenBuild: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      coargsttax[vertlocnum1] = multlocnum + multlocadj; /* Don't care if single multinode */
    }
    else {
      Gnum                edgelocnum;
      Gnum                vertglbnum1;
      Gnum                vertgstnum1;
      int                 coarsndidx;
      int                 procngbnum;

      edgelocnum = -2 - vertlocnum1;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((edgelocnum < grafptr->baseval) ||
          (edgelocnum >= (grafptr->edgelocsiz + grafptr->baseval))) {
        errorPrint ("dgraphCoarsenBuild: internal error (3)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      vertglbnum1 = edgeloctax[edgelocnum];
      vertgstnum1 = edgegsttax[edgelocnum];

      procngbnum = procgsttax[vertgstnum1];       /* Find neighbor owner process       */
      if (procngbnum < 0) {                       /* If neighbor had not been computed */
        errorPrint ("dgraphCoarsenBuild: internal error (4)");
        return     (1);
      }

      coarsndidx = nsndidxtab[procngbnum] ++;     /* Get position of message in send array */
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (coarsndidx >= coarptr->vsnddsptab[procngbtab[procngbnum] + 1]) {
        errorPrint ("dgraphCoarsenBuild: internal error (5)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      vsnddattab[coarsndidx].datatab[0] = vertglbnum1;
      vsnddattab[coarsndidx].datatab[1] = multlocnum + multlocadj; /* Send multinode value */
    }
  }

  if (coarptr->multloctmp != NULL) {              /* If we allocated the multinode array */
    coarptr->multloctmp =
    coarptr->multloctab = memRealloc (coarptr->multloctab, coarptr->multlocnbr * sizeof (DgraphCoarsenMulti)); /* In the mean time, resize multinode array */
  }

  if ((((grafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphCoarsenBuildPtop : dgraphCoarsenBuildColl) (coarptr) != 0)
    return (1);

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  ercvcnttab = coarptr->coargrafptr->procrcvtab;  /* TRICK: re-use some private coarse graph arrays after vertex exchange phase */
  ercvdsptab = coarptr->coargrafptr->procsndtab;
  for (procnum = 0, ercvdspval = 0; procnum < grafptr->procglbnbr; procnum ++) { /* TRICK: dcntglbtab array no longer needed afterwards; can be freed */
    ercvdsptab[procnum] = ercvdspval;
    ercvcnttab[procnum] = coarptr->dcntglbtab[procnum].vertsndnbr * ((veloloctax != NULL) ? 2 : 1) +
                          coarptr->dcntglbtab[procnum].edgesndnbr * ((edloloctax != NULL) ? 2 : 1);
    ercvdspval += ercvcnttab[procnum];
  }

  memFree (coarptr->nsndidxtab);                  /* Free now useless work memory  */
  coarptr->nsndidxtab = NULL;                     /* This block won't be reclaimed */

#ifdef SCOTCH_DEBUG_DGRAPH2
  for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) {
    if (coargsttax[vertlocnum] < 0) {
      errorPrint ("dgraphCoarsenBuild: invalid matching");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (dgraphHaloSync (grafptr, coargsttax + grafptr->baseval, GNUM_MPI) != 0) {
    errorPrint ("dgraphCoarsenBuild: cannot propagate multinode indices");
    return     (1);
  }

  edgelocnbr = coarptr->edgekptnbr + coarptr->edgercvnbr; /* Upper bound on number of edges     */
  ercvdatsiz = coarptr->vertrcvnbr + coarptr->edgercvnbr; /* Basic size: degrees plus edge data */
  esnddatsiz = coarptr->vertsndnbr + coarptr->edgesndnbr;
  if (grafptr->veloloctax != NULL) {              /* Add vertex loads if necessary */
    ercvdatsiz += coarptr->vertrcvnbr;
    esnddatsiz += coarptr->vertsndnbr;
  }
  if (grafptr->edloloctax != NULL) {              /* Add edge loads if necessary */
    ercvdatsiz += coarptr->edgercvnbr;
    esnddatsiz += coarptr->edgesndnbr;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (ercvdspval != ercvdatsiz) {
    errorPrint ("dgraphCoarsenBuild: internal error (6)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (coarhashmsk = 31; coarhashmsk < grafptr->degrglbmax; coarhashmsk = coarhashmsk * 2 + 1) ;
  coarhashmsk = coarhashmsk * 4 + 3;
  coarhashnbr = coarhashmsk + 1;

  cheklocval = 0;
  coargrafptr->flagval = DGRAPHFREETABS | DGRAPHFREEPRIV | DGRAPHVERTGROUP; /* Coarse graph is not yet based */
  coarptr->coarprvptr = NULL;                     /* Transfer ownership of private arrays to coarse graph    */
  if (memAllocGroup ((void **) (void *)
                     &coargrafptr->vertloctax, (size_t) ((coarptr->multlocnbr + 1) * sizeof (Gnum)),
                     &coargrafptr->veloloctax, (size_t) ( coarptr->multlocnbr      * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphCoarsenBuild: out of memory (1)");
    cheklocval = 1;
  }
  else if ((coargrafptr->edgeloctax = memAlloc (edgelocnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphCoarsenBuild: out of memory (2)");
    cheklocval = 1;
  }
  else if ((coargrafptr->edloloctax = memAlloc (edgelocnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphCoarsenBuild: out of memory (3)");
    cheklocval = 1;
  }
  else if ((coarptr->nsndidxtab = memAllocGroup ((void **) (void *) /* TRICK: allow data array to be released on error */
                                                 &esndcnttab,  (size_t) (grafptr->procglbnbr * sizeof (int)),
                                                 &esnddsptab,  (size_t) (grafptr->procglbnbr * sizeof (int)),
                                                 &esnddattab,  (size_t) (esnddatsiz * sizeof (Gnum)),
                                                 &ercvdattab,  (size_t) (ercvdatsiz * sizeof (Gnum)),
#ifdef SCOTCH_DEBUG_DGRAPH2
                                                 &ercvdbgtab,  (size_t) (grafptr->procglbnbr * sizeof (int)),
#endif /* SCOTCH_DEBUG_DGRAPH2 */
                                                 &coarhashtab, (size_t) (coarhashnbr * sizeof (DgraphCoarsenHash)), NULL)) == NULL) {
    errorPrint ("dgraphCoarsenBuild: out of memory (4)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be overlapped by a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (2)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0) {
    dgraphFree (coargrafptr);
    return     (1);
  }

  memSet (coarhashtab, ~0, coarhashnbr * sizeof (DgraphCoarsenHash));

  coargrafptr->baseval     = grafptr->baseval;
  coargrafptr->vertlocnnd  = coargrafptr->baseval + coargrafptr->vertlocnbr;
  coargrafptr->vertloctax -= coargrafptr->baseval;
  coargrafptr->vendloctax  = coargrafptr->vertloctax + 1; /* Graph is compact */
  coargrafptr->veloloctax -= coargrafptr->baseval;
  coargrafptr->edgeloctax -= coargrafptr->baseval;
  coargrafptr->edloloctax -= coargrafptr->baseval;

  for (procngbnum = procnum = 0, esnddspval = 0; procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;
    int                 vrcvidxnnd;
    int                 vrcvidxnum;

    procglbnum = procngbtab[procngbnum];
    while (procnum < procglbnum) {                /* Fill empty slots */
      esnddsptab[procnum] = esnddspval;
      esndcnttab[procnum] = 0;
      procnum ++;
    }
    esnddsptab[procnum] = esnddspval;

    for (vrcvidxnum = coarptr->vrcvdsptab[procglbnum], vrcvidxnnd = coarptr->nrcvidxtab[procngbnum]; /* For all multinode requests received, in order */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                vertlocnum;
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      vertlocnum = vrcvdattab[vrcvidxnum].datatab[0] - vertlocadj;
      edgelocnum = vertloctax[vertlocnum];
      edgelocnnd = vendloctax[vertlocnum];
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((esnddspval + ((veloloctax != NULL) ? 2 : 1) + ((edloloctax != NULL) ? 2 : 1) * (edgelocnnd - edgelocnum)) > esnddatsiz) {
        errorPrint ("dgraphCoarsenBuild: internal error (7)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      esnddattab[esnddspval ++] = (edgelocnnd - edgelocnum); /* Write degree */
      if (veloloctax != NULL)
        esnddattab[esnddspval ++] = veloloctax[vertlocnum];
      if (edloloctax != NULL) {
        for ( ; edgelocnum < edgelocnnd; edgelocnum ++) {
          esnddattab[esnddspval ++] = coargsttax[edgegsttax[edgelocnum]];
          esnddattab[esnddspval ++] = edloloctax[edgelocnum];
        }
      }
      else {
        for ( ; edgelocnum < edgelocnnd; edgelocnum ++)
          esnddattab[esnddspval ++] = coargsttax[edgegsttax[edgelocnum]];
      }
    }
    esndcnttab[procnum] = esnddspval - esnddsptab[procnum];
    procnum ++;
  }
  while (procnum < grafptr->procglbnbr) {         /* Complete fill-in of empty slots */
    esnddsptab[procnum] = esnddspval;
    esndcnttab[procnum] = 0;
    procnum ++;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (esnddspval != esnddatsiz) {
    errorPrint ("dgraphCoarsenBuild: internal error (8)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  while (procnum < grafptr->procglbnbr) {         /* Complete edge data send displacement array */
    esnddsptab[procnum] = esnddspval;
    esndcnttab[procnum] = 0;
    procnum ++;
  }

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Alltoall (esndcnttab, 1, MPI_INT, ercvdbgtab, 1, MPI_INT, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (3)");
    return     (1);
  }
  for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
    if (ercvdbgtab[procnum] != ercvcnttab[procnum]) {
      errorPrint ("dgraphCoarsenBuild: internal error (9)");
      return     (1);
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  if (MPI_Alltoallv (esnddattab, esndcnttab, esnddsptab, GNUM_MPI,
                     ercvdattab, ercvcnttab, ercvdsptab, GNUM_MPI, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (4)");
    return     (1);
  }

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++)
    ercvdsptab[procngbnum] = ercvdsptab[procngbtab[procngbnum]];

  multloctax = coarptr->multloctab - grafptr->baseval;

  edlolocval = 1;
  coarvelolocsum = 0;
  coardegrlocmax = 0;
  coarvertloctax = coargrafptr->vertloctax;
  coarveloloctax = coargrafptr->veloloctax;
  coaredgeloctax = coargrafptr->edgeloctax;
  coaredloloctax = coargrafptr->edloloctax;
  for (coarvertlocnum = coaredgelocnum = grafptr->baseval, coarvertglbnum = multlocadj, coarvertlocnnd = coarvertlocnum + coargrafptr->vertlocnbr;
       coarvertlocnum < coarvertlocnnd; coarvertlocnum ++, coarvertglbnum ++) {
    Gnum                coarvelolocval;
    Gnum                vertlocnum;
    int                 i;

    coarvertloctax[coarvertlocnum] = coaredgelocnum;

    i = 0;
    coarvelolocval = 0;
    vertlocnum = multloctax[coarvertlocnum].vertglbnum[0] - vertlocadj;
    while (1) {                                   /* Pseudo-infinite loop on both vertices of the multinode */
      Gnum                  vertglbnum;
      Gnum                  edgelocnum;
      Gnum                  edgelocnnd;
      Gnum                  degrlocval;
      int                   procngbnum;
      int                   ercvidxnum;

      coarvelolocval += (veloloctax != NULL) ? veloloctax[vertlocnum] : 1;
      for (edgelocnum = vertloctax[vertlocnum], edgelocnnd = vendloctax[vertlocnum]; /* Loop on edges of first (and sometimes second) local mate */
           edgelocnum < edgelocnnd; edgelocnum ++) {
        Gnum                  coarvertglbend;
        Gnum                  h;

        coarvertglbend = coargsttax[edgegsttax[edgelocnum]];
        if (coarvertglbend == coarvertglbnum)     /* If end of collapsed edge */
          continue;

        if (edloloctax != NULL)
          edlolocval = edloloctax[edgelocnum];
        for (h = (coarvertglbend * COARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
          if (coarhashtab[h].vertorgnum != coarvertglbnum) { /* If old slot           */
            coarhashtab[h].vertorgnum = coarvertglbnum; /* Mark it in reference array */
            coarhashtab[h].vertendnum = coarvertglbend;
            coarhashtab[h].edgelocnum = coaredgelocnum;
#ifdef SCOTCH_DEBUG_DGRAPH2
            if (coaredgelocnum >= (edgelocnbr + coargrafptr->baseval)) {
              errorPrint ("dgraphCoarsenBuild: internal error (10)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
            coaredgeloctax[coaredgelocnum] = coarvertglbend; /* One more edge created */
            coaredloloctax[coaredgelocnum] = edlolocval;
            coaredgelocnum ++;
            break;                                /* Give up hashing */
          }
          if (coarhashtab[h].vertendnum == coarvertglbend) { /* If coarse edge already exists */
            coaredloloctax[coarhashtab[h].edgelocnum] += edlolocval;
            break;                                /* Give up hashing */
          }
        }
      }

      if (i ++ > 0)                               /* If second local vertex has been processed, exit */
        break;

      vertglbnum = multloctax[coarvertlocnum].vertglbnum[1];

      if (vertglbnum >= 0) {                      /* If second multinode vertex is local */
        if ((vertglbnum - vertlocadj) == vertlocnum) /* If single multinode              */
          break;
        vertlocnum = (vertglbnum - vertlocadj);
        continue;
      }

      edgelocnum = -2 - vertglbnum;
      multloctax[coarvertlocnum].vertglbnum[1] = edgeloctax[edgelocnum]; /* Set second vertex of multinode */
      procngbnum = procgsttax[edgegsttax[edgelocnum]];
      ercvidxnum = ercvdsptab[procngbnum];
      degrlocval = ercvdattab[ercvidxnum ++];
      coarvelolocval += (veloloctax != NULL) ? ercvdattab[ercvidxnum ++] : 1;

      while (degrlocval -- > 0) {
        Gnum                  coarvertglbend;
        Gnum                  h;

        coarvertglbend = ercvdattab[ercvidxnum ++];
        if (edloloctax != NULL)
          edlolocval = ercvdattab[ercvidxnum ++];
        if (coarvertglbend == coarvertglbnum)     /* If end of collapsed edge */
          continue;

        for (h = (coarvertglbend * COARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
          if (coarhashtab[h].vertorgnum != coarvertglbnum) { /* If old slot           */
            coarhashtab[h].vertorgnum = coarvertglbnum; /* Mark it in reference array */
            coarhashtab[h].vertendnum = coarvertglbend;
            coarhashtab[h].edgelocnum = coaredgelocnum;
#ifdef SCOTCH_DEBUG_DGRAPH2
            if (coaredgelocnum >= (edgelocnbr + coargrafptr->baseval)) {
              errorPrint ("dgraphCoarsenBuild: internal error (11)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
            coaredgeloctax[coaredgelocnum] = coarvertglbend; /* One more edge created */
            coaredloloctax[coaredgelocnum] = edlolocval;
            coaredgelocnum ++;
            break;                                /* Give up hashing */
          }
          if (coarhashtab[h].vertendnum == coarvertglbend) { /* If coarse edge already exists */
            coargrafptr->edloloctax[coarhashtab[h].edgelocnum] += edlolocval;
            break;                                /* Give up hashing */
          }
        }
      }

      ercvdsptab[procngbnum] = ercvidxnum;        /* Write back updated receive index       */
      break;                                      /* Exit loop after processing remote mate */
    }
    coarvelolocsum += coarvelolocval;
    coarveloloctax[coarvertlocnum] = coarvelolocval;
    if (coardegrlocmax < (coaredgelocnum - coarvertloctax[coarvertlocnum]))
      coardegrlocmax = (coaredgelocnum - coarvertloctax[coarvertlocnum]);
  }
  coarvertloctax[coarvertlocnum] = coaredgelocnum; /* Set end of compact edge array */
  coargrafptr->velolocsum = coarvelolocsum;
  coargrafptr->veloglbsum = grafptr->veloglbsum;
  coargrafptr->edgelocnbr =
  coargrafptr->edgelocsiz = coaredgelocnum - coargrafptr->baseval;

  coargrafptr->edgeloctax  = memRealloc (coaredgeloctax + coargrafptr->baseval, coargrafptr->edgelocnbr * sizeof (Gnum));
  coargrafptr->edgeloctax -= coargrafptr->baseval;
  coargrafptr->edloloctax  = memRealloc (coaredloloctax + coargrafptr->baseval, coargrafptr->edgelocnbr * sizeof (Gnum));
  coargrafptr->edloloctax -= coargrafptr->baseval;
 
  reduloctab[0] = coargrafptr->vertlocnbr;        /* Get maximum over all processes */
  reduloctab[1] = coargrafptr->edgelocnbr;
  reduloctab[2] = coardegrlocmax;                 /* Get local maximum degree */
  reduloctab[3] = coargrafptr->edgelocnbr;

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 3, 1, proccomm) != 0) {
    errorPrint ("dgraphCoarsenBuild: communication error (5)");
    return     (1);
  }

  coargrafptr->vertglbmax = reduglbtab[0];
  coargrafptr->edgeglbmax = reduglbtab[1];
  coargrafptr->degrglbmax = reduglbtab[2];        /* It is now a real global maximum degree */
  coargrafptr->edgeglbnbr = reduglbtab[3];
  coargrafptr->edgeglbsmx = coargrafptr->edgeglbmax;

  return (0);
}

/***************************/
/*                         */
/* The coarsening routine. */
/*                         */
/***************************/

/* This routine coarsens the given fine distributed
** graph, as long as the coarsening ratio remains
** below some threshold value and the coarsened graph
** is not too small.
** If a multinode array is provided (*multlocptr != NULL),
** it must be of a size sufficient to hold multinode data
** in any configuration, including in the case of folding
** with duplication, where folded data is spread across
** floor(P/2) processes.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

int
dgraphCoarsen (
Dgraph * restrict const               finegrafptr, /*+ Graph to coarsen                   +*/
Dgraph * restrict const               coargrafptr, /*+ Coarse graph to build              +*/
DgraphCoarsenMulti * restrict * const multlocptr, /*+ Pointer to un-based multinode array +*/
const Gnum                            passnbr,    /*+ Number of coarsening passes to go   +*/
const Gnum                            coarnbr,    /*+ Minimum number of coarse vertices   +*/
const double                          coarrat,    /*+ Maximum contraction ratio           +*/
const int                             flagval)    /*+ Flag value                          +*/
{
  DgraphMatchData           matedat;              /* Matching state data; includes coarsening handling data   */
  Gnum                      vertrcvnbr;           /* Overall number of vertices to be received from neighbors */
  Gnum                      edgercvnbr;           /* Overall number of edges to be received from neighbors    */
  Gnum                      vertsndnbr;           /* Overall number of vertices to be sent to neighbors       */
  Gnum                      edgesndnbr;           /* Overall number of edges to be sent to neighbors          */
  int                       cheklocval;
  int                       chekglbval;
  Gnum                      coarvertmax;
  Gnum                      passnum;
  int                       procnum;
  int                       o;

#ifdef SCOTCH_DEBUG_DGRAPH1
  if (coarrat < 0.5L)                             /* If impossible coarsening ratio wanted */
    return (1);                                   /* We will never succeed                 */
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  coarvertmax = (Gnum) ((double) finegrafptr->vertglbnbr * coarrat); /* Maximum number of coarse vertices */
  if (coarvertmax < coarnbr)                      /* If there are too few vertices in graph               */
    return (1);                                   /* It is useless to go any further                      */

  if (dgraphGhst (finegrafptr) != 0) {            /* Compute ghost edge array of fine graph if not already present */
    errorPrint ("dgraphCoarsen: cannot compute ghost edge array");
    return     (2);
  }

  matedat.c.flagval    = flagval;
  matedat.c.multloctab = *multlocptr;             /* Propagate the provided multinode array or NULL if it has to be allocated */
  cheklocval  = dgraphCoarsenInit (&matedat.c, finegrafptr, coargrafptr);
  cheklocval |= dgraphMatchInit   (&matedat, 0.5F);

#ifdef SCOTCH_DEBUG_DGRAPH1                       /* This communication cannot be covered by a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsen: communication error (1)");
    return     (2);
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0)
    return (2);

  for (passnum = 0; passnum < passnbr; passnum ++) {
    ((passnum == 0) ? dgraphMatchHl : dgraphMatchHy) (&matedat); /* If first pass, process lightest vertices first */

    if ((((finegrafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphMatchSyncPtop : dgraphMatchSyncColl) (&matedat) != 0) {
      errorPrint        ("dgraphCoarsen: cannot perform matching");
      dgraphMatchExit   (&matedat);
      dgraphCoarsenExit (&matedat.c);
      return            (2);
    }
  }
  dgraphMatchLy (&matedat);                       /* All remaining vertices are matched locally */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphMatchCheck (&matedat) != 0) {
    errorPrint        ("dgraphCoarsen: invalid matching");
    dgraphMatchExit   (&matedat);
    dgraphCoarsenExit (&matedat.c);
    return            (2);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  dgraphMatchExit (&matedat);

  vertsndnbr = 
  edgesndnbr = 0;
  for (procnum = 0; procnum < finegrafptr->procglbnbr; procnum ++) {
    vertsndnbr += matedat.c.dcntloctab[procnum].vertsndnbr;
    edgesndnbr += matedat.c.dcntloctab[procnum].edgesndnbr;
    matedat.c.dcntloctab[procnum].vertlocnbr = matedat.c.multlocnbr;
  }
  matedat.c.vertsndnbr = vertsndnbr;
  matedat.c.edgesndnbr = edgesndnbr;

  if (MPI_Alltoall (matedat.c.dcntloctab, 3, GNUM_MPI, matedat.c.dcntglbtab, 3, GNUM_MPI, finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsen: communication error (2)");
    return     (2);
  }

  vertrcvnbr = 
  edgercvnbr = 0;
  coargrafptr->procdsptab[0] = finegrafptr->baseval; /* Build vertex-to-process array */
  for (procnum = 0; procnum < finegrafptr->procglbnbr; procnum ++) {
    Gnum                proccntval;

    vertrcvnbr += matedat.c.dcntglbtab[procnum].vertsndnbr;
    edgercvnbr += matedat.c.dcntglbtab[procnum].edgesndnbr;
    proccntval  = matedat.c.dcntglbtab[procnum].vertlocnbr;
    coargrafptr->proccnttab[procnum] = proccntval;
    coargrafptr->procdsptab[procnum + 1] = coargrafptr->procdsptab[procnum] + proccntval;
  }
  coargrafptr->vertlocnbr = matedat.c.multlocnbr;
  coargrafptr->vertglbnbr = coargrafptr->procdsptab[finegrafptr->procglbnbr] - finegrafptr->baseval;
  matedat.c.vertrcvnbr = vertrcvnbr;
  matedat.c.edgercvnbr = edgercvnbr;

  if (coargrafptr->vertglbnbr > coarvertmax) {    /* If coarsening ratio not met */
    dgraphCoarsenExit (&matedat.c);
    return            (1);
  }

  if (dgraphCoarsenBuild (&matedat.c) != 0) {     /* Build coarse graph */
    dgraphCoarsenExit (&matedat.c);
    return            (2);
  }

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (coargrafptr) != 0) {           /* Check graph consistency */
    errorPrint ("dgraphCoarsen: inconsistent graph data");
    dgraphFree (coargrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  matedat.c.multloctmp = NULL;                    /* So that it will not be freed    */
  dgraphCoarsenExit (&matedat.c);                 /* Free all other temporary arrays */

  o = 0;                                          /* Assume everything is now all right */
#ifdef PTSCOTCH_FOLD_DUP
  if (((flagval & DGRAPHCOARSENFOLDDUP) != 0) &&  /* If some form of folding is requested */
      (coargrafptr->procglbnbr >= 2)) {           /* And if there is need to it           */
    Dgraph                coargrafdat;            /* Coarse graph data before folding     */
    DgraphCoarsenMulti *  coarmultptr;            /* Pointer to folded multinode array    */
    MPI_Datatype          coarmultype;

    MPI_Type_contiguous (2, GNUM_MPI, &coarmultype); /* Define type for MPI transfer */
    MPI_Type_commit (&coarmultype);               /* Commit new type                 */

    coargrafdat = *coargrafptr;                   /* Copy unfolded coarse graph data to save area */
    coarmultptr = NULL;                           /* Assume we will not get a multinode array     */
    if ((flagval & DGRAPHCOARSENFOLDDUP) == DGRAPHCOARSENFOLD) { /* Do a simple folding           */
      memSet (coargrafptr, 0, sizeof (Dgraph));   /* Also reset procglbnbr for unused processes   */
      o = dgraphFold (&coargrafdat, 0, coargrafptr, (void *) matedat.c.multloctab, (void **) (void *) &coarmultptr, coarmultype);
    }
    else {                                        /* Do a duplicant-folding */
      int               loopval;
      int               dumyval;

      o = dgraphFoldDup (&coargrafdat, coargrafptr, (void *) matedat.c.multloctab, (void **) (void *) &coarmultptr, coarmultype);
      loopval = intRandVal (finegrafptr->proclocnum + intRandVal (finegrafptr->proclocnum * 2 + 1) + 1);
      while (loopval --)                          /* Desynchronize pseudo-random generator across processes */
        dumyval = intRandVal (2);
    }
    dgraphExit    (&coargrafdat);                 /* Free unfolded graph */
    MPI_Type_free (&coarmultype);
    if (*multlocptr == NULL)                      /* If unfolded multinode array was not user-provided, free it */
      memFree (matedat.c.multloctab);
    *multlocptr = coarmultptr;                    /* Return folded multinode array or NULL */
  }
  else                                            /* No folding at all */
#endif /* PTSCOTCH_FOLD_DUP */
    *multlocptr = matedat.c.multloctab;           /* Return un-based pointer (maybe the same as initially user-provided) */

  return (o);
}
