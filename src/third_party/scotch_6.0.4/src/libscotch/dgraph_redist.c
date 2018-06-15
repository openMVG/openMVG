/* Copyright 2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_redist.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file implements the distributed    **/
/**                graph redistribution method.            **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 10 may 2010     **/
/**                                 to   : 13 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_REDIST 

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_redist.h"

/**********************************/
/*                                */
/* Graph redistribution routines. */
/*                                */
/**********************************/

/* This routine creates a redistributed
** destination graph by redistributing
** the contents of the given source graph
** according to the provided information.
** It returns:
** - 0   : if the redistributed graph has been created.
** - !0  : on error.
*/

int
dgraphRedist (
Dgraph * restrict const     srcgrafptr,           /* Source distributed graph         */
const Gnum * restrict const srcpartloctax,        /* Array of process destinations    */
const Gnum * restrict const srcpermgsttax,        /* Redistribution permutation array */
const Gnum                  dstvertlocdlt,        /* Extra size of local vertex array */
const Gnum                  dstedgelocdlt,        /* Extra size of local edge array   */
Dgraph * restrict const     dstgrafptr)           /* Destination distributed graph    */
{
  Gnum * restrict       permgsttax;
  const Gnum * restrict permgsttmp;
  Gnum                  permgstnbr;
  Gnum * restrict       procdsptab;
  Gnum * restrict       procvrttab;
  Gnum * restrict       vadjloctab;
  Gnum * restrict       vadjglbtab;
  Gnum                  vadjglbnbr;
  Gnum                  vertlocnum;
  int                   cheklocval;
  int                   chekglbval;
  Gnum                  procdspval;
  Gnum                  procvrtval;
  int                   procglbnbr;
  int                   procnum;
  int                   o;

  if (srcpartloctax == NULL) {
    errorPrint ("dgraphRedist: part array must be provided");
    return     (1);
  }

  cheklocval = 0;
  procglbnbr = srcgrafptr->procglbnbr;
  if (srcpermgsttax != NULL) {                    /* Do not allocate permutation array if already provided */
    permgstnbr =
    vadjglbnbr = 0;
  }
  else {
    if (dgraphGhst (srcgrafptr) != 0) {             /* Compute ghost edge array if not already present */
      errorPrint ("dgraphRedist: cannot compute ghost edge array");
      return     (1);
    }

    permgstnbr = srcgrafptr->vertgstnbr;
    vadjglbnbr = procglbnbr;
  }
  if (memAllocGroup ((void **) (void *)
                     &procvrttab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                     &procdsptab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                     &vadjloctab, (size_t) (procglbnbr       * sizeof (Gnum)),
                     &vadjglbtab, (size_t) (vadjglbnbr       * sizeof (Gnum)),
                     &permgsttax, (size_t) (permgstnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphRedist: out of memory");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DGRAPH2 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  if (chekglbval != 0) {
    if (procvrttab != NULL)
      memFree (procvrttab);
    return (1);
  }

  memSet (vadjloctab, 0, procglbnbr * sizeof (Gnum));

  for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) /* Count number of vertices for each processor */
    vadjloctab[srcpartloctax[vertlocnum]] ++;

  if (MPI_Allreduce (vadjloctab, procdsptab, procglbnbr, GNUM_MPI, MPI_SUM, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist: communication error (2)");
    return     (1);
  }

  for (procnum = 0, procdspval = procvrtval = srcgrafptr->baseval; procnum < procglbnbr; procnum ++) { /* Create vertex range arrays */
    Gnum                      procglbval;

    procglbval = procdsptab[procnum];
    procdsptab[procnum] = procdspval;             /* Build displacement array */
    procdspval += procglbval;
    procvrttab[procnum] = procvrtval;             /* Build vertex index array by adding vertlocdlt */
    procvrtval += procglbval + dstvertlocdlt;
  }
  procdsptab[procnum] = procdspval;               /* Set end of vertex range arrays */
  procvrttab[procnum] = procvrtval;

  if (srcpermgsttax == NULL) {
    permgsttax -= srcgrafptr->baseval;

    if (MPI_Scan (vadjloctab, vadjglbtab, procglbnbr, GNUM_MPI, MPI_SUM, srcgrafptr->proccomm) != MPI_SUCCESS) { /* Compute permutation start indices */
      errorPrint ("dgraphRedist: communication error (3)");
      return     (1);
    }
    for (procnum = 0; procnum < procglbnbr; procnum ++) /* Finalize permutation start indices */
      vadjglbtab[procnum] -= vadjloctab[procnum] - procvrttab[procnum];

    for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) /* Renumber local vertices */
      permgsttax[vertlocnum] = vadjglbtab[srcpartloctax[vertlocnum]] ++;

    if (dgraphHaloSync (srcgrafptr, permgsttax + srcgrafptr->baseval, GNUM_MPI) != 0) {
      errorPrint ("dgraphRedist: cannot compute halo");
      memFree    (procvrttab);                      /* Free group leader */
      return     (1);
    }

    permgsttmp = permgsttax;
  }
  else
    permgsttmp = srcpermgsttax;

  o = dgraphRedist2 (srcgrafptr, srcpartloctax, permgsttmp, procdsptab, procvrttab, 0, dstedgelocdlt, dstgrafptr);

  memFree (procvrttab);                           /* Free group leader */

  return (o);
}

static
int
dgraphRedist2 (
Dgraph * restrict const     srcgrafptr,           /* Source distributed graph           */
const Gnum * restrict const srcpartloctax,        /* Array of process destinations      */
const Gnum * restrict const srcpermgsttax,        /* Redistribution permutation array   */
const Gnum * const          dstprocdsptab,        /* New distribution of graph vertices */
const Gnum * const          dstprocvrttab,        /* New distribution of graph vertices */
const Gnum                  dstvertlocdlt,        /* Extra size of local vertex array   */
const Gnum                  dstedgelocdlt,        /* Extra size of local edge array     */
Dgraph * restrict const     dstgrafptr)           /* Destination distributed graph      */
{
  Gnum                      baseval;
  int                       flveval;              /* Number of data to send per vertex                   */
  int                       fledval;              /* Number of data to send per edge                     */
  Gnum *                    drcvdattab;           /* Receive array for vertex and edge data [norestrict] */
  Gnum *                    dsnddattab;           /* Send array for vertex and edge data [norestrict]    */
  int * restrict            drcvcnttab;           /* Count array for received data                       */
  int * restrict            dsndcnttab;           /* Count array for sent data                           */
  int * restrict            drcvdsptab;           /* Displacement array for received data                */
  int * restrict            dsnddsptab;           /* Displacement array for sent data                    */
  int                       drcvdatnbr;           /* Amount of data to allocate                          */
  int                       dsnddatnbr;
  int                       drcvdatidx;
  int                       dsnddatidx;
  Gnum                      srcvertlocnum;
  Gnum                      srcvertlocnnd;
  Gnum                      srcvertlocadj;
  Gnum * restrict           dstvertloctax;
  Gnum                      dstvertlocadj;
  Gnum                      dstvertlocnbr;
  Gnum                      dstvertlocnnd;
  Gnum                      dstvertlocnum;
  Gnum * restrict           dstveloloctax;
  Gnum                      dstvelolocsiz;
  Gnum                      dstvelolocsum;
  Gnum * restrict           dstvlblloctax;
  Gnum *                    dstedgeloctax;        /* Pointer to destination edge array [norestrict]      */
  Gnum                      dstedgelocnbr;
  Gnum                      dstedgelocsiz;
  Gnum                      dstedgelocnum;
  Gnum *                    dstedloloctax;        /* Pointer to destination edge load array [norestrict] */
  Gnum                      dstedlolocsiz;
  int                       dstvertloctmp;        /* Vertex and edge numbers, as (int)s                  */
  int                       dstedgeloctmp;
  Gnum                      procdspval;
  int                       procglbnbr;
  int                       cheklocval;
  int                       chekglbval;
  int                       procnum;

  const Gnum * restrict const srcvertloctax = srcgrafptr->vertloctax;
  const Gnum * restrict const srcvendloctax = srcgrafptr->vendloctax;
  const Gnum * restrict const srcveloloctax = srcgrafptr->veloloctax;
  const Gnum * restrict const srcvlblloctax = srcgrafptr->vlblloctax;
  const Gnum * restrict const srcedgegsttax = srcgrafptr->edgegsttax;
  const Gnum * restrict const srcedloloctax = srcgrafptr->edloloctax;

  dstgrafptr->flagval |= (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;

  cheklocval = 0;
  procglbnbr = srcgrafptr->procglbnbr;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &dstgrafptr->procdsptab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                     &dstgrafptr->procvrttab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                     &dstgrafptr->proccnttab, (size_t) (procglbnbr       * sizeof (Gnum)),
                     &dstgrafptr->procngbtab, (size_t) (procglbnbr       * sizeof (int)),
                     &dstgrafptr->procrcvtab, (size_t) (procglbnbr       * sizeof (int)),
                     &dstgrafptr->procsndtab, (size_t) (procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphRedist2: out of memory (1)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist2: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DGRAPH2 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  if (chekglbval != 0) {
    dgraphFree (dstgrafptr);
    return     (1);
  }

  dsndcnttab = (int *) dstgrafptr->procdsptab;    /* TRICK: use procdsptab and procvrttab as paired send count arrays */
  drcvcnttab = (int *) dstgrafptr->proccnttab;    /* TRICK: use proccnttab and procngbtab as paired receive arrays    */
  dsnddsptab = (int *) dstgrafptr->procrcvtab;
  drcvdsptab = (int *) dstgrafptr->procsndtab;

  memSet (dsndcnttab, 0, procglbnbr * 2 * sizeof (int)); /* TRICK: Pairs of vertex and edge counts will be exchanged */

  baseval = srcgrafptr->baseval;
  for (srcvertlocnum = baseval, srcvertlocnnd = srcgrafptr->vertlocnnd;
       srcvertlocnum < srcvertlocnnd; srcvertlocnum ++) {
    Gnum                procngbnum;

    procngbnum = srcpartloctax[srcvertlocnum];
    dsndcnttab[2 * procngbnum] ++;                /* One more vertex */
    dsndcnttab[2 * procngbnum + 1] += (int) (srcvendloctax[srcvertlocnum] - srcvertloctax[srcvertlocnum]); /* More edges */
  }

  if (MPI_Alltoall (dsndcnttab, 2, MPI_INT,       /* Get amounts of vertex and edge data to receive */
                    drcvcnttab, 2, MPI_INT, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist2: communication error (2)");
    return     (1);
  }

  fledval = ((srcgrafptr->edloloctax != NULL) ? 1 : 0) + 1; /* Amount of data to exchange per edge  */
  flveval = ((srcgrafptr->veloloctax != NULL) ? 1 : 0) + 3; /* Number, degree and label count for 3 */

  for (procnum = 0, drcvdatidx = dsnddatidx = 0, dstvertloctmp = dstedgeloctmp = 0;
       procnum < procglbnbr; procnum ++) {        /* Compute start indices for data send and receive arrays */
    int                 dsndcntval;
    int                 drcvcntval;

    dsndcntval = dsndcnttab[2 * procnum] * flveval + dsndcnttab[2 * procnum + 1] * fledval;
    drcvcntval = drcvcnttab[2 * procnum] * flveval + drcvcnttab[2 * procnum + 1] * fledval;
    dstvertloctmp += drcvcnttab[2 * procnum];     /* Accumulate number of vertices and edges */
    dstedgeloctmp += drcvcnttab[2 * procnum + 1];

    dsnddsptab[procnum] = dsnddatidx;
    dsnddatidx += dsndcntval;
    drcvdsptab[procnum] = drcvdatidx;
    drcvdatidx += drcvcntval;
  }

  dsnddatnbr = dsnddatidx;                        /* Preserve amount of data to allocate for sending and receiving */
  drcvdatnbr = drcvdatidx;

  for (procnum = procglbnbr - 1; procnum >= 0; procnum --) { /* Compute count arrays for data send and receive arrays */
    int                 dsnddspval;
    int                 drcvdspval;

    dsnddspval = dsnddsptab[procnum];
    dsndcnttab[procnum] = dsnddatidx - dsnddspval;
    dsnddatidx = dsnddspval;
    drcvdspval = drcvdsptab[procnum];
    drcvcnttab[procnum] = drcvdatidx - drcvdspval;
    drcvdatidx = drcvdspval;
  }

  dstvertlocnbr = (Gnum) dstvertloctmp;
  dstedgelocnbr = (Gnum) dstedgeloctmp;
  dstedgelocsiz = dstedgelocnbr + dstedgelocdlt;
  dstvelolocsiz = (srcgrafptr->veloloctax != NULL) ? dstvertlocnbr + dstvertlocdlt : 0;
  dstedlolocsiz = (srcgrafptr->edloloctax != NULL) ? dstedgelocsiz : 0;
  if (memAllocGroup ((void **) (void *)
                     &dstvertloctax, (size_t) ((dstvertlocnbr + dstvertlocdlt + 1) * sizeof (Gnum)), /* Create compact array */
                     &dstveloloctax, (size_t) ( dstvelolocsiz                      * sizeof (Gnum)),
                     &dstvlblloctax, (size_t) ((dstvertlocnbr + dstvertlocdlt)     * sizeof (Gnum)), NULL) == NULL) { /* Vertex labels always present */
    errorPrint ("dgraphRedist2: out of memory (2)");
    cheklocval = 1;
  }
  else if (dstvertloctax -= baseval,
           dstveloloctax  = ((srcgrafptr->veloloctax != NULL) ? dstveloloctax - baseval : NULL),
           dstvlblloctax -= baseval,
           memAllocGroup ((void **) (void *)
                          &dstedgeloctax, (size_t) ((dstedlolocsiz +  /* TRICK: extra space required only if edge loads                                */
                                                     MAX (dstedgelocdlt, srcgrafptr->degrglbmax)) * sizeof (Gnum)), /* TRICK: degrmax to avoid overlap */
                          &dsnddattab,    (size_t) (dsnddatnbr * sizeof (Gnum)), /* TRICK: send space will be edgeloctab or edloloctab                 */
                          &drcvdattab,    (size_t) (drcvdatnbr * sizeof (Gnum)), NULL) == NULL) { /* TRICK: Remaining space will be freed              */
    errorPrint ("dgraphRedist2: out of memory (3)");
    cheklocval = 1;
  }
  else {
    dstedgeloctax -= baseval;
    dstedloloctax  = (srcgrafptr->edloloctax != NULL) ? (dstedgeloctax + dstedgelocsiz) : NULL;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist2: communication error (3)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DGRAPH2 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  if (chekglbval != 0) {
    dgraphFree (dstgrafptr);
    return     (1);
  }

  srcvertlocadj = srcgrafptr->procvrttab[srcgrafptr->proclocnum] - baseval;
  for (srcvertlocnum = baseval; srcvertlocnum < srcvertlocnnd; srcvertlocnum ++) { /* Record data to send */
    Gnum                procngbnum;
    int                 dsnddatidx;
    Gnum                srcedgelocnum;
    Gnum                srcedgelocnnd;
    Gnum                srcdegrval;

    procngbnum = srcpartloctax[srcvertlocnum];    /* Retrieve destination process number */
    dsnddatidx = dsnddsptab[procngbnum];

    srcedgelocnum = srcvertloctax[srcvertlocnum];
    srcedgelocnnd = srcvendloctax[srcvertlocnum];
    srcdegrval    = srcedgelocnnd - srcedgelocnum;

    dsnddattab[dsnddatidx ++] = srcpermgsttax[srcvertlocnum]; /* Record destination vertex global number */
    dsnddattab[dsnddatidx ++] = srcdegrval;       /* Record number of edges                              */
    dsnddattab[dsnddatidx ++] = (srcvlblloctax != NULL) /* Record source vertex global number or label   */
                                ? srcvlblloctax[srcvertlocnum]
                                : srcvertlocnum + srcvertlocadj;
    if (srcveloloctax != NULL)
      dsnddattab[dsnddatidx ++] = srcveloloctax[srcvertlocnum]; /* Record vertex load if needed */

    if (srcedloloctax != NULL) {                  /* If edge loads have to be sent too */
      memCpy (dsnddattab + dsnddatidx, srcedloloctax + srcedgelocnum, srcdegrval * sizeof (Gnum)); /* Copy edge loads */
      dsnddatidx += srcdegrval;
    }

    for ( ; srcedgelocnum < srcedgelocnnd; srcedgelocnum ++) /* Record translated edge array */
      dsnddattab[dsnddatidx ++] = srcpermgsttax[srcedgegsttax[srcedgelocnum]];

    dsnddsptab[procngbnum] = dsnddatidx;
  }

  for (procnum = 0, dsnddatidx = 0;               /* Recompute dsnddsptab */
       procnum < procglbnbr; procnum ++) {
    dsnddsptab[procnum] = dsnddatidx;
    dsnddatidx += dsndcnttab[procnum];
  }

  if (MPI_Alltoallv (dsnddattab, dsndcnttab, dsnddsptab, GNUM_MPI, /* Exchange graph data */
                     drcvdattab, drcvcnttab, drcvdsptab, GNUM_MPI, srcgrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphRedist2: communication error (4)");
    return     (1);
  }

  dstvertlocadj = dstprocvrttab[srcgrafptr->proclocnum] - baseval;
  for (drcvdatidx = 0; drcvdatidx < drcvdatnbr; ) {
    Gnum                dstvertlocnum;
    Gnum                dstdegrval;

    dstvertlocnum = drcvdattab[drcvdatidx ++] - dstvertlocadj; /* Get vertex index */
    dstdegrval    = drcvdattab[drcvdatidx ++];    /* Get number of edges           */

    dstvertloctax[dstvertlocnum] = dstdegrval;    /* Record vertex degree to compute index array */
    dstvlblloctax[dstvertlocnum] = drcvdatidx;    /* TRICK: record data position in label array  */

    drcvdatidx += dstdegrval * fledval + (flveval - 2); /* Increase index by proper value */
  }

  dstvelolocsum = (dstveloloctax != NULL) ? 0 : dstvertlocnbr; /* Set local vertex load sum if no vertex loads present         */
  for (dstvertlocnum = dstedgelocnum = baseval, dstvertlocnnd = dstvertlocnbr + baseval; /* Copy edge information in due place */
       dstvertlocnum < dstvertlocnnd; dstvertlocnum ++) {
    int                 drcvdatidx;
    Gnum                dstdegrval;

    drcvdatidx = dstvlblloctax[dstvertlocnum];
    dstdegrval = dstvertloctax[dstvertlocnum];
    dstvertloctax[dstvertlocnum] = dstedgelocnum;
    dstvlblloctax[dstvertlocnum] = drcvdattab[drcvdatidx ++]; /* Set vertex label */
    if (dstveloloctax != NULL) {
      dstvelolocsum +=
      dstveloloctax[dstvertlocnum] = drcvdattab[drcvdatidx ++]; /* Set vertex load */
    }

    if (dstedloloctax != NULL) {
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (abs ((dstedloloctax + dstedgelocnum) - (drcvdattab + drcvdatidx)) < dstdegrval) { /* Memory areas should never overlap */
        errorPrint ("dgraphRedist2: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      memCpy (dstedloloctax + dstedgelocnum, drcvdattab + drcvdatidx, dstdegrval * sizeof (Gnum));
      drcvdatidx += dstdegrval;
    }
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (abs ((dstedgeloctax + dstedgelocnum) - (drcvdattab + drcvdatidx)) < dstdegrval) { /* TRICK: memory areas should never overlap because of degrmax */
      errorPrint ("dgraphRedist2: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    memCpy (dstedgeloctax + dstedgelocnum, drcvdattab + drcvdatidx, dstdegrval * sizeof (Gnum)); /* TRICK: will never overlap */

    dstedgelocnum += dstdegrval;
  }
  dstvertloctax[dstvertlocnum] = dstedgelocnum;   /* Set end of compact vertex array */

  dstedgeloctax = memRealloc (dstedgeloctax + baseval, dstedgelocsiz * fledval * sizeof (Gnum));
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dstedgeloctax == NULL) {                    /* Shrinking should never fail */
    errorPrint ("dgraphRedist2: out of memory (4)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  dstedgeloctax -= baseval;
  if (dstedloloctax != NULL)
    dstedloloctax = dstedgeloctax + dstedgelocsiz;

  dstgrafptr->procglbnbr = procglbnbr;
  dstgrafptr->proclocnum = srcgrafptr->proclocnum;
  memCpy (dstgrafptr->procvrttab, dstprocvrttab, (procglbnbr + 1) * sizeof (Gnum)); /* Set vertex range array (possibly with holes) */
  memCpy (dstgrafptr->procdsptab, dstprocdsptab, (procglbnbr + 1) * sizeof (Gnum)); /* Set vertex displacement array                */
  for (procnum = procglbnbr - 1, procdspval = dstprocdsptab[procglbnbr]; /* Set vertex count array                                  */
       procnum >= 0; procnum --) {
    Gnum                procdsptmp;

    procdsptmp = dstprocdsptab[procnum];
    dstgrafptr->proccnttab[procnum] = procdspval - procdsptmp;
    procdspval = procdsptmp;
  }

  if (dgraphBuild3 (dstgrafptr, baseval,
                    dstvertlocnbr, dstvertloctax, dstvertloctax + 1, dstveloloctax, dstvelolocsum, NULL, NULL,
                    dstedgelocnbr, dstedgelocsiz, dstedgeloctax, NULL, dstedloloctax, srcgrafptr->degrglbmax) != 0) {
    errorPrint ("dgraphRedist2: cannot build redistributed graph");
    dgraphFree (dstgrafptr);
    return     (1);
  }

  dstgrafptr->vlblloctax = dstvlblloctax;         /* Set label array after building so that labels not taken into account at build time */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (dstgrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("dgraphRedist2: inconsistent graph data");
    dgraphFree (dstgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  return (0);
}
