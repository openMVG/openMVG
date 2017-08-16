/* Copyright 2007-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_band_grow.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module flags vertices according    **/
/**                to a breadth-first search traversal of  **/
/**                the distributed graph. It is used both  **/
/**                by dgraphBand() and dgraphGrow().       **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 20 feb 2011     **/
/**                # Version 6.0  : from : 03 apr 2012     **/
/**                                 to   : 26 sep 2012     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  vdgraph_separate_bd.c in version      **/
/**                  5.0. It was first moved to            **/
/**                  dgraph_band.c, then to here to be     **/
/**                  mutualized between dgraphBand and     **/
/**                  dgraphGrow.                           **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPHBANDGROWNSTR          STRINGIFY (DGRAPHBANDGROWNAME)

/**********************************/
/*                                */
/* Distance computation routines. */
/*                                */
/**********************************/

/* This routine computes a distributed index array
** of given width around the current separator,
** using collective communication.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
DGRAPHBANDGROWNAMECOLL (
Dgraph * restrict const           grafptr,        /*+ Distributed graph                                        +*/
const Gnum                        queulocnbr,     /*+ Number of frontier vertices, start size for vertex queue +*/
Gnum * restrict const             queuloctab,     /*+ Array of frontier vertices, re-used as queue array       +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices                 +*/
Gnum * restrict const             vnumgsttax,     /*+ Flag or index array to fill                              +*/
Gnum * restrict const             bandvertlvlptr, /*+ Pointer to based start index of last level               +*/
Gnum * restrict const             bandvertlocptr, /*+ Pointer to bandvertlocnnd                                +*/
Gnum * restrict const             bandedgelocptr) /*+ Pointer to bandedgelocnbr                                +*/
{
  Gnum                    queulocnum;
  Gnum                    vertlocnnd;
  Gnum                    vnumgstsiz;             /* Size of vnumgsttax; TRICK: re-use     */
  Gnum                    vrcvdatsiz;             /* Sizes of data send and receive arrays */
  Gnum                    vsnddatsiz;
  Gnum *                  vrcvdattab;             /* Data arrays [norestrict:async]        */
  Gnum * restrict         vsnddattab;
  Gnum * restrict         procvgbtab;             /* Array of neighbor bounds [+1]         */
  int                     procngbnbr;
  int                     procngbnum;
  int * restrict          vrcvcnttab;
  int * restrict          vsndcnttab;
  int * restrict          vrcvdsptab;
  int * restrict          vsnddsptab;
  int * restrict          nsndidxtab;
  int                     vrcvdspnum;
  int                     vsnddspnum;
  Gnum                    queuheadidx;            /* Index of head of queue                */
  Gnum                    queutailidx;            /* Index of tail of queue                */
  Gnum                    bandvertlvlnum;
  Gnum                    bandvertlocnnd;
  Gnum                    bandedgelocnbr;
  Gnum                    distval;
#ifdef SCOTCH_DEBUG_DGRAPH1
  Gnum                    reduloctab[3];
  Gnum                    reduglbtab[3];
#else /* SCOTCH_DEBUG_DGRAPH1 */
  Gnum                    reduglbtab[1];
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  const Gnum * restrict const vertloctax = grafptr->vertloctax;
  const Gnum * restrict const vendloctax = grafptr->vendloctax;
  const Gnum * restrict const edgegsttax = grafptr->edgegsttax;
  const Gnum * restrict const edgeloctax = grafptr->edgeloctax;

  procngbnbr = grafptr->procngbnbr;

  reduglbtab[0] = 0;                             /* Assume everything is all right */

  vrcvdatsiz = DGRAPHBANDGROWSMUL (grafptr->procsndnbr); /* Senders and receivers inverted because we send local, not halo vertices */
  vsnddatsiz = DGRAPHBANDGROWSMUL (grafptr->vertgstnbr - grafptr->vertlocnbr);
  if (memAllocGroup ((void **) (void *)
                     &procvgbtab, (size_t) ((procngbnbr + 1)    * sizeof (Gnum)),
                     &nsndidxtab, (size_t) (procngbnbr          * sizeof (int)),
                     &vrcvcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &vsndcnttab, (size_t) (grafptr->procglbnbr * sizeof (int)), /* TRICK: vsndcnttab, vrcvdsptab, vrcvdattab, vrcvdattab joined */
                     &vrcvdsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &vsnddsptab, (size_t) (grafptr->procglbnbr * sizeof (int)),
                     &vrcvdattab, (size_t) (vrcvdatsiz          * sizeof (Gnum)),
                     &vsnddattab, (size_t) (vsnddatsiz          * sizeof (Gnum)), NULL) == NULL) {
    errorPrint (DGRAPHBANDGROWNSTR "Coll: out of memory (1)");
    reduglbtab[0] = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication not needed and not absorbable by the algorithm */
  reduloctab[0] = reduglbtab[0];
  reduloctab[1] =   distmax;
  reduloctab[2] = - distmax;
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint (DGRAPHBANDGROWNSTR "Coll: communication error (1)");
    return     (1);
  }
  if (reduglbtab[1] != - reduglbtab[2]) {
    errorPrint (DGRAPHBANDGROWNSTR "Coll: invalid parameters");
    reduglbtab[0] = 1;
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (reduglbtab[0] != 0) {
    if (vnumgsttax != NULL) {
      if (procvgbtab != NULL)
        memFree (procvgbtab);                     /* Free group leader */
      memFree (vnumgsttax);
    }
    return (1);
  }

  memSet (vsndcnttab, 0, (((int *) vrcvdattab) - vsndcnttab) * sizeof (int)); /* TRICK: vsndcnttab, vrcvdsptab, vrcvdattab, vrcvdattab joined */

  for (procngbnum = 0, vrcvdspnum = vsnddspnum = 0; /* Build communication index arrays */
       procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = grafptr->procngbtab[procngbnum];
    procvgbtab[procngbnum] = grafptr->procvrttab[procglbnum];
    vrcvdsptab[procglbnum] = vrcvdspnum;
    vsnddsptab[procglbnum] = vsnddspnum;
    vrcvdspnum += DGRAPHBANDGROWSMUL (grafptr->procsndtab[procglbnum]); /* Senders and receivers are reversed */
    vsnddspnum += DGRAPHBANDGROWSMUL (grafptr->procrcvtab[procglbnum]);
  }
  procvgbtab[procngbnum] = grafptr->procvrttab[grafptr->procglbnbr];

  bandvertlvlnum =                                /* Start index of last level is start index */
  bandvertlocnnd = grafptr->baseval;              /* Reset number of band vertices, plus base */
  bandedgelocnbr = 0;                             /* No edges accounted for yet               */

  DGRAPHBANDGROWENQ1;                             /* Record queued vertices in index array and account for their edges */

  vertlocnnd = grafptr->vertlocnnd;

  queuheadidx = 0;                                /* No queued vertex read yet                  */
  queutailidx = queulocnbr;                       /* All frontier vertices are already in queue */
  for (distval = 0; ++ distval <= distmax; ) {
    Gnum              queunextidx;                /* Tail index for enqueuing vertices of next band */
    int               procngbnum;

    bandvertlvlnum = bandvertlocnnd;              /* Save start index of current level, based */
    *bandvertlvlptr = bandvertlvlnum;

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) /* Build communication index arrays */
      nsndidxtab[procngbnum] = vsnddsptab[grafptr->procngbtab[procngbnum]];

    for (queunextidx = queutailidx; queuheadidx < queutailidx; ) { /* For all vertices in queue */
      Gnum              vertlocnum;
      Gnum              edgelocnum;

      vertlocnum = queuloctab[queuheadidx ++];    /* Dequeue vertex */
      for (edgelocnum = vertloctax[vertlocnum]; edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
        Gnum              vertlocend;

        vertlocend = edgegsttax[edgelocnum];
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        if (vertlocend < vertlocnnd) {            /* If end vertex is local        */
          vnumgsttax[vertlocend] = DGRAPHBANDGROWENQ2; /* Label vertex as enqueued */
          queuloctab[queunextidx ++] = vertlocend; /* Enqueue vertex for next pass */
          DGRAPHBANDGROWEDGE (vertlocend);        /* Account for its edges         */
        }
        else {                                    /* End vertex is a ghost */
          Gnum              vertglbend;
          int               procngbnum;
          int               procngbmax;
          int               nsndidxnum;

          vnumgsttax[vertlocend] = 0;             /* Label ghost vertex as enqueued  */
          vertglbend = edgeloctax[edgelocnum];    /* Get global number of end vertex */

          procngbnum = 0;
          procngbmax = procngbnbr;
          while ((procngbmax - procngbnum) > 1) { /* Find owner process by dichotomy on procvgbtab */
            int                 procngbmed;

            procngbmed = (procngbmax + procngbnum) / 2;
            if (procvgbtab[procngbmed] > vertglbend)
              procngbmax = procngbmed;
            else
              procngbnum = procngbmed;
          }
#ifdef SCOTCH_DEBUG_DGRAPH2
          if ((grafptr->procvrttab[grafptr->procngbtab[procngbnum]]     >  vertglbend) ||
              (grafptr->procvrttab[grafptr->procngbtab[procngbnum] + 1] <= vertglbend) ||
              (nsndidxtab[procngbnum] >= (vsnddsptab[grafptr->procngbtab[procngbnum]] + grafptr->procrcvtab[grafptr->procngbtab[procngbnum]]))) {
            errorPrint (DGRAPHBANDGROWNSTR "Coll: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          nsndidxnum = nsndidxtab[procngbnum];
          vsnddattab[nsndidxnum ++] = vertglbend - procvgbtab[procngbnum] + grafptr->baseval; /* Buffer local value on neighbor processor */
          DGRAPHBANDGROWENQ3;                     /* Send additional data if needed */
          nsndidxtab[procngbnum] = nsndidxnum;
        }
      }
    }

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
      int                 procglbnum;

      procglbnum = grafptr->procngbtab[procngbnum];
      vsndcnttab[procglbnum] = nsndidxtab[procngbnum] - vsnddsptab[procglbnum];
    }

    if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, vrcvcnttab, 1, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint (DGRAPHBANDGROWNSTR "Coll: communication error (2)");
      return     (1);
    }
    if (MPI_Alltoallv (vsnddattab, vsndcnttab, vsnddsptab, GNUM_MPI,
                       vrcvdattab, vrcvcnttab, vrcvdsptab, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint (DGRAPHBANDGROWNSTR "Coll: communication error (3)");
      return     (1);
    }

    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* For all receive buffers */
      Gnum * restrict   vrcvdatptr;
      int               vertrcvnum;
      int               procglbnum;
      int               statsiz;

      procglbnum = grafptr->procngbtab[procngbnum];
      vrcvdatptr = vrcvdattab + vrcvdsptab[procglbnum];
      statsiz    = vrcvcnttab[procglbnum];
      for (vertrcvnum = 0; vertrcvnum < statsiz; vertrcvnum += DGRAPHBANDGROWSMUL (1)) {
        Gnum              vertlocend;

        vertlocend = vrcvdatptr[vertrcvnum];      /* Get local vertex from message */
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocend < grafptr->baseval) || (vertlocend >= vertlocnnd)) {
          errorPrint (DGRAPHBANDGROWNSTR "Coll: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        vnumgsttax[vertlocend] = DGRAPHBANDGROWENQ4; /* Label vertex as enqueued  */
        queuloctab[queunextidx ++] = vertlocend;  /* Enqueue vertex for next pass */
        DGRAPHBANDGROWEDGE (vertlocend);          /* Account for its edges        */
      }
    }

    queutailidx = queunextidx;                    /* Prepare queue for next sweep */
  }

  memFree (procvgbtab);                           /* Free group leader */

  *bandvertlocptr = bandvertlocnnd - grafptr->baseval;
  *bandedgelocptr = bandedgelocnbr;

  return (0);
}

/* This routine computes a distributed index array
** of given width around the current separator,
** using point-to-point communication.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
DGRAPHBANDGROWNAMEPTOP (
Dgraph * restrict const           grafptr,        /*+ Distributed graph                                        +*/
const Gnum                        queulocnbr,     /*+ Number of frontier vertices, start size for vertex queue +*/
Gnum * restrict const             queuloctab,     /*+ Array of frontier vertices, re-used as queue array       +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices                 +*/
Gnum * restrict const             vnumgsttax,     /*+ Flag or index array to fill                              +*/
Gnum * restrict const             bandvertlvlptr, /*+ Pointer to based start index of last level               +*/
Gnum * restrict const             bandvertlocptr, /*+ Pointer to bandvertlocnnd                                +*/
Gnum * restrict const             bandedgelocptr) /*+ Pointer to bandedgelocnbr                                +*/
{
  Gnum                    queulocnum;
  Gnum                    vertlocnnd;
  Gnum                    vnumgstsiz;             /* Size of vnumgsttax; TRICK: re-use     */
  Gnum                    vrcvdatsiz;             /* Sizes of data send and receive arrays */
  Gnum                    vsnddatsiz;
  Gnum *                  vrcvdattab;             /* Data arrays [norestrict:async]        */
  Gnum *                  vsnddattab;
  Gnum * restrict         procvgbtab;             /* Array of neighbor bounds [+1]         */
  int                     procngbnbr;
  int                     procngbnum;
  int                     procngbnxt;
  Gnum * restrict         nrcvdsptab;
  Gnum * restrict         nsnddsptab;
  Gnum                    nrcvdspnum;
  Gnum                    nsnddspnum;
  Gnum * restrict         nsndidxtab;
  Gnum                    queuheadidx;            /* Index of head of queue                */
  Gnum                    queutailidx;            /* Index of tail of queue                */
  MPI_Request *           nrcvreqtab;             /* Array of receive requests             */
  MPI_Request *           nsndreqtab;             /* Array of receive requests             */
  Gnum                    bandvertlvlnum;
  Gnum                    bandvertlocnnd;
  Gnum                    bandedgelocnbr;
  Gnum                    distval;
#ifdef SCOTCH_DEBUG_DGRAPH1
  Gnum                    reduloctab[3];
  Gnum                    reduglbtab[3];
#else /* SCOTCH_DEBUG_DGRAPH1 */
  Gnum                    reduglbtab[1];
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  const Gnum * restrict const vertloctax = grafptr->vertloctax;
  const Gnum * restrict const vendloctax = grafptr->vendloctax;
  const Gnum * restrict const edgegsttax = grafptr->edgegsttax;
  const Gnum * restrict const edgeloctax = grafptr->edgeloctax;

  procngbnbr = grafptr->procngbnbr;

  reduglbtab[0] = 0;                             /* Assume everything is all right */

  vrcvdatsiz = DGRAPHBANDGROWSMUL (grafptr->procsndnbr); /* Senders and receivers inverted because we send local, not halo vertices */
  vsnddatsiz = DGRAPHBANDGROWSMUL (grafptr->vertgstnbr - grafptr->vertlocnbr);
  if (memAllocGroup ((void **) (void *)
                     &procvgbtab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)),
                     &nrcvdsptab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)), /* +1 to check against end of array */
                     &nsnddsptab, (size_t) ((procngbnbr + 1) * sizeof (Gnum)), /* +1 to check against end of array */
                     &nsndidxtab, (size_t) (procngbnbr       * sizeof (Gnum)), /* Here Gnum's since point-to-point */
                     &nrcvreqtab, (size_t) (procngbnbr       * sizeof (MPI_Request)),
                     &nsndreqtab, (size_t) (procngbnbr       * sizeof (MPI_Request)),
                     &vrcvdattab, (size_t) (vrcvdatsiz       * sizeof (Gnum)),
                     &vsnddattab, (size_t) (vsnddatsiz       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint (DGRAPHBANDGROWNSTR "Ptop: out of memory (1)");
    reduglbtab[0] = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication not needed and not absorbable by the algorithm */
  reduloctab[0] = reduglbtab[0];
  reduloctab[1] =   distmax;
  reduloctab[2] = - distmax;
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, GNUM_MPI, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (1)");
    return     (1);
  }
  if (reduglbtab[1] != - reduglbtab[2]) {
    errorPrint (DGRAPHBANDGROWNSTR "Ptop: invalid parameters");
    reduglbtab[0] = 1;
  }
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (reduglbtab[0] != 0) {
    if (vnumgsttax != NULL) {
      if (procvgbtab != NULL)
        memFree (procvgbtab);                     /* Free group leader */
      memFree (vnumgsttax);
    }
    return (1);
  }

  for (procngbnum = 0, nrcvdspnum = nsnddspnum = procngbnxt = 0; /* Build communication index arrays */
       procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = grafptr->procngbtab[procngbnum];
    if ((procngbnxt == 0) && (procglbnum > grafptr->proclocnum)) /* Find index of first neighbor of higher rank */
      procngbnxt = procngbnum;
    procvgbtab[procngbnum] = grafptr->procvrttab[procglbnum];
    nrcvdsptab[procngbnum] = nrcvdspnum;          /* Arrays are indexed per neighbor since we are doing point-to-point communication */
    nsnddsptab[procngbnum] = nsnddspnum;
    nrcvdspnum += DGRAPHBANDGROWSMUL (grafptr->procsndtab[procglbnum]); /* Senders and receivers are reversed */
    nsnddspnum += DGRAPHBANDGROWSMUL (grafptr->procrcvtab[procglbnum]);
  }
  procvgbtab[procngbnum] = grafptr->procvrttab[grafptr->procglbnbr];
  nrcvdsptab[procngbnum] = nrcvdspnum;            /* Mark end of communication index arrays */
  nsnddsptab[procngbnum] = nsnddspnum;

  procngbnum = procngbnxt;                        /* Create receive requests in descending order */
  if (procngbnbr != 0) {
    do {
      procngbnum = (procngbnum + (procngbnbr - 1)) % procngbnbr; /* Pre-decrement neighbor rank */
      if (MPI_Recv_init (vrcvdattab + nrcvdsptab[procngbnum], (int) (nrcvdsptab[procngbnum + 1] - nrcvdsptab[procngbnum]), GNUM_MPI,
                         grafptr->procngbtab[procngbnum], TAGBAND,
                         grafptr->proccomm, nrcvreqtab + procngbnum) != MPI_SUCCESS) {
        errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (2)");
        return     (1);
      }
    } while (procngbnum != procngbnxt);
  }

  bandvertlvlnum =                                /* Start index of last level is start index */
  bandvertlocnnd = grafptr->baseval;              /* Reset number of band vertices, plus base */
  bandedgelocnbr = 0;                             /* No edges accounted for yet               */

  DGRAPHBANDGROWENQ1;                             /* Record queued vertices in index array and account for their edges */

  vertlocnnd = grafptr->vertlocnnd;

  queuheadidx = 0;                                /* No queued vertex read yet                  */
  queutailidx = queulocnbr;                       /* All frontier vertices are already in queue */
  for (distval = 0; ++ distval <= distmax; ) {
    Gnum              queunextidx;                /* Tail index for enqueuing vertices of next band */
    int               vrcvreqnbr;

    if (MPI_Startall (procngbnbr, nrcvreqtab) != MPI_SUCCESS) { /* Start all receive operations from neighbors */
      errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (3)");
      return     (1);
    }

    bandvertlvlnum = bandvertlocnnd;              /* Save start index of current level, based */
    *bandvertlvlptr = bandvertlvlnum;

    memCpy (nsndidxtab, nsnddsptab, procngbnbr * sizeof (Gnum)); /* Reset send buffer indices */

    for (queunextidx = queutailidx; queuheadidx < queutailidx; ) { /* For all vertices in queue */
      Gnum              vertlocnum;
      Gnum              edgelocnum;

      vertlocnum = queuloctab[queuheadidx ++];    /* Dequeue vertex */
      for (edgelocnum = vertloctax[vertlocnum]; edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
        Gnum              vertlocend;

        vertlocend = edgegsttax[edgelocnum];
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        if (vertlocend < vertlocnnd) {            /* If end vertex is local        */
          vnumgsttax[vertlocend] = DGRAPHBANDGROWENQ2; /* Label vertex as enqueued */
          queuloctab[queunextidx ++] = vertlocend; /* Enqueue vertex for next pass */
          DGRAPHBANDGROWEDGE (vertlocend);        /* Account for its edges         */
        }
        else {                                    /* End vertex is a ghost */
          Gnum              vertglbend;
          int               procngbnum;
          int               procngbmax;
          int               nsndidxnum;

          vnumgsttax[vertlocend] = 0;             /* Label ghost vertex as enqueued  */
          vertglbend = edgeloctax[edgelocnum];    /* Get global number of end vertex */

          procngbnum = 0;
          procngbmax = procngbnbr;
          while ((procngbmax - procngbnum) > 1) { /* Find owner process by dichotomy on procvgbtab */
            int                 procngbmed;

            procngbmed = (procngbmax + procngbnum) / 2;
            if (procvgbtab[procngbmed] > vertglbend)
              procngbmax = procngbmed;
            else
              procngbnum = procngbmed;
          }
#ifdef SCOTCH_DEBUG_DGRAPH2
          if ((grafptr->procvrttab[grafptr->procngbtab[procngbnum]]     >  vertglbend) ||
              (grafptr->procvrttab[grafptr->procngbtab[procngbnum] + 1] <= vertglbend) ||
              (nsndidxtab[procngbnum] >= nsnddsptab[procngbnum + 1])) {
            errorPrint (DGRAPHBANDGROWNSTR "Ptop: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          nsndidxnum = nsndidxtab[procngbnum];
          vsnddattab[nsndidxnum ++] = vertglbend - procvgbtab[procngbnum] + grafptr->baseval; /* Buffer local value on neighbor processor */
          DGRAPHBANDGROWENQ3;                     /* Send additional data if needed */
          nsndidxtab[procngbnum] = nsndidxnum;
        }
      }
    }

    procngbnum = procngbnxt;                      /* Send all buffers to neighbors */
    if (procngbnbr != 0) {
      do {
        int               procglbnum;

        procglbnum = grafptr->procngbtab[procngbnum];

        if (MPI_Isend (vsnddattab + nsnddsptab[procngbnum], nsndidxtab[procngbnum] - nsnddsptab[procngbnum],
                       GNUM_MPI, grafptr->procngbtab[procngbnum], TAGBAND, grafptr->proccomm,
                       nsndreqtab + procngbnum) != MPI_SUCCESS) {
          errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (4)");
          return     (1);
        }
        procngbnum = (procngbnum + 1) % procngbnbr; /* Post-increment neighbor rank */
      } while (procngbnum != procngbnxt);
    }

    for (vrcvreqnbr = procngbnbr; vrcvreqnbr > 0; vrcvreqnbr --) { /* For all pending receive requests */
      Gnum * restrict   vrcvdatptr;
      int               vertrcvnum;
      MPI_Status        statdat;
      int               statsiz;
      int               o;

#ifdef SCOTCH_DETERMINISTIC
      procngbnum = vrcvreqnbr - 1;
      o = MPI_Wait (&nrcvreqtab[procngbnum], &statdat);
#else /* SCOTCH_DETERMINISTIC */
      o = MPI_Waitany (procngbnbr, nrcvreqtab, &procngbnum, &statdat);
#endif /* SCOTCH_DETERMINISTIC */
      if ((o != MPI_SUCCESS) ||
          (MPI_Get_count (&statdat, GNUM_MPI, &statsiz) != MPI_SUCCESS)) {
        errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (5)");
        return     (1);
      }
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (statdat.MPI_SOURCE != grafptr->procngbtab[procngbnum]) {
        errorPrint (DGRAPHBANDGROWNSTR "Ptop: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      vrcvdatptr = vrcvdattab + nrcvdsptab[procngbnum];
      for (vertrcvnum = 0; vertrcvnum < statsiz; vertrcvnum += DGRAPHBANDGROWSMUL (1)) {
        Gnum              vertlocend;

        vertlocend = vrcvdatptr[vertrcvnum];      /* Get local vertex from message */
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocend < grafptr->baseval) || (vertlocend >= vertlocnnd)) {
          errorPrint (DGRAPHBANDGROWNSTR "dgraphBandPtop: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        if (vnumgsttax[vertlocend] != ~0)         /* If end vertex has already been processed */
          continue;                               /* Skip to next vertex                      */

        vnumgsttax[vertlocend] = DGRAPHBANDGROWENQ4; /* Label vertex as enqueued  */
        queuloctab[queunextidx ++] = vertlocend;  /* Enqueue vertex for next pass */
        DGRAPHBANDGROWEDGE (vertlocend);          /* Account for its edges        */
      }
    }

    queutailidx = queunextidx;                    /* Prepare queue for next sweep */

    if (MPI_Waitall (procngbnbr, nsndreqtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait until all send operations completed */
      errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (6)");
      return     (1);
    }
  }
  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* Free persistent receive requests */
    if (MPI_Request_free (nrcvreqtab + procngbnum) != MPI_SUCCESS) {
      errorPrint (DGRAPHBANDGROWNSTR "Ptop: communication error (7)");
      return     (1);
    }
  }

  memFree (procvgbtab);                           /* Free group leader */

  *bandvertlocptr = bandvertlocnnd - grafptr->baseval;
  *bandedgelocptr = bandedgelocnbr;

  return (0);
}
