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
/**   NAME       : dgraph_band.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a distributed band **/
/**                graph from the given frontier array.    **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 20 feb 2011     **/
/**                # Version 6.0  : from : 03 apr 2012     **/
/**                                 to   : 26 sep 2012     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  vdgraph_separate_bd.c in version 5.0. **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPH_BAND

#include "module.h"
#include "common.h"
#include "dgraph.h"

/**********************************/
/*                                */
/* Distance computation routines. */
/*                                */
/**********************************/

#define DGRAPHBANDGROWNAME          dgraphBand
#define DGRAPHBANDGROWEDGE(n)       bandedgelocnbr += vendloctax[n] - vertloctax[n]
#define DGRAPHBANDGROWENQ1          for (queulocnum = 0; queulocnum < queulocnbr; queulocnum ++) { \
                                      Gnum              vertlocnum;                                \
                                                                                                   \
                                      vertlocnum = queuloctab[queulocnum];                         \
                                      vnumgsttax[vertlocnum] = bandvertlocnnd ++;                  \
                                      DGRAPHBANDGROWEDGE (vertlocnum);                             \
                                    }
#define DGRAPHBANDGROWENQ2          bandvertlocnnd ++
#define DGRAPHBANDGROWENQ3                        /* Nothing more to send */
#define DGRAPHBANDGROWENQ4          bandvertlocnnd ++
#define DGRAPHBANDGROWSMUL(n)       (n)
#include "dgraph_band_grow.h"
#include "dgraph_band_grow.c"
#undef DGRAPHBANDGROWNAME
#undef DGRAPHBANDGROWEDGE
#undef DGRAPHBANDGROWENQU
#undef DGRAPHBANDGROWSMUL

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed index array
** of given width around the current separator.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
dgraphBand (
Dgraph * restrict const             grafptr,      /*+ Distributed graph                                      +*/
const Gnum                          fronlocnbr,   /*+ Number of frontier vertices                            +*/
Gnum * restrict const               fronloctab,   /*+ Array of frontier vertices, re-used as queue array     +*/
const GraphPart * restrict const    partgsttax,   /*+ Part array for original graph ({0,1} or {0,1,2})       +*/
const Gnum                          complocload0, /*+ Load in part 0 or {0,2}                                +*/
const Gnum                          complocload1, /*+ Load in part 1                                         +*/
Gnum                                distmax,      /*+ Maximum distance from separator vertices               +*/
Dgraph * restrict const             bandgrafptr,  /*+ Pointer to band graph structure to fill                +*/
Gnum * restrict * const             bandfronlocptr, /*+ Pointer to bandfronloctab                            +*/
GraphPart * restrict * const        bandpartgstptr, /*+ Pointer to bandpartgsttax                            +*/
Gnum * const                        bandvertlvlptr, /*+ Pointer to based start index of last level           +*/
Gnum * const                        bandvertlocptr1, /*+ Pointer to number of band vertices in part 1        +*/
Gnum * const                        bandvertlocancptr) /*+ Pointer to flag set if anchor vertices overloaded +*/
{
  Gnum                    bandvertlocnnd;         /* End of local band vertex array, (without anchor vertices) */
  Gnum                    bandvertlocnbr;         /* Number of local band vertices (including anchor vertices) */
  Gnum                    bandvertlocnbr1;        /* Number of band graph vertices in part 1 except anchor 1   */
  Gnum                    bandvertlvlnum;         /* Index of first band vertex belonging to last level        */
  Gnum * restrict         bandvertloctax;
  Gnum                    bandvertlocadj;         /* Ajust value for local-to-global band vertex indices       */
  Gnum                    bandvertlocancadj;      /* Flag set when anchor(s) represent unexistent vertices     */
  Gnum                    bandvertlocnum;
  Gnum                    bandvelolocsum;
  Gnum                    bandvelolocsum1;
  Gnum * restrict         bandedgeloctax;
  Gnum                    bandedgelocnum;
  Gnum                    bandedgeloctmp;
  Gnum                    bandedgelocnbr;         /* Number of local edges in band graph                       */
  Gnum * restrict         bandedloloctax;
  Gnum                    bandedlolocnbr;         /* Size of local band edge load array                        */
  Gnum * restrict         bandfronloctab;
  GraphPart * restrict    bandpartgsttax;
  Gnum                    bandvnumgstsiz;
  Gnum * restrict         bandvnumgsttax;         /* Indices of selected band vertices in band graph           */
  Gnum                    banddegrlocmax;
  Gnum                    degrval;
  Gnum                    veloval;
  const Gnum * restrict   edgegsttax;
  Gnum                    fronlocnum;
  int                     cheklocval;
  int                     chekglbval;
  int                     procngbnum;

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint ("dgraphBand: cannot compute ghost edge array");
    return     (1);
  }

  cheklocval = 0;
  bandvnumgstsiz = MAX ((grafptr->vertgstnbr * sizeof (Gnum)), (grafptr->procglbnbr * sizeof (int))); /* TRICK: re-use array for further error collective communications */
  if ((bandvnumgsttax = memAlloc (bandvnumgstsiz)) == NULL) {
    errorPrint ("dgraphBand: out of memory (1)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* This communication cannot be covered by a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphBand: communication error (1)");
    return     (1);
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0)
    return (1);

  memSet (bandvnumgsttax, ~0, grafptr->vertgstnbr * sizeof (Gnum)); /* Reset part array */
  bandvnumgsttax -= grafptr->baseval;

  if ((((grafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphBandPtop : dgraphBandColl)
      (grafptr, fronlocnbr, fronloctab, distmax, bandvnumgsttax, &bandvertlvlnum, &bandvertlocnbr, &bandedgelocnbr) != 0)
    return (1);

  if (bandvertlvlptr != NULL)
    *bandvertlvlptr = bandvertlvlnum;

  bandedgelocnbr += 2 * ((bandvertlocnbr + grafptr->baseval - bandvertlvlnum) + (grafptr->procglbnbr - 1)); /* Add edges to and from anchors */
  bandvertlocnbr += 2;                            /* Add anchor vertices */
  bandedlolocnbr  = (grafptr->edloloctax != NULL) ? bandedgelocnbr : 0;

  dgraphInit (bandgrafptr, grafptr->proccomm);
  bandgrafptr->flagval = (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* Arrays created by the routine itself */
  bandgrafptr->baseval = grafptr->baseval;

  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &bandgrafptr->procdsptab, (size_t) ((grafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &bandgrafptr->proccnttab, (size_t) (grafptr->procglbnbr       * sizeof (Gnum)),
                     &bandgrafptr->procngbtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procrcvtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procsndtab, (size_t) (grafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphBand: out of memory (2)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &bandgrafptr->vertloctax, (size_t) ((bandvertlocnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &bandgrafptr->vnumloctax, (size_t) (bandvertlocnbr       * sizeof (Gnum)),
                          &bandgrafptr->veloloctax, (size_t) (bandvertlocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphBand: out of memory (3)");
    cheklocval = 1;
  }
  else if (bandgrafptr->vertloctax -= bandgrafptr->baseval,
           bandgrafptr->veloloctax -= bandgrafptr->baseval,
           bandgrafptr->vnumloctax -= bandgrafptr->baseval,
           (memAllocGroup ((void **) (void *)
                           &bandedgeloctax, (size_t) (bandedgelocnbr * sizeof (Gnum)),
                           &bandedloloctax, (size_t) (bandedlolocnbr * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("dgraphBand: out of memory (4)");
    cheklocval = 1;
  }
  else {
    bandedgeloctax -= bandgrafptr->baseval;
    bandedloloctax  = (grafptr->edloloctax != NULL) ? (bandedloloctax - bandgrafptr->baseval) : NULL;

    if ((bandfronloctab = memAlloc (bandvertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("dgraphBand: out of memory (5)");
      cheklocval = 1;
    }
    else if ((bandpartgsttax = memAlloc ((bandvertlocnbr + bandedgelocnbr) * sizeof (GraphPart))) == NULL) { /* Upper bound on number of ghost vertices */
      errorPrint ("dgraphBand: out of memory (6)");
      cheklocval = 1;
    }
    else
      bandpartgsttax -= bandgrafptr->baseval;
  }

  if (cheklocval != 0) {                          /* In case of memory error */
    bandgrafptr->procdsptab[0] = -1;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI, /* Send received data to dummy array */
                       bandvnumgsttax + bandgrafptr->baseval, 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBand: communication error (2)");
      return     (1);
    }
    if (bandfronloctab != NULL)
      memFree (bandfronloctab);
    dgraphFree (bandgrafptr);
    memFree    (bandvnumgsttax + bandgrafptr->baseval);
    return     (1);
  }
  else {
    bandgrafptr->procdsptab[0] = bandvertlocnbr;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &bandgrafptr->procdsptab[1], 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphBand: communication error (3)");
      return     (1);
    }
  }
  bandgrafptr->procdsptab[0] = bandgrafptr->baseval; /* Build vertex-to-process array */
#ifdef SCOTCH_DEBUG_DGRAPH2
  memSet (bandgrafptr->vnumloctax + bandgrafptr->baseval, ~0, (bandvertlocnbr * sizeof (Gnum)));
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (procngbnum = 1; procngbnum <= grafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
    if (bandgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                       */
      if (bandpartgsttax != NULL)
        memFree (bandpartgsttax + bandgrafptr->baseval);
      if (bandfronloctab != NULL)
        memFree (bandfronloctab);
      dgraphFree (bandgrafptr);
      memFree    (bandvnumgsttax + bandgrafptr->baseval);
      return     (1);
    }
    bandgrafptr->procdsptab[procngbnum]    += bandgrafptr->procdsptab[procngbnum - 1];
    bandgrafptr->proccnttab[procngbnum - 1] = bandgrafptr->procdsptab[procngbnum] - bandgrafptr->procdsptab[procngbnum - 1];
  }

  for (fronlocnum = 0, bandvertlocnum = bandgrafptr->baseval, bandvertlocadj = bandgrafptr->procdsptab[grafptr->proclocnum] - bandgrafptr->baseval;
       fronlocnum < fronlocnbr; fronlocnum ++, bandvertlocnum ++) { /* Turn all graph frontier vertices into band frontier vertices */
    Gnum              vertlocnum;

    bandfronloctab[fronlocnum] = bandvertlocnum;  /* All frontier vertices are first vertices of band graph */
    vertlocnum = fronloctab[fronlocnum];
    bandgrafptr->vnumloctax[bandvertlocnum] = vertlocnum;
    bandvnumgsttax[vertlocnum] += bandvertlocadj; /* Turn local indices in band graph into global indices */
  }
  for (bandvertlocnnd = bandvertlocnbr + bandgrafptr->baseval - 2; /* Pick selected band vertices from rest of frontier array without anchors */
       bandvertlocnum < bandvertlocnnd; fronlocnum ++, bandvertlocnum ++) {
    Gnum              vertlocnum;

    vertlocnum = fronloctab[fronlocnum];
    bandgrafptr->vnumloctax[bandvertlocnum] = vertlocnum;
    bandvnumgsttax[vertlocnum] += bandvertlocadj; /* Turn local indices in band graph into global indices */
  }
  bandgrafptr->vnumloctax[bandvertlocnnd]     =   /* Prevent Valgrind from yelling when centralizing band graphs */
  bandgrafptr->vnumloctax[bandvertlocnnd + 1] = -1;

  if (dgraphHaloSync (grafptr, (byte *) (bandvnumgsttax + bandgrafptr->baseval), GNUM_MPI) != 0) { /* Share global indexing of halo vertices */
    errorPrint ("dgraphBand: cannot perform halo exchange");
    return     (1);
  }

  edgegsttax = grafptr->edgegsttax;

  veloval = 1;
  bandvertloctax  = bandgrafptr->vertloctax;
  bandvertlocnbr1 = 0;
  bandvelolocsum  = 0;
  bandvelolocsum1 = 0;
  banddegrlocmax  = 0;
  for (bandvertlocnum = bandedgelocnum = bandgrafptr->baseval; /* Build global vertex array of band graph             */
       bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { /* For all vertices that do not belong to the last level */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;
    GraphPart         partval;
    Gnum              partval1;

    vertlocnum = bandgrafptr->vnumloctax[bandvertlocnum];
    partval    = partgsttax[vertlocnum];
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (partval > 2) {
      errorPrint ("dgraphBand: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    partval1 = partval & 1;
    bandvertlocnbr1 += partval1;                  /* Count vertices in part 1 */
    bandpartgsttax[bandvertlocnum] = partval;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (grafptr->veloloctax != NULL) {
      veloval = grafptr->veloloctax[vertlocnum];
      bandvelolocsum  += veloval;
      bandvelolocsum1 += veloval & (- partval1);  /* Sum vertex load if (partval == 1) */
    }
    bandgrafptr->veloloctax[bandvertlocnum] = veloval;

    degrval = grafptr->vendloctax[vertlocnum] - grafptr->vertloctax[vertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;

    for (edgelocnum = grafptr->vertloctax[vertlocnum]; /* For all original edges */
         edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (bandvnumgsttax[edgegsttax[edgelocnum]] == ~0) { /* All ends should belong to the band graph too */
        errorPrint ("dgraphBand: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      bandedgeloctax[bandedgelocnum ++] = bandvnumgsttax[edgegsttax[edgelocnum]];
    }
  }
  for ( ; bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) { /* For all vertices that belong to the last level except anchors */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;
    GraphPart         partval;
    Gnum              partval1;

    vertlocnum = bandgrafptr->vnumloctax[bandvertlocnum];
    partval    = partgsttax[vertlocnum];
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (partval > 2) {
      errorPrint ("dgraphBand: internal error (3)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    partval1 = partval & 1;
    bandvertlocnbr1 += partval1;                  /* Count vertices in part 1 */
    bandpartgsttax[bandvertlocnum] = partval;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (grafptr->veloloctax != NULL) {
      veloval = grafptr->veloloctax[vertlocnum];
      bandvelolocsum  += veloval;
      bandvelolocsum1 += veloval & (- partval1);  /* Sum vertex load if (partval == 1) */
    }
    bandgrafptr->veloloctax[bandvertlocnum] = veloval;

    for (edgelocnum = grafptr->vertloctax[vertlocnum]; /* For all original edges */
         edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
      Gnum              bandvertlocend;

      bandvertlocend = bandvnumgsttax[edgegsttax[edgelocnum]];
      if (bandvertlocend != ~0) {                 /* If end vertex belongs to band graph  */
        if (bandedloloctax != NULL)               /* If graph has edge weights, copy load */
          bandedloloctax[bandedgelocnum] = grafptr->edloloctax[edgelocnum];
        bandedgeloctax[bandedgelocnum ++] = bandvertlocend;
      }
    }
    if (bandedloloctax != NULL)                   /* If graph has edge weights */
      bandedloloctax[bandedgelocnum] = 1;         /* Edge to anchor has load 1 */
    bandedgeloctax[bandedgelocnum ++] = bandvertlocnnd + bandvertlocadj + partval1; /* Add edge to anchor of proper part */

    degrval = bandedgelocnum - bandvertloctax[bandvertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;
  }

  memFree (bandvnumgsttax + bandgrafptr->baseval); /* Free useless space */

  bandpartgsttax[bandvertlocnnd]     = 0;         /* Set parts of anchor vertices */
  bandpartgsttax[bandvertlocnnd + 1] = 1;
  bandvertloctax[bandvertlocnum] = bandedgelocnum; /* Process anchor vertex in part 0                          */
  for (procngbnum = 0; procngbnum < grafptr->proclocnum; procngbnum ++) /* Build clique with anchors of part 0 */
    bandedgeloctax[bandedgelocnum ++] = bandgrafptr->procdsptab[procngbnum + 1] - 2;
  for (procngbnum ++; procngbnum < grafptr->procglbnbr; procngbnum ++) /* Build clique with anchors of part 0 */
    bandedgeloctax[bandedgelocnum ++] = bandgrafptr->procdsptab[procngbnum + 1] - 2;
  bandedgeloctmp = bandedgelocnum + (bandvertlocnnd - bandvertlvlnum);
  for (procngbnum = 0; procngbnum < grafptr->proclocnum; procngbnum ++) /* Build clique with anchors of part 1 */
    bandedgeloctax[bandedgeloctmp ++] = bandgrafptr->procdsptab[procngbnum + 1] - 1;
  for (procngbnum ++; procngbnum < grafptr->procglbnbr; procngbnum ++) /* Build clique with anchors of part 1 */
    bandedgeloctax[bandedgeloctmp ++] = bandgrafptr->procdsptab[procngbnum + 1] - 1;
  bandvertloctax[bandvertlocnnd + 2] = bandedgeloctmp;
  bandedgelocnbr = bandedgeloctmp - bandgrafptr->baseval; /* Set real number of edges */
  for (bandvertlocnum = bandvertlvlnum, bandedgeloctmp = bandedgelocnum + (bandvertlocnnd - bandvertlvlnum); /* Link vertices of last level to anchors */
       bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) {
    if (bandpartgsttax[bandvertlocnum] == 0)
      bandedgeloctax[bandedgelocnum ++] = bandvertlocnum + bandvertlocadj;
    else
      bandedgeloctax[-- bandedgeloctmp] = bandvertlocnum + bandvertlocadj;
  }
  bandvertloctax[bandvertlocnnd + 1] = bandedgeloctmp;
  degrval = bandvertloctax[bandvertlocnnd + 1] - bandvertloctax[bandvertlocnnd];
  if (banddegrlocmax < degrval)
    banddegrlocmax = degrval;
  degrval = bandvertloctax[bandvertlocnnd + 2] - bandvertloctax[bandvertlocnnd + 1];
  if (banddegrlocmax < degrval)
    banddegrlocmax = degrval;

  if (bandedloloctax != NULL) {                   /* If graph has edge weights */
    Gnum              edgelocnum;
    Gnum              edgelocnnd;

    for (bandvertlocnum = bandgrafptr->baseval;   /* For all vertices that do not belong to the last level */
         bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { 
      Gnum              vertlocnum;
      Gnum              bandedgelocnum;

      vertlocnum     = bandgrafptr->vnumloctax[bandvertlocnum];
      bandedgelocnum = bandvertloctax[bandvertlocnum];
      memCpy (bandedloloctax + bandedgelocnum,    /* Copy edge load array */
              &grafptr->edloloctax[grafptr->vertloctax[vertlocnum]],
              (bandvertloctax[bandvertlocnum + 1] - bandedgelocnum) * sizeof (Gnum));
    }                                             /* Vertices of last level have been processed before */
    for (edgelocnum = bandvertloctax[bandvertlocnnd], /* Loads of anchor edges are all 1's too         */
         edgelocnnd = bandvertloctax[bandvertlocnnd + 2];
         edgelocnum < edgelocnnd; edgelocnum ++)
      bandedloloctax[edgelocnum] = 1;
  }

  if (grafptr->veloloctax == NULL) {              /* If original graph is not weighted */
    bandgrafptr->veloloctax[bandvertlocnnd]     = complocload0 + bandvertlocnbr1 - bandvertlocnbr + 2; /* Plus 2 for anchors */
    bandgrafptr->veloloctax[bandvertlocnnd + 1] = complocload1 - bandvertlocnbr1;
  }
  else {
    bandgrafptr->veloloctax[bandvertlocnnd]     = complocload0 + bandvelolocsum1 - bandvelolocsum;
    bandgrafptr->veloloctax[bandvertlocnnd + 1] = complocload1 - bandvelolocsum1;
  }
  bandvertlocancadj = 0;
  if ((bandgrafptr->veloloctax[bandvertlocnnd]     == 0) || /* If at least one anchor is empty */
      (bandgrafptr->veloloctax[bandvertlocnnd + 1] == 0)) {
    bandvertlocancadj = 1;
    bandgrafptr->veloloctax[bandvertlocnnd] ++;   /* Increase weight of both anchors to keep balance */
    bandgrafptr->veloloctax[bandvertlocnnd + 1] ++;
  }

  bandgrafptr->procvrttab = bandgrafptr->procdsptab; /* Graph does not have holes */
  bandgrafptr->vertlocnbr = bandvertlocnbr;
  bandgrafptr->vertlocnnd = bandvertlocnbr + bandgrafptr->baseval;
  bandgrafptr->vendloctax = bandvertloctax + 1;   /* Band graph is compact */
  bandgrafptr->velolocsum = grafptr->velolocsum + 2 * bandvertlocancadj;
  bandgrafptr->edgeloctax = bandedgeloctax;
  bandgrafptr->edloloctax = bandedloloctax;
  bandgrafptr->edgelocnbr = bandedgelocnbr;
  bandgrafptr->edgelocsiz = bandedgelocnbr;
  bandgrafptr->degrglbmax = banddegrlocmax;       /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (bandgrafptr) != 0) {
    errorPrint ("dgraphBand: cannot build band graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (bandgrafptr) != 0) {
    errorPrint ("dgraphBand: internal error (4)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  *bandfronlocptr    = bandfronloctab;
  *bandpartgstptr    = bandpartgsttax;
  *bandvertlocptr1   = bandvertlocnbr1;
  *bandvertlocancptr = bandvertlocancadj;

  return (0);
}
