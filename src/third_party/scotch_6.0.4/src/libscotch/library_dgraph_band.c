/* Copyright 2011,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_band.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted band graph building routine of    **/
/**                the libScotch library.                  **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 28 oct 2011     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_halo.h"
#include "ptscotch.h"

#define DGRAPHBANDGROWNAME          dgraphBand
#include "dgraph_band_grow.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mapping routines.            */
/*                                  */
/************************************/

/*+ This routine builds a distributed
*** band graph, without anchors, from the
*** given distributed graph.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphBand (
SCOTCH_Dgraph * const       orggrafptr,
const SCOTCH_Num            fronlocnbr,
SCOTCH_Num * const          fronloctab,
const SCOTCH_Num            distval,
SCOTCH_Dgraph * const       bndgrafptr)
{
  MPI_Comm              bandproccomm;
  Gnum                  bandvertlocnnd;           /* End of local band vertex array                      */
  Gnum                  bandvertlocnbr;           /* Number of local band vertices                       */
  Gnum                  bandvertlvlnum;           /* Index of first band vertex belonging to last level  */
  Gnum * restrict       bandvertloctax;
  Gnum                  bandvertlocadj;           /* Ajust value for local-to-global band vertex indices */
  Gnum                  bandvertlocnum;
  Gnum * restrict       bandveloloctax;
  Gnum                  bandvelolocnbr;
  Gnum                  bandvelolocsum;
  Gnum * restrict       bandedgeloctax;
  Gnum                  bandedgelocnum;
  Gnum                  bandedgelocsiz;           /* Number of local edges in band graph                 */
  Gnum * restrict       bandedloloctax;
  Gnum                  bandedlolocsiz;           /* Size of local band edge load array                  */
  Gnum                  bandvnumgstsiz;
  Gnum * restrict       bandvnumgsttax;           /* Indices of selected band vertices in band graph     */
  Gnum * restrict       bandvlblloctax;
  Gnum                  banddegrlocmax;
  Gnum                  degrval;
  Gnum                  veloval;
  Gnum                  vertlocadj;
  const Gnum * restrict edgegsttax;
  SCOTCH_Num *          fronloctax;
  Gnum                  fronlocnum;
  int                   cheklocval;
  int                   chekglbval;
  int                   procngbnum;

  Dgraph * restrict const grafptr     = (Dgraph *) orggrafptr;
  Dgraph * restrict const bandgrafptr = (Dgraph *) bndgrafptr;
  const Gnum * restrict const vertloctax = grafptr->vertloctax;
  const Gnum * restrict const vendloctax = grafptr->vendloctax;
  const Gnum * restrict const vlblloctax = grafptr->vlblloctax;
  const Gnum * restrict const veloloctax = grafptr->veloloctax;
  const Gnum * restrict const edloloctax = grafptr->edloloctax;

#ifdef SCOTCH_DEBUG_LIBRARY1
  int                   o;

  MPI_Comm_compare (((Dgraph * restrict const) orggrafptr)->proccomm,
                    ((Dgraph * restrict const) bndgrafptr)->proccomm, &o);
  if ((o != MPI_IDENT) && (o != MPI_CONGRUENT)) {
    errorPrint ("SCOTCH_dgraphBand: communicators are not congruent");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint ("SCOTCH_dgraphBand: cannot compute ghost edge array");
    return     (1);
  }

  cheklocval = 0;
  bandvnumgstsiz = MAX ((grafptr->vertgstnbr * sizeof (Gnum)), (grafptr->procglbnbr * sizeof (int))); /* TRICK: re-use array for further error collective communications */
  if ((bandvnumgsttax = memAlloc (bandvnumgstsiz)) == NULL) {
    errorPrint ("SCOTCH_dgraphBand: out of memory (1)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* This communication cannot be covered by a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("SCOTCH_dgraphBand: communication error (1)");
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
      (grafptr, fronlocnbr, fronloctab, distval, bandvnumgsttax, &bandvertlvlnum, &bandvertlocnbr, &bandedgelocsiz) != 0)
    return (1);

  bandvelolocnbr = (veloloctax != NULL) ? bandvertlocnbr : 0;
  bandedlolocsiz = (edloloctax != NULL) ? bandedgelocsiz : 0;

  bandgrafptr->flagval |= (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* Arrays created by the routine itself */
  bandgrafptr->baseval  = grafptr->baseval;

  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &bandgrafptr->procdsptab, (size_t) ((grafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &bandgrafptr->proccnttab, (size_t) (grafptr->procglbnbr       * sizeof (Gnum)),
                     &bandgrafptr->procngbtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procrcvtab, (size_t) (grafptr->procglbnbr       * sizeof (int)),
                     &bandgrafptr->procsndtab, (size_t) (grafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("SCOTCH_dgraphBand: out of memory (2)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &bandgrafptr->vertloctax, (size_t) ((bandvertlocnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &bandvlblloctax,          (size_t) (bandvertlocnbr       * sizeof (Gnum)),
                          &bandveloloctax,          (size_t) (bandvelolocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("SCOTCH_dgraphBand: out of memory (3)");
    cheklocval = 1;
  }
  else if (bandgrafptr->vertloctax -= bandgrafptr->baseval,
           bandvlblloctax          -= bandgrafptr->baseval,
           bandveloloctax = (veloloctax != NULL) ? (bandveloloctax - bandgrafptr->baseval) : NULL,
           (memAllocGroup ((void **) (void *)
                           &bandedgeloctax, (size_t) (bandedgelocsiz * sizeof (Gnum)),
                           &bandedloloctax, (size_t) (bandedlolocsiz * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("SCOTCH_dgraphBand: out of memory (4)");
    cheklocval = 1;
  }
  else {
    bandedgeloctax -= bandgrafptr->baseval;
    bandedloloctax  = (edloloctax != NULL) ? (bandedloloctax - bandgrafptr->baseval) : NULL;
  }

  if (cheklocval != 0) {                          /* In case of memory error */
    bandgrafptr->procdsptab[0] = -1;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI, /* Send received data to dummy array */
                       bandvnumgsttax + bandgrafptr->baseval, 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("SCOTCH_dgraphBand: communication error (2)");
      return     (1);
    }
    dgraphExit (bandgrafptr);
    memFree    (bandvnumgsttax + bandgrafptr->baseval);
    return     (1);
  }
  else {
    bandgrafptr->procdsptab[0] = bandvertlocnbr;
    if (MPI_Allgather (&bandgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &bandgrafptr->procdsptab[1], 1, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("SCOTCH_dgraphBand: communication error (3)");
      return     (1);
    }
  }
  bandgrafptr->procdsptab[0] = bandgrafptr->baseval; /* Build vertex-to-process array */
#ifdef SCOTCH_DEBUG_DGRAPH2
  memSet (bandvlblloctax + bandgrafptr->baseval, ~0, (bandvertlocnbr * sizeof (Gnum)));
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (procngbnum = 1; procngbnum <= grafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
    if (bandgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                       */
      dgraphExit (bandgrafptr);
      memFree    (bandvnumgsttax + bandgrafptr->baseval);
      return     (1);
    }
    bandgrafptr->procdsptab[procngbnum]    += bandgrafptr->procdsptab[procngbnum - 1];
    bandgrafptr->proccnttab[procngbnum - 1] = bandgrafptr->procdsptab[procngbnum] - bandgrafptr->procdsptab[procngbnum - 1];
  }

  fronloctax = fronloctab - bandgrafptr->baseval;
  for (vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval,
       bandvertlocnum = bandgrafptr->baseval,
       bandvertlocnnd = bandvertlocnbr + bandgrafptr->baseval,
       bandvertlocadj = bandgrafptr->procdsptab[grafptr->proclocnum] - bandgrafptr->baseval;
       bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) { /* Turn all kept graph vertices into band graph vertices */
    Gnum              vertlocnum;

    vertlocnum = fronloctax[bandvertlocnum];
    bandvlblloctax[bandvertlocnum] = (vlblloctax == NULL) ? (vertlocnum + vertlocadj) : vlblloctax[vertlocnum];
    bandvnumgsttax[vertlocnum] += bandvertlocadj; /* Turn local indices in band graph into global indices */
  }

  if (dgraphHaloSync (grafptr, (byte *) (bandvnumgsttax + bandgrafptr->baseval), GNUM_MPI) != 0) { /* Share global indexing of halo vertices */
    errorPrint ("SCOTCH_dgraphBand: cannot perform halo exchange");
    return     (1);
  }

  edgegsttax = grafptr->edgegsttax;

  veloval = 1;
  bandvertloctax = bandgrafptr->vertloctax;
  bandvelolocsum = 0;
  banddegrlocmax = 0;
  for (bandvertlocnum = bandedgelocnum = bandgrafptr->baseval; /* Build global vertex array of band graph */
       bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { /* For all vertices save for the last level  */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;

    vertlocnum = bandvlblloctax[bandvertlocnum] - vertlocadj;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (veloloctax != NULL) {
      veloval = veloloctax[vertlocnum];
      bandvelolocsum += veloval;
      bandveloloctax[bandvertlocnum] = veloval;
    }

    degrval = vendloctax[vertlocnum] - vertloctax[vertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;

    for (edgelocnum = vertloctax[vertlocnum];     /* For all original edges */
         edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (bandvnumgsttax[edgegsttax[edgelocnum]] == ~0) { /* All ends should belong to the band graph too */
        errorPrint ("SCOTCH_dgraphBand: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      bandedgeloctax[bandedgelocnum ++] = bandvnumgsttax[edgegsttax[edgelocnum]];
    }
  }
  for ( ; bandvertlocnum < bandvertlocnnd; bandvertlocnum ++) { /* For all vertices that belong to the last level */
    Gnum              vertlocnum;
    Gnum              edgelocnum;
    Gnum              degrval;

    vertlocnum = bandvlblloctax[bandvertlocnum] - vertlocadj;
    bandvertloctax[bandvertlocnum] = bandedgelocnum;
    if (veloloctax != NULL) {
      veloval = veloloctax[vertlocnum];
      bandvelolocsum += veloval;
      bandveloloctax[bandvertlocnum] = veloval;
    }

    for (edgelocnum = vertloctax[vertlocnum];     /* For all original edges */
         edgelocnum < vendloctax[vertlocnum]; edgelocnum ++) {
      Gnum              bandvertlocend;

      bandvertlocend = bandvnumgsttax[edgegsttax[edgelocnum]];
      if (bandvertlocend != ~0) {                 /* If end vertex belongs to band graph  */
        if (bandedloloctax != NULL)               /* If graph has edge weights, copy load */
          bandedloloctax[bandedgelocnum] = edloloctax[edgelocnum];
        bandedgeloctax[bandedgelocnum ++] = bandvertlocend;
      }
    }

    degrval = bandedgelocnum - bandvertloctax[bandvertlocnum];
    if (banddegrlocmax < degrval)
      banddegrlocmax = degrval;
  }
  bandvertloctax[bandvertlocnnd] = bandedgelocnum; /* Set end of vertex array */

  memFree (bandvnumgsttax + bandgrafptr->baseval);  /* Free useless space */

  if (bandedloloctax != NULL) {                   /* If graph has edge weights */
    Gnum              edgelocnum;
    Gnum              edgelocnnd;

    for (bandvertlocnum = bandgrafptr->baseval;   /* For all vertices that do not belong to the last level */
         bandvertlocnum < bandvertlvlnum; bandvertlocnum ++) { 
      Gnum              vertlocnum;
      Gnum              bandedgelocnum;

      vertlocnum     = bandvlblloctax[bandvertlocnum] - vertlocadj;
      bandedgelocnum = bandvertloctax[bandvertlocnum];
      memCpy (bandedloloctax + bandedgelocnum,    /* Copy edge load array */
              &edloloctax[vertloctax[vertlocnum]],
              (bandvertloctax[bandvertlocnum + 1] - bandedgelocnum) * sizeof (Gnum));
    }                                             /* Vertices of last level have been processed before */
  }

  bandgrafptr->procvrttab = bandgrafptr->procdsptab; /* Graph does not have holes */
  bandgrafptr->vertlocnbr = bandvertlocnbr;
  bandgrafptr->vertlocnnd = bandvertlocnbr + bandgrafptr->baseval;
  bandgrafptr->vendloctax = bandvertloctax + 1;   /* Band graph is compact */
  bandgrafptr->veloloctax = bandveloloctax;
  bandgrafptr->velolocsum = bandvelolocsum;
  bandgrafptr->vlblloctax = bandvlblloctax;
  bandgrafptr->edgeloctax = bandedgeloctax;
  bandgrafptr->edloloctax = bandedloloctax;
  bandgrafptr->edgelocnbr = bandedgelocnum - bandgrafptr->baseval;
  bandgrafptr->edgelocsiz = bandedgelocsiz;
  bandgrafptr->degrglbmax = banddegrlocmax;       /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (bandgrafptr) != 0) {
    errorPrint ("SCOTCH_dgraphBand: cannot build band graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (bandgrafptr) != 0) {
    errorPrint ("SCOTCH_dgraphBand: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}
