/* Copyright 2008-2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_map_view.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted mapping routines of the libSCOTCH **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 26 jul 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY
#define LIBRARY_DGRAPH_MAP_VIEW

#include "module.h"
#include "common.h"
#include "parser.h"
#include "dgraph.h"
#include "dgraph_halo.h"
#include "arch.h"
#include "dmapping.h"
#include "kdgraph.h"
#include "library_dmapping.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mapping routines.            */
/*                                  */
/************************************/

/*+ This routine writes distributed mapping
*** statistics to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphMapView (
SCOTCH_Dgraph * const         libgrafptr,
const SCOTCH_Dmapping * const libmappptr,
FILE * const                  stream)
{
  Dgraph * restrict             grafptr;
  const LibDmapping * restrict  mappptr;
  ArchDom                       domnfrst;         /* Largest domain in architecture          */
  unsigned int * restrict       nmskloctab;       /* Local neighbor bitfield                 */
  unsigned int * restrict       nmskglbtab;       /* Local neighbor bitfield                 */
  int                           nmskidxnbr;       /* Size of bitfield; int since sent by MPI */
  Gnum * restrict               tgloloctab;       /* Local array of terminal domain loads    */
  Gnum * restrict               tgloglbtab;       /* Global array of terminal domain loads   */
  Gnum * restrict               termgsttax;       /* Terminal domain ghost mapping array     */
  Anum                          tgtnbr;           /* Number of processors in target topology */
  Anum                          tgtnum;
  Anum                          mapnbr;           /* Number of processors effectively used   */
  double                        mapavg;           /* Average mapping weight                  */
  Gnum                          mapmin;
  Gnum                          mapmax;
  Gnum                          mapsum;           /* (Partial) sum of vertex loads           */
  double                        mapdlt;
  double                        mapmmy;           /* Maximum / average ratio                 */
  Anum                          ngbsum;
  Anum                          ngbmin;
  Anum                          ngbmax;
  Gnum                          vertlocnum;
  Gnum                          veloval;
  Gnum                          edloval;
  Gnum                          commlocdist[256 + 3]; /* Array of local load distribution    */
  Gnum                          commglbdist[256 + 3];
  Gnum                          commlocload;      /* Total local edge load (edge sum)        */
  Gnum                          commlocdilat;     /* Total edge dilation                     */
  Gnum                          commlocexpan;     /* Total edge expansion                    */
  Anum                          distmax;
  Anum                          distval;
  int                           cheklocval;
  int                           chekglbval;
  DgraphHaloRequest             requdat;

  grafptr = (Dgraph *) libgrafptr;
  mappptr = (LibDmapping *) libmappptr;

  if ((grafptr->vertglbnbr == 0) ||               /* Return if nothing to do */
      (grafptr->edgeglbnbr == 0))
    return (0);

  archDomFrst (&mappptr->m.archdat, &domnfrst);   /* Get architecture domain      */
  tgtnbr = archDomSize (&mappptr->m.archdat, &domnfrst); /* Get architecture size */

  if (archVar (&mappptr->m.archdat)) {
    errorPrint ("SCOTCH_dgraphMapView: not implemented");
    return     (1);
  }

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint ("SCOTCH_dgraphMapView: cannot compute ghost edge array");
    return     (1);
  }

  nmskidxnbr = (tgtnbr + 1 + ((sizeof (int) << 3) - 1)) / (sizeof (int) << 3); /* Size of neighbor subdomain bitfield; TRICK: "+1" to have a "-1" cell for unmapped vertices */

  cheklocval = 0;
  if (memAllocGroup ((void **) (void *)
                     &nmskloctab, (size_t) (nmskidxnbr          * sizeof (unsigned int)),
                     &nmskglbtab, (size_t) (nmskidxnbr          * sizeof (unsigned int)),
                     &tgloloctab, (size_t) ((tgtnbr + 1)        * sizeof (Gnum)), /* TRICK: "+1" to have a "-1" cell for unmapped vertices */
                     &tgloglbtab, (size_t) (tgtnbr              * sizeof (Gnum)),
                     &termgsttax, (size_t) (grafptr->vertgstnbr * sizeof (Gnum)), NULL) == NULL) {
    cheklocval = 1;
  }
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("SCOTCH_dgraphMapView: communication error (1)");
    return     (1);
  }
  if (chekglbval != 0) {
    if (nmskloctab != NULL)
      memFree (nmskloctab);
    errorPrint ("SCOTCH_dgraphMapView: out of memory");
    return     (1);
  }

  if (dmapTerm (&mappptr->m, grafptr, termgsttax) != 0) {
    errorPrint ("SCOTCH_dgraphMapView: cannot build local terminal array");
    memFree    (nmskloctab);
    return     (1);
  }
  dgraphHaloAsync (grafptr, termgsttax, GNUM_MPI, &requdat);
  termgsttax -= grafptr->baseval;

  memSet (tgloloctab, 0, (tgtnbr + 1) * sizeof (Gnum));
  tgloloctab ++;                                  /* TRICK: trim array for "-1" cell */

  veloval = 1;
  for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) {
#ifdef SCOTCH_DEBUG_DMAP2
    if ((termgsttax[vertlocnum] < -1) || (termgsttax[vertlocnum] >= tgtnbr)) {
      errorPrint ("SCOTCH_dgraphMapView: invalid local terminal array");
      memFree    (nmskloctab);                      /* Free group leader */
      return     (1);
    }
#endif /* SCOTCH_DEBUG_DMAP2 */
    if (grafptr->veloloctax != NULL)
      veloval = grafptr->veloloctax[vertlocnum];
    tgloloctab[termgsttax[vertlocnum]] += veloval; /* One more vertex of given weight assigned to this target */
  }

  if (MPI_Allreduce (tgloloctab, tgloglbtab, tgtnbr, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("SCOTCH_dgraphMapView: communication error (2)");
    memFree    (nmskloctab);                      /* Free group leader */
    return     (1);
  }

  mapmin = GNUMMAX;
  mapmax = 0;
  mapsum = 0;
  mapnbr = 0;
  for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++) {
    Gnum                tgtsum;

    tgtsum = tgloglbtab[tgtnum];
    if (tgtsum != 0) {
      mapnbr ++;
      mapsum += tgtsum;
      if (tgtsum < mapmin)
        mapmin = tgtsum;
      if (tgtsum > mapmax)
        mapmax = tgtsum;
    }
  }
  mapavg = (mapnbr == 0) ? 0.0L : ((double) mapsum / (double) mapnbr);

  mapdlt = 0.0L;
  for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++)
    mapdlt += fabs ((double) tgloglbtab[tgtnum] - mapavg);
  mapdlt = (mapnbr != 0) ? mapdlt / ((double) mapnbr * mapavg) : 0.0L;
  mapmmy = (mapnbr != 0) ? (double) mapmax / (double) mapavg : 0.0L;

  if (stream != NULL) {
    fprintf (stream, "M\tProcessors " GNUMSTRING "/" GNUMSTRING "(%g)\n",
             (Gnum) mapnbr,
             (Gnum) tgtnbr,
             (double) mapnbr / (double) tgtnbr);
    fprintf (stream, "M\tTarget min=" GNUMSTRING "\tmax=" GNUMSTRING "\tavg=%g\tdlt=%g\tmaxavg=%g\n",
             (Gnum) mapmin,
             (Gnum) mapmax,
             mapavg,
             mapdlt,
             mapmmy);
  }

  if (dgraphHaloWait (&requdat) != 0) {           /* Wait for ghost terminal data to be exchanged */
    errorPrint ("SCOTCH_dgraphMapView: cannot complete asynchronous halo exchange");
    memFree    (nmskloctab);                      /* Free group leader */
    return     (1);
  }

  ngbmin = ANUMMAX;
  ngbmax = 0;
  ngbsum = 0;
  for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++) {  /* For all subdomain indices */
    int                 nmskidxnum;
    Gnum                vertlocnum;
    Anum                ngbnbr;

    if (tgloglbtab[tgtnum] <= 0)                  /* If empty subdomain, skip it */
      continue;

    memSet (nmskloctab, 0, nmskidxnbr * sizeof (int)); /* Reset neighbor bit mask */

    for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) { /* For all local vertices */
      Gnum                termnum;
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      termnum = termgsttax[vertlocnum];
      if (termnum != tgtnum)                      /* If vertex does not belong to current part or is not mapped, skip it */
        continue;

      for (edgelocnum = grafptr->vertloctax[vertlocnum], edgelocnnd = grafptr->vendloctax[vertlocnum];
           edgelocnum < edgelocnnd; edgelocnum ++) {
        Gnum                termend;

        termend = termgsttax[grafptr->edgegsttax[edgelocnum]];
        if (termend != tgtnum) {                  /* If edge is not internal             */
          termend ++;                             /* TRICK: turn unmapped to 0 and so on */
          nmskloctab[termend / (sizeof (int) << 3)] |= 1 << (termend & ((sizeof (int) << 3) - 1)); /* Flag neighbor in bit array */
        }
      }
    }
    nmskloctab[0] &= ~1;                          /* Do not account for unmapped vertices (terminal domain 0 because of "+1") */

    if (MPI_Allreduce (nmskloctab, nmskglbtab, nmskidxnbr, MPI_INT, MPI_BOR, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("SCOTCH_dgraphMapView: communication error (3)");
      memFree    (nmskloctab);                    /* Free group leader */
      return     (1);
    }

    for (nmskidxnum = 0, ngbnbr = 0; nmskidxnum < nmskidxnbr; nmskidxnum ++) {
      unsigned int        nmskbitval;

      for (nmskbitval = nmskglbtab[nmskidxnum]; nmskbitval != 0; nmskbitval >>= 1)
        ngbnbr += nmskbitval & 1;
    }

    ngbsum += ngbnbr;
    if (ngbnbr < ngbmin)
      ngbmin = ngbnbr;
    if (ngbnbr > ngbmax)
      ngbmax = ngbnbr;
  }

  if (stream != NULL) {
    fprintf (stream, "M\tNeighbors min=" GNUMSTRING "\tmax=" GNUMSTRING "\tsum=" GNUMSTRING "\n",
             (Gnum) ngbmin,
             (Gnum) ngbmax,
             (Gnum) ngbsum);
  }

  memSet (commlocdist, 0, 256 * sizeof (Gnum));   /* Initialize the data */
  commlocload  =
  commlocdilat =
  commlocexpan = 0;

  edloval = 1;
  for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) { /* For all local vertices */
    Gnum                termlocnum;
    ArchDom             termdomdat;
    Gnum                edgelocnum;
    Gnum                edgelocnnd;

    termlocnum = termgsttax[vertlocnum];
    if (termlocnum == ~0)                         /* Skip unmapped vertices */
      continue;

    archDomTerm (&mappptr->m.archdat, &termdomdat, termlocnum);

    for (edgelocnum = grafptr->vertloctax[vertlocnum], edgelocnnd = grafptr->vendloctax[vertlocnum];
         edgelocnum < edgelocnnd; edgelocnum ++) {
      ArchDom             termdomend;
      Gnum                termgstend;
      Anum                distval;

      termgstend = termgsttax[grafptr->edgegsttax[edgelocnum]];
      if (termgstend == ~0)                       /* Skip unmapped end vertices */
        continue;

      distval = 0;
      if (grafptr->edloloctax != NULL)            /* Get edge weight if any */
        edloval = grafptr->edloloctax[edgelocnum];
      if (termgstend != termlocnum) {             /* If not same domain, compute distance */
        archDomTerm (&mappptr->m.archdat, &termdomend, termgstend);
        distval = archDomDist (&mappptr->m.archdat, &termdomdat, &termdomend);
      }
      commlocdist[(distval > 255) ? 255 : distval] += edloval;
      commlocload  += edloval;
      commlocdilat += distval;
      commlocexpan += distval * edloval;
    }
  }
  commlocdist[256]     = commlocload;
  commlocdist[256 + 1] = commlocdilat;
  commlocdist[256 + 2] = commlocexpan;

  if (MPI_Allreduce (commlocdist, commglbdist, 256 + 3, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("SCOTCH_dgraphMapView: communication error (4)");
    memFree    (nmskloctab);                      /* Free group leader */
    return     (1);
  }

  if (stream != NULL) {
    Gnum                commglbload;

    commglbload = commglbdist[256];
    fprintf (stream, "M\tCommDilat=%f\t(" GNUMSTRING ")\n", /* Print expansion parameters */
           (double) commglbdist[256 + 1] / grafptr->edgeglbnbr,
           (Gnum) (commglbdist[256 + 1] / 2));
    fprintf (stream, "M\tCommExpan=%f\t(" GNUMSTRING ")\n",
             ((commglbload == 0) ? (double) 0.0L
                                 : (double) commglbdist[256 + 2] / (double) commglbload),
             (Gnum) (commglbdist[256 + 2] / 2));
    fprintf (stream, "M\tCommCutSz=%f\t(" GNUMSTRING ")\n",
             ((commglbload == 0) ? (double) 0.0L
                                 : (double) (commglbload - commglbdist[0]) / (double) commglbload),
             (Gnum) ((commglbload - commglbdist[0]) / 2));
    fprintf (stream, "M\tCommDelta=%f\n",
             (((double) commglbload  * (double) commglbdist[256 + 1]) == 0.0L)
             ? (double) 0.0L
             : ((double) commglbdist[256 + 2] * (double) grafptr->edgeglbnbr) /
               ((double) commglbload * (double) commglbdist[256 + 2]));

    for (distmax = 255; distmax != -1; distmax --)  /* Find longest distance */
      if (commglbdist[distmax] != 0)
        break;
    for (distval = 0; distval <= distmax; distval ++) /* Print distance histogram */
      fprintf (stream, "M\tCommLoad[" ANUMSTRING "]=%f\n",
               (Anum) distval,
               (double) commglbdist[distval] / (double) commglbload);
  }

  memFree (nmskloctab);                           /* Free group leader */

  return (0);
}
