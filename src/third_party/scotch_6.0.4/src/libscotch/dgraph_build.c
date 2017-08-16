/* Copyright 2007,2010,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_build.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                Nicolas GICQUEL (P0.1)                  **/
/**                Jerome LACOSTE (P0.1)                   **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : These lines are the distributed source  **/
/**                graph building routines.                **/
/**                                                        **/
/**   DATES      : # Version P0.1 : from : 01 apr 1997     **/
/**                                 to   : 20 jun 1997     **/
/**                # Version P0.2 : from : 02 feb 2000     **/
/**                                 to   : 02 feb 2000     **/
/**                # Version 5.0  : from : 22 jul 2005     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to   : 03 nov 2010     **/
/**                # Version 6.0  : from : 23 dec 2013     **/
/**                                 to   : 05 jun 2014     **/
/**                                                        **/
/************************************************************/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_allreduce.h"
#include "dgraph_build.h"

/* This routine builds a distributed graph from
** the local arrays that are passed to it. If
** a vertex label array is given, it is assumed
** that edge ends are given with respect to these
** labels, and thus they are updated so as to be
** given with respect to the implicit (based)
** global numbering.
** As for all routines that build graphs, the private
** fields of the Dgraph structure have to be initialized
** if they are not already.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphBuild (
Dgraph * restrict const     grafptr,              /* Graph                                */
const Gnum                  baseval,              /* Base for indexing                    */
const Gnum                  vertlocnbr,           /* Number of local vertices             */
const Gnum                  vertlocmax,           /* Maximum number of local vertices     */
Gnum * const                vertloctax,           /* Local vertex begin array             */
Gnum * const                vendloctax,           /* Local vertex end array               */
Gnum * const                veloloctax,           /* Local vertex load array (if any)     */
Gnum * const                vnumloctax,           /* Local vertex number array (if any)   */
Gnum * const                vlblloctax,           /* Local vertex label array (if any)    */
const Gnum                  edgelocnbr,           /* Number of local edges                */
const Gnum                  edgelocsiz,           /* Size of local edge array             */
Gnum * const                edgeloctax,           /* Local edge array                     */
Gnum * const                edgegsttax,           /* Ghost edge array (if any); not const */
Gnum * const                edloloctax)           /* Local edge load array (if any)       */
{
  Gnum                vertlocnum;
  Gnum                vertlocnnd;
  Gnum                velolocsum;
  Gnum                degrlocmax;                 /* Local maximum degree */

  for (vertlocnum = baseval, vertlocnnd = vertlocnbr + baseval, degrlocmax = 0;
       vertlocnum < vertlocnnd; vertlocnum ++) {
    Gnum                degrval;

    degrval = vendloctax[vertlocnum] - vertloctax[vertlocnum];
    if (degrlocmax < degrval)
      degrlocmax = degrval;
  }

  if (veloloctax == NULL)                         /* Get local vertex load sum */
    velolocsum = vertlocnbr;
  else {
    Gnum                vertlocnum;

    for (vertlocnum = baseval, velolocsum = 0;
         vertlocnum < vertlocnnd; vertlocnum ++)
      velolocsum += veloloctax[vertlocnum];
  }

  return (dgraphBuild2 (grafptr, baseval,
                        vertlocnbr, vertlocmax, vertloctax, vendloctax, veloloctax, velolocsum, vnumloctax, vlblloctax,
                        edgelocnbr, edgelocsiz, edgeloctax, edgegsttax, edloloctax, degrlocmax));
}

/* This routine builds a distributed graph from
** the local arrays that are passed to it. If
** a vertex label array is given, it is assumed
** that edge ends are given with respect to these
** labels, and thus they are updated so as to be
** given with respect to the implicit (based)
** global numbering.
** As for all routines that build graphs, the private
** fields of the Dgraph structure have to be initialized
** if they are not already.
** These graphs do not have holes, since procvrttab
** points to procdsptab.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphBuild2 (
Dgraph * restrict const     grafptr,              /* Graph                                */
const Gnum                  baseval,              /* Base for indexing                    */
const Gnum                  vertlocnbr,           /* Number of local vertices             */
const Gnum                  vertlocmax,           /* Maximum number of local vertices     */
Gnum * const                vertloctax,           /* Local vertex begin array             */
Gnum * const                vendloctax,           /* Local vertex end array               */
Gnum * const                veloloctax,           /* Local vertex load array (if any)     */
const Gnum                  velolocsum,           /* Local sum of vertex loads            */
Gnum * const                vnumloctax,           /* Local vertex number array (if any)   */
Gnum * const                vlblloctax,           /* Local vertex label array (if any)    */
const Gnum                  edgelocnbr,           /* Number of local edges                */
const Gnum                  edgelocsiz,           /* Size of local edge array             */
Gnum * const                edgeloctax,           /* Local edge array                     */
Gnum * const                edgegsttax,           /* Ghost edge array (if any); not const */
Gnum * const                edloloctax,           /* Local edge load array (if any)       */
const Gnum                  degrlocmax)
{
  Gnum                  procnum;
  int                   reduloctab[2];
  int                   cheklocval;               /* Local consistency flag */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((vertlocmax < vertlocnbr) ||
      (edgelocsiz < edgelocnbr)) {
    errorPrint ("dgraphBuild2: invalid parameters");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  cheklocval = 0;
  if (grafptr->procdsptab == NULL) {              /* If private data not yet allocated */
    int                 procglbnbr;

    procglbnbr = grafptr->procglbnbr;
    if (memAllocGroup ((void **) (void *)         /* Allocate distributed graph private data */
                       &grafptr->procdsptab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                       &grafptr->procvrttab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                       &grafptr->proccnttab, (size_t) (procglbnbr       * sizeof (Gnum)),
                       &grafptr->procngbtab, (size_t) (procglbnbr       * sizeof (int)),
                       &grafptr->procrcvtab, (size_t) (procglbnbr       * sizeof (int)),
                       &grafptr->procsndtab, (size_t) (procglbnbr       * sizeof (int)), NULL) == NULL) {
      int *               dummtab;

      errorPrint ("dgraphBuild2: out of memory");
      if ((dummtab = memAlloc ((procglbnbr * 2) * sizeof (int))) != NULL) {
        reduloctab[0] =
        reduloctab[1] = -1;
        if (MPI_Allgather (reduloctab, 2, MPI_INT, /* Use dummy receive array (if can be allocated too) */
                           dummtab,    2, MPI_INT, grafptr->proccomm) != MPI_SUCCESS)
          errorPrint ("dgraphBuild2: communication error (1)");

        memFree (dummtab);
      }
      return (1);
    }
  }

  reduloctab[0] = (int) vertlocnbr;
  reduloctab[1] = (int) vertlocmax;
  if (MPI_Allgather (reduloctab,          2, MPI_INT, /* Use procngbtab and procrcvtab as a joint allreduce receive array */
                     grafptr->procngbtab, 2, MPI_INT, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphBuild2: communication error (2)");
    return     (1);
  }

  grafptr->procdsptab[0] =                        /* Build vertex-to-process array */
  grafptr->procvrttab[0] = baseval;
  for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
    if (grafptr->procngbtab[procnum] < 0) {       /* If error notified by another process during memory allocation */
      memFree (grafptr->procdsptab);
      grafptr->procdsptab = NULL;                 /* Free memory that has just been allocated */
      return (1);
    }

    grafptr->procdsptab[procnum + 1] = grafptr->procdsptab[procnum] + (Gnum) grafptr->procngbtab[2 * procnum];
    grafptr->procvrttab[procnum + 1] = grafptr->procvrttab[procnum] + (Gnum) grafptr->procngbtab[2 * procnum + 1];
    grafptr->proccnttab[procnum]     = grafptr->procdsptab[procnum + 1] - grafptr->procdsptab[procnum];
  }

  grafptr->flagval |= DGRAPHFREEPRIV;
  return (dgraphBuild3 (grafptr, baseval,
                        vertlocnbr, vertloctax, vendloctax, veloloctax, velolocsum, vnumloctax, vlblloctax,
                        edgelocnbr, edgelocsiz, edgeloctax, edgegsttax, edloloctax, degrlocmax));
}

/* This routine builds a distributed graph from
** the local arrays that are passed to it. If
** a vertex label array is given, it is assumed
** that edge ends are given with respect to these
** labels, and thus they are updated so as to be
** given with respect to the implicit (based)
** global numbering.
** This alternate interface assumes that the private
** fields of the Dgraph structure have already been
** initialized.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (17, 3)

int
dgraphBuild3 (
Dgraph * restrict const     grafptr,              /* Graph                                   */
const Gnum                  baseval,              /* Base for indexing                       */
const Gnum                  vertlocnbr,           /* Number of local vertices                */
Gnum * const                vertloctax,           /* Local vertex begin array                */
Gnum * const                vendloctax,           /* Local vertex end array                  */
Gnum * const                veloloctax,           /* Local vertex load array (if any)        */
const Gnum                  velolocsum,           /* Local sum of vertex loads               */
Gnum * const                vnumloctax,           /* Local vertex number array (if any)      */
Gnum * const                vlblloctax,           /* Local vertex label array (if any)       */
const Gnum                  edgelocnbr,           /* Number of local edges                   */
const Gnum                  edgelocsiz,           /* Minimum useful size of local edge array */
Gnum * const                edgeloctax,           /* Local edge array                        */
Gnum * const                edgegsttax,           /* Ghost edge array (if any); not const    */
Gnum * const                edloloctax,           /* Local edge load array (if any)          */
const Gnum                  degrlocmax)
{
  int                   procglbnbr;               /* Number of processes sharing graph data       */
  int                   procrcvnum;               /* Number of process from which to receive      */
  int                   procsndnum;               /* Number of process to which to send           */
  int                   procngbnbr;               /* Number of neighbors processed                */
  int                   procngbnum;               /* Number of current neighbor process           */
  int                   procngbsel;               /* Value of the currently used neighbor buffers */
  Gnum                  vertngbmin;               /* Smallest vertex number of neighbor process   */
  Gnum                  vertlocnum;
  Gnum                  edgelocnum;
  const Gnum *          vlbllocptr;               /* Pointer to current vertex label                  */
  DgraphLablSortVert *  vesongbptr;               /* Pointer to current sort cell                     */
  DgraphLablSortVert *  vesongbtnd;               /* Pointer to end of current sort array             */
  DgraphLablSortVert *  vesongbtab[2];            /* Neighbor vertex sorting array [norestrict:async] */
  int                   vesongbnbr[2];            /* Sizes of both vertex sort arrays                 */
  DgraphLablSortEdge *  edsoloctab;               /* Local edge sorting array                         */
  DgraphLablSortEdge *  edsoloctnd;               /* Pointer to end of edge sort array                */
  DgraphLablSortEdge *  edsolocptr;               /* Pointer to current sort edge                     */
  MPI_Request           requloctab[2];            /* Arrays for pipelined communication               */
  MPI_Status            statloctab[2];
  int                   cheklocval;               /* Local consistency flag                           */
  int                   chekglbval;               /* Global consistency flag                          */
  Gnum                  reduloctab[20];           /* Arrays for reductions                            */
  Gnum                  reduglbtab[20];

  reduloctab[0]  =   baseval;                     /* Check argument consistency */
  reduloctab[1]  = - baseval;
  reduloctab[2]  =   (veloloctax != NULL) ? 1 : 0;
  reduloctab[3]  = - reduloctab[2];
  reduloctab[4]  =   (vnumloctax != NULL) ? 1 : 0;
  reduloctab[5]  = - reduloctab[4];
  reduloctab[6]  =   (vlblloctax != NULL) ? 1 : 0;
  reduloctab[7]  = - reduloctab[6];
  reduloctab[8]  =   (edloloctax != NULL) ? 1 : 0;
  reduloctab[9]  = - reduloctab[8];
  reduloctab[10] =   (edgegsttax != NULL) ? 1 : 0;
  reduloctab[11] = - reduloctab[10];
  reduloctab[12] =   vertlocnbr;                  /* Get maximum number of local vertices */
  reduloctab[13] =   edgelocnbr;
  reduloctab[14] =   edgelocsiz;
  reduloctab[15] =   degrlocmax;
  reduloctab[16] =   (grafptr->procdsptab == NULL) ? 1 : 0; /* Error if private data not yet allocated */

  reduloctab[17] = vertlocnbr;                    /* Sum local sizes */
  reduloctab[18] = velolocsum;
  reduloctab[19] = edgelocnbr;

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 17, 3, grafptr->proccomm) != 0) {
    errorPrint ("dgraphBuild3: cannot compute reductions");
    return     (1);
  }
  if (reduglbtab[16] != 0) {
    errorPrint ("dgraphBuild3: no private data");
    return (1);
  }
  if ((reduglbtab[1]  != - reduglbtab[0]) ||
      (reduglbtab[3]  != - reduglbtab[2]) ||
      (reduglbtab[5]  != - reduglbtab[4]) ||
      (reduglbtab[7]  != - reduglbtab[6]) ||
      (reduglbtab[9]  != - reduglbtab[8]) ||
      (reduglbtab[11] != - reduglbtab[10])) {
    errorPrint ("dgraphBuild3: inconsistent parameters");
    return     (1);
  }

  grafptr->vertglbmax = reduglbtab[12];           /* Set maximum number of local vertices  */
  grafptr->edgeglbmax = reduglbtab[13];           /* Set maximum number of local edges     */
  grafptr->edgeglbsmx = reduglbtab[14];           /* Set maximum size of local edge arrays */
  grafptr->degrglbmax = reduglbtab[15];           /* Set maximum degree                    */

  grafptr->baseval    = baseval;
  grafptr->vertglbnbr = reduglbtab[17];           /* Set global and local data */
  grafptr->vertlocnbr = vertlocnbr;
  grafptr->vertlocnnd = vertlocnbr + baseval;
  grafptr->velolocsum = velolocsum;
  grafptr->veloglbsum = reduglbtab[18];
  grafptr->vertloctax = vertloctax;
  grafptr->vendloctax = vendloctax;
  grafptr->veloloctax = veloloctax;
  grafptr->vnumloctax = vnumloctax;
  grafptr->vlblloctax = vlblloctax;
  grafptr->edgeglbnbr = reduglbtab[19];
  grafptr->edgelocnbr = edgelocnbr;
  grafptr->edgelocsiz = edgelocsiz;
  grafptr->edgegsttax = edgegsttax;
  grafptr->edgeloctax = edgeloctax;
  grafptr->edloloctax = edloloctax;
#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((grafptr->procdsptab[grafptr->procglbnbr] - baseval) < grafptr->vertglbnbr) {
    errorPrint ("dgraphBuild3: invalid process vertex array");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (vlblloctax != NULL) {                       /* If vertex labels given */
    procglbnbr = grafptr->procglbnbr;

    if (memAllocGroup ((void **) (void *)
                       &vesongbtab[0], (size_t) (grafptr->vertglbmax * sizeof (DgraphLablSortVert)),
                       &vesongbtab[1], (size_t) (grafptr->vertglbmax * sizeof (DgraphLablSortVert)),
                       &edsoloctab,    (size_t) (grafptr->edgeglbmax * sizeof (DgraphLablSortEdge)),
                       NULL) == NULL) {
      errorPrint ("dgraphBuild3: out of memory");
      return     (1);
    }

    for (vertlocnum = 0, vesongbptr = vesongbtab[0], vlbllocptr = vlblloctax + baseval;
         vertlocnum < vertlocnbr;
         vertlocnum ++, vesongbptr ++, vlbllocptr ++) {
      vesongbptr->vlblglbnum = *vlbllocptr;       /* Build vertex sort array  */
      vesongbptr->vertlocnum =  vertlocnum;       /* Local index is not based */
    }
    intSort2asc1 (vesongbtab[0], vertlocnbr);
    vesongbnbr[0] = vertlocnbr;                   /* Set array size */

    cheklocval = 0;

    for (vesongbptr = vesongbtab[0] + 1, vesongbtnd = vesongbtab[0] + vertlocnbr;
         vesongbptr < vesongbtnd; vesongbptr ++) {
      if (vesongbptr[0].vlblglbnum == vesongbptr[-1].vlblglbnum) {
        cheklocval = 1;
        break;
      }
    }
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be merged with a useful one */
    MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm);
#else /* SCOTCH_DEBUG_DGRAPH1 */
    chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
    if (chekglbval != 0) {
      errorPrint ("dgraphBuild3: duplicate vertex label (1)");
      memFree    (vesongbtab[0]);
      return     (1);
    }

    for (edsolocptr = edsoloctab, edsoloctnd = edsoloctab + edgelocnbr, edgelocnum = baseval;
         edsolocptr < edsoloctnd; edsolocptr ++, edgelocnum ++) {
      edsolocptr->vlblglbnum = edgeloctax[edgelocnum];
      edsolocptr->edgelocnum = edgelocnum;
    }
    intSort2asc2 (edsoloctab, grafptr->edgelocnbr);

    procrcvnum = (grafptr->proclocnum + 1) % procglbnbr; /* Compute indices of neighbors */
    procsndnum = (grafptr->proclocnum - 1 + procglbnbr) % procglbnbr;

    for (procngbnbr = 0, procngbsel = 0;          /* For all processes */
         procngbnbr < procglbnbr;
         procngbnbr ++, procngbsel = 1 - procngbsel) {
      procngbnum = (grafptr->proclocnum + procngbnbr) % procglbnbr; /* Get neighbor process */
      vertngbmin = grafptr->procvrttab[procngbnum]; /* Get neighbor vertex number range     */

      if (procngbnbr < (procglbnbr - 1)) {        /* If not last iteration */
        MPI_Irecv (vesongbtab[1 - procngbsel], 2 * grafptr->vertglbmax, GNUM_MPI, procrcvnum, TAGVLBLLOCTAB, grafptr->proccomm, &requloctab[0]);
        MPI_Isend (vesongbtab[procngbsel], 2 * vesongbnbr[procngbsel], GNUM_MPI, procsndnum, TAGVLBLLOCTAB, grafptr->proccomm, &requloctab[1]);
      }

      if (vesongbnbr[procngbsel] > 0) {           /* If neighbor vertex sort array not empty */
        for (edsolocptr = edsoloctab,             /* Replace label by global vertex number   */
             vesongbptr = vesongbtab[procngbsel], vesongbtnd = vesongbptr + vesongbnbr[procngbsel];
             edsolocptr < edsoloctnd; ) {
          if (edsolocptr->vlblglbnum == vesongbptr->vlblglbnum) {
            if (edsolocptr->edgelocnum == -1)     /* If edge label already replaced */
              cheklocval = 1;                     /* Set error flag                 */
            else {
              edgeloctax[edsolocptr->edgelocnum] = vertngbmin + vesongbptr->vertlocnum;
              edsolocptr->edgelocnum = -1;        /* Edge has been processed */
            }
            edsolocptr ++;                        /* One more edge processed      */
            continue;                             /* Go on as quickly as possible */
          }
          if (edsolocptr->vlblglbnum < vesongbptr->vlblglbnum) {
            edsolocptr ++;                        /* One more edge processed      */
            continue;                             /* Go on as quickly as possible */
          }
          while (edsolocptr->vlblglbnum > vesongbptr->vlblglbnum) {
            if (++ vesongbptr >= vesongbtnd) {    /* Break if all labels processed */
              edsolocptr = edsoloctnd;
              break;
            }
          }
        }
      }

      if (procngbnbr < (procglbnbr - 1)) {        /* If not last iteration             */
        MPI_Waitall (2, requloctab, statloctab);  /* Wait for communication completion */
        MPI_Get_count (&statloctab[0], GNUM_MPI, &vesongbnbr[1 - procngbsel]);
        vesongbnbr[1 - procngbsel] /= 2;          /* Count items, not fields */
      }
    }

    memFree (vesongbtab[0]);

#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be merged with a useful one */
    MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm);
#else /* SCOTCH_DEBUG_DGRAPH1 */
    chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
    if (chekglbval != 0) {
      errorPrint ("dgraphBuild3: duplicate vertex label (2)");
      return     (1);
    }
  }

  return (0);
}

/* This subroutine computes the reduced values
** of all of the distributed graph fields.
** It does not deal with vertex labels, nor with
** the ghost edge array.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (4, 3)

int
dgraphBuild4 (
Dgraph * restrict const     grafptr)              /* Distributed graph */
{
  Gnum                  reduloctab[7];            /* Arrays for reductions */
  Gnum                  reduglbtab[7];

  reduloctab[0] = grafptr->vertlocnbr;            /* Get maximum over all processes */
  reduloctab[1] = grafptr->edgelocnbr;
  reduloctab[2] = grafptr->edgelocsiz;
  reduloctab[3] = grafptr->degrglbmax;            /* Here, degrglbmax may store only a local maximum degree before calling */

  reduloctab[4] = grafptr->vertlocnbr;            /* Sum local sizes */
  reduloctab[5] = grafptr->velolocsum;
  reduloctab[6] = grafptr->edgelocnbr;

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 4, 3, grafptr->proccomm) != 0) {
    errorPrint ("dgraphBuild4: cannot compute reductions");
    return     (1);
  }

  grafptr->vertglbmax = reduglbtab[0];            /* Set maximum number of local vertices  */
  grafptr->edgeglbmax = reduglbtab[1];            /* Set maximum number of local edges     */
  grafptr->edgeglbsmx = reduglbtab[2];            /* Set maximum size of local edge arrays */
  grafptr->degrglbmax = reduglbtab[3];            /* Set maximum degree                    */

  grafptr->vertglbnbr = reduglbtab[4];
  grafptr->veloglbsum = reduglbtab[5];
  grafptr->edgeglbnbr = reduglbtab[6];

#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((grafptr->procdsptab[grafptr->procglbnbr] - grafptr->baseval) < grafptr->vertglbnbr) {
    errorPrint ("dgraphBuild4: invalid process vertex array");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}
