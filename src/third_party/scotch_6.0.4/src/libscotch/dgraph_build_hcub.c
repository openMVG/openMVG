/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_build_hcub.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : These lines are the distributed source  **/
/**                graph building routines for hypercube   **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version P0.1 : from : 19 may 1999     **/
/**                                 to   : 19 may 1999     **/
/**                # Version P0.2 : from : 02 feb 2000     **/
/**                                 to   : 02 feb 2000     **/
/**                # Version 5.0  : from : 20 jul 2005     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 12 nov 2008     **/
/**                                 to   : 12 nov 2008     **/
/**                                                        **/
/************************************************************/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"

/*******************************************/
/*                                         */
/* The distributed graph building routine. */
/*                                         */
/*******************************************/

/* This routine builds a distributed hypercube of
** given dimension.
** Since this routine calls dgraphBuild, the private
** data of the Dgraph structure will be initialized
** by the latter, if needed.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphBuildHcub (
Dgraph * restrict const     grafptr,              /* Graph               */
const Gnum                  hcubdim,              /* Hypercube dimension */
const Gnum                  baseval,              /* Base value          */
const Gnum                  flagval)              /* Flags               */
{
  Gnum                procngbnum;
  Gnum                vertglbnbr;                 /* Total number of vertices        */
  Gnum                vertglbnum;                 /* Global number of current vertex */
  Gnum                vertlocnbr;                 /* Number of local vertices        */
  Gnum                velolocnbr;                 /* Number of local vertex loads    */
  Gnum                vertlocnnd;
  Gnum                vertlocnum;
  Gnum * restrict     vertloctax;
  Gnum * restrict     veloloctax;
#ifdef SCOTCH_DEBUG_DGRAPH3
  Gnum * restrict     vlblloctax;
#endif /* SCOTCH_DEBUG_DGRAPH3 */
  Gnum                edgelocnbr;
  Gnum * restrict     edgeloctax;
  Gnum                edgelocnum;
  Gnum                edlolocnbr;
  Gnum * restrict     edloloctax;
  int                 cheklocval;
  Gnum                reduloctab[7];
  Gnum                reduglbtab[7];

  vertglbnbr = 1 << hcubdim;
  vertlocnbr = DATASIZE (vertglbnbr, grafptr->procglbnbr, grafptr->proclocnum);
  velolocnbr = ((flagval & 1) != 0) ? vertlocnbr : 0;
  edgelocnbr = vertlocnbr * hcubdim;              /* Set local number of arcs */
  edlolocnbr = ((flagval & 2) != 0) ? edgelocnbr : 0;

  for (procngbnum = 0, vertglbnum = 0;            /* Compute index of first local vertex */
       procngbnum < grafptr->proclocnum; procngbnum ++)
    vertglbnum += DATASIZE (vertglbnbr, grafptr->procglbnbr, procngbnum);

  cheklocval = 0;
  vertloctax =
  edgeloctax = NULL;
  if (memAllocGroup ((void **) (void *)
                     &vertloctax, (size_t) ((vertlocnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
#ifdef SCOTCH_DEBUG_DGRAPH3
                     &vlblloctax, (size_t) (vertlocnbr       * sizeof (Gnum)),
#endif /* SCOTCH_DEBUG_DGRAPH3 */
                     &veloloctax, (size_t) (vertlocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphBuildHcub: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)
                          &edgeloctax, (size_t) (edgelocnbr * sizeof (Gnum)),
                          &edloloctax, (size_t) (edlolocnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphBuildHcub: out of memory (2)");
    cheklocval = 1;
  }
  reduloctab[0] =   hcubdim;
  reduloctab[1] = - hcubdim;
  reduloctab[2] =   baseval;
  reduloctab[3] = - baseval;
  reduloctab[4] =   flagval;
  reduloctab[5] = - flagval;
  reduloctab[6] =   cheklocval;

  if (MPI_Allreduce (reduloctab, reduglbtab, 7, GNUM_MPI, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphBuildHcub: communication error");
    return     (1);
  }
  if (reduglbtab[6] != 0) {
    if (vertloctax != NULL) {
      if (edgeloctax != NULL)
        memFree (edgeloctax);
      memFree (vertloctax);
    }
    return (1);
  }
  if ((reduglbtab[1] != - reduglbtab[0]) ||
      (reduglbtab[3] != - reduglbtab[2]) ||
      (reduglbtab[5] != - reduglbtab[4])) {
    errorPrint ("dgraphBuildHcub: inconsistent parameters");
    return     (1);
  }
  vertloctax -= baseval;
  veloloctax  = ((flagval & 1) != 0) ? (veloloctax - baseval) : NULL;
  edgeloctax -= baseval;
  edloloctax  = ((flagval & 2) != 0) ? (edloloctax - baseval) : NULL;
#ifdef SCOTCH_DEBUG_DGRAPH3
  vlblloctax -= baseval;
#endif /* SCOTCH_DEBUG_DGRAPH3 */

  for (vertlocnum = baseval, vertlocnnd = vertlocnbr + baseval, edgelocnum = baseval;
       vertlocnum < vertlocnnd; vertlocnum ++, vertglbnum ++) {
    Gnum                vertngbbit;               /* Bit that differs between neighbors */

    if (veloloctax != NULL)
      veloloctax[vertlocnum] = 1 + (vertglbnum & 3);  /* Pseudo random weight (1 to 5) */
#ifdef SCOTCH_DEBUG_DGRAPH3
    vlblloctax[vertlocnum] = ((vertglbnum * COARHASHPRIME) % vertglbnbr) + baseval; /* Hash vertices to spread labels */
#endif /* SCOTCH_DEBUG_DGRAPH3 */
    vertloctax[vertlocnum] = edgelocnum;

    for (vertngbbit = 1; vertngbbit < vertglbnbr; vertngbbit <<= 1, edgelocnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH3
      edgeloctax[edgelocnum] = (((vertglbnum ^ vertngbbit) * COARHASHPRIME) % vertglbnbr) + baseval;
#else /* SCOTCH_DEBUG_DGRAPH3 */
      edgeloctax[edgelocnum] = (vertglbnum ^ vertngbbit) + baseval;
#endif /* SCOTCH_DEBUG_DGRAPH3 */
      if (edloloctax != NULL)
        edloloctax[edgelocnum] = ((vertglbnum + edgeloctax[edgelocnum]) % 16) + 1; /* Pseudo random weight (1 to 16) */
    }
  }
  vertloctax[vertlocnum] = edgelocnum;            /* Mark end of local vertex array */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (edgelocnum != edgelocnbr + baseval) {
    errorPrint ("dgraphBuildHcub: internal error");
    memFree    (vertloctax + baseval);           /* Free memory group leader */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (dgraphBuild2 (grafptr, baseval,             /* Build the distributed graph */
#ifdef SCOTCH_DEBUG_DGRAPH3
                    vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, NULL, vertlocnbr, NULL, vlblloctax,
#else /* SCOTCH_DEBUG_DGRAPH3 */
                    vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, NULL, vertlocnbr, NULL, NULL,
#endif /* SCOTCH_DEBUG_DGRAPH3 */
                    edgelocnbr, edgelocnbr, edgeloctax, NULL, edloloctax, hcubdim) != 0) {
    memFree (edgeloctax + baseval);           /* Free memory group leaders */
    memFree (vertloctax + baseval);
    return  (1);
  }

  grafptr->flagval |= DGRAPHFREETABS | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* Arrays created by the routine itself */

  return (0);
}
