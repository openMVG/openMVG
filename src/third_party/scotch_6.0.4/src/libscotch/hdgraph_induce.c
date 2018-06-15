/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hdgraph_induce.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the halo source     **/
/**                graph subgraph-making functions.        **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 19 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 27 jun 2008     **/
/**                                 to   : 22 oct 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HDGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "hdgraph.h"

/******************************/
/*                            */
/* These routines handle halo */
/* distributed source graphs. */
/*                            */
/******************************/

/* This routine builds the graph induced
** by the original graph and the list of
** selected vertices.
** The induced vnumtab array is the global
** translation of the list array if the
** original graph does not have a vnumtab,
** or the proper subset of the original
** vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hdgraphInduceList (
Hdgraph * restrict const    orggrafptr,
const Gnum                  indlistnbr,
const Gnum * restrict const indlisttab,           /* Local list of kept vertices */
Hdgraph * restrict const    indgrafptr)
{
  const Gnum * restrict orgvertloctax;
  const Gnum * restrict orgvendloctax;
  const Gnum * restrict orgveloloctax;
  const Gnum * restrict orgedgegsttax;
  const Gnum * restrict orgedgeloctax;
  Gnum * restrict       orgindxgsttax;            /* Based access to vertex translation array       */
  Gnum * restrict       orgindxhaltax;            /* Based access to halo vertex translation array  */
  Gnum * restrict       indvertloctax;
  Gnum * restrict       indveloloctax;
  Gnum                  indvertlocnnd;            /* Based index of end of local vertex array       */
  Gnum                  indvertlocnum;            /* Number of current vertex in induced graph      */
  Gnum                  indvertglbnum;            /* Number of current vertex in global ordering    */
  Gnum * restrict       indvendloctax;
  Gnum                  indvelolocnbr;            /* Size of local vertex load array                */
  Gnum                  indvelolocsum;            /* Sum of vertex loads                            */
  Gnum                  indvhallocnum;            /* Number of halo vertex to be declared           */
  Gnum * restrict       indedgeloctax;
  Gnum                  indedgelocmax;            /* (Approximate) number of edges in induced graph */
  Gnum                  indedgelocsiz;            /* Real size of edge array, including halo        */
  Gnum                  indedgelocnbr;            /* Real number of edges in induced graph          */
  Gnum                  indedgelocnum;
  Gnum                  inddegrlocmax;            /* Local maximum degree over non-halo vertices    */
  const Gnum * restrict indlisttax;               /* Based access to list of kept vertices          */
  int                   cheklocval;
  int                   chekglbval;

  if (dgraphGhst (&orggrafptr->s) != 0) {         /* Compute ghost edge array if not already present */
    errorPrint ("hdgraphInduceList: cannot compute ghost edge array");
    return     (1);
  }

  memSet (indgrafptr, 0, sizeof (Hdgraph));       /* Pre-initialize graph fields */

  indgrafptr->s.proccomm   = orggrafptr->s.proccomm;
  indgrafptr->s.procglbnbr = orggrafptr->s.procglbnbr;
  indgrafptr->s.proclocnum = orggrafptr->s.proclocnum;
  indgrafptr->s.flagval    = (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP; /* For premature freeing on error; do not free vhndloctab as it is grouped with vertloctab */

  if (orggrafptr->s.veloloctax != NULL) {
    indvelolocnbr = indlistnbr;
    indvelolocsum = 0;
  }
  else {
    indvelolocnbr = 0;
    indvelolocsum = indlistnbr;
  }
  indedgelocmax = orggrafptr->s.edgelocnbr;       /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->s.degrglbmax > 0) && (indlistnbr < (indedgelocmax / orggrafptr->s.degrglbmax)))
    indedgelocmax = indlistnbr * orggrafptr->s.degrglbmax;
  indedgelocmax += orggrafptr->ehallocnbr;

  cheklocval =
  chekglbval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &indgrafptr->s.procdsptab, (size_t) ((orggrafptr->s.procglbnbr + 1) * sizeof (Gnum)),
                     &indgrafptr->s.proccnttab, (size_t) (orggrafptr->s.procglbnbr       * sizeof (Gnum)),
                     &indgrafptr->s.procngbtab, (size_t) (orggrafptr->s.procglbnbr       * sizeof (int)),
                     &indgrafptr->s.procrcvtab, (size_t) (orggrafptr->s.procglbnbr       * sizeof (int)),
                     &indgrafptr->s.procsndtab, (size_t) (orggrafptr->s.procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("hdgraphInduceList: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &indgrafptr->s.vertloctax, (size_t) ((indlistnbr + 1) * sizeof (Gnum)), /* Compact vertex arrays                  */
                          &indgrafptr->s.vendloctax, (size_t) (indlistnbr       * sizeof (Gnum)), /* Vertex end array for non-halo vertices */
                          &indgrafptr->s.vnumloctax, (size_t) (indlistnbr       * sizeof (Gnum)),
                          &indgrafptr->s.veloloctax, (size_t) (indvelolocnbr    * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("hdgraphInduceList: out of memory (2)");
    cheklocval = 1;
  }
  else if (indgrafptr->s.vertloctax -= orggrafptr->s.baseval,
           indgrafptr->s.vendloctax -= orggrafptr->s.baseval,
           indgrafptr->s.vnumloctax -= orggrafptr->s.baseval,
           indgrafptr->s.veloloctax  = ((orggrafptr->s.veloloctax != NULL) ? indgrafptr->s.veloloctax - orggrafptr->s.baseval : NULL),
           memAllocGroup ((void **) (void *)
                          &indgrafptr->s.edgeloctax, (size_t) (indedgelocmax            * sizeof (Gnum)), /* Pre-allocate space for edgeloctab              */
                          &orgindxgsttax,            (size_t) (orggrafptr->s.vertgstnbr * sizeof (Gnum)), /* orgindxgsttab and orgindxhaltab are at the end */
                          &orgindxhaltax,            (size_t) (orggrafptr->vhallocnbr   * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("hdgraphInduceList: out of memory (3)");
    cheklocval = 1;
  }
  else
    indgrafptr->s.edgeloctax -= orggrafptr->s.baseval;

  if (cheklocval != 0) {                          /* In case of memory error */
    Gnum                procngbnum;
    int                 dummyval;

    dummyval   = -1;
    chekglbval = 1;
    if (MPI_Allgather (&dummyval, 1, GNUM_MPI,    /* Use proccnttab of orggraf as dummy receive array (will be regenerated) */
                       orggrafptr->s.proccnttab, 1, GNUM_MPI, indgrafptr->s.proccomm) != MPI_SUCCESS)
      errorPrint ("hdgraphInduceList: communication error (1)");

    for (procngbnum = 1; procngbnum <= orggrafptr->s.procglbnbr; procngbnum ++) /* Rebuild proccnttab of orggraf */
      orggrafptr->s.proccnttab[procngbnum - 1] = orggrafptr->s.procdsptab[procngbnum] - orggrafptr->s.procdsptab[procngbnum - 1];
  }
  else {
    indgrafptr->s.procvrttab = indgrafptr->s.procdsptab; /* Graph does not have holes */
    indgrafptr->s.procdsptab[0] = indlistnbr;
    if (MPI_Allgather (&indgrafptr->s.procdsptab[0], 1, GNUM_MPI,
                       &indgrafptr->s.proccnttab[0], 1, GNUM_MPI, indgrafptr->s.proccomm) != MPI_SUCCESS) {
      errorPrint ("hdgraphInduceList: communication error (2)");
      chekglbval = 1;
    }
    else {
      Gnum                procngbnum;

      indgrafptr->s.procdsptab[0] = orggrafptr->s.baseval; /* Build vertex-to-process array                                              */
      for (procngbnum = 0; procngbnum < indgrafptr->s.procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
        if (indgrafptr->s.proccnttab[procngbnum] < 0) { /* If error notified by another process                                          */
          chekglbval = 1;
          break;
        }
        indgrafptr->s.procdsptab[procngbnum + 1] = indgrafptr->s.procdsptab[procngbnum] + indgrafptr->s.proccnttab[procngbnum];
      }
    }
  }
  if (chekglbval != 0) {                          /* If something went wrong in all of the above */
    hdgraphExit (indgrafptr);
    return      (1);
  }

  memSet (orgindxgsttax, ~0, orggrafptr->s.vertgstnbr * sizeof (Gnum)); /* Preset index arrays */
  orgindxgsttax -= orggrafptr->s.baseval;
  memSet (orgindxhaltax, ~0, orggrafptr->vhallocnbr * sizeof (Gnum));
  orgindxhaltax -= orggrafptr->s.baseval;

  indlisttax    = indlisttab - orggrafptr->s.baseval;
  indvertlocnnd = indlistnbr + orggrafptr->s.baseval;
  for (indvertlocnum = orggrafptr->s.baseval, indvertglbnum = indgrafptr->s.procdsptab[indgrafptr->s.proclocnum]; /* Set adjustment for global ordering */
       indvertlocnum < indvertlocnnd; indvertlocnum ++, indvertglbnum ++)
    orgindxgsttax[indlisttax[indvertlocnum]] = indvertglbnum; /* Mark selected vertices */

  if (dgraphHaloSync (&orggrafptr->s, (byte *) (orgindxgsttax + orggrafptr->s.baseval), GNUM_MPI) != 0) { /* Share global indexing of subgraph vertices */
    errorPrint  ("hdgraphInduceList: cannot perform halo exchange");
    hdgraphExit (indgrafptr);
    return      (1);
  }

  orgvertloctax = orggrafptr->s.vertloctax;
  orgvendloctax = orggrafptr->s.vendloctax;
  orgveloloctax = orggrafptr->s.veloloctax;
  orgedgegsttax = orggrafptr->s.edgegsttax;
  orgedgeloctax = orggrafptr->s.edgeloctax;
  indvertloctax = indgrafptr->s.vertloctax;
  indvendloctax = indgrafptr->s.vendloctax;
  indveloloctax = indgrafptr->s.veloloctax;
  indedgeloctax = indgrafptr->s.edgeloctax;
  inddegrlocmax = 0;
  for (indvertlocnum = indedgelocnum = indvhallocnum = orggrafptr->s.baseval, indedgelocnbr = 0;
       indvertlocnum < indvertlocnnd; indvertlocnum ++) {
    Gnum                orgvertlocnum;
    Gnum                orgedgelocnum;
    Gnum                orgdegrlocval;
    Gnum                inddegrlocval;
    Gnum                indedgelocnnd;

    orgvertlocnum = indlisttax[indvertlocnum];
    orgdegrlocval = orgvendloctax[orgvertlocnum] - orgvertloctax[orgvertlocnum];

    indvertloctax[indvertlocnum] = indedgelocnum;
    if (orgveloloctax != NULL) {                  /* If graph has vertex weights */
      indvelolocsum +=                            /* Accumulate vertex loads     */
      indveloloctax[indvertlocnum] = orgveloloctax[orgvertlocnum];
    }
    indedgelocnnd = indedgelocnum + orgdegrlocval;
#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (indedgelocnnd > (indedgelocmax + orggrafptr->s.baseval)) {
      errorPrint  ("hdgraphInduceList: internal error (1)");
      return      (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
    for (orgedgelocnum = orgvertloctax[orgvertlocnum]; /* Process local and ghost non-halo vertices */
         orgedgelocnum < orgvendloctax[orgvertlocnum]; orgedgelocnum ++) {
      Gnum                orgvertlocend;
      Gnum                indvertgstend;

      orgvertlocend = orgedgegsttax[orgedgelocnum];
#ifdef SCOTCH_DEBUG_HDGRAPH2
      if ((orgvertlocend < orggrafptr->s.baseval) || (orgvertlocend > orggrafptr->s.vertgstnnd)) {
        errorPrint  ("hdgraphInduceList: internal error (2)");
        return      (1);
      }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
      indvertgstend = orgindxgsttax[orgvertlocend];
      if (indvertgstend >= 0)                    /* If edge is local or halo                           */
        indedgeloctax[indedgelocnum ++] = indvertgstend; /* Keep it as regular edge                    */
      else {                                     /* If edge is halo edge                               */
        if (indvertgstend == ~0)                 /* If halo vertex not assigned yet                    */
          orgindxgsttax[orgvertlocend] = indvertgstend = -2 - indvhallocnum ++; /* Set new halo number */
        indedgeloctax[-- indedgelocnnd] = -2 - indvertgstend;
      }
    }
#ifdef SCOTCH_DEBUG_HDGRAPH2
    if (indedgelocnnd != indedgelocnum) {
      errorPrint  ("hdgraphInduceList: internal error (3)");
      return      (1);
    }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
    indvendloctax[indvertlocnum] = indedgelocnum;
    inddegrlocval  = indedgelocnum - indvertloctax[indvertlocnum];
    indedgelocnbr += inddegrlocval;
    if (inddegrlocmax < inddegrlocval)
      inddegrlocmax = inddegrlocval;

    for (indedgelocnum = indvertloctax[indvertlocnum] + orgdegrlocval; /* Process local halo vertices */
         orgedgelocnum < orggrafptr->vhndloctax[orgvertlocnum]; orgedgelocnum ++) {
      Gnum                orgvhallocend;
      Gnum                indvhallocend;

      orgvhallocend = orgedgeloctax[orgedgelocnum]; /* Halo vertices only exist in the edgeloctab array */
#ifdef SCOTCH_DEBUG_HDGRAPH2
      if ((orgvhallocend < orggrafptr->s.baseval) || (orgvhallocend >= (orggrafptr->vhallocnbr + orggrafptr->s.baseval))) {
        errorPrint  ("hdgraphInduceList: inconsistent halo vertex numbers");
        hdgraphExit (indgrafptr);
        return      (1);
      }
      if (indedgelocnum >= (indedgelocmax + orggrafptr->s.baseval)) {
        errorPrint  ("hdgraphInduceList: internal error (4)");
        return      (1);
      }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
      indvhallocend = orgindxhaltax[orgvhallocend];
      if (indvhallocend == ~0)                    /* If halo vertex not assigned yet            */
        orgindxhaltax[orgvhallocend] = indvhallocend = indvhallocnum ++; /* Set new halo number */
      indgrafptr->s.edgeloctax[indedgelocnum ++] = indvhallocend;
    }
  }
  indvertloctax[indvertlocnum] = indedgelocnum;   /* Mark end of edge array for vhndloctax                 */
  indedgelocsiz = indedgelocnum - orggrafptr->s.baseval; /* Global number of edges, both non-halo and halo */

  indgrafptr->s.edgeloctax = memRealloc (indgrafptr->s.edgeloctax + orggrafptr->s.baseval,
                                         (size_t) (indedgelocsiz * sizeof (Gnum)));
  indgrafptr->s.edgeloctax -= orggrafptr->s.baseval;

  if (orggrafptr->s.vnumloctax != NULL) {         /* Adjust vnumloctax */
    for (indvertlocnum = orggrafptr->s.baseval; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->s.vnumloctax[indvertlocnum] = orggrafptr->s.vnumloctax[indlisttax[indvertlocnum]];
  }
  else {
    Gnum                orgvertglbadj;

    orgvertglbadj = orggrafptr->s.procvrttab[orggrafptr->s.proclocnum] - orggrafptr->s.baseval; /* Set adjustement for global ordering */
    for (indvertlocnum = orggrafptr->s.baseval; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indgrafptr->s.vnumloctax[indvertlocnum] = indlisttax[indvertlocnum] + orgvertglbadj;
  }
  indgrafptr->vhallocnbr = indvhallocnum - orggrafptr->s.baseval;
  indgrafptr->vhndloctax = indgrafptr->s.vertloctax + 1; /* Compact edge array with halo vertices   */
  indgrafptr->ehallocnbr = indedgelocsiz - indedgelocnbr; /* Get number of halo edges by difference */
  indgrafptr->levlnum    = orggrafptr->levlnum + 1; /* Induced subgraph is one level below          */

  indgrafptr->s.baseval    = orggrafptr->s.baseval;
  indgrafptr->s.vertlocnbr = indlistnbr;
  indgrafptr->s.vertlocnnd = indlistnbr + orggrafptr->s.baseval;
  indgrafptr->s.velolocsum = indvelolocsum;
  indgrafptr->s.edgelocnbr = indedgelocnbr;
  indgrafptr->s.edgelocsiz = indedgelocsiz;
  indgrafptr->s.degrglbmax = orggrafptr->s.degrglbmax;
  if (dgraphBuild4 (&indgrafptr->s) != 0) {
    errorPrint ("hdgraphInduceList: cannot build induced graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_HDGRAPH2
  if (hdgraphCheck (indgrafptr) != 0) {           /* Check graph consistency */
    errorPrint  ("hdgraphInduceList: internal error (5)");
    hdgraphExit (indgrafptr);
    return      (1);
  }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */

  return (0);
}
