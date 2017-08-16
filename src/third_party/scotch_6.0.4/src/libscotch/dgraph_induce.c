/* Copyright 2007-2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_induce.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                subgraph-making functions.              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 apr 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 31 mar 2008     **/
/**                                 to   : 30 jul 2010     **/
/**                # Version 6.0  : from : 29 aug 2012     **/
/**                                 to   : 13 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DGRAPH
#define DGRAPH_INDUCE

#include "module.h"
#include "common.h"
#include "dgraph.h"

/****************************************/
/*                                      */
/* These routines handle source graphs. */
/*                                      */
/****************************************/

int
dgraphInduce2 (
Dgraph * restrict const       orggrafptr,
Gnum                       (* orgfuncptr) (Dgraph * restrict const, Dgraph * restrict const, const void * restrict const, Gnum * restrict const),
const void * const            orgdataptr,         /* Pointer to routine-specific data                 */
const Gnum                    indvertlocnbr,      /* Number of vertices in induced subgraph           */
Gnum *                        indvnumloctmp,      /* Pointer to temporary index array; TRICK: [alias] */
Dgraph * restrict const       indgrafptr)
{
  Gnum * restrict       orgindxgsttax;            /* Based access to vertex translation array       */
  Gnum                  indvertlocnnd;            /* Based index of end of local vertex array       */
  Gnum                  indvertlocnum;            /* Number of current vertex in induced graph      */
  Gnum                  indvertglbnum;            /* Number of current vertex in global ordering    */
  Gnum                  indvelolocnbr;            /* Size of local vertex load array                */
  Gnum                  indvelolocsum;            /* Sum of vertex loads                            */
  Gnum *                indvnumloctax;            /* TRICK: maybe alias of indvnumloctmp            */
  Gnum                  indvlbllocnbr;            /* Size of local vertex label array               */
  Gnum                  indedgelocmax;            /* (Approximate) number of edges in induced graph */
  Gnum                  indedgelocnbr;            /* Real number of edges in induced graph          */
  Gnum                  indedgelocnum;
  Gnum * restrict       indedloloctax;
  Gnum                  inddegrlocmax;            /* Local maximum degree                           */
  const Gnum * restrict indlisttax;
  Gnum                  baseval;
  int                   cheklocval;
  int                   chekglbval;

  const Gnum * restrict const orgvertloctax = orggrafptr->vertloctax;
  const Gnum * restrict const orgvendloctax = orggrafptr->vendloctax;
  const Gnum * restrict const orgvnumloctax = orggrafptr->vnumloctax;
  const Gnum * restrict const orgvlblloctax = orggrafptr->vlblloctax;
  const Gnum * restrict const orgveloloctax = orggrafptr->veloloctax;
  const Gnum * restrict const orgedloloctax = orggrafptr->edloloctax;

  if (dgraphGhst (orggrafptr) != 0) {             /* Compute ghost edge array if not already present */
    errorPrint ("dgraphInduce2: cannot compute ghost edge array");
    return     (1);
  }

  baseval                = orggrafptr->baseval;
  indgrafptr->flagval   |= (DGRAPHFREEALL ^ DGRAPHFREECOMM) | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;
  indgrafptr->baseval    = baseval;
  indgrafptr->vertlocnbr = indvertlocnbr;         /* Must be set before orgfuncptr() is called */
  indgrafptr->vertlocnnd = indvertlocnbr + baseval;

  if (orgveloloctax != NULL) {
    indvelolocnbr = indvertlocnbr;
    indvelolocsum = 0;
  }
  else {
    indvelolocnbr = 0;
    indvelolocsum = indvertlocnbr;
  }
  indvlbllocnbr = (orgvlblloctax != NULL) ? indvertlocnbr : 0;
  indedgelocmax = orggrafptr->edgelocnbr;         /* Choose best upper bound on number of edges (avoid multiply overflow) */
  if ((orggrafptr->degrglbmax > 0) && (indvertlocnbr < (indedgelocmax / orggrafptr->degrglbmax)))
    indedgelocmax = indvertlocnbr * orggrafptr->degrglbmax;
  if (orggrafptr->edloloctax != NULL)             /* If graph has edge weights */
    indedgelocmax *= 2;                           /* Account for edge weights  */

  cheklocval =
  chekglbval = 0;
  if (memAllocGroup ((void **) (void *)           /* Allocate distributed graph private data */
                     &indgrafptr->procdsptab, (size_t) ((orggrafptr->procglbnbr + 1) * sizeof (Gnum)),
                     &indgrafptr->proccnttab, (size_t) (orggrafptr->procglbnbr       * sizeof (Gnum)),
                     &indgrafptr->procngbtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procrcvtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)),
                     &indgrafptr->procsndtab, (size_t) (orggrafptr->procglbnbr       * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphInduce2: out of memory (1)");
    cheklocval = 1;
  }
  else if (memAllocGroup ((void **) (void *)      /* Allocate distributed graph public data */
                          &indgrafptr->vertloctax, (size_t) ((indvertlocnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                          &indgrafptr->vnumloctax, (size_t) (indvertlocnbr       * sizeof (Gnum)),
                          &indgrafptr->veloloctax, (size_t) (indvelolocnbr       * sizeof (Gnum)),
                          &indgrafptr->vlblloctax, (size_t) (indvlbllocnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphInduce2: out of memory (2)");
    cheklocval = 1;
  }
  else if (indgrafptr->vertloctax -= baseval,
           indgrafptr->vnumloctax -= baseval,
           indgrafptr->veloloctax  = (orgveloloctax != NULL) ? indgrafptr->veloloctax - baseval : NULL,
           indgrafptr->vlblloctax  = indgrafptr->vlblloctax - baseval, /* If no vertex labels, vlblloctax will point to vnumloctax afterward */
           memAllocGroup ((void **) (void *)
                          &indgrafptr->edgeloctax, (size_t) (indedgelocmax          * sizeof (Gnum)), /* Pre-allocate space for edgetab (and edlotab) */
                          &orgindxgsttax,          (size_t) (orggrafptr->vertgstnbr * sizeof (Gnum)), NULL) == NULL) { /* orgindxgsttab is at the end */
    errorPrint ("dgraphInduce2: out of memory (3)");
    cheklocval = 1;
  }
  else
    indgrafptr->edgeloctax -= baseval;

  if (cheklocval != 0) {                          /* In case of memory error */
    Gnum                procngbnum;
    Gnum                dummyval;

    dummyval   = -1;
    chekglbval = 1;
    if (MPI_Allgather (&dummyval, 1, GNUM_MPI,    /* Use proccnttab of orggraf as dummy receive array (will be regenerated) */
                       orggrafptr->proccnttab, 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS)
      errorPrint ("dgraphInduce2: communication error (1)");

    for (procngbnum = 1; procngbnum <= orggrafptr->procglbnbr; procngbnum ++) /* Rebuild proccnttab of orggraf */
      orggrafptr->proccnttab[procngbnum - 1] = orggrafptr->procdsptab[procngbnum] - orggrafptr->procdsptab[procngbnum - 1];
  }
  else {
    indgrafptr->procdsptab[0] = indvertlocnbr;
    if (MPI_Allgather (&indgrafptr->procdsptab[0], 1, GNUM_MPI,
                       &indgrafptr->proccnttab[0], 1, GNUM_MPI, indgrafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphInduce2: communication error (2)");
      chekglbval = 1;
    }
    else {
      Gnum                procngbnum;

      indgrafptr->procdsptab[0] = baseval;        /* Build vertex-to-process array                                                     */
      for (procngbnum = 0; procngbnum < indgrafptr->procglbnbr; procngbnum ++) { /* Process potential error flags from other processes */
        if (indgrafptr->procdsptab[procngbnum] < 0) { /* If error notified by another process                                          */
          chekglbval = 1;
          break;
        }
        indgrafptr->procdsptab[procngbnum + 1] = indgrafptr->procdsptab[procngbnum] + indgrafptr->proccnttab[procngbnum];
      }
    }
    indgrafptr->procvrttab = indgrafptr->procdsptab; /* Graph does not have holes */
  }
  if (chekglbval != 0) {                          /* If something went wrong in all of the above */
    dgraphFree (indgrafptr);
    return     (1);
  }

  memSet (orgindxgsttax, ~0, orggrafptr->vertlocnbr * sizeof (Gnum)); /* Preset index array */
  orgindxgsttax -= baseval;

  indedgelocmax = orgfuncptr (indgrafptr, orggrafptr, orgdataptr, orgindxgsttax); /* Call flagging subroutine */

  if (dgraphHaloSync (orggrafptr, (byte *) (orgindxgsttax + baseval), GNUM_MPI) != 0) { /* Share global indexing of subgraph vertices */
    errorPrint ("dgraphInduce2: cannot perform halo exchange");
    dgraphFree (indgrafptr);
    return     (1);
  }

  if (indvnumloctmp == NULL)                      /* indgrafptr->vnumloctax did not exist when function was called */
    indvnumloctmp = indgrafptr->vnumloctax;

  indedloloctax = (orggrafptr->edloloctax != NULL) ? indgrafptr->edgeloctax + indedgelocmax : NULL;
  inddegrlocmax = 0;
  for (indvertlocnum = indedgelocnum = baseval, indvertlocnnd = indvertlocnbr + baseval;
       indvertlocnum < indvertlocnnd; indvertlocnum ++) {
    Gnum                orgvertlocnum;
    Gnum                orgedgelocnum;

    orgvertlocnum = indvnumloctmp[indvertlocnum];
    indgrafptr->vertloctax[indvertlocnum] = indedgelocnum;
    if (orgveloloctax != NULL) {                  /* If graph has vertex weights */
      indvelolocsum +=                            /* Accumulate vertex loads     */
      indgrafptr->veloloctax[indvertlocnum] = orgveloloctax[orgvertlocnum];
    }
    if (orgvlblloctax != NULL)                    /* If graph has vertex labels */
      indgrafptr->vlblloctax[indvertlocnum] = orgvlblloctax[orgvertlocnum];

    for (orgedgelocnum = orgvertloctax[orgvertlocnum];
         orgedgelocnum < orgvendloctax[orgvertlocnum]; orgedgelocnum ++) {
      Gnum                indvertgstend;

      indvertgstend = orgindxgsttax[orggrafptr->edgegsttax[orgedgelocnum]];
      if (indvertgstend != ~0) {                  /* If edge should be kept */
        indgrafptr->edgeloctax[indedgelocnum] = indvertgstend;
        if (indedloloctax != NULL)
          indedloloctax[indedgelocnum] = orgedloloctax[orgedgelocnum];
        indedgelocnum ++;
      }
    }
    if (inddegrlocmax < (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]))
      inddegrlocmax = (indedgelocnum - indgrafptr->vertloctax[indvertlocnum]);
  }
  indedgelocnbr = indedgelocnum - baseval;
  indgrafptr->vertloctax[indvertlocnum] = indedgelocnum; /* Mark end of edge array */
  indgrafptr->vendloctax = indgrafptr->vertloctax + 1; /* Induced graph is compact */
  indgrafptr->velolocsum = indvelolocsum;
  indgrafptr->edgelocnbr = indedgelocnbr;
  indgrafptr->edgelocsiz = indedgelocnbr;
  if (orgvlblloctax == NULL)                      /* If we didn't have vertex labels, use vertex index array as vertex label array */
    indgrafptr->vlblloctax = indgrafptr->vnumloctax;

  if (indedloloctax != NULL) {                    /* Re-allocate arrays and delete orgindxtab             */
    size_t              indedlooftval;            /* Offset of edge load array with respect to edge array */

    indedlooftval = indedloloctax - indgrafptr->edgeloctax;
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + baseval, (indedlooftval + indedgelocnbr) * sizeof (Gnum));
    indgrafptr->edgeloctax -= baseval;
    indedloloctax = indgrafptr->edgeloctax + indedlooftval; /* Use old index into old array as new index to avoid alignment problems */
  }
  else {
    indgrafptr->edgeloctax  = memRealloc (indgrafptr->edgeloctax + baseval, indedgelocnbr * sizeof (Gnum));
    indgrafptr->edgeloctax -= baseval;
  }

  indvertlocnum = baseval;
  indvnumloctax = indgrafptr->vnumloctax;         /* TRICK: maybe alias */
  if (orgvnumloctax != NULL) {                    /* Adjust vnumloctax  */
    for ( ; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indvnumloctax[indvertlocnum] = orgvnumloctax[indvnumloctmp[indvertlocnum]]; /* TRICK: indvnumloctmp and indgrafptr->vnumloctax may be aliases */
  }
  else {
    Gnum                orgvertglbadj;

    orgvertglbadj = orggrafptr->procvrttab[orggrafptr->proclocnum] - baseval; /* Set adjustement for global indexing */
    for ( ; indvertlocnum < indvertlocnnd; indvertlocnum ++)
      indvnumloctax[indvertlocnum] = indvnumloctmp[indvertlocnum] + orgvertglbadj; /* TRICK: indvnumloctmp and indgrafptr->vnumloctax may be aliases */
  }

  indgrafptr->edloloctax = indedloloctax;
  indgrafptr->degrglbmax = inddegrlocmax;         /* Local maximum degree will be turned into global maximum degree */
  if (dgraphBuild4 (indgrafptr) != 0) {
    errorPrint ("dgraphInduce2: cannot build induced graph");
    return     (1);
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (indgrafptr) != 0) {            /* Check graph consistency */
    errorPrint ("dgraphInduce2: inconsistent graph data");
    dgraphFree (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}

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

static
Gnum
dgraphInduceList2 (
Dgraph * restrict const     indgrafptr,
Dgraph * restrict const     orggrafptr,
const void * restrict const orgdataptr,           /* Data pointer is list array pointer */
Gnum * restrict const       orgindxgsttax)
{
  Gnum                orglistnbr;
  Gnum                orglistnum;
  Gnum                indvertglbnum;
  Gnum                indedgelocmax;

  const Gnum * restrict const orgvertloctax = orggrafptr->vertloctax;
  const Gnum * restrict const orgvendloctax = orggrafptr->vendloctax;
  const Gnum * restrict const orglisttab    = (const Gnum * restrict) orgdataptr;

  for (orglistnum = 0, indvertglbnum = indgrafptr->procvrttab[indgrafptr->proclocnum], /* Fill index array while recomputing tighter upper bound on arcs */
       orglistnbr = indgrafptr->vertlocnbr, indedgelocmax = 0;
       orglistnum < orglistnbr; orglistnum ++, indvertglbnum ++) {
    Gnum                orgvertlocnum;

    orgvertlocnum = orglisttab[orglistnum];
    orgindxgsttax[orgvertlocnum] = indvertglbnum; /* Mark selected vertices */
    indedgelocmax += orgvendloctax[orgvertlocnum] - orgvertloctax[orgvertlocnum];
  }

  return (indedgelocmax);
}

int
dgraphInduceList (
Dgraph * restrict const       orggrafptr,
const Gnum                    orglistnbr,
const Gnum * const            orglisttab,         /* Local list of kept vertices */
Dgraph * restrict const       indgrafptr)
{
  return (dgraphInduce2 (orggrafptr, dgraphInduceList2, (const void * const) orglisttab, orglistnbr, (Gnum * const) (orglisttab - orggrafptr->baseval), indgrafptr));
}

/* This routine builds the graph induced
** by the original graph and the vector of
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

typedef struct DgraphInducePartData_ {
  const GraphPart *         orgpartloctax;
  GraphPart                 indpartval;
} DgraphInducePartData;

static
Gnum
dgraphInducePart2 (
Dgraph * restrict const     indgrafptr,
Dgraph * restrict const     orggrafptr,
const void * restrict const orgdataptr,
Gnum * restrict const       orgindxgsttax)
{
  Gnum                orgvertlocnnd;
  Gnum                orgvertlocnum;
  Gnum                indvertlocnum;
  Gnum                indvertglbnum;
  Gnum                indedgelocmax;

  const Gnum * restrict const       orgvertloctax = orggrafptr->vertloctax;
  const Gnum * restrict const       orgvendloctax = orggrafptr->vendloctax;
  const GraphPart * restrict const  orgpartloctax = ((const DgraphInducePartData * restrict const) orgdataptr)->orgpartloctax;
  const GraphPart                   indpartval    = ((const DgraphInducePartData * restrict const) orgdataptr)->indpartval;
  Gnum * restrict const             indvnumloctax = indgrafptr->vnumloctax;

  for (orgvertlocnum = indvertlocnum = orggrafptr->baseval, indvertglbnum = indgrafptr->procvrttab[indgrafptr->proclocnum], /* Fill index array while recomputing tighter upper bound on arcs */
       orgvertlocnnd = orggrafptr->vertlocnnd, indedgelocmax = 0;
       orgvertlocnum < orgvertlocnnd; orgvertlocnum ++) {
    if (orgpartloctax[orgvertlocnum] == indpartval) {
      orgindxgsttax[orgvertlocnum] = indvertglbnum; /* Mark selected vertices */
      indvnumloctax[indvertlocnum] = orgvertlocnum;
      indedgelocmax += orgvendloctax[orgvertlocnum] - orgvertloctax[orgvertlocnum];
      indvertlocnum ++, indvertglbnum ++;
    }
    else
      orgindxgsttax[orgvertlocnum] = ~0;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((indvertlocnum - orggrafptr->baseval) != indgrafptr->vertlocnbr) {
    errorPrint ("dgraphInducePart2: inconsistent data");
    dgraphFree (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (indedgelocmax);
}

int
dgraphInducePart (
Dgraph * restrict const           orggrafptr,     /* Pointer to original distributed graph       */
const GraphPart * restrict const  orgpartloctax,  /* Based array of local vertex partition flags */  
const Gnum                        indvertlocnbr,  /* Number of local vertices in selected part   */
const GraphPart                   indpartval,
Dgraph * restrict const           indgrafptr)
{
  DgraphInducePartData  orgdatadat;

  orgdatadat.orgpartloctax = orgpartloctax;
  orgdatadat.indpartval    = indpartval;

  return (dgraphInduce2 (orggrafptr, dgraphInducePart2, &orgdatadat, indvertlocnbr, NULL, indgrafptr));
}
