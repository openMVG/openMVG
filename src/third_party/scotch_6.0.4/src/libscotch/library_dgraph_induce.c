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
/**   NAME       : library_dgraph_induce.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the          **/
/**                distributed graph induction routine of  **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 30 aug 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/**   NOTES      : # This code is directly derived from    **/
/**                  the code of dgraphInducePart() and    **/
/**                  of its subroutines. The only change   **/
/**                  is that it uses Gnum's instead of     **/
/**                  GraphPart's as part values.           **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "ptscotch.h"

/***********************************/
/*                                 */
/* This routine is the C API for   */
/* the distributed graph induction */
/* routine.                        */
/*                                 */
/***********************************/

/*+ This routine creates a distributed induced graph
*** from the given graph, according to the partition
*** map that is passed to the routine.
*** It returns:
*** - 0   : if the induced graph has been created.
*** - !0  : on error.
+*/

typedef struct _SCOTCHDgraphInducePartData_ {
  const Gnum *              orgpartloctax;        /* In the public interface, parts are represented as Gnum's */
  Gnum                      indpartval;
} _SCOTCHDgraphInducePartData;

static
Gnum
_SCOTCHdgraphInducePart2 (
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

  const Gnum * restrict const orgvertloctax = orggrafptr->vertloctax;
  const Gnum * restrict const orgvendloctax = orggrafptr->vendloctax;
  const Gnum * restrict const orgpartloctax = ((const _SCOTCHDgraphInducePartData * restrict const) orgdataptr)->orgpartloctax;
  const Gnum                  indpartval    = ((const _SCOTCHDgraphInducePartData * restrict const) orgdataptr)->indpartval;
  Gnum * restrict const       indvnumloctax = indgrafptr->vnumloctax;

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
    dgraphExit (indgrafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (indedgelocmax);
}

int
SCOTCH_dgraphInducePart (
SCOTCH_Dgraph * const       orggrafptr,           /* Original graph                   */
const SCOTCH_Num * const    orgpartloctab,        /* Partition array                  */
const SCOTCH_Num            indpartval,           /* Part value of induced subgraph   */
const SCOTCH_Num            indvertlocnbr,        /* Number of local vertices in part */
SCOTCH_Dgraph * const       indgrafptr)           /* Induced subgraph                 */
{
  _SCOTCHDgraphInducePartData orgdatadat;
  Gnum                        indvertloctmp;
  int                         o;

#ifdef SCOTCH_DEBUG_LIBRARY1
  MPI_Comm_compare (((Dgraph * restrict const) orggrafptr)->proccomm,
                    ((Dgraph * restrict const) indgrafptr)->proccomm, &o);
  if ((o != MPI_IDENT) && (o != MPI_CONGRUENT)) {
    errorPrint ("SCOTCH_dgraphInducePart: communicators are not congruent");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  if (indvertlocnbr < 0) {                        /* If number of kept vertices is not known, compute it */
    Gnum                orgvertlocnum;
    Gnum                orgvertlocnbr;

    for (orgvertlocnum = indvertloctmp = 0, orgvertlocnbr = ((Dgraph * restrict const) orggrafptr)->vertlocnbr;
         orgvertlocnum < orgvertlocnbr; orgvertlocnum ++) {
      if (orgpartloctab[orgvertlocnum] == indpartval)
        indvertloctmp ++;
    }
  }
  else
    indvertloctmp = indvertlocnbr;

  orgdatadat.orgpartloctax = orgpartloctab - ((Dgraph *) orggrafptr)->baseval;
  orgdatadat.indpartval    = indpartval;

  o = dgraphInduce2 ((Dgraph *) orggrafptr, _SCOTCHdgraphInducePart2, &orgdatadat, indvertloctmp, NULL, (Dgraph *) indgrafptr);
  ((Dgraph *) indgrafptr)->vnumloctax = NULL;     /* Do not impact subsequent inductions */
  return (o);
}
