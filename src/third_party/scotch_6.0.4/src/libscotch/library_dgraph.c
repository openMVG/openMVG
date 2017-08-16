/* Copyright 2007,2009,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source graph handling routines of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 26 apr 2006     **/
/**                                 to     14 apr 2008     **/
/**                # Version 5.1  : from : 26 mar 2009     **/
/**                                 to     17 nov 2010     **/
/**                # Version 6.0  : from : 27 nov 2012     **/
/**                                 to     29 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"                                /* For graphPtscotch() */
#include "dgraph.h"
#include "ptscotch.h"

/****************************************/
/*                                      */
/* These routines are the C API for the */
/* distributed graph handling routines. */
/*                                      */
/****************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** distributed graph structure.
*** It returns:
*** - !NULL  : if the initialization succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Dgraph *
SCOTCH_dgraphAlloc ()
{
  return ((SCOTCH_Dgraph *) memAlloc (sizeof (SCOTCH_Dgraph)));
}

/*+ This routine initializes the opaque
*** distributed graph structure used to
*** handle distributed graphs in the
*** Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphInit (
SCOTCH_Dgraph * const       grafptr,
MPI_Comm                    proccomm)             /* Communicator to be used for all communications */
{
#ifdef SCOTCH_PTHREAD
  int                 thrdlvlval;
#endif /* SCOTCH_PTHREAD */

#ifdef SCOTCH_PTHREAD
  MPI_Query_thread (&thrdlvlval);
  if (thrdlvlval < MPI_THREAD_MULTIPLE) {
    errorPrint ("SCOTCH_dgraphInit: Scotch compiled with SCOTCH_PTHREAD and program not launched with MPI_THREAD_MULTIPLE");
    return     (1);
  }
#endif /* SCOTCH_PTHREAD */

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint ("SCOTCH_dgraphInit: internal error (1)");
    return     (1);
  }
  if (sizeof (SCOTCH_Dgraph) < sizeof (Dgraph)) {
    errorPrint ("SCOTCH_dgraphInit: internal error (2)");
    return     (1);
  }

  return (dgraphInit ((Dgraph *) grafptr, proccomm));
}

/*+ This routine frees the contents of the
*** given opaque graph structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dgraphExit (
SCOTCH_Dgraph * const       grafptr)
{
  dgraphExit ((Dgraph *) grafptr);
}

/*+ This routine frees the contents of the
*** given opaque graph structure but does
*** not free its private data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dgraphFree (
SCOTCH_Dgraph * const       grafptr)
{
  dgraphFree ((Dgraph *) grafptr);
}

/*+ This routine accesses graph size data.
*** NULL pointers on input indicate unwanted
*** data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dgraphSize (
const SCOTCH_Dgraph * const grafptr,
SCOTCH_Num * const          vertglbnbr,
SCOTCH_Num * const          vertlocnbr,
SCOTCH_Num * const          edgeglbnbr,
SCOTCH_Num * const          edgelocnbr)
{
  const Dgraph *      srcgrafptr;

  srcgrafptr = (Dgraph *) grafptr;

  if (vertglbnbr != NULL)
    *vertglbnbr = (SCOTCH_Num) (srcgrafptr->vertglbnbr);
  if (vertlocnbr != NULL)
    *vertlocnbr = (SCOTCH_Num) (srcgrafptr->vertlocnbr);
  if (edgeglbnbr != NULL)
    *edgeglbnbr = (SCOTCH_Num) srcgrafptr->edgeglbnbr;
  if (edgelocnbr != NULL)
    *edgelocnbr = (SCOTCH_Num) srcgrafptr->edgelocnbr;
}

/*+ This routine accesses all of the graph data.
*** NULL pointers on input indicate unwanted
*** data. NULL pointers on output indicate
*** unexisting arrays.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dgraphData (
const SCOTCH_Dgraph * const grafptr,              /* Graph structure to read          */
SCOTCH_Num * const          baseptr,              /* Base value                       */
SCOTCH_Num * const          vertglbptr,           /* Number of global vertices        */
SCOTCH_Num * const          vertlocptr,           /* Number of local vertices         */
SCOTCH_Num * const          vertlocptz,           /* Maximum number of local vertices */
SCOTCH_Num * const          vertgstptr,           /* Number of local + ghost vertices */
SCOTCH_Num ** const         vertloctab,           /* Vertex array [vertnbr+1]         */
SCOTCH_Num ** const         vendloctab,           /* Vertex array [vertnbr]           */
SCOTCH_Num ** const         veloloctab,           /* Vertex load array                */
SCOTCH_Num ** const         vlblloctab,           /* Vertex label array               */
SCOTCH_Num * const          edgeglbptr,           /* Number of global edges (arcs)    */
SCOTCH_Num * const          edgelocptr,           /* Number of local edges (arcs)     */
SCOTCH_Num * const          edgelocptz,           /* Size of local edge array         */
SCOTCH_Num ** const         edgeloctab,           /* Local edge array [edgelocsiz]    */
SCOTCH_Num ** const         edgegsttab,           /* Ghost edge array [edgelocsiz]    */
SCOTCH_Num ** const         edloloctab,           /* Edge load array [edgelocsiz]     */
MPI_Comm * const            comm)                 /* MPI Communicator                 */
{
  const Dgraph *      srcgrafptr;                 /* Pointer to source graph structure */

  srcgrafptr = (const Dgraph *) grafptr;

  if (baseptr != NULL)
    *baseptr = srcgrafptr->baseval;
  if (vertglbptr != NULL)
    *vertglbptr = srcgrafptr->vertglbnbr;
  if (vertlocptr != NULL)
    *vertlocptr = srcgrafptr->vertlocnbr;
  if (vertlocptz != NULL)
    *vertlocptz = srcgrafptr->procvrttab[srcgrafptr->proclocnum + 1] - srcgrafptr->procvrttab[srcgrafptr->proclocnum];
  if (vertgstptr != NULL)
    *vertgstptr = ((srcgrafptr->flagval & DGRAPHHASEDGEGST) != 0) ? srcgrafptr->vertgstnbr : -1;
  if (vertloctab != NULL)
    *vertloctab = srcgrafptr->vertloctax + srcgrafptr->baseval;
  if (vendloctab != NULL)
    *vendloctab = srcgrafptr->vendloctax + srcgrafptr->baseval;
  if (veloloctab != NULL)
    *veloloctab = (srcgrafptr->veloloctax != NULL) ? srcgrafptr->veloloctax + srcgrafptr->baseval : NULL;
  if (vlblloctab != NULL)
    *vlblloctab = (srcgrafptr->vlblloctax != NULL) ? srcgrafptr->vlblloctax + srcgrafptr->baseval : NULL;
  if (edgeglbptr != NULL)
    *edgeglbptr = srcgrafptr->edgeglbnbr;
  if (edgelocptr != NULL)
    *edgelocptr = srcgrafptr->edgelocnbr;
  if (edgelocptz != NULL)
    *edgelocptz = srcgrafptr->edgelocsiz;
  if (edgeloctab != NULL)
    *edgeloctab = srcgrafptr->edgeloctax + srcgrafptr->baseval;
  if (edgegsttab != NULL)
    *edgegsttab = (srcgrafptr->edgegsttax != NULL) ? srcgrafptr->edgegsttax + srcgrafptr->baseval : NULL;
  if (edloloctab != NULL)
    *edloloctab = (srcgrafptr->edloloctax != NULL) ? srcgrafptr->edloloctax + srcgrafptr->baseval : NULL;
  if (comm != NULL)
    *comm = srcgrafptr->proccomm;
}
