/* Copyright 2007,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_build.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source graph handling routines of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 feb 2007     **/
/**                                 to     18 jul 2007     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
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
#include "ptscotch.h"

/****************************************/
/*                                      */
/* These routines are the C API for the */
/* distributed graph handling routines. */
/*                                      */
/****************************************/

/*+ This routine fills the contents of the given
*** opaque distributed graph structure with the
*** data provided by the user. The base value
*** allows the user to set the graph base to 0 or 1.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphBuild (
SCOTCH_Dgraph * const       grafptr,              /* Distributed graph structure to fill  */
const Gnum                  baseval,              /* Base for indexing                    */
const Gnum                  vertlocnbr,           /* Number of local vertices             */
const Gnum                  vertlocmax,           /* Maximum number of local vertices     */
Gnum * const                vertloctab,           /* Local vertex begin array             */
Gnum * const                vendloctab,           /* Local vertex end array               */
Gnum * const                veloloctab,           /* Local vertex load array (if any)     */
Gnum * const                vlblloctab,           /* Local vertex label array (if any)    */
const Gnum                  edgelocnbr,           /* Number of local edges                */
const Gnum                  edgelocsiz,           /* Size of local edge array             */
Gnum * const                edgeloctab,           /* Local edge array                     */
Gnum * const                edgegsttab,           /* Ghost edge array (if any); not const */
Gnum * const                edloloctab)           /* Local edge load array (if any)       */
{
  Dgraph *                    srcgrafptr;         /* Pointer to source graph structure */
  Gnum *                      vertloctax;
  Gnum *                      vendloctax;
  Gnum *                      veloloctax;
  Gnum *                      vlblloctax;
  Gnum *                      edgeloctax;
  Gnum *                      edgegsttax;
  Gnum *                      edloloctax;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Dgraph) < sizeof (Dgraph)) {
    errorPrint ("SCOTCH_dgraphBuild: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
  if ((baseval < 0) || (baseval > 1)) {
    errorPrint ("SCOTCH_dgraphBuild: invalid base parameter");
    return     (1);
  }

  srcgrafptr = (Dgraph *) grafptr;                /* Use structure as source graph */
  
  vertloctax = (Gnum *) vertloctab - baseval;
  vendloctax = ((vendloctab == NULL) || (vendloctab == vertloctab + 1)) ? vertloctax + 1 : (Gnum *) vendloctab - baseval;
  veloloctax = ((veloloctab == NULL) || (veloloctab == vertloctab)) ? NULL : (Gnum *) veloloctab - baseval;
  vlblloctax = ((vlblloctab == NULL) || (vlblloctab == vertloctab)) ? NULL : (Gnum *) vlblloctab - baseval;
  edgeloctax = (Gnum *) edgeloctab - baseval;
  edgegsttax = ((edgegsttab == NULL) || (edgegsttab == edgeloctab)) ? NULL : (Gnum *) edgegsttab - baseval;
  edloloctax = ((edloloctab == NULL) || (edloloctab == edgeloctab)) ? NULL : (Gnum *) edloloctab - baseval;

  return (dgraphBuild (srcgrafptr, baseval,
                       vertlocnbr, vertlocmax, vertloctax, vendloctax, veloloctax, NULL, vlblloctax,
                       edgelocnbr, edgelocsiz, edgeloctax, edgegsttax, edloloctax));
}
