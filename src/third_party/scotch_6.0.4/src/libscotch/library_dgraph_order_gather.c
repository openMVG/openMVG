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
/**   NAME       : library_dgraph_order_gather.c           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted ordering gathering routines of    **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 21 jul 2007     **/
/**                                 to     04 aug 2007     **/
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
#include "order.h"
#include "dorder.h"
#include "library_order.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the distributed ordering         */
/* handling routines.               */
/*                                  */
/************************************/

/*+ This routine initializes an API centralized
*** ordering with respect to the given distributed
*** source graph and the locations of output parameters.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphCorderInit (
const SCOTCH_Dgraph * const grafptr,              /*+ Distributed graph to order         +*/
SCOTCH_Ordering * const     cordptr,              /*+ Ordering structure to initialize   +*/
SCOTCH_Num * const          permtab,              /*+ Direct permutation array           +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Separator tree array               +*/
{
  Dgraph *            srcgrafptr;
  LibOrder *          libcordptr;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Ordering) < sizeof (LibOrder)) {
    errorPrint ("SCOTCH_dgraphCorderInit: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  srcgrafptr = (Dgraph *) grafptr;                /* Use structure as distributed source graph */
  libcordptr = (LibOrder *) cordptr;
  libcordptr->permtab = ((permtab == NULL) || ((void *) permtab == (void *) grafptr)) ? NULL : (Gnum *) permtab;
  libcordptr->peritab = ((peritab == NULL) || ((void *) peritab == (void *) grafptr)) ? NULL : (Gnum *) peritab;
  libcordptr->cblkptr = ((cblkptr == NULL) || ((void *) cblkptr == (void *) grafptr)) ? NULL : (Gnum *) cblkptr;
  libcordptr->rangtab = ((rangtab == NULL) || ((void *) rangtab == (void *) grafptr)) ? NULL : (Gnum *) rangtab;
  libcordptr->treetab = ((treetab == NULL) || ((void *) treetab == (void *) grafptr)) ? NULL : (Gnum *) treetab;

  return (orderInit (&libcordptr->o, srcgrafptr->baseval, srcgrafptr->vertglbnbr, libcordptr->peritab));
}

/*+ This routine frees an API centralized ordering.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dgraphCorderExit (
const SCOTCH_Dgraph * const grafptr,
SCOTCH_Ordering * const     cordptr)
{
  orderExit (&((LibOrder *) cordptr)->o);
}

/*+ This routine gathers the contents of
*** the given distributed ordering into the
*** given centralized ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphOrderGather (
const SCOTCH_Dgraph * const     grafptr,          /*+ Not used             +*/
const SCOTCH_Dordering * const  dordptr,          /*+ Distributed ordering +*/
SCOTCH_Ordering * const         cordptr)          /*+ Centralized ordering +*/
{
  LibOrder *          libcordptr;                 /* Pointer to ordering */

  if ((cordptr != NULL) && ((void *) cordptr != (void *) dordptr)) { /* If potential root process */
    libcordptr = (LibOrder *) cordptr;            /* Get centralized ordering                     */

    if (dorderGather ((Dorder *) dordptr, &libcordptr->o) != 0)
      return (1);

    if (libcordptr->permtab != NULL)              /* Build direct permutation if wanted */
      orderPeri (libcordptr->o.peritab, libcordptr->o.baseval, libcordptr->o.vnodnbr, libcordptr->permtab, libcordptr->o.baseval);
    if (libcordptr->rangtab != NULL)              /* Build range array if column block data wanted */
      orderRang (&libcordptr->o, libcordptr->rangtab);
    if (libcordptr->treetab != NULL)              /* Build separator tree array if wanted */
      orderTree (&libcordptr->o, libcordptr->treetab);
    if (libcordptr->cblkptr != NULL)              /* Set number of column blocks if wanted */
      *(libcordptr->cblkptr) = libcordptr->o.cblknbr;

    return (0);
  }
  else
    return (dorderGather ((Dorder *) dordptr, NULL));
}
