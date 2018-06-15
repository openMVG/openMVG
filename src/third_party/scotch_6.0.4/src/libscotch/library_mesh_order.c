/* Copyright 2004,2007,2008,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_mesh_order.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the mesh     **/
/**                ordering routines of the libSCOTCH      **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 aug 1998     **/
/**                                 to     22 aug 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     27 mar 1999     **/
/**                # Version 4.0  : from : 29 jan 2002     **/
/**                                 to     20 dec 2005     **/
/**                # Version 5.0  : from : 04 aug 2007     **/
/**                                 to     31 may 2008     **/
/**                # Version 5.1  : from : 29 mar 2010     **/
/**                                 to     14 aug 2010     **/
/**                # Version 6.0  : from : 14 nov 2012     **/
/**                                 to     14 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "mesh.h"
#include "order.h"
#include "hmesh.h"
#include "hmesh_order_st.h"
#include "library_order.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mesh ordering routines.      */
/*                                  */
/************************************/

/*+ This routine initializes an API ordering
*** with respect to the given source graph
*** and the locations of output parameters.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderInit (
const SCOTCH_Mesh * const   meshptr,              /*+ Mesh to order                      +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering structure to initialize   +*/
SCOTCH_Num * const          permtab,              /*+ Direct permutation array           +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Separator tree array               +*/
{
  Mesh *              srcmeshptr;
  LibOrder *          libordeptr;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Ordering) < sizeof (LibOrder)) {
    errorPrint ("SCOTCH_meshOrderInit: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  srcmeshptr = (Mesh *) meshptr;                  /* Use structure as source mesh */
  libordeptr = (LibOrder *) ordeptr;
  libordeptr->permtab = ((permtab == NULL) || ((void *) permtab == (void *) meshptr)) ? NULL : (Gnum *) permtab;
  libordeptr->peritab = ((peritab == NULL) || ((void *) peritab == (void *) meshptr)) ? NULL : (Gnum *) peritab;
  libordeptr->cblkptr = ((cblkptr == NULL) || ((void *) cblkptr == (void *) meshptr)) ? NULL : (Gnum *) cblkptr;
  libordeptr->rangtab = ((rangtab == NULL) || ((void *) rangtab == (void *) meshptr)) ? NULL : (Gnum *) rangtab;
  libordeptr->treetab = ((treetab == NULL) || ((void *) treetab == (void *) meshptr)) ? NULL : (Gnum *) treetab;

  return (orderInit (&libordeptr->o, srcmeshptr->baseval, srcmeshptr->vnodnbr, libordeptr->peritab));
}

/*+ This routine frees an API ordering.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_meshOrderExit (
const SCOTCH_Mesh * const   meshptr,
SCOTCH_Ordering * const     ordeptr)
{
  orderExit (&((LibOrder *) ordeptr)->o);
}

/*+ This routine saves the contents of
*** the given ordering to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderSave (
const SCOTCH_Mesh * const     meshptr,            /*+ Mesh to order    +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSave (&((LibOrder *) ordeptr)->o, ((Mesh *) meshptr)->vlbltax, stream));
}

/*+ This routine saves the mapping data
*** associated with the given ordering
*** to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderSaveMap (
const SCOTCH_Mesh * const     meshptr,            /*+ Mesh to order    +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSaveMap (&((LibOrder *) ordeptr)->o, ((Mesh *) meshptr)->vlbltax, stream));
}

/*+ This routine saves to the given stream
*** the separator tree data associated with
*** the given ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderSaveTree (
const SCOTCH_Mesh * const     meshptr,            /*+ Mesh to order    +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSaveTree (&((LibOrder *) ordeptr)->o, ((Mesh *) meshptr)->vlbltax, stream));
}

/*+ This routine computes an ordering
*** of the API ordering structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderCompute (
SCOTCH_Mesh * const         meshptr,              /*+ Mesh to order       +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering to compute +*/
SCOTCH_Strat * const        stratptr)             /*+ Ordering strategy   +*/
{
  return (SCOTCH_meshOrderComputeList (meshptr, ordeptr, 0, NULL, stratptr));
}

/*+ This routine computes a partial ordering
*** of the listed nodes of the API ordering
*** structure mesh with respect to the given
*** strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderComputeList (
SCOTCH_Mesh * const         meshptr,              /*+ Mesh to order                   +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering to compute             +*/
const SCOTCH_Num            listnbr,              /*+ Number of vertices in list      +*/
const SCOTCH_Num * const    listtab,              /*+ List of vertex indices to order +*/
SCOTCH_Strat * const        stratptr)             /*+ Ordering strategy               +*/
{
  LibOrder *          libordeptr;                 /* Pointer to ordering             */
  Mesh *              srcmeshptr;                 /* Pointer to source mesh          */
  Hmesh               srcmeshdat;                 /* Halo source mesh structure      */
  VertList            srclistdat;                 /* Subgraph vertex list            */
  VertList *          srclistptr;                 /* Pointer to subgraph vertex list */
  const Strat *       ordstratptr;                /* Pointer to ordering strategy    */

  srcmeshptr = (Mesh *) meshptr;

#ifdef SCOTCH_DEBUG_MESH2
  if (meshCheck (srcmeshptr) != 0) {
    errorPrint ("SCOTCH_meshOrderComputeList: invalid input mesh");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  if (*((Strat **) stratptr) == NULL)             /* Set default ordering strategy if necessary */
    SCOTCH_stratMeshOrderBuild (stratptr, SCOTCH_STRATQUALITY, 0.1);

  ordstratptr = *((Strat **) stratptr);
  if (ordstratptr->tabl != &hmeshorderststratab) {
    errorPrint ("SCOTCH_meshOrderComputeList: not a mesh ordering strategy");
    return     (1);
  }

  memCpy (&srcmeshdat.m, srcmeshptr, sizeof (Mesh)); /* Copy non-halo mesh data  */
  srcmeshdat.m.flagval &= ~MESHFREETABS;          /* Do not allow to free arrays */
  srcmeshdat.vehdtax    = srcmeshdat.m.vendtax;   /* End of non-halo vertices    */
  srcmeshdat.veihnbr    = 0;                      /* No halo isolated elements   */
  srcmeshdat.vnohnbr    = srcmeshdat.m.vnodnbr;   /* All nodes are non-halo      */
  srcmeshdat.vnohnnd    = srcmeshdat.m.vnodnnd;   /* No halo present             */
  srcmeshdat.vnhlsum    = srcmeshdat.m.vnlosum;   /* Sum of node vertex weights  */
  srcmeshdat.enohnbr    = srcmeshdat.m.edgenbr;   /* All edges are non-halo      */
  srcmeshdat.levlnum    = 0;                      /* Start from level zero       */

  libordeptr         = (LibOrder *) ordeptr;      /* Get ordering      */
  srclistdat.vnumnbr = (Gnum)   listnbr;          /* Build vertex list */
  srclistdat.vnumtab = (Gnum *) listtab;
  srclistptr = ((srclistdat.vnumnbr == 0) ||
                (srclistdat.vnumnbr == srcmeshdat.m.vnodnbr))
               ? NULL : &srclistdat;              /* Is the list really necessary */
  if (srclistptr != NULL) {
    errorPrint ("SCOTCH_meshOrderComputeList: node lists not yet implemented");
    return     (1);
  }

  intRandInit ();                                 /* Check that random number generator is initialized */

  hmeshOrderSt (&srcmeshdat, &libordeptr->o, 0, &libordeptr->o.cblktre, ordstratptr);

#ifdef SCOTCH_DEBUG_LIBRARY2
  orderCheck (&libordeptr->o);
#endif /* SCOTCH_DEBUG_LIBRARY2 */

  if (libordeptr->permtab != NULL)                 /* Build direct permutation if wanted */
    orderPeri (libordeptr->o.peritab, libordeptr->o.baseval, libordeptr->o.vnodnbr, libordeptr->permtab, libordeptr->o.baseval);
  if (libordeptr->rangtab != NULL)                /* Build range array if column block data wanted */
    orderRang (&libordeptr->o, libordeptr->rangtab);
  if (libordeptr->treetab != NULL)                /* Build separator tree array if wanted */
      orderTree (&libordeptr->o, libordeptr->treetab);
  if (libordeptr->cblkptr != NULL)                /* Set number of column blocks if wanted */
    *(libordeptr->cblkptr) = libordeptr->o.cblknbr;

  meshExit (&srcmeshdat.m);                       /* Free in case mesh had been reordered */

  return (0);
}

/*+ This routine computes an ordering
*** of the API ordering structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrder (
SCOTCH_Mesh * const         meshptr,              /*+ Mesh to order                      +*/
SCOTCH_Strat * const        stratptr,             /*+ Ordering strategy                  +*/
SCOTCH_Num * const          permtab,              /*+ Ordering permutation               +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Separator tree array               +*/
{
  SCOTCH_Ordering     ordedat;
  int                 o;

  SCOTCH_meshOrderInit (meshptr, &ordedat, permtab, peritab, cblkptr, rangtab, treetab);
  o = SCOTCH_meshOrderCompute (meshptr, &ordedat, stratptr);
  SCOTCH_meshOrderExit (meshptr, &ordedat);

  return (o);
}

/*+ This routine computes an ordering
*** of the submesh of the API ordering
*** structure mesh induced by the given
*** vertex list, with respect to the given
*** strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderList (
SCOTCH_Mesh * const         meshptr,              /*+ Mesh to order                      +*/
const SCOTCH_Num            listnbr,              /*+ Number of vertices in list         +*/
const SCOTCH_Num * const    listtab,              /*+ List of vertex indices to order    +*/
SCOTCH_Strat * const        stratptr,             /*+ Ordering strategy                  +*/
SCOTCH_Num * const          permtab,              /*+ Ordering permutation               +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Column block range array           +*/
{
  SCOTCH_Ordering     ordedat;
  int                 o;

  SCOTCH_meshOrderInit (meshptr, &ordedat, permtab, peritab, cblkptr, rangtab, treetab);
  o = SCOTCH_meshOrderComputeList (meshptr, &ordedat, listnbr, listtab, stratptr);
  SCOTCH_meshOrderExit (meshptr, &ordedat);

  return (o);
}

/*+ This routine checks the consistency
*** of the given mesh ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshOrderCheck (
const SCOTCH_Mesh * const     meshptr,
const SCOTCH_Ordering * const ordeptr)            /*+ Ordering to check +*/
{
  return (orderCheck (&((LibOrder *) ordeptr)->o));
}

/*+ This routine parses the given
*** mesh ordering strategy.
*** It returns:
*** - 0   : if string successfully scanned.
*** - !0  : on error.
+*/

int
SCOTCH_stratMeshOrder (
SCOTCH_Strat * const        stratptr,
const char * const          string)
{
  if (*((Strat **) stratptr) != NULL)
    stratExit (*((Strat **) stratptr));

  if ((*((Strat **) stratptr) = stratInit (&hmeshorderststratab, string)) == NULL) {
    errorPrint ("SCOTCH_stratMeshOrder: error in ordering strategy");
    return     (1);
  }

  return (0);
}

/*+ This routine provides predefined
*** ordering strategies.
*** It returns:
*** - 0   : if string successfully initialized.
*** - !0  : on error.
+*/

int
SCOTCH_stratMeshOrderBuild (
SCOTCH_Strat * const        stratptr,             /*+ Strategy to create      +*/
const SCOTCH_Num            flagval,              /*+ Desired characteristics +*/
const double                balrat)               /*+ Desired imbalance ratio +*/
{
  char                bufftab[8192];              /* Should be enough */
  char                bbaltab[32];

  strcpy (bufftab, "c{rat=0.7,cpr=n{sep=/(vnod>120)?m{vnod=100,low=h{pass=10},asc=f{bal=<BBAL>}}:;,ole=v{strat=d{cmin=0,cmax=10000000,frat=0}},ose=g},unc=n{sep=/(vnod>120)?m{vnod=100,low=h{pass=10},asc=f{bal=<BBAL>}}:;,ole=v{strat=d{cmin=0,cmax=10000000,frat=0}},ose=g}}");

  sprintf (bbaltab, "%lf", balrat);
  stringSubst (bufftab, "<BBAL>", bbaltab);

  if (SCOTCH_stratMeshOrder (stratptr, bufftab) != 0) {
    errorPrint ("SCOTCH_stratMeshOrderBuild: error in sequential ordering strategy");
    return     (1);
  }

  return (0);
}
