/* Copyright 2004,2007,2008,2010,2012-2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_order.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the graph    **/
/**                ordering routines of the libSCOTCH      **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 aug 1998     **/
/**                                 to     22 aug 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     27 mar 1999     **/
/**                # Version 4.0  : from : 29 jan 2002     **/
/**                                 to     08 sep 2006     **/
/**                # Version 5.0  : from : 19 dec 2006     **/
/**                                 to     04 aug 2007     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to     14 aug 2010     **/
/**                # Version 6.0  : from : 08 jan 2012     **/
/**                                 to     15 nov 2014     **/
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
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_st.h"
#include "library_order.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph ordering routines.     */
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
SCOTCH_graphOrderInit (
const SCOTCH_Graph * const  grafptr,              /*+ Graph to order                     +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering structure to initialize   +*/
SCOTCH_Num * const          permtab,              /*+ Direct permutation array           +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Separator tree array               +*/
{
  const Graph *       srcgrafptr;
  LibOrder *          libordeptr;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Ordering) < sizeof (LibOrder)) {
    errorPrint ("SCOTCH_graphOrderInit: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  srcgrafptr = (Graph *) grafptr;                 /* Use structure as source graph */
  libordeptr = (LibOrder *) ordeptr;
  libordeptr->permtab = ((permtab == NULL) || ((void *) permtab == (void *) grafptr)) ? NULL : (Gnum *) permtab;
  libordeptr->peritab = ((peritab == NULL) || ((void *) peritab == (void *) grafptr)) ? NULL : (Gnum *) peritab;
  libordeptr->cblkptr = ((cblkptr == NULL) || ((void *) cblkptr == (void *) grafptr)) ? NULL : (Gnum *) cblkptr;
  libordeptr->rangtab = ((rangtab == NULL) || ((void *) rangtab == (void *) grafptr)) ? NULL : (Gnum *) rangtab;
  libordeptr->treetab = ((treetab == NULL) || ((void *) treetab == (void *) grafptr)) ? NULL : (Gnum *) treetab;

  return (orderInit (&libordeptr->o, srcgrafptr->baseval, srcgrafptr->vertnbr, libordeptr->peritab));
}

/*+ This routine frees an API ordering.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphOrderExit (
const SCOTCH_Graph * const  grafptr,
SCOTCH_Ordering * const     ordeptr)
{
  orderExit (&((LibOrder *) ordeptr)->o);
}

/*+ This routine loads the contents of
*** the given ordering from the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderLoad (
const SCOTCH_Graph * const        grafptr,        /*+ Graph to order   +*/
SCOTCH_Ordering * restrict const  ordeptr,        /*+ Ordering to load +*/
FILE * restrict const             stream)         /*+ Output stream    +*/
{
  const Graph *       srcgrafptr;
  LibOrder *          libordeptr;

  srcgrafptr = (Graph *) grafptr;
  libordeptr = (LibOrder *) ordeptr;

  if (orderLoad (&libordeptr->o, srcgrafptr->vlbltax, stream) != 0)
    return (1);

  if (libordeptr->permtab != NULL)                /* Build inverse permutation if wanted */
    orderPeri (libordeptr->o.peritab, srcgrafptr->baseval, libordeptr->o.vnodnbr, libordeptr->permtab, srcgrafptr->baseval);

  return (0);
}

/*+ This routine saves the contents of
*** the given ordering to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderSave (
const SCOTCH_Graph * const    grafptr,            /*+ Graph to order   +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSave (&((LibOrder *) ordeptr)->o, ((Graph *) grafptr)->vlbltax + ((Graph *) grafptr)->baseval, stream));
}

/*+ This routine saves to the given stream
*** the mapping data associated with the
*** given ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderSaveMap (
const SCOTCH_Graph * const    grafptr,            /*+ Graph to order   +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSaveMap (&((LibOrder *) ordeptr)->o, ((Graph *) grafptr)->vlbltax + ((Graph *) grafptr)->baseval, stream));
}

/*+ This routine saves to the given stream
*** the separator tree data associated with
*** the given ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderSaveTree (
const SCOTCH_Graph * const    grafptr,            /*+ Graph to order   +*/
const SCOTCH_Ordering * const ordeptr,            /*+ Ordering to save +*/
FILE * const                  stream)             /*+ Output stream    +*/
{
  return (orderSaveTree (&((LibOrder *) ordeptr)->o, ((Graph *) grafptr)->vlbltax + ((Graph *) grafptr)->baseval, stream));
}

/*+ This routine computes an ordering
*** of the API ordering structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderCompute (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order      +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering to compute +*/
SCOTCH_Strat * const        stratptr)             /*+ Ordering strategy   +*/
{
  return (SCOTCH_graphOrderComputeList (grafptr, ordeptr, ((Graph *) grafptr)->vertnbr, NULL, stratptr));
}

/*+ This routine computes a partial ordering
*** of the listed vertices of the API ordering
*** structure graph with respect to the given
*** strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderComputeList (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order                  +*/
SCOTCH_Ordering * const     ordeptr,              /*+ Ordering to compute             +*/
const SCOTCH_Num            listnbr,              /*+ Number of vertices in list      +*/
const SCOTCH_Num * const    listtab,              /*+ List of vertex indices to order +*/
SCOTCH_Strat * const        stratptr)             /*+ Ordering strategy               +*/
{
  const Graph * restrict  srcgrafptr;
  LibOrder *              libordeptr;             /* Pointer to ordering             */
  Hgraph                  halgrafdat;             /* Halo source graph structure     */
  Hgraph                  halgraftmp;             /* Halo source graph structure     */
  Hgraph *                halgrafptr;             /* Pointer to halo graph structure */
  const Strat *           ordstratptr;            /* Pointer to ordering strategy    */
  OrderCblk *             cblkptr;

  srcgrafptr = (Graph *) grafptr;
  libordeptr = (LibOrder *) ordeptr;              /* Get ordering */

#ifdef SCOTCH_DEBUG_LIBRARY1
  if ((listnbr < 0) || (listnbr > srcgrafptr->vertnbr)) {
    errorPrint ("SCOTCH_graphOrderComputeList: invalid parameters (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
#ifdef SCOTCH_DEBUG_LIBRARY2
  if (graphCheck (srcgrafptr) != 0) {
    errorPrint ("SCOTCH_graphOrderComputeList: invalid input graph");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY2 */

  if (listnbr == 0) {                             /* If empty list, return identity peremutation */
    intAscn (libordeptr->o.peritab, srcgrafptr->vertnbr, srcgrafptr->baseval);
    return  (0);
  }

  if (*((Strat **) stratptr) == NULL)             /* Set default ordering strategy if necessary */
    SCOTCH_stratGraphOrderBuild (stratptr, SCOTCH_STRATQUALITY, 0, 0.2);

  ordstratptr = *((Strat **) stratptr);
  if (ordstratptr->tabl != &hgraphorderststratab) {
    errorPrint ("SCOTCH_graphOrderComputeList: not an ordering strategy");
    return     (1);
  }

  memCpy (&halgrafdat.s, grafptr, sizeof (Graph)); /* Copy non-halo graph data   */
  halgrafdat.s.flagval &= ~GRAPHFREETABS;         /* Do not allow to free arrays */
  halgrafdat.s.edlotax  = NULL;                   /* Don't mind about edge loads */
  halgrafdat.vnohnbr    = halgrafdat.s.vertnbr;   /* All vertices are non-halo   */
  halgrafdat.vnohnnd    = halgrafdat.s.vertnnd;   /* No halo present             */
  halgrafdat.vnhdtax    = halgrafdat.s.vendtax;   /* End of non-halo vertices    */
  halgrafdat.vnlosum    = halgrafdat.s.velosum;   /* Sum of node vertex weights  */
  halgrafdat.enohnbr    = halgrafdat.s.edgenbr;   /* No halo present             */
  halgrafdat.enohsum    = halgrafdat.s.edlosum;
  halgrafdat.levlnum    = 0;                      /* No nested dissection yet */

  if (listnbr == srcgrafptr->vertnbr) {           /* If work on full graph */
    halgrafptr = &halgrafdat;
    cblkptr    = &libordeptr->o.cblktre;
  }
  else {
    VertList              listdat;
    Gnum * restrict       peritax;
    Gnum                  listnum;
    Gnum                  vertnum;
    Gnum                  halonum;

    if ((cblkptr = (OrderCblk *) memAlloc (2 * sizeof (OrderCblk))) == NULL) {
      errorPrint ("SCOTCH_graphOrderComputeList: out of memory");
      return     (1);
    }
    libordeptr->o.treenbr = 3;
    libordeptr->o.cblknbr = 2;
    libordeptr->o.cblktre.typeval = ORDERCBLKNEDI; /* Node becomes a (fake) nested dissection node */
    libordeptr->o.cblktre.vnodnbr = srcgrafptr->vertnbr;
    libordeptr->o.cblktre.cblknbr = 2;
    libordeptr->o.cblktre.cblktab = cblkptr;

    cblkptr[0].typeval = ORDERCBLKOTHR;           /* Build column blocks */
    cblkptr[0].vnodnbr = listnbr;
    cblkptr[0].cblknbr = 0;
    cblkptr[0].cblktab = NULL;
    cblkptr[1].typeval = ORDERCBLKOTHR;
    cblkptr[1].vnodnbr = srcgrafptr->vertnbr - listnbr;
    cblkptr[1].cblknbr = 0;
    cblkptr[1].cblktab = NULL;

    memSet (libordeptr->o.peritab, 0, srcgrafptr->vertnbr * sizeof (Gnum)); /* Fill inverse permutation with dummy values */
    for (listnum = 0, peritax = libordeptr->o.peritab - srcgrafptr->baseval;
         listnum < listnbr; listnum ++) {
#ifdef SCOTCH_DEBUG_LIBRARY2
      if ((listtab[listnum] <  srcgrafptr->baseval) ||
          (listtab[listnum] >= srcgrafptr->vertnnd)) {
        errorPrint ("SCOTCH_graphOrderComputeList: invalid parameters (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_LIBRARY2 */
      peritax[listtab[listnum]] = ~0;             /* TRICK: use peritab as flag array to mark used vertices */
    }
    for (vertnum = halonum = srcgrafptr->vertnnd - 1; vertnum >= srcgrafptr->baseval; vertnum --) {
      if (peritax[vertnum] == 0)
        peritax[halonum --] = vertnum;
    }
#ifdef SCOTCH_DEBUG_LIBRARY2
    if (halonum != (listnbr + srcgrafptr->baseval - 1)) { 
      errorPrint ("SCOTCH_graphOrderComputeList: internal error");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_LIBRARY2 */

    listdat.vnumnbr = listnbr;
    listdat.vnumtab = (Gnum * const) listtab;
    if (hgraphInduceList (&halgrafdat, &listdat, srcgrafptr->vertnbr - listnbr, &halgraftmp) != 0) {
      errorPrint ("SCOTCH_graphOrderComputeList: cannot create induced subgraph");
      return     (1);
    }
    halgrafptr = &halgraftmp;
  }

  intRandInit ();                                 /* Check that random number generator is initialized */

  hgraphOrderSt (halgrafptr, &libordeptr->o, 0, cblkptr, ordstratptr);

  if (halgrafptr != &halgrafdat)                  /* If induced subgraph created */
    hgraphExit (halgrafptr);                      /* Free it                     */

#ifdef SCOTCH_DEBUG_LIBRARY2
  if (orderCheck (&libordeptr->o) != 0)
    return (1);
#endif /* SCOTCH_DEBUG_LIBRARY2 */

  if (libordeptr->permtab != NULL)                /* Build direct permutation if wanted */
    orderPeri (libordeptr->o.peritab, srcgrafptr->baseval, libordeptr->o.vnodnbr, libordeptr->permtab, srcgrafptr->baseval);
  if (libordeptr->rangtab != NULL)                /* Build range array if column block data wanted */
    orderRang (&libordeptr->o, libordeptr->rangtab);
  if (libordeptr->treetab != NULL)                /* Build separator tree array if wanted */
      orderTree (&libordeptr->o, libordeptr->treetab);
  if (libordeptr->cblkptr != NULL)                /* Set number of column blocks if wanted */
    *(libordeptr->cblkptr) = libordeptr->o.cblknbr;

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
SCOTCH_graphOrder (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order                     +*/
SCOTCH_Strat * const        stratptr,             /*+ Ordering strategy                  +*/
SCOTCH_Num * const          permtab,              /*+ Ordering permutation               +*/
SCOTCH_Num * const          peritab,              /*+ Inverse permutation array          +*/
SCOTCH_Num * const          cblkptr,              /*+ Pointer to number of column blocks +*/
SCOTCH_Num * const          rangtab,              /*+ Column block range array           +*/
SCOTCH_Num * const          treetab)              /*+ Separator tree array               +*/
{
  SCOTCH_Ordering     ordedat;
  int                 o;

  if (SCOTCH_graphOrderInit (grafptr, &ordedat, permtab, peritab, cblkptr, rangtab, treetab) != 0)
    return (1);

  o = SCOTCH_graphOrderCompute (grafptr, &ordedat, stratptr);
  SCOTCH_graphOrderExit (grafptr, &ordedat);

  return (o);
}

/*+ This routine computes an ordering
*** of the subgraph of the API ordering
*** structure graph induced by the given
*** vertex list, with respect to the given
*** strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderList (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order                     +*/
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

  SCOTCH_graphOrderInit (grafptr, &ordedat, permtab, peritab, cblkptr, rangtab, treetab);
  o = SCOTCH_graphOrderComputeList (grafptr, &ordedat, listnbr, listtab, stratptr);
  SCOTCH_graphOrderExit (grafptr, &ordedat);

  return (o);
}

/*+ This routine checks the consistency
*** of the given graph ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphOrderCheck (
const SCOTCH_Graph * const    grafptr,
const SCOTCH_Ordering * const ordeptr)            /*+ Ordering to check +*/
{
  return (orderCheck (&((LibOrder *) ordeptr)->o));
}

/*+ This routine parses the given
*** graph ordering strategy.
*** It returns:
*** - 0   : if string successfully scanned.
*** - !0  : on error.
+*/

int
SCOTCH_stratGraphOrder (
SCOTCH_Strat * const        stratptr,
const char * const          string)
{
  if (*((Strat **) stratptr) != NULL)
    stratExit (*((Strat **) stratptr));

  if ((*((Strat **) stratptr) = stratInit (&hgraphorderststratab, string)) == NULL) {
    errorPrint ("SCOTCH_stratGraphOrder: error in ordering strategy");
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
SCOTCH_stratGraphOrderBuild (
SCOTCH_Strat * const        stratptr,             /*+ Strategy to create                 +*/
const SCOTCH_Num            flagval,              /*+ Desired characteristics            +*/
const SCOTCH_Num            levlnbr,              /*+ Number of nested dissection levels +*/
const double                balrat)               /*+ Desired imbalance ratio            +*/
{
  char                bufftab[8192];              /* Should be enough */
  char                levltab[32];
  char                bbaltab[32];
  char *              sepaptr;
  char *              tstsptr;
  char *              oleaptr;
  char *              osepptr;

  sprintf (bbaltab, "%lf", balrat);
  sprintf (levltab, GNUMSTRING, levlnbr);

  strcpy (bufftab, "c{rat=0.7,cpr=n{sep=/(<TSTS>)?m{rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=<BBAL>},org=(|h{pass=10})f{bal=<BBAL>}}}<SEPA>;,ole=<OLEA>,ose=<OSEP>},unc=n{sep=/(<TSTS>)?m{rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=<BBAL>},org=(|h{pass=10})f{bal=<BBAL>}}}<SEPA>;,ole=<OLEA>,ose=<OSEP>}}");

  switch (flagval & (SCOTCH_STRATLEVELMIN | SCOTCH_STRATLEVELMAX)) {
    case SCOTCH_STRATLEVELMIN :
      tstsptr = "(levl<<LEVL>)|(vert>240)";
      break;
    case SCOTCH_STRATLEVELMAX :
      tstsptr = "(levl<<LEVL>)&(vert>240)";
      break;
    case (SCOTCH_STRATLEVELMIN | SCOTCH_STRATLEVELMAX) :
      tstsptr = "levl<<LEVL>";
      break;
    default :
      tstsptr = "vert>240";
      break;
  }

  sepaptr = ((flagval & SCOTCH_STRATSPEED) != 0)
            ? ""
            : "|m{rat=0.7,vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=<BBAL>},org=(|h{pass=10})f{bal=<BBAL>}}}";

  oleaptr = ((flagval & SCOTCH_STRATLEAFSIMPLE) != 0)
            ? "s"
            : "f{cmin=15,cmax=100000,frat=0.0}";

  osepptr = ((flagval & SCOTCH_STRATSEPASIMPLE) != 0)
            ? "s"
            : "g";

  stringSubst (bufftab, "<SEPA>", sepaptr);
  stringSubst (bufftab, "<TSTS>", tstsptr);
  stringSubst (bufftab, "<LEVL>", levltab);
  stringSubst (bufftab, "<OLEA>", oleaptr);
  stringSubst (bufftab, "<OSEP>", osepptr);
  stringSubst (bufftab, "<BBAL>", bbaltab);

  if (SCOTCH_stratGraphOrder (stratptr, bufftab) != 0) {
    errorPrint ("SCOTCH_stratGraphOrderBuild: error in sequential ordering strategy");
    return     (1);
  }

  return (0);
}
