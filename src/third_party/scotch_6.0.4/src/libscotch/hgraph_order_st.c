/* Copyright 2004,2007,2008,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_st.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the generic call to the  **/
/**                graph ordering module, using a given    **/
/**                strategy.                               **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 oct 1996     **/
/**                                 to     09 sep 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     07 sep 2001     **/
/**                # Version 4.0  : from : 27 dec 2001     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.0  : from : 31 may 2008     **/
/**                                 to     31 may 2008     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to     17 oct 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_ST

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_bl.h"
#include "hgraph_order_cp.h"
#include "hgraph_order_gp.h"
#include "hgraph_order_hd.h"
#include "hgraph_order_hf.h"
#include "hgraph_order_kp.h"
#include "hgraph_order_nd.h"
#include "hgraph_order_si.h"
#include "hgraph_order_st.h"
#include "kgraph.h"
#include "kgraph_map_st.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"

/*
**  The static and global variables.
*/

static Hgraph               hgraphorderstgraphdummy; /* Dummy graph for offset computations */

static union {                                    /* Default parameters for block splitting method */
  HgraphOrderBlParam        param;                /* Parameter zone                                */
  StratNodeMethodData       padding;              /* To avoid reading out of structure             */
} hgraphorderstdefaultbl = { { &stratdummy, 8 } };

static union {
  HgraphOrderCpParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaultcp = { { 0.70L, &stratdummy, &stratdummy } };

static union {
  HgraphOrderGpParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaultgp = { { 3 } };

static union {
  HgraphOrderHdParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaulthd = { { 1, 10000, 0.08L } };

static union {
  HgraphOrderHfParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaulthf = { { 1, 1000000, 0.08L } };

static union {
  HgraphOrderKpParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaultkp = { { 1, &stratdummy } };

static union {                                    /* Default parameters for nested dissection method */
  HgraphOrderNdParam        param;
  StratNodeMethodData       padding;
} hgraphorderstdefaultnd = { { &stratdummy, &stratdummy, &stratdummy } };

static StratMethodTab       hgraphorderstmethtab[] = { /* Graph ordering methods array */
                              { HGRAPHORDERSTMETHBL, "b",  hgraphOrderBl, &hgraphorderstdefaultbl },
                              { HGRAPHORDERSTMETHCP, "c",  hgraphOrderCp, &hgraphorderstdefaultcp },
                              { HGRAPHORDERSTMETHGP, "g",  hgraphOrderGp, &hgraphorderstdefaultgp },
                              { HGRAPHORDERSTMETHHD, "d",  hgraphOrderHd, &hgraphorderstdefaulthd },
                              { HGRAPHORDERSTMETHHF, "f",  hgraphOrderHf, &hgraphorderstdefaulthf },
                              { HGRAPHORDERSTMETHKP, "k",  hgraphOrderKp, &hgraphorderstdefaultkp },
                              { HGRAPHORDERSTMETHND, "n",  hgraphOrderNd, &hgraphorderstdefaultnd },
                              { HGRAPHORDERSTMETHSI, "s",  hgraphOrderSi, NULL },
                              { -1,                  NULL, NULL,          NULL } };

static StratParamTab        hgraphorderstparatab[] = { /* The method parameter list */
                              { HGRAPHORDERSTMETHBL,  STRATPARAMSTRAT,  "strat",
                                (byte *) &hgraphorderstdefaultbl.param,
                                (byte *) &hgraphorderstdefaultbl.param.strat,
                                (void *) &hgraphorderststratab },
                              { HGRAPHORDERSTMETHBL,  STRATPARAMINT,    "cmin",
                                (byte *) &hgraphorderstdefaultbl.param,
                                (byte *) &hgraphorderstdefaultbl.param.cblkmin,
                                NULL },
                              { HGRAPHORDERSTMETHCP,  STRATPARAMDOUBLE, "rat",
                                (byte *) &hgraphorderstdefaultcp.param,
                                (byte *) &hgraphorderstdefaultcp.param.comprat,
                                NULL },
                              { HGRAPHORDERSTMETHCP,  STRATPARAMSTRAT,  "cpr",
                                (byte *) &hgraphorderstdefaultcp.param,
                                (byte *) &hgraphorderstdefaultcp.param.stratcpr,
                                (void *) &hgraphorderststratab },
                              { HGRAPHORDERSTMETHCP,  STRATPARAMSTRAT,  "unc",
                                (byte *) &hgraphorderstdefaultcp.param,
                                (byte *) &hgraphorderstdefaultcp.param.stratunc,
                                (void *) &hgraphorderststratab },
                              { HGRAPHORDERSTMETHGP,  STRATPARAMINT,    "pass",
                                (byte *) &hgraphorderstdefaultgp.param,
                                (byte *) &hgraphorderstdefaultgp.param.passnbr,
                                NULL },
                              { HGRAPHORDERSTMETHHD,  STRATPARAMINT,    "cmin",
                                (byte *) &hgraphorderstdefaulthd.param,
                                (byte *) &hgraphorderstdefaulthd.param.colmin,
                                NULL },
                              { HGRAPHORDERSTMETHHD,  STRATPARAMINT,    "cmax",
                                (byte *) &hgraphorderstdefaulthd.param,
                                (byte *) &hgraphorderstdefaulthd.param.colmax,
                                NULL },
                              { HGRAPHORDERSTMETHHD,  STRATPARAMDOUBLE, "frat",
                                (byte *) &hgraphorderstdefaulthd.param,
                                (byte *) &hgraphorderstdefaulthd.param.fillrat,
                                NULL },
                              { HGRAPHORDERSTMETHHF,  STRATPARAMINT,    "cmin",
                                (byte *) &hgraphorderstdefaulthf.param,
                                (byte *) &hgraphorderstdefaulthf.param.colmin,
                                NULL },
                              { HGRAPHORDERSTMETHHF,  STRATPARAMINT,    "cmax",
                                (byte *) &hgraphorderstdefaulthf.param,
                                (byte *) &hgraphorderstdefaulthf.param.colmax,
                                NULL },
                              { HGRAPHORDERSTMETHHF,  STRATPARAMDOUBLE, "frat",
                                (byte *) &hgraphorderstdefaulthf.param,
                                (byte *) &hgraphorderstdefaulthf.param.fillrat,
                                NULL },
                              { HGRAPHORDERSTMETHKP,  STRATPARAMINT,    "siz",
                                (byte *) &hgraphorderstdefaultkp.param,
                                (byte *) &hgraphorderstdefaultkp.param.partsiz,
                                NULL },
                              { HGRAPHORDERSTMETHKP,  STRATPARAMSTRAT,  "strat",
                                (byte *) &hgraphorderstdefaultkp.param,
                                (byte *) &hgraphorderstdefaultkp.param.strat,
                                (void *) &kgraphmapststratab },
                              { HGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "sep",
                                (byte *) &hgraphorderstdefaultnd.param,
                                (byte *) &hgraphorderstdefaultnd.param.sepstrat,
                                (void *) &vgraphseparateststratab },
                              { HGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "ole",
                                (byte *) &hgraphorderstdefaultnd.param,
                                (byte *) &hgraphorderstdefaultnd.param.ordstratlea,
                                (void *) &hgraphorderststratab },
                              { HGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "ose",
                                (byte *) &hgraphorderstdefaultnd.param,
                                (byte *) &hgraphorderstdefaultnd.param.ordstratsep,
                                (void *) &hgraphorderststratab },
                              { HGRAPHORDERSTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        hgraphorderstcondtab[] = { /* Graph condition parameter table */
                              { STRATNODECOND,        STRATPARAMINT,    "edge",
                                (byte *) &hgraphorderstgraphdummy,
                                (byte *) &hgraphorderstgraphdummy.s.edgenbr,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,     "levl",
                                (byte *) &hgraphorderstgraphdummy,
                                (byte *) &hgraphorderstgraphdummy.levlnum,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "load",
                                (byte *) &hgraphorderstgraphdummy,
                                (byte *) &hgraphorderstgraphdummy.vnlosum,
                                NULL },
                              { STRATNODECOND,        STRATPARAMDOUBLE, "mdeg",
                                (byte *) &hgraphorderstgraphdummy,
                                (byte *) &hgraphorderstgraphdummy.s.degrmax,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "vert",
                                (byte *) &hgraphorderstgraphdummy,
                                (byte *) &hgraphorderstgraphdummy.vnohnbr, /* Only consider non-halo vertices */
                                NULL },
                              { STRATNODENBR,         STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    hgraphorderststratab = { /* Strategy tables for graph ordering methods */
                              hgraphorderstmethtab,
                              hgraphorderstparatab,
                              hgraphorderstcondtab };

/************************************/
/*                                  */
/* This routine is the entry point  */
/* for the graph ordering routines. */
/*                                  */
/************************************/

/* This routine computes an ordering
** with respect to a given strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hgraphOrderSt (
const Hgraph * restrict const   grafptr,          /*+ Subgraph to order          +*/
Order * restrict const          ordeptr,          /*+ Ordering to complete       +*/
const Gnum                      ordenum,          /*+ Index to start ordering at +*/
OrderCblk * restrict const      cblkptr,          /*+ Current column block       +*/
const Strat * restrict const    strat)            /*+ Graph ordering strategy    +*/
{
  StratTest           val;
  int                 o;

  if (grafptr->vnohnbr == 0)                      /* Return immediately if nothing to do */
    return (0);

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      errorPrint ("hgraphOrderSt: concatenation operator not available for graph ordering strategies");
      return     (1);
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_HGRAPH2
        if ((val.typetest != STRATTESTVAL) &&
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("hgraphOrderSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_HGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                                              */
          o = hgraphOrderSt (grafptr, ordeptr, ordenum, cblkptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                                        */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement                                  */
            o = hgraphOrderSt (grafptr, ordeptr, ordenum, cblkptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      hgraphOrderSi (grafptr, ordeptr, ordenum, cblkptr); /* Always maintain a consistent ordering */
      break;
    case STRATNODESELECT :
      errorPrint ("hgraphOrderSt: selection operator not available for graph ordering strategies");
      return     (1);
#ifdef SCOTCH_DEBUG_HGRAPH2
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_HGRAPH2 */
    default :
#endif /* SCOTCH_DEBUG_HGRAPH2 */
      return (strat->tabl->methtab[strat->data.method.meth].func (grafptr, ordeptr, ordenum, cblkptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_HGRAPH2
    default :
      errorPrint ("hgraphOrderSt: invalid parameter");
      return     (1);
#endif /* SCOTCH_DEBUG_HGRAPH2 */
  }
  return (o);
}
