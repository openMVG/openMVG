/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hmesh_order_st.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the generic call to the  **/
/**                halo mesh ordering module, using a      **/
/**                given strategy.                         **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 28 sep 2002     **/
/**                                 to     05 jan 2005     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_ORDER_ST

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "hgraph.h"
#include "mesh.h"
#include "hmesh.h"
#include "order.h"
#include "hgraph_order_st.h"
#include "hmesh_order_bl.h"
#include "hmesh_order_cp.h"
#include "hmesh_order_gp.h"
#include "hmesh_order_gr.h"
#include "hmesh_order_hd.h"
#include "hmesh_order_hf.h"
#include "hmesh_order_nd.h"
#include "hmesh_order_si.h"
#include "hmesh_order_st.h"
#include "vmesh.h"
#include "vmesh_separate_st.h"

/*
**  The static and global variables.
*/

static Hmesh                hmeshorderstmeshdummy; /* Dummy mesh for offset computations */

static union {                                    /* Default parameters for block splitting method */
  HmeshOrderBlParam         param;                /* Parameter zone                                */
  StratNodeMethodData       padding;              /* To avoid reading out of structure             */
} hmeshorderstdefaultbl = { { &stratdummy, 8 } };

static union {
  HmeshOrderCpParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaultcp = { { 0.70L, &stratdummy, &stratdummy } };

static union {                                    /* Default parameters for nested dissection method */
  HmeshOrderGpParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaultgp = { { 3 } };

static union {                                    /* Default parameters for nested dissection method */
  HmeshOrderGrParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaultgr = { { &stratdummy } };

static union {
  HmeshOrderHdParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaulthd = { { 1, 1000000, 0.08L } };

static union {
  HmeshOrderHfParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaulthf = { { 1, 1000000, 0.08L } };

static union {                                    /* Default parameters for nested dissection method */
  HmeshOrderNdParam         param;
  StratNodeMethodData       padding;
} hmeshorderstdefaultnd = { { &stratdummy, &stratdummy, &stratdummy } };

static StratMethodTab       hmeshorderstmethtab[] = { /* Mesh ordering methods array */
                              { HMESHORDERSTMETHBL, "b",  hmeshOrderBl, &hmeshorderstdefaultbl },
                              { HMESHORDERSTMETHCP, "c",  hmeshOrderCp, &hmeshorderstdefaultcp },
                              { HMESHORDERSTMETHGP, "g",  hmeshOrderGp, &hmeshorderstdefaultgp },
                              { HMESHORDERSTMETHGR, "v",  hmeshOrderGr, &hmeshorderstdefaultgr },
                              { HMESHORDERSTMETHHD, "d",  hmeshOrderHd, &hmeshorderstdefaulthd },
                              { HMESHORDERSTMETHHF, "f",  hmeshOrderHf, &hmeshorderstdefaulthf },
                              { HMESHORDERSTMETHND, "n",  hmeshOrderNd, &hmeshorderstdefaultnd },
                              { HMESHORDERSTMETHSI, "s",  hmeshOrderSi, NULL },
                              { -1,                 NULL, NULL,         NULL } };

static StratParamTab        hmeshorderstparatab[] = { /* The method parameter list */
                              { HMESHORDERSTMETHBL,   STRATPARAMSTRAT,  "strat",
                                (byte *) &hmeshorderstdefaultbl.param,
                                (byte *) &hmeshorderstdefaultbl.param.strat,
                                (void *) &hmeshorderststratab },
                              { HMESHORDERSTMETHBL,   STRATPARAMINT,    "cmin",
                                (byte *) &hmeshorderstdefaultbl.param,
                                (byte *) &hmeshorderstdefaultbl.param.cblkmin,
                                NULL },
                              { HMESHORDERSTMETHCP,   STRATPARAMDOUBLE, "rat",
                                (byte *) &hmeshorderstdefaultcp.param,
                                (byte *) &hmeshorderstdefaultcp.param.comprat,
                                NULL },
                              { HMESHORDERSTMETHCP,   STRATPARAMSTRAT,  "cpr",
                                (byte *) &hmeshorderstdefaultcp.param,
                                (byte *) &hmeshorderstdefaultcp.param.stratcpr,
                                (void *) &hmeshorderststratab },
                              { HMESHORDERSTMETHCP,   STRATPARAMSTRAT,  "unc",
                                (byte *) &hmeshorderstdefaultcp.param,
                                (byte *) &hmeshorderstdefaultcp.param.stratunc,
                                (void *) &hmeshorderststratab },
                              { HMESHORDERSTMETHGP,   STRATPARAMINT,    "pass",
                                (byte *) &hmeshorderstdefaultgp.param,
                                (byte *) &hmeshorderstdefaultgp.param.passnbr,
                                NULL },
                              { HMESHORDERSTMETHGR,   STRATPARAMSTRAT,  "strat",
                                (byte *) &hmeshorderstdefaultgr.param,
                                (byte *) &hmeshorderstdefaultgr.param.stratptr,
                                (void *) &hgraphorderststratab },
                              { HMESHORDERSTMETHHD,   STRATPARAMINT,    "cmin",
                                (byte *) &hmeshorderstdefaulthd.param,
                                (byte *) &hmeshorderstdefaulthd.param.colmin,
                                NULL },
                              { HMESHORDERSTMETHHD,   STRATPARAMINT,    "cmax",
                                (byte *) &hmeshorderstdefaulthd.param,
                                (byte *) &hmeshorderstdefaulthd.param.colmax,
                                NULL },
                              { HMESHORDERSTMETHHD,   STRATPARAMDOUBLE, "frat",
                                (byte *) &hmeshorderstdefaulthd.param,
                                (byte *) &hmeshorderstdefaulthd.param.fillrat,
                                NULL },
                              { HMESHORDERSTMETHHF,   STRATPARAMINT,    "cmin",
                                (byte *) &hmeshorderstdefaulthf.param,
                                (byte *) &hmeshorderstdefaulthf.param.colmin,
                                NULL },
                              { HMESHORDERSTMETHHF,   STRATPARAMINT,    "cmax",
                                (byte *) &hmeshorderstdefaulthf.param,
                                (byte *) &hmeshorderstdefaulthf.param.colmax,
                                NULL },
                              { HMESHORDERSTMETHHF,   STRATPARAMDOUBLE, "frat",
                                (byte *) &hmeshorderstdefaulthf.param,
                                (byte *) &hmeshorderstdefaulthf.param.fillrat,
                                NULL },
                              { HMESHORDERSTMETHND,   STRATPARAMSTRAT,  "sep",
                                (byte *) &hmeshorderstdefaultnd.param,
                                (byte *) &hmeshorderstdefaultnd.param.sepstrat,
                                (void *) &vmeshseparateststratab },
                              { HMESHORDERSTMETHND,   STRATPARAMSTRAT,  "ole",
                                (byte *) &hmeshorderstdefaultnd.param,
                                (byte *) &hmeshorderstdefaultnd.param.ordstratlea,
                                (void *) &hmeshorderststratab },
                              { HMESHORDERSTMETHND,   STRATPARAMSTRAT,  "ose",
                                (byte *) &hmeshorderstdefaultnd.param,
                                (byte *) &hmeshorderstdefaultnd.param.ordstratsep,
                                (void *) &hmeshorderststratab },
                              { HMESHORDERSTMETHNBR,  STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        hmeshorderstcondtab[] = { /* Mesh condition parameter table */
                              { STRATNODECOND,        STRATPARAMINT,    "edge",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.m.edgenbr,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "levl",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.levlnum,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "load",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.vnhlsum,
                                NULL },
                              { STRATNODECOND,        STRATPARAMDOUBLE, "mdeg",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.m.degrmax,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "vnod",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.vnohnbr,
                                NULL },
                              { STRATNODECOND,        STRATPARAMINT,    "velm",
                                (byte *) &hmeshorderstmeshdummy,
                                (byte *) &hmeshorderstmeshdummy.m.velmnbr,
                                NULL },
                              { STRATNODENBR,         STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    hmeshorderststratab = { /* Strategy tables for mesh ordering methods */
                              hmeshorderstmethtab,
                              hmeshorderstparatab,
                              hmeshorderstcondtab };

/***********************************/
/*                                 */
/* This routine is the entry point */
/* for the mesh ordering routines. */
/*                                 */
/***********************************/

/* This routine computes an ordering
** with respect to a given strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hmeshOrderSt (
const Hmesh * restrict const    meshptr,          /*+ Submesh to which list apply +*/
Order * restrict const          ordeptr,          /*+ Ordering to complete        +*/
const Gnum                      ordenum,          /*+ Index to start ordering at  +*/
OrderCblk * restrict const      cblkptr,          /*+ Current column block        +*/
const Strat * restrict const    strat)            /*+ Mesh ordering strategy      +*/
{
  StratTest           val;
  int                 o;

  if (meshptr->vnohnbr == 0)                      /* Return immediately if nothing to do */
    return (0);

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      errorPrint ("hmeshOrderSt: concatenation operator not implemented for ordering strategies");
      return     (1);
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) meshptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_HMESH2
        if ((val.typetest != STRATTESTVAL) &&
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("hmeshOrderSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_HMESH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                                             */
          o = hmeshOrderSt (meshptr, ordeptr, ordenum, cblkptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                                       */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement                                 */
            o = hmeshOrderSt (meshptr, ordeptr, ordenum, cblkptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      hmeshOrderSi (meshptr, ordeptr, ordenum, cblkptr); /* Always maintain a coherent ordering */
      break;
    case STRATNODESELECT :
      errorPrint ("hmeshOrderSt: selection operator not available for mesh ordering strategies");
      return     (1);
#ifdef SCOTCH_DEBUG_HMESH2
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_HMESH2 */
    default :
#endif /* SCOTCH_DEBUG_HMESH2 */
      return (strat->tabl->methtab[strat->data.method.meth].func (meshptr, ordeptr, ordenum, cblkptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_HMESH2
    default :
      errorPrint ("hmeshOrderSt: invalid parameter");
      return     (1);
#endif /* SCOTCH_DEBUG_HMESH2 */
  }
  return (o);
}
