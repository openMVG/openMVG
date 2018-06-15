/* Copyright 2007-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : wgraph_part_st.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the global         **/
/**                vertex overlapped graph partitioning    **/ 
/**                strategy and method tables.             **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to     16 apr 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_ST

#include "module.h"
#include "common.h"
#include "gain.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "graph_coarsen.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "wgraph.h"
#include "wgraph_part_fm.h"
#include "wgraph_part_gg.h"
#include "wgraph_part_gp.h"
#include "wgraph_part_ml.h"
#include "wgraph_part_rb.h"
#include "wgraph_part_st.h"
#include "wgraph_part_zr.h"

/*
**  The static and global variables.
*/

static Wgraph               wgraphdummy;          /* Dummy overlap graph for offset computations */

static union {
  WgraphPartFmParam         param;
  StratNodeMethodData       padding;
} wgraphpartdefaultfm = { { 10, 40, 0.1L } };

static union {
  WgraphPartGgParam         param;
  StratNodeMethodData       padding;
} wgraphpartdefaultgg = { { 10 } };

static union {
  WgraphPartGpParam         param;
  StratNodeMethodData       padding;
} wgraphpartdefaultgp = { { 5 } };

static union {
  WgraphPartMlParam         param;
  StratNodeMethodData       padding;
} wgraphpartdefaultml = { { 20, 0.8L, &stratdummy, &stratdummy } };

static union {
  WgraphPartRbParam         param;
  StratNodeMethodData       padding;
} wgraphpartdefaultrb = { { &stratdummy } };

static StratMethodTab       wgraphpartstmethtab[] = { /* Graph overlap partitioning methods array */
                              { WGRAPHSEPASTMETHGG, "h",  wgraphPartGg, &wgraphpartdefaultgg },
                              { WGRAPHSEPASTMETHGP, "g",  wgraphPartGp, &wgraphpartdefaultgp },
                              { WGRAPHSEPASTMETHFM, "f",  wgraphPartFm, &wgraphpartdefaultfm },
                              { WGRAPHSEPASTMETHML, "m",  wgraphPartMl, &wgraphpartdefaultml },
                              { WGRAPHSEPASTMETHRB, "r",  wgraphPartRb, &wgraphpartdefaultrb },
                              { WGRAPHSEPASTMETHZR, "z",  wgraphPartZr, NULL },
                              { -1,                NULL,  NULL,         NULL } };

static StratParamTab        wgraphpartstparatab[] = { /* Method parameter list */
                              { WGRAPHSEPASTMETHFM,  STRATPARAMINT,    "pass",
                                (byte *) &wgraphpartdefaultfm.param,
                                (byte *) &wgraphpartdefaultfm.param.passnbr,
                                NULL },
                              { WGRAPHSEPASTMETHFM,  STRATPARAMINT,    "move",
                                (byte *) &wgraphpartdefaultfm.param,
                                (byte *) &wgraphpartdefaultfm.param.movenbr,
                                NULL },
                              { WGRAPHSEPASTMETHFM,  STRATPARAMDOUBLE, "bal",
                                (byte *) &wgraphpartdefaultfm.param,
                                (byte *) &wgraphpartdefaultfm.param.deltrat,
                                NULL },
                              { WGRAPHSEPASTMETHGG,  STRATPARAMINT,    "pass",
                                (byte *) &wgraphpartdefaultgg.param,
                                (byte *) &wgraphpartdefaultgg.param.passnbr,
                                NULL },
                              { WGRAPHSEPASTMETHGP,  STRATPARAMINT,    "pass",
                                (byte *) &wgraphpartdefaultgp.param,
                                (byte *) &wgraphpartdefaultgp.param.passnbr,
                                NULL },
                              { WGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "asc",
                                (byte *) &wgraphpartdefaultml.param,
                                (byte *) &wgraphpartdefaultml.param.stratasc,
                                (void *) &wgraphpartststratab },
                              { WGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "low",
                                (byte *) &wgraphpartdefaultml.param,
                                (byte *) &wgraphpartdefaultml.param.stratlow,
                                (void *) &wgraphpartststratab },
                              { WGRAPHSEPASTMETHML,  STRATPARAMINT,    "vert",
                                (byte *) &wgraphpartdefaultml.param,
                                (byte *) &wgraphpartdefaultml.param.coarnbr,
                                NULL },
                              { WGRAPHSEPASTMETHML,  STRATPARAMDOUBLE, "rat",
                                (byte *) &wgraphpartdefaultml.param,
                                (byte *) &wgraphpartdefaultml.param.coarval,
                                NULL },
                              { WGRAPHSEPASTMETHRB,  STRATPARAMSTRAT,  "sep",
                                (byte *) &wgraphpartdefaultrb.param,
                                (byte *) &wgraphpartdefaultrb.param.stratptr,
                                (void *) &vgraphseparateststratab },
                              { WGRAPHSEPASTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        wgraphpartstcondtab[] = { /* Overlap graph condition parameter table*/
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    wgraphpartststratab = { /* Strategy tables for overlap partitioning methods */
                              wgraphpartstmethtab,
                              wgraphpartstparatab,
                              wgraphpartstcondtab };

/*********************************************/
/*                                           */
/* This is the generic partitioning routine. */
/*                                           */
/*********************************************/

/* This routine computes the separation of
** the given graph according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
wgraphPartSt (
Wgraph * restrict const      grafptr,            /*+ Overlap partitioning graph    +*/
const Strat * restrict const strat)              /*+ Overlap partitioning strategy +*/
{
  StratTest           val;                        /* Result of condition evaluation */
  WgraphStore         savetab[2];                 /* Results of the two strategies  */
  int                 o;
  int                 o2;

#ifdef SCOTCH_DEBUG_WGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("wgraphPartSt: invalid type specification for parser variables");
    return     (1);
  }
  if ((sizeof (WgraphPartFmParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (WgraphPartGgParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (WgraphPartGpParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (WgraphPartMlParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (WgraphPartRbParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("wgraphPartSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_WGRAPH2 */
#ifdef SCOTCH_DEBUG_WGRAPH1
  if ((strat->tabl != &wgraphpartststratab) &&
      (strat       != &stratdummy)) {
    errorPrint ("wgraphPartSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_WGRAPH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = wgraphPartSt (grafptr, strat->data.concat.strat[0]); /* Apply the first strategy      */
      if (o == 0)                                 /* If it worked all right                     */
        o |= wgraphPartSt (grafptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_WGRAPH2
        if ((val.typetest != STRATTESTVAL) ||
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("wgraphPartSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_WGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                  */
          o = wgraphPartSt (grafptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false            */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement      */
            o = wgraphPartSt (grafptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      if (((wgraphStoreInit (grafptr, &savetab[0])) != 0) || /* Allocate save areas */
          ((wgraphStoreInit (grafptr, &savetab[1])) != 0)) {
        errorPrint ("wgraphPartSt: out of memory");
        wgraphStoreExit (&savetab[0]);
        return          (1);
      }

      wgraphStoreSave   (grafptr, &savetab[1]);   /* Save initial partition              */
      o = wgraphPartSt  (grafptr, strat->data.select.strat[0]); /* Apply first strategy  */
      wgraphStoreSave   (grafptr, &savetab[0]);   /* Save its result                     */
      wgraphStoreUpdt   (grafptr, &savetab[1]);   /* Restore initial partition           */
      o2 = wgraphPartSt (grafptr, strat->data.select.strat[1]); /* Apply second strategy */

      if ((o == 0) || (o2 == 0)) {                /* If at least one method make a k-partition */
        if (savetab[0].fronload < grafptr->fronload) /* If first strategy is better            */
          wgraphStoreUpdt (grafptr, &savetab[0]); /* Restore its result                        */
      }

      wgraphStoreExit (&savetab[0]);              /* Free both save areas */
      wgraphStoreExit (&savetab[1]);
      break;
#ifdef SCOTCH_DEBUG_WGRAPH2
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_WGRAPH2 */
    default :
#endif /* SCOTCH_DEBUG_WGRAPH2 */
      return (strat->tabl->methtab[strat->data.method.meth].func (grafptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_WGRAPH2
    default :
      errorPrint ("wgraphPartSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_WGRAPH2 */
  }
  return (o);
}
