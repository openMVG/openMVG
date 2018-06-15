/* Copyright 2008-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kdgraph_map_st.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the strategy and   **/
/**                method tables for the parallel          **/
/**                multi-way static mapping routines.      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 jun 2008     **/
/**                                 to     14 apr 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KDGRAPH_MAP_ST

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "arch.h"
#include "mapping.h"
#include "dmapping.h"
#include "bdgraph.h"
#include "bdgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_st.h"
#include "kdgraph.h"
#include "kdgraph_map_rb.h"
#include "kdgraph_map_st.h"

/*
**  The static and global variables.
*/

static union {
  KdgraphMapRbParam         param;
  StratNodeMethodData       padding;
} kdgraphmapstdefaultrb = { { &stratdummy, &stratdummy, 0.05 } };

static StratMethodTab       kdgraphmapstmethtab[] = { /* Mapping methods array */
                              { KDGRAPHMAPSTMETHRB, "r",  kdgraphMapRb, &kdgraphmapstdefaultrb },
                              { -1,                 NULL, NULL,         NULL } };

static StratParamTab        kdgraphmapstparatab[] = { /* Method parameter list */
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMSTRAT,  "sep",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.stratsep,
                                (void *) &bdgraphbipartststratab },
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMSTRAT,  "seq",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.stratseq,
                                (void *) &kgraphmapststratab },
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMDOUBLE, "bal",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.kbalval,
                                NULL },
                              { KDGRAPHMAPSTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        kdgraphmapstcondtab[] = { /* Graph condition parameter table */
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    kdgraphmapststratab = { /* Strategy tables for graph mapping methods */
                              kdgraphmapstmethtab,
                              kdgraphmapstparatab,
                              kdgraphmapstcondtab };

/****************************************/
/*                                      */
/* This is the generic mapping routine. */
/*                                      */
/****************************************/

/* This routine computes the given
** mapping according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kdgraphMapSt (
Kdgraph * restrict const      grafptr,            /*+ Mapping graph    +*/
Kdmapping * restrict const    mappptr,            /*+ Dynamic mapping  +*/
const Strat * restrict const  strat)              /*+ Mapping strategy +*/
{
  StratTest           val;
  int                 o;

#ifdef SCOTCH_DEBUG_KDGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("kdgraphMapSt: invalid type specification for parser variables");
    return     (1);
  }
  if ((sizeof (KdgraphMapRbParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("kdgraphMapSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KDGRAPH2 */
#ifdef SCOTCH_DEBUG_KDGRAPH1
  if ((strat->tabl != &kdgraphmapststratab) &&
      (strat       != &stratdummy)) {
    errorPrint ("kdgraphMapSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KDGRAPH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = kdgraphMapSt (grafptr, mappptr, strat->data.concat.strat[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                              */
        o |= kdgraphMapSt (grafptr, mappptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_KDGRAPH2
        if ((val.typetest != STRATTESTVAL) ||
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("kdgraphMapSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_KDGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                           */
          o = kdgraphMapSt (grafptr, mappptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                     */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement               */
            o = kdgraphMapSt (grafptr, mappptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      errorPrint ("kdgraphMapSt: selection operator not implemented for k-way strategies");
      return      (1);
#ifdef SCOTCH_DEBUG_KDGRAPH1
    case STRATNODEMETHOD :
#else  /* SCOTCH_DEBUG_KDGRAPH1 */
    default :
#endif /* SCOTCH_DEBUG_KDGRAPH1 */
      return (strat->tabl->methtab[strat->data.method.meth].func (grafptr, mappptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_KDGRAPH1
    default :
      errorPrint ("kdgraphMapSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_KDGRAPH1 */
  }
  return (o);
}
