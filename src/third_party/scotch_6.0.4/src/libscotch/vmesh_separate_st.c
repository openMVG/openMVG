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
/**   NAME       : vmesh_separate_st.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the global mesh    **/
/**                separation strategy and method tables.  **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 20 sep 2002     **/
/**                                 to     08 feb 2003     **/
/**                # Version 5.0  : from : 04 aug 2007     **/
/**                                 to     04 aug 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_SEPARATE_ST

#include "module.h"
#include "common.h"
#include "gain.h"
#include "parser.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "mesh.h"
#include "mesh_coarsen.h"
#include "vmesh.h"
#include "vmesh_separate_fm.h"
#include "vmesh_separate_gg.h"
#include "vmesh_separate_gr.h"
#include "vmesh_separate_ml.h"
#include "vmesh_separate_zr.h"
#include "vmesh_separate_st.h"

/*
**  The static and global variables.
*/

static Vmesh                vmeshdummy;           /* Dummy separator mesh for offset computations */

static union {
  VmeshSeparateFmParam      param;
  StratNodeMethodData       padding;
} vmeshseparatedefaultfm = { { 200, 1000, 0.1L } };

static union {
  VmeshSeparateGgParam      param;
  StratNodeMethodData       padding;
} vmeshseparatedefaultgg = { { 5 } };

static union {
  VmeshSeparateGrParam      param;
  StratNodeMethodData       padding;
} vmeshseparatedefaultgr = { { &stratdummy } };

static union {
  VmeshSeparateMlParam      param;
  StratNodeMethodData       padding;
} vmeshseparatedefaultml = { { 1000, 0.8L, MESHCOARSENNGB, &stratdummy, &stratdummy } };

static StratMethodTab       vmeshseparatestmethtab[] = { /* Mesh separation methods array */
                              { VMESHSEPASTMETHFM, "f",  vmeshSeparateFm, &vmeshseparatedefaultfm },
                              { VMESHSEPASTMETHGG, "h",  vmeshSeparateGg, &vmeshseparatedefaultgg },
#ifdef SCOTCH_DEBUG_VMESH2
                              { VMESHSEPASTMETHGR, "v",  vmeshSeparateGr, &vmeshseparatedefaultgr },
#endif /* SCOTCH_DEBUG_VMESH2 */
                              { VMESHSEPASTMETHML, "m",  vmeshSeparateMl, &vmeshseparatedefaultml },
                              { VMESHSEPASTMETHZR, "z",  vmeshSeparateZr, NULL },
                              { -1,                NULL, NULL,            NULL } };

static StratParamTab        vmeshseparatestparatab[] = { /* Mesh separation method parameter list */
                              { VMESHSEPASTMETHFM,  STRATPARAMINT,    "move",
                                (byte *) &vmeshseparatedefaultfm.param,
                                (byte *) &vmeshseparatedefaultfm.param.movenbr,
                                NULL },
                              { VMESHSEPASTMETHFM,  STRATPARAMINT,    "pass",
                                (byte *) &vmeshseparatedefaultfm.param,
                                (byte *) &vmeshseparatedefaultfm.param.passnbr,
                                NULL },
                              { VMESHSEPASTMETHFM,  STRATPARAMDOUBLE, "bal",
                                (byte *) &vmeshseparatedefaultfm.param,
                                (byte *) &vmeshseparatedefaultfm.param.deltrat,
                                NULL },
                              { VMESHSEPASTMETHGG,  STRATPARAMINT,    "pass",
                                (byte *) &vmeshseparatedefaultgg.param,
                                (byte *) &vmeshseparatedefaultgg.param.passnbr,
                                NULL },
#ifdef SCOTCH_DEBUG_VMESH2
                              { VMESHSEPASTMETHGR,  STRATPARAMSTRAT,  "strat",
                                (byte *) &vmeshseparatedefaultgr.param,
                                (byte *) &vmeshseparatedefaultgr.param.stratptr,
                                (void *) &vgraphseparateststratab },
#endif /* SCOTCH_DEBUG_VMESH2 */
                              { VMESHSEPASTMETHML,  STRATPARAMSTRAT,  "asc",
                                (byte *) &vmeshseparatedefaultml.param,
                                (byte *) &vmeshseparatedefaultml.param.stratasc,
                                (void *) &vmeshseparateststratab },
                              { VMESHSEPASTMETHML,  STRATPARAMSTRAT,  "low",
                                (byte *) &vmeshseparatedefaultml.param,
                                (byte *) &vmeshseparatedefaultml.param.stratlow,
                                (void *) &vmeshseparateststratab },
                              { VMESHSEPASTMETHML,  STRATPARAMCASE,   "type",
                                (byte *) &vmeshseparatedefaultml.param,
                                (byte *) &vmeshseparatedefaultml.param.coartype,
                                (void *) "hsn" },
                              { VMESHSEPASTMETHML,  STRATPARAMINT,    "vnod",
                                (byte *) &vmeshseparatedefaultml.param,
                                (byte *) &vmeshseparatedefaultml.param.vnodnbr,
                                NULL },
                              { VMESHSEPASTMETHML,  STRATPARAMDOUBLE, "rat",
                                (byte *) &vmeshseparatedefaultml.param,
                                (byte *) &vmeshseparatedefaultml.param.coarrat,
                                NULL },
                              { VMESHSEPASTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        vmeshseparatestcondtab[] = { /* Mesh condition parameter table */
                              { STRATNODECOND,      STRATPARAMINT,    "edge",
                                (byte *) &vmeshdummy,
                                (byte *) &vmeshdummy.m.edgenbr,
                                NULL },
                              { STRATNODECOND,      STRATPARAMINT,    "levl",
                                (byte *) &vmeshdummy,
                                (byte *) &vmeshdummy.levlnum,
                                NULL },
                              { STRATNODECOND,      STRATPARAMINT,    "load",
                                (byte *) &vmeshdummy,
                                (byte *) &vmeshdummy.m.vnlosum,
                                NULL },
                              { STRATNODECOND,      STRATPARAMINT,    "velm",
                                (byte *) &vmeshdummy,
                                (byte *) &vmeshdummy.m.velmnbr,
                                NULL },
                              { STRATNODECOND,      STRATPARAMINT,    "vnod",
                                (byte *) &vmeshdummy,
                                (byte *) &vmeshdummy.m.vnodnbr,
                                NULL },
                              { STRATNODENBR,       STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                      vmeshseparateststratab = { /* Strategy tables for mesh separation methods */
                                vmeshseparatestmethtab,
                                vmeshseparatestparatab,
                                vmeshseparatestcondtab };

/*******************************************/
/*                                         */
/* This is the generic separation routine. */
/*                                         */
/*******************************************/

/* This routine computes the separation of
** the given graph according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
vmeshSeparateSt (
Vmesh * restrict const        meshptr,            /*+ Separation mesh     +*/
const Strat * restrict const  strat)              /*+ Separation strategy +*/
{
  StratTest           val;
  VmeshStore          save[2];                    /* Results of the two strategies */
  int                 o;

#ifdef SCOTCH_DEBUG_VMESH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("vmeshSeparateSt: invalid type specification for parser variables");
    return     (1);
  }
  if ((sizeof (VmeshSeparateFmParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VmeshSeparateGgParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VmeshSeparateGrParam) > sizeof (StratNodeMethodData)) ||
      (sizeof (VmeshSeparateMlParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("vmeshSeparateSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */
#ifdef SCOTCH_DEBUG_VMESH1
  if (strat->tabl != &vmeshseparateststratab) {
    errorPrint ("vmeshSeparateSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = vmeshSeparateSt (meshptr, strat->data.concat.strat[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                        */
        o |= vmeshSeparateSt (meshptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) meshptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_VMESH2
        if ((val.typetest != STRATTESTVAL) &&
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("vmeshSeparateSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_VMESH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                     */
          o = vmeshSeparateSt (meshptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false               */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement         */
            o = vmeshSeparateSt (meshptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      if (((vmeshStoreInit (meshptr, &save[0])) != 0) || /* Allocate save areas */
          ((vmeshStoreInit (meshptr, &save[1])) != 0)) {
        errorPrint     ("vmeshSeparateSt: out of memory");
        vmeshStoreExit (&save[0]);
        return         (1);
      }

      vmeshStoreSave  (meshptr, &save[1]);        /* Save initial bipartition          */
      vmeshSeparateSt (meshptr, strat->data.select.strat[0]); /* Apply first strategy  */
      vmeshStoreSave  (meshptr, &save[0]);        /* Save its result                   */
      vmeshStoreUpdt  (meshptr, &save[1]);        /* Restore initial bipartition       */
      vmeshSeparateSt (meshptr, strat->data.select.strat[1]); /* Apply second strategy */

      if ( (save[0].fronnbr <  meshptr->fronnbr) || /* If first strategy is better */
          ((save[0].fronnbr == meshptr->fronnbr) &&
           (abs (save[0].ncmploaddlt) < abs (meshptr->ncmploaddlt))))
        vmeshStoreUpdt (meshptr, &save[0]);       /* Restore its result */

      vmeshStoreExit (&save[0]);                  /* Free both save areas */
      vmeshStoreExit (&save[1]);
      break;
#ifdef SCOTCH_DEBUG_VMESH1
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_VMESH1 */
    default :
#endif /* SCOTCH_DEBUG_VMESH1 */
      return (strat->tabl->methtab[strat->data.method.meth].func (meshptr, (void *) &strat->data.method.data));
#ifdef SCOTCH_DEBUG_VMESH1
    default :
      errorPrint ("vmeshSeparateSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_VMESH1 */
  }
  return (o);
}
