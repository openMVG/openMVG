/* Copyright 2007-2009,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vdgraph_separate_st.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : This module contains the global         **/
/**                distributed separation strategy and     **/
/**                method tables.                          **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 16 feb 2006     **/
/**                                 to     01 aug 2007     **/
/**                # Version 5.1  : from : 05 nov 2007     **/
/**                                 to     26 may 2009     **/
/**                # Version 6.0  : from : 01 may 2014     **/
/**                                 to     30 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_SEPARATE_ST

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "vdgraph.h"
#include "vdgraph_separate_bd.h"
#include "vdgraph_separate_df.h"
#include "vdgraph_separate_ml.h"
#include "vdgraph_separate_sq.h"
#include "vdgraph_separate_st.h"
#include "vdgraph_separate_zr.h"

/*
**  The static and global variables.
*/

static Vdgraph              vdgraphdummy;         /* Dummy distributed separator graph for offset computations */

static union {
  VdgraphSeparateBdParam    param;
  StratNodeMethodData       padding;
} vdgraphseparatedefaultbd = { { 3, &stratdummy } };

static union {
  VdgraphSeparateDfParam    param;
  StratNodeMethodData       padding;
} vdgraphseparatedefaultdf = { { 0, 300, 1.0, 0.0, 0.2 } };

static union {
  VdgraphSeparateMlParam    param;
  StratNodeMethodData       padding;
} vdgraphseparatedefaultml = { { 5, 1000, 2, 10000, 0.8L, &stratdummy, &stratdummy, &stratdummy } };

static union {
  VdgraphSeparateSqParam    param;
  StratNodeMethodData       padding;
} vdgraphseparatedefaultsq = { { &stratdummy } };


static StratMethodTab       vdgraphseparatestmethtab[] = { /* Distributed graph separation methods array */
                             { VDGRAPHSEPASTMETHBD, "b",  vdgraphSeparateBd, &vdgraphseparatedefaultbd },
                             { VDGRAPHSEPASTMETHDF, "d",  vdgraphSeparateDf, &vdgraphseparatedefaultdf },
                             { VDGRAPHSEPASTMETHML, "m",  vdgraphSeparateMl, &vdgraphseparatedefaultml },
                             { VDGRAPHSEPASTMETHSQ, "q",  vdgraphSeparateSq, &vdgraphseparatedefaultsq },
                             { VDGRAPHSEPASTMETHZR, "z",  vdgraphSeparateZr, NULL },
                             { -1,                  NULL, NULL,              NULL } };


static StratParamTab        vdgraphseparatestparatab[] = { /* Distributed graph separation method parameter list */
                              { VDGRAPHSEPASTMETHBD,  STRATPARAMINT,    "width",
                                (byte *) &vdgraphseparatedefaultbd.param,
                                (byte *) &vdgraphseparatedefaultbd.param.distmax,
                                NULL },
                              { VDGRAPHSEPASTMETHBD,  STRATPARAMSTRAT,  "strat",
                                (byte *) &vdgraphseparatedefaultbd.param,
                                (byte *) &vdgraphseparatedefaultbd.param.strat,
                                (void *) &vdgraphseparateststratab },
                              { VDGRAPHSEPASTMETHDF,  STRATPARAMINT,    "part",
                                (byte *) &vdgraphseparatedefaultdf.param,
                                (byte *) &vdgraphseparatedefaultdf.param.partval,
                                NULL },
                              { VDGRAPHSEPASTMETHDF,  STRATPARAMINT,    "pass",
                                (byte *) &vdgraphseparatedefaultdf.param,
                                (byte *) &vdgraphseparatedefaultdf.param.passnbr,
                                NULL },
                              { VDGRAPHSEPASTMETHDF,  STRATPARAMDOUBLE, "bal",
                                (byte *) &vdgraphseparatedefaultdf.param,
                                (byte *) &vdgraphseparatedefaultdf.param.deltval,
                                NULL },
                              { VDGRAPHSEPASTMETHDF,  STRATPARAMDOUBLE, "dif",
                                (byte *) &vdgraphseparatedefaultdf.param,
                                (byte *) &vdgraphseparatedefaultdf.param.cdifval,
                                NULL },
                              { VDGRAPHSEPASTMETHDF,  STRATPARAMDOUBLE, "rem",
                                (byte *) &vdgraphseparatedefaultdf.param,
                                (byte *) &vdgraphseparatedefaultdf.param.cremval,
                                NULL },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "asc",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.stratasc,
                                (void *) &vdgraphseparateststratab },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "low",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.stratlow,
                                (void *) &vdgraphseparateststratab },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMSTRAT,  "seq",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.stratseq,
                                (void *) &vdgraphseparateststratab },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMINT,    "pass",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.passnbr,
                                NULL },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMINT,    "vert",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.coarnbr,
                                NULL },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMINT,    "dvert",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.foldmax,
                                NULL },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMCASE,   "fold",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.foldval,
                                (void *) "nfd" },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMDOUBLE, "rat",
                                (byte *) &vdgraphseparatedefaultml.param,
                                (byte *) &vdgraphseparatedefaultml.param.coarrat,
                                NULL },
                              { VDGRAPHSEPASTMETHML,  STRATPARAMDEPRECATED | STRATPARAMINT, "dlevl", NULL, NULL, NULL }, /* Wait until MUMPS 5.0 */
                              { VDGRAPHSEPASTMETHML,  STRATPARAMDEPRECATED | STRATPARAMINT, "proc",  NULL, NULL, NULL },
                              { VDGRAPHSEPASTMETHSQ,  STRATPARAMSTRAT,  "strat",
                                (byte *) &vdgraphseparatedefaultsq.param,
                                (byte *) &vdgraphseparatedefaultsq.param.strat,
                                (void *) &vgraphseparateststratab },
                              { VDGRAPHSEPASTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        vdgraphseparatestcondtab[] = { /* Distributed graph condition parameter table */
                              { STRATNODECOND,       STRATPARAMINT,    "edge",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.s.edgeglbnbr,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "levl",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.levlnum,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "load",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.s.veloglbsum,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "proc",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.s.procglbnbr,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "rank",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.s.proclocnum,
                                NULL },
                              { STRATNODECOND,       STRATPARAMINT,    "vert",
                                (byte *) &vdgraphdummy,
                                (byte *) &vdgraphdummy.s.vertglbnbr,
                                NULL },
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    vdgraphseparateststratab = { /* Strategy tables for distributed vertex separation methods */
                              vdgraphseparatestmethtab,
                              vdgraphseparatestparatab,
                              vdgraphseparatestcondtab };

/*******************************************/
/*                                         */
/* This is the generic separation routine. */
/*                                         */
/*******************************************/

/* This routine computes the separation of the
** given distributed graph according to the given
** strategy.
** All distributed vertex separation routines must
** be collective, that is, they must all return
** the same success or failure return value on all
** of the processors onto which they are run. Else,
** the behavior of the software is unpredictable.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
vdgraphSeparateSt (
Vdgraph * restrict const      grafptr,            /*+ Distributed separation graph +*/
const Strat * restrict const  strat)              /*+ Separation strategy          +*/
{
  StratTest           val;
  VdgraphStore        savetab[2];                 /* Results of the two strategies */
  Gnum                compglbload2;               /* Saved global separator load   */
  int                 o;
#ifdef SCOTCH_DEBUG_VDGRAPH2
  MPI_Comm            proccommold;                /* Save area for old communicator */
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("vdgraphSeparateSt: invalid type specification for parser variables");
    return     (1);
  }
  if (sizeof (VdgraphSeparateSqParam) > sizeof (StratNodeMethodData)) {
    errorPrint ("vdgraphSeparateSt: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
#ifdef SCOTCH_DEBUG_VDGRAPH1
  if ((strat->tabl != &vdgraphseparateststratab) &&
      (strat       != &stratdummy)) {
    errorPrint ("vdgraphSeparateSt: invalid parameter (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH1 */

  o = 0;
  switch (strat->type) {
    case STRATNODECONCAT :
      o = vdgraphSeparateSt (grafptr, strat->data.concat.strat[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                          */
        o |= vdgraphSeparateSt (grafptr, strat->data.concat.strat[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (strat->data.cond.test, &val, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct                  */
#ifdef SCOTCH_DEBUG_VDGRAPH2
        if ((val.typetest != STRATTESTVAL) ||
            (val.typenode != STRATPARAMLOG)) {
          errorPrint ("vdgraphSeparateSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
        if (val.data.val.vallog == 1)             /* If expression is true                       */
          o = vdgraphSeparateSt (grafptr, strat->data.cond.strat[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                 */
          if (strat->data.cond.strat[1] != NULL)  /* And if there is an else statement           */
            o = vdgraphSeparateSt (grafptr, strat->data.cond.strat[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :                        /* TODO: Can be multithreaded!     */
      if (((vdgraphStoreInit (grafptr, &savetab[0])) != 0) || /* Allocate save areas */
          ((vdgraphStoreInit (grafptr, &savetab[1])) != 0)) {
        errorPrint       ("vdgraphSeparateSt: out of memory");
        vdgraphStoreExit (&savetab[0]);
        return           (1);
      }

      vdgraphStoreSave (grafptr, &savetab[1]);    /* Save initial bipartition                               */
      if (vdgraphSeparateSt (grafptr, strat->data.select.strat[0]) != 0) { /* If first strategy didn't work */
        vdgraphStoreUpdt (grafptr, &savetab[1]);  /* Restore initial bipartition                            */
        vdgraphStoreSave (grafptr, &savetab[0]);  /* Save it as result                                      */
      }
      else {                                      /* First strategy worked       */
        vdgraphStoreSave (grafptr, &savetab[0]);  /* Save its result             */
        vdgraphStoreUpdt (grafptr, &savetab[1]);  /* Restore initial bipartition */
      }
      if (vdgraphSeparateSt (grafptr, strat->data.select.strat[1]) != 0) /* If second strategy didn't work */
        vdgraphStoreUpdt (grafptr, &savetab[1]);  /* Restore initial bipartition as its result             */

      compglbload2 = grafptr->s.veloglbsum - savetab[0].compglbload[0] - savetab[0].compglbload[1]; /* Compute saved separator load */
      if ( (compglbload2 <  grafptr->compglbload[2]) || /* If first strategy is better */
          ((compglbload2 == grafptr->compglbload[2]) &&
           (abs (savetab[0].compglbloaddlt) < abs (grafptr->compglbloaddlt))))
        vdgraphStoreUpdt (grafptr, &savetab[0]);  /* Restore its result */

      vdgraphStoreExit (&savetab[0]);             /* Free both save areas */
      vdgraphStoreExit (&savetab[1]);
      break;
#ifdef SCOTCH_DEBUG_VDGRAPH1
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_VDGRAPH1 */
    default :
#endif /* SCOTCH_DEBUG_VDGRAPH1 */
#ifdef SCOTCH_DEBUG_VDGRAPH2
      proccommold = grafptr->s.proccomm;          /* Create new communicator to isolate method communications */
      MPI_Comm_dup (proccommold, &grafptr->s.proccomm); 
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
      o = strat->tabl->methtab[strat->data.method.meth].func (grafptr, (void *) &strat->data.method.data);
#ifdef SCOTCH_DEBUG_VDGRAPH2
      MPI_Comm_free (&grafptr->s.proccomm);       /* Restore old communicator */
      grafptr->s.proccomm = proccommold;
#endif /* SCOTCH_DEBUG_VDGRAPH2 */
#ifdef SCOTCH_DEBUG_VDGRAPH1
      break;
    default :
      errorPrint ("vdgraphSeparateSt: invalid parameter (2)");
      return     (1);
#endif /* SCOTCH_DEBUG_VDGRAPH1 */
  }
  return (o);
}
