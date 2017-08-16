/* Copyright 2004,2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_fm.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module refines a k-way mapping of  **/
/**                the given mapping graph by applying a   **/
/**                Fiduccia-Mattheyses-like gradient       **/
/**                method.                                 **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 03 mar 2011     **/ 
/**                                 to     23 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_FM

#define SCOTCH_TABLE_GAIN

#include "module.h"
#include "common.h"
#include "gain.h"
#include "fibo.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "parser.h"
#include "kgraph.h"
#include "kgraph_map_fm.h"
#include "kgraph_map_st.h"

/*
**  The static variables.
*/

static const Gnum           kgraphloadone  = 1;
static const Gnum           kgraphnotfixed = -1;

/*********************************/
/*                               */
/* Gain table handling routines. */
/*                               */
/*********************************/

/* This routine returns the vertex of best gain
** whose swap will keep the balance correct.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no more vertices available.
*/

#ifdef SCOTCH_TABLE_GAIN

static
KgraphMapFmEdge *
kgraphMapFmTablGet (
KgraphMapFmTabl * restrict const    tablptr,      /*+ Gain table                       +*/
KgraphMapFmVertex * restrict const  vexxtab,      /*+ Extended vertex hash table array +*/
const Gnum * restrict const         comploaddlt,  /*+ Current imbalance array          +*/
const Gnum * restrict const         comploadmax,  /*+ Maximum imbalance array          +*/
Gnum * restrict                     comploaddiff,
Gnum * restrict                     flagval)
{
  GainTabl *                      gaintab;
  KgraphMapFmEdge *               edxxptr;
  KgraphMapFmEdge *               edxxbest;
  Gnum                            gainbest;
  const GainEntr * restrict       tablbest;
  Gnum                            deltbest;

  gaintab  = *tablptr;
  tablbest = gaintab->tend;                       /* Assume no candidate vertex found yet */
  gainbest = GAINMAX;
  edxxbest = NULL;
  deltbest = GNUMMAX;

  for (edxxptr = (KgraphMapFmEdge *) gainTablFrst (gaintab); /* Select candidate edges */
       (edxxptr != NULL) && (edxxptr->gainlink.tabl < tablbest);
       edxxptr = (KgraphMapFmEdge *) gainTablNext (gaintab, &edxxptr->gainlink)) {
    Gnum                                vexxidx;
    Gnum                                veloval;
    Anum                                domnnumold;
    Anum                                domnnumnew;
    Gnum                                deltvalold;
    Gnum                                deltvalnew;
    Gnum                                deltnew;
    Gnum                                commgain;

    vexxidx = edxxptr->vexxidx;
    veloval = vexxtab[vexxidx].veloval;
#ifdef SCOTCH_DEBUG_KGRAPH2
    if (vexxtab[vexxidx].lockptr != NULL) {       /* If vertex is locked */
      errorPrint ("kgraphMapFmTablGet: internal error (1)");
      return     (NULL); 
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    domnnumold = vexxtab[vexxidx].domnnum;
    domnnumnew = edxxptr->domnnum;
    deltvalold = abs (comploaddlt[domnnumold] - veloval);
    deltvalnew = abs (comploaddlt[domnnumnew] + veloval);
    if (((deltvalold > comploadmax[domnnumold]) && (deltvalold >= abs (comploaddlt[domnnumold]))) || /* If vertex does not enforce or improve balance, skip it */
        ((deltvalnew > comploadmax[domnnumnew]) && (deltvalnew >= abs (comploaddlt[domnnumnew])))) {
      if (edxxptr->cmigmask == ~0) {
        edxxptr->cmigmask = 0;
        kgraphMapFmTablDel (&gaintab, edxxptr);
        kgraphMapFmTablAdd (&gaintab, edxxptr);
      }
      continue;
    }

    deltnew = deltvalold - abs (comploaddlt[domnnumold]) + /* Compute difference in imbalance load sum */
              deltvalnew - abs (comploaddlt[domnnumnew]);
    commgain = edxxptr->commgain + edxxptr->cmiggain;
    if ((commgain < gainbest) ||                  /* And if it gives better gain */
        ((commgain == gainbest) &&                /* Or if it gives better load  */
         (deltnew < deltbest))) {
      tablbest = edxxptr->gainlink.tabl;          /* Select it */
      gainbest = commgain;
      edxxbest = edxxptr;
      deltbest = deltnew;
      if ((abs (comploaddlt[domnnumold]) > comploadmax[domnnumold]) ||
          (abs (comploaddlt[domnnumnew]) > comploadmax[domnnumnew]))
        *flagval = 1;
      else
        *flagval = 0;
    }
  }

  if (edxxbest != NULL)
    *comploaddiff += deltbest;

  return (edxxbest);
}

#else /* SCOTCH_TABLE_GAIN */

/* kgraphMapFmCmpFunc(a,b) must return a negative
** number if a is "better" than b. The smaller, the
** better.
*/

static
int
kgraphMapFmCmpFunc (
const FiboNode * const      data0ptr,             /* TRICK: KgraphMapFmLink is FIRST in KgraphMapFmEdge */
const FiboNode * const      data1ptr)
{
  const KgraphMapFmEdge * const  node0ptr = (KgraphMapFmEdge *) data0ptr;
  const KgraphMapFmEdge * const  node1ptr = (KgraphMapFmEdge *) data1ptr;
  Gnum                          node0val;
  Gnum                          node1val;

  node0val = node0ptr->commgain + node0ptr->cmiggain;
  node1val = node1ptr->commgain + node1ptr->cmiggain;

  if (node0val < node1val)
    return (-1);
  if (node0val > node1val)
    return (1);
  return (0);
}

static
KgraphMapFmEdge *
kgraphMapFmTablGet (
KgraphMapFmTabl * restrict const    tablptr,      /*+ Gain table                       +*/
KgraphMapFmVertex * restrict const  vexxtab,      /*+ Extended vertex hash table array +*/
const Gnum * restrict const         comploaddlt,  /*+ Current imbalance array          +*/
const Gnum * restrict const         comploadmax,  /*+ Maximum imbalance array          +*/
Gnum * restrict                     comploaddiff,
Gnum * restrict                     flagval)
{
  FiboNode *                      remoptr;        /* List of removed links        */
  KgraphMapFmEdge *               edxxptr;
  Gnum                            deltnew;
  FiboNode *                      linkptr;        /* Pointer to current gain link */

  edxxptr = NULL;
  remoptr = NULL;

  while ((linkptr = fiboTreeMin (tablptr)) != NULL) { /* Select candidate vertices */
    Gnum                      vexxidx;
    Gnum                      veloval;
    Anum                      domnnumold;
    Anum                      domnnumnew;
    Gnum                      deltvalold;
    Gnum                      deltvalnew;

    edxxptr = (KgraphMapFmEdge *) linkptr;
    vexxidx = edxxptr->vexxidx;
    veloval = vexxtab[vexxidx].veloval;
#ifdef SCOTCH_DEBUG_KGRAPH2
    if (vexxtab[vexxidx].lockptr != NULL) {       /* If vertex is locked */
      errorPrint ("kgraphMapFmTablGet: internal error (2)");
      return     (NULL); 
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

    fiboTreeDel (tablptr, linkptr);               /* Remove vertex link from table         */
    linkptr->linkdat.prevptr = remoptr;           /* Node has been removed but is not kept */
    remoptr = linkptr;                            /* It will be chained back afterwards    */

    domnnumold = vexxtab[vexxidx].domnnum;
    domnnumnew = edxxptr->domnnum;
    deltvalold = abs (comploaddlt[domnnumold] - veloval);
    deltvalnew = abs (comploaddlt[domnnumnew] + veloval);
    if (((deltvalold <= comploadmax[domnnumold]) || (deltvalold < abs (comploaddlt[domnnumold]))) && /* If vertex does enforce or improve balance, keep it */
        ((deltvalnew <= comploadmax[domnnumnew]) || (deltvalnew < abs (comploaddlt[domnnumnew])))) {
      deltnew = deltvalold - abs (comploaddlt[domnnumold]) + /* Compute difference in imbalance load sum */
                deltvalnew - abs (comploaddlt[domnnumnew]);
      if ((abs (comploaddlt[domnnumold]) > comploadmax[domnnumold]) ||
          (abs (comploaddlt[domnnumnew]) > comploadmax[domnnumnew]))
        *flagval = 1;
      break;
    }
    else {
      if (edxxptr->cmigmask == ~0) {
        edxxptr->cmigmask = 0;
        kgraphMapFmTablDel (&gaintab, edxxptr);
        kgraphMapFmTablAdd (&gaintab, edxxptr);
      }
    }
  }

  while (remoptr != NULL) {                       /* Put back all temporarily removed nodes */
    FiboNode *          tempptr;

    tempptr = remoptr;                            /* Get pointer to node */
    remoptr = remoptr->linkdat.prevptr;           /* Find next node      */
    fiboTreeAdd (tablptr, tempptr);               /* Re-link node        */
  }

  if (linkptr == NULL)
    return (NULL);
      
  *comploaddiff += deltnew;
  return (edxxptr);
}

#endif /* SCOTCH_TABLE_GAIN */

/* This routine checks the consistency of
** the hash structures.
** It returns:
** - 0   : in case of success.
** - !0  : in case of error.
*/

#ifdef SCOTCH_DEBUG_KGRAPH3

static
int
kgraphMapFmCheck (
KgraphMapFmTabl * restrict const            tablptr,        /*+ Gain table          +*/
const Kgraph * restrict const               grafptr,
const KgraphMapFmVertex * restrict const    vexxtab,
const KgraphMapFmEdge * const               edxxtab,        /*+ Extended edge array +*/
const Gnum                                  hashmsk,
const Gnum                                  commload,
Gnum * const                                chektab)
{
  Gnum                  vexxidx;
  Gnum                  edxxidx;
  Gnum                  commloadtmp;
  Anum                  domnnum;

  Anum * restrict const       parttax = grafptr->m.parttax;
  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edlotax = grafptr->s.edlotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * restrict const vmlotax = grafptr->r.vmlotax;

  Gnum * const                edlosumtab  = chektab;
  Gnum * const                edgenbrtab  = chektab + grafptr->m.domnnbr;
  Gnum * const                commgaintab = chektab + (grafptr->m.domnnbr * 2);

  commloadtmp = 0;
  for (vexxidx = 0; vexxidx <= hashmsk; vexxidx ++) { /* For all vertex slots */
    Gnum                vertnum;
    Gnum                edgenum;
    Anum                domnorg;
    Anum                domnlst;                  /* Domain of last vertex for which a distance was computed */
    Anum                distlst;                  /* Last distance computed                                  */
    Gnum                commloadloctmp;
    Gnum                edlosum;
    Gnum                edgenbr;

    commloadloctmp = 0;
    vertnum = vexxtab[vexxidx].vertnum;
    if ((vertnum == ~0))                          /* If unallocated    */
      continue;                                   /* Skip to next slot */

    domnorg = vexxtab[vexxidx].domnnum;
    if ((domnorg < 0) || (domnorg >= grafptr->m.domnnbr)) {
      errorPrint ("kgraphMapFmCheck: invalid vertex part value");
      return     (1);
    }
    if (domnorg != parttax[vertnum]) {
      errorPrint ("kgraphMapFmCheck: invalid extended vertex part value");
      return     (1);
    }
    edxxidx = vexxtab[vexxidx].edxxidx;
    if ((edxxidx >= 0) &&
        (edxxtab[edxxidx].vexxidx != vexxidx)) {
      errorPrint ("kgraphMapFmCheck: invalid extended vertex pointer in extended edge");
      return     (1);
    }

    edlosum =
    edgenbr = 0;
    domnlst = -1;                                 /* Invalid domnain to recompute distance                         */
    memSet (edlosumtab, 0, grafptr->m.domnnbr * 3 * sizeof(Gnum)); /* Reset edlosumtab, edgenbrtab and commgaintab */
    for (edgenum = verttax[vertnum];              /* For all neighbors */
         edgenum < vendtax[vertnum]; edgenum ++) {
      Gnum                vertend;
      Anum                domnend;
      Gnum                edloval;

      vertend = edgetax[edgenum];
      domnend = parttax[vertend];
      edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;

      if (domnorg != domnend) {
        Anum                distval;

        distval = (domnend != domnlst) ? archDomDist (grafptr->m.archptr, &grafptr->m.domntab[domnorg], &grafptr->m.domntab[domnend]) : distlst;
        distlst = distval;
        domnlst = domnend;

        edlosumtab[domnend] += edloval;
        edgenbrtab[domnend] ++;

        commloadloctmp += (Gnum) distval * edloval * grafptr->r.crloval;

        for (edxxidx = vexxtab[vexxidx].edxxidx; (edxxidx != -1) && (edxxtab[edxxidx].domnnum != domnend) ;
             edxxidx = edxxtab[edxxidx].edxxidx) ; /* Search if edge slot exists */
        if (edxxidx == -1) {
          Gnum                vertndd;
          Gnum                vexxndd;

          vertndd = edgetax[edgenum];
          for (vexxndd = 0; (vexxtab[vexxndd].vertnum != vertndd) && (vexxndd <= hashmsk); vexxndd ++) ; /* For all vertex slots */
          errorPrint ("kgraphMapFmCheck: no link for migration of vertex %d to domain %d", vexxtab[vexxidx].vertnum, domnend);
          return     (1);
        }
      }
      else {
        edlosum += edloval;
        edgenbr ++;
      }
    }
    if (edlosum != vexxtab[vexxidx].edlosum) {
      errorPrint ("kgraphMapFmCheck: invalid edge load sum for vertex %d", vexxtab[vexxidx].vertnum);
      return     (1);
    }
    if (edgenbr != vexxtab[vexxidx].edgenbr) {
      errorPrint ("kgraphMapFmCheck: invalid edge number for vertex %d", vexxtab[vexxidx].vertnum);
      return     (1);
    }
    commloadtmp += commloadloctmp;

    for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) {
      Gnum                        domncur;
      Gnum                        domnflg;
  
      domnflg = 0;
      domncur = edxxtab[edxxidx].domnnum;
  
      commgaintab[domncur] = -commloadloctmp;
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                      vertend;
        Anum                      domnend;
        Gnum                      edloval;
  
        vertend = edgetax[edgenum];
        domnend = parttax[vertend];
        edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;
  
        if (domnend == domncur) {
          domnflg = 1;
          continue;
        }
      
        edloval *= grafptr->r.crloval;
        commgaintab[domncur] += edloval        /* Add edge contribution to target domain */
                              * archDomDist (&grafptr->a, &grafptr->m.domntab[domncur], &grafptr->m.domntab[domnend]);
      }
      if (domnflg == 0) {
        errorPrint ("kgraphMapFmCheck: extra link for migration of vertex %d to domain %d", vexxtab[vexxidx].vertnum, domncur);
        return     (1);
      }
    }
    if ((vexxtab[vexxidx].domoptr != NULL) &&
        (vexxtab[vexxidx].cmigload != (archDomIncl (grafptr->m.archptr, &grafptr->m.domntab[domnorg], vexxtab[vexxidx].domoptr) == 1) ? 0
                                       : grafptr->r.cmloval * ((vmlotax != NULL) ? vmlotax[vertnum] : 1)
                                       * archDomDist (grafptr->m.archptr, &grafptr->m.domntab[domnorg], vexxtab[vexxidx].domoptr))) {
      errorPrint ("kgraphMapFmCheck: invalid migration communication load for extended vertex");
      return     (1);
    }
    for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) { /* For vertex links */
      if (edlosumtab[edxxtab[edxxidx].domnnum] != edxxtab[edxxidx].edlosum) {
        errorPrint ("kgraphMapFmCheck: invalid extended edge edlosum value");
        return     (1);
      }
      if (edgenbrtab[edxxtab[edxxidx].domnnum] != edxxtab[edxxidx].edgenbr) {
        errorPrint ("kgraphMapFmCheck: invalid extended edge edgenbr value");
        return     (1);
      }
      if (commgaintab[edxxtab[edxxidx].domnnum] != edxxtab[edxxidx].commgain) {
        errorPrint ("kgraphMapFmCheck: invalid extended edge commgain value");
        return     (1);
      }
      if ((vexxtab[vexxidx].domoptr != NULL) &&
          (edxxtab[edxxidx].cmiggain + vexxtab[vexxidx].cmigload != (archDomIncl (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr) == 1) ? 0
                                                                        : grafptr->r.cmloval * ((vmlotax != NULL) ? vmlotax[vertnum] : 1)
                                                                        * archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr))) {
        errorPrint ("kgraphMapFmCheck: invalid migration communication gain for extended edge");
        return     (1);
      }
    }
  }

  if (grafptr->pfixtax != NULL) {                 /* We have fixed vertices                                   */
    Anum                  domnlst;                /* Domnain of last vertex for which a distance was computed */
    Anum                  distlst;                /* Last distance computed                                   */
    Gnum                  vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      if (grafptr->pfixtax[vertnum] != -1) {      /* If vertex is fixed */
        Gnum                domnorg;
        Gnum                edgenum;

        domnorg = parttax[vertnum];
        domnlst = -1;                             /* Invalid part to recompute distance */
        for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
          Gnum                vertend;
          Anum                domnend;
          Gnum                edloval;

          vertend = edgetax[edgenum];
          domnend = parttax[vertend];
          edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;
          edloval *= grafptr->r.crloval;

          if (domnorg != domnend) {
            Anum                distval;
  
            distval = (domnend != domnlst) ? archDomDist (grafptr->m.archptr, &grafptr->m.domntab[domnorg], &grafptr->m.domntab[domnend]) : distlst;
            distlst = distval;
            domnlst = domnend;

            commloadtmp += (Gnum) distval * edloval;
          }
        }
      }
    }
  }

  if (commload != commloadtmp / 2) {
    errorPrint ("kgraphMapFmCheck: invalid communication load");
    return     (1);
  }

  return (0);
}

#endif /* SCOTCH_DEBUG_KGRAPH3 */

/*****************************/
/*                           */
/* These routines handle the */
/* insertion of vertices     */
/* into the gain arrays.     */
/*                           */
/*****************************/

static
int
kgraphMapFmEdgeResize (
KgraphMapFmVertex * restrict const  vexxtab,      /*+ Extended vertex hash table array +*/
Gnum                                vexxidx,
KgraphMapFmEdge * restrict * const  edxxtabptr,
Gnum * restrict const               edxxsizptr,
Gnum const                          edxxnbr,
KgraphMapFmTabl * restrict const    tablptr)
{
  KgraphMapFmEdge * restrict    edxxtmp;
  KgraphMapFmEdge * restrict    edxxtab;
  Gnum                          edxxsiz;

  edxxtab = *edxxtabptr;
  edxxsiz = *edxxsizptr;

  edxxsiz *= 2;                                 /* Compute new array size  */
  *edxxsizptr = edxxsiz;                        /* Propagate back new size */

  if ((edxxtmp = memRealloc (edxxtab, edxxsiz * sizeof (KgraphMapFmEdge))) == NULL) {
    errorPrint ("kgraphMapFmEdgeResize: out of memory");
    return     (1);
  }
  if (edxxtmp != edxxtab) {                     /* If array has been reallocated, re-link all edges */
    Gnum                edxxidx;

    edxxtab = *edxxtabptr = edxxtmp;            /* Point to new array location */

    kgraphMapFmTablFree (tablptr);              /* Free all edges in gain structure */

    for (edxxidx = 0; edxxidx < edxxnbr; edxxidx ++) {
      if ((vexxtab[edxxtab[edxxidx].vexxidx].lockptr == NULL) && /* Vertex is not locked         */
          (edxxtab[edxxidx].vexxidx != vexxidx) && /* Edge of current vertex will be added later */
          (edxxtab[edxxidx].edxxidx != -2))     /* Skip deprecated extended edges                */
        kgraphMapFmTablAdd (tablptr, &edxxtab[edxxidx]);
    }
  }

  return (0);
}

static
int
kgraphMapFmPartAdd2 (
const Kgraph * restrict const       grafptr,
KgraphMapFmVertex * restrict const  vexxtab,      /*+ Extended vertex hash table array   +*/
Gnum                                vexxidx,
KgraphMapFmEdge * restrict * const  edxxtabptr,
Gnum * restrict const               edxxsizptr,
Gnum * restrict const               edxxnbrptr,
Anum                                domnnum,      /*+ Vertex domain                      +*/
Anum                                domnend,      /*+ Domain of the extended edge to add +*/
Gnum                                edloval,      /*+ Load of the edge to domnend        +*/
KgraphMapFmTabl * restrict const    tablptr)
{
  KgraphMapFmEdge * restrict    edxxtab;
  Gnum                          edxxidx;
  Gnum                          edgenum;
  Gnum                          edlosum;
  Gnum                          edgenbr;
  Gnum                          edxxtmp;
  Gnum                          commgain;

  if (*edxxnbrptr >= *edxxsizptr)                 /* If new slot would not fit  */
    kgraphMapFmEdgeResize (vexxtab, -1, edxxtabptr, edxxsizptr, *edxxnbrptr, tablptr); /* No vexxidx because vertex extended edges will be readd later */

  edxxtab = *edxxtabptr;
  edxxidx = (*edxxnbrptr) ++;                     /* Allocate new slot */

  edxxtab[edxxidx].domnnum = domnend;             /* Set extended edge data */
  edxxtab[edxxidx].distval = archDomDist (&grafptr->a, &grafptr->m.domntab[domnnum], &grafptr->m.domntab[domnend]);
  edxxtab[edxxidx].edlosum = edloval;
  edxxtab[edxxidx].edgenbr = 1;
  edxxtab[edxxidx].vexxidx = vexxidx;
  edxxtab[edxxidx].mswpnum = 0;

  commgain = 0;                                   /* Compute commgain */
  for (edxxtmp = vexxtab[vexxidx].edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) {
    commgain += edxxtab[edxxtmp].edlosum *
                (archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxtmp].domnnum], &grafptr->m.domntab[domnend])
                - edxxtab[edxxtmp].distval);
  }
  commgain += (vexxtab[vexxidx].edlosum - edloval) * edxxtab[edxxidx].distval;
  edxxtab[edxxidx].commgain = commgain * grafptr->r.crloval;

  edxxtab[edxxidx].edxxidx  = vexxtab[vexxidx].edxxidx; /* Link edge to vertex  */
  vexxtab[vexxidx].edxxidx  = edxxidx;

  edxxtab[edxxidx].cmiggain = 0;                  /* Compute migration commgain */
  edxxtab[edxxidx].cmigmask = 0;
  if (vexxtab[vexxidx].domoptr != NULL) {
    Gnum                migcoef;                  /* Equal to -migedloval if vertex was mapped in old mapping */

    migcoef = grafptr->r.cmloval * ((grafptr->r.vmlotax != NULL) ? grafptr->r.vmlotax[vexxtab[vexxidx].vertnum] : 1);

    edxxtab[edxxidx].cmiggain = (archDomIncl (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr) == 1) ? 0
                                : migcoef * archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr);
    edxxtab[edxxidx].cmiggain -= vexxtab[vexxidx].cmigload;
    edxxtab[edxxidx].cmigmask = ~0;
  }

  if (vexxtab[vexxidx].lockptr == NULL)           /* If value has to be linked */
    kgraphMapFmTablAdd (tablptr, &edxxtab[edxxidx]);

  return (0);
}

/* This routine adds a vertex to hash table.
*/

static
int
kgraphMapFmPartAdd (
const Kgraph * restrict const               grafptr,
const Gnum                                  vertnum,
const Gnum                                  vexxidx,  /* Hash value for insertion in vexxtab */
KgraphMapFmVertex * restrict const          vexxtab,
KgraphMapFmEdge **                          edxxtabptr,
Gnum * restrict const                       edxxsizptr,
Gnum * restrict const                       edxxnbrptr, 
KgraphMapFmTabl * restrict const            tablptr)
{
  Gnum                          oldvertnum;       /* Number of current vertex */
  KgraphMapFmEdge * restrict    edxxtab;
  Gnum                          edxxidx;
  Gnum                          edgenum;
  Gnum                          edlosum;
  Gnum                          edgenbr;
  Anum                          domnnum;
  Gnum                          commload;         /* Communication load for local domain */

  const Anum * restrict const parttax = grafptr->m.parttax;
  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * restrict const edlotax = grafptr->s.edlotax;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (vexxtab[vexxidx].vertnum != ~0) {
    errorPrint ("kgraphMapFmPartAdd: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  domnnum = parttax[vertnum];

  vexxtab[vexxidx].vertnum = vertnum;
  vexxtab[vexxidx].domnnum = domnnum;
  vexxtab[vexxidx].veloval = (grafptr->s.velotax != NULL) ? grafptr->s.velotax[vertnum] : 1; /* Vertex will be linked since value is positive */
  vexxtab[vexxidx].mswpnum = 0;                   /* Implicitly set slot as used */
  vexxtab[vexxidx].edxxidx = -1;                  /* No target domains yet       */
  vexxtab[vexxidx].lockptr = NULL;                /* No locked yet               */

  oldvertnum = ((grafptr->s.vnumtax != NULL) &&   /* If there is ancestor graph vertex numbers           */
                (grafptr->s.flagval & KGRAPHHASANCHORS) == 0) /* That are not the ones of the band graph */
             ? grafptr->s.vnumtax[vertnum] : vertnum; /* Get vertex number in original graph             */

  if ((grafptr->r.m.parttax != NULL) &&           /* If we are doing a repartitioning                                     */
      (grafptr->r.m.parttax[oldvertnum] != -1))   /* And if vertex was mapped to an old domain                            */
    vexxtab[vexxidx].domoptr = mapDomain (&grafptr->r.m, oldvertnum); /* Domain in which the vertex was previously mapped */
  else
    vexxtab[vexxidx].domoptr = NULL;

  edxxtab = *edxxtabptr;                          /* Compute and link edges */
  edgenbr = 0;

  for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) {
    Gnum                domnend;

    domnend = edxxtab[edxxidx].domnnum;
    edxxtab[edxxidx].edlosum = 0;
    edxxtab[edxxidx].edgenbr = 0;
    edxxtab[edxxidx].distval = archDomDist (&grafptr->a, &grafptr->m.domntab[domnnum], &grafptr->m.domntab[domnend]);
  }

  commload = 0;                                   /* Load associated with vertex edges */
  edlosum  = 0;
  for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
    Gnum                edxxidx;
    Gnum                vertend;
    Anum                domnend;
    Gnum                edloval;

    vertend = edgetax[edgenum];
    domnend = parttax[vertend];

    edloval = (edlotax != NULL) ? edlotax[edgenum] : 1; /* Do not account for crloval yet */

    if (domnend == domnnum) {                     /* If end vertex belongs to same domain */
      edlosum += edloval;                         /* Record local edge load sum           */
      edgenbr ++;                                 /* Record local edge                    */
      continue;                                   /* Skip further processing              */
    }

    for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) { /* Search for edge */
      if (edxxtab[edxxidx].domnnum == domnend)    /* If edge slot found */
        break;
    }

    if (edxxidx == -1) {                          /* If edge slot not found    */
      if (*edxxnbrptr >= *edxxsizptr)             /* If new slot would not fit */
        kgraphMapFmEdgeResize (vexxtab, vexxidx, edxxtabptr, edxxsizptr, *edxxnbrptr, tablptr);

      edxxidx = (*edxxnbrptr) ++;                 /* Allocate new slot */
      edxxtab = *edxxtabptr;                      /* Update edxxtab    */

      edxxtab[edxxidx].commgain = 0;
      edxxtab[edxxidx].cmiggain = 0;
      edxxtab[edxxidx].cmigmask = (grafptr->r.m.parttax != NULL) ? ~0 : 0;
      edxxtab[edxxidx].domnnum  = domnend;
      edxxtab[edxxidx].distval  = archDomDist (&grafptr->a, &grafptr->m.domntab[domnnum], &grafptr->m.domntab[domnend]);
      edxxtab[edxxidx].edlosum  = 0;
      edxxtab[edxxidx].edgenbr  = 0;
      edxxtab[edxxidx].vexxidx  = vexxidx;
      edxxtab[edxxidx].edxxidx  = vexxtab[vexxidx].edxxidx; /* Link edge to vertex */
      vexxtab[vexxidx].edxxidx  = edxxidx;
      edxxtab[edxxidx].mswpnum  = 0;
    }

    commload += edloval * edxxtab[edxxidx].distval;
    edxxtab[edxxidx].edlosum += edloval;
    edxxtab[edxxidx].edgenbr ++;
  }
  commload *= grafptr->r.crloval;                 /* Multiply all local loads by crloval */
  vexxtab[vexxidx].edlosum = edlosum;
  vexxtab[vexxidx].edgenbr = edgenbr;

  for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) {
    Gnum                domncur;
    Gnum                edxxtmp;
    Gnum                commgain;

    domncur = edxxtab[edxxidx].domnnum;

    commgain = 0;
    for (edxxtmp = vexxtab[vexxidx].edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) {
      Anum                domnend;

      if (edxxtmp == edxxidx)
        continue;

      domnend = edxxtab[edxxtmp].domnnum;
      commgain += edxxtab[edxxtmp].edlosum *      /* Add edge contribution to target domain */
                  archDomDist (&grafptr->a, &grafptr->m.domntab[domncur], &grafptr->m.domntab[domnend]);
    }
    commgain += vexxtab[vexxidx].edlosum * edxxtab[edxxidx].distval;
    edxxtab[edxxidx].commgain = commgain * grafptr->r.crloval - commload;
  }

  vexxtab[vexxidx].cmigload = 0;
  if (vexxtab[vexxidx].domoptr != NULL) {
    Gnum                migcoef;                  /* Equal to -migedloval if vertex was mapped in old mapping */

    migcoef = grafptr->r.cmloval * ((grafptr->r.vmlotax != NULL) ? grafptr->r.vmlotax[vertnum] : 1);

    vexxtab[vexxidx].cmigload = (archDomIncl (&grafptr->a, &grafptr->m.domntab[domnnum], vexxtab[vexxidx].domoptr) == 1) ? 0
                                 : migcoef * archDomDist (&grafptr->a, &grafptr->m.domntab[domnnum], vexxtab[vexxidx].domoptr);
    for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) {
      edxxtab[edxxidx].cmiggain = (archDomIncl (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr) == 1) ? 0
                                   : migcoef * archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxidx].domnnum], vexxtab[vexxidx].domoptr);
      edxxtab[edxxidx].cmiggain -= vexxtab[vexxidx].cmigload;
      edxxtab[edxxidx].cmigmask = ~0;
    }
  }

  if (vexxtab[vexxidx].lockptr == NULL) {         /* If value has to be (re)linked */
    for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) /* Insert edges to neighbors in gain arrays */
      kgraphMapFmTablAdd (tablptr, &edxxtab[edxxidx]);
  }

  return (0);
}

/* This routine doubles the size all of the arrays
** involved in handling the hash table and hash
** vertex arrays.
** It returns:
** - 0   : if resizing succeeded.
** - !0  : if out of memory.
*/

static
int
kgraphMapFmResize (
KgraphMapFmVertex * restrict *    vexxtabptr,     /*+ Extended vertex array      +*/
Gnum * restrict const             hashmaxptr,     /*+ Size of vertex array       +*/
Gnum * const                      hashmskptr,     /*+ Pointer to hash table mask +*/
KgraphMapFmSave * restrict        savetab,        /*+ Move array                 +*/
const Gnum                        savenbr,        /*+ Number of moves recorded   +*/
KgraphMapFmTabl * const           tablptr,        /*+ Gain table                 +*/
KgraphMapFmEdge *                 edxxtab,        /*+ Extended edge array        +*/
KgraphMapFmVertex ** const        lockptr)        /*+ Pointer to locked list     +*/
{
  KgraphMapFmVertex * restrict    vexxtab;        /* Extended vertex array                */
  KgraphMapFmSave * restrict      saveold;        /* Pointer to translated old save array */
  Gnum                            savenum;
  Gnum                            hashold;        /* Size of old hash table (half of new) */
  Gnum                            hashsiz;
  Gnum                            hashmax;
  Gnum                            hashmsk;
  Gnum                            vexxidx;

  hashmax = *hashmaxptr << 1;                     /* Compute new sizes */
  hashold = *hashmaxptr << 2;
  hashsiz = *hashmaxptr << 3;
  hashmsk = hashsiz - 1;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (sizeof (KgraphMapFmVertex) < sizeof (KgraphMapFmSave)) { /* Should always be true */
    errorPrint ("kgraphMapFmResize: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  if ((vexxtab = memRealloc (*vexxtabptr, (size_t) hashsiz * sizeof (KgraphMapFmVertex))) == NULL) {
    errorPrint ("kgraphMapFmResize: out of memory");
    return (1);
  }

  for (savenum = savenbr - 1; savenum >= 0; savenum --) { /* Move save array, in reverse order */
    if (savetab[savenum].type == KGRAPHMAPPFMSAVEVEXX)
      savetab[savenum].u.vexxdat.vexxidx = vexxtab[savetab[savenum].u.vexxdat.vexxidx].vertnum; /* Temporarily translate from hash index to number */
    else if ((savetab[savenum].type & KGRAPHMAPPFMSAVELINK) != 0)
      savetab[savenum].u.linkdat.vexxidx = vexxtab[savetab[savenum].u.linkdat.vexxidx].vertnum; /* Temporarily translate from hash index to number */
  }

  *vexxtabptr = vexxtab;
  *hashmaxptr = hashmax;
  *hashmskptr = hashmsk;

  memSet (vexxtab + hashold, ~0, hashold * sizeof (KgraphMapFmVertex)); /* Set new slots to ~0 */
 
  kgraphMapFmTablFree (tablptr);                  /* Reset gain table  */
  *lockptr = (KgraphMapFmVertex *) -1;            /* Rebuild lock list */

  if (vexxtab[0].vertnum != ~0) {                 /* If vertex overflowing may have occured in old hash table */
    Gnum                        hashtmp;
    Gnum                        hashnew;

    for (hashtmp = hashold - 1, hashnew = 0;      /* Temporarily move vertices away from end of old table to prevent overflowing */
         vexxtab[hashtmp].vertnum != ~0; hashtmp --) {
      while (vexxtab[++ hashnew].vertnum != ~0) ; /* Find an empty slot to receive moved vertex */
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (hashnew >= hashtmp) {
        errorPrint ("kgraphMapFmResize: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

      vexxtab[hashnew] = vexxtab[hashtmp];        /* Move vertex from end of table */
      vexxtab[hashtmp].vertnum = ~0;              /* Set old slot as free          */
    }
  }

  for (vexxidx = 0; vexxidx < hashold; vexxidx ++) { /* Re-compute position of vertices in new table */
    Gnum                        vertnum;

    vertnum = vexxtab[vexxidx].vertnum;
    if (vertnum != ~0) {                          /* If hash slot used */
      Gnum                        hashnew;
      Gnum                        edxxtmp;

      for (hashnew = (vertnum * KGRAPHMAPFMHASHPRIME) & hashmsk; ; hashnew = (hashnew + 1) & hashmsk) {
        if (hashnew == vexxidx)                   /* If hash slot is the same      */
          break;                                  /* There is nothing to do        */
        if (vexxtab[hashnew].vertnum == ~0) {     /* If new slot is empty          */
          vexxtab[hashnew] = vexxtab[vexxidx];    /* Copy data to new slot         */
          vexxtab[vexxidx].mswpnum = ~0;          /* TRICK: not tested at creation */
          vexxtab[vexxidx].vertnum = ~0;          /* Make old slot empty           */
          break;
        }
      }
      if ((hashnew > vexxidx) && (hashnew < hashold)) /* If vertex was an overflowed vertex which will be replaced at end of old table */
        continue;                                 /* It will be re-processed again and re-linked once for good at the end of the loop  */

      for (edxxtmp = vexxtab[hashnew].edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) /* Insert edges to neighbors in gain arrays */
        edxxtab[edxxtmp].vexxidx = hashnew;

      edxxtmp = vexxtab[hashnew].edxxidx;         /* vertices lock is done through their first links           */
      if (vexxtab[hashnew].lockptr == NULL) {     /* If vertex was linked, re-link its edges                   */
        for ( ; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) /* Insert edges to neighbors in gain arrays */
          kgraphMapFmTablAdd (tablptr, &edxxtab[edxxtmp]);
      }
      else                                        /* Re-lock vertices                                          */
        kgraphMapFmLock (*lockptr, &vexxtab[hashnew]);
    }
  }

  for (savenum = 0; savenum < savenbr; savenum ++) {
    if (savetab[savenum].type == KGRAPHMAPPFMSAVEVEXX) {
      Gnum                  vertnum;
      Gnum                  vexxidx;
  
      vertnum = savetab[savenum].u.vexxdat.vexxidx; /* Get vertex number */
      for (vexxidx = (vertnum * KGRAPHMAPFMHASHPRIME) & hashmsk; vexxtab[vexxidx].vertnum != vertnum; vexxidx = (vexxidx + 1) & hashmsk) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (vexxtab[vexxidx].vertnum == ~0) {
          errorPrint ("kgraphMapFmResize: internal error (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      }
      savetab[savenum].u.vexxdat.vexxidx = vexxidx; /* Set new hash table index */
    }
    else if ((savetab[savenum].type & KGRAPHMAPPFMSAVELINK) != 0) {
      Gnum                  vertnum;
      Gnum                  vexxidx;
  
      vertnum = savetab[savenum].u.linkdat.vexxidx; /* Get vertex number */
      for (vexxidx = (vertnum * KGRAPHMAPFMHASHPRIME) & hashmsk; vexxtab[vexxidx].vertnum != vertnum; vexxidx = (vexxidx + 1) & hashmsk) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (vexxtab[vexxidx].vertnum == ~0) {
          errorPrint ("kgraphMapFmResize: internal error (4)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      }
      savetab[savenum].u.linkdat.vexxidx = vexxidx; /* Set new hash table index */
    }
  }

  return (0);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the k-way partitioning.
** It returns:
** - 0 : if k-partition could be computed.
** - 1 : on error.
*/

int
kgraphMapFm (
Kgraph * restrict const           grafptr,        /*+ Active graph      +*/
const KgraphMapFmParam * const    paraptr)        /*+ Method parameters +*/
{
  KgraphMapFmSave *               savetab;        /* Pointer to move array                          */
  KgraphMapFmVertex * restrict    vexxtab;        /* Extended vertex hash table array               */
  KgraphMapFmEdge   *             edxxtab;        /* Edge extended array                            */
  Gnum                            edxxnbr;        /* Current number of links in link array          */
  Gnum                            edxxsiz;
  Gnum                            edxxnum;
  Gnum                            edxxcdx;        /* Index for extended edge copy                   */
  Gnum                            edcunbr;        /* Current number of unused extended edge slots   */
  Gnum                            edxunbr;        /* Number of unused extended edge slots           */
  Gnum *                          comploadmax;    /* Array of maximum imbalances                    */
  Gnum *                          comploaddlt;
  Gnum                            commload;
  Gnum                            cmigload;
  KgraphMapFmTabl * restrict      tablptr;        /* Pointer to gain table for easy access          */
  KgraphMapFmTabl                 tabldat;        /* Gain table                                     */
  KgraphMapFmVertex *             lockptr;
  Gnum                            fronnum;
  Gnum                            fronnbr;
  Gnum                            commloadbst;
  Gnum                            cmigloadbst;
  Gnum                            moveflag;       /* Flag set if useful moves made                  */
  Gnum                            edcpflag;       /* Extended edge array compacting flag            */
  Gnum                            comploaddiff;
  Gnum                            flagval;
  Anum                            domnnum;
  Gnum                            vexxidx;
  Gnum                            passnbr;        /* Maximum number of passes to go                 */
  Gnum                            movenbr;        /* Number of uneffective moves done               */
  Gnum                            savenbr;        /* Number of recorded backtrack moves             */
  Gnum                            savesiz;        /* Size of save array                             */
  Gnum                            mswpnum;        /* Current number of recording sweep              */
  Gnum                            vertnum;
  Gnum                            hashsiz;        /* Size of hash table                             */
  Gnum                            hashmsk;        /* Mask for access to hash table                  */
  Gnum                            hashnum;        /* Hash value                                     */
  Gnum                            hashmax;        /* Maximum number of entries in vertex hash table */
  Gnum                            hashnbr;        /* Current number of entries in vertex hash table */
#ifdef SCOTCH_DEBUG_KGRAPH3
  Gnum *                          chektab;        /* Extra memory needed for the check routine      */
#endif /* SCOTCH_DEBUG_KGRAPH3 */

  Anum * restrict const           parttax = grafptr->m.parttax;
  const Gnum * restrict const     verttax = grafptr->s.verttax;
  const Gnum * restrict const     vendtax = grafptr->s.vendtax;
  const Gnum * restrict const     edgetax = grafptr->s.edgetax;
  const Gnum * restrict const     edlotax = grafptr->s.edlotax;
  const Gnum * restrict const     pfixtax = grafptr->pfixtax;

#ifdef SCOTCH_DEBUG_KGRAPH3                       /* Allocation of extra memory needed for the check routine */
  if ((chektab = memAlloc (grafptr->m.domnnbr * 3 * sizeof(Gnum))) == NULL) {
    errorPrint ("kgraphMapFm: out of memory (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH3 */
  tablptr = &tabldat;

  grafptr->kbalval = paraptr->deltval;            /* Store last k-way imbalance ratio */
  kgraphCost (grafptr);  
  grafptr->commload *= grafptr->r.crloval;        /* crloval must be 1 if we are not doing a repartitioning of a no-band graph */
  if (memAllocGroup ((void **) (void *)           /* Allocation and initialization of imbalance arrays                         */
                    &comploadmax, (size_t) (grafptr->m.domnnbr * sizeof (Gnum)),
                    &comploaddlt, (size_t) (grafptr->m.domnnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("kgraphMapFm: out of memory (2)");
    return     (1);
  }
  for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++) {
    comploadmax[domnnum] = (Gnum) ((double) grafptr->comploadavg[domnnum] * paraptr->deltval);
    comploaddlt[domnnum] = grafptr->comploaddlt[domnnum];
  }

  if (grafptr->fronnbr == 0) {                    /* If no current frontier */
    Anum               domnnum;
    
    for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++) {
      if (abs (grafptr->comploaddlt[domnnum]) > comploadmax[domnnum])
        break;
    }
    if (domnnum == grafptr->m.domnnbr) {          /* If balance is correct */
      memFree (comploadmax);                      /* Nothing to do         */
      return  (0);
    }
    else {                                        /* Imbalance must be fought */
      const Strat *         strat;

      strat = stratInit (&kgraphmapststratab, "r{sep=h{pass=10}}"); /* Use a standard algorithm */
      kgraphMapSt (grafptr, strat);               /* Perform mapping */

      if (grafptr->fronnbr == 0) {                /* If new partition has no frontier */
        memFree (comploadmax);
        return  (0);                              /* This algorithm is still useless  */
      }

      if (memReallocGroup ((void *) comploadmax,  /* domnnbr has changed after mapping */
                           &comploadmax, (size_t) (grafptr->m.domnnbr * sizeof (Gnum)),
                           &comploaddlt, (size_t) (grafptr->m.domnnbr * sizeof (Gnum)), NULL) == NULL) {
        errorPrint ("kgraphMapFm: out of memory (3)");
        return     (1);
      }

      for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++) { /* Else update compload{max,dlt} according to new partition */
        comploadmax[domnnum] = (Gnum) ((double) grafptr->comploadavg[domnnum] * paraptr->deltval);
        comploaddlt[domnnum] = grafptr->comploaddlt[domnnum];
      }
    }
  }

#ifdef SCOTCH_DEBUG_KGRAPH2                       /* Allocation of extended vertex hash table and extended edge array */
  hashnbr = 2 * grafptr->fronnbr + 1;             /* Ensure resizing will be performed, for maximum code coverage     */
  savesiz = 2 * grafptr->fronnbr + 1;             /* Ensure resizing will be performed, for maximum code coverage     */
  edxxsiz = 2 * grafptr->fronnbr + 1;             /* Ensure resizing will be performed, for maximum code coverage     */
#else /* SCOTCH_DEBUG_KGRAPH2 */ 
  hashnbr = 4 * (grafptr->fronnbr + paraptr->movenbr + grafptr->s.degrmax);
  savesiz = 4 * (grafptr->fronnbr + paraptr->movenbr + grafptr->s.degrmax) * 2;
  edxxsiz = 4 * (grafptr->fronnbr + paraptr->movenbr + grafptr->s.degrmax) * 4;
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  if (hashnbr > grafptr->s.vertnbr)
    hashnbr = grafptr->s.vertnbr;
  if (edxxsiz > grafptr->s.edgenbr)
    edxxsiz = grafptr->s.edgenbr;

  for (hashsiz = 256; hashsiz < hashnbr; hashsiz <<= 1) ; /* Get upper power of two */
  hashmsk = hashsiz - 1;
  hashmax = hashsiz >> 2;

  if (kgraphMapFmTablInit (tablptr) != 0) {
    errorPrint ("kgraphMapFm: internal error (1)"); /* Unable to do proper initialization */
    kgraphMapFmTablExit (tablptr);
    return (1);  
  }
  else {
    if (((vexxtab = memAlloc ((size_t) hashsiz * sizeof (KgraphMapFmVertex))) == NULL) ||
        ((savetab = memAlloc ((size_t) savesiz * sizeof (KgraphMapFmSave)))   == NULL) ||
        ((edxxtab = memAlloc ((size_t) edxxsiz * sizeof (KgraphMapFmEdge)))   == NULL)) {
      errorPrint ("kgraphMapFm: out of memory (4)");
      kgraphMapFmTablExit (tablptr);
      return (1);
    }
  }
  memSet (vexxtab, ~0, hashsiz * sizeof (KgraphMapFmVertex)); /* Set all vertex numbers to ~0 */
  memSet (edxxtab, ~0, edxxsiz * sizeof (KgraphMapFmEdge));   /* Set all edge numbers to ~0   */

  hashnbr = grafptr->fronnbr;
  while (hashnbr >= hashmax) {
    if (kgraphMapFmResize (&vexxtab, &hashmax, &hashmsk, savetab, 0, tablptr, edxxtab, &lockptr) != 0) {
      errorPrint ("kgraphMapFm: out of memory (5)");
      memFree    (vexxtab);                       /* Free group leader */
      kgraphMapFmTablExit (tablptr);
      return (1);
    }
  }

  edxxnbr = 0;
  for (fronnum = 0; fronnum < hashnbr; fronnum ++) { /* Set initial gains */
    Gnum                        vertnum;

    vertnum = grafptr->frontab[fronnum];
    if ((pfixtax == NULL) || (pfixtax[vertnum] == -1)) { /* Add only not fixed vertices */
      for (hashnum = (vertnum * KGRAPHMAPFMHASHPRIME) & hashmsk; vexxtab[hashnum].vertnum != ~0; hashnum = (hashnum + 1) & hashmsk) ;

      kgraphMapFmPartAdd (grafptr, vertnum, hashnum, vexxtab, &edxxtab, &edxxsiz, &edxxnbr, tablptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
      if (vexxtab[hashnum].edxxidx == -1) {       /* If vertex does not have any neighbor */
        errorPrint ("kgraphMapFm: vertex does not belong to frontier");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    }
  }
#ifdef SCOTCH_DEBUG_KGRAPH2
    if (hashnbr >= hashmax) {                     /* Hash table too small (must not occur) */
      errorPrint ("kgraphMapFm: hash table to small");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  commloadbst = grafptr->commload;                /* Start from initial situation                              */
  cmigloadbst = 0;                                /* Do not take initial migration cost situation into account */

#ifdef SCOTCH_DEBUG_KGRAPH3
  if (kgraphMapFmCheck (tablptr, grafptr, vexxtab, edxxtab, hashmsk, commloadbst, chektab) != 0) {
    errorPrint ("kgraphMapFm: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH3 */

  passnbr = paraptr->passnbr;                     /* Set remaining number of passes    */
  savenbr = 0;                                    /* For empty backtrack of first pass */
  mswpnum = 0;                                    /* Will be incremented afterwards    */
  lockptr = (KgraphMapFmVertex *) -1;             /* Locked list is empty              */
  edxunbr = 0;

  do {                                            /* As long as there are improvements */
    KgraphMapFmEdge *     edxxptr;
    Gnum                  oldsavenbr;

    oldsavenbr = savenbr;

    while (savenbr -- > 0) {                      /* Delete exceeding moves */
      /* Restore last correct state of the graph 
       * All unlocked vertices have all there valid extended edges in the table
       * Any of the deprecated edges are in the table
       * Any of locked vertices have extended edges in the table */
      Gnum                vexxidx;  
      Gnum                oldveloval;
      Anum                domnnum;
      Anum                domnorg;
      Gnum                edxxidx;
      Gnum *              edxiptr;

      switch (savetab[savenbr].type) {
        case KGRAPHMAPPFMSAVEVEXX:
          vexxidx    = savetab[savenbr].u.vexxdat.vexxidx;
          domnnum    = savetab[savenbr].u.vexxdat.domnnum;
          domnorg    = vexxtab[vexxidx].domnnum;
          oldveloval = vexxtab[vexxidx].veloval;
          vexxtab[vexxidx].domnnum  = domnnum;    /* Restore vertex data */
          vexxtab[vexxidx].veloval  = savetab[savenbr].u.vexxdat.veloval;
          vexxtab[vexxidx].cmigload = savetab[savenbr].u.vexxdat.cmigload;
          vexxtab[vexxidx].edlosum  = savetab[savenbr].u.vexxdat.edlosum;
          vexxtab[vexxidx].edgenbr  = savetab[savenbr].u.vexxdat.edgenbr;
          parttax[vexxtab[vexxidx].vertnum] = domnnum;
          comploaddlt[domnorg] -= oldveloval;     /* Update domain load delta */
          comploaddlt[domnnum] += oldveloval;
          break;

        case KGRAPHMAPPFMSAVEEDXX:
          edxxidx                   = savetab[savenbr].u.edxxdat.edxxidx; /* Restore extended edge data */
          edxxtab[edxxidx].domnnum  = savetab[savenbr].u.edxxdat.domnnum;
          edxxtab[edxxidx].commgain = savetab[savenbr].u.edxxdat.commgain;
          edxxtab[edxxidx].cmiggain = savetab[savenbr].u.edxxdat.cmiggain;
          edxxtab[edxxidx].cmigmask = (grafptr->r.m.parttax != NULL) ? ~0 : 0;
          edxxtab[edxxidx].edlosum  = savetab[savenbr].u.edxxdat.edlosum;
          edxxtab[edxxidx].edgenbr  = savetab[savenbr].u.edxxdat.edgenbr;
          edxxtab[edxxidx].distval  = savetab[savenbr].u.edxxdat.distval;
          break;

        case KGRAPHMAPPFMSAVELINKDEL:
          edxxidx = savetab[savenbr].u.linkdat.edxxidx;
          vexxidx = savetab[savenbr].u.linkdat.vexxidx;
          if (edxxtab[edxxidx].vexxidx != vexxidx) /* Restore correct vexxidx after resize */
            edxxtab[edxxidx].vexxidx = vexxidx;
          edxxtab[edxxidx].edxxidx = vexxtab[vexxidx].edxxidx; /* Add it back to vertex list and set it as used */
          vexxtab[vexxidx].edxxidx = edxxidx;
          edxunbr --;                             /* One more used edge slot                                    */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (edxunbr < 0) {
            errorPrint ("kgraphMapFm: internal error (3)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

          if (vexxtab[vexxidx].lockptr == NULL)
            kgraphMapFmTablAdd (tablptr, &edxxtab[edxxidx]); /* Link it back */
          break;

        case KGRAPHMAPPFMSAVELINKADD:
          edxxidx = savetab[savenbr].u.linkdat.edxxidx;
          vexxidx = savetab[savenbr].u.linkdat.vexxidx;
          if (edxxtab[edxxidx].vexxidx != vexxidx) /* Restore correct vexxidx after resize */
            edxxtab[edxxidx].vexxidx = vexxidx;

          if (vexxtab[vexxidx].lockptr == NULL)
            kgraphMapFmTablDel (tablptr, &edxxtab[edxxidx]); /* Unlink it */

          for (edxiptr = &vexxtab[vexxidx].edxxidx; (*edxiptr != edxxidx) && (*edxiptr != -1); edxiptr = &edxxtab[*edxiptr].edxxidx) ;
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (*edxiptr == -1) {                   /* Since it has been added to the list, extended edge must be in it */
            errorPrint ("kgraphMapFm: internal error (4)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          *edxiptr = edxxtab[edxxidx].edxxidx;    /* Remove it from the vertex list */
          edxxtab[edxxidx].edxxidx = -2;          /* Set slot as unused             */
          edxunbr ++;                             /* One more unused edge slot      */
          break;
      }
    }
    edcpflag = 0;                                /* Assume that extended edge array not need to be compacted   */
    if (edxunbr > (edxxnbr / KGRAPHMAPFMEDXXCOMP)) { /* If more than 20% of edxxtab is unused, compact edxxtab */
      edcpflag = 1;
      for (vexxidx = 0; vexxidx <= hashmsk; vexxidx ++) /* Extended edge list must be recomputed */
        vexxtab[vexxidx].edxxidx = -1;
      edxxcdx = 0;
      edcunbr = 0;
      kgraphMapFmTablFree (tablptr);              /* Free all edges in gain structure */
      for (edxxnum = 0; edxxnum < edxxnbr; edxxnum ++) { /* Remove unused slots of edxxtab */
        if (edxxtab[edxxnum].edxxidx == -2) {
          if (edcunbr > 0) {
            memmove (&edxxtab[edxxcdx], &edxxtab[edxxcdx + edcunbr], (edxxnum - edxxcdx - edcunbr) * sizeof (KgraphMapFmEdge)); /* Since there is overlapping, use memmove */
            edxxcdx = edxxnum - edcunbr;
          }
          else
            edxxcdx = edxxnum;
          edcunbr ++;
        }
        else { 
          vexxidx = edxxtab[edxxnum].vexxidx;

          edxxtab[edxxnum].edxxidx  = vexxtab[vexxidx].edxxidx; /* Link edge to vertex  */
          vexxtab[vexxidx].edxxidx  = edxxnum - edcunbr;        /* Set to new index     */
        }
      }
      if ((edcunbr > 0) && (edxxtab[edxxnbr - 1].edxxidx != -2))
        memmove (&edxxtab[edxxcdx], &edxxtab[edxxcdx + edcunbr], (edxxnum - edxxcdx - edcunbr) * sizeof (KgraphMapFmEdge));
      edxxnbr -= edcunbr;
      edxunbr = 0;
    } 
    else {
      while (oldsavenbr -- > 0) {                 /* Must be sure that all parttax is correct before recompute vertices gains */
        if (savetab[oldsavenbr].type == KGRAPHMAPPFMSAVEVEXX) {
          Gnum                vexxidx;
          Gnum                edxxtmp;

          vexxidx = savetab[oldsavenbr].u.vexxdat.vexxidx;
          if (vexxtab[vexxidx].lockptr == NULL) {
            for (edxxtmp = vexxtab[vexxidx].edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) { /* Relink all vertex links */
              kgraphMapFmTablDel (tablptr, &edxxtab[edxxtmp]); 
              kgraphMapFmTablAdd (tablptr, &edxxtab[edxxtmp]); 
            }
          }
        }
      }
    }
    while (lockptr != (KgraphMapFmVertex *) -1) { /* Unlock locked vertices */ 
      KgraphMapFmVertex *           vexxptr;
      Gnum                          edxxtmp;

      vexxptr = lockptr;                          /* Get vertex associated with lock list */
      lockptr = kgraphMapFmLockNext (lockptr);    /* Point to next vertex to unlock       */ 

#ifdef SCOTCH_DEBUG_KGRAPH2
      if (vexxptr->lockptr == NULL) {
        errorPrint ("kgraphMapFm: internal error (5)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      vexxptr->lockptr = NULL;                    /* Set vertex as unlocked */

      if (edcpflag == 0)                          /* It has not been done during edxxtab compacting                    */
        for (edxxtmp = vexxptr->edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) /* Relink all vertex links */
          kgraphMapFmTablAdd (tablptr, &edxxtab[edxxtmp]);
    }
    if (edcpflag == 1)
      for (edxxnum = 0; edxxnum < edxxnbr; edxxnum ++)
        kgraphMapFmTablAdd (tablptr, &edxxtab[edxxnum]); /* Add all used links, all vertices will be unlocked */

    commload = commloadbst;
    cmigload = cmigloadbst;
    mswpnum ++;                                   /* Forget all recorded moves */

#ifdef SCOTCH_DEBUG_KGRAPH3
    if (kgraphMapFmCheck (tablptr, grafptr, vexxtab, edxxtab, hashmsk, commload, chektab) != 0) {
      errorPrint ("kgraphMapFm: internal error (6)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH3 */

    moveflag     = 0;                              /* No useful moves made              */
    movenbr      = 0;                              /* No ineffective moves recorded yet */
    savenbr      = 0;                              /* Back up to beginning of table     */
    comploaddiff = 0;
    flagval      = 0;

    while ((movenbr < paraptr->movenbr) &&
           ((edxxptr = (KgraphMapFmEdge *) kgraphMapFmTablGet (tablptr, vexxtab, comploaddlt, comploadmax, &comploaddiff, &flagval)) != NULL)) { 
      /* Move one vertex */
      Gnum                vexxidx;
      Gnum                edxxtmp;
      Gnum                veloval;
      Gnum                edgenum;
      Anum                domnnum;
      Gnum                edlosum;
      Gnum                edgenbr;
      Anum                domnend;
      Gnum                edxxidx;
      Gnum *              edxpptr;
      Gnum *              vpexptr;
      Gnum                vexdflg;

      vexxidx = edxxptr->vexxidx;                 /* Get relevant information */
      vertnum = vexxtab[vexxidx].vertnum;
      veloval = vexxtab[vexxidx].veloval;
      domnnum = vexxtab[vexxidx].domnnum;
      edlosum = vexxtab[vexxidx].edlosum;
      edgenbr = vexxtab[vexxidx].edgenbr;
      domnend = edxxptr->domnnum;

      /* Save moved vertex information */
      if (vexxtab[vexxidx].mswpnum != mswpnum) {  /* If extended vertex data not yet recorded */ 
        vexxtab[vexxidx].mswpnum            = mswpnum;
        savetab[savenbr].type               = KGRAPHMAPPFMSAVEVEXX; /* Save extended vertex data */
        savetab[savenbr].u.vexxdat.vexxidx  = vexxidx;
        savetab[savenbr].u.vexxdat.veloval  = veloval;
        savetab[savenbr].u.vexxdat.domnnum  = domnnum;
        savetab[savenbr].u.vexxdat.cmigload = vexxtab[vexxidx].cmigload;
        savetab[savenbr].u.vexxdat.edlosum  = edlosum;
        savetab[savenbr].u.vexxdat.edgenbr  = edgenbr;
        savenbr ++;                               /* One more data recorded */
      } 
      for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxxidx = edxxtab[edxxidx].edxxidx) { /* Save vertex links */
        if (edxxtab[edxxidx].mswpnum != mswpnum) { /* If extended edge data not yet recorded */ 
          edxxtab[edxxidx].mswpnum            = mswpnum;
          savetab[savenbr].type               = KGRAPHMAPPFMSAVEEDXX; /* Save extended edge data */
          savetab[savenbr].u.edxxdat.edxxidx  = edxxidx;
          savetab[savenbr].u.edxxdat.domnnum  = edxxtab[edxxidx].domnnum;
          savetab[savenbr].u.edxxdat.commgain = edxxtab[edxxidx].commgain;
          savetab[savenbr].u.edxxdat.cmiggain = edxxtab[edxxidx].cmiggain;
          savetab[savenbr].u.edxxdat.edlosum  = edxxtab[edxxidx].edlosum;
          savetab[savenbr].u.edxxdat.edgenbr  = edxxtab[edxxidx].edgenbr;
          savetab[savenbr].u.edxxdat.distval  = edxxtab[edxxidx].distval;
          savenbr ++;                             /* One more data recorded */
        }
      }
      movenbr ++;                                 /* One more move done */

      commload += edxxptr->commgain; 
      cmigload += edxxptr->cmiggain;

#ifdef SCOTCH_DEBUG_KGRAPH2
      if (vexxtab[vexxidx].lockptr != NULL) {     /* Vertex is locked */
        errorPrint ("kgraphMapFm: internal error (7)");
        return     (1);
      }
      if (vexxtab[vexxidx].edxxidx == -1) {       /* Vertex not in the frontier */
        errorPrint ("kgraphMapFm: internal error (8)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      kgraphMapFmLock (lockptr, &vexxtab[vexxidx]); /* Set part as having changed (lock vertex) */

      for (edxxtmp = vexxtab[vexxidx].edxxidx; edxxtmp != -1; edxxtmp = edxxtab[edxxtmp].edxxidx) /* Unlink all vertex links from gain arrays */
        kgraphMapFmTablDel (tablptr, &edxxtab[edxxtmp]); 
      edxxtmp = vexxtab[vexxidx].edxxidx;         /* Lock vertex through its first link */

      /* Switch information with corresponding extended edge
       * No add or del of extended edges needed. */
      parttax[vertnum] = domnend;                 /* Change vertex part */

      comploaddlt[domnnum] -= veloval;            /* Account for move */
      comploaddlt[domnend] += veloval;

      vexxtab[vexxidx].domnnum   = domnend;       /* Swap edges */
      edxxptr->domnnum           = domnnum;
      vexxtab[vexxidx].edlosum   = edxxptr->edlosum;
      edxxptr->edlosum           = edlosum;
      vexxtab[vexxidx].edgenbr   = edxxptr->edgenbr;
      edxxptr->edgenbr           = edgenbr;
      vexxtab[vexxidx].cmigload += edxxptr->cmiggain;


      edxpptr = &vexxtab[vexxidx].edxxidx;
      for (edxxidx = vexxtab[vexxidx].edxxidx; edxxidx != -1; edxpptr = &edxxtab[edxxidx].edxxidx, edxxidx = edxxtab[edxxidx].edxxidx) { /* Update vertex links */
        Gnum              domncur;

        domncur = edxxtab[edxxidx].domnnum;
        if (domncur == domnnum) {
          vpexptr = edxpptr;
          continue;
        }
        edxxtab[edxxidx].commgain -= edxxptr->commgain;    
        edxxtab[edxxidx].cmiggain -= edxxptr->cmiggain; 
        edxxtab[edxxidx].distval = archDomDist (&grafptr->a, &grafptr->m.domntab[domnend], &grafptr->m.domntab[domncur]);
      }
      edxxptr->commgain = - edxxptr->commgain;
      edxxptr->cmiggain = - edxxptr->cmiggain;

      vexdflg = 0;
      if (edgenbr == 0) {
        Gnum              edxxidx;  
 
        edxxidx = *vpexptr;
        savetab[savenbr].type = KGRAPHMAPPFMSAVELINKDEL; /* Save it */
        savetab[savenbr].u.linkdat.edxxidx = edxxidx;
        savetab[savenbr].u.linkdat.vexxidx = vexxidx;
        savenbr ++;                               /* One more data recorded              */
        *vpexptr = edxxtab[edxxidx].edxxidx;      /* Remove it from extended vertex list */
        edxxtab[edxxidx].edxxidx = -2;            /* Set extended edge slot as unused    */
        edxunbr ++;                               /* One more unused edge slot           */
      }
      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) { /* (Re-)link neighbors */
        /* - Add the vertex in vexxtab if not yet inserted
         * - Del the edge to the vertnum old domain if it was the only vertex
         *   linked to this domain.
         * - Add and edge to the vertnum new domain if it was the first vertex
         *   linked to this domain.
         * - Update commgain of other edges.
         * - Relink extended edges
         */  
        Gnum                edxxend;
        Gnum                vexxend;
        Gnum                edxoidx;              /* Index of extended edge to old domain                    */ 
        Gnum *              edxcptr;              /* Pointer to index of current extended edge               */
        Gnum *              edxoptr;              /* Pointer to index of extended edge to old domain         */
        Gnum                edxnidx;              /* index of extended edge to new domain                    */
        Gnum                vertend;              /* Number of current end neighbor vertex                   */
        Gnum                edxfidx;              /* Index of first extended edge to update                  */ 
        Gnum                edgeend;
        Gnum                edodnbr;
        Gnum                edndnbr;
        Anum                divoval;              /* Distance between current neighbor domain and old domain */
        Anum                divnval;              /* Distance between current neighbor domain and new domain */
        Gnum                edloval;

        vertend = edgetax[edgenum];
        edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;

        if (parttax[edgetax[edgenum]] == domnnum)
          vexdflg = 1;

        if ((pfixtax != NULL) && (pfixtax[vertend] != -1)) /* Do not link fixed vertices */
          continue;
         
        if (savenbr >= (savesiz - (grafptr->m.domnnbr + 4) * 4)) {
          KgraphMapFmSave *               saveptr; /* Pointer to move array */

          while (savenbr >= (savesiz - (grafptr->m.domnnbr + 4) * 4))
            savesiz += savesiz / 2;

          if ((saveptr = memRealloc (savetab, savesiz * sizeof (KgraphMapFmSave))) == NULL) {
            errorPrint ("kgraphMapFm: out of memory (6)");
            memFree    (savetab);                 /* Free group leader */
            return     (1);
          }
          savetab = saveptr;
        }
        if (hashnbr >= hashmax) {                 /* If extended vertex table is already full */
          if (kgraphMapFmResize (&vexxtab, &hashmax, &hashmsk, savetab, savenbr, tablptr, edxxtab, &lockptr) != 0) { 
            errorPrint ("kgraphMapFm: out of memory (7)");
            memFree    (vexxtab);                 /* Free group leader */
            kgraphMapFmTablExit (tablptr);
            return       (1);
          }
        }
        for (vexxend = (vertend * KGRAPHMAPFMHASHPRIME) & hashmsk; /* Search for vertex or first free slot */
             (vexxtab[vexxend].vertnum != vertend) && (vexxtab[vexxend].vertnum != ~0); vexxend = (vexxend + 1) & hashmsk) ;

        if (vexxtab[vexxend].vertnum == ~0) {     /* If neighbor vertex not yet inserted, create it */
          kgraphMapFmPartAdd (grafptr, vertend, vexxend, vexxtab, &edxxtab, &edxxsiz, &edxxnbr, tablptr);
          hashnbr ++;                             /* One more vertex in hash table */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (vexxtab[vexxend].edxxidx == -1) {
            errorPrint ("kgraphMapFm: internal error (9)");
            return     (1);
          }
          if (edxxtab[vexxtab[vexxend].edxxidx].domnnum != domnend) {
            errorPrint ("kgraphMapFm: internal error (10)");
            return     (1);
          }
          if (edxxtab[vexxtab[vexxend].edxxidx].edxxidx != -1) {
            errorPrint ("kgraphMapFm: internal error (11)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          vexxtab[vexxend].mswpnum            = mswpnum;
          savetab[savenbr].type               = KGRAPHMAPPFMSAVEVEXX; /* Save extended vertex data */
          savetab[savenbr].u.vexxdat.vexxidx  = vexxend;
          savetab[savenbr].u.vexxdat.veloval  = vexxtab[vexxend].veloval;
          savetab[savenbr].u.vexxdat.domnnum  = vexxtab[vexxend].domnnum;
          savetab[savenbr].u.vexxdat.cmigload = vexxtab[vexxend].cmigload;
          savetab[savenbr].u.vexxdat.edlosum  = vexxtab[vexxend].edlosum + edloval; /* Save state before vertex move */
          savetab[savenbr].u.vexxdat.edgenbr  = vexxtab[vexxend].edgenbr + 1;
          savenbr ++;                             /* One more data saved                  */
          savetab[savenbr].type = KGRAPHMAPPFMSAVELINKADD; /* Save extended edge creation */
          savetab[savenbr].u.linkdat.edxxidx = vexxtab[vexxend].edxxidx;
          savetab[savenbr].u.linkdat.vexxidx = vexxend;
          savenbr ++;                             /* One more data recorded */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (savenbr > savesiz) {
            errorPrint ("kgraphMapFm: save array error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          continue;
        }

        if (vexxtab[vexxend].mswpnum != mswpnum) { /* If vertex data not yet recorded */
          vexxtab[vexxend].mswpnum            = mswpnum;
          savetab[savenbr].type               = KGRAPHMAPPFMSAVEVEXX; /* Save extended vertex data */
          savetab[savenbr].u.vexxdat.vexxidx  = vexxend; 
          savetab[savenbr].u.vexxdat.veloval  = vexxtab[vexxend].veloval;
          savetab[savenbr].u.vexxdat.domnnum  = vexxtab[vexxend].domnnum;
          savetab[savenbr].u.vexxdat.cmigload = vexxtab[vexxend].cmigload;
          savetab[savenbr].u.vexxdat.edlosum  = vexxtab[vexxend].edlosum;
          savetab[savenbr].u.vexxdat.edgenbr  = vexxtab[vexxend].edgenbr;
          savenbr ++;                             /* One more data saved */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (savenbr > savesiz) {
            errorPrint ("kgraphMapFm: save array error (2)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        }
        edxoidx =
        edxnidx = -1;                             /* Assume there are no extended edges to domnnum or domnend */
        for (edxcptr = &vexxtab[vexxend].edxxidx, edxxidx = vexxtab[vexxend].edxxidx; edxxidx != -1;
             edxcptr = &edxxtab[edxxidx].edxxidx, edxxidx = edxxtab[edxxidx].edxxidx) { /* Loop on domains */
          Gnum            domncur;                /* Save vertex links */

          domncur = edxxtab[edxxidx].domnnum;
          if (edxxtab[edxxidx].mswpnum != mswpnum) { /* If extended edge data not yet recorded */ 
            edxxtab[edxxidx].mswpnum            = mswpnum;
            savetab[savenbr].type               = KGRAPHMAPPFMSAVEEDXX; /* Save extended edge data */
            savetab[savenbr].u.edxxdat.edxxidx  = edxxidx;
            savetab[savenbr].u.edxxdat.domnnum  = domncur;
            savetab[savenbr].u.edxxdat.commgain = edxxtab[edxxidx].commgain;
            savetab[savenbr].u.edxxdat.cmiggain = edxxtab[edxxidx].cmiggain;
            savetab[savenbr].u.edxxdat.edlosum  = edxxtab[edxxidx].edlosum;
            savetab[savenbr].u.edxxdat.edgenbr  = edxxtab[edxxidx].edgenbr;
            savetab[savenbr].u.edxxdat.distval  = edxxtab[edxxidx].distval;
            savenbr ++;                           /* One more data recorded */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (savenbr > savesiz) {
            errorPrint ("kgraphMapFm: save array error (3)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          } 
          if (domncur == domnnum) {
            edxoidx = edxxidx;
            edxoptr = edxcptr;
            edxxtab[edxxidx].edlosum -= edloval;
            edxxtab[edxxidx].edgenbr --;
            divoval = edxxtab[edxxidx].distval;
          }
          else if (domncur == domnend) {
            edxnidx = edxxidx;
            edxxtab[edxxidx].edlosum += edloval;
            edxxtab[edxxidx].edgenbr ++;
            divnval = edxxtab[edxxidx].distval;
          }
        }
#ifdef SCOTCH_DEBUG_KGRAPH2 
        if ((edxoidx == -1) && (vexxtab[vexxend].domnnum != domnnum)) {
          errorPrint ("kgraphMapFm: internal error (12)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

        if (vexxtab[vexxend].domnnum == domnend) {
          vexxtab[vexxend].edlosum += edloval;
          vexxtab[vexxend].edgenbr ++;
          divnval = 0;
        }

        if (vexxtab[vexxend].domnnum == domnnum) { /* Remove edge from neighbor domain */
          vexxtab[vexxend].edlosum -= edloval;
          vexxtab[vexxend].edgenbr --;
          divoval = 0;
        }
        else {
          if (edxxtab[edxoidx].edgenbr == 0) {    /* If it was the last edge in the end domain, save it */
            savetab[savenbr].type = KGRAPHMAPPFMSAVELINKDEL;
            savetab[savenbr].u.linkdat.edxxidx = edxoidx;
            savetab[savenbr].u.linkdat.vexxidx = vexxend;
            savenbr ++;                           /* One more data recorded */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (savenbr > savesiz) {
            errorPrint ("kgraphMapFm: save array error (4)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
            *edxoptr = edxxtab[edxoidx].edxxidx;  /* Remove it from extended vertex list   */
            if (vexxtab[vexxend].lockptr == NULL) /* If edge is of use (vertex not locked) */
              kgraphMapFmTablDel (tablptr, &edxxtab[edxoidx]); /* Remove it                */
            edxxtab[edxoidx].edxxidx = -2;        /* Set extended edge slot as unused      */
            edxunbr ++;                           /* One more unused edge slot             */
          }
        }

        if ((edxnidx == -1) && (vexxtab[vexxend].domnnum != domnend)) { /* If was first vertex linked to this domain, add edge to new domain */
          Gnum        edxxidx;

          kgraphMapFmPartAdd2 (grafptr, vexxtab, vexxend, &edxxtab, &edxxsiz, &edxxnbr, vexxtab[vexxend].domnnum, domnend, edloval, tablptr); /* Add new extended edge */
#ifdef SCOTCH_DEBUG_KGRAPH2 
          for (edxxidx = vexxtab[vexxend].edxxidx; (edxxidx != -1) && (edxxtab[edxxidx].domnnum != domnend); edxxidx = edxxtab[edxxidx].edxxidx) ;
          if (edxxidx == -1) {
            errorPrint ("kgraphMapFm: internal error (13)");
            return     (1);
          }
          if (edxxidx != edxxnbr - 1) {
            errorPrint ("kgraphMapFm: internal error (14)");
            return     (1);
          }
          if (edxxtab[edxxidx].domnnum != domnend) {
            errorPrint ("kgraphMapFm: internal error (15)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          savetab[savenbr].type = KGRAPHMAPPFMSAVELINKADD; /* Save extended edge creation */
          savetab[savenbr].u.linkdat.edxxidx = edxxnbr - 1;
          savetab[savenbr].u.linkdat.vexxidx = vexxend;
          savenbr ++;                             /* One more data recorded */
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (savenbr > savesiz) {
            errorPrint ("kgraphMapFm: save array error (5)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          divnval = edxxtab[edxxnbr - 1].distval;
          edxfidx = edxxtab[edxxnbr - 1].edxxidx; /* Skip update of the newly added extended edge */
        }
        else
          edxfidx = vexxtab[vexxend].edxxidx;

        edloval *= grafptr->r.crloval;
        for (edxxend = edxfidx; edxxend != -1; edxxend = edxxtab[edxxend].edxxidx) /* Update vertex links */
           edxxtab[edxxend].commgain -= edloval * (divnval - archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxend].domnnum], &grafptr->m.domntab[domnend])
                                                 - divoval + archDomDist (&grafptr->a, &grafptr->m.domntab[edxxtab[edxxend].domnnum], &grafptr->m.domntab[domnnum]));
        if (vexxtab[vexxend].lockptr == NULL) { /* If vertex is not locked */
          for (edxxend = edxfidx; edxxend != -1; edxxend = edxxtab[edxxend].edxxidx) { /* Relink its extended edges */
            kgraphMapFmTablDel (tablptr, &edxxtab[edxxend]); /* Remove it and re-link it                            */
            kgraphMapFmTablAdd (tablptr, &edxxtab[edxxend]);
          }
        }
      }
      if (flagval == 1) {                         /* If move improves balance and we do not respect it */
        commloadbst  = commload;                  /* This move was effective */
        cmigloadbst  = cmigload;
        moveflag     = 1;
        movenbr      =
        savenbr      = 0;
        flagval      = 0;
        comploaddiff = 0;
        mswpnum ++;
      }
      else if ((commload + cmigload) < (commloadbst + cmigloadbst)) { /* If move improves the cost */ 
        commloadbst  = commload;                  /* This move was effective                       */
        cmigloadbst  = cmigload;
        moveflag     = 1;
        movenbr      =
        savenbr      = 0;
        flagval      = 0;
        comploaddiff = 0;
        mswpnum ++;
      } 
      else if (((commload + cmigload) == (commloadbst + cmigloadbst)) && (comploaddiff < 0)) { /* If move improves balance and cut does not decrease */
        commloadbst  = commload;                  /* This move was effective */
        cmigloadbst  = cmigload;
        moveflag     = 1;
        movenbr      =
        savenbr      = 0;
        flagval      = 0;
        comploaddiff = 0;
        mswpnum ++;
      }
      else if (((commload + cmigload) == (commloadbst + cmigloadbst)) && (comploaddiff == 0)) {
        commloadbst = commload;                   /* Forget backtracking */
        cmigloadbst = cmigload;
        movenbr     =
        savenbr     = 0;
        flagval     = 0;
        mswpnum ++;
      }
#ifdef SCOTCH_DEBUG_KGRAPH3
      if (kgraphMapFmCheck (tablptr, grafptr, vexxtab, edxxtab, hashmsk, commload, chektab) != 0) {
        errorPrint ("kgraphMapFm: internal error (16)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH3 */
    }
#ifdef SCOTCH_DEBUG_KGRAPH3
    if (kgraphMapFmCheck (tablptr, grafptr, vexxtab, edxxtab, hashmsk, commload, chektab) != 0) {
      errorPrint ("kgraphMapFm: internal error (17)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH3 */
  } while ((moveflag != 0) &&                     /* As long as vertices are moved                          */
           (-- passnbr != 0));                    /* And we are allowed to loop (TRICK for negative values) */

#ifdef SCOTCH_DEBUG_KGRAPH3
  if (kgraphMapFmCheck (tablptr, grafptr, vexxtab, edxxtab, hashmsk, commload, chektab) != 0) {
    errorPrint ("kgraphMapFm: internal error (18)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH3 */

  while (savenbr -- > 0) {                        /* Delete exceeding moves */
    Gnum                vexxidx;
    Anum                domnnum;
    Gnum                veloval;

    if (savetab[savenbr].type == KGRAPHMAPPFMSAVEVEXX) {
      vexxidx = savetab[savenbr].u.vexxdat.vexxidx;
      domnnum = savetab[savenbr].u.vexxdat.domnnum;
      veloval = savetab[savenbr].u.vexxdat.veloval;

      comploaddlt[vexxtab[vexxidx].domnnum] -= veloval;
      comploaddlt[domnnum] += veloval;
      vexxtab[vexxidx].domnnum = domnnum;         /* Restore vertex data */
      parttax[vexxtab[vexxidx].vertnum] = domnnum;
    }
  }

  commload = 0;
  for (vexxidx = fronnbr = 0;                     /* Build new frontier, compute commload, update parttax */
       vexxidx <= hashmsk; vexxidx ++) {          /* hashsiz no longer valid after resizing, so use hashmsk */
    Gnum                vertnum;
    Gnum                edgenum;
    Anum                domnnum;
    Anum                domnlst;                  /* Domain of last vertex for which a distance was computed */
    Anum                distlst;                  /* Last distance computed                                  */
    Gnum                commcut;

    vertnum = vexxtab[vexxidx].vertnum;           /* Get vertex data from slot */
    if (vertnum != ~0) {
      commcut = 0;
      domnnum = parttax[vertnum];
      domnlst = -1;                              /* Invalid domnain to recompute distance */
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (vexxtab[vexxidx].domnnum != parttax[vertnum]) {
        errorPrint ("kgraphMapFm: internal error (19)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      for (edgenum = verttax[vertnum];
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;
        Gnum                domnend;

        vertend = edgetax[edgenum];
        domnend = parttax[vertend];
        if (domnend != domnnum) {
          Anum              distval;
          Gnum              edloval;

          distval = (domnend != domnlst) ? archDomDist (grafptr->m.archptr, &grafptr->m.domntab[domnnum], &grafptr->m.domntab[domnend]) : distlst;
          distlst = distval;
          domnlst = domnend;
          edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;

          commload += (Gnum) distval * edloval;
          commcut   = 1;
        }
      }
      if (commcut != 0)
        grafptr->frontab[fronnbr ++] = vertnum;
    }
  }

  if (grafptr->pfixtax != NULL) {                 /* We have fixed vertices */
    Gnum                  vertnum;

    for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++) {
      if ((grafptr->pfixtax != NULL) && (grafptr->pfixtax[vertnum] != -1)) { /* If vertex is fixed */
        Gnum                edgenum;
        Anum                domnnum;
        Gnum                commcut;
        Anum                domnlst;              /* Domain of last vertex for which a distance was computed */
        Anum                distlst;              /* Last distance computed                                  */

        commcut = 0;
        domnnum = parttax[vertnum];
        domnlst = -1;                             /* Invalid domnain to recompute distance */
        for (edgenum = verttax[vertnum];          /* For all neighbors                     */
             edgenum < vendtax[vertnum]; edgenum ++) {
          Gnum                vertend;
          Gnum                domnend;

          vertend = edgetax[edgenum];
          domnend = parttax[vertend];
          if (domnend != domnnum) {
            Anum              distval;
            Gnum              edloval;

            distval = (domnend != domnlst) ? archDomDist (grafptr->m.archptr, &grafptr->m.domntab[domnnum], &grafptr->m.domntab[domnend]) : distlst;
            distlst = distval;
            domnlst = domnend;
            edloval = (edlotax != NULL) ? edlotax[edgenum] : 1;

            commload += (Gnum) distval * edloval;
            commcut   = 1;
          }
        }
        if (commcut != 0)
          grafptr->frontab[fronnbr ++] = vertnum;
      }
    }
  }
  grafptr->fronnbr  = fronnbr;
  grafptr->commload = commload / 2;

  for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++)  /* Update graph information */
    grafptr->comploaddlt[domnnum] = comploaddlt[domnnum];

#ifdef SCOTCH_DEBUG_KGRAPH3
  memFree (chektab);                              /* Free group leader */
#endif /* SCOTCH_DEBUG_KGRAPH3 */
  memFree (comploadmax);                          /* Free group leader */
  memFree (vexxtab);
  memFree (savetab);
  memFree (edxxtab);
  kgraphMapFmTablExit (tablptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapFm: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}


