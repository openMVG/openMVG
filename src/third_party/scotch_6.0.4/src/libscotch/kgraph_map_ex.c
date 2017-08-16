/* Copyright 2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_ex.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tries to balance the        **/
/**                subgraphs of the partition as best as   **/
/**                it can.                                 **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 27 may 2011     **/
/**                                 to     21 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_EX

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "kgraph.h"
#include "kgraph_map_ex.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the balanced mapping.
** It returns:
** - 0 : if mapping could be computed.
** - 1 : on error.
*/

int
kgraphMapEx (
Kgraph * const                  grafptr,          /*+ Graph to map      +*/
const KgraphMapExParam * const  paraptr)          /*+ Method parameters +*/
{
  KgraphMapExDom * restrict   doextab;
  KgraphMapExSort * restrict  sorttab;
  KgraphMapExTerm * restrict  termtab;
  KgraphMapExTree * restrict  treetab;
  Anum * restrict             parttax;
  Anum                        treenbr;            /* Number of nodes in tree structure       */
  const Arch * restrict       archptr;
  ArchDom                     domndat;            /* Root domain                             */
  Anum                        domnnbr;
  Anum                        domnnum;
  Gnum                        sortnbr;            /* Number of non-fixed vertices            */
  Gnum                        sortnum;
  Anum                        termnbr;            /* Number of terminal domains in mapping   */
  Gnum                        vertnum;
  Gnum                        vertnnd;
  double                      velosum;            /* Sum of vertex weights                   */
  double                      wghtsum;            /* Sum of architecture weights             */
  Anum                        wghttmp;            /* Sum of architecture weights for archVar */
  int                         flagval;            /* Flag unset if load imbalance to fix     */

  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Anum * restrict const pfixtax = grafptr->pfixtax;

  grafptr->kbalval = paraptr->kbalval;            /* Store last k-way imbalance ratio */
  domnnbr = grafptr->m.domnnbr;
  sortnbr = grafptr->s.vertnbr - grafptr->vfixnbr; /* Only sort non-fixed vertices */
  if (memAllocGroup ((void **) (void *)
                     &doextab, (size_t) (domnnbr * sizeof (KgraphMapExDom)),
                     &sorttab, (size_t) (sortnbr * sizeof (KgraphMapExSort)),
                     &termtab, (size_t) (domnnbr * sizeof (KgraphMapExTerm)),
                     &treetab, (size_t) (domnnbr * sizeof (KgraphMapExTree) * 2), NULL) == NULL) {
    errorPrint ("kgraphMapEx: out of memory");
    return     (1);
  }

  archptr = grafptr->m.archptr;
  archDomFrst (archptr, &domndat);
  wghtsum = (double) archDomWght (archptr, &domndat);
  velosum = (double) grafptr->s.velosum;

  for (domnnum = 0, termnbr = 0, wghttmp = 0, flagval = 1; domnnum < domnnbr; domnnum ++) {
    const ArchDom * restrict  domnptr;

    domnptr = &grafptr->m.domntab[domnnum];
    if (archDomSize (archptr, domnptr) <= 1) {    /* If domain is a terminal (even variable-sized) */
      Anum                termnum;

      wghttmp                     +=              /* Accumulate subdomain loads in case of variable-sized architectures */
      doextab[domnnum].domnwght    = archDomWght (archptr, domnptr);
      doextab[domnnum].compload    = 0;
      doextab[domnnum].comploadmax = ((double) doextab[domnnum].domnwght * velosum * (1.0 + paraptr->kbalval)) / wghtsum;

      termnum = archDomNum (archptr, domnptr);
      termtab[termnbr].termnum = termnum;         /* Record domain in terminal domain array */
      termtab[termnbr].domnnum = domnnum;
      termnbr ++;                                 /* One more terminal domain */

      if ((grafptr->comploadavg[domnnum] + grafptr->comploaddlt[domnnum]) > doextab[domnnum].comploadmax)
        flagval = 0;                              /* Set flag if at least one domain is imbalanced */
    }
  }
  if (archVar (archptr)) {                        /* If architecture is variable-sized */
    Anum                termnum;

    wghtsum = (double) wghttmp / wghtsum;         /* Recompute real load sum */

    for (termnum = 0; termnum < termnbr; termnum ++) {
      Anum                domnnum;

      domnnum = termtab[termnum].domnnum;
      doextab[domnnum].comploadmax = ((double) doextab[domnnum].domnwght * velosum * (1.0 + paraptr->kbalval)) / wghtsum;
      
      if ((grafptr->comploadavg[domnnum] + grafptr->comploaddlt[domnnum]) > doextab[domnnum].comploadmax)
        flagval = 0;                              /* Set flag if at least one domain is imbalanced */
    }
  }

  if (flagval != 0) {                             /* If nothing to do  */
    memFree (doextab);                            /* Free group leader */
    return  (0);
  }

  intSort2asc1 (termtab, termnbr);                /* Sort terminal domains to allow for dichotomy */

  treenbr = 0;                                    /* Prepare to fill tree array; next slot to fill                  */
  kgraphMapExTree (archptr, termtab, termnbr, doextab, treetab, &treenbr, &domndat); /* Recursively fill tree array */

  parttax = grafptr->m.parttax;
  for (vertnum = grafptr->s.baseval, vertnnd = grafptr->s.vertnnd, sortnbr = 0; /* Get vertex weights */
       vertnum < vertnnd; vertnum ++) {
    Gnum                veloval;

    veloval = (velotax != NULL) ? velotax[vertnum] : 1;
    if ((pfixtax == NULL) || (pfixtax[vertnum] < 0)) { /* If vertex is not fixed */
      sorttab[sortnbr].veloval = veloval;         /* Record it for sorting       */
      sorttab[sortnbr].vertnum = vertnum;
      sortnbr ++;
    }
    else
      doextab[parttax[vertnum]].comploadmax -= veloval; /* Reduce available room in domain for non-fixed vertices */
  }
#ifdef SCOTCH_DEBUG_KGRAPH2
  if (sortnbr != (grafptr->s.vertnbr - grafptr->vfixnbr)) {
    errorPrint ("kgraphMapEx: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  if (velotax != NULL)                            /* If vertices are weighted, sort them in ascending order */
    intSort2asc1 (sorttab, sortnbr);

  for (sortnum = sortnbr - 1; sortnum >= 0; sortnum --) { /* For all sorted vertex indices, by descending weights */
    Gnum                  vertnum;
    Gnum                  veloval;
    Anum                  domnnum;

    vertnum = sorttab[sortnum].vertnum;
    veloval = sorttab[sortnum].veloval;
    domnnum = parttax[vertnum];

    if ((doextab[domnnum].compload + veloval) > doextab[domnnum].comploadmax) { /* If leaving vertex in place would cause imbalance */
      domnnum = kgraphMapExFind (archptr, treetab, doextab, domnnum, veloval); /* Try to find better location for vertex load       */
      if (parttax[vertnum] != domnnum) {          /* If vertex moved to another part */
        parttax[vertnum] = domnnum;               /* Set vertex to new part          */
        flagval = 0;                              /* Record change                   */
      }
    }
    doextab[domnnum].compload += veloval;
  }

  memFree (doextab);                              /* Free group leader */

  if (flagval == 0) {                             /* If something changed */
    kgraphFron (grafptr);                         /* Recompute frontier   */
    kgraphCost (grafptr);
  }

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapEx: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}

/* This routine fills the tree structure
** with the relevant node information.
** It returns:
** - void  : in all cases.
*/

static
Anum
kgraphMapExTree (
const Arch * restrict const             archptr,
const KgraphMapExTerm * restrict const  termtab,
const Anum                              termnbr,
KgraphMapExDom * restrict const         doextab,  /*+ Extended domain array, for adding link to tree +*/
KgraphMapExTree * restrict const        treetab,
Anum * restrict const                   treeptr,
const ArchDom * restrict const          domnptr)  /*+ Pointer to subdomain to consider for this node +*/
{
  Anum                treenum;
  int                 o;

  if (archDomSize (archptr, domnptr) > 1) {       /* If not variable-sized architecture and can bipartition */
    ArchDom             domntab[2];               /* Temporary area to store subdomains                     */
    Anum                sonstab[2];
    int                 i, j;

    o = archDomBipart (archptr, domnptr, &domntab[0], &domntab[1]);
#ifdef SCOTCH_DEBUG_KGRAPH2
    if (o != 0) {
      errorPrint ("kgraphMapExTree: internal error");
      return     (-1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    sonstab[0] = kgraphMapExTree (archptr, termtab, termnbr, doextab, treetab, treeptr, &domntab[0]);
    sonstab[1] = kgraphMapExTree (archptr, termtab, termnbr, doextab, treetab, treeptr, &domntab[1]);
    if (sonstab[0] + sonstab[1] < -1)             /* If both sub-branches do not exist      */
      return (-1);                                /* Return that this branch does not exist */

    treenum = (*treeptr) ++;                      /* Reserve slot for node */

    treetab[treenum].fathnum    =                 /* Assume node has no father (yet)   */
    treetab[treenum].sonstab[1] = -1;             /* Assume second node does not exist */
    for (i = j = 0; i < 2; i ++) {                /* For both prospective sons         */
      Anum                sonsnum;

      sonsnum = sonstab[i];
      if (sonsnum == -1)                          /* If this son does not exist, skip it */
        continue;

      treetab[treenum].sonstab[j] = sonsnum;      /* Link son to current node */
      treetab[sonsnum].fathnum = treenum;
      j ++;                                       /* One more son created */
    }
    treetab[treenum].domndat = *domnptr;
  }
  else {                                          /* If domain is terminal */
    Anum                termmin;
    Anum                termmax;
    Anum                termnum;
    Anum                domnnum;

#ifdef SCOTCH_DEBUG_KGRAPH2
    if (archVar (archptr)) {
      errorPrint ("kgraphMapExTree: not implemented");
      return     (-1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

    termnum = archDomNum (archptr, domnptr);      /* Get number of terminal domain */
    for (termmin = 0, termmax = termnbr; (termmax - termmin) > 1; ) {
      Anum                termmed;

      termmed = (termmax + termmin) / 2;
      if (termtab[termmed].termnum <= termnum)
        termmin = termmed;
      else
        termmax = termmed;
    }
    if (termtab[termmin].termnum != termnum)      /* If terminal not found */
      return (-1);                                /* This branch is dead   */

    treenum = (*treeptr) ++;                      /* Reserve slot for terminal            */
    domnnum = termtab[termmin].domnnum;           /* Get domain location                  */
    treetab[treenum].sonstab[0] = -1;             /* Slot is terminal                     */
    treetab[treenum].sonstab[1] = domnnum;        /* Record domain number of tree node    */
    treetab[treenum].domndat = *domnptr;          /* Record domain data                   */
    doextab[domnnum].treenum = treenum;           /* Record entry point in tree structure */
  }

  return (treenum);
}

/* This routine tries to find a destination
** target vertex that creates the least imbalance.
** It returns:
** - 0 : if a suitable terminal vertex has been found.
** - 1 : if process has to be continued..
*/

static
Anum
kgraphMapExFind (
const Arch * restrict const             archptr,  /*+ Target architecture      +*/
const KgraphMapExTree * restrict const  treetab,  /*+ Subdomain tree structure +*/
const KgraphMapExDom * restrict const   doextab,  /*+ Extended domain array    +*/
const Anum                              domnnum,  /*+ Initial domain number    +*/
const Gnum                              veloval)  /*+ Weight of vertex to map  +*/
{
  KgraphMapExFind     bestdat;
  Anum                treenum;
  int                 o;

  bestdat.comploaddlt = (doextab[domnnum].compload + veloval - doextab[domnnum].comploadmax) / doextab[domnnum].domnwght; /* Compute weighted imbalance */
  bestdat.domnnum     = domnnum;

  treenum = doextab[domnnum].treenum;             /* Start from leaf of subdomain tree */
  do {                                            /* Traverse nodes up to the root     */
    Anum                nodenum;                  /* Number of the son we come from    */
    Anum                othrnum;                  /* Number of the other son           */

    nodenum = treenum;                            /* Record position of current node */
    treenum = treetab[treenum].fathnum;           /* Get father node                 */
    if (treenum == -1)                            /* If already reached root of tree */
      break;

    othrnum = treetab[treenum].sonstab[(treetab[treenum].sonstab[0] == nodenum) ? 1 : 0]; /* Don't consider the branch we come from */

    if (othrnum == -1)                            /* If parent node has only one son */
      continue;                                   /* Skip to upper level             */

    o = kgraphMapExFind2 (archptr, treetab, doextab, &bestdat, treenum, othrnum, veloval);
  }
  while (o != 0);                                 /* As long as proper candidate not found */

  return (bestdat.domnnum);                       /* Return best candidate found */
}

/* This routine tries to find a destination
** target vertex that creates the least imbalance.
** It returns:
** - 0 : if a suitable terminal vertex has been found.
** - 1 : if process has to be continued.
*/

static
int
kgraphMapExFind2 (
const Arch * restrict const             archptr,
const KgraphMapExTree * restrict const  treetab,
const KgraphMapExDom * restrict const   doextab,
KgraphMapExFind * restrict const        bestptr,  /*+ Pointer to structure that keeps best terminal found */
const Anum                              treenum,
const Anum                              nodenum,
const Gnum                              veloval)
{
  Anum                son0num;
  Anum                son1num;

  son0num = treetab[nodenum].sonstab[0];
  son1num = treetab[nodenum].sonstab[1];
  if (son0num != -1) {                            /* If node is not a terminal */
    int                 i;
    int                 o;

    if (son1num == -1)                            /* If node has only one son */
      return (kgraphMapExFind2 (archptr, treetab, doextab, bestptr, treenum, son0num, veloval)); /* Process it directly */

    i = (archDomDist (archptr, &treetab[treenum].domndat, &treetab[son0num].domndat) <= /* Get closest subdomain */
         archDomDist (archptr, &treetab[treenum].domndat, &treetab[son1num].domndat)) ? 0 : 1;

    o = kgraphMapExFind2 (archptr, treetab, doextab, bestptr, treenum, treetab[nodenum].sonstab[i], veloval); /* Process closest branch */
    if (o != 0)                                   /* If didn't find suitable terminal in closest branch */
      o = kgraphMapExFind2 (archptr, treetab, doextab, bestptr, treenum, treetab[nodenum].sonstab[i ^ 1], veloval); /* Process farthest one */
    return (o);
  }
  else {                                          /* If current node is terminal */
    Anum                domnnum;
    Gnum                comploaddlt;

    domnnum = son1num;                            /* Second son records domain number */
    comploaddlt = (doextab[domnnum].compload + veloval - doextab[domnnum].comploadmax) / doextab[domnnum].domnwght; /* Compute weighted imbalance */

    if (comploaddlt < bestptr->comploaddlt) {     /* If found vertex that potentially improves balance */
      bestptr->comploaddlt = comploaddlt;
      bestptr->domnnum     = domnnum;
    }

    return ((comploaddlt <= 0) ? 0 : 1);          /* Return immediatly or go on whether found proper terminal to host vertex or not */
  }
}
