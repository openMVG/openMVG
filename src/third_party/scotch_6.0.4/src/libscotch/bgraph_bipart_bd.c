/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_bd.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module builds a band graph around  **/
/**                the frontier in order to decrease       **/
/**                problem size for the strategy to be     **/
/**                applied.                                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 27 nov 2006     **/
/**                                 to   : 23 dec 2007     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 26 mar 2011     **/
/**                # Version 6.0  : from : 07 nov 2011     **/
/**                                 to   : 08 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BGRAPH_BIPART_BD

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_bd.h"
#include "bgraph_bipart_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

int
bgraphBipartBd (
Bgraph * restrict const           orggrafptr,     /*+ Active graph      +*/
const BgraphBipartBdParam * const paraptr)        /*+ Method parameters +*/
{
  Gnum * restrict             queutab;
  Gnum                        queuheadval;
  Gnum                        queutailval;
  Gnum                        distmax;            /* Maximum distance allowed                                     */
  Gnum * restrict             orgindxtax;         /* Based access to index array for original graph               */
  Gnum                        orgfronnum;
  Gnum                        ancfronnum;
  Gnum                        bndfronnum;
  Bgraph                      bndgrafdat;         /* Band graph structure                                         */
  Gnum                        bndvertnbr;         /* Number of regular vertices in band graph (without anchors)   */
  Gnum                        bndvertnnd;
  const Gnum * restrict       bndvnumtax;         /* Band vertex number array, recycling queutab                  */
  Gnum * restrict             bndveextax;         /* External gain array of band graph, if present                */
  Gnum                        bndveexnbr;         /* Number of external array vertices                            */
  Gnum                        bndvelosum;         /* Load of regular vertices in band graph                       */
  Gnum                        bndedlosum;         /* Sum of edge loads                                            */
  Gnum                        bndcompsize1;       /* Number of regular vertices in part 1 of band graph           */
  Gnum                        bndcompload1;       /* Load of regular vertices in part 1                           */
  Gnum                        bndvlvlnum;         /* Index of first band graph vertex to belong to the last layer */
  Gnum                        bndvertnum;
  Gnum                        bndeancnbr;         /* Number of anchor edges                                       */
  Gnum                        bndedgenbr;         /* Upper bound on the number of edges, including anchor edges   */
  Gnum                        bndedgenum;
  Gnum * restrict             bndedgetax;
  Gnum * restrict             bndedlotax;
  Gnum                        bndedgetmp;
  Gnum                        bnddegrmax;
  Gnum                        bndcommgainextn;    /* Sum of all external gains in band graph                      */
  Gnum                        bndcommgainextn1;   /* Sum of external gains accounted for in load, since in part 1 */
  size_t                      bndedlooftval;      /* Offset of edge load array with respect to edge array         */
  const Gnum * restrict const orgverttax = orggrafptr->s.verttax; /* Fast accesses                                */
  const Gnum * restrict const orgvendtax = orggrafptr->s.vendtax;
  const Gnum * restrict const orgvelotax = orggrafptr->s.velotax;
  const Gnum * restrict const orgedgetax = orggrafptr->s.edgetax;
  const Gnum * restrict const orgedlotax = orggrafptr->s.edlotax;
  GraphPart * restrict const  orgparttax = orggrafptr->parttax;
  Gnum * restrict const       orgfrontab = orggrafptr->frontab;

  if (orggrafptr->fronnbr == 0)                   /* If no separator vertices, apply strategy to full (original) graph */
    return (bgraphBipartSt (orggrafptr, paraptr->stratorg));

  distmax = (Gnum) paraptr->distmax;
  if (distmax < 1)                                /* Always at least one layer of vertices around separator */
    distmax = 1;

  if (memAllocGroup ((void **) (void *)
                     &queutab,    (size_t) (orggrafptr->s.vertnbr * sizeof (Gnum)),
                     &orgindxtax, (size_t) (orggrafptr->s.vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("bgraphBipartBd: out of memory (1)");
    return     (1);
  }
  memSet (orgindxtax, ~0, orggrafptr->s.vertnbr * sizeof (Gnum)); /* Initialize index array */
  orgindxtax -= orggrafptr->s.baseval;

  queuheadval = orggrafptr->fronnbr;              /* First layer is vertices in frontier array                     */
  for (orgfronnum = 0, bndvertnnd = orggrafptr->s.baseval; /* Flag vertices belonging to frontier as band vertices */
       orgfronnum < queuheadval; orgfronnum ++) {
    Gnum                orgvertnum;

    orgvertnum = orgfrontab[orgfronnum];
    orgindxtax[orgvertnum] = bndvertnnd ++;
    queutab[orgfronnum] = orgvertnum;             /* Copy frontier array in queue array */
  }

  bndvelosum   = 0;
  bndedgenbr   = 0;                               /* Estimate upper bound on the number of edges */
  bndcompsize1 = 0;
  bndcompload1 = 0;
  queutailval  = 0;
  bndvlvlnum   = 0;                               /* Assume first layer is last layer   */
  while (distmax -- > 0) {                        /* For all passes except the last one */
    Gnum                orgvertnum;
    Gnum                orgdistval;

    bndvlvlnum = queuheadval;                     /* Record start of last layer */
    while (queutailval < bndvlvlnum) {            /* For all vertices in queue  */
      Gnum                orgvertnum;
      Gnum                orgedgenum;
      Gnum                orgpartval;

      orgvertnum = queutab[queutailval ++];
#ifdef SCOTCH_DEBUG_BGRAPH2
      if ((orgvertnum < orggrafptr->s.baseval) || (orgvertnum >= orggrafptr->s.vertnnd)) {
        errorPrint ("bgraphBipartBd: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
      bndedgenbr += orgvendtax[orgvertnum] - orgverttax[orgvertnum]; /* Exact number of edges */
      orgpartval = orgparttax[orgvertnum];
      bndcompsize1 += orgpartval;                 /* Count vertices in part 1 */
      if (orgvelotax != NULL) {
        Gnum                orgveloval;

        orgveloval    = orgvelotax[orgvertnum];
        bndvelosum   += orgveloval;
        bndcompload1 += orgveloval * orgpartval;
      }

      for (orgedgenum = orgverttax[orgvertnum];
           orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
        Gnum                orgvertend;

        orgvertend = orgedgetax[orgedgenum];
        if (orgindxtax[orgvertend] == ~0) {       /* If vertex not visited yet */
          orgindxtax[orgvertend] = bndvertnnd ++; /* Flag it as enqueued       */
          queutab[queuheadval ++] = orgvertend;   /* Enqueue it                */
        }
      }
    }
  }
  bndedgenbr += queuheadval - queutailval;        /* As many edges from anchors as remaining vertices */
  while (queutailval < queuheadval) {             /* Process vertices in last layer                   */
    Gnum                orgvertnum;
    Gnum                orgpartval;

    orgvertnum = queutab[queutailval ++];
#ifdef SCOTCH_DEBUG_BGRAPH2
    if ((orgvertnum < orggrafptr->s.baseval) || (orgvertnum >= orggrafptr->s.vertnnd)) {
      errorPrint ("bgraphBipartBd: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
    bndedgenbr += orgvendtax[orgvertnum] - orgverttax[orgvertnum]; /* Upper bound on number of edges, including anchor edge */
    orgpartval = orgparttax[orgvertnum];
    bndcompsize1 += orgpartval;                   /* Count vertices in part 1 */
    if (orgvelotax != NULL) {
      Gnum                orgveloval;

      orgveloval    = orgvelotax[orgvertnum];
      bndvelosum   += orgveloval;
      bndcompload1 += orgveloval * orgpartval;
    }
  }
  bndvertnbr = bndvertnnd - orggrafptr->s.baseval;
  if (orgvelotax == NULL) {
    bndvelosum   = bndvertnbr;
    bndcompload1 = bndcompsize1;
  }

  if ((bndcompsize1 >= (orggrafptr->s.vertnbr - orggrafptr->compsize0)) || /* If either part has all of its vertices in band, use plain graph instead */
      ((bndvertnbr - bndcompsize1) >= orggrafptr->compsize0)) {
    memFree (queutab);
    return  (bgraphBipartSt (orggrafptr, paraptr->stratorg));
  }                                               /* TRICK: since always at least one missing vertex per part, there is room for anchor vertices */

  queutab[bndvertnbr]     =                       /* Anchor vertices do not have original vertex numbers */
  queutab[bndvertnbr + 1] = -1;

  memSet (&bndgrafdat, 0, sizeof (Bgraph));
  bndgrafdat.s.flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP | BGRAPHHASANCHORS; /* All Bgraph arrays are non-freeable by bgraphExit() */
  bndgrafdat.s.baseval = orggrafptr->s.baseval;
  bndgrafdat.s.vertnbr = bndvertnbr += 2;         /* "+ 2" for anchor vertices */
  bndgrafdat.s.vertnnd = bndvertnnd + 2;

  bndveexnbr = (orggrafptr->veextax != NULL) ? bndvertnbr : 0;
  if (memAllocGroup ((void **) (void *)           /* Do not allocate vnumtax but keep queutab instead */
                     &bndgrafdat.s.verttax, (size_t) ((bndvertnbr + 1) * sizeof (Gnum)),
                     &bndgrafdat.s.velotax, (size_t) (bndvertnbr       * sizeof (Gnum)),
                     &bndveextax,           (size_t) (bndveexnbr       * sizeof (Gnum)),
                     &bndgrafdat.frontab,   (size_t) (bndvertnbr       * sizeof (Gnum)),
                     &bndgrafdat.parttax,   (size_t) (bndvertnbr       * sizeof (GraphPart)), NULL) == NULL) {
    errorPrint ("bgraphBipartBd: out of memory (2)");
    memFree    (queutab);
    return     (1);
  }
  bndgrafdat.parttax   -= orggrafptr->s.baseval;  /* From now on we should free a Bgraph and not a Graph */
  bndgrafdat.s.verttax -= orggrafptr->s.baseval;
  bndgrafdat.s.vendtax  = bndgrafdat.s.verttax + 1; /* Band graph is compact */
  bndgrafdat.s.velotax -= orggrafptr->s.baseval;
  bndgrafdat.s.vnumtax  = queutab - orggrafptr->s.baseval; /* TRICK: re-use queue array as vertex number array since vertices taken in queue order; will not be freed as graph vertex arrays are said to be grouped */
  bndgrafdat.s.velosum  = orggrafptr->s.velosum;
  bndgrafdat.s.velotax[bndvertnnd]     = orggrafptr->compload0 - (bndvelosum - bndcompload1); /* Set loads of anchor vertices */
  bndgrafdat.s.velotax[bndvertnnd + 1] = orggrafptr->s.velosum - orggrafptr->compload0 - bndcompload1;
  if (bndveexnbr != 0) {
    bndveextax -= orggrafptr->s.baseval;
    bndgrafdat.veextax = bndveextax;
  }
  else
    bndveextax = NULL;

  if (memAllocGroup ((void **) (void *)
                     &bndgrafdat.s.edgetax, (size_t) (bndedgenbr * sizeof (Gnum)),
                     &bndgrafdat.s.edlotax, (size_t) (bndedgenbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("bgraphBipartBd: out of memory (3)");
    bgraphExit (&bndgrafdat);
    memFree    (queutab);
    return     (1);
  }
  bndgrafdat.s.edgetax -= orggrafptr->s.baseval;
  bndgrafdat.s.edlotax -= orggrafptr->s.baseval;
  bndedgetax = bndgrafdat.s.edgetax;
  bndedlotax = bndgrafdat.s.edlotax;
  bndvnumtax = bndgrafdat.s.vnumtax;

  for (bndvertnum = bndedgenum = orggrafptr->s.baseval, bnddegrmax = bndedlosum = bndcommgainextn = bndcommgainextn1 = 0;
       bndvertnum < bndvlvlnum; bndvertnum ++) {  /* Fill index array for vertices not belonging to last level */
    Gnum                orgvertnum;
    GraphPart           orgpartval;
    Gnum                orgedgenum;
    Gnum                orgedloval;
    Gnum                bnddegrval;

    orgvertnum = bndvnumtax[bndvertnum];
    orgpartval = orgparttax[orgvertnum];
    bndgrafdat.s.verttax[bndvertnum] = bndedgenum;
    bndgrafdat.s.velotax[bndvertnum] = (orgvelotax != NULL) ? orgvelotax[orgvertnum] : 1;
    bndgrafdat.parttax[bndvertnum] = orgpartval;
    if (bndveextax != NULL) {
      Gnum                orgveexval;

      orgveexval = orggrafptr->veextax[orgvertnum];
      bndveextax[bndvertnum] = orgveexval;
      bndcommgainextn       += orgveexval;
      bndcommgainextn1      += orgveexval * (Gnum) orgpartval;
    }
    orgedloval = 1;                               /* Assume unity edge loads if not present */
    for (orgedgenum = orgverttax[orgvertnum];     /* All edges of first levels are kept     */
         orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++, bndedgenum ++) {

#ifdef SCOTCH_DEBUG_BGRAPH2
      if ((bndedgenum >= (bndedgenbr + orggrafptr->s.baseval)) ||
          (orgindxtax[orgedgetax[orgedgenum]] < 0)) {
        errorPrint ("bgraphBipartBd: internal error (3)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
      if (orgedlotax != NULL)
        orgedloval = orgedlotax[orgedgenum];
      bndedlosum += orgedloval;
      bndedgetax[bndedgenum] = orgindxtax[orgedgetax[orgedgenum]];
      bndedlotax[bndedgenum] = orgedloval;
    }

    bnddegrval = bndedgenum - bndgrafdat.s.verttax[bndvertnum];
    if (bnddegrmax < bnddegrval)
      bnddegrmax = bnddegrval;
  }
  bndeancnbr = 0;
  for ( ; bndvertnum < bndvertnnd; bndvertnum ++) { /* Fill index array for vertices belonging to last level */
    Gnum                orgvertnum;
    Gnum                orgedgenum;
    GraphPart           orgpartval;
    Gnum                bnddegrval;
    Gnum                orgedloval;
    Gnum                ancedloval;               /* Accumulated edge load for anchor edge */

    orgvertnum = bndvnumtax[bndvertnum];
    orgpartval = orgparttax[orgvertnum];
    bndgrafdat.s.verttax[bndvertnum] = bndedgenum;
    bndgrafdat.s.velotax[bndvertnum] = (orgvelotax != NULL) ? orgvelotax[orgvertnum] : 1;
    bndgrafdat.parttax[bndvertnum] = orgpartval;  /* Record part for vertices of last level */
    if (bndveextax != NULL) {
      Gnum                orgveexval;

      orgveexval = orggrafptr->veextax[orgvertnum];
      bndveextax[bndvertnum] = orgveexval;
      bndcommgainextn       += orgveexval;
      bndcommgainextn1      += orgveexval * (Gnum) orgpartval;
    }

    ancedloval = 0;
    orgedloval = 1;                               /* Assume unity edge loads if not present */
    for (orgedgenum = orgverttax[orgvertnum];     /* Keep only band edges                   */
         orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
      Gnum                bndvertend;

#ifdef SCOTCH_DEBUG_BGRAPH2
      if (bndedgenum >= (bndedgenbr + orggrafptr->s.baseval)) {
        errorPrint ("bgraphBipartBd: internal error (4)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
      if (orgedlotax != NULL)
        orgedloval = orgedlotax[orgedgenum];
      bndedlosum += orgedloval;                   /* Normal arcs are accounted for twice; anchor arcs only once */
      bndvertend  = orgindxtax[orgedgetax[orgedgenum]];
      if (bndvertend != ~0) {
        bndedgetax[bndedgenum]    = bndvertend;
        bndedlotax[bndedgenum ++] = orgedloval;
      }
      else
        ancedloval += orgedloval;                 /* Accumulate loads of edges linking to anchor vertex */
    }

    bndedlosum += ancedloval;                     /* Account for anchor edges a second time */
    if (ancedloval > 0) {                         /* If vertex is connected to rest of part */
      bndedlotax[bndedgenum]    = ancedloval;
      bndedgetax[bndedgenum ++] = bndvertnnd + (Gnum) orgpartval; /* Add anchor edge to proper anchor vertex */
      bndeancnbr ++;
    }
    bnddegrval = bndedgenum - bndgrafdat.s.verttax[bndvertnum];
    if (bnddegrmax < bnddegrval)
      bnddegrmax = bnddegrval;
  }
  bndgrafdat.parttax[bndvertnnd]     = 0;         /* Set parts of anchor vertices */
  bndgrafdat.parttax[bndvertnnd + 1] = 1;

  bndgrafdat.s.edlosum = bndedlosum;
  bndgrafdat.s.verttax[bndvertnnd] = bndedgenum;  /* Mark end of regular edge array and start of first anchor edge array */
  bndedgetmp = bndedgenum + bndeancnbr;
#ifdef SCOTCH_DEBUG_BGRAPH2
  if ((bndedgetmp - 1) >= (bndedgenbr + orggrafptr->s.baseval)) {
    errorPrint ("bgraphBipartBd: internal error (5)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
  bndgrafdat.s.edgenbr = bndedgetmp - orggrafptr->s.baseval;
  bndgrafdat.s.verttax[bndvertnnd + 2] = bndedgetmp; /* Mark end of edge array with anchor vertices  */
  for (bndvertnum = bndvlvlnum; bndvertnum < bndvertnnd; bndvertnum ++) { /* Fill anchor edge arrays */
    Gnum                orgvertnum;

    orgvertnum = bndvnumtax[bndvertnum];
    if (bndgrafdat.s.verttax[bndvertnum + 1] > bndgrafdat.s.verttax[bndvertnum]) { /* If vertex is not isolated */
      Gnum                bndedgelst;             /* Number of last edge */
      Gnum                bndvertend;

      bndedgelst = bndgrafdat.s.verttax[bndvertnum + 1] - 1;
      bndvertend = bndedgetax[bndedgelst];        /* Get last neighbor of its edge sub-array */
      if (bndvertend >= bndvertnnd) {             /* If it is an anchor                      */
        Gnum                bndedloval;

        bndedloval  = bndedlotax[bndedgelst];
        bndedlosum += bndedloval;
        if (bndvertend == bndvertnnd) {           /* Add edge from proper anchor */
          bndedgetax[bndedgenum]    = bndvertnum;
          bndedlotax[bndedgenum ++] = bndedloval;
        }
        else {
          bndedgetax[-- bndedgetmp] = bndvertnum;
          bndedlotax[bndedgetmp]    = bndedloval;
        }
      }
    }
  }
  bndgrafdat.s.verttax[bndvertnnd + 1] = bndedgenum; /* Mark end of edge array of first anchor and start of second */
#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bndedgenum != bndedgetmp) {
    errorPrint ("bgraphBipartBd: internal error (6)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  if ((bndedgenum == bndgrafdat.s.verttax[bndvertnnd]) || /* If any of the anchor edges is isolated */
      (bndedgenum == bndgrafdat.s.verttax[bndvertnnd + 2])) {
    bgraphExit (&bndgrafdat);                     /* Free all band graph related data */
    memFree    (queutab);
    return     (bgraphBipartSt (orggrafptr, paraptr->stratorg)); /* Work on original graph */
  }

  if (bnddegrmax < (bndgrafdat.s.verttax[bndvertnnd + 1] - bndgrafdat.s.verttax[bndvertnnd]))
    bnddegrmax = (bndgrafdat.s.verttax[bndvertnnd + 1] - bndgrafdat.s.verttax[bndvertnnd]);
  if (bnddegrmax < (bndgrafdat.s.verttax[bndvertnnd + 2] - bndgrafdat.s.verttax[bndvertnnd + 1]))
    bnddegrmax = (bndgrafdat.s.verttax[bndvertnnd + 2] - bndgrafdat.s.verttax[bndvertnnd + 1]);
  bndgrafdat.s.degrmax = bnddegrmax;

  bndedlooftval = bndedlotax - bndedgetax;
  bndgrafdat.s.edgetax = (Gnum *) memRealloc (bndedgetax + bndgrafdat.s.baseval, (bndedlooftval + bndgrafdat.s.edgenbr) * sizeof (Gnum)) - bndgrafdat.s.baseval;
  bndgrafdat.s.edlotax = bndgrafdat.s.edgetax + bndedlooftval; /* Use old index into old array as new index */
  bndedgetax = bndgrafdat.s.edgetax;
  bndedlotax = bndgrafdat.s.edlotax;

  for (bndfronnum = 0, bndvertnum = orggrafptr->s.baseval; /* Fill band frontier array with first vertex indices as they make the separator */
       bndfronnum < orggrafptr->fronnbr; bndfronnum ++, bndvertnum ++)
    bndgrafdat.frontab[bndfronnum] = bndvertnum;

  if (bndveextax != NULL) {
    Gnum                bndcommloadintn;
    Gnum                bndfronnnd;
    Gnum                bndvertnum;
    Gnum                bndedgenum;
    Gnum                bndedloval;

    bndedloval      = 1;                          /* Assume unity edge weights */
    bndcommloadintn = 0;
    for (bndvertnum = orggrafptr->s.baseval, bndfronnnd = bndvertnum + orggrafptr->fronnbr; /* Compute communication load at frontier */
         bndvertnum < bndfronnnd; bndvertnum ++) {
      Gnum                bndpartval;

      bndpartval = (Gnum) bndgrafdat.parttax[bndvertnum];
      if (bndpartval != 0)                        /* Process only frontier vertices in part 0 */
        continue;

      for (bndedgenum = bndgrafdat.s.verttax[bndvertnum]; bndedgenum < bndgrafdat.s.vendtax[bndvertnum]; bndedgenum ++) {
        Gnum                bndpartend;

        bndpartend = (Gnum) bndgrafdat.parttax[bndedgetax[bndedgenum]];
        bndedloval = bndedlotax[bndedgenum];
        bndcommloadintn += bndedloval * bndpartend;
      }
    }

    bndcommloadintn *= orggrafptr->domndist;
    bndveextax[bndvertnnd + 1] = (orggrafptr->commload - orggrafptr->commloadextn0 - bndcommloadintn) - bndcommgainextn1;
    bndveextax[bndvertnnd]     = (orggrafptr->commload - orggrafptr->commloadextn0 - bndcommloadintn) - bndcommgainextn + bndcommgainextn1 + orggrafptr->commgainextn;
  }

  bndgrafdat.fronnbr       = orggrafptr->fronnbr;
  bndgrafdat.compload0     = orggrafptr->compload0;
  bndgrafdat.compload0min  = orggrafptr->compload0min;
  bndgrafdat.compload0max  = orggrafptr->compload0max;
  bndgrafdat.compload0avg  = orggrafptr->compload0avg;
  bndgrafdat.compload0dlt  = orggrafptr->compload0dlt;
  bndgrafdat.compsize0     = bndvertnbr - bndcompsize1 - 1; /* "- 1" for anchor vertex in part 0 */
  bndgrafdat.commload      = orggrafptr->commload;
  bndgrafdat.commloadextn0 = orggrafptr->commloadextn0;
  bndgrafdat.commgainextn  = orggrafptr->commgainextn;
  bndgrafdat.commgainextn0 = orggrafptr->commgainextn0;
  bndgrafdat.domndist      = orggrafptr->domndist;
  bndgrafdat.domnwght[0]   = orggrafptr->domnwght[0];
  bndgrafdat.domnwght[1]   = orggrafptr->domnwght[1];
  bndgrafdat.vfixload[0]   = orggrafptr->vfixload[0];
  bndgrafdat.vfixload[1]   = orggrafptr->vfixload[1];
  bndgrafdat.bbalval       = orggrafptr->bbalval;
  bndgrafdat.levlnum       = orggrafptr->levlnum;

#ifdef SCOTCH_DEBUG_BGRAPH2
  if ((graphCheck (&bndgrafdat.s) != 0) ||        /* Check band graph consistency */
      (bgraphCheck (&bndgrafdat)  != 0)) {
    errorPrint ("bgraphBipartBd: inconsistent band graph data");
    bgraphExit (&bndgrafdat);
    memFree    (queutab);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  if (bgraphBipartSt (&bndgrafdat, paraptr->stratbnd) != 0) { /* Apply strategy to band graph */
    errorPrint ("bgraphBipartBd: cannot bipartition band graph");
    bgraphExit (&bndgrafdat);
    memFree    (queutab);
    return     (1);
  }
  if (bndgrafdat.parttax[bndvertnnd] ==           /* If band graph was too small and anchors went to the same part, apply strategy on full graph */
      bndgrafdat.parttax[bndvertnnd + 1]) {
    bgraphExit (&bndgrafdat);
    memFree    (queutab);
    return     (bgraphBipartSt (orggrafptr, paraptr->stratorg));
  }

  orggrafptr->compload0    = bndgrafdat.compload0;
  orggrafptr->compload0dlt = bndgrafdat.compload0dlt;
  orggrafptr->commload     = bndgrafdat.commload;
  orggrafptr->commgainextn = bndgrafdat.commgainextn;
  orggrafptr->bbalval      = bndgrafdat.bbalval;

  if (bndgrafdat.parttax[bndvertnnd] != 0) {      /* If anchors swapped parts, swap all parts of original vertices */
    Gnum                orgvertnum;

    orggrafptr->compsize0 = orggrafptr->s.vertnbr - orggrafptr->compsize0 - bndcompsize1 + bndgrafdat.compsize0 - 1; /* "- 1" for anchor 0 */

    for (orgvertnum = orggrafptr->s.baseval; orgvertnum < orggrafptr->s.vertnnd; orgvertnum ++)
      orgparttax[orgvertnum] ^= 1;
  }
  else
    orggrafptr->compsize0 = orggrafptr->compsize0 - (bndvertnbr - bndcompsize1) + bndgrafdat.compsize0 + 1; /* "+ 1" for anchor 0 */

  for (bndvertnum = bndgrafdat.s.baseval; bndvertnum < bndvertnnd; bndvertnum ++) /* Update part array of full graph */
    orgparttax[bndvnumtax[bndvertnum]] = bndgrafdat.parttax[bndvertnum];

  for (bndfronnum = orgfronnum = ancfronnum = 0;  /* Update frontier array of full graph */
       bndfronnum < bndgrafdat.fronnbr; bndfronnum ++) {
    Gnum                bndvertnum;
    Gnum                orgvertnum;

    bndvertnum = bndgrafdat.frontab[bndfronnum];
    orgvertnum = bndvnumtax[bndvertnum];
    if (orgvertnum != -1)                         /* If frontier vertex is not an anchor vertex */
      orgfrontab[orgfronnum ++] = orgvertnum;     /* Record it as original frontier vertex      */
    else
      bndgrafdat.frontab[ancfronnum ++] = bndvertnum; /* Else record it for future processing */
  }

  while (ancfronnum > 0) {                        /* For all recorded frontier anchor vertices     */
    Gnum                bndvertnum;               /* Index of frontier anchor vertex in band graph */
    GraphPart           ancpartval;

    bndvertnum = bndgrafdat.frontab[-- ancfronnum];
    ancpartval = bndgrafdat.parttax[bndvertnum];

    for (bndedgenum = bndgrafdat.s.verttax[bndvertnum];
         bndedgenum < bndgrafdat.s.vendtax[bndvertnum]; bndedgenum ++) {
      Gnum                bndvertend;             /* Index of neighbor of anchor vertex in band graph     */
      Gnum                orgvertnum;             /* Index of neighbor of anchor vertex in original graph */
      Gnum                orgedgenum;

      bndvertend = bndedgetax[bndedgenum];
      if (bndgrafdat.parttax[bndvertend] == ancpartval) /* If neighbor is in same part as anchor, skip to next */
        continue;

      orgvertnum = bndvnumtax[bndvertend];

      for (orgedgenum = orgverttax[orgvertnum];   /* For all neighbors of neighbor */
           orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
        Gnum                orgvertend;

        orgvertend = orgedgetax[orgedgenum];      /* Get end vertex in original graph  */
        if (orgindxtax[orgvertend] == ~0) {       /* If vertex never considered before */
#ifdef SCOTCH_DEBUG_BGRAPH2
          if (orgparttax[orgvertend] != ancpartval) { /* Original vertex should always be in same part as anchor */
            errorPrint ("bgraphBipartBd: internal error (7)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_BGRAPH2 */
          orggrafptr->frontab[orgfronnum ++] = orgvertend; /* Add vertex to frontier array */
          orgindxtax[orgvertend] = 0;             /* Flag vertex as already enqueued       */
        }
      }
    }
  }
  orggrafptr->fronnbr = orgfronnum;

  bgraphExit (&bndgrafdat);                       /* Free band graph structures */
  memFree    (queutab);

#ifdef SCOTCH_DEBUG_BGRAPH2
  if (bgraphCheck (orggrafptr) != 0) {
    errorPrint ("bgraphBipartBd: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BGRAPH2 */

  return (0);
}
