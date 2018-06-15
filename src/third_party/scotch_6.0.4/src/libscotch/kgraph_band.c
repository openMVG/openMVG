/* Copyright 2009-2011,2013-2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_band.c                           **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)		   **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a k-way band       **/
/**                graph from the given frontier           **/
/**                array.                                  **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 jan 2009     **/
/**                                 to   : 01 mar 2015     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  kdgraph_band.c in version 5.2 for     **/
/**                  direct k-way partitioning.            **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_BAND

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "mapping.h"
#include "kgraph.h"
#include "kgraph_band.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a index array of given
** width around the current separator.
** It returns:
** - 0   : if the index array could be computed.
** - !0  : on error.
*/

int
kgraphBand (
Kgraph * restrict const           grafptr,        /*+ Pointer to original graph                  +*/
const Gnum                        distmax,        /*+ Maximum distance from separator vertices   +*/
Kgraph * restrict const           bandgrafptr,    /*+ Pointer to band graph structure to fill    +*/
Gnum * const                      bandvertlvlptr, /*+ Pointer to based start index of last level +*/
Gnum * restrict * restrict const  bandvnumptr)    /*+ Pointer to bandvnumtax                     +*/
{
  Gnum                      bandvertnbr;          /* Number of band vertices (including anchor vertices)                                    */
  Gnum                      bandvertnnd;          /* Based end of band vertex array (without anchor vertices)                               */
  Gnum                      bandvertnum;
  Gnum * restrict           bandverttax;
  Gnum * restrict           bandvelotax;
  Gnum * restrict           bandvmlotax;
  Gnum * restrict           bandvnumtax;          /* Original numbers of vertices in band graph                                             */
#define bandedgetab                 bandcompload  /* TRICK: use delta array to compute edge offsets                                         */
  Gnum * restrict           bandedgetax;
  Gnum * restrict           bandedlotax;
  Gnum * restrict           bandeeextab;
  Gnum                      bandedgenbr;
  Gnum                      bandedgenum;
  Gnum                      banddegrmax;
  Gnum * restrict           bandfrontab;
  Anum * restrict           bandparttax;
  Anum * restrict           bandparotax;
  Gnum                      bandvfixnbr;
  Gnum                      bandvertlvlnum;       /* Index of first band vertex belonging to last level                                     */
  Gnum                      bandvertancadj;       /* Flag set when anchor(s) represent inexistent vertices                                  */
  Gnum * restrict         vnumotbdtax;            /* Original to band graph vertex numbers (~0 if not in band graph, -2 for fixed vertices) */
  Gnum * restrict         bandanlotab;            /* Temporary array to store loads to anchors                                              */
  Gnum                      bandedlonbr;          /* Size of local band edge load array                                                     */
  Gnum                      bandedlosum;
  Gnum * restrict           bandcompload;
  Gnum                      bandcommload;
  Gnum * restrict         compload;               /* Load of parts in original graph                                                        */
  Gnum                    fronnum;
  Anum                    domnnbr;
  Anum                    domnnum;
  Gnum                    veloval;
  Gnum                    vertnum;
  Gnum                    vfixnum;
  Gnum                    vfixflag;
  KgraphBandHash * restrict termhashtab;
  Anum                      termhashmsk;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  const Gnum * restrict const edlotax = grafptr->s.edlotax;
  const Gnum * restrict const frontab = grafptr->frontab;
  const Gnum * restrict const pfixtax = grafptr->pfixtax;
  const Anum * restrict const parttax = grafptr->m.parttax;
  const Anum * restrict const parotax = grafptr->r.m.parttax;
  const Gnum * restrict const vmlotax = grafptr->r.vmlotax;

  if (graphBand (&grafptr->s, grafptr->fronnbr, grafptr->frontab, distmax,
                 &vnumotbdtax, &bandvertlvlnum, &bandvertnbr, &bandedgenbr,
                 pfixtax, &bandvfixnbr) != 0) {   /* Get vertices to keep in band graph */
    errorPrint ("kgraphBand: cannot number graph vertices");
    return     (1);
  }

  if (bandvertlvlptr != NULL)
    *bandvertlvlptr = bandvertlvlnum;

  domnnbr = grafptr->m.domnnbr;
  termhashtab = NULL;
  bandanlotab = NULL;
  bandeeextab = NULL;
  if (pfixtax != NULL) {                          /* Fixed vertices may be neighbors of band graph vertices and may not belong to the band graph */
    Anum                termhashsiz;

    const Arch * restrict const tgtarchptr = grafptr->m.archptr;

    for (termhashsiz = 0, termhashmsk = domnnbr; termhashmsk != 0; termhashsiz ++, termhashmsk >>= 1) ; /* Get upper power of two */
    termhashsiz = 1 << (termhashsiz + 2);         /* Fill hash table at 25% maximum */
    termhashmsk = termhashsiz - 1;

    if (memAllocGroup ((void **) (void *)         /* Allocation and initialization of fixed vertices temporary arrays */
                       &termhashtab, (size_t) (termhashsiz * sizeof (KgraphBandHash)),
                       &bandanlotab, (size_t) (domnnbr     * sizeof (Gnum)),
                       &bandeeextab, (size_t) (domnnbr     * sizeof (Gnum)), NULL) == NULL) {
      errorPrint   ("kgraphBand: out of memory (1)");
      return       (1);
    }
    memSet (termhashtab, ~0, termhashsiz * sizeof (KgraphBandHash));
    memSet (bandanlotab,  0, domnnbr     * sizeof (Gnum)); /* Assume there are no extra loads to anchors */
    memSet (bandeeextab,  0, domnnbr     * sizeof (Gnum));
    for (domnnum = 0; domnnum < domnnbr; domnnum ++) {
      ArchDom *                 domnptr;

      domnptr = &grafptr->m.domntab[domnnum];
      if (archDomSize (tgtarchptr, domnptr) == 1) { /* If domain is terminal */
        Gnum                termhashnum;
        Anum                termnum;

        termnum = archDomNum (tgtarchptr, domnptr); /* Get terminal domain number */

        for (termhashnum = (termnum * KGRAPHBANDHASHPRIME) & termhashmsk; ; termhashnum = (termhashnum + 1) & termhashmsk) {
          if (termhashtab[termhashnum].termnum == ~0) { /* If hash slot empty */
            termhashtab[termhashnum].termnum = termnum; /* Create slot        */
            termhashtab[termhashnum].domnnum = domnnum;
            break;
          }
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (termhashtab[termhashnum].termnum == termnum) { /* If hash slot found */
            errorPrint ("kgraphBand: duplicate terminal domain in domain array");
            memFree    (termhashtab);
            kgraphExit (bandgrafptr);
            return     (1);
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        }
      }
    }
  }

  bandedgenbr += 2 * (bandvertnbr + grafptr->s.baseval - bandvertlvlnum) + /* Add edges to and from anchors */
                 grafptr->s.degrmax * grafptr->vfixnbr; /* A band graph vertex that is the neighbour of a fixed vertex will get an extra edge, even if the fixed vertex is not in the band graph */
  bandvertnbr += domnnbr;                         /* Add anchor vertices */
  bandedlonbr  = ((edlotax != NULL) || (pfixtax != NULL)) ? bandedgenbr : 0;

  graphInit (&bandgrafptr->s);
  bandgrafptr->s.flagval   = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP | /* Arrays created by the routine itself */
                             KGRAPHFREEFRON | KGRAPHFREECOMP | KGRAPHHASANCHORS;
  bandgrafptr->s.baseval   = grafptr->s.baseval;
  bandgrafptr->s.vertnbr   = bandvertnbr;
  bandgrafptr->s.vertnnd   = bandvertnbr + bandgrafptr->s.baseval; /* With anchor vertices */
  bandgrafptr->a           = grafptr->a;
  mapInit2 (&bandgrafptr->m,   &bandgrafptr->s, &bandgrafptr->a, &grafptr->m.domnorg, grafptr->m.domnmax,   grafptr->m.domnnbr);
  mapInit2 (&bandgrafptr->r.m, &bandgrafptr->s, &bandgrafptr->a, &grafptr->m.domnorg, grafptr->r.m.domnmax, grafptr->r.m.domnnbr);
  bandgrafptr->r.m.domntab = grafptr->r.m.domntab; /* Band old mapping domain array is a clone of old mapping (no freeing) */
  bandgrafptr->r.m.domnnbr = grafptr->r.m.domnnbr;
  bandgrafptr->r.crloval   = grafptr->r.crloval;
  bandgrafptr->r.cmloval   = grafptr->r.cmloval;
  bandgrafptr->r.vmlotax   = NULL;
  bandgrafptr->vfixnbr     = 0;                   /* Band graphs do not have fixed vertices */
  bandgrafptr->pfixtax     = NULL;
  bandgrafptr->frontab     = NULL;                /* Frontier array not yet allocated          */
  bandgrafptr->comploadavg = NULL;                /* Computation load arrays not yet allocated */
  bandgrafptr->comploaddlt = NULL;
  bandgrafptr->commload    = grafptr->commload;   /* Communication load is preserved */
  bandgrafptr->kbalval     = grafptr->kbalval;
  bandgrafptr->levlnum     = grafptr->levlnum;

  if (memAllocGroup ((void **) (void *)           /* Allocate graph data */
                     &bandgrafptr->s.verttax, (size_t) ((bandvertnbr + 1) * sizeof (Gnum)), /* Compact vertex array */
                     &bandgrafptr->s.velotax, (size_t) (bandvertnbr       * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("kgraphBand: out of memory (2)");
    return     (1);
  }
  if ((bandvnumtax = memAlloc ((bandvertnbr) * sizeof (Gnum))) == NULL) { /* Allocate alone since it is an output */
    errorPrint ("kgraphBand: out of memory (3)");
    return     (1);
  }
  if (vmlotax != NULL) {
    if ((bandvmlotax = memAlloc (bandvertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("kgraphBand: out of memory (4)");
      return     (1);
    }
    memSet (bandvnumtax + bandvertnbr - domnnbr, ~0, domnnbr * sizeof (Gnum)); /* Prevent Valgrind from yelling when centralizing band graphs */
    bandvmlotax -= bandgrafptr->s.baseval;
    bandgrafptr->r.vmlotax  = bandvmlotax;
    bandgrafptr->s.flagval |= KGRAPHFREEVMLO;
  }
  if (parotax != NULL) {
    if ((bandparotax = memAlloc (bandvertnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("kgraphBand: out of memory (5)");
      return     (1);
    }
    memSet (bandparotax + bandvertnbr - bandgrafptr->r.m.domnnbr, ~0, bandgrafptr->r.m.domnnbr * sizeof (Gnum)); /* Prevent Valgrind from yelling when centralizing band graphs */
    bandparotax -= bandgrafptr->s.baseval;
    bandgrafptr->r.m.parttax = bandparotax;
    bandgrafptr->r.m.flagval |= MAPPINGFREEPART;
  }

  bandgrafptr->s.verttax -= bandgrafptr->s.baseval;
  bandvnumtax            -= bandgrafptr->s.baseval;
  bandgrafptr->s.velotax -= bandgrafptr->s.baseval;

  if ((bandgrafptr->s.edgetax = memAlloc ((bandedgenbr + bandedlonbr) * sizeof (Gnum))) == NULL) {
    errorPrint ("kgraphBand: out of memory (6)");
    return     (1);
  }
  bandedlotax             = NULL;
  bandedgetax             =
  bandgrafptr->s.edgetax -= bandgrafptr->s.baseval;
  if ((edlotax != NULL) || (pfixtax != NULL)) {
    bandgrafptr->s.edlotax =
    bandedlotax            = bandedgetax + bandedgenbr;
  }

  if (((bandgrafptr->frontab = memAlloc (bandvertnbr * sizeof (Gnum))) == NULL) || /* Allocation and initialization of imbalance arrays */
      (memAllocGroup ((void **) (void *)
                      &bandgrafptr->comploadavg, (size_t) ((domnnbr + 2) * sizeof (Gnum)), /* TRICK: always keep two slots for collective communication */
                      &bandgrafptr->comploaddlt, (size_t) ((domnnbr + 2) * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("kgraphBand: out of memory (7)");
    if (bandgrafptr->frontab != NULL)
      memFree (bandgrafptr->frontab);
    return (1);
  }
  bandfrontab = bandgrafptr->frontab;

  if ((bandparttax = memAlloc (bandvertnbr * sizeof (Anum))) == NULL) {
    errorPrint ("kgraphBand: out of memory (8)");
    return     (1);
  }
  bandgrafptr->m.parttax =
  bandparttax           -= bandgrafptr->s.baseval;

  if ((bandgrafptr->m.domntab = memAlloc (domnnbr * sizeof (ArchDom))) == NULL) {
    errorPrint ("kgraphBand: out of memory (9)");
    return     (1);
  }
  bandgrafptr->m.flagval |= MAPPINGFREEDOMN;
#ifdef SCOTCH_DEBUG_KGRAPH2
  memSet (bandvnumtax + bandgrafptr->s.baseval, ~0, (bandvertnbr * sizeof (Gnum)));
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  vfixnum = 0;
  for (fronnum = 0, bandvertnum = bandgrafptr->s.baseval;
       fronnum < grafptr->fronnbr; fronnum ++) {  /* Turn all graph frontier vertices into band frontier vertices */
    Gnum              vertnum;

    vertnum = frontab[fronnum];
    if ((pfixtax != NULL) && (pfixtax[vertnum] != -1)) /* It is a fixed vertex */
      vfixnum ++;
    else {
      bandfrontab[bandvertnum - bandgrafptr->s.baseval] = bandvertnum; /* All frontier vertices are first vertices of band graph */
      bandvnumtax[bandvertnum] = vertnum;
      bandvertnum ++;
    }
  }
  bandgrafptr->fronnbr = grafptr->fronnbr - vfixnum; /* Remove fixed vertices from frontier */
  for (bandvertnnd = bandvertnbr + bandgrafptr->s.baseval - domnnbr; /* Pick selected band vertices from rest of frontier array without anchors */
       bandvertnum < bandvertnnd + bandvfixnbr - vfixnum; fronnum ++) {
    Gnum              vertnum;

    vertnum = frontab[fronnum];
    if ((pfixtax != NULL) && (pfixtax[vertnum] != -1)) /* It is a fixed vertex */
      vfixnum ++;
    else {
      bandvnumtax[bandvertnum] = vertnum;
      bandvertnum ++;
    }
  }
#ifdef SCOTCH_DEBUG_KGRAPH2
  if (vfixnum != bandvfixnbr) {
    errorPrint ("kgraphBand: internal error (1)"); /* All fixed vertices indices must be at the beginning of frontab */
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  memSet (bandvnumtax + bandvertnnd, ~0, domnnbr * sizeof (Gnum)); /* Prevent Valgrind from yelling when centralizing band graphs */

  bandverttax  = bandgrafptr->s.verttax;
  bandvelotax  = bandgrafptr->s.velotax;
  banddegrmax  = 0;
  bandcommload = 0;
  bandcompload = bandgrafptr->comploaddlt;        /* TRICK: use delta array to compute load sums */
  memSet (bandcompload, 0, domnnbr * sizeof (Gnum));
  bandedlosum = 0;
  vfixflag = 0;
  for (bandvertnum = bandedgenum = bandgrafptr->s.baseval; /* Build vertex array of band graph                  */
       bandvertnum < bandvertlvlnum; bandvertnum ++) { /* For all vertices that do not belong to the last level */
    Gnum              vertnum;
    Gnum              edgenum;
    Anum              partval;
    Gnum              degrval;

    vertnum = bandvnumtax[bandvertnum];
    if (vfixflag == 1) {                          /* Last vertex had neighbours fixed vertices      */
      memSet (bandanlotab, 0, domnnbr * sizeof (Gnum)); /* Reset loads to anchors                   */
      vfixflag = 0;                               /* Guess that these are no extra loads to anchors */
    }
    if (vmlotax != NULL)
      bandvmlotax[bandvertnum] = vmlotax[vertnum];
    if (parotax != NULL)
      bandparotax[bandvertnum] = parotax[vertnum];
    partval = parttax[vertnum];
#ifdef SCOTCH_DEBUG_KGRAPH2
    if ((partval < 0) || (partval >= domnnbr)) {
      errorPrint ("kgraphBand: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    bandparttax[bandvertnum] = partval;
    bandverttax[bandvertnum] = bandedgenum;
    veloval = (velotax != NULL) ? velotax[vertnum] : 1;
    bandcompload[partval]   += veloval;           /* Sum vertex load for each part */
    bandvelotax[bandvertnum] = veloval;

    degrval = vendtax[vertnum] - verttax[vertnum];
    if (banddegrmax < degrval)
      banddegrmax = degrval;

    for (edgenum = verttax[vertnum];              /* For all original edges */
         edgenum < vendtax[vertnum]; edgenum ++) {
      Gnum                vertend;

      vertend = edgetax[edgenum];
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (vnumotbdtax[vertend] == -1) {           /* All ends should belong to the band graph too */
        errorPrint ("kgraphBand: internal error (3)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      if (bandedlotax != NULL) {                  /* If graph has edge weights (always true when fixed vertices present) */
        Gnum                edloval;
        Gnum                termnum;

        edloval = (edlotax == NULL) ? 1 : edlotax[edgenum];
        if ((pfixtax != NULL) &&
            ((termnum = pfixtax[vertend]) >= 0)) { /* If end vertex is fixed */
          Gnum                termhashnum;

          for (termhashnum = (termnum * KGRAPHBANDHASHPRIME) & termhashmsk; ; termhashnum = (termhashnum + 1) & termhashmsk) {
            if (termhashtab[termhashnum].termnum == termnum) { /* If hash slot found */
              bandanlotab[termhashtab[termhashnum].domnnum] += edloval;
              vfixflag = 1;                       /* Vertex have some end vertices fixed */
              break;
            }
#ifdef SCOTCH_DEBUG_KGRAPH2
            if (termhashtab[termhashnum].termnum == ~0) { /* If hash slot not found */
              errorPrint ("kgraphBand: missing terminal domain in domain array (1)");
              memFree    (termhashtab);
              kgraphExit (bandgrafptr);
              return     (1);
            }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          }
        }
        else {                                    /* End vertex is not fixed */
          bandedlosum += edloval;
          bandedlotax[bandedgenum] = edloval;
          bandedgetax[bandedgenum ++] = vnumotbdtax[vertend];
        }
      }
      else
        bandedgetax[bandedgenum ++] = vnumotbdtax[vertend];
    }
    
    if (vfixflag == 1) {                          /* If vertex has at least one neighbours that is fixed */
      Gnum                    edloval;

      edloval = 0;
      for (domnnum = 0; domnnum < domnnbr; domnnum ++) { /* Traverse bandanlotab to handle loads to anchors */
        if (bandanlotab[domnnum] != 0)            /* Fixed neighbours are linked to this domain             */
          edloval += bandanlotab[domnnum];        /* Add load induced by edges to fixed vertices            */
        if (edloval != 0) {                       /* We have to add an edge to the anchor                   */
          Gnum              degrval;

          bandedlotax[bandedgenum] = edloval;
          bandedlosum += edloval;
          bandedgetax[bandedgenum ++] = bandvertnnd + domnnum; /* Add edge to anchor of proper part */
          bandeeextab[domnnum] ++;                /* One more extra edge to the anchor              */

          degrval = bandedgenum - bandverttax[bandvertnum];
          if (banddegrmax < degrval)
            banddegrmax = degrval;
          edloval = 0;
        }
      }
    }
  }
  for ( ; bandvertnum < bandvertnnd; bandvertnum ++) { /* For all vertices that belong to the last level except anchors */
    Gnum              vertnum;
    Gnum              edgenum;
    Anum              partval;

    vertnum = bandvnumtax[bandvertnum];
    if (vfixflag == 1) {                          /* Last vertex had neighbours fixed vertices      */ 
      memSet (bandanlotab, 0, domnnbr * sizeof (Gnum)); /* Reset loads to anchors                   */
      vfixflag = 0;                               /* Guess that these are no extra loads to anchors */
    }
    if (vmlotax != NULL)
      bandvmlotax[bandvertnum] = vmlotax[vertnum];
    if (parotax != NULL)
      bandparotax[bandvertnum] = parotax[vertnum];
    partval = parttax[vertnum];
#ifdef SCOTCH_DEBUG_KGRAPH2
    if ((partval < 0) || (partval >= domnnbr)) {
      errorPrint ("kgraphBand: internal error (4)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
    bandparttax[bandvertnum] = partval;
    bandverttax[bandvertnum] = bandedgenum;
    veloval = (velotax != NULL) ? velotax[vertnum] : 1;
    bandcompload[partval]   += veloval;           /* Sum vertex load for each part */
    bandvelotax[bandvertnum] = veloval;

    for (edgenum = verttax[vertnum];              /* For all original edges */
         edgenum < vendtax[vertnum]; edgenum ++) {
      Gnum              vertend;
      Gnum              bandvertend;

      vertend     = edgetax[edgenum];
      bandvertend = vnumotbdtax[vertend];
      if (bandedlotax != NULL) {                  /* If graph has edge weights, copy load */
        Gnum                edloval;
        Gnum                termnum;

        edloval = (edlotax == NULL) ? 1 : edlotax[edgenum];
        if ((pfixtax != NULL) &&
            ((termnum = pfixtax[vertend]) >= 0)) { /* If end vertex is fixed */
          Gnum                termhashnum;

          for (termhashnum = (termnum * KGRAPHBANDHASHPRIME) & termhashmsk; ; termhashnum = (termhashnum + 1) & termhashmsk) {
            if (termhashtab[termhashnum].termnum == termnum) { /* If hash slot found */
              bandanlotab[termhashtab[termhashnum].domnnum] += edloval;
              vfixflag = 1;                       /* Vertex have some end vertices fixed */
              break;
            }
#ifdef SCOTCH_DEBUG_KGRAPH2
            if (termhashtab[termhashnum].termnum == ~0) { /* If hash slot not found */
              errorPrint ("kgraphBand: missing terminal domain in domain array (2)");
              memFree    (termhashtab);
              kgraphExit (bandgrafptr);
              return     (1);
            }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          }
        }
        else if (bandvertend >= 0) {              /* If end vertex is not fixed and in the band graph */
          bandedlosum += edloval;
          bandedlotax[bandedgenum] = edloval;
        }
      }
      if (bandvertend >= 0)                       /* If end vertex is not fixed and in the band graph */
        bandedgetax[bandedgenum ++] = vnumotbdtax[vertend];
    }
    if (vfixflag == 1) {                          /* If vertex has at least one fixed neighbor */
      Gnum                    edloval;

      for (domnnum = 0; domnnum < domnnbr; domnnum ++) { /* Traverse bandanlotab to handle loads to anchors */

        edloval = 0;
        if (domnnum == partval)                   /* The current vertex is mapped to current domain */
          edloval += 1;                           /* Add basic edge load to anchor                  */
        if (bandanlotab[domnnum] != 0)            /* Fixed neighbours are linked to this domain     */
          edloval += bandanlotab[domnnum];        /* Add load induced by edges to fixed vertices    */
        if (edloval != 0) {                       /* We have to add an edge to the anchor           */
          Gnum                degrval;

          bandedlotax[bandedgenum] = edloval;
          bandedlosum += edloval;
          bandedgetax[bandedgenum ++] = bandvertnnd + domnnum; /* Add edge to anchor of proper part */
          if (domnnum != partval)
            bandeeextab[domnnum] ++;              /* One more extra edge to the anchor */

          degrval = bandedgenum - bandverttax[bandvertnum];
          if (banddegrmax < degrval)
            banddegrmax = degrval;
          edloval = 0;
        }
      }
    }
    else {
      Gnum              degrval;

      if (bandedlotax != NULL) {                    /* If graph has edge weights */
        bandedlotax[bandedgenum] = 1;               /* Edge to anchor has load 1 */
        bandedlosum ++;
      }
      bandedgetax[bandedgenum ++] = bandvertnnd + partval; /* Add edge to anchor of proper part */
      degrval = bandedgenum - bandverttax[bandvertnum];
      if (banddegrmax < degrval)
        banddegrmax = degrval;
    }
  }

  memFree (vnumotbdtax + bandgrafptr->s.baseval); /* Free useless space */

  compload = bandgrafptr->comploadavg;            /* Use average array to store actual part loads */
  memSet (compload, 0, domnnbr * sizeof (Gnum));
  for (vertnum = grafptr->s.baseval; vertnum < grafptr->s.vertnnd; vertnum ++)
    compload[parttax[vertnum]] += (velotax != NULL) ? velotax[vertnum] : 1;

  for (domnnum = 0, bandvertancadj = 0; domnnum < domnnbr; domnnum ++) { /* For all anchors */
    Gnum                bandveloval;

    bandveloval = compload[domnnum] - bandcompload[domnnum]; /* Get load of anchor */
    bandparttax[bandvertnnd + domnnum] = domnnum; /* Set parts of anchor vertices  */
    bandvelotax[bandvertnnd + domnnum] = bandveloval;
    if (bandveloval == 0)
      bandvertancadj = 1;
  }

  if (bandvertancadj == 1)                        /* Anchors have to be adjusted                       */
    for (domnnum = 0; domnnum < domnnbr; domnnum ++) /* Increase weight of all anchors to keep balance */
      bandvelotax[bandvertnnd + domnnum] ++;   

  bandverttax[bandvertnum] = bandedgenum;         /* Fill last element without anchors */
  if (pfixtax != NULL)
    memCpy (bandedgetab, bandeeextab, domnnbr * sizeof (Gnum));
  else
    memSet (bandedgetab, 0, domnnbr * sizeof (Gnum));

  for (bandvertnum = bandvertlvlnum; bandvertnum < bandvertnnd; bandvertnum ++)
    bandedgetab[bandparttax[bandvertnum]] ++;  /* Fill array of anchors' degrees */

  for (domnnum = 0; domnnum < domnnbr; domnnum ++) { /* Set bandverttax for anchors vertices and pre-set bandedgetab */
    Gnum                degrval;                     /* to be able to quickly fill bandedgetax in next loop          */
    Gnum                dispval;

    degrval = bandedgetab[domnnum];
    dispval = bandverttax[bandvertnnd + domnnum];
    if (banddegrmax < degrval)                    /* Update maximum degree value */
      banddegrmax = degrval;
    bandverttax[bandvertnnd + domnnum + 1] = dispval + degrval;
    bandedgetab[domnnum] = dispval;               /* Start index for edges to vertices of last layer */
  }

  if (pfixtax != NULL) {                          /* We have fixed vertices */
    memFree (termhashtab);                        /* Free group leader      */
    for (bandvertnum = bandgrafptr->s.baseval, bandedgenum = 0; bandvertnum < bandvertnnd; 
        bandvertnum ++) {                         /* Link anchors to vertices */
      for ( ; bandedgenum < bandverttax[bandvertnum + 1]; bandedgenum ++) {
        Gnum                  bandvertend;

        bandvertend = bandedgetax[bandedgenum];
        if (bandvertend >= bandvertnnd) {         /* If it is an edge to an anchor          */
          Gnum                partval;            /* Add the symmetric edge from the anchor */
          Gnum                edloval;

          partval = bandvertend - bandvertnnd;
          edloval = bandedlotax[bandedgenum];

          bandedlotax[bandedgetab[partval]] = edloval;
          bandedlosum += edloval;
          bandedgetax[bandedgetab[partval] ++] = bandvertnum; 
        }
      }
    }
  }
  else {
    if (bandedlotax != NULL) {                    /* If graph has edge weights */
      Gnum              edgenum;
      Gnum              edgennd;

      for (bandvertnum = bandgrafptr->s.baseval;  /* For all vertices not belonging to last level */
           bandvertnum < bandvertlvlnum; bandvertnum ++) { 
        Gnum              vertnum;
        Gnum              bandedgenum;

        vertnum     = bandvnumtax[bandvertnum];
        bandedgenum = bandverttax[bandvertnum];
        memCpy (&bandedlotax[bandedgenum], &edlotax[verttax[vertnum]], /* Copy edge load array */
                (bandverttax[bandvertnum + 1] - bandedgenum) * sizeof (Gnum));
      }                                           /* Vertices of last level have been processed before */
      for (edgenum = bandverttax[bandvertnnd],    /* Loads of anchor edges are all 1's too       */
         edgennd = bandverttax[bandvertnnd + domnnbr];
         edgenum < edgennd; edgenum ++)
        bandedlotax[edgenum] = 1;
    }
    for (bandvertnum = bandvertlvlnum; bandvertnum < bandvertnnd; /* We do not have fixed vertices */
         bandvertnum ++) {                        /* Link anchors to vertices of last level        */
      Anum              partval;
      Gnum              vertnum;

      vertnum = bandvnumtax[bandvertnum];
      partval = bandparttax[bandvertnum];
      bandedgetax[bandedgetab[partval] ++] = bandvertnum;

      if (bandedlotax != NULL)
        bandedlotax[bandedgetab[partval] - 1] = 1;
      bandedlosum ++;
    }
#ifdef SCOTCH_DEBUG_KGRAPH2
    for (domnnum = 0; domnnum < domnnbr; domnnum ++) {
      if (bandedgetab[domnnum] != bandverttax[bandvertnnd + 1 + domnnum]) {
        errorPrint ("kgraphBand: internal error (6)");
        return     (1);
      }
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
  }  
  bandedgenbr = bandgrafptr->s.verttax[bandvertnnd + domnnbr] - bandgrafptr->s.baseval; /* Set real number of edges */

  bandgrafptr->s.vendtax = bandgrafptr->s.verttax + 1; /* Band graph is compact */
  bandgrafptr->s.velosum = grafptr->s.velosum + domnnbr * bandvertancadj;
  bandgrafptr->s.edgenbr = bandedgenbr;
  if (bandedlotax == NULL) 
    bandedlosum = bandedgenbr;
  bandgrafptr->s.edlosum = bandedlosum;
  bandgrafptr->s.degrmax = banddegrmax;           /* Local maximum degree will be turned into global maximum degree */

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (graphCheck (&bandgrafptr->s) != 0) {
    errorPrint ("kgraphBand: internal error (7)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  memCpy (bandgrafptr->m.domntab, grafptr->m.domntab, domnnbr * sizeof (ArchDom));

  if (pfixtax != NULL)
    kgraphFron (bandgrafptr);

  kgraphCost (bandgrafptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (bandgrafptr) != 0) {
    errorPrint ("kgraphBand: internal error (8)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  *bandvnumptr = bandvnumtax;

  return (0);
}

