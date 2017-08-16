/* Copyright 2004,2007,2008,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_rb.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module performs the Dual Recursive **/
/**                Bipartitioning mapping algorithm.       **/
/**                It is now a branching routine.          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     08 dec 1998     **/
/**                # Version 3.4  : from : 01 jun 2001     **/
/**                                 to     07 nov 2001     **/
/**                # Version 4.0  : from : 12 jan 2004     **/
/**                                 to     06 mar 2005     **/
/**                # Version 5.1  : from : 22 nov 2007     **/
/**                                 to     07 oct 2008     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     28 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_RB

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_dist.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_rb.h"
#include "kgraph_map_rb_map.h"
#include "kgraph_map_rb_part.h"

/********************************************/
/*                                          */
/* This is the entry point for the Dual     */
/* Recursive Bipartitioning mapping method. */
/*                                          */
/********************************************/

/* This routine runs the Dual Recursive
** Bipartitioning algorithm.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRb (
Kgraph * const                          grafptr,
const KgraphMapRbParam * restrict const paraptr)
{
  KgraphMapRbData               datadat;          /* Data passed to each bipartitioning job            */
  Graph                         indgrafdat;       /* Induced graph without fixed vertices              */
  Graph * restrict              indgrafptr;       /* Pointer to top-level graph without fixed vertices */
  KgraphMapRbVflo * restrict    vflotab;          /* Array of fixed vertex load slots                  */
  Anum                          vflonbr;          /* Number of fixed vertex load slots                 */
  Gnum                          vertnum;
  int                           o;

  Gnum * const                frontab = grafptr->frontab;
  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;

  grafptr->kbalval = paraptr->kbalval;            /* Store last k-way imbalance ratio */

  datadat.grafptr = &grafptr->s;
  datadat.mappptr = &grafptr->m;

  datadat.r.mappptr   = (grafptr->r.m.parttax != NULL) ? &grafptr->r.m : NULL;
  datadat.r.vmlotax   = grafptr->r.vmlotax;
  datadat.r.cmloval   = grafptr->r.cmloval;
  datadat.r.crloval   = grafptr->r.crloval;
  datadat.pfixtax     = grafptr->pfixtax;
  datadat.paraptr     = paraptr;
  datadat.comploadrat = grafptr->comploadrat;
  datadat.comploadmin = (1.0 - paraptr->kbalval) * grafptr->comploadrat; /* Ratio can have been tilted when working on subgraph */
  datadat.comploadmax = (1.0 + paraptr->kbalval) * grafptr->comploadrat;

  if (grafptr->pfixtax == NULL) {
    indgrafptr = &grafptr->s;                     /* Work on the original graph */
    vflonbr    = 0;                               /* No fixed vertex load slots */
    vflotab    = NULL;
  }
  else {
    if (kgraphMapRbVfloBuild (grafptr->m.archptr, &grafptr->s, grafptr->vfixnbr, grafptr->pfixtax,
                              &indgrafdat, &vflonbr, &vflotab) != 0) {
      errorPrint ("kgraphMapRb: cannot create induced graph");
      return     (1);
    }
    indgrafptr = &indgrafdat;
  }

  o = ((archPart (grafptr->m.archptr) != 0) ? kgraphMapRbPart : kgraphMapRbMap) (&datadat, indgrafptr, vflonbr, vflotab); /* Compute recursive bipartitioning */

  if (grafptr->pfixtax != NULL) {                 /* If fixed vertices   */
    memFree (vflotab);                            /* Not used any longer */

    graphExit (&indgrafdat);
    if (kgraphMapRbVfloMerge (&grafptr->m, grafptr->vfixnbr, grafptr->pfixtax, vflonbr) != 0) {
      errorPrint ("kgraphMapRb: cannot merge fixed vertex domains");
      return     (1);
    }
  }

  if (memReallocGroup (grafptr->comploadavg,      /* Reallocate cost array according to potential new size                                        */
                       &grafptr->comploadavg, (size_t) (grafptr->m.domnmax * sizeof (Gnum)), /* TRICK: can send both compload arrays in one piece */
                       &grafptr->comploaddlt, (size_t) (grafptr->m.domnmax * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("kgraphMapRb: out of memory (3)");
    return     (1);
  }
  kgraphFron (grafptr);
  kgraphCost (grafptr);                           /* Compute cost of full k-way partition */

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapRb: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (o);
}

/*******************************************/
/*                                         */
/* These routines handle fixed vertex load */
/* arrays and fixed vertices.              */
/*                                         */
/*******************************************/

/* If the given graph has fixed vertices, this
** routine builds simultaneously an induced subgraph
** from which such vertices have been removed, and
** a list of these fixed vertices.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbVfloBuild (
const Arch * restrict const                 archptr, /*+ Target architecture                           +*/
const Graph * restrict const                orggrafptr, /*+ Original graph with fixed vertices         +*/
const Gnum                                  orgvfixnbr, /*+ Number of fixed vertices in graph          +*/
const Anum * restrict const                 orgpfixtax, /*+ Array of fixed vertex terminal domains     +*/
Graph * restrict const                      indgrafptr, /*+ Induced subgraph without fixed vertices    +*/
Anum * restrict const                       vflonbrptr, /*+ Pointer to number of fixed vertex slots    +*/
KgraphMapRbVflo * restrict * restrict const vflotabptr) /*+ Pointer to fixed vertex load array pointer +*/
{
  ArchDom                     domndat;
  Gnum                        orgvertnum;
  GraphPart * restrict        orgparttax;         /* Part array for induced graph            */
  KgraphMapRbVflo * restrict  hashtab;            /* Hash table for merging fixed vertices   */
  Gnum                        hashnbr;            /* Prospective number of cells in table    */
  Gnum                        hashsiz;            /* Size of hash table                      */
  Gnum                        hashnum;
  Gnum                        hashmsk;            /* Mask for access to hash table           */
  Gnum                        velomsk;            /* Zero if all fixed vertex loads are zero */
  Anum                        vflonbr;

  const Gnum * restrict const orgvelotax = orggrafptr->velotax;

  if (archVar (archptr) == 0) {                   /* If fixed size architecture */
    archDomFrst (archptr, &domndat);
    hashnbr = archDomSize (archptr, &domndat);    /* Get maximum size of distinct terminal domains */
    if (orgvfixnbr < hashnbr)
      hashnbr = orgvfixnbr;
  }
  else
    hashnbr = orgvfixnbr;
  for (hashsiz = 0, hashmsk = hashnbr; hashmsk != 0; hashsiz ++, hashmsk >>= 1) ; /* Get upper power of two */
  hashsiz = 1 << (hashsiz + 2);                   /* Fill hash table at 25% maximum                         */
  hashmsk = hashsiz - 1;

  if (memAllocGroup ((void **) (void *)
                     &hashtab,    (size_t) (hashsiz             * sizeof (KgraphMapRbVflo)), /* Use fixed vertex load slots as hash slots */
                     &orgparttax, (size_t) (orggrafptr->vertnbr * sizeof (GraphPart)), NULL) == NULL) {
    errorPrint ("kgraphMapRbVfloBuild: out of memory");
    return     (1);
  }
  orgparttax -= orggrafptr->baseval;

  memSet (hashtab, ~0, hashsiz * sizeof (KgraphMapRbVflo)); /* Set all vertex numbers to ~0 */

  velomsk = 0;                                    /* Assume all fixed vertex loads are zero */
  for (orgvertnum = orggrafptr->baseval; orgvertnum < orggrafptr->vertnnd; orgvertnum ++) {
    Anum                orgpfixval;

    orgpfixval = orgpfixtax[orgvertnum];
    if (orgpfixval >= 0) {                        /* If vertex is a fixed vertex */
      Gnum                hashnum;
      Gnum                veloval;

      veloval  = (orgvelotax == NULL) ? 1 : orgvelotax[orgvertnum]; /* Get fixed vertex load */
      velomsk |= veloval;                         /* See if all fixed vertex loads are zero  */

      for (hashnum = (orgpfixval * KGRAPHMAPRBVFLOHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
        if (hashtab[hashnum].termnum == orgpfixval) { /* If hash slot found   */
          hashtab[hashnum].veloval += veloval;    /* Add contribution to slot */
          break;
        }
        if (hashtab[hashnum].termnum == ~0) {     /* If hash slot empty */
          hashtab[hashnum].termnum = orgpfixval;  /* Create slot        */
          hashtab[hashnum].veloval = veloval;
          break;
        }
      }
      orgparttax[orgvertnum] = 1;                 /* Fixed vertex will not be kept in induced subgraph */
    }
    else
      orgparttax[orgvertnum] = 0;                 /* Keep non-fixed vertex in induced subgraph */
  }

  if (graphInducePart (orggrafptr, orgparttax, orggrafptr->vertnbr - orgvfixnbr, 0, indgrafptr) != 0) { /* Keep non-fixed vertices in induced graph */
    errorPrint ("kgraphMapRbVfloBuild: cannot build induced subgraph");
    memFree    (hashtab);
    return     (1);
  }

  if (velomsk == 0) {                             /* If all fixed vertex loads are zero */
    memFree (hashtab);                            /* No need to allocate a table        */
    *vflonbrptr = 0;
    *vflotabptr = NULL;
    return (0);
  }

  for (hashnum = vflonbr = 0; hashnum < hashsiz; hashnum ++) { /* Recycle hash table into fixed vertex load table */
    if (hashtab[hashnum].termnum != ~0) {
      hashtab[vflonbr] = hashtab[hashnum];
      vflonbr ++;
    }
  }

  *vflonbrptr = vflonbr;
  *vflotabptr = memRealloc (hashtab, vflonbr * sizeof (KgraphMapRbVflo));

  return (0);
}

/* This routine splits the fixed vertex load
** array in two parts, one for each of the two
** provided subdomains.
** It returns:
** - void  : in all cases.
*/

void
kgraphMapRbVfloSplit (
const Arch * restrict const       archptr,        /*+ Target architecture                 +*/
const ArchDom * restrict const    domnsubtab,     /*+ Array of the two subdomains         +*/
const Anum                        vflonbr,        /*+ Number of fixed vertex load slots   +*/
KgraphMapRbVflo * restrict const  vflotab,        /*+ Fixed vertex load array             +*/
Anum * restrict const             vflonbrtab,     /*+ Number of slots in each subdomain   +*/
Gnum * restrict const             vflowgttab)     /*+ Fixed vertex load in each subdomain +*/
{
  KgraphMapRbVflo     vflodat;                    /* Temporary swap area                      */
  ArchDom             domndat;                    /* Terminal domain attached to fixed vertex */
  Gnum                compload0;                  /* Load of slots in first subdomain         */
  Gnum                compload1;                  /* Load of slots in second subdomain        */
  Gnum                vflomax;
  Gnum                vflonum;
  Gnum                vflonnd;

  compload0 =
  compload1 = 0;

  vflomax = vflonbr;
  if (archVar (archptr) == 0) {
    for (vflonum = 0, vflonnd = vflonbr - 1; vflonum < vflonnd; ) {
      while (1) {
#ifdef SCOTCH_DEBUG_KGRAPH2
        int                 o;

        o =
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        archDomTerm (archptr, &domndat, vflotab[vflonum].termnum);
#ifdef SCOTCH_DEBUG_KGRAPH2
        if (o != 0) {
          errorPrint ("kgraphMapRbVfloSplit: internal error (1)");
          return;
        }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
        if (archDomIncl (archptr, &domnsubtab[0], &domndat) == 1) { /* If terminal vertex subdomain included in first subdomain */
          compload0 += vflotab[vflonum].veloval;  /* Fixed vertex belongs to first subdomain                                    */
          if (++ vflonum > vflonnd)               /* If passed beyond the limit of the second subdomain                         */
            goto quit;
        }
        else {
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (archDomIncl (archptr, &domnsubtab[1], &domndat) != 1) { /* If terminal vertex subdomain not included in second subdomain */
            errorPrint ("kgraphMapRbVfloSplit: internal error (2)");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          break;                                  /* We have found a candidate for swapping */
        }
      }
      while (1) {
        archDomTerm (archptr, &domndat, vflotab[vflonnd].termnum);
        if (archDomIncl (archptr, &domnsubtab[1], &domndat) == 1) { /* If terminal vertex subdomain included in second subdomain */
          compload1 += vflotab[vflonnd].veloval;  /* Fixed vertex belongs to second subdomain                                    */
          if (-- vflonnd <= vflonum) {            /* If matched the location of a slot that also belongs to second subdomain     */
            compload1 += vflotab[vflonnd].veloval; /* Add load of reached slot to second subdomain                               */
            goto quit;
          }
        }
        else {
#ifdef SCOTCH_DEBUG_KGRAPH2
          if (archDomIncl (archptr, &domnsubtab[0], &domndat) != 1) { /* If terminal vertex subdomain not included in first subdomain */
            errorPrint ("kgraphMapRbVfloSplit: internal error (3)");
            return;
          }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
          break;                                  /* We have found a candidate for swapping */
        }
      }

      vflodat          = vflotab[vflonum];        /* Swap slots */
      vflotab[vflonum] = vflotab[vflonnd];
      vflotab[vflonnd] = vflodat;
      compload0 += vflotab[vflonum ++].veloval;
      compload1 += vflotab[vflonnd --].veloval;
    }
  }
  else {                                          /* If variable-sized architecture, pseudo-terminals may not always be included */
    for (vflonum = 0, vflonnd = vflonbr - 1; vflonum <= vflonnd; ) {
#ifdef SCOTCH_DEBUG_KGRAPH2
      int                 o;

      o =
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      archDomTerm (archptr, &domndat, vflotab[vflonum].termnum);
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (o != 0) {
        errorPrint ("kgraphMapRbVfloSplit: internal error (4)");
        return;
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      if (archDomIncl (archptr, &domnsubtab[0], &domndat) == 1) { /* If vertex subdomain included in first subdomain */
        compload0 += vflotab[vflonum].veloval;    /* Fixed vertex belongs to first subdomain                         */
        vflonum ++;
        continue;                                 /* Keep it in place */
      }
      if (archDomIncl (archptr, &domnsubtab[1], &domndat) == 1) { /* If vertex subdomain included in second subdomain */
        compload1 += vflotab[vflonum].veloval;    /* Fixed vertex belongs to second subdomain                         */

        vflodat          = vflotab[vflonum];      /* Swap slots */
        vflotab[vflonum] = vflotab[vflonnd];
        vflotab[vflonnd] = vflodat;
      }
      else {                                      /* Fixed vertex is more generic than the two subdomains */
        vflomax --;                               /* One less slot to consider in the future              */
        vflodat          = vflotab[vflonum];      /* Swap slots */
        vflotab[vflonum] = vflotab[vflonnd];
        vflotab[vflonnd] = vflotab[vflomax];
        vflotab[vflomax] = vflodat;
      }
      vflonnd --;                                 /* One less slot to consider */
    }
  }

quit:
  vflonbrtab[0] = vflonum;
  vflonbrtab[1] = vflomax - vflonum;
  vflowgttab[0] = compload0;
  vflowgttab[1] = compload1;
}

/* This routine prolongs the given mapping
** with fixed vertices, by merging back with
** regular domains.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbVfloMerge (
Mapping * restrict const    mappptr,              /*+ Mapping to prolong                     +*/
const Gnum                  vfixnbr,              /*+ Number of fixed vertices in graph      +*/
const Anum * restrict const pfixtax,              /*+ Array of fixed vertex terminal domains +*/
const Anum                  vflonbr)              /*+ Number of fixed vertex load slots      +*/
{
  ArchDom                         domndat;
  Anum                            domnmax;
  Anum                            domnnum;
  Gnum                            vertnum;
  Gnum                            vertnnd;
  KgraphMapRbVfloHash * restrict  hashtab;        /* Hash table for merging fixed vertices   */
  Gnum                            hashnbr;        /* Prospective number of cells in table    */
  Gnum                            hashsiz;        /* Size of hash table                      */
  Gnum                            hashnum;
  Gnum                            hashmsk;        /* Mask for access to hash table           */

  const Arch * restrict const archptr = mappptr->archptr;
  Anum * restrict const       parttax = mappptr->parttax;

  hashnbr = mappptr->domnnbr + vflonbr;
  for (hashsiz = 0, hashmsk = hashnbr; hashmsk != 0; hashsiz ++, hashmsk >>= 1) ; /* Get upper power of two */
  hashsiz = 1 << (hashsiz + 2);                   /* Fill hash table at 25% maximum                         */
  hashmsk = hashsiz - 1;

  if ((hashtab = memAlloc (hashsiz * sizeof (KgraphMapRbVfloHash))) == NULL) { /* Use fixed vertex load slots as hash slots */
    errorPrint ("kgraphMapRbVfloMerge: out of memory (1)");
    return     (1);
  }
  memSet (hashtab, ~0, hashsiz * sizeof (KgraphMapRbVfloHash)); /* Set all vertex numbers to ~0 */

  for (domnnum = 0; domnnum < mappptr->domnnbr; domnnum ++) { /* Load all existing domains into hash table */
    Anum                termnum;

    termnum = archDomNum (archptr, &mappptr->domntab[domnnum]);
    for (hashnum = (termnum * KGRAPHMAPRBVFLOHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
      if (hashtab[hashnum].termnum == termnum)    /* If domain found */
        break;
      if (hashtab[hashnum].termnum == ~0) {       /* If empty slot found */
        hashtab[hashnum].termnum = termnum;       /* Fill it             */
        hashtab[hashnum].domnnum = domnnum;
        break;
      }
    }
  }

  domnmax = mappptr->domnmax;
  for (vertnum = mappptr->grafptr->baseval, vertnnd = mappptr->grafptr->vertnnd; vertnum < vertnnd; vertnum ++) {
    Anum                pfixval;

    pfixval = pfixtax[vertnum];
    if (pfixval < 0) {                            /* If vertex is not a fixed vertex */
#ifdef SCOTCH_DEBUG_KGRAPH2
      if (mappptr->parttax[vertnum] < 0) {        /* If vertex has not been mapped */
        errorPrint ("kgraphMapRbVfloMerge: internal error (1)");
        return;
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      continue;                                   /* Skip to next vertex */
    }

#ifdef SCOTCH_DEBUG_KGRAPH2
    if (mappptr->parttax[vertnum] >= 0) {         /* If fixed vertex has been mapped */
      errorPrint ("kgraphMapRbVfloMerge: internal error (2)");
      return;
    }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

    for (hashnum = (pfixval * KGRAPHMAPRBVFLOHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
      if (hashtab[hashnum].termnum == pfixval)    /* If hash slot found */
        break;
      if (hashtab[hashnum].termnum == ~0) {       /* If hash slot empty */
        if (domnnum >= mappptr->domnmax) {
          if (mapResize (mappptr, mappptr->domnmax + (mappptr->domnmax >> 2) + 8) != 0) { /* Increase size by 25% */
            errorPrint ("kgraphMapRbVfloMerge: out of memory (2)");
            return     (1);
          }
        }
        archDomTerm (archptr, &mappptr->domntab[domnnum], pfixval); /* Add new domain to domain array */

        hashtab[hashnum].termnum = pfixval;       /* Create slot */
        hashtab[hashnum].domnnum = domnnum;
        domnnum ++;                               /* One more domain created */
        break;
      }
    }
    parttax[vertnum] = hashtab[hashnum].domnnum;  /* Assign fixed vertex to existing domain */
  }
  mappptr->domnnbr = domnnum;

  memFree (hashtab);

  return (0);
}

/*****************************************/
/*                                       */
/* This routine computes external gains. */
/*                                       */
/*****************************************/

/* This routines computes the three kinds of
** external loads and gains that can be associated
** with a bipartitioning job: regular external
** gains arising from former bipartitions, fixed
** vertex gains and remapping gains.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kgraphMapRbBgraph (
const KgraphMapRbData * restrict const  dataptr,  /*+ Global data                     +*/
Bgraph * restrict const                 actgrafptr, /*+ Graph to build                +*/
const Graph * restrict const            srcgrafptr, /*+ Source graph                  +*/
const Mapping * restrict const          srcmappptr, /*+ Current mapping               +*/
const ArchDom * restrict const          domnsubtab, /*+ Array of the two subdomains   +*/
const Gnum * restrict const             vflowgttab) /*+ Array of vertex weight biases +*/
{
  Gnum                  actvertnum;               /* Number of current active vertex   */
  Gnum                  commloadextn0;            /* External communication load       */
  Gnum                  commgainextn0;            /* External communication gain       */
  Gnum * restrict       veextax;                  /* External gain array               */
  Gnum                  veexmsk;                  /* Flag set if external array useful */
  int                   flagval;
  int                   o;

  const Arch * restrict const     archptr    = dataptr->mappptr->archptr;
  const Gnum * restrict const     orgverttax = dataptr->grafptr->verttax;
  const Gnum * restrict const     orgvendtax = dataptr->grafptr->vendtax;
  const Gnum * restrict const     orgvelotax = dataptr->grafptr->velotax;
  const Gnum * restrict const     orgedgetax = dataptr->grafptr->edgetax;
  const Gnum * restrict const     orgedlotax = dataptr->grafptr->edlotax;
  const Mapping * restrict const  oldmappptr = dataptr->r.mappptr;
  const Gnum * restrict const     orgvmlotax = dataptr->r.vmlotax;
  const Anum * restrict const     orgpfixtax = dataptr->pfixtax;
  const Gnum * restrict const     actverttax = srcgrafptr->verttax; /* Get pointers from source graph before bgraphInit() */
  const Gnum * restrict const     actvendtax = srcgrafptr->vendtax;
  const Gnum * restrict const     actedgetax = srcgrafptr->edgetax;
  const Gnum * restrict const     actvnumtax = srcgrafptr->vnumtax;

  if (bgraphInit (actgrafptr, srcgrafptr, srcmappptr->archptr, domnsubtab, vflowgttab) != 0) {
    errorPrint ("kgraphMapRbBgraph: cannot create bipartition graph");
    return     (1);
  }

  flagval = KGRAPHMAPRBVEEXNONE;                  /* Assume no processing */
  if ((! archPart (archptr)) && (actvnumtax != NULL))
    flagval |= KGRAPHMAPRBVEEXMAPP;
  if (orgpfixtax != NULL)                         /* Fixed vertices always imply (actvnumtax != NULL) */
    flagval |= KGRAPHMAPRBVEEXVFIX;
  if (dataptr->r.mappptr != NULL)
    flagval |= KGRAPHMAPRBVEEXREMA;

  if (flagval == KGRAPHMAPRBVEEXNONE)             /* If nothing to do */
    return (0);

  if ((veextax = (Gnum *) memAlloc (actgrafptr->s.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("kgraphMapRbBgraph: out of memory");
    return     (1);
  }
  veextax -= actgrafptr->s.baseval;

  o = 1;                                          /* Assume failure                */
  veexmsk = 0;                                    /* No useful array entry yet     */
  commloadextn0 =                                 /* No external communication yet */
  commgainextn0 = 0;
  for (actvertnum = actgrafptr->s.baseval;        /* Compute external loads */
       actvertnum < actgrafptr->s.vertnnd; actvertnum ++) {
    Gnum                commloadextn;             /* External communication load for current vertex */
    Gnum                commgainextn;             /* External communication gain for current vertex */
    Gnum                orgvertnum;               /* Number of current original vertex              */

    commloadextn =                                /* Assume no external loads */
    commgainextn = 0;

    if (actvnumtax == NULL)                       /* If active graph is not a subgraph      */
      orgvertnum = actvertnum;                    /* Its number is that of original graph   */
    else {                                        /* Else we have external edges to process */
      orgvertnum = actvnumtax[actvertnum];        /* Get vertex number in original graph    */

      if ((flagval & KGRAPHMAPRBVEEXEDGE) != 0) { /* If edge-based gains have to be computed */
        Gnum                orgedgenum;
        Gnum                orgedgennd;
        Gnum                actedgenum;
        Gnum                actedgennd;

        orgedgenum = orgverttax[orgvertnum];
        orgedgennd = orgvendtax[orgvertnum];
        actedgenum = actverttax[actvertnum];
        actedgennd = actvendtax[actvertnum];
        if ((orgedgennd - orgedgenum) != (actedgennd - actedgenum)) { /* If vertex has external edges */
          Gnum                orgedloval;
          Gnum                actvertend;         /* Next internal end vertex index to find in original edge array */

          orgedloval = 1;                         /* Assume no edge loads */
          actvertend = (actedgenum >= actedgennd) ? ~0 : actvnumtax[actedgetax[actedgenum]];
          for ( ; orgedgenum < orgedgennd; orgedgenum ++) {
            Gnum                orgvertend;

            orgvertend = orgedgetax[orgedgenum];
            if (orgvertend == actvertend) {       /* If internal edge found    */
              actedgenum ++;                      /* Skip internal active edge */
              actvertend = (actedgenum >= actedgennd) ? ~0 : actvnumtax[actedgetax[actedgenum]]; /* Set next internal end vertex index to fetch */
              continue;                           /* Skip internal original edge */
            }

            if (orgedlotax != NULL)
              orgedloval = orgedlotax[orgedgenum];

            if (orgpfixtax != NULL) {
              ArchDom             domndat;
              Anum                pfixval;

              pfixval = orgpfixtax[orgvertend];
              if (pfixval >= 0) {                 /* If end vertex is fixed */
#ifdef SCOTCH_DEBUG_KGRAPH2
                if (dataptr->mappptr->parttax[orgvertend] != ~0) { /* If original vertex has a current part */
                  errorPrint ("kgraphMapRbBgraph: internal error");
                  goto fail;
                }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

                if (archDomTerm (archptr, &domndat, pfixval) != 0) { /* Get its domain */
                  errorPrint ("kgraphMapRbBgraph: invalid fixed part array");
                  goto fail;
                }

                commloadextn += (archDomIncl (archptr, &domnsubtab[0], &domndat) != 0)
                                ? 0 : (orgedloval * archDomDist (archptr, &domnsubtab[0], &domndat));
                commgainextn += (archDomIncl (archptr, &domnsubtab[1], &domndat) != 0)
                                ? 0 : (orgedloval * archDomDist (archptr, &domnsubtab[1], &domndat));

                continue;                         /* Process next edge */
              }
            }

            if ((flagval & KGRAPHMAPRBVEEXMAPP) != 0) { /* If mapping */
              ArchDom *           domnptr;

              domnptr = mapDomain (srcmappptr, orgvertend);
              commloadextn += orgedloval * archDomDist (archptr, &domnsubtab[0], domnptr);
              commgainextn += orgedloval * archDomDist (archptr, &domnsubtab[1], domnptr);
            }
          }

          commloadextn *= dataptr->r.crloval;     /* Scale regular, non-remapping gains */
          commgainextn *= dataptr->r.crloval;
        }
      }
    }

    if (oldmappptr != NULL) {                     /* If remapping gains have to be computed */
      ArchDom *           domnptr;
      Gnum                edloval;

      edloval = dataptr->r.cmloval;
      if (orgvmlotax != NULL)
        edloval *= orgvmlotax[orgvertnum];

      domnptr = mapDomain (oldmappptr, orgvertnum);
      commloadextn += (archDomIncl (archptr, &domnsubtab[0], domnptr) != 0)
                      ? 0 : (edloval * archDomDist (archptr, &domnsubtab[0], domnptr));
      commgainextn += (archDomIncl (archptr, &domnsubtab[1], domnptr) != 0)
                      ? 0 : (edloval * archDomDist (archptr, &domnsubtab[1], domnptr));
    }

    commgainextn  -= commloadextn;                /* Compute vertex gain        */
    commloadextn0 += commloadextn;                /* Account for external edges */
    commgainextn0 += commgainextn;

    veextax[actvertnum] = commgainextn;           /* Record external gain value */
    veexmsk            |= commgainextn;           /* Accumulate non-zero values */
  }
  o = 0;                                          /* Computations succeeded */

fail:
  if ((o != 0) || (veexmsk == 0)) {               /* If external gain array is useless */
    memFree (veextax + actgrafptr->s.baseval);    /* Forget about it                   */
    return  (o);                                  /* Return error code                 */
  }

  actgrafptr->s.flagval |= BGRAPHFREEVEEX;        /* Keep external gain array */
  actgrafptr->veextax    = veextax;

  actgrafptr->commload      = commloadextn0;      /* Account for external gains in future computations */
  actgrafptr->commgainextn  = commgainextn0;
  actgrafptr->commloadextn0 = commloadextn0;
  actgrafptr->commgainextn0 = commgainextn0;

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (bgraphCheck (actgrafptr) != 0) {
    errorPrint ("kgraphMapRbBgraph: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}
