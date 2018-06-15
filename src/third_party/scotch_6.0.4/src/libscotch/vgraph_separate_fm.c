/* Copyright 2004,2007,2008,2010,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vgraph_separate_fm.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates an active         **/
/**                graph using a vertex-oriented version   **/
/**                of our improved Fiduccia-Mattheyses     **/
/**                heuristics, similar in principle to     **/
/**                the algorithm of Ashcraft and Liu 1994. **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 02 nov 1997     **/
/**                                 to     17 jul 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     31 dec 1998     **/
/**                # Version 4.0  : from : 07 jan 2002     **/
/**                                 to     18 aug 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     22 may 2008     **/
/**                # Version 5.1  : from : 10 nov 2008     **/
/**                                 to     01 jun 2010     **/
/**                # Version 6.0  : from : 31 mar 2014     **/
/**                                 to     01 apr 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VGRAPH_SEPARATE_FM

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_gg.h"
#include "vgraph_separate_fm.h"

/*
**  The static definitions.
*/

static VgraphSeparateFmVertex vexxdat;            /* Dummy structure for computing offsets */

/*********************************/
/*                               */
/* Gain table handling routines. */
/*                               */
/*********************************/

/* This routine returns the vertex of best gain
** whose swap will keep the balance correct.
** It returns:
** - !NULL  : pointer to the vertex gainlink.
** - NULL   : if no more vertices available.
*/

static
GainLink *
vgraphSeparateFmTablGet (
GainTabl * const            tablptr,              /* Gain table        */
const Gnum                  deltcur,              /* Current imbalance */
const Gnum                  deltmax,              /* Maximum imbalance */
const int                   partval)              /* Current preferred */
{
  const VgraphSeparateFmVertex *  vexxptr;        /* Pointer to vertex of current link */
  const GainLink *                linkptr;        /* Pointer to current gain link      */
  const GainLink *                linkbest;       /* Pointer to best link found        */
  const GainEntr *                tablbest;       /* Gain table entry of best link     */
  Gnum                            gaincur;        /* Separator gain of current link    */
  Gnum                            gainbest;       /* Separator gain of best link       */

  linkbest = NULL;                                /* Assume no candidate vertex found yet */
  tablbest = tablptr->tend;
  gainbest = GAINMAX;

  for (linkptr = gainTablFrst (tablptr);          /* Select candidate vertices */
       (linkptr != NULL) && (linkptr->tabl <= tablbest);
       linkptr = gainTablNext (tablptr, linkptr)) {
    int                       vertpart;           /* Part of current vertex */

    vertpart = 0;                                 /* Assume we point to gainlink0      */
    vexxptr  = (VgraphSeparateFmVertex *) linkptr; /* TRICK: gainlink0 is at beginning */
    if (vexxptr->veloval >= 0) {                  /* If in fact we point to gainlink1  */
      vertpart = 1;                               /* Then point to vertex structure    */
      vexxptr  = (VgraphSeparateFmVertex *) ((byte *) vexxptr - ((byte *) &vexxdat.gainlink1 - (byte *) &vexxdat));
    }
    gaincur = vexxptr->compgain[vertpart];        /* Get separator gain and vertex balance */

    if (gaincur == vexxptr->veloval)              /* If vertex is isolated separator vertex */
      return ((GainLink *) linkptr);              /* Select it immediatly                   */

    if (abs (deltcur + (1 - 2 * vertpart) * (gaincur - 2 * vexxptr->veloval)) <= deltmax)  { /* If vertex enforces balance */
      if ((gaincur < gainbest) ||                 /* And if it gives better gain */
          ((gaincur == gainbest) &&               /* Or is in preferred part     */
           (partval == vertpart))) {
        linkbest = linkptr;                       /* Select it */
        tablbest = linkptr->tabl;
        gainbest = gaincur;
      }
    }
  }

  return ((GainLink *) linkbest);                 /* Return best link found */
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

int
vgraphSeparateFm (
Vgraph * restrict const             grafptr,      /*+ Active graph      +*/
const VgraphSeparateFmParam * const paraptr)      /*+ Method parameters +*/
{
  GainTabl * restrict             tablptr;        /* Pointer to gain tables                  */
  INT                             passnbr;        /* Maximum number of passes to go          */
  Gnum                            movenbr;        /* Number of uneffective moves done        */
  int                             moveflag;       /* Flag set if useful moves made           */
  Gnum                            savenbr;        /* Position of current move                */
  VgraphSeparateFmSave * restrict savetab;        /* Pointer to move array                   */
  Gnum                            hashmax;        /* Maximum number of elements in table     */
  Gnum                            hashsiz;        /* Size of hash and save tables            */
  Gnum                            hashmsk;        /* Mask for access to hash table           */
  Gnum                            hashnbr;        /* Number of elements in hash table        */
  VgraphSeparateFmVertex *        hashtab;        /* Hash vertex table                       */
  GainLink                        lockdat;        /* Double linked list of locked vertices   */
  VgraphSeparateFmVertex *        vexxptr;        /* Pointer to current vertex               */
  VgraphSeparateFmVertex *        sepaptr;        /* Pointer to current vertex in table      */
  Gnum                            fronnum;        /* Current index of frontier vertex        */
  Gnum                            comploaddlt;    /* Current load imbalance                  */
  Gnum                            comploaddltmat; /* Theoretical maximum unbalance           */
  Gnum                            comploaddltmax; /* Largest unbalance allowed               */
  Gnum                            comploaddltbst; /* Unbalance of best solution to date      */
  Gnum                            compload2;      /* Current load of separator               */
  Gnum                            compload2bst;   /* Separator load of best solution to date */
  Gnum                            mswpnum;        /* Number of current move sweep            */
  Gnum                            compsize1add;   /* Number of vertices to add to counters   */
  Gnum                            compsize1sub;

  const Gnum * restrict const verttax = grafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;
  GraphPart * restrict const  parttax = grafptr->parttax;

  comploaddltmat = (paraptr->deltrat > 0.0L)
                   ? MAX ((Gnum) ((grafptr->compload[0] + grafptr->compload[1]) * paraptr->deltrat),
                          ((2 * grafptr->s.velosum) / grafptr->s.vertnbr))
                   : 0;

  if (grafptr->fronnbr == 0) {                    /* If no frontier defined     */
    if (abs (grafptr->comploaddlt) <= comploaddltmat) /* If balance is achieved */
      return (0);                                 /* This algorithm is useless  */
    else {                                        /* Imbalance must be fought   */
      VgraphSeparateGgParam paradat;

      paradat.passnbr = 4;                        /* Use a standard algorithm */
      vgraphSeparateGg (grafptr, &paradat);
      if (grafptr->fronnbr == 0)                  /* If new partition has no frontier */
        return (0);                               /* This algorithm is still useless  */
    }
  }

  hashnbr = 16 * (grafptr->fronnbr + paraptr->movenbr + grafptr->s.degrmax) + 1;
#ifdef SCOTCH_DEBUG_VGRAPH2
  hashnbr /= 8;                                   /* Ensure resizing routine will be called */
#endif /* SCOTCH_DEBUG_VGRAPH2 */
  if (hashnbr > grafptr->s.vertnbr)
    hashnbr = grafptr->s.vertnbr;
  for (hashsiz = 512; hashsiz < hashnbr; hashsiz <<= 1) ; /* Get upper power of two */
  hashmsk = hashsiz - 1;
  hashmax = hashsiz >> 2;                         /* Use hash table at 1/4 of its capacity */

  if (((tablptr = gainTablInit (GAINMAX, VGRAPHSEPAFMGAINBITS)) == NULL) || /* Use logarithmic array only */
      (memAllocGroup ((void **) (void *)
                      &hashtab, (size_t) (hashsiz * sizeof (VgraphSeparateFmVertex)),
                      &savetab, (size_t) (hashsiz * sizeof (VgraphSeparateFmSave)), NULL) == NULL)) {
    errorPrint ("vgraphSeparateFm: out of memory (1)");
    if (tablptr != NULL)
      gainTablExit (tablptr);
    return (1);
  }
  memSet (hashtab, ~0, hashsiz * sizeof (VgraphSeparateFmVertex)); /* Set all vertex numbers to ~0 */

  for (fronnum = 0, hashnbr = grafptr->fronnbr;   /* Set initial gains */
       fronnum < hashnbr; fronnum ++) { 
    Gnum                vertnum;
    Gnum                hashnum;

    vertnum = grafptr->frontab[fronnum];
#ifdef SCOTCH_DEBUG_VGRAPH2
    if (parttax[vertnum] != 2) {
      errorPrint   ("vgraphSeparateFm: vertex not in separator");
      return       (1);
    }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

    for (hashnum = (vertnum * VGRAPHSEPAFMHASHPRIME) & hashmsk; hashtab[hashnum].vertnum != ~0; hashnum = (hashnum + 1) & hashmsk) ;

    if (velotax != NULL) {                        /* If vertex loads present */
      Gnum                edgenum;
      Gnum                veloval;
      Gnum                compgain0;
      Gnum                compgain01;

      for (edgenum = verttax[vertnum], compgain0 = compgain01 = 0;
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;
        Gnum                partend;
        Gnum                veloend;

        vertend = edgetax[edgenum];
        partend = (Gnum) parttax[vertend];
        veloend = velotax[vertend];

        compgain0  += (partend & 1)       * veloend;
        compgain01 += (2 - (partend & 2)) * veloend;
      }
      veloval = velotax[vertnum];
      hashtab[hashnum].veloval     = - veloval;   /* TRICK: -veloval: stored value is opposite of load */
      hashtab[hashnum].compgain[0] = compgain0 - veloval;
      hashtab[hashnum].compgain[1] = (compgain01 >> 1) - compgain0 - veloval;
    }
    else {                                        /* No vertex loads */
      Gnum                edgenum;
      Gnum                compgain0;
      Gnum                compgain2;

      for (edgenum = verttax[vertnum], compgain0 = compgain2 = 0;
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;
        Gnum                partend;

        vertend = edgetax[edgenum];
        partend = (Gnum) parttax[vertend];

        compgain0 += (partend & 1);
        compgain2 += (partend & 2);
      }
      hashtab[hashnum].veloval     = -1;          /* TRICK: -veloval */
      hashtab[hashnum].compgain[0] = compgain0 - 1;
      hashtab[hashnum].compgain[1] = vendtax[vertnum] - verttax[vertnum] - (compgain2 >> 1) - compgain0 - 1;
    }
    hashtab[hashnum].partval = 2;
    hashtab[hashnum].vertnum = vertnum;

    gainTablAdd (tablptr, &hashtab[hashnum].gainlink0, hashtab[hashnum].compgain[0]); /* Link both directions of separator vertex */
    gainTablAdd (tablptr, &hashtab[hashnum].gainlink1, hashtab[hashnum].compgain[1]);
  }

  comploaddltmax = MAX (comploaddltmat, abs (grafptr->comploaddlt)); /* Set current maximum distance */
  comploaddltbst = grafptr->comploaddlt;
  compload2bst   = grafptr->compload[2];

#ifdef SCOTCH_DEBUG_VGRAPH3
  if (vgraphSeparateFmCheck (grafptr, hashtab, hashmsk, compload2bst, comploaddltbst) != 0) {
    errorPrint ("vgraphSeparateFm: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH3 */

  passnbr = paraptr->passnbr;                     /* Set remaining number of passes    */
  savenbr = 0;                                    /* For empty backtrack of first pass */
  mswpnum = -1;                                   /* Will be incremented afterwards    */
  lockdat.next =                                  /* List of locked vertices is empty  */
  lockdat.prev = &lockdat;
  do {                                            /* As long as there is improvement */
    Gnum                comploadabsdltbst;

    while (savenbr -- > 0) {                      /* Delete exceeding moves */
      Gnum                hashnum;
      int                 partval;

      hashnum = savetab[savenbr].hashnum;
      partval = savetab[savenbr].partval;
      hashtab[hashnum].partval     = partval;     /* Restore vertex data */
      hashtab[hashnum].compgain[0] = savetab[savenbr].compgain[0];
      hashtab[hashnum].compgain[1] = savetab[savenbr].compgain[1];

      if (hashtab[hashnum].gainlink0.next >= VGRAPHSEPAFMSTATELINK) { /* If vertex is linked */
        gainTablDel (tablptr, &hashtab[hashnum].gainlink0); /* Unlink it                     */
        gainTablDel (tablptr, &hashtab[hashnum].gainlink1);
        hashtab[hashnum].gainlink0.next = VGRAPHSEPAFMSTATEFREE; /* Set it as free */
      }
      if ((hashtab[hashnum].gainlink0.next == VGRAPHSEPAFMSTATEFREE) && (partval == 2)) { /* If vertex not locked and in separator */
        gainTablAdd (tablptr, &hashtab[hashnum].gainlink0, hashtab[hashnum].compgain[0]); /* Re-link it                            */
        gainTablAdd (tablptr, &hashtab[hashnum].gainlink1, hashtab[hashnum].compgain[1]);
      }
    }
    compload2   = compload2bst;                   /* Restore best separator parameters */
    comploaddlt = comploaddltbst;
    comploadabsdltbst = abs (comploaddltbst);
    if (comploadabsdltbst > comploaddltmax)       /* If the former state had a higher maximum imbalance ratio */
      comploaddltmax = comploadabsdltbst;         /* Restore this maximum imbalance ratio                     */

    mswpnum ++;                                   /* Forget all recorded moves */

#ifdef SCOTCH_DEBUG_VGRAPH3
    if (vgraphSeparateFmCheck (grafptr, hashtab, hashmsk, compload2, comploaddlt) != 0) {
      errorPrint ("vgraphSeparateFm: internal error (2)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VGRAPH3 */

    while (lockdat.next != &lockdat) {            /* For all vertices in locked list */
      VgraphSeparateFmVertex *  vexxptr;

      vexxptr      = (VgraphSeparateFmVertex *) ((byte *) lockdat.next - ((byte *) &vexxdat.gainlink1 - (byte *) &vexxdat));
      lockdat.next = (GainLink *) vexxptr->gainlink1.next; /* Unlink vertex from list */

      if (vexxptr->partval == 2) {                /* If vertex belongs to separator  */
#ifdef SCOTCH_DEBUG_VGRAPH2
        if (vexxptr->gainlink0.next != VGRAPHSEPAFMSTATEUSED) {
          errorPrint ("vgraphSeparateFm: linked non-used vertex");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
        gainTablAdd (tablptr, &vexxptr->gainlink0, vexxptr->compgain[0]); /* Link it */
        gainTablAdd (tablptr, &vexxptr->gainlink1, vexxptr->compgain[1]);
      }
      else                                        /* Vertex does not belong to separator */
        vexxptr->gainlink0.next = VGRAPHSEPAFMSTATEFREE; /* Set it as free for this run  */
    }
    lockdat.prev = &lockdat;                      /* Restore backward chaining */

    moveflag = 0;                                 /* No moves to date                          */
    movenbr  =                                    /* No uneffective moves yet                  */
    savenbr  = 0;                                 /* No recorded moves yet                     */
    while ((movenbr < paraptr->movenbr) &&        /* As long as we can find effective vertices */
           ((vexxptr = (VgraphSeparateFmVertex *) vgraphSeparateFmTablGet (tablptr, comploaddlt, comploaddltmax, (passnbr & 1))) != NULL)) {
      Gnum                comploadabsdlt;
      int                 partval;                /* Part of current vertex */
      Gnum                vertnum;
      Gnum                edgenum;

      partval = 0;                                /* Assume we point to gainlink0     */
      if (vexxptr->veloval >= 0) {                /* If in fact we point to gainlink1 */
        partval = 1;                              /* Then point to vertex structure   */
        vexxptr = (VgraphSeparateFmVertex *) ((byte *) vexxptr - ((byte *) &vexxptr->gainlink1 - (byte *) vexxptr));
      }
#ifdef SCOTCH_DEBUG_VGRAPH2
      if (vexxptr->partval != 2) {
        errorPrint ("vgraphSeparateFm: linked non-separator vertex (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
      gainTablDel (tablptr, &vexxptr->gainlink0); /* Remove it from table */
      gainTablDel (tablptr, &vexxptr->gainlink1);
      vexxptr->gainlink0.next = VGRAPHSEPAFMSTATESUCH; /* Mark it as used and avoid chaining */
      vexxptr->gainlink1.prev = &lockdat;         /* Lock it                                 */
      vexxptr->gainlink1.next = lockdat.next;
      lockdat.next->prev      = &vexxptr->gainlink1;
      lockdat.next            = &vexxptr->gainlink1;

      vertnum      = vexxptr->vertnum;            /* Get vertex number */
      compload2   += vexxptr->compgain[partval];
      comploaddlt -= (2 * partval - 1) * (vexxptr->compgain[partval] - 2 * vexxptr->veloval); /* TRICK: -veloval */

      if (vexxptr->mswpnum != mswpnum) {          /* If vertex data not yet recorded */
        vexxptr->mswpnum = mswpnum;
        savetab[savenbr].hashnum     = vexxptr - hashtab;
        savetab[savenbr].partval     = 2;
        savetab[savenbr].compgain[0] = vexxptr->compgain[0];
        savetab[savenbr].compgain[1] = vexxptr->compgain[1];
        savenbr ++;                               /* One more move recorded */
      }
      movenbr ++;                                 /* One more move done */

      sepaptr = NULL;                             /* No separator vertices to relink yet */
      for (edgenum = verttax[vertnum];            /* Update neighbors                    */
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;
        Gnum                hashnum;

        vertend = edgetax[edgenum];
        for (hashnum = (vertend * VGRAPHSEPAFMHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
          VgraphSeparateFmVertex *            vexxend; /* Pointer to neighbor of current vertex */

          vexxend = hashtab + hashnum;            /* Point to neighbor              */
          if (vexxend->vertnum == ~0) {           /* If neighbor does not exist yet */
            if (parttax[vertend] == partval)      /* If no use to create it         */
              break;                              /* Skip to next vertex            */
#ifdef SCOTCH_DEBUG_VGRAPH2
            if (parttax[vertend] != (1 - partval)) {
              errorPrint ("vgraphSeparateFm: undeclared separator vertex");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

            vexxend->vertnum        = vertend;    /* Set its number (TRICK: mswpnum assumed to be always -1) */
            vexxend->partval        = 1 - partval; /* Vertex will be in separator                            */
            vexxend->veloval        = - ((velotax != NULL) ? velotax[vertend] : 1);
            vexxend->gainlink0.next = VGRAPHSEPAFMSTATEFREE; /* Vertex will be linked */
            hashnbr ++;                           /* One more vertex in hash table    */
#ifdef SCOTCH_DEBUG_VGRAPH2
            if (hashnbr > hashmsk) {
              errorPrint ("vgraphSeparateFm: hash table overflow");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
          }
          if (vexxend->vertnum == vertend) {      /* If end vertex has been found       */
            if (vexxend->partval == 2) {          /* If already in separator or chained */
              if (vexxend->mswpnum != mswpnum) {  /* If vertex data not yet recorded    */
                vexxend->mswpnum = mswpnum;
                savetab[savenbr].hashnum     = hashnum;
                savetab[savenbr].partval     = 2;
                savetab[savenbr].compgain[0] = vexxend->compgain[0];
                savetab[savenbr].compgain[1] = vexxend->compgain[1];
                savenbr ++;                       /* One more move recorded */
              }
              vexxend->compgain[1 - partval] -= vexxptr->veloval; /* TRICK: -veloval                              */
              if (vexxend->gainlink0.next >= VGRAPHSEPAFMSTATELINK) { /* If vertex is linked                      */
                gainTablDel (tablptr, &vexxend->gainlink0); /* Unlink it temporarily                              */
                gainTablDel (tablptr, &vexxend->gainlink1); /* TRICK: gainlink1.next != NULL                      */
                vexxend->gainlink0.next = VGRAPHSEPAFMSTATEFREE; /* Mark separator vertex as temporarily unlinked */
                vexxend->gainlink0.prev = (GainLink *) sepaptr; /* Chain it for relinking                         */
                sepaptr                 = vexxend;
              }
              else if (vexxend->gainlink0.next == VGRAPHSEPAFMSTATEUSED) {
                vexxend->gainlink0.next = VGRAPHSEPAFMSTATESUCH; /* Mark separator vertex as chained-used */
                vexxend->gainlink0.prev = (GainLink *) sepaptr; /* Chain it for relinking                 */
                sepaptr                 = vexxend;
              }
            }
            else if (vexxend->partval == (1 - partval)) { /* Vertex is in other part */
              Gnum                edgeend;
              Gnum                compgainp;      /* Gain to be added to gain of part partval */

              if (vexxend->mswpnum != mswpnum) {  /* If vertex data not yet recorded */
                vexxend->mswpnum = mswpnum;
                savetab[savenbr].hashnum     = hashnum;
                savetab[savenbr].partval     = 1 - partval;
                savetab[savenbr].compgain[0] =    /* Vertex not in separator so gains are not relevant */
                savetab[savenbr].compgain[1] = 0;
                savenbr ++;                       /* One more move recorded */
              }

              vexxend->partval               = 2; /* Vertex will be in separator                       */
              vexxend->compgain[partval]     = vexxend->veloval; /* Moved vertex still in separator    */
              vexxend->compgain[1 - partval] = vexxend->veloval - vexxptr->veloval; /* TRICK: -veloval */

              for (edgeend = verttax[vertend], compgainp = 0;
                   edgeend < vendtax[vertend]; edgeend ++) {
                Gnum                vertent;
                Gnum                hashnum;

                vertent = edgetax[edgeend];
                for (hashnum = (vertent * VGRAPHSEPAFMHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
                  VgraphSeparateFmVertex *    vexxent; /* Pointer to neighbor of neighbor of current vertex */

                  vexxent = hashtab + hashnum;
                  if (vexxent->vertnum == ~0) {   /* If neighbor does not exist */
#ifdef SCOTCH_DEBUG_VGRAPH2
                    if (parttax[vertent] != (1 - partval)) {
                      errorPrint ("vgraphSeparateFm: broken separator (1)");
                      return     (1);
                    }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
                    compgainp += (velotax != NULL) ? velotax[vertent] : 1;
                    break;                        /* Skip to next vertex */
                  }
                  if (vexxent->vertnum == vertent) { /* If end vertex found */
#ifdef SCOTCH_DEBUG_VGRAPH2
                    if (vexxent->partval == partval) {
                      errorPrint ("vgraphSeparateFm: broken separator (2)");
                      return     (1);
                    }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
                    if (vexxent->partval == 2) {  /* If vertex is in separator (or is vexxptr) */
                      if (vexxent->mswpnum != mswpnum) { /* If vertex data not yet recorded    */
                        vexxent->mswpnum = mswpnum;
                        savetab[savenbr].hashnum     = hashnum;
                        savetab[savenbr].partval     = 2;
                        savetab[savenbr].compgain[0] = vexxent->compgain[0];
                        savetab[savenbr].compgain[1] = vexxent->compgain[1];
                        savenbr ++;               /* One more move recorded */
                      }

                      vexxent->compgain[partval] += vexxend->veloval; /* TRICK: -veloval                                  */
                      if (vexxent->gainlink0.next >= VGRAPHSEPAFMSTATELINK) { /* If not already chained                   */
                        gainTablDel (tablptr, &vexxent->gainlink0); /* Unlink it temporarily                              */
                        gainTablDel (tablptr, &vexxent->gainlink1); /* TRICK: gainlink1.next != NULL                      */
                        vexxent->gainlink0.next = VGRAPHSEPAFMSTATEFREE; /* Mark separator vertex as temporarily unlinked */
                        vexxent->gainlink0.prev = (GainLink *) sepaptr; /* Chain it                                       */
                        sepaptr                 = vexxent;
                      }
                      else if (vexxent->gainlink0.next == VGRAPHSEPAFMSTATEUSED) {
                        vexxent->gainlink0.next = VGRAPHSEPAFMSTATESUCH; /* Mark separator vertex as chained-used */
                        vexxent->gainlink0.prev = (GainLink *) sepaptr; /* Chain it for relinking                 */
                        sepaptr                 = vexxent;
                      }
                    }
                    else                          /* Vertex is in same part as vexxend */
                      compgainp -= vexxent->veloval; /* TRICK: -veloval                */
                    break;
                  }
                }
              }
              vexxend->compgain[partval] += compgainp;
              if (vexxend->gainlink0.next == VGRAPHSEPAFMSTATEUSED) /* If vertex was already used    */
                vexxend->gainlink0.next = VGRAPHSEPAFMSTATESUCH; /* Set it as separator-used-chained */
              vexxend->gainlink0.prev = (GainLink *) sepaptr; /* Chain it for relinking              */
              sepaptr                 = vexxend;
            }
            break;                                /* If in same part, ignore */
          }
        }
      }
      vexxptr->gainlink0.next = VGRAPHSEPAFMSTATEUSED; /* Mark it as used and not chained */
      vexxptr->partval = partval;                 /* Set vertex part last                 */
      while (sepaptr != NULL) {                   /* For all vertices in chain list       */
        vexxptr = sepaptr;                        /* Unlink vertex from list              */
        sepaptr = (VgraphSeparateFmVertex *) vexxptr->gainlink0.prev;
#ifdef SCOTCH_DEBUG_VGRAPH2
        if (vexxptr->partval != 2) {
          errorPrint ("vgraphSeparateFm: linked non-separator vertex (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

        if (vexxptr->gainlink0.next == VGRAPHSEPAFMSTATEFREE) { /* If vertex is not used */
          gainTablAdd (tablptr, &vexxptr->gainlink0, vexxptr->compgain[0]); /* Link it   */
          gainTablAdd (tablptr, &vexxptr->gainlink1, vexxptr->compgain[1]);
        }
        else {
          vexxptr->gainlink0.next = VGRAPHSEPAFMSTATEUSED;
          if (vexxptr->compgain[partval] == vexxptr->veloval) { /* If immediate gain                 */
            vexxptr->gainlink1.next->prev = vexxptr->gainlink1.prev; /* Remove vertex from lock list */
            vexxptr->gainlink1.prev->next = vexxptr->gainlink1.next;
            gainTablAdd (tablptr, &vexxptr->gainlink0, vexxptr->compgain[0]); /* Link it */
            gainTablAdd (tablptr, &vexxptr->gainlink1, vexxptr->compgain[1]);
          }
        }
      }
#ifdef SCOTCH_DEBUG_VGRAPH3
      if (vgraphSeparateFmCheck (grafptr, hashtab, hashmsk, compload2, comploaddlt) != 0) {
        errorPrint ("vgraphSeparateFm: internal error (3)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VGRAPH3 */

      if (hashnbr >= hashmax) {
        if (vgraphSeparateFmResize (&hashtab, &hashmax, &hashmsk, &savetab, savenbr, tablptr, &lockdat) != 0) {
          errorPrint ("vgraphSeparateFm: out of memory (2)");
          return     (1);
        }
      }

      comploadabsdltbst = abs (comploaddltbst);
      comploadabsdlt    = abs (comploaddlt);
      if ((comploadabsdlt    < comploaddltmat) || /* Record move only if it is within bounds */
          (comploadabsdltbst > comploaddltmat)) { /* Or if we have always been out of bounds */
        if (compload2 < compload2bst) {           /* If move improves the cost               */
          compload2bst   = compload2;             /* This move was effective                 */
          comploaddltbst = comploaddlt;
          movenbr  =
          savenbr  = 0;
          moveflag = 1;
          mswpnum ++;
        } else if (compload2 == compload2bst) {
          if (comploadabsdlt < comploadabsdltbst) {
            comploaddltbst = comploaddlt;         /* This move was effective */
            movenbr  =
            savenbr  = 0;
            moveflag = 1;
            mswpnum ++;
          }
          else if (comploadabsdlt == comploadabsdltbst) {
            comploaddltbst = comploaddlt;         /* Might be the opposite, so record */
            savenbr = 0;                          /* Forget backtracking              */
            mswpnum ++;
          }
        }
      }

      if (comploadabsdlt > comploaddltmax)        /* If an isolated vertex unbalanced the partition */
        comploaddltmax = comploadabsdlt;          /* Record that we degraded maximum load imbalance */
      else if (comploaddltmax > comploaddltmat) { /* Else if we must restrict distance bounds       */
        Gnum                comploaddlttmp;

        comploaddlttmp = comploaddltmax;          /* Save old working compdeltmax value           */
        comploaddltmax = MAX (comploaddltmat, comploadabsdlt); /* Restrict at most to the maximum */
        if ((comploadabsdltbst > comploaddltmat) &&  /* If we have never achieved balance yet     */
            (comploaddltmax < comploaddlttmp)) {  /* And if we have done something useful         */
          compload2bst   = compload2;             /* Then record best move done                   */
          comploaddltbst = comploaddlt;
          movenbr =                               /* Never set moveflag so as not to create an infinite loop */
          savenbr = 0;
          mswpnum ++;
        }
      }
    }
  } while ((moveflag != 0) &&                     /* As long as vertices are moved                          */
           (-- passnbr != 0));                    /* And we are allowed to loop (TRICK for negative values) */

  while (savenbr -- > 0) {                        /* Delete exceeding moves */
    Gnum                hashnum;
    int                 partval;

    hashnum = savetab[savenbr].hashnum;
    partval = savetab[savenbr].partval;
    hashtab[hashnum].partval = partval;           /* Restore vertex part only for update computation */
  }
  compload2    = compload2bst;                    /* Restore best separator parameters */
  comploaddlt  = comploaddltbst;
  compsize1add =                                  /* Variables for superscalar update */
  compsize1sub = 0;
  for (vexxptr = hashtab, fronnum = 0;            /* Build new frontier                */
       vexxptr < hashtab + (hashmax << 2); vexxptr ++) { /* From all vertices in table */
    Gnum                vertnum;

    vertnum = vexxptr->vertnum;
    if (vertnum != ~0) {                          /* If vertex slot is used     */
      int                 partval;                /* New part of current vertex */
      int                 partold;                /* Old part of current vertex */

      partval = vexxptr->partval;
      partold = parttax[vexxptr->vertnum];        /* Get old part value from array */
      if (partval != partold) {                   /* If vertex part changed        */
        parttax[vertnum] = partval;               /* Set new part value            */
        compsize1add += (partval & 1);            /* Superscalar update            */
        compsize1sub += (partold & 1);
      }
      if (partval == 2)                           /* If vertex belongs to cut */
        grafptr->frontab[fronnum ++] = vertnum;   /* Add vertex to frontier   */
    }
  }
  grafptr->compload[0] = ((grafptr->s.velosum - compload2) + comploaddlt) / 2;
  grafptr->compload[1] = ((grafptr->s.velosum - compload2) - comploaddlt) / 2;
  grafptr->compload[2] = compload2;
  grafptr->comploaddlt = comploaddlt;
  grafptr->compsize[1] = grafptr->compsize[1] + compsize1add - compsize1sub;
  grafptr->compsize[0] = grafptr->s.vertnbr - grafptr->compsize[1] - fronnum;
  grafptr->fronnbr     = fronnum;

#ifdef SCOTCH_DEBUG_VGRAPH2
  if (vgraphCheck (grafptr) != 0) {
    errorPrint ("vgraphSeparateFm: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VGRAPH2 */

  memFree      (hashtab);                         /* Free group leader */
  gainTablExit (tablptr);

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
vgraphSeparateFmResize (
VgraphSeparateFmVertex * restrict * hashtabptr,   /*+ Pointer to hash vertex table                   +*/
Gnum * const                        hashmaxptr,   /*+ Pointer to maximum number of elements in table +*/
Gnum * const                        hashmskptr,   /*+ Pointer to hash table mask                     +*/
VgraphSeparateFmSave * restrict *   savetabptr,   /*+ Pointer to move array                          +*/
const Gnum                          savenbr,      /*+ Current number of active slots in move array   +*/
GainTabl * const                    tablptr,      /*+ Gain table                                     +*/
GainLink * const                    lockptr)
{
  VgraphSeparateFmVertex * restrict hashtab;      /* Pointer to new hash table                    */
  VgraphSeparateFmSave *            savetab;      /* Pointer to new save array                    */
  VgraphSeparateFmSave *            saveold;      /* Pointer to translated old save array         */
  Gnum                              savenum;
  Gnum                              hashold;      /* Size of old hash table (half of new)         */
  Gnum                              hashsiz;
  Gnum                              hashmax;
  Gnum                              hashmsk;
  Gnum                              hashsta;      /* Start index of range of hash indices to move */
  Gnum                              hashend;      /* End index of range of hash indices to move   */
  Gnum                              hashnum;

  hashmax = *hashmaxptr << 1;                     /* Compute new sizes */
  hashold = *hashmaxptr << 2;
  hashsiz = *hashmaxptr << 3;
  hashmsk = hashsiz - 1;

  if (memReallocGroup ((void *) *hashtabptr,
                       &hashtab, (size_t) (hashsiz * sizeof (VgraphSeparateFmVertex)),
                       &savetab, (size_t) (hashsiz * sizeof (VgraphSeparateFmSave)), NULL) == NULL) {
    errorPrint ("vgraphSeparateFmResize: out of memory");
    return (1);
  }

  saveold = (VgraphSeparateFmSave *) ((byte *) hashtab + ((byte *) *savetabptr - (byte *) *hashtabptr));
  for (savenum = savenbr - 1; savenum >= 0; savenum --) { /* Move save array, in reverse order */
    savetab[savenum].compgain[1] = saveold[savenum].compgain[1];
    savetab[savenum].compgain[0] = saveold[savenum].compgain[0];
    savetab[savenum].partval     = saveold[savenum].partval;
    savetab[savenum].hashnum     = hashtab[saveold[savenum].hashnum].vertnum; /* Temporarily translate from hash index to number */
  }

  *hashtabptr = hashtab;
  *hashmaxptr = hashmax;
  *hashmskptr = hashmsk;
  *savetabptr = savetab;

  memSet (hashtab + hashold, ~0, hashold * sizeof (VgraphSeparateFmVertex));

  gainTablFree (tablptr);                         /* Reset gain table  */
  lockptr->next =                                 /* Rebuild lock list */
  lockptr->prev = lockptr;

  for (hashsta = hashold - 1; hashtab[hashsta].vertnum != ~0; hashsta --) ; /* Start index of first segment to reconsider is last empty slot */
  hashend = hashold;                              /* First segment to reconsider ends at the end of the old array                            */
  while (hashend != hashsta) {                    /* For each of the two segments to consider                                                */
    for (hashnum = hashsta; hashnum < hashend; hashnum ++) { /* Re-compute position of vertices in new table                                 */
      Gnum                        vertnum;

      vertnum = hashtab[hashnum].vertnum;
      if (vertnum != ~0) {                        /* If hash slot used */
        Gnum                        hashnew;

        for (hashnew = (vertnum * VGRAPHSEPAFMHASHPRIME) & hashmsk; ; hashnew = (hashnew + 1) & hashmsk) {
          if (hashnew == hashnum)                 /* If hash slot is the same */
            break;                                /* There is nothing to do   */
          if (hashtab[hashnew].vertnum == ~0) {   /* If new slot is empty     */
#ifdef SCOTCH_DEBUG_VGRAPH2
            if ((hashnew > hashnum) && (hashnew < hashend)) { /* If vertex is not moved either before its old position or after the end of the segment */
              errorPrint ("vgraphSeparateFmResize: internal error (1)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
            hashtab[hashnew] = hashtab[hashnum];  /* Copy data to new slot         */
            hashtab[hashnum].mswpnum = ~0;        /* TRICK: not tested at creation */
            hashtab[hashnum].vertnum = ~0;        /* Make old slot empty           */
            break;
          }
        }

        if (hashtab[hashnew].gainlink0.next >= VGRAPHSEPAFMSTATELINK) { /* If vertex was linked, re-link it */
          gainTablAdd (tablptr, &hashtab[hashnew].gainlink0, hashtab[hashnew].compgain[0]);
          gainTablAdd (tablptr, &hashtab[hashnew].gainlink1, hashtab[hashnew].compgain[1]);
        }
        else if (hashtab[hashnew].gainlink0.next == VGRAPHSEPAFMSTATEUSED) { /* Re-lock used vertices */
          hashtab[hashnew].gainlink1.prev = lockptr; /* Lock it */
          hashtab[hashnew].gainlink1.next = lockptr->next;
          lockptr->next->prev = &hashtab[hashnew].gainlink1;
          lockptr->next       = &hashtab[hashnew].gainlink1;
        }
      }
    }

    hashend = hashsta;                            /* End of second segment to consider is start of first one    */
    hashsta = 0;                                  /* Start of second segment is beginning of array              */
  }                                               /* After second segment, hashsta = hashend = 0 and loop stops */

  for (savenum = 0; savenum < savenbr; savenum ++) {
    Gnum                  vertnum;
    Gnum                  hashnum;

    vertnum = savetab[savenum].hashnum;           /* Get vertex number temporarily saved */
    for (hashnum = (vertnum * VGRAPHSEPAFMHASHPRIME) & hashmsk; hashtab[hashnum].vertnum != vertnum; hashnum = (hashnum + 1) & hashmsk) {
#ifdef SCOTCH_DEBUG_VGRAPH2
      if (hashtab[hashnum].vertnum == ~0) {
        errorPrint ("vgraphSeparateFmResize: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VGRAPH2 */
    }
    savetab[savenum].hashnum = hashnum;           /* Set new hash table index */
  }

  return (0);
}

/* This routine checks the consistency of
** the hash structures.
** It returns:
** - 0   : in case of success.
** - !0  : in case of error.
*/

#ifdef SCOTCH_DEBUG_VGRAPH3
static
int
vgraphSeparateFmCheck (
const Vgraph * restrict const                 grafptr,
const VgraphSeparateFmVertex * restrict const hashtab,
const Gnum                                    hashmsk,
const Gnum                                    compload2,
const Gnum                                    comploaddlt)
{
  Gnum                  hashnum;
  Gnum                  comploadtmp[3];

  const Gnum * restrict const       verttax = grafptr->s.verttax; /* Fast accesses */
  const Gnum * restrict const       vendtax = grafptr->s.vendtax;
  const Gnum * restrict const       velotax = grafptr->s.velotax;
  const Gnum * restrict const       edgetax = grafptr->s.edgetax;
  const GraphPart * restrict const  parttax = grafptr->parttax;

  comploadtmp[0] = grafptr->compload[0];
  comploadtmp[1] = grafptr->compload[1];
  comploadtmp[2] = grafptr->compload[2];
  for (hashnum = 0; hashnum <= hashmsk; hashnum ++) { /* For all vertex slots */
    Gnum                vertnum;
    int                 partval;

    vertnum = hashtab[hashnum].vertnum;
    if (vertnum == ~0)                            /* If unallocated slot */
      continue;                                   /* Skip to next slot   */

    if (hashtab[hashnum].veloval != - ((velotax == NULL) ? 1 : velotax[vertnum])) {
      errorPrint ("vgraphSeparateFmCheck: invalid vertex load (1)");
      return     (1);
    }
    partval = hashtab[hashnum].partval;
    if ((partval < 0) || (partval > 2)) {
      errorPrint ("vgraphSeparateFmCheck: invalid part value");
      return     (1);
    }

    if (partval != parttax[vertnum]) {
      comploadtmp[parttax[vertnum]] += hashtab[hashnum].veloval; /* TRICK: -veloval */
      comploadtmp[partval]          -= hashtab[hashnum].veloval;
    }

    if (partval < 2) {                            /* If not separator vertex */
      if (hashtab[hashnum].gainlink0.next >= VGRAPHSEPAFMSTATELINK) {
        errorPrint ("vgraphSeparateFmCheck: linked non-separator vertex");
        return     (1);
      }
    }
    else {                                        /* Separator vertex */
      Gnum                compload[3];
      Gnum                edgenum;

      if (hashtab[hashnum].gainlink0.next == VGRAPHSEPAFMSTATEFREE) {
        errorPrint ("vgraphSeparateFmCheck: free separator vertex");
        return     (1);
      }

      compload[0] =
      compload[1] =
      compload[2] = 0;
      for (edgenum = verttax[vertnum];            /* For all element neighbors */
           edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;
        Gnum                hashnum;
        int                 partend;
        Gnum                veloend;

        vertend = edgetax[edgenum];
        for (hashnum = (vertend * VGRAPHSEPAFMHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
          if (hashtab[hashnum].vertnum == vertend) { /* If end vertex found */
            partend = hashtab[hashnum].partval;
            veloend = hashtab[hashnum].veloval;

            if (veloend != - ((velotax == NULL) ? 1 : velotax[vertend])) {
              errorPrint ("vgraphSeparateFmCheck: invalid vertex load (2)");
              return     (1);
            }
            break;
          }
          if (hashtab[hashnum].vertnum == ~0) {     /* If element not present */
            partend = parttax[vertend];
            veloend = - ((velotax == NULL) ? 1 : velotax[vertend]);
            break;
          }
        }
        compload[partend] += veloend;
      }

      if ((hashtab[hashnum].compgain[0] != (hashtab[hashnum].veloval - compload[1])) ||
          (hashtab[hashnum].compgain[1] != (hashtab[hashnum].veloval - compload[0]))) {
        errorPrint ("vgraphSeparateFmCheck: invalid vertex gains");
        return     (1);
      }
    }
  }
  if (compload2 != comploadtmp[2]) {
    errorPrint ("vgraphSeparateFmCheck: invalid frontier load");
    return     (1);
  }
  if (comploaddlt != (comploadtmp[0] - comploadtmp[1])) {
    errorPrint ("vgraphSeparateFmCheck: invalid separator balance");
    return     (1);
  }

  return (0);
}
#endif /* SCOTCH_DEBUG_VGRAPH3 */
