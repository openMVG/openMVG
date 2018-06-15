/* Copyright 2007-2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_fm.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module is the improved Fiduccia-   **/
/**                Mattheyses refinement routine for the   **/
/**                vertex overlapped graph partitioning.   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to     24 dec 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define WGRAPH_PART_FM

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "wgraph.h"
#include "wgraph_part_gg.h"
#include "wgraph_part_fm.h"

/*
**  The static variables.
*/

static const Gnum           wgraphpartfmloadone = 1;


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
WgraphPartFmLink *
wgraphPartFmTablGet (
Wgraph * restrict const     grafptr,          /*+ Active graph      +*/
WgraphPartFmVertex *        hashtab,
GainTabl * const            tablptr,          /*+ Gain table        +*/
Gnum * restrict             diffload)
{
  Gnum                        gainbest;       /* Separator gain of best link       */
  Gnum                        bestdiffload;
  const GainEntr *            tablbest;       /* Gain table entry of best link     */
  WgraphPartFmLink *          linkptr;        /* Pointer to current gain link      */
  WgraphPartFmLink *          linkbest;       /* Pointer to best link found        */

  linkbest = NULL;                            /* Assume no candidate vertex found yet */
  tablbest = tablptr->tend;
  gainbest = GAINMAX;

  for (linkptr = (WgraphPartFmLink *) gainTablFrst (tablptr);          /* Select candidate vertices */
       (linkptr != NULL) && (linkptr->gainlink.tabl <= tablbest);
       linkptr = (WgraphPartFmLink *) gainTablNext (tablptr, (GainLink *) linkptr)) {
    Gnum                      vertpart;       /* Part of current vertex */
    Gnum                      gaincur;        /* Separator gain of current link    */

    vertpart = linkptr->partval;
    gaincur  = linkptr->gain;                 /* Get separator gain and vertex balance */

    if (linkptr->minloadpartval == -1) {
      return (linkptr);                       /* Return best link found */
    }

    if ((gaincur < gainbest) ||               /* And if it gives better gain than gain max */
        ((gaincur == gainbest) &&             /* Or is in preferred part     */
         ((grafptr->compload[vertpart] - grafptr->compload[linkptr->minloadpartval]) <= bestdiffload))) {
      linkbest     = linkptr;                 /* Select it */
      gainbest     = gaincur;
      tablbest     = linkptr->gainlink.tabl;
      bestdiffload = grafptr->compload[vertpart] - grafptr->compload[linkptr->minloadpartval];
    }
  }

  *diffload += bestdiffload;
  return (linkbest);                          /* Return best link found */
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

int
wgraphPartFm (
Wgraph * restrict const         grafptr,    /*+ Active graph      +*/
const WgraphPartFmParam * const paraptr)    /*+ Method parameters +*/
{
  int                             passnbr;    /* Maximum number of passes to go       */
  int                             moveflag;   /* Flag set if useful moves made        */
  Gnum                            fronnum;    /* Current index of frontier vertex     */
  Gnum                            partval;
  Gnum                            partnbr;
  Gnum                            vertnum;
  Gnum                            hashnum;
  Gnum                            velosum;
  Gnum                            minload;
  Gnum                            maxload;
  Gnum                            frlobst;
  Gnum                            linknbr;
  Gnum                            savenum;
  Gnum                            hashmax;    /* Maximum number of elements in table  */
  Gnum                            hashsiz;    /* Size of hash and save tables         */
  Gnum                            hashmsk;    /* Mask for access to hash table        */
  Gnum                            hashnbr;    /* Number of elements in hash table     */
  Gnum                            velomsk;
  Gnum * restrict                 movetab;
  GainTabl * restrict             tablptr;    /* Pointer to gain tables               */
  const Gnum * restrict           velobax;    /* Data for handling of optional arrays */
  WgraphPartFmVertex *            hashtab;    /* Hash vertex table                    */
  WgraphPartFmVertex *            vexxptr;
  WgraphPartFmVertex *            vertlist;
  WgraphPartFmVertex *            locklist;
  WgraphPartFmSave * restrict     savetab;
  WgraphPartFmLink * restrict     linktab;
  WgraphPartFmPartList * restrict partlist;
  WgraphPartFmPartList * restrict partlistptr;
  WgraphPartFmPartList * restrict partlistloadptr;

  partnbr = grafptr->partnbr;
  for (partval = 0, velosum = 0; partval < partnbr; partval ++)
    velosum += grafptr->compload[partval];

  if (grafptr->s.velotax == NULL) {               /* Set accesses to optional arrays             */
    velobax = &wgraphpartfmloadone;  /* In case vertices not weighted (least often) */
    velomsk = 0;
  }
  else {
    velobax = grafptr->s.velotax;
    velomsk = ~((Gnum) 0);
  }

  minload = velosum * (1. - paraptr->deltrat) / (float) (partnbr);
  maxload = velosum * (1. + paraptr->deltrat) / (float) (partnbr);

  if (grafptr->fronnbr == 0) {                    /* If no frontier defined    */
    for (partval = 0; partval < partnbr; partval ++) {
      if (grafptr->compload[partval] < minload      /* If balance is not achieved */
          || grafptr->compload[partval] > maxload) { /* Imbalance must be fought  */
        WgraphPartGgParam         paradat;

        paradat.passnbr = 5;                      /* Use a standard algorithm */

        wgraphPartGg (grafptr, &paradat);

        if (grafptr->fronnbr == 0)                /* If new partition has no frontier */
          return (0);                             /* This algorithm is still useless  */

      }
      return (0);                                 /* This algorithm is still useless  */
    }
  }

  hashnbr = 16 * (grafptr->fronnbr + paraptr->movenbr + grafptr->s.degrmax) + 1;

  if (hashnbr > grafptr->s.vertnbr)
    hashnbr = 4 * grafptr->s.vertnbr;

  hashnbr *= 2;

  for (hashsiz = 512; hashsiz < hashnbr; hashsiz <<= 1) ; /* Get upper power of two */
  hashmsk = hashsiz - 1;
  hashmax = hashsiz >> 2;                         /* Use hash table at 1/4 of its capacity */
  hashnbr = 0;

  if (((tablptr = gainTablInit (GAINMAX, WGRAPHSEPAFMGAINBITS)) == NULL) || /* Use logarithmic array only */
      (memAllocGroup ((void **) (void *)
                      &hashtab, (size_t) (hashsiz * sizeof (WgraphPartFmVertex)),
                      &linktab, (size_t) (hashsiz * sizeof (WgraphPartFmLink)),
                      &partlist,(size_t) (partnbr * sizeof (WgraphPartFmPartList)),
                      &savetab, (size_t) (hashsiz * sizeof (WgraphPartFmSave)),
                      &movetab, (size_t) (hashsiz * sizeof (Gnum)), NULL) == NULL)) {
    errorPrint ("wgraphPartFm: out of memory (1)");

    if (tablptr != NULL)
      gainTablExit (tablptr);

    return (1);
  }
  memSet (hashtab, ~0, hashsiz * sizeof (WgraphPartFmVertex));
  memSet (linktab, ~0, hashsiz * sizeof (WgraphPartFmLink));

  linknbr  = 0;
  locklist = NULL;
  passnbr  = paraptr->passnbr;
  frlobst  = grafptr->fronload;
  for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {       /* Set initial gains */
    Gnum                        edgenum;
    Gnum                        compgain;
    Gnum                        minloadpartval;
    Gnum                        minloadpartload;
    Gnum                        secondloadpartval;
    WgraphPartFmLink * restrict linklist;
    WgraphPartFmLink * restrict linkptr;

    minloadpartval    =
    secondloadpartval = -1;
    partlistptr       = NULL;
    vertnum           = grafptr->frontab[fronnum];

    compgain = - velobax[vertnum & velomsk];    /* Gain initialised as negative value for the frontier */
    memSet (partlist, 0, partnbr * sizeof (WgraphPartFmPartList));
    for (edgenum = grafptr->s.verttax[vertnum];
         edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
      Gnum                              vertnum2;
      WgraphPartFmVertex *              vexxptr;

      vertnum2 = grafptr->s.edgetax[edgenum];
      for (hashnum = (vertnum2 * WGRAPHSEPAFMHASHPRIME) & hashmsk;
           (hashtab[hashnum].vertnum != vertnum2) && (hashtab[hashnum].vertnum != ~0);
           hashnum = (hashnum + 1) & hashmsk);

      vexxptr = hashtab + hashnum;
      if (vexxptr->vertnum == ~0) {          /* If vertex not found add it in the hash table */
        vexxptr->vertnum = vertnum2;
        vexxptr->partval = grafptr->parttax[vertnum2];
        vexxptr->linklist = NULL;
        hashnbr ++;
        if (hashnbr >= hashmax) {
          if (wgraphPartFmResize () != 0) {
            errorPrint ("wgraphPartFm: out of memory (2)");
            return     (1);
          }
        }
      }
      if (vexxptr->partval != -1) {                      /* If its part is not the separator */
        if (partlist[vexxptr->partval].gain == 0) {      /* and not yet linked               */
          partlist[vexxptr->partval].prev = partlistptr; /* link it                          */
          partlistptr = partlist + vexxptr->partval;

          if (minloadpartval == -1 ||
              minloadpartload > grafptr->compload[vexxptr->partval]) {
            secondloadpartval = minloadpartval;
            minloadpartval    = vexxptr->partval;
            minloadpartload   = grafptr->compload[minloadpartval];
          }
	  else if (secondloadpartval != vexxptr->partval) {
            secondloadpartval = vexxptr->partval;
          }
        }
        partlist[vexxptr->partval].gain -= velobax[vertnum2 & velomsk]; /* Store the gain of this vertex move for this part */
        compgain += velobax[vertnum2 & velomsk];                        /* Store the global gain of this vertex move        */
      }
    }

    for (hashnum = (vertnum * WGRAPHSEPAFMHASHPRIME) & hashmsk;
         hashtab[hashnum].vertnum != vertnum && hashtab[hashnum].vertnum != ~0;
         hashnum = (hashnum + 1) & hashmsk);
    if (hashtab[hashnum].vertnum == ~0) {         /* If vertex not found      */
      hashtab[hashnum].vertnum = vertnum;         /* Add it in the hash table */
      hashtab[hashnum].partval = -1;
      hashtab[hashnum].linklist = NULL;
      hashnbr ++;

      if (hashnbr >= hashmax) {
        if (wgraphPartFmResize () != 0) {
          errorPrint ("wgraphPartFm: out of memory (2)");
          return     (1);
        }
      }
    }

    hashtab[hashnum].linklist = (partlistptr != NULL) /* If selected vertex is not isolated         */
                                ? linktab + linknbr /* Then first link will next element in linktab */
                                : NULL;           /* Else, no link                                  */
    linklist =
    linkptr  = NULL;                              /* Assume list is empty */
    while (partlistptr != NULL) {                 /* For each linked part */
      partval          = partlistptr - partlist;
      partlistptr      = partlistptr->prev;
      linkptr          = linktab + linknbr;       /* Create link at the end of linktab */
      linkptr->partval = partval;
      linkptr->hashnum = hashnum;
      linkptr->gain    = compgain + partlist[partval].gain;

      if (partval != minloadpartval)
        linkptr->minloadpartval = minloadpartval;
      else
        linkptr->minloadpartval = secondloadpartval;
      linkptr->minloadpartload = minloadpartload;

      if (linklist != NULL)                       /* Add link to the list */
        linklist->next = linkptr;
      linklist = linkptr;

      gainTablAdd (tablptr, (GainLink *) (linktab + linknbr), linkptr->gain); /* add the link in the gain table */
      linknbr ++;

      partlist[partval].prev = NULL;
      partlist[partval].gain = 0;
    }
    if (linkptr != NULL)
      linkptr->next = NULL;                         /* close the end of the list */
  }

  movetab[0] = -1;

  do {                                            /* As long as there is improvement */
    Gnum                        partval2;
    Gnum                        movenbr;          /* Number of uneffective moves done */
    Gnum                        savenbr;          /* Position of current move         */
    Gnum                        diffload;
    WgraphPartFmLink * restrict linkptr;

    movenbr  = 0;                                 /* No uneffective moves yet */
    savenbr  = 0;                                 /* No recorded moves yet    */
    moveflag = 0;                                 /* No moves to date         */
    diffload = 0;

    while ((movenbr < paraptr->movenbr) &&        /* As long as we can find effective vertices */
           ((linkptr = wgraphPartFmTablGet (grafptr, hashtab, tablptr, &diffload)) != NULL)) {
      Gnum                          edgenum;
      WgraphPartFmVertex *          vexxptr2; /* Pointer to current vertex */
      WgraphPartFmVertex *          vertlist2;

      movenbr ++;
      minload  = velosum * (1. - paraptr->deltrat) / (float) (partnbr);
      maxload  = velosum * (1. + paraptr->deltrat) / (float) (partnbr);
      partval  = linkptr->partval;                 /* Get data from the selected link */
      hashnum  = linkptr->hashnum;
      vexxptr2 = hashtab + hashnum;
      vertnum  = vexxptr2->vertnum;

      vertlist        = NULL;                      /* and a list of vertex             */
      partlistptr     = NULL;                      /* initialise an empty list of part */
      partlistloadptr = NULL;

      vexxptr2->linked   = 1;                /* link the vertex */
      vexxptr2->partval2 = partval;        /* the vertex will move to the part partval */
      vexxptr2->prev     = vertlist;
      vertlist           = hashtab + hashnum;

      for (edgenum = grafptr->s.verttax[vertnum]; /* for neighbours vertices */
           edgenum < grafptr->s.vendtax[vertnum]; edgenum ++) {
        WgraphPartFmVertex *            vexxptr3;    /* Pointer to current vertex */
        Gnum                                vertnum2;
        Gnum                                edgenum2;

        vertnum2 = grafptr->s.edgetax[edgenum];

        for (hashnum = (vertnum2 * WGRAPHSEPAFMHASHPRIME) & hashmsk; /* search the vertex in the hash table */
             hashtab[hashnum].vertnum != vertnum2 && hashtab[hashnum].vertnum != ~0;
             hashnum = (hashnum + 1) & hashmsk);
        vexxptr3 = hashtab + hashnum;

        if (vexxptr3->vertnum == ~0) {     /* if vertex not found      */
          vexxptr3->vertnum = vertnum2;    /* add it in the hash table */
          vexxptr3->partval = grafptr->parttax[vertnum2];
          vexxptr3->linklist = NULL;
          hashnbr ++;
          if (hashnbr >= hashmax) {
            if (wgraphPartFmResize () != 0) {
              errorPrint ("wgraphPartFm: out of memory (2)");
              return     (1);
            }
          }
        }

        if (vexxptr3->partval == -1) {     /* if vertex is in the separator */
          if (vexxptr3->linked == ~0) {    /* link the vertex               */
            vexxptr3->linked = 1;
            vexxptr3->prev = vertlist;
            vertlist = hashtab + hashnum;
            vexxptr3->partval2 = -1;       /* the vertex will stay in the separator */
          }
        }
        else if ((vexxptr3->partval != partval)) { /* or if it's not in the part of the selected vertex */
          if (vexxptr3->linked == ~0) {            /* link the vertex                                   */
            vexxptr3->linked = 1;
            vexxptr3->prev = vertlist;
            vertlist = hashtab + hashnum;
            vexxptr3->partval2 = -1;       /* the vertex will move to the separator */
          }

          for (edgenum2 = grafptr->s.verttax[vertnum2]; /* for neighbours vertices */
               edgenum2 < grafptr->s.vendtax[vertnum2]; edgenum2 ++) {
            WgraphPartFmVertex *            vexxptr4; /* Pointer to current vertex */
            Gnum                                vertnum3;

            vertnum3 = grafptr->s.edgetax[edgenum2];

            for (hashnum = (vertnum3 * WGRAPHSEPAFMHASHPRIME) & hashmsk; /* search the vertex in the hash table */
                 hashtab[hashnum].vertnum != vertnum3 && hashtab[hashnum].vertnum != ~0;
                 hashnum = (hashnum + 1) & hashmsk);
            vexxptr4 = hashtab + hashnum;

            if (vexxptr4->vertnum == ~0) {  /* if vertex not found      */
              vexxptr4->vertnum = vertnum3; /* add it in the hash table */
              vexxptr4->partval = grafptr->parttax[vertnum3];
              vexxptr4->linklist = NULL;
              hashnbr ++;
              if (hashnbr >= hashmax) {
                if (wgraphPartFmResize () != 0) {
                  errorPrint ("wgraphPartFm: out of memory (2)");
                  return     (1);
                }
              }
            }

            if (vexxptr4->partval == -1) {  /* if vertex in separator */
              if (vexxptr4->linked == ~0) { /* link the vertex        */
                vexxptr4->linked = 1;
                vexxptr4->prev = vertlist;
                vexxptr4->partval2 = -1;
                vertlist = hashtab + hashnum;
              }
            }
          }
        }
      }
      vertlist2 = vertlist;
      grafptr->fronload -= velobax[vertnum & velomsk]; /* Decrease the frontier load (due to the selected vertex) */
      grafptr->fronnbr  --;                            /* Decrease the frontier size (due to the selected vertex) */

      while (vertlist != NULL) {                       /* For each vertex in the list of linked vertices */
        Gnum                    newpart;               /* Move vertex from part parval to part partval2  */
        Gnum                    oldpart;
        Gnum                    vertnum2;

        vertnum2 = vertlist->vertnum;
        oldpart  = vertlist->partval;
        newpart  = vertlist->partval2;
        if (oldpart != newpart)
          moveflag = 1;

        if ((newpart == -1) && (oldpart != -1)) {       /* If the vertex will move into separator */
          grafptr->fronload += velobax[vertnum2 & velomsk]; /* Increase the frontier load             */
          grafptr->fronnbr  ++;                             /* Increase the frontier size */
        }

        if (oldpart != -1) {
          grafptr->compload[oldpart] -= velobax[vertnum2 & velomsk]; /* Decrease the load of the selected part */
          grafptr->compsize[oldpart] --;                             /* Decrease the size of the selected part */

          partlist[oldpart].loadgain -= velobax[vertnum2 & velomsk];
          partlist[oldpart].sizegain --;
          if (partlist[oldpart].isinloadlist != 1) {
            partlist[oldpart].loadprev = partlistloadptr;
            partlistloadptr = partlist + oldpart;
            partlist[oldpart].isinloadlist = 1;
          }
        }
        if (newpart != -1) {
          grafptr->compload[newpart] += velobax[vertnum2 & velomsk]; /* increase the load of the selected part */
          grafptr->compsize[newpart] ++;           /* increase the size of the selected part */
          partlist[newpart].loadgain += velobax[vertnum2 & velomsk];
          partlist[newpart].sizegain ++;
          if (partlist[newpart].isinloadlist != 1) {
            partlist[newpart].loadprev = partlistloadptr;
            partlistloadptr = partlist + newpart;
            partlist[newpart].isinloadlist = 1;
          }
        }

        vertlist->partval2 = oldpart;             /* exchange the old and the new parts */
        vertlist->partval = newpart;
        if (savenbr >= hashsiz) {
          if (wgraphPartFmResize () != 0) {
            errorPrint ("wgraphPartFm: out of memory (2)");
            return     (1);
          }
        }

        savetab[savenbr].type = WGRAPHSEPAFMSAVEMOVE; /* save the move */
        savetab[savenbr].u.movedata.hashnum = vertlist - hashtab;
        savetab[savenbr].u.movedata.partval = vertlist->partval2;
        movetab[movenbr] = savenbr;
        savenbr ++;

        vertlist = vertlist->prev;
      }

      while (vertlist2 != NULL) {                 /* for each vertex in the same list of linked vertices */
        Gnum                                vertnum2;

        vertnum2 = vertlist2->vertnum;

        if (vertlist2->partval2 == -1) {          /* If vertex was in separator */
          WgraphPartFmLink * restrict linkptr2;

          for (linkptr2 = vertlist2->linklist;    /* For each link of the vertex */
               linkptr2 != NULL; linkptr2 = linkptr2->next) {
            if (linkptr2->gainlink.next != NULL) {
              grafptr->compload[linkptr2->partval] -= velobax[vertnum2 & velomsk]; /* decrease the load of this part */
              grafptr->compsize[linkptr2->partval]--; /* decrease the size of this part                              */
              partlist[linkptr2->partval].loadgain -= velobax[vertnum2 & velomsk];
              partlist[linkptr2->partval].sizegain --;

              if (partlist[linkptr2->partval].isinloadlist != 1) {
                partlist[linkptr2->partval].loadprev = partlistloadptr;
                partlistloadptr = partlist + linkptr2->partval;
                partlist[linkptr2->partval].isinloadlist = 1;
              }
            }

            if (linkptr2->gainlink.next != NULL) {
              if (linkptr2->gainlink.next != (GainLink *) 1) /* If link is in gain table */
                gainTablDel (tablptr, (GainLink *) linkptr2); /* Remove link from table  */

              if (savenbr >= hashsiz) {
                if (wgraphPartFmResize () != 0) { /* TODO: Check if table resizing changes linkptr2 ? */
                  errorPrint ("wgraphPartFm: out of memory (2)");
                  return     (1);
                }
              }

              savetab[savenbr].type = WGRAPHSEPAFMSAVELINKDEL; /* Save link removal from table */
              savetab[savenbr].u.linkdata.linknum = linkptr2 - linktab;
              savetab[savenbr].u.linkdata.gain = linkptr2->gain;
              movetab[movenbr] = savenbr;
              savenbr ++;
            }

            linkptr2->gainlink.next = NULL;
          }
        }

        if (vertlist2->partval == -1) {           /* if vertex move (or stay) in the separator */
          Gnum                        gain;
          Gnum                        minloadpartval;
          Gnum                        minloadpartload;
          Gnum                        secondloadpartval;
          WgraphPartFmLink * restrict linklist;
          WgraphPartFmLink * restrict linkptr2;

          vertnum2 = vertlist2->vertnum;
          partlistptr = NULL;
          gain = -1;
          minloadpartval = -1;
          secondloadpartval = -1;

          for (hashnum = (vertnum2 * WGRAPHSEPAFMHASHPRIME) & hashmsk; /* search the vertex in the hash table */
               hashtab[hashnum].vertnum != vertnum2;
               hashnum = (hashnum + 1) & hashmsk);

          vexxptr = hashtab + hashnum;

          for (edgenum = grafptr->s.verttax[vertnum2]; /* recompute the gain */
               edgenum < grafptr->s.vendtax[vertnum2]; edgenum ++) {
            WgraphPartFmVertex *            vexxptr4; /* Pointer to current vertex */
            Gnum                                vertnum3;

            vertnum3 = grafptr->s.edgetax[edgenum];

            for (hashnum = (vertnum3 * WGRAPHSEPAFMHASHPRIME) & hashmsk; /* search the vertex in the hash table */
                 hashtab[hashnum].vertnum != vertnum3 && hashtab[hashnum].vertnum != ~0;
                 hashnum = (hashnum + 1) & hashmsk);
            vexxptr4 = hashtab + hashnum;

            if (vexxptr4->vertnum == ~0) {  /* if vertex not found      */
              vexxptr4->vertnum = vertnum3; /* add it in the hash table */
              vexxptr4->partval = grafptr->parttax[vertnum3];
              vexxptr4->linklist = NULL;
              hashnbr ++;
              if (hashnbr >= hashmax) {
                if (wgraphPartFmResize () != 0) {
                  errorPrint ("wgraphPartFm: out of memory (2)");
                  return     (1);
                }
              }
            }
            if (vexxptr4->partval != -1) {                      /* if part is not the separator */
              if (partlist[vexxptr4->partval].gain == 0) {      /* and not yet linked           */
                partlist[vexxptr4->partval].prev = partlistptr; /* link it                      */
                partlistptr = partlist + vexxptr4->partval;

                if (minloadpartval == -1 ||
                    minloadpartload > grafptr->compload[vexxptr4->partval]) {
                  secondloadpartval = minloadpartval;
                  minloadpartval = vexxptr4->partval;
                  minloadpartload = grafptr->compload[minloadpartval];
                } else if (secondloadpartval != vexxptr4->partval) {
                  secondloadpartval = vexxptr4->partval;
                }
              }

              partlist[vexxptr4->partval].gain --;
              gain ++;
            }
          }

          for (linkptr2 = vexxptr->linklist, linklist = NULL; /* For each vertex link in list */
               linkptr2 != NULL; linkptr2 = linkptr2->next) {
            linklist = linkptr2;
            partval2 = linkptr2->partval;

            if (partlist[partval2].gain != 0) {   /* If part is linked */
              hashtab[linkptr2->hashnum].vertnum = vertnum2;
              linkptr2->partval = partval2;
              linkptr2->gain = gain + partlist[partval2].gain;

              if (partval2 != minloadpartval)
                linkptr2->minloadpartval = minloadpartval;
              else
                linkptr2->minloadpartval = secondloadpartval;
              linkptr2->minloadpartload = minloadpartload;

              if (hashtab[linkptr2->hashnum].lockprev == (WgraphPartFmVertex *) ~0) { /* If vertex is not locked */
                if (savenbr >= hashsiz) {
                  if (wgraphPartFmResize () != 0) {
                    errorPrint ("wgraphPartFm: out of memory (2)");
                    return     (1);
                  }
                }
                savetab[savenbr].type               = WGRAPHSEPAFMSAVELINKADD;
                savetab[savenbr].u.linkdata.linknum = linkptr2 - linktab;
                savetab[savenbr].u.linkdata.gain    = linkptr2->gain;
                movetab[movenbr]                    = savenbr;
                savenbr ++;

                gainTablAdd(tablptr, (GainLink *) linkptr2, linkptr2->gain); /* add the link in the gain table of the part */
              }
              else {
                if (savenbr >= hashsiz) {
                  if (wgraphPartFmResize () != 0) {
                    errorPrint ("wgraphPartFm: out of memory (2)");
                    return     (1);
                  }
                }

                savetab[savenbr].type               = WGRAPHSEPAFMSAVELINKADD;
                savetab[savenbr].u.linkdata.linknum = linkptr2 - linktab;
                savetab[savenbr].u.linkdata.gain    = linkptr2->gain;
                movetab[movenbr]                    = savenbr;
                savenbr ++;

                linkptr2->gainlink.next = (GainLink *) 1;
              }

              partlist[partval2].gain = 0;
              grafptr->compload[linkptr2->partval] += velobax[hashtab[linkptr2->hashnum].vertnum & velomsk]; /* increase the load of this part */
              grafptr->compsize[linkptr2->partval] ++; /* increase the size of this part */
              partlist[linkptr2->partval].loadgain += velobax[hashtab[linkptr2->hashnum].vertnum & velomsk];
              partlist[linkptr2->partval].sizegain ++;

              if (partlist[linkptr2->partval].isinloadlist != 1) {
                partlist[linkptr2->partval].isinloadlist = 1;
                partlist[linkptr2->partval].loadprev     = partlistloadptr;
                partlistloadptr                          = partlist + linkptr2->partval;
              }
            }
          }

          while (partlistptr != NULL) {           /* for each part in the linked list */
            partval2 = partlistptr - partlist;
            if (partlist[partval2].gain != 0) {
              WgraphPartFmLink * restrict linkptr3;

              linkptr3 = linktab + linknbr;
              if (linklist != NULL)               /* If vertex has a list of link */
                linklist->next = linkptr3;        /* Add link to the list         */
              else
                vexxptr->linklist = linkptr3;     /* Else create the list */

              linknbr ++;

              linkptr3->hashnum                  = vexxptr - hashtab;
              hashtab[linkptr3->hashnum].vertnum = vertnum2;
              linkptr3->partval                  = partval2;
              linkptr3->gain                     = gain + partlist[partval2].gain;

              if (partval2 != minloadpartval)
                linkptr3->minloadpartval = minloadpartval;
              else
                linkptr3->minloadpartval = secondloadpartval;
              linkptr3->minloadpartload = minloadpartload;

              if (vexxptr->lockprev == (WgraphPartFmVertex *) ~0/*  || linkptr3->gain == -1 */)
                gainTablAdd(tablptr, (GainLink *) linkptr3, linkptr3->gain); /* add the link in the gain table of the part */
              else
                linkptr3->gainlink.next = (GainLink *) 1;

              if (savenbr >= hashsiz) {
                if (wgraphPartFmResize () != 0) {
                  errorPrint ("wgraphPartFm: out of memory (2)");
                  return     (1);
                }
              }
              savetab[savenbr].type               = WGRAPHSEPAFMSAVELINKADD;
              savetab[savenbr].u.linkdata.linknum = linkptr3 - linktab;
              savetab[savenbr].u.linkdata.gain    = linkptr3->gain;
              movetab[movenbr]                    = savenbr;
              savenbr ++;

              linklist = linkptr3;

              grafptr->compload[partval2] += velobax[vertnum2 & velomsk];
              partlist[partval2].loadgain += velobax[vertnum2 & velomsk];
              grafptr->compsize[partval2] ++;
              partlist[partval2].sizegain ++;

              if (partlist[partval2].isinloadlist != 1) { /* link the part in the load list */
                partlist[partval2].isinloadlist = 1;
                partlist[partval2].loadprev     = partlistloadptr;
                partlistloadptr                 = partlist + partval2;
              }

              partlist[partval2].gain = 0;
            }

            partlistptr = partlistptr->prev;
            partlist[partval2].prev = NULL;
          }
          if (linklist != NULL)
	    linklist->next = NULL;
        }
        vertlist2->linked = ~0;
        vertlist2 = vertlist2->prev;
      }

      while (partlistloadptr != NULL) {           /* For each part in the load list */
        if (partlistloadptr->loadgain != 0) {
          if (savenbr >= hashsiz) {
            if (wgraphPartFmResize () != 0) {
              errorPrint ("wgraphPartFm: out of memory (2)");
              return     (1);
            }
          }
          savetab[savenbr].type                = WGRAPHSEPAFMSAVELOAD; /* Save load variation */
          savetab[savenbr].u.loaddata.partval  = partlistloadptr - partlist;
          savetab[savenbr].u.loaddata.loaddiff = partlistloadptr->loadgain;
          savetab[savenbr].u.loaddata.sizediff = partlistloadptr->sizegain;
          movetab[movenbr]                     = savenbr;
          savenbr ++;
          velosum += partlistloadptr->loadgain;
        }
        partlistloadptr->loadgain     =
        partlistloadptr->sizegain     =
        partlistloadptr->isinloadlist = 0;
        partlistloadptr               = partlistloadptr->loadprev;
      }

      vexxptr2->lockprev = locklist;              /* Lock the selected vertex */
      locklist           = vexxptr2;
      if (grafptr->fronload < frlobst) {
        movenbr  =
        savenbr  =
        diffload = 0;
        moveflag = 1;
        frlobst  = grafptr->fronload;
      }
      else if ((grafptr->fronload == frlobst) && (diffload < -1)) {
        movenbr  =
        savenbr  =
        diffload = 0;
        moveflag = 1;
      }
      else if (linkptr->minloadpartval == -1) { 
	movenbr  =
	savenbr  =
	diffload = 0;
      }
      else if ((grafptr->compload[linkptr->partval] < minload) ||
	       (grafptr->compload[linkptr->partval] < grafptr->compload[linkptr->minloadpartval])) { 
	movenbr  =
	savenbr  =
	diffload = 0;
      }
    }

    while (locklist != NULL) {                    /* For each locked vertex */
      WgraphPartFmLink * restrict linkptr2;

      vexxptr           = locklist;
      locklist          = locklist->lockprev;     /* Unlock it */
      vexxptr->lockprev = (WgraphPartFmVertex *) ~0;

      for (linkptr2 = vexxptr->linklist; linkptr2 != NULL; linkptr2 = linkptr2->next) {
        if (linkptr2->gainlink.next == (GainLink *) 1)
          gainTablAdd(tablptr, (GainLink *) linkptr2, linkptr2->gain); /* Add link to part gain table */
      }
    }
    locklist = NULL;

    for ( ; movenbr > 0; movenbr --) {             /* For each move to undo */
      for (savenum = movetab[movenbr]; savenum > movetab[movenbr - 1]; savenum --) {
        WgraphPartFmLink * restrict linkptr2;

        switch (savetab[savenum].type) {
          case WGRAPHSEPAFMSAVEMOVE :
            if (savetab[savenum].u.movedata.partval == -1)
              grafptr->fronload += velobax[hashtab[savetab[savenum].u.movedata.hashnum].vertnum & velomsk];
            if (hashtab[savetab[savenum].u.movedata.hashnum].partval == -1)
              grafptr->fronload -= velobax[hashtab[savetab[savenum].u.movedata.hashnum].vertnum & velomsk];

            hashtab[savetab[savenum].u.movedata.hashnum].partval = savetab[savenum].u.movedata.partval;
            break;

          case WGRAPHSEPAFMSAVELINKDEL :
            linkptr2 = linktab + savetab[savenum].u.linkdata.linknum;
            linkptr2->gain = savetab[savenum].u.linkdata.gain;
            gainTablAdd (tablptr, (GainLink *) linkptr2, linkptr2->gain); /* Add link into part gain table */
            break;

          case WGRAPHSEPAFMSAVELINKADD :
            linkptr2 = linktab + savetab[savenum].u.linkdata.linknum;
            if (linkptr2->gainlink.next != (GainLink *) 1) {
              gainTablDel (tablptr, (GainLink *) linkptr2); /* Remove link from table */
              linkptr2->gainlink.next = NULL;
            }
            break;

          case WGRAPHSEPAFMSAVELOAD:
            grafptr->compload[savetab[savenum].u.loaddata.partval] -= savetab[savenum].u.loaddata.loaddiff;
            grafptr->compsize[savetab[savenum].u.loaddata.partval] -= savetab[savenum].u.loaddata.sizediff;
            break;
        }
      }
    }
  } while ((moveflag != 0) &&                     /* As long as vertices are moved */
           (-- passnbr > 0));                     /* and we are allowed to loop    */

  grafptr->fronload = 0;

  for (vexxptr = hashtab, fronnum = 0;            /* Build new frontier                */
       vexxptr < hashtab + (hashmax << 2); vexxptr ++) { /* from all vertices in table */
    Gnum                vertnum;

    vertnum = vexxptr->vertnum;
    if (vertnum != ~0) {                          /* If vertex slot is used     */
      Gnum                 partval;               /* New part of current vertex */
      Gnum                 partold;               /* Old part of current vertex */

      partval = vexxptr->partval;
      partold = grafptr->parttax[vertnum];        /* Get old part value from array */
      if (partval != partold)                     /* If vertex part changed        */
        grafptr->parttax[vertnum] = partval;      /* Set new part value            */

      if (partval == -1)
      {
        grafptr->fronload += velobax[vertnum & velomsk];
        grafptr->frontab[fronnum ++] = vexxptr->vertnum;
      }
    }
  }

  grafptr->fronnbr = fronnum;
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
wgraphPartFmResize ()

{
  errorPrint ("wgraphPartFmResize: not implemented");

  return (0);
}
