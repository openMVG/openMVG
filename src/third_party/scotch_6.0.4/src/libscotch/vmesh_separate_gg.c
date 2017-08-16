/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vmesh_separate_gg.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates a node separation **/
/**                mesh using an element-oriented version  **/
/**                of the Greedy Graph Growing algorithm.  **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 16 sep 2002     **/
/**                                 to     18 aug 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     24 mar 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     09 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_SEPARATE_GG

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "mesh.h"
#include "vmesh.h"
#include "vmesh_separate_gg.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
vmeshSeparateGg (
Vmesh * restrict const                      meshptr, /*+ Node separation mesh +*/
const VmeshSeparateGgParam * restrict const paraptr) /*+ Method parameters    +*/
{
  GainTabl * restrict               tablptr;      /* Pointer to gain table                                              */
  byte * restrict                   vexxtab;      /* Start of auxiliary arrays                                          */
  Gnum                              vexxsiz;      /* Size of auxiliary arrays to be reset every pass                    */
  VmeshSeparateGgElem * restrict    velxtax;      /* Based auxiliary element array                                      */
  VmeshSeparateGgNode * restrict    vnoxtax;      /* Based auxiliary node array                                         */
  Gnum * restrict                   velitax;      /* Array of sums of weights of isolated neighboring nodes of elements */
  Gnum * restrict                   velstax;      /* Array of sums of weights of neighboring nodes of elements          */
  Gnum                              velssiz;      /* Size of element neighboring node load sum array                    */
  VmeshSeparateGgElem *             sepaptr;      /* Head of chained list of elements to re-link                        */
  Gnum * restrict                   permtab;      /* Element permutation table for finding new roots                    */
  Gnum *                            permptr;      /* Pointer to current permutation index                               */
  INT                               passnum;      /* Number of current pass                                             */
  Gnum                              ecmpsize0;    /* Number of elements in part 0         */
  Gnum                              ncmploaddlt;  /* Current imbalance of bipartition     */
  Gnum                              ncmpload2;    /* Current number of nodes in separator */
  Gnum                              vnodnum;
  Gnum                              fronnum;
  Gnum                              ncmpsize1;
  Gnum                              ncmpsize2;

  if (meshptr->m.velmnbr == 0) {                  /* If only a single node or disconnected nodes */
    vmeshZero (meshptr);                          /* Don't bother with parts                     */
    return    (0);
  }

  velssiz = (meshptr->m.vnlotax == NULL) ? 0 : meshptr->m.velmnbr; /*  Compute size of vetex load sum array */
  if (((tablptr = gainTablInit (GAINMAX, VMESHSEPAGGSUBBITS)) == NULL) || /* Use logarithmic array only     */
      ((vexxtab = (byte *) memAllocGroup ((void **) (void *)
                             &velxtax, (size_t) (meshptr->m.velmnbr * sizeof (VmeshSeparateGgElem)),
                             &vnoxtax, (size_t) (meshptr->m.vnodnbr * sizeof (VmeshSeparateGgNode)),
                             &velitax, (size_t) (meshptr->m.velmnbr * sizeof (Gnum)),
                             &velstax, (size_t) (velssiz            * sizeof (Gnum)), NULL)) == NULL)) { /* Indicates end of group allocated array */
    if (tablptr != NULL)
      gainTablExit (tablptr);
    errorPrint ("vmeshSeparateGg: out of memory (1)");
    return     (1);
  }
  vexxsiz  = (byte *) velitax - vexxtab;          /* Size of arrays that must be reset at each pass */
  velxtax -= meshptr->m.velmbas;                  /* Base access to auxiliary arrays                */
  vnoxtax -= meshptr->m.vnodbas;
  velitax -= meshptr->m.velmbas;

  if (velssiz == 0) {                             /* If no vertex load array */
    Gnum                velmnum;

    for (velmnum = meshptr->m.velmbas; velmnum < meshptr->m.velmnnd; velmnum ++) {
      Gnum                eelmnum;
      Gnum                velisum;

      for (eelmnum = meshptr->m.verttax[velmnum], velisum = 0;
           eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
        Gnum                vnodnum;

        vnodnum = meshptr->m.edgetax[eelmnum];
        if ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]) == 1)
          velisum --;
      }
      velitax[velmnum] = velisum;
    }
  }
  else {
    Gnum                velmnum;

    velstax -= meshptr->m.velmbas;

    for (velmnum = meshptr->m.velmbas; velmnum < meshptr->m.velmnnd; velmnum ++) {
      Gnum                eelmnum;
      Gnum                velisum;
      Gnum                velssum;

      for (eelmnum = meshptr->m.verttax[velmnum], velisum = velssum = 0;
           eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
        Gnum                vnodnum;
        Gnum                vnloval;

        vnodnum  = meshptr->m.edgetax[eelmnum];
        vnloval  = meshptr->m.vnlotax[vnodnum];
        velssum += vnloval;
        if ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]) == 1)
          velisum -= vnloval;
      }
      velitax[velmnum] = velisum;
      velstax[velmnum] = velssum;
    }
  }

  permtab = NULL;                                 /* Do not allocate permutation array yet */
  for (passnum = 0; passnum < paraptr->passnbr; passnum ++) { /* For all passes            */
    VmeshSeparateGgElem *             velxptr;    /* Pointer to selected element           */

    memSet (vexxtab, 0, vexxsiz);                 /* All vertices to part 0      */
    gainTablFree (tablptr);                       /* Reset gain table            */
    permptr     = NULL;                           /* No permutation built yet    */
    ecmpsize0   = meshptr->m.velmnbr;             /* All elements to part 0      */
    ncmpload2   = 0;                              /* Reset separation parameters */
    ncmploaddlt = meshptr->m.vnlosum;

    velxptr = (VmeshSeparateGgElem *) vexxtab + intRandVal (meshptr->m.velmnbr); /* Randomly select first root element vertex */

    do {                                          /* Loop on root element vertices        */
      Gnum                velmnum;                /* Number of current element to process */

      velxptr->gainlink.next =                    /* TRICK: allow deletion of root vertex */
      velxptr->gainlink.prev = (GainLink *) velxptr;

      velmnum = velxptr - velxtax;                /* Get root element number */
      {
        Gnum                ncmpgain1;            /* Gain (2->1) */
        Gnum                ncmpgain2;            /* Gain (0->2) */

        ncmpgain2 = (meshptr->m.vnlotax == NULL)  /* Set gains */
                    ? meshptr->m.vendtax[velmnum] - meshptr->m.verttax[velmnum]
                    : velstax[velmnum];
        ncmpgain1 = velitax[velmnum];

        velxptr->ncmpgain2   = ncmpgain1 + ncmpgain2;
        velxptr->ncmpgaindlt = ncmpgain1 - ncmpgain2;
      }

      do {                                        /* While element vertices can be retrieved */
        Gnum                eelmnum;              /* Number of current element edge          */

        velmnum = velxptr - velxtax;              /* Get based number of selected element */

        if (ncmploaddlt < abs (ncmploaddlt + velxtax[velmnum].ncmpgaindlt)) { /* If swapping would cause imbalance */
          permptr = permtab + meshptr->m.velmnbr; /* Terminate swapping process                                    */
          velxptr = NULL;
          break;
        }
        ecmpsize0 --;                             /* One less element in part 0   */
        gainTablDel (tablptr, (GainLink *) velxptr); /* Remove element from table */
        velxptr->gainlink.next = VMESHSEPAGGSTATEPART1; /* Move element to part 1 */
        ncmpload2   += velxptr->ncmpgain2;        /* Update partition parameters  */
        ncmploaddlt += velxptr->ncmpgaindlt;           

        sepaptr = NULL;                           /* No frontier elements to relink yet */
        for (eelmnum = meshptr->m.verttax[velmnum]; /* For all neighbor node vertices   */
             eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
          Gnum                vnodnum;            /* Number of current node neighbor */

          vnodnum = meshptr->m.edgetax[eelmnum];  /* Get number of neighbor node */

#ifdef SCOTCH_DEBUG_VMESH2
          if (vnoxtax[vnodnum].partval == 1) {
            errorPrint ("vmeshSeparateGg: internal error (1)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_VMESH2 */
          if (vnoxtax[vnodnum].partval == 0) {    /* If yet untouched neighbor node        */
            Gnum                enodnum;          /* Current egde of current neighbor node */
            Gnum                vnloval;

            vnloval = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];
            vnoxtax[vnodnum].partval   = 2;       /* Node now belongs to separator */
            vnoxtax[vnodnum].commsize0 = meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum] - 1;
            vnoxtax[vnodnum].velmisum0 = - velmnum;

            for (enodnum = meshptr->m.verttax[vnodnum]; /* For all its elements */
                 enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
              Gnum                velmend;

              velmend = meshptr->m.edgetax[enodnum]; /* Get neighbor element */

              vnoxtax[vnodnum].velmisum0 += velmend; /* Sum-up element indices for neighbor node */
#ifdef SCOTCH_DEBUG_VMESH2
              if ((velxtax[velmend].gainlink.next == VMESHSEPAGGSTATEPART1) &&
                  (velmend != velmnum)) {
                errorPrint ("vmeshSeparateGg: internal error (2)");
                return     (1);
              }
#endif /* SCOTCH_DEBUG_VMESH2 */
              if (velxtax[velmend].gainlink.next == VMESHSEPAGGSTATEPART0) { /* If untouched element */
                Gnum                ncmpgain1;    /* Gain (2->1) */
                Gnum                ncmpgain2;    /* Gain (0->2) */
#ifdef SCOTCH_DEBUG_VMESH2
                Gnum                eelmend;
#endif /* SCOTCH_DEBUG_VMESH2 */

                velxtax[velmend].gainlink.next = VMESHSEPAGGSTATEPART2; /* Move element to frontier */
                velxtax[velmend].gainlink.prev = (GainLink *) sepaptr; /* Chain vertex              */
                sepaptr                        = &velxtax[velmend];

                ncmpgain2 = (meshptr->m.vnlotax == NULL) /* Set gains */
                            ? meshptr->m.vendtax[velmend] - meshptr->m.verttax[velmend] - 1
                            : velstax[velmend] - vnloval;
                ncmpgain1 = velitax[velmend];

#ifdef SCOTCH_DEBUG_VMESH2
                for (eelmend = meshptr->m.verttax[velmend]; /* For all its neighboring nodes */
                     eelmend < meshptr->m.vendtax[velmend]; eelmend ++) {
                  Gnum                vnodend;

                  vnodend = meshptr->m.edgetax[eelmend];

                  if ((vnoxtax[vnodend].partval == 1) ||
                      ((vnoxtax[vnodend].partval == 2) && (vnodend != vnodnum))) {
                    errorPrint ("vmeshSeparateGg: internal error (3)");
                    return     (1);
                  }
                  if (meshptr->m.vendtax[vnodend] - meshptr->m.verttax[vnodend] == 1) {
                    if (vnoxtax[vnodend].partval != 0) {
                      errorPrint ("vmeshSeparateGg: internal error (4)");
                      return     (1);
                    }
                  }
                }
#endif /* SCOTCH_DEBUG_VMESH2 */
                velxtax[velmend].ncmpgain2   = ncmpgain1 + ncmpgain2;
                velxtax[velmend].ncmpgaindlt = ncmpgain1 - ncmpgain2;
              }
              else {                              /* Neighbor element belongs to frontier                 */
                velxtax[velmend].ncmpgain2   -= vnloval; /* One less node to add to separator for element */
                velxtax[velmend].ncmpgaindlt += vnloval; /* One less node to remove from part 0           */

                if (velxtax[velmend].gainlink.next >= VMESHSEPAGGSTATELINK) {
                  gainTablDel (tablptr, (GainLink *) &velxtax[velmend]); /* Unlink vertex */
                  velxtax[velmend].gainlink.next = VMESHSEPAGGSTATEPART2; /* Chain vertex */
                  velxtax[velmend].gainlink.prev = (GainLink *) sepaptr;
                  sepaptr                        = &velxtax[velmend];
                }
              }
            }
          }
          else {                                  /* Neighbor node already in separator  */
            vnoxtax[vnodnum].commsize0 --;        /* One less neighbor element in part 0 */
            vnoxtax[vnodnum].velmisum0 -= velmnum; /* Subtract index of removed element  */
          }

          if (vnoxtax[vnodnum].commsize0 == 0)    /* If node no longer has neighbors in part 0  */
            vnoxtax[vnodnum].partval = 1;         /* Node moves from separator to part 1        */
          else if (vnoxtax[vnodnum].commsize0 == 1) { /* If only one neighbor element in part 0 */
            Gnum                velmend;          /* Index of remaining element in part 0       */
            Gnum                vnloval;

            velmend = vnoxtax[vnodnum].velmisum0; /* Get neighbor element from remaining index */

#ifdef SCOTCH_DEBUG_VMESH2
            if (velxtax[velmend].gainlink.next < VMESHSEPAGGSTATEPART2) { /* Element should have been declared in part 0 at this stage */
              errorPrint ("vmeshSeparateGg: internal error (5)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VMESH2 */

            vnloval = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];
            velxtax[velmend].ncmpgain2   -= vnloval;
            velxtax[velmend].ncmpgaindlt -= vnloval;

            if (velxtax[velmend].gainlink.next >= VMESHSEPAGGSTATELINK) {
              gainTablDel (tablptr, (GainLink *) &velxtax[velmend]); /* Unlink vertex */
              velxtax[velmend].gainlink.next = VMESHSEPAGGSTATEPART2; /* Chain vertex */
              velxtax[velmend].gainlink.prev = (GainLink *) sepaptr;
              sepaptr                        = &velxtax[velmend];
            }
          }
        }

        while (sepaptr != NULL) {                 /* For all vertices in chain list */
          velxptr = sepaptr;                      /* Unlink vertex from list        */
          sepaptr = (VmeshSeparateGgElem *) velxptr->gainlink.prev;
          gainTablAdd (tablptr, (GainLink *) velxptr, velxptr->ncmpgain2); /* Relink it */
        }

#ifdef SCOTCH_DEBUG_VMESH3
        if (vmeshSeparateGgCheck (meshptr, ncmpload2, ncmploaddlt, velxtax, vnoxtax) != 0) {
          errorPrint ("vmeshSeparateGg: internal error (6)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VMESH3 */
      } while ((velxptr = (VmeshSeparateGgElem *) gainTablFrst (tablptr)) != NULL);

      if (permptr == NULL) {                      /* If element permutation not yet built   */
        if (permtab == NULL) {                    /* If permutation array not yet allocated */
          if ((permtab = (Gnum *) memAlloc (meshptr->m.velmnbr * sizeof (Gnum))) == NULL) {
            errorPrint   ("vmeshSeparateGg: out of memory (2)");
            memFree      (vexxtab);
            gainTablExit (tablptr);
            return (1);
          }
          intAscn (permtab, meshptr->m.velmnbr, meshptr->m.baseval); /* Initialize permutation array */
        }
        intPerm (permtab, meshptr->m.velmnbr);    /* Build random permutation          */
        permptr = permtab;                        /* Start at beginning of permutation */
      }
      for ( ; permptr < permtab + meshptr->m.velmnbr; permptr ++) { /* Find next root vertex */
#ifdef SCOTCH_DEBUG_VMESH2
        if (velxtax[*permptr].gainlink.next >= VMESHSEPAGGSTATEPART2) {
          errorPrint ("vmeshSeparateGg: internal error (7)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VMESH2 */
        if (velxtax[*permptr].gainlink.next == VMESHSEPAGGSTATEPART0) {
          velxptr = velxtax + (*permptr ++);
          break;
        }
      }
    } while (velxptr != NULL);

    if ((passnum == 0) ||                         /* If it is the first try         */
        ( (meshptr->ncmpload[2] >  ncmpload2) ||  /* Or if better solution reached */
         ((meshptr->ncmpload[2] == ncmpload2) &&
          (abs (meshptr->ncmploaddlt) > abs (ncmploaddlt))))) {
      Gnum                vertnum;

      meshptr->ecmpsize[0] = ecmpsize0;           /* Set graph parameters */
      meshptr->ncmpload[2] = ncmpload2;
      meshptr->ncmploaddlt = ncmploaddlt;

      for (vertnum = meshptr->m.velmbas; vertnum < meshptr->m.velmnnd; vertnum ++) /* Copy element bipartition state */
        meshptr->parttax[vertnum] = (velxtax[vertnum].gainlink.next == VMESHSEPAGGSTATEPART1) ? 1 : 0;
      for (vertnum = meshptr->m.vnodbas; vertnum < meshptr->m.vnodnnd; vertnum ++) /* Copy node bipartition state */
        meshptr->parttax[vertnum] = vnoxtax[vertnum].partval;
    }
  }

  meshptr->ecmpsize[1] = meshptr->m.velmnbr - meshptr->ecmpsize[0];
  meshptr->ncmpload[1] = ((meshptr->m.vnlosum - meshptr->ncmpload[2]) - meshptr->ncmploaddlt) >> 1;
  meshptr->ncmpload[0] = (meshptr->m.vnlosum - meshptr->ncmpload[2]) - meshptr->ncmpload[1];

  for (vnodnum = meshptr->m.vnodbas, fronnum = 0, ncmpsize1 = ncmpsize2 = 0;
       vnodnum < meshptr->m.vnodnnd; vnodnum ++) {
    Gnum                partval;

    partval = meshptr->parttax[vnodnum];
    ncmpsize1 += (partval &  1);                  /* Superscalar update */
    ncmpsize2 += (partval >> 1);
    if (partval == 2)
      meshptr->frontab[fronnum ++] = vnodnum;     /* Vertex belongs to frontier */
  }
  meshptr->ncmpsize[0] = meshptr->m.vnodnbr - (ncmpsize1 + ncmpsize2);
  meshptr->ncmpsize[1] = ncmpsize1;
  meshptr->fronnbr     = ncmpsize2;

#ifdef SCOTCH_DEBUG_VMESH2
  if (vmeshCheck (meshptr) != 0) {
    errorPrint ("vmeshSeparateGg: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  if (permtab != NULL)
    memFree (permtab);
  memFree      (vexxtab);                         /* Free group leader */
  gainTablExit (tablptr);

/* printf ("GG Sepa\tsize=%ld\tload=%ld\tbal=%ld\n", (long) meshptr->fronnbr, (long) meshptr->ncmpload[2], (long) meshptr->ncmploaddlt); */

  return (0);
}

/* This routine checks the consistency
** of the current bipartition.
** It returns:
** - 0   : if the bipartition is consistent.
** - !0  : on error.
*/

#ifdef SCOTCH_DEBUG_VMESH3
static
int
vmeshSeparateGgCheck (
Vmesh * restrict const                      meshptr,
const Gnum                                  ncmpload2,
const Gnum                                  ncmploaddlt,
const VmeshSeparateGgElem * restrict const  velxtax,
const VmeshSeparateGgNode * restrict const  vnoxtax)
{
  Gnum                vnodnum;
  Gnum                vnloval;
  Gnum                ncmpsize0c;
  Gnum                ncmpsize1c;
  Gnum                ncmpsize2c;
  Gnum                ncmpload0c;
  Gnum                ncmpload1c;
  Gnum                ncmpload2c;

  ncmpsize1c =
  ncmpsize2c =
  ncmpload1c =
  ncmpload2c = 0;
  vnloval    = 1;
  for (vnodnum = meshptr->m.vnodbas; vnodnum < meshptr->m.vnodnnd; vnodnum ++) {
    int                 partval;
    Gnum                partval1;
    Gnum                partval2;

    partval  = vnoxtax[vnodnum].partval;
    partval1 = partval &  1;
    partval2 = partval >> 1;
    if (meshptr->m.vnlotax != NULL)
      vnloval = meshptr->m.vnlotax[vnodnum];
    if (partval > 2) {
      errorPrint ("vmeshSeparateGgCheck: invalid node part value");
      return     (1);
    }
    ncmpsize1c += partval1;
    ncmpsize2c += partval2;
    ncmpload1c += partval1 * vnloval;
    ncmpload2c += partval2 * vnloval;
  }
  ncmpsize0c = meshptr->m.vnodnbr - ncmpsize1c - ncmpsize2c;
  ncmpload0c = meshptr->m.vnlosum - ncmpload1c - ncmpload2c;

  if (ncmpload2c != ncmpload2) {
    errorPrint ("vmeshSeparateGgCheck: invalid separator size");
    return     (1);
  }
  if (ncmploaddlt != (ncmpload0c - ncmpload1c)) {
    errorPrint ("vmeshSeparateGgCheck: invalid separator balance");
    return     (1);
  }

  return (0);
}
#endif /* SCOTCH_DEBUG_VMESH3 */
