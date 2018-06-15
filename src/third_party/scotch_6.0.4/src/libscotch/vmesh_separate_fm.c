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
/**   NAME       : vmesh_separate_fm.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates an active         **/
/**                mesh using an element-oriented version  **/
/**                of our improved Fiduccia-Mattheyses     **/
/**                heuristics.                             **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 26 feb 2003     **/
/**                                 to     06 may 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     22 may 2008     **/
/**                # Version 5.1  : from : 12 nov 2008     **/
/**                                 to     12 nov 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_SEPARATE_FM

/* #define SCOTCH_DEBUG_VMESH3 */              /* For intensive debugging */

#include "module.h"
#include "common.h"
#include "gain.h"
#include "graph.h"
#include "mesh.h"
#include "vmesh.h"
#include "vmesh_separate_gg.h"
#include "vmesh_separate_fm.h"

/* This routine resizes the hash arrays
** as well as the associated structures.
** In the group of allocated arrays,
** the element array must be put before
** the node array, because the element
** structure is larger than the node
** structure, such that the old and new
** node arrays can never overlap after
** the doubling in size of the element
** array. The same for the move array.
** It returns:
** - 0   : if resize succeeded.
** - !0  : in case of error.
*/

static
int
vmeshSeparateFmResize (
GainTabl * restrict const                 tablptr, /*+ Pointer to gain table                         +*/
VmeshSeparateFmElement * restrict * const helmptr, /*+ Pointer to pointer to element hash table      +*/
VmeshSeparateFmNode * restrict * const    hnodptr, /*+ Pointer to pointer to node hash table         +*/
VmeshSeparateFmSave * restrict * const    saveptr, /*+ Pointer to pointer to move array              +*/
const Gnum                                savenbr, /*+ Current number of items in save array         +*/
VmeshSeparateFmElement **                 lockptr, /*+ Pointer to list of locked elements            +*/
VmeshSeparateFmElement **                 sepaptr, /*+ Pointer to list of separator elements, if any +*/
const Gnum                                hashold) /*+ Maximum number of vertices in hash structures +*/
{
  Gnum                              hashsiz;      /* Size of hash table                   */
  Gnum                              hashmsk;      /* Mask for access to hash table        */
  Gnum                              hashmax;      /* Maximum number of objects in tables  */
  VmeshSeparateFmSave * restrict    movetab;      /* Pointer to move array                */
  VmeshSeparateFmElement * restrict helmtab;      /* Element hash table                   */
  VmeshSeparateFmNode *             hnodtab;      /* Node hash table                      */
  size_t                            addradj;      /* Address adjustment                   */
  Gnum                              helmold;
  VmeshSeparateFmNode *             hnodtld;
  Gnum                              hnodold;
  VmeshSeparateFmSave * restrict    savetab;
  Gnum                              savenum;

  hashmax = 2 * hashold;                          /* Set new number */
  hashsiz = 4 * hashmax;                          /* Set new size   */
  hashmsk = hashsiz - 1;

  savetab = *saveptr;                             /* Point to old move array                  */
  for (savenum = 0; savenum < savenbr; savenum ++) { /* Turn hash indices into vertex indices */
    Gnum                hertnum;

    hertnum = savetab[savenum].hertnum;
    savetab[savenum].hertnum = (hertnum >= 0) ? (*helmptr)[hertnum].velmnum : (-1 - (*hnodptr)[-1 - hertnum].vnodnum);
  }

  if (memReallocGroup ((void *) *helmptr,         /* Get old group leader */
        &helmtab, (size_t) (hashsiz * sizeof (VmeshSeparateFmElement)),
        &hnodtab, (size_t) (hashsiz * sizeof (VmeshSeparateFmNode)),
        &movetab, (size_t) (hashmax * sizeof (VmeshSeparateFmSave)), NULL) == NULL) {
    errorPrint ("vmeshSeparateFmResize: cannot resize arrays");
    return     (1);                               /* If cannot reallocate */
  }
#ifdef SCOTCH_DEBUG_VMESH2
  if (((byte *) hnodtab - (byte *) helmtab) < ((byte *) (*saveptr) - (byte *) (*helmptr))) { /* If cannot simply copy node hash array */
    errorPrint ("vmeshSeparateFmResize: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  memMov (movetab, ((byte *) helmtab) + ((byte *) *saveptr - (byte *) *helmptr), savenbr * sizeof (VmeshSeparateFmSave)); /* Old array may have moved but arrays cannot overlap */

  memSet (hnodtab, ~0, hashsiz * sizeof (VmeshSeparateFmNode)); /* Cannot overlap */
  hnodtld = (VmeshSeparateFmNode *) ((byte *) helmtab) + ((byte *) *hnodptr - (byte *) *helmptr); /* Point to old node array */
  for (hnodold = 0; hnodold < (hashold * 4); hnodold ++) { /* For all old allocated nodes */
    Gnum                              hnodnew;

    if (hnodtld[hnodold].vnodnum == ~0)           /* If unallocated slot */
      continue;                                   /* Skip to next slot   */

    for (hnodnew = (hnodtld[hnodold].vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk;
         hnodtab[hnodnew].vnodnum != ~0; hnodnew = (hnodnew + 1) & hashmsk) ;
    hnodtab[hnodnew] = hnodtld[hnodold];          /* Move node data to new position */
  }

/* TODO */
  fprintf (stderr, "hertnum no longer valid !\n");
  exit    (1);

  addradj = (byte *) helmtab - (byte *) (*helmptr); /* Compute address difference */

  gainTablFree (tablptr);                         /* Reset gain table */
  memSet (helmtab + hashold, ~0, hashold * 2 * sizeof (VmeshSeparateFmElement));
  for (helmold = 0; helmold < (hashold * 4); helmold ++) { /* For all old allocated elements */
    Gnum                              helmnew;

    if (helmtab[helmold].velmnum == ~0)           /* If unallocated slot */
      continue;                                   /* Skip to next slot   */

    for (helmnew = (helmtab[helmold].velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; helmnew = (helmnew + 1) & hashmsk) {
      if (helmtab[helmnew].velmnum == ~0) {
        helmtab[helmnew].velmnum     = helmtab[helmold].velmnum;
        helmtab[helmnew].vertpart    = helmtab[helmold].vertpart;
        helmtab[helmnew].ncmpcut2    = helmtab[helmold].ncmpcut2;
        helmtab[helmnew].ncmpgain2   = helmtab[helmold].ncmpgain2;
        helmtab[helmnew].ncmpgaindlt = helmtab[helmold].ncmpgaindlt;
        helmtab[helmnew].mswpnum     = helmtab[helmold].mswpnum;
        helmtab[helmold].velmnum     = ~0;        /* Free old slot      */
        helmtab[helmold].mswpnum     = ~0;        /* Reset sweep number */
        break;
      }
      if (helmtab[helmnew].velmnum != helmtab[helmold].velmnum) /* If element not found */
        continue;                                 /* Go on searching                    */
    }
    if (helmtab[helmold].gainlink.next >= VMESHSEPAFMSTATELINK) /* If element was linked               */
      gainTablAdd (tablptr, (GainLink *) &helmtab[helmnew], helmtab[helmnew].ncmpgain2); /* Re-link it */
    else {                                        /* Element may be chained in some list               */
      helmtab[helmnew].gainlink.next = helmtab[helmold].gainlink.next; /* Save it                      */
      helmtab[helmnew].gainlink.prev = (GainLink *) ((byte *) helmtab[helmold].gainlink.prev + addradj);
    }
  }

  if (*lockptr != NULL)
    *lockptr = (VmeshSeparateFmElement *) ((byte *) (*lockptr) + addradj);
  if (sepaptr != NULL) {
    if (*sepaptr != NULL)
      *sepaptr = (VmeshSeparateFmElement *) ((byte *) (*sepaptr) + addradj);
  }

  for (savenum = 0; savenum < savenbr; savenum ++) { /* Turn vertex indices back into hash indices */
    Gnum                vertnum;

    vertnum = movetab[savenum].hertnum;           /* Read vertex index */
    if (vertnum >= 0) {                           /* If element vertex */
      Gnum                helmnum;

      for (helmnum = (vertnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; helmnum = (helmnum + 1) & hashmsk) {
#ifdef SCOTCH_DEBUG_VMESH2
        if (helmtab[helmnum].velmnum == ~0) {     /* We should always find the elements */
          errorPrint ("vmeshSeparateFmResize: internal error (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VMESH2 */
        if (helmtab[helmnum].velmnum == vertnum)  /* If element found */
          break;
      }
      movetab[savenum].hertnum = helmnum;         /* Save element hash index */
    }
    else {                                        /* If node vertex */
      Gnum                hnodnum;

      vertnum = -1 - vertnum;
      for (hnodnum = (vertnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; hnodnum = (hnodnum + 1) & hashmsk) {
#ifdef SCOTCH_DEBUG_VMESH2
        if (hnodtab[hnodnum].vnodnum == ~0) {     /* We should always find the nodes */
          errorPrint ("vmeshSeparateFmResize: internal error (3)");
          return     (1);
        }
        if (hnodtab[hnodnum].vnodnum == vertnum)  /* If element found */
          break;
#endif /* SCOTCH_DEBUG_VMESH2 */
      }
      movetab[savenum].hertnum = -1 - hnodnum;    /* Save node hash index */
    }
  }

  fprintf (stderr, "########### vmeshSeparateFmResize (%ld) !!!\n", (long) hashold);

  return (0);
}

/* This routine returns the vertex of best gain
** whose swap will keep the balance correct.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no more vertices available.
*/

static
VmeshSeparateFmElement *
vmeshSeparateFmTablGet (
GainTabl * const            tablptr,              /*+ Gain table        +*/
const Gnum                  deltcur,              /*+ Current imbalance +*/
const Gnum                  deltmax)              /*+ Maximum imbalance +*/
{
  const VmeshSeparateFmElement *  velmptr;
  VmeshSeparateFmElement *        vertbest;
  Gnum                            gainbest;
  const GainEntr *                tablbest;
  Gnum                            deltbest;
  Gnum                            deltnew;

  tablbest = tablptr->tend;                       /* Assume no candidate vertex found yet */
  gainbest = GAINMAX;
  vertbest = NULL;
  deltbest = deltmax;

  for (velmptr = (VmeshSeparateFmElement *) gainTablFrst (tablptr); /* Select candidate vertices */
       (velmptr != NULL) && (velmptr->gainlink.tabl < tablbest);
       velmptr = (VmeshSeparateFmElement *) gainTablNext (tablptr, &velmptr->gainlink)) {
    deltnew = abs (deltcur + velmptr->ncmpgaindlt);
    if (deltnew <= deltmax) {                     /* If vertex enforces balance  */
      if ((velmptr->ncmpgain2 < gainbest) ||      /* And if it gives better gain */
          ((velmptr->ncmpgain2 == gainbest) &&    /* Or if it gives better load  */
           (deltnew < deltbest))) {
        tablbest = velmptr->gainlink.tabl;        /* Select it */
        gainbest = velmptr->ncmpgain2;
        vertbest = (VmeshSeparateFmElement *) velmptr;
        deltbest = deltnew;
      }
    }
  }

  return (vertbest);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes the
** separation of the given mesh.
** It returns:
** - 0   : if bipartitioning could be computed.
** - !0  : on error.
*/

int
vmeshSeparateFm (
Vmesh * restrict const                      meshptr, /*+ Node separation mesh +*/
const VmeshSeparateFmParam * restrict const paraptr) /*+ Method parameters    +*/
{
  GainTabl * restrict               tablptr;      /* Pointer to gain table                   */
  long                              passnbr;      /* Maximum number of passes to go          */
  VmeshSeparateFmSave * restrict    movetab;      /* Pointer to move array                   */
  Gnum                              movenbr;      /* Number of uneffective moves done        */
  Gnum                              savenbr;      /* Number of recorded backtrack moves      */
  Gnum                              mswpnum;      /* Current number of recording sweep       */
  int                               moveflag;     /* Flag set if useful moves made           */
  Gnum                              fronnum;      /* Current index in frontier array         */
  Gnum                              vertnbr;
  Gnum                              hashsiz;      /* Size of hash table                      */
  Gnum                              hashmsk;      /* Mask for access to hash table           */
  Gnum                              hashmax;      /* Maximum number of objects in tables     */
  Gnum                              hashnbr;      /* Estimated number of onjects in hash     */
  Gnum                              helmnbr;      /* Number of elements in hash table        */
  Gnum                              helmnum;
  Gnum                              hnodnum;
  Gnum                              hnodnbr;      /* Number of nodes in hash table           */
  VmeshSeparateFmElement * restrict helmtab;      /* Element hash table                      */
  VmeshSeparateFmNode * restrict    hnodtab;      /* Node hash table                         */
  VmeshSeparateFmElement *          lockptr;      /* Linked list of locked elements          */
  VmeshSeparateFmElement *          velmptr;      /* Pointer to current element              */
  Gnum                              ncmploaddltmat; /* Theoretical latgest imbalance allowed */
  Gnum                              ncmploaddltmax; /* Largest imbalance allowed             */
  Gnum                              ncmploaddlt;  /* Current imbalance                       */
  Gnum                              ncmpload2;    /* Current size of separator               */
  Gnum                              ncmpload2bst;
  Gnum                              ncmploaddltbst;
  Gnum                              ecmpload1;
  Gnum                              ncmpload1;
  Gnum                              ncmpsize1;
  Gnum                              ncmpsize2;

  if (paraptr->deltrat > 0.0L) {
    Gnum                ncmploaddlttmp;

    ncmploaddltmat = (Gnum) (paraptr->deltrat * meshptr->m.vnlosum) + 1;
    ncmploaddlttmp = (Gnum) (((float) meshptr->m.edgenbr * (float) meshptr->m.vnlosum) /
                             ((float) meshptr->m.velmnbr   * (float) meshptr->m.vnodnbr));
    if (ncmploaddltmat < ncmploaddlttmp)
      ncmploaddltmat = ncmploaddlttmp;
  }
  else
    ncmploaddltmat = 0;

  ncmploaddltmat = (paraptr->deltrat > 0.0L) ? ((Gnum) (paraptr->deltrat * meshptr->m.vnlosum) + 1) : 0;

/* printf ("FM Mbal=%ld\n", (long) ncmploaddltmat); */

  if ((meshptr->fronnbr == 0) &&                  /* If imbalance in graph with no frontier */
      (abs (meshptr->ncmploaddlt) > ncmploaddltmat)) {
    VmeshSeparateGgParam  paramdat;

    paramdat.passnbr = 3;
    vmeshSeparateGg (meshptr, &paramdat);         /* Compute a balanced initial partition */
  }

  vertnbr = meshptr->m.velmnbr + meshptr->m.vnodnbr;
  hashnbr = 2 * ((meshptr->fronnbr + paraptr->movenbr) * (1 + (Gnum) ((float) meshptr->m.edgenbr / (float) vertnbr)));
  if (hashnbr > vertnbr)                          /* Set bound on hash table */
    hashnbr = vertnbr;
  for (hashmax = 256; hashmax < hashnbr; hashmax <<= 1) ; /* Get upper power of two */
/* TODO */ hashmax *= 4;
  hashsiz = 4 * hashmax;
  hashmsk = hashsiz - 1;

  if (((tablptr = gainTablInit (meshptr->m.vnlosum, VMESHSEPAFMGAINBITS)) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &helmtab, (size_t) (hashsiz * sizeof (VmeshSeparateFmElement)),
                      &hnodtab, (size_t) (hashsiz * sizeof (VmeshSeparateFmNode)),
                      &movetab, (size_t) (hashmax * sizeof (VmeshSeparateFmSave)), NULL) == NULL)) {
    if (tablptr != NULL) {
      errorPrint   ("vmeshSeparateFm: out of memory (1)");
      gainTablExit (tablptr);
    }
    return (1);
  }

  passnbr        = paraptr->passnbr;              /* Set remaining number of passes */
  ncmpload2      = meshptr->ncmpload[2];          /* Set current partition loads    */
  ncmploaddlt    = meshptr->ncmploaddlt;

  memSet (helmtab, ~0, (byte *) &hnodtab[hashsiz] - (byte *) helmtab); /* Set all vertex numbers to ~0 */

  helmnbr =
  hnodnbr = 0;
  savenbr = 0;                                    /* No recorded moves yet         */
  lockptr = NULL;                                 /* Set locked list as empty      */
  for (fronnum = 0; fronnum < meshptr->fronnbr; fronnum ++) { /* Set initial gains */
    Gnum                vnloval;
    Gnum                vnodnum;
    Gnum                enodnum;
    Gnum                hnodnum;
    Gnum                ecmpsize1;

    vnodnum = meshptr->frontab[fronnum];
#ifdef SCOTCH_DEBUG_VMESH2
    if (meshptr->parttax[vnodnum] != 2) {
      errorPrint ("vmeshSeparateFm: internal error (1)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VMESH2 */
    vnloval = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];

    for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; hnodnum = (hnodnum + 1) & hashmsk) {
      if (hnodtab[hnodnum].vnodnum == ~0)         /* If node slot found */
        break;                                    /* No need to go on   */
#ifdef SCOTCH_DEBUG_VMESH2
      if (hnodtab[hnodnum].vnodnum == vnodnum) {  /* If node already present in frontier array */
        errorPrint ("vmeshSeparateFm: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH2 */
    }

    hnodnbr ++;                                   /* One more node in hash table */
    hnodtab[hnodnum].vnodnum = vnodnum;           /* Insert node in hash table   */
    hnodtab[hnodnum].vnloval = vnloval;

    if ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]) == 0) { /* If single node */
      int                 vnodpart;

      vnodpart     = (ncmploaddlt > 0) ? 1 : 0;
      ncmpload2   -= vnloval;                     /* Node cannot belong to the separator */
      ncmploaddlt += (1 - 2 * vnodpart) * vnloval;
      hnodtab[hnodnum].vertpart  = vnodpart;
      hnodtab[hnodnum].ecmpsize0 = 0;
      continue;                                   /* No need to process elements of node vertex */
    }

    for (enodnum = meshptr->m.verttax[vnodnum], ecmpsize1 = 0; /* For all (at least one) element neighbors of node */
         enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
      Gnum                velmnum;
      Gnum                helmnum;

      velmnum = meshptr->m.edgetax[enodnum];
      for (helmnum = (velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; helmnum = (helmnum + 1) & hashmsk) {
        if (helmtab[helmnum].velmnum == ~0) {     /* If element not yet inserted   */
          if (helmnbr >= hashmax) {               /* If element hash table is full */
            if (vmeshSeparateFmResize (tablptr, &helmtab, &hnodtab, &movetab, savenbr, &lockptr, NULL, hashmax) != 0) {
              errorPrint   ("vmeshSeparateFm: cannot resize arrays (1)");
              memFree      (helmtab);             /* Free group leader */
              gainTablExit (tablptr);
              return       (1);
            }
            hashmax <<= 1;
            hashsiz <<= 1;
            hashmsk   = (hashmsk << 1) | 1;
#ifdef SCOTCH_DEBUG_VMESH3
            if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
              errorPrint ("vmeshSeparateFm: internal error (3)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VMESH3 */

            for (helmnum = (velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; /* Re-compute position in table */
                 helmtab[helmnum].velmnum != ~0; helmnum = (helmnum + 1) & hashmsk) ;
          }
          helmtab[helmnum].gainlink.prev = (GainLink *) lockptr; /* Link it */
          lockptr = &helmtab[helmnum];
          helmtab[helmnum].velmnum  = velmnum;    /* Insert it */
          helmtab[helmnum].vertpart = meshptr->parttax[velmnum] & 1; /* Separator elements to part 0 */
          helmnbr ++;
          break;
        }
        if (helmtab[helmnum].velmnum == velmnum)  /* If element already or newly inserted */
          break;                                  /* It will already be processed later   */
      }
      ecmpsize1 += helmtab[helmnum].vertpart;     /* Account for its (possibly modified) part */
    }
    hnodtab[hnodnum].vertpart  = 2;               /* Assume node is in separator */
    hnodtab[hnodnum].ecmpsize0 = meshptr->m.vendtax[vnodnum] -
                                 meshptr->m.verttax[vnodnum] - ecmpsize1;
    if (hnodtab[hnodnum].ecmpsize0 == 0) {        /* If all neighboring elements are in part 1 */
      ncmpload2   -= vnloval;                     /* Node moves from separator to part 1 too   */
      ncmploaddlt -= vnloval;
      hnodtab[hnodnum].vertpart = 1;
    }
    else if (ecmpsize1 == 0) {                    /* If all neighboring elements are in part 0 */
      ncmpload2   -= vnloval;                     /* Node moves from separator to part 0 too   */
      ncmploaddlt += vnloval;
      hnodtab[hnodnum].vertpart = 0;
    }
  }
  for (velmptr = lockptr; velmptr != NULL;        /* Process all frontier elements */
       velmptr = (VmeshSeparateFmElement *) velmptr->gainlink.prev) {
    Gnum                velmnum;
    Gnum                eelmnum;
    Gnum                ncmpcut2;
    Gnum                ncmpgain2;
    Gnum                ncmpgaindlt;

    velmnum     = velmptr->velmnum;
    ncmpcut2    =
    ncmpgain2   =
    ncmpgaindlt = 0;
    for (eelmnum = meshptr->m.verttax[velmnum];   /* For all neighbors of element */
         eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
      Gnum                vnodnum;
      Gnum                hnodnum;
      Gnum                vnoddeg;

      vnodnum = meshptr->m.edgetax[eelmnum];
      vnoddeg = meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum];
      for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; hnodnum = (hnodnum + 1) & hashmsk) {
        if (hnodtab[hnodnum].vnodnum == vnodnum) { /* If node exists (can be in same part or separator) */
          int                 vnodpart;
          Gnum                vnloval;

          vnodpart = hnodtab[hnodnum].vertpart;
          vnloval  = hnodtab[hnodnum].vnloval;
          if (vnodpart == 2) {                    /* If vertex is in separator  */
            ncmpcut2 ++;                          /* One more node in separator */
            if ((hnodtab[hnodnum].ecmpsize0 - 1) == ((vnoddeg - 2) * velmptr->vertpart)) { /* If element is only neighbor in its part */
              ncmpgain2   -= vnloval;
              ncmpgaindlt += vnloval * (2 * velmptr->vertpart - 1);
            }
          }
          else {                                  /* Vertex not in separator */
#ifdef SCOTCH_DEBUG_VMESH2
            if (vnodpart != velmptr->vertpart) {  /* Node should be in same part as element */
              errorPrint ("vmeshSeparateFm: internal error (4)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VMESH2 */
            if (vnoddeg <= 1)                     /* If element is node's sole neighbor */
              ncmpgaindlt += (vnloval * (2 * vnodpart - 1)) * 2;
            else {
              ncmpgain2   += vnloval;
              ncmpgaindlt += (vnloval * (2 * vnodpart - 1));
            }
          }
          break;
        }
        if (hnodtab[hnodnum].vnodnum == ~0) {     /* If node does not exist */
          int                 vnodpart;
          Gnum                vnloval;

          vnodpart = velmptr->vertpart;           /* Get its part */
          vnloval  = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];
#ifdef SCOTCH_DEBUG_VMESH2
          if (vnodpart != meshptr->parttax[vnodnum]) { /* Node should be in same part as element */
            errorPrint ("vmeshSeparateFm: internal error (5)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_VMESH2 */
          if (vnoddeg > 1) {                      /* If node will move to separator */
            ncmpgain2   += vnloval;               /* Increase size of separator     */
            ncmpgaindlt += (2 * vnodpart - 1) * vnloval; /* Account for imbalance   */
          }
          else                                    /* Node will move with element */
            ncmpgaindlt += (2 * vnodpart - 1) * 2 * vnloval; /* Double imbalance */
          break;
        }
      }
    }
    velmptr->ncmpcut2    = ncmpcut2;
    velmptr->ncmpgain2   = ncmpgain2;
    velmptr->ncmpgaindlt = ncmpgaindlt;
  }
  ncmploaddltmax = MAX (ncmploaddltmat, abs (ncmploaddlt)); /* Set current maximum imbalance after cleaning */

#ifdef SCOTCH_DEBUG_VMESH3
  if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
    errorPrint ("vmeshSeparateFm: internal error (6)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH3 */

  mswpnum        = 0;                             /* First sweep       */
  ncmpload2bst   = ncmpload2;                     /* Record best state */
  ncmploaddltbst = ncmploaddlt;
  do {                                            /* As long as there are improvements */
    VmeshSeparateFmElement *  velmptr;
    Gnum                      velmgain2;
    Gnum                      velmgaindlt;

    while (lockptr != NULL) {                     /* For all elements in locked list */
      VmeshSeparateFmElement *  velmptr;

      velmptr = lockptr;                          /* Unlink element from list */
      lockptr = (VmeshSeparateFmElement *) velmptr->gainlink.prev;

      velmptr->gainlink.next = VMESHSEPAFMSTATEFREE;/* Set it free anyway                   */
      if (velmptr->ncmpcut2 > 0)                  /* If element has nodes in separator      */
        gainTablAdd (tablptr, (GainLink *) velmptr, velmptr->ncmpgain2); /* Put it in table */
    }

/* fprintf (stderr, "LOOP %ld\t(%ld,\t%ld)\n", (long) passnbr, (long) ncmpload2bst, (long) ncmploaddltbst); */

    moveflag = 0;                                 /* No moves to date                          */
    movenbr  = 0;                                 /* No uneffective moves yet                  */
    while ((movenbr < paraptr->movenbr) &&        /* As long as we can find effective elements */
           ((velmptr = vmeshSeparateFmTablGet (tablptr, ncmploaddlt, ncmploaddltmax)) != NULL)) {
      VmeshSeparateFmElement *  sepaptr;          /* Linked list of separator frontier elements */
      Gnum                      velmnum;          /* Number of current element                  */
      Gnum                      velmpart;         /* Old part of current element                */
      Gnum                      eelmnum;

      gainTablDel (tablptr, &velmptr->gainlink);  /* Remove it from table */
      velmptr->gainlink.next = VMESHSEPAFMSTATEUSED; /* Mark it as used   */
      velmptr->gainlink.prev = (GainLink *) lockptr; /* Lock it           */
      lockptr                = velmptr;

      if (velmptr->mswpnum != mswpnum) {          /* If element data not yet recorded */
        movetab[savenbr].hertnum               = velmptr - helmtab; /* Record them    */
        movetab[savenbr].data.elem.vertpart    = velmptr->vertpart;
        movetab[savenbr].data.elem.ncmpcut2    = velmptr->ncmpcut2;
        movetab[savenbr].data.elem.ncmpgain2   = velmptr->ncmpgain2;
        movetab[savenbr].data.elem.ncmpgaindlt = velmptr->ncmpgaindlt;
        velmptr->mswpnum = mswpnum;
        savenbr ++;                               /* One more move recorded */
      }
      movenbr ++;                                 /* One more assumed uneffective move performed */

      velmgain2    = velmptr->ncmpgain2;          /* Save old gains for this vertex */
      velmgaindlt  = velmptr->ncmpgaindlt;
      ncmpload2   += velmgain2;                   /* Account for gains */
      ncmploaddlt += velmgaindlt;

      velmnum           = velmptr->velmnum;       /* Move element to other part */
      velmpart          = velmptr->vertpart;
      velmptr->vertpart = velmpart ^ 1;

      sepaptr = NULL;                             /* No frontier elements to relink yet */
      for (eelmnum = meshptr->m.verttax[velmnum]; /* (Re-)link neighbors              */
           eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
        Gnum                vnoddeg;
        Gnum                vnodnum;
        Gnum                hnodnum;
        Gnum                enodnum;
        Gnum                vnloval;              /* Load of current node     */
        int                 vnodpartold;          /* Old part of current node */
        int                 vnodpartnew;          /* New part of current node */
        Gnum                ecmpsize0old;
        Gnum                ecmpsize0new;

        vnodnum = meshptr->m.edgetax[eelmnum];
        vnoddeg = meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum];
        for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; hnodnum = (hnodnum + 1) & hashmsk) {
          if (hnodtab[hnodnum].vnodnum == vnodnum) /* If node found */
            break;
          if (hnodtab[hnodnum].vnodnum == ~0) {   /* If node not yet inserted */
#ifdef SCOTCH_DEBUG_VMESH2
            if (meshptr->parttax[vnodnum] != velmpart) {
              errorPrint ("vmeshSeparateFm: internal error (7)");
              return     (1);
            }
#endif /* SCOTCH_DEBUG_VMESH2 */
            if (hnodnbr >= hashmax) {             /* If node hash table is full */
              if (vmeshSeparateFmResize (tablptr, &helmtab, &hnodtab, &movetab, savenbr, &lockptr, &sepaptr, hashmax) != 0) {
                errorPrint   ("vmeshSeparateFm: cannot resize arrays (2)");
                memFree      (helmtab);           /* Free group leader */
                gainTablExit (tablptr);
                return       (1);
              }
              hashmax <<= 1;
              hashsiz <<= 1;
              hashmsk   = (hashmsk << 1) | 1;
#ifdef SCOTCH_DEBUG_VMESH3
              if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
                errorPrint ("vmeshSeparateFm: internal error (8)");
                return     (1);
              }
#endif /* SCOTCH_DEBUG_VMESH3 */
              for (helmnum = (velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; /* Re-compute positions in tables */
                   helmtab[helmnum].velmnum != velmnum; helmnum = (helmnum + 1) & hashmsk) ;
              velmptr = helmtab + helmnum;
              for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk;
                   hnodtab[hnodnum].vnodnum != ~0; hnodnum = (hnodnum + 1) & hashmsk) ;
            }

            hnodtab[hnodnum].vnodnum   = vnodnum;  /* Insert node in hash table */
            hnodtab[hnodnum].vnloval   = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];
            hnodtab[hnodnum].ecmpsize0 = vnoddeg * (1 - velmpart);
            hnodtab[hnodnum].vertpart  = velmpart; /* Node belongs to old part */
            hnodnbr ++;                           /* One more node created     */
            break;
          }
        }

        if (hnodtab[hnodnum].mswpnum != mswpnum) { /* If node data not yet recorded */
          movetab[savenbr].hertnum             = -1 - hnodnum;
          movetab[savenbr].data.node.vertpart  = hnodtab[hnodnum].vertpart;
          movetab[savenbr].data.node.ecmpsize0 = hnodtab[hnodnum].ecmpsize0;
          hnodtab[hnodnum].mswpnum = mswpnum;
          savenbr ++;                             /* One more move recorded */
        }

        vnloval = hnodtab[hnodnum].vnloval;
        if (vnoddeg <= 1) {                       /* If node only has one neighbor */
#ifdef SCOTCH_DEBUG_VMESH2
          if (hnodtab[hnodnum].vertpart != velmpart) {
            errorPrint ("vmeshSeparateFm: internal error (9)");
            return     (1);
          }
#endif /* SCOTCH_DEBUG_VMESH2 */
          hnodtab[hnodnum].vertpart  = 1 - velmpart; /* Directly move node to other part */
          hnodtab[hnodnum].ecmpsize0 = velmpart;
          continue;                               /* Skip to next node */
        }

        ecmpsize0old = hnodtab[hnodnum].ecmpsize0;
        ecmpsize0new =
        hnodtab[hnodnum].ecmpsize0 += (2 * velmpart - 1); /* One less neighbor element for this node */
#ifdef SCOTCH_DEBUG_VMESH2
        if ((hnodtab[hnodnum].ecmpsize0 < 0) ||
            (hnodtab[hnodnum].ecmpsize0 > vnoddeg)) {
          errorPrint ("vmeshSeparateFm: internal error (10)");
          return     (1);
        }
        if (hnodtab[hnodnum].vertpart == (1 - velmpart)) {
          errorPrint ("vmeshSeparateFm: internal error (11)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_VMESH2 */

        vnodpartold = hnodtab[hnodnum].vertpart;  /* Get old node part value */
        vnodpartnew = 2;                          /* Assume new node part    */
        if (vnodpartold != vnodpartnew)           /* If belonged to old part */
          hnodtab[hnodnum].vertpart = vnodpartnew; /* Move to separator      */
        else if ((ecmpsize0old - 1) == (vnoddeg - 2) * velmpart) {
          vnodpartnew =                           /* Belonged to separator and last element in this part */
          hnodtab[hnodnum].vertpart = 1 - velmpart;
        }

        for (enodnum = meshptr->m.verttax[vnodnum]; /* For all element neighbors of node */
             enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
          Gnum                velmend;
          Gnum                helmend;
          int                 vendpart;
          Gnum                ncmpcut2;
          Gnum                ncmpgain2;
          Gnum                ncmpgaindlt;

          velmend = meshptr->m.edgetax[enodnum];
          for (helmend = (velmend * VMESHSEPAFMHASHPRIME) & hashmsk; ; helmend = (helmend + 1) & hashmsk) {
            if (helmtab[helmend].velmnum == velmend) /* If element found */
              break;
            if (helmtab[helmend].velmnum == ~0) { /* If element not yet inserted */
              Gnum                ncmpgain2;
              Gnum                ncmpgaindlt;

#ifdef SCOTCH_DEBUG_VMESH2
              if (vnodpartold == 2) {             /* Elements neighboring the frontier should exist */
                errorPrint ("vmeshSeparateFm: internal error (12)");
                return     (1);
              }
              if (vnodpartold != meshptr->parttax[velmend]) { /* Unexisting elements should be in same part as their neighboring nodes */
                errorPrint ("vmeshSeparateFm: internal error (13)");
                return     (1);
              }
#endif /* SCOTCH_DEBUG_VMESH2 */
              if (helmnbr >= hashmax) {           /* If element hash table is full */
                if (vmeshSeparateFmResize (tablptr, &helmtab, &hnodtab, &movetab, savenbr, &lockptr, &sepaptr, hashmax) != 0) {
                  errorPrint   ("vmeshSeparateFm: cannot resize arrays (3)");
                  memFree      (helmtab);         /* Free group leader */
                  gainTablExit (tablptr);
                  return       (1);
                }
                hashmax <<= 1;
                hashsiz <<= 1;
                hashmsk   = (hashmsk << 1) | 1;
#ifdef SCOTCH_DEBUG_VMESH3
                if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
                  errorPrint ("vmeshSeparateFm: internal error (14)");
                  return     (1);
                }
#endif /* SCOTCH_DEBUG_VMESH3 */
                for (helmnum = (velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; /* Re-compute positions in tables */
                     helmtab[helmnum].velmnum != velmnum; helmnum = (helmnum + 1) & hashmsk) ;
                velmptr = helmtab + helmnum;
                for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk;
                     hnodtab[hnodnum].vnodnum != vnodnum; hnodnum = (hnodnum + 1) & hashmsk) ;
                for (helmend = (velmend * VMESHSEPAFMHASHPRIME) & hashmsk;
                     helmtab[helmend].velmnum != ~0; helmend = (helmend + 1) & hashmsk) ;
              }

              helmtab[helmend].gainlink.next = VMESHSEPAFMSTATEFREE;
              helmtab[helmend].velmnum       = velmend;
              helmtab[helmend].vertpart      = vnodpartold;
              helmtab[helmend].ncmpcut2      = 0;

              if (meshptr->m.vnlotax == NULL) {
                Gnum                eelmend;

                ncmpgain2   = meshptr->m.vendtax[velmend] - meshptr->m.verttax[velmend];
                ncmpgaindlt = (2 * vnodpartold - 1) * ncmpgain2;

                for (eelmend = meshptr->m.verttax[velmend]; /* For all neighboring nodes */
                     eelmend < meshptr->m.vendtax[velmend]; eelmend ++) {
                  Gnum                vnodend;

                  vnodend = meshptr->m.edgetax[eelmend];
                  if ((meshptr->m.vendtax[vnodend] - meshptr->m.verttax[vnodend]) <= 1) { /* If node linked to element only */
                    ncmpgain2   --;               /* Node will directly move to other part                                  */
                    ncmpgaindlt += (2 * velmpart - 1);
                  }
                }
              }
              else {
                Gnum                eelmend;
                Gnum                veloend;

                for (eelmend = meshptr->m.verttax[velmend], ncmpgain2 = ncmpgaindlt = veloend = 0; /* For all neighboring nodes */
                     eelmend < meshptr->m.vendtax[velmend]; eelmend ++) {
                  Gnum                vnodend;
                  Gnum                vnloend;

                  vnodend  = meshptr->m.edgetax[eelmend];
                  vnloend  = meshptr->m.vnlotax[vnodend];
                  veloend += vnloend;
                  if ((meshptr->m.vendtax[vnodend] - meshptr->m.verttax[vnodend]) <= 1) { /* If node linked to element only */
                    ncmpgain2   -= vnloend;
                    ncmpgaindlt += vnloend * (2 * velmpart - 1);
                  }
                }
                ncmpgain2   += veloend;
                ncmpgaindlt += (2 * vnodpartold - 1) * veloend;
              }
              helmtab[helmend].ncmpgain2   = ncmpgain2;
              helmtab[helmend].ncmpgaindlt = ncmpgaindlt;

              helmnbr ++;
              break;
            }
          }

          if (helmtab[helmend].mswpnum != mswpnum) { /* If element data not yet recorded */
            movetab[savenbr].hertnum               = helmend;
            movetab[savenbr].data.elem.vertpart    = helmtab[helmend].vertpart;
            movetab[savenbr].data.elem.ncmpcut2    = helmtab[helmend].ncmpcut2;
            movetab[savenbr].data.elem.ncmpgain2   = helmtab[helmend].ncmpgain2;
            movetab[savenbr].data.elem.ncmpgaindlt = helmtab[helmend].ncmpgaindlt;
            helmtab[helmend].mswpnum = mswpnum;
            savenbr ++;                           /* One more move recorded */
          }

          if (helmtab[helmend].gainlink.next != VMESHSEPAFMSTATEUSED) { /* If element available */
            if (helmtab[helmend].gainlink.next >= VMESHSEPAFMSTATELINK) /* If element linked    */
              gainTablDel (tablptr, &helmtab[helmend].gainlink); /* Unlink element              */
            helmtab[helmend].gainlink.next = VMESHSEPAFMSTATEUSED; /* Chain neighbor elements   */
            helmtab[helmend].gainlink.prev = (GainLink *) sepaptr;
            sepaptr                        = &helmtab[helmend];
          }

          vendpart    = helmtab[helmend].vertpart;
          ncmpcut2    = helmtab[helmend].ncmpcut2; /* Get element values */
          ncmpgain2   = helmtab[helmend].ncmpgain2;
          ncmpgaindlt = helmtab[helmend].ncmpgaindlt;
          if (vnodpartold != 2) {                 /* If node was in same part as the element */
            ncmpgain2   -= vnloval;
            ncmpgaindlt -= (2 * vendpart - 1) * vnloval;
          }
          else {                                  /* If node was in separator */
            ncmpcut2 --;
            if ((ecmpsize0old - 1) == ((vnoddeg - 2) * vendpart)) { /* If element was the only one in its part */
              ncmpgain2   += vnloval;
              ncmpgaindlt -= (2 * vendpart - 1) * vnloval;
            }
          }
          if (vnodpartnew != 2) {                 /* If node is now in same part as the element */
            ncmpgain2   += vnloval;
            ncmpgaindlt += (2 * vendpart - 1) * vnloval;
          }
          else {                                  /* If node is now in separator */
            ncmpcut2 ++;
            if ((ecmpsize0new - 1) == ((vnoddeg - 2) * vendpart)) { /* If element is the only one in its part */
              ncmpgain2   -= vnloval;
              ncmpgaindlt += (2 * vendpart - 1) * vnloval;
            }
          }
          helmtab[helmend].ncmpcut2    = ncmpcut2; /* Adjust element values */
          helmtab[helmend].ncmpgain2   = ncmpgain2;
          helmtab[helmend].ncmpgaindlt = ncmpgaindlt;
        }
      }
      velmptr->ncmpgain2   = - velmgain2;         /* Set new gains of element */
      velmptr->ncmpgaindlt = - velmgaindlt;

      while (sepaptr != NULL) {                   /* As long as there are element to re-link */
        VmeshSeparateFmElement *  velmptr;

        velmptr = sepaptr;                        /* Get element to re-link */
        sepaptr = (VmeshSeparateFmElement *) velmptr->gainlink.prev;
        velmptr->gainlink.next = VMESHSEPAFMSTATEFREE;
        if (velmptr->ncmpcut2 != 0)               /* If element belongs to frontier      */
          gainTablAdd (tablptr, (GainLink *) velmptr, velmptr->ncmpgain2); /* Re-link it */
      }

#ifdef SCOTCH_DEBUG_VMESH3
      if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
        errorPrint ("vmeshSeparateFm: internal error (15)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH3 */

      if (ncmpload2 < ncmpload2bst) {             /* If move improves separator size */
        ncmpload2bst   = ncmpload2;               /* This move was effective         */
        ncmploaddltbst = ncmploaddlt;
        movenbr  =
        savenbr  = 0;
        moveflag = 1;
        mswpnum ++;
      } else if (ncmpload2 == ncmpload2bst) {
        if (abs (ncmploaddlt) < abs (ncmploaddltbst)) {
          ncmploaddltbst = ncmploaddlt;           /* This move was effective */
          movenbr  =
          savenbr  = 0;
          moveflag = 1;
          mswpnum ++;
        }
        else if (abs (ncmploaddlt) == abs (ncmploaddltbst)) {
          ncmploaddltbst = ncmploaddlt;           /* Might be the opposite, so record */
          savenbr = 0;                            /* Forget backtracking              */
          mswpnum ++;
        }
      }
      if (ncmploaddltmax > ncmploaddltmat) {      /* If must restrict distance bounds */
        Gnum                ncmploaddlttmp;

        ncmploaddlttmp = ncmploaddltmax;          /* Save old working ncmpdltmax value */
        ncmploaddltmax = MAX (ncmploaddltmat,     /* Restrict at most to maximum       */
                              abs (ncmploaddlt));
        if (ncmploaddltmax < ncmploaddlttmp) {    /* If we have done something useful */
          ncmpload2bst   = ncmpload2;             /* Then record best move done       */
          ncmploaddltbst = ncmploaddlt;
          movenbr =
          savenbr = 0;
          mswpnum ++;
        }
      }
    }

    while (savenbr > 0) {                         /* Delete exceeding moves */
      Gnum                hertnum;
      
      hertnum = movetab[-- savenbr].hertnum;      /* Get vertex hash number */
      if (hertnum >= 0) {                         /* If vertex is element   */
        helmtab[hertnum].vertpart    = movetab[savenbr].data.elem.vertpart;
        helmtab[hertnum].ncmpcut2    = movetab[savenbr].data.elem.ncmpcut2;
        helmtab[hertnum].ncmpgain2   = movetab[savenbr].data.elem.ncmpgain2;
        helmtab[hertnum].ncmpgaindlt = movetab[savenbr].data.elem.ncmpgaindlt;
        if (helmtab[hertnum].gainlink.next != VMESHSEPAFMSTATEUSED) { /* If element not already removed  */
          if (helmtab[hertnum].gainlink.next >= VMESHSEPAFMSTATELINK) /* If vertex is still linked       */
            gainTablDel (tablptr, &helmtab[hertnum].gainlink); /* Remove it from table                   */
          helmtab[hertnum].gainlink.next = VMESHSEPAFMSTATEUSED;
          helmtab[hertnum].gainlink.prev = (GainLink *) lockptr; /* Lock it */
          lockptr                        = &helmtab[hertnum];
        }
      }
      else {                                      /* Vertex is node */
        hertnum = -1 - hertnum;                   /* Get hash index */
        hnodtab[hertnum].vertpart  = movetab[savenbr].data.node.vertpart;
        hnodtab[hertnum].ecmpsize0 = movetab[savenbr].data.node.ecmpsize0;
      }
    }
    ncmpload2   = ncmpload2bst;                   /* Restore best separator parameters */
    ncmploaddlt = ncmploaddltbst;
    mswpnum ++;                                   /* Forget all recorded moves */

#ifdef SCOTCH_DEBUG_VMESH3
    if (vmeshSeparateFmCheck (meshptr, helmtab, hnodtab, hashmsk, ncmpload2, ncmploaddlt) != 0) {
      errorPrint ("vmeshSeparateFm: internal error (16)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_VMESH3 */
  } while ((moveflag != 0) &&                     /* As long as vertices are moved                          */
           (-- passnbr != 0));                    /* And we are allowed to loop (TRICK for negative values) */

  ecmpload1 = 0;                                  /* Assume no change in elements */
  for (helmnum = 0; helmnum < hashsiz; helmnum ++) {
    Gnum                velmnum;

    velmnum = helmtab[helmnum].velmnum;
    if ((velmnum != ~0) && (helmtab[helmnum].vertpart != meshptr->parttax[velmnum])) {
#ifdef SCOTCH_DEBUG_VMESH2
      if ((helmtab[helmnum].vertpart < 0) ||      /* Separator elements should have been removed */
          (helmtab[helmnum].vertpart > 1)) {
        errorPrint ("vmeshSeparateFm: internal error (17)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH2 */
      ecmpload1 += helmtab[helmnum].vertpart - (meshptr->parttax[velmnum] & 1);
      meshptr->parttax[velmnum] = helmtab[helmnum].vertpart;
    }
  }
  meshptr->ecmpsize[1] += ecmpload1;              /* No longer elements in separator */
  meshptr->ecmpsize[0]  = meshptr->m.velmnbr - meshptr->ecmpsize[1];

#ifdef SCOTCH_DEBUG_VMESH2
  if ((meshptr->ecmpsize[0] + meshptr->ecmpsize[1]) != meshptr->m.velmnbr) { /* Separator elements should have been removed */
    errorPrint ("vmeshSeparateFm: internal error (18)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  ncmpsize1 =
  ncmpsize2 = 0;
  ncmpload1 =
  ncmpload2 = 0;
  for (hnodnum = 0, fronnum = 0; hnodnum < hashsiz; hnodnum ++) {
    Gnum                vnodnum;

    vnodnum = hnodtab[hnodnum].vnodnum;
    if (vnodnum != ~0) {
#ifdef SCOTCH_DEBUG_VMESH2
      if ((hnodtab[hnodnum].vertpart < 0) ||
          (hnodtab[hnodnum].vertpart > 2)) {
        errorPrint ("vmeshSeparateFm: internal error (19)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH2 */

      if (hnodtab[hnodnum].vertpart == 2)         /* If node belongs to separator */
        meshptr->frontab[fronnum ++] = vnodnum;   /* Add it to separator array    */
      if (hnodtab[hnodnum].vertpart != meshptr->parttax[vnodnum]) {
        Gnum                diffpart1;
        Gnum                diffpart2;

        diffpart1  = (hnodtab[hnodnum].vertpart &  1) - (meshptr->parttax[vnodnum] &  1);
        diffpart2  = (hnodtab[hnodnum].vertpart >> 1) - (meshptr->parttax[vnodnum] >> 1);
        ncmpsize1 += diffpart1;
        ncmpsize2 += diffpart2;
        ncmpload1 += hnodtab[hnodnum].vnloval * diffpart1;
        ncmpload2 += hnodtab[hnodnum].vnloval * diffpart2;
        meshptr->parttax[vnodnum] = hnodtab[hnodnum].vertpart;
      }
    }
  }
#ifdef SCOTCH_DEBUG_VMESH2
  if ((meshptr->fronnbr + ncmpsize2) != fronnum) {
    errorPrint ("vmeshSeparateFm: internal error (20)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */
  meshptr->ncmpload[1] += ncmpload1;
  meshptr->ncmpload[2] += ncmpload2;
  meshptr->ncmpload[0]  = meshptr->m.vnlosum - meshptr->ncmpload[1] - meshptr->ncmpload[2];
  meshptr->ncmploaddlt  = meshptr->ncmpload[0] - meshptr->ncmpload[1];
  meshptr->fronnbr      = fronnum;
  meshptr->ncmpsize[1] += ncmpsize1;
  meshptr->ncmpsize[0]  = meshptr->m.vnodnbr - fronnum - meshptr->ncmpsize[1];

  memFree      (helmtab);                         /* Free group leader */
  gainTablExit (tablptr);

#ifdef SCOTCH_DEBUG_VMESH2
  if (vmeshCheck (meshptr) != 0) {
    errorPrint ("vmeshSeparateFm: internal error (21)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

/* printf ("FM Sepa\tsize=%ld\tload=%ld\tbal=%ld\n", (long) meshptr->fronnbr, (long) meshptr->ncmpload[2], (long) meshptr->ncmploaddlt); */

  return (0);
}

/* This routine checks the consistency of
** the hash structures.
** It returns:
** - 0   : in case of success.
** - !0  : in case of error.
*/

#ifdef SCOTCH_DEBUG_VMESH3
static
int
vmeshSeparateFmCheck (
const Vmesh * const                     meshptr,
const VmeshSeparateFmElement * restrict helmtab,
const VmeshSeparateFmNode * restrict    hnodtab,
const Gnum                              hashmsk,
const Gnum                              ncmpload2,
const Gnum                              ncmploaddlt)
{
  Gnum                  vertnbr;
  Gnum                  helmnum;
  Gnum                  hnodnum;
  Gnum                  ncmploadtmp[3];
  GraphPart * restrict  parttax;

  vertnbr = meshptr->m.velmnbr + meshptr->m.vnodnbr;
  if ((parttax = (GraphPart *) memAlloc (vertnbr * sizeof (GraphPart))) == NULL) {
    errorPrint ("vmeshSeparateFmCheck: out of memory");
    return     (1);
  }
  memCpy (parttax, meshptr->parttax + meshptr->m.baseval, vertnbr * sizeof (GraphPart));
  parttax -= meshptr->m.baseval;

  ncmploadtmp[0] = meshptr->ncmpload[0];
  ncmploadtmp[1] = meshptr->ncmpload[1];
  ncmploadtmp[2] = meshptr->ncmpload[2];
  for (hnodnum = 0; hnodnum <= hashmsk; hnodnum ++) { /* For all node slots */
    Gnum                vnodnum;
    Gnum                enodnum;
    Gnum                ecmpsize0;
    int                 vnodpart;

    vnodnum = hnodtab[hnodnum].vnodnum;
    if (vnodnum == ~0)                            /* If unallocated slot */
      continue;                                   /* Skip to next slot   */

    if (hnodtab[hnodnum].vnloval != ((meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum])) {
      errorPrint ("vmeshSeparateFmCheck: invalid node load");
      return     (1);
    }
    if ((hnodtab[hnodnum].ecmpsize0 < 0) ||
        (hnodtab[hnodnum].ecmpsize0 > (meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]))) {
      errorPrint ("vmeshSeparateFmCheck: invalid node neighbors in part 0");
      return     (1);
    }
    vnodpart = hnodtab[hnodnum].vertpart;
    if (vnodpart != meshptr->parttax[vnodnum]) {
      ncmploadtmp[meshptr->parttax[vnodnum]] -= hnodtab[hnodnum].vnloval;
      ncmploadtmp[vnodpart]                  += hnodtab[hnodnum].vnloval;
      parttax[vnodnum]                        = vnodpart;
    }

    ecmpsize0 = 0;
    for (enodnum = meshptr->m.verttax[vnodnum];   /* For all element neighbors */
         enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
      Gnum                velmnum;
      Gnum                helmnum;
      int                 velmpart;

      velmnum = meshptr->m.edgetax[enodnum];
      for (helmnum = (velmnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; helmnum = (helmnum + 1) & hashmsk) {
        if (helmtab[helmnum].velmnum == velmnum) { /* If element found */
          velmpart         = helmtab[helmnum].vertpart;
          parttax[velmnum] = velmpart;
          break;
        }
        if (helmtab[helmnum].velmnum == ~0) {     /* If element not present */
          velmpart = meshptr->parttax[velmnum];
          break;
        }
      }
      if (velmpart == 0)
        ecmpsize0 ++;
    }

    if (ecmpsize0 != hnodtab[hnodnum].ecmpsize0) {
      errorPrint ("vmeshSeparateFmCheck: invalid node neighbor count");
      return     (1);
    }
  }
  if (ncmpload2 != ncmploadtmp[2]) {
    errorPrint ("vmeshSeparateFmCheck: invalid frontier load");
    return     (1);
  }
  if (ncmploaddlt != (ncmploadtmp[0] - ncmploadtmp[1])) {
    errorPrint ("vmeshSeparateFmCheck: invalid separator balance");
    return     (1);
  }

  for (helmnum = 0; helmnum <= hashmsk; helmnum ++) { /* For all element slots */
    Gnum                velmnum;
    Gnum                eelmnum;
    Gnum                ncmpcut2;
    Gnum                ncmpgain2;
    Gnum                ncmpgaindlt;
    Gnum                ncmpsize[3];
    int                 velmpart;

    velmnum = helmtab[helmnum].velmnum;
    if (velmnum == ~0)                            /* If unallocated slot */
      continue;                                   /* Skip to next slot   */

    ncmpcut2    =
    ncmpgain2   =
    ncmpgaindlt = 0;
    ncmpsize[0] =
    ncmpsize[1] =
    ncmpsize[2] = 0;
    velmpart    = helmtab[helmnum].vertpart;
    for (eelmnum = meshptr->m.verttax[velmnum];   /* For all node neighbors */
         eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
      Gnum                vnodnum;
      Gnum                hnodnum;
      int                 vnodpart;
      Gnum                vnloval;

      vnodnum = meshptr->m.edgetax[eelmnum];
      vnloval = (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];
      for (hnodnum = (vnodnum * VMESHSEPAFMHASHPRIME) & hashmsk; ; hnodnum = (hnodnum + 1) & hashmsk) {
        if (hnodtab[hnodnum].vnodnum == vnodnum) { /* If element found */
          vnodpart = hnodtab[hnodnum].vertpart;
          if (vnodpart != 2) {
            if (vnodpart != velmpart) {
              errorPrint ("vmeshSeparateFmCheck: invalid separator node (1)");
              return     (1);
            }
            if ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum] - 1) == 0)
              ncmpgaindlt += (2 * vnodpart - 1) * 2 * vnloval;
            else {
              ncmpgain2   += vnloval;
              ncmpgaindlt += (2 * vnodpart - 1) * vnloval;
            }
          }
          else if (((hnodtab[hnodnum].ecmpsize0 == 1) && (velmpart == 0)) ||
                   ((hnodtab[hnodnum].ecmpsize0 == (meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum] - 1)) && (velmpart == 1))) {
            ncmpgain2   -= vnloval;
            ncmpgaindlt += (2 * velmpart - 1) * vnloval;
          }
          break;
        }
        if (hnodtab[hnodnum].vnodnum == ~0) {     /* If element not present */
          vnodpart = meshptr->parttax[vnodnum];
          if (vnodpart != velmpart) {
            errorPrint ("vmeshSeparateFmCheck: invalid separator node (2)");
            return     (1);
          }
          if ((meshptr->m.vendtax[vnodnum] - meshptr->m.verttax[vnodnum]) == 1) {
            if (vnodpart == 2) {
              errorPrint ("vmeshSeparateFmCheck: invalid separator node (3)");
              return     (1);
            }
            ncmpgaindlt += (2 * vnodpart - 1) * 2 * vnloval;
          }
          else {
            ncmpgain2   += vnloval;
            ncmpgaindlt += (2 * vnodpart - 1) * vnloval;
          }
          break;
        }
      }
      ncmpsize[vnodpart] ++;
    }
    if ((ncmpsize[0] != 0) && (ncmpsize[1] != 0)) {
      errorPrint ("vmeshSeparateFmCheck: invalid element nodes");
      return     (1);
    }
    if (ncmpsize[2] != helmtab[helmnum].ncmpcut2) {
      errorPrint ("vmeshSeparateFmCheck: invalid element separator count");
      return     (1);
    }
    if ((ncmpgain2   != helmtab[helmnum].ncmpgain2) ||
        (ncmpgaindlt != helmtab[helmnum].ncmpgaindlt)) {
      errorPrint ("vmeshSeparateFmCheck: invalid element gains");
      return     (1);
    }
  }

  memFree (parttax + meshptr->m.baseval);

  return (0);
}
#endif /* SCOTCH_DEBUG_VMESH3 */
