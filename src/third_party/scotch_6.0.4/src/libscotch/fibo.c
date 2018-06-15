/* Copyright 2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : fibo.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles Fibonacci trees.    **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 may 2010     **/
/**                                 to     12 may 2010     **/
/**                # Version 6.0  : from : 22 oct 2011     **/
/**                                 to     22 oct 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define FIBO

#include "module.h"
#include "common.h"
#include "fibo.h"

/*********************************************/
/*                                           */
/* These routines deal with Fibonacci trees. */
/*                                           */
/*********************************************/

/* This routine initializes a Fibonacci
** tree structure.
** It returns:
** - 0   : in case of success.
** - !0  : on error.
*/

int
fiboTreeInit (
FiboTree * const            treeptr,
int                      (* cmpfptr) (const FiboNode * const, const FiboNode * const))
{
  if ((treeptr->degrtab = (FiboNode **) memAlloc ((sizeof (INT) << 3) * sizeof (FiboNode *))) == NULL) /* As many cells as there are bits in an INT */
    return (1);

  memSet (treeptr->degrtab, 0, (sizeof (INT) << 3) * sizeof (FiboNode *)); /* Make degree array ready for consolidation: all cells set to NULL */

  treeptr->rootdat.linkdat.prevptr =              /* Link root node to itself */
  treeptr->rootdat.linkdat.nextptr = &treeptr->rootdat;
  treeptr->cmpfptr = cmpfptr;

  return (0);
}

/* This routine flushes the contents of
** the given Fibonacci tree.
** It returns:
** - VOID  : in all cases.
*/

void
fiboTreeExit (
FiboTree * const            treeptr)
{
  if (treeptr->degrtab != NULL)
    memFree (treeptr->degrtab);
}

/* This routine flushes the contents of
** the given Fibonacci tree.
** It returns:
** - VOID  : in all cases.
*/

void
fiboTreeFree (
FiboTree * const            treeptr)
{
  treeptr->rootdat.linkdat.prevptr =              /* Link root node to itself */
  treeptr->rootdat.linkdat.nextptr = &treeptr->rootdat;
}

/* This routine perform the consolidation
** of roots per degree. It returns the best
** element found because this element is not
** recorded in the data structure itself.
** It returns:
** - !NULL  : pointer to best element found.
** - NULL   : Fibonacci tree is empty.
*/

FiboNode *
fiboTreeConsolidate (
FiboTree * const            treeptr)
{
  FiboNode ** restrict  degrtab;
  int                   degrmax;
  int                   degrval;
  FiboNode *            rootptr;
  FiboNode *            nextptr;
  FiboNode *            bestptr;

  degrtab = treeptr->degrtab;

  for (rootptr = treeptr->rootdat.linkdat.nextptr, nextptr = rootptr->linkdat.nextptr, degrmax = 0; /* For all roots in root list */
       rootptr != &treeptr->rootdat; ) {
    degrval = rootptr->deflval >> 1;              /* Get degree, getting rid of flag part */
#ifdef SCOTCH_DEBUG_FIBO2
    if (degrval >= (sizeof (INT) << 3))
      errorPrint ("fiboTreeConsolidate: invalid node degree");
#endif /* SCOTCH_DEBUG_FIBO2 */
    if (degrtab[degrval] == NULL) {               /* If no tree with same degree already found */
      if (degrval > degrmax)                      /* Record highest degree found               */
        degrmax = degrval;

      degrtab[degrval] = rootptr;                 /* Record tree as first tree with this degree      */
      rootptr = nextptr;                          /* Process next root in list during next iteration */
      nextptr = rootptr->linkdat.nextptr;
    }
    else {
      FiboNode *            oldrptr;              /* Root which will no longer be a root */
      FiboNode *            chldptr;

      oldrptr = degrtab[degrval];                 /* Assume old root is worse           */
      if (treeptr->cmpfptr (oldrptr, rootptr) <= 0) { /* If old root is still better    */
        oldrptr = rootptr;                        /* This root will be be linked to it  */
        rootptr = degrtab[degrval];               /* We will go on processing this root */
      }

      degrtab[degrval] = NULL;                    /* Remaining root changes degree so leaves this cell */
      fiboTreeUnlink (oldrptr);                   /* Old root is no longer a root                      */
      oldrptr->deflval &= ~1;                     /* Whatever old root flag was, it is reset to 0      */
      oldrptr->pareptr = rootptr;                 /* Remaining root is now father of old root          */

      chldptr = rootptr->chldptr;                 /* Get first child of remaining root                                    */
      if (chldptr != NULL) {                      /* If remaining root had already some children, link old root with them */
        rootptr->deflval += 2;                    /* Increase degree by 1, that is, by 2 with left shift in deflval       */
        fiboTreeLinkAfter (chldptr, oldrptr);
      }
      else {                                      /* Old root becomes first child of remaining root */
        rootptr->deflval = 2;                     /* Real degree set to 1, and flag set to 0        */
        rootptr->chldptr = oldrptr;
        oldrptr->linkdat.prevptr =                /* Chain old root to oneself as only child */
        oldrptr->linkdat.nextptr = oldrptr;
      }
    }                                             /* Process again remaining root as its degree has changed */
  }

  bestptr = NULL;
  for (degrval = 0; degrval <= degrmax; degrval ++) {
    if (degrtab[degrval] != NULL) {               /* If some tree is found           */
      bestptr = degrtab[degrval];                 /* Record it as potential best     */
      degrtab[degrval] = NULL;                    /* Clean-up used part of array     */
      degrval ++;                                 /* Go on at next cell in next loop */
      break;
    }
  }
  for ( ; degrval <= degrmax; degrval ++) {       /* For remaining roots once a potential best root has been found */
    if (degrtab[degrval] != NULL) {
      if (treeptr->cmpfptr (degrtab[degrval], bestptr) < 0) /* If new root is better */
        bestptr = degrtab[degrval];               /* Record new root as best root    */
      degrtab[degrval] = NULL;                    /* Clean-up used part of array     */
    }
  }

  return (bestptr);
}

/* This routine returns the node of minimum
** key in the given tree. The node is searched
** for each time this routine is called, so this
** information should be recorded if needed.
** This is the non-macro version, for testing
** and setting up breakpoints.
** It returns:
** - !NULL  : pointer to best element found.
** - NULL   : Fibonacci tree is empty.
*/

#ifndef fiboTreeMinIsMacro

FiboNode *
fiboTreeMin (
FiboTree * const            treeptr)
{
  FiboNode *            bestptr;

  bestptr = fiboTreeMinMacro (treeptr);

#ifdef SCOTCH_DEBUG_FIBO3
  fiboTreeCheck (treeptr);
#endif /* SCOTCH_DEBUG_FIBO3 */

  return (bestptr);
}

#endif /* fiboTreeMinIsMacro */

/* This routine adds the given node to the
** given tree. This is the non-macro version,
** for testing and setting up breakpoints.
** It returns:
** - void  : in all cases.
*/

#ifndef fiboTreeAddIsMacro

void
fiboTreeAdd (
FiboTree * const            treeptr,
FiboNode * const            nodeptr)
{
  fiboTreeAddMacro (treeptr, nodeptr);

#ifdef SCOTCH_DEBUG_FIBO3
  fiboTreeCheck (treeptr);
#endif /* SCOTCH_DEBUG_FIBO3 */
}

#endif /* fiboTreeAddIsMacro */

/* This routine deletes the given node from
** the given tree, whatever ths node is (root
** or non root). This is the non-macro version,
** for testing and setting up breakpoints.
** It returns:
** - void  : in all cases.
*/

#ifndef fiboTreeDelIsMacro

void
fiboTreeDel (
FiboTree * const            treeptr,
FiboNode * const            nodeptr)
{
  fiboTreeDelMacro (treeptr, nodeptr);

#ifdef SCOTCH_DEBUG_FIBO3
  nodeptr->pareptr =
  nodeptr->chldptr =
  nodeptr->linkdat.prevptr =
  nodeptr->linkdat.nextptr = NULL;

  fiboTreeCheck (treeptr);
#endif /* SCOTCH_DEBUG_FIBO3 */
}

#endif /* fiboTreeDelIsMacro */

/* This routine checks the consistency of the
** given linked list.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no such vertex available.
*/

#ifdef SCOTCH_DEBUG_FIBO3

static
int
fiboTreeCheck2 (
const FiboNode * const      nodeptr)
{
  FiboNode *            chldptr;
  int                   degrval;

  degrval = 0;
  chldptr = nodeptr->chldptr;
  if (chldptr != NULL) {
    do {
      if (chldptr->linkdat.nextptr->linkdat.prevptr != chldptr) {
        errorPrint ("fiboTreeCheck: bad child linked list");
        return     (1);
      }

      if (chldptr->pareptr != nodeptr) {
        errorPrint ("fiboTreeCheck: bad child parent");
        return (1);
      }

      if (fiboTreeCheck2 (chldptr) != 0)
        return (1);

      degrval ++;
      chldptr = chldptr->linkdat.nextptr;
    } while (chldptr != nodeptr->chldptr);
  }

  if (degrval != (nodeptr->deflval >> 1)) {       /* Real node degree is obtained by discarding lowest bit */
    errorPrint ("fiboTreeCheck2: invalid child information");
    return     (1);
  }

  return (0);
}

int
fiboTreeCheck (
const FiboTree * const      treeptr)
{
  FiboNode *            nodeptr;

  for (nodeptr = treeptr->rootdat.linkdat.nextptr;
       nodeptr != &treeptr->rootdat; nodeptr = nodeptr->linkdat.nextptr) {
    if (nodeptr->linkdat.nextptr->linkdat.prevptr != nodeptr) {
      errorPrint ("fiboTreeCheck: bad root linked list");
      return     (1);
    }

    if (nodeptr->pareptr != NULL) {
      errorPrint ("fiboTreeCheck: bad root parent");
      return (1);
    }

    if (fiboTreeCheck2 (nodeptr) != 0)
      return (1);
  }

  return (0);
}

#endif /* SCOTCH_DEBUG_FIBO3 */
