/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : gain.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles logarithmic gain    **/
/**                table structures.                       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 26 oct 1996     **/
/**                                 to     30 nov 1996     **/
/**                # Version 0.1  : from : 10 may 1999     **/
/**                                 to     10 may 1999     **/
/**                # Version 4.0  : from : 10 jan 2004     **/
/**                                 to     18 mar 2005     **/
/**                # Version 5.0  : from : 24 mar 2008     **/
/**                                 to     24 mar 2008     **/
/**                                                        **/
/**   NOTES      : # Most of the contents of this module   **/
/**                  comes from "map_b_fm" of the SCOTCH   **/
/**                  project (v3.2).                       **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GAIN

#include "module.h"
#include "common.h"
#include "gain.h"

/*
**  The static variables.
*/

static GainLink             gainLinkDummy;        /*+ Dummy link space for fast linked list operations +*/

/*************************************************/
/*                                               */
/* These routines deal with generic gain tables. */
/*                                               */
/*************************************************/

/* This routine allocates and initializes
** a gain table structure with the proper
** number of subbits.
** It returns:
** - !NULL  : pointer to the gain table;
** - NULL   : on error.
*/

GainTabl *
gainTablInit (
const INT                   gainmax,
const INT                   subbits)
{
  GainEntr *          entrptr;
  GainTabl *          tablptr;
  INT                 totsize;

  if (gainmax >= GAIN_LINMAX) {                   /* If logarithmic indexing           */
    totsize = ((sizeof (INT) << 3) - subbits) << (subbits + 1); /* Allocate gain table */

    if ((tablptr = (GainTabl *) memAlloc (sizeof (GainTabl) + (totsize - 1) * sizeof (GainEntr))) == NULL)
      return (NULL);

    tablptr->tablAdd = gainTablAddLog;

    tablptr->subbits = subbits;                   /* Fill gain table fields                    */
    tablptr->submask = (1 << (subbits + 1)) - 1;  /* Mask with all subbits, plus one, set to 1 */
  }
  else {                                          /* Linear indexing     */
    totsize = 2 * GAIN_LINMAX;                    /* Allocate gain table */

    if ((tablptr = (GainTabl *) memAlloc (sizeof (GainTabl) + (totsize - 1) * sizeof (GainEntr))) == NULL)
      return (NULL);

    tablptr->tablAdd = gainTablAddLin;
    tablptr->subbits = 0;                         /* Fill gain table fields */
    tablptr->submask = 0;
  }

  tablptr->totsize = totsize;
  tablptr->tabl    = tablptr->tabk + (totsize / 2);
  tablptr->tend    = tablptr->tabk + (totsize - 1); /* End of gain entry array */
  tablptr->tmin    = tablptr->tend;               /* Entries of extremal gain  */
  tablptr->tmax    = tablptr->tabk;

  for (entrptr  = tablptr->tabk;                  /* Initialize gain table entries */
       entrptr <= tablptr->tend;
       entrptr ++)
    entrptr->next = &gainLinkDummy;               /* Point to dummy link area */

  return (tablptr);
}

/* This routine deletes a gain list
** It returns:
** - VOID  : in all cases.
*/

void
gainTablExit (
GainTabl * const         tablptr)
{
  memFree (tablptr);                              /* Free table structure itself */
}

/* This routine flushes the contents of
** the given gain table.
** It returns:
** - VOID  : in all cases.
*/

void
gainTablFree (
GainTabl * const            tablptr)
{
  GainEntr *          entrptr;

  for (entrptr  = tablptr->tmin;                  /* Flush only used area */
       entrptr <= tablptr->tmax;
       entrptr ++)
    entrptr->next = &gainLinkDummy;               /* Point to dummy link area */

  tablptr->tmin = tablptr->tend;                  /* Entries of extremal gain */
  tablptr->tmax = tablptr->tabk;
}

/* This routine adds a vertex to the table
** and table gain indicated in the vertex
** fields.
** It returns:
** - VOID  : in all cases.
*/

void
gainTablAddLin (
GainTabl * const            tablptr,              /*+ Pointer to gain table   +*/
GainLink * const            linkptr,              /*+ Pointer to entry to add +*/
const INT                   gain)                 /*+ Gain value              +*/
{
  GainEntr *          entrptr;                    /* Pointer to gain entry   */
  GainLink *          headptr;                    /* Pointer to head of list */

#ifdef SCOTCH_DEBUG_GAIN2
  if (tablptr->tablAdd != gainTablAddLin) {
    errorPrint ("gainTablAddLin: table type mismatch");
    return;
  }
#endif /* SCOTCH_DEBUG_GAIN2 */

  entrptr = tablptr->tabl + gain;
  if (entrptr < tablptr->tabk)
    entrptr = tablptr->tabk;
  else if (entrptr > tablptr->tend)
    entrptr = tablptr->tend;

  if (entrptr < tablptr->tmin)
    tablptr->tmin = entrptr;
  if (entrptr > tablptr->tmax)
    tablptr->tmax = entrptr;

  headptr = (GainLink *) entrptr;                 /* TRICK: assume gain entry is a link */
  linkptr->tabl       = entrptr;                  /* Set table position                 */
  headptr->next->prev = linkptr;                  /* Link vertex in gain list: TRICK    */
  linkptr->prev       = headptr;
  linkptr->next       = headptr->next;
  headptr->next       = linkptr;
}

/* This routine adds a vertex to the table
** and table gain indicated in the vertex
** fields.
** It returns:
** - VOID  : in all cases.
*/

void
gainTablAddLog (
GainTabl * const            tablptr,              /*+ Pointer to gain table   +*/
GainLink * const            linkptr,              /*+ Pointer to entry to add +*/
const INT                   gain)                 /*+ Gain value              +*/
{
  GainEntr *          entrptr;                    /* Pointer to gain entry   */
  INT                 i, j;

#ifdef SCOTCH_DEBUG_GAIN2
  if (tablptr->tablAdd != gainTablAddLog) {
    errorPrint ("gainTablAddLog: table type mismatch");
    return;
  }
#endif /* SCOTCH_DEBUG_GAIN2 */

  if (gain >= 0) {                                /* Compute table entry for gain */
    for (i = 0, j = gain; j > tablptr->submask; i ++, j >>= 1) ;
    i = (i << tablptr->subbits) + j;
  }
  else {
    for (i = 0, j = - (gain + 1); j > tablptr->submask; i ++, j >>= 1) ;
    i = - ((i << tablptr->subbits) + j + 1);
  }
  entrptr = tablptr->tabl + i;

  if (entrptr < tablptr->tmin)
    tablptr->tmin = entrptr;
  if (entrptr > tablptr->tmax)
    tablptr->tmax = entrptr;

#ifdef SCOTCH_DEBUG_GAIN3
  if ((entrptr->next != &gainLinkDummy) &&
      (entrptr->next->prev != (GainLink *) entrptr)) {
    errorPrint ("gainTablAddLog: bad first element");
    return;
  }
  if (gainTablCheck (entrptr) != 0) {
    errorPrint ("gainTablAddLog: bad chaining");
  }
#endif /* SCOTCH_DEBUG_GAIN3 */

  entrptr->next->prev = linkptr;                  /* Link vertex in gain list: TRICK */
  linkptr->prev       = (GainLink *) entrptr;
  linkptr->next       = entrptr->next;
  linkptr->tabl       = entrptr;                  /* Set table position */
  entrptr->next       = linkptr;
}

/* This routine removes a link
** from the table.
** It returns:
** - VOID  : in all cases.
*/

#ifdef SCOTCH_DEBUG_GAIN1                         /* Compiled only in debug mode */
void
gainTablDel (
GainTabl * const            tablptr,
GainLink * const            linkptr)              /*+ Pointer to link to delete +*/
{
#ifdef SCOTCH_DEBUG_GAIN3
  if (linkptr->tabl != NULL) {
    if ((linkptr->tabl->next != &gainLinkDummy) &&
        (linkptr->tabl->next->prev != (GainLink *) linkptr->tabl)) {
      errorPrint ("gainTablDel: bad first element");
      return;
    }
    if (gainTablCheck (linkptr->tabl) != 0) {
      errorPrint ("gainTablDel: bad chaining");
      return;
    }
  }
#endif /* SCOTCH_DEBUG_GAIN3 */

  linkptr->next->prev = linkptr->prev;            /* TRICK: may write in dummy link area */
  linkptr->prev->next = linkptr->next;
}
#endif /* SCOTCH_DEBUG_GAIN1 */

/* This routine returns the link of best
** gain in the table structure.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no such vertex available.
*/

GainLink *
gainTablFrst (
GainTabl * const            tablptr)
{
  GainEntr *          entrptr;

  for (entrptr  = tablptr->tmin;
       entrptr <= tablptr->tend;
       entrptr ++) {
    if (entrptr->next != &gainLinkDummy) {
      tablptr->tmin = entrptr;
#ifdef SCOTCH_DEBUG_GAIN3
      if (gainTablCheck (entrptr) != 0) {
        errorPrint ("gainTablFrst: bad chaining");
        return     (NULL);
      }
#endif /* SCOTCH_DEBUG_GAIN3 */
      return (entrptr->next);
    }
  }
  tablptr->tmin = tablptr->tend;
  tablptr->tmax = tablptr->tabk;

  return (NULL);
}

/* This routine returns the next best vertex
** following the given vertex.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no such vertex available.
*/

GainLink *
gainTablNext (
GainTabl * const            tablptr,
const GainLink * const      linkptr)
{
  GainEntr *          entrptr;

  if (linkptr->next != &gainLinkDummy) {
#ifdef SCOTCH_DEBUG_GAIN3
    if (gainTablCheck (linkptr->tabl) != 0) {
      errorPrint ("gainTablNext: bad chaining (1)");
      return     (NULL);
    }
#endif /* SCOTCH_DEBUG_GAIN3 */
    return (linkptr->next);
  }

  for (entrptr = linkptr->tabl + 1;
       entrptr < tablptr->tend;
       entrptr ++) {
    if (entrptr->next != &gainLinkDummy) {
#ifdef SCOTCH_DEBUG_GAIN3
      if (gainTablCheck (entrptr) != 0) {
        errorPrint ("gainTablNext: bad chaining (2)");
        return     (NULL);
      }
#endif /* SCOTCH_DEBUG_GAIN3 */
      return (entrptr->next);
    }
  }

  return (NULL);
}

/* This routine checks the consistency of the
** given linked list.
** It returns:
** - !NULL  : pointer to the vertex.
** - NULL   : if no such vertex available.
*/

#ifdef SCOTCH_DEBUG_GAIN3

static
int
gainTablCheck2 (
GainEntr * const            entrptr,
GainLink * const            linkptr)
{
  GainLink *          nextptr;
  int                 o;

  if (linkptr->next == &gainLinkDummy)            /* If end of list successfully reached */
    return (0);
  if (linkptr->next == NULL)                      /* If loop discovered */
    return (1);

  if (linkptr->tabl != entrptr)
    return (1);

  nextptr       = linkptr->next;                  /* Save pointer to next cell    */
  linkptr->next = NULL;                           /* Flag cell as already visited */
  o             = gainTablCheck2 (entrptr, nextptr); /* Process next cell         */
  linkptr->next = nextptr;                        /* Restore cell state           */

  return (o);
}

int
gainTablCheck (
GainEntr * const            entrptr)
{
  if (entrptr->next == &gainLinkDummy)            /* If end of list successfully reached */
    return (0);

  if ((entrptr->next == NULL) ||                  /* If problem with link */
      (gainTablCheck2 (entrptr, entrptr->next) != 0)) {
    errorPrint ("gainTablCheck: bad linked list");
    return     (1);
  }

  return (0);
}

#endif /* SCOTCH_DEBUG_GAIN3 */
