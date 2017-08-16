/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : graph_list.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                vertex lists functions.                 **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to     18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     31 oct 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to     08 jun 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     15 sep 1998     **/
/**                # Version 4.0  : from : 10 dec 2001     **/
/**                                 to     10 dec 2001     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_LIST

#include "module.h"
#include "common.h"
#include "graph.h"

/********************************/
/*                              */
/* These routines handle vertex */
/* number lists.                */
/*                              */
/********************************/

/* This routine creates a vertex number list.
** It returns:
** - 0  : in all cases.
*/

int
listInit (
VertList *          listptr)
{
  listptr->vnumnbr = 0;                           /* Initialize list fields */
  listptr->vnumtab = NULL;

  return (0);
}

/* This routine deletes the given vertex number list.
** It returns:
** - VOID  : in all cases.
*/

void
listExit (
VertList *          listptr)
{
  if (listptr->vnumtab != NULL)
    memFree (listptr->vnumtab);                   /* Free vertex list array */

#ifdef SCOTCH_DEBUG_GRAPH2
  memSet (listptr, 0, sizeof (VertList));         /* Purge list fields */
#endif /* SCOTCH_DEBUG_GRAPH2 */
}

/* This routine allocates a vertex
** number list array.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
listAlloc (
VertList *          listptr,
Gnum                vnumnbr)
{
  if (vnumnbr == listptr->vnumnbr)                /* If array is already dimensioned */
    return (0);                                   /* Keep it as it is                */

  listFree (listptr);                             /* Free vertex array       */
  if (vnumnbr > 0) {                              /* Reallocate vertex space */
    if ((listptr->vnumtab = (Gnum *) memAlloc (vnumnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("listAlloc: out of memory");
      return     (1);
    }
    listptr->vnumnbr = vnumnbr;
  }

  return (0);
}

/* This routine frees a vertex
** number list array.
** It returns:
** - 0  : in all cases.
*/

int
listFree (
VertList *          listptr)
{
  if (listptr->vnumtab != NULL)                   /* Free vertex list array */
    memFree (listptr->vnumtab);

  listptr->vnumnbr = 0;                           /* Reset list values */
  listptr->vnumtab = NULL;

  return (0);
}

/* These routines load a vertex number list
** from the given stream. Because of the search
** for duplicates, the list read is always
** sorted by ascending order.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
listLoad (
VertList *          listptr,
FILE *              stream)
{
  Gnum                vnumnbr;
  Gnum                vnumnum;

  if (intLoad (stream, &vnumnbr) != 1) {          /* Read number of vertices */
    errorPrint ("listLoad: bad input (1)");
    return     (1);
  }

  if (listAlloc (listptr, vnumnbr) != 0) {        /* Allocate vertex space */
    errorPrint ("listLoad: out of memory");
    return     (1);
  }

  for (vnumnum = 0; vnumnum < vnumnbr; vnumnum ++) { /* Read vertex list contents */
    if (intLoad (stream, &listptr->vnumtab[vnumnum]) != 1) {
      errorPrint ("listLoad: bad input (2)");
      return     (1);
    }
  }

  listSort (listptr);                             /* Sort vertex list by ascending order */
  for (vnumnum = 1; vnumnum < vnumnbr; vnumnum ++) { /* Search list for duplicates       */
    if (listptr->vnumtab[vnumnum] == listptr->vnumtab[vnumnum - 1]) {
      errorPrint ("listLoad: duplicate vertex numbers");
      return     (1);
    }
  }

  return (0);
}

/* This routine sorts a vertex list
** by ascending order.
** It returns:
** - VOID  : in all cases.
*/

void
listSort (
VertList *          listptr)
{
  intSort1asc1 (listptr->vnumtab, listptr->vnumnbr);
}

/* This routine saves a vertex number list
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
listSave (
VertList *          listptr,
FILE *              stream)
{
  Gnum                vnumnum;
  int                 o;

  o = (intSave (stream, listptr->vnumnbr) == 0);  /* Write number of vertices */
  for (vnumnum = 0; (o == 0) && (vnumnum < listptr->vnumnbr); vnumnum ++) {
    o = (fprintf (stream, "%c" GNUMSTRING,
                  ((vnumnum % 8) == 0) ? '\n' : '\t',
                  (Gnum) listptr->vnumtab[vnumnum]) == EOF);
  }
  o |= (fprintf (stream, "\n") == EOF);

  if (o != 0)
    errorPrint ("listSave: bad output");

  return (o);
}

/* This routine copies the contents
** of a vertex list into another.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
listCopy (
VertList *          dstlistptr,                   /* Destination list */
VertList *          srclistptr)                   /* Source list      */
{
  if (listAlloc (dstlistptr, dstlistptr->vnumnbr) != 0) { /* Allocate vertex space */
    errorPrint ("listCopy: out of memory");
    return     (1);
  }
  memCpy (dstlistptr->vnumtab,                    /* Copy list data */
          srclistptr->vnumtab,
          srclistptr->vnumnbr * sizeof (Gnum));

  return (0);
}
