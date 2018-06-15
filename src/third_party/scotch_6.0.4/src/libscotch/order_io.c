/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : order_io.c                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles generic orderings.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 oct 1996     **/
/**                                 to     27 aug 1998     **/
/**                # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     28 jun 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     27 feb 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER_IO

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine loads an ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
orderLoad (
Order * restrict const      ordeptr,
const Gnum * restrict const vlbltab,
FILE * restrict const       stream)
{
  Gnum * restrict     permtab;
  Gnum                vertnum;

  if (vlbltab != NULL) {
    errorPrint ("orderLoad: vertex labels not yet supported");
    return     (1);
  }

  if ((permtab = memAlloc (ordeptr->vnodnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("orderLoad: out of memory");
    return     (1);
  }

  if (intLoad (stream, &ordeptr->vnodnbr) != 1) {
    errorPrint ("orderLoad: bad input (1)");
    memFree    (permtab);
    return     (1);
  }

  if (vlbltab == NULL) {                          /* If ordering does not have label array */
    for (vertnum = 0; vertnum < ordeptr->vnodnbr; vertnum ++) {
      Gnum                vertval;

      if ((intLoad (stream, &vertval)          != 1) || /* Read item data */
          (intLoad (stream, &permtab[vertnum]) != 1)) {
        errorPrint ("orderLoad: bad input (2)");
        memFree    (permtab);
        return     (1);
      }
      if (vertval != (vertnum + ordeptr->baseval)) { /* Read item data */
        errorPrint ("orderLoad: bad input (3)");
        memFree    (permtab);
        return     (1);
      }
    }
  }

  orderPeri (permtab, ordeptr->baseval, ordeptr->vnodnbr, ordeptr->peritab, ordeptr->baseval); /* Compute inverse permutation */

  memFree (permtab);
  return  (0);
}

/* This routine saves an ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
orderSave (
const Order * restrict const  ordeptr,
const Gnum * restrict const   vlbltab,
FILE * restrict const         stream)
{
  const Gnum * restrict vlbltax;
  Gnum * restrict       permtab;
  Gnum                  vertnum;

  vlbltax = (vlbltab != NULL) ? (vlbltab - ordeptr->baseval) : NULL;

  if ((permtab = memAlloc (ordeptr->vnodnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("orderSave: out of memory");
    return     (1);
  }

  if (fprintf (stream, GNUMSTRING "\n",
               (Gnum) ordeptr->vnodnbr) == EOF) {
    errorPrint ("orderSave: bad output (1)");
    memFree    (permtab);
    return     (1);
  }

  orderPeri (ordeptr->peritab, ordeptr->baseval, ordeptr->vnodnbr, permtab, ordeptr->baseval); /* Compute direct permutation */

  if (vlbltax != NULL) {                          /* If ordering has label array */
    for (vertnum = 0; vertnum < ordeptr->vnodnbr; vertnum ++) {
      if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                   (Gnum) vlbltax[vertnum + ordeptr->baseval],
                   (Gnum) vlbltax[permtab[vertnum]]) == EOF) {
        errorPrint ("orderSave: bad output (2)");
        memFree    (permtab);
        return     (1);
      }
    }
  }
  else {
    for (vertnum = 0; vertnum < ordeptr->vnodnbr; vertnum ++) {
      if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                   (Gnum) (vertnum + ordeptr->baseval),
                   (Gnum) permtab[vertnum]) == EOF) {
        errorPrint ("orderSave: bad output (3)");
        memFree    (permtab);
        return     (1);
      }
    }
  }

  memFree (permtab);
  return  (0);
}

/* This routine saves the column block
** mapping data of the given ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
orderSaveMap (
const Order * restrict const  ordeptr,
const Gnum * restrict const   vlbltab,
FILE * const                  stream)
{
  const Gnum * restrict vlbltax;
  const Gnum * restrict peritax;
  Gnum * restrict       rangtab;
  Gnum * restrict       cblktax;
  Gnum                  vnodnnd;
  Gnum                  vnodnum;
  Gnum                  cblknum;
  int                   o;

  if (fprintf (stream, GNUMSTRING "\n",
               (Gnum) ordeptr->vnodnbr) == EOF) {
    errorPrint ("orderSaveMap: bad output (1)");
    return     (1);
  }

  if (memAllocGroup ((void **) (void *)
        &rangtab, (size_t) ((ordeptr->vnodnbr + 1) * sizeof (Gnum)),
        &cblktax, (size_t) ( ordeptr->vnodnbr      * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("orderSaveMap: out of memory");
    return     (1);
  }
  cblktax -= ordeptr->baseval;

  orderRang (ordeptr, rangtab);
  peritax = ordeptr->peritab - ordeptr->baseval;
  for (vnodnum = ordeptr->baseval, vnodnnd = vnodnum + ordeptr->vnodnbr, cblknum = 0;
       vnodnum < vnodnnd; vnodnum ++) {
    if (vnodnum >= rangtab[cblknum + 1])
      cblknum ++;
    cblktax[peritax[vnodnum]] = cblknum;
  }

  vlbltax = (vlbltab != NULL) ? (vlbltab - ordeptr->baseval) : NULL;
  for (vnodnum = ordeptr->baseval, o = 0; vnodnum < vnodnnd; vnodnum ++) {
    if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                 (Gnum) ((vlbltax != NULL) ? vlbltax[vnodnum] : vnodnum),
                 (Gnum) cblktax[vnodnum]) == EOF) {
      errorPrint ("orderSaveMap: bad output (2)");
      o = 1;
      break;
    }
  }

  memFree (rangtab);                              /* Free memory group leader */
  return  (o);
}

/* This routine saves the separator
** tree data of the given ordering.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
orderSaveTree (
const Order * restrict const  ordeptr,
const Gnum * restrict const   vlbltab,
FILE * const                  stream)
{
  const Gnum * restrict vlbltax;
  const Gnum * restrict peritax;
  Gnum * restrict       rangtab;
  Gnum * restrict       treetab;
  Gnum * restrict       cblktax;
  Gnum                  vnodnnd;
  Gnum                  vnodnum;
  Gnum                  cblknum;
  int                   o;

  if (fprintf (stream, GNUMSTRING "\n",
               (Gnum) ordeptr->vnodnbr) == EOF) {
    errorPrint ("orderSaveTree: bad output (1)");
    return     (1);
  }

  if (memAllocGroup ((void **) (void *)
        &rangtab, (size_t) ((ordeptr->vnodnbr + 1) * sizeof (Gnum)),
        &treetab, (size_t) ( ordeptr->vnodnbr      * sizeof (Gnum)),
        &cblktax, (size_t) ( ordeptr->vnodnbr      * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("orderSaveTree: out of memory");
    return     (1);
  }
  cblktax -= ordeptr->baseval;

  orderRang (ordeptr, rangtab);
  orderTree (ordeptr, treetab);
  peritax = ordeptr->peritab - ordeptr->baseval;
  for (vnodnum = ordeptr->baseval, vnodnnd = vnodnum + ordeptr->vnodnbr, cblknum = 0;
       vnodnum < vnodnnd; vnodnum ++) {
    if (vnodnum >= rangtab[cblknum + 1])
      cblknum ++;
    cblktax[peritax[vnodnum]] = treetab[cblknum];
  }

  vlbltax = (vlbltab != NULL) ? (vlbltab - ordeptr->baseval) : NULL;
  for (vnodnum = ordeptr->baseval, o = 0; vnodnum < vnodnnd; vnodnum ++) {
    if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                 (Gnum) ((vlbltax != NULL) ? vlbltax[vnodnum] : vnodnum),
                 (Gnum) cblktax[vnodnum]) == EOF) {
      errorPrint ("orderSaveMap: bad output (2)");
      o = 1;
      break;
    }
  }

  memFree (rangtab);                              /* Free memory group leader */
  return  (o);
}
