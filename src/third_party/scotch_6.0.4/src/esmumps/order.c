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
/**   NAME       : order.c                                 **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module computes orderings.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 20 aug 1998     **/
/**                                 to     24 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ORDER

#include "common.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/* This routine initializes the given
** ordering structure.
** It returns:
** - 0  : in all cases.
*/

int
orderInit (
Order * const               ordeptr)
{
  memSet (ordeptr, 0, sizeof (Order));

  return (0);
}

/* This routine frees the contents
** of the given ordering.
** It returns:
** - VOID  : in all cases.
*/

void
orderExit (
Order * const               ordeptr)
{
  if (ordeptr->rangtab != NULL)
    memFree (ordeptr->rangtab);
  if (ordeptr->permtab != NULL)
    memFree (ordeptr->permtab);
  if (ordeptr->peritab != NULL)
    memFree (ordeptr->peritab);

#ifdef ORDER_DEBUG
  memSet (ordeptr, ~0, sizeof (Order));
#endif /* ORDER_DEBUG */
}
