/* Copyright 2004,2007,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_order_si.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders halo graph vertices  **/
/**                using a simple method.                  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 01 nov 1996     **/
/**                                 to     21 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     02 oct 1998     **/
/**                # Version 4.0  : from : 19 dec 2001     **/
/**                                 to     11 dec 2002     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to   : 04 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HGRAPH_ORDER_SI

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_si.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the ordering.
** It returns:
** - 0   : if the ordering could be computed.
** - !0  : on error.
*/

int
hgraphOrderSi (
const Hgraph * restrict const grafptr,
Order * restrict const        ordeptr,
const Gnum                    ordenum,            /*+ Zero-based ordering number +*/
OrderCblk * restrict const    cblkptr)            /*+ Single column-block        +*/
{
  Gnum                vnohnnd;
  Gnum                vertnum;
  Gnum                vnumnum;

  Gnum * restrict const       peritab = ordeptr->peritab;
  const Gnum * restrict const vnumtax = grafptr->s.vnumtax;

  vnohnnd = grafptr->vnohnnd;
  if (vnumtax == NULL) {                          /* If graph is original graph */
    for (vertnum = grafptr->s.baseval, vnumnum = ordenum;
         vertnum < vnohnnd; vertnum ++, vnumnum ++)
      peritab[vnumnum] = vertnum;
  }
  else {                                          /* Graph is not original graph */
    for (vertnum = grafptr->s.baseval, vnumnum = ordenum;
         vertnum < vnohnnd; vertnum ++, vnumnum ++)
      peritab[vnumnum] = vnumtax[vertnum];
  }

  return (0);
}
