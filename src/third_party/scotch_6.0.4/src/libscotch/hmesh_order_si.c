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
/**   NAME       : hmesh_order_si.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders halo mesh vertices   **/
/**                using the natural order.                **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 01 jan 2002     **/
/**                                 to     27 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_ORDER_SI

#include "module.h"
#include "common.h"
#include "graph.h"
#include "order.h"
#include "mesh.h"
#include "hmesh.h"
#include "hmesh_order_si.h"

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
hmeshOrderSi (
const Hmesh * restrict const  meshptr,
Order * restrict const        ordeptr,
const Gnum                    ordenum,
OrderCblk * restrict const    cblkptr)            /*+ Single column-block +*/
{
  Gnum                vnodnum;
  Gnum                ordeval;

  if (meshptr->m.vnumtax == NULL) {               /* If mesh is original mesh (no halo) */
#ifdef SCOTCH_DEBUG_ORDER2
    if (meshptr->m.vnodnbr != ordeptr->vnodnbr) {
      errorPrint ("hmeshOrderSi: invalid permutation bounds");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_ORDER2 */
    for (vnodnum = ordeptr->baseval, ordeval = ordenum;
         vnodnum < ordeptr->baseval + ordeptr->vnodnbr; vnodnum ++, ordeval ++) {
      ordeptr->peritab[ordeval] = vnodnum;
    }
  }
  else {                                          /* Mesh is not original mesh */
    for (vnodnum = meshptr->m.vnodbas, ordeval = ordenum;
         vnodnum < meshptr->vnohnnd; vnodnum ++, ordeval ++) {
      ordeptr->peritab[ordeval] = meshptr->m.vnumtax[vnodnum];
#ifdef SCOTCH_DEBUG_ORDER2
      if ((ordeptr->peritab[ordeval] <   ordeptr->baseval) ||
          (ordeptr->peritab[ordeval] >= (ordeptr->baseval + ordeptr->vnodnbr))) {
        errorPrint ("hmeshOrderSi: invalid permutation index");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_ORDER2 */
    }
  }

  return (0);
}
