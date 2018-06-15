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
/**   NAME       : hmesh_check.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the halo source     **/
/**                mesh functions.                         **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 12 sep 2002     **/
/**                                 to     11 may 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "hmesh.h"

/****************************************/
/*                                      */
/* These routines handle source meshes. */
/*                                      */
/****************************************/

/* This routine checks the consistency
** of the given halo mesh.
** It returns:
** - 0   : if halo mesh data are consistent.
** - !0  : on error.
*/

int
hmeshCheck (
const Hmesh * const         meshptr)
{
  Gnum                vnhlsum;

  if ((meshptr->vnohnnd < meshptr->m.vnodbas) ||
      (meshptr->vnohnnd > meshptr->m.vnodnnd)) {
    errorPrint ("hmeshCheck: invalid halo node numbers");
    return     (1);
  }
 
  if (meshCheck (&meshptr->m) != 0) {
    errorPrint ("hmeshCheck: invalid non-halo mesh structure");
    return     (1);
  }

  if (meshptr->vehdtax != meshptr->m.vendtax) {
    Gnum                veihnbr;
    Gnum                velmnum;

    for (velmnum = meshptr->m.velmbas, veihnbr = 0; /* For all element vertices */
         velmnum < meshptr->m.velmnnd; velmnum ++) {
      if ((meshptr->vehdtax[velmnum] < meshptr->m.verttax[velmnum]) ||
          (meshptr->vehdtax[velmnum] > meshptr->m.vendtax[velmnum])) {
        errorPrint ("hmeshCheck: invalid non-halo end vertex array");
        return     (1);
      }
      if (meshptr->vehdtax[velmnum] == meshptr->m.verttax[velmnum])
        veihnbr ++;
    }
    if (veihnbr != meshptr->veihnbr) {
      errorPrint ("hmeshCheck: invalid number of halo-isolated element vertices (1)");
      return     (1);
    }
  }
  else {
    if (meshptr->veihnbr != 0) {
      errorPrint ("hmeshCheck: invalid number of halo-isolated element vertices (2)");
      return     (1);
    }
  }

  if (meshptr->m.vnlotax == NULL)                 /* Recompute non-halo node vertex load sum */
    vnhlsum = meshptr->vnohnnd - meshptr->m.vnodbas;
  else {
    Gnum                vnodnum;

    for (vnodnum = meshptr->m.vnodbas, vnhlsum = 0;
         vnodnum < meshptr->vnohnnd; vnodnum ++)
      vnhlsum += meshptr->m.vnlotax[vnodnum];
  }
  if (vnhlsum != meshptr->vnhlsum) {
    errorPrint ("hmeshCheck: invalid non-halo vertex load sum");
    return     (1);
  }

  return (0);
}
