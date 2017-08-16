/* Copyright 2004,2007,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vmesh_check.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the consistency    **/
/**                checker for separation meshes.          **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 21 mar 2003     **/
/**                                 to     11 may 2004     **/
/**                # Version 6.0  : from : 02 jun 2014     **/
/**                                 to     02 jun 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "vmesh.h"

/*************************/
/*                       */
/* These routines handle */
/* separator meshes.     */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given separator mesh.
** It returns:
** - 0   : if mesh data are consistent.
** - !0  : on error.
*/

int
vmeshCheck (
const Vmesh * const         meshptr)
{
  Gnum                velmnum;                    /* Number of current element vertex  */
  Gnum                vnodnum;                    /* Number of current node vertex     */
  Gnum                fronnum;                    /* Number of current frontier vertex */
  int * restrict      frontax;                    /* Frontier flag array               */
  Gnum                ecmpsize[2];                /* Elements never in separator       */
  Gnum                ncmpsize[3];
  Gnum                ncmpload[3];
  int                 o;

  if ((meshptr->ecmpsize[0] + meshptr->ecmpsize[1]) > meshptr->m.velmnbr) {
    errorPrint ("vmeshCheck: invalid element balance");
    return     (1);
  }
  if (meshptr->ncmploaddlt != (meshptr->ncmpload[0] - meshptr->ncmpload[1])) {
    errorPrint ("vmeshCheck: invalid node balance");
    return     (1);
  }

  ecmpsize[0] =
  ecmpsize[1] = 0;
  for (velmnum = meshptr->m.velmbas; velmnum < meshptr->m.velmnnd; velmnum ++) {
    Gnum                edgecut[3];               /* Array of cut edges     */
    Gnum                partnum;                  /* Part of current vertex */
    Gnum                eelmnum;                  /* Number of current edge */

    partnum = meshptr->parttax[velmnum];
    if ((partnum < 0) || (partnum > 1)) {
      errorPrint ("vmeshCheck: invalid part array (1)");
      return     (1);
    }
    ecmpsize[partnum] ++;

    if ((partnum != 0) &&
        (meshptr->m.verttax[velmnum] == meshptr->m.vendtax[velmnum])) {
      errorPrint ("vmeshCheck: isolated element not in part 0");
      return     (1);
    }

    edgecut[0] =
    edgecut[1] =
    edgecut[2] = 0;
    for (eelmnum = meshptr->m.verttax[velmnum];
         eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++)
      edgecut[meshptr->parttax[meshptr->m.edgetax[eelmnum]]] ++;

    if (partnum == 2) {
      if ((edgecut[0] != 0) || (edgecut[1] != 0)) {
        errorPrint ("vmeshCheck: separator element not surrounded by separator nodes");
        return     (1);
      }
    }
    else {
      if (edgecut[1 - partnum] != 0) {
        errorPrint ("vmeshCheck: element should be in separator (%ld)", (long) velmnum);
        return     (1);
      }
    }
  }
  if ((meshptr->ecmpsize[0] != ecmpsize[0]) ||
      (meshptr->ecmpsize[1] != ecmpsize[1])) {
    errorPrint ("vmeshCheck: invalid element parameters");
    return     (1);
  }

  ncmpload[0] =
  ncmpload[1] =
  ncmpload[2] = 0;
  ncmpsize[0] =
  ncmpsize[1] =
  ncmpsize[2] = 0;
  for (vnodnum = meshptr->m.vnodbas; vnodnum < meshptr->m.vnodnnd; vnodnum ++) {
    Gnum                edgecut[3];               /* Array of cut edges     */
    Gnum                partnum;                  /* Part of current vertex */
    Gnum                enodnum;                  /* Number of current edge */

    partnum = meshptr->parttax[vnodnum];
    if ((partnum < 0) || (partnum > 2)) {
      errorPrint ("vmeshCheck: invalid part array (2)");
      return     (1);
    }

    ncmpsize[partnum] ++;
    ncmpload[partnum] += (meshptr->m.vnlotax == NULL) ? 1 : meshptr->m.vnlotax[vnodnum];

    edgecut[0] =
    edgecut[1] =
    edgecut[2] = 0;
    for (enodnum = meshptr->m.verttax[vnodnum];
         enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++)
      edgecut[meshptr->parttax[meshptr->m.edgetax[enodnum]]] ++;

#ifdef SCOTCH_DEBUG_VMESH3
    if (partnum == 2) {
      if ((edgecut[0] == 0) ||
          (edgecut[1] == 0))
        errorPrint ("vmeshCheck: no-use separator vertex%s (%ld)", /* Warning only */
                    ((meshptr->levlnum == 0) ? " at level 0" : ""),
                    (long) vnodnum);
    }
    else {
#else
    if (partnum != 2) {
#endif /* SCOTCH_DEBUG_VMESH3 */
      if (edgecut[1 - partnum] != 0) {
        errorPrint ("vmeshCheck: node should be in separator (%ld)", (long) vnodnum);
        return     (1);
      }
    }
  }
  if ((meshptr->ncmpload[0] != ncmpload[0]) ||
      (meshptr->ncmpload[1] != ncmpload[1]) ||
      (meshptr->ncmpload[2] != ncmpload[2]) ||
      (meshptr->ncmpsize[0] != ncmpsize[0]) ||
      (meshptr->ncmpsize[1] != ncmpsize[1]) ||
      (meshptr->fronnbr     != ncmpsize[2])) {
    errorPrint ("vmeshCheck: invalid node parameters");
    return     (1);
  }

  if ((meshptr->fronnbr < 0) ||
      (meshptr->fronnbr > meshptr->m.vnodnbr)) {
    errorPrint ("vmeshCheck: invalid number of frontier vertices");
    return     (1);
  }
  if ((frontax = memAlloc (meshptr->m.vnodnbr * sizeof (int))) == NULL) {
    errorPrint ("vmeshCheck: out of memory");
    return     (1);
  }
  memSet (frontax, 0, meshptr->m.vnodnbr * sizeof (int));
  frontax -= meshptr->m.vnodbas;

  o = 1;                                          /* Assume failure when checking */
  for (fronnum = 0; fronnum < meshptr->fronnbr; fronnum ++) {
    Gnum                vnodnum;

    vnodnum = meshptr->frontab[fronnum];

    if ((vnodnum <  meshptr->m.vnodbas) ||
        (vnodnum >= meshptr->m.vnodnnd)) {
      errorPrint ("vmeshCheck: invalid vertex in frontier array");
      goto fail;
    }
    if (meshptr->parttax[vnodnum] != 2) {
      errorPrint ("vmeshCheck: invalid frontier array");
      goto fail;
    }
    if (frontax[vnodnum] != 0) {
      errorPrint ("vmeshCheck: duplicate node in frontier array");
      goto fail;
    }
    frontax[vnodnum] = 1;
  }

  o = 0;                                          /* Everything turned well */

fail :
  memFree (frontax + meshptr->m.vnodbas);

  return (o);
}
