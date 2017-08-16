/* Copyright 2010,2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_check.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the mapping graph  **/
/**                consistency checking routine.           **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 13 jul 2010     **/
/**                                 to     13 jul 2010     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     14 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "kgraph.h"

/*************************/
/*                       */
/* These routines handle */
/* mapping graphs.       */
/*                       */
/*************************/

/* This routine checks the consistency
** of the given mapping graph.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
kgraphCheck (
const Kgraph * restrict const grafptr)
{
  int * restrict      flagtax;                    /* Frontier flag array       */
  Gnum                vertnum;                    /* Number of current vertex  */
  Gnum                fronnum;                    /* Number of frontier vertex */
  Gnum                vfixnbr;                    /* Number of fixed vertices  */
  Gnum * restrict     comploadtab;
  Gnum                commload;
  Gnum                edloval;
  Anum                domnnum;
  int                 o;

  const Gnum                      baseval = grafptr->s.baseval;
  const Gnum                      vertnnd = grafptr->s.vertnnd;
  const Gnum * restrict const     verttax = grafptr->s.verttax;
  const Gnum * restrict const     vendtax = grafptr->s.vendtax;
  const Gnum * restrict const     velotax = grafptr->s.velotax;
  const Gnum * restrict const     edgetax = grafptr->s.edgetax;
  const Gnum * restrict const     edlotax = grafptr->s.edlotax;
  const Anum * restrict const     parttax = grafptr->m.parttax;
  const ArchDom * restrict const  domntab = grafptr->m.domntab;
  const Arch * const              archptr = grafptr->m.archptr;
  const Anum * restrict const     pfixtax = grafptr->pfixtax;
  const Gnum * restrict const     frontab = grafptr->frontab;

  if (&grafptr->s != grafptr->m.grafptr) {
    errorPrint ("kgraphCheck: invalid mapping graph");
    return     (1);
  }

  if ((grafptr->m.domnmax <= 0)                 ||
      (grafptr->m.domnnbr > grafptr->m.domnmax) ||
      (grafptr->m.domnnbr > grafptr->s.vertnbr)) {
    errorPrint ("kgraphCheck: invalid number of domains");
    return     (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &comploadtab, (size_t) (grafptr->m.domnnbr * sizeof (Gnum)),
                     &flagtax,     (size_t) (grafptr->s.vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("kgraphCheck: out of memory");
    return     (1);
  }
  memSet (comploadtab, 0, grafptr->m.domnnbr * sizeof (Gnum));
  memSet (flagtax,    ~0, grafptr->s.vertnbr * sizeof (Gnum));
  flagtax -= baseval;

  o = 1;                                          /* Assume failure when checking */
  for (vertnum = baseval, vfixnbr = 0;
       vertnum < vertnnd; vertnum ++) {
    Anum                partval;

    partval = parttax[vertnum];
    if ((partval < 0) ||
        (partval >= grafptr->m.domnnbr)) {
      errorPrint ("kgraphCheck: invalid part array");
      goto fail;
    }
    if (pfixtax != NULL) {
      Anum                pfixval;

      pfixval = pfixtax[vertnum];
      if (pfixval != ~0) {
        if (pfixval < 0) {
          errorPrint ("kgraphCheck: invalid fixed part value");
          goto fail;
        }
        if (pfixval != archDomNum (archptr, &grafptr->m.domntab[partval])) {
          errorPrint ("kgraphCheck: part index does not match fixed array");
          goto fail;
        }
        vfixnbr ++;
      }
    }
  }
  if (vfixnbr != grafptr->vfixnbr) {
    errorPrint ("kgraphCheck: invalid number of fixed vertices");
    goto fail;
  }

  if ((grafptr->fronnbr < 0) ||
      (grafptr->fronnbr > grafptr->s.vertnbr)) {
    errorPrint ("kgraphCheck: invalid number of frontier vertices");
    goto fail;
  }
  for (fronnum = 0; fronnum < grafptr->fronnbr; fronnum ++) {
    Gnum                vertnum;
    Gnum                edgenum;
    Anum                partval;
    Anum                flagval;

    vertnum = frontab[fronnum];
    if ((vertnum < baseval) || (vertnum >= vertnnd)) {
      errorPrint ("kgraphCheck: invalid vertex index in frontier array");
      goto fail;
    }
    if (flagtax[vertnum] != ~0) {
      errorPrint ("kgraphCheck: duplicate vertex in frontier array");
      goto fail;
    }
    flagtax[vertnum] = 0;
    partval = parttax[vertnum];

    for (edgenum = verttax[vertnum], flagval = 0;
         edgenum < vendtax[vertnum]; edgenum ++)
      flagval |= parttax[edgetax[edgenum]] ^ partval; /* Flag set if neighbor part differs from vertex part */

    if (flagval == 0) {
      errorPrint ("kgraphCheck: invalid vertex in frontier array");
      goto fail;
    }
  }

  commload = 0;
  edloval  = 1;                                   /* Assume edges are not weighted */
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    Anum                partval;                  /* Part of current vertex                         */
    Anum                tdomnum;                  /* Terminal domain number of current vertex       */
    Anum                tdfinum;                  /* Terminal domain number of current fixed vertex */
    Gnum                edgenum;                  /* Number of current edge                         */
    Gnum                commcut;

    partval = parttax[vertnum];
    tdomnum = archDomNum (&grafptr->a, &grafptr->m.domntab[partval]);
    tdfinum = (grafptr->pfixtax != NULL) ? grafptr->pfixtax[vertnum] : -1;

    if ((tdfinum != -1) &&                        /* If it is a fixed vertex */
        (tdfinum != tdomnum)) {
      errorPrint ("kgraphCheck: invalid fixed vertex part");
      goto fail;
    }

    comploadtab[partval] += (velotax == NULL) ? 1 : velotax[vertnum];

    commcut = 0;
    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
      Anum                partend;

      if (edlotax != NULL)
        edloval = edlotax[edgenum];
      partend = parttax[edgetax[edgenum]];

      if (partend == partval)                     /* If same part, no communication load */
        continue;

      commcut += edloval * archDomDist (archptr, &domntab[partval], &domntab[partend]); /* Loads are accounted for twice */
    }

    if ((commcut != 0) && (flagtax[vertnum] != 0)) { /* If vertex should be in frontier array */
      errorPrint ("kgraphCheck: vertex should be in frontier array");
      goto fail;
    }

    commload += commcut;
  }
  commload /= 2;
  if (commload != grafptr->commload) {
    errorPrint ("kgraphCheck: invalid communication load");
    goto fail;
  }

  for (domnnum = 0; domnnum < grafptr->m.domnnbr; domnnum ++) {
    if (comploadtab[domnnum] != (grafptr->comploadavg[domnnum] + grafptr->comploaddlt[domnnum])) {
      errorPrint ("kgraphCheck: invalid computation load");
      goto fail;
    }
  }

  o = 0;                                          /* Everything turned well */

fail :
  memFree (comploadtab);                          /* Free group leader */

  return (o);
}
