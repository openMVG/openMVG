/* Copyright 2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_cp.c                         **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This method copies a given old mapping  **/
/**                as a mapping result.                    **/ 
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 16 jan 2012     **/
/**                                 to     23 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define KGRAPH_MAP_CP

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"
#include "parser.h"
#include "kgraph.h"
#include "kgraph_map_cp.h"
#include "kgraph_map_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the k-way partitioning.
** It returns:
** - 0 : if k-partition could be computed.
** - 1 : on error.
*/

/* TODO handle the case when the old and the new architectures are differents. */

int
kgraphMapCp (
Kgraph * restrict const     grafptr)              /*+ Graph +*/
{
  Gnum                          baseval;
  Anum                          domnnbr;

  const Anum * restrict const pfixtax = grafptr->pfixtax;

  if (grafptr->r.m.parttax == NULL) {             /* If we do not have an old partition */
    errorPrint ("kgraphMapCp: inconsistent old mapping data");
    return     (1);
  }
  baseval = grafptr->s.baseval;

  if (mapCopy (&grafptr->m, &grafptr->r.m) != 0) {
    errorPrint ("kgraphMapCp: cannot copy old mapping");
    return     (1);
  }

  if (pfixtax != NULL) {                          /* If we have fixed vertices */
    if (mapMerge (&grafptr->m, pfixtax) != 0) {
      errorPrint ("kgraphMapCp: cannot merge with fixed vertices");
      return     (1);
    }
  }

  kgraphFron (grafptr);
  kgraphCost (grafptr);

#ifdef SCOTCH_DEBUG_KGRAPH2
  if (kgraphCheck (grafptr) != 0) {
    errorPrint ("kgraphMapCp: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_KGRAPH2 */

  return (0);
}


