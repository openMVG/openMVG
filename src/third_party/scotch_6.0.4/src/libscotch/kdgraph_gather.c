/* Copyright 2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : kdgraph_gather.c                        **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the routine which  **/
/**                builds a centralized mapping graph by   **/
/**                gathering the pieces of a distributed   **/
/**                mapping graph.                          **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 17 jun 2008     **/
/**                                 to   : 17 jun 2008     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define KDGRAPH_GATHER

#include "module.h"
#include "common.h"
#include "arch.h"
#include "graph.h"
#include "mapping.h"
#include "kgraph.h"
#include "dgraph.h"
#include "dmapping.h"
#include "kdgraph.h"

/******************************/
/*                            */
/* These routines handle      */
/* distributed source graphs. */
/*                            */
/******************************/

/* This function gathers the pieces of
** a distributed mapping graph to build a
** centralized mapping graph.
** There is no gathered vnumtab array if
** the original graph did not have one, as
** vertices are gathered in global order, or
** else the original vnumloctab is gathered.
** It returns:
** - 0   : if graph data are consistent.
** - !0  : on error.
*/

int
kdgraphGather (
Kdgraph * restrict const    dgrfptr,              /* Distributed halo graph */
Kgraph * restrict const     cgrfptr)              /* Centralized halo graph */
{
  return (dgraphGather (&dgrfptr->s, &cgrfptr->s));
}
