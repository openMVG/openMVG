/* Copyright 2007,2008,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph.c                               **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                bipartitioning graph data structure     **/
/**                handling routines.                      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 10 sep 2007     **/
/**                                 to     14 apr 2011     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to     31 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH

#include "module.h"
#include "common.h"
#include "arch.h"
#include "dgraph.h"
#include "dmapping.h"
#include "bdgraph.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* bipartition graphs.               */
/*                                   */
/*************************************/

/* This routine builds the active graph
** corresponding to the given bipartitioning
** job parameters.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
bdgraphInit (
Bdgraph * restrict const        actgrafptr,       /* Active graph                     */
const Dgraph * restrict const   indgrafptr,       /* Induced source subdgraph         */
const Dgraph * restrict const   srcgrafptr,       /* Original source graph            */
const Arch * restrict const     archptr,          /* Current mapping of halo vertices */
const ArchDom                   domnsubtab[])     /* Subdomains                       */
{
  Anum                domndist;                   /* Distance between both subdomains   */
  Anum                domnwght0;                  /* Processor workforce in each domain */
  Anum                domnwght1;

  domndist  = archDomDist (archptr, &domnsubtab[0], &domnsubtab[1]); /* Get distance between subdomains */
  domnwght0 = archDomWght (archptr, &domnsubtab[0]); /* Get weights of subdomains                       */
  domnwght1 = archDomWght (archptr, &domnsubtab[1]);
  actgrafptr->s            = *indgrafptr;            /* Get source graph data                        */
  actgrafptr->s.flagval   &= ~DGRAPHFREEALL;         /* Do not free contents of separation graph     */
  actgrafptr->s.vlblloctax = NULL;                   /* Never mind about vertex labels in the future */
  actgrafptr->veexloctax   = NULL;                   /* No external gain (yet)                       */
  actgrafptr->veexglbsum   = 0;
  actgrafptr->partgsttax   = NULL;                   /* Do not allocate frontier arrays yet */
  actgrafptr->fronloctab   = NULL;
 
  bdgraphInit2 (actgrafptr, domndist, domnwght0, domnwght1);

/* TODO: Compute external gains */
  
#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (actgrafptr) != 0) {
    errorPrint ("bdgraphInit: inconsistent graph data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  return (0);
}

void
bdgraphInit2 (
Bdgraph * restrict const        actgrafptr,       /* Active graph                       */
const Anum                      domndist,         /* Distance between both subdomains   */
const Anum                      domnwght0,        /* Processor workforce in each domain */
const Anum                      domnwght1)
{
  actgrafptr->fronlocnbr       =                  /* No frontier vertices */
  actgrafptr->fronglbnbr       = 0;
  actgrafptr->complocload0     = actgrafptr->s.velolocsum;
  actgrafptr->compglbload0     = actgrafptr->s.veloglbsum;
  actgrafptr->compglbload0min  = 0;               /* No external constraints on bipartition (yet) */
  actgrafptr->compglbload0max  = actgrafptr->s.veloglbsum;
  actgrafptr->compglbload0avg  = (Gnum) (((double) actgrafptr->s.veloglbsum * (double) domnwght0) / (double) (domnwght0 + domnwght1));
  actgrafptr->compglbload0dlt  = actgrafptr->s.veloglbsum - actgrafptr->compglbload0avg;
  actgrafptr->complocsize0     = actgrafptr->s.vertlocnbr;
  actgrafptr->compglbsize0     = actgrafptr->s.vertglbnbr;
  actgrafptr->commglbload      = 0;
  actgrafptr->commglbloadextn0 = 0;
  actgrafptr->commglbgainextn  = 0;
  actgrafptr->commglbgainextn0 = 0; 
  actgrafptr->bbalglbval       = (double) actgrafptr->compglbload0dlt / (double) actgrafptr->compglbload0avg;
  actgrafptr->domndist         = domndist;
  actgrafptr->domnwght[0]      = domnwght0;
  actgrafptr->domnwght[1]      = domnwght1;
  actgrafptr->levlnum          = 0;
}

/* This routine frees the contents
** of the given distributed active graph.
** It returns:
** - VOID  : in all cases.
*/

void
bdgraphExit (
Bdgraph * const             grafptr)
{
  if (grafptr->partgsttax != NULL)
    memFree (grafptr->partgsttax + grafptr->s.baseval);
  if (grafptr->fronloctab != NULL)
    memFree (grafptr->fronloctab);
  if (grafptr->veexloctax != NULL)
    memFree (grafptr->veexloctax + grafptr->s.baseval);

  dgraphExit (&grafptr->s);                       /* Free distributed source graph and its private data (flagval may be corrupted afterwards) */

#ifdef SCOTCH_DEBUG_BDGRAPH2
  memSet (grafptr, ~0, sizeof (Bdgraph));
#endif /* SCOTCH_DEBUG_BDGRAPH2 */
}

/* This routine moves all of the graph
** vertices to the first part.
** It returns:
** - VOID  : in all cases.
*/

void
bdgraphZero (
Bdgraph * const             grafptr)
{
  if (grafptr->partgsttax != NULL)
    memSet (grafptr->partgsttax + grafptr->s.baseval, 0, grafptr->s.vertgstnbr * sizeof (GraphPart)); /* Set all local and ghost vertices to part 0 */

  grafptr->fronlocnbr      =                      /* No frontier vertices */
  grafptr->fronglbnbr      = 0;
  grafptr->complocload0    = grafptr->s.velolocsum;
  grafptr->compglbload0    = grafptr->s.veloglbsum;
  grafptr->compglbload0dlt = grafptr->s.veloglbsum - grafptr->compglbload0avg;
  grafptr->complocsize0    = grafptr->s.vertlocnbr;
  grafptr->compglbsize0    = grafptr->s.vertglbnbr;
  grafptr->commglbload     = grafptr->commglbloadextn0;
  grafptr->commglbgainextn = grafptr->commglbgainextn0;
}
