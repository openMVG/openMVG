/* Copyright 2007-2009,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_match.h                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : This file provides the structure        **/
/**                definitions for the coasening phase of  **/
/**                the multi-level method.                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 04 may 2006     **/
/**                                 to   : 21 jun 2007     **/
/**                # Version 5.1  : from : 23 nov 2008     **/
/**                                 to   : 04 apr 2009     **/
/**                # Version 6.0  : from : 03 oct 2012     **/
/**                                 to   : 03 oct 2012     **/
/**                                                        **/
/************************************************************/

/*
** The type and structure definitions.
*/

/* This structure contains all the tables describing the matching state */

typedef struct DgraphMatchData_ {
  DgraphCoarsenData         c;                    /*+ Coarsening structure                                        +*/
  Gnum *                    mategsttax;           /*+ Mating table for local and ghost vertices                   +*/
  Gnum                      matelocnbr;           /*+ Number of local matchings                                   +*/
  Gnum *                    queuloctab;           /*+ Queue of unmated local vertex                               +*/
  Gnum                      queulocnbr;           /*+ Number of enqueued unmated vertices                         +*/
  Gnum *                    procvgbtab;           /*+ Global vertex number bounds for neighboring processors [+1] +*/
  float                     probval;              /*+ Vertex mating probability (1.0 is certain)                  +*/
} DgraphMatchData;

/*
** The function prototypes.
*/

int                         dgraphMatchInit     (DgraphMatchData * restrict const, const float);
void                        dgraphMatchExit     (DgraphMatchData * restrict const);
int                         dgraphMatchSync     (DgraphMatchData * restrict const);
int                         dgraphMatchSyncColl (DgraphMatchData * restrict const);
int                         dgraphMatchSyncPtop (DgraphMatchData * restrict const);
int                         dgraphMatchCheck    (DgraphMatchData * restrict const);

void                        dgraphMatchHl       (DgraphMatchData * restrict const);
void                        dgraphMatchSc       (DgraphMatchData * restrict const);
void                        dgraphMatchHy       (DgraphMatchData * restrict const);
void                        dgraphMatchLc       (DgraphMatchData * restrict const);
void                        dgraphMatchLy       (DgraphMatchData * restrict const);
