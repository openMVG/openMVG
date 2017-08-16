/* Copyright 2009-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_ml.h                         **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module contains the function       **/
/**                declarations for the diffusion scheme   **/
/**                k-partitioning method.                  **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 22 dec 2009     **/
/**                                 to     02 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/* Small non-zero float value. */

#define KGRAPHMAPDFEPSILON       (1.0F / (float) (GNUMMAX))

/*+ Sign masking operator. +*/

#define KGRAPHMAPDFGNUMSGNMSK(i) ((Gnum) 0 - (((Gunum) (i)) >> (sizeof (Gnum) * 8 - 1)))

/*
**  The type and structure definitions.
*/

/*+ Method parameters. +*/

typedef struct KgraphMapDfParam_ {
  INT                       passnbr;              /*+ Number of passes to do        +*/
  double                    cdifval;              /*+ Coefficient of diffused load  +*/
  double                    cremval;              /*+ Coefficient of remaining load +*/
} KgraphMapDfParam;

/*+ The complementary vertex structure. +*/

typedef struct KgraphMapDfVertex_ {
  Anum                      partval;              /*+ Type of liquid in barrel                                       +*/
  float                     diffval;              /*+ Value to be diffused to everybody                              +*/
  float                     fdifval;              /*+ Value to be diffused to other parts for mapping                +*/
  float                     mdisval;              /*+ Value to be diffused to vertnum if parotax[vertnum] == partval +*/
  float                     mdidval;              /*= Value to be diffused to vertnum else                           +*/
} KgraphMapDfVertex;

/*+ The sort structure. +*/

typedef struct KgraphMapDfSort_ {
  Anum                      partval;              /*+ Type of liquid in barrel +*/
  float                     diffval;              /*+ Value to be diffused     +*/
  Anum                      distval;              /*+ Distance value           +*/
  Gnum                      edlosum;              /*+ Sum of edge loads        +*/
} KgraphMapDfSort;

/*+ The loop routine parameter
    structure. It contains the
    thread-independent data.   +*/

typedef struct KgraphMapDfData_ {
  ThreadGroupHeader         thrddat;
  const Kgraph *            grafptr;              /*+ Graph to work on           +*/
  float *                   vanctab;
  float *                   valotab;              /*+ Fraction of load to leak   +*/
  Gnum *                    velstax;              /*+ Vertex edge load sum array +*/
  KgraphMapDfVertex *       difntax;              /*+ New diffusion value array  +*/
  KgraphMapDfVertex *       difotax;              /*+ Old diffusion value array  +*/
  int                       passnbr;              /*+ Number of passes           +*/
#ifdef KGRAPHMAPDFTHREAD
  int                       abrtval;              /*+ Abort value                +*/
#endif /* KGRAPHMAPDFTHREAD */
} KgraphMapDfData;

/*+ The thread-specific data block. +*/

typedef struct KgraphMapDfThread_ {
  ThreadHeader              thrddat;              /*+ Thread management data          +*/
  Gnum                      vertbas;              /*+ Minimum regular vertex index    +*/
  Gnum                      vertnnd;              /*+ After-last regular vertex index +*/
  Anum                      domnbas;              /*+ Minimum anchor vertex index     +*/
  Anum                      domnnnd;              /*+ After-last anchor vertex index  +*/
} KgraphMapDfThread;

/*
**  The function prototypes.
*/

#ifndef KGRAPH_MAP_DF
#define static
#endif

int                         kgraphMapDf      (Kgraph * restrict const, const KgraphMapDfParam * const);

static void                 kgraphMapDfSort  (void * const, const INT);

#undef static
