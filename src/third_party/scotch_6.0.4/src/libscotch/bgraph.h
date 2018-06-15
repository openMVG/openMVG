/* Copyright 2004,2007,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the job control routines of the     **/
/**                Dual Recursive Bipartitioning method.   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 23 mar 1993     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     09 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     04 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     28 sep 1995     **/
/**                # Version 3.1  : from : 15 nov 1995     **/
/**                                 to     15 nov 1995     **/
/**                # Version 3.2  : from : 24 aug 1996     **/
/**                                 to     03 nov 1997     **/
/**                # Version 3.3  : from : 01 dec 1998     **/
/**                                 to     02 dec 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     05 may 2006     **/
/**                # Version 5.0  : from : 30 nov 2006     **/
/**                                 to     30 nov 2006     **/
/**                # Version 5.1  : from : 08 jan 2008     **/
/**                                 to     18 mar 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to     08 aug 2014     **/
/**                                                        **/
/************************************************************/

#define BGRAPH_H

/*
**  The type and structure definitions.
*/

/*+ Graph option flags. +*/

#define BGRAPHFREEFRON              (GRAPHBITSNOTUSED) /* Free frontier array               */
#define BGRAPHFREEPART              (GRAPHBITSNOTUSED << 1) /* Free part array              */
#define BGRAPHFREEVEEX              (GRAPHBITSNOTUSED << 2) /* Free external gain array     */
#define BGRAPHHASANCHORS            (GRAPHBITSNOTUSED << 3) /* If graph has anchor vertices */

/*+ The bipartition graph structure. +*/

typedef struct Bgraph_ {
  Graph                     s;                    /*+ Source graph data                                  +*/
  Gnum *                    veextax;              /*+ Array of vertex external gain if moved to part 1   +*/
  GraphPart *               parttax;              /*+ Array of parts for every vertex                    +*/
  Gnum *                    frontab;              /*+ Array of frontier vertex numbers                   +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices                        +*/
  Gnum                      compload0min;         /*+ Minimum allowed load in part 0 (strategy variable) +*/
  Gnum                      compload0max;         /*+ Maximum allowed load in part 0 (strategy variable) +*/
  Gnum                      compload0avg;         /*+ Average load of part 0                             +*/
  Gnum                      compload0dlt;         /*+ Difference from the average                        +*/
  Gnum                      compload0;            /*+ Load in part 0 (strategy variable)                 +*/
  Gnum                      compsize0;            /*+ Number of vertices in part 0                       +*/
  Gnum                      commload;             /*+ Communication load                                 +*/
  Gnum                      commloadextn0;        /*+ Communication load if all moved to part 0          +*/
  Gnum                      commgainextn0;        /*+ External gain if all swapped from part 0           +*/
  Gnum                      commgainextn;         /*+ External gain if all swapped                       +*/
  double                    bbalval;              /*+ Bipartitioning imbalance ratio (strategy variable) +*/
  Anum                      domndist;             /*+ Distance between subdomains                        +*/
  Gnum                      domnwght[2];          /*+ Weights of the two subdomains                      +*/
  Gnum                      vfixload[2];          /*+ Vertex load biases of the two subdomains           +*/
  INT                       levlnum;              /*+ Coarsening level                                   +*/
} Bgraph;

/*+ The save graph structure. +*/

typedef struct BgraphStore_ {
  Gnum                      fronnbr;              /*+ Number of frontier nodes      +*/
  Gnum                      compload0dlt;         /*+ Difference from the average   +*/
  Gnum                      compsize0;            /*+ Number of vertices in part 0  +*/
  Gnum                      commload;             /*+ Communication load            +*/
  Gnum                      commgainextn;         /*+ External gain if all swapped  +*/
  byte *                    datatab;              /*+ Variable-sized data array     +*/
} BgraphStore;

/*
**  The function prototypes.
*/

#ifndef BGRAPH
#define static
#endif

int                         bgraphInit          (Bgraph * restrict const, const Graph * restrict const, const Arch * restrict const, const ArchDom * restrict const, const Gnum * restrict const);
void                        bgraphInit2         (Bgraph * restrict const, const Anum, const Anum, const Anum, const Gnum, const Gnum);
void                        bgraphExit          (Bgraph * restrict const);
void                        bgraphSwal          (Bgraph * restrict const);
void                        bgraphZero          (Bgraph * restrict const);
int                         bgraphCheck         (const Bgraph * restrict const);

int                         bgraphStoreInit     (const Bgraph * const, BgraphStore * const);
void                        bgraphStoreExit     (BgraphStore * const);
void                        bgraphStoreSave     (const Bgraph * const, BgraphStore * const);
void                        bgraphStoreUpdt     (Bgraph * const, const BgraphStore * const);

#undef static
