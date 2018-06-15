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
/**   NAME       : vgraph_separate_es.h                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the edge-separation-based node      **/
/**                separation module.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 24 oct 1996     **/
/**                                 to   : 07 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 18 aug 2004     **/
/**                                 to     19 aug 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Separator type. +*/

typedef enum VgraphSeparateEsWidth_ {
  VGRAPHSEPAESWIDTHTHIN,                          /*+ Thin vertex separator +*/
  VGRAPHSEPAESWIDTHFAT                            /*+ Fat vertex separator  +*/
} VgraphSeparateEsWidth;

/*+ This structure holds the method parameters. +*/

typedef struct VgraphSeparateEsParam_ {
  Strat *                   strat;                /*+ Edge bipartitioning strategy used +*/
  VgraphSeparateEsWidth     widtval;              /*+ Separator width                   +*/
} VgraphSeparateEsParam;

/*+ These are the type of subgraphs vertices
    belong to, used to represent the Dulmage-
    Mendelsohn decomposition.
    TRICK: item numbers have been carefully
    chosen, so that one can easily sort out
    vertices that belong to (HR u SC u VC)
    and to (HR u SR u VC).                    +*/

typedef enum VgraphSeparateEsType_ {
  VGRAPHSEPAESTYPEHC     = 0x0000,
  VGRAPHSEPAESTYPEVR     = 0x0001,
  VGRAPHSEPAESTYPEHRSCVC = 0x0002,                /* Bit mask for testing */
  VGRAPHSEPAESTYPESC     = 0x0003,
  VGRAPHSEPAESTYPEHRSRVC = 0x0004,                /* Bit mask for testing */
  VGRAPHSEPAESTYPESR     = 0x0005,
  VGRAPHSEPAESTYPEHR     = 0x0006,
  VGRAPHSEPAESTYPEVC     = 0x0007
} VgraphSeparateEsType;

#define VGRAPHSEPAESTYPEBITC        1             /* Bit index for VGRAPHSEPAESTYPEHRSCVC */
#define VGRAPHSEPAESTYPEBITR        2             /* Bit index for VGRAPHSEPAESTYPEHRSRVC */

/*+ Vertex traversal flag. +*/

typedef enum VgraphSeparateEsTrav_ {
  VGRAPHSEPAESTRAVFREE = 0,                       /*+ Vertex not traversed                                     +*/
  VGRAPHSEPAESTRAVUSED,                           /*+ Vertex traversed by search for free rows                 +*/
  VGRAPHSEPAESTRAVDRTY                            /*+ Vertex traversed by backtracking search for free columns +*/
} VgraphSeparateEsTrav;

/*
**  The function prototypes.
*/

#ifndef VGRAPH_SEPARATE_ES
#define static
#endif

static int                  vgraphSeparateEsCover (const Graph * const, const Gnum, Gnum * const, Gnum * const);
static int                  vgraphSeparateEsCoverAugment (const Gnum * restrict const, const Gnum, Gnum * restrict const, VgraphSeparateEsTrav * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum);
static void                 vgraphSeparateEsCoverCol (const Gnum * restrict const, VgraphSeparateEsType * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum);
static void                 vgraphSeparateEsCoverRow (const Gnum * restrict const, VgraphSeparateEsType * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum * restrict const, const Gnum);

int                         vgraphSeparateEs    (Vgraph * const, const VgraphSeparateEsParam * const);

#undef static
