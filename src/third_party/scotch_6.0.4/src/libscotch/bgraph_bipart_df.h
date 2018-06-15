/* Copyright 2004,2007,2011,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bgraph_bipart_df.h                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the function       **/
/**                declarations for the diffusion scheme   **/
/**                bipartitioning method.                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 09 jan 2007     **/
/**                                 to     28 may 2007     **/
/**                # Version 5.1  : from : 29 oct 2007     **/
/**                                 to     28 mar 2011     **/
/**                # Version 6.0  : from : 08 nov 2011     **/
/**                                 to     20 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/* Small non-zero float value. */

#define BGRAPHBIPARTDFEPSILON       (1.0F / (float) (GNUMMAX))

/*+ Sign masking operator. +*/

#define BGRAPHBIPARTDFGNUMSGNMSK(i) (- (Gnum) (((Gunum) (i)) >> (sizeof (Gnum) * 8 - 1)))

/*
**  The type and structure definitions.
*/

/*+ Job selection policy types. +*/

typedef enum BgraphBipartDfType_ {
  BGRAPHBIPARTDFTYPEBAL = 0,                      /*+ Balance to average         +*/
  BGRAPHBIPARTDFTYPEKEEP                          /*+ Preserve current imbalance +*/
} BgraphBipartDfType;

/*+ Method parameters. +*/

typedef struct BgraphBipartDfParam_ {
  INT                       passnbr;              /*+ Number of passes to do   +*/
  BgraphBipartDfType        typeval;              /*+ Type of balance to reach +*/
} BgraphBipartDfParam;

/*+ The loop routine parameter
    structure. It contains the
    thread-independent data.   +*/

typedef struct BgraphBipartDfData_ {
  ThreadGroupHeader         thrddat;
  Bgraph *                  grafptr;              /*+ Graph to work on          +*/
  float *                   difntax;              /*+ New diffusion value array +*/
  float *                   difotax;              /*+ Old diffusion value array +*/
  INT                       passnbr;              /*+ Number of passes          +*/
  Gnum                      vanctab[2];           /*+ Anchor load arrays        +*/
#ifdef BGRAPHBIPARTDFTHREAD
  int                       abrtval;              /*+ Abort value               +*/
#endif /* BGRAPHBIPARTDFTHREAD */
} BgraphBipartDfData;

/*+ The thread-specific data block. +*/

typedef struct BgraphBipartDfThread_ {
  ThreadHeader              thrddat;              /*+ Thread management data                       +*/
  Gnum                      vertbas;              /*+ Minimum regular vertex index                 +*/
  Gnum                      vertnnd;              /*+ After-last regular vertex index              +*/
  Gnum                      fronnnd;              /*+ After-last frontier vertex index             +*/
  Gnum                      compload1;            /*+ State return values to aggregate             +*/
  Gnum                      compsize1;
  Gnum                      commloadextn;
  Gnum                      commloadintn;
  Gnum                      commgainextn;
  float                     vanctab[2];           /*+ Area for (reducing) contributions to anchors +*/
#ifdef BGRAPHBIPARTDFTHREAD
  Gnum                      veexsum;              /*+ Area for reducing sums of external gains     +*/
  Gnum                      veexsum1;
#endif /* BGRAPHBIPARTDFTHREAD */
} BgraphBipartDfThread;

/*
**  The function prototypes.
*/

#ifndef BGRAPH_BIPART_DF
#define static
#endif

int                         bgraphBipartDf      (Bgraph * restrict const, const BgraphBipartDfParam * const);

#undef static
