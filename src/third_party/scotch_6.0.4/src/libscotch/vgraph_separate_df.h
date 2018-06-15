/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vgraph_separate_df.h                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the diffusion vertex separation     **/
/**                method.                                 **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 29 oct 2007     **/
/**                                 to     24 may 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/* Small non-zero float value. */

#define VGRAPHSEPARATEDFEPSILON     (1.0F / (float) (GNUMMAX))

/*
**  The type and structure definitions.
*/

/*+ Method parameters. +*/

typedef struct VgraphSeparateDfParam_ {
  INT                       partval;              /*+ Part to aggregate to separator +*/
  INT                       movenbr;              /*+ Number of moves to do          +*/
  INT                       passnbr;              /*+ Number of passes to do         +*/
  double                    cdifval;              /*+ Coefficient of diffused load   +*/
  double                    cremval;              /*+ Coefficient of remaining load  +*/
} VgraphSeparateDfParam;

/*
**  The function prototypes.
*/

#ifndef VGRAPH_SEPARATE_DF
#define static
#endif

int                         vgraphSeparateDf    (Vgraph * restrict const, const VgraphSeparateDfParam * restrict const);

#undef static
