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
/**   NAME       : vmesh_separate_gg.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the the greedy mesh growing node    **/
/**                separation method.                      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 16 sep 2002     **/
/**                                 to     07 apr 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ System-defined constants. +*/

#define VMESHSEPAGGSUBBITS          4

#define VMESHSEPAGGSTATEPART0       ((GainLink *) 0) /*+ Element vertex in part 0 (initial state)  +*/
#define VMESHSEPAGGSTATEPART1       ((GainLink *) 1) /*+ Element vertex in part 1                  +*/
#define VMESHSEPAGGSTATEPART2       ((GainLink *) 2) /*+ Element vertex in part 2, chained         +*/
#define VMESHSEPAGGSTATELINK        ((GainLink *) 3) /*+ Currently in gain table if higher         +*/

/*
**  The type and structure definitions.
*/

/*+ Method parameters. +*/

typedef struct VmeshSeparateGgParam_ {
  INT                       passnbr;              /*+ Number of passes to perform +*/
} VmeshSeparateGgParam;

/*+ The complementary element vertex structure.
    For trick reasons, the gain table data structure
    must be the first field of the structure.        +*/

typedef struct VmeshSeparateGgElem_ {
  GainLink                  gainlink;             /*+ Gain link: FIRST                               +*/
  Gnum                      ncmpgain2;            /*+ Computation gain in separator: (0->2) - (2->1) +*/
  Gnum                      ncmpgaindlt;          /*+ Overall computation delta: - (0->2) - (2->1)   +*/
} VmeshSeparateGgElem;

/*+ The complementary vertex structure. Only
    partval is always valid. Other fields are
    valid only when vertex belongs to separator. +*/

typedef struct VmeshSeparateGgNode_ {
  int                       partval;              /*+ Part to which node vertex belongs                                   +*/
  Gnum                      commsize0;            /*+ Number of neighbors in part 0                                       +*/
  Gnum                      velmisum0;            /*+ Sum of all element indices in part 0; the last one is the right one +*/
} VmeshSeparateGgNode;

/*
**  The function prototypes.
*/

#ifndef VMESH_SEPARATE_GG
#define static
#endif

int                         vmeshSeparateGg     (Vmesh * restrict const, const VmeshSeparateGgParam * restrict const);

#ifdef SCOTCH_DEBUG_VMESH3
static int                  vmeshSeparateGgCheck (Vmesh * restrict const, const Gnum, const Gnum, const VmeshSeparateGgElem * restrict const, const VmeshSeparateGgNode * restrict const  vnoxtax);
#endif /* SCOTCH_DEBUG_VMESH3 */

#undef static
