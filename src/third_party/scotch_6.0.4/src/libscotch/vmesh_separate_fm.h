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
/**   NAME       : vmesh_separate_fm.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the improved Fiduccia-Mattheyses    **/
/**                mesh element separation routine.        **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 26 feb 2003     **/
/**                                 to     06 may 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ Gain table subbits. +*/

#define VMESHSEPAFMGAINBITS         4

/*+ Prime number for hashing vertex numbers. +*/

#define VMESHSEPAFMHASHPRIME        11            /*+ Prime number for hashing +*/

/*+ Gain table vertex status. +*/

#define VMESHSEPAFMSTATEFREE        ((GainLink *) 0) /*+ Element is free or separator-chained +*/
#define VMESHSEPAFMSTATEUSED        ((GainLink *) 1) /*+ Element already swapped              +*/
#define VMESHSEPAFMSTATELINK        ((GainLink *) 2) /*+ Currently in gain table if higher    +*/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct VmeshSeparateFmParam_ {
  INT                       movenbr;              /*+ Maximum number of uneffective moves that can be done +*/
  INT                       passnbr;              /*+ Number of passes to be performed (-1 : infinite)     +*/
  double                    deltrat;              /*+ Maximum weight imbalance ratio                       +*/
} VmeshSeparateFmParam;

/*+ The hash element structure. The goal
    of both hash arrays is to record partition
    data that supercedes the one contained in
    the calling Vmesh structure, until the newly
    computed partition is written back to the
    Vmesh.                                       +*/

typedef struct VmeshSeparateFmElement_ {
  GainLink                  gainlink;             /*+ Gain link if moved to other part; FIRST +*/
  Gnum                      velmnum;              /*+ Number of vertex in hash table          +*/
  int                       vertpart;             /*+ Vertex part                             +*/
  Gnum                      ncmpcut2;             /*+ Number of neighbor nodes in separator   +*/
  Gnum                      ncmpgain2;            /*+ Separator gain if moved to given part   +*/
  Gnum                      ncmpgaindlt;          /*+ Node load imbalance if element swapped  +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded +*/
} VmeshSeparateFmElement;

/*+ The hash node structure. +*/

typedef struct VmeshSeparateFmNode_ {
  Gnum                      vnodnum;              /*+ Number of vertex in hash table          +*/
  Gnum                      vnloval;              /*+ Vertex load                             +*/
  int                       vertpart;             /*+ Vertex part                             +*/
  Gnum                      ecmpsize0;            /*+ Number of element neighbors in part 0   +*/
  Gnum                      mswpnum;              /*+ Number of move sweep when data recorded +*/
} VmeshSeparateFmNode;

/*+ The move recording structure. +*/

typedef struct VmeshSeparateFmSave_ {
  Gnum                      hertnum;              /*+ Hash index of vertex, (helmnum) or (-1 - hnodnum) +*/
  union {                                         /*+ Stored data to recover                            +*/
    struct {                                      /*+ Recovery data for element                         +*/
      int                   vertpart;             /*+ Vertex part                                       +*/
      Gnum                  ncmpcut2;             /*+ Number of neighbor nodes in separator             +*/
      Gnum                  ncmpgain2;            /*+ Separator gain if moved to given part             +*/
      Gnum                  ncmpgaindlt;          /*+ Node load imbalance if element swapped            +*/
    } elem;
    struct {                                      /*+ Recovery data for node                            +*/
      int                   vertpart;             /*+ Vertex part                                       +*/
      Gnum                  ecmpsize0;            /*+ Number of element neighbors in part 0             +*/
    } node;
  } data;
} VmeshSeparateFmSave;

/*
**  The function prototypes.
*/

#ifndef VMESH_SEPARATE_FM
#define static
#endif

int                         vmeshSeparateFm     (Vmesh * restrict const, const VmeshSeparateFmParam * restrict const);

static VmeshSeparateFmElement * vmeshSeparateFmTablGet (GainTabl * const, const Gnum, const Gnum);
static int                  vmeshSeparateFmResize (GainTabl * restrict const, VmeshSeparateFmElement * restrict * const, VmeshSeparateFmNode * restrict * const, VmeshSeparateFmSave * restrict * const, const Gnum, VmeshSeparateFmElement **, VmeshSeparateFmElement **, const Gnum);
#ifdef SCOTCH_DEBUG_VMESH3
static int                  vmeshSeparateFmCheck (const Vmesh * const, const VmeshSeparateFmElement * restrict, const VmeshSeparateFmNode * restrict, const Gnum, const Gnum, const Gnum);
#endif /* SCOTCH_DEBUG_VMESH3 */

#undef static
