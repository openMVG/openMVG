/* Copyright 2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dmapping.h                              **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                                                        **/
/**   FUNCTION   : These lines are the declarations for    **/
/**                the parallel mapping handling routines. **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 31 mar 2008     **/
/**                                 to     04 nov 2010     **/
/**                                                        **/
/************************************************************/

#define DMAPPING_H

/*
**  The type definitions.
*/

/*+ This structure defines a mapping fragment. +*/

typedef struct DmappingFrag_ {
  struct DmappingFrag_ *    nextptr;              /*+ Pointer to next fragment            +*/
  Gnum                      vertnbr;              /*+ Number of local vertices in mapping +*/
  Gnum *                    vnumtab;              /*+ Vertex index array                  +*/
  Anum *                    parttab;              /*+ Mapping array [vertlocnbr]          +*/
  Anum                      domnnbr;              /*+ Local number of domains             +*/
  ArchDom *                 domntab;              /*+ Array of domains [domnnbr]          +*/
} DmappingFrag;

/*+ This structure defines an (eventually
    partial) mapping of a source graph to
    a target architecture.                +*/

typedef struct Dmapping_ {
  struct DmappingFrag_ *    fragptr;              /*+ Pointer to first mapping fragment   +*/
  Gnum                      fragnbr;              /*+ Number of local fragments           +*/
  Gnum                      vertlocmax;           /*+ Size of biggest local fragment      +*/
  Gnum                      vertlocnbr;           /*+ Number of local vertices in mapping +*/
  Arch                      archdat;              /*+ Architecture data                   +*/
#ifdef SCOTCH_PTHREAD
  pthread_mutex_t           mutelocdat;           /*+ Local mutex for updates             +*/
#endif /* SCOTCH_PTHREAD */
} Dmapping;

/*+ The sort structure, used to sort mapped vertices.
    Field vertnum is first and field termnum is a Gnum
    and not an Anum because of intSort2asc1.           +*/

typedef struct DmappingTermSort_ {
  Gnum                      vertnum;              /*+ Vertex number: FIRST     +*/
  Gnum                      termnum;              /*+ Direct permutation index +*/
} DmappingTermSort;

/*
**  The function prototypes.
*/

#ifndef DMAPPING
#define static
#endif

int                         dmapInit            (Dmapping * restrict const, const Arch * restrict const);
void                        dmapExit            (Dmapping * const);
void                        dmapAdd             (Dmapping * restrict const, DmappingFrag * restrict const);
int                         dmapSave            (const Dmapping * restrict const, const Dgraph * restrict const, FILE * restrict const);
int                         dmapTerm            (const Dmapping * restrict const, const Dgraph * restrict const, Gnum * restrict const);

#undef static
