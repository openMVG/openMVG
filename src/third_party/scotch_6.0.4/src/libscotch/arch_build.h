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
/**   NAME       : arch_build.h                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the architecture building routine.  **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 29 may 1997     **/
/**                                 to     01 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     29 nov 2003     **/
/**                                                        **/
/**   NOTES      : # This file contains pieces of code     **/
/**                  extracted from release 3.1 of         **/
/**                  "amk_src.c".                          **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Job to process. +*/

typedef struct ArchBuildJob_ {
  struct ArchBuildJob_ *    joblink;              /*+ Link to job pool                        +*/
  ArchDomNum                domnum;               /*+ Mapping domain to which vertices belong +*/
  Graph                     grafdat;              /*+ Job graph data                          +*/
} ArchBuildJob;

/*+ Vertex distance information. +*/

typedef struct ArchBuildDistElem_ {
  int                       queued;               /*+ Flag set if vertex queued  +*/
  Anum                      distval;              /*+ Distance to initial vertex +*/
} ArchBuildDistElem;

/*+ Queue element. +*/

typedef struct ArchBuildQueuElem_ {
  Gnum                      vertnum;              /*+ Vertex number in source graph +*/
  Anum                      distval;              /*+ Distance reached              +*/
} ArchBuildQueuElem;

/*
**  The function prototypes.
*/

#ifndef ARCH_BUILD
#define static
#endif

static void                 archBuildJobExit    (ArchBuildJob *);

int                         archBuild           (Arch * const, const Graph * const, const VertList * const, const Strat * const);

#undef static
