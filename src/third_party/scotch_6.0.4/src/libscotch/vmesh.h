/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vmesh.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for mesh vertex separation        **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 10 sep 2002     **/
/**                                 to   : 10 sep 2002     **/
/**                # Version 5.1  : from : 04 nov 2010     **/
/**                                 to   : 04 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Active graph structure. +*/

typedef struct Vmesh_ {
  Mesh                      m;                    /*+ Source mesh                                        +*/
  GraphPart *               parttax;              /*+ Based part array: 0,1: part; 2: separator          +*/
  Gnum                      ecmpsize[2];          /*+ Number of elements in each part (not in separator) +*/
  Gnum                      ncmpload[3];          /*+ Loads of nodes in both parts and separator         +*/
  Gnum                      ncmploaddlt;          /*+ Node load difference between both parts            +*/
  Gnum                      ncmpsize[2];          /*+ Number of nodes in parts (separator is fronnbr)    +*/
  Gnum                      fronnbr;              /*+ Number of frontier nodes; TRICK: ncmpsize[2]       +*/
  Gnum *                    frontab;              /*+ Array of frontier node numbers                     +*/
  Gnum                      levlnum;              /*+ Nested dissection or coarsening level              +*/
} Vmesh;

/*+ The graph separator storing structure. +*/

typedef struct VmeshStore_ {
  Gnum                      ecmpsize[2];          /*+ Number of elements in each part                 +*/
  Gnum                      ncmpload[3];          /*+ Loads of nodes in both parts and separator      +*/
  Gnum                      ncmploaddlt;          /*+ Node load difference between both parts         +*/
  Gnum                      ncmpsize[2];          /*+ Number of nodes in parts (separator is fronnbr) +*/
  Gnum                      fronnbr;              /*+ Number of frontier nodes; TRICK: ncmpsize[2]    +*/
  byte *                    datatab;              /*+ Variable-sized data array                       +*/
} VmeshStore;

/*
**  The function prototypes.
*/

#ifndef VMESH
#define static
#endif

void                        vmeshExit           (Vmesh * const);
void                        vmeshZero           (Vmesh * const);
int                         vmeshCheck          (const Vmesh * const);

int                         vmeshStoreInit      (const Vmesh * const, VmeshStore * const);
void                        vmeshStoreExit      (VmeshStore * const);
void                        vmeshStoreSave      (const Vmesh * const , VmeshStore * const);
void                        vmeshStoreUpdt      (Vmesh * const, const VmeshStore * const);

#undef static
