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
/**   NAME       : order.h                                 **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the graph ordering routine.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 aug 1998     **/
/**                                 to     01 may 1999     **/
/**                # Version 2.0  : from : 25 oct 2003     **/
/**                                 to     02 jul 2010     **/
/**                                                        **/
/************************************************************/

#define ORDER_H

/*
**  The type and structure definitions.
*/

/*+ Ordering structure. vnodbas holds the base
    value for node indexings. vnodbas is equal
    to baseval for graphs, and to vnodbas for
    meshes. The same holds for rangtab, with
    rangtab[0] = vnodbas.                      +*/

typedef struct Order_ {
  INT                       cblknbr;              /*+ Number of column blocks             +*/
  INT *                     rangtab;              /*+ Column block range array [based,+1] +*/
  INT *                     permtab;              /*+ Permutation array [based]           +*/
  INT *                     peritab;              /*+ Inverse permutation array [based]   +*/
} Order;

/*
**  The function prototypes.
*/

#ifndef ORDER
#define static
#endif

int                         orderInit           (Order * const ordeptr);
void                        orderExit           (Order * const ordeptr);
int                         orderLoad           (Order * const ordeptr, FILE * const stream);
int                         orderSave           (const Order * const ordeptr, FILE * const stream);
void                        orderBase           (Order * restrict const ordeptr, const INT baseval);

int                         orderCheck          (const Order * const ordeptr);

int                         orderGrid2          (Order * const ordeptr, const INT xnbr, const INT ynbr, const INT baseval, const INT xlim, const INT ylim);
int                         orderGrid2C         (Order * const ordeptr, const INT xnbr, const INT ynbr, const INT baseval, const INT xlim, const INT ylim);
int                         orderGrid3          (Order * const ordeptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const INT xlim, const INT ylim, const INT zlim);
int                         orderGrid3C         (Order * const ordeptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const INT xlim, const INT ylim, const INT zlim);

#ifdef GRAPH_H
int                         orderGraph          (Order * restrict const ordeptr, Graph * restrict const grafptr);
int                         orderGraphList      (Order * restrict const ordeptr, Graph * restrict const grafptr, const INT listnbr, const INT * restrict const listtab);
int                         orderGraphStrat     (Order * restrict const ordeptr, Graph * restrict const grafptr, const char * restrict const);
int                         orderGraphListStrat (Order * restrict const ordeptr, Graph * restrict const grafptr, const INT listnbr, const INT * restrict const listtab, const char * const);
#endif /* GRAPH_H */

#ifdef MESH_H
int                         orderMesh           (Order * restrict const ordeptr, Mesh * restrict const meshptr);
int                         orderMeshList       (Order * restrict const ordeptr, Mesh * restrict const meshptr, const INT listnbr, const INT * restrict const listtab);
int                         orderMeshStrat      (Order * restrict const ordeptr, Mesh * restrict const meshptr, const char * const);
int                         orderMeshListStrat  (Order * restrict const ordeptr, Mesh * restrict const meshptr, const INT listnbr, const INT * restrict const listtab, const char * const);
#endif /* MESH_H */

#undef static
