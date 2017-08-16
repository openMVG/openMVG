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
/**   NAME       : fax.h                                   **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Jean ROMAN (v0.0)                       **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic factorization routine. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     24 sep 1998     **/
/**                # Version 0.1  : from : 04 apr 1999     **/
/**                                 to     01 may 1999     **/
/**                # Version 1.0  : from : 01 jun 2002     **/
/**                                 to     25 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     25 sep 2002     **/
/**                # Version 1.3  : from : 17 jun 2003     **/
/**                                 to     17 jul 2003     **/
/**                # Version 2.0  : from : 21 mar 2003     **/
/**                                 to     29 oct 2003     **/
/**                # Version 2.0  : from : 03 mar 2004     **/
/**                                 to     03 mar 2004     **/
/**                # Version 3.0  : from : 23 nov 2004     **/
/**                                 to     03 mar 2005     **/
/**                                                        **/
/************************************************************/

#define FAX_H

/*
**  The function prototypes.
*/

#ifndef FAX
#define static
#endif

int                         symbolCompact       (SymbolMatrix * const symbptr);
int                         symbolFax           (SymbolMatrix * const symbptr, const INT vertnbr, const INT edgenbr, const INT baseval, void * const ngbdptr, INT ngbfrst (void * const, const INT), INT ngbnext (void * const), INT ngbdegr (void * const, const INT), const Order * const ordeptr);
#ifdef GRAPH_H
int                         symbolFaxGraph      (SymbolMatrix * const symbptr, const Graph * const grafptr, const Order * const ordeptr);
#endif /* GRAPH_H */
int                         symbolFaxGrid2C     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT baseval, const Order * const ordeptr);
int                         symbolFaxGrid2D     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT baseval, const Order * const ordeptr);
int                         symbolFaxGrid2E     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT baseval, const Order * const ordeptr);
int                         symbolFaxGrid3C     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const Order * const ordeptr);
int                         symbolFaxGrid3D     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const Order * const ordeptr);
int                         symbolFaxGrid3E     (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const Order * const ordeptr);
#ifdef MESH_H
int                         symbolFaxMesh       (SymbolMatrix * const symbptr, const Mesh * const meshptr, const Order * const ordeptr);
#endif /* MESH_H */

int                         symbolFaxi          (SymbolMatrix * const symbptr, const INT vertnbr, const INT edgenbr, const INT baseval, void * const ngbdptr, INT ngbfrst (void * const, const INT), INT ngbnext (void * const), INT ngbdegr (void * const, const INT), const Order * const ordeptr, const INT levfmax);
#ifdef GRAPH_H
int                         symbolFaxiGraph     (SymbolMatrix * const symbptr, const Graph * const grafptr, const Order * const ordeptr, const INT levfmax);
#endif /* GRAPH_H */
int                         symbolFaxiGrid2D    (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT baseval, const Order * const ordeptr, const INT levfmax);
int                         symbolFaxiGrid2E    (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT baseval, const Order * const ordeptr, const INT levfmax);
int                         symbolFaxiGrid3D    (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const Order * const ordeptr, const INT levfmax);
int                         symbolFaxiGrid3E    (SymbolMatrix * const symbptr, const INT xnbr, const INT ynbr, const INT znbr, const INT baseval, const Order * const ordeptr, const INT levfmax);

#undef static
