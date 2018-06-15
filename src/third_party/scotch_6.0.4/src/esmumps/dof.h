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
/**   NAME       : dof.h                                   **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the DOF handling structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 oct 1998     **/
/**                                 to     16 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

#define DOF_H

#define DOF_CONSTANT                              /* Constant DOFs for ESMUMPS */

/*
**  The type and structure definitions.
*/

/*+ The DOF structure. This structure is
    always associated to a Graph structure,
    which holds the base value.             +*/

typedef struct Dof_ {
  INT                       baseval;              /*+ Base value for indexing                                       +*/
  INT                       nodenbr;              /*+ Number of nodes in DOF array                                  +*/
  INT                       noddval;              /*+ DOF value for every node (if noddtab == NULL, 0 else)         +*/
  INT * restrict            noddtab;              /*+ Array of node->first DOF indexes (if noddval == 0) [+1,based] +*/
} Dof;

/*
**  The function prototypes.
*/

#ifndef DOF
#define static
#endif

int                         dofInit             (Dof * const deofptr);
void                        dofExit             (Dof * const deofptr);
int                         dofLoad             (Dof * const deofptr, FILE * const stream);
int                         dofSave             (const Dof * const deofptr, FILE * const stream);
void                        dofConstant         (Dof * const deofptr, const INT baseval, const INT nodenbr, const INT noddval);
#ifdef GRAPH_H
int                         dofGraph            (Dof * const deofptr, const Graph * grafptr, const INT, const INT * const peritab);
#endif /* GRAPH_H */

#undef static

/*
**  The macro definitions.
*/

#ifdef DOF_CONSTANT
#define noddVal(deofptr,nodenum)    ((deofptr)->baseval + (deofptr)->noddval * ((nodenum) - (deofptr)->baseval))
#define noddDlt(deofptr,nodenum)    ((deofptr)->noddval)
#else /* DOF_CONSTANT */
#define noddVal(deofptr,nodenum)    (((deofptr)->noddtab != NULL) ? (deofptr)->noddtab[(deofptr)->baseval + (nodenum)] : ((deofptr)->baseval + (deofptr)->noddval * ((nodenum) - (deofptr)->baseval)))
#define noddDlt(deofptr,nodenum)    (((deofptr)->noddtab != NULL) ? ((deofptr)->noddtab[(deofptr)->baseval + (nodenum) + 1] - (deofptr)->noddtab[(deofptr)->baseval + (nodenum)]) : (deofptr)->noddval)
#endif /* DOF_CONSTANT */
