/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dof.c                                   **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the general purpose     **/
/**                routines for the DOF structure.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 oct 1998     **/
/**                                 to     14 oct 1998     **/
/**                # Version 3.0  : from : 28 feb 2004     **/
/**                                 to     03 feb 2006     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DOF

#include "common.h"
#ifdef SCOTCH_PTSCOTCH
#include "ptscotch.h"
#else /* SCOTCH_PTSCOTCH */
#include "scotch.h"
#endif /* SCOTCH_PTSCOTCH */
#include "graph.h"
#include "dof.h"

/******************************/
/*                            */
/* The DOF handling routines. */
/*                            */
/******************************/

/*+ This routine initializes
*** the given DOF structure.
*** It returns:
*** - 0  : in all cases.
+*/

int
dofInit (
Dof * const                 deofptr)
{
  deofptr->baseval = 0;
  deofptr->nodenbr = 0;
  deofptr->noddval = 1;                           /* Set constant, non zero, number of DOFs */
  deofptr->noddtab = NULL;

  return (0);
}

/*+ This routine frees the contents
*** of the given DOF structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
dofExit (
Dof * const                 deofptr)
{
  if (deofptr->noddtab != NULL)
    memFree (deofptr->noddtab);

#ifdef DOF_DEBUG
  dofInit (deofptr);
#endif /* DOF_DEBUG */
}

/*+ This routine sets the number of DOFs
*** per node to a constant value.
*** It returns:
*** - VOID  : in all cases.
+*/

void
dofConstant (
Dof * const                 deofptr,
const INT                   baseval,
const INT                   nodenbr,
const INT                   noddval)
{
  deofptr->baseval = baseval;
  deofptr->nodenbr = nodenbr;
  if (deofptr->noddtab != NULL) {                 /* If DOF array already allocated */
    memFree (deofptr->noddtab);                   /* It is no longer of use         */
    deofptr->noddtab = NULL;
  }
  deofptr->noddval = noddval;
}

/*+ This routine builds the DOF index
*** array from the graph vertex array.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
dofGraph (
Dof * const                 deofptr,              /*+ DOF index array to build [based]            +*/
const Graph * const         grafptr,              /*+ Matrix adjacency structure [based]          +*/
const INT                   deofval,              /*+ DOFs per node if no graph vertex load array +*/
const INT * const           peritab)              /*+ Inverse vertex->node permutation array      +*/
{
  INT                 baseval;
  INT                 vertnbr;
  INT *               velotab;
  INT                 edgenbr;

  SCOTCH_graphData (grafptr, &baseval, &vertnbr, NULL, NULL, &velotab, NULL, &edgenbr, NULL, NULL);

  deofptr->baseval = baseval;
  deofptr->nodenbr = vertnbr;
  if (velotab == NULL) {                          /* If no vertex weight (i.e. DOF) array */
    deofptr->noddtab = NULL;                      /* No DOF array                         */
    deofptr->noddval = deofval;                   /* Get node DOF value                   */
  }
  else {                                          /* Vertex load array present */
#ifdef DOF_CONSTANT
    deofptr->noddtab = NULL;                      /* No DOF array */
    deofptr->noddval = deofval;
#else /* DOF_CONSTANT */
    const INT * restrict  velotax;                /* Based access to grafptr->velotab  */
    INT                   nodenum;                /* Number of current node            */
    INT *                 noddtnd;                /* Pointer to end of DOF index array */
    INT *                 noddptr;                /* Pointer to current DOF index      */
    const INT *           periptr;

    deofptr->noddval = 0;                         /* DOF values are not constant */
    if ((deofptr->noddtab = (INT *) memAlloc ((vertnbr + 1) * sizeof (INT))) == NULL) {
      errorPrint ("dofGraph: out of memory");
      return     (1);
    }
    for (noddptr = deofptr->noddtab, noddtnd = noddptr + vertnbr,
         periptr = peritab, nodenum = baseval,
         velotax = velotab - baseval;
         noddptr < noddtnd; noddptr ++, periptr ++) {
      *noddptr = nodenum;                         /* Set index to DOF array        */
      nodenum += velotax[*periptr];               /* Add number of DOFs for vertex */
    }
    *noddptr = nodenum;                           /* Set end of DOF array */
#endif /* DOF_CONSTANT */
  }

  return (0);
}
