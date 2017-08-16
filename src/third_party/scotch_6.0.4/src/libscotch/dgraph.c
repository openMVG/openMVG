/* Copyright 2007,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Francois CHATENET (P0.0)                **/
/**                Sebastien FOUCAULT (P0.0)               **/
/**                Nicolas GICQUEL (P0.1)                  **/
/**                Jerome LACOSTE (P0.1)                   **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                graph data structure handling           **/
/**                routines.                               **/
/**                                                        **/
/**    DATES     : # Version P0.0 : from : 01 apr 1997     **/
/**                                 to     01 apr 1997     **/
/**                # Version P0.1 : from : 12 apr 1998     **/
/**                                 to     20 jun 1998     **/
/**                # Version 5.0  : from : 16 feb 2005     **/
/**                                 to   : 17 jul 2008     **/
/**                # Version 5.1  : from : 21 jun 2008     **/
/**                                 to   : 30 jul 2010     **/
/**                # Version 6.0  : from : 12 sep 2012     **/
/**                                 to   : 12 sep 2012     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH

#include "module.h"
#include "common.h"
#include "dgraph.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* source graphs.                    */
/*                                   */
/*************************************/

/* This routine initializes a distributed graph
** structure. In order to avoid collective
** communication whenever possible, the allocation
** of send and receive index arrays is not performed
** in the routine itself, but rather delegated to
** subsequent routines such as dgraphBuild.
** However, these arrays will not be freed by
** dgraphFree, but by dgraphExit.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphInit (
Dgraph * restrict const     grafptr,              /* Distributed graph structure                    */
MPI_Comm                    proccomm)             /* Communicator to be used for all communications */
{
  memSet (grafptr, 0, sizeof (Dgraph));           /* Clear public and private graph fields */

  grafptr->proccomm = proccomm;                   /* Set private fields    */
  MPI_Comm_size (proccomm, &grafptr->procglbnbr); /* Get communicator data */
  MPI_Comm_rank (proccomm, &grafptr->proclocnum);

  return (0);
}

/* This routine frees the public and private data
** of the given distributed graph, but not its
** communicator.
** Private data could have been kept and freed only
** in dgraphExit(). Yet, freeing it along with the
** graph public data is a way to avoid memory
** fragmentation. Moreover, it allows to use compact
** private data if some graphs are known not to have
** holes.
** It is not a collective routine, as no communication
** is needed to perform the freeing of memory structures.
** It returns:
** - VOID  : in all cases.
*/

static
void
dgraphFree2 (
Dgraph * restrict const     grafptr)
{
  if ((grafptr->flagval & DGRAPHFREETABS) != 0) { /* If local arrays must be freed */
    if (grafptr->vertloctax != NULL)
      memFree (grafptr->vertloctax + grafptr->baseval);
    if ((grafptr->flagval & DGRAPHVERTGROUP) == 0) { /* If vertex arrays not grouped */
      if (grafptr->vendloctax != (grafptr->vertloctax + 1))
	memFree (grafptr->vendloctax + grafptr->baseval);
      if (grafptr->veloloctax != NULL)
	memFree (grafptr->veloloctax + grafptr->baseval);
      if (grafptr->vnumloctax != NULL)
        memFree (grafptr->vnumloctax + grafptr->baseval);
      if (grafptr->vlblloctax != NULL)
	memFree (grafptr->vlblloctax + grafptr->baseval);
    }
    if (grafptr->edgeloctax != NULL)
      memFree (grafptr->edgeloctax + grafptr->baseval);
    if ((grafptr->flagval & DGRAPHEDGEGROUP) == 0) { /* If edge arrays not grouped */
      if (grafptr->edloloctax != NULL)
	memFree (grafptr->edloloctax + grafptr->baseval);
    }
  }
  if ((grafptr->flagval & DGRAPHFREEPSID) != 0) { /* If process send arrays must be freed */
    if (grafptr->procsidtab != NULL)
      memFree (grafptr->procsidtab);
  }
  if ((grafptr->flagval & DGRAPHFREEEDGEGST) != 0) { /* If ghost array must be freed */
    if (grafptr->edgegsttax != NULL)
      memFree (grafptr->edgegsttax + grafptr->baseval);
  }
  if ((grafptr->flagval & DGRAPHFREEPRIV) != 0)   /* If private data has to be freed */
    if (grafptr->procdsptab != NULL)
      memFree (grafptr->procdsptab);              /* Free group leader of graph private data */
}

void
dgraphFree (
Dgraph * restrict const     grafptr)
{
  DgraphFlag          flagval;
  MPI_Comm            proccomm;                   /* Data for temporarily saving private data */
  int                 procglbnbr;
  int                 proclocnum;

  dgraphFree2 (grafptr);                          /* Free all user fields */

  flagval    = grafptr->flagval & DGRAPHFREECOMM;
  proccomm   = grafptr->proccomm;                 /* Save private fields only */
  procglbnbr = grafptr->procglbnbr;
  proclocnum = grafptr->proclocnum;

  memSet (grafptr, 0, sizeof (Dgraph));           /* Reset graph structure */

  grafptr->flagval    = flagval;                  /* Restore private fields */
  grafptr->proccomm   = proccomm;
  grafptr->procglbnbr = procglbnbr;
  grafptr->proclocnum = proclocnum;

  return;
}

/* This routine destroys a distributed graph structure.
** It is not a collective routine, as no communication
** is needed to perform the freeing of memory structures.
** Private data are always destroyed. If this is not
** wanted, use dgraphFree() instead.
** It returns:
** - VOID  : in all cases.
*/

void
dgraphExit (
Dgraph * restrict const     grafptr)
{
  if ((grafptr->flagval & DGRAPHFREECOMM) != 0)   /* If communicator has to be freed */
    MPI_Comm_free (&grafptr->proccomm);           /* Free it                         */

  dgraphFree2 (grafptr);

#ifdef SCOTCH_DEBUG_DGRAPH1
  memSet (grafptr, 0, sizeof (Dgraph));
#endif /* SCOTCH_DEBUG_DGRAPH1 */
}
