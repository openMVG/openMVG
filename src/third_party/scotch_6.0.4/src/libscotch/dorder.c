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
/**   NAME       : dorder.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 18 apr 2006     **/
/**                                 to     28 jul 2006     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dorder.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine initializes a distributed
** ordering with respect to the given parameters.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dorderInit (
Dorder * restrict const     ordeptr,
const Gnum                  baseval,
const Gnum                  vnodglbnbr,
MPI_Comm                    proccomm)
{
  ordeptr->baseval    = baseval;
  ordeptr->vnodglbnbr = vnodglbnbr;
  ordeptr->cblklocnbr = 0;

  ordeptr->linkdat.nextptr = &ordeptr->linkdat;   /* Loop double-chained list */
  ordeptr->linkdat.prevptr = &ordeptr->linkdat;

  MPI_Comm_dup  (proccomm, &ordeptr->proccomm);   /* Duplicate communicator to avoid lifespan problems */
  MPI_Comm_rank (ordeptr->proccomm, &ordeptr->proclocnum);

#ifdef SCOTCH_PTHREAD
  pthread_mutex_init (&ordeptr->mutelocdat, NULL); /* Initialize local mutex */
#endif /* SCOTCH_PTHREAD */

  return (0);
}

/* This routine frees the column blocks
** of the given distributed ordering.
** It returns:
** - void  : in all cases.
*/

static
void
dorderFreeCblk (
DorderCblk * restrict const cblkptr)
{
#ifdef SCOTCH_DEBUG_DORDER2
  if ((cblkptr->typeval < DORDERCBLKNEDI) ||
      (cblkptr->typeval > (DORDERCBLKNEDI | DORDERCBLKLEAF)))
    errorPrint ("dorderFreeCblk: invalid column block type");
#endif /* SCOTCH_DEBUG_DORDER2 */

  if ((cblkptr->typeval & DORDERCBLKLEAF) != 0) {
    memFree (cblkptr->data.leaf.periloctab);
    if (cblkptr->data.leaf.nodeloctab != NULL)
      memFree (cblkptr->data.leaf.nodeloctab);
  }

  memFree (cblkptr);                              /* Free column block structure */
}

void
dorderFree (
Dorder * restrict const     ordeptr)
{
  DorderCblk *          cblkptr;
  DorderLink *          linkptr;

  for (linkptr = ordeptr->linkdat.nextptr; linkptr != &ordeptr->linkdat; ) {
    cblkptr = (DorderCblk *) linkptr;             /* TRICK: FIRST */
    linkptr = linkptr->nextptr;

    dorderFreeCblk (cblkptr);
  }

  ordeptr->linkdat.nextptr =                      /* Loop double-chained list */
  ordeptr->linkdat.prevptr = &ordeptr->linkdat;
}

/* This routine frees the contents
** of the given ordering.
** It returns:
** - void  : in all cases.
*/

void
dorderExit (
Dorder * restrict const     ordeptr)
{
  dorderFree (ordeptr);

  MPI_Comm_free (&ordeptr->proccomm);             /* Free duplicated communicator */

#ifdef SCOTCH_PTHREAD
  pthread_mutex_destroy (&ordeptr->mutelocdat);   /* Destroy local mutex */
#endif /* SCOTCH_PTHREAD */

#ifdef SCOTCH_DEBUG_DORDER2
  memSet (ordeptr, ~0, sizeof (Dorder));
#endif /* SCOTCH_DEBUG_DORDER2 */
}

/* This routine creates the root column
** block slot in the given distributed
** ordering structure.
** It returns:
** - !NULL  : root column block.
** - NULL   : on error.
*/

DorderCblk *
dorderFrst (
Dorder * const              ordeptr)
{
  DorderCblk            cblkdat;
  DorderCblk *          cblkptr;

  cblkdat.ordelocptr         = ordeptr;           /* Fake father node                                       */
  cblkdat.cblknum.proclocnum = 0;                 /* Belongs to process 0 to ease displacement computations */
  cblkdat.cblknum.cblklocnum = -1;

  if ((cblkptr = dorderNew (&cblkdat, ordeptr->proccomm)) == NULL)
    return (NULL);

  cblkptr->ordeglbval = 0;                        /* Un-based inverse permutation index */
  cblkptr->vnodglbnbr = ordeptr->vnodglbnbr;
  cblkptr->cblkfthnum = 0;

  return (cblkptr);
}

/* This routine gives back a new distributed
** column block slot in the same ordering
** structure as the given column block. 
** It returns:
** - !NULL  : new column block.
** - NULL   : on error.
*/

DorderCblk *
dorderNew (
DorderCblk * const          cblkptr,              /* One of the column blocks       */
MPI_Comm                    proccomm)             /* Communicator sharing the block */
{
  Dorder * restrict     ordeptr;
  DorderCblk * restrict cblknewptr;
  Gnum                  reduloctab[3];
  Gnum                  reduglbtab[3];
  int                   proclocnum;

  MPI_Comm_rank (proccomm, &proclocnum);

  ordeptr = cblkptr->ordelocptr;

  reduloctab[1] =                                 /* Assume process is not root for this column block */
  reduloctab[2] = 0;
  if ((cblknewptr = (DorderCblk *) memAlloc (sizeof (DorderCblk))) == NULL) {
    errorPrint ("dorderNew: out of memory");
    reduloctab[0] = 2;                            /* Indicate error without doubt */
  }
  else {
    reduloctab[0] = 0;                            
    if (proclocnum == 0) {                        /* If root of sub-tree                 */
      reduloctab[0] = 1;                          /* Indicate it is the root             */
      reduloctab[1] = ordeptr->proclocnum;        /* Broadcast global rank of block root */
#ifdef SCOTCH_PTHREAD
      pthread_mutex_lock (&ordeptr->mutelocdat);  /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */
      reduloctab[2] = ordeptr->cblklocnbr ++;     /* One more root block in local ordering */
#ifdef SCOTCH_PTHREAD
      pthread_mutex_unlock (&ordeptr->mutelocdat); /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */
    }
  }
  if (MPI_Allreduce (&reduloctab, &reduglbtab, 3, GNUM_MPI, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderNew: communication error");
    return     (NULL);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("dorderNew: cannot create new node");
    if (cblknewptr != NULL)
      memFree (cblknewptr);
    return (NULL);
  }

  cblknewptr->ordelocptr         = ordeptr;
  cblknewptr->typeval            = DORDERCBLKNONE;
  cblknewptr->fathnum            = cblkptr->cblknum;
  cblknewptr->cblknum.proclocnum = (int) reduglbtab[1];
  cblknewptr->cblknum.cblklocnum = reduglbtab[2];

#ifdef SCOTCH_PTHREAD
  pthread_mutex_lock (&ordeptr->mutelocdat);      /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */
  cblknewptr->linkdat.nextptr = &ordeptr->linkdat; /* Link new block at end of local ordering node list */
  cblknewptr->linkdat.prevptr = ordeptr->linkdat.prevptr;
  ordeptr->linkdat.prevptr->nextptr = &cblknewptr->linkdat;
  ordeptr->linkdat.prevptr = &cblknewptr->linkdat;
#ifdef SCOTCH_PTHREAD
  pthread_mutex_unlock (&ordeptr->mutelocdat);    /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */

  return (cblknewptr);
}

/* This routine gives back a new centralized
** column block slot in the same ordering
** structure as the given column block.
** It returns:
** - !NULL  : new column block.
** - NULL   : on error.
*/

DorderCblk *
dorderNewSequ (
DorderCblk * const          cblkptr)              /* One of the column blocks */
{
  Dorder * restrict     ordeptr;
  DorderCblk * restrict cblknewptr;

  if ((cblknewptr = (DorderCblk *) memAlloc (sizeof (DorderCblk))) == NULL) {
    errorPrint ("dorderNewSequ: out of memory");
    return     (NULL);
  }

  ordeptr = cblkptr->ordelocptr;

  cblknewptr->ordelocptr         = ordeptr;
  cblknewptr->typeval            = DORDERCBLKNONE;
  cblknewptr->fathnum            = cblkptr->cblknum;
  cblknewptr->cblknum.proclocnum = ordeptr->proclocnum; /* Node belongs to this process */
#ifdef SCOTCH_PTHREAD
  pthread_mutex_lock (&ordeptr->mutelocdat);      /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */
  cblknewptr->cblknum.cblklocnum = ordeptr->cblklocnbr ++; /* One more locally-rooted block in ordering */

  cblknewptr->linkdat.nextptr = &ordeptr->linkdat; /* Link new block at end of local ordering node list */
  cblknewptr->linkdat.prevptr = ordeptr->linkdat.prevptr;
  ordeptr->linkdat.prevptr->nextptr = &cblknewptr->linkdat;
  ordeptr->linkdat.prevptr = &cblknewptr->linkdat;
#ifdef SCOTCH_PTHREAD
  pthread_mutex_unlock (&ordeptr->mutelocdat);    /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */

  return (cblknewptr);
}

/* This routine gives back a new centralized
** index range in the same ordering structure
** as the given column block.
** It returns:
** - !NULL  : new column block.
** - NULL   : on error.
*/

Gnum
dorderNewSequIndex (
DorderCblk * const          cblkptr,              /* One of the column blocks     */
const Gnum                  cblknbr)              /* Number of indices to reserve */
{
  Dorder * restrict     ordeptr;
  Gnum                  cblklocnum;

  ordeptr = cblkptr->ordelocptr;

#ifdef SCOTCH_PTHREAD
  pthread_mutex_lock (&ordeptr->mutelocdat);      /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */
  cblklocnum = ordeptr->cblklocnbr;               /* Get current local index number           */
  ordeptr->cblklocnbr += cblknbr;                 /* These more root blocks in local ordering */
#ifdef SCOTCH_PTHREAD
  pthread_mutex_unlock (&ordeptr->mutelocdat);    /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */

  return (cblklocnum);
}

/* This routine removes a no longer used
** column block from the given distributed
** ordering. Leaves or locally-rooted column
** blocks are kept, others are removed.
** It returns:
** - void  : in all cases.
*/

void
dorderDispose (
DorderCblk * const          cblkptr)              /* Column block to consider */
{
  Dorder * restrict     ordeptr;

  ordeptr = cblkptr->ordelocptr;

  if (cblkptr->cblknum.proclocnum == ordeptr->proclocnum) /* If node is local root of column block, keep it */
    return;

  if ((cblkptr->typeval & DORDERCBLKLEAF) == 0) { /* If node is not non-rooted leaf of distributed ordering */
#ifdef SCOTCH_PTHREAD
    pthread_mutex_lock (&ordeptr->mutelocdat);    /* Lock local mutex */
#endif /* SCOTCH_PTHREAD */
    cblkptr->linkdat.nextptr->prevptr = cblkptr->linkdat.prevptr; /* Unchain node from double-chained list  */
    cblkptr->linkdat.prevptr->nextptr = cblkptr->linkdat.nextptr;
#ifdef SCOTCH_PTHREAD
    pthread_mutex_unlock (&ordeptr->mutelocdat);  /* Unlock local mutex */
#endif /* SCOTCH_PTHREAD */

    memFree (cblkptr);
  }
}
