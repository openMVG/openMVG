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
/**   NAME       : dorder_io.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 27 apr 2006     **/
/**                                 to     13 jun 2008     **/
/**                # Version 5.1  : from : 30 jul 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "comm.h"
#include "dgraph.h"
#include "dorder.h"
#include "order.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine saves a distributed ordering.
** The distributed graph structure is provided
** to access the distribution of vertex labels,
** whenever present.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dorderSave (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
FILE * restrict const         stream)
{
  Gnum * restrict       peritab;
  Gnum * restrict       permtab;
  Gnum * restrict       vlbltax;
  int                   procglbnbr;
  int                   reduloctab[3];
  int                   reduglbtab[3];
  int                   protnum;

  if (stream != NULL) {                           /* If file provided         */
    reduloctab[0] = 1;                            /* This process is the root */
    reduloctab[1] = ordeptr->proclocnum;          /* Get its rank             */
  }
  else {
    reduloctab[0] =                               /* This process is not the root */
    reduloctab[1] = 0;
  }
  reduloctab[2] = (grafptr->vlblloctax != NULL) ? 1 : 0; /* See if vertex labels provided */
  if (MPI_Allreduce (reduloctab, reduglbtab, 3, MPI_INT, MPI_SUM, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderSave: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 1) {
    errorPrint ("dorderSave: should have only one root");
    return     (1);
  }
  MPI_Comm_size (ordeptr->proccomm, &procglbnbr);
  if ((reduglbtab[2] != 0) && (reduglbtab[2] != procglbnbr)) {
    errorPrint ("dorderSave: inconsistent parameters");
    return     (1);
  }
  protnum = (int) reduglbtab[1];                  /* Get rank of root process */

  reduloctab[0] = 0;
  permtab       = NULL;
  if (protnum == ordeptr->proclocnum) {
    Gnum                  vlblnbr;

    vlblnbr = (grafptr->vlblloctax != NULL) ? ordeptr->vnodglbnbr : 0;
    if (memAllocGroup ((void **) (void *)
                       &permtab, (size_t) (ordeptr->vnodglbnbr * sizeof (Gnum)),
                       &peritab, (size_t) (ordeptr->vnodglbnbr * sizeof (Gnum)),
                       &vlbltax, (size_t) (vlblnbr             * sizeof (Gnum)), NULL) == NULL) {
      errorPrint ("dorderSave: out of memory");
#ifdef SCOTCH_DEBUG_DORDER1
      reduloctab[0] = 1;
#else /* SCOTCH_DEBUG_DORDER1 */
      return (1);
#endif /* SCOTCH_DEBUG_DORDER1 */
    }
  }
#ifdef SCOTCH_DEBUG_DORDER1                       /* This communication cannot be covered by a useful one */
  if (MPI_Bcast (&reduloctab[0], 1, MPI_INT, protnum, ordeptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dorderSave: communication error (2)");
    if (permtab != NULL);
      memFree (permtab);                          /* Free group leader */
    return (1);
  }
  if (reduloctab[0] != 0)
    return (1);
#endif /* SCOTCH_DEBUG_DORDER1 */

  if (grafptr->vlblloctax != NULL) {
    if (commGatherv (grafptr->vlblloctax + grafptr->baseval, grafptr->vertlocnbr, GNUM_MPI,
                     vlbltax, grafptr->proccnttab, grafptr->procdsptab, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dorderSave: communication error (3)");
      return     (1);
    }
  }

  if (protnum == ordeptr->proclocnum) {
    Gnum                  vertnum;

    for (vertnum = 0; vertnum < ordeptr->vnodglbnbr; ) { /* Till all inverse permutation indices collected */
      const DorderLink *    linkptr;

      for (linkptr = ordeptr->linkdat.nextptr; linkptr != &ordeptr->linkdat; linkptr = linkptr->nextptr) {
        const DorderCblk *    cblkptr;

        cblkptr = (DorderCblk *) linkptr;         /* TRICK: FIRST */
        if (((cblkptr->typeval & DORDERCBLKLEAF) != 0) &&
            (cblkptr->data.leaf.ordelocval == vertnum) && /* If column block fragment starts at proper index       */
            (cblkptr->data.leaf.vnodlocnbr > 0)) { /* And is not an empty local block with relevent data elsewhere */
          memCpy (peritab + vertnum, cblkptr->data.leaf.periloctab, cblkptr->data.leaf.vnodlocnbr * sizeof (Gnum));
          vertnum += cblkptr->data.leaf.vnodlocnbr;
          break;
        }
      }
      if (linkptr == &ordeptr->linkdat) {         /* If fragment not found locally */
        MPI_Status            statdat;
        int                   recvnbr;

        if (MPI_Bcast (&vertnum, 1, GNUM_MPI, protnum, ordeptr->proccomm) != MPI_SUCCESS) { /* Broadcast missing fragment */
          errorPrint ("dorderSave: communication error (4)");
          memFree    (permtab);                       /* Free group leader */
          return     (1);
        }
        if (MPI_Recv (peritab + vertnum, ordeptr->vnodglbnbr - vertnum, GNUM_MPI,
                      MPI_ANY_SOURCE, DORDERTAGPERI, ordeptr->proccomm, &statdat) != MPI_SUCCESS) {
          errorPrint ("dorderSave: communication error (5)");
          return     (1);
        }
        MPI_Get_count (&statdat, GNUM_MPI, &recvnbr);
        vertnum += recvnbr;
      }
    }
    vertnum = -1;                                 /* Indicate termination */
    if (MPI_Bcast (&vertnum, 1, GNUM_MPI, protnum, ordeptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dorderSave: communication error (6)");
      memFree    (permtab);                       /* Free group leader */
      return     (1);
    }

    if (fprintf (stream, GNUMSTRING "\n", (Gnum) ordeptr->vnodglbnbr) == EOF) {
      errorPrint ("dorderSave: bad output (1)");
      memFree    (permtab);
      return     (1);
    }

    orderPeri (peritab, ordeptr->baseval, ordeptr->vnodglbnbr, permtab, ordeptr->baseval); /* Compute direct permutation */

    if (grafptr->vlblloctax != NULL) {            /* If ordering has label array */
      vlbltax -= ordeptr->baseval;                /* Base label array            */

      for (vertnum = 0; vertnum < ordeptr->vnodglbnbr; vertnum ++) {
        if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                     (Gnum) vlbltax[vertnum + ordeptr->baseval],
                     (Gnum) vlbltax[permtab[vertnum]]) == EOF) {
          errorPrint ("dorderSave: bad output (2)");
          memFree    (permtab);
          return     (1);
        }
      }
    }
    else {
      for (vertnum = 0; vertnum < ordeptr->vnodglbnbr; vertnum ++) {
        if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                     (Gnum) (vertnum + ordeptr->baseval),
                     (Gnum) permtab[vertnum]) == EOF) {
          errorPrint ("dorderSave: bad output (3)");
          memFree    (permtab);
          return     (1);
        }
      }
    }

    memFree (permtab);                            /* Free group leader */
  }
  else {
    while (1) {
      const DorderLink *    linkptr;
      Gnum                  vertnum;

      if (MPI_Bcast (&vertnum, 1, GNUM_MPI, protnum, ordeptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dorderSave: communication error (7)");
        return     (1);
      }
      if (vertnum == -1)                          /* If asked to quit */
        break;                                    /* Finish           */

      for (linkptr = ordeptr->linkdat.nextptr; linkptr != &ordeptr->linkdat; linkptr = linkptr->nextptr) {
        const DorderCblk *    cblkptr;

        cblkptr = (DorderCblk *) linkptr;         /* TRICK: FIRST */

        if (((cblkptr->typeval & DORDERCBLKLEAF) != 0) && /* If matching column block fragment found */
            (cblkptr->data.leaf.ordelocval == vertnum) &&
            (cblkptr->data.leaf.vnodlocnbr > 0)) { /* And is not an empty local block with relevent data elsewhere */
          if (MPI_Send (cblkptr->data.leaf.periloctab, cblkptr->data.leaf.vnodlocnbr,
                        GNUM_MPI, protnum, DORDERTAGPERI, ordeptr->proccomm) != MPI_SUCCESS) {
            errorPrint ("dorderSave: communication error (8)");
            return     (1);
          }
          break;
        }
      }
    }
  }

  return  (0);
}
