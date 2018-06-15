/* Copyright 2008,2010,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dmapping_io.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                mappings.                               **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 13 jun 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 29 oct 2014     **/
/**                                 to     29 oct 2014     **/
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
#include "dgraph_allreduce.h"
#include "arch.h"
#include "dmapping.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

DGRAPHALLREDUCEMAXSUMOP (1, 5)

/* This routine saves a distributed mapping.
** The distributed graph structure is provided
** to access the distribution of vertex labels,
** whenever present.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dmapSave (
const Dmapping * restrict const dmapptr,
const Dgraph * restrict const   grafptr,
FILE * restrict const           stream)
{
  const DmappingFrag * restrict fragptr;
  Gnum                          fragglbnbr;
  Gnum * restrict               termloctab;
  Gnum * restrict               termrcvtab;
  Gnum                          vertrcvmax;
  Gnum                          vertglbnbr;
  Gnum * restrict               vlbltax;
  Gnum                          reduloctab[6];
  Gnum                          reduglbtab[6];
  int                           protnum;

  reduloctab[0] = dmapptr->vertlocmax;
  reduloctab[1] = dmapptr->vertlocnbr;
  reduloctab[2] = dmapptr->fragnbr;
  if (stream != NULL) {                           /* If file provided         */
    reduloctab[3] = 1;                            /* This process is the root */
    reduloctab[4] = grafptr->proclocnum;          /* Get its rank             */
  }
  else {
    reduloctab[3] =                               /* This process is not the root */
    reduloctab[4] = 0;
  }
  reduloctab[5] = (grafptr->vlblloctax != NULL) ? 1 : 0; /* See if vertex labels provided */

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 1, 5, grafptr->proccomm) != 0) {
    errorPrint ("dmapSave: communication error (1)");
    return     (1);
  }
  if (reduglbtab[3] != 1) {
    errorPrint ("dmapSave: should have only one root");
    return     (1);
  }
  if ((reduglbtab[5] != 0) && (reduglbtab[5] != grafptr->procglbnbr)) {
    errorPrint ("dmapSave: inconsistent parameters");
    return     (1);
  }
  if ((reduglbtab[1] < 0) && (reduglbtab[1] > grafptr->procglbnbr)) {
    errorPrint ("dmapSave: invalid mapping (1)");
    return     (1);
  }
  vertrcvmax = reduglbtab[0];                     /* Size of largest fragment to receive */
  vertglbnbr = reduglbtab[1];
  fragglbnbr = reduglbtab[2];
  protnum    = (int) reduglbtab[4];               /* Get rank of root process */

  reduloctab[0] = 0;
  if (protnum == grafptr->proclocnum) {
    Gnum                  vlblnbr;

    vlblnbr = (grafptr->vlblloctax != NULL) ? grafptr->vertglbnbr : 0;
    if ((termloctab = memAllocGroup ((void **) (void *) /* termloctab not used on root processor, but used only for freeing the block              */
                                     &termrcvtab, (size_t) (vertrcvmax * 2 * sizeof (Gnum)), /* TRICK: "*2" as vnumrcvtab is sent after termrcvtab */
                                     &vlbltax,    (size_t) (vlblnbr        * sizeof (Gnum)), NULL)) == NULL) {
      errorPrint ("dmapSave: out of memory (1)");
      reduloctab[0] = 1;
    }
    else if (fprintf (stream, GNUMSTRING "\n", (Gnum) vertglbnbr) == EOF) {
      errorPrint ("dmapSave: bad output (1)");
      reduloctab[0] = 1;
    }
  }
  else {
    vlbltax = NULL;                               /* Prevent Valgrind from yelling */
    if ((termloctab = memAlloc (dmapptr->vertlocmax * sizeof (Gnum))) == NULL) {
      errorPrint ("dmapSave: out of memory (2)");
      reduloctab[0] = 1;
    }
  }
#ifdef SCOTCH_DEBUG_DMAP1                         /* This communication cannot be covered by a useful one */
  if (MPI_Allreduce (reduloctab, reduglbtab, 1, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmapSave: communication error (2)");
    reduglbtab[0] = 1;
  }
#else /* SCOTCH_DEBUG_DMAP1 */
  reduglbtab[0] = reduloctab[0];
#endif /* SCOTCH_DEBUG_DMAP1 */
  if (reduglbtab[0] != 0) {
    if (termloctab != NULL)
      memFree (termloctab);                       /* Free group leader */
    return (1);
  }

  if (grafptr->vlblloctax != NULL) {
    if (commGatherv (grafptr->vlblloctax + grafptr->baseval, grafptr->vertlocnbr, GNUM_MPI,
                     vlbltax, grafptr->proccnttab, grafptr->procdsptab, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dmapSave: communication error (3)");
      return     (1);
    }
    vlbltax -= grafptr->baseval;                  /* Base label array */
  }

  if (protnum == grafptr->proclocnum) {
    Gnum                  vertrcvnbr;
    Gnum * restrict       vnumrcvptr;
    Gnum * restrict       termrcvptr;

    for (fragptr = dmapptr->fragptr; fragptr != NULL; fragptr = fragptr->nextptr) { /* Output local fragments */
      Gnum                  fraglocnum;

      for (fraglocnum = 0; fraglocnum < fragptr->vertnbr; fraglocnum ++) {
        Gnum                  vnumnum;
        Gnum                  termnum;

        vnumnum = fragptr->vnumtab[fraglocnum];
#ifdef SCOTCH_DEBUG_DMAP2
        if ((vnumnum < 0) || (vnumnum >= (grafptr->vertglbnbr + grafptr->baseval))) {
          errorPrint ("dmapSave: invalid mapping (2)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DMAP2 */
        termnum = archDomNum (&dmapptr->archdat, &fragptr->domntab[fragptr->parttab[fraglocnum]]);

        if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                     (Gnum) ((grafptr->vlblloctax != NULL) ? vlbltax[vnumnum] : vnumnum),
                     (Gnum) termnum) == EOF) {
          errorPrint ("dmapSave: bad output (2)");
          reduloctab[0] = 1;
          break;
        }
      }
    }

    for (fragglbnbr -= dmapptr->fragnbr; fragglbnbr > 0; fragglbnbr --) { /* For all non-local fragments */
      Gnum * restrict     termrcvnnd;
      MPI_Status          statdat;
      int                 recvnbr;

      if (MPI_Recv (termrcvtab, (int) (vertrcvmax * 2), GNUM_MPI, MPI_ANY_SOURCE, MPI_ANY_TAG, grafptr->proccomm, &statdat) != MPI_SUCCESS) {
        errorPrint ("dmapSave: communication error (4)"); /* TRICK: "*2" as vnumrcvtab is sent after termrcvtab */
        return     (1);
      }

      if (reduloctab[0] != 0)
        continue;

      MPI_Get_count (&statdat, GNUM_MPI, &recvnbr);
      vertrcvnbr = (Gnum) (recvnbr / 2);          /* We received a composite message made of both vectors   */
      vnumrcvptr = termrcvtab + vertrcvnbr;       /* Vertex index array is just after terminal number array */

      for (termrcvptr = termrcvtab, termrcvnnd = termrcvtab + vertrcvnbr; termrcvptr < termrcvnnd; termrcvptr ++, vnumrcvptr ++) {
        if (fprintf (stream, GNUMSTRING "\t" GNUMSTRING "\n",
                     (Gnum) ((grafptr->vlblloctax != NULL) ? vlbltax[*vnumrcvptr] : *vnumrcvptr),
                     (Gnum) *termrcvptr) == EOF) {
          errorPrint ("dmapSave: bad output (3)");
          reduloctab[0] = 1;
          break;
        }
      }
    }
  }
  else {
    int                 typecnttab[2];
    MPI_Aint            typedsptab[2];
    MPI_Datatype        typedat;

    for (fragptr = dmapptr->fragptr; fragptr != NULL; fragptr = fragptr->nextptr) { /* Output local fragments */
      Gnum                fraglocnum;

      for (fraglocnum = 0; fraglocnum < fragptr->vertnbr; fraglocnum ++) {
#ifdef SCOTCH_DEBUG_DMAP2
        Gnum                vnumnum;

        vnumnum = fragptr->vnumtab[fraglocnum];
        if ((vnumnum < 0) || (vnumnum >= (grafptr->vertglbnbr + grafptr->baseval))) {
          errorPrint ("dmapSave: invalid mapping (3)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_DMAP2 */
        termloctab[fraglocnum] = archDomNum (&dmapptr->archdat, &fragptr->domntab[fragptr->parttab[fraglocnum]]);
      }

#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
      MPI_Address (termloctab, &typedsptab[0]);
      MPI_Address (fragptr->vnumtab, &typedsptab[1]);
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
      MPI_Get_address (termloctab, &typedsptab[0]);
      MPI_Get_address (fragptr->vnumtab, &typedsptab[1]);
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
      typedsptab[1] -= typedsptab[0];
      typedsptab[0] = 0;
      typecnttab[0] =
      typecnttab[1] = (int) fragptr->vertnbr;
#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
      MPI_Type_hindexed (2, typecnttab, typedsptab, GNUM_MPI, &typedat);
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
      MPI_Type_create_hindexed (2, typecnttab, typedsptab, GNUM_MPI, &typedat);
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
      MPI_Type_commit   (&typedat);

      if (MPI_Send (termloctab, 1, typedat, protnum, 0, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("dmapSave: communication error (5)");
        return     (1);
      }

      MPI_Type_free (&typedat);
    }
  }

  memFree (termloctab);                           /* Free group leader */

#ifdef SCOTCH_DEBUG_DMAP1                         /* This communication cannot be covered by a useful one */
  if (MPI_Allreduce (reduloctab, reduglbtab, 1, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmapSave: communication error (6)");
    reduglbtab[0] = 1;
  }
#else /* SCOTCH_DEBUG_DMAP1 */
  reduglbtab[0] = reduloctab[0];
#endif /* SCOTCH_DEBUG_DMAP1 */

  return  ((int) reduglbtab[0]);
}
