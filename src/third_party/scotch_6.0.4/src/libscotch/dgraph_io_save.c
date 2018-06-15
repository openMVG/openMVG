/* Copyright 2007,2010,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_io.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the input/output routines for       **/
/**                distributed graphs.                     **/
/**                                                        **/
/**                # Version P0.2 : from : 11 may 1999     **/
/**                                 to     12 may 1999     **/
/**                # Version 5.0  : from : 22 jul 2005     **/
/**                                 to   : 22 apr 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to   : 11 aug 2010     **/
/**                # Version 6.0  : from : 18 sep 2013     **/
/**                                 to   : 18 sep 2013     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_IO

#include "module.h"
#include "common.h"
#include "dgraph.h"

/* This routine saves a distributed source
** graph to the given streams.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphSave (
Dgraph * restrict const     grafptr,              /* Not const since halo may update structure */
FILE * const                stream)
{
  Gnum * restrict     vlblgsttax;                 /* Based index to ghost label array */
  Gnum                vertlocnum;
  char                propstr[4];                 /* Property string                  */
  int                 o;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Barrier (grafptr->proccomm) != MPI_SUCCESS) { /* Synchronize for debugging */
    errorPrint ("dgraphSave: communication error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  vlblgsttax = NULL;                              /* Assume ghost label array is free         */
  if ((grafptr->vlblloctax != NULL) ||            /* If graph has vertex labels or            */
      (grafptr->edgeloctax == NULL) ||            /* If no global index edge array present or */
      (grafptr->procvrttab[grafptr->procglbnbr] != grafptr->procdsptab[grafptr->procglbnbr])) { /* If graph may have holes in its numbering */
    if (dgraphGhst (grafptr) != 0) {              /* Compute ghost edge array */
      errorPrint ("dgraphSave: cannot compute ghost edge array");
      return     (1);
    }
    if ((vlblgsttax = (Gnum *) memAlloc (grafptr->vertgstnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("dgraphSave: out of memory");
      return     (1);
    }

    if (grafptr->vlblloctax != NULL)
      memCpy (vlblgsttax, grafptr->vlblloctax + grafptr->baseval, grafptr->vertlocnbr * sizeof (Gnum));
    else {
      for (vertlocnum = 0; vertlocnum < grafptr->vertlocnbr; vertlocnum ++) /* vlblgsttax is not based yet at this time */
        vlblgsttax[vertlocnum] = grafptr->procvrttab[grafptr->proclocnum] + vertlocnum;
    }

    if (dgraphHaloSync (grafptr, (byte *) vlblgsttax, GNUM_MPI) != 0) { /* vlblgsttax is not based yet at this time */
      errorPrint ("dgraphSave: cannot halo labels");
      memFree    (vlblgsttax);
      return     (1);
    }
    vlblgsttax -= grafptr->baseval;
  }

  propstr[0] = (vlblgsttax != NULL) ? '1' : '0';  /* Set property string */
  propstr[1] = (grafptr->edloloctax != NULL) ? '1' : '0';
  propstr[2] = (grafptr->veloloctax != NULL) ? '1' : '0';
  propstr[3] = '\0';

  if (fprintf (stream, "2\n" GNUMSTRING "\t" GNUMSTRING "\n" GNUMSTRING "\t" GNUMSTRING "\n" GNUMSTRING "\t" GNUMSTRING "\n" GNUMSTRING "\t%3s\n", /* Write file header */
               (Gnum) grafptr->procglbnbr,
               (Gnum) grafptr->proclocnum,
               (Gnum) grafptr->vertglbnbr,
               (Gnum) grafptr->edgeglbnbr,
               (Gnum) grafptr->vertlocnbr,
               (Gnum) grafptr->edgelocnbr,
               (Gnum) grafptr->baseval,
               propstr) == EOF) {
    errorPrint ("dgraphSave: bad output (1)");
    return     (1);
  }

  o = 0;
  for (vertlocnum = grafptr->baseval; (vertlocnum < grafptr->vertlocnnd) && (o == 0); vertlocnum ++) {
    Gnum                edgelocnum;

    if (vlblgsttax != NULL)                       /* Write vertex label if necessary */
      o  = (fprintf (stream, GNUMSTRING "\t", (Gnum) vlblgsttax[vertlocnum]) == EOF);
    if (grafptr->veloloctax != NULL)              /* Write vertex load if necessary */
      o |= (fprintf (stream, GNUMSTRING "\t", (Gnum) grafptr->veloloctax[vertlocnum]) == EOF);

    o |= (fprintf (stream, GNUMSTRING, (Gnum) (grafptr->vendloctax[vertlocnum] - grafptr->vertloctax[vertlocnum])) == EOF); /* Write vertex degree */

    for (edgelocnum = grafptr->vertloctax[vertlocnum];
         edgelocnum < grafptr->vendloctax[vertlocnum]; edgelocnum ++) {
      o |= (putc ('\t', stream) == EOF);
      if (grafptr->edloloctax != NULL)            /* Write edge load if necessary */
        o |= (fprintf (stream, "\t" GNUMSTRING " ", (Gnum) grafptr->edloloctax[edgelocnum]) == EOF);
      o |= (fprintf (stream, GNUMSTRING, (Gnum) ((vlblgsttax != NULL) /* Write edge end */
                                                   ? vlblgsttax[grafptr->edgegsttax[edgelocnum]]
                                                   : grafptr->edgeloctax[edgelocnum])) == EOF);
    }
    o |= (putc ('\n', stream) == EOF);
  }

  if (o != 0)
    errorPrint ("dgraphSave: bad output (2)");

  if (vlblgsttax != NULL)                         /* Free ghost label array if used */
    memFree (vlblgsttax + grafptr->baseval);

  return (o);
}
