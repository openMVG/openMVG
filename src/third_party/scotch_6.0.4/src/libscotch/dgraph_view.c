/* Copyright 2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_view.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed graph general       **/
/**                purpose routines.                       **/
/**                                                        **/
/**    DATES     : # Version P0.0 : from : 01 apr 1997     **/
/**                                 to     01 apr 1997     **/
/**                # Version P0.1 : from : 12 apr 1998     **/
/**                                 to     20 jun 1998     **/
/**                # Version 5.0  : from : 16 feb 2005     **/
/**                                 to   : 15 aug 2006     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to   : 12 aug 2010     **/
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

/* This routine displays the contents
** of the given graph structure.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphView (
const Dgraph * restrict const grafptr,
FILE * const                  stream)
{
  MPI_Comm            proccomm;                   /* Graph communicator                     */
  int                 procglbnbr;                 /* Number of processes sharing graph data */
  int                 proclocnum;                 /* Number of this process                 */
  int                 procngbnbr;
  int                 procngbnum;
  Gnum                vertlocnum;
  Gnum                edgelocnum;
  Gnum *              edgelocptr;

  proccomm = grafptr->proccomm;                   /* Simplify                  */
  MPI_Comm_size (proccomm, &procglbnbr);          /* Rely on communicator data */
  MPI_Comm_rank (proccomm, &proclocnum);

  fflush (stream);                                /* Flush previous data */
  for (procngbnbr = 0; procngbnbr < procglbnbr; procngbnbr ++) {
    MPI_Barrier (proccomm);
    if (procngbnbr == proclocnum) {
      fprintf (stream, "Process %d:\n",
	       proclocnum);
      fprintf (stream, "  vertglbnbr: " GNUMSTRING "\n  vertgstnbr: " GNUMSTRING "\n vertgstnnd: " GNUMSTRING "\n  vertlocnbr: " GNUMSTRING "\n vertlocnnd: " GNUMSTRING "\n",
	       (Gnum) grafptr->vertglbnbr,
	       (Gnum) grafptr->vertgstnbr,
	       (Gnum) grafptr->vertgstnnd,
	       (Gnum) grafptr->vertlocnbr,
	       (Gnum) grafptr->vertlocnnd);
      fprintf (stream, "  vertloctax:");
      if (grafptr->vendloctax == grafptr->vertloctax + 1) {
	for (vertlocnum = grafptr->baseval; vertlocnum <= grafptr->vertlocnnd; vertlocnum ++)/**/
	  fprintf (stream, " " GNUMSTRING,
		   (Gnum) grafptr->vertloctax[vertlocnum]);
	fprintf (stream, " x\n  vendloctax: = vertloctax + 1");
      }
      else {
	for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++)
	  fprintf (stream, " " GNUMSTRING,
		   (Gnum) grafptr->vertloctax[vertlocnum]);
	fprintf (stream, "  vendloctax: x");
	for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++)
	  fprintf (stream, " " GNUMSTRING,
		   (Gnum) grafptr->vendloctax[vertlocnum]);
      }
      fprintf (stream, "\n  edgeglbnbr: " GNUMSTRING "\n  edgelocnbr: " GNUMSTRING "\n",
	       (Gnum) grafptr->edgeglbnbr,
	       (Gnum) grafptr->edgelocnbr);
      fprintf (stream, "  edgeloctax:");
      for (edgelocnum = grafptr->baseval, edgelocptr = grafptr->edgeloctax;
	   edgelocnum < grafptr->edgelocnbr + grafptr->baseval;
	   edgelocnum ++, edgelocptr ++)
	fprintf (stream, " " GNUMSTRING,
		 (Gnum) *edgelocptr);
      if ((grafptr->flagval & DGRAPHHASEDGEGST) != 0) {
	fprintf (stream, "\n  edgegsttax:");
	for (edgelocnum = grafptr->baseval, edgelocptr = grafptr->edgegsttax;
	     edgelocnum < grafptr->edgelocnbr + grafptr->baseval;
	     edgelocnum ++, edgelocptr ++)
	  fprintf (stream, " " GNUMSTRING,
		   (Gnum) *edgelocptr);
      }
      fprintf (stream, "\n  procdsptab:");
      for (procngbnum = 0; procngbnum <= procglbnbr ; procngbnum ++)
	fprintf (stream, " " GNUMSTRING,
		 (Gnum) grafptr->procdsptab[procngbnum]);
      fprintf (stream, "\n  procngbnbr: %d",
	       grafptr->procngbnbr);
      fprintf (stream, "\n  procngbtab:");
      for (procngbnum = 0; procngbnum < grafptr->procngbnbr; procngbnum ++)
	fprintf (stream, " %d",
		 grafptr->procngbtab[procngbnum]);
      fprintf (stream, "\n  procrcvtab:");
      for (procngbnum = 0; procngbnum < grafptr->procglbnbr; procngbnum ++)
	fprintf (stream, " %d",
		 grafptr->procrcvtab[procngbnum]);
      fprintf (stream, "\n  procsndnbr: %d",
	       grafptr->procsndnbr);
      fprintf (stream, "\n  procsndtab:");
      for (procngbnum = 0; procngbnum < grafptr->procglbnbr; procngbnum ++)
	fprintf (stream, " %d",
		 grafptr->procsndtab[procngbnum]);
      fprintf (stream, "\n  degrglbmax: " GNUMSTRING,
	       (Gnum) grafptr->degrglbmax);
      fprintf (stream, "\n");
      fflush  (stream);                           /* Flush data */
    }
  }
  MPI_Barrier (proccomm);

  return (0);
}
