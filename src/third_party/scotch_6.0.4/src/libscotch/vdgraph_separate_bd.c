/* Copyright 2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : vdgraph_separate_bd.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER                        **/
/**                                                        **/
/**   FUNCTION   : This module computes a separator of the **/
/**                given distributed separator graph by    **/
/**                creating a band graph of given witdh    **/
/**                around the current separator, computing **/
/**                an improved separator of the band       **/
/**                graph, and projecting back the obtained **/
/**                separator in the original graph.        **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 04 mar 2006     **/
/**                                 to   : 07 nov 2007     **/
/**                # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 01 mar 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VDGRAPH_SEPARATE_BD

#include "module.h"
#include "common.h"
#include "parser.h"
#include "dgraph.h"
#include "vdgraph.h"
#include "vdgraph_separate_bd.h"
#include "vdgraph_separate_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed band graph
** of given width around the current separator and
** applies distributed separation routines to it.
** The distributed graph is not guaranteed to be
** balanced at at all.
** It returns:
** - 0   : if the distributed band graph could be computed.
** - !0  : on error.
*/

int
vdgraphSeparateBd (
Vdgraph * const                       grafptr,    /*+ Distributed graph +*/
const VdgraphSeparateBdParam * const  paraptr)    /*+ Method parameters +*/
{
  Vdgraph                 bandgrafdat;            /* Vertex separator band graph structure                   */
  Gnum                    bandvertancnnd;         /* End of local vertex array, without anchors              */
  Gnum                    bandvertlocnbr1;        /* Number of band graph vertices in part 1 except anchor 1 */
  Gnum                    bandvertlocnum;
  Gnum                    bandvertlocancadj;      /* Flag set when anchor(s) represent unexistent vertices   */
  Gnum                    bandvertglbancadj;      /* Global adjustment of anchor vertices                    */
  Gnum                    complocsizeadj0;
  Gnum                    complocsizeadj1;
  Gnum                    reduloctab[3];
  Gnum                    reduglbtab[3];
  Gnum * restrict         edloloctax;             /* Save value for edge loads while we pretend we don't have them */
  Gnum                    fronlocnum;

  if (grafptr->compglbsize[2] == 0)               /* If no frontier to base on  */
    return (0);                                   /* Then do nothing            */
  if (paraptr->distmax < 1)                       /* If distance is 0 (or less) */
    return (0);                                   /* Then do nothing            */

  edloloctax = grafptr->s.edloloctax;             /* Fake no edge loads on original graph as we do not need them */
  grafptr->s.edloloctax = NULL;
  if (dgraphBand (&grafptr->s, grafptr->complocsize[2], grafptr->fronloctab, grafptr->partgsttax,
                  grafptr->complocload[0] + grafptr->complocload[2], grafptr->complocload[1], paraptr->distmax,
                  &bandgrafdat.s, &bandgrafdat.fronloctab, &bandgrafdat.partgsttax,
                  NULL, &bandvertlocnbr1, &bandvertlocancadj) != 0) {
    grafptr->s.edloloctax = edloloctax;
    errorPrint ("vdgraphSeparateBd: cannot create band graph");
    return     (1);
  }
  grafptr->s.edloloctax = edloloctax;             /* Restore edge loads, if any */

  bandgrafdat.complocsize[0] = bandgrafdat.s.vertlocnbr - (bandvertlocnbr1 + 1) - grafptr->complocsize[2]; /* Add 1 for anchor vertex 1 */
  bandgrafdat.complocsize[1] = bandvertlocnbr1 + 1; /* Add 1 for anchor vertex 1 */
  complocsizeadj0 = grafptr->complocsize[0] - bandgrafdat.complocsize[0];
  complocsizeadj1 = grafptr->complocsize[1] - bandgrafdat.complocsize[1];
  reduloctab[0] = bandgrafdat.complocsize[0];
  reduloctab[1] = bandgrafdat.complocsize[1];
  reduloctab[2] = bandvertlocancadj;              /* Sum increases in size and load */
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 3, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateBd: communication error (1)");
    return     (1);
  }
  bandvertglbancadj = reduglbtab[2];
  bandgrafdat.compglbload[0] = grafptr->compglbload[0] + bandvertglbancadj; /* All loads are kept in band graph */
  bandgrafdat.compglbload[1] = grafptr->compglbload[1] + bandvertglbancadj;
  bandgrafdat.compglbload[2] = grafptr->compglbload[2];
  bandgrafdat.compglbloaddlt = grafptr->compglbloaddlt; /* Balance is not changed by anchor vertices */
  bandgrafdat.complocload[0] = grafptr->complocload[0] + bandvertlocancadj;
  bandgrafdat.complocload[1] = grafptr->complocload[1] + bandvertlocancadj;
  bandgrafdat.complocload[2] = grafptr->complocload[2];
  bandgrafdat.compglbsize[0] = reduglbtab[0];
  bandgrafdat.compglbsize[1] = reduglbtab[1];
  bandgrafdat.compglbsize[2] = grafptr->compglbsize[2]; /* All separator vertices are kept in band graph */
  bandgrafdat.complocsize[2] = grafptr->complocsize[2];
  bandgrafdat.levlnum        = grafptr->levlnum;

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (&bandgrafdat) != 0) {
    errorPrint ("vdgraphSeparateBd: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  if (vdgraphSeparateSt (&bandgrafdat, paraptr->strat) != 0) { /* Separate distributed band graph */
    errorPrint  ("vdgraphSeparateBd: cannot separate band graph");
    vdgraphExit (&bandgrafdat);
    return      (1);
  }

  bandvertancnnd = bandgrafdat.s.vertlocnnd - 2;
  reduloctab[0] = ((bandgrafdat.partgsttax[bandvertancnnd]     == 0) && /* Check if anchor vertices remain in their parts */
                   (bandgrafdat.partgsttax[bandvertancnnd + 1] == 1)) ? 0 : 1;
  reduloctab[1] = bandgrafdat.complocsize[0] + complocsizeadj0;
  reduloctab[2] = bandgrafdat.complocsize[1] + complocsizeadj1;
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 3, GNUM_MPI, MPI_SUM, grafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("vdgraphSeparateBd: communication error (2)");
    return     (1);
  }
  if (reduglbtab[0] != 0) {                       /* If at least one anchor changed of part */
    vdgraphExit (&bandgrafdat);                   /* Then keep original partition           */
    return      (0);
  }

  grafptr->compglbload[0] = bandgrafdat.compglbload[0] - bandvertglbancadj;
  grafptr->compglbload[1] = bandgrafdat.compglbload[1] - bandvertglbancadj;
  grafptr->compglbload[2] = bandgrafdat.compglbload[2];
  grafptr->compglbloaddlt = bandgrafdat.compglbloaddlt;
  grafptr->compglbsize[0] = reduglbtab[1];
  grafptr->compglbsize[1] = reduglbtab[2];
  grafptr->compglbsize[2] = bandgrafdat.compglbsize[2];
  grafptr->complocload[0] = bandgrafdat.complocload[0] - bandvertlocancadj;
  grafptr->complocload[1] = bandgrafdat.complocload[1] - bandvertlocancadj;
  grafptr->complocload[2] = bandgrafdat.complocload[2];
  grafptr->complocsize[0] = reduloctab[1];
  grafptr->complocsize[1] = reduloctab[2];
  grafptr->complocsize[2] = bandgrafdat.complocsize[2];
  for (fronlocnum = 0; fronlocnum < bandgrafdat.complocsize[2]; fronlocnum ++) /* Project back separator */
    grafptr->fronloctab[fronlocnum] = bandgrafdat.s.vnumloctax[bandgrafdat.fronloctab[fronlocnum]];
  for (bandvertlocnum = bandgrafdat.s.baseval; bandvertlocnum < bandvertancnnd; bandvertlocnum ++) /* For all vertices except anchors */
    grafptr->partgsttax[bandgrafdat.s.vnumloctax[bandvertlocnum]] = bandgrafdat.partgsttax[bandvertlocnum];

#ifdef SCOTCH_DEBUG_VDGRAPH2
  if (vdgraphCheck (grafptr) != 0) {
    errorPrint ("vdgraphSeparateBd: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VDGRAPH2 */

  vdgraphExit (&bandgrafdat);

  return (0);
}
