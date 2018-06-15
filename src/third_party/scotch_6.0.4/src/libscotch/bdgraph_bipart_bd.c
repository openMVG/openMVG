/* Copyright 2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : bdgraph_bipart_bd.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a bipartition of   **/
/**                the given distributed separator graph   **/
/**                by creating a band graph of given witdh **/
/**                around the current frontier, computing  **/
/**                an improved bipartition of the band     **/
/**                graph, and projecting back the obtained **/
/**                frontier to the original graph.         **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 nov 2007     **/
/**                                 to   : 14 apr 2011     **/
/**                # Version 6.0  : from : 11 sep 2011     **/
/**                                 to   : 31 aug 2014     **/
/**                                                        **/
/**   NOTES      : # Since only edges from local vertices  **/
/**                  to local anchors are created in       **/
/**                  dgraphBand(), the communication const **/
/**                  might be wrong if a local vertex of   **/
/**                  the last layer is linked to a remote  **/
/**                  vertex of different part which was    **/
/**                  not in the band graph. Hence, commun- **/
/**                  ication costs have to be recomputed   **/
/**                  from scratch.                         **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define BDGRAPH_BIPART_BD

#include "module.h"
#include "common.h"
#include "parser.h"
#include "arch.h"
#include "dgraph.h"
#include "dgraph_halo.h"
#include "bdgraph.h"
#include "bdgraph_bipart_bd.h"
#include "bdgraph_bipart_st.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine computes a distributed band graph
** of given width around the current frontier and
** applies distributed bipartitioning routines to it.
** The distributed graph is not guaranteed to be
** balanced at at all.
** It returns:
** - 0   : if the distributed band graph could be computed.
** - !0  : on error.
*/

int
bdgraphBipartBd (
Bdgraph * const                     orggrafptr,   /*+ Distributed graph +*/
const BdgraphBipartBdParam * const  paraptr)      /*+ Method parameters +*/
{
  Bdgraph                 bndgrafdat;             /* Bipartitioning band graph structure                          */
  Gnum                    bndvertancnnd;          /* End of local vertex array, without anchors                   */
  Gnum                    bndvertlocnbr1;         /* Number of band graph vertices in part 1 except anchor 1      */
  Gnum                    bndvertlocnum;
  Gnum                    bndvertlvlnum;          /* Based number of first band vertex in last layer              */
  Gnum                    bndvertlocancadj;       /* Flag set when anchor(s) represent unexistent vertices        */
  Gnum                    bndvertglbancadj;       /* Global adjustment of anchor vertices                         */
  Gnum                    bndveexlocsum;          /* Local sum of veexloctax array cells for band graph           */
  Gnum                    bndveexlocsum0;         /* Local sum of veexloctax array cells in part 0 for band graph */
  Gnum                    bndedlolocval;
  Gnum                    bndfronlocnum;
  Gnum                    orgfronlocnum;
  int * restrict          orgflagloctab;
  Gnum                    orgvertlocnum;
  Gnum                    orgedlolocval;
  const int * restrict    orgprocsidtab;
  int                     orgprocsidnbr;
  int                     orgprocsidnum;
  int                     orgprocsidval;
  Gnum                    complocsizeadj0;
  Gnum                    commlocloadintn;
  Gnum                    commlocloadintn2;       /* Twice twice (4 times) the internal communication load of last layer */
  Gnum                    commlocloadextn;
  Gnum                    commlocgainextn;
  Gnum                    reduloctab[7];
  Gnum                    reduglbtab[7];
  DgraphHaloRequest       requdat;

  if (orggrafptr->fronglbnbr == 0)                /* If no separator vertices, apply strategy to full (original) graph */
    return (bdgraphBipartSt (orggrafptr, paraptr->stratorg));

  if (dgraphBand (&orggrafptr->s, orggrafptr->fronlocnbr, orggrafptr->fronloctab, orggrafptr->partgsttax,
                  orggrafptr->complocload0, orggrafptr->s.velolocsum - orggrafptr->complocload0, paraptr->distmax,
                  &bndgrafdat.s, &bndgrafdat.fronloctab, &bndgrafdat.partgsttax,
                  &bndvertlvlnum, &bndvertlocnbr1, &bndvertlocancadj) != 0) {
    errorPrint ("bdgraphBipartBd: cannot create band graph");
    return     (1);
  }
  bndvertancnnd = bndgrafdat.s.vertlocnnd - 2;

  reduloctab[0]  = 0;                             /* Assume no memory allocation problem */
  bndveexlocsum  =
  bndveexlocsum0 = 0;
  bndgrafdat.veexloctax = NULL;                   /* Assume no external gains */
  if (orggrafptr->veexloctax != NULL) {
    if ((bndgrafdat.veexloctax = memAlloc (bndgrafdat.s.vertlocnbr * sizeof (Gnum))) == NULL) {
      errorPrint ("bdgraphBipartBd: out of memory (1)");
      reduloctab[0] = 1;                         /* Memory error */
    }
    else {
      Gnum                bndvertlocnum;

      bndgrafdat.veexloctax -= bndgrafdat.s.baseval;

      for (bndvertlocnum = bndgrafdat.s.baseval; bndvertlocnum < bndvertancnnd; bndvertlocnum ++) {
        Gnum                veexval;

        veexval = orggrafptr->veexloctax[bndgrafdat.s.vnumloctax[bndvertlocnum]];
        bndgrafdat.veexloctax[bndvertlocnum] = veexval;
        bndveexlocsum  += veexval;
        bndveexlocsum0 += veexval & (((Gnum) bndgrafdat.partgsttax[bndvertlocnum]) - 1);
      }
    }
  }
  reduloctab[1] = bndgrafdat.s.vendloctax[bndvertancnnd]     - bndgrafdat.s.vertloctax[bndvertancnnd]     - (orggrafptr->s.procglbnbr - 1); /* Anchor degrees */
  reduloctab[2] = bndgrafdat.s.vendloctax[bndvertancnnd + 1] - bndgrafdat.s.vertloctax[bndvertancnnd + 1] - (orggrafptr->s.procglbnbr - 1);

  bndgrafdat.complocsize0 = bndgrafdat.s.vertlocnbr - (bndvertlocnbr1 + 1); /* Add 1 for anchor vertex 1 */
  complocsizeadj0 = orggrafptr->complocsize0 - bndgrafdat.complocsize0; /* -1 less because of anchor 0   */
  reduloctab[3] = bndgrafdat.complocsize0;
  reduloctab[4] = bndvertlocancadj;               /* Sum increases in size and load */
  reduloctab[5] = bndveexlocsum;
  reduloctab[6] = bndveexlocsum0;
  if (MPI_Allreduce (reduloctab, reduglbtab, 7, GNUM_MPI, MPI_SUM, orggrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartBd: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 0) {
    bdgraphExit (&bndgrafdat);
    return      (1);
  }
  if ((reduglbtab[1] == 0) ||                     /* If graph is too small to have any usable anchors */
      (reduglbtab[2] == 0)) {
    bdgraphExit (&bndgrafdat);
    return      (bdgraphBipartSt (orggrafptr, paraptr->stratorg));
  }

  bndvertglbancadj            = reduglbtab[4];
  bndgrafdat.veexglbsum       = orggrafptr->veexglbsum; /* All external gains preserved                  */
  bndgrafdat.fronlocnbr       = orggrafptr->fronlocnbr; /* All separator vertices are kept in band graph */
  bndgrafdat.fronglbnbr       = orggrafptr->fronglbnbr;
  bndgrafdat.complocload0     = orggrafptr->complocload0    + bndvertlocancadj; /* All loads are kept in band graph */
  bndgrafdat.compglbload0     = orggrafptr->compglbload0    + bndvertglbancadj;
  bndgrafdat.compglbload0min  = orggrafptr->compglbload0min + bndvertglbancadj; /* Tilt extrema loads according to adjustments */
  bndgrafdat.compglbload0max  = orggrafptr->compglbload0max + bndvertglbancadj;
  bndgrafdat.compglbload0avg  = orggrafptr->compglbload0avg + bndvertglbancadj; /* Tilt average load according to adjustments */
  bndgrafdat.compglbload0dlt  = orggrafptr->compglbload0dlt;
  bndgrafdat.compglbsize0     = reduglbtab[3];
  bndgrafdat.commglbload      = orggrafptr->commglbload;
  bndgrafdat.commglbgainextn  = orggrafptr->commglbgainextn;
  bndgrafdat.commglbloadextn0 = orggrafptr->commglbloadextn0;
  bndgrafdat.commglbgainextn0 = orggrafptr->commglbgainextn0;
  bndgrafdat.bbalglbval       = orggrafptr->bbalglbval;
  bndgrafdat.domndist         = orggrafptr->domndist;
  bndgrafdat.domnwght[0]      = orggrafptr->domnwght[0];
  bndgrafdat.domnwght[1]      = orggrafptr->domnwght[1];
  bndgrafdat.levlnum          = orggrafptr->levlnum;

  if (bndgrafdat.veexloctax != NULL) {
    Gnum                bndveexglbanc0;
    Gnum                bndveexglbanc1;

    bndveexglbanc0 = (orggrafptr->veexglbsum + orggrafptr->commglbgainextn) / 2 - reduglbtab[6]; /* Compute global external gains of anchors */
    bndveexglbanc1 = (orggrafptr->veexglbsum - bndveexglbanc0) - reduglbtab[5];

    bndgrafdat.veexloctax[bndvertancnnd]     = DATASIZE (bndveexglbanc0, bndgrafdat.s.procglbnbr, bndgrafdat.s.proclocnum); /* Spread gains across local anchors */
    bndgrafdat.veexloctax[bndvertancnnd + 1] = DATASIZE (bndveexglbanc1, bndgrafdat.s.procglbnbr, bndgrafdat.s.proclocnum);
  }

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (&bndgrafdat) != 0) {
    errorPrint ("bdgraphBipartBd: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  if (bdgraphBipartSt (&bndgrafdat, paraptr->stratbnd) != 0) { /* Separate distributed band graph */
    errorPrint  ("bdgraphBipartBd: cannot separate band graph");
    bdgraphExit (&bndgrafdat);
    return      (1);
  }

  reduloctab[0] = (Gnum) bndgrafdat.partgsttax[bndvertancnnd]; /* Check if anchor vertices remain in their parts */
  reduloctab[1] = (Gnum) bndgrafdat.partgsttax[bndvertancnnd + 1];
  reduloctab[2] = complocsizeadj0;
  reduloctab[3] = 0;                              /* Assume memory allocation is all right */
  if ((orgflagloctab = memAlloc (flagSize (orggrafptr->s.vertlocnnd) * sizeof (int))) == NULL) { /* Eventually keep space for based indices */
    errorPrint ("bdgraphBipartBd: out of memory (2)");
    reduloctab[3] = 1;
  }

  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 4, GNUM_MPI, MPI_SUM, orggrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartBd: communication error (2)");
    return     (1);
  }

  if (((reduglbtab[0] + reduglbtab[1]) != orggrafptr->s.procglbnbr)         || /* If not all anchors of initial same parts in same parts */
      ((reduglbtab[0] != 0) && (reduglbtab[0] != orggrafptr->s.procglbnbr)) ||
      (reduglbtab[3] != 0)) {
    if (orgflagloctab != NULL)
      memFree (orgflagloctab);
    bdgraphExit (&bndgrafdat);                    /* Apply original strategy to full graph */
    return      (bdgraphBipartSt (orggrafptr, paraptr->stratorg));
  }

  if (dgraphGhst (&bndgrafdat.s) != 0) {          /* Compute ghost edge array if not already present */
    errorPrint ("bdgraphBipartBd: cannot compute ghost edge array");
    return     (1);
  }

  if (reduglbtab[0] == orggrafptr->s.procglbnbr) { /* If all anchors swapped parts, swap all parts of original vertices */
    Gnum                orgvertnum;

    orggrafptr->complocsize0 = orggrafptr->s.vertlocnbr - reduloctab[2] - bndgrafdat.s.vertlocnbr + bndgrafdat.complocsize0;
    orggrafptr->compglbsize0 = orggrafptr->s.vertglbnbr - reduglbtab[2] - bndgrafdat.s.vertglbnbr + bndgrafdat.compglbsize0;

    for (orgvertnum = orggrafptr->s.baseval; orgvertnum < orggrafptr->s.vertlocnnd; orgvertnum ++)
      orggrafptr->partgsttax[orgvertnum] ^= 1;
  }
  else {
    orggrafptr->complocsize0 = reduloctab[2] + bndgrafdat.complocsize0;
    orggrafptr->compglbsize0 = reduglbtab[2] + bndgrafdat.compglbsize0;
  }

  for (bndvertlocnum = bndgrafdat.s.baseval; bndvertlocnum < bndvertancnnd; bndvertlocnum ++) /* Update part array of all vertices except anchors */
    orggrafptr->partgsttax[bndgrafdat.s.vnumloctax[bndvertlocnum]] = bndgrafdat.partgsttax[bndvertlocnum];

  dgraphHaloAsync (&orggrafptr->s, (byte *) (orggrafptr->partgsttax + orggrafptr->s.baseval), GRAPHPART_MPI, &requdat); /* Share part array of full graph */

  commlocloadintn =
  commlocloadextn =
  commlocgainextn = 0;
  bndedlolocval   = 1;                            /* Assume no edge loads */
  for (bndvertlocnum = bndgrafdat.s.baseval; bndvertlocnum < bndvertlvlnum; bndvertlocnum ++) { /* For all vertices of band graph save for last layer */
    Gnum                bndedgelocnum;
    Gnum                bndedgelocnnd;
    Gnum                bndpartval;

    bndpartval = (Gnum) bndgrafdat.partgsttax[bndvertlocnum];
    if (bndgrafdat.veexloctax != NULL) {
      commlocloadextn += bndgrafdat.veexloctax[bndvertlocnum] * bndpartval;
      commlocgainextn += bndgrafdat.veexloctax[bndvertlocnum] * (1 - bndpartval * 2);
    }
    for (bndedgelocnum = bndgrafdat.s.vertloctax[bndvertlocnum], bndedgelocnnd = bndgrafdat.s.vendloctax[bndvertlocnum];
         bndedgelocnum < bndedgelocnnd; bndedgelocnum ++) {
      Gnum                bndvertlocend;
      Gnum                bndpartend;

      bndvertlocend = bndgrafdat.s.edgegsttax[bndedgelocnum];
      bndpartend    = bndgrafdat.partgsttax[bndvertlocend];

      if (bndgrafdat.s.edloloctax != NULL)
        bndedlolocval = bndgrafdat.s.edloloctax[bndedgelocnum];
      commlocloadintn += (bndpartval ^ bndpartend) * bndedlolocval; /* Internal load is accounted for twice */
    }
  }
  for ( ; bndvertlocnum < bndvertancnnd; bndvertlocnum ++) { /* For all vertices of last layer, remove internal loads to band vertices once */
    Gnum                bndedgelocnum;
    Gnum                bndedgelocnnd;
    Gnum                bndpartval;

    bndpartval = (Gnum) bndgrafdat.partgsttax[bndvertlocnum];
    if (bndgrafdat.veexloctax != NULL) {
      commlocloadextn += bndgrafdat.veexloctax[bndvertlocnum] * bndpartval;
      commlocgainextn += bndgrafdat.veexloctax[bndvertlocnum] * (1 - bndpartval * 2);
    }
    for (bndedgelocnum = bndgrafdat.s.vertloctax[bndvertlocnum], bndedgelocnnd = bndgrafdat.s.vendloctax[bndvertlocnum] - 1; /* "-1" to avoid anchor edges */
         bndedgelocnum < bndedgelocnnd; bndedgelocnum ++) {
      Gnum                bndvertlocend;
      Gnum                bndpartend;

      bndvertlocend = bndgrafdat.s.edgegsttax[bndedgelocnum];
      bndpartend    = bndgrafdat.partgsttax[bndvertlocend];

      if (bndgrafdat.s.edloloctax != NULL)
        bndedlolocval = bndgrafdat.s.edloloctax[bndedgelocnum];
      commlocloadintn -= (bndpartval ^ bndpartend) * bndedlolocval; /* Remove internal loads to band graph vertices once because afterwards they will be accounted for twice */
    }
  }

  memSet (orgflagloctab, 0, flagSize (orggrafptr->s.vertlocnnd) * sizeof (int)); /* Set vertices as not already considered */

  for (bndfronlocnum = orgfronlocnum = 0; bndfronlocnum < bndgrafdat.fronlocnbr; bndfronlocnum ++) { /* Project back separator except for last layer */
    Gnum                bndvertlocnum;

    bndvertlocnum = bndgrafdat.fronloctab[bndfronlocnum];
    if (bndvertlocnum < bndvertlvlnum) {          /* If vertex does not belong to last layer */
      Gnum                orgvertlocnum;

      orgvertlocnum = bndgrafdat.s.vnumloctax[bndvertlocnum];
      flagSet (orgflagloctab, orgvertlocnum);     /* Set vertex as processed */
      orggrafptr->fronloctab[orgfronlocnum ++] = orgvertlocnum;
    }
  }

  if (dgraphHaloWait (&requdat) != 0) {
    errorPrint ("bdgraphBipartBd: cannot complete asynchronous halo exchange");
    return     (1);
  }

  orgedlolocval    = 1;                           /* Assume no edge loads */
  commlocloadintn2 = 0;
  for (bndvertlocnum = bndvertlvlnum; bndvertlocnum < bndvertancnnd; bndvertlocnum ++) { /* For all vertices of last layer */
    Gnum                orgedgelocnum;
    Gnum                orgedgelocnnd;
    Gnum                orgvertlocnum;
    GraphPart           orgpartval;
    Gnum                orgflagval;

    orgvertlocnum = bndgrafdat.s.vnumloctax[bndvertlocnum];
    orgpartval    = bndgrafdat.partgsttax[bndvertlocnum];

    orgflagval = 0;                               /* Assume vertex does not belong to the frontier */
    for (orgedgelocnum = orggrafptr->s.vertloctax[orgvertlocnum], orgedgelocnnd = orggrafptr->s.vendloctax[orgvertlocnum];
         orgedgelocnum < orgedgelocnnd; orgedgelocnum ++) {
      Gnum                orgvertlocend;
      Gnum                orgpartend;
      Gnum                orgflagtmp;

      orgvertlocend = orggrafptr->s.edgegsttax[orgedgelocnum];
      orgpartend    = orggrafptr->partgsttax[orgvertlocend];

      orgflagtmp = orgpartval ^ orgpartend;
      if (bndgrafdat.s.edloloctax != NULL)
        orgedlolocval = orggrafptr->s.edloloctax[orgedgelocnum];
      orgflagval       |= orgflagtmp;
      commlocloadintn2 += orgflagtmp * orgedlolocval; /* Internal load to band and original graph vertices are accounted for twice */
      if ((orgflagtmp != 0) && (orgvertlocend < orggrafptr->s.vertlocnnd) && (flagVal (orgflagloctab, orgvertlocend) == 0)) {
        orggrafptr->fronloctab[orgfronlocnum ++] = orgvertlocend;
        flagSet (orgflagloctab, orgvertlocend);
      }
    }
    if ((orgflagval != 0) && (flagVal (orgflagloctab, orgvertlocnum) == 0))
      orggrafptr->fronloctab[orgfronlocnum ++] = orgvertlocnum;

    flagSet (orgflagloctab, orgvertlocnum);       /* Set vertex as processed anyway */
  }
  commlocloadintn += 2 * commlocloadintn2;        /* Add twice the internal load of original graph edges and once the one of band edges (one removed before) */

  orggrafptr->complocload0    = bndgrafdat.complocload0 - bndvertlocancadj;
  orggrafptr->compglbload0    = bndgrafdat.compglbload0 - bndvertglbancadj;
  orggrafptr->compglbload0dlt = orggrafptr->compglbload0 - orggrafptr->compglbload0avg;

  orgprocsidnbr = orggrafptr->s.procsidnbr;
  if (orgprocsidnbr == 0)
    goto loop_exit;
  orgvertlocnum = orggrafptr->s.baseval;
  orgprocsidnum = 0;
  orgprocsidtab = orggrafptr->s.procsidtab;
  orgprocsidval = orgprocsidtab[orgprocsidnum ++];  
  while (1) {             /* Scan all vertices which have foreign neighbors */
    while (orgprocsidval < 0) {
      orgvertlocnum -= (Gnum) orgprocsidval;
      orgprocsidval  = orgprocsidtab[orgprocsidnum ++];  
    }

    if (flagVal (orgflagloctab, orgvertlocnum) == 0) { /* If vertex not already processed */
      Gnum                orgedgelocnum;
      Gnum                orgedgelocnnd;
      GraphPart           orgpartval;

      orgpartval = orggrafptr->partgsttax[orgvertlocnum];
      for (orgedgelocnum = orggrafptr->s.vertloctax[orgvertlocnum], orgedgelocnnd = orggrafptr->s.vendloctax[orgvertlocnum];
           orgedgelocnum < orgedgelocnnd; orgedgelocnum ++) {
        if (orggrafptr->partgsttax[orggrafptr->s.edgegsttax[orgedgelocnum]] != orgpartval) {
          orggrafptr->fronloctab[orgfronlocnum ++] = orgvertlocnum;
          break;
        }
      }
    }

    do {
      if (orgprocsidnum >= orgprocsidnbr)
        goto loop_exit;
    } while ((orgprocsidval = orgprocsidtab[orgprocsidnum ++]) >= 0);
  }
loop_exit :
  memFree (orgflagloctab);

  reduloctab[0] = commlocloadintn;                /* Twice the internal load; sum globally before dividing by two */
  reduloctab[1] = commlocloadextn;
  reduloctab[2] = commlocgainextn;
  reduloctab[3] = orgfronlocnum;
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 4, GNUM_MPI, MPI_SUM, orggrafptr->s.proccomm) != MPI_SUCCESS) {
    errorPrint ("bdgraphBipartBd: communication error (3)");
    return     (1);
  }
  orggrafptr->fronlocnbr      = orgfronlocnum;
  orggrafptr->fronglbnbr      = reduglbtab[3];
  orggrafptr->commglbload     = (reduglbtab[0] / 2) * orggrafptr->domndist + reduglbtab[1];
  orggrafptr->commglbgainextn = reduglbtab[2];
  orggrafptr->bbalglbval      = (double) ((orggrafptr->compglbload0dlt < 0) ? (- orggrafptr->compglbload0dlt) : orggrafptr->compglbload0dlt) / (double) orggrafptr->compglbload0avg;

#ifdef SCOTCH_DEBUG_BDGRAPH2
  if (bdgraphCheck (orggrafptr) != 0) {
    errorPrint ("bdgraphBipartBd: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_BDGRAPH2 */

  bdgraphExit (&bndgrafdat);

  return (0);
}
