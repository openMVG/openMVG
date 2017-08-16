/* Copyright 2008-2010,2012,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_match.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed graph matching      **/
/**                routines.                               **/
/**                                                        **/
/**    DATES     : # Version 5.1  : from : 01 dec 2008     **/
/**                                 to   : 30 jul 2010     **/
/**                # Version 6.0  : from : 03 oct 2012     **/
/**                                 to   : 10 oct 2013     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_MATCH

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "dgraph_match.h"

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
dgraphMatchInit (
DgraphMatchData * restrict const    mateptr,
const float                         probval)
{
  Gnum * restrict     procvgbtab;
  int                 procngbnum;
  Gnum                vertlocnbr;
  Gnum                vertgstnbr;

  Dgraph * restrict const     grafptr    = mateptr->c.finegrafptr;
  const int * restrict const  procngbtab = grafptr->procngbtab;
  const Gnum * restrict const procvrttab = grafptr->procvrttab;

  vertlocnbr = grafptr->vertlocnbr;
  vertgstnbr = grafptr->vertgstnbr;

  if (memAllocGroup ((void **) (void *)
                     &mateptr->procvgbtab, (size_t) ((grafptr->procngbnbr + 1) * sizeof (Gnum)),
                     &mateptr->queuloctab, (size_t) (vertlocnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("dgraphMatchInit: out of memory");
    return     (1);
  }

  mateptr->c.multlocnbr = 0;
  mateptr->mategsttax = mateptr->c.coargsttax;    /* TRICK: re-use array               */
  mateptr->matelocnbr = 0;                        /* All vertices need to be processed */
  mateptr->queulocnbr = 0;
  mateptr->probval = (grafptr->procngbnbr == 0) ? 1.0F : probval;

  memSet (mateptr->mategsttax + grafptr->vertlocnnd, ~0, (vertgstnbr - vertlocnbr) * sizeof (Gnum)); /* No ghost vertices matched to date */

  for (procngbnum = 0, procvgbtab = mateptr->procvgbtab; procngbnum < grafptr->procngbnbr; procngbnum ++)
    procvgbtab[procngbnum] = (Gnum) procvrttab[procngbtab[procngbnum]];
  procvgbtab[procngbnum] = (Gnum) procvrttab[grafptr->procglbnbr]; /* Mark end */

  return (0);
}

/* This routine frees the contents of a matching
** data structure.
** It returns:
** - VOID  : in all cases.
*/

void
dgraphMatchExit (
DgraphMatchData * restrict const  mateptr)
{
  memFree (mateptr->procvgbtab);
}

/* These routines perform a round of computations
** among enqueued vertices to produce matching requests.
** They return:
** - 0   : on success.
** - !0  : on error.
*/

#define DGRAPHMATCHSCANNAME         dgraphMatchSc /* Scan matching (no edge weights) */
#define DGRAPHMATCHSCANINIT                      \
  probmax = (Gnum) (mateptr->probval * 32768.0);  /* Compute integer threshold of random value */
#define DGRAPHMATCHSCANCOUNTDECL                 ;
#define DGRAPHMATCHSCANCOUNTINIT                 \
      probval = intRandVal (32768);               /* Get proba for this vertex */
#define DGRAPHMATCHSCANCOUNTSELECT               \
          edgefrenbr ++;
#define DGRAPHMATCHSCANFINDSELECT                \
          (edgefrenbr -- == 0)
#include "dgraph_match_scan.c"
#undef DGRAPHMATCHSCANFINDSELECT
#undef DGRAPHMATCHSCANCOUNTSELECT
#undef DGRAPHMATCHSCANCOUNTINIT
#undef DGRAPHMATCHSCANCOUNTDECL
#undef DGRAPHMATCHSCANINIT
#undef DGRAPHMATCHSCANNAME

#define DGRAPHMATCHSCANNAME         dgraphMatchHy /* Heavy edge matching */
#define DGRAPHMATCHSCANINIT                                     \
  const Gnum * restrict const edloloctax = grafptr->edloloctax; \
  if (edloloctax == NULL) {                                     \
    dgraphMatchSc (mateptr);                                    \
    return;                                                     \
  }                                                             \
  probmax = (Gnum) (mateptr->probval * 32768.0);  /* Compute integer threshold of random value */
#define DGRAPHMATCHSCANCOUNTDECL                                \
      Gnum                edlolocmax;
#define DGRAPHMATCHSCANCOUNTINIT                                \
      edlolocmax = 0;                                           \
      probval = intRandVal (32768);               /* Get proba for this vertex */
#define DGRAPHMATCHSCANCOUNTSELECT                              \
          Gnum                edlolocval;                       \
          edlolocval = edloloctax[edgelocnum];                  \
          if (edlolocval > edlolocmax) {                        \
            edlolocmax = edlolocval;                            \
            edgefrenbr = 1;                                     \
          }                                                     \
          else if (edlolocval == edlolocmax)                    \
            edgefrenbr ++;
#define DGRAPHMATCHSCANFINDSELECT                               \
          ((edloloctax[edgelocnum] == edlolocmax) &&            \
           (edgefrenbr -- == 0))
#include "dgraph_match_scan.c"
#undef DGRAPHMATCHSCANFINDSELECT
#undef DGRAPHMATCHSCANCOUNTSELECT
#undef DGRAPHMATCHSCANCOUNTINIT
#undef DGRAPHMATCHSCANCOUNTDECL
#undef DGRAPHMATCHSCANINIT
#undef DGRAPHMATCHSCANNAME

#define DGRAPHMATCHSCANNAME         dgraphMatchHl /* Heavy edge matching of lightest vertices */
#define DGRAPHMATCHSCANINIT                                        \
  const Gnum * restrict const veloloctax = grafptr->veloloctax;    \
  const Gnum * restrict const edloloctax = grafptr->edloloctax;    \
  if ((veloloctax == NULL) || (edloloctax == NULL)) {              \
    dgraphMatchHy (mateptr);                                       \
    return;                                                        \
  }                                                                \
  probmax = (1 * grafptr->veloglbsum) / (5 * grafptr->vertglbnbr);
#define DGRAPHMATCHSCANCOUNTDECL                                   \
      Gnum                edlolocmax;
#define DGRAPHMATCHSCANCOUNTINIT                                   \
      edlolocmax = 0;                                              \
      probval = veloloctax[vertlocnum];           /* Process vertex if vertex weight smaller than threshold */
#define DGRAPHMATCHSCANCOUNTSELECT                                 \
          Gnum                edlolocval;                          \
          edlolocval = edloloctax[edgelocnum];                     \
          if (edlolocval > edlolocmax) {                           \
            edlolocmax = edlolocval;                               \
            edgefrenbr = 1;                                        \
          }                                                        \
          else if (edlolocval == edlolocmax)                       \
            edgefrenbr ++;
#define DGRAPHMATCHSCANFINDSELECT                                  \
          ((edloloctax[edgelocnum] == edlolocmax) &&               \
           (edgefrenbr -- == 0))
#include "dgraph_match_scan.c"
#undef DGRAPHMATCHSCANFINDSELECT
#undef DGRAPHMATCHSCANCOUNTSELECT
#undef DGRAPHMATCHSCANCOUNTINIT
#undef DGRAPHMATCHSCANCOUNTDECL
#undef DGRAPHMATCHSCANINIT
#undef DGRAPHMATCHSCANNAME

#define DGRAPHMATCHSCANNAME         dgraphMatchLc /* Local scan matching */
#define DGRAPHMATCHSCANINIT             \
  probmax = 0;                                    /* Vertices will always be active */
#define DGRAPHMATCHSCANCOUNTDECL        ;
#define DGRAPHMATCHSCANCOUNTINIT        \
      probval = 0;                                /* Vertices will always be active */
#define DGRAPHMATCHSCANCOUNTSELECT      \
          if (vertgstend < vertlocnnd)  \
            edgefrenbr ++;              \
          else                          \
            edgeendnbr --;
#define DGRAPHMATCHSCANFINDSELECT       \
          ((vertgstend < vertlocnnd) && \
           (edgefrenbr -- == 0))
#include "dgraph_match_scan.c"
#undef DGRAPHMATCHSCANFINDSELECT
#undef DGRAPHMATCHSCANCOUNTSELECT
#undef DGRAPHMATCHSCANCOUNTINIT
#undef DGRAPHMATCHSCANCOUNTDECL
#undef DGRAPHMATCHSCANINIT
#undef DGRAPHMATCHSCANNAME

#define DGRAPHMATCHSCANNAME         dgraphMatchLy /* Local heavy edge matching */
#define DGRAPHMATCHSCANINIT                                                    \
  const Gnum * restrict const edloloctax = mateptr->c.finegrafptr->edloloctax; \
  if (edloloctax == NULL) {                                                    \
    dgraphMatchLc (mateptr);                                                   \
    return;                                                                    \
  }                                                                            \
  probmax = 0;                                    /* Vertices will always be active */
#define DGRAPHMATCHSCANCOUNTDECL                                               \
      Gnum                edlolocmax;
#define DGRAPHMATCHSCANCOUNTINIT                                               \
      edlolocmax = 0;                                                          \
      probval = 0;                                /* Vertices will always be active */
#define DGRAPHMATCHSCANCOUNTSELECT                                             \
          if (vertgstend < vertlocnnd) {                                       \
            Gnum                edlolocval;                                    \
            edlolocval = edloloctax[edgelocnum];                               \
            if (edlolocval > edlolocmax) {                                     \
              edlolocmax = edlolocval;                                         \
              edgefrenbr = 1;                                                  \
            }                                                                  \
            else if (edlolocval == edlolocmax)                                 \
              edgefrenbr ++;                                                   \
          }                                                                    \
          else                                                                 \
            edgeendnbr --;
#define DGRAPHMATCHSCANFINDSELECT                                              \
          ((vertgstend < vertlocnnd) &&                                        \
           (edloloctax[edgelocnum] == edlolocmax) &&                           \
           (edgefrenbr -- == 0))
#include "dgraph_match_scan.c"
#undef DGRAPHMATCHSCANFINDSELECT
#undef DGRAPHMATCHSCANCOUNTSELECT
#undef DGRAPHMATCHSCANCOUNTINIT
#undef DGRAPHMATCHSCANCOUNTDECL
#undef DGRAPHMATCHSCANINIT
#undef DGRAPHMATCHSCANNAME
