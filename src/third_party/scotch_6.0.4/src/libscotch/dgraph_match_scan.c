/* Copyright 2008,2012,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_match_scan.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines define a generic scanning   **/
/**                framework for distributed graph         **/
/**                matching routines.                      **/
/**                                                        **/
/**    DATES     : # Version 6.0  : from : 04 dec 2008     **/
/**                                 to   : 10 oct 2013     **/
/**                                                        **/
/************************************************************/

void
DGRAPHMATCHSCANNAME (
DgraphMatchData * restrict const  mateptr)
{
  int                 flagval;
  Gnum                vertlocnnd;
  Gnum                vertlocadj;
  Gnum                edgekptnbr;
  Gnum                queulocnbr;
  Gnum                matelocnbr;                 /* TRICK: Initial number of local mated vertices, from which are subtracted single multinodes */
  Gnum                multlocnbr;
  Gnum                probmax;

  const Dgraph * restrict const       grafptr    = mateptr->c.finegrafptr;
  const Gnum * restrict const         vertloctax = grafptr->vertloctax;
  const Gnum * restrict const         vendloctax = grafptr->vendloctax;
  const Gnum * restrict const         edgegsttax = grafptr->edgegsttax;
  Gnum * restrict const               queuloctab = mateptr->queuloctab;
  Gnum * restrict const               mategsttax = mateptr->mategsttax;
  DgraphCoarsenMulti * restrict const multloctab = mateptr->c.multloctab;
  DGRAPHMATCHSCANINIT

  flagval = mateptr->c.flagval;                   /* Get flag value */
  vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;
  vertlocnnd = grafptr->vertlocnnd;
  matelocnbr = mateptr->matelocnbr;
  multlocnbr = mateptr->c.multlocnbr;
  edgekptnbr = mateptr->c.edgekptnbr;

  if (matelocnbr == 0) {                          /* If initial pass or nothing useful done */
    Gnum                vertlocnnt;               /* End of useable local vertices          */
    Gnum                vertlocnum;

    memSet (mategsttax + grafptr->baseval, ~0, grafptr->vertlocnbr * sizeof (Gnum)); /* No local vertices matched to date: wipe all unsatisfied queries */

    queulocnbr = 0;                               /* Build queue from scratch */
    for (vertlocnum = grafptr->baseval, vertlocnnt = vertlocnnd;
         vertlocnum < vertlocnnt; vertlocnum ++) {
      Gnum                edgelocnum;
      Gnum                edgelocnnd;
      Gnum                edgeendnbr;
      Gnum                edgefrenbr;
      Gnum                probval;
      DGRAPHMATCHSCANCOUNTDECL

      if (mategsttax[vertlocnum] >= 0)            /* If vertex has been matched by one of the previous ones, skip it */
        continue;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (mategsttax[vertlocnum] < -1) {          /* Vertex must not be requesting yet */
        errorPrint ("dgraphMatchSc: internal error (1)");
        return;
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      DGRAPHMATCHSCANCOUNTINIT

      if (probval > probmax) {                    /* If vertex not active this turn */
        queuloctab[queulocnbr ++] = vertlocnum;   /* Enqueue it for next time       */
        continue;                                 /* Skip to next vertex            */
      }

      edgelocnum = vertloctax[vertlocnum];
      edgelocnnd = vendloctax[vertlocnum];
      if (((flagval & DGRAPHCOARSENNOMERGE) == 0) && /* If merging isolated vertices is allowed  */
          ((edgelocnnd - edgelocnum) == 0)) {     /* And if vertex is isolated                   */
        while (mategsttax[-- vertlocnnt] != ~0) ; /* Search for first matchable local "neighbor" */

        mategsttax[vertlocnum] = (vertlocnnt + vertlocadj); /* At worst we will stop at vertlocnum */
        mategsttax[vertlocnnt] = (vertlocnum + vertlocadj);
        multloctab[multlocnbr].vertglbnum[0] = (vertlocnum + vertlocadj);
        multloctab[multlocnbr].vertglbnum[1] = (vertlocnnt + vertlocadj);
        multlocnbr ++;                            /* One more coarse vertex created (two more local mates) */
        edgekptnbr += vendloctax[vertlocnnt] - vertloctax[vertlocnnt]; /* Add edges of other vertex only   */
        continue;
      }

      for (edgeendnbr = edgefrenbr = 0; edgelocnum < edgelocnnd; edgelocnum ++) { /* For all edges, count yet unmatched ends */
        Gnum                vertgstend;

        vertgstend = edgegsttax[edgelocnum];
        if (mategsttax[vertgstend] == -1) {       /* Count relevant free end vertices */
          DGRAPHMATCHSCANCOUNTSELECT
        }
        if (mategsttax[vertgstend] < 0)           /* Count not yet assigned end vertices */
          edgeendnbr ++;
      }
      if (edgeendnbr <= 0) {                      /* If vertex has no possible neighbor */
        mategsttax[vertlocnum] =                  /* Create single multinode            */
        multloctab[multlocnbr].vertglbnum[0] =
        multloctab[multlocnbr].vertglbnum[1] = vertlocnum + vertlocadj;
        multlocnbr ++;                            /* One more coarse vertex created  */
        matelocnbr --;                            /* TRICK: But with only one vertex */
        edgekptnbr += edgelocnnd - vertloctax[vertlocnum];
        continue;
      }

      if (edgefrenbr > 0) {                       /* If vertex has some free neighbor */
        Gnum                vertgstend;

        edgefrenbr = intRandVal (edgefrenbr);     /* Select one of them randomly */

        for (edgelocnum = vertloctax[vertlocnum]; ; edgelocnum ++) { /* Loop again on edges */
#ifdef SCOTCH_DEBUG_DGRAPH2
          if (edgelocnum >= edgelocnnd) {
            errorPrint ("dgraphMatchSc: internal error (2)");
            return;
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          vertgstend = edgegsttax[edgelocnum];
          if (mategsttax[vertgstend] == -1) {     /* If free end vertex found */
            if (DGRAPHMATCHSCANFINDSELECT)        /* If it is the one we want */
              break;                              /* Exit loop                */
          }
        }

        if (vertgstend >= vertlocnnd) {           /* If end vertex is a ghost vertex             */
          queuloctab[queulocnbr ++] = vertlocnum; /* Enqueue vertex for communication processing */
          mategsttax[vertlocnum] = -2 - edgelocnum; /* Local vertex index codes edge number      */
        }
        else {                                    /* Perform local matching */
          mategsttax[vertlocnum] = (vertgstend + vertlocadj);
          mategsttax[vertgstend] = (vertlocnum + vertlocadj);
          multloctab[multlocnbr].vertglbnum[0] = (vertlocnum + vertlocadj);
          multloctab[multlocnbr].vertglbnum[1] = (vertgstend + vertlocadj);
          multlocnbr ++;                          /* One more coarse vertex created (two more local mates) */
          edgekptnbr += (edgelocnnd - vertloctax[vertlocnum]) + (vendloctax[vertgstend] - vertloctax[vertgstend]) - 2; /* "-2" for collapsed arcs */
        }
      }
      else
        queuloctab[queulocnbr ++] = vertlocnum;   /* Enqueue vertex for next time */
    }
  }
  else {                                          /* Vertices to consider are already enqueued */
    Gnum                queulocnum;
    Gnum                queulocnew;

    for (queulocnum = queulocnew = 0, queulocnbr = mateptr->queulocnbr; queulocnum < queulocnbr; queulocnum ++) { /* For all vertices in queue */
      Gnum                vertlocnum;
      Gnum                mategstnum;

      vertlocnum = queuloctab[queulocnum];        /* Get current vertex */
      mategstnum = mategsttax[vertlocnum];
      if (mategstnum > -1)                        /* If already mated */
        continue;                                 /* Find another one */
      queuloctab[queulocnew ++] = vertlocnum;
      if (mategstnum < -1)
        mategsttax[vertlocnum] = -1;
    }
    queulocnbr = queulocnew;

    for (queulocnum = 0; queulocnum < queulocnbr; queulocnum ++) { /* For all vertices in queue */
      Gnum                vertlocnum;
      Gnum                edgelocnum;
      Gnum                edgelocnnd;
      Gnum                edgeendnbr;
      Gnum                edgefrenbr;
      Gnum                probval;
      DGRAPHMATCHSCANCOUNTDECL

      vertlocnum = queuloctab[queulocnum];        /* Get current vertex */
      if (mategsttax[vertlocnum] >= 0)            /* If already mated   */
        continue;                                 /* Find another one   */

#ifdef SCOTCH_DEBUG_DGRAPH2
      if (mategsttax[vertlocnum] < -1) {          /* Vertex must not be requesting yet */
        errorPrint ("dgraphMatchSc: internal error (3)");
        return;
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      DGRAPHMATCHSCANCOUNTINIT

      if (probval > probmax)                      /* If vertex not active this turn               */
        continue;                                 /* Keep vertex enqueued and skip to next vertex */

      edgelocnum = vertloctax[vertlocnum];
      edgelocnnd = vendloctax[vertlocnum];        /* No need to test for isolated vertices this turn */

      for (edgeendnbr = edgefrenbr = 0; edgelocnum < edgelocnnd; edgelocnum ++) { /* For all edges, count yet unmatched ends */
        Gnum                vertgstend;

        vertgstend = edgegsttax[edgelocnum];
        if (mategsttax[vertgstend] == -1) {       /* Count free end vertices */
          DGRAPHMATCHSCANCOUNTSELECT
        }
        if (mategsttax[vertgstend] < 0)           /* Count not yet assigned end vertices */
          edgeendnbr ++;
      }
      if (edgeendnbr <= 0) {                      /* If vertex has no possible neighbor */
        mategsttax[vertlocnum] =                  /* Create single multinode            */
        multloctab[multlocnbr].vertglbnum[0] =
        multloctab[multlocnbr].vertglbnum[1] = vertlocnum + vertlocadj;
        multlocnbr ++;                            /* One more coarse vertex created  */
        matelocnbr --;                            /* TRICK: But with only one vertex */
        edgekptnbr += edgelocnnd - vertloctax[vertlocnum];
        continue;
      }

      if (edgefrenbr > 0) {                       /* If vertex has some free neighbor */
        Gnum                vertgstend;

        edgefrenbr = intRandVal (edgefrenbr);     /* Select one of them randomly */

        for (edgelocnum = vertloctax[vertlocnum]; ; edgelocnum ++) { /* Loop again on edges */
#ifdef SCOTCH_DEBUG_DGRAPH2
          if (edgelocnum >= edgelocnnd) {
            errorPrint ("dgraphMatchSc: internal error (4)");
            return;
          }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
          vertgstend = edgegsttax[edgelocnum];
          if (mategsttax[vertgstend] == -1) {     /* If free end vertex found */
            if (DGRAPHMATCHSCANFINDSELECT)        /* If it is the one we want */
              break;                              /* Exit loop                */
          }
        }

        if (vertgstend >= vertlocnnd)             /* If end vertex is a ghost vertex        */
          mategsttax[vertlocnum] = -2 - edgelocnum; /* Local vertex index codes edge number */
        else {                                    /* Perform local matching */
          mategsttax[vertlocnum] = (vertgstend + vertlocadj);
          mategsttax[vertgstend] = (vertlocnum + vertlocadj);
          multloctab[multlocnbr].vertglbnum[0] = (vertlocnum + vertlocadj);
          multloctab[multlocnbr].vertglbnum[1] = (vertgstend + vertlocadj);
          multlocnbr ++;                          /* One more coarse vertex created (two more local mates) */
          edgekptnbr += (edgelocnnd - vertloctax[vertlocnum]) + (vendloctax[vertgstend] - vertloctax[vertgstend]) - 1;
        }
      }                                           /* Else vertex stays enqueued */
    }
  }

  mateptr->matelocnbr   = matelocnbr + 2 * (multlocnbr - mateptr->c.multlocnbr); /* TRICK: Two times new multinodes, minus single multinode adjustment */
  mateptr->queulocnbr   = queulocnbr;
  mateptr->c.multlocnbr = multlocnbr;
  mateptr->c.edgekptnbr = edgekptnbr;
}
