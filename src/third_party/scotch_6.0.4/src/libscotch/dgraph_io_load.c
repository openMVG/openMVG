/* Copyright 2007-2009,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_io_load.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the input/output routines for       **/
/**                distributed graphs.                     **/
/**                                                        **/
/**                # Version  5.0 : from : 28 apr 2007     **/
/**                                 to   : 24 mar 2008     **/
/**                # Version  5.1 : from : 23 jun 2008     **/
/**                                 to   : 27 jan 2009     **/
/**                # Version  6.0 : from : 25 aug 2012     **/
/**                                 to   : 18 nov 2012     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define DGRAPH_IO_LOAD

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph.h"
#include "dgraph_allreduce.h"
#include "dgraph_io_load.h"

/* This routine loads a distributed source
** graph from the given stream(s). Either
** one processor holds a non-NULL stream
** of a centralized graph, or all of them
** hold valid streams to either a centralized
** or a distributed graph.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (6, 3)
DGRAPHALLREDUCEMAXSUMOP (10, 2)

int
dgraphLoad (
Dgraph * restrict const     grafptr,              /* Not const since halo may update structure         */
FILE * const                stream,               /* One single centralized stream or distributed ones */
const Gnum                  baseval,              /* Base value (-1 means keep file base)              */
const DgraphFlag            flagval)              /* Graph loading flags                               */
{
  Gnum                reduloctab[12];
  Gnum                reduglbtab[12];
  Gnum                versval;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Barrier (grafptr->proccomm) != MPI_SUCCESS) { /* Synchronize for debugging */
    errorPrint ("dgraphLoad: communication error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  reduloctab[0] = baseval;                        /* Exchange baseval to check it is the same for all */
  reduloctab[1] = - baseval;
  reduloctab[2] = flagval;                        /* Exchange flagval to check it is the same for all */
  reduloctab[3] = - flagval;
  reduloctab[4] = 0;                              /* Set uneffective values for versval */
  reduloctab[5] = -2;
  reduloctab[6] =                                 /* Assume everything will be fine */
  reduloctab[7] =                                 /* Assume does not have a stream  */
  reduloctab[8] = 0;
  if (stream != NULL) {
    if (intLoad (stream, &versval) != 1) {        /* Read version number */
      errorPrint ("dgraphLoad: bad input (1)");
      versval       = 0;
      reduloctab[6] = 1;
    }
    else if ((versval != 0) && (versval != 2)) {  /* If not a graph format */
      errorPrint ("dgraphLoad: not a graph format");
      reduloctab[6] = 1;
    }
    reduloctab[4] = versval;
    reduloctab[5] = - versval;
    reduloctab[7] = 1;                            /* One more process involved in loading */
    reduloctab[8] = grafptr->proclocnum;
  }

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 6, 3, grafptr->proccomm) != 0) {
    errorPrint ("dgraphLoad: communication error (2)");
    return     (1);
  }

  if (reduglbtab[6] != 0)                         /* Return from previous errors */
    return (1);

  if ((reduglbtab[0] != - reduglbtab[1])) {
    errorPrint ("dgraphLoad: inconsistent base value");
    return     (1);
  }
  if ((reduglbtab[2] != - reduglbtab[3])) {
    errorPrint ("dgraphLoad: inconsistent flag value");
    return     (1);
  }
  if ((reduglbtab[7] != 0) &&
      (reduglbtab[4] != - reduglbtab[5])) {
    errorPrint ("dgraphLoad: inconsistent graph file version value");
    return     (1);
  }

  if (reduglbtab[4] == 2) {                       /* If distributed graph format             */
    if (reduglbtab[7] == grafptr->procglbnbr)     /* If as many input streams as processors  */
      return (dgraphLoadDist (grafptr, stream, baseval, flagval)); /* Read distributed graph */
  }
  else {                                          /* If centralized graph format */
    if (reduglbtab[7] == 1)                       /* If only one reader stream   */
      return (dgraphLoadCent (grafptr, stream, baseval, flagval, reduglbtab[8])); /* Distribute centralized graph from known root */
    else if (reduglbtab[7] == grafptr->procglbnbr)
      return (dgraphLoadMulti (grafptr, stream, baseval, flagval)); /* Read multi-centralized graph */
  }

  errorPrint ((reduglbtab[7] == 0)
              ? "dgraphLoad: no input stream provided"
              : "dgraphLoad: invalid number of input streams");
  return     (1);
}

/* This routine loads a centralized source
** graph from a single stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dgraphLoadCent (
Dgraph * restrict const     grafptr,              /* Distributed graph to load            */
FILE * const                stream,               /* One single centralized stream        */
Gnum                        baseval,              /* Base value (-1 means keep file base) */
const DgraphFlag            flagval,              /* Graph loading flags                  */
const int                   protnum)              /* Root process number                  */
{
  Gnum                vertglbnbr;
  Gnum                vertglbmax;
  Gnum                vertlocnbr;
  Gnum *              vertloctax;                 /* [norestrict:async] */
  Gnum *              vertlocptr;
  Gnum * restrict     vertredtax;
  Gnum                velolocnbr;
  Gnum                velolocsum;
  Gnum *              veloloctax;
  Gnum * restrict     veloredtax;
  Gnum                vlbllocnbr;
  Gnum *              vlblloctax;
  Gnum * restrict     vlblredtax;
  Gnum                edgelocnbr;
  Gnum *              edgeloctax;
  Gnum                edgeredmnd;
  Gnum * restrict     edgeredtax;
  Gnum *              edloloctax;
  Gnum * restrict     edloredtax;
  Gnum                degrglbmax;
  Gnum                baseadj;
  Gnum                reduglbtab[5];
  char                proptab[4];                 /* Property string array */
  int                 cheklocval;
  int                 chekglbval;
  int                 o;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (((stream != NULL) && (protnum != grafptr->proclocnum)) || /* Enforce single stream */
      ((stream == NULL) && (protnum == grafptr->proclocnum))) {
    errorPrint ("dgraphLoadCent: invalid parameter (1)");
    return     (1);
  }
  if ((protnum < 0) || (protnum >= grafptr->procglbnbr)) {
    errorPrint ("dgraphLoadCent: invalid parameter (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  reduglbtab[0] = 0;                              /* Assume everything will be fine */
  if (stream != NULL) {
    if ((intLoad (stream, &reduglbtab[1]) != 1) || /* Read rest of header */
        (intLoad (stream, &reduglbtab[2]) != 1) ||
        (intLoad (stream, &reduglbtab[3]) != 1) ||
        (intLoad (stream, &reduglbtab[4]) != 1) ||
        (reduglbtab[4] < 0)                     ||
        (reduglbtab[4] > 111)) {
      errorPrint ("dgraphLoadCent: bad input (1)");
      cheklocval = 1;
    }
  }

  if (MPI_Bcast (&reduglbtab[0], 5, GNUM_MPI, protnum, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphLoadCent: communication error (1)");
    return     (1);
  }
  if (reduglbtab[0] != 0)
    return (1);

  if (baseval == -1) {                            /* If keep file graph base     */
    baseval = reduglbtab[3];                      /* Set graph base as file base */
    baseadj = 0;                                  /* No base adjustment needed   */
  }
  else                                            /* If set graph base  */
    baseadj = baseval - reduglbtab[3];            /* Update base adjust */

  vertglbnbr = reduglbtab[1];
  vertglbmax = DATASIZE (vertglbnbr, grafptr->procglbnbr, 0);
  vertlocnbr = DATASIZE (vertglbnbr, grafptr->procglbnbr, grafptr->proclocnum);

  sprintf (proptab, "%3.3d", (int) reduglbtab[4]); /* Compute file properties */
  proptab[0] -= '0';                              /* Vertex labels flag       */
  proptab[1] -= '0';                              /* Edge weights flag        */
  proptab[2] -= '0';                              /* Vertex loads flag        */

  velolocnbr = ((proptab[2] != 0) && ((flagval & GRAPHIONOLOADVERT) == 0)) ? vertglbmax : 0;
  vlbllocnbr = (proptab[0] != 0) ? vertglbmax : 0;
  vlblloctax =
  veloloctax =
  vertloctax =
  edgeloctax =                                    /* Send arrays not allocated yet for root process */
  edgeredtax = NULL;                              /* No read edge array yet                         */
  cheklocval = 0;
  if ((vertlocptr = memAlloc ((vertglbmax + 2 + velolocnbr + vlbllocnbr) * sizeof (Gnum))) == NULL) { /* TRICK: "+2" for space for velolocsum */
    errorPrint ("dgraphLoadCent: out of memory (1)");
    cheklocval = 1;
  }
  else {
    vertloctax  =
    vertlocptr -= baseval;
    vertlocptr += vertglbmax + 2;                 /* TRICK: "+2" for space for velolocsum */
    if (proptab[2] != 0) {
      veloloctax  = vertlocptr;
      vertlocptr += vertglbmax;
    }
    if (proptab[0] != 0) {
      vlblloctax = vertlocptr;
      baseadj    = 0;                             /* No vertex adjustments if vertex labels */
    }

    if (stream != NULL) {                         /* Allocate read edge array */
      Gnum                edgeredmax;
      Gnum                edloredmax;

      edgeredmax  = reduglbtab[2] / grafptr->procglbnbr + 1;
      edgeredmax += (edgeredmax >> 2) + 4;        /* Add 25% more space for edges than average */
      edloredmax  = ((proptab[1] != 0) && ((flagval & GRAPHIONOLOADEDGE) == 0)) ? edgeredmax : 0;

      if ((edgeredtax = memAlloc ((edgeredmax + edloredmax) * sizeof (Gnum))) == NULL) {
        errorPrint ("dgraphLoadCent: out of memory (2)");
        cheklocval = 1;
      }
      else {
        edgeredtax -= baseval;
        edloredtax = (edloredmax != 0) ? (edgeredtax + edgeredmax) : NULL;

        vertredtax = vertloctax;                  /* Prepare read vertex arrays, which will never move */
        veloredtax = veloloctax;
        vlblredtax = vlblloctax;
      }
      edgeredmnd = edgeredmax + baseval;
    }
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphLoadCent: communication error (2)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    if (edgeredtax != NULL)
      memFree (edgeredtax + baseval);
    if (vertloctax != NULL)
      memFree (vertloctax + baseval);
    return  (1);
  }

  degrglbmax = 0;                                 /* No maximum degree yet */

  if (stream != NULL) {
    Gnum                procnum;

    for (procnum = 0; procnum < grafptr->procglbnbr; procnum ++) {
      Gnum                vertrednnd;
      Gnum                vertrednum;
      Gnum                edgerednum;
      Gnum                veloredsum;

      for (vertrednum = edgerednum = baseval, veloredsum = 0,
           vertrednnd = DATASIZE (vertglbnbr, grafptr->procglbnbr, procnum) + baseval;
           vertrednum < vertrednnd; vertrednum ++) {
        Gnum                degrredval;

        if (vlblredtax != NULL) {                 /* If must read label            */
          Gnum                vlblredval;         /* Vertex label value to be read */

          if (intLoad (stream, &vlblredval) != 1) {  /* Read label data */
            errorPrint ("dgraphLoadCent: bad input (2)");
            cheklocval = 1;
            break;
          }
          vlblredtax[vertrednum] = vlblredval;
        }
        if (proptab[2] != 0) {                    /* If must read vertex load */
          Gnum                veloredval;

          if (intLoad (stream, &veloredval) != 1) { /* Read vertex load data */
            errorPrint ("dgraphLoadCent: bad input (3)");
            cheklocval = 1;
            break;
          }
          if (veloredtax != NULL)
            veloredsum            +=
            veloredtax[vertrednum] = veloredval;
        }
        if (intLoad (stream, &degrredval) != 1) { /* Read vertex degree */
          errorPrint ("dgraphLoadCent: bad input (4)");
          cheklocval = 1;
          break;
        }
        if (degrglbmax < degrredval)              /* Set maximum degree */
          degrglbmax = degrredval;

        vertredtax[vertrednum] = edgerednum;      /* Set index in edge array */
        degrredval += edgerednum;
        if (degrredval > edgeredmnd) {            /* Check if edge array overflows */
          Gnum                edgeredmax;
          Gnum                edgenewmax;
          Gnum                edgenewsiz;
          Gnum * restrict     edgenewtab;

          edgenewmax =
          edgeredmax = edgeredmnd - baseval;
          do                                      /* Increase edge array size by 25 % */
            edgenewmax += (edgenewmax >> 2) + 4;
          while (edgenewmax < (degrredval - baseval));
          edgenewsiz = (edloredtax != NULL) ? (2 * edgenewmax) : edgenewmax;
          if ((edgenewtab = memRealloc (edgeredtax + baseval, edgenewsiz * sizeof (Gnum))) == NULL) {
            errorPrint ("dgraphLoadCent: out of memory (3)");
            cheklocval = 1;
            break;
          }
          edgeredtax = edgenewtab - baseval;
          edgeredmnd = edgenewmax + baseval;
          if (edloredtax != NULL) {               /* Move edge load array if present */
            memMov (edgenewtab + edgenewmax, edgenewtab + edgeredmax, (edgerednum - baseval) * sizeof (Gnum));
            edloredtax = edgeredtax + edgenewmax;
          }
        }

        for ( ; edgerednum < degrredval; edgerednum ++) {
          Gnum                edgeredval;         /* Value where to read edge end */

          if (proptab[1] != 0) {                  /* If must read edge load        */
            Gnum                edloredval;       /* Value where to read edge load */

            if (intLoad (stream, &edloredval) != 1) { /* Read edge load data */
              errorPrint ("dgraphLoadCent: bad input (5)");
              cheklocval = 1;
              break;
            }
            if (edloredtax != NULL)
              edloredtax[edgerednum] = edloredval;
          }
          if (intLoad (stream, &edgeredval) != 1) { /* Read edge data */
            errorPrint ("dgraphLoadCent: bad input (6)");
            cheklocval = 1;
            break;
          }
          edgeredtax[edgerednum] = edgeredval + baseadj;
        }
        if (cheklocval != 0)
          break;
      }
      vertredtax[vertrednum ++] = edgerednum;     /* Set end of edge array */

      if (cheklocval == 0) {
        if (procnum != grafptr->proclocnum) {     /* If arrays have to be sent */
          MPI_Request         requtab[5];
          MPI_Status          stattab[5];
          int                 requnbr;

          vertredtax[baseval]    = edgerednum - baseval; /* First slot is number of edges */
          vertredtax[vertrednum] = (veloredtax != NULL) ? veloredsum : (vertrednnd - baseval); /* Add vertex load sum to send vertex array */
          if (MPI_Isend (vertredtax + baseval, vertrednnd - baseval + 2, /* TRICK: "+2" and not "+1" because of space for velolocsum       */
                         GNUM_MPI, procnum, TAGVERTLOCTAB, grafptr->proccomm, &requtab[0]) != MPI_SUCCESS) {
            errorPrint ("dgraphLoadCent: communication error (5)");
            return     (1);                       /* Dirty exit as we can do nothing */
          }
          requnbr = 1;
          if (veloredtax != NULL) {
            if (MPI_Isend (veloredtax + baseval, vertrednnd - baseval, GNUM_MPI,
                           procnum, TAGVELOLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
              errorPrint ("dgraphLoadCent: communication error (6)");
              return     (1);
            }
          }
          if (vlblredtax != NULL) {
            if (MPI_Isend (vlblredtax + baseval, vertrednnd - baseval, GNUM_MPI,
                           procnum, TAGVLBLLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
              errorPrint ("dgraphLoadCent: communication error (7)");
              return     (1);
            }
          }
          if (MPI_Recv (&reduglbtab[0], 0, MPI_INT, procnum, MPI_ANY_TAG, grafptr->proccomm, &stattab[0]) != MPI_SUCCESS) {
            errorPrint ("dgraphLoadCent: communication error (8)");
            return     (1);
          }
          if (stattab[0].MPI_TAG != TAGOK)        /* If receiver could not allocate memory for edge arrays */
            cheklocval = 1;
          else {
            if (MPI_Isend (edgeredtax + baseval, edgerednum - baseval, GNUM_MPI,
                           procnum, TAGEDGELOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
              errorPrint ("dgraphLoadCent: communication error (9)");
              return     (1);
            }
            if (edloredtax != NULL) {
              if (MPI_Isend (edloredtax + baseval, edgerednum - baseval, GNUM_MPI,
                             procnum, TAGEDLOLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
                errorPrint ("dgraphLoadCent: communication error (10)");
                return     (1);
              }
            }
            MPI_Waitall (requnbr, &requtab[0], &stattab[0]);
          }
        }
        else {                                    /* Arrays are local */
          velolocsum = (veloredtax != NULL) ? veloredsum : vertlocnbr; /* Save accumulated values as local data */
          edgelocnbr = edgerednum - baseval;
          if (edgeredmnd - edgerednum > 10000) {  /* If array can be compacted */
            if (edloredtax != NULL) {
              memMov (edgeredtax + edgerednum, edloredtax + baseval, edgelocnbr * sizeof (Gnum));
              edgeredtax  = memRealloc (edgeredtax + baseval, edgelocnbr * 2 * sizeof (Gnum));
              edgeredtax -= baseval;
              edloredtax  = edgeredtax + edgelocnbr;
            }
            else {
              edgeredtax  = memRealloc (edgeredtax + baseval, edgelocnbr * sizeof (Gnum));
              edgeredtax -= baseval;
            }
          }
          edgeloctax = edgeredtax;                /* Keep read edge array as local edge array */
          edloloctax = edloredtax;

          if (grafptr->proclocnum == (grafptr->procglbnbr - 1)) { /* If root process is last process */
            vertredtax =                          /* No need to reallocate read arrays               */
            edgeredtax = NULL;
            break;                                /* And we can exit now */
          }

          if ((vertlocptr = memAlloc ((vertglbmax + 2 + velolocnbr + vlbllocnbr) * sizeof (Gnum))) == NULL) { /* TRICK: "+2" for space for velolocsum */
            errorPrint ("dgraphLoadCent: out of memory (4)");
            cheklocval = 1;
          }
          else {
            Gnum                edgeredmax;
            Gnum                edloredmax;

            vertredtax  =
            vertlocptr -= baseval;
            vertlocptr += vertglbmax + 2;         /* TRICK: "+2" for space for velolocsum */
            if (veloredtax != NULL) {
              veloredtax  = vertlocptr;
              vertlocptr += vertglbmax;
            }
            if (vlblredtax != NULL)
              vlblredtax = vertlocptr;

            edgeredmax = edgeredmnd - baseval;
            edloredmax = (edloloctax != NULL) ? edgeredmax : 0;
            if ((edgeredtax = memAlloc ((edgeredmax + edloredmax) * sizeof (Gnum))) == NULL) {
              errorPrint ("dgraphLoadCent: out of memory (5)");
              cheklocval = 1;
            }
            else {
              edgeredtax -= baseval;
              if (edloredtax != NULL)
                edloredtax = edgeredtax + edgeredmax;
            }
          }
        }
      }

      if (cheklocval != 0) {                      /* If error encountered                           */
        for ( ; procnum < grafptr->procglbnbr; procnum ++) { /* Send abortion messages              */
          if (procnum != grafptr->proclocnum) {   /* Abortion messages complete vertloctab receives */
            if (MPI_Send (vertredtax + baseval, 0, GNUM_MPI, procnum, TAGVERTLOCTAB, grafptr->proccomm) != MPI_SUCCESS)
              errorPrint ("dgraphLoadCent: communication error (11)");
          }
        }
        break;
      }
    }

    if (vertredtax != NULL) {                     /* Free reader arrays if reallocated                   */
      if (vertredtax != vertloctax)               /* If equal, vertloctax will be deallocated afterwards */
        memFree (vertredtax + baseval);
      memFree (edgeredtax + baseval);
    }
  }
  else {                                          /* Process is not reader */
    MPI_Request         requtab[5];
    MPI_Status          stattab[5];
    int                 requnbr;
    int                 vertrcvnbr;               /* int because of the MPI API */

    if (MPI_Irecv (vertloctax + baseval, vertlocnbr + 2, GNUM_MPI, /* TRICK: "+2" and not "+1" because of velolocsum                      */
                   protnum, TAGVERTLOCTAB, grafptr->proccomm, &requtab[2]) != MPI_SUCCESS) { /* requtab[2] is first surely available slot */
      errorPrint ("dgraphLoadCent: communication error (10)");
      return     (1);                             /* Dirty exit as we can do nothing */
    }
    requnbr = 0;
    if (veloloctax != NULL) {
      if (MPI_Irecv (veloloctax + baseval, vertlocnbr, GNUM_MPI,
                     protnum, TAGVELOLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("dgraphLoadCent: communication error (11)");
        return     (1);
      }
    }
    if (vlblloctax != NULL) {
      if (MPI_Irecv (vlblloctax + baseval, vertlocnbr, GNUM_MPI,
                     protnum, TAGVLBLLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
        errorPrint ("dgraphLoadCent: communication error (12)");
        return     (1);
      }
    }

    MPI_Wait (&requtab[2], &stattab[2]);          /* Wait until vertloctab received   */
    MPI_Get_count (&stattab[2], GNUM_MPI, &vertrcvnbr); /* Get size of received array */
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (((Gnum) vertrcvnbr != 0) &&               /* vertrcvnbr == 0 in the case of abortion message */
        ((Gnum) vertrcvnbr != (vertlocnbr + 2))) { /* TRICK: "+2" and not "+1" because of velolocsum */
      errorPrint ("dgraphLoadCent: invalid vertex array size");
      vertrcvnbr = 0;
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    if (vertrcvnbr == 0) {                        /* Empty message means abortion wanted */
      while (requnbr > 0) {                       /* Cancel all pending requests         */
        requnbr --;
        MPI_Cancel       (&requtab[requnbr]);
        MPI_Request_free (&requtab[requnbr]);
      }                                           /* No more pending requests */
    }
    else {
      Gnum                edlolocnbr;

      edgelocnbr = vertloctax[baseval];           /* edgelocnbr is first cell     */
      vertloctax[baseval] = baseval;              /* First cell is always baseval */
      velolocsum = vertloctax[vertlocnbr + baseval + 1]; /* TRICK: get velolocsum */
      edlolocnbr = ((proptab[1] != 0) && ((flagval & GRAPHIONOLOADEDGE) == 0)) ? edgelocnbr : 0;
      edloloctax = NULL;                          /* Assume no edge load array */

      if ((edgeloctax = memAlloc ((edgelocnbr + edlolocnbr) * sizeof (Gnum))) == NULL) {
        errorPrint ("dgraphLoadCent: out of memory (6)");
        MPI_Send   (&cheklocval, 0, MPI_INT, protnum, TAGBAD, grafptr->proccomm); /* Memory could not be allocated */
        cheklocval = 1;
      }
      else {
        if (MPI_Irecv (edgeloctax, edgelocnbr, GNUM_MPI, protnum, TAGEDGELOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
          errorPrint ("dgraphLoadCent: communication error (13)");
          return     (1);
        }
        if (edlolocnbr != 0) {
          edloloctax = edgeloctax + edgelocnbr;
          if (MPI_Irecv (edloloctax, edgelocnbr, GNUM_MPI, protnum, TAGEDLOLOCTAB, grafptr->proccomm, &requtab[requnbr ++]) != MPI_SUCCESS) {
            errorPrint ("dgraphLoadCent: communication error (14)");
            return     (1);
          }
          edloloctax -= baseval;
        }
        edgeloctax -= baseval;
        MPI_Isend (&cheklocval, 0, MPI_INT, protnum, TAGOK, grafptr->proccomm, &requtab[requnbr ++]); /* Send ready to receive */
      }
    }

    MPI_Waitall (requnbr, &requtab[0], &stattab[0]); /* Wait until all pending communications completed and all arrays received */
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphLoadCent: communication error (15)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    if (edgeloctax != NULL)
      memFree (edgeloctax + baseval);
    memFree (vertloctax + baseval);
    return  (1);
  }

  o = dgraphBuild2 (grafptr, baseval,             /* Build distributed graph */
                    vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, veloloctax, velolocsum, NULL, vlblloctax,
                    edgelocnbr, edgelocnbr, edgeloctax, NULL, edloloctax, degrglbmax); /* Non-readers will have degrglbmax set to 0 */
  grafptr->flagval |= DGRAPHFREETABS | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;

  return (o);
}

/* This routine loads a distributed source
** graph from a distributed source graph
** file spread across all of the streams.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dgraphLoadDist (
Dgraph * restrict const     grafptr,              /* Distributed graph to load            */
FILE * const                stream,               /* One single centralized stream        */
Gnum                        baseval,              /* Base value (-1 means keep file base) */
const DgraphFlag            flagval)              /* Graph loading flags                  */
{
  Gnum                proclocnum;
  Gnum                vertlocnbr;
  Gnum                vertlocnnd;
  Gnum                vertlocnum;
  Gnum * restrict     vertloctax;
  Gnum *              vertlocptr;
  Gnum                velolocnbr;
  Gnum                velolocsum;
  Gnum * restrict     veloloctax;
  Gnum                vlbllocnbr;
  Gnum * restrict     vlblloctax;
  Gnum                edgelocnbr;
  Gnum                edgelocnnd;
  Gnum                edgelocnum;
  Gnum * restrict     edgeloctax;
  Gnum * restrict     edloloctax;
  Gnum                degrlocmax;
  Gnum                baseadj;
  Gnum                reduloctab[12];
  Gnum                reduglbtab[12];
  char                proptab[4];                 /* Property string array */
  int                 cheklocval;
  int                 chekglbval;
  int                 o;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (stream == NULL) {
    errorPrint ("dgraphLoadDist: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  reduloctab[0] = 0;                              /* Assume everything will be fine */
  o  = intLoad (stream, &reduloctab[1]);          /* Read rest of header            */
  o += intLoad (stream, &proclocnum);
  o += intLoad (stream, &reduloctab[3]);
  o += intLoad (stream, &reduloctab[5]);
  o += intLoad (stream, &reduloctab[10]);
  o += intLoad (stream, &reduloctab[11]);
  o += intLoad (stream, &reduloctab[7]);
  o += intLoad (stream, &reduloctab[9]);
  if ((o != 8)            ||
      (reduloctab[9] < 0) ||
      (reduloctab[9] > 111)) {
    errorPrint ("dgraphLoadDist: bad input (1)");
    reduloctab[0] = 2;                            /* Immediate abort has maximum value so as to be propagated by MAX reduce */
  }
  reduloctab[2]  = - reduloctab[1];
  reduloctab[4]  = - reduloctab[3];
  reduloctab[6]  = - reduloctab[5];
  reduloctab[8]  = - reduloctab[7];
  if ((int) proclocnum != grafptr->proclocnum)    /* If fragment is not read by proper process */
    reduloctab[0] |= 1;
  if ((int) reduloctab[1] != grafptr->procglbnbr) {
    errorPrint ("dgraphLoadDist: wrong number of processors to read distributed graph");
    reduloctab[0] = 2;
  }

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 10, 2, grafptr->proccomm) != 0) {
    errorPrint ("dgraphLoadDist: communication error (1)");
    reduglbtab[0] = 2;
  }
  if (reduglbtab[0] >= 2)                         /* If has to abort immediately */
    return (1);

  if ((reduglbtab[2] != - reduglbtab[1]) ||
      (reduglbtab[4] != - reduglbtab[3]) ||
      (reduglbtab[6] != - reduglbtab[5]) ||
      (reduglbtab[8] != - reduglbtab[7])) {
    errorPrint ("dgraphLoadDist: inconsistent distributed graph headers");
    return     (1);
  }
  if (reduloctab[0] == 1)
    errorPrint ("dgraphLoadDist: distributed graph file not read by proper process");
  if (reduglbtab[0] != 0)                         /* If cannot go on anyway */
    return (1);

  if ((reduglbtab[10] != reduloctab[3]) ||
      (reduglbtab[11] != reduloctab[5]))
    errorPrint ("dgraphLoadDist: bad input (2)");
  if ((reduglbtab[10] != reduglbtab[3]) ||
      (reduglbtab[11] != reduglbtab[5]))
    return (1);

  if (baseval == -1) {                            /* If keep file graph base     */
    baseval = reduglbtab[7];                      /* Set graph base as file base */
    baseadj = 0;                                  /* No base adjustment needed   */
  }
  else                                            /* If set graph base  */
    baseadj = baseval - reduglbtab[7];            /* Update base adjust */

  vertlocnbr = reduloctab[10];
  edgelocnbr = reduloctab[11];

  sprintf (proptab, "%3.3d", (int) reduglbtab[9]); /* Compute file properties */
  proptab[0] -= '0';                              /* Vertex labels flag       */
  proptab[1] -= '0';                              /* Edge weights flag        */
  proptab[2] -= '0';                              /* Vertex loads flag        */

  velolocnbr = ((proptab[2] != 0) && ((flagval & GRAPHIONOLOADVERT) == 0)) ? vertlocnbr : 0;
  vlbllocnbr = (proptab[0] != 0) ? vertlocnbr : 0;
  vlblloctax =
  veloloctax =
  vertloctax =
  edgeloctax = NULL;
  cheklocval = 0;
  if ((vertlocptr = memAlloc ((vertlocnbr + 1 + velolocnbr + vlbllocnbr) * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphLoadDist: out of memory (1)");
    cheklocval = 1;
  }
  else {
    Gnum                edlolocnbr;

    vertloctax  =
    vertlocptr -= baseval;
    vertlocptr += vertlocnbr + 1;
    if ((proptab[2] != 0) && ((flagval & GRAPHIONOLOADVERT) == 0)) {
      veloloctax  = vertlocptr;
      vertlocptr += vertlocnbr;
    }
    if (proptab[0] != 0) {
      vlblloctax = vertlocptr;
      baseadj    = 0;                             /* No vertex adjustments if vertex labels */
    }

    edlolocnbr = ((proptab[1] != 0) && ((flagval & GRAPHIONOLOADEDGE) == 0)) ? edgelocnbr : 0;
    if ((edgeloctax = memAlloc ((edgelocnbr + edlolocnbr) * sizeof (Gnum))) == NULL) {
      errorPrint ("dgraphLoadDist: out of memory (2)");
      cheklocval = 1;
    }
    else {
      edgeloctax -= baseval;
      edloloctax  = ((proptab[1] != 0) && ((flagval & GRAPHIONOLOADEDGE) == 0)) ? (edgeloctax + edgelocnbr) : NULL;
    }
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphLoadDist: communication error (2)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    if (edgeloctax != NULL)
      memFree (edgeloctax + baseval);
    if (vertloctax != NULL)
      memFree (vertloctax + baseval);
    return  (1);
  }

  degrlocmax = 0;                                 /* No maximum degree yet */
  velolocsum = (veloloctax != NULL) ? 0 : vertlocnbr;
  for (vertlocnum = edgelocnum = baseval, vertlocnnd = vertlocnbr + baseval, edgelocnnd = edgelocnbr + baseval;
       vertlocnum < vertlocnnd; vertlocnum ++) {
    Gnum                degrlocval;

    if (vlblloctax != NULL) {                     /* If must read label            */
      Gnum                vlbllocval;             /* Vertex label value to be read */

      if (intLoad (stream, &vlbllocval) != 1) {  /* Read label data */
        errorPrint ("dgraphLoadDist: bad input (2)");
        cheklocval = 1;
        break;
      }
      vlblloctax[vertlocnum] = vlbllocval;
    }
    if (proptab[2] != 0) {                    /* If must read vertex load */
      Gnum                velolocval;

      if (intLoad (stream, &velolocval) != 1) { /* Read vertex load data */
        errorPrint ("dgraphLoadDist: bad input (3)");
        cheklocval = 1;
        break;
      }
      if (veloloctax != NULL)
        velolocsum            +=
        veloloctax[vertlocnum] = velolocval;
    }
    if (intLoad (stream, &degrlocval) != 1) {     /* Read vertex degree */
      errorPrint ("dgraphLoadDist: bad input (4)");
      cheklocval = 1;
      break;
    }
    if (degrlocmax < degrlocval)                  /* Set maximum degree */
      degrlocmax = degrlocval;

    vertloctax[vertlocnum] = edgelocnum;          /* Set index in edge array */
    degrlocval += edgelocnum;
    if (degrlocval > edgelocnnd) {                /* Check if edge array overflows */
      errorPrint ("dgraphLoadDist: invalid arc count (1)");
      cheklocval = 1;
      break;
    }

    for ( ; edgelocnum < degrlocval; edgelocnum ++) {
      Gnum                edgelocval;         /* Value where to read edge end */

      if (proptab[1] != 0) {                  /* If must read edge load        */
        Gnum                edlolocval;       /* Value where to read edge load */

        if (intLoad (stream, &edlolocval) != 1) { /* Read edge load data */
          errorPrint ("dgraphLoadDist: bad input (5)");
          cheklocval = 1;
          break;
        }
        if (edloloctax != NULL)
          edloloctax[edgelocnum] = edlolocval;
      }
      if (intLoad (stream, &edgelocval) != 1) { /* Read edge data */
        errorPrint ("dgraphLoadDist: bad input (6)");
        cheklocval = 1;
        break;
      }
      edgeloctax[edgelocnum] = edgelocval + baseadj;
    }
    if (cheklocval != 0)
      break;
  }
  vertloctax[vertlocnum] = edgelocnum;            /* Set end of edge array             */
  if (edgelocnum != edgelocnnd) {                 /* Check if number of edges is valid */
    errorPrint ("dgraphLoadDist: invalid arc count (2)");
    cheklocval = 1;
  }

  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphLoadDist: communication error (17)");
    chekglbval = 1;
  }
  if (chekglbval != 0) {
    memFree (edgeloctax + baseval);
    memFree (vertloctax + baseval);
    return  (1);
  }

  o = dgraphBuild2 (grafptr, baseval,             /* Build distributed graph */
                    vertlocnbr, vertlocnbr, vertloctax, vertloctax + 1, veloloctax, velolocsum, NULL, vlblloctax,
                    edgelocnbr, edgelocnbr, edgeloctax, NULL, edloloctax, degrlocmax);
  grafptr->flagval |= DGRAPHFREETABS | DGRAPHVERTGROUP | DGRAPHEDGEGROUP;

  return (o);
}

/* This routine loads a distributed source
** graph from a centralized source graph
** file replicated on all of the streams.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dgraphLoadMulti (
Dgraph * restrict const     grafptr,              /* Distributed graph to load            */
FILE * const                stream,               /* One single centralized stream        */
Gnum                        baseval,              /* Base value (-1 means keep file base) */
const DgraphFlag            flagval)              /* Graph loading flags                  */
{
  errorPrint ("dgraphLoadMulti: not implemented");
  return     (1);
}
