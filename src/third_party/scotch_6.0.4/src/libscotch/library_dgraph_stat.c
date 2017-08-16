/* Copyright 2007,2008,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_stat.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the source   **/
/**                graph handling routines of the          **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 jun 2007     **/
/**                                 to     03 apr 2008     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to     29 oct 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "library_dgraph_stat.h"
#include "ptscotch.h"

/*
**  The static variables.
*/

static int                  dgraphstatblentab[2] = { 7, 3 };
static MPI_Datatype         dgraphstattypetab[2] = { GNUM_MPI, MPI_DOUBLE };

/************************************/
/*                                  */
/* These routines are the C API for */
/* the graph handling routines.     */
/*                                  */
/************************************/

/* This routine is the reduction-loc operator which
** returns in inout[2] the rank of the process which
** holds the best partition.
** It returns:
** - void  : in all cases.
*/

static
void
dgraphStatReduceAll (
const DgraphStatData * const  in,                 /* First operand                              */
DgraphStatData * const        inout,              /* Second and output operand                  */
const int * const             len,                /* Number of instances; should be 1, not used */
const MPI_Datatype * const    typedat)            /* MPI datatype; not used                     */
{
  if (inout->velomin > in->velomin)
    inout->velomin = in->velomin;
  if (inout->velomax < in->velomax)
    inout->velomax = in->velomax;
  if (inout->degrmin > in->degrmin)
    inout->degrmin = in->degrmin;
  if (inout->degrmax < in->degrmax)
    inout->degrmax = in->degrmax;
  if (inout->edlomin > in->edlomin)
    inout->edlomin = in->edlomin;
  if (inout->edlomax < in->edlomax)
    inout->edlomax = in->edlomax;
  inout->edlosum += in->edlosum;
  inout->velodlt += in->velodlt;
  inout->degrdlt += in->degrdlt;
  inout->edlodlt += in->edlodlt;
}

/*+ This routine computes statistics
*** on the given distributed graph.
*** It returns:
*** - VOID  : in all cases.
+*/

int
SCOTCH_dgraphStat (
const SCOTCH_Dgraph * const grafptr,
SCOTCH_Num * const          velominptr,
SCOTCH_Num * const          velomaxptr,
SCOTCH_Num * const          velosumptr,
double *                    veloavgptr,
double *                    velodltptr,
SCOTCH_Num * const          degrminptr,
SCOTCH_Num * const          degrmaxptr,
double *                    degravgptr,
double *                    degrdltptr,
SCOTCH_Num * const          edlominptr,
SCOTCH_Num * const          edlomaxptr,
SCOTCH_Num * const          edlosumptr,
double *                    edloavgptr,
double *                    edlodltptr)
{
  const Dgraph *      srcgrafptr;
  DgraphStatData      srcgstadat;
  DgraphStatData      srclstadat;
  MPI_Datatype        srctypedat;
  MPI_Aint            srcdisptab[2];
  MPI_Op              srcoperdat;
  Gnum                vertlocnum;
  double              veloglbavg;
  double              velolocdlt;
  Gnum                degrlocmin;
  Gnum                degrlocmax;
  double              degrglbavg;
  double              degrlocdlt;
  Gnum                edloglbsum;
  double              edloglbavg;
  double              edlolocdlt;
  int                 o;

  srcgrafptr = (Dgraph *) grafptr;

  velolocdlt = 0.0L;
  if (srcgrafptr->vertglbnbr > 0) {
    if (srcgrafptr->veloloctax != NULL) {         /* If graph has vertex loads */
      const Gnum * restrict veloloctax;
      Gnum                  velolocmin;
      Gnum                  velolocmax;

      veloloctax = srcgrafptr->veloloctax;
      velolocmin = GNUMMAX;
      velolocmax = 0;
      veloglbavg = (double) srcgrafptr->veloglbsum / (double) srcgrafptr->vertglbnbr;

      for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) {
        Gnum                velolocval;

        velolocval = veloloctax[vertlocnum];
        if (velolocval < velolocmin)
          velolocmin = velolocval;
        if (velolocval > velolocmax)
          velolocmax = velolocval;
        velolocdlt += fabs ((double) velolocval - veloglbavg);
      }

      srclstadat.velomin = velolocmin;
      srclstadat.velomax = velolocmax;
    }
    else {
      srclstadat.velomin =
      srclstadat.velomax = 1;
      veloglbavg         = 1.0L;
    }
  }
  else {
    srclstadat.velomin =
    srclstadat.velomax = 0;
    veloglbavg         = 0.0L;
  }
  srclstadat.velodlt = velolocdlt;

  degrlocmax = 0;
  degrlocdlt = 0.0L;
  if (srcgrafptr->vertglbnbr > 0) {
    degrlocmin = GNUMMAX;
    degrglbavg = (double) srcgrafptr->edgeglbnbr / (double) srcgrafptr->vertglbnbr;
    for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) {
      Gnum                degrlocval;

      degrlocval = srcgrafptr->vendloctax[vertlocnum] - srcgrafptr->vertloctax[vertlocnum]; /* Get vertex degree */
      if (degrlocval < degrlocmin)
        degrlocmin = degrlocval;
      if (degrlocval > degrlocmax)
        degrlocmax = degrlocval;
      degrlocdlt += fabs ((double) degrlocval - degrglbavg);
    }
  }
  else {
    degrlocmin = 0;
    degrglbavg = 0.0L;
  }
  srclstadat.degrmin = degrlocmin;
  srclstadat.degrmax = degrlocmax;
  srclstadat.degrdlt = degrlocdlt;

  edlolocdlt = 0.0L;
  if (srcgrafptr->edgeglbnbr > 0) {
    if (srcgrafptr->edloloctax != NULL) {         /* If graph has edge loads */
      Gnum                edlolocmin;
      Gnum                edlolocmax;
      Gnum                edlolocsum;

      edlolocmin = GNUMMAX;
      edlolocmax = 0;
      edlolocsum = 0;

      for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) {
        Gnum                edgelocnum;

        for (edgelocnum = srcgrafptr->vertloctax[vertlocnum];
             edgelocnum < srcgrafptr->vendloctax[vertlocnum]; edgelocnum ++) {
          Gnum                edlolocval;

          edlolocval  = srcgrafptr->edloloctax[edgelocnum];
          edlolocsum += edlolocval;
          if (edlolocval < edlolocmin)            /* Account for edge load */
            edlolocmin = edlolocval;
          if (edlolocval > edlolocmax)
            edlolocmax = edlolocval;
        }
      }

      if (MPI_Allreduce (&edlolocsum, &edloglbsum, 1, GNUM_MPI, MPI_SUM, srcgrafptr->proccomm) != MPI_SUCCESS) {
        errorPrint ("SCOTCH_dgraphStat: communication error (1)");
        return     (1);
      }
      edloglbavg = (double) edloglbsum / (double) (2 * srcgrafptr->edgeglbnbr);

      for (vertlocnum = srcgrafptr->baseval; vertlocnum < srcgrafptr->vertlocnnd; vertlocnum ++) {
        Gnum                edgelocnum;

        for (edgelocnum = srcgrafptr->vertloctax[vertlocnum];
             edgelocnum < srcgrafptr->vendloctax[vertlocnum]; edgelocnum ++)
          edlolocdlt += fabs ((double) srcgrafptr->edloloctax[edgelocnum] - edloglbavg);
      }
    }
    else {
      srclstadat.edlomin =
      srclstadat.edlomax = 1;
      edloglbsum         = srcgrafptr->edgeglbnbr / 2;
      edloglbavg         = 1.0L;
    }
  }
  else {
    srclstadat.edlomin =
    srclstadat.edlomax = 0;
    edloglbsum         = 0;
    edloglbavg         = 0.0L;
  }
  srclstadat.edlodlt = edlolocdlt;

#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
  MPI_Address (&srclstadat.velomin, &srcdisptab[0]);
  MPI_Address (&srclstadat.velodlt, &srcdisptab[1]);
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  MPI_Get_address (&srclstadat.velomin, &srcdisptab[0]);
  MPI_Get_address (&srclstadat.velodlt, &srcdisptab[1]);
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  srcdisptab[1] -= srcdisptab[0];
  srcdisptab[0] -= srcdisptab[0];

  o = 1;                                          /* Assume something will go wrong */
#if ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100))
  if ((MPI_Type_struct (2, dgraphstatblentab, srcdisptab, dgraphstattypetab, &srctypedat) == MPI_SUCCESS) &&
#else /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
  if ((MPI_Type_create_struct (2, dgraphstatblentab, srcdisptab, dgraphstattypetab, &srctypedat) == MPI_SUCCESS) &&
#endif /* ((defined COMMON_MPI_VERSION) && (COMMON_MPI_VERSION <= 100)) */
      (MPI_Type_commit (&srctypedat) == MPI_SUCCESS)) {
    if (MPI_Op_create ((MPI_User_function *) dgraphStatReduceAll, 0, &srcoperdat) == MPI_SUCCESS) {
      if (MPI_Allreduce (&srclstadat, &srcgstadat, 1, srctypedat, srcoperdat, srcgrafptr->proccomm) == MPI_SUCCESS)
        o = 0;

      MPI_Op_free (&srcoperdat);
    }
    MPI_Type_free (&srctypedat);
  }
  if (o != 0) {
    errorPrint ("SCOTCH_dgraphStat: communication error (2)");
    return     (1);
  }

  if (velominptr != NULL)
    *velominptr = (SCOTCH_Num) srcgstadat.velomin;
  if (velomaxptr != NULL)
    *velomaxptr = (SCOTCH_Num) srcgstadat.velomax;
  if (velosumptr != NULL)
    *velosumptr = (SCOTCH_Num) srcgrafptr->veloglbsum;
  if (veloavgptr != NULL)
    *veloavgptr = (double) veloglbavg;
  if (velodltptr != NULL)
    *velodltptr = srcgstadat.velodlt / (double) srcgrafptr->vertglbnbr;

  if (degrminptr != NULL)
    *degrminptr = (SCOTCH_Num) srcgstadat.degrmin;
  if (degrmaxptr != NULL)
    *degrmaxptr = (SCOTCH_Num) srcgstadat.degrmax;
  if (degravgptr != NULL)
    *degravgptr = (double) degrglbavg;
  if (degrdltptr != NULL)
    *degrdltptr = srcgstadat.degrdlt / (double) srcgrafptr->vertglbnbr;

  if (edlominptr != NULL)
    *edlominptr = (SCOTCH_Num) srcgstadat.edlomin;
  if (edlomaxptr != NULL)
    *edlomaxptr = (SCOTCH_Num) srcgstadat.edlomax;
  if (edlosumptr != NULL)
    *edlosumptr = (SCOTCH_Num) edloglbsum;
  if (edloavgptr != NULL)
    *edloavgptr = (double) edloglbavg;
  if (edlodltptr != NULL)
    *edlodltptr = srcgstadat.edlodlt / (double) srcgrafptr->edgeglbnbr;

  return (0);
}
