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
/**   NAME       : library_mesh.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the source   **/
/**                mesh handling routines of the           **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 23 sep 2002     **/
/**                                 to     11 may 2004     **/
/**                # Version 5.1  : from : 17 nov 2010     **/
/**                                 to     17 nov 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mesh handling routines.      */
/*                                  */
/************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** centralized mesh structure.
*** It returns:
*** - !NULL  : if the initialization succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Mesh *
SCOTCH_meshAlloc ()
{
  return ((SCOTCH_Mesh *) memAlloc (sizeof (SCOTCH_Mesh)));
}

/*+ This routine initializes the opaque
*** mesh structure used to handle meshes
*** in the Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_meshInit (
SCOTCH_Mesh * const         meshptr)
{
  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint ("SCOTCH_meshInit: internal error (1)");
    return     (1);
  }
  if (sizeof (SCOTCH_Mesh) < sizeof (Mesh)) {
    errorPrint ("SCOTCH_meshInit: internal error (2)");
    return     (1);
  }

  return (meshInit ((Mesh *) meshptr));
}

/*+ This routine frees the contents of the
*** given opaque mesh structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_meshExit (
SCOTCH_Mesh * const         meshptr)
{
  meshExit ((Mesh *) meshptr);
}

/*+ This routine loads the given opaque mesh
*** structure with the data of the given stream.
*** The base value allows the user to set the
*** mesh base to 0 or 1, or to the base value
*** of the stream if the base value is equal
*** to -1.
*** It returns:
*** - 0   : if the loading succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_meshLoad (
SCOTCH_Mesh * const         meshptr,
FILE * const                stream,
const SCOTCH_Num            baseval)
{
  if ((baseval < -1) || (baseval > 1)) {
    errorPrint ("SCOTCH_meshLoad: invalid base parameter");
    return     (1);
  }

  return (meshLoad ((Mesh * const) meshptr, stream, (Gnum) baseval));
}

/*+ This routine saves the contents of the given
*** opaque mesh structure to the given stream.
*** It returns:
*** - 0   : if the saving succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_meshSave (
const SCOTCH_Mesh * const   meshptr,
FILE * const                stream)
{
  return (meshSave ((const Mesh * const) meshptr, stream));
}

/*+ This routine fills the contents of the given
*** opaque mesh structure with the data provided
*** by the user. The base value allows the user to
*** set the mesh base to 0 or 1.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshBuild (
SCOTCH_Mesh * const         meshptr,              /* Mesh structure to fill               */
const SCOTCH_Num            velmbas,              /* Base index for element vertices      */
const SCOTCH_Num            vnodbas,              /* Base index for node vertices         */
const SCOTCH_Num            velmnbr,              /* Number of elements in mesh graph     */
const SCOTCH_Num            vnodnbr,              /* Number of vertices in mesh graph     */
const SCOTCH_Num * const    verttab,              /* Vertex array [vertnbr or vertnbr+1]  */
const SCOTCH_Num * const    vendtab,              /* Vertex end array [vertnbr]           */
const SCOTCH_Num * const    velotab,              /* Element vertex load array            */
const SCOTCH_Num * const    vnlotab,              /* Node vertex load array               */
const SCOTCH_Num * const    vlbltab,              /* Vertex label array                   */
const SCOTCH_Num            edgenbr,              /* Number of edges (arcs)               */
const SCOTCH_Num * const    edgetab)              /* Edge array [edgenbr]                 */
{
  Mesh *              srcmeshptr;                 /* Pointer to source mesh structure    */
  Gnum                degrmax;                    /* Maximum degree                      */
  Gnum                veisnbr;                    /* Number of isolated element vertices */
  Gnum                vertnum;                    /* Current vertex number               */

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (sizeof (SCOTCH_Mesh) < sizeof (Mesh)) {
    errorPrint ("SCOTCH_meshBuild: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
  if ((velmbas < 0) ||
      (vnodbas < 0) ||
      ((velmbas > 1) && (vnodbas > 1))) {
    errorPrint ("SCOTCH_meshBuild: invalid base parameters");
    return     (1);
  }
  if (((velmbas + velmnbr) != vnodbas) &&
      ((vnodbas + vnodnbr) != velmbas)) {
    errorPrint ("SCOTCH_meshBuild: invalid element or node range");
    return     (1);
  }

  srcmeshptr = (Mesh *) meshptr;                  /* Use structure as source mesh */

  srcmeshptr->flagval = MESHNONE;
  srcmeshptr->baseval = MIN (velmbas, vnodbas);
  srcmeshptr->velmnbr = velmnbr;
  srcmeshptr->velmbas = velmbas;
  srcmeshptr->velmnnd = velmbas + velmnbr;
  srcmeshptr->vnodnbr = vnodnbr;
  srcmeshptr->vnodbas = vnodbas;
  srcmeshptr->vnodnnd = vnodbas + vnodnbr;
  srcmeshptr->verttax = (Gnum *) verttab - srcmeshptr->baseval;
  srcmeshptr->vendtax = ((vendtab == NULL) || (vendtab == verttab) || (vendtab == verttab + 1))
                          ? srcmeshptr->verttax + 1 : (Gnum *) vendtab - srcmeshptr->baseval;
  srcmeshptr->velotax = ((velotab == NULL) || (velotab == verttab)) ? NULL : (Gnum *) velotab - srcmeshptr->velmbas;
  srcmeshptr->vnlotax = ((vnlotab == NULL) || (vnlotab == verttab)) ? NULL : (Gnum *) vnlotab - srcmeshptr->vnodbas;
  srcmeshptr->vlbltax = ((vlbltab == NULL) || (vlbltab == verttab)) ? NULL : (Gnum *) vlbltab - srcmeshptr->baseval;
  srcmeshptr->edgenbr = edgenbr;
  srcmeshptr->edgetax = (Gnum *) edgetab - srcmeshptr->baseval;

  if (srcmeshptr->velotax == NULL)                /* Compute element vertex load sum */
    srcmeshptr->velosum = srcmeshptr->velmnbr;
  else {
    Gnum                velosum;                  /* Sum of element vertex loads */

    for (vertnum = srcmeshptr->velmbas, velosum = 0;
         vertnum < srcmeshptr->velmnnd; vertnum ++)
      velosum += srcmeshptr->velotax[vertnum];

    srcmeshptr->velosum = velosum;
  }
  if (srcmeshptr->vnlotax == NULL)                /* Compute node vertex load sum */
    srcmeshptr->vnlosum = srcmeshptr->vnodnbr;
  else {
    Gnum                vnlosum;                  /* Sum of node vertex loads */

    for (vertnum = srcmeshptr->vnodbas, vnlosum = 0;
         vertnum < srcmeshptr->vnodnnd; vertnum ++)
      vnlosum += srcmeshptr->vnlotax[vertnum];

    srcmeshptr->vnlosum = vnlosum;
  }

  for (vertnum = srcmeshptr->velmbas, veisnbr = degrmax = 0; /* Compute maximum degree */
       vertnum < srcmeshptr->velmnnd; vertnum ++) {
    Gnum                degrval;                  /* Degree of current vertex */

    degrval = srcmeshptr->vendtax[vertnum] - srcmeshptr->verttax[vertnum];
    if (degrval > degrmax)
      degrmax = degrval;
    else if (degrval == 0)                        /* Count number of isolated element vertices */
      veisnbr ++;
  }
  srcmeshptr->veisnbr = veisnbr;

  for (vertnum = srcmeshptr->vnodbas;             /* Compute maximum degree */
       vertnum < srcmeshptr->vnodnnd; vertnum ++) {
    Gnum                degrval;                  /* Degree of current vertex */

    degrval = srcmeshptr->vendtax[vertnum] - srcmeshptr->verttax[vertnum];
    if (degrval > degrmax)
      degrmax = degrval;
  }
  srcmeshptr->degrmax = degrmax;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (meshCheck (srcmeshptr) != 0) {
    errorPrint ("SCOTCH_meshBuild: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  return (0);
}

/*+ This routine checks the consistency
*** of the given mesh.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_meshCheck (
const SCOTCH_Mesh * const   meshptr)
{
  return (meshCheck ((const Mesh * const) meshptr));
}

/*+ This routine accesses mesh size data.
*** NULL pointers on input indicate unwanted
*** data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_meshSize (
const SCOTCH_Mesh * const   meshptr,
SCOTCH_Num * const          velmnbr,
SCOTCH_Num * const          vnodnbr,
SCOTCH_Num * const          edgenbr)
{
  const Mesh *        srcmeshptr;

  srcmeshptr = (Mesh *) meshptr;

  if (velmnbr != NULL)
    *velmnbr = (SCOTCH_Num) srcmeshptr->velmnbr;
  if (vnodnbr != NULL)
    *vnodnbr = (SCOTCH_Num) srcmeshptr->vnodnbr;
  if (edgenbr != NULL)
    *edgenbr = (SCOTCH_Num) srcmeshptr->edgenbr;
}

/*+ This routine accesses all of the mesh data.
*** NULL pointers on input indicate unwanted
*** data. NULL pointers on output indicate
*** unexisting arrays.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_meshData (
const SCOTCH_Mesh * const   meshptr,              /* Mesh structure to read    */
SCOTCH_Num * const          velmbas,              /* Base index for elements   */
SCOTCH_Num * const          vnodbas,              /* Base index for nodes      */
SCOTCH_Num * const          velmnbr,              /* Number of elements        */
SCOTCH_Num * const          vnodnbr,              /* Number of nodes           */
SCOTCH_Num ** const         verttab,              /* Vertex array [vertnbr+1]  */
SCOTCH_Num ** const         vendtab,              /* Vertex array [vertnbr]    */
SCOTCH_Num ** const         velotab,              /* Element vertex load array */
SCOTCH_Num ** const         vnlotab,              /* Vertex load array         */
SCOTCH_Num ** const         vlbltab,              /* Vertex label array        */
SCOTCH_Num * const          edgenbr,              /* Number of edges (arcs)    */
SCOTCH_Num ** const         edgetab,              /* Edge array [edgenbr]      */
SCOTCH_Num * const          degrnbr)              /* Maximum degree            */
{
  const Mesh *       srcmeshptr;                  /* Pointer to source mesh structure */

  srcmeshptr = (const Mesh *) meshptr;

  if (velmnbr != NULL)
    *velmnbr = srcmeshptr->velmnbr;
  if (vnodnbr != NULL)
    *vnodnbr = srcmeshptr->vnodnbr;
  if (velmbas != NULL)
    *velmbas = srcmeshptr->velmbas;
  if (vnodbas != NULL)
    *vnodbas = srcmeshptr->vnodbas;
  if (verttab != NULL)
    *verttab = srcmeshptr->verttax + srcmeshptr->baseval;
  if (vendtab != NULL)
    *vendtab = srcmeshptr->vendtax + srcmeshptr->baseval;
  if (velotab != NULL)
    *velotab = (srcmeshptr->velotax != NULL) ? (srcmeshptr->velotax + srcmeshptr->velmbas) : NULL;
  if (vnlotab != NULL)
    *vnlotab = (srcmeshptr->vnlotax != NULL) ? (srcmeshptr->vnlotax + srcmeshptr->vnodbas) : NULL;
  if (vlbltab != NULL)
    *vlbltab = (srcmeshptr->vlbltax != NULL) ? (srcmeshptr->vlbltax + srcmeshptr->baseval) : NULL;
  if (edgenbr != NULL)
    *edgenbr = srcmeshptr->edgenbr;
  if (edgetab != NULL)
    *edgetab = srcmeshptr->edgetax + srcmeshptr->baseval;
  if (degrnbr != NULL)
    *degrnbr = srcmeshptr->degrmax;
}

/*+ This routine computes statistics
*** on the given graph.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_meshStat (
const SCOTCH_Mesh * const   meshptr,
SCOTCH_Num * const          vnlominptr,           /* Vertex loads only for nodes */
SCOTCH_Num * const          vnlomaxptr,
SCOTCH_Num * const          vnlosumptr,
double * const              vnloavgptr,
double * const              vnlodltptr,
SCOTCH_Num * const          edegminptr,           /* Element degree data */
SCOTCH_Num * const          edegmaxptr,
double * const              edegavgptr,
double * const              edegdltptr,
SCOTCH_Num * const          ndegminptr,           /* Node degree data */
SCOTCH_Num * const          ndegmaxptr,
double * const              ndegavgptr,
double * const              ndegdltptr)
{
  const Mesh * restrict srcmeshptr;
  Gnum                  vnlomin;
  Gnum                  vnlomax;
  double                vnloavg;
  double                vnlodlt;
  Gnum                  degrmin;
  Gnum                  degrmax;
  double                degravg;
  double                degrdlt;

  srcmeshptr = (Mesh *) meshptr;

  vnlodlt = 0.0L;
  if (srcmeshptr->vnodnbr > 0) {
    if (srcmeshptr->vnlotax != NULL) {            /* If mesh has node vertex loads */
      Gnum                  vnodnum;

      vnlomin = GNUMMAX;
      vnlomax = 0;
      vnloavg = (double) srcmeshptr->vnlosum / (double) srcmeshptr->vnodnbr;

      for (vnodnum = srcmeshptr->vnodbas; vnodnum < srcmeshptr->vnodnnd; vnodnum ++) {
        if (srcmeshptr->vnlotax[vnodnum] < vnlomin) /* Account for vertex loads */
          vnlomin = srcmeshptr->vnlotax[vnodnum];
        if (srcmeshptr->vnlotax[vnodnum] > vnlomax)
          vnlomax = srcmeshptr->vnlotax[vnodnum];
        vnlodlt += fabs ((double) srcmeshptr->vnlotax[vnodnum] - vnloavg);
      }
      vnlodlt /= (double) srcmeshptr->vnodnbr;
    }
    else {
      vnlomin =
      vnlomax = 1;
      vnloavg = 1.0L;
    }
  }
  else {
    vnlomin =
    vnlomax = 0;
    vnloavg = 0.0L;
  }

  if (vnlominptr != NULL)
    *vnlominptr = (SCOTCH_Num) vnlomin;
  if (vnlomaxptr != NULL)
    *vnlomaxptr = (SCOTCH_Num) vnlomax;
  if (vnlosumptr != NULL)
    *vnlosumptr = (SCOTCH_Num) srcmeshptr->vnlosum;
  if (vnloavgptr != NULL)
    *vnloavgptr = (double) vnloavg;
  if (vnlodltptr != NULL)
    *vnlodltptr = (double) vnlodlt;

  degrmax = 0;
  degrdlt = 0.0L;
  if (srcmeshptr->velmnbr > 0) {
    Gnum                  velmnum;

    degrmin = GNUMMAX;
    degravg = (double) srcmeshptr->edgenbr / (double) (2 * srcmeshptr->velmnbr);
    for (velmnum = srcmeshptr->velmbas; velmnum < srcmeshptr->velmnnd; velmnum ++) {
      Gnum                  degrval;

      degrval = srcmeshptr->vendtax[velmnum] - srcmeshptr->verttax[velmnum]; /* Get element degree */
      if (degrval < degrmin)
        degrmin = degrval;
      if (degrval > degrmax)
        degrmax = degrval;
      degrdlt += fabs ((double) degrval - degravg);
    }
    degrdlt /= (double) srcmeshptr->velmnbr;
  }
  else {
    degrmin = 0;
    degravg = 0.0L;
  }

  if (edegminptr != NULL)
    *edegminptr = (SCOTCH_Num) degrmin;
  if (edegmaxptr != NULL)
    *edegmaxptr = (SCOTCH_Num) degrmax;
  if (edegavgptr != NULL)
    *edegavgptr = (double) degravg;
  if (edegdltptr != NULL)
    *edegdltptr = (double) degrdlt;

  degrmax = 0;
  degrdlt = 0.0L;
  if (srcmeshptr->vnodnbr > 0) {
    Gnum                  vnodnum;

    degrmin = GNUMMAX;
    degravg = (double) srcmeshptr->edgenbr / (double) (2 * srcmeshptr->vnodnbr);
    for (vnodnum = srcmeshptr->vnodbas; vnodnum < srcmeshptr->vnodnnd; vnodnum ++) {
      Gnum                  degrval;

      degrval = srcmeshptr->vendtax[vnodnum] - srcmeshptr->verttax[vnodnum]; /* Get element degree */
      if (degrval < degrmin)
        degrmin = degrval;
      if (degrval > degrmax)
        degrmax = degrval;
      degrdlt += fabs ((double) degrval - degravg);
    }
    degrdlt /= (double) srcmeshptr->vnodnbr;
  }
  else {
    degrmin = 0;
    degravg = 0.0L;
  }

  if (ndegminptr != NULL)
    *ndegminptr = (SCOTCH_Num) degrmin;
  if (ndegmaxptr != NULL)
    *ndegmaxptr = (SCOTCH_Num) degrmax;
  if (ndegavgptr != NULL)
    *ndegavgptr = (double) degravg;
  if (ndegdltptr != NULL)
    *ndegdltptr = (double) degrdlt;
}
