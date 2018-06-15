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
/**   NAME       : library_geom.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the geom     **/
/**                graph handling routines of the          **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 10 oct 1999     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     19 jan 2004     **/
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
#include "geom.h"
#include "graph.h"
#include "scotch.h"

/****************************************/
/*                                      */
/* These routines are the C API for the */
/* graph geometry handling routines.    */
/*                                      */
/****************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** geometry structure.
*** It returns:
*** - !NULL  : if the initialization succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Geom *
SCOTCH_geomAlloc ()
{
  return ((SCOTCH_Geom *) memAlloc (sizeof (SCOTCH_Geom)));
}

/*+ This routine initializes the opaque
*** geom structure used to handle graph
*** geometry in the Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_geomInit (
SCOTCH_Geom * const         geomptr)
{
  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint ("SCOTCH_geomInit: internal error (1)");
    return     (1);
  }
  if (sizeof (SCOTCH_Geom) < sizeof (Geom)) {
    errorPrint ("SCOTCH_geomInit: internal error (2)");
    return     (1);
  }

  return (geomInit ((Geom *) geomptr));
}

/*+ This routine frees the contents of the
*** given opaque geometry structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_geomExit (
SCOTCH_Geom * const         geomptr)
{
  geomExit ((Geom *) geomptr);
}

/*+ This routine accesses all of the geometry data.
*** NULL pointers on input indicate unwanted
*** data. NULL pointers on output indicate
*** unexisting arrays.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_geomData (
const SCOTCH_Geom * const   geomptr,              /* Geometry structure to read */
SCOTCH_Num * const          dimnptr,              /* Number of dimensions       */
double ** const             geomtab)              /* Geometry array [vertnbr]   */
{
  const Geom *        srcgeomptr;                 /* Pointer to source geometry structure */

  srcgeomptr = (const Geom *) geomptr;

  if (dimnptr != NULL)
    *dimnptr = srcgeomptr->dimnnbr;
  if (geomtab != NULL)
    *geomtab = (double *) srcgeomptr->geomtab;
}
