/* Copyright 2004,2007,2008,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : arch_mesh.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the mesh graph      **/
/**                target architectures.                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to   : 24 mar 1993     **/
/**                # Version 1.2  : from : 04 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 1.3  : from : 20 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 23 dec 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to   : 29 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     08 sep 1995     **/
/**                # Version 3.1  : from : 22 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 16 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     14 fev 2011     **/
/**                                                        **/
/**   NOTES      : # The vertices of the (dX,dY) mesh are  **/
/**                  numbered as terminals so that         **/
/**                  t(0,0) = 0, t(1,0) = 1,               **/
/**                  t(dX - 1, 0) = dX - 1, t(0,1) = dX,   **/
/**                  and t(x,y) = (y * dX) + x.            **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_MESH

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_mesh.h"

/***********************************************/
/*                                             */
/* These are the 2-dimensional mesh routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** bidimensional mesh architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archMesh2ArchLoad (
ArchMesh2 * restrict const  archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchMesh2)    > sizeof (ArchDummy)) ||
      (sizeof (ArchMesh2Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archMesh2ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1)) {
    errorPrint ("archMesh2ArchLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves the
** bidimensional mesh architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archMesh2ArchSave (
const ArchMesh2 * const     archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchMesh2)    > sizeof (ArchDummy)) ||
      (sizeof (ArchMesh2Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archMesh2ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0],
               (Anum) archptr->c[1]) == EOF) {
    errorPrint ("archMesh2ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archMesh2DomNum (
const ArchMesh2 * const     archptr,
const ArchMesh2Dom * const  domptr)
{
  return ((domptr->c[1][0] * archptr->c[0]) + domptr->c[0][0]); /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archMesh2DomTerm (
const ArchMesh2 * const     archptr,
ArchMesh2Dom * const        domptr,
const ArchDomNum            domnum)
{
  if (domnum < (archptr->c[0] * archptr->c[1])) { /* If valid label */
    domptr->c[0][0] =                             /* Set the domain */
    domptr->c[0][1] = domnum % archptr->c[0];
    domptr->c[1][0] =
    domptr->c[1][1] = domnum / archptr->c[0];

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the rectangular domain.
*/

Anum 
archMesh2DomSize (
const ArchMesh2 * const     archptr,
const ArchMesh2Dom * const  domptr)
{
  return ((domptr->c[0][1] - domptr->c[0][0] + 1) *
          (domptr->c[1][1] - domptr->c[1][0] + 1));
}

/* This function returns the average
** distance between two rectangular
** domains (in fact the distance between
** the centers of the domains).
*/

Anum 
archMesh2DomDist (
const ArchMesh2 * const    archptr,
const ArchMesh2Dom * const dom0ptr,
const ArchMesh2Dom * const dom1ptr)
{
  return (((abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
                 dom1ptr->c[0][0] - dom1ptr->c[0][1]) + 1) / 2) +
          ((abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
                 dom1ptr->c[1][0] - dom1ptr->c[1][1]) + 1) / 2));
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh2DomFrst (
const ArchMesh2 * const       archptr,
ArchMesh2Dom * restrict const domptr)
{
  domptr->c[0][0] =
  domptr->c[1][0] = 0;
  domptr->c[0][1] = archptr->c[0] - 1;
  domptr->c[1][1] = archptr->c[1] - 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh2DomLoad (
const ArchMesh2 * const       archptr,
ArchMesh2Dom * restrict const domptr,
FILE * restrict const         stream)
{
  if ((intLoad (stream, &domptr->c[0][0]) != 1) ||
      (intLoad (stream, &domptr->c[1][0]) != 1) ||
      (intLoad (stream, &domptr->c[0][1]) != 1) ||
      (intLoad (stream, &domptr->c[1][1]) != 1)) {
    errorPrint ("archMesh2DomLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves domain information
** to the given stream.
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh2DomSave (
const ArchMesh2 * const     archptr,
const ArchMesh2Dom * const  domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->c[0][0], (Anum) domptr->c[1][0],
               (Anum) domptr->c[0][1], (Anum) domptr->c[1][1]) == EOF) {
    errorPrint ("archMesh2DomSave: bad output");
    return     (1);
  }

  return (0);
}

/* These functions try to split a rectangular
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archMesh2DomBipart (
const ArchMesh2 * const       archptr,
const ArchMesh2Dom * const    domptr,
ArchMesh2Dom * restrict const dom0ptr,
ArchMesh2Dom * restrict const dom1ptr)
{
  Anum                dimsiz[2];
  int                 dimval;                     /* Dimension along which to split */

  dimsiz[0] = domptr->c[0][1] - domptr->c[0][0];
  dimsiz[1] = domptr->c[1][1] - domptr->c[1][0];

  if ((dimsiz[0] | dimsiz[1]) == 0)               /* Return if cannot bipartition more */
    return (1);

  dimval = 1;
  if ((dimsiz[0] > dimsiz[1]) ||                  /* Split domain in two along largest dimension */
      ((dimsiz[0] == dimsiz[1]) && (archptr->c[0] > archptr->c[1])))
    dimval = 0;

  if (dimval == 0) {                              /* Split across the X dimension */
    dom0ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];
  }
  else {                                          /* Split across the Y dimension */
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];
  }

  return (0);
}

int
archMesh2DomBipartO (
const ArchMesh2 * const       archptr,
const ArchMesh2Dom * const    domptr,
ArchMesh2Dom * restrict const dom0ptr,
ArchMesh2Dom * restrict const dom1ptr)
{
  if ((domptr->c[0][0] == domptr->c[0][1]) &&     /* Return if cannot bipartition more */
      (domptr->c[1][0] == domptr->c[1][1]))
    return (1);

  if (domptr->c[1][1] == domptr->c[1][0]) {       /* If the Y dimension cannot be cut */
    dom0ptr->c[0][0] = domptr->c[0][0];           /* Cut in the X dimension           */
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];
  }
  else {                                          /* If the Y dimension can be cut, cut it */
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];
  }

  return (0);
}

int
archMesh2DomBipartU (
const ArchMesh2 * const       archptr,
const ArchMesh2Dom * const    domptr,
ArchMesh2Dom * restrict const dom0ptr,
ArchMesh2Dom * restrict const dom1ptr)
{
  if ((domptr->c[0][0] == domptr->c[0][1]) &&     /* Return if cannot bipartition more */
      (domptr->c[1][0] == domptr->c[1][1]))
    return (1);

  if ((domptr->c[0][1] - domptr->c[0][0]) >       /* Split domain unevenly along largest dimension */
      (domptr->c[1][1] - domptr->c[1][0])) {
    dom0ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1] + domptr->c[0][1]) / 3;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];
  }
  else {
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1] + domptr->c[1][1]) / 3;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];
  }

  return (0);
}

/* This function checks if dom1 is
** included in dom0.
** It returns:
** - 0  : if dom1 is not included in dom0.
** - 1  : if dom1 is included in dom0.
** - 2  : on error.
*/

int
archMesh2DomIncl (
const ArchMesh2 * const     archptr,
const ArchMesh2Dom * const  dom0ptr,
const ArchMesh2Dom * const  dom1ptr)
{
  if ((dom0ptr->c[0][0] <= dom1ptr->c[0][0]) &&
      (dom0ptr->c[0][1] >= dom1ptr->c[0][1]) &&
      (dom0ptr->c[1][0] <= dom1ptr->c[1][0]) &&
      (dom0ptr->c[1][1] >= dom1ptr->c[1][1]))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** 2D mesh domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archMesh2DomMpiType (
const ArchMesh2 * const       archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (4, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */

/***********************************************/
/*                                             */
/* These are the 3-dimensional mesh routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** tridimensional mesh architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archMesh3ArchLoad (
ArchMesh3 * restrict const  archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchMesh3)    > sizeof (ArchDummy)) ||
      (sizeof (ArchMesh3Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archMesh3ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (intLoad (stream, &archptr->c[2]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1) || (archptr->c[2] < 1)) {
    errorPrint ("archMesh3ArchLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves the
** tridimensional mesh architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archMesh3ArchSave (
const ArchMesh3 * const     archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchMesh3)    > sizeof (ArchDummy)) ||
      (sizeof (ArchMesh3Dom) > sizeof (ArchDomDummy))) {
    errorPrint ("archMesh3ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0], (Anum) archptr->c[1], (Anum) archptr->c[2]) == EOF) {
    errorPrint ("archMesh3ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archMesh3DomNum (
const ArchMesh3 * const     archptr,
const ArchMesh3Dom * const  domptr)
{
  return ((((domptr->c[2][0]  * archptr->c[1]) +  /* Return the vertex number */
             domptr->c[1][0]) * archptr->c[0]) +
             domptr->c[0][0]);
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archMesh3DomTerm (
const ArchMesh3 * const     archptr,
ArchMesh3Dom * const        domptr,
const ArchDomNum            domnum)
{
  if (domnum < (archptr->c[0] * archptr->c[1] * archptr->c[2])) { /* If valid label */
    domptr->c[0][0] =                             /* Set the domain                 */
    domptr->c[0][1] = domnum % archptr->c[0];
    domptr->c[1][0] =
    domptr->c[1][1] = (domnum / archptr->c[0]) % archptr->c[1];
    domptr->c[2][0] =
    domptr->c[2][1] = domnum / (archptr->c[0] * archptr->c[1]);

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the cubic domain.
*/

Anum 
archMesh3DomSize (
const ArchMesh3 * const     archptr,
const ArchMesh3Dom * const  domptr)
{
  return ((domptr->c[0][1] - domptr->c[0][0] + 1) *
          (domptr->c[1][1] - domptr->c[1][0] + 1) *
          (domptr->c[2][1] - domptr->c[2][0] + 1));
}

/* This function returns the average distance
** between two cubic domains (in fact the
** distance between the centers of the domains).
*/

Anum 
archMesh3DomDist (
const ArchMesh3 * const     archptr,
const ArchMesh3Dom * const  dom0ptr,
const ArchMesh3Dom * const  dom1ptr)
{
  return (((abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
                 dom1ptr->c[0][0] - dom1ptr->c[0][1]) + 1) / 2) +
          ((abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
                 dom1ptr->c[1][0] - dom1ptr->c[1][1]) + 1) / 2) +
          ((abs (dom0ptr->c[2][0] + dom0ptr->c[2][1] -
                 dom1ptr->c[2][0] - dom1ptr->c[2][1]) + 1) / 2));
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh3DomFrst (
const ArchMesh3 * const       archptr,
ArchMesh3Dom * restrict const domptr)
{
  domptr->c[0][0] =
  domptr->c[1][0] =
  domptr->c[2][0] = 0;
  domptr->c[0][1] = archptr->c[0] - 1;
  domptr->c[1][1] = archptr->c[1] - 1;
  domptr->c[2][1] = archptr->c[2] - 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh3DomLoad (
const ArchMesh3 * const       archptr,
ArchMesh3Dom * restrict const domptr,
FILE * restrict const         stream)
{
  if ((intLoad (stream, &domptr->c[0][0]) != 1) ||
      (intLoad (stream, &domptr->c[1][0]) != 1) ||
      (intLoad (stream, &domptr->c[2][0]) != 1) ||
      (intLoad (stream, &domptr->c[0][1]) != 1) ||
      (intLoad (stream, &domptr->c[1][1]) != 1) ||
      (intLoad (stream, &domptr->c[2][1]) != 1)) {
    errorPrint ("archMesh3DomLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archMesh3DomSave (
const ArchMesh3 * const     archptr,
const ArchMesh3Dom * const  domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->c[0][0], (Anum) domptr->c[1][0], (Anum) domptr->c[2][0],
               (Anum) domptr->c[0][1], (Anum) domptr->c[1][1], (Anum) domptr->c[2][1]) == EOF) {
    errorPrint ("archMesh3DomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a cubic
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archMesh3DomBipart (
const ArchMesh3 * const       archptr,
const ArchMesh3Dom * const    domptr,
ArchMesh3Dom * restrict const dom0ptr,
ArchMesh3Dom * restrict const dom1ptr)
{
  Anum                dimsiz[3];
  int                 dimtmp;
  int                 dimval;

  dimsiz[0] = domptr->c[0][1] - domptr->c[0][0];
  dimsiz[1] = domptr->c[1][1] - domptr->c[1][0];
  dimsiz[2] = domptr->c[2][1] - domptr->c[2][0];

  if ((dimsiz[0] | dimsiz[1] | dimsiz[2]) == 0)   /* Return if cannot bipartition more */
    return (1);

  dimval = (archptr->c[1] > archptr->c[0]) ? 1 : 0; /* Assume all subdomain dimensions are equal */
  if (archptr->c[2] > archptr->c[dimval])         /* Find priviledged dimension                  */
    dimval = 2;

  dimtmp = dimval;                                /* Find best dimension */
  if (dimsiz[(dimtmp + 1) % 3] > dimsiz[dimval])
    dimval = (dimtmp + 1) % 3;
  if (dimsiz[(dimtmp + 2) % 3] > dimsiz[dimval])
    dimval = (dimtmp + 2) % 3;

  if (dimval == 0) {                              /* Split domain in two along largest dimension */
    dom0ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = dom1ptr->c[2][1] = domptr->c[2][1];
  }
  else if (dimval == 1) {
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = dom1ptr->c[2][1] = domptr->c[2][1];
  }
  else {
    dom0ptr->c[0][0] = dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = (domptr->c[2][0] + domptr->c[2][1]) / 2;
    dom1ptr->c[2][0] = dom0ptr->c[2][1] + 1;
    dom1ptr->c[2][1] = domptr->c[2][1];
  }

  return (0);
}

/* This function checks if dom1 is
** included in dom0.
** It returns:
** - 0  : if dom1 is not included in dom0.
** - 1  : if dom1 is included in dom0.
** - 2  : on error.
*/

int
archMesh3DomIncl (
const ArchMesh3 * const     archptr,
const ArchMesh3Dom * const  dom0ptr,
const ArchMesh3Dom * const  dom1ptr)
{
  if ((dom0ptr->c[0][0] <= dom1ptr->c[0][0]) &&
      (dom0ptr->c[0][1] >= dom1ptr->c[0][1]) &&
      (dom0ptr->c[1][0] <= dom1ptr->c[1][0]) &&
      (dom0ptr->c[1][1] >= dom1ptr->c[1][1]) &&
      (dom0ptr->c[2][0] <= dom1ptr->c[2][0]) &&
      (dom0ptr->c[2][1] >= dom1ptr->c[2][1]))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** 3D mesh domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archMesh3DomMpiType (
const ArchMesh3 * const       archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (6, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
