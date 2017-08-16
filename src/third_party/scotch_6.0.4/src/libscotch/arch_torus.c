/* Copyright 2004,2007,2008,2010,2011,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_torus.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the torus graph     **/
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
/**                # Version 3.1  : from : 07 may 1996     **/
/**                                 to     22 jul 1996     **/
/**                # Version 3.2  : from : 16 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     30 nov 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_TORUS

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_torus.h"

/***********************************************/
/*                                             */
/* These are the 2-dimensional torus routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** bidimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTorus2ArchLoad (
ArchTorusX * restrict const archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((ARCHTORUSDIMMAX        < 2)                  ||
      (sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus2ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1)) {
    errorPrint ("archTorus2ArchLoad: bad input");
    return     (1);
  }

  archptr->dimmax = 2;

  return (0);
}

/* This routine saves the
** bidimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTorus2ArchSave (
const ArchTorusX * const    archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((ARCHTORUSDIMMAX        < 2)                  ||
      (sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus2ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0],
               (Anum) archptr->c[1]) == EOF) {
    errorPrint ("archTorus2ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archTorus2DomNum (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
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
archTorus2DomTerm (
const ArchTorusX * const    archptr,
ArchTorusXDom * const       domptr,
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
archTorus2DomSize (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
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
archTorus2DomDist (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
{
  Anum               dc0, dc1;
  Anum               ds0, ds1;

  dc0 = abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
             dom1ptr->c[0][0] - dom1ptr->c[0][1]);
  ds0 = (dc0 > archptr->c[0]) ? (2 * archptr->c[0] - dc0) : dc0;

  dc1 = abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
             dom1ptr->c[1][0] - dom1ptr->c[1][1]);
  ds1 = (dc1 > archptr->c[1]) ? (2 * archptr->c[1] - dc1) : dc1;

  return ((ds0 + ds1) >> 1);
}

/* This function tries to split a rectangular
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archTorus2DomBipart (
const ArchTorusX * const        archptr,
const ArchTorusXDom * const     domptr,
ArchTorusXDom * restrict const  dom0ptr,
ArchTorusXDom * restrict const  dom1ptr)
{
  Anum                dimsiz0;
  Anum                dimsiz1;
  int                 dimnum;                     /* Dimension along which to split */

  dimsiz0 = domptr->c[0][1] - domptr->c[0][0];
  dimsiz1 = domptr->c[1][1] - domptr->c[1][0];

  if ((dimsiz0 | dimsiz1) == 0)                   /* Return if cannot bipartition more */
    return (1);

  dimnum = 1;
  if ((dimsiz0 > dimsiz1) ||                      /* Split domain in two along largest dimension */
      ((dimsiz0 == dimsiz1) && (archptr->c[0] > archptr->c[1])))
    dimnum = 0;

  dom0ptr->c[0][0] = domptr->c[0][0];
  dom1ptr->c[1][1] = domptr->c[1][1];
  if (dimnum == 0) {                              /* Split across the X dimension */
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = domptr->c[1][1];
  }
  else {                                          /* Split across the Y dimension */
    dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];
    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
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
archTorus2DomIncl (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
{
  if ((dom0ptr->c[0][0] <= dom1ptr->c[0][0]) &&
      (dom0ptr->c[0][1] >= dom1ptr->c[0][1]) &&
      (dom0ptr->c[1][0] <= dom1ptr->c[1][0]) &&
      (dom0ptr->c[1][1] >= dom1ptr->c[1][1]))
    return (1);

  return (0);
}

/***********************************************/
/*                                             */
/* These are the 3-dimensional torus routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTorus3ArchLoad (
ArchTorusX * restrict const archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((ARCHTORUSDIMMAX        < 3)                  ||
      (sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus3ArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->c[0]) != 1) ||
      (intLoad (stream, &archptr->c[1]) != 1) ||
      (intLoad (stream, &archptr->c[2]) != 1) ||
      (archptr->c[0] < 1) || (archptr->c[1] < 1) || (archptr->c[2] < 1)) {
    errorPrint ("archTorus3ArchLoad: bad input");
    return     (1);
  }

  archptr->dimmax = 3;

  return (0);
}

/* This routine saves the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTorus3ArchSave (
const ArchTorusX * const    archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorus3ArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) archptr->c[0], (Anum) archptr->c[1], (Anum) archptr->c[2]) == EOF) {
    errorPrint ("archTorus3ArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archTorus3DomNum (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
{
  return ((((domptr->c[2][0]  * archptr->c[1]) +  /* Return vertex number */
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
archTorus3DomTerm (
const ArchTorusX * const    archptr,
ArchTorusXDom * const       domptr,
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
archTorus3DomSize (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
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
archTorus3DomDist (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
{
  Anum               dc0, dc1, dc2;
  Anum               ds0, ds1, ds2;

  dc0 = abs (dom0ptr->c[0][0] + dom0ptr->c[0][1] -
             dom1ptr->c[0][0] - dom1ptr->c[0][1]);
  ds0 = (dc0 > archptr->c[0]) ? (2 * archptr->c[0] - dc0) : dc0;

  dc1 = abs (dom0ptr->c[1][0] + dom0ptr->c[1][1] -
             dom1ptr->c[1][0] - dom1ptr->c[1][1]);
  ds1 = (dc1 > archptr->c[1]) ? (2 * archptr->c[1] - dc1) : dc1;

  dc2 = abs (dom0ptr->c[2][0] + dom0ptr->c[2][1] -
             dom1ptr->c[2][0] - dom1ptr->c[2][1]);
  ds2 = (dc2 > archptr->c[2]) ? (2 * archptr->c[2] - dc2) : dc2;

  return ((ds0 + ds1 + ds2) >> 1);
}

/* This function tries to split a cubic
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archTorus3DomBipart (
const ArchTorusX * const        archptr,
const ArchTorusXDom * const     domptr,
ArchTorusXDom * restrict const  dom0ptr,
ArchTorusXDom * restrict const  dom1ptr)
{
  Anum                dimsiz[3];
  int                 dimnum;

  dimsiz[0] = domptr->c[0][1] - domptr->c[0][0];
  dimsiz[1] = domptr->c[1][1] - domptr->c[1][0];
  dimsiz[2] = domptr->c[2][1] - domptr->c[2][0];

  if ((dimsiz[0] | dimsiz[1] | dimsiz[2]) == 0)   /* Return if cannot bipartition more */
    return (1);

  dimnum = ((dimsiz[1] > dimsiz[2]) ||            /* Find largest or priviledged subdomain dimension */
            ((dimsiz[1] == dimsiz[2]) &&
             (archptr->c[1] > archptr->c[2]))) ? 1 : 2;
  if ((dimsiz[0] > dimsiz[dimnum]) ||
      ((dimsiz[0] == dimsiz[dimnum]) &&
       (archptr->c[0] > archptr->c[dimnum])))
    dimnum = 0;

  dom0ptr->c[0][0] = domptr->c[0][0];
  dom1ptr->c[2][1] = domptr->c[2][1];
  if (dimnum == 0) {                              /* Split domain in two along largest dimension */
    dom0ptr->c[0][1] = (domptr->c[0][0] + domptr->c[0][1]) / 2;
    dom1ptr->c[0][0] = dom0ptr->c[0][1] + 1;
    dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = domptr->c[2][1];
  }
  else if (dimnum == 1) {
    dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = (domptr->c[1][0] + domptr->c[1][1]) / 2;
    dom1ptr->c[1][0] = dom0ptr->c[1][1] + 1;
    dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = dom1ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = domptr->c[2][1];
  }
  else {
    dom1ptr->c[0][0] = domptr->c[0][0];
    dom0ptr->c[0][1] = dom1ptr->c[0][1] = domptr->c[0][1];

    dom0ptr->c[1][0] = dom1ptr->c[1][0] = domptr->c[1][0];
    dom0ptr->c[1][1] = dom1ptr->c[1][1] = domptr->c[1][1];

    dom0ptr->c[2][0] = domptr->c[2][0];
    dom0ptr->c[2][1] = (domptr->c[2][0] + domptr->c[2][1]) / 2;
    dom1ptr->c[2][0] = dom0ptr->c[2][1] + 1;
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
archTorus3DomIncl (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
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

/***********************************************/
/*                                             */
/* These are the x-dimensional torus routines. */
/*                                             */
/***********************************************/

/* This routine loads the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTorusXArchLoad (
ArchTorusX * restrict const archptr,
FILE * restrict const       stream)
{
  Anum                dimnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorusXArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->dimmax) != 1) ||
      (archptr->dimmax > ARCHTORUSDIMMAX)) {
    errorPrint ("archTorusXArchLoad: bad input (1)");
    return     (1);
  }

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    if ((intLoad (stream, &archptr->c[dimnum]) != 1) ||
        (archptr->c[dimnum] < 1)) {
      errorPrint ("archTorusXArchLoad: bad input (2)");
      return     (1);
    }
  }

  return (0);
}

/* This routine saves the
** tridimensional torus architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTorusXArchSave (
const ArchTorusX * const    archptr,
FILE * restrict const       stream)
{
  Anum                dimnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTorusX)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTorusXDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTorusXArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " ",
               (Anum) archptr->dimmax) == EOF) {
    errorPrint ("archTorusXArchSave: bad output (1)");
    return     (1);
  }

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    if (fprintf (stream, ANUMSTRING " ",
               (Anum) archptr->c[dimnum]) == EOF) {
      errorPrint ("archTorusXArchSave: bad output (2)");
      return     (1);
    }
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archTorusXDomNum (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
{
  Anum                dimnum;
  Anum                domnum;

  for (dimnum = archptr->dimmax - 2, domnum = domptr->c[archptr->dimmax - 1][0]; dimnum >= 0; dimnum --)
    domnum = (domnum * archptr->c[dimnum]) + domptr->c[dimnum][0];

  return (domnum);                                /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archTorusXDomTerm (
const ArchTorusX * const    archptr,
ArchTorusXDom * const       domptr,
const ArchDomNum            domnum)
{
  Anum                dimnum;
  Anum                domtmp;

  for (dimnum = 0, domtmp = domnum; dimnum < archptr->dimmax; dimnum ++) { /* Set the domain */
    domptr->c[dimnum][0] =
    domptr->c[dimnum][1] = domtmp % archptr->c[dimnum];
    domtmp /= archptr->c[dimnum];
  }

  if (domtmp > 0)                                 /* If residual is not zero, terminal domain number is invalid since too high */
    return (1);

  return (0);
}

/* This function returns the number of
** elements in the cubic domain.
*/

Anum 
archTorusXDomSize (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr)
{
  Anum                dimnum;
  Anum                domsiz;

  for (dimnum = 0, domsiz = 1; dimnum < archptr->dimmax; dimnum ++)
    domsiz *= domptr->c[dimnum][1] - domptr->c[dimnum][0] + 1;

  return (domsiz);
}

/* This function returns the average distance
** between two cubic domains (in fact the
** distance between the centers of the domains).
*/

Anum 
archTorusXDomDist (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
{
  Anum                dimnum;
  Anum                distval;
  Anum                disttmp;

  for (dimnum = 0, distval = 0; dimnum < archptr->dimmax; dimnum ++) {
    disttmp = abs (dom0ptr->c[dimnum][0] + dom0ptr->c[dimnum][1] -
                   dom1ptr->c[dimnum][0] - dom1ptr->c[dimnum][1]);
    distval += (disttmp > archptr->c[dimnum]) ? (2 * archptr->c[dimnum] - disttmp) : disttmp;
  }

  return (distval >> 1);
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorusXDomFrst (
const ArchTorusX * const        archptr,
ArchTorusXDom * restrict const  domptr)
{
  Anum                dimnum;

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    domptr->c[dimnum][0] = 0;
    domptr->c[dimnum][1] = archptr->c[dimnum] - 1;
  }

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTorusXDomLoad (
const ArchTorusX * const        archptr,
ArchTorusXDom * restrict const  domptr,
FILE * restrict const           stream)
{
  Anum                dimnum;

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    if ((intLoad (stream, &domptr->c[dimnum][0]) != 1) ||
        (intLoad (stream, &domptr->c[dimnum][1]) != 1) ||
        (domptr->c[dimnum][0] > domptr->c[dimnum][1])  ||
        (domptr->c[dimnum][0] < 0)) {
      errorPrint ("archTorusXDomLoad: bad input");
      return     (1);
    }
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
archTorusXDomSave (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const domptr,
FILE * restrict const       stream)
{
  Anum                dimnum;

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
                 (Anum) domptr->c[dimnum][0],
                 (Anum) domptr->c[dimnum][1]) == EOF) {
      errorPrint ("archTorusXDomSave: bad output");
      return     (1);
    }
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
archTorusXDomBipart (
const ArchTorusX * const        archptr,
const ArchTorusXDom * const     domptr,
ArchTorusXDom * restrict const  dom0ptr,
ArchTorusXDom * restrict const  dom1ptr)
{
  Anum                archdimsizmax;              /* Maximum span on largest architecture dimension */
  Anum                domndimsizmax;              /* Maximum span on largest domain dimension       */
  Anum                domndimval;                 /* Dimension to be split                          */
  Anum                domndimflg;                 /* Flag set if subdomain can be bipartitioned     */
  Anum                domndimtmp;
  Anum                dimnum;

  for (dimnum = domndimval = archptr->dimmax - 1, archdimsizmax = domndimflg = 0, domndimsizmax = -1;
       dimnum >= 0; dimnum --) {
    Anum                archdimsiz;
    Anum                domndimsiz;
    Anum                domndim0;
    Anum                domndim1;

    dom0ptr->c[dimnum][0] =                       /* Set up subdomain data as copy of original domain data */
    dom1ptr->c[dimnum][0] = domndim0 = domptr->c[dimnum][0];
    dom0ptr->c[dimnum][1] =
    dom1ptr->c[dimnum][1] = domndim1 = domptr->c[dimnum][1];

    domndimsiz  = domndim1 - domndim0;            /* Span on current dimension            */
    domndimflg |= domndimsiz;                     /* Flag set if at least one is not zero */

    if (domndimsiz < domndimsizmax)               /* If dimension is too small, skip it */
      continue;
    archdimsiz = archptr->c[dimnum];
    if ((domndimsiz == domndimsizmax) &&          /* If dimension to split is not priviledged, skip it */
        (archdimsiz <= archdimsizmax))
      continue;

    archdimsizmax = archdimsiz;                   /* Record dimension to split */
    domndimsizmax = domndimsiz;
    domndimval    = dimnum;
  }

  if (domndimflg == 0)                            /* Return if cannot bipartition more */
    return (1);

  domndimtmp = (domptr->c[domndimval][0] + domptr->c[domndimval][1]) / 2;
  dom0ptr->c[domndimval][1] = domndimtmp;
  dom1ptr->c[domndimval][0] = domndimtmp + 1;

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
archTorusXDomIncl (
const ArchTorusX * const    archptr,
const ArchTorusXDom * const dom0ptr,
const ArchTorusXDom * const dom1ptr)
{
  Anum                dimnum;

  for (dimnum = 0; dimnum < archptr->dimmax; dimnum ++) {
    if ((dom1ptr->c[dimnum][0] < dom0ptr->c[dimnum][0]) ||
        (dom1ptr->c[dimnum][1] > dom1ptr->c[dimnum][1]))
      return (0);
  }

  return (1);
}

/* This function creates the MPI_Datatype for
** xD torus domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archTorusXDomMpiType (
const ArchTorusX * const      archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (2 * archptr->dimmax, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
