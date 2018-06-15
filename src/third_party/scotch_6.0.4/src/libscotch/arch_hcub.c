/* Copyright 2004,2007,2008,2010-2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : arch_hcub.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the hypercube       **/
/**                target architecture.                    **/
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
/**                # Version 3.1  : from : 11 jun 1996     **/
/**                                 to     11 jun 1996     **/
/**                # Version 3.2  : from : 21 sep 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 4.0  : from : 11 nov 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     14 fev 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_HCUB

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_hcub.h"

/********************************************/
/*                                          */
/* These are the binary hypercube routines. */
/*                                          */
/********************************************/

/* This routine loads the binary
** hypercube architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archHcubArchLoad (
ArchHcub * restrict const   archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchHcub)    > sizeof (ArchDummy)) ||
      (sizeof (ArchHcubDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archHcubArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &archptr->dimmax) != 1) ||
      (archptr->dimmax < 1)                     ||
      (archptr->dimmax > (sizeof (archptr->dimmax) << 3))) {
    errorPrint ("archHcubArchLoad: bad input");
    return     (1);
  }

  return (0);
}

/* This routine saves the
** binary hypercube architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archHcubArchSave (
const ArchHcub * const      archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchHcub)    > sizeof (ArchDummy)) ||
      (sizeof (ArchHcubDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archHcubArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " ", (Anum) archptr->dimmax) == EOF) {
    errorPrint ("archHcubArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archHcubDomNum (
const ArchHcub * const      archptr,
const ArchHcubDom * const   domptr)
{
  return (domptr->bitset);                        /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archHcubDomTerm (
const ArchHcub * const      archptr,
ArchHcubDom * const         domptr,
const ArchDomNum            domnum)
{
  if (domnum < (1 << archptr->dimmax)) {          /* If valid label */
    domptr->dimcur = 0;                           /* Set the domain */
    domptr->bitset = domnum;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the hypercube domain.
*/

Anum 
archHcubDomSize (
const ArchHcub * const      archptr,
const ArchHcubDom * const   domptr)
{
  return (1 << domptr->dimcur);
}

/* This function returns the average distance
** between two sub-hypercubes.
*/

Anum
archHcubDomDist (
const ArchHcub * const      archptr,
const ArchHcubDom * const   dom0ptr,
const ArchHcubDom * const   dom1ptr)
{
  Anum                i, j, k;

  if (dom0ptr->dimcur > dom1ptr->dimcur) {        /* Get smallest set dimension value */
    i = dom0ptr->dimcur;
    j = i - dom1ptr->dimcur;
  }
  else {
    i = dom1ptr->dimcur;
    j = i - dom0ptr->dimcur;
  }
  j /= 2;                                         /* For set/unset bits, assume 1/2 difference */

  for (k = (dom0ptr->bitset ^ dom1ptr->bitset) >> i, i = archptr->dimmax - i;
       i > 0;
       k >>= 1, i --)
    j += (k & 1);                                 /* Add Hamming difference on set dimensions */

  return (j);
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archHcubDomFrst (
const ArchHcub * const        archptr,
ArchHcubDom * restrict const  domptr)
{
  domptr->dimcur = archptr->dimmax;
  domptr->bitset = 0;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archHcubDomLoad (
const ArchHcub * const        archptr,
ArchHcubDom * restrict const  domptr,
FILE * restrict const         stream)
{
  if ((intLoad (stream, &domptr->dimcur) != 1) ||
      (intLoad (stream, &domptr->bitset) != 1) ||
      (domptr->dimcur > archptr->dimmax)) {
    errorPrint ("archHcubDomLoad: bad input");
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
archHcubDomSave (
const ArchHcub * const      archptr,
const ArchHcubDom * const   domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->dimcur,
               (Anum) domptr->bitset) == EOF) {
    errorPrint ("archHcubDomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a hypercube
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archHcubDomBipart (
const ArchHcub * const        archptr,
const ArchHcubDom * const     domptr,
ArchHcubDom * restrict const  dom0ptr,
ArchHcubDom * restrict const  dom1ptr)
{
  if (domptr->dimcur <= 0)                        /* Return if cannot bipartition more */
    return (1);

  dom0ptr->dimcur =
  dom1ptr->dimcur = domptr->dimcur - 1;
  dom0ptr->bitset = domptr->bitset;
  dom1ptr->bitset = domptr->bitset | (1 << dom1ptr->dimcur);

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
archHcubDomIncl (
const ArchHcub * const      archptr,
const ArchHcubDom * const   dom0ptr,
const ArchHcubDom * const   dom1ptr)
{
  if ((dom0ptr->dimcur >= dom1ptr->dimcur) &&
      (((dom0ptr->bitset ^ dom1ptr->bitset) >> dom0ptr->dimcur) == 0))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** hypercube domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archHcubDomMpiType (
const ArchHcub * const        archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (2, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
