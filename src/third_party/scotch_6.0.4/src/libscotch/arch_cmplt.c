/* Copyright 2004,2007,2010,2011 ENSEIRB, INRIA & CNRS
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
/**   NAME       : arch_cmplt.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the complete graph  **/
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
/**                                 to     13 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 4.0  : from : 09 jan 2004     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 19 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     14 fev 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_CMPLT

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_cmplt.h"

/******************************************/
/*                                        */
/* These are the complete graph routines. */
/*                                        */
/******************************************/

/* This routine loads the complete
** graph architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archCmpltArchLoad (
ArchCmplt * restrict const  archptr,
FILE * restrict const       stream)
{
  long                numnbr;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmplt)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((fscanf (stream, "%ld", &numnbr) != 1) ||
      (numnbr < 1)) {
    errorPrint ("archCmpltArchLoad: bad input");
    return     (1);
  }
  archptr->numnbr = (Anum) numnbr;

  return (0);
}

/* This routine saves the
** complete graph architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archCmpltArchSave (
const ArchCmplt * const     archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmplt)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING " ", (Anum) archptr->numnbr) == EOF) {
    errorPrint ("archCmpltArchSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archCmpltDomNum (
const ArchCmplt * const     archptr,
const ArchCmpltDom * const  domptr)
{
  return (domptr->nummin);                        /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archCmpltDomTerm (
const ArchCmplt * const     archptr,
ArchCmpltDom * const        domptr,
const ArchDomNum            domnum)
{
  if (domnum < archptr->numnbr) {                 /* If valid label */
    domptr->nummin = domnum;                      /* Set the domain */
    domptr->numnbr = 1;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the complete domain.
*/

Anum 
archCmpltDomSize (
const ArchCmplt * const     archptr,
const ArchCmpltDom * const  domptr)
{
  return (domptr->numnbr);
}

/* This function returns the average
** distance between two complete
** subdomains.
*/

Anum 
archCmpltDomDist (
const ArchCmplt * const     archptr,
const ArchCmpltDom * const  dom0ptr,
const ArchCmpltDom * const  dom1ptr)
{
  return (((dom0ptr->nummin == dom1ptr->nummin) && /* All domains are at distance 1 */
           (dom0ptr->numnbr == dom1ptr->numnbr)) ? 0 : 1); /* If they are different */
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltDomFrst (
const ArchCmplt * const         archptr,
ArchCmpltDom * restrict const   domptr)
{
  domptr->nummin = 0;
  domptr->numnbr = archptr->numnbr;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltDomLoad (
const ArchCmplt * const       archptr,
ArchCmpltDom * restrict const domptr,
FILE * const                  stream)
{
  long                nummin;
  long                numnbr;

  if ((fscanf (stream, "%ld%ld",
               &nummin,
               &numnbr) != 2) ||
      (numnbr < 1)            ||
      (numnbr + nummin > (long) archptr->numnbr)) {
    errorPrint ("archCmpltDomLoad: bad input");
    return     (1);
  }
  domptr->nummin = (Anum) nummin;
  domptr->numnbr = (Anum) numnbr;

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltDomSave (
const ArchCmplt * const     archptr,
const ArchCmpltDom * const  domptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->nummin,
               (Anum) domptr->numnbr) == EOF) {
    errorPrint ("archCmpltDomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a complete
** graph domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archCmpltDomBipart (
const ArchCmplt * const       archptr,
const ArchCmpltDom * const    domptr,
ArchCmpltDom * restrict const dom0ptr,
ArchCmpltDom * restrict const dom1ptr)
{
  if (domptr->numnbr <= 1)                        /* Return if cannot bipartition more */
    return (1);

  dom0ptr->nummin = domptr->nummin;               /* Bipartition vertices */
  dom0ptr->numnbr = domptr->numnbr / 2;
  dom1ptr->nummin = domptr->nummin + dom0ptr->numnbr;
  dom1ptr->numnbr = domptr->numnbr - dom0ptr->numnbr;

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
archCmpltDomIncl (
const ArchCmplt * const     archptr,
const ArchCmpltDom * const  dom0ptr,
const ArchCmpltDom * const  dom1ptr)
{
  if ((dom1ptr->nummin >= dom0ptr->nummin) &&
      ((dom1ptr->nummin + dom1ptr->numnbr) <= (dom0ptr->nummin + dom0ptr->numnbr)))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** complete graph domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archCmpltDomMpiType (
const ArchCmplt * const       archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (2, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
