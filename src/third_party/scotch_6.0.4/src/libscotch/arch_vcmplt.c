/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_vcmplt.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the variable-sized  **/
/**                complete graph target architecture.     **/
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
/**                                 to     16 aug 1995     **/
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     20 jul 1996     **/
/**                # Version 3.2  : from : 15 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 14 sep 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 05 nov 2003     **/
/**                                 to     05 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     26 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_VCMPLT

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_vcmplt.h"

/*****************************************/
/*                                       */
/* These are the variable-sized complete */
/* graph handling routines.              */
/*                                       */
/*****************************************/

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archVcmpltDomNum (
const ArchVcmplt * const    archptr,
const ArchVcmpltDom * const domptr)
{
  return (domptr->termnum);                       /* Return terminal number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archVcmpltDomTerm (
const ArchVcmplt * const    archptr,
ArchVcmpltDom * const       domptr,
const ArchDomNum            domnum)
{
  Anum                termnum;
  Anum                termlvl;

  if (domnum != ARCHDOMNOTTERM) {                 /* If valid label     */
    if (domnum == 0)                              /* Not a legal domain */
      return (2);

    domptr->termnum = domnum;                     /* Set the domain */
    for (termnum = domnum, termlvl = 0; termnum > 1; termnum >>= 1, termlvl ++) ; /* Compute level */
    domptr->termlvl = termlvl;                    /* Set level */

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the domain.
*/

Anum 
archVcmpltDomSize (
const ArchVcmplt * const    archptr,
const ArchVcmpltDom * const domptr)
{
  return (1);                                     /* All domains have same size for bipartitioning */
}

/* This function returns the average
** distance between two subdomains.
*/

Anum 
archVcmpltDomDist (
const ArchVcmplt * const    archptr,
const ArchVcmpltDom * const dom0ptr,
const ArchVcmpltDom * const dom1ptr)
{
  return ((dom0ptr->termnum == dom1ptr->termnum) ? 0 : 1); /* All distinct terminals are at distance 1 */
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archVcmpltDomFrst (
const ArchVcmplt * const        archptr,
ArchVcmpltDom * restrict const  domptr)
{
  domptr->termlvl = 0;                            /* First terminal number */
  domptr->termnum = 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archVcmpltDomLoad (
const ArchVcmplt * const        archptr,
ArchVcmpltDom * restrict const  domptr,
FILE * const                    stream)
{
  Anum                termnum;
  Anum                termlvl;

  if (intLoad (stream, &domptr->termnum) != 1) {
    errorPrint ("archVcmpltDomLoad: bad input");
    return     (1);
  }

  for (termnum = domptr->termnum, termlvl = 0; termnum > 1; termnum >>= 1, termlvl ++) ; /* Compute level */
  domptr->termlvl = termlvl;

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archVcmpltDomSave (
const ArchVcmplt * const    archptr,
const ArchVcmpltDom * const domptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " ",
               (Anum) domptr->termnum) == EOF) {
    errorPrint ("archVcmpltDomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function splits a domain
** into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 2  : on error.
*/

int
archVcmpltDomBipart (
const ArchVcmplt * const        archptr,
const ArchVcmpltDom * const     domptr,
ArchVcmpltDom * restrict const  dom0ptr,
ArchVcmpltDom * restrict const  dom1ptr)
{
  dom0ptr->termlvl =                              /* Bipartition the domain */
  dom1ptr->termlvl = domptr->termlvl + 1;
  dom0ptr->termnum = domptr->termnum << 1;
  dom1ptr->termnum = dom0ptr->termnum + 1;

  return ((dom1ptr->termnum < domptr->termnum) ? 2 : 0); /* Return error on overflow */
}

/* This function checks if dom1 is
** included in dom0.
** It returns:
** - 0  : if dom1 is not included in dom0.
** - 1  : if dom1 is included in dom0.
** - 2  : on error.
*/

int
archVcmpltDomIncl (
const ArchVcmplt * const    archptr,
const ArchVcmpltDom * const dom0ptr,
const ArchVcmpltDom * const dom1ptr)
{
  if ((dom1ptr->termlvl >= dom0ptr->termlvl) &&
      ((dom1ptr->termnum >> (dom1ptr->termlvl - dom0ptr->termlvl)) == dom0ptr->termnum))
      return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** variable-sized complete graph domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archVcmpltDomMpiType (
const ArchVcmplt * const      archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (2, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
