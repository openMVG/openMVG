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
/**   NAME       : arch_vhcub.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the variable-sized  **/
/**                hypercube target architecture.          **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 08 nov 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 04 nov 2003     **/
/**                                 to     04 nov 2003     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     27 feb 2008     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     26 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_VHCUB

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_vhcub.h"

/********************************/
/*                              */
/* These are the variable-sized */
/* hypercube handling routines. */
/*                              */
/********************************/

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archVhcubDomNum (
const ArchVhcub * const     archptr,
const ArchVhcubDom * const  domptr)
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
archVhcubDomTerm (
const ArchVhcub * const     archptr,
ArchVhcubDom * const        domptr,
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
archVhcubDomSize (
const ArchVhcub * const     archptr,
const ArchVhcubDom * const  domptr)
{
  return (1);                                     /* All domains have same size for bipartitioning */
}

/* This function returns the average
** distance between two subdomains.
*/

Anum 
archVhcubDomDist (
const ArchVhcub * const     archptr,
const ArchVhcubDom * const  dom0ptr,
const ArchVhcubDom * const  dom1ptr)
{
  Anum                dom0num;
  Anum                dom1num;
  Anum                distval;

  if (dom0ptr->termlvl > dom1ptr->termlvl) {
    dom0num = dom0ptr->termnum >> (dom0ptr->termlvl - dom1ptr->termlvl);
    dom1num = dom1ptr->termnum;
    distval = (dom0ptr->termlvl - dom1ptr->termlvl) >> 1; /* One half of unknown bits */
  }
  else {
    dom0num = dom0ptr->termnum;
    dom1num = dom1ptr->termnum >> (dom1ptr->termlvl - dom0ptr->termlvl);
    distval = (dom1ptr->termlvl - dom0ptr->termlvl) >> 1; /* One half of unknown bits */
  }

  for (dom0num ^= dom1num; dom0num != 0;          /* Compute Hamming distance */
       distval += (dom0num & 1), dom0num >>= 1) ;

  return (distval);
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archVhcubDomFrst (
const ArchVhcub * const       archptr,
ArchVhcubDom * restrict const domptr)
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
archVhcubDomLoad (
const ArchVhcub * const       archptr,
ArchVhcubDom * restrict const domptr,
FILE * const                  stream)
{
  Anum                termnum;
  Anum                termlvl;

  if (intLoad (stream, &domptr->termnum) != 1) {
    errorPrint ("archVhcubDomLoad: bad input");
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
archVhcubDomSave (
const ArchVhcub * const     archptr,
const ArchVhcubDom * const  domptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " ",
               (Anum) domptr->termnum) == EOF) {
    errorPrint ("archVhcubDomSave: bad output");
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
archVhcubDomBipart (
const ArchVhcub * const       archptr,
const ArchVhcubDom * const    domptr,
ArchVhcubDom * restrict const dom0ptr,
ArchVhcubDom * restrict const dom1ptr)
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
archVhcubDomIncl (
const ArchVhcub * const     archptr,
const ArchVhcubDom * const  dom0ptr,
const ArchVhcubDom * const  dom1ptr)
{
  if ((dom1ptr->termlvl >= dom0ptr->termlvl) &&
      ((dom1ptr->termnum >> (dom1ptr->termlvl - dom0ptr->termlvl)) == dom0ptr->termnum))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** variable-sized hypercube domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archVhcubDomMpiType (
const ArchVhcub * const       archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (2, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
