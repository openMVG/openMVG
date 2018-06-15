/* Copyright 2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_cmpltw.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the weighted        **/
/**                complete graph target architecture.     **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 dec 2007     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     23 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_CMPLTW

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "arch_cmpltw.h"

/******************************************/
/*                                        */
/* These are the complete graph routines. */
/*                                        */
/******************************************/

/* This routine builds a complete weighted
** graph architecture from the given load array.
** It returns:
** - 0   : if the architecture has been successfully built.
** - !0  : on error.
*/

static
void
archCmpltwArchBuild3 (
ArchCmpltwLoad * restrict const velotab,
ArchCmpltwLoad * restrict const vesotab,
Anum                            vertnbr,
Anum                            velosum)
{
  Anum                velosum0;
  Anum                velosum1;
  Anum                vertnbr0;
  Anum                vertnbr1;
  Anum                vertnum0;
  Anum                vertnum1;
  Anum                vertnum;

  vertnum0 =
  vertnum1 = vertnbr - 1;
  velosum0 = velotab[vertnum0 --].veloval;
  velosum1 = 0;
  for (vertnum = vertnum0; vertnum >= 0; vertnum --) {
    if (velosum1 < velosum0) {
      velosum1            += velotab[vertnum].veloval;
      vesotab[vertnum1 --] = velotab[vertnum];
    }
    else {
      velosum0            += velotab[vertnum].veloval;
      velotab[vertnum0 --] = velotab[vertnum];
    }
  }

  if (velosum0 >= velosum1) {
    vertnbr0 = vertnbr - vertnum0 - 1;
    vertnbr1 = vertnbr - vertnbr0;
    memMov (velotab,            velotab + vertnbr1, vertnbr0 * sizeof (ArchCmpltwLoad));
    memCpy (velotab + vertnbr0, vesotab + vertnbr0, vertnbr1 * sizeof (ArchCmpltwLoad));
  }
  else {
    Anum                velotmp;

    vertnbr0 = vertnbr - vertnum1 - 1;
    vertnbr1 = vertnbr - vertnbr0;
    memCpy (velotab, vesotab + vertnbr1, vertnbr0 * sizeof (ArchCmpltwLoad));
    velotmp  = velosum0;
    velosum0 = velosum1;
    velosum1 = velotmp;
  }

  if (vertnbr0 > 2)
    archCmpltwArchBuild3 (velotab, vesotab, vertnbr0, velosum0);
  if (vertnbr1 > 2)
    archCmpltwArchBuild3 (velotab + vertnbr0, vesotab + vertnbr0, vertnbr1, velosum1);
}

static
int
archCmpltwArchBuild2 (
ArchCmpltw * restrict const archptr)
{
  ArchCmpltwLoad * restrict vesotab;              /* Auxiliary sort array for weighted vertices */

  if (archptr->vertnbr < 3)                       /* No need to sort if less than 3 vertices */
    return (0);

  if ((vesotab = (ArchCmpltwLoad *) memAlloc (archptr->vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchBuild2: out of memory");
    memFree (archptr->velotab);
    archptr->velotab = NULL;
    return (1);
  }

  intSort2asc2 (archptr->velotab, archptr->vertnbr); /* Sort load array by both keys to be portable across sorting implementations */
  
  archCmpltwArchBuild3 (archptr->velotab, vesotab, archptr->vertnbr, archptr->velosum);

  memFree (vesotab);

  return (0);
}

int
archCmpltwArchBuild (
ArchCmpltw * restrict const archptr,
const Anum                  vertnbr,
const Anum * restrict const velotab)
{
  Anum                vertnum;
  Anum                velosum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchBuild: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (vertnbr <= 0) {
    errorPrint ("archCmpltwArchBuild: invalid parameters");
    return     (1);
  }

  archptr->vertnbr = (Anum) vertnbr;

  if ((archptr->velotab = (ArchCmpltwLoad *) memAlloc (archptr->vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchBuild: out of memory");
    return     (1);
  }

  for (vertnum = 0, velosum = 0; vertnum < archptr->vertnbr; vertnum ++) { /* Fill vertex load array */
    Anum                veloval;

    veloval  = velotab[vertnum];
    velosum += veloval;
    archptr->velotab[vertnum].veloval = veloval;
    archptr->velotab[vertnum].vertnum = vertnum;
  }
  archptr->velosum = (Anum) velosum;

  return (archCmpltwArchBuild2 (archptr));
}

/* This routine loads the weighted complete
** graph architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archCmpltwArchLoad (
ArchCmpltw * restrict const  archptr,
FILE * restrict const       stream)
{
  long                vertnbr;
  Anum                velosum;
  Anum                vertnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((fscanf (stream, "%ld", &vertnbr) != 1) ||
      (vertnbr < 1)) {
    errorPrint ("archCmpltwArchLoad: bad input (1)");
    return     (1);
  }
  archptr->vertnbr = (Anum) vertnbr;

  if ((archptr->velotab = (ArchCmpltwLoad *) memAlloc (archptr->vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchLoad: out of memory");
    return     (1);
  }

  for (vertnum = 0, velosum = 0; vertnum < archptr->vertnbr; vertnum ++) {
    long                veloval;
    Anum                velotmp;

    if ((fscanf (stream, "%ld", &veloval) != 1) ||
        (veloval < 1)) {
      errorPrint ("archCmpltwArchLoad: bad input (2)");
      return     (1);
    }

    velotmp  = (Anum) veloval;
    velosum += velotmp;
    archptr->velotab[vertnum].veloval = velotmp;
    archptr->velotab[vertnum].vertnum = vertnum;
  }
  archptr->velosum = (Anum) velosum;

  return (archCmpltwArchBuild2 (archptr));
}

/* This routine saves the weighted
** complete graph architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archCmpltwArchSave (
const ArchCmpltw * const    archptr,
FILE * restrict const       stream)
{
  Anum                vertnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING, (Anum) archptr->vertnbr) == EOF) {
    errorPrint ("archCmpltwArchSave: bad output (1)");
    return     (1);
  }

  for (vertnum = 0; vertnum < archptr->vertnbr; vertnum ++) { /* For all weights to output */
    Anum                verttmp;

    for (verttmp = 0; verttmp < archptr->vertnbr; verttmp ++) { /* For all vertex indices: O(n^2) loop but we don't really care */
      if (archptr->velotab[verttmp].vertnum == vertnum) {
        if (fprintf (stream, " " ANUMSTRING, (Anum) archptr->velotab[verttmp].veloval) == EOF) {
          errorPrint ("archCmpltwArchSave: bad output (2)");
          return     (1);
        }
        break;
      }
      if (verttmp == archptr->vertnbr) {
        errorPrint ("archCmpltwArchSave: internal error");
        return     (1);
      }
    }
  }

  return (0);
}

/* This routine frees the weighted complete
** graph architecture data structures.
** It returns:
** - 0   : if the architecture has been successfully freed.
** - !0  : on error.
*/

int
archCmpltwArchFree (
ArchCmpltw * const          archptr)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchFree: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (archptr->velotab != NULL) {
    memFree (archptr->velotab);
    archptr->velotab = NULL;
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archCmpltwDomNum (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domptr)
{
  return (archptr->velotab[domptr->vertmin].vertnum); /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** 
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archCmpltwDomTerm (
const ArchCmpltw * const    archptr,
ArchCmpltwDom * const       domptr,
const ArchDomNum            domnum)
{
  if (domnum < archptr->vertnbr) {                /* If valid label */
    Anum                vertnum;

    for (vertnum = 0; vertnum < archptr->vertnbr; vertnum ++) { /* Search for terminal domain index matching vertex label */
      if (archptr->velotab[vertnum].vertnum == domnum)
        break;
    }
#ifdef SCOTCH_DEBUG_ARCH2
    if (vertnum == archptr->vertnbr) {            /* If index not found */
      errorPrint ("archCmpltwDomTerm: internal error");
      return     (2);
    }
#endif /* SCOTCH_DEBUG_ARCH2 */

    domptr->vertmin = vertnum;                    /* Set the domain */
    domptr->vertnbr = 1;
    domptr->veloval = archptr->velotab[vertnum].veloval;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the complete domain.
*/

Anum 
archCmpltwDomSize (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domptr)
{
  return (domptr->vertnbr);
}

/* This function returns the weight of
** the complete domain.
*/

Anum 
archCmpltwDomWght (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domptr)
{
  return (domptr->veloval);
}

/* This function returns the average
** distance between two complete
** subdomains.
*/

Anum 
archCmpltwDomDist (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const dom0ptr,
const ArchCmpltwDom * const dom1ptr)
{
  return (((dom0ptr->vertmin == dom1ptr->vertmin) && /* All domains are at distance 1 */
           (dom0ptr->vertnbr == dom1ptr->vertnbr)) ? 0 : 1); /* If they are different */
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomFrst (
const ArchCmpltw * const        archptr,
ArchCmpltwDom * restrict const  domptr)
{
  domptr->vertmin = 0;
  domptr->vertnbr = archptr->vertnbr;
  domptr->veloval = archptr->velosum;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomLoad (
const ArchCmpltw * const        archptr,
ArchCmpltwDom * restrict const  domptr,
FILE * const                    stream)
{
  long                vertmin;
  long                vertnbr;
  Anum                vertnum;
  Anum                vertnnd;
  Anum                velosum;

  if ((fscanf (stream, "%ld%ld",
               &vertmin,
               &vertnbr) != 2) ||
      (vertnbr < 1)            ||
      (vertnbr + vertmin > (long) archptr->vertnbr)) {
    errorPrint ("archCmpltwDomLoad: bad input");
    return     (1);
  }
  domptr->vertmin = (Anum) vertmin;
  domptr->vertnbr = (Anum) vertnbr;

  for (vertnum = domptr->vertmin, vertnnd = vertnum + domptr->vertnbr, velosum = 0;
       vertnum < vertnnd; vertnum ++)
    velosum += archptr->velotab[vertnum].veloval;

  domptr->veloval += velosum;

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomSave (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) domptr->vertmin,
               (Anum) domptr->vertnbr) == EOF) {
    errorPrint ("archCmpltwDomSave: bad output");
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
archCmpltwDomBipart (
const ArchCmpltw * const        archptr,
const ArchCmpltwDom * const     domptr,
ArchCmpltwDom * restrict const  dom0ptr,
ArchCmpltwDom * restrict const  dom1ptr)
{
  Anum                vertnum;
  Anum                velosum1;
  Anum                velosum2;                   /* Half of overall load sum */

  if (domptr->vertnbr <= 1)                       /* Return if cannot bipartition more */
    return (1);

  vertnum  = domptr->vertmin + domptr->vertnbr - 1;
  velosum1 = (Anum) archptr->velotab[vertnum].veloval;
  velosum2 = domptr->veloval / 2;
  for (vertnum --; vertnum > domptr->vertmin; vertnum --) {
    Anum                velotmp;

    velotmp = velosum1 + (Anum) archptr->velotab[vertnum].veloval;
    if (velotmp > velosum2)                       /* Domain 1 is always the least loaded */
      break;
    velosum1 = velotmp;
  }

  dom0ptr->vertmin = domptr->vertmin;             /* Bipartition vertices */
  dom1ptr->vertmin = vertnum + 1;
  dom0ptr->vertnbr = dom1ptr->vertmin - domptr->vertmin;
  dom1ptr->vertnbr = domptr->vertnbr - dom0ptr->vertnbr;
  dom0ptr->veloval = domptr->veloval - velosum1;
  dom1ptr->veloval = velosum1;

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
archCmpltwDomIncl (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const dom0ptr,
const ArchCmpltwDom * const dom1ptr)
{
  if ((dom1ptr->vertmin >= dom0ptr->vertmin) &&
      ((dom1ptr->vertmin + dom1ptr->vertnbr) <= (dom0ptr->vertmin + dom0ptr->vertnbr)))
    return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** weighted complete graph domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archCmpltwDomMpiType (
const ArchCmpltw * const      archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (3, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
