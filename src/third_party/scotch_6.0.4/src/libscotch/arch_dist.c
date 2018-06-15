/* Copyright 2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_dist.c                             **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the distance        **/
/**                multiplicator pseudo-architecture       **/
/**                functions. This pseudo-architecture is  **/
/**                used by graph repartitioning routines   **/
/**                to handle floating-point migration      **/
/**                costs.                                  **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 14 fev 2011     **/
/**                                 to   : 30 jun 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_DIST

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_dist.h"

/**************************************/
/*                                    */
/* These are the entry points for the */
/* distance graph routines. They are  */
/* used only in debugging mode, to    */
/* provide breakpoints for routines   */
/* which are else implemented as      */
/* macros for the sake of efficiency. */
/*                                    */
/**************************************/

/* This routine loads the distance
** graph architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archDistArchLoad (
ArchDist * restrict const   archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if (sizeof (ArchDist) > sizeof (ArchDummy)) {
    errorPrint ("archDistArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (intLoad (stream, &archptr->crloval) != 1) {
    errorPrint ("archDistArchLoad: bad input");
    return     (1);
  }

  return (archLoad (archptr->archptr, stream));   /* Load sub-architecture */
}

/* This routine saves the
** distance graph architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archDistArchSave (
const ArchDist * const      archptr,
FILE * restrict const       stream)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if (sizeof (ArchDist) > sizeof (ArchDummy)) {
    errorPrint ("archDistArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING "\t",
               (Anum) archptr->crloval) == EOF) {
    errorPrint ("archDistArchSave: bad output");
    return     (1);
  }

  return (archSave (archptr->archptr, stream));   /* Save sub-architecture */
}

/* This routine build the
** distance graph architecture of
** an original one.
** It returns:
** - 0   : if the architecture has been successfully built.
** - !0  : on error.
*/

int
archDistArchBuild (
Arch * const                archptr,
Arch * const                orgarchptr,
const Anum                  crloval)
{
  ArchDist *          archdataptr;

  archInit (archptr);                             /* Initialize architecture body */
  archptr->class   = archClass ("dist");          /* Set architecture class       */
  archptr->flagval = orgarchptr->flagval;         /* Set architecture flag        */
  archdataptr = (ArchDist *) (void *) &archptr->data;
  archdataptr->archptr = orgarchptr;
  archdataptr->crloval = crloval;

  return (0); 
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archDistDomNum (
const ArchDist * const      archptr,
const ArchDom * const       domptr)
{
  return (archDomNum (archptr->archptr, domptr)); /* Call proper routine */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archDistDomTerm (
const ArchDist * const      archptr,
ArchDom * const             domptr,
const ArchDomNum            domnum)
{
  return (archDomTerm (archptr->archptr, domptr, domnum)); /* Call proper routine */
}

/* This function returns the number of
** elements in the distance domain.
*/

Anum 
archDistDomSize (
const ArchDist * const      archptr,
const ArchDom * const       domptr)
{
  return (archDomSize (archptr->archptr, domptr)); /* Call proper routine */
}

/* This function returns the weight of
** the given distance domain.
*/

Anum 
archDistDomWght (
const ArchDist * const      archptr,
const ArchDom * const       domptr)
{
  return (archDomWght (archptr->archptr, domptr)); /* Call proper routine */
}

/* This function returns the average
** distance between two distance
** subdomains.
*/

Anum 
archDistDomDist (
const ArchDist * const      archptr,
const ArchDom * const       dom0ptr,
const ArchDom * const       dom1ptr)
{
  return (archptr->crloval * archDomDist (archptr->archptr, dom0ptr, dom1ptr));
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archDistDomFrst (
const ArchDist * const      archptr,
ArchDom * restrict const    domptr)
{
  return (archDomFrst (archptr->archptr, domptr)); /* Call proper routine */
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archDistDomLoad (
const ArchDist * const      archptr,
ArchDom * restrict const    domptr,
FILE * const                stream)
{
  return (archDomLoad (archptr->archptr, domptr, stream)); /* Call proper routine */
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archDistDomSave (
const ArchDist * const      archptr,
const ArchDom * const       domptr,
FILE * const                stream)
{
  return (archDomSave (archptr->archptr, domptr, stream)); /* Call proper routine */
}

/* This function tries to split a distance
** graph domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archDistDomBipart (
const ArchDist * const      archptr,
const ArchDom * const       domptr,
ArchDom * restrict const    dom0ptr,
ArchDom * restrict const    dom1ptr)
{
  return (archDomBipart (archptr->archptr, domptr, dom0ptr, dom1ptr)); /* Call proper routine */
}

/* This function checks if dom1 is
** included in dom0.
** It returns:
** - 0  : if dom1 is not included in dom0.
** - 1  : if dom1 is included in dom0.
** - 2  : on error.
*/

int
archDistDomIncl (
const ArchDist * const      archptr,
const ArchDom * const       dom0ptr,
const ArchDom * const       dom1ptr)
{
  return (archDomIncl (archptr->archptr, dom0ptr, dom1ptr)); /* Call proper routine */
}

/* This function creates the MPI_Datatype for
** distance graph domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archDistDomMpiType (
const ArchDist * const        archptr,
MPI_Datatype * const          typeptr)
{
  return (archDomMpiType (archptr->archptr, typeptr)); /* Call proper routine as we don't add any parameter */
}
#endif /* SCOTCH_PTSCOTCH */

