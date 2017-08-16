/* Copyright 2004,2007,2008,2010-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_tleaf.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the tree-leaf       **/
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
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     20 jul 1996     **/
/**                # Version 3.2  : from : 10 oct 1996     **/
/**                                 to     14 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 07 jun 2001     **/
/**                                 to     07 jun 2001     **/
/**                # Version 4.0  : from : 10 dec 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 jul 2012     **/
/**                                                        **/
/**   NOTES      : # The ltleaf architecture was proposed  **/
/**                  by Emmanuel Jeannot and Francois      **/
/**                  Tessier.                              **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_TLEAF

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_tleaf.h"

/*******************************************/
/*                                         */
/* These are the tree-leaf graph routines. */
/*                                         */
/*******************************************/

/* This routine loads the
** tree leaf architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archTleafArchLoad (
ArchTleaf * restrict const  archptr,
FILE * restrict const       stream)
{
  Anum                sizeval;
  Anum                levlnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTleaf)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTleafDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTleafArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (intLoad (stream, &archptr->levlnbr) != 1) {
    errorPrint ("archTleafArchLoad: bad input (1)");
    return     (1);
  }

  if ((archptr->sizetab = memAlloc ((archptr->levlnbr * 2 + 1) * sizeof (Anum))) == NULL) { /* TRICK: One more slot for linktab[-1] */
    errorPrint ("archTleafArchLoad: out of memory");
    return     (1);
  }
  archptr->linktab     = archptr->sizetab + archptr->levlnbr + 1; /* TRICK: One more slot     */
  archptr->linktab[-1] = 0;                       /* Dummy slot for for level-0 communication */
  archptr->permtab     = NULL;                    /* Assume no permutation array              */

  for (levlnum = 0, sizeval = 1; levlnum < archptr->levlnbr; levlnum ++) {
    if ((intLoad (stream, &archptr->sizetab[levlnum]) != 1) ||
        (intLoad (stream, &archptr->linktab[levlnum]) != 1) ||
        (archptr->sizetab[levlnum] < 2)                     ||
        (archptr->linktab[levlnum] < 1)) {
      errorPrint ("archTleafArchLoad: bad input (2)");
      return (1);
    }
    sizeval *= archptr->sizetab[levlnum];
  }
  archptr->termnbr = sizeval;

  return (0);
}

/* This routine frees the tree
** leaf architecture structures.
** It returns:
** - 0   : if the architecture has been successfully freed.
** - !0  : on error.
*/

int
archTleafArchFree (
ArchTleaf * const           archptr)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTleaf)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTleafDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTleafArchFree: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  memFree (archptr->sizetab);                     /* Free group leader */
  if (archptr->permtab != NULL)
    memFree (archptr->permtab);                   /* Free group leader */

#ifdef SCOTCH_DEBUG_ARCH2
  archptr->sizetab =
  archptr->linktab =
  archptr->permtab = NULL;
#endif /* SCOTCH_DEBUG_ARCH2 */

  return (0);
}

/* This routine saves the
** tree leaf architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archTleafArchSave (
const ArchTleaf * const     archptr,
FILE * restrict const       stream)
{
  Anum                levlnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchTleaf)    > sizeof (ArchDummy)) ||
      (sizeof (ArchTleafDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archTleafArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING,
               (Anum) archptr->levlnbr) == EOF) {
    errorPrint ("archTleafArchSave: bad output (1)");
    return     (1);
  }

  for (levlnum = 0; levlnum < archptr->levlnbr; levlnum ++) {
    if (fprintf (stream, " " ANUMSTRING " " ANUMSTRING,
                 (Anum) archptr->sizetab[levlnum],
                 (Anum) archptr->linktab[levlnum]) == EOF) {
      errorPrint ("archTleafArchSave: bad output (2)");
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
archTleafDomNum (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  domnptr)
{
  Anum                levlnum;
  Anum                sizeval;

  sizeval = 1;                                    /* Compute size of blocks below */
  for (levlnum = domnptr->levlnum; levlnum < archptr->levlnbr; levlnum ++)
    sizeval *= archptr->sizetab[levlnum];

  return (domnptr->indxmin * sizeval);
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archTleafDomTerm (
const ArchTleaf * const     archptr,
ArchTleafDom * const        domnptr,
const ArchDomNum            domnnum)
{
#ifdef SCOTCH_DEBUG_ARCH2
  if (domnnum < 0) {
    errorPrint ("archTleafDomTerm: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH2 */

  if (domnnum < archptr->termnbr) {               /* If valid label */
    domnptr->levlnum = archptr->levlnbr;          /* Set the domain */
    domnptr->indxmin = domnnum;
    domnptr->indxnbr = 1;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the subtree domain.
*/

Anum 
archTleafDomSize (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  domnptr)
{
  Anum                levlnum;
  Anum                sizeval;

  sizeval = 1;                                    /* Compute size of blocks below */
  for (levlnum = domnptr->levlnum; levlnum < archptr->levlnbr; levlnum ++)
    sizeval *= archptr->sizetab[levlnum];

  return (sizeval * domnptr->indxnbr);
}

/* This function returns the average
** distance between two tree leaf
** subdomains.
*/

Anum 
archTleafDomDist (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  dom0ptr,
const ArchTleafDom * const  dom1ptr)
{
  Anum                lev0num;
  Anum                lev1num;
  Anum                idx0min;
  Anum                idx1min;
  Anum                idx0nbr;
  Anum                idx1nbr;
  Anum                distval;

  const Anum * const  sizetab = archptr->sizetab;

  lev0num = dom0ptr->levlnum;
  lev1num = dom1ptr->levlnum;
  idx0min = dom0ptr->indxmin;
  idx1min = dom1ptr->indxmin;
  idx0nbr = dom0ptr->indxnbr;
  idx1nbr = dom1ptr->indxnbr;

  if (lev0num != lev1num) {
    if (lev0num > lev1num) {
      idx0nbr = 1;
      do {
        lev0num --;
        idx0min /= sizetab[lev0num];
      } while (lev0num > lev1num);
    }
    else {
      idx1nbr = 1;
      do {
        lev1num --;
        idx1min /= sizetab[lev1num];
      } while (lev1num > lev0num);
    }
  }

  distval = archptr->linktab[lev0num - 1];        /* Get cost at this level */

  return (((idx0min >= (idx1min + idx1nbr)) ||    /* If inclusion, only half of the distance */
           (idx1min >= (idx0min + idx0nbr))) ? distval : ((idx0nbr == idx1nbr) ? 0 : (distval >> 1)));
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTleafDomFrst (
const ArchTleaf * const       archptr,
ArchTleafDom * restrict const domnptr)
{
  domnptr->levlnum = 0;
  domnptr->indxmin = 0;
  domnptr->indxnbr = 1;                           /* The root vertex is unique */

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archTleafDomLoad (
const ArchTleaf * const       archptr,
ArchTleafDom * restrict const domnptr,
FILE * const                  stream)
{
  if ((intLoad (stream, &domnptr->levlnum) != 1) ||
      (intLoad (stream, &domnptr->indxmin) != 1) ||
      (intLoad (stream, &domnptr->indxnbr) != 1) ||
      (domnptr->levlnum < 0)                     ||
      (domnptr->levlnum > archptr->levlnbr)) {
    errorPrint ("archTleafDomLoad: bad input");
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
archTleafDomSave (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  domnptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " " ANUMSTRING " ",
               (Anum) domnptr->levlnum,
               (Anum) domnptr->indxmin,
               (Anum) domnptr->indxnbr) == EOF) {
    errorPrint ("archTleafDomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a tree leaf
** domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archTleafDomBipart (
const ArchTleaf * const       archptr,
const ArchTleafDom * const    domnptr,
ArchTleafDom * restrict const dom0ptr,
ArchTleafDom * restrict const dom1ptr)
{
  Anum                sizeval;

  if (domnptr->indxnbr <= 1) {                    /* If dubdomain has only one node at this level */
    if (domnptr->levlnum >= archptr->levlnbr)     /* Return if cannot bipartition more            */
      return (1);

    sizeval = archptr->sizetab[domnptr->levlnum]; /* Partition all the vertices of a new level */

    dom0ptr->levlnum =
    dom1ptr->levlnum = domnptr->levlnum + 1;
    dom0ptr->indxmin = domnptr->indxmin * sizeval;
  }
  else {                                          /* Subdomain has several indices */
    sizeval = domnptr->indxnbr;                   /* Base on existing block size   */

    dom0ptr->levlnum =                            /* Stay at same level */
    dom1ptr->levlnum = domnptr->levlnum;
    dom0ptr->indxmin = domnptr->indxmin;          /* Start from the existing start index */
  }

  dom0ptr->indxnbr = (sizeval + 1) >> 1;          /* Subdomain 0 is always the largest one */
  dom1ptr->indxmin = dom0ptr->indxmin + dom0ptr->indxnbr;
  dom1ptr->indxnbr = sizeval - dom0ptr->indxnbr;

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
archTleafDomIncl (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  dom0ptr,
const ArchTleafDom * const  dom1ptr)
{
  Anum                lev0num;
  Anum                lev1num;
  Anum                idx0min;
  Anum                idx1min;
  Anum                idx0nbr;
  Anum                idx1nbr;

  const Anum * const  sizetab = archptr->sizetab;

  lev0num = dom0ptr->levlnum;
  lev1num = dom1ptr->levlnum;
  idx0min = dom0ptr->indxmin;
  idx1min = dom1ptr->indxmin;
  idx0nbr = dom0ptr->indxnbr;
  idx1nbr = dom1ptr->indxnbr;

  if (lev0num != lev1num) {
    if (lev1num > lev0num) {
      idx1nbr = 1;
      do {
        lev1num --;
        idx1min /= sizetab[lev1num];
      } while (lev1num > lev0num);
    }
    else
      return (0);
  }

  return (((idx0min >= (idx1min + idx1nbr)) ||
           (idx1min >= (idx0min + idx0nbr))) ? 0 : 1);
}

/* This function creates the MPI_Datatype for
** tree-leaf domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archTleafDomMpiType (
const ArchTleaf * const       archptr,
MPI_Datatype * const          typeptr)
{
  MPI_Type_contiguous (3, ANUM_MPI, typeptr);

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */

/***********************************/
/*                                 */
/* These are the labeled tree-leaf */
/* graph routines.                 */
/*                                 */
/***********************************/

/* This routine loads the labeled
** tree leaf architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archLtleafArchLoad (
ArchTleaf * restrict const  archptr,
FILE * restrict const       stream)
{
  Anum                sizeval;
  Anum                levlnum;
  Anum                permnum;

  if (archTleafArchLoad (archptr, stream) != 0)   /* Read tree part */
    return (1);

  if ((intLoad (stream, &archptr->permnbr) != 1) ||
      (archptr->permnbr <= 0)) {
    errorPrint ("archLtleafArchLoad: bad input (1)");
    return     (1);
  }

#ifdef SCOTCH_DEBUG_ARCH2
  if (archptr->permnbr != 1) {                    /* Valid empty permutation is of size 1 */
    for (levlnum = archptr->levlnbr - 1, sizeval = archptr->sizetab[levlnum];
         sizeval != archptr->permnbr; levlnum --, sizeval *= archptr->sizetab[levlnum]) {
      if (levlnum < 0) {
        errorPrint ("archLtleafArchLoad: permutation size does not match level boundaries");
        return     (1);
      }
    }
  }
#endif /* SCOTCH_DEBUG_ARCH2 */

  if ((archptr->permtab = memAlloc (archptr->permnbr * 2 * sizeof (Anum))) == NULL) { /* TRICK: space for peritab too */
    errorPrint ("archLtleafArchLoad: out of memory");
    return     (1);
  }

  for (permnum = 0; permnum < archptr->permnbr; permnum ++) {
#ifdef SCOTCH_DEBUG_ARCH2
    Anum                permtmp;
#endif /* SCOTCH_DEBUG_ARCH2 */

    if ((intLoad (stream, &archptr->permtab[permnum]) != 1) ||
        (archptr->permtab[permnum] < 0)                     ||
        (archptr->permtab[permnum] >= archptr->permnbr)) {
      errorPrint ("archLtleafArchLoad: bad input (2)");
      return (1);
    }
#ifdef SCOTCH_DEBUG_ARCH2
    for (permtmp = 0; permtmp < permnum; permtmp ++) {
      if (archptr->permtab[permtmp] == archptr->permtab[permnum]) {
        errorPrint ("archLtleafArchLoad: duplicate permutation index");
        return     (1);
      }
    }
#endif /* SCOTCH_DEBUG_ARCH2 */
  }

  archptr->peritab = archptr->permtab + archptr->permnbr;
  for (permnum = 0; permnum < archptr->permnbr; permnum ++) /* Build inverse permutation */
    archptr->peritab[archptr->permtab[permnum]] = permnum;

  return (0);
}

/* This routine saves the labeled
** tree leaf architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archLtleafArchSave (
const ArchTleaf * const     archptr,
FILE * restrict const       stream)
{
  Anum                permnum;

  if (archTleafArchSave (archptr, stream) != 0)   /* Save tree part */
    return (1);

  if (fprintf (stream, ANUMSTRING,
               (Anum) archptr->permnbr) == EOF) {
    errorPrint ("archLtleafArchSave: bad output (1)");
    return     (1);
  }

  for (permnum = 0; permnum < archptr->permnbr; permnum ++) {
    if (fprintf (stream, " " ANUMSTRING,
                 (Anum) archptr->permtab[permnum]) == EOF) {
      errorPrint ("archLtleafArchSave: bad output (2)");
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
archLtleafDomNum (
const ArchTleaf * const     archptr,
const ArchTleafDom * const  domnptr)
{
  Anum                levlnum;
  Anum                sizeval;
  Anum                domnnum;
  Anum                permnum;

  sizeval = 1;                                    /* Compute size of blocks below */
  for (levlnum = domnptr->levlnum; levlnum < archptr->levlnbr; levlnum ++)
    sizeval *= archptr->sizetab[levlnum];

  domnnum = domnptr->indxmin * sizeval;
  permnum = domnnum % archptr->permnbr;           /* Get non permuted index as terminal domain */

  return (domnnum - permnum + archptr->permtab[permnum]); /* Return permuted index */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archLtleafDomTerm (
const ArchTleaf * const     archptr,
ArchTleafDom * const        domnptr,
const ArchDomNum            domnnum)
{
#ifdef SCOTCH_DEBUG_ARCH2
  if (domnnum < 0) {
    errorPrint ("archLtleafDomTerm: invalid parameter");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH2 */

  if (domnnum < archptr->termnbr) {               /* If valid label */
    Anum                permnum;

    permnum = domnnum % archptr->permnbr;         /* Get permuted index as terminal domain */

    domnptr->levlnum = archptr->levlnbr;          /* Set the domain */
    domnptr->indxmin = domnnum - permnum + archptr->peritab[permnum];
    domnptr->indxnbr = 1;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}
