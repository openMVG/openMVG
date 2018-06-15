/* Copyright 2004,2007-2009,2011,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : mapping.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles (partial) mappings. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 31 mar 1993     **/
/**                                 to     31 mar 1993     **/
/**                # Version 1.0  : from : 04 oct 1993     **/
/**                                 to     06 oct 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 09 apr 1994     **/
/**                                 to     11 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     19 oct 1995     **/
/**                # Version 3.1  : from : 30 oct 1995     **/
/**                                 to     14 jun 1996     **/
/**                # Version 3.2  : from : 23 aug 1996     **/
/**                                 to     07 sep 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to     30 mar 1999     **/
/**                # Version 3.4  : from : 11 sep 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to     05 jan 2005     **/
/**                # Version 5.1  : from : 25 jun 2008     **/
/**                                 to     28 apr 2009     **/
/**                # Version 6.0  : from : 04 mar 2011     **/
/**                                 to     16 sep 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MAPPING

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "mapping.h"

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine builds a mapping.
** It returns:
** - 0   : if mapping successfully initialized.
** - !0  : on error.
*/

void
mapInit (
Mapping * restrict const        mappptr,          /*+ Mapping structure                  +*/
const Graph * restrict const    grafptr,          /*+ Graph data                         +*/
const Arch * restrict const     archptr,          /*+ Architecture data                  +*/
const ArchDom * restrict const  domnptr)          /*+ Target architecture initial domain +*/
{
  Anum                domnmax;                    /* Maximum number of domains */

  domnmax = (archVar (archptr) != 0)              /* If target architecture is variable-sized    */
            ? MIN (1023, grafptr->vertnbr)        /* Pre-set number of domains                   */
            : archDomSize (archptr, domnptr);     /* Else get architecture size                  */
  domnmax ++;                                     /* +1 for empty domain in mapBuild()/mapLoad() */

  mapInit2 (mappptr, grafptr, archptr, domnptr, domnmax, 0);
}

void
mapInit2 (
Mapping * restrict const        mappptr,          /*+ Mapping structure                  +*/
const Graph * restrict const    grafptr,          /*+ Graph data                         +*/
const Arch * restrict const     archptr,          /*+ Architecture data                  +*/
const ArchDom * restrict const  domnptr,          /*+ Target architecture initial domain +*/
const Anum                      domnmax,
const Anum                      domnnbr)
{
  mappptr->flagval = MAPPINGNONE;
  mappptr->grafptr = grafptr;
  mappptr->archptr = archptr;
  mappptr->parttax = NULL;
  mappptr->domntab = NULL;
  mappptr->domnnbr = domnnbr;
  mappptr->domnmax = domnmax;
  mappptr->domnorg = *domnptr;                    /* Use provided domain as original domain (e.g. when running a piece of parallel partitioning) */
}

/* This routine allocates the contents of a mapping.
** It returns:
** - 0   : if mapping successfully allocated.
** - !0  : on error.
*/

int
mapAlloc (
Mapping * restrict const    mappptr)              /*+ Mapping structure to fill +*/
{
  if ((mappptr->flagval & MAPPINGFREEPART) == 0) { /* If no private partition array yet */
    Anum * restrict     parttab;
    
    if ((parttab = (Anum *) memAlloc (mappptr->grafptr->vertnbr * sizeof (Anum))) == NULL) {
      errorPrint ("mapAlloc: out of memory (1)");
      return     (1);
    }
    mappptr->flagval |= MAPPINGFREEPART;
    mappptr->parttax  = parttab - mappptr->grafptr->baseval;
  }

  if ((mappptr->flagval & MAPPINGFREEDOMN) == 0) { /* If no private domain array yet */
    if ((mappptr->domntab = (ArchDom *) memAlloc (mappptr->domnmax * sizeof (ArchDom))) == NULL) {
      errorPrint ("mapAlloc: out of memory (2)");
      return     (1);
    }
    mappptr->flagval |= MAPPINGFREEDOMN;
  }

  return (0);
}

/* This routine resizes the domain array of a
** mapping and preserves its existing contents.
** It returns:
** - 0   : if mapping successfully allocated.
** - !0  : on error.
*/

int
mapResize (
Mapping * restrict const    mappptr,              /*+ Mapping structure to fill +*/
const Anum                  domnmax)
{
  int                       flagval;
  const ArchDom * restrict  domntab;

  flagval = mappptr->flagval;                     /* Save old flag value              */
  domntab = mappptr->domntab;                     /* Save pointer to old domain array */

  if (mapResize2 (mappptr, domnmax) != 0)         /* Resize array */
    return (1);

  if (flagval != mappptr->flagval)                /* If a new private array has been created */
    memCpy (mappptr->domntab, domntab, mappptr->domnnbr * sizeof (ArchDom));

  return (0);
}

/* This routine resizes the domain array of a
** mapping without preserving its existing contents.
** It returns:
** - 0   : if mapping successfully allocated.
** - !0  : on error.
*/

int
mapResize2 (
Mapping * restrict const    mappptr,              /*+ Mapping structure to fill +*/
const Anum                  domnmax)
{
  ArchDom *           domntab;

  domntab = ((mappptr->flagval & MAPPINGFREEDOMN) != 0) /* If it was a privately owned array */
            ? memRealloc (mappptr->domntab, domnmax * sizeof (ArchDom)) /* Reallocate it     */
            : memAlloc (domnmax * sizeof (ArchDom)); /* Else allocate it privately           */
  if (domntab == NULL) {
    errorPrint ("mapResize2: out of memory");
    return     (1);
  }

  mappptr->domntab  = domntab;
  mappptr->domnmax  = domnmax;
  mappptr->flagval |= MAPPINGFREEDOMN;            /* Array is now private anyway */

  return (0);
}

/* This routine builds an initial mapping.
** It returns:
** - void  : in all cases.
*/

void
mapFrst (
Mapping * restrict const    mappptr)              /*+ Mapping structure to fill +*/
{
  mappptr->domnnbr = 1;                           /* One domain in mapping to date */
  mappptr->domntab[0] = mappptr->domnorg;         /* Set first domain              */
  memSet (mappptr->parttax + mappptr->grafptr->baseval, 0, mappptr->grafptr->vertnbr * sizeof (Anum)); /* Set parttax to first domain */
}

/* This routine frees the contents
** of the given mapping.
** It returns:
** - void  : in all cases.
*/

void
mapFree (
Mapping * const             mappptr)
{
  if (((mappptr->flagval & MAPPINGFREEDOMN) != 0) && /* If domntab must be freed */
      (mappptr->domntab != NULL))                 /* And if exists               */
    memFree (mappptr->domntab);                   /* Free it                     */
  if (((mappptr->flagval & MAPPINGFREEPART) != 0) && /* If parttax must be freed */
      (mappptr->parttax != NULL))                 /* And if exists               */
    memFree (mappptr->parttax + mappptr->grafptr->baseval); /* Free it           */

  mappptr->parttax = NULL;
  mappptr->domntab = NULL;
}

/* This routine frees the contents
** of the given mapping.
** It returns:
** - void  : in all cases.
*/

void
mapExit (
Mapping * const             mappptr)
{
  mapFree (mappptr);

#ifdef SCOTCH_DEBUG_MAP2
  memSet (mappptr, ~0, sizeof (Mapping));
#endif /* SCOTCH_DEBUG_MAP2 */
}

/* This routine copies a mapping onto another.
** It returns:
** - 0   : if mapping successfully copied.
** - !0  : on error.
*/

int
mapCopy (
Mapping * restrict const       mappptr,           /*+ Mapping to set +*/
const Mapping * restrict const mapoptr)           /*+ Old mapping    +*/
{
  ArchDom *           domntab;
  Anum                domnnbr;
  Gnum                baseval;

#ifdef SCOTCH_DEBUG_MAP2
  if (mappptr->grafptr->vertnbr != mapoptr->grafptr->vertnbr) {
    errorPrint ("mapCopy: mappings do not match");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MAP2 */

  baseval = mapoptr->grafptr->baseval;
  domnnbr = mapoptr->domnnbr;
  if (domnnbr > mappptr->domnmax) {               /* If we have to resize domain array */
    if (mapResize2 (mappptr, domnnbr) != 0)       /* Resize it                         */
      return (1);
  }

  mappptr->domnnbr = domnnbr;
  memCpy (mappptr->domntab, mapoptr->domntab, domnnbr * sizeof (ArchDom));
  memCpy (mappptr->parttax + baseval, mapoptr->parttax + baseval, mapoptr->grafptr->vertnbr * sizeof (Anum));

  return (0);
}

/* This routine builds a mapping from a
** terminal domain partition array.
** It returns:
** - 0   : if mapping successfully filled.
** - !0  : on error.
*/

static
int
mapBuild2 (
Mapping * restrict const        mappptr,          /*+ Mapping to set                  +*/
MappingHash * restrict * const  hashtabptr,       /*+ Pointer to hash table to set up +*/
Gnum * const                    hashsizptr)       /*+ Size of hash table              +*/
{
  ArchDom             domndat;
  MappingHash *       hashtab;
  Gnum                hashnbr;                    /* Prospective number of cells in table */
  Gnum                hashsiz;                    /* Size of hash table                   */

  const Arch * restrict const archptr = mappptr->archptr;

#ifdef SCOTCH_DEBUG_MAP2
  if (mappptr->domnmax < 1) {
    errorPrint ("mapBuild2: domain array is too small");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MAP2 */

  archDomFrst (archptr, &domndat);
  hashnbr = (archVar (archptr) == 0)              /* If fixed size architecture                     */
            ? archDomSize (archptr, &domndat)     /* Get maximum size of distinct terminal domains  */
            : mappptr->grafptr->vertnbr;          /* Else take upper bound as number of vertices    */
  hashnbr ++;                                     /* Add one extra slot for unknown terminal domain */

  for (hashsiz = 32; hashsiz < hashnbr; hashsiz <<= 1) ; /* Get upper power of two  */
  hashsiz <<= 2;                                  /* Fill hash table at 25% maximum */

  if ((hashtab = (MappingHash *) memAlloc (hashsiz * sizeof (MappingHash))) == NULL) {
    errorPrint ("mapBuild2: out of memory");
    return (1);
  }
  memSet (hashtab, ~0, hashsiz * sizeof (MappingHash)); /* Set all vertex numbers to ~0 */

  *hashtabptr = hashtab;
  *hashsizptr = hashsiz;

  return (0);
}

static
int
mapBuild3 (
Mapping * restrict const      mappptr,            /*+ Mapping to fill        +*/
MappingHash * restrict const  hashtab,            /*+ Hash table             +*/
const Gnum                    hashsiz,            /*+ Hash table size        +*/
const Anum * restrict const   termtax)            /*+ Terminal array to load +*/
{
  ArchDom * restrict  domntab;
  Anum                domnnbr;
  Anum                domnmax;
  Gnum                hashmsk;
  Gnum                vertnnd;
  Gnum                vertnum;
  int                 o;

  const Arch * restrict const archptr = mappptr->archptr;
  Anum * restrict const       parttax = mappptr->parttax;

  o = 1;                                          /* Assume loop will fail */
  hashmsk = hashsiz - 1;
  domntab = mappptr->domntab;
  domnnbr = mappptr->domnnbr;
  domnmax = mappptr->domnmax;
  for (vertnum = mappptr->grafptr->baseval, vertnnd = mappptr->grafptr->vertnnd;
       vertnum < vertnnd; vertnum ++) {
    Gnum                hashnum;
    Anum                termnum;
    Anum                domnnum;

    termnum = termtax[vertnum];
    if (termnum == ~0)                            /* If unknown part, skip it */
      continue;

    for (hashnum = (termnum * MAPPINGHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
      if (hashtab[hashnum].termnum == termnum) {  /* If hash slot found  */
        domnnum = hashtab[hashnum].domnnum;       /* Domain number found */
        break;
      }
      if (hashtab[hashnum].termnum == ~0) {       /* If hash slot empty */
        hashtab[hashnum].termnum = termnum;       /* Create slot        */
        hashtab[hashnum].domnnum = domnnbr;

        if (domnnbr == domnmax) {
          domnmax += (domnmax >> 2) + 8;          /* Increase size by 25% */
          if (mapResize (mappptr, domnmax) != 0)
            goto fail;
          domntab = mappptr->domntab;             /* Re-read pointer in case it changed */
        }
        archDomTerm (archptr, &domntab[domnnbr], termnum); /* Create slot with terminal number domain */
        domnnum = domnnbr ++;                     /* Get position of new slot; one more slot created  */
        break;
      }
    }
    parttax[vertnum] = domnnum;                   /* Refer to the proper domain */
  }
  o = 0;                                          /* Success */

fail:
  mappptr->domnnbr = domnnbr;                     /* Set updated number of domains */

  memFree (hashtab);                              /* Free hash table */

  return (o);
}

int
mapBuild (
Mapping * restrict const    mappptr,              /*+ Mapping to set         +*/
const Anum * restrict const termtax)              /*+ Terminal array to load +*/
{
  MappingHash * restrict  hashtab;
  Gnum                    hashsiz;                /* Size of hash table */

  if (mapBuild2 (mappptr, &hashtab, &hashsiz) != 0)
    return (1);
  return (mapBuild3 (mappptr, hashtab, hashsiz, termtax));
}

/* This routine updates a mapping with a
** terminal domain partition array.
** It returns:
** - 0   : if mapping successfully updated.
** - !0  : on error.
*/

int
mapMerge (
Mapping * restrict const    mappptr,              /*+ Mapping to set         +*/
const Anum * restrict const termtab)              /*+ Terminal array to load +*/
{
  MappingHash * restrict  hashtab;
  Gnum                    hashsiz;                /* Size of hash table */
  Gnum                    hashmsk;
  Anum                    domnnbr;
  Anum                    domnnum;

  const Arch * restrict const     archptr = mappptr->archptr;
  const ArchDom * restrict const  domntab = mappptr->domntab;

  if (mapBuild2 (mappptr, &hashtab, &hashsiz) != 0)
    return (1);

  hashmsk = hashsiz - 1;
  for (domnnum = 0, domnnbr = mappptr->domnnbr; domnnum < domnnbr; domnnum ++) {
    const ArchDom *     domnptr;
    Gnum                hashnum;
    Anum                termnum;

    domnptr = &domntab[domnnum];
    if (archDomSize (archptr, domnptr) != 1)      /* If domain is not terminal, skip it */
      continue;

    termnum = archDomNum (archptr, domnptr);

    for (hashnum = (termnum * MAPPINGHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) { /* Fill hash table with existing domains */
#ifdef SCOTCH_DEBUG_MAP2
      if (hashtab[hashnum].termnum == termnum) {  /* If hash slot found                         */
        errorPrint ("mapMerge: internal error");  /* Multiple domains with same terminal number */
        return     (1);
      }
#endif /* SCOTCH_DEBUG_MAP2 */
      if (hashtab[hashnum].termnum == ~0) {       /* If hash slot empty */
        hashtab[hashnum].termnum = termnum;       /* Create slot        */
        hashtab[hashnum].domnnum = domnnum;
        break;
      }
    }
  }

  return (mapBuild3 (mappptr, hashtab, hashsiz, termtab)); /* Add new domains to existing domain array */
}

/* This routine propagates back mapping
** information to a terminal part array.
** It returns:
** - void  : in all cases.
*/

void
mapTerm (
const Mapping * restrict const  mappptr,
Anum * restrict const           termtax)
{
  Gnum                vertnnd;
  Gnum                vertnum;

  const Arch * restrict const     archptr = mappptr->archptr;
  const ArchDom * restrict const  domntab = mappptr->domntab;
  const Anum * restrict const     parttax = mappptr->parttax;

  vertnum = mappptr->grafptr->baseval;
  if (domntab != NULL) {
    for (vertnnd = mappptr->grafptr->vertnnd;
         vertnum < vertnnd; vertnum ++)
      termtax[vertnum] = archDomNum (archptr, &domntab[parttax[vertnum]]);
  }
  else
    memSet (termtax + vertnum, ~0, mappptr->grafptr->vertnbr * sizeof (Anum));
}
