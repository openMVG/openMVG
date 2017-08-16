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
/**   NAME       : arch_deco.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the decomposition-  **/
/**                defined target architecture.            **/
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
/**                                 to     14 sep 1995     **/
/**                # Version 3.1  : from : 20 jul 1996     **/
/**                                 to     23 jul 1996     **/
/**                # Version 3.2  : from : 11 sep 1996     **/
/**                                 to     28 sep 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     17 may 1999     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.0  : from : 10 sep 2007     **/
/**                                 to     28 feb 2008     **/
/**                # Version 5.1  : from : 21 jan 2008     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     14 fev 2011     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_DECO

#include "module.h"
#include "common.h"
#include "arch.h"
#include "arch_deco.h"

/***************************************/
/*                                     */
/* These are the decomposition-defined */
/* architecture routines.              */
/*                                     */
/***************************************/

/* This routine builds a compiled
** decomposition-defined architecture
** from the raw terminal tables that are
** passed to it.
** It returns:
** - 0   : if the decomposition has been successfully computed.
** - !0  : on error.
*/

int
archDecoArchBuild (
ArchDeco * restrict const       archptr,          /*+ Architecture to build                            +*/
const Anum                      termdomnbr,       /*+ Number of terminal domains (ie processors)       +*/
const Anum                      termdommax,       /*+ Maximum domain number given to a terminal domain +*/
const ArchDecoTermVert * const  termverttab,      /*+ Terminal vertex array                            +*/
const Anum * const              termdisttab)      /*+ Terminal distance map                            +*/
{
  Anum                i, j, k;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchDeco)    > sizeof (ArchDummy)) ||
      (sizeof (ArchDecoDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archDecoArchBuild: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (memAllocGroup ((void **) (void *)
                     &archptr->domverttab, (size_t) (termdommax * sizeof (ArchDecoVert)),
                     &archptr->domdisttab, (size_t) ((((termdommax * (termdommax - 1)) / 2) + 1) * sizeof (Anum)), NULL) == NULL) {
    errorPrint ("archDecoArchBuild: out of memory");
    return     (1);
  }
  archptr->flagval    = ARCHDECOFREE;
  archptr->domtermnbr = termdomnbr;
  archptr->domvertnbr = termdommax;

  for (i = 0; i < termdommax; i ++) {
    archptr->domverttab[i].labl = ARCHDOMNOTTERM; /* Assume domain is not a terminal */
    archptr->domverttab[i].size = 0;              /* Assume domain is not used (yet)  */
    archptr->domverttab[i].wght = 0;              /* Assume domain is not used (yet)  */
  }

  for (i = 0; i < termdomnbr; i ++) {             /* Set terminal data of all declared terminals */
#ifdef SCOTCH_DEBUG_ARCH1
    if (termverttab[i].num > termdommax) {        /* If incorrect maximum terminal number */
      errorPrint       ("archDecoArchBuild: bad maximum terminal");
      archDecoArchFree (archptr);
      return           (1);
    }
#endif /* SCOTCH_DEBUG_ARCH1 */

    archptr->domverttab[termverttab[i].num - 1].labl = termverttab[i].labl;
    archptr->domverttab[termverttab[i].num - 1].size = 1;
    archptr->domverttab[termverttab[i].num - 1].wght = termverttab[i].wght;
  }

  for (i = termdommax - 1; i > 0; i --) {         /* Accumulate data from terminals to root */
    j = ((i - 1) >> 1);
    if (archptr->domverttab[i].labl != ARCHDOMNOTTERM) {
      if ((archptr->domverttab[j].labl == ARCHDOMNOTTERM) || /* Get smallest label */
          (archptr->domverttab[j].labl > archptr->domverttab[i].labl))
        archptr->domverttab[j].labl = archptr->domverttab[i].labl;
      archptr->domverttab[j].size += archptr->domverttab[i].size;
      archptr->domverttab[j].wght += archptr->domverttab[i].wght;
    }
  }
#ifdef SCOTCH_DEBUG_ARCH1
  if (archptr->domverttab[0].size != termdomnbr) { /* If incorrect accumulation */
    errorPrint ("archDecoArchBuild: bad terminal count");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  memSet (archptr->domdisttab, 0, termdommax * (termdommax - 1) * sizeof (Anum) / 2);
                                                  /* Assume distance is not known   */
  for (i = 1, k = 0; i < termdomnbr; i ++) {      /* Read the terminal distance map */
    for (j = 0; j < i; j ++, k ++)
      archDecoArchDist (archptr, termverttab[i].num, termverttab[j].num) = termdisttab[k];
  }

  for (j = termdommax; j > 0; j --) {             /* Loop on domains              */
    if (archDecoArchSize (archptr, j) == 0)       /* If domain is unused, skip it */
      continue;
    for (i = termdommax; i > j; i --) {           /* Double loop on distance array values */
      if (archDecoArchSize (archptr, i) == 0)     /* If domain is unused, skip it         */
        continue;
      if (archDecoArchSize (archptr, i) > 1) {    /* If domain i has subdomains */
        if (archDecoArchSize (archptr, j) > 1)    /* If domain j has subdomains */
          archDecoArchDist (archptr, i, j) = (archDecoArchDistE (archptr, 2 * i,     2 * j)     +
                                              archDecoArchDistE (archptr, 2 * i,     2 * j + 1) +
                                              archDecoArchDistE (archptr, 2 * i + 1, 2 * j)     +
                                              archDecoArchDistE (archptr, 2 * i + 1, 2 * j + 1) + 2) / 4;
        else                                      /* If domain j is a terminal */
          archDecoArchDist (archptr, i, j) = (archDecoArchDistE (archptr, 2 * i,     j) +
                                              archDecoArchDistE (archptr, 2 * i + 1, j) + 1) / 2;
      }
      else {                                      /* If domain i is a terminal  */
        if (archDecoArchSize (archptr, j) > 1)    /* If domain j has subdomains */
          archDecoArchDist (archptr, i, j) = (archDecoArchDistE (archptr, i, 2 * j)     +
                                              archDecoArchDistE (archptr, i, 2 * j + 1) + 1) / 2;
#ifdef SCOTCH_DEBUG_ARCH1
        else {                                    /* If both domain are terminals                  */
          if (archDecoArchDist (archptr, i, j) == 0) { /* Distance value must be greater than zero */
            errorPrint       ("archDecoArchBuild: invalid null distance");
            archDecoArchFree (archptr);
            return           (1);
          }
        }
#endif /* SCOTCH_DEBUG_ARCH1 */
      }
    }
  }

  return (0);
}

/* This routine loads and computes the
** decomposition-defined architecture
** tables.
** It returns:
** - 0   : if the decomposition has been successfully read.
** - !0  : on error.
*/

int
archDecoArchLoad (
ArchDeco * restrict const   archptr,
FILE * restrict const       stream)
{
  INT                         decotype;           /* Type of decomposition                            */
  INT                         termdomnbr;         /* Number of terminal domains (ie processors)       */
  INT                         termdommax;         /* Maximum domain number given to a terminal domain */
  ArchDecoTermVert * restrict termverttab;        /* Table of terminal vertex data                    */
  Anum * restrict             termdisttab;        /* Table of terminal-to-terminal distances          */
  INT                         i, j;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchDeco)    > sizeof (ArchDummy)) ||
      (sizeof (ArchDecoDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archDecoArchLoad: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((intLoad (stream, &decotype)   != 1) ||     /* Read header */
      (intLoad (stream, &termdomnbr) != 1) ||
      (intLoad (stream, &termdommax) != 1) ||
      (decotype   < 0)                     ||
      (decotype   > 1)                     ||
      (termdommax < termdomnbr)            ||
      (termdomnbr < 1)) {
    errorPrint ("archDecoArchLoad: bad input (1)");
    return     (1);
  }

  if (decotype == 0) {                            /* If raw decomposition */
    if (memAllocGroup ((void **) (void *)
                       &termverttab, (size_t) (termdomnbr * sizeof (ArchDecoTermVert)),
                       &termdisttab, (size_t) ((((termdommax * (termdommax - 1)) / 2) + 1) * sizeof (Anum)), NULL) == NULL) {
      errorPrint ("archDecoArchLoad: out of memory (1)");
      return     (1);
    }

    for (i = 0; i < termdomnbr; i ++) {           /* For all declared terminals  */
      INT                 termvertlabl;
      INT                 termvertwght;
      INT                 termvertnum;

      if ((intLoad (stream, &termvertlabl) != 1) || /* Read terminal data */
          (intLoad (stream, &termvertwght) != 1) ||
          (intLoad (stream, &termvertnum)  != 1) ||
          (termvertnum < 1)                      ||
          (termvertnum > termdommax)) {
        errorPrint       ("archDecoArchLoad: bad input (2)");
        memFree          (termverttab);           /* Free group leader */
        return           (1);
      }
      termverttab[i].labl = (ArchDomNum) termvertlabl;
      termverttab[i].wght = (Anum)       termvertwght;
      termverttab[i].num  = (Anum)       termvertnum;
    }

    for (i = 0, j = (termdomnbr * (termdomnbr - 1)) / 2; i < j; i ++) { /* Read terminal distance map */
      INT                 termdistval;

      if ((intLoad (stream, &termdistval) != 1) ||
          (termdistval < 1)) {
        errorPrint       ("archDecoArchLoad: bad input (3)");
        memFree          (termverttab);           /* Free group leader */
        return           (1);
      }
      termdisttab[i] = (Anum) termdistval;
    }

    archDecoArchBuild (archptr, termdomnbr, termdommax, termverttab, termdisttab);

    memFree (termverttab);                        /* Free group leader */
  }
  else {                                          /* If it is a compiled decomposition */
    if (memAllocGroup ((void **) (void *)
                       &archptr->domverttab, (size_t) (termdommax * sizeof (ArchDecoVert)),
                       &archptr->domdisttab, (size_t) ((((termdommax * (termdommax - 1)) / 2) + 1) * sizeof (Anum)), NULL) == NULL) {
      errorPrint       ("archDecoArchLoad: out of memory (2)");
      return           (1);
    }
    archptr->flagval    = ARCHDECOFREE;
    archptr->domtermnbr = (Anum) termdomnbr;
    archptr->domvertnbr = (Anum) termdommax;

    for (i = 0; i < termdommax; i ++) {           /* Read domain array */
      INT                 domvertlabl;
      INT                 domvertsize;
      INT                 domvertwght;

      if ((intLoad (stream, &domvertlabl) != 1) ||
	  (intLoad (stream, &domvertsize) != 1) ||
	  (intLoad (stream, &domvertwght) != 1)) {
        errorPrint       ("archDecoArchLoad: bad input (4)");
        archDecoArchFree (archptr);
        return           (1);
      }
      archptr->domverttab[i].labl = (ArchDomNum) domvertlabl;
      archptr->domverttab[i].size = (Anum)       domvertsize;
      archptr->domverttab[i].wght = (Anum)       domvertwght;
    }

    for (i = 0; i < (termdommax * (termdommax - 1)) / 2; i ++) { /* Read distance array */
      INT                 domdistval;

      if (intLoad (stream, &domdistval) != 1) {
        errorPrint       ("archDecoArchLoad: bad input (5)");
        archDecoArchFree (archptr);
        return           (1);
      }
      archptr->domdisttab[i] = domdistval;
    }
  }

  return (0);
}

/* This routine frees the decomposition
** architecture structures.
** It returns:
** - 0   : if the decomposition has been successfully freed.
** - !0  : on error.
*/

int
archDecoArchFree (
ArchDeco * const            archptr)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchDeco)    > sizeof (ArchDummy)) ||
      (sizeof (ArchDecoDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archDecoArchFree: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (((archptr->flagval & ARCHDECOFREE) != 0) &&
      (archptr->domverttab != NULL))
    memFree (archptr->domverttab);                /* Free group leader */

  archptr->domtermnbr =
  archptr->domvertnbr = 0;
  archptr->domverttab = NULL;
  archptr->domdisttab = NULL;

  return (0);
}

/* This routine saves the given target architecture
** as compiled decomposition tables.
** It returns:
** - 0   : if the decomposition has been successfully written.
** - !0  : on error.
*/

int
archDecoArchSave (
const ArchDeco * const      archptr,
FILE * restrict const       stream)
{
  Anum                i, j;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchDeco)    > sizeof (ArchDummy)) ||
      (sizeof (ArchDecoDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archDecoArchSave: invalid type specification");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, "1\n" ANUMSTRING "\t" ANUMSTRING "\n", /* Write number of domains */
               (Anum) archptr->domtermnbr,
               (Anum) archptr->domvertnbr) == EOF) {
    errorPrint ("archDecoArchSave: bad output (1)");
    return     (1);
  }

  for (i = 0; i < archptr->domvertnbr; i ++) {    /* Write domain array */
    if (fprintf (stream, ANUMSTRING "\t" ANUMSTRING "\t" ANUMSTRING "\n",
                 (Anum) archptr->domverttab[i].labl,
                 (Anum) archptr->domverttab[i].size,
                 (Anum) archptr->domverttab[i].wght) == EOF) {
      errorPrint ("archDecoArchSave: bad output (2)");
      return     (1);
    }
  }

  j = (archptr->domvertnbr * (archptr->domvertnbr - 1)) / 2;
  for (i = 0; i < j; i ++) {                      /* Write distance array */
    if (fprintf (stream, ANUMSTRING "%c",
                 (Anum) archptr->domdisttab[i],
                 (((i % 8) == 7) && (i != (j - 1))) ? '\n' : '\t') == EOF) {
      errorPrint ("archDecoArchSave: bad output (3)");
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
archDecoDomNum (
const ArchDeco * const      archptr,
const ArchDecoDom * const   domptr)
{
  return (archptr->domverttab[domptr->num - 1].labl);
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archDecoDomTerm (
const ArchDeco * const      archptr,
ArchDecoDom * const         domptr,
const ArchDomNum            domnum)
{
  Anum                domtermnum;
  Anum                domvertnum;

  for (domtermnum = archptr->domtermnbr, domvertnum = archptr->domvertnbr - 1;
       (domtermnum > 0) && (domvertnum != (Anum) (-1)); domvertnum --) {
    if (archptr->domverttab[domvertnum].size == 1) { /* If terminal vertex                     */
      domtermnum --;                              /* One more terminal scanned                 */
      if (archptr->domverttab[domvertnum].labl == domnum) { /* If terminal domain number found */
        domptr->num = domvertnum;                 /* Set domain number                         */
        return (0);
      }
    }
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the given domain.
*/

Anum
archDecoDomSize (
const ArchDeco * const      archptr,
const ArchDecoDom * const   domptr)
{
  return (archptr->domverttab[domptr->num - 1].size);
}

/* This function returns the weight of
** the given domain.
*/

Anum
archDecoDomWght (
const ArchDeco * const      archptr,
const ArchDecoDom * const   domptr)
{
  return (archptr->domverttab[domptr->num - 1].wght);
}

/* This function returns the average distance
** between two domains, which is extracted
** from the table.
*/

Anum
archDecoDomDist (
const ArchDeco * const      archptr,
const ArchDecoDom * const   dom0ptr,
const ArchDecoDom * const   dom1ptr)
{
  return (archDecoArchDistE (archptr, dom0ptr->num, dom1ptr->num));
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archDecoDomFrst (
const ArchDeco * const        archptr,
ArchDecoDom * restrict const  domptr)
{
  domptr->num = 1;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archDecoDomLoad (
const ArchDeco * const        archptr,
ArchDecoDom * restrict const  domptr,
FILE * restrict const         stream)
{
  if ((intLoad (stream, &domptr->num) != 1) ||
      (domptr->num < 1) || (domptr->num > archptr->domvertnbr)) {
    errorPrint ("archDecoDomLoad: bad input");
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
archDecoDomSave (
const ArchDeco * const      archptr,
const ArchDecoDom * const   domptr,
FILE * restrict const       stream)
{
  if (fprintf (stream, ANUMSTRING " ",
               (Anum) domptr->num) == EOF) {
    errorPrint ("archDecoDomSave: bad output");
    return     (1);
  }

  return (0);
}

/* This function tries to split a
** decomposition domain into two
** subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archDecoDomBipart (
const ArchDeco * const        archptr,
const ArchDecoDom * const     domptr,
ArchDecoDom * restrict const  dom0ptr,
ArchDecoDom * restrict const  dom1ptr)
{
  if (archptr->domverttab[domptr->num - 1].size <= 1) /* Return if cannot bipartition more */
    return (1);

  dom0ptr->num = domptr->num << 1;                /* Compute subdomain numbers from domain number */
  dom1ptr->num = dom0ptr->num + 1;

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
archDecoDomIncl (
const ArchDeco * const      archptr,
const ArchDecoDom * const   dom0ptr,
const ArchDecoDom * const   dom1ptr)
{
  Anum          dom0num;
  Anum          dom1num;

  for (dom1num = dom1ptr->num, dom0num = dom0ptr->num; dom1num != 0; dom1num >>= 1)
    if (dom1num == dom0ptr->num)
      return (1);

  return (0);
}

/* This function creates the MPI_Datatype for
** decomposition-described domains.
** It returns:
** - 0  : if type could be created.
** - 1  : on error.
*/

#ifdef SCOTCH_PTSCOTCH
int
archDecoDomMpiType (
const ArchDeco * const        archptr,
MPI_Datatype * const          typeptr)
{
  *typeptr = ANUM_MPI;

  return (0);
}
#endif /* SCOTCH_PTSCOTCH */
