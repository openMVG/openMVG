/* Copyright 2004,2007 ENSEIRB, INRIA & CNRS
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
/**   NAME       : symbol.c                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the general purpose     **/
/**                routines for the symbolic matrix.       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 0.1  : from : 03 dec 1998     **/
/**                                 to     03 dec 1998     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL

#include "common.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine initializes the given
*** symbolic block matrix structure.
*** It returns:
*** - 0  : in all cases.
+*/

int
symbolInit (
SymbolMatrix * const        symbptr)
{
  memSet (symbptr, 0, sizeof (SymbolMatrix));

  return (0);
}

/*+ This routine frees the contents
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolExit (
SymbolMatrix * const        symbptr)
{
  if (symbptr->cblktab != NULL)
    memFree (symbptr->cblktab);
  if (symbptr->bloktab != NULL)
    memFree (symbptr->bloktab);

#ifdef SYMBOL_DEBUG
  symbolInit (symbptr);
#endif /* SYMBOL_DEBUG */
}

/*+ This routine reallocates the arrays
*** of the given symbolic block matrix.
*** It returns:
*** - VOID  : in all cases.
+*/

void
symbolRealloc (
SymbolMatrix * const        symbptr)
{
  SymbolCblk *        cblktab;
  SymbolBlok *        bloktab;

  if ((cblktab = (SymbolCblk *) memAlloc ((symbptr->cblknbr + 1) * sizeof (SymbolCblk))) == NULL)
    return;                                       /* Cannot move smallest array */
  memCpy  (cblktab, symbptr->cblktab, (symbptr->cblknbr + 1) * sizeof (SymbolCblk));
  memFree (symbptr->cblktab);                     /* Move column block array */
  symbptr->cblktab = cblktab;

  if ((bloktab = (SymbolBlok *) memAlloc (symbptr->bloknbr * sizeof (SymbolBlok))) == NULL)
    return;                                       /* Cannot move array */
  memCpy  (bloktab, symbptr->bloktab, symbptr->bloknbr * sizeof (SymbolBlok));
  memFree (symbptr->bloktab);                     /* Move column block array */
  symbptr->bloktab = bloktab;
}
