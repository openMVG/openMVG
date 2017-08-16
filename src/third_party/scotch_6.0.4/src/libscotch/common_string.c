/* Copyright 2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : common_string.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are common routines used    **/
/**                by all modules.                         **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 23 jul 2010     **/
/**                                 to     23 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include <time.h>
#include "common.h"

/********************************/
/*                              */
/* String substitution routine. */
/*                              */
/********************************/

static
void
stringSubst2 (
char * const                bsrcptr,
char * const                bdstptr,
const char * const          pattstr,
const char * const          replstr,
const int                   pattsiz,
const int                   replsiz)
{
  char *              pattptr;
  int                 pattidx;

  pattptr = strstr (bsrcptr, pattstr);            /* Search for the pattern in the remaining source string   */
  pattidx = (pattptr == NULL) ? (strlen (bsrcptr) + 1): (pattptr - bsrcptr); /* Get length of unchanged part */

  if (replsiz < pattsiz)                          /* If replacement is smaller, pre-move unchanged part */
    memMov (bdstptr, bsrcptr, pattidx * sizeof (char));

  if (pattptr != NULL)                            /* If remaining part of string has to be processed */
    stringSubst2 (pattptr + pattsiz, bdstptr + pattidx + replsiz, pattstr, replstr, pattsiz, replsiz);

  if (replsiz > pattsiz)                          /* If replacement is longer, post-move unchanged part */
    memMov (bdstptr, bsrcptr, pattidx * sizeof (char));

  if (pattptr != NULL)                            /* If there is something to replace         */
    memCpy (bdstptr + pattidx, replstr, replsiz * sizeof (char)); /* Write replacement string */

  return;
}

void
stringSubst (
char * const                buffptr,              /* String to search into */
const char * const          pattstr,              /* Pattern to search for */
const char * const          replstr)              /* Replacement string    */
{
  int                 pattsiz;
  int                 replsiz;

  pattsiz = strlen (pattstr);
  replsiz = strlen (replstr);

  stringSubst2 (buffptr, buffptr, pattstr, replstr, pattsiz, replsiz);
}
