/* Copyright 2004,2007,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_parser.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the generic  **/
/**                strategy handling routines of the       **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 aug 1998     **/
/**                                 to     19 aug 1998     **/
/**   DATES      : # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     31 may 1999     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 23 dec 2001     **/
/**                                 to   : 23 dec 2001     **/
/**                # Version 6.0  : from : 07 jan 2014     **/
/**                                 to     07 jan 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "parser.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the strategy handling routines.  */
/*                                  */
/************************************/

/* This routine initializes a strategy
** structure.
** It returns:
** - 0  : in all cases.
*/

int
SCOTCH_stratInit (
SCOTCH_Strat * const        stratptr)
{
  if (sizeof (SCOTCH_Strat) < sizeof (Strat *)) {
    errorPrint ("SCOTCH_stratInit: internal error (1)");
    return     (1);
  }

  *((Strat **) stratptr) = NULL;                  /* Initialize pointer to strategy */

  return (0);
}

/* This routine frees a strategy structure.
** It returns:
** - VOID  : in all cases.
*/

void
SCOTCH_stratExit (
SCOTCH_Strat * const        stratptr)
{
  if (*((Strat **) stratptr) != NULL)             /* If strategy is not null */
    stratExit (*((Strat **) stratptr));           /* Free strategy structure */
}

/* This routine frees the contents of a
** strategy structure and cleans it for
** future use.
** It returns:
** - VOID  : in all cases.
*/

void
SCOTCH_stratFree (
SCOTCH_Strat * const        stratptr)
{
  if (*((Strat **) stratptr) != NULL) {           /* If strategy is not null        */
    stratExit (*((Strat **) stratptr));           /* Free strategy structure        */
    *((Strat **) stratptr) = NULL;                /* Initialize pointer to strategy */
  }
}

/* This routine outputs the contents of the
** given strategy to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
SCOTCH_stratSave (
const SCOTCH_Strat * const  stratptr,
FILE * const                stream)
{
  return (stratSave (*((Strat **) stratptr), stream));
}
