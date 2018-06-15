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
/**   NAME       : common_error.c                          **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module handles errors.             **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 08 may 1998     **/
/**                                 to     02 oct 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON_ERROR

#include "common.h"

/********************************/
/*                              */
/* The error handling routines. */
/*                              */
/********************************/

static char                 errorProgName[32] = "";

/* This routine sets the program name for
** error reporting.
** It returns:
** - VOID  : in all cases.
*/

void
errorProg (
const char * const          progstr)              /*+ Program name +*/
{
  strncpy (errorProgName, progstr, 29);
  errorProgName[29] = '\0';
  strcat  (errorProgName, ": ");

}

/* This routine prints an error message with
** a variable number of arguments, as printf ()
** does, and exits.
** It returns:
** - EXIT  : in all cases.
*/

void
errorPrint (
const char * const          errstr,               /*+ printf-like variable argument list */
...)
{
  va_list             errlist;                    /* Argument list of the call */

  fprintf  (stderr, "%sERROR: ", errorProgName);
  va_start (errlist, errstr);
  vfprintf (stderr, errstr, errlist);             /* Print arguments */
  va_end   (errlist);
  fprintf  (stderr, "\n");
  fflush   (stderr);                              /* In case it has been set to buffered mode */
}

/* This routine prints a warning message with
** a variable number of arguments, as printf ()
** does.
** It returns:
** - VOID  : in all cases.
*/

void
errorPrintW (
const char * const          errstr,               /*+ printf-like variable argument list */
...)
{
  va_list             errlist;                    /* Argument list of the call */

  fprintf  (stderr, "%sWARNING: ", errorProgName);
  va_start (errlist, errstr);
  vfprintf (stderr, errstr, errlist);             /* Print arguments */
  va_end   (errlist);
  fprintf  (stderr, "\n");
  fflush   (stderr);                              /* In case it has been set to buffered mode */
}
