/* Copyright 2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : common_stub.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are common stub routines    **/
/**                for several systems, used by all        **/
/**                modules.                                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 11 may 2008     **/
/**                                 to     11 may 2008     **/
/**                # Version 5.1  : from : 27 jun 2010     **/
/**                                 to     01 jul 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"

/*
**  The static and global variables.
*/

int commonStubDummy;                              /* Dummy variable so that the object file will never be empty */

/******************/
/*                */
/* Stub routines. */
/*                */
/******************/

#ifdef COMMON_STUB_FORK
int
fork ()
{
  errorPrint ("fork() not implemented on this system");
  return     (-1);
}
#endif /* COMMON_STUB_FORK */
