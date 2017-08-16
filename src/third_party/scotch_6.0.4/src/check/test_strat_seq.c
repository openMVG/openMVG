/* Copyright 2012,2013 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_strat_seq.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the sequential        **/
/**                strategy building routines.             **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 08 jan 2012     **/
/**                                 to     11 oct 2013     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>

#include "scotch.h"

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                 argc,
char *              argv[])
{
  SCOTCH_Strat        stradat;

  SCOTCH_errorProg (argv[0]);

  printf ("Sequential mapping strategy, SCOTCH_STRATDEFAULT\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphMapBuild (&stradat, SCOTCH_STRATDEFAULT, 16, 0.03);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential mapping strategy, SCOTCH_STRATRECURSIVE\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphMapBuild (&stradat, SCOTCH_STRATRECURSIVE, 16, 0.03);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential mapping strategy, SCOTCH_STRATREMAP\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphMapBuild (&stradat, SCOTCH_STRATREMAP, 16, 0.03);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential mapping strategy, SCOTCH_STRATRECURSIVE | SCOTCH_STRATREMAP\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphMapBuild (&stradat, SCOTCH_STRATRECURSIVE | SCOTCH_STRATREMAP, 16, 0.03);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATDEFAULT\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATDEFAULT, 0, 0.2);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATLEVELMAX\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATLEVELMAX, 3, 0.2);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATLEVELMIN\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATLEVELMIN, 3, 0.2);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN, 3, 0.2);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATLEAFSIMPLE\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATLEAFSIMPLE, 3, 0.2);
  SCOTCH_stratExit (&stradat);

  printf ("Sequential ordering strategy, SCOTCH_STRATSEPASIMPLE\n");

  SCOTCH_stratInit (&stradat);
  SCOTCH_stratGraphOrderBuild (&stradat, SCOTCH_STRATSEPASIMPLE, 3, 0.2);
  SCOTCH_stratExit (&stradat);

  return (0);
}
