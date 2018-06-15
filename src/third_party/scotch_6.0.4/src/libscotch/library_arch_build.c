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
/**   NAME       : library_arch_build.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the target   **/
/**                architecture building routine of the    **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.3  : from : 02 oct 1998     **/
/**                                 to     29 mar 1999     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 08 mar 2005     **/
/**                                 to     17 mar 2005     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_build.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "scotch.h"

/**************************************/
/*                                    */
/* These routines are the C API for   */
/* the architecture building routine. */
/*                                    */
/**************************************/

/*+ This routine parses the given
*** bipartitioning strategy.
*** It returns:
*** - 0   : if string successfully scanned.
*** - !0  : on error.
+*/

int
SCOTCH_stratGraphBipart (
SCOTCH_Strat * const        stratptr,
const char * const          string)
{
  if (*((Strat **) stratptr) != NULL)
    stratExit (*((Strat **) stratptr));

  if ((*((Strat **) stratptr) = stratInit (&bgraphbipartststratab, string)) == NULL) {
    errorPrint ("SCOTCH_stratBipart: error in bipartitioning strategy");
    return     (1);
  }

  return (0);
}

/*+ This routine fills the contents of the given
*** opaque target structure with the data provided
*** by the user. The source graph provided on input
*** is turned into a decomposition-defined target
*** architecture.
*** It returns:
*** - 0   : if the computation succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_archBuild (
SCOTCH_Arch * const         archptr,              /*+ Target architecture to build    +*/
const SCOTCH_Graph * const  grafptr,              /*+ Graph to turn into architecture +*/
const SCOTCH_Num            listnbr,              /*+ Number of elements in sublist   +*/
const SCOTCH_Num * const    listptr,              /*+ Pointer to sublist              +*/
const SCOTCH_Strat * const  stratptr)             /*+ Bipartitoning strategy          +*/
{
  Strat *             bipstratptr;
  VertList            graflistdat;
  VertList *          graflistptr;
  int                 o;

  if ((sizeof (SCOTCH_Num) != sizeof (Gnum)) ||
      (sizeof (SCOTCH_Num) != sizeof (Anum))) {
    errorPrint ("SCOTCH_archBuild: internal error (1)");
    return     (1);
  }

  if (*((Strat **) stratptr) == NULL)             /* Set default mapping strategy if necessary */
    *((Strat **) stratptr) = stratInit (&bgraphbipartststratab, "(m{vert=50,low=h{pass=10},asc=f{move=100,bal=0.1}}f{move=100,bal=0.05})(/((load0=load)|(load0=0))?x;)");
  bipstratptr = *((Strat **) stratptr);
  if (bipstratptr->tabl != &bgraphbipartststratab) {
    errorPrint ("SCOTCH_archBuild: not a bipartitioning strategy");
    return     (1);
  }

  if ((listnbr == ((Graph *) grafptr)->vertnbr) || (listnbr == 0) || (listptr == NULL))
    graflistptr = NULL;
  else {
    graflistptr = &graflistdat;
    graflistdat.vnumnbr = (Gnum)   listnbr;
    graflistdat.vnumtab = (Gnum *) listptr;
  }

  o = archBuild ((Arch * const) archptr, (const Graph * const) grafptr, graflistptr, bipstratptr);

  return (o);
}
