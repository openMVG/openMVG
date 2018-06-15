/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dorder_io_block.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles distributed         **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 26 may 2008     **/
/**                                 to     26 may 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define DORDER

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "order.h"
#include "dorder.h"

/************************************/
/*                                  */
/* These routines handle orderings. */
/*                                  */
/************************************/

/* This routine saves a distributed ordering on
** a combined block ordering format.
** The distributed graph structure is provided
** to access the distribution of vertex labels,
** whenever present.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dorderSaveBlock2 (
const Order * const           cordptr,
const Gnum * const            vlbltab,
FILE * const                  stream)
{
  Gnum              vertnum;
  Gnum * restrict   rangtab;
  Gnum              cblknum;
  int               o;

  if ((rangtab = memAlloc ((cordptr->vnodnbr + 1) * sizeof (Gnum))) == NULL) {
    errorPrint ("dorderSaveBlock2: out of memory");
    return     (1);
  }
  orderRang (cordptr, rangtab);

  if (fprintf (stream, "0\n" GNUMSTRING "\t" GNUMSTRING "\n",
               (Gnum) cordptr->cblknbr,
               (Gnum) cordptr->vnodnbr) < 0) {
    errorPrint ("dorderSaveBlock2: bad output (1)");
    return     (1);
  }

  for (cblknum = 0, o = 1; (o == 1) && (cblknum < cordptr->cblknbr); cblknum ++) { /* Save column-block range array */
    o = intSave (stream, rangtab[cblknum]);
    putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, rangtab[cblknum]);
  putc ('\n', stream);

  orderPeri (cordptr->peritab, cordptr->baseval, cordptr->vnodnbr, rangtab, cordptr->baseval); /* TRICK: re-use rangtab as permtab */

  for (vertnum = 0; (o == 1) && (vertnum < (cordptr->vnodnbr - 1)); vertnum ++) { /* Save direct permutation */
    o = intSave (stream, rangtab[vertnum]);
    putc (((vertnum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, rangtab[vertnum]);
  putc ('\n', stream);

  if (o != 1)
    errorPrint ("dorderSaveBlock2: bad output (2)");

  return (1 - o);
}

int
dorderSaveBlock (
const Dorder * restrict const ordeptr,
const Dgraph * restrict const grafptr,
FILE * restrict const         stream)
{
  return (dorderSaveTree2 (ordeptr, grafptr, stream, dorderSaveBlock2));
}
