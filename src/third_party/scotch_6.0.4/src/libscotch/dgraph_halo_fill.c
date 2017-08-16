/* Copyright 2007-2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : dgraph_halo_fill.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel static mapper.       **/
/**                This module contains the halo update    **/
/**                routines.                               **/
/**                                                        **/
/**                # Version 5.0  : from : 31 dec 2006     **/
/**                                 to     05 feb 2008     **/
/**                # Version 5.1  : from : 28 aug 2008     **/
/**                                 to     29 aug 2010     **/
/**                                                        **/
/************************************************************/

/* This function fills the send array used by
** all of the halo routines.
** It returns:
** - void  : in all cases.
*/

static
void
DGRAPHHALOFILLNAME (
const Dgraph * restrict const grafptr,
const void * restrict const   attrgsttab,         /* Attribute array to diffuse           */
const int                     attrglbsiz,         /* Type extent of attribute             */
byte ** restrict const        attrdsptab)         /* Temporary address displacement array */
{
  byte * restrict       attrgstptr;
  const int * restrict  procsidptr;
  const int * restrict  procsidnnd;

  for (procsidptr = grafptr->procsidtab, procsidnnd = procsidptr + grafptr->procsidnbr, attrgstptr = (byte *) attrgsttab;
       procsidptr < procsidnnd; procsidptr ++) {
    int                 procsidval;

    procsidval = *procsidptr;
    if (procsidval < 0)
      attrgstptr -= ((Gnum) procsidval) * DGRAPHHALOFILLSIZE;
    else {
      byte *              attrdspptr;

      attrdspptr = attrdsptab[procsidval];
      attrdsptab[procsidval] = attrdspptr + DGRAPHHALOFILLSIZE; /* Skip to next position in send buffer */
      DGRAPHHALOFILLCOPY (attrdspptr, attrgstptr, DGRAPHHALOFILLSIZE);
    }
  }
}
