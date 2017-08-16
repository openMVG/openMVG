/* Copyright 2004,2007,2010,2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hgraph_induce_edge.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This commodity file contains the edge   **/
/**                arrays building subroutine which is     **/
/**                duplicated, with minor modifications,   **/
/**                into hgraph_induce.c                    **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 10 jan 2002     **/
/**                                 to     17 jan 2003     **/
/**                # Version 5.0  : from : 19 dec 2006     **/
/**                                 to     19 dec 2006     **/
/**                # Version 5.1  : from : 24 oct 2010     **/
/**                                 to     24 oct 2010     **/
/**                # Version 6.0  : from : 22 mar 2012     **/
/**                                 to     07 nov 2014     **/
/**                                                        **/
/************************************************************/

static
void
HGRAPHINDUCE2NAME (
const Hgraph * restrict const orggrafptr,         /* Pointer to original halo graph        */
Gnum * restrict const         orgindxtax,         /* Array of numbers of selected vertices */
Hgraph * restrict const       indgrafptr)         /* Pointer to induced halo graph         */
{
  Gnum                indvertnum;                 /* Number of current induced vertex                */
  Gnum                indvertnnd;                 /* Number of after-last induced (halo) vertex      */
  Gnum                indvelosum;                 /* Overall induced vertex load                     */
  Gnum                indedgenum;                 /* Number of current edge in induced halo subgraph */
  Gnum                indenohnbr;                 /* Number of non-halo edges in halo subgraph       */
  Gnum                inddegrmax;                 /* Maximum degree                                  */
#ifdef SCOTCH_DEBUG_HGRAPH2
  Gnum                indedgenbs;                 /* Revised number of edges in halo subgraph        */
#endif /* SCOTCH_DEBUG_HGRAPH2 */
#ifdef HGRAPHINDUCE2L                             /* If edge loads present */
  Gnum                indedlosum;
  Gnum                indenohsum;
#endif /* HGRAPHINDUCE2L */

  const Gnum * restrict const orgverttax = orggrafptr->s.verttax;
  const Gnum * restrict const orgvendtax = orggrafptr->s.vendtax;
  const Gnum * restrict const orgvelotax = orggrafptr->s.velotax;
  const Gnum * restrict const orgvnumtax = orggrafptr->s.vnumtax;
  const Gnum * restrict const orgedgetax = orggrafptr->s.edgetax;
  Gnum * restrict const       indvnhdtax = indgrafptr->vnhdtax;
  Gnum * restrict const       indverttax = indgrafptr->s.verttax;
  Gnum * restrict const       indvelotax = indgrafptr->s.velotax;
  Gnum * restrict const       indvnumtax = indgrafptr->s.vnumtax;
  Gnum * restrict const       indedgetax = indgrafptr->s.edgetax;
#ifdef HGRAPHINDUCE2L                             /* If edge loads present */
  const Gnum * restrict const orgedlotax = orggrafptr->s.edlotax;
  Gnum * restrict             indedlotax = indgrafptr->s.edlotax; /* Not const because location will change */

  indedlosum =
  indenohsum = 0;
#endif /* HGRAPHINDUCE2L */
  inddegrmax = 0;
  for (indvertnum = indedgenum = indgrafptr->s.baseval, indvelosum = indenohnbr = 0, indvertnnd = indgrafptr->vnohnnd; /* For all non-halo vertices */
       indvertnum < indgrafptr->vnohnnd; indvertnum ++) {
    Gnum                orgvertnum;               /* Number of current vertex in original halo graph       */
    Gnum                orgedgenum;               /* Number of current edge in original halo graph         */
    Gnum                indedgennd;               /* Index of after-last edge position in edge array       */
    Gnum                indedhdnum;               /* Index of after-last edge linking to non-halo vertices */
    Gnum                inddegrval;

    orgvertnum = indvnumtax[indvertnum];
    indverttax[indvertnum] = indedgenum;
    indenohnbr -= indedgenum;                     /* Subtract base of non-halo edges */
    if (indvelotax != NULL) {                     /* If graph has vertex weights     */
      indvelosum +=                               /* Accumulate vertex loads         */
      indvelotax[indvertnum] = orgvelotax[orgvertnum];
    }

    inddegrval = orgvendtax[orgvertnum] - orgverttax[orgvertnum]; /* Get degree of non-halo node */
    if (inddegrmax < inddegrval)                  /* Keep maximum degree                         */
      inddegrmax = inddegrval;

    for (orgedgenum = orgverttax[orgvertnum], indedhdnum = indedgennd = indedgenum + inddegrval;
         orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
      Gnum                orgvertend;             /* Number of current end vertex in original halo graph   */
      Gnum                indvertend;             /* Number of current end vertex in induced halo subgraph */

      orgvertend = orgedgetax[orgedgenum];
      indvertend = orgindxtax[orgvertend];
      if (indvertend == ~0) {                     /* If neighbor is yet undeclared halo vertex */
        indvnumtax[indvertnnd] = orgvertend;      /* Add number of halo vertex to array        */
        indvertend = orgindxtax[orgvertend] = indvertnnd ++; /* Get induced number of vertex   */
      }
      if (indvertend >= indgrafptr->vnohnnd) {    /* If neighbor is halo vertex            */
        indedhdnum --;                            /* Add neighbor at end of edge sub-array */
        indedgetax[indedhdnum] = indvertend;
        HGRAPHINDUCE2EDLOINIT (indedhdnum);
      }
      else {                                      /* If heighbor is non-halo vertex              */
        indedgetax[indedgenum] = indvertend;      /* Add neighbor at beginning of edge sub-array */
        HGRAPHINDUCE2EDLOINIT (indedgenum);
        HGRAPHINDUCE2ENOHINIT;
        indedgenum ++;
      }
    }
#ifdef SCOTCH_DEBUG_HGRAPH2
    if (indedgenum != indedhdnum) {
      errorPrint (STRINGIFY (HGRAPHINDUCE2NAME) ": internal error (1)");
      return;
    }
#endif /* SCOTCH_DEBUG_HGRAPH2 */
    indenohnbr += indedhdnum;                     /* Add position to number of non-halo edges */
    indvnhdtax[indvertnum] = indedhdnum;          /* Set end of non-halo sub-array            */
    indedgenum = indedgennd;                      /* Point to next free space in edge array   */
  }
  indgrafptr->vnlosum = (indvelotax != NULL) ? indvelosum : indgrafptr->vnohnbr;
  indgrafptr->enohnbr = indenohnbr;

#ifdef SCOTCH_DEBUG_HGRAPH2
  indedgenbs = 2 * (indedgenum - indgrafptr->s.baseval) - indenohnbr; /* Compute total number of edges */
#endif /* SCOTCH_DEBUG_HGRAPH2 */
#ifdef HGRAPHINDUCE2L                             /* If edge loads present */
  {
    Gnum *              indedgetab;               /* Dummy area to recieve un-based edgetab */
    Gnum *              indedlotab;               /* Save of old position of edgetab array  */
#ifndef SCOTCH_DEBUG_HGRAPH2
    Gnum                indedgenbs;               /* Revised number of edges in halo subgraph */

    indedgenbs = 2 * (indedgenum - indgrafptr->s.baseval) - indenohnbr; /* Compute total number of edges */
#endif /* SCOTCH_DEBUG_HGRAPH2 */

    indedlotab = indedlotax + indgrafptr->s.baseval; /* Save old offset of move area */
    memOffset (indedgetax + indgrafptr->s.baseval, /* Compute new offsets            */
               &indedgetab, (size_t) (indedgenbs * sizeof (Gnum)),
               &indedlotax, (size_t) (indedgenbs * sizeof (Gnum)), NULL);
    memMov (indedlotax, indedlotab, (indedgenum - indgrafptr->s.baseval) * sizeof (Gnum)); /* Move already existing edge load array */
    indgrafptr->s.edlotax =                       /* Record new position of edge load array */
    indedlotax           -= indgrafptr->s.baseval;
  }
#endif /* HGRAPHINDUCE2L */

  for ( ; indvertnum < indvertnnd; indvertnum ++) { /* For all halo vertices found during first pass */
    Gnum                orgvertnum;               /* Number of current vertex in original halo graph */
    Gnum                orgedgenum;               /* Number of current edge in original halo graph   */

    orgvertnum = indvnumtax[indvertnum];
    indverttax[indvertnum] = indedgenum;
    if (indvelotax != NULL) {                     /* If graph has vertex weights */
      indvelosum +=                               /* Accumulate vertex loads     */
      indvelotax[indvertnum] = orgvelotax[orgvertnum];
    }

    for (orgedgenum = orgverttax[orgvertnum];
         orgedgenum < orgvendtax[orgvertnum]; orgedgenum ++) {
      Gnum                orgvertend;             /* Number of current end vertex in original halo graph   */
      Gnum                indvertend;             /* Number of current end vertex in induced halo subgraph */

      orgvertend = orgedgetax[orgedgenum];
      indvertend = orgindxtax[orgvertend];
      if ((indvertend != ~0) &&                   /* If end vertex in induced halo subgraph */
          (indvertend < indgrafptr->vnohnnd)) {   /* And in its non-halo part only          */
        indedgetax[indedgenum] = indvertend;
        HGRAPHINDUCE2EDLOINIT (indedgenum);
        indedgenum ++;
      }
    }
    if (inddegrmax < (indedgenum - indverttax[indvertnum]))
      inddegrmax = (indedgenum - indverttax[indvertnum]);
  }
#ifdef SCOTCH_DEBUG_HGRAPH2
  if ((indedgenum - indgrafptr->s.baseval) != indedgenbs) {
    errorPrint (STRINGIFY (HGRAPHINDUCE2NAME) ": internal error (2)");
    return;
  }
#endif /* SCOTCH_DEBUG_HGRAPH2 */
  indverttax[indvertnnd] = indedgenum; /* Set end of compact vertex array */

  indgrafptr->s.vertnbr = indvertnnd - indgrafptr->s.baseval;
  indgrafptr->s.vertnnd = indvertnnd;
  indgrafptr->s.velosum = (indvelotax != NULL) ? indvelosum : indgrafptr->s.vertnbr;
  indgrafptr->s.edgenbr = indedgenum - indgrafptr->s.baseval; /* Set actual number of edges */
  indgrafptr->s.edlosum = HGRAPHINDUCE2EDLOSUM;
  indgrafptr->s.degrmax = inddegrmax;
  indgrafptr->enohsum   = HGRAPHINDUCE2ENOHSUM;
}
