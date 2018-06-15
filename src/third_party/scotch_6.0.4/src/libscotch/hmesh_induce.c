/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
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
/**   NAME       : hmesh_induce.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the halo source     **/
/**                mesh subgraph-making functions.         **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 07 jan 2002     **/
/**                                 to     11 may 2004     **/
/**                # Version 5.0  : from : 22 dec 2006     **/
/**                                 to     11 jun 2007     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH
#define HMESH_INDUCE

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "hmesh.h"

/*********************************************/
/*                                           */
/* These routines handle halo source meshes. */
/*                                           */
/*********************************************/

/* This routine builds the halo mesh induced
** by the original halo mesh and the array of
** selected node vertices. Elements which are
** adjacent to the selected nodes are themselves
** selected, as well as their adjacent node
** vertices, which comprise the new halo.
** In the induced halo mesh, elements are
** placed first, then non-halo nodes, then
** halo nodes. This order is quite important
** as it eases the building of vertex separation
** meshes from halo meshes, just by ignoring
** halo nodes.
** The induced vnumtab array is a baseval-based
** list of the selected node vertices if the
** original halo mesh does not have a vnumtab, or
** the proper subset of the original vnumtab else.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hmeshInducePart (
const Hmesh * restrict const      orgmeshptr,     /* Pointer to original graph                                    */
const GraphPart * restrict const  orgparttax,     /* Array of vertex partition flags                              */
const GraphPart                   orgpartval,     /* Partition value of vertices to keep (0 or 1)                 */
const Gnum                        orgvelmnbr,     /* Number of (maybe isolated) element vertices in selected part */
const Gnum                        orgvnodnbr,     /* Number of node vertices in selected part                     */
const Gnum                        orgvnspnbr,     /* Number of node vertices in separator                         */
Hmesh * restrict const            indmeshptr)     /* Pointer to induced halo submesh                              */
{
  Gnum                orgvelmnum;                 /* Number of current element vertex in original halo mesh          */
  Gnum                orgvnodnum;                 /* Number of current node vertex in original halo mesh             */
  Gnum * restrict     orgindxtax;                 /* Original to induced vertex number translation array             */
  Gnum                indvertnbr;                 /* Upper bound on the number of vertices in induced halo mesh      */
  Gnum                indvertnum;                 /* Number of current vertex in induced mesh graph                  */
  Gnum                indvelmnbr;                 /* Number of element vertices in induced halo mesh                 */
  Gnum                indveihnbr;                 /* Number of newly created halo-isolated elements                  */
  Gnum                indvnodnbr;                 /* Upper bound on the number of node vertices in induced halo mesh */
  Gnum                indvnodnum;                 /* Number of current node vertex in induced halo mesh              */
  Gnum * restrict     indedgetax;                 /* Based access to induced mesh graph edge arrays                  */
  Gnum                indedgenbr;                 /* (Approximate) number of edges in induced halo mesh              */
  Gnum                indedgenum;                 /* Number of current edge in induced halo mesh                     */
  Gnum * restrict     indvnuhtax;                 /* Array of vertex numbers for halo nodes (aka vnumtab)            */
  Gnum                indvelonbr;
  Gnum                indvelosum;
  Gnum                indvnlonbr;
  Gnum                indvnlosum;
  Gnum                indvnhlsum;

  indvelmnbr = orgvelmnbr - orgmeshptr->veihnbr;  /* Remove known halo-isolated elements               */
  if (orgpartval == 0)                            /* If submesh is part zero                           */
    indvelmnbr -= orgmeshptr->m.veisnbr;          /* Also remove isolated elements, which belong to it */

  indvnodnbr = orgvnodnbr + orgvnspnbr + orgmeshptr->m.vnodnnd - orgmeshptr->vnohnnd; /* Compute upper bound on number of node vertices */
  indvertnbr = indvnodnbr + indvelmnbr;

  indedgenbr = ((orgmeshptr->m.degrmax > 0) && (indvertnbr < (orgmeshptr->m.edgenbr / orgmeshptr->m.degrmax))) /* Choose best upper bound on number of edges (avoid multiply overflow) */
               ? (indvertnbr * orgmeshptr->m.degrmax) : orgmeshptr->m.edgenbr;

  memSet (indmeshptr, 0, sizeof (Hmesh));         /* Initialize halo mesh fields */
  indmeshptr->m.baseval = orgmeshptr->m.baseval;  /* Inherit mesh properties     */
  indmeshptr->m.flagval = MESHFREETABS | MESHVERTGROUP;
  indmeshptr->m.velmnbr = indvelmnbr;
  indmeshptr->m.velmbas = indmeshptr->m.baseval;  /* Elements are placed first */
  indmeshptr->m.velmnnd =
  indmeshptr->m.vnodbas = orgvelmnbr + indmeshptr->m.baseval; /* Node vertices are placed after elements */
  indmeshptr->m.veisnbr = 0;                      /* All isolated elements will be removed in submesh    */

  indvelonbr = (orgmeshptr->m.velotax != NULL) ? indmeshptr->m.velmnbr : 0;
  indvnlonbr = (orgmeshptr->m.vnlotax != NULL) ? indvnodnbr            : 0;

  if (memAllocGroup ((void **) (void *)
                     &indmeshptr->m.verttax, (size_t) ((indvertnbr + 1) * sizeof (Gnum)),
                     &indmeshptr->vehdtax,   (size_t) ( orgvelmnbr      * sizeof (Gnum)), /* vehdtab is limited to elements */
                     &indmeshptr->m.velotax, (size_t) ( indvelonbr      * sizeof (Gnum)),
                     &indmeshptr->m.vnlotax, (size_t) ( indvnlonbr      * sizeof (Gnum)),
                     &indmeshptr->m.vnumtax, (size_t) ( orgvnodnbr      * sizeof (Gnum)), NULL) == NULL) {  /* vnumtab is of size vnohnbr */
    errorPrint ("hmeshInducePart: out of memory (1)"); /* Allocate induced mesh graph structure */
    return     (1);
  }
  indmeshptr->m.verttax -= indmeshptr->m.baseval;
  indmeshptr->m.vendtax  = indmeshptr->m.verttax + 1; /* Compact array */
  indmeshptr->m.velotax  = (indvelonbr != 0) ? (indmeshptr->m.velotax - indmeshptr->m.velmbas) : NULL;
  indmeshptr->m.vnlotax  = (indvnlonbr != 0) ? (indmeshptr->m.vnlotax - indmeshptr->m.vnodbas) : NULL;
  indmeshptr->m.vnumtax -= indmeshptr->m.vnodbas; /* Only for non-halo nodes */
  indmeshptr->m.degrmax  = orgmeshptr->m.degrmax;
  indmeshptr->vnohnbr    = orgvnodnbr;
  indmeshptr->vnohnnd    = indmeshptr->m.vnodbas + orgvnodnbr;
  indmeshptr->vehdtax   -= indmeshptr->m.velmbas;
  indmeshptr->vnhlsum    = orgvnodnbr;            /* Assume no vertex loads */

  if (memAllocGroup ((void **) (void *)
                     &indedgetax, (size_t) (indedgenbr                                                 * sizeof (Gnum)),
                     &orgindxtax, (size_t) ((orgmeshptr->m.velmnbr + orgmeshptr->m.vnodnbr)            * sizeof (Gnum)),
                     &indvnuhtax, (size_t) ((orgvnspnbr + orgmeshptr->m.vnodnnd - orgmeshptr->vnohnnd) * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("hmeshInducePart: out of memory (2)"); /* Allocate induced mesh graph structure */
    hmeshExit  (indmeshptr);
    return     (1);
  }
  indedgetax -= indmeshptr->m.baseval;
  orgindxtax -= orgmeshptr->m.baseval;
  indvnuhtax -= indmeshptr->vnohnnd;              /* Base array so as to catch halo nodes only */

  indvnhlsum = 0;
  for (orgvnodnum = orgmeshptr->m.vnodbas,        /* For all original non-halo node vertices */
       indvnodnum = indmeshptr->m.vnodbas;        /* And assuming all elements will be kept  */
       orgvnodnum < orgmeshptr->vnohnnd; orgvnodnum ++) {
    if (orgparttax[orgvnodnum] == orgpartval) {   /* If in right part      */
      orgindxtax[orgvnodnum] = indvnodnum;        /* Set new index of node */
      indmeshptr->m.vnumtax[indvnodnum] = orgvnodnum - (orgmeshptr->m.vnodbas - orgmeshptr->m.baseval);
      if (orgmeshptr->m.vnlotax != NULL)
        indvnhlsum += (indmeshptr->m.vnlotax[indvnodnum] = orgmeshptr->m.vnlotax[orgvnodnum]);
      indvnodnum ++;
    }
    else if (orgparttax[orgvnodnum] == 2)         /* If node belongs to separator      */
      orgindxtax[orgvnodnum] = ~0;                /* Pre-set array for separator nodes */
  }
#ifdef SCOTCH_DEBUG_HMESH2
  if ((indvnodnum - indmeshptr->m.vnodbas) != orgvnodnbr) {
    errorPrint ("hmeshInducePart: internal error (1)");
    memFree    (indedgetax + indmeshptr->m.baseval);
    hmeshExit  (indmeshptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HMESH2 */
  memSet (orgindxtax + orgmeshptr->vnohnnd, ~0, (orgmeshptr->m.vnodnnd - orgmeshptr->vnohnnd) * sizeof (Gnum)); /* Pre-set halo node vertices */

  indveihnbr = 0;
  indvelosum = 0;
  indvnlosum = 0;
  for (orgvelmnum = orgmeshptr->m.velmbas,        /* For all elements of original graph              */
       indvertnum = indedgenum = indmeshptr->m.baseval; /* Elements are placed first in vertex array */
       orgvelmnum < orgmeshptr->m.velmnnd; orgvelmnum ++) {
    if (orgparttax[orgvelmnum] == orgpartval) {   /* If element belongs to right part */
      Gnum                orgedgenum;
      Gnum                indedgennd;             /* Index of after-last edge position in edge array       */
      Gnum                indedhdnum;             /* Index of after-last edge linking to non-halo vertices */

      orgedgenum = orgmeshptr->m.verttax[orgvelmnum];
      if (orgedgenum == orgmeshptr->vehdtax[orgvelmnum]) /* If (halo-)isolated element vertex */
        continue;                                 /* Discard element in induced submesh       */

#ifdef SCOTCH_DEBUG_HMESH2
      if (indvertnum >= indmeshptr->m.velmnnd) {  /* If too many element vertices kept                 */
        errorPrint ("hmeshInducePart: internal error (2)"); /* Maybe a problem with veisnbr or veihnbr */
        memFree    (indedgetax + indmeshptr->m.baseval);
        hmeshExit  (indmeshptr);
        return     (1);
      }
#endif /* SCOTCH_DEBUG_HMESH2 */

      indmeshptr->m.verttax[indvertnum] = indedgenum;
      indedhdnum = orgmeshptr->m.vendtax[orgvelmnum] - orgmeshptr->m.verttax[orgvelmnum] + indedgenum;
      indedgennd = indedhdnum;
      if (orgmeshptr->m.velotax != NULL) {
        Gnum                orgveloval;

        orgveloval = orgmeshptr->m.velotax[orgvelmnum];
        indmeshptr->m.velotax[indvertnum] = orgveloval;
        indvelosum += orgveloval;
      }

      for ( ; orgedgenum < orgmeshptr->m.vendtax[orgvelmnum]; orgedgenum ++) {
        Gnum                orgvertend;

        orgvertend = orgmeshptr->m.edgetax[orgedgenum];

        if (orgindxtax[orgvertend] == ~0) {       /* If found yet un-numbered halo node */
#ifdef SCOTCH_DEBUG_HMESH2
          if ((orgvertend < orgmeshptr->vnohnnd) &&
              (orgparttax[orgvertend] != 2)) {
            errorPrint ("hmeshInducePart: internal error (3)");
            memFree    (indedgetax + indmeshptr->m.baseval);
            hmeshExit  (indmeshptr);
            return     (1);
          }
#endif /* SCOTCH_DEBUG_HMESH2 */
          if (orgmeshptr->m.vnlotax != NULL) {
            Gnum                orgvnloval;

            orgvnloval = orgmeshptr->m.vnlotax[orgvertend];
            indmeshptr->m.vnlotax[indvnodnum] = orgvnloval;
            indvnlosum += orgvnloval;
          }
          orgindxtax[orgvertend] = indvnodnum;    /* Set number of halo node  */
          indvnuhtax[indvnodnum] = orgvertend;    /* Keep number of halo node */
          indvnodnum ++;
          indedgetax[-- indedhdnum] = orgindxtax[orgvertend];
          continue;
        }
        if (orgindxtax[orgvertend] < indmeshptr->vnohnnd) /* If non-halo vertex */
          indedgetax[indedgenum ++] = orgindxtax[orgvertend];
        else                                      /* Else if halo vertex */
          indedgetax[-- indedhdnum] = orgindxtax[orgvertend];
      }
#ifdef SCOTCH_DEBUG_HMESH2
      if (indedgenum != indedhdnum) {
        errorPrint ("hmeshInducePart: internal error (4)");
        memFree    (indedgetax + indmeshptr->m.baseval);
        hmeshExit  (indmeshptr);
        return     (1);
      }
#endif /* SCOTCH_DEBUG_HMESH2 */
      if (indedhdnum == indmeshptr->m.verttax[indvertnum]) /* If element has halo nodes only */
        indveihnbr ++;                            /* One more halo-isolated element created  */

      indmeshptr->vehdtax[indvertnum] = indedhdnum;
      indedgenum = indedgennd;
      orgindxtax[orgvelmnum] = indvertnum;
      indvertnum ++;                              /* One more element created */
    }
    else
      orgindxtax[orgvelmnum] = ~0;
  }
#ifdef SCOTCH_DEBUG_HMESH2
  if (indvertnum != indmeshptr->m.velmnnd) {
    errorPrint ("hmeshInducePart: internal error (5)"); /* Maybe a problem with veisnbr or veihnbr */
    memFree    (indedgetax + indmeshptr->m.baseval);
    hmeshExit  (indmeshptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HMESH2 */

  indmeshptr->veihnbr = indveihnbr;

  indmeshptr->m.vnodnbr = indvnodnum - indmeshptr->m.vnodbas;
  indmeshptr->m.vnodnnd = indvertnum + indmeshptr->m.vnodnbr;
  indmeshptr->m.velosum = (indmeshptr->m.velotax != NULL) ? indvelosum : indmeshptr->m.velmnbr;
  if (indmeshptr->m.vnlotax != NULL) {            /* If vertex loads wanted */
    indmeshptr->m.vnlosum = indvnhlsum + indvnlosum;
    indmeshptr->vnhlsum   = indvnhlsum;
  }
  else {
    indmeshptr->m.vnlosum = indmeshptr->m.vnodnbr;
    indmeshptr->vnhlsum   = indmeshptr->vnohnbr;
  }

  indedgenbr = 2 * (indedgenum - indmeshptr->m.baseval); /* Twice as many arcs as element arcs         */
  for ( ; indvertnum < indmeshptr->vnohnnd; indvertnum ++) { /* For all non-halo induced node vertices */
    Gnum                orgvnodnum;
    Gnum                orgedgenum;

    orgvnodnum = indmeshptr->m.vnumtax[indvertnum] + (orgmeshptr->m.vnodbas - orgmeshptr->m.baseval); /* Get number of original node */

    indmeshptr->m.verttax[indvertnum] = indedgenum;

    for (orgedgenum = orgmeshptr->m.verttax[orgvnodnum];
         orgedgenum < orgmeshptr->m.vendtax[orgvnodnum]; orgedgenum ++) {
      Gnum                orgvertend;

      orgvertend = orgmeshptr->m.edgetax[orgedgenum];
#ifdef SCOTCH_DEBUG_HMESH2
      if (orgindxtax[orgvertend] == ~0) {
        errorPrint ("hmeshInducePart: internal error (6)");
        memFree    (indedgetax + indmeshptr->m.baseval);
        hmeshExit  (indmeshptr);
        return     (1);
      }
#endif /* SCOTCH_DEBUG_HMESH2 */
      indedgetax[indedgenum ++] = orgindxtax[orgvertend];
    }
  }

  indmeshptr->enohnbr = indedgenum - indmeshptr->m.baseval;

  for ( ; indvertnum < indmeshptr->m.vnodnnd; indvertnum ++) { /* For all halo induced node vertices */
    Gnum                orgvnodnum;
    Gnum                orgedgenum;

    orgvnodnum = indvnuhtax[indvertnum];        /* Get number of original node */

    indmeshptr->m.verttax[indvertnum] = indedgenum;

    for (orgedgenum = orgmeshptr->m.verttax[orgvnodnum];
         orgedgenum < orgmeshptr->m.vendtax[orgvnodnum]; orgedgenum ++) {
      Gnum                orgvertend;

      orgvertend = orgmeshptr->m.edgetax[orgedgenum];
      if (orgindxtax[orgvertend] != ~0) {       /* If end element belongs to right part */
        indedgetax[indedgenum ++] = orgindxtax[orgvertend];
      }
    }
  }
  indmeshptr->m.verttax[indvertnum] = indedgenum; /* Set end of edge array */
  indmeshptr->m.edgenbr = indedgenum - indmeshptr->m.baseval;
  indmeshptr->m.vnodnnd = indvertnum;             /* Record number of induced non-element vertices */
  indmeshptr->m.vnodnbr = indvertnum - indmeshptr->m.vnodbas;

  if (orgmeshptr->m.vnumtax != NULL) {          /* If source mesh is not original mesh */
    for (indvnodnum = indmeshptr->m.vnodbas; indvnodnum < indmeshptr->vnohnnd; indvnodnum ++)
      indmeshptr->m.vnumtax[indvnodnum] = orgmeshptr->m.vnumtax[indmeshptr->m.vnumtax[indvnodnum] + (orgmeshptr->m.vnodbas - orgmeshptr->m.baseval)];
  }

  indmeshptr->m.edgetax  = memRealloc (indedgetax + indmeshptr->m.baseval, indedgenbr * sizeof (Gnum));
  indmeshptr->m.edgetax -= indmeshptr->m.baseval;

#ifdef SCOTCH_DEBUG_HMESH2
  if (hmeshCheck (indmeshptr) != 0) {             /* Check halo mesh consistency */
    errorPrint ("hmeshInducePart: inconsistent halo mesh data");
    hmeshExit  (indmeshptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HMESH2 */

  return (0);
}
