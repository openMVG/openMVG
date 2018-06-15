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
/**   NAME       : hmesh_hgraph.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source halo    **/
/**                mesh to halo graph conversion function. **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 30 nov 2003     **/
/**                                 to     05 may 2004     **/
/**                # Version 5.0  : from : 10 sep 2007     **/
/**                                 to   : 10 sep 2007     **/
/**                                                        **/
/**   NOTES      : # From a given halo mesh is created a   **/
/**                  halo graph, such that all vertices of **/
/**                  the graph represent the nodes of the  **/
/**                  mesh, and there exists an edge        **/
/**                  between two vertices if there exists  **/
/**                  at least one element to which the two **/
/**                  associated nodes belong.              **/
/**                  While all non-halo nodes become non-  **/
/**                  halo vertices, some halo nodes may    **/
/**                  disappear from the graph if their     **/
/**                  elements are only connected to other  **/
/**                  halo nodes.                           **/
/**                  Since the contents of vnumtab are     **/
/**                  based with respect to s.baseval and   **/
/**                  not to vnodbas, the vnumtab array can **/
/**                  simply be shared by the graph.        **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HMESH_HGRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "hgraph.h"
#include "mesh.h"
#include "hmesh.h"
#include "hmesh_hgraph.h"

/************************************/
/*                                  */
/* The halo graph building routine. */
/*                                  */
/************************************/

/* This routine builds a halo graph from
** the given halo mesh.
** It returns:
** - 0  : if the halo graph has been successfully built.
** - 1  : on error.
*/

int
hmeshHgraph (
const Hmesh * restrict const  meshptr,            /*+ Original mesh  +*/
Hgraph * restrict const       grafptr)            /*+ Graph to build +*/
{
  Gnum                        hashnbr;            /* Number of vertices in hash table       */
  Gnum                        hashsiz;            /* Size of hash table                     */
  Gnum                        hashmsk;            /* Mask for access to hash table          */
  HmeshHgraphHash * restrict  hashtab;            /* Table of edges to other node vertices  */
  Gnum                        edgemax;            /* Upper bound of number of edges in mesh */
  Gnum                        edgennd;            /* Based upper bound on number of edges   */
  Gnum                        enohnbr;            /* Number of non-halo edges               */
  Gnum                        edgenum;            /* Number of current graph edge           */
  Gnum                        vertnum;            /* Number of current graph vertex         */
  Gnum                        degrmax;

#ifdef SCOTCH_DEBUG_HMESH2
  if (hmeshCheck (meshptr) != 0) {
    errorPrint ("hmeshHgraph: invalid input halo mesh");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HMESH2 */

  grafptr->s.flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  grafptr->s.baseval = meshptr->m.baseval;
  grafptr->s.vertnbr = meshptr->m.vnodnbr;
  grafptr->s.vertnnd = meshptr->m.vnodnbr + meshptr->m.baseval;
  grafptr->vnohnbr   = meshptr->vnohnbr;
  grafptr->vnohnnd   = meshptr->vnohnbr + grafptr->s.baseval;
  grafptr->vnlosum   = meshptr->vnhlsum;

  for (hashsiz = 2, hashnbr = meshptr->m.degrmax * meshptr->m.degrmax * 2; /* Compute size of hash table */
       hashsiz < hashnbr; hashsiz <<= 1) ;
  hashmsk = hashsiz - 1;

  if (memAllocGroup ((void **) (void *)
                     &grafptr->s.verttax, (size_t) ((grafptr->s.vertnbr + 1) * sizeof (Gnum)),
                     &grafptr->vnhdtax,   (size_t) ( grafptr->vnohnbr        * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("hmeshHgraph: out of memory (1)");
    return     (1);
  }
  if ((hashtab = memAlloc (hashsiz * sizeof (HmeshHgraphHash))) == NULL) {
    errorPrint ("hmeshHgraph: out of memory (2)");
    memFree    (grafptr->s.verttax);
    return     (1);
  }
  grafptr->s.verttax -= grafptr->s.baseval;
  grafptr->s.vendtax  = grafptr->s.verttax + 1;
  grafptr->vnhdtax   -= grafptr->s.baseval;
  if (meshptr->m.vnumtax != NULL)                 /* If (node) vertex index array present, point to its relevant part                     */
    grafptr->s.vnumtax = meshptr->m.vnumtax + (meshptr->m.vnodbas - grafptr->s.baseval); /* Since GRAPHVERTGROUP, no problem on graphFree */

  if (meshptr->m.vnlotax != NULL)                 /* Keep node part of mesh (node) vertex load array as graph vertex load array           */
    grafptr->s.velotax = meshptr->m.vnlotax + (meshptr->m.vnodbas - grafptr->s.baseval); /* Since GRAPHVERTGROUP, no problem on graphFree */

  grafptr->s.velosum = meshptr->m.vnlosum;

  edgemax = ((meshptr->m.degrmax * meshptr->m.degrmax) / 2 + 1) * meshptr->m.vnodnbr; /* Compute estimated number of edges in graph */
#ifdef SCOTCH_DEBUG_HMESH2
  edgemax = meshptr->m.degrmax + 4;               /* Test dynamic reallocation of edge array; 4 guarantees that 25% > 0 */
#endif /* SCOTCH_DEBUG_HMESH2 */

  if ((grafptr->s.edgetax = memAlloc (edgemax * sizeof (Gnum))) == NULL) {
    errorPrint ("hmeshHgraph: out of memory (3)");
    hgraphFree (grafptr);
    return     (1);
  }
  grafptr->s.edgetax -= grafptr->s.baseval;

  memSet (hashtab, ~0, hashsiz * sizeof (HmeshHgraphHash)); /* Initialize hash table */

  for (vertnum = edgenum = grafptr->s.baseval, edgennd = edgemax + grafptr->s.baseval, enohnbr = degrmax = 0; /* Build graph edges for non-halo vertices */
       vertnum < grafptr->vnohnnd; vertnum ++) {
    Gnum                        vnodnum;
    Gnum                        hnodnum;
    Gnum                        enodnum;
    Gnum                        enhdnum;          /* Index of first non-halo neighbor in edge array for current vertex */

    grafptr->s.verttax[vertnum] = edgenum;

    vnodnum = vertnum + (meshptr->m.vnodbas - meshptr->m.baseval);
    hnodnum = (vnodnum * HMESHHGRAPHHASHPRIME) & hashmsk; /* Prevent adding loop edge */
    hashtab[hnodnum].vertnum = vnodnum;
    hashtab[hnodnum].vertend = vnodnum;

    for (enodnum = meshptr->m.verttax[vnodnum], enhdnum = edgenum;
         enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
      Gnum                        velmnum;
      Gnum                        eelmnum;

      velmnum = meshptr->m.edgetax[enodnum];

      for (eelmnum = meshptr->m.verttax[velmnum];
           eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
        Gnum                        vnodend;
        Gnum                        hnodend;

        vnodend = meshptr->m.edgetax[eelmnum];

        for (hnodend = (vnodend * HMESHHGRAPHHASHPRIME) & hashmsk; ; hnodend = (hnodend + 1) & hashmsk) {
          if (hashtab[hnodend].vertnum != vnodnum) { /* If edge not yet created */
            Gnum                        vertend;

            if (edgenum == edgennd) {             /* If edge array already full */
              Gnum                        edgemax;
              Gnum * restrict             edgetmp;

              edgemax = edgennd - grafptr->s.baseval; /* Increase size by 25 % */
              edgemax = edgemax + (edgemax >> 2);

              if ((edgetmp = memRealloc (grafptr->s.edgetax + grafptr->s.baseval, edgemax * sizeof (Gnum))) == NULL) {
                errorPrint ("hmeshHgraph: out of memory (4)");
                hgraphFree (grafptr);
                memFree    (hashtab);
                return     (1);
              }

              grafptr->s.edgetax = edgetmp - grafptr->s.baseval;
              edgennd            = edgemax + grafptr->s.baseval;
            }

            hashtab[hnodend].vertnum = vnodnum;   /* Record new edge */
            hashtab[hnodend].vertend = vnodend;
            vertend = vnodend - (meshptr->m.vnodbas - grafptr->s.baseval);
            if (vnodend >= meshptr->vnohnnd)      /* If halo edge                       */
              grafptr->s.edgetax[edgenum ++] = vertend; /* Build at end of array        */
            else {                                /* If non-halo edge                   */
              if (edgenum != enhdnum)             /* If already halo edges              */
                grafptr->s.edgetax[edgenum] = grafptr->s.edgetax[enhdnum]; /* Make room */
              grafptr->s.edgetax[enhdnum ++] = vertend; /* Record new edge              */
              edgenum ++;                         /* One more edge created              */
            }
            break;
          }
          if (hashtab[hnodend].vertend == vnodend) /* If edge already exists */
            break;                                /* Skip to next neighbor   */
        }
      }
    }
    grafptr->vnhdtax[vertnum] = enhdnum;          /* Set end of non-halo edge array */
    enohnbr += enhdnum - grafptr->s.verttax[vertnum];

    if ((edgenum - grafptr->s.verttax[vertnum]) > degrmax) /* Compute maximum degree */
      degrmax = (edgenum - grafptr->s.verttax[vertnum]);
  }
  grafptr->enohnbr = enohnbr;                     /* All other edges will be halo edges */

  for ( ; vertnum < grafptr->s.vertnnd; vertnum ++) { /* Build graph edges for halo vertices */
    Gnum                        vnodnum;
    Gnum                        enodnum;

    vnodnum = vertnum + (meshptr->m.vnodbas - meshptr->m.baseval);
    grafptr->s.verttax[vertnum] = edgenum;

    for (enodnum = meshptr->m.verttax[vnodnum];
         enodnum < meshptr->m.vendtax[vnodnum]; enodnum ++) {
      Gnum                        velmnum;
      Gnum                        eelmnum;

      velmnum = meshptr->m.edgetax[enodnum];

      for (eelmnum = meshptr->m.verttax[velmnum]; /* Only consider non-halo edges of elements this time */
           eelmnum < meshptr->vehdtax[velmnum]; eelmnum ++) {
        Gnum                        vnodend;
        Gnum                        hnodend;

        vnodend = meshptr->m.edgetax[eelmnum];

#ifdef SCOTCH_DEBUG_HMESH2
        if (vnodend >= meshptr->vnohnnd) {        /* Not visiting halo edges should prevent this */
          errorPrint ("hmeshHgraph: internal error (1)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_HMESH2 */

        for (hnodend = (vnodend * HMESHHGRAPHHASHPRIME) & hashmsk; ; hnodend = (hnodend + 1) & hashmsk) {
          if (hashtab[hnodend].vertnum != vnodnum) { /* If edge not yet created */
            Gnum                        vertend;

            if (edgenum == edgennd) {             /* If edge array already full */
              Gnum                        edgemax;
              Gnum * restrict             edgetmp;

              edgemax = edgennd - grafptr->s.baseval; /* Increase size by 25 % */
              edgemax = edgemax + (edgemax >> 2);

              if ((edgetmp = memRealloc (grafptr->s.edgetax + grafptr->s.baseval, edgemax * sizeof (Gnum))) == NULL) {
                errorPrint ("hmeshHgraph: out of memory (5)");
                hgraphFree (grafptr);
                memFree    (hashtab);
                return     (1);
              }

              grafptr->s.edgetax = edgetmp - grafptr->s.baseval;
              edgennd            = edgemax + grafptr->s.baseval;
            }

            hashtab[hnodend].vertnum = vnodnum;   /* Record new edge */
            hashtab[hnodend].vertend = vnodend;
            vertend = vnodend - (meshptr->m.vnodbas - grafptr->s.baseval);
            grafptr->s.edgetax[edgenum ++] = vertend; /* Build halo edge */
            break;
          }
          if (hashtab[hnodend].vertend == vnodend) /* If edge already exists */
            break;                                /* Skip to next neighbor   */
        }
      }
    }

    if ((edgenum - grafptr->s.verttax[vertnum]) > degrmax) /* Compute maximum degree */
      degrmax = (edgenum - grafptr->s.verttax[vertnum]);
  }
  grafptr->s.verttax[vertnum] = edgenum;          /* Set end of vertex array */
  grafptr->s.edgenbr = edgenum - grafptr->s.baseval;
  grafptr->s.degrmax = degrmax;

  memFree (hashtab);

#ifdef SCOTCH_DEBUG_HMESH2
  if (hgraphCheck (grafptr) != 0) {
    errorPrint ("hmeshHgraph: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_HMESH2 */

  return (0);
}
