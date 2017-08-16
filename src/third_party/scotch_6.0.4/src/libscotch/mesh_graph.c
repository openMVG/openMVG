/* Copyright 2004,2007,2009 ENSEIRB, INRIA & CNRS
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
/**   NAME       : mesh_graph.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source mesh    **/
/**                to graph conversion function.           **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 11 oct 2003     **/
/**                                 to     05 may 2004     **/
/**                # Version 5.1  : from : 19 nov 2009     **/
/**                                 to     19 nov 2009     **/
/**                                                        **/
/**   NOTES      : # From a given mesh is created a graph, **/
/**                  such that all vertices of the graph   **/
/**                  represent the nodes of the mesh, and  **/
/**                  there exists an edge between two      **/
/**                  vertices if there exists at least one **/
/**                  element to which the two associated   **/
/**                  nodes belong.                         **/
/**                  In order to extract mesh vertex       **/
/**                  partitions from graph vertex          **/
/**                  partitions easily, the vertices of    **/
/**                  the graph are numbered in the same    **/
/**                  order as the nodes of the mesh.       **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH_GRAPH

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "mesh_graph.h"

/*******************************/
/*                             */
/* The graph building routine. */
/*                             */
/*******************************/

/* This routine builds a graph from the
** given mesh.
** It returns:
** - 0  : if the graph has been successfully built.
** - 1  : on error.
*/

int
meshGraph (
const Mesh * restrict const   meshptr,            /*+ Original mesh  +*/
Graph * restrict const        grafptr)            /*+ Graph to build +*/
{
  Gnum                        hashnbr;            /* Number of vertices in hash table       */
  Gnum                        hashsiz;            /* Size of hash table                     */
  Gnum                        hashmsk;            /* Mask for access to hash table          */
  MeshGraphHash * restrict    hashtab;            /* Table of edges to other node vertices  */
  Gnum                        edgemax;            /* Upper bound of number of edges in mesh */
  Gnum                        edgennd;            /* Based upper bound on number of edges   */
  Gnum                        edgenum;            /* Number of current graph edge           */
  Gnum                        vertnum;            /* Number of current graph vertex         */
  Gnum                        degrmax;

  grafptr->flagval = GRAPHFREETABS | GRAPHVERTGROUP | GRAPHEDGEGROUP;
  grafptr->baseval = meshptr->baseval;
  grafptr->vertnbr = meshptr->vnodnbr;
  grafptr->vertnnd = meshptr->vnodnbr + meshptr->baseval;

  for (hashsiz = 32, hashnbr = meshptr->degrmax * meshptr->degrmax * 2; /* Compute size of hash table */
       hashsiz < hashnbr; hashsiz <<= 1) ;
  hashmsk = hashsiz - 1;

  if (((grafptr->verttax = memAlloc ((meshptr->vnodnbr + 1) * sizeof (Gnum)))          == NULL) ||
      ((hashtab          = memAlloc (hashsiz                * sizeof (MeshGraphHash))) == NULL)) {
    errorPrint ("meshGraph: out of memory (1)");
    if (grafptr->verttax != NULL)
      memFree (grafptr->verttax);
    return (1);
  }
  grafptr->verttax -= grafptr->baseval;
  grafptr->vendtax  = grafptr->verttax + 1;

  if (meshptr->vnlotax != NULL)                   /* Keep node part of mesh vertex load array as graph vertex load array        */
    grafptr->velotax = meshptr->vnlotax + meshptr->vnodbas - grafptr->baseval; /* Since GRAPHVERTGROUP, no problem on graphFree */

  grafptr->velosum = meshptr->vnlosum;

  edgemax = 2 * meshptr->edgenbr;                 /* Compute lower bound on number of edges in graph */
#ifdef SCOTCH_DEBUG_MESH2
  edgemax = meshptr->degrmax + 1;                 /* Allow testing dynamic reallocation of edge array */
#endif /* SCOTCH_DEBUG_MESH2 */

  if ((grafptr->edgetax = memAlloc (edgemax * sizeof (Gnum))) == NULL) {
    errorPrint ("meshGraph: out of memory (2)");
    graphFree  (grafptr);
    return     (1);
  }
  grafptr->edgetax -= grafptr->baseval;

  memSet (hashtab, ~0, hashsiz * sizeof (MeshGraphHash)); /* Initialize hash table */

  for (vertnum = edgenum = grafptr->baseval, edgennd = edgemax + grafptr->baseval, degrmax = 0; /* Build graph edges */
       vertnum < grafptr->vertnnd; vertnum ++) {
    Gnum                        vnodnum;
    Gnum                        hnodnum;
    Gnum                        enodnum;

    grafptr->verttax[vertnum] = edgenum;

    vnodnum = vertnum + (meshptr->vnodbas - meshptr->baseval);
    hnodnum = (vnodnum * MESHGRAPHHASHPRIME) & hashmsk; /* Prevent adding loop edge */
    hashtab[hnodnum].vertnum = vnodnum;
    hashtab[hnodnum].vertend = vnodnum;

    for (enodnum = meshptr->verttax[vnodnum]; enodnum < meshptr->vendtax[vnodnum]; enodnum ++) {
      Gnum                        velmnum;
      Gnum                        eelmnum;

      velmnum = meshptr->edgetax[enodnum];

      for (eelmnum = meshptr->verttax[velmnum]; eelmnum < meshptr->vendtax[velmnum]; eelmnum ++) {
        Gnum                        vnodend;
        Gnum                        hnodend;

        vnodend = meshptr->edgetax[eelmnum];

        for (hnodend = (vnodend * MESHGRAPHHASHPRIME) & hashmsk; ; hnodend = (hnodend + 1) & hashmsk) {
          if (hashtab[hnodend].vertnum != vnodnum) { /* If edge not yet created */
            if (edgenum == edgennd) {             /* If edge array already full */
              Gnum                        edgemax;
              Gnum * restrict             edgetmp;

              edgemax = edgennd - grafptr->baseval; /* Increase size by 25 % */
              edgemax = edgemax + (edgemax >> 2);

              if ((edgetmp = memRealloc (grafptr->edgetax + grafptr->baseval, edgemax * sizeof (Gnum))) == NULL) {
                errorPrint ("meshGraph: out of memory (3)");
                graphFree  (grafptr);
                memFree    (hashtab);
                return     (1);
              }

              grafptr->edgetax = edgetmp - grafptr->baseval;
              edgennd          = edgemax + grafptr->baseval;
            }

            hashtab[hnodend].vertnum = vnodnum;   /* Record new edge */
            hashtab[hnodend].vertend = vnodend;
            grafptr->edgetax[edgenum ++] = vnodend - (meshptr->vnodbas - grafptr->baseval); /* Build new edge */
            break;
          }
          if (hashtab[hnodend].vertend == vnodend) /* If edge already exists */
            break;                                /* Skip to next neighbor   */
        }
      }
    }

    if ((edgenum - grafptr->verttax[vertnum]) > degrmax) /* Compute maximum degree */
      degrmax = (edgenum - grafptr->verttax[vertnum]);
  }
  grafptr->verttax[vertnum] = edgenum;            /* Set end of vertex array */

  grafptr->edgenbr = edgenum - grafptr->baseval;
  grafptr->degrmax = degrmax;

#ifdef SCOTCH_DEBUG_MESH2
  if (graphCheck (grafptr) != 0) {
    errorPrint ("meshGraph: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  memFree(hashtab);

  return (0);
}
