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
/**   NAME       : vmesh_separate_gr.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates a node separation **/
/**                mesh by turning the mesh into a graph   **/
/**                and using a graph separation strategy.  **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 13 oct 2003     **/
/**                                 to     13 oct 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_SEPARATE_GR

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "vgraph.h"
#include "vgraph_separate_st.h"
#include "mesh.h"
#include "vmesh.h"
#include "vmesh_separate_gr.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine performs the bipartitioning.
** It returns:
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
vmeshSeparateGr (
Vmesh * restrict const                      meshptr, /*+ Node separation mesh +*/
const VmeshSeparateGrParam * restrict const paraptr) /*+ Method parameters    +*/
{
  Vgraph                            grafdat;
  Gnum                              fronnum;
  Gnum                              velmnum;
  Gnum                              ecmpsize1;

  graphInit (&grafdat.s);
  if (meshGraph (&meshptr->m, &grafdat.s) != 0) {
    errorPrint ("vmeshSeparateGr: cannot build graph");
    return     (1);
  }
  grafdat.parttax     = meshptr->parttax + (meshptr->m.vnodbas - grafdat.s.baseval); /* Get node area of part array */
  grafdat.compload[0] = meshptr->ncmpload[0];
  grafdat.compload[1] = meshptr->ncmpload[1];
  grafdat.compload[2] = meshptr->ncmpload[2];
  grafdat.comploaddlt = meshptr->ncmploaddlt;
  grafdat.compsize[0] = meshptr->ncmpsize[0];
  grafdat.compsize[1] = meshptr->ncmpsize[1];
  grafdat.fronnbr     = meshptr->fronnbr;
  grafdat.frontab     = meshptr->frontab;         /* Re-use frontier array */
  grafdat.levlnum     = meshptr->levlnum;

  for (fronnum = 0; fronnum < grafdat.fronnbr; fronnum ++)
    grafdat.frontab[fronnum] -= (meshptr->m.vnodbas - grafdat.s.baseval);

#ifdef SCOTCH_DEBUG_VMESH2
  if (vgraphCheck (&grafdat) != 0) {
    errorPrint ("vmeshSeparateGr: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  if (vgraphSeparateSt (&grafdat, paraptr->stratptr) != 0) {
    errorPrint ("vmeshSeparateGr: cannot separate graph");
    return     (1);
  }

  for (fronnum = 0; fronnum < grafdat.fronnbr; fronnum ++) /* Restore mesh-based frontier array */
    grafdat.frontab[fronnum] += (meshptr->m.vnodbas - grafdat.s.baseval);
  meshptr->ncmpload[0] = grafdat.compload[0];
  meshptr->ncmpload[1] = grafdat.compload[1];
  meshptr->ncmpload[2] = grafdat.compload[2];
  meshptr->ncmploaddlt = grafdat.comploaddlt;
  meshptr->ncmpsize[0] = grafdat.compsize[0];
  meshptr->ncmpsize[1] = grafdat.compsize[1];
  meshptr->fronnbr     = grafdat.fronnbr;

  for (velmnum = meshptr->m.velmbas, ecmpsize1 = 0;
       velmnum < meshptr->m.velmnnd; velmnum ++) { /* Compute part of all elements */
    Gnum                              eelmnum;
    GraphPart                         partval;

    partval = 0;                                  /* Empty elements move to part 0 */
    for (eelmnum = meshptr->m.verttax[velmnum];
         eelmnum < meshptr->m.vendtax[velmnum]; eelmnum ++) {
      Gnum                              vnodnum;

      vnodnum = meshptr->m.edgetax[eelmnum];
      partval = meshptr->parttax[vnodnum];
      if (partval != 2)
        break;
    }
    partval   &= 1;                               /* In case all nodes in separator */
    ecmpsize1 += (Gnum) partval;                  /* Count elements in part 1       */
    meshptr->parttax[velmnum] = partval;          /* Set part of element            */
  }
  meshptr->ecmpsize[0] = meshptr->m.velmnbr - ecmpsize1;
  meshptr->ecmpsize[1] = ecmpsize1;

#ifdef SCOTCH_DEBUG_VMESH2
  if (vmeshCheck (meshptr) != 0) {
    errorPrint ("vmeshSeparateGr: internal error (2)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  return (0);
}
