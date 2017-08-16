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
/**   NAME       : vmesh_separate_ml.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module separates an active         **/
/**                mesh using a multi-level scheme.        **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 feb 2003     **/
/**                                 to     31 aug 2005     **/
/**                # Version 5.0  : from : 30 jan 2008     **/
/**                                 to     30 jan 2008     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define VMESH_SEPARATE_ML

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "mesh.h"
#include "mesh_coarsen.h"
#include "vmesh.h"
#include "vmesh_separate_ml.h"
#include "vmesh_separate_st.h"

/*********************************************/
/*                                           */
/* The coarsening and uncoarsening routines. */
/*                                           */
/*********************************************/

/* This routine builds a coarser mesh from the
** mesh that is given on input. The coarser
** meshes differ at this stage from classical
** active meshes as their internal gains are not
** yet computed.
** It returns:
** - 0  : if the coarse mesh has been built.
** - 1  : if threshold achieved.
** - 2  : on error.
*/

static
int
vmeshSeparateMlCoarsen (
const Vmesh * restrict const        finemeshptr,  /*+ Finer mesh                          +*/
Vmesh * restrict const              coarmeshptr,  /*+ Coarser mesh to build               +*/
Gnum * restrict * const             finecoarptr,  /*+ Pointer to multinode table to build +*/
const VmeshSeparateMlParam * const  paraptr)      /*+ Method parameters                   +*/
{
  int                 o;

  if (finemeshptr->m.vnodnbr <= (Gnum) paraptr->vnodnbr)
    return (1);

  if ((o = meshCoarsen (&finemeshptr->m, &coarmeshptr->m, finecoarptr, (Gnum) paraptr->vnodnbr, paraptr->coarrat, paraptr->coartype)) != 0)
    return (o);                                   /* Return if coarsening failed */

  coarmeshptr->parttax = NULL;                    /* Do not allocate partition data yet     */
  coarmeshptr->frontab = finemeshptr->frontab;    /* Re-use frontier array for coarser mesh */
  coarmeshptr->levlnum = finemeshptr->levlnum + 1; /* Mesh level is coarsening level        */

  return (0);
}

/* This routine propagates the separation of the
** coarser mesh back to the finer mesh, according
** to the multinode table of collapsed elements.
** After the separation is propagated, it finishes
** to compute the parameters of the finer mesh that
** were not computed at the coarsening stage.
** It returns:
** - 0   : if coarse mesh data has been propagated to fine mesh.
** - !0  : on error.
*/

static
int
vmeshSeparateMlUncoarsen (
Vmesh * restrict const        finemeshptr,        /*+ Finer mesh      +*/
const Vmesh * restrict const  coarmeshptr,        /*+ Coarser mesh    +*/
const Gnum * restrict const   finecoartax)        /*+ Multinode array +*/
{
  if (finemeshptr->parttax == NULL) {             /* If partition array not yet allocated */
    if ((finemeshptr->parttax = (GraphPart *) memAlloc ((finemeshptr->m.velmnbr + finemeshptr->m.vnodnbr) * sizeof (GraphPart))) == NULL) {
      errorPrint ("vmeshSeparateMlUncoarsen: out of memory");
      return     (1);                             /* Allocated data will be freed along with mesh structure */
    }
    finemeshptr->parttax -= finemeshptr->m.baseval;
  }

  if (coarmeshptr != NULL) {                      /* If coarser mesh provided */
    Gnum                finevelmnum;
    Gnum                fineecmpsize1;            /* Number of fine elements */
    Gnum                fineecmpsize2;
    Gnum                finevnodnum;
    Gnum                finencmpsize1;            /* Number of fine nodes                     */
    Gnum                finefronnbr;              /* Number of frontier vertices in fine mesh */

    for (finevelmnum = finemeshptr->m.velmbas, fineecmpsize1 = fineecmpsize2 = 0;
         finevelmnum < finemeshptr->m.velmnnd; finevelmnum ++) {
      Gnum                partval;

#ifdef SCOTCH_DEBUG_VMESH2
      if ((finecoartax[finevelmnum] < coarmeshptr->m.baseval) ||
          (finecoartax[finevelmnum] >= (coarmeshptr->m.velmnbr + coarmeshptr->m.vnodnbr + coarmeshptr->m.baseval))) {
        errorPrint ("vmeshSeparateMlUncoarsen: internal error (1)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH2 */
      partval = (Gnum) coarmeshptr->parttax[finecoartax[finevelmnum]];
      finemeshptr->parttax[finevelmnum] = partval;

      fineecmpsize1 += (partval & 1);             /* Superscalar update of counters */
      fineecmpsize2 += (partval & 2);
    }
    finemeshptr->ecmpsize[0] = finemeshptr->m.velmnbr - fineecmpsize1 - (fineecmpsize2 >> 1);
    finemeshptr->ecmpsize[1] = fineecmpsize1;

    for (finevnodnum = finemeshptr->m.vnodbas, finencmpsize1 = finefronnbr = 0;
         finevnodnum < finemeshptr->m.vnodnnd; finevnodnum ++) {
      Gnum                partval;

#ifdef SCOTCH_DEBUG_VMESH2
      if ((finecoartax[finevnodnum] <  coarmeshptr->m.vnodbas) || /* Sons of nodes are always nodes */
          (finecoartax[finevnodnum] >= coarmeshptr->m.vnodnnd)) {
        errorPrint ("vmeshSeparateMlUncoarsen: internal error (2)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_VMESH2 */
      partval = coarmeshptr->parttax[finecoartax[finevnodnum]];
      finemeshptr->parttax[finevnodnum] = partval;

      if ((partval & 2) != 0)                     /* If node is in separator   */
        finemeshptr->frontab[finefronnbr ++] = finevnodnum; /* Add to frontier */

      finencmpsize1 += (partval & 1);
    }

    finemeshptr->ncmpload[0] = coarmeshptr->ncmpload[0];
    finemeshptr->ncmpload[1] = coarmeshptr->ncmpload[1];
    finemeshptr->ncmpload[2] = coarmeshptr->ncmpload[2];
    finemeshptr->ncmploaddlt = coarmeshptr->ncmploaddlt;
    finemeshptr->ncmpsize[0] = finemeshptr->m.vnodnbr - finencmpsize1 - finefronnbr;
    finemeshptr->ncmpsize[1] = finencmpsize1;
    finemeshptr->fronnbr     = finefronnbr;
  }
  else                                            /* No coarse mesh provided       */
    vmeshZero (finemeshptr);                      /* Assign all vertices to part 0 */

#ifdef SCOTCH_DEBUG_VMESH2
  if (vmeshCheck (finemeshptr) != 0) {
    errorPrint ("vmeshSeparateMlUncoarsen: internal error (3)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_VMESH2 */

  return (0);
}

/* This routine recursively performs the
** coarsening recursion.
** It returns:
** - 0   : if recursion proceeded well.
** - !0  : on error.
*/

static
int
vmeshSeparateMl2 (
Vmesh * restrict const                      finemeshptr, /* Vertex-separation mesh */
const VmeshSeparateMlParam * restrict const paraptr) /* Method parameters          */
{
  Vmesh               coarmeshdat;
  Gnum * restrict     finecoartax;
  int                 o;

  o = 1;                                          /* Assume an error if "switch...case" returns a strange value in debug mode */
  switch (vmeshSeparateMlCoarsen (finemeshptr, &coarmeshdat, &finecoartax, paraptr)) {
    case 0 :
      if (((o = vmeshSeparateMl2         (&coarmeshdat, paraptr))                  == 0) &&
          ((o = vmeshSeparateMlUncoarsen (finemeshptr, &coarmeshdat, finecoartax)) == 0) &&
          ((o = vmeshSeparateSt          (finemeshptr, paraptr->stratasc))         != 0)) /* Apply ascending strategy */
        errorPrint ("vmeshSeparateMl2: cannot apply ascending strategy");
      coarmeshdat.frontab = NULL;                 /* Prevent frontab of fine mesh from being freed */
      vmeshExit (&coarmeshdat);
      memFree (finecoartax + finemeshptr->m.baseval); /* Free finecoartab as not part of coarse mesh vertex group (unlike for graphCoarsen) */
      break;
#ifdef SCOTCH_DEBUG_VMESH2
    case 1 :
    case 2 :                                      /* Cannot coarsen due to lack of memory */
      finecoartax = NULL;                         /* Prevent Valgrind from yelling */
#else /* SCOTCH_DEBUG_VMESH2 */
    default :
#endif /* SCOTCH_DEBUG_VMESH2 */
      if (((o = vmeshSeparateMlUncoarsen (finemeshptr, NULL, finecoartax)) == 0) && /* Finalize mesh    */
          ((o = vmeshSeparateSt          (finemeshptr, paraptr->stratlow)) != 0)) /* Apply low strategy */
        errorPrint ("vmeshSeparateMl2: cannot apply low strategy");
      break;
  }

  return (o);
}

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

int
vmeshSeparateMl (
Vmesh * restrict const              meshptr,      /*+ Node-separation mesh +*/
const VmeshSeparateMlParam * const  paraptr)      /*+ Method parameters    +*/
{
  Gnum                levlnum;                    /* Save value for mesh level */
  int                 o;

  levlnum = meshptr->levlnum;                     /* Save mesh level                */
  meshptr->levlnum = 0;                           /* Initialize coarsening level    */
  o = vmeshSeparateMl2 (meshptr, paraptr);        /* Perform multi-level separation */
  meshptr->levlnum = levlnum;                     /* Restore mesh level             */

  return (o);
}
