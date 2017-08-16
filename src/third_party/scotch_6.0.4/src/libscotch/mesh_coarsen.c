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
/**   NAME       : mesh_coarsen.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source mesh    **/
/**                coarsening functions.                   **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 30 jan 2004     **/
/**                                 to     05 may 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     12 sep 2007     **/
/**                                                        **/
/**   NOTES      : # The coarsening process is as follows. **/
/**                  First, node collapsing is performed,  **/
/**                  such that pairs of matching nodes are **/
/**                  created, or kept as single nodes.     **/
/**                  Then, elements are built, and merged  **/
/**                  whenever possible.                    **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH_COARSEN

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "mesh_coarsen.h"

/*
**  The static variables.
*/

static void              (* meshCoarFuncTab[MESHCOARSENNBR]) () = { /* Tables of matching routines */
                             meshCoarsenMatchNg };

/***************************/
/*                         */
/* The coarsening routine. */
/*                         */
/***************************/

/* This routine coarsens the given "finemesh" into
** "coarmesh", as long as the coarsening ratio remains
** below some threshold value and the coarsened mesh
** is not too small.
** It returns:
** - 0  : if the mesh has been coarsened.
** - 1  : if the mesh could not be coarsened.
** - 2  : on error.
*/

int
meshCoarsen (
const Mesh * restrict const   finemeshptr,        /*+ Mesh to coarsen                +*/
Mesh * restrict const         coarmeshptr,        /*+ Coarse mesh to build           +*/
Gnum * restrict * const       finecoarptr,        /*+ Pointer to multinode data      +*/
const Gnum                    coarnbr,            /*+ Minimum number of coarse nodes +*/
const double                  coarrat,            /*+ Maximum contraction ratio      +*/
const MeshCoarsenType         coartype)           /*+ Matching type                  +*/
{
  Gnum                        coarhashsiz;        /* Size of the hash table                      */
  Gnum                        coarhashmsk;        /* Mask for access to hash table               */
  MeshCoarsenHngb * restrict  coarhngbtab;        /* Table of edges to other multinodes          */
  MeshCoarsenHbdg * restrict  coarhbdgtab;        /* Table of bridge nodes to other multinodes   */
  Gnum * restrict             coarverttax;        /* Pointer to coarse vertex array              */
  Gnum * restrict             coarvelotax;        /* Pointer to coarse vertex load array         */
  Gnum * restrict             coaredgetax;        /* Pointer to coarse edge array                */
  Gnum                        coaredgenbr;        /* (Upper bound of) number of edges in mesh    */
  Gnum                        coaredgenum;        /* Number of current coarse edge               */
  Gnum                        coarvertnbr;        /* Number of vertices in coarse mesh           */
  Gnum                        coarvelmnbr;        /* Number of coarse element vertices           */
  Gnum                        coarvnodnbr;        /* Number of coarse node vertices              */
  Gnum                        finevertnbr;        /* Number of vertices in fine graph            */
  Gnum * restrict             finecoartax;        /* Based access to finecoartab                 */
  Gnum                        coarvelmnum;        /* Number of currently selected coarse element */
  Gnum                        coareelmnum;
  Gnum                        coarvnodnum;
  Gnum                        coardegrmax;
  MeshCoarsenMult * restrict  finemulttax;
  Gnum                        coaredgetmp;
  size_t                      coarvelooftval;
  size_t                      coaredgeoftval;

#ifdef SCOTCH_DEBUG_MESH2
  if (coartype >= MESHCOARSENNBR) {
    errorPrint ("meshCoarsen: invalid parameter");
    return     (2);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  memSet (coarmeshptr, 0, sizeof (Mesh));         /* Initialize coarse mesh */
  coarmeshptr->flagval = GRAPHFREEVERT;
  coarmeshptr->baseval = finemeshptr->baseval;

  finevertnbr = finemeshptr->velmnbr + finemeshptr->vnodnbr;
  if ((finecoartax = (Gnum *) memAlloc (finevertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("meshCoarsen: out of memory (1)"); /* Allocate coarse mesh uncoarsening array */
    return     (2);
  }
  memSet (finecoartax, ~0, finevertnbr * sizeof (Gnum));
  finecoartax -= finemeshptr->baseval;            /* Set based access to finecoartax */

  for (coarhashmsk = 31, coarhashsiz = finemeshptr->degrmax * finemeshptr->degrmax - 1; /* Compute size of hash table */
       coarhashmsk < coarhashsiz; coarhashmsk = coarhashmsk * 2 + 1) ;
  coarhashsiz = coarhashmsk + 1;

  if (memAllocGroup ((void **) (void *)
        &coarverttax, (size_t) ((finevertnbr + 1)     * sizeof (Gnum)), /* Upper bound on number of coarse vertices */
        &coarvelotax, (size_t) ( finevertnbr          * sizeof (Gnum)), /* Upper bound on number of coarse vertices */
        &coaredgetax, (size_t) ( finemeshptr->edgenbr * sizeof (Gnum)),
        &coarhngbtab, (size_t) ( coarhashsiz          * sizeof (MeshCoarsenHngb)),
        &coarhbdgtab, (size_t) ( coarhashsiz          * sizeof (MeshCoarsenHbdg)),
        &finemulttax, (size_t) ( finemeshptr->velmnbr * sizeof (MeshCoarsenMult)), NULL) == NULL) {
    errorPrint ("meshCoarsen: out of memory (2)"); /* Allocate coarser mesh structure */
    memFree    (finecoartax + finemeshptr->baseval);
    return     (2);
  }
  memSet (coarhngbtab, ~0, coarhashsiz * sizeof (MeshCoarsenHngb));
  memSet (coarhbdgtab, ~0, coarhashsiz * sizeof (MeshCoarsenHbdg));
  finemulttax -= coarmeshptr->baseval;

#define SCOTCH_DEBUG_MESH3
#ifdef SCOTCH_DEBUG_MESH3
fprintf (stderr, "-------- ENTERING COARSENING ---------\n");

fprintf (stderr, "finenodenbr=%ld, fineelemnbr=%ld, fineedgenbr=%ld, finedegrmax=%ld\n", (long) finemeshptr->vnodnbr, (long) finemeshptr->velmnbr, (long) finemeshptr->edgenbr, (long) finemeshptr->degrmax);
#endif /* SCOTCH_DEBUG_MESH3 */

  meshCoarFuncTab[coartype] (finemeshptr, finemulttax, finecoartax, &coarvelmnbr, &coarvnodnbr, &coaredgenbr); /* Call proper matching function */

#ifndef DEAD_CODE
  coarvnodnbr = finemeshptr->vnodnbr;             /* TODO : coarvnodnbr estimator is wrong : too tight */
  coaredgenbr = finemeshptr->edgenbr;
#endif /* DEAD_CODE */
  coarvertnbr = coarvelmnbr + coarvnodnbr;

  memOffset ((void *) coarverttax,
             &coarverttax, (size_t) ((coarvertnbr + 1) * sizeof (Gnum)),
             &coarvelotax, (size_t) ( coarvertnbr      * sizeof (Gnum)),
             &coaredgetax, (size_t) ( coaredgenbr      * sizeof (Gnum)), NULL); /* Hast tables and finemulttax stay in place */
  coarverttax -= coarmeshptr->baseval;
  coarvelotax -= coarmeshptr->baseval;
  coaredgetax -= coarmeshptr->baseval;

  coarmeshptr->velmbas  = coarmeshptr->baseval;
  coarmeshptr->velmnbr  = coarvelmnbr;
  coarmeshptr->velmnnd  =
  coarmeshptr->vnodbas  = coarvelmnbr + coarmeshptr->velmbas;

  for (coarvelmnum = coaredgenum = coarmeshptr->baseval, coarvnodnum = coarmeshptr->vnodbas, coardegrmax = 0; /* For all coarse elements */
       coarvelmnum < coarmeshptr->velmnnd; coarvelmnum ++) {
    Gnum                coarveloval;              /* Weight of coarsened element    */
    Gnum                coarvnisnum;              /* Number of coarse isolated node */
    Gnum                finevelmnum;              /* Number of current element      */
    int                 i;

    coarverttax[coarvelmnum] = coaredgenum;

    coarvnisnum = ~0;                             /* No isolated node yet for this element pair */
    coarveloval = 0;
    i = 0;
    do {                                          /* For both elements of element pair (if they are different) */
      Gnum                fineeelmnum;

      finevelmnum = finemulttax[coarvelmnum].finevelmnum[i]; /* Get number of current element */
      coarveloval += ((finemeshptr->velotax != NULL) ? finemeshptr->velotax[finevelmnum] : 1);
      for (fineeelmnum = finemeshptr->verttax[finevelmnum];
           fineeelmnum < finemeshptr->vendtax[finevelmnum]; fineeelmnum ++) {
        Gnum                finevnodnum;          /* Number of current node neighbor */
        Gnum                fineenodnum;
        Gnum                finevdegval;
        Gnum                finevnloval;
        Gnum                finevelmend;
        Gnum                coarvnodtmp;
        Gnum                coarhnodtmp;

        finevnodnum = finemeshptr->edgetax[fineeelmnum];
        fineenodnum = finemeshptr->verttax[finevnodnum];
        finevdegval = finemeshptr->vendtax[finevnodnum] - fineenodnum;
        finevnloval = (finemeshptr->vnlotax != NULL) ? finemeshptr->vnlotax[finevnodnum] : 1;

        if ((finevdegval == 2) &&                 /* If node is an external bridge to another coarse element */
            ((finevelmend = (finemeshptr->edgetax[fineenodnum] + finemeshptr->edgetax[fineenodnum + 1] - finevelmnum)) != finemulttax[coarvelmnum].finevelmnum[1 - i])) {
          Gnum                coarvelmend;
          Gnum                coarhelmend;

          coarvelmend = finecoartax[finevelmend]; /* Get coarse index of end element */
          coarvnodtmp = finecoartax[finevnodnum]; /* Get coarse number of fine node  */
          for (coarhelmend = (coarvelmend * MESHCOARSENHASHPRIME) & coarhashmsk; ; coarhelmend = (coarhelmend + 1) & coarhashmsk) {
            if (coarhbdgtab[coarhelmend].coarvelmnum != coarvelmnum) { /* If bridge not yet considered */
              coarhbdgtab[coarhelmend].coarvelmnum = coarvelmnum; /* Add it to element neighbor list   */
              coarhbdgtab[coarhelmend].coarvelmend = coarvelmend;
              if (coarvnodtmp == -1) {            /* If bridge nodes not considered before by other element */
                coarhbdgtab[coarhelmend].coarvnodnum   = /* Assign it                                       */
                finecoartax[finevnodnum] = coarvnodtmp = coarvnodnum ++;
                coarverttax[coarvnodtmp] = 2;     /* Prepare the fact that another element will see the node */
                coarvelotax[coarvnodtmp] = finevnloval;
              }
              coaredgetax[coaredgenum ++] = coarvnodtmp; /* Directly add coarse node to element neighborhood */
              break;
            }
            if (coarhbdgtab[coarhelmend].coarvelmend == coarvelmend) { /* If bridge already present               */
              if (coarvnodtmp == -1) {            /* If we are the first element to see the bridge node           */
                finecoartax[finevnodnum]  = coarvnodtmp = coarhbdgtab[coarhelmend].coarvnodnum; /* Assign it      */
                coarvelotax[coarvnodtmp] += finevnloval; /* Update the weight of the node                         */
              }                                   /* Else node already processed with full load, so nothing to do */
              break;
            }
          }
          continue;                               /* Edge has been added or will not be */
        }
        else if (finevdegval < 3) {               /* Else if node is isolated or is an internal bridge */
          if ((finevdegval == 2) &&               /* Process bridge edges only once                    */
              (finevelmnum >= finemulttax[coarvelmnum].finevelmnum[1 - i]))
            continue;

          if (coarvnisnum == ~0) {                /* If no isolated node for this element pair */
            coarvnisnum = coarvnodnum ++;         /* Create isolated node                      */
            coarverttax[coarvnisnum]    = 1;
            coarvelotax[coarvnisnum]    = finevnloval;
            coaredgetax[coaredgenum ++] = coarvnisnum;
          }
          else                                    /* If isolated node already exists */
            coarvelotax[coarvnisnum] += finevnloval; /* Add node contribution to it  */
          finecoartax[finevnodnum] = coarvnisnum; /* Map fine node to isolated node  */
          continue;
        }
        else {
          coarvnodtmp = finecoartax[finevnodnum]; /* Get coarse number of fine node    */
          if (coarvnodtmp == ~0) {                /* If coarse number not yet assigned */
            finecoartax[finevnodnum] = coarvnodtmp = coarvnodnum ++; /* Assign it      */
            coarverttax[coarvnodtmp] = 0;         /* No connections to the node yet    */
            coarvelotax[coarvnodtmp] = finevnloval;
          }
        }

        for (coarhnodtmp = (coarvnodtmp * MESHCOARSENHASHPRIME) & coarhashmsk; ; coarhnodtmp = (coarhnodtmp + 1) & coarhashmsk) {
          if (coarhngbtab[coarhnodtmp].coarvelmnum != coarvelmnum) { /* If node neighbor not yet considered */
            coarhngbtab[coarhnodtmp].coarvelmnum = coarvelmnum; /* Add it to element neighbor list          */
            coarhngbtab[coarhnodtmp].coarvnodnum = coarvnodtmp;
            coaredgetax[coaredgenum ++]          = coarvnodtmp;
            coarverttax[coarvnodtmp] ++;          /* One more edge referencing the node */
            break;
          }
          if (coarhngbtab[coarhnodtmp].coarvnodnum == coarvnodtmp) /* If node already present, nothing to do */
            break;
        }
      }
    } while (i ++, finevelmnum != finemulttax[coarvelmnum].finevelmnum[1]);

    coarvelotax[coarvelmnum] = coarveloval;       /* Lose initial weights of elements, if any, to keep coarsening weights */
    if ((coaredgenum - coarverttax[coarvelmnum]) > coardegrmax)
      coardegrmax = (coaredgenum - coarverttax[coarvelmnum]);
  }
  coarmeshptr->vnodnnd = coarvnodnum;
  coarmeshptr->vnodnbr = coarvnodnum - coarmeshptr->vnodbas;
  coarmeshptr->velosum = finemeshptr->velosum;
  coarmeshptr->vnlosum = finemeshptr->vnlosum;
  coarmeshptr->edgenbr = 2 * (coaredgenum - coarmeshptr->baseval);

  for (coarvnodnum = coarmeshptr->vnodbas, coaredgetmp = coaredgenum; /* Build start indices for node edge sub-arrays */
       coarvnodnum < coarmeshptr->vnodnnd; coarvnodnum ++) {
    Gnum                coardegrval;

    coardegrval = coarverttax[coarvnodnum];
    coarverttax[coarvnodnum] = coaredgetmp;
    coaredgetmp += coardegrval;

    if (coardegrval > coardegrmax)
      coardegrmax = coardegrval;
  }
  coarmeshptr->degrmax = coardegrmax;
  for (coarvelmnum = coareelmnum = coarmeshptr->baseval;
       coarvelmnum < coarmeshptr->velmnnd; coarvelmnum ++) {
    Gnum                coareelmnnd;

    coareelmnnd = (coarvelmnum < (coarmeshptr->velmnnd - 1)) ? coarverttax[coarvelmnum + 1] : coaredgenum;
    while (coareelmnum < coareelmnnd) {
      Gnum                coarvnodnum;

      coarvnodnum = coaredgetax[coareelmnum ++];
      coaredgetax[coarverttax[coarvnodnum] ++] = coarvelmnum;
    }
  }
  memMov (&coarverttax[coarmeshptr->vnodbas + 1], /* Re-build start indices for node edge sub-arrays */
          &coarverttax[coarmeshptr->vnodbas],
          coarmeshptr->vnodnbr * sizeof (Gnum));
  coarverttax[coarmeshptr->vnodbas] = coaredgenum;

  coarvelooftval = coarvelotax - coarverttax;
  coaredgeoftval = coaredgetax - coarverttax;
  coarverttax = memRealloc (coarverttax + coarmeshptr->baseval, (coaredgeoftval + coarmeshptr->edgenbr) * sizeof (Gnum)); /* Re-allocate array to save space */
  coarmeshptr->verttax = coarverttax - coarmeshptr->baseval;
  coarmeshptr->vendtax = coarmeshptr->verttax + 1;
  coarmeshptr->velotax = coarmeshptr->verttax + coarvelooftval;
  coarmeshptr->vnlotax = coarmeshptr->velotax;    /* Same array for both vertex load sub-arrays */
  coarmeshptr->edgetax = coarmeshptr->verttax + coaredgeoftval;

#ifdef SCOTCH_DEBUG_MESH2
  if (meshCheck (coarmeshptr) != 0) {             /* Check mesh consistency */
    errorPrint ("meshCoarsen: internal error (7)");
    return     (2);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  *finecoarptr = finecoartax;                     /* Return multinode array */

#ifdef SCOTCH_DEBUG_MESH3
fprintf (stderr, "coarvnodnbr=%ld\tcoarvelmnbr=%ld\tcoaredgenbr=%ld, coardegrmax=%ld\n", (long) coarmeshptr->vnodnbr, (long) coarmeshptr->velmnbr, (long) coarmeshptr->edgenbr, (long) coarmeshptr->degrmax);
fprintf (stderr, "-------- EXITING COARSENING ---------\n"); /* TODO REMOVE */
#endif /* SCOTCH_DEBUG_MESH3 */

  return (0);
}

/********************************************/
/*                                          */
/* The matching subroutines. In fact, these */
/* are merging routines, which merge        */
/* elements of the fine mesh to form larger */
/* elements in the coarse mesh.             */
/* New elements are ordered in increasing   */
/* order from baseval, while nodes are      */
/* ordered in decreasing order from -2, as  */
/* -1 is a reserved flag value used         */
/* for labelling non yet considered         */
/* vertices.                                */
/*                                          */
/********************************************/


/* This routine performs elements matching by
** selecting the elements that share most nodes
** with the first element.
*/

static
void
meshCoarsenMatchNg (
const Mesh * restrict const       finemeshptr,    /* Fine mesh to perform matching on                           */
MeshCoarsenMult * restrict const  finemulttax,    /* Array of fine multielements                                */
Gnum * restrict const             finecoartax,    /* Fine to coarse vertex array                                */
Gnum * restrict const             coarvelmptr,    /* Pointer to number of coarse element vertices               */
Gnum * restrict const             coarvnodptr,    /* Pointer to (upper bound on) number of coarse node vertices */
Gnum * restrict const             coaredgeptr)    /* Pointer to (upper bound on) number of edges                */
{
  Gnum                          coarvelmnum;      /* Number of current coarse element vertex */
  Gnum                          finepertbas;      /* Index of base of perturbation area      */
  Gnum                          finepertnbr;      /* Size of perturbation area               */
  MeshCoarsenNgHash * restrict  finehashtab;      /* Hash table of neighbor elements         */
  Gnum                          finehashsiz;
  Gnum                          finehashmsk;
  Gnum                          coarvnodnbr;
  Gnum                          coaredgenbr;

  for (finehashmsk = 31, finehashsiz = finemeshptr->degrmax * finemeshptr->degrmax - 1; /* Compute size of hash table */
       finehashmsk < finehashsiz; finehashmsk = finehashmsk * 2 + 1) ;
  finehashsiz = finehashmsk + 1;
  if ((finehashtab = (MeshCoarsenNgHash *) memAlloc (finehashsiz * sizeof (MeshCoarsenNgHash))) == NULL) {
    *coarvelmptr = finemeshptr->velmnbr;          /* Indicate no coarsening occured */
    return;
  }
  memSet (finehashtab, ~0, finehashsiz * sizeof (MeshCoarsenNgHash));
  finehashmsk = finehashsiz - 1;

  coarvelmnum = finemeshptr->baseval;             /* Start numbering elements in ascending order */
  coarvnodnbr = finemeshptr->vnodnbr;
  coaredgenbr = finemeshptr->edgenbr;

  if (finemeshptr->velotax != NULL) {             /* If fine mesh has element coarsening vertex weights, perform first pass */
    Gnum                          finevelomin;
    Gnum                          finevelomax;
    Gnum                          finevelmnum;

    finevelomin = (3 * finemeshptr->velosum) / (5 * finemeshptr->velmnbr);
    finevelomax = (5 * finemeshptr->velosum) /      finemeshptr->velmnbr;

    for (finevelmnum = finemeshptr->velmbas; finevelmnum < finemeshptr->velmnnd; finevelmnum ++) {
      Gnum                fineeelmnum;
      Gnum                finehelmnum;
      Gnum                finevnisnbr;            /* Number of isolated node vertices         */
      Gnum                finehebsnum;            /* Hash number of best matching element     */
      Gnum                finevebsnum;            /* Number of best matching element          */
      Gnum                finevnbsnbr;            /* Number of nodes shared with best element */

      if (finecoartax[finevelmnum] != ~0)         /* If element already selected */
        continue;
      if (finemeshptr->velotax[finevelmnum] >= finevelomin) { /* If element is large enough, leave it for the second pass */
        if (finemeshptr->velotax[finevelmnum] > finevelomax) { /* Except if it is too large, as then it is not matched    */
          finecoartax[finevelmnum] = coarvelmnum;
          finemulttax[coarvelmnum].finevelmnum[0] =
          finemulttax[coarvelmnum].finevelmnum[1] = finevelmnum;
fprintf (stderr, "++ %ld %ld\n", (long) finevelmnum, (long) finemeshptr->velotax[finevelmnum]); /* TODO REMOVE */
          coarvelmnum ++;                         /* One more single vertex created */
        }
        continue;
      }

      finecoartax[finevelmnum] = coarvelmnum;     /* Set vertex as used so that it will not be considered as an end vertex */

      finehelmnum = (finevelmnum * MESHCOARSENHASHPRIME) & finehashmsk;
      finehashtab[finehelmnum].velmnum = finevelmnum; /* Put element in hash table so that number of end vertex is right even for uncoarsened elements */
      finehashtab[finehelmnum].velmend = finevelmnum;
      finehebsnum = finehelmnum;                  /* Mate is element itself */
      finevnbsnbr = 0;                            /* Will never be selected */

      finevnisnbr = 0;                            /* No isolated node vertices yet                       */
      for (fineeelmnum = finemeshptr->verttax[finevelmnum]; /* For all node neighbors of current element */
           fineeelmnum < finemeshptr->vendtax[finevelmnum]; fineeelmnum ++) {
        Gnum                finevnodnum;
        Gnum                fineenodnum;
        Gnum                fineenodnnd;
        Gnum                finevdegval;
        Gnum                finevnbgval;

        finevnodnum = finemeshptr->edgetax[fineeelmnum];
        fineenodnum = finemeshptr->verttax[finevnodnum];
        fineenodnnd = finemeshptr->vendtax[finevnodnum];
        finevdegval = fineenodnnd - fineenodnum;
        if (finevdegval == 1) {                   /* If node is isolated */
          finevnisnbr ++;
          continue;                               /* Directly skip to next node */
        }
        finevnbgval = (finevdegval == 2) ? 1 : 0; /* If node is a bridge which connects the element to only one other element */

        for ( ; fineenodnum < fineenodnnd; fineenodnum ++) { /* For all elements which are neighbors of current node */
          Gnum                finevelmend;
          Gnum                finehelmend;
          Gnum                finevnngnbr;        /* Current number of neigoboring nodes that connect the two elements */

          finevelmend = finemeshptr->edgetax[fineenodnum];
          if (finecoartax[finevelmend] != ~0)     /* If end element vertex already matched, do not consider it */
            continue;

          for (finehelmend = (finevelmend * MESHCOARSENHASHPRIME) & finehashmsk; ; finehelmend = (finehelmend + 1) & finehashmsk) {
            if (finehashtab[finehelmend].velmnum != finevelmnum) { /* If element neighbor not yet considered */
              finevnngnbr = 1;
              finehashtab[finehelmend].velmnum = finevelmnum;
              finehashtab[finehelmend].velmend = finevelmend;
              finehashtab[finehelmend].vnngnbr = finevnngnbr;
              finehashtab[finehelmend].vnbgnbr = finevnbgval;
            }
            else if (finehashtab[finehelmend].velmend == finevelmend) { /* Else if element found */
              finevnngnbr = ++ finehashtab[finehelmend].vnngnbr;
              finehashtab[finehelmend].vnbgnbr += finevnbgval;
            }
            else                                  /* Else go on searching */
              continue;

            if (finevnngnbr > finevnbsnbr) {
              finehebsnum = finehelmend;
              finevnbsnbr = finevnngnbr;
            }
            break;
          }
        }
      }

      finevebsnum = finehashtab[finehebsnum].velmend;
      finemulttax[coarvelmnum].finevelmnum[0] = finevelmnum; /* Set matching pair */
      finemulttax[coarvelmnum].finevelmnum[1] = finevebsnum;
      if (finevelmnum != finevebsnum) {           /* If a matching element has been found */
        finecoartax[finevebsnum] = coarvelmnum;
        if (finevnisnbr > 0)
          finevnisnbr --;
        coarvnodnbr -= finehashtab[finehebsnum].vnbgnbr + finevnisnbr;
        coaredgenbr -= 2 * finevnisnbr + 4 * finehashtab[finehebsnum].vnbgnbr;
      }
      coarvelmnum ++;                             /* Number nodes in ascending order */
    }
  }

  for (finepertbas = finemeshptr->velmbas,        /* Run cache-friendly perturbation on elements  */
       finepertnbr = 2 + intRandVal (MESHCOARSENPERTPRIME - 2); /* Compute perturbation area size */
       finepertbas < finemeshptr->velmnnd; finepertbas += finepertnbr) {
    Gnum                finepertval;              /* Current index in perturbation area */

    if (finepertbas + finepertnbr > finemeshptr->velmnnd)
      finepertnbr = finemeshptr->velmnnd - finepertbas;

    finepertval = 0;                              /* Start from first perturbation element vertex     */
    do {                                          /* Loop on perturbation element vertices            */
      Gnum                finevelmnum;            /* Number of currently selected fine element vertex */
      Gnum                fineeelmnum;
      Gnum                finehelmnum;
      Gnum                finevnisnbr;            /* Number of isolated node vertices         */
      Gnum                finehebsnum;            /* Hash number of best matching element     */
      Gnum                finevebsnum;            /* Number of best matching element          */
      Gnum                finevnbsnbr;            /* Number of nodes shared with best element */

      finevelmnum = finepertbas + finepertval;    /* Compute corresponding elemennt number */
      if (finecoartax[finevelmnum] != ~0)         /* If element already selected           */
        continue;

      finecoartax[finevelmnum] = coarvelmnum;     /* Set vertex as used so that it will not be considered as an end vertex */

      finehelmnum = (finevelmnum * MESHCOARSENHASHPRIME) & finehashmsk;
      finehashtab[finehelmnum].velmnum = finevelmnum; /* Put element in hash table so that number of end vertex is right even for uncoarsened elements */
      finehashtab[finehelmnum].velmend = finevelmnum;
      finehebsnum = finehelmnum;                  /* Mate is element itself */
      finevnbsnbr = 0;                            /* Will never be selected */

      finevnisnbr = 0;                            /* No isolated node vertices yet                         */
      for (fineeelmnum = finemeshptr->verttax[finevelmnum]; /* For all node neighbors of current element */
           fineeelmnum < finemeshptr->vendtax[finevelmnum]; fineeelmnum ++) {
        Gnum                finevnodnum;
        Gnum                fineenodnum;
        Gnum                fineenodnnd;
        Gnum                finevdegval;
        Gnum                finevnbgval;

        finevnodnum = finemeshptr->edgetax[fineeelmnum];
        fineenodnum = finemeshptr->verttax[finevnodnum];
        fineenodnnd = finemeshptr->vendtax[finevnodnum];
        finevdegval = fineenodnnd - fineenodnum;
        if (finevdegval == 1) {                   /* If node is isolated */
          finevnisnbr ++;
          continue;                               /* Directly skip to next node */
        }
        finevnbgval = (finevdegval == 2) ? 1 : 0; /* If node is a bridge which connects the element to only one other element */

        for ( ; fineenodnum < fineenodnnd; fineenodnum ++) { /* For all elements which are neighbors of current node */
          Gnum                finevelmend;
          Gnum                finehelmend;
          Gnum                finevnngnbr;        /* Current number of neigoboring nodes that connect the two elements */

          finevelmend = finemeshptr->edgetax[fineenodnum];
          if (finecoartax[finevelmend] != ~0)     /* If end element vertex already matched, do not consider it */
            continue;

          for (finehelmend = (finevelmend * MESHCOARSENHASHPRIME) & finehashmsk; ; finehelmend = (finehelmend + 1) & finehashmsk) {
            if (finehashtab[finehelmend].velmnum != finevelmnum) { /* If element neighbor not yet considered */
              finevnngnbr = 1;
              finehashtab[finehelmend].velmnum = finevelmnum;
              finehashtab[finehelmend].velmend = finevelmend;
              finehashtab[finehelmend].vnngnbr = finevnngnbr;
              finehashtab[finehelmend].vnbgnbr = finevnbgval;
            }
            else if (finehashtab[finehelmend].velmend == finevelmend) { /* Else if element found */
              finevnngnbr = ++ finehashtab[finehelmend].vnngnbr;
              finehashtab[finehelmend].vnbgnbr += finevnbgval;
            }
            else                                  /* Else go on searching */
              continue;

            if (finevnngnbr > finevnbsnbr) {
              finehebsnum = finehelmend;
              finevnbsnbr = finevnngnbr;
            }
            break;
          }
        }
      }

      finevebsnum = finehashtab[finehebsnum].velmend;
      finemulttax[coarvelmnum].finevelmnum[0] = finevelmnum; /* Set matching pair */
      finemulttax[coarvelmnum].finevelmnum[1] = finevebsnum;
      if (finevelmnum != finevebsnum) {           /* If a matching element has been found */
        finecoartax[finevebsnum] = coarvelmnum;
        if (finevnisnbr > 0)
          finevnisnbr --;
        coarvnodnbr -= finehashtab[finehebsnum].vnbgnbr + finevnisnbr;
        coaredgenbr -= 2 * finevnisnbr + 4 * finehashtab[finehebsnum].vnbgnbr;
      }
      coarvelmnum ++;                             /* Number nodes in ascending order */
    } while ((finepertval = (finepertval + MESHCOARSENPERTPRIME) % finepertnbr) != 0); /* Compute next perturbation index */
  }

  memFree (finehashtab);

  *coarvelmptr = coarvelmnum - finemeshptr->velmbas;
  *coarvnodptr = coarvnodnbr;
  *coaredgeptr = coaredgenbr;

  return;
}
