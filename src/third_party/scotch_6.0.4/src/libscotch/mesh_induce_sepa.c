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
/**   NAME       : mesh_induce_sepa.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the source         **/
/**                mesh separator subgraph-making          **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 26 jan 2003     **/
/**                                 to     11 may 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     03 apr 2008     **/
/**                                                        **/
/**   NOTES      : # This routine differs from the         **/
/**                  standard mesh induction routine by    **/
/**                  the fact that neighboring elements do **/
/**                  not bear the same part number as the  **/
/**                  vertices to keep and are found on the **/
/**                  fly.                                  **/
/**                                                        **/
/**                # This routine is used in hmeshInduce-  **/
/**                  Nd to build the induced separator     **/
/**                  mesh.                                 **/
/**                  Since separator vertices are likely   **/
/**                  to be surrounded by  elements that    **/
/**                  belong to different parts but connect **/
/**                  the same vertices in the separator,   **/
/**                  a simple detection mechanism is used  **/
/**                  to try to reduce the number of such   **/
/**                  redundant elements.                   **/
/**                  Also, since this mesh is to be        **/
/**                  ordered with the highest available    **/
/**                  numbers, halo is not really needed    **/
/**                  so this "mesh*" routine is used. If   **/
/**                  halo was needed, a "hmesh*" routine   **/
/**                  should be used instead.               **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH
#define MESH_INDUCE_SEPA

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "mesh_induce_sepa.h"

/***************************************/
/*                                     */
/* This routine handles source meshes. */
/*                                     */
/***************************************/

/* This routine builds the separator mesh
** induced by the original mesh, the list of
** selected separator node vertices, and the
** vertex part array. Elements which are
** adjacent to the selected nodes are themselves
** selected.
** The induced vnumtab array is the baseval-based
** list of remaining node vertices if the original
** mesh does not have a vnumtab, or the proper
** subset of the original vnumtab else.
** This routine builds a mesh with nodes first
** and elements last. If a halo version of this
** algorithm had to be created, it should be
** re-worked such that elements are put first
** and nodes last, to enforce the structure of
** halo meshes.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
meshInduceSepa (
const Mesh * restrict const       orgmeshptr,     /* Pointer to original mesh             */
const GraphPart * restrict const  orgparttax,     /* Array of vertex partition flags      */
const Gnum                        orgsepanbr,     /* Number of node vertices in separator */
const Gnum * restrict const       orgsepatab,     /* Array of node indices                */
Mesh * restrict const             indmeshptr)     /* Pointer to induced submesh           */
{
  Gnum                          orgvertnbr;       /* Number of vertices in original mesh                        */
  Gnum                          orgbitssiz;       /* Size of int array holding as many bits as element vertices */
  unsigned int * restrict       orgbitstab;       /* Array of bit flags for mesh elements                       */
  Gnum                          indvnumnum;       /* Current vertex number index in separator node list         */
  Gnum * restrict               indvneltax;       /* Array of original numbers for induced elements             */
  Gnum                          indvelmnbr;       /* (Approximate) number of element vertices in induced mesh   */
  Gnum                          indvertnbr;       /* (Approximate) number of vertices in induced mesh           */
  Gnum                          indvnodnum;       /* Number of current node vertex in induced mesh              */
  Gnum                          indvnlosum;       /* Sum of node vertex load array                              */
  Gnum                          indvelmnum;       /* Number of current element vertex in induced mesh           */
  Gnum * restrict               indedgetax;       /* Based array to edge array of induced mesh                  */
  Gnum                          indedgenum;       /* Number of current edge to be created in induced mesh       */
  Gnum                          indedgenbr;       /* (Approximate) number of edges in induced mesh graph        */
  Gnum * restrict               orgindxtax;       /* Original to induced vertex number translation              */
  MeshInduceSepaHash * restrict hashtab;          /* Node hash array                                            */
  Gnum                          hashmsk;          /* Mask for hash array                                        */
  void * restrict               p;

  orgbitssiz = (orgmeshptr->velmnbr + (INTSIZEBITS - 1)) / INTSIZEBITS; /* Compute size of element bit array */
  if ((orgbitstab = (unsigned int *) memAlloc (orgbitssiz * sizeof (unsigned int))) == NULL) {
    errorPrint ("meshInduceSepa: out of memory (1)");
    return     (1);
  }
  memSet (orgbitstab, ~0, orgbitssiz * sizeof (unsigned int));

  for (indvnumnum = 0, indvelmnbr = 0, indedgenbr = 0; /* For all separator nodes in list */
       indvnumnum < orgsepanbr; indvnumnum ++) {
    Gnum                orgvnodnum;               /* Number of current node in original mesh */
    Gnum                orgenodnum;

    orgvnodnum = orgsepatab[indvnumnum];          /* Get number of separator node */

    for (orgenodnum = orgmeshptr->verttax[orgvnodnum], indedgenbr += (orgmeshptr->vendtax[orgvnodnum] - orgenodnum);
         orgenodnum < orgmeshptr->vendtax[orgvnodnum]; orgenodnum ++) {
      Gnum                orgvelmidx;             /* Index of incident element in bit array */
      Gnum                orgvelmbit;             /* Bit to update in element bit array     */

      orgvelmidx  = orgmeshptr->edgetax[orgenodnum] - orgmeshptr->velmbas;
      orgvelmbit  = orgvelmidx & (INTSIZEBITS - 1);
      orgvelmidx /= INTSIZEBITS;
      indvelmnbr += (orgbitstab[orgvelmidx] >> orgvelmbit) & 1; /* Count element vertex if not yet accounted for */
      orgbitstab[orgvelmidx] &= ~(1 << orgvelmbit);
    }
  }
  indedgenbr *= 2;                                /* Number of arcs is twice number of edges */
  memFree (orgbitstab);

  memSet (indmeshptr, 0, sizeof (Mesh));          /* Initialize mesh fields */
  indmeshptr->baseval = orgmeshptr->baseval;
  indmeshptr->flagval = MESHFREETABS | MESHVERTGROUP;
  indmeshptr->vnodnbr = orgsepanbr;               /* Nodes first */
  indmeshptr->vnodbas = indmeshptr->baseval;
  indmeshptr->vnodnnd =
  indmeshptr->velmbas = orgsepanbr + indmeshptr->baseval; /* Elements last */
  indmeshptr->veisnbr = 0;                        /* No isolated elements  */

  for (hashmsk = 15; hashmsk < orgmeshptr->degrmax; hashmsk = hashmsk * 2 + 1) ;
  hashmsk = hashmsk * 4 + 3;                      /* Compute size of node hash table */

  indvertnbr = indvelmnbr + orgsepanbr;           /* Upper bound on number of all vertices */
  if (orgmeshptr->velotax != NULL) {
    p = memAllocGroup ((void **) (void *)
          &indmeshptr->verttax, (size_t) (indvertnbr * sizeof (Gnum)), /* verttab is not compact */
          &indmeshptr->vendtax, (size_t) (indvertnbr * sizeof (Gnum)),
          &indmeshptr->vnumtax, (size_t) (orgsepanbr * sizeof (Gnum)), /* vnumtab is of size indnodenbr                */
          &indmeshptr->velotax, (size_t) (orgsepanbr * sizeof (Gnum)), NULL); /* Element vertices assumed not weighted */
    indmeshptr->velotax -= indmeshptr->vnodbas;
  }
  else
    p = memAllocGroup ((void **) (void *)
          &indmeshptr->verttax, (size_t) (indvertnbr * sizeof (Gnum)),
          &indmeshptr->vendtax, (size_t) (indvertnbr * sizeof (Gnum)),
          &indmeshptr->vnumtax, (size_t) (orgsepanbr * sizeof (Gnum)), NULL); /* vnumtab is of size indnodenbr */
  orgvertnbr = orgmeshptr->velmnbr + orgmeshptr->vnodnbr;
  if (p != 0) {
    indmeshptr->verttax -= indmeshptr->baseval;
    indmeshptr->vendtax -= indmeshptr->baseval;
    indmeshptr->vnumtax -= indmeshptr->vnodbas; /* vnumtab is of size indmeshptr->vnodnbr */

    p = memAllocGroup ((void **) (void *)
          &indedgetax, (size_t) (indedgenbr    * sizeof (Gnum)),
          &indvneltax, (size_t) (indvelmnbr    * sizeof (Gnum)),
          &orgindxtax, (size_t) (orgvertnbr    * sizeof (Gnum)),
          &hashtab,    (size_t) ((hashmsk + 1) * sizeof (MeshInduceSepaHash)), NULL);
  }
  if (p == 0) {
    errorPrint ("meshInduceSepa: out of memory (2)"); /* Allocate induced mesh graph structure */
    return     (1);
  }
  memSet (orgindxtax, ~0, orgvertnbr    * sizeof (Gnum));
  memSet (hashtab,    ~0, (hashmsk + 1) * sizeof (MeshInduceSepaHash));
  indedgetax -= indmeshptr->baseval;
  indvneltax -= indmeshptr->velmbas;
  orgindxtax -= orgmeshptr->baseval;

  for (indvnumnum = 0, indvnodnum = indedgenum = indmeshptr->baseval, indvelmnum = orgsepanbr + indvnodnum, indedgenbr = 0, indvnlosum = 0;
       indvnumnum < orgsepanbr; indvnumnum ++) {  /* For all node vertices in separator      */
    Gnum                orgvnodnum;               /* Number of current node in original mesh */
    Gnum                orgenodnum;
    Gnum                indenodnum;               /* Current index in node edge sub-array */

    orgvnodnum = orgsepatab[indvnumnum];          /* Get number of separator node */

    if (orgindxtax[orgvnodnum] == ~0)             /* If separator node not yet numbered          */
      orgindxtax[orgvnodnum] = indvnodnum ++;     /* One more node vertex created                */
    indmeshptr->verttax[orgindxtax[orgvnodnum]] = indedgenum; /* Set beginning of edge sub-array */
    indmeshptr->vnumtax[orgindxtax[orgvnodnum]] = orgvnodnum - (orgmeshptr->vnodbas - orgmeshptr->baseval);
    if (indmeshptr->velotax != NULL) {
      Gnum                indvnloval;

      indvnloval  = orgmeshptr->velotax[orgvnodnum];
      indvnlosum += indvnloval;
      indmeshptr->velotax[orgindxtax[orgvnodnum]] = indvnloval;
    }

    indenodnum  = indedgenum;                     /* Record position of node edge sub-array */
    indedgenum += orgmeshptr->vendtax[orgvnodnum] - /* Reserve space for node edges         */
                  orgmeshptr->verttax[orgvnodnum];

    for (orgenodnum = orgmeshptr->verttax[orgvnodnum]; /* For all element vertices neighboring the separator node */
         orgenodnum < orgmeshptr->vendtax[orgvnodnum]; orgenodnum ++) {
      Gnum                orgvelmnum;             /* Number of element node in original mesh */

      orgvelmnum = orgmeshptr->edgetax[orgenodnum]; /* Get number of element */

      if (orgindxtax[orgvelmnum] == -2)           /* If discarded element       */
        continue;                                 /* Skip to next element       */
      if (orgindxtax[orgvelmnum] == ~0) {         /* If element not yet created */
        Gnum                indedgetmp;           /* Save value for edge index  */
        Gnum                orgeelmnum;
        Gnum                indenodtmp;

        indmeshptr->verttax[indvelmnum] = indedgenum; /* Set beginning of element edge sub-array */
        indedgetmp = indedgenum;

        for (orgeelmnum = orgmeshptr->verttax[orgvelmnum]; /* For all nodes neighboring the element vertex */
             orgeelmnum < orgmeshptr->vendtax[orgvelmnum]; orgeelmnum ++) {
          Gnum                orgvnodend;         /* Number of current end node vertex in original mesh */

          orgvnodend = orgmeshptr->edgetax[orgeelmnum]; /* Get number of end node */

          if (orgparttax[orgvnodend] == 2) {      /* If end node is in separator */
            Gnum                hashnum;

            if (orgindxtax[orgvnodend] == ~0)     /* If end node not yet numbered  */
              orgindxtax[orgvnodend] = indvnodnum ++; /* Assign number to end node */

            indedgetax[indedgenum ++] = orgindxtax[orgvnodend]; /* Add node to element edge sub-array */

            for (hashnum = (orgindxtax[orgvnodend] * MESHINDUCESEPAHASHPRIME) & hashmsk; 
                 hashtab[hashnum].orgvelmnum == orgvelmnum; hashnum = (hashnum + 1) & hashmsk) ;
            hashtab[hashnum].orgvelmnum = orgvelmnum; /* Add vertex to hash table */
            hashtab[hashnum].indvnodnum = orgindxtax[orgvnodend];
          }
        }
        indmeshptr->vendtax[indvelmnum] = indedgenum; /* Set end of element edge sub-array */

        for (indenodtmp = indenodnum - 1; indenodtmp >= indmeshptr->verttax[orgindxtax[orgvnodnum]]; indenodtmp --) {
          Gnum                indvelmtmp;         /* Number of current already declared element for current node         */
          Gnum                indeelmtmp;         /* Number of current edge of already declared element for current node */
          Gnum                mtchnbr;            /* Number of matches between new element and current element           */

          indvelmtmp = indedgetax[indenodtmp];

          for (indeelmtmp = indmeshptr->verttax[indvelmtmp], mtchnbr = 0;
               indeelmtmp < indmeshptr->vendtax[indvelmtmp]; indeelmtmp ++) {
            Gnum                indvnodend;
            Gnum                hashnum;

            indvnodend = indedgetax[indeelmtmp];  /* Get number of current end node of already declared element */

            for (hashnum = (indvnodend * MESHINDUCESEPAHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
              if (hashtab[hashnum].orgvelmnum != orgvelmnum) /* If end node not present */
                break;
              if (hashtab[hashnum].indvnodnum == indvnodend) { /* If end node found */
                mtchnbr ++;                       /* One more matching node         */
                break;
              }
            }
          }

          if (mtchnbr == (indedgenum - indmeshptr->verttax[indvelmnum])) { /* If we are useless */
            orgindxtax[orgvelmnum] = -2;          /* Never consider this element again          */
            indedgenum = indedgetmp;              /* Recycle edge sub-array of new element      */
            break;                                /* No need to go any further                  */
          }
          if (mtchnbr == (indmeshptr->vendtax[indvelmtmp] - indmeshptr->verttax[indvelmtmp])) { /* If other element is useless */
            indedgenbr -= mtchnbr;                /* Remove its edges */
            indmeshptr->verttax[indvelmtmp] = indmeshptr->verttax[indvelmnum]; /* Recycle it into our element */
            indmeshptr->vendtax[indvelmtmp] = indmeshptr->vendtax[indvelmnum];
            orgindxtax[indvneltax[indvelmtmp]] = -2;
            indvneltax[indvelmtmp] = orgvelmnum;
          }
        }
        if (indenodtmp < indmeshptr->verttax[orgindxtax[orgvnodnum]]) { /* If new element distinct from all other neighboring elements */
          indedgetax[indenodnum ++] = indvelmnum; /* Record element in edge sub-array of current node                                  */
          indedgenbr += indedgenum - indmeshptr->verttax[indvelmnum];
          indvneltax[indvelmnum] = orgvelmnum;    /* Record original number of element in case we have to remove it afterwards */
          orgindxtax[orgvelmnum] = indvelmnum ++;
        }
      }
      else                                        /* Element already exists                                     */
        indedgetax[indenodnum ++] = orgindxtax[orgvelmnum]; /* Record element in edge sub-array of current node */
    }
    indmeshptr->vendtax[orgindxtax[orgvnodnum]] = indenodnum; /* Set end of edge sub-array for current node */
    indedgenbr += (indenodnum - indmeshptr->verttax[orgindxtax[orgvnodnum]]); /* Account for thes edges     */
  }
#ifdef SCOTCH_DEBUG_MESH2
  if (indvnodnum != (orgsepanbr + indmeshptr->baseval)) {
    errorPrint ("meshInduceSepa: internal error");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  indmeshptr->velmnbr   = indvelmnum - indmeshptr->velmbas;
  indmeshptr->velmnnd   = indvelmnum;
  indmeshptr->edgenbr = indedgenbr;

  if ((indmeshptr->edgetax = (Gnum *) memRealloc (indedgetax + indmeshptr->baseval, (indedgenum - indmeshptr->baseval) * sizeof (Gnum))) == NULL) {
    errorPrint ("meshInduceSepa: out of memory (3)");
    memFree    (indedgetax + indmeshptr->baseval); /* Free group leader */
    meshFree   (indmeshptr);
    return     (1);
  }
  indmeshptr->edgetax -= indmeshptr->baseval;

  indmeshptr->velosum = indmeshptr->velmnbr;      /* TODO: account for element weights */
  indmeshptr->vnlosum = (indmeshptr->velotax == NULL) ? orgsepanbr : indvnlosum;
  indmeshptr->degrmax = orgmeshptr->degrmax;      /* Rough estimate */

  if (orgmeshptr->vnumtax != NULL) {              /* If mesh is a submesh */
    for (indvnodnum = indmeshptr->vnodbas; indvnodnum < indmeshptr->vnodnnd; indvnodnum ++)
      indmeshptr->vnumtax[indvnodnum] = orgmeshptr->vnumtax[indmeshptr->vnumtax[indvnodnum] + (orgmeshptr->vnodbas - orgmeshptr->baseval)];
  }

  return (0);
}
