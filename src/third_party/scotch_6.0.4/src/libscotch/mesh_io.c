/* Copyright 2004,2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : mesh_io.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the source graph    **/
/**                input/output functions.                 **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 05 nov 2002     **/
/**                                 to     06 may 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to     27 feb 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to     11 aug 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH_IO

#include "module.h"
#include "common.h"
#include "graph.h"
#include "graph_io.h"
#include "mesh.h"
#include "mesh_io.h"

/******************************************/
/*                                        */
/* These routines handle source mesh I/O. */
/*                                        */
/******************************************/

/* This routine loads a source mesh from
** the given stream.
** Edge loads, whenever present, are
** always discarded.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
meshLoad (
Mesh * restrict const       meshptr,              /* Mesh structure to fill               */
FILE * restrict const       stream,               /* Stream from which to read graph data */
const Gnum                  baseval)              /* Base value (-1 means keep file base) */
{
  Gnum                edgenum;                    /* Number of edges really allocated */
  Gnum                edgennd;
  Gnum                vertnbr;
  Gnum                velmnbr;                    /* Number of elements in mesh       */
  Gnum                velmbas;                    /* Base index for element vertices  */
  Gnum                vnodnbr;                    /* Number of nodes in mesh          */
  Gnum                vnodbas;                    /* Base index for node vertices     */
  Gnum                vlblmax;                    /* Maximum vertex label number      */
  Gnum                vlblval;                    /* Value where to read vertex label */
  Gnum                velonbr;                    /* Size of vertex load array        */
  Gnum                veloval;                    /* Value where to read vertex load  */
  Gnum                vlblnbr;                    /* Size of vertex label array       */
  Gnum                edloval;                    /* Value where to read edge load    */
  Gnum                edgeval;                    /* Value where to read edge end     */
  Gnum                baseadj;
  Gnum                versval;
  Gnum                degrmax;
  Gnum                propval;
  char                proptab[4];
  Gnum                vertbastab[2];
  Gnum                vertnndtab[2];
  Gnum                edgeadjtab[2];
  int                 i;

  memSet (meshptr, 0, sizeof (Mesh));

  if ((intLoad (stream, &versval) != 1) ||        /* Read version number */
      (versval != 1)) {
    errorPrint ("meshLoad: bad input (1)");
    return     (1);
  }

  if ((intLoad (stream, &velmnbr)          != 1) || /* Read rest of header */
      (intLoad (stream, &vnodnbr)          != 1) ||
      (intLoad (stream, &meshptr->edgenbr) != 1) ||
      (intLoad (stream, &velmbas)          != 1) ||
      (intLoad (stream, &vnodbas)          != 1) ||
      (intLoad (stream, &propval)          != 1) ||
      (velmnbr < 0)                              ||
      (vnodnbr < 0)                              ||
      (velmbas < 0)                              ||
      (vnodbas < 0)                              ||
      (propval < 0)                              ||
      (propval > 111)                            ||
      (((velmbas + velmnbr) != vnodbas) &&
       ((vnodbas + vnodnbr) != velmbas))) {
    errorPrint ("meshLoad: bad input (2)");
    return     (1);
  }
  sprintf (proptab, "%3.3d", (int) propval);      /* Compute file properties */
  proptab[0] -= '0';                              /* Vertex labels flag      */
  proptab[1] -= '0';                              /* Edge weights flag       */
  proptab[2] -= '0';                              /* Vertex loads flag       */

  baseadj = MIN (velmbas, vnodbas);               /* Get file graph base value   */
  if (baseval == -1) {                            /* If keep file graph base     */
    meshptr->baseval = baseadj;                   /* Set graph base as file base */
    baseadj            = 0;                       /* No base adjustment needed   */
  }
  else {                                          /* If set graph base     */
    meshptr->baseval = baseval;                   /* Set wanted graph base */
    baseadj          = baseval - baseadj;         /* Update base adjust    */
  }
  meshptr->flagval = MESHFREEVERT | MESHVERTGROUP; /* Edge array grouped with vertex array */
  meshptr->velmnbr = velmnbr;
  meshptr->velmbas = velmbas + baseadj;
  meshptr->velmnnd = velmnbr + (velmbas + baseadj);
  meshptr->vnodnbr = vnodnbr;
  meshptr->vnodbas = vnodbas + baseadj;
  meshptr->vnodnnd = vnodnbr + (vnodbas + baseadj);
  vertnbr          = velmnbr + vnodnbr;

  velonbr = (proptab[2] != 0) ? vertnbr : 0;
  vlblnbr = (proptab[0] != 0) ? vertnbr : 0;

  if (memAllocGroup ((void **) (void *)
      &meshptr->verttax, (size_t) ((vertnbr + 1)     * sizeof (Gnum)),
      &meshptr->vlbltax, (size_t) ( vlblnbr          * sizeof (Gnum)),
      &meshptr->velotax, (size_t) ( velonbr          * sizeof (Gnum)), /* Allocate single array for both element and node vertices */
      &meshptr->edgetax, (size_t) ( meshptr->edgenbr * sizeof (Gnum)), NULL) == NULL) { /* Edge array grouped with vertex arrays   */
    errorPrint ("meshLoad: out of memory (1)");
    meshFree   (meshptr);
    return     (1);
  }
  meshptr->verttax -= meshptr->baseval;
  meshptr->vendtax  = meshptr->verttax + 1;       /* Use compact vertex array */
  meshptr->velotax  = (velonbr != 0) ? (meshptr->velotax - meshptr->baseval) : NULL; /* Store based load array access in velotax */
  meshptr->vnlotax  = meshptr->velotax;
  meshptr->vlbltax  = (vlblnbr != 0) ? (meshptr->vlbltax - meshptr->baseval) : NULL;
  meshptr->velosum  = meshptr->velmnbr;           /* Assume element and node vertices not weighted */
  meshptr->vnlosum  = meshptr->vnodnbr;
  meshptr->edgetax -= meshptr->baseval;

  degrmax = 0;
  edgennd = meshptr->edgenbr + meshptr->baseval;
  edgenum = meshptr->baseval;                     /* No edges allocated yet */
  vlblmax = vertnbr + meshptr->baseval - 1;       /* No vertex labels known */

  if (meshptr->velmbas <= meshptr->vnodbas) {     /* If elements first */
    vertbastab[0] = meshptr->velmbas;
    vertnndtab[0] = meshptr->velmnnd;
    edgeadjtab[0] = meshptr->vnodbas - meshptr->baseval;
    vertbastab[1] = meshptr->vnodbas;
    vertnndtab[1] = meshptr->vnodnnd;
    edgeadjtab[1] = meshptr->velmbas - meshptr->baseval;
  }
  else {
    vertbastab[0] = meshptr->vnodbas;
    vertnndtab[0] = meshptr->vnodnnd;
    edgeadjtab[0] = meshptr->velmbas - meshptr->baseval;
    vertbastab[1] = meshptr->velmbas;
    vertnndtab[1] = meshptr->velmnnd;
    edgeadjtab[1] = meshptr->vnodbas - meshptr->baseval;
  }
  for (i = 0; i < 2; i ++) {                      /* For both kinds of vertices */
    Gnum                vertbas;
    Gnum                vertnnd;
    Gnum                edgeadj;
    Gnum                vertnum;
    Gnum                velosum;
    Gnum                velomax;

    vertbas = vertbastab[i];
    vertnnd = vertnndtab[i];
    edgeadj = edgeadjtab[i];
    velosum = 0;
    velomax = 1;                                  /* Assume vertex loads all equal to 1 */

    for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) { /* For all vertices of same kind */
      Gnum                degrval;

      if (meshptr->vlbltax != NULL) {             /* If must read label */
        if (intLoad (stream, &vlblval) != 1) {    /* Read label data    */
          errorPrint ("meshLoad: bad input (3)");
          meshFree   (meshptr);
          return     (1);
        }
        meshptr->vlbltax[vertnum] = vlblval + vertbas + baseadj; /* Adjust vertex label */
        if (meshptr->vlbltax[vertnum] > vlblmax)  /* Get maximum vertex label           */
          vlblmax = meshptr->vlbltax[vertnum];
      }
      if (proptab[2] != 0) {                      /* If must read vertex load */
        if ((intLoad (stream, &veloval) != 1) ||  /* Read vertex load data    */
            (veloval < 1)) {
          errorPrint ("meshLoad: bad input (4)");
          meshFree   (meshptr);
          return     (1);
        }
        if (veloval > velomax)
          velomax = veloval;
        meshptr->velotax[vertnum] = veloval;
        velosum += veloval;
      }
      if (intLoad (stream, &degrval) != 1) {      /* Read vertex degree */
        errorPrint ("meshLoad: bad input (5)");
        meshFree   (meshptr);
        return     (1);
      }
      if (degrmax < degrval)                      /* Set maximum degree */
        degrmax = degrval;

      meshptr->verttax[vertnum] = edgenum;        /* Set index in edge array */
      degrval += edgenum;
      if (degrval > edgennd) {                    /* Check if edge array overflows */
        errorPrint ("meshLoad: invalid arc count (1)");
        meshFree   (meshptr);
        return     (1);
      }

      for ( ; edgenum < degrval; edgenum ++) {
        if (proptab[1] != 0) {                    /* If must read edge load        */
          if (intLoad (stream, &edloval) != 1) {  /* Read edge load data (useless) */
            errorPrint ("meshLoad: bad input (6)");
            meshFree   (meshptr);
            return     (1);
          }
        }
        if (intLoad (stream, &edgeval) != 1) {    /* Read edge data */
          errorPrint ("meshLoad: bad input (7)");
          meshFree   (meshptr);
          return     (1);
        }
        meshptr->edgetax[edgenum] = edgeval + edgeadj;
      }
    }

    if (vertbastab[i] == meshptr->velmbas) {      /* If elements are processed        */
      if (velomax == 1)                           /* If element loads not significant */
        meshptr->velotax = NULL;
      else
        meshptr->velosum = velosum;
    }
    else {
      if (velomax == 1)                           /* If node loads not significant */
        meshptr->vnlotax = NULL;
      else
        meshptr->vnlosum = velosum;
    }
  }
  meshptr->verttax[vertnbr + meshptr->baseval] = meshptr->edgenbr + meshptr->baseval; /* Set end of edge array */

  if (edgenum != edgennd) {                       /* Check if number of edges is valid */
    errorPrint ("meshLoad: invalid arc count (2)");
    meshFree   (meshptr);
    return     (1);
  }

  meshptr->degrmax = degrmax;

  if (meshptr->vlbltax != NULL) {                 /* If vertex label renaming necessary                 */
    if (graphLoad2 (meshptr->baseval, vertnbr + meshptr->baseval, meshptr->verttax, /* Rename edge ends */
                    meshptr->vendtax, meshptr->edgetax, vlblmax, meshptr->vlbltax) != 0) {
      errorPrint ("meshLoad: cannot relabel vertices");
      meshFree   (meshptr);
      return     (1);
    }
  }

#ifdef SCOTCH_DEBUG_MESH2
  if (meshCheck (meshptr) != 0) {                 /* Check mesh consistency */
    errorPrint  ("meshLoad: inconsistent mesh data");
    meshFree    (meshptr);
    return      (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  return (0);
}

/* This routine saves a source mesh to
** the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
meshSave (
const Mesh * restrict const meshptr,
FILE * restrict const       stream)
{
  char                propstr[4];                 /* Property string */
  Gnum                vertbastab[2];
  Gnum                vertnndtab[2];
  Gnum *              velotabtab[2];
  Gnum                edgeadjtab[2];
  int                 i;
  int                 o;

  propstr[0] = (meshptr->vlbltax != NULL) ? '1' : '0'; /* Set property string */
  propstr[1] = '0';                               /* No edge loads written    */
  propstr[2] = ((meshptr->velotax != NULL) || (meshptr->vnlotax != NULL)) ? '1' : '0';
  propstr[3] = '\0';

  if (fprintf (stream, "1\n" GNUMSTRING "\t" GNUMSTRING "\t" GNUMSTRING "\n" GNUMSTRING "\t" GNUMSTRING "\t%3s\n", /* Write file header */
               (Gnum) meshptr->velmnbr,
               (Gnum) meshptr->vnodnbr,
               (Gnum) meshptr->edgenbr,
               (Gnum) meshptr->velmbas,
               (Gnum) meshptr->vnodbas,
               propstr) == EOF) {
    errorPrint ("meshSave: bad output (1)");
    return     (1);
  }

  vertbastab[0] = meshptr->baseval;
  vertnndtab[1] = meshptr->velmnbr + meshptr->vnodnbr + meshptr->baseval;
  if (meshptr->velmbas <= meshptr->vnodbas) {     /* If elements first */
    vertnndtab[0] = meshptr->velmnnd;
    velotabtab[0] = meshptr->velotax;
    edgeadjtab[0] = meshptr->vnodbas - meshptr->baseval;
    vertbastab[1] = meshptr->vnodbas;
    velotabtab[0] = meshptr->vnlotax;
    edgeadjtab[1] = meshptr->velmbas - meshptr->baseval;
  }
  else {
    vertnndtab[0] = meshptr->vnodnnd;
    velotabtab[0] = meshptr->vnlotax;
    edgeadjtab[0] = meshptr->velmbas - meshptr->baseval;
    vertbastab[1] = meshptr->velmbas;
    velotabtab[1] = meshptr->velotax;
    edgeadjtab[1] = meshptr->vnodbas - meshptr->baseval;
  }
  for (i = 0; i < 2; i ++) {                      /* For both kinds of vertices */
    Gnum                vertbas;
    Gnum                vertnnd;
    Gnum * restrict     velotax;
    Gnum                edgeadj;
    Gnum                vertnum;

    vertbas = vertbastab[i];
    vertnnd = vertnndtab[i];
    velotax = velotabtab[i];
    edgeadj = edgeadjtab[i];

    for (vertnum = vertbas, o = 0; (vertnum < vertnnd) && (o == 0); vertnum ++) {
      Gnum                edgenum;

      if (meshptr->vlbltax != NULL)               /* Write vertex label if necessary */
        o  = (fprintf (stream, GNUMSTRING "\t", (Gnum) meshptr->vlbltax[vertnum]) == EOF);
      if (propstr[2] != '0')                      /* Write vertex load if necessary */
        o |= (fprintf (stream, GNUMSTRING "\t", (Gnum) ((velotax != NULL) ? velotax[vertnum] : 1)) == EOF);
      o |= (fprintf (stream, GNUMSTRING, (Gnum) (meshptr->vendtax[vertnum] - meshptr->verttax[vertnum])) == EOF); /* Write vertex degree */

      for (edgenum = meshptr->verttax[vertnum];
           (edgenum < meshptr->vendtax[vertnum]) && (o == 0); edgenum ++) {
        o |= (putc ('\t', stream) == EOF);
        o |= (intSave (stream,                    /* Write edge end */
                       (meshptr->vlbltax != NULL) ? meshptr->vlbltax[meshptr->edgetax[edgenum]] : meshptr->edgetax[edgenum] - edgeadj) != 1);
      }
      o |= (putc ('\n', stream) == EOF);
    }
  }

  if (o != 0)
    errorPrint ("meshSave: bad output (2)");

  return (o);
}
