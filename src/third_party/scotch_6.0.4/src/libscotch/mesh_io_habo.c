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
/**   NAME       : mesh_io_habo.c                          **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the I/O routines   **/
/**                for handling the Harwell-Boeing         **/
/**                elemental matrix format.                **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 jan 2004     **/
/**                                 to     20 jan 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define MESH_IO_HABO

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "mesh.h"

/* This routine loads the geometrical mesh
** in the Harwell-Boeing matrix format, and
** allocates the proper structures.
** - 0   : on success.
** - !0  : on error.
*/

int
meshGeomLoadHabo (
Mesh * restrict const       meshptr,              /* Mesh to load     */
Geom * restrict const       geomptr,              /* Geometry to load */
FILE * const                filesrcptr,           /* Topological data */
FILE * const                filegeoptr,           /* No use           */
const char * const          dataptr)              /* No use           */
{
  long              habmattag;                    /* Matrix tag number in file                  */
  long              habmatnum;                    /* Current matrix number                      */
  char              habmatbuf[4][84];             /* Matrix header line buffers                 */
  char              habmattype[4];                /* Matrix type                                */
  long              habcrdnbr;                    /* Total number of data lines                 */
  long              habrhsnbr;                    /* Number of right hand side lines            */
  int               habcolsiz;                    /* Size of column type                        */
  int               habnzrsiz;                    /* Size of non-zero type                      */
  Gnum              habevlnbr;                    /* Number of element matrix values (not uesd) */
  Gnum              habvixnbr;                    /* Number of element matrix variable indices  */
  Gnum              habeltnbr;                    /* Number of elements                         */
  Gnum              habvarnbr;                    /* Number of variables, that is, nodes        */
  Gnum              velmnum;                      /* Number of current element vertex           */
  Gnum              vnodnum;                      /* Number of current node vertex              */
  Gnum              edgenum;                      /* Number of current edge (arc)               */
  Gnum              edgesum;                      /* Accumulated number of edges in node part   */
  Gnum              degrmax;                      /* Maximum degree                             */
  int               c;

  if (((habmattag = atol (dataptr)) == 0) &&      /* Get tag value */
      (dataptr[0] != '0') && (dataptr[0] != '\0')) {
    errorPrint ("meshGeomLoadHabo: bad input (1)");
    return     (1);
  }

  for (habmatnum = 0; habmatnum <= habmattag; habmatnum ++) { /* Read headers and skip if necessary */
    memSet (habmatbuf[0], ' ', &habmatbuf[3][83] - &habmatbuf[0][0]); /* Initialize header buffers  */
    if ((fgets (habmatbuf[0], 83, filesrcptr) == NULL) || /* Read mesh header                       */
        (fgets (habmatbuf[1], 83, filesrcptr) == NULL) ||
        (fgets (habmatbuf[2], 83, filesrcptr) == NULL) ||
        (fgets (habmatbuf[3], 83, filesrcptr) == NULL)) {
      errorPrint ("meshGeomLoadHabo: bad input (2)");
      return     (1);
    }
    habmatbuf[1][70] = '\0';                      /* Extract header values */
    habrhsnbr = atol (&habmatbuf[1][56]);
    habmatbuf[1][14] = '\0';
    habcrdnbr = atol (&habmatbuf[1][00]);
    habmattype[0] = toupper (habmatbuf[2][0]);
    habmattype[1] = toupper (habmatbuf[2][1]);
    habmattype[2] = toupper (habmatbuf[2][2]);
    habmatbuf[2][70] = '\0';
    habevlnbr = (Gnum) atol (&habmatbuf[2][57]);  /* NELTVL */
    habmatbuf[2][56] = '\0';
    habvixnbr = (Gnum) atol (&habmatbuf[2][43]);  /* VARIX */
    habmatbuf[2][42] = '\0';
    habeltnbr = (Gnum) atol (&habmatbuf[2][29]);  /* NELT */
    habmatbuf[2][28] = '\0';
    habvarnbr = (Gnum) atol (&habmatbuf[2][14]);  /* NVAR */
    habmatbuf[2][14] = '\0';

    if ((c = sscanf (habmatbuf[3], "(%*d%*[Ii]%d) (%*d%*[Ii]%d)",
                     &habcolsiz, &habnzrsiz)) != 2) {
      errorPrint ("meshGeomLoadHabo: bad input (3, %d)", c);
      return     (1);
    }

    if (habrhsnbr != 0) {
      while ((c = getc (filesrcptr)) != '\n'){    /* Skip RHS format line */
        if (c == EOF) {
          errorPrint ("meshGeomLoadHabo: bad input (4)");
          return     (1);
        }
      }
    }

    if (habmatnum < habmattag) {                 /* If we have to skip file */
      while (habcrdnbr -- > 0) {                 /* Skip all of file lines  */
        while ((c = getc (filesrcptr)) != '\n') { /* Skip line              */
          if (c == EOF) {
            errorPrint ("meshGeomLoadHabo: bad input (5)");
            return     (1);
          }
        }
      }
    }
  }

  if (habmattype[2] != 'E') {
    errorPrint ("meshGeomLoadHabo: only elemental matrices supported");
    return     (1);
  }
  if (habmattype[1] == 'R') {
    errorPrint ("meshGeomLoadHabo: rectangular matrices not supported");
    return     (1);
  }

  if (((meshptr->verttax = (Gnum *) memAlloc ((habeltnbr + habvarnbr + 1) * sizeof (Gnum))) == NULL) ||
      ((meshptr->edgetax = (Gnum *) memAlloc (habvixnbr * 2               * sizeof (Gnum))) == NULL)) {
    errorPrint ("meshGeomLoadHabo: out of memory (1)");
    if (meshptr->verttax != NULL) {
      memFree (meshptr->verttax);
      meshptr->verttax = NULL;
    }
    return (1);
  }

  meshptr->flagval = MESHFREETABS;                /* Totally new mesh structure       */
  meshptr->baseval = 1;                           /* Harwell-Boeing meshs have base 1 */
  meshptr->vendtax = meshptr->verttax;
  meshptr->verttax --;
  meshptr->edgenbr = habvixnbr * 2;
  meshptr->edgetax --;
  meshptr->velmnbr = habeltnbr;
  meshptr->velmbas = habvarnbr + 1;
  meshptr->velmnnd = habeltnbr + habvarnbr + 1;
  meshptr->vnodnbr = habvarnbr;
  meshptr->vnodbas = 1;
  meshptr->vnodnnd = meshptr->velmbas;
  meshptr->vnlosum = habvarnbr;

  for (velmnum = meshptr->velmbas; velmnum <= meshptr->velmnnd; velmnum ++) { /* Read ELTPTR as second part of vertex array */
    Gnum              habcolval;                  /* Current column value                                                   */
    int               habcolidx;                  /* Current index in column value                                          */

    while (((c = getc (filesrcptr)) == '\n') || (c == '\r')) ;
    habcolval = (c == ' ') ? 0 : (c - '0');
    for (habcolidx = 1; habcolidx < habcolsiz; habcolidx ++) {
      if ((c = getc (filesrcptr)) != ' ')
        habcolval = habcolval * 10 + c - '0';
    }
    if (c == EOF) {
      errorPrint ("meshGeomLoadHabo: bad input (6)");
      meshFree   (meshptr);
      return     (1);
    }
#ifdef SCOTCH_DEBUG_MESH2
    if ((habcolval < 1) || (habcolval > (habvixnbr + 1))) {
      errorPrint ("meshGeomLoadHabo: bad input (7)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_MESH2 */

    meshptr->verttax[velmnum] = habcolval + habvixnbr;
  }
  if (meshptr->verttax[velmnum - 1] != (habvixnbr + habvixnbr + 1)) {
    errorPrint ("meshGeomLoadHabo: bad input (8)");
    meshFree   (meshptr);
    return     (1);
  }

  memSet (meshptr->verttax + 1, 0, habvarnbr * sizeof (Gnum)); /* Pre-set node adjacency array */

  for (edgenum = habvixnbr + 1; edgenum <= meshptr->edgenbr; edgenum ++) { /* Read VARIND as second part of edge array */
    Gnum              habnodval;                  /* Current non-zero value                                              */
    int               habnzridx;                  /* Current index in non-zero value                                     */

    while (((c = getc (filesrcptr)) == '\n') || (c == '\r')) ;
    habnodval = (c == ' ') ? 0 : (c - '0');
    for (habnzridx = 1; habnzridx < habnzrsiz; habnzridx ++) {
      if ((c = getc (filesrcptr)) != ' ')
        habnodval = habnodval * 10 + c - '0';
    }
    if (c == EOF) {
      errorPrint ("meshGeomLoadHabo: bad input (9)");
      meshFree   (meshptr);
      return     (1);
    }
#ifdef SCOTCH_DEBUG_MESH2
    if ((habnodval < 1) || (habnodval > habvarnbr)) {
      errorPrint ("meshGeomLoadHabo: bad input (10)");
      return     (1);
    }
#endif /* SCOTCH_DEBUG_MESH2 */
    meshptr->edgetax[edgenum] = habnodval;
    meshptr->verttax[habnodval] ++;
  }

  degrmax = 1;
  for (vnodnum = edgesum = 1; vnodnum < meshptr->vnodnnd; vnodnum ++) { /* Accumulate start indices for node part of vertex array */
    Gnum                degrval;

    degrval = meshptr->verttax[vnodnum];
    if (degrval > degrmax)
      degrmax = degrval;

    meshptr->verttax[vnodnum] = edgesum;
    edgesum += degrval;
  }
#ifdef SCOTCH_DEBUG_MESH2
  if (edgesum != meshptr->verttax[meshptr->velmbas]) {
    errorPrint ("meshGeomLoadHabo: internal error (1)");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  for (velmnum = meshptr->velmbas, edgenum = habvixnbr + 1; /* Re-scan element part of edge array to build node part of edge array */
       velmnum < meshptr->velmnnd; velmnum ++) {
    if ((meshptr->vendtax[velmnum] - edgenum) > degrmax)
      degrmax = (meshptr->vendtax[velmnum] - edgenum);

    for ( ; edgenum < meshptr->vendtax[velmnum]; edgenum ++) {
      Gnum                vnodnum;

      vnodnum = meshptr->edgetax[edgenum];
#ifdef SCOTCH_DEBUG_MESH2
      if ((vnodnum < 1) || (vnodnum > habvarnbr)) {
        errorPrint ("meshGeomLoadHabo: internal error (2)");
        return     (1);
      }
      if (meshptr->verttax[vnodnum] > habvixnbr) {
        errorPrint ("meshGeomLoadHabo: internal error (3)");
        return     (1);
      }
#endif /* SCOTCH_DEBUG_MESH2 */
      meshptr->edgetax[meshptr->verttax[vnodnum] ++] = velmnum;
    }
  }

  meshptr->degrmax = degrmax;

  for (vnodnum = edgesum = 1; vnodnum < meshptr->vnodnnd; vnodnum ++) { /* Restore start indices for node part of vertex array */
    Gnum                edgenum;

    edgenum = meshptr->verttax[vnodnum];
    meshptr->verttax[vnodnum] = edgesum;
    edgesum = edgenum;
  }


#ifdef SCOTCH_DEBUG_MESH2
  if (meshCheck (meshptr) != 0) {
    errorPrint ("meshGeomLoadHabo: inconsistent mesh data");
    return     (1);
  }
#endif /* SCOTCH_DEBUG_MESH2 */

  return (0);
}
