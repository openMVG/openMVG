/* Copyright 2004,2007,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : library_mesh_f.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the source mesh handling routines of    **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 23 sep 2002     **/
/**                                 to     13 dec 2005     **/
/**                # Version 5.1  : from : 27 mar 2010     **/
/**                                 to     15 apr 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include "module.h"
#include "common.h"
#include "scotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the mesh handling routines.    */
/*                                    */
/**************************************/

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHINIT, scotchfmeshinit, (             \
SCOTCH_Mesh * const         meshptr,            \
int * const                 revaptr),           \
(meshptr, revaptr))
{
  *revaptr = SCOTCH_meshInit (meshptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHEXIT, scotchfmeshexit, (             \
SCOTCH_Mesh * const         meshptr),           \
(meshptr))
{
  SCOTCH_meshExit (meshptr);
}

/* When an input stream is built from the given
** file handle, it is set as unbuffered, so as to
** allow for multiple stream reads from the same
** file handle. If it were buffered, too many
** input characters would be read on the first
** block read.
*/

FORTRAN (                                       \
SCOTCHFMESHLOAD, scotchfmeshload, (             \
SCOTCH_Mesh * const         meshptr,            \
int * const                 fileptr,            \
const SCOTCH_Num * const    baseptr,            \
int * const                 revaptr),           \
(meshptr, fileptr, baseptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFMESHLOAD: cannot duplicate handle");
    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "r")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFMESHLOAD: cannot open input stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }
  setbuf (stream, NULL);                          /* Do not buffer on input */

  o = SCOTCH_meshLoad (meshptr, stream, *baseptr);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHSAVE, scotchfmeshsave, (             \
const SCOTCH_Mesh * const   meshptr,            \
int * const                 fileptr,            \
int * const                 revaptr),           \
(meshptr, fileptr, revaptr))
{
  FILE *              stream;                     /* Stream to build from handle */
  int                 filenum;                    /* Duplicated handle           */
  int                 o;

  if ((filenum = dup (*fileptr)) < 0) {           /* If cannot duplicate file descriptor */
    errorPrint ("SCOTCHFMESHSAVE: cannot duplicate handle");

    *revaptr = 1;                                 /* Indicate error */
    return;
  }
  if ((stream = fdopen (filenum, "w")) == NULL) { /* Build stream from handle */
    errorPrint ("SCOTCHFMESHSAVE: cannot open output stream");
    close      (filenum);
    *revaptr = 1;
    return;
  }

  o = SCOTCH_meshSave (meshptr, stream);

  fclose (stream);                                /* This closes filenum too */

  *revaptr = o;
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHBUILD, scotchfmeshbuild, (           \
SCOTCH_Mesh * const         meshptr,            \
const SCOTCH_Num * const    velmbas,            \
const SCOTCH_Num * const    vnodbas,            \
const SCOTCH_Num * const    velmnbr,            \
const SCOTCH_Num * const    vnodnbr,            \
const SCOTCH_Num * const    verttab,            \
const SCOTCH_Num * const    vendtab,            \
const SCOTCH_Num * const    velotab,            \
const SCOTCH_Num * const    vnlotab,            \
const SCOTCH_Num * const    vlbltab,            \
const SCOTCH_Num * const    edgenbr,            \
const SCOTCH_Num * const    edgetab,            \
int * const                 revaptr),           \
(meshptr, velmbas, vnodbas, velmnbr, vnodnbr,   \
 verttab, vendtab, velotab, vnlotab, vlbltab,   \
 edgenbr, edgetab, revaptr))
{
  *revaptr = SCOTCH_meshBuild (meshptr, *velmnbr, *vnodnbr, *velmbas, *vnodbas,
                               verttab, vendtab, velotab, vnlotab, vlbltab,
                               *edgenbr, edgetab);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHCHECK, scotchfmeshcheck, (           \
const SCOTCH_Mesh * const   meshptr,            \
int * const                 revaptr),           \
(meshptr, revaptr))
{
  *revaptr = SCOTCH_meshCheck (meshptr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHSIZE, scotchfmeshsize, (             \
const SCOTCH_Mesh * const  meshptr,             \
SCOTCH_Num * const          velmnbr,            \
SCOTCH_Num * const          vnodnbr,            \
SCOTCH_Num * const          edgenbr),           \
(meshptr, velmnbr, vnodnbr, edgenbr))
{
  SCOTCH_meshSize (meshptr, velmnbr, vnodnbr, edgenbr);
}

/*
**
*/

FORTRAN (                                       \
SCOTCHFMESHDATA, scotchfmeshdata, (             \
const SCOTCH_Mesh * const   meshptr,            \
const SCOTCH_Num * const    indxptr,            \
SCOTCH_Num * const          velmbas,            \
SCOTCH_Num * const          vnodbas,            \
SCOTCH_Num * const          velmnbr,            \
SCOTCH_Num * const          vnodnbr,            \
SCOTCH_Idx * const          vertidx,            \
SCOTCH_Idx * const          vendidx,            \
SCOTCH_Idx * const          veloidx,            \
SCOTCH_Idx * const          vnloidx,            \
SCOTCH_Idx * const          vlblidx,            \
SCOTCH_Num * const          edgenbr,            \
SCOTCH_Idx * const          edgeidx,            \
SCOTCH_Num * const          degrnbr),           \
(meshptr, indxptr, velmbas, vnodbas, velmnbr, vnodnbr, \
 vertidx, vendidx, veloidx, vnloidx, vlblidx,          \
 edgenbr, edgeidx, degrnbr))
{
  SCOTCH_Num *        verttab;                    /* Pointer to mesh arrays */
  SCOTCH_Num *        vendtab;
  SCOTCH_Num *        velotab;
  SCOTCH_Num *        vnlotab;
  SCOTCH_Num *        vlbltab;
  SCOTCH_Num *        edgetab;

  SCOTCH_meshData (meshptr, velmbas, vnodbas, velmnbr, vnodnbr,
                   &verttab, &vendtab, &velotab, &vnlotab, &vlbltab,
                   edgenbr, &edgetab, degrnbr);
  *vertidx = (verttab - indxptr) + 1;             /* Add 1 since Fortran indices start at 1 */
  *vendidx = (vendtab - indxptr) + 1;
  *veloidx = (velotab != NULL) ? (velotab - indxptr) + 1 : *vertidx;
  *vnloidx = (vnlotab != NULL) ? (vnlotab - indxptr) + 1 : *vertidx;
  *vlblidx = (vlbltab != NULL) ? (vlbltab - indxptr) + 1 : *vertidx;
  *edgeidx = (edgetab - indxptr) + 1;
}
