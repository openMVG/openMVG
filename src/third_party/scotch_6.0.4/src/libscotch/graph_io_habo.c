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
/**   NAME       : graph_io_habo.c                         **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the I/O routines   **/
/**                for handling the Harwell-Boeing matrix  **/
/**                format.                                 **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 06 nov 1997     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 13 dec 1998     **/
/**                                 to     24 dec 1998     **/
/**                # Version 4.0  : from : 18 dec 2001     **/
/**                                 to     21 mar 2005     **/
/**                # Version 5.0  : from : 06 jun 2007     **/
/**                                 to     31 aug 2007     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to     27 apr 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define GRAPH_IO_HABO

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "graph_io_habo.h"

/* This routine loads the geometrical graph
** in the Harwell-Boeing matrix format, and
** allocates the proper structures.
** - 0   : on success.
** - !0  : on error.
*/

int
graphGeomLoadHabo (
Graph * restrict const      grafptr,              /* Graph to load    */
Geom * restrict const       geomptr,              /* Geometry to load */
FILE * const                filesrcptr,           /* Topological data */
FILE * const                filegeoptr,           /* No use           */
const char * const          dataptr)              /* Tag value        */
{
  Gnum                          habmattag;        /* Matrix tag number in file       */
  Gnum                          habmatnum;        /* Current matrix number           */
  char                          habmatbuf[4][84]; /* Matrix header line buffers      */
  char                          habmattype[3];    /* Matrix type                     */
  Gnum                          habcrdnbr;        /* Total number of data lines      */
  Gnum                          habrhsnbr;        /* Number of right hand side lines */
  Gnum                          habrownbr;        /* Number of rows                  */
  GraphGeomHaboLine             habcolfmt;        /* Format of column line           */
  int                           habvalnum;        /* Number of value in line         */
  Gnum                          habcolnbr;        /* Number of columns               */
  Gnum                          habcolnum;        /* Number of current column index  */
  Gnum * restrict               habcoltab;        /* Index array                     */
  GraphGeomHaboLine             habnzrfmt;        /* Format of non-zero type         */
  Gnum                          habnzrnbr;        /* Number of non-zero indices      */
  Gnum                          habnzrnum;        /* Number of current row data      */
  Gnum * restrict               habnzrtab;        /* Row data array                  */
  GraphGeomHaboHash * restrict  hashtab;          /* Neighbor hash table             */
  Gnum                          hashmsk;          /* Mask for access to hash table   */
  Gnum                          vertnum;          /* Number of current vertex        */
  Gnum                          edgenum;          /* Number of current edge (arc)    */
  Gnum                          edgeold;          /* Number of non-purged edge       */
  Gnum                          edgetmp;          /* Temporary edge number           */
  Gnum                          degrmax;          /* Maximum degree                  */
  int                           c;

  if ((dataptr != NULL)                          && /* If tag value provided */
      (dataptr[0] != '\0')                       &&
      ((habmattag = (Gnum) atol (dataptr)) == 0) && /* Get tag value */
      (dataptr[0] != '0')) {
    errorPrint ("graphGeomLoadHabo: invalid parameter");
    return     (1);
  }

  habmattype[0] =
  habmattype[1] =
  habmattype[2] = '\0';
  for (habmatnum = 0; habmatnum <= habmattag; habmatnum ++) { /* Read headers and skip if necessary */
    memSet (habmatbuf[0], ' ', &habmatbuf[3][83] - &habmatbuf[0][0]); /* Initialize header buffers  */
    if ((fgets (habmatbuf[0], 83, filesrcptr) == NULL) || /* Read graph header                      */
        (fgets (habmatbuf[1], 83, filesrcptr) == NULL) ||
        (fgets (habmatbuf[2], 83, filesrcptr) == NULL) ||
        (fgets (habmatbuf[3], 83, filesrcptr) == NULL)) {
      errorPrint ("graphGeomLoadHabo: bad input (1)");
      return     (1);
    }
    habmatbuf[1][70] = '\0';                      /* Extract header values */
    habrhsnbr = (Gnum) atol (&habmatbuf[1][56]);
    habmatbuf[1][14] = '\0';
    habcrdnbr = (Gnum) atol (&habmatbuf[1][00]);
    habmattype[0] = toupper (habmatbuf[2][0]);
    habmattype[1] = toupper (habmatbuf[2][1]);
    habmattype[2] = toupper (habmatbuf[2][2]);
    habmatbuf[2][56] = '\0';
    habnzrnbr = (Gnum) atol (&habmatbuf[2][43]);
    habmatbuf[2][42] = '\0';
    habcolnbr = (Gnum) atol (&habmatbuf[2][29]);
    habmatbuf[2][28] = '\0';
    habrownbr = (Gnum) atol (&habmatbuf[2][14]);

    habmatbuf[3][32] = '\0';
    if (graphGeomLoadHaboFormat (&habnzrfmt, &habmatbuf[3][16]) != 0) {
      errorPrint ("graphGeomLoadHabo: bad input (2)");
      return     (1);
    }
    habmatbuf[3][16] = '\0';
    if (graphGeomLoadHaboFormat (&habcolfmt, &habmatbuf[3][0]) != 0) {
      errorPrint ("graphGeomLoadHabo: bad input (3)");
      return     (1);
    }

    if (habrhsnbr != 0) {
      while ((c = getc (filesrcptr)) != '\n'){    /* Skip RHS format line */
        if (c == EOF) {
          errorPrint ("graphGeomLoadHabo: bad input (4)");
          return     (1);
        }
      }
    }

    if (habmatnum < habmattag) {                 /* If we have to skip file */
      while (habcrdnbr -- > 0) {                 /* Skip all of file lines  */
        while ((c = getc (filesrcptr)) != '\n') { /* Skip line              */
          if (c == EOF) {
            errorPrint ("graphGeomLoadHabo: bad input (5)");
            return     (1);
          }
        }
      }
    }
  }

  if (habmattype[2] != 'A') {
    errorPrint ("graphGeomLoadHabo: only assembled matrices supported; for unassembled matrices, use the mesh version of the tools");
    return     (1);
  }
  if (habmattype[1] == 'R') {
    errorPrint ("graphGeomLoadHabo: rectangular matrices not supported");
    return     (1);
  }

  if (((grafptr->verttax = (Gnum *) memAlloc ((habcolnbr + 1) * sizeof (Gnum))) == NULL) ||
      ((grafptr->edgetax = (Gnum *) memAllocGroup ((void **) (void *)
                                                   &grafptr->edgetax, (size_t) (habnzrnbr * 2   * sizeof (Gnum)),
                                                   &habcoltab,        (size_t) ((habcolnbr + 1) * sizeof (Gnum)),
                                                   &habnzrtab,        (size_t) (habnzrnbr       * sizeof (Gnum)), NULL)) == NULL)) {
    errorPrint ("graphGeomLoadHabo: out of memory (1)");
    if (grafptr->verttax != NULL) {
      memFree (grafptr->verttax);
      grafptr->verttax = NULL;
    }
    return (1);
  }
  grafptr->flagval = GRAPHFREETABS;               /* Totally new graph structure       */
  grafptr->baseval = 1;                           /* Harwell-Boeing graphs have base 1 */
  grafptr->vertnbr = (Gnum) habcolnbr;
  grafptr->vertnnd = grafptr->vertnbr + 1;
  grafptr->velosum = grafptr->vertnbr;
  grafptr->vendtax = grafptr->verttax;            /* Use compact representation for array based at 1     */
  grafptr->verttax --;                            /* Base verttab array at 1, with vendtab = verttab + 1 */
  grafptr->edgetax --;

  ungetc ('\n', filesrcptr);                      /* Create fake previous line     */
  for (habcolnum = 0, habvalnum = habcolfmt.datanbr; /* Eat up fake previous line  */
       habcolnum <= habcolnbr; habcolnum ++) {    /* Read column indices           */
    Gnum              habcolval;                  /* Current column value          */
    int               habcolidx;                  /* Current index in column value */

    c = getc (filesrcptr);
    if (habvalnum ++ >= habcolfmt.datanbr) {      /* If all useful data read from line */
      habvalnum = 1;                              /* Start at beginning of new line    */

      while ((c != '\n') && (c != '\r'))          /* Eat up all remaining spaces */
        c = getc (filesrcptr);
      while (((c = getc (filesrcptr)) == '\n') || (c == '\r')) ; /* Read till end of line */

      for (habcolidx = 0; habcolidx < habcolfmt.strtnbr; habcolidx ++) /* Get start of line */
        c = getc (filesrcptr);
    }

    habcolval = (c == ' ') ? 0 : (c - '0');
    for (habcolidx = 1; habcolidx < habcolfmt.datalen; habcolidx ++) {
      if ((c = getc (filesrcptr)) != ' ')
        habcolval = habcolval * 10 + c - '0';
    }
    if (c == EOF) {
      errorPrint ("graphGeomLoadHabo: bad input (6)");
      graphFree  (grafptr);
      return     (1);
    }
    habcoltab[habcolnum] = habcolval;
  }
  if (habcoltab[habcolnbr] != (Gnum) habnzrnbr + 1) {
    errorPrint ("graphGeomLoadHabo: bad input (7)");
    graphFree  (grafptr);
    return     (1);
  }

  memSet (grafptr->vendtax, 0, habcolnbr * sizeof (Gnum)); /* Here, vendtax = verttab */

  for (vertnum = 1, habnzrnum = 0, habvalnum = habnzrfmt.datanbr; /* Start by eating end of previous line  */
       vertnum < grafptr->vertnnd; vertnum ++) {  /* Read matrix pattern                                   */
    for ( ; habnzrnum < (habcoltab[vertnum] - 1); habnzrnum ++) { /* All right since vertnum is based at 1 */
      Gnum              habnzrval;                /* Current non-zero value                                */
      int               habnzridx;                /* Current index in non-zero value                       */

      c = getc (filesrcptr);
      if (habvalnum ++ >= habnzrfmt.datanbr) {    /* If all useful data read from line */
        habvalnum = 1;                            /* Start at beginning of new line    */

        while ((c != '\n') && (c != '\r'))        /* Eat up all remaining spaces */
          c = getc (filesrcptr);
        while (((c = getc (filesrcptr)) == '\n') || (c == '\r')) ; /* Read till end of line */

        for (habnzridx = 0; habnzridx < habnzrfmt.strtnbr; habnzridx ++) /* Get start of line */
          c = getc (filesrcptr);
      }
      habnzrval = (c == ' ') ? 0 : (c - '0');
      for (habnzridx = 1; habnzridx < habnzrfmt.datalen; habnzridx ++) {
        if ((c = getc (filesrcptr)) != ' ')
          habnzrval = habnzrval * 10 + c - '0';
      }
      if (c == EOF) {
        errorPrint ("graphGeomLoadHabo: bad input (8)");
        graphFree  (grafptr);
        return     (1);
      }
      habnzrtab[habnzrnum] = habnzrval;
      if (habnzrval != vertnum) {                 /* If not loop edge      */
        grafptr->verttax[vertnum] ++;             /* Account for arc       */
        grafptr->verttax[habnzrval] ++;           /* Add arc to symmetrize */
      }
    }
  }

  degrmax = 1;
  for (vertnum = edgenum = 1; vertnum < grafptr->vertnnd; vertnum ++) { /* Build (superset of) vertex array */
    Gnum              edgetmp;

    edgetmp = grafptr->verttax[vertnum];
    grafptr->verttax[vertnum] = edgenum;
    edgenum += edgetmp;
    if (edgetmp > degrmax)                        /* Update bound on maximum degree */
      degrmax = edgetmp;
  }
  grafptr->verttax[vertnum] = edgenum;            /* Set end of vertex array */

  for (vertnum = 1, habnzrnum = 0; vertnum < grafptr->vertnnd; vertnum ++) { /* Build (superset of) edge array */
    for ( ; habnzrnum < (habcoltab[vertnum] - 1); habnzrnum ++) {
      Gnum              vertend;                  /* Number of end vertex */

      vertend = habnzrtab[habnzrnum];
      if (vertend != vertnum) {                   /* If not loop edge           */
        grafptr->edgetax[grafptr->verttax[vertnum] ++] = vertend; /* Build arc  */
        grafptr->edgetax[grafptr->verttax[vertend] ++] = vertnum; /* Symmetrize */
      }
    }
  }

  for (hashmsk = 31; hashmsk < degrmax; hashmsk = hashmsk * 2 + 1) ; /* Set neighbor hash table size */
  hashmsk = hashmsk * 4 + 3;

  if ((hashtab = (GraphGeomHaboHash *) memAlloc ((hashmsk + 1) * sizeof (GraphGeomHaboHash))) == NULL) {
    errorPrint ("graphGeomLoadHabo: out of memory (2)");
    graphFree  (grafptr);
    return     (1);
  }
  memSet (hashtab, ~0, (hashmsk + 1) * sizeof (GraphGeomHaboHash)); /* Pre-set hash table */

  degrmax = 1;
  for (vertnum = edgetmp = edgenum = 1; vertnum < grafptr->vertnnd; vertnum ++) { /* Remove duplicates from edge array */
    for (edgeold = edgetmp, edgetmp = grafptr->verttax[vertnum], grafptr->verttax[vertnum] = edgenum;
         edgeold < edgetmp; edgeold ++) {
      Gnum              vertend;                  /* Number of end vertex */
      Gnum              hashnum;                  /* Current hash index   */

      vertend = grafptr->edgetax[edgeold];
      for (hashnum = (vertend * GRAPHGEOMHABOHASHPRIME) & hashmsk; ; hashnum = (hashnum + 1) & hashmsk) {
        if (hashtab[hashnum].vertnum != vertnum) { /* If edge not found */
          hashtab[hashnum].vertnum = vertnum;
          hashtab[hashnum].vertend = vertend;
          grafptr->edgetax[edgenum ++] = vertend;
          break;
        }
        if (hashtab[hashnum].vertend == vertend)  /* Do not add duplicate edges */
          break;
      }
    }
    if ((edgenum - grafptr->verttax[vertnum]) > degrmax) /* Set real maximum degree */
      degrmax = edgenum - grafptr->verttax[vertnum];
  }
  grafptr->verttax[vertnum] = edgenum;            /* Set end of vertex array */

  grafptr->edgenbr = edgenum - 1;
  grafptr->edlosum = grafptr->edgenbr;
  grafptr->degrmax = degrmax;

  memFree (hashtab);

  grafptr->edgetax = ((Gnum *) memRealloc (grafptr->edgetax + 1, grafptr->edgenbr * sizeof (Gnum))) - 1;

#ifdef SCOTCH_DEBUG_GRAPH2
  if (graphCheck (grafptr) != 0) {                /* Check graph consistency */
    errorPrint ("graphGeomLoadHabo: internal error");
    graphFree  (grafptr);
    return     (1);
  }
#endif /* SCOTCH_DEBUG_GRAPH2 */

  return (0);
}

/* This routine reads a Fortran format structure
** and returns the size of the integers to read.
*/

static
int
graphGeomLoadHaboFormat (
GraphGeomHaboLine * restrict const  lineptr,      /* Line format to fill */
const char * const                  dataptr)      /* Format string       */
{
  const char * restrict   charptr;
  int                     number;                 /* Number to read */

  lineptr->strtnbr =
  lineptr->datanbr =
  lineptr->datalen = 0;
  for (charptr = dataptr; ; charptr ++) {         /* Skip to first '(' */
    if (*charptr == '(')
      break;
    if (*charptr == '\0')                         /* Error if end of string */
      return (1);
  }

  number = 0;                                     /* Read number */
  for (charptr ++ ; ; charptr ++) {
    if (*charptr == '\0')                         /* Error if end of string */
      return (1);
    if (! isdigit ((int) (unsigned char) *charptr))
      break;
    number = number * 10 + *charptr - '0';
  }

  if ((*charptr == 'x') || (*charptr == 'X')) {   /* If dummy characters at beginning of line */
    lineptr->strtnbr = number;
    for (charptr ++; ; charptr ++) {              /* Skip to first ','      */
      if (*charptr == '\0')                       /* Error if end of string */
        return (1);
      if (*charptr == ',')
        break;
    }

    number = 0;                                   /* Read number */
    for (charptr ++; *charptr != '\0'; charptr ++) {
      if (*charptr == '\0')                       /* Error if end of string */
        return (1);
      if (! isdigit (*charptr < '0'))
        break;
      number = number * 10 + *charptr - '0';
    }
  }

  if ((*charptr != 'I') && (*charptr != 'i'))     /* If not integer specification */
    return (1);

  lineptr->datanbr = number;

  number = 0;                                     /* Read size of integer */
  for (charptr ++; ; charptr ++) {
    if (*charptr == '\0')                         /* Error if end of string */
      return (1);
    if (! isdigit ((int) (unsigned char) *charptr))
      break;
    number = number * 10 + *charptr - '0';
  }

  if (number == 0)
    return (1);

  lineptr->datalen = number;

  return (0);
}
