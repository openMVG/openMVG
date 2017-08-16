/* Copyright 2007,2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : common_file.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles files and file      **/
/**                names.                                  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 21 may 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 27 jun 2010     **/
/**                                 to     27 jun 2010     **/
/**                # Version 6.0  : from : 10 nov 2014     **/
/**                                 to     12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON_FILE

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"
#include "common_file.h"

/*********************************/
/*                               */
/* Basic routines for filenames. */
/*                               */
/*********************************/

/* This routine expands distributed filenames
** according to process numbers and the root
** process number.
** It returns:
** - 0  : on success.
** - !0 : on error.
*/

int
fileNameDistExpand (
char ** const               nameptr,              /*+ Pointer to name string pointer +*/
const int                   procnbr,              /*+ Number of processes            +*/
const int                   procnum,              /*+ Number of current process      +*/
const int                   protnum)              /*+ Root process number            +*/
{
  int                 namemax;
  int                 namenum;
  char *              naexptr;
  int                 naexmax;
  int                 naexnum;
  int                 flagval;                    /* Flag set if expansion took place */

  namemax = strlen (*nameptr);
  naexmax = namemax + FILENAMEDISTEXPANDNBR * 2;

  if ((naexptr = memAlloc ((naexmax + 1) * sizeof (char))) == NULL) /* "+ 1" for terminating character */
    return (1);

#ifdef COMMON_DEBUG
  sprintf (naexptr, FILENAMEDISTEXPANDSTR, procnbr); /* TRICK: Test if FILENAMEDISTEXPANDNBR is a size large enough */
  if (atoi (naexptr) != procnbr) {
    errorPrint ("fileNameDistExpand: undersized integer string size");
    memFree    (naexptr);
    return     (1);
  }
#endif /* COMMON_DEBUG */

  for (namenum = naexnum = flagval = 0; namenum < namemax; ) {
    char                charval;
    int                 dataval = 0;
    int                 datasiz;

    charval = (*nameptr)[namenum ++];             /* Get current characted                */
    datasiz = 1;                                  /* Assume individual expanded character */
    if (charval == '%') {
      char                chnxval;                /* Value of next char */

      switch (chnxval = (*nameptr)[namenum ++]) {
        case 'p' :                                /* "%p" translates into number of processes */
          flagval = 1;
          datasiz = FILENAMEDISTEXPANDNBR;
          dataval = procnbr;
          break;
        case 'r' :                                /* "%r" translates into process rank */
          flagval = 1;
          datasiz = FILENAMEDISTEXPANDNBR;
          dataval = procnum;
          break;
        case '-' :                                /* "%-" translates into nothing but indicates distributed file */
          datasiz = 0;
          flagval = 1;
          break;
        case '%' :                                /* "%%" translates into '%' */
          break;
        default :
          charval = chnxval;                      /* Unrecognized character */
      }
    }
    if (datasiz > 0) {
      if ((naexnum + datasiz) > naexmax) {
        char *              nanwptr;

        naexmax += datasiz + FILENAMEDISTEXPANDNBR;
        if ((nanwptr = memRealloc (naexptr, (naexmax + 1) * sizeof (char))) == NULL) { /* "+ 1" for terminating character */
          memFree (naexptr);
          return  (1);
        }
        naexptr = nanwptr;
      }
      if (datasiz == 1)
        naexptr[naexnum ++] = charval;
      else {
        sprintf (&naexptr[naexnum], FILENAMEDISTEXPANDSTR, dataval); /* TRICK: Change format string if FILENAMEDISTEXPANDNBR changes */
        naexptr[naexnum + FILENAMEDISTEXPANDNBR] = ' ';
        naexnum = strchr (&naexptr[naexnum], ' ') - naexptr;
      }
    }
  }
  naexptr[naexnum] = '\0';                        /* Set terminating character as there is always room for it */

  if ((flagval != 0) || (procnum == protnum))     /* If file name is a distributed one or we are the root process */
    *nameptr = naexptr;                           /* Replace new string by old one                                */
  else {
    memFree (naexptr);                            /* No need for the expanded string */
    *nameptr = NULL;
  }

  return (0);
}

/* This routine initializes a block of
** file descriptor structures.
** It returns:
** - void  : in all cases.
*/

void
fileBlockInit (
File * const                filetab,
const int                   filenbr)
{
  int                 i;

  for (i = 0; i < filenbr; i ++) {                /* For all file names     */
    filetab[i].nameptr = "-";                     /* Assume standard stream */
    filetab[i].fileptr = (filetab[i].modeptr[0] == 'r') ? stdin : stdout;
    filetab[i].dataptr = NULL;                    /* Assume nothing to free */
  }
}

/* This routine opens a block of file
** descriptor structures.
** It returns:
** - 0  : on success.
** - !0 : on error.
*/

int
fileBlockOpen (
File * const                filetab,
const int                   filenbr)
{
  int                 i, j;

  for (i = 0; i < filenbr; i ++) {                /* For all file names             */
    if (filetab[i].fileptr == NULL)               /* If unwanted stream, do nothing */
      continue;

    for (j = 0; j < i; j ++) {
      if ((filetab[i].modeptr[0] == filetab[j].modeptr[0]) && /* If very same name with same opening mode */
          (filetab[j].nameptr != NULL)                     &&
          (strcmp (filetab[i].nameptr, filetab[j].nameptr) == 0)) {
        filetab[i].fileptr = filetab[j].fileptr;  /* Share pointer to already processed stream */
        filetab[i].nameptr = NULL;                /* Do not close this stream multiple times   */
        break;
      }
    }
    if (j == i) {                                 /* If original stream                */
      int                 compval;                /* Compression type                  */
      FILE *              compptr;                /* Processed ((un)compressed) stream */

      if (filetab[i].nameptr[0] != '-') {         /* If not standard stream, open it */
        if ((filetab[i].fileptr = fopen (filetab[i].nameptr, filetab[i].modeptr)) == NULL) { /* Open the file */
          errorPrint ("fileBlockOpen: cannot open file (%d)", i);
          return     (1);
        }
      }
      compval = (filetab[i].modeptr[0] == 'r') ? fileUncompressType (filetab[i].nameptr) : fileCompressType (filetab[i].nameptr);
      if (compval < 0) {
        errorPrint ("fileBlockOpen: (un)compression type not implemented");
        return     (1);
      }
      compptr = (filetab[i].modeptr[0] == 'r') ? fileUncompress (filetab[i].fileptr, compval) : fileCompress (filetab[i].fileptr, compval);
      if (compptr == NULL) {
        errorPrint ("fileBlockOpen: cannot create (un)compression subprocess");
        return     (1);
      }
      filetab[i].fileptr = compptr;               /* Use processed stream instead of original stream */
    }
  }

  return (0);
}

/* This routine opens a block of eventually
** distributed file descriptor structures.
** It returns:
** - 0  : on success.
** - !0 : on error.
*/

int
fileBlockOpenDist (
File * const                filetab,
const int                   filenbr,
const int                   procglbnbr,
const int                   proclocnum,
const int                   protglbnum)
{
  int                 i, j;

  for (i = 0; i < filenbr; i ++) {                /* For all file names             */
    if (filetab[i].fileptr == NULL) {             /* If unwanted stream, do nothing */
      filetab[i].nameptr = NULL;                  /* Do nothing on closing, too     */
      continue;
    }

    if (fileNameDistExpand (&filetab[i].nameptr, procglbnbr, proclocnum, protglbnum) != 0) { /* If cannot allocate new name */
      errorPrint ("fileBlockOpenDist: cannot create file name (%d)", i);
      return     (1);
    }
    if (filetab[i].nameptr == NULL) {             /* If inexisting stream because not root process and centralized stream */
      filetab[i].fileptr = NULL;
      continue;
    }
    filetab[i].dataptr = filetab[i].nameptr;      /* From now on, we will have to free it anyway */

    for (j = 0; j < i; j ++) {
      if ((filetab[i].modeptr[0] == filetab[j].modeptr[0]) && /* If very same name with same opening mode */
          (filetab[j].nameptr != NULL)                     &&
          (strcmp (filetab[i].nameptr, filetab[j].nameptr) == 0)) {
        filetab[i].fileptr = filetab[j].fileptr;  /* Share pointer to already processed stream */
        filetab[i].nameptr = NULL;                /* Do not close this stream multiple times   */
        break;
      }
    }
    if (j == i) {                                 /* If original stream                */
      int                 compval;                /* Compression type                  */
      FILE *              compptr;                /* Processed ((un)compressed) stream */

      if (filetab[i].nameptr[0] != '-') {         /* If not standard stream, open it */
        if ((filetab[i].fileptr = fopen (filetab[i].nameptr, filetab[i].modeptr)) == NULL) { /* Open the file */
          errorPrint ("fileBlockOpenDist: cannot open file (%d)", i);
          return     (1);
        }
      }
      compval = (filetab[i].modeptr[0] == 'r') ? fileUncompressType (filetab[i].nameptr) : fileCompressType (filetab[i].nameptr);
      if (compval < 0) {
        errorPrint ("fileBlockOpenDist: (un)compression type not implemented");
        return     (1);
      }
      compptr = (filetab[i].modeptr[0] == 'r') ? fileUncompress (filetab[i].fileptr, compval) : fileCompress (filetab[i].fileptr, compval);
      if (compptr == NULL) {
        errorPrint ("fileBlockOpenDist: cannot create (un)compression subprocess");
        return     (1);
      }
      filetab[i].fileptr = compptr;               /* Use processed stream instead of original stream */
    }
  }

  return (0);
}

/* This routine opens a block of file
** descriptor structures.
** It returns:
** - 0  : on success.
** - !0 : on error.
*/

void
fileBlockClose (
File * const                filetab,
const int                   filenbr)
{
  int                i;

  for (i = 0; i < filenbr; i ++) {                /* For all file names   */
    if ((filetab[i].fileptr != NULL) &&           /* If stream exists     */
        (filetab[i].nameptr != NULL) &&           /* And it has a name    */
        (filetab[i].nameptr[0] != '-'))           /* Which is not default */
      fclose (filetab[i].fileptr);                /* Close the stream     */

    if (filetab[i].dataptr != NULL)               /* If file name data to free */
      memFree (filetab[i].dataptr);
  }
}
