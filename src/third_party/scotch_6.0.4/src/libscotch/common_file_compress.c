/* Copyright 2008,2010 ENSEIRB, INRIA & CNRS
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
/**   NAME       : common_file_compress.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles compressed streams  **/
/**                for compression.                        **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 13 mar 2008     **/
/**                                 to   : 15 may 2008     **/
/**                # Version 5.1  : from : 27 jun 2010     **/
/**                                 to     27 jun 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON_FILE
#define COMMON_FILE_COMPRESS

#ifndef COMMON_NOMODULE
#include "module.h"
#endif /* COMMON_NOMODULE */
#include "common.h"
#include "common_file.h"
#include "common_file_compress.h"
#ifdef COMMON_FILE_COMPRESS_BZ2
#include "bzlib.h"
#endif /* COMMON_FILE_COMPRESS_BZ2 */
#ifdef COMMON_FILE_COMPRESS_GZ
#include "zlib.h"
#endif /* COMMON_FILE_COMPRESS_GZ */

/*
**  The static definitions.
*/

static FileCompressTab      filetab[] = {
#ifdef COMMON_FILE_COMPRESS_BZ2
                                          { ".bz2",  FILECOMPRESSTYPEBZ2,    },
#else /* COMMON_FILE_COMPRESS_BZ2 */
                                          { ".bz2",  FILECOMPRESSTYPENOTIMPL },
#endif /* COMMON_FILE_COMPRESS_BZ */
#ifdef COMMON_FILE_COMPRESS_GZ
                                          { ".gz",   FILECOMPRESSTYPEGZ,     },
#else /* COMMON_FILE_COMPRESS_GZ */
                                          { ".gz",   FILECOMPRESSTYPENOTIMPL },
#endif /* COMMON_FILE_COMPRESS_GZ */
/* #ifdef COMMON_FILE_COMPRESS_LZMA */ /* TODO: Current lzmadec library does not offer compression to date */
/*                                          { ".lzma", FILECOMPRESSTYPELZMA    }, */
/* #else COMMON_FILE_COMPRESS_LZMA */
                                          { ".lzma", FILECOMPRESSTYPENOTIMPL },
/* #endif /\* COMMON_FILE_COMPRESS_LZMA *\/ */
                                          { NULL,    FILECOMPRESSTYPENOTIMPL } };

/*********************************/
/*                               */
/* Basic routines for filenames. */
/*                               */
/*********************************/

/* This routine searches the given file name
** for relevant extensions and returns the
** corresponding code if it is the case.
** It returns:
** - FILECOMPRESSTYPENONE     : no recognized file extension.
** - FILECOMPRESSTYPENOTIMPL  : compression algorithm not implemented.
** - FILECOMPRESSTYPExxxx     : implemented compression algorithm.
*/

int
fileCompressType (
const char * const          nameptr)              /*+ Name string +*/
{
  int                 namelen;
  int                 i;

  namelen = strlen (nameptr);
  for (i = 0; filetab[i].name != NULL; i ++) {
    int                 extnlen;                  /* Name of extension string */

    extnlen = strlen (filetab[i].name);
    if ((namelen >= extnlen) && (strncmp (filetab[i].name, nameptr + (namelen - extnlen), extnlen) == 0))
      return (filetab[i].type);
  }

  return (FILECOMPRESSTYPENONE);
}

/* This routine creates a thread to compress the
** given stream according to the given compression
** algorithm.
** If threads are available, compression will be
** performed by an auxiliary thread. Else, a child process
** will be fork()'ed, and after completion this process
** will remain a zombie until the main process terminates.
** It returns:
** - !NULL  : stream holding compressed data.
** - NULL   : on error.
*/

static
void *                                            /* (void *) to comply to the Posix pthread API */
fileCompress2 (
FileCompressData * const  dataptr)
{
  switch (dataptr->typeval) {
#ifdef COMMON_FILE_COMPRESS_BZ2
    case FILECOMPRESSTYPEBZ2 :
      fileCompressBz2 (dataptr);
      break;
#endif /* COMMON_FILE_COMPRESS_BZ2 */
#ifdef COMMON_FILE_COMPRESS_GZ
    case FILECOMPRESSTYPEGZ :
      fileCompressGz (dataptr);
      break;
#endif /* COMMON_FILE_COMPRESS_GZ */
/* #ifdef COMMON_FILE_COMPRESS_LZMA /\* TODO: Current lzmadec library does not offer compression to date *\/ */
/*     case FILECOMPRESSTYPELZMA : */
/*       fileCompressLzma (dataptr); */
/*       break; */
/* #endif /\* COMMON_FILE_COMPRESS_LZMA *\/ */
    default :
      errorPrint ("fileCompress2: method not implemented");
  }

  close   (dataptr->innerfd);                     /* Close writer's end */
  memFree (dataptr);                              /* Free buffers       */

  return ((void *) 0);                            /* Don't care anyway */
}

FILE *
fileCompress (
FILE * const                stream,               /*+ Uncompressed stream       +*/
const int                   typeval)              /*+ (Un)compression algorithm +*/
{
  int                 filetab[2];
  FILE *              writptr;
  FileCompressData *  dataptr;
#ifdef COMMON_PTHREAD
  pthread_t           thrdval;
#endif /* COMMON_PTHREAD */

  if (typeval <= FILECOMPRESSTYPENONE)            /* If uncompressed stream, return original stream pointer */
    return (stream);

  if (pipe (filetab) != 0) {
    errorPrint ("fileCompress: cannot create pipe");
    return (NULL);
  }

  if ((writptr = fdopen (filetab[1], "w")) == NULL) {
    errorPrint ("fileCompress: cannot create stream");
    close  (filetab[0]);
    close  (filetab[1]);
    return (NULL);
  }

  if ((dataptr = memAlloc (sizeof (FileCompressData) + FILECOMPRESSDATASIZE)) == NULL) {
    errorPrint ("fileCompress: out of memory");
    close  (filetab[0]);
    fclose (writptr);
    return (NULL);
  }

  dataptr->typeval     = typeval;                 /* Fill structure to be passed to compression thread/process */
  dataptr->innerfd     = filetab[0];
  dataptr->outerstream = stream;                  /* Stream to write to */

#ifdef COMMON_PTHREAD
  if (pthread_create (&thrdval, NULL, (void * (*) (void *)) fileCompress2, (void *) dataptr) != 0) { /* If could not create thread */
    errorPrint ("fileCompress: cannot create thread");
    memFree (dataptr);
    close   (filetab[0]);
    fclose  (writptr);
    return  (NULL);
  }
#else /* COMMON_PTHREAD */
  switch (fork ()) {
    case -1 :                                     /* Error */
      errorPrint ("fileCompress: cannot create child process");
      memFree (dataptr);
      close   (filetab[0]);
      fclose  (writptr);
      return  (NULL);
    case 0 :                                      /* We are the son process    */
      fclose (writptr);                           /* Close writer pipe stream  */
      fileCompress2 (dataptr);                    /* Perform compression       */
      exit (0);                                   /* Exit gracefully           */
    default :                                     /* We are the father process */
      close (filetab[0]);                         /* Close the reader pipe end */
  }
#endif /* COMMON_PTHREAD */

  return (writptr);
}

/* This routine compresses a stream compressed in the
** gzip format.
** It returns:
** - void  : in all cases. Compression stops immediately
**           in case of error.
*/

#ifdef COMMON_FILE_COMPRESS_BZ2
static
void
fileCompressBz2 (
FileCompressData * const  dataptr)
{
  BZFILE *              bzfile;
  int                   bzsize;
  int                   bzerror;

  if ((bzfile = BZ2_bzWriteOpen (&bzerror, dataptr->outerstream, 9, 0, 0)) == NULL) {
    errorPrint ("fileCompressBz2: cannot start compression");
    BZ2_bzWriteClose (&bzerror, bzfile, 1, NULL, NULL);
    return;
  }

  while ((bzsize = read (dataptr->innerfd, &dataptr->datatab, FILECOMPRESSDATASIZE)) > 0) { /* Read from pipe */
    BZ2_bzWrite (&bzerror, bzfile, &dataptr->datatab, bzsize);
    if (bzerror != BZ_OK) {
      errorPrint ("fileCompressBz2: cannot write");
      break;
    }
  }
  if (bzsize < 0) {
    errorPrint ("fileCompressBz2: cannot read");
    bzerror = BZ_STREAM_END;                      /* Will set abandon flag to 1 in BZ2_bzWriteClose */
  }

  BZ2_bzWriteClose (&bzerror, bzfile, (bzerror != BZ_OK) ? 1 : 0, NULL, NULL);
  fclose (dataptr->outerstream);                  /* Do as zlib does */
}
#endif /* COMMON_FILE_COMPRESS_BZ2 */

/* This routine compresses a stream compressed in the
** gzip format.
** It returns:
** - void  : in all cases. Compression stops immediately
**           in case of error.
*/

#ifdef COMMON_FILE_COMPRESS_GZ
static
void
fileCompressGz (
FileCompressData * const  dataptr)
{
  gzFile                gzfile;
  int                   gzsize;

  if ((gzfile = gzdopen (fileno (dataptr->outerstream), "wb")) == NULL) {
    errorPrint ("fileCompressGz: cannot start compression");
    return;
  }
  gzsetparams (gzfile, 9, Z_DEFAULT_STRATEGY);    /* Maximum compression */

  while ((gzsize = read (dataptr->innerfd, &dataptr->datatab, FILECOMPRESSDATASIZE)) > 0) { /* Read from pipe */
    if (gzwrite (gzfile, &dataptr->datatab, gzsize) != gzsize) {
      errorPrint ("fileCompressGz: cannot write");
      break;
    }
  }
  if (gzsize < 0)
    errorPrint ("fileCompressGz: cannot read");

  gzclose (gzfile);
}
#endif /* COMMON_FILE_COMPRESS_GZ */
