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
/**   NAME       : common_file_uncompress.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles compressed streams  **/
/**                for uncompression.                      **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 11 mar 2008     **/
/**                                 to   : 15 may 2008     **/
/**                # Version 5.1  : from : 27 jun 2010     **/
/**                                 to     27 jun 2010     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define COMMON_FILE
#define COMMON_FILE_UNCOMPRESS

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
#ifdef COMMON_FILE_COMPRESS_LZMA
#include "lzmadec.h"                              /* TODO: Temporary interface */
#endif /* COMMON_FILE_COMPRESS_LZMA */

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
#ifdef COMMON_FILE_COMPRESS_LZMA
                                          { ".lzma", FILECOMPRESSTYPELZMA    },
#else /* COMMON_FILE_COMPRESS_LZMA */
                                          { ".lzma", FILECOMPRESSTYPENOTIMPL },
#endif /* COMMON_FILE_COMPRESS_LZMA */
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
fileUncompressType (
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

/* This routine creates a thread to uncompress the
** given stream according to the given (un)compression
** algorithm.
** If threads are available, uncompression will be
** performed by an auxiliary thread. Else, a child process
** will be fork()'ed, and after completion this process
** will remain a zombie until the main process terminates.
** It returns:
** - !NULL  : stream holding uncompressed data.
** - NULL   : on error.
*/

static
void *                                            /* (void *) to comply to the Posix pthread API */
fileUncompress2 (
FileCompressData * const  dataptr)
{
  switch (dataptr->typeval) {
#ifdef COMMON_FILE_COMPRESS_BZ2
    case FILECOMPRESSTYPEBZ2 :
      fileUncompressBz2 (dataptr);
      break;
#endif /* COMMON_FILE_COMPRESS_BZ2 */
#ifdef COMMON_FILE_COMPRESS_GZ
    case FILECOMPRESSTYPEGZ :
      fileUncompressGz (dataptr);
      break;
#endif /* COMMON_FILE_COMPRESS_GZ */
#ifdef COMMON_FILE_COMPRESS_LZMA
    case FILECOMPRESSTYPELZMA :
      fileUncompressLzma (dataptr);
      break;
#endif /* COMMON_FILE_COMPRESS_LZMA */
    default :
      errorPrint ("fileUncompress2: method not implemented");
  }

  close   (dataptr->innerfd);                     /* Close writer's end */
  memFree (dataptr);                              /* Free buffers       */

  return ((void *) 0);                            /* Don't care anyway */
}

FILE *
fileUncompress (
FILE * const                stream,               /*+ Compressed stream         +*/
const int                   typeval)              /*+ (Un)compression algorithm +*/
{
  int                 filetab[2];
  FILE *              readptr;
  FileCompressData *  dataptr;
#ifdef COMMON_PTHREAD
  pthread_t           thrdval;
#endif /* COMMON_PTHREAD */

  if (typeval <= FILECOMPRESSTYPENONE)            /* If uncompressed stream, return original stream pointer */
    return (stream);

  if (pipe (filetab) != 0) {
    errorPrint ("fileUncompress: cannot create pipe");
    return (NULL);
  }

  if ((readptr = fdopen (filetab[0], "r")) == NULL) {
    errorPrint ("fileUncompress: cannot create stream");
    close  (filetab[0]);
    close  (filetab[1]);
    return (NULL);
  }

  if ((dataptr = memAlloc (sizeof (FileCompressData) + FILECOMPRESSDATASIZE)) == NULL) {
    errorPrint ("fileUncompress: out of memory");
    fclose (readptr);
    close  (filetab[1]);
    return (NULL);
  }

  dataptr->typeval     = typeval;                 /* Fill structure to be passed to uncompression thread/process */
  dataptr->innerfd     = filetab[1];
  dataptr->outerstream = stream;

#ifdef COMMON_PTHREAD
  if (pthread_create (&thrdval, NULL, (void * (*) (void *)) fileUncompress2, (void *) dataptr) != 0) { /* If could not create thread */
    errorPrint ("fileUncompress: cannot create thread");
    memFree (dataptr);
    fclose  (readptr);
    close   (filetab[1]);
    return  (NULL);
  }
  pthread_detach (thrdval);                       /* Detach thread so that it will end up gracefully by itself */
#else /* COMMON_PTHREAD */
  switch (fork ()) {
    case -1 :                                     /* Error */
      errorPrint ("fileUncompress: cannot create child process");
      memFree (dataptr);
      fclose  (readptr);
      close   (filetab[1]);
      return  (NULL);
    case 0 :                                      /* We are the son process    */
      fclose (readptr);                           /* Close reader pipe stream  */
      fileUncompress2 (dataptr);                  /* Perform uncompression     */
      exit (0);                                   /* Exit gracefully           */
    default :                                     /* We are the father process */
      close (filetab[1]);                         /* Close the writer pipe end */
  }
#endif /* COMMON_PTHREAD */

  return (readptr);
}

/* This routine uncompresses a stream compressed in the
** bzip2 format.
** It returns:
** - void  : in all cases. Uncompression stops immediately
**           in case of error.
*/

#ifdef COMMON_FILE_COMPRESS_BZ2
static
void
fileUncompressBz2 (
FileCompressData * const  dataptr)
{
  BZFILE *              bzfile;
  int                   bzsize;
  int                   bzerror;

  if (FILECOMPRESSDATASIZE < (BZ_MAX_UNUSED)) {
    errorPrint ("fileUncompressBz2: cannot start decompression (1)");
    return;
  }
  if ((bzfile = BZ2_bzReadOpen (&bzerror, dataptr->outerstream, 0, 0, NULL, 0)) == NULL) {
    errorPrint ("fileUncompressBz2: cannot start decompression (2)");
    BZ2_bzReadClose (&bzerror, bzfile);
    return;
  }

  while ((bzsize = BZ2_bzRead (&bzerror, bzfile, &dataptr->datatab, FILECOMPRESSDATASIZE), bzerror) >= BZ_OK) { /* If BZ_OK or BZ_STREAM_END */
    if (write (dataptr->innerfd, &dataptr->datatab, bzsize) != bzsize) {
      errorPrint ("fileUncompressBz2: cannot write");
      bzerror = BZ_STREAM_END;                    /* Avoid other error message */
      break;
    }
    if (bzerror == BZ_STREAM_END) {               /* If end of compressed stream */
      void *                bzunusptr;
      int                   bzunusnbr;

      BZ2_bzReadGetUnused (&bzerror, bzfile, &bzunusptr, &bzunusnbr); /* Get remaining chars in stream   */
      if ((bzunusnbr == 0) && (feof (dataptr->outerstream) != 0)) { /* If end of uncompressed stream too */
        bzerror = BZ_STREAM_END;
        break;
      }
      memMov (&dataptr->datatab, bzunusptr, bzunusnbr);
      BZ2_bzReadClose (&bzerror, bzfile);
      if ((bzfile = BZ2_bzReadOpen (&bzerror, dataptr->outerstream, 0, 0, &dataptr->datatab, bzunusnbr)) == NULL) {
        errorPrint ("fileUncompressBz2: cannot start decompression (3)");
        bzerror = BZ_STREAM_END;
        break;
      }
    }
  }
  if (bzerror != BZ_STREAM_END)
    errorPrint ("fileUncompressBz2: cannot read");

  BZ2_bzReadClose (&bzerror, bzfile);
  fclose (dataptr->outerstream);                  /* Do as zlib does */
}
#endif /* COMMON_FILE_COMPRESS_BZ2 */

/* This routine uncompresses a stream compressed in the
** gzip format.
** It returns:
** - void  : in all cases. Uncompression stops immediately
**           in case of error.
*/

#ifdef COMMON_FILE_COMPRESS_GZ
static
void
fileUncompressGz (
FileCompressData * const  dataptr)
{
  gzFile                gzfile;
  int                   gzsize;

  if ((gzfile = gzdopen (fileno (dataptr->outerstream), "rb")) == NULL) {
    errorPrint ("fileUncompressGz: cannot start decompression");
    return;
  }

  while ((gzsize = gzread (gzfile, &dataptr->datatab, FILECOMPRESSDATASIZE)) > 0) {
    if (write (dataptr->innerfd, &dataptr->datatab, gzsize) != gzsize) {
      errorPrint ("fileUncompressGz: cannot write");
      break;
    }
  }
  if (gzsize < 0)
    errorPrint ("fileUncompressGz: cannot read");

  gzclose (gzfile);
}
#endif /* COMMON_FILE_COMPRESS_GZ */

/* This routine uncompresses a stream compressed in the
** lzma format.
** It returns:
** - void  : in all cases. Uncompression stops immediately
**           in case of error.
*/

#ifdef COMMON_FILE_COMPRESS_LZMA
static
void
fileUncompressLzma (
FileCompressData * const  dataptr)
{
  lzmadec_FILE *        lzmafile;
  ssize_t               lzmasize;

  if ((lzmafile = lzmadec_dopen (fileno (dataptr->outerstream))) == NULL) {
    errorPrint ("fileUncompressLzma: cannot start decompression");
    return;
  }

  while ((lzmasize = lzmadec_read (lzmafile, (void *) &dataptr->datatab, FILECOMPRESSDATASIZE)) > 0) {
    if (write (dataptr->innerfd, &dataptr->datatab, lzmasize) != lzmasize) {
      errorPrint ("fileUncompressLzma: cannot write");
      break;
    }
  }
  if (lzmasize < 0)
    errorPrint ("fileUncompressLzma: cannot read");

  lzmadec_close (lzmafile);
}
#endif /* COMMON_FILE_COMPRESS_LZMA */
