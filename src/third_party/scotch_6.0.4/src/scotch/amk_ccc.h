/* Copyright 2004,2007,2008,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : amk_ccc.h                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the distance map for CCC        **/
/**                graphs, to be used to build the         **/
/**                architecture description files for      **/
/**                these graphs.                           **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 24 apr 1994     **/
/**                                 to   : 24 apr 1994     **/
/**                # Version 2.0  : from : 13 jul 1994     **/
/**                                 to   : 18 jul 1994     **/
/**                # Version 3.0  : from : 18 sep 1995     **/
/**                                 to   : 18 sep 1995     **/
/**                # Version 3.2  : from : 07 may 1997     **/
/**                                 to   : 07 may 1997     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 02 oct 1998     **/
/**                # Version 6.0  : from : 12 nov 2014     **/
/**                                 to   : 12 nov 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

/*+ File name aliases. +*/

#define C_FILENBR                   1             /* Number of files in list */

#define C_filenamearcout            fileBlockName (C_fileTab, 0) /* Architecture output file name */

#define C_filepntrarcout            fileBlockFile (C_fileTab, 0) /* Architecture output file */

/*
**  The type and structure definitions.
*/

/*+ This structure defines a CCC vertex. +*/

typedef struct C_Vertex_ {
  SCOTCH_Num                lvl;                  /*+ Vertex level    +*/
  SCOTCH_Num                pos;                  /*+ Vertex position +*/
} C_Vertex;

/*+ This structure defines a vertex distance information. +*/

typedef struct C_VertDist_ {
  int                       queued;               /*+ Flag set if vertex queued  +*/
  SCOTCH_Num                dist;                 /*+ Distance to initial vertex +*/
} C_VertDist;

/*+ This is a neighbor queue element. +*/

typedef struct C_QueueElem_ {
  C_Vertex                  vert;                 /*+ Vertex number    +*/
  SCOTCH_Num                dist;                 /*+ Distance reached +*/
} C_QueueElem;

/*+ This is the distance queue. +*/

typedef struct C_Queue_ {
  C_QueueElem *             tab;                  /*+ Pointer to queue data    +*/
  SCOTCH_Num                min;                  /*+ Pointer to first element +*/
  SCOTCH_Num                max;                  /*+ Pointer to last element  +*/
} C_Queue;

/*
**  The macro definitions.
*/

#define C_vertLabl(v)               (((v)->lvl << ccdim) | (v)->pos)

#define C_queueInit(q,n)            ((((q)->tab = (C_QueueElem *) memAlloc ((n) * sizeof (C_QueueElem))) == NULL) ? 1 : 0)
#define C_queueExit(q)              memFree ((q)->tab)
#define C_queueFlush(q)             (q)->min = \
                                    (q)->max = 0
#define C_queuePut(q,v,d)           ((q)->tab[(q)->max].vert    = *(v),  \
                                     (q)->tab[(q)->max ++].dist = (d))
#define C_queueGet(q,v,d)           (((q)->min < (q)->max) ? (*(v) = (q)->tab[(q)->min].vert,    \
                                                              *(d) = (q)->tab[(q)->min ++].dist, \
                                                              1)                                 \
                                                           : 0)

#define C_distaRoot(v)              (C_queueFlush (&C_distaQueue),         \
                                     C_queuePut   (&C_distaQueue, (v), 0), \
                                     C_distaTab[C_vertLabl (v)].queued = 1)
#define C_distaGet(v,d)             (C_queueGet (&C_distaQueue, (v), (d)))
#define C_distaPut(v,d)             ((C_distaTab[C_vertLabl (v)].queued == 0) \
                                     ? C_queuePut (&C_distaQueue, (v), d),    \
                                       C_distaTab[C_vertLabl (v)].queued = 1  \
                                     : 0)
