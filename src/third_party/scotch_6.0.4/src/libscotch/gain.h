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
/**   NAME       : gain.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the definitions of **/
/**                the generic gain tables.                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 26 oct 1996     **/
/**                                 to     17 nov 1997     **/
/**                # Version 0.1  : from : 10 may 1999     **/
/**                                 to     18 mar 2005     **/
/**                # Version 5.0  : from : 24 mar 2008     **/
/**                                 to     01 jun 2008     **/
/**                                                        **/
/**   NOTES      : # Most of the contents of this module   **/
/**                  comes from "map_b_fm" of the SCOTCH   **/
/**                  project.                              **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

#define GAINMAX                  ((INT) (((UINT) 1 << ((sizeof (INT) << 3) - 1)) - 2))

#define GAIN_LINMAX              1024

/*
**  The type and structure definitions.
*/

/* The gain link data structure. This must be the
   first item of objects that are linked into gain
   tables.                                         */

typedef struct GainLink_ {
  struct GainLink_ *        next;                 /*+ Pointer to next element: FIRST +*/
  struct GainLink_ *        prev;                 /*+ Pointer to previous element    +*/
  struct GainEntr_ *        tabl;                 /*+ Index into the gain table      +*/
} GainLink;

/* Gain table entry structure. */

typedef struct GainEntr_ {
  GainLink *                next;                 /*+ Pointer to first element: FIRST +*/
} GainEntr;

/* The gain table structure, built from table entries.
   For trick reasons, the pointer to the first entry
   must be the first field of the structure.           */

typedef struct GainTabl_ {
  void                   (* tablAdd)  (struct GainTabl_ * const, GainLink * const, const INT); /*+ Add method +*/
  INT                       subbits;              /*+ Number of subbits                      +*/
  INT                       submask;              /*+ Subbit mask                            +*/
  INT                       totsize;              /*+ Total table size                       +*/
  GainEntr *                tmin;                 /*+ Non-empty entry of minimum gain        +*/
  GainEntr *                tmax;                 /*+ Non-empty entry of maximum gain        +*/
  GainEntr *                tend;                 /*+ Point after last valid gain entry      +*/
  GainEntr *                tabl;                 /*+ Gain table structure is.. [SIZE - ADJ] +*/
  GainEntr                  tabk[1];              /*+ Split in two for relative access [ADJ] +*/
} GainTabl;

/*
**  The function prototypes.
*/

#ifndef GAIN
#define static
#endif

GainTabl *                  gainTablInit        (const INT, const INT);
void                        gainTablExit        (GainTabl * const);
void                        gainTablFree        (GainTabl * const);
void                        gainTablAddLin      (GainTabl * const, GainLink * const, const INT);
void                        gainTablAddLog      (GainTabl * const, GainLink * const, const INT);
void                        gainTablDel         (GainTabl * const, GainLink * const);
GainLink *                  gainTablFrst        (GainTabl * const);
GainLink *                  gainTablNext        (GainTabl * const, const GainLink * const);
#ifdef SCOTCH_DEBUG_GAIN3
int                         gainTablCheck       (GainEntr * const);
static int                  gainTablCheck2      (GainEntr * const, GainLink * const);
#endif /* SCOTCH_DEBUG_GAIN3 */

#undef static

/*
**  The marco definitions.
*/

#define gainTablEmpty(tabl)         ((tabl)->tmin == (tabl)->tend)
#define gainTablAdd(tabl,link,gain) ((tabl)->tablAdd  ((tabl), (link), (gain)))
#if ((! defined GAIN) && (! defined SCOTCH_DEBUG_GAIN1))
#define gainTablDel(tabl,link)      (((GainLink *) (link))->next->prev = ((GainLink *) (link))->prev, \
                                     ((GainLink *) (link))->prev->next = ((GainLink *) (link))->next)
#endif /* ((! defined GAIN) && (! defined SCOTCH_DEBUG_GAIN1)) */
