/* Copyright 2004,2007,2012 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hall_order_hd.c                         **/
/**                                                        **/
/**   AUTHOR     : Patrick AMESTOY                         **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module orders a halo graph or mesh **/
/**                structure using the block-oriented Halo **/
/**                Approximate (Multiple) Minimum Degree   **/
/**                algorithm, with super-variable          **/
/**                accounting (HaloAMD v2.0).              **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 09 aug 1998     **/
/**                                 to     18 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 05 jan 1999     **/
/**                # Version 4.0  : from : 14 jan 2003     **/
/**                                 to   : 29 aug 2007     **/
/**                # Version 6.0  : from : 08 mar 2012     **/
/**                                 to   : 08 mar 2012     **/
/**                                                        **/
/**   NOTES      : # This module contains pieces of code   **/
/**                  that belong to other people; see      **/
/**                  below.                                **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define HALL_ORDER_HD

#include "module.h"
#include "common.h"
#include "graph.h"
#include "hall_order_hd.h"

/*  -- translated by f2c (version 19970219). */

/** -------------------------------------------------------------------- **/
/** December 8th 2003                                                    **/
/** Unique version for both graph of variables and graphs of elements    **/
/** Let us refer to as                                                   **/
/**       Gv a graph with only variables                                 **/
/**       Ge a graph with both variables and elements                    **/
/**                                                                      **/
/** Notations used:                                                      **/
/**                                                                      **/
/**     Let V be the set of nodes                                        **/
/**       V = Ve + V0 + V1                                               **/
/**           V0 = Set of variable nodes (not in halo)                   **/
/**           V1 = Set of variable nodes (in halo)                       **/
/**           Ve = Set of element nodes                                  **/
/**                                                                      **/
/**       All 3 sets are disjoint, Ve and V1 can be empty                **/
/**                                                                      **/
/**  Modifications w.r.t. previous version :                             **/
/**                                                                      **/  
/**  New Input:                                                          **/
/**  ---------                                                           **/
/**         nbelts : integer holding size of Ve                          **/
/**                            =0 if Gv (graph of variables)             **/
/**                            >0 if Ge                                  **/
/**                                                                      **/
/**  Extension of the meaning of input entry len for nodes in Ve         **/
/**  ---------                                                           **/
/**         len(i) = | Adj(i) | if i \in V0 U Ve                         **/
/**           ( Note that in the case of a GE graph                      **/
/**              if v\in V0 then len(v) = nb of elements adjacent to v ) **/
/**         len(i) = - | Adj(i) | if i \in V1                            **/
/**                  or -N -1 if  | Adj(i) | = 0 and i \in V1            **/
/**                                                                      **/
/**  Modified the meaning of input entry elen                            **/
/**  ---------                                                           **/
/**         if e \in Ve then elen (e) = -N-1                             **/
/**         if v \in V0 then elen (v) = External degree of v             **/
/**                             Gv : elen (v) = len(v)                   **/
/**                             Ge : elen (v)                            **/
/**                                  should be computed in SCOTCH        **/
/**         if v \in V1 then elen (v) = 0                                **/
/**                                                                      **/
/**                                                                      **/
/**  Output is unchanged                                                 **/
/**  ---------                                                           **/
/**                                                                      **/
/**                                                                      **/
/** End remarks done on December 8th 2003                                **/
/** ---------------------------------------------------------------------**/

void
hallOrderHdHalmd (
const Gnum          n,                            /* Matrix order                             */
const Gnum          nbelts,                       /* Number of elements                       */
const Gnum          iwlen,                        /* Length of array iw                       */
Gnum * restrict     pe,                           /* Array of indexes in iw of start of row i */
Gnum                pfree,                        /* Useful size in iw                        */
Gnum * restrict     len,                          /* Array of lengths of adjacency lists      */
Gnum * restrict     iw,                           /* Adjacency list array                     */
Gnum * restrict     nv,                           /* Array of element degrees                 */
Gnum * restrict     elen,                         /* Array that holds the inverse permutation */
Gnum * restrict     last,                         /* Array that holds the permutation         */
Gnum * restrict     ncmpa,                        /* Number of times array iw was compressed  */
Gnum * restrict     degree,                       /* Array that holds degree data             */
Gnum * restrict     head,                         /* Linked list structure                    */
Gnum * restrict     next,                         /* Linked list structure                    */
Gnum * restrict     w)                            /* Flag array                               */
{
  Gnum                deg, degme, dext, dmax, e, elenme, eln, hash, hmod, i,
                      ilast, inext, j, jlast, jnext, k, knt1, knt2, knt3,
                      lenj, ln, me = 0, mem, mindeg, nel, newmem,
                      nleft, nvi, nvj, nvpiv, slenme, we, wflg, wnvi, x,
                      nbflag, nreal, lastd, nelme;
  Gnum                p, p1, p2, p3, pdst, pend, pj, pme, pme1, pme2, pn, psrc;

/** -------------------------------------------------------------------- **/
/** HALOAMD_V6: (January 1999, P. Amestoy)                               **/
/** ***********                                                          **/
/**  1/ ERROR 2 detection followed by stop statement suppressed.         **/
/**  2/ Pb 1  identified in V5 was not correctly solved.                 **/
/**                                                                      **/
/** HALOAMD_V5: (December 1998, P. Amestoy)                              **/
/** ***********                                                          **/
/**  1/ Solved problem with matrix psmigr 1, because upper bound degree  **/
/**     DEG>N was considered as a node of V1.                            **/
/**                                                                      **/
/** HALOAMD_V4: (October 1998, P. Amestoy)                               **/
/** ***********                                                          **/
/**  Only MA41 interface (ok for both scotch and MA41) is included in    **/
/**  this file.                                                          **/
/**                                                                      **/
/** HALOAMD_V3: (August 1998, P. Amestoy)                                **/
/** **********                                                           **/
/**  Solved problem in version 2: variables of V1 with len(i)=0 were not **/
/**  well processed. See modification of the input to characterize those **/
/**  variables.                                                          **/
/**  Problem detected by Jacko Koster while experimenting with C version **/
/**  2 of haloAMD in the context of multiple front method based on       **/
/**  MA27: "if for an interface variable i, row i in the matrix has only **/
/**  a nonzero entry on the diagonal, we first remove this entry and     **/
/**  len(i) is set to zero on input to HALOAMD. However, this means that **/
/**  HALOAMD will treat variable i as an interior variable (in V0)       **/
/**  instead as an interface variable (in V1). It is indeed a bit        **/
/**  strange to have such interface variables but we encountered some    **/
/**  in our debugging experiments with some random partitionings.        **/
/**  Solution:                                                           **/
/**  IF on input i \in V1 and len(i)=0 (that is adjlist(i)={}) THEN      **/
/**  len(i) must be set on input to -N-1.                                **/
/**  ENDIF                                                               **/
/**  Therefore, all variables i / len(i) < 0 and only those are in V1.   **/
/**  Variables with len(i) = -N-1 are then processed differently at the  **/
/**  beginning of the code.                                              **/
/**                                                                      **/
/** HALOAMD_V2: (April 1998)                                             **/
/** **********                                                           **/
/**  The end of the tree (including links to block of flagged indices    **/
/**  is built) . The list of flagged indices is considered as a dense    **/
/**  amalgamated node.                                                   **/
/**  Tested on rosanna: ~amestoy/MA41_NEW/SUN_RISC_dbl/SOFT              **/
/**                                                                      **/
/**  Comments on the OUTPUT:                                             **/
/**  ----------------------                                              **/
/**                                                                      **/
/**  Let V= V0 U V1 the nodes of the initial graph (|V|=n).              **/
/**  The assembly tree corresponds to the tree of the supernodes (or     **/
/**  supervariables). Each node of the assembly tree is then composed of **/
/**  one principal variable and a list of secondary variables. The list  **/
/**  of variable of a node (principal + secondary variables) then        **/
/**  describes the structure of the diagonal bloc of the supernode.      **/
/**  The elimination tree denotes the tree of all the variables(=nodes)  **/
/**  and is therefore of order n. The arrays NV(N) and PE(N) give a      **/
/**  description of the assembly tree.                                   **/
/**                                                                      **/
/**   1/ Description of array nv(N) (on OUPUT)                           **/
/**    nv(i)=0 i is a secondary variable.                                **/
/**    N+1> nv(i) >0 i is a principal variable, nv(i) holds the number   **/
/**    of elements in column i of L (true degree of i)                   **/
/**    nv(i) = N+1 then i is a flagged variable (belonging to V1)        **/
/**                                                                      **/
/**   2/ Description of array PE(N) (on OUPUT)                           **/
/**    pe(i) = -(father of variable/node i) in the elimination tree.     **/
/**    If nv (i) .gt. 0, then i represents a node in the assembly tree,  **/
/**    and the parent of i is -pe (i), or zero if i is a root.           **/
/**    If nv (i) = 0, then (i,-pe (i)) represents an edge in a           **/
/**    subtree, the root of which is a node in the assembly tree.        **/
/**                                                                      **/
/**   3/ Example:                                                        **/
/**    Let If be a root node father of Is in the assembly tree.          **/
/**    If is the principal variable of the node If and let If1, If2, If3 **/
/**    be the secondary variables of node If. Is is the principal        **/
/**    variable of the node Is and let Is1, Is2 be the secondary         **/
/**    variables of node Is.                                             **/
/**    Then:                                                             **/
/**        NV(If1)=NV(If2)=NV(If3) = 0  (secondary variables)            **/
/**        NV(Is1)=NV(Is2) = 0  (secondary variables)                    **/
/**        NV(If) > 0  (principal variable)                              **/
/**        NV(Is) > 0  (principal variable)                              **/
/**        PE(If)  = 0 (root node)                                       **/
/**        PE(Is)  = -If (If is the father of Is in the assembly tree)   **/
/**        PE(If1)=PE(If2)=PE(If3)= -If  (If is the principal variable)  **/
/**        PE(Is1)=PE(Is2)= -Is  (Is is the principal variable)          **/
/**                                                                      **/
/** HALOAMD_V1: (September 1997)                                         **/
/** **********                                                           **/
/**  Initial version designed to experiment the numerical (fill-in)      **/
/**  impact of taking into account the halo. This code should be able to **/
/**  experiment no-halo, partial halo, complete halo.                    **/
/** -------------------------------------------------------------------- **/
/** HALOAMD is designed to process a graph composed of two types         **/
/**            of nodes, V0 and V1, extracted from a larger gragh.       **/
/**            V0^V1 = {},                                               **/
/**            We used Min. degree heuristic to order only               **/
/**            nodes in V0, but the adjacency to nodes                   **/
/**            in V1 is taken into account during ordering.              **/
/**            Nodes in V1 are odered at last.                           **/
/**            Adjacency between nodes of V1 need not be provided,       **/
/**            however |len(i)| must always corresponds to the number of **/
/**            edges effectively provided in the adjacency list of i.    **/
/**          On input :                                                  **/
/**          ********                                                    **/
/**            Nodes INODE in V1 are flagged with len(INODE) = -degree   **/
/**            Update version HALO V3 (August 1998):                     **/
/**            if len(i)=0 and i \in V1 then len(i) must be set          **/
/**            on input to -N-1.                                         **/
/**          ERROR return :                                              **/
/**          ************                                                **/
/**            Negative value in ncmpa indicates an error detected       **/
/**            by HALOAMD.                                               **/
/**                                                                      **/
/**            The graph provided MUST follow the rule:                  **/
/**             if (i,j) is an edge in the gragh then                    **/
/**             j must be in the adjacency list of i AND                 **/
/**             i must be in the adjacency list of j.                    **/
/**                                                                      **/
/**    REMARKS :                                                         **/
/**    -------                                                           **/
/**        1/  Providing edges between nodes of V1 should not            **/
/**            affect the final ordering, only the amount of edges       **/
/**            of the halo should effectively affect the solution.       **/
/**            This code should work in the following cases:             **/
/**              1/ halo not provided                                    **/
/**              2/ halo partially provided                              **/
/**              3/ complete halo                                        **/
/**              4/ complete halo+interconnection between nodes of V1.   **/
/**                                                                      **/
/**             1/ should run and provide identical results (w.r.t to    **/
/**                current implementation of AMD in SCOTCH).             **/
/**             3/ and 4/ should provide identical results.              **/
/**                                                                      **/
/**        2/ All modifications of the MC47 initial code are indicated   **/
/**           with begin HALO .. end HALO                                **/
/**                                                                      **/
/** Ordering of nodes in V0 is based on Approximate Minimum Degree       **/
/** ordering algorithm, with aggressive absorption:                      **/
/** Given a representation of the nonzero pattern of a symmetric matrix, **/
/**       A, (excluding the diagonal) perform an approximate minimum     **/
/**       (UMFPACK/MA38-style) degree ordering to compute a pivot order  **/
/**       such that fill-in in the Cholesky                              **/
/**       factors A = LL^T is kept low.  At each step, the pivot         **/
/**       selected is the one with the minimum UMFPACK/MA38-style        **/
/**       upper-bound on the external degree.  Aggresive absorption is   **/
/**       used to tighten the bound on the degree.  This can result an   **/
/**       significant improvement in the quality of the ordering for     **/
/**       some matrices.                                                 **/
/**       The approximate degree algorithm implemented here is the       **/
/**       symmetric analogue of the degree update algorithm in MA38, by  **/
/**       Davis and Duff, also in the Harwell Subroutine Library.        **/
/**                                                                      **/
/** **** CAUTION:  ARGUMENTS ARE NOT CHECKED FOR ERRORS ON INPUT.  ***** **/
/** ** If you want error checking, a more versatile input format, and ** **/
/** ** a simpler user interface, then use MC47A/AD in the Harwell     ** **/
/** ** Subroutine Library, which checks for errors, transforms the    ** **/
/** ** input, and calls MC47B/BD.                                     ** **/
/** ******************************************************************** **/
/**       References:  (UF Tech Reports are available via anonymous ftp  **/
/**       to ftp.cis.ufl.edu:cis/tech-reports).                          **/
/**       [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern    **/
/**               multifrontal method for sparse LU factorization",      **/
/**               SIAM J. Matrix Analysis and Applications, to appear.   **/
/**               also Univ. of Florida Technical Report TR-94-038.      **/
/**               Discuss UMFPACK / MA38.                                **/
/**       [2] Patrick Amestoy, Timothy A. Davis, and Iain S. Duff,       **/
/**               "An approximate minimum degree ordering algorithm,"    **/
/**               SIAM J. Matrix Analysis and Applications (to appear),  **/
/**               also Univ. of Florida Technical Report TR-94-039.      **/
/**               Discusses this routine.                                **/
/**       [3] Alan George and Joseph Liu, "The evolution of the          **/
/**               minimum degree ordering algorithm," SIAM Review, vol.  **/
/**               31, no. 1, pp. 1-19, March 1989.  We list below the    **/
/**               features mentioned in that paper that this code        **/
/**               includes:                                              **/
/**       mass elimination:                                              **/
/**               Yes.  MA27 relied on supervariable detection for mass  **/
/**               elimination.                                           **/
/**       indistinguishable nodes:                                       **/
/**               Yes (we call these "supervariables").  This was also   **/
/**               in the MA27 code - although we modified the method of  **/
/**               detecting them (the previous hash was the true degree, **/
/**               which we no longer keep track of).  A supervariable is **/
/**               a set of rows with identical nonzero pattern.  All     **/
/**               variables in a supervariable are eliminated together.  **/
/**               Each supervariable has as its numerical name that of   **/
/**               one of its variables (its principal variable).         **/
/**       quotient graph representation:                                 **/
/**               Yes.  We use the term "element" for the cliques formed **/
/**               during elimination.  This was also in the MA27 code.   **/
/**               The algorithm can operate in place, but it will work   **/
/**               more efficiently if given some "elbow room."           **/
/**       element absorption:                                            **/
/**               Yes.  This was also in the MA27 code.                  **/
/**       external degree:                                               **/
/**               Yes.  The MA27 code was based on the true degree.      **/
/**       incomplete degree update and multiple elimination:             **/
/**               No.  This was not in MA27, either.  Our method of      **/
/**               degree update within MC47B/BD is element-based, not    **/
/**               variable-based.  It is thus not well-suited for use    **/
/**               with incomplete degree update or multiple elimination. **/
/** -------------------------------------------------------------------- **/
/** Authors, and Copyright (C) 1995 by:                                  **/
/**       Timothy A. Davis, Patrick Amestoy, Iain S. Duff, &             **/
/**       John K. Reid.                                                  **/
/** Modified (V1) by P.R. Amestoy ENSEEIHT (1997)                        **/
/** Modified (V2) by P.R. Amestoy ENSEEIHT (1998)                        **/
/** Modified (V3) by P.R. Amestoy ENSEEIHT (1998)                        **/
/** Modified (V4) by P.R. Amestoy ENSEEIHT (1998)                        **/
/** Modified (V5) by P.R. Amestoy ENSEEIHT (1998)                        **/
/** Modified (V6) by P.R. Amestoy ENSEEIHT (1999)                        **/
/**                                                                      **/
/** Dates: September, 1995                                               **/
/**        September, 1997 (halo AMD V1)                                 **/
/**        April, 1998 (halo AMD V2)                                     **/
/**        August, 1998 (halo AMD V3)                                    **/

  -- w;                                           /* Parameter adjustments */
  -- next;
  -- head;
  -- degree;
  -- last;
  -- elen;
  -- nv;
  -- len;
  -- pe;
  -- iw;

  wflg = 2;
  mindeg = 1;
  *ncmpa = 0;
  nel = 0;
  hmod = MAX (1, (n - 1));
  dmax = 0;
  mem = pfree - 1;
  nbflag = 0;
  lastd = 0;

  memSet (last + 1, 0, n * sizeof (Gnum));
  memSet (head + 1, 0, n * sizeof (Gnum));
  
  if (nbelts == 0) {                              /* Patch 8/12/03 <PA> */
    memSet (elen + 1, 0, n * sizeof (Gnum));      
    for (i = 1; i <= n; i ++) {
      nv[i] = 1;
      w[i]  = 1;
      if (len[i] < 0) {
        degree[i] = n + 1;
        nbflag ++;
        if (len[i] == - (n + 1)) {                /* Patch 09/08/98 <PA+FP> */
          len[i] = 0;
          pe[i]  = 0;                             /* Patch 12/12/03 <PA>: Because of compress, we force skipping those entries (which are anyway empty) */
        }
        else
          len[i] = - len[i];
      }
      else
        degree[i] = len[i];
    }
  }
  else  {                                         /* Patch 08/12/03 <PA>: Duplicate part of previous loop to avoid sytematic testing for elements */
    for (i = 1; i <= n; i ++) {
      nv[i] = 1;
      w[i]  = 1;
      if (len[i] < 0) {                           /* i \in V1 */
        degree[i] = n + 1;
        nbflag ++;
        if (len[i] == - (n + 1)) {                /* Patch 09/08/98 <PA+FP> */
          len[i]  = 0;
          pe[i]   = 0;                            /* Patch 12/12/03 <PA>: because of compress, we force skipping those entries (which are anyway empty) */
          elen[i] = 0;                            /* Patch 16/12/03 <PA> */
        }
        else {
          len[i]  = - len[i];
          elen[i] = len[i];                       /* Patch 16/12/03 <PA>: only elements are adjacent to a variable */
        }
      }
      else {                                      /* i \in Ve or V0 */
        if (elen[i] < 0) {                        /* i \in Ve       */
          nel ++;
          degree[i] = len[i];
          elen[i]   = - nel;
          dmax      = MAX (dmax, degree[i]);      /* Patch 11/03/04 <PA> */
        }
        else {
          degree[i] = elen[i];
          elen[i]   = len[i];                     /* Patch 16/12/03 <PA>: only elements are adjacent to a variable */
        }
      }
    }
  }

#ifdef SCOTCH_DEBUG_ORDER2
  if (nbelts != nel)                              /* Temporary Patch 8/12/03 <PA> */
    printf ("error 8Dec2003\n");
#endif /* SCOTCH_DEBUG_ORDER2 */

  nreal = n - nbflag;

  for (i = 1; i <= n; i ++) {
    if (elen[i] < 0 )                             /* Patch 16/12/03 <PA>: Skip elements */
      continue;

    deg = degree[i];
    if (deg == (n + 1)) {
      deg = n;
      if (lastd == 0) {
        lastd     = i;
        head[deg] = i;
        next[i]   = 0;
        last[i]   = 0;
      }
      else {
        next[lastd] = i;
        last[i]     = lastd;
        lastd       = i;
        next[i]     = 0;
      }
    }
    else if (deg > 0) {
      inext = head[deg];
      if (inext != 0)
        last[inext] = i;
      next[i]   = inext;
      head[deg] = i;
    }
    else {
      nel ++;
      elen[i] = - nel;
      pe[i]   = 0;
      w[i]    = 0;
    }
  }                                               /* L20: */

  nleft = n - nel;                                /* Patch v5 12/12/98 <PA+FP> */

  while (nel < nreal) {                           /* WHILE (selecting pivots) DO */
    for (deg = mindeg; deg <= n; deg ++) {        /* Patch 17/11/97 <PA+FP>      */
       me = head[deg];
       if (me > 0)
         break;                                   /* GO to 50 */
    }                                             /* L40:     */
    mindeg = deg;
    if (me <= 0) {                                /* Error 1 */
      *ncmpa = -n;
      return;
    }

    inext = next[me];
    if (inext != 0)
      last[inext] = 0;
    head[deg] = inext;

    elenme   = elen[me];
    elen[me] = - (nel + 1);
    nvpiv    = nv[me];
    nel     += nvpiv;

    nv[me] = - nvpiv;
    degme  = 0;
    if (elenme == 0) {
      pme1 = pe[me];
      pme2 = pme1 - 1;

      for (p = pme1; p <= pme1 + len[me] - 1; p ++) {
        i   = iw[p];
        nvi = nv[i];
        if (nvi > 0) {
          degme +=   nvi;
          nv[i]  = - nvi;
          pme2 ++;
          iw[pme2] = i;

          if (degree[i] <= n) {
            ilast = last[i];
            inext = next[i];
            if (inext != 0)
              last[inext] = ilast;
            if (ilast != 0)
              next[ilast] = inext;
            else
              head[degree[i]] = inext;
          }
        }
      }                                           /* L60: */

      newmem = 0;
    }
    else {
      p    = pe[me];
      pme1 = pfree;
      slenme = len[me] - elenme;
      for (knt1 = 1; knt1 <= elenme + 1; knt1 ++) {
        if (knt1 > elenme) {
          e  = me;
          pj = p;
          ln = slenme;
        }
        else {
          e  = iw[p ++];
          pj = pe[e];
          ln = len[e];
        }

        for (knt2 = 1; knt2 <= ln; knt2 ++) {
          i   = iw[pj ++];
          nvi = nv[i];
          if (nvi > 0) {
            if (pfree > iwlen) {
              pe[me]   = p;
              len[me] -= knt1;
              if (len[me] == 0)
                pe[me] = 0;
              pe[e]  = pj;
              len[e] = ln - knt2;
              if (len[e] == 0)
                pe[e] = 0;
              (*ncmpa) ++;

              for (j = 1; j <= n; j ++) {
                pn = pe[j];
                if (pn > 0) {
                  pe[j]  = iw[pn];
                  iw[pn] = - j;
                }
              }                                   /* L70: */

              pdst = 1;
              psrc = 1;
              pend = pme1 - 1;

              while (psrc <= pend) {              /* L80: */
                j = - iw[psrc ++];
                if (j > 0) {
                  iw[pdst] = pe[j];
                  pe[j]    = pdst ++;
                  lenj     = len[j];
                  for (knt3 = 0; knt3 <= lenj - 2; knt3 ++)
                    iw[pdst + knt3] = iw[psrc + knt3];
                  pdst = pdst + (lenj - 1);
                  psrc = psrc + (lenj - 1);
                }
              }

              p1 = pdst;
              for (psrc = pme1; psrc <= pfree - 1; psrc ++, pdst ++) /* L100: */
                iw[pdst] = iw[psrc];
              pme1 = p1;
              pfree = pdst;
              pj = pe[e];
              p  = pe[me];
            }

            degme +=   nvi;
            nv[i]  = - nvi;
            iw[pfree] = i;
            pfree ++;

            if (degree[i] <= n) {
              ilast = last[i];
              inext = next[i];
              if (inext != 0)
                last[inext] = ilast;
              if (ilast != 0)
                next[ilast] = inext;
              else
                head[degree[i]] = inext;
            }
          }
        }                                         /* L110: */

        if (e != me) {
          pe[e] = -me;
          w[e]  = 0;
        }
      }                                           /* L120: */
      pme2 = pfree - 1;

      newmem = pfree - pme1;
      mem   += newmem;
    }

    degree[me] = degme;
    pe[me]     = pme1;
    len[me]    = pme2 - pme1 + 1;

    if (wflg + n <= wflg) {
      for (x = 1; x <= n; x ++) {
        if (w[x] != 0)
          w[x] = 1;
      }                                           /* L130: */
      wflg = 2;
    }

    for (pme = pme1; pme <= pme2; pme ++) {
      i   = iw[pme];
      eln = elen[i];
      if (eln > 0) {
        nvi  = - nv[i];
        wnvi = wflg - nvi;
        for (p = pe[i]; p < pe[i] + eln; p ++) {
          e  = iw[p];
          we = w[e];
          if (we >= wflg)
            we -= nvi;
          else if (we != 0)
            we = degree[e] + wnvi;
          w[e] = we;
        }                                         /* L140: */
      }
    }                                             /* L150: */

    for (pme = pme1; pme <= pme2; pme ++) {
      i  = iw[pme];
      p1 = pe[i];
      p2 = p1 + elen[i] - 1;
      pn = p1;
      hash = 0;
      deg  = 0;

      for (p = p1; p <= p2; p ++) {
        e    = iw[p];
        dext = w[e] - wflg;
        if (dext > 0) {
          deg      += dext;
          iw[pn ++] = e;
          hash     += e;
        }
        else if (dext == 0) {
          pe[e] = -me;
          w[e]  = 0;
        }
      }                                           /* L160: */
      elen[i] = pn - p1 + 1;

      p3 = pn;
      for (p = p2 + 1; p < p1 + len[i]; p ++) {
        j   = iw[p];
        nvj = nv[j];
        if (nvj > 0) {
          deg += nvj;
          iw[pn ++] = j;
          hash += j;
        }
      }                                           /* L170: */

      if (degree[i] == (n + 1))
        deg = n + 1;
      if (deg == 0) {
        pe[i]   = - me;
        nvi     = - nv[i];
        degme  -= nvi;
        nvpiv  += nvi;
        nel    += nvi;
        nv[i]   = 0;
        elen[i] = 0;
      }
      else {
        if (degree[i] != (n + 1)) {               /* Patch v6 05/01/99 <PA+FP> */
          deg       = MIN (nleft,     deg);       /* Patch v5 12/12/98 <PA+FP> */
          degree[i] = MIN (degree[i], deg);
        }

        iw[pn] = iw[p3];
        iw[p3] = iw[p1];
        iw[p1] = me;
        len[i] = pn - p1 + 1;

        if (deg <= n) {
          hash = (hash % hmod) + 1;
          j = head[hash];
          if (j <= 0) {
            next[i]    = - j;
            head[hash] = - i;
          }
          else {
            next[i] = last[j];
            last[j] = i;
          }
          last[i] = hash;
        }
      }
    }                                             /* L180: */
    degree[me] = degme;

    dmax  = MAX (dmax, degme);
    wflg += dmax;

    if (wflg + n <= wflg) {
      for (x = 1; x <= n; x ++) {
        if (w[x] != 0)
          w[x] = 1;
      }
      wflg = 2;
    }

    for (pme = pme1; pme <= pme2; pme ++) {
      i = iw[pme];
      if ((nv[i] < 0) && (degree[i] <= n)) {
        hash = last[i];
        j    = head[hash];
        if (j == 0)
          continue;
        if (j < 0) {
          i = - j;
          head[hash] = 0;
        }
        else {
          i       = last[j];
          last[j] = 0;
        }
        if (i == 0)
          continue;

L200:                                             /* WHILE LOOP: */
        if (next[i] != 0) {
          ln  = len[i];
          eln = elen[i];
          for (p = pe[i] + 1; p < pe[i] + ln; p ++)
            w[iw[p]] = wflg;

          jlast = i;
          j = next[i];

L220:                                             /* WHILE LOOP: */
          if (j != 0) {
            if (len[j] != ln)
              goto L240;
            if (elen[j] != eln)
              goto L240;

            for (p = pe[j] + 1; p < pe[j] + ln; p ++) {
              if (w[iw[p]] != wflg)
                goto L240;
            }                                     /* L230: */

            pe[j]   = -i;
            nv[i]  += nv[j];
            nv[j]   = 0;
            elen[j] = 0;

            j           = next[j];
            next[jlast] = j;
            goto L220;

L240:
            jlast = j;
            j     = next[j];
            goto L220;
          }

          wflg ++;
          i = next[i];
          if (i != 0)
            goto L200;
        }
      }
    }

    p     = pme1;
    nleft = n - nel;
    for (pme = pme1; pme <= pme2; pme ++) {
      i   = iw[pme];
      nvi = - nv[i];
      if (nvi > 0) {
        nv[i] = nvi;
        if (degree[i] <= n) {
          deg = MIN (degree[i] + degme, nleft) - nvi;

          inext = head[deg];
          if (inext != 0)
            last[inext] = i;
          next[i]   = inext;
          last[i]   = 0;
          head[deg] = i;

          mindeg    = MIN (mindeg, deg);
          degree[i] = deg;
        }

        iw[p ++] = i;
      }
    } /* L260: */

    nv[me]  = nvpiv + degme;
    len[me] = p - pme1;
    if (len[me] == 0) {
      pe[me] = 0;
      w[me]  = 0;
    }
    if (newmem != 0) {
      pfree = p;
      mem   = mem - newmem + len[me];
    }
  }                                             /* END WHILE (selecting pivots) */

  if (nel < n) {                                /* Patch 12/12/98 <PA+FP> (old: nreal < n) */
    for (deg = mindeg; deg <= n; deg ++) {
      me = head[deg];
      if (me > 0)
        break;
    }

    mindeg = deg;
    nelme  = - (nel + 1);
    for (x = 1; x <= n; x ++) {
      if ((pe[x] > 0) && (elen[x] < 0))
        pe[x] = - me;
      else if (degree[x] == (n + 1)) {
        nel    += nv[x];
        pe[x]   = - me;
        elen[x] = 0;
        nv[x]   = 0;                              /* Patch 12/12/98 <PA+FP> (old: n + 1) */
      }
    }

    elen[me] = nelme;
    nv[me]   = n - nreal;                         /* Patch 12/12/98 <PA+FP> (old: n + 1) */
    pe[me]   = 0;
    if (nel != n) {                               /* Error 2 */
      *ncmpa = - (n + 1);
      return;
    }
  }

  for (i = 1; i <= n; i ++) {
    if (elen[i] == 0) {
      j = - pe[i];

      while (elen[j] >= 0)                        /* L270: */
        j = - pe[j];
      e = j;

      k = - elen[e];
      j = i;

      while (elen[j] >= 0) {                      /* L280: */
        jnext = - pe[j];
        pe[j] = - e;
        if (elen[j] == 0)
          elen[j] = k ++;
        j = jnext;
      }
      elen[e] = - k;
    }
  }                                               /* L290: */

#ifdef DEAD_CODE                                  /* No need for permutations */
  for (i = 1; i <= n; i ++) {                     /* Patch 19/10/98 <PA+FP>   */
    k = abs (elen[i]);
    last[k] = i;
    elen[i] = k;
  }                                               /* L300: */
#endif /* DEAD_CODE */
}
