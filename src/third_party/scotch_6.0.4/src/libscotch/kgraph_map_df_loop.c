/* Copyright 2010-2012,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kgraph_map_df_loop.c                    **/
/**                                                        **/
/**   AUTHOR     : Sebastien FOURESTIER (v6.0)             **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes a k-way partition  **/
/**                of the given mapping graph by applying  **/
/**                a diffusion method to what is assumed   **/
/**                to be a band graph.                     **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 05 jan 2010     **/
/**                                 to   : 23 aug 2014     **/
/**                                                        **/
/************************************************************/

/* Tests flags for mapping TODO remove it after performance tests */
/* #define KGRAPHDIFFMAPPNONE */                    /* No special code for mapping                      */
/* #define KGRAPHDIFFMAPPMORE */                    /* Give more liquid on expensive architecture edges */
#define KGRAPHDIFFMAPPLESS                        /* Give less liquid on expensive architecture edges */

/*
**  The defines and includes.
*/

/****************************/
/*                          */
/* The diffusion subroutine */
/* pattern.                 */
/*                          */
/****************************/

/* This routine computes the diffusion of liquids
** on the given part of the k-way band graph.
** It returns:
** - void  : in all cases
*/

static
int
KGRAPHMAPDFLOOPNAME (
KgraphMapDfThread * restrict  thrdptr)            /* Thread-dependent data */
{
  KgraphMapDfVertex * restrict  difotax;          /* Old diffusion value array  */
  KgraphMapDfVertex * restrict  difntax;          /* New diffusion value array  */
  KgraphMapDfSort * restrict    sorttab;          /* Liquid sort array          */
  Gnum                          vancnnd;
  Gnum                          vertnum;
  Anum                          domnnum;
  Gnum                          passnum;
  int                           velsmsk;
  int                           mappflag;         /* We are computing a mapping */

  KgraphMapDfData * restrict const  loopptr = (KgraphMapDfData *) thrdptr->thrddat.grouptr;
  const Kgraph * restrict const     grafptr = loopptr->grafptr;
  float * restrict const            vanctab = loopptr->vanctab;
  float * restrict const            valotab = loopptr->valotab; /* Fraction of load to leak */
  Gnum * restrict const             velstax = loopptr->velstax;
  const Gnum                        vertbas = thrdptr->vertbas; /* Range of non-anchor vertices to process */
  const Gnum                        vertnnd = thrdptr->vertnnd;
  const Anum                        domnbas = thrdptr->domnbas; /* Range of anchor vertices to process */
  const Anum                        domnnnd = thrdptr->domnnnd;
  const Anum                        domnnbr = grafptr->m.domnnbr;
  const Gnum                        crloval = grafptr->r.crloval;
  const Gnum                        cmloval = grafptr->r.cmloval;
  Anum * restrict const             parttax = grafptr->m.parttax;
  const Anum * restrict const       parotax = grafptr->r.m.parttax;
  const Gnum * restrict const       verttax = grafptr->s.verttax;
  const Gnum * restrict const       vendtax = grafptr->s.vendtax;
  const Gnum * restrict const       velotax = grafptr->s.velotax;
  const Gnum * const                edgetax = grafptr->s.edgetax;
  const Gnum * const                edlotax = grafptr->s.edlotax;

  vancnnd = grafptr->s.vertnnd - domnnbr;
  velsmsk = 1;                                    /* Assume no anchors are isolated */
  if (edlotax != NULL) {
    for (domnnum = domnbas; domnnum < domnnnd; domnnum ++) { /* For all local anchor vertices */
      Gnum                edgenum;
      Gnum                edgennd;
      Gnum                velssum;

      for (edgenum = verttax[vancnnd + domnnum], edgennd = vendtax[vancnnd + domnnum], velssum = 0;
           edgenum < edgennd; edgenum ++)
        velssum += edlotax[edgenum];
      velstax[vancnnd + domnnum] = velssum;
      velsmsk &= (velssum != 0);
    }
  }
  else {
    for (domnnum = domnbas; domnnum < domnnnd; domnnum ++) {
      Gnum                velssum;

      velssum = vendtax[vancnnd + domnnum] - verttax[vancnnd + domnnum]; /* Local degree of anchor vertices */
      velstax[vancnnd + domnnum] = velssum;
      velsmsk &= (velssum != 0);
    }
  }
  if (velsmsk == 0) {                             /* If graph is too small to have any usable anchors */
#ifdef KGRAPHMAPDFLOOPTHREAD
    loopptr->abrtval == 1;                        /* We will leave during the first iteration */
#else /* KGRAPHMAPDFLOOPTHREAD */
    return (1);
#endif /* KGRAPHMAPDFLOOPTHREAD */
  }

  if ((sorttab = memAlloc (domnnbr * sizeof (KgraphMapDfSort))) == NULL) { /* Allocate here for memory affinity as it is a private array */
    errorPrint (STRINGIFY (KGRAPHMAPDFLOOPNAME) ": out of memory");
#ifdef KGRAPHMAPDFLOOPTHREAD
    loopptr->abrtval = 1;
#else /* KGRAPHMAPDFLOOPTHREAD */
    return (1);
#endif /* KGRAPHMAPDFLOOPTHREAD */
  }
    
  if (velotax == NULL) {
    for (domnnum = domnbas; domnnum < domnnnd; domnnum ++)
      valotab[domnnum] = 1.0F;
  }
  else {
    for (domnnum = domnbas; domnnum < domnnnd; domnnum ++)
      valotab[domnnum] = (float) velotax[vancnnd + domnnum];
  }

  difntax = loopptr->difntax;
  difotax = loopptr->difotax;

  if (edlotax != NULL) {
    for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) {
      Gnum                velssum;
      Gnum                edgenum;
      Gnum                edgennd;

#ifdef SCOTCH_DEBUG_KGRAPH2
      if ((vendtax[vertnum] - verttax[vertnum]) == 0) {
        errorPrint (STRINGIFY (KGRAPHMAPDFLOOPNAME) ": internal error (1)");
#ifdef KGRAPHMAPDFLOOPTHREAD
        loopptr->abrtval = 1;
#else /* KGRAPHMAPDFLOOPTHREAD */
        return (1);
#endif /* KGRAPHMAPDFLOOPTHREAD */
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      difotax[vertnum].partval = parttax[vertnum]; /* Set initial part by default */
      difotax[vertnum].diffval =
      difotax[vertnum].fdifval =
      difotax[vertnum].mdisval =
      difotax[vertnum].mdidval =
      difntax[vertnum].fdifval =
      difntax[vertnum].mdisval =
      difntax[vertnum].mdidval = 0.0F;

      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum], velssum = 0;
           edgenum < edgennd; edgenum ++)
        velssum += edlotax[edgenum];
      velstax[vertnum] = velssum;
    }
  }
  else {                                          /* Graph has no edge loads */
    for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) {
#ifdef SCOTCH_DEBUG_KGRAPH2
      if ((vendtax[vertnum] - verttax[vertnum]) == 0) {
        errorPrint (STRINGIFY (KGRAPHMAPDFLOOPNAME) ": internal error (2)");
#ifdef KGRAPHMAPDFLOOPTHREAD
        loopptr->abrtval = 1;
#else /* KGRAPHMAPDFLOOPTHREAD */
        return (1);
#endif /* KGRAPHMAPDFLOOPTHREAD */
      }
#endif /* SCOTCH_DEBUG_KGRAPH2 */
      difotax[vertnum].partval = parttax[vertnum]; /* Set initial part by default */
      difotax[vertnum].diffval =
      difotax[vertnum].fdifval = 
      difotax[vertnum].mdisval =
      difotax[vertnum].mdidval =
      difntax[vertnum].fdifval = 
      difntax[vertnum].mdisval =
      difntax[vertnum].mdidval = 0.0F;
      velstax[vertnum] = vendtax[vertnum] - verttax[vertnum];
    }
  }
  for (domnnum = domnbas, vertnum = vancnnd + domnbas; /* For all the subset of anchor vertices */
       domnnum < domnnnd; domnnum ++, vertnum ++) {
    float               vancval;
    Gnum                comploadbal;              /* Compload to reach to get wished balance */

    if (velstax[vancnnd + domnnum] <= 0) {
      vancval          =
      vanctab[domnnum] = 0.0F; 
      velstax[vertnum] = -1;
    }
    else {
      comploadbal = grafptr->comploadavg[domnnum];
      vancval = ((float) comploadbal - valotab[domnnum]) / (float) velstax[vancnnd + domnnum]; /* Amount of liquid to be added at each step */
      vanctab[domnnum] = comploadbal;
    }
    difotax[vertnum].diffval = vancval;           /* Load anchor vertices for first pass */
    difotax[vertnum].partval =
    difntax[vertnum].partval = domnnum;
    difntax[vertnum].diffval =                    /* In case of isolated anchors, do not risk overflow because of NaN */
    difotax[vertnum].fdifval =
    difotax[vertnum].mdisval =
    difotax[vertnum].mdidval = 
    difntax[vertnum].fdifval =
    difntax[vertnum].mdisval =
    difntax[vertnum].mdidval = 0.0F;              /* Do not consider migration costs for anchors */
  }

#ifdef KGRAPHMAPDFLOOPTHREAD
  threadBarrier (thrdptr);

  if (loopptr->abrtval == 1) {                    /* If process alone or some decided to quit */
    memFree (sorttab);                            /* Free local array                         */
    return  (1);
  }
#endif /* KGRAPHMAPDFLOOPTHREAD */

#ifndef KGRAPHDIFFMAPPNONE
  if (archPart (grafptr->m.archptr))   
    mappflag = 0;
  else
    mappflag = 1;
#else
  mappflag = 0;
#endif /* KGRAPHDIFFMAPPNONE */

  passnum = loopptr->passnbr;
  for ( ; passnum > 0; passnum --) {              /* For all passes       */
    KgraphMapDfVertex * difttax;                  /* Temporary swap value */
    Gnum                vertnum;
    float               veloval;

    veloval = 1.0F;                               /* Assume no vertex loads */

    for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) { /* For all local regular vertices */
      Gnum                edgenum;
      Gnum                edgennd;
      Gnum                soplval;                /* Load sum of edges going to vertex old part                 */
      Gnum                sfplval;                /* Load sum of edges going to vertex of other parts           */
      Gnum                dfplval;                /* Load sum of edges going to vertex of other parts * distval */
      Gnum                migrval;
      Anum                partnbr;                /* Number of active parts */
      Anum                partnum;
      float               diffval;
      float               signval;
      Anum                paronum;
      Anum                partcur;

      partnbr = 1;                                /* Keep vertex in first place to preserve its part */
      partcur            =
      sorttab[0].partval = difotax[vertnum].partval; /* Always keep old part value                */
      sorttab[0].diffval = 0.0F;                  /* Assume at first it is null                   */
      sorttab[0].edlosum = 0;                     /* Assume at first there are no loads           */
      sorttab[0].distval = 1;                     /* Do not take distval of our part into account */

      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum];
           edgenum < edgennd; edgenum ++) {
        Gnum                vertend;
        float               diffval;
        float               fdifval;
        float               mdisval;
        float               mdidval;
        Anum                partval;
        Anum                partnum;
        Gnum                edloval;
        Anum                distval;

        vertend = edgetax[edgenum];
        edloval = (float) ((edlotax != NULL) ? edlotax[edgenum] : 1);

        partval = difotax[vertend].partval;
        diffval = difotax[vertend].diffval;       /* Value is not yet scaled with respect to diffusion coefficient */
        fdifval = difotax[vertend].fdifval;
        mdisval = difotax[vertend].mdisval;
        mdidval = difotax[vertend].mdidval;

        if ((mappflag == 1) && (partval != partcur))
          diffval = fdifval;

        diffval *= edloval * crloval;
        if (parotax != NULL) {
          if (difotax[vertnum].partval == parotax[vertend])
            diffval += mdisval;
          else
            diffval += mdidval;
        }

        for (partnum = 0; partnum < partnbr; partnum ++) {
          if (sorttab[partnum].partval == partval) {
            sorttab[partnum].diffval += diffval;  /* Accumulate contribution in slot */
            sorttab[partnum].edlosum += edloval;
            goto endloop1;                        /* Do not consider creating a new slot */
          }
        }
        sorttab[partnbr].partval = partval;       /* Create new slot */
        sorttab[partnbr].distval = ((mappflag == 1) && (partcur != partval)) ? archDomDist (&grafptr->a, &grafptr->m.domntab[partcur], &grafptr->m.domntab[partval]) : 1;
        sorttab[partnbr].diffval = diffval;
        sorttab[partnbr].edlosum = edloval;
        partnbr ++;
endloop1 : ;
      }

      if (mappflag == 1)
        for (partnum = 0; partnum < partnbr; partnum ++)
#ifdef KGRAPHDIFFMAPPMORE
          sorttab[partnum].diffval *= sorttab[partnum].distval;
#else /* KGRAPHDIFFMAPPLESS */
          sorttab[partnum].diffval /= sorttab[partnum].distval;
#endif /* KGRAPHDIFFMAPPMORE */

      if (partnbr > 1)                            /* If not interior vertex           */
        kgraphMapDfSort (sorttab, partnbr);       /* Sort array by descending amounts */

      soplval = 0;
      if (parotax != NULL) {
        for (partnum = 0; partnum < partnbr; partnum ++) {
          if (sorttab[partnum].partval == parotax[vertnum]) {
            soplval = sorttab[partnum].edlosum;
            break;
          }
        }
      }

      sfplval = 0;
      dfplval = 0;
      if (mappflag == 1) {                        /* We are doing a mapping */
        for (partnum = 1; partnum < partnbr; partnum ++) {
          sfplval += sorttab[partnum].edlosum;
#ifdef KGRAPHDIFFMAPPMORE
          dfplval += sorttab[partnum].edlosum * sorttab[partnum].distval;
#else /* KGRAPHDIFFMAPPLESS */
          dfplval += sorttab[partnum].edlosum / sorttab[partnum].distval;
#endif /* KGRAPHDIFFMAPPMORE */
        }
      }

      difntax[vertnum].partval = sorttab[0].partval; /* New part is part of most abundant liquid */

      diffval = sorttab[0].diffval;               /* Get amount of most abundant liquid    */

      if (velotax != NULL)                        /* Account for capacity of barrel        */
        veloval = (float) velotax[vertnum];
      diffval -= veloval;                         /* Leak liquid from barrel               */
      if (diffval <= 0.0F)                        /* Amount of liquid cannot be negative   */
        diffval = 0.0F;
      migrval = ((soplval == 0) || (soplval == velstax[vertnum])) ? 0 : grafptr->r.cmloval * ((grafptr->r.vmlotax != NULL) ? grafptr->r.vmlotax[vertnum] : 1);
      if (migrval > diffval) {
        migrval = diffval;
        diffval = 0;
      }
      else
        diffval -= migrval;
      diffval = diffval / (velstax[vertnum] * crloval);
      if (isnan (diffval)) {                      /* If overflow occured */
#ifdef SCOTCH_DEBUG_KGRAPH2
        errorPrintW (STRINGIFY (KGRAPHMAPDFLOOPNAME) ": overflow (1)");
#endif /* SCOTCH_DEBUG_KGRAPH2 */
#ifdef KGRAPHMAPDFLOOPTHREAD
        loopptr->abrtval = 1;                     /* Threads need to halt              */
        goto abort1;                              /* Skip computations but synchronize */
#else /* KGRAPHMAPDFLOOPTHREAD */
        goto abort2;                              /* Exit this loop without swapping arrays */
#endif /* KGRAPHMAPDFLOOPTHREAD */
      }

      if (parotax != NULL) {
        if (migrval == 0) {
          difntax[vertnum].mdisval = 
          difntax[vertnum].mdidval = 0;
        }
        else {
          if (parotax[vertnum] == sorttab[0].partval) {
            difntax[vertnum].mdisval = migrval / soplval;
            difntax[vertnum].mdidval = 0;
          }  
          else {
            difntax[vertnum].mdisval = 0;
            difntax[vertnum].mdidval = migrval / (velstax[vertnum] - soplval);
          }
        }
      }

      difntax[vertnum].diffval = diffval;
      if (dfplval != 0)
        difntax[vertnum].fdifval = diffval * sfplval / dfplval;
      else
        difntax[vertnum].fdifval = 0;
    }

    for (domnnum = domnbas, vertnum = vancnnd + domnbas; /* For all the subset of anchor vertices */
         domnnum < domnnnd; domnnum ++, vertnum ++) {
      Gnum                edgenum;
      Gnum                edgennd;
      Anum                partnbr;                /* Number of active parts */
      Anum                partnum;
      float               diffval;
      float               signval;

      partnbr = 1;                                /* Keep vertex in first place to preserve its part */
      sorttab[0].partval = domnnum;               /* Always keep initial part value                  */
      sorttab[0].diffval = 0.0F;                  /* Assume at first it is null                      */

      edgenum = verttax[vertnum];
      edgennd = vendtax[vertnum];
      if (edgenum == edgennd)                     /* If isolated anchor */
        continue;                                 /* Barrel is empty    */

      for ( ; edgenum < edgennd; edgenum ++) {    /* For all edges except anchors */
        Gnum                vertend;
        float               diffval;
        Anum                partval;
        Anum                partnum;

        vertend = edgetax[edgenum];

        partval = difotax[vertend].partval;
        diffval = difotax[vertend].diffval;       /* Value is not yet scaled with respect to diffusion coefficient */

        diffval *= (float) ((edlotax != NULL) ? edlotax[edgenum] : 1);
        diffval *= crloval;

        for (partnum = 0; partnum < partnbr; partnum ++) {
          if (sorttab[partnum].partval == partval) {
            sorttab[partnum].diffval += diffval;  /* Accumulate contribution in slot     */
            goto endloop2;                        /* Do not consider creating a new slot */
          }
        }
        sorttab[partnbr].partval = partval;       /* Create new slot */
        sorttab[partnbr].diffval = diffval;
        partnbr ++;
endloop2 : ;
      }

      if (partnbr > 1)                            /* If not interior vertex           */
        kgraphMapDfSort (sorttab, partnbr);       /* Sort array by descending amounts */

      diffval = sorttab[0].diffval;               /* Get amount of most abundant liquid */

      if (sorttab[0].partval != domnnum)          /* Add liquid from tap to barrel */
        diffval = vanctab[domnnum] - diffval;
      else
        diffval += vanctab[domnnum];

      diffval = (diffval - valotab[domnnum]) / (velstax[vertnum] * crloval); /* Add input and leak liquid from barrel */

      if (diffval <= 0.0F)                        /* Amount of liquid cannot be negative   */
        diffval = 0.0F;
      if (isnan (diffval)) {                      /* If overflow occured */
#ifdef SCOTCH_DEBUG_KGRAPH2
        errorPrintW (STRINGIFY (KGRAPHMAPDFLOOPNAME) ": overflow (2)");
#endif /* SCOTCH_DEBUG_KGRAPH2 */
#ifdef KGRAPHMAPDFLOOPTHREAD
        loopptr->abrtval = 1;                     /* Threads need to halt              */
        goto abort1;                              /* Skip computations but synchronize */
#else /* KGRAPHMAPDFLOOPTHREAD */
        goto abort2;                              /* Exit this loop without swapping arrays */
#endif /* KGRAPHMAPDFLOOPTHREAD */
      }

      difntax[vertnum].partval = domnnum;         /* Anchor part is always domain part */
      difntax[vertnum].diffval = diffval;
    }

    difttax = (KgraphMapDfVertex *) difntax;      /* Swap old and new diffusion arrays          */
    difntax = (KgraphMapDfVertex *) difotax;      /* Casts to prevent IBM compiler from yelling */
    difotax = (KgraphMapDfVertex *) difttax;
abort1 : ;                                        /* If overflow occured, resume here */
#ifdef KGRAPHMAPDFLOOPTHREAD
    threadBarrier (thrdptr);

    if (loopptr->abrtval == 1)                    /* If all threads need to abort */
      break;
#endif /* KGRAPHMAPDFLOOPTHREAD */
  }
abort2 : ;

  for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) /* Set new part distribution of local vertices */
    parttax[vertnum] = difntax[vertnum].partval;

  memFree (sorttab);                              /* Free local array */

  return (0);
}

