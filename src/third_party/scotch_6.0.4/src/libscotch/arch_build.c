/* Copyright 2004,2007,2008,2010,2011,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_build.c                            **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module builds a decomposition-     **/
/**                based architecture from a source graph. **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 29 may 1997     **/
/**                                 to     30 aug 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     01 oct 1998     **/
/**                # Version 3.4  : from : 30 oct 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 29 nov 2003     **/
/**                                 to     10 mar 2005     **/
/**                # Version 5.0  : from : 10 sep 2007     **/
/**                                 to     03 apr 2008     **/
/**                # Version 5.1  : from : 28 sep 2008     **/
/**                                 to     28 jun 2011     **/
/**                # Version 6.0  : from : 28 jun 2011     **/
/**                                 to     05 aug 2014     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_BUILD

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_deco.h"
#include "arch_vcmplt.h"
#include "mapping.h"
#include "bgraph.h"
#include "bgraph_bipart_st.h"
#include "arch_build.h"

/************************************/
/*                                  */
/* These routines handle job pools. */
/*                                  */
/************************************/

/* This routine frees the contents of
** the given job pool.
** It returns:
** - VOID  : in all cases.
*/

static
void
archBuildJobExit (
ArchBuildJob * const        jobtab)
{
  ArchBuildJob *          jobptr;

  jobptr = jobtab;
  do {
    graphExit (&jobptr->grafdat);
    jobptr = jobptr->joblink;
  } while (jobptr != NULL);
}

/********************************************/
/*                                          */
/* The main routine, which computes the     */
/* decomposition-based target architecture. */
/*                                          */
/********************************************/

/*
** This routine builds a target architecture from
** the given source graph and the optional vertex
** list.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archBuild (
Arch * restrict const       tgtarchptr,           /*+ Decomposition architecture to build    +*/
const Graph * const         tgtgrafptr,           /*+ Source graph modeling the architecture +*/
const VertList * const      tgtlistptr,           /*+ Subset of source graph vertices        +*/
const Strat * const         mapstrat)             /*+ Bipartitioning strategy                +*/
{
  Gnum * restrict               mapparttax;       /* Based access to mapping part array             */
  Arch                          archdat;          /* Variable-sized architecture for bipartitioning */
  ArchDom                       domsub0;          /* Temporary space for subdomain 0                */
  Gnum                          termdomnbr;       /* Number of terminal domains                     */
  Gnum                          termdommax;       /* Maximum terminal number                        */
  ArchDecoTermVert * restrict   termverttab;      /* Terminal vertex table                          */
  Anum * restrict               termdisttab;      /* Vertex distance table                          */
  ArchBuildDistElem * restrict  disttax;          /* Distance table                                 */
  ArchBuildQueuElem * restrict  queutab;          /* Distance queue table                           */
  Gnum                          queuhead;         /* Head-of-queue index                            */
  Gnum                          queutail;         /* Tail-of-queue index                            */
  Mapping                       mappdat;          /* Partial and final mapping data                 */
  ArchBuildJob * restrict       jobtab;           /* Job array                                      */
  ArchBuildJob *                joblink;          /* Linked list of jobs to process                 */
  ArchBuildJob *                joborgptr;        /* Pointer to original job and first subjob       */
  ArchBuildJob *                jobsubptr;        /* Pointer to second subjob                       */
  Bgraph                        actgrafdat;       /* Active graph for bipartitioning                */
  Gnum                          invedlosiz;       /* Size of inversed edge load array               */
  Gnum * restrict               invedlotax;       /* Inversed edge load array for cutting           */
  Gnum * restrict               actfrontab;       /* Frontier array for all bipartitionings         */
  GraphPart * restrict          actparttax;       /* Part array for all bipartitionings             */
  GraphPart                     actpartval;       /* Part value to put to subjob                    */
  Gnum                          actpartnbr;       /* Size of part value to put to subjob            */
  Gnum                          termdomnum;

  const Gnum * restrict const tgtverttax = tgtgrafptr->verttax;
  const Gnum * restrict const tgtvendtax = tgtgrafptr->vendtax;
  const Gnum * restrict const tgtedgetax = tgtgrafptr->edgetax;
  const Gnum * restrict const tgtedlotax = tgtgrafptr->edlotax;

  archInit (tgtarchptr);                          /* Initialize architecture body */
  tgtarchptr->class = archClass ("deco");         /* Set architecture class       */

  termdomnbr = (tgtlistptr != NULL) ? tgtlistptr->vnumnbr : tgtgrafptr->vertnbr;
  if (termdomnbr == 0)                            /* If nothing to do */
    return (0);

  intRandInit ();                                 /* Initialize random generator */

  invedlosiz = (tgtedlotax != NULL) ? tgtgrafptr->edgenbr : 0;
  if ((memAllocGroup ((void **) (void *)
                      &jobtab,     (size_t) (termdomnbr * sizeof (ArchBuildJob)),
                      &actfrontab, (size_t) (termdomnbr * sizeof (Gnum)),
                      &actparttax, (size_t) (termdomnbr * sizeof (GraphPart)),
                      &invedlotax, (size_t) (invedlosiz * sizeof (Gnum)), NULL) == NULL) ||
      ((mappdat.parttax = memAlloc (tgtgrafptr->vertnbr * sizeof (ArchDomNum))) == NULL) || /* Final mapping array is for all graph vertices */
      ((mappdat.domntab = memAlloc (termdomnbr * sizeof (ArchDom)))             == NULL)) {
    errorPrint ("archBuild: out of memory (1)");
    if (jobtab != NULL) {
      memFree (jobtab);
      if (mappdat.parttax != NULL)
        memFree (mappdat.parttax);
    }
    return (1);
  }
  memSet (mappdat.parttax, 0, termdomnbr * sizeof (ArchDomNum));
  actparttax      -= tgtgrafptr->baseval;
  mappdat.flagval  = MAPPINGFREEPART | MAPPINGFREEDOMN;
  mappdat.grafptr  = tgtgrafptr;
  mappdat.archptr  = &archdat;
  mappdat.parttax -= tgtgrafptr->baseval;
  mappdat.domnmax  = termdomnbr;

  archInit (&archdat);                            /* Initialize terminal architecture */
  archdat.class = archClass ("varcmplt");         /* Set architecture class           */
  archDomFrst (&archdat, &mappdat.domntab[0]);    /* Get initial domain               */
  mappdat.domnnbr = 1;

  jobtab[0].domnum = 0;                           /* All vertices mapped to first domain  */
  if ((tgtlistptr != NULL) && (tgtlistptr->vnumtab != NULL)) /* If vertex list given      */
    graphInduceList (tgtgrafptr, tgtlistptr, &jobtab[0].grafdat); /* Restrict initial job */
  else {                                          /* If no vertex list given              */
    memCpy (&jobtab[0].grafdat, tgtgrafptr, sizeof (Graph)); /* Job takes whole graph     */
    jobtab[0].grafdat.flagval &= ~GRAPHFREETABS;  /* Graph is a clone                     */
    jobtab[0].grafdat.vnumtax  = NULL;            /* Assume we have no vertex index array */
  }

  if (tgtedlotax != NULL) {                       /* If architecture graph has edge loads */
    Gnum                  vertnum;
    Gnum                  vertnnd;
    Gnum                  edlomin;
    Gnum                  edlomax;
    float                 prodval;

    const Gnum * restrict const indverttax = jobtab[0].grafdat.verttax;
    const Gnum * restrict const indvendtax = jobtab[0].grafdat.vendtax;
    const Gnum * restrict const indedlotax = jobtab[0].grafdat.edlotax; /* Point to possibly induced original edge array */

    invedlotax -= tgtgrafptr->baseval;            /* Base inversed edge load array */

    edlomin = GNUMMAX;
    edlomax = 1;
    for (vertnum = jobtab[0].grafdat.baseval, vertnnd = jobtab[0].grafdat.vertnnd; /* Handle non-compact graphs as well as compact graphs */
         vertnum < vertnnd; vertnum ++) {
      Gnum                  edgenum;
      Gnum                  edgennd;

      for (edgenum = indverttax[vertnum], edgennd = indvendtax[vertnum];
           edgenum < edgennd; edgenum ++) {
        Gnum                edloval;

        edloval = indedlotax[edgenum];
        if (edloval < edlomin)
          edlomin = edloval;
        if (edloval > edlomax)
          edlomax = edloval;
      }
    }

    prodval = (float) edlomin * (float) edlomax;

    for (vertnum = jobtab[0].grafdat.baseval;
         vertnum < vertnnd; vertnum ++) {
      Gnum                  edgenum;
      Gnum                  edgennd;

      for (edgenum = indverttax[vertnum], edgennd = indvendtax[vertnum];
           edgenum < edgennd; edgenum ++) {
        Gnum                edloval;

        edloval = indedlotax[edgenum];
        if (edloval == edlomin)
          edloval = edlomax;
        else if (edloval == edlomax)
          edloval = edlomin;
        else
          edloval = (Gnum) (prodval / (float) edloval + 0.49F);
#ifdef SCOTCH_DEBUG_ARCH2
        if ((edloval < edlomin) || (edloval > edlomax)) {
          errorPrint ("archBuild: internal error (1)");
          return     (1);
        }
#endif /* SCOTCH_DEBUG_ARCH2 */
        invedlotax[edgenum] = edloval;            /* Write inversed cost in working array */
      }
    }

    jobtab[0].grafdat.edlotax = invedlotax;       /* Replace potentially induced edge array with inversed one */
  }                                               /* Edge array will be freed along with jobtab group leader  */

  mapparttax = mappdat.parttax;

  actgrafdat.veextax = NULL;                      /* No external gain array      */
  actgrafdat.parttax = actparttax;                /* Set global auxiliary arrays */
  actgrafdat.frontab = actfrontab;
  joblink = NULL;                                 /* Initialize job list            */
  if (jobtab[0].grafdat.vertnbr > 1) {            /* If job is worth bipartitioning */
    jobtab[0].joblink = joblink;                  /* Add initial job to list        */
    joblink = &jobtab[0];
  }
  while (joblink != NULL) {                       /* For all jobs in list */
    joborgptr          = joblink;                 /* Get job              */
    joblink            = joblink->joblink;        /* Remove job from list */
    joborgptr->joblink = NULL;                    /* In case of freeing   */

    memCpy (&actgrafdat.s, &joborgptr->grafdat, sizeof (Graph));
    actgrafdat.s.flagval = joborgptr->grafdat.flagval & ~GRAPHFREETABS;
    bgraphInit2 (&actgrafdat, 1, 1, 1, 0, 0);     /* Create active graph         */
    if (bgraphBipartSt (&actgrafdat, mapstrat) != 0) { /* Perform bipartitioning */
      errorPrint       ("archBuild: internal error (2)");
      archBuildJobExit (joborgptr);
      archBuildJobExit (joblink);
      archExit         (&archdat);
      mapExit          (&mappdat);
      memFree          (jobtab);
      return           (1);
    }
    if ((actgrafdat.compsize0 == 0) ||            /* If one of the jobs is empty */
        (actgrafdat.compsize0 == actgrafdat.s.vertnbr)) {
      errorPrint       ("archBuild: strategy leads to empty domains");
      graphExit        (&actgrafdat.s);           /* Only free graph part, global arrays kept */
      archBuildJobExit (joborgptr);
      archBuildJobExit (joblink);
      archExit         (&archdat);
      mapExit          (&mappdat);
      memFree          (jobtab);
      return           (1);
    }

    archVcmpltDomBipart ((const ArchVcmplt * const) (void *) &archdat, /* Update mapping domains */
                         (const ArchVcmpltDom * const) (void *) &mappdat.domntab[joborgptr->domnum],
                         (ArchVcmpltDom * const) (void *) &domsub0,
                         (ArchVcmpltDom * const) (void *) &mappdat.domntab[mappdat.domnnbr]);
    mappdat.domntab[joborgptr->domnum] = domsub0;
    actpartval = actgrafdat.parttax[actgrafdat.s.baseval]; /* Always keep first vertex in sub0 */
    actpartnbr = (actpartval == 0) ? actgrafdat.compsize0 : (actgrafdat.s.vertnbr - actgrafdat.compsize0);
    if (actgrafdat.s.vnumtax != NULL) {           /* Update mapping fraction */
      Gnum                actvertnum;

      for (actvertnum = actgrafdat.s.baseval; actvertnum < actgrafdat.s.vertnnd; actvertnum ++) {
        if (actgrafdat.parttax[actvertnum] != actpartval)
          mappdat.parttax[actgrafdat.s.vnumtax[actvertnum]] = mappdat.domnnbr;
      }
    }
    else {
      Gnum                actvertnum;

      for (actvertnum = actgrafdat.s.baseval; actvertnum < actgrafdat.s.vertnnd; actvertnum ++) {
        if (actgrafdat.parttax[actvertnum] != actpartval)
          mappdat.parttax[actvertnum] = mappdat.domnnbr;
      }
    }

    jobsubptr = jobtab + mappdat.domnnbr;         /* Point to new subjob          */
    jobsubptr->domnum = mappdat.domnnbr ++;       /* Build subjobs                */
    actgrafdat.s.flagval = joborgptr->grafdat.flagval; /* Active is now main copy */

    if (actpartnbr < (actgrafdat.s.vertnbr - 1)) { /* If part 1 splittable */
      graphInducePart (&actgrafdat.s, actgrafdat.parttax, actgrafdat.s.vertnbr - actpartnbr,
                       1 - actpartval, &jobsubptr->grafdat);
      jobsubptr->joblink = joblink;               /* Link subjobs in list */
      joblink            = jobsubptr;
    }
    if (actpartnbr > 1) {                         /* If part 0 splittable */
      graphInducePart (&actgrafdat.s, actgrafdat.parttax, actpartnbr,
                       actpartval, &joborgptr->grafdat);
      joborgptr->joblink = joblink;               /* Link subjobs in list */
      joblink            = joborgptr;
    }
    graphExit (&actgrafdat.s);                    /* Only free graph part, global arrays kept */
  }

  memFree (jobtab);                               /* Free group leader */

  if (memAllocGroup ((void **) (void *)
                     &termverttab, (size_t) (termdomnbr                            * sizeof (ArchDecoTermVert)),
                     &termdisttab, (size_t) (((termdomnbr * (termdomnbr - 1)) / 2) * sizeof (Anum)),
                     &disttax,     (size_t) (tgtgrafptr->vertnbr                   * sizeof (ArchBuildDistElem)), 
                     &queutab,     (size_t) (tgtgrafptr->vertnbr                   * sizeof (ArchBuildQueuElem)), NULL) == NULL) {
    errorPrint ("archBuild: out of memory (2)");
    mapExit    (&mappdat);
    archExit   (&archdat);
    return     (1);
  }

  for (termdomnum = 0, termdommax = 0; termdomnum < termdomnbr; termdomnum ++) { /* Set terminal vertex array */
    Gnum                tgtvertnum;

    tgtvertnum = (tgtlistptr != NULL) ? tgtlistptr->vnumtab[termdomnum] : (termdomnum + tgtgrafptr->baseval);
    termverttab[termdomnum].labl = tgtvertnum;
    termverttab[termdomnum].wght = (tgtgrafptr->velotax != NULL) ? tgtgrafptr->velotax[tgtvertnum] : 1;
    termverttab[termdomnum].num  = archDomNum (&archdat, mapDomain (&mappdat, tgtvertnum - tgtgrafptr->baseval));
    if (termverttab[termdomnum].num > termdommax) /* Find maximum terminal number */
      termdommax = termverttab[termdomnum].num;
  }

  disttax -= tgtgrafptr->baseval;
  for (termdomnum = 1; termdomnum < termdomnbr; termdomnum ++) { /* For all active terminal vertices except the first */
    Gnum                termdomend;
    Gnum                tgtvertnum;

    for (tgtvertnum = tgtgrafptr->baseval; tgtvertnum < (tgtgrafptr->vertnbr + tgtgrafptr->baseval); tgtvertnum ++) {
      disttax[tgtvertnum].queued  = 0;            /* Vertex not queued       */
      disttax[tgtvertnum].distval = INTVALMAX;    /* Assume maximum distance */
    }

    queuhead =                                    /* Reset the queue */
    queutail = 0;
    tgtvertnum = termverttab[termdomnum].labl;
    queutab[queutail].vertnum    = tgtvertnum;    /* Insert root vertex */
    queutab[queutail ++].distval = 0;
    disttax[tgtvertnum].queued  = 1;              /* Mark vertex as queued */
    disttax[tgtvertnum].distval = 0;

    while (queuhead < queutail) {                 /* As long as there are vertices in queue */
      Gnum                vertnum;                /* Number of current vertex               */
      Gnum                vertdist;               /* Current distance value                 */
      Gnum                edgenum;

      vertnum  = queutab[queuhead].vertnum;       /* Retrieve vertex from queue */
      vertdist = queutab[queuhead ++].distval;

      for (edgenum = tgtverttax[vertnum];         /* For all vertex edges */
           edgenum < tgtvendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = tgtedgetax[edgenum];
        if (disttax[vertend].queued == 0) {       /* If end vertex not queued */
          queutab[queutail].vertnum    = vertend; /* Queue the vertex         */
          queutab[queutail ++].distval =
          disttax[vertend].distval = vertdist + ((tgtedlotax != NULL) ? tgtedlotax[edgenum] : 1);
          disttax[vertend].queued  = 1;           /* Mark vertex as queued */
        }
      }
    }

    for (termdomend = 0; termdomend < termdomnum; termdomend ++) /* For all previous terminal domains */
      termdisttab[((termdomnum * (termdomnum - 1)) / 2) + termdomend] = /* Retrieve distance          */
        disttax[termverttab[termdomend].labl].distval;
  }

  archDecoArchBuild ((ArchDeco *) (void *) &tgtarchptr->data, termdomnbr, termdommax, termverttab, termdisttab);

  memFree  (termverttab);                         /* Free group leader */
  mapExit  (&mappdat);
  archExit (&archdat);

  return (0);
}
