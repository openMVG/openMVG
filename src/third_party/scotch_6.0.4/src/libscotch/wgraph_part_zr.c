/*
**  The defines and includes.
*/

#define WGRAPH_PART_ZR

#include "module.h"
#include "common.h"
#include "graph.h"
#include "wgraph.h"
#include "wgraph_part_zr.h"

/*****************************/
/*                           */
/* This is the main routine. */
/*                           */
/*****************************/

/* This routine moves all of the graph vertices
** to the first part of the partition.
** It returns:
** - 0   : if the bipartitioning could be computed.
** - !0  : on error.
*/

int
wgraphPartZr (
Wgraph * const              grafptr)              /*+ Active graph +*/
{
  
  if (grafptr->compload[0] != grafptr->s.velosum) /* If not all vertices already in part zero */
    wgraphZero (grafptr);

  return (0);
}
