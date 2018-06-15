/* Copyright 2004,2007-2011,2013,2014 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch.h                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the generic target architecture     **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 dec 1992     **/
/**                                 to   : 24 mar 1993     **/
/**                # Version 1.2  : from : 04 feb 1994     **/
/**                                 to   : 11 feb 1994     **/
/**                # Version 1.3  : from : 20 apr 1994     **/
/**                                 to   : 20 apr 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 12 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to   : 30 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     08 sep 1995     **/
/**                # Version 3.1  : from : 02 may 1996     **/
/**                                 to     20 jul 1996     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to     13 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     07 oct 1998     **/
/**                # Version 3.4  : from : 08 nov 2001     **/
/**                                 to     08 nov 2001     **/
/**                # Version 4.0  : from : 01 jan 2002     **/
/**                                 to     07 dec 2004     **/
/**                # Version 5.1  : from : 11 dec 2007     **/
/**                                 to     11 aug 2010     **/
/**                # Version 6.0  : from : 14 fev 2011     **/
/**                                 to     01 aug 2014     **/
/**                                                        **/
/************************************************************/

#define ARCH_H

/*
** The defines.
*/

/*+ Architecture flags. +*/

#define ARCHNONE                    0x0000        /*+ No options set                      +*/
#define ARCHPART                    0x0001        /*+ Architecture without external gains +*/
#define ARCHVAR                     0x0002        /*+ Variable-sized architecture         +*/

/*
**  The type and structure definitions.
*/

typedef INT Anum;                                 /*+ Generic integer +*/

#define ANUMMAX                     INTVALMAX
#define ANUMSTRING                  INTSTRING
#define ANUM_MPI                    COMM_INT      /*+ MPI type for Gnum is MPI type for INT +*/

/*+ The domain number type. +*/

typedef Anum ArchDomNum;                          /*+ Domain number +*/

#define ARCHDOMNOTTERM              ((ArchDomNum) ~0) /*+ Not-terminal number +*/

/*+ The sub-includes for structure size computations. +*/

#define ARCH_NOPROTO 
#include "arch_cmplt.h"
#include "arch_cmpltw.h"
#include "arch_deco.h"
#include "arch_dist.h"
#include "arch_hcub.h"
#include "arch_tleaf.h"
#include "arch_mesh.h"
#include "arch_torus.h"
#include "arch_vcmplt.h"
#include "arch_vhcub.h"
#undef ARCH_NOPROTO

/*+ The architecture class type. +*/

typedef struct ArchClass_ {
  char *                    archname;             /*+ Architecture name                 +*/
  int                       flagval;              /*+ Architecture flags of the class   +*/
  int                    (* archLoad)   ();       /*+ Architecture loading function     +*/
  int                    (* archSave)   ();       /*+ Architecture saving function      +*/
  int                    (* archFree)   ();       /*+ Architecture freeing function     +*/
  ArchDomNum             (* domNum)     ();       /*+ Domain labeling function          +*/
  int                    (* domTerm)    ();       /*+ Terminal domain building function +*/
  Anum                   (* domSize)    ();       /*+ Domain size function              +*/
  Anum                   (* domWght)    ();       /*+ Domain weight function            +*/
  Anum                   (* domDist)    ();       /*+ Distance computation function     +*/
  int                    (* domFrst)    ();       /*+ Compute biggest domain            +*/
  int                    (* domLoad)    ();       /*+ Domain loading routine            +*/
  int                    (* domSave)    ();       /*+ Domain saving routine             +*/
  int                    (* domBipart)  ();       /*+ Domain bipartitioning routine     +*/
  int                    (* domIncl)    ();       /*+ Domain inclusion routine          +*/
#ifdef SCOTCH_PTSCOTCH
  int                    (* domMpiType) ();       /*+ Domain MPI type building routine  +*/
#endif /* SCOTCH_PTSCOTCH */
  int                       domsizeof;            /*+ Size in bytes of domain data      +*/
} ArchClass;

/*+ The architecture union type. +*/

typedef union {                                   /*+ Architecture data                          +*/
  ArchCmplt                 cmplt;                /*+ Complete graph architecture                +*/
  ArchCmpltw                cmpltw;               /*+ Weighted complete graph architecture       +*/
  ArchDeco                  deco;                 /*+ Decomposition-described architecture       +*/
  ArchDist                  dist;                 /*+ Distance multiplicator pseudo-architecture +*/
  ArchHcub                  hcub;                 /*+ Hypercube architecture                     +*/
  ArchMesh2                 mesh2;                /*+ 2D-mesh architecture                       +*/
  ArchMesh3                 mesh3;                /*+ 3D-mesh architecture                       +*/
  ArchTleaf                 tleaf;                /*+ Tree-leaf architecture                     +*/
  ArchTorusX                torusx;               /*+ xD-torus architecture (includes 2D and 3D) +*/
  ArchVcmplt                vcmplt;               /*+ Variable-sized complete graph architecture +*/
  ArchVhcub                 vhcub;                /*+ Variable-sized hypercube architecture      +*/
} ArchDummy;

/*+ The architecture type. +*/

typedef struct Arch_ {
  const ArchClass *         class;                /*+ Pointer to architecture class         +*/
  int                       flagval;              /*+ (Possibly updated) architecture flags +*/
  ArchDummy                 data;                 /*+ Architecture data                     +*/
} Arch;

/*+ The architecture domain union type. +*/

typedef union {                                   /*+ The domain data                           +*/
  ArchCmpltDom              cmplt;                /*+ Complete graph domain type                +*/
  ArchCmpltwDom             cmpltw;               /*+ Weighted complete graph domain type       +*/
  ArchDecoDom               deco;                 /*+ Decomposition-descripted domain type      +*/
/*+ ArchDistDom             dist;                  *+ Distance multiplicator domain is ArchDom  +*/
  ArchHcubDom               hcub;                 /*+ Hypercube domain type                     +*/
  ArchMesh2Dom              mesh2;                /*+ 2D-mesh domain type                       +*/
  ArchMesh3Dom              mesh3;                /*+ 3D-mesh domain type                       +*/
  ArchTleafDom              tleaf;                /*+ Tree-leaf domain type                     +*/
  ArchTorusXDom             torusx;               /*+ xD-torus domain type (includes 2D and 3D) +*/
  ArchVcmpltDom             vcmplt;               /*+ Variable-sized complete graph domain type +*/
  ArchVhcubDom              vhcub;                /*+ Variable-sized hypercube domain type      +*/
} ArchDomDummy;

/*+ The domain structure type. +*/

typedef struct ArchDom_ {
  ArchDomDummy              data;                 /*+ The domain data +*/
} ArchDom;

/*
**  The function prototypes.
*/

#ifndef ARCH
#define static
#endif

int                         archInit            (Arch * restrict const);
int                         archExit            (Arch * restrict const);
int                         archFree            (Arch * restrict const);
int                         archLoad            (Arch * restrict const, FILE * const);
int                         archSave            (const Arch * const, FILE * const);
char *                      archName            (const Arch * const);
const ArchClass *           archClass           (const char * const);

ArchDomNum                  archDomNum          (const Arch * const, const ArchDom * const);
int                         archDomTerm         (const Arch * const, ArchDom * const, const ArchDomNum);
Anum                        archDomSize         (const Arch * const, const ArchDom * const);
Anum                        archDomWght         (const Arch * const, const ArchDom * const);
Anum                        archDomDist         (const Arch * const, const ArchDom * const, const ArchDom * const);
int                         archDomFrst         (const Arch * const, ArchDom * const);
int                         archDomLoad         (const Arch * const, ArchDom * const, FILE * const);
int                         archDomSave         (const Arch * const, const ArchDom * const, FILE * const);
int                         archDomBipart       (const Arch * const, const ArchDom * const, ArchDom * const, ArchDom * const);
int                         archDomIncl         (const Arch * const, const ArchDom * const, const ArchDom * const);
#ifdef SCOTCH_PTSCOTCH
int                         archDomMpiType      (const Arch * const, MPI_Datatype * const);
#endif /* SCOTCH_PTSCOTCH */

#undef static

/*
**  The macro definitions.
*/

#define archDomSizeof(a)            ((a)->class->domsizeof)
#define archName(a)                 (((a)->class == NULL) ? "" : (a)->class->archname)
#define archPart(a)                 ((((a)->flagval) & ARCHPART) != 0)
#define archVar(a)                  ((((a)->flagval) & ARCHVAR) != 0)

#if ((! defined SCOTCH_DEBUG_ARCH2) || (defined ARCH))
#define archDomNum2(arch,dom)       (((ArchDomNum (*) (const void * const, const void * const)) (arch)->class->domNum) ((const void * const) &(arch)->data, (const void * const) &(dom)->data))
#define archDomTerm2(arch,dom,num)  (((int (*) (const void * const, void * const, const ArchDomNum)) (arch)->class->domTerm) ((void *) &(arch)->data, (void *) &(dom)->data, (num)))
#define archDomSize2(arch,dom)      (((Anum (*) (const void * const, const void * const)) (arch)->class->domSize) ((void *) &(arch)->data, (void *) &(dom)->data))
#define archDomWght2(arch,dom)      (((Anum (*) (const void * const, const void * const)) (arch)->class->domWght) ((void *) &(arch)->data, (void *) &(dom)->data))
#define archDomDist2(arch,dom0,dom1) (((Anum (*) (const void * const, const void * const, const void * const)) (arch)->class->domDist) ((const void *) &(arch)->data, (const void *) &(dom0)->data, (const void *) &(dom1)->data))
#define archDomFrst2(arch,dom)      (((int (*) (const void * const, void * const)) (arch)->class->domFrst) ((const void * const) &(arch)->data, (void * const) &(dom)->data))
#define archDomBipart2(arch,dom,dom0,dom1) (((int (*) (const void * const, const void * const, void * const, void * const)) (arch)->class->domBipart) ((const void * const) &(arch)->data, (const void * const) &(dom)->data, (void * const) &(dom0)->data, (void * const) &(dom1)->data))
#define archDomIncl2(arch,dom0,dom1) (((int (*) (const void * const, const void * const, void * const)) (arch)->class->domIncl) ((const void * const) &(arch)->data, (void * const) &(dom0)->data, (void * const) &(dom1)->data))
#endif
#ifndef SCOTCH_DEBUG_ARCH2
#define archDomNum                  archDomNum2
#define archDomTerm                 archDomTerm2
#define archDomSize                 archDomSize2
#define archDomWght                 archDomWght2
#define archDomDist                 archDomDist2
#define archDomFrst                 archDomFrst2
#define archDomBipart               archDomBipart2
#define archDomIncl                 archDomIncl2
#endif /* SCOTCH_DEBUG_ARCH2 */

#ifdef SCOTCH_PTSCOTCH
#define ARCHCLASSBLOCK(s,n,f)       { s, f,		   \
                                      arch##n##ArchLoad,   \
                                      arch##n##ArchSave,   \
                                      arch##n##ArchFree,   \
                                      arch##n##DomNum,     \
                                      arch##n##DomTerm,    \
                                      arch##n##DomSize,    \
                                      arch##n##DomWght,    \
                                      arch##n##DomDist,    \
                                      arch##n##DomFrst,    \
                                      arch##n##DomLoad,    \
                                      arch##n##DomSave,    \
                                      arch##n##DomBipart,  \
                                      arch##n##DomIncl,    \
                                      arch##n##DomMpiType, \
                                      sizeof (Arch##n##Dom) }
#else /* SCOTCH_PTSCOTCH */
#define ARCHCLASSBLOCK(s,n,f)       { s, f,		  \
                                      arch##n##ArchLoad,  \
                                      arch##n##ArchSave,  \
                                      arch##n##ArchFree,  \
                                      arch##n##DomNum,    \
                                      arch##n##DomTerm,   \
                                      arch##n##DomSize,   \
                                      arch##n##DomWght,   \
                                      arch##n##DomDist,   \
                                      arch##n##DomFrst,   \
                                      arch##n##DomLoad,   \
                                      arch##n##DomSave,   \
                                      arch##n##DomBipart, \
                                      arch##n##DomIncl,   \
                                      sizeof (Arch##n##Dom) }
#endif /* SCOTCH_PTSCOTCH */

#define ARCHCLASSBLOCKNULL          { NULL, ARCHNONE }

#define ARCH_H_END
