/* Copyright 2004,2007-2015 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : module.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This is the global configuration file   **/
/**                for the whole libSCOTCH library module. **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 22 jun 1998     **/
/**                                 to     13 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     03 oct 1998     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 12 dec 2001     **/
/**                                 to     24 nov 2005     **/
/**                # Version 5.0  : from : 24 feb 2007     **/
/**                                 to     24 jul 2007     **/
/**                # Version 5.1  : from : 25 oct 2007     **/
/**                                 to     20 feb 2011     **/
/**                # Version 6.0  : from : 12 sep 2008     **/
/**                                 to     01 mar 2015     **/
/**                                                        **/
/************************************************************/

#define MODULE_H

/*
** Version string.
*/

#define SCOTCH_VERSION_STRING       SCOTCH_VERSION_STRING2(SCOTCH_VERSION) "." SCOTCH_VERSION_STRING2(SCOTCH_RELEASE) "." SCOTCH_VERSION_STRING2(SCOTCH_PATCHLEVEL)
#define SCOTCH_VERSION_STRING2(x)   SCOTCH_VERSION_STRING3(x)
#define SCOTCH_VERSION_STRING3(x)   #x

/*
** Debug values.
*/

#ifdef SCOTCH_DEBUG_ALL
#ifndef SCOTCH_DEBUG
#define SCOTCH_DEBUG
#endif /* SCOTCH_DEBUG */

#define COMMON_DEBUG
#define SCOTCH_DEBUG_ARCH2
#define SCOTCH_DEBUG_FIBO2
#define SCOTCH_DEBUG_GAIN2
#define SCOTCH_DEBUG_PARSER2
#define SCOTCH_DEBUG_BDGRAPH2
#define SCOTCH_DEBUG_BGRAPH2
#define SCOTCH_DEBUG_DGRAPH2
#define SCOTCH_DEBUG_DMAP2
#define SCOTCH_DEBUG_DORDER2
#define SCOTCH_DEBUG_GEOM2
#define SCOTCH_DEBUG_GRAPH2
#define SCOTCH_DEBUG_HDGRAPH2
#define SCOTCH_DEBUG_HGRAPH2
#define SCOTCH_DEBUG_HMESH2
#define SCOTCH_DEBUG_KDGRAPH2
#define SCOTCH_DEBUG_KDMAP2
#define SCOTCH_DEBUG_KGRAPH2
#define SCOTCH_DEBUG_LIBRARY2
#define SCOTCH_DEBUG_MAP2
#define SCOTCH_DEBUG_MESH2
#define SCOTCH_DEBUG_ORDER2
#define SCOTCH_DEBUG_PARSER2
#define SCOTCH_DEBUG_VDGRAPH2
#define SCOTCH_DEBUG_VGRAPH2
#define SCOTCH_DEBUG_VMESH2
#define SCOTCH_DEBUG_WGRAPH2
#endif /* SCOTCH_DEBUG_ALL */

#ifdef SCOTCH_DEBUG
#define SCOTCH_DEBUG_ARCH1
#define SCOTCH_DEBUG_FIBO1
#define SCOTCH_DEBUG_GAIN1
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_BDGRAPH1
#define SCOTCH_DEBUG_BGRAPH1
#define SCOTCH_DEBUG_DGRAPH1
#define SCOTCH_DEBUG_DMAP1
#define SCOTCH_DEBUG_DORDER1
#define SCOTCH_DEBUG_GEOM1
#define SCOTCH_DEBUG_GRAPH1
#define SCOTCH_DEBUG_HDGRAPH1
#define SCOTCH_DEBUG_HGRAPH1
#define SCOTCH_DEBUG_HMESH1
#define SCOTCH_DEBUG_KDGRAPH1
#define SCOTCH_DEBUG_KDMAP1
#define SCOTCH_DEBUG_KGRAPH1
#define SCOTCH_DEBUG_LIBRARY1
#define SCOTCH_DEBUG_MAP1
#define SCOTCH_DEBUG_MESH1
#define SCOTCH_DEBUG_ORDER1
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_VDGRAPH1
#define SCOTCH_DEBUG_VGRAPH1
#define SCOTCH_DEBUG_VMESH1
#define SCOTCH_DEBUG_WGRAPH1
#endif /* SCOTCH_DEBUG */

/*
** Function renaming.
*/

#if ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME))
#define memCur                      SCOTCH_memCur
#define memMax                      SCOTCH_memMax

#define clockGet                    _SCOTCHclockGet

#define commonStubDummy             _SCOTCHcommonStubDummy

#define errorPrint                  SCOTCH_errorPrint
#define errorPrintW                 SCOTCH_errorPrintW
#define errorProg                   SCOTCH_errorProg

#define fileBlockClose              _SCOTCHfileBlockClose
#define fileBlockOpen               _SCOTCHfileBlockOpen
#define fileBlockOpenDist           _SCOTCHfileBlockOpenDist
#define fileCompress                _SCOTCHfileCompress
#define fileCompressType            _SCOTCHfileCompressType
#define fileUncompress              _SCOTCHfileUncompress
#define fileUncompressType          _SCOTCHfileUncompressType
#define fileNameDistExpand          _SCOTCHfileNameDistExpand 

#define intLoad                     _SCOTCHintLoad
#define intSave                     _SCOTCHintSave
#define intAscn                     _SCOTCHintAscn
#define intGcd                      _SCOTCHintGcd
#define intPerm                     _SCOTCHintPerm
#define intRandInit                 _SCOTCHintRandInit
#define intRandProc                 _SCOTCHintRandProc
#define intRandReset                _SCOTCHintRandReset
#define intRandSeed                 _SCOTCHintRandSeed
#ifndef COMMON_RANDOM_SYSTEM
#define intRandVal                  _SCOTCHintRandVal
#endif /* COMMON_RANDOM_SYSTEM */
#define intSort1asc1                _SCOTCHintSort1asc1
#define intSort2asc1                _SCOTCHintSort2asc1
#define intSort2asc2                _SCOTCHintSort2asc2
#define intSort3asc1                _SCOTCHintSort3asc1
#define intSort3asc2                _SCOTCHintSort3asc2

#define memAllocGroup               _SCOTCHmemAllocGroup
#define memAllocRecord              _SCOTCHmemAllocRecord
#define memCheck                    _SCOTCHmemCheck
#define memCheckExists              _SCOTCHmemCheckExists
#define memCheckSize                _SCOTCHmemCheckSize
#define memCheckToggle              _SCOTCHmemCheckToggle
#define memCheckWatch               _SCOTCHmemCheckWatch
#define memFreeRecord               _SCOTCHmemFreeRecord
#define memReallocGroup             _SCOTCHmemReallocGroup
#define memReallocRecord            _SCOTCHmemReallocRecord
#define memOffset                   _SCOTCHmemOffset

#define stringSubst                 _SCOTCHstringSubst

#define usagePrint                  _SCOTCHusagePrint
#endif /* ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME)) */

#ifdef SCOTCH_RENAME
#define archInit                    _SCOTCHarchInit
#define archExit                    _SCOTCHarchExit
#define archFree                    _SCOTCHarchFree
#define archLoad                    _SCOTCHarchLoad
#define archSave                    _SCOTCHarchSave
/* #define archName                 _SCOTCHarchName Already a macro */
#define archClass                   _SCOTCHarchClass
#define archClassTab                _SCOTCHarchClassTab
#define archDomLoad                 _SCOTCHarchDomLoad
#define archDomSave                 _SCOTCHarchDomSave
#ifdef SCOTCH_DEBUG_ARCH2                         /* If already redefined */
#define archDomNum                  _SCOTCHarchDomNum
#define archDomDist                 _SCOTCHarchDomDist
#define archDomFrst                 _SCOTCHarchDomFrst
#define archDomIncl                 _SCOTCHarchDomIncl
#define archDomSize                 _SCOTCHarchDomSize
#define archDomTerm                 _SCOTCHarchDomTerm
#define archDomWght                 _SCOTCHarchDomWght
#define archDomBipart               _SCOTCHarchDomBipart
#endif /* SCOTCH_DEBUG_ARCH2 */
#define archDomMpiType              _SCOTCHarchDomMpiType
#define archBuild                   _SCOTCHarchBuild
#define archCmpltArchLoad           _SCOTCHarchCmpltArchLoad
#define archCmpltArchSave           _SCOTCHarchCmpltArchSave
#define archCmpltDomNum             _SCOTCHarchCmpltDomNum
#define archCmpltDomTerm            _SCOTCHarchCmpltDomTerm
#define archCmpltDomSize            _SCOTCHarchCmpltDomSize
/* #define archCmpltDomWght            _SCOTCHarchCmpltDomWght Already a macro */
#define archCmpltDomDist            _SCOTCHarchCmpltDomDist
#define archCmpltDomFrst            _SCOTCHarchCmpltDomFrst
#define archCmpltDomIncl            _SCOTCHarchCmpltDomIncl
#define archCmpltDomLoad            _SCOTCHarchCmpltDomLoad
#define archCmpltDomSave            _SCOTCHarchCmpltDomSave
#define archCmpltDomBipart          _SCOTCHarchCmpltDomBipart
#define archCmpltDomMpiType         _SCOTCHarchCmpltDomMpiType
#define archCmpltwArchBuild         _SCOTCHarchCmpltwArchBuild
#define archCmpltwArchFree          _SCOTCHarchCmpltwArchFree
#define archCmpltwArchLoad          _SCOTCHarchCmpltwArchLoad
#define archCmpltwArchSave          _SCOTCHarchCmpltwArchSave
#define archCmpltwDomNum            _SCOTCHarchCmpltwDomNum
#define archCmpltwDomTerm           _SCOTCHarchCmpltwDomTerm
#define archCmpltwDomSize           _SCOTCHarchCmpltwDomSize
#define archCmpltwDomWght           _SCOTCHarchCmpltwDomWght
#define archCmpltwDomDist           _SCOTCHarchCmpltwDomDist
#define archCmpltwDomFrst           _SCOTCHarchCmpltwDomFrst
#define archCmpltwDomIncl           _SCOTCHarchCmpltwDomIncl
#define archCmpltwDomLoad           _SCOTCHarchCmpltwDomLoad
#define archCmpltwDomSave           _SCOTCHarchCmpltwDomSave
#define archCmpltwDomBipart         _SCOTCHarchCmpltwDomBipart
#define archCmpltwDomMpiType        _SCOTCHarchCmpltwDomMpiType
#define archDecoArchBuild           _SCOTCHarchDecoArchBuild
#define archDecoArchFree            _SCOTCHarchDecoArchFree
#define archDecoArchLoad            _SCOTCHarchDecoArchLoad
#define archDecoArchSave            _SCOTCHarchDecoArchSave
#define archDecoDomNum              _SCOTCHarchDecoDomNum
#define archDecoDomTerm             _SCOTCHarchDecoDomTerm
#define archDecoDomSize             _SCOTCHarchDecoDomSize
#define archDecoDomWght             _SCOTCHarchDecoDomWght
#define archDecoDomDist             _SCOTCHarchDecoDomDist
#define archDecoDomFrst             _SCOTCHarchDecoDomFrst
#define archDecoDomIncl             _SCOTCHarchDecoDomIncl
#define archDecoDomLoad             _SCOTCHarchDecoDomLoad
#define archDecoDomSave             _SCOTCHarchDecoDomSave
#define archDecoDomBipart           _SCOTCHarchDecoDomBipart
#define archDecoDomMpiType          _SCOTCHarchDecoDomMpiType
#define archDistArchLoad            _SCOTCHarchDistArchLoad
#define archDistArchSave            _SCOTCHarchDistArchSave
#define archDistArchBuild           _SCOTCHarchDistArchBuild
#define archDistDomNum              _SCOTCHarchDistDomNum
#define archDistDomTerm             _SCOTCHarchDistDomTerm
#define archDistDomSize             _SCOTCHarchDistDomSize
#define archDistDomWght             _SCOTCHarchDistDomWght
#define archDistDomDist             _SCOTCHarchDistDomDist
#define archDistDomFrst             _SCOTCHarchDistDomFrst
#define archDistDomIncl             _SCOTCHarchDistDomIncl
#define archDistDomLoad             _SCOTCHarchDistDomLoad
#define archDistDomSave             _SCOTCHarchDistDomSave
#define archDistDomBipart           _SCOTCHarchDistDomBipart
#define archDistDomMpiType          _SCOTCHarchDistDomMpiType
#define archHcubArchLoad            _SCOTCHarchHcubArchLoad
#define archHcubArchSave            _SCOTCHarchHcubArchSave
#define archHcubDomNum              _SCOTCHarchHcubDomNum
#define archHcubDomTerm             _SCOTCHarchHcubDomTerm
#define archHcubDomSize             _SCOTCHarchHcubDomSize
/* #define archHcubDomWght             _SCOTCHarchHcubDomWght Already a macro */
#define archHcubDomDist             _SCOTCHarchHcubDomDist
#define archHcubDomFrst             _SCOTCHarchHcubDomFrst
#define archHcubDomIncl             _SCOTCHarchHcubDomIncl
#define archHcubDomLoad             _SCOTCHarchHcubDomLoad
#define archHcubDomSave             _SCOTCHarchHcubDomSave
#define archHcubDomBipart           _SCOTCHarchHcubDomBipart
#define archHcubDomMpiType          _SCOTCHarchHcubDomMpiType
#define archLtleafArchLoad          _SCOTCHarchLtleafArchLoad
#define archLtleafArchSave          _SCOTCHarchLtleafArchSave
#define archLtleafDomNum            _SCOTCHarchLtleafDomNum
#define archLtleafDomTerm           _SCOTCHarchLtleafDomTerm
#define archTleafArchLoad           _SCOTCHarchTleafArchLoad
#define archTleafArchFree           _SCOTCHarchTleafArchFree
#define archTleafArchSave           _SCOTCHarchTleafArchSave
#define archTleafDomNum             _SCOTCHarchTleafDomNum
#define archTleafDomTerm            _SCOTCHarchTleafDomTerm
#define archTleafDomSize            _SCOTCHarchTleafDomSize
/* #define archTleafDomWght            _SCOTCHarchTleafDomWght Already a macro */
#define archTleafDomDist            _SCOTCHarchTleafDomDist
#define archTleafDomFrst            _SCOTCHarchTleafDomFrst
#define archTleafDomIncl            _SCOTCHarchTleafDomIncl
#define archTleafDomLoad            _SCOTCHarchTleafDomLoad
#define archTleafDomSave            _SCOTCHarchTleafDomSave
#define archTleafDomBipart          _SCOTCHarchTleafDomBipart
#define archTleafDomMpiType         _SCOTCHarchTleafDomMpiType
#define archMesh2ArchLoad           _SCOTCHarchMesh2ArchLoad
#define archMesh2ArchSave           _SCOTCHarchMesh2ArchSave
#define archMesh2DomNum             _SCOTCHarchMesh2DomNum
#define archMesh2DomTerm            _SCOTCHarchMesh2DomTerm
#define archMesh2DomSize            _SCOTCHarchMesh2DomSize
/* #define archMesh2DomWght            _SCOTCHarchMesh2DomWght Already a macro */
#define archMesh2DomDist            _SCOTCHarchMesh2DomDist
#define archMesh2DomFrst            _SCOTCHarchMesh2DomFrst
#define archMesh2DomIncl            _SCOTCHarchMesh2DomIncl
#define archMesh2DomLoad            _SCOTCHarchMesh2DomLoad
#define archMesh2DomSave            _SCOTCHarchMesh2DomSave
#define archMesh2DomBipart          _SCOTCHarchMesh2DomBipart
#define archMesh2DomBipartO         _SCOTCHarchMesh2DomBipartO
#define archMesh2DomBipartU         _SCOTCHarchMesh2DomBipartU
#define archMesh2DomMpiType         _SCOTCHarchMesh2DomMpiType
#define archMesh3ArchLoad           _SCOTCHarchMesh3ArchLoad
#define archMesh3ArchSave           _SCOTCHarchMesh3ArchSave
#define archMesh3DomNum             _SCOTCHarchMesh3DomNum
#define archMesh3DomTerm            _SCOTCHarchMesh3DomTerm
#define archMesh3DomSize            _SCOTCHarchMesh3DomSize
/* #define archMesh3DomWght            _SCOTCHarchMesh3DomWght Already a macro */
#define archMesh3DomDist            _SCOTCHarchMesh3DomDist
#define archMesh3DomFrst            _SCOTCHarchMesh3DomFrst
#define archMesh3DomIncl            _SCOTCHarchMesh3DomIncl
#define archMesh3DomLoad            _SCOTCHarchMesh3DomLoad
#define archMesh3DomSave            _SCOTCHarchMesh3DomSave
#define archMesh3DomBipart          _SCOTCHarchMesh3DomBipart
#define archMesh3DomMpiType         _SCOTCHarchMesh3DomMpiType
#define archTermArchLoad            _SCOTCHarchTermArchLoad
#define archTermArchSave            _SCOTCHarchTermArchSave
#define archTermDomNum              _SCOTCHarchTermDomNum
#define archTermDomTerm             _SCOTCHarchTermDomTerm
#define archTermDomSize             _SCOTCHarchTermDomSize
/* #define archTermDomWght             _SCOTCHarchTermDomWght Already a macro */
#define archTermDomDist             _SCOTCHarchTermDomDist
#define archTermDomFrst             _SCOTCHarchTermDomFrst
#define archTermDomIncl             _SCOTCHarchTermDomIncl
#define archTermDomLoad             _SCOTCHarchTermDomLoad
#define archTermDomSave             _SCOTCHarchTermDomSave
#define archTermDomBipart           _SCOTCHarchTermDomBipart
#define archTermDomMpiType          _SCOTCHarchTermDomMpiType
#define archTorus2ArchLoad          _SCOTCHarchTorus2ArchLoad
#define archTorus2ArchSave          _SCOTCHarchTorus2ArchSave
#define archTorus2DomNum            _SCOTCHarchTorus2DomNum
#define archTorus2DomTerm           _SCOTCHarchTorus2DomTerm
#define archTorus2DomSize           _SCOTCHarchTorus2DomSize
/* #define archTorus2DomWght           _SCOTCHarchTorus2DomWght Already a macro */
#define archTorus2DomDist           _SCOTCHarchTorus2DomDist
/* #define archTorus2DomFrst           _SCOTCHarchTorus2DomFrst Already a macro */
#define archTorus2DomIncl           _SCOTCHarchTorus2DomIncl
#define archTorus2DomBipart         _SCOTCHarchTorus2DomBipart
/* #define archTorus2DomLoad           _SCOTCHarchTorus2DomLoad Already a macro */
/* #define archTorus2DomSave           _SCOTCHarchTorus2DomSave Already a macro */
#define archTorus2DomBipart         _SCOTCHarchTorus2DomBipart
/* #define archTorus2DomMpiType        _SCOTCHarchTorus2DomMpiTypeA lready a macro */
#define archTorus3ArchLoad          _SCOTCHarchTorus3ArchLoad
#define archTorus3ArchSave          _SCOTCHarchTorus3ArchSave
#define archTorus3DomNum            _SCOTCHarchTorus3DomNum
#define archTorus3DomTerm           _SCOTCHarchTorus3DomTerm
#define archTorus3DomSize           _SCOTCHarchTorus3DomSize
/* #define archTorus3DomWght           _SCOTCHarchTorus3DomWght Already a macro */
#define archTorus3DomDist           _SCOTCHarchTorus3DomDist
/* #define archTorus3DomFrst           _SCOTCHarchTorus3DomFrst Already a macro */
#define archTorus3DomIncl           _SCOTCHarchTorus3DomIncl
/* #define archTorus3DomLoad           _SCOTCHarchTorus3DomLoad Already a macro */
/* #define archTorus3DomSave           _SCOTCHarchTorus3DomSave Already a macro */
#define archTorus3DomBipart         _SCOTCHarchTorus3DomBipart
/* #define archTorus3DomMpiType        _SCOTCHarchTorus3DomMpiType Already a macro */
#define archTorusXArchLoad          _SCOTCHarchTorusXArchLoad
#define archTorusXArchSave          _SCOTCHarchTorusXArchSave
#define archTorusXDomNum            _SCOTCHarchTorusXDomNum
#define archTorusXDomTerm           _SCOTCHarchTorusXDomTerm
#define archTorusXDomSize           _SCOTCHarchTorusXDomSize
/* #define archTorusXDomWght           _SCOTCHarchTorusXDomWght Already a macro */
#define archTorusXDomDist           _SCOTCHarchTorusXDomDist
#define archTorusXDomFrst           _SCOTCHarchTorusXDomFrst
#define archTorusXDomIncl           _SCOTCHarchTorusXDomIncl
#define archTorusXDomLoad           _SCOTCHarchTorusXDomLoad
#define archTorusXDomSave           _SCOTCHarchTorusXDomSave
#define archTorusXDomBipart         _SCOTCHarchTorusXDomBipart
#define archTorusXDomMpiType        _SCOTCHarchTorusXDomMpiType
/* #define archVcmpltArchLoad          _SCOTCHarchVcmpltArchLoad Already a macro */
/* #define archVcmpltArchSave          _SCOTCHarchVcmpltArchSave Already a macro */
#define archVcmpltDomNum            _SCOTCHarchVcmpltDomNum
#define archVcmpltDomTerm           _SCOTCHarchVcmpltDomTerm
#define archVcmpltDomSize           _SCOTCHarchVcmpltDomSize
/* #define archVcmpltDomWght           _SCOTCHarchVcmpltDomWght Already a macro */
#define archVcmpltDomDist           _SCOTCHarchVcmpltDomDist
#define archVcmpltDomFrst           _SCOTCHarchVcmpltDomFrst
#define archVcmpltDomIncl           _SCOTCHarchVcmpltDomIncl
#define archVcmpltDomBipart         _SCOTCHarchVcmpltDomBipart
#define archVcmpltDomLoad           _SCOTCHarchVcmpltDomLoad
#define archVcmpltDomSave           _SCOTCHarchVcmpltDomSave
#define archVcmpltDomBipart         _SCOTCHarchVcmpltDomBipart
#define archVcmpltDomMpiType        _SCOTCHarchVcmpltDomMpiType
/* #define archVhcubArchLoad           _SCOTCHarchVhcubArchLoad Already a macro */
/* #define archVhcubArchSave           _SCOTCHarchVhcubArchSave Already a macro */
#define archVhcubDomNum             _SCOTCHarchVhcubDomNum
#define archVhcubDomTerm            _SCOTCHarchVhcubDomTerm
#define archVhcubDomSize            _SCOTCHarchVhcubDomSize
/* #define archVhcubDomWght            _SCOTCHarchVhcubDomWght Already a macro */
#define archVhcubDomDist            _SCOTCHarchVhcubDomDist
#define archVhcubDomFrst            _SCOTCHarchVhcubDomFrst
#define archVhcubDomIncl            _SCOTCHarchVhcubDomIncl
#define archVhcubDomLoad            _SCOTCHarchVhcubDomLoad
#define archVhcubDomSave            _SCOTCHarchVhcubDomSave
#define archVhcubDomBipart          _SCOTCHarchVhcubDomBipart
#define archVhcubDomMpiType         _SCOTCHarchVhcubDomMpiType

#define bdgraphInit                 _SCOTCHbdgraphInit
#define bdgraphInit2                _SCOTCHbdgraphInit2
#define bdgraphExit                 _SCOTCHbdgraphExit
#define bdgraphZero                 _SCOTCHbdgraphZero
#define bdgraphbipartststratab      _SCOTCHbdgraphbipartststratab
#define bdgraphCheck                _SCOTCHbdgraphCheck
#define bdgraphGatherAll            _SCOTCHbdgraphGatherAll
#define bdgraphBipartBd             _SCOTCHbdgraphBipartBd
#define bdgraphBipartDf             _SCOTCHbdgraphBipartDf
#define bdgraphBipartEx             _SCOTCHbdgraphBipartEx
#define bdgraphBipartMl             _SCOTCHbdgraphBipartMl
#define bdgraphBipartSq             _SCOTCHbdgraphBipartSq
#define bdgraphBipartSt             _SCOTCHbdgraphBipartSt
#define bdgraphBipartZr             _SCOTCHbdgraphBipartZr
#define bdgraphStoreInit            _SCOTCHbdgraphStoreInit
#define bdgraphStoreExit            _SCOTCHbdgraphStoreExit
#define bdgraphStoreSave            _SCOTCHbdgraphStoreSave
#define bdgraphStoreUpdt            _SCOTCHbdgraphStoreUpdt

#define bgraphbipartststratab       _SCOTCHbgraphbipartststratab
#define bgraphInit                  _SCOTCHbgraphInit
#define bgraphInit2                 _SCOTCHbgraphInit2
#define bgraphInit3                 _SCOTCHbgraphInit3
#define bgraphInit4                 _SCOTCHbgraphInit4
#define bgraphInit5                 _SCOTCHbgraphInit5
#define bgraphExit                  _SCOTCHbgraphExit
#define bgraphCheck                 _SCOTCHbgraphCheck
#define bgraphSwal                  _SCOTCHbgraphSwal
#define bgraphZero                  _SCOTCHbgraphZero
#define bgraphBipartBd              _SCOTCHbgraphBipartBd
#define bgraphBipartDf              _SCOTCHbgraphBipartDf
#define bgraphBipartDf2             _SCOTCHbgraphBipartDf2
#define bgraphBipartDfJoin          _SCOTCHbgraphBipartDfJoin
#define bgraphBipartEx              _SCOTCHbgraphBipartEx
#define bgraphBipartFm              _SCOTCHbgraphBipartFm
#define bgraphBipartGg              _SCOTCHbgraphBipartGg
#define bgraphBipartGp              _SCOTCHbgraphBipartGp
#define bgraphBipartMl              _SCOTCHbgraphBipartMl
#define bgraphBipartSt              _SCOTCHbgraphBipartSt
#define bgraphBipartZr              _SCOTCHbgraphBipartZr
#define bgraphStoreInit             _SCOTCHbgraphStoreInit
#define bgraphStoreExit             _SCOTCHbgraphStoreExit
#define bgraphStoreSave             _SCOTCHbgraphStoreSave
#define bgraphStoreUpdt             _SCOTCHbgraphStoreUpdt

#if ((defined INTSIZE64) || (defined COMM))
#define commAllgatherv              _SCOTCHcommAllgatherv
#define commGatherv                 _SCOTCHcommGatherv
#define commScatterv                _SCOTCHcommScatterv
#endif /* ((defined INTSIZE64) || (defined COMM)) */

#define dgraphAllreduceMaxSum2      _SCOTCHdgraphAllreduceMaxSum2
#define dgraphBuild                 _SCOTCHdgraphBuild
#define dgraphBuild2                _SCOTCHdgraphBuild2
#define dgraphBuild3                _SCOTCHdgraphBuild3
#define dgraphBuild4                _SCOTCHdgraphBuild4
#define dgraphBuildGrid3D           _SCOTCHdgraphBuildGrid3D
#define dgraphBuildHcub             _SCOTCHdgraphBuildHcub
#define dgraphCheck                 _SCOTCHdgraphCheck
#define dgraphBand                  _SCOTCHdgraphBand
#define dgraphBandColl              _SCOTCHdgraphBandColl
#define dgraphBandPtop              _SCOTCHdgraphBandPtop
#define dgraphCoarsen               _SCOTCHdgraphCoarsen
#define dgraphExit                  _SCOTCHdgraphExit
#define dgraphFold                  _SCOTCHdgraphFold
#define dgraphFold2                 _SCOTCHdgraphFold2
#define dgraphFoldComm              _SCOTCHdgraphFoldComm
#define dgraphFoldDup               _SCOTCHdgraphFoldDup
#define dgraphFree                  _SCOTCHdgraphFree
#define dgraphGather                _SCOTCHdgraphGather
#define dgraphGatherAll             _SCOTCHdgraphGatherAll
#define dgraphGatherAll2            _SCOTCHdgraphGatherAll2
/* #define dgraphGhst                  _SCOTCHdgraphGhst Already a macro        */
/* #define dgraphGhstReplace           _SCOTCHdgraphGhstReplace Already a macro */
#define dgraphGhst2                 _SCOTCHdgraphGhst2
#define dgraphGrow                  _SCOTCHdgraphGrow /* Used before macro replacement */
#define dgraphGrowColl              _SCOTCHdgraphGrowColl
#define dgraphGrowPtop              _SCOTCHdgraphGrowPtop
#define dgraphHaloSync              _SCOTCHdgraphHaloSync
#define dgraphHaloAsync             _SCOTCHdgraphHaloAsync
#define dgraphHaloWait              _SCOTCHdgraphHaloWait
#define dgraphHaloCheck             _SCOTCHdgraphHaloCheck
#define dgraphInduceList            _SCOTCHdgraphInduceList
#define dgraphInducePart            _SCOTCHdgraphInducePart
#define dgraphInduce2               _SCOTCHdgraphInduce2
#define dgraphInit                  _SCOTCHdgraphInit
#define dgraphLoad                  _SCOTCHdgraphLoad
#define dgraphMatchInit             _SCOTCHdgraphMatchInit
#define dgraphMatchExit             _SCOTCHdgraphMatchExit
#define dgraphMatchSync             _SCOTCHdgraphMatchSync
#define dgraphMatchSyncColl         _SCOTCHdgraphMatchSyncColl
#define dgraphMatchSyncPtop         _SCOTCHdgraphMatchSyncPtop
#define dgraphMatchCheck            _SCOTCHdgraphMatchCheck
#define dgraphMatchHl               _SCOTCHdgraphMatchHl
#define dgraphMatchHy               _SCOTCHdgraphMatchHy
#define dgraphMatchLc               _SCOTCHdgraphMatchLc
#define dgraphMatchLy               _SCOTCHdgraphMatchLy
#define dgraphMatchSc               _SCOTCHdgraphMatchSc
#define dgraphRedist                _SCOTCHdgraphRedist
#define dgraphSave                  _SCOTCHdgraphSave
#define dgraphScatter               _SCOTCHdgraphScatter
#define dgraphView                  _SCOTCHdgraphView

#define dmapInit                    _SCOTCHdmapInit
#define dmapExit                    _SCOTCHdmapExit
#define dmapAdd                     _SCOTCHdmapAdd
#define dmapTerm                    _SCOTCHdmapTerm
#define dmapSave                    _SCOTCHdmapSave

#define dorderDispose               _SCOTCHdorderDispose
#define dorderExit                  _SCOTCHdorderExit
#define dorderFree                  _SCOTCHdorderFree
#define dorderFrst                  _SCOTCHdorderFrst
#define dorderGather                _SCOTCHdorderGather
#define dorderGatherTree            _SCOTCHdorderGatherTree
#define dorderInit                  _SCOTCHdorderInit
#define dorderNew                   _SCOTCHdorderNew
#define dorderNewSequ               _SCOTCHdorderNewSequ
#define dorderNewSequIndex          _SCOTCHdorderNewSequIndex
#define dorderPerm                  _SCOTCHdorderPerm
#define dorderSave                  _SCOTCHdorderSave
#define dorderSaveBlock             _SCOTCHdorderSaveBlock
#define dorderSaveMap               _SCOTCHdorderSaveMap
#define dorderSaveTree              _SCOTCHdorderSaveTree
#define dorderSaveTree2             _SCOTCHdorderSaveTree2
#define dorderCblkDist              _SCOTCHdorderCblkDist
#define dorderTreeDist              _SCOTCHdorderTreeDist

#define fiboTreeCheck               _SCOTCHfiboTreeCheck
#define fiboTreeConsolidate         _SCOTCHfiboTreeConsolidate
/* #define fiboTreeAdd              _SCOTCHfiboTreeAdd Already a macro */
#define fiboTreeDel                 _SCOTCHfiboTreeDel
#define fiboTreeExit                _SCOTCHfiboTreeExit
#define fiboTreeFree                _SCOTCHfiboTreeFree
#define fiboTreeInit                _SCOTCHfiboTreeInit
#define fiboTreeMin                 _SCOTCHfiboTreeMin

#define gainTablAddLin              _SCOTCHgainTablAddLin
#define gainTablAddLog              _SCOTCHgainTablAddLog
#define gainTablCheck               _SCOTCHgainTablCheck
#ifdef SCOTCH_DEBUG_GAIN1                         /* If not already redefined as accelerated macro */
#define gainTablDel                 _SCOTCHgainTablDel
#endif /* SCOTCH_DEBUG_GAIN1 */
#define gainTablExit                _SCOTCHgainTablExit
#define gainTablFree                _SCOTCHgainTablFree
#define gainTablFrst                _SCOTCHgainTablFrst
#define gainTablInit                _SCOTCHgainTablInit
#define gainTablNext                _SCOTCHgainTablNext

#define geomExit                    _SCOTCHgeomExit
#define geomInit                    _SCOTCHgeomInit

#define graphInit                   _SCOTCHgraphInit
#define graphExit                   _SCOTCHgraphExit
#define graphFree                   _SCOTCHgraphFree
#define graphLoad                   _SCOTCHgraphLoad
#define graphLoad2                  _SCOTCHgraphLoad2
#define graphSave                   _SCOTCHgraphSave
#define graphBand                   _SCOTCHgraphBand
#define graphBase                   _SCOTCHgraphBase
#define graphCheck                  _SCOTCHgraphCheck
#define graphCoarsen                _SCOTCHgraphCoarsen
#define graphInduceList             _SCOTCHgraphInduceList
#define graphInducePart             _SCOTCHgraphInducePart
#define graphMatch                  _SCOTCHgraphMatch
#define graphMatchInit              _SCOTCHgraphMatchInit
#define graphGeomLoadChac           _SCOTCHgraphGeomLoadChac
#define graphGeomLoadHabo           _SCOTCHgraphGeomLoadHabo
#define graphGeomLoadMmkt           _SCOTCHgraphGeomLoadMmkt
#define graphGeomLoadScot           _SCOTCHgraphGeomLoadScot
#define graphGeomSaveChac           _SCOTCHgraphGeomSaveChac
#define graphGeomSaveScot           _SCOTCHgraphGeomSaveScot
#define graphGeomSaveMmkt           _SCOTCHgraphGeomSaveMmkt
#define graphPtscotch               _SCOTCHgraphPtscotch

#define hallOrderHdHalmd            _SCOTCHhallOrderHdHalmd
#define hallOrderHfR2hamdf4         _SCOTCHhallOrderHfR2hamdf4
#define hallOrderHxBuild            _SCOTCHhallOrderHxBuild
#define hallOrderHxTree             _SCOTCHhallOrderHxTree

#define hdgraphorderststratab       _SCOTCHhdgraphorderststratab
#define hdgraphInit                 _SCOTCHhdgraphInit
#define hdgraphExit                 _SCOTCHhdgraphExit
#define hdgraphCheck                _SCOTCHhdgraphCheck
#define hdgraphFold                 _SCOTCHhdgraphFold
#define hdgraphFold2                _SCOTCHhdgraphFold2
#define hdgraphGather               _SCOTCHhdgraphGather
#define hdgraphInduceList           _SCOTCHhdgraphInduceList
#define hdgraphOrderNd              _SCOTCHhdgraphOrderNd
#define hdgraphOrderSi              _SCOTCHhdgraphOrderSi
#define hdgraphOrderSq              _SCOTCHhdgraphOrderSq
#define hdgraphOrderSq2             _SCOTCHhdgraphOrderSq2
#define hdgraphOrderSt              _SCOTCHhdgraphOrderSt

#define hgraphorderststratab        _SCOTCHhgraphorderststratab
#define hgraphInit                  _SCOTCHhgraphInit
#define hgraphExit                  _SCOTCHhgraphExit
#define hgraphFree                  _SCOTCHhgraphFree
#define hgraphInduceList            _SCOTCHhgraphInduceList
#define hgraphCheck                 _SCOTCHhgraphCheck
#define hgraphOrderBl               _SCOTCHhgraphOrderBl
#define hgraphOrderCp               _SCOTCHhgraphOrderCp
#define hgraphOrderGp               _SCOTCHhgraphOrderGp
#define hgraphOrderHd               _SCOTCHhgraphOrderHd
#define hgraphOrderHf               _SCOTCHhgraphOrderHf
#define hgraphOrderHxFill           _SCOTCHhgraphOrderHxFill
#define hgraphOrderKp               _SCOTCHhgraphOrderKp
#define hgraphOrderNd               _SCOTCHhgraphOrderNd
#define hgraphOrderSi               _SCOTCHhgraphOrderSi
#define hgraphOrderSt               _SCOTCHhgraphOrderSt
#define hgraphUnhalo                _SCOTCHhgraphUnhalo

#define hmeshorderststratab         _SCOTCHhmeshorderststratab
#define hmeshExit                   _SCOTCHhmeshExit
#define hmeshBase                   _SCOTCHhmeshBase
#define hmeshCheck                  _SCOTCHhmeshCheck
#define hmeshInducePart             _SCOTCHhmeshInducePart
#define hmeshHgraph                 _SCOTCHhmeshHgraph
#define hmeshMesh                   _SCOTCHhmeshMesh
#define hmeshOrderBl                _SCOTCHhmeshOrderBl
#define hmeshOrderCp                _SCOTCHhmeshOrderCp
#define hmeshOrderGp                _SCOTCHhmeshOrderGp
#define hmeshOrderGr                _SCOTCHhmeshOrderGr
#define hmeshOrderHd                _SCOTCHhmeshOrderHd
#define hmeshOrderHf                _SCOTCHhmeshOrderHf
#define hmeshOrderHxFill            _SCOTCHhmeshOrderHxFill
#define hmeshOrderNd                _SCOTCHhmeshOrderNd
#define hmeshOrderSi                _SCOTCHhmeshOrderSi
#define hmeshOrderSt                _SCOTCHhmeshOrderSt

#define kdgraphmapststratab         _SCOTCHkdgraphmapststratab
#define kdgraphInit                 _SCOTCHkdgraphInit
#define kdgraphExit                 _SCOTCHkdgraphExit
#define kdgraphGather               _SCOTCHkdgraphGather
#define kdgraphMapRb                _SCOTCHkdgraphMapRb
#define kdgraphMapRbAdd2            _SCOTCHkdgraphMapRbAdd2
#define kdgraphMapRbAddBoth         _SCOTCHkdgraphMapRbAddBoth
#define kdgraphMapRbAddOne          _SCOTCHkdgraphMapRbAddOne
#define kdgraphMapRbAddPart         _SCOTCHkdgraphMapRbAddPart
#define kdgraphMapRbMap             _SCOTCHkdgraphMapRbMap
#define kdgraphMapRbPart            _SCOTCHkdgraphMapRbPart
#define kdgraphMapSt                _SCOTCHkdgraphMapSt

#define kgraphmapststratab          _SCOTCHkgraphmapststratab
#define kgraphInit                  _SCOTCHkgraphInit
#define kgraphExit                  _SCOTCHkgraphExit
#define kgraphCheck                 _SCOTCHkgraphCheck
#define kgraphBand                  _SCOTCHkgraphBand
#define kgraphCost                  _SCOTCHkgraphCost
#define kgraphFron                  _SCOTCHkgraphFron
#define kgraphFrst                  _SCOTCHkgraphFrst
#define kgraphMapBd                 _SCOTCHkgraphMapBd
#define kgraphMapCp                 _SCOTCHkgraphMapCp
#define kgraphMapDf                 _SCOTCHkgraphMapDf
#define kgraphMapEx                 _SCOTCHkgraphMapEx
#define kgraphMapFm                 _SCOTCHkgraphMapFm
#define kgraphMapMl                 _SCOTCHkgraphMapMl
#define kgraphMapRb                 _SCOTCHkgraphMapRb
#define kgraphMapRbMap              _SCOTCHkgraphMapRbMap
#define kgraphMapRbBgraph           _SCOTCHkgraphMapRbBgraph
#define kgraphMapRbPart             _SCOTCHkgraphMapRbPart
#define kgraphMapRbVfloBuild        _SCOTCHkgraphMapRbVfloBuild
#define kgraphMapRbVfloMerge        _SCOTCHkgraphMapRbVfloMerge
#define kgraphMapRbVfloSplit        _SCOTCHkgraphMapRbVfloSplit
#define kgraphMapSt                 _SCOTCHkgraphMapSt
#define kgraphStoreInit             _SCOTCHkgraphStoreInit
#define kgraphStoreExit             _SCOTCHkgraphStoreExit
#define kgraphStoreSave             _SCOTCHkgraphStoreSave
#define kgraphStoreUpdt             _SCOTCHkgraphStoreUpdt

#define listInit                    _SCOTCHlistInit
#define listExit                    _SCOTCHlistExit
#define listAlloc                   _SCOTCHlistAlloc
#define listFree                    _SCOTCHlistFree
#define listLoad                    _SCOTCHlistLoad
#define listSave                    _SCOTCHlistSave
#define listSort                    _SCOTCHlistSort
#define listCopy                    _SCOTCHlistCopy

#define mapInit                     _SCOTCHmapInit
#define mapInit2                    _SCOTCHmapInit2
#define mapExit                     _SCOTCHmapExit
#define mapAlloc                    _SCOTCHmapAlloc
#define mapBuild                    _SCOTCHmapBuild
#define mapCopy                     _SCOTCHmapCopy
#define mapFree                     _SCOTCHmapFree
#define mapFrst                     _SCOTCHmapFrst
#define mapLoad                     _SCOTCHmapLoad
#define mapMerge                    _SCOTCHmapMerge
#define mapResize                   _SCOTCHmapResize
#define mapResize2                  _SCOTCHmapResize2
#define mapSave                     _SCOTCHmapSave
#define mapTerm                     _SCOTCHmapTerm

#define meshInit                    _SCOTCHmeshInit
#define meshExit                    _SCOTCHmeshExit
#define meshFree                    _SCOTCHmeshFree
#define meshLoad                    _SCOTCHmeshLoad
#define meshSave                    _SCOTCHmeshSave
#define meshBase                    _SCOTCHmeshBase
#define meshGraph                   _SCOTCHmeshGraph
#define meshCoarsen                 _SCOTCHmeshCoarsen
#define meshInduceList              _SCOTCHmeshInduceList
#define meshInducePart              _SCOTCHmeshInducePart
#define meshInduceSepa              _SCOTCHmeshInduceSepa
#define meshCheck                   _SCOTCHmeshCheck
#define meshGeomLoadHabo            _SCOTCHmeshGeomLoadHabo
#define meshGeomLoadScot            _SCOTCHmeshGeomLoadScot
#define meshGeomSaveScot            _SCOTCHmeshGeomSaveScot

#define orderInit                   _SCOTCHorderInit
#define orderExit                   _SCOTCHorderExit
#define orderLoad                   _SCOTCHorderLoad
#define orderSave                   _SCOTCHorderSave
#define orderSaveMap                _SCOTCHorderSaveMap
#define orderSaveTree               _SCOTCHorderSaveTree
#define orderCheck                  _SCOTCHorderCheck
#define orderPeri                   _SCOTCHorderPeri
#define orderRang                   _SCOTCHorderRang
#define orderTree                   _SCOTCHorderTree

#define parsermethtokentab          _SCOTCHparsermethtokentab
#define parserparamcurr             _SCOTCHparserparamcurr
#define parserstratcurr             _SCOTCHparserstratcurr
#define parserstrattab              _SCOTCHparserstrattab

#define stratdummy                  _SCOTCHstratdummy
#define stratInit                   _SCOTCHstratInit
#define stratExit                   _SCOTCHstratExit
#define stratSave                   _SCOTCHstratSave
#define stratCondEval               _SCOTCHstratCondEval
#define stratCondExit               _SCOTCHstratCondExit
#define stratCondSave               _SCOTCHstratCondSave
#define stratParserInit             _SCOTCHstratParserInit
#define stratParserInput            _SCOTCHstratParserInput
#define stratParserLex              _SCOTCHstratParserLex
#define stratParserRemain           _SCOTCHstratParserRemain
#define stratParserSelect           _SCOTCHstratParserSelect
#define stratParserError            _SCOTCHstratParserError
#define stratParserParse            _SCOTCHstratParserParse
#define stratParserParse2           _SCOTCHstratParserParse2
#define stratTestEval               _SCOTCHstratTestEval
#define stratTestExit               _SCOTCHstratTestExit
#define stratTestSave               _SCOTCHstratTestSave

#define threadLaunch                _SCOTCHthreadLaunch
#define threadReduce                _SCOTCHthreadReduce
#define threadScan                  _SCOTCHthreadScan

#define vdgraphseparateststratab    _SCOTCHvdgraphseparateststratab
#define vdgraphCheck                _SCOTCHvdgraphCheck
#define vdgraphExit                 _SCOTCHvdgraphExit
#define vdgraphGatherAll            _SCOTCHvdgraphGatherAll
#define vdgraphInit                 _SCOTCHvdgraphInit
#define vdgraphSeparateBd           _SCOTCHvdgraphSeparateBd
#define vdgraphSeparateDf           _SCOTCHvdgraphSeparateDf
#define vdgraphSeparateMl           _SCOTCHvdgraphSeparateMl
#define vdgraphSeparateSq           _SCOTCHvdgraphSeparateSq
#define vdgraphSeparateSt           _SCOTCHvdgraphSeparateSt
#define vdgraphSeparateZr           _SCOTCHvdgraphSeparateZr
#define vdgraphStoreExit            _SCOTCHvdgraphStoreExit
#define vdgraphStoreInit            _SCOTCHvdgraphStoreInit
#define vdgraphStoreSave            _SCOTCHvdgraphStoreSave
#define vdgraphStoreUpdt            _SCOTCHvdgraphStoreUpdt
#define vdgraphZero                 _SCOTCHvdgraphZero

#define vgraphseparateststratab     _SCOTCHvgraphseparateststratab
#define vgraphInit                  _SCOTCHvgraphInit
#define vgraphExit                  _SCOTCHvgraphExit
#define vgraphCheck                 _SCOTCHvgraphCheck
#define vgraphZero                  _SCOTCHvgraphZero
#define vgraphSeparateBd            _SCOTCHvgraphSeparateBd
#define vgraphSeparateDf            _SCOTCHvgraphSeparateDf
#define vgraphSeparateEs            _SCOTCHvgraphSeparateEs
#define vgraphSeparateFm            _SCOTCHvgraphSeparateFm
#define vgraphSeparateGg            _SCOTCHvgraphSeparateGg
#define vgraphSeparateGp            _SCOTCHvgraphSeparateGp
#define vgraphSeparateMl            _SCOTCHvgraphSeparateMl
#define vgraphSeparateMt            _SCOTCHvgraphSeparateMt
#define vgraphSeparateSt            _SCOTCHvgraphSeparateSt
#define vgraphSeparateTh            _SCOTCHvgraphSeparateTh
#define vgraphSeparateVw            _SCOTCHvgraphSeparateVw
#define vgraphSeparateZr            _SCOTCHvgraphSeparateZr
#define vgraphStoreInit             _SCOTCHvgraphStoreInit
#define vgraphStoreExit             _SCOTCHvgraphStoreExit
#define vgraphStoreSave             _SCOTCHvgraphStoreSave
#define vgraphStoreUpdt             _SCOTCHvgraphStoreUpdt

#define vmeshseparateststratab      _SCOTCHvmeshseparateststratab
#define vmeshExit                   _SCOTCHvmeshExit
#define vmeshCheck                  _SCOTCHvmeshCheck
#define vmeshZero                   _SCOTCHvmeshZero
#define vmeshSeparateFm             _SCOTCHvmeshSeparateFm
#define vmeshSeparateGg             _SCOTCHvmeshSeparateGg
#define vmeshSeparateGr             _SCOTCHvmeshSeparateGr
#define vmeshSeparateMl             _SCOTCHvmeshSeparateMl
#define vmeshSeparateSt             _SCOTCHvmeshSeparateSt
#define vmeshSeparateZr             _SCOTCHvmeshSeparateZr
#define vmeshStoreInit              _SCOTCHvmeshStoreInit
#define vmeshStoreExit              _SCOTCHvmeshStoreExit
#define vmeshStoreSave              _SCOTCHvmeshStoreSave
#define vmeshStoreUpdt              _SCOTCHvmeshStoreUpdt

#define wgraphpartststratab         _SCOTCHwgraphpartststratab
#define wgraphAlloc                 _SCOTCHwgraphAlloc
#define wgraphInit                  _SCOTCHwgraphInit
#define wgraphExit                  _SCOTCHwgraphExit
#define wgraphCheck                 _SCOTCHwgraphCheck
#define wgraphZero                  _SCOTCHwgraphZero
#define wgraphPartFm                _SCOTCHwgraphPartFm
#define wgraphPartGg                _SCOTCHwgraphPartGg
#define wgraphPartGp                _SCOTCHwgraphPartGp
#define wgraphPartMl                _SCOTCHwgraphPartMl
#define wgraphPartRb                _SCOTCHwgraphPartRb
#define wgraphPartSt                _SCOTCHwgraphPartSt
#define wgraphPartZr                _SCOTCHwgraphPartZr
#define wgraphStoreInit             _SCOTCHwgraphStoreInit
#define wgraphStoreExit             _SCOTCHwgraphStoreExit
#define wgraphStoreSave             _SCOTCHwgraphStoreSave
#define wgraphStoreUpdt             _SCOTCHwgraphStoreUpdt
#endif /* SCOTCH_RENAME */
