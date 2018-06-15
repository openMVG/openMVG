/* Copyright 2009 ENSEIRB, INRIA & CNRS
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
/**                for the ESMUMPS library module.         **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                                                        **/
/************************************************************/

#define MODULE_H

/*
** Function renaming.
*/

#if ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME))
#define clockGet                    _SCOTCHclockGet

#define fileNameDistExpand          _SCOTCHfileNameDistExpand 

#define usagePrint                  _SCOTCHusagePrint

#define errorPrint                  SCOTCH_errorPrint
#define errorPrintW                 SCOTCH_errorPrintW
#define errorProg                   SCOTCH_errorProg

#define intLoad                     _SCOTCHintLoad
#define intSave                     _SCOTCHintSave
#define intAscn                     _SCOTCHintAscn
#define intPerm                     _SCOTCHintPerm
#define intRandReset                _SCOTCHintRandReset
#define intRandInit                 _SCOTCHintRandInit
/* #define intRandVal               _SCOTCHintRandVal Already a macro */
#define intSearchDicho              _SCOTCHintSearchDicho
#define intSort1asc1                _SCOTCHintSort1asc1
#define intSort2asc1                _SCOTCHintSort2asc1
#define intSort2asc2                _SCOTCHintSort2asc2
#define intSort3asc1                _SCOTCHintSort3asc1

#define memAllocGroup               _SCOTCHmemAllocGroup
#define memReallocGroup             _SCOTCHmemReallocGroup
#define memOffset                   _SCOTCHmemOffset
#endif /* ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME)) */
