!* Copyright 2004,2007,2009,2010,2012 IPB, Universite de Bordeaux, INRIA & CNRS
!*
!* This file is part of the Scotch software package for static mapping,
!* graph partitioning and sparse matrix ordering.
!*
!* This software is governed by the CeCILL-C license under French law
!* and abiding by the rules of distribution of free software. You can
!* use, modify and/or redistribute the software under the terms of the
!* CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
!* URL: "http://www.cecill.info".
!* 
!* As a counterpart to the access to the source code and rights to copy,
!* modify and redistribute granted by the license, users are provided
!* only with a limited warranty and the software's author, the holder of
!* the economic rights, and the successive licensors have only limited
!* liability.
!* 
!* In this respect, the user's attention is drawn to the risks associated
!* with loading, using, modifying and/or developing or reproducing the
!* software by the user in light of its specific status of free software,
!* that may mean that it is complicated to manipulate, and that also
!* therefore means that it is reserved for developers and experienced
!* professionals having in-depth computer knowledge. Users are therefore
!* encouraged to load and test the software's suitability as regards
!* their requirements in conditions enabling the security of their
!* systems and/or data to be ensured and, more generally, to use and
!* operate it in the same conditions as regards security.
!* 
!* The fact that you are presently reading this means that you have had
!* knowledge of the CeCILL-C license and that you accept its terms.
!*
!***********************************************************
!*                                                        **
!*   NAME       : library_pt_f.h                          **
!*                                                        **
!*   AUTHOR     : Francois PELLEGRINI                     **
!*                                                        **
!*   FUNCTION   : FORTRAN declaration file for the        **
!*                LibPtscotch parallel static mapping and **
!*                sparse matrix block ordering sequential **
!*                library.                                **
!*                                                        **
!*   DATES      : # Version 3.4  : from : 04 feb 2000     **
!*                                 to     22 oct 2001     **
!*                # Version 4.0  : from : 16 jan 2004     **
!*                                 to     16 jan 2004     **
!*                # Version 5.0  : from : 26 apr 2006     **
!*                                 to     26 apr 2006     **
!*                # Version 5.1  : from : 26 mar 2009     **
!*                                 to     12 feb 2011     **
!*                # Version 6.0  : from : 22 oct 2011     **
!*                                 to     27 nov 2012     **
!*                                                        **
!***********************************************************

!* Flag definitions for the coarsening
!* routines.

        INTEGER SCOTCH_COARSENNOMERGE
        PARAMETER (SCOTCH_COARSENNOMERGE = 16384)

!* Flag definitions for the strategy
!* string selection routines.

        INTEGER SCOTCH_STRATDEFAULT
        INTEGER SCOTCH_STRATQUALITY
        INTEGER SCOTCH_STRATSPEED
        INTEGER SCOTCH_STRATBALANCE
        INTEGER SCOTCH_STRATSAFETY
        INTEGER SCOTCH_STRATSCALABILITY
        INTEGER SCOTCH_STRATRECURSIVE
        INTEGER SCOTCH_STRATREMAP
        INTEGER SCOTCH_STRATLEVELMAX
        INTEGER SCOTCH_STRATLEVELMIN
        INTEGER SCOTCH_STRATLEAFSIMPLE
        INTEGER SCOTCH_STRATSEPASIMPLE
        PARAMETER (SCOTCH_STRATDEFAULT     = 0)
        PARAMETER (SCOTCH_STRATQUALITY     = 1)
        PARAMETER (SCOTCH_STRATSPEED       = 2)
        PARAMETER (SCOTCH_STRATBALANCE     = 4)
        PARAMETER (SCOTCH_STRATSAFETY      = 8)
        PARAMETER (SCOTCH_STRATSCALABILITY = 16)
        PARAMETER (SCOTCH_STRATRECURSIVE   = 256)
        PARAMETER (SCOTCH_STRATREMAP       = 512)
        PARAMETER (SCOTCH_STRATLEVELMAX    = 4096)
        PARAMETER (SCOTCH_STRATLEVELMIN    = 8192)
        PARAMETER (SCOTCH_STRATLEAFSIMPLE  = 16384)
        PARAMETER (SCOTCH_STRATSEPASIMPLE  = 32768)

!* Size definitions for the SCOTCH opaque
!* structures. These structures must be
!* allocated as arrays of DOUBLEPRECISION
!* values for proper padding. The dummy
!* sizes are computed at compile-time by
!* program "dummysizes".

        INTEGER SCOTCH_ARCHDIM
        INTEGER SCOTCH_DGRAPHDIM
        INTEGER SCOTCH_DGRAPHHALOREQDIM
        INTEGER SCOTCH_DORDERDIM
        INTEGER SCOTCH_GEOMDIM
        INTEGER SCOTCH_GRAPHDIM
        INTEGER SCOTCH_MAPDIM
        INTEGER SCOTCH_MESHDIM
        INTEGER SCOTCH_ORDERDIM
        INTEGER SCOTCH_STRATDIM
        PARAMETER (SCOTCH_ARCHDIM          = DUMMYSIZEARCH)
        PARAMETER (SCOTCH_DGRAPHDIM        = DUMMYSIZEDGRAPH)
        PARAMETER (SCOTCH_DGRAPHHALOREQDIM = DUMMYSIZEDGRAPHHALOREQ)
        PARAMETER (SCOTCH_DORDERDIM        = DUMMYSIZEDORDER)
        PARAMETER (SCOTCH_GEOMDIM          = DUMMYSIZEGEOM)
        PARAMETER (SCOTCH_GRAPHDIM         = DUMMYSIZEGRAPH)
        PARAMETER (SCOTCH_MAPDIM           = DUMMYSIZEMAP)
        PARAMETER (SCOTCH_MESHDIM          = DUMMYSIZEMESH)
        PARAMETER (SCOTCH_ORDERDIM         = DUMMYSIZEORDER)
        PARAMETER (SCOTCH_STRATDIM         = DUMMYSIZESTRAT)
