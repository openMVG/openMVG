************************************************************
**                                                        **
**   NAME       : scotch_example.f                        **
**                                                        **
**   AUTHOR     : Francois PELLEGRINI                     **
**                                                        **
**   FUNCTION   : FORTRAN testbed for the LibSCOTCH       **
**                library routines.                       **
**                                                        **
**   DATES      : # Version 3.4  : from : 04 feb 2000     **
**                                 to     07 feb 2000     **
**                # Version 4.0  : from : 13 mar 2005     **
**                                 to     13 mar 2005     **
**                                                        **
*234567*****************************************************

        PROGRAM SCOTCH_TEST
        IMPLICIT NONE
        INCLUDE "scotchf.h"
          DOUBLEPRECISION   SCOTCHGRAPH (SCOTCH_GRAPHDIM)
          INTEGER           VERTNBR
          DATA              VERTNBR / 3 /
          INTEGER           EDGENBR
          DATA              EDGENBR / 4 /
          INTEGER           VERTTAB (4)
          DATA              VERTTAB / 1, 2, 4, 5 /
          INTEGER           EDGETAB (4)
          DATA              EDGETAB / 2, 1, 3, 2 /
          INTEGER           INDXTAB (1)
          INTEGER           IDXVERTNBR
          INTEGER           IDXVERTTABIDX, IDXVENDTABIDX
          INTEGER           IDXVELOTABIDX, IDXVLBLTABIDX
          INTEGER           IDXEDGENBR
          INTEGER           IDXEDGETABIDX, IDXEDLOTABIDX
          INTEGER           IDXBASEVAL, IDXFLAGVAL
          INTEGER           IERR

          PRINT *, 'Starting'

          CALL SCOTCHFGRAPHINIT (SCOTCHGRAPH (1), IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Cannot initialize graph'
            STOP
          ENDIF

          CALL SCOTCHFGRAPHBUILD (SCOTCHGRAPH (1), 1, VERTNBR,
     *                            VERTTAB (1), VERTTAB (2),
     *                            VERTTAB (1), VERTTAB (1),
     *                            EDGENBR,
     *                            EDGETAB (1), EDGETAB (1), IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Cannot build graph'
            STOP
          ENDIF

          CALL SCOTCHFGRAPHCHECK (SCOTCHGRAPH (1), IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Invalid check'
            STOP
          ENDIF

          PRINT *, 'Outputing original graph'

          CALL SCOTCHFGRAPHSAVE (SCOTCHGRAPH (1), 1, IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Invalid graph output'
            STOP
          ENDIF

          CALL SCOTCHFGRAPHDATA (SCOTCHGRAPH (1), INDXTAB (1),
     *                           IDXBASEVAL, IDXVERTNBR,
     *                           IDXVERTTABIDX, IDXVENDTABIDX,
     *                           IDXVELOTABIDX, IDXVLBLTABIDX,
     *                           IDXEDGENBR,
     *                           IDXEDGETABIDX, IDXEDLOTABIDX,
     *                           IDXFLAGVAL, IERR);
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Cannot get graph data'
            STOP
          ENDIF

          PRINT *, 'Number of vertices : ', IDXVERTNBR
          PRINT *, 'Index of verttab   : ', IDXVERTTABIDX
          PRINT *, 'Index of vendtab   : ', IDXVENDTABIDX
          PRINT *, 'Index of velotab   : ', IDXVELOTABIDX
          PRINT *, 'Index of vlbltab   : ', IDXVLBLTABIDX
          PRINT *, 'Number of edges    : ', IDXEDGENBR
          PRINT *, 'Index of edgetab   : ', IDXEDGETABIDX
          PRINT *, 'Index of edlotab   : ', IDXEDLOTABIDX

          PRINT *, 'Updating vertex and edge arrays'
          INDXTAB (IDXVERTTABIDX + 1) = 3
          INDXTAB (IDXEDGETABIDX)     = 2
          INDXTAB (IDXEDGETABIDX + 1) = 3
          INDXTAB (IDXEDGETABIDX + 2) = 1
          INDXTAB (IDXEDGETABIDX + 3) = 1

          PRINT *, 'Outputting updated graph'

          CALL SCOTCHFGRAPHCHECK (SCOTCHGRAPH (1), IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Invalid check'
            STOP
          ENDIF

          CALL SCOTCHFGRAPHSAVE (SCOTCHGRAPH (1), 1, IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Invalid graph output'
            STOP
          ENDIF

          CALL SCOTCHFGRAPHEXIT (SCOTCHGRAPH (1), IERR)
          IF (IERR .NE. 0) THEN
            PRINT *, 'ERROR : MAIN : Cannot destroy graph'
            STOP
          ENDIF

          PRINT *, 'Test complete'

          RETURN
        END
