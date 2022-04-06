      subroutine getbmat(bmati, bmat)
      
      integer bmati
      character bmat

      select case (bmati)
        case (0)
            bmat = 'I'
        case (1)
            bmat = 'G'
        case default
            bmat = 'I'
      end select

      end



      subroutine getwhich(whichi, which)

      integer whichi
      character which*2

      select case (whichi)
        case (0)
            which = 'LM'
        case (1)
            which = 'SM'
        case (2)
            which = 'LR'
        case (3)
            which = 'SR'
        case (4)
            which = 'LI'
        case (5)
            which = 'SI'
        case (6)
            which = 'LA'
        case (7)
            which = 'SA'
        case (8)
            which = 'BE'
        case default
            which = 'LM'
      end select

      end



      subroutine gethowmny(howmnyi, howmny)

      integer howmnyi
      character howmny

      select case (howmnyi)
        case (0)
            howmny = 'A'
        case (1)
            howmny = 'P'
        case (2)
            howmny = 'S'
        case default
            howmny = 'A'
      end select

      end



      subroutine dsaupdwr
     &   ( ido, bmati, n, whichi, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )
      
      implicit none

      integer ido, n, nev, ncv, ldv, iparam(11), ipntr(11), lworkl, info
      integer bmati, whichi
      character bmat, which*2
      double precision tol, resid(n), v(ldv, ncv), workd(3*n),
     &                 workl(lworkl)
      
      external dsaupd, getbmat, getwhich
      
      call getbmat(bmati, bmat)
      call getwhich(whichi, which)

      call dsaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )

      end



      subroutine dseupdwr(rvec  , howmnyi, select, d    ,
     &                    z     , ldz    , sigma , bmati,
     &                    n     , whichi , nev   , tol  ,
     &                    resid , ncv    , v     , ldv  ,
     &                    iparam, ipntr  , workd , workl,
     &                    lworkl, info )
      
      implicit none

      integer ldz, n, nev, ncv, ldv, iparam(7), ipntr(11), lworkl, info
      double precision d(nev), z(ldz, nev), sigma, tol, resid(n)
      double precision v(ldv, ncv), workd(2*n), workl(lworkl)
      integer howmnyi, bmati, whichi
      character howmny, bmat, which*2
      logical rvec, select(ncv)

      external dseupd, gethowmny, getbmat, getwhich

      call gethowmny(howmnyi, howmny)
      call getbmat(bmati, bmat)
      call getwhich(whichi, which)

      call dseupd(rvec  , howmny, select, d    ,
     &            z     , ldz   , sigma , bmat ,
     &            n     , which , nev   , tol  ,
     &            resid , ncv   , v     , ldv  ,
     &            iparam, ipntr , workd , workl,
     &            lworkl, info )

      end



      subroutine dnaupdwr
     &   ( ido, bmati, n, whichi, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )

      implicit none

      integer ido, n, nev, ncv, ldv, iparam(11), ipntr(14), lworkl, info
      integer bmati, whichi
      character bmat, which*2
      double precision tol, resid(n), v(ldv, ncv), workd(3*n),
     &                 workl(lworkl)
      
      external dnaupd, getbmat, getwhich
      
      call getbmat(bmati, bmat)
      call getwhich(whichi, which)
      
      call dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
     &     ipntr, workd, workl, lworkl, info )

      end



      subroutine dneupdwr(rvec , howmnyi, select, dr    , di,    
     &                    z    , ldz    , sigmar, sigmai, workev,
     &                    bmati, n      , whichi, nev   , tol,
     &                    resid, ncv    , v     , ldv   , iparam,
     &                    ipntr, workd  , workl , lworkl, info)

      implicit none

      integer ldz, n, nev, ncv, ldv, iparam(11), ipntr(14), lworkl, info
      double precision dr(nev+1), di(nev+1), z(ldz,*), sigmar, sigmai
      double precision workev(3*ncv), tol, resid(n), v(ldv, ncv)
      double precision workd(3*n), workl(lworkl)
      integer howmnyi, bmati, whichi
      character howmny, bmat, which*2
      logical rvec, select(ncv)

      external dneupd, gethowmny, getbmat, getwhich

      call gethowmny(howmnyi, howmny)
      call getbmat(bmati, bmat)
      call getwhich(whichi, which)
      
      call dneupd (rvec , howmny, select, dr    , di,    
     &             z    , ldz   , sigmar, sigmai, workev,
     &             bmat , n     , which , nev   , tol,
     &             resid, ncv   , v     , ldv   , iparam,
     &             ipntr, workd , workl , lworkl, info)

      end
