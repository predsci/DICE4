      subroutine gensingle(sh,school,nparam,
     $     curpars,ndata, tps, rtn, ndays, epi_model) 

      implicit none

      integer iday_per_week
      parameter(iday_per_week=7)
      integer nstep
      parameter (nstep = 6)
      integer ndata, nparam, ndays
      integer epi_model
      real*8 sh(ndata), school(ndata)
      real*8 curpars(nparam)
      real*8 pC, e_bckgrnd
      real*8 tps(ndata)
      real*8 rtn(ndata)
      real*8 dsdt((ndays)*nstep)
      integer imid((ndata+1))


!
! Need to find the mid-point
! This will need to be further cleaned - the CDC data is absolute day numbers..
!

      if (tps(1) > (iday_per_week * 10)) Then
         call BuildIMIDforCDC(ndata, nstep, tps, imid)
      else
         call BuildIMID(ndata, nstep, tps, imid)
      endif

      pC = curpars(5)
      e_bckgrnd = curpars(8)

      rtn = 0.0d0
      dsdt = 0.0d0


      select case (epi_model) 
         case (1) 
            call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (2) 
            call RK4SEIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
        case (3) 
           call RK4VecSIROneD(ndata,ndays,nstep,tps,
     $          nparam,curpars,dsdt)
        case (4) 
           call RK4VecSEIROneD(ndata,ndays,nstep,tps,
     $          nparam,curpars,dsdt)              
        case default 
           call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $            nparam,curpars,dsdt)
      end select 


      call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd, 
     $     rtn)


      
      return
      end subroutine gensingle


