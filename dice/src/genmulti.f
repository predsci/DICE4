      subroutine genmulti(n,sh,school,nparam,curpars,
     $     parCPL,fRij,ndata,tps,rtn,ndays,epi_model) 


      implicit none

      integer iday_per_week
      parameter(iday_per_week=7)
      integer nCPL
      parameter(nCPL=2)
      integer nstep
      parameter (nstep = 6)
      integer n, ndata, ndays, nparam
      integer epi_model
      real*8 sh(ndata,n), school(ndata,n)
      real*8 curpars(nparam,n)
      real*8 rtn(ndata,n)
      real*8 fN(n)
      real*8 pC(n)
      real*8 tps(ndata)
      real*8 e_bckgrnd(n)
      real*8 dsdt((ndays)*nstep,n)
! Unique to the coupling between regions
      real*8 fRij(n,n)
      real*8 parCPL(nCPL)
      real*8 fMij(n,n)
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


! integrate the ODEs to get the profiles

      rtn = 0.0d0
      dsdt = 0.0d0
    
      fN = curpars(1, 1:n)
      pC = curpars(5, 1:n)
      e_bckgrnd = curpars(8, 1:n)

      fMij = 0.0d0
 
      call buildMij(n,nCPL,fN,fRij,parCPL,fMij)
      

      select case (epi_model) 
         case (1) 
             call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
         case (2) 
             call RK4SEIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
          case default 
             call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
      end select 
    
      call weekly(n,ndata, ndays, nstep, imid, dsdt, pC,
     $     e_bckgrnd, rtn)

      return
      end subroutine genmulti
