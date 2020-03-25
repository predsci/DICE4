!----------------------------------------------------------------

      subroutine fnProposeParamUpdates(nparam,curval,
     $     valmin,valmax,step,logflag,
     $     parupdt,iseed,nopt,iopt)

      implicit none
      integer nparam, nopt, iopt(nopt)
      real*8 curval(nparam),valmin(nparam),valmax(nparam)
      real*8 parupdt(nparam),step(nparam)
      real*8 x, rv, rtn
      real*8 ran1
      real*8 SR_to_unit, SR_from_unit
      integer iseed
      integer i,j
      
      logical logflag(nparam)
      external SR_to_unit, SR_from_unit, ran1

       do j = 1, nopt
          i = iopt(j)
          parupdt(i)=curval(i)

          rv = ran1(iseed)

          rv = (rv - 0.50d0)*step(i)

! convert to a zero - one scale

          x = SR_to_unit(curval(i),valmin(i),valmax(i),
     $         logflag(i))

          x = x + rv

      if (x .lt. 0.0d0) x = 1.0d0 + x
      if (x .gt. 1.0d0) x = x - 1.0d0

! Do not use period boundary conditions here but rather re-smaple
c$$$      if (x .le. 0.0d0 .or. x .ge. 1.0d0) go to 101


! bring value back to original scale
      
         rtn = SR_from_unit(x,valmin(i),valmax(i),
     $        logflag(i))

         parupdt(i) = rtn

      enddo

      return
      end subroutine fnProposeParamUpdates


c----------------------------------------------------------------

      function SR_to_unit(y,ymin,ymax,logflag)

      implicit none
      real*8 y, ymin,ymax,rtn,SR_to_Unit
      logical logflag

      if (logflag) Then
         rtn = (log10(y) - log10(ymin)) /
     $        (log10(ymax)-log10(ymin))
      else
         rtn = (y - ymin)/(ymax - ymin)
      endif

      SR_to_unit = rtn
 
      return
      end function SR_to_unit

c----------------------------------------------------------------

      function SR_from_unit(x,ymin,ymax,logflag)
      
      implicit none
   
      real*8 x,ymin,ymax,rtn,SR_from_unit
      logical logflag


      if (logflag) Then
         rtn = ymin * 
     $        10.0**(x*(log10(ymax)-log10(ymin)))

      else
         rtn = ymin + (ymax-ymin)*x
      endif

      SR_from_unit = rtn

      return
      end function SR_from_unit

c----------------------------------------------------------------

      subroutine BuildIMIDforCDC(ndata,nstep,tps,imid)

      implicit none

      integer ndata,nstep
      integer imid((ndata+1))
      real*8 tps(ndata)
      real*8 t0
      integer i
      
! Build the vector needed for calculating the incidence
! days is needed only to calculate imid
! tps is the 'cumulative' day number
! This works correctly provided tps holds cumulative days numbers
! that start with the first month or week and not some arbitrary start

! this is for the CDC where we begin with EW27 - days start ~182
! If we move to start with first value being 7 days we can move 
! to using the cleaner dengue version..

      t0 =  (tps(2) - tps(1) )
      imid(1) = int(t0) * int(nstep * 0.5)

      do i = 1, ndata
         imid((i+1)) = imid(1) + int(tps(i) - tps(1) + t0)* nstep
        ! imid((i+1)) = imid(1) + int(tps(i) * nstep)
      enddo

      return
      end subroutine BuildIMIDforCDC

C----------------------------------------------------------------------

!--------------------------------------------------------------------------------
           
        subroutine weekly1D(ndata,ndays,nstep,imid,dsdt,pC,e_nonflu,x)
        
        implicit none
        integer ndata,nstep, ndays
        real*8 dsdt((ndays)*nstep)
        real*8 x(ndata), pC, e_nonflu
        integer i
        integer istart, iend
        integer imid((ndata+1))

        do i=1,ndata
c$$$           x(i) = dsdt(i*nstep*iday_per_week) - 
c$$$     $     dsdt(1+(i-1)*nstep*iday_per_week)

           iend   = imid(i+1)
           istart = imid(i)
           x(i) = dsdt(iend) - dsdt(istart)

        enddo

        x= x* pC +e_nonflu

        return
        end subroutine weekly1D

C--------------------------------------------------------------------------------
        function calcFit1D(y,gamay,x,wght,ndata)

        implicit none

        integer ndata, i
        real*8 y(ndata),x(ndata), gamay(ndata)
        real*8 wght(ndata)
        real*8 xi, yi,sum,val,CalcFit1D


c x is the simulated data
c y is the base profile

C calculate the P(yi,xi)


        sum = 0.0
        do i=1,ndata
           
           yi = y(i)
           xi = x(i)
           
           val = yi * log(xi) - xi - gamay(i)
           sum = sum + val  * wght(i) 
           
        enddo
        sum = -sum   
        
        CalcFit1D = sum 
           
        return
        end function CalcFit1D

C--------------------------------------------------------------------------------
        function Ratios(nparam,iopt,parnew,parold,ymu,sigma,Temp,
     $     shaveg,imodel)

        implicit none

        integer nparam,iopt(nparam), imodel
        real*8 parnew(nparam),parold(nparam)
        real*8 ymu(nparam),sigma(nparam), Temp, shaveg
        real*8 Ratios
        real*8 dx2New, dx2Old, del
        real*8 tmpOLD, tmpNEW
        integer i, j, i1, i2
        real*8 denom

        
        Ratios = 1.0d0
        del = 0.0d0

! The first two spots are pC and R0 which are true for all models

! c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd", "deltaR", "aparam", "alpha", "delta", "ts", "dur")

        do j = 1, 2
           i = iopt(j)
           dx2New = (parnew(i) - ymu(j))*(parnew(i) - ymu(j))
           dx2Old = (parold(i) - ymu(j))*(parold(i) - ymu(j))
           denom = (sigma(j) * sigma(j) * 2.0d0 * Temp)
           del = del  -(dx2New - dx2Old )/ denom

        enddo

! now we have to become model specific 

        if (imodel .eq. 1 ) Then
           i1 = iopt(3)
           i2 = iopt(4)
           j = 3
           tmpNEW = parnew(i1) * exp(-parnew(i2) * shaveg)
           tmpOLD = parold(i1) * exp(-parold(i2) * shaveg)
           dx2New = (tmpNEW - ymu(j))*(tmpNEW - ymu(j))
           dx2Old = (tmpOLD - ymu(j))*(tmpOLD - ymu(j))
           denom = (sigma(j) * sigma(j) * 2.0d0 * Temp)
           del = del -(dx2New - dx2Old )/ denom

        end if

! Simple SV term only 

        if (imodel .eq. 2) Then 
           i = iopt(3)
           j = 3
           dx2New = (parnew(i) - ymu(j))*(parnew(i) - ymu(j))
           dx2Old = (parold(i) - ymu(j))*(parold(i) - ymu(j))
           denom = (sigma(j) * sigma(j) * 2.0d0 * Temp)
           del = del -(dx2New - dx2Old )/ denom

        endif

! First handle the SH term and then the SV term
        if (imodel .eq. 3) Then
        
           i1 = iopt(3)
           i2 = iopt(4)
           j = 3
           tmpNEW = parnew(i1) * exp(-parnew(i2) * shaveg)
           tmpOLD = parold(i1) * exp(-parold(i2) * shaveg)
           dx2New = (tmpNEW - ymu(j))*(tmpNEW - ymu(j))
           dx2Old = (tmpOLD - ymu(j))*(tmpOLD - ymu(j))
           denom = (sigma(j) * sigma(j) * 2.0d0 * Temp)
           del = del -(dx2New - dx2Old )/ denom
           i = iopt(5)
           j = 4
           dx2New = (parnew(i) - ymu(j))*(parnew(i) - ymu(j))
           dx2Old = (parold(i) - ymu(j))*(parold(i) - ymu(j))
           denom = (sigma(j) * sigma(j) * 2.0d0 * Temp)
           del = del -(dx2New - dx2Old )/ denom

        endif

        Ratios = exp(del)

        return
        end function Ratios

C--------------------------------------------------------------------------------
