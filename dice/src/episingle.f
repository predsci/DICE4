      subroutine episingle(epi, gamaEpi,sh,school,wght,nparam,
     $     par,parmin,parmax,step,ilog,ymu,sigma,Temp,
     $     imask,iseed,nsamps,ithin,ndata,tps,rtn,accept_rate,
     $     curMin, tab, imodel, ndays, profiles,
     $     nRnd, epi_model) 


      implicit none

      integer iday_per_week, nblock
      parameter(iday_per_week=7,nblock=3)
      integer nstep
      parameter (nstep = 6)
! Y and gamay hold the synthetic profiles
      integer iseed, ndata,  nsamps, nparam
      integer ithin, imodel, ndays, nRnd
      integer epi_model
      integer ilog(nparam),imask(nparam)
      real*8 epi(ndata), gamaEpi(ndata)
      real*8 sh(ndata), school(ndata)
      real*8 wght(ndata)
      real*8 par(nparam)
      real*8 parmin(nparam), parmax(nparam),step(nparam)
      real*8 ymu(nparam),sigma(nparam),Temp
      real*8 shaveg
      real*8 tps(ndata)
      real*8 pC,e_bckgrnd
      real*8 rtn(ndata)
      real  tab(nsamps/ithin,nparam+1)
      real  profiles(nRnd, ndata)
      real*8 rtnBest(ndata),rtnNew(ndata)
c$$$      real*8 dsdt(ndata*nstep*iday_per_week)
      real*8 dsdt((ndays)*nstep)
      real*8 curpars(nparam),copypar(nparam),savepar(nparam)
      real*8 savestep(nparam),parupdt(nparam),parBest(nparam)
      real*8 curLLK, curMin
      real*8 fnewLLK
      real*8 ran1
      real*8 accept_rate(nblock)
      logical plog(nparam)
      integer i, k, icount
      integer noptv, ioptv(nparam)
      integer nopt(nblock),iopt(nparam,nblock)
      integer iaccept(nblock)
      integer iadapt(nblock)
      integer imid((ndata+1))
      real*8 AUX, Ratios
      real*8  range_min, range_max, scale, myaccept
      real*8 step_max, step_min
      real*8 calcFit1D
      real*8 scaleNdata
      integer ionep
      integer nlines, iline, iburn
      integer iweeks_per_year
      external calcFit1D,ran1


! Number of lines in the tab array

      nlines = nsamps/ithin

! calculate the yearly average value of SH
      shaveg = 0.0d0
      shaveg = sum(sh) / dble(ndata)
      
! For an adaptive size MCMC - decide if we need to update step size every 1% 
! or 1,000 steps
      scale = 2.0d0
      ionep = int(nsamps * 0.01)
      ionep = min(1000, ionep)
      range_min = 0.20d0
      range_max = 0.30d0
      step_max = 1.d0
      step_min = 1e-5
      iadapt = 0

      iweeks_per_year = 52 

!! Special case of daily data..
c$$$      if (sum(wght) .gt. iweeks_per_year ) Then
c$$$         scaleNdata = Temp
c$$$      else
         scaleNdata =sum(wght) * Temp
c$$$      endif
!!         scaleNData = 1.0d0

! ran1 works with -iseed

      iseed = -abs(iseed)

! ilog can have values of 
! 0 - uniform 
! 1 - log uniform
! 2 - Gaussian 

! set an auxilary array with the indices of the parameters we are optimizing 

      noptv = 0
      ioptv = 0
      nopt  = 0
      iopt = 0

      do i=1,nparam
         if (imask(i) > 0) Then
            noptv = noptv + 1
            ioptv(noptv) = i
         endif
      enddo
 ! It just helps to separate pC and R0 from everyone else.

      do i = 1,nparam
         if (imask(i) > 0) Then
            k = ilog(i) + 1
            nopt(k) = nopt(k) + 1
            iopt(nopt(k),k) = i
            if (par(i) .le. parmin(i) .or. par(i) .ge. parmax(i)) Then
               par(i) = parmin(i) + 0.50d0 * (parmax(i) - parmin(i))
            endif

         endif
      enddo



! convert 0 and 1 to true and false

      do i=1,nparam
         if (ilog(i) .eq. 0) Then
            plog(i) = .false.
         else
            plog(i) = .true.
         endif
      enddo

!
! Need to find the mid-point
! This will need to be further cleaned - the CDC data is absolute day numbers..
!

      if (tps(1) > (iday_per_week * 10)) Then
         call BuildIMIDforCDC(ndata, nstep, tps, imid)
      else
         call BuildIMID(ndata, nstep, tps, imid)
      endif

      curpars = par

! integrate the ODEs to get the profiles
! This is the parameter order
! c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd", "deltaR", "aparam", "alpha", "delta", "ts", "dur")
! The vector parameters follow
! c("bite","vec_k","muV","T_VH", "T_HV", "sigmaV")


      pC = curpars(5)
      e_bckgrnd = curpars(8)


      dsdt = 0.0d0
      rtn  = 0.0d0

!      
! Initial solution 
!
      select case (epi_model) 
         case (1) 
            call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (2) 
            call RK4SEIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (3) 
            call RK4VecSIROneD(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (4) 
            call RK4VecSEIROneD(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)            
         case default 
           call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $            nparam,curpars,dsdt)
      end select 


      call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd, 
     $     rtn)

!
! calculate Likelihood of the solution
!

      curLLK = calcFit1D(epi,gamaEpi,rtn,wght,ndata)

      curLLK = curLLK / (scaleNdata)

!
! MCMC loop starts here 
!
      curMIN = curLLK
      rtnBest = rtn

      iaccept = 0
      icount = 0

      savestep = step
     
      parupdt = curpars

      do i =1,nsamps
 
         do k = 1, nblock

            if (nopt(k) .eq. 0) go to 102
            savepar = curpars
            copypar = curpars
    
            
! Propose a new value for each parameter we are optimizing 
! half the steps will be small and half will be the given size
            call fnProposeParamUpdates(nparam,copypar,
     $           parmin,parmax,step,
     $           plog,parupdt,iseed,nopt(k),iopt(1:nopt(k),k))
            
!     update all the parameters and calculate a new solution and LLK

            curpars(iopt(1:nopt(k),k))= parupdt(iopt(1:nopt(k),k))

            pC = curpars(5)
            e_bckgrnd = curpars(8)

            rtnNew = 0.0d0
            dsdt   = 0.0d0

            select case (epi_model) 
            case (1) 
               call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $              nparam,curpars,dsdt)
            case (2) 
               call RK4SEIROneD(sh,school,ndata,ndays,nstep,tps,
     $              nparam,curpars,dsdt)
            case (3) 
               call RK4VecSIROneD(ndata,ndays,nstep,tps,
     $              nparam,curpars,dsdt)
            case (4) 
               call RK4VecSEIROneD(ndata,ndays,nstep,tps,
     $              nparam,curpars,dsdt)               
            case default 
               call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
            end select 

            rtnNew = 0.0d0
            call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd,
     $           rtnNew)

            fnewLLK = calcFit1D(epi,gamaEpi,rtnNew,wght,ndata)
            fnewLLK = fnewLLK / (scaleNdata)

            
            AUX = 1.0d0
            if (k .eq. 3) Then
               AUX = Ratios(nparam,iopt(1:nopt(k),k),
     $              curpars,savepar,ymu,sigma,Temp,shaveg,
     $              imodel)

            endif
            
            call MH1D(fnewLLK,curLLK,AUX,curMin,nparam,ndata,
     $           iseed,iaccept(k),iadapt(k),curpars,savepar,
     $           parBest,rtnNew,rtnBest)

 102        continue
            
         enddo  ! End of loop over blocks
 
         if (mod(i,ithin) .eq. 0) Then
            icount = icount + 1
            tab(icount,1:nparam)   = real(curpars)
            tab(icount,(nparam+1)) = real(curLLK)
         endif


         if (mod(i,ionep) .eq. 0) Then
             do k = 1,nblock
               if (nopt(k) .eq. 0) go to 202
               myaccept = dble(iadapt(k))/dble(ionep)
               if (myaccept > range_max .and.
     $         all(step(iopt(1:nopt(k),k))*scale < step_max)) Then
                 step(iopt(1:nopt(k),k))= step(iopt(1:nopt(k),k))*scale
               endif
               if (myaccept < range_min .and.
     $         all(step(iopt(1:nopt(k),k))/scale > step_min)) Then
                 step(iopt(1:nopt(k),k))= step(iopt(1:nopt(k),k))/scale
               endif
               iadapt(k) = 0    !reset acceptance number
 202           continue
               
            enddo
! Print information to the screen:
            call intpr('MCMC Step Number:',-1,i,1)
            call dblepr('curLLK',-1,curLLK,1)
            call dblepr('accept%',-1,dble(iaccept)/dble(i)*100.0,nblock)

         endif
         
      enddo

!update change information only every ithin iterations 

      accept_rate = dble(iaccept)/dble(nsamps) *100.0d0
               
! generate the profiles _ could really be done on the fly in the MCMC procedure
! see aove 

      iburn = nlines/5
      do i = 1, nRnd
         iline = iburn + floor(ran1(iseed)*(nlines - iburn)) 
         curpars = tab(iline, 1:nparam)
         pC = curpars(5)
         e_bckgrnd = curpars(8)
         select case (epi_model) 
         case (1) 
            call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (2) 
           call RK4SEIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (3) 
            call RK4VecSIROneD(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         case (4) 
            call RK4VecSEIROneD(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)            
         case default 
           call RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)
         end select 
         
         call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd, 
     $        rtn)
         profiles(i, 1:ndata) = real(rtn)
      enddo

!! End of generating profiles


      rtn = rtnBest
      par = parBest
      
      
      return
         
      end subroutine episingle

!
! -------------------------------------------------------------------
!
      subroutine MH1D(fnewLLK,curLLK,AUX,curMin,nparam,ndata,
     $              iseed,iaccept,iadapt,curpars,savepar,parBest,
     $              rtn,rtnBest)


      implicit none

      integer nparam,ndata,iseed,iaccept,iadapt
      real*8 fnewLLK, curLLK, AUX,curMin
      real*8 curpars(nparam),savepar(nparam)
      real*8 parBest(nparam)
      real*8 rtn(ndata), rtnBest(ndata)
      real*8 rnd,ran1,diff_LLK
      external ran1
      logical accept

      diff_LLK = fnewLLK - curLLK
            
      accept = .false.
 
      if (diff_LLK .le. 0.0d0 .and. AUX .ge. 1.0d0) Then
         accept = .true.
      else
         rnd = ran1(iseed)
         if ((exp(-diff_LLK)*AUX) .gt.rnd) Then
            accept = .true.
         else
            accept = .false.
         endif                 
      endif
               
! if step is accepted 

      if (accept) Then
         iaccept = iaccept + 1
         iadapt = iadapt + 1 
         savepar = curpars
         curLLK = fnewLLK
 
         if (curLLK .le. curMin) then !keep track of the best profile we have
            parBest = curpars
            curMin = curLLK
            rtnBest = rtn
         endif
!     if step is rejected - restore saved values of parameters
      else
         
         curpars = savepar
              
      endif

      return
      end subroutine MH1D

! ------------------------------------------------------------------------

      subroutine RK4SIROneD(sh,school,ndata,ndays,nstep,tps,
     $     nparam,curpars,dsdt)


      implicit none

      integer iday_per_week, np, nc
      parameter (iday_per_week = 7, np = 3, nc = 4)
      integer ndata, nstep, ndays, nparam
      real*8 sh(ndata), school(ndata)
      real*8 curpars(nparam)
      real*8 deltaR, Aparam, alpha
      real*8 delta, ts, dur
      real*8 dt
      real*8 t0, R0, Tg, sigma
      real*8 fN
c$$$      real*8 dsdt(ndata*nstep*iday_per_week)
      real*8 dsdt((ndays)*nstep)
      real*8 seed !initial number of infectious
      real*8 tps(ndata) !time-series
      real*8 tps2(0:(ndata+1))
      real*8 y(nc), tmpY(nc)
      real*8 dy1(nc), dy2(nc), dy3(nc), dy4(nc)
      real*8 t_cur, Rt
      real*8 p5, p3, p6
      real*8 pars(np)

      integer istep, iday, iweek
      integer icount, istart, iend

      dt = 1.0d0 /dble(nstep)

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0


      fN = curpars(1)
      Tg  = curpars(2)
      R0 = curpars(3)
      sigma = curpars(4)
      t0 = curpars(6)
      seed = curpars(7)
      deltaR = curpars(9)
      Aparam = curpars(10)
      alpha  = curpars(11)
      delta = curpars(12)
      ts = curpars(13) 
      dur = curpars(14)
 

!! Pack what we need - N, beta(t) and sigma - for the ODEs
      pars(1) = fN
      pars(2) = Tg
      pars(3) = R0/Tg ! will be updated below to beta(t)


! initialize things - this is the Y vector

      y(1) = fN - seed !fN * fraction-1.0
      y(2) = seed
      y(3) = 0.0d0
      y(4) = 0.0d0

      dsdt = 0.0d0

      tps2 = 0.0d0
      tps2(1:ndata) = tps
      tps2(0) = tps2(1) - (tps(2) - tps(1)) 
      tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 

c$$$      t_cur = tps(1)
      t_cur = tps2(0)

      icount = 0

!     Loop over ndays and in each day loop over nsteps
!     need to asign each day to a week because the sh and school data are weekly
!     tps has an arbitrary start (could be week 27 for example) but sh/school 
!     go from week 1 to week ndata hence we substract the start

       istart = int(tps2(0))
       iend   = int(tps2((ndata+1)))

       do iday = istart, iend
          iweek = (iday-int(tps2(0)))/iday_per_week + 1
          
          iweek = min(iweek,ndata)
          iweek = max(1,iweek)

          do istep=1,nstep
             t_cur = t_cur+dt

             if (t_cur <= t0) then
                go to 101
             end if
             
! calculate the time dependent R(t)
             Rt = R0 * 
     $            * (1.0d0+deltaR*exp(-Aparam*sh(iweek)))
     $            * (1.0d0-alpha*school(iweek))
             if (t_cur >= ts .and. t_cur <= (ts+dur)) Then
                Rt = Rt * (1.0d0 + delta)
             endif
             
             Rt = Rt / Tg


             pars(3) = Rt

             call derivSIR(pars, Y, dY1, np, nc)

             tmpY = Y + dt * dY1 * p5

             call derivSIR(pars, tmpY, dY2, np, nc)

             tmpY = Y + dt * dY2 * p5

             call derivSIR(pars, tmpY, dY3, np, nc)
               
             tmpY = Y + dt * dY3

             call derivSIR(pars, tmpY, dY4, np, nc)
               
             tmpY = Y + dt * ( dY1 * p6 + dY2 * p3 +
     $            dY3 * p3 + dY4 * p6 )
             Y = tmpY 
 101         continue
!     Update -dS/dt       
             icount= icount + 1
               
             dsdt(icount) = y(4)

          enddo                 ! End of loop over days
      enddo
       return
       end subroutine RK4SIROneD


!------------------------------------
C The Main Differential Equations- SIR model

      SUBROUTINE derivSIR(pars, Pop, dPop, np, nc)

      implicit none

      integer np, nc
      REAL*8 pars(np)
      REAL*8 Pop(nc), dPop(nc)

C
C     The SIR differential equations
C

C     dS/dt =
      dPop(1) = - pars(3) * Pop(1) * Pop(2)/pars(1)
C     dI/dt =
      dPop(2) =   pars(3) * Pop(1) * Pop(2)/pars(1)
     $          - Pop(2)/pars(2)
C     dR/dt =
      dPop(3) =   Pop(2)/pars(2)
c cummulative dsdt
      dPop(4) =  - dPop(1)

      RETURN
      END subroutine derivSIR



!
! -------------------------------------------------------------------
!

      subroutine RK4SEIROneD(sh,school,ndata,ndays,nstep,tps,
     $     nparam,curpars,dsdt) 


      implicit none
      integer iday_per_week, nc, np
      parameter (iday_per_week=7, nc = 5, np=4)
      integer iseed
      integer nparam, ndata, ndays, nstep
      real*8 sh(ndata), school(ndata)
      real*8 curpars(nparam)
      real*8 fN, R0, t0, Tg, sigma
      real*8 alpha, Aparam, delta, deltaR, dur, ts, seed
      real*8 Rt, sSh, sEh, sIh, sRh
      real*8 dsdt((ndays)*nstep)
      real*8 dt,tps(ndata)     !time-series
      real*8 tps2(0:(ndata+1))
      real*8 t_cur, p2, p3, p6
      integer istep, iday, iweek
      integer icount, istart, iend

      real*8 param(np)
      real*8 Y(nc), tmpY(nc)
      real*8 dPop1(nc), dPop2(nc), dPop3(nc), dPop4(nc)


      dt = 1.0/ dble(nstep)

      p2 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0
!     
!     ran1 works with -iseed
!     
      iseed = -abs(iseed)

      fN = curpars(1)
      Tg  = curpars(2)
      R0 = curpars(3)
      sigma = curpars(4)
      t0 = curpars(6)
      seed = curpars(7)
      deltaR = curpars(9)
      Aparam = curpars(10)
      alpha  = curpars(11)
      delta = curpars(12)
      ts = curpars(13) 
      dur = curpars(14)

      tps2 = 0.0d0
      tps2(1:ndata) = tps
      tps2(0) = tps2(1) - (tps(2) - tps(1)) 
      tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 

      t_cur = tps2(0)

      istart = int(tps2(0))
      iend   = int(tps2((ndata+1)))


      sSh = fN - seed
      sEh = seed
      sIh = 0.0
      sRh = 0.0


!     initialize the vectors SEIR and cumulative

      Y(1) = sSh
      Y(2) = sEh
      Y(3) = sIh
      Y(4) = sRh
      y(5) = 0.0d0

      dsdt = 0.0d0
 
      t_cur = tps2(0)

!     The equations use beta and not R(t)

      param(1) = fN 
      param(2) = Tg
      param(3) = R0 / Tg
      param(4) = sigma

      icount = 0
      do iday=istart, iend
          iweek = (iday-int(tps2(0)))/iday_per_week + 1
          
          iweek = min(iweek,ndata)
          iweek = max(1,iweek)

          do istep=1,nstep

            t_cur = t_cur+dt
            if (t_cur < t0) go to 101

! calculate the time dependent R(t)
             Rt = R0 * 
     $            * (1.0d0+deltaR*exp(-Aparam*sh(iweek)))
     $            * (1.0d0-alpha*school(iweek))
             if (t_cur >= ts .and. t_cur <= (ts+dur)) Then
                Rt = Rt * (1.0d0 + delta)
             endif
             
             param(3) = Rt/Tg  

            call derivSEIR(Param, Y, dpop1, np, nc)

            tmpY=Y+dt*dPop1 * p2
               
            call derivSEIR(Param, tmpY, dpop2, np, nc)

            tmpY=Y+dt*dPop2 * p2

            call derivSEIR(Param, tmpY, dpop3, np, nc)

            tmpY=Y+dt*dPop3

            call derivSEIR(Param, tmpY ,dpop4, np, nc)

            tmpY=Y+dt*(dPop1*p6 + dPop2*p3 +
     $           dPop3*p3 + dPop4*p6)
            Y = tmpY

!     Update -dS/dt
 101        continue
            icount = icount + 1
            dsdt(icount) = y(nc)


         enddo                  ! End of a single day         
      enddo                     ! End of loop over all days
 
      return
      end subroutine RK4SEIROneD

!------------------------------------
C The Main Differential Equations.

      SUBROUTINE derivSEIR(Param, Pop, dPop, np, nc)

      implicit real*8(a-h,o-z)

      integer np, nc
      REAL*8 Param(np), Pop(nc), dPop(nc)

C
C     The differential equations-including source/sink term
C

C     dS/dt =
      dPop(1) = - Param(3)*Pop(1)*Pop(3)/Param(1)

C     dE/dt 
      dPop(2) = Param(3)*Pop(1)*Pop(3)/Param(1) - param(4) * Pop(2)
C     dI/dt =
      dPop(3) = param(4) * Pop(2) - Pop(3)/Param(2)
C     dR/dt =
      dPop(4) = Pop(3)/Param(2) 
c cummulative sigma * E 
      dPop(5) = Param(4)*Pop(2)
      RETURN
      END SUBROUTINE derivSEIR

C--------------------------------------------------------------------------------

      subroutine BuildIMID(ndata,nstep,tps,imid)

      implicit none

      integer ndata,nstep
      integer imid((ndata+1))
      real*8 tps(ndata)
      integer i
      
! Build the vector needed for calculating the incidence
! days is needed only to calculate imid
! tps is the 'cumulative' day number
! This works correctly provided tps holds cumulative days numbers
! that start with the first month or week and not some arbitrary start

! This line works for the dengue monthly/weekly data but not for the cdc data
! For now we have this somewhat dirty solution


      imid(1) = int(tps(1)) * int(nstep * 0.50)

      do i = 1, ndata
         imid((i+1)) = imid(1) + int(tps(i) * nstep)
      enddo

      return
      end subroutine BuildIMID

C--------------------------------------------------------------------------------


      subroutine RK4VecSIROneD(ndata,ndays,nstep,tps,
     $     nparam,curpars,dsdt)


      implicit none

      integer iday_per_week, np, nc
      parameter (iday_per_week = 7, np = 7, nc = 6)
      integer ndata, nstep, ndays, nparam
      real*8 curpars(nparam)
      real*8 deltaR, Aparam, alpha
      real*8 delta, ts, dur
      real*8 dt
      real*8 t0, R0, Tg, sigma
      real*8 fN
      real*8 bite, vec_k, muV, t_HV, t_VH, sigmaV
      real*8 sSh, sIh, sRh, sSv, sIv
      real*8 dsdt((ndays)*nstep)
      real*8 seed !initial number of infectious
      real*8 tps(ndata) !time-series
      real*8 tps2(0:(ndata+1))
      real*8 y(6), tmpY(6)
      real*8 dy1(6), dy2(6), dy3(6), dy4(6)
      real*8 t_cur
      real*8 p5, p3, p6
      real*8 pars(np)

      integer istep, iday, iweek
      integer icount, istart, iend

      dt = 1.0d0 /dble(nstep)

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0


      fN = curpars(1)
      Tg  = curpars(2)
      R0 = curpars(3)
      sigma = curpars(4)
      t0 = curpars(6)
      seed = curpars(7)
      deltaR = curpars(9)
      Aparam = curpars(10)
      alpha  = curpars(11)
      delta = curpars(12)
      ts = curpars(13) 
      dur = curpars(14)



! Number of Bites per day
      bite  = curpars(15) 
! ratio of vector/human 
      vec_k = curpars(16)
! vector birth/death rate
      muV   = curpars(17)
! Probability that an infected vector when biting a human will transmit the infection
      t_HV  = curpars(18)
! Probability of transmission in reverse direction
      t_VH  = curpars(19)
! vector rate of moving from exposed to infectious 
      sigmaV = curpars(20) 

!! Pack what we need for the ODEs
      pars(1) = fN
      pars(2) = Tg
      pars(3) = bite
      pars(4) = vec_k
      pars(5) = muV
      pars(6) = t_HV
      pars(7) = t_VH

! initialize things - this is the Y vector
! Order is S_h, I_h, R_h, cum state, S_v, I_v

      sSh = fN - seed 
      sIh = seed
      sRh = 0.0d0
      sSv = vec_k * fN - seed
      sIv = seed 

      y(1) = sSh
      y(2) = sIh
      y(3) = sRh
      y(4) = 0.0d0
      y(5) = sSv  
      y(6) = sIv

      dsdt = 0.0d0

      tps2 = 0.0d0
      tps2(1:ndata) = tps
      tps2(0) = tps2(1) - (tps(2) - tps(1)) 
      tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 

c$$$      t_cur = tps(1)
      t_cur = tps2(0)

      icount = 0

!     Loop over ndays and in each day loop over nsteps
!     need to asign each day to a week because the sh and school data are weekly
!     tps has an arbitrary start (could be week 27 for example) but sh/school 
!     go from week 1 to week ndata hence we substract the start

       istart = int(tps2(0))
       iend   = int(tps2((ndata+1)))

       do iday = istart, iend
          iweek = (iday-int(tps2(0)))/iday_per_week + 1
          
          iweek = min(iweek,ndata)
          iweek = max(1,iweek)

          do istep=1,nstep
             t_cur = t_cur+dt

             if (t_cur <= t0) then
                go to 101
             end if
             

             call derivVecSIR(pars, Y, dY1, np, nc)

             tmpY = Y + dt * dY1 * p5

             call derivVecSIR(pars, tmpY, dY2, np, nc)

             tmpY = Y + dt * dY2 * p5

             call derivVecSIR(pars, tmpY, dY3, np, nc)
               
             tmpY = Y + dt * dY3

             call derivVecSIR(pars, tmpY, dY4, np, nc)
               
             tmpY = Y + dt * ( dY1 * p6 + dY2 * p3 +
     $            dY3 * p3 + dY4 * p6 )
             Y = tmpY 
 101         continue
!     Update -dS/dt       
             icount= icount + 1
               
             dsdt(icount) = y(4)

          enddo                 ! End of loop over days
       enddo

      return
      end subroutine RK4VecSIROneD


!------------------------------------
C The Main Differential Equations- (SIR) Host / (SI) Vector  model

      subroutine derivVECSIR(pars, Pop, dPop, np, nc)

      implicit none

      integer np, nc
      REAL*8 pars(np)
      REAL*8 Pop(nc), dPop(nc)

C
C     The SIR differential equations - Human Equations
C

C     dS_H/dt =
      dPop(1) = - pars(3) * pars(6) * Pop(1) * Pop(6)/pars(1)
C     dI_H/dt =
      dPop(2) =   pars(3) * pars(6) * Pop(1) * Pop(6)/pars(1)
     $          - Pop(2)/pars(2)
C     dR_H/dt =
      dPop(3) =   Pop(2)/pars(2)
c cummulative dsdt- Human
      dPop(4) =  - dPop(1)


C
C  The SI - Vector Equations

C     dS_V/dt = 
      dPop(5) = -pars(3) * pars(7) * Pop(5) * Pop(2)/pars(1)
     $     - pars(5) * (pop(5) - pars(4) * pars(1) )
C     dI_V/dt =
      dPop(6) =   pars(3) * pars(7) * Pop(5) * Pop(2)/pars(1)
     $          - pars(5) * pop(6) 
            
      RETURN
      END subroutine derivVECSIR


C--------------------------------------------------------------------------------


      subroutine RK4VecSEIROneD(ndata,ndays,nstep,tps,
     $     nparam,curpars,dsdt)


      implicit none

      integer iday_per_week, np, nc
      parameter (iday_per_week = 7, np = 9, nc = 8)
      integer ndata, nstep, ndays, nparam
      real*8 curpars(nparam)
      real*8 deltaR, Aparam, alpha
      real*8 delta, ts, dur
      real*8 dt
      real*8 t0, R0, Tg, sigmaH
      real*8 fN
      real*8 bite, vec_k, muV, t_HV, t_VH, sigmaV
      real*8 sSh, sEh, sIh, sRh, sSv, sEv, sIv
      real*8 dsdt((ndays)*nstep)
      real*8 seed !initial number of infectious
      real*8 tps(ndata) !time-series
      real*8 tps2(0:(ndata+1))
      real*8 y(nc), tmpY(nc)
      real*8 dy1(nc), dy2(nc), dy3(nc), dy4(nc)
      real*8 t_cur
      real*8 p5, p3, p6
      real*8 pars(np)

      integer istep, iday, iweek
      integer icount, istart, iend

      dt = 1.0d0 /dble(nstep)

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0


      fN = curpars(1)
      Tg  = curpars(2)
      R0 = curpars(3)
      sigmaH = curpars(4)
      t0 = curpars(6)
      seed = curpars(7)
      deltaR = curpars(9)
      Aparam = curpars(10)
      alpha  = curpars(11)
      delta = curpars(12)
      ts = curpars(13) 
      dur = curpars(14)



! Number of Bites per day
      bite  = curpars(15) 
! ratio of vector/human 
      vec_k = curpars(16)
! vector birth/death rate
      muV   = curpars(17)
! Probability that an infected vector when biting a human will transmit the infection
      t_HV  = curpars(18)
! Probability of transmission in reverse direction
      t_VH  = curpars(19)
! vector rate of moving from exposed to infectious 
      sigmaV = curpars(20) 

!! Pack what we need for the ODEs
      pars(1) = fN
      pars(2) = Tg
      pars(3) = bite
      pars(4) = vec_k
      pars(5) = muV
      pars(6) = t_HV
      pars(7) = t_VH
      pars(8) = sigmaH
      pars(9) = sigmaV 


      sSh = fN - seed
      sEh = seed
      sIh = 0.0
      sRh = 0.0
      
      sSv = vec_k * fN - seed
      sEv = seed
      sIv = 0.0d0

! initialize things - this is the Y vector
! Order is S_h, I_h, R_h, cum state, S_v, I_v
      y(1) = sSh
      y(2) = sEh
      y(3) = sIh
      y(4) = sRh
      y(5) = 0.0d0
      y(6) = sSv 
      y(7) = sEv
      y(8) = sIv

      dsdt = 0.0d0

      tps2 = 0.0d0
      tps2(1:ndata) = tps
      tps2(0) = tps2(1) - (tps(2) - tps(1)) 
      tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 

c$$$      t_cur = tps(1)
      t_cur = tps2(0)

      icount = 0

!     Loop over ndays and in each day loop over nsteps
!     need to asign each day to a week because the sh and school data are weekly
!     tps has an arbitrary start (could be week 27 for example) but sh/school 
!     go from week 1 to week ndata hence we substract the start

       istart = int(tps2(0))
       iend   = int(tps2((ndata+1)))

       do iday = istart, iend
          iweek = (iday-int(tps2(0)))/iday_per_week + 1
          
          iweek = min(iweek,ndata)
          iweek = max(1,iweek)

          do istep=1,nstep
             t_cur = t_cur+dt

             if (t_cur <= t0) then
                go to 101
             end if
             

             call derivVecSEIR(pars, Y, dY1, np, nc)

             tmpY = Y + dt * dY1 * p5

             call derivVecSEIR(pars, tmpY, dY2, np, nc)

             tmpY = Y + dt * dY2 * p5

             call derivVecSEIR(pars, tmpY, dY3, np, nc)
               
             tmpY = Y + dt * dY3

             call derivVecSIR(pars, tmpY, dY4, np, nc)
               
             tmpY = Y + dt * ( dY1 * p6 + dY2 * p3 +
     $            dY3 * p3 + dY4 * p6 )
             Y = tmpY 
 101         continue
!     Update -dS/dt       
             icount= icount + 1
               
             dsdt(icount) = y(5)

          enddo                 ! End of loop over days
      enddo
       return
       end subroutine RK4VecSEIROneD

C
C-----------------------------------------------------------
C The Main Differential Equations- (SEIR) Host / (SEI) Vector  model

      subroutine derivVECSEIR(pars, Pop, dPop, np, nc)

      implicit none

      integer np, nc
      REAL*8 pars(np)
      REAL*8 Pop(nc), dPop(nc)

C
C     The SEIR differential equations - Human Equations
C

C     dS_H/dt =
      dPop(1) = - pars(3) * pars(6) * Pop(1) * Pop(8)/pars(1)
C     dE_H/dt = 
      dPop(2) = pars(3) * pars(6) * Pop(1) * Pop(8)/pars(1) 
     $     - pars(8) * pop(2)
C     dI_H/dt =
      dPop(3) =   pars(8) * pop(2) - pop(3) / pars(2)
C     dR_H/dt =
      dPop(4) =   Pop(3)/pars(2)
c cummulative dsdt- Human
      dPop(5) =  pars(8) * pop(2)

C
C  The SEI - Vector Equations

C     dS_V/dt = 
      dPop(6) = -pars(3) * pars(7) * Pop(6) * Pop(3)/pars(1)
     $     - pars(5) * (pop(6) - pars(4) * pars(1) )
C     dE_V/dt =
      dpop(7) = pars(3) * pars(7) * Pop(6) * Pop(3)/pars(1)
     $    - pop(7) * (pars(9) + pars(5)) 
C     dI_V/dt =
      dPop(8) = pars(9) * pop(7) - pars(5) * pop(8)
      

      RETURN
      END subroutine derivVECSEIR


