      subroutine fitplague(epi, gamaEpi,wght,nparam,
     $     par,parmin,parmax,step,ilog,Temp,
     $     imask,iseed,nsamps,ithin,ndata,tps,rtn,accept_rate,
     $     curMin, tab, ndays, profiles, nRnd) 


      implicit none

      integer iday_per_week, nblock
      parameter(iday_per_week=7,nblock=3)
      integer nstep
      parameter (nstep = 6)
! Y and gamay hold the synthetic profiles
      integer iseed, ndata,  nsamps, nparam
      integer ithin, ndays, nRnd
      integer ilog(nparam),imask(nparam)
      real*8 epi(ndata), gamaEpi(ndata)
      real*8 wght(ndata)
      real*8 par(nparam)
      real*8 parmin(nparam), parmax(nparam),step(nparam)
      real*8 Temp
      real*8 tps(ndata)
      real*8 pC,e_bckgrnd
      real*8 rtn(ndata)
      real  tab(nsamps/ithin,nparam+1)
      real  profiles(nRnd, ndata)
      real*8 rtnBest(ndata),rtnNew(ndata)
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


      scaleNdata = Temp
! Number of lines in the tab array

      nlines = nsamps/ithin 

!     For an adaptive size MCMC - decide if we need to update step size ever 1% of steps
      scale = 2.0d0
      ionep = nsamps * 0.01
      ionep = min(1000,ionep)
      range_min = 0.20d0
      range_max = 0.30d0
      step_max = 1.d0
      step_min = 1e-5
      iadapt = 0

!     ran1 works with -iseed

      iseed = -abs(iseed)

!     ilog can have values of 
!     0 - uniform 
!     1 - log uniform
!     2 - Gaussian 

!     set an auxilary array with the indices of the parameters we are optimizing 

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
      
      icount = 0


!     convert 0 and 1 to true and false

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

      call BuildIMID(ndata, nstep, tps, imid)


!     ! the paramters are in the following order


      curpars = par

!     integrate the ODEs to get the profiles

      pC = curpars(5)
      e_bckgrnd = curpars(8)


      call RK4Plague(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)

      
      call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd, 
     $     rtn)

!     calculate Likelihood of the solution

      curLLK = calcFit1D(epi,gamaEpi,rtn,wght,ndata)

      curLLK = (curLLK)/scaleNdata
      call dblepr('curLLK',-1,curLLK,1)

!     MCMC loop starts here 

      curMIN = curLLK
      rtnBest = rtn

      iaccept = 0
      icount = 0

      savestep = step
     
      parupdt = curpars

      do i =1, nsamps
         
         savepar = curpars
         copypar = curpars
 
         do k = 1, nblock

            if (nopt(k) .eq. 0) go to 102
            savepar = curpars
            copypar = curpars
!     Propose a new value for each parameter we are optimizing 
!     half the steps will be small and half will be the given size
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

            call RK4Plague(ndata,ndays,nstep,tps,
     $           nparam,curpars,dsdt)

         
            rtnNew = 0.0d0
            call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd,
     $           rtnNew)

            fnewLLK = calcFit1D(epi,gamaEpi,rtnNew,wght,ndata)
            fnewLLK = fnewLLK / (scaleNdata)

            call MH1D(fnewLLK,curLLK,1.0d0,curMin,nparam,ndata,
     $        iseed,iaccept(k),iadapt(k),curpars,savepar,
     $        parBest,rtnNew,rtnBest)

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
         
         call RK4Plague(ndata,ndays,nstep,tps,
     $        nparam,curpars,dsdt)
         call weekly1D(ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd, 
     $        rtn)
         profiles(i, 1:ndata) = real(rtn)
      enddo

!! End of generating profiles


      rtn = rtnBest
      par = parBest
      
      return
      
      end subroutine fitplague



!
! -------------------------------------------------------------------
!

       subroutine RK4Plague(ndata,ndays,nstep,tps,
     $     nparam,curpars,dsdt) 


       implicit none
       integer  nc, np
       parameter (nc = 5, np=4)
       integer iseed
       integer nparam, ndata, ndays, nstep
       real*8 curpars(nparam)
       real*8 fN, R0, t0, Tg, sigma
       real*8 seed
       real*8 Rt, sSh, sEh, sIh, sRh
       real*8 dsdt((ndays)*nstep)
       real*8 dt,tps(ndata)     !time-series
       real*8 tps2(0:(ndata+1))
       real*8 t_cur, p2, p3, p6
       integer istep, iday
       integer icount, istart, iend
       
       real*8 param(np)
       real*8 Y(nc), tmpY(nc)
       real*8 dPop1(nc), dPop2(nc), dPop3(nc), dPop4(nc)
       

       dt = 1.0d0 /dble(nstep-1)
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

       tps2 = 0.0d0
       tps2(1:ndata) = tps
       tps2(0) = tps2(1) - (tps(2) - tps(1)) 
       tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 
       
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

          do istep=1,nstep
             
             t_cur = t_cur+dt
             if (t_cur < t0) go to 101
             
             call diffPlague(Param, Y, dpop1, np, nc)
             
             tmpY=Y+dt*dPop1 * p2
             
             call diffPlague(Param, tmpY, dpop2, np, nc)

             tmpY=Y+dt*dPop2 * p2
               
             call diffPlague(Param, tmpY, dpop3, np, nc)

             tmpY=Y+dt*dPop3

             call diffPlague(Param, tmpY ,dpop4, np, nc)
               
             tmpY=Y+dt*(dPop1*p6 + dPop2*p3 +
     $            dPop3*p3 + dPop4*p6)
             Y = tmpY

!     Update -dS/dt
 101         continue
               
             icount= icount + 1
               
             dsdt(icount) = y(nc)

          enddo                 ! End of a single day         
       enddo                    ! End of loop over all days
 
       return
       end subroutine RK4Plague


!------------------------------------
C The Main Differential Equations.

      SUBROUTINE diffPlague(Param, Pop, dPop, np, nc)

      implicit none
      integer np, nc
      REAL*8 Param(np), Pop(nc), dPop(nc)

C
C     The differential equations-including source/sink term
C

C     dS/dt =
      dPop(1) = - Param(3)*Pop(1)*Pop(2)/Param(1)

C     dE/dt 
      dPop(2) = Param(3)*Pop(1)*Pop(2)/Param(1) - param(4) * Pop(2)
C     dI/dt =
      dPop(3) = param(4) * Pop(2) - Pop(3)/Param(2)
C     dR/dt =
      dPop(4) = Pop(3)/Param(2) 
c cummulative sigma * E 
      dPop(5) = Param(4)*Pop(2)
      RETURN
      END




