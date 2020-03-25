      subroutine fitsirb(ndata,tps,epi,gamaEpi,wght,nparam,
     $     par,parmin,parmax,step,ilog,imask,iseed,
     $     nsamps,ithin,tab,profiles, profilesB,
     $     nreals, Temp, ndays) 


      implicit none

      integer iday_per_week, day_per_year
      parameter (iday_per_week=7)
      parameter (day_per_year = 365)

      integer nstep
      parameter(nstep = 10)
      integer ndata, ndays, iseed

      real*8 pC, e_bckgrnd
      real*8 dsdt(ndays*nstep)
      real*8 sSh, sIh, sRh
      real*8 dt,tps(ndata)     !time-series
      real*8 rtn(ndata), wght(ndata)
!     Y and gamay hold the synthetic profiles
      integer nsamps, nparam
      integer ithin, imodel
      integer ilog(nparam),imask(nparam)
      real*8 epi(ndata), gamaEpi(ndata)
      real*8 par(nparam)
      real*8 parmin(nparam), parmax(nparam),step(nparam)
      real  tab(nsamps/ithin,nparam+1)
      integer nreals
      real  profiles( ndata,nreals)
      real  profilesB(ndata,nreals)
      real*8 rtnBest(ndata),rtnNew(ndata)
      real*8 Bact(ndata)
      real*8 curpars(nparam),copypar(nparam),savepar(nparam)
      real*8 savestep(nparam),parupdt(nparam),parBest(nparam)
      real*8 curLLK, curMin, bestLLK
      real*8 fnewLLK, diff_LLK
      real*8 calcFit, ran1, rnd
      real*8 fnProposeParamUpdatesSingle
      real*8 accept_rate
      logical plog(nparam), verbose, accept
      integer i, j, m, icount
      integer nopt, iopt(nparam)
      integer iaccept
      integer iadapt
      real*8  range_min, range_max, scale, myaccept
      real*8 step_max, step_min
      real*8 calcFit1D
      real*8 scaleLLK, Temp, Aux
      integer ionep, nsamps5, ifreq, isave
      integer imid((ndata+1))
      external calcFit1D,ran1



      dt = 1.0d0 / dble(nstep - 1)

      scaleLLK = 0.1 / dble(nsamps)
!     For saving nreals profiles from the second half of the MCMC chain
      nsamps5 = nsamps * 3 / 4
      ifreq = (nsamps-nsamps5)/(nreals-1)

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

      nopt = 0
      iopt = 0
      do i=1,nparam
         if (imask(i) > 0) Then
            nopt = nopt + 1
            iopt(nopt) = i
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


!     ! the paramters are in the following order
!     !

      curpars = par

!     integrate the ODEs to get the profiles
! c("NH", "Tg", "birth", 't0', 'pC', "e_bckgrnd",'seed', 'a', 'lamda','growth_loss','e_shedd')
!
      pC = curpars(5)
      e_bckgrnd = curpars(6)

      rtn = 0.0d0
      rtnNew = 0.0d0
      dsdt = 0.0d0
!     

      call BuildIMID(ndata, nstep, tps, imid)

      call RK4SIRBOneD(nparam, curpars,
     $     ndata, ndays,nstep,dt,tps,dsdt, Bact)

      call weekly1D(ndata, ndays, nstep, imid, dsdt, pC, 
     $     e_bckgrnd, rtn)

!     calculate Likelihood of the solution

      curLLK = calcFit1D(epi,gamaEpi,rtn,wght,ndata)

      curLLK = (curLLK)/Temp

      call dblepr('LLK', -1, curLLK, 1)

!     MCMC loop starts here 

      curMIN = curLLK
      rtnBest = rtn

      iaccept = 0
      icount = 0
      isave = 1
      savestep = step
      
      parupdt = curpars

      AUX  = 1.0d0

      do i =1,nsamps
         
         savepar = curpars
         copypar = curpars
 
!     Propose a new value for each parameter we are optimizing 
!     half the steps will be small and half will be the given size
         call fnProposeParamUpdates(nparam,copypar,parmin,
     $        parmax,step,plog,parupdt,iseed,nopt,iopt(1:nopt))
         
!     update all the parameters and calculate a new solution and LLK

         curpars(iopt(1:nopt))= parupdt(iopt(1:nopt))

         pC    = curpars(5)
         e_bckgrnd = curpars(6)

         rtnNew = 0.0d0
         dsdt = 0.0d0
!     

         call RK4SIRBOneD(nparam, curpars,
     $        ndata, ndays, nstep,dt,tps,dsdt, Bact)
         
         rtnNew = 0.0d0
 
         call weekly1D(ndata, ndays, nstep, imid, dsdt, pC, 
     $        e_bckgrnd, rtnNew)

         fnewLLK = calcFit1D(epi,gamaEpi,rtnNew,wght, ndata)
         fnewLLK = (fnewLLk)/Temp

         call MH1D(fnewLLK,curLLK,AUX,curMin,nparam,ndata,
     $        iseed,iaccept,iadapt,curpars,savepar,
     $        parBest,rtnNew,rtnBest)

         if (mod(i,ithin) .eq. 0) Then
            icount = icount + 1
            tab(icount,1:nparam)   = sngl(curpars)
            tab(icount,(nparam+1)) = sngl(curLLK)
         endif

         if (mod(i,ifreq) .eq. 0 .and. i .ge. nsamps5 .and. 
     $        isave .le. nreals) Then
            profiles( 1:ndata,isave) = sngl(rtnNew)
            profilesB(1:ndata,isave) = sngl(Bact)
             isave = isave + 1
         endif


         if (mod(i,ionep) .eq. 0) Then
            myaccept = dble(iadapt)/dble(ionep)
            if (myaccept > range_max .and.
     $           all(step(iopt(1:nopt))*scale < step_max)) Then
               step(iopt(1:nopt))= step(iopt(1:nopt))*scale
            endif
            if (myaccept < range_min .and.
     $           all(step(iopt(1:nopt))/scale > step_min)) Then
               step(iopt(1:nopt))= step(iopt(1:nopt))/scale
            endif
            iadapt = 0          !reset acceptance number

!     Print information to the screen:
            call intpr('MCMC Step Number:',-1,i,1)
            call dblepr('curLLK',-1,curLLK,1)
            call dblepr('accept%',-1,dble(iaccept)/dble(i)*100.0,1)

         endif
         
      enddo

!     update change information only every ithin iterations 

      accept_rate = dble(iaccept)/dble(nsamps) *100.0d0
      
      rtn = rtnBest
      par = parBest
       
      return
      
      end subroutine fitsirb


!
!
! -------------------------------------------------------------------
!


      subroutine RK4SIRBOneD(nparam, curpars,
     $     ndata,ndays, nstep,dt,tps,dsdt,Bact)


      implicit none

      integer np, nc
      parameter (np = 7, nc = 4)
      integer iseed
      integer nparam
      integer ndata, ndays, nstep
      real*8 curpars(nparam)
      real*8 dsdt((ndays)*nstep)

      real*8 dt,tps(ndata)     !time-series
      real*8 tps2(0:(ndata+2))
      real*8 t_cur, p2, p3, p6
      real*8 ran1
      integer i, iweek, istep, iday, j
      integer icount
      logical debug
      real*8 R0
      real*8 Nh, Tg, birth
      real*8 t0, pC, e_bckgrnd, seed
      real*8 a_par, K, lamda, growth_loss, e_shedd
      real*8 delta, ts, tend
      real*8 bact(ndata)
      real*8 param(np)
      real*8 Y(nc), tmpY(nc)
      real*8 dPop1(nc), dPop2(nc), dPop3(nc), dPop4(nc)
      integer istart, iend
      
      debug = .false.

      p2 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0
!     
!     ran1 works with -iseed
!     
      iseed = -abs(iseed)

    
      nH = curpars(1)
      Tg = curpars(2)
      birth = curpars(3)
      t0 = curpars(4)
      pC = curpars(5)
      e_bckgrnd = curpars(6)
      seed      = curpars(7)
      a_par     = curpars(8)
      K         = curpars(9) * NH 
      growth_loss = curpars(10)
      e_shedd     = curpars(11)

      delta       = curpars(12)
      ts          = curpars(13)
      tend        = ts + 7 !curpars(14)

!     initialize the vectors

      Y(1) = Nh - seed
      Y(2) = seed
      Y(3) = 0.0d0
      Y(4) = 0.0d0

      tmpY = Y

      dsdt = 0.0d0


      tps2 = 0.0d0
      tps2(1:ndata) = tps
      tps2(0) = tps2(1) - (tps(2) - tps(1)) 
      tps2((ndata+1)) = tps2(ndata) +  (tps(ndata) - tps((ndata-1))) 

      t_cur = tps2(0)

      istart = tps2(0)
      iend   = tps2((ndata+1))

      icount = 0
      
      param(1) = birth
      param(2) = Nh
      param(3) = a_par
      param(4) = y(3) / (K + y(3))
      param(5) = Tg
      param(6) = growth_loss
      param(7) = e_shedd

      j = 1

      do iday = istart, iend


            do istep=1,nstep


               t_cur = t_cur+dt
               
               if (t_cur < t0) go to 101

               param(3) = a_par 
               if (t_cur .ge. ts .and. t_cur .le. tend) Then
                  param(3) = a_par * (1.0d0+delta)
               end if

               param(4) = y(3) / (K + y(3))
              
               call derivSIRB(Param,   Y,dpop1,np,nc)
               
               tmpY=Y+dt*dPop1 * p2
               
               call derivSIRB(Param,tmpY,dpop2,np,nc)

               tmpY=Y+dt*dPop2 * p2

               call derivSIRB(Param,tmpY,dpop3,np,nc)

               tmpY=Y+dt*dPop3

               call derivSIRB(Param,tmpY,dpop4,np,nc)

               tmpY=Y+dt*(dPop1*p6 + dPop2*p3 +
     $              dPop3*p3 + dPop4*p6)
               Y = tmpY

!     Update -dS/dt
 101           continue
               
               icount= icount + 1
               
               dsdt(icount) = y(4)

            enddo               ! End of a single day  
            if (iday .eq. tps(j) ) Then
               Bact(j) = y(3)
               j = j + 1
            endif
      enddo                      ! End of loop over all days

       return
      end subroutine RK4SIRBOneD


!------------------------------------
C The Main Differential Equations.

      SUBROUTINE DerivSIRB(Param, Pop, dPop, np, nc)

      implicit real*8(a-h,o-z)

      integer np, nc
      REAL*8 Param(np), Pop(nc), dPop(nc)

C
C     The differential equations-including source/sink term
C

C     dS/dt =
      dPop(1) = Param(1)*(Param(2) - Pop(1))
     $     - Param(3) * Param(4) * Pop(1)
C     dI/dt =
      dPop(2) =    Param(3) * Param(4) * Pop(1)
     $          - Pop(2)/Param(5)
     $          - param(1) * pop(2) 
C     dB/dt =
      dPop(3) =  - Pop(3) * Param(6) + Param(7) * Pop(2)
c cummulative dsdt
      dPop(4) =  Param(3) * Param(4) * Pop(1)
      RETURN
      END
