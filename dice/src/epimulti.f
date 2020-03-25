      subroutine epimulti(n,epi, gamaEpi, sh,school,wght,nparam,
     $     par,parmin,parmax, step,ilog,imask,iseed,nsamps ,
     $     ithin,ndata,tps,rtn,curMin,coef,
     $     parCPL,pminCPL,pmaxCPL,stepCPL,ilogCPL,imaskCPL,fRij,tab,
     $     tabCPL,profiles,nRnd, ndays,epi_model, Temp) 


      implicit none


      integer iday_per_week, nblock
      parameter(iday_per_week=7,nblock=3)

      integer nCPL
      parameter(nCPL=2)
      integer nstep
      parameter (nstep = 6)
! Y and gamay hold the synthetic profiles
      integer n, nparam, iseed
      integer ndays, ithin, ndata, nsamps
      integer epi_model, nRnd
      integer ilog(nparam),imask(nparam)
      real*8 epi(ndata,n), gamaEpi(ndata,n)
      real*8 sh(ndata,n), school(ndata,n)
      real*8 wght(ndata,n)
      real*8 par(nparam,n)
      real*8 parmin(nparam,n), parmax(nparam,n),step(nparam,n)
      real*8 Temp
!!    real*8 ymu(nparam,n), sigma(nparam,n)      
      real*8 fN(n)
      real*8 pC(n)
      real*8 tps(ndata)
      real*8 e_bckgrnd(n)
      real*8 rtn(ndata,n)
      real  tab(nsamps/ithin,nparam+1,n)
      real  tabCPL(nsamps/ithin,nCPL)
      real  profiles(nRnd, ndata, n)
      real*8 rtnBest(ndata,n),rtnNew(ndata,n)
      real*8 dsdt(((ndays)*nstep),n)
      real*8 curpars(nparam,n),copypar(nparam,n),savepar(nparam,n)
      real*8 parupdt(nparam,n),parBest(nparam,n)
      real*8 curLLK(n), curMin(n)
      real*8 fnewLLK(n)
      real*8 ran1
      real*8 accept_rate(nblock)
      logical plog(nparam)
      integer i, j, k, m, icount
      integer noptv, ioptv(nparam)
      integer nopt(nblock),iopt(nparam,nblock)
      integer iaccept(nblock)
      integer iadapt(nblock)
      real*8 AUX
      real*8  range_min, range_max, scale, myaccept
      real*8 step_max, step_min
      real*8 tmp((nblock+1))
      integer ionep
      integer imid((ndata+1))
! Unique to the coupling between regions
      real*8 fRij(n,n), fMij(n,n), coef(n)
      real*8 parCPL(nCPL),pminCPL(nCPL)
      real*8 pmaxCPL(nCPL),stepCPL(nCPL)
      real*8 savefMij(n,n),parBestCPL(nCPL)
      real*8 saveparCPL(nCPL), copyparCPL(nCPL)
      real*8 curparsCPL(nCPL),parupdtCPL(nCPL)
      real*8 fnewLLK1D, curLLK1D,curMin1D
!      real*8 rtn1D(ndata),rtnNew1D(ndata)

      integer imaskCPL(nCPL), ilogCPL(nCPL)
      integer ioptCPL(nCPL), noptCPL
      integer iacceptCPL, iadaptCPL
      logical plogCPL(nCPL)

      real*8 tempCPL(sum(imaskCPL))
      integer iburn, iline, nlines
      external calcFit,ran1


      nlines = nsamps / ithin
! For an adaptive size MCMC - decide if we need to update step size every 1% 
! or 1,000 steps
      scale = 2.0d0
      ionep = int(nsamps * 0.01)
      ionep = min(1000, ionep)
      range_min = 0.20d0
      range_max = 0.30d0
      step_max = 1.d0
      step_min = 1e-5

 
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

!
! Separate R0 and pC into their own block.  It helps..
!
      ilog(3) = 2
      ilog(5) = 2

      do i = 1,nparam
         if (imask(i) > 0) Then
            k = ilog(i) + 1
            nopt(k) = nopt(k) + 1
            iopt(nopt(k),k) = i
            do j = 1, n
            if (par(i,j) .le. parmin(i,j) .or. par(i,j) 
     $              .ge. parmax(i,j)) Then
               par(i,j) = parmin(i,j) +
     $              0.50d0 * (parmax(i,j) - parmin(i,j))

            endif
         enddo
         endif
      enddo

      icount = 0

! convert 0 and 1 to true and false

      do i=1,nparam
         if (ilog(i) .eq. 0) Then
            plog(i) = .false.
         else
            plog(i) = .true.
         endif
      enddo


! For the coupling parameters
      noptCPL = 0
      do i =1,nCPL
         if(imaskCPL(i) > 0) Then
            noptCPL=noptCPL+1
            ioptCPL(noptCPL) = i
           if (parCPL(i) .le. pminCPL(i) .or. parCPL(i) .ge. pmaxCPL(i))
     $           Then
               parCPL(i) = pminCPL(i) + 
     $         0.50d0 * (pmaxCPL(i) - pminCPL(i))
            endif
         endif
         if(ilogCPL(i) .eq. 0) Then
            plogCPL(i) = .false.
         else
            plogCPL(i) = .true.
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


! the paramters are in the following order
! t0, pC, R0, seed, 

      curpars = par

! integrate the ODEs to get the profiles

      fN = curpars(1, 1:n)
      pC = curpars(5, 1:n)
      e_bckgrnd = curpars(8, 1:n)


      fMij = 0.0d0
      tab = 0.0d0
      tabCPL = 0.0d0

 
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



! calculate Likelihood of the solution

      call calcFit(n,epi,gamaEpi,rtn,wght,ndata,curLLK)

      curLLK1D = dot_product(curLLK, coef) /dble(n *Temp)
     
      curMIN = curLLK
      curMIN1D = curLLK1D
      rtnBest = rtn

      iaccept = 0
      iadapt = 0
      icount = 0

      iacceptCPL = 0
      iadaptCPL = 0
      
      parupdt = curpars
      parupdtCPL = parCPL
      curparsCPL = parCPL
      
      tab = 0.0d0
      tabCPL = 0.0d0

! MCMC loop starts here 
      do i =1, nsamps
 
         do k = 1, nblock

            if (nopt(k) .eq. 0) go to 102
            savepar = curpars
            copypar = curpars
    

            call fnProposeParamUpdatesBlock(n,nparam,copypar,
     $           parmin,parmax,step,plog,parupdt,
     $           iseed,nopt(k),iopt(1:nopt(k),k))
            
!     update all the parameters and calculate a new solution and LLK

            curpars(iopt(1:nopt(k),k),:)=parupdt(iopt(1:nopt(k),k),:)
            
            fN = curpars(1, 1:n)
            pC = curpars(5, 1:n)
            e_bckgrnd = curpars(8, 1:n) 
           
            select case (epi_model) 
            case (1) 
              call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
            case (2)                
               call RK4SEIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
            case default 
               call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $              curpars,fMij,dsdt)
            end select 
    

            rtnNew = 0.0d0
            call weekly(n, ndata, ndays, nstep, imid, dsdt, pC,
     $           e_bckgrnd,rtnNew)

            call calcFit(n, epi,gamaEpi,rtnNew,wght,ndata,
     $           fnewLLK)

            fnewLLK1D = dot_product(fnewLLK, coef)/dble(n * Temp)

            AUX = 1.0d0
c$$$            if (k .eq. 3) Then
c$$$               AUX = Ratios(nparam,iopt(1:nopt(k),k),
c$$$     $              curpars,savepar,ymu,sigma,Temp(1))
c$$$            endif
            
            call MH(n,fnewLLK1D,curLLK1D,AUX,curMin1D,nparam,ndata,
     $           iseed,iaccept(k),iadapt(k),curpars,
     $           savepar,parBest,rtnNew,rtnBest,fnewLLK,curLLK,curMin)

 102        continue

         enddo ! End of loop over blocks 
         

!! Now updated the coupling parameters 
         if (noptCPL .eq. 0) go to 103
         savefMij = fMij
         saveparCPL = curparsCPL
         copyparCPL = curparsCPL

         call fnProposeParamUpdates(nCPL,copyparCPL,
     $        pminCPL,pmaxCPL,stepCPL,plogCPL,
     $        parupdtCPL,iseed,noptCPL,ioptCPL(1:noptCPL))


         curparsCPL(ioptCPL(1:noptCPL))=
     $        parupdtCPL(ioptCPL(1:noptCPL))

         fN = curpars(1, 1:n)
         pC = curpars(5, 1:n)
         e_bckgrnd = curpars(8, 1:n)

         call buildMij(n,nCPL,fN,fRij,curparsCPL,fMij)
         select case (epi_model)
         case (1) 
             call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $            curpars,fMij,dsdt)
         case (2)                
            call RK4SEIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $           curpars,fMij,dsdt)
         case default 
            call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $           curpars,fMij,dsdt)
         end select 

         rtnNew = 0.0d0
         call weekly(n, ndata, ndays, nstep, imid, dsdt, pC,e_bckgrnd,
     $        rtnNew)
            
         call calcFit(n, epi,gamaEpi,rtnNew,wght,ndata,
     $        fnewLLK)

         fnewLLK1D = dot_product(fnewLLK , coef) /dble(n * Temp)

         call MH_CPL(n,fnewLLK1D,curLLK1D,AUX,curMin1D,nCPL,ndata,
     $        iseed,iacceptCPL,iadaptCPL,curparsCPL,
     $        saveparCPL,parBestCPL,fMij,savefMij,rtnNew,rtnBest,
     $        fnewLLK,curLLK,curMin)
         
 103     continue

         
!! End of adjusting coupling parameters
         
         if (mod(i,ithin) .eq. 0) Then
            icount = icount + 1
            tab(icount,1:nparam,:)   = real(curpars)
            tab(icount,(nparam+1),:) = real(curLLK)
            tabCPL(icount,:) = real(curparsCPL)
         endif
         
         if (mod(i,ionep) .eq. 0) Then
            m = 0
            do k = 1,nblock
               if (nopt(k) .eq. 0) go to 202
               myaccept = dble(iadapt(k))/dble(ionep)

               if (myaccept > range_max .and.
     $         all(step(iopt(1:nopt(k),k),:)*scale < step_max)) Then
                 step(iopt(1:nopt(k),k),:) =
     $                 step(iopt(1:nopt(k),k),:) *scale
               endif
               if (myaccept < range_min .and.
     $         all(step(iopt(1:nopt(k),k),:)/scale > step_min)) Then
                 step(iopt(1:nopt(k),k),:) =
     $                 step(iopt(1:nopt(k),k),:) /scale
               endif
               iadapt(k) = 0    !reset acceptance number
               m = m + 1
               tmp(m) = dble(iaccept(k))/dble(i)*100.0
 202           continue
               
            enddo
            
!! Repeat for coupling
            if (noptCPL. eq. 0) go to 302
            myaccept = dble(iadaptCPL)/dble(ionep)

            if (myaccept > range_max) Then
              tempCPL = stepCPL(ioptCPL)
              where (tempCPL*scale < step_max)
               tempCPL = tempCPL*scale
              end where
              stepCPL(ioptCPL) = tempCPL
            endif
            if (myaccept < range_min ) Then
              tempCPL = stepCPL(ioptCPL)
              where (tempCPL/scale > step_min)
                tempCPL = tempCPL/scale
              end where
              stepCPL(ioptCPL) = tempCPL
            endif

            iadaptCPL = 0       !reset acceptance number

            m = m + 1
            tmp(m) = dble(iacceptCPL)/dble(i)*100.0

 302        continue
            call intpr('MCMC Step Number:',-1,i,1)
            call dblepr('curLLK1D',-1,curLLK1D,1)
            call dblepr('accept%',-1,tmp,m)
         endif

      enddo
      
!update change information only every ithin iterations 

      accept_rate = dble(iaccept)/dble(nsamps) *100.0d0
            
      call dblepr('acceptence %',-1,accept_rate,nblock)

! generate the profiles _ could really be done on the fly in the MCMC procedure
! see aove 

      iburn = nlines/5
      do i = 1, nRnd
         iline = iburn + floor(ran1(iseed)*(nlines - iburn)) 
         curpars = tab(iline, 1:nparam, 1:n)
         parCPL  = tabCPL(iline,:)

         fN = curpars(1, 1:n)
         pC = curpars(5, 1:n)
         e_bckgrnd = curpars(8, 1:n)


         fMij = 0.0d0
 
         call buildMij(n,nCPL,fN,fRij,parCPL,fMij)
      

         select case (epi_model) 
         case (1) 
            call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $           curpars,fMij,dsdt)
         case (2) 
            call RK4SEIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $           curpars,fMij,dsdt)
         case default 
            call RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $           curpars,fMij,dsdt)
         end select 
      
         call weekly(n,ndata, ndays, nstep, imid, dsdt, pC,
     $        e_bckgrnd, rtn)
         
         profiles(i, 1:ndata, 1:n) = real(rtn)
      enddo


      rtn = rtnBest
      par = parBest
      
      return
         
      end subroutine epimulti


! ------------------------------------------------------------------------

      subroutine RK4SIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $     curpars,fMij,dsdt)


      implicit none
      integer iday_per_week 
      parameter(iday_per_week=7)
      integer n, ndata, ndays, nstep, nparam
      real*8 sh(ndata,n), school(ndata,n)
      real*8 curpars(nparam,n)
      real*8 deltaR(n), Aparam(n), alpha(n)
      real*8 delta(n), ts(n), dur(n)
      real*8 dt
      real*8 t0(n), R0(n), Tg, sigma(n)
      real*8 fN(n), fMij(n,n)
      real*8 dsdt((ndays)*nstep,n)
      real*8 seed(n) !initial number of infectious
      real*8 tps(ndata) !time-series
      real*8 tps2(0:(ndata+1))
      real*8 y(4,n), tmpY(4,n)
      real*8 dy1(4,n), dy2(4,n), dy3(4,n), dy4(4,n)
      real*8 lambda(n)
      real*8 t_cur, Rt(n)
      real*8 p5, p3, p6
      integer np(n), index
      integer k, iweek, istep, iday
      integer icount, istart, iend
      
      dt = 1.0d0 / dble(nstep)

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0
  
      fN = curpars(1,:)
      Tg  = curpars(2,1)
      R0 = curpars(3,:)
      sigma = curpars(4,:)
!pC = curpars(5,:)
      t0 = curpars(6,:)
      seed = curpars(7,:)
!e_bckgrnd = curpars(8,:)
      deltaR = curpars(9,:)
      Aparam = curpars(10,:)
      alpha  = curpars(11,:)
      delta = curpars(12,:)
      ts = curpars(13,:) 
      dur = curpars(14,:)

! initialize things - this is the Y vector

      y(1,1:n) = fN - seed !fN * fraction-1.0
      y(2,1:n) = seed
      y(3,1:n) = fN   !* (1.0d0 - fraction)
      y(4,1:n) = 0.0d0

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
             index = 0
             dY1 = 0.0d0
             dY2 = 0.0d0
             dY3 = 0.0d0
             dY4 = 0.0d0
             tmpY = 0.0d0
! calculate the time dependent R(t)
             Rt = R0  
     $            * (1.0d0+deltaR*exp(-Aparam*sh(iweek,1:n)))
     $            * (1.0d0-alpha*school(iweek,1:n))

             do k = 1,n
                if (t_cur > t0(k)) then
                   index = index+1
                   np(index) = k
                end if

                if (t_cur >= ts(k) .and. t_cur <= (ts(k)+dur(k))) Then
                   Rt(k) = Rt(k) * (1.0d0 + delta(k))
                endif
             enddo

             Rt = Rt / Tg

             call calcForce(n,y(2,:),fN,fMij,Rt,lambda)

             call diffSIR(n, np,index, lambda, Tg, Y, dY1)

             tmpY = Y + dt * dY1 * p5
 
             call calcForce(n,tmpY(2,:),fN,fMij,Rt,lambda)

             call diffSIR(n, np,index, lambda, Tg, tmpY, dY2)

             tmpY = Y + dt * dY2 * p5
             
             call calcForce(n,tmpY(2,:),fN,fMij,Rt,lambda)

             call diffSIR(n, np,index, lambda, Tg, tmpY, dY3)

             tmpY = Y + dt * dY3
 
             call calcForce(n,tmpY(2,:),fN,fMij,Rt,lambda)

             call diffSIR(n, np,index, lambda, Tg, tmpY, dY4)

             tmpY = Y + dt * ( dY1 * p6 + dY2 * p3 +
     $              dY3 * p3 + dY4 * p6 )
             Y = tmpY 

!     Update -dS/dt       
             icount= icount + 1
               
             dsdt(icount,1:n) = y(4,1:n)
          enddo                 ! End of a single day 
        enddo                    ! End of all days

       return
      end subroutine RK4SIR


! ------------------------------------------------------------------------

      subroutine RK4SEIR(n,sh,school,ndata,ndays,nstep,tps,nparam,
     $     curpars,fMij,dsdt)


      implicit none
      integer iday_per_week 
      parameter(iday_per_week=7)
      integer n, ndata, ndays, nstep, nparam
      real*8 sh(ndata,n), school(ndata,n)
      real*8 curpars(nparam,n)
      real*8 deltaR(n), Aparam(n), alpha(n)
      real*8 delta(n), ts(n), dur(n)
      real*8 dt
      real*8 t0(n), R0(n), Tg, sigma(n)
      real*8 fN(n), fMij(n,n)
      real*8 dsdt((ndays)*nstep,n)
      real*8 seed(n) !initial number of infectious
      real*8 tps(ndata) !time-series
      real*8 tps2(0:(ndata+1))
      real*8 y(5,n), tmpY(5,n)
      real*8 dy1(5,n), dy2(5,n), dy3(5,n), dy4(5,n)
      real*8 lambda(n)
      real*8 t_cur, Rt(n)
      real*8 p5, p3, p6
      integer np(n), index
      integer k, iweek, istep, iday
      integer icount, istart, iend
      
      dt = 1.0d0 / dble(nstep)

      p5 = 1.0d0 / 2.0d0
      p3 = 1.0d0 / 3.0d0
      p6 = 1.0d0 / 6.0d0
  
      fN = curpars(1,:)
      Tg  = curpars(2,1)
      R0 = curpars(3,:)
      sigma = curpars(4,:)
!pC = curpars(5,:)
      t0 = curpars(6,:)
      seed = curpars(7,:)
!e_bckgrnd = curpars(8,:)
      deltaR = curpars(9,:)
      Aparam = curpars(10,:)
      alpha  = curpars(11,:)
      delta = curpars(12,:)
      ts = curpars(13,:) 
      dur = curpars(14,:)

! initialize things - this is the Y vector S-E-I-R-cummulative

      y(1,1:n) = fN - seed !fN * fraction-1.0
      y(2,1:n) = seed
      y(3,1:n) = 0.0d0   !* (1.0d0 - fraction)
      y(4,1:n) = 0.0d0
      y(5,1:n) = 0.0d0
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
             index = 0
             dY1 = 0.0d0
             dY2 = 0.0d0
             dY3 = 0.0d0
             dY4 = 0.0d0
             tmpY = 0.0d0
! calculate the time dependent R(t)
             Rt = R0  
     $            * (1.0d0+deltaR*exp(-Aparam*sh(iweek,1:n)))
     $            * (1.0d0-alpha*school(iweek,1:n))

             do k = 1,n
                if (t_cur > t0(k)) then
                   index = index+1
                   np(index) = k
                end if

                if (t_cur >= ts(k) .and. t_cur <= (ts(k)+dur(k))) Then
                   Rt(k) = Rt(k) * (1.0d0 + delta(k))
                endif
             enddo

             Rt = Rt / Tg

             call calcForce(n,y(3,:),fN,fMij,Rt,lambda)

             call diffSEIR(n,np,index, lambda, sigma, Tg, Y, dY1)

             tmpY = Y + dt * dY1 * p5
 
             call calcForce(n,tmpY(3,:),fN,fMij,Rt,lambda)

             call diffSEIR(n,np,index, lambda, sigma, Tg, tmpY, dY2)

             tmpY = Y + dt * dY2 * p5
             
             call calcForce(n,tmpY(3,:),fN,fMij,Rt,lambda)

             call diffSEIR(n,np,index, lambda, sigma, Tg, tmpY, dY3)

             tmpY = Y + dt * dY3
 
             call calcForce(n,tmpY(3,:),fN,fMij,Rt,lambda)

             call diffSEIR(n,np,index, lambda, sigma, Tg, tmpY, dY4)

             tmpY = Y + dt * ( dY1 * p6 + dY2 * p3 +
     $              dY3 * p3 + dY4 * p6 )
             Y = tmpY 

!     Update -dS/dt       
             icount= icount + 1
               
             dsdt(icount,1:n) = y(5,1:n)
          enddo                 ! End of a single day         
        enddo                    ! End of all days

       return
      end subroutine RK4SEIR


           
!-------------------------------------------------------
           
!-------------------------------------------------------
C Calcluate the force of infection for each compartment

      subroutine calcForce(n,fI,fN,fMij,beta,lambda) 

      implicit none

      integer n
      real*8 fI(n), fN(n), fMij(n,n)
      real*8 beta(n), lambda(n)
      real*8 sum_numer,sum_denom, ratio(n)
      integer i, j, k

      lambda = 0.0d0


      do j=1,n
         sum_numer = 0.0d0
         sum_denom = 0.0d0

         do k=1,n

            sum_numer = sum_numer + fMij(k,j) * fI(k)
            sum_denom = sum_denom + fMij(k,j) * fN(k)
         enddo
         ratio(j) = sum_numer/sum_denom
      enddo

      do i = 1,n
         do j = 1,n 
            lambda(i) = lambda(i) + beta(j) * fMij(i,j) *
     $                 ratio(j)
         enddo

      enddo

      return
      end subroutine calcForce



!------------------------------------
C The Main Differential Equations- SIR

      SUBROUTINE DiffSIR(n, np,index, lambda, Tg, Pop, dPop)

      implicit none

      integer n, np(n), index
      REAL*8 lambda(n), Tg
      REAL*8 Pop(4,n), dPop(4,n)
      integer i, j
C
C     The differential equations
C

C     dS/dt =
      do i = 1,index
         j = np(i)
         dPop(1,j) = - lambda(j)*Pop(1,j) 
C     dI/dt =
         dPop(2,j) =   lambda(j)*Pop(1,j)
     $        - Pop(2,j)/Tg
C     dR/dt =
         dPop(3,j) =   Pop(2,j)/Tg 
c cummulative dsdt
         dPop(4,j) = lambda(j)*Pop(1,j)
      enddo
      RETURN
      END subroutine DiffSIR


!------------------------------------
C The Main Differential Equations- SIR

      SUBROUTINE DiffSEIR(n,np,index, lambda, sigma, Tg, Pop, dPop)

      implicit none

      integer n, np(n), index
      REAL*8 lambda(n), sigma(n), Tg
      REAL*8 Pop(5,n), dPop(5,n)
      integer i, j
C
C     The differential equations
C

C     dS/dt =
      do i = 1,index
         j = np(i)
         dPop(1,j) = -lambda(j)*Pop(1,j)
C     dE/dt = 
         dPop(2,J) = lambda(j)*Pop(1,j)
     $        - sigma(j) * Pop(2,j)
C     dI/dt =
         dPop(3,j) = sigma(j)*Pop(2,j)
     $        - Pop(3,j)/Tg
C     dR/dt =
         dPop(4,j) = Pop(3,j)/Tg 
c cummulative dsdt
         dPop(5,j) = sigma(j) * Pop(2,j)
      enddo
      RETURN
      END subroutine DiffSEIR


!--------------------------------------------------------------------------------
           
        subroutine weekly(n,ndata,ndays,nstep,imid,dsdt,pC,e_bckgrnd,x)
        
        implicit none

        integer n,ndata, ndays, nstep
        real*8 dsdt((ndays)*nstep,n)
        real*8 x(ndata,n), pC(n), e_bckgrnd(n)
        integer imid(ndata+1)
        integer i, istart, iend 
        integer iday_per_week, ishift, j

c$$$        do i=1,ndata
c$$$           iend   = imid(i+1)
c$$$           istart = imid(i)
c$$$
c$$$           x(i,1:n) = dsdt(iend,1:n) - dsdt(istart,1:n)
c$$$
c$$$           x(i,1:n) = x(i,1:n)* pC(1:n) +e_bckgrnd(1:n)
c$$$        enddo

        iday_per_week = 7

        ishift = nstep * iday_per_week * 0.50d0
        do i=1,ndata
           do j = 1,n
              x(i,j) = dsdt(i*nstep*iday_per_week+ishift,j) - 
     $             dsdt(i*nstep*iday_per_week-ishift,j)
           enddo
           x(i,1:n) = x(i,1:n)* pC(1:n) +e_bckgrnd(1:n)
        enddo


        return
        end subroutine weekly

C--------------------------------------------------------------------------------
       subroutine calcFit(n,y,gamay,x,wght,ndata,curLLK)

        implicit none

        integer n, ndata, i, j

        real*8 y(ndata,n),x(ndata,n), gamay(ndata,n)
        real*8 wght(ndata,n)
        real*8 xi, yi,val,CurLLK(n)
        

c x is the simulated data
c y is the base profile

C calculate the P(yi,xi)

           do j = 1,n
              curLLK(j) = 0.0d0
              do i=1,ndata
                 
                 yi = y(i,j)
                 xi = x(i,j)
                 
                 val = yi * log(xi) - xi - gamay(i,j)
                 curLLK(j) = curLLK(j) + val * wght(i,j)
                 
              enddo
              curLLK(j) = -curLLK(j) /sum(wght(:,j))
           enddo           
   
        return

        end subroutine CalcFit

!
! -------------------------------------------------------------------
!
      subroutine MH(n,fnewLLK1D,curLLK1D,AUX,curMin1D,nparam,
     $     ndata,iseed,iaccept,iadapt,curpars,
     $     savepar,parBest,rtn,rtnBest,fnewLLK,curLLK,curMin)


      implicit none
      integer n,nparam,ndata,iseed,iaccept,iadapt
      real*8 fnewLLK1D, curLLK1D, AUX,curMin1D
      real*8 curMin(n), curLLK(n), fnewLLK(n)
      real*8 curpars(nparam,n),savepar(nparam,n)
      real*8 parBest(nparam,n)
      real*8 rtn(ndata,n), rtnBest(ndata,n)
      real*8 rnd,ran1,diff_LLK
      external ran1
      logical accept

      diff_LLK = fnewLLK1D - curLLK1D
      
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
         iadapt= iadapt + 1 
         savepar = curpars
         curLLK1D = fnewLLK1D
         curLLK   = fnewLLK

         if (curLLK1D .le. curMin1D) then !keep track of the best profile we have
            parBest = curpars
            curMin1D = curLLK1D
            curMin   = curLLK
            rtnBest = rtn
         endif
!     if step is rejected - restore saved values of parameters
      else
         
         curpars = savepar
         
      endif
         
      return
      end subroutine MH

! ------------------------------------------------------------------------


      subroutine MH_CPL(n,fnewLLK1D,curLLK1D,AUX,curMin1D,nparam,
     $     ndata,iseed,iaccept,iadapt,curpars,savepar,
     $     parBest,curfMij,savefMij,rtn,rtnBest,fnewLLK,curLLK,curMin)


      implicit none
      integer n,nparam,ndata,iseed,iaccept,iadapt
      real*8 fnewLLK1D, curLLK1D, AUX,curMin1D
      real*8 fnewLLK(n), curLLK(n), curMin(n)
      real*8 curpars(nparam),savepar(nparam)
      real*8 parBest(nparam)
      real*8 curfMij(n,n), savefMij(n,n)
      real*8 rtn(ndata,n), rtnBest(ndata,n)
      real*8 rnd,ran1,diff_LLK
      external ran1
      logical accept

      diff_LLK = fnewLLK1D - curLLK1D
            
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
         curLLK1D = fnewLLK1D
         curLLK   = fnewLLK
         savefMij = curfMij

         if (curLLK1D .le. curMin1D) then !keep track of the best profile we have
            parBest= curpars
            curMin1D = curLLK1D
            curMin = fnewLLK
            rtnBest = rtn
         endif
!     if step is rejected - restore saved values of parameters
      else
         
         curpars = savepar
         curfMij = savefMij
         
      endif
         

      return
      end subroutine MH_CPL

! ------------------------------------------------------------------------


      subroutine SR_from_unit_Block_Mult(n,nopt,x,ymin,ymax,
     $     logflag,rtn)
      
      implicit none

      integer n, nopt
      real*8 x(nopt,n),ymin(nopt,n),ymax(nopt,n),rtn(nopt,n)
      integer i, j
      logical logflag


! Block version
! Only allow for log base 10

      if (logflag) Then
            do i = 1,nopt
               do j = 1,n
                  rtn(i,j) = ymin(i,j) * 
     $           10.0**(x(i,j)*(log10(ymax(i,j))-log10(ymin(i,j))))
               enddo
            enddo
      else
         do i=1,nopt
            do j = 1,n
               rtn(i,j) = ymin(i,j) + (ymax(i,j)-ymin(i,j))*x(i,j)
            enddo
         enddo
      endif

      return
      end subroutine SR_from_unit_Block_Mult

c----------------------------------------------------------------

      subroutine SR_to_unit_Block_Mult(n,nopt,y,ymin,ymax,
     $     logflag,rtn)

      implicit none
      integer n,nopt
      real*8 y(nopt,n), ymin(nopt,n),ymax(nopt,n),rtn(nopt,n)
      integer i,j
      logical logflag

! Block version
! Only allow for log base 10

      if (logflag) Then

         do i = 1,nopt
            do j = 1,n
               rtn(i,j) = (log10(y(i,j)) - log10(ymin(i,j))) /
     $              (log10(ymax(i,j))-log10(ymin(i,j)))
            enddo
         enddo
      else
         do i = 1, nopt
            do j = 1, n
               rtn(i,j) = (y(i,j) - ymin(i,j))/(ymax(i,j) - ymin(i,j))
            enddo
         enddo
      endif
 
      return
      end subroutine SR_to_unit_Block_Mult

c----------------------------------------------------------------

      subroutine fnProposeParamUpdatesBlock(n,nparam,curval,
     $     valmin,valmax,step,logflag,
     $     parupdt,iseed,nopt,iopt)

      implicit none
      integer n,nparam, nopt, iopt(nopt)
      real*8 curval(nparam,n),valmin(nparam,n),valmax(nparam,n)
      real*8 parupdt(nparam,n),step(nparam,n)
      real*8 x(nopt,n), rv(nopt,n), rtn(nopt,n)
      real*8 ran1
      integer iseed
      integer j,k
      
      logical logflag(nparam)
      external ran1


      parupdt = curval

      do j = 1, nopt
         do k = 1,n
            rv(j,k) = ran1(iseed)
            rv(j,k) = (rv(j,k) - 0.50d0) * step(iopt(j),k)
         enddo
      enddo


      call SR_to_unit_Block_Mult(n,nopt,curval(iopt,:),valmin(iopt,:),
     $         valmax(iopt,:),logflag(1),x)

      x = x + rv

      do j = 1, nopt
         do k = 1,n
            if (x(j,k) .le. 0.0d0) x(j,k) = 1.0d0 + x(j,k)
            if (x(j,k) .ge. 1.0d0) x(j,k) = x(j,k) - 1.0d0         
         enddo
      enddo

      call SR_from_unit_Block_Mult(n,nopt,x,valmin(iopt,:),
     $     valmax(iopt,:),logflag(1),rtn)

      parupdt(iopt,:) = rtn
      
      return
      end subroutine fnProposeParamUpdatesBlock

!
! -----------------------------------------------------------------
!

      subroutine buildMij(n,nCPL,fN,fRij,parCPL,fMij)

      implicit none

      integer n,nCPL
      real*8 parCPL(nCPL)
      real*8 fMij(n,n),fN(n),fRij(n,n)
      real*8 tmp(n,n),sd, gamma, denom
      integer i, j

      sd    = parCPL(1)
      gamma = parCPL(2)

      tmp = 0.0d0
      tmp = fRij/sd
! The Mills & Riley Kernel 
      do i=1,n
         do j=i,n
            denom = 1.0d0 + tmp(i,j)**gamma
            fMij(i,j) = fN(j)/denom
            fMij(j,i) = fN(i)/denom
         enddo
         fMij(i,:) = fMij(i,:) /sum(fMIJ(i,:))
      enddo


c$$$      fMij = 0.0d0
c$$$      do i = 1, n
c$$$         fMij(i,i) = 1.0d0
c$$$      enddo

      return
      end subroutine buildMij

! ------------------------------------------------------------------

      subroutine checkParam(n,nparam,pmin,pmax,pval)

      implicit none

      integer n, nparam
      real*8 pmin(nparam,n), pmax(nparam,n)
      real*8 pval(nparam,n)
      integer i, j

      do i=1,nparam
         do j = 1,n
            if (pval(i,j) .ge. pmax(i,j)) Then
               call dblepr('Greater than Pmax',-1,pval(i,j),1)
               call intpr('i',-1,i,1)
               call intpr('j',-1,j,1)
               stop
            endif

            if (pval(i,j) .le. pmin(i,j)) Then
               call dblepr('Smaller than Pmin',-1,pval(i,j),1)
               call intpr('i',-1,i,1)
               call intpr('j',-1,j,1)
               stop
            endif

         enddo
      enddo

      return
      end subroutine checkParam
