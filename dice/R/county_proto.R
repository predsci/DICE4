#' A function to simulate a flu-like epidemic and be able to output 
#' which sub population infected which other sub-population
county_proto <- function(
    df_pop,
    R0 = 1.8,
    Tg = 2.6,
    p_c = 1,
    trickle = 1 / 300000000,
    kernPower = 0,
    d0 = 41,        ## days after the 1st day of week 26 to start
    ws = 40,        ## week of recording to start
    wf = 3,         ## week of recording to finish
    ## seedName = "CA Los Angeles County", ## up to here
    seedName = "First",
    seedSize = 50,
    dt = 1.0,
    ## ndt = 60,
    nreals = 1,
    verbose=FALSE,
    beta_fudge = 1.0,
    week1 = 26
) {
  
  ## Define the parmeters
  beta <- beta_fudge * R0 / Tg
  
  ## Derived parameters
  if (seedName == "First") {
    seedIndex <- 1
  } else {
    seedIndex <- match(seedName,df_pop$scname)
    if (is.na(seedIndex)) stop("name lookup failed for seed population")
  }
  
  ## Set the time scale for ctual days and reported weeks
  ifelse(ws < 53,ws <-ws-26+1,ws-ws+52-26+1) 
  ifelse(wf < 53,wf <-wf-26+1,wf- wf+52-26+1)
  nwk <- wf - ws + 1
  w0 <- d0 / 7
  ndt <- wf * 7 - d0
  ## TODO finish these edits here
  
  ## Setup vectors of S, I, R
  nPop <- dim(df_pop)[1]
  vecS <- vector(mode="numeric",length=nPop)
  vecI <- vector(mode="numeric",length=nPop)
  vecR <- vector(mode="numeric",length=nPop)
  vecN <- as.numeric(df_pop$size)
  
  ## Define the summary variables for each run
  runWhoWhom <- array(dim=c(nreals, nPop))
  runTimeSeed <- array(dim=c(nreals, nPop))
  runDistSeed <- array(dim=c(nreals, nPop))
  runInc <- matrix(nrow=nreals, ncol=ndt)
  
  ## Define data structures used for each simulation
  tsWhoInfectedWhom <- array(dim=c(ndt, nPop, nPop))
  tsFOI <- matrix(nrow = nPop, ncol = nPop)
  tsRec <- vector(mode="numeric", length = nPop)
  tsInf <- vector(mode="numeric", length = nPop)
  
  browser()
  ## TODO Finish defining the timescale
  
  weekoffset <- week1 + d0 %/% 7
  dayoffset <- d0 %% 7
  
  ## - number of infectious at each timestep in each place
  ## nwks <- ndt %/% 7 
  tsMatI <- matrix(nrow=ndt,ncol=nPop)
  raTsMatIInf <- array(dim=c(nreals,nPop,ndt))
  cumRaTsMatI <- array(dim=c(nreals,nPop,ndt))
  wksubinc <- array(dim=c(nreals,nPop,nwks))
  
  ## Calculate the connectivity matrix
  ## Loop should be made with apply funcs
  kappa <- matrix(nrow=nPop,ncol=nPop)
  matDist <- distm(cbind(df_pop$long,df_pop$lat))
  matDist <- matDist / 1000
  for (i in 1:nPop) {
    for (j in 1:nPop) {
      dist <- matDist[i,j]
      if (dist < 1) {rawkern <- 1}
      ## else {rawkern <- exp(-dist/kernPower)}
      else {rawkern <- (1/dist)^kernPower}
      kappa[i,j] <- rawkern
    }
    kappa[i,] <- kappa[i,] / sum(kappa[i,])
  }
  
  ## browser()
  ## image(log(kappa))
  
  ## Initiate each realization loop
  for (i in 1:nreals) {
    
    ## Up to here. Error below
    vecS[] <- vecN
    vecR[] <- 0
    vecI[] <- 0
    ## seedIndex <- match(seedName,df_pop$scname)
    vecS[seedIndex] <- vecS[seedIndex] - seedSize
    vecI[seedIndex] <- seedSize
    
    ## The timestep loop
    if (verbose) {cat("\n")}
    for (j in 1:ndt) {
      
      ## Need to give a bit of feedback here about progress
      if (verbose) {
        cat("\rTimestep ", j ," of ", ndt,". Realization ", i, " of ", nreals, ".")
      }
      
      ## Calculate the force of infection matrix
      ## Might be worth playing with different elements here
      ## This should be swapped out using native R matrix fuctions
      ## Force of infection by k from l
      ## Should be swapped out for c or while loops (maybe)
      for (k in 1:nPop) {
        for (l in 1:nPop) {
          tmpTsFOI <- 0
          for (m in 1:nPop) {
            destN <- vecN[m]
            if (destN > 0) {
              tmpTsFOI <- tmpTsFOI + beta * (kappa[k,m] * kappa[l,m] * vecI[l]) / vecN[m]
            } else {
              tmpTsFOI <- 0
            }
          }
          tsFOI[k,l] <- tmpTsFOI
        }
      }
      
      ## Total expected number of infection events for each
      ## susceptible population
      vecTotalHazInf <- rowSums(tsFOI)
      expInfEvents <- vecS*(1-exp(-dt*vecTotalHazInf))
      
      ## Calcuate the rate of recovery events
      probRecEvents <- 1-exp(-dt/Tg)
      
      ## Draw from multinomial for numebrs of infections and recoveries
      ## This could be replaced by a native R matrix function
      for (k in 1:nPop) { 
        if (expInfEvents[k]>0) { 
          tsWhoInfectedWhom[j,k,] <- rmultinom(1,expInfEvents[k],tsFOI[k,])
        } else {
          tsWhoInfectedWhom[j,k,] <- rep(0,nPop)
        }
        if (vecI[k] > 0) {
          tsRec[k] <- rbinom(1,vecI[k],probRecEvents)
        } else {
          tsRec[k] <- 0
        }
      }
      
      ## Update the state variables
      tsInf <- rowSums(tsWhoInfectedWhom[j,,])
      vecS <- vecS - tsInf
      vecI <- vecI + tsInf - tsRec
      vecR <- vecR + tsRec
      
      ## Record the number infectious at each timepoint
      ## Note that its the number at the end of the timestep
      ## to matchup with the main use of this structure
      tsMatI[j,] <- vecI
      raTsMatIInf[i,,j] <- tsInf    
      if (j==1) {
        cumRaTsMatI[i,,j] <- tsInf
      } else {
        cumRaTsMatI[i,,j] <- cumRaTsMatI[i,,j-1] + tsInf
      }
    }
    
    ## Record an overall epicurve for each run
    ## The obvious apply function needed a massive tmp vector
    ## still included as a comment here
    ## runInc[i,] <- apply(tsWhoInfectedWhom,1,sum)
    ## memory.size()
    this_ts <- 1
    while (this_ts <= ndt) {
      curtot <- 0
      this_pop_a <- 1
      while (this_pop_a <= nPop) {
        this_pop_b <- 1
        while (this_pop_b <= nPop) {
          curtot <- curtot + tsWhoInfectedWhom[this_ts,this_pop_a,this_pop_b]
          this_pop_b <- this_pop_b + 1
        }
        this_pop_a <- this_pop_a + 1
      }
      runInc[i,this_ts] <- curtot
      this_ts <- this_ts + 1
    }
    
    ## Identify seeds for each location
    vecRunLength <- vector(length=ndt,mode="numeric")
    for (k in 1:nPop) {
      
      ## First find the first timestep of the longest
      ## Note the reverse ordering of the loop
      ## Ensures that the max run count is the seed time
      maxRun <- 0
      for (j in ndt:1) {
        if (tsMatI[j,k] > 0) {
          maxRun <- maxRun + 1
        } else {
          maxRun <- 0
        }
        vecRunLength[j] <- maxRun
      }
      ## Not quite working
      intMaxLength <- max(vecRunLength)
      if (k==seedIndex) {
        runWhoWhom[i,k] <- -2
        runTimeSeed[i,k] <- -2
        runDistSeed[i,k] <- -2
      } else if (intMaxLength==0) {
        runWhoWhom[i,k] <- -1
        runTimeSeed[i,k] <- -1
        runDistSeed[i,k] <- -1
      } else {
        timeSeed <- match(intMaxLength,vecRunLength)
        runTimeSeed[i,k] <- timeSeed
        vecPossSeeds <- tsWhoInfectedWhom[timeSeed,k,]
        vecPossSeeds[k] <- 0
        infector <- match(1,rmultinom(1,1,vecPossSeeds))
        runWhoWhom[i,k] <- infector
        runDistSeed[i,k] <- matDist[k,infector]
      }
    }
  }

  if (verbose) {cat("\n")}

  browser()
  
  ## Assign a day 0
  
  ## Switch this into weekly via a cumulative version if needed
  wksubinc[,,] <- cumRaTsMatI[,,(1:nwks)*7]
  wksubinc[,,2:nwks] <- wksubinc[,,2:nwks] - wksubinc[,,1:(nwks-1)] 
  wksubinc <- wksubinc * p_c + trickle*7
  
  ## Wrangle some of the required objects to be returned
  list(
    inc = runInc,
    subinc=raTsMatIInf, 
    cumsubinc=cumRaTsMatI,
    wksubinc=wksubinc
  )
  
}

#' Function to give a rough overall log likelihood for likleihood for 
#' To be efficient, needs to assume that all counties in the data are
#' in the model and that all season weeks in the data are in the model.
#' Assumes strictly that the ordering of the state county names is the 
#' same for the 
county_lnlike <-function(
    v_theta,              ## Vector of names parameters used in the function
    mat_mod,              ## Matrix of model results 1st D no pops, 2nd D no weeks
    df_dat,               ## Database of observations
    df_pop,               ## Population database with list of state county names
    seas=2017,            ## Season that we want to fit
    result = "num_spec",  ## The result field from df_dat that we are fitting
    pseudo = FALSE        ## Indicator if seudo data should be created 
) {
  
  ## Add a few checks that all the observations are in the model results
  ## and also that the datasets are correctly ordered
    
  ## set return objects
  noobs <- dim(df_dat)[1]
  lnlike <- 0
  obs_inc <- 0
  obs_mat <- matrix(nrow=nrow(mat_mod),ncol=ncol(mat_mod))
  if (pseudo) {
    df_dat <- cbind(df_dat,pseudo=rep(-999,noobs))
  } 
  
  ## setup loop limits
  noweeks <- dim(mat_mod)[2]
  weeks <- (1:noweeks) + v_theta["w_1"]
  nodat <- dim(df_dat)[1]
  nopops <- dim(df_pop)[1]
  
  dfscname <- as.character(df_pop$scname)
  vecdatweek <- df_dat$seasweek - v_theta["w_1"] + 1
  
  ## Main loop to generate a likelihood
  i <- 1
  mod_sc_ind <- 1  
  while (i <= nodat) {
    dat_scn <- df_dat[i,"scname"]
    dat_seas <- df_dat[i,"season"]
    dat_week <- vecdatweek[i]
  if (dat_seas == seas) {
    while (dfscname[mod_sc_ind] != dat_scn) { 
      mod_sc_ind <- mod_sc_ind + 1
      if (mod_sc_ind > nopops) error("While loop overrun wnpzrd")
    }
    if ((dat_week >= 1) && (dat_week <= noweeks)) {
      obs_inf <- df_dat[i,result]
      if (!is.na(obs_inf)) {
        if (obs_inf < 0) {
          stop("Non NA value sbelow zero are probably bad pseudo data") 
        }
      }
      obs_mat[mod_sc_ind,dat_week] <- obs_inf
      exp_inf <- mat_mod[mod_sc_ind,dat_week] 
      lnlike <- lnlike + dpois(obs_inf,exp_inf,log=TRUE)
      df_dat[i,"pseudo"] <- round(exp_inf)
      obs_inc <- obs_inc + 1
    }
  }
  i <- i + 1 
  }
  
  if (pseudo) {
    datrtn <- df_dat
  } else {
    datrtn <- NULL
  }
  
  list(lnlike=lnlike,obsinc=obs_inc,obsmat=obs_mat,pseudodat = datrtn)
  
}

mask_counties <- function(
    dfdat,
    dfpop,
    lat=39.83,
    long=-98.58,
    dist=200) {
  
  require(geosphere)
    
  vecdist <- geosphere::distHaversine(
      c(long,lat),
      as.matrix(cbind(dfpop$long,dfpop$lat),
      r=6378137)
  )
  
  mask <- (vecdist < (dist * 1000))
  
  dfpop_new <- dfpop[mask,]
  dfdat_new <- dfdat[dfdat$scname %in% dfpop_new$scname,]  
  
  list(dat=dfdat_new, pop=dfpop_new)
  
}

load_quidel_proto_data <- function(
    dataDir = "~/Dropbox/shares/pete_LEPR03/data/quidel/",
    week1=26) {
  
  #' So had a look at the data that Jamie extracted
  
  
#' Construct the state county combo
  dat_tmp1 <- read.csv(paste(dataDir,"county_weekly.csv",sep=""))
  dat_tmp1$scname <- paste(dat_tmp1$state,dat_tmp1$county)
  head(dat_tmp1$scname)
  
#' Have a quick look at the data 
  head(dat_tmp1)
  tail(sort(table(dat_tmp1$county)),n=20)
  
#' Load up the US county data and create the right table
  dfStateLU <- read.table(paste(dataDir,"../population/state_lookups_v2.csv",sep=""),
      sep="\t",row.names=NULL,header=TRUE)
  dfUSGaz <- read.table(
      paste(dataDir,"../population/2015_Gaz_counties_national.txt",sep=""),
      header=TRUE,
      sep="\t",
      colClasses=c("NULL","numeric","NULL","character",
          "NULL","NULL","NULL","NULL","numeric","numeric")
  )
  dfPopEst <- read.csv(paste(dataDir,"../population/co-est2016-alldata.csv",sep=""))
  dfPopEst <- dfPopEst[,c("STATE","COUNTY","CTYNAME","CENSUS2010POP")]
  dfPopEst$GEOID <- paste(
      dfPopEst$STATE,sprintf("%03d",dfPopEst$COUNTY),sep="")
  dfUSRaw <- merge(dfPopEst,dfUSGaz,by=c("GEOID"),all.y=TRUE)
  sum(dfUSRaw$CTYNAME[1:2720]==dfUSRaw$NAME[1:2720])
  sum(dfUSRaw$CTYNAME[1:2721]==dfUSRaw$NAME[1:2721])
  dfUSRaw <- dfUSRaw[1:2720,]
  sum(dfUSRaw$CTYNAME==dfUSRaw$NAME)
  dfUSRaw <- dfUSRaw[dfUSRaw$CTYNAME==dfUSRaw$NAME,]
  sum(dfUSRaw$CTYNAME==dfUSRaw$NAME)
  dim(dfUSRaw)
  dfPop = data.frame(
      scname=paste(
          dfStateLU$usps[match(as.numeric(dfUSRaw$STATE),dfStateLU$fips)],
          dfUSRaw$NAME),
      ## name=dfUSRaw$NAME,
      lat=dfUSRaw$INTPTLAT,
      long=dfUSRaw$INTPTLONG,
      size=dfUSRaw$CENSUS2010POP
  )
  dfPop <- dfPop[order(dfPop$scname),]
  head(dfPop)
  
#' Check for counties in the data that are not in the county data table
#' and then extract them from the quidel data so we have at least a consistent 
#' pair of datasets. Need to look at the bad names to see if they are obviously 
#' in Alaska or something like that
  sum(dat_tmp1$scname %in% dfPop$scname) - length(dat_tmp1$scname)
  badnames <- unique(dat_tmp1[!(dat_tmp1$scname %in% dfPop$scname),"county"])
  length(badnames)
  
#' Reduce the quidel data to only those with county names that are in our list of counties
  dat <- dat_tmp1[!(dat_tmp1$scname %in% badnames),]
  sum(dat$scname %in% dfPop$scname) - length(dat$scname)
  dim(dat)
  names(dat)
  head(dat)
  table(dat$num_sofias)
  
  ##  Define seasons and season weeks for the data
  ## no week 53s for some reason
  dat$seasweek <- ifelse(
      dat$week >= week1, 
      dat$week - week1 + 1, 
      dat$week + (52 - week1) +1)
  dat$season <- ifelse(
      dat$week >= week1, 
      dat$year, 
      dat$year-1)
  
#' Merge on the population sizes from the county database and sort the data
  dat <- merge(dat,data.frame(scname=dfPop$scname,popsize=dfPop$size))
  dat <- dat[order(dat$scname,dat$season,dat$seasweek),]
  
  list(dat=dat,dfPop=dfPop,week1=week1)
  
}
