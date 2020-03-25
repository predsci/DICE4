#' # Steven's test script to play a little with the quidel data
#' 
#' Designed to make a report by running spin("this_file_name") from package knitr.
#'  
#' ## Outline ideas
#' 
#' The general flow of a paper could be to introduce and describe the data and then
#' make epidemic models out of it with the best possible forecast models that uswe only
#' regional averages of the previous year. There could be a step mapping those models
#' onto the current year.
#' 
#' ## Data wrangling
#' 
#' irst set the current wd if not using this in batch mode and then remove
#' all current objects. ANd make sure we can use quite a bit of memory
## setwd("~/Dropbox/git/DICE-private/examples")
memory.limit(size=10000)
rm(list = ls(all=TRUE))

#' Now require DICE and check that it loads with all its dependencies
## library(DICE)
source("~/Dropbox/git/DICE-private/dice/R/county_proto.R")
library(geosphere)

tmp <- load_quidel_proto_data()

dat <- tmp$dat
dfPop <- tmp$dfPop
season_start <- tmp$week1


## Should have a function that ends here by producing the incidence data and pop data


#' Reduce the problem set to a smaller number of counties if needed
tmp_2 <- mask_counties(dat,dfPop,lat=33.75,long=-84.39,dist=50)

#' Need to sort variable naming below here
dat_main <- tmp_2$dat
dfPop_main <- tmp_2$pop
dim(dfPop_main)


#' ## Get some results from a prototype model
#' 
#' First setup some pseudo data
x1 <- county_proto(
    dfPop_main,
    nreals = 1,
    verbose=TRUE,
    kernPower = 0.4,
    beta_fudge = 0.6,
    p_c = 0.001,
    trickle = 1/3000000)

x1$wksubinc[1,,]

theta_base <- c(
    w_1 = 40-season_start
)

tmp <- county_lnlike(
    theta_base,
    as.matrix(x1$wksubinc[1,,]),
    dat_main,
    dfPop_main,
    seas=2018,
    pseudo=TRUE)
dat_pseudo <- tmp$pseudo

vecp1 <- seq(0,3,0.2) 
vecp2 <- 10^(seq(-6,0,0.25))
nop1 <- length(vecp1)
nop2 <- length(vecp2)
nosamples <- nop1*nop2
scanres <- matrix(nrow=length(vecp1),ncol=length(vecp2))
rownames(scanres) <- vecp1
colnames(scanres) <- vecp2

cat("\n")
for (p1 in vecp1) {
  for (p2 in vecp2) {
    x2 <- county_proto(
        dfPop_main,
        nreals = 1,
        verbose=FALSE,
        kernPower = p1,
        p_c = p2,
        beta_fudge = 0.6)
    wksubinc <- x2$wksubinc[1,,]
    
    #' Some a paramter vector that we will use for the likelihood
    theta <- c(
        w_1 = 40-season_start
    )
    
    #' Get a single value of the likelihood and the number of observations used
    modmat <- as.matrix(wksubinc)
    tmp <- county_lnlike(
        theta,
        modmat,
        dat_pseudo,
        dfPop_main,
        seas=2018,
        result="pseudo")
        ## result="num_spec")
        countp1 <- match(p1,vecp1)
    countp2 <- match(p2,vecp2)
    scanres[countp1,countp2] <- tmp$lnlike
    
    cat("\rParam sample ", countp2 + (countp1-1)*nop2 ," of ", nosamples,".")
    
  }
}

#' So this seems to work pretty well. Lets try it at a slightly larger size
image(log10(abs(scanres)),x=vecp1,y=log10(vecp2))

maskmodel <- dat_pseudo$pseudo > -0.1
hist(dat_pseudo$num_spec[maskmodel])
hist(dat_pseudo$pseudo[maskmodel])
hist(dat_pseudo$num_spec[maskmodel] - dat_pseudo$pseudo[maskmodel])

#' Plot variation by sub population 
noreals <- dim(x2$inc)[1]
ntp <- dim(x2$inc)[2]
plot(x2$inc[1,],type="l")
if (noreals > 1) {
  for (i in 2:noreals) {
    points(x2$inc[i,],type="l")
  }
}

#' Plot the county specific incidences and they look too small
#' large at the moment
nosubs <- dim(x2$subinc)[2]
chosen_inc <- x2$subinc[1,,]
plot(chosen_inc[1,],type="l",ylim=c(0,max(x2$subinc)))
for (i in 2:nosubs){
  points(chosen_inc[i,],type="l")
}
hist(rowSums(chosen_inc) / dfPop_main$size)

image(wksubinc)
head(wksubinc)
sum(wksubinc) / sum(dfPop_main$size)