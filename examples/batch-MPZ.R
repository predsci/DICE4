##
## DRIVER FOR base level runs for military data. Currently for internal use only
##



rm(list=ls())

require(DICE)

##
## Number of MCMC chains
nreal = 1

##
## Number of steps/trials in each MCMC chain
nMCMC = 1e5

## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model

model  = 5

##
## Number of weeks of data that are fitted - default is to fit all of the available data - for the 2009-2010 pandemic it is 66 weeks long!
## if nperiodsFit is set to 0, the code will reset it to nperiods 

nperiodsFit = 52

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions, here it must be 1 since this is a single patch
isingle     = 1


## 
## device - device type for plots: currently supports pdf, png and x11 and can accpet all in a vector

device = 'pdf'

## Set default values for all parameters - these are being reset when the script is called
## with difefernt values

# data type to be modeled - cdc or gft

data_source = 'misccases'

##
## Plotting - 0 = NO, 1 = Base Plot, 2 = ggplot2 

plot = 2


## Start year of the flu season
year  = 2009

start_week = 13

## Choose a zipcode and then let DICE determine the Country, state and if relevant the HHS region this MPZ is in:
##

MPZind = diceData$attr$ABBV_5<"AAAAA" & diceData$attr$ABBV_5 != ""
MPZattr = diceData$attr[MPZind,]

## To see a vector of base names use: MPZattr$NAME_5
##

## To see a vector of zip codes use: as.numeric(MPZattr$ABBV_5)
##

## If you know the zip code just set it, for example
## Portsmouth, VA 23708
NAME_5 = '23708'

index = which(MPZattr$ABBV_5 == as.numeric(NAME_5))

NAME_4 = MPZattr$ABBV_4[index]

NAME_3 = MPZattr$NAME_3[index]

NAME_2 = MPZattr$ABBV_2[index]

# ## Fort Carson, CO 80913 set:

# NAME_5 = '80913'

library('batch')

## args will be a data frame 
args = parseCommandArgs()

if(length(args) == 0) {
	cat("
	data_source      type of data: cdc or gft    
	year      starting year for ILI season
	nMCMC     Number of MCMC steps in chain
	nreal     Number of MCMC chains/realizations. Script runs one at a time
	nperiodsFit Number of weeks to fit.  Default is the entire season which is 52 weeks
	isingle    If set to 1==TRUE each region's profile is fitted individually \n
	Example: \n
	Rscript batch-example.R data_source cdc year 2015 real 1 isingle 0  nMCMC 1e5 nperiods 45  
	\n\n
	")
	q(save='no')
} else {
	if ('NAME_3'  %in% names(args))  NAME_3 = args$NAME_3
	if ('NAME_4'  %in% names(args))  NAME_4 = args$NAME_4	
	if ('NAME_5'  %in% names(args))  NAME_5 = args$NAME_5	
	if ('year'   	%in% names(args))  year = args$year	
	if ('nreal'     %in% names(args))  nreal = args$nreal
	if ('nMCMC'     %in% names(args))  nMCMC = args$nMCMC
	if ('model'  	%in% names(args))  model= args$model
	if ('nperiodsFit' %in% names(args))  nperiodsFit = args$nperiodsFit
	if ('subDir'    %in% names(args))  subDir = args$subDir
	if ('device'    %in% names(args))  device = args$device	
 }
 
mydata = get.subset(start.year=year,end.year=(year+1),start.week=start_week,mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)
 

## start the clock 
start.time <- proc.time()

##
## Name of sub-directory where all the output files/plots will be saved. 

subDir  <- paste(mydata$model$attr$ABBV_4,'-data-',year,'-',(year+1),'-v',model,sep="")


output <- runSinglePatch(mydata=mydata, year=year, nreal=nreal, nMCMC = nMCMC, model= model, nperiodsFit=nperiodsFit, isingle = isingle, subDir = subDir, device = device, plot = plot)
 

# stop the clock        
cat("\n\nElapsed Time: \n\n")
end.time = proc.time()
run.time = end.time-start.time
print(end.time - start.time)


