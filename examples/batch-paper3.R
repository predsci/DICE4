##
## DRIVER FOR county level runs for paper-3. Currently for internal use only
##



rm(list=ls())


require(DICE)

##
## Number of MCMC chains
nreal = 1

##
## Number of steps/trials in each MCMC chain
nMCMC = 1e4

## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model
model  = 5

##
## Number of weeks of data that are fitted - default is to fit all of the available data
nperiodsFit = 52

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions
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

plot = 1


##
## 
prior = 0
Temp  = 1

## Start year of the flu season
year  = 2010

NAME_2 = 'USA'

## To retrieve the data for any of the five counties use:
##


## AZ - MARICOPA

NAME_3 = "Region9"

NAME_4 = 'Arizona'

NAME_5 = 'Maricopa'

# mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)

# ## CA _ SAN DIEGO

# NAME_3 = "Region9"

# NAME_4 = "California"

# NAME_5 = "San.Diego"

# mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)

# ## MO - EASTERN

# NAME_3 = "Region7"

# NAME_4 = "Missouri"

# NAME_5 = 'Eastern'

# mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)

# ## VA - EASTERN

# NAME_3 = "Region3"

# NAME_4 = "Virginia"

# NAME_5 = "Eastern"

# mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)

# ## TN - NDR

# NAME_3 = 'Region4'

# NAME_4 = "Tennessee"

# NAME_5 = "Nashville-Davidson"

# mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)


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
	if ('prior'     %in% names(args))  prior = args$prior
	if ('Temp '     %in% names(args))  Temp  = args$Temp	
	if ('RegState'  %in% names(args))  RegState = args$RegState		
 }
 
mydata = get.subset(start.year=year,end.year=(year+1),mod_level = 5,fit_level = 5,name=c(NAME_2=NAME_2,NAME_3=NAME_3,NAME_4=NAME_4,NAME_5=NAME_5),data_source=data_source)
 

## start the clock 
start.time <- proc.time()

##
## Name of sub-directory where all the output files/plots will be saved. 

subDir  <- paste(mydata$model$attr$ABBV_4,'-data-',year,'-',(year+1),'-v',model,sep="")


output <- runSinglePatch(mydata=mydata, year=year, nreal=nreal, nMCMC = nMCMC, model= model, nperiodsFit=nperiodsFit, isingle = isingle, subDir = subDir, device = device, plot = plot, Temp = Temp, prior = prior)
 

# stop the clock        
cat("\n\nElapsed Time: \n\n")
end.time = proc.time()
run.time = end.time-start.time
print(end.time - start.time)


