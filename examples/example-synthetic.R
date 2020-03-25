##
## Driver for modeling  of Synthetic Data - 
## The script first retrives all the data for a 'cdc' type run modeling the 
## national data using the 10 HHS regions (either coupled or not).  It then
## replaces the incidence (%ILI) and number of cases with simple synthetic profiles
## which it then fits.  This example has been set up with:
## data_source = 'cdc', RegState = 'usa', mod_level = 2, fit_level = 3, and model = 4
## It does not really matter what season is chosen
## 

rm(list=ls())

## Set default values for all parameters - these are being reset when the script is called
## with difefernt values

# data type to be modeled - cdc or gft

data_source = 'cdc'


## Start year of the flu season
year  = 2015

## Spatial level of model/fit data - 2 entire usa, 3 HHS regions, 4 States
## For CDC we only support level 2  
## For GFT we support levels 2 (entire usa) and  3 (a single HHS region)
mod_level = 2

## To fit HHS regions using state level data use: mod_level=2 and fit_level = 3
## To fit a single HHS region using state level data use mod_level = 3 and fit_level = 4
## In this case you will also need to say which of the 1-10 HHS regions you would like to fit.
## This is done by setting the parametr RegState below to a number between 1-10


## Spatial level of data used to fit the model data >= mod_level
##

fit_level = 3

## 
## The Name of the Model Data
## If mod_level = 2 this should be set to 'usa'
## RegState - if mod_level = 3 choose the HHS region number, 1-10.
## 

RegState = 'usa'

##
## Number of MCMC chains
nreal = 1


## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model
model  = 4


## -------------------------------------------
## Parameters that the User may want to change
##
## Number of steps/trials in each MCMC chain
nMCMC = 1e6

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions
isingle     = 1

##
## Name of sub-directory where all the output files/plots will be saved. 

subDir  <- "output" 

##
## plot - 0, 1 or 2 (no, basic plotting or ggplot2)
plot = 1

## 
## device - device type for plots: currently supports pdf, png and x11 and can accpet all in a vector: device = c('png','pdf')

device = 'pdf'

## 
## Seed for RNG if NULL DICE will seed it

iseed = NULL


## End of parameters that the User may want to change
## ----------------------------------------------------


## load DICE
## to see what data has been loaded with DICE use: data(package='dice')

require(DICE)

## Retrive data using the User's choices

mydata <- get.DICE.data(data_source = data_source, mod_level = mod_level, fit_level = fit_level, year = year, model = model, RegState = RegState, isingle=isingle)


## augment the mydata with some more information that the user has chosen
Temp = 1

prior = 0

mydata$single = isingle

mydata$prior = prior

mydata$Temp = Temp

##
## Pack the information for the 'synthetic' run 

opt.list <- set.opt.list(model = mydata$imodel, isingle = isingle)

misc.list <- set.misc.list(model = mydata$imodel)

run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, device = device, subDir = subDir, plot = plot)

##
## Fit the synthetic data and generate all the plots/files.

if (isingle == 1) {
	output = runSyntheticSingle(mydata = mydata, opt.list = opt.list, run.list = run.list, misc.list = misc.list, ireal = 1, iseed = iseed)	
} else {
	output = runSyntheticMulti( mydata = mydata, opt.list = opt.list, run.list = run.list, misc.list = misc.list, ireal = 1, iseed = iseed)	
}


