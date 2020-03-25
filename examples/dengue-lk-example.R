##
## Dengue Example for a DICE run - Using Rscript or R CMD BATCH from the Command Line
##

rm(list=ls())

library('batch')

## Set default values for all parameters - these are being reset when the script is called
## with differnt values


##
## Fitting/Forecasting Method either compartmental models ('mech') or statistical ('stat')
##
method = 'mech' 


## 
## Define the SARIMA model if method = 'stat' Default is NULL and code uses auto.arima
##
arima_model = list(p = 1, d = 0, q = 0, P = 3, D = 1, Q = 0) # NULL

##
## Optional co-variate for SARIMA. Used ONLY if 'arima_model' is defined by the User and method = 'stat'
##

covar = 'precip' # or 'temp' or 'precip' default is FALSE

##
## Lag time (in cadence of data units) for co-variate
## 
covar_lag = 1 

##
## Disease to model 
##

disease ='dengue' 

##
## An integer or string-abbreviation selecting which data source should be used.  A list of data sources appears in DICE_data_sources.  If value is NULL, ## DICE will attempt to auto-choose a data source.
##

data_source="Sri_weekly"
## 
## Use the MySQl database (dengue) or not (flu)
##

db_opts=list(DICE_db="predsci", CDC_server=FALSE)

##
## Define mod_name
##

RegState = NULL

mod_name = c(NAME_2="LK")

fit_names = 'all'

##
## Start year of the disease season
##

year  = 2017

##
## Spatial level of model/fit data 
##

mod_level = 2

##
## Spatial level of data used to fit the model data >= mod_level
##

fit_level = 2

##
## Number of data points that are fitted - default is to fit all of the available data
##
nfit = 52

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions
isingle     = 1

## 
## SIR (1) SEIR (2) can also accept string: sir, seir, vsir , vseir or: SIR, SEIR, VSIR, VSEIR
## 
epi_model = 2

##
## Generation time  (in days)
##
Tg = 10.0

## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model
## Models 2 and 3 should not used for dengue
## Model 1 also makes little sense as does model 5..
## The code will reset the model to 4 if the user chose 2 or 3 
##
model  = 4

## Select Method for Data Augmentation
## 0 No data augmentation
## 1 Use historic monthly average, NULL model
## 2 Use the most similar Season
da = 0
##
## Use a prior (prior = 1) or no (prior = 0). Currently relevant only for flu
##
prior = 0
##
## Temperature for LLK in MCMC procedure
##
Temp = 10

##
## Number of steps/trials in MCMC chain
##
nMCMC = 2e6

##
## Number of MCMC chains
##
nreal = 1

plot = TRUE

##
## Optional- Name of sub-directory where all the output files/plots will be saved. If set to NULL code will build the directory
## name based on the parameters selection
##
subDir = NULL 

## args will be a data frame
args = parseCommandArgs()

if(length(args) == 0) {
	cat("
	RegState      name of region/state we want to model (case insensitive BR, MX, TH)
	method        String mech or stat
	DICE_db       Determine which SQL database to use 'PredSci' or 'BSVE' 
  CDC_server    Logical determining if CDC data should be retrieved from the CDC server
	disease       flu or dengue
	epi_model     SIR or SEIR (can also use 1 or 2 or lower case)
	arima_model   If method = 'stat' a list with p, d, q, P, D, Q values (Default NULL)
	covar         Optional co-variate to use in SARIMA modeling.  'sh', 'temp', 'precip'
	covar_lag     optional time-lag for covariate variable. Time units are the same as cadence of data
	year          starting year for disease season
	nfit          Number of data points to fit.  Default is to fit all available data for the season	
	mod_level     Spatial level of model data 
	fit_level     Spatial level of data used in fitting of mod_level.  fit_level >= mod_level	
	isingle       If set to 1==TRUE each region's profile is fitted individually. If isingle = 0 the fit is coupled	
	model         Model for the force of infection. 1-5 for flu 4-5 for dengue
	Tg            Recovery time in days, relevant for SIR and SEIR.  Default is 3 and 10 days for flu/dengue
	da            Data augmentation option. default is da = 0, no augmentation. da=1/2 uses the average/most similar season
	prior         Only for flu. Use a prior (prior=1) or not(prior = 0). Default is prior = 0
	Temp          Temperature for LLK in MCMC procedure. Default is Temp = 1. 
	nMCMC         Number of MCMC steps in chain
	nreal         Number of MCMC chains/realizations. Script runs one at a time \n
	Example: \n
	Rscript dengue-example.R method mech epi_model SEIR RegState MX mod_level 2 fit_level 3 year 2010 Tg 8  isingle 0 nMCMC 1e6
	Rscript dengue-example.R method mech epi_model SIR  RegState BR mod_level 2 fit_level 3 year 2010 Tg 10 isingle 1 nMCMC 1e6 
	Rscript dengue-example.R method stat RegState TH year 2012 mod_level 2 fit_level 3 `list(p = 1, d = 0, q = 0, P = 3, D = 1, Q = 0)` 
	Rscript dengue-example.R method stat RegState MX year 2012 mod_level 2 fit_level 3 `list(p = 1, d = 0, q = 0, P = 3, D = 1, Q = 0)` covar precip covar_lag 1
	
	\n\n
	")
	q(save='no')
} else {
	if ('RegState'    %in% names(args))  RegState = args$RegState
	if ('method'      %in% names(args))  method = args$method
	if ('disease'     %in% names(args))  disease = args$disease
	if ('epi_model'   %in% names(args))  epi_model = args$epi_model
	if ('arima_model' %in% names(args))  arima_model = args$arima_model
	if ('covar'       %in% names(args))  covar = args$covar
	if ('covar_lag'   %in% names(args))  covar_lag = args$covar_lag
	if ('data_source'    %in% names(args))  data_source = args$data_source
	if ('year'   	  %in% names(args))  year = args$year	
	if ('nreal'       %in% names(args))  nreal = args$nreal
	if ('nMCMC'       %in% names(args))  nMCMC = args$nMCMC
	if ('model'  	  %in% names(args))  model= args$model
	if ('nfit'        %in% names(args))  nfit = args$nfit
	if ('isingle' 	  %in% names(args))  isingle = args$isingle
	if ('Tg'          %in% names(args))  Tg = args$Tg
	if ('prior'       %in% names(args))  prior = args$prior
	if ('Temp '       %in% names(args))  Temp  = args$Temp	
	if ('da'          %in% names(args))  da    = args$da
	if ('DICE_db'     %in% names(args))  db_opts$DICE_db = args$DICE_db
	if ('CDC_server'  %in% names(args))  db_opts$CDC_server = args$CDC_server
 }


## load vectorDICE
## to see what data has been loaded with DICE use: data(package='vectorDICE')

require(DICE)

## start the clock
start.time <- proc.time()

## Start a DICE Run
##


output <- runDICE(data_source=data_source,year=year, nreal=nreal, nMCMC = nMCMC, mod_level = mod_level, fit_level = fit_level, RegState = RegState,mod_name = mod_name, model= model, nfit = nfit, isingle = isingle, subDir = subDir, device = 'pdf', plot = plot, Tg = Tg, prior = prior, Temp = Temp, da = da, epi_model=epi_model, db_opts=db_opts, disease = disease, method = method, arima_model = arima_model, covar = covar, covar_lag = covar_lag)



# stop the clock
cat("\n\nElapsed Time: \n\n")
end.time = proc.time()
run.time = end.time-start.time
print(end.time - start.time)



