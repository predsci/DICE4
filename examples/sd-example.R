##
## Driver for modeling San Diego County Data
## 

rm(list=ls())

library('batch')

##
## Set default values for all parameters - these can be reset when the script is called
## with difefernt values

##
## Fitting/Forecasting Method either compartmental models ('mech') or statistical ('stat')
##
method = 'mech' 

## 
## Define the SARIMA model if method = 'stat' Default is NULL and code uses auto.arima
##
arima_model = list(p = 1, d = 1, q = 0, P = 0, D = 0, Q = 0) # or NULL

##
## Optional co-variate for SARIMA. Used ONLY if 'arima_model' is defined by the User and method = 'stat'
##

covar = 'sh' # or 'temp' or 'precip' default is FALSE

##
## Lag time (in cadence of data units) for co-variate
## 
covar_lag = 1 

##
## Disease to model 
##

disease ='flu' 

## 
## Use the MySQl database 
##
db_opts=list(DICE_db="predsci", CDC_server=FALSE)

##
## data source
##
data_source = NULL

##
## Start year of the disease season
##

year  = 2018

## 
## Define the name - complicated in this case
## 

mod_name=c(NAME_2="US", NAME_3="R9", NAME_4="CA", NAME_5="CHD1", NAME_6="SD")

##
## Spatial level of model data - 6 is county level 
##

mod_level = fit_level = 6

##
## Number of data points that are fitted - default is to fit all of the available data
##
nfit = 52

##
## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the fit_level spatial regions, relevant only when method='mech'
##
isingle = 1

## 
## SIR (1) SEIR (2) can also accept string: sir, seir or SIR, SEIR
## 

epi_model = 1

##
## Generation time
##
Tg = 3.0

## 
## Model Number (1-5), 1- SH Only, 2 - School Vacation only, 3 - Both SH and SV, 4 - Constant, 5 -  A two value model
##
model  = 4

##
## Number of steps/trials in each MCMC chain
##

nMCMC = 1e6

##
## Number of MCMC chains
##
nreal = 1
##
## Use a prior (prior = 1) or no (prior = 0). Currently relevant only for flu
##
prior = 0

##
## Temperature for LLK (1, 10, 100)
##
Temp = 1

## 
## Data Augmentation - 0 (No), 1 (Historic Null Model), 2 (Most Similar Season)
##

da = 0

##
## plot - TRUE, FALSE or EXTERNAL (or 0, 1 and 2 respectively) No reason to ever set to FALSE
##
plot = TRUE

## 
## Optional - device - device type for plots: currently supports pdf, png and x11 and can accpet all in a vector: device = c('png','pdf')
##
device = 'pdf'

##
## Optional- Name of sub-directory where all the output files/plots will be saved. If set to NULL code will build the directory
## name based on the parameters selection
##
subDir = NULL 

##
## args will be a data frame 
##
args = parseCommandArgs()

if(length(args) == 0) {
	cat("
	mod_name      name of region/state/county we want to model (case insensitive)
	method        String mech or stat
  	DICE_db       Determine which SQL database to use 'PredSci' or 'BSVE' 
  	CDC_server    Logical determining if CDC data should be retrieved from the CDC server
	disease       flu or dengue
	epi_model     SIR or SEIR (can also use 1 or 2 or lower case)
	arima_model   If method = 'stat' a list with p, d, q, P, D, Q values (Default NULL)
	covar         Optional co-variate to use in SARIMA modeling.  'sh', 'temp', 'precip'
	covar_lag     optional time-lag for covariate variable. Time units are the same as cadence of data
	data_source      type of data: cdc, gft, dengue    
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
	Examples: \n
	Rscript batch-example.R disease flu method mech epi_model SIR data_source cdc year 2016 model 3 nreal 1 isingle 0  nMCMC 1e5 nfit 45 \n
	Rscript batch-example.R disease flu method stat 
	
	\n\n
	")
	q(save='no')
} else {
	if ('mod_name'    %in% names(args))  mod_name = args$mod_name
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
	if ('subDir'      %in% names(args))  subDir = args$subDir
	if ('device'      %in% names(args))  device = args$device
	if ('plot'        %in% names(args))  plot = args$plot	
	if ('Tg'          %in% names(args))  Tg = args$Tg
	if ('prior'       %in% names(args))  prior = args$prior
	if ('Temp '       %in% names(args))  Temp  = args$Temp	
	if ('da'          %in% names(args))  da    = args$da
	if ('db_opts'     %in% names(args))  db_opts = args$db_opts
 }

##
## load DICE
## to see what data has been loaded with DICE use: data(package='dice')
##

require(DICE)

# ## First let's just get the data 

# all_years_data = get.DICE.data(data_source=data_source, mod_level=mod_level, fit_level=fit_level, mod_name=mod_name, disease=disease, year=year, db_opts=db_opts)

# mydata = all_years_data$mydata
# all_years_epi = all_years_data$all_years_epi
# all_years_clim = all_years_data$all_years_clim

## start the clock 
start.time <- proc.time()

output <- runDICE(data_source=data_source,year=year, nreal=nreal, nMCMC = nMCMC, mod_level = mod_level, fit_level = fit_level, mod_name = mod_name,model= model, nfit = nfit, isingle = isingle, subDir = subDir, device = device, plot = plot, Tg = Tg, prior = prior, Temp = Temp, da = da, epi_model=epi_model, db_opts=db_opts, disease = disease, method = method, arima_model = arima_model, covar = covar, covar_lag = covar_lag)

# stop the clock        
cat("\n\nElapsed Time: \n\n")
end.time = proc.time()
run.time = end.time-start.time
print(end.time - start.time)





 

