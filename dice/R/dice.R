##
## Main DICE driver functions
##

runDICE <- function(data_source = NULL, year = 2016, mod_level = 2, fit_level = 3, nfit = 52, model = 5, isingle = 0, nMCMC = 1e+05, nreal = 1, device = "pdf", prior = 0, Temp = 1, da = 0, mod_name=c(NAME_2="US"), RegState=NULL, fit_names="all", subDir = NULL, plot = 1, iseed = NULL, Tg = NULL, epi_model = NULL, disease='flu', db_opts=list(DICE_db="predsci", CDC_server=TRUE), raw_col=NULL, arima_model = NULL, method = 'mech', covar = FALSE, covar_lag = 1, all_data=NULL) {

  #' Main Driver for the \code{DICE} package
  #'
  #' The main driver for \code{DICE} which grabs the data for the requested dataset, season and model/fit spatial regions combinations.
  #' After data is retrieved, the simulation is setup and the function calls either the single or multi options - depending on the user's
  #' request for an uncoupled or coupled run. In either case the fit begins with an MCMC procedure on the model-level data.
  #' In the case of statistical modeling - the code fits the model level data both directly and as a weighted sum of the region level
  #' data fits.
  #' @param db_opts A list of database options.  $DICE_db Determines which SQL database the data is retrieved from.  'PredSci' is the default SQL database, 'BSVE' is in development.  Additional flags are for outside sources of data (currently only the CDC Influenza-Like_Illness (ILI) is supported: $CDC_server=TRUE).
  #' @param disease String - disease name. Options for modeling are: flu, dengue, yellow$\_$fever, ebola, zika, cholera, chik, plague. To graphically explore the data see: \url{predsci.com/id$\_$data/}. A full list of diseases in the DICE database can be found from an R-prompt by following one of the examples below.
  #' @param data_source Describes the data source for the incidence data. Default is 'cdc' (for \code{disease = 'flu'}). It can be selected by source_key (integer) or source abbreviation (string). Most disease/location combinations have only one data source.  In this case, it may be easier to set data_source=NULL.  However, when multiple data sources exist, setting data_source=NULL will essentially choose from the available sources at random.  To determine a data source by graphical interface, see: \url{predsci.com/id_data/}.  Looking-up the disease and location will result in a list of data sources that can be entered into DICE.  Alternatively, all country/disease/data_source combinations are listed in the `Data Sources Table' tab at the same url.  To access the list of sources directly from an R-prompt, see the examples below.
  #' @param year Integer - start year of the disease season
  #' @param mod_level Integer - Spatial level of the model data. For CDC can only be 2, 3 and 4. For Dengue - country dependent
  #' @param fit_level Integer - Spatial level of data used for a coupled or uncoupled fit of the model data, fit_level = 2,3,4 for flu
  #'@param mod_name Named vector of strings specifying the model-level spatial patch.  If \code{is.null(mod_name)}, the code reverts to using \code{RegState} (see next entry).  To specify New York state, set \code{mod_name=c(NAME_2="United States", NAME_3="R1", NAME_4="New York"). Here NAME_X is either the full name or abbreviation of the level-X patch. Replacing 'United States' with 'US' or 'R1' with 'Region 1' would result in the same outcome.  Also, vector entries for mod_name should go from NAME_2,....,NAME_n where mod_level=n.}
  #' @param fit_name A character vector indicating which fit-regions to use.  If \code{fit_name='all'}, then DICE uses all child-regions of the model region with level equal to \code{fit_level}. The other mode for fit_name is to specifiy a subset of the fit regions to construct an aggregate representation of the model region.  For example if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=3}, and \code{fit_names=c("R1", "R2", "R3")}, DICE will create an Atlantic super-region to model (as opposed to using all 10 HHS regions).  Similarly, if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=4}, and \code{fit_names=c("WA", "OR", "CA")}, DICE will create and model a super-state of Pacific states.
  #' @param RegState  Single element: determines which single region from \code{mod_level} is to be modeled. Depending on the model level, \code{RegState} should adhere to the following format: \code{mod_level = 2} - 3-letter ISO3 RegState code, \code{mod_level=3} - an integer describing the HHS region, \code{mod_level=4} - a 2-letter state code. Where possible, RegState should be replaced by mod_name.  RegState is limited to country-level data and US regions/states.
  #' @param raw_col A string specifying which data column to use for modeling.  Data column names are associated with data sources and are listed in the data_sources table under 'col_name'.  When raw_col is NULL, DICE uses the default (first) data column.
  #' @param nfit Integer - Number of data points that will be fitted.  Default is to fit all the data.  This will be reset if nfit > nperiodsData
  #' @param method String either 'mech' for compartmental mechanistic models or 'stat' for SARIMA models
  #' @param epi_model String , name of mechanistic compartmental model: SIR, SEIR, (SIR)H/(SI)V,  (SEIR)H/(SEI)/V SIRB integer 1,2,3,4,5 (case insensitive)
  #' @param model Integer - The model number, see manual for more details (1-4 are supported for flu 4 for dengue). Relevant only when method = 'mech'
  #' @param Tg - recovery time in days.  If NULL it is set to three/eight days for flu/dengue. Relevant only when method = 'mech'
  #' @param isingle Integer - 0: couple the fit spatial regions; 1: no coupling. Relevant only when method = 'mech'
  #' @param prior Integer - if greater than zero use a prior for the MCMC procedure. Relevant only when method = 'mech'
  #' @param Temp Integer 1, 5, 10, 100 - Temperature for the MCMC procedure.    Relevant only when method = 'mech'
  #' @param da Integer 0, 1 or 2. Data augmentation options: 0-none, 1-using historic average and 2-using the most similar season.
  #' Relevant only when method = 'mech'
  #' @param nMCMC Integer - number of steps/trials in the MCMC procedure. Relevant only when method = 'mech'
  #' @param nreal Integer - number of MCMC chains. Relevant only when method = 'mech'
  #' @param arima_model  - A List of ARIMA model parameters: list(p=, d=, q=, P=, D, Q=) can be set to NULL to trigger the
  #' \code{auto.arima} process
  #' @param covar String, optional.  Covariate for use in ARIMA fitting.  Options are: 'sh', 'precip', 'temp'
  #' @param covar_lag Numeric lag time for optional covariate variable in time units of cadence of the data
  #' @param device  Either 'pdf' (default) or 'x11'
  #' @param plot TRUE, FALSE or EXTERNAL (or 0, 1, 2) allows the Users to implement their own plotting routines
  #' @param subDir Name of output sub-directory where all plots and files will be written.  Default is NULL -let the code build it.
  #' reproducible.
  #' @param all_data The complete DICE data structure returned by get.DICE.data().  This input allows the user to input the data rather than making a call to the DICE database.  If this input is not NULL, DICE will attempt to use the data assigned to all_data.


  #' @return solution a list with the input and entire output of the run.
  #' @examples
  #' For a run of the 2015-2016 cdc national mydata using the ten HHS regions with coupling between the regions use:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0)
  #'
  #' For a run of the 2015-2016 cdc national mydata using the ten HHS regions without coupling between the regions use:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 1)
  #'
  #' For a run of the 2014-2015 GFT mydata for HHS region number 9, using state level mydata with coupling between the states in region 9 use:
  #' output <- runDICE(data_source='gft', year = 2014, mod_level = 3, fit_level = 4, RegState = 9, isingle = 0)
  #'
  #' To control which model is used for the basic reproduction number, set the parameter model in your call. Default value is 5:
  #' output <- runDICE(data_source='gft', year = 2014, mod_level = 3, fit_level = 4, RegState = 9, isingle = 0, model = 3)
  #'
  #' To control the number of MCMC chains that the code will run set the parameter  nreal in your call, default is 1:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nreal = 3)
  #'
  #' To control the number of MCMC steps/trial in each chain set the parameter nMCMC in your call, default is 1e5:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6)
  #'
  #' To control the name of the sub-directory where all the output files and plots are saved use the keyword subDir, default is output:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6, subDir = 'test')
  #'
  #' To control the file format for the plots (pdf, png or x11) set the parameter device:
  #' output <- runDICE(data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6, device = 'pdf')
  #' (The package can accept an array of file formats, i.e. device = c('pdf','png'), in which case more both 'png' and 'pdf' files will be created.)
  #'
  #' To run in a forecast or predictive mode you can set the number of weeks the code uses in the fit to be lower than the number of weeks in the season.
  #' (Note that for the current season it is always running in a predictive mode because the season is not yet completed.)
  #' output <- runDICE(data_source='gft', year = 2013, mod_level = 3, fit_level = 4, isingle = 1, nMCMC = 1e6, nfit = 35)
  #'
  #' To select only a few HHS regions and run them coupled (for example the Eastern Regions 1, 2 and 3) use:
  #' output <- runDICE(data_source='cdc', year=2015, mod_level=2, fit_level=3, RegState=c('Region1','Region2','Region3'), isingle = 0)
  #'
  #'To select only a few states and run them coupled  use for example:
  #' output <- runDICE(data_source='gft',y ear=2014, mod_level=3, fit_level=4, RegState=c('WA','OR','CA'), isingle = 0)
  #'
  #' -- Data diseases and data_sources -------
  #' Access the database and list all available diseases:
  #' library(DICE)
  #' myDB = OpenCon()
  #' data_sources = dbReadTable(con=myDB, name="data_sources")
  #' unique(data_sources$disease)
  #' # then list all data sources
  #' str(data_sources)
  #' data_sources$source_abbv
  #' dbDisconnect(myDB)



  ## Method - Statistical or Mechanistic
  if(is.null(method)) {
    method = 'mech'
  }

  ## In the interface to the USer we call it nfit because nperiodsFit is too long
  nperiodsFit = nfit


  if (is.null(all_data)) {
    ## if data_source has been entered as a number, convert to abbreviation
    if (!is.null(data_source)) {
      if (is.numeric(data_source)) {
        if (!db_opts$CDC_server) {
          myDB = OpenCon(db_opts)
          data_sources = dbReadTable(myDB, "data_sources")
          dbDisconnect(myDB)
          if (data_source %in% data_sources$source_key) {
            data_source = tolower(data_sources$source_abbv[data_source])
          } else {
            cat("Warning: number entered for data_source does not match a DICE source_key.  DICE will now attempt to choose a default data source.")
            data_source = NULL
          }
        } else {
          if (data_source %in% DICE_data_sources$source_key) {
            data_source = tolower(DICE_data_sources$source_abbv[data_source])
          } else {
            cat("Warning: number entered for data_source does not match a DICE source_key.  Setting data_source='cdc'.")
            data_source = "cdc"
          }
        }
      } else {
        # check that source abbreviation matches a DICE source
        if (!db_opts$CDC_server) {
          myDB = OpenCon(db_opts)
          data_sources = dbReadTable(myDB, "data_sources")
          dbDisconnect(myDB)
          if (tolower(data_source) %in% tolower(data_sources$source_abbv)) {
            data_source = tolower(data_source)
          } else {
            cat("Warning: string entered for data_source does not match a DICE source_abbv.  DICE will now attempt to choose a default data source.")
            data_source = NULL
          }
        }
      }
    }
  } else {
    disease = all_data$mydata$disease
    mod_level = all_data$mydata$model$level
    fit_level = all_data$mydata$fit$level
    year = all_data$mydata$season
    data_source = all_data$mydata$data_source
  }

  ## data_source: if sql_db==F, cdc is the only current option
  if (is.null(data_source) & db_opts$CDC_server) {
    if (disease == 'flu') data_source = "cdc"
  }

  ## Tg in days
 if (is.null(Tg)) {
 	Tg = 3
 	if (tolower(disease) == "dengue")
 		Tg = 8

 	if (tolower(disease) == "yellow_fever")
 		Tg = 5

 	if (tolower(disease) == "plague")
 		Tg = 5

 	if (tolower(disease) == "chik")
 		Tg = 8

 	if (tolower(disease) == "zika")
 		Tg = 7

 	if (tolower(disease) == "sars")
 		Tg = 4
 }


  if (is.null(epi_model)) {
    epi_model = 1
  } else if (tolower(epi_model) == "sir" | epi_model == 1) {
    epi_model = 1
  } else if (tolower(epi_model) == "seir" | epi_model == 2) {
    epi_model = 2
  } else if (tolower(epi_model) == "vsir" | epi_model == 3) {
    epi_model = 3
  } else if (tolower(epi_model) == "vseir" | epi_model == 4) {
    epi_model = 4
  } else if (tolower(epi_model) == "sirb" | epi_model == 5) {
    epi_model = 5
  } else {
    epi_model = 1
  }

  ## Takes care of the user giving a NULL value from the command line.
  ## R interprets is as a string and NOT as an object of type NULL
  ## This will convert the NULL string to a NULL object

  if (!is.null(arima_model) && (!is.list(arima_model) | arima_model[1] == "NULL")) {
    arima_model = NULL
    covar = FALSE
  }

  ## Number of MCMC chains
  if (is.null(nreal))
    nreal = 1

  ## Number of steps/trials in each MCMC chain

  if (is.null(nMCMC))
    nMCMC = 10000

  ## Number of times the history of the MCMC chain is saved.

  nlines = round(nMCMC/100)
  nlines = min(10000, nlines)
  nlines = max(100, nlines)

  if (is.null(mod_level))
    mod_level = 2

  if (is.null(fit_level))
    fit_level = 3
  ## Start year of the flu season

  if (is.null(year))  {
    year = 2010
    if (disease == 'flu') year = 2018
  }


  ## Model Number
  if (is.null(model)) {
    model = 4
    if(disease == 'flu') model = 3
  }


  ## Number of periods of mydata that are fitted
  if (is.null(nperiodsFit))  {
    nperiodsFit = 12
    if (disease == 'flu') nperiodsFit = 52
  }

  ## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the spatial regions

  if (is.null(isingle))
    isingle = 1

  ## device name for plotting the results - pdf or X11
  if (is.null(device))
    device = "pdf"

  ## Use (prior > 0) a prior in the MCMC procedure, or not (prior = 0)
  if (is.null(prior))
    prior = 0

  ## 'Temperature' for the prior width. Used only if prior > 0

  if (is.null(Temp))
    Temp = 1

  ## mydata augmentation
  if (is.null(da))
    da = 0

  # One of these should be filled.  If not, set to something that will probably work.
  if (is.null(RegState) & is.null(mod_name) & is.null(all_data)) {
    if (mod_level == 2) {
      RegState = "usa"
      if (tolower(disease) == 'dengue' | tolower(disease) == 'yellow_fever') RegState = 'BR'
      if (tolower(disease) == 'zika') RegState = 'PR'
    }
    if (mod_level == 3) {
      RegState = 9
      if (tolower(disease) == 'dengue') RegState = 1
    }
  }

  if (is.null(plot))
    plot = 1

  if (tolower(disease) != "flu" ) {

    if (model == 2 | model == 3) {
      cat("\nSchool Vacation Models can not be used for Dengue Data \n\n  \tResetting the model to 4 - Fixed Force of Infection !!\n\n")
      model = 4
    }


    if (isingle == 0 & epi_model > 2 & fit_level > mod_level) {
      cat("\n Explicit Modeling of Vector States is Supported only for Uncoupled Runs, Resetting to SIR or SEIR with quasi-equilibrium \n
		assumption for Vector States\n")
      if (epi_model == 3) epi_model = 1
      if (epi_model == 4) epi_model = 2
    }
  }

  if (is.null(all_data)) {
    # Get data. Default behavior is to return the specified season in $mydata as well as all available data in $all_years_epi.
    complete_mydata <- get.DICE.data(data_source = data_source, mod_level = mod_level, fit_level = fit_level, year = year, model = model, nperiodsFit = nperiodsFit, mod_name=mod_name, RegState = RegState, fit_names=fit_names, isingle = isingle, db_opts=db_opts, disease = disease, epi_model = epi_model, method = method, all_years_flag=T, all_cad_clim=T, raw_col=raw_col)
  } else {
    complete_mydata = all_data
  }

  mydata = complete_mydata$mydata # season for modeling
  all_years_epi = complete_mydata$all_years_epi

  nmydata = length(all_years_epi$years)
  all_years_epi$FY = paste0(all_years_epi$years[1], "-", all_years_epi$years[nmydata])

  all_years_epi$mydataName = paste0(mydata$model$name, "-", all_years_epi$FY)

  mydata$sql_db = !db_opts$CDC_server

  mydata$db_opts = db_opts

  if (is.null(data_source)) {
    data_source = mydata$data_source
  }

  ## SIR (1) or SEIR (2) Vector-SIR (3) Vector SEIR (4) SIRB (5)

  mydata$epi_model = epi_model

  mydata$single = isingle


  ## prior, da and temperature for MCMC procedure

  if (data_source == 'cdc') {
    mydata$prior = prior # We only have a prior in the case of cdc data
  }  else {
    mydata$prior = 0
  }

  mydata$da = da
  mydata$Temp = Temp
  mydata$fit_level = fit_level
  mydata$mod_level = mod_level

  # if not specified, auto-generate a directory name
  if (is.null(subDir)) {
  	myName = mydata$model$name
   	myName = gsub(" ","",myName)
    if (method == 'mech') {
      subDir = paste0(myName, "_", mydata$season, "_level_",mod_level,'_prior_',prior,'_Temp_',Temp,'_da_',da, "_Tg_",Tg,"/")
    } else {
      subDir = paste0(myName, "_", mydata$season, "_level_",mod_level)
    }
  }

  mydata$subDir = subDir
  # if necessary create directory
  if (!dir.exists(subDir)) {
    dir.create(subDir)
  }

  mydata$method = method
  mydata$Tg = Tg

  ##
  ## covariate (in SARIMA case)
  ##

  if (is.null(covar))
    covar = FALSE

  if (is.null(covar_lag))
    covar_lag = 0

  mydata$covar = covar
  mydata$covar_lag = covar_lag


  ##
  ## mydata Augmentation - 0 (no), 1 (historic NULL model), 2 (Most Similar Season) In sql_db case this is done later
  ##

  if (disease == "flu" & data_source == 'cdc' & method == 'mech') {
    
    mydata <- get.cdc.prior(mydata = mydata)

    mydata <- get.cdc.da(mydata = mydata)

  }

  # retrieve flu point-of-care units
  flu_poc_units = get_flu_poc_units()

  if (disease == "flu" & mydata$model$raw_units%in%flu_poc_units & method == 'mech') {
    mydata <- get.flu_poc.da(mydata = mydata, all_years_epi = all_years_epi)
  }

  ## DA for San Diego
  if ( length(mod_name) == 5 && mod_name[5] == "SD") {
  	mydata <- get.sd.da(mydata = mydata)
  }


  ##
  ## Plot the entire incidence data
  ##

  if (plot == TRUE || plot == 1 || plot == 2) {
  	err <- plotDisease(mydata = mydata, all_years_epi = all_years_epi, device = device)
  }

  ## SARIMA CASE
  if (tolower(method) == 'stat') {
  	if (tolower(disease) == 'ebola' || tolower(disease) == 'cholera' || tolower(disease) == 'plague' || tolower(disease) == 'sars')
  	stop("SARIMA Statistical Modeling is NOT supported for Ebola/Cholera/Plague/SARS")

    auto_arima_model = NULL

    if (is.null(arima_model)) {  ## If User has not provided a model these are the maximum values we allow..

      if(disease == 'flu') auto_arima_model = list(p = 0, d = 1, q = 2, P = 0, D = 1, Q = 2)
      if (tolower(disease) == 'dengue' || tolower(disease) == 'yellow_fever') {
        auto_arima_model = list(p = 1, d = 0, q = 0, P = 2, D = 1, Q = 0)
        if (mydata$cadence == "Monthly")
          auto_arima_model = list(p = 1, d = 0, q = 0, P = 3, D = 1, Q = 0)
      }


    }
    # Special Code for SD county for now - may change this
    if (length(mod_name) == 5 && mod_name[5] == "SD") {
 		solution = fitSARIMASD(mydata = mydata, all_years_epi = all_years_epi, arima_model = arima_model, auto_arima_model = auto_arima_model, device = device)
    } else {
    	solution = fitSARIMA(mydata = mydata, all_years_epi = all_years_epi, arima_model = arima_model, auto_arima_model = auto_arima_model, device = device)
    }

    return(solution)
  }

  ##
  ## Mechanistic Modeling
  ## Pack the information for the run

  par_names <- set.param.list(epi_model = mydata$epi_model)

  opt.list <- set.opt.list(mydata = mydata, isingle = isingle)

  run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, nlines = nlines, device = device, subDir = subDir, plot = plot)


  ##
  ## Fitting a SINGLE region/patch
  ##


  if (mydata$fit$nregions == 1) {

    cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
    cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$disease), " data", "\n")
    cat("\n\n Fitting ", mydata$nperiodsFit, " weeks out of ", mydata$nperiodsData, " periods of data", "\n\n")
    cat("\n\n Fitting ", toupper(mydata$model$name), " data", "\n")

    for (ireal in 1:nreal) {
      if (tolower(mydata$data_source) == "who_flu") {
        solution = fitOnePatchWHO(mydata = mydata, all_years_epi = all_years_epi, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
      } else if (tolower(mydata$disease) == 'cholera') {
        solution = fitCholera(mydata = mydata, all_years_epi = all_years_epi, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
      } else if (tolower(mydata$disease) == 'plague') {
        solution = fitPlague(mydata = mydata, all_years_epi = all_years_epi, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
      } else if (tolower(mydata$disease) == 'sars') {
       solution = fitSARS(mydata = mydata, all_years_epi = all_years_epi, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)

      }
      else {
        solution = fitOnePatch(mydata = mydata, all_years_epi = all_years_epi, opt.list = opt.list, run.list = run.list, ireal = ireal,
                               iseed = iseed)
      }
    }

    return(solution)
  }

  ## ALL OTHER CASES - number of regions/patches of fit_level > 1

  cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
  cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$disease), " data", "\n")
  cat("\n\n Fitting ", mydata$nperiodsFit, " weeks out of ", mydata$nperiodsData, " periods of data", "\n\n")
  cat("\n\n Fitting ", toupper(mydata$model$name), " using ", toupper(mydata$fit$name), " data", "\n")

  ## UNCOUPLED CASE

  if (isingle == TRUE) {
    for (ireal in 1:nreal) {

      if (disease == "flu" & tolower(mydata$data_source) == 'who_flu') {
        solution = fitSingleWHO(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal,
                                iseed = iseed)
      } else {
        solution = fitSingle(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal,
                             iseed = iseed)
      }

    }
    ## COUPLED CASE

  } else {
    for (ireal in 1:nreal) {
      if (disease == "flu" & tolower(mydata$data_source) == 'who_flu') {
        solution = fitMultiWHO(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal,
                               iseed = iseed)
      } else {
        solution = fitMulti(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal,
                            iseed = iseed)

      }

    }
  }


  return(solution)
}


fitOnePatch <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine - Fitting a Single Incidence Profile
    #'
    #' A driver for fitting a single incidence profile of a region/patch. The R code calls a Fortran routine which
    #' uses an MCMC procedure to maximize the likelihood of the solution using SIR, SEIR, (SIR)_H/(SI)_V or (SEIR)_H/(SEI)_V models.
    #' @param mydata A list containing all available information for the region: incidence, specific-humidity, school schedule, population etc.
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    #' These values are set based on the user chosen model for the compartmental model and the force of infection.
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param ireal Integer - the MCMC chain number.  Default is 1.
    #' @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly.
    #' Setting the seed to a known integer allows an MCMC chain to be reproducible.
    #' @return A list with the following arguments:
    #' \describe{
    #' \item{model_rtn}{The best result of the MCMC procedure for the data}
    #' \item{model_profile}{Randomly selected results from the MCMC procedure of fitting the data}
    #' \item{tab.model}{The MCMC history of fitting the data}
    #' }
    #' @examples
    #' fitSingle{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}


	if (mydata$cadence == "Weekly") {
		cadence = 7
	} else if (mydata$cadence == "Monthly") {
		cadence = 31
	} else if (mydata$cadence == "Daily") {
		cadence = 1
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}

    epi_model = mydata$epi_model

    # Take only the first list from opt.list, the second is optional and is needed only for a coupled run-irrelevant in this case

    opt.list = opt.list[[1]]

    device = run.list$device
    subDir = run.list$subDir

    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    prior = mydata$prior
    Temp = mydata$Temp
    n = mydata$fit$nregions

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed() * ireal%%.Machine$integer.max
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generation seed for R
    set.seed(iseed)


	mod_level = mydata$mod_level
	fit_level = mydata$fit_level


    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    mod_name = mydata$model$name
    fit_name = mydata$fit$name

    tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

	##  Observed number of cases

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)

    ## Historic NULL model

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod

	tables.mod$epi.null[nperiodsFit, , mod_id] = epi.null.mod[, mod_id]

	if (tolower(mydata$disease) == "flu" | tolower(mydata$disease) == "chik" ||
		tolower(mydata$disease) == "ebola") {
		week0 = mydata$weeks[1]
		day0 = cadence * (week0 - 1)
		tps = seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}
	if (tolower(mydata$disease) == "zika" || tolower(mydata$disease) == "sars") {
		tps = (1:mydata$nperiods) * cadence
	}

	##
	## Missing data is given a zero weight
	##
	if (tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == "zika" || tolower(mydata$disease) == "ebola") {
		wght = rep(0, mydata$nperiods)
		wght[1:nperiodsFit] = 1
		mydata$model$wght = wght
		mydata$model$wght[is.na(mydata$model$raw)] = 0
		mydata$fit$wght = mydata$fit$epi
		for (iregion in 1:mydata$fit$nregions) {
			mydata$fit$wght[,iregion] = wght
			mydata$fit$wght[is.na(mydata$fit$raw[,iregion]), iregion] = 0
		}
	}

	if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == "yellow_fever") {
		## Fit the mod_level
		## Start by calculating correlation/distance with past years

		corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		epi.da.df = get.sql.da(nfit = nperiodsFit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
			corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

		mydata$model$wght = epi.da.df$wght

		mydata$model$epirun = epi.da.df$epi.new

		gamaepi = rep(0, nperiods)
		for (i in 1:nperiods) {
			gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
		}

		mydata$model$gamaepirun = gamaepi

		tps = mydata$ndays

	}

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    if (tolower(mydata$disease) == "flu" || tolower(mydata$disease) == 'sars' || tolower(mydata$disease) == 'ebola') {
    	setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    } else if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == "yellow_fever" || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'zika') {
    	setup = setup.vector.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    }


    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    nparam = setup$nparam
    nopt = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec

    tab = setup$tab

    nreal = run.list$nreal
    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

	setup.prior = setup.single.prior(mydata = mydata, logvec = logvec)

	n1 = length(setup.prior$ymu[, 1])

	ymu = setup.prior$ymu[n1, ]
	sigma = setup.prior$sigma[n1, ]

	npad = nparam - length(ymu)
	if (npad > 0) {
		ymu = c(ymu, rep(0, npad))
		sigma = c(sigma, rep(1, npad))
		ymu = as.numeric(ymu)
		sigma = as.numeric(sigma)
	}

	if (mydata$prior == 1 | mydata$prior == 2 & tolower(mydata$data_source) == "cdc")
		logvec = setup.prior$logvec

	imask = set.imask(par_names = par_names, opt.list = opt.list)

	# The number of paramters that are optimized
	nopt = length(which(imask == 1))

	# Here loop on each region and call a single-region fitting routine
	tab.model = NULL

	pois = 0
	nblock = 3
	accept.vec = array(0, c(n, nblock))

	# retrieve flu point-of-care units
	flu_poc_units = get_flu_poc_units()

    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units || (mydata$model$attr$level==6 && mydata$model$attr$ABBV_4=="CA" && mydata$model$attr$ABBV_6 == "SD")) {

    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 2.0

    }

	if (tolower(mydata$disease) == "dengue" || (mydata$model$level==6 && mydata$model$attr$ABBV_4=="CA" && mydata$model$attr$ABBV_6=="SD")) {
		cases = mydata$model$epirun
		gamaepi = mydata$model$gamaepirun
	} else {
		cases = mydata$model$epi
		gamaepi = mydata$model$gamaepi
	}

	sh     = mydata$model$sh
	school = mydata$model$school
	wght   = mydata$model$wght

	nRnd = 1000
	model_profile = array(0, c(nRnd, nperiods))
	## Need nperiods + 2 because we pathc the weeks at the beginning and end

	nmydata = length(tps)
	if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever') {
		time = cumsum(tps)
		ndays = sum(tps) + tps[1] + tps[nmydata]
	} else {
		time = tps[1:nperiods]
		ndays = (nperiods + 2) * cadence
	}

	cat("\n ****** Fitting ", mydata$model$name, " Using an MCMC Procedure ****** \n")

	if (is.null(school)) school = rep(0,mydata$nperiods)

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(time),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))


     model_profile = array(solution$profile, c(nRnd, nperiods))
     model_rtn = solution$rtn

    tab.model =  matrix(solution$tab, ncol = (nparam + 1))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

 	tables = calc.null.err(tables = tables.mod, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

	model_factor = mydata$model$factor
	model_profile_ili = model_profile/model_factor
	model_rtn_ili = model_rtn/model_factor

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile_ili, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

    # # write the MCMC statistics

    success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, mydata = mydata, imask = imask, ireal = ireal)


    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

    ## Write CSV files: (1) with the ILI profiles and (2) a CSV file with onset week, week maximum, maximum value etc.

	if (mydata$data_source == 'cdc' || mydata$model$raw_units%in%flu_poc_units) {
	    success = writecsvOnePatch(mydata = mydata, run.list = run.list, tab = tab.model, model_rtn = model_rtn, model_profile = model_profile,
        ireal = ireal)
    }
    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

     tables = calc.null.err(tables = tables.mod, nfit = nperiodsFit, state_id = mod_id)
	 tables.mod = tables

	 tables = calc.mech.err(tables = tables.mod, profiles = model_profile_ili, nfit = nperiodsFit, state_id = mod_id)
	 tables.mod = tables

    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
    	if (as.character(run.list$plot) == "TRUE" || as.character(run.list$plot) == "1") {
    		for (k in 1:length(run.list$device)) success = plotFitOnePatch(model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
    			ireal = ireal, run.list = run.list, idevice = k)
    	}

    	if (tolower(as.character(run.list$plot)) == "external" || as.character(run.list$plot) == "2") {
    		for (k in 1:length(run.list$device)) success = plotFitOnePatch.external(model_rtn = model_rtn, model_profile = model_profile,
    			mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)
    	}

    } else {
    	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id,
    		ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)

    	err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)

    ##
    ## now dump all the profiles we have to a file    - in the case of the CDC this is done in the ploting routine
    ##
    	err <- saveProfilesOnePatch( mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list)

    }

	##
	## Plot the posterior distribution of the parameters
	##
    err <- plotMCMC(mydata = mydata, tab.model = tab.model, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)


    results = list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model)

    return(results)

}

fitMulti <- function(mydata = NULL,  all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine - Coupled Spatial Model
    #'
    #' A spatially coupled MCMC fit of the model data. The code first fits the model data directly and then fits it
    #'   as a weighted sum of the coupled fit level data. This fit uses a coupling matrix to describe the interaction between
    #'   different spatial regions and it generates all the fit level profiles at once and minimizes their weighted
    #'   likelihood with the weights given by the relative population of each region. The data can be either cdc or gft data,
    #'   and the model/fit data should have different spatial scales. For example in the case of cdc data: the model is
    #'   national and the fit are the ten HHS regions. Or the model can be an HHS region and the fit fit data is state level data.
    #' @param mydata A complex list with all available data for a given disease season model/fit spatial levels and data type
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    #' These values are based on the user's chosen compartmental model (SIR/SEIR) and the force of infection.
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param ireal Integer - the MCMC chain number.  Default is 1.
    #' @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly.
    #' Setting the seed to a known integer allows an MCMC chain to be reproducible.
    #' @return A list with the following arguments:
    #' \describe{
    #' \item{model_rtn}{The best result of the MCMC procedure for the model data}
    #' \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    #' \item{tab.model}{The MCMC history of the direct fit to the model data}
    #' \item{rtn}{The best result for indirectly fitting the model data using the fit regions}
    #' \item{profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    #' \item{tab}{The MCMC history of indirectly fitting the model data using the fit regions}
    #' }
    #' @examples
    #' fitMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}

 	epi_model = mydata$epi_model

	if (mydata$cadence == "Weekly") {
		cadence = 7
		weeks = mydata$weeks
	} else if (mydata$cadence == "Monthly") {
		cadence = 31
		months = mydata$months
	} else if (mydata$cadence == "Daily") {
		cadence = 1
		days = mydata$days
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed() * ireal%%.Machine$integer.max
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generator seed for R
	set.seed(iseed)


    # Take both the first and the second elements for opt.list.  The first is always present The second is needed only in case of a
    # coupled run
    opt.cpl = opt.list[[2]]
    opt.list = opt.list[[1]]

	nreal  = run.list$nreal
	nMCMC  = run.list$nMCMC
	nlines = run.list$nlines
	ithin = run.list$ithin
	device = run.list$device
	subDir = run.list$subDir


   	nperiods = mydata$nperiods
   	nperiodsFit = mydata$nperiodsFit

   	prior = mydata$prior
   	Temp = mydata$Temp
   	n = mydata$fit$nregions

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    mod_name = mydata$model$name
    fit_name = mydata$fit$name

    mod_level = mydata$mod_level
    fit_level = mydata$fit_level

    tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)

   if (mod_level < fit_level) {

   	tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)
   }

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.numeric(epi.obsrv)
	epi.obsrv = mydata$fit$raw

	tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod
	epi.null.fit = epi.null$cases.ave.fit

	tables.mod$epi.null[nperiodsFit, , mod_id] = epi.null.mod[, mod_id]

	if (mod_level < fit_level) {

		tables.fit$epi.null[nperiodsFit, , fit_id] = epi.null.fit[, fit_id]
	}

	if (tolower(mydata$disease) == "dengue") {
		## alculating correlation/distance with past years and augmenting appropriately

		corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		epi.da.df = get.sql.da(nfit = nperiodsFit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
			corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

		mydata$model$wght = epi.da.df$wght

		mydata$model$epirun = epi.da.df$epi.new

		gamaepi = rep(0, nperiods)
		for (i in 1:nperiods) {
			gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
		}

		mydata$model$gamaepirun = gamaepi

		tps = mydata$ndays

		## Repeat for fit_level


		if (fit_level > mod_level) {


			mydata$fit$epirun = mydata$fit$epi

			mydata$fit$gamaepirun = mydata$fit$gamaepi

			mydata$fit$wght = mydata$fit$epi

			for (iregion in 1:mydata$fit$nregions) {

				corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nperiodsFit)

				distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nperiodsFit)

				epi.da.df = get.sql.da(nfit = nperiodsFit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$fit$epi[, iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

				mydata$fit$wght[, iregion] = epi.da.df$wght

				mydata$fit$epirun[, iregion] = epi.da.df$epi.new

				for (i in 1:length(mydata$fit$epirun[, iregion])) {
					gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
				}
				mydata$fit$gamaepirun[, iregion] = gamaepi

			}
		}

	}


	if (tolower(mydata$disease) == "flu" || tolower(mydata$disease) == "chik" || tolower(mydata$disease) == "ebola") {
    	week0 = mydata$weeks[1]
    	day0 = cadence * (week0 - 1)
	    tps =  seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}

	if (tolower(mydata$disease) == "zika") {
		tps = (1:mydata$nperiods) * cadence
	}

	## In the case of Chikungunya/Zika data have to take care of the weights and tps
	## Missing data is given a zero weight
	##
	if (tolower(mydata$disease) == "chik" || tolower(mydata$disease) == "zika" || tolower(mydata$disease) == "ebola") {
		wght = rep(0, mydata$nperiods)
		wght[1:nperiodsFit] = 1
		mydata$model$wght = wght
		mydata$model$wght[is.na(mydata$model$raw)] = 0
		mydata$fit$wght = mydata$fit$epi
		for (iregion in 1:mydata$fit$nregions) {
			mydata$fit$wght[,iregion] = wght
			mydata$fit$wght[is.na(mydata$fit$raw[,iregion]), iregion] = 0
		}
	}

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    par_names_cpl <- set.param.cpl.list()

    if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever' || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'zika') {
   	setup = setup.vector.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

   } else {

   	setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

   }

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    ## Just for the CDC prediction limit R0 nd pC
    # retrieve flu poc_units
    flu_poc_units = get_flu_poc_units()
    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {

    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 2.0

    }

     # these are the same for both model and fit

    nparam = setup$nparam
    logbase = setup$logbase
    logvec = setup$logvec

    tab = setup$tab
    nparam = length(model_par)

	## Setting up for Fit level

    if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever' || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'zika') {
    	setup = setup.vector.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    } else {
    	setup = setup.fit.mcmc(       mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    }

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par

   ## Just for the CDC prediction limit pC and R0

   if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
   	fit_pmax["R0", 1:n] = 2.0
   	fit_pmax["pC", 1:n] = 0.2

   }


       ## And let's sample pC and R0 from the prior

    ## We get the mean and sigma of the prior but we actually use it only with prior = 1 or 2 This is controlled by logvec see if
    ## statement below

    setup.prior = setup.single.prior(mydata = mydata, logvec = logvec)

    n1 = length(setup.prior$ymu[, 1])

    fit_ymu = setup.prior$ymu[1:n, ]
    fit_sigma = setup.prior$sigma[1:n, ]
    model_ymu = setup.prior$ymu[n1, ]
    model_sigma = setup.prior$sigma[n1, ]
	#logvec  = setup.prior$logvec


    if (mydata$fit_level < 4 & mydata$mod_level < 4 & mydata$imodel != 5 & mydata$data_source == "cdc") {

    	prior.ini <- sample.from.prior(mydata = mydata, model_par = model_par, fit_par = fit_par, model_pmin = model_pmin, model_pmax = model_pmax, fit_pmin = fit_pmin, fit_pmax = fit_pmax, model_ymu = model_ymu, model_sigma = model_sigma, fit_ymu = fit_ymu, fit_sigma = fit_sigma)

		model_par = prior.ini$model_par
		fit_par   = prior.ini$fit_par

    }

    if (mydata$prior == 1 || mydata$prior == 2 & tolower(mydata$data_source) == "cdc")
        logvec = setup.prior$logvec

    imask = set.imask(par_names = par_names, opt.list = opt.list)

   # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    nfit.prior = dim(fit_ymu)[2]

    setup = setup.coupling.mcmc(opt.list = opt.cpl)

    cpl_pmin = setup$parmin
    cpl_pmax = setup$parmax
    cpl_dx = setup$pardx
    cpl_par = setup$par
    cpl_nparam = setup$nparam
    cpl_nopt = setup$nopt
    cpl_logvec = setup$logvec

    ## Build the Rij matrix - there are zero's on the diagonal

    Rij = distance.matrix(mydata = mydata)

    cpl_imask = set.cpl.imask(nparam = cpl_nparam, opt.list = opt.cpl)

    # Here loop on each region and call a single-region fitting routine

    pois = 0
    nblock = 3
    accept.vec = array(0, c(n, nblock))
  	### Get the incidence, gamaepi, sh and school

 	if (tolower(mydata$disease) == "dengue") {
		cases = mydata$model$epirun
		gamaepi = mydata$model$gamaepirun
	} else {
		cases = mydata$model$epi
		gamaepi = mydata$model$gamaepi
	}

	sh = mydata$model$sh
	school = mydata$model$school
	wght = mydata$model$wght

    mypar = model_par

    ## Select the prior for this model region
    ymu = model_ymu
    sigma = model_sigma
    ## pad with zeros
    npad = nparam - length(ymu)
    if (npad > 0) {
        ymu = c(ymu, rep(0, npad))
        sigma = c(sigma, rep(1, npad))
        ymu = as.numeric(ymu)
        sigma = as.numeric(sigma)
    }

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

    profile = array(data = 0, dim = c(nRnd, nperiods, n))

	nmydata = length(tps)

	if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever') {
		time = cumsum(tps)
		ndays = sum(tps) + tps[1] + tps[nmydata]
	} else {
		time = tps[1:nperiods]
		ndays = (nperiods + 2) * cadence
	}

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

# c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd", "deltaR", "aparam", "alpha", "delta", "ts", "dur")

	if (is.null(school)) school = rep(0,mydata$nperiods)

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), Temp = as.double(Temp),
        imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(time[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd),epi_model=as.integer(epi_model))

    cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")


    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile, c(nRnd, nperiods))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

	tables = calc.null.err(tables = tables.mod, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

    ## Now fit the coupled regions

	if (tolower(mydata$disease)  == 'dengue') {
		 cases = as.matrix(mydata$fit$epirun)
   		 gamaepi = as.matrix(mydata$fit$gamaepirun)
   		 model_gamaepi = as.numeric(mydata$model$gamaepirun)
	} else {
    	cases = as.matrix(mydata$fit$epi)
    	gamaepi = as.matrix(mydata$fit$gamaepi)
    	model_gamaepi = as.numeric(mydata$model$gamaepi)
	}

    sh = as.matrix(mydata$fit$sh)

    if (is.null(mydata$fit$school)) {
    	school = matrix(data=0, nrow=mydata$nperiods, ncol=mydata$fit$nregions)
    } else {
    	school = as.matrix(mydata$fit$school)
    }

    mypar = fit_par
    wght = as.matrix(mydata$fit$wght)

    rtn = array(0, c(nperiods, n))
    tab = array(0, c(nlines, (nparam + 1), n))
    tab.cpl = array(0, c(nlines, 2))

    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")


	for (iregion in 1:n) {
		# Also check for NAs
		school[is.na(school[, iregion]), iregion] <- 0
		sh[is.na(sh[, iregion]), iregion]         <- 0
	}

    solution = .Fortran("epimulti", n = as.integer(n), epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin), parmax = as.double(fit_pmax),
        dx = as.double(fit_dx), ilog = as.integer(logvec), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
        ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(time[1:nperiods]),
        rtn = as.double(rtn), pois = as.double(rep(0, n)), coef = as.double(mydata$fit$coef),
        parCPL = as.double(cpl_par), pminCPL = as.double(cpl_pmin),
        pmaxCPL = as.double(cpl_pmax), stepCPL = as.double(cpl_dx), ilogCPL = as.integer(cpl_logvec), imaskCPL = as.integer(cpl_imask),
        Rij = as.double(Rij), tab = as.single(tab), tabCPL = as.single(tab.cpl), profiles = as.single(profile), nRnd = as.integer(nRnd), ndays = as.integer(ndays), epi_model=as.integer(epi_model), scales = as.double(mydata$Temp))


    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fit$name, " Data Completed \n")

    rtn = matrix(solution$rtn, ncol = n)

    tab.cpl = matrix(solution$tabCPL, ncol = 2)

    tab = array(solution$tab, c(nlines, (nparam + 1), n))


    profile = array(solution$profiles, c(nRnd, nperiods, n))


	## Convert LLK to AICc and change tab to a list
	tab.fit = list()
	for (iregion in 1:n) {
		tab.tmp = tab[,,iregion]
		colnames(tab.tmp) = c(myColName, "AICc")
		tab.tmp[, "AICc"] <- 2 * tab.tmp[, "AICc"] + 2 * nopt
		tab.tmp[, "AICc"] <- tab.tmp[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)
		tab.fit[[iregion]] = tab.tmp

	}

	## Build the coupled fit to model level mydata. It is a simple sum since these are cases and not %ILI

    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {

        fit_model[i] = sum(rtn[i, 1:n])

        for (irnd in 1:nRnd) {
            for (k in 1:n) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile[irnd, i, k]
        }

	}

	## update table
	tables.agg.mod = tables.mod
	tables = calc.mech.err(tables = tables.agg.mod, profiles = fit_model_profile, nfit = nperiodsFit, state_id = mod_id)
	tables.agg.mod = tables

    success = mcmc.multi.write(tab.model = tab.model, tab.fit = tab.fit, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list,
        mydata = mydata, imask = imask, cpl_imask = cpl_imask, ireal = ireal)


 	tables = calc.null.err(tables = tables.fit, nfit = nperiodsFit, state_id = fit_id)
	tables.fit = tables

 	## Update the results table

 	for (iregion in 1:n) {
 		tables = calc.mech.err(tables = tables.fit, profiles = profile[,,iregion], nfit = nperiodsFit, state_id = fit_id[iregion])
		tables.fit = tables
 	}

     # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

	if (tolower(mydata$disease) == 'flu') {
    success = writeCSV(mydata = mydata, run.list = run.list, model_rtn = model_rtn, model_profile = model_profile,
        rtn = rtn, profile = profile, ireal = ireal)
    }

	##
    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
    	if (as.character(run.list$plot) == "TRUE" || as.character(run.list$plot) == "1") {
    		for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
    			ireal = ireal, run.list = run.list, idevice = k)

    		for (k in 1:length(run.list$device)) success = plotHists(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
    			ireal = ireal, run.list = run.list, idevice = k)
    	}

    	if (tolower(as.character(run.list$plot)) == "external" || as.character(run.list$plot) == "2") {
    		for (k in 1:length(run.list$device)) {
    			success = plotFitCDCPercentILI.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata, ireal = ireal,
    				run.list = run.list, idevice = k)
    		}
    		for (k in 1:length(run.list$device)) {
    			success = plotHists.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata, ireal = ireal, run.list = run.list,
    				idevice = k)
    		}


    	}

    } else {
    	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL,
    		ireal = ireal, run.list = run.list, idevice = 1)
    	err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)
    	err <- saveProfiles(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, rtn = rtn,  profile = profile,  ireal = ireal, run.list = run.list)
    }

    err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit= tab.fit, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)

    results = list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model, rtn = rtn, profile = profile, tab = tab)

    return(results)

}



fitSingle <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine - Uncoupled Spatial Model
    #'
    #' A spatially uncoupled MCMC fit of the model data. The code first fits the model data directly and then
    #'   fits each of the sub-regions sequentially - minimizing the likelihood of each one. The final indirect
    #'   model fit is obtained as a weighted sum of these individual fits with the weights given by the relative
    #'   population of each region. The data can be either cdc or gft data, and the model/fit data should have
    #'   different spatial scales. For example in the case of cdc/gft data: the model can be national and the fit are
    #'   the ten HHS regions. Or the model can be an HHS region and the fit are the states in that region.
    #' @param mydata A complex list with all available data for a given season model/fit spatial levels and data type
    #' @param all_years_epi A complex list with all available data for all seasons model/fit spatial levels and data type
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.   These values are set based on the user chosen model
    #' for the basic reproduction value.
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param ireal Integer - the MCMC chain number.  Default is 1.
    #' #' @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    #' @return A list with the following arguments:
    #' \describe{
    #' \item{model_rtn}{The best result of the MCMC procedure for the model data}
    #' \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    #' \item{tab.model}{The MCMC history of the direct fit to the model data}
    #' \item{fit_rtn}{The best result for indirectly fitting the model data using the fit regions}
    #' \item{fit_profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    #' \item{tab.fit}{The MCMC history of indirectly fitting the model data using the fit regions}
    #' }
    #' @examples
    #' fitSingle{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}

 	epi_model = mydata$epi_model

	if (mydata$cadence == "Weekly") {
		cadence = 7
		weeks = mydata$weeks
	} else if (mydata$cadence == "Monthly") {
		cadence = 31
		months = mydata$months
	} else if (mydata$cadence == "Daily") {
		cadence = 1
		days = mydata$days
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed() * ireal%%.Machine$integer.max
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generator seed for R
	set.seed(iseed)

    # Take only the first list from opt.list, the second is optional and is needed only for a coupled run

    opt.list = opt.list[[1]]

	nreal  = run.list$nreal
	nMCMC  = run.list$nMCMC
	nlines = run.list$nlines
	ithin = run.list$ithin
	device = run.list$device
	subDir = run.list$subDir


   	nperiods = mydata$nperiods
   	nperiodsFit = mydata$nperiodsFit

   	prior = mydata$prior
   	Temp = mydata$Temp
   	n = mydata$fit$nregions

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    mod_name = mydata$model$name
    fit_name = mydata$fit$name

    tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)

   if (mod_level < fit_level) {

   	tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)
   }

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.numeric(epi.obsrv)
	epi.obsrv = mydata$fit$raw
	tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod
	epi.null.fit = epi.null$cases.ave.fit

	tables.mod$epi.null[nperiodsFit, , mod_id] = epi.null.mod[, mod_id]

	if (mod_level < fit_level) {

		tables.fit$epi.null[nperiodsFit, , fit_id] = epi.null.fit[, fit_id]
	}

	if (tolower(mydata$disease) == "dengue") {
		## alculating correlation/distance with past years and augmenting appropriately

		corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nperiodsFit)

		epi.da.df = get.sql.da(nfit = nperiodsFit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
			corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

		mydata$model$wght = epi.da.df$wght

		mydata$model$epirun = epi.da.df$epi.new

		gamaepi = rep(0, nperiods)
		for (i in 1:nperiods) {
			gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
		}

		mydata$model$gamaepirun = gamaepi

		tps = mydata$ndays

		## Repeat for fit_level


		if (fit_level > mod_level) {


			mydata$fit$epirun = mydata$fit$epi

			mydata$fit$gamaepirun = mydata$fit$gamaepi

			mydata$fit$wght = mydata$fit$epi

			for (iregion in 1:mydata$fit$nregions) {

				corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nperiodsFit)

				distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nperiodsFit)

				epi.da.df = get.sql.da(nfit = nperiodsFit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$fit$epi[, iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

				mydata$fit$wght[, iregion] = epi.da.df$wght

				mydata$fit$epirun[, iregion] = epi.da.df$epi.new

				for (i in 1:length(mydata$fit$epirun[, iregion])) {
					gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
				}
				mydata$fit$gamaepirun[, iregion] = gamaepi

			}
		}

	}

	if (tolower(mydata$disease) == "flu" || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == "ebola") {
    	week0 = mydata$weeks[1]
    	day0 = cadence * (week0 - 1)
	    tps =  seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}

	if (tolower(mydata$disease) == "zika") {
		tps = (1:mydata$nperiods) * cadence
	}

	## In the case of Chikungunya/Zika/Ebola data have to take care of the weights
	## Missing data is given a zero weight
	##
	if (tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == "zika" || tolower(mydata$disease) == "ebola") {
		wght = rep(0, mydata$nperiods)
		wght[1:nperiodsFit] = 1
		mydata$model$wght = wght
		mydata$model$wght[is.na(mydata$model$raw)] = 0
		mydata$fit$wght = mydata$fit$epi
		for (iregion in 1:mydata$fit$nregions) {
			mydata$fit$wght[,iregion] = wght
			mydata$fit$wght[is.na(mydata$fit$raw[,iregion]), iregion] = 0
		}
	}

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever' || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'zika') {
   	setup = setup.vector.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

   } else {

   	setup = setup.model.mcmc(       mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

   }

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    ## Just for the CDC prediction limit R0 nd pC
    # retrieve flu point-of-care units
    flu_poc_units = get_flu_poc_units()
    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {

    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 2.0

    }


    # these are the same for both model and fit

    nparam = setup$nparam
    logbase = setup$logbase
    logvec = setup$logvec

    tab = setup$tab
    nparam = length(model_par)

	## Setting up for Fit level

    if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever' || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'zika') {
    	setup = setup.vector.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    } else {
    	setup = setup.fit.mcmc(       mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    }

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par

    ## Just for the CDC prediction limit pC and R0

   if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
   	fit_pmax["R0", 1:n] = 1.4
   	fit_pmax["pC", 1:n] = 0.2

   }

    ## And let's sample pC and R0 from the prior

    ## We get the mean and sigma of the prior but we actually use it only with prior = 1 or 2 This is controlled by logvec see if
    ## statement below

    setup.prior = setup.single.prior(mydata = mydata, logvec = logvec)

    n1 = length(setup.prior$ymu[, 1])

    fit_ymu = setup.prior$ymu[1:n, ]
    fit_sigma = setup.prior$sigma[1:n, ]
    model_ymu = setup.prior$ymu[n1, ]
    model_sigma = setup.prior$sigma[n1, ]

    if (mydata$fit_level < 4 & mydata$mod_level < 4 & mydata$imodel != 5 & mydata$data_source == "cdc") {

    	prior.ini <- sample.from.prior(mydata = mydata, model_par = model_par, fit_par = fit_par, model_pmin = model_pmin, model_pmax = model_pmax, fit_pmin = fit_pmin, fit_pmax = fit_pmax, model_ymu = model_ymu, model_sigma = model_sigma, fit_ymu = fit_ymu, fit_sigma = fit_sigma)

		model_par = prior.ini$model_par
		fit_par   = prior.ini$fit_par

    }


    if (mydata$prior == 1 || mydata$prior == 2 & tolower(mydata$data_source) == "cdc")
        logvec = setup.prior$logvec


    imask = set.imask(par_names = par_names, opt.list = opt.list)

    # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL
    tab.fit = list()

    pois = 0
    nblock = 3
    accept.vec = array(0, c(n, nblock))


 	### Get the incidence, gamaepi, sh and school

 	if (tolower(mydata$disease) == "dengue") {
		cases = mydata$model$epirun
		gamaepi = mydata$model$gamaepirun
	} else {
		cases = mydata$model$epi
		gamaepi = mydata$model$gamaepi
	}

	sh = mydata$model$sh
	school = mydata$model$school
	wght = mydata$model$wght

    mypar = model_par

    ## Select the prior for this model region
    ymu = model_ymu
    sigma = model_sigma
    ## pad with zeros
    npad = nparam - length(ymu)
    if (npad > 0) {
        ymu = c(ymu, rep(0, npad))
        sigma = c(sigma, rep(1, npad))
        ymu = as.numeric(ymu)
        sigma = as.numeric(sigma)
    }

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

    profile = array(data = 0, dim = c(nRnd, nperiods, n))

	nmydata = length(tps)

	if (tolower(mydata$disease) == "dengue" || tolower(mydata$disease) == 'yellow_fever') {
		time = cumsum(tps)
		ndays = sum(tps) + tps[1] + tps[nmydata]
	} else {
		time = tps[1:nperiods]
		ndays = (nperiods + 2) * cadence
	}

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

	if (is.null(school)) school = rep(0,mydata$nperiods)

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(time[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))

    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile, c(nRnd, nperiods))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

	tables = calc.null.err(tables = tables.mod, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables

	rtn = array(0, c(nperiods, n))
	profile = array(0, c(nRnd, nperiods, n))

    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")
    for (iregion in 1:n) {

        cat("\n Direct Uncoupled Fitting of: ", mydata$fit$name[iregion], "\n\n")

        # Set the seed iseed = set.iseed()

        mypar = fit_par[, iregion]

        if (tolower(mydata$disease) == 'dengue') {
        	cases = mydata$fit$epirun[, iregion]
        	gamaepi = mydata$fit$gamaepirun[, iregion]
        } else {
        	cases = mydata$fit$epi[, iregion]
        	gamaepi = mydata$fit$gamaepi[, iregion]
        }

        sh = mydata$fit$sh[, iregion]
        school = mydata$fit$school[, iregion]
        wght = mydata$fit$wght[, iregion]

        ## Select the prior for this region
        ymu = fit_ymu[iregion, ]
        sigma = fit_sigma[iregion, ]
        ## pad with zeros
        npad = nparam - length(ymu)
        if (npad > 0) {
            ymu = c(ymu, rep(0, npad))
            sigma = c(sigma, rep(1, npad))
        }

		if (is.null(school)) school = rep(0,mydata$nperiods)

		# Also check for NAs
		school[is.na(school)] <- 0
		sh[is.na(sh)]         <- 0

        solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
            wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin[, iregion]), parmax = as.double(fit_pmax[,
                iregion]), dx = as.double(fit_dx[, iregion]), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma),
            scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
            ithin = as.integer(ithin),nperiods = as.integer(nperiods),
            tps = as.double(time[1:nperiods]), rtn = as.double(rtn[, iregion]), accept = as.double(rep(0, nblock)), pois = as.double(pois),
            tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays),profile = 	    as.single(profile[,,iregion]), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))
      
        rtn[, iregion] = solution$rtn

        tab.tmp = matrix(solution$tab, ncol = (nparam + 1))

		# Convert LLK to AICc

		colnames(tab.tmp) = c(myColName, "AICc")
		tab.tmp[, "AICc"] <- 2 * tab.tmp[, "AICc"] + 2 * nopt
		tab.tmp[, "AICc"] <- tab.tmp[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)
		tab.fit[[iregion]] = tab.tmp


        accept.vec[iregion, 1:nblock] = solution$accept

        profile[, , iregion] = array(solution$profile, c(nRnd, nperiods))

		tables = calc.null.err(tables = tables.fit, nfit = nperiodsFit, state_id = fit_id[iregion])
		tables.fit = tables

		tables = calc.mech.err(tables = tables.fit, profiles = profile[,,iregion], nfit = nperiodsFit, state_id = fit_id[iregion])
		tables.fit = tables

    }

    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed\n")

    # # write the MCMC statistics

    success = mcmc.single.write(tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, mydata = mydata,
        imask = imask, ireal = ireal)

     # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

   if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
   	success = writeCSV(mydata = mydata, run.list = run.list, model_rtn = model_rtn, model_profile = model_profile,
   		rtn = rtn, profile = profile, ireal = ireal)
   }

	## Build the aggregate fit to model level mydata. It is a simple sum since these are cases and not %ILI

    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {

        fit_model[i] = sum(rtn[i, 1:n])

        for (irnd in 1:nRnd) {
            for (k in 1:n) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile[irnd, i, k]
        }

	}

	## update table
	tables.agg.mod = tables.mod
	tables = calc.mech.err(tables = tables.agg.mod, profiles = fit_model_profile, nfit = nperiodsFit, state_id = mod_id)
	tables.agg.mod = tables
	##
    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

    if (mydata$data_source == "cdc" || mydata$model$raw_units%in%flu_poc_units) {
    	if (as.character(toupper(run.list$plot)) == "TRUE" || as.character(run.list$plot) == "1") {
    		for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    			mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

    		for (k in 1:length(run.list$device)) success = plotHists(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    			mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)
    	}

    	if (tolower(as.character(run.list$plot)) == "external" || as.character(run.list$plot) == "2") {
    		for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    			mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

    		for (k in 1:length(run.list$device)) success = plotHists.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    			mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

    	}

    } else {

    	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)

    	err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)

    	err <- saveProfiles(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, rtn = rtn,  profile = profile,  ireal = ireal, run.list = run.list)
    }

    err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)

    results = list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)

    return(results)
}

calcRandomProfiles <- function(sh = NULL, school = NULL, pop = NULL, tab = NULL, tps = NULL, nRnd = 100, cadence = 7, epi_model=1) {

    #' Calculate Incidence Disease Profiles
    #'
    #' Given an array with the MCMC history of an uncoupled \code{DICE} run
    #'   randomly choose sets of parameters and use them to create disease profiles.
    #'   The code works on a level of single spatial region and supports both SIR and SEIR models.
    #' @param sh Numeric array with averaged specific humidity for the region
    #' @param school Numeric array with the school vacation schedule for the region
    #' @param pop Numeric - the region's total population
    #' @param tab A numeric array with the MCMC history of all the parameters (both optimized and not)
    #' @param tps A numeric array of days for the flu season - the day numbers correspond to the week numbers
    #' @param nRnd Integer  - the number of random profiles that will be calculated.  (Default is 100)
    #' @param cadence Integer - the cadence of the disease data
    #' @param epi_model Integer 1 = SIR 2 = SEIR
    #' @return profile a 2D array with each row holding a single random disease incidence profile for the region
    #' @examples
    #' calcRandomProfiles{sh = mydata$model$sh, school = mydata$model$school,pop = mydata$model$pop, tab = tab,
    #' tps = tps, nRnd = 100, cadence = 7}

    nperiods = length(sh)
    profile = array(data = 0, dim = c(nRnd, nperiods))

    nlines = dim(tab)[1]
    nparam = dim(tab)[2] - 1
    iburn = nlines/5

    ndays = (nperiods + 2) * cadence

    for (irnd in 1:nRnd) {
        iline = iburn + floor(runif(1, 0, 1) * (nlines - iburn))
        mypar = tab[iline, 1:nparam]

        solution = .Fortran("gensingle", sh = as.double(sh), school = as.double(school), nparam = as.integer(length(mypar)), par = as.double(mypar),
            nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]), rtn = as.double(rep(0, nperiods)), ndays = as.integer(ndays), epi_model=as.integer(epi_model))
        profile[irnd, ] = solution$rtn

    }
    return(profile)

}


calcCPLRandomProfiles <- function(mydata = NULL, tab = NULL, tab.cpl = NULL, Rij = NULL, tps = NULL, nRnd = 100) {

    #' Calculate Coupled Incidence Disease Profiles
    #'
    #' Given an array with the MCMC history of a coupled \code{DICE} run
    #'   randomly choose sets of parameters and use them to create disease incidence profiles for the fit regions.
    #'   The code works on the level of the coupled fit regions.
    #' @param mydata A list with the entire data available for this \code{DICE} run
    #' @param tab The MCMC history of an indirect fit of the model using a coupled model
    #'   This  array includes all the parameters except the two that define the coupling matrix.
    #' @param Rij A 2D array with the Euclidean distances between the centroids of regions \emph{i} and \emph{j}
    #' @param tps A numeric array of days for the disease season
    #' @param nRnd Integer  - the number of random profiles that will be calculated.  (Default is 100)
    #' @return profile a 3D  array with each row holding nRnd profiles for each of the fit regions.
    #' @examples
    #' calcRandomProfiles{mydata = mydata, tab = tab, tab.cpl = tab.cpl, Rij = Rij,
    #' tps = tps, nRnd = 100}



    if (mydata$cadence == "Weekly") {
		cadence = 7
	} else if (mydata$cadence == "Monthly") {
		cadence = 31
	} else if (mydata$cadence == "Daily") {
		cadence = 1
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}


    nperiods = mydata$nperiods
    sh = as.matrix(mydata$fit$sh)
    school = as.matrix(mydata$fit$school)
    pop = mydata$fit$pop
    n = mydata$fit$nregions
    coef = mydata$fit$coef
    epi_model = mydata$epi_model

    profile = array(data = 0, dim = c(nRnd, nperiods, n))
    rtn = array(data = 0, dim = c(nperiods, n))

    nlines = dim(tab)[1]
    nparam = dim(tab)[2] - 1
    nparamCPL = dim(tab.cpl)[2]
    iburn = nlines/5
    mypar = array(0, c(nparam, n))

    ## Need nperiods + 2 because we pathc the weeks at the beginning and end
    ndays = (nperiods + 2) * cadence

    for (irnd in 1:nRnd) {
        iline = iburn + floor(runif(1, 0, 1) * (nlines - iburn))
        iline = min(iline, nlines)

        mypar[1:nparam, 1:n] = tab[iline, 1:nparam, 1:n]
        myparCPL = tab.cpl[iline, ]
        solution = .Fortran("genmulti", n = as.integer(n), sh = as.double(sh), school = as.double(school), nparam = as.integer(nparam),
            par = as.double(mypar), parCPL = as.double(myparCPL), Rij = as.double(Rij), nperiods = as.integer(nperiods),
            tps = as.double(tps[1:nperiods]), rtn = as.double(array(data = 0, dim = c(nperiods, n))), ndays = as.integer(ndays), epi_model=as.integer(epi_model))

        profile[irnd, , ] = matrix(solution$rtn, ncol = n)
        # profile = matrix(solution$rtn, ncol=n)

    }

    return(profile)

}
