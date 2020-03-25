##
## All the functions in this file have been used/worked in the past but they are no longer supported or used by the DICE code
##
##

runDICE2 <- function(dataType = "cdc", year = NULL, mod_level = NULL, mod_name = NULL, fit_level = NULL, fit_name = NULL, nperiodsFit = 52,
    model = 5, isingle = 0, nMCMC = 1e+05, nreal = 1, device = "pdf", prior = 0, Temp = 1, subDir = NULL, plot = 1, movie = 0, emcee = FALSE,
    nwalk = 1, iseed = NULL, Tg = NULL) {


    # The main driver for the DICE package
    #
    # The main driver for dice which grabs the mydata for the requested mydata set, season and model/fit spatial regions combinations. After mydata is retrived, the simulation is setup and the function calls either the single or multi options - depending on the user's request for an uncoupled or coupled run. In either case the fit begins with an MCMC procedure on the model-level mydata.
    # @param mod_level Integer - Spatial level of the model data. Spatial levels for the United States are defined as follows:
    # \tabular{rl}{
    # 2 \tab Country \cr
    # 3 \tab HHS Region \cr
    # 4 \tab State \cr
    # 5 \tab County \cr
    # 6 \tab City \cr}
    # @param mod_name \itemize{
    #   \item{Character Vector that describes the model population.  For example \code{mod_name=c(NAME_2='USA',NAME_3='Region2',NAME_4='NY')} combined with \code{mod_level=4} specifies the state of New York.  \code{NAME_X} may be the fullname or abbreviation of the \code{X}-level population.  For instance \code{mod_name=c(NAME_2='United States',NAME_3='R2',NAME_4='New York')} can also be used to specify New York.}
    #   \item{For \code{mod_name='create'} the populations specified by fit_level/fit_name will be used to create a model population.  Currently dataTypes 'CDC', 'GFT', and 'MiscILI' are supported for the 'create' functionality.}}
    # @param fit_level If \code{mod_name='create'}, then \code{fit_level} is a vector of integers indicating the spatial level of each region to be used. Otherwise \code{fit_level} is a single integer.
    # @param fit_name \itemize{
    #   \item{If \code{mod_name='create'}, then \code{fit_name} is a list of character vectors: \code{fit_name = list(c(NAME_2='USA',NAME_3='Region10'),c(NAME_2='USA',NAME_3='Region9',NAME_4='AZ'),c(NAME_2='USA',NAME_3='Region9',NAME_4='CA'))}. This example combined with \code{mod_name='create'} and \code{fit_level=c(3,4,4)} creates a super-region from HHS Region 10, California, and Arizona.}
    #   \item{When mod_name/mod_level is not in 'create' mode, the only option currently supported is \code{fit_name='all'}. This will pull all populations at fit_level that are subpopulations of the mod_name population.}}
    # @param dataType Select from mydata-types: 'CDC', 'GFT', 'MiscILI', 'MiscCases', 'ZikaS', or 'ZikaC' (NOT case sensitive).  These indicate the source and type of mydata: Centers for Disease Control (CDC) mydata is \%ILI, Google Flu Trends (GFT) is \%ILI, 'MiscILI' is \%ILI mydata from miscelaneous sources, 'MiscCases' is absolute-cases mydata from miscelaneous sources, 'ZikaS' is Zika suspected cases (PAHO), and 'ZikaC' is Zika confirmed cases (PAHO). \itemize{
    # \item{If \code{mod_name='create'}, then dataType is a vector of strings indicating the mydata-type to be used for each fit_name population. \code{dataType = c('CDC','GFT','MiscCases')}}
    # \item{Otherwise dataType is a single string that determines the mydata-type for all model/fit populations.}}
    # @param year Integer - start year of the flu season
    # @param nperiodsFit Integer - Number of weeks that will be fitted.  Default is to fit all the mydata.  This will be reset if nperiodsFit > nperiodsData
    # @param model Integer - The model number, see manual for more details (1-5 are supported)
    # @param isingle Integer - 0: couple the fit spatial regions; 1: no coupling
    # @param nMCMC Integer - number of steps/trials in the MCMC procedure
    # @param nreal Integer - number of chains
    # @param device  Either 'pdf' (default) or 'x11'
    # @param plot TRUE, FALSE or EXTERNAL (or 0, 1, 2) allows the Users to implement their own plotting routines
    # @param movie Integer 1/ 0 save or not png frames that can be stringed to a movie. Default is NO = 0 (movie frame generation is only supported for USA CDC/HHS Regions)
    # @param subDir Name of output sub-directory where all plots and files will be written.  Default is 'output'
    # @param prior Integer - if greater than zero use a prior for the MCMC procedure
    # @param Temp Integer 1, 2, ..10 Relevant only if using a prior - controls the width of the Gaussian prior
    # @param emcee Logical - relevant only when running a sinlge region/patch slects between MCMC and EMCEE procedures
    # @param nwalk Integer - Number of walkers, relevant ONLY for an EMCEE procedure
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @param Tg - recovery time in days, reasonable range is 2-3 days
    # @return solution a list with the input and entire output of the run.
    # @examples
    # For a run of the 2015-2016 cdc national mydata using the ten HHS regions with coupling between the regions use:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0)
    #
    # For a run of the 2015-2016 cdc national mydata using the ten HHS regions without coupling between the regions use:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 1)
    #
    # For a run of the 2014-2015 GFT mydata for HHS region number 9, using state level mydata with coupling between the states in region 9 use:
    # output <- runDICE2(dataType='gft', year = 2014, mod_level = 3, mod_name=c(NAME_2='USA',NAME_3='Region9'), fit_level = 4, fit_name='all', isingle = 0)
    #
    # To control which model is used for the basic reproduction number, set the parameter model in your call. Default value is 5:
    # output <- runDICE2(dataType='gft', year = 2014, mod_level = 3, mod_name=c(NAME_2='USA',NAME_3='Region9'), fit_level = 4, fit_name='all', isingle = 0, model = 3)
    #
    # To control the number of MCMC chains that the code will run set the parameter  nreal in your call, default is 1:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0, nreal = 3)
    #
    # To control the number of MCMC steps/trial in each chain set the parameter nMCMC in your call, default is 1e5:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0, nMCMC = 1e6)
    #
    # To control the name of the sub-directory where all the output files and plots are saved use the keyword subDir, default is output:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0, nMCMC = 1e6, subDir = 'test')
   #
    # To control the file format for the plots (pdf, png or x11) set the parameter device:
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0, nMCMC = 1e6, device = 'pdf')
    # (The package can accept an array of file formats, i.e. device = c('pdf','png'), in which case more both 'png' and 'pdf' files will be created.)
    #
    # By default png frames are created for each epidemic week and they can be used to create a movie. The keyword movie = 1 (yes), or 0 (no) controls this
    # output <- runDICE2(dataType='cdc', year = 2015, mod_level = 2, mod_name=c(NAME_2='USA'), fit_level = 3, fit_name='all', isingle = 0, nMCMC = 1e6, movie=1)
    # To run in a forecast or predictive mode you can set the number of weeks the code uses in the fit to be lower than the number of weeks in the season.
    # (Note that for the current season it is always running in a predictive mode because the season is not yet completed.)
    # output <- runDICE2(dataType='gft', year = 2013, mod_level = 3, mod_name=c(NAME_2='USA',NAME_3='Region9'), fit_level = 4, fit_name='all', isingle = 1, nMCMC = 1e6, nperiodsFit = 35)
    #
    # To select only a few HHS regions and run them coupled (for example the Eastern Regions 1, 2 and 3) use:
    # output <- runDICE2(dataType=c('cdc','cdc','cdc'), year=2015, mod_level=2, mod_name='create', fit_level=c(3,3,3), fit_name=list(c(NAME_2='USA',NAME_3='Region1'),c(NAME_2='USA',NAME_3='Region2'),c(NAME_2='USA',NAME_3='R3')), isingle = 0)
    #
    #To select only a few states and run them coupled  use for example:
    # output <- runDICE2(dataType=c('gft','gft','gft'), year=2014, mod_level=3, mod_name='create', fit_level=c(4,4,4), fit_name=list(c(NAME_2='USA',NAME_3='Region10',NAME_4='WA'),c(NAME_2='USA',NAME_3='Region10',NAME_4='OR'),c(NAME_2='USA',NAME_3='R9',NAME_4='CA')), isingle = 0)

    ## Data Type: cdc or gft

    if (is.null(dataType))
        dataType = "cdc"

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

    if (is.null(year))
        year = 2015

    ## Model Number
    if (is.null(model))
        model = 5

    ## Number of weeks of mydata that are fitted
    if (is.null(nperiodsFit))
        nperiodsFit = 52

    ## Recovery time in days
    if (is.null(Tg))
        Tg = 3

    ## Coupled (isingle = 0) or Uncoupled (isingle = 1) fit for the spatial regions

    if (is.null(isingle))
        isingle = 0

    ## device name for plotting the results - pdf or X11
    if (is.null(device))
        device = "pdf"

    ## Use (prior > 0) a prior in the MCMC procedure, or not (prior = 0)
    if (is.null(prior))
        prior = 0

    ## 'Temperature' for the prior width. Used only if prior > 0

    if (is.null(Temp))
        Temp = 1

    if (is.null(subDir))
        subDir = "output"

    if (is.null(plot))
        plot = 1
    if (is.null(movie))
        movie = 1

    if (is.null(emcee))
        emcee = FALSE
    if (is.null(nwalk))
        nwalk = 1

    # this will get us the mydata at the region level we requested for the years and weeks we asked for

    ## This get's state level mydata For the SH it has both the state and region level mydata

    # Set start week and end week
    if (year == 2009) {
        week.start = 13
        week.end = 26
    } else {
        # This is the standard CDC definition of start/end of EPI year
        week.start = 27
        week.end = 26
    }

    # Recover ili,sh,school,pop mydata at both model and fit levels
    mydata = get.subset2(mod_level = mod_level, mod_name = mod_name, fit_level = fit_level, fit_name = fit_name, start.year = year, start.week = week.start,
        end.year = year + 1, end.week = week.end, dataType = dataType)

    if (nperiodsFit <= 0)
        nperiodsFit = mydata$nperiodsData
    if (nperiodsFit > mydata$nperiodsData)
        nperiodsFit = mydata$nperiodsData
    if (nperiodsFit < 10) {
        cat("\n\n WARNING: Resetting nperiodsFit to its minimal value of 10 \n\n")
        nperiodsFit = 10
    }
    if (nperiodsFit < 15) {
        cat("\n\n WARNING: Using a small Number of weeks in fitting procedure
        \n
        Results will have limited information \n\n")
    }
    mydata$nperiodsFit = nperiodsFit

    mydata$imodel = model
    if (length(unique(fit_level)) == 1) {
        mydata$fitLevelName = get.FitLevelName(iFit = fit_level[1])
    } else {
        mydata$fitLevelName = "mixed"
    }

    # mydata$modelLevelName = get.ModelLevelName(mod_level, RegState)
    mydata$modelLevelName = mydata$model$name

    if (length(unique(tolower(dataType))) == 1) {
        dataType = tolower(dataType)
        mydata$dataType = dataType

        if (isingle == 0)
            mydata$dataName = paste(dataType, "-", mydata$model$name, "-cpl-", mydata$FY, "-", mydata$imodel, sep = "")
        if (isingle == 1)
            mydata$dataName = paste(dataType, "-", mydata$model$name, "-uncpl-", mydata$FY, "-", mydata$imodel, sep = "")
        if (all(fit_level == mod_level))
            mydata$dataName = paste(dataType, "-", mydata$model$name, "-", mydata$FY, "-", mydata$imodel, sep = "")
    } else {
        dataType_temp = "mixed"
        mydata$dataType = dataType_temp

        if (isingle == 0)
            mydata$dataName = paste(dataType_temp, "-", mydata$model$name, "-cpl-", mydata$FY, "-", mydata$imodel, sep = "")
        if (isingle == 1)
            mydata$dataName = paste(dataType_temp, "-", mydata$model$name, "-uncpl-", mydata$FY, "-", mydata$imodel, sep = "")
        if (all(fit_level == mod_level))
            mydata$dataName = paste(dataType_temp, "-", mydata$model$name, "-", mydata$FY, "-", mydata$imodel, sep = "")
    }

    mydata$FY = paste(year, "-", (year + 1), sep = "")
    mydata$season = year

    ## augment the mydata with some more information that the user has chosen

    mydata$single = isingle
    ## TRUE or FALSE - default is FALSE

    mydata$prior = prior
    mydata$Temp = Temp

    ############ get.subset2() returns NAs where mydata is unavailable. For now, we replace them with 0s to avoid errors in the MCMC procedure
    if (any(is.na(mydata$fit$school))) {
        cat("\n\n WARNING: NAs in school fit_level mydata are being set to 0. \n\n")
        mydata$fit$school[is.na(mydata$fit$school)] = 0
    }
    if (any(is.na(mydata$fit$sh))) {
        cat("\n\n WARNING: NAs in specific humidity fit_level mydata are being set to 0. \n\n")
        mydata$fit$sh[is.na(mydata$fit$sh)] = 0
    }
    if (any(is.na(mydata$fit$raw))) {
        cat("\n\n WARNING: NAs in epi fit_level mydata are being set to 0. \n\n")
        mydata$fit$epi[is.na(mydata$fit$epi)] = 0
        mydata$fit$gamaepi[is.na(mydata$fit$gamaepi)] = 0
        mydata$fit$raw[is.na(mydata$fit$raw)] = 0
    }

    if (any(is.na(mydata$model$school))) {
        cat("\n\n WARNING: NAs in school model_level mydata are being set to 0. \n\n")
        mydata$model$school[is.na(mydata$model$school)] = 0
    }
    if (any(is.na(mydata$model$sh))) {
        cat("\n\n WARNING: NAs in specific humidity model_level mydata are being set to 0. \n\n")
        mydata$model$sh[is.na(mydata$model$sh)] = 0
    }
    if (any(is.na(mydata$model$raw))) {
        cat("\n\n WARNING: NAs in epi model_level mydata are being set to 0. \n\n")
        mydata$model$epi[is.na(mydata$model$epi)] = 0
        mydata$model$gamaepi[is.na(mydata$model$gamaepi)] = 0
        mydata$model$raw[is.na(mydata$model$raw)] = 0
    }

    ## Pack the information for the run
    opt.list <- set.opt.list(model = mydata$imodel, isingle = isingle)

    run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, nlines = nlines, device = device, subDir = subDir, plot = plot, movie = movie,
        nwalk = nwalk)

    ## Fitting a SINGLE region/patch
    if (mydata$fit$nregions == 1) {

        cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
        cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$dataType), " mydata", "\n")
        cat("\n\n Fitting ", mydata$nperiodsFit, " weeks out of ", mydata$nperiodsData, " weeks of mydata", "\n\n")
        cat("\n\n Fitting ", toupper(mydata$modelLevelName), " mydata", "\n")

        for (ireal in 1:nreal) {
            solution = fitOnePatch(mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
        }

        return(solution)
    }

    ## ALL OTHER CASES - number of regions/patches of fit_level > 1

    cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
    cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$dataType), " mydata", "\n")
    cat("\n\n Fitting ", mydata$nperiodsFit, " weeks out of ", mydata$nperiodsData, " weeks of mydata", "\n\n")
    cat("\n\n Fitting ", toupper(mydata$modelLevelName), " using ", toupper(mydata$fitLevelName), " mydata", "\n")

    ## UNCOUPLED CASAE

    if (isingle == TRUE) {
        for (ireal in 1:nreal) {

            solution = fitSingle(mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)

        }
        ## COUPLED CASE

    } else {
        for (ireal in 1:nreal) {
            solution = fitMulti(mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
        }
    }

    return(solution)
}


get.subset <- function(mod_level=2, fit_level=3, mod_name=c(NAME_2 = "USA"), start.year=2016, start.week=27, end.year=2017, end.week=26, data_source="cdc") {

  # Load the full \pkg{dice} dataset and then subset by date and region
  #
  # This function is generally called by \code{\link{get.DICE.data}}.
  # @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has mydata at levels 2-3 for CDC and 2-4 for GFT.
  # @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  # @param mod_name A named-vector of character strings that specify which region is to be \strong{model}ed.  In other words, \code{mod_name} specifies the country, region, state, etc. of the mod_level region.  \code{mod_name} should be of the form \code{mod_name = c(NAME_2='a', NAME_3='b',..., NAME_i='x'} where i=mod_level and 'a', 'b',...,'x' are the appropriate level names. NAME_i='x' also accepts abbreviations.  Choose appropriate names from \code{\link{diceData}}.  For example, mod_name=c(NAME_2='United.States',NAME_3='Region4',NAME_4='North.Carolina') and mod_level=4 specifies North Carolina. To achieve the same result, use all abbreviations mod_name=c(NAME_2='USA',NAME_3='R4',NAME_4='NC') or a mix of names and abbreviations mod_name=c(NAME_2='USA',NAME_3='Region4',NAME_4='NC')
  # @param start.year An integer - start year of the flu season
  # @param end.year  An integer  - end year of the flu season (default is end.year = start.year + 1)
  # @param start.week An integer - starting CDC week of the flu season (default is 27)
  # @param end.week   An integer - ending CDC week of the flu season (default is 26)
  # @param data_source A string specifying the type of mydata the user wants to model, currently we support `cdc' and `gft' (default is cdc)
  # @return mydata  A list ILI, SH, School, and Census mydata for both the model- and fit-level region(s).
  # @examples
  # require(DICE)
  # mydata = get.subset(mod_level = 3, fit_level = 4, mod_name = c(NAME_2='USA',NAME_3='Region4'), start.year = 2014, data_source = 'gft')

  #require(countrycode) # this should be removed once mod_name["NAME_2"] is converted to ISO2 package has been added to DESCRIPTION file

  cat("Warning: The function get.subset() has been deprecated and will not be supported going forward.  Please use get.cdc.data() or get.mysql() instead. \n")

  # This function loads the full DICE datasetand then subsets by date and region
  mydata = list()
  # Current year. Used to determine if mydata should be downloaded from CDC curYear = 2016
  today = Sys.Date()
  CDCweek = Date2CDCweek(today)
  if (CDCweek$CDCweek > 26) {
    curSeason = CDCweek$CDCyear
  } else {
    curSeason = CDCweek$CDCyear - 1
  }

  # End year. Evaluate which season end.year/end.week falls in
  if (end.week>26) {
    endSeason = end.year
  } else {
    endSeason = end.year - 1
  }

  # load the dataset(later this will be a package-compatable call OR the subsetting may be done at the SQL level) no need for this
  # because LazyData was set to yes in the DESCRIPTION file

  # mydata(DICE_dataset, package='DICE') returns a mydataset/list called 'diceData'

  # Determine index of mod_level and fit_level patch(es) in 'attr'
  if (mod_level >= 2) {
    mod_ind = (toupper(mod_name["NAME_2"]) == diceData$attr$ABBV_2 | mod_name["NAME_2"] == diceData$attr$NAME_2) & mod_level == diceData$attr$level
    fit_ind = (toupper(mod_name["NAME_2"]) == diceData$attr$ABBV_2 | mod_name["NAME_2"] == diceData$attr$NAME_2) & fit_level == diceData$attr$level
    if (mod_level > 2) {
      for (lev in 3:mod_level) {
        mod_ind = mod_ind & (mod_name[paste0("NAME_", lev)] == diceData$attr[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_",
                                                                                                                               lev)]) == diceData$attr[, paste0("ABBV_", lev)])
        fit_ind = fit_ind & (mod_name[paste0("NAME_", lev)] == diceData$attr[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_",
                                                                                                                               lev)]) == diceData$attr[, paste0("ABBV_", lev)])
      }
    }
  } else {
    # Global/Continent level mydata
  }
  ModObj = which(mod_ind)
  FitObj = which(fit_ind)

  # Check that model-patch exists in mydata set
  if (length(ModObj) == 0) {
    stop("model name/level combination not present in mydata.")
  }
  if (length(FitObj) == 0) {
    stop("fit-level mydata is not present for model name/level.")
  }

  # Create a mod_level identifier
  if (mod_level >= 2) {
    mod_ident = diceData$attr$ABBV_2[ModObj]
    if (mod_level > 2) {
      for (lev in 3:mod_level) {
        mod_ident = paste0(mod_ident, ".", diceData$attr[ModObj, paste0("ABBV_", lev)])
      }
    }
  } else {
    # Global/Continent level mydata
  }

  # Create a vector of fit_level identifiers
  fit_ident = character(length(FitObj))
  for (ii in 1:length(FitObj)) {
    obj = FitObj[ii]
    if (fit_level >= 2) {
      fit_ident[ii] = diceData$attr$ABBV_2[obj]
      if (fit_level > 2) {
        for (lev in 3:fit_level) {
          fit_ident[ii] = paste0(fit_ident[ii], ".", diceData$attr[obj, paste0("ABBV_", lev)])
        }
      }
    } else {
      # Global/Continent level mydata
    }
  }

  # Check that model and fit data exist in the specified data_source
  if (tolower(data_source) == "cdc") {
    dataNames = colnames(diceData$CDCili)
  } else if (tolower(data_source) == "gft") {
    dataNames = colnames(diceData$GFTili)
  } else if (tolower(data_source) == "miscili") {
    dataNames = colnames(diceData$MiscIli)
  } else if (tolower(data_source) == "misccases") {
    dataNames = colnames(diceData$MiscCases)
  }

  if (!any(dataNames == mod_ident)) {
    stop("DICE does not have the specified mod_level mydata for mydata-type: ", data_source)
  }
  if (!all(fit_ident %in% dataNames)) {
    stop("DICE does not have the specified fit_level mydata for mydata-type: ", data_source)
  }

  # Check that non-NA ILI mydata exists for the date range


  # Determine what time indices are included
  start.ind = which(diceData$school$year == start.year & diceData$school$week == start.week)
  end.ind = which(diceData$school$year == end.year & diceData$school$week == end.week)

  mydata$dates = diceData$school$date[start.ind:end.ind]
  mydata$weeks = diceData$school$week[start.ind:end.ind]
  mydata$years = diceData$school$year[start.ind:end.ind]
  mydata$nperiods = length(mydata$weeks)


  # Pull attr and school mydata for fit_level
  mydata$fit = list()
  mydata$fit$level = fit_level
  mydata$fit$name = diceData$attr[fit_ind, paste0("NAME_", fit_level)]
  mydata$fit$attr = diceData$attr[fit_ind, ]
  mydata$fit$school = diceData$school[start.ind:end.ind, fit_ident]
  mydata$fit$school[is.na(mydata$fit$school)] = 0
  mydata$fit$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, fit_ident])
  mydata$fit$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, fit_ident]
  mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)


  # Repeat for model data
  mydata$model = list()
  mydata$model$level = mod_level
  mydata$model$name = diceData$attr[mod_ind, paste0("NAME_", mod_level)]
  mydata$model$attr = diceData$attr[mod_ind, ]
  mydata$model$school = diceData$school[start.ind:end.ind, mod_ident]
  mydata$model$school[is.na(mydata$model$school)] = 0
  mydata$model$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, mod_ident])
  mydata$model$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, mod_ident]
  mydata$model$coef = mydata$model$pop/sum(mydata$model$pop)

  # Get ili mydata
  if (tolower(data_source) == "gft") {
    start.ind = which(diceData$GFTili$year == start.year & diceData$GFTili$week == start.week)
    end.ind = which(diceData$GFTili$year == end.year & diceData$GFTili$week == end.week)
    mydata$fit$raw = diceData$GFTili[start.ind:end.ind, fit_ident]
    mydata$model$raw = diceData$GFTili[start.ind:end.ind, mod_ident]

    mydata$model$raw_units = diceData$MetaILI$GFTili_raw_units[mod_ident]
    mydata$model$epi_units = diceData$MetaILI$GFTili_epi_units[mod_ident]
    mydata$model$factor = diceData$MetaILI$GFTili_factor[diceData$MetaILI$GFTili_factor$year == start.year, mod_ident]

    mydata$fit$raw_units = diceData$MetaILI$GFTili_raw_units[fit_ident]
    mydata$fit$epi_units = diceData$MetaILI$GFTili_epi_units[fit_ident]
    mydata$fit$factor = diceData$MetaILI$GFTili_factor[diceData$MetaILI$GFTili_factor$year == start.year, fit_ident]

  } else if (tolower(data_source == "miscili")) {
    start.ind = which(diceData$MiscIli$year == start.year & diceData$MiscIli$week == start.week)
    end.ind = which(diceData$MiscIli$year == end.year & diceData$MiscIli$week == end.week)
    mydata$fit$raw = diceData$MiscIli[start.ind:end.ind, fit_ident]
    mydata$model$raw = diceData$MiscIli[start.ind:end.ind, mod_ident]

    mydata$model$raw_units = diceData$MetaILI$MiscIli_raw_units[mod_ident]
    mydata$model$epi_units = diceData$MetaILI$MiscIli_epi_units[mod_ident]
    mydata$model$factor = diceData$MetaILI$MiscIli_factor[diceData$MetaILI$MiscIli_factor$year == start.year, mod_ident]

    mydata$fit$raw_units = diceData$MetaILI$MiscIli_raw_units[fit_ident]
    mydata$fit$epi_units = diceData$MetaILI$MiscIli_epi_units[fit_ident]
    mydata$fit$factor = diceData$MetaILI$MiscIli_factor[diceData$MetaILI$MiscIli_factor$year == start.year, fit_ident]

  } else if (tolower(data_source == "misccases")) {
    start.ind = which(diceData$MiscCases$year == start.year & diceData$MiscCases$week == start.week)
    end.ind = which(diceData$MiscCases$year == end.year & diceData$MiscCases$week == end.week)
    mydata$fit$raw = diceData$MiscCases[start.ind:end.ind, fit_ident]
    mydata$model$raw = diceData$MiscCases[start.ind:end.ind, mod_ident]

    mydata$model$raw_units = diceData$MetaILI$MiscCases_raw_units[mod_ident]
    mydata$model$epi_units = diceData$MetaILI$MiscCases_epi_units[mod_ident]
    mydata$model$factor = diceData$MetaILI$MiscCases_factor[diceData$MetaILI$MiscCases_factor$year == start.year, mod_ident]

    mydata$fit$raw_units = diceData$MetaILI$MiscCases_raw_units[fit_ident]
    mydata$fit$epi_units = diceData$MetaILI$MiscCases_epi_units[fit_ident]
    mydata$fit$factor = diceData$MetaILI$MiscCases_factor[diceData$MetaILI$MiscCases_factor$year == start.year, fit_ident]

    # 'cases' type mydata should not use the CDC-onset values
    mydata$fit$onset = numeric()
    mydata$model$onset = numeric()

  } else if (tolower(data_source) == "cdc") {
    GoodPull_mod = FALSE
    GoodPull_fit = FALSE
    if (endSeason >= curSeason)
    {
      # Attempt to download 'model' ILI mydata from CDC server
      if (mod_level == 2) {
        CDCreg = "national"
      } else if (mod_level == 3) {
        CDCreg = "hhs"
        sub_region = mydata$model$attr$ID_3
      } else if (mod_level == 4) {
        CDCreg = "state"
        sub_region = gsub(".", " ", mydata$model$attr$NAME_4, fixed = TRUE)
      }

      # Pull ili mydata from CDC server (make 5 attempts)
      n = 0
      while (n < 5) {
        CDCmodel = try(ilinet(region = CDCreg, years = (start.year-1):(end.year-1)), silent = TRUE)
        n = n + 1
        if (!is(CDCmodel, "try-error")) {
          GoodPull_mod = TRUE
          if (mod_level == 3) {
            # reduce mydata to specified regions
            reg_ind = CDCmodel$region %in% paste0("Region ", sub_region)
            CDCmodel = CDCmodel[reg_ind, ]
          } else if (mod_level == 4) {
            state_ind = CDCmodel$region %in% sub_region
            CDCmodel = CDCmodel[state_ind, ]
          }
          break
        } else Sys.sleep(5)
      }
      # If mydata download is unsuccessful, print warning and revert to saved mydata.
      if (!GoodPull_mod) {
        print("WARNING: Download from CDC server was unsuccessful for model-level mydata. Reverting to ili mydata in local DICE-package mydatabase.  DICE package mydata is likely less up-to-date than the CDC server.")
      }

      # Attempt to download 'fit' ILI mydata from CDC server
      if (fit_level == 2) {
        CDCreg = "national"
      } else if (fit_level == 3) {
        CDCreg = "hhs"
        sub_region = mydata$fit$attr$ID_3
      } else if (fit_level == 4) {
        CDCreg = "state"
      }
      # Pull ili mydata from CDC server (make 5 attempts)
      n = 0
      while (n < 5) {
        CDCfit = try(ilinet(region = CDCreg, years = (start.year-1):(end.year-1)), silent = TRUE)
        n = n + 1
        if (!is(CDCfit, "try-error")) {
          GoodPull_fit = TRUE
          if (fit_level == 3) {
            # reduce mydata to specified regions
            reg_ind = CDCfit$region %in% paste0("Region ", sub_region)
            CDCfit = CDCfit[reg_ind, ]
          }
          break
        } else Sys.sleep(5)
      }
      # If mydata download is unsuccessful, print warning and revert to saved mydata.
      if (!GoodPull_fit) {
        print("WARNING: Download from CDC server was unsuccessful for fit-level mydata. Reverting to ili mydata in local DICE-package mydatabase.  DICE package mydata is likely less up-to-date than the CDC server.")
      }
    }  # End if CurYear

    if (!GoodPull_mod) {
      # Pull mydata from 'diceData'
      start.ind = which(diceData$CDCili$year == start.year & diceData$CDCili$week == start.week)
      end.ind = which(diceData$CDCili$year == end.year & diceData$CDCili$week == end.week)
      # if the datasetdoes not contain the full year, pull the weeks that are available
      if (length(end.ind) == 0) {
        end.ind = length(diceData$CDCili$week)
      }
      mydata$model$raw = diceData$CDCili[start.ind:end.ind, mod_ident]

    } else {
      # Process CDC download
      start.ind = which(CDCmodel$year == start.year & CDCmodel$week == start.week)
      end.ind = which(CDCmodel$year == end.year & CDCmodel$week == end.week)
      if (length(end.ind) == 0) {
        end.ind = length(CDCmodel$year)
      }
      if (mod_level == 4) {
        mydata$model$raw = CDCmodel$unweighted_ili[start.ind:end.ind]
      } else {
        mydata$model$raw = CDCmodel$weighted_ili[start.ind:end.ind]
      }
    }  # end if GoodPull_mod

    # pad with zeros to fill-out a complete year
    padweeks = mydata$nperiods - length(start.ind:end.ind)
    mydata$model$raw = c(mydata$model$raw, rep(0, padweeks))

    # Fill-in MetaILI mydata for 'model'
    mydata$model$raw_units = diceData$MetaILI$CDCili_raw_units[mod_ident]
    mydata$model$epi_units = diceData$MetaILI$CDCili_epi_units[mod_ident]
    mydata$model$factor = diceData$MetaILI$CDCili_factor[diceData$MetaILI$CDCili_factor$year == start.year, mod_ident]

    if (!GoodPull_fit) {
      # Pull fit-mydata from DICE mydataset
      start.ind = which(diceData$CDCili$year == start.year & diceData$CDCili$week == start.week)
      end.ind = which(diceData$CDCili$year == end.year & diceData$CDCili$week == end.week)
      # if the datasetdoes not contain the full year, pull the weeks that are available
      if (length(end.ind) == 0) {
        end.ind = length(diceData$CDCili$week)
      }
      mydata$fit$raw = diceData$CDCili[start.ind:end.ind, fit_ident]
    } else {
      if (fit_level == 4) {
        # sort out only the states requested
        req_names = gsub(".", " ", mydata$fit$attr$NAME_4, fixed = TRUE)
        keep_ind = CDCfit$region %in% req_names
        CDCfit = CDCfit[keep_ind, ]
      }
      nregions = length(unique(CDCfit$region))
      # Process CDC 'fit' download
      if (nregions > 1) {
        tempYear = CDCfit$year[seq(from = 1, to = length(CDCfit$year), by = nregions)]
        tempWeek = CDCfit$week[seq(from = 1, to = length(CDCfit$year), by = nregions)]
        start.ind = which(tempYear == start.year & tempWeek == start.week)
        end.ind = which(tempYear == end.year & tempWeek == end.week)
        if (length(end.ind) == 0) {
          end.ind = length(tempYear)
        }
        # reshape CDC mydata to DICE-matrix form
        if (fit_level == 4) {
          FitILI = array(data = CDCfit$unweighted_ili, dim = c(nregions, length(CDCfit$year)/nregions), dimnames = list(unique(CDCfit$region),
                                                                                                                        NULL))[, start.ind:end.ind]
          # Then re-order to match mydata$fit$attr
          FitILI = FitILI[gsub(".", " ", mydata$fit$attr$NAME_4, fixed = TRUE), ]
        } else {
          FitILI = array(data = CDCfit$weighted_ili, dim = c(nregions, length(CDCfit$year)/nregions))[, start.ind:end.ind]
        }
        mydata$fit$raw = t(FitILI)
      } else {
        # only one fit region
        start.ind = which(CDCfit$year == start.year & CDCfit$week == start.week)
        end.ind = which(CDCfit$year == end.year & CDCfit$week == end.week)
        if (length(end.ind) == 0) {
          end.ind = length(CDCfit$year)
        }
        mydata$fit$raw = CDCfit$weighted_ili[start.ind:end.ind]
      }

      # If working with state-level, check for columns with NAs
      if (fit_level == 4)
      {
        NA_ind = apply(mydata$fit$raw, 2, function(x) any(is.na(x)))
        # if any states are missing mydata, attempt to reconstruct their ILI from region totals
        if (any(NA_ind))
        {
          regions = unique(mydata$fit$attr$ID_3[NA_ind])

          n = 0
          while (n < 5) {
            CDCregs = try(ilinet(region = "hhs", years = (start.year-1):(end.year-1)), silent = TRUE)
            n = n + 1
            if (!is(CDCregs, "try-error")) {
              GoodPull_regs = TRUE
              # reduce mydata to specified regions
              reg_ind = CDCregs$region %in% paste0("Region ", regions)
              CDCregs = CDCregs[reg_ind, ]
              break
            } else Sys.sleep(5)
          }

          # grab state ILITOTAL and TOTAL PATIENTS
          StateILITOTAL = array(data = CDCfit$ilitotal, dim = c(nregions, length(CDCfit$year)/nregions), dimnames = list(unique(CDCfit$region), NULL))[, start.ind:end.ind]
          StateILITOTAL = t(StateILITOTAL[gsub(".", " ", mydata$fit$attr$NAME_4, fixed = TRUE), ])
          StateTotal = array(data = CDCfit$total_patients, dim = c(nregions, length(CDCfit$year)/nregions), dimnames = list(unique(CDCfit$region), NULL))[, start.ind:end.ind]
          StateTotal = t(StateTotal[gsub(".", " ", mydata$fit$attr$NAME_4, fixed = TRUE), ])
          # reshape region mydata to a 'raw' type array
          nregions = length(regions)
          if (nregions == 1) {
            RegILITOTAL = array(data = CDCregs$ilitotal[start.ind:end.ind], dim = c(end.ind - start.ind + 1, nregions))
            RegTotal = array(data = CDCregs$total_patients[start.ind:end.ind], dim = c(end.ind - start.ind + 1, nregions))
          } else {
            RegILITOTAL = t(array(data = CDCregs$ilitotal, dim = c(nregions, length(CDCregs$year)/nregions))[, start.ind:end.ind])
            RegTotal = t(array(data = CDCregs$total_patients, dim = c(nregions, length(CDCregs$year)/nregions))[, start.ind:end.ind])
          }


          # go through each week and attempt to reconstruct state ILI
          for (week in 1:nrow(mydata$fit$raw)) {
            NA_ind = is.na(mydata$fit$raw[week, ])
            if (any(NA_ind)) {
              NA_reg = mydata$fit$attr$ID_3[NA_ind]
              # if more than one state is missing from a region, do not attempt to reconstruct
              dup_ind = duplicated(NA_reg)
              if (any(dup_ind)) {
                dup_reg = NA_reg[dup_ind]
                # cat('Warning!!!! Region(s) ',dup_reg,' are missing mydata from more that one state.  Individual state ILI cannot be reconstructed
                # from region totals.')
                remove_ind = NA_reg %in% dup_reg
                NA_reg = NA_reg[!remove_ind]
                NA_ind = NA_ind & mydata$fit$attr$ID_3 %in% NA_reg
              }
            }

            if (any(NA_ind)) {
              # use region totals to reconstruct state ILI
              for (state_num in which(NA_ind)) {
                region = mydata$fit$attr$ID_3[state_num]
                state_tot = RegTotal[week, regions == region] - sum(StateTotal[week, mydata$fit$attr$ID_3 == region], na.rm = TRUE)
                state_ilitot = RegILITOTAL[week, regions == region] - sum(StateILITOTAL[week, mydata$fit$attr$ID_3 == region],
                                                                          na.rm = TRUE)
                mydata$fit$raw[week, state_num] = state_ilitot/state_tot * 100
              }
            }
          }
        }  # end if any(state NA)
      }  # end if state NA fill
    }  # end if GoodPull_fit

    # padweeks if necessary
    padweeks = mydata$nperiods - length(start.ind:end.ind)
    if (length(fit_ident) > 1) {
      # pad with zeros to fill-out a complete year
      mydata$fit$raw = rbind(as.matrix(mydata$fit$raw), matrix(0, padweeks, length(fit_ident)))
    } else {
      # pad with zeros to fill-out a complete year
      mydata$fit$raw = c(mydata$fit$raw, rep(0, padweeks))
    }

    # Fill-in MetaILI mydata for 'fit'
    mydata$fit$raw_units = diceData$MetaILI$CDCili_raw_units[fit_ident]
    mydata$fit$epi_units = diceData$MetaILI$CDCili_epi_units[fit_ident]
    mydata$fit$factor = diceData$MetaILI$CDCili_factor[diceData$MetaILI$CDCili_factor$year == start.year, fit_ident]

    # download data for 'all_years_epi'
    # always grab MySQL data so NOAA precip and temp can be appended to NASA sh
    sh_mod_name = mod_name
    sh_mod_name["NAME_2"] = "US"
    if (mod_level>2) sh_mod_name["NAME_3"] = mydata$model$attr$ABBV_3
    if (mod_level>3) sh_mod_name["NAME_4"] = mydata$model$attr$ABBV_4
    start_date = CDCweek2date(CDCweek=27, year=2003)
    mysql_dat = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=sh_mod_name, start.date=start_date, end.date=NULL, disease="flu", sql_data_source=7)

    # Attempt to download 'model' ILI all_years_epi from CDC server
    if (mod_level == 2) {
      CDCreg = "national"
      sub_region = NULL
      sub_name   = "National"
    } else if (mod_level == 3) {
      CDCreg = "hhs"
      sub_region = mydata$model$attr$ID_3
      sub_name   = paste0("Region ", sub_region)
    } else if (mod_level == 4) {
      CDCreg = "state"
      # sub_region = gsub(".", " ", mydata$model$attr$NAME_4, fixed = TRUE)
      sub_region = mydata$model$attr$NAME_4
      sub_name   = gsub(".", " ", mydata$model$attr$NAME_4, fixed=TRUE)
    }

    download_attempt = cdc_download(CDCreg=CDCreg, years=2002:year(Sys.Date()))
    if (download_attempt$GoodPull) {
      CDCmodel = as.data.frame(download_attempt$CDCdata)
      # extract appropriate sub_region
      CDCmodel = CDCmodel[CDCmodel$region==sub_name, ]
      if (mod_level<4) {
        temp_epi = CDCmodel[, c("year", "week", "weighted_ili")]
        temp_epi$weighted_ili[temp_epi$weighted_ili==0] = NA
      } else {
        temp_epi = CDCmodel[, c("year", "week", "unweighted_ili")]
      }

      # generate master dates
      min_year = min(CDCmodel$year)
      min_week = min(CDCmodel$week[CDCmodel$year==min_year])
      min_date = CDCweek2date(CDCweek=min_week, year=min_year)
      min_date = max(min_date, CDCweek2date(27, 2003))
      season_dates = get.StartEndDates(country="US", model_year=year(min_date), disease="flu")
      min_date = min(min_date, season_dates$start_date)
      max_year = max(CDCmodel$year)
      max_week = max(CDCmodel$week[CDCmodel$year==max_year])
      max_date = CDCweek2date(CDCweek=max_week, year=max_year)
      all_years_epi = BuildDateVecs(cadence=2, start.date=min_date, end.date=max_date)
      # merge master dates and data
      temp_epi = merge(x=all_years_epi, y=temp_epi, all.x=T, sort=F)
      temp_epi = temp_epi[order(temp_epi$date), ]
      # build all_years_epi
      all_years_epi$model = list()
      all_years_epi$model$raw = temp_epi[, 4]

    } else {

      # pull all available years from datafile
      # sh_mod_name = mod_name
      # sh_mod_name["NAME_2"] = "US"
      # if (mod_level>2) sh_mod_name["NAME_3"] = mydata$model$attr$ABBV_3
      # if (mod_level>3) sh_mod_name["NAME_4"] = mydata$model$attr$ABBV_4
      # start_date = CDCweek2date(CDCweek=27, year=2003)
      # mysql_dat = get.mysql(mod_level=mod_level, fit_level=mod_level, mod_name=sh_mod_name, start.date=start_date, disease="flu", sql_data_source=7)
      # copy master dates
      all_years_epi = as.list(add_ndays(data.frame(mysql_dat$mydata[1:3]), cadence=2))
      all_years_epi$model = list()
      all_years_epi$model$raw = mysql_dat$mydata$model[[mysql_dat$mydata$model$attr$identifier]]$ILI

    }

    # calc $epi and $gammaepi
    all_years_epi$model$epi = all_years_epi$model$raw*mydata$model$factor
    all_years_epi$model$epi[is.na(all_years_epi$model$epi)] = 0
    all_years_epi$model$gamaepi = lgamma(all_years_epi$model$epi + 1)

    # retrieve SH for model region (from NASA)
    sh_df = get.sh.data(region=CDCreg, sub_region=sub_region, start_year=all_years_epi$year[1], start_week=all_years_epi$week[1], end_year=all_years_epi$year[length(all_years_epi$year)], end_week=all_years_epi$week[length(all_years_epi$week)])
    all_years_epi$model$sh = sh_df[, 4]

    # retrieve temp and precip for model region (from NOAA)
    temp_clim = cbind(data.frame(mysql_dat$mydata[1:3]), mysql_dat$mydata$model[[mysql_dat$mydata$model$attr$identifier]][, c("temp", "precip")])
    temp_clim = merge(x=pluralize(all_years_epi[1:3]), y=temp_clim, all.x=T)
    all_years_epi$model$temp = temp_clim$temp
    all_years_epi$model$precip = temp_clim$precip

    # Attempt to download 'fit' ILI all_years_epi from CDC server
    if (fit_level == 2) {
      CDCreg = "national"
      sub_region = "national"
      sub_name = "National"
      col_names = c("year", "week", "weighted_ili")
    } else if (fit_level == 3) {
      CDCreg = "hhs"
      sub_region = mydata$fit$attr$ID_3
      sub_name   = paste0("Region ", mydata$fit$attr$ID_3)
      col_names = c("year", "week", "weighted_ili")
    } else if (fit_level == 4) {
      CDCreg = "state"
      sub_region = mydata$fit$attr$NAME_4
      sub_name   = gsub(".", " ", mydata$fit$attr$NAME_4, fixed = TRUE)
      col_names = c("year", "week", "unweighted_ili")
    }

    download_attempt = cdc_download(CDCreg=CDCreg, years=2002:year(Sys.Date()))
    if (download_attempt$GoodPull) {
      CDCfit = as.data.frame(download_attempt$CDCdata)
      if (fit_level<4) CDCfit$weighted_ili[CDCfit$weighted_ili==0] = NA
      # extract appropriate sub_region(s)
      temp_epi = as.data.frame(all_years_epi[1:3])
      for (reg in sub_name) {
        reg_ind = CDCfit$region==reg
        reg_epi = CDCfit[reg_ind, col_names]
        names(reg_epi)[3] = reg
        temp_epi = merge(x=temp_epi, y=reg_epi, all.x=T, sort=F)
      }
      # re-order temp_epi
      temp_epi = temp_epi[order(temp_epi$date), ]
      # build all_years_epi
      all_years_epi$fit = list()
      all_years_epi$fit$raw = temp_epi[, 4:length(temp_epi)]

      all_years_epi = pluralize(all_years_epi)

    } else {
      # pull all available years from MySQL
      # sh_mod_name = mod_name
      # sh_mod_name["NAME_2"] = "US"
      # if (mod_level>2) sh_mod_name["NAME_3"] = mydata$model$attr$ABBV_3
      # if (mod_level>3) sh_mod_name["NAME_4"] = mydata$model$attr$ABBV_4
      # start_date = CDCweek2date(CDCweek=27, year=2003)
      # mysql_dat = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=sh_mod_name, start.date=start_date, disease="flu", sql_data_source=7)
      # extract fit regions
      all_years_epi$fit = list()
      temp_epi = data.frame(mysql_dat$all_years_epi[1:3])
      temp_epi = cbind(temp_epi, mysql_dat$all_years_epi$fit$raw)
      if (length(temp_epi)==4) {
        names(temp_epi)[4] = mysql_dat$mydata$fit$attr$identifier
      }
      # match dates with master dates
      temp_epi = merge(x=all_years_epi[1:3], y=temp_epi, all.x=TRUE, sort=F)
      temp_epi = temp_epi[order(temp_epi$date), ]
      # remove date columns
      temp_epi = temp_epi[, 4:length(temp_epi)]
      if (!is.data.frame(temp_epi)) {
        all_years_epi$fit$raw = temp_epi
      } else {
        # match regions with mydata$fit
        # temp_epi = temp_epi[, names(temp_epi) %in% mydata$fit$attr$identifier]
        all_years_epi$fit$raw = temp_epi
      }

    }
    # calc $epi and $gammaepi
    all_years_epi$fit$epi = all_years_epi$fit$raw
    all_years_epi$fit$gamaepi = all_years_epi$fit$raw
    if (!is.null(ncol(all_years_epi$fit$epi))) {
      for (ii in 1:ncol(all_years_epi$fit$epi)) {
        all_years_epi$fit$epi[, ii] = all_years_epi$fit$raw[, ii] * as.numeric(mydata$fit$factor[ii])
        all_years_epi$fit$epi[is.na(all_years_epi$fit$epi[, ii]), ii] = 0
        all_years_epi$fit$gamaepi[, ii] = lgamma(all_years_epi$fit$epi[, ii]+1)
      }
    } else {
      all_years_epi$fit$epi     = all_years_epi$fit$raw * mydata$fit$factor
      all_years_epi$fit$epi[is.na(all_years_epi$fit$epi)] = 0
      all_years_epi$fit$gamaepi = lgamma(all_years_epi$fit$epi+1)
    }

    # retrieve SH for all_years_epi fit regions (from NASA)
    sh_df = get.sh.data(region=CDCreg, sub_region=sub_region, start_year=all_years_epi$year[1], start_week=all_years_epi$week[1], end_year=all_years_epi$year[length(all_years_epi$year)], end_week=all_years_epi$week[length(all_years_epi$week)])
    all_years_epi$fit$sh = sh_df[, 4:ncol(sh_df)]

    # retrieve temp and precip for fit regions (from NOAA)
    temp_precip = cbind(data.frame(mysql_dat$all_years_epi[1:3]), mysql_dat$all_years_epi$fit$precip)
    temp_precip = merge(x=all_years_epi[1:3], y=temp_precip, all.x=T)
    all_years_epi$fit$precip = temp_precip[4:length(temp_precip)]
    temp_temp = cbind(data.frame(mysql_dat$all_years_epi[1:3]), mysql_dat$all_years_epi$fit$temp)
    temp_temp = merge(x=all_years_epi[1:3], y=temp_temp, all.x=T)
    all_years_epi$fit$temp = temp_temp[4:length(temp_temp)]
    if (length(all_years_epi$fit$temp)==1) {
      all_years_epi$fit$precip = unlist(all_years_epi$fit$precip, use.names=F)
      all_years_epi$fit$temp   = unlist(all_years_epi$fit$temp, use.names=F)
    }

    # pluralize date vecs
    all_years_epi = pluralize(all_years_epi)

    mydata$data_source = "cdc"
    mydata$data_desc = DICE_data_sources[tolower(DICE_data_sources$source_abbv)=="cdc", ]
  }  # end if data_source == 'GFT' .... 'CDC'

  mydata$nperiodsData = length(start.ind:end.ind)

  # Use raw to calc epi, gamaepi for model profile
  mydata$model$epi = round(mydata$model$raw * mydata$model$factor)
  mydata$model$epi[is.na(mydata$model$epi)] = 0
  mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

  # Use raw to calc epi, gamaepi for fit profile(s)
  mydata$fit$epi = as.matrix(mydata$fit$raw)
  if (length(fit_ident) > 1) {
    for (ii in 1:length(fit_ident)) {
      mydata$fit$epi[, ii] = round(mydata$fit$raw[, ii] * mydata$fit$factor[[ii]])
    }
  } else {
    mydata$fit$epi = round(mydata$fit$raw * mydata$fit$factor)
  }
  mydata$fit$epi[is.na(mydata$fit$epi)] = 0
  mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)


  # Pull SH for fit profile(s)
  Good_Pull = FALSE
  if (endSeason >= curSeason & fit_level < 5) {
    # if mydata comes from current season, pull directly from the mydata server
    if (fit_level == 2) {
      region = "national"
      sub_region = ""
    } else if (fit_level == 3) {
      region = "hhs"
      sub_region = as.integer(sub("^.*\\R", "", fit_ident))
    } else if (fit_level == 4) {
      region = "state"
      sub_region = mydata$fit$name
    }

    sh = try(get.sh.data(region = region, sub_region = sub_region, start_year = 2000, start_week = 1, end_year = end.year, end_week = end.week),
             silent = TRUE)
    if (!is(sh, "try-error")) {
      Good_Pull = TRUE
      print("Fit-level specific humidity mydata successfully downloaded from PSI server.")
      # Check if all weeks are present
      nsh = length(sh$date)
      if (sh$week[nsh] < end.week | sh$year[nsh] < end.year) {
        # if not, add appropriate number of weeks to sh
        start.ind = which(sh$year == start.year & sh$week == start.week)
        mydata$fit$sh = sh[start.ind:nsh, 4:ncol(sh)]
        shEndWeek = sh$week[nsh]
        shEndYear = sh$year[nsh]
        date_ind = which(mydata$weeks == shEndWeek & mydata$years == shEndYear)
        predWeeks = (date_ind + 1):length(mydata$weeks)
        # pad mydata$fit$sh with NAs if more than one mydata column
        if (length(sh) > 4) {
          mydata$fit$sh[predWeeks, ] = NA
        } else {
          # assume mydata$fit$sh is a vector
          mydata$fit$sh[predWeeks] = NA
        }

      } else {
        start.ind = which(sh$year == start.year & sh$week == start.week)
        end.ind = which(sh$year == end.year & sh$week == end.week)
        mydata$fit$sh = sh[start.ind:end.ind, 4:ncol(sh)]
      }

      # If any NAs appear in the sh mydata, replace with historic average
      if (length(sub_region) > 1) {
        for (kk in 1:length(sub_region)) {
          NA_ind = which(is.na(mydata$fit$sh[, kk]))
          for (ind in NA_ind) {
            week = mydata$weeks[ind]
            WeekInd = sh$week == week
            mydata$fit$sh[ind, kk] = mean(sh[WeekInd, kk + 3], na.rm = TRUE)
          }
        }
      } else {
        NA_ind = which(is.na(mydata$fit$sh))
        for (ind in NA_ind) {
          week = mydata$weeks[ind]
          WeekInd = sh$week == week
          mydata$fit$sh[ind] = mean(sh[WeekInd, 4], na.rm = TRUE)
        }
      }

    } else {
      print("WARNING: Specific Humidity download from PSI server was unsuccessful for fit-level mydata. Reverting to SH mydata in local DICE-package mydatabase.  DICE package mydata is likely less up-to-date than the PSI server.")
    }
  }
  if (!Good_Pull) {
    # older mydata is pulled from DICE mydataset
    start.ind = which(diceData$sh$year == start.year & diceData$sh$week == start.week)
    end.ind = which(diceData$sh$year == end.year & diceData$sh$week == end.week)
    mydata$fit$sh = diceData$sh[start.ind:end.ind, fit_ident]
  }
  names(mydata$fit$sh) = fit_ident

  # Pull SH for model profile
  Good_Pull = FALSE
  if (endSeason >= curSeason & mod_level < 5) {
    # if mydata comes from current season, pull directly from the mydata server
    if (mod_level == 2) {
      region = "national"
    } else if (mod_level == 3) {
      region = "hhs"
      sub_region = as.integer(sub("^.*\\R", "", mod_ident))
    } else if (mod_level == 4) {
      region = "state"
      sub_region = mydata$model$name
      # state_Abbs = sub('^.*\\.','',fit_ident) sub_region = c() for (ii in 1:length(state_Abbs)) { if (state_Abbs[ii] == 'DC') {
      # sub_region[ii] = 'District.of.Columbia' } else { sub_region[ii] = state.name[state_Abbs[ii]==state.abb] } } sub_region = gsub('
      # ','.',sub_region)
    }

    sh = try(get.sh.data(region = region, sub_region = sub_region, start_year = 2000, start_week = 1, end_year = end.year, end_week = end.week),
             silent = TRUE)
    if (!is(sh, "try-error")) {
      Good_Pull = TRUE
      print("Model-level specific humidity mydata successfully downloaded from PSI server.")

      # Check if all weeks are present
      nsh = length(sh$date)
      if (sh$week[nsh] < end.week | sh$year[nsh] < end.year) {
        # if not, add appropriate number of weeks to sh
        start.ind = which(sh$year == start.year & sh$week == start.week)
        mydata$model$sh = sh[start.ind:nsh, 4:ncol(sh)]
        shEndWeek = sh$week[nsh]
        shEndYear = sh$year[nsh]
        date_ind = which(mydata$weeks == shEndWeek & mydata$years == shEndYear)
        predWeeks = (date_ind + 1):length(mydata$weeks)
        # pad mydata$model$sh with NAs
        mydata$model$sh[predWeeks] = NA
      } else {
        start.ind = which(sh$year == start.year & sh$week == start.week)
        end.ind = which(sh$year == end.year & sh$week == end.week)
        mydata$model$sh = sh[start.ind:end.ind, 4:ncol(sh)]
      }

      # If any NAs appear in the sh mydata, replace with historic average
      NA_ind = which(is.na(mydata$model$sh))
      for (ind in NA_ind) {
        week = mydata$weeks[ind]
        WeekInd = sh$week == week
        mydata$model$sh[ind] = mean(sh[WeekInd, 4], na.rm = TRUE)
      }
    } else {
      print("WARNING: Specific Humidity download from PSI server was unsuccessful for model-level mydata. Reverting to SH mydata in local DICE-package mydatabase.  DICE package mydata is likely less up-to-date than the PSI server.")
    }
  }
  if (!Good_Pull) {
    # older mydata is pulled from DICE mydataset
    start.ind = which(diceData$sh$year == start.year & diceData$sh$week == start.week)
    end.ind = which(diceData$sh$year == end.year & diceData$sh$week == end.week)
    mydata$model$sh = diceData$sh[start.ind:end.ind, mod_ident]
  }

  if (tolower(data_source)=="cdc") {
    # pull temp and precip data from mysql for model
    temp_clim = cbind(data.frame(mysql_dat$mydata[1:3]), mysql_dat$mydata$model[[mysql_dat$mydata$model$attr$identifier]][, c("temp", "precip")])
    temp_clim = merge(x=mydata[1:3], y=temp_clim, all.x=T)
    mydata$model$temp   = temp_clim$temp
    mydata$model$precip = temp_clim$precip
    # fill any NAs with historic average
    for (var in c("temp", "precip")) {
      na_ind = which(is.na(mydata$model[[var]]))
      for (ind in na_ind) {
        week = mydata$weeks[ind]
        week_ind = all_years_epi$weeks==week
        mydata$model[[var]][ind] = mean(all_years_epi$model[[var]][week_ind], na.rm=T)
      }
    }
    # pull temp and precip data from mysql for fit
    temp_precip = cbind(data.frame(all_years_epi[1:3]), all_years_epi$fit$precip)
    temp_precip = merge(x=mydata[1:3], y=temp_precip, all.x=T)
    mydata$fit$precip = temp_precip[4:length(temp_precip)]
    temp_temp = cbind(data.frame(all_years_epi[1:3]), all_years_epi$fit$temp)
    temp_temp = merge(x=mydata[1:3], y=temp_temp, all.x=T)
    mydata$fit$temp = temp_temp[4:length(temp_temp)]
    # fill any NAs with historic average
    for (var in c("temp", "precip")) {
      for (ii in 1:length(mydata$fit[[var]])) {
        na_ind = which(is.na(mydata$fit[[var]][, ii]))
        for (ind in na_ind) {
          week = mydata$weeks[ind]
          week_ind = all_years_epi$weeks==week
          mydata$fit[[var]][ind, ii] = mean(all_years_epi$model[[var]][week_ind], na.rm=T)
        }
      }
    }
    if (length(mydata$fit$temp)==1) {
      mydata$fit$precip = unlist(mydata$fit$precip, use.names=F)
      mydata$fit$temp   = unlist(mydata$fit$temp, use.names=F)
    }
  }


  mydata$fit$nregions = length(mydata$fit$pop)

  mydata$model$attr$identifier = mod_ident
  mydata$fit$attr$identifier = fit_ident

  if (tolower(data_source)=="cdc") {
    return(list(mydata=mydata, all_years_epi=all_years_epi))
  } else {
    return(mydata)
  }

}


get.subset2 <- function(mod_level = 2, fit_level = 3, mod_name = NULL, fit_name = NULL, start.year = 2015, start.week = 27, end.year = 2016,
    end.week = 26, dataType = "cdc") {
    # Load the full \pkg{dice} mydataset and then subset by date and region
    #
    # This function populates the $model and $fit entries in \code{mydata}.  It requires a more complicated input than get.DICE.mydata/get.subset, but allows for much more flexibility when pulling mydata for DICE.
    # @param mod_level Integer - Spatial level of the model data. Spatial levels for the United States are defined as follows:
    # \tabular{rl}{
    # 2 \tab Country \cr
    # 3 \tab HHS Region \cr
    # 4 \tab State \cr
    # 5 \tab County \cr
    # 6 \tab City \cr}
    # @param mod_name - two options:
    # \itemize{
    #   \item{Character Vector that describes the model population.  For example \code{mod_name=c(NAME_2='USA',NAME_3='Region2',NAME_4='NY')} combined with \code{mod_level=4} specifies the state of New York.  \code{NAME_X} may be the fullname or abbreviation of the \code{X}-level population.  For instance \code{mod_name=c(NAME_2='United States',NAME_3='R2',NAME_4='New York')} can also be used to specify New York.}
    #   \item{For \code{mod_name='create'} the populations specified by fit_level/fit_name will be used to create a model population.  Currently dataTypes 'CDC', 'GFT', and 'MiscILI' are supported for the 'create' functionality.}}
    # @param fit_level If \code{mod_name='create'}, then \code{fit_level} is a vector of integers indicating the spatial level of each region to be used. Otherwise \code{fit_level} is a single integer.
    # @param fit_name \itemize{
    #   \item{If \code{mod_name='create'}, then \code{fit_name} is a list of character vectors: \code{fit_name = list(c(NAME_2='USA',NAME_3='Region10'),c(NAME_2='USA',NAME_3='Region9',NAME_4='AZ'),c(NAME_2='USA',NAME_3='Region9',NAME_4='CA'))}. This example combined with \code{mod_name='create'} and \code{fit_level=c(3,4,4)} creates a super-region from HHS Region 10, California, and Arizona.}
    #   \item{When mod_name/mod_level is not in 'create' mode, the only option currently supported is \code{fit_name='all'}. This will pull all populations at fit_level that are subpopulations of the mod_name population.}}
    # @param dataType Select from mydata-types: 'CDC', 'GFT', 'MiscILI', 'MiscCases', 'ZikaS', or 'ZikaC' (NOT case sensitive).  These indicate the source and type of mydata: Centers for Disease Control (CDC) mydata is \%ILI, Google Flu Trends (GFT) is \%ILI, 'MiscILI' is \%ILI mydata from miscelaneous sources, 'MiscCases' is absolute-cases mydata from miscelaneous sources, 'ZikaS' is Zika suspected cases (PAHO), and 'ZikaC' is Zika confirmed cases (PAHO). \itemize{
    # \item{If \code{mod_name='create'}, then dataType is a vector of strings indicating the mydata-type to be used for each fit_name population. \code{dataType = c('CDC','GFT','MiscCases')}}
    # \item{Otherwise dataType is a single string that determines the mydata-type for all model/fit populations.}}
    # @param start.year An integer - start year of the flu season
    # @param end.year  An integer  - end year of the flu season (default is end.year = start.year + 1)
    # @param start.week An integer - starting CDC week of the flu season (default is 27)
    # @param end.week   An integer - ending CDC week of the flu season (default is 26)
    # @return mydata  A list ILI, SH, School, and Census mydata for both the model- and fit-level region(s).
    # @examples
    # require(DICE)
    # To pull mydata for a single region:
    # mydata = get.subset2(mod_level = 3, mod_name=c(NAME_2='USA',NAME_3='Region4'),fit_level = 3, fit_name = 'all', start.year = 2014, end.year=2015, dataType = 'cdc')
    #
    # OR
    #
    # mydata = get.subset2(mod_level = 5, mod_name=c(NAME_2='USA',NAME_3='Region9',NAME_4='CA',NAME_5='SD'),fit_level = 5, fit_name = 'all', start.year = 2014, end.year=2015, dataType = 'MiscCases')
    #
    # To pull mydata for all states in a region:
    # mydata = get.subset2(mod_level = 3, mod_name=c(NAME_2='USA',NAME_3='Region4'),fit_level = 4, fit_name = 'all', start.year = 2014, end.year=2015, dataType = 'gft')
    #
    # To create a custom region from a few states:
    # mydata = get.subset2(dataType=c('gft','gft','gft'), mod_level=3, mod_name='create', fit_level=c(4,4,4), fit_name=list(c(NAME_2='USA',NAME_3='Region10',NAME_4='WA'),c(NAME_2='USA',NAME_3='Region10',NAME_4='OR'),c(NAME_2='USA',NAME_3='R9',NAME_4='CA')), start.year=2014, end.year=2015)

    # Initialize 'mydata'
    mydata = list()
    # Current year/season. Used to determine if mydata should be downloaded from CDC
    today = Sys.Date()
    CDCweek = Date2CDCweek(today)
    if (CDCweek$CDCweek > 26) {
        curYear = CDCweek$CDCyear
    } else {
        curYear = CDCweek$CDCyear - 1
    }

    # load the mydataset
    mydata(DICE_mydataset, package = "DICE")
    # returns a mydataset/list called 'diceData'

    # Determine what time indices are included
    start.ind = which(diceData$school$year == start.year & diceData$school$week == start.week)
    end.ind = which(diceData$school$year == end.year & diceData$school$week == end.week)

    mydata$dates = diceData$school$date[start.ind:end.ind]
    mydata$weeks = diceData$school$week[start.ind:end.ind]
    mydata$years = diceData$school$year[start.ind:end.ind]
    mydata$nperiods = length(mydata$weeks)

    # Index all (sub)populations in $attr Start with the 'model' population. In this step, also find all sub-populations of 'model' and
    # index in 'fit_ind'.
    if (length(mod_name) == 1 && tolower(mod_name) == "create") {
        # Do nothing. 'model' info will be generated by combining the fit populations
        mod_ind = rep(FALSE, length(diceData$attr$object))
    } else {
        if (mod_level >= 2) {
            fit_ind = rep(TRUE, length(diceData$attr$object))
            for (lev in 2:mod_level) {
                fit_ind = fit_ind & (toupper(mod_name[paste0("NAME_", lev)]) == toupper(diceData$attr[, paste0("NAME_", lev)]) | toupper(mod_name[paste0("NAME_",
                  lev)]) == diceData$attr[, paste0("ABBV_", lev)])
            }
            mod_ind = fit_ind & mod_level == diceData$attr$level
        } else {
            # Global/Continent
            mod_ind = (toupper(mod_name[paste0("NAME_", mod_level)]) == diceData$attr[, paste0("ABBV_", mod_level)] | toupper(mod_name[paste0("NAME_",
                mod_level)]) == toupper(diceData$attr[, paste0("NAME_", mod_level)])) & mod_level == diceData$attr$level
            fit_ind = (toupper(mod_name[paste0("NAME_", mod_level)]) == diceData$attr[, paste0("ABBV_", mod_level)] | toupper(mod_name[paste0("NAME_",
                mod_level)]) == toupper(diceData$attr[, paste0("NAME_", mod_level)])) & mod_level > diceData$attr$level
        }
    }

    # Reduce 'fit_ind' to the sub-populations specified in fit_level and fit_name
    if (length(fit_name) == 1 && fit_name == "best") {
        # do something clever to fuse the smallest-scale available mydata to a representative 'fit' set in this step, simply return all
        # populations in attr that are within the mod population (already done)
    } else {
        if (length(fit_name) == 1 && fit_name == "all") {
            # Find all subpopulations at fit_level within mod_name
            if (length(fit_level) != 1) {
                stop("For fit_name option 'all', fit_level must be a single integer.")
            }
            fit_ind = fit_ind & fit_level == diceData$attr$level
        } else {
            # Find fit_level entries one at a time
            if (length(fit_level) != length(fit_name) || length(fit_level) != length(dataType)) {
                stop("For multiple fit_name specification: fit_level,  fit_name, and dataType must have the same number of elements.")
            }
            reorder_Dtype = integer(length = length(diceData$attr$object))
            fit_ind = rep(FALSE, length(diceData$attr$object))
            for (ii in 1:length(fit_level)) {
                temp_ind = diceData$attr$level == fit_level[ii]
                for (lev in 2:fit_level[ii]) {
                  temp_ind = temp_ind & (toupper(fit_name[[ii]][paste0("NAME_", lev)]) == toupper(diceData$attr[, paste0("NAME_", lev)]) |
                    toupper(fit_name[[ii]][paste0("NAME_", lev)]) == diceData$attr[, paste0("ABBV_", lev)])
                }
                if (!any(temp_ind)) {
                  stop("fit_name entry number ", ii, " does not match any DICE populations.")
                } else {
                  reorder_Dtype[temp_ind] = ii
                  fit_ind = fit_ind | temp_ind
                }
            }
            dataType = dataType[reorder_Dtype[reorder_Dtype != 0]]
        }
    }

    # Import attr info
    ModObj = which(mod_ind)
    FitObj = which(fit_ind)

    if (length(mod_name) == 1 && tolower(mod_name) == "create") {
        # do not attempt to populate $mod$attr
        mod_ident = character()
    } else {
        if (any(mod_ind)) {
            mydata$model$attr = diceData$attr[ModObj, ]
            mydata$model$level = mydata$model$attr$level
            mydata$model$name = diceData$attr[ModObj, paste0("NAME_", mod_level)]
        } else {
            stop("The mod_level/mod_name combination entered did not match any mydata entries.")
        }
        # Create a mod_level identifier
        if (mod_level >= 2) {
            mod_ident = diceData$attr$ABBV_2[ModObj]
            if (mod_level > 2) {
                for (lev in 3:mod_level) {
                  mod_ident = paste0(mod_ident, ".", diceData$attr[ModObj, paste0("ABBV_", lev)])
                }
            }
        } else {
            mod_ident = diceData$attr[, paste0("ABBV_", mod_level)][ModObj]
        }
    }

    mydata$fit$attr = diceData$attr[fit_ind, ]
    mydata$fit$level = mydata$fit$attr$level
    # Create a vector of fit_level identifiers
    fit_ident = character(length(FitObj))
    for (ii in 1:length(FitObj)) {
        obj = FitObj[ii]
        level = mydata$fit$attr$level[ii]
        if (level >= 2) {
            fit_ident[ii] = diceData$attr$ABBV_2[obj]
            if (level > 2) {
                for (lev in 3:level) {
                  fit_ident[ii] = paste0(fit_ident[ii], ".", diceData$attr[obj, paste0("ABBV_", lev)])
                }
            }
        } else {
            fit_ident[ii] = diceData$attr[, paste0("ABBV_", level)][obj]
        }
        mydata$fit$name[ii] = diceData$attr[, paste0("NAME_", level)][obj]
    }

    # Pull school mydata, onset, pop
    mydata$fit$school = diceData$school[start.ind:end.ind, fit_ident]
    # Should change code to understand NA as 'no-mydata' mydata$fit$school[is.na(mydata$fit$school)] = 0
    mydata$fit$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, fit_ident])
    mydata$fit$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, fit_ident]
    mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)

    # Import mod_levelschool, population, and onset mydata
    mydata$model$school = diceData$school[start.ind:end.ind, mod_ident]
    # mydata$model$school[is.na(mydata$model$school)] = 0
    mydata$model$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, mod_ident])
    mydata$model$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, mod_ident]


    # Associate dataType input with appropriate mydataframe
    FrameName = c()
    for (ii in 1:length(dataType)) {
        FrameName[ii] = switch(tolower(dataType[ii]), cdc = "CDCili", gft = "GFTili", miscili = "MiscIli", misccases = "MiscCases", zikac = "ZikaConf",
            zikas = "ZikaSusp")
    }


    # Check that model and fit data exist in the specified dataType
    if (!(length(mod_name) == 1 && tolower(mod_name) == "create")) {
        mydataNames = colnames(diceData[[FrameName]])
        if (!(mod_ident %in% mydataNames)) {
            stop("DICE does not have the specified mod_level mydata for mydata-type: ", dataType)
        }
    } else {
        for (ii in 1:length(FrameName)) {
            mydataNames = colnames(diceData[[FrameName[ii]]])
            if (!(fit_ident[ii] %in% mydataNames)) {
                stop(mydata$fit$name[ii], " does not have an entry in the ", FrameName[ii], " mydatabase.")
            }
        }
    }

    # Get ILI mydata
    if (!(length(mod_name) == 1 && tolower(mod_name) == "create")) {
        # Get model incidence mydata
        temp_dat = get.ili(FrameName = FrameName, ident = mod_ident, start.year = start.year, end.year = end.year, start.week = start.week,
            end.week = end.week, curYear = curYear, nperiods = mydata$nperiods)
        mydata$model$raw = temp_dat$raw
        # Get fit incidence mydata
        temp_dat = get.ili(FrameName = FrameName, ident = fit_ident, start.year = start.year, end.year = end.year, start.week = start.week,
            end.week = end.week, curYear = curYear, nperiods = mydata$nperiods)
        mydata$fit$raw = temp_dat$raw
        nperiodsData = temp_dat$nperiodsData


        # Pull model-incidence Meta-mydata
        temp_meta = get.meta(FrameName = FrameName, ident = mod_ident, start.year = start.year)
        mydata$model$raw_units = temp_meta$raw_units
        mydata$model$epi_units = temp_meta$epi_units
        mydata$model$factor = temp_meta$factor

        # Pull fit-incidence Meta-mydata
        temp_meta = get.meta(FrameName = FrameName, ident = fit_ident, start.year = start.year)
        mydata$fit$raw_units = temp_meta$raw_units
        mydata$fit$epi_units = temp_meta$epi_units
        mydata$fit$factor = temp_meta$factor

    } else {
        # if mod_name=='create' Determine which dataTypes are needed
        Types = unique(FrameName)

        ii = 1
        # for each mydata type, pull and add to mydata$fit$raw
        temp_ident = fit_ident[FrameName == Types[ii]]
        temp_dat = get.ili(FrameName = Types[ii], ident = temp_ident, start.year = start.year, end.year = end.year, start.week = start.week,
            end.week = end.week, curYear = curYear, nperiods = mydata$nperiods)
        mydata$fit$raw = data.frame(temp_dat$raw)
        names(mydata$fit$raw) = temp_ident
        nperiodsData = temp_dat$nperiodsData

        if (length(Types) > 1) {
            for (ii in 2:length(Types)) {
                # for each mydata type, pull and add to mydata$fit$raw
                temp_ident = fit_ident[FrameName == Types[ii]]
                temp_dat = get.ili(FrameName = Types[ii], ident = temp_ident, start.year = start.year, end.year = end.year, start.week = start.week,
                  end.week = end.week, curYear = curYear, nperiods = mydata$nperiods)
                mydata$fit$raw[temp_ident] = temp_dat$raw
                # Need a better procedure when different dataTypes have different nperiodsData
                nperiodsData = min(nperiodsData, temp_dat$nperiodsData)
            }
        }
        # Reorder raw
        mydata$fit$raw = mydata$fit$raw[fit_ident]

        # Pull fit-incidence Meta-mydata
        temp_meta = get.meta(FrameName = FrameName, ident = fit_ident, start.year = start.year)
        mydata$fit$raw_units = temp_meta$raw_units
        mydata$fit$epi_units = temp_meta$epi_units
        mydata$fit$factor = temp_meta$factor
    }



    mydata$nperiodsData = length(start.ind:end.ind)

    # Use raw to calc epi, gamaepi for model profile
    if (!(length(mod_name) == 1 && tolower(mod_name) == "create")) {
        mydata$model$epi = round(mydata$model$raw * mydata$model$factor)
        # mydata$model$epi[is.na(mydata$model$epi)] = 0
        mydata$model$gamaepi = lgamma(mydata$model$epi + 1)
    }

    # Use raw to calc epi, gamaepi for fit profile(s)
    mydata$fit$epi = as.matrix(mydata$fit$raw)
    if (length(fit_ident) > 1) {
        for (ii in 1:length(fit_ident)) {
            mydata$fit$epi[, ii] = round(mydata$fit$raw[, ii] * mydata$fit$factor[[ii]])
        }
    } else {
        mydata$fit$epi = round(mydata$fit$raw * mydata$fit$factor)
    }
    # mydata$fit$epi[is.na(mydata$fit$epi)] = 0
    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)

    # Pull SH for fit profile(s)
    start.ind = which(diceData$sh$year == start.year & diceData$sh$week == start.week)
    end.ind = which(diceData$sh$year == end.year & diceData$sh$week == end.week)
    mydata$fit$sh = diceData$sh[start.ind:end.ind, fit_ident]
    mydata$fit$nregions = length(mydata$fit$pop)

    # Pull SH for model profile
    mydata$model$sh = diceData$sh[start.ind:end.ind, mod_ident]

    # For the 'create' option, construct mydata$model from mydata$fit
    if (length(mod_name) == 1 && tolower(mod_name) == "create") {
        # first update the coefficients in fit
        mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)

        # Create a name for the new model patch
        mydata$model$level = mod_level
        mydata$model$name = mydata$fit$attr[, paste0("ABBV_", mydata$fit$level[1])][1]
        if (length(mydata$fit$level) > 1) {
            for (ii in 2:length(mydata$fit$level)) {
                mydata$model$name = paste0(mydata$model$name, ".", mydata$fit$attr[, paste0("ABBV_", mydata$fit$level[ii])][ii])
            }
        }

        # Port some values from fit$attr to model$attr
        mydata$model$attr$object = NA
        mydata$model$attr$level = mod_level
        for (lev in 1:(mod_level - 1)) {
            mydata$model$attr[[paste0("NAME_", lev)]] = mydata$fit$attr[[paste0("NAME_", lev)]][1]
            mydata$model$attr[[paste0("ABBV_", lev)]] = mydata$fit$attr[[paste0("ABBV_", lev)]][1]
            mydata$model$attr[[paste0("ID_", lev)]] = mydata$fit$attr[[paste0("ID_", lev)]][1]
        }
        mydata$model$attr[[paste0("NAME_", mod_level)]] = mydata$model$name
        mydata$model$attr[[paste0("ABBV_", mod_level)]] = ""
        mydata$model$attr[[paste0("ID_", mod_level)]] = NA
        # with population-weighted lat lon
        mydata$model$attr$lat = sum(mydata$fit$attr$lat * mydata$fit$pop)/sum(mydata$fit$pop)
        mydata$model$attr$lon = sum(mydata$fit$attr$lon * mydata$fit$pop)/sum(mydata$fit$pop)

        # population-weighted average fields
        mydata$model$school = (mydata$fit$pop %*% t(as.matrix(mydata$fit$school)))/sum(mydata$fit$pop)
        mydata$model$sh = (mydata$fit$pop %*% t(as.matrix(mydata$fit$sh)))/sum(mydata$fit$pop)
        mydata$model$onset = sum(mydata$fit$onset * mydata$fit$pop)/sum(mydata$fit$pop)
        mydata$model$raw = (mydata$fit$pop %*% t(as.matrix(mydata$fit$raw)))/sum(mydata$fit$pop)
        # summed fields
        mydata$model$pop = sum(mydata$fit$pop)
        mydata$model$factor = sum(mydata$fit$factor)
        # newly calculated fields
        mydata$model$epi = round(apply(mydata$fit$epi, 1, sum, na.rm = TRUE))
        # mydata$model$epi = round(mydata$model$raw*mydata$model$factor)
        mydata$model$gamaepi = lgamma(mydata$model$epi + 1)
    }

    return(mydata)
}



get_disease_mydata_old <- function(countries="all", years="all", disease="dengue", daily_clim=FALSE) {

  # open mydatabase connection
  myDB = OpenCon()
  # retrieve data_sources table
  data_sources = dbReadTable(myDB, name="data_sources")

  available_diseases = unique(data_sources$disease)

  clim_by_disease = dbReadTable(myDB, name="clim_by_disease")
  clim_names = clim_by_disease$clim_names[clim_by_disease$disease==disease]
  clim_names = strsplit(clim_names, ";")[[1]]

  if (!any(disease==available_diseases)) {
    stop("The disease entered c(", disease, "), does not match available diseases c(", paste(available_diseases, collapse=", "), ")." )
  }

  dis_lut = dbReadTable(myDB, name=paste0(disease,"_lut"))

  # construct the WHERE string(s) for MySQL query
  where = character(2)
  if (length(countries)==1 && countries=="all") {
    where[1] = ""
  } else {
    master_keys = dis_lut$master_key[dis_lut$NAME_2 %in% countries]
    where[1] = paste0("master_key IN('",paste(master_keys, collapse="','"),"')")
  }

  if (length(years)==1 && years=="all") {
    where[2] = ""
  } else {
    min_date = as.Date(paste0(min(years)-1, "-12-28"))
    max_date = as.Date(paste0(max(years), "-12-31"))
    where[2] = paste0("date BETWEEN '", min_date, "' AND '", max_date, "'")
  }

  if (any(where!="")) {
    where = where[where!=""]
    where_string = paste0("WHERE ", paste(where, collapse=" AND "), " ")
  } else {
    where_string = ""
  }

  cat("Downloading climate data......")
  query_string = paste0("SELECT * FROM ",disease,"_data ", where_string, "ORDER BY master_key, source_key, date")
  inc_data = dbGetQuery(myDB, statement=query_string)
  cat("Complete\n")

  inc_data$date = as.Date(inc_data$date)

  # add auxilary columns to inc_data
  source_ind = match(inc_data$source_key, data_sources$source_key)
  inc_data$cadence = data_sources$cadence[source_ind]
  # add year, month, week, columns
  inc_data = Date2CadLabel(query_result=inc_data)

  if (length(years)>1 || years!="all") {
    # trim mydata to appropriate years
    inc_data = inc_data[inc_data$year %in% years, ]
  }

  master_keys = unique(inc_data$master_key)
  cads = unique(inc_data$cadence)

  if (daily_clim) {
    # retrieve clim_by_disease table
    clim_by_disease = dbReadTable(conn=myDB, name="clim_by_disease")
    clim_cols = strsplit(clim_by_disease$clim_names[clim_by_disease$disease==disease], ";")[[1]]
    # retrieve NOAA climate data from MySQL
    min_date = min(inc_data$date)
    max_date = max(inc_data$date)
    # dates are the begining of the cadence cycle, and climate data is daily. So adjust max_date such that all possible days are pulled from noaa_clim
    for (cad in cads) {
      temp_max = max(inc_data$date[inc_data$cadence==cad])
      if (cad==2) {
        temp_max = temp_max + 6
      } else if (cad==3) {
        temp_max = temp_max + 30
      } else if (cad==4) {
        temp_max = temp_max + 366
      } else {
        temp_max = temp_max
      }
      max_date = max(max_date, temp_max)
    }

    cat("Downloading daily climate data......")
    query_string = paste0("SELECT master_key, date, sh, temp, precip FROM noaa_clim WHERE date BETWEEN '", min_date, "' AND '", max_date, "' AND master_key IN('", paste(master_keys, collapse="','"), "') ORDER BY master_key, date")
    noaa_clim = dbGetQuery(myDB, statement=query_string)
    cat("Complete\n")
  }

  # convert mydata into useable and user-friendly format
  dis_lut = dis_lut[dis_lut$master_key %in% master_keys, ]
  out = list()
  out$lut = dis_lut

  countries = dis_lut$NAME_2[dis_lut$level==2]
  cadences = unique(inc_data$cadence)
  cad_names = CadenceNum2Name(cadences)

  cat("Formatting mydata......")
  for (ii in 1:length(cadences)) {
    out[[cad_names[ii]]] = list()
    cad_list_ind = which(names(out)==cad_names[ii])
    cad_data_ind = inc_data$cadence==cadences[ii]

    # establish which columns of inc_data are needed
    data_names = names(inc_data)[substr(names(inc_data), start=1, stop=4)=="mydata"]
    n_data_names = length(data_names)
    if (cadences[ii] %in% c(0,1)) {
      data_cols = c("date", "year", "day", data_names, clim_names)
      max_day_fun <- function(max_date) {
        return(max_date)
      }
    } else if (cadences[ii]==2) {
      data_cols = c("date", "year", "week", data_names, clim_names)
      max_day_fun <- function(max_date) {
        return(max_date+6)
      }
    } else if (cadences[ii]==3) {
      data_cols = c("date", "year", "month", data_names, clim_names)
      max_day_fun <- function(max_date) {
        month(max_date) = month(max_date)+1
        return(max_date-1)
      }
    } else if (cadences[ii]==4) {
      data_cols = c("date", "year", data_names, clim_names)
      max_day_fun <- function(max_date) {
        year(max_date) = year(max_date)+1
        return(max_date-1)
      }
    }

    # determine cadence patches
    master_keys = unique(inc_data$master_key[inc_data$cadence==cadences[ii]])
    for (jj in 1:length(master_keys)) {
      patch_ind = which(dis_lut$master_key==master_keys[jj])
      ident = dis_lut$identifier[patch_ind]
      out[[cad_list_ind]][[ident]] = list()

      sources = unique(inc_data$source_key)
      for (source_key in sources) {
        data_ind = inc_data$master_key==master_keys[jj] & cad_data_ind & inc_data$source_key==source_key
        source_abbv = data_sources$source_abbv[data_sources$source_key==source_key]

        # pull climate data
        out[[cad_list_ind]][[ident]][[source_abbv]] = inc_data[data_ind, data_cols]
        # set 'mydata' column names from source info
        col_names = data_sources$col_names[data_sources$source_key==source_key]
        col_names = strsplit(col_names, ";")[[1]]
        for (kk in 1:length(col_names)) {
          names(out[[cad_list_ind]][[ident]][[source_abbv]])[names(out[[cad_list_ind]][[ident]][[source_abbv]])==paste0("data",kk)] = col_names[kk]
        }
        # clean-up any extra mydata columns
        if (n_data_names > data_sources$data_cols[data_sources$source_key==source_key]) {
          out[[cad_list_ind]][[ident]][[source_abbv]][,paste0("data",(data_sources$data_cols[data_sources$source_key==source_key]+1):n_data_names)] = NULL
        }
      }

      if (daily_clim) {
        # pull daily climate data
        min_date = out[[cad_list_ind]][[ident]]$mydata$date[1]
        max_date = out[[cad_list_ind]][[ident]]$mydata$date[length(out[[cad_list_ind]][[ident]]$mydata$date)]
        max_date = max_day_fun(max_date)

        if (dis_lut$identifier[patch_ind]==dis_lut$clim_ident[patch_ind]) {
          clim_key = master_keys[jj]
        } else {
          # patches that are smaller than the smallest GADM level are assigned the climate data of the smallest encompassing GADM patch
          clim_key = dis_lut$master_key[dis_lut$identifier==dis_lut$clim_ident[patch_ind]]
        }
        out[[cad_list_ind]][[ident]]$noaa_daily = noaa_clim[noaa_clim$master_key==clim_key, c("date", clim_cols)]
        # reduce to applicable dates
        date_ind = out[[cad_list_ind]][[ident]]$noaa_daily$date>=min_date & out[[cad_list_ind]][[ident]]$noaa_daily$date<=max_date
        out[[cad_list_ind]][[ident]]$noaa_daily = out[[cad_list_ind]][[ident]]$noaa_daily[date_ind, ]

        if (max_date>max(out[[cad_list_ind]][[ident]]$noaa_daily$date)) {
          # when incidence is ahead of climate data, pad climate data with NAs
          na.rows = data.frame(date=seq.Date(from=as.Date(max(out[[cad_list_ind]][[ident]]$noaa_daily$date))+1, to=max_date, by=1))
          for(kk in 1:length(clim_cols)) {
            na.rows[[clim_cols[kk]]] = NA
          }
          out[[cad_list_ind]][[ident]]$noaa_daily = rbind(out[[cad_list_ind]][[ident]]$noaa_daily, na.rows)
        }
      }

    }
  }
  cat("Complete\n")

  dbDisconnect(myDB)
  return(out)
}


get.StartEndDates <- function(country = "TH", model_year = 2010, disease = "dengue") {
  # Set Start/End of Season for a Given Country
  #
  # \code{build.cadence.vector} Given a country name and a model year determine the start/end month of a season,
  # and create a list which contains the months/weeks names, etc.
  # @param country The two letter (ISO2) country code
  # @param model_year Integer, the dengue season we want to model (using calendar year at start of season)
  # @param cadence String with cadence of mydata as inferred by the code 'month', 'week', 'day'
  # @return cadence.year.info.for.subset A list with the start and end month of the season, a vector with the
  # month numbers, the month names and the calendar year of each month.
  # @examples
  # build.cadence.vector(country = 'TH', model_year = 2010)
  #
  #
  if (disease=="dengue") {
    switch(as.character(country), TH = {
      # season follows calendar year and cadence is monthly
      start_date = as.Date(paste0(model_year,"-01-01"))
      end_date = as.Date(paste0(model_year,"-12-31"))
    }, BR = {
      # season starts Sep 1
      # Use EW to set date so they are compatible with both Monthly and Weekly cadence
      sep1_date = as.Date(paste0(model_year-1,"-09-01"))
      start_date = sep1_date
      end_date = sep1_date + 364 - wday(sep1_date)
      # start_week = Date2CDCweek(sep1_date)
      # start_date = CDCweek2date(CDCweek=start_week$CDCweek, year=start_week$CDCyear)
      # if (month(start_date)==9) {
      #   end_date = CDCweek2date(CDCweek=start_week$CDCweek-1, year=model_year)
      # } else {
      #   start_date = start_date+7
      #   end_date   = CDCweek2date(CDCweek=start_week$CDCweek, year=model_year)
      # }

    }, MX = {
      # season starts Apr 1
      start_date = as.Date(paste0(model_year,"-04-01"))
      end_date = as.Date(paste0(model_year+1,"-03-31"))
    }, SG = {
      # season follows calendar year.  dates should work for monthly or weekly cadence
      start_date = CDCweek2date(CDCweek=1, year=model_year)
      jan1date = as.Date(paste0(model_year, "-01-01"))
      if (start_date >= jan1date) {
        end_date = as.Date(paste0(model_year, "-12-31"))
      } else {
        end_date = start_date-1
        year(end_date) = year(end_date) + 1
        start_date = jan1date
      }
    },  {
      # for all other countries assume the season follows calendar year.  dates should work for monthly or weekly cadence
      # start_date = as.Date(paste0(model_year-1,"-12-29"))
      # end_date = as.Date(paste0(model_year,"-12-28"))
      start_date = CDCweek2date(CDCweek=1, year=model_year)
      jan1date = as.Date(paste0(model_year, "-01-01"))
      if (start_date >= jan1date) {
        end_date = as.Date(paste0(model_year, "-12-31"))
      } else {
        end_date = start_date-1
        year(end_date) = year(end_date) + 1
        start_date = jan1date
      }
    })
  } else if (disease=="flu") {
    if (country %in% c("AF", "SA", "OC")) {
      start_date = CDCweek2date(CDCweek=1, year=model_year)
      end_date  = CDCweek2date(CDCweek=52, year=model_year)
    } else {
      # Until we have more countries, assume season starts EW27 and ends EW26
      start_date = CDCweek2date(CDCweek=27, year=model_year)
      end_date  = CDCweek2date(CDCweek=26, year=model_year+1) + 1
    }
  }

  out = list(start_date=start_date, end_date=end_date)

  return(out)
}



get.mydata <- function(region = "world", nregion = NULL) {

    # Get Data for a DICE run
    # @section Warning:
    # \bold{This function is no longer supported}

    if (region == "world") {
        regionDF = get.world.mydata()
    } else if (region == "wonderland") {
        regionDF = make.wonderland.mydata(nregion = nregion)
    } else {
        regionDF = get.usa.mydata(region = region)

    }
    regionDF
}



get.usa.mydata <- function(region = "cdc") {

    # Get Data for a DICE run  on the continental US
    # @section Warning:
    # \bold{This function is no longer supported}

    mydata("usa-centrdiceData", package = "DICE")

    # remove the first two columns which are just some number - not clear even what
    usamydata = usamydata[, -c(1:2)]

    # these are the columns we now have: NAME ISO3 ID_1 NAME_1 ENGTYPE_1 long lat poplon poplat pop

    # for now: remove the district of columbia, alaska and hawaii
    usamydata = usamydata[-which(as.character(usamydata$NAME_1) == "Alaska"), ]  #c(\District of Columbia\',\'Alaska\',\'Hawaii\')'
    usamydata = usamydata[-which(as.character(usamydata$NAME_1) == "Hawaii"), ]

    # for now: remove the district of columbia, alaska and hawaii
    usamydata = usamydata[-which(as.character(usamydata$NAME_1) == "District of Columbia"), ]
    # now we are left with 48 states - CONUS

    statelist = as.character(usamydata$NAME_1)


    # same for all states

    iso3list = as.character(usamydata$ISO3)

    usapop = data.frame(usamydata$pop)

    rownames(usapop) = statelist

    # calculate the distances between the centroids - these are population weighted centroids

    usalonlat = data.frame(lon = usamydata$poplon, lat = usamydata$poplat)

    rownames(usalonlat) = as.character(statelist)

    rij = distance.matrix(as.matrix(usalonlat))

    mydata = list()
    # mydata$NAME=usamydata$NAME
    mydata$ISO3 = as.character(usamydata$ISO3)
    mydata$ID_1 = usamydata$ID_1
    mydata$NAME = as.character(usamydata$NAME_1)
    mydata$ENGTYPE_1 = usamydata$ENGTYPE_1
    mydata$long = usamydata$lon
    mydata$lat = usamydata$lat
    mydata$poplon = usamydata$poplon
    mydata$poplat = usamydata$poplat
    mydata$pop = usamydata$pop

    # now get the SH mydata on state or CDC level
    mydata("sh_region_state")

    usaSH = sh_region_state

    usaDF = list(mydata = mydata, rij = rij, usaSH = usaSH)

    usaDF

}



get.ili <- function(FrameName = NULL, ident = NULL, start.year = NULL, end.year = NULL, start.week = 27, end.week = 26, curYear = 2017,
    nperiods = 52) {

    # Extract incidence mydata from the mydata frame \code{diceData$FrameName}
    #
    # This function uses the input dates to pull incidence mydata from mydata frame \code{FrameName}.  If 'CDC' mydata from the current year is requested, the function will attempt to pull the mydata directly from CDC servers. \code{ident} is used to specify which columns are pulled from the mydata frame.
    #
    #  @param FrameName A string specifying the incidence mydata frame name: 'CDCili', 'GFTili', 'MiscCases', 'MiscILI', 'ZikaSusp', or 'ZikaConf'.
    #  @param ident A character vector specifying the identifiers (mydata frame column names) of the target population(s).  For example \code{ident=c('USA','USA.R9.CA')} specifies both the United States and the state of California.
    # @param start.year An integer - start year of the flu season
    # @param end.year  An integer  - end year of the flu season
    # @param start.week An integer - starting CDC week of the flu season (default is 27)
    # @param end.week   An integer - ending CDC week of the flu season (default is 26)
    # @param curYear if \code{start.year==curYear} and \code{FrameName=='CDCili'} the function will attempt to pull mydata directly from the CDC server.
    # @param nperiods Integer indicating the length of the returned incidence vectors.  If the available mydata is shorter than \code{nperiods}, the vectors will be padded with NAs.
    # @return a list with the incidence mydata $raw and the number of weeks of mydata $nperiodsData.
    # @examples
    # require(DICE)
    # myILI = get.ili(FrameName='GFTili',ident=c('USA','USA.R9.CA'),start.year=2014,end.year=2015,start.week=27,end.week=26,curYear=2015, nperiods=52)

    raw = list()
    raw[ident] = 0  # Initialize raw list names

    if (FrameName == "CDCili" && start.year == curYear) {
        # Attempt to download 'fit' ILI mydata from CDC server National mydata
        if ("USA" %in% ident) {
            CDCdat = get.cdc.CrYr(CDCreg = "national", sub_region = 1, start.year = start.year, end.year = end.year, start.week = start.week,
                end.week = end.week, ident = "USA", nperiods = nperiods)
            raw$USA = CDCdat$raw
            nperiodsData = CDCdat$nperiodsData
        }
        # Regional mydata
        RegIdent = paste0("USA.R", 1:10)
        ident.ind = ident %in% RegIdent
        if (any(ident.ind)) {
            CDCreg = "hhs"
            sub_region = as.numeric(sub("USA.R", "", ident[ident.ind]))
            CDCdat = get.cdc.CrYr(CDCreg = CDCreg, sub_region = sub_region, start.year = start.year, end.year = end.year, start.week = start.week,
                end.week = end.week, ident = ident[ident.ind], nperiods = nperiods)
            raw[ident.ind] = as.data.frame(CDCdat$raw)
            nperiodsData = CDCdat$nperiodsData
        }
        raw = as.data.frame(raw)
    } else {
        # if !(recent CDC) pull incidence from diceData
        start.ind = which(diceData[[FrameName]]$year == start.year & diceData[[FrameName]]$week == start.week)
        end.ind = which(diceData[[FrameName]]$year == end.year & diceData[[FrameName]]$week == end.week)
        raw = diceData[[FrameName]][start.ind:end.ind, ident]
        nperiodsData = length(start.ind:end.ind)
    }
    return(list(raw = raw, nperiodsData = nperiodsData))
}

get.meta <- function(FrameName = NULL, ident = NA, start.year = 2014) {
    # Get meta mydata for DICE incidence mydata
    #
    # Pulls the mydata units, epi/cases units, and a factor for converting from mydata(raw) to cases(epi).
    #
    #  @param FrameName A string specifying the incidence mydata frame name: 'CDCili', 'GFTili', 'MiscCases', 'MiscILI', 'ZikaSusp', or 'ZikaConf'.
    #  @param ident A character vector specifying the identifiers (mydata frame column names) of the target population(s).  For example \code{ident=c('USA','USA.R9.CA')} specifies both the United States and the state of California.
    # @param start.year An integer - start year of the flu season
    # @return meta a list with the mydata units $raw_units, the cases units $epi_units, and the conversion factor $factor.
    # @examples
    # require(DICE)
    # meta = get.meta(FrameName='GFTili',ident=c('USA','USA.R9.CA'),start.year=2014)

    meta = list()
    Types = unique(FrameName)

    ii = 1
    temp_ident = ident[FrameName == Types[ii]]
    # Pull incidence Meta-mydata
    meta$raw_units = diceData$MetaILI[[paste0(Types[ii], "_raw_units")]][temp_ident]
    meta$epi_units = diceData$MetaILI[[paste0(Types[ii], "_epi_units")]][temp_ident]
    meta$factor = diceData$MetaILI[[paste0(Types[ii], "_factor")]][diceData$MetaILI[[paste0(Types[ii], "_factor")]]$year == start.year,
        temp_ident]
    names(meta$factor) = temp_ident

    if (length(Types) > 1) {
        for (ii in 2:length(Types)) {
            temp_ident = ident[FrameName == Types[ii]]
            # Pull incidence Meta-mydata
            meta$raw_units = c(meta$raw_units, diceData$MetaILI[[paste0(Types[ii], "_raw_units")]][temp_ident])
            meta$epi_units = c(meta$epi_units, diceData$MetaILI[[paste0(Types[ii], "_epi_units")]][temp_ident])
            meta$factor[temp_ident] = diceData$MetaILI[[paste0(Types[ii], "_factor")]][diceData$MetaILI[[paste0(Types[ii], "_factor")]]$year ==
                start.year, temp_ident]
        }
    }
    # reorder each list then, convert to vector/matrix
    meta$raw_units = unlist(meta$raw_units[ident])
    meta$epi_units = unlist(meta$epi_units[ident])
    meta$factor = unlist(meta$factor[ident])
    return(meta)
}

get.cdc.CrYr <- function(CDCreg = "national", sub_region = 1:10, start.year = 2017, end.year = 2018, start.week = 27, end.week = 26,
    ident = "USA", nperiods = 52) {
    # Get current-year CDC incidence mydata
    #
    # Attempt to pull CDC mydata from the CDC server. If that fails, pull from diceData.
    #
    #  @param CDCreg String specifying USA ('national') or CDC Region ('hhs').
    #  @param sub_region Integer(s) indicating which of the 10 CDC regions are to be retrieved.
    # @param start.year An integer - start year of the flu season
    # @param end.year  An integer  - end year of the flu season
    # @param start.week An integer - starting CDC week of the flu season (default is 27)
    # @param end.week   An integer - ending CDC week of the flu season (default is 26)
    # @param ident A character vector specifying diceData-type identifiers.  In the case that the CDC serer cannot be contacted, \code{ident} is used to pull mydata from the DICE mydatabase (diceData).
    # @param nperiods Integer indicating the length of the returned incidence vectors.  If the available mydata is shorter than \code{nperiods}, the vectors will be padded with NAs.
    # @return a list with the incidence mydata $raw and the number of weeks of mydata $nperiodsData.
    # @examples
    # require(DICE)
    # myILI = get.cdc.CrYr(CDCreg='hhs',sub_region=c(1,9),ident=c('USA.R1','USA.R9'),start.year=2015,end.year=2016,start.week=27,end.week=26, nperiods=52)

    GoodPull = FALSE

    # Pull ili mydata from CDC server (make 5 attempts)
    n = 0
    while (n < 5) {
        CDCfit = try(get_flu_mydata(region = CDCreg, sub_region = sub_region, mydata_source = "ilinet", years = c(start.year - 1, start.year)),
            silent = TRUE)
        n = n + 1
        if (!is(CDCfit, "try-error")) {
            GoodPull = TRUE
            break
        } else Sys.sleep(5)
    }
    # If mydata download is unsuccessful, print warning and revert to saved mydata.
    if (!GoodPull) {
        print("WARNING: Download from CDC server was unsuccessful. Reverting to ili mydata in local DICE-package mydatabase.  DICE package mydata is likely less up-to-date than the CDC server.")

        # Pull fit-mydata from DICE mydataset
        start.ind = which(diceData$CDCili$year == start.year & diceData$CDCili$week == start.week)
        end.ind = which(diceData$CDCili$year == end.year & diceData$CDCili$week == end.week)
        # if the mydataset does not contain the full year, pull the weeks that are available
        if (length(end.ind) == 0) {
            end.ind = length(diceData$CDCili$week)
        }
        raw = diceData$CDCili[start.ind:end.ind, ident]
    } else {
        # Process CDC 'fit' download
        if (length(sub_region) > 1) {
            tempYear = CDCfit$YEAR[seq(from = 1, to = length(CDCfit$YEAR), by = length(sub_region))]
            tempWeek = CDCfit$WEEK[seq(from = 1, to = length(CDCfit$YEAR), by = length(sub_region))]
            start.ind = which(tempYear == start.year & tempWeek == start.week)
            end.ind = which(tempYear == (start.year + 1) & tempWeek == end.week)
            if (length(end.ind) == 0) {
                end.ind = length(tempYear)
            }

            # reshape CDC mydata to DICE-matrix form
            FitILI = array(mydata = CDCfit$`% WEIGHTED ILI`, dim = c(length(sub_region), length(CDCfit$YEAR)/length(sub_region)))[, start.ind:end.ind]
            raw = t(FitILI)

        } else {
            # only one fit region
            start.ind = which(CDCfit$YEAR == start.year & CDCfit$WEEK == start.week)
            end.ind = which(CDCfit$YEAR == (start.year + 1) & CDCfit$WEEK == end.week)
            if (length(end.ind) == 0) {
                end.ind = length(CDCfit$YEAR)
            }
            raw = CDCfit$`% WEIGHTED ILI`[start.ind:end.ind]
        }
    }  # end if GoodPull_fit

    # padweeks if necessary
    padweeks = nperiods - length(start.ind:end.ind)
    if (length(sub_region) > 1) {
        # pad with NAs to fill-out a complete year
        raw = rbind(as.matrix(raw), matrix(NA, padweeks, length(sub_region)))
    } else {
        # pad with zeros to fill-out a complete year
        raw = c(raw, rep(NA, padweeks))
    }

    nperiodsData = length(start.ind:end.ind)

    return(list(raw = raw, nperiodsData = nperiodsData))
}


get.FitLevelName <- function(iFit = 3) {

    # Given a spatial Fit level return the corresponding spatial name
    #
    # This function is called by \code{\link{get.DICE.data}}. It assigns the correct spatial name to the spatial level
    # @param iFit An integer (Currently support 2-4 for gft mydata and 2-3 for cdc mydata)
    # @return A character string with the corresponding spatial Name: earth, continent, country, hhs_region, state, county and city for iFit between 0 and 6
    # @examples
    # get.FitLevelName(iFit=3)
    #
    # get.FitLevelName(iFit=4)

    if (iFit == 0)
        fitLevelName = "earth"
    if (iFit == 1)
        fitLevelName = "continent"
    if (iFit == 2)
        fitLevelName = "country"
    if (iFit == 3)
        fitLevelName = "hhs_region"
    if (iFit == 4)
        fitLevelName = "state"
    if (iFit == 5)
        fitLevelName = "county"
    if (iFit == 6)
        fitLevelName = "city"

    return(fitLevelName)
}

get.ModelLevelName <- function(mod_level = 2, RegState = "usa", NAME_4 = NULL, NAME_5 = NULL) {

    # Given a spatial Model level and the corresponding RegState value return the corresponding spatial name
    # @param mod_level An integer 2-4 for gft mydata and 2-3 for cdc mydata
    # @param RegState A string, USA, or in case of mod_level=2, or region number of mod_level=3
    # @param NAME_4 Relevant for mod_level = 4 and up.  State name
    # @param NAME_5 Relevant for mod_level = 5 - county name
    if (mod_level == 2)
        modelLevelName = tolower(RegState)
    if (mod_level == 3)
        modelLevelName = paste("Region", RegState, sep = "")
    if (mod_level == 4)
        modelLevelName = tolower(RegState)
    if (mod_level == 5)
        modelLevelName = paste(tolower(NAME_4), "-", tolower(NAME_5), sep = "")
    return(modelLevelName)
}


calc.tmax <- function(loc = NULL) {

    # set the time of year where the force of infection is maximal
    #
    # Given the latitude of the country(ies) set the time of year when the force of infection is maximal.
    #   Currently this is not really used anywhere but will become useful if we model more than one season
    #   at a time and/or if we model both the northern and southern hemispheres.
    # @param loc Numeric - an array (or a single number) with the latitude of the region(s).
    # @return Numeric an array (or a single number) with the time of year of the maximum in the force of infection
    #   January 15 for the northern hemisphere and July 15 for the southern. The results appear as a day number.
    # @examples
    # calc.tmax{loc=45.}
    # calc.tmax{loc=c(-25,-35,40)}

    tmax.north = "2009-01-15 19:12:04 BST"
    tmax.south = "2009-07-15 19:12:04 BST"
    doy.north = strftime(tmax.north, format = "%j")
    doy.south = strftime(tmax.south, format = "%j")
    tmax.north = as.numeric(doy.north)
    tmax.south = as.numeric(doy.south)

    ncountry = length(loc)
    tmax <- rep(NA, ncountry)
    for (i in 1:ncountry) {
        if (loc[i] <= -20) {
            tmax[i] = tmax.south
        } else if (loc[i] >= 20) {
            tmax[i] = tmax.north
        } else {
            tmax[i] = 0  #The tropical region
        }
    }
    return(tmax)
}


runSinglePatch <- function(mydata = NULL, year = NULL, nperiodsFit = 52, model = 5, isingle = 1, nMCMC = 1e+05, nreal = 1, device = "pdf",
    prior = 0, Temp = 1, subDir = NULL, plot = 1, emcee = FALSE, nwalk = 10, iseed = NULL, Tg = NULL) {

    # Driver for running a single patch
    #
    # Run a single calculation on a specific country/region/state in a given state for given season using a given model and number of MCMC steps
    # @param mydata - dataframe with all the data for the run
    # @param year Integer - start year of the flu season
    # @param nperiodsFit Integer - Number of data points that will be fitted.  Default is to fit all the data.  This will be reset if nperiodsFit > nperiodsData or nperiodsFit = 0
    # @param model Integer - The model number, see manual for more details (1-5 are supported)
    # @param isingle Integer - should always be 1 since we are running a single patch
    # @param nMCMC Integer - number of steps/trials in the MCMC procedure
    # @param nreal Integer - number of chains
    # @param device  Either 'pdf' (default) 'png' or 'x11'. Accepts also an array with more than one device
    # @param subDir Name of output sub-directory where all plots and files will be written.  Default is 'output'
    # @param plot Character TRUE, FALSE or EXTRENAL (or 0, 1 or 2) Allows Users to write their own plotting routines
    # @param prior Integer - if greater than zero use a prior for the MCMC procedure - not used for now so set to zero
    # @param Temp Integer 1, 2, ..10 Relevant only if using a prior - controls the width of the Gaussian prior
    # @param emcee Logical - relevant only when running a sinlge region/patch selects between MCMC and EMCEE procedures
    # @param nwalk Integer - Number of walkers, relevant ONLY for an EMCEE procedure
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @param Tg - recovery time in days.  Reasonable range is 2-3 days
    # @return solution a list with the input and entire output of the run.
    # @examples
    # # Use model 3 and basic plotting routine
    #
    # runSinglePatch{mydata=mydata,year=2009,nperiodsFit=52,model=3,nMCMC=1e6,nreal=1,device='pdf',subDir='test',plot=1}
    #
    # # Use model 5, ggplot2 plotting and fit only the first 40 weeks
    #
    # runSinglePatch{mydata=mydata,year=2009,nperiodsFit=40,model=5,nMCMC=1e6,nreal=1,device='pdf',subDir='test',plot=2}

    ## Data Type: cdc or gft

    if (is.null(dataType))
        dataType = "misccases"

    ## Number of MCMC chains
    if (is.null(nreal))
        nreal = 1

    ## Number of steps/trials in each MCMC chain

    if (is.null(nMCMC))
        nMCMC = 10000

    ## Number of times the history of the MCMC chain is saved.

    nlines = round(nMCMC/100)
    nlines = min(10000, nlines)

    ## Start year of the flu season

    if (is.null(year))
        year = 2010

    ## Model Number
    if (is.null(model))
        model = 3

    ## Number of weeks of mydata that are fitted
    if (is.null(nperiodsFit))
        nperiodsFit = 52

    ## Recovery time in days - reasonable range is 2-3 days
    if (is.null(Tg))
        Tg = 3

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

    if (is.null(subDir))
        subDir = "output"

    if (is.null(plot))
        plot = 1

    if (is.null(emcee))
        emcee = FALSE

    ## We already have thee mydata but some cleaning is needed and we need to build some names

    ## For the SH it has both the state and region level mydata

    ## mydata <- get.DICE.data(dataType = dataType, mod_level = mod_level, fit_level = fit_level, year = year, model = model, nperiodsFit =
    ## nperiodsFit, RegState = RegState, isingle=isingle)

    ## augment the mydata with some more information that the user has chosen

    mydata$single = isingle

    ## TRUE or FALSE - default is FALSE

    mydata$prior = prior
    mydata$Temp = Temp


    if (nperiodsFit <= 0)
        nperiodsFit = mydata$nperiodsData
    if (nperiodsFit > mydata$nperiodsData)
        nperiodsFit = mydata$nperiodsData
    if (nperiodsFit < 10) {
        cat("\n\n WARNING: Resetting nperiodsFit to its minimal value of 10 \n\n")
        nperiodsFit = 10
    }
    if (nperiodsFit < 15) {
        cat("\n\n WARNING: Using a small Number of weeks in fitting procedure Results will have limited information \n\n")
    }

    mydata$nperiodsFit = nperiodsFit

    mydata$imodel = model

   # mydata$fitLevelName = get.FitLevelName(iFit = 5)

   # mydata$modelLevelName = get.ModelLevelName(mod_level = 5, NAME_4 = NAME_4, NAME_5 = NAME_5)

    dataType = tolower(dataType)
    mydata$dataType = dataType

    mydata$FY = paste(year, "-", (year + 1), sep = "")

    mydata$season = year

    mydata$dataName = paste(mydata$modelLevelName, "-", mydata$FY, "-", mydata$imodel, sep = "")

    ## Pack the information for the run

    opt.list <- set.opt.list(model = mydata$imodel, isingle = isingle)

    run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, nlines = nlines, device = device, subDir = subDir, plot = plot,
        nwalk = nwalk)

    ## Fitting a SINGLE region/patch

    cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
    cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$dataType), " mydata", "\n")
    cat("\n\n Fitting ", mydata$nperiodsFit, " weeks out of ", mydata$nperiodsData, " weeks of mydata", "\n\n")
    cat("\n\n Fitting ", toupper(mydata$model$name), " mydata", "\n")


    for (ireal in 1:nreal) {
        solution = fitOnePatch(mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal)
    }

    return(solution)

}

##
## Functions for checking that country/season is in the dengue mydatabase, getting the dengue mydata, subsetting, cadence etc.
##

get.dengue.mydata <- function(country = "TH", level = 2) {


    # Retrieve Dengue Data for Selected Country and Level
    #
    # \code{get.mydata} Retrieves Dengue mydata for country at requested spatial level
    # @param country Two letter (ISO2) country name
    # @param level Integer, requested spatial level of mydata (default = 2 country level)
    # @return mydata A list containing:
    # cases A mydata frame with year, month/week, ndays and cases information for country/states/regions \
    #
    # the country/states/regions full name(s) (state_names) \
    #
    # the country/states/regions abbreavation(s) (state_id) \
    #
    # the country/state/regions population (pop) \
    #
    # the country/states/regions population weighted long/lat levels (lon/lat)

    # @examples
    # get.mydata(country = 'TH', level = 2)
    #

    # filename = '~/Dropbox/LEPR03/mydata/dengue/DICE_dengue.Rds' DENVdat = readRDS(file = filename) level depends from country to country
    attr_ind = DENVdat$attr$ABBV_2 == country & DENVdat$attr$level == level
    identifiers = DENVdat$attr$identifier[attr_ind]

    ## Check to see that the spatial level exists for this country

    if (length(identifiers) == 0) {
        cat("\n The User Selected Level Does Not Exist For This Country \n\n")
        cat("\n Available Levels are: \n")
        cat("\n Br - levels 2,3,4 \n")
        cat("\n MX - levels 2 and 3\n")
        cat("\n TH - levels 2,3,4,5 \n")
        cat("\n Code will Stop! \n")
        quit(save = "no")
    }

    # This depends on the level the user is asking for BR - level 4 TH - level 5 MX - level 3 LK - level 3 SG - level 2 TW - level 2
    # state_names is used somehwat loosely here..

    my_level = level
    switch(as.character(my_level), `2` = {
        state_names <- DENVdat$attr$NAME_2[attr_ind]
    }, `3` = {
        state_names <- DENVdat$attr$NAME_3[attr_ind]
    }, `4` = {
        state_names <- DENVdat$attr$NAME_4[attr_ind]
    }, `5` = {
        state_names <- DENVdat$attr$NAME_5[attr_ind]
    }, {
        state_names <- DENVdat$attr$NAME_2[attr_ind]
    })

    # population-weighted centroid lat/lon
    my_lat = DENVdat$attr[attr_ind, "sedac_lat"]
    my_lon = DENVdat$attr[attr_ind, "sedac_lon"]
    my_pop = DENVdat$attr[attr_ind, "sedac_pop"]
    my_pop = round(my_pop)


    ## will need to be able to choose between weekly and monthly here..

    ident_monthly = identifiers[identifiers %in% names(DENVdat$DENG_monthly)]
    ident_weekly = identifiers[identifiers %in% names(DENVdat$DENG_weekly)]
    if (length(ident_monthly) > 0) {
        cases = DENVdat$DENG_monthly[, c("year", "month", "ndays", identifiers)]
    } else {
        cases = DENVdat$DENG_weekly[, c("year", "week", "ndays", identifiers)]
    }

    if (length(identifiers) > 1) {
        NA_ind = apply(is.na(cases[, identifiers]), 1, all)
    } else {
        NA_ind = is.na(cases[, identifiers])
    }

    cases = cases[!NA_ind, ]

    noTPts = dim(cases)[1]

    nstates = length(state_names)

    ## We can generalize this to all countries if we want

    if (level > 2 && country == "MX") {

        ikeep = NULL
        istate_keep = NULL
        for (istate in 1:nstates) {
            id = identifiers[istate]
            nzero = sum(cases[, id] == 0, na.rm = TRUE)
            if (nzero < (noTPts/2) || length(nzero) == 0) {
                ikeep = cbind(ikeep, id)
                istate_keep = cbind(istate_keep, istate)
            } else {

            }
        }
        nstates = length(ikeep)
        cases = cases[, c("year", "month", "ndays", ikeep)]
        identifiers = identifiers[istate_keep]
        my_lon = my_lon[istate_keep]
        my_lat = my_lat[istate_keep]
        my_pop = my_pop[istate_keep]
        state_names = state_names[istate_keep]


    }



    ## Now remove states with less than half the mydata missing This does not happen on the national level

    mydata = list(state_names = state_names, pop = my_pop, cases = cases, lat = my_lat, lon = my_lon, state_id = identifiers)
    return(mydata)

}

get.dengue.years <- function(cases = NULL, cadence = NULL) {

    # Year and Month Information for the Dengue Data
    #
    # \code{get.years} Determines the start/end year and month/week for the dengue cases mydata frame, the unique
    # years and the unique number of years
    # @param cases  A mydata frame with year, month/week, ndays and cases information
    # @param cadence A string with the cadence, month, week, day. It is inferred from the mydata by the code
    # @return years.info A list with start/end information for the year and month (or week/day), a vector of unique years
    # and the number of unique years
    # @examples
    # get.years(cases = cases, cadence = cadence)
    #
    #
    nlines = dim(cases)[1]

    start.year = cases[1, "year"]

    end.year = cases[nlines, "year"]

    if (cadence == "month") {
        start.month = cases[1, "month"]
        end.month = cases[nlines, "month"]
        start.week = end.week = NA
    } else {
        start.week = cases[1, "week"]
        end.week = cases[nlines, "week"]
        start.month = end.month = NA
    }


    years = unique(cases[, "year"])

    nyears = length(years)

    years.info = list(start.year = start.year, end.year = end.year, start.month = start.month, end.month = end.month, start.week = start.week,
        end.week = end.week, years = years, nyears = nyears)

    return(years.info)
}

check.dengue.country <- function(country = NULL) {
    # Check that the user has selected a country that is in the mydata set and can be modeled
    # @param country Two letter (ISO2) country name
    # @return err = 0 if country is in the mydata base
    # @examples check.country(country = country)

    country.vec = c("BR", "MX", "TH", "SG")
    if (country %in% country.vec) {
        return(err = 0)
    } else {
        cat("\n\n vectorDICE Currently Supports Only the Following Countries: \n")
        cat("\n BR (Brazil)\n")
        cat("\n MX (Mexico)\n")
        cat("\n TH (Thailand)\n")
        cat("\n SG (Singapore)\n")
        cat("\n Please Select One of the Above and Resubmit \n\n")
        cat("\n Code Will Stop!\n")
        quit(save = "no")

    }

}


build.dengue.cadence.vector <- function(country = "TH", model_year = 2010, cadence = "month") {

    # Set Start/End of Season for a Given Country
    #
    # \code{build.cadence.vector} Given a country name and a model year determine the start/end month of a season,
    # and create a list which contains the months/weeks names, etc.
    # @param country The two letter (ISO2) country code
    # @param model_year Integer, the dengue season we want to model (using calendar year at start of season)
    # @param cadence String with cadence of mydata as inferred by the code 'month', 'week', 'day'
    # @return cadence.year.info.for.subset A list with the start and end month of the season, a vector with the
    # month numbers, the month names and the calendar year of each month.
    # @examples
    # build.cadence.vector(country = 'TH', model_year = 2010)
    #
    #


    if (cadence == "month") {


        month_per_year = 12

        switch(as.character(country), TH = {
            cadence.start = 1
            month.end = 12
            cadence.vec = seq(cadence.start, month.end, by = 1)
            cadence.names = paste0("month", cadence.vec)
            year.vec = rep(model_year, length(cadence.vec))
        }, BR = {
            cadence.start = 9
            month.end = 8
            cadence.vec = cadence.start:month_per_year
            cadence.vec = c(cadence.vec, 1:month.end)
            cadence.names = paste0("month", cadence.vec)
            year.vec = rep((model_year - 1), length(cadence.start:month_per_year))
            year.vec = c(year.vec, rep(model_year, length(1:month.end)))
        }, MX = {
            cadence.start = 4
            month.end = 3
            cadence.vec = cadence.start:month_per_year
            cadence.vec = c(cadence.vec, 1:month.end)
            cadence.names = paste0("month", cadence.vec)
            year.vec = rep(model_year, length(cadence.start:month_per_year))
            year.vec = c(year.vec, rep((model_year + 1), length(1:month.end)))
        }, {
            cadence.start = 1
            month.end = 12
            cadence.vec = seq(cadence.start, month.end, by = 1)
            cadence.names = paste0("month", cadence.vec)
            year.vec = rep(model_year, length(cadence.vec))
        })


    } else if (cadence == "week") {
        week_per_year = 52  # Will need to see what to do if we have 53 weeks..
        switch(as.character(country), SG = {
            week.start = 1
            week.end = week_per_year
            cadence.vec = seq(week.start, week.end, by = 1)
            cadence.names = paste0(cadence, cadence.vec)
            year.vec = rep(model_year, length(cadence.vec))
        }, {
            week.start = 1
            week.end = week_per_year
            cadence.vec = seq(week.start, week.end, by = 1)
            cadence.names = paste0(cadence, cadence.vec)
            year.vec = rep(model_year, length(cadence.vec))

        })



    } else {
        cat("\n\n Currently Supporting only Monthly and Weekly Cadence \n\n")
        cat("\n\n Code Will Stop!! \n\n")
        quit(save = "no")
    }


    cadence.year.info.for.subset = list(cadence.vec = cadence.vec, cadence.names = cadence.names, year.vec = year.vec)

    return(cadence.year.info.for.subset)

}



check.dengue.season <- function(years.info = NULL, cadence.year = NULL) {
	#  Check that Selected Season is in the mydatabase
	#
# @param year.info A list with the start/end year and month for the incidence.
# @param cadence.year A list with the Year and Month/Week information for the Season the User wants to model
# @return err = 0 If season is in the mydatabase
# @examples
# build.cadence.vector(years.info = years.info, cadence.year = cadence.year)
#
#

	cadence.vec = cadence.year$cadence.vec
	year.vec = cadence.year$year.vec

	ntps = length(cadence.vec)

	end.year = years.info$end.year
	end.month = years.info$end.month

	if (end.year == year.vec[ntps] && end.month < cadence.vec[ntps] | end.year < year.vec[ntps]) {
		cat("\n vectorDICE Does Not Have Data For This Country,Level,Season Selection\n\n")
		cat("\n Please Choose a Previous Season \n")
		cat("\n Code Will Stop!\n")
		quit(save = "no")
	} else {
		return(err = 0)
	}

}



mydata.dengue.augment_old <- function(da = 0, nfit = NULL, cases = NULL, epi = NULL, state_id = NULL, cadence.year = NULL, cadence = "month",
    corrMAT = NULL, distMAT = null, deng.null = NULL) {
    # Data Augmentaton for SEIR Model
    #
    # \code{mydata.augment} Allows the User to augment the dengue time series beyond the \code{nfit} mydata points
    # using either the historic NULL model \code{da = 0}, or the most similar past season \code{da=1}.
    # If \code{da=-1}, the mydata is not augmented.
    # @param da Integer 1 - mydata augmentation using the most similar past season
    #  0 mydata augmentation using the historic NULL model
    #  -1 NO mydata augmentation
    # @param cases A mydata frame with year, month, ndays and cases information
    # @param epi A mydata frame with year, month/week, ndays and cases information, only for the season chosen
    # for modeling
    # @param nfit An integer - the number of mydata points used for the curent dengue season
    # @param state_id The abbreviation of the state/region
    # @param cadence.year A list with year and month/week information for the time series that we are modeling
    # @param cadence String, the cadence of the mydata: month, week
    # @param corrMAT A 2D array with the Pearson corrletion values for all available dengue seasons
    # @param distMAT A 2D array with the Euler distance values for all available dengue seasons
    # @param deng.null the Historic NULL model - average monthly number of cases
    # @return epi An updated mydata frame with the original and augmented mydata (if da=0,1)
    # @examples
    # mydata.augment(da = 0, cases = cases, state_id = state_id, nfit = nfit, cadence.year=cadence.year,
    # cadence = cadence, tables=tables, corrMAT = corrMAT, distMAT = distMAT, deng.null=deng.null)
    #
    if (is.null(nfit) | is.null(cases) | is.null(epi))
        return


    if (is.null(da))
        da = 0

    da = FALSE

    if (da >= 0)
        da = TRUE

    epi.state = epi[, state_id]

    nmydata = length(epi.state)

    nda = nmydata - nfit

    if (nda == 0)
        da = FALSE


    if (da) {
        cat("\n Entered da == TRUE Option \n\n")
        epi.save = epi
        epi.state.save = epi.state
        my.year = cadence.year$year.vec[1]
        iyear = which(rownames(corrMAT) == my.year)

        if (all(is.na(corrMAT[iyear, ])) || da == 0) {
            corr.max = -1
        } else {
            corr.max = max(corrMAT[iyear, ], na.rm = TRUE)
        }


        if (corr.max > 0) {
            jyear = which.max(corrMAT[iyear, ])
            jyear = colnames(corrMAT)[jyear]

            jyear = as.numeric(jyear)

            futur = subset.sql.mydata(year.start = jyear, cadence.start = cadence.year$cadence.vec[1], cadence = cadence, ntps = length(cadence.year$cadence.vec),
                cases = cases)
            futur = futur[, state_id]
            shft = epi.state[nfit] - futur[nfit]

            x = futur[1:nfit]
            y = epi.state[1:nfit]

            b1 = sum((x - mean(x)) * (y - mean(y)))/sum((x - mean(x)) * (x - mean(x)))

            b0 = mean(y) - b1 * mean(x)

            z = b0 + b1 * futur

            cat("\n Using Pearson Correlation to Choose Season for Data Augmentation \n")
            cat("\n Season Selected is: ", jyear, "\n")


        } else {

            cat(" Using Historic NULL model to Augment Data \n")

            futur = deng.null[, state_id]
            x = futur[1:nfit]
            y = epi.state[1:nfit]

            shft = epi.state[nfit] - futur[nfit]

            b1 = sum((x - mean(x)) * (y - mean(y)))/sum((x - mean(x)) * (x - mean(x)))

            b0 = mean(y) - b1 * mean(x)

            z = b0 + b1 * futur



        }

        epi.new = epi.state
        wght = rep(1, nmydata)

        if (nfit <= 4 || all(is.na(corrMAT[iyear, ]))) {
            epi.new[(nfit + 1):nmydata] = futur[(nfit + 1):nmydata] + shft
        } else {
            epi.new[(nfit + 1):nmydata] = z[(nfit + 1):nmydata]
            wght[(nfit + 1):nmydata] = cor(x[1:nfit], y[1:nfit])
        }

        epi.new[(nfit + 1):nmydata] = round(epi.new[(nfit + 1):nmydata])

        for (i in (nfit + 1):nmydata) {
            epi.new[i] = max(1, epi.new[i])
        }

    } else {

        cat("\n\n NO Data Augmentation \n\n")
        epi.save = epi
        epi.state.save = epi.state
        epi.new = epi.state
        wght = rep(0, nmydata)
        wght[1:nfit] = 1

    }

    return(list(nmydata = nmydata, nfit = nfit, epi.new = epi.new, wght = wght, epi.state.save = epi.state.save, epi.save = epi.save))

}


saveDengueRData <- function(mydata = NULL, input = NULL, country = NULL, mod_level = NULL, nfit = NULL, tab.fit=NULL,tab.model=NULL,  ireal = 1) {

    # Save a Binary RData file for a Dengue Run
    #
    # Saves a binary RData file with all the input parameters and the data for this run.
    #   Each MCMC chain saves its own RData file and the file name includes the MCMC chain number
    #   as the last digit before the .RData extension of the file name.
    # @param mydata a list with all the subset mydata for the run
    # @param model An integer, the model number for the force of infection
    # @param input a list of all other input  needed for the run
    # @param country The Model country or region abbrevaiation
    # @param mod_level Numeric - The spatial model level
    # @param nfit Numeric - the number of data points used in the fit
    # @param ireal Integer indicating the MCMC chain number
    # @return err   \eqn{err =0} if successful.
    # @examples
    # saveDengueRData(mydata=mydata,input=input, country = mydata$model$name, mod_level = mod_level, nfit = nfit,
    # tab.fit = tab.fit, tab.model = tab.model, ireal=1)

    model = mydata$imodel
    subDir = input$run.list$subDir
    year.vec = cadence.year$year.vec
	Tg = mydata$Tg

    # check to see if 'mydata' sub-directory exists, if not create it
    if (is.null(subDir))
        subDir = "output"

    err = makeDir(subDir = subDir)
    filename = paste0(subDir, "mcmc_", country, "_year_", year.vec[1], "_level_", mod_level, "_Tg_", as.character(Tg),"_nfit_",nfit,"_model_",model,".RData")
    cat("\n Saving Input Run Information to: ", filename, "\n")

    save(mydata, input, tab.model, tab.fit, nfit, mod_level, file = filename)

    err = 0
    return(err)

}


run.vectorDICE <- function(RegState = "BR", mod_level = 2, fit_level = 3, year = 2010, Tg = 10, nMCMC = 1e+06, nfit = 12, model = 4, isingle = 0, da = 0,
    plot = TRUE, nreal = 1, subDir = NULL, prior = 0, Temp = 1.0, epi_model = SEIR, sql_db = TRUE, disease ='dengue', device = 'pdf', movie = NULL, emcee = NULL, nwalk = NULL, iseed = NULL) {


    # Driver for the Dengue SEIR/SARIMA Fitting Package
    #
    # \code{run.vectorDICE} Runs eight SARIMA models and a mechanistic compartmental S-E-I-R model
    # fitting and forecasting dengue monthly (or weekly) incidence profiles.
    # @param RegState Two letter (ISO2) RegState name
    # @param level Integer, requested spatial level of mydata (default = 2 RegState level)
    # @param year Integer, Dengue season User chooses to fit and/or forecast
    # @param epi_model String (SIR/SEIR) or integer (1/2)
    # @param Tg Numeric, recovery rate in days used in the S-E-I-R Eqs.
    # @param nMCMC, Integer, number of steps used in the MCMC procedure (SEIR model, default =1e6)
    # @param nfit Integer, number of mydata points in season we want to fit.
    # @param da Inetegr, type of mydata augmentation for calculation
    # @param Temp - Temperature for MCMC procedure, default is 1.
    # If set to NULL the code will set this to: using the historic null model (early in th season),
    # most similar season (middle of season), and finally no mydata augmentation (late in the season).
    # @param new Logical TRUE (default) / FALSE. Is this a new calculation or are we loading a past calculation
    # and plotting/analyzing the fit/forecast
    # @param plot Logical TRUE (default) / FALSE.  Plot the dengue incidence and the results (or not)
    # @return mydata  A list containing:
    # tables - the tables with both the observed and fit/forecast for the SEIR and SARIMA models.
    # The tables contain detailed analysis of the errors of the models: mean absolute error (mae),
    # mean relative error (mre), and mean error relative to the historic NULL model (rel).
    # profiles - 10,000 randomly selected profiles from the MCMC fitting prcedure
    # @examples
    # Fit the 2010 Dengue season in Thailand at the RegState level
    # run.vectorDICE(RegState = 'TH', level = 2, model_year = 2010, Tg = 10.,
    #nMCMC = 1e+06, nfit.vec = 12, da = NULL, plot = TRUE, Temp =1)
    #

 	nperiodsFit = nfit

  if (is.null(Tg)) {
  	Tg = 3.0
  	if (disease == "dengue")
  		Tg = 10.0
  }


   if (is.null(epi_model)) {
   	epi_model = 1
   }

   if (tolower(epi_model) == "sir") {
   	epi_model = 1
   }

   if (tolower(epi_model) == "seir") {
   	epi_model = 2
   }

   if (epi_model != 1 & epi_model != 2) {
   	epi_model = 1
   }

    ## Number of MCMC chains
    if (is.null(nreal))
        nreal = 1

    ## Number of steps/trials in each MCMC chain

    if (is.null(nMCMC))
        nMCMC = 10000

    if (is.null(historic))
        historic = FALSE

   ## Number of steps/trials in each MCMC chain

   if (is.null(nMCMC))
   	nMCMC = 10000


   ## Number of times the history of the MCMC chain is saved.

   nlines = round(nMCMC/100)
   nlines = min(10000, nlines)
   nlines = max(100, nlines)

   ithin = nMCMC / nlines

    if (is.null(mod_level))
        mod_level = 2

    if (is.null(fit_level))
        fit_level = 3
    ## Start year of the flu season

    if (is.null(year))
        year = 2010

    ## Model Number
    if (is.null(model))
        model = 4

    ## Number of periods of mydata that are fitted
    if (is.null(nperiodsFit))
        nperiodsFit = 12

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

    if (is.null(RegState)) {
        if (mod_level == 2)
            RegState = "usa"
            if (disease == 'dengue') RegState = 'BR'
        if (mod_level == 3)
            RegState = 1
    }

   if (is.null(subDir))
   	subDir = paste0(RegState, "_", year, "_level_",mod_level,'_prior_',prior,'_Temp_',Temp,'_da_',da,"/")

   if (!dir.exists(subDir)) {
   	dir.create(subDir)
   }


    if (is.null(plot))
        plot = 1
    if (is.null(movie))
        movie = 1

    if (is.null(emcee))
        emcee = FALSE
    if (is.null(nwalk))
        nwalk = 1


  if (model == 2 || model == 3) {
  	cat("\nSchool Vacation Models can not be used for Dengue Data \n
  	Resetting the model to 4 - Fixed Force of Infection !!\n\n")
  	model = 4
  }




   ## get the mydata - all years and for the RegState and level determined by the User,


	if (disease == "dengue") {


		complete_mydata = get.DICE.data(mod_level = mod_level, fit_level = fit_level, year = year, nperiodsFit = nfit, model = model, isingle = isingle,
			disease = disease, RegState = RegState)


		mydata = complete_mydata$mydata # season for modeling

		all_years_epi = complete_mydata$all_years_epi

		nmydata = length(all_years_epi$year)

		all_years_epi$FY = paste0(all_years_epi$year[1], "-", all_years_epi$year[nmydata])

		all_years_epi$mydataName = paste0(mydata$model$name, "-", all_years_epi$FY)


	}

 	## SIR/ or SEIR case insensetive

  	mydata$epi_model = epi_model

    mydata$isingle = isingle

    ## TRUE or FALSE - default is FALSE

    if (disease == 'dengue') {
    	mydata$prior = 0 # We do not have any priors in the case od dengue
    }  else {
    	mydata$prior = prior
    }

    mydata$Temp = Temp

    mydata$fit_level = fit_level

    mydata$mod_level = mod_level

	mydata$subDir = subDir

	mydata$da = da

   ## determine if incidence is weekly or monthly

   if(mydata$cadence == 'Weekly') {
   	cadence = 'week'
   } else if (mydata$cadence == 'Monthly') {
   	cadence = 'month'
   } else if (mydata$cadence == 'Daily') {
   	cadence = 'day'
   } else {
   	cadence = 'Unknown'
   	cat('\n\n Cadence of Data is Unknown, Code will Stop !!\n\n')
   	q()
   }


	if (disease == "dengue") {

		if (is.null(mydata$model$school)) {
			nperiods = mydata$nperiods
			mydata$model$school = rep(0, nperiods)
		}

		if (is.null(mydata$fit$school)) {
			nperiods = mydata$nperiods
			nregions = mydata$fit$nregions
			mydata$fit$school = mydata$fit$sh
			mydata$fit$school[1:nperiods, 1:nregions] = 0
		}


		err <- plotDisease(mydata = mydata, all_years_epi = all_years_epi, PDF = TRUE)

		cadence.year <- build.cadence.vector.sql(mydata = mydata)

	}

  ## Pack the information for the run, isingle=1 uncoupled run

     par_names <- set.param.list()

      run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, nlines = nlines, device = device, subDir = subDir, plot = plot, movie = movie, nwalk = nwalk)


	opt.list <- set.opt.list(model = mydata$imodel, isingle = isingle)

	   # ## Fitting a SINGLE region/patch

    if (mydata$fit$nregions == 1) {

        cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
        cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$disease), " mydata", "\n")
        cat("\n\n Fitting ", mydata$nperiodsFit, " periods out of ", mydata$nperiodsData, " periods of mydata", "\n\n")
        cat("\n\n Fitting ", toupper(mydata$modelLevelName), " mydata", "\n")

        for (ireal in 1:nreal) {
            solution = fitOnePatchDengue(mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed)
        }

        return(solution)
    }

    ## ALL OTHER CASES - number of regions/patches of fit_level > 1

    cat("\n\n ********  DICE: Dynamics of Interacting Community Epidemics  ********\n\n")
    cat("\n\n Fitting the ", mydata$FY, " Season using ", toupper(mydata$disease), " mydata", "\n")
    cat("\n\n Fitting ", mydata$nperiodsFit, " periods out of ", mydata$nperiodsData, " periods of mydata", "\n\n")
    cat("\n\n Fitting ", toupper(mydata$mode$name), " using ", toupper(mydata$fit$name), " mydata", "\n")

    ## UNCOUPLED CASE

    if (isingle == TRUE) {
        for (ireal in 1:nreal) {

            solution = fitSingleDengue(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal, iseed = iseed)

        }
        ## COUPLED CASE

    } else {

        for (ireal in 1:nreal) {
            solution = fitMultiDengue(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, opt.list = opt.list, ireal = ireal, iseed = iseed)
        }
    }

	return(solution)
}


###
### This file holds ALL the R functions used for the EMCEE procedure -  THEY ARE NO LONGER SUPPORTED
###


setup.emcee <- function(nMCMC = 1e+05, nlines = 10000, reals = 1, nparam = NULL, nwalk = 10) {

    # Setup for an EMCEE Procedure
    #
    # @section Warning:
    # \bold{This function is no longer supported}


    ithin <- round(nMCMC/nlines)
    if (ithin < 1)
        ithin = 1

    # accept <- 0.0 # accept rate for MCMC procedure

    tab <- array(data = 0, dim = c(nlines, (nparam + 1), nwalk))

    list(ithin = ithin, nMCMC = nMCMC, tab = tab, accept = 0)
}

setup.model.emcee <- function(n = 10, nwalk = 10, nperiods = 52, cadence = 7, tps = NULL, t0 = 1, nR0 = 1) {

    # Setup all the min/max and initial value parameters for an EMCEE \code{DICE} run
    #
    # @section Warning:
    # \bold{This function is no longer supported}

    parmin = rep(tps[1], n)
    parmax = rep(tps[(nperiods/2)], n)
    iparam = n
    vecopt = rep("t0", n)

    parmin = append(parmin, rep(1e-06, n))
    parmax = append(parmax, rep(1, n))

    nparam = length(par)
    vecopt = append(vecopt, rep("pC", n))
    iparam = append(iparam, n)

    # LET'S ADD R0 TO THE FITTING PROCEDURE

    parmin = append(parmin, rep(1.001, nR0))
    parmax = append(parmax, rep(3, nR0))
    vecopt = append(vecopt, rep("R0", nR0))
    iparam = append(iparam, nR0)

    # now seed
    parmin = append(parmin, 10)
    parmax = append(parmax, 1e+05)
    vecopt = append(vecopt, "seed")
    iparam = append(iparam, 1)

    # now the parameters a, b, c and d but we shall start with 'd' only
    parmin = append(parmin, 0.5)
    parmax = append(parmax, 1)
    vecopt = append(vecopt, "d")
    iparam = append(iparam, 1)

    # now add c - the power of the distance
    parmin = append(parmin, 0.5)
    parmax = append(parmax, 8)
    vecopt = append(vecopt, "c")
    iparam = append(iparam, 1)

    # and the power of the populations - apow for n_i and bpow for n_j give them the same range and the same initial guess of 1.0

    parmin = append(parmin, rep(0.1, 2))
    parmax = append(parmax, rep(5, 2))

    vecopt = append(vecopt, c("a", "b"))
    iparam = append(iparam, c(1, 1))

    nparam = length(parmin)

    par = matrix(data = 0, nrow = nparam, ncol = nwalk)

    for (iwalk in 1:nwalk) {
        par[1:n, iwalk] = runif(n, t0[1], (t0[1] + cadence * 4))  #t0
        par[(n + 1):(2 * n), iwalk] = runif(n, 0.001, 0.01)  #pC
        par[(2 * n + 1):(2 * n + nR0), iwalk] = runif(nR0, 1.1, 1.5)  #R0
        par[(2 * n + nR0 + 1), iwalk] = runif(1, min = 100, max = 1000)  #seed
        par[(2 * n + nR0 + 2), iwalk] = runif(1, 0.9, 1)  #diag
        par[(2 * n + nR0 + 3), iwalk] = runif(1, 2, 4)  #cpow
        par[(2 * n + nR0 + 4), iwalk] = runif(1, 1, 2)  #apow
        par[(2 * n + nR0 + 5), iwalk] = runif(1, 1, 2)  #bpow

    }


    emcee_list = list(parmin = parmin, parmax = parmax, iparam = iparam, par = par, vecopt = vecopt)
    return(emcee_list)
}

fitEMCEE <- function(region = "cdc", nwalk = 10, regionDF = NULL, e_nonflu = 1, iseed = NULL, w = 0.01, Tg = 3, pC = 0.005, one_R0 = FALSE,
    R0 = 1.4, R1 = 0, nperiods = 52, cadence = 7, dt = 0.2, nMCMC = 1000, nlines = 1000, device = "x11", restart = FALSE) {


    # Driver function for an EMCEE run
    #
    # @section Warning:
    # \bold{This function is no longer supported}


    country = regionDF$data
    rij = regionDF$rij
    # make sure the populations are integers
    country$pop = round(country$pop)
    ncountry <- length(country$pop)
    pop = country$pop

    country$pop = round(country$pop)
    ncountry <- length(country$pop)

    # in case there are NA must replace them with zero

    regionDF$epi[is.na(regionDF$epi)] = 0

    regionDF$gamaepi = gammaEpi(epi = regionDF$epi)

    nR0 = ncountry
    if (one_R0 == TRUE)
        nR0 = 1

    weeks = regionDF$weeks
    nperiods = regionDF$nperiods
    nperiodsFit = regionDF$nperiodsFit

    out <- setupODEforCDC(region = region, regionDF = regionDF, R0 = R0, Tg = Tg, pC = pC, iseed = iseed, nperiods = nperiods,
        cadence = cadence, dt = dt)

    tps = out$tps

    R0 = out$R0
    Tg = out$Tg
    pC = out$pC

    # The index of the region where ILI started

    ip = epiStart$ip

    pop = regionDF$data$pop

    n = length(pop)

    mij = regionDF$mij

    rij = regionDF$rij

    rtn = out$rtn

    sh = as.matrix(regionDF$sh)

    school = matrix(0, nr = nperiods, nc = ncountry)

    time = out$t0

    t0 = rep(out$t0, n)

    # t0 for regions 1-n will be fitted let's limit the max to be half the season

    if (restart == FALSE) {

        setup = setup.model.emcee(n = n, nwalk = nwalk, nperiods = nperiods, t0 = t0[1], cadence = cadence, tps = tps, nR0 = nR0)

        parmin = setup$parmin
        parmax = setup$parmax
        iparam = setup$iparam
        par = setup$par
        vecopt = setup$vecopt

    } else {

        cat("\n Restart of Fitting Loading data/parameters from: output/restart.RData\n\n")
        mainDir = getwd()
        outDir = "output"
        subDir = paste(mainDir, "/", outDir, sep = "")
        fileRestart = paste("restart-emcee-", regionDF$FY, ".RData", sep = "")
        filePath = paste(subDir, "/", fileRestart, sep = "")
        if (file.exists(filePath)) {
            load(filePath)
        } else {
            cat("\n\n No Restart File Found \n\n")
            cat("\n Set restart = FALSE and rerun \n")
            q(save = "no")
        }

        nwalk = input$nwalk

        e_nonflu = input$e_nonflu
        Tg = rep(input$Tg, n)
        iseed = set.iseed()
        w = rep(input$w, n)
        cadence = input$cadence
        dt = input$dt
        one_R0 = input$one_R0
        par = solution$par
        parmin = solution$parmin
        parmax = solution$parmax
        vecopt = solution$vecopt
        iparam = solution$iparam

        nR0 = ncountry
        if (one_R0 == TRUE)
            nR0 = 1

    }

    nparam = length(parmin)
    nopt = length(vecopt)

    # all of this will be moved to a function once it is working
    model = list()  #setup.Epimodel(param=t0,vecopt=vecopt)
    model$ilog = rep(0, nopt)
    model$vecopt = vecopt

    logbase = 10  #use log base 10 when needed
    logvec = rep(1, nparam)
    # logvec <- rep(0,nparam) #1 use log uniform sampling , 0 use uniform sampling - must use this for any parameter that can be
    # non-positive logvec[which(pmin <= 0)] <- 0

    imask = rep(1, nparam)

    # names(imask) <- vecnames imask[vecopt] <- 1 #vecopt holds the sub-list of what we want to optimize

    mcmc = setup.emcee(nMCMC = nMCMC, nlines = nlines, nparam = length(par), nwalk = nwalk)

    iverbose = 1  # TRUE, 0 == FALSE

    cat("Fitting ", nperiodsFit, " out of ", nperiods, " of data", "\n\n")

    solution = .Fortran("emcee", n = as.integer(n), epi = as.double(regionDF$epi), gamaepi = as.double(regionDF$gamaepi), sh = as.double(sh),
        school = as.double(school), nparam = as.integer(nparam), nwalk = as.integer(nwalk), iparam = as.integer(iparam), par = as.double(par),
        parmin = as.double(parmin), parmax = as.double(parmax), iverbose = as.integer(iverbose), imask = as.integer(imask), iseed = as.integer(iseed),
        nsamps = as.integer(nMCMC), ithin = as.integer(mcmc$ithin), pop = as.double(pop), mij = as.double(mij), rij = as.double(rij),
        Tg = as.double(Tg), e_nonflu = as.double(e_nonflu), nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]), rtn = as.double(rtn),
        accept = double(0), tab = as.single(mcmc$tab), iBest = as.integer(1), nperiodsFit = as.integer(nperiodsFit))

    rtn = matrix(solution$rtn, nr = nperiods, nc = n)

    tab = array(solution$tab, c(nlines, (nparam + 1), nwalk))

    iBest = solution$iBest

    # write the MCMC statistics

    success = emcee.write(region = region, tab = tab, regionDF = regionDF, model = model, mcmc = mcmc, nwalk = nwalk, iBest = iBest,
        accept = solution$accept)

    # calculate the error in the CAR for each region - this uses the best result we have

    carErr = calcCARErr(rtn = rtn, mydata = mydata)

    # call the routine that will plot the results This routine will also calculate randomly chosen profiles
    nRnd = 100
    profile = array(data = 0, dim = c(nRnd, nperiods, n))
    iburn = 1

    if (restart == FALSE)
        iburn = nlines/5

    for (irnd in 1:nRnd) {
        iline = iburn + floor(runif(1, 0, 1) * (nlines - iburn))
        par = tab[iline, 1:nparam, iBest]

        solution = .Fortran("genprofiles", n = as.integer(n), epi = as.double(regionDF$epi), gamaepi = as.double(regionDF$gamaepi), sh = as.double(sh),
            school = as.double(school), nparam = as.integer(length(par)), iparam = as.integer(iparam), par = as.double(par), pop = as.double(pop),
            mij = as.double(mij), rij = as.double(rij), Tg = as.double(Tg),
            e_nonflu = as.double(e_nonflu), nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
            rtn = as.double(matrix(data = 0, nr = nperiods, nc = n)), nperiodsFit = as.integer(nperiodsFit))

        profile[irnd, , ] = matrix(solution$rtn, nr = nperiods, nc = n)
    }
    success = plotFitCDC(rtn = rtn, profile = profile, regionDF = regionDF, ireal = iBest, carErr = carErr[52, ], device = device)

    # save a binary RData file

    input = list(region = region, regionDF = regionDF, e_nonflu = e_nonflu, iseed = iseed, w = w[1], Tg = Tg[1], one_R0 = one_R0, cadence = cadence,
        dt = dt, nMCMC = nMCMC, nlines = nlines, nwalk = nwalk)

    # use the last line of tab for the value of the restart parameters

    success = saveRData.emcee(input = input, nwalk = nwalk, solution = solution, par = tab[nlines, 1:nparam, 1:nwalk], parmax = parmax,
        parmin = parmin, vecopt = vecopt, iparam = iparam)

    success = printParam(regionDF = regionDF, input = input, rtn = rtn, par = solution$par, tab = tab)

    out <- list(rtn = rtn, par = solution$par, parmin = parmin, parmax = parmax, pois = solution$pois, tab = tab)

    return(out)

}


emcee.write <- function(region = region, tab = NULL, regionDF = NULL, model = model, mcmc = mcmc, nwalk = 10, iBest = 1, accept = NULL) {


    # Write an RData file of the chains for an EMCEE  \code{DICE} run
    #
    # @section Warning:
    # \bold{This function is no longer supported}


    saveTab = tab
    tab = saveTab[, , iBest]

    myName = paste(region, "-", regionDF$FY, sep = "")

    nperiods = regionDF$nperiods
    weeks = regionDF$weeks

    # convert MCMC output back to matrix and to an MCMC object
    nparam <- length(model$veopt)

    colnames(tab) <- c(model$vecopt, "AICc")
    # convert the negLLK to actual AICc score

    nopt = length(model$vecopt)

    # For Statistics use the best chain
    tab[, "AICc"] <- 2 * tab[, "AICc"] + 2 * nopt
    tab[, "AICc"] <- tab[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

    imax = dim(tab)[1]
    # how many steps to burn - here we set it to half which is very rigid
    iburn <- imax/5

    # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics
    results <- mcmc(data = tab, start = (mcmc$ithin * iburn), end = mcmc$nMCMC, thin = mcmc$ithin)

    # check to see if 'data' sub-directory exists, if not create it
    subDir <- getwd()
    subDir = paste(subDir, "/output", sep = "")

    if (!file.exists(subDir)) {
        dir.create(file.path(subDir))
        cat(" Created ", subDir, "Directory for all the Data of MCMC chain \n")
        text <- paste(" Created ", subDir, " Directory for all the Data of MCMC chain", "\n", sep = "")
    }
    # here write an R data file print the chains statistics to the screen-only for optimized variables
    print(summary(results, quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 4))


    # create csv file with mcmc statistics
    mcmc.summary <- summary(results, quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE))
    # this is a bug in R - when a parameter is fixed sometimes NA comes in the Time-series SE

    mcmc.quantiles <- mcmc.summary$quantiles
    mcmc.stats <- mcmc.summary$statistics

    # retrive the mean AICc of this chain
    mean.llk <- mcmc.stats["AICc", "Mean"]

    # randomly select 100 sets of parameters from the chain
    myparam <- matrix(nr = 100, nc = nopt)
    colnames(myparam) <- model$vecopt
    for (i in 1:100) {
        irnd <- runif(1, min = iburn, max = dim(tab)[1])
        irnd <- round(irnd)
        myparam[i, ] <- tab[irnd, 1:nopt]
    }


    # give the parameter names to the rows
    rownames(mcmc.quantiles) <- c(model$vecopt, "AICc")
    rownames(mcmc.stats) <- c(model$vecopt, "AICc")

    # create two csv files with the stats and quantiles
    filename <- paste(subDir, "/param-stats-", myName, "-", iBest, ".csv", sep = "")
    write.csv(mcmc.stats, file = filename)
    cat("\n Writing MCMC Statistics for this Chain to File: ", filename, "\n")

    filename <- paste(subDir, "/param-quantiles-", myName, "-", iBest, ".csv", sep = "")
    write.csv(mcmc.quantiles, file = filename)
    cat("\n Writing MCMC Quantiles for this Chain to File: ", filename, "\n")

    # Now write a more condensed Table with only the optimized parameters and Mean, SD and quantiles
    mcmc.both <- cbind(mcmc.stats, mcmc.quantiles)
    colnames(mcmc.both) <- c(colnames(mcmc.stats), colnames(mcmc.quantiles))
    rownames(mcmc.both) <- rownames(mcmc.stats)

    mcmc.both <- mcmc.both[c(model$vecopt, "AICc"), ]

    mcmc.both <- mcmc.both[, c("Mean", "SD", colnames(mcmc.quantiles))]
    filename <- paste(subDir, "/param-table-", myName, "-", iBest, ".csv", sep = "")
    write.csv(mcmc.both, file = filename)
    cat("\n Writing MCMC Condensed Statistics for this Chain to File: ", filename, "\n")

    cat("\n  Acceptance rate for Chain: ", accept, "\n")

    filename <- paste(subDir, "/emcee-", myName, ".RData", sep = "")
    cat("\n Writing R object Data file for this EMCEE Run: ", filename, "\n")

    # save All the chains for all the walkers here
    emcee.list = list()
    for (iwalk in 1:nwalk) {
        emcee.list[[iwalk]] = mcmc(data = saveTab[, , iwalk], start = (mcmc$ithin), end = mcmc$nMCMC, thin = mcmc$ithin)
    }
    results <- mcmc.list(emcee.list)
    # results <- mcmc.list(saveTab,start=(mcmc$ithin),end=mcmc$nMCMC,thin=mcmc$ithin)
    save(results, model, file = filename)
    save.image()

    list(mcmc = mcmc, mean.llk = mean.llk, chain.param = myparam)

}

saveRData.emcee <- function(nwalk = 10, input = NULL, solution = NULL, par = NULL, parmax = NULL, parmin = NULL, vecopt = NULL, iparam = NULL) {

    # Write an RData file with all the input parameters of an EMCEE \code{DICE} run
    #
    # @section Warning:
    # \bold{This function is no longer supported}


    FY = input$regionDF$FY
    solution$par = par
    solution$parmin = parmin
    solution$parmax = parmax
    solution$vecopt = vecopt
    solution$iparam = iparam

    mainDir = getwd()
    outDir = "output"
    subDir = paste(mainDir, "/", outDir, sep = "")
    if (!file.exists(subDir))
        dir.create(file.path(subDir))
    fileRestart = paste("restart-emcee-", FY, ".RData", sep = "")
    filePath = paste(subDir, "/", fileRestart, sep = "")
    if (file.exists(filePath)) {
        filePathOld = paste(filePath, ".old", sep = "")
        file.rename(from = filePath, to = filePathOld)
        cat("\n Renaming restart file from: ", filePath, " to ", filePathOld, "\n")
        cat("\n Saving Fitting Information to: ", filePath, "\n")
        save(input, solution, file = filePath)

    } else {
        cat("\n Saving Fitting Information to: ", filePath, "\n")
        save(input, solution, file = filePath)
    }


    err = 0

    return(err)

}


calc.cdc.null <- function(mydata = NULL, all_years_epi = NULL, mod_id = NULL, fit_id = NULL) {
    # Calculate the Historic NULL Model
    #
    # \code{calc.null} Creates a NULL model for CDC data using the average weekly mydata
    # @param mydata A complex list with all available mydata for a given disease season model/fit spatial levels and data type
    # for all the regions in the country/level selection
    # @param all_years_epi A dataframe with all available incidence data for all years
    # @param mod_id the abbreviated name for the model region
    # @param fit_id the abbreviated name for the fit regions
    # @examples
    # calc.cdc.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
    # @return A list with a NULL model for the model and fit regions
    #

    cadence = mydata$cadence
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    if (cadence == 'Monthly') {
    	my_row = which(all_years_epi$years == mydata$years[nperiodsFit] & all_years_epi$months == mydata$months[nperiodsFit])
    	cadence.vec = mydata$months
    }

    if (cadence == 'Weekly') {

    	my_row = which(all_years_epi$years == mydata$years[nperiodsFit] & all_years_epi$weeks == mydata$weeks[nperiodsFit])
    	cadence.vec = mydata$weeks

    }

    year.vec = mydata$years

    nregions = mydata$fit$nregions

    attr = c("year", cadence)

    cases.ave.mod  = array(0, c(nperiods, (length(attr)+1)))

    colnames(cases.ave.mod) = c(attr, mod_id)
    cases.ave.mod[, cadence] = cadence.vec
	cases.ave.mod[, 'year' ] = year.vec
 	cases.ave.fit  = array(0, c(nperiods, (length(attr)+nregions)))

    colnames(cases.ave.fit)  = c(attr, fit_id)  # just allocates the right number of month for the mydata frame of the null model
    cases.ave.fit[, cadence] = cadence.vec
    cases.ave.fit[, 'year' ] = year.vec

    istart = my_row - nperiods * 10
    istart = max(istart, 1)
    past.cases.mod = all_years_epi$model$raw[istart:(my_row)]
    past.cases.fit = all_years_epi$fit$raw[istart:(my_row), ]

    for (j in 1:nperiods) {
    	if (mydata$cadence == "Monthly")
    		ind = which(all_years_epi$months[istart:my_row] == mydata$months[j])
    	if (mydata$cadence == "Weekly")
    		ind = which(all_years_epi$weeks[istart:my_row] == mydata$weeks[j])
    	cases.ave.mod[j, mod_id] = mean(past.cases.mod[ind])

    }

    for (i in 1:nregions) {


    	for (j in 1:nperiods) {
    		if (mydata$cadence == "Monthly")
    			ind = which(all_years_epi$months[istart:my_row] == mydata$months[j])
    		if (mydata$cadence == "Weekly")
    			ind = which(all_years_epi$weeks[istart:my_row] == mydata$weeks[j])
    		cases.ave.fit[j, fit_id[i]] = mean(past.cases.fit[ind, i])

    	}
    }


    cases.ave = list(cases.ave.mod = cases.ave.mod, cases.ave.fit = cases.ave.fit)

    return(cases.ave)

}



fitOnePatchDengue <- function(mydata = mydata, all_years_epi = NULL, run.list = NULL, opt.list = NULL, ireal = 1, iseed = NULL) {
    # Driver for running a single patch - Dengue
    #
    # Run a single calculation on a specific country/region/state for a given dengue season using a given method and model
    # @param mydata - dataframe with all the data for the run
	# @param all_years_epi the epi data for all years
    # @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    #  These values are set based on the user chosen model
    # for the basic reproduction value.
    # @param run.list A list with parameters needed for the MCMC procedure
    # @param ireal Integer - the MCMC chain number.  Default is 1.
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is
    # generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @return A list with the following arguments:
    # \describe{
    # \item{model_rtn}{The best result of the MCMC procedure for the data}
    # \item{model_profile}{Randomly selected results from the MCMC procedure of fitting the data}
    # \item{tab.model}{The MCMC history of fitting the data}
    # }
    # fitOnePatchDengue{mydata=mydata,all_years_epi = all_years_epi, run.list = run.list,
    #               	opt.list = opt.list, ireal = ireal, iseed = iseed}
    #


	if (mydata$cadence == "Weekly") {
		cadence = "week"
	} else if (mydata$cadence == "Monthly") {
		cadence = "month"
	} else if (mydata$cadence == "Daily") {
		cadence = "day"
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
    # Here we set the random number generation seed for R
    set.seed(iseed)

 	epi_model = mydata$epi_model


    opt.list = opt.list[[1]]

	nreal  = run.list$nreal
	nMCMC  = run.list$nMCMC
	nlines = run.list$nlines
	ithin = run.list$ithin
	device = run.list$device
	subDir = run.list$subDir

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

	mod_id = mydata$model$attr[[paste0("ABBV_", mod_level)]]
	fit_id = mydata$fit$attr[[  paste0("ABBV_", fit_level)]]

	mod_name = mydata$model$name
	fit_name = mydata$fit$name

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

	## First fit the mod_level

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)


	cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")
	## mbn changed until here
	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod

	tables.mod$epi.null[nfit, , mod_id] = epi.null.mod[, mod_id]


	## Fit the mod_level
	## Start by calculating correlation/distance with past years

	corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)

	distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)

	epi.da.df = get.sql.da(nfit = nfit, years = all_years_epi$years, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
		corrMAT = corrMAT, distMAT = distMAT, deng.null = epi.null.mod, my_id = mod_id)

	mydata$model$wght = epi.da.df$wght

	mydata$model$epirun = epi.da.df$epi.new

	gamaepi = rep(0, mydata$nperiods)
	for (i in 1:length(mydata$model$epirun)) {
		gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
	}
	mydata$model$gamaepirun = gamaepi

	tps = mydata$ndays

	nmydata = length(tps)

	nperiods = mydata$nperiods

	wght = mydata$model$wght

	par_names <- set.param.list(epi_model = mydata$epi_model)

    setup = setup.dengue.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    # these are the same for both model and fit

    nparam = setup$nparam
	nopt    = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec


    tab = setup$tab
    nparam = length(model_par)

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC


	imask = set.imask(par_names = par_names, opt.list = opt.list)

	# The number of paramters that are optimized
    nopt = length(which(imask == 1))

    iverbose = 1  # TRUE, 0 == FALSE
    iupdate = 1

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL

    parBest = par
    pois = 0
    nblock = 3
    n = dim(mydata$fit$epirun)[2]

    accept.vec = array(0, c(n, nblock))

	cases = mydata$model$epirun
	sh = mydata$model$sh
	school = mydata$model$school
	gamaepi = mydata$model$gamaepirun
    wght = mydata$model$wght
    mypar = model_par

	accept_rate = 0

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

	rtn = rep(0, nperiods)

	    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")
	## set the seed again - not really required
	iseed <- set.iseed()
	## Number of days - and patch it with a month or a week before and after
	ndays = sum(tps) + tps[1] + tps[nmydata]
	ymu = rep(0, nparam)
	sigma = rep(1, nparam)

	out <- .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school), wght = as.double(wght),
		nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax), step = as.double(model_dx),
		ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scale = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed),
		nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nmydata = as.integer(nmydata), tps = as.double(cumsum(tps)), rtn = as.double(rtn), accept_rate = as.double(rep(0,nblock)),
		curMin = as.double(1e+05), tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profiles = as.single(model_profile),
		nRnd = as.integer(nRnd), epi_model = as.integer(mydata$epi_model))

	cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")

	model_rtn = out$rtn
	tab.model = matrix(out$tab, ncol = (nparam + 1))
	model_profile = array(out$profiles, c(nRnd, nmydata))

	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nfit, state_id = mod_id)
	tables.mod = tables


    success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, mydata = mydata,
        imask = imask, ireal = ireal)

     # # save a binary RData file of the input of this run

    input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, input = input, ireal = ireal, subDir = run.list$subDir)

   # success = writeCSVOnePatch(mydata = mydata, run.list = run.list, tab = tab, tab.model = tab.model, model_rtn = model_rtn, model_profile = model_profile,
   #ireal = ireal)


	## Plot results

	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id,
		ymax.input = NULL, ireal = ireal)

	err <- plotMCMC(mydata = mydata, tab.model = tab.model, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal)

    list(fit_rtn = NULL, model_rtn = model_rtn, fit_profile = NULL, model_profile = model_profile, tab.model = tab.model, tab.fit = NULL)

}



fitSingleDengue <- function(mydata = mydata, all_years_epi = NULL, run.list = NULL, opt.list = NULL, ireal = 1, iseed = NULL) {
    # Driver for an Uncoupled Dengue Run
    #
    #
    # A spatailly uncoupled MCMC fit of the model data. The code first fits the model data directly and then
    #   fits each of the sub-regions sequentially - minimizing the likelihood of each one. The final indirect
    #   model fit is obtained as a weighted sum of these individual fits with the weights given by the relative
    #   population of each region. The data can be either cdc or gft data, and the model/fit data should have
    #   different spatial scales. For example in the case of cdc/gft data: the model can be national and the fit are
    #   the ten HHS regions. Or the model can be an HHS regionand the fit are the states in that region.
    # @param mydata - dataframe with all the data for the run
	# @param all_years_epi the epi data for all years
    # @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    # These values are set based on the user chosen model
    # for the basic reproduction value.
    # @param run.list A list with parameters needed for the MCMC procedure
    # @param ireal Integer - the MCMC chain number.  Default is 1.
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is
    # generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
  # @return A list with the following arguments:
    # \describe{
    # \item{model_rtn}{The best result of the MCMC procedure for the model data}
    # \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    # \item{tab.model}{The MCMC history of the direct fit to the model data}
    # \item{fit_rtn}{The best result for indirectly fitting the model data using the fit regions}
    # \item{fit_profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    # \item{tab.fit}{The MCMC history of indirectly fitting the model data using the fit regions}
    # }
    # fitOnePatchDengue{mydata=mydata,all_years_epi = all_years_epi, run.list = run.list,
    #               	opt.list = opt.list, ireal = ireal, iseed = iseed}
    #

 	epi_model = mydata$epi_model

	if (mydata$cadence == "Weekly") {
		cadence = "week"
	} else if (mydata$cadence == "Monthly") {
		cadence = "month"
	} else if (mydata$cadence == "Daily") {
		cadence = "day"
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
    # Here we set the random number generation seed for R
    set.seed(iseed)

    opt.cpl = opt.list[[2]]
    opt.list = opt.list[[1]]

	nreal  = run.list$nreal
	nMCMC  = run.list$nMCMC
	nlines = run.list$nlines
	ithin = run.list$ithin
	device = run.list$device
	subDir = run.list$subDir

	mod_id = mydata$model$attr[[paste0("ABBV_", mod_level)]]
	fit_id = mydata$fit$attr[[paste0("ABBV_", fit_level)]]

	mod_name = mydata$model$name
	fit_name = mydata$fit$name


	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)


	if (mod_level < fit_level) {

		tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
	}

	## First fit the mod_level

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)
	epi.obsrv = mydata$fit$raw

	tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

	cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod
	epi.null.fit = epi.null$cases.ave.fit


	tables.mod$epi.null[nfit, , mod_id] = epi.null.mod[, mod_id]

	if (mod_level < fit_level) {

		tables.fit$epi.null[nfit, , fit_id] = epi.null.fit[, fit_id]
	}


	## Fit the mod_level
	## Start by calculating correlation/distance with past years

	corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)

	distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)


	epi.da.df = get.sql.da(nfit = nfit, years = all_years_epi$years, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
		corrMAT = corrMAT, distMAT = distMAT, deng.null = epi.null.mod, my_id = mod_id)

	mydata$model$wght = epi.da.df$wght

	mydata$model$epirun = epi.da.df$epi.new

	gamaepi = rep(0, mydata$nperiods)
	for (i in 1:length(mydata$model$epirun)) {
		gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
	}
	mydata$model$gamaepirun = gamaepi

	## Repeat for fit_level


	if (fit_level > mod_level) {


		mydata$fit$epirun = mydata$fit$epi

		mydata$fit$gamaepirun = mydata$fit$gamaepi

		mydata$fit$wght = mydata$fit$epi

		for (iregion in 1:mydata$fit$nregion) {

			corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nfit)

			distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nfit)

			epi.da.df = get.sql.da(nfit = nfit, years = all_years_epi$years, mydata = mydata, cases = all_years_epi$fit$epi[,
				iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, distMAT = distMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

			mydata$fit$wght[, iregion] = epi.da.df$wght

			mydata$fit$epirun[, iregion] = epi.da.df$epi.new

			for (i in 1:length(mydata$fit$epirun[, iregion])) {
				gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
			}
			mydata$fit$gamaepirun[, iregion] = gamaepi

		}

	}


	tps = mydata$ndays

	nmydata = length(tps)

	nperiods = mydata$nperiods

	wght = mydata$model$wght

	par_names <- set.param.list(epi_model = mydata$epi_model)

    setup = setup.dengue.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = tps, par_names = par_names)

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    # these are the same for both model and fit

    nparam = setup$nparam
	nopt    = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec


    tab = setup$tab
    nparam = length(model_par)

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

     setup = setup.dengue.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
         tps = tps, par_names = par_names)

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par


	imask = set.imask(par_names = par_names, opt.list = opt.list)

	# The number of paramters that are optimized
    nopt = length(which(imask == 1))

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL
    tab.fit = list()

    pois = 0
    nblock = 3
    n = dim(mydata$fit$epirun)[2]

    accept.vec = array(0, c(n, nblock))

	cases = mydata$model$epirun
	sh = mydata$model$sh
	school = mydata$model$school
	gamaepi = mydata$model$gamaepirun
    wght = mydata$model$wght
    mypar = model_par

	accept_rate = 0

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

	rtn = rep(0, nperiods)

	cat("\n\nSEIR Fitting to the Data\n\n")
	## set the seed again - not really required
	iseed <- set.iseed()
	## Number of days - and patch it with a month or a week before and after
	ndays = sum(tps) + tps[1] + tps[nmydata]
	ymu = rep(0, nparam)
	sigma = rep(1, nparam)

	out <- .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school), wght = as.double(wght),
		nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax), step = as.double(model_dx),
		ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scale = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed),
		nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nmydata = as.integer(nmydata), tps = as.double(cumsum(tps)), rtn = as.double(rtn), accept_rate = as.double(rep(0,nblock)),
		curMin = as.double(1e+05), tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profiles = as.single(model_profile),
		nRnd = as.integer(nRnd), epi_model = as.integer(mydata$epi_model))


	model_rtn = out$rtn
	tab.model = matrix(out$tab, ncol = (nparam + 1))
	model_profile = array(out$profiles, c(nRnd, nmydata))

	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nfit, state_id = mod_id)
	tables.mod = tables


	rtn = array(0, c(nperiods, n))
	profile = array(0, c(nRnd, nperiods, n))

    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")

    for (iregion in 1:n) {

        cat("\n Direct Uncoupled Fitting of: ", mydata$fit$name[iregion], "\n\n")

        # Set the seed iseed = set.iseed()

        mypar = fit_par[, iregion]

        cases = mydata$fit$epirun[, iregion]
        sh = mydata$fit$sh[, iregion]
        school = mydata$fit$school[, iregion]
        gamaepi = mydata$fit$gamaepirun[, iregion]
        wght = mydata$fit$wght[, iregion]

        solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
            wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin[, iregion]), parmax = as.double(fit_pmax[,
                iregion]), dx = as.double(fit_dx[, iregion]), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma),
            scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
            ithin = as.integer(ithin),nmydata = as.integer(nmydata), tps = as.double(cumsum(tps)), rtn = as.double(rtn[, iregion]), accept = as.double(rep(0, nblock)),
		curMin = as.double(1e+05), tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays),profile = 	    as.single(profile[,,iregion]), nRnd = as.integer(nRnd), epi_model = as.integer(mydata$epi_model))

        rtn[, iregion] = solution$rtn

        tab.fit[[iregion]] = matrix(solution$tab, ncol = (nparam + 1))
        accept.vec[iregion, 1:nblock] = solution$accept


        profile[, , iregion] = array(solution$profile, c(nRnd, nperiods))

		tables = calc.null.err(tables = tables.fit, nfit = nfit, state_id = fit_id[iregion])
		tables.fit = tables

		tables = calc.mech.err(tables = tables.fit, profiles = profile[,,iregion], nfit = nfit, state_id = fit_id[iregion])
		tables.fit = tables

    }


    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed\n")

    # # write the MCMC statistics

    success = mcmc.single.write(tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, mydata = mydata,
        imask = imask, ireal = ireal)

     # # save a binary RData file of the input of this run

    input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, input = input, ireal = ireal, subDir = run.list$subDir)

    #success = writeCSV(mydata = mydata, run.list = run.list, tab = tab, tab.model = tab.model, model_rtn = model_rtn, model_profile = model_profile,
    #    rtn = rtn, profile = profile, ireal = ireal)


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
	tables = calc.mech.err(tables = tables.agg.mod, profiles = fit_model_profile, nfit = nfit, state_id = mod_id)
	tables.agg.mod = tables

	## Plot results

	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal)

	err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal)

    list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)

}

fitSingle_old <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    # Driver Routine for an Uncoupled Spatial Model
    #
    # A spatailly uncoupled MCMC fit of the model data. The code first fits the model data directly and then
    #   fits each of the sub-regions sequentially - minimizing the likelihood of each one. The final indirect
    #   model fit is obtained as a weighted sum of these individual fits with the weights given by the relative
    #   population of each region. The data can be either cdc or gft data, and the model/fit data should have
    #   different spatial scales. For example in the case of cdc/gft data: the model can be national and the fit are
    #   the ten HHS regions. Or the model can be an HHS regionand the fit are the states in that region.
    # @param mydata A complex list with all available data for a given season model/fit spatial levels and data type
    # @param all_years_epi A complex list with all available data for all seasons model/fit spatial levels and data type
    # @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.   These values are set based on the user chosen model
    # for the basic reproduction value.
    # @param run.list A list with parameters needed for the MCMC procedure
    # @param ireal Integer - the MCMC chain number.  Default is 1.
    # # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @return A list with the following arguments:
    # \describe{
    # \item{model_rtn}{The best result of the MCMC procedure for the model data}
    # \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    # \item{tab.model}{The MCMC history of the direct fit to the model data}
    # \item{fit_rtn}{The best result for indirectly fitting the model data using the fit regions}
    # \item{fit_profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    # \item{tab.fit}{The MCMC history of indirectly fitting the model data using the fit regions}
    # }
    # @examples
    # fitSingle{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}

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


    # Take only the first list from opt.list, the second is optional and is needed only for a coupled run

    opt.list = opt.list[[1]]

   nperiods = mydata$nperiods
   nperiodsFit = mydata$nperiodsFit

   prior = mydata$prior
   Temp = mydata$Temp
   n = mydata$fit$nregion


	if (mydata$disease == "flu") {

		out <- setupODEforSingleCDC(mydata = mydata)

		tps = out$tps

		rtn = out$rtn

		time = out$t0

		t0 = rep(out$t0, n)

	}


    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    mod_name = mydata$model$name
    fit_name = mydata$fit$name

    if (mydata$disease == "dengue") {

    	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

    	if (mod_level < fit_level) {

    		tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
    	}

    	setup = setup.dengue.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

    } else {

    	setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)
    }

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    ## Just for the CDC prediction limit R0

    if (mydata$data_source == "cdc") {

    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 1.4

    }


    # these are the same for both model and fit

    nparam = setup$nparam
	nopt    = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec


    tab = setup$tab
    nparam = length(model_par)

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    setup = setup.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = out$tps, par_names = par_names)

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par

    ## Just for the CDC prediction limit pC and R0

    if (mydata$data_source == "cdc") {
	fit_pmax['R0', 1:n] = 1.4
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

    if (mydata$fit_level < 4 & mydata$mod_level < 4 & mydata$imodel != 5) {

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

    cases = mydata$model$epi
    sh = mydata$model$sh
    school = mydata$model$school
    gamaepi = mydata$model$gamaepi
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

    ## Need nperiods + 2 because we patch the weeks at the beginning and end
    ndays = (nperiods + 2) * cadence

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

    profile = array(data = 0, dim = c(nRnd, nperiods, n))

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))

    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile, c(nRnd, nperiods))


    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")
    for (iregion in 1:n) {

        cat("\n Direct Uncoupled Fitting of: ", mydata$fit$name[iregion], "\n\n")

        # Set the seed iseed = set.iseed()

        mypar = fit_par[, iregion]

        cases = mydata$fit$epi[, iregion]
        sh = mydata$fit$sh[, iregion]
        school = mydata$fit$school[, iregion]
        gamaepi = mydata$fit$gamaepi[, iregion]
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
	    ## Need nperiods + 2 because we pathc the weeks at the beginning and end
    	ndays = (nperiods + 2) * cadence
        solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
            wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin[, iregion]), parmax = as.double(fit_pmax[,
                iregion]), dx = as.double(fit_dx[, iregion]), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma),
            scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
            ithin = as.integer(ithin),nperiods = as.integer(nperiods),
            tps = as.double(tps[1:nperiods]), rtn = as.double(rtn[, iregion]), accept = as.double(rep(0, nblock)), pois = as.double(pois),
            tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays),profile = 	    as.single(profile[,,iregion]), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))

        rtn[, iregion] = solution$rtn

        tab.fit[[iregion]] = matrix(solution$tab, ncol = (nparam + 1))
        accept.vec[iregion, 1:nblock] = solution$accept


        profile[, , iregion] = array(solution$profile, c(nRnd, nperiods))

    }
    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed\n")

    # # write the MCMC statistics

    success = mcmc.single.write(tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, mydata = mydata,
        imask = imask, ireal = ireal)

    # # calculate the error in the CAR for each region

    carErr = calcCARErr(rtn = rtn, mydata = mydata)

     # # save a binary RData file of the input of this run

    input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, input = input, ireal = ireal, subDir = run.list$subDir)

    success = writeCSV(mydata = mydata, run.list = run.list, tab = tab, tab.model = tab.model, model_rtn = model_rtn, model_profile = model_profile,
        rtn = rtn, profile = profile, ireal = ireal)
    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

    if (as.character(toupper(run.list$plot)) == "TRUE" || as.character(run.list$plot) == "1") {
        for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
            mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

        for (k in 1:length(run.list$device)) success = plotHists(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
            mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)
    }

    if (tolower(as.character(run.list$plot)) == "external" || as.character(run.list$plot) == "2") {
        for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI.external(rtn = rtn, profile = profile, model_rtn = model_rtn,
            model_profile = model_profile, mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

        for (k in 1:length(run.list$device)) success = plotHists.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
            mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

    }


    ## Create frames for movie if requested by the User - currently works only for CONUS and fit_level = 3, mod_level = 2

    if (run.list$movie == 1 && (mydata$model$level == 2) && (mydata$fit$level == 3)) {
        cat("\n Creating PNG frames \n")
        success = makeMovieFrames.ggplot2(rtn = rtn, profile = profile, mydata = mydata, ireal = ireal, run.list = run.list, convert = TRUE)
    }

    err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal)

    results = list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)

    return(results)
}

fitMulti_old <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    # Driver Routine for a Coupled Spatial Model
    #
    # A spatailly coupled MCMC fit of the model data. The code first fits the model data directly and then fits it
    #   as a weighted sum of the coupled fit level data. This fit uses a coupling matrix to describe the interaction between
    #   different spatial regions and it generates all the fit level profiles at once and minimizes their weighted
    #   likelihood with the weights given by the relative population of each region. The data can be either cdc or gft data,
    #   and the model/fit data should have different spatial scales. For example in the case of cdc data: the model is
    #   national and the fit are the ten HHS regions. Or the model can be an HHS region and the fit fit data is state level data.
    # @param mydata A complex list with all available data for a given disease season model/fit spatial levels and data type
    # @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    # These values are based on the user's chosen compartmental model (SIR/SEIR) and the force of infection.
    # @param run.list A list with parameters needed for the MCMC procedure
    # @param ireal Integer - the MCMC chain number.  Default is 1.
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly.
    # Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @return A list with the following arguments:
    # \describe{
    # \item{model_rtn}{The best result of the MCMC procedure for the model data}
    # \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    # \item{tab.model}{The MCMC history of the direct fit to the model data}
    # \item{rtn}{The best result for indirectly fitting the model data using the fit regions}
    # \item{profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    # \item{tab}{The MCMC history of indirectly fitting the model data using the fit regions}
    # }
    # @examples
    # fitMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}


    # Take both the first and the second elements for opt.list.  The first is always present The second is needed only in case of a
    # coupled run
    opt.cpl = opt.list[[2]]
    opt.list = opt.list[[1]]

    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    prior = mydata$prior
    Temp = mydata$Temp

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

    epi_model = mydata$epi_model
    if (is.null(epi_model)) {
    	epi_model = 1
    }

    if (tolower(epi_model) == "sir") {
    	epi_model = 1
    }
    if (tolower(epi_model) == "seir") {
    	epi_model = 2
    }

    if (epi_model != 1 & epi_model != 2) {
    	epi_model = 1
    }

	mydata$epi_model = epi_model

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed() * ireal%%.Machine$integer.max
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generation seed for R
    set.seed(iseed)

    out <- setupODEforSingleCDC(mydata = mydata)

    tps = out$tps

    n = mydata$fit$nregion

    rtn = out$rtn

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    par_names_cpl <- set.param.cpl.list()

    setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = out$tps, par_names = par_names)

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    ## Just for the CDC prediction limit pC and R0
    model_pmax[2] = 0.2
    model_pmax[3] = 2

    # these are the same for both model and fit

    nparam = setup$nparam
    nopt = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec
    nparam = length(model_par)
    tab.model = setup$tab

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    setup = setup.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = out$tps, par_names = par_names)

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par

    ## Just for the CDC prediction limit pC and R0
    fit_pmax[2, 1:n] = 0.2
    fit_pmax[3, 1:n] = 2


    setup.prior = setup.single.prior(mydata = mydata, logvec = logvec)

    n1 = length(setup.prior$ymu[, 1])

    fit_ymu = setup.prior$ymu[1:n, ]
    fit_sigma = setup.prior$sigma[1:n, ]
    model_ymu = setup.prior$ymu[n1, ]
    model_sigma = setup.prior$ymu[n1, ]
    logvec = setup.prior$logvec

    nmodel.prior = length(model_ymu)
    nfit.prior = dim(fit_ymu)[2]

    imask = set.imask(par_names = par_names, opt.list = opt.list)

    # set up initial guess / min/max values etc for the coupling elements - there are three of them

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

    cases = mydata$model$epi
    sh = mydata$model$sh
    school = mydata$model$school
    gamaepi = mydata$model$gamaepi
    mypar = model_par

    ## Select the prior for this model region
    ymu = model_ymu
    sigma = model_sigma
    ## pad with zeros
    npad = nparam - length(ymu)
    ymu = c(ymu, rep(0, npad))
    sigma = c(sigma, rep(1, npad))

    wght = mydata$model$wght

	nRnd = 1000
	model_profile = array(0,c(nRnd, nperiods))
	profile       = array(0,c(nRnd, nperiods, n))

    ## Need nperiods + 2 because we patch the weeks at the beginning and end
    ndays = (nperiods + 2) * cadence

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")
    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), Temp = as.double(Temp),
        imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab.model),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd),epi_model=as.integer(epi_model))


    cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")

    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile,c(nRnd, nperiods))
    # nRnd = 1000

    # if (!is.null(tab.model)) {
        # model_profile <- calcRandomProfiles(sh = mydata$model$sh, school = mydata$model$school, pop = mydata$model$pop, tab = tab.model,
            # tps = tps, nRnd = nRnd, cadence = mydata$cadence)
    # }



    ## Now fit the coupled regions

    cases = as.matrix(mydata$fit$epi)
    sh = as.matrix(mydata$fit$sh)
    school = as.matrix(mydata$fit$school)
    gamaepi = as.matrix(mydata$fit$gamaepi)
    mypar = fit_par
    wght = mydata$fit$wght
    tab = array(0, c(nlines, (nparam + 1), n))
    tab.cpl = array(0, c(nlines, 2))

    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")

    solution = .Fortran("epimulti", n = as.integer(n), epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin), parmax = as.double(fit_pmax),
        dx = as.double(fit_dx), ilog = as.integer(logvec), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
        ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rtn), pois = as.double(rep(0, n)), coef = as.double(mydata$fit$coef),
        parCPL = as.double(cpl_par), pminCPL = as.double(cpl_pmin),
        pmaxCPL = as.double(cpl_pmax), stepCPL = as.double(cpl_dx), ilogCPL = as.integer(cpl_logvec), imaskCPL = as.integer(cpl_imask),
        Rij = as.double(Rij), tab = as.single(tab), tabCPL = as.single(tab.cpl), profiles = as.single(profile), nRnd = as.integer(nRnd), ndays = as.integer(ndays), epi_model=as.integer(epi_model), scales = as.double(mydata$Temp))


    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed \n")
    rtn = matrix(solution$rtn, ncol = n)

    tab.cpl = matrix(solution$tabCPL, ncol = 2)

    tab = array(solution$tab, c(nlines, (nparam + 1), n))

    profile = array(solution$profiles, c(nRnd, nperiods, n))

    success = mcmc.multi.write(tab.model = tab.model, tab = tab, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list,
        mydata = mydata, imask = imask, cpl_imask = cpl_imask, ireal = ireal)

    # # call the routine that will plot the results # This routine will also calculate randomly chosen profiles

    ##profile <- calcCPLRandomProfiles(mydata = mydata, tab = tab, tab.cpl = tab.cpl, Rij = Rij, tps = tps, nRnd = nRnd)

    # # save a binary RData file of the input of this run

    input = list(run.list = run.list, opt.list = opt.list, opt.cpl = opt.cpl)

    success = saveRData(mydata = mydata, input = input, ireal = ireal, subDir = run.list$subDir)

    success = writeCSV(mydata = mydata, run.list = run.list, tab = tab, tab.model = tab.model, model_rtn = model_rtn, model_profile = model_profile,
        rtn = rtn, profile = profile, ireal = ireal)

    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

    if (as.character(run.list$plot) == "TRUE" || as.character(run.list$plot) == "1") {
        for (k in 1:length(run.list$device)) success = plotFitCDCPercentILI(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
            mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)

        for (k in 1:length(run.list$device)) success = plotHists(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
            mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)
    }

    if (tolower(as.character(run.list$plot)) == "external" || as.character(run.list$plot) == "2") {
        for (k in 1:length(run.list$device)) {
            success = plotFitCDCPercentILI.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
                mydata = mydata, ireal = ireal, run.list = run.list, idevice = k)
        }
        for (k in 1:length(run.list$device)) {
            success = plotHists.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
                ireal = ireal, run.list = run.list, idevice = k)
        }


    }

    err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit= tab, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, ireal = ireal)

    results = list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model, rtn = rtn, profile = profile, tab = tab)

    return(results)

}



fitMultiDengue_old <- function(mydata = mydata, all_years_epi = NULL, run.list = NULL, opt.list = NULL, ireal = 1, iseed = NULL) {

    # Driver Routine for a Coupled Spatial Model - Dengue
    #
    # A spatailly coupled MCMC fit of the model data. The code first fits the model data directly and then fits it
    #   as a weighted sum of the coupled fit level data. This fit uses a coupling matrix to describe the interaction between
    #   different spatial regions and it generates all the fit level profiles at once and minimizes their weighted
    #   likelihood with the weights given by the relative population of each region. The data can be either cdc or gft data,
    #   and the model/fit data should have different spatial scales. For example in the case of cdc data: the model is
    #   national and the fit are the ten HHS regions. Or the model can be an HHS region and the fit fit data is state level data.
    # @param mydata A complex list with all available data for a given disease season model/fit spatial levels and data type
    # @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
    # These values are based on the user's chosen compartmental model (SIR/SEIR) and the force of infection.
    # @param run.list A list with parameters needed for the MCMC procedure
    # @param ireal Integer - the MCMC chain number.  Default is 1.
    # @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly.
    # Setting the seed to a known integer allows an MCMC chain to be reproducible.
    # @return A list with the following arguments:
    # \describe{
    # \item{model_rtn}{The best result of the MCMC procedure for the model data}
    # \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
    # \item{tab.model}{The MCMC history of the direct fit to the model data}
    # \item{rtn}{The best result for indirectly fitting the model data using the fit regions}
    # \item{profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    # \item{tab}{The MCMC history of indirectly fitting the model data using the fit regions}
    # }
    # @examples
    # fitMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed =iseed}


	if (mydata$cadence == "Weekly") {
		cadence = "week"
	} else if (mydata$cadence == "Monthly") {
		cadence = "month"
	} else if (mydata$cadence == "Daily") {
		cadence = "day"
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
    # Here we set the random number generation seed for R
    set.seed(iseed)


 	epi_model = mydata$epi_model

    opt.cpl = opt.list[[2]]
    opt.list = opt.list[[1]]

	nreal  = run.list$nreal
	nMCMC  = run.list$nMCMC
	nlines = run.list$nlines
	ithin = run.list$ithin
	device = run.list$device
	subDir = run.list$subDir

	nfit      = mydata$nperiodsFit
	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

	mod_id = mydata$model$attr[[paste0("ABBV_", mod_level)]]
	fit_id = mydata$fit$attr[[paste0("ABBV_", fit_level)]]

	mod_name = mydata$model$name
	fit_name = mydata$fit$name

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)


	if (mod_level < fit_level) {

		tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
	}


	## First fit the mod_level

	epi.obsrv = mydata$model$epi
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)
	epi.obsrv = mydata$fit$epi

	tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

	cat("\n Calculating NULL Model (Historic Monthly or Weekly Average)\n\n")
	## mbn changed until here
	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod
	epi.null.fit = epi.null$cases.ave.fit


	tables.mod$epi.null[nfit, , mod_id] = epi.null.mod[, mod_id]

	if (mod_level < fit_level) {

		tables.fit$epi.null[nfit, , fit_id] = epi.null.fit[, fit_id]
	}


	## Fit the mod_level
	## Start by calculating correlation/distance with past years

	corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)

	distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)


	epi.da.df = get.sql.da(nfit = nfit, years = all_years_epi$years, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
		corrMAT = corrMAT, distMAT = distMAT, deng.null = epi.null.mod, my_id = mod_id)

	mydata$model$wght = epi.da.df$wght

	mydata$model$epirun = epi.da.df$epi.new

	gamaepi = rep(0, mydata$nperiods)
	for (i in 1:length(mydata$model$epirun)) {
		gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
	}
	mydata$model$gamaepirun = gamaepi

	## Repeat for fit_level


	if (fit_level > mod_level) {


		mydata$fit$epirun = mydata$fit$epi

		mydata$fit$gamaepirun = mydata$fit$gamaepi

		mydata$fit$wght = mydata$fit$epi

		for (iregion in 1:mydata$fit$nregion) {

			corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nfit)

			distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion], nfit = nfit)

			epi.da.df = get.sql.da(nfit = nfit, years = all_years_epi$years, mydata = mydata, cases = all_years_epi$fit$epi[,
				iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, distMAT = distMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

			mydata$fit$wght[, iregion] = epi.da.df$wght

			mydata$fit$epirun[, iregion] = epi.da.df$epi.new

			for (i in 1:length(mydata$fit$epirun[, iregion])) {
				gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
			}
			mydata$fit$gamaepirun[, iregion] = gamaepi

		}

	}

	tps = mydata$ndays

	nmydata = length(tps)

	nperiods = mydata$nperiods

	wght = mydata$model$wght

	par_names <- set.param.list(epi_model = mydata$epi_model)

	par_names_cpl <- set.param.cpl.list()

    setup = setup.dengue.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = tps, par_names = par_names)

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx
    model_par = setup$par

    # these are the same for both model and fit

    nparam = setup$nparam
	nopt    = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec


    tab = setup$tab
    nparam = length(model_par)

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

     setup = setup.dengue.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
         tps = tps, par_names = par_names)

    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par


	imask = set.imask(par_names = par_names, opt.list = opt.list)

    # set up initial guess / min/max values etc for the coupling elements - there are three of them

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

	# The number of paramters that are optimized
    nopt = length(which(imask == 1))

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL
    tab.fit = list()

    pois = 0
    nblock = 3
    n = dim(mydata$fit$epirun)[2]

    accept.vec = array(0, c(n, nblock))

	cases = mydata$model$epirun
	sh = mydata$model$sh
	school = mydata$model$school
	gamaepi = mydata$model$gamaepirun
    wght = mydata$model$wght
    mypar = model_par

	accept_rate = 0

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

	rtn = rep(0, nperiods)

	cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")
	## set the seed again - not really required
	iseed <- set.iseed()
	## Number of days - and patch it with a month or a week before and after
	ndays = sum(tps) + tps[1] + tps[nmydata]
	ymu = rep(0, nparam)
	sigma = rep(1, nparam)


	out <- .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school), wght = as.double(wght),
		nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax), step = as.double(model_dx),
		ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scale = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed),
		nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nmydata = as.integer(nmydata), tps = as.double(cumsum(tps)), rtn = as.double(rtn), accept_rate = as.double(rep(0,nblock)),
		curMin = as.double(1e+05), tab = as.single(tab), imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profiles = as.single(model_profile),
		nRnd = as.integer(nRnd), epi_model = as.integer(mydata$epi_model))

	cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")


	model_rtn = out$rtn
	tab.model = matrix(out$tab, ncol = (nparam + 1))
	model_profile = array(out$profiles, c(nRnd, nmydata))

	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	rtn = array(0, c(nperiods, n))
	profile = array(0, c(nRnd, nperiods, n))

    ## Now fit the coupled regions


    cases = as.matrix(mydata$fit$epirun)

    sh = as.matrix(mydata$fit$sh)

    school = as.matrix(mydata$fit$school)

    gamaepi = as.matrix(mydata$fit$gamaepirun)

    mypar = fit_par
    wght =  as.matrix(mydata$fit$wght)

    tab = array(0, c(nlines, (nparam + 1), n))
    tab.cpl = array(0, c(nlines, 2))


    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")


    solution = .Fortran("epimulti", n = as.integer(n) , epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin), parmax = as.double(fit_pmax),
        dx = as.double(fit_dx), ilog = as.integer(logvec), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
        ithin = as.integer(ithin), nmydata = as.integer(nmydata), tps = as.double(cumsum(tps)),
        rtn = as.double(rtn), pois = as.double(rep(0, n)), coef = as.double(mydata$fit$coef),
        parCPL = as.double(cpl_par), pminCPL = as.double(cpl_pmin),
        pmaxCPL = as.double(cpl_pmax), stepCPL = as.double(cpl_dx), ilogCPL = as.integer(cpl_logvec), imaskCPL = as.integer(cpl_imask),
        Rij = as.double(Rij), tab = as.single(tab), tabCPL = as.single(tab.cpl), profiles = as.single(profile), nRnd = as.integer(nRnd), ndays = as.integer(ndays), epi_model=as.integer(mydata$epi_model), scales = as.double(mydata$Temp))


    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fit$name, " Data Completed \n")
    rtn = matrix(solution$rtn, ncol = n)

    tab.cpl = matrix(solution$tabCPL, ncol = 2)

    tab = array(solution$tab, c(nlines, (nparam + 1), n))

    profile = array(solution$profiles, c(nRnd, nmydata, n))

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
	tables = calc.mech.err(tables = tables.agg.mod, profiles = fit_model_profile, nfit = nfit, state_id = mod_id)
	tables.agg.mod = tables

	 # # write the MCMC statistics

    success = mcmc.multi.write(tab.model = tab.model, tab = tab, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list,
        mydata = mydata, imask = imask, cpl_imask = cpl_imask, ireal = ireal)


   # # call the routine that will plot the results # This routine will also calculate randomly chosen profiles

    #profile <- calcCPLRandomProfiles(mydata = mydata, tab = tab, tab.cpl = tab.cpl, Rij = Rij, tps = cumsum(tps), nRnd = nRnd)

 	tables = calc.null.err(tables = tables.fit, nfit = nfit, state_id = fit_id)
	tables.fit = tables

 	## Update the results table

 	for (iregion in 1:n) {
 		tables = calc.mech.err(tables = tables.fit, profiles = profile[,,iregion], nfit = nfit, state_id = fit_id[iregion])
		tables.fit = tables
 	}

     # # save a binary RData file of the input of this run

    input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, input = input, ireal = ireal, subDir = run.list$subDir)

    #success = writeCSV(mydata = mydata, run.list = run.list, tab = tab, tab.model = tab.model, model_rtn = model_rtn, model_profile = model_profile,
    #    rtn = rtn, profile = profile, ireal = ireal)

	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal)

	err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit= tab, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, ireal = ireal)

    list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)

}

makeMovieFrames.ggplot2 <- function(rtn = NULL, profile = NULL, mydata = NULL, ireal = 1, run.list = NULL, convert = TRUE) {

    # For a CDC or GFT run of USA create frames of weekly incidence
    #
    # For a fit_level=3 and mod_level =3 run of the US create weekly maps of ILI incidence using ggplot2.  The function creates a sub-directory where all the frames are written. Default name is frames-ireal where ireal is the MCMC chain number.  A frame is written for each week.  The user can then string these frmaes to create a movie.  Each frames shows a map of CONUS colord by the ten HHS regions and overlayed on it are circles whos radius is proportional to the \%ILI.  At each HHS region at a that week.  The circles are drawn at the (population density weighted) centroid of each HHS region.  We also show a circle for the national results - drawn at the (population density weighted) centroid of the continental US (CONUS).  Each circle show the best result along with the 75\% CI. The thin black circle denotes the data.
    #
    # @param rtn A 1D numeric array with the best in-direct prediction to the model region
    # @param mydata A list with the entire data set of this \code{DICE} run
    # @param ireal Integer - the MCMC chain number
    # @param run.list  A list with various run parameters
    # @param convert Logical - tells the code if we need to convert the profiles to \%ILI (TRUE)  or not (FALSE)
    #
    # @return  Returns \eqn{err = 0} if successful
    #
    # @examples
    # makeMovieFrames.ggplot2{rtn = rtn, profile = profile, mydata = mydata, ireal = ireal}


    FY = mydata$FY
    model = mydata$imodel
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$nperiodsData
    reg.fit.name = mydata$fit$name
    reg.model.name = mydata$model$name

    nRnd = dim(profile)[1]

    n.model = 1
    n.fit = mydata$fit$nregions
    nRnd = dim(profile)[1]

    fit_coef = mydata$fit$coef
    fit_onset = mydata$fit$onset
    model_coef = mydata$model$coef
    model_onset = mydata$model$onset
    weeks = mydata$weeks

    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName


    # Make a movie

    ## Get map-data from GADM and simplify it port1 = suppressMessages(getData('GADM', country = 'USA', level = 1))
    port1$NAME_1 = as.factor(as.character(port1$NAME_1))
    name = port1$NAME_1
    port1 = gSimplify(port1, tol = 0.01, topologyPreserve = TRUE)

    ## fills the map based on different regions
    region_list = list(Region1 = c("Maine", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "Vermont"), Region2 = c("New York",
        "New Jersey", "Puerto Rico"), Region3 = c("Pennsylvania", "Delaware", "Maryland", "West Virginia", "Virginia", "District of Columbia"),
        Region4 = c("Kentucky", "Tennessee", "North Carolina", "South Carolina", "Georgia", "Florida", "Alabama", "Mississippi"), Region5 = c("Minnesota",
            "Wisconsin", "Illinois", "Indiana", "Michigan", "Ohio"), Region6 = c("New Mexico", "Texas", "Oklahoma", "Arkansas", "Louisiana"),
        Region7 = c("Nebraska", "Kansas", "Iowa", "Missouri"), Region8 = c("Utah", "Colorado", "Wyoming", "Montana", "South Dakota",
            "North Dakota"), Region9 = c("California", "Nevada", "Arizona", "Hawaii"), Region10 = c("Oregon", "Washington", "Idaho",
            "Alaska"))
    col = rainbow(n.fit)

    # labels = mydata$fit$attr$NAME_3
    labels = mydata$fit$name

    state.list = name
    nstate = length(state.list)
    for (i in nstate) {
        k = regexpr(":", state.list[i])
        if (k != -1)
            state.list[i] = substr(state.list[i], start = 1, stop = (k - 1))
    }
    state.col = NULL
    state.label = NULL
    nregion = length(region_list)
    for (j in 1:nregion) {
        myregion = region_list[[j]]
        index = which(as.character(state.list) %in% myregion)
        state.col[index] = col[j]
        state.label[index] = labels[j]
    }
    ## change to data.frame
    map1 = fortify(port1)
    map1$id = as.integer(map1$id)
    dat = data.frame(id = 1:(length(name)), state = name)
    map1.df = inner_join(map1, dat, by = "id")


    labels = mydata$fit$attr$NAME_3

    map1.df$col = state.col[map1.df$id]
    map1.df$col = factor(map1.df$col, levels = unique(state.col)[order(unique(state.label))], labels = labels)

    longitude = mydata$fit$attr$lon
    latitude = mydata$fit$attr$lat

    factor = as.numeric(mydata$fit$factor)
    model_factor = mydata$model$factor

    ## This is the fit and model %ILI data
    fit_ili = mydata$fit$raw
    model_ili = mydata$model$raw

    ## If needed Convert the fit results to % ILI When this routine is running after a run is complete the loaded profiles are already
    ## %ILI and there is no need to convert

    rtn_ili = rtn
    profile_ili = profile

    if (convert == TRUE) {
        for (i in 1:n.fit) {
            rtn_ili[, i] = rtn[, i]/factor[i]
            profile_ili[, , i] = profile[, , i]/factor[i]
        }

    }


    ## Build the indirect model results
    fit_model = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {
        fit_model[i] = sum(rtn_ili[i, 1:n.fit] * fit_coef[1:n.fit])
        for (irnd in 1:nRnd) {
            for (k in 1:n.fit) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile_ili[irnd, i, k] * fit_coef[k]
        }

    }

    region_r = rtn_ili

    ## Quantiles for regions and nation

    quan.out = apply(profile_ili, c(2, 3), quanci)
    region.lb = quan.out[1, , ]
    region.ub = quan.out[2, , ]


    nation_long = mydata$model$attr$lon
    nation_lat = mydata$model$attr$lat

    nation_r = fit_model

    nation.out = apply(fit_model_profile, 2, quanci)
    nation.lb = nation.out[1, ]
    nation.ub = nation.out[2, ]

    maxrtn = max(region_r, nation_r, na.rm = TRUE)


    ## create a directory for the frames
    tmp <- getwd()  # get the current working directory
    setwd(subDir)  # set the working directory to the output dir
    fname <- paste0("frames-", ireal)
    # create the directory name
    err <- makeDir(subDir = fname)  # create the new directory
    setwd(tmp)  # return to the application working directory
    subDir2 <- file.path(subDir, fname)

    for (i in 1:nperiods) {
        cat("making frame number: ", i, "\n")

        if (i < 10) {
            name = paste(subDir2, "/000", i, "plot.png", sep = "")
        }
        if (i >= 10) {
            name = paste(subDir2, "/00", i, "plot.png", sep = "")
        }
        p = ggplot() + geom_map(data = map1.df, map = map1.df, aes(map_id = id, x = long, y = lat, group = group, fill = col), size = 0.25) +
            coord_map() + geom_point(aes_string(x = nation_long, y = nation_lat), size = 8 * (nation.ub[i]/maxrtn), col = "mediumpurple",
            show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = nation_long, y = nation_lat), size = 8 * (nation.lb[i]/maxrtn),
            col = "khaki", show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = nation_long, y = nation_lat), size = 8 * (nation_r[i]/maxrtn),
            col = "grey", show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = nation_long, y = nation_lat), size = 8 * (model_ili[i]/maxrtn),
            shape = 1, col = "black", show.legend = FALSE, alpha = 0.6, na.rm = TRUE) + scale_x_continuous(name = "", limits = c(-130,
            -60)) + scale_y_continuous(name = "", limits = c(25, 50)) + ggtitle(paste(FY, " Seaon - Epidemic Week ", weeks[i], sep = "")) +
            theme(legend.title = element_blank(), plot.title = element_text(size = 15, family = "serif", face = "bold"))
        for (j in 1:n.fit) {
            p = p + geom_point(aes_string(x = longitude[j], y = latitude[j]), size = 8 * (region.ub[i, j]/maxrtn), col = "forestgreen",
                show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = longitude[j], y = latitude[j]), size = 8 * (region.lb[i,
                j]/maxrtn), col = "deeppink", show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = longitude[j], y = latitude[j]),
                size = 8 * (region_r[i, j]/maxrtn), col = "#24576D", show.legend = FALSE, alpha = 0.5) + geom_point(aes_string(x = longitude[j],
                y = latitude[j]), size = 8 * (fit_ili[i, j]/maxrtn), shape = 1, col = "black", show.legend = FALSE, alpha = 0.6, na.rm = TRUE)
        }
        suppressMessages(ggsave(filename = name))
    }

    cat("  png frames were saved at: ", subDir2, "\n")
    cat(" They can be used to make a movie using for example: convert -delay 20 *png my_movie.mov", "\n")


    err = 0
    return(err)


}


## Not clear why this is here but I was not sure where to put it.  This function is called but the information is not really used

calcCARErr <- function(rtn = NULL, mydata = NULL) {

    # Calculate Weekly Error in the Cumulative Attack Rate
    #
    # Calculate the weekly error in the cumulative attack rate for each spatial region using the best result from the MCMC chain and the mydata.
    # @param rtn Numeric matrix with \emph{nperiods} rows and \emph{nregions} columns with our best prediction for the weekly number of cases.
    # @param mydata - dataframe with all the data for this \code{DICE} run
    # @return   car.err the Absolute monthly/weekly/daily \% error for each of the \emph{nregions}
    # @examples
    # calcCARErr{rtn=rtn, mydata=mydata}

    epi = mydata$fit$epi
    n = dim(epi)[2]
    nperiods = mydata$nperiods
    # # calculate the cummulative attack rate-for the mydata and the estimate

    rtn.cumsum = rtn
    epi.cumsum = epi

    car.err = matrix(0, nr = nperiods, nc = n)

    for (i in 1:n) {
        rtn.cumsum[, i] = cumsum(rtn[, i])
        epi.cumsum[, i] = cumsum(epi[, i])
        car.err[, i] = (epi.cumsum[, i] - rtn.cumsum[, i])/epi.cumsum[, i] * 100
    }
    car.err = abs(car.err)
    return(car.err)
}


build.cadence.vector.sql <- function(mydata=NULL) {

    # Given a dataframe of dengue data Build a cadence array
    #
    # \code{build.cadence.vector.sql} Given the dengue data dataframe
    # create a list which contains the months/weeks names, etc.
    # @param mydata The dengue mydata set for the run
    # @examples
    # build.cadence.vector.sql(mydata = mydata)
    #
    #


    year.vec = mydata$years
    if (mydata$cadence == "Monthly") {
    	cadence.vec = mydata$months
    	cadence.names = month.abb[cadence.vec]
    } else if (mydata$cadence == "Weekly") {
    	cadence.vec = mydata$weeks
    	cadence.names = paste0("EW", cadence.vec)
    } else {
    	cat("\n\n Currently Supporting only Monthly and Weekly Cadence \n\n")
    	cat("\n\n Code Will Stop!! \n\n")
    	quit(save = "no")
    }

    cadence.year.info.for.subset = list(cadence.vec = cadence.vec, cadence.names = cadence.names, year.vec = year.vec)

    return(cadence.year.info.for.subset)

}

setupODEforCDC <- function(region = "cdc", regionDF = NULL, R0=1.2, Tg = 2.6, pC = 1, iseed = NULL, 
    nperiods = 52, cadence = 7, dt = 0.1) {
    
    # Setup Parameters and Arrays for the ODEs
    # 
    # @section Warning:
    # \bold{This function is no longer supported}
    
    country = regionDF$mydata
    mij = regionDF$mij
    
    ncountry <- length(country$pop)
    
    # make sure the populations are integers
    country$pop = round(country$pop)
    

    week0 = regionDF$weeks[1]
    
    day0 = cadence * (week0 - 1)
    
    # if iseed is NULL we need to see the RNG
    if (is.null(iseed)) 
        iseed = set.iseed()
    
    rtn = regionDF$epi
    
    # for now just turn R0 and Tg to vectors with the same value for all countries in the future we may change this
    
    Tg <- rep(Tg, ncountry)
    R0 <- rep(R0, ncountry)
    pC <- rep(pC, ncountry)
    
    
    days_per_week = 7
    tps <- seq(from = day0, to = (day0 + nperiods * days_per_week), by = cadence)
    noPts = length(tps)
    rtn <- array(0, c(noPts - 1, ncountry))
    
    list(R0 = R0, pC = pC, Tg = Tg, tps = tps, nperiods = nperiods, rtn = rtn, t0 = day0)
    
}
