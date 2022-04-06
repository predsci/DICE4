##
## All Functions related to getting the mydata from internal or external mydatabase
##

get.DICE.data <- function(data_source = "cdc", mod_level = 2, mod_name=c(NAME_2="US"), fit_names="all", fit_level = 3, RegState=NULL, year = 2015, nperiodsFit = 52, model = 4, isingle = 0, db_opts=list(DICE_db="predsci", CDC_server=TRUE), disease="flu", epi_model = 1, method ='mech', all_years_flag=T, all_cad_clim=T, raw_col=NULL) {

  #' Retrieve all available data from the DICE database.
  #'
  #' \code{get.DICE.data} retrieves all the information for the model and fit regions from the \pkg{DICE} data base.  \pkg{DICE} currently has Google Flu Trends (GFT) and Centers for Disease Control (CDC) for the United States and Dengue data for a large number of countries.  It is assumed that one might be using finer resolution data (fit_level) to create a forecast for a larger area (mod_level).
  #' @param data_source Describes the data source for the incidence data. Default is 'cdc' (for \code{disease = 'flu'}). It can be selected by source_key (integer) or source abbreviation (string). Most disease/location combinations have only one data source.  In this case, it may be easier to set data_source=NULL.  However, when multiple data sources exist, setting data_source=NULL will essentially choose from the available sources at random.  To determine a data source by graphical interface, see: \url{predsci.com/id_data/}.  Looking-up the disease and location will result in a list of data sources that can be entered into DICE.  Alternatively, all country/disease/data_source combinations are listed in the `Data Sources Table' tab at the same url.  To access the list of sources directly from an R-prompt, see the examples below.
  #' @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has data at levels 2-4 for CDC and GFT.
  #' @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  #' @param mod_name Named vector of strings specifying the model-level spatial patch.  If \code{is.null(mod_name)}, the code reverts to using \code{RegState} (see next entry).  To specify New York state, set \code{mod_name=c(NAME_2="United States", NAME_3="R1", NAME_4="New York"). Here NAME_X is either the full name or abbreviation of the level-X patch. Replacing 'United States' with 'US' or 'R1' with 'Region 1' would result in the same outcome.  Also, vector entries for mod_name should go from NAME_2,....,NAME_n where mod_level=n.}
  #' @param fit_name A character vector indicating which fit-regions to use.  If \code{fit_name='all'}, then DICE uses all child-regions of the model region with level equal to \code{fit_level}. The other mode for fit_name is to specifiy a subset of the fit regions to construct an aggregate representation of the model region.  For example if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=3}, and \code{fit_names=c("R1", "R2", "R3")}, DICE will create an Atlantic super-region to model (as opposed to using all 10 HHS regions).  Similarly, if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=4}, and \code{fit_names=c("WA", "OR", "CA")}, DICE will create and model a super-state of Pacific states.
  #' @param RegState  \strong{Single element}: determines which single region from \code{mod_level} is to be modeled. Depending on the model level, \code{RegState} should adhere to the following format: \code{mod_level = 2} - 3-letter ISO3 country code, \code{mod_level=3} - an integer describing the HHS region, \code{mod_level=4} - a 2-letter state code. \cr
  #' @param year A Number - The starting year of the flu season  (Default value is 2017). \pkg{dice} currently has data for years 2003-2015 for CDC and 2003-2014 for GFT.
  #' @param nperiodsFit A number - the number of data periods the user wants to include in the fit.  (Default is to include all available data)
  #' @param model A number - the model number (currently we support models 1-5 for flu and 4-5 for dengue. Default is model 4 )
  #' @param isingle Integer 0 couple the fit-level regions/patched 1 do NOT couple.  Default is couple
  #' @param disease String - disease name. Options for modeling are: flu, dengue, yellow$\_$fever, ebola, zika, cholera, chik, plague. To graphically explore the data see: \url{predsci.com/id$\_$data/}. A full list of diseases in the DICE database can be found from an R-prompt by following one of the examples below.
  #' @param db_opts A list of database options.  $DICE_db Determines which SQL database the data is retrieved from.  'PredSci' is the default SQL database, 'BSVE' is in development.  Additional flags are for outside sources of data (currently only the CDC Influenza-Like_Illness (ILI) is supported: $CDC_server=TRUE).
  #' @param epi_model Numeric, 1 == sir (default) 2 == seir.  Used to build a filename for output
  #' @param method String either 'mech' for compartmental mechanistic models or 'stat' for SARIMA models. Used to build a filename for output
  #' @param all_years_flag TRUE/FALSE, grab all years of incidence data in addition to the specified season.
  #' @param all_cad_clim TRUE/FALSE, grab all years of climate data that are available.
  #' @param raw_col A string specifying which data column to use for modeling.  Data column names are associated with data sources and are listed in the data_sources table under 'col_name'.  When raw_col is NULL, DICE uses the default (first) data column.
  #  #' @param cdc_NOAA_clim TRUE/FALSE, specific to data_source='cdc' and sql_opts$CDC_server=TRUE.  get.cdc.data() automatically returns specific humidity compiled from NASA data.  If cdc_NOAA_clim is TRUE, population-weighted temperature and precipitation compiled from NOAA data will be returned as well.
  #' @return mydata - a list with all available data and auxillary information for both the model and fit data sets.
  #'
  #' For both we provide the percent weighted ILI, the number of cases, the weekly averaged specific humidity, precipitation and temerature and the school vacation schedule.
  #' For dengue - most of the data is monthly and almost all the data is number of cases.  We also provide averaged specific humidity, precipitation and temperature on the same cadence as the dengue data.
  #'
  #' The auxillary information, for both data sets, includes the populations, the lon/lat values and all the names describing the region.
  #' @examples
  #' require(DICE)
  #' # Get national and regional CDC mydata
  #' get.DICE.data(data_source = 'cdc', mod_level = 2, fit_level = 3, RegState = 'usa', year = 2015, nperiodsFit = 45, mode = 5)
  #'
  #' # Get Region9 and state GFT mydata
  #' get.DICE.data(data_source = 'gft', mod_level = 3, fit_level = 4, RegState = 9    , year = 2013, nperiodsFit = 45, mode = 5)
  #'
  #' # Create a 'west coast' region from California, Oregon, Washington
  #' get.DICE.data(data_source = 'gft', mod_level = 3, fit_level = 4, RegState = c('CA','OR','WA'), year = 2013, nperiodsFit = 45, mode = 5)
  #'
  #' Dengue data
  #' get.DICE.data(mod_level = 3, fit_level = 3, year = 1010, nperiodsFit = 12, model = 4, isingle = 0,
  #'              sql_db = TRUE, disease = 'dengue', RegState = 'BR')
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


  # Check that the input mod and fit levels are appropriate
  if (fit_level < mod_level)
    stop("Error in DICE inputs: fit-level must be >= mod_level")

    # decide between mod_name and RegState
    if (!is.null(RegState)) {
      if (is.character(RegState) && nchar(RegState)<=3) {
        if (mod_level==1) {
          mod_name=c(NAME_1=RegState)
        } else if (mod_level==2) {
          mod_name=c(NAME_2=RegState)
        } else {
          stop("Error in DICE inputs: RegState only supports mod_level = 1 or 2 for sql_opts$CDC_server=FALSE.")
        }
      } else {
        stop("Error in DICE inputs: RegState only supports integer and character abbreviations (nchar(RegState)<4) values for sql_opts$CDC_server=TRUE.")
      }
      # get.DICE.data works with ISO3, get.mysql works with ISO2.  if NAME_2 is ISO3, convert to ISO2
      if (mod_level>1 && nchar(mod_name["NAME_2"])==3) {
        mod_name["NAME_2"] = countrycode(mod_name["NAME_2"], origin="iso3c", destination="iso2c")
      }
    } else {
      # attempt to check if mod_name is properly formatted/constructed
      if (!is.character(mod_name)) {
        stop("Error in DICE inputs: mod_name must be a character vector.")
      }
      if (mod_level > 1) {
        if (!all(paste0("NAME_", 2:mod_level) %in% names(mod_name))) {
          stop("Error in DICE inputs: names(mod_name) must contain ", paste0("NAME_", 2:mod_level, collapse=", "), " for mod_level=", mod_level, ". \n")
        }
      } else if (mod_level==1) {
        if (!("NAME_1" %in% names(mod_name))) {
          stop("Error in DICE inputs: names(mod_name) must contain 'NAME_1' when mod_level=1.\n ")
        }
      }
    }

    if (db_opts$CDC_server & (is.null(data_source) || (tolower(data_source)!='cdc' &  data_source!=7))) {
      stop("If db_opts$CDC_server is set to TRUE, data_source must be either 'cdc' or 7.  To use a different data_source, please set db_opts$CDC_server to FALSE.")
    }

    # retrieve mydata from MySQL mydatabase
    mysql_out = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=mod_name, fit_names=fit_names, disease=disease, sql_data_source=data_source, season=year, all_years_flag=all_years_flag, db_opts=db_opts, cad_clim=TRUE, raw_col=raw_col)
    mydata = mysql_out$mydata
    all_years_epi = mysql_out$all_years_epi
    all_years_clim = mysql_out$all_years_clim

    # convert new 'mydata' structure to old version
    mydata = mydata_New2Old(mydata)
    # if any NAs introduced to mydata$fit$epi, remove now
    mydata$fit$epi[is.na(mydata$fit$epi)] = 0

    # check if 'school' mydata streams exist
    if (is.null(mydata$model$school)) {
      cat("Warning:  No school mydata for model-level incidence.  Inserting all zeros for school schedule. \n")
      mydata$model$school = rep(0,mydata$nperiods)
    } else {
      mydata$model$school[is.na(mydata$model$school)] = 0
    }

    if (is.null(mydata$fit$school)) {
      cat("Warning:  No school mydata for fit-level incidence.  Inserting all zeros for school schedule. \n")
      mydata$fit$school = as.data.frame(matrix(data=0, nrow=mydata$nperiods, ncol=mydata$fit$nregions))
      names(mydata$fit$school) = mydata$fit$attr$identifier
    } else {
      mydata$fit$school[is.na(mydata$fit$school)] = 0
    }

    nperiodsDataUse = mydata$nperiodsData

    mydata$nperiodsDataUse = nperiodsDataUse

    if (nperiodsFit <= 0)
      nperiodsFit = mydata$nperiodsData
    if (nperiodsFit > mydata$nperiodsData)
      nperiodsFit = mydata$nperiodsData
    mydata$nperiodsFit = nperiodsFit

    if (mydata$years[length(mydata$years)] > mydata$years[1]) {
      mydata$FY = paste0(mydata$years[1], "-", mydata$years[length(mydata$years)])
    } else {
      mydata$FY = paste0(mydata$years[1])
    }

    # mydata$season = year(start.date)
    mydata$season = year
    # set data_source to string/abbreviation to be used in filenames
    data_source = mydata$data_source
  # }

  # if data_source=="cdc", convert names to old form
  if (mydata$data_source=="cdc") {
    # first mydata$model
    if (mydata$model$level==3) {
      mydata$model$name = gsub(pattern=" ", replacement="", x=mydata$model$name, fixed=TRUE)
    } else {
      mydata$model$name = gsub(pattern=" ", replacement=".", x=mydata$model$name, fixed=TRUE)
    }
    mydata$model$attr$NAME_2 = gsub(pattern=" ", replacement=".", x=mydata$model$attr$NAME_2, fixed=TRUE)
    mydata$model$attr$NAME_3 = gsub(pattern=" ", replacement="", x=mydata$model$attr$NAME_3, fixed=TRUE)
    mydata$model$attr$NAME_4 = gsub(pattern=" ", replacement=".", x=mydata$model$attr$NAME_4, fixed=TRUE)

    # now mydata$fit
    if (mydata$fit$level==3) {
      mydata$fit$name = gsub(pattern=" ", replacement="", x=mydata$fit$name, fixed=TRUE)
    } else {
      mydata$fit$name = gsub(pattern=" ", replacement=".", x=mydata$fit$name, fixed=TRUE)
    }
    mydata$fit$attr$NAME_2 = gsub(pattern=" ", replacement=".", x=mydata$fit$attr$NAME_2, fixed=TRUE)
    mydata$fit$attr$NAME_3 = gsub(pattern=" ", replacement="", x=mydata$fit$attr$NAME_3, fixed=TRUE)
    mydata$fit$attr$NAME_4 = gsub(pattern=" ", replacement=".", x=mydata$fit$attr$NAME_4, fixed=TRUE)
  }

  mydata$imodel = model

  mydata$isingle = isingle

  mydata$method = method

  ## include compartmental model name in filename
  if(epi_model == 1) {
  	epi_name = 'sir'
  } else if (epi_model == 2) {
  	epi_name = 'seir'
  } else if (epi_model == 3) {
  	epi_name = 'vsir'
  } else if (epi_model == 4) {
  	epi_name = 'vseir'
  } else if (epi_model == 5) {
  	epi_name = 'sirb'
  } else {
  	epi_name = 'unknown'
  }

  ## if input column raw_col is specified by user, include in filename
  if (!is.null(raw_col)) {
    data_source_desc = paste0(data_source, "-", raw_col)
  } else {
    data_source_desc = data_source
  }

  if (isingle==0) {
    cpl_desc = "cpl"
  } else if (isingle==1) {
    cpl_desc = "uncpl"
  }

  if (method == "mech") {
    mydata$dataName = paste0(data_source_desc, "-", epi_name, "-", mydata$model$name, "-", cpl_desc, "-", mydata$FY, "-", mydata$imodel)
  } else {
    mydata$dataName = paste0(data_source_desc, "-", mydata$model$name, "-", mydata$FY)
  }

  # if (!db_opts$CDC_server | tolower(data_source)=="cdc") {
  if (all_years_flag | all_cad_clim) {
    complete = list(mydata=mydata)
    if (all_years_flag) {
      complete$all_years_epi = all_years_epi
      complete$all_season_dates = mysql_out$all_SE_dates
    } else {
      complete$all_years_epi = NA
    }
    if (all_cad_clim) {
      complete$all_years_clim = all_years_clim
    } else {
      complete$all_years_clim = NA
    }
    return(complete)
  } else {
    return(mydata)
  }

}



get.cdc.data <- function(mod_level=2, fit_level=3, mod_name=c(NAME_2 = "USA"), start.year=2016, start.week=27, end.year=2017, end.week=26, data_source="cdc", all_years_flag=T, NOAA_clim=T, DICE_db='predsci') {

  #' Load the full \pkg{DICE} dataset and subset by date and region.
  #'
  #' This function is generally called by \code{\link{get.DICE.data}}.
  #' For the current season, access the CDC server for most up-to-date data.
  #'
  #' @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has mydata at levels 2-3 for CDC and 2-4 for GFT.
  #' @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  #' @param mod_name A named-vector of character strings that specify which region is to be \strong{model}ed.  In other words, \code{mod_name} specifies the country, region, state, etc. of the mod_level region.  \code{mod_name} should be of the form \code{mod_name = c(NAME_2='a', NAME_3='b',..., NAME_i='x'} where i=mod_level and 'a', 'b',...,'x' are the appropriate level names. NAME_i='x' also accepts abbreviations.  Choose appropriate names from \code{\link{diceData}}.  For example, mod_name=c(NAME_2='United.States',NAME_3='Region4',NAME_4='North.Carolina') and mod_level=4 specifies North Carolina. To achieve the same result, use all abbreviations mod_name=c(NAME_2='USA',NAME_3='R4',NAME_4='NC') or a mix of names and abbreviations mod_name=c(NAME_2='USA',NAME_3='Region4',NAME_4='NC')
  #' @param start.year An integer - start year of the flu season
  #' @param end.year  An integer  - end year of the flu season (default is end.year = start.year + 1)
  #' @param start.week An integer - starting CDC week of the flu season (default is 27)
  #' @param end.week   An integer - ending CDC week of the flu season (default is 26)
  #' @param data_source A string specifying the type of mydata the user wants to model, currently we support `cdc' and `gft' (default is cdc)
  #' @param all_years_flag Logical that determines if all years of incidence are returned in addition to the specified time period.
  #' @param NOAA_clim Logical that determines if population-weighted precipitation and temperature (NOAA) are returned.  NASA specific humidity is always returned.
  #' @return mydata  A list ILI, SH, School, and Census mydata for both the model- and fit-level region(s).
  #' @examples
  #' require(DICE)
  #' mydata = get.CDC.data(mod_level = 3, fit_level = 4, mod_name = c(NAME_2='USA',NAME_3='Region4'), start.year = 2014)


  # This function loads the full DICE dataset and then subsets by date and region
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
  # if (tolower(data_source) == "cdc") {
    dataNames = colnames(diceData$CDCili)
  # } else if (tolower(data_source) == "gft") {
  #   dataNames = colnames(diceData$GFTili)
  # } else if (tolower(data_source) == "miscili") {
  #   dataNames = colnames(diceData$MiscIli)
  # } else if (tolower(data_source) == "misccases") {
  #   dataNames = colnames(diceData$MiscCases)
  # }

  if (!any(dataNames == mod_ident)) {
    stop("DICE does not have the specified mod_level mydata for mydata-type: ", data_source)
  }
  if (!all(fit_ident %in% dataNames)) {
    stop("DICE does not have the specified fit_level mydata for mydata-type: ", data_source)
  }

  # Check that non-NA ILI mydata exists for the date range

  # Build date vectors
  mydata = BuildDateVecs(cadence=2, start.date=CDCweek2date(CDCweek=start.week, year=start.year), end.date=CDCweek2date(CDCweek=end.week, year=end.year))
  mydata = pluralize(mydata)
  mydata$nperiods = length(mydata$weeks)


  # Pull attr and school mydata for fit_level
  mydata$fit = list()
  mydata$fit$level = fit_level
  mydata$fit$name = diceData$attr[fit_ind, paste0("NAME_", fit_level)]
  mydata$fit$attr = diceData$attr[fit_ind, ]
  # mydata$fit$school = diceData$school[start.ind:end.ind, fit_ident]
  # mydata$fit$school[is.na(mydata$fit$school)] = 0
  mydata$fit$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, fit_ident])
  mydata$fit$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, fit_ident]
  mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)


  # Repeat for model data
  mydata$model = list()
  mydata$model$level = mod_level
  mydata$model$name = diceData$attr[mod_ind, paste0("NAME_", mod_level)]
  mydata$model$attr = diceData$attr[mod_ind, ]
  # mydata$model$school = diceData$school[start.ind:end.ind, mod_ident]
  # mydata$model$school[is.na(mydata$model$school)] = 0
  mydata$model$pop = as.numeric(diceData$pop[diceData$pop$year == start.year, mod_ident])
  mydata$model$onset = diceData$CDCbaseline[diceData$CDCbaseline$year == start.year, mod_ident]
  mydata$model$coef = mydata$model$pop/sum(mydata$model$pop)


  # Get ili mydata
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
      } else Sys.sleep(.1)
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
      } else Sys.sleep(.1)
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
          } else Sys.sleep(.1)
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
  if (NOAA_clim) {
    mysql_dat = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=sh_mod_name, start.date=start_date, end.date=NULL, disease="flu", sql_data_source=7, db_opts=list(DICE_db=DICE_db, CDC_server=T))
  }

  # Generate all_years_epi data
  if (all_years_flag) {

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
      # season_dates = get.StartEndDates(country="US", model_year=year(min_date), disease="flu")
      myDB = OpenCon(list(DICE_db=DICE_db))
      all_season_dates = get_season_dates(myDB=myDB, continent="NA", country="US", cadence=2, disease="flu")
      season_dates = all_season_dates[all_season_dates$season==year(min_date), ]
      dbDisconnect(myDB)
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

      # copy master dates
      all_years_epi = as.list(add_ndays(data.frame(mysql_dat$mydata[1:3]), cadence=2))
      all_years_epi$model = list()
      if(!exists(mysql_dat)) {
        mysql_dat = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=sh_mod_name, start.date=start_date, end.date=NULL, disease="flu", sql_data_source=7, db_opts=list(DICE_db=DICE_db, CDC_server=T))
      }
      all_years_epi$model$raw = mysql_dat$mydata$model[[mysql_dat$mydata$model$attr$identifier]]$ILI

    }

    # calc $epi and $gammaepi
    all_years_epi$model$epi = all_years_epi$model$raw*mydata$model$factor
    all_years_epi$model$epi[is.na(all_years_epi$model$epi)] = 0
    all_years_epi$model$gamaepi = lgamma(all_years_epi$model$epi + 1)

    # retrieve SH for model region (from NASA)
    sh_df = get.sh.nasa(region=CDCreg, sub_region=sub_region, start_year=all_years_epi$year[1], start_week=all_years_epi$week[1], end_year=all_years_epi$year[length(all_years_epi$year)], end_week=all_years_epi$week[length(all_years_epi$week)])
    all_years_epi$model$sh = sh_df[, 4]

    if (NOAA_clim) {
      # retrieve temp and precip for model region (from NOAA)
      temp_clim = cbind(data.frame(mysql_dat$mydata[1:3]), mysql_dat$mydata$model[[mysql_dat$mydata$model$attr$identifier]][, c("temp", "precip")])
      temp_clim = merge(x=pluralize(all_years_epi[1:3]), y=temp_clim, all.x=T)
      all_years_epi$model$temp = temp_clim$temp
      all_years_epi$model$precip = temp_clim$precip
    }


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

      # extract fit regions
      all_years_epi$fit = list()
      if(!exists(mysql_dat)) {
        mysql_dat = get.mysql(mod_level=mod_level, fit_level=fit_level, mod_name=sh_mod_name, start.date=start_date, end.date=NULL, disease="flu", sql_data_source=7, db_opts=list(DICE_db=DICE_db, CDC_server=T))
      }
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
    sh_df = get.sh.nasa(region=CDCreg, sub_region=sub_region, start_year=all_years_epi$year[1], start_week=all_years_epi$week[1], end_year=all_years_epi$year[length(all_years_epi$year)], end_week=all_years_epi$week[length(all_years_epi$week)])
    all_years_epi$fit$sh = sh_df[, 4:ncol(sh_df)]

    if (NOAA_clim) {
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
    }

    # pluralize date vecs
    all_years_epi = pluralize(all_years_epi)
  } else {
    all_years_epi = NULL
  }

  mydata$data_source = "cdc"
  mydata$data_desc = DICE_data_sources[tolower(DICE_data_sources$source_abbv)=="cdc", ]


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

    sh = try(get.sh.nasa(region = region, sub_region = sub_region, start_year = 2000, start_week = 1, end_year = end.year, end_week = end.week),
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

    sh = try(get.sh.nasa(region = region, sub_region = sub_region, start_year = 2000, start_week = 1, end_year = end.year, end_week = end.week),
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

  if (tolower(data_source)=="cdc" & NOAA_clim) {
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
        week_ind = mysql_dat$all_years_epi$weeks==week
        mydata$model[[var]][ind] = mean(mysql_dat$all_years_epi$model[[var]][week_ind], na.rm=T)
      }
    }
    # pull temp and precip data from mysql for fit
    temp_precip = cbind(data.frame(mysql_dat$all_years_epi[1:3]), mysql_dat$all_years_epi$fit$precip)
    temp_precip = merge(x=mydata[1:3], y=temp_precip, all.x=T)
    mydata$fit$precip = temp_precip[4:length(temp_precip)]
    temp_temp = cbind(data.frame(mysql_dat$all_years_epi[1:3]), mysql_dat$all_years_epi$fit$temp)
    temp_temp = merge(x=mydata[1:3], y=temp_temp, all.x=T)
    mydata$fit$temp = temp_temp[4:length(temp_temp)]
    # fill any NAs with historic average
    for (var in c("temp", "precip")) {
      for (ii in 1:length(mydata$fit[[var]])) {
        na_ind = which(is.na(mydata$fit[[var]][, ii]))
        for (ind in na_ind) {
          week = mydata$weeks[ind]
          week_ind = mysql_dat$all_years_epi$weeks==week
          mydata$fit[[var]][ind, ii] = mean(mysql_dat$all_years_epi$model[[var]][week_ind], na.rm=T)
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


get.mysql <- function(mod_level=2, fit_level=3, mod_name=c(NAME_2 = "BR"), fit_names="all", start.date=as.Date("2010-01-01"), end.date=as.Date("2010-12-31"), disease="dengue", sql_data_source=NULL, season=NULL, all_years_flag=T, db_opts=list(DICE_db='predsci', CDC_server=F), cad_clim=TRUE, raw_col=NULL) {

  #' Retrieve incidence and climate data from MySQL database.
  #'
  #' This function is generally called by \code{\link{get.DICE.data}}. After retrieving the data is is formatted for use by \pkg{DICE}
  #' @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has mydata at levels 2-3 for CDC and 2-4 for GFT.
  #' @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  #' @param mod_name A named-vector of character strings that specify which region is to be \strong{model}ed.  In other words, \code{mod_name} specifies the country, region, state, etc. of the mod_level region.  \code{mod_name} should be of the form \code{mod_name = c(NAME_2='a', NAME_3='b',..., NAME_i='x'} where i=mod_level and 'a', 'b',...,'x' are the appropriate level names. NAME_i='x' also accepts abbreviations.  Choose appropriate names from \code{\link{diceData}}.  For example, mod_name=c(NAME_2='United.States',NAME_3='Region4',NAME_4='North.Carolina') and mod_level=4 specifies North Carolina. To achieve the same result, use all abbreviations mod_name=c(NAME_2='US',NAME_3='R4',NAME_4='NC') or a mix of names and abbreviations mod_name=c(NAME_2='US',NAME_3='Region4',NAME_4='NC'). Unlike get.cdc.data(), get.mysql() expects ISO2 country abbreviations.
  #' @param fit_name A character vector indicating which fit-regions to use.  If \code{fit_name='all'}, then DICE uses all child-regions of the model region with level equal to \code{fit_level}. The other mode for fit_name is to specifiy a subset of the fit regions to construct an aggregate representation of the model region.  For example if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=3}, and \code{fit_names=c("R1", "R2", "R3")}, DICE will create an Atlantic super-region to model (as opposed to using all 10 HHS regions).  Similarly, if \code{mod_level=c(NAME_2="US")}, \code{mod_level=2}, \code{fit_level=4}, and \code{fit_names=c("WA", "OR", "CA")}, DICE will create and model a super-state of Pacific states.
  #' @param start.date A Date-class variable - start date of fitting period.  Passing a NULL value causes the earliest available mydata to be returned.
  #' @param end.date  A Date-class variable - end date of forecasting period.  Passing a NULL value causes the latest available mydata to be returned.
  #' @param disease String - disease name. Options for modeling are: flu, dengue, yellow$\_$fever, ebola, zika, cholera, chik, plague. To graphically explore the data see: \url{predsci.com/id$\_$data/}. A full list of diseases in the DICE database can be found from an R-prompt by following one of the examples below.
  #' @param sql_data_source Describes the data source for the incidence data. Default is 'cdc' (for \code{disease = 'flu'}). It can be selected by source_key (integer) or source abbreviation (string). Most disease/location combinations have only one data source.  In this case, it may be easier to set data_source=NULL.  However, when multiple data sources exist, setting data_source=NULL will essentially choose from the available sources at random.  To determine a data source by graphical interface, see: \url{predsci.com/id_data/}.  Looking-up the disease and location will result in a list of data sources that can be entered into DICE.  Alternatively, all country/disease/data_source combinations are listed in the `Data Sources Table' tab at the same url.  To access the list of sources directly from an R-prompt, see the examples below.
  #' @param season An integer (year) specifying the season.  When only one year is needed and/or start/end dates are not certain, season can be set.  If !is.null(season), it overrides the start./end.date inputs.
  #' @param all_years_flag Logic flag indicating if all years of incidence should be returned in addition to the specified season/date-range.
  #' @param raw_col A string specifying which data column to use for modeling.  Data column names are associated with data sources and are listed in the data_sources table under 'col_name'.  When raw_col is NULL, DICE uses the default (first) data column.
  #' @return mydata  A list ILI, SH, School, and Census mydata for both the model- and fit-level region(s).
  #' @examples
  #' require(DICE)
  #' mydata = get.mysql(mod_level = 2, fit_level = 4, mod_name = c(NAME_2='BR'), start.date=as.Date("2010-01-01"), end.date=as.Date("2010-12-31"), disease="dengue", sql_data_source=1)
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

  if (is.null(season) && !is.null(start.date) && !is.Date(start.date)) {
    stop("In get.mysql(), start.date should of type Date or NULL.")
  }
  if (is.null(season) && !is.null(end.date) && !is.Date(end.date)) {
    stop("In get.mysql(), end.date should of type Date or NULL.")
  }
  if (is.null(season) && is.Date(start.date) & is.Date(end.date)) {
    if (end.date<start.date) {
      stop("In get.mysql(), end.date must be greater-than or equal to start.date.")
    }
  }

  # Initiate the output strucure 'mydata'
  mydata = list()

  # download data_sources tables
  myDB = OpenCon(db_opts)
  data_sources = dbReadTable(myDB, name="data_sources")
  clim_by_disease = dbReadTable(myDB, name="clim_by_disease")
  col_units = dbReadTable(myDB, name="col_units")
  unit_types = dbReadTable(myDB, "unit_types")

  # check if sql_data_source is numeric. If not, convert to numeric
  if (!is.null(sql_data_source)) {
    if (is.character(sql_data_source)) {
      temp_source = data_sources$source_key[tolower(sql_data_source)==tolower(data_sources$source_abbv)]
      if (length(temp_source)==0) {
        cat("Warning: Character string entered into data_source does not match any DICE source abbreviations.  DICE will now attempt to choose a default data source.")
        sql_data_source = NULL
      } else {
        sql_data_source = temp_source
      }
    } else if (is.numeric(sql_data_source)){
      if (!(sql_data_source %in% data_sources$source_key)) {
        cat("Warning: Numeric value entered into data_source does not match any DICE source_keys.  DICE will now attempt to choose a default data source.")
        sql_data_source = NULL
      }
    } else {
      cat("Warning: data_source is neither a single string or an integer.  DICE will now attempt to choose a default data source.")
      sql_data_source = NULL
    }
  }

  # currently downloads all years and all fit regions
  all_years = get_disease_data(mod_level=mod_level, fit_level=fit_level, mod_name=mod_name, years="all", disease=disease, sql_data_source=sql_data_source, daily_clim=FALSE, db_opts=db_opts, cad_clim=cad_clim)
  # if alternate data sources are selected, check them now to overwrite/augment all_years
  if (any(c(db_opts$CDC_server))) {
    if ((tolower(sql_data_source)=="cdc" | sql_data_source==7) &db_opts$CDC_server) {
      all_years = get.cdc(all_years=all_years, mod_level=mod_level, fit_level=fit_level, mod_name=mod_name, myDB=myDB)
    }
  }
  # use sql_data_source to determine source_abbv and cadence. if user did not specify a data source, pick one for them
  if (length(sql_data_source)!=0) {
    if (is.numeric(sql_data_source)) {
      data_sources_ind = sql_data_source==data_sources$source_key
    } else {
      data_sources_ind = tolower(sql_data_source) == tolower(data_sources$source_abbv)
    }
  } else {
    # search all_years for available data_sources and then pick one
    all_sources = character()
    for(ii in 2:length(all_years)) {  # loop through cadences
      cadence_sources = data_sources$source_abbv
      for(jj in 1:length(all_years[[ii]])) {  # loop through spatial identifiers
        cadence_sources = cadence_sources[cadence_sources %in% names(all_years[[ii]][[jj]])]
        all_sources = union(all_sources, cadence_sources)
      }
    }
    if (length(all_sources)==0) {
      stop("There are no data sources for the spatial patch entered.\n")
    } else if (length(all_sources)==1) {
      sql_data_source = data_sources$source_key[all_sources==data_sources$source_abbv]
      data_sources_ind = sql_data_source==data_sources$source_key
    } else if (length(all_sources)>1) {
      sql_data_source = min(which(data_sources$source_abbv %in% all_sources))
      data_sources_ind = sql_data_source==data_sources$source_key
      cat("Warning: No sql_data_source was specified and mulitple data sources exist (", paste0(data_sources$source_key[which(data_sources$source_abbv %in% all_sources)], collapse=", "), ") for the spatial patch entered.  DICE is defaulting to data source ", sql_data_source, ". \n")
    }


  }
  # use source to determine cadence and source_abbv
  cadence  = data_sources$cadence[data_sources_ind]
  cad_name = CadenceNum2Name(cadences=cadence)
  source_abbv = data_sources$source_abbv[data_sources_ind]

  # extract column names, units, and aggregation methods from source
  data_col_names = strsplit(data_sources$col_names[data_sources_ind], ";")[[1]]
  data_unit_names = strsplit(data_sources$col_units[data_sources_ind], ";")[[1]]
  data_DICE_units = col_units$unit_type[match(x=data_unit_names, table=col_units$col_unit)]
  data_agg_method = unit_types$aggregate_method[match(x=data_DICE_units, table=unit_types$unit_type)]

  # determine which climate columns are needed
  clim_cols = clim_by_disease$clim_names[clim_by_disease$disease==disease]
  clim_cols = strsplit(clim_cols, ";")[[1]]

  # determine mod_level master_key
  if (mod_level>=2) {
    mod_ind = (toupper(mod_name["NAME_2"]) == all_years$lut$ABBV_2 | mod_name["NAME_2"] == all_years$lut$NAME_2) & mod_level == all_years$lut$level
    fit_ind = (toupper(mod_name["NAME_2"]) == all_years$lut$ABBV_2 | mod_name["NAME_2"] == all_years$lut$NAME_2) & fit_level == all_years$lut$level
    if (mod_level > 2) {
      for (lev in 3:mod_level) {
        mod_ind = mod_ind & (mod_name[paste0("NAME_", lev)] == all_years$lut[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_", lev)]) == all_years$lut[, paste0("ABBV_", lev)])
        fit_ind = fit_ind & (mod_name[paste0("NAME_", lev)] == all_years$lut[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_", lev)]) == all_years$lut[, paste0("ABBV_", lev)])
      }
    }
    # determine if an aggregate has been specified
    if (!is.null(fit_names) && (length(fit_names)>1 || fit_names!="all")) {
      mod_sum = TRUE
      # allow fit regions to be specified by NAME_x, ABBV_x, or identifier
      fit_ind = fit_ind & (tolower(all_years$lut[, paste0("NAME_", fit_level)]) %in% tolower(fit_names) | tolower(all_years$lut[, paste0("ABBV_", fit_level)]) %in% tolower(fit_names) | tolower(all_years$lut$identifier) %in% tolower(fit_names))
    } else {
      mod_sum = FALSE
    }
  } else if (mod_level==1) {
    # narrow fit_ind to mod_name children at fit_level
    fit_ind = (toupper(mod_name["NAME_1"]) == all_years$lut$ABBV_1 | mod_name["NAME_1"] == all_years$lut$NAME_1) & fit_level == all_years$lut$level
    # Continent level $model will be constructed as a sum of $fit regions
    mod_sum = TRUE
    if (!is.null(fit_names) && (length(fit_names)>1 || fit_names!="all")) {
      # fit regions can be specified by NAME_x, ABBV_x, or identifier
      fit_ind = fit_ind & (tolower(all_years$lut[, paste0("NAME_", fit_level)]) %in% tolower(fit_names) | tolower(all_years$lut[, paste0("ABBV_", fit_level)]) %in% tolower(fit_names) | tolower(all_years$lut$identifier) %in% tolower(fit_names))
    } else if (fit_level==1) {
      # sum all available countries
      fit_ind = (toupper(mod_name["NAME_1"]) == all_years$lut$ABBV_1 | mod_name["NAME_1"] == all_years$lut$NAME_1) & 2 == all_years$lut$level
    } else {
      # determine indices of fit_level indices
      fit_ind = (toupper(mod_name["NAME_1"]) == all_years$lut$ABBV_1 | mod_name["NAME_1"] == all_years$lut$NAME_1) & fit_level == all_years$lut$level
    }

  } else {
    # Global level mydata
    # no reason to support this at this time
  }

  # Check that model-patch exists in MySQL database
  if (!mod_sum && all(!mod_ind)) {
    stop("model name/level combination not present in MySQL '", disease,"_data' table.")
  }
  if (all(!fit_ind)) {
    stop("fit-level mydata is not present for model name/level in MySQL '", disease,"_data' table.")
  }


  # Now that the source is concretely determined, check if 'season' was defined by the user
  if (mod_sum) {
    continent = all_years$lut$ABBV_1[1]
    if (mod_level<=1) {
      country = NULL
    } else if (mod_level>1) {
      country = all_years$lut$ABBV_2[1]
    }
  } else {
    continent = all_years$lut$ABBV_1[mod_ind]
    country   = all_years$lut$ABBV_2[mod_ind]
  }

  all_SE_dates = get_season_dates(myDB=myDB, season=season, continent=continent, country=country, cadence=cadence, disease=disease)
  if (!is.null(season) & mod_level>0) {
    # if so, use season to determine start/end dates
    SE_dates = all_SE_dates[all_SE_dates$season==season, ]
    start.date = as.Date(SE_dates$start_date)
    end.date   = as.Date(SE_dates$end_date)
  } else if (is.null(season)) {
    season = year(start.date)
  }


  if (!mod_sum) {
    # build 'model' list
    mod_ident = all_years$lut$identifier[mod_ind]
    mydata$model = list()
    # index dates
    if (is.null(start.date)) {
      if (is.null(end.date)) {
        # if start and end dates are null, grab all mydata
        date_ind = rep(TRUE, length(all_years[[cad_name]][[]][[source_abbv]]$date))
        # Generate date/year/etc vectors
        mydata = BuildDateVecs(cadence=cadence, start.date=min(all_years[[cad_name]][[mod_ident]][[source_abbv]]$date), end.date=max(all_years[[cad_name]][[mod_ident]][[source_abbv]]$date))
      } else {
        # if start.date is null, grab all mydata before end.date
        date_ind = all_years[[cad_name]][[mod_ident]][[source_abbv]]$date<=end.date
        # Generate date/year/etc vectors
        mydata = BuildDateVecs(cadence=cadence, start.date=min(all_years[[cad_name]][[mod_ident]][[source_abbv]]$date), end.date=end.date)
      }
    } else if (is.null(end.date)) {
      # if end.date is null, grab all mydata after start.date
      date_ind = all_years[[cad_name]][[mod_ident]][[source_abbv]]$date>=start.date
      # Generate date/year/etc vectors
      mydata = BuildDateVecs(cadence=cadence, start.date=start.date, end.date=max(all_years[[cad_name]][[mod_ident]][[source_abbv]]$date))
    } else {
      # if neither start or end dates are null, grab the interval between
      date_ind = all_years[[cad_name]][[mod_ident]][[source_abbv]]$date>=start.date & all_years[[cad_name]][[mod_ident]][[source_abbv]]$date<=end.date
      # Generate date/year/etc vectors
      mydata = BuildDateVecs(cadence=cadence, start.date=start.date, end.date=end.date)
    }
    # If using WHO_flu data, shift date to ISO weeks (week starts on Monday)
    if (sql_data_source==15) {
      mydata$date = mydata$date + 1
    }
    mydata$model$level = mod_level
    mydata$model$name  = all_years$lut[[paste0("NAME_",mod_level)]][mod_ind]
    mydata$model$attr  = all_years$lut[mod_ind, ]

    # grab mydata master date vecs
    if (cadence>3) {
      date_df = as.data.frame(mydata[1:2])
    } else {
      date_df = as.data.frame(mydata[1:3])
    }

    period_cols = names(date_df)[names(date_df) %in% c("year", "month", "week", "day")]

    # match available data with master date vecs
    temp_epi = all_years[[cad_name]][[mod_ident]][[source_abbv]]
    temp_epi = temp_epi[, names(temp_epi)[names(temp_epi)!="date"]]
    mydata$model[[mod_ident]] = merge(x=date_df, y=temp_epi, all.x=TRUE)

    # first add period columns to noaa_xxxx
    all_years[[cad_name]][[mod_ident]][[paste0("noaa_", tolower(cad_name))]] =  add_periods_cols(ident_dat=all_years[[cad_name]][[mod_ident]][[paste0("noaa_", tolower(cad_name))]], cadence=cadence)
    # then match climate data
    temp_clim = all_years[[cad_name]][[mod_ident]][[paste0("noaa_", tolower(cad_name))]]
    # remove date column so that matching is done only on year and day/week/month
    temp_clim = temp_clim[, names(temp_clim)!="date"]
    mydata$model[[mod_ident]] = merge(x=mydata$model[[mod_ident]], y=temp_clim, all.x=TRUE)
    # re-order by date
    mydata$model[[mod_ident]] = mydata$model[[mod_ident]][order(mydata$model[[mod_ident]]$date), ]
    # wrap-up
    mydata$model[[mod_ident]] = pluralize(ident_dat=mydata$model[[mod_ident]])
  }


  # build 'fit' list
  mydata$fit = list()
  mydata$fit$level = fit_level
  mydata$fit$name  = all_years$lut[[paste0("NAME_",fit_level)]][fit_ind]
  mydata$fit$attr  = all_years$lut[fit_ind, ]

  if (mod_sum) {
    # determine begining and end dates if start.date is NULL by finding the earliest and latest incidence data from fit-regions
    max_fit_date = as.Date(NA)
    min_fit_date = as.Date(NA)
    remove_ind = rep(FALSE, nrow(mydata$fit$attr))

    if (is.null(start.date)) {
      # determine earliest data point
      for (kk in 1:nrow(mydata$fit$attr)) {
        min_fit_date = min(min_fit_date, mydata$fit[[mydata$fit$attr$identifier[kk]]]$date, na.rm=T)
      }
      start.date = min_fit_date
    }

    if (is.null(end.date)) {
      # determine latest data point
      for (kk in 1:nrow(mydata$fit$attr)) {
        max_fit_date = max(max_fit_date, mydata$fit[[mydata$fit$attr$identifier[kk]]]$date, na.rm=T)
      }
      end.date = max_fit_date
    }
    # generate master date vector
    date_df = BuildDateVecs(cadence=cadence, start.date=start.date, end.date=end.date)
    mydata = c(as.list(date_df), mydata)
  }
  
  # Import fit region incidence
  for (kk in 1:nrow(mydata$fit$attr)) {
    # match available data with master date vecs
    fit_ident = mydata$fit$attr$identifier[kk]
    # take dates from date_df.  Match year and month/week/day columns
    temp_epi = all_years[[cad_name]][[fit_ident]][[source_abbv]]
    temp_epi = temp_epi[, names(temp_epi)[names(temp_epi)!="date"]]
    mydata$fit[[fit_ident]] = merge(x=date_df, y=temp_epi, all.x=TRUE)

    # first add period columns to noaa_xxxx
    if (mod_level!=fit_level) {
      all_years[[cad_name]][[fit_ident]][[paste0("noaa_", tolower(cad_name))]] =  add_periods_cols(ident_dat=all_years[[cad_name]][[fit_ident]][[paste0("noaa_", tolower(cad_name))]], cadence=cadence)
    }  # else this process was already performed for mod_level calcs

    # then match climate data
    temp_clim = all_years[[cad_name]][[fit_ident]][[paste0("noaa_", tolower(cad_name))]]
    # remove date column so that matching is done only on year and day/week/month
    temp_clim = temp_clim[, names(temp_clim)!="date"]
    
    mydata$fit[[fit_ident]] = merge(x=mydata$fit[[fit_ident]], y=temp_clim, all.x=TRUE)
    # re-order by date
    mydata$fit[[fit_ident]] = mydata$fit[[fit_ident]][order(mydata$fit[[fit_ident]]$date), ]
    # wrap-up
    mydata$fit[[fit_ident]] = pluralize(ident_dat=mydata$fit[[fit_ident]])
    mydata$fit[[fit_ident]] = add_ndays(ident_dat=mydata$fit[[fit_ident]], cadence=cadence)
  }

  ## This section tries to fill $pop with
  #     1. population from the incidence report
  #     2. census data
  #     3. then with the sedac_pop approximation

  # try SQL-table 'pop-yearly' first
  all_years_pop = get.census.data(master_keys=c(mydata$model$attr$master_key, mydata$fit$attr$master_key), myDB=myDB)
  # population retrieval for $fit
  mydata$fit$pop = rep(NA, nrow(mydata$fit$attr))
  for (ii in 1:nrow(mydata$fit$attr)) {
    ident = mydata$fit$attr$identifier[ii]
    master_key = mydata$fit$attr$master_key[ii]
    pop_ind = all_years_pop$master_key==master_key
    if ("pop" %in% names(mydata$fit[[ident]])) {
      # pull population from incidence data
      not_na_ind = which(!is.na(mydata$fit[[ident]]$pop))
      if (any(not_na_ind)) {
        use_ind = min(not_na_ind)
        mydata$fit$pop[ii] = mydata$fit[[ident]]$pop[use_ind]
      }
    }
    if (is.na(mydata$fit$pop[ii]) && any(pop_ind)) {
      # pull population from census data
      match_ind = season==year(all_years_pop$date[pop_ind])
      if (any(match_ind)) {
        mydata$fit$pop[ii] = all_years_pop$pop[pop_ind][match_ind]
      }
    }
    if (is.na(mydata$fit$pop[ii])) {
      # if first two options fail, use sedac_pop
      mydata$fit$pop[ii] = mydata$fit$attr$sedac_pop[ii]
    }
  } # end fit-level population loop

  if (!mod_sum) {
    mydata$model$pop = NA
    pop_ind = all_years_pop$master_key==mydata$model$attr$master_key
    if ("pop" %in% names(mydata$model[[mod_ident]])) {
      # pull population from incidence data
      not_na_ind = which(!is.na(mydata$model[[mod_ident]]$pop))
      if (any(not_na_ind)) {
        use_ind = min(not_na_ind)
        mydata$model$pop = mydata$model[[mod_ident]]$pop[use_ind]
      }
    }
    if (is.na(mydata$model$pop) && any(pop_ind)) {
      match_ind = season==year(all_years_pop$date[pop_ind])
      if (any(match_ind)) {
        mydata$model$pop = all_years_pop$pop[pop_ind][match_ind]
      } else {
        mydata$model$pop = NA
      }
    } else {
      mydata$model$pop = NA
    }
    # else recover SEDAC_pop from $attr
    if (is.na(mydata$model$pop)) {
      mydata$model$pop = mydata$model$attr$sedac_pop
    }
  } else {
    # sum fit-level populations
    model_pop = sum(mydata$fit$pop)
  }
  
  # use data_source to determine which data column will be 'raw'
  mydata = Proc_IncCol_BySource_mydata(mydata=mydata, data_source=sql_data_source, cadence=cadence, mod_sum=mod_sum, col_units=col_units, unit_types=unit_types, raw_col=raw_col, data_sources=data_sources)

  # check for regions with sum(raw)==0
  fit_idents = mydata$fit$attr$identifier
  remove_ind = rep(FALSE, nrow(mydata$fit$attr))
  for (kk in 1:nrow(mydata$fit$attr)) {
    if (sum(mydata$fit[[fit_idents[kk]]]$raw, na.rm=T)==0) {
      remove_ind[kk] = TRUE
      mydata$fit[[fit_idents[kk]]] = NULL
    }
  }
  if (all(remove_ind)) {
    stop("Fit_level region(s) either do not have data for the date-range or incidence data time-series is(are) zero-sum.  Please choose a different level/region/season combination.")
  }

  mydata$fit$attr = mydata$fit$attr[!remove_ind, ]
  mydata$fit$name = mydata$fit$name[!remove_ind]
  mydata$fit$factor = mydata$fit$factor[!remove_ind]
  mydata$fit$pop  = mydata$fit$pop[!remove_ind]
  fit_idents = mydata$fit$attr$identifier

  # if model will be a sum of fit regions, calc now
  if (mod_sum) {
    # build 'model' list
    mydata$model = list()
    mydata$model$level = mod_level
    mydata$model$name  = paste0(all_years$lut[[paste0("NAME_",mod_level)]][1], "-aggregate")
    # aggregate model_lut from remaining fit entries
    mydata$model$attr  = mydata$fit$attr[1,]
    mydata$model$attr[, paste0("NAME_", (mod_level+1):fit_level)] = NA
    mydata$model$attr[, paste0("ID_", (mod_level+1):fit_level)] = NA
    mydata$model$attr[, paste0("ABBV_", (mod_level+1):fit_level)] = NA

    mydata$model$attr$identifier = paste(mydata$model$attr[, paste0("ABBV_", 1:mod_level)], collapse=".")
    mydata$model$attr[, c("master_key", "inc_key", "gadm_name", "gadm_lvl",  "clim_ident", " gadm_noaa_sedac_ident")] = NA
    # sum area and population; average lat/lon
    mydata$model$attr$gadm_area = sum(mydata$fit$attr$gadm_area)
    mydata$model$attr$sedac_pop = sum(mydata$fit$attr$sedac_pop)
    mydata$model$attr$gadm_lat  = sum(mydata$fit$attr$gadm_area*mydata$fit$attr$gadm_lat)/mydata$model$attr$gadm_area
    mydata$model$attr$gadm_lon  = sum(mydata$fit$attr$gadm_area*mydata$fit$attr$gadm_lon)/mydata$model$attr$gadm_area
    mydata$model$attr$sedac_lat  = sum(mydata$fit$attr$sedac_pop*mydata$fit$attr$sedac_lat)/mydata$model$attr$sedac_pop
    mydata$model$attr$sedac_lon  = sum(mydata$fit$attr$sedac_pop*mydata$fit$attr$sedac_lon)/mydata$model$attr$sedac_pop

    # aggregate fit patches
    # first determine all dates needed
    const_cols = c("years", "weeks", "months", "days")
    # sum_cols   = c("cases", "ili_out", "pop", "cases_B", "cases_dhf", "deaths")
    sum_cols     = data_col_names[data_agg_method=="sum"]
    # pop_avg_cols = c("ILI", "ili_rate", "sh", "precip", "temp")
    pop_avg_cols = c(data_col_names[data_agg_method=="pop_avg"], clim_cols)
    fit_agg_method = data_agg_method[mydata$fit$raw_units==data_unit_names]
    if (fit_agg_method=="sum") {
      sum_cols = c(sum_cols, "raw")
    } else if (fit_agg_method=="pop_avg") {
      pop_avg_cols = c(pop_avg_cols, "raw")
    } else {
      stop("Aggregation method '", fit_agg_method, "' is not yet supported.\n", sep="")
    }

    dates = mydata$date

    mod_ident = mydata$model$attr$identifier
    mydata$model[[mod_ident]] = data.frame(dates=dates, years=NA)

    fit_ind = integer(length(mydata$fit$attr$identifier))

    col_names = names(mydata$fit[[fit_idents[1]]])
    # const_cols
    const_names = col_names[col_names %in% const_cols]
    sum_names   = col_names[col_names %in% sum_cols]
    pop_avg_names = col_names[col_names %in% pop_avg_cols]
    for (ii in 1:length(dates)) {
      date = dates[ii]
      for(jj in 1:length(fit_ind)) {
        temp_ind = which(mydata$fit[[fit_idents[jj]]]$dates==date)
        if (length(temp_ind)==1) {
          fit_ind[jj] = temp_ind
        } else {
          fit_ind[jj] = NA
        }
      }
      # aggregate constant columns (choose first region available)
      ident_ind = which(!is.na(fit_ind))[1]
      for (name in const_names) {
        mydata$model[[mod_ident]][[name]][ii] = mydata$fit[[fit_idents[ident_ind]]][[name]][fit_ind[ident_ind]]
      }
      # aggreagate summable columns (cases, deaths, etc)
      ident_ind = which(!is.na(fit_ind))
      for (name in sum_names) {
        temp_sum = NA
        for (ind in ident_ind) {
          ident = fit_idents[ind]
          add_x = mydata$fit[[ident]][[name]][fit_ind[ind]]
          if (!is.na(add_x)) {
            temp_sum = sum(temp_sum, add_x, na.rm=T)
          }
        }
        mydata$model[[mod_ident]][[name]][ii] = temp_sum
      }
      # aggreagate pop_avg columns (ILI rate, climate, etc)
      for (name in pop_avg_names) {
        temp_sum = NA
        for (ind in ident_ind) {
          ident = fit_idents[ind]
          add_x = mydata$fit[[ident]][[name]][fit_ind[ind]] * mydata$fit$attr$sedac_pop[ind]
          if (!is.na(add_x)) {
            temp_sum = sum(temp_sum, add_x, na.rm=T)
          }
        }
        mydata$model[[mod_ident]][[name]][ii] = temp_sum/sum(mydata$fit$attr$sedac_pop[ident_ind])
      }
    }
    mydata$model$pop       = model_pop
    mydata$model$raw_units = mydata$fit$raw_units
    # mydata$model$factor    = sum(mydata$fit$factor*mydata$fit$attr$sedac_pop)/sum(mydata$fit$attr$sedac_pop)
    mydata$model$factor    = sum(mydata$fit$factor*mydata$fit$pop)/sum(mydata$fit$pop)

  }

  # change date column names to plurals
  mydata = pluralize(ident_dat=mydata)

  # calculate 'ndays' for all time-series
  mydata$model[[mod_ident]] = add_ndays(ident_dat=mydata$model[[mod_ident]], cadence=cadence)
  if (cadence %in% 1:3) {
    mydata = c(mydata[1:3], list(ndays=mydata$model[[mod_ident]]$ndays), mydata[4:length(mydata)])
  } else {
    mydata = c(mydata[1:2], list(ndays=mydata$model[[mod_ident]]$ndays), mydata[4:length(mydata)])
  }

  # for (fit_ident in mydata$fit$attr$identifier) {
  #   mydata$fit[[fit_ident]] = add_ndays(ident_dat=mydata$fit[[fit_ident]], cadence=cadence)
  # }

  mydata$nperiods = length(mydata$dates)

  ##mydata$nperiodsData = sum(!is.na(mydata$model[[mod_ident]]$raw))
  mydata$nperiodsData =  max(which(!is.na(mydata$model[[mod_ident]]$raw)))

  # check if any NAs appear in climate data and fill with historic average
      # fit
  for (fit_ident in mydata$fit$attr$identifier) {
    # then fill NA climate entries
    mydata$fit[[fit_ident]] = AvgClim(ident_data=mydata$fit[[fit_ident]], clim_cols=clim_cols, all_years_ident=all_years[[cad_name]][[fit_ident]][[paste0("noaa_", tolower(cad_name))]], cadence=cadence)
  }
    # model
  if (!mod_sum) {
    if (fit_level>mod_level) {
      # then fill NA climate entries
      mydata$model[[mod_ident]] = AvgClim(ident_data=mydata$model[[mod_ident]], clim_cols=clim_cols, all_years_ident=all_years[[cad_name]][[mod_ident]][[paste0("noaa_", tolower(cad_name))]], cadence=cadence)
    }
  } else {
    # population average fit regions
    clim_names = col_names[col_names %in% clim_cols]
    for (col_name in clim_names) {
      na_ind = which(is.na(mydata$model[[mod_ident]][col_name]))
      for (ii in na_ind) {
        temp_sum = NA
        for (jj in 1:length(fit_idents)) {
          ident = fit_idents[jj]
          fit_ind = mydata$fit[[ident]]$dates==mydata$model[[mod_ident]]$dates[ii]
          add_x = mydata$fit[[ident]][[col_name]][fit_ind]*mydata$fit$attr$sedac_pop[jj]
          if (length(add_x)>0 && !is.na(add_x)) {
            temp_sum = sum(temp_sum, add_x, na.rm=T)
          }
        }
        mydata$model[[mod_ident]][[col_name]][ii] = temp_sum/sum(mydata$fit$attr$sedac_pop)
      }
    }
  }
  
  mydata$fit$nregions = length(mydata$fit$attr$identifier)

  data_units = data_sources$col_units[data_sources_ind]
  data_units = strsplit(data_units, split=";")[[1]]
  # determine which 'mydataX' column to use, then convert to 'epi' values as input to DICE


  # retrieve all years of school data
  all_years_school = get.school.data(master_keys=c(mydata$model$attr$master_key, mydata$fit$attr$master_key), cadence=cadence, myDB=myDB)

  # assign fit-level school data
  for (ii in 1:mydata$fit$nregions) {
    ident = mydata$fit$attr$identifier[ii]
    master_key = mydata$fit$attr$master_key[ii]
    school_ind = all_years_school$master_key==master_key
    if (any(school_ind)) {
      match_ind = match(mydata$dates, all_years_school$date[school_ind])
      mydata$fit[[ident]]$school = all_years_school$school[school_ind][match_ind]
    } else {
      mydata$fit[[ident]]$school = rep(NA, mydata$nperiods)
    }
  }

  # assign model-level school data
  if (!mod_sum) {
    school_ind = all_years_school$master_key==mydata$model$attr$master_key
    if (any(school_ind)) {
      # import school schedule into mydata
      match_ind = match(mydata$dates, all_years_school$date[school_ind])
      mydata$model[[mydata$model$attr$identifier]]$school = all_years_school$school[school_ind][match_ind]
    } else {
      mydata$model[[mydata$model$attr$identifier]]$school = rep(NA, mydata$nperiods)
    }
  } else {
    # aggregate model school schedule from fit schedules
    temp_school = matrix(NA, nrow=mydata$fit$nregions, ncol=mydata$nperiods)
    pop_date    = integer(mydata$nperiods)
    for (ii in 1:mydata$fit$nregions) {
      ident = mydata$fit$attr$identifier[ii]
      temp_school[ii, ] = mydata$fit[[ident]]$school * mydata$fit$pop[ii]
      na_ind = is.na(temp_school[ii, ])
      pop_date[!na_ind] = pop_date[!na_ind] + mydata$fit$pop[ii]
    }
    mydata$model[[mydata$model$attr$identifier]]$school = apply(temp_school, 2, sum, na.rm=T)/pop_date
  }

  # incorporate school data into all_years
  for (ii in 1:length(all_years$lut$identifier)) {
    ident = all_years$lut$identifier[ii]
    key   = all_years$lut$master_key[ii]
    # grab data from epi/climate structure for master_key
    temp_epi = all_years[[cad_name]][[ident]][[source_abbv]]
    # grab school data for master_key
    if (nrow(all_years_school)>0) {
      temp_school = all_years_school[all_years_school$master_key==key, c("date", "school")]
    } else {
      # PostGRES empty queries return an empty data.frame without columns, so create an empty data.frame with columns
      temp_school = data.frame(date=character(), school=integer())
    }
    # merge, keeping all dates from both data types
    temp_epi = merge(x=temp_epi, y=temp_school, all.x=TRUE, sort=F)
    # re-order by date
    temp_epi = temp_epi[order(temp_epi$date), ]
    # re-assign to all_years
    all_years[[cad_name]][[ident]][[source_abbv]] = temp_epi
  }


  # Now attempt to retrieve onset value
  all_onset = get.onset.data(master_keys=c(mydata$model$attr$master_key, mydata$fit$attr$master_key), start.year=season, end.year=season, myDB=myDB)

  if (nrow(all_onset)>0) {
    if (!mod_sum) {
      onset_ind = mydata$model$attr$master_key==all_onset$master_key
      if (any(onset_ind)) {
        mydata$model$onset = all_onset$onset[onset_ind]
      } else {
        mydata$model$onset = NA
      }
    } else {
      # model-region is the result of aggregation, so there has been no onset value assigned to it.
      mydata$model$onset = NA
    }

    # repeat onset for $fit
    mydata$fit$onset = numeric(mydata$fit$nregions)
    for (ii in 1:mydata$fit$nregions) {
      onset_ind = mydata$fit$attr$master_key[ii]==all_onset$master_key
      if (any(onset_ind)) {
        mydata$fit$onset[ii] = all_onset$onset[onset_ind]
      } else {
        mydata$fit$onset[ii] = NA
      }
    }
  }


  # determine if coef is 1 or population-weighted
  fit_agg_method = data_agg_method[mydata$fit$raw_units==data_unit_names]
  if (fit_agg_method=="pop_avg") {
    mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)
  } else if (fit_agg_method=="sum") {
    mydata$fit$coef = rep(1, mydata$fit$nregions)
  } else {
    stop("Aggregation method '", fit_agg_method, "' is not yet supported.\n", sep="")
  }

  # Convert 'raw' to epi and gamaepi
  mydata$model[[mod_ident]]$epi = mydata$model[[mod_ident]]$raw*mydata$model$factor
  # convert NAs to zero and floats to integers
  mydata$model[[mod_ident]]$epi[is.na(mydata$model[[mod_ident]]$epi)] = 0
  mydata$model[[mod_ident]]$epi = as.integer(round(mydata$model[[mod_ident]]$epi))
  mydata$model[[mod_ident]]$gamaepi = lgamma(mydata$model[[mod_ident]]$epi+1)

  for (ii in 1:length(mydata$fit$attr$identifier)) {
    fit_ident = mydata$fit$attr$identifier[ii]
    mydata$fit[[fit_ident]]$epi = mydata$fit[[fit_ident]]$raw*mydata$fit$factor[ii]
    # convert NAs to zero and floats to integers
    mydata$fit[[fit_ident]]$epi[is.na(mydata$fit[[fit_ident]]$epi)] = 0
    mydata$fit[[fit_ident]]$epi = as.integer(round(mydata$fit[[fit_ident]]$epi))
    mydata$fit[[fit_ident]]$gamaepi = lgamma(mydata$fit[[fit_ident]]$epi+1)
  }

  dbDisconnect(myDB)


  mydata$disease = disease
  mydata$cadence = cad_name
  mydata$data_source = tolower(data_sources$source_abbv[data_sources_ind])
  mydata$data_desc   = data_sources[data_sources_ind, ]

  # if all_years_flag, compile all_years_epi
  if (all_years_flag) {
    # reformat all_years and export for use in mydata-augmentation routines
    if (cadence %in% 1:3) {
      date_cols = 1:3
    } else if (cadence==4) {
      date_cols = 1:2
    }

    # create master vector of dates for all_years_epi
    min_date = as.Date(NA)
    max_date = as.Date(NA)

    if (!mod_sum) {
      min_date = min(min_date, all_years[[cad_name]][[mod_ident]][[source_abbv]]$date, na.rm=T)
      max_date = max(max_date, all_years[[cad_name]][[mod_ident]][[source_abbv]]$date, na.rm=T)
    }
    for(ii in 1:mydata$fit$nregions) {
      fit_ident = mydata$fit$attr$identifier[ii]
      min_date = min(min_date, all_years[[cad_name]][[fit_ident]][[source_abbv]]$date, na.rm=T)
      max_date = max(max_date, all_years[[cad_name]][[fit_ident]][[source_abbv]]$date, na.rm=T)
    }
    all_years_epi = BuildDateVecs(cadence=cadence, start.date=min_date, end.date=max_date)
    date_names = names(all_years_epi)[date_cols]
    all_master_dates = as.data.frame(all_years_epi)

    keep_cols = c("date", "year", "month", "ndays", "week", "day","raw", "epi", "sh", "temp", "precip", "press", "rh", "school")
    
    temp_all_years = list()
    for (ii in 1:mydata$fit$nregions) {
      fit_ident = mydata$fit$attr$identifier[ii]
      all_years[[cad_name]][[fit_ident]][[source_abbv]] = Proc_IncCol_BySource(ident_data=all_years[[cad_name]][[fit_ident]][[source_abbv]], data_source=sql_data_source, factor=mydata$fit$factor[ii], raw_col=mydata$fit$raw_name)
      all_years_cols = keep_cols[keep_cols %in% names(all_years[[cad_name]][[fit_ident]][[source_abbv]])]
      for (jj in (max(date_cols)+1):length(all_years_cols)) {
        # prepare a temp_epi structure to merge with all_years_epi$fit.  This automates date-matching and data.frame creation for regions with different length histories
        temp_epi  = cbind(all_years[[cad_name]][[fit_ident]][[source_abbv]][, date_cols], all_years[[cad_name]][[fit_ident]][[source_abbv]][, all_years_cols[jj]])
        names(temp_epi)[max(date_cols)+1] = fit_ident
        exist_epi = temp_all_years[[all_years_cols[jj]]]
        if (!is.null(exist_epi)) {
          temp_all_years[[all_years_cols[jj]]] = merge(exist_epi, temp_epi, all=TRUE)
        } else {
          temp_all_years[[all_years_cols[jj]]] = temp_epi
        }

      }
    }

    # ensure that fit data frames contain all master dates
    for (var_name in names(temp_all_years)) {
      temp_all_years[[var_name]] = merge(x=all_master_dates[, ], y=temp_all_years[[var_name]][, 2:length(temp_all_years[[var_name]])], all.x=T, sort=FALSE)
      temp_all_years[[var_name]] = temp_all_years[[var_name]][order(temp_all_years[[var_name]]$date), ]
    }

    # append appropriate climate data from noaa_xxxxxx
    for (clim in clim_cols) {
      temp_all_years[[clim]] = all_master_dates
      for (ident in mydata$fit$attr$identifier) {
        # separate each identifier
        temp_clim = all_years[[cad_name]][[ident]][[paste0("noaa_", tolower(cad_name))]]
        temp_clim = temp_clim[, names(temp_clim) %in% c(date_names, clim)]
        # remove 'date' column.  This dictates that merge() matches by year and period rather than date.
        temp_clim = temp_clim[, names(temp_clim)!="date"]
        # rename for merging
        names(temp_clim)[names(temp_clim)==clim] = ident
        # merge with data.frame
        temp_all_years[[clim]] = merge(x=temp_all_years[[clim]], y=temp_clim, all.x=TRUE, sort=FALSE)
      }
      # re-sort by date
      temp_all_years[[clim]] = temp_all_years[[clim]][order(temp_all_years[[clim]]$date), ]
    }


    all_years_epi$fit = list()
    for (var_name in names(temp_all_years)) {
      all_years_epi$fit[[var_name]] = temp_all_years[[var_name]][, (max(date_cols)+1):length(temp_all_years[[var_name]])]
    }
    # if any additional NAs are introduced to $epi, replace now with 0
    all_years_epi$fit$epi[is.na(all_years_epi$fit$epi)] = 0

    if (!mod_sum) {
      # add epi vector
      all_years[[cad_name]][[mod_ident]][[source_abbv]] = Proc_IncCol_BySource(ident_data=all_years[[cad_name]][[mod_ident]][[source_abbv]], data_source=sql_data_source, factor=mydata$model$factor, raw_col=mydata$model$raw_name)
      # ensure that all master dates exist
      all_years[[cad_name]][[mod_ident]][[source_abbv]] = merge(x=all_master_dates, y=all_years[[cad_name]][[mod_ident]][[source_abbv]][, names(all_years[[cad_name]][[mod_ident]][[source_abbv]])!="date"], all.x=T)
      all_years[[cad_name]][[mod_ident]][[source_abbv]] = all_years[[cad_name]][[mod_ident]][[source_abbv]][order(all_years[[cad_name]][[mod_ident]][[source_abbv]]$date), ]
      # columns to keep for model
      all_years_cols = keep_cols[keep_cols %in% names(all_years[[cad_name]][[mod_ident]][[source_abbv]])]
      all_years_cols = all_years_cols[-date_cols]
      # move into all_years_epi
      all_years_epi$model = all_years[[cad_name]][[mod_ident]][[source_abbv]][, all_years_cols]
      # append appropriate climate data
      temp_clim = all_master_dates
      temp_noaa = all_years[[cad_name]][[mod_ident]][[paste0("noaa_", tolower(cad_name))]]
      # remove 'date' column.  Force merge to use year and period for matching
      temp_noaa = temp_noaa[, names(temp_noaa)!="date"]
      temp_clim = merge(x=temp_clim, y=temp_noaa, all.x=T, sort=F)
      # ensure that rows are ordered by date
      temp_clim = temp_clim[order(temp_clim$date), ]
      # combine epi and climate
      all_years_epi$model = cbind(all_years_epi$model, temp_clim[, clim_cols])
      # all_years_epi$model = pluralize(all_years_epi$model)
      # all_years_epi$model = add_ndays(ident_dat=all_years_epi$model, cadence=cadence)
    } else {
      # aggregate from fit
      all_years_epi$model = data.frame(raw=numeric(nrow(all_years_epi$fit$raw)))
      # calc $raw aggregate based on raw_units
      if (fit_agg_method=="sum") {
        # sum fit regions
        all_years_epi$model$raw = apply(all_years_epi$fit$raw, MARGIN=1, FUN=sum, na.rm=T)
      } else if (fit_agg_method=="pop_avg") {
        # population-weighted average fit regions
        all_years_epi$model$raw = apply(all_years_epi$fit$raw, MARGIN=1, FUN=function(x,y,na.rm) sum(x*y, na.rm=na.rm)/sum(y), na.rm=T, y=mydata$fit$pop)
      } else {
        stop("Aggregation method '", fit_agg_method, "' is not yet supported.\n", sep="")
      }

      # check for rows where all fit regions have NA
      NA_rows = apply(all_years_epi$fit$raw, MARGIN=1, FUN=function(x) all(is.na(x)))
      all_years_epi$model$raw[NA_rows] = NA

      # sum epi vectors from fit
      all_years_epi$model$epi = apply(all_years_epi$fit$epi, MARGIN=1, sum, na.rm=T)
      # population average climate and school factors
      for (col_name in names(temp_all_years)[names(temp_all_years) %in% c("sh", "temp", "precip", "school")]) {
        for (ii in 1:nrow(all_years_epi$model)) {
          use_ind = !is.na(all_years_epi$fit[[col_name]][ii,])
          if (any(use_ind)) {
            all_years_epi$model[[col_name]][ii] = sum(all_years_epi$fit[[col_name]][ii, use_ind]*mydata$fit$pop[use_ind], na.rm=T)/sum(mydata$fit$pop[use_ind])
          } else {
            all_years_epi$model[[col_name]][ii] = NA
          }
        }
      }
    }
    all_years_epi = pluralize(all_years_epi)
  } else {
    all_years_epi = NULL
  }

  if (cad_clim) {
    all_years_clim = list()
    if (cadence %in% 1:3) {
      date_cols = 1:3
    } else if (cadence==4) {
      date_cols = 1:2
    }

    # import all fit regions
    noaa_cad = paste0("noaa_", tolower(cad_name))
    fit_ident = mydata$fit$attr$identifier[1]
    total_df = list()
    for (clim in clim_cols) {
      total_df[[clim]] = cbind(all_years[[cad_name]][[fit_ident]][[noaa_cad]][, date_cols], all_years[[cad_name]][[fit_ident]][[noaa_cad]][, clim])
      names(total_df[[clim]])[max(date_cols)+1] = fit_ident
    }

    if (nrow(mydata$fit$attr)>1) {
      for (ii in 2:nrow(mydata$fit$attr)) {
        fit_ident = mydata$fit$attr$identifier[ii]
        for (clim in clim_cols) {
          temp_clim = cbind(all_years[[cad_name]][[fit_ident]][[noaa_cad]][, date_cols], all_years[[cad_name]][[fit_ident]][[noaa_cad]][, clim])
          names(temp_clim)[max(date_cols)+1] = fit_ident
          total_df[[clim]] = merge(x=total_df[[clim]], y=temp_clim, all=TRUE)
        }
      }
    }

    mod_ident = mydata$model$attr$identifier
    if (mod_sum) {
      # population-average fit regions
      for (clim in clim_cols) {
        temp_clim = total_df[[clim]][, (max(date_cols)+1):ncol(total_df[[clim]])]
        for (ii in 1:ncol(temp_clim)) {
          temp_clim[, ii] = temp_clim[, ii]*mydata$fit$pop[ii]
        }
        total_df[[clim]][[mod_ident]] = apply(temp_clim, MARGIN=1, FUN=sum, na.rm=TRUE)/sum(mydata$fit$pop)
      }
    } else {
      # pull from all_years
      for (clim in clim_cols) {
        temp_clim = cbind(all_years[[cad_name]][[mod_ident]][[noaa_cad]][, date_cols], all_years[[cad_name]][[mod_ident]][[noaa_cad]][, clim])
        names(temp_clim)[max(date_cols)+1] = mod_ident
        total_df[[clim]] = merge(x=total_df[[clim]], y=temp_clim, all=TRUE)
      }
    }

    # Now that all climate data is consolidated (and all dates are included), create all_years_clim to match the format of all_years_epi.
    date_col_names = names(total_df[[1]])[date_cols]
    for(col_name in date_col_names) {
      all_years_clim[[col_name]] = total_df[[1]][[col_name]]
    }
    # construct $model data.frame
    temp_clim = data.frame(total_df[[clim_cols[1]]][[mod_ident]])
    names(temp_clim) = clim_cols[1]
    if (length(clim_cols)>1) {
      for (ii in 2:length(clim_cols)) {
        temp_clim[[clim_cols[ii]]] = total_df[[clim_cols[ii]]][[mod_ident]]
      }
    }
    all_years_clim$model = temp_clim

    # construct $fit list
    all_years_clim$fit = list()
    for (clim in clim_cols) {
      all_years_clim$fit[[clim]] = total_df[[clim]][, mydata$fit$attr$identifier]
    }
    all_years_clim = pluralize(all_years_clim)
  } else {
    all_years_clim = NULL
  }

  mysql_out = list(mydata=mydata, all_years_epi=all_years_epi, all_years_clim=all_years_clim, all_SE_dates=all_SE_dates)
  return(mysql_out)
}


get.sh.nasa <- function(region = "all", sub_region = 1, start_year = 2000, start_week = 1, end_year = 2016, end_week = 10) {
    #' Retrieve specific humidity data from data-server
    #'
    #' Data is retrieved for the United States, HHS Regions, and States.
    #' Attempt to pull CDC data from the CDC server. If that fails, pull from diceData.
    #'
    #' @param region String specifying USA ('national'), CDC Region ('hhs'), US state ('state'), or all mydata regions ('all').
    #' @param sub_region For region='hhs': Integer(s) indicating which of the 10 CDC regions are to be retrieved.  For region='state': a vector of state names to be returned.
    #' @param start_year An integer - start year of the flu season
    #' @param end_year  An integer  - end year of the flu season
    #' @param start_week An integer - starting CDC week of the flu season (default is 27)
    #' @param end_week   An integer - ending CDC week of the flu season (default is 26)
    #' @return a data frame with columns for $date, $week, $year, and whichever regions were requested.
    #' @examples
    #' require(DICE)
    #' sh = get.sh.nasa(region='hhs',sub_region=c(1,9),start_year=2015,end_year=2016,start_week=27,end_week=26)
    #' sh2 = get.sh.nasa(region='state',sub_region=c('District.of.Columbia','Kansas'),start_year=2015,end_year=2016,start_week=27,end_week=26)

    if (length(region) != 1)
        stop("'region' must have exactly one string")
    if (!any(region == c("all", "national", "hhs", "state")))
        stop("'region' must be one of the following strings: 'all', 'national', 'hhs', or 'state'")
    if (region == "hhs" && length(setdiff(sub_region, 1:10)) != 0)
        stop("for region='hhs', sub_region must be an integer or vector of integers between 1 and 10.")

    # require(RMySQL) create connection to MySQL mydatabase user = 'epi_guest' userp = 'UY5GE2kfUa' host = 'shadow.predsci.com'
    # myDB = dbConnect(MySQL(), user = "epi_guest", password = "UY5GE2kfUa", dbname = "epi_data", host = "shadow.predsci.com")
    myDB = OpenCon(list(DICE_db="PredSci"))

    # Determine date range
    startdate = CDCweek2date(start_week, start_year)
    enddate = CDCweek2date(end_week, end_year)

    if (region == "all") {
        SHquery = dbSendQuery(myDB, paste0("SELECT * from SHData where date>='", startdate, "' and date<='", enddate, "';"))
        sh = dbFetch(SHquery, n = -1)
        dbClearResult(SHquery)
    } else if (region == "national") {
        SHquery = dbSendQuery(myDB, paste0("SELECT date,week,year,national from SHData where date>='", startdate, "' and date<='", enddate,
            "';"))
        sh = dbFetch(SHquery, n = -1)
        dbClearResult(SHquery)
    } else if (region == "hhs") {
        SHquery = dbSendQuery(myDB, paste0("SELECT date,week,year,", paste0("Region", sub_region, collapse = ","), " from SHData where date>='",
            startdate, "' and date<='", enddate, "';"))
        sh = dbFetch(SHquery, n = -1)
        dbClearResult(SHquery)
    } else if (region == "state") {
        SHquery = dbSendQuery(myDB, paste0("SELECT date,week,year,", paste0("`", sub_region, "`", collapse = ","), " from SHData where date>='",
            startdate, "' and date<='", enddate, "';"))
        sh = dbFetch(SHquery, n = -1)
        dbClearResult(SHquery)
    }
    dbDisconnect(myDB)

    sh$date = as.Date(sh$date)
    sh$week = as.integer(sh$week)
    sh$year = as.integer(sh$year)

    return(sh)
}


create_patch <- function(mydata = NULL, RegState = NULL, fit_level = 3, mod_level = 2) {
    #' Create a New Region
    #'
    #' Using the populations specified in `RegState` create a new region. If the user lists multiple entries to `RegState', this function will aggregate those regions into a new super-region.
    #' @param mydata A mydata structure resulting from get.cdc.data() that contains all the regions listed in `RegState' in the $fit sub structure.
    #' @param RegState A vector of strings specifying either the name or abbreviation of US regions or states.
    #' @param fit_level An integer specifying the spatial-level of the elements of RegState. All elements of RegState must be from the same level.
    #' @param mod_level An integer that defines the spatial-level of the created region.
    #' @return An updated version of `mydata' with only the RegState elements in $fit and the aggregated super-region in $model
    #' @examples
    #' require(DICE)
    #' mydata = get.DICE.data()
    #' mydata2 = create_patch(mydata=mydata,RegState=c('R1','R2','R3'),fit_level=mydata$fit$level,mod_level=2)

    # Delete entries to $fit that are not contained in RegState
    keep_ind = (mydata$fit$attr[, paste0("NAME_", fit_level)] %in% RegState) | (mydata$fit$attr[, paste0("ABBV_", fit_level)] %in% toupper(RegState))

    if (!any(keep_ind)) {
        stop("No level ", fit_level, " entries match the names provided in RegState.")
    }

    mydata$fit$name = mydata$fit$name[keep_ind]
    mydata$fit$attr = mydata$fit$attr[keep_ind, ]
    mydata$fit$school = mydata$fit$school[, keep_ind]
    mydata$fit$pop = mydata$fit$pop[keep_ind]
    mydata$fit$onset = mydata$fit$onset[keep_ind]
    mydata$fit$coef = mydata$fit$pop/sum(mydata$fit$pop)
    mydata$fit$raw = mydata$fit$raw[, keep_ind]
    mydata$fit$raw_units = mydata$fit$raw_units[keep_ind]
    mydata$fit$factor = mydata$fit$factor[keep_ind]
    mydata$fit$epi = mydata$fit$epi[, keep_ind]
    mydata$fit$gamaepi = mydata$fit$gamaepi[, keep_ind]
    mydata$fit$sh = mydata$fit$sh[, keep_ind]
    mydata$fit$nregions = sum(keep_ind)

    # Use $fit data to construct a new $model patch
    mydata$model$level = mod_level
    mydata$model$name = paste0(mydata$fit$attr[, paste0("ABBV_", mydata$fit$level)], collapse = ".")
    # Currently just returning the input $model$attr
    mydata$model$attr = mydata$model$attr
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
    # newly calculated fields mydata$model$epi = apply(mydata$fit$epi,1,sum)
    mydata$model$epi = round(mydata$model$raw * mydata$model$factor)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)


    return(mydata)
}


get.FitLevelName <- function(iFit = 3) {
    #' Spatial Name for Fit Level
    #'
    #' Given a spatial Fit level this functions returns the corresponding spatial name. It is called by \code{\link{get.DICE.data}}. It assigns the correct spatial name to the spatial level
    #' @param iFit An integer (Currently support 2-4 for cdc and gft data)
    #' @return A character string with the corresponding spatial Name: earth, continent, country, hhs_region, state, county and city for iFit between 0 and 6
    #' @examples
    #' get.FitLevelName(iFit=3)
    #'
    #' get.FitLevelName(iFit=4)

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
    #' Spatial Name for Model Level
    #'
    #' Given a spatial Model level and the corresponding RegState value this function returns the corresponding spatial name
    #' @param mod_level An integer 2-4 for cdc and gft data
    #' @param RegState A string, USA, or in case of mod_level=2, or region number of mod_level=3
    #' @param NAME_4 Relevant for mod_level = 4 and up.  State name
    #' @param NAME_5 Relevant for mod_level = 5 - county name
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


get_disease_data <- function(mod_level=2, fit_level=3, mod_name=c(NAME_2="BR"), years="all", disease="dengue", sql_data_source=integer(), daily_clim=FALSE, cad_clim=TRUE, mod_sum=FALSE, db_opts=list(DICE_db="predsci", CDC_server=F)) {
  #' Retrieve Incidence and Climate Data
  #'
  #' This function retrieves incidence and climate data from the MySQL database. It is generally called by \code{\link{get.DICE.data}}.
  #' @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has mydata at levels 2-3 for CDC and 2-4 for GFT.
  #' @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  #' @param mod_name A named-vector of character strings that specify which region is to be \strong{model}ed.  In other words, \code{mod_name} specifies the country, region, state, etc. of the mod_level region.  \code{mod_name} should be of the form \code{mod_name = c(NAME_2='a', NAME_3='b',..., NAME_i='x'} where i=mod_level and 'a', 'b',...,'x' are the appropriate level names. NAME_i='x' also accepts abbreviations.  Choose appropriate names from \code{\link{diceData}}.  For example, mod_name=c(NAME_2='United.States',NAME_3='Region4',NAME_4='North.Carolina') and mod_level=4 specifies North Carolina. To achieve the same result, use all abbreviations mod_name=c(NAME_2='US',NAME_3='R4',NAME_4='NC') or a mix of names and abbreviations mod_name=c(NAME_2='US',NAME_3='Region4',NAME_4='NC'). Unlike get.cdc.data(), get_disease_data() expects ISO2 country abbreviations.
  #' @param years An integer vector - all years to be included.  If a value of "all" is passed, the function returns all available dates.
  #' @param disease A character string - for example 'flu' or 'dengue'
  #' @param sql_data_source An integer - designates which mydata source is to be used.
  #' @param daily_clim Logic value determining if daily climate data is retrieved.  (Note: climate data averaged to the incidence-cadence is always returned.)
  #' @param mod_sum Logic value indicating if model data will be aggregated from fit data.  If so, get_disease_data will not attempt to return model-level data.
  #' @return out  A list with a lookup table 'attr' and then climate data organized by cadence, spatial patch, and then data_source.
  #' @examples
  #' require(DICE)
  #' out = get_disease_data(mod_level = 2, fit_level = 4, mod_name = c(NAME_2='BR'), years="all", disease="dengue", sql_data_source=1)

  # open mydatabase connection
  myDB = OpenCon(db_opts)
  # retrieve data_sources table
  data_sources = dbReadTable(myDB, name="data_sources")

  available_diseases = unique(data_sources$disease)

  clim_by_disease = dbReadTable(myDB, name="clim_by_disease")
  clim_names = clim_by_disease$clim_names[clim_by_disease$disease==disease]
  clim_names = strsplit(clim_names, ";")[[1]]

  if (!any(disease==available_diseases)) {
    stop("The disease entered c(", disease, "), does not match available diseases c(", paste(available_diseases, collapse=", "), ")." )
  }

  # download the disease-specific lookup table
  # dis_lut = dbReadTable(myDB, name=paste0(disease,"_lut"))
  dis_lut = dbGetQuery(myDB, statement=paste0("SELECT * FROM ", disease, "_lut ORDER BY master_key;"))
  dis_lut = StandardizeLUT(dis_lut=dis_lut, DICE_db=db_opts$DICE_db)

  # construct the WHERE string(s) for MySQL query
  where = character(3)
  # start with conditional for master_key
  if (length(mod_name)==1 && mod_name=="all") {
    where[1] = ""
  } else {
    # determine model- and fit-level master_keys
    if (mod_level>=2) {
      mod_ind = (toupper(mod_name["NAME_2"]) == dis_lut$ABBV_2 | mod_name["NAME_2"] == dis_lut$NAME_2) & mod_level == dis_lut$level
      fit_ind = (toupper(mod_name["NAME_2"]) == dis_lut$ABBV_2 | mod_name["NAME_2"] == dis_lut$NAME_2) & fit_level == dis_lut$level
      if (mod_level > 2) {
        for (lev in 3:mod_level) {
          mod_ind = mod_ind & (mod_name[paste0("NAME_", lev)] == dis_lut[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_", lev)]) == dis_lut[, paste0("ABBV_", lev)])
          fit_ind = fit_ind & (mod_name[paste0("NAME_", lev)] == dis_lut[, paste0("NAME_", lev)] | toupper(mod_name[paste0("NAME_", lev)]) == dis_lut[, paste0("ABBV_", lev)])
        }
      }
    } else if (mod_level==1) {
      # Continent level $model will be constructed as a sum of $fit regions
      if (fit_level==1) {
        # sum all available countries
        fit_ind = (toupper(mod_name["NAME_1"]) == dis_lut$ABBV_1 | mod_name["NAME_1"] == dis_lut$NAME_1) & 2 == dis_lut$level
      } else {
        # determine indices of fit_level indices
        fit_ind = (toupper(mod_name["NAME_1"]) == dis_lut$ABBV_1 | mod_name["NAME_1"] == dis_lut$NAME_1) & fit_level == dis_lut$level
      }
      mod_ind = rep(FALSE, length(fit_ind))
    } else {
      # Global/Continent level mydata
    }

    if (mod_sum) {
      mod_ind = rep(FALSE, length(fit_ind))
    }
    master_keys = dis_lut$master_key[mod_ind | fit_ind]
    where[1] = paste0("master_key IN('",paste(master_keys, collapse="','"),"')")
  }

  # create date condition for MySQL query
  if (length(years)==1 && years=="all") {
    where[2] = ""
  } else {
    min_date = as.Date(paste0(min(years)-1, "-12-28"))
    max_date = as.Date(paste0(max(years), "-12-31"))
    where[2] = paste0("date BETWEEN '", min_date, "' AND '", max_date, "'")
  }

  # create source_key condition for MySQL query
  if (length(sql_data_source)==0 || sql_data_source==26) {
    where[3] = ""
  } else {
    where[3] = paste0("source_key IN(",paste(sql_data_source, collapse=","),")")
  }

  if (any(where!="")) {
    where = where[where!=""]
    where_string = paste0("WHERE ", paste(where, collapse=" AND "), " ")
  } else {
    where_string = ""
  }

  cat("Downloading incidence and climate data......")
  # REMOVE_PUBLIC_START
  if (length(sql_data_source)>0 && sql_data_source==26) {
    qDB = OpenCon(db_opts)
    query_string = paste0("SELECT * FROM quidel_data ", where_string, "ORDER BY master_key, date")
    inc_data = dbGetQuery(qDB, statement=query_string)
    dbDisconnect(qDB)
  } else {
    # REMOVE_PUBLIC_END
    query_string = paste0("SELECT * FROM ",disease,"_data ", where_string, "ORDER BY master_key, source_key, date")
    inc_data = dbGetQuery(myDB, statement=query_string)
    # REMOVE_PUBLIC_START
  }
  # REMOVE_PUBLIC_END

  # check if any data was returned
  if (length(inc_data$master_key)==0) {
    stop("No data in MySQL mydatabase matching the mod_name/level/dates/disease/mydata source combination entered into get_disease_data().")
  }
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

  if (daily_clim | cad_clim) {
    # retrieve clim_by_disease table
    clim_by_disease = dbReadTable(conn=myDB, name="clim_by_disease")
    clim_cols = strsplit(clim_by_disease$clim_names[clim_by_disease$disease==disease], ";")[[1]]
  }

  if (daily_clim) {
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
    query_string = paste0("SELECT * FROM noaa_daily WHERE date BETWEEN '", min_date, "' AND '", max_date, "' AND master_key IN('", paste(master_keys, collapse="','"), "') ORDER BY master_key, date")
    noaa_clim = dbGetQuery(myDB, statement=query_string)
    cat("Complete\n")
  }

  if (cad_clim) {
    cat("Downloading climate data at appropriate cadence. \n")
    noaa_cad = list()
    # determine which master_keys are needed
    clim_idents = dis_lut$clim_ident[dis_lut$master_key %in% master_keys]
    clim_keys   = dis_lut$master_key[dis_lut$identifier %in% clim_idents]
    for (cad in cads) {
      cad_name = tolower(CadenceNum2Name(cadences=cad))

      query_string = paste0("SELECT * FROM noaa_", cad_name, " WHERE master_key IN('", paste(clim_keys, collapse="','"), "') ORDER BY master_key, date")
      noaa_cad[[paste0("noaa_", cad_name)]] = dbGetQuery(myDB, statement=query_string)
    }
  }

  # convert data into useable and user-friendly format
  out = list()
  out$lut = dis_lut[dis_lut$master_key %in% master_keys, ]

  countries = out$lut$NAME_2[out$lut$level==2]
  cadences = unique(inc_data$cadence)
  cad_names = CadenceNum2Name(cadences)

  cat("Formatting data......")
  for (ii in 1:length(cadences)) {
    out[[cad_names[ii]]] = list()
    cad_list_ind = which(names(out)==cad_names[ii])
    cad_data_ind = inc_data$cadence==cadences[ii]
    sources_cadii = data_sources$source_key[data_sources$cadence==cadences[ii]]

    # establish which columns of inc_data are needed
    data_names = names(inc_data)[substr(names(inc_data), start=1, stop=4)=="data"]
    n_data_names = length(data_names)
    if (cadences[ii] %in% c(0,1)) {
      data_cols = c("date", "year", "day", data_names)
      max_day_fun <- function(max_date) {
        return(max_date)
      }
    } else if (cadences[ii]==2) {
      data_cols = c("date", "year", "week", data_names)
      max_day_fun <- function(max_date) {
        return(max_date+6)
      }
    } else if (cadences[ii]==3) {
      data_cols = c("date", "year", "month", data_names)
      max_day_fun <- function(max_date) {
        month(max_date) = month(max_date)+1
        return(max_date-1)
      }
    } else if (cadences[ii]==4) {
      data_cols = c("date", "year", data_names)
      max_day_fun <- function(max_date) {
        year(max_date) = year(max_date)+1
        return(max_date-1)
      }
    }

    # determine cadence patches
    master_keys = unique(inc_data$master_key[inc_data$cadence==cadences[ii]])
    for (jj in 1:length(master_keys)) {
      patch_ind = which(out$lut$master_key==master_keys[jj])
      dis_patch_ind = which(dis_lut$master_key==master_keys[jj])
      ident = out$lut$identifier[patch_ind]
      out[[cad_list_ind]][[ident]] = list()

      # loop through data sources with the same cadence
      sources = unique(inc_data$source_key)
      sources = sources[sources %in% sources_cadii]
      for (source_key in sources) {
        data_ind = inc_data$master_key==master_keys[jj] & cad_data_ind & inc_data$source_key==source_key
        source_abbv = data_sources$source_abbv[data_sources$source_key==source_key]

        # pull incidence data
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
        if (source_key==sources[1]) {
          min_date = out[[cad_list_ind]][[ident]][[source_abbv]]$date[1]
          max_date = out[[cad_list_ind]][[ident]][[source_abbv]]$date[length(out[[cad_list_ind]][[ident]][[source_abbv]]$date)]
          max_date = max_day_fun(max_date)
        } else {
          min_date = min(min_date, out[[cad_list_ind]][[ident]][[source_abbv]]$date[1])
          max_date = max(max_date, max_day_fun(out[[cad_list_ind]][[ident]][[source_abbv]]$date[length(out[[cad_list_ind]][[ident]][[source_abbv]]$date)]))
        }
      }

      if (daily_clim) {
        # pull daily climate data
        if (dis_lut$identifier[dis_patch_ind]==dis_lut$clim_ident[dis_patch_ind]) {
          clim_key = master_keys[jj]
        } else {
          # patches that are smaller than the smallest GADM level are assigned the climate data of the smallest encompassing GADM patch
          clim_key = dis_lut$master_key[dis_lut$identifier==dis_lut$clim_ident[dis_patch_ind]]
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

      if (cad_clim) {
        # pull climate data with proper cadence
        if (dis_lut$identifier[dis_patch_ind]==dis_lut$clim_ident[dis_patch_ind]) {
          clim_key = master_keys[jj]
        } else {
          # patches that are smaller than the smallest GADM level are assigned the climate data of the smallest encompassing GADM patch
          clim_key = dis_lut$master_key[dis_lut$identifier==dis_lut$clim_ident[dis_patch_ind]]
        }
        out[[cad_list_ind]][[ident]][[paste0("noaa_", tolower(cad_names[ii]))]] = noaa_cad[[paste0("noaa_", tolower(cad_names[ii]))]][noaa_cad[[paste0("noaa_", tolower(cad_names[ii]))]]$master_key==clim_key, c("date", clim_cols)]
        # Reformat dates to Date class
        out[[cad_list_ind]][[ident]][[paste0("noaa_", tolower(cad_names[ii]))]]$date = as.Date(out[[cad_list_ind]][[ident]][[paste0("noaa_", tolower(cad_names[ii]))]]$date)
      }

    }
  }
  cat("Complete\n")
  
  dbDisconnect(myDB)
  return(out)
}


OpenCon <- function(db_opts=list(DICE_db="predsci")) {
  #' Open a Connection to MySQL Database.
  #'
  #' This function uses guest credentials to open a read-only connection to the MySQL database.
  #' @param sql_db Generally, data has been consolidated to 'PredSci'.
  #' @return myDB  An S4 object that inherits from DBIConnection.
  #' @examples
  #' require(DICE)
  #' myDB = OpenCon()
  #' myDB
  #' dbDisconnect(myDB)

  sql_db = db_opts$DICE_db
  
  if (tolower(sql_db)=="predsci") {
    drv  = MySQL()
    user = "epi_guest"
    password="UY5GE2kfUa"
    dbname='epi_data'
    port = 3306
    if (Sys.info()["nodename"]=="Q") {
      host="shadow"
    } else {
      host="shadow.predsci.com"
    }
  } else if (tolower(sql_db)=="predsci2") {
    drv  = MySQL()
    user = "epi_guest"
    password="UY5GE2kfUa"
    dbname='epi_data'
    port = 3306
    if (Sys.info()["nodename"]=="Q") {
      host="shadow2"
    } else {
      host="shadow2.predsci.com"
    }
  } else if (tolower(sql_db)=="quidel") {
    drv = MySQL()
    if (!(all(names(db_opts) %in% c("quidel_user", "quidel_pwd")))) {
      stop("Accessing Quidel data requires username and password.")
    }
    user = db_opts$quidel_user
    password = db_opts$quidel_pwd
    dbname = "quidel_data"
    port = 3306
    if (Sys.info()["nodename"]=="Q") {
      host="shadow"
    } else {
      host="shadow.predsci.com"
    }
  }

  for (ii in 1:10) {
    # attempt to establish the database connection
    myDB=try(dbConnect(drv=drv, user=user, password=password, dbname=dbname, host=host, port=port), silent=TRUE)
    if (!is(myDB, "try-error")) {
      # if successful, continue
      break
    } else {
      cat("Failed to connect to ", sql_db, " database on attempt ", ii, ".\n", sep="")
      # if not, wait 0.1s and try again
      Sys.sleep(0.1)
    }
  }

  # if loop exited without making a connection, throw error
  if (is(myDB, "try-error")) {
    stop("After ", ii, " attempts, could not establish a connection to ", sql_db, " database.\n")
  }

  return(myDB)
}


CadenceNum2Name <- function(cadences=NULL) {
  #' Number to String Cadence Conversion
  #'
  #' Convert number-based cadence representations to descriptive strings.
  #'
  #' @param cadences An integer vector of cadences.
  #' @return cadence_names  A a vector of descriptive strings.
  #' @examples
  #' require(DICE)
  #' cadence_names = CadenceNum2Name(cadences=0:4)
  cadence_names = character(length(cadences))
  cadence_names[cadences==0] = "By Case"
  cadence_names[cadences==1] = "Daily"
  cadence_names[cadences==2] = "Weekly"
  cadence_names[cadences==3] = "Monthly"
  cadence_names[cadences==4] = "Annual"
  return(cadence_names)
}


Date2CadLabel <- function(query_result=NULL) {
  #' Add Cadence to Data
  #'
  #' Take the date/cadence output from a MySQL query and add appropriate 'year', 'month', 'week', or 'day' columns.
  #'
  #' @param query_result A data frame resulting from a MySQL incidence query.  A cadence column should be added after the query, but before calling this function.
  #' @return query_result  The same dataframe, but with added 'year', 'month', 'week', or 'day' columns appropriate to the incidence cadence.
  #' @examples
  #' require(DICE)
  #' myDB = OpenCon()
  #' data_sources = dbReadTable(myDB, name="data_sources")
  #' dbDisconnect(myDB)
  #' query_result = get_disease_data()
  #' query_result$date = as.Date(query_result$date)
  #' # add auxilary columns to query_result
  #' source_ind = match(query_result$source_key, data_sources$source_key)
  #' query_result$cadence = data_sources$cadence[source_ind]
  #' # add year, month, week, columns
  #' query_result = Date2CadLabel(query_result=query_result)

  yearly_ind = query_result$cadence==4
  if (any(yearly_ind)) {
    query_result$year[yearly_ind] = year(query_result$date[yearly_ind])
  }

  monthly_ind = query_result$cadence==3
  if (any(monthly_ind)) {
    query_result$year[monthly_ind] = year(query_result$date[monthly_ind])
    query_result$month[monthly_ind] = month(query_result$date[monthly_ind])
  }

  weekly_ind = query_result$cadence==2
  if (any(weekly_ind)) {
    # find all unique dates
    dates = unique(query_result$date[weekly_ind])
    # map to full list
    dates_ind = match(query_result$date[weekly_ind], dates)

    # convert dates to MMWR week and year
    weeks = integer(length(dates))
    years = weeks
    for (ii in 1:length(dates)) {
      CDCweek = Date2CDCweek(dates[ii])
      weeks[ii] = CDCweek$CDCweek
      years[ii] = CDCweek$CDCyear
    }

    query_result$week[weekly_ind] = weeks[dates_ind]
    query_result$year[weekly_ind] = years[dates_ind]
  }

  daily_ind = query_result$cadence==1
  if (any(daily_ind)) {
    query_result$year[daily_ind] = year(query_result$date[daily_ind])
    query_result$day[daily_ind] = yday(query_result$date[daily_ind])
  }
  return(query_result)
}


BuildDateVecs <- function(cadence=cadence, start.date=start.date, end.date=end.date) {
  #' Build A Date Vector
  #'
  #' Given a cadence, start and end date build a vector of dates
  #' @param cadence Numeric (1 = days, 2 = weeks, 3 = months and 4 = years)
  #' @param start.date Start Date - Type Date
  #' @param end.date   End Date   - Type Date
  #' @return mydata A list with two/three elements: date, year and day/week/month
  #' @examples
  #' mydata <- BuildDateVecs(cadence = 2, start.date = "01-01-2015", end.date = "12-31-2015")
  mydata = list()
  if (cadence==1) { # days
    mydata$date = seq.Date(from=start.date, to=end.date, by=1)
    mydata$year = year(mydata$date)
    mydata$day  = yday(mydata$date)
  } else if (cadence==2) { # weeks
    start.week = Date2CDCweek(start.date)
    end.week   = Date2CDCweek(end.date)
    start.date = CDCweek2date(CDCweek=start.week$CDCweek, year=start.week$CDCyear)
    end.date   = CDCweek2date(CDCweek=end.week$CDCweek, year=end.week$CDCyear)
    mydata$date = seq.Date(from=start.date, to=end.date, by=7)
    mydata$year = integer(length(mydata$date))
    mydata$week = integer(length(mydata$date))
    for (ii in 1:length(mydata$date)) {
      CDCweek = Date2CDCweek(date=mydata$date[ii])
      mydata$year[ii] = CDCweek$CDCyear
      mydata$week[ii] = CDCweek$CDCweek
    }
  } else if (cadence==3) { # months
    mday(start.date) = 1
    mday(end.date)   = 1
    mydata$date      = seq.Date(from=start.date, to=end.date, by="month")
    mydata$year      = year(mydata$date)
    mydata$month     = month(mydata$date)
  } else if (cadence==4) { # years
    start.year = year(start.date)
    end.year   = year(end.date)
    years      = start.year:end.year
    mydata$date = as.Date(paste0(years, "-01-01"))
    mydata$year = years
  }
  return(mydata)
}



PadNAs <- function(ident_data=NULL, cadence=NULL, end.date=NULL, nmod=NULL) {
  #' Pad Time Series wit NA's
  #'
  #' The function \code{PadNAs} takes a single identifier time-series dataframe from a 'data' \pkg{DICE} structure and add dates through 'end.date'.
  #' This function is generally called by get.mysql() to add forecast dates beyond the last climate data point.
  #' @param ident_data A dataframe from 'mydata' for a single identifier's time-series mydata.
  #' @param cadence An integer that specifies incidence cadence.
  #' @param end.date A Date class variable specifying the last date to be forecast.
  #' @return ident_data  The same mydataframe, but with added cadence periods until 'end.date'.
  #' @examples
  #' require(DICE)
  #'
  nmod = length(ident_data$dates)
  if (cadence==2) {
    CDCweek = Date2CDCweek(end.date)
    # just in case 'end.date' falls in the middle of a week
    lastweek_date = CDCweek2date(CDCweek=CDCweek$CDCweek, year=CDCweek$CDCyear)
    first_padweek_date = ident_data$dates[nmod]+7
    if (first_padweek_date<=lastweek_date) {
      addDates = seq.Date(from=first_padweek_date, to=lastweek_date, by=7)
      new_rows = data.frame(dates=addDates, years=NA, weeks=NA)
      for (kk in 1:length(addDates)) {
        CDCweek = Date2CDCweek(date=addDates[kk])
        new_rows$years[kk] = CDCweek$CDCyear
        new_rows$weeks[kk] = CDCweek$CDCweek
      }
    } else {
      new_rows = data.frame(dates=NULL, years=NULL, months=NULL)
    }

  } else if (cadence==3) {
    lastmonth_date = end.date
    mday(lastmonth_date) = 1
    new_month = ident_data$dates[nmod]
    if (new_month < lastmonth_date) {
      month(new_month) = month(new_month) + 1
      add_month = new_month
      while(new_month < lastmonth_date) {
        month(new_month) = month(new_month) + 1
        add_month = c(add_month, new_month)
      }
      new_rows = data.frame(dates=add_month, years=year(add_month), months=month(add_month))
    } else {
      new_rows = data.frame(dates=NULL, years=NULL, months=NULL)
    }

  } else if (cadence==4) {
    lastyear = year(end.date)
    new_year = (year(ident_data$dates[nmod])+1):lastyear
    new_rows = data.frame(dates=as.Date(paste0(new_year,"-01-01")), years=new_year)
  }
  # fill non-date columns with NA
  if (length(new_rows$dates)>0) {
    model_names = names(ident_data)
    for (kk in 4:length(model_names)) {
      new_rows[[model_names[kk]]] = NA
    }
    ident_data = rbind(ident_data, new_rows)
  }

  return(ident_data)
}


AvgClim <- function(ident_data=NULL, clim_cols=NULL, all_years_ident=NULL, cadence=NULL) {
  #' Climate Data Processing.
  #'
  #' This function is generally called by \code{\link{get.mysql()}} to fill future climate dates with appropriate historic averages.
  #' After padding incidence with NAs, check climate data for NAs and replace with historic average.  Historic average to be calculated from all years available.
  #'
  #' @param ident_data A mydataframe from 'mydata' for a single identifier's time-series mydata.
  #' @param clim_cols The names of climate columns in ident_data (should match with column names in all_years_ident).
  #' @param all_years_ident The output of get_disease_data containing all available years of mydata.
  #' @param cadence An integer that specifies incidence cadence.
  #' @return ident_data  The same mydataframe, but with added historic averages for future climate dates.
  #' @examples
  #' require(DICE)
  #'
  for (clim in clim_cols) {
    na_ind = is.na(ident_data[[clim]])
    if (any(na_ind)) {
      # average past periods to estimate future
      if (cadence==2) {
        for(ii in which(na_ind)) {
          week = ident_data$weeks[ii]
          # index week in all_years mydata
          week_ind = all_years_ident$week==week
          temp_mean = mean(all_years_ident[[clim]][week_ind], na.rm=TRUE)
          if (!is.nan(temp_mean)) {
            ident_data[[clim]][ii] = temp_mean
          } else {
            ident_data[[clim]][ii] = NA
          }
        }
      } else if (cadence==3) {
        for(ii in which(na_ind)) {
          month = ident_data$months[ii]
          # index month in all_years mydata
          month_ind = all_years_ident$month==month
          temp_mean = mean(all_years_ident[[clim]][month_ind], na.rm=TRUE)
          if (!is.nan(temp_mean)) {
            ident_data[[clim]][ii] = temp_mean
          } else {
            ident_data[[clim]][ii] = NA
          }
        }
      } else if (cadence==1) {
        for(ii in which(na_ind)) {
          day = ident_data$days[ii]
          # index month in all_years mydata
          day_ind = all_years_ident$day==day
          temp_mean = mean(all_years_ident[[clim]][day_ind], na.rm=TRUE)
          if (!is.nan(temp_mean)) {
            ident_data[[clim]][ii] = temp_mean
          } else {
            ident_data[[clim]][ii] = NA
          }
        }
      } else if (cadence==4) {
        # do nothing, running mechanistic forecasts on a yearly cadence doesn't make sense
      }

    }
  }
  return(ident_data)
}



mydata_New2Old <- function(mydata_new=NULL) {
  #' Data Format Conversion
  #'
  #' Convert the dataframe format from the new one to the old one.
  #' Can only be used by expert Users
  #'
  mydata_old = mydata_new

  # convert 'model' structure
  mod_ident = mydata_new$model$attr$identifier
  temp_ident = mydata_old$model[[mod_ident]]
  var_names = names(mydata_new$model[[mod_ident]])
  date_names= var_names[1:length(var_names) <= which(var_names=="ndays")]
  var_names = var_names[1:length(var_names) > which(var_names=="ndays")]
  var_names = var_names[var_names!="pop"]
  mydata_old$model = c(mydata_old$model, mydata_old$model[[mod_ident]][, var_names])
  # remove identifier-time-series structure
  mydata_old$model[[mod_ident]] = NULL

  # convert 'fit' structure
  ii = 1
  fit_ident = mydata_new$fit$attr$identifier[ii]
  var_names = names(mydata_new$fit[[fit_ident]])
  date_names= var_names[1:length(var_names) <= which(var_names=="ndays")]
  var_names = var_names[1:length(var_names) > which(var_names=="ndays")]
  var_names = var_names[var_names!="pop"]
  for (vname in var_names) {
    mydata_old$fit[[vname]] = data.frame(mydata_new$fit[[fit_ident]][, vname])
    names(mydata_old$fit[[vname]]) = fit_ident
  }
  # remove 'new' structure
  mydata_old$fit[[fit_ident]] = NULL

  # process the remaining regions
  if (mydata_old$fit$nregions>1) {
    for (ii in 2:mydata_old$fit$nregions) {
      fit_ident = mydata_new$fit$attr$identifier[ii]
      for (vname in var_names) {
        if (vname %in% names(mydata_new$fit[[fit_ident]])) {
          mydata_old$fit[[vname]][[fit_ident]] = mydata_new$fit[[fit_ident]][, vname]
        } else {
          mydata_old$fit[[vname]][[fit_ident]] = as.numeric(NA)
        }
        # names(mydata_old$fit[[vname]])[ii] = fit_ident
      }
      mydata_old$fit[[fit_ident]] = NULL
    }
  }
  return(mydata_old)
}


pluralize <- function(ident_dat=NULL) {
  #' Pluralize names
  #'
  #' Pluralize names such as: date, day, week, month and year
  #'
  old_names = names(ident_dat)
  new_names = old_names
  new_names[new_names=="date"]  = "dates"
  new_names[new_names=="day"]   = "days"
  new_names[new_names=="week"]  = "weeks"
  new_names[new_names=="month"] = "months"
  new_names[new_names=="year"]  = "years"

  names(ident_dat) = new_names
  return(ident_dat)
}


add_ndays <- function(ident_dat=NULL, cadence=NULL) {
  #'
  #' Add Days
  #'
  ncols = length(ident_dat)
  if (cadence==1) {
    ndays = rep(1, length(ident_dat$dates))
    if (length(ident_dat)>3) {
      ident_dat = cbind(ident_dat[,1:3], data.frame(ndays=ndays), ident_dat[,4:ncols, drop=FALSE])
    }
  } else if (cadence==2) {
    ndays = rep(7, length(ident_dat$dates))
    if (length(ident_dat)>3) {
      ident_dat = cbind(ident_dat[,1:3], data.frame(ndays=ndays), ident_dat[,4:ncols, drop=FALSE])
    } else {
      ident_dat = cbind(ident_dat[,1:3], data.frame(ndays=ndays))
    }

  } else if (cadence==3) {
    ndays = integer(length(ident_dat$dates))
    for (ii in 1:length(ident_dat$dates)) {
      ndays[ii] = days_in_month(ident_dat$dates[ii])
    }
    if (length(ident_dat)>3) {
      ident_dat = cbind(ident_dat[,1:3], data.frame(ndays=ndays), ident_dat[,4:ncols, drop=FALSE])
    } else {
      ident_dat = cbind(ident_dat[,1:3], data.frame(ndays=ndays))
    }

  } else if (cadence==4) {
    ndays = integer(length(ident_dat$dates))
    for (ii in 1:length(ident_dat$dates)) {
      next_year = ident_dat$dates[ii]
      year(next_year) = year(next_year) + 1
      ndays[ii] = as.integer(next_year-ident_dat$dates[ii])
    }
    if (length(ident_dat)>2) {
      ident_dat = cbind(ident_dat[,1:2], data.frame(ndays=ndays), ident_dat[,3:ncols, drop=FALSE])
    } else {
      ident_dat = cbind(ident_dat[,1:2], data.frame(ndays=ndays))
    }

  }
  # ident_dat$ndays = ndays
  return(ident_dat)
}


add_periods_cols <- function(ident_dat=NULL, cadence=NULL) {
  #'
  #' Add a year column as well as the appropriate day/week/month column
  #'
  ncols = length(ident_dat)
  if (nrow(ident_dat) > 0) {
    if (cadence==1) {
      year = year(ident_dat$date)
      day  = yday(ident_dat$date)

      ident_dat = cbind(date=ident_dat$date, data.frame(year=year, day=day), ident_dat[2:ncols])
    } else if (cadence==2) {
      year = epiyear(ident_dat$date)
      week = epiweek(ident_dat$date)

      ident_dat = cbind(date=ident_dat$date, data.frame(year=year, week=week), ident_dat[2:ncols])

    } else if (cadence==3) {
      year  = year(ident_dat$date)
      month = month(ident_dat$date)

      ident_dat = cbind(date=ident_dat$date, data.frame(year=year, month=month), ident_dat[2:ncols])

    } else if (cadence==4) {
      year  = year(ident_dat$date)

      ident_dat = cbind(date=ident_dat$date, data.frame(year=year), ident_dat[2:ncols])

    }
  } else {
    # add columns to empty data.frame
    ident_names = names(ident_dat)
    if (cadence==1) {
      ident_dat = data.frame(matrix(ncol=length(ident_names)+2, nrow=0))
      names(ident_dat) = c("date", "year", "day", ident_names[2:length(ident_names)])
    } else if (cadence==2) {
      ident_dat = data.frame(matrix(ncol=length(ident_names)+2, nrow=0))
      names(ident_dat) = c("date", "year", "week", ident_names[2:length(ident_names)])
    } else if (cadence==3) {
      ident_dat = data.frame(matrix(ncol=length(ident_names)+2, nrow=0))
      names(ident_dat) = c("date", "year", "month", ident_names[2:length(ident_names)])
    } else if (cadence==4) {
      ident_dat = data.frame(matrix(ncol=length(ident_names)+1, nrow=0))
      names(ident_dat) = c("date", "year", ident_names[2:length(ident_names)])
    }

  }


  return(ident_dat)
}



get_season_dates <- function(myDB=NULL, season=NULL, continent=NULL, country=NULL, cadence=NULL, disease=NULL) {

  if (is.null(country)) {
    # look for continent default dates
    SE_dates = dbGetQuery(conn=myDB, statement=paste0("SELECT * FROM season_se_dates WHERE disease='", disease, "' AND cadence=", cadence," AND abbv_1='", continent, "' AND abbv_2='default'"))
  } else {
    # Retrieve all seasons, both default and country-specific
    query_dates = dbGetQuery(conn=myDB, statement=paste0("SELECT * FROM season_se_dates WHERE disease='", disease, "' AND cadence=", cadence," AND abbv_1='", continent, "' AND abbv_2 IN ('default','", country,"')"))
    # where country-specific dates exist, remove defaults
    country_ind = query_dates$abbv_2!="default"
    SE_dates    = query_dates[country_ind, ]
    # add defaults for other years
    default_use = query_dates$abbv_2=="default" & !(query_dates$season %in% query_dates$season[country_ind])
    SE_dates = rbind(SE_dates, query_dates[default_use, ])
    SE_dates = SE_dates[order(SE_dates$season), ]
  }
  return(SE_dates)
}



SummaryByCountry <- function(diseases="all", countries="all", DICE_db="predsci") {
  #' Disease Summary By Country
  #'
  #' Generate a summary of each disease/country/level available in the database.  Attributes returned include number of regions per level, data cadence, and min/max dates by country/level/cadence.
  #' @param diseases Can be a vector of diseases c("dengue", "flu") or "all"
  #' @param countries Can be a vector of countries (or ISO-2 abbreviations) c("Brazil", "MX", "Thailand") or "all".  Not case sensitive.
  #' @return A list of diseases.  Each disease entry contains a sub-list of countries.  Each country contains total regions, regions-by-level, and min/max data-dates by cadence/level.

  # open database connection
  myDB = OpenCon(list(DICE_db=DICE_db))
  # retrieve data_sources table
  data_sources = dbReadTable(myDB, name="data_sources")

  available_diseases = unique(data_sources$disease)

  if (tolower(diseases)=="all") {
    diseases = available_diseases
  } else if (!all(diseases %in% available_diseases)) {
    stop("The diseases entered c(", paste(diseases, collapse=", "), "), do not match available diseases c(", paste(available_diseases, collapse=", "), ")." )
  }

  summary = list()
  for(disease in diseases) {
    summary[[disease]] = list()
    # retrieve look-up table
    # dis_lut = dbReadTable(myDB, name=paste0(disease, "_lut"))
    dis_lut = dbGetQuery(myDB, statement=paste0("SELECT * FROM ", disease, "_lut ORDER BY master_key;"))
    dis_lut = StandardizeLUT(dis_lut=dis_lut, DICE_db=DICE_db)
    # retrieve min/max dates for each patch
    date_ranges = dbGetQuery(myDB, statement=paste0("SELECT master_key, source_key, MAX(date), MIN(date) FROM ",disease,"_data GROUP BY master_key, source_key"))
    names(date_ranges) = c("master_key", "source_key", "max_date", "min_date")
    date_ranges$min_date = as.Date(date_ranges$min_date, tz = "UTC")
    date_ranges$max_date = as.Date(date_ranges$max_date, tz = "UTC")

    # augment with info from lookup table and source table
    lut_ind = match(date_ranges$master_key, dis_lut$master_key)
    date_ranges$country = dis_lut$NAME_2[lut_ind]
    date_ranges$ABBV_2  = dis_lut$ABBV_2[lut_ind]
    date_ranges$level   = dis_lut$level[lut_ind]
    source_ind = match(date_ranges$source_key, data_sources$source_key)
    date_ranges$cadence = data_sources$cadence[source_ind]

    data_countries = tolower(unique(date_ranges$country))
    data_abbv2s    = tolower(unique(date_ranges$ABBV_2))
    if (length(countries)==1 && tolower(countries)=="all") {
      dis_countries = data_countries
    } else {
      dis_countries = data_countries[data_countries %in% tolower(countries) | data_abbv2s %in% tolower(countries)]
    }


    for (country in dis_countries) {
      summary[[disease]][[country]] = list()
      country_ind = tolower(date_ranges$country)==country
      country_patches = unique(date_ranges$master_key[country_ind])
      # country_sources = unique(date_ranges$source_key[date_ranges$master_key %in% country_patches])
      cadences = sort(unique(date_ranges$cadence[country_ind]))
      cadence_names = CadenceNum2Name(cadences)
      summary[[disease]][[country]]$tot_regions = length(country_patches)

      country_levels = unique(date_ranges$level[country_ind])
      num_regions = country_levels
      for (ii in 1:length(country_levels)) {
        lev = country_levels[ii]
        lev_patches = unique(date_ranges$master_key[date_ranges$level==lev & country_ind])
        num_regions[ii] = length(lev_patches)
      }
      summary[[disease]][[country]]$RegionsByLevel = data.frame(level=country_levels, num_regions=num_regions)

      date_ranges_ind = date_ranges$master_key %in% country_patches
      country_sources = unique(date_ranges$source_key[date_ranges_ind])
      for (ii in 1:length(cadences)) {
        cadence_ind = date_ranges$cadence==cadences[ii] & country_ind
        levels = unique(date_ranges$level[cadence_ind])
        summary[[disease]][[country]][[cadence_names[ii]]] = data.frame(level=levels, num_regions=0, min_date=rep(as.Date("2000-01-01"), length(levels)), max_date=rep(as.Date("2000-01-01"), length(levels)))
        for (jj in 1:length(levels)) {
          lev_ind = cadence_ind & date_ranges$level==levels[jj]
          num_regions = length(unique(date_ranges$master_key[lev_ind]))
          min_date = min(date_ranges$min_date[lev_ind])
          max_date = max(date_ranges$max_date[lev_ind])
          summary[[disease]][[country]][[cadence_names[ii]]][jj, 2:4] = list(num_regions, min_date, max_date)
        }
      }
    }
  }
  dbDisconnect(myDB)
  return(summary)
}


SummaryBySource <- function(diseases="all", countries="all", DICE_db="predsci") {
  #' Disease Summary By Source
  #'
  #' Add text here
  #'
  # open database connection
  myDB = OpenCon(list(DICE_db=DICE_db))
  # retrieve data_sources table
  data_sources = dbReadTable(myDB, name="data_sources")

  available_diseases = unique(data_sources$disease)

  if (tolower(diseases)=="all") {
    diseases = available_diseases
  } else if (!all(diseases %in% available_diseases)) {
    stop("The diseases entered c(", paste(diseases, collapse=", "), "), do not match available diseases c(", paste(available_diseases, collapse=", "), ")." )
  }

  summary = list()
  for(disease in diseases) {
    summary[[disease]] = list()
    # retrieve look-up table
    # dis_lut = dbReadTable(myDB, name=paste0(disease, "_lut"))
    dis_lut = dbGetQuery(myDB, statement=paste0("SELECT * FROM ", disease, "_lut ORDER BY master_key;"))
    dis_lut = StandardizeLUT(dis_lut=dis_lut, DICE_db=DICE_db)
    # retrieve min/max dates for each patch
    date_ranges = dbGetQuery(myDB, statement=paste0("SELECT master_key, source_key, MAX(date), MIN(date) FROM ",disease,"_data GROUP BY master_key, source_key"))
    names(date_ranges) = c("master_key", "source_key", "max_date", "min_date")
    date_ranges$min_date = as.Date(date_ranges$min_date, tz = "UTC")
    date_ranges$max_date = as.Date(date_ranges$max_date, tz = "UTC")

    # augment with info from lookup table and source table
    lut_ind = match(date_ranges$master_key, dis_lut$master_key)
    date_ranges$country = dis_lut$NAME_2[lut_ind]
    date_ranges$ABBV_2  = dis_lut$ABBV_2[lut_ind]
    date_ranges$level   = dis_lut$level[lut_ind]
    source_ind = match(date_ranges$source_key, data_sources$source_key)
    date_ranges$cadence = data_sources$cadence[source_ind]

    data_countries = unique(date_ranges$country)
    data_abbv2s    = tolower(unique(date_ranges$ABBV_2))
    if (length(countries)==1 && tolower(countries)=="all") {
      dis_countries = data_countries
    } else {
      dis_countries = data_countries[tolower(data_countries) %in% tolower(countries) | data_abbv2s %in% tolower(countries)]
    }

    # ignore any NA entries
    dis_countries = dis_countries[!is.na(dis_countries)]
    for (country in dis_countries) {
      summary[[disease]][[country]] = list()
      country_ind = date_ranges$country==country

      # determine root location of each data source
      source_info = data.frame(source_keys=unique(date_ranges$source_key[country_ind]), level=NA, loc=NA, master_key=NA, cadence=NA)
      for (ii in 1:nrow(source_info)) {
        temp_ind = source_info$source_keys[ii]==date_ranges$source_key & country_ind
        # find common parent of all date_ranges$master_key[temp_ind]
        loc_lut_ind = dis_lut$master_key %in% date_ranges$master_key[temp_ind]

        if (sum(loc_lut_ind)==1) {
          parent_ind = loc_lut_ind
          parent_lev = dis_lut$level[parent_ind]
        } else {
          parent_lev = 0
          parent_ind = dis_lut$ABBV_2==date_ranges$ABBV_2[country_ind][1]
          if (min(dis_lut$level[loc_lut_ind])==2) {
            parent_lev = 2
            parent_ind = parent_ind & dis_lut$level==2
          } else {
            for (lev in 3:(min(dis_lut$level[loc_lut_ind])+1)) {
              if (length(unique(dis_lut[[paste0("NAME_", lev)]][loc_lut_ind]))!=1) {
                parent_lev = lev-1
                parent_ind = parent_ind & dis_lut$level==parent_lev
                break
              } else {
                parent_ind = parent_ind & dis_lut[[paste0("ABBV_", lev)]]==dis_lut[[paste0("ABBV_", lev)]][loc_lut_ind][1] & !is.na(dis_lut[[paste0("ABBV_", lev)]])
              }
            }
          }
        }
        source_info$level[ii] = parent_lev
        if (sum(parent_ind)!=1) {
          # browser()
          stop("Data source ", source_info$source_keys[ii], " does not have a root location.")
        }

        # record location name and master key of root location
        source_info$master_key[ii] = dis_lut$master_key[parent_ind]
        source_info$loc[ii] = dis_lut[[paste0("NAME_", parent_lev)]][parent_ind]
        source_info$cadence[ii] = date_ranges$cadence[temp_ind][1]
      }
      source_info$source_abbv = data_sources$source_abbv[match(source_info$source_keys, data_sources$source_key)]

      locations = unique(source_info$loc)
      for (loc in locations) {
        summary[[disease]][[country]][[loc]] = list()
        loc_ind = loc==source_info$loc
        for (source in which(loc_ind)) {
          source_key = source_info$source_keys[source]
          s_abbv     = source_info$source_abbv[source]
          summary[[disease]][[country]][[loc]][[s_abbv]] = list()
          summary[[disease]][[country]][[loc]][[s_abbv]]$source_info = data_sources[data_sources$source_key==source_key, ]
          source_ind = country_ind & date_ranges$source_key==source_key
          levels = unique(date_ranges$level[source_ind])
          summary[[disease]][[country]][[loc]][[s_abbv]]$data_summary = data.frame(levels=levels, num_regions=0, min_date=as.Date("2000-01-01"), max_date=as.Date("2000-01-01"))
          for (jj in 1:length(levels)) {
            lev_ind = source_ind & date_ranges$level==levels[jj]
            num_regions = sum(lev_ind)
            min_date = min(date_ranges$min_date[lev_ind])
            max_date = max(date_ranges$max_date[lev_ind])
            summary[[disease]][[country]][[loc]][[s_abbv]]$data_summary$num_regions[jj] = num_regions
            summary[[disease]][[country]][[loc]][[s_abbv]]$data_summary$min_date[jj] = min_date
            summary[[disease]][[country]][[loc]][[s_abbv]]$data_summary$max_date[jj] = max_date
          }
        } # source loop
      } # location loop
    } # country loop
  } # disease loop
  dbDisconnect(myDB)
  return(summary)
}


Proc_IncCol_BySource <- function(ident_data=NULL, data_source=NULL, factor=NULL, raw_col=NULL) {
  # for each data source, write a procedure that determines which data column to use and converts it to either 'cases' or 'percent'.  Then designate the appropriate column with a matching column name.
  if (is.null(raw_col)) {
    if (data_source %in% c(7,11)) {
      raw_name = "ILI"
    } else if (data_source %in% c(26)) {
      raw_name = "total_specimens"
    } else {
      raw_name = "cases"
    }

    if (data_source == 15) {
      # for each patch, choose between ILI_CASES and ILI_OUT by determining which has more periods of data
      n_ili_cases = sum(!is.na(ident_data$cases))
      n_ili_out   = sum(!is.na(ident_data$ili_out))
      if (n_ili_cases>n_ili_out) {
        ident_data$raw = ident_data$cases
      } else if (n_ili_cases<n_ili_out) {
        ident_data$raw = ident_data$ili_out
      } else {
        # sum both columns and use the one with more cases
        n_ili_cases = sum(ident_data$cases, na.rm=T)
        n_ili_out   = sum(ident_data$ili_out, na.rm=T)
        if (n_ili_cases>=n_ili_out) {
          ident_data$raw = ident_data$cases
        } else {
          ident_data$raw = ident_data$ili_out
        }
      }
    } else {
      # for all other data sources
      # grab the appropriate column
      ident_data$raw = ident_data[[raw_name]]
    }
  } else {
    # if raw_col has been specified, check that it matches a column of ident_dat
    if (raw_col %in% names(ident_data)) {
      ident_data$raw = ident_data[[raw_col]]
    } else {
      stop("The raw_col '", raw_col, "' is not a valid column name for data_source=", data_source, ".\n", sep="")
    }

  }
  # convert raw to epi
  ident_data$epi = ident_data$raw*factor
  # set NAs to 0
  ident_data$epi[is.na(ident_data$epi)] = 0
  return(ident_data)
}


Proc_IncCol_BySource_mydata <- function(mydata=NULL, data_source=NULL, cadence=NULL, mod_sum=NULL, col_units=NULL, unit_types=NULL, raw_col=NULL, data_sources=NULL) {
  # for each data source, write a procedure that determines which data column to use and converts it to either 'cases' or 'percent'.  Then designate the appropriate column with a matching column name.
  nperiodsPerYear = switch(cadence, 365, 52, 12, 1)

  if (!mod_sum) mod_ident = mydata$model$attr$identifier
  fit_ident = mydata$fit$attr$identifier

  data_col_string = data_sources$col_names[data_sources$source_key==data_source]
  data_col_names = strsplit(x=data_col_string, split=";")[[1]]
  col_units_string = data_sources$col_units[data_sources$source_key==data_source]
  col_units_names = strsplit(x=col_units_string, split=";")[[1]]

  if (!(data_source==15 & is.null(raw_col))) {
    # in all other cases

    if (is.null(raw_col)) {
      # if user has not specified the data metric, use first one
      raw_name = data_col_names[1]
      raw_unit = col_units_names[1]
    } else {
      # use raw_col to determine which data column to use
      if (raw_col %in% data_col_names) {
        name_ind = which(raw_col==data_col_names)
        raw_name = raw_col
        raw_unit = col_units_names[name_ind]
      } else {
        stop("The raw_col '", raw_col, "' is not a valid column name for data_source=", data_source, ".\n", sep="")
      }
    }

    DICE_unit = col_units$unit_type[col_units$col_unit==raw_unit]
    factor_string = unit_types$factor[unit_types$unit_type==DICE_unit]

    # use raw_units to assign raw vector and calculate raw->epi factor
    if (!mod_sum) {
      mydata$model[[mod_ident]]$raw = mydata$model[[mod_ident]][[raw_name]]
      # Assign variables needed for factor calculation
      pop = mydata$model$pop
      consultPerCapita = 2  # For now this is fixed.  Future this will be region/year specific
      mydata$model$factor = eval(parse(text=factor_string))
    }
    mydata$model$raw_units = raw_unit

    # repeat for fit-level
    mydata$fit$factor = rep(1, length(fit_ident))
    consultPerCapita = 2  # For now this is fixed.  Future this will be region/year specific
    for (ii in 1:length(fit_ident)) {
      mydata$fit[[fit_ident[ii]]]$raw = mydata$fit[[fit_ident[ii]]][[raw_name]]
      # Assign variables needed for factor calculation
      pop = mydata$fit$pop[ii]
      mydata$fit$factor[ii] = eval(parse(text=factor_string))
    }
    mydata$fit$raw_units = raw_unit

  } else {
    # auto-choose WHO metric to use
    if (!mod_sum) {
      # for each patch, choose between ILI_CASES and ILI_OUT by determining which has more periods of data
      n_ili_cases = sum(!is.na(mydata$model[[mod_ident]]$cases))
      n_ili_out   = sum(!is.na(mydata$model[[mod_ident]]$ili_out))
      if (n_ili_cases>n_ili_out) {
        raw_name = "cases"
      } else if (n_ili_cases<n_ili_out) {
        raw_name = "ili_out"
      } else {
        # sum both columns and use the one with more cases
        n_ili_cases = sum(mydata$model[[mod_ident]]$cases, na.rm=T)
        n_ili_out   = sum(mydata$model[[mod_ident]]$ili_out, na.rm=T)
        if (n_ili_cases>=n_ili_out) {
          raw_name = "cases"
        } else {
          raw_name = "ili_out"
        }
      }
      mydata$model[[mod_ident]]$raw = mydata$model[[mod_ident]][[raw_name]]
      mydata$model$raw_units = col_units_names[data_col_names==raw_name]
    }
    mydata$model$factor = 1

    # Repeat for each patch in $fit
    for (ii in 1:length(fit_ident)) {
      # for each patch, choose between ILI_CASES and ILI_OUT by determining which has more periods of data
      n_ili_cases = sum(!is.na(mydata$fit[[fit_ident[ii]]]$cases))
      n_ili_out   = sum(!is.na(mydata$fit[[fit_ident[ii]]]$ili_out))
      if (n_ili_cases>n_ili_out) {
        raw_name = "cases"
      } else if (n_ili_cases<n_ili_out) {
        raw_name = "ili_out"
      } else {
        # if equal, sum both columns and use the one with more cases
        n_ili_cases = sum(mydata$fit[[fit_ident[ii]]]$cases, na.rm=T)
        n_ili_out   = sum(mydata$fit[[fit_ident[ii]]]$ili_out, na.rm=T)
        if (n_ili_cases>=n_ili_out) {
          raw_name = "cases"
        } else {
          raw_name = "ili_out"
        }
      }
      mydata$fit[[fit_ident[ii]]]$raw = mydata$fit[[fit_ident[ii]]][[raw_name]]
    }
    # determine raw_units and factor
    mydata$fit$raw_units = col_units_names[data_col_names==raw_name]
    mydata$fit$factor = rep(1, length(fit_ident))
  }
  # record which column was used for raw
  mydata$fit$raw_name = raw_name
  mydata$model$raw_name = raw_name
  
  return(mydata)
}


cdc_download <- function (CDCreg=NULL, years=NULL) {
  #'
  #' Pull ili mydata from CDC server (make 5 attempts)
  #'
  GoodPull = FALSE
  n = 0
  while (n < 5) {
    CDCmodel = try(ilinet(region = CDCreg, years = years), silent = TRUE)
    n = n + 1
    if (!is(CDCmodel, "try-error")) {
      GoodPull = TRUE
      break
    } else Sys.sleep(.1)
  }
  return(list(CDCdata=CDCmodel, GoodPull=GoodPull))
}



StandardizeLUT <- function(dis_lut=dis_lut, DICE_db=db_opts$DICE_db) {
  # if LUT is coming from PostGreSQL, capitalize some column names
  if (tolower(DICE_db)=="bsve") {
    col_names = names(dis_lut)
    for (key in c("name_", "abbv_", "id_")) {
      key_length = nchar(key)
      name_sub   = substr(col_names, start=1, stop=key_length)
      key_ind    = name_sub==key
      col_names[key_ind] = toupper(col_names[key_ind])
    }
    names(dis_lut) = col_names
  }

  return(dis_lut)
}


get.school.data <- function(master_keys=NULL, cadence=NULL, start.date=NULL, end.date=NULL, myDB=NULL) {

  cad_name = tolower(CadenceNum2Name(cadences=cadence))

  where_string = character()

  master_keys = unique(master_keys)
  if (!is.null(master_keys)) {
    where_string = paste0("master_key IN('",paste(master_keys, collapse="','"),"')")
  }
  if (!is.null(start.date)) {
    where_string = paste0(where_string, " AND date >= '", start.date, "'")
  }
  if (!is.null(end.date)) {
    where_string = paste0(where_string, " AND date <= '", end.date, "'")
  }

  if (length(where_string)>0) {
    query_string = paste0("SELECT * FROM school_",cad_name," WHERE ", where_string, "ORDER BY master_key, date")
  } else {
    query_string = paste0("SELECT * FROM ",cad_name,"_data ", "ORDER BY master_key, date")
  }

  school_data = dbGetQuery(myDB, statement=query_string)
  if (nrow(school_data)>0) {
    school_data$date = as.Date(school_data$date)
  }

  return(school_data)
}


get.census.data <- function(master_keys=NULL, start.date=NULL, end.date=NULL, myDB=NULL) {

  where_string = character()

  master_keys = unique(master_keys)
  if (!is.null(master_keys)) {
    where_string = paste0("master_key IN('",paste(master_keys, collapse="','"),"')")
  }
  if (!is.null(start.date)) {
    where_string = paste0(where_string, " AND date >= '", start.date, "'")
  }
  if (!is.null(end.date)) {
    where_string = paste0(where_string, " AND date <= '", end.date, "'")
  }

  if (length(where_string)>0) {
    query_string = paste0("SELECT * FROM pop_yearly WHERE ", where_string, "ORDER BY master_key, date")
  } else {
    query_string = paste0("SELECT * FROM pop_yearly ORDER BY master_key, date")
  }

  pop_data = dbGetQuery(myDB, statement=query_string)
  if (nrow(pop_data)>0) {
    pop_data$date = as.Date(pop_data$date)
  }

  return(pop_data)
}


get.onset.data <- function(master_keys=NULL, start.year=NULL, end.year=NULL, myDB=NULL) {

  where_string = character()

  master_keys = unique(master_keys)
  if (!is.null(master_keys)) {
    where_string = paste0("master_key IN('",paste(master_keys, collapse="','"),"')")
  }
  if (!is.null(start.year)) {
    if (length(where_string)>0) {
      where_string = paste0(where_string, " AND year >= '", start.year, "'")
    } else {
      where_string = paste0(where_string, " year >= '", start.year, "'")
    }
  }
  if (!is.null(end.year)) {
    if (length(where_string)>0) {
      where_string = paste0(where_string, " AND year <= '", end.year, "'")
    } else {
      where_string = paste0(where_string, " year <= '", end.year, "'")
    }
  }

  if (length(where_string)>0) {
    query_string = paste0("SELECT * FROM season_onset WHERE ", where_string, "ORDER BY master_key, year")
  } else {
    query_string = paste0("SELECT * FROM season_onset ORDER BY master_key, year")
  }

  onset_data = dbGetQuery(myDB, statement=query_string)
  return(onset_data)
}


get.cdc <- function(all_years=NULL, mod_level=2, fit_level=3, mod_name=c(NAME_2 = "USA"), myDB=NULL) {

  #' Retrieve incidence data from CDC server.
  #'
  #' This function is generally called by \code{\link{get.DICE.data}}.
  #' Augment/overwrite all_years with incidence data retrieved from the CDC server.  Also augment/overwrite specific humidity with NASA data retrieved from the Predictive Science server.
  #'
  #' @param all_years Data structure resulting from a call to \code{\link{get_disease_data}}
  #' @param mod_level An integer describing the spatial level of the model data.(Default value is 2)  Levels: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City.  \pkg{dice} currently has data at levels 2-4 for CDC and GFT.
  #' @param fit_level An integer describing the spatial level of the fits used to construct the model-level profile/forecast (Default value is 3, must be >= mod_level).
  #' @param mod_name A named-vector of character strings that specify which region is to be \strong{model}ed.  In other words, \code{mod_name} specifies the country, region, state, etc. of the mod_level region.  \code{mod_name} should be of the form \code{mod_name = c(NAME_2='a', NAME_3='b',..., NAME_i='x'} where i=mod_level and 'a', 'b',...,'x' are the appropriate level names. NAME_i='x' also accepts abbreviations.  Choose appropriate names from \code{\link{diceData}}.  For example, mod_name=c(NAME_2='United.States',NAME_3='Region4',NAME_4='North.Carolina') and mod_level=4 specifies North Carolina. To achieve the same result, use all abbreviations mod_name=c(NAME_2='USA',NAME_3='R4',NAME_4='NC') or a mix of names and abbreviations mod_name=c(NAME_2='USA',NAME_3='Region4',NAME_4='NC')
  #' @param myDB A connection to one of the databases.  Generally opened using \code{\link{OpenCon}}.

  # Get ili mydata
  GoodPull_mod = FALSE
  GoodPull_fit = FALSE

  # Attempt to download 'model' ILI mydata from CDC server
  if (mod_level == 2) {
    CDCreg = "national"
  } else if (mod_level == 3) {
    CDCreg = "hhs"
    sub_region = all_years$lut$ID_3[all_years$lut$level==mod_level]
  } else if (mod_level == 4) {
    CDCreg = "state"
    sub_region = all_years$lut$NAME_4[all_years$lut$level==mod_level]
  }

  cat("Replacing database incidence with CDC-server incidence.....")

  # Pull ili mydata from CDC server (make 5 attempts)
  n = 0
  while (n < 5) {
    CDCmodel = try(ilinet(region = CDCreg), silent = TRUE)
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
    } else Sys.sleep(0.1)
  }
  # If mydata download is unsuccessful, print warning and revert to SQL server data.
  if (!GoodPull_mod) {
    print("\nWARNING: Download from CDC server was unsuccessful for model-level ILI. Reverting to ILI data in SQL database.  SQL data is likely less up-to-date than the CDC server.")
  }

  # Attempt to download 'fit' ILI mydata from CDC server
  if (fit_level == 2) {
    CDCreg = "national"
  } else if (fit_level == 3) {
    CDCreg = "hhs"
    sub_region = all_years$lut$ID_3[all_years$lut$level==fit_level]
  } else if (fit_level == 4) {
    CDCreg = "state"
  }
  # Pull ili mydata from CDC server (make 5 attempts)
  n = 0
  while (n < 5) {
    CDCfit = try(ilinet(region = CDCreg), silent = TRUE)
    n = n + 1
    if (!is(CDCfit, "try-error")) {
      GoodPull_fit = TRUE
      if (fit_level == 3) {
        # reduce CDCfit to specified regions
        reg_ind = CDCfit$region %in% paste0("Region ", sub_region)
        CDCfit = CDCfit[reg_ind, ]
      }
      break
    } else Sys.sleep(0.1)
  }
  # If mydata download is unsuccessful, print warning and revert to saved mydata.
  if (!GoodPull_fit) {
    print("WARNING: Download from CDC server was unsuccessful for fit-level ILI. Reverting to ILI data in SQL database.  SQL data is likely less up-to-date than the CDC server.")
  }


  if (GoodPull_mod) {
    # Overwrite/augment all_years
    if (mod_level<4) {
      mod_df = CDCmodel[, c("year", "week", "weighted_ili", "ilitotal")]
    } else {
      mod_df = CDCmodel[, c("year", "week", "unweighted_ili", "ilitotal")]
    }
    names(mod_df) = c("year", "week", "ILI2", "cases2")

    ident = all_years$lut$identifier[all_years$lut$level==mod_level]

    temp_epi = merge(mod_df, all_years$Weekly[[ident]]$CDC, all=T, sort=F)

    # fill empty dates
    na_ind = is.na(temp_epi$date)
    if (any(na_ind)) {
      for (ii in which(na_ind)) {
        temp_epi$date[ii] = CDCweek2date(CDCweek=temp_epi$week[ii], year=temp_epi$year[ii])
      }
    }
    # if any gaps in CDC-server data, fill with SQL data
    na_ind = is.na(temp_epi$ILI2)
    if (any(na_ind)) {
      temp_epi$ILI2[na_ind] = temp_epi$ILI[na_ind]
      temp_epi$cases2[na_ind] = temp_epi$cases[na_ind]
    }
    # remove SQL data columns
    temp_epi = temp_epi[, c("date", "year", "week", "ILI2", "cases2")]
    names(temp_epi) = c("date", "year", "week", "ILI", "cases")
    # re-sort by date
    temp_epi = temp_epi[order(temp_epi$date), ]
    # update all_years
    all_years$Weekly[[ident]]$CDC = temp_epi
  }  # end if GoodPull_mod


  if (GoodPull_fit) {
    if (fit_level == 4) {
      # sort out only the states requested
      req_names = all_years$lut$NAME_4[all_years$lut$level==fit_level]
      keep_ind = CDCfit$region %in% req_names
      CDCfit = CDCfit[keep_ind, ]
    }
    nregions = length(unique(CDCfit$region))
    # Process CDC 'fit' download
    if (nregions > 1) {
      if (fit_level==4) {
        # choose a sample state that is not Puerto Rico
        cdc_states = unique(CDCfit$region)
        sample_state = cdc_states[1]
        if (sample_state=="Puerto Rico") {
          sample_state = cdc_states[2]
        }
        sample_ind = CDCfit$region==sample_state
        tempYear = CDCfit$year[sample_ind]
        tempWeek = CDCfit$week[sample_ind]
      } else {
        tempYear = CDCfit$year[seq(from=1, to=length(CDCfit$year), by=nregions)]
        tempWeek = CDCfit$week[seq(from=1, to=length(CDCfit$year), by=nregions)]
      }

      # reshape CDC download to DICE-matrix form
      if (fit_level == 4) {
        FitILI = array(data = NA, dim = c(nregions, length(tempYear)), dimnames = list(req_names, NULL))
        FitCases = array(data = NA, dim = c(nregions, length(tempYear)), dimnames = list(req_names, NULL))
        StateTotal = FitCases
        # match dates for each row of FitILI
        for (ii in 1:nregions) {
          temp_fit = CDCfit[CDCfit$region==req_names[ii], ]
          match_df = merge(x=temp_fit, y=data.frame(year=tempYear, week=tempWeek), all.y=T, sort=F)
          match_df = match_df[order(match_df$year, match_df$week), ]
          FitILI[ii, ]     = match_df$unweighted_ili
          FitCases[ii, ]   = match_df$ilitotal
          StateTotal[ii, ] = match_df$total_patients
        }
        StateILITOTAL = t(FitCases)
        StateTotal    = t(StateTotal)
      } else {
        FitILI = array(data = CDCfit$weighted_ili, dim = c(nregions, length(CDCfit$year)/nregions))
        FitCases = array(data = CDCfit$ilitotal, dim = c(nregions, length(CDCfit$year)/nregions))
      }
    } else {
      if (fit_level==4) {
        FitILI = matrix(CDCfit$unweighted_ili, nrow=1)
      } else {
        FitILI = matrix(CDCfit$weighted_ili, nrow=1)
      }
      FitCases = matrix(CDCfit$ilitotal, nrow=1)
      tempYear = CDCfit$year
      tempWeek = CDCfit$week
    }


    # If working with state-level, check for columns with NAs
    if (fit_level == 4) {
      # check for states that are missing completely (Florida)
      if ((mod_level==2) | (mod_level==3 && all_years$lut$ABBV_3[all_years$lut$level==3]=="R4")) {
        if (!("Florida") %in% row.names(FitILI)) {
          FitILI = rbind(FitILI, matrix(NA, nrow=1, ncol=ncol(FitILI), dimnames=list(c("Florida"), NULL)))
          FitCases = rbind(FitCases, matrix(NA, nrow=1, ncol=ncol(FitCases), dimnames=list(c("Florida"), NULL)))
          # flu_lut = dbReadTable(myDB, "flu_lut")
          flu_lut = dbGetQuery(myDB, statement=paste0("SELECT * FROM flu_lut ORDER BY master_key;"))
          flu_lut = StandardizeLUT(dis_lut=flu_lut, DICE_db="bsve")
          all_years$lut = rbind(all_years$lut, flu_lut[flu_lut$level==4 & flu_lut$NAME_4=="Florida", ])
        }
      }
      # determine which states have NAs
      NA_ind = apply(FitILI, 1, function(x) any(is.na(x)))
      # if any states are missing ILI data, attempt to reconstruct their ILI from region totals
      if (any(NA_ind)) {
        regions = unique(all_years$lut$ID_3[all_years$lut$level==fit_level][NA_ind])

        n = 0
        while (n < 5) {
          CDCregs = try(ilinet(region = "hhs"), silent = TRUE)
          n = n + 1
          if (!is(CDCregs, "try-error")) {
            GoodPull_regs = TRUE
            # reduce mydata to specified regions
            reg_ind = CDCregs$region %in% paste0("Region ", regions)
            CDCregs = CDCregs[reg_ind, ]
            break
          } else Sys.sleep(5)
        }

        # grab state ILITOTAL and TOTAL PATIENTS (now done previously)
        # StateILITOTAL = array(data = CDCfit$ilitotal, dim = c(nregions, length(CDCfit$year)/nregions), dimnames = list(unique(CDCfit$region), NULL))
        # StateILITOTAL = t(StateILITOTAL[req_names, ])
        # StateTotal = array(data = CDCfit$total_patients, dim = c(nregions, length(CDCfit$year)/nregions), dimnames = list(unique(CDCfit$region), NULL))
        # StateTotal = t(StateTotal[req_names, ])
        state_reg  = all_years$lut$ID_3[match(req_names, all_years$lut$NAME_4)]
        # reshape region mydata to a 'raw' type array
        nregions = length(regions)
        if (nregions == 1) {
          nweeks = nrow(unique(CDCregs[, c("week", "year")]))
          RegYear = CDCregs$year
          RegWeek = CDCregs$week
          RegILITOTAL = array(data = CDCregs$ilitotal, dim = c(nweeks, nregions))
          RegTotal = array(data = CDCregs$total_patients, dim = c(nweeks, nregions))
        } else {
          RegILITOTAL = t(array(data = CDCregs$ilitotal, dim = c(nregions, length(CDCregs$year)/nregions)))
          RegTotal = t(array(data = CDCregs$total_patients, dim = c(nregions, length(CDCregs$year)/nregions)))
          RegYear = CDCregs$year[seq(from = 1, to = length(CDCregs$year), by = nregions)]
          RegWeek = CDCregs$week[seq(from = 1, to = length(CDCregs$year), by = nregions)]
        }


        # go through each week and attempt to reconstruct state ILI
        state_list = row.names(FitILI)
        for (week in 1:ncol(FitILI)) {
          NA_ind = is.na(FitILI[, week])
          epi_week = tempWeek[week]
          epi_year = tempYear[week]
          if (any(NA_ind)) {
            NA_reg = all_years$lut$ID_3[NA_ind]
            # if more than one state is missing from a region, do not attempt to reconstruct
            dup_ind = duplicated(NA_reg)
            if (any(dup_ind)) {
              dup_reg = NA_reg[dup_ind]
              # cat('Warning!!!! Region(s) ',dup_reg,' are missing mydata from more than one state.  Individual state ILI cannot be reconstructed
              # from region totals.')
              remove_ind = NA_reg %in% dup_reg
              NA_reg = NA_reg[!remove_ind]
              NA_ind = NA_ind & all_years$lut$ID_3 %in% NA_reg
            }
          }

          if (any(NA_ind)) {
            # use region totals to reconstruct state ILI
            reg_ind = RegWeek==epi_week & RegYear==epi_year
            for (state_num in which(NA_ind)) {
              state = state_list[state_num]
              region = all_years$lut$ID_3[state==all_years$lut$NAME_4 & all_years$lut$level==4]
              state_tot = RegTotal[reg_ind, regions==region] - sum(StateTotal[week, state_reg==region], na.rm = TRUE)
              state_ilitot = RegILITOTAL[reg_ind, regions == region] - sum(StateILITOTAL[week, state_reg==region], na.rm = TRUE)
              FitILI[state_num, week] = state_ilitot/state_tot * 100
              FitCases[state_num, week] = state_ilitot
            }
          }
        }
      }  # end if any(state NA)
    }  # end if state NA fill

    # overwrite/augment all_years incidence data
    reg_ind = which(all_years$lut$level==fit_level)
    for (ii in 1:length(reg_ind)) {
      ident    = all_years$lut$identifier[reg_ind[ii]]
      fit_df   = data.frame(year=tempYear, week=tempWeek, ILI2=FitILI[ii, ], cases2=FitCases[ii, ])
      if (ident %in% names(all_years$Weekly)) {
        temp_epi = merge(fit_df, all_years$Weekly[[ident]]$CDC, all=T, sort=F)
      } else {
        temp_epi = fit_df
      }


      # fill any empty dates
      if ("date" %in% names(temp_epi)) {
        na_ind = is.na(temp_epi$date)
      } else {
        na_ind = rep(TRUE, nrow(temp_epi))
        temp_epi$date = as.Date("1900-01-01")
      }
      if (any(na_ind)) {
        for (jj in which(na_ind)) {
          temp_epi$date[jj] = CDCweek2date(CDCweek=temp_epi$week[jj], year=temp_epi$year[jj])
        }
      }

      # if any gaps in CDC-server data, fill with SQL data
      if (all(c("ILI", "cases") %in% names(temp_epi))) {
        na_ind = is.na(temp_epi$ILI2)
        if (any(na_ind)) {
          temp_epi$ILI2[na_ind] = temp_epi$ILI[na_ind]
          temp_epi$cases2[na_ind] = temp_epi$cases[na_ind]
        }
      }

      # remove SQL data columns
      keep_names = c("date", "year", "week", "ILI2", "cases2")
      keep_ind = keep_names %in% names(temp_epi)
      temp_epi = temp_epi[, c("date", "year", "week", "ILI2", "cases2")[keep_ind]]
      names(temp_epi) = c("date", "year", "week", "ILI", "cases")[keep_ind]
      # re-sort by date
      temp_epi = temp_epi[order(temp_epi$date), ]
      # update all_years
      all_years$Weekly[[ident]]$CDC = temp_epi
    }
  }  # end if GoodPull_fit

  cat("Complete. \n")
  cat("Replacing database NOAA-sh with NASA-sh.....")

  # retrieve SH for model region (from NASA)
  if (mod_level == 2) {
    CDCreg = "national"
    sub_region = "national"
  } else if (mod_level == 3) {
    CDCreg = "hhs"
    sub_region = all_years$lut$ID_3[all_years$lut$level==mod_level]
  } else if (mod_level == 4) {
    CDCreg = "state"
    sub_region = gsub(" ",".", all_years$lut$NAME_4[all_years$lut$level==mod_level], fixed=TRUE)
  }
  ident = all_years$lut$identifier[all_years$lut$level==mod_level]
  temp_week = Date2CDCweek(Sys.Date())
  sh_df = get.sh.nasa(region=CDCreg, sub_region=sub_region, start_year=all_years$Weekly[[ident]]$CDC$year[1], start_week=all_years$Weekly[[ident]]$CDC$week[1], end_year=temp_week$CDCyear, end_week=temp_week$CDCweek)
  # match dates from all_years and sh_df
  names(sh_df)[4] = "sh2"
  temp_clim = merge(x=sh_df, y=all_years$Weekly[[ident]]$noaa_weekly, all=T, sort=F)

  # if any gaps in NASA SH data, fill with NOAA SH data
  if ("sh" %in% names(temp_clim)) {
    na_ind = is.na(temp_clim$sh2)
    if (any(na_ind)) {
      temp_clim$sh2[na_ind] = temp_clim$sh[na_ind]
    }
  }

  # remove NOAA sh column
  keep_names = c("date", "sh2", "temp", "precip", "press", "rh")
  keep_ind  = keep_names %in% names(temp_clim)
  temp_clim = temp_clim[, keep_names[keep_ind]]
  names(temp_clim) = c("date", "sh", "temp", "precip", "press", "rh")[keep_ind]
  # re-sort by date
  temp_clim = temp_clim[order(temp_clim$date), ]
  # update all_years
  all_years$Weekly[[ident]]$noaa_weekly = temp_clim


  # retrieve SH for all_years_epi fit regions (from NASA)
  if (fit_level == 2) {
    CDCreg = "national"
    sub_region = "national"
  } else if (fit_level == 3) {
    CDCreg = "hhs"
    sub_region = all_years$lut$ID_3[all_years$lut$level==fit_level]
  } else if (fit_level == 4) {
    CDCreg = "state"
    sub_region = gsub(" ", ".", all_years$lut$NAME_4[all_years$lut$level==fit_level], fixed=TRUE)
  }
  sh_df_all = get.sh.nasa(region=CDCreg, sub_region=sub_region, start_year=all_years$Weekly[[ident]]$CDC$year[1], start_week=all_years$Weekly[[ident]]$CDC$week[1], end_year=temp_week$CDCyear, end_week=temp_week$CDCweek)

  # match dates from all_years and sh_df
  reg_ind = all_years$lut$level==fit_level
  for (ii in 1:sum(reg_ind)) {
    sh_df = sh_df_all[, c(1:3, 3+ii)]
    ident = all_years$lut$identifier[reg_ind][ii]
    if ("noaa_weekly" %in% names(all_years$Weekly[[ident]])) {
      names(sh_df)[4] = "sh2"
      temp_clim = merge(x=sh_df, y=all_years$Weekly[[ident]]$noaa_weekly, all=T, sort=F)

      # if any gaps in NASA SH data, fill with NOAA SH data
      if ("sh" %in% names(temp_clim)) {
        na_ind = is.na(temp_clim$sh2)
        if (any(na_ind)) {
          temp_clim$sh2[na_ind] = temp_clim$sh[na_ind]
        }
      }

      # remove SQL data columns
      keep_names = c("date", "sh2", "temp", "precip", "press", "rh")
      keep_ind  = keep_names %in% names(temp_clim)
      temp_clim = temp_clim[, keep_names[keep_ind]]
      names(temp_clim) = c("date", "sh", "temp", "precip", "press", "rh")[keep_ind]
      # re-sort by date
      temp_clim = temp_clim[order(temp_clim$date), ]
    } else {
      # rename sh column
      names(sh_df)[4] = "sh"
      # keep only date and sh
      temp_clim = sh_df[, c(1,4)]
      # add placeholder columns for temp and precip
      temp_clim$temp = as.numeric(NA)
      temp_clim$precip = as.numeric(NA)
    }

    # update all_years
    all_years$Weekly[[ident]]$noaa_weekly = temp_clim
  }
  cat("Complete. \n")

  return(all_years)
}


PullSeason_AllYears <- function(all_years_epi=NULL, season_start_date=NULL, season_end_date=NULL) {

  # index dates from all_years_epi
  date_ind = all_years_epi$dates>=season_start_date & all_years_epi$dates<=season_end_date

  out_epi = list()
  names1 = names(all_years_epi)

  for (list_name in names1) {
    if (list_name=="model") {
      out_epi$model = all_years_epi$model[date_ind, ]
    } else if (list_name=="fit") {
      names2 = names(all_years_epi$fit)
      out_epi$fit = list()
      if (mydata$fit$nregions>1) {
        for (fit_name in names2) {
          out_epi$fit[[fit_name]] = all_years_epi$fit[[fit_name]][date_ind, ]
        }
      } else {
        for (fit_name in names2) {
          out_epi$fit[[fit_name]] = all_years_epi$fit[[fit_name]][date_ind]
        }
      }
    } else {
      out_epi[[list_name]] = all_years_epi[[list_name]][date_ind]
    }
  }
  return(out_epi)
}


get_flu_poc_units <- function() {
  # define flu point-of-care units
  flu_poc_units = c("specimens tested", "valid specimens", "specimens positive for A", "specimens positive for B", "normalized total specimens", "normalized positive A", "normalized positive B")
  return(flu_poc_units)
}


