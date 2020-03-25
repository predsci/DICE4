##
## Driver and Functions for fitting and forecasting User provided incidence data
##

runDICEUserData <-
  function(filename = 'data.csv',
           pop = 1e4,
           epi_model = 1,
           Tg = 3,
           sigma = NULL) {
    #' Main Driver for using \code{DICE} for Fitting User Provided Incidence
    #'
    #' \code{runDICEUserData} reads in a user provided incidence file and uses that along with the user
    #' provided population size and value for the generation time, Tg, and if relevant the latent period
    #' sigma to fit, and if desired, forecast the incidence.
    #' Data cadence is arbitrary but at most can be monthly
    #' We support only S-I-R and S-E-I-R models for a single population with a fixed force of infection
    #' The incidence file must have two columns: dates and cases
    #' All other parameters are set by the code
    #' @param filenmae String - input csv file name.  Default 'data.csv'
    #' @param Tg - Numeric, generation time in days. Default is 3 days
    #' @param pop - Integer population of the region for which incidence is provided
    #' @param sigma - inverse of of the latent period in days. Needed only for an SEIR model. Default NULL
    #' @param epi_model - integer 1 (SIR) or 2 (SEIR), default is
    #' @return results A list with the input and entire output of the run.
    #' @examples
    #' Run an SEIR model using the incidence file and assuming a population of 1 million people.
    #' The generation time and latent period are set to 2.6 days and 3 days respectively
    #'
    #' output <- runDICEUserData(filename = 'data.csv', pop = 1e6, epi_model = 2, Tg = 2.6, sigma = 3.)
    #'
    #' Run an SIR model using the incidence file and assuming a population of 10,000 people.
    #' The generation time is set to 3 days. (No need to define a latent period.)
    #' output <- runDICEUserData(filename = 'data.csv', pop = 1e5, epi_model = 1, Tg = 3.)
    #'
    
    ## Set the MCMC parameters
    
    nreal = 1
    nMCMC = 1e6
    nlines = 1e4
    plot = 1
    
    device = 'pdf'
 
   # Build the mydata list
   
   mydata <- build.mydata(filename = filename, pop = pop, epi_model = epi_model, Tg = Tg, sigma = sigma)

   #
   # Build a name for the sub-directory that will hold the results
   #
   myName = Sys.time()
   myName = gsub(" ", "_", myName)
   subDir = paste0(mydata$model$name, "_", myName, "_Temp_", mydata$Temp, "_Tg_", Tg, "/")

   mydata$subDir = subDir

   # if necessary create this directory
   if (!dir.exists(subDir)) {
   	dir.create(subDir)
   }

   # Plot the incidence - there is no historic data 
   
   err <- plotUserDisease(mydata = mydata, device = device)

   ##
   ## Mechanistic Modeling
## Pack the information for the run

   par_names <- set.user.data.param.list()

   opt.list <- set.user.data.opt.list(mydata = mydata)

   run.list <- set.run.list(nreal = nreal, nMCMC = nMCMC, nlines = nlines, device = device, subDir = mydata$subDir, plot = plot)
   ## Fit the data
   
   output <- fit.user.data(mydata = mydata, par_names = par_names, opt.list = opt.list, run.list = run.list, iseed = NULL)
   return(output)
   
  }


plotUserDisease <-

  function(mydata = NULL,
           run.list = NULL,
           device = 'png') {
    #' Plot The Time Series of the User Provided Incidence
    #'
    #' \code{plotUserDisease} Plots the user provided disease time series 
    #' @param mydata A list with all the available data for this \pkg{DICE} run
    #' @param run.list a list with various parameters for the run
    #' @param device - 'pdf or 'png'.  Default is 'png'
    #' plotUserDisease(mydata = mydata, run.list = run.list, device = 'pdf')
    #' @return err=0 if plot was created
    #'
    
    if (is.null(device))
      device = "png"
    
    subDir = mydata$subDir
    
    if (is.null(subDir)) {
      subDir = getwd()
    }
    
    if (!dir.exists(subDir)) {
      dir.create(subDir)
    }
    
    nperiods = mydata$nperiods
    
    dates = mydata$dates
    x.axis.label = format(dates, "%b-%d")
    
    myName = mydata$model$name
    myName = gsub(" ", "", myName)
    
    if (tolower(device) == "pdf") {
   	filename = paste0(subDir, "/", myName, "-", mydata$disease, "-incidence", ".pdf")
   	pdf(file = filename, onefile = TRUE, width = 15, height = 9)
   } else if (tolower(device) == "png") {
   	filename = paste0(subDir, "/", myName, "-", mydata$disease, "-incidence", ".png")
   	png(file = filename, width = 1500, height = 900)
   } else {
   	dev.next()
   	dev.new()
   }

   cat("\n For a Plot of Incidence See: ", filename, "\n\n")

   nrow = ncol = 1
   par(mfrow = c(nrow, 1), mar = c(5, 5, 3, 2))

   FY = mydata$FY

   ylab = " # Cases"

   tps = cumsum(mydata$ndays)

   plot(tps, mydata$model$raw, type = "l", col = "red", xaxt = "n", xlab = "", ylab = ylab, lwd = 2)
   lines(tps, mydata$model$raw, type = "p", col = "red", pch = 19)
   axis(1, at = tps, labels = x.axis.label, las = 2)

   legend("topleft", c(mydata$model$name, paste0("Year ", FY)), bty = "n")

   if (tolower(device) == "pdf" || tolower(device) == "png") 
   	dev.off()

   err = 0
   return(err)
   
  }


build.mydata <-
  function(filename = "data.csv",
           pop = 1e4,
           epi_model = 1,
           Tg = 3,
           sigma = NULL) {
    #' Read incidence data file and Build the DICE data list, mydata
    #'
    #' \code{build.mydata} Reads the user provided csv file with two columns: date and cases
    #' The date format is: Year-month-day.
    #' cases - integer number of cases
    #' cadence can be anything including irregular
    #' The code builds and populates the mydata DICE list
    #' @param filenmae String - input csv file name.  Default 'data.csv'
    #' @param pop - Integer, population of the region for which incidence is provided 
    #' @param epi_model - integer 1 (SIR) or 2 (SEIR), default is 1
    #' @param Tg - Numeric, generation time in days. Default is 3 days
    #' @param sigma - inverse of of the latent period. Needed only for an SEIR model. Default 5 days
    #' @return mydata - a DICE list
    #' @examples
    #' mydata <- build.mydata(filename = "data.csv", pop = 1e5, epi_model = 2, Tg = 3., sigma = 5.)
    
    user.data = read.csv(file = filename, sep = ",")
    
    raw = user.data$cases
    
    dates = as.Date(user.data$date, format = '%m/%d/%y')
    
    cases = raw
    
    cases[is.na(cases)] <- 0
    
    nperiods = length(dates)
    
    ## Find how many data points we have and how many we may need to forecast
    
    nperiodsData = trimdata.in(longvec = raw)
    
    nperiodsDataUse = nperiodsFit = nperiodsData
    
    mydata = list()
    
    model = list()
    
    mydata$dates = as.Date(dates, format = '%Y-%m-%d')
    
    mydata$years = year(dates)
    
    mydata$days  = yday(dates)  # day of year
    
    mydata$months = month(dates)
    
    mydata$weeks = epiweek(dates)
    
    mydata$nperiods = nperiods
    
    mydata$nperiodsData = nperiodsData
    
    mydata$nperiodsDataUse = nperiodsDataUse
    
    mydata$nperiodsFit = nperiodsData
    
    # build the cadence between the dates
    # the first one is the median of all the other ones
    
    ndays = as.numeric(diff(dates, lag = 1))
    
    ndays0 = median(ndays)
    
    mydata$ndays = c(ndays0, ndays)
    
    mydata$season = year(dates)[1]
    
    year.end =  year(dates)[nperiods]
    
    year.start = year(dates)[1]
    
    if (year.end > year.start) {
      mydata$FY = paste0(year.start, '-', year.end)
    } else {
      mydata$FY = year.start
    }
    
    # create an array of zeros - will be used later
    
    zero.vec = rep(0, nperiods)
    
    # determine cadence
    if (all(mydata$ndays == 1)) {
      cadence = "Daily"
    } else if (all(mydata$ndays == 7)) {
      cadence = 'Weekly'
    } else if (all(mydata$ndays >= 28 && mydata$ndays <= 31)) {
      cadence = 'Monthly'
    } else {
      cadence = 'nonuniform'
    }
    
    mydata$cadence = cadence
    
    mydata$disease = 'unknown'
    
    mydata$dataName = 'user_data'
    
    mydata$data_source = 'user'
    
    mydata$data_desc = NULL
    
    mydata$method = 'mech'
    
    if ((nperiods * max(cases)) > 1e6) {
      Temp = 100
    } else if ((nperiods * max(cases)) >= 1e4 &&
               (nperiods * max(cases)) < 1e6) {
      Temp = 10
    } else {
      Temp = 1
    }
    mydata$Temp = Temp
    
    mydata$epi_model = epi_model
    
    mydata$single = 1
    
    mydata$imodel = 4
    
    mydata$prior = 0
    
    mydata$da = 0
    
    mydata$fit_level = mydata$mod_level = 'unknown'
    
    mydata$Tg = Tg
    
    mydata$sigma = sigma
    
    ## Build the model list
    
    model$level = 'unknown'
    
    model$name = 'user_data'
    
    model$attr = NULL
    
    model$pop = pop
    
    model$raw_units = '# of cases'
    
    model$factor = 1
    
    model$wght = zero.vec
    
    model$wght[1:nperiodsFit] = 1.0
    
    model$raw   = raw
    
    model$cases = cases
    
    model$epi   = cases
    
    model$gamaepi = lgamma((mydata$model$epi + 1))
    
    model$sh = zero.vec
    
    model$temp = zero.vec
    
    model$precip = zero.vec
    
    model$school = zero.vec
    
    model$attr = NULL
    
    mydata$model = model
    
    fit = list()
    
    fit$names = 'user_data'
    
    fit$nregions = 1
    
    mydata$fit = fit
    
    
    return(mydata)
  }

trimdata.in <- function(longvec) {
#' Trim an array of raw cases to find the last observed data point
#'
#' @param longvec array with raw cases numbers
#'
#' @return nperiodsData the number of last observed data point 
#' @export
#'
#' @examples 
#' nperiodsData <- trimdata.in(longvec = cases)
#' 
  
  noobs = length(longvec)
  last  = noobs
  
  while (is.na(longvec[last]))
    last = last - 1
  
  return(last)
}


##
## General setup of parameter lists routines
##

set.user.data.param.list <- function() {
  #' Short Parameter List - Fixed Force of Infection Case
  #'
  #'\code{set.user.data.param.list} creates a list with parameters that the \pkg{DICE} code recognizes.
  #' @return An array with parameter names - the order of parameters will set the order for the min/max arrays also
  #' and the 'mask' for which parameters are optimized (or not)
  #' @examples
  #' set.user.data.param.list()
  
  par_names <-
    c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd")
  
  return(par_names)
  
}

set.user.data.opt.list <- function(mydata = NULL) {
  #' Create a Logical List  with TRUE/FALSE values for parameter optimization
  #'
  #' \pkg{DICE} has a list of model parameters it recognizes, for both uncoupled and coupled runs.
  #' This function sets the values of these parameters to either TRUE or FALSE for the simple
  #' case of user provided data
  #' @param mydata - The \code{DICE}  data list
  #' @examples
  #' set.opt.list{mydata = mydata}
  #'
 
 opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE, pC = TRUE, t0 = TRUE, seed = TRUE, 
 	e_bckgrnd = TRUE)

 return(opt.list)
  
}


set.user.data.model <-
  function(mydata = NULL,
           par_names = NULL,
           opt.list = NULL,
           run.list = NULL) {
    #' Setup of Parameters for an MCMC procedure
    #'
    #' \code{set.user.data.model} Creates the arrays for modeling of user 
    #' provided incidence data with min/max values, step size,
    #' and initial guess/default values for all the parameters.
    #' @param par_names A character array with parameter names
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    #'   supported by \pkg{DICE}. 
    #' @return A list with min, max, step size and initial guess/default values for the sir parameters
    #' @examples
    #' set.user.data.model(mydata = mydata, par_names = par_names,
    #' run.list = run.list, opt.list = opt.list)
    
 
   nparam = length(par_names)

   parmin = parmax = pardx = par = logvec = rep(0, length = nparam)
   names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names
   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

   NH = mydata$model$pop
   R0 = 1.4
   t0 = 1
   pC = 0.005
   seed = 1
   sigma = mydata$sigma
   Tg = mydata$Tg

   if (is.null(Tg)) 
   	Tg = 3
   if (is.null(sigma)) 
   	sigma = 5
   ## Population, Tg and sigma
   
   par["NH"] = parmin["NH"] = parmax["NH"] = NH
   par["Tg"] = parmin["Tg"] = parmax["Tg"] = Tg
   par["sigma"] = parmin["sigma"] = parmax["sigma"] = sigma
   ## R0
   
   parmin["R0"] = 1.1
   parmax["R0"] = 7
   par["R0"] = runif(1, 1.2, 2)

   ## pC
   parmin["pC"] = 1e-06
   parmax["pC"] = 1
   par["pC"] = pC * runif(1, 0.8, 1.2)

   nperiodsData = mydata$nperiodsData
   nperiodsData2 = round(nperiodsData/2)
   sum.days = sum(mydata$ndays[1:(nperiodsData2)])
   ## Time of first infection
   parmin["t0"] = 0.001
   parmax["t0"] = sum.days
   par["t0"] = t0

   ## Initial Number of cases
   
   parmin["seed"] = 1
   parmax["seed"] = 0.001 * mydata$model$pop
   par["seed"] = 10 * runif(1, 0.8, 1.2)

   est_bckgrnd = mean(mydata$model$cases[1:3], na.rm = TRUE)
   parmin["e_bckgrnd"] = round(est_bckgrnd * 1)
   parmax["e_bckgrnd"] = round(est_bckgrnd * 10)
   par["e_bckgrnd"] = est_bckgrnd

   ##
   dx = 0.01
   pardx = c(NH = dx, Tg = dx, R0 = dx, sigma = dx, pC = dx, t0 = dx, seed = dx, e_bckgrnd = dx)


   logbase = 10 #use log base 10 when needed
   logvec[1:nparam] = 1
   logvec["t0"] = 0 #linear for t0

   nlines = run.list$nlines

   tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))

   setup = list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam, nopt = nopt, par_names = par_names, tab = tab)

   return(setup)
    
  }


plotUserDataMECH <-
  function(mydata = NULL,
           tables.mod = NULL,
           mod_id = NULL,
           ireal = 1,
           run.list = NULL,
           idevice = 1) {
    #' Plot Mechanistic Results of Fitting/Forecasting User Incidence Data
    #'
    #' \code{plotUserDataMECH} Creates a PDF or PNG file with the results of the MCMC fits to 
    #' incidence data
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param tables.mod - a table with the data and the results for the model fit. 
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param mod_id The abbreviation of the states/regions-model level
    #' @param ireal - Numeric, the number of the MCMC chain, default is 1
    #' @param idevice Integer the index in the device array, default is 1
    #' @examples
    #' plotUserDataMECH(mydata = mydata, tables.mod = tables.mod,
    #' mod_id = mod_id, ireal = ireal, run.list = run.list, idevice = 1)
    #' @return err=0 if plots were created
    #'
    #'
    
    device = run.list$device[idevice]
    if (is.null(device))
      device = "png"
    
    subDir = mydata$subDir
    
    Tg     = mydata$Tg
    
    isingle = mydata$isingle
    
    nperiods = mydata$nperiods
    
    nperiodsFit = mydata$nperiodsFit
    
    if (mydata$epi_model == 1 || mydata$epi_model == "SIR") {
      model_name = 'SIR'
    } else if (mydata$epi_model == 2 ||
               mydata$epi_model == "SEIR") {
      model_name = 'SEIR'
    } else {
      model_name = 'UNKNOWN'
    }
    
    cadence = mydata$cadence
    las = 1
    
    if (cadence == "Monthly") {
      month.names = month.abb[mydata$months]
      dates = month.names  #  paste0(month.name,'-',year.vec)
      ind = seq(from = 1, to = length(dates))
      xlab = ""
    } else if (cadence == "Weekly") {
      dates = mydata$weeks
      dates = sub("week", "EW", dates)
      ind = seq(1, length(dates), by = 4)
      dates = dates[ind]
      xlab = 'EW #'
    } else if (cadence == "Daily") {
      dates = mydata$dates
      ind = seq(1, length(dates), by = 7)
      dates = dates[ind]
      dates = as.Date(dates, format = '%Y-%m-%d')
      dates = format(dates, "%b-%d")
      dates = dates[ind]
      xlab = 'Date'
      las = 2
    } else if (cadence == "nonuniform") {
      dates = mydata$dates
      dates = as.Date(dates, format = '%Y-%m-%d')
      dates = format(dates, "%b-%d")
      ind = seq(1, length(dates), by = 1)
      xlab = 'Date'
      las = 2
    } else {
      dates = mydata$dates
      ind = seq(1, length(dates), by = 7)
      dates = dates[ind]
      dates = as.Date(dates, format = '%Y-%m-%d')
      dates = format(dates, "%b-%d")
      xlab = 'Date'
      las = 2
    }
    
    tps = cumsum(mydata$ndays)
    
    nperiods = mydata$nperiods
    
    nperiodsFit = mydata$nperiodsFit
    myName = mydata$dataName
    myName = gsub(" ", "", myName)
    
    if (tolower(device) == "pdf") {
   	filename = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
   	pdf(file = filename, onefile = TRUE, width = 18, height = 11)
   } else if (tolower(device) == "png") {
   	filename = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".png", sep = "")
   	png(file = filename, width = 1800, height = 1100)
   } else {
   	dev.next()
   	dev.new()
   }

   cat("\nFor a Plot of the Results See: ", filename, "\n")

   nrow = 1
   ncol = 1
   par(mar = c(6, 6, 3, 2), mfrow = c(nrow, ncol))

   ## Plot direct fit to Data
   
   ## subset only the state we are plotting
   observed = tables.mod$epi.obsrv[, mod_id]
   epi.null = tables.mod$epi.null[, , mod_id]

   model = tables.mod$seir.frcst[, , mod_id]
   model.lower = tables.mod$seir.lower[, , mod_id]
   model.upper = tables.mod$seir.upper[, , mod_id]

   if (mydata$model$name == "user_data") {
   	title = paste0("User Provided Data: ", mydata$FY, " Season")
   } else {
   	title = paste0(mydata$model$name, " - ", mydata$FY, " Season")
   }

   ymin = 0
   ymax = max(observed, model.upper[nperiodsFit, ], na.rm = TRUE)


   ylab = "# Cases"

   plot(tps, observed, type = "n", col = "black", lwd = 1, xlab = xlab, ylab = ylab, ylim = c(ymin, 
   	ymax), xlim = c(tps[1], tps[nperiods]), main = title, xaxt = "n")

   polygon(c(tps[1:nperiods], rev(tps[1:nperiods])), c(model.upper[nperiodsFit, ], rev(model.lower[nperiodsFit, ])), col = "lightblue", border = FALSE) #rgb(0, 0, 0.6, 0.2)
   lines(tps[1:nperiods], model[nperiodsFit, ], type = "l", col = "red", xlab = "", ylab = "")
   lines(tps[1:nperiods], observed[1:nperiods], type = "p", col = "red", lwd = 1, lty = 1, xlab = "", 
   	ylab = "")
   lines(tps[1:nperiodsFit], observed[1:nperiodsFit], type = "l", col = "black", lwd = 1, lty = 1, 
   	xlab = "", ylab = "")
   lines(tps[1:nperiodsFit], observed[1:nperiodsFit], type = "p", col = "black", lwd = 1, lty = 1, 
   	xlab = "", ylab = "", pch = 19)

   lines(tps[1:nperiods], epi.null[nperiodsFit, ], type = "l", col = "darkgrey", lwd = 1, lty = 1)
   rect(xleft = tps[nperiodsFit], ybottom = ymin, xright = tps[(nperiods)], ytop = ymax * 1.2, col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
   if (all(is.na(epi.null[nperiodsFit, ]))) {
   	legend("topleft", c("Observed", paste0(model_name, " Direct")), text.col = c("black", "red"), 
   		bty = "n")
   } else {
   	legend("topleft", c("Observed", "NULL", paste0(model_name, " Direct")), text.col = c("black", 
   		"darkgrey", "red"), bty = "n")
   }
   legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
   axis(1, at = tps, labels = FALSE)
   axis(1, at = tps[ind], labels = dates, las = las)

   dev.off()

   return(err = 0)
  }


fit.user.data <-
  function(mydata = NULL,
           par_names = NULL,
           opt.list = NULL,
           run.list = NULL,
           iseed = NULL) {

    #' Fit User Provided Incidence
    #'
    #' \code{fit.user.data} fits and forecasts user provided incidence using an MCMC procedure
    #' and a compartmental S-I-R orr S-E-I-R model
    #'
    #' @param mydata The \pkg{DICE} data list
    #' @param par_names Array with the \pkg{DICE} parameter names
    #' @param opt.list Array with TRUE/FALSE for each f the parameters
    #' @param run.list A list with MCMC parameters
    #' @param iseed Integer - seed for RNG. Default is NULL 
    #'
    #' @return results A list with all the MCMC fitting/forecasting results and the input mydata list
    #' @export
    #'
    #' @examples
    #' results <- fit.user.data(mydata = mydata, par_names = par_names, 
    #' opt.list = opt.list, run.list = run.list, iseed = 12345)	
	tps = cumsum(mydata$ndays)
	
	device = run.list$device
	subDir = run.list$subDir
	
	nperiods = mydata$nperiods
	nperiodsFit = mydata$nperiodsFit
	
	prior = mydata$prior
	Temp = mydata$Temp
	
	# if iseed is not specified, create it. Else generate a unique, reproducible seed for each
	# realization.
	ireal = 1
	if (is.null(iseed)) {
		iseed = set.iseed() * ireal%%.Machine$integer.max
	} else {
		iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
	}
	# Here we set the random number generation seed for R
	set.seed(iseed)
	
	setup <- set.user.data.model(par_names = par_names, mydata = mydata, opt.list = opt.list, run.list = run.list)
	
	par_names = setup$par_names
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
	
	imask = set.imask(par_names = par_names, opt.list = opt.list)
	
	# The number of paramters that are optimized
	nopt = length(which(imask == 1))
	
	pois = 0
	
	cases = mydata$model$epi
	gamaepi = mydata$model$gamaepi
	wght = mydata$model$wght
	
	nRnd = 1000
	model_profile = array(0, c(nRnd, nperiods))
	
	## Need nperiods + 2 because we patch the data at the beginning and end
	
	ndays = tps[nperiods] + mydata$ndays[1] + mydata$ndays[nperiods]
	
	
	## MCMC Fitting procedure 
	##
	if (mydata$epi_model == 1) {
		cat("\n\n Fitting S-I-R model to User Data \n\n")
		out <- .Fortran("fitsir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), 
			parmax = as.double(model_pmax), step = as.double(model_dx), ilog = as.integer(logvec), 
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd), 
			tab = as.single(tab), profile = as.single(model_profile))
	
	} else if (mydata$epi_model == 2) {
		cat("\n\n Fitting S-E-I-R model to User Data \n\n")
		out <- .Fortran("fitseir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), 
			parmax = as.double(model_pmax), step = as.double(model_dx), ilog = as.integer(logvec), 
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd), 
			tab = as.single(tab), profile = as.single(model_profile))
	
	} else {
		cat("\n\n Fitting S-I-R model to User Data \n\n")
		out <- .Fortran("fitsir", epi = as.double(cases), gamaepi = as.double(gamaepi), wght = as.double(wght), nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), 
			parmax = as.double(model_pmax), step = as.double(model_dx), ilog = as.integer(logvec), 
			Temp = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), nperiods = as.integer(nperiods), tps = as.double(tps), rtn = as.double(rep(0, nperiods)), pois = as.double(pois), ndays = as.integer(ndays), nRrnd = as.integer(nRnd), 
			tab = as.single(tab), profile = as.single(model_profile))
	}
	
	model_profile = array(out$profile, c(nRnd, nperiods))
	model_rtn = out$rtn
	
	
	# # 	for (i in 1:nperiods) {
	# model_profile[1:nRnd, i] = rpois(n = nRnd, lambda = model_rtn[i])
	# }
	
	tab.model = matrix(out$tab, ncol = (nparam + 1))
	
	mod_name = mydata$model$name
	mod_id = paste0("ABBV_", mydata$mod_level)
	
	
	## Build the results table
	##
	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
	
	##  Observed number of cases
	##
	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)
	
	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables
	
	## Plot the fit/forecast along with the data
	##
	err <- plotUserDataMECH(mydata = mydata, tables.mod = tables.mod, mod_id = mod_id, ireal = ireal, 
		run.list = run.list, idevice = 1)
	
	
	## write the MCMC statistics
	##
	success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, 
		mydata = mydata, imask = imask, ireal = ireal)
	
	
	## save a binary RData file of the input of this run
	##
	run.input = list(run.list = run.list, opt.list = opt.list)
	
	success = saveRData(mydata = mydata, all_years_epi = NULL, run.input = run.input, ireal = ireal)
	
	## Dump all the profiles we have to a file  
	##
	err <- saveProfilesOnePatch(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, 
		ireal = ireal, run.list = run.list)
	
	## Plot the posterior distribution of the parameters
	##
	
	err <- plotMCMC(mydata = mydata, tab.model = tab.model, opt.list = opt.list, run.list = run.list, 
		imask = imask, ireal = ireal, idevice = 1)
	
	## Build and return the results list
	##
	results = list(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model, tables.mod = tables.mod)
	
	return(results)
	    
  }
