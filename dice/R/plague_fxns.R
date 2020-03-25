set.plague.model <- function(mydata = NULL, par_names = NULL, est_bckgrnd = 1, opt.list = NULL, run.list = NULL) {
  #' Setup of Parameters for an MCMC SEIR procedure
  #'
  #' \code{set.plague.model} Creates the arrays for SEIR modeling of the Plague data with min/max values, step size,
  #' and initial guess/default values for all the parameters.
  #' @param Tg Numeric, recovery rate in days
  #' @param par_names A character array with seir parameter names
  #' @param pop Numeric, the population of the state/region
  #' @param est_bckgrnd Numeric, estimated value of background cholera cases
  #' @param run.list A list with parameters needed for the MCMC procedure  
  #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
  #'   supported by \pkg{DICE}. These values are based on the model used for the basic
  #'   reproduction number  
  #' @return A list with min, max, step size and initial guess/default values for the sir parameters
  #' @examples
  #' set.plague.model(mydata = mydata, par_names = par_names, 
  #' est_bckgrnd=1, run.list = run.list, opt.list = opt.list)

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

	if (is.null(par_names)) {
		par_names <- c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd", "deltaR", "aparam", "alpha", "delta", "ts", "dur")
	}
		
	nparam = length(par_names)
		
	parmin = parmax = pardx = par = logvec = rep(0.0, length = nparam)
	names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names
   	par_opt = opt.list

   	nopt = length(par_opt[par_opt == TRUE])
	
	days_per_week = 7
	NH = mydata$model$pop
	R0 = 2.0
	t0 = 1 * days_per_week
	pC = 0.0001
	seed = 1
	sigma = 1/7.
	Tg = mydata$Tg
	
	if (is.null (Tg)) Tg = 5.

   	## Population
	par["NH"] = parmin["NH"] = parmax["NH"] = NH	
  	par['Tg'] = parmin['Tg'] = parmax['Tg'] = Tg 
  	
   ## R0

   parmin["R0"] = 1.1
   parmax["R0"] = 3.0
   par["R0"]   = R0 * runif(1, 0.8, 1.2)

   ## Sigma 
   parmin['sigma'] = 1/14
   parmax['sigma'] = 1/2
   par['sigma']    = sigma * runif(1, 0.8, 1.2)
   
   ## pC 
   parmin['pC'] = 1e-06
   parmax['pC'] = 1e-03
   par['pC']    = pC * runif(1, 0.8, 1.2)
   
   ## Time of first infection
   parmin['t0'] = 1
   parmax['t0'] = 40
   par['t0']    = t0 * runif(1, 0.8, 1.2)
	
   ## Initial Number of cases 
   
   parmin['seed'] = 1
   parmax['seed'] = 1
   par['seed']    = 1
   
   parmin['e_bckgrnd'] = round(est_bckgrnd*1)
   parmax['e_bckgrnd'] = round(est_bckgrnd*2)
   par['e_bckgrnd']    = est_bckgrnd
   
   par['deltaR'] = par['aprarm'] = par['alpha'] = par['ts'] = par['dur'] = 0
   parmin['deltaR'] = parmin['aprarm'] = parmin['alpha'] = parmin['ts'] = parmin['dur'] = 0
   parmax['deltaR'] = parmax['aprarm'] = parmax['alpha'] = parmax['ts'] = parmax['dur'] = 0
      
	##
	dx = 0.01
	pardx = c(NH = dx, Tg = dx, R0 = dx, sigma = dx, pC = dx, t0 = dx, seed = dx, e_bckgrnd = dx, deltaR = dx, aparam = dx, alpha = dx, ts = dx, dur = dx)
	
	
   	logbase = 10  #use log base 10 when needed
	logvec[1:nparam] = 1
    logvec['t0'] = 0  #linear for t0
    
    nlines = run.list$nlines
    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))


    setup = list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)
        
    return(setup)
		
	
}


fitPlague <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine - Fitting a Single Plague Incidence Profile
    #'
    #' A driver for fitting a single cholera incidence profile of a region/patch. The R code calls a Fortran routine which
    #' uses an MCMC procedure to maximize the likelihood of the solution using and seir mode.
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
    #' fitPlague{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}


	if (mydata$cadence == "Weekly") {
		cadence = 7
		tps = mydata$weeks
	} else if (mydata$cadence == "Monthly") {
		cadence = 31
		tps = mydata$months
	} else if (mydata$cadence == "Daily") {
		cadence = 1
		tps = mydata$days
		ndays = length(tps)
		tps  = 1:ndays
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

	##  Observed number of cases

	epi.obsrv = mydata$model$raw

	cases.model = mydata$model$epi

	est_bckgrnd = 1
	
	par_names <- set.param.list(epi_model = mydata$epi_model)	
	
	setup  <- set.plague.model(par_names = par_names, mydata = mydata, est_bckgrnd = est_bckgrnd, opt.list = opt.list, run.list = run.list)
	
	
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

	imask = set.imask(par_names = par_names, opt.list = opt.list)

	# The number of paramters that are optimized
	nopt = length(which(imask == 1))


	# Here loop on each region and call a single-region fitting routine
	tab.model = NULL

	pois = 0
	nblock = 3
	accept.vec = array(0, c(n, nblock))

	cases = mydata$model$epi
	gamaepi = mydata$model$gamaepi

	sh = mydata$model$sh
	school = mydata$model$school
	
	wght = rep(0, nperiods)
	wght[1:nperiodsFit] = 1.
	mydata$model$wght = wght

	nRnd = 1000
	model_profile = array(0, c(nRnd, nperiods))
	
	## Need nperiods + 2 because we pathc the weeks at the beginning and end

	nmydata = length(tps)

	time = tps
	
    ndays = length(tps) + cadence * 2
	
	cat("\n ****** Fitting ", mydata$model$name, " Using an MCMC Procedure ****** \n")
	
    solution = .Fortran("fitplague", epi = as.double(cases), gamaepi = as.double(gamaepi),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), scales = as.double(Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(time),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd))


     model_profile = array(solution$profile, c(nRnd, nperiods))
     model_rtn = solution$rtn

    tab.model =  matrix(solution$tab, ncol = (nparam + 1))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)
    

   tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

	##  Observed number of cases

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)

 
	model_factor = mydata$model$factor
	model_profile_ili = model_profile/model_factor
	model_rtn_ili = model_rtn/model_factor

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile_ili, nfit = nperiodsFit, state_id = mod_id)
	tables.mod = tables
	
	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, 
    		ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)
    		
 	err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)   		
	
    # # write the MCMC statistics

    success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, mydata = mydata, imask = imask, ireal = ireal)


    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)
    ## now dump all the profiles we have to a file    - in the case of the CDC this is done in the ploting routine
    ## 
    	err <- saveProfilesOnePatch( mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list)
    		
	
	##
	## Plot the posterior distribution of the parameters
	##
	
    err <- plotMCMC(mydata = mydata, tab.model = tab.model, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)


    results = list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model)

    return(results)

}
