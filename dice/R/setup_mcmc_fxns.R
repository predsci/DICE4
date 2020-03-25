##
## Functions that set up parameters for the Fortran MCMC routines (including for the ODEs)
##

setup.model.mcmc <- function(mydata = NULL, run.list = NULL, opt.list = NULL, tps = NULL, par_names = NULL) {

    #' Setup Parameters for MCMC Procedure - Model Data
    #'
    #' Setup all the required parameters for the MCMC procedure on the model data.
    #'   These include: The min/max and step size values for all the parameters,
    #'   the initial values for all the parameters, a vector with the list of parameters
    #'   that will be optimized, the total number of parameters and the number of
    #'   parameters that will be optimized.  The code also allocates an array, tab,
    #'   where the history of the MCMC chain is recorded.
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    #'   supported by \pkg{DICE}. These values are based on the model used for the basic
    #'   reproduction number.
    #' @param tps  A numeric array of days for the disease season - the day numbers are
    #'  consistent with the weeks/months
    #' @param par_names - An array with the all parameters ordered as required by DICE
    #' @return
    #' A list with the following  arguments
    #' \describe{
    #' \item{par_min}{Minimum values for all the parameters supported by \pkg{DICE}}
    #' \item{par_max}{Maximum values for all the parameters supported by \pkg{DICE}}
    #' \item{pardx}{Step-size for MCMC for all the parameters supported by \pkg{DICE}}
    #' \item{par}{Initial values for all the parameters supported by \pkg{DICE}}
    #' \item{logbase}{Base for log values - currently code assumes base 10}
    #' \item{logvec}{Array of integers with 1, use log base, or 0 - do not}
    #' \item{nparam}{Integer-the total number of parameters recognized by the \pkg{DICE} code}
    #' \item{nopt}{Integer-the number of parameters that will be optimized}
    #' }
    #' @examples
    #' setup.model.mcmc{mydata = mydata,
    #' run.list = run.list, opt.list = opt.list, tps = tps}

   dataType = mydata$data_source
   cases = mydata$model$epi
   nperiods = mydata$nperiods
   cadence = mydata$cadence
   model = mydata$imodel
   disease = tolower(mydata$disease)
	
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


   pop = mydata$model$pop


   nparam = length(par_names)

   parmin = parmax = pardx = par = logvec = rep(0.0, length = nparam)

   names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names

   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

   ## Population
   par["NH"] = parmin["NH"] = parmax["NH"] = pop

   ## Tg - in days
   Tg = mydata$Tg
   if (is.null(mydata$Tg)) {
   	Tg = 3
   	if (tolower(mydata$disease) == "dengue")
   		Tg = 10
	if (tolower(mydata$disease) == 'ebola'){
		if (mydata$epi_model == 1) Tg = 12.0
		if (mydata$epi_model == 1) Tg = 10.0		
	}
		
   }

   par["Tg"] = parmin["Tg"] = parmax["Tg"] = Tg

   ## R0

   parmin["R0"] = 1.0
   parmax["R0"] = 3
   pardx["R0"] = 0.05

   if (opt.list$R0) {
   	R0.rnd = runif(1, 1.1, 1.4)
   } else {
   	R0.rnd = 1.2
   }

   if(disease == 'ebola') {
   	parmax["R0"] = 4.0
   	R0.rnd = runif(1, 1.1, 2.)
   	 par["R0"] = R0.rnd
   	}
   	
   if (disease == "dengue") {
   	parmin["R0"] = 1.1
   	parmax["R0"] = 8
   	pardx["R0"] = 0.05
   	R0.rnd = 3
   	if (opt.list$R0)
   		R0.rnd = runif(1, 2, 3)
   	par["R0"] = R0.rnd
   }

   
   parmin["sigma"] = 1/5.
   parmax["sigma"] = 1/1.
   pardx["sigma"] = 0.05

	if (disease == "ebola") {
		parmin["sigma"] = 1/10
		parmax["sigma"] = 1/2
		pardx["sigma"] = 0.05
	}
	
   if (opt.list$sigma) {
   	sigma.rnd = runif(1, parmin["sigma"], parmax["sigma"])
   } else {
   	sigma.rnd = mean(parmin["sigma"] + parmax["sigma"])
   }

   par['sigma'] = sigma.rnd

   parmin["pC"] = 1e-05
   parmax["pC"] = 1
   pardx["pC"] = 0.05

   if (opt.list$pC) {
   	pC.rnd = runif(1, 0.02, 0.1)
   } else {
   	pC.rnd = 0.05
   }

   par['pC'] = pC.rnd

	## Time of first case

	t0 = tps[1]
	parmin["t0"] = tps[1]
	parmax["t0"] = tps[(nperiods/2)]
	pardx["t0"] = 0.05
	par["t0"] = t0 + round(runif(1, min = 1 * cadence, max = 4 * cadence))

	# Initial number of cases - seed
	parmin["seed"] = 1
	parmax["seed"] = 1000
	pardx["seed"] = 0.05

	if (opt.list$seed) {
		seed.rnd = runif(1, min = 1, max = 10)
		seed.rnd = round(seed.rnd)
	} else {
		seed.rnd = 1
	}

	par["seed"] = seed.rnd

	## background cases
	Baseline = mean(cases[1:5])
	if (dataType == "dengue")
		Baseline = mean(cases[1:2]) #, cases[(nperiods - 4):nperiods]))
	Baseline = max(Baseline, 1)
	parmin["e_bckgrnd"] = 0.1 * Baseline
	parmax["e_bckgrnd"] = 2 * Baseline
	pardx["e_bckgrnd"] = 0.05

    if (opt.list$e_bckgrnd) {
   	e_bckgrnd.rnd = runif(1, parmin["e_bckgrnd"], parmax["e_bckgrnd"])
   } else {
   	e_bckgrnd.rnd = 1
   }

   par["e_bckgrnd"] = e_bckgrnd.rnd

   ## for SH Term

   if (any(model == c(1, 3))) {
   	deltaR = runif(1, 0.1, 0.3)
   	aparam = runif(1, 100, 200)
   } else {
   	deltaR = 0
   	aparam = 0
   }

    par['deltaR'] = deltaR
    par['aparam'] = aparam

    parmin['deltaR'] = 1e-06
    parmin['aparam'] = 1

    parmax['deltaR'] = 4
    parmax['aparam'] = 1000

    pardx['deltaR'] = 0.05
    pardx['aparam'] = 0.05

    ## For school term
    if (any(model == c(2, 3))) {
        alpha = runif(1, 0.05, 0.2)
	} else {
        alpha = 0
    }

    par["alpha"] = alpha
    parmin["alpha"] = 1e-06
    parmax["alpha"] = 1
    pardx["alpha"] = 0.05

    ## Two-value model: delta, ts, dur
    parmin["delta"] = -1
    parmin["ts"] = t0
    parmin["dur"] = 0.01

    parmax["delta"] = 1
    parmax["ts"] = t0 + 30 * cadence
    parmax["dur"] = 30 * cadence

    pardx["delta"] = 0.05
    pardx["ts"] = 0.05
    pardx["dur"] = 0.05

    if (model == 5) {
    	delta = runif(1, min = -0.2, max = 0.2)
    	ts = runif(1, min = (t0 + cadence * 10), max = (t0 + cadence * 20))
    	dur = runif(1, min = (cadence * 2), max = (cadence * 10))
    } else {
    	delta = 0
    	ts = 0
    	dur = 0
    }

    par['delta'] = delta
    par['ts'] = ts
    par['dur'] = dur

    logbase = 10  #use log base 10 when needed
	logvec[1:nparam] = 1
    logvec['t0'] = 0  #linear for t0
    logvec['delta'] = 0  # and linear for delta since it can be negative

    nlines = run.list$nlines
    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))

    list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)

}

setup.fit.mcmc <- function(mydata = NULL, run.list = NULL, opt.list = NULL, tps = NULL, par_names = par_names) {

    #' Setup parameters for MCMC Procedure - Fit Data
    #'
    #' Setup all the required parameters for the MCMC procedure on the fit data.
    #'   These include: The min/max and step size values for all the parameters,
    #'   the initial values for all the parameters, a vector with the list of parameters
    #'   that will be optimized, the total number of parameters and the number of
    #'   parameters that will be optimized.  The code also allocates an array, tab,
    #'   where the history of the MCMC chain is recorded.
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    #'   and find the proper values needed for this prior.
    #'   supported by \pkg{DICE}. These values are based on the model used for the basic
    #'   reproduction number.
    #' @param tps  A numeric array of days for the disease season - the day numbers are
    #'   consistent with the weeks/months
    #' @param par_names - - An array with the all parameters ordered as required by DICE
    #' @return
    #' A list with the following  arguments:
    #' \describe{
    #' \item{par_min}{Minimum values for all the parameters supported by \pkg{DICE} for each region}
    #' \item{par_max}{Maximum values for all the parameters supported by \pkg{DICE} for each region}
    #' \item{pardx}{Step-size for MCMC for all the parameters supported by \pkg{DICE} for each region}
    #' \item{par}{Initial values for all the parameters for each region}
    #' \item{logbase}{Base for log values - currently code assumes base 10}
    #' \item{nopt}{Integer-the number of parameters that will be optimized}
    #' \item{tab}{A 2D numeric array with nlines and (nparam+1) columns used to store the MCMC history
    #'   of all the parameters and the likelihood}
    #' }
    #' @examples
    #' setup.fit.mcmc{mydata = mydata,
    #' run.list = run.list, opt.list = opt.list, tps = tps, par_names = par_names}

    dataType = mydata$data_source
    nperiods = mydata$nperiods
    cases  = mydata$fit$epi
    cadence = mydata$cadence
    n = dim(cases)[2]
    model = mydata$imodel
	disease = tolower(mydata$disease)
	
	
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
	
    pop = mydata$fit$pop

    nparam = length(par_names)

    par_opt = opt.list

   	nopt = length(par_opt[par_opt == TRUE])

    nlines = run.list$nlines
    t0 = tps[1]

    parmin = parmax = pardx = par = array(0, dim = c(nparam, n))

	rownames(parmin) = rownames(parmax) = rownames(pardx) = rownames(par) = par_names

	logvec = rep(0, nparam)

    names(logvec) = par_names

    # Population

   par["NH", 1:n] = parmin['NH', 1:n] = parmax['NH', 1:n] = pop

	Tg = mydata$Tg
   if (is.null(mydata$Tg)) {
   	Tg = 3
   	if (disease == "dengue")
   		Tg = 10
	if (tolower(mydata$disease) == 'ebola'){
		if (mydata$epi_model == 1) Tg = 12.0
		if (mydata$epi_model == 1) Tg = 10.0		
	}   	
   }

    par['Tg', 1:n] = parmin['Tg',1:n] = parmax['Tg', 1:n] = Tg

    ## R0

   parmin["R0", 1:n] = 1.0
   parmax["R0", 1:n] = 3
   pardx["R0", 1:n] = 0.05
   

   if (opt.list$R0) {
   	R0.rnd = runif(n, 1.1, 1.4)
   } else {
   	R0.rnd = rep(1.2, n)
   }

   if(disease == 'ebola') {
   	parmax["R0", 1:n] = 4.0
   	R0.rnd = runif(n, 1.1, 2.)
   	}
   	
  if (disease == "dengue") {
   	parmin["R0",1:n] = 1.1
   	parmax["R0",1:n] = 8
   	pardx["R0",1:n] = 0.05
   	R0.rnd = rep(3, n)
   	if (opt.list$R0)
   		R0.rnd = runif(n, 2, 3)
   }

    par["R0", 1:n] = R0.rnd

   parmin["sigma", 1:n] = 1./5.
   parmax["sigma",1:n] = 1./1.
   pardx["sigma",1:n] = 0.05

	if (disease == "ebola") {
		parmin["sigma", 1:n] = 1/10
		parmax["sigma", 1:n] = 1/2
		pardx["sigma", 1:n]  = 0.05
	}
	
   if (opt.list$sigma) {
   	sigma.rnd = runif(n, parmin["sigma", 1], parmax["sigma", 1])
   } else {
   	sigma.rnd = rep(mean(parmin["sigma", 1] + parmax["sigma", 1]), n)
   }

   par['sigma', 1:n] = sigma.rnd

   parmin["pC", 1:n] = 1e-05
   parmax["pC", 1:n] = 1
   pardx["pC", 1:n] = 0.05

   if (opt.list$pC) {
   	pC.rnd = runif(n, 0.02, 0.1)
   } else {
   	pC.rnd = rep(0.05, n)
   }

   par['pC',1:n] = pC.rnd

	## Time of first case

	t0 = tps[1]
	parmin["t0", 1:n] = tps[1]
	parmax["t0", 1:n] = tps[(nperiods/2)]
	pardx["t0", 1:n] = 0.05
	par["t0", 1:n] = t0 + round(runif(n, min = 1 * cadence, max = 4 * cadence))

	# Initial number of cases - seed
	parmin["seed", 1:n] = 1
	parmax["seed", 1:n] = 1000
	pardx[ "seed", 1:n] = 0.05
	if(disease == 'ebola') parmax["seed", 1:n] = 100
	
	if (opt.list$seed) {
		seed.rnd = runif(n, min = 1, max = 10)
		seed.rnd = round(seed.rnd)
	} else {
		seed.rnd = rep(1, n)
	}

	par["seed", 1:n] = seed.rnd

	## background cases
	Baseline = e_bckgrnd.rnd = rep(0, n)

	for(i in 1:n) Baseline[i] = mean(cases[1:5, i])
	if (dataType == "dengue")
		for(i in 1:n) Baseline[i] = mean(cases[1:2, i]) #, cases[(nperiods - 4):nperiods]))
	for (i in 1:n) Baseline[i] = max(Baseline[i], 1)
	parmin["e_bckgrnd", 1:n] = 0.1 * Baseline
	parmax["e_bckgrnd", 1:n] = 2 * Baseline
	pardx["e_bckgrnd", 1:n] = 0.05

    if (opt.list$e_bckgrnd) {
    	for (i in 1:n ) e_bckgrnd.rnd[i] = runif(1, parmin["e_bckgrnd", i], parmax["e_bckgrnd", i])
    } else {
    	e_bckgrnd.rnd = rep(1, n)
    }

   par["e_bckgrnd", 1:n] = e_bckgrnd.rnd

    # now add deltaR, aparam and alpha for SH (first two) and School (third)

    if (model == 1 | model == 3) {
        deltaR = runif(n, 0.05, 0.1)
        aparam = runif(n, 100, 200)
     } else {
        deltaR = rep(0, n)
        aparam = rep(0, n)
    }

    par['deltaR', ] = deltaR
    par['aparam', ] = aparam
    parmin['deltaR', ] = rep(1e-06, n)
    parmin['aparam', ] = rep(1, n)

    parmax['deltaR', ] = rep(4, n)
    parmax['aparam', ] = rep(1000, n)
    pardx['deltaR', ] = rep(0.05, n)
    pardx['aparam', ] = rep(0.05, n)

    if (model == 2 | model == 3) {
        alpha = runif(n, 0.05, 0.1)
     } else {
        alpha = rep(0, n)
    }

    par['alpha', 1:n] = alpha
    parmin['alpha', 1:n] = rep(1e-06, n)
    parmax['alpha', 1:n] = rep(1, n)
    pardx['alpha', 1:n] = rep(0.05, n)

    parmin['delta', 1:n] = rep(-1, n)
    parmin['ts', 1:n] = rep(t0, n)
    parmin['dur', 1:n] = rep((0.01), n)
    parmax['delta', 1:n] = rep(1, n)
    parmax['ts', 1:n] = rep((t0 + 30 * cadence), n)
    parmax['dur', 1:n] = rep((30 * cadence), n)

    pardx['delta', 1:n] = rep(0.05, n)
    pardx['ts', 1:n] = rep(0.05, n)
    pardx['dur', 1:n] = rep(0.05, n)

    if (model == 5) {
        delta = runif(n, min = -0.2, max = 0.2)
        ts = runif(n, min = (t0 + cadence * 10), max = (t0 + cadence * 20))
        dur = runif(n, min = (cadence * 2), max = (cadence * 10))
    } else {
        delta = rep(0, n)
        ts = rep(0, n)
        dur = rep(0, n)
    }

    par['delta', 1:n] = delta
    par['ts'   , 1:n] = ts
    par['dur'  , 1:n] = dur


    logbase = 10  #use log base 10 when needed
    logvec[1:nparam] # use log uniform for everything
    logvec['t0'] = 0  #linear for t0
    logvec['delta'] = 0  # and linear for delta since it can be negative

    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))

    list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)

}

setup.coupling.mcmc <- function(par_names_cpl = NULL, opt.list = NULL) {

    #' Setup For Coupling Parameters
    #'
    #' For a coupled \pkg{DICE} run we setup all the required values for
    #' the two parameters that help define the coupling matrix.  These
    #'   include min/max and step size values, initial values, logscale
    #' and number of parameters. The two parameters are: the saturation distance
    #'   \eqn{s_d} in km, and the distance power \eqn{\gamma}.
    #' @param par_names_cpl An array with the names of the coupling parameters
    #' @param opt.list a logical list with TRUE or FALSE values for each
    #' @return coupling.list A named list with the following arguments:
    #' \describe{
    #' \item{parmin}{Minimum values for the two parameters}
    #' \item{parmax}{Maximum values for the two parameters}
    #' \item{pardx}{Step-size for the two parameters}
    #' \item{par}{Initial values for the two parameters}
    #' \item{logbase}{Base for log values - currently code assumes base 10}
    #' \item{logvec}{Array of integers with 1, use log base, or 0 - do not}
    #' \item{nparam}{Integer-the total number of parameters recognized by the \pkg{DICE} code}
    #' \item{nopt}{Integer-the number of parameters that will be optimized}
    #' }
    #' @examples
    #' setup.coupling.mcmc{par_names_cpl = par_names_cpl, opt.list=opt.cpl}



   nparam = length(par_names_cpl)

   parmin = parmax = pardx = par = logvec = rep(0.0, length = nparam)

   names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names_cpl

   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

    # # Now populate the saturation distance sd (in km)
    parmin['sd'] = 5
    parmax['sd'] = 500
    pardx['sd'] = 0.01

    if (opt.list$sd) {
        sd.rnd = runif(1, 100, 300)
    } else {
        sd.rnd = 200
    }
    par['sd'] = sd.rnd

    parmin['gamma'] = 1
    parmax['gamma'] = 6
    pardx['gamma'] = 0.01


    if (opt.list$gamma) {
        gamma.rnd = runif(1, 2, 5)
    } else {
        gamma.rnd = 4
    }

    par['gamma'] = gamma.rnd

    nparam = length(par)

    logbase = 10  #use log base 10 when needed
    logvec[1:nparam] = 1 # use log uniform for everything

    coupling.list = list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec,
        nparam = nparam, nopt = nopt)

    return(coupling.list)
}


set.imask <- function(par_names = NULL, opt.list = NULL) {

    #' Set a Mask for Parameters
    #'
    #' Given a named logical list of all the parameters that \pkg{DICE} supports
    #'   prepare an integer list of the same length with +1/-1 for parameters
    #'   that are/are not optimized
    #' @param par_names - A list with all the \pkg{DICE} parameter names
    #' @param opt.list A named logical list with TRUE/FALSE values for each parameter
    #' @return imask An integer list of length nparam with +1 or -1 values
    #' @examples
    #' set.imask{par_names = par_names, opt.list = opt.list}

    nparam = length(par_names)

    imask = rep(-1, nparam)
    names(imask) = par_names

    for (i in 1:nparam) {
    	if(opt.list[i] == TRUE) imask[i] = 1
    }

    return(imask)
}

set.cpl.imask <- function(nparam = NULL, opt.list = NULL) {

    #' Set a Mask for the Coupling Parameters
    #'
    #' Set +1 or -1 values for the two parameter that help define the
    #' coupling matrix.    These are needed and called for only in the case
    #' of a coupled \pkg{DICE} run.
    #' @param nparam Integer - The number of parameters that help define the coupling matrix.   Currently two.
    #' @param opt.list A named logical list with TRUE or FALSE values for the coupling parameters.
    #' @return imask An integer list of length nparam with +1 or -1 values
    #' examples
    #' set.cpl.imask{nparam=cpl_nparam,opt.list=opt.cpl}

    imask = rep(1, nparam)

    imask[1:nparam] = -1
    if (opt.list$sd)
        imask[1] = 1  #optimize the saturation distance
    if (opt.list$gamma)
        imask[2] = 1  #optimize the distance power

    return(imask)
}

setup.vector.model.mcmc <- function(mydata = NULL, opt.list = NULL, run.list = NULL, tps = NULL, par_names = NULL) {


    #' Setup Parameters for MCMC Procedure - Model Dengue Data
    #'
    #' Setup all the required parameters for the MCMC procedure on the model data.
    #'   These include: The min/max and step size values for all the parameters,
    #'   the initial values for all the parameters, a vector with the list of parameters
    #'   that will be optimized, the total number of parameters and the number of
    #'   parameters that will be optimized.  The code also allocates an array, tab,
    #'   where the history of the MCMC chain is recorded.
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    #'   supported by \pkg{DICE}. These values are based on the model used for the basic
    #'   reproduction number.
    #' @param tps  A numeric array of days for the disease season - the day numbers are
    #'  consistent with the weeks/months
    #' @param par_names - An array with the all parameters ordered as required by DICE
    #' @return A list with the following  arguments
    #' \describe{
    #' \item{par_min}{Minimum values for all the parameters supported by \pkg{DICE}}
    #' \item{par_max}{Maximum values for all the parameters supported by \pkg{DICE}}
    #' \item{pardx}{Step-size for MCMC for all the parameters supported by \pkg{DICE}}
    #' \item{par}{Initial values for all the parameters supported by \pkg{DICE}}
    #' \item{logbase}{Base for log values - currently code assumes base 10}
    #' \item{logvec}{Array of integers with 1, use log base, or 0 - do not}
    #' \item{nparam}{Integer-the total number of parameters recognized by the \pkg{DICE} code}
    #' \item{nopt}{Integer-the number of parameters that will be optimized}
    #' }
    #' @examples
    #' setup.model.mcmc{mydata = mydata,
    #' run.list = run.list, opt.list = opt.list, tps = tps, par_names = par_names}

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
	
    if (is.null(par_names))
        par_names <- set.param.list(epi_model = mydata$epi_model)

    nparam = length(par_names)

   parmin = parmax = pardx = par = logvec = rep(0.0, length = nparam)

   names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names

   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

   cases = mydata$model$epi

   ## Population
   par["NH"] = parmin["NH"] = parmax["NH"] = mydata$model$pop

   ## Tg - in days
   Tg = mydata$Tg
   if (is.null(mydata$Tg)) {
   	Tg = 3
   	if (tolower(mydata$disease) == "dengue") {Tg = 8. }
   	if (tolower(mydata$disease) == "yellow_fever") {Tg = 5.}
   	if (tolower(mydata$disease) == "chik") {Tg = 9.}
   	if (tolower(mydata$disease) == "zika") {Tg = 7.}
   }
	
   par["Tg"] = parmin["Tg"] = parmax["Tg"] = Tg

   ## R0

   parmin["R0"] =  1.1
   parmax["R0"] = 4
   if (tolower(mydata$disease) == "yellow_fever" || tolower(mydata$disease) == 'zika') {parmax["R0"] = 12.}
   if (tolower(mydata$disease) == "chik")		  {parmin["R0"] = 1.5} 			
   pardx["R0"] = 0.05
   R0.rnd = 3
   if (opt.list$R0)
   	R0.rnd = runif(1, 2, 3)

   par["R0"] = R0.rnd

   parmin["sigma"] = 1/7
   parmax["sigma"] = 1/4
   if (tolower(mydata$disease) == "yellow_fever") {
    parmin["sigma"] = 1/8
   	parmax["sigma"] = 1/2 	 	
   }
   if (tolower(mydata$disease) == "chik") {
    parmin["sigma"] = 1/8
   	parmax["sigma"] = 1/2 	 	
   }
   if (tolower(mydata$disease) == "zika") {
    parmin["sigma"] = 1/10
   	parmax["sigma"] = 1/1	 	
   }  
      
   pardx["sigma"] = 0.05

   if (opt.list$sigma) {
   	sigma.rnd = runif(1, parmin["sigma"], parmax["sigma"])
   } else {
   	sigma.rnd = mean(parmin["sigma"] + parmax["sigma"])
   }

   par["sigma"] = sigma.rnd


   parmin["pC"] = 1e-06
   parmax["pC"] = 1   
   pardx["pC"] = 0.05

   if (opt.list$pC) {
   	pC.rnd = runif(1, 0.001, 0.005)
   } else {
   	pC.rnd = 0.01
   }

   par['pC'] = pC.rnd

	if(mydata$mod_level == 2 & mydata$model$attr$ABBV_2 == 'LK') {
		parmax['pC'] = 0.01
	}
	day_per_week = 7

	parmin["t0"] = 0.1
	parmax["t0"] =  200.
	if (tolower(mydata$disease) == "yellow_fever") {
		parmax["t0"] =  90
	}
	if (tolower(mydata$disease) == "chik") {
		t0 = tps[1]
		parmin["t0"] = tps[1]
		parmax["t0"] = tps[(nperiods/4)]
		pardx["t0"] = 0.05
		par["t0"] = t0 + round(runif(1, min = 1 * cadence, max = 4 * cadence))
	}
		
	pardx["t0"] = 0.05
	par["t0"] = round(runif(1, min = 1, max = day_per_week))

	# Initial number of cases - seed
	parmin["seed"] = 1
	parmax["seed"] = 100

   if (tolower(mydata$disease) == "yellow_fever") {
   		parmin["seed"] = 50
   		parmax["seed"] = 1000
   }
   if (tolower(mydata$disease) == "chik") {
   		parmin["seed"] = 1
   		parmax["seed"] = 100
   }   		
	pardx["seed"] = 0.05

	if (opt.list$seed) {
		seed.rnd = runif(1, min = 2, max = parmax["seed"])
		seed.rnd = round(seed.rnd)
	} else {
		seed.rnd = 10
	}

	par["seed"] = seed.rnd

	## background cases
	Baseline = mean(cases[1:2])
	Baseline = max(Baseline, 1)
	parmin["e_bckgrnd"] = 0.1 * Baseline
	parmin["e_bckgrnd"] = max(1,parmin["e_bckgrnd"])
	parmax["e_bckgrnd"] = 2 * Baseline
	parmin["e_bckgrnd"] = round(parmin["e_bckgrnd"])
	parmax["e_bckgrnd"] = round(parmax["e_bckgrnd"])
   if (tolower(mydata$disease) == "yellow_fever") {
   		parmin["e_bckgrnd"] = 1
   		parmax["e_bckgrnd"] = 5
   }	
	pardx["e_bckgrnd"] = 0.05

    if (opt.list$e_bckgrnd) {
   	e_bckgrnd.rnd = runif(1, parmin["e_bckgrnd"], parmax["e_bckgrnd"])
   } else {
   	e_bckgrnd.rnd = 1
   }

   par["e_bckgrnd"] = e_bckgrnd.rnd


  ## for SH Term

   if (any(mydata$imodel == c(1, 3))) {
   	deltaR = runif(1, 0.1, 0.3)
   	aparam = runif(1, 100, 200)
   } else {
   	deltaR = 0
   	aparam = 0
   }

    par['deltaR'] = deltaR
    par['aparam'] = aparam

    parmin['deltaR'] = 1e-06
    parmin['aparam'] = 1

    parmax['deltaR'] = 4
    parmax['aparam'] = 1000

    pardx['deltaR'] = 0.05
    pardx['aparam'] = 0.05

    ## For school term
    if (any(mydata$imodel == c(2, 3))) {
        alpha = runif(1, 0.05, 0.2)
	} else {
        alpha = 0
    }

    par["alpha"] = alpha
    parmin["alpha"] = 1e-06
    parmax["alpha"] = 1
    pardx["alpha"] = 0.05

    ## Two-value model: delta, ts, dur
    parmin["delta"] = -1
    parmin["ts"] = tps[1]
    parmin["dur"] = 0.01

    parmax["delta"] = 1
    parmax["ts"] = tps[1] + 10 * day_per_week
    parmax["dur"] = 30 * day_per_week

    pardx["delta"] = 0.05
    pardx["ts"] = 0.05
    pardx["dur"] = 0.05

    if (mydata$imodel == 5) {
    	t0 = tps[1]
    	delta = runif(1, min = -0.2, max = 0.2)
    	ts = runif(1, min = (t0 + cadence * 10), max = (t0 + cadence * 20))
    	dur = runif(1, min = (cadence * 2), max = (cadence * 10))
    } else {
    	delta = 0
    	ts = 0
    	dur = 0
    }

    par['delta'] = delta
    par['ts'] = ts
    par['dur'] = dur

	##
	## If Vector States are included explicitly
	##

	if (mydata$epi_model > 2) {
		pardx["bite"] = pardx["vec_k"] = pardx["muV"] = pardx["T_HV"] = pardx["T_VH"] = pardx["sigmaV"] = 0.1

		parmin["bite"] = 0.1
		parmax["bite"] = 1.
		par["bite"]    = 0.5

		if (opt.list$bite)
			par["bite"] = runif(1, parmin["bite"], 1.0)

		parmin["vec_k"] = 0.3
		parmax["vec_k"] = 4
		par["vec_k"] = 2

		if (opt.list$vec_k)
			par["vec_k"] = runif(1, parmin["vec_k"], parmax["vec_k"])

		parmin["muV"] = 1/42
		parmax["muV"] = 1/8
		par["muV"] = 1/25
		if (opt.list$muV)
			par["muV"] = runif(1, parmin["muV"], parmax["muV"])

		parmin["T_HV"] = parmin["T_VH"] = 0.33
		parmax["T_HV"] = parmax["T_VH"] = 1.0
		par["T_HV"] = par["T_VH"] = 0.5
		if (opt.list$T_HV)
			par["T_HV"] = runif(1, parmin["T_HV"], parmax["T_HV"])
		if (opt.list$T_VH)
			par["T_VH"] = runif(1, parmin["T_VH"], parmax["T_VH"])

		parmin["sigmaV"] = 1/14
		parmax["sigmaV"] = 1/7
		par["sigmaV"] = 1/10

		if (opt.list$sigmaV)
			par["sigmaV"] = runif(1, parmin["sigmaV"], parmax["sigmaV"])

	}


    logbase = 10  #use log base 10 when needed
	logvec[1:nparam] = 1
    logvec['t0'] = 0  #linear for t0
    logvec['delta'] = 0  # and linear for delta since it can be negative

    nlines = run.list$nlines
    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))

    list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)


}


setup.vector.fit.mcmc <- function(mydata = NULL, opt.list = NULL, run.list = NULL, tps = NULL, par_names = NULL) {

    #' Setup Parameters for MCMC Procedure - Fit Dengue Data
    #'
    #' Setup all the required parameters for the MCMC procedure on the fit data.
    #'   These include: The min/max and step size values for all the parameters,
    #'   the initial values for all the parameters, a vector with the list of parameters
    #'   that will be optimized, the total number of parameters and the number of
    #'   parameters that will be optimized.  The code also allocates an array, tab,
    #'   where the history of the MCMC chain is recorded.
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
    #'   and find the proper values needed for this prior.
    #'   supported by \pkg{DICE}. These values are based on the model used for the basic
    #'   reproduction number.
    #' @param tps  A numeric array of days for the disease season - the day numbers are
    #'   consistent with the weeks/months
    #' @param par_names - - An array with the all parameters ordered as required by DICE
    #' @return A list with the following  arguments:
    #' \describe{
    #' \item{par_min}{Minimum values for all the parameters supported by \pkg{DICE} for each region}
    #' \item{par_max}{Maximum values for all the parameters supported by \pkg{DICE} for each region}
    #' \item{pardx}{Step-size for MCMC for all the parameters supported by \pkg{DICE} for each region}
    #' \item{par}{Initial values for all the parameters for each region}
    #' \item{logbase}{Base for log values - currently code assumes base 10}
    #' \item{nopt}{Integer-the number of parameters that will be optimized}
    #' \item{tab}{A 2D numeric array with nlines and (nparam+1) columns used to store the MCMC history
    #'   of all the parameters and the likelihood}
    #' }
    #' @examples
    #' setup.fit.mcmc{mydata = mydata,
    #' run.list = run.list, opt.list = opt.list, tps = tps, par_names = par_names}

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
	
    if (is.null(par_names))
        par_names <- set.param.list(epi_model = mydata$epi_model)

    cases = mydata$fit$epi

    n = dim(cases)[2]

	nparam = length(par_names)

    parmin = parmax = pardx = par = array(0, dim = c(nparam, n))

	rownames(parmin) = rownames(parmax) = rownames(pardx) = rownames(par) = par_names


   par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])

   ## Population
   par["NH",1:n] = parmin["NH",1:n] = parmax["NH",1:n] = mydata$fit$pop

   ## Tg - in days
   Tg = mydata$Tg
   if (is.null(mydata$Tg)) {
   	Tg = 3
   	if (tolower(mydata$disease) == "dengue")
   		Tg = 8.
   	if (tolower(mydata$disease) == "yellow_fever")
   		Tg = 5. 
   	if (tolower(mydata$disease) == "chik")
   		Tg = 9.    		  		
   	if (tolower(mydata$disease) == "zika")
   		Tg = 7.   		
   }


   par["Tg",1:n] = parmin["Tg",1:n] = parmax["Tg",1:n] = Tg

   ## R0

   parmin["R0", 1:n] = 1.1
   parmax["R0", 1:n] = 4.0
   if (tolower(mydata$disease) == "yellow_fever" || tolower(mydata$disease) == "zika") {parmax["R0", 1:n] = 12.0}
   if (tolower(mydata$disease) == "chik" ) {parmin["R0", 1:n] = 1.5}
   pardx["R0", 1:n] = 0.05
   R0.rnd = 3.
   if (opt.list$R0)
   	R0.rnd = runif(n, 1.5, 2.)

   par["R0", 1:n] = R0.rnd

   parmin["sigma", 1:n] = 1/7
   parmax["sigma", 1:n] = 1/4
   pardx["sigma" , 1:n] = 0.05

   if (tolower(mydata$disease) == "yellow_fever") {
    parmin["sigma", 1:n] = 1/8
   	parmax["sigma", 1:n] = 1/2 	 	
   }
   if (tolower(mydata$disease) == "chik") {
    parmin["sigma", 1:n] = 1/8
   	parmax["sigma", 1:n] = 1/2 	 	
   } 
   if (tolower(mydata$disease) == "zika") {
    parmin["sigma", 1:n] = 1/10
   	parmax["sigma", 1:n] = 1/1	 	
   }         
   if (opt.list$sigma) {
   	sigma.rnd = runif(n, parmin["sigma",1], parmax["sigma",1])
   } else {
   	sigma.rnd = rep(mean(parmin["sigma", 1] + parmax["sigma", 1]), n)
   }

   par["sigma", 1:n ] = sigma.rnd


   parmin["pC", 1:n] = 1e-06
   parmax["pC", 1:n] = 1
   pardx["pC" , 1:n] = 0.05

   if (opt.list$pC) {
   	pC.rnd = runif(n, 0.0005, 0.002)
   } else {
   	pC.rnd = rep(0.01, n)
   }

   par['pC', 1:n] = pC.rnd

	if(mydata$mod_level == 2 & mydata$model$attr$ABBV_2 == 'LK') {
		parmax['pC', 1:n] = 0.01
	}
	day_per_week = 7

	parmin["t0", 1:n] = 0.1
	parmax["t0", 1:n] =  200.
	pardx["t0", 1:n] = 0.05
   if (tolower(mydata$disease) == "yellow_fever") {
   	parmin["t0", 1:n] = 1.
    parmax["t0", 1:n] = 90.
   }
	par["t0", 1:n] = round(runif(n, min = 1, max = day_per_week))	

	if (tolower(mydata$disease) == "chik") {
		t0 = tps[1]
		parmin["t0", 1:n] = tps[1]
		parmax["t0", 1:n] = tps[(nperiods/4)]
		pardx["t0", 1:n] = 0.05
		par["t0", 1:n] = t0 + round(runif(n, min = 1 * cadence, max = 4 * cadence))
	}
		


	# Initial number of cases - seed
	parmin["seed", 1:n] = 1
	parmax["seed", 1:n] = 100
	if (tolower(mydata$disease) == 'yellow_fever') {
		parmin["seed", 1:n] = 50.0
		parmax["seed", 1:n] = 1000.0
	}

   if (tolower(mydata$disease) == "chik") {
   		parmin["seed", 1:n] = 1
   		parmax["seed", 1:n] = 10
   }   		
   			
	pardx["seed" , 1:n] = 0.05

	if (opt.list$seed) {
		seed.rnd = runif(n, min = 2, max = parmax["seed", 1])
		seed.rnd = round(seed.rnd)
	} else {
		seed.rnd = rep(10, n)
	}

	par["seed", 1:n ] = seed.rnd

	## background cases
	Baseline = rep(0, n)
	for (i in 1:n) {
		Baseline[i] = mean(cases[1:2, i])
		Baseline[i] = max(Baseline[i], 1)
	}

	parmin["e_bckgrnd", 1:n] = 0.1 * Baseline
	parmax["e_bckgrnd", 1:n] = 2.0 * Baseline
	for (i in 1:n ) parmin["e_bckgrnd", i] = max(1, parmin["e_bckgrnd", i])
	parmin["e_bckgrnd", 1:n] = round(parmin["e_bckgrnd", 1:n])
	parmax["e_bckgrnd", 1:n] = round(parmax["e_bckgrnd", 1:n])	
	
	if (tolower(mydata$disease) == 'yellow_fever') {
		parmin["e_bckgrnd", 1:n] = 1.
		parmax["e_bckgrnd", 1:n] = 5.
	}	
	pardx["e_bckgrnd" , 1:n] = 0.05

	e_bckgrnd.rnd = rep(0, n)
    if (opt.list$e_bckgrnd) {
    	for (i in 1:n )
   		e_bckgrnd.rnd[i] = runif(1, parmin["e_bckgrnd", i], parmax["e_bckgrnd", i])
   } else {
   	e_bckgrnd.rnd = rep(1, n)
   }

   par["e_bckgrnd", 1:n] = e_bckgrnd.rnd


  ## for SH Term

   if (any(mydata$imodel == c(1, 3))) {
   	deltaR = runif(n, 0.1, 0.3)
   	aparam = runif(n, 100, 200)
   } else {
   	deltaR = rep(0, n)
   	aparam = rep(0, n)
   }

    par['deltaR', 1:n] = deltaR
    par['aparam', 1:n] = aparam

    parmin['deltaR', 1:n] = 1e-06
    parmin['aparam', 1:n] = 1

    parmax['deltaR', 1:n] = 4
    parmax['aparam', 1:n] = 1000

    pardx['deltaR', 1:n] = 0.05
    pardx['aparam', 1:n] = 0.05

    ## For school term
    if (any(mydata$imodel == c(2, 3))) {
        alpha = runif(n, 0.05, 0.2)
	} else {
        alpha = rep(0, n)
    }

    par["alpha", 1:n] = alpha
    parmin["alpha", 1:n] = 1e-06
    parmax["alpha", 1:n] = 1
    pardx["alpha" , 1:n] = 0.05

    ## Two-value model: delta, ts, dur
    parmin["delta", 1:n] = -1
    parmin["ts", 1:n] = tps[1]
    parmin["dur", 1:n] = 0.01

    parmax["delta", 1:n] = 1
    parmax["ts"   , 1:n] = tps[1] + 10 * day_per_week
    parmax["dur"  , 1:n] = 30 * day_per_week

    pardx["delta", 1:n] = 0.05
    pardx["ts"   , 1:n] = 0.05
    pardx["dur"  , 1:n] = 0.05

    if (mydata$imodel == 5) {
    	delta = runif(n, min = -0.2, max = 0.2)
    	ts = runif(n, min = (t0 + cadence * 10), max = (t0 + cadence * 20))
    	dur = runif(n, min = (cadence * 2), max = (cadence * 10))
    } else {
    	delta = rep(0, n)
    	ts    = rep(0, n)
    	dur   = rep(0, n)
    }

    par['delta', 1:n] = delta
    par['ts'   , 1:n] = ts
    par['dur'  , 1:n] = dur

	if (mydata$epi_model > 2) {
		pardx["bite", 1:n] = pardx["vec_k", 1:n] = pardx["muV", 1:n] = pardx["T_HV", 1:n] = pardx["T_VH", 1:n] = pardx["sigmaV", 1:n] = 0.05

		parmin["bite", 1:n] = 0.1
		parmax["bite", 1:n] = 10.
		par["bite", 1:n] = 0.5

		if (opt.list$bite)
			par["bite", 1:n] = runif(n, parmin["bite", 1], parmax["bite", 1])

		parmin["vec_k", 1:n] = 0.3
		parmax["vec_k", 1:n] = 4
		par["vec_k", 1:n] = 2

		if (opt.list$vec_k)
			par["vec_k", 1:n] = runif(n, parmin["vec_k", 1], parmax["vec_k", 1])

		parmin["muV", 1:n] = 1/42
		parmax["muV", 1:n] = 1/8
		par["muV", 1:n] = 1/25
		if (opt.list$muV)
			par["muV", 1:n] = runif(n, parmin["muV", 1], parmax["muV", 1])

		parmin["T_HV", 1:n] = parmin["T_VH", 1:n] = 0.33
		parmax["T_HV", 1:n] = parmax["T_VH", 1:n] = 1
		par["T_HV", 1:n]    = par["T_VH", 1:n] = 0.5
		if (opt.list$T_HV)
			par["T_HV", 1:n] = runif(n, parmin["T_HV", 1:n], parmax["T_HV", 1:n])
		if (opt.list$T_VH)
			par["T_VH", 1:n] = runif(n, parmin["T_VH", 1:n], parmax["T_VH", 1:n])

		parmin["sigmaV", 1:n] = 1/14
		parmax["sigmaV", 1:n] = 1/7
		par["sigmaV"   , 1:n] = 1/10

		if (opt.list$sigmaV)
			par["sigmaV", 1:n] = runif(n, parmin["sigmaV", 1], parmax["sigmaV", 1])

	}



	logvec = rep(0, length = nparam)
    logbase = 10  #use log base 10 when needed
	logvec[1:nparam] = 1
    logvec['t0'] = 0  #linear for t0
    logvec['delta'] = 0  # and linear for delta since it can be negative

    nlines = run.list$nlines
    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))


    list(parmin = parmin, parmax = parmax, pardx = pardx, par = par, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)


}


setupODEforSingleCDC <- function(mydata = NULL) {
    
    # Pack and setup some parameters for the ODEs
    #
    # @param mydata The complete list of all the data for the run
    # @return A named list with the following arguments:
    # \describe{
    # \item{tmax}{A numeric array with the day number when the basic reproduction   number is at its maximum.  Currently the code does not use this.}
    # \item{tps}{A numeric array with the day number for each week in the flu season}
    # \item{nperiods}{The number of weeks in the flu season (typically 52)}
    # \item{rtn}{A numeric 2D array with with nperiods rows and nregions columns}
    # \item{t0}{Numeric - the first day of the first week in the flu season}
    # }
    # @examples 
    # setupODEforSingleCDC{mydata = mydata}
    
    # set up the time-dependence of R(t) using the following formula set R(t) using the formula: R(t) = ((1-R0min/R0max)*
    # sin[2*pi/365*(t-t_max)+pi/2]+1+R0min/R0max)/2 where tmax is January 15 in the northern hemisphere and July 15 in the southern. To
    # change this go into calc.tmax using the latitude information for each pixel we determine tmax not really needed for region =
    # usa/cdc- it will come up as January 15 for all states or regions
    
    
    
    nregion = mydata$fit$nregions
    nregion1 = nregion + 1
    weeks = mydata$weeks
    nperiods = length(weeks)
    
   if (mydata$cadence == "Weekly") {
        cadence = 7
        tps = mydata$weeks
   } else if (mydata$cadence == "Monthly") {
        cadence = 31
        tps = mydata$months
   } else if (mydata$cadence == "Daily") {
        cadence = 1
        tps = mydata$days
   } else {
        tps = mydata$weeks
   }

    week0 = weeks[1]
        
    days_per_week =  7
    day0 =days_per_week * (week0 - 1)
    tps <- seq(from = day0, to = (day0 + nperiods * days_per_week), by = cadence)

    noPts = length(tps)
    rtn <- array(0, c(noPts - 1, nregion))

    ode_list = list(tps = tps, nperiods = nperiods, rtn = rtn, t0 = day0)
    
    return(ode_list)
    
}





