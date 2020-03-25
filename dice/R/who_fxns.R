## All functions realted to fitting WHO data

fitSingleWHO <- function(mydata = NULL, all_years_epi = all_years_epi, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine for an Uncoupled Spatial Model - WHO Data
    #'
    #' A spatially uncoupled MCMC fit of the model data. The code first fits the model data directly and then
    #' fits each of the sub-regions sequentially - minimizing the likelihood of each one. The final indirect
    #' model fit is obtained as a weighted sum of these individual fits with the weights given by the relative
    #' population of each region. The data can be either cdc or gft data, and the model/fit data should have
    #' different spatial scales. For example in the case of cdc/gft data: the model can be national and the fit are
    #' the ten HHS regions. Or the model can be an HHS region and the fit are the states in that region.
    #' @param mydata A complex list with all available mydata for a given flu season model/fit spatial levels and mydata type
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
	nfit = mydata$nperiodsFit
	
	prior = mydata$prior
	n = mydata$fit$nregion
	
	mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
	fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]
	
	mod_name = mydata$model$name
	fit_name = mydata$fit$name

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

	if (mod_level < fit_level) {

		tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
	}

	if (tolower(mydata$disease) == "flu") {
    	week0 = mydata$weeks[1]
    	day0 = cadence * (week0 - 1)
	    tps =  seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}

	## First fit the mod_level

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)
	epi.obsrv = mydata$fit$raw

	tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

	cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")
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


	epi.da.df = get.sql.da(nfit = nfit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
		corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

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

			epi.da.df = get.sql.da(nfit = nfit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$fit$epi[,
				iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

			mydata$fit$wght[, iregion] = epi.da.df$wght

			mydata$fit$epirun[, iregion] = epi.da.df$epi.new

			for (i in 1:length(mydata$fit$epirun[, iregion])) {
				gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
			}
			mydata$fit$gamaepirun[, iregion] = gamaepi

		}

	}

	nmydata = length(tps)

	nperiods = mydata$nperiods

	wght = mydata$model$wght

	par_names <- set.param.list(epi_model = mydata$epi_model)

    setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

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

    setup = setup.fit.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = tps, par_names = par_names)


    fit_pmin = setup$parmin
    fit_pmax = setup$parmax
    fit_dx = setup$pardx
    fit_par = setup$par

    if (mydata$data_source == 15) {
    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 1.4

    	fit_pmax["pC", 1:n] = 0.2
    	fit_pmax["R0", 1:n] = 1.4
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

    imask = set.imask(par_names = par_names, opt.list = opt.list)

    # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL
    tab.fit = list()

    pois = 0
    nblock = 3
    accept.vec = array(0, c(n, nblock))

    if (mydata$sql_db == TRUE) {
    	cases = mydata$model$epirun
    	gamaepi = mydata$model$gamaepirun
    	wght = mydata$model$wght
    } else {
   		cases = mydata$model$epi
    	gamaepi = mydata$model$gamaepi
    	wght = mydata$model$wght
    }

    sh = mydata$model$sh
    school = mydata$model$school

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

    ## Need nperiods + 2 because we pathc the weeks at the beginning and end
    ndays = (nperiods + 2) * cadence

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

    profile = array(data = 0, dim = c(nRnd, nperiods, n))

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
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

        if (mydata$sql_db == TRUE) {
        	cases = mydata$fit$epirun[, iregion]
        	gamaepi = mydata$fit$gamaepirun[, iregion]
        	wght = mydata$fit$wght[, iregion]
        } else {
        	cases = mydata$fit$epi[, iregion]
        	gamaepi = mydata$fit$gamaepi[, iregion]
        	wght = mydata$fit$wght[, iregion]
        }

        sh = mydata$fit$sh[, iregion]
        school = mydata$fit$school[, iregion]

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
            scales = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
            ithin = as.integer(ithin),nperiods = as.integer(nperiods),
            tps = as.double(tps[1:nperiods]), rtn = as.double(rtn[, iregion]), accept = as.double(rep(0, nblock)), pois = as.double(pois),
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

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

    success = writeCSV(mydata = mydata, run.list = run.list, model_rtn = model_rtn, model_profile = model_profile,
        rtn = rtn, profile = profile, ireal = ireal)
    list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)



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

	err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)

 err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)
 
    err <- saveProfiles(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, rtn = rtn,  profile = profile,  ireal = ireal, run.list = run.list)

	err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)

    list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)



}

fitOnePatchWHO <- function(mydata = NULL, all_years_epi = all_years_epi, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine for a single region - WHO Data
    #'
    #' An MCMC fit of the model WHO data.
    #' @param mydata A complex list with all available mydata for a given flu season model
    #' @param all_years_epi - the entire history of incidence for this region
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
    #' }
    #' @examples
    #' fitOnePatchWHO{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}

 	mydata$data_source = 'who_flu'

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
    nfit = mydata$nperiodsFit

    prior = mydata$prior

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    mod_name = mydata$model$name
    fit_name = mydata$fit$name

	mod_name = mydata$model$name

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)

	if (tolower(mydata$disease) == "flu") {
    	week0 = mydata$weeks[1]
    	day0 = cadence * (week0 - 1)
	    tps =  seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}

	## mod level

	epi.obsrv = mydata$model$raw
	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)

	cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

	## Fit the mod_level
	## Start by calculating correlation/distance with past years

	corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)

	distMAT = calc.sqldata.dist(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = nfit)


	epi.da.df = get.sql.da(nfit = nfit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
		corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

	mydata$model$wght = epi.da.df$wght

	mydata$model$epirun = epi.da.df$epi.new

	gamaepi = rep(0, mydata$nperiods)
	for (i in 1:length(mydata$model$epirun)) {
		gamaepi[i] = lgamma((mydata$model$epirun[i] + 1))
	}
	mydata$model$gamaepirun = gamaepi

	nmydata = length(tps)

	wght = mydata$model$wght

	par_names <- set.param.list(epi_model = mydata$epi_model)

    setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list, tps = tps, par_names = par_names)

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

    if (mydata$data_source == "who") {
    	model_pmax["pC"] = 0.2
    	model_pmax["R0"] = 1.4
    }

    ## And let's sample pC and R0 from the prior

    ## We get the mean and sigma of the prior but we actually use it only with prior = 1 or 2 This is controlled by logvec see if
    ## statement below

    setup.prior = setup.single.prior(mydata = mydata, logvec = logvec)

    n1 = length(setup.prior$ymu[, 1])

    model_ymu = setup.prior$ymu[n1, ]
    model_sigma = setup.prior$sigma[n1, ]

    imask = set.imask(par_names = par_names, opt.list = opt.list)

    # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    tab.model = NULL

    pois = 0
    nblock = 3

    if (mydata$sql_db == TRUE) {
    	cases = mydata$model$epirun
    	gamaepi = mydata$model$gamaepirun
    	wght = mydata$model$wght
    } else {
   		cases = mydata$model$epi
    	gamaepi = mydata$model$gamaepi
    	wght = mydata$model$wght
    }

    sh = mydata$model$sh
    school = mydata$model$school

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

    ## Need nperiods + 2 because we pathc the weeks at the beginning and end
    ndays = (nperiods + 2) * cadence

    nRnd = 1000

    model_profile = array(0, c(nRnd, nperiods))

    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(mydata$Temp), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd), epi_model = as.integer(epi_model))

	model_factor = mydata$model$factor

    model_rtn = solution$rtn
    model_rtn_ili = model_rtn/model_factor
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile, c(nRnd, nperiods))
	model_profile_ili = model_profile/model_factor

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	tables = calc.mech.err(tables = tables.mod, profiles = model_profile, nfit = nfit, state_id = mod_id)
	tables.mod = tables

    # # write the MCMC statistics

    success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, mydata = mydata, imask = imask, ireal = ireal)

    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

    ## Plot all of the results - both profiles and histograms - the device list can have more than one element

     tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	 tables.mod = tables

	 tables = calc.mech.err(tables = tables.mod, profiles = model_profile_ili, nfit = nfit, state_id = mod_id)
	 tables.mod = tables

    err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)
    
 err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)    
    ##
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

fitMultiWHO <- function(mydata = mydata, all_years_epi = NULL, run.list = NULL, opt.list = NULL, ireal = 1, iseed = NULL) {

  #' Driver Routine Coupled Spatial Model - WHO Data
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
  #' fitMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed =iseed}


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
  # Here we set the random number generation seed for R
  set.seed(iseed)

  nregion = n = mydata$fit$nregion

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

  opt.cpl = opt.list[[2]]
  opt.list = opt.list[[1]]

  nperiods = mydata$nperiods
  nfit = mydata$nperiodsFit

  nreal  = run.list$nreal
  nMCMC  = run.list$nMCMC
  nlines = run.list$nlines
  ithin = run.list$ithin
  device = run.list$device
  subDir = run.list$subDir

  mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
  fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

  mod_name = mydata$model$name
  fit_name = mydata$fit$name

  tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)


  if (mod_level < fit_level) {

    tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = mydata$nperiods)
  }

	if (tolower(mydata$disease) == "flu") {
    	week0 = mydata$weeks[1]
    	day0 = cadence * (week0 - 1)
	    tps =  seq(from = day0, to = (day0 + nperiods * cadence), by = cadence)
	}

  ## First fit the mod_level

  epi.obsrv = mydata$model$raw
  tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)
  epi.obsrv = mydata$fit$raw

  tables.fit$epi.obsrv[,fit_id]= as.matrix(epi.obsrv)

  cat("\n Calculating NULL Model (Historic Monthly or Weekly Average)\n\n")

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

  epi.da.df = get.sql.da(nfit = nfit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$model$epi, epi = mydata$model$epi,
                         corrMAT = corrMAT, deng.null = epi.null.mod, my_id = mod_id)

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

      epi.da.df = get.sql.da(nfit = nfit, all_years_epi = all_years_epi, mydata = mydata, cases = all_years_epi$fit$epi[,
                                                                                                                      iregion], epi = mydata$fit$epi[, iregion], corrMAT = corrMAT, deng.null = epi.null.fit, my_id = fit_id[iregion])

      mydata$fit$wght[, iregion] = epi.da.df$wght

      mydata$fit$epirun[, iregion] = epi.da.df$epi.new

      for (i in 1:length(mydata$fit$epirun[, iregion])) {
        gamaepi[i] = lgamma((mydata$fit$epirun[i, iregion] + 1))
      }
      mydata$fit$gamaepirun[, iregion] = gamaepi

    }

  }


  nperiodsFit = mydata$nperiodsFit

  wght = mydata$model$wght

  par_names <- set.param.list(epi_model = mydata$epi_model)

  par_names_cpl <- set.param.cpl.list()

  setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
                           tps = tps, par_names = par_names)

  model_pmin = setup$parmin
  model_pmax = setup$parmax

  model_dx = setup$pardx
  model_par = setup$par

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
                         tps = tps, par_names = par_names)

  fit_pmin = setup$parmin
  fit_pmax = setup$parmax
  fit_dx = setup$pardx
  fit_par = setup$par

  ## Just for the flu prediction limit pC and R0
  model_pmax['pC'] = 0.2
  model_pmax['R0'] = 2
  fit_pmax['pC', 1:n] = 0.2
  fit_pmax['R0', 1:n] = 2

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

  pois = 0
  nblock = 3

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

  ## Number of days - and patch it with a month or a week before and after

  ndays = (nperiods + 2) * cadence

  ymu = rep(0, nparam)
  sigma = rep(1, nparam)

	# Also check for NAs
	school[is.na(school)] <- 0
	sh[is.na(sh)]         <- 0

  out <- .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
                  wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
                  dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), Temp = as.double(mydata$Temp),
                  imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
                  nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
                  rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(1e5), tab.model = as.single(tab.model),
                  imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd),epi_model=as.integer(mydata$epi_model))

  cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")

  model_rtn = out$rtn
  tab.model = matrix(out$tab.model, ncol = (nparam + 1))
  model_profile = array(out$profile, c(nRnd, nperiods))


	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

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

  cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using Level ", mydata$fit_level, " Data\n")

	if (is.null(school)) school = matrix(data=0, nrow=mydata$nperiods, ncol=mydata$fit$nregions)

	for (iregion in 1:nregion) {
		# Also check for NAs
		school[is.na(school[, iregion]), iregion] <- 0
		sh[is.na(sh[, iregion]), iregion]         <- 0
	}

  solution = .Fortran("epimulti", n = as.integer(n) , epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
                      wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin), parmax = as.double(fit_pmax),
                      dx = as.double(fit_dx), ilog = as.integer(logvec), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
                      ithin = as.integer(ithin), nmydata = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
                      rtn = as.double(rtn), pois = as.double(rep(0, n)),
                      coef = as.double(mydata$fit$coef), parCPL = as.double(cpl_par), pminCPL = as.double(cpl_pmin),
                      pmaxCPL = as.double(cpl_pmax), stepCPL = as.double(cpl_dx), ilogCPL = as.integer(cpl_logvec), imaskCPL = as.integer(cpl_imask),
                      Rij = as.double(Rij), tab = as.single(tab), tabCPL = as.single(tab.cpl), profiles = as.single(profile), nRnd = as.integer(nRnd), ndays = as.integer(ndays), epi_model=as.integer(epi_model), scales = as.double(mydata$Temp))


  cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using Level ", mydata$fit_level, " Data Completed \n")
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
  tables = calc.mech.err(tables = tables.agg.mod, profiles = fit_model_profile, nfit = nfit, state_id = mod_id)
  tables.agg.mod = tables

  ## write the MCMC statistics

  success = mcmc.multi.write(tab.model = tab.model, tab.fit = tab.fit, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list,
                             mydata = mydata, imask = imask, cpl_imask = cpl_imask, ireal = ireal)


  ## Calculate the NULL model
  tables = calc.null.err(tables = tables.fit, nfit = nfit, state_id = fit_id)
  tables.fit = tables

  ## Update the results table

  for (iregion in 1:n) {
    tables = calc.mech.err(tables = tables.fit, profiles = profile[,,iregion], nfit = nfit, state_id = fit_id[iregion])
    tables.fit = tables
  }

  ## save a binary RData file of the input of this run

  run.input = list(run.list = run.list, opt.list = opt.list)

  success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)

  ## Plot the results and the posterior

  err <- plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal, run.list = run.list, idevice = 1)

 err <- saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal)
 
  err <- plotMCMC(mydata = mydata, tab.model = tab.model, tab.fit= tab, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)

    err <- saveProfiles(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, rtn = rtn,  profile = profile,  ireal = ireal, run.list = run.list)

  list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab)

}

