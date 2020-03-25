set.sirb.model <- function(mydata = NULL, par_names = NULL, est_bckgrnd = 10, opt.list = NULL, run.list = NULL) {
  #' Setup of Parameters for an MCMC sirb procedure
  #'
  #' \code{set.sirb.model} Creates the arrays for the SIRB model with min/max values, step size,
  #' and initial guess/default values for all the parameters.
  #' @param mydata A list containing all available information for the region: incidence, specific-humidity, school schedule, population etc.
  #' @param par_names A character array with SIRB parameter names
  #' @param est_bckgrnd Numeric, estimated value of background cholera cases
  #' @param run.list A list with parameters needed for the MCMC procedure  
  #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters
  #'   supported by \pkg{DICE}. These values are based on the model used for the basic
  #'   reproduction number  
  #' @return A list with min, max, step size and initial guess/default values for the sir parameters
  #' @examples
  #' set.sirb.model(mydata = mydata, par_names = par_names,
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
		par_names <- c("NH", "Tg", "birth", 't0', 'pC', "e_bckgrnd",'seed', 'a', 'K','growth_loss', 'e_shedd', 'delta', 'ts')

	}
	
	nparam = length(par_names)
		
	parmin = parmax = pardx = par = logvec = rep(0.0, length = nparam)
	names(parmin) = names(parmax) = names(pardx) = names(par) = names(logvec) = par_names
   	par_opt = opt.list

   nopt = length(par_opt[par_opt == TRUE])
	
	pop = mydata$model$pop
		
	Tg = mydata$Tg 
	
	growth_loss = 0.33
	e_shedd = 10. 
		
	if (mydata$data_source == tolower("SOM_MH")) {
		NH = 2e5
		a = 0.5
		K = 5.5
		delta = 0.0
		ts = 1
		t0 = 7.
		pC = 1 # 0.1
		e_bckgrnd = est_bckgrnd
		seed = 2
		dx = 0.001
	} else if (mydata$data_source == tolower("ZMB_MH")) {
		NH = 2e4
		a = 1
		K = 5.5
		delta = 1.5
		index = which(mydata$weeks == 9)
		index1 = index + 1
		ts   = index *  cadence - 3.5
		t0 = 7.
		pC = 1 # 0.1
		e_bckgrnd = est_bckgrnd
		seed = 5.
		growth_loss = 0.7
		e_shedd = 2.
		dx = 0.01
	} else {
		NH = 2e5
		a = 0.5
		K = 5.5
		delta = 0.0
		ts = 1
		t0 = 7.
		pC = 1 # 0.1
		e_bckgrnd = est_bckgrnd
		seed = 2
		dx = 0.01
	}
	

	birth = 1/(60 * 365)
	
	
	par.input <- c(NH = NH, Tg = Tg, birth = birth, t0 = t0, pC = pC, e_bckgrnd = e_bckgrnd , seed = seed, a = a, K=K, growth_loss = growth_loss, e_shedd = e_shedd, delta = delta, ts = ts)

	##
	## min/max values for each parameter
	##
	parmin = c(NH = 1e5, Tg = Tg, birth = birth, t0 = 0, pC = 1e-5, e_bckgrnd = 1 , seed = 1, a = 0.1, K = 1, growth_loss = 0.01, e_shedd = 1, delta = 0.01, ts = 100)

	parmax = c(NH = 1e6, Tg = Tg, birth = birth, t0 = 21., pC = 1, e_bckgrnd = 100, seed = 50, a = 2, K = 10, growth_loss = 5., e_shedd = 50, delta = 1, ts= 200)

	if (mydata$data_source == tolower("ZMB_MH")) {
		parmin = c(NH = 1e3, Tg = Tg, birth = birth, t0 = 0, pC = 1e-5, e_bckgrnd = 1 , seed = 1, a = 0.9, K = 1, growth_loss = 0.01, e_shedd = 1, delta = 0.01, ts = 100)

		parmax = c(NH = 1e5, Tg = Tg, birth = birth, t0 = 21., pC = 1, e_bckgrnd = 20, seed = 20, a = 2, K = 10, growth_loss = 1., e_shedd = 20, delta = 2, ts = 200)
	}
	pardx = c(NH = dx, Tg = dx, birth = dx, t0 = dx, pC = dx, e_bckgrnd = dx, seed = dx, a = dx, lamda = dx, growth_loss = dx, e_shedd = dx, delta = dx, ts = dx)
	
	
   	logbase = 10  #use log base 10 when needed
	logvec[1:nparam] = 1
    logvec['t0'] = 0  #linear for t0
    
    nlines = run.list$nlines
    tab <- matrix(data = 0, nr = nlines, nc = (nparam + 1))


    setup = list(parmin = parmin, parmax = parmax, pardx = pardx, par = par.input, logbase = logbase, logvec = logvec, nparam = nparam,
        nopt = nopt, tab = tab)
        
    return(setup)
        	
}


fitCholera <- function(mydata = NULL, all_years_epi = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine - Fitting a Single Cholera Incidence Profile
    #'
    #' A driver for fitting a single cholera incidence profile of a region/patch. The R code calls a Fortran routine which
    #' uses an MCMC procedure to maximize the likelihood of the solution using and SIRB mode.
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
    #' fitCholera{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}


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
	
	## Lets padd it with 4 weeks of zero cases
	
	if (mydata$data_source == tolower("SOM_MH")) {

		nadd = sum(is.na(mydata$model$raw))

		tps_add = 1:nadd 

		tps = tps + nadd 
		
		tps = c(tps_add, tps)

		tps = tps * cadence 
		
		cases.model = c(rep(0, nadd), cases.model)
		
		est_bckgrnd = 10
		
	}

	if (mydata$data_source == tolower("ZMB_MH")) {

		nadd = 0

		tps = 1:length(tps) * cadence 
				
		est_bckgrnd = 10
		
	}

	par_names <- set.param.list(epi_model = mydata$epi_model)	
	
	setup  <- set.sirb.model(par_names = par_names, mydata = mydata, est_bckgrnd = est_bckgrnd, opt.list = opt.list, run.list = run.list)
	
	
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

	gamaepi = cases.model

	ndata = length(gamaepi)
	
	for (i in 1:ndata) {
		gamaepi[i] = lgamma((cases.model[i] + 1))
	}	

	nfit = mydata$nperiodsFit
	wght = rep(0, ndata)
	wght[1:(nfit+nadd)] = 1
	if ((nfit+nadd) < length(tps)) wght[(nfit+nadd+1):ndata] = 0
	
	nRnd = 1000
	profiles = profilesB = array(0, c(ndata, nRnd))

	## Need nperiods + 2 because we pathc the weeks at the beginning and end

	nmydata = length(tps)

	ndays = (ndata + 2) * cadence
	
	cat("\n ****** Fitting ", mydata$model$name, " Using an MCMC Procedure ****** \n")

solution <- .Fortran("fitsirb", nperiods = as.integer(ndata), tps = as.double(tps), epi = as.double(cases.model), gamaepi = as.double(gamaepi), wght = as.double(wght), 
	nparam = as.integer(nparam), par = as.double(model_par), parmin = as.double(model_pmin), parmax = as.double(model_pmax), step = as.double(model_dx), ilog = as.integer(logvec), 
	imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin), tab = as.single(tab), profiles = as.single(profiles), profilesB = as.single(profiles), nRnd = as.integer(nRnd), Temp = as.double(mydata$Temp), ndays = as.integer(ndays))


	profiles = array(solution$profiles, c(ndata, nRnd))

	profilesB = array(solution$profilesB, c(ndata, nRnd))

	clean.profiles = profiles
	clean.profilesB = profilesB
	
	if (mydata$data_source == tolower("SOM_MH")) {
		irmv = which(profiles[1, 1:nRnd] > 1000) # max(cases.model/2)
	
		if (length(irmv) > 0) {
	
			clean.profiles = profiles[, -irmv]
			clean.profilesB = profilesB[, -irmv]
	
		}
	}
		
	if (mydata$data_source == tolower("ZMB_MH")) {
		irmv = which(profiles[1, 1:nRnd] > 200) # max(cases.model/2)
		irmv2 = which(profiles[2, 1:nRnd] > 200)
		irmv3 = which(profiles[3, 1:nRnd] > 200)
		print(irmv)
		print(irmv2)
		print(irmv3)

		irmv = c(irmv, irmv2, irmv3)
		irmv = unique(irmv)

		if (length(irmv) > 0) {

			clean.profiles = profiles[, -irmv]
			clean.profilesB = profilesB[, -irmv]
		}

	}

	clean.profiles = clean.profiles[(nadd+1):ndata,]
	clean.profilesB = clean.profilesB[(nadd+1):ndata,]
	
	model_rtn = rep(NA, mydata$nperiods)
	for (i in 1:mydata$nperiods) model_rtn[i] = mean(clean.profiles[i,])

    tab.model =  matrix(solution$tab, ncol = (nparam + 1))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

	tables = calc.sirb.err(profiles = clean.profiles, profilesB = clean.profilesB, mydata = mydata)
	tables.mod = tables

	model_factor = mydata$model$factor
	
	model_profile = t(clean.profiles) 
	model_profile_ili = model_profile/model_factor
	
	model_profileB = t(clean.profilesB) 
	
	err <- plotSIRB(tables.mod = tables, profiles = clean.profiles, profilesB = clean.profilesB, mydata = mydata, run.list = run.list, ireal = ireal, idevice = 1)

    		
 	err <- saveTables(mydata = mydata, tables.mod = tables, tables.fit = NULL, tables.agg.mod = NULL, mod_id = mod_id, fit_id = fit_id, run.list = run.list, ireal = ireal) 
 		
    # # write the MCMC statistics

    success = mcmc.onepatch.write(tab.model = tab.model, opt.list = opt.list, run.list = run.list, mydata = mydata, imask = imask, ireal = ireal)


    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = all_years_epi, run.input = run.input, ireal = ireal)
    ## now dump all the profiles we have to a file    - in the case of the CDC this is done in the ploting routine
    ## 
    	err <- saveProfilesOnePatch( mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, model_profileB = model_profileB, ireal = ireal, run.list = run.list)
    		
	
	##
	## Plot the posterior distribution of the parameters
	##
	
    err <- plotMCMC(mydata = mydata, tab.model = tab.model, opt.list = opt.list, run.list = run.list, imask = imask, ireal = ireal, idevice = 1)


    results = list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model)

    return(results)

}


calc.sirb.err <- function(profiles=NULL, profilesB =NULL, mydata = NULL) {

	ntraj = dim(profiles)[2]
	ndays = dim(profiles)[1]

	sirb.lower = sirb.upper = sirb.frcst = sirb.bact = sirb.lower.bact = sirb.upper.bact = rep(0, ndays)
		for (iday in 1:ndays) {
			sirb.lower[iday] = quantile(profiles[iday, ], prob = 0.05)
			sirb.upper[iday] = quantile(profiles[iday, ], prob = 0.95)
			sirb.frcst[iday] = mean(profiles[iday,])
			sirb.bact[iday]  = mean(profilesB[iday,])
			sirb.lower.bact[iday] = quantile(profilesB[iday, ], prob = 0.05)
			sirb.upper.bact[iday] = quantile(profilesB[iday, ], prob = 0.95)
			
		}


	tables = list()
	tables$observed   = mydata$model$raw
	tables$sirb.lower = sirb.lower
	tables$sirb.upper = sirb.upper
	tables$sirb.frcst = sirb.frcst
	tables$sirb.bact  = sirb.bact
	tables$sirb.lower.bact = sirb.lower.bact
	tables$sirb.upper.bact = sirb.upper.bact	
	return(tables)

}


plotSIRB <- function(tables.mod = NULL, profiles = NULL, profilesB = NULL, mydata = NUll, ymax.input = NULL, ireal = NULL, run.list = NULL, idevice = 1) {

    device = run.list$device[idevice]

    if (is.null(device))
        device = "png"
    if (is.null(tables.mod))
        return

        
    subDir = mydata$subDir

    Tg     = mydata$Tg

    model = mydata$imodel

	isingle = mydata$isingle

    nperiods = mydata$nperiods

    nperiodsFit = mydata$nperiodsFit

    if (mydata$epi_model == 1) {
    	epi_model = 'SIR'
    } else if (mydata$epi_model == 2) {
    	epi_model = 'SEIR'
    } else if (mydata$epi_model == 3) {
    	epi_model = 'V-SIR'
    } else if (mydata$epi_model == 4) {
    	epi_model = 'V-SEIR'
    } else if (mydata$epi_model == 5) {
    	epi_model = 'SIRB'
    } else {
    	epi_model = 'UNKNOWN'
    }


   ## subset only the state we are plotting
   
    cadence = mydata$cadence

    if (cadence == "Monthly") {
        month.names = month.abb[mydata$months]
        dates = month.names  #  paste0(month.name,'-',year.vec)
        ind = seq(from = 1, to = length(dates))
    } else if (cadence == "Weekly"){
        dates = mydata$weeks
        dates = sub("week", "EW", dates)
        ind = seq(1, length(dates), by = 4)
        dates = dates[ind]
    } else {
    	quit("Unknown cadence in plotSIRB")
    }

    nregions = mydata$fit$nregions
    nregions1 = nregions + 1
    nregions2 = nregions + 2

    ndata = mydata$nperiods

    nperiodsFit = mydata$nperiodsFit
    myName = mydata$dataName
    myName = gsub(" ","",myName)
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

    if (nregions > 1) {
    	ncol = 3
    	nrow = round(nregions2 / ncol)
    	if ((nrow * ncol) < nregions2) nrow = nrow + 1
    	if (nrow > 3) nrow = min(nrow,3)
    	if (nregions == 10) {
    		nrow = 3
    		ncol = 4
    	}
    	par(mar = c(4, 4, 1, 1), mfrow = c(nrow, ncol))
    } else {
    	nrow = 1
    	ncol = 1
    	par(mar = c(5, 5, 3, 2), mfrow = c(nrow, ncol))
    }


	
	sirb.lower = tables.mod$sirb.lower
	sirb.upper = tables.mod$sirb.upper
	sirb.frcst = tables.mod$sirb.frcst
	sirb.bact  = tables.mod$sirb.bact
	observed   = tables.mod$observed

	sirb.bact  = tables.mod$sirb.bact
	sirb.lower.bact  = tables.mod$sirb.lower.bact
	sirb.upper.bact  = tables.mod$sirb.upper.bact
	
	
	title = paste0(mydata$model$name, " - ", mydata$FY, " season")

	xlab = ""

	xat = 1:ndata

	xlab = dates
	
	ymax = max(sirb.upper, observed, na.rm=TRUE)
	ymin = 0
	plot(1:ndata, observed, type = "n", col = "red", lwd = 1, xlab = "Epi Week", ylab = "Suspected Cases", ylim = c(ymin, ymax), main = title, xlim = c(1, ndata), xaxt = "n", cex=1., cex.lab=1., cex.axis=1)
	polygon(c(1:ndata, rev(1:ndata)), c(sirb.upper, rev(sirb.lower)), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
	lines(1:ndata, sirb.frcst, type = "l", col = "blue", xlab = "", ylab = "", lwd = 2)
	lines(1:nperiodsFit, observed[1:nperiodsFit], type = "h", col = rgb(1,0,0, alpha=0.7), lwd = 10, lty = 1, xlab = "", ylab = "")
	lines((nperiodsFit:ndata), observed[nperiodsFit:ndata], type = "h", col = rgb(1,0,0, alpha=0.4), lwd = 10, lty = 1, xlab = "", ylab = "")	
    rect(xleft = nperiodsFit, ybottom = 0, xright = (ndata), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)	
	legend("topleft", c("Observed", "SIR-B", "V. Cholera Concentration"), text.col = c("red", "blue", "green"), bty = "n")
	axis(1, at = 1:ndata, labels = FALSE, las = 2)
	axis(1, at = ind, labels = xlab, las = 1, cex.lab = 1., cex.axis=1., cex.lab=1.)
	par(new = TRUE)
	ymax = max(sirb.upper.bact, na.rm=TRUE)
	plot(1:ndata, sirb.bact, type = "l", col = "green", lwd = 2, xlab = " ", ylab = " ",xlim = c(1, ndata), ylim = c(ymin, ymax), xaxt = "n", yaxt='n')	
	polygon(c(1:ndata, rev(1:ndata)), c(sirb.upper.bact, rev(sirb.lower.bact)), col = rgb(0, 1., 0., 0.2), border = FALSE)
	axis(4, col.axis = 'green')		
	axis(4, col.axis = 'green')	
	mtext("V. Cholera Concentraion (cells/ml) ", side = 4, line = 3, cex.lab = 1, col = "green")
	dev.off()
	
	return(err=0)
}

