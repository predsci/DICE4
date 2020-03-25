###
### Functions related to Running DICE with Synthetic Data
###

psi.sim.modelA <- function(R0 = 1.4, Tg = 2.6, t0 = 20, seed = 1, N = 10000, tps = (0:52) * 7, reals = 20, dt = 0.1, stochastic = TRUE,
    plot = TRUE) {

    #' A Compartmental SIR Model
    #'
    #' An R code for generating SIR profile(s) using either a stochastic or deterministic event-driven procedure.
    #'   The routine calculates the weekly number of cases and if requested plots them.
    #'   When using a stochastic procedure the user should set the number of realization to be > 50.
    #'   When using a deterministic procedure the user should only ask for only one realization, see examples below.
    #' @param R0 Numeric - the basic reproduction number
    #' @param Tg Numeric - the recovery rate in days
    #' @param t0 Numeric - the onset day of the epidemic
    #' @param seed Integer - the initial number of cases
    #' @param N Integer - the total population
    #' @param tps a numeric array of time points in days
    #' @param reals Integer - the number of requested realizations.
    #' @param dt Numeric - the time-step for the event-driven procedure
    #' @param stochastic Logical - if TRUE (default) use a stochastic event-driven procedure, otherwise use a deterministic procedure
    #' @param plot Logical - if TRUE plot the results
    #' @return rtn   a numeric array with the weekly incidence of all requested realizations
    #' @examples
    #' psi.sim.modelA(R0 = 1.4, Tg = 2.6, t0 = 20, seed = 1, N = 1e4, tps = (0:52) * 7, reals = 100, dt = 0.2, stochastic = TRUE , plot = TRUE)
    #' psi.sim.modelA(R0 = 1.4, Tg = 2.6, t0 = 20, seed = 1, N = 1e4, tps = (0:52) * 7, reals = 1  , dt = 0.2, stochastic = FALSE, plot = TRUE)

    noTPts <- length(tps)
    epsilon <- 1e-10
    if (seed >= N)
        stop("psi.sim.modelA seed greater than N")

    if (noTPts < 3)
        stop("simSIR: tps must be of length 2 or greater")
    if (reals > 1 && stochastic == FALSE)
        stop("simSIR: why make multiple realisations of a deterministic model?")

    rtn <- array(0, c(noTPts - 1, reals))

    for (i in 1:reals) {

        t_cur <- tps[1]
        ind_t <- 2
        t_next <- tps[ind_t]

        sS <- N
        sI <- 0
        sR <- 0

        blNotYetSeed <- TRUE

        while (ind_t <= noTPts) {

            if (t_cur >= t0 && blNotYetSeed) {
                seedInf <- seed
                blNotYetSeed <- FALSE
            } else {
                seedInf <- 0
            }

            lambda <- R0 * sI/Tg/N
            pInf <- 1 - exp(-lambda * dt)
            pRec <- 1 - exp(-dt/Tg)


            if (stochastic) {
                if (lambda > 0 || seedInf > 0) {
                  nInf <- rbinom(1, sS, pInf) + seedInf
                } else {
                  nInf <- 0
                }

                nRec <- rbinom(1, sI, pRec)
            } else {
                if (lambda > 0) {
                  nInf <- sS * pInf + seedInf
                } else {
                  nInf <- 0
                }

                nRec <- sI * pRec

            }

            sS <- sS - nInf
            sI <- sI + nInf - nRec
            sR <- sR - nRec


            rtn[ind_t - 1, i] <- rtn[ind_t - 1, i] + nInf

            t_cur <- t_cur + dt

            if (t_cur > (t_next - epsilon)) {
                t_next <- tps[ind_t]
                ind_t <- ind_t + 1

            }

        }
    }

    if (plot) {
        plot(rtn[, 1], type = "l", col = "grey", ylim = c(0, max(rtn)), xlab = "Week Number", ylab = "Number of Cases")
        if (reals > 1) {
            for (j in 2:reals) points(rtn[, j], type = "l", col = "grey")
        }
    }

    return(allruns = rtn)

}

runSyntheticSingle <- function(mydata = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine Uncoupled Spatial Model - Synthetic Data
    #'
    #' A spatailly uncoupled MCMC fit of synthetic data. The code first fits the model data directly and then
    #'   fits each of the sub-regions sequentially - minimizing the likelihood of each one. The final indirect
    #'   model fit is obtained as a weighted sum of these individual fits with the weights given by the relative
    #'   population of each region. The data is synthetic (although it has the dataType 'cdc'). This synthetic example
    #'   has been set up for the case of fit_level = 3 and mod_level = 2, that is fit the national data using the ten
    #'   uncoupled HHS regions
    #' @param mydata A complex list with all the data for a single  season. Everything except the number of cases and \%ILI
    #'  is from the 'cdc' dataType and will not be replaced.  Only the number of cases and the \%ILI will be replaced
    #'  by syntehtic profiles
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \pkg{DICE}.
    #'  The synthetic example has been set-up for model = 4
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' #' @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    #' @return A list with the following arguments:
    #' \describe{
    #' \item{model_rtn}{The best result of the MCMC procedure for the model synthetic data}
    #' \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model synthetic data}
    #' \item{tab.model}{The MCMC history of the direct fit to the model synthetic data}
    #' \item{fit_rtn}{The best result for indirectly fitting the model synthetic data using the synthetic fit regions}
    #' \item{fit_profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
    #' \item{tab.fit}{The MCMC history of indirectly fitting the model data using the fit regions}
    #' }
    #' @examples
    #' runSyntheticSingle{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = 1}


    # Take only the first list from opt.list, the second is optional and is needed only for a coupled run

    opt.list = opt.list[[1]]
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    prior = mydata$prior
    season = mydata$season
    FY = mydata$FY
    Temp = mydata$Temp

    cadence = mydata$cadence

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed()
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generation seed for R
    set.seed(iseed)


    out <- setupODEforSingleCDC(mydata = mydata)

    tps = out$tps

    n = mydata$fit$nregions
    rtn = out$rtn

    time = out$t0

    t0 = rep(out$t0, n)

    tps = out$tps

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    # Random initial guess for the model parameters
    setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = out$tps, par_names = par_names)


    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx

    ## These random guesses for the initial parameters will be used in the optimization process
    model_par = setup$par


    # these are the same for both model and fit

    nparam = setup$nparam
    nopt = setup$nopt
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


    # replace some values

    fit_par[1, 1:n] = round(runif(n, min = t0[1] + 1, max = t0[1] + cadence * 1))  #t0
    fit_par[2, 1:n] = runif(n, min = 0.02, max = 0.1)  #pC
    fit_par[3, 1:n] = runif(n, min = 1.2, max = 1.3)  #R0
    for (j in 1:n) fit_par[8, j] = min(mydata$fit$epi[1, j], mydata$fit$epi[nperiods, j])

    ## Save a copy of the input parameters - they will be used when plotting the results
    save_fit_par = fit_par
    colnames(save_fit_par) = mydata$fit$name
    rownames(save_fit_par) = names(opt.list)

    # generate a synthetic profile

    fit_epi = array(0, c(nperiods, n))
    for (i in 1:n) {
        mypar = fit_par[, i]
        sh = mydata$fit$sh[, i]
        school = mydata$fit$school[, i]
        pop = mydata$fit$pop[i]
        solution = .Fortran("gensingle", sh = as.double(sh), school = as.double(school), nparam = as.integer(length(mypar)), par = as.double(mypar),
            nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]), rtn = as.double(rep(0, nperiods)), ndays = as.integer(tps[nperiods]))
        fit_epi[, i] = solution$rtn

    }

    # Build the national number of cases as the sum of the HHS regions
    model_epi = rep(0, nperiods)
    for (i in 1:nperiods) model_epi[i] = sum(fit_epi[i, 1:n])

    ## Create %ILI for model and fit using cases and factor
    fit_ili = fit_epi
    model_ili = model_epi
    for (i in 1:nperiods) {
        for (j in 1:n) fit_ili[i, j] = fit_ili[i, j]/as.numeric(mydata$fit$factor[j])
        model_ili[i] = model_ili[i]/mydata$model$factor
    }

    model_epi = round(model_epi)
    fit_epi = round(fit_epi)
    fit_gamma = gammaEpi(epi = fit_epi)
    model_gamma = lgamma((model_epi + 1))

    synData = mydata
    synData$model$raw = model_ili
    synData$model$epi = model_epi
    synData$model$gamaepi = model_gamma

    synData$fit$raw = fit_ili
    synData$fit$epi = fit_epi
    synData$fit$gamaepi = fit_gamma

    saveData = mydata
    mydata = synData

    ## These are not used in the case of synthetic data but we need to allocate the arrays
    ymu = rep(0, nparam)
    sigma = rep(1, nparam)

    imask = set.imask(par_names = par_namses, opt.list = opt.list)

    # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    # Here loop on each region and call a single-region fitting routine
    tab.model = NULL
    tab.fit = list()
    parBest = par
    pois = 0
    nblock = 3
    accept.vec = array(0, c(n, nblock))

    cases = mydata$model$epi
    sh = mydata$model$sh
    school = mydata$model$school
    gamaepi = mydata$model$gamaepi
    wght = rep(0, nperiods)
    wght[1:nperiodsFit] = 1


    # In the case of model 4 - the minimum value for R0 should be 1
    if (mydata$imodel == 4) {
        fit_pmin[3, 1:n] = 1
        model_pmin[3] = 1
    }

    mypar = model_par
    ndays = tps[nperiods]
    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma), scales = as.double(c(n,
            nperiodsFit)), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays))

    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))


	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

    # replace some values

    fit_par[1, 1:n] = round(runif(n, min = t0[1] + 1, max = t0[1] + cadence * 1))  #t0
    fit_par[2, 1:n] = runif(n, min = 0.02, max = 0.1)  #pC
    fit_par[3, 1:n] = runif(n, min = 1.2, max = 1.3)  #R0
    for (j in 1:n) fit_par[8, j] = fit_par[8, j] = runif(1, fit_par[8, j] * 0.75, fit_par[8, j] * 1.5)


    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")
    for (iregion in 1:n) {

        cat("\n Direct Uncoupled Fitting of: ", mydata$fit$name[iregion], "\n\n")

        # Set the seed

        iseed = set.iseed()

        mypar = fit_par[, iregion]
        cases = mydata$fit$epi[, iregion]
        sh = mydata$fit$sh[, iregion]
        school = mydata$fit$school[, iregion]
        gamaepi = mydata$fit$gamaepi[, iregion]
        ndays = tps[nperiods]

        solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
            wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin[, iregion]), parmax = as.double(fit_pmax[,
                iregion]), dx = as.double(fit_dx[, iregion]), ilog = as.integer(logvec), ymu = as.double(ymu), sigma = as.double(sigma),
            scales = as.double(c(1, 1)), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC), ithin = as.integer(ithin),
            nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
            rtn = as.double(rtn[, iregion]), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab), 				
            imodel = as.integer(mydata$imodel), ndays = as.integer(ndays))

        rtn[, iregion] = solution$rtn
        tab.tmp = matrix(solution$tab, ncol = (nparam + 1))

		# Convert LLK to AICc
		
		colnames(tab.tmp) = c(myColName, "AICc")
		tab.tmp[, "AICc"] <- 2 * tab.tmp[, "AICc"] + 2 * nopt
		tab.tmp[, "AICc"] <- tab.tmp[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)
		tab.fit[[iregion]] = tab.tmp
		
        accept.vec[iregion, 1:nblock] = solution$accept

    }
    cat("\n Uncoupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed\n")
    # # write the MCMC statistics

    success = mcmc.single.write(tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list, run.list = run.list, mydata = mydata,
        imask = imask, ireal = ireal)

    # # calculate the error in the CAR for each region

    ##carErr = calcCARErr(rtn = rtn, mydata = mydata)

    # # call the routine that will plot the results # This routine will also calculate randomly chosen profiles
    nRnd = 1000

    profile = array(data = 0, dim = c(nRnd, nperiods, n))
    for (i in 1:n) {
        if (!is.null(tab.fit[[i]])) {
            profile[, , i] <- calcRandomProfiles(sh = mydata$fit$sh[, i], school = mydata$fit$school[, i], pop = mydata$fit$pop[i], tab = tab.fit[[i]],
                tps = tps, nRnd = nRnd, cadence = mydata$cadence)
        }
    }
    if (!is.null(tab.model)) {
        model_profile <- calcRandomProfiles(sh = mydata$model$sh, school = mydata$model$school, pop = mydata$model$pop, tab = tab.model,
            tps = tps, nRnd = nRnd, cadence = mydata$cadence)
    }

    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list)

    success = saveRData(mydata = mydata, all_years_epi = NULL, run.input = run.input, ireal = ireal)

    success = writeCSV(mydata = mydata, run.list = run.list, model_rtn = model_rtn, model_profile = model_profile,
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

    ## create a graphic files that compares the input values for pC and R0 with their posterior distribution

    success = synthetic.compare(tab.fit = tab.fit, mydata = mydata, fit_par = save_fit_par, opt.list = opt.list, run.list = run.list,
        imask = imask, ireal = ireal)

    list(fit_rtn = rtn, model_rtn = model_rtn, fit_profile = profile, model_profile = model_profile, tab.model = tab.model, tab.fit = tab.fit)

}



runSyntheticMulti <- function(mydata = NULL, opt.list = NULL, run.list = NULL, ireal = 1, iseed = NULL) {

    #' Driver Routine Coupled Spatial Model - Synthetic Data
    #'
    #' A spatailly coupled MCMC fit of synthetic data. The code first fits the model data directly and then
    #'   fits it using a weighted coupled model of the regions. The data is synthetic (although it has the dataType 'cdc').
    #'   This synthetic example has been set up for the case of fit_level = 3 and mod_level = 2, that is fit the national data using the ten
    #'   coupled HHS regions
    #' @param mydata A complex list with all the data for a single  season. Everything except the number of cases and \%ILI
    #'  is from the 'cdc' dataType and will not be replaced.  Only the number of cases and the \%ILI will be replaced
    #'  by syntethic profiles
    #' @param opt.list A Logical list with TRUE or FALSE values for all the parameters supported by \pkg{DICE}.
    #'  The synthetic example has been set-up for model = 4
    #' @param run.list A list with parameters needed for the MCMC procedure
    #' @param ireal Integer - the MCMC chain number.  Default is 1.
    #' @param iseed An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly. Setting the seed to a known integer allows an MCMC chain to be reproducible.
    #' @return A list with the following arguments:
    #' \describe{
    #' \item{model_rtn}{The best result of the MCMC procedure for the model synthetic data}
    #' \item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model synthetic data}
    #' \item{tab.model}{The MCMC history of the direct fit to the model synthetic data}
    #' \item{fit_rtn}{The best result for fitting the model synthetic data using the coupled synthetic fit regions}
    #' \item{fit_profile}{Randomly selected results from the MCMC procedure of fitting the model data using the coupled fit regions}
    #' \item{tab.fit}{The MCMC history of fitting the model data using the coupled fit regions}
    #' }
    #' @examples
    #' runSyntheticMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = 1}


    # Take only the first list from opt.list, the second is optional and is needed only for a coupled run
    opt.cpl = opt.list[[2]]
    opt.list = opt.list[[1]]
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    prior = mydata$prior
    season = mydata$season
    FY = mydata$FY
    Temp = mydata$Temp

    cadence = mydata$cadence

    # if iseed is not specified, create it. Else generate a unique, reproducible seed for each
    # realization.
    if (is.null(iseed)) {
        iseed = set.iseed()
    } else {
        iseed = as.integer((iseed * ireal)%%.Machine$integer.max)
    }
    # Here we set the random number generation seed for R
    set.seed(iseed)


    out <- setupODEforSingleCDC(mydata = mydata)

    tps = out$tps

    n = mydata$fit$nregions
    rtn = out$rtn

    time = out$t0

    t0 = rep(out$t0, n)

    ## Set up list of parameter names

    par_names <- set.param.list(epi_model = mydata$epi_model)

    par_names_cpl <- set.param.cpl.list()

    setup = setup.model.mcmc(mydata = mydata, opt.list = opt.list, run.list = run.list,
        tps = out$tps, par_names = par_names)

    model_pmin = setup$parmin
    model_pmax = setup$parmax

    model_dx = setup$pardx

    ## These random values in model_par wil be used as the initial guess for teh direct fit of the model fata

    model_par = setup$par


    # these are the same for both model and fit

    nparam = setup$nparam
    nopt = setup$nopt
    logbase = setup$logbase
    logvec = setup$logvec
    tab.model = setup$tab
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
    # replace some values - and use them to generate the synthetic profiles
    fit_par[1, 1:n] = round(runif(n, min = t0[1] + 1, max = t0[1] + cadence * 1))  #t0
    fit_par[2, 1:n] = runif(n, min = 0.02, max = 0.1)  #pC
    fit_par[3, 1:n] = runif(n, min = 1.2, max = 1.3)  #R0

    for (j in 1:n) fit_par[8, j] = min(mydata$fit$epi[1, j], mydata$fit$epi[nperiods, j])

    save_fit_par = fit_par
    colnames(save_fit_par) = mydata$fit$name
    rownames(save_fit_par) = names(opt.list)


    # set up initial guess / min/max values etc for the coupling elements - there are three of them

    setup = setup.coupling.mcmc(par_names_cpl = par_names_cpl, opt.list = opt.cpl)

    cpl_pmin = setup$parmin
    cpl_pmax = setup$parmax
    cpl_dx = setup$pardx
    cpl_par = setup$par
    cpl_nparam = setup$nparam
    cpl_nopt = setup$nopt
    cpl_logvec = setup$logvec
    cpl_par["sd"] = 200
    cpl_par["gamma"] = 4

    save_cpl_par = cpl_par

    ## Build the Rij matrix - there are zero's on the diagonal

    Rij = distance.matrix(mydata = mydata)

    cpl_imask = set.cpl.imask(nparam = cpl_nparam, opt.list = opt.cpl)

    imask = set.imask(par_names = par_names, opt.list = opt.list)

    # The number of paramters that are optimized
    nopt = length(which(imask == 1))

    ## generate a synthetic profile

    fit_epi = array(0, c(nperiods, n))

    mypar = fit_par
    sh = as.matrix(mydata$fit$sh)
    school = as.matrix(mydata$fit$school)
    pop = mydata$fit$pop
    coef = mydata$fit$coef

    myparCPL = cpl_par

    solution = .Fortran("genmulti", n = as.integer(n), sh = as.double(sh), school = as.double(school), nparam = as.integer(nparam), par = as.double(mypar),
        parCPL = as.double(myparCPL), Rij = as.double(Rij), nperiods = as.integer(nperiods),
        tps = as.double(tps[1:nperiods]), rtn = as.double(array(data = 0, dim = c(nperiods, n))), ndays = as.integer(tps[nperiods]))


    fit_epi = array(0, c(nperiods, n))
    fit_epi = matrix(solution$rtn, ncol = n)


    # Build the national number of cases as the sum of the HSH regions
    model_epi = rep(0, nperiods)
    for (i in 1:nperiods) model_epi[i] = sum(fit_epi[i, 1:n])


    ## Create %ILI for model and fit using cases and factor
    fit_ili = fit_epi
    model_ili = model_epi
    for (i in 1:nperiods) {
        for (j in 1:n) fit_ili[i, j] = fit_ili[i, j]/as.numeric(mydata$fit$factor[j])
        model_ili[i] = model_ili[i]/mydata$model$factor
    }

    model_epi = round(model_epi)
    fit_epi = round(fit_epi)
    fit_gamma = gammaEpi(epi = fit_epi)
    model_gamma = lgamma((model_epi + 1))

    synData = mydata
    synData$model$raw = model_ili
    synData$model$epi = model_epi
    synData$model$gamaepi = model_gamma

    synData$fit$raw = fit_ili
    synData$fit$epi = fit_epi
    synData$fit$gamaepi = fit_gamma

    saveData = mydata
    mydata = synData

    # Here loop on each region and call a single-region fitting routine

    parBest = par
    pois = 0
    nblock = 3
    accept.vec = array(0, c(n, nblock))

    cases = mydata$model$epi
    sh = mydata$model$sh
    school = mydata$model$school
    gamaepi = mydata$model$gamaepi
    mypar = model_par

    ## These are not used in the synthetic case but we need to allocate the array
    model_ymu = rep(0, nparam)
    model_sigma = rep(0, nparam)
    wght = rep(0, nperiods)
    wght[1:nperiodsFit] = 1
    ndays = (nperiods + 2) * 7

 	nRnd = 1000
	model_profile = array(0,c(nRnd, nperiods))
	profile       = array(0,c(nRnd, nperiods, n))


    cat("\n ****** Fitting ", mydata$model$name, " Directly ****** \n")

    solution = .Fortran("episingle", epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.doubel(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(model_pmin), parmax = as.double(model_pmax),
        dx = as.double(model_dx), ilog = as.integer(logvec), ymu = as.double(model_ymu), sigma = as.double(model_sigma), Temp = as.double(c(n, nperiodsFit)),
        nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rep(0, nperiods)), accept = as.double(rep(0, nblock)), pois = as.double(pois), tab = as.single(tab.model),
        imodel = as.integer(mydata$imodel), ndays = as.integer(ndays), profile = as.single(model_profile), nRnd = as.integer(nRnd),epi_model=as.integer(1))

    cat("\n ****** Direct Fitting of ", mydata$model$name, " Completed ****** \n")

    model_rtn = solution$rtn
    tab.model = array(solution$tab, c(nlines, (nparam + 1)))
    model_profile = array(solution$profile,c(nRnd, nperiods))

	## Convert LLK to AICc
 	myColName = names(opt.list)
    colnames(tab.model) = c(myColName, "AICc")
    tab.model[, "AICc"] <- 2 * tab.model[, "AICc"] + 2 * nopt
    tab.model[, "AICc"] <- tab.model[, "AICc"] + (2 * nopt * (nopt + 1))/(nperiods - nopt - 1)

    #if (!is.null(tab.model)) {
    #    model_profile <- calcRandomProfiles(sh = mydata$model$sh, school = mydata$model$school, pop = mydata$model$pop, tab = tab.model,
    #        tps = tps, nRnd = nRnd, cadence = mydata$cadence)
    #}

    ## Now fit the coupled regions

    cases = as.matrix(mydata$fit$epi)
    sh = as.matrix(mydata$fit$sh)
    school = as.matrix(mydata$fit$school)
    gamaepi = as.matrix(mydata$fit$gamaepi)
    mypar = fit_par
    tab = array(0, c(nlines, (nparam + 1), n))
    tab.cpl = array(0, c(nlines, 2))


    ## Randomize the initial guess for the parameters we will be fitting


    cpl_par[1] = runif(1, 100, 300)
    cpl_par[2] = runif(1, cpl_par[2] - 1, cpl_par[2] + 1)

    # replace some values

    fit_par[1, 1:n] = round(runif(n, min = t0[1] + 1, max = t0[1] + cadence * 1))  #t0
    fit_par[2, 1:n] = runif(n, min = 0.02, max = 0.1)  #pC
    fit_par[3, 1:n] = runif(n, min = 1.2, max = 1.3)  #R0
    for (j in 1:n) fit_par[8, j] = fit_par[8, j] = runif(1, fit_par[8, j] * 0.75, fit_par[8, j] * 1.5)

    mypar = fit_par

    # weights - just one for all the weeks that are being fitted

    wght = array(0, c(nperiods, n))
    wght[1:nperiodsFit, 1:n] = 1



    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data\n")

    solution = .Fortran("epimulti", n = as.integer(n), epi = as.double(cases), gamaepi = as.double(gamaepi), sh = as.double(sh), school = as.double(school),
        wght = as.double(wght), nparam = as.integer(nparam), par = as.double(mypar), parmin = as.double(fit_pmin), parmax = as.double(fit_pmax),
        dx = as.double(fit_dx), ilog = as.integer(logvec), imask = as.integer(imask), iseed = as.integer(iseed), nsamps = as.integer(nMCMC),
        ithin = as.integer(ithin),nperiods = as.integer(nperiods), tps = as.double(tps[1:nperiods]),
        rtn = as.double(rtn), pois = as.double(rep(0, n)), 
        coef = as.double(mydata$fit$coef), parCPL = as.double(cpl_par), pminCPL = as.double(cpl_pmin),
        pmaxCPL = as.double(cpl_pmax), stepCPL = as.double(cpl_dx), ilogCPL = as.integer(cpl_logvec), imaskCPL = as.integer(cpl_imask),
        Rij = as.double(Rij), tab = as.single(tab), tabCPL = as.single(tab.cpl),  profiles = as.single(profile), nRnd = as.integer(nRnd), ndays = as.integer(ndays), epi_model=as.integer(1))

    cat("\n Coupled Indirect Fitting of ", mydata$model$name, " Using ", mydata$fitLevelName, " Data Completed \n")
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

    success = mcmc.multi.write(tab.model = tab.model, tab.fit = tab.fit, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list,
        mydata = mydata, imask = imask, cpl_imask = cpl_imask, ireal = ireal)

    # # call the routine that will plot the results # This routine will also calculate randomly chosen profiles

    #profile <- calcCPLRandomProfiles(mydata = mydata, tab = tab, tab.cpl = tab.cpl, Rij = Rij, tps = tps, nRnd = nRnd)

    # # save a binary RData file of the input of this run

    run.input = list(run.list = run.list, opt.list = opt.list, opt.cpl = opt.cpl)

    success = saveRData(mydata = mydata, all_years_epi = NULL, run.input = run.input, ireal = ireal)

    success = writeCSV(mydata = mydata, run.list = run.list, model_rtn = model_rtn, model_profile = model_profile,
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

    success = synthetic.compare(tab.fit = tab.fit, tab.cpl = tab.cpl, mydata = mydata, fit_par = save_fit_par, cpl_par = save_cpl_par, opt.list = opt.list,
        run.list = run.list, imask = imask, ireal = ireal)


    list(model_rtn = model_rtn, model_profile = model_profile, tab.model = tab.model, rtn = rtn, profile = profile, tab = tab)


}


synthetic.compare <- function(tab.fit = NULL, tab.cpl = NULL, fit_par = NULL, cpl_par = NULL, mydata = mydata, opt.list = NULL, run.list = NULL,
    imask = NULL, ireal = 1) {

    #' Plot Posterior Distributions for a Synthetic Example
    #' For each region we plot in blue the posterior density of R0 and of pC.  In red we show the input value
    #' The entire history of the chain is used when showing the density.
    #' In the case of a coupled model we also compare the values for the saturation distance and the
    #' distance power (s_d and gamma respectively).
    #' @param tab.fit  The MCMC history of the fit
    #' @param tab.cpl MCMC history of the coupling parameters (if present)
    #' @param fit_par Parameter values used to create the synthetic profiles
    #' @param cpl_par Coupling parameter values used to create the synthetic profile
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param mydata a list with all the data available for this run
    #' @param imask An array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param ireal Integer - the MCMC chain number
    #' @return err   Returns \eqn{err = 0}
    #' @examples
    #' synthetic.compare{tab.fit = tab.fit, tab.cpl = tab.cpl, fit_par = save_fit_par, cpl_par = save_cpl_par,opt.list = opt.list,
    #' run.list = run.list, mydata = mydata, imask = imask, ireal = 1}

    myName = mydata$dataName
    nperiodsFit = mydata$nperiodsFit
    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC
    device = run.list$device

    myColName = names(opt.list)
    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]


    n.fit = mydata$fit$nregions
    iburn <- nlines/2

    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName

    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/syn-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n   For a plot Comparing Posterior Distributions to Input Parameters See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste(subDir, "/syn-", myName, "-", nperiodsFit, "-", ireal, ".png", sep = "")
        cat("\n   For a plot Comparing Posterior Distributions to Input Parameters See: ", pngName, "\n\n")
        png(file = pngName, width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }

    if (mydata$imodel == 1) {


        nc = n.fit/2
        nr = 4
        var2plot = c("R0", "pC")

        par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

        for (ivar in 1:length(var2plot)) {
            myvar = var2plot[ivar]
            irow = which(rownames(fit_par) == myvar)
            for (i in 1:n.fit) {
                tab = tab.fit[[i]]
                colnames(tab) = c(myColName, "llk")
                icol = which(colnames(tab) == myvar)
                tab <- mcmc(data = tab, start = (ithin * iburn), end = nMCMC, thin = ithin)
                densplot(tab[, icol], ylab = "Frequency", xlab = myvar, main = mydata$fit$name[i], col = "blue", show.obs = FALSE)
                abline(v = fit_par[irow, i], col = "red", lwd = 2)
            }

        }

    } else {

        nc = n.fit/2
        nr = 5
        var2plot = c("R0", "pC")

        par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

        for (ivar in 1:length(var2plot)) {
            myvar = var2plot[ivar]
            irow = which(rownames(fit_par) == myvar)
            for (i in 1:n.fit) {
                tab = tab.fit[, , i]
                colnames(tab) = c(myColName, "llk")
                icol = which(colnames(tab) == myvar)
                tab <- mcmc(data = tab, start = (ithin * iburn), end = nMCMC, thin = ithin)
                densplot(tab[, icol], ylab = "Frequency", xlab = myvar, main = mydata$fit$name[i], col = "blue", show.obs = FALSE)
                abline(v = fit_par[irow, i], col = "red", lwd = 2)
            }

        }

        var2plotCPL = names(cpl_par)
        title = c(expression("s"[d]), expression(gamma))
        for (ivar in 1:length(var2plot)) {
            myvar = var2plotCPL[ivar]
            icol = ivar
            tab <- mcmc(data = tab.cpl, start = (ithin * iburn), end = nMCMC, thin = ithin)
            densplot(tab[, icol], ylab = "Frequency", xlab = myvar, main = title[ivar], xlim = range(tab[, icol], cpl_par[ivar]), col = "blue",
                show.obs = FALSE)
            abline(v = cpl_par[ivar], col = "red", lwd = 2)
        }


    }
    if (tolower(device) == "pdf" || tolower(device) == "png")
        dev.off()

    err = 0
    return(err)
}


