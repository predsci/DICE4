##
## Functions Related to writing or Plotting the MCMC output
##

mcmc.onepatch.write <- function(tab.model = NULL, opt.list = NULL, run.list = NULL, mydata = NULL, imask = NULL, ireal = 1) {

    #' Write the MCMC History of a Run - Single Region
    #'
    #' Writes an RData file with the MCMC history for a \pkg{DICE} run on a single region/patch.
    #'   The function also calcualtes and prints to the screen the statistics for all the
    #'   parameters that were optimized and does a Gaussian fit to these parameters.
    #'   The results of these fits are written to separate a csv file.
    #' @param tab.model The MCMC history of the direct fit of the mydata
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param mydata A dataframe with all the data available for this run
    #' @param imask An array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param ireal Integer - the MCMC chain number
    #' @return err   Returns \eqn{err = 0}
    #' @examples
    #' mcmc.onepatch.write{tab.model = tab.model, opt.list = opt.list,
    #' run.list = run.list, mydata = mydata, imask = imask, ireal = 1}

    myName = mydata$dataName
	myName = gsub(" ","",myName)
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    weeks = mydata$weeks
    epi_model = mydata$model$epi
    nperiodsFit = mydata$nperiodsFit

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    myColName = names(opt.list)
    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]

    # how many steps to burn - here we set it to half which is very rigid
    iburn <- nlines/5

    # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics


    if (!is.null(tab.model))
        results = mcmc(data = tab.model, start = (ithin * iburn), end = nMCMC, thin = ithin)

    # check to see if 'mydata' sub-directory exists, if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)

    # here write an R mydata file print the chains statistics to the screen-only for optimized variables
    mcmc.summary.list = list()

    if (!is.null(tab.model)) {
        print(summary(results[,c(names.opt,"AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
        mcmc.summary.list[[1]] <- summary(results[,c(names.opt,"AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE))
    }


    filename <- paste(subDir, "/mcmc-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    cat("\n Writing R object Data file for this Chain: ", filename, "\n")
    
    # save all the information needed for this function 
    dump = list()
    dump = list(mydata = mydata, tab.model = tab.model, tab.fit = NULL, tab.cpl = NULL, opt.list = opt.list, opt.cpl = NULL, run.list = run.list, imask = imask, cpl_imask = NULL, ireal = ireal, idevice = 1)
    save(dump, file = filename)


    opt.normal = array(0, c(((length(names.opt) + 1)), 2))

    rownames(opt.normal) = c(names.opt, "AICc")
    colnames(opt.normal) = c("mean", "sd")
    # for pC and R0 - do a maximum LLK fit to a normal coefficient
    op = options(digits = 3)
    k = 0
    sumAICc = 0

    if (!is.null(tab.model)) {
        tab = tab.model
        colnames(tab) = c(myColName, "AICc")
        results = tab[iburn:nlines, ]
        ncol = dim(results)[2]
        for (j in 1:(ncol - 1)) {
            if (imask[j] != 1)
                next
            k = k + 1
            z = fitdistr(results[, j], "normal")
            opt.normal[k, 1] = z$estimate["mean"]
            opt.normal[k, 2] = z$estimate["sd"]
        }
        llk = results[, ncol]

        k = k + 1
        aicc = 2 * llk + 2 * (nopt) + (2 * (nopt) * ((nopt) + 1))/(nperiodsFit - (nopt) - 1)


        opt.normal[k, 1] = mean(aicc)
        opt.normal[k, 2] = sd(aicc)

    }

    filename = paste(subDir, "/gaussian-fit-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(opt.normal, file = filename)

    success = 0
    return(success)

}


mcmc.onepatch.plot <- function(tab.model = NULL, opt.list = NULL, run.list = NULL, mydata = NULL, imask = NULL, ireal = 1, idevice = 1) {

    #' Plot the Posterior Distribution of the Parameters - Single Region
    #'
    #' Plots the posterior distribution of the MCMC parameters from a  \pkg{DICE} run on a single region/patch.
    #' @param tab.model The MCMC history of the direct fit of the data
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param mydata A dataframe with all the data available for this run
    #' @param imask An array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param ireal Integer - the MCMC chain number
    #' @return err   Returns \eqn{err = 0}
    #' @examples
    #' mcmc.onepatch.plot{tab.model = tab.model, opt.list = opt.list,
    #' run.list = run.list, mydata = mydata, imask = imask, ireal = 1}

    myName = mydata$dataName
	myName = gsub(" ","",myName)
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    weeks = mydata$weeks
    epi_model = mydata$model$epi

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC
    device = run.list$device[idevice]

    myColName = names(opt.list)
    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]
	nparam <- dim(tab.model)[2] - 1
	
    # how many steps to burn - set to 1/5
    iburn <- nlines/5

    # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics


    if (!is.null(tab.model))
        results = mcmc(data = tab.model, start = (ithin * iburn), end = nMCMC, thin = ithin)

    # check to see if 'mydata' sub-directory exists, if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)


    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/posterior-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n   For a Plot of the Posterior Distributions See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste(subDir, "/posterior-", myName, "-", nperiodsFit, "-", ireal, ".png", sep = "")
        cat("\n   For a Plot of the Posterior Distributions See: ", pngName, "\n\n")
        png(file = pngName, width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }


    nc = 4
    nr = 4
    var2plot = names.opt

    par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

    results = tab.model[iburn:nlines, ]
    legend = rep(0, 2)

    for (ivar in 1:length(var2plot)) {
        myvar = var2plot[ivar]
        icol = which(colnames(tab.model) == myvar)
        tab <- mcmc(data = tab.model, start = (ithin * iburn), end = nMCMC, thin = ithin)
        densplot(tab[, icol], ylab = "Frequency", xlab = myvar, main = myvar, col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

        z = fitdistr(results[, icol], "normal")
        val.mean = z$estimate["mean"]
        val.sd = z$estimate["sd"]
        legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
        legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
        legend("topleft", legend = legend, bty = "n")
    }

    if (tolower(device) == "pdf" || tolower(device[1]) == "png")
        dev.off()

    success = 0
    return(success)

}


	mcmc.single.write <- function(tab.model = NULL, tab.fit = NULL, opt.list = NULL, run.list = NULL, mydata = NULL, imask = NULL, ireal = 1) {

    #' Write the MCMC History of an Uncoupled \pkg{DICE} Run
    #'
    #' Writes an RData file with the MCMC history for an uncoupled \pkg{DICE} run.
    #'   The function also calculates and prints to the screen the statistics for all the
    #'   parameters that were optimized and does a Gaussian fit to these parameters.
    #'   The results of these fits are written to separate a csv file.
    #' @param tab.model The MCMC history of the direct fit of the model data
    #' @param tab.fit  The MCMC history of an indirect fit of the model using the fit data
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param mydata A dataframe with all the data available for this run
    #' @param imask An array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param ireal Integer - the MCMC chain number
    #' @return err   Returns \eqn{err = 0}
    #' @examples
    #' mcmc.single.write{tab.model = tab.model, tab.fit = tab.fit, opt.list = opt.list,
    #' run.list = run.list, mydata = mydata, imask = imask, ireal = 1}

    myName = mydata$dataName
	myName = gsub(" ","",myName)
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    weeks = mydata$weeks
    epi_fit = mydata$fit$epi
    epi_model = mydata$model$epi

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    myColName = names(opt.list)
    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]

    n.fit = mydata$fit$nregions
    n.fit1 = n.fit + 1


	nparam <- dim(tab.fit[[1]])[2] - 1
	nparam1 = nparam + 1
	
    # how many steps to burn - set to 1/5 here
    iburn <- nlines/5

    # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics
    results.list = list()

    for (i in 1:n.fit) {
        tab = tab.fit[[i]]
        results <- mcmc(data = tab, start = (ithin * iburn), end = nMCMC, thin = ithin)
        results.list[[i]] = results
    }
    if (!is.null(tab.model))
        results.list[[n.fit1]] = mcmc(data = tab.model, start = (ithin * iburn), end = nMCMC, thin = ithin)

    # check to see if 'mydata' sub-directory exists, if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)


    # print the chains statistics to the screen-only for optimized variables and AICc

    point = which(imask == 1)
    point = append(point, nparam1)
    for (i in 1:n.fit) {
        tmp = results.list[[i]]
        colnames(tmp) = c(myColName, "AICc")
        print(summary(tmp[,c(names.opt,"AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
    }


    if (!is.null(tab.model)) {
        tmp = results.list[[i]]
        colnames(tmp) = c(myColName, "AICc")
        print(summary(tmp[,c(names.opt,"AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
    }


    # save the complete chain here
    filename <- paste(subDir, "/mcmc-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    cat("\n Writing R object Data file for this Chain: ", filename, "\n")

    dump = list()
    dump = list(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, tab.cpl = NULL, opt.list = opt.list, opt.cpl = NULL, run.list = run.list, imask = imask, cpl_imask = NULL, ireal = ireal, idevice = 1)
    save(dump, file = filename)


    opt.normal = array(0, c(((length(names.opt) + 1) * n.fit1), 2))

    rownames(opt.normal) = c(rep(c(names.opt, "AICc"), n.fit1))
    colnames(opt.normal) = c("mean", "sd")
    # for pC and R0 - do a maximum LLK fit to a normal coefficient
    op = options(digits = 3)
    k = 0
    sumAICc = 0


    for (i in 1:n.fit) {
        tab = tab.fit[[i]]
        colnames(tab) = c(myColName, "AICc")
        results = tab[iburn:nlines, ]
        ncol = dim(results)[2]
        for (j in 1:(ncol - 1)) {
            if (imask[j] != 1)
                next
            k = k + 1
            z = fitdistr(results[, j], "normal")
            opt.normal[k, 1] = z$estimate["mean"]
            opt.normal[k, 2] = z$estimate["sd"]
        }
        llk = results[, ncol]

        k = k + 1
        aicc = 2 * llk + 2 * (nopt) + (2 * (nopt) * ((nopt) + 1))/(nperiodsFit - (nopt) - 1)


        opt.normal[k, 1] = mean(aicc)
        opt.normal[k, 2] = sd(aicc)


    }

    if (!is.null(tab.model)) {
        tab = tab.model
        colnames(tab) = c(myColName, "AICc")
        results = tab[iburn:nlines, ]
        ncol = dim(results)[2]
        for (j in 1:(ncol - 1)) {
            if (imask[j] != 1)
                next
            k = k + 1
            z = fitdistr(results[, j], "normal")
            opt.normal[k, 1] = z$estimate["mean"]
            opt.normal[k, 2] = z$estimate["sd"]
        }
        llk = results[, ncol]

        k = k + 1
        aicc = 2 * llk + 2 * (nopt) + (2 * (nopt) * ((nopt) + 1))/(nperiodsFit - (nopt) - 1)


        opt.normal[k, 1] = mean(aicc)
        opt.normal[k, 2] = sd(aicc)

    }

    filename = paste(subDir, "/gaussian-fit-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(opt.normal, file = filename)

    success = 0
    return(success)

}

mcmc.multi.write <- function(tab.model = NULL, tab.fit = NULL, tab.cpl = NULL, opt.list = NULL, opt.cpl = NULL, run.list = NULL, mydata = NULL,
    imask = NULL, cpl_imask = NULL, ireal = 1) {

    #' Write the MCMC History of a Coupled \pkg{DICE} Run
    #'
    #' Writes an RData file with the MCMC history for an coupled \pkg{DICE} run.
    #'   The function also calculates and prints to the screen the statistics for all the
    #'   parameters that were optimized and does a Gaussian fit to these parameters.
    #'   The results of these fits are written to separate a csv file.
    #' @param tab.model The MCMC history of the direct fit of the model data
    #' @param tab.fit  The MCMC history of an indirect fit of the model using a coupled model.
    #'   This  array includes all the parameters except the two that define the coupling matrix.
    #' @param tab.cpl The MCMC history of the two parameters that help define the coupling matrix:
    #'   the saturation distance and the distance power.
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param opt.cpl A logical list for the two parameters that help define the coupling matrix.
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param mydata A dataframe with all the data available for this run
    #' @param imask AN array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param cpl_imask An array of length two with +1/-1 for the two parameters that help define
    #'   the coupling matrix
    #' @param ireal Integer - the MCMC chain number
    #' @return err   Returns \eqn{err = 0}
    #' @examples
    #' mcmc.single.write{tab.model = tab.model, tab= tab , tab.cpl = tab.cpl,
    #' opt.list = opt.list,run.list = run.list, mydata = mydata, imask = imask,
    #' cpl_imask = cpl_imask, ireal = 1}


    myName = mydata$dataName
	myName = gsub(" ","",myName)
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$weeksData
    weeks = mydata$weeks
    epi_fit = mydata$fit$epi
    epi_model = mydata$model$epi

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    myColName = names(opt.list)
    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]
	nparam <- dim(tab.model)[2] - 1
	nparam1 = nparam + 1
	
    n.fit = mydata$fit$nregions
    n.fit1 = n.fit + 1

    ##
    myColNameCPL = names(opt.cpl)
    noptCPL = length(which(cpl_imask == 1))
    names.optCPL = myColNameCPL[which(cpl_imask == 1)]

    
    # how many steps to burn - set to 1/5
    iburn <- nlines/5

    # This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics
    results.list = list()

    point = which(imask == 1)
    point = append(point, nparam1)
    pointCPL = which(cpl_imask == 1)
    pointFit = append(point, (pointCPL + nparam1))
    for (i in 1:n.fit) {
        tmp = cbind(tab.fit[[i]], tab.cpl)
        colnames(tmp) = c(myColName, "AICc", myColNameCPL)
        results <- mcmc(data = tmp, start = (ithin * iburn), end = nMCMC, thin = ithin)
        results.list[[i]] = results
    }


    if (!is.null(tab.model))
        results.list[[n.fit1]] = mcmc(data = tab.model, start = (ithin * iburn), end = nMCMC, thin = ithin)

    # check to see if 'mydata' sub-directory exists, if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)

    # here write an R mydata file print the chains statistics to the screen-only for optimized variables and AICc

    for (i in 1:n.fit) {
        tmp = results.list[[i]]
        colnames(tmp) = c(myColName, "AICc", myColNameCPL)
        print(summary(tmp[,c(names.opt, "AICc", myColNameCPL)], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
    }


    if (!is.null(tab.model)) {
        tmp = results.list[[n.fit1]]
         colnames(tmp) = c(myColName, "AICc")
        print(summary(tmp[,c(names.opt, "AICc")], quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95, na.rm = TRUE), digits = 2))
    }


    filename <- paste(subDir, "/mcmc-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    cat("\n Writing R object Data file for this Chain: ", filename, "\n")

	dump = list()
    dump = list(mydata = mydata, tab.model = tab.model, tab.fit = tab.fit, tab.cpl = tab.cpl, opt.list = opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, cpl_imask = cpl_imask, ireal = ireal, idevice = 1)
    
    save(dump, file = filename)

    opt.normal = array(0, c(((length(names.opt) + 1) * n.fit1), 2))

    rownames(opt.normal) = c(rep(c(names.opt, "AICc"), n.fit1))
    colnames(opt.normal) = c("mean", "sd")
    # for pC and R0 - do a maximum LLK fit to a normal coefficient
    op = options(digits = 3)
    k = 0
    sumAICc = 0

    for (i in 1:n.fit) {        
        results = tab.fit[[i]]
        colnames(results) = c(myColName, "AICc")
        ncol = dim(results)[2]
        for (j in 1:(ncol - 1)) {
            if (imask[j] != 1)
                next
            k = k + 1
            z = fitdistr(results[, j], "normal")
            opt.normal[k, 1] = z$estimate["mean"]
            opt.normal[k, 2] = z$estimate["sd"]
        }
        llk = results[, ncol]

        k = k + 1
        aicc = 2 * llk + 2 * (nopt) + (2 * (nopt) * ((nopt) + 1))/(nperiodsFit - (nopt) - 1)


        opt.normal[k, 1] = mean(aicc)
        opt.normal[k, 2] = sd(aicc)


    }

    if (!is.null(tab.model)) {
        tab = tab.model
        colnames(tab) = c(myColName, "AICc")
        results = tab[iburn:nlines, ]
        ncol = dim(results)[2]
        for (j in 1:(ncol - 1)) {
            if (imask[j] != 1)
                next
            k = k + 1
            z = fitdistr(results[, j], "normal")
            opt.normal[k, 1] = z$estimate["mean"]
            opt.normal[k, 2] = z$estimate["sd"]
        }
        llk = results[, ncol]

        k = k + 1
        aicc = 2 * llk + 2 * (nopt) + (2 * (nopt) * ((nopt) + 1))/(nperiodsFit - (nopt) - 1)


        opt.normal[k, 1] = mean(aicc)
        opt.normal[k, 2] = sd(aicc)

    }

    filename = paste(subDir, "/gaussian-fit-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(opt.normal, file = filename)


    success = 0
    return(success)

}


