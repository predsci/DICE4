##
## All Plotting functions using Generic plot Command and the wrapper external function
##

plotFitOnePatch <- function(model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1, run.list = NULL, idevice = 1) {

    #' Plot the Results of a \pkg{DICE} Run - Single Region
    #'
    #' Plot the results of \pkg{DICE} run for a single region/patch. We show the observed incidence along with our fits and
    #' if appropriate predictions. We show the best result and randomly selected results from the MCMC chain.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the mydata
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotFitOnePatch{model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, device = device, idevice = 1}

    device = run.list$device[idevice]

    if (is.null(device))
        device = "png"
    if (is.null(model_profile))
        return

    FY = mydata$FY
    model = mydata$imodel
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$nperiodsData
    reg.model.name = mydata$model$name

    nRnd = dim(model_profile)[1]
    step = max(1, nRnd/100)
	irnd.set = seq(from = 1, to = nRnd, by = step)
    n.model = 1

    model_onset = mydata$model$onset
    tps = mydata$weeks

    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName
	myName = gsub(" ","",myName)
    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n\n For a Plot of the Results See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".png", sep = "")
        cat("\n\n For a Plot of the Results See: ", pngName, "\n\n")
        png(file = pngName, width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }


    # convert the fit and model results to %ILI calculate the national

    # convert the model to %ILI calculate the national

    model_factor = mydata$model$factor

    ## This model raw mydata

    model_ili = mydata$model$raw

    model_rtn_ili = NULL
    model_profile_ili = NULL

    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile
        model_profile_ili = model_profile_ili/model_factor
    }


    ## for plotting - replace zero's with NA

    if (nperiodsFit < nperiods) {
        index <- which(model_ili == 0)
        if (length(index) >= 1)
            model_ili[index] = NA

    }

    lwd = 2

    # Plot the mydata and fits/predictions

    par(mfrow = c(3, 3), mar = c(4, 4, 1, 1))

    ## Determine ylabel based on dataType
    if (mydata$data_source == "cdc" | mydata$data_source == "gft") {
        ylab = "%ILI"
    } else {
        ylab = "# Cases"
    }

    if (!is.null(model_rtn) && !is.null(model_profile)) {
        model_mean = rep(0, nperiods)
        for (iweek in 1:nperiods) model_mean[iweek] = mean(model_profile_ili[, iweek])
        ymax = max(model_rtn_ili[1:nperiodsData], model_profile_ili[, 1:nperiodsData], model_ili[1:nperiodsData], na.rm = TRUE)

        plot(1:nperiods, model_ili, type = "l", col = "black", xlim = c(1, nperiods), ylim = c(0, ymax), xaxt = "n", xlab = "EW #", ylab = ylab,
            lwd = 2, lty = 1)
        points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

        for (irnd in irnd.set) {
            lines(1:nperiodsFit, model_profile_ili[irnd, 1:nperiodsFit], col = "lightgrey", lwd = 1, lty = 1)
            lines(nperiodsFit:nperiods, model_profile_ili[irnd, nperiodsFit:nperiods], col = "lightgrey", lwd = 1, lty = 2)
        }

        lines(1:nperiodsFit, model_rtn_ili[1:nperiodsFit], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        lines(nperiodsFit:nperiods, model_rtn_ili[nperiodsFit:nperiods], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)

        lines(1:nperiodsFit, model_mean[1:nperiodsFit], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        lines(nperiodsFit:nperiods, model_mean[nperiodsFit:nperiods], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)
        lines(1:nperiods, model_ili, type = "l", col = "black", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

        if (length(model_onset) > 0)
            abline(h = model_onset, col = "grey", lty = 2, lwd = 2)

        rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
        reg.name = paste(mydata$model$name, c("-Data", "-Model-Best", "-Model-Mean", "-Model-Random"), sep = "")

        legend("topleft", c(mydata$FY, reg.name), text.col = c("black", "black", "red", "blue", "grey"), bty = "n")
        axis(1, at = 1:nperiods, labels = weeks)

        if (any(model == c(2, 3, 101, 102, 104, 105))) {
            school = mydata$model$school
            school[school == 0] = NA
            par(new = TRUE)
            plot(school, ylim = c(0, 6), xaxt = "n", yaxt = "n", xlab = "n", ylab = "n", col = "grey", type = "p", pch = 22, lwd = 4)

        }

        if (any(model == c(1, 3, 102, 103, 104, 105))) {
            sh = mydata$model$sh
            par(new = TRUE)
            plot(sh, xaxt = "n", yaxt = "n", xlab = "n", ylab = "n", col = "black", type = "l", lwd = 1)

        }

    }


    # maximum week in mydata

    dat_model_wk_max = which.max(model_ili)

    # maximum week in model

    drct_model_wk_max = rep(0, nRnd)

    if (!is.null(model_profile_ili))
        for (i in 1:nRnd) drct_model_wk_max[i] = which.max(model_profile_ili[i, ])

    if (!is.null(model_profile_ili)) {
        wk.min = round(min(dat_model_wk_max, drct_model_wk_max))
        wk.max = round(max(dat_model_wk_max, drct_model_wk_max))
    }
    wk.min = round(0.5 * wk.min)
    wk.max = round(1.5 * wk.max)
    wk.max = min(wk.max, nperiods)

    wk.min = max(1, wk.min)
    breaks = seq(from = wk.min, to = wk.max, by = 1)
    ylab = "Probability Density"
    xlab = "EW #"

    if (!is.null(model_profile_ili)) {
        hist(dat_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab,
            xaxt = "n", main = "", freq = FALSE)
        hist(drct_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = "", ylab = "",
            xaxt = "n", main = "", freq = FALSE, add = TRUE)
    }

    axis(1, at = wk.min:wk.max, labels = weeks[wk.min:wk.max])
    leg.text = c(mydata$FY)
    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    leg.text = c(leg.text, model.name, "Data", "Model")
    legend("topleft", legend = c("Observed/Predicted", "Peak Week"), text.col = "black", bty = "n")
    legend("topright", legend = leg.text, text.col = c("black", "black", rgb(0, 0, 1, 1/2), rgb(0, 1, 0, 1/2)), bty = "n")

    ## Histogram plots - of %ILI binned
    ylab = "Probability Density"
    if (mydata$data_source == "cdc" | mydata$data_source == "gft") {
        xlab = "%ILI"
    } else {
        xlab = "# Cases"
    }

    if (!is.null(model_profile_ili)) {
        for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
            min.val = 0
            if (iweek <= nperiodsData) {

                max.val = ceiling(max(model_ili[iweek], model_profile_ili[, iweek], na.rm = TRUE))
                max.val = 2 * max.val
	      if (mydata$model$factor != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
                if (!is.na(model_ili[iweek])) {
                  hist(model_ili[iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab,
                    ylab = ylab, main = "", freq = F)
                  hist(model_profile_ili[, iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = "",
                    ylab = "", main = "", freq = F, add = TRUE)
                } else {
                  hist(model_profile_ili[, iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = xlab,
                    ylab = "", main = "", freq = F)
                }

                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, "Data", "Model")
                legend("topright", legend = leg.text, text.col = c("black", "black", rgb(0, 0, 1, 1/2), rgb(0, 1, 0, 1/2)), bty = "n")
                legend("topleft", legend = c("Observed/Predicted", paste("%ILI for ", my.week, sep = "")), text.col = "black", bty = "n")
            } else {
                max.val = ceiling(max(model_profile_ili[, iweek], na.rm = TRUE))
                max.val = 2 * max.val

                breaks = seq(from = min.val, to = max.val, by = 0.5)
                hist(model_profile_ili[, iweek], breaks = breaks, xlim = range(breaks), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = xlab,
                  ylab = "", main = "", freq = F)

                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, "Model")
                legend("topright", legend = leg.text, text.col = c("black", "black", rgb(0, 1, 0, 1/2)), bty = "n")
                legend("topleft", legend = c("Predicted", paste("%ILI for ", my.week, sep = "")), text.col = "black", bty = "n")
            }
        }
    }

    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list, idevice = idevice, model_profile_ili = model_profile_ili, model_rtn_ili = model_rtn_ili)

    filename = paste(subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)

    # success = writecsvOnePatch(mydata=mydata,run.list = run.list, tab = tab, model_rtn_ili = model_rtn_ili,
    # model_profile_ili=model_profile_ili, ireal= ireal)

    err = 0
    return(err)
}

plotFitCDCPercentILI <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1, run.list = NULL, idevice = 1) {

    #' Plot the Results of a \pkg{DICE} Run
    #'
    #' Plot the results of an a coupled or uncoupled  \pkg{DICE} run. For each of the fit regions we plot the incidence of the region along with
    #' our predictions for it based on randomly selected results from the  history of the MCMC chain of each region.  Using the predictions
    #' for the fit regions we then show the results for the model region as a weighted sum of the fit regions.  The last panel
    #' shows our prediction for the model region using a direct fit to the model data.  The function also writes a binary RData file with
    #' all the profile predictions for the model and fit regions. Note that in the case of a coupled run the fit regions are never individually
    #' optimized. It is their weighted sum that is optimized, with the weights given by the relative population of each fit region.
    #' @param rtn An nweeks x nregion 2D  numeric array with the best fit for each region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotFitCDCPercentILI{ rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, device = device, idevice = 1}

    device = run.list$device[idevice]

    if (is.null(device))
        device = "png"
    if (is.null(profile))
        return
    if (is.null(rtn))
        return

    FY = mydata$FY
    model = mydata$imodel
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$nperiodsData
    reg.fit.name = mydata$fit$name
    reg.model.name = mydata$model$name

    nRnd = dim(profile)[1]
    step = max(1, nRnd/100)
    irnd.set = seq(from = 1, to = nRnd, by = step)

    n.model = 1
    nregions = mydata$fit$nregions

    fit_coef = mydata$fit$coef
    fit_onset = mydata$fit$onset
    model_coef = mydata$model$coef
    model_onset = mydata$model$onset
    tps = mydata$weeks

    ##
    ## check to see if output directory exists and if not create it
    ##

    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName
	myName = gsub(" ","",myName)
    if (tolower(device) == "pdf") {
        pdfName = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".pdf")
        cat("\n\n For a plot of %ILI See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, "_pg")
        png(file = paste0(pngName, "%01d.png"), width = 1200, height = 900)
        cat("\n\n For a plot of  %ILI  See: ", pngName, "\n\n")
    } else {
        dev.next()
        dev.new()
    }


    ##
    ##convert the fit and model results to %ILI calculate the national
    ## convert the model to %ILI calculate the national
    ##

    factor = as.numeric(mydata$fit$factor)
    model_factor = mydata$model$factor

    ## This is the fit and model %ILI mydata
    fit_ili = mydata$fit$raw
    model_ili = mydata$model$raw

    rtn_ili = rtn
    profile_ili = profile
    for (i in 1:nregions) {
        rtn_ili[, i] = rtn[, i]/factor[i]
        profile_ili[, , i] = profile[, , i]/factor[i]
    }

    model_rtn_ili = NULL
    model_profile_ili = NULL

    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile
        model_profile_ili = model_profile_ili/model_factor
    }

    ##
    ## for plotting - replace zero's with NA
    ##

	if (nperiodsFit < nperiods) {
		for (iregion in 1:nregions) {
			index <- which(fit_ili[, 1] == 0)

			if (length(index) >= 1) {
				fit_ili[index, iregion] = NA
			}

		}
		index <- which(fit_ili[, 1] == 0)
		if (length(index) >= 1) {
			model_ili[index] = NA
		}

	}

    colnames(fit_ili) = mydata$fit$attr$NAME_3
    colnames(rtn_ili) = mydata$fit$attr$NAME_3

    colvec = rainbow(nregions)
    lwd = rep(2, nregions)

    ##
    ## Plot one region at a time
    ##

    par(mfrow = c(3, 4), mar = c(4, 4, 1, 1))

    for (i in 1:nregions) {

        # ymax = max(fit_ili[1:nperiodsData, i], rtn_ili[1:nperiodsData, i], profile_ili[,1:nperiodsData , i], unlist(fit_onset[i]), na.rm = TRUE)

        ymax = max(fit_ili[1:nperiods, i], rtn_ili[1:nperiods, i], profile_ili[, 1:nperiods, i], unlist(fit_onset[i]), na.rm = TRUE)

        plot(1:nperiods, fit_ili[1:nperiods, i], type = "l", col = "black", xlim = c(1, nperiods), ylim = c(0, ymax), xaxt = "n", xlab = "EW #",
            ylab = mydata$fit$raw_units, lwd = lwd[1], lty = 1)

        for (irnd in irnd.set) {
            lines(1:nperiodsFit, profile_ili[irnd, 1:nperiodsFit, i], col = colvec[i], lwd = lwd[i], lty = 1)
            lines(nperiodsFit:nperiods, profile_ili[irnd, nperiodsFit:nperiods, i], col = colvec[i], lwd = lwd[i], lty = 2)
        }

        lines(nperiodsFit:nperiods, rtn_ili[nperiodsFit:nperiods, i], type = "l", col = colvec[i], lwd = lwd[i], lty = 2)

        lines(1:nperiodsFit, rtn_ili[1:nperiodsFit, i], type = "l", col = colvec[i], lwd = lwd[i], lty = 1)

        lines(1:nperiodsFit, fit_ili[1:nperiodsFit, i], type = "b", col = "black", pch = 20, lty = 2)

        if (length(fit_onset[i]) > 0)
            abline(h = fit_onset[i], col = "grey", lty = 2, lwd = 2)

        rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

        reg.name = reg.fit.name[i]
        rel.pop = fit_coef[i]
        rel.pop = round(rel.pop, digits = 3)
        legend("topleft", c(mydata$FY, reg.name, rel.pop), text.col = c("black", colvec[i], colvec[i]), bty = "n")

        axis(1, at = 1:nperiods, labels = weeks)

        if (any(model == c(2, 3, 101, 102, 104, 105))) {
            school = mydata$fit$school[, i]
            school[school == 0] = NA
            par(new = TRUE)
            plot(school, ylim = c(0, 6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "grey", type = "p", pch = 22, lwd = 4)

        }
        if (any(model == c(1, 3, 102, 103, 104, 105))) {
            sh = mydata$fit$sh[, i]
            par(new = TRUE)
            plot(sh, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "black", type = "l", lwd = 1)

        }
    }

    ##
    ## now do the national
    ##

    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {
        fit_model[i] = sum(rtn_ili[i, 1:nregions] * fit_coef[1:nregions])
        tmp = rep(0, nregions)
        for (k in 1:nregions) tmp[k] = mean(profile_ili[, i, k])
        fit_model_mean[i] = sum(tmp[1:nregions] * fit_coef[1:nregions])
        for (irnd in 1:nRnd) {
            for (k in 1:nregions) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile_ili[irnd, i, k] * fit_coef[k]
        }

    }

    ymax = max(fit_model[1:nperiods], fit_model_mean[1:nperiods], fit_model_profile[1:nperiods], model_ili[1:nperiods], na.rm = TRUE)

    plot(1:nperiods, model_ili, type = "l", col = "black", xlim = c(1, nperiods), ylim = c(0, ymax), xaxt = "n", xlab = "EW #", ylab = mydata$model$raw_units,
        lwd = 2, lty = 1)
    points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

    for (irnd in irnd.set) {
        lines(1:nperiodsFit, fit_model_profile[irnd, 1:nperiodsFit], col = "lightgrey", lwd = 3, lty = 1)
        lines(nperiodsFit:nperiods, fit_model_profile[irnd, nperiodsFit:nperiods], col = "lightgrey", lwd = 3, lty = 2)
    }

    lines(1:nperiodsFit, fit_model[1:nperiodsFit], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
    lines(nperiodsFit:nperiods, fit_model[nperiodsFit:nperiods], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)

    lines(1:nperiodsFit, fit_model_mean[1:nperiodsFit], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
    lines(nperiodsFit:nperiods, fit_model_mean[nperiodsFit:nperiods], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)
    lines(1:nperiods, model_ili, type = "l", col = "black", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
    points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

    if (length(model_onset) > 0)
        abline(h = model_onset, col = "grey", lty = 2, lwd = 2)

    rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
    reg.name = paste(mydata$model$name, c("-Data", "-Best", "-Mean", "-Random"), sep = "")

    legend("topleft", c(mydata$FY, reg.name), text.col = c("black", "black", "red", "blue", "grey"), bty = "n")

    axis(1, at = 1:nperiods, labels = weeks)

    if (model == 2 || model == 3) {
        school = mydata$model$school
        school[school == 0] = NA
        par(new = TRUE)
        plot(school, ylim = c(0, 6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "grey", type = "p", pch = 22, lwd = 4)

    }

    if (model == 1 || model == 3) {
        sh = mydata$model$sh
        par(new = TRUE)
        plot(sh, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "black", type = "l", lwd = 1)

    }

    ##
    ## Repeat with the direct fit of the national mydata - if it was done!
    ##

    if (!is.null(model_rtn) && !is.null(model_profile)) {
        model_mean = rep(0, nperiods)
        for (iweek in 1:nperiods) model_mean[iweek] = mean(model_profile_ili[, iweek])
        ymax = max(model_rtn_ili[1:nperiodsData], model_profile_ili[, 1:nperiodsData], model_ili[1:nperiodsData], na.rm = TRUE)

        plot(1:nperiods, model_ili, type = "l", col = "black", xlim = c(1, nperiods), ylim = c(0, ymax), xaxt = "n", xlab = "EW #", ylab = mydata$model$raw_units,
            lwd = 2, lty = 1)
        points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

        for (irnd in irnd.set) {
            lines(1:nperiodsFit, model_profile_ili[irnd, 1:nperiodsFit], col = "lightgrey", lwd = 3, lty = 1)
            lines(nperiodsFit:nperiods, model_profile_ili[irnd, nperiodsFit:nperiods], col = "lightgrey", lwd = 3, lty = 2)
        }

        lines(1:nperiodsFit, model_rtn_ili[1:nperiodsFit], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        lines(nperiodsFit:nperiods, model_rtn_ili[nperiodsFit:nperiods], col = "red", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)

        lines(1:nperiodsFit, model_mean[1:nperiodsFit], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        lines(nperiodsFit:nperiods, model_mean[nperiodsFit:nperiods], col = "blue", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 2)
        lines(1:nperiods, model_ili, type = "l", col = "black", xaxt = "n", xlab = "", ylab = "", lwd = 2, lty = 1)
        points(1:nperiodsFit, model_ili[1:nperiodsFit], type = "p", col = "black", pch = 20, xaxt = "n", xlab = "", ylab = "")

        if (length(model_onset) > 0)
            abline(h = model_onset, col = "grey", lty = 2, lwd = 2)

        rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
        reg.name = paste(mydata$model$name, c("-Data", "-Direct-Best", "-Direct-Mean", "-Direct-Random"), sep = "")

        legend("topleft", c(mydata$FY, reg.name), text.col = c("black", "black", "red", "blue", "grey"), bty = "n")
        axis(1, at = 1:nperiods, labels = weeks)

        if (model == 2 || model == 3) {
            school = mydata$model$school
            school[school == 0] = NA
            par(new = TRUE)
            plot(school, ylim = c(0, 6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "grey", type = "p", pch = 22, lwd = 4)

        }

        if (model == 1 || model == 3) {
            sh = mydata$model$sh
            par(new = TRUE)
            plot(sh, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "black", type = "l", lwd = 1)

        }

    }

    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list, idevice = idevice, fit_ili = fit_ili, rtn_ili = rtn_ili, profile_ili = profile_ili,
        model_ili = model_ili, model_rtn_ili = model_rtn_ili, model_profile_ili = model_profile_ili, fit_model = fit_model, fit_model_mean = fit_model_mean, fit_model_profile = fit_model_profile)

    filename = pdfName = paste(subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)

    err = 0
    return(err)
}

plotHists <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1, run.list = NULL,
    idevice = 1) {

    #' Plot Histograms of Predicted and Observed Peak Value and Week
    #'
    #' Plots two histogram files one for the observed and predicted peak week and the other for the
    #' observed and predicted \% ILI value for all weeks included in the range of the number of weeks fitted
    #' to the number of weeks of mydata.  The \% ILI is presented in bins of 0.5\%.The default is to
    #' fit all the available mydata.  In each case we first plot the results for all the fit regions, followed by
    #' the results for an indirect (coupled or uncoupled) fit to the model region, and the last panel is for
    #' the direct fit to the model region.
    #' @param rtn A 1D numeric array with the best in-direct prediction to the model region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param mydata A list with the entire mydata set of this \code{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list a list with various parameters for the run
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotHists(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, run.list = run.list, idevice = 1)

    device = run.list$device[idevice]
    if (is.null(device))
        device = "pdf"
    if (is.null(mydata))
        return
    if (is.null(profile))
        return

    FY = mydata$FY
    model = mydata$imodel
    weeks = mydata$weeks
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$nperiodsData
    reg.fit.name = mydata$fit$name
    reg.model.name = mydata$model$name

    n.model = 1
    nregions = mydata$fit$nregions
    nRnd = dim(profile)[1]
    nRndiceData = nRnd + 1
    nregions1 = nregions + 1
    # the coefficients are given by the relative populations for safety normalize it

    fit_coef = mydata$fit$coef
    fit_onset = mydata$fit$onset
    model_coef = mydata$model$coef
    model_onset = mydata$model$onset

    tps = mydata$weeks

    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName
	myName = gsub(" ","",myName)

    factor = as.numeric(mydata$fit$factor)
    model_factor = mydata$model$factor
    model_ili = mydata$model$raw
    fit_ili = mydata$fit$raw

    rtn_ili = rtn
    profile_ili = profile
    for (i in 1:nregions) {
        rtn_ili[, i] = rtn[, i]/factor[i]
        profile_ili[, , i] = profile[, , i]/factor[i]
    }
    model_rtn_ili = NULL
    model_profile_ili = NULL
    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile/model_factor
    }
    ## for plotting - replace zero's with NA
    if (nperiodsFit < nperiods) {
        index <- which(fit_ili[, 1] == 0)

        if (length(index) >= 1) {
            fit_ili[index, 1:nregions] = NA
            model_ili[index] = NA
        }
    }

    colnames(fit_ili) = mydata$fit$attr$NAME_3
    colnames(rtn_ili) = mydata$fit$attr$NAME_3


    model_profile_indr = array(data = 0, dim = c(nRnd, nperiods))

    # this is the indirect modeling
    for (i in 1:nperiods) {
        for (j in 1:nRnd) {
            model_profile_indr[j, i] = sum(profile_ili[j, i, 1:nregions] * fit_coef[1:nregions])
        }
    }

    # direct fitting of nregions regions and indirect of the model
    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/hist-prfl-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n\n For a Histogram Plot of Predicted Profiles See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste0(subDir, "/hist-prfl-", myName, "-", nperiodsFit, "-", ireal, "_pg")
        cat("\n\n For a Histogram Plot of Predicted Profiles See: ", pngName, "\n\n")
        png(file = paste0(pngName, "%01d.png"), width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }
    nrow = 3
    ncol = length(nperiodsFit:min(nperiods, nperiodsFit + 4))

    par(mfrow = c(nrow, ncol), mar = c(4, 4, 1, 2))
    ## a histogram plot of binned mydata
    ylab = "Probability Density"
    xlab = "% ILI"

    for (iregion in 1:nregions) {

        for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
            min.val = 0
            if (iweek <= nperiodsData) {
                max.val = ceiling(max(fit_ili[iweek, iregion], profile_ili[, iweek, iregion], na.rm = TRUE))
                max.val = 2 * max.val
	      if (mydata$fit$factor[iregion] != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }

                hist(profile_ili[, iweek, iregion], breaks = breaks, xlim = c(min.val,max.val), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "",
                  ylab = "", main = "", freq = F)
                if (!is.na(fit_ili[iweek, iregion]))
                  hist(fit_ili[iweek, iregion], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab,
                    ylab = ylab, main = "", freq = F, add = TRUE)
                leg.text = c(mydata$FY)
                reg.name = reg.fit.name[iregion]
                rel.pop = round(fit_coef[iregion], digits = 3)
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                leg.text = c(leg.text, reg.name, rel.pop, my.week, "Data", "Direct")
                legend("topright", legend = leg.text, text.col = c("black", "black", "black", "black", rgb(0, 0, 1, 1/2), rgb(1, 0, 0,
                  1/2)), bty = "n")
            } else {
                max.val = ceiling(max(profile_ili[, iweek, iregion], na.rm = TRUE))
                max.val = 2 * max.val
	      if (mydata$fit$factor[iregion] != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
                hist(profile_ili[, iweek, iregion], breaks = breaks, xlim = c(min.val,max.val), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "",
                  ylab = "", main = "", freq = F)
                leg.text = c(mydata$FY)
                reg.name = reg.fit.name[iregion]
                rel.pop = round(fit_coef[iregion], digits = 3)
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                leg.text = c(leg.text, reg.name, rel.pop, my.week, "Direct")
                legend("topright", legend = leg.text, text.col = c("black", "black", "black", "black", rgb(1, 0, 0, 1/2)), bty = "n")

            }

        }
    }
    for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
        min.val = 0
        if (iweek <= nperiodsData) {
            max.val = ceiling(max(model_ili[iweek], model_profile_indr[, iweek], na.rm = TRUE))
            max.val = 2 * max.val
	      if (mydata$model$factor != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
            hist(model_ili[iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab,
                main = "", freq = F)
            hist(model_profile_indr[, iweek], breaks = breaks,xlim = c(min.val,max.val), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "",
                ylab = "", main = "", freq = F, add = TRUE)
            leg.text = c(mydata$FY)
            model.name = mydata$model$name
            my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
            gsub(model.name, ".", " ", model.name)

            leg.text = c(leg.text, model.name, my.week, "Data", "Indirect")
            legend("topright", legend = leg.text, text.col = c("black", "black", "black", rgb(0, 0, 1, 1/2), rgb(1, 0, 0, 1/2)), bty = "n")

        } else {

            max.val = ceiling(max(model_profile_indr[, iweek], na.rm = TRUE))
            max.val = 2 * max.val
	      if (mydata$model$factor != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
            hist(model_profile_indr[, iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "",
                ylab = "", main = "", freq = F)
            leg.text = c(mydata$FY)
            model.name = mydata$model$name
            my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
            gsub(model.name, ".", " ", model.name)

            leg.text = c(leg.text, model.name, my.week, "Indirect")
            legend("topright", legend = leg.text, text.col = c("black", "black", "black", rgb(1, 0, 0, 1/2)), bty = "n")

        }
    }
    if (!is.null(model_profile_ili)) {
        for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {

            if (iweek <= nperiodsData) {
                max.val = ceiling(max(model_ili[iweek], model_profile_ili[, iweek], na.rm = TRUE))
                max.val = 2 * max.val
	      if (mydata$model$factor != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
                 hist(model_ili[iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab,
                  main = "", freq = F)
                hist(model_profile_ili[, iweek], breaks = breaks, xlim = c(min.val,max.val), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = "",
                  ylab = "", main = "", freq = F, add = TRUE)

                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, my.week, "Data", "Direct")
                legend("topright", legend = leg.text, text.col = c("black", "black", "black", rgb(0, 0, 1, 1/2), rgb(0, 1, 0, 1/2)),
                  bty = "n")
            } else {
                max.val = ceiling(max(model_profile_ili[, iweek], na.rm = TRUE))
                max.val = 2 * max.val
	      if (mydata$model$factor != 1){
                	breaks = seq(from = min.val, to = max.val, by = 0.5)
                } else {
                	breaks = "Sturges"
                }
                hist(model_profile_ili[, iweek], breaks = breaks, xlim = c(min.val, max.val), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = "",
                  ylab = "", main = "", freq = F)
                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, my.week, "Direct")
                legend("topright", legend = leg.text, text.col = c("black", "black", "black", rgb(0, 1, 0, 1/2)), bty = "n")
            }
        }
    }


    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()


    colvec = rainbow(nregions)
    lwd = rep(2, nregions)

    # Plot the individual fits along with the model fit
    ipage <- 1
    iplots <- 0
    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/hist-week-max-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n\n For a Histogram Plot of Peak Week See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 16, height = 12)

    } else if (tolower(device) == "png") {
        pngName = paste(subDir, "/hist-week-max-", myName, "-", nperiodsFit, "-", ireal, "-pg", ipage, ".png", sep = "")
        cat("\n\n For a Histogram Plot of Peak Week See: ", pngName, "\n\n")
        png(file = pngName, width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }

    nrow = 3
    ncol = 4

    # maximum week for model and fit in the mydata
    dat_model_wk_max = which.max(model_ili)
    # if ((mydata$fit$level > mydata$model$level)) {
    if (nregions > 1) {
        dat_fit_wk_max = rep(0, nregions)
        for (i in 1:nregions) dat_fit_wk_max[i] = which.max(fit_ili[, i])
    }

    # maximum fit in direct modeling
    drct_model_wk_max = rep(0, nRnd)

    if (!is.null(model_profile_ili))
        for (i in 1:nRnd) drct_model_wk_max[i] = which.max(model_profile_ili[i, ])


    indrct_model_wk_max = rep(0, nRnd)
    for (i in 1:nRnd) indrct_model_wk_max[i] = which.max(model_profile_indr[i, ])

    sim_fit_wk_max = array(0, c(nRnd, nregions))

    for (i in 1:nregions) {
        for (j in 1:nRnd) {
            sim_fit_wk_max[j, i] = which.max(profile_ili[j, , i])
        }
    }

    par(mfrow = c(nrow, ncol), mar = c(4, 4, 1, 4))

    # Continue plotting if fitted at a higher resolution
    wk.min = round(min(dat_fit_wk_max, sim_fit_wk_max))
    wk.max = round(max(dat_fit_wk_max, sim_fit_wk_max))

    wk.min = round(0.5 * wk.min)
    wk.max = round(1.5 * wk.max)
    wk.max = min(wk.max, nperiods)
    wk.min = max(1, wk.min)

    breaks = seq(from = wk.min, to = wk.max, by = 1)
    ylab = "Probability Density"
    xlab = "EW #"

    for (iregion in 1:nregions) {
        hist(dat_fit_wk_max[iregion], breaks = breaks, xlim = range(breaks), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab,
            xaxt = "n", main = "", freq = FALSE)
        hist(sim_fit_wk_max[, iregion], breaks = breaks, xlim = range(breaks), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "", ylab = "",
            xaxt = "n", main = "", freq = FALSE, add = TRUE)
        axis(1, at = wk.min:wk.max, labels = weeks[wk.min:wk.max])
        leg.text = c(mydata$FY)
        reg.name = reg.fit.name[iregion]
        rel.pop = round(fit_coef[iregion], digits = 3)
        legend("topleft", c(leg.text, reg.name, rel.pop), bty = "n")
        leg.text = c("Peak Week", "Data", "Direct")
        legend("topright", legend = leg.text, text.col = c("black", rgb(0, 0, 1, 1/2), rgb(1, 0, 0, 1/2)), bty = "n")
        iplots <- iplots + 1  # count plots
        if ((device == "png") & (iplots%%12 == 0)) {
            dev.off()  # close plot
            ipage <- ipage + 1  # incurument page number
            pngName = paste(subDir, "/hist-week-max-", myName, "-", nperiodsFit, "-", ireal, "-pg", ipage, ".png", sep = "")
            cat("\n\n For a Histogram Plot of Peak Week See: ", pngName, "\n\n")
            png(file = pngName, width = 1200, height = 900)
            par(mfrow = c(nrow, ncol), mar = c(4, 4, 1, 4))
        }
    }


    # Plot the national
    wk.min = round(min(dat_model_wk_max, indrct_model_wk_max))
    wk.max = round(max(dat_model_wk_max, indrct_model_wk_max))
    if (!is.null(model_profile_ili)) {
        wk.min = round(min(wk.min, drct_model_wk_max))
        wk.max = round(max(wk.max, drct_model_wk_max))
    }
    wk.min = round(0.5 * wk.min)
    wk.max = round(1.5 * wk.max)
    wk.max = min(wk.max, nperiods)
    wk.min = max(1, wk.min)
    breaks = seq(from = wk.min, to = wk.max, by = 1)


    hist(dat_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab, xaxt = "n",
        main = "", freq = FALSE)
    hist(indrct_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(1, 0, 0, 1/4), border = NULL, xlab = "", ylab = "", xaxt = "n",
        main = "", freq = FALSE, add = TRUE)
    axis(1, at = wk.min:wk.max, labels = weeks[wk.min:wk.max])
    leg.text = c(mydata$FY)
    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    legend("topleft", c(leg.text, model.name), bty = "n")
    leg.text = c("Peak Week", "Data", "Indirect")
    legend("topright", legend = leg.text, text.col = c("black", rgb(0, 0, 1, 1/2), rgb(1, 0, 0, 1/2)), bty = "n")
    if (!is.null(model_profile_ili)) {
        hist(dat_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(0, 0, 1, 1/4), border = NULL, xlab = xlab, ylab = ylab,
            xaxt = "n", main = "", freq = FALSE)
        hist(drct_model_wk_max, breaks = breaks, xlim = range(breaks), col = rgb(0, 1, 0, 1/4), border = NULL, xlab = "", ylab = "",
            xaxt = "n", main = "", freq = FALSE, add = TRUE)
    }
    axis(1, at = wk.min:wk.max, labels = weeks[wk.min:wk.max])
    leg.text = c(mydata$FY)
    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    legend("topleft", c(leg.text, model.name), bty = "n")

    leg.text = c("Peak Week", "Data", "Direct")
    legend("topright", legend = leg.text, text.col = c("black", rgb(0, 0, 1, 1/2), rgb(0, 1, 0, 1/2)), bty = "n")
    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    return(err = 0)
}


## place holders for extrenal plotting routines that the User can implement

plotFitOnePatch.external <- function(model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = NULL, run.list = NULL, idevice = 1) {

    #' Plot the Results of a \pkg{DICE} Run - Using an External User Provided Routine
    #'
    #' Plot the results of \pkg{DICE} run for a single region/patch. We show the observed %ILI (or number of cases) along with our fits and
    #' if appropriate predictions. We show the best result and randomly selected results from the MCMC chain.
    #' For now this is just a wrapper that calls the \pkg{DICE} \code{\link{plotFitOnePatch.ggplot2}} function
    #' @param model_rtn A 1D numeric array with the best direct prediction to the region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the mydata
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotFitOnePatch{model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, device = device, idevice = 1}

    success = plotFitOnePatch.ggplot2(model_rtn = model_rtn, model_profile = model_profile, mydata = mydata, ireal = ireal,
        run.list = run.list, idevice = idevice)

    success

}


plotFitCDCPercentILI.external <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = NULL,
    run.list = NULL, idevice = 1) {
    #'
    #' Plot the Results of a \pkg{DICE} Run - Using an External User Provided Routine
    #'
    #' For now this is just a wrapper that calls the \pkg{DICE} \code{\link{plotFitCDCPercentILI.ggplot2}} function
    #' Plot the results of an a coupled or uncoupled  \pkg{DICE} run. For each of the fit regions we plot  the  \% ILI of the region along with
    #' our predictions for it based on randomly selected results from the  history of the MCMC chain of each region.  Using the predictions
    #' for the fit regions we then show the results for the model  region as a weighted sum of the fit regions.  The last panel
    #' shows our prediction for the model region using a direct fit to the model data.  The function also writes a binary RData file with
    #' all the profile predictions for the model and fit regions. Note that in the case of a coupled run the fit regions are never individually
    #' optimized. It is their weighted sum that is optimized, with the weights given by the relative population of each fit region.
    #'
    #' @param rtn A 1D numeric array with the best in-direct prediction to the model region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param mydata A dataframe with all the data available for this \code{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotFitCDCPercentILI{ rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, device = device, idevice = 1}


    success = plotFitCDCPercentILI.ggplot2(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
        ireal = ireal, run.list = run.list, idevice = idevice)
    success
}


plotHists.external <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = NULL, run.list = NULL,
    idevice = 1) {

    #' Plot Histograms of Predicted and Observed Peak Value and Week - Using an External User Provided Routine
    #'
    #' For now this is just a wrapper that calls the \pkg{DICE} \code{\link{plotHists.ggplot2}} function
    #' Plots two histogram files one for the observed and predicted peak week and the other for the
    #' observed and predicted \% ILI value for all weeks included in the range of the number of weeks fitted
    #' to the number of weeks of mydata.  The \% ILI is presented in bins of 0.5\%.The default is to
    #' fit all the available data.  In each case we first plot the results for all the fit regions, followed by
    #' the results for an indirect (coupled or uncoupled) fit to the model region, and the last panel is for
    #' the direct fit to the model region.
    #' @param rtn A 1D numeric array with the best in-direct prediction to the model region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list a list with various parameters for the run
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotHists.external(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, run.list = run.list, idevice = 1)

    success = plotHists.ggplot2(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, mydata = mydata,
        ireal = ireal, run.list = run.list, idevice = idevice)

    success
}

plotDisease <- function(mydata=NULL, all_years_epi = NULL, run.list = NULL, device = NULL) {
    #' Plot The entire Time Series of the Disease Data
    #'
    #' \code{plotDisease} Plots the disease time series for the country/level data along with the historic
    #' average
    #' @param mydata A dataframe with all the available data for this \pkg{DICE} run
    #' @param all_years_epi the epi data for all years
    #' @param run.list a list with various parameters for the run
    #' @param device - 'pdf or 'png'.  Default is 'png'
    #' plotDisease(mydata = mydata, all_years_epi = all_years_epi, run.list = run.list, device = 'pdf')
    #' @return err=0 if plots were created
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

   cadence.names = 1:nperiods
   if (mydata$cadence == "Monthly") {
   	cadence.names = month.abb[all_years_epi$months]
   	index.year.start = which(all_years_epi$months == 1)
   } else if (mydata$cadence == "Weekly") {
   	cadence.names = paste0("EW", all_years_epi$weeks)
   	index.year.start = which(all_years_epi$weeks == 1)
   } else if (mydata$cadence == "Daily") {
   	index.year.start = which(all_years_epi$days == 1)
    dates = all_years_epi$dates
 	dates = as.Date(dates, format='%Y-%m-%d')
 	cadence.names = format(dates, "%m-%d")
   } else {
	index.year.start = 1
   }
	if (is.null(index.year.start)) index.year.start = 1

   years = all_years_epi$years[index.year.start]

   nregions = mydata$fit$nregions
   nregions1 = nregions + 1

	myName = all_years_epi$mydataName
	if(is.null(myName)) myName = mydata$model$name
	myName = gsub(" ","",myName)

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
   if (nregions > 1) {
   	nrow = min(nregions1, 5)
   	ncol = 2
   	par(mfrow = c(nrow, ncol), mar = c(4, 5, 2, 1))
   }

   FY = all_years_epi$FY

   ylab = paste0(mydata$cadence, " # Cases")

   nts = length(all_years_epi$model$raw)
   lower = upper = cases.ave = rep(0, nts)


   for (j in 1:nperiods) {
   	if (mydata$cadence == "Monthly")
   		ind = which(all_years_epi$months == mydata$months[j])

   	if (mydata$cadence == "Weekly")
   		ind = which(all_years_epi$weeks == mydata$weeks[j])

   	if (mydata$cadence == "Daily")
   		ind = which(all_years_epi$days == mydata$days[j])

  	if (mydata$cadence == "nonuniform")
   		ind = which(all_years_epi$days == mydata$days[j])

   	cases.ave[ind] = mean(all_years_epi$model$raw[ind], na.rm = TRUE)
   	lower[ind] = quantile(all_years_epi$model$raw[ind], na.rm = TRUE, probs = c(0.05, 0.95))[1]
   	upper[ind] = quantile(all_years_epi$model$raw[ind], na.rm = TRUE, probs = c(0.05, 0.95))[2]
   }

   nts = length(lower)

   plot(all_years_epi$model$raw, type = "l", col = "red", xaxt = "n", xlab = paste0("Year ", FY), ylab = ylab, lwd = 2)
   polygon(c(1:nts, rev(1:nts)), c(upper, rev(lower)), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
   if (nts <= length(mydata$model$raw)) {
   	axis(1, at = 1:nts, labels = cadence.names) #cadence.names[1:nts] #all_years_epi$dates[1:nts]
   } else {
   	axis(1, at = index.year.start, labels = years)
   }
   abline(v = index.year.start, col = "blue", lwd = 1, lty = 2)
   lines(cases.ave, type = "l", col = "black")
   legend("topleft", mydata$model$name, bty = "n")
   legend("topright", legend = c("Data", "Historic Average", "95% CI"), text.col = c("red", "black", rgb(0, 0, 0.6, 0.2)), bty = "n")


   if (nregions > 1) {

   	for (i in 1:nregions) {
   		nts = length(all_years_epi$fit$raw[, i])
   		lower = upper = cases.ave = rep(0, nts)
   		for (j in 1:nperiods) {
   			if (mydata$cadence == "Monthly")
   				ind = which(all_years_epi$months == mydata$months[j])
   			if (mydata$cadence == "Weekly")
   				ind = which(all_years_epi$weeks == mydata$weeks[j])
   			cases.ave[ind] = mean(all_years_epi$fit$raw[ind, i], na.rm = TRUE)
   			lower[ind] = quantile(all_years_epi$fit$raw[ind, i], na.rm = TRUE, probs = c(0.05, 0.95))[1]
   			upper[ind] = quantile(all_years_epi$fit$raw[ind, i], na.rm = TRUE, probs = c(0.05, 0.95))[2]
   		}

   		plot(all_years_epi$fit$raw[, i], type = "l", col = "red", xaxt = "n", xlab = "Year", ylab = ylab)
   		polygon(c(1:nts, rev(1:nts)), c(upper, rev(lower)), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
   		if (nts <= length(mydata$fit$raw[,i])) {
   			axis(1, at = 1:nts, labels = cadence.names)
   		} else {
   			axis(1, at = index.year.start, labels = years)
   		}
    	abline(v = index.year.start, col = "blue", lwd = 1, lty = 2)
   		lines(cases.ave, type = "l", col = "black")
   		legend("topleft", mydata$fit$name[i], bty = "n")
   		legend("topright", legend = c("Data", "Historic Average", "95% CI"), text.col = c("red", "black", rgb(0, 0, 0.6, 0.2)), bty = "n")

   		if (i%%nrow == 0) {
   			if ((tolower(device) != 'pdf') & (tolower(device) != 'png'))
   				dev.new()
   		}
   	}
   }

	if (tolower(device) == 'pdf' || tolower(device) == 'png') dev.off()

	err = 0
	return(err)
}

plotMECH <- function(mydata = NULL, tables.mod = NULL, tables.fit = NULL, tables.agg.mod = NULL,  mod_id = NULL, fit_id = NULL, ymax.input = NULL, ireal = NULL, run.list = NULL, idevice = 1) {

    #' Plot Mechanistic Results of Fitting/Forecasting the Disease Data
    #'
    #' \code{plotMECH} Creates a PDF file with the results of the MCMC fits to the model and
    #' the fit level data
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param tables.mod - a table with the data and the results for the model fit. It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.fit - Optional a table with the data and the results for the fits at the fit_level.
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.agg.mod - Optional. A table with the aggregate of the uncoupled fit_level results present only in the case of an uncoupled run
    #' @param mod_id The abbreviation of the states/regions-model level
    #' @param fit_id The abbreviation of the states/regions-fit level
    #' @param ymax.input Optional, Maximum value for y-axis in plots (numeric)
    #' @param ireal - Numeric, the number of the MCMC chain
    #' @examples
    #' plotMECH(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL,
    #' mod_id = mod_id, fit_id = fit_id, ymax.input = NULL, ireal = ireal)
    #' @return err=0 if plots were created
    #'
    #'

    device = run.list$device[idevice]
    if (is.null(device))
        device = "png"

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
    	epi_model = 'SIR-B'
    } else {
    	epi_model = 'UNKNOWN'
    }

    ## subset only the state we are plotting
    observed = tables.mod$epi.obsrv[, mod_id]
    epi.null = tables.mod$epi.null[, , mod_id]
    #arima = tables$arima.frcst[, , state_id]
    #arima.lower = tables$arima.lower[, , state_id]
    #arima.upper = tables$arima.upper[, , state_id]

    seir = tables.mod$seir.frcst[, , mod_id]
    seir.lower = tables.mod$seir.lower[, , mod_id]
    seir.upper = tables.mod$seir.upper[, , mod_id]

    # Data including augmentation - may be the same as original if no augmentation used
    if (mydata$da > 0) {
    	epi.mod.agmnt  = as.numeric(mydata$model$epirun)
    } else {
    	epi.mod.agmnt = mydata$model$raw
    }

    cadence = mydata$cadence

    if (cadence == "Monthly") {
        month.names = month.abb[mydata$months]
        dates = month.names  #  paste0(month.name,'-',year.vec)
        ind = seq(from = 1, to = length(dates))
    } else if (cadence == "Weekly") {
        dates = mydata$weeks
        dates = sub("week", "EW", dates)
        ind = seq(1, length(dates), by = 4)
        dates = dates[ind]
    } else if (cadence == "Daily") {
        dates = mydata$dates
        ind = seq(1, length(dates), by = 7)
        dates = dates[ind]
 		dates = as.Date(dates, format='%Y-%m-%d')
 		dates = format(dates, "%b-%d")
        dates = dates[ind]
    } else if (cadence == "nonuniform") {
        dates = mydata$dates
        ind = seq(1, length(dates), by = 1)
        dates = dates[ind]
 		dates = as.Date(dates, format='%Y-%m-%d')
 		dates = format(dates, "%b-%d")
    } else {
        dates = mydata$dates
        ind = seq(1, length(dates), by = 7)
        dates = dates[ind]
 		dates = as.Date(dates, format='%Y-%m-%d')
 		dates = format(dates, "%b-%d")
    }

    nregions = mydata$fit$nregions
    nregions1 = nregions + 1
    nregions2 = nregions + 2

    my.dims = dim(seir)

    nmydata = mydata$nperiods

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

	ymin = 0.0

	if (tolower(mydata$disease) == 'dengue' || tolower(mydata$disease) == 'yellow_fever' || tolower(mydata$disease) == 'plague' || tolower(mydata$disease) == 'chik' || tolower(mydata$disease) == 'sars') {
		ylab = 'Cases'
	} else if (tolower(mydata$disease) == 'flu') {
		ylab = ' % ILI'
	} else {
		ylab = 'Incidence'
	}

  		if (!is.null(tables.fit)) {

 			## subset only the states we are plotting
 			observed = tables.fit$epi.obsrv[, fit_id]
 			epi.null = tables.fit$epi.null[, , fit_id]


 			seir = tables.fit$seir.frcst[, , fit_id]
 			seir.lower = tables.fit$seir.lower[, , fit_id]
 			seir.upper = tables.fit$seir.upper[, , fit_id]

    		# Data including augmentation - may be the same as original if no augmentation used
    		if (mydata$da > 0){
    			epi.fit.agmnt  = as.matrix(mydata$fit$epirun)
    		}	else {
    			epi.fit.agmnt  = as.matrix(mydata$fit$raw)
    		}

    		colnames(epi.fit.agmnt) = fit_id

 			nstates = length(fit_id)

 			colvec_save = rainbow(nstates)
 			colvec = col2rgb(colvec_save)

 			for (istate in 1:nstates) {

 				id = fit_id[istate]

 				title = paste0(mydata$fit$name[istate], " - ", mydata$FY, " Season")

 				xlab = ""

 				ymax = ymax.input

 				if (length(ymax.input) == 0) {

 					ymax = max(observed[,id], seir.upper[nperiodsFit, , id], na.rm = TRUE)
 				}


 				xlab = paste0("Time (", cadence, ")")

 				plot(1:nmydata, observed[, id], type = "n", col = "red", lwd = 1, xlab = xlab, ylab = ylab, ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata),
 					xaxt = "n")

				poly_col = as.numeric(colvec[,istate])/255

 				polygon(c(1:nmydata, rev(1:nmydata)), c(seir.upper[nperiodsFit, , id], rev(seir.lower[nperiodsFit, , id])), col = rgb(poly_col[1], poly_col[2], poly_col[3], 0.2), border = FALSE)
 				lines(1:nmydata, seir[nperiodsFit, , id], type = "l", col = colvec_save[istate], xlab = "", ylab = "")
 				lines(1:nmydata, observed[, id], type = "p", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "")

 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "l", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "")
 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "p", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "", pch = 19)
 				if(mydata$da >= 0) lines(1:nperiods, epi.fit.agmnt[1:nperiods, id], type = 'l', col = 'black', lty = 2, lwd = 1, xlab = "", ylab = "")
 				lines(1:nmydata, epi.null[nperiodsFit, , id], type = "l", col = "darkgrey", lwd = 2, lty = 1)
 				rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

 				if (all(is.na(epi.null[nperiodsFit, , id]))) {
 					legend("topleft", c("Observed", epi_model), text.col = c("black", colvec_save[istate]), bty = "n")
 				} else {
 					legend("topleft", c("Observed", "NULL", epi_model), text.col = c("black", "darkgrey",colvec_save[istate]), bty = "n")
 				}

 				axis(1, at = 1:nmydata, labels = FALSE)
 				axis(1, at = ind, labels = dates, las = 2)


 			}
 		}

		if (!is.null(tables.agg.mod)) {
			## subset only the state we are plotting
			observed = tables.agg.mod$epi.obsrv[, mod_id]
			epi.null = tables.agg.mod$epi.null[, , mod_id]


			seir = tables.agg.mod$seir.frcst[, , mod_id]
			seir.lower = tables.agg.mod$seir.lower[, , mod_id]
			seir.upper = tables.agg.mod$seir.upper[, , mod_id]

			if (length(ymax.input) == 0) {

				ymax = max(observed, seir.upper[nperiodsFit, ], na.rm = TRUE)
			}

			title = paste0(mydata$model$name, " - ", mydata$FY, " season")

			xlab = paste0("Time (", cadence, ")")

			plot(1:nmydata, observed, type = "n", col = "black", lwd = 1, xlab = xlab, ylab = ylab, ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

			polygon(c(1:nmydata, rev(1:nmydata)), c(seir.upper[nperiodsFit, ], rev(seir.lower[nperiodsFit, ])), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
			lines(1:nmydata, seir[nperiodsFit, ], type = "l", col = "red", lwd = 2, xlab = "", ylab = "")
			lines(1:nmydata, observed, type = "p", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "")
			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "")
			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "black", lwd = 2, lty = 1, xlab = "", ylab = "", pch = 19)
			if(mydata$da >= 0) lines(1:nperiods, epi.mod.agmnt[1:nperiods], type = 'l', col = 'red', lty = 2, lwd = 1, xlab = "", ylab = "")
			lines(1:nmydata, epi.null[nperiodsFit, ], type = "l", col = "darkgrey", lwd = 1, lty = 1)
			rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

			if (all(is.na(epi.null[nperiodsFit, ]))) {
				if (isingle == 1) { #Uncoupled
					legend("topleft", c("Observed", paste0(epi_model, " aggregate")), text.col = c("black", "red"), bty = "n")
				} else { #Coupled
					legend("topleft", c("Observed", paste0(epi_model, " coupled")), text.col = c("black", "red"), bty = "n")
				}
			} else {
				if (isingle == 1) { #Uncoupled
					legend("topleft", c("Observed", "NULL", paste0(epi_model, " aggregate")), text.col = c("black", "darkgrey",
						"red"), bty = "n")
				} else { #Coupled
					legend("topleft", c("Observed", "NULL", paste0(epi_model, " coupled")), text.col = c("black", "darkgrey", "red"),
						bty = "n")
				}
			}

			legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
			axis(1, at = 1:nmydata, labels = FALSE)
			axis(1, at = ind, labels = dates, las = 2)


		}

		## Plot direct fit to National

		## subset only the state we are plotting
		observed = tables.mod$epi.obsrv[, mod_id]
		epi.null = tables.mod$epi.null[, , mod_id]


		seir = tables.mod$seir.frcst[, , mod_id]
		seir.lower = tables.mod$seir.lower[, , mod_id]
		seir.upper = tables.mod$seir.upper[, , mod_id]

    	title = paste0(mydata$model$name, " - ", mydata$FY, " season")

    	xlab = ""

    	ymax = ymax.input

    	if (length(ymax.input) == 0) {

    		ymax = max(observed, seir.upper[nperiodsFit, ], na.rm = TRUE)
    	}


    	xlab = paste0("Time (", cadence, ")")

    	plot(1:nmydata, observed, type = "n", col = "black", lwd = 1, xlab = xlab, ylab = ylab, ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

    	polygon(c(1:nmydata, rev(1:nmydata)), c(seir.upper[nperiodsFit, ], rev(seir.lower[nperiodsFit, ])), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
    	lines(1:nmydata, seir[nperiodsFit, ], type = "l", col = "red", xlab = "", ylab = "")
    	lines(1:nmydata, observed, type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
    	lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")
    	lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)
    	if(mydata$da >= 0) lines(1:nperiods, epi.mod.agmnt[1:nperiods], type = 'l', col = 'black', lty = 2, lwd = 1, xlab = "", ylab = "")
    	lines(1:nmydata, epi.null[nperiodsFit, ], type = "l", col = "darkgrey", lwd = 1, lty = 1)
    	rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
    	if (all(is.na(epi.null[nperiodsFit, ]))) {
   		legend("topleft", c("Observed", paste0(epi_model, " direct")), text.col = c("black", "red"), bty = "n")
   	} else {
   		legend("topleft", c("Observed", "NULL", paste0(epi_model, " direct")), text.col = c("black", "darkgrey", "red"),
   			bty = "n")
   	}
    	legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
    	axis(1, at = 1:nmydata, labels = FALSE)
    	axis(1, at = ind, labels = dates, las = 2)

    dev.off()

    return(err = 0)
}


plotARIMA <- function(mydata = NULL, tables.mod = NULL, tables.fit = NULL, tables.agg.mod = NULL, arima_model = NULL, arima_model_all = NULL, covar = NULL, covar_lag = 0, device = 'png') {

    #' Plot SARIMA Results of Fitting/Forecasting the Disease Data
    #'
    #' \code{plotARIMA} Creates a PDF file with the results of the SARIMA fits
    #' to the model and fit data.  Both direct fitting and the aggregate results are plotted.
    #' @param all_years_epi the entire time series for the disease
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param tables.mod - a table with the data and the results for the direct Model fit.
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.fit - Optional a table with the data and the results for the fits at the fit_level.
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.agg.mod - A table with aggregate results of the fits
    #' @param arima_model - A list with the user selection of the ARIMA model.
    #' If the user has not selected a model this will be NULL and arima_model_all will have the models selected by auto.arima
    #' @param arima_model_all - An array with details of the ARIMA model used for
    #' model region, the fits regions and the aggregate (the last is only if the User
    #' has chosen an ARIMA model for the run)
    #' @param covar - String. Optional covariate variable for ARIMA fit.
    #' 'sh', 'precip' and 'temp' are currently supported. Default is NULL - no covariate
    #' @param covar_lag - Numeric. Lag time for covariate variable in units of the cadence
    #' @param device - 'pdf' or 'png'. Default is 'png'
    #' of the data
    #' @examples
    #' plotARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit,
    #' tables.agg.mod = tables.agg.mod,arima_model = arima_model,
    #' arima_model_all = NULL, covar = covar, covar_lag = covar_lag, device = device)
    #' @return err=0 if plots were created
    #'
    #'

    if (is.null(device))
        device = "png"

    disease = toupper(mydata$disease)

    subDir = mydata$subDir

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    nregions = mydata$fit$nregions
    nregions1 = nregions + 1
    nregions2 = nregions + 2

	covar = mydata$covar

    if (!is.null(arima_model)) {
    	p = rep(as.numeric(arima_model[["p"]]), nregions2)
    	d = rep(as.numeric(arima_model[["d"]]), nregions2)
    	q = rep(as.numeric(arima_model[["q"]]), nregions2)

    	P = rep(as.numeric(arima_model[["P"]]), nregions2)
    	D = rep(as.numeric(arima_model[["D"]]), nregions2)
    	Q = rep(as.numeric(arima_model[["Q"]]), nregions2)

    	epi_model = rep(0, nregions2)
    	for (iregion in 1:nregions2) epi_model[iregion] = paste0("(", p[iregion], ",", d[iregion], ",", q[iregion], ")(", P[iregion], ",", D[iregion],
    		",", Q[iregion], ")", nperiods)

    	epi.model = paste0(p[1], d[1], q[1], "-", P[1], D[1], Q[1])

    	if (covar != FALSE) {
    		epi.model = paste0(epi.model,'-',covar,'-lag-',covar_lag)
    	}
    } else {
    	p = d = q = P = D = Q = epi_model = rep(0, nregions2)
    	istart = 1
    	if (is.null(tables.fit)) istart = 2
    	for (iregion in istart:nregions1) {
    		tmp.arima <- arima_model_all[[iregion]]

    		p[iregion] = as.numeric(tmp.arima[["p"]])
    		d[iregion] = as.numeric(tmp.arima[["d"]])
    		q[iregion] = as.numeric(tmp.arima[["q"]])

    		P[iregion] = as.numeric(tmp.arima[["P"]])
    		D[iregion] = as.numeric(tmp.arima[["D"]])
    		Q[iregion] = as.numeric(tmp.arima[["Q"]])
    		epi_model[iregion] = paste0("(", p[iregion], ",", d[iregion], ",", q[iregion], ")(", P[iregion], ",", D[iregion], ",", Q[iregion], ")",
    			nperiods)
    	}
    	epi_model[nregions2] = "aggregate-arima"
    	epi.model = "auto-arima"
    }


    ##
    ## If a covar was used - it will be plotted
    ##

    if(covar != FALSE) {
    	if (tolower(covar) == "sh") {
			model_covar = mydata$model$sh
			fit_covar   = mydata$fit$sh
		} else if (tolower(covar) == "precip") {
			model_covar = mydata$model$precip
			fit_covar   = mydata$fit$precip
		} else if (tolower(covar) == "temp") {
			model_covar = mydata$model$temp
			fit_covar   = mydata$fit$temp
		} else if (tolower(covar) == "school") {
			model_covar = mydata$model$school
			fit_covar   = mydata$fit$school
		} else {
			model_covar = mydata$model$sh
			fit_covar   = mydata$fit$sh
		}
    }

    cadence = mydata$cadence

    if (cadence == "Monthly") {
        month.names = month.abb[mydata$months]
        dates = month.names  #  paste0(month.name,'-',year.vec)
        ind = seq(from = 1, to = length(dates))
    } else if (cadence == 'Weekly'){
        dates = mydata$weeks
        dates = sub("week", "EW", dates)
        ind = seq(1, length(dates), by = 4)
        dates = dates[ind]
    } else if (cadence == "Daily") {
        dates = mydata$dates
        ind = seq(1, length(dates), by = 7)
        dates = dates[ind]
 		dates = as.Date(dates, format='%Y-%m-%d')
 		dates = format(dates, "%m-%d")
    } else {
    	nperiods = mydata$nperiods
    	dates = 1:nperiods
    	ind = seq(1, length(dates), by = 1)
    }

    my.dims = dim(arima)

    nmydata = mydata$nperiods

    nperiodsFit = mydata$nperiodsFit
    myName = mydata$dataName
	myName = gsub(" ","",myName)

    if (tolower(device) == "pdf") {
    	filename = paste0(subDir, "/results-sarima-", myName, "-", epi.model,'-',nperiodsFit, ".pdf", sep = "")
        pdf(file = filename, onefile = TRUE, width = 18, height = 11)
    } else if (tolower(device) == "png") {
        filename = paste0(subDir, "/results-sarima-", myName, "-", epi.model,'-',nperiodsFit, ".png", sep = "")
        png(file = filename, width = 1800, height = 1100)
    } else {
        dev.next()
        dev.new()
    }


    cat("\nFor a plot of SARIMA results see: ", filename, "\n")

    if (nregions > 1) {
    	ncol = 3
    	nrow = round(nregions2 / ncol)
    	if ((nrow * ncol) < nregions2) nrow = nrow + 1
    	if (nrow > 3) nrow = min(nrow,3)
		if (nregions == 10) {
			nrow = 3
			ncol = 4
		}
    	par(mar = c(4,4,1,1), mfrow = c(nrow, ncol))
    } else {
    	nrow = 1
    	ncol = 1
    	par(mar = c(4,4,1,1), mfrow = c(nrow, ncol))
    }

	ymin = 0.0

  		if (!is.null(tables.fit)) {

 			## subset only the states we are plotting
 			observed = tables.fit$epi.obsrv[, fit_id]
 			epi.null = tables.fit$epi.null[, , fit_id]

 			arima = tables.fit$arima.frcst[, , fit_id]
 			arima.lower = tables.fit$arima.lower[, , fit_id]
 			arima.upper = tables.fit$arima.upper[, , fit_id]

 			nstates = length(fit_id)
 			for (istate in 1:nstates) {

 				#observed[observed[,istate] == 0, istate]  <- NA
 				epi.null[nperiodsFit,   !is.finite(epi.null[nperiodsFit,,istate]), istate]  <- NA
 				arima[  nperiodsFit,    !is.finite(arima[nperiodsFit,,istate])   , istate]  <- NA
 				arima.lower[nperiodsFit,!is.finite(arima.lower[nperiodsFit,,istate]), istate]  <- NA
 				arima.upper[nperiodsFit,!is.finite(arima.upper[nperiodsFit,,istate]), istate]  <- NA


 				id = fit_id[istate]

 				title = paste0(mydata$fit$name[istate], ": ", mydata$FY, " season")

 				xlab = ""

 				ymax = max(observed[,id], arima.upper[nperiodsFit, , id], na.rm = TRUE)

 				xlab = paste0("Time (", cadence, ")")

 				plot(1:nmydata, observed[, id], type = "n", col = "red", lwd = 1, xlab = xlab, ylab = "Cases", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata),
 					xaxt = "n")

 				polygon(c(1:nmydata, rev(1:nmydata)), c(arima.upper[nperiodsFit, , id], rev(arima.lower[nperiodsFit, , id])), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
 				lines(1:nmydata, arima[nperiodsFit, , id], type = "l", col = "blue", xlab = "", ylab = "")
 				lines(1:nmydata, observed[, id], type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")

 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "l", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)
 				lines(1:nmydata, epi.null[nperiodsFit, , id], type = "l", col = "black", lwd = 1, lty = 1)
 				rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
 				legend("topleft", c("Observed", "NULL", epi_model[istate]), text.col = c("red", "black", "blue"), bty = "n")
 				legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
 				axis(1, at = 1:nmydata, labels = FALSE)
 				axis(1, at = ind, labels = dates, las = 2)
 				if(covar != FALSE) {
 					par(new = TRUE)
 					plot(1:nmydata, fit_covar[,istate], type='l', lwd = 1, col = 'cyan', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
 					legend('topleft', c("","","",covar),text.col='cyan', bty='n')
				}
 			}
 		}

		if (!is.null(tables.agg.mod)) {

			## subset only the state we are plotting
			observed = tables.agg.mod$epi.obsrv[, mod_id]
			epi.null = tables.agg.mod$epi.null[, , mod_id]


			arima = tables.agg.mod$arima.frcst[, , mod_id]
			arima.lower = tables.agg.mod$arima.lower[, , mod_id]
			arima.upper = tables.agg.mod$arima.upper[, , mod_id]

			#observed[observed == 0] <- NA
			epi.null[nperiodsFit,!is.finite(epi.null[nperiodsFit,])] <- NA
			arima[      nperiodsFit,   !is.finite(arima[nperiodsFit,])   ] <- NA
			arima.lower[nperiodsFit,!is.finite(arima.lower[nperiodsFit,])] <- NA
			arima.upper[nperiodsFit,!is.finite(arima.upper[nperiodsFit,])] <- NA

			ymax = max(observed, arima.upper[nperiodsFit, ], na.rm = TRUE)

			title = paste0(mydata$model$name, ": ", mydata$FY, " season")

			xlab = paste0("Time (", cadence, ")")

			plot(1:nmydata, observed, type = "n", col = "red", lwd = 1, xlab = xlab, ylab = "Cases", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

			polygon(c(1:nmydata, rev(1:nmydata)), c(arima.upper[nperiodsFit, ], rev(arima.lower[nperiodsFit, ])), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
			lines(1:nmydata, arima[nperiodsFit, ], type = "l", col = "blue", xlab = "", ylab = "")
			lines(1:nmydata, observed, type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)
			lines(1:nmydata, epi.null[nperiodsFit, ], type = "l", col = "black", lwd = 1, lty = 1)
			rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

			legend("topleft", c("Observed", "NULL", paste0(epi_model[nregions2], " aggregate")), text.col = c("red", "black", "blue"), bty = "n")
			legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
			axis(1, at = 1:nmydata, labels = FALSE)
			axis(1, at = ind, labels = dates, las = 2)
			if(covar != FALSE) {
				par(new = TRUE)
				plot(1:nmydata, model_covar, type='l', lwd = 1, col = 'cyan', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
				legend('topleft', c("","","",covar),text.col='cyan', bty='n')
			}


		}

		## Plot direct fit to National

		## subset only the state we are plotting
		observed = tables.mod$epi.obsrv[, mod_id]
		epi.null = tables.mod$epi.null[, , mod_id]


		arima = tables.mod$arima.frcst[, , mod_id]
		arima.lower = tables.mod$arima.lower[, , mod_id]
		arima.upper = tables.mod$arima.upper[, , mod_id]


		#observed[              observed == 0] <- NA
		epi.null[nperiodsFit   ,!is.finite(epi.null[nperiodsFit,])] <- NA
		arima[   nperiodsFit   ,!is.finite(arima[nperiodsFit,])   ] <- NA
		arima.lower[nperiodsFit,!is.finite(arima.lower[nperiodsFit,])] <- NA
		arima.upper[nperiodsFit,!is.finite(arima.upper[nperiodsFit,])] <- NA

    	title = paste0(mydata$model$name, " - ", mydata$FY, " season")

    	xlab = ""

    	ymax = max(observed, arima.upper[nperiodsFit, ], na.rm = TRUE)

    	xlab = paste0("Time (", cadence, ")")

    	plot(1:nmydata, observed, type = "n", col = "red", lwd = 1, xlab = xlab, ylab = "Cases", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

    	polygon(c(1:nmydata, rev(1:nmydata)), c(arima.lower[nperiodsFit, ], rev(arima.upper[nperiodsFit, ])), col = rgb(0, 0, 0.6, 0.2), border = FALSE)
    	lines(1:nmydata, arima[nperiodsFit, ], type = "l", col = "blue", xlab = "", ylab = "")
    	lines(1:nmydata, observed, type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
    	lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "")
    	lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "red", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)
    	lines(1:nmydata, epi.null[nperiodsFit, ], type = "l", col = "black", lwd = 1, lty = 1)
    	rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
    	legend("topleft", c("Observed", "NULL", paste0(epi_model[nregions1], " direct")), text.col = c("red", "black", "blue"), bty = "n")
    	legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
    	axis(1, at = 1:nmydata, labels = FALSE)
    	axis(1, at = ind, labels = dates, las = 2)
		if(covar != FALSE) {
			par(new = TRUE)
			plot(1:nmydata, model_covar, type='l', lwd = 1, col = 'cyan', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
			legend('topleft', c("","","",covar),text.col='cyan', bty='n')
		}
    dev.off()

    ##
    ## Repeate for the Error of NULL and SARIMA models
    ##

    if (tolower(device) == "pdf") {
    	filename = paste0(subDir, "/error-sarima-", myName, "-", epi.model,'-',nperiodsFit, ".pdf", sep = "")
        pdf(file = filename, onefile = TRUE, width = 18, height = 11)
    } else if (tolower(device) == "png") {
        filename = paste0(subDir, "/error-sarima-", myName, "-", epi.model,'-',nperiodsFit, ".png", sep = "")
        png(file = filename, width = 1800, height = 1100)
    } else {
        dev.next()
        dev.new()
    }

    cat("\nFor a plot of SARIMA and NULL Errors See: ", filename, "\n")

    par(mar = c(4, 4, 1, 1), mfrow = c(nrow, ncol))

	nperiodsData = mydata$nperiodsData

  		if (!is.null(tables.fit)) {

 			## subset only the states we are plotting

			observed = tables.fit$epi.obsrv[, fit_id]
			epi.null = tables.fit$epi.nul[  , fit_id]
 			arima = tables.fit$arima.frcst[, , fit_id]

 			null.mae = tables.fit$epi.null.mae[, , fit_id]
			null.mre = tables.fit$epi.null.mre[, , fit_id]

 			arima.mae = tables.fit$arima.mae[, , fit_id]
 			arima.mre = tables.fit$arima.mre[, , fit_id]
			arima.rel = tables.fit$arima.rel[, , fit_id]

 			for (istate in 1:nstates) {

				observed[observed[, istate] == 0, istate] <- NA
				if (nperiodsData < nperiods) {
					null.mae[nperiodsFit, (nperiodsData + 1):nperiods, istate] <- NA
					null.mre[nperiodsFit, (nperiodsData + 1):nperiods, istate] <- NA

					arima.mae[nperiodsFit, (nperiodsData + 1):nperiods, istate] <- NA
					arima.mre[nperiodsFit, (nperiodsData + 1):nperiods, istate] <- NA
					arima.rel[nperiodsFit, (nperiodsData + 1):nperiods, istate] <- NA
				}
				index <- is.nan(null.mre[nperiodsFit, , istate]	)
				null.mre[nperiodsFit, index , istate] <- NA
				index <- is.infinite(null.mre[nperiodsFit, , istate])
				null.mre[nperiodsFit, index , istate] <- NA

				index <- is.nan(arima.mre[nperiodsFit, , istate]	)
				arima.mre[nperiodsFit, index , istate] <- NA
				index <- is.infinite(arima.mre[nperiodsFit, , istate])
				arima.mre[nperiodsFit, index , istate] <- NA

 				id = fit_id[istate]

 				title = paste0(mydata$fit$name[istate], ": ", mydata$FY, " season")

 				xlab = ""

    			ymax = max(arima.mre[nperiodsFit, ,id],null.mre[nperiodsFit, ,id], na.rm=TRUE)
				ymin = min(arima.mre[nperiodsFit, ,id],null.mre[nperiodsFit, ,id], na.rm=TRUE)

 				xlab = paste0("Time (", cadence, ")")

 				plot(1:nmydata, null.mre[nperiodsFit, , istate], type = "n", col = "red", lwd = 2, xlab = xlab, ylab = "Mean Relative Error", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata),
 					xaxt = "n")
 				lines(1:nmydata, null.mre[nperiodsFit, ,istate], type='l', col='red', xlab='',ylab='', lty = 2, lwd =2)
 				lines(1:nperiodsFit, null.mre[nperiodsFit,1:nperiodsFit , istate], type='l', col='red', xlab='',ylab='',lty = 1, lwd =2)

 				lines(1:nmydata, arima.mre[nperiodsFit, , istate], type = "l", col = "blue", xlab = "", ylab = "", lty = 2, lwd=2)
 				lines(1:nperiodsFit, arima.mre[nperiodsFit,1:nperiodsFit , istate], type = "l", col = "blue", xlab = "", ylab = "", lty = 1, lwd=2)

 				rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)
 				legend("topleft", c("NULL", epi_model[iregion], 'Incidence'), text.col = c("red", "blue", 'black'), bty = "n")
 				legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
 				axis(1, at = 1:nmydata, labels = FALSE)
 				axis(1, at = ind, labels = dates, las = 2)
   				par(new=TRUE)
    			plot(1:nmydata, observed[,id], type = "n", col = "black", lwd = 2, xlab = '', ylab = '', main = '', xlim = c(1, nmydata), xaxt = "n", yaxt='n')
				lines(1:nmydata, observed[,id], type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")

 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "l", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")
 				lines(1:nperiodsFit, observed[1:nperiodsFit, id], type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)

 			}
 		}

		if (!is.null(tables.agg.mod)) {
 			## subset only the states we are plotting

			observed = tables.agg.mod$epi.obsrv[, mod_id]
			epi.null = tables.agg.mod$epi.nul[  , mod_id]
 			arima = tables.agg.mod$arima.frcst[, , mod_id]

 			null.mae = tables.agg.mod$epi.null.mae[, , mod_id]
			null.mre = tables.agg.mod$epi.null.mre[, , mod_id]

 			arima.mae = tables.agg.mod$arima.mae[, , mod_id]
 			arima.mre = tables.agg.mod$arima.mre[, , mod_id]
			arima.rel = tables.agg.mod$arima.rel[, , mod_id]

			observed[observed == 0] <- NA
			if (nperiodsData < nperiods) {
				null.mae[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
				null.mre[nperiodsFit, (nperiodsData + 1):nperiods] <- NA

				arima.mae[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
				arima.mre[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
				arima.rel[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
			}
			index <- is.nan(null.mre[nperiodsFit, ]	)
			null.mre[nperiodsFit, index ] <- NA
			index <- is.infinite(null.mre[nperiodsFit, ])
			null.mre[nperiodsFit, index ] <- NA

			index <- is.nan(arima.mre[nperiodsFit, ]	)
			arima.mre[nperiodsFit, index ] <- NA
			index <- is.infinite(arima.mre[nperiodsFit, ])
			arima.mre[nperiodsFit, index ] <- NA


    		ymax = max(arima.mre[nperiodsFit, ],null.mre[nperiodsFit, ], na.rm=TRUE)
			ymin = min(arima.mre[nperiodsFit, ],null.mre[nperiodsFit, ], na.rm=TRUE)


			title = paste0(mydata$model$name, ": ", mydata$FY, " season")

			xlab = paste0("Time (", cadence, ")")

			plot(1:nmydata, null.mre[nperiodsFit, ], type = "n", col = "red", xlab = xlab, ylab = "Mean Relative Error", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

			lines(1:nmydata, null.mre[nperiodsFit, ], type = "l", col = "red", lwd = 2, lty = 2)
			lines(1:nperiodsFit, null.mre[nperiodsFit, 1:nperiodsFit], type = "l", col = "red", lwd = 2, lty = 1)
			lines(1:nmydata, arima.mre[nperiodsFit, ], type = "l", col = "blue", xlab = "", ylab = "", lwd=2, lty = 2)
			lines(1:nperiodsFit, arima.mre[nperiodsFit, 1:nperiodsFit], type = "l", col = "blue", xlab = "", ylab = "", lwd=2, lty = 1)

			rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

			legend("topleft", c("NULL", paste0(epi_model[nregions2], " aggregate"), 'Incidence'), text.col = c("red", "blue", 'black'), bty = "n")

			legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
			axis(1, at = 1:nmydata, labels = FALSE)
			axis(1, at = ind, labels = dates, las = 2)
   			par(new=TRUE)
    		plot(1:nmydata, observed, type = "n", col = "black", lwd = 2, xlab = '', ylab = '', main = '', xlim = c(1, nmydata), xaxt = "n", yaxt='n')
			lines(1:nmydata, observed, type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")

 			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")
 			lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)


		}

		## Plot direct fit to National

		observed = tables.mod$epi.obsrv[, mod_id]
		epi.null = tables.mod$epi.nul[  , mod_id]
 		arima = tables.mod$arima.frcst[, , mod_id]

 		null.mae = tables.mod$epi.null.mae[, , mod_id]
		null.mre = tables.mod$epi.null.mre[, , mod_id]

 		arima.mae = tables.mod$arima.mae[, , mod_id]
 		arima.mre = tables.mod$arima.mre[, , mod_id]
		arima.rel = tables.mod$arima.rel[, , mod_id]


		observed[observed == 0] <- NA
		if (nperiodsData < nperiods) {
			null.mae[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
			null.mre[nperiodsFit, (nperiodsData + 1):nperiods] <- NA

			arima.mae[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
			arima.mre[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
			arima.rel[nperiodsFit, (nperiodsData + 1):nperiods] <- NA
		}

		index <- is.nan(null.mre[nperiodsFit, ]	)
		null.mre[nperiodsFit, index ] <- NA
		index <- is.infinite(null.mre[nperiodsFit, ])
		null.mre[nperiodsFit, index ] <- NA

		index <- is.nan(arima.mre[nperiodsFit, ]	)
		arima.mre[nperiodsFit, index ] <- NA
		index <- is.infinite(arima.mre[nperiodsFit, ])
		arima.mre[nperiodsFit, index ] <- NA


    	title = paste0(mydata$model$name, ": ", mydata$FY, " season")

    	xlab = ""

    	ymax = max(arima.mre[nperiodsFit, ],null.mre[nperiodsFit, ], na.rm=TRUE)
		ymin = min(arima.mre[nperiodsFit, ],null.mre[nperiodsFit, ], na.rm=TRUE)

    	xlab = paste0("Time (", cadence, ")")

		plot(1:nmydata, null.mre[nperiodsFit, ], type = "n", col = "red", lwd = 2, xlab = xlab, ylab = "Mean Relative Error", ylim = c(ymin, ymax), main = title, xlim = c(1, nmydata), xaxt = "n")

		lines(1:nmydata, null.mre[nperiodsFit, ], type = "l", col = "red", lwd = 2, lty = 2)
		lines(1:nperiodsFit, null.mre[nperiodsFit, 1:nperiodsFit], type = "l", col = "red", lwd = 2, lty = 1)
		lines(1:nmydata, arima.mre[nperiodsFit, ], type = "l", col = "blue", xlab = "", ylab = "", lwd=2, lty = 2)
		lines(1:nperiodsFit, arima.mre[nperiodsFit, 1:nperiodsFit], type = "l", col = "blue", xlab = "", ylab = "", lwd=2, lty = 1)

		rect(xleft = nperiodsFit, ybottom = 0, xright = (nperiodsFit + 4), ytop = (ymax * 1.2), col = rgb(0.1, 0.1, 0.1, alpha = 0.1), border = NA)

    	legend("topleft", c("NULL", paste0(epi_model[nregions1], " direct"), "Incidence"), text.col = c("red", "blue", 'black'), bty = "n")
    	legend("topright", paste0("nperiodsFit=", nperiodsFit), bty = "n")
    	axis(1, at = 1:nmydata, labels = FALSE)
    	axis(1, at = ind, labels = dates, las = 2)
    	par(new=TRUE)
    	plot(1:nmydata, observed, type = "n", col = "black", lwd = 2, xlab = '', ylab = '', main = '', xlim = c(1, nmydata), xaxt = "n", yaxt='n')
		lines(1:nmydata, observed, type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")

 		lines(1:nperiodsFit, observed[1:nperiodsFit], type = "l", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "")
 		lines(1:nperiodsFit, observed[1:nperiodsFit], type = "p", col = "black", lwd = 1, lty = 1, xlab = "", ylab = "", pch = 19)

    dev.off()

    return(err = 0)
}



plotMCMC <- function(mydata=NULL, tab.model = NULL, tab.fit=NULL, tab.cpl = NULL, opt.list=NULL, opt.cpl = NULL, run.list = NULL, imask = NULL, ireal = 1, idevice = 1) {
    #' Plot Posterior Distribution from MCMC Procedure
    #'
    #' \code{plotMCMC} Creates a PDF file with plots of the posterior distributions of all the mechanistic model \
    #' parameters the User has chosen to fit and of the AICc score
    #' @param mydata A dataframe with all the available data for this \pkg{DICE} run
    #' @param tab.model The MCMC history of the direct fit of the model data
    #' @param tab  The MCMC history of an indirect fit of the model using a coupled model.
    #'   This  array includes all the parameters except the two that define the coupling matrix.
    #' @param tab.cpl The MCMC history of the two parameters that help define the coupling matrix:
    #'   the saturation distance and the distance power.
    #' @param opt.list A logical list of all the parameters \pkg{DICE} recognizes and
    #'   can optimize with TRUE/FALSE
    #' @param opt.cpl A logical list for the two parameters that help define the coupling matrix.
    #' @param run.list a list with parameters used for the MCMC procedure
    #' @param imask AN array of integers with +1/-1 values for parameters that are optimized (or not)
    #' @param ireal - Integer, the MCMC chain number
    #' @examples
    #' plotMCMC(mydata=mydata, Tg = Tg, tab.model = tab.model, tab.fit= tab.fit, tab.cpl = tab.cpl,
    #' opt.list=opt.list, opt.cpl = opt.cpl, run.list = run.list, imask = imask, ireal = ireal)
    #' @return err=0 if plots were created
    #'
    #'
    #'

    device = run.list$device[idevice]
    if (is.null(device))
        device = "png"
    if (is.null(tab.model) && is.null(tab.fit))
        return

    subDir = run.list$subDir

    if (!dir.exists(subDir)) {
        dir.create(subDir)
    }

    ithin = run.list$ithin
    nlines = run.list$nlines
    nMCMC = run.list$nMCMC

    myColName = names(opt.list)

    nopt = length(which(imask == 1))
    names.opt = myColName[which(imask == 1)]

    model = mydata$imodel

	FY = mydata$FY

	nperiods = mydata$nperiods
	nperiodsFit = mydata$nperiodsFit
	nregions = mydata$fit$nregions

    # how many steps to burn - here we set it to half which is very rigid

    ## First operate on tab.mod
    if (!is.null(tab.model)) {
    	tab = tab.model
		colnames(tab) = c(myColName, "AICc")
    	nlines <- dim(tab)[1]

    	nparam <- dim(tab)[2] - 1

    	iburn <- nlines/2

    	nMCMC = nlines * ithin
    	# This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics

    	results = mcmc(data = tab[iburn:nlines, ], start = (iburn * ithin), end = nMCMC, thin = ithin)

    	nr = nc = 3

    	myName = mydata$dataName
		myName = gsub(" ","",myName)
    	if (tolower(device) == "pdf") {
    		filename = paste0(subDir, "/params-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
    		pdf(file = filename, onefile = TRUE, width = 11, height = 11)
    	} else if (tolower(device) == "png") {
    		filename = paste0(subDir, "/params-", myName, "-", nperiodsFit, "-", ireal, ".png", sep = "")
    		png(file = filename, width = 1100, height = 1100)
    	} else {
    		dev.next()
    		dev.new()
    	}

    	cat("\nFor plot of MCMC posterior densities see: ", filename, "\n")

    	par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

    	legend = rep(0, 2)


    	for (ivar in 1:nparam) {
    		if (imask[ivar] < 0)
    			next
    		myvar = myColName[ivar]
    		icol = which(colnames(tab) == myvar)

    		densplot(results[, icol], ylab = "Frequency", xlab = myvar, main = mydata$model$name, col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

    		z = fitdistr(results[, icol], "normal")
    		val.mean = z$estimate["mean"]
    		val.sd = z$estimate["sd"]
    		legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
    		legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
    		legend("topleft", legend = legend, bty = "n")
    		legend("topright",myvar, bty = 'n')
    	}

    	## repeat for the AICc

    	icol = nparam + 1
    	myvar = "AICc"
    	densplot(results[, icol], ylab = "Frequency", xlab = myvar, main = mydata$model$name, col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

    	z = fitdistr(results[, icol], "normal")
    	val.mean = z$estimate["mean"]
    	val.sd = z$estimate["sd"]
    	legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
    	legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
    	legend("topleft", legend = legend, bty = "n")
		legend("topright", legend = myvar, bty = 'n')

    }


	if (is.null(tab.fit)) {
		dev.off()
		return(err = 0)
	}

	nskip = (nopt+1) %% (nc)
	nskip = (nc - nskip)

	for (i in 1:nskip) plot.new()


    if (!is.null(tab.fit)) {

    	for (iregion in 1:nregions) {
    		
    		tab = tab.fit[[iregion]]

			colnames(tab) = c(myColName, "AICc")
    		# This mcmc object is only from 'iburn' and up and used only for the purpose of the statistics

    		results = mcmc(data = tab[iburn:nlines, ], start = (iburn * ithin), end = nMCMC, thin = ithin)

    		for (ivar in 1:nparam) {
    			if (imask[ivar] < 0)
    				next
    			myvar = myColName[ivar]
    			icol = which(colnames(tab) == myvar)

    			densplot(results[, icol], ylab = "Frequency", xlab = myvar, main = mydata$fit$name[iregion], col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

    			z = fitdistr(results[, icol], "normal")
    			val.mean = z$estimate["mean"]
    			val.sd = z$estimate["sd"]
    			legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
    			legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
    			legend("topleft", legend = legend, bty = "n")
    			legend("topright", legend = myvar, bty = 'n')
    		}

    		## repeat for the AICc

    		icol = nparam + 1
    		myvar = "AICc"
    		densplot(results[, icol], ylab = "Frequency", xlab = myvar, main = mydata$fit$name[iregion], col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

    		z = fitdistr(results[, icol], "normal")
    		val.mean = z$estimate["mean"]
    		val.sd = z$estimate["sd"]
    		legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
    		legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
    		legend("topleft", legend = legend, bty = "n")
			legend("topright", legend = myvar, bty= 'n')

			for (i in 1:nskip) plot.new()

    	}
    }

   if (!is.null(tab.cpl)) {
   	tab = tab.cpl
   	ncpl = dim(tab)[2]
   	colnames(tab) = names(opt.cpl)

   	results = mcmc(data = tab[iburn:nlines, ], start = (iburn * ithin), end = nMCMC, thin = ithin)

   	for (ivar in 1:ncpl) {
   		myvar = names(opt.cpl)[ivar]
   		icol = which(colnames(tab) == myvar)

   		densplot(results[, icol], ylab = "Frequency", xlab = myvar, main = "", col = "blue", show.obs = FALSE, xlim = range(tab[, icol]))

   		z = fitdistr(results[, icol], "normal")
   		val.mean = z$estimate["mean"]
   		val.sd = z$estimate["sd"]
   		legend[1] = paste("Mean", formatC(signif(val.mean, digits = 3), digits = 3, format = "fg", flag = "#"))
   		legend[2] = paste("SD", formatC(signif(val.sd, digits = 3), digits = 3, format = "fg", flag = "#"))
   		legend("topleft", legend = legend, bty = "n")
   		legend("topright", legend = myvar, bty = "n")
   	}


   }
    dev.off()

    return(err = 0)

}

