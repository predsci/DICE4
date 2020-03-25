###
### This file holds ALL the ggplot2 R functions
###

plotFitOnePatch.ggplot2 <- function(model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1, run.list = NULL, idevice = 1) {

    #' Plot the results of a \pkg{DICE} Run - Single Region
    #'
    #' Plot the results of \pkg{DICE} run for a single region/patch. We show the incidence along with our fits and
    #' if appropriate predictions. We show the best result and randomly selected results from the MCMC chain. This is the ggplot2 version.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the data
    #' @param mydata A dataframe with all the data available for this \code{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param idevice Integer - the index of the device in the device array. Default is 1 - make only one format of plot results
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' plotFitOnePatch{model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, idevice = device}



    device = run.list$device[idevice]
    if (is.null(device))
        device = "png"
    if (is.null(model_profile))
        return
    if (is.null(device))
        device = "x11"
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

    ## This model raw data
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

    ## Determine ylabel based on dataType
    if (mydata$data_source == "cdc" | mydata$data_source == "gft") {
        ylab = "%ILI"
    } else {
        ylab = "# Cases"
    }

    # Plot the data and fits/predictions
    plotlist = list()  # For saving all the plots
    if (!is.null(model_rtn) && !is.null(model_profile)) {
        model_mean = rep(0, nperiods)
        for (iweek in 1:nperiods) model_mean[iweek] = mean(model_profile_ili[, iweek])
        ymax = max(model_rtn_ili[1:nperiodsData], model_profile_ili[, 1:nperiodsData], model_ili[1:nperiodsData], na.rm = TRUE)

        breaks = seq(from = 1, to = nperiods, by = 4)
        labels = weeks[breaks]

        plotlist[[1]] = ggplot(data = NULL) + scale_x_continuous(name = "EW #", limits = c(1, nperiods), breaks = breaks, labels = labels) +
            scale_y_continuous(name = ylab, limits = c(0, ymax)) + theme(text = element_text(size = 10, color = "gray20", face = "italic"),
            axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

        step = max(1, nRnd/100)
        irnd.set = seq(from = 1, to = nRnd, by = step)
        dat.rnd = t(model_profile[irnd.set, 1:nperiodsFit])
        dat.rnd.pred = t(model_profile_ili[irnd.set, nperiodsFit:nperiods])
        data.rnd = melt(dat.rnd)
        data.rnd.pred = melt(dat.rnd.pred)

        plotlist[[1]] = plotlist[[1]] + geom_line(aes(x = data.rnd[, 1], y = data.rnd[, 3], group = data.rnd[, 2]), col = "#E495A5",
            size = 2, alpha = 0.4) + geom_line(aes(x = 1:nperiodsFit, y = model_rtn_ili[1:nperiodsFit]), col = "#39BEB1", size = 1) + geom_line(aes(x = 1:nperiodsFit,
            y = model_mean[1:nperiodsFit]), col = "#099DD7", size = 0.8) + geom_line(aes(x = 1:nperiods, y = model_ili), col = "black", na.rm = TRUE) +
            geom_point(aes(x = 1:nperiodsFit, y = model_ili[1:nperiodsFit]), col = "#24576D", size = 1, na.rm = TRUE)

        if (nperiodsFit < nperiods) {
            plotlist[[1]] = plotlist[[1]] + geom_line(aes(x = (data.rnd.pred[, 1] + nperiodsFit - 1), y = data.rnd.pred[, 3], group = data.rnd.pred[,
                2]), col = "#E495A5", size = 2, linetype = 2, alpha = 0.4) + geom_line(aes(x = nperiodsFit:nperiods, y = model_mean[nperiodsFit:nperiods]),
                col = "#099DD7", size = 0.8, linetype = 2) + geom_line(aes(x = nperiodsFit:nperiods, y = model_rtn_ili[nperiodsFit:nperiods]),
                col = "#39BEB1", size = 1, linetype = 2) + geom_rect(aes(xmin = nperiodsFit, xmax = min(nperiodsFit + 4, nperiods), ymin = 0,
                ymax = ymax), fill = "#D497D3", alpha = 0.7)
        }

        if (length(model_onset) > 0) {
            plotlist[[1]] = plotlist[[1]] + geom_hline(yintercept = model_onset, col = "#D497D3", size = 1, linetype = 2)

        }

        reg.name = paste("   ", mydata$model$name, c("-Data", "-Model-Best", "-Model-Mean", "-Model-Random"), sep = "")
        plotlist[[1]] = plotlist[[1]] + annotate("text", x = rep(-Inf, 5), y = rep(Inf, 5), label = c(paste("   ", mydata$FY, sep = ""),
            reg.name), hjust = rep(0, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black", "black", "#39BEB1", "#099DD7",
            "#E495A5"), family = "serif", size = 3.5)
        # The following two blocks have not been tested
        if (any(model == c(2, 3, 101, 102, 104, 105))) {
            school = mydata$model$school
            school[school == 0] = NA
            plotlist[[1]] = plotlist[[1]] + geom_point(aes(x = 1:nperiods, y = school * (ymax/5)), fill = "grey50", col = "grey", size = 3,
                pch = 22, na.rm = TRUE)
        }
        if (any(model == c(1, 3, 102, 103, 104, 105))) {
            sh = mydata$model$sh
            plotlist[[1]] = plotlist[[1]] + geom_line(aes(x = 1:nperiods, y = sh * (ymax/max(sh))), col = "black")
        }
    }

    # maximum week in data
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
    ylab = "Probability Density"
    xlab = "EW #"

    # Plot histograms of maximum week in model
    if (!is.null(model_profile_ili)) {

        breaks = seq(from = wk.min, to = wk.max, by = 1)
        breaks_x = seq(from = wk.min, to = wk.max, by = 4)
        labels = weeks[breaks_x]
        data = data.frame(x = c(dat_model_wk_max, drct_model_wk_max), y = c(rep(0, length(dat_model_wk_max)), rep(1, length(drct_model_wk_max))))
        plotlist[[2]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..), breaks = breaks,
            fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..), breaks = breaks,
            fill = "deeppink", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks), breaks = breaks_x,
            labels = labels) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10, color = "gray20", face = "italic"),
            axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))
    }

    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    leg.text = c(mydata$FY, model.name, "Data", "Model")
    plotlist[[2]] = plotlist[[2]] + annotate("text", x = c(rep(-Inf, 2), rep(Inf, 4)), y = rep(Inf, 6), label = c(paste("   ", "Observed/Predicted",
        sep = ""), paste("   ", "Peak Week", sep = ""), leg.text), hjust = c(0, 0, 1, 1, 1, 1), vjust = c(2.5, 4, seq(from = 2.5, to = 7,
        by = 1.5)), col = c("black", "black", "black", "black", "dodgerblue", "deeppink"), family = "serif", size = 3.5)

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
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                if (max.val <= 20)
                  step = 2 else if (max.val <= 40)
                  step = 4 else if (max.val <= 150)
                  step = 10 else step = 50
                data = data.frame(x = c(model_ili[iweek], model_profile_ili[, iweek]), y = c(rep(0, length(model_ili[iweek])), rep(1,
                  length(model_profile_ili[, iweek]))))
                if (!is.na(model_ili[iweek])) {
                  plotlist[[iweek - nperiodsFit + 3]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..),
                    breaks = breaks, fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..),
                    breaks = breaks, fill = "deeppink", col = "black", alpha = 0.7)

                } else {
                  plotlist[[iweek - nperiodsFit + 3]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..),
                    breaks = breaks, fill = "deeppink", col = "black", alpha = 0.7)
                }
                plotlist[[iweek - nperiodsFit + 3]] = plotlist[[iweek - nperiodsFit + 3]] + scale_x_continuous(name = xlab, limits = range(breaks),
                  breaks = seq(from = min.val, to = max.val, by = step)) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10,
                  color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

                model.name = mydata$model$name
                my.week = paste("EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(mydata$FY, model.name, "Data", "Model")
                plotlist[[iweek - nperiodsFit + 3]] = plotlist[[iweek - nperiodsFit + 3]] + annotate("text", x = c(rep(-Inf, 2), rep(Inf,
                  4)), y = rep(Inf, 6), label = c(paste("   ", "Observed/Predicted", sep = ""), paste("   %ILI for ", my.week, sep = ""),
                  leg.text), hjust = c(0, 0, 1, 1, 1, 1), vjust = c(2.5, 4, seq(from = 2.5, to = 7, by = 1.5)), col = c("black", "black",
                  "black", "black", "dodgerblue", "deeppink"), family = "serif", size = 3.5)
            } else {
                max.val = ceiling(max(model_profile_ili[, iweek], na.rm = TRUE))
                max.val = 2 * max.val
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                plotlist[[iweek - nperiodsFit + 3]] = ggplot(data = NULL) + geom_histogram(aes(model_profile_ili[, iweek], y = ..density..),
                  fill = "dodgerblue", col = "black", breaks = breaks, alpha = 0.7) + scale_x_continuous(name = xlab, limits = c(breaks[1],
                  breaks[length(breaks)]), breaks = seq(min.val, max.val, by = 2)) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10,
                  color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

                model.name = mydata$model$name
                my.week = paste("EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(mydata$FY, model.name, "Model")
                plotlist[[iweek - nperiodsFit + 3]] = plotlist[[iweek - nperiodsFit + 3]] + annotate("text", x = c(rep(-Inf, 2), rep(Inf,
                  3)), y = rep(Inf, 5), label = c(paste("   ", "Observed/Predicted", sep = ""), paste("   %ILI for ", my.week, sep = ""),
                  leg.text), hjust = c(0, 0, 1, 1, 1), vjust = c(2.5, 4, 2.5, 4, 5.5), col = c("black", "black", "black", "black", "dodgerblue"),
                  family = "serif", size = 3.5)
            }
        }
    }
    layout = matrix(seq(1, 9, 1), nrow = 3, byrow = TRUE)
    multiplot(plotlist = plotlist, layout = layout)

    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list, idevice = idevice, model_profile_ili = model_profile_ili, model_rtn_ili = model_rtn_ili)
    
    filename = paste(subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)



    filename = paste(subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)



    # Make a movie
    if ((mydata$model$level == 2) && (mydata$fit$level == 3)) {
        ## Get map-data from GADM and simplify it
        port1 = suppressMessages(getData("GADM", country = "USA", level = 1))
        port1$NAME_1 = as.factor(as.character(port1$NAME_1))
        name = port1$NAME_1
        port1 = gSimplify(port1, tol = 0.01, topologyPreserve = TRUE)

        ## fills the map basing on different region
        region_list = list(Region0 = c("Maine", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "Vermont"), Region1 = c("New York",
            "New Jersey", "Puerto Rico"), Region2 = c("Pennsylvania", "Delaware", "Maryland", "West Virginia", "Virginia", "District of Columbia"),
            Region3 = c("Kentucky", "Tennessee", "North Carolina", "South Carolina", "Georgia", "Florida", "Alabama", "Mississippi"),
            Region4 = c("Minnesota", "Wisconsin", "Illinois", "Indiana", "Michigan", "Ohio"), Region5 = c("New Mexico", "Texas", "Oklahoma",
                "Arkansas", "Louisiana"), Region6 = c("Nebraska", "Kansas", "Iowa", "Missouri"), Region7 = c("Utah", "Colorado", "Wyoming",
                "Montana", "South Dakota", "North Dakota"), Region8 = c("California", "Nevada", "Arizona", "Hawaii"), Region9 = c("Oregon",
                "Washington", "Idaho", "Alaska"))
        col = rainbow(length(region_list))
        labels = c("Region0", "Region1", "Region2", "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Region9")
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

        labels = c("Region1", "Region2", "Region3", "Region4", "Region5", "Region6", "Region7", "Region8", "Region9", "Region10")
        map1.df$col = state.col[map1.df$id]
        map1.df$col = factor(map1.df$col, levels = unique(state.col)[order(unique(state.label))], labels = labels)

        # mydata = get.DICE.data(dataType = dataType, year = year, mod_level = 2, fit_level = 3, model = model, nperiodsFit = nperiodsFit,
        # isingle=isingle)
        longitude = mydata$fit$attr$lon
        latitude = mydata$fit$attr$lat
        # output = runDICE(dataType = dataType, year = year, mod_level = 2, fit_level = 3, model = model, isingle = isingle, nMCMC = nMCMC,
        # nreal = nreal) onset = output$rtn
        factor = as.numeric(mydata$fit$factor)
        region_r = rtn
        for (i in 1:ncol(rtn)) {
            region_r[, i] = rtn[, i]/factor[i]
        }

        nation_long = suppressMessages(geocode("USA")$lon)
        nation_lat = suppressMessages(geocode("USA")$lat)
        # region_r = onset

        ## Normailization: divided by minimum value
        for (i in 1:ncol(region_r)) {
            region_r[, i] = region_r[, i]/min(region_r[, i])
        }
        nation_r = apply(region_r, 1, sum)

        week = mydata$weeks
        # ## Make a movie mainDir = getwd() err = makeDir(subDir = subDir) setwd(paste(mainDir,subDir,sep = '/'))
        # frame = nrow(region_r)

        # for (i in 1:frame) {
            # cat("making frame number: ", i, "\n")
            # if (i < 10) {
                # name = paste(subDir, "/000", i, "plot.png", sep = "")
            # }
            # if (i >= 10) {
                # name = paste(subDir, "/00", i, "plot.png", sep = "")
            # }
            # p = ggplot() + geom_map(data = map1.df, map = map1.df, aes(map_id = id, x = long, y = lat, group = group, fill = col), size = 0.25) +
                # coord_map() + geom_point(aes(x = nation_long, y = nation_lat), size = 2 * nation_r[i], col = "grey", show.legend = FALSE,
                # alpha = 0.5) + scale_x_continuous(name = "", limits = c(-130, -60)) + scale_y_continuous(name = "", limits = c(25, 50)) +
                # ggtitle(paste("USA map colored by CDC region - Epidemic Week ", weeks[i], sep = "")) + theme(legend.title = element_blank(),
                # plot.title = element_text(size = 15, family = "serif", face = "bold"))
            # for (j in 1:10) {
                # p = p + geom_point(aes_string(x = longitude[j], y = latitude[j]), size = 4 * region_r[i, j], col = "#24576D", show.legend = FALSE,
                  # alpha = 0.8)
            # }
            # suppressMessages(ggsave(filename = name))
        # }

        # moviename = paste("map-", myName, "-", nperiodsFit, "-", ireal, ".mov", sep = "")
        # command_to_convert = paste("convert -delay 20 *.png ", moviename, sep = "")
        # system('convert -delay 20 *.png map.gif') system(command_to_convert) file.remove(list.files(pattern = '*.png'))

        # cat(' For movie see: ',subDir,'/',moviename, sep = '')
    }



    err = 0
    return(err)


}

plotFitCDCPercentILI.ggplot2 <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1,
    run.list = NULL, idevice = 1) {

    #' Plot the results of a \pkg{DICE} Run Using ggplot2
    #'
    #' Plot the results of an a coupled or uncoupled  \pkg{DICE} run. For each of the fit regions we plot  the disease incidence for
    #' the region along with our predictions for it based on randomly selected results from the  history of the MCMC chain of each region.
    #' Using the predictions
    #' for the fit regions we then show the results for the model  region as a weighted sum of the fit regions.  The last panel
    #' shows our prediction for the model region using a direct fit to the model data.  The function also writes a binary RData file with
    #' all the profile predictions for the model and fit regions. Note that in the case of a coupled run the fit regions are never individually
    #' optimized. It is their weighted sum that is optimized, with the weights given by the relative population of each fit region.
    #' @param rtn A 1D numeric array with the best indirect prediction to the model region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @param device String with format for output of plots - pdf or png
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
    n.model = 1
    n.fit = mydata$fit$nregions
    nRnd = dim(profile)[1]

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
    if (tolower(device) == "pdf") {
        pdfName = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, ".pdf")
        cat("\n\n For a plot of %ILI See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 15, height = 9)
    } else if (tolower(device) == "png") {
        pngName = paste0(subDir, "/results-", myName, "-", nperiodsFit, "-", ireal, "_pg")
        cat("\n\n For a plot of  %ILI  See: ", pngName, "\n\n")
        png(file = paste0(pngName, "%01d.png"), width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }

    # convert the fit and model results to %ILI calculate the national

    # convert the model to %ILI calculate the national
    factor = as.numeric(mydata$fit$factor)
    model_factor = mydata$model$factor

    ## This is the fit and model %ILI data
    fit_ili = mydata$fit$raw
    model_ili = mydata$model$raw

    rtn_ili = rtn
    profile_ili = profile
    for (i in 1:n.fit) {
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

    ## for plotting - replace zero's with NA
    if (nperiodsFit < nperiods) {
        index <- which(fit_ili[, 1] == 0)

        if (length(index) >= 1) {
            fit_ili[index, 1:n.fit] = NA
            model_ili[index] = NA
        }
    }

    colnames(fit_ili) = mydata$fit$attr$NAME_3
    colnames(rtn_ili) = mydata$fit$attr$NAME_3

    colvec = rainbow(n.fit)
    lwd = rep(2, n.fit)

    # Plot one region at a time
    plot_region = list()
    for (i in 1:n.fit) {
        ymax = max(fit_ili[1:nperiodsData, i], rtn_ili[1:nperiodsData, i], profile_ili[, 1:nperiodsData, i], unlist(fit_onset[i]), na.rm = TRUE)
        breaks = seq(from = 1, to = nperiods, by = 4)
        labels = weeks[breaks]
        plot_region[[i]] = ggplot(data = NULL) + scale_x_continuous(name = "EW #", limit = c(1, nperiods), breaks = breaks, labels = labels) +
            scale_y_continuous(name = "% ILI") + coord_cartesian(ylim = c(0, ymax)) + theme(text = element_text(size = 10, color = "gray20",
            face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

        step = max(1, nRnd/100)
        irnd.set = seq(from = 1, to = nRnd, by = step)
        dat.rnd = t(profile_ili[irnd.set, 1:nperiodsFit, i])
        dat.rnd.pred = t(profile_ili[irnd.set, nperiodsFit:nperiods, i])
        data.rnd = melt(dat.rnd)
        data.rnd.pred = melt(dat.rnd.pred)
        plot_region[[i]] = plot_region[[i]] + geom_line(aes_string(x = data.rnd[, 1], y = data.rnd[, 3], group = data.rnd[, 2]), col = colvec[i],
            size = 1.5, alpha = 0.4) + geom_line(aes_string(x = 1:nperiodsFit, y = rtn_ili[1:nperiodsFit, i]), col = colvec[i]) + geom_line(aes_string(x = 1:nperiodsFit,
            y = fit_ili[1:nperiodsFit, i]), col = "black", linetype = 2, na.rm = TRUE) + geom_point(aes_string(x = 1:nperiodsFit, y = fit_ili[1:nperiodsFit,
            i]), col = "black", pch = 20)

        if (nperiodsFit < nperiods) {
            plot_region[[i]] = plot_region[[i]] + geom_line(aes_string(x = data.rnd.pred[, 1] + nperiodsFit - 1, y = data.rnd.pred[, 3],
                group = data.rnd.pred[, 2]), col = colvec[i], size = 1.5, linetype = 2, alpha = 0.4) + geom_line(aes_string(x = nperiodsFit:nperiods,
                y = rtn_ili[nperiodsFit:nperiods, i]), col = colvec[i], linetype = 2) + geom_rect(aes_string(xmin = nperiodsFit, xmax = min(nperiodsFit +
                4, nperiods), ymin = 0, ymax = Inf), fill = "#D497D3", alpha = 0.7)
        }
        if (length(fit_onset[i]) > 0) {
            y = rep(as.numeric(fit_onset[i]), nperiods)
            plot_region[[i]] = plot_region[[i]] + geom_line(aes_string(x = 1:nperiods, y = y), col = "#D497D3", size = 1, linetype = 2)
        }
        reg.name = reg.fit.name[i]
        rel.pop = fit_coef[i]
        rel.pop = round(rel.pop, digits = 3)
        plot_region[[i]] = plot_region[[i]] + annotate("text", x = rep(-Inf, 3), y = rep(Inf, 3), label = c(paste("   ", mydata$FY, sep = ""),
            paste("   ", reg.name, sep = ""), paste("   ", rel.pop, sep = "")), hjust = rep(0, 3), vjust = c(2.5, 4, 5.5), col = c("black",
            colvec[i], colvec[i]), family = "serif", size = 3.5)
        if (any(model == c(2, 3, 101, 102, 104, 105))) {
            school = mydata$fit$school[, i]
            school[school == 0] = NA
            plot_region[[i]] = plot_region[[i]] + geom_point(aes(x = 1:nperiods, y = school * (ymax/5)), fill = "grey50", col = "grey",
                size = 3, pch = 22, na.rm = TRUE)
        }
        if (any(model == c(1, 3, 102, 103, 104, 105))) {
            sh = mydata$fit$sh[, i]
            plot_region[[i]] = plot_region[[i]] + geom_line(aes(x = 1:nperiods, y = sh * (ymax/max(sh))), col = "black")
        }
    }

    ## now do the national
    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {
        fit_model[i] = sum(rtn_ili[i, 1:n.fit] * fit_coef[1:n.fit])
        tmp = rep(0, n.fit)
        for (k in 1:n.fit) tmp[k] = mean(profile_ili[, i, k])
        fit_model_mean[i] = sum(tmp[1:n.fit] * fit_coef[1:n.fit])
        for (irnd in 1:nRnd) {
            for (k in 1:n.fit) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile_ili[irnd, i, k] * fit_coef[k]
        }
    }

    ymax = max(fit_model[1:nperiodsData], fit_model_mean[1:nperiodsData], fit_model_profile[1:nperiodsData], model_ili[1:nperiodsData], na.rm = TRUE)
    breaks = seq(from = 1, to = nperiods, by = 4)
    labels = weeks[breaks]
    plot_national = ggplot(data = NULL) + scale_x_continuous(name = "EW #", limits = c(1, nperiods), breaks = breaks, labels = labels) +
        scale_y_continuous(name = "% ILI") + coord_cartesian(ylim = c(0, ymax)) + theme(text = element_text(size = 10, color = "gray20",
        face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))
    step = max(1, nRnd/100)
    irnd.set = seq(from = 1, to = nRnd, by = step)
    dat.rnd = t(fit_model_profile[irnd.set, 1:nperiodsFit])
    dat.rnd.pred = t(fit_model_profile[irnd.set, nperiodsFit:nperiods])
    data.rnd = melt(dat.rnd)
    data.rnd.pred = melt(dat.rnd.pred)
    plot_national = plot_national + geom_line(aes(x = data.rnd[, 1], y = data.rnd[, 3], group = data.rnd[, 2]), col = "#E495A5", size = 2,
        alpha = 0.4) + geom_line(aes(x = 1:nperiodsFit, y = fit_model[1:nperiodsFit]), col = "#39BEB1", size = 1) + geom_line(aes(x = 1:nperiodsFit,
        y = fit_model_mean[1:nperiodsFit]), col = "#099DD7", size = 0.8) + geom_line(aes(x = 1:nperiods, y = model_ili), col = "black", na.rm = TRUE) +
        geom_point(aes(x = 1:nperiodsFit, y = model_ili[1:nperiodsFit]), col = "#24576D", size = 1)

    if (nperiodsFit < nperiods) {
        plot_national = plot_national + geom_line(aes(x = data.rnd.pred[, 1] + nperiodsFit - 1, y = data.rnd.pred[, 3], group = data.rnd.pred[,
            2]), col = "#E495A5", size = 2, linetype = 2, alpha = 0.4) + geom_line(aes(x = nperiodsFit:nperiods, y = fit_model[nperiodsFit:nperiods]),
            col = "#39BEB1", size = 1, linetype = 2) + geom_line(aes(x = nperiodsFit:nperiods, y = fit_model_mean[nperiodsFit:nperiods]), col = "#099DD7",
            size = 0.8, linetype = 2) + geom_rect(aes(xmin = nperiodsFit, xmax = min(nperiodsFit + 4, nperiods), ymin = 0, ymax = Inf), fill = "#D497D3",
            alpha = 0.7)
    }

    if (length(model_onset) > 0) {
        plot_national = plot_national + geom_hline(yintercept = model_onset, col = "#D497D3", size = 1, linetype = 2)
    }
    reg.name = paste("   ", mydata$model$name, c("-Data", "-Best", "-Mean", "-Random"), sep = "")
    plot_national = plot_national + annotate("text", x = rep(-Inf, 5), y = rep(Inf, 5), label = c(paste("   ", mydata$FY, sep = ""),
        reg.name), hjust = rep(0, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black", "black", "#39BEB1", "#099DD7", "#E495A5"),
        family = "serif", size = 3.5)

    if (model == 2 || model == 3) {
        school = mydata$model$school
        school[school == 0] = NA
        plot_national = plot_national + geom_point(aes(x = 1:nperiods, y = school * (ymax/5)), fill = "grey50", col = "grey", size = 3,
            pch = 22, na.rm = TRUE)
    }

    if (model == 1 || model == 3) {
        sh = mydata$model$sh
        plot_national = plot_national + geom_line(aes(x = 1:nperiods, y = sh * (ymax/max(sh))), col = "black")
    }
    plot_region[[(n.fit + 1)]] = plot_national

    ## Repeat with the direct fit of the national data - if it was done!

    if (!is.null(model_rtn) && !is.null(model_profile)) {
        model_mean = rep(0, nperiods)
        for (iweek in 1:nperiods) model_mean[iweek] = mean(model_profile_ili[, iweek])
        ymax = max(model_rtn_ili[1:nperiodsData], model_profile_ili[, 1:nperiodsData], model_ili[1:nperiodsData], na.rm = TRUE)
        breaks = seq(from = 1, to = nperiods, by = 4)
        labels = weeks[breaks]
        plot_direct = ggplot(data = NULL) + scale_x_continuous(name = "EW #", limits = c(1, nperiods), breaks = breaks, labels = labels) +
            scale_y_continuous(name = "% ILI") + coord_cartesian(ylim = c(0, ymax)) + theme(text = element_text(size = 10, color = "gray20",
            face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

        step = max(1, nRnd/100)
        irnd.set = seq(from = 1, to = nRnd, by = step)
        dat_rnd = t(model_profile_ili[irnd.set, 1:nperiodsFit])
        dat_rnd_pred = t(model_profile_ili[irnd.set, nperiodsFit:nperiods])
        data_rnd = melt(dat_rnd)
        data_rnd_pred = melt(dat_rnd_pred)
        plot_direct = plot_direct + geom_line(aes(x = data_rnd[, 1], y = data_rnd[, 3], group = data_rnd[, 2]), col = "#E495A5", size = 2,
            alpha = 0.4) + geom_line(aes(x = 1:nperiodsFit, y = model_rtn_ili[1:nperiodsFit]), col = "#39BEB1", size = 1) + geom_line(aes(x = 1:nperiodsFit,
            y = model_mean[1:nperiodsFit]), col = "#099DD7", size = 0.8) + geom_line(aes(x = 1:nperiods, y = model_ili), col = "black", na.rm = TRUE) +
            geom_point(aes(x = 1:nperiodsFit, y = model_ili[1:nperiodsFit]), col = "#24576D", size = 1)

        if (nperiodsFit < nperiods) {
            plot_direct = plot_direct + geom_line(aes(x = data_rnd_pred[, 1] + nperiodsFit - 1, y = data_rnd_pred[, 3], group = data_rnd_pred[,
                2]), col = "#E495A5", size = 2, linetype = 2, alpha = 0.4) + geom_line(aes(x = nperiodsFit:nperiods, y = model_rtn_ili[nperiodsFit:nperiods]),
                col = "#39BEB1", size = 1, linetype = 2) + geom_line(aes(x = nperiodsFit:nperiods, y = model_mean[nperiodsFit:nperiods]), col = "#099DD7",
                size = 0.8, linetype = 2) + geom_rect(aes(xmin = nperiodsFit, xmax = min(nperiodsFit + 4, nperiods), ymin = 0, ymax = Inf),
                fill = "#D497D3", alpha = 0.7)
        }

        if (length(model_onset) > 0) {
            plot_direct = plot_direct + geom_hline(yintercept = model_onset, col = "#D497D3", size = 1, linetype = 2)
        }
        reg.name = paste(mydata$model$name, c("-Data", "-Direct-Best", "-Direct-Mean", "-Direct-Random"), sep = "")
        plot_direct = plot_direct + annotate("text", x = rep(-Inf, 5), y = rep(Inf, 5), label = c(paste("   ", mydata$FY, sep = ""),
            paste("   ", reg.name, sep = "")), hjust = rep(0, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black", "black",
            "#39BEB1", "#099DD7", "#E495A5"), family = "serif", size = 3.5)

        if (model == 2 || model == 3) {
            school = mydata$model$school
            school[school == 0] = NA
            plot_direct = plot_direct + geom_point(aes(x = 1:nperiods, y = school * (ymax/5)), fill = "grey50", col = "grey", size = 3,
                pch = 22, na.rm = TRUE)
        }

        if (model == 1 || model == 3) {
            sh = mydata$model$sh
            plot_direct = plot_direct + geom_line(aes(x = 1:nperiods, y = sh * (ymax/max(sh))), col = "black")
        }
    }
    plot_region[[(n.fit + 2)]] = plot_direct
    layout = matrix(seq(1, 12, 1), nrow = 3, byrow = TRUE)
    multiplot(plotlist = plot_region, layout = layout)

    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile, ireal = ireal, run.list = run.list, idevice = idevice, fit_ili = fit_ili, rtn_ili = rtn_ili, profile_ili = profile_ili,
        model_ili = model_ili, model_rtn_ili = model_rtn_ili, model_profile_ili = model_profile_ili, fit_model = fit_model, fit_model_mean = fit_model_mean, fit_model_profile = fit_model_profile)

    filename = paste(subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)

    err = 0
    return(err)
}

plotHists.ggplot2 <- function(rtn = NULL, profile = NULL, model_rtn = NULL, model_profile = NULL, mydata = NULL, ireal = 1, run.list = NULL,
    idevice = 1) {

    #' Plot Histograms of Predicted and Observed Peak Week and Value
    #'
    #' Plots two histogram files one for the observed and predicted peak week and the other for the
    #' observed and predicted \% ILI value for all weeks included in the range of the number of weeks fitted
    #' to the number of weeks of data.  The \% ILI is presented in bins of 0.5\%.The default is to
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
    #' plotHists.ggplot2(rtn = rtn, profile = profile, model_rtn = model_rtn, model_profile = model_profile,
    #' mydata = mydata, ireal = ireal, run.list = run.list, idevice = 1)

    device = run.list$device[idevice]

    if (is.null(device))
        device = "png"
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
    n.fit = mydata$fit$nregions
    nRnd = dim(profile)[1]
    nRndiceData = nRnd + 1
    n.fit1 = n.fit + 1

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
    for (i in 1:n.fit) {
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
            fit_ili[index, 1:n.fit] = NA
            model_ili[index] = NA
        }
    }

    colnames(fit_ili) = mydata$fit$attr$NAME_3
    colnames(rtn_ili) = mydata$fit$attr$NAME_3
    model_profile_indr = array(data = 0, dim = c(nRnd, nperiods))

    # this is the indirect modeling
    for (i in 1:nperiods) {
        for (j in 1:nRnd) {
            model_profile_indr[j, i] = sum(profile_ili[j, i, 1:n.fit] * fit_coef[1:n.fit])
        }
    }

    # direct fitting of n.fit regions and indirect of the model
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
    ## a histogram plot of binned data
    ylab = "Probability Density"
    if (mydata$data_source == "cdc" | mydata$data_source == "gft") {
        xlab = "%ILI"
    } else {
        xlab = "# Cases"
    }
    plotlist = list()
    plot_hist = list()
    for (iregion in 1:n.fit) {
        plot_hist[[iregion]] = list()
        for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
            min.val = 0
            if (iweek <= nperiodsData) {
                max.val = ceiling(max(fit_ili[iweek, iregion], profile_ili[, iweek, iregion], na.rm = TRUE))
                max.val = 2 * max.val
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                data = data.frame(x = c(fit_ili[iweek, iregion], profile_ili[, iweek, iregion]), y = c(rep(0, length(fit_ili[iweek, iregion])),
                  rep(1, length(profile_ili[, iweek, iregion]))))
                if (!is.na(fit_ili[iweek, iregion])) {
                  plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data,
                    y == 0), aes(y = ..density..), breaks = breaks, fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data,
                    y == 1), aes(y = ..density..), breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7)
                } else {
                  plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data,
                    y == 1), aes(y = ..density..), breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7)
                }
                plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] + scale_x_continuous(name = xlab,
                  limits = range(breaks), breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) +
                  theme(text = element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"),
                    axis.text.y = element_text(face = "plain"))
                leg.text = c(mydata$FY)
                reg.name = reg.fit.name[iregion]
                rel.pop = round(fit_coef[iregion], digits = 3)
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                leg.text = c(leg.text, reg.name, rel.pop, my.week, "Data", "Direct")
                plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] + annotate("text",
                  x = rep(Inf, 6), y = rep(Inf, 6), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 6), vjust = seq(from = 2.5,
                    to = 10, by = 1.5), col = c("black", "black", "black", "black", "dodgerblue", "forestgreen"), family = "serif", size = 3.5)
            } else {
                max.val = ceiling(max(profile_ili[, iweek, iregion], na.rm = TRUE))
                max.val = 2 * max.val
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = ggplot(data = NULL, aes_string(profile_ili[, iweek, iregion])) + geom_histogram(aes(y = ..density..),
                  breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
                  breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) +
                  theme(text = element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"),
                    axis.text.y = element_text(face = "plain"))
                leg.text = c(mydata$FY)
                reg.name = reg.fit.name[iregion]
                rel.pop = round(fit_coef[iregion], digits = 3)
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                leg.text = c(leg.text, reg.name, rel.pop, my.week, "Direct")
                plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] = plot_hist[[iregion]][[(iweek - nperiodsFit + 1)]] + annotate("text",
                  x = rep(Inf, 5), y = rep(Inf, 5), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 5), vjust = seq(from = 2.5,
                    to = 8.5, by = 1.5), col = c("black", "black", "black", "black", "forestgreen"), family = "serif", size = 3.5)
            }
        }
        plotlist = append(plotlist, plot_hist[[iregion]])
    }

    plot_dr = list()
    for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
        min.val = 0
        if (iweek <= nperiodsData) {
            max.val = ceiling(max(model_ili[iweek], model_profile_indr[, iweek]))
            max.val = 2 * max.val
            breaks = seq(from = min.val, to = max.val, by = 0.5)
            data = data.frame(x = c(model_ili[iweek], model_profile_indr[, iweek]), y = c(rep(0, length(model_ili[iweek])), rep(1, length(model_profile_indr[,
                iweek]))))
            plot_dr[[(iweek - nperiodsFit + 1)]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..),
                breaks = breaks, fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..),
                breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
                breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10,
                color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))
            leg.text = c(mydata$FY)
            model.name = mydata$model$name
            my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
            gsub(model.name, ".", " ", model.name)
            leg.text = c(leg.text, model.name, my.week, "Data", "Indirect")
            plot_dr[[(iweek - nperiodsFit + 1)]] = plot_dr[[(iweek - nperiodsFit + 1)]] + annotate("text", x = rep(Inf, 5), y = rep(Inf,
                5), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black",
                "black", "black", "dodgerblue", "forestgreen"), family = "serif", size = 3.5)
        } else {
            max.val = ceiling(max(model_profile_indr[, iweek]))
            max.val = 2 * max.val
            breaks = seq(from = min.val, to = max.val, by = 0.5)
            plot_dr[[(iweek - nperiodsFit + 1)]] = ggplot(data = NULL, aes_string(model_profile_indr[, iweek])) + geom_histogram(aes(y = ..density..),
                breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
                breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10,
                color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))
            leg.text = c(mydata$FY)
            model.name = mydata$model$name
            my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
            gsub(model.name, ".", " ", model.name)

            leg.text = c(leg.text, model.name, my.week, "Indirect")
            plot_dr[[(iweek - nperiodsFit + 1)]] = plot_dr[[(iweek - nperiodsFit + 1)]] + annotate("text", x = rep(Inf, 4), y = rep(Inf,
                4), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 4), vjust = seq(from = 2.5, to = 7, by = 1.5), col = c("black",
                "black", "black", "forestgreen"), family = "serif", size = 3.5)
        }
    }
    plotlist = append(plotlist, plot_dr)

    if (!is.null(model_profile_ili)) {
        plot_idr = list()
        for (iweek in nperiodsFit:min(nperiods, nperiodsFit + 4)) {
            if (iweek <= nperiodsData) {
                max.val = ceiling(max(model_ili[iweek], model_profile_ili[, iweek]))
                max.val = 2 * max.val
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                data = data.frame(x = c(model_ili[iweek], model_profile_ili[, iweek]), y = c(rep(0, length(model_ili[iweek])), rep(1,
                  length(model_profile_ili[, iweek]))))
                plot_idr[[(iweek - nperiodsFit + 1)]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..),
                  breaks = breaks, fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..),
                  breaks = breaks, fill = "mediumpurple", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
                  breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) +
                  theme(text = element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"),
                    axis.text.y = element_text(face = "plain"))

                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, my.week, "Data", "Direct")
                plot_idr[[(iweek - nperiodsFit + 1)]] = plot_idr[[(iweek - nperiodsFit + 1)]] + annotate("text", x = rep(Inf, 5), y = rep(Inf,
                  5), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 5), vjust = seq(from = 2.5, to = 8.5, by = 1.5), col = c("black",
                  "black", "black", "dodgerblue", "mediumpurple"), family = "serif", size = 3.5)
            } else {
                max.val = ceiling(max(model_profile_ili[, iweek]))
                max.val = 2 * max.val
                breaks = seq(from = min.val, to = max.val, by = 0.5)
                plot_idr[[(iweek - nperiodsFit + 1)]] = ggplot(data = NULL, aes_string(model_profile_ili[, iweek])) + geom_histogram(aes(y = ..density..),
                  breaks = breaks, fill = "mediumpurple", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
                  breaks = seq(min.val, max.val, by = (max.val - min.val)/4), oob = rescale_none) + scale_y_continuous(name = ylab) +
                  theme(text = element_text(size = 10, color = "gray20", face = "italic"), axis.text.x = element_text(face = "plain"),
                    axis.text.y = element_text(face = "plain"))

                leg.text = c(mydata$FY)
                model.name = mydata$model$name
                my.week = paste("%ILI for EW # ", weeks[iweek], sep = "")
                gsub(model.name, ".", " ", model.name)
                leg.text = c(leg.text, model.name, my.week, "Direct")
                plot_idr[[(iweek - nperiodsFit + 1)]] = plot_idr[[(iweek - nperiodsFit + 1)]] + annotate("text", x = rep(Inf, 4), y = rep(Inf,
                  4), label = paste(leg.text, "   ", sep = ""), hjust = rep(1, 4), vjust = seq(from = 2.5, to = 7, by = 1.5), col = c("black",
                  "black", "black", "mediumpurple"), family = "serif", size = 3.5)
            }
        }
    }

    plotlist = append(plotlist, plot_idr)
    s = ncol * nrow
    for (i in 1:ceiling(length(plotlist)/s)) {
        layout = matrix(seq(1, s, 1), nrow = 3, byrow = TRUE)
        if (i == ceiling(length(plotlist))/s) {
            multiplot(plotlist = plotlist[((i - 1) * s + 1):length(plotlist)], layout = layout)
        } else {
            multiplot(plotlist = plotlist[((i - 1) * s + 1):(i * s)], layout = layout)
        }
    }

    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    colvec = rainbow(n.fit)
    lwd = rep(2, n.fit)

    # Plot the individual fits along with the model fit

    if (tolower(device) == "pdf") {
        pdfName = paste(subDir, "/hist-week-max-", myName, "-", nperiodsFit, "-", ireal, ".pdf", sep = "")
        cat("\n\n For a Histogram Plot of Peak Week See: ", pdfName, "\n\n")
        pdf(file = pdfName, onefile = TRUE, width = 16, height = 12)

    } else if (tolower(device) == "png") {
        pngName = paste0(subDir, "/hist-week-max-", myName, "-", nperiodsFit, "-", ireal, "_pg")
        cat("\n\n For a Histogram Plot of Peak Week See: ", pngName, "\n\n")
        png(file = paste0(pngName, "%01d.png"), width = 1200, height = 900)
    } else {
        dev.next()
        dev.new()
    }

    nrow = 3
    ncol = 4

    # maximum week for model and fit in the data

    dat_model_wk_max = which.max(model_ili)

    if (mydata$fit$level > mydata$model$level) {
        dat_fit_wk_max = rep(0, n.fit)
        for (i in 1:n.fit) dat_fit_wk_max[i] = which.max(fit_ili[, i])
    }

    # maximum fit in direct modeling

    drct_model_wk_max = rep(0, nRnd)

    if (!is.null(model_profile_ili))
        for (i in 1:nRnd) drct_model_wk_max[i] = which.max(model_profile_ili[i, ])


    indrct_model_wk_max = rep(0, nRnd)
    for (i in 1:nRnd) indrct_model_wk_max[i] = which.max(model_profile_indr[i, ])

    sim_fit_wk_max = array(0, c(nRnd, n.fit))

    for (i in 1:n.fit) {
        for (j in 1:nRnd) {
            sim_fit_wk_max[j, i] = which.max(profile_ili[j, , i])
        }
    }

    # Continue plotting if fitted at a higher resolution
    wk.min = round(min(dat_fit_wk_max, sim_fit_wk_max))
    wk.max = round(max(dat_fit_wk_max, sim_fit_wk_max))

    wk.min = round(0.5 * wk.min)
    wk.max = round(1.5 * wk.max)
    wk.max = min(wk.max, nperiods)
    wk.min = max(1, wk.min)

    breaks = seq(from = wk.min, to = wk.max, by = 1)
    breaks_x = seq(from = wk.min, to = wk.max, by = 4)
    labels = weeks[breaks_x]
    ylab = "Probability Density"
    xlab = "EW #"

    plotlist2 = list()
    plot_region = list()
    for (iregion in 1:n.fit) {
        data = data.frame(x = c(dat_fit_wk_max[iregion], sim_fit_wk_max[, iregion]), y = c(rep(0, length(dat_fit_wk_max[iregion])), rep(1,
            length(sim_fit_wk_max[, iregion]))))
        plot_region[[iregion]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..),
            breaks = breaks, fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..),
            breaks = breaks, fill = "forestgreen", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks),
            breaks = breaks_x, labels = labels) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10, color = "gray20",
            face = "italic"), axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

        reg.name = reg.fit.name[iregion]
        rel.pop = round(fit_coef[iregion], digits = 3)
        leg.text = c("Peak Week", "Data", "Direct")
        plot_region[[iregion]] = plot_region[[iregion]] + annotate("text", x = c(rep(-Inf, 3), rep(Inf, 3)), y = rep(Inf, 6), label = c(paste("   ",
            mydata$FY, sep = ""), paste("   ", reg.name, sep = ""), paste("   ", rel.pop, sep = ""), paste(leg.text, "   ", sep = "")),
            hjust = c(0, 0, 0, 1, 1, 1), vjust = c(2.5, 4, 5.5, 2.5, 4, 5.5), col = c("black", "black", "black", "black", "dodgerblue",
                "forestgreen"), family = "serif", size = 3.5)
    }
    plotlist2 = append(plotlist2, plot_region)

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
    breaks_x = seq(from = wk.min, to = wk.max, by = 4)
    labels = weeks[breaks_x]
    plot_nation = list()
    data = data.frame(x = c(dat_model_wk_max, indrct_model_wk_max), y = c(rep(0, length(dat_model_wk_max)), rep(1, length(indrct_model_wk_max))))
    plot_nation[[1]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..), breaks = breaks,
        fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..), breaks = breaks,
        fill = "forestgreen", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks), breaks = breaks_x,
        labels = labels) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10, color = "gray20", face = "italic"),
        axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))

    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    leg.text = c("Peak Week", "Data", "Indirect")
    plot_nation[[1]] = plot_nation[[1]] + annotate("text", x = c(rep(-Inf, 2), rep(Inf, 3)), y = rep(Inf, 5), label = c(paste("   ",
        mydata$FY, sep = ""), paste("   ", model.name, sep = ""), paste(leg.text, "   ", sep = "")), hjust = c(0, 0, 1, 1, 1), vjust = c(2.5,
        4, 2.5, 4, 5.5), col = c("black", "black", "black", "dodgerblue", "forestgreen"), family = "serif", size = 3.5)

    if (!is.null(model_profile_ili)) {
        data = data.frame(x = c(dat_model_wk_max, drct_model_wk_max), y = c(rep(0, length(dat_model_wk_max)), rep(1, length(drct_model_wk_max))))
        plot_nation[[2]] = ggplot(data = data, aes(x = x)) + geom_histogram(data = subset(data, y == 0), aes(y = ..density..), breaks = breaks,
            fill = "dodgerblue", col = "black", alpha = 0.7) + geom_histogram(data = subset(data, y == 1), aes(y = ..density..), breaks = breaks,
            fill = "mediumpurple", col = "black", alpha = 0.7) + scale_x_continuous(name = xlab, limits = range(breaks), breaks = breaks_x,
            labels = labels) + scale_y_continuous(name = ylab) + theme(text = element_text(size = 10, color = "gray20", face = "italic"),
            axis.text.x = element_text(face = "plain"), axis.text.y = element_text(face = "plain"))
    }

    model.name = mydata$model$name
    gsub(model.name, ".", " ", model.name)
    leg.text = c("Peak Week", "Data", "Direct")
    plot_nation[[2]] = plot_nation[[2]] + annotate("text", x = c(rep(-Inf, 2), rep(Inf, 3)), y = rep(Inf, 5), label = c(paste("   ",
        mydata$FY, sep = ""), paste("   ", model.name, sep = ""), paste(leg.text, "   ", sep = "")), hjust = c(0, 0, 1, 1, 1), vjust = c(2.5,
        4, 2.5, 4, 5.5), col = c("black", "black", "black", "dodgerblue", "mediumpurple"), family = "serif", size = 3.5)

    plotlist2 = append(plotlist2, plot_nation)
    s = nrow * ncol
    layout = matrix(seq(1, s, 1), nrow = nrow, byrow = TRUE)
    multiplot(plotlist = plotlist2, layout = layout)


    if (tolower(device) == "pdf" | tolower(device) == "png")
        dev.off()

    return(err = 0)

}


# This function is used by the plotting routines that are based on ggplot
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
	#'
	#' Setup a grid for multi-plots on one page using ggplot2
	#'
	#' @param plotlist A list with a ggplot in each element
	#' @param file    filename for output
	#' @param cols Integer - the number of columns for the grid
	#' layout optional a 2D matrix with the plot number in each location
	#' matrix has ncols and nrow
	#'
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots == 1) {

    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
        }
    }
}




