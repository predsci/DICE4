	##
## Functions Related to Priors and Data Augmentation for the CDC Data: either Posterior of a Past year or Data Augmentation
##

get.cdc.da    <- function(mydata = NULL) {

	#' Augment a CDC Time Series
	#'
	#' \code{get.cdc.da} is a wrapper for the data augmentation procedure
	#' Based on the value of the parameter 'da' it calls the function that augments
	#' with the historic null model (da=1) or the most similar season (da=2).
	#' If da=0 there is no data augmentation
  #' @param mydata A dataframe with all the available for this \pkg{DICE} run
	#' @return - A mydata data frame with the additional augmented data and its weight
	#' @examples
	#' mydata <- get.cdc.da(mydata = mydata)
	#'

	if(is.null(mydata)) return

    if (mydata$da == 1 & tolower(mydata$data_source) == "cdc") { ## Augement with historic average
    	cat("\n Getting Future Data: Historic Average \n")
    	mydata = get.future.average.data(mydata = mydata)
    } else if (mydata$da == 2 & tolower(mydata$data_source) == "cdc") {

    	cat("\n Getting Future Data: Most Similar Season \n") ## Augment with most similar season
    	mydata = get.future.similar.data(mydata = mydata)
    } else { ## No Augmentation
    	cat("\n No Data Augmentation \n")
    	mydata$model$wght = rep(0, mydata$nperiods)
    	mydata$model$wght[1:mydata$nperiodsFit] = 1
    	mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
    	mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
    }

	return(mydata)

}

get.cdc.prior <- function(mydata = NULL) {

	#' Retrieve Prior Parameters for a CDC Run
	#'
	#' \code{get.cdc.prior} is a wrapper that calls the \code{\link{get.prior}} function
  #' @param mydata A dataframe with all the available for this \pkg{DICE} run
	#' @return mydata dataframe with additional information about the prior
  #'
  #' @examples
  #' mydata <- get.cdc.prior(mydata = mydata)
  #'

	if(is.null(mydata)) return
    ## Let's get prior for all cases for CDC mydata-type and current season, but not for state level mydata!!
    if (tolower(mydata$data_source) == "cdc") {
        if (mydata$fit_level > 3 | mydata$mod_level > 3 | mydata$imodel == 5) {
            mydata$fit$prior.start = NA
            mydata$fit$prior.end = NA
            mydata$model$prior.start = NA
            mydata$model$prior.end = NA
        } else {
            cat("\n Getting Prior \n")
            my.prior = get.prior(mydata = mydata)
            mydata$fit$prior.start = my.prior$fit.start
            mydata$fit$prior.end = my.prior$fit.end
            mydata$model$prior.start = my.prior$model.start
            mydata$model$prior.end = my.prior$model.end
        }
    }

	return(mydata)
}

get.prior <- function(mydata = NULL) {

    #' Find the Most Similar Season
    #'
    #' Given the fit and model ILI data find the most similar season for each of the fit and model regions
    #'
    #' This function is needed when running in a forecast mode using a prior
    #'
    #' @param mydata A dataframe with all the available for this \pkg{DICE} run
    #' @return  A list with the most similar season year for each of the fit and model regions
    #' @examples
    #' prior = get.prior(mydata = mydata)

    dataType = mydata$data_source
    nperiodsData = min(mydata$nperiodsData, 52)
    nperiods = mydata$nperiods

    my.year = mydata$years[1]
    my.year1 = mydata$years[nperiods]

    mod_level = mydata$mod_level
    fit_level = mydata$fit_level

    model_raw = mydata$model$raw
    fit_raw = mydata$fit$raw
    nregion = mydata$fit$nregions

    nregion1 = nregion + 1

	## Once get.first.year is fixed uncomment the next three lines and comment the following two lines
    #firstYear = get.first.year(mydata = mydata)
    #year.start <- seq(from = firstYear$model[1], to = (my.year - 1), by = 1)
    #year.end <- seq(from = (firstYear$model[1] + 1), to = my.year, by = 1)

    year.start <- seq(from = 2004, to = (my.year - 1), by = 1)
    year.end <- seq(from = 2005, to = my.year, by = 1)

    nyears = length(year.start)

    corrMat = array(0, c(nyears, nregion1))
    colnames(corrMat) = c(mydata$fit$name, mydata$model$name)
    rownames(corrMat) = year.start
    if (mod_level == 2) {
        name = c(NAME_2 = "USA")

    } else {
        name3 = mydata$model$attr$ABBV_3
        name = c(NAME_2 = "USA", NAME_3 = name3)
    }


    for (iyear in 1:nyears) {
        year = year.start[iyear]
        year1 = year.end[iyear]

        #pastdata = get.cdc.data(mod_level = mod_level, fit_level = fit_level, mod_name = name, start.year = year, end.year = year1, data_source = mydata$data_source, all_years_flag=F, NOAA_clim=F)


 pastdata <- get.DICE.data(data_source = 'cdc', mod_level = mod_level, fit_level = fit_level, year = year, mod_name=mydata$model$name, RegState = 'usa', fit_names='all', db_opts=mydata$db_opts, disease = 'flu',  all_years_flag=F, all_cad_clim=F)
        #pastdata = get.subset(start.year = year, end.year = year1, mod_level = mod_level, fit_level = fit_level, data_source = dataType)
        #pastdata = pastdata$mydata

        if (nregion == 1) {
            corrMat[iyear, nregion] = cor(fit_raw[1:nperiodsData, 1], pastdata$fit$raw[1:nperiodsData, 1], method = "pearson")
        } else {
            for (iregion in 1:nregion) {
                corrMat[iyear, iregion] = cor(fit_raw[1:nperiodsData, iregion], pastdata$fit$raw[1:nperiodsData, iregion], method = "pearson")
            }
        }

        corrMat[iyear, nregion1] = cor(model_raw[1:nperiodsData], pastdata$model$raw[1:nperiodsData], method = "pearson")

    }

    prior.start = rep(0, nregion1)
    prior.end = rep(0, nregion1)
    names(prior.start) = names(prior.end) = c(mydata$fit$name, mydata$model$name)

    for (iregion in 1:nregion1) {
        j = which.max(corrMat[, iregion])
        prior.start[iregion] = year.start[j]
        prior.end[iregion] = year.end[j]
    }

    prior = list(fit.start = prior.start[1:nregion], fit.end = prior.end[1:nregion], model.start = prior.start[nregion1], model.end = prior.end[nregion1])

    return(prior)
}

get.future.similar.data <- function(mydata = NULL) {
    #' Data Augmentation Using the Most Similar Season
    #'
    #' given the fit and model ILI data find the most similar season for each and grab the data from
    #' these years for future data points > nperiodsData
    #'
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata A dataframe with all the available for this \pkg{DICE} run
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' prior = get.future.similar.data(mydata = mydata)
    #'

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsFit
	nperiodsFit     = mydata$nperiodsFit
    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has ILI Data for this Entire Season\n\nData Will NOT be Augmented\n\nInternally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your da Value of 2 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

    cat("\n\nGetting Future Data\n\n")
    my.year = mydata$years[1]
    my.year1 = mydata$years[nperiods]

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

    if (mod_level == 2) {
    	name2 = "USA" #mydata$model$attr$ABBV_2
        name = c(NAME_2 = name2)
 	    if (fit_level == 2 || fit_level == 3) {
 	    	year.start = setdiff(2003:(my.year-1), my.year)
 	    } else {
 	    	year.start = setdiff(2010:(my.year-1), my.year)
 	    }
    	year.end = year.start + 1
    } else if (mod_level == 3){
    	name2 = "USA" # mydata$model$attr$ABBV_2
        name3 = mydata$model$attr$ABBV_3
        name = c(NAME_2 = name2, NAME_3 = name3)
    	year.start = setdiff(2010:(my.year-1), my.year)
    	year.end = year.start + 1
    } else if (mod_level == 4){
   		name2 = "USA" # mydata$model$attr$ABBV_2
        name3 = mydata$model$attr$ABBV_3
        name4 = mydata$model$attr$ABBV_4
        name = c(NAME_2 = name2, NAME_3 = name3, NAME_4 = name4)
    	year.start = setdiff(2010:(my.year-1), my.year)
    	year.end = year.start + 1
    } else {
    	cat("\n\n Can NOT run DA for mod_level > 4\n Code will QUIT !!\n")
    	quit()

    }

    model_raw = mydata$model$raw
    fit_raw = mydata$fit$raw
    nregion = mydata$fit$nregions
    nregion1 = nregion + 1

    nyears = length(year.start)

    corrMat = array(0, c(nyears, nregion1))
    colnames(corrMat) = c(mydata$fit$name, mydata$model$name)
    rownames(corrMat) = year.start

    nperiods52 = 52

    for (iyear in 1:nyears) {
        year = year.start[iyear]
        year1 = year.end[iyear]
        # exclude the H1N1 pandemic year
        if (year == 2009)
            next
		#pastdata = get.cdc.data(mod_level = mod_level, fit_level = fit_level, mod_name = name, start.year = year, end.year = year1, data_source = mydata$data_source, all_years_flag=F, NOAA_clim=F)
        #pastdata = pastdata$mydata
		pastdata = get.DICE.data(data_source = 'cdc', mod_level = mod_level, fit_level = fit_level, year = year,  mod_name=mydata$model$name, RegState = 'usa', fit_names='all', db_opts=mydata$db_opts, disease = 'flu', all_years_flag=F, all_cad_clim=F)


        for (iregion in 1:nregion) {
            corrMat[iyear, iregion] = cor(fit_raw[1:nperiodsDataUse, iregion], pastdata$fit$raw[1:nperiodsDataUse, iregion], method = "pearson", use = 'complete.obs')
        }
        corrMat[iyear, nregion1] = cor(model_raw[1:nperiodsDataUse], pastdata$model$raw[1:nperiodsDataUse], method = "pearson", use = 'complete.obs')

    }

    smlr.start = rep(0, nregion1)
    smlr.end = rep(0, nregion1)
    names(smlr.start) = names(smlr.end) = c(mydata$fit$name, mydata$model$name)
    wght_ftr = rep(0, nregion1)
    # The weight is determined by the correlation
    wght_ftr = rep(0, nregion1)
    for (iregion in 1:nregion1) {
        j = which.max(corrMat[, iregion])
        smlr.start[iregion] = year.start[j]
        smlr.end[iregion] = year.end[j]
        wght_ftr[iregion] = corrMat[j, iregion]
    }
    fit.epi = mydata$fit$epi
    model.epi = mydata$model$epi
    # Now we need to grab the data from these years and put it as the future data

    for (iregion in 1:nregion) {
        year = smlr.start[iregion]
        year1 = smlr.end[iregion]
        #future.data = get.cdc.data(mod_level = mod_level, fit_level = fit_level, mod_name = name, start.year = year, end.year = year1, data_source = mydata$data_source, all_years_flag=F, NOAA_clim=F)

        #future.data = future.data$mydata

        future.data = get.DICE.data(data_source = 'cdc', mod_level = mod_level, fit_level = fit_level, year = year,  mod_name=mydata$model$name, RegState = 'usa', fit_names='all', db_opts=mydata$db_opts, disease = 'flu', all_years_flag=F, all_cad_clim=F)


        future = future.data$fit$epi[nperiodsDataUse, iregion]
        current = mydata$fit$epi[nperiodsDataUse, iregion]
        shft = current - future
        mydata$fit$epi[(nperiodsDataUse + 1):nperiods52, iregion] = future.data$fit$epi[(nperiodsDataUse + 1):nperiods52, iregion] + shft
        if (nperiods > nperiods52)
            mydata$fit$epi[nperiods, iregion] = mydata$fit$epi[nperiods52, iregion]
    }

    # Now for the mod_level
    year = smlr.start[nregion1]
    year1 = smlr.end[nregion1]
    #future.data = get.cdc.data(mod_level = mod_level, fit_level = fit_level, mod_name = name, start.year = year, end.year = year1, data_source = mydata$data_source, all_years_flag=F, NOAA_clim=F)

    #future.data = future.data$mydata
	future.data = get.DICE.data(data_source = 'cdc', mod_level = mod_level, fit_level = fit_level, year = year, mod_name=mydata$model$name, RegState = 'usa', fit_names='all', db_opts=mydata$db_opts, disease = 'flu', all_years_flag=F, all_cad_clim=F)

    future = future.data$model$epi[nperiodsDataUse]
    current = mydata$model$epi[nperiodsDataUse]
    shft = current - future
    mydata$model$epi[(nperiodsDataUse + 1):nperiods52] = future.data$model$epi[(nperiodsDataUse + 1):nperiods52] + shft
    if (nperiods > nperiods52)
        mydata$model$epi[nperiods] = mydata$model$epi[nperiods52]

    # avoid having negative values - this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gama function

    for (iweek in (nperiodsDataUse + 1):nperiods) {
        for (iregion in 1:nregion) {
            mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1)
        }
        mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1)
    }

    wght_data = 1

    cat("\n Giving Future Data Points a weight of", wght_ftr, "\n\n")
    mydata$fit$wght = array(0, c(nperiods, nregion))
    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsDataUse] = wght_data
    mydata$model$wght[(nperiodsDataUse + 1):nperiods] = wght_ftr[nregion1]

    for (iregion in 1:nregion) {
        mydata$fit$wght[1:nperiodsDataUse, iregion] = wght_data
        mydata$fit$wght[(nperiodsDataUse + 1):nperiods, iregion] = wght_ftr[iregion]
    }

    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

    return(mydata)
}

get.future.average.data <- function(mydata = NULL) {
    #' Data Augmentation Using the Historic NULL Model
    #'
    #' Use the average weekly incidence value for each region to augment the data
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata - A dataframe with all the available data for this \pkg{DICE} run
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' prior = get.future.average.data(mydata = mydata)

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsFit
    nperiodsFit     = mydata$nperiodsFit

    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has ILI Data for this Entire Season\n\nData Will NOT be Augmented\n\nInternally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your da Value of 1 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

    my.year = mydata$years[1]
    my.year1 = mydata$years[nperiods]

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

     if (mod_level == 2) {
    	name2 = "USA" #mydata$model$attr$ABBV_2
        name = c(NAME_2 = name2)
 	    if (fit_level == 2 || fit_level == 3) {
 	    	year.start = setdiff(2003:(my.year-1), my.year)
 	    } else {
 	    	year.start = setdiff(2010:(my.year-1), my.year)
 	    }
    	year.end = year.start + 1
    } else if (mod_level == 3){
    	name2 = "USA" # mydata$model$attr$ABBV_2
        name3 = mydata$model$attr$ABBV_3
        name = c(NAME_2 = name2, NAME_3 = name3)
    	year.start = setdiff(2010:(my.year-1), my.year)
    	year.end = year.start + 1
    } else if (mod_level == 4){
   		name2 = "USA" # mydata$model$attr$ABBV_2
        name3 = mydata$model$attr$ABBV_3
        name4 = mydata$model$attr$ABBV_4
        name = c(NAME_2 = name2, NAME_3 = name3, NAME_4 = name4)
    	year.start = setdiff(2010:(my.year-1), my.year)
    	year.end = year.start + 1
    } else {
    	cat("\n\n Can NOT run DA for mod_level > 4\n Code will QUIT !!\n")
    	quit()

    }

    model_raw = mydata$model$raw
    fit_raw = mydata$fit$raw
    nregion = mydata$fit$nregions
    nregion1 = nregion + 1

    nyears = length(year.start)

    hstrc.ave = array(0, c(nperiods, nregion1))
    shft = rep(0, nregion1)
    my.cor = rep(0, nregion1)

    nperiods52 = 52

    jj = 1
    for (iyear in 1:nyears) {
        year = year.start[iyear]
        year1 = year.end[iyear]
        # exclude the H1N1 pandemic year
        if (year == 2009)
            next
        #pastdata = get.cdc.data(mod_level = mod_level, fit_level = fit_level, mod_name = name, start.year = year, end.year = year1, data_source = mydata$data_source, all_years_flag=F, NOAA_clim=F)
        #pastdata = pastdata$mydata

        pastdata = get.DICE.data(data_source = 'cdc', mod_level = mod_level, fit_level = fit_level, year = year, mod_name=mydata$model$name, RegState = 'usa', fit_names='all', db_opts=mydata$db_opts, disease = 'flu', all_years_flag=F, all_cad_clim=F)


        hstrc.ave[1:nperiods52, 1:nregion] = hstrc.ave[1:nperiods52, 1:nregion] + as.matrix(pastdata$fit$raw[1:nperiods52, 1:nregion], ncol = nregion)
        hstrc.ave[1:nperiods52, nregion1] = hstrc.ave[1:nperiods52, nregion1] + pastdata$model$raw[1:nperiods52]

        jj = jj + 1
    }
    if (nperiods > nperiods52)
        hstrc.ave[nperiods, ] = hstrc.ave[nperiods52, ]
    hstrc.ave = hstrc.ave/(jj - 1)

    for (iregion in 1:nregion) {
        shft[iregion] = mydata$fit$raw[nperiodsDataUse, iregion] - hstrc.ave[nperiodsDataUse, iregion]
        my.cor[iregion] = cor(mydata$fit$raw[1:nperiodsDataUse, iregion], hstrc.ave[1:nperiodsDataUse, iregion], method = "pearson", use = 'complete.obs')
    }

    shft[nregion1] = mydata$model$raw[nperiodsDataUse] - hstrc.ave[nperiodsDataUse, nregion1]
    my.cor[nregion1] = cor(mydata$model$raw[1:nperiodsDataUse], hstrc.ave[1:nperiodsDataUse, nregion1], method = "pearson", use = 'complete.obs')


    fit.epi = mydata$fit$epi
    model.epi = mydata$model$epi

    ## We augment the epi data not the raw using the factor to go between raw (% ILI) and epi (number of cases)

    for (iregion in 1:nregion) {
        mydata$fit$epi[(nperiodsDataUse + 1):nperiods, iregion] = round((hstrc.ave[(nperiodsDataUse + 1):nperiods, iregion] + shft[iregion]) *
            as.numeric(mydata$fit$factor[iregion]))

    }

    mydata$model$epi[(nperiodsDataUse + 1):nperiods] = round((hstrc.ave[(nperiodsDataUse + 1):nperiods, nregion1] + shft[nregion1]) * as.numeric(mydata$model$factor))

    # avoid having negative values - this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gamma function
    for (iweek in (nperiodsDataUse + 1):nperiods) {
        for (iregion in 1:nregion) {
            mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1)
        }
        mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1)
    }

    wght_data = 1

    wght_ftr = my.cor
    cat("\n Giving Future Data Points a weight of", wght_ftr, "\n\n")
    mydata$fit$wght = array(0, c(nperiods, nregion))
    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsDataUse] = wght_data
    mydata$model$wght[(nperiodsDataUse + 1):nperiods] = wght_ftr[nregion1]

    for (iregion in 1:nregion) {
        mydata$fit$wght[1:nperiodsDataUse, iregion] = wght_data
        mydata$fit$wght[(nperiodsDataUse + 1):nperiods, iregion] = wght_ftr[iregion]
    }

    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

    return(mydata)
}

get.sd.da <- function(mydata = NULL) {

	#' Augment a San Diego Time Series
	#'
	#' \code{get.sd.da} is a wrapper for the data augmentation procedure
	#' Based on the value of the parameter 'da' it calls the function that augments
	#' with the historic null model (da=1) or the most similar season (da=2).
	#' If da=0 there is no data augmentation
  #' @param mydata A dataframe with all the available for this \pkg{DICE} run
	#' @return - A mydata data frame with the additional augmented data and its weight
	#' @examples
	#' mydata <- get.cdc.da(mydata = mydata)
	#'

	if(is.null(mydata)) return
	## mod_name is built here

    if (mydata$da == 1) { ## Augement with historic average
    	cat("\n Getting Future Data: Historic Average \n")
    	mydataNew = get.sd.future.average.data(mydata = mydata)
    	mydata$model$epirun = mydataNew$model$epi
    	mydata$fit$epirun = mydataNew$fit$epi
    	mydata$model$wght = mydataNew$model$wght
    	mydata$fit$wght = mydataNew$fit$wght
    	mydata$model$gamaepirun = mydataNew$model$gamaepi
    	mydata$fit$gamaepirun = mydataNew$fit$gamaepi
    } else if (mydata$da == 2) {

    	cat("\n Getting Future Data: Most Similar Season \n") ## Augment with most similar season
    	mydataNew = get.sd.future.similar.data(mydata = mydata)
   		mydata$model$epirun = mydataNew$model$epi
    	mydata$fit$epirun = mydataNew$fit$epi
    	mydata$model$wght = mydataNew$model$wght
    	mydata$fit$wght = mydataNew$fit$wght
    	mydata$model$gamaepirun = mydataNew$model$gamaepi
    	mydata$fit$gamaepirun = mydataNew$fit$gamaepi

    } else { ## No Augmentation
    	cat("\n No Data Augmentation \n")
  		mydata$model$epirun = mydata$model$epi
    	mydata$fit$epirun = mydata$fit$epi
    	mydata$model$gamaepirun = mydata$model$gamaepi
    	mydata$fit$gamaepirun = mydata$fit$gamaepi
    	mydata$model$wght = rep(0, mydata$nperiods)
    	mydata$model$wght[1:mydata$nperiodsFit] = 1
    	mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
    	mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
    }

	return(mydata)

}

get.sd.future.similar.data <- function(mydata = NULL) {
    #' Data Augmentation Using the Most Similar Season-San Diego
    #'
    #' given the fit and model ILI data find the most similar season for each and grab the data from
    #' these years for future data points > nperiodsData
    #'
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata A dataframe with all the available for this \pkg{DICE} run
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' prior = get.sd.future.similar.data(mydata = mydata)
    #'
	## To avoid problems due to non-reporting at different weeks for different years
	## this fxn works with mydata$model$epi/mydata$fit$epi and NOT mydata$model$raw/mydata$fit$raw
	## Hard coded for San Diego

	mod_name=c(NAME_2="US", NAME_3="R9", NAME_4="CA", NAME_5="CHD1", NAME_6="SD")

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsDataUse
	nperiodsFit     = mydata$nperiodsFit
    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has ILI Data for this Entire Season\n\nData Will NOT be Augmented\n\nInternally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your Initial Prior Value of 3 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

    cat("\n\nGetting Future Data\n\n")
    my.year = mydata$years[1]
    my.year1 = mydata$years[nperiods]

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level


    model_raw = mydata$model$epi
    fit_raw = mydata$fit$epi
    nregion = mydata$fit$nregions
    nregion1 = nregion + 1

    year.start = setdiff(2008:(my.year-1), my.year)
    year.end = year.start + 1

    nyears = length(year.start)

    corrMat = array(0, c(nyears, nregion1))
    colnames(corrMat) = c(mydata$fit$name, mydata$model$name)
    rownames(corrMat) = year.start

    nperiods52 = 52

    for (iyear in 1:nyears) {
   	year = year.start[iyear]
   	# exclude the H1N1 pandemic year
   	if (year == 2009)
   		next
   	pastdata = get.DICE.data(data_source = mydata$data_source, mod_level = mydata$mod_level,
   		fit_level = mydata$fit_level, mod_name = mod_name, disease = mydata$disease, year = year,
   		db_opts = mydata$db_opts)
   	pastdata = pastdata$mydata
   	for (iregion in 1:nregion) {
   		corrMat[iyear, iregion] = cor(fit_raw[1:nperiodsFit, iregion], pastdata$fit$epi[1:nperiodsFit,
   			iregion], method = "pearson")
   	}
   	corrMat[iyear, nregion1] = cor(model_raw[1:nperiodsFit], pastdata$model$epi[1:nperiodsFit],
   		method = "pearson")

   }

    smlr.start = rep(0, nregion1)
    smlr.end = rep(0, nregion1)
    names(smlr.start) = names(smlr.end) = c(mydata$fit$name, mydata$model$name)
    wght_ftr = rep(0, nregion1)
    # The weight is determined by the correlation
    wght_ftr = rep(0, nregion1)
    for (iregion in 1:nregion1) {
        j = which.max(corrMat[, iregion])
        smlr.start[iregion] = year.start[j]
        smlr.end[iregion] = year.end[j]
        wght_ftr[iregion] = corrMat[j, iregion]
    }
    fit.epi = mydata$fit$epi
    model.epi = mydata$model$epi
    # Now we need to grab the data from these years and put it as the future data

    for (iregion in 1:nregion) {

        year = smlr.start[iregion]

        future.data = get.DICE.data(data_source=mydata$data_source, mod_level=mydata$mod_level, fit_level=mydata$fit_level, mod_name=mod_name, disease=mydata$disease, year=year, db_opts=mydata$db_opts)
        future.data = future.data$mydata
        future = future.data$fit$epi[nperiodsFit, iregion]
        current = mydata$fit$epi[nperiodsFit, iregion]
        shft = current - future
        mydata$fit$epi[(nperiodsFit + 1):nperiods52, iregion] = future.data$fit$epi[(nperiodsFit + 1):nperiods52, iregion] + shft
        if (nperiods > nperiods52)
            mydata$fit$epi[nperiods, iregion] = mydata$fit$epi[nperiods52, iregion]
    }

    # Now for the mod_level
    year  = smlr.start[nregion1]
    year1 = smlr.end[nregion1]

    future.data = get.DICE.data(data_source=mydata$data_source, mod_level=mydata$mod_level, fit_level=mydata$fit_level, mod_name=mod_name, disease=mydata$disease, year=year, db_opts=mydata$db_opts)

	future.data = future.data$mydata

    future = future.data$model$epi[nperiodsFit]
    current = mydata$model$epi[nperiodsFit]
    shft = current - future

    mydata$model$epi[(nperiodsFit + 1):nperiods52] = future.data$model$epi[(nperiodsFit + 1):nperiods52] + shft
    if (nperiods > nperiods52)
        mydata$model$epi[nperiods] = mydata$model$epi[nperiods52]

    # avoid having negative values - this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gama function

    for (iweek in (nperiodsFit + 1):nperiods) {
        for (iregion in 1:nregion) {
            mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1, na.rm = TRUE)
        }
        mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1, na.rm = TRUE)
    }

    wght_data = 1

    cat("\n Giving Future Data Points a weight of", wght_ftr, "\n\n")
    mydata$fit$wght = array(0, c(nperiods, nregion))
    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsFit] = wght_data
    mydata$model$wght[(nperiodsFit + 1):nperiods] = wght_ftr[nregion1]

    for (iregion in 1:nregion) {
        mydata$fit$wght[1:nperiodsFit, iregion] = wght_data
        mydata$fit$wght[(nperiodsFit + 1):nperiods, iregion] = wght_ftr[iregion]
    }

    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

    return(mydata)
}

get.sd.future.average.data <- function(mydata = NULL) {
    #' Data Augmentation Using the Historic NULL Model
    #'
    #' Use the average weekly incidence value for each region to augment the data
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata - A dataframe with all the available data for this \pkg{DICE} run
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' future = get.sd.future.average.data(mydata = mydata)

	## To avoid problems due to non-reporting at different weeks for different years
	## this fxn works with mydata$model$epi/mydata$fit$epi and NOT mydata$model$raw/mydata$fit$raw
	## Hard coded for San Diego

	mod_name=c(NAME_2="US", NAME_3="R9", NAME_4="CA", NAME_5="CHD1", NAME_6="SD")

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsDataUse
    nperiodsFit     = mydata$nperiodsFit

    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has ILI Data for this Entire Season\n\nData Will NOT be Augmented\n\nInternally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your Initial Prior Value of 3 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

    my.year = mydata$years[1]
    my.year1 = mydata$years[nperiods]

	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

    model_raw = mydata$model$epi
    fit_raw = mydata$fit$epi
    nregion = mydata$fit$nregions
    nregion1 = nregion + 1

    year.start <- seq(from = 2008, to = (my.year - 1), by = 1)
    year.end <- seq(from = 2009, to = my.year, by = 1)

    nyears = length(year.start)
    hstrc.ave = array(0, c(nperiods, nregion1))
    shft = rep(0, nregion1)
    my.cor = rep(0, nregion1)

    nperiods52 = 52

    jj = 1
    for (iyear in 1:nyears) {
        year = year.start[iyear]
        year1 = year.end[iyear]
        # exclude the H1N1 pandemic year
        if (year == 2009)
            next
        pastdata = get.DICE.data(data_source=mydata$data_source, mod_level=mydata$mod_level, fit_level=mydata$fit_level, mod_name=mod_name, disease=mydata$disease, year=year, db_opts=mydata$db_opts)
        pastdata = pastdata$mydata
        hstrc.ave[1:nperiods52,1:nregion] = hstrc.ave[1:nperiods52, 1:nregion] + pastdata$fit$epi[1:nperiods52,1:nregion]

        hstrc.ave[1:nperiods52, nregion1] = hstrc.ave[1:nperiods52, nregion1] + pastdata$model$epi[1:nperiods52]
        jj = jj + 1
    }
    if (nperiods > nperiods52)
   	hstrc.ave[nperiods, ] = hstrc.ave[nperiods52, ]
   hstrc.ave = hstrc.ave/(jj - 1)

   for (iregion in 1:nregion) {
   	shft[iregion] = mydata$fit$epi[nperiodsFit, iregion] - hstrc.ave[nperiodsFit,
   		iregion]
   	my.cor[iregion] = cor(mydata$fit$raw[1:nperiodsFit, iregion], hstrc.ave[1:nperiodsFit,
   		iregion], method = "pearson")
   }

    shft[nregion1] = mydata$model$epi[nperiodsFit] - hstrc.ave[nperiodsFit, nregion1]
    my.cor[nregion1] = cor(mydata$model$epi[1:nperiodsFit], hstrc.ave[1:nperiodsFit, nregion1], method = "pearson")

    fit.epi = mydata$fit$epi
    model.epi = mydata$model$epi

    ## This routine works with epi now raw to avoid the missing NA data

    for (iregion in 1:nregion) {
        mydata$fit$epi[(nperiodsFit + 1):nperiods, iregion] = round(hstrc.ave[(nperiodsFit + 1):nperiods, iregion] + shft[iregion])

    }

    mydata$model$epi[(nperiodsFit + 1):nperiods] = round(hstrc.ave[(nperiodsFit + 1):nperiods, nregion1] + shft[nregion1])

    # avoid having negative values - this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gamma function
    for (iweek in (nperiodsFit + 1):nperiods) {
        for (iregion in 1:nregion) {
            mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1, na.rm = TRUE)
        }
        mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1, na.rm = TRUE)
    }

    wght_data = 1

    wght_ftr = my.cor
    cat("\n Giving Future Data Points a weight of", wght_ftr, "\n\n")
    mydata$fit$wght = array(0, c(nperiods, nregion))
    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsFit] = wght_data
    mydata$model$wght[(nperiodsFit + 1):nperiods] = wght_ftr[nregion1]

    for (iregion in 1:nregion) {
        mydata$fit$wght[1:nperiodsFit, iregion] = wght_data
        mydata$fit$wght[(nperiodsFit + 1):nperiods, iregion] = wght_ftr[iregion]
    }

    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

    return(mydata)
}


setup.single.prior <- function(mydata = NULL, logvec = NULL) {

    #' Obtain Prior parameters for an Uncoupled CDC \pkg{DICE} Run
    #'
    #' Use the pre-computed data base of mean and sd values for the parameters
    #' to obtain them for each region
    #' @param mydata - A dataframe with all the information for this \pkg{DICE} run
    #' @param logvec - a vector of zero and 1 need to be updated to value of 2 for prior paramters
    #' @return
    #' A list with the mean and sd values for the fit regions followed by the model region
    #' The number of columns depends on the model
    #' an updated logvec
    #' model 1 (SH): pC, R0, SH
    #' model 2 (SV): pC, R0, SV
    #' model 3 (SH+SV): pC, R0, SH, SV
    #' model 4: pC, R0
    #' model 5: pC, R0
    #' example:
    #' setup.prior = setup.single.prior(mydata=mydata,logvec=logvec)

    if (mydata$mod_level > 3 | mydata$fit_level > 3 | mydata$imodel == 5 | mydata$data_source != 'cdc') {

        ymu = array(0, c((mydata$fit$nregions + 1), 4))

        sigma = array(1, c((mydata$fit$nregions + 1), 4))

        setup.prior = list(ymu = ymu, sigma = sigma, logvec = logvec)

         return(setup.prior)
    } else {



        n = mydata$fit$nregions

        # n1 = 1

        # if (n > 1)
            n1 = n + 1

        ## order of parameters as in par_names

        prior.list = c("R0", "pC") # c("pC", "R0")
        logvec['pC'] = 2
        logvec['R0'] = 2

        if (mydata$imodel == 1) {
            ## deltaR, aparam
            prior.list = c(prior.list, "SH")
            logvec['deltaR'] = 2
            logvec['aparam'] = 2
        }
        if (mydata$imodel == 2) {
            ## alpha
            prior.list = c(prior.list, "SV")
            logvec['alpha'] = 2
        }
        if (mydata$imodel == 3) {
            ## deltaR, aparam, alpha
            prior.list = c(prior.list, "SH", "SV")
            logvec['deltaR'] = 2
            logvec['aparam'] = 2
            logvec['alpha'] = 2
        }

        nprior = length(prior.list)
        ymu = array(0, dim = c(n1, nprior))
        sigma = array(1, dim = c(n1, nprior))

        prior.list_mean = paste(prior.list, "-Mean", sep = "")
        prior.list_sd = paste(prior.list, "-SD", sep = "")

        colnames(ymu) = prior.list
        colnames(sigma) = prior.list

        rownames(ymu) = c(mydata$fit$name, mydata$model$name)
        rownames(sigma) = c(mydata$fit$name, mydata$model$name)

        # subset the prior Table to the model we are using
        myPrior = priorTable[as.numeric(priorTable[, "Model"]) == mydata$imodel, ]

        for (i in 1:n) {
            if (!is.na(mydata$fit$prior.start[i])) {
                start = mydata$fit$prior.start[i]
                end = mydata$fit$prior.end[i]
                start = as.numeric(start)
                end = as.numeric(end)
                myName = mydata$fit$name[i]
                #myName = gsub(" ","",myName)
                myName = gsub(" ",".",myName)
                cat(i, start, end, myName, '\n')
                index = which(as.numeric(myPrior[, "start"]) == start & myPrior[, "Region"] == myName)
               # cat(index,'\n')
               # cat(myPrior[index, prior.list_mean],'\n')
               # cat(myPrior[index, prior.list_sd],'\n')
                ymu[i, 1:nprior] = myPrior[index, prior.list_mean]
                sigma[i, 1:nprior] = myPrior[index, prior.list_sd]
            } else {
                ymu[i, 1:nprior] = 0
                sigma[i, 1:nprior] = 1
            }

        }
        # Sanity check
        for (j in 1:nprior) {
        	for (i in 1:n) {
        		if (sigma[i, j] <= 0) {
        			tmp = sigma[1:n,j]
        			tmp = tmp[-c(i)]
        			sigma[i, j] = mean(as.numeric(tmp))
        		}
        	}
        }

        if ('SV' %in% prior.list) {
        	for(i in 1:n) {
        		if (ymu[i,'SV'] <= 0) {
        			tmp = ymu[1:n,'SV']
        			tmp = tmp[-c(i)]
        			ymu[i,'SV'] = mean(as.numeric(tmp))
        		}
        	}
        }
        # repeat for the model
        if (!is.na(mydata$model$prior.start)) {
            start = mydata$model$prior.start
            end = mydata$model$prior.end
            start = as.numeric(start)
            end = as.numeric(end)
            myName = mydata$model$name
            myName = gsub(" ",".",myName)
            	cat(i, start, end, myName, '\n')
            index = which(as.numeric(myPrior[, "start"]) == start & myPrior[, "Region"] == myName)
            ymu[n1, 1:nprior] = myPrior[index, prior.list_mean]
            sigma[n1, 1:nprior] = myPrior[index, prior.list_sd]
        } else {
            ymu[n1, 1:nprior] = 0
            sigma[n1, 1:nprior] = 1
        }

        setup.prior = list(ymu = ymu, sigma = sigma, logvec = logvec)

    }


    return(setup.prior)


}

sample.from.prior <- function(mydata = NULL, model_par = NULL, fit_par = NULL, model_pmin = NULL, model_pmax = NULL, fit_pmin = NULL, fit_pmax = NULL, model_ymu = NULL, model_sigma = NULL, fit_ymu = NULL, fit_sigma = NULL) {

	#' Sample from a Prior Gaussian Distribution
	#'
	#' \code{sample.from.prior} samples initial guesses for the MCMC procedure using
	#' a prior Gaussian distribution for: pC, R0 and the parameters that control
	#' the force of infection for the SH and/or SV models
	#' @param mydata A dataframe with all the information for this \pkg{DICE} run
	#' @param model_par An array with initial guesses for the parameters for the model region
	#' @param fit_par   A matrix with initial guesses for the parameters for the fit regions
	#' @param model_pmin An array with minimum values for all the model region parameters
	#' @param model_pmax An array with maximum values for all the model region parameters
	#' @param fit_pmin A matrix with minimum values for all the fit regions parameters
	#' @param fit_pmax A matrix with maximum values for all the fit regions parameters
	#' @param model_ymu An array with the prior mean values for the model region
	#' @param model_sigma An array with the prior standard deviation values for the model region
	#' @param fit_ymu A matrix with the prior mean values for the fit regions
	#' @param fit_sigma An array with the prior standard deviation values for the model region
	#' @return An array and matrix of initial parameters for the MCMC procedure
	#' @examples
	#' sample.from.prior(mydata = mydata, model_par = model_par, fit_par = fit_par,
	#' model_pmin = model_pmin, model_pmax = model_pmax, fit_pmin = fit_pmin, fit_pmax = fit_pmax,
	#' model_ymu = model_ymu, model_sigma = model_sigma, fit_ymu = fit_ymu, fit_sigma = fit_sigma)
	#'

	n = dim(fit_par)[2]

	## First and second places are pC and R0 - which are 2nd and 3rd in parameter list

	for (i in 1:n) {
		success = FALSE
		while (!success) {
			xmean = fit_ymu[i, "pC"]
			xsd = fit_sigma[i, "pC"]
			fit_par["pC", i] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
			if (fit_par["pC", i] > fit_pmin["pC", i] & fit_par["pC", i] < fit_pmax["pC", i])
				success = TRUE
		}

		success = FALSE
		while (!success) {
			xmean = fit_ymu[i, "R0"]
			xsd = fit_sigma[i, "R0"]
			fit_par["R0", i] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
			if (fit_par["R0", i] > fit_pmin["R0", i] & fit_par["R0", i] < fit_pmax["R0", i])
				success = TRUE
		}

	}

	success = FALSE
	while (!success) {
		xmean = model_ymu["pC"]
		xsd = model_sigma["pC"]
		model_par["pC"] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
		if (model_par["pC"] > model_pmin["pC"] & model_par["pC"] < model_pmax["pC"])
			success = TRUE

	}

	success = FALSE
	while (!success) {
		xmean = model_ymu["R0"]
		xsd = model_sigma["R0"]
		model_par["R0"] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
		if (model_par["R0"] > model_pmin["R0"] & model_par["R0"] < model_pmax["R0"])
			success = TRUE
	}

	#' model 1 (SH): pC, R0, SH
	#' model 2 (SV): pC, R0, SV
#' model 3 (SH+SV): pC, R0, SH, SV
#' model 4: pC, R0
#' model 5: pC, R0

	## SV term
	if (mydata$imodel == 2)
		ind = 3
	if (mydata$imodel == 3)
		ind = 4

	if (mydata$imodel == 2 || mydata$imodel == 3 & mydata$prior > 0) {
		for (i in 1:n) {
			success = FALSE
			while (!success) {
				xmean = fit_ymu[i, "SV"]
				xsd = fit_sigma[i, "SV"]
				fit_par["alpha", i] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
				if (fit_par["alpha", i] > fit_pmin["alpha", i] & fit_par["alpha", i] < fit_pmax["alpha", i])
					success = TRUE
			}
		}


		xmean = model_ymu["SV"]
		xsd = model_sigma["SV"]
		cat("xmean", xmean, "xsd", xsd, "\n")
		success = FALSE
		while (!success) {
			model_par["alpha"] = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
			if (model_par["alpha"] > model_pmin["alpha"] & model_par["alpha"] < model_pmax["alpha"])
				success = TRUE
		}

	}

	## SH term is more difficult
	if (mydata$imodel == 1 || mydata$imodel == 3) {

		for (i in 1:n) {
			success = FALSE
			xmean = fit_ymu[i, "SH"]
			xsd = fit_sigma[i, "SH"]
			sh.aveg = mean(mydata$fit$sh[, i])
			while (!success) {
				deltaR = runif(1, min = fit_pmin["deltaR", i], max = fit_pmax["deltaR", i])
				aparam = runif(1, min = fit_pmin["aparam", i], max = fit_pmax["aparam", i])
				myProb = deltaR * exp(-aparam * sh.aveg)
				myPrior = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
				if (myProb >= myPrior) {
					fit_par["deltaR", i] = deltaR
					fit_par["aparam", i] = aparam
					success = TRUE
				}
			}
		}

		success = FALSE
		xmean = model_ymu["SH"]
		xsd = model_sigma["SH"]
		sh.aveg = mean(mydata$model$sh)
		while (!success) {
			deltaR = runif(1, min = model_pmin["deltaR"], max = model_pmax["deltaR"])
			aparam = runif(1, min = model_pmin["aparam"], max = model_pmax["aparam"])
			myProb = deltaR * exp(-aparam * sh.aveg)
			myPrior = rnorm(1, mean = as.numeric(xmean), sd = as.numeric(xsd))
			if (myProb >= myPrior) {
				model_par["deltaR"] = deltaR
				model_par["aparam"] = aparam
				success = TRUE
			}
		}

	}


	prior.ini = list(model_par=model_par, fit_par=fit_par)

	return(prior.ini)
}

get.first.year <- function(mydata = NULL) {
    #' Find the First Year for the Data
    #'
    #' This function is called after mydata was obtained to determine the first year of data for this data type
    #' @param mydata - A dataframe with all the information for this \pkg{DICE} run
    #' @return
    #' StartYear - a list with first year of data for this data type: both model level and fit level
    #' @examples
    #' FirstYear = get.first.year(mydata)

    StartYear = list()
    modelName = names(mydata$model$raw_units)
    fitNames = names(mydata$fit$raw_units)
    # determine first year for model region/dataType
    if (tolower(mydata$data_source) == "cdc") {
        dataName = "CDCili"
    } else if (tolower(mydata$data_source) == "gft") {
        dataName = "GFTili"
    } else if (tolower(mydata$data_source) == "miscili") {
        dataName = "MiscIli"
    } else if (tolower(mydata$data_source) == "misccases") {
        dataName = "MiscCases"
    }
    iliIndex = which(!is.na(diceData[[dataName]][, modelName]))
    StartIndex = min(iliIndex)
    StartYear[["model"]] = diceData[[dataName]][StartIndex, "year"]
    # determine first year for fit region(s)/dataType(s)
    StartYear[["fit"]] = integer(length = length(fitNames))
    for (ii in 1:length(fitNames)) {
        iliIndex = which(!is.na(diceData[[dataName]][, fitNames[ii]]))
        StartIndex = min(iliIndex)
        StartYear[["fit"]][ii] = diceData[[dataName]][StartIndex, "year"]
    }
    return(StartYear)
}


calc.epi.null <- function(mydata = NULL, all_years_epi = NULL, mod_id = NULL, fit_id = NULL) {
    #' Calculate the Historic NULL Model
    #'
    #' \code{calc.null} Creates an epidemic NULL model using average monthly/weekly values from the previous ten years
    #' @param mydata - A dataframe with all the information for this \pkg{DICE} run
    #' @param all_years_data the epi data for all years
    #' @param model_id the abbreviated name for model level region
    #' @param fit_id the abbreviated name for fit level regions
    #' @examples
    #' calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
    #' @return
    #' A list with the NULL model values for the model and fit level data
    #'

    cadence = mydata$cadence
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    if (cadence == 'Monthly') {
    	my_row = which(all_years_epi$years == mydata$years[nperiodsFit] & all_years_epi$months == mydata$months[nperiodsFit])
    	cadence.vec = mydata$months
    }

    if (cadence == 'Weekly') {
    	my_row = which(all_years_epi$years == mydata$years[nperiodsFit] & all_years_epi$weeks == mydata$weeks[nperiodsFit])
    	cadence.vec = mydata$weeks
    }

    if (cadence == 'Daily') {
    	my_row = which(all_years_epi$years == mydata$years[nperiodsFit] & all_years_epi$days == mydata$days[nperiodsFit])
    	cadence.vec = mydata$days
    }

    year.vec = mydata$years

    nregions = mydata$fit$nregions

    attr = c("year", cadence)

    cases.ave.mod  = array(0, c(nperiods, (length(attr)+1)))

    colnames(cases.ave.mod) = c(attr, mod_id)
    cases.ave.mod[, cadence] = cadence.vec

 	cases.ave.fit  = array(0, c(nperiods, (length(attr)+nregions)))

    colnames(cases.ave.fit)  = c(attr, fit_id)  # just allocates the right number of month for the mydata frame of the null model
    cases.ave.fit[, cadence] = cadence.vec

    istart = my_row - nperiods * 10
    istart = max(istart, 1)
    past.cases.mod = all_years_epi$model$raw[istart:(my_row)]

    if (nregions > 1) {
    	past.cases.fit = all_years_epi$fit$raw[istart:(my_row), ]
    } else {
    	past.cases.fit = all_years_epi$fit$raw[istart:(my_row)]
    }



    for (j in 1:nperiods) {
    	if (mydata$cadence == "Monthly")
    		ind = which(all_years_epi$months[istart:my_row] == mydata$months[j])
    	if (mydata$cadence == "Weekly")
    		ind = which(all_years_epi$weeks[istart:my_row] == mydata$weeks[j])
   		if (mydata$cadence == "Daily")
    		ind = which(all_years_epi$days[istart:my_row] == mydata$days[j])

    	cases.ave.mod[j, mod_id] = mean(all_years_epi$model$raw[ind], na.rm = TRUE)

    }
	cases.ave.mod[is.nan(cases.ave.mod[,mod_id]),mod_id] <- NA

    for (i in 1:nregions) {

    	for (j in 1:nperiods) {
    		if (mydata$cadence == "Monthly")
    			ind = which(all_years_epi$months[istart:my_row] == mydata$months[j])
    		if (mydata$cadence == "Weekly")
    			ind = which(all_years_epi$weeks[istart:my_row] == mydata$weeks[j])
   			if (mydata$cadence == "Daily")
    			ind = which(all_years_epi$days[istart:my_row] == mydata$days[j])
    		if (nregions > 1) {
    			cases.ave.fit[j, fit_id[i]] = mean(all_years_epi$fit$raw[ind, i], na.rm = TRUE)
    		} else {
    			cases.ave.fit[j, fit_id[i]] = mean(all_years_epi$fit$raw[ind], na.rm = TRUE)
    		}
    	}
    	cases.ave.fit[is.nan(cases.ave.fit[,fit_id[i]]),fit_id[i]] <- NA
    }


    cases.ave = list(cases.ave.mod = cases.ave.mod, cases.ave.fit = cases.ave.fit)

    return(cases.ave)

}

get.flu_poc.da <- function(mydata = NULL, all_years_epi = NULL) {

	#' Augment an influenza point-of-care data Time Series
	#'
	#' \code{get.flu_poc.da} is a wrapper for the data augmentation procedure
	#' Based on the value of the parameter 'da' it calls the function that augments
	#' with the historic null model (da=1) or the most similar season (da=2).
	#' If da=0 there is no data augmentation
  #' @param mydata A dataframe with all the available for this \pkg{DICE} run
	#' @return - A mydata data frame with the additional augmented data and its weight
	#' @examples
	#' mydata <- get.flu_poc.da(mydata = mydata, all_years_epi=all_years_epi)
	#'

	if(is.null(mydata)) return

    if (mydata$da == 1) { ## Augement with historic average
    	cat("\n Getting Future Data: Historic Average \n")
    	mydata = get.flu_poc.future.average.data(mydata = mydata, all_years_epi = all_years_epi)
    } else if (mydata$da == 2) {

    	cat("\n Getting Future Data: Most Similar Season \n") ## Augment with most similar season
    	mydata = get.flu_poc.future.similar.data(mydata = mydata, all_years_epi = all_years_epi)
    } else { ## No Augmentation
    	cat("\n No Data Augmentation \n")
    	mydata$model$wght = rep(0, mydata$nperiods)
    	mydata$model$wght[1:mydata$nperiodsFit] = 1
    	mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
    	mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
    }

	return(mydata)

}

get.flu_poc.future.similar.data <- function(mydata = NULL, all_years_epi = NULL) {
    #' Data Augmentation Using the Most Similar Season
    #'
    #' given the fit and model flu point-of-care data find the most similar season for each and grab the data from
    #' these years for future data points > nperiodsData
    #'
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata A dataframe with all the available for this \pkg{DICE} run
    #' @param all_years_epi A dataframe with the entire history of disease incidence
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' prior = get.flu_poc.future.similar.data(mydata = mydata)
    #'

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsDataUse
	nperiodsFit     = mydata$nperiodsFit

    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has Point-of-Care Incidence  Data for this Entire Season\n\nData Will NOT be Augmented\n\n
        Internally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your Initial Prior Value of 3 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

    cat("\n\nGetting Future Data Using Most Similar Season\n\n")
    my.year = mydata$years[1]


	corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$model$epi, nfit = mydata$nperiodsFit)
    # This calculated a year but we only need the row or correlation of the current season will all other seasons
    index = which(rownames(corrMAT) == my.year)
    my.cor = corrMAT[index, ]
    index = which.max(my.cor)
    smlr.start = names(my.cor)[index]
    smlr.start = as.numeric(smlr.start)
    smlr.end   = smlr.start + 1
    model_wght_ftr = as.numeric(my.cor[index])

    wght_data = 1.0

    ## Assume Point-of-Care data is weekly - this will subset it properly

    my_row = which(all_years_epi$years == smlr.start & all_years_epi$weeks == mydata$weeks[1])
		future.data = all_years_epi$model$epi[my_row:length(all_years_epi$model$epi)]
    future = future.data[nperiodsFit]
    current = mydata$model$epi[nperiodsFit]
    shft = current - future

    mydata$model$epi[(nperiodsFit + 1):nperiods] = future.data[(nperiodsFit + 1):nperiods] + shft

    # avoid having negative values -
    # this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gama function

    for (iweek in (nperiodsFit + 1):nperiods) {
    	mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1)
    	mydata$model$epi[iweek] = round(mydata$model$epi[iweek])
    }


    cat("\n Giving Future Data Points a weight of", model_wght_ftr, "\n\n")

    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsFit] = wght_data
    mydata$model$wght[(nperiodsFit + 1):nperiods] = model_wght_ftr

    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)


    ## repeat for fit level- if needed
    nregions = mydata$fit$nregions

    mydata$fit$wght = array(data = 0, c(nperiods, nregions))
	print(all_years_epi$fit$epi)
    if (nregions == 1) {
    	mydata$fit$epi[, 1] = mydata$model$epi
    	mydata$fit$wght[, 1] = mydata$model$wght
    	mydata$fit$gamaepi[, 1] = mydata$model$gamaepi
    } else {
    	for (iregion in 1:nregions) {
    		corrMAT = calc.sqldata.cor(mydata = mydata, all_years_epi = all_years_epi, cases = all_years_epi$fit$epi[, iregion],
    			nfit = mydata$nperiodsFit)
    		# This calculated a year but we only need the row or correlation of the current season will all other seasons
    		index = which(rownames(corrMAT) == my.year)
    		my.cor = corrMAT[index, ]
    		index = which.max(my.cor)
    		smlr.start = names(my.cor)[index]
    		smlr.start = as.numeric(smlr.start)
    		smlr.end = smlr.start + 1
    		wght_ftr = as.numeric(my.cor[index])

    		## Assume Point-of-Care data is weekly - this will subset it properly

    		my_row = which(all_years_epi$years == smlr.start & all_years_epi$weeks == mydata$weeks[1])
    		future.data = all_years_epi$fit$epi[my_row:length(all_years_epi$model$epi), iregion]
    		future = future.data[nperiodsFit]
    		current = mydata$fit$epi[nperiodsFit, iregion]
    		shft = current - future

    		mydata$fit$epi[(nperiodsFit + 1):nperiods, iregion] = future.data[(nperiodsFit + 1):nperiods] + shft

    		# avoid having negative values -
    		# this does not affect the fit because it happens at the end of the season, so is really done just for
# the gama function

    		for (iweek in (nperiodsFit + 1):nperiods) {
    			mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1)
    			mydata$fit$epi[iweek, iregion] = round(mydata$fit$epi[iweek, iregion])
    		}

    		cat("\n Giving Future Data Points a weight of", wght_ftr, "\n\n")

    		mydata$fit$wght[1:nperiodsFit, iregion] = wght_data
    		mydata$fit$wght[(nperiodsFit + 1):nperiods, iregion] = wght_ftr
    		mydata$fit$gamaepi[, iregion] = lgamma(mydata$fit$epi[, iregion] + 1)

    	}

    }

    return(mydata)
}

get.flu_poc.future.average.data <- function(mydata = NULL, all_years_epi = NULL) {
    #' Data Augmentation Using the Most Similar Season
    #'
    #' given the model Point-of-Care data, use the historic average to augment the data
    #' This function is needed when running in a forecast mode using data augmentation
    #'
    #' @param mydata A dataframe with all the available for this \pkg{DICE} run
    #' @param all_years_epi A dataframe with the entire history of disease incidence
    #' @return An updated mydata structure that has the future data and has weights for each week of data
    #' @examples
    #' prior = get.flu_poc.future.average.data(mydata = mydata, all_years_epi = all_years_epi)
    #'

    dataType = mydata$data_source
    nperiodsData = mydata$nperiodsData
    nperiods = mydata$nperiods
    nperiodsDataUse = mydata$nperiodsDataUse
	nperiodsFit     = mydata$nperiodsFit

    ## Sanity check - if there is data for the entire ILI season there is no reason to 'augment' the data Let the user know that we are
    ## using the available data for this year Weights are set to 1.0 as is the case when prior != 3

    if (nperiodsFit == nperiods) {
        cat("\n-------------- WARNING --------------\n")
        cat("\nDICE Has Point-of-Care Incidence  Data for this Entire Season\n\nData Will NOT be Augmented\n\n
        Internally Resetting Prior in Code to Zero\n")
        cat("\nDirectory and File Names will Still Show Your Initial Prior Value of 3 !!\n")
        cat("\n--------- END OF WARNING -----------\n\n")
        mydata$model$wght = rep(0, mydata$nperiods)
        mydata$model$wght[1:mydata$nperiodsFit] = 1
        mydata$fit$wght = array(0, c(mydata$nperiods, mydata$fit$nregions))
        mydata$fit$wght[1:mydata$nperiodsFit, 1:mydata$fit$nregions] = 1
        mydata$prior = 0
        return(mydata)
    }

 	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

    mod_id = mydata$model$attr[[paste0("ABBV_", mydata$mod_level)]]
    fit_id = mydata$fit$attr[[paste0("ABBV_", mydata$fit_level)]]

    ## Historic NULL model

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)
	epi.null.mod = epi.null$cases.ave.mod[,mod_id]
	epi.null.fit = epi.null$cases.ave.fit[,fit_id]


	my.cor = cor(epi.null.mod[1:nperiodsFit], mydata$model$epi[1:nperiodsFit])
    model_wght_ftr = my.cor

    nregions = mydata$fit$nregions
    fit_wght_ftr = rep(0, nregions)

    if (nregions > 1) {
    	for (iregion in 1:nregions) {
    	  non_na_ind = !is.na(epi.null.fit[1:nperiodsFit, iregion]) & !is.na(mydata$fit$epi[1:nperiodsFit, iregion])
    		fit_wght_ftr[iregion] = cor(epi.null.fit[1:nperiodsFit, iregion][non_na_ind], mydata$fit$epi[1:nperiodsFit, iregion][non_na_ind])
    	}
    } else {
      non_na_ind = !is.na(epi.null.fit[1:nperiodsFit]) & !is.na(mydata$fit$epi[1:nperiodsFit,1])
    	fit_wght_ftr[1] = cor(epi.null.fit[1:nperiodsFit][non_na_ind], mydata$fit$epi[1:nperiodsFit,1][non_na_ind])
    }


    ## Point-of-Care data is weekly - this will subset it properly

	future.data = epi.null.mod
    future = future.data[nperiodsFit]
    current = mydata$model$epi[nperiodsFit]
    shft = current - future
    mydata$model$epi[(nperiodsFit + 1):nperiods] = future.data[(nperiodsFit + 1):nperiods] + shft

	wght_data = 1.0

	## repeat for all fit regions

	mydata$fit$wght = array(data=0, c(nperiods, nregions))

	if (nregions > 1) {
		for (iregion in 1:nregions) {
			future.data = epi.null.fit[, iregion]
			future = future.data[nperiodsFit]
			current = mydata$fit$epi[nperiodsFit, iregion]
			shft = round(current - future)
			mydata$fit$epi[(nperiodsFit + 1):nperiods, iregion] = future.data[(nperiodsFit + 1):nperiods] + shft
			mydata$fit$wght[1:nperiodsFit, iregion] = wght_data
			mydata$fit$wght[(nperiodsFit + 1):nperiods, iregion] = fit_wght_ftr[iregion]
		}
	} else {
		future.data = epi.null.fit
		future = future.data[nperiodsFit]
		current = mydata$fit$epi[nperiodsFit,1]
		shft = round(current - future)
		mydata$fit$epi[(nperiodsFit + 1):nperiods, 1] = future.data[(nperiodsFit + 1):nperiods] + shft
		mydata$fit$wght[1:nperiodsFit, 1] = wght_data
		mydata$fit$wght[(nperiodsFit + 1):nperiods, 1] = fit_wght_ftr[1]

	}

    # avoid having negative values -
    # this does not affect the fit because it happens at the end of the season, so is really done just for
    # the gama function

    for (iweek in (nperiodsFit + 1):nperiods) {
        mydata$model$epi[iweek] = max(mydata$model$epi[iweek], 1)
        mydata$model$epi[iweek] = round(mydata$model$epi[iweek])
        for (iregion in 1:nregions) {
        	mydata$fit$epi[iweek, iregion] = max(mydata$fit$epi[iweek, iregion], 1)
        	mydata$fit$epi[iweek, iregion] = round(mydata$fit$epi[iweek, iregion])
        }
    }

    cat("\n Giving Future Data Points a weight of", model_wght_ftr, "\n\n")
	cat('\n Data Shifted by: ', shft, '\n\n')

    mydata$model$wght = rep(0, nperiods)

    mydata$model$wght[1:nperiodsFit] = wght_data
    mydata$model$wght[(nperiodsFit + 1):nperiods] = model_wght_ftr

    mydata$fit$gamaepi = lgamma(mydata$fit$epi + 1)
    mydata$model$gamaepi = lgamma(mydata$model$epi + 1)

    return(mydata)
}
