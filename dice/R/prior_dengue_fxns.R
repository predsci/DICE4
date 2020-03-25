##
## Functions related to Data Augmentation for Dengue: calculate correlation with past season, distance matrix, historic NULL model
## and actual aumentation
##


calc.epi.cor   <- function(mydata = NULL, all_years_epi = NULL) {
	# This is the season
	my.year = mydata$year
	nregion = mydata$fit$nregion

   	mod_level = mydata$mod_level
    fit_level = mydata$fit_level

    model_raw = mydata$model$raw
    fit_raw = mydata$fit$raw
    nregion = mydata$fit$nregions

    nregion1 = nregion + 1

	
	 
}
calc.sqldata.cor <- function(mydata= NULL, all_years_epi = NULL, cases = NULL, nfit = NULL) {

    #' Calculate Pearson Correlation Matrix Between Seasons
    #'
    #' \code{calc.sqldata.cor} Calculates an nyears by nyears correlation matrix for dengue incidence
    #' @param mydata A dataframe with all the available data for this \pkg{DICE} run
    #' @param all_years_epi the epi data for all years
    #' @param epi the epi data for current season
    #' @param nfit  The number of data points used when calculating the correlation
    #' @param cases  AN array with the entire time series of incidence
    #' @return corrMAT An nyear x nyears matrix with the Pearson correlation values
    #' @examples
    #' calc.cor(mydata= mydata, all_years_epi = al_years_epi, cases = all_years_epi$model$epi, nfit = nfit)
    #' calc.cor(mydata= mydata, all_years_epi = al_years_epi, cases = all_years_epi$fit$epi[,iregion], nfit = nfit)

    year.vec = all_years_epi$years
    years = unique(year.vec)

 	nyears = length(years)

    nperiods = mydata$nperiods
    nperiodsData = mydata$nperiodsData
    nperiodsFit  = mydata$nperiodsFit

    if (mydata$cadence == 'Monthly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$months == mydata$months[1])
    	cadence.vec = all_years_epi$months
    	cadence.start = mydata$months[1]
    	cadence = 'month'

    }

    if (mydata$cadence == 'Weekly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$weeks[1])
    	cadence.vec = all_years_epi$weeks
    	cadence.start = mydata$weeks[1]
    	cadence = 'week'
    }

   if (mydata$cadence == 'Daily') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$days[1])
    	cadence.vec = all_years_epi$days
    	cadence.start = mydata$days[1]
    	cadence = 'day'
    }
    
    year.start = mydata$years[1]
    istart = which(all_years_epi$years == year.start & cadence.vec == cadence.start)
    nlines = length(all_years_epi$years)
    # This is the epi data we will calculate correlation with
    
    # Now walk back in time 
    
    corrMAT = array(NA, c(nyears, nyears))

    colnames(corrMAT) = rownames(corrMAT) = years
	
	for (iyear in rev(1:(nyears))) {
		istart = which(all_years_epi$years == (year.start-nyears+iyear) & cadence.vec == cadence.start)
		if(length(istart) == 0) next
		iend   = istart + nperiods -1
		if(iend > nlines) iend = nlines
		x = cases[istart:iend]
		if (all(is.na(x))) next

		for (jyear in rev(1:nyears)) {
			if (iyear == jyear) next

			istart = which(all_years_epi$years == (year.start-nyears+jyear) & cadence.vec == cadence.start)
			if(length(istart) == 0) next
			iend   = istart + nperiods -1
			if(iend > nlines) iend = nlines
			y = cases[istart:iend]
			if(all(is.na(y))) next

			corrMAT[iyear,jyear] = 	cor(x[1:nfit], y[1:nfit], method = "pearson", use = "pairwise.complete.obs")		
		}		
	}


    return(corrMAT)

}

calc.sqldata.dist <- function(mydata= NULL, all_years_epi = NULL, cases = NULL, nfit = NULL) {


    #' Calculate Euler Distance Matrix Between Seasons
    #'
    #' \code{calc.dist} Calculates an nyears X nyears Euler distance matrix for the incidence profiles
     #' @param mydata A dataframe with all the available data for this \pkg{DICE} run
    #' @param all_years_epi the epi data for all years
    #' @param nfit  The number of data points used when calculating the correlation
    #' @param cases  AN array with the entire time series of incidence
    #' @return distMAT An nyear x nyears matrix with Euler distance values
    #' @examples
    #' calc.sqldata.dist(mydata= mydata, all_years_epi = al_years_epi, cases = all_years_epi$model$epi, nfit = nfit)
	#' calc.sqldata.distr(mydata= mydata, all_years_epi = al_years_epi, cases = all_years_epi$fit$epi[,iregion], nfit = nfit)


 
    year.vec = all_years_epi$years
    years = unique(year.vec)

 	nyears = length(years)

    nperiods = mydata$nperiods
    nperiodsData = mydata$nperiodsData
    nperiodsFit  = mydata$nperiodsFit

    if (mydata$cadence == 'Monthly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$months == mydata$months[1])
    	cadence.vec = all_years_epi$months
    	cadence.start = mydata$months[1]
    	cadence = 'month'

    }

    if (mydata$cadence == 'Weekly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$weeks[1])
    	cadence.vec = all_years_epi$weeks
    	cadence.start = mydata$weeks[1]
    	cadence = 'week'
    }

   if (mydata$cadence == 'Daily') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$days[1])
    	cadence.vec = all_years_epi$days
    	cadence.start = mydata$days[1]
    	cadence = 'day'
    }
    
    year.start = mydata$years[1]
    istart = which(all_years_epi$years == year.start & cadence.vec == cadence.start)
    nlines = length(all_years_epi$years)
    # This is the epi data we will calculate correlation with
    
    # Now walk back in time 

    distMAT = array(NA, c(nyears, nyears))

    colnames(distMAT) = rownames(distMAT) = years

	for (iyear in rev(1:(nyears))) {
		istart = which(all_years_epi$years == (year.start-nyears+iyear) & cadence.vec == cadence.start)
		if(length(istart) == 0) next
		iend   = istart + nperiods -1
		if(iend > nlines) iend = nlines
		x = cases[istart:iend]
		if (all(is.na(x))) next

		for (jyear in rev(1:nyears)) {
			
			istart = which(all_years_epi$years == (year.start-nyears+jyear) & cadence.vec == cadence.start)
			if(length(istart) == 0) next
			iend   = istart + nperiods -1
			if(iend > nlines) iend = nlines
			y = cases[istart:iend]
			if(all(is.na(y))) next

			distMAT[iyear,jyear] = 	sqrt(sum((as.numeric(x[1:nfit] - y[1:nfit])) * (as.numeric(x[1:nfit] - y[1:nfit]))))/nfit	
		}		
	}


    return(distMAT)

}


get.sql.da <- function(nfit = NULL, mydata = NULL, all_years_epi = NULL, cases = NULL, epi = NULL,
    corrMAT = NULL, deng.null = NULL, my_id = NULL) {
    #'
    #' Data Augmentation for Mechanistic Model
    #'
    #' \code{get.sql.da} Allows the User to augment the incidence time series beyond the \code{nfit} data points
    #' using either the historic NULL model \code{da = 1}, or the most similar past season \code{da=2}.
    #' If \code{da=0}, the data is not augmented.
    #' @param mydata A dataframe with all the available data for this \pkg{DICE} run
    #' @param all_years_epi the epi data for all years
    #' @param cases An array with the entire history of incidence information for one region/patch
    #' @param epi An array with the incidence information for one season and one region/patch
    #' @param nfit An integer - the number of data points use in the run
    #' @param my_id The abbreviation of the state/region
    #' @param corrMAT A 2D array with the Pearson corrletion values for all available disease seasons
    #' @param deng.null the Historic NULL model - average monthly/weekly number of cases
    #' @return epi A list with the original and augmented mydata (if da=1,2) and the weight for the augmented
    #' data points
    #' @examples
    #' get.sql.da(nfit = nfit, mydata = mydata, years = all_years_epi$years, cases = all_years_epi$model$epi,
    #'  epi = mydata$model$epi,corrMAT = corrMAT, distMAT = distMAT, deng.null=deng.null)
    #'
    #' get.sql.da(nfit = nfit, mydata = mydata, years = all_years_epi$years, cases = all_years_epi$fit$epi[,iregion],
    #'  epi = mydata$fit$epi[,iregion],corrMAT = corrMAT, distMAT = distMAT, deng.null=deng.null)

   	if (is.null(nfit) | is.null(cases) | is.null(epi))
        return


    da = mydata$da

    if (is.null(da))
        da = 0

    cadence = mydata$cadence
    nperiods = mydata$nperiods
	nperiodsFit = mydata$nperiodsFit
	
    if (mydata$cadence == 'Monthly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$months == mydata$months[1])
    	cadence.vec = all_years_epi$months
    	cadence.start = mydata$months[1]
    	cadence = 'month'

    }

    if (mydata$cadence == 'Weekly') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$weeks[1])
    	cadence.vec = all_years_epi$weeks
    	cadence.start = mydata$weeks[1]
    	cadence = 'week'
    }

   if (mydata$cadence == 'Daily') {
    	my_row = which(all_years_epi$years == mydata$years[1] & all_years_epi$weeks == mydata$days[1])
    	cadence.vec = all_years_epi$days
    	cadence.start = mydata$days[1]
    	cadence = 'day'
    }

    augment = FALSE

    if (da > 0)
        augment = TRUE

    epi.state =  epi

    nperiods = mydata$nperiods

    nda = nperiods - nfit

    if (nda == 0)
        augment = FALSE


    if (augment) {
        cat("\n Entered da == TRUE Option \n\n")
        epi.save = epi
        epi.state.save = epi.state
        my.year = mydata$season
        print(my.year)
        iyear = which(rownames(corrMAT) == my.year)
		iyear = as.numeric(iyear)
		
		
        if (all(is.na(corrMAT[iyear, ])) | da == 1) {
            corr.max = -1
        } else {
            corr.max = max(corrMAT[iyear, ], na.rm = TRUE)
        }


        if (corr.max > 0) {  ## Looking for the Most similar year 
        	
            jyear = which.max(corrMAT[iyear, ])
        	
    		year.start = mydata$years[1]

    		nyears = length(unique(all_years_epi$years))
 
 			istart = which(all_years_epi$years == (year.start-nyears+jyear) & cadence.vec == cadence.start)
			#if(length(istart) == 0) 
			iend   = istart + mydata$nperiods -1
			
			nlines = length(cases)
			if(iend > nlines) iend = nlines
			y = cases[istart:iend]
 
			futur = cases[istart:iend]
			
            shft = epi.state[nfit] - futur[nfit]

            x = futur[1:nfit]

            y = epi.state[1:nfit]

			shft = y[nfit] - x[nfit]
			
			prob = cor(x, y, method = "pearson", use = "pairwise.complete.obs")
			
            #b1 = sum((x - mean(x)) * (y - mean(y)))/sum((x - mean(x)) * (x - mean(x)))

            #b0 = mean(y) - b1 * mean(x)

            #z = b0 + b1 * futur

            cat("\n Using Pearson Correlation to Choose Season for Data Augmentation \n")
            cat("\n Season Selected is: ", colnames(corrMAT)[jyear], "\n")


        } else {

            cat(" Using Historic NULL model to Augment Data \n")

            futur = deng.null[,my_id]

            x = futur[1:nfit]
            y = epi.state[1:nfit]

            shft = epi.state[nfit] - futur[nfit]

			prob = cor(x, y, method = "pearson", use = "pairwise.complete.obs")
			
            #b1 = sum((x - mean(x)) * (y - mean(y)))/sum((x - mean(x)) * (x - mean(x)))

            #b0 = mean(y) - b1 * mean(x)

            #z = b0 + b1 * futur

        }

        epi.new = epi.state

        wght = rep(1, nperiods)

        if (nfit < nperiods) epi.new[(nfit + 1):nperiods] = futur[(nfit + 1):nperiods] + shft
           	
        if (prob <= 0 | length(prob) == 0) prob = (nperiods-nfit)/nperiods

        if (nfit < nperiods)  {
        	wght[(nfit + 1):nperiods] = prob ##cor(x[1:nfit], y[1:nfit])

        	epi.new[(nfit + 1):nperiods] = round(epi.new[(nfit + 1):nperiods])
        }

        for (i in (nfit + 1):nperiods) {
            epi.new[i] = max(1, epi.new[i])
        }

    } else {

        cat("\n\n NO Data Augmentation \n\n")
        epi.save = epi
        epi.state.save = epi.state
        epi.new = epi.state
        wght = rep(0, nperiods)
        wght[1:nperiods] = 1

    }

    return(list(nfit = nfit, epi.new = epi.new, wght = wght))

}
