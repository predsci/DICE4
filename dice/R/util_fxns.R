###
### Utility functions for DICE
###

print.error.dice <- function(err, ...) {

    #' Print an Error Message
    #'
    #' Prints the error message and quits \pkg{DICE} with save option set to yes
    #' @param err A character string with the error message
    #' @return No return value. Quit \pkg{DICE} with the save option set to yes
    #' @examples
    #' print.error.dice('For GFT mydata sets DICE only supports flu seasons 2003-2014')

    cat("Error in Dice: \n")
    print(err)
    cat("\n")
    quit(save = "yes")

}

set.iseed <- function() {

    #' Set the seed for the RNG
    #'
    #' Create the seed in the R code for the Random Number Generator
    #' @param none required
    #' @return An integer - seed for the RNG
    #' @examples
    #' set.iseed()

    set.seed(seed = NULL)
    iseed <- runif(1) * 10000
    iseed <- round(iseed)

    return(iseed)
}


makeDir <- function(subDir = NULL) {
    #' Create a sub-directory
    #'
    #' Create a sub-directory for a \pkg{DICE} run under the current working directory.   If this sub-directory already exists it does NOT overwrite it.
    #' @param subDir String - the sub-directory name for all output files of the run. Default it 'output'
    #' @examples
    #' makeDir{subDir='output'}

    if (is.null(subDir))
        subDir = "output"
    mainDir = getwd()
    if (!file.exists(subDir)) {
        dir.create(file.path(mainDir, subDir))
    }

    return(err = 0)
}

gammaEpi <- function(epi = NULL) {

    #' Calculate the \eqn{Gamma} Function of Disease Number of Cases
    #'
    #' Given a matrix or vector of integers with monthly/weekly/daily number of disease cases calculate for each spatial region, the \eqn{\Gamma} function of this matrix/vector.
    #' This is done to save computational time during the MCMC procedure: the \eqn{\Gamma} function of the number of cases is needed when evaluating the maximum likelihood.
    #' @param epi A nperiods x nregions matrix of integers  with monthly/weekly/daily number of disease cases for each region
    #' @return gamaepi A matrix  of nperiods x nregions with the \eqn{\Gamma} function of monthly/weekly/daily number of disease cases in each region
    #' @examples
    #' gammaEpi{epi = mydata$model$epi$}
    #' gammaEpi{epi = mydata$fit$epi$  }

    copyepi = epi
    gamaepi = epi

    n = dim(epi)[2]
    nperiods = dim(epi)[1]

    # to save computational time we calculate the gamma(epi+1)
    for (i in 1:nperiods) {
        for (j in 1:n) {
            gamaepi[i, j] <- lgamma((epi[i, j] + 1))
        }
    }
    return(gamaepi)
}

distance.matrix <- function(mydata = NULL) {

    #' Calculate a Euclidean Distance Matrix
    #'
    #' Given a data frame with the centroid location (in lat/lon) of all the regions
    #'   construct the corresponding Euclidian distance matrix.
    #' @param mydata The complete list with all the information for the run.
    #' @return A 2D matrix with the Euclidian distance between the centroids of the fit
    #'   regions. (Units are km.)
    #' @examples
    #' distance.matrix{mydata=mydata}

    n = mydata$fit$nregions
    mydistance <- array(0, c(n, n))
    fit_level = mydata$fit_level

    rownames(mydistance) <- mydata$fit$attr[[paste0("ABBV_", fit_level)]]
    colnames(mydistance) <- mydata$fit$attr[[paste0("ABBV_", fit_level)]]

    # # spDistsN1 with longlat=TRUE will give us the great sphere distance from 'pt' to the matrix of long/lat units are km the diagonal
    # cntr.matrix = cbind(mydata$fit$attr$sedac_lon, mydata$fit$attr$sedac_lat)
    # colnames(cntr.matrix) = c("sedac_lon", "sedac_lat")
    # # elements of dist0 are zero.
    # for (j in 1:n) {
    #     mydistance[j, ] <- spDistsN1(pts = cntr.matrix, pt = cntr.matrix[j, ], longlat = TRUE)
    # }

    # Hard-code the great sphere calc to remove sp dependency
    cntr.matrix = cbind(mydata$fit$attr$sedac_lat, mydata$fit$attr$sedac_lon)
    colnames(cntr.matrix) = c("sedac_lat", "sedac_lon")
    # convert lat/lon to spherical radians
    cntr.rads = pi*cntr.matrix/180
    mean_earth_radius_km = 6371
    # calc great circle distance - by haversine formula
    for (jj in 1:n) {
      for (kk in 1:jj) {
        del_lon = cntr.rads[jj, 2] - cntr.rads[kk, 2]
        del_lat = cntr.rads[jj, 1] - cntr.rads[kk, 1]
        intermidiate_a = sin(del_lat/2)^2 + cos(cntr.rads[jj, 1])*cos(cntr.rads[kk, 1])*sin(del_lon/2)^2
        central_angle = 2*asin(min(1,sqrt(intermidiate_a)))
        mydistance[jj, kk] = mean_earth_radius_km*central_angle
      }
    }
    mydistance = mydistance + t(mydistance)

    return(mydistance)
}

Date2CDCweek <- function(date=NULL) {
    #' Convert a Date to CDC Week and Year
    #'
    #' CDC week numbering: first Wednesday of the year is in week 1.
    #' @param date A `Date' class object
    #' @return A list containing the CDC week and year that `date' falls in.
    #' @examples
    #' require(DICE)
    #' Date2CDCweek(as.Date('2011-03-15'))
    #'

    # extract year and week from date
    year = as.integer(format(date, "%Y"))
    Rweek = as.integer(format(date, "%U"))
    DayOfYear = as.integer(format(date, "%j"))
    # shift R-weeks to CDCweeks based off first day of the year
    Jan1 = weekdays(as.Date(paste0(year, "-01-01")))

    # check if date belongs to previous year
    if (DayOfYear < 4) {
        if ((Jan1 == "Thursday") || (Jan1 == "Friday" && DayOfYear < 3) || (Jan1 == "Saturday" && DayOfYear == 1)) {
            year = year - 1
            Jan1 = weekdays(as.Date(paste0(year, "-01-01")))
            Rweek = as.integer(format(date - 7, "%U")) + 1
        }
    }

    # convert Rweeks to CDCweeks
    if (any(Jan1 == c("Monday", "Tuesday", "Wednesday"))) {
        CDCweek = Rweek + 1
    } else CDCweek = Rweek

    # check if date belongs to the following year
    if (DayOfYear > 362) {
        NumDays = as.integer(format(as.Date(paste0(year, "-12-31")), "%j"))
        Jan1next = weekdays(as.Date(paste0(year + 1, "-01-01")))
        if ((Jan1next == "Wednesday" && NumDays - DayOfYear < 3) || (Jan1next == "Tuesday" && NumDays - DayOfYear < 2) || (Jan1next ==
            "Monday" && NumDays - DayOfYear == 0)) {
            year = year + 1
            CDCweek = 1
        }
    }

    list(CDCweek = CDCweek, CDCyear = year)
}


CDCweek2date <- function(CDCweek, year) {
    #' Convert a CDC week and year to Sunday's Date.
    #'
    #' CDC week numbering: first Wednesday of the year is in week 1.
    #' @param CDCweek An integer indicating the CDC week to be converted
    #' @param year An integer specifying the year
    #' @return A `Date' class object containing the date of Sunday for the week and year provided.
    #' @examples
    #' require(DICE)
    #' CDCweek2date(10,2011)

    # Determine if CDCweeks and Rweeks are the same using 1-Jan-year
    Jan1 = weekdays(as.Date(paste0(year, "-01-01")))
    if (any(Jan1 == c("Monday", "Tuesday", "Wednesday"))) {
        # They are not, shift by 1
        Rweek = CDCweek - 1
        if (Rweek == 0) {
            date = as.Date(paste0("0-", 1, "-", year), format = "%w-%U-%Y") - 7
        } else {
            date = as.Date(paste0("0-", Rweek, "-", year), format = "%w-%U-%Y")
        }
    } else {
        # this year Rweeks and CDCweeks are the same
        Rweek = CDCweek
        if (Rweek == 53) {
            date = as.Date(paste0("0-", 52, "-", year), format = "%w-%U-%Y") + 7
        } else {
            date = as.Date(paste0("0-", Rweek, "-", year), format = "%w-%U-%Y")
        }
    }
    return(date)
}

calcLLK <- function(mod, dat) {
    #' Calculate the Likelihood of each Data Point
    #'
    #' Calculate the likelihood of each data point given two time-series: the observed one and predicted (dat and mod respectively).
    #' @param mod Numeric time series with the predicted number of cases
    #' @param dat Integer time series with the observed number of cases
    #' @return myLLK  A vector with the likelihood value of each data point
    #' @examples
    #' calcLLK(mod=model_rtn, dat = mydata$model$epi)
    #' calcLLK(mod=model_rtn, dat = mydata$fit$epi[,iregion])
    myLLK <- dpois(dat, mod, log = TRUE)
    return(myLLK)
}

psi.lnlike <- function(mod, dat) {
    #' Calculate the Likelihood of a Model
    #'
    #' Calculate the total likelihood given two time-series: the observed one and predicted (dat and mod respectively).
    #' @param mod Numeric time series with the predicted number of cases
    #' @param dat Integer time series with the observed number of cases
    #' @return myLLK  The total likelihood of all the data points
    #' @examples
    #' psi.lnlike(mod=model_rtn, dat = mydata$model$epi)
    #' psi.lnlike(mod=model_rtn, dat = mydata$fit$epi[,iregion])
    #'
    myLLK <- dpois(dat, mod, log = TRUE)
    sum(myLLK)
}


psi.lnlike.deviance <- function(mod, dat) {
    #' Calculate the Deviance of a Model
    #'
    #' Calculate the deviance given two time-series: the observed one and predicted (dat and mod respectively).
    #' @param mod Numeric time series with the predicted number of cases
    #' @param dat Integer time series with the observed number of cases
    #' @return devi the deviance of the model
    #' @examples
    #' psi.lnlike.deviance(mod=model_rtn, dat = mydata$model$epi)
    #' psi.lnlike.deviance(mod=model_rtn, dat = mydata$fit$epi[,iregion])
    #'
    fit <- dpois(dat, mod, log = TRUE)
    sat <- dpois(dat, dat, log = TRUE)

    devi = -2 * (sum(fit) - sum(sat))
    devi
}

psi.lnlike.deldeviance <- function(modI, modJ, dat) {

    #' Calculate the Difference Deviance of Two models
    #'
    #' Calculate the difference in deviance given three time-series: the observed one and two predicted ones (dat, modI and modJ respectively).
    #' @param modI Numeric time series with the predicted number of cases using model I
    #' @param modJ Numeric time series with the predicted number of cases using model J
    #' @param dat Integer time series with the observed number of cases
    #' @return deldevi the difference in deviance of the two models

    fitJ <- dpois(dat, modJ, log = TRUE)
    fitI <- dpois(dat, modI, log = TRUE)

    deldevi = -2 * (sum(fitI) - sum(fitJ))
    deldevi
}


quanci <- function(x) {

    #' Calculate the 0.75 confidence interval of an a distribution
    #'
    #' @param x - Numeric, a one dimensional vector of numbers
    #' @return the 0.75 lower and upper bounds

    ub = quantile(x, 0.875)
    lb = quantile(x, 0.125)
    return(c(lb, ub))
}

