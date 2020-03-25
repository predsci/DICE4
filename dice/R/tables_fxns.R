##
## Functions that set up tables for results and calculates the error in the fit/forecast
##
results.table <- function(state_names = NULL, state_id = NULL, cadence.names = NULL, cadence = "month", nperiods = NULL) {
    #' Create a List of Tables for Observations and Results
    #'
    #' \code{results.table} Creates all the tables used to store the observed data and
    #' the results/errors of the SIR/SEIR and SARIMA fits/forecasts
    #' @param state_names The full name of the states/regions
    #' @param state_id The abbreviation of the states/regions
    #' @param cadence.names An array with month or week names
    #' @param cadence String, the cadence of the data: month, week
    #' @param nperiods Integer - the number of data points in the season
    #' @return tables A list with a table for each observation/model/error measure.
    #' @examples
    #' results.table(state_names = state_names, state_id = state_id, cadence.names = cadence.names, cadennce = cadence)
    #'



    Period = cadence_per_year = nperiods

    nstates = length(state_names)

    epi.obsrv = array(NA, c(cadence_per_year, nstates))

    epi.null = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    epi.null.mae = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    epi.null.mre = array(NA, c(cadence_per_year, cadence_per_year, nstates))


    seir.frcst = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    seir.mae = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    seir.mre = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    seir.rel = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    seir.lower = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    seir.upper = array(NA, c(cadence_per_year, cadence_per_year, nstates))

    arima.frcst = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    arima.mae = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    arima.mre = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    arima.rel = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    arima.lower = array(NA, c(cadence_per_year, cadence_per_year, nstates))
    arima.upper = array(NA, c(cadence_per_year, cadence_per_year, nstates))

    dimnames(epi.obsrv)[[1]] = cadence.names
    dimnames(epi.obsrv)[[2]] = as.list(state_id)

    dimnames(arima.frcst)[[1]] = dimnames(arima.mae)[[1]] = dimnames(arima.mre)[[1]] = dimnames(arima.rel)[[1]] = dimnames(arima.lower)[[1]] = dimnames(arima.upper)[[1]] = dimnames(epi.null)[[1]] = dimnames(epi.null.mae)[[1]] = dimnames(epi.null.mre)[[1]] = dimnames(seir.frcst)[[1]] = dimnames(seir.mae)[[1]] = dimnames(seir.mre)[[1]] = dimnames(seir.rel)[[1]] = dimnames(seir.lower)[[1]] = dimnames(seir.upper)[[1]] = cadence.names

    dimnames(arima.frcst)[[2]] = dimnames(arima.mae)[[2]] = dimnames(arima.mre)[[2]] = dimnames(arima.rel)[[2]] = dimnames(arima.lower)[[2]] = dimnames(arima.upper)[[2]] = dimnames(epi.null)[[2]] = dimnames(epi.null.mae)[[2]] = dimnames(epi.null.mre)[[2]] = dimnames(seir.frcst)[[2]] = dimnames(seir.mae)[[2]] = dimnames(seir.mre)[[2]] = dimnames(seir.rel)[[2]] = dimnames(seir.lower)[[2]] = dimnames(seir.upper)[[2]] = cadence.names

    dimnames(arima.frcst)[[3]] = dimnames(arima.mae)[[3]] = dimnames(arima.mre)[[3]] = dimnames(arima.rel)[[3]] = dimnames(arima.lower)[[3]] = dimnames(arima.upper)[[3]] = dimnames(epi.null)[[3]] = dimnames(epi.null.mae)[[3]] = dimnames(epi.null.mre)[[3]] = dimnames(seir.frcst)[[3]] = dimnames(seir.mae)[[3]] = dimnames(seir.mre)[[3]] = dimnames(seir.rel)[[3]] = dimnames(seir.lower)[[3]] = dimnames(seir.upper)[[3]] = as.list(state_id)


    return(list(epi.obsrv = epi.obsrv, arima.frcst = arima.frcst, arima.mae = arima.mae, arima.mre = arima.mre,
        arima.rel = arima.rel, arima.lower = arima.lower, arima.upper = arima.upper, epi.null = epi.null,
        epi.null.mae = epi.null.mae, epi.null.mre = epi.null.mre, seir.frcst = seir.frcst, seir.mae = seir.mae,
        seir.mre = seir.mre, seir.rel = seir.rel, seir.lower = seir.lower, seir.upper = seir.upper))

}

calc.mech.err <- function(tables = NULL, profiles = NULL, nfit = NULL, state_id = NULL) {
    #' Calculate Error in Mechanistic Fit/Forecast
    #'
    #' \code{calc.mech.err} calculates the error:
    #' mean absoulte error, mean relative error and mean error relative to the NULL model for each
    #' month/week for which a fit/forecast was calculated
    #' @param tables A list with all the data tables - both observed and fitted/forecasted
    #' @param profiles An array with 10,000 selected dengue profiles from the MCMC fitting procedure
    #' @param nfit An integer - the number of data points used in the fit
    #' @param state_id The abbreviation of the state/region
    #' @return tables An updated tables list with the \code{nfit} row of the SEIR errors filled
    #' @examples
    #' calc.mech.err(tables=tables, profiles=profiles, nfit= nfit, state_id = state_id)


    epi.obsrv = tables$epi.obsrv

    epi.null = tables$epi.null
    epi.null.mae = tables$epi.null.mae
    epi.null.mre = tables$epi.null.mre

    seir.frcst = tables$seir.frcst
    seir.mae = tables$seir.mae
    seir.mre = tables$seir.mre
    seir.rel = tables$seir.rel
    seir.lower = tables$seir.lower
    seir.upper = tables$seir.upper

    dims = dim(seir.frcst)
    cadence_per_year = dims[2]

    ntraj = dim(profiles)[2]

	nstates = length(state_id)

	for (istate in 1:nstates) {
        id = state_id[istate]
        seir.mae[nfit, 1:cadence_per_year, id] = 0
        for (itraj in 1:ntraj) {
            seir.mae[nfit, 1:cadence_per_year, id] = seir.mae[nfit, 1:cadence_per_year, id] + abs(profiles[itraj, 1:cadence_per_year] - 	epi.obsrv[1:cadence_per_year, id])
        }

        seir.mae[nfit, 1:cadence_per_year, id] = seir.mae[nfit, 1:cadence_per_year, id]/ntraj
        seir.mre[nfit, 1:cadence_per_year, id] = seir.mae[nfit, 1:cadence_per_year, id]/epi.obsrv[1:cadence_per_year, id]
        seir.rel[nfit, 1:cadence_per_year, id] = seir.mae[nfit, 1:cadence_per_year, id]/epi.null.mae[nfit, 1:cadence_per_year,
            id]

        for (imonth in 1:cadence_per_year) {
            seir.lower[nfit, imonth, id] = quantile(profiles[,imonth], prob = 0.05)
            seir.upper[nfit, imonth, id] = quantile(profiles[,imonth], prob = 0.95)
            seir.frcst[nfit, imonth, id] = mean(profiles[,imonth])
        }

    }
        
    tables$seir.frcst = seir.frcst
    tables$seir.mae = seir.mae
    tables$seir.mre = seir.mre
    tables$seir.rel = seir.rel
    tables$seir.lower = seir.lower
    tables$seir.upper = seir.upper

    return(tables)

}

calc.stat.err <- function(tables = NULL, nfit = NULL, state_id = NULL) {
    #' Calculate Error in SARIMA Fit/Forecast
    #'
    #' \code{calc.stat.err} calculates the error:
    #' mean absolute error, mean relative error and mean error relative to the NULL model for each
    #' month/week for which a fit/forecast was calculated
    #' @param tables A list with all the mydata tables - both observed and fitted/forecasted
    #' @param nfit An integer - the number of data points used in the fit
    #' @param state_id The abbreviation of the state/region
    #' @return tables An updated tables list with the \code{nfit} row of the SEIR errors filled
    #' @examples
    #' calc.stat.err(tables=tables, nfit= nfit, state_id = state_id)


    epi.obsrv = tables$epi.obsrv

    epi.null = tables$epi.null
    epi.null.mae = tables$epi.null.mae
    epi.null.mre = tables$epi.null.mre

    arima.frcst = tables$arima.frcst
    arima.mae = tables$arima.mae
    arima.mre = tables$arima.mre
    arima.rel = tables$arima.rel

    dims = dim(arima.frcst)
    cadence_per_year = dims[2]

    nstates = length(state_id)

    for (istate in 1:nstates) {
    	id = state_id[istate]
    	arima.mae[nfit, 1:cadence_per_year, id] = abs(arima.frcst[nfit, 1:cadence_per_year, id] - epi.obsrv[1:cadence_per_year, id])
    	arima.mre[nfit, 1:cadence_per_year, id] = arima.mae[nfit, 1:cadence_per_year, id]/epi.obsrv[1:cadence_per_year, id]
    	arima.rel[nfit, 1:cadence_per_year, id] = arima.mae[nfit, 1:cadence_per_year, id]/epi.null.mae[nfit, 1:cadence_per_year, id]
    }

    tables$arima.frcst[nfit,, state_id] = arima.frcst[nfit,, state_id]
    tables$arima.mae[nfit,, state_id]   = arima.mae[nfit,, state_id]
    tables$arima.mre[nfit,, state_id]   = arima.mre[nfit,, state_id]
    tables$arima.rel[nfit,, state_id]   = arima.rel[nfit,, state_id]

    return(tables)

}

calc.null.err <- function(tables = NULL, nfit = NULL, state_id = NULL) {
    #' Calculate Error of Historic Null Model
    #'
    #' \code{calc.null.err} calculates the error of the historic NULL model:
    #' the weekly/monthly value is given by the ten year average value
    #' @param tables A list with all the data tables - both observed and fitted/forecasted
    #' @param nfit An integer - the number of data points used in the fit
    #' @param state_id The abbreviation of the state/region
    #' @return tables An updated tables list with the \code{nfit} row of the NULL errors
    #' - absolute and relative error
    #' filled
    #' @examples
    #' calc.null.err(tables=tables, nfit= nfit, state_id = state_id)


    epi.obsrv = tables$epi.obsrv

    epi.null = tables$epi.null
    epi.null.mae = tables$epi.null.mae
    epi.null.mre = tables$epi.null.mre

    dims = dim(epi.null)
    cadence_per_year = dims[2]

    nstates = length(state_id)

    for (istate in 1:nstates) {
    	id = state_id[istate]
    	epi.null.mae[nfit, 1:cadence_per_year, id] = abs(epi.null[nfit, 1:cadence_per_year, id] - epi.obsrv[1:cadence_per_year, id])
    	epi.null.mre[nfit, 1:cadence_per_year, id] = epi.null.mae[nfit, 1:cadence_per_year, id]/epi.obsrv[1:cadence_per_year, id]
    	epi.null.mae[nfit,is.nan(epi.null.mae[nfit,, id]), id] <- NA
    	epi.null.mre[nfit,is.nan(epi.null.mae[nfit,, id]), id] <- NA
    }


    tables$epi.null.mae = epi.null.mae
    tables$epi.null.mre = epi.null.mre


    return(tables)

}


