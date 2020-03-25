fitSARIMA <- function(mydata = NULL, all_years_epi = NULL, arima_model = NULL, auto_arima_model = NULL, device = 'png') {
	#'
	#' SARIMA Fitting to Disease Data
	#'
	#'\code{fitSARIMA} fits the model and fit level regions with a SARIMA model (either chosen by the user or using \pkg{auto.arima})
	#' and calls the plotting function
	#' The model level is fitted both directly and as the weighted average of the fit level regions
	#' @param mydata - dataframe with all the data for this \pkg{DICE} run
	#' @param all_years_epi the epi data for all years
	#' @param arima_model List of ARIMA model parameters: list(p=, d=, q=, P=, D, Q=) can be set to NULL to trigger the
	#' @param device - 'pdf' or 'png'. Default is 'png'
	#' \code{auto.arima} process
	#' @param auto_arima_model A list of upper limit values for the ARIMA parameters.  This list is created by the
	#' \code{DICE} code ONLY if arima_model is not set by the user
	#' @examples
	#' fitSARIMA((mydata = mydata, all_years_epi = all_years_epi, arima_model = arima_model,
	#' auto_arima_model = NULL, covar = covar, covar_lag = covar_lag, run.list = run.list)
	#' @return
	#' A list of tables with the results for the model (direct and aggregate) and fit regions
	#'


	mod_level = mydata$mod_level
	fit_level = mydata$fit_level

	cadence = mydata$cadence

	nfit = mydata$nperiodsFit

	mod_id = mydata$model$attr[[paste0("ABBV_", mod_level)]]
	fit_id = mydata$fit$attr[[paste0("ABBV_", fit_level)]]

	mod_name = mydata$model$name
	fit_name = mydata$fit$name
	nregion = mydata$fit$nregions
	nregion1 = nregion + 1
	nregion2 = nregion1 + 1

	covar = mydata$covar
	covar_lag = mydata$covar_lag

	if (mydata$cadence == "Weekly") {
		cadence = "week"
		nperiods = mydata$nperiods
	} else if (mydata$cadence == "Monthly") {
		cadence = "month"
		nperiods = mydata$nperiods
	} else if (mydata$cadence == "Daily") {
		cadence = "day"
		nperiods = mydata$nperiods
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}

	## Prepare a list for the arima models - really needed only when arima_model = NULL and each region can have a  different
	## arima model coming out of auto.arima

	arima_model_all <- list()

	for (iregion in 1:nregion2) {
		arima_model_all[[iregion]] <- arima_model
	}

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)

	epi.obsrv = mydata$model$raw

	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)

	epi.null.mod = epi.null$cases.ave.mod

	tables.mod$epi.null[nfit, , mod_id] = epi.null.mod[, mod_id]


	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables


	if (mydata$mod_level < mydata$fit_level) {

		tables.fit <- results.table(state_names = fit_name, state_id = fit_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)

		epi.obsrv = mydata$fit$raw

		tables.fit$epi.obsrv[, fit_id] = as.matrix(epi.obsrv)

		epi.null.fit = epi.null$cases.ave.fit
		tables.fit$epi.null[nfit, , fit_id] = epi.null.fit[, fit_id]

		tables = calc.null.err(tables = tables.fit, nfit = nfit, state_id = fit_id)
		tables.fit = tables

	}

	if (mydata$cadence == "Monthly") {
		my_row = which(all_years_epi$years == mydata$years[nfit] & all_years_epi$months == mydata$months[nfit])
		cadence.vec = all_years_epi$month
		my.freq = 12
	}

	if (mydata$cadence == "Weekly") {
		my_row = which(all_years_epi$years == mydata$years[nfit] & all_years_epi$weeks == mydata$weeks[nfit])
		cadence.vec = all_years_epi$weeks
		my.freq = 52
	}

	years.vec = all_years_epi$years
	year.vec = mydata$years

	attr = c("year", cadence, "ndays")

	istart = my_row - my.freq * 10
	istart = max(istart, 1)
	ntrn = length(istart:my_row)

	if (!is.null(arima_model)) {
		p = as.numeric(arima_model[["p"]])
		d = as.numeric(arima_model[["d"]])
		q = as.numeric(arima_model[["q"]])

		P = as.numeric(arima_model[["P"]])
		D = as.numeric(arima_model[["D"]])
		Q = as.numeric(arima_model[["Q"]])
	} else {
		max.p = as.numeric(auto_arima_model[["p"]])
		max.d = as.numeric(auto_arima_model[["d"]])
		max.q = as.numeric(auto_arima_model[["q"]])

		max.P = as.numeric(auto_arima_model[["P"]])
		max.D = as.numeric(auto_arima_model[["D"]])
		max.Q = as.numeric(auto_arima_model[["Q"]])
		max.order = max.p + max.q + max.P + max.Q
	}


	cases.train.mod = all_years_epi$model$raw[istart:(my_row)]

	cat("\n SARIMA Fitting to:", mydata$model$name, " Data\n\n")
	##
	## Now the covar if the User asked for it:
	tiny = 1e-06
	if (!is.null(arima_model) & covar != FALSE) {
		cat("\n Using ", covar, "For SARIMA Model With a Lag of ", covar_lag,"\n\n")
		ifirst = istart - covar_lag
		nadd = NULL
		if (ifirst <= 0) nadd = abs(ifirst) + 1
		ifirst = max(1,ifirst)
		# Lagged predictors. Test 0, 1 or 2 lags.
		if (tolower(covar) == "sh") {
			y <- all_years_epi$model$sh[(ifirst):(my_row - covar_lag)]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
				z <- mydata$model$sh[(nfit+1 - covar_lag):(nperiods - covar_lag)]
			} else {
				z= NULL
			}
		} else if (tolower(covar) == "precip") {

			y <- all_years_epi$model$precip[(ifirst):(my_row - covar_lag)]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
				z <- mydata$model$precip[(nfit+1 - covar_lag):(nperiods - covar_lag)]
			} else {
				z= NULL
			}
		} else if (tolower(covar) == "temp") {

			y <- all_years_epi$model$temp[(ifirst):(my_row - covar_lag)]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
				z <- mydata$model$temp[(nfit+1 - covar_lag):(nperiods - covar_lag)]
			} else {
				z= NULL
			}
		} else if (tolower(covar) == "school") {

			y <- all_years_epi$model$school[(ifirst):(my_row - covar_lag)]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
				z <- mydata$model$school[(nfit+1 - covar_lag):(nperiods - covar_lag)]
			} else {
				z= NULL
			}
		} else {
			y <- all_years_epi$model$sh[(ifirst):(my_row - covar_lag)]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
				z <- mydata$model$sh[(nfit+1 - covar_lag):(nperiods - covar_lag)]
			} else {
				z= NULL
			}
		}

		if (!is.null(nadd)) y = c(rep(0, nadd),y)
		y <- na.interp(y)
		y <- log(y + tiny)
		mean.y = mean(y)
		y = y - mean.y
		sd.y = sd(y)
		y = y/sd.y
		if (!is.null(z)) {
			z <- na.interp(z)
			z = log(z + tiny)
			z = z - mean(z)
			if (length(z) > 1) 
				z = z/sd(z)
		}

	}
	## Train the model on the mod_level data
	tiny = 0.001
	x = cases.train.mod
	x <- na.interp(x)
	x = log(x + tiny)
	mean.x = mean(x)
	x = x - mean.x
	sd.x = sd(x)
	x = x/sd.x

	if (is.null(arima_model)) {
		ts.x <- ts(x, frequency = my.freq)
		mod.sarima <- auto.arima(ts.x, d = 1, D = 1, max.p = max.p, max.q = max.q, max.P = max.P, max.Q = max.Q, max.d = max.d, max.D = max.D)
		fit.order <- arimaorder(mod.sarima)
		p = fit.order[1]
		d = fit.order[2]
		q = fit.order[3]
		P = D = Q = 0
		if (length(fit.order > 3)) {
			P = fit.order[4]
			D = fit.order[5]
			Q = fit.order[6]
		}
		arima_model_all[[nregion1]] = list(p = p, d = d, q = q, P = P, D = D, Q = Q)
		arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")

	} else {

		if (covar != FALSE) {
			mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), xreg = y, method = "ML"))
		} else {
			mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), method = "ML"))
		}
		arima_model_all[[nregion1]] = arima_model
		arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
		if (isTRUE(class(mod.sarima) == "try-error")) {
			mod.sarima <- try(arima(x, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = my.freq), method = "ML"))
			arima_model_all[[nregion1]] = list(p = 0, d = 1, q = 1, P = 0, D = 1, Q = 1)
			arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
		}

	}
	cat(" ARIMA MODEL: ", arima_model.chosen,'\n')

	## these are the fitted results
	x.fit = fitted(mod.sarima)

	nlast = length(x.fit)

	result = result.lower = result.upper = rep(NA, nperiods)

	result[1:nfit] = exp(x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x)
	
	lower = x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x - 1.96 * sqrt(mod.sarima$sigma2)
	upper = x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x + 1.96 * sqrt(mod.sarima$sigma2)
	
	result.lower[1:nfit] = exp(lower[1:nfit])
	result.upper[1:nfit] = exp(upper[1:nfit])

	if (nfit < nperiods) {
		if (!is.null(arima_model) & covar != FALSE) {
			future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), newxreg = z, se.fit = TRUE))
		} else {
			future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), se.fit = TRUE))
		}


		if (length(future$pred) == (nperiods - nfit)) {
			## Need to convert back from log to linear

			x.pred = as.numeric(future$pred)
			result[(nfit + 1):nperiods] = exp(x.pred * sd.x + mean.x)

			# This is an approximation.  We actually know that exp(sd(log[y])) != sd(y)

			lower = x.pred * sd.x + mean.x - 1.96 * sqrt(mod.sarima$sigma2)
			upper = x.pred * sd.x + mean.x + 1.96 * sqrt(mod.sarima$sigma2)

			result.upper[(nfit + 1):nperiods] = exp(upper)
			result.lower[(nfit + 1):nperiods] = exp(lower)
			
		}
	}

	tables.mod$arima.frcst[nfit, 1:nperiods, mod_id] = result
	tables.mod$arima.lower[nfit, 1:nperiods, mod_id] = result.lower
	tables.mod$arima.upper[nfit, 1:nperiods, mod_id] = result.upper

	tables = calc.stat.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	sarima.results = list(tables.mod = tables.mod, arima_model_all = arima_model_all)

	if (mydata$mod_level == mydata$fit_level) {
		err <- plotARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, arima_model = arima_model, arima_model_all = 		arima_model_all, covar = covar, covar_lag = covar_lag, device = device)

		err <- saveARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, arima_model = arima_model,
		arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag)

		return(sarima.results)
	}

	cases.train.fit = all_years_epi$fit$raw[istart:(my_row), ]
	colnames(cases.train.fit) = colnames(all_years_epi$fit$raw)

	fit.result = fit.result.lower = fit.result.upper = array(NA, c(nperiods, nregion))

	for (iregion in 1:nregion) {

		cat("\n SARIMA Fitting to: ", mydata$fit$name[iregion], " Data\n\n")
		x = cases.train.fit[, iregion]
		if(any(is.na(x))) x <- na.interp(x)
		x = log(x + tiny)
		mean.x = mean(x)
		x = x - mean.x
		sd.x = sd(x)
		x = x/sd.x

		if (!is.null(arima_model) & covar != FALSE) {

			ifirst = istart - covar_lag
			nadd = NULL
			if (ifirst <= 0) nadd = abs(ifirst) + 1
			ifirst = max(1,ifirst)
			# Lagged predictor is done here
			if (tolower(covar) == "sh") {
				y <- all_years_epi$fit$sh[(ifirst):(my_row - covar_lag), iregion]
				if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
					z <- mydata$fit$sh[(nfit+1 - covar_lag):(nperiods - covar_lag), iregion]
				} else {
					z = NULL
				}
			} else if (tolower(covar) == "precip") {
				y <- all_years_epi$fit$precip[(ifirst):(my_row - covar_lag), iregion]
				if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
					z <- mydata$fit$precip[(nfit+1 - covar_lag):(nperiods - covar_lag), iregion]
				} else {
					z = NULL
				}
			} else if (tolower(covar) == "temp") {
				y <- all_years_epi$fit$temp[(ifirst):(my_row - covar_lag), iregion]
				if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
					z <- mydata$fit$temp[(nfit+1 - covar_lag):(nperiods - covar_lag), iregion]
				} else {
					z = NULL
				}
			} else if (tolower(covar) == "school") {
				y <- all_years_epi$fit$school[(ifirst):(my_row - covar_lag), iregion]
			if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
					z <- mydata$fit$school[(nfit+1 - covar_lag):(nperiods - covar_lag), iregion]
				} else {
					z = NULL
				}
			} else {
				y <- all_years_epi$fit$sh[(ifirst):(my_row - covar_lag), iregion]
				if((nfit+1 - covar_lag) <= (nperiods - covar_lag)) {
					z <- mydata$fit$sh[(nfit+1 - covar_lag):(nperiods - covar_lag), iregion]
				} else {
					z = NULL
				}
			}
			if (!is.null(nadd)) y = c(rep(0, nadd),y)
			y <- na.interp(y)
			y <- log(y + tiny)
			mean.y = mean(y)
			y = y - mean.y
			sd.y = sd(y)
			y = y/sd.y
			if (!is.null(z)) {
				z <- na.interp(z)
				z = log(z + tiny)
				z = z - mean(z)
				if (length(z) > 1) 
					z = z/sd(z)

			}
		}

		if (is.null(arima_model)) {
			ts.x <- ts(x, frequency = my.freq)
			mod.sarima <- auto.arima(ts.x, d = 1, D = 1, max.p = max.p, max.q = max.q, max.P = max.P, max.Q = max.Q, max.d = max.d, max.D = max.D)
			fit.order <- arimaorder(mod.sarima)
			p = fit.order[1]
			d = fit.order[2]
			q = fit.order[3]
			P = D = Q = 0
			if (length(fit.order > 3)) {
				P = fit.order[4]
				D = fit.order[5]
				Q = fit.order[6]
			}
			arima_model_all[[iregion]] = list(p = p, d = d, q = q, P = P, D = D, Q = Q)
			arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
		} else {
			if (covar != FALSE) {
				mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), xreg = y, method = "ML"))
			} else {
				mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), method = "ML"))
			}
			arima_model_all[[iregion]] = arima_model
			arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
			if (isTRUE(class(mod.sarima) == "try-error")) {
				mod.sarima <- try(arima(x, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = my.freq), method = "ML"))
				arima_model_all[[iregion]] = list(p = 0, d = 1, q = 1, P = 0, D = 1, Q = 1)
				arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
			}

		}

		cat("\n ARIMA MODEL: ", arima_model.chosen,'\n')

		## these are the fitted results
		x.fit = fitted(mod.sarima)

		nlast = length(x.fit)

		fitted <- x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x
		fit.result[1:nfit, iregion] = exp(fitted)

		lower = x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x - 1.96 * sqrt(mod.sarima$sigma2)
		upper = x.fit[(nlast - nfit + 1):nlast] * sd.x + mean.x + 1.96 * sqrt(mod.sarima$sigma2)

		fit.result.lower[1:nfit, iregion] = exp(lower[1:nfit])
		fit.result.upper[1:nfit, iregion] = exp(upper[1:nfit])

		if (nfit < nperiods) {
			if (!is.null(arima_model) & covar != FALSE) {
				future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), newxreg = z, se.fit = TRUE))
			} else {
				future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), se.fit = TRUE))
			}

			if (length(future$pred) == (nperiods - nfit)) {
				## Need to convert back from log to linear

				x.pred = as.numeric(future$pred)
				fit.result[(nfit + 1):nperiods, iregion] = exp(x.pred * sd.x + mean.x)

				# This is an approximation.  We actually know that exp(sd(log[y])) != sd(y)

				lower = x.pred * sd.x + mean.x - 1.96 * sqrt(mod.sarima$sigma2)
				upper = x.pred * sd.x + mean.x + 1.96 * sqrt(mod.sarima$sigma2)

				fit.result.upper[(nfit + 1):nperiods, iregion] = exp(upper)
				fit.result.lower[(nfit + 1):nperiods, iregion] = exp(lower)

			}
		}

	}


	tables.fit$arima.frcst[nfit, 1:nperiods, 1:nregion] = fit.result
	tables.fit$arima.lower[nfit, 1:nperiods, 1:nregion] = fit.result.lower
	tables.fit$arima.upper[nfit, 1:nperiods, 1:nregion] = fit.result.upper

	## Build the aggregate

	tables.agg.mod = tables.mod
	for (iperiod in 1:nperiods) {
		tables.agg.mod$arima.frcst[nfit, iperiod, mod_id] = sum(tables.fit$arima.frcst[nfit, iperiod, 1:nregion] * mydata$fit$coef[1:nregion])
		tables.agg.mod$arima.lower[nfit, iperiod, mod_id] = sum(tables.fit$arima.lower[nfit, iperiod, 1:nregion] * mydata$fit$coef[1:nregion])
		tables.agg.mod$arima.upper[nfit, iperiod, mod_id] = sum(tables.fit$arima.upper[nfit, iperiod, 1:nregion] * mydata$fit$coef[1:nregion])

	}


	tables = calc.stat.err(tables = tables.agg.mod, nfit = nfit, state_id = mod_id)
	tables.agg.mod = tables

	tables = calc.stat.err(tables = tables.fit, nfit = nfit, state_id = fit_id)
	tables.fit = tables


	err <- plotARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, arima_model = arima_model,
		arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag, device = device)

	sarima.results = list(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, arima_model = arima_model,
		arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag)

	## Save and RData file with all the results

	err <- saveARIMA(all_years_epi = all_years_epi, mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, arima_model = arima_model,
		arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag)


	return(sarima.results)
}

fitSARIMASD <- function(mydata = NULL, all_years_epi = NULL, arima_model = NULL, auto_arima_model = NULL, device = 'png') {
	#'
	#' SARIMA Fitting to San Diego ILI Data
	#'
	#'\code{fitSARIMASD} fits the model data with a SARIMA model (either chosen by the user or using \pkg{auto.arima})
	#' and calls the plotting function
	#' @param mydata - dataframe with all the data for this \pkg{DICE} run
	#' @param all_years_epi the epi data for SD county all years
	#' @param arima_model List of ARIMA model parameters: list(p=, d=, q=, P=, D, Q=) can be set to NULL to trigger the
	#' @param device - 'pdf' or 'png'. Default is 'png'
	#' \code{auto.arima} process
	#' @param auto_arima_model A list of upper limit values for the ARIMA parameters.  This list is created by the
	#' \code{DICE} code ONLY if arima_model is not set by the user
	#' @examples
	#' fitSARIMASD((mydata = mydata, all_years_epi = all_years_epi, arima_model = arima_model,
	#' auto_arima_model = NULL, covar = covar, covar_lag = covar_lag, run.list = run.list)
	#' @return
	#' A list of tables with the results for the model region
	#'


	## Remove the pandemic year from all_years_epi - it messes things for the training!
	## 
	
	istart =  which(all_years_epi$years == 2010 & all_years_epi$weeks == 27)
	iend = length(all_years_epi$years)
	all_years_epi$years = all_years_epi$years[istart:iend]
	all_years_epi$weeks = all_years_epi$weeks[istart:iend]
	all_years_epi$dates = all_years_epi$dates[istart:iend]
	all_years_epi$model = all_years_epi$model[istart:iend,]
	all_years_epi$fit$raw    = all_years_epi$fit$raw[   istart:iend]
	all_years_epi$fit$epi    = all_years_epi$fit$epi[   istart:iend]
	all_years_epi$fit$school = all_years_epi$fit$school[istart:iend]
	all_years_epi$fit$sh     = all_years_epi$fit$sh[    istart:iend]
	all_years_epi$fit$temp   = all_years_epi$fit$temp[  istart:iend]
	all_years_epi$fit$precip = all_years_epi$fit$precip[istart:iend]
	
	mod_level = mydata$mod_level
	fit_level = mydata$fit_level
	
	cadence = mydata$cadence

	nfit = mydata$nperiodsFit

	mod_id = mydata$model$attr[[paste0("ABBV_", mod_level)]]
	fit_id = mydata$fit$attr[[paste0("ABBV_", fit_level)]]

	mod_name = mydata$model$name
	fit_name = mydata$fit$name
	nregion = mydata$fit$nregions
	nregion1 = nregion + 1
	nregion2 = nregionnregion1 + 1

	covar = mydata$covar
	covar_lag = mydata$covar_lag

	if (mydata$cadence == "Weekly") {
		cadence = "week"
		nperiods = mydata$nperiods
	} else if (mydata$cadence == "Monthly") {
		cadence = "month"
		nperiods = mydata$nperiods
	} else if (mydata$cadence == "Daily") {
		cadence = "day"
		nperiods = mydata$nperiods
	} else {
		cadence = "Unknown"
		cat("\n\n Cadence of Data is Unknown, Code will Stop !!\n\n")
		q()
	}

	## Prepare a list for the arima models - really needed only when arima_model = NULL and each region can have a  different
	## arima model coming out of auto.arima

	arima_model_all <- list()

	for (iregion in 1:nregion2) {
		arima_model_all[[iregion]] <- arima_model
	}

	tables.mod <- results.table(state_names = mod_name, state_id = mod_id, cadence.names = mydata$cadence.names, cadence = cadence, nperiods = nperiods)

	epi.obsrv = mydata$model$raw

	tables.mod$epi.obsrv[, mod_id] = as.matrix(epi.obsrv)

	epi.null = calc.epi.null(mydata = mydata, all_years_epi = all_years_epi, mod_id = mod_id, fit_id = fit_id)

	epi.null.mod = epi.null$cases.ave.mod

	tables.mod$epi.null[nfit, , mod_id] = epi.null.mod[, mod_id]


	tables = calc.null.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	if (mydata$cadence == "Monthly") {
		my_row = which(all_years_epi$years == mydata$years[nfit] & all_years_epi$months == mydata$months[nfit])
		cadence.vec = all_years_epi$month
		my.freq = 12
	}

	if (mydata$cadence == "Weekly") {
		my_row = which(all_years_epi$years == mydata$years[nfit] & all_years_epi$weeks == mydata$weeks[nfit])
		cadence.vec = all_years_epi$weeks
		my.freq = 52
	}

	years.vec = all_years_epi$years
	year.vec = mydata$years

	attr = c("year", cadence, "ndays")

	istart = my_row - my.freq * 10
	istart = max(istart, 1)
	ntrn = length(istart:my_row)

	if (!is.null(arima_model)) {
		p = as.numeric(arima_model[["p"]])
		d = as.numeric(arima_model[["d"]])
		q = as.numeric(arima_model[["q"]])

		P = as.numeric(arima_model[["P"]])
		D = as.numeric(arima_model[["D"]])
		Q = as.numeric(arima_model[["Q"]])
	} else {
		max.p = as.numeric(auto_arima_model[["p"]])
		max.d = as.numeric(auto_arima_model[["d"]])
		max.q = as.numeric(auto_arima_model[["q"]])

		max.P = as.numeric(auto_arima_model[["P"]])
		max.D = as.numeric(auto_arima_model[["D"]])
		max.Q = as.numeric(auto_arima_model[["Q"]])
		max.order = max.p + max.q + max.P + max.Q
	}


	cases.train.mod = all_years_epi$model$raw[istart:(my_row)]

	cat("\n SARIMA Fitting to:", mydata$model$name, " Data\n\n")
	##
	## Now the covar if the User asked for it:
	tiny = 1e-06
	if (!is.null(arima_model) & covar != FALSE) {
		cat("\n Using ", covar, "For SARIMA Model With a Lag of ", covar_lag,"\n\n")
		ifirst = istart - covar_lag
		nadd = NULL
		if (ifirst <= 0) nadd = abs(ifirst) + 1
		ifirst = max(1,ifirst)
		# Lagged predictors. Test 0, 1 or 2 lags.
		if (tolower(covar) == "sh") {
			y <- all_years_epi$model$sh[(ifirst):(my_row - covar_lag)]
			z <- mydata$model$sh[(nfit + 1 - covar_lag):(nperiods - covar_lag)]
		} else if (tolower(covar) == "precip") {

			y <- all_years_epi$model$precip[(ifirst):(my_row - covar_lag)]
			z <- mydata$model$precip[(nfit + 1 - covar_lag):(nperiods - covar_lag)]
		} else if (tolower(covar) == "temp") {

			y <- all_years_epi$model$temp[(ifirst):(my_row - covar_lag)]
			z <- mydata$model$temp[(nfit + 1 - covar_lag):(nperiods - covar_lag)]
		} else if (tolower(covar) == "school") {

			y <- all_years_epi$model$school[(ifirst):(my_row - covar_lag)]
			z <- mydata$model$school[(nfit + 1 - covar_lag):(nperiods - covar_lag)]
		} else {
			y <- all_years_epi$model$sh[(ifirst):(my_row - covar_lag)]
			z <- mydata$model$sh[(nfit+1 - covar_lag):(nperiods - covar_lag)]
		}

		if (!is.null(nadd)) y = c(rep(0, nadd),y)
		y <- na.interp(y)
		y <- log(y + tiny)
		mean.y = mean(y)
		y = y - mean.y
		sd.y = sd(y)
		y = y/sd.y
		z <- na.interp(z)
		z = log(z + tiny)
		z = z - mean(z)
		if(length(z) > 1) z = z/sd(z)
	}
	## Train the model on the mod_level data
	tiny = 0.001
	x = cases.train.mod
	x <- na.interp(x)

	if (is.null(arima_model)) {
		ts.x <- ts(x, frequency = my.freq)
		mod.sarima <- auto.arima(ts.x, d = 1, D = 1, max.p = max.p, max.q = max.q, max.P = max.P, max.Q = max.Q, max.d = max.d, max.D = max.D)
		fit.order <- arimaorder(mod.sarima)
		p = fit.order[1]
		d = fit.order[2]
		q = fit.order[3]
		P = D = Q = 0
		if (length(fit.order > 3)) {
			P = fit.order[4]
			D = fit.order[5]
			Q = fit.order[6]
		}
		arima_model_all[[nregion1]] = list(p = p, d = d, q = q, P = P, D = D, Q = Q)
		arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")

	} else {

		if (covar != FALSE) {
			mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), xreg = y, method = "ML"))
		} else {
			mod.sarima <- try(arima(x, order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = my.freq), method = "ML"))
		}
		arima_model_all[[nregion1]] = arima_model
		arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
		if (isTRUE(class(mod.sarima) == "try-error")) {
			mod.sarima <- try(arima(x, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = my.freq), method = "ML"))
			arima_model_all[[nregion1]] = list(p = 0, d = 1, q = 1, P = 0, D = 1, Q = 1)
			arima_model.chosen = paste0("(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")")
		}

	}
	cat(" ARIMA MODEL: ", arima_model.chosen,'\n')

	## these are the fitted results
	x.fit = fitted(mod.sarima)

	nlast = length(x.fit)

	result = result.lower = result.upper = rep(NA, nperiods)
	result[1:nfit] = x.fit[(nlast - nfit + 1):nlast]
	
	lower = x.fit[(nlast - nfit + 1):nlast] - 1.96 * sqrt(mod.sarima$sigma2)
	upper = x.fit[(nlast - nfit + 1):nlast] + 1.96 * sqrt(mod.sarima$sigma2)
	
	result.lower[1:nfit] = lower[1:nfit]
	result.upper[1:nfit] = upper[1:nfit]

	if (nfit < nperiods) {
		if (!is.null(arima_model) & covar != FALSE) {
			future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), newxreg = z, se.fit = TRUE))
		} else {
			future <- try(predict(mod.sarima, n.ahead = (nperiods - nfit), se.fit = TRUE))
		}


		if (length(future$pred) == (nperiods - nfit)) {

			x.pred = as.numeric(future$pred)
			result[(nfit + 1):nperiods] = x.pred 

			# This is an approximation.  We actually know that exp(sd(log[y])) != sd(y)

			lower = x.pred - 1.96 * sqrt(mod.sarima$sigma2)
			upper = x.pred + 1.96 * sqrt(mod.sarima$sigma2)

			result.upper[(nfit + 1):nperiods] = upper
			result.lower[(nfit + 1):nperiods] = lower
			

		}
	}

	tables.mod$arima.frcst[nfit, 1:nperiods, mod_id] = result
	tables.mod$arima.lower[nfit, 1:nperiods, mod_id] = result.lower
	tables.mod$arima.upper[nfit, 1:nperiods, mod_id] = result.upper

	tables = calc.stat.err(tables = tables.mod, nfit = nfit, state_id = mod_id)
	tables.mod = tables

	sarima.results = list(tables.mod = tables.mod, arima_model_all = arima_model_all)

	err <- plotARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, 
		arima_model = arima_model, arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag, 
		device = device)

	err <- saveARIMA(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL, tables.agg.mod = NULL, 
		arima_model = arima_model, arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag)

	return(sarima.results)
}


saveARIMA <- function(all_years_epi = NULL, mydata = NULL, tables.mod = NULL, tables.fit = NULL, tables.agg.mod = NULL, arima_model = NULL,
	arima_model_all = NULL, covar = NULL, covar_lag = NULL) {

	#' Save the Results of a SARIMA Fit
	#'
	#' \code{saveARIMA} saves the results of a SARIMA fit to the data into an RData file
	#' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param tables.mod - a table with the data and the results for the direct Model fit
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.fit - Optional a table with the data and the results for the fits at the fit_level. It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.agg.mod - A table with aggregate results of the fits
    #' @param arima_model - A list with the user selection of the ARIMA model
    #' If the user has not selected a model this will be NULL and arima_model_all
    #' will have the models selected by auto.arima
    #' @param arima_model_all - An array with details of the ARIMA model used for
    #' model region, the fits regions and the aggregate (the last is only if the User
    #' has chosen an ARIMA model for the run)
    #' @param covar - String. Optional covariate variable for ARIMA fit.
    #' 'sh', 'precip' and 'temp' are currently supported. Default is NULL - no covariate
    #' @param covar_lag - Numeric. Lag time for covariate variable in units of the cadence
    #' of the data
	#' @examples
    #' saveARIMA(all_years_epi = all_years_epi, mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit,
    #' tables.agg.mod = tables.agg.mod,arima_model = arima_model,
    #' arima_model_all = NULL, covar = covar, covar_lag = covar_lag)
	#' @return
	#' A list of tables with the results for the model (direct and aggregate) and fit regions
	#'

	subDir = mydata$subDir
	myName = mydata$dataName
	nperiodsFit = mydata$nperiodsFit
	myName = gsub(" ","",myName)
	
	if (!is.null(arima_model)) {
		p = as.numeric(arima_model[["p"]])
		d = as.numeric(arima_model[["d"]])
		q = as.numeric(arima_model[["q"]])

		P = as.numeric(arima_model[["P"]])
		D = as.numeric(arima_model[["D"]])
		Q = as.numeric(arima_model[["Q"]])

		epi.model = paste0(p, d, q, "-", P, D, Q)

		if (covar != FALSE) {
			epi.model = paste0(epi.model, "-", covar, "-lag-", covar_lag)
		}
	} else {

		epi.model = "auto-arima"
	}

	sarima.results = list(all_years_epi = all_years_epi, mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, arima_model = arima_model,
		arima_model_all = arima_model_all, covar = covar, covar_lag = covar_lag)

	filename = paste0(subDir, "/sarima-", myName, "-", epi.model, "-", nperiodsFit, ".RData", sep = "")

	save(sarima.results, file = filename)
	save.image()

	cat("\n Writing R object Data file for this Calculation:",filename,'\n\n')

	return(err = 0)
}


