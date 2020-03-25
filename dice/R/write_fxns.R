##
## Output writing functions routines.  Not all of them are here.
## In some cases the plotting routines (plot_fxns.R) also do the writing
##

writecsvOnePatch <- function(mydata = NULL, run.list = NULL, tab = NULL, model_rtn = NULL, model_profile = NULL, ireal = 1) {

    #' Write the Results of a  \pkg{DICE} Run - Single Region
    #'
    #' Write the results of a  \pkg{DICE} run for a single region/patch to csv files. We write the observed %ILI (or number of cases) along with our fits and
    #' if appropriate predictions. We write the best result along with the mean
    #' @param tab An array with the history of the MCMC chain
    #' @param model_rtn A 1D numeric array with the best direct prediction to the region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the mydata
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' writecsvOnePatch{mydata=mydata,run.list = run.list, tab = tab.model, model_rtn = model_rtn,
    #'                  model_profile = model_profile, ireal = ireal}


    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName
	myName = gsub(" ","",myName)
	
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
   	tps = mydata$weeks
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

    # write a file with week number, date, mydata, model - mean, model best, model sd

    model_ili_sd = rep(0, nperiods)
    model_ili_mean = rep(0, nperiods)
    for (i in 1:nperiods) {
        model_ili_mean[i] = mean(model_profile_ili[, i])
        model_ili_sd[i] = sd(model_profile_ili[, i])
    }
    weeks_fitted = rep(0, nperiods)
    weeks_fitted[1:nperiodsFit] = 1

    table_ili = data.frame(EW = weeks, YEARS = mydata$years, WKENDAT = mydata$dates, WK_FITTED = weeks_fitted, DATA = mydata$model$raw,
        MODEL_ILI_BEST = model_rtn_ili, MODEL_ILI_MEAN = model_ili_mean, MODEL_ILI_SD = model_ili_sd)

    filename = paste(subDir, "/ili-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(table_ili, file = filename)
    cat("\n For a CSV Table with the ILI Data and Model Results See: ", filename, "\n")

    # statistics about onset, peak week, peak value, percent clinical and the force of infection

    # Need to add the onset here - three consequative weeks of values >= onset value
    iline = 0
    if (length(model_onset) > 0 && !is.na(model_onset)) {

        onset = FALSE
        onset_mydata = NA

        for (i in 1:(nperiodsData - 2)) {
            if (onset == FALSE) {
                val_i1 = mydata$model$raw[i]
                val_i2 = mydata$model$raw[(i + 1)]
                val_i3 = mydata$model$raw[(i + 2)]
                if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                  onset = TRUE
                  onset_mydata = weeks[i]
                }
            }
        }


        onset = FALSE
        onset_best = NA
        for (i in 1:(nperiods - 2)) {
            if (onset == FALSE) {
                val_i1 = model_rtn_ili[i]
                val_i2 = model_rtn_ili[(i + 1)]
                val_i3 = model_rtn_ili[(i + 2)]
                if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                  onset = TRUE
                  onset_best = weeks[i]
                }
            }
        }

        # now do the statistics on the onset
        onset_vec = rep(NA, nRnd)
        for (irnd in 1:nRnd) {
            onset = FALSE
            for (i in 1:(nperiods - 2)) {
                if (onset == FALSE) {
                  val_i1 = model_profile_ili[irnd, i]
                  val_i2 = model_profile_ili[irnd, (i + 1)]
                  val_i3 = model_profile_ili[irnd, (i + 2)]
                  if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                    onset = TRUE
                    onset_vec[irnd] = weeks[i]
                  }
                }
            }
        }
        onset_mean = mean(onset_vec, na.rm = TRUE)
        onset_sd = sd(onset_vec, na.rm = TRUE)
        iline = 1
    } else {
        onset_mydata = NA
        onset_best = NA
        onset_mean = NA
        onset_sd = NA
    }


    if (any(model == c(4, 101, 103, 105))) {
        table_results = array(0, c(4 + iline, 4))

    } else if (model == 3 || model == 5) {
        table_results = array(0, c(6 + iline, 4))
    } else if (any(model == c(1, 2, 102, 104))) {
        table_results = array(0, c(5 + iline, 4))
    }



    mydata_pk_val = max(mydata$model$raw, na.rm = TRUE)
    mydata_pk_wk = weeks[which.max(model_rtn_ili)]
    model_pk_val_best = max(model_rtn_ili, na.rm = TRUE)
    model_pk_wk_best = weeks[which.max(model_rtn_ili)]

    model_max = rep(0, nRnd)
    model_wk = rep(0, nRnd)
    for (i in 1:nRnd) {
        model_max[i] = max(model_profile_ili[i, ], na.rm = TRUE)
        model_wk[i] = weeks[which.max(model_profile_ili[i, ])]
    }
    model_pk_val_mean = mean(model_max)
    model_pk_val_sd = sd(model_max)
    model_pk_wk_mean = mean(model_wk)
    model_pk_wk_sd = sd(model_wk)

    # percent clinical using second half of teh chain
    nlines = dim(tab)[1]
    nparam1 = dim(tab)[2]
    nlines2 = nlines/2
    index = 2

    pC_mean = mean(tab[nlines2:nlines, index])
    pC_sd = sd(tab[nlines2:nlines, index])
    imin = which.min(tab[nlines2:nlines, nparam1])
    pC_best = tab[imin, index]

    # for R0

    index = 3
    R0_mean = mean(tab[nlines2:nlines, index])
    R0_sd = sd(tab[nlines2:nlines, index])
    R0_best = tab[imin, index]



    if (any(model == c(1, 3, 102))) {
        index1 = 4
        index2 = 5
        sh = mydata$model$sh
        tmp = rep(0, length(tab[nlines2:nlines, 1]))
        for (i in nlines2:nlines) tmp[i] = mean(tab[i, index1] * exp(-tab[i, index2] * sh))
        del2_mean = mean(tmp)
        del2_sd = sd(tmp)
        del2_best = mean(tab[imin, index1] * exp(-tab[imin, index2] * sh))
    }

    if (model == 1 || model == 102) {
        if (length(model_onset) > 0 && !is.na(model_onset))
            rownames(table_results) = c("Onset_Wk", "Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "SH_Term") else rownames(table_results) = c("Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "SH_Term")

        table_results[5 + iline, ] = c(NA, del2_best, del2_mean, del2_sd)

    } else if (model == 2 || model == 104) {

        if (length(model_onset) > 0 && !is.na(model_onset))
            rownames(table_results) = c("Onset_Wk", "Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "School_Term") else rownames(table_results) = c("Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "School_Term")
        index = 6
        del_best = tab[imin, index]
        del_mean = mean(tab[nlines2:nlines, index])
        del_sd = sd(tab[nlines2:nlines, index])
        table_results[5 + iline, ] = c(NA, del_best, del_mean, del_sd)
    } else if (model == 3) {


        if (length(model_onset) > 0 && !is.na(model_onset))
            rownames(table_results) = c("Onset_Wk", "Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "SH_Term", "School_Term") else rownames(table_results) = c("Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "SH_Term", "School_Term")
        index = 6
        del_best = tab[imin, index]
        del_mean = mean(tab[nlines2:nlines, index])
        del_sd = sd(tab[nlines2:nlines, index])
        table_results[5 + iline, ] = c(NA, del2_best, del2_mean, del2_sd)
        table_results[6 + iline, ] = c(NA, del_best, del_mean, del_sd)
    } else if (model == 5) {


        if (length(model_onset) > 0 && !is.na(model_onset))
            rownames(table_results) = c("Onset_Wk", "Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "delta_R0", "Duration") else rownames(table_results) = c("Peak_EW", "Peak_Value", "Percent_Clinical", "R0", "delta_R0", "Duration")
        index = 9
        del_best = tab[imin, index]
        del_mean = mean(tab[nlines2:nlines, index])
        del_sd = sd(tab[nlines2:nlines, index])
        table_results[5, ] = c(NA, del_best, del_mean, del_sd)
        index = 11
        # duration - converted to week
        del_best = tab[imin, index]/7
        del_mean = mean(tab[nlines2:nlines, index])/7
        del_sd = sd(tab[nlines2:nlines, index]/7)
        table_results[6 + iline, ] = c(NA, del_best, del_mean, del_sd)

    } else if (any(model == c(4, 101, 103, 105))) {

        if (length(model_onset) > 0 && !is.na(model_onset))
            rownames(table_results) = c("Onset_Wk", "Peak_EW", "Peak_Value", "Percent_Clinical", "R0") else rownames(table_results) = c("Peak_EW", "Peak_Value", "Percent_Clinical", "R0")
    }

    colnames(table_results) = c("DATA", "MODEL_BEST", "MODEL_MEAN", "MODEL_SD")

    if (length(model_onset) > 0 && !is.na(model_onset)) {
        table_results[iline, ] = c(onset_mydata, onset_best, onset_mean, onset_sd)
    }
    table_results[1 + iline, ] = c(mydata_pk_wk, model_pk_wk_best, model_pk_wk_mean, model_pk_wk_sd)
    table_results[2 + iline, ] = c(mydata_pk_val, model_pk_val_best, model_pk_val_mean, model_pk_val_sd)
    table_results[3 + iline, ] = c(NA, pC_best, pC_mean, pC_sd)
    table_results[4 + iline, ] = c(NA, R0_best, R0_mean, R0_sd)

    filename = paste(subDir, "/stats-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(table_results, file = filename)
    cat("\n For a CSV Table with  Onset Week, Week Maximum, Peak value ...: ", filename, "\n")

    cat("\n Summary of Results \n")
    print(table_results)
    err = 0
    return(err)
}

writeCSV <- function(mydata = NULL, run.list = NULL, model_rtn = NULL, model_profile = NULL, rtn = NULL,
    profile = NULL, ireal = 1) {

    #' Write the Results of a  \pkg{DICE} Run - Uncoupled Multiple Region
    #'
    #' Write the results of a \pkg{DICE} run for multiple uncoupled regions into csv files. We write the observed %ILI (or number of cases) along with our fits and
    #' if appropriate predictions. We write the best result along with the mean
    #' @param model_rtn A 1D numeric array with the best direct prediction to the region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the mydata
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' writecsvOnePatch{mydata = mydata,run.list = run.list, model_rtn = model_rtn,
    #'                  model_profile=model_profile, ireal= ireal}


    ## check to see if output directory exists and if not create it
    subDir = run.list$subDir
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    myName = mydata$dataName
	myName = gsub(" ","",myName)
    FY = mydata$FY
    model = mydata$imodel
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit
    nperiodsData = mydata$nperiodsData
    reg.model.name = mydata$model$name
    nRnd = dim(model_profile)[1]
    n.model = 1
    n.fit = mydata$fit$nregions
    model_onset = mydata$model$onset
    fit_onset = as.numeric(mydata$fit$onset)
    fit_coef = mydata$fit$coef

   if (mydata$cadence == "Weekly") {
   	cadence = 7
   	tps = mydata$weeks
   	weeks = mydata$weeks
   } else if (mydata$cadence == "Monthly") {
   	cadence = 31
   	tps = mydata$months
   } else if (mydata$cadence == "Daily") {
   	cadence = 1
   	tps = mydata$days
   } else {
   	tps = mydata$weeks
   }


    # convert the fit and model results to %ILI calculate the national

    # convert the model to %ILI calculate the national

    model_factor = mydata$model$factor
    fit_factor = as.numeric(mydata$fit$factor)

    ## These are the model and fit raw mydata

    model_ili = mydata$model$raw
    fit_ili = mydata$fit$raw

    model_rtn_ili = NULL
    model_profile_ili = NULL

    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile/model_factor
    }

    # convert the rtn and profile to %ILI
    rtn_ili = rtn
    profile_ili = profile

    for (i in 1:n.fit) {
        rtn_ili[, i] = rtn[, i]/fit_factor[i]
        profile_ili[, , i] = profile[, , i]/fit_factor[i]
    }


    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_sd = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {
        fit_model[i] = sum(rtn_ili[i, 1:n.fit] * fit_coef[1:n.fit])
        tmp = rep(0, n.fit)
        for (k in 1:n.fit) tmp[k] = mean(profile_ili[, i, k])
        fit_model_mean[i] = sum(tmp[1:n.fit] * fit_coef[1:n.fit])

        for (irnd in 1:nRnd) {
            for (k in 1:n.fit) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile_ili[irnd, i, k] * fit_coef[k]
            fit_model_sd[i] = sd(fit_model_profile[, i])
        }

    }

    # write a file with week number, date, mydata, model - mean, model best, model sd

    model_ili_sd = rep(0, nperiods)
    model_ili_mean = rep(0, nperiods)
    fit_ili_mean = array(0, c(nperiods, n.fit))
    fit_ili_sd = array(0, c(nperiods, n.fit))
    for (i in 1:nperiods) {
        model_ili_mean[i] = mean(model_profile_ili[, i])
        model_ili_sd[i] = sd(model_profile_ili[, i])
        for (j in 1:n.fit) {
            fit_ili_mean[i, j] = mean(profile_ili[, i, j])
            fit_ili_sd[i, j] = sd(profile_ili[, i, j])
        }
    }
    weeks_fitted = rep(0, nperiods)
    weeks_fitted[1:nperiodsFit] = 1

    table_ili = data.frame(EW = weeks, YEARS = mydata$years, WKENDAT = mydata$dates, WK_FITTED = weeks_fitted, DATA = mydata$model$raw,
        DIRECT_ILI_BEST = model_rtn_ili, DIRECT_ILI_MEAN = model_ili_mean, DIRECT_ILI_SD = model_ili_sd, INDIRECT_ILI_BEST = fit_model,
        INDIRECT_ILI_MEAN = fit_model_mean, INDIRECT_ILI_SD = fit_model_sd, DATA_fit = mydata$fit$raw, FIT_ILI_BEST = rtn_ili, FIT_DATA_MEAN = fit_ili_mean,
        FIT_DATA_SD = fit_ili_sd)

    filename = paste(subDir, "/ili-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(table_ili, file = filename)
    cat("\n For a CSV Table with the ILI Data and Model Results See: ", filename, "\n")

    # # # statistics about onset, peak week, peak value, percent clinical and the force of infection


    table_row_names = c("ONSET_MODEL_DIRECT", "ONSET_MODEL_INDIRECT", paste("ONSET_REGION", 1:n.fit, sep = ""), "PEAK_VAL_DIRECT", "PEAK_WEEK_DIRECT",
        "PEAK_VAL_INDIRECT", "PEAK_WEEK_INDIRECT", paste(c("PEAK_VAL_REGION", "PEAK_WK_REGION"), rep(1:n.fit, each = 2), sep = ""))


    nrow = length(table_row_names)

    table_results = array(0, c(nrow, 4))
    colnames(table_results) = c("DATA", "BEST", "MEAN", "SD")
    rownames(table_results) = table_row_names
    # Need to add the onset here - three consequative weeks of values >= onset value
    iline = 0
    if (length(model_onset) > 0 && !is.na(model_onset)) {

        onset = FALSE
        onset_mydata = NA
        for (i in 1:(nperiodsData - 2)) {
            if (onset == FALSE) {
                val_i1 = mydata$model$raw[i]
                val_i2 = mydata$model$raw[(i + 1)]
                val_i3 = mydata$model$raw[(i + 2)]
                if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                  onset = TRUE
                  onset_mydata = weeks[i]
                }
            }
        }


        onset = FALSE
        onset_best_direct = NA
        for (i in 1:(nperiods - 2)) {
            if (onset == FALSE) {
                val_i1 = model_rtn_ili[i]
                val_i2 = model_rtn_ili[(i + 1)]
                val_i3 = model_rtn_ili[(i + 2)]
                if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                  onset = TRUE
                  onset_best_direct = weeks[i]
                }
            }
        }

        # repeat for indirect
        onset = FALSE
        onset_best_indirect = NA
        for (i in 1:(nperiods - 2)) {
            if (onset == FALSE) {
                val_i1 = fit_model[i]
                val_i2 = fit_model[(i + 1)]
                val_i3 = fit_model[(i + 2)]
                if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                  onset = TRUE
                  onset_best_indirect = weeks[i]
                }
            }
        }

        # now do the statistics on the onset - direct case
        onset_vec = rep(NA, nRnd)
        for (irnd in 1:nRnd) {
            onset = FALSE
            for (i in 1:(nperiods - 2)) {
                if (onset == FALSE) {
                  val_i1 = model_profile_ili[irnd, i]
                  val_i2 = model_profile_ili[irnd, (i + 1)]
                  val_i3 = model_profile_ili[irnd, (i + 2)]
                  if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                    onset = TRUE
                    onset_vec[irnd] = weeks[i]
                  }
                }
            }
        }
        onset_mean_direct = mean(onset_vec, na.rm = TRUE)
        onset_sd_direct = sd(onset_vec, na.rm = TRUE)
        # now do the statistics on the onset - indirect case
        onset_vec = rep(NA, nRnd)
        for (irnd in 1:nRnd) {
            onset = FALSE
            for (i in 1:(nperiods - 2)) {
                if (onset == FALSE) {
                  val_i1 = fit_model_profile[irnd, i]
                  val_i2 = fit_model_profile[irnd, (i + 1)]
                  val_i3 = fit_model_profile[irnd, (i + 2)]
                  if (val_i1 >= model_onset && val_i2 >= model_onset && val_i3 >= model_onset) {
                    onset = TRUE
                    onset_vec[irnd] = weeks[i]
                  }
                }
            }
        }
        onset_mean_indirect = mean(onset_vec, na.rm = TRUE)
        onset_sd_indirect = sd(onset_vec, na.rm = TRUE)

    } else {
        onset_mydata = NA
        onset_best_direct = NA
        onset_mean_direct = NA
        onset_sd_direct = NA
        onset_best_indirect = NA
        onset_mean_indirect = NA
        onset_sd_indirect = NA
    }

    table_results[1, ] = c(onset_mydata, onset_best_direct, onset_mean_direct, onset_sd_direct)
    table_results[2, ] = c(onset_mydata, onset_best_indirect, onset_mean_indirect, onset_sd_indirect)
    onset_mydata = rep(NA, n.fit)
    onset_best_fit = rep(NA, n.fit)
    onset_mean = rep(NA, n.fit)
    onset_sd = rep(NA, n.fit)
    iline = 2

    for (j in 1:n.fit) {
        if (!is.na(fit_onset[j])) {
            onset = FALSE
            for (i in 1:(nperiodsData - 2)) {
                if (onset == FALSE) {
                  val_i1 = mydata$fit$raw[i, j]
                  val_i2 = mydata$fit$raw[(i + 1), j]
                  val_i3 = mydata$fit$raw[(i + 2), j]
                  if (!any(is.na(c(val_i1, val_i2, val_i3)))) {
                    if (val_i1 >= fit_onset[j] && val_i2 >= fit_onset[j] && val_i3 >= fit_onset[j]) {
                      onset = TRUE
                      onset_mydata[j] = weeks[i]
                    }
                  }
                }
            }  # end of loop on weeks

            onset = FALSE
            for (i in 1:(nperiods - 2)) {
                if (onset == FALSE) {
                  val_i1 = rtn_ili[i, j]
                  val_i2 = rtn_ili[(i + 1), j]
                  val_i3 = rtn_ili[(i + 2), j]
                  if (val_i1 >= fit_onset[j] && val_i2 >= fit_onset[j] && val_i3 >= fit_onset[j]) {
                    onset = TRUE
                    onset_best_fit[j] = weeks[i]
                  }
                }
            }  # end of loop on weeks

            onset_vec = rep(NA, nRnd)
            for (irnd in 1:nRnd) {
                onset = FALSE
                for (i in 1:(nperiods - 2)) {
                  if (onset == FALSE) {
                    val_i1 = profile_ili[irnd, i, j]
                    val_i2 = profile_ili[irnd, (i + 1), j]
                    val_i3 = profile_ili[irnd, (i + 2), j]
                    if (val_i1 >= fit_onset[j] && val_i2 >= fit_onset[j] && val_i3 >= fit_onset[j]) {
                      onset = TRUE
                      onset_vec[irnd] = weeks[i]
                    }
                  }
                }
            }
            onset_mean[j] = mean(onset_vec, na.rm = TRUE)
            onset_sd[j] = mean(onset_vec, na.rm = TRUE)
        }
        table_results[(iline + j), ] = c(onset_mydata[j], onset_best_fit[j], onset_mean[j], onset_sd[j])

    }  # end of loop on regions
    iline = iline + j



    # peak wk and peak value for mydata
    mydata_pk_val = max(mydata$model$raw, na.rm = TRUE)
    mydata_pk_wk = weeks[which.max(mydata$model$raw)]

    # peak week and peak value for direct model

    model_pk_val_best = max(model_rtn_ili, na.rm = TRUE)
    model_pk_wk_best = weeks[which.max(model_rtn_ili)]

    model_max = rep(0, nRnd)
    model_wk = rep(0, nRnd)
    for (i in 1:nRnd) {
        model_max[i] = max(model_profile_ili[i, ], na.rm = TRUE)
        model_wk[i] = weeks[which.max(model_profile_ili[i, ])]
    }
    model_pk_val_mean = mean(model_max)
    model_pk_val_sd = sd(model_max)
    model_pk_wk_mean = mean(model_wk)
    model_pk_wk_sd = sd(model_wk)

    table_results[(iline + 1), ] = c(mydata_pk_val, model_pk_val_best, model_pk_val_mean, model_pk_val_sd)
    table_results[(iline + 2), ] = c(mydata_pk_wk, model_pk_wk_best, model_pk_wk_mean, model_pk_wk_sd)
    # repeat for indirect modeling
    model_pk_val_best = max(fit_model, na.rm = TRUE)
    model_pk_wk_best = weeks[which.max(fit_model)]

    model_max = rep(0, nRnd)
    model_wk = rep(0, nRnd)
    for (i in 1:nRnd) {
        model_max[i] = max(fit_model_profile[i, ], na.rm = TRUE)
        model_wk[i] = weeks[which.max(fit_model_profile[i, ])]
    }
    model_pk_val_mean = mean(model_max)
    model_pk_val_sd = sd(model_max)
    model_pk_wk_mean = mean(model_wk)
    model_pk_wk_sd = sd(model_wk)

    table_results[(iline + 3), ] = c(mydata_pk_val, model_pk_val_best, model_pk_val_mean, model_pk_val_sd)
    table_results[(iline + 4), ] = c(mydata_pk_wk, model_pk_wk_best, model_pk_wk_mean, model_pk_wk_sd)
    icount = iline + 4 + 1

    # now repeart for each region
    for (j in 1:n.fit) {

        mydata_pk_val = max(mydata$fit$raw[, j], na.rm = TRUE)
        mydata_pk_wk = weeks[which.max(mydata$fit$raw[, j])]

        # peak week and peak value for model - best

        model_pk_val_best = max(rtn_ili[, j], na.rm = TRUE)
        model_pk_wk_best = weeks[which.max(rtn_ili[, j])]

        for (i in 1:nRnd) {
            model_max[i] = max(profile_ili[i, , j], na.rm = TRUE)
            model_wk[i] = weeks[which.max(profile_ili[i, , j])]
        }
        model_pk_val_mean = mean(model_max)
        model_pk_val_sd = sd(model_max)
        model_pk_wk_mean = mean(model_wk)
        model_pk_wk_sd = sd(model_wk)

        table_results[icount, ] = c(mydata_pk_val, model_pk_val_best, model_pk_val_mean, model_pk_val_sd)
        table_results[(icount + 1), ] = c(mydata_pk_wk, model_pk_wk_best, model_pk_wk_mean, model_pk_wk_sd)

        icount = icount + 2

    }

    filename = paste(subDir, "/stats-", myName, "-", nperiodsFit, "-", ireal, ".csv", sep = "")
    write.csv(table_results, file = filename)
    cat("\n For a CSV Table with  Onset Week, Week Maximum, Peak Value See: ", filename, "\n")

    err = 0
    return(err)
}

saveRData <- function(mydata = NULL, all_years_epi = NULL, run.input = NULL, ireal = 1) {

    #' Save a Binary RData file
    #'
    #' Saves a binary RData file with all the input parameters and the data for this run.
    #'   Each MCMC chain saves its own RData file and the file name includes the MCMC chain number
    #'   as the last digit before the .RData extension of the file name.
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param all_years_epi - a dataframe with the entire incidence and climate history
    #' @param run.input a list of all other input  needed for the run
    #' @param ireal Integer indicating the MCMC chain number

    #' @return err   \eqn{err =0} if successful.
    #' @examples
    #' saveRData(mydata = mydata, run.input= run.input,ireal = 1)


	subDir = run.input$run.list$subDir
	
    myName = mydata$dataName
    myName = gsub(" ","",myName)
    nperiodsFit = mydata$nperiodsFit
    # check to see if 'mydata' sub-directory exists, if not create it
    if (is.null(subDir))
        subDir = "output"
    err = makeDir(subDir = subDir)
    filename <- paste(subDir, "/input-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    cat("\n Saving Input Run Information to: ", filename, "\n")

    save(mydata, all_years_epi, run.input, file = filename)

    err = 0
    return(err)

}

saveProfilesOnePatch <- function(mydata=NULL, model_rtn = NULL, model_profile = NULL, model_profileB = NULL, ireal = 1, run.list = NULL) {
	#'
	#' Save the profiles for Model region
	#' Save the profiles for Model region so that incidence plot can be re-created
	#' 
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run	
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param model_profileB A 2D numeric array with randomly chosen predicted profiles for the Bacteria obtained by fitting the model region directly. Relevant only in the case of an SIRB model
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' saveProfilesOnePatch{ mydata = mydata, model_rtn = model_rtn, model_profile = model_profile,
    #' model_profileB = model_profileB, ireal = ireal, run.list = run.list}	
    
    myName = mydata$dataName
    myName = gsub(" ","",myName)
    nperiods = mydata$nperiods
    nperiodsFit = mydata$nperiodsFit

    model_coef = mydata$model$coef
    model_onset = mydata$model$onset
    
    model_factor = mydata$model$factor

    ## This is the model raw data

    model_ili = mydata$model$raw

    model_rtn_ili = NULL
    model_profile_ili = NULL

    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile
        model_profile_ili = model_profile_ili/model_factor
    }

    # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, model_profileB = model_profileB, model_ili = model_ili, model_rtn_ili = model_rtn_ili, model_profile_ili = model_profile_ili, run.list = run.list, ireal = ireal, idevice = 1)

    filename = pdfName = paste(mydata$subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename)    
 
	cat("\nSaving Profiles to a Binary File : ", filename, "\n")
	    
    return(err = 0)
}

saveProfiles <- function(mydata=NULL, model_rtn = NULL, model_profile = NULL, rtn = NULL, profile = NULL, ireal = 1, run.list = NULL) {
	#'
	#' Save the profiles for Model region
	#' Save the profiles for Model region so that incidence plot can be re-created
	#' 
    #' @param mydata A dataframe with all the data available for this \pkg{DICE} run	
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC chains.
    #' @param model_rtn A 1D numeric array with the best direct prediction to the model region
    #' @param model_profile A 2D numeric array with randomly chosen predicted profiles obtained by fitting the model region directly.
    #' @param rtn An nweeks x nregion 2D  numeric array with the best fit for each region
    #' @param profile A 3D numeric array holding random predictions for each of the fit regions based on the history of their MCMC 
    #' chains.
    #' @param ireal Integer - the MCMC chain number
    #' @param run.list  A list with various run parameters
    #' @return  Returns \eqn{err = 0} if successful
    #' @examples
    #' saveProfiles{ mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, rtn = rtn, 
    #' profile = profile,  ireal = ireal, run.list = run.list}	
 
     
    myName = mydata$dataName
    myName = gsub(" ","",myName)
   
	nperiods = mydata$nperiods	
    nperiodsFit = mydata$nperiodsFit

    nRnd = dim(profile)[1]

    nregions = mydata$fit$nregions

    model_coef = mydata$model$coef
    model_onset = mydata$model$onset
    
    model_factor = mydata$model$factor
    
    fit_coef = mydata$fit$coef
    fit_onset = mydata$fit$onset
    
    fit_factor = mydata$fit$factor    

    ## This is the model and fit  raw data

    model_ili = mydata$model$raw
	fit_ili   = mydata$fit$raw
		
    model_rtn_ili = NULL
    model_profile_ili = NULL

    if (!is.null(model_rtn))
        model_rtn_ili = model_rtn/model_factor
    if (!is.null(model_profile)) {
        model_profile_ili = model_profile
        model_profile_ili = model_profile_ili/model_factor
    }

    rtn_ili = rtn
    profile_ili = profile
    for (i in 1:nregions) {
        rtn_ili[, i] = rtn[, i]/fit_factor[i]
        profile_ili[, , i] = profile[, , i]/fit_factor[i]
    }
    
    fit_model = rep(NA, nperiods)
    fit_model_mean = rep(NA, nperiods)
    fit_model_profile = array(0, dim = c(nRnd, nperiods))
    for (i in 1:nperiods) {

        fit_model[i] = sum(rtn[i, 1:nregions])
        tmp = rep(0, nregions)
        for (k in 1:nregions) tmp[k] = mean(profile_ili[, i, k])
        fit_model_mean[i] = sum(tmp[1:nregions] * fit_coef[1:nregions])
        for (irnd in 1:nRnd) {
            for (k in 1:nregions) fit_model_profile[irnd, i] = fit_model_profile[irnd, i] + profile_ili[irnd, i, k] * fit_coef[k]
        }
	}

   # now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, model_rtn = model_rtn, model_profile = model_profile, model_ili = model_ili, model_rtn_ili = model_rtn_ili, model_profile_ili = model_profile_ili, fit_ili = fit_ili, rtn_ili = rtn_ili, profile_ili = profile_ili, fit_model = fit_model, fit_model_mean = fit_model_mean, fit_model_profile = fit_model_profile, run.list = run.list, ireal = ireal, idevice = 1)
        
    filename = pdfName = paste(mydata$subDir, "/profiles-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename) 
    
	cat("\nSaving Profiles to a Binary File : ", filename, "\n")    
         
    return(err = 0)
}

saveTables <- function(mydata = NULL, tables.mod = NULL, tables.fit = NULL, tables.agg.mod = NULL, mod_id = NULL, fit_id = NULL, run.list = NULL, ireal = 1) {


    #' Save Tables with Results of Fitting/Forecasting the Disease Data
    #'
    #' \code{saveTables} Saves and RData binary files with tables of the forecast/fit results as well as the observations
    #' the fit level data
    #' @param mydata - dataframe with all the data for this \pkg{DICE} run
    #' @param tables.mod - a table with the data and the results for the model fit. It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.fit - Optional a table with the data and the results for the fits at the fit_level.
    #' It includes the mean and the 5-95\% percentile (as well as the error)
    #' @param tables.agg.mod - Optional. A table with the aggregate of the uncoupled fit_level results present only in the case of an uncoupled run
    #' @param mod_id The abbreviation of the states/regions-model level
    #' @param fit_id The abbreviation of the states/regions-fit level
    #' @param ireal - Numeric, the number of the MCMC chain
    #' @param run.list  A list with various run parameters
    #' @examples
    #' saveTables(mydata = mydata, tables.mod = tables.mod, tables.fit = NULL,
    #' mod_id = mod_id, fit_id = fit_id, ireal = ireal)
    #' @return err=0 if file was saved
    #'
    #'
        
    
    myName = mydata$dataName
    myName = gsub(" ","",myName)
 
    nperiodsFit = mydata$nperiodsFit

	# now dump all the profiles we have to a file
    dump = list()
    dump = list(mydata = mydata, tables.mod = tables.mod, tables.fit = tables.fit, tables.agg.mod = tables.agg.mod, mod_id = mod_id, fit_id = fit_id, ireal = ireal, run.list = run.list, idevice = 1)

    filename = pdfName = paste(mydata$subDir, "/tables-", myName, "-", nperiodsFit, "-", ireal, ".RData", sep = "")
    save(dump, file = filename) 
    
	cat("\nSaving Tables to a Binary File : ", filename, "\n")    
         
    return(err = 0)

}
