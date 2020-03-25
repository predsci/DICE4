##
## General setup of parameter lists routines
##
set.param.list <- function(epi_model = 1) {

	#' Create a List of all Parameters the Code Recognizes
	#'
	#'\code{set.param.list} creates a list with all the parameters that the \pkg{DICE} code recognizes. the list is `epi_model' specific.  If the vector states are included it will have additional parameters.
  #' @param epi_model Integer Mechanistic Model Number: 1=SIR, 2=SEIR, 3=(SIR)_H/(SI)_V, 4=(SEIR)_H/(SEI)_V
  #' @return An array with parameter names - the order of parameters will set the order for the min/max arrays also
	#' and the 'mask' for which parameters are optimized (or not)
	#' @examples
	#' set.param.list(epi_model = 2)

	par_names <- c("NH", "Tg", "R0", "sigma",  "pC", "t0", "seed", "e_bckgrnd", "deltaR", "aparam", "alpha", "delta", "ts", "dur")

	## If Vector states are explicitly represented - add the needed parameters

	if (epi_model > 2) {
		par_names <- append(par_names, c("bite", "vec_k", "muV", "T_VH", "T_HV", "sigmaV"))
	}
	
	## SIRB model
	
	if (epi_model == 5)  {
		par_names <- c("NH", "Tg", "birth", 't0', 'pC', "e_bckgrnd",'seed', 'a', 'K','growth_loss','e_shedd', 'delta', 'ts')
	}
	
	return(par_names)

}

set.param.cpl.list <- function( ) {

	#' Create a List of Coupling Parameters that the Code Recognizes
	#'
	#' Create a list with the coupling paramter names
	#' @return
	#' An array with parameter names - the order of parameters will set the order for the min/max arrays also
	#' and the 'mask' for which parameters are optimized (or not)
	#' @examples
	#' set.param.cpl.list()

	par_names_cpl <- c("sd", 'gamma')

	return(par_names_cpl)

}

set.opt.list <- function(mydata = NULL, isingle = 1) {
    #' Create a Logical List  with all the \pkg{DICE} model parameters
    #'
    #' \pkg{DICE} has a list of model parameters it recognizes, for both uncoupled and coupled runs.
    #' This function sets the values of these parameters to either TRUE or FALSE depending on the model number
    #' Given a model number TRUE or FALSE values are assigned to all of these parameters.
    #' In the case of a coupled run, the coupled logical list is also constructed.
    #' @param Compartmental Model Number (1-5) \itemize{
    #'	\item1 - SIR
    #'	\item 2 - SEIR
    #'	\item 3 - (SIR)_H / (SI)_V
    #'	\item 4 - (SEIR)_H / (SEI)_V
    #'  \item 5 - SIRB}
    #' @param model Number - The model Number (1-5, default is 5) \itemize{
    #'   \item 1 - specific humidity (SH) only
    #'   \item 2 - school vacation (SV) only
    #'   \item 3 - both SH and SV
    #'   \item 4 - Fixed value of R(t)
    #'   \item 5 - A two-value model for R(t) }
    #' @param isingle A Number:
    #'   0 - coupled run
    #'   1 - uncoupled run
    #' @return A logical list with all the \pkg{DICE} parameters properly set to TRUE or FALSE.
    #' @examples
    #' set.opt.list{mydata = mydata, isingle = 0}
    #' set.opt.list{mydata = mydata, isingle = 1}
    #'

	epi_model = mydata$epi_model
	
	model = mydata$imodel
	
    opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE,  pC = TRUE, t0 = TRUE, seed = TRUE, e_bckgrnd = TRUE, deltaR = FALSE, aparam = FALSE, alpha = FALSE, delta = FALSE, ts = FALSE, dur = FALSE)

    opt.list.cpl = list()
    if (isingle == 0) {
        opt.list.cpl = list(sd = TRUE, gamma = TRUE)
    } else {
        opt.list.cpl = list(sd = FALSE, gamma = FALSE)
    }

    ## For models with a Human Exposed (E) state
    if (epi_model == 2 || epi_model == 4) opt.list$sigma = TRUE
	
    # In case of a model that includes SH we need to include deltaR and aparam

    if (epi_model == 1 || epi_model == 2) {
    	if (model == 1) {
    		opt.list$deltaR = TRUE
    		opt.list$aparam = TRUE
    		cat("\n Running Model 1- Effect of SH \n")
    	} else if (model == 2) {
    		opt.list$alpha = TRUE
    		cat("\n Running Model 2 - Effect of School Schedule \n")
    	} else if (model == 3) {
    		opt.list$deltaR = TRUE
    		opt.list$aparam = TRUE
    		opt.list$alpha  = TRUE
    		cat("\n Running Model 3 - Effect of SH & School Schedule \n")
    	} else if (model == 4) {
    		cat("\n Running Model 4 - Fixed R(t) \n")
    	} else if (model == 5) {
    		opt.list$delta = TRUE
    		opt.list$ts    = TRUE
    		opt.list$dur   = TRUE
    		cat("\n Running Model 5 - Two value R(t) \n")
    	} else {
    		# default to model 4
    		cat("\n Running Model 5 - Two value R(t) \n")
    	}
    } else if (epi_model == 3 || epi_model == 4)  { ## For explicit Modeling of vector states

    	opt.list$R0 = FALSE

    	opt.list[["bite"]]  = TRUE
    	opt.list[["vec_k"]] = FALSE
    	opt.list[["muV"]]   = FALSE
    	opt.list[["T_HV"]]  = FALSE
    	opt.list[["T_VH"]]  = FALSE
    	opt.list[["sigmaV"]] = FALSE
    	if (epi_model == 4) { ## SEIR (H) - SEI (V)
    		opt.list[["sigmaV"]] = TRUE
    	}
     } else {
     	opt.list = list('NH' = TRUE, 'Tg' = FALSE, 'birth' = FALSE, 't0' = TRUE,  'pC' = FALSE, "e_bckgrnd" = TRUE, 'seed' = TRUE, "a" = TRUE, "K" = TRUE,  'growth_loss' = TRUE, 'e_shedd' = TRUE, 'delta' = FALSE, 'ts' = FALSE)
     	if (mydata$data_source == tolower("ZMB_MH")) {
     		opt.list['delta'] = opt.list['ts'] = TRUE
     		
     	}
     
    }

	##
	## Special case of the Plague, optimize onl: "R0", "sigma", "pC", "t0"
	##
	
	if(tolower(mydata$disease) == 'plague')
	    opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = TRUE,  pC = TRUE, t0 = TRUE, seed = FALSE, e_bckgrnd = FALSE, deltaR = FALSE, aparam = FALSE, alpha = FALSE, delta = FALSE, ts = FALSE, dur = FALSE)
	    
	##
	## Special case of SARS
	## 
	if(tolower(mydata$disease) == 'sars') {
	    opt.list = list(NH = FALSE, Tg = FALSE, R0 = TRUE, sigma = FALSE,  pC = TRUE, t0 = TRUE, seed = TRUE, e_bckgrnd = TRUE) 	
	}
	
    optimize = list(opt.list, opt.list.cpl)
    return(optimize)

}


set.run.list <- function(nreal = 1, nMCMC = 1e+05, nlines = NULL, device = "pdf", subDir = "output", plot = 1) {

    #' Create a List of Parameters For MCMC Procedure
    #'
    #' Pack various parameters related to the MCMC procedure and plotting into a single list.
    #'   These include: The number of MCMC chains, the number of MCMC steps in each chain, the number
    #'   of instances saved for each chain and the device name for plotting the results.
    #' @param nreal Integer  - The number of MCMC chains (default is 1)
    #' @param nMCMC Integer - the number of steps in each MCMC chain (default is 1e5)
    #' @param nlines Integer - the number of instances saved for the history of the chain (default is to save every 100 steps and not to exceed 1e4)
    #' @param device String  - the device name for plotting.  Default is 'pdf' but we also support 'png'
    #' @param subDir String  - the sub-directory name for all output files of the run. Default it 'output'
    #' @param plot Character - TRUE, FALSE or EXTERNAL allows users to use their own plotting routines for the results (can also use 0, 1 and 2)
    #' @return A list packed with these parameters
    #' @examples
    #' set.run.list{nreal = 1, nMCMC = 1e+05, nlines = 1e3, device = 'pdf'}
    #' set.run.list{nreal = 3, nMCMC = 5e+06, nlines = 1e4, device = 'pdf'}

    run.list = list()
    run.list$nreal = nreal
    run.list$nMCMC = nMCMC

    if (is.null(nlines)) {
        nlines = round(nMCMC/100)
        nlines = max(nlines, 100)
        nlines = min(nlines, 10000)
        if (nlines < nMCMC)
            nlines = nMCMC
    }

    run.list$nlines = nlines
    run.list$device = device
    run.list$subDir = subDir
    ithin <- round(nMCMC/nlines)
    ithin = max(1, ithin)

    run.list$ithin = ithin
    run.list$plot = as.character(plot)

    return(run.list)
}
