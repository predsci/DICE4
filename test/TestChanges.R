rm(list=ls())
#detach("package:pmedds.core",unload=TRUE)
require(DICE)

# Function not needed
# PrintError = function(dataType=NULL,job.year=NULL,wflag=NULL,wweight=NULL, myMPZ=NULL, iregion=NULL, national=NULL, job.name=NULL, imodel=NULL, Tg=NULL, optTg=NULL, RandInit=NULL, nMCMC=NULL, nlines=NULL, dr=NULL, reals=NULL, walkers=NULL, iseed=NULL, debug=NULL, verbose=NULL, plot=NULL, device=NULL, MCMC=NULL) {
# 	
# 	
# }


testDICE = function(data_source="cdc",year=2014, disease="flu", mod_name=c(NAME_2="US"), fit_names="all", mod_level=2, fit_level=3, nfit=41, model=5, epi_model=NULL, prior=0, Temp=1, da=0, db_opts=list(DICE_db="predsci", CDC_server=TRUE), arima_model=NULL, covar="sh", covar_lag=1, method='mech', isingle=1, nMCMC=100, nreal=1, plot=1, device="pdf", RegState=NULL, Tg=Tg) {
	
  # browser()
  # test = list(data_source=data_source, year=year, disease=disease, mod_name=mod_name, fit_names=fit_names, mod_level=mod_level, fit_level=fit_level, nfit=nfit, model=model, epi_model=epi_model, prior=prior, Temp=Temp, da=da, db_opts=db_opts, arima_model=arima_model, covar=covar, covar_lag=covar_lag, method=method, isingle=isingle, nMCMC=nMCMC, nreal=nreal, plot=plot, device=device, RegState=RegState)
  # str(test)
  
	stdout = vector('character')
	output = textConnection('stdout','wr',local=TRUE)
	sink(output)
  
	test = try(runDICE(data_source=data_source, year=year, disease=disease, mod_name=mod_name, fit_names=fit_names, mod_level=mod_level, fit_level=fit_level, nfit=nfit, model=model, epi_model=epi_model, prior=prior, Temp=Temp, da=da, db_opts=db_opts, arima_model=arima_model, covar=covar, covar_lag=covar_lag, method=method, isingle=isingle, nMCMC=nMCMC, nreal=nreal, plot=plot, device=device, RegState=RegState, Tg=Tg))

	sink()
	close(output)
	
	return(list(stdout=stdout,test=test))	
}



TestFunction = function(db_name="predsci") {
	# Set default values				
	# nMCMC    = 100
	plot     = 1

	
	data_sources = list('cdc',NULL,'cdc', 'cdc', NULL, "SOM_MH", NULL, "who_flu", NULL, NULL)
	years        = c(2017,2009,2016, 2016, 2012, 2017, 2016, 2013, 2014, 2017)
	diseases     = c("flu", "flu", "flu", "flu", "dengue", "cholera", "yellow_fever", "flu", "chik", "plague")
	mod_name     = list(c(NAME_2="US"), c(NAME_2="US"), c(NAME_2="US"), c(NAME_2="US", NAME_3="R4"), c(NAME_2="BR"), c(NAME_2="SO"), c(NAME_2="BR"), c(NAME_2="US"), c(NAME_2="JM"), c(NAME_2="MG"))
	mod_levels   = c(2,2,2,3,2,2,2,2,2,2)
	fit_levels   = c(3,2,4,4,3,2,2,2,3,2)
	fit_names    = list("all", "all", c("CA", "OR", "WA"), "all", "all", "all", "all","all",c("HA", "KSA", "SC", "SE", "SJ", "TR", "WE"), "all")
	nFits        = c(52,60,41,52,12,52,52,41,52,104)
	models       = c(1,2,3,4,4,4,4,3,4,4)
	Tg           = c(3,3,3,3,10,5,5,3,10,5)
	epi_model    = c(1,1,1,1,2,5,2,1,2,2)
	prior        = c(0,1,0,0,0,0,0,0,0,0)
	Temp         = c(10,1,1,10,1,1,1,1,1,10)
	da           = c(0,0,1,2,0,0,0,0,0,0)
	db_opts      = list(list(DICE_db=db_name, CDC_server=TRUE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=TRUE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE), list(DICE_db=db_name, CDC_server=FALSE))
	arima_model  = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
	covar        = list(FALSE, FALSE, FALSE, "sh", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
	covar_lag    = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
	method       = list('mech', 'mech', 'mech', 'stat', 'mech', 'mech', 'mech', 'mech', 'mech', 'mech')
	isingles     = c(0,1,0,1,0,1,1,1,1,1)
	nMCMC        = c(1e2, 1e2, 1e2, 1e2, 1e2, 5e3, 1e2, 1e2, 1e2, 1e2)
	nreals       = c(1,2,1,1,1,1,1,1,1,1)
	devices      = c("pdf","pdf","pdf","pdf","pdf","pdf","pdf","pdf","pdf","pdf")
	RegStates    = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,NULL,NULL)

	
	NumOpts    = length(data_sources)
	for (ii in 1:NumOpts) {
	# for (ii in 7) {
	  cat("Running test ",ii," of ",NumOpts,"....... ", sep="")

	  # Test runPMEDDS() function
	  test1 = testDICE(data_source=data_sources[[ii]],year=years[ii], disease=diseases[ii], mod_name=mod_name[[ii]], fit_names=fit_names[[ii]], mod_level=mod_levels[ii], fit_level=fit_levels[ii], db_opts=db_opts[[ii]], nfit=nFits[ii], model=models[ii], epi_model=epi_model[ii], prior=prior[ii], Temp=Temp[ii], da=da[ii], arima_model=arima_model[[ii]], covar=covar[[ii]], covar_lag=covar_lag[[ii]], method=method[[ii]], isingle=isingles[ii], nMCMC=nMCMC[ii], nreal=nreals[ii], plot=plot, device=devices[ii], RegState=RegStates[[ii]], Tg=Tg[ii])
	  
	  if (is.character(test1$test)) {
	    cat("Error occurred in execution of runDICE() for the following parameters: \n")
	    dataPar = list(data_source=data_sources[[ii]],year=years[ii], disease=diseases[ii], mod_name=mod_name[[ii]], fit_names=fit_names[[ii]], mod_level=mod_levels[ii], fit_level=fit_levels[ii], db_opts=db_opts[[ii]], RegState=RegStates[[ii]])
	    numericPar = list(nfit=nFits[ii], model=models[ii], Tg=Tg[ii], epi_model=epi_model[ii], prior=prior[ii], Temp=Temp[ii], da=da[ii], arima_model=arima_model[[ii]], covar=covar[[ii]], covar_lag=covar_lag[[ii]], method=method[[ii]], isingle=isingles[ii], nMCMC=nMCMC[ii], nreal=nreals[ii], plot=plot, device=devices[ii])
	    ParList = list(dataPar=dataPar,DICEpar=numericPar)
	    str(ParList)
	    return(list(stdout=test1$stdout,errmsg=test1$test, ParList=ParList))		
	  }
	  cat("no errors found \n")
	}
	
	cat("\n Test complete. \n\n")
	# would be nice to erase the directories created by this script, but would need 'subdir' to be passed from runDICE
	# system("rm -rf ./output/*")
	# system("rm -d output")
}





#----- Main Script -------------------------

db_name="predsci"
# db_name="bsve"

start.time <- proc.time()

output = TestFunction(db_name=db_name)

# stop the clock	
      cat("\n\nElapsed Time: \n\n")
      end.time = proc.time()
      run.time = end.time-start.time
      text.label = "user    system   elapsed"
      text.time  = paste(run.time[[1]],run.time[[2]],run.time[[3]],sep="   ")
      print(end.time - start.time)

#print("Erasing all test files:")
#system("find ./output/* -maxdepth 1 -delete")

#print("post-test operation")
