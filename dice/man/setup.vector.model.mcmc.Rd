% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_mcmc_fxns.R
\name{setup.vector.model.mcmc}
\alias{setup.vector.model.mcmc}
\title{Setup Parameters for MCMC Procedure - Model Dengue Data}
\usage{
setup.vector.model.mcmc(mydata = NULL, opt.list = NULL,
  run.list = NULL, tps = NULL, par_names = NULL)
}
\arguments{
\item{mydata}{- dataframe with all the data for this \pkg{DICE} run}

\item{opt.list}{A Logical list with TRUE or FALSE values for all the parameters
supported by \pkg{DICE}. These values are based on the model used for the basic
reproduction number.}

\item{run.list}{A list with parameters needed for the MCMC procedure}

\item{tps}{A numeric array of days for the disease season - the day numbers are
consistent with the weeks/months}

\item{par_names}{- An array with the all parameters ordered as required by DICE}
}
\value{
A list with the following  arguments
\describe{
\item{par_min}{Minimum values for all the parameters supported by \pkg{DICE}}
\item{par_max}{Maximum values for all the parameters supported by \pkg{DICE}}
\item{pardx}{Step-size for MCMC for all the parameters supported by \pkg{DICE}}
\item{par}{Initial values for all the parameters supported by \pkg{DICE}}
\item{logbase}{Base for log values - currently code assumes base 10}
\item{logvec}{Array of integers with 1, use log base, or 0 - do not}
\item{nparam}{Integer-the total number of parameters recognized by the \pkg{DICE} code}
\item{nopt}{Integer-the number of parameters that will be optimized}
}
}
\description{
Setup all the required parameters for the MCMC procedure on the model data.
  These include: The min/max and step size values for all the parameters,
  the initial values for all the parameters, a vector with the list of parameters
  that will be optimized, the total number of parameters and the number of
  parameters that will be optimized.  The code also allocates an array, tab,
  where the history of the MCMC chain is recorded.
}
\examples{
setup.model.mcmc{mydata = mydata,
run.list = run.list, opt.list = opt.list, tps = tps, par_names = par_names}
}
