% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dice.R
\name{fitMulti}
\alias{fitMulti}
\title{Driver Routine - Coupled Spatial Model}
\usage{
fitMulti(mydata = NULL, all_years_epi = NULL, opt.list = NULL,
  run.list = NULL, ireal = 1, iseed = NULL)
}
\arguments{
\item{mydata}{A complex list with all available data for a given disease season model/fit spatial levels and data type}

\item{opt.list}{A Logical list with TRUE or FALSE values for all the parameters supported by \code{DICE}.
These values are based on the user's chosen compartmental model (SIR/SEIR) and the force of infection.}

\item{run.list}{A list with parameters needed for the MCMC procedure}

\item{ireal}{Integer - the MCMC chain number.  Default is 1.}

\item{iseed}{An integer used to set the random-number-generator seed. When not specified, the seed is generated randomly.
Setting the seed to a known integer allows an MCMC chain to be reproducible.}
}
\value{
A list with the following arguments:
\describe{
\item{model_rtn}{The best result of the MCMC procedure for the model data}
\item{model_profile}{Randomly selected results from the MCMC procedure of directly fitting the model data}
\item{tab.model}{The MCMC history of the direct fit to the model data}
\item{rtn}{The best result for indirectly fitting the model data using the fit regions}
\item{profile}{Randomly selected results from the MCMC procedure of indirectly fitting the model data using the fit regions}
\item{tab}{The MCMC history of indirectly fitting the model data using the fit regions}
}
}
\description{
A spatially coupled MCMC fit of the model data. The code first fits the model data directly and then fits it
  as a weighted sum of the coupled fit level data. This fit uses a coupling matrix to describe the interaction between
  different spatial regions and it generates all the fit level profiles at once and minimizes their weighted
  likelihood with the weights given by the relative population of each region. The data can be either cdc or gft data,
  and the model/fit data should have different spatial scales. For example in the case of cdc data: the model is
  national and the fit are the ten HHS regions. Or the model can be an HHS region and the fit fit data is state level data.
}
\examples{
fitMulti{mydata = mydata, opt.list = opt.list, run.list = run.list, ireal = ireal, iseed = iseed}
}
