% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_output_fxns.R
\name{mcmc.onepatch.plot}
\alias{mcmc.onepatch.plot}
\title{Plot the Posterior Distribution of the Parameters - Single Region}
\usage{
mcmc.onepatch.plot(tab.model = NULL, opt.list = NULL,
  run.list = NULL, mydata = NULL, imask = NULL, ireal = 1,
  idevice = 1)
}
\arguments{
\item{tab.model}{The MCMC history of the direct fit of the data}

\item{opt.list}{A logical list of all the parameters \pkg{DICE} recognizes and
can optimize with TRUE/FALSE}

\item{run.list}{a list with parameters used for the MCMC procedure}

\item{mydata}{A dataframe with all the data available for this run}

\item{imask}{An array of integers with +1/-1 values for parameters that are optimized (or not)}

\item{ireal}{Integer - the MCMC chain number}
}
\value{
err   Returns \eqn{err = 0}
}
\description{
Plots the posterior distribution of the MCMC parameters from a  \pkg{DICE} run on a single region/patch.
}
\examples{
mcmc.onepatch.plot{tab.model = tab.model, opt.list = opt.list,
run.list = run.list, mydata = mydata, imask = imask, ireal = 1}
}
