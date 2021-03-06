% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_mcmc_fxns.R
\name{set.imask}
\alias{set.imask}
\title{Set a Mask for Parameters}
\usage{
set.imask(par_names = NULL, opt.list = NULL)
}
\arguments{
\item{par_names}{- A list with all the \pkg{DICE} parameter names}

\item{opt.list}{A named logical list with TRUE/FALSE values for each parameter}
}
\value{
imask An integer list of length nparam with +1 or -1 values
}
\description{
Given a named logical list of all the parameters that \pkg{DICE} supports
  prepare an integer list of the same length with +1/-1 for parameters
  that are/are not optimized
}
\examples{
set.imask{par_names = par_names, opt.list = opt.list}
}
