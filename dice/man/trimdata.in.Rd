% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_user_data_fxns.R
\name{trimdata.in}
\alias{trimdata.in}
\title{Trim an array of raw cases to find the last observed data point}
\usage{
trimdata.in(longvec)
}
\arguments{
\item{longvec}{array with raw cases numbers}
}
\value{
nperiodsData the number of last observed data point
}
\description{
Trim an array of raw cases to find the last observed data point
}
\examples{
nperiodsData <- trimdata.in(longvec = cases)

}
