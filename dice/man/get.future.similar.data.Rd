% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_fxns.R
\name{get.future.similar.data}
\alias{get.future.similar.data}
\title{Data Augmentation Using the Most Similar Season}
\usage{
get.future.similar.data(mydata = NULL)
}
\arguments{
\item{mydata}{A dataframe with all the available for this \pkg{DICE} run}
}
\value{
An updated mydata structure that has the future data and has weights for each week of data
}
\description{
given the fit and model ILI data find the most similar season for each and grab the data from
these years for future data points > nperiodsData
}
\details{
This function is needed when running in a forecast mode using data augmentation
}
\examples{
prior = get.future.similar.data(mydata = mydata)

}
