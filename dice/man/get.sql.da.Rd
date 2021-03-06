% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_dengue_fxns.R
\name{get.sql.da}
\alias{get.sql.da}
\title{Data Augmentation for Mechanistic Model}
\usage{
get.sql.da(nfit = NULL, mydata = NULL, all_years_epi = NULL,
  cases = NULL, epi = NULL, corrMAT = NULL, deng.null = NULL,
  my_id = NULL)
}
\arguments{
\item{nfit}{An integer - the number of data points use in the run}

\item{mydata}{A dataframe with all the available data for this \pkg{DICE} run}

\item{all_years_epi}{the epi data for all years}

\item{cases}{An array with the entire history of incidence information for one region/patch}

\item{epi}{An array with the incidence information for one season and one region/patch}

\item{corrMAT}{A 2D array with the Pearson corrletion values for all available disease seasons}

\item{deng.null}{the Historic NULL model - average monthly/weekly number of cases}

\item{my_id}{The abbreviation of the state/region}
}
\value{
epi A list with the original and augmented mydata (if da=1,2) and the weight for the augmented
data points
}
\description{
\code{get.sql.da} Allows the User to augment the incidence time series beyond the \code{nfit} data points
using either the historic NULL model \code{da = 1}, or the most similar past season \code{da=2}.
If \code{da=0}, the data is not augmented.
}
\examples{
get.sql.da(nfit = nfit, mydata = mydata, years = all_years_epi$years, cases = all_years_epi$model$epi,
 epi = mydata$model$epi,corrMAT = corrMAT, distMAT = distMAT, deng.null=deng.null)

get.sql.da(nfit = nfit, mydata = mydata, years = all_years_epi$years, cases = all_years_epi$fit$epi[,iregion],
 epi = mydata$fit$epi[,iregion],corrMAT = corrMAT, distMAT = distMAT, deng.null=deng.null)
}
