% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables_fxns.R
\name{results.table}
\alias{results.table}
\title{Create a List of Tables for Observations and Results}
\usage{
results.table(state_names = NULL, state_id = NULL,
  cadence.names = NULL, cadence = "month", nperiods = NULL)
}
\arguments{
\item{state_names}{The full name of the states/regions}

\item{state_id}{The abbreviation of the states/regions}

\item{cadence.names}{An array with month or week names}

\item{cadence}{String, the cadence of the data: month, week}

\item{nperiods}{Integer - the number of data points in the season}
}
\value{
tables A list with a table for each observation/model/error measure.
}
\description{
\code{results.table} Creates all the tables used to store the observed data and
the results/errors of the SIR/SEIR and SARIMA fits/forecasts
}
\examples{
results.table(state_names = state_names, state_id = state_id, cadence.names = cadence.names, cadennce = cadence)

}
