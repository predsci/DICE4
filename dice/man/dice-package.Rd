\name{DICE-package}
\alias{DICE-package}
\alias{DICE}
\docType{package}
\title{
DICE: Dynamics of Interacting Community Epidemics
}
\description{
An arbitrary scale global Modeling Framework for modeling the spread of influenza and other infectious diseases.

The \code{DICE} package implements a global model of influenza within civilian populations with arbitrary spatial resolution, with the objective of making best use of available data at different spatial scales. \code{DICE} uses a deterministic compartmental S-I-R model, applied on various interacting spatial scales, and a robust Markov-Chain-Monte-Carlo (MCMC) fitting procedure that can quickly characterize Influenza-Like-Illness (ILI) incidence profiles, of either coupled or uncoupled regions, in real time. Upon completing a run \code{DICE} provides estimates for the epidemic transmissibility and severity.

Publicly available weekly ILI data is included in the package (from the CDC and Google Flu Trends), along with weekly averaged specific humidity and school opening/closing schedules.  The last two sets of data are needed for three of the five different models for the time dependence of the basic reproduction number that the user can choose from.  These data sets are provided on spatial scales of various, increasing, resolutions.

Upon completion of the MCMC fitting procedure, the \code{DICE} package analyzes the results and produces an extensive set of publication-quality plots (in a png or pdf format) and tables (in csv format).
The complete history of the run and the fitting procedure is saved as an \emph{RData} file which the user can later load and further analyze.


}
\details{
\tabular{ll}{
Package: \tab DICE\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2017-12-24\cr
License: \tab Currently Used by PSI and Leidos Only\cr
}

}
\author{
Predictive Science Inc.

Maintainer:  Michal Ben-Nun <mbennun@predsci.com>
}
\references{

Ben-Nun, M., Riley, P., Turtle, J., Bacon, D. P., and Riley, S., DICE: Dynamics of Interacting Community Epidemics, to be submitted to PLOS – Computational Biology, 2018.

Riley P, Ben-Nun M, Armenta R, Linker JA, Cost AA, et al. (2015) Early Characterization of the Severity and Transmissibility of Pandemic Influenza Using Clinical Episode Data from Multiple Populations. PLoS Comput Biol \url{http://dx.doi.org/10.1371/journal.pcbi.1004392}


Riley P, Ben-Nun M, Armenta R, Linker JA, Eick AA, et al. (2013) Multiple Estimates of Transmissibility for the 2009 Influenza Pandemic Based on Influenza-like-Illness Data from Small US Military Populations. PLoS Comput Biol 9: e1003064. doi:10.1371/journal.pcbi.1003064. \url{http://dx.doi.org/10.1371/journal.pcbi.1003064}


Centers for Disease Control and Prevention Influenza Surveillance: \url{http://www.cdc.gov/flu/weekly/fluactivitysurv.htm}

Google Flu Trends: \url{http://www.google.org/flutrends/us/}

}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
The \code{DICE} manual
}
\examples{

For a \code{DICE} run of the 2015-2016 cdc national data using the ten HHS regions \emph{with} coupling between the regions use:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0)

For a \code{DICE} run of the 2015-2016 cdc national data using the ten HHS regions \emph{without} coupling between the regions use:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 1)

For a \code{DICE} run of the 2014-2015 GFT data for HHS region number 9 \emph{RegState = 9} using state level data with coupling between the states in region 9 use:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2014, mod_level = 3, fit_level = 4, RegState = 9, isingle = 0)

To control which model is used for the basic reproduction number \emph{R0}, set the parameter \emph{model} in your call. Default value is \emph{5}:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2014, mod_level = 3, fit_level = 4, RegState = 9, isingle = 0, model = 3)

To control the number of MCMC chains that the code will run set the parameter \emph{nreal} in your call, default is \emph{1}:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nreal = 3)

To control the number of MCMC steps/trial in each chain set the parameter \emph{nMCMC} in your call, default is \emph{1e5}:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6)

To control the name of the sub-directory where all the output files and plots are saved use the keyword \emph{subDir}, default is \emph{output}:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6, subDir = 'test')

To control the file format for the plots (pdf, png or x11) set the parameter \emph{device}:

output <- runDICE(diseae = 'flu', data_source='cdc', year = 2015, mod_level = 2, fit_level = 3, isingle = 0, nMCMC = 1e6, device = 'pdf')
(The package can accept an array of file formats, i.e. device = c('pdf','png'), in which case more both 'png' and 'pdf' files will be created.)

To run \code{DICE} in a forecast or predictive mode you can set the number of weeks the code uses in the fit to be lower than the number of weeks in the season.
(Note that for the current season it is always running in a predictive mode because the season is not yet completed.)

output <- runDICE(disease = 'flu', data_source='gft', year = 2013, mod_level = 3, fit_level = 4, isingle = 1, nMCMC = 1e6, nfit = 35)

To select only  a few HHS regions and run them coupled (for example using the three east coast regions: 1, 2 and 3) use:
output <- runDICE(disease = 'flu', data_source="cdc", year=2015, mod_name=c(NAME_2="US"), mod_level=2, fit_level=3, fit_names=c("R1","R2","R3"), isingle = 0)

To select only a few states and run them coupled  use for example:
output <- runDICE(disease = 'flu', data_source="cdc",year=2014, mod_name=c(NAME_2="US"), mod_level=2, fit_level=4, fit_names=c("WA","OR","CA"), isingle = 0)

}

