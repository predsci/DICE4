
#' Pearson Correlation Table for \% ILI 
#'
#' A numeric correlation table with 132 rows and 6 columns showing for each season and each HHS region the season that has the highest correlation with it.
#' This season is referred to as the 'prior season'
#' Correlation is also shown for the national data- which is denoted as region 11.
#' \describe{
#'   \item{Region}{Starting year for the prior flu season for regions 1-10 and national - 11}
#'   \item{year.start}{Starting year for the flu season}
#'   \item{year.end}{Ending year for the flu season}
#'   \item{year.start.prior}{Starting year for the prior flu season}
#'   \item{year.end.prior}{Ending year for the prior flu season}
#'   \item{Cor}{The Pearson correlation between the two seasons}
#'   ...
#' }

"corTable"

#' Pearson Correlation Table for \% ILI - New version
#'
#' A numeric correlation table with 132 rows and 6 columns showing for each season and each HHS region the season that has the highest correlation with it.
#' This season is referred to as the 'prior season'
#' Correlation is also shown for the national data- which is denoted as region 11.
#' \describe{
#'   \item{Region}{Region name,  RegionX and national - United.States}
#'   \item{start}{Starting year for the flu season}
#'   \item{end}{Ending year for the flu season}
#'   \item{R0-mean}{Mean value for R0}
#'   \item{R0-SD}{Standard deviation for R0}
#'   \item{pC-mean}{Mean value for percent clinical, pC}
#'   \item{pC-SD}{Standard deviation for percent clinical, pC}
#'   \item{SH-mean}{Mean value for SH term}
#'   \item{SH-SD}{Standard deviation for SH term}
#'   \item{SV-mean}{Mean value for SV term}
#'   \item{SV-SD}{Standard deviation for SV term}
#'   \item{Cor}{The Pearson correlation between the two seasons}
#'   ...
#' }

"priorTable"


