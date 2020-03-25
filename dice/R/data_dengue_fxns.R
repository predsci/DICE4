subset.sql.mydata <- function(years = NULL, year.start = 2010, cadence = "month", cadence.start = 1, ntps = 12, cases = NULL) {
    
    #' Subsets the Data
    #'
    #' \code{subset.sql.mydata} selects \code{ntps} time points from the entire data set starting
    #' with the year and month or week selected by the user
    #' @param years A vector of years 
    #' @param year.start Start year for subseting data
    #' @param cadence String, the cadence of the data: month, week, day     
    #' @param cadence.start Start month/week for subseting data
    #' @param ntps Length of subset data    
    #' @param cases  A data frame with year, month/week, ndays and cases information
    #' @return cases A data frame of length ntps starting with the requested month/week and year
    #' @examples
    #' subset.mydata(year.start = 2010, cadence.start = 1, ntps = 12, cases = mydata)
    #'
    #'
    if (is.null(cases)) 
        return
    nlines = length(cases)
    istart = which(years == year.start & cadence == cadence.start)
    iend = (istart + ntps - 1)

	
    if (length(istart) == 0 | length(iend) == 0) {
    	 cases = rep(NA, ntps)
    } else if (iend > nlines) {
    	cases = cases[istart:nlines]
    	while(length(cases) < ntps) cases = c(cases, NA)        
    } else {
        cases = cases[istart:iend]
    }
    
    return(cases)
    
}




