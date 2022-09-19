

#' @title Cluster-specific sample quantiles for clustered data
#' 
#' @description
#' For a user-supplied sequence of percentiles, calculates sample quantiles in each independent cluster of observations.
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param endpoint \link[base]{character} scalar, column name for the \link[survival]{Surv} endpoint
#' 
#' @param subjID \link[base]{character} scalar, column name for subject/patient ID
#' 
#' @param sampleID (optional) \link[base]{character} scalar, 
#' column name for \code{sampleID}, which is nested within \code{subjID}.
#' When \code{sampleID} is missing, the analysis is performed with only one-level cluster of \code{subjID}
#' 
#' @param Qpredictor \link[base]{character} scalar, column name of the predictor variable
#' 
#' @param probs \link[base]{numeric} vector, probabilities of the \link[stats]{quantile}s 
#' 
#' @param aggregate_sampleID_quantile \link[base]{character} vector, will be implemented in the next release
#' 
#' @param ... additional parameters, currently not in use ..
#' 
#' @return 
#' 
#' \link{sampleQp} returns 
#' data set with one row per \code{subjID}, Np columns with Q(p), p=1,…, Np, and additional
#' subject level variables designated to be kept in the data set.
#' The code should allow for multiple clusters (\code{clusterID}) per subject (\code{subjID}) and output the
#' mean or other summary statistic (according to user input for type of summary statistic) Q(p) for
#' every p.
#' (*) Let us first make use of \link[base]{summary.default}, which cover all options I can think about for now
#' 
#' @examples 
#' # see ?BBC_dichotom
#' 
#' @export
sampleQp <- function(
    data,
    endpoint = 'RECURRENCE',
    subjID = 'PATIENT_ID',
    sampleID,
    Qpredictor = 'Marker',
    probs = seq(from = .05, to = .95, by = .05),
    aggregate_sampleID_quantile = c('median', 'mean', 'min', 'Q1', 'Q3', 'max'),
    ...
) {
  
  if (anyNA(data)) stop('first release of Qindex package does not allow missningness in the data')
  
  if (missing(sampleID)) {
    fom <- eval(call('~', as.symbol(Qpredictor), Reduce(f = function(e1, e2) call('+', e1, e2), lapply(setdiff(names(data), y = c(Qpredictor)), FUN = as.symbol))))
    ret <- aggregate(fom, data = data, FUN = 'quantile', probs = probs)
    
  } else {
    stop('`sampleID` will be implemented in the next release')
    # aggregate quantiles over tissueID
  }
  
  #attr(ret, 'formula') <- 
  return(ret)
  
}



