
#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' This function calculates vectors of \code{nQ} sample quantiles in each independent cluster of observations (sample or subject)
#' Equidistant probabilities between 0 and 1 are used (including both ends).
#' 
#' @param data \link[base]{data.frame}
#' 
# @param endpoint \link[base]{character} scalar, column name of the endpoint in \code{data}.
# Currently only \link[survival]{Surv} endpoint is supported.
#' 
#' @param subjID \link[base]{character} scalar, column name of the subject/patient index in \code{data}
#'
# @param sampleID (optional) \link[base]{character} scalar, 
# column name of \code{sampleID}, which is nested within \code{subjID}.
# When \code{sampleID} is missing, the analysis is performed with only one-level cluster of \code{subjID}
#' 
#' @param Qpredictor \link[base]{character} scalar, column name of the predictor variable
#' 
#' @param FUN \link[base]{function} for summary statistics, currently only supports \link[stats]{quantile}
#' 
#' @param include \link[base]{character} \link{vector}, predictors to be included
#' 
#' @param exclude \link[base]{character} \link{vector}, predictors to be excluded
#' 
#' @param nQ \link[base]{integer} scalar, number of equidistant quantiles between 0 and 1 (inclusive). 
#' Default is \code{101L}.
#' 
#' @param type \link[base]{integer} scalar, see \link[stats]{quantile}
#' 
# @param aggregate_sampleID_quantile \link[base]{character} vector, will be implemented in the next release
#' 
#' @param ... additional parameters, currently not in use.
#' 
#' @details 
#' \link{sampleQp} calculates \eqn{N_p} sample quantiles in each independent cluster of observations
#' defined by \code{subjID} or \code{sampleID}. 
#' Sample quantiles are stored in \eqn{N_p} columns for \eqn{Q(p)}, \eqn{p=1, \cdots, Np}  
#' Additional subject-level predictors may be designated to be kept in the data set.
# The function allows for multiple clusters (\code{sampleID}) per subject (\code{subjID}) and output the
# mean or other summary statistic (according to user input for type of summary statistic) 
# for multiple cluster-specific Q(p) per \code{subjID}.
#' 
# (*) Let us first make use of \link[base]{summary.default}, which cover all options I can think about for now
#' 
#' @return 
#' \link{sampleQp} returns a \link[base]{data.frame}, with aggregated \link[survival]{Surv} endpoint and 
#' \link[base]{numeric}, \link[base]{character}, or \link[base]{logical} predictors, 
#' in addition to a \link[base]{matrix} of quantiles.
#'
#' @examples 
#' Ki67_Qps = sampleQp(data = Ki67, subjID = 'PATIENT_ID', Qpredictor = 'Marker')
#' head(Ki67_Qps)
#' sapply(Ki67_Qps, FUN = class)   
#' head(Ki67_Qps$Marker)
#'    
#' PR_Qps = sampleQp(data = PR, subjID = 'PATIENT_ID',  Qpredictor = 'Marker')
#' 
#' @export
sampleQp <- function(
    data,
    subjID = 'PATIENT_ID',
    #sampleID, # to be implemented in future
    Qpredictor = 'Marker',
    FUN = 'quantile', #
    include = setdiff(names(data), y = c(Qpredictor)),
    exclude = NULL,
    nQ = 101L,
    type = 7,
    # aggregate_sampleID_quantile = c('median', 'mean', 'min', 'Q1', 'Q3', 'max'),
    ...
) {
  
  if (anyNA(data)) stop('first release of Qindex package does not allow missningness in the data')
  
  probs <- seq.int(from = 0, to = 1, length.out = nQ)
  
#  if (missing(sampleID)) {
    fom <- eval(call('~', as.symbol(Qpredictor), Reduce(f = function(e1, e2) call('+', e1, e2), lapply(setdiff(include, exclude), FUN = as.symbol))))
    # to create the formula for aggregate
    # left-hand-side is `Qpredictor` (which is 'Marker' in our example)
    # right-hand-side is all other columns in `data`, except for `Qpredictor`
    # ... i.e., for current version of Ki67, they are 'PATIENT_ID' 'RECURRENCE' 'AGE_AT_DX' 'Tstage' 'NodeSt'
    # ... these variables must be unique for each `subjID`
    
    ret <- aggregate(fom, data = data, FUN = FUN, probs = probs, type = type)
    # when we use ?stats::aggregate and provide a `FUN` which returns a vector
    # we will get a generated new data column of class 'matrix'
    if (anyDuplicated.default(ret[[subjID]])) stop('wrong')
    
#  } else {
#    stop('`sampleID` will be implemented in the next release')
    # aggregate quantiles over tissueID
#  }
  
  #attr(ret, 'formula') <- 
  return(ret)
  
}



