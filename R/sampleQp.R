
#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' This function calculates vectors of sample quantiles in each independent cluster of observations (sample or subject)
#' Equidistant probabilities between user provided p_{min} and p_{max} are used (including both ends).
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param subjID \link[base]{character} scalar, column name of the subject/patient index in \code{data}
#'
# @param sampleID (optional) \link[base]{character} scalar, 
# column name of \code{sampleID}, which is nested within \code{subjID}.
# When \code{sampleID} is missing, the analysis is performed with only one-level cluster of \code{subjID}
#' 
#' @param Qpredictor \link[base]{character} scalar, column name of the predictor variable
#' 
# @param FUN \link[base]{character} scalar, name of
# \link[base]{function} for summary statistics, 
# currently only supports \link[stats]{quantile} (\code{'quantile'} default) 
# \link[np]{npquantile} (\code{'npquantile'}), or 
# using \link[np]{npquantile} only if sample size is less or equal to 100 (\code{'npquantile100'}).
#' 
#' @param include \link[base]{character} \link{vector}, predictors to be included in output data set
#' 
#' @param exclude \link[base]{character} \link{vector}, predictors to be excluded from output data set
#' 
#' @param from,to,by see \link[base]{seq} a sequence of probabilities with starting, ending values, and interval
#' to calculate corresponding quantiles 
#' 
#' @param type \link[base]{integer} scalar, type of the formula for quantiles,see \link[stats]{quantile}
#' 
# @param aggregate_sampleID_quantile \link[base]{character} vector, will be implemented in the next release
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' \link{sampleQp} calculates \eqn{N_p} sample quantiles in each independent cluster of observations
#' defined by \code{subjID} or \code{sampleID}. 
#' Sample quantiles are stored in \eqn{N_p} columns for \eqn{Q(p)}, \eqn{p=1, \cdots, Np}  
#' Additional subject-level predictors may be designated to be kept in the output data set.
#' Observation-level predictors should be excluded from the input data set.
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
#' Ki67_Qps = sampleQp(data = Ki67, subjID = 'PATIENT_ID', 
#'   exclude = c('tissueID','inner_x','inner_y'), Qpredictor = 'Marker')
#' head(Ki67_Qps)
#' sapply(Ki67_Qps, FUN = class)   
#' head(Ki67_Qps$Marker)
#'
#' 
# @importFrom np npquantile
#' @export
sampleQp <- function(
    data,
    subjID = 'PATIENT_ID',
    #sampleID, # to be implemented in future
    Qpredictor = 'Marker',
    # FUN = c('npquantile100', 'quantile', 'npquantile'), #
    include = setdiff(names(data), y = c(Qpredictor)),
    exclude = NULL,
    
    from = .01,
    to = .99,
    by = .01,
    
    type = 7,
    # aggregate_sampleID_quantile = c('median', 'mean', 'min', 'Q1', 'Q3', 'max'),
    ...
) {
  
 # if (anyNA(data)) stop('first release of Qindex package does not allow missningness in the data')
  
  #FUN <- match.arg(FUN)
  
  probs <- seq.int(from = from, to = to, by = by)
  
  
#  if (missing(sampleID)) {
    fom <- eval(call('~', as.symbol(Qpredictor), Reduce(f = function(e1, e2) call('+', e1, e2), lapply(setdiff(include, exclude), FUN = as.symbol))))
    # to create the formula for aggregate
    # left-hand-side is `Qpredictor` (which is 'Marker' in our example)
    # right-hand-side is all other columns in `data`, except for `Qpredictor`
    # ... i.e., for current version of Ki67, they are 'PATIENT_ID' 'RECURRENCE' 'AGE_AT_DX' 'Tstage' 'NodeSt'
    # ... these variables must be unique for each `subjID`
    
    ret <- #switch(FUN, quantile = {
      aggregate(fom, data = data, FUN = quantile, probs = probs, type = type)
    #}, npquantile = {
    #  aggregate(fom, data = data, FUN = FUN, tau = probs)
    #}, npquantile100 = {
    #  aggregate(fom, data = data, FUN = function(x) {
    #    if (length(x) <= 100L) {
    #      npquantile(x, tau = probs)
    #    } else quantile(x, probs = probs, type = type)
    #  })
    #})
    # when we use ?stats::aggregate and provide a `FUN` which returns a vector
    # we will get a generated new data column of class 'matrix'
    if (anyDuplicated.default(ret[[subjID]])) stop('wrong')
    attr(ret, which = 'xarg') <- probs
    
#  } else {
#    stop('`sampleID` will be implemented in the next release')
    # aggregate quantiles over tissueID
#  }
  
  #attr(ret, 'formula') <- 
  return(ret)
  
}



