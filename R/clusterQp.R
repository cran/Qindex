
#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' Obtain vectors of sample quantiles in each cluster of observations
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param contX \link[base]{character} scalar, 
#' column name of the variable, for which the quantiles will be calculated
#' 
# @param FUN \link[base]{character} scalar, name of
# \link[base]{function} for summary statistics, 
# currently only supports \link[stats]{quantile} (\code{'quantile'} default) 
# \link[np]{npquantile} (\code{'npquantile'}), or 
# using \link[np]{npquantile} only if sample size is less or equal to 100 (\code{'npquantile100'}).
#' 
#' @param include (optional) \link[base]{character} \link[base]{vector}, 
#' names of columns in \code{data} to be \link[stats]{aggregate}d by.
#' Default values are all columns in \code{data} (except for \code{contX})
#' 
#' @param exclude (optional) \link[base]{character} \link[base]{vector}, 
#' names of columns in \code{data} to be excluded from aggregation.
#' Default value is \code{NULL}
#' 
#' @param from,to,by \link[base]{numeric} scalars,
#' the starting, end, and increment values
#' to specify a sequence (via \link[base]{seq.int}) of probabilities 
#' \eqn{\bold{p} = (p_1,\cdots,p_N)'}
#' for the sample quantiles \eqn{\bold{q} = (q_1,\cdots,q_N)'}
#' 
#' @param type \link[base]{integer} scalar, type of quantile algorithm, see \link[stats]{quantile}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' \link{clusterQp} calculates \eqn{N} sample quantiles 
#' in each aggregated cluster of observations.
#' The aggregated clusters are specified by parameters \code{include} and/or \code{exclude} via \link[stats]{aggregate}. 
#' Sample quantiles \eqn{\bold{q}}, for all aggregated clusters, are stored in a \link[base]{matrix} with \eqn{N} columns.  
#' 
#' @returns 
#' \link{clusterQp} returns an aggregated \link[base]{data.frame}, 
#' with a \link[base]{matrix} of sample quantiles named \code{contX}.
#' The column names of the quantile \link[base]{matrix} are the probabilities \eqn{\bold{p}}.
#'
#' @examples 
#' # see ?`Qindex-package`
#' 
# @importFrom np npquantile
#' @importFrom stats aggregate quantile
#' @export
clusterQp <- function(
    data,
    contX = 'Marker',
    # FUN = c('npquantile100', 'quantile', 'npquantile'),
    include = names(data),
    exclude = NULL,
    from = .01, to = .99, by = .01,
    type = 7,
    ...
) {
  
  # if (anyNA(data)) # okay to have NA in `data`
  
  #FUN <- match.arg(FUN)
  
  probs <- seq.int(from = from, to = to, by = by)

  fom <- eval(call('~', 
                   as.symbol(contX), 
                   Reduce(f = function(e1, e2) call('+', e1, e2), 
                          x = lapply(setdiff(include, y = c(contX, exclude)), FUN = as.symbol))))
  # to create the formula for aggregate
  # left-hand-side is `contX` (which is 'Marker' in our example)
  # right-hand-side is all other columns in `data`, except for `contX`
  # ... i.e., for current version of Ki67, they are 'PATIENT_ID' 'RECURRENCE' 'AGE_AT_DX' 'Tstage' 'NodeSt'
  
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
    
  #attr(ret[[contX]], which = 'xarg') <- probs
  colnames(ret[[contX]]) <- probs
  return(ret)
  
}



