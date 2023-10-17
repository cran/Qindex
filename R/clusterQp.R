
#' @title Cluster-Specific Sample Quantiles
#'  
#' @description
#' Obtain vectors of sample \link[stats]{quantile}s in each cluster of observations
#' 
#' @param formula \link[stats]{formula} passed to \link[stats]{aggregate.formula}.
#' To calculate the cluster-specific statistics
#' for response \eqn{y}, the user may use
#' \describe{
#' \item{`y ~ id`}{to retain only the cluster `id` in the returned value}
#' \item{`y ~ id + x1 + x2`}{to retain the cluster `id` and cluster-specific variables \eqn{x_1} and \eqn{x_2} in the returned value}
#' \item{`y ~ .`}{to retain all (supposedly cluster-specific) variables from `data` in the returned value}
#' }
#' 
#' @param data \link[base]{data.frame}
#' 
#' 
# @param FUN \link[base]{character} scalar, name of
# \link[base]{function} for summary statistics, 
# currently only supports \link[stats]{quantile} (`'quantile'` default) 
# \link[np]{npquantile} (`'npquantile'`), or 
# using \link[np]{npquantile} only if sample size is less or equal to 100 (`'npquantile100'`).
#' 
#' @param exclude (optional) \link[stats]{formula} or \link[base]{character} \link[base]{vector}, 
#' (supposedly non-cluster-specific) variables to be excluded from aggregation.
#' To remove variables \eqn{z_1} and \eqn{z_2}, the user may use either
#' \itemize{
#' \item `exclude = c('z1', 'z2')`; or
#' \item `exclude = . ~ . - z1 - z2`
#' }
#' 
#' 
#' 
#' @param from,to,by \link[base]{double} scalars,
#' the starting, end, and increment values
#' to specify a \link[base]{seq}uence of probabilities 
#' \eqn{p = (p_1,\cdots,p_N)'}
#' for the sample \link[stats]{quantile}s \eqn{q = (q_1,\cdots,q_N)'}
#' 
#' @param type \link[base]{integer} scalar, type of \link[stats]{quantile} algorithm
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' Function [clusterQp()] calculates \eqn{N} sample \link[stats]{quantile}s 
#' in each \link[stats]{aggregate}d cluster of observations.
#' The aggregation is specified by parameters `formula` and `exclude`. 
#' 
#' @returns 
#' Function [clusterQp()] returns an \link[stats]{aggregate}d \link[base]{data.frame}.
#' A \link[base]{double} \link[base]{matrix} of \eqn{N} columns is created to store  
#' the sample \link[stats]{quantile}s \eqn{q} of each \link[stats]{aggregate}d cluster.
#' The column names of this quantile \link[base]{matrix} are the probabilities \eqn{p}.
#' 
#' @examples 
#' Ki67q = clusterQp(Marker ~ ., data = Ki67, exclude = c('tissueID','inner_x','inner_y'))
#' tmp = clusterQp(Marker ~ ., data = Ki67, exclude = . ~ . - tissueID - inner_x - inner_y)
#' # stopifnot(identical(Ki67q, tmp))
#' # stopifnot(!anyDuplicated.default(Ki67q$subjID))
#' head(Ki67q)
#' sapply(Ki67q, FUN = class)
#' 
# @importFrom np npquantile
#' @importFrom stats aggregate quantile update.formula terms.formula
#' @export
clusterQp <- function(
    formula,
    data,
    # FUN = c('npquantile100', 'quantile', 'npquantile'),
    exclude,
    from = .01, to = .99, by = .01,
    type = 7,
    ...
) {
  
  # if (anyNA(data)) # okay to have NA in `data`
  
  if (!is.symbol(X <- formula[[2L]])) stop('variable to be aggregated must be `symbol`')
  
  # to create the formula for aggregate
  fom <- if (missing(exclude) || !length(exclude)) {
    formula 
  } else {
    newfom <- if (is.character(exclude) && !anyNA(exclude)) {
      eval(call(
        name = '~', 
        as.symbol('.'),
        Reduce(f = function(e1, e2) call('-', e1, e2), x = lapply(c('.', exclude), FUN = as.symbol))))
    } else if (is.call(exclude) && exclude[[1L]] == '~') {
      exclude
    } else stop('illegal `exclude`')
    update.formula(terms.formula(formula, data = data), new = newfom)
  }

  probs <- seq.int(from = from, to = to, by = by)
  #FUN <- match.arg(FUN)
  ret <- #switch(FUN, quantile = {
    aggregate(fom, data = data, FUN = quantile, probs = probs, type = type, names = FALSE)
  #}, npquantile = {
  #  aggregate(fom, data = data, FUN = FUN, tau = probs)
  #}, npquantile100 = {
  #  aggregate(fom, data = data, FUN = function(x) {
  #    if (length(x) <= 100L) {
  #      npquantile(x, tau = probs)
  #    } else quantile(x, probs = probs, type = type, names = FALSE)
  #  })
  #})
  # when we use ?stats::aggregate and provide a `FUN` which returns a vector
  # we will get a generated new data column of class 'matrix'
  
  colnames(ret[[X]]) <- probs
  return(ret)
  
}






