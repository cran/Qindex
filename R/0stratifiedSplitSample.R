
#' @title Stratified Split Sampling
#' 
#' @description 
#' Split sampling, stratified based on the type of the endpoint.
#' 
#' @param x a \link[base]{numeric} \link[base]{vector}, 
#' a \link[base]{logical} \link[base]{vector},
#' a \link[base]{factor}, 
#' or a \link[survival]{Surv} object, 
#' endpoint of the analysis, on which the stratification will be based
#' 
#' @param nsplit \link[base]{integer} scalar, copies of sample split to be performed
#' 
#' @param trainFrac \link[base]{numeric} scalar between 0 and 1, 
#' fraction of the training set, default \code{.8}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' Split a \link[survival]{Surv} endpoint, stratified by the censoring status.
#' Specifically, 
#' split subjects with observed event into a training and a test set with training set fraction \code{trainFrac},
#' and split the censored subjects into a training and a test set with training set fraction \code{trainFrac}.
#' Then combine the training sets from subjects with observed events and censored subjects,
#' and combine the test sets from subjects with observed events and censored subjects.
#' 
#' Split a \link[base]{logical} endpoint, stratified by the endpoint itself.
#' Specifically, 
#' split the subjects with \code{TRUE} endpoint into a training and a test set with training set fraction \code{trainFrac},
#' and split the subjects with \code{FALSE} endpoint into a training and a test set with training set fraction \code{trainFrac}.
#' Then combine the training sets, and the test sets, in a similar fashion as described above.
#' 
#' Split a \link[base]{factor} endpoint, stratified by the level of the endpoint.
#' Specifically, 
#' split the subjects in each \link[base]{levels} of the endpoint into a training and a test set by ratio \code{trainFrac}.
#' Then combine the training sets, and the test sets, from all levels of the endpoint.
#' 
#' Split a \link[base]{numeric} endpoint into a training and a test set by ratio \code{trainFrac}, without stratification.
#' 
#' @returns 
#' \link{stratifiedSplitSample} returns a \link[base]{list} of length \code{nsplit}, each element of which 
#' is a \link[base]{list} of two elements
#' \describe{
#' \item{\code{train}}{\link[base]{integer} \link[base]{vector}, \link[base]{sort}ed indexes of the training set}
#' \item{\code{test}}{\link[base]{integer} \link[base]{vector}, \link[base]{sort}ed indexes of the test set}
#' }
#' 
#' @note
#' \code{caTools::sample.split} is not what we need.
#' 
#' @examples
#' stratifiedSplitSample(x = rep(c(TRUE, FALSE), times = c(20, 30)), nsplit = 5L)
#' 
#' @export 
stratifiedSplitSample <- function(x, nsplit, trainFrac = .8, ...) {
  
  if (anyNA(x)) stop('do not allow missingness in the response, for now')
  
  if (inherits(x, what = 'Surv')) {
    # stratified by censoring status
    if (dim(x)[2L] == 3L) stop('3-col Surv endpoint not supported yet')
    endpoint <- 'Surv'
    n <- dim(x)[1L]
    xevent <- as.logical(x[,2L])
    idx <- list(
      which(!xevent), # 'integer' indexes of censored events
      which(xevent) # 'integer' indexes of observed events
    )

  } else if (is.logical(x) || all(x %in% c(0, 1))) {
    # stratified by the binary response
    n <- length(x)
    endpoint <- 'logical'
    x <- as.logical(x)
    idx <- list(
      which(!x), # 'integer' indexes of non-responder
      which(x) # 'integer' indexes of responder
    )

  } else if (is.factor(x)) {
    n <- length(x)
    endpoint <- 'factor'
    idx <- lapply(seq_along(attr(x, which = 'levels', exact = TRUE)), FUN = function(i) {
      which(unclass(x) == i)
    })
    
  } else if (is.vector(x, mode = 'numeric')) {
    n <- length(x)
    endpoint <- 'numeric'
    idx <- seq_len(n)
    
  } else stop('unsupported class of `x`: ', class(x)[1L])
  
  ret <- replicate(n = nsplit, expr = {
    idx_train <- lapply(idx, FUN = function(iidx) {
      sample(iidx, size = floor(length(iidx) * trainFrac), replace = FALSE)
    })
    train <- sort.int(unlist(idx_train, use.names = FALSE))
    list(train = train, test = sort.int(setdiff(seq_len(n), y = train)))
  }, simplify = FALSE)
  
  attr(ret, which = 'endpoint') <- endpoint
  return(ret)
  
}

