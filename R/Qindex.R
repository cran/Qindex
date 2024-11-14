
setOldClass('gam')
setOldClass('gam.prefit')

#' @title Sign-Adjusted Quantile Indices
#' 
#' @description
#' Sign-adjusted quantile indices 
#' based on linear and/or nonlinear functional predictors.
#' 
#' @slot .Data \link[base]{double} \link[base]{vector},
#' sign-adjusted quantile indices, see section **Details** of function [integrandSurface]
#' 
#' @slot formula see section **Arguments**, parameter `formula`
#'   
# @slot xgrid strictly increasing \link[base]{double} \link[base]{vector}.
# In package \pkg{Qindex},
# the functional predictor is the \link[stats]{quantile} function, 
# and the input argument for parameter `data` is the returned object of function [clusterQp],
# thus `xgrid`, 
# the common probability grid on which the functional predictor values \eqn{X} are tabulated,
# are the column names of the \link[base]{matrix} \eqn{X}.
#' 
#' @slot gam a \link[mgcv]{gam} object
#' 
#' @slot gpf a `'gam.prefit'` object, which is the returned object 
#' from function \link[mgcv]{gam} with argument `fit = FALSE`
#' 
#' @slot p.value \link[base]{numeric} scalar, 
#' \eqn{p}-value for the test of significance of the functional predictor, 
#' based on slot `@gam`
#' 
#' @slot sign \link[base]{double} scalar of either 1 or -1, 
#' \link[base]{sign}-adjustment, see section **Details** of function [integrandSurface]
#' 
#' @slot sign_prob \link[base]{double} scalar, section **Arguments**, parameter `sign_prob`
#' 
#' @name Qindex
#' @aliases Qindex-class
#' @importFrom methods setClass
#' @export
setClass(Class = 'Qindex', contains = 'numeric', slots = c(
  formula = 'formula',
  gam = 'gam', gpf = 'gam.prefit',
  p.value = 'numeric',
  sign = 'numeric', # scalar
  sign_prob = 'numeric'
  #xgrid = 'numeric' # vector
))


#' @rdname Qindex
#' 
#' @param formula \link[stats]{formula}, e.g., `y~X`. 
#' Response \eqn{y} may be \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.
#' Functional predictor \eqn{X} is a tabulated \link[base]{double} \link[base]{matrix};
#' the rows of \eqn{X} correspond to the subjects, 
#' while the columns of \eqn{X} correspond to a *common tabulating grid* shared by all subjects.
#' The \link[base]{numeric} values of the grid are in the \link[base]{colnames} of \eqn{X}
#' 
#' @param data \link[base]{data.frame}, must be a returned object from function [clusterQp]
#' 
#' @param sign_prob \link[base]{double} scalar between 0 and 1,
#' user-specified probability \eqn{\tilde{p}} 
#' for the nearest-even \link[stats]{quantile} in the grid, 
#' which is used to determine the \link[base]{sign}-adjustment.
#' Default is `.5`, i.e., the nearest-even \link[stats]{median} of the grid
#' 
#' @param family \link[stats]{family} object, 
#' see function \link[mgcv]{gam}.
#' Default values are
#' \itemize{
#' \item `mgcv::cox.ph()` for \link[survival]{Surv} response \eqn{y};
#' \item `binomial(link = 'logit')` for \link[base]{logical} response \eqn{y};
#' \item `gaussian(link = 'identity')` for \link[base]{double} response \eqn{y}
#' }
#' 
#' @param nonlinear \link[base]{logical} scalar, 
#' whether to use nonlinear or linear functional model.
#' Default `FALSE`
#' 
#' @param ... additional parameters for functions \link[mgcv]{s} and \link[mgcv]{ti},
#' most importantly `k`
#' 
# @details 
# Function [Qindex] calculates the sign-adjusted quantile indices in the following steps.
# \enumerate{
# \item Fit a functional model (via \link[mgcv]{gam}) 
# of response \eqn{y} with functional predictor \eqn{X};
# \item Obtain the \link[base]{sign}-adjustment, see section **Details** of function [integrandSurface];
# }
# 
# *Sign-adjusted quantile indices* (slot `@@.Data` of returned object)
# are the product of 
# `sign` (from Step 2) and `gam(.)$linear.predictors` (from Step 1).
# Multiplication by `sign` ensures
# that the sign-adjusted quantile indices
# are positively correlated with the user-selected \eqn{X_{\cdot,j}}.
#' 
#' 
#' @returns 
#' Function [Qindex] returns an \linkS4class{Qindex} object, 
#' which is an instance of an \link[base]{S4} class.
#' See section **Slots** for details.
#' 
#' 
#' @examples 
#' # see ?`Qindex-package`
#' @importFrom stats cor quantile binomial gaussian
#' @importFrom methods new
#' @importFrom mgcv gam cox.ph s ti summary.gam
#' @name Qindex
#' @export
Qindex <- function(
    formula, data,
    sign_prob = .5,
    ...
) {
  
  rhs <- formula[[3L]] # right-hand-side
  X <- data[[rhs]]
  if (!is.symbol(rhs) || !is.matrix(X)) stop('Right-hand-side of `formula` must be a symbol, indicating a matrix column in `data`')
  
  dm <- dim(X)
  xgrid <- as.double(colnames(X))
  if (!is.numeric(xgrid) || anyNA(xgrid) || 
      is.unsorted(xgrid, strictly = TRUE)) {
    stop('`data` needs to be a returned object of function clusterQp()')
  }

  gpf_obj <- Qindex_prefit_(formula = formula, data = data, ...)
  
  gam_obj <- gam(G = gpf_obj, data = data, control = list(keepData = TRUE))
  
  sign_id <- which(xgrid == quantile(xgrid, probs = sign_prob, type = 3L))[1L]
  
  cor_ <- cor(
    x = X[, sign_id], # quantiles of `marker` (at selected quantile) 
    y = gam_obj$linear.predictors, # integration
    use = 'complete.obs'
  ) # scalar
  
  # we use `sign_` to make easy the interpretation
  sign_ <- sign(cor_)
  
  return(new(
    Class = 'Qindex', 
    sign_ * gam_obj$linear.predictors, # sign-adjusted quantile index
    formula = formula,
    #xgrid = xgrid,
    gam = gam_obj, gpf = gpf_obj,
    p.value = summary.gam(gam_obj)$s.table[, 'p-value'],
    sign = sign_,
    sign_prob = sign_prob
  ))
  
}







#' @rdname Qindex
#' @export
Qindex_prefit_ <- function(
    formula, data,
    family,
    nonlinear = FALSE,
    ...
) {
  
  rhs <- formula[[3L]] # right-hand-side
  X <- data[[rhs]]
  if (!is.symbol(rhs) || !is.matrix(X)) stop('Right-hand-side of `formula` must be a symbol, indicating a matrix column in `data`')
  
  dm <- dim(X)
  xgrid <- as.double(colnames(X))
  nxgrid <- length(xgrid)
  if (!is.numeric(xgrid) || anyNA(xgrid) || 
      is.unsorted(xgrid, strictly = TRUE)) {
    stop('`data` needs to be a returned object of function clusterQp()')
  }
  
  y <- eval(formula[[2L]], envir = data)
  
  xgrid_ <- tcrossprod(rep(1, times = dm[1L]), xgrid)
  
  # for numeric integration of the functional term
  L <- array(1/nxgrid, dim = dm)
  
  trm_ <- if (nonlinear) as.call(list(
    quote(ti), 
    quote(xgrid_), rhs, by = quote(L), 
    bs = 'cr', # cubic regression spline
    mc = c( # see ?mgcv::ti; which marginals should have centering constraints applied
      FALSE, 
      TRUE
    ), 
    ...
  )) else as.call(list(
    quote(s),
    quote(xgrid_), by = call('*', quote(L), rhs), bs = 'cr', ...
  ))
  
  gam_cl <- if (inherits(y, what = 'Surv')) call(
    name = 'gam', 
    formula = call(name = '~', quote(y[,1L]), trm_),
    weights = quote(y[,2L]), 
    family = if (missing(family)) quote(cox.ph()) else substitute(family),
    fit = FALSE, data = quote(data)
  ) else call(
    name = 'gam', 
    formula = call(name = '~', quote(y), trm_),
    family = if (!missing(family)) {
      substitute(family) 
    } else if (is.logical(y) || all(y %in% c(0, 1))) {
      quote(binomial(link = 'logit'))
    } else if (is.numeric(y)) {
      quote(gaussian(link = 'identity'))
    } else stop('not supported yet'), 
    fit = FALSE, data = quote(data)
  )
  
  return(eval(gam_cl)) # class 'gam.prefit'

}












#' @title Show \linkS4class{Qindex} Object
#' 
#' @description
#' Show \linkS4class{Qindex} object.
#' 
#' @param object an \linkS4class{Qindex} object
#' 
#' @returns 
#' The \link[base]{S4} \link[methods]{show} method of \linkS4class{Qindex} object 
#' does not have a returned value.
#' 
#' @keywords internal
#' @importFrom methods show
#' @importFrom utils head
#' @export
setMethod(f = show, signature = 'Qindex', definition = function(object) {
  
  cat('p-value from gam: test significance of `marker` as a functional predictor\n')
  print(object@p.value)

  cat('Sign-adjusted quantile index\n')
  print(head(c(object)))
    
  return(invisible())
  
})



