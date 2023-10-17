
setOldClass('gam')


#' @title Functional Regression Indices & Weights
#' 
#' @description
#' 
#' Functions explained in this documentation are,
#' 
#' \describe{
#' 
#' \item{[FRindex()]}{
#' to compute the functional regression indices and weights based on the functional predictors.}
#' 
#' \item{[predict.FRindex()]}{
#' to compute the predicted values based on functional regression indices and weights model.}
#' 
#' 
#' \item{[FR_gam()]}{
#' a helper function to fit a functional regression model 
#' using generalized additive models with integrated smoothness estimation (\link[mgcv]{gam}).}
#' 
#' }
#' 
#' 
#' @slot formula,data,xarg see explanations in section **Arguments**
#'  
#' @slot gam \link[mgcv]{gam} object, the returned object of helper function [FR_gam()]
#' 
#' @slot sign \link[base]{double} scalar of either 1 or -1, 
#' see Step 5 in section **Details** on function [FRindex()]
#' 
#' @slot index,weight \link[base]{double} \link[base]{vector}s, 
#' functional regression indices and functional regression weights, respectively.
#' See section **Details** on function [FRindex()]
#' 
#' @name FRindex
#' @aliases FRindex-class
#' @export
setClass(Class = 'FRindex', slots = c(
  formula = 'formula',
  data = 'data.frame',
  gam = 'gam',
  sign = 'numeric', # scalar
  index = 'numeric', # vector
  weight = 'numeric', # vector
  xarg = 'numeric' # vector
))






#' @param formula a two-sided \link[stats]{formula}.
#' \describe{
#' \item{Left-hand-side}{is the \link[base]{name} of the response \eqn{y}. 
#' Supported types of responses are \link[base]{double}, \link[base]{logical} and \link[survival]{Surv}.}
#' \item{Right-hand-side}{is the \link[base]{name} 
#' of the tabulated \link[base]{double} \link[base]{matrix} \eqn{X} of functional predictor values.
#' Each row of \eqn{X} represents the tabulated values for a subject.
#' All rows/subjects are tabulated on a common grid `xarg`.
#' Each column of \eqn{X} represents the tabulated values at a point on the common grid for each subject.}
#' }
#' 
#' @param data \link[base]{data.frame}, with 
#' the response \eqn{y} and the tabulated functional predictor values \eqn{X}
#' specified in `formula`.
#' If the functional predictor is the \link[stats]{quantile} function,
#' then `data` is preferably the returned object of [clusterQp()].
#' 
#' @param sign_prob \link[base]{double} scalar between 0 and 1,
#' probability corresponding to the selected nearest-even quantile in `xarg`, 
#' which is used to define the \link[base]{sign} of the functional regression weights.
#' Default is `.5`, i.e., the nearest-even \link[stats]{median} of `xarg`
#' 
#' @param ... for function [predict.FRindex()] and helper function [FR_gam()], 
#' these are currently not in use.
#' For function [FRindex()], see a detailed explanation in section **Using `...` in `FRindex()`**
#' 
#' 
#' @details 
#' 
#' ## Functional regression indices & weights model
#' 
#' Function [FRindex()] defines and calculates 
#' the functional regression indices and weights in the following steps.
#' 
#' \enumerate{
#' 
#' \item Fit a functional regression model to the response \eqn{y} 
#' using the functional predictor \eqn{X}, 
#' with tabulated tabulated on a same grid `xarg` for all subjects,
#' using helper function [FR_gam()]
#' 
#' \item Select one point in the tabulating grid `xarg`.
#' For one-dimensional domain, 
#' we select the nearest-even \link[stats]{quantile} of the tabulating grid `xarg`,
#' corresponding to the user-specified probability `sign_prob`.
#' Default `sign_prob = .5` indicates the \link[stats]{median} of `xarg`.
#' 
#' \item Obtain the fitted coefficient function \eqn{\hat\beta(x)},
#' tabulated on the grid `xarg`, 
#' using internal helper function [gam2beta()]
#' 
#' \item Calculate the integral of the product of 
#' the fitted coefficient function \eqn{\hat\beta(x)} (from Step 3) and
#' the functional predictor values \eqn{X}, 
#' using the \link[pracma]{trapz}oid rule
#' 
#' \item Obtain the \link[base]{sign} of the \link[stats]{cor}relation between 
#' \itemize{
#' \item the subject-specific functional predictor *values*, 
#' at the selected quantile of `xarg` (from Step 2), and
#' \item the subject-specific integrals from Step 4
#' }
#' 
#' }
#' 
#' *Functional regression weights* (slot `@@weight`)
#' are the tabulated weight function on the grid `xarg`.
#' These weights are defined as the product of 
#' `sign` (from Step 5) and \eqn{\hat\beta(x)} (from Step 3).
#' 
#' *Functional regression indices* (slot `@@index`)
#' are defined as the product of 
#' `sign` (from Step 5) and `intg` (from Step 4).
#' Multiplication by `sign` is required to ensure
#' that the resulting functional regression indices
#' are positively associated with the functional predictor values
#' at the selected quantile of `xarg` (from Step 2).
#' 
#' 
#' 
#' @section Using `...` in `FRindex()`:
#' 
#' Function [FRindex()] passes the parameters 
#' `xarg`, `family`, `knot_pct` and `knot.value` 
#' into helper function [FR_gam()] through three dots `...`.
#' 
#' The most important parameter among them is `xarg`.
#' The default argument of the parameter `xarg` comes 
#' from the column names of the \link[base]{matrix} of 
#' tabulated functional predictor values \eqn{X}.
#' This is particularly convenient when 
#' the functional predictor is the \link[stats]{quantile} function, 
#' and `data` is the returned object of function [clusterQp()].
#' 
#' Both [FRindex()] and helper function [FR_gam()] accept user-provided `xarg`.
#' In such case, the provided values will be checked such that
#' \enumerate{
#' \item `xarg` is a \link[base]{numeric} \link[base]{vector} without missingness
#' \item \link[base]{length} of `xarg` is the same as the number of columns of \link[base]{matrix} \eqn{X}
#' \item `xarg` must be strictly sorted (see \link[base]{is.unsorted})
#' }
#' Otherwise, an error message will be returned.
#' 
#' 
#' @section Details of Helper Function: 
#'
#' Helper function [FR_gam()] uses \link[mgcv]{gam} to estimate the functional coefficient by fitting functional regression model. 
#' 
#' 
#' @returns 
#' 
#' ## Functional regression indices & weights model
#' 
#' Function [FRindex()] returns an \link[base]{S4} \linkS4class{FRindex} object. 
#' The slots of \link[base]{S4} class \linkS4class{FRindex} are described in section **Slots**.
#' 
#' 
#' @section Returns of Helper Functions: 
#' Helper function [FR_gam()] returns a \link[mgcv]{gam} object, with additional \link[base]{attributes}
#' 
#' \describe{
#' \item{`attr(,'X')`}{\link[base]{double} \link[base]{matrix} of tabulated functional predictor values \eqn{X}}
#' \item{`attr(,'xarg')`}{\link[base]{double} \link[base]{vector}, see explanation of parameter `xarg`}
#' }
#' 
#' 
#' 
#' @examples 
#' library(survival)
#' 
#' pt = unique(Ki67$PATIENT_ID)
#' length(pt) # 622
#' # set.seed if necessary
#' train_pt = sample(pt, size = 500L)
#' Ki67q = clusterQp(Marker ~ ., data = Ki67, exclude = c('tissueID','inner_x','inner_y'))
#' train_q = subset(Ki67q, PATIENT_ID %in% train_pt)
#' test_q = subset(Ki67q, !(PATIENT_ID %in% train_pt))
#' train_q$Marker = log1p(train_q$Marker)
#' test_q$Marker = log1p(test_q$Marker)
#' 
#' FRi = FRindex(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = train_q)
#' FRi@@index # functional regression index
#' FRi@@weight # functional regression weights
#' head(show(FRi)) # append `FRi` to the data
#' 
#' (FRi_test = predict(FRi, newdata = test_q))
#' 
#' FRi_train = predict(FRi)
#' # stopifnot(identical(FRi@@index, c(FRi_train)), 
#' #  identical(FRi@@weight, attr(FRi_train, 'weight')))
#' 
#' # set.seed if necessary
#' Ki67bbc_v2 = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, 
#'   data = data.frame(train_q, FRi_std = std_IQR(FRi_train)), 
#'   dichotom = ~ FRi_std)
#' summary(Ki67bbc_v2)
#' 
#' @importFrom pracma cumtrapz
#' @importFrom stats cor quantile
#' @rdname FRindex
#' @export
FRindex <- function(
    formula, data,
    sign_prob = .5,
    ... # includes `xarg`
) {
  
  gam_obj <- FR_gam(formula = formula, data = data, ...)
  
  X <- attr(gam_obj, which = 'X', exact = TRUE)
  xarg <- attr(gam_obj, which = 'xarg', exact = TRUE)
  N <- length(xarg)
  
  sign_id <- which(xarg == quantile(xarg, probs = sign_prob, type = 3L))[1L]
  
  beta_ <- gam2beta(gam_obj, n = N)
  
  intg_ <- cumtrapz(x = xarg, y = t.default(X) * beta_)[N, ]
  
  sign_ <- sign(cor(X[, sign_id], intg_, use = 'complete.obs')) # scalar
  
  return(new(
    Class = 'FRindex', 
    formula = formula, data = data,
    gam = gam_obj,
    sign = sign_,
    index = sign_ * intg_, # functional regression index
    weight = sign_ * beta_, # functional regression weight
    xarg = xarg
  ))
  
}



#' @title Show \linkS4class{FRindex} Object
#' 
#' @description
#' Show \linkS4class{FRindex} object.
#' 
#' @param object an \linkS4class{FRindex} object, returned from function [FRindex()]
#' 
#' @returns 
#' The \link[base]{S4} \link[methods]{show} method of \linkS4class{FRindex} object returns a \link[base]{data.frame}.
#' This is the input `data` with the functional predictor values \eqn{X} removed, 
#' and the functional regression indices `@@index` appended.
#' 
#' @keywords internal
#' @export
setMethod(f = show, signature = signature(object = 'FRindex'), definition = function(object) {
  fom <- object@formula
  data <- object@data
  data[[fom[[3L]]]] <- NULL # remove functional predictor values from printing
  data[['FRi']] <- object@index
  return(data)
})





















#' @param xarg strictly increasing \link[base]{double} \link[base]{vector},
#' the common grid on which the functional predictor values \eqn{X} are tabulated
#' 
#' @param family \link[stats]{family} object, the distribution 
#' and link function to be used in \link[mgcv]{gam}.
#' Default family for \link[survival]{Surv} response is `mgcv::cox.ph()`,
#' for \link[base]{logical} response is `binomial(link = 'logit')`,
#' for \link[base]{double} response is `gaussian(link = 'identity')`.
#' 
#' @param knot_pct positive \link[base]{double} scalar, 
#' percentage of the number of columns of \eqn{X},
#' to be used as `knot.value`.  
#' Default is \eqn{40\%}.
#' If `knot.value` is provided by the end-user, then `knot_pct` is ignored.
#' 
#' @param knot.value positive \link[base]{integer} scalar, number of knots 
#' (i.e., parameter `k` in the spline smooth function \link[mgcv]{s})
#' used in \link[mgcv]{gam}.
#' Default is the \link[base]{ceiling} of `knot_pct` of
#' the column dimension of \eqn{X}
#' 
#' @references 
#' 
#' Cui, E., Crainiceanu, C. M., & Leroux, A. (2021). 
#' Additive Functional Cox Model. Journal of Computational and Graphical Statistics. 
#' \doi{10.1080/10618600.2020.1853550}
#' 
#' Gellar, J. E., Colantuoni, E., Needham, D. M., & Crainiceanu, C. M. (2015). 
#' Cox regression models with functional covariates for survival data. Statistical Modelling.
#' \doi{10.1177/1471082X14565526}
#' 
#' @examples 
#' Ki67q = clusterQp(Marker ~ ., data = Ki67, exclude = c('tissueID','inner_x','inner_y'))
#' Ki67q$Marker = log1p(Ki67q$Marker)
#' 
#' library(survival)
#' FR_gam(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = Ki67q)
#' 
#' @importFrom mgcv gam cox.ph s
#' @importFrom stats binomial gaussian
#' @rdname FRindex
#' @export
FR_gam <- function(
    formula, data,
    xarg = as.double(colnames(X)),
    family,
    knot_pct = .4,
    knot.value = ceiling(length(xarg) * knot_pct), 
    ...
) {
  
  if (!is.symbol(rhs <- formula[[3L]]) || !is.matrix(X <- data[[rhs]])) {
    stop('Right-hand-side of `formula` must be a symbol, indicating a matrix column in `data`')
  }
  
  dm <- dim(X)
  force(xarg)
  if (!length(xarg) || !is.numeric(xarg) || anyNA(xarg) || 
      (length(xarg) != dm[2L]) ||
      is.unsorted(xarg, strictly = TRUE)) {
    stop('Provide a strictly increasing numeric grid of length', dm[2L], 'as argument for parameter `xarg` of FR_gam() or FRindex()')
  }
  
  y <- eval(formula[[2L]], envir = data)
  
  ### `Tm`: time indices
  Tm <- tcrossprod(rep(1, times = dm[1L]), xarg)
  # identical to `t.default(array(xarg, dim = c(dm[2L], dm[1L])))`, which is much slower
  
  L <- array(1/dm[2L], dim = dm)
  s_term <- quote(s(Tm, by = L * X, bs = 'cr', k = knot.value))

  gam_cl <- if (inherits(y, what = 'Surv')) {
    call('gam', 
         formula = call('~', quote(y[,1L]), s_term),
         weights = y[,2L], 
         family = if (missing(family)) quote(cox.ph()) else substitute(family))
  } else {
    call('gam', formula = call('~', quote(y), s_term), 
         family = if (!missing(family)) {
           substitute(family) 
         } else if (is.logical(y) || all(y %in% c(0, 1))) {
           quote(binomial(link = 'logit'))
         } else if (is.numeric(y)) {
           quote(gaussian(link = 'identity'))
         } else stop('not supported yet'))
  }
  
  ret <- eval(gam_cl)
  attr(ret, which = 'X') <- X
  attr(ret, which = 'xarg') <- xarg
  return(ret)
}



# future: nonlinear functional regression gam
# nFR_gam





if (FALSE) { # only for developer
  # why no element carries `length(xarg)` ??
  #id_mat = vapply(xx, inherits, what = 'matrix', FUN.VALUE = NA)
  #lapply(xx[id_mat], dim)
  #id_vec = vapply(xx, is.vector, FUN.VALUE = NA)
  #lengths(xx[id_vec])
} # only for developer












#' @param object an \linkS4class{FRindex} object for the \link[stats]{predict} method, 
#' the returned object from function [FRindex()]
#' 
#' @param newdata \link[base]{data.frame}, with at least 
#' the tabulated functional predictor values \eqn{X^{new}}
#' based on `object@@formula`
#' 
#' @param newX \link[base]{double} \link[base]{matrix}, 
#' functional predictor values \eqn{X^{new}} for a set of new subjects.
#' Each row of \eqn{X^{new}} represents the tabulated values for a new subject.
#' All rows/subjects are tabulated on a common grid `new_xarg`.
#' Each column of \eqn{X^{new}} represents the tabulated values at a point on the common grid for each new subject.
#' 
#' @param new_xarg strictly increasing \link[base]{double} \link[base]{vector},
#' the common grid on which the functional predictor values \eqn{X^{new}} are tabulated.
#' The length of `new_xarg` 
#' does not need to be the same as the length of `object@@xarg`, 
#' but they must share the same range.
#' 
#' @details 
#' 
#' ## Predict method for functional regression indices & weights
#' 
#' Function [predict.FRindex()] computes functional regression indices and weights 
#' based on the tabulated functional predictors \eqn{X^{new}} in a new sets of subjects.
#' It's important that the new tabulation grid `new_xarg` must have the same \link[base]{range} 
#' as the model tabulation grid `object@@xarg`.
#' Then,
#' 
#' \enumerate{
#' 
#' \item Obtain the fitted coefficient function \eqn{\hat\beta(x^{new})}
#' of the existing generalized additive model `object@@gam`,
#' but tabulated on the new grid `new_xarg`, 
#' using internal helper function [gam2beta()]
#' 
#' \item Calculate the integral of the product of 
#' the fitted coefficient function \eqn{\hat\beta(x^{new})} (from Step 1) and
#' the new functional predictor values \eqn{X^{new}}, 
#' using the \link[pracma]{trapz}oid rule
#' 
#' }
#' 
#' Predicted functional regression weights
#' are the tabulated weight function on the new grid `new_xarg`.
#' These weights are defined as the product of 
#' `object@@sign` and \eqn{\hat\beta(x^{new})} (from Step 1).
#' 
#' Predicted functional regression indices
#' are defined as the product of 
#' `object@@sign` and `intg` (from Step 2).
#' Multiplication by `object@@sign` is required to ensure
#' that the resulting functional regression indices
#' are positively associated with the functional predictor values
#' at the selected quantile of `object@@xarg`.
#' 
#' 
#' @returns 
#' 
#' ## Predict method for functional regression indices & weights
#' 
#' Function [predict.FRindex()] returns a 
#' \link[base]{double} \link[base]{vector}, 
#' which is the predicted functional regression indices.
#' The returned object contains an \link[base]{attributes}
#' \describe{
#' \item{`attr(,'weight')`}{\link[base]{double} \link[base]{vector}, 
#' the predicted functional regression weights}
#' }
#'  
#' 
#' @rdname FRindex
#' @importFrom stats predict
#' @importFrom pracma cumtrapz
#' @export predict.FRindex
#' @export
predict.FRindex <- function(
    object, newdata = object@data,
    newX = newdata[[object@formula[[3L]]]],
    new_xarg = as.double(colnames(newX)),
    ...
) {
  
  if (!is.matrix(newX)) stop('`newdata` does not contain a matrix column of functional predictor values')
  
  xarg <- object@xarg
  force(new_xarg)
  new_dm <- dim(newX)
  newN <- length(new_xarg)
  
  if (!newN || !is.numeric(new_xarg) || anyNA(new_xarg) ||
      (length(new_xarg) != new_dm[2L]) ||
      is.unsorted(new_xarg, strictly = TRUE)) {
    stop('Provide a strictly increasing numeric grid of length', new_dm[2L], 'as argument for parameter `new_xarg` of predict.FRindex()')
  }
  
  if (!isTRUE(all.equal.numeric(min(xarg), min(new_xarg))) ||
      !isTRUE(all.equal.numeric(max(xarg), max(new_xarg)))) {
    stop('`new_xarg` and `object@xarg` should share the same domain (i.e., range)')
  }
  
  new_beta <- gam2beta(object@gam, n = newN)
  
  new_intg <- cumtrapz(x = new_xarg, y = t.default(newX) * new_beta)[newN, ]
  
  new_FRi <- object@sign * new_intg
  attr(new_FRi, which = 'weight') <- object@sign * new_beta
  
  return(new_FRi)
  
}



if (FALSE) {
  
  # FRindex(y ~ x1 + x2 + x3, xarg = c(1, 2, 3), data = data)
  
  FRindex(y ~ X, xarg = c(1, 2, 3), data = data)
  
  FRindex(endpoint = y, Xs = 101:103, xarg = c(1, 2, 3), data = data)
  
  
  
}


