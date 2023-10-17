
setOldClass('gam')

#' @title Nonlinear Functional Regression Indices
#' 
#' @description 
#' 
#' Functions explained in this documentation are,
#' 
#' \describe{
#' 
#' \item{[nlFRindex()]}{
#' to compute the non-linear functional regression indices based on the functional predictors.}
#' 
#' \item{[predict.FRindex()]}{
#' to compute the predicted values based on functional regression indices model.}
#' 
#' }
#' 
#' 
#' @slot formula,data,xarg see explanations in section **Arguments**
#' 
#' @slot gam \link[mgcv]{gam} object
#' 
#' @slot p.value \link[base]{numeric} scalar, 
#' \eqn{p}-value for the test of significance of the functional predictor
#' 
#' @slot index \link[base]{double} \link[base]{vector}, 
#' functional regression indices.
#' 
#' @name nlFRindex
#' @aliases nlFRindex-class
#' @export
setClass(Class = 'nlFRindex', slots = c(
  formula = 'formula',
  data = 'data.frame',
  gam = 'gam',
  p.value = 'numeric',
  index = 'numeric', 
  xarg = 'numeric'
))





#' @rdname nlFRindex
#' 
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
#' @param xarg \link[base]{numeric} \link[base]{vector}.
#' The default argument comes from the column names of the \link[base]{matrix} of 
#' tabulated functional predictor values \eqn{X}.
#' This is particularly convenient when 
#' the functional predictor is the \link[stats]{quantile} function, 
#' and `data` is the returned object of function [clusterQp()].
#' The user-provided `xarg` will be checked such that
#' \enumerate{
#' \item `xarg` is a \link[base]{numeric} \link[base]{vector} without missingness
#' \item \link[base]{length} of `xarg` is the same as the number of columns of \link[base]{matrix} \eqn{X}
#' \item `xarg` must be strictly sorted (see \link[base]{is.unsorted})
#' }
#' Otherwise, an error message will be returned.
#' 
#' @param family ..
#' 
#' @param fit \link[base]{logical} scalar, see \link[mgcv]{gam}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' ## Functional regression indices & weights model
#' 
#' Function [nlFRindex()] fits a non-linear functional regression model to the response \eqn{y} 
#' using the functional predictor \eqn{X}, 
#' with values tabulated on a same grid `xarg` for all subjects (Cui et al, 2021).
# *Functional regression indices* (slot `@@index`)
# are defined as **`gam$linear.predictors`**.
#' 
#' 
#' 
#' @returns 
#' 
#' ## Functional regression indices & weights model
#' 
#' Function [nlFRindex()] returns an \link[base]{S4} \linkS4class{nlFRindex} object.
#' The slots of \link[base]{S4} class \linkS4class{nlFRindex} are described in section **Slots**.
#' 
#' 
#' @references 
#' 
#' Cui, E., Crainiceanu, C. M., & Leroux, A. (2021). 
#' Additive Functional Cox Model. Journal of Computational and Graphical Statistics. 
#' \doi{10.1080/10618600.2020.1853550}
#' 
#'  
#' @examples
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
#' # using Cox model
#' m = nlFRindex(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = train_q)
#' m@@p.value # test significance of `Marker` as a functional predictor
#' train_index = predict(m, newdata = train_q) # non-linear FR index of training data
#' # stopifnot(identical(train_index, m@@index))
#' predict(m, newdata = test_q) # non-linear FR index of test data
#' 
#' # using logistic regression model
#' nlFRindex(RECURRENCE ~ Marker, data = train_q)
#' 
#' # using Gaussian model
#' nlFRindex(RECFREESURV_MO ~ Marker, data = train_q)
#' 
#' @importFrom mgcv gam cox.ph ti summary.gam
#' @importFrom stats binomial gaussian
#' @export
nlFRindex <- function(
    formula, data,
    xarg = as.double(colnames(X)),
    family,
    fit = TRUE, # either FALSE/TRUE
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
  
  S <- length(xarg) ## number of functional observations per subject
  
  ## Add variables related to numerical integration
  ### lmat: numeric integration of the functional term
  L <- array(1/dm[2L], dim = dm)
  
  ### tmat: time indices, xarg
  Tm <- tcrossprod(rep(1, times = dm[1L]), xarg)
  
  ti_term <- #quote(ti(Tm, Marker, by = L, bs = c('cr', 'cr'), mc = c(FALSE, TRUE)))
    call('ti', quote(Tm), rhs, by = quote(L), bs = c('cr', 'cr'), mc = c(FALSE, TRUE))
  
  gam_cl <- if (inherits(y, what = 'Surv')) {
    call('gam', 
         formula = call('~', quote(y[,1L]), ti_term),
         weights = quote(y[,2L]), 
         family = if (missing(family)) quote(cox.ph()) else substitute(family),
         fit = fit,
         data = quote(data))
  } else {
    call('gam', formula = call('~', quote(y), ti_term), 
         family = if (!missing(family)) {
           substitute(family) 
         } else if (is.logical(y) || all(y %in% c(0, 1))) {
           quote(binomial(link = 'logit'))
         } else if (is.numeric(y)) {
           quote(gaussian(link = 'identity'))
         } else stop('not supported yet'),
         fit = fit,
         data = data)
  }
  
  gam_obj <- eval(gam_cl)
  
  if (fit) {
    return(new(
      Class = 'nlFRindex', 
      formula = formula, data = data,
      gam = gam_obj,
      index = gam_obj$linear.predictors,
      p.value = summary.gam(gam_obj)$s.table[, 'p-value'],
      xarg = xarg
    ))
  } else {
    return(list(
      xarg = xarg, 
      gam.prefit = gam_obj
    )) # temperary return for [predict.nlFRindex]
    #new(Class = 'nlFRindex.prefit', 
    #    formula = formula, data = data,
    #    gam.prefit = eval(gam_cl),
    #    xarg = xarg)
  }
  
}



#' @param object an \linkS4class{nlFRindex} object for the \link[stats]{predict} method, 
#' the returned object from function [nlFRindex()]
#' 
#' @param newdata \link[base]{data.frame}, with at least 
#' the tabulated functional predictor values \eqn{X^{new}}
#' based on `object@@formula`
#' 
#' @details 
#' 
#' ## Predict method for non-linear functional regression indices
#' 
#' Function [predict.nlFRindex()] computes non-linear functional regression indices
#' based on the tabulated functional predictors \eqn{X^{new}} in a new sets of subjects.
#' It's important that the new tabulation grid must be exactly the same  
#' as the model tabulation grid `object@@xarg`.
#' 
#' 
#' 
#' 
#' @returns 
#' 
#' ## Predict method for non-linear functional regression indices
#' 
#' Function [predict.nlFRindex()] returns a 
#' \link[base]{double} \link[base]{vector}, 
#' which is the predicted non-linear functional regression indices.
#'  
#' 
#' @rdname nlFRindex
#' @importFrom stats predict
#' @export predict.FRindex
#' @export
predict.nlFRindex <- function(object, newdata, ...) {
  if (identical(object@data, newdata)) return(object@index)
  
  new <- nlFRindex(formula = object@formula, data = newdata, fit = FALSE)
  if (!identical(new$xarg, object@xarg)) stop('`xarg` of new and old data not identical')
  as.numeric(new$gam.prefit$X %*% object@gam$coefficients)
}



