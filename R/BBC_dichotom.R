

#' @title Bootstrap Bias Correction for Dichotomization 
#' 
#' @description Bootstrap-based optimism correction for dichotomizing selected continuous predictor(s)
#' 
#' @param formula \link[stats]{formula}, with a \link[survival]{Surv}, \link[base]{logical}, 
#' or \link[base]{numeric} endpoint 
#' and covariates to be included as is (i.e., without dichotomization).
#' 
#' @param data \link[base]{data.frame}
#'  
#' @param contX \link[base]{character} scalar, 
#' name of the \link[base]{matrix} 
#' containing the continuous predictor(s) to be dichotomized optimally
#' 
#' @param R \link[base]{integer} scalar, number of bootstrap samples, see \link[boot]{boot}. 
#' Default \code{200L}.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' The bootstrap optimism correction procedure
#' is performed as described for a general model selection. 
#' First, \code{R} bootstrap samples are drawn with replacement from the main sample. 
#' In each bootstrap sample, the recursive partitioning tree model is used
#' to establish an objective data-driven optimal cut-point for selected continuous predictors
#' The cut-point from the current bootstrap sample is used to compute the effect size 
#' (log hazards ratio (HR), odds ratio (OR), or coefficient) for
#' each dichotomized predictor in the current bootstrap sample (\strong{bootstrap performance}) 
#' and in the main sample (\strong{test performance}), 
#' and the optimism in log OR/HR or coefficient estimation is computed as the difference between 
#' log OR/HR or coefficient for "Bootstrap performance" and for "Test performance". 
#' The median optimism estimate is computed as the median of optimism estimates over all bootstrap samples. 
#' The cutpoint for dichotomizing each selected continuous predictor is also established in the 
#' main sample and its "apparent performance" is computed as the log OR/HR or coefficient for 
#' dichotomized quantile in the univariate Cox models, logistic regression, or linear regression. 
#' Finally, the optimism-corrected performance estimate is computed by subtracting the median
#' optimism estimate from the apparent performance estimate.
#' 
#' @return 
#' 
#' \link{BBC_dichotom} returns a \link[base]{numeric} bootstrap-based bias adjusted  
#' coefficients of model, corresponding median optimism, \code{R} bootstrap resampling based thresholds, 
#' and a single apparent performance threshold
#' 
#' @references 
#' Ewout W. Steyerberg (2009) Clinical Prediction Models.
#' \doi{10.1007/978-0-387-77244-8}
#' 
#' Frank E. Harrell Jr., Kerry L. Lee, Daniel B. Mark. (1996) Multivariable prognostic models: issues in developing models, evaluating
#' assumptions and adequacy, and measuring and reducing errors.
#' \doi{10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4} 
#' 
#' @examples 
#' # see ?FRindex
#' # see ?eval_dichotom
#' 
#' @export
BBC_dichotom <- function(formula, data, contX, R = 200L, ...) {
  
  if (!inherits(formula, 'formula') || length(formula) == 2L) stop('`formula` must be 2-sided formula')
  
  if (!is.matrix(data[[contX]])) {
    data[[contX]] <- matrix(data[[contX]], ncol = 1L, dimnames = list(NULL, contX))
  }
  
  q <- dim(data[[contX]])[2L] # number of quantiles
  
  #############################################
  ## Apparent performance #
  #############################################
  
  apparent_dichotom <- dichotom_int(edp = eval(formula[[2L]], envir = data), contX = data[[contX]])
  apparent_thresholds <- attr(apparent_dichotom, which = 'thres', exact = TRUE)
  apparent_coef <- model_dichotom(formula = formula, data = data, contX = contX, dQ = apparent_dichotom)
  
  #############################################
  ## Bootstrap  #
  #############################################
  
  thresholds.fun <- function(data, indices) {
    d <- data[indices, ]
    boot_dich <- dichotom_int(edp = eval(formula[[2L]], envir = d), contX = d[[contX]])
    boot_coef <- model_dichotom(formula = formula, data = d, contX = contX, dQ = boot_dich)
    return(c(
      attr(boot_dich, which = 'thres', exact = TRUE), # length `q`
      boot_coef # length `q`
    ))
  }
  
  b_ret <- boot(data = data, statistic = thresholds.fun, R = R)
  b_t <- b_ret$t
  # dim(b_t)
  # head(b_t)
  # dimnames(b_t) # boot::boot removes the dimension names 

  BB_thresholds <- b_t[, 1L:q, drop = FALSE] # threshold
  # dimnames(BB_thresholds) # nothing

  #############################################################
  ## Test performance using BB_thresholds       			#
  ## kth bootstrap cutoffs #
  #############################################################	
  
  test_t <- do.call(rbind, args = lapply(1:R, FUN = function(r) { # (r = 1)
    model_dichotom(formula = formula, data = data, contX = contX, dQ = t.default(t.default(data[[contX]]) > BB_thresholds[r,]))
     # vector of length `q`
  })) # to mimic the behavior of ?boot::boot, which ?base::rbind the bootstrap copies
  # dim(test_t) # R * q
  # dimnames(test_t)
  
  M_Optimism <- b_t[, -(1L:q), drop = FALSE] - test_t # optimistically biased matrix of coefficients
  # dim(M_Optimism) # # R * q
  dimnames(M_Optimism) <- dimnames(test_t)
  
  ## later: maybe trimmed mean may also work instead of median.
  median_optimism <- colMedians(M_Optimism, useNames = TRUE, na.rm = TRUE)
  
  ret <- attr(apparent_coef, which = 'model', exact = TRUE)
  ncf <- length(ret$coefficients)
  ret$coefficients[(ncf-q+1L):ncf] <- ret$coefficients[(ncf-q+1L):ncf] - median_optimism
  ## Subtract the mean optimism estimates from the apparent performance 
  ## estimates to obtain the optimism-corrected performance estimates.
  
  # Tingting: we update the `$coefficients` of `ret`
  # so that the Wald-type z-statistics and p-values can be automatically calculated using ?summary
  # We need to update
  # ret$var
  # we still need cov(ret$coefficients$coefficients, median.optimism.Mbeta)
  # !!!! for now, just leave the variance/covariance as it was !!!
  # end of Tingting
  
  attr(ret, which = 'median_optimism') <- median_optimism
  attr(ret, which = 'BB_thresholds') <- BB_thresholds
  attr(ret, which = 'apparent_thresholds') <- apparent_thresholds
  class(ret) <- c('BBC_dichotom', class(ret))
  return(ret)
  
}






#' @title Dichotomizing
#' 
#' @description 
#' 
#' Generate dichotimized variables by a threshold calculated from \link{rpart} 
#' using \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint and each quantile  
#' 
#' @param edp \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} object, the endpoint
#' 
#' @param contX see \link{BBC_dichotom}
#' 
#' @return 
#' \link{dichotom_int} returns a \link[base]{logical} \link[base]{matrix} of 
#' the same dimension and dimension names as argument \code{contX}, 
#' as dichotomized using the first node of \link[rpart]{rpart} as threshold. 
#' The thresholds (per column) are returned as an \link[base]{numeric} vector 
#' in attribute \code{'thres'}.
#' 
#' @export
dichotom_int <- function(edp, contX) {
  dmQ <- dim(contX)
  dnmQ <- dimnames(contX)
  thres <- setNames(rep(NA_real_, times = dmQ[2L]), nm = dnmQ[[2L]])
  ret <- array(NA, dim = dmQ, dimnames = dnmQ)
  for (i in seq_len(dmQ[2L])) {
    dd <- data.frame(y = edp, x = contX[, i])
    irp <- rpart(formula = y ~ x, data = dd, na.action = na.exclude, control = rpart.control(
      cp = .Machine$double.eps, # to force a split even if the overall lack of fit is not decreased
      maxdepth = 2L
    ))
    if (length(irp$splits)) {
      thres[i] <- irp$splits[1,4]
      ret[, i] <- (contX[, i] > thres[i])
    } # else do nothing
  }
  attr(ret, which = 'thres') <- thres
  return(ret)
}





#' @title Model with Dichotomized Predictors
#' 
#' @description fit Cox proportional hazard model, logistic regression, or linear regression with 
#' dichotomized marker and/or other covariates built as formula in \link{BBC_dichotom} 
#' 
#' @param formula same formula in \link{BBC_dichotom} 
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param contX see \link{BBC_dichotom}
#' 
#' @param dQ returned object from \link{dichotom_int}
#' 
#' @return 
#' \link{model_dichotom} returns a \link[base]{numeric} vector of the concatenated 
#' coefficients.
#' 
#'
#' @export
model_dichotom <- function(
    formula, data, contX, 
    dQ = dichotom_int(edp = eval(formula[[2L]], envir = data), contX = data[[contX]])
) {
  
  #data <- data.frame(data, dQ)
  dQ_nms <- make.names(dimnames(dQ)[[2L]])
  data[dQ_nms] <- dQ
  
  # Removed 2023-01-18
  #ufom <- call('~', quote(.), call('-', quote(.)))
  #for (i in lapply(dQ_nms, FUN = as.symbol)) {
  #  ufom[[3]] <- call('+', ufom[[3]], i)
  #}
  #ucox <- coxph(update.formula(old = formula, new = ufom), data = data)
  # end of removal: end user could select the 'univariable' model by using `formula = edp ~ 1`
  
  y <- formula[[2L]] # endpoint (in language)
  yval <- eval(y, envir = data) # this was programmed for 'Surv' value.  Now to be generalized to *any type* of endpoint
  if (anyNA(yval)) stop('do not allow missingness in the response, for now')
  
  mfom <- call('~', quote(.), quote(.))
  for (i in lapply(dQ_nms, FUN = as.symbol)) {
    mfom[[3]] <- call('+', mfom[[3]], i)
  }
  new_fom <- update.formula(old = formula, new = mfom)
  
  suppressWarnings(mod <- if (inherits(yval, what = 'Surv')) {
    coxph(new_fom, data = data)
  } else if (is.logical(yval) || all(yval %in% c(0, 1))) {
    glm(new_fom, data = data, family = binomial(link = 'logit'))
  } else {
    lm(new_fom, data = data)
  })
  
  coef_ <- tail(mod$coefficients, n = length(dQ_nms))
  attr(coef_, which = 'model') <- mod
  return(coef_)

}




