


# @references 
# \doi{10.1007/978-0-387-77244-8}  Clinical Prediction Models by Steyerberg
### to add in the next release 




#' @title Bootstrap Bias Correction for Dichotomization 
#' 
#' @description Bootstrap-based optimism correction for dichotomizing quantiles
#' 
#' @param formula \link[stats]{formula}
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param optQ \link[base]{character} scalar
#' 
#' @param R see \link[boot]{boot}
#' 
#' @param seed \link[base]{integer} scalar, random seeds for generating repeated samples, see \link[base]{set.seed}
#' 
#' @param ... ..
#' 
#' @details 
#' 
#' The bootstrap optimism correction procedure (1) 
#' is performed as described for a general model selection in (2). 
#' First, \code{R} bootstrap samples are drawn with replacement from the main sample. 
#' In each bootstrap sample, the recursive partitioning tree model is used
#' to establish an objective data-driven optimal cut-point for a specific optimal quantile. 
#' The cut-point from the current bootstrap sample is used to compute the log odds-ratio (OR) and/or hazards ratio (HR) for
#' dichotomized quantile predictor in the current bootstrap sample (\strong{bootstrap performance}) 
#' and in the main sample (\strong{test performance}), 
#' and the optimism in log OR/HR estimation is computed as the difference between 
#' log OR/HR for “Bootstrap performance” and for “Test performance”. 
#' The median optimism estimate is computed as the median of optimism estimates over all bootstrap samples. 
#' The cutpoint for dichotomizing each selected optimal quantile is also established in the 
#' main sample and its “apparent performance” is computed as the log OR/HR for 
#' dichotomized quantile in the univariate logistic regression or Cox models. 
#' Finally, the optimism-corrected performance estimate is computed by subtracting the mean
#' optimism estimate from the apparent performance estimate.
#' 
#' @return ..
#' 
#' @examples 
#' 
#' \donttest{
#' Ki67_Qps = sampleQp(data = Ki67, endpoint = 'RECURRENCE', subjID = 'PATIENT_ID', 
#'   Qpredictor = 'Marker', probs = seq(from = .05, to = .95, by = .05))
#' sapply(Ki67_Qps, FUN = class)   
#' head(Ki67_Qps$Marker)
#' (Ki67_eval = evalQp(RECURRENCE ~ Marker, data = Ki67_Qps, seeds = 1:100))
#' lapply(Ki67_eval, FUN = dim)
#' summary(Ki67_eval)
#' Ki67_opt = optQp(Ki67_eval, n = 2L)
#' colnames(Ki67_opt$Marker) = paste0('Ki67_', colnames(Ki67_opt$Marker))
#'     
#' PR_Qps = sampleQp(data = PR, endpoint = 'RECURRENCE', subjID = 'PATIENT_ID', 
#'   Qpredictor = 'Marker', probs = seq(from = .05, to = .95, by = .05))
#' head(PR_Qps$Marker)
#' PR_eval = evalQp(RECURRENCE ~ Marker, data = PR_Qps, seeds = 1:100)
#' PR_opt = optQp(PR_eval, n = 3L)
#' colnames(PR_opt$Marker) = paste0('PR_', colnames(PR_opt$Marker))
#' 
#' cellOpt0 = merge(Ki67_opt, PR_opt, by = setdiff(names(Ki67_opt), 'Marker'),
#'   suffixes = c('.Ki67', '.PR'))
#' cellOpt = within(cellOpt0, expr = {
#'   Marker = cbind(Marker.Ki67, Marker.PR)
#'   Marker.Ki67 = Marker.PR = NULL
#' })
#' dim(cellOpt)
#' names(cellOpt)
#' head(cellOpt$Marker)
#' 
#' mod = BBC_dichotom(RECURRENCE ~ NodeSt + Tstage, data = cellOpt, 
#'   optQ = 'Marker', R = 100, seed = 1)
#' names(mod)
#' mod$ucox_dictom
#' mod$mcox_dictom
#' mod$median.optimism
#' }  
#' 
#' @export
BBC_dichotom <- function(formula, data, optQ, R, seed, ...) {
  
  if (!inherits(formula, 'formula') || length(formula) == 2L) stop('`formula` must be 2-sided formula')
  
  if (!is.matrix(data[[optQ]])) stop('optimal markers should be a \'matrix\'')
  
  q <- dim(data[[optQ]])[2L]
  
  #############################################
  ## Apparent performance #
  #############################################
  
  apparent_dQ <- dichotom_optQ(edp = eval(formula[[2L]], envir = data), optQ = data[[optQ]])
  apparent_cox <- coxph_optQ(formula = formula, data = data, optQ = optQ, dQ = apparent_dQ)
  if (FALSE) {
    attr(apparent_dQ, which = 'thres', exact = TRUE) # `thres`, vector of length `q`
    #apparent_cox$u # univariable coefs # 'matrix' of dimension `q` by 1
    #apparent_cox$m # multivariable coefs # 'matrix' of dimension `q` by `x+1` # `x` is the number of other predictors
    apparent_cox$u2 # vector of length `q`
    apparent_cox$m2 # vector of length `x + q`
  }
  
  
  #############################################
  ## Bootstrap  #
  #############################################
  
  thresholds.fun <- function(data, indices) {
    d <- data[indices, ] # allows boot to select sample
    
    boot_dQ <- dichotom_optQ(edp = eval(formula[[2L]], envir = d), optQ = d[[optQ]])
    boot_cox <- coxph_optQ(formula = formula, data = d, optQ = optQ, dQ = boot_dQ)
    
    ret <- c(
      attr(boot_dQ, which = 'thres', exact = TRUE), # length `q`
      #boot_cox$u,
      #boot_cox$m,
      boot_cox$u2, # length `q`
      boot_cox$m2 # length `x + q`
    )
    # vector length q + q + (x+q) = 5 + 5 + (5+2) = 17
    return(ret)
  }
  
  set.seed(seed = seed)
  
  b_ret <- boot(data = data, statistic = thresholds.fun, R = R)
  b_t <- b_ret$t
  # dim(b_t)
  # dimnames(b_t) # boot::boot removes the dimension names 

  thresholds.B <- b_t[, 1L:q] # threshold
  # dimnames(thresholds.B) # nothing

  #############################################################
  ## Test performance using thresholds.B       			#
  ## kth bootstrap cutoffs #
  #############################################################	
  
  test_t <- do.call(rbind, args = lapply(1:R, FUN = function(r) {
    icox <- coxph_optQ(formula = formula, data = data, optQ = optQ, dQ = t.default(t.default(data[[optQ]]) > thresholds.B[r,]))
    c(icox$u2, icox$m2) # vector of length `q` + `x + q`
  })) # to mimic the behavior of ?boot::boot, which ?base::rbind the bootstrap copies
  # dim(test_t) # 100 * (5+(5+2))
  # dimnames(test_t)
  
  M_Optimism <- b_t[, -(1L:q)] - test_t # optimistically biased matrix of coefficients
  # dim(M_Optimism) # # 100 * (5+(5+2))
  dimnames(M_Optimism) <- dimnames(test_t)
  
  ## later: maybe trimmed mean may also work instead of median.
  
  median.optimism.Mbeta <- colMedians(M_Optimism, useNames = TRUE, na.rm = TRUE)
  
  ## Subtract the mean optimism estimates from the apparent performance 
  ## estimates to obtain the optimism-corrected performance estimates.
  ucox <- attr(apparent_cox, which = 'ucox', exact = TRUE)
  ucox$coefficients <- ucox$coefficients - median.optimism.Mbeta[1:q]

  mcox <- attr(apparent_cox, which = 'mcox', exact = TRUE)
  mcox$coefficients <- mcox$coefficients - median.optimism.Mbeta[-(1:q)]
  
  # Tingting: we update `ucox` and `mcox` with its $coefficients
  # so that the Wald-type z-statistics and p-values can be automatically calculated using ?summary
  # We need to update
  # ucox$var and mcox$var
  # we still need cov(mcox$coefficients$coefficients, median.optimism.Mbeta)
  # !!!! for now, just leave the variance/covariance as it was !!!
  # end of Tingting
  
  class(ucox) <- c('dictom.coxph', class(ucox))
  class(mcox) <- c('dictom.coxph', class(mcox))

  return(list(
    'ucox_dictom' = ucox, 
    'mcox_dictom' = mcox, 
    'median.optimism' = list(u = median.optimism.Mbeta[1:q], m = median.optimism.Mbeta[-(1:q)])))
}



#' @title Helper function \link{dichotom_optQ}
#' 
#' @description ..
#' 
#' @param edp \link[survival]{Surv} object, the endpoint
#' 
#' @param optQ \link[base]{matrix} of optimal quantiles
#' 
#' @param control see \link[rpart]{rpart}
#' 
#' @return 
#' \link{dichotom_optQ} returns a \link[base]{logical} \link[base]{matrix} of 
#' the same dimension and dimension names as argument \code{optQ}, 
#' as dichotomized using the first node of \link[rpart]{rpart} as threshold. 
#' The thresholds (per column) are returned as an \link[base]{numeric} vector 
#' in attribute \code{'thres'}.
#' 
#' @export
dichotom_optQ <- function(edp, optQ, control = rpart.control(minbucket = 40)) {
  dmQ <- dim(optQ)
  dnmQ <- dimnames(optQ)
  thres <- setNames(rep(NA_real_, times = dmQ[2L]), nm = dnmQ[[2L]])
  ret <- array(NA, dim = dmQ, dimnames = dnmQ)
  for (i in seq_len(dmQ[2L])) {
    dd <- data.frame(y = edp, x = optQ[, i])
    irp <- rpart(formula = y ~ x, data = dd, na.action = na.exclude, control = control)
    if (length(irp$splits)) {
      thres[i] <- irp$splits[1,4]
      ret[, i] <- (optQ[, i] > thres[i])
    } # else do nothing
  }
  attr(ret, which = 'thres') <- thres
  return(ret)
}


#' @title ..
#' 
#' @description ..
#' 
#' @param formula ..
#' 
#' @param data ..
#' 
#' @param optQ \link[base]{character} scalar
#' 
#' @param dQ returned object from \link{dichotom_optQ}
#' 
#' @return 
#' \link{coxph_optQ} returns a \link[base]{numeric} vector of the concatenated 
#' coefficients of \strong{univariable model} (which is returned in attribute \code{'ucox'}) 
#' and coefficients \strong{multivariable model} (which is returned in attribute \code{'mcox'}).
#' 
#'
#' @export
coxph_optQ <- function(
    formula, data, optQ, 
    dQ = dichotom_optQ(edp = eval(formula[[2L]], envir = data), optQ = data[[optQ]])
) {
  
  #tmp <- lapply(setNames(seq_len(dim(dQ)[2L]), nm = dimnames(dQ)[[2L]]), FUN = function(i) {
  #  idat <- data.frame(data, .dQ = dQ[, i])
  #  ufom <- update(formula, . ~ - . + .dQ)
  #  mfom <- update(formula, . ~ . + .dQ)
  #  list(u = coxph(ufom, data = idat), m = coxph(mfom, data = idat))
  #})
  #mods <- setNames(.mapply(list, dots = tmp, MoreArgs = NULL), nm = c('u', 'm')) # models
  
  dat2 <- data.frame(data, dQ)
  dQ_nms <- make.names(dimnames(dQ)[[2L]])
  ufom2 <- call('~', as.symbol('.'), call('-', as.symbol('.')))
  for (i in lapply(dQ_nms, FUN = as.symbol)) {
    ufom2[[3]] <- call('+', ufom2[[3]], i)
  }
  ucox <- coxph(update(formula, ufom2), data = dat2)
  mfom2 <- call('~', as.symbol('.'), as.symbol('.'))
  for (i in lapply(dQ_nms, FUN = as.symbol)) {
    mfom2[[3]] <- call('+', mfom2[[3]], i)
  }
  mcox <- coxph(update(formula, mfom2), data = dat2)
  
  ret <- #c(lapply(mods, FUN = function(i) {
    #  do.call(rbind, args = lapply(i, FUN = `[[`, 'coefficients'))
    #}), 
    list(
      u2 = ucox$coefficients,
      m2 = mcox$coefficients
    )#)
  # attr(ret, which = 'mods') <- mods
  attr(ret, which = 'ucox') <- ucox
  attr(ret, which = 'mcox') <- mcox
  return(ret)
}




