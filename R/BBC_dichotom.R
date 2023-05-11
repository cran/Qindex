

#' @title Bootstrap-based Optimism Corrected Optimal Dichotomization 
#' 
#' @description 
#' Bootstrap-based optimism correction for optimal dichotomization of selected \link[base]{numeric} predictor(s).
#' 
#' @param formula \link[stats]{formula}, left-hand-side being the endpoint and 
#' right-hand-side being the predictors *in addition to* the \link[base]{numeric} predictor(s) to be dichotomized.
#' If there is no additional predictors, use \code{y ~ 1}
#' 
#' @param data \link[base]{data.frame}
#'  
#' @param contX \link[base]{character} scalar, 
#' name of the \link[base]{matrix} column in \code{data}
#' which contains the \link[base]{numeric} predictor(s) to be dichotomized
#' 
#' @param ... additional parameters of \link[boot]{boot}, most importantly the number of bootstrap replicates \code{R}
#' 
#' @details 
#' 
#' The *apparent performance estimate* are the coefficients 
#' of the regression model with dichotomized predictors (via \link{coef_dichotom})
#' fitted to the entire data set.
#' 
#' The *optimism-corrected performance estimate* is computed by subtracting the 
#' median optimism estimate (via \link{boot_optim_dichotom}) from the apparent performance estimate.
#' 
#' @returns 
#' 
#' \link{BBC_dichotom} returns a \link[base]{numeric} bootstrap-based bias adjusted  
#' coefficients of a \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} regression model,
#' with \link[base]{attributes}
#' \describe{
#' \item{\code{attr(,'median_optimism')}}{the returned object from \link{boot_optim_dichotom}}
#' \item{\code{attr(,'apparent_branch')}}{a \link[base]{list} of \link[base]{language} objects, branches of the apparent model}
#' } 
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
#' # see ?`Qindex-package`
#' 
#' @export
BBC_dichotom <- function(formula, data, contX, ...) {
  
  if (!inherits(formula, 'formula') || length(formula) == 2L) stop('`formula` must be 2-sided formula')
  
  if (!is.matrix(data[[contX]])) {
    data[[contX]] <- matrix(data[[contX]], ncol = 1L, dimnames = list(NULL, contX))
  }
  
  apparent_cf <- coef_dichotom(formula = formula, data = data, contX = contX) # Apparent performance 
  
  median_optimism <- boot_optim_dichotom(formula = formula, data = data, contX = contX, ...) # Bootstrap

  ret <- attr(apparent_cf, which = 'model', exact = TRUE)
  ncf <- length(ret$coefficients)
  q <- dim(data[[contX]])[2L] # number of predictors to be dichotomized
  ret$coefficients[(ncf-q+1L):ncf] <- ret$coefficients[(ncf-q+1L):ncf] - median_optimism
  ## Subtract the mean optimism estimates from the apparent performance 
  ## estimates to obtain the optimism-corrected performance estimates.
  
  # Tingting: we update the `$coefficients` of `ret`
  # so that the Wald-type z-statistics and p-values can be automatically calculated using ?summary
  # We need to update
  # ret$var
  # we still need cov(ret$coefficients, median_optimism)
  # !!!! for now, just leave the variance/covariance as it was !!!
  # end of Tingting
  
  attr(ret, which = 'median_optimism') <- median_optimism
  attr(ret, which = 'apparent_branch') <- attr(apparent_cf, which = 'branch', exact = TRUE)
  class(ret) <- c('BBC_dichotom', class(ret))
  return(ret)
  
}









#' @title Bootstrap-based Optimism Correction
#' 
#' @description 
#' Computes optimism correction for effect sizes corresponding to dichotomized \link[base]{numeric} predictor(s).
#' 
#' @param formula \link[stats]{formula}, left-hand-side being the endpoint and 
#' right-hand-side being the predictors *in addition to* the \link[base]{numeric} predictor(s) to be dichotomized.
#' If there is no additional predictors, use \code{y ~ 1}
#' 
#' @param data \link[base]{data.frame}
#'  
#' @param contX \link[base]{character} scalar, 
#' name of the \link[base]{matrix} column in \code{data}
#' which contains the \link[base]{numeric} predictor(s) to be dichotomized
#' 
#' @param ... additional parameters of \link[boot]{boot}, most importantly the number of bootstrap replicates \code{R}
#' 
#' @details 
#' 
#' \eqn{R} bootstrap samples are generated. For each bootstrap sample, 
#' dichotomizing branches for the \link[base]{numeric} predictors 
#' are optimized in the bootstrap sample using \link{m_rpartD}
#' 
#' \describe{
#' 
#' \item{Bootstrap performance}{Regression coefficients for the dichotomized predictors are estimated 
#' in the bootstrap sample via \link{coef_dichotom}}
#' 
#' \item{Test performance}{Regression coefficients for the dichotomized predictors are estimated  
#' in the entire data set via \link{coef_dichotom}}. The \link[base]{numeric} predictors in the entire data set are dichotomized 
#' using dichotomizing branches determined in the bootstrap sample.
#' 
#' }
#' 
#' The *optimism or optimistic bias* is estimated by the difference between regression coefficients 
#' from bootstrap performance and test performance, for each of the \eqn{R} bootstrap sample.
#' 
#' The *median optimism estimate* is the median of optimism estimates over all \eqn{R} bootstrap samples. 
#' 
#' @returns 
#' 
#' \link{boot_optim_dichotom} returns a \link[base]{numeric} \link[base]{vector} of 
#' median optimism estimate, with \link[base]{attributes}
#' \describe{
#' \item{\code{attr(,'boot_branch')}}{a \link[base]{list} of length \eqn{R}, each element of which is 
#' a \link[base]{list} of \link[base]{language} objects (see attribute \code{'branch'} of \link{m_rpartD})}
#' }
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
#' # see ?`Qindex-package`
#' 
#' @importFrom matrixStats colMedians
#' @export
boot_optim_dichotom <- function(formula, data, contX, ...) {
  
  yval <- eval(formula[[2L]], envir = data)
  
  bootID <- bootIDX(N = length(yval), ...)
  
  boot_dich <- lapply(bootID, FUN = function(i) {
    m_rpartD(y = yval[i], X = data[[contX]][i, , drop = FALSE])
  })
  boot_branch <- lapply(boot_dich, FUN = attr, which = 'branch', exact = TRUE)
  
  boot_cf <- lapply(seq_along(bootID), FUN = function(i) {
    coef_dichotom(formula = formula, data = data[bootID[[i]], ], dX = boot_dich[[i]])
  })
  
  test_cf <- lapply(boot_branch, FUN = function(b) { # (b = boot_branch[[1L]])
    dX <- do.call(cbind, args = lapply(seq_along(b), FUN = function(ib) { # (ib = 1L)
      cl <- b[[ib]]
      cl[[2L]] <- quote(data[[contX]][,ib])
      eval(cl)
    }))
    dimnames(dX) <- dimnames(data[[contX]])
    coef_dichotom(formula = formula, data = data, dX = dX)
  })
  
  # optimistically biased matrix of coefficients
  # .mapply(FUN = `-`, dots = list(boot_cf, test_cf), MoreArgs = NULL) # slower
  M_optim <- do.call(rbind, args = boot_cf) - do.call(rbind, args = test_cf) 
  
  ## later: maybe trimmed mean may also work instead of median.
  ret <- colMedians(M_optim, useNames = TRUE, na.rm = TRUE)
  #dim(M_optim)
  #do.call(rbind, args = test_cf)
  # cov(A, B) - var(B)
  
  attr(ret, which = 'boot_branch') <- boot_branch
  return(ret)
  
}





