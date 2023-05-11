
#' @title Dichotomize Multiple Predictors via Recursive Partitioning
#' 
#' @description 
#' To dichotomize multiple \link[base]{numeric} predictors 
#' via recursive partitioning \link[rpart]{rpart}.
#' 
#' @param y a \link[survival]{Surv} object, 
#' a \link[base]{logical} \link[base]{vector}, 
#' or \link[base]{numeric} \link[base]{vector}, the endpoint
#' 
#' @param X \link[base]{numeric} \link[base]{matrix}, 
#' columns being the predictors to be dichotomized
#' 
#' @param ... additional parameters of \link{rpartD}, currently not in use
#' 
#' @details
#' 
#' For each column of \code{X},
#' 
#' \itemize{
#' 
#' \item Find the dichotomizing branch of the endpoint \code{y} using the predictor \code{X[,i]}, 
#' via \link{rpartD}
#' 
#' \item Dichotomize \code{X[,i]} using this branch
#' 
#' }
#' 
#' @returns 
#' \link{m_rpartD} returns a \link[base]{logical} \link[base]{matrix} of 
#' the same dimension and dimension names as the argument \code{X}, 
#' with \link[base]{attributes}
#' \describe{
#' \item{\code{attr(,'branch')}}{\link[base]{list} of \link[base]{language} objects 
#' (i.e., return of \link{rpartD}), 
#' dichotomozing branches for each predictor}
#' }
#' 
#' @examples
#' y = rnorm(30)
#' X = array(rnorm(30*4), dim = c(30, 4), dimnames = list(NULL, letters[1:4]))
#' m_rpartD(y, X)
#' 
#' @export
m_rpartD <- function(y, X, ...) {
  dm <- dim(X)
  dnm <- dimnames(X)
  
  cseq <- seq_len(dm[2L]) # column sequence
  names(cseq) <- dnm[[2L]]
  
  branch <- lapply(cseq, FUN = function(i) {
    rpartD(y = y, x = X[,i], ...)
  })
  
  ret <- do.call(cbind, lapply(cseq, FUN = function(i) {
    cl <- branch[[i]]
    cl[[2L]] <- quote(X[,i])
    eval(cl)
  }))
  
  attr(ret, which = 'branch') <- branch
  return(ret)
}





#' @title Regression Coefficients of Dichotomized Predictors
#' 
#' @description 
#' Regression coefficients of dichotomized predictors
#' 
#' @param formula \link[stats]{formula}, left-hand-side being the endpoint and 
#' right-hand-side being the predictors *in addition to* the \link[base]{numeric} predictor(s) to be dichotomized.
#' If there is no additional predictors, use \code{y ~ 1}
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param contX \link[base]{character} scalar, 
#' name of the \link[base]{matrix} column in \code{data}
#' which contains the \link[base]{numeric} predictor(s) to be dichotomized.
#' This parameter is not needed if parameter \code{dX} is provided
#' 
#' @param dX (optional) \link[base]{logical} \link[base]{matrix}, dichotomized \code{data[[contX]]}.
#' Default is the returned object from \link{m_rpartD}, 
#' 
#' @details
#' ..
#' 
#' 
#' @returns 
#' \link{coef_dichotom} returns a \link[base]{numeric} \link[base]{vector} of the
#' coefficients of dichotomized predictor(s), with \link[base]{attributes}
#' \describe{
#' \item{\code{attr(,'model')}}{the \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} regression model}
#' \item{\code{attr(,'branch')}}{see \link{m_rpartD}}
#' }
#' 
#'
#' @importFrom survival coxph
#' @importFrom stats lm glm binomial
#' @importFrom utils tail
#' @export
coef_dichotom <- function(
    formula, data, contX, 
    dX = m_rpartD(y = eval(formula[[2L]], envir = data), X = data[[contX]])
) {
  
  nms <- make.names(dimnames(dX)[[2L]])
  data[nms] <- dX
  
  yval <- eval(formula[[2L]], envir = data)
  if (anyNA(yval)) stop('do not allow missingness in the response, for now')
  
  new_fom <- formula
  for (i in nms) new_fom[[3]] <- call('+', new_fom[[3]], as.symbol(i))
  
  suppressWarnings(mod <- if (inherits(yval, what = 'Surv')) {
    coxph(formula = new_fom, data = data)
  } else if (is.logical(yval) || all(yval %in% c(0, 1))) {
    glm(formula = new_fom, data = data, family = binomial(link = 'logit'))
  } else {
    lm(formula = new_fom, data = data)
  })
  
  # if (anyNA(mod$coefficients)) stop('why having NA coefficient estimates?') # future work
  coef_ <- tail(mod$coefficients, n = length(nms))
  attr(coef_, which = 'branch') <- attr(dX, which = 'branch', exact = TRUE)
  attr(coef_, which = 'model') <- mod
  return(coef_)
  
}




