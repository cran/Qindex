
#' @title Predicted Sign-Adjusted Quantile Indices
#' 
#' @description
#' To predict sign-adjusted quantile indices of a test set.
#' 
#' @param object an \linkS4class{Qindex} object based on the training set.
#' 
#' @param newdata test \link[base]{data.frame}, with at least 
#' the response \eqn{y^{\text{new}}} and
#' the \link[base]{double} \link[base]{matrix} of 
#' functional predictor values \eqn{X^{\text{new}}}
#' of the test set, tabulated on the same \eqn{x}-grid as the training set \eqn{X}.
#' If missing, the training set `object@gam$data` will be used.
#' 
#' @param ... additional parameters, currently not in use.
#' 
#' @details 
#' 
#' Function [predict.Qindex] computes 
#' the predicted sign-adjusted quantile indices on the test set, 
#' which is 
#' the product of function \link[mgcv]{predict.gam} return
#' and the correlation sign based on training set
#' (`object@@sign`, see Step 3 of section **Details** of function [Qindex]).
#' Multiplication by `object@@sign` is required to ensure
#' that the predicted sign-adjusted quantile indices
#' are positively associated with the **training** functional predictor values
#' at the selected tabulating grid.
#' 
#' 
#' @returns 
#' Function [predict.Qindex] returns a 
#' \link[base]{double} \link[base]{vector}, 
#' which is the predicted sign-adjusted quantile indices on the test set.
#' 
#' @importFrom mgcv predict.gam
#' @importFrom stats predict
#' @export predict.Qindex
#' @export
predict.Qindex <- function(
    object, 
    newdata = object@gam$data,
    ...
) {
  
  rhs <- object@formula[[3L]]
  olddata <- object@gam$data
  X <- olddata[[rhs]] 
  xgrid <- as.double(colnames(X))
  
  newX <- newdata[[rhs]]
  if (!is.matrix(newX)) stop('`newdata` does not contain a matrix column of functional predictor values')
  new_xgrid <- as.double(colnames(newX))
  if (!all.equal.numeric(new_xgrid, xgrid)) stop('grid of training and test data must be exactly the same')
  
  dm <- dim(newX)
  newdata$xgrid_ <- tcrossprod(rep(1, times = dm[1L]), xgrid)
  newdata$L <- array(1/length(xgrid), dim = dm)
  
  fv <- predict.gam(object = object@gam, newdata = newdata)
  # fitted value `fv` is 'array'
  return(as.double(fv) * object@sign)
  #?base::as.double much faster than ?base::c
}






