
#' @title Dichotomizing Split Sample
#' 
#' @description
#' 
#' For a given predictor and a random sample split, 
#' finds the optimal dichotomizing branch from the training set,
#' applies this branch to dichotomize the predictor in the test set,
#' and runs a regression analysis in the test set using the dichotomized predictor.
#' 
#' @param y a \link[survival]{Surv} object, a \link[base]{logical} \link[base]{vector}, 
#' or a \link[base]{numeric} \link[base]{vector}, the endpoint
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, the predictor
#' 
#' @param splitID a named length-2 \link[base]{list} of \link[base]{integer} \link[base]{vector}s,
#' the indexes of \code{'train'} and \code{'test'} set. 
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' Given a sample split in \code{splitID},
#' 
#' \enumerate{
#' 
#' \item finds the optimal dichotomizing branch of endpoint \code{y} 
#' using predictor \code{x}, via \link{rpartD}, in the training set
#' 
#' \item dichotomizes the predictor \code{x} using the branch identified in Step 1, 
#' in the test set.
#' The dichotomized \code{x} is denoted as \code{TRUE} 
#' if \code{x} is greater than, or greater-than-or-equal-to, the threshold
#' (see \link{rpartD}).
#' 
#' \item runs a regression model on the endpoint \code{y} 
#' using the dichotomized predictor \code{x} in Step 2, in the test set.
#' Currently the Cox proportional hazards (\link[survival]{coxph}) regression for \link[survival]{Surv} endpoint, 
#' logistic (\link[stats]{glm}) regression for \link[base]{logical} endpoint and 
#' linear (\link[stats]{lm}) regression for \link[stats]{gaussian} endpoint
#' are supported.
#' 
#' }
#' 
#' @note
#' \link[rpart]{rpart} in \link{rpartD} is the greatest computation cost in \link{splitSample_dichotom}.
#' 
#' @returns 
#' 
#' \link{splitSample_dichotom} returns an object of S3 class \link{splitSample_dichotom},
#' which is essentially
#' a Cox proportional hazards (\link[survival]{coxph}), 
#' or a logistic (\link[stats]{glm}), 
#' or a linear (\link[stats]{lm}) 
#' regression model, 
#' with three (3) additional \link[base]{attributes}
#' 
#' \describe{
#' 
#' \item{\code{attr(,'branch')}}{\link[base]{language} object, 
#' the branch identified in the training set, via \link{rpartD}}
#' 
#' \item{\code{attr(,'highX')}}{\link[base]{numeric} scalar, 
#' percentage of \code{x}, in the test set, which is greater than, or greater-than-or-equal-to, 
#' the threshold represented in \code{attr(, 'branch')}}
#' 
#' \item{\code{attr(,'coef')}}{\link[base]{numeric} scalar, 
#' the regression coefficient corresponding to the dichotomized \code{x} in the test set}
#' 
#' }
#' 
#' 
#' @importFrom survival coxph
#' @importFrom stats lm glm binomial
#' @export
splitSample_dichotom <- function(y, x, splitID, ...) {
  
  id_train <- splitID[['train']]
  id_test <- splitID[['test']]
  
  branch <- branch_ <- rpartD(y = y[id_train], x = x[id_train])
  branch_[[2L]] <- quote(x[id_test])
  
  # then use the branch to dichotomize the predictors of the test data 
  dtest <- data.frame(y = y[id_test], high = eval(branch_))
  
  if (!any(dtest$high) || all(dtest$high)) {
    # exception
    mtest <- logical()
    attr(mtest, which = 'coef') <- NA_real_
    return(mtest)
  }
  
  mtest <- if (inherits(y, what = 'Surv')) {
    suppressWarnings(coxph(formula = y ~ high, data = dtest))
  } else if (is.logical(y) || all(y %in% c(0, 1))) {
    suppressWarnings(glm(formula = y ~ high, family = binomial(link = 'logit'), data = dtest))
  } else if (is.vector(y, mode = 'numeric')) {
    suppressWarnings(lm(formula = y ~ high, data = dtest))
  }
  
  cf_test <- mtest$coefficients[length(mtest$coefficients)]
  attr(mtest, which = 'coef') <- if (is.finite(cf_test)) unname(cf_test) else NA_real_
  attr(mtest, which = 'highX') <- mean.default(dtest$high)
  attr(mtest, which = 'branch') <- branch
  class(mtest) <- c('splitSample_dichotom', class(mtest))
  return(mtest)
  
}










#' @title Median Effect Size for Multiple Sample Splits
#' 
#' @description 
#' Defining the median of multiple \link{splitSample_dichotom} objects, 
#' as the \link{splitSample_dichotom} object with 
#' the median of the regression coefficient of dichotomized predictor (\code{attr(,'coef')}).
#' 
#' @param y a \link[survival]{Surv} object, a \link[base]{logical} \link[base]{vector}, 
#' or a \link[base]{numeric} \link[base]{vector}, the endpoint
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, the predictor
#' 
#' @param splitIDs (optional) a \link[base]{list}.
#' Each element of \code{splitIDs} is a length-2 \link[base]{list} of \link[base]{integer} \link[base]{vector}s,
#' the indexes of \code{'train'} and \code{'test'} set. 
#' Default value is a *series of* sample splits (via \link{stratifiedSplitSample}) 
#' based on the endpoint \code{y}.
#' 
#' @param ... additional parameters for \link{stratifiedSplitSample} if \code{splitIDs} is missing, 
#' most importantly the copies of sample splits \code{nsplit}
#' 
#' @details
#' 
#' Given a *series of* sample splits \code{splitIDs}, 
#' 
#' \enumerate{
#' 
#' \item {for each sample split, use \link{splitSample_dichotom}}
#' 
#' \item {finds the nearest-even median (\code{type = 3} of \link[stats]{quantile})
#' of the regression coefficient of dichotomized predictor \code{attr(,'coef')},
#' from all \link{splitSample_dichotom} objects obtained in Step 1}
#' 
#' \item {returns the \link{splitSample_dichotom} object 
#' with the median regression coefficient of dichotomized predictor, identified in Step 2}
#' 
#' }
#' 
#' 
#' @returns 
#' \link{median_splitSample_dichotom} returns an object of S3 class \link{splitSample_dichotom}.
#' 
#' @importFrom stats quantile
#' @export
median_splitSample_dichotom <- function(y, x, splitIDs = stratifiedSplitSample(y, ...), ...) {
  
  tmp <- lapply(splitIDs, FUN = function(splitID) splitSample_dichotom(y, x, splitID = splitID, ...))
  
  cf <- vapply(tmp, FUN = attr, which = 'coef', exact = TRUE, FUN.VALUE = NA_real_)
  
  medianID <- which(cf == quantile(cf, probs = .5, type = 3L, na.rm = TRUE))[1L]
  return(tmp[[medianID]])
  
}






