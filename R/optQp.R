
#' @title Optimal Quantile Predictor
#' 
#' @description
#' From a given set of sample quantiles, this function selects the optimal quantile 
#' with the largest effect size for predicting given \link[survival]{Surv}, 
#' \link[base]{logical}, or \link[base]{numeric} outcome 
#' 
#' @param formula \link[stats]{formula}, supports \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint 
#' and one \link[base]{matrix} predictor of the descriptive statistics of markers per cluster.
#' See details of parameter \code{data}.
#' 
#' @param data \link[base]{data.frame}, with at least 
#' \itemize{
#' \item {two \link[survival]{Surv} columns including time-to-event and event indicator,
#'        one \link[base]{logical} column, 
#'        or one \link[base]{numeric} column as the outcome}
#' \item {one \link[base]{matrix} column as the descriptive statistics of markers per cluster.  
#' Currently only a quantile sequence is supported.}
#' }
#' 
#' @param seeds \link[base]{integer} vector of random seeds for generating repeated random split samples, see \link[base]{set.seed}
#' 
#' @param pct_train \link[base]{numeric} scalar, proportion of the training set, default \code{.8}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' Optimal quantile selection algorithm
#' For a sample \eqn{x_i ,1 \leq i \leq n} of repeated measures of independent variable \eqn{X}
#' and probability \eqn{p, 0 < p < 1}, 
#' the empirical quantile function \eqn{Q_n(p)} is defined as the k th order statistic of the sample,
#' where \eqn{k} is such that \eqn{(k-1)/n < p < k/n}. For \eqn{p = .01,\cdots,.99} and \eqn{k = 100p},
#' \eqn{Q_n(p)} is also known as \eqn{k}th percentile.
#' The following algorithm is proposed to identify the optimal \eqn{Q(p)} predictor of
#' survival outcome in a screening data set:
#' 
#' \enumerate{
#' 
#' \item Select the set of quantiles to be evaluated as predictors and the desired ratio
#' for training/test sets.
#' 
#' \item Split the data into training and test sets
#' \itemize{
#' \item In case of \link[survival]{Surv} outcomes, split the group of subjects with event into a training set and a test set randomly
#' with desired ratio. Similarly, split the group of subjects without event into a training
#' set and a test set randomly with desired ratio. Combine training sets with and
#' without event and test sets with and without event.
#' \item In case of \link[base]{logical} outcomes, split the group of subjects with one level of outcome into a training set and a test set randomly
#' with desired ratio. Similarly, split the group of subjects with the other level of outcome into a training
#' set and a test set randomly with desired ratio. Combine training sets with both levels of logical outcome and test sets with both levels of logical outcome.
#' \item In case of \link[base]{numeric} outcomes, split the entire subjects into a training set and a test set randomly
#' with desired ratio. 
#' }
#' 
#' \item For each training/test set pair and each considered quantile,
#' \itemize{
#' \item Determine the optimal cutoff (e.g. using the R package rpart) in the combined
#' training set.
#' \item Apply the optimal cutoff to the combined test set and estimate the effect size
#' (hazard ratio, odds ratio, or exponentiated coefficient).
#' }
#' 
#' \item Repeat steps 2 and 3 for 100 training/test splits, compute the median log effect size
#' (log hazard ratio(HR), log odds ratio(OR), or coefficient) for each quantile.
#' 
#' \item Rank the effect sizes (absolute value of log OR, HR, or coefficient) for all considered quantiles and
#' select the optimal quantile with the highest effect size.
#' 
#' \item Perform bootstrap-based optimism correction for the selected optimal quantile(s).
#' In this work, for each random split, \eqn{80\%} of subjects were assigned to the training
#' set and \eqn{20\%} of subjects were assigned to the test set. We considered every fifth
#' quantile starting from the 5th quantile to the 95th quantile plus 99th quantile as
#' candidate predictors of PFS (a total of 20 quantile predictors). Also, we identified
#' quantiles with the second and third highest effect sizes to compare them to the
#' optimal ones.
#' }
#' 
#' @return 
#' \link{eval_dichotom} returns a \link{eval_dichotom} object, which is a \link[base]{list} with elements  elements
#' \code{thresholds} and \code{coefs}.
#' \describe{
#' \item{\code{thresholds}}{\link[base]{matrix}, data set with all cut points for all candidate quantiles and Np columns and Nsplit rows }
#' \item{\code{coefs}}{\link[base]{matrix}, coefficients corresponding to thresholds}
#' \item{\code{data}}{\link[base]{data.frame}, this is the unmodified input \code{data}}
#' }
#' 
#' 
#' @references 
#' Selection of optimal quantile protein biomarkers based on cell-level immunohistochemistry data,
#' Misung Yi, Tingting Zhan , Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva, under review
#' 
#' @examples 
#' 
#' if (FALSE) { # masked to save time
#' library(survival)
#' 
#' Ki67_Qps = sampleQp(data = Ki67, subjID = 'PATIENT_ID', 
#'   exclude = c('tissueID','inner_x','inner_y'), Qpredictor = 'Marker')
#' Ki67c = eval_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = Ki67_Qps, 
#'   seeds = 1:20)
#' # summary(Ki67c) # works, but not needed
#' head(Ki67_opt <- optQp(Ki67c, n = 2L))
#' 
#' set.seed(1)
#' mod_c = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, 
#'   data = optQp(Ki67c, n = 2L), contX = 'Marker', R = 100)
#' summary(mod_c)
#' }
#' 
#' if (FALSE) { # mask to save time
#' Ki67a = eval_dichotom(RECFREESURV_MO ~ Marker, data = Ki67_Qps, seeds = 1:20)
#' set.seed(1)
#' mod_a = BBC_dichotom(RECFREESURV_MO ~ NodeSt + Tstage, data = optQp(Ki67a, n = 2L), 
#'   contX = 'Marker', R = 100)
#' summary(mod_a)
#' 
#' Ki67b = eval_dichotom(RECURRENCE ~ Marker, data = Ki67_Qps, seeds = 1:20)
#' set.seed(1)
#' mod_b = BBC_dichotom(RECURRENCE ~ NodeSt + Tstage, data = optQp(Ki67b, n = 2L), 
#'   contX = 'Marker', R = 100)
#' summary(mod_b)
#' }
#' 
#' @export
eval_dichotom <- function(formula, data, seeds, pct_train = .8, ...) {
  
  y <- formula[[2L]] # endpoint (in language)
  yval <- eval(y, envir = data) # this was programmed for 'Surv' value.  Now to be generalized to *any type* of endpoint
  if (anyNA(yval)) stop('do not allow missingness in the response, for now')
  
  marker <- formula[[3L]] # marker (in language)
  markerval <- eval(marker, envir = data) # 'matrix' of markers.  To be processed one column at a time
  if (anyNA(markerval)) stop('do not allow missingness in the \'marker\' matrix; for now')
  pseq <- colnames(markerval)
  nq <- dim(markerval)[2L] # number of quantiles
  n <- dim(markerval)[1L] # total number of subjects
  
  if (inherits(yval, what = 'Surv')) {
    # Split into training/test, stratified by censoring status
    if (y[[1L]] != 'Surv') stop('endpoint must be specified like `Surv(time, event)`')
    yevent <- as.logical(eval(y[[3L]], envir = data))
    idx1 <- which(yevent) # 'integer' indexes of observed events
    idx0 <- which(!yevent) # 'integer' indexes of censored events
    n1 <- sum(yevent)
    
  } else if (is.logical(yval) || all(yval %in% c(0, 1))) {
    # Split into training/test, stratified by the binary response
    yval <- as.logical(yval)
    idx1 <- which(yval) # 'integer' indexes of response
    idx0 <- which(!yval) # 'integer' indexes of non-response
    n1 <- sum(yval)
    
  } else {
    # Split into training/test, without stratification
    idx1 <- seq_len(n)
    idx0 <- integer()
    n1 <- n
  }
  
  coefs <- thresholds <- array(NA_real_, dim = c(length(seeds), nq), dimnames = list(NULL, pseq))
  
  # random splits
  for (k in seq_along(seeds)) { #k=1
    
    set.seed(seed = seeds[k])
    
    idx1_train <- sort.int(sample(idx1, size = floor(length(idx1) * pct_train), replace = FALSE))
    idx1_test <- setdiff(idx1, idx1_train)
    idx0_train <- sort.int(sample(idx0, size = floor(length(idx0) * pct_train), replace = FALSE))
    idx0_test <- setdiff(idx0, idx0_train)
    
    trainData <- data[c(idx1_train, idx0_train), ]
    testData <- data[c(idx1_test, idx0_test), ]
    
    for (p in 1:nq) { # for each column of 'marker'
      tmp_train <- data.frame(y = eval(y, envir = trainData), x = trainData[[marker]][,p])
      train_tree <- rpart(formula = y ~ x, data = tmp_train, 
                          # method = ???,
                          control = rpart.control(
                            cp = .Machine$double.eps, # to force a split even if the overall lack of fit is not decreased
                            maxdepth = 2L
                          ))
      if (!length(train_tree$splits)) stop('we must force a split')
      thresholds[k,p] <- train_tree$splits[1L, 4L]
        
      testData$high <- (testData[[marker]][,p] > thresholds[k,p]) # a new binary column, indicating if the selected marker is higher than the cutoff
      if (!any(testData$high) || all(testData$high)) next
      test_model <- if (inherits(yval, what = 'Surv')) {
        suppressWarnings(coxph(formula = eval(call('~', y, quote(high))), data = testData))
      } else if (is.logical(yval)) {
        suppressWarnings(glm(formula = eval(call('~', y, quote(high))), family = binomial(link = 'logit'), data = testData))
      } else {
        suppressWarnings(lm(formula = eval(call('~', y, quote(high))), data = testData))
      }
      cf_test <- test_model$coefficients[length(test_model$coefficients)]
      #if (is.finite(cf_test) && (abs(cf_test) < 2)) coefs[k,p] <- cf_test # else do nothing
      if (is.finite(cf_test)) coefs[k,p] <- cf_test # else do nothing
    }
    
  }#end of k loop
  
  ret <- list(
    thresholds = thresholds,
    coefs = coefs,
    endpoint = if (inherits(yval, what = 'Surv')) {
      'Surv' 
    } else if (is.logical(yval)) {
      'logical'
    } else 'numeric',
    data = data
  )
  class(ret) <- 'eval_dichotom'
  return(ret)
  
}



#' @title Summary Information of \link{eval_dichotom} Object
#' 
#' @description To summarize \link{eval_dichotom} object
#' 
#' @param object \link{eval_dichotom} object
#' 
#' @param FUN summarizing \link[base]{function}, either \link[stats]{median} (default) or \link[base]{mean} 
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' \link{summary.eval_dichotom} present effect sizes (absolute value of log hazard ratio, log odds ratio, or coefficient), corresponding hazard ratio or odds ratio, and threshold to dichotomize quantile.
#' And then sort by \code{abs.coef}
#' 
#' @return 
#' \link{summary.eval_dichotom} returns a \link[base]{data.frame} with three (3) columns
#' \describe{
#'   \item{\code{abs.coef}}{effect sizes (absolute value of log hazard ratio, log odds ratio, or coefficient)}
#'   \item{\code{HR}}{corresponding hazard ratio or odds ratio}
#'   \item{\code{threshold}}{threshold for dichotomizing the quantile}
#' }
#' 
#' @seealso \link[base]{summary}
#' 
#' @export
summary.eval_dichotom <- function(object, FUN = median, ...) {
  if (!inherits(object, 'eval_dichotom')) stop('input must be `eval_dichotom`')
  ret <- data.frame(coef = apply(object$coefs, MARGIN = 2L, FUN = FUN, na.rm = TRUE))
  ret <- within.data.frame(ret, expr = {
    abs.coef <- abs(coef)
  })
  
  switch(object$endpoint, Surv = {
    ret <- within.data.frame(ret, expr = {
      HazardsRatio <- exp(coef)
    })
  }, logical = {
    ret <- within.data.frame(ret, expr = {
      OddsRatio <- exp(coef)
    })
  }, numeric = {
    # do nothing
  })
  return(ret[order(ret$abs.coef, decreasing = TRUE), ]  )
  
}


#' @title Print \link{eval_dichotom} Object
#' 
#' @description To print the optimal quantile that has the largest absolute value of coefficients  
#' for a given \link[survival]{Surv} outcome, \link[base]{logical}, or \link[base]{numeric} outcome 
#' @param x \link{eval_dichotom} object
#' 
#' @param n \link[base]{integer} scalar, number of optimal quantiles to be printed
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' \link{print.eval_dichotom} simply calls \link{summary.eval_dichotom} and
#'  print three largest absolute value of coefficients, hazard ratio or odds ratio, and threshold
#' 
#' @return 
#' \link{print.eval_dichotom} does not have a returned value.
#' 
#' @export
print.eval_dichotom <- function(x, n = 3L, ...) {
  ret <- summary.eval_dichotom(x)
  cat(sprintf('by median, first %d rows\n', n))
  print(head(ret, n = n))
  return(invisible())
}


#' @title Optimal Quantile Predictor
#' 
#' @description 
#' 
#' Wrapper of \link{eval_dichotom}.
#' 
#' 
#' To create a \link[base]{data.frame} with designated number of optimal quantiles returned from \link{eval_dichotom} 
#' in addition to \code{subjID}, \link[survival]{Surv} endpoint and 
#' \link[base]{numeric}, \link[base]{character}, or \link[base]{logical} predictors
#' 
#' @param x \link{eval_dichotom} object
#' 
#' @param n \link[base]{integer} scalar, number of optimal quantiles to be printed
#' 
#' @param ... additional parameters, not currently in use
#' 
#' @return 
#' \link{optQp} returns a \link[base]{data.frame} with designated number of optimal quantiles
#' 
#' @examples 
#' # see ?BBC_dichotom
#' 
#' @export
optQp <- function(x, n = 3L, ...) {
  if (!inherits(x, what = 'eval_dichotom')) stop('input must be eval_dichotom object')
  tmp <- summary.eval_dichotom(x)
  ret <- x$data
  # to replace the full grid (of probabilities) of 'Marker', by the top-selected markers
  ret$Marker <- ret$Marker[, head(rownames(tmp), n = n)]
  return(ret)
}



# only for developers

# if (FALSE) {
# # eval_dichotom(TStage ~ Marker, data = Ki67_Qps, seeds = 1:20) # next task
# }
# 
# if (FALSE) {
# colnames(Ki67_opt$Marker) = paste0('Ki67_', colnames(Ki67_opt$Marker))
#
# PR_Qps = sampleQp(data = PR, subjID = 'PATIENT_ID', Qpredictor = 'Marker')
# PR_eval = eval_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = PR_Qps, seeds = 1:20)
# PR_opt = optQp(PR_eval, n = 3L)
# colnames(PR_opt$Marker) = paste0('PR_', colnames(PR_opt$Marker))
# 
# cellOpt0 = merge(Ki67_opt, PR_opt, by = setdiff(names(Ki67_opt), 'Marker'),
#   suffixes = c('.Ki67', '.PR'))
# cellOpt = within(cellOpt0, expr = {
#   Marker = cbind(Marker.Ki67, Marker.PR)
#   Marker.Ki67 = Marker.PR = NULL
# })
# dim(cellOpt)
# names(cellOpt)
# head(cellOpt$Marker)
# 
# set.seed(1)
# mod0 = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, data = cellOpt, 
#   contX = 'Marker', R = 100)
# names(mod0)
# }
# 
