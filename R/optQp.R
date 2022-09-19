

## output# (i) medianHR: data set with Np rows including median abs(logHR), corresponding HR
# (ii) thresholds: data set with all cut points for all candidate quantiles and Np columns and Nsplit rows
# (iii) medianHR.Ntop: short list of Ntop best quantile predictors
# (iv) thresholds.Ntop: short list of Ntop median cut points

#' @title Optimal Quantile Predictor 
#' 
#' @description 
#' From the sequence of sample quantiles, selects the optimal quantile that has highest predictive value for a given survival outcome 
#' 
#' @param formula \link[stats]{formula}, currently only supports one \link[survival]{Surv} endpoint 
#' and one predictor of the markers.
#' See details of parameter \code{data}.
#' 
#' @param data \link[base]{data.frame}, with at least 
#' \itemize{
#' \item one \link[survival]{Surv} column as the survival outcome
#' \item one \link[base]{matrix} column as the predictor of quantile sequence
#' }
#' 
#' @param seeds \link[base]{integer} vector of random seeds for generating repeated samples, see \link[base]{set.seed}
#' 
#' @param pct_train \link[base]{numeric} scalar, proportion of the training set, default \code{.8}
#' 
#' @param minbucket.p (for \link[rpart]{rpart} function) the minimum number of observations in any terminal \strong{leaf} node. If only
#'  one of \code{minbucket} or \code{minsplit} is specified, the code either sets \code{minsplit} to
#'  \code{minbucket*3} or \code{minbucket} to \code{minsplit/3}, as appropriate.
#' 
#' @return 
#' \link{evalQp} returns a \link{evalQp} object, which is a \link[base]{list} with two \link[base]{matrix} elements
#' \code{thresholds} and \code{coefs}.
#' 
#' \itemize{
#' \item {\code{thresholds} is ...}
#' \item {\code{coefs} is ...}
#' }
#' 
#' 
#' @examples 
#' # see ?optQp
#' 
#' 
#' @export
evalQp <- function(formula, data, seeds, pct_train = .8, minbucket.p = 0.5) {
  
  y <- formula[[2L]]
  marker <- formula[[3L]]
  
  if (anyNA(data[[y]])) stop('do not allow missingness in the response; remove manually')
  if (anyNA(data[[marker]])) stop('do not allow missingness in the marker predictor; remove manually')
  pseq <- colnames(data[[marker]])
  
  #a. Split data by event
  id_obs <- (data[[y]][,2L] == 1) # 'logical' indexes of observed events
  idx_obs <- which(id_obs) # 'integer' indexes of observed events
  idx_cens <- which(!id_obs) # 'integer' indexes of censored events
  n_obs <- sum(id_obs)
  n <- length(id_obs) # total number of patients
  nq <- dim(data[[marker]])[2L] # number of quantiles
  
  # set default threshold as 0
  # set default coefficient as 0
  coefs <- thresholds <- array(0, dim = c(length(seeds), nq), dimnames = list(NULL, pseq))
  
  # random splits
  for (k in seq_along(seeds)) { #k=1
    
    #b.Randomly shuffle the data
    set.seed(seed = seeds[k])
    idx_obs_shuffled <- sample(idx_obs)
    idx_cens_shuffled <- sample(idx_cens)
    
    folds_obs <- cut.default(seq(1,n_obs), breaks=1/(1-pct_train))
    folds_cens <- cut.default(seq(1,n - n_obs),breaks=1/(1-pct_train))
    # Tingting: really wierd programming!!!
    
    #d. Model fitting and evaluation
    #Creating empty object to hold fit information
    
    #firstQp = dim(data)[2]-length(pseq)+1
    #lastQp = dim(data)[2]
    
    #Segment data by fold using the which() function 
    
    testID_obs <- (unclass(folds_obs) == 1L)
    testID_cens <- (unclass(folds_cens) == 1L)
    dim(testData <- data[c(idx_obs_shuffled[testID_obs], idx_cens_shuffled[testID_cens]), ])
    dim(trainData <- data[c(idx_obs_shuffled[!testID_obs], idx_cens_shuffled[!testID_cens]), ])
    
    minbucket.value <- quantile(trainData[[y]][,1L][trainData[[y]][,2L]==1], probs = minbucket.p)
    
    #Use the test and train data partitions
    #Model fitting and evaluation
    for (p in 1:nq) { #p=1
      stree <- rpart(formula = trainData[[y]] ~ trainData[[marker]][,p], control = rpart.control(minbucket = minbucket.value))
      if (length(stree$splits)) thresholds[k,p] <- stree$splits[1L, 4L]
      
      #evaluating fit on the test fold
      testData$high <- (testData[[marker]][,p] > thresholds[k,p])
      if (!any(testData$high) || all(testData$high)) next
      #suppressWarnings(cox_test <- coxph(formula = eval(call('~', y, quote(high))), data = testData))
      cox_test <- coxph(formula = eval(call('~', y, quote(high))), data = testData)
      cf_test <- cox_test$coefficients
      if (!is.na(cf_test) && (-2 < cf_test) && (cf_test < 2)) coefs[k,p] <- cf_test # else do nothing
    }#end of p loop
    
  }#end of k loop
  
  ret <- list(
    thresholds = thresholds,
    coefs = coefs,
    data = data
  )
  class(ret) <- 'evalQp'
  return(ret)
  
}


# write some document to say `FUN` can only be stats::median or base::mean

#' @title Summary Information of \link{evalQp} Object
#' 
#' @description Summary method for \link{optQp} Object
#' 
#' @param object \link{evalQp} object
#' 
#' @param FUN \link[stats]{median} (default) or \link[base]{mean} 
#' 
#' @param ... ..
#' 
#' @return 
#' \link{summary.evalQp} summarize ...
#' 
#' 
# ?base::summary
#' @export
summary.evalQp <- function(object, FUN = median, ...) {
  ret <- data.frame(
    abs.logHR = abs(apply(object$coefs, MARGIN = 2L, FUN = FUN)),
    HR = exp(apply(object$coefs, MARGIN = 2L, FUN = FUN)),
    threshold = apply(object$thresholds, MARGIN = 2L, FUN = FUN)
  )
  ret[order(ret$abs.logHR, decreasing = TRUE), ]
}


#' @title Print \link{evalQp} Object
#' 
#' @description ..
#' 
#' @param x ..
#' 
#' @param n ..
#' 
#' @param ... ..
#' 
#' @details 
#' \link{print.evalQp} simply calls \link{summary.evalQp} ..
#' 
#' @return 
#' \link{print.evalQp} does not have a returned value.
#' 
#' @export
print.evalQp <- function(x, n = 3L, ...) {
  ret <- summary.evalQp(x)
  cat(sprintf('by median, first %d rows\n', n))
  print(head(ret, n = n))
  return(invisible())
}




#' @title \link{optQp}
#' 
#' @description \link{optQp}
#' 
#' @param x \link{evalQp} object
#' 
#' @param n \link[base]{integer} scalar
#' 
#' @param ... additional parameters, not currently in use
#' 
#' @return 
#' \link{optQp} returns ..
#' 
#' @examples 
#' # see ?BBC_dichotom
#' 
#' @export
optQp <- function(x, n = 3L, ...) {
  if (!inherits(x, 'evalQp')) stop('input must be evalQp object')
  tmp <- summary.evalQp(x)
  head(rownames(tmp), n = n)
  ret <- x$data
  ret$Marker <- ret$Marker[, head(rownames(tmp), n = n)]
  return(ret)
}

