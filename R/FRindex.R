
#' @title Functional Regression Index (FRindex)
#' 
#' @description 
#' Computes a scalar functional index as a predictor in the functional regression model.
#' The training data are used to estimate the functional coefficient by fitting functional regression model 
#' with any response supported by \CRANpkg{gam} package. 
#' The functional regression index is computed for training and (optional) test data 
#' using the functional coefficient from the model fitted to the training data.
#' 
#' @param trainData \link[base]{data.frame}, training data including 
#' a response variable and a \link[base]{matrix} of tabulated functional predictor.
#' If the functional predictor is the \link[stats]{quantile} function,
#' then \code{trainData} is preferably the returned object of \link{sampleQp}.
#' 
#' @param testData \link[base]{data.frame}, test data including 
#' a response variable and a \link[base]{matrix} of tabulated functional predictor.
#' The number of values in the grid for tabulating functional predictor 
#' does not need to be the same as that number for the training data.
#' If the functional predictor is the \link[stats]{quantile} function,
#' then \code{testData} is preferably the returned object of \link{sampleQp}.
#' 
#' @param trainArgs \link[base]{numeric} \link[base]{vector}
#' 
#' @param testArgs \link[base]{numeric} \link[base]{vector}
#' 
#' @param response \link[base]{language} of the response variable.
#' All response variable types supported by \CRANpkg{gam} package are allowed, 
#' including continuous, binary or survival outcome.
#' 
#' @param family \link[stats]{family} object specifying the distribution 
#' and link to use in \link[mgcv]{gam}
#' 
#' @param log \link[base]{logical} scalar, whether to perform \link[base]{log} transformation on quantiles (default \code{FALSE}).
#'  
#' @param predictor \link[base]{character} scalar, 
#' name of the \link[base]{matrix} column variable in \code{trainData} and \code{testData}, 
#' with each row representing the tabulated, on a common grid, functional predictor values for each subject.
#' If the functional predictor is the \link[stats]{quantile} function,
#' then the number of columns in the \link[base]{matrix} column identified by \code{predictor}
#' equals to the number of quantiles used.
#' 
#' @param knot_pct \link[base]{numeric} scalar, 
#' percentage of the column dimension of the \link[base]{matrix} column \code{predictor},
#' to be used as \code{knot.value}.  
#' Default is \eqn{40\%}.
#' If \code{knot.value} is provided by the end-user, then \code{knot_pct} is not used.
#' 
#' @param knot.value \link[base]{integer} scalar, number of knots 
#' (i.e., parameter \code{k} in the spline smooth function \link[mgcv]{s})
#' used in \link[mgcv]{gam} function.
#' Default is the \link[base]{ceiling} of \code{knot_pct} of
#' the column dimension of the \link[base]{matrix} column \code{predictor}.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' The Functional Regression Index (\link{FRindex}) is defined
#' as the integral of the functional predictor multiplied by common weight function.
#' The weight function is proportional to the functional coefficient in the functional
#' regression model for a given response variable fitted to training data.
#' The \link{FRindex} values computed for the training data and test data, if provided. 
#' 
#' @references 
#' 
# \url{https://link.springer.com/book/10.1007/b98888#about-this-book}
#' J. O. Ramsay, B. W. Silverman (2005). Functional Data Analysis, ed 2. Springer New York, NY
#' \doi{10.1007/b98888}  
#' 
#' Cui, E., Crainiceanu, C. M., & Leroux, A. (2021). Additive Functional Cox Model. Journal of Computational and Graphical Statistics. 2021;30(3):780-793.
#' \doi{10.1080/10618600.2020.1853550}
#' 
#' Gellar, J. E., Colantuoni, E., Needham, D. M., & Crainiceanu, C. M. (2015). Cox regression models with functional covariates for survival data. Statistical Modelling, 15(3), 256-278.
#' \doi{10.1177/1471082X14565526}
#' 
#' @return 
#' \link{FRindex} returns a \link[base]{data.frame} with the training data and FRindex, 
#' and the vectors of weights on the grid of the functional predictor
#' 
#' 
#' @examples 
#' pt = unique(Ki67$PATIENT_ID)
#' length(pt) # 622
#' train = subset(Ki67, PATIENT_ID %in% pt[1:500])
#' test = subset(Ki67, PATIENT_ID %in% pt[501:622])
#' train_Qps = sampleQp(data = train, subjID = 'PATIENT_ID', 
#'   exclude = c('tissueID','inner_x','inner_y'), Qpredictor = 'Marker')
#' test_Qps = sampleQp(data = test, subjID = 'PATIENT_ID',
#'   exclude = c('tissueID','inner_x','inner_y'), Qpredictor = 'Marker')
#' 
#' if (FALSE) { # masked to save time
#' FRQI = FRindex(trainData = train_Qps, testData = test_Qps, 
#'   response = Surv(RECFREESURV_MO, RECURRENCE), predictor = 'Marker', log = TRUE)
#' set.seed(1)
#' mod = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, data = FRQI$trainData, 
#'   contX = 'FRindex_std', R = 100)
#' }
#' 
#' if (FALSE) { # masked to save time
#' FRindex(trainData = train_Qps, testData = test_Qps, response = RECURRENCE)
#' FRindex(trainData = train_Qps, testData = test_Qps, response = RECFREESURV_MO)
#' 
#' set.seed(1)
#' mod = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ 1, data = FRQI$trainData, 
#'   contX = 'FRindex_std', R = 200)
#' summary(mod)
#' names(attributes(mod))
#' attr(mod, 'median_optimism')
#' }
#' 
#' @export
FRindex <- function(
    trainData, testData,
    trainArgs = attr(trainData, which = 'xarg', exact = TRUE),
    testArgs = if (!missing(testData)) attr(testData, which = 'xarg', exact = TRUE) else trainArgs,
    response, family,
    predictor = 'Marker', 
    log = TRUE,
    knot_pct = .4,
    knot.value = ceiling(ncol(trainData[[predictor]]) * knot_pct), 
    ...
) {
  
  #r_pred
  #s_pred
  
  if (!length(trainArgs) || !is.numeric(trainArgs) || anyNA(trainArgs)) stop('`trainArgs` needs to be a numeric sequence')
  if (!length(testArgs) || !is.numeric(testArgs) || anyNA(testArgs)) stop('`testArgs` needs to be a numeric sequence')
  if (!isTRUE(all.equal.numeric(min(trainArgs), min(testArgs))) ||
      !isTRUE(all.equal.numeric(max(trainArgs), max(testArgs)))) stop('Training and test data should share the same domain')

  train_pred <- trainData[[predictor]]
  test_pred <- if (missing(testData)) train_pred else testData[[predictor]]
  
  train_pred <- if (log) log1p(train_pred) else train_pred
  test_pred <- if (missing(testData)) train_pred else if (log) log1p(test_pred) else test_pred
  
  train_dm <- dim(train_pred)
  # stopifnot(train_dm[2L] == length(trainArgs))
  test_dm <- dim(test_pred)
  # stopifnot(test_dm[2L] == length(testArgs))
  
  L <- array(1/train_dm[2L], dim = train_dm)
  
  ### `Tm`: time indices, we assume an equally-spaced grid of (0, 1) for example
  Tm <- tcrossprod(rep(1, times = train_dm[1L]), trainArgs)
  
  response <- substitute(response)
  s_term <- quote(s(Tm, by = L * train_pred, bs = 'cr', k = knot.value))
  train_gam <- if (is.call(response) && (response[[1L]] == 'Surv')) {
    family <- if (missing(family)) quote(cox.ph()) else substitute(family)
    eval(call('gam', formula = call('~', response[[2L]], s_term),
         weights = response[[3L]], data = quote(trainData), family = family))
  } else if (is.name(response)) {
    yval <- trainData[[response]]
    family <- if (!missing(family)) {
      substitute(family) 
    } else if (is.logical(yval) || all(yval %in% c(0, 1))) {
      quote(binomial(link = 'logit'))
    } else if (is.numeric(yval)) {
      quote(gaussian(link = 'identity'))
    } else stop('not supported yet')
    eval(call('gam', formula = call('~', response, s_term), data = quote(trainData), family = family))
  } else stop('illegal `response`; not supported yet')
  
  # extract beta function from ?mgcv::gam return, using ?mgcv::plot.gam
  # suppress the figure from printing; inspired by ?oddsratio::no_plot
  png(tmp_figure <- file.path(tempdir(), 'tmp.png'))
  train_beta <- c(plot.gam(train_gam, n = train_dm[2L])[[1]]$fit) # length train_pred
  test_beta <- if (missing(testData)) train_beta else c(plot.gam(train_gam, n = test_dm[2L])[[1]]$fit) # length `test_pred`
  dev.off()
  file.remove(tmp_figure)
  
  # ?stats::integrate: might be more difficult than Tingting originally thought..
  tmp <- train_pred * tcrossprod(rep(1, times = train_dm[1L]), train_beta) # 500 * 99
  
  #train_FRi <- ((1-0)/(train_dm[2L]+1L)) * (train_dm[2L] * rowMeans(tmp) - tmp[,1L]/2 - tmp[,train_dm[2L]]/2) # Newton integral
  train_FRi <- ((max(trainArgs)-min(trainArgs))/(train_dm[2L]+1L)) * (train_dm[2L] * rowMeans(tmp) - tmp[,1L]/2 - tmp[,train_dm[2L]]/2) # Newton integral
  
  # 2023-02-07: we are thinking about two questions
  # 1. if we could change the domain from (0, 1) to an arbitrary (x1, x2)
  # 2. how do we deal with the two extreme values of the domain
  # To solve these problems, we may add two more parameters
  # @param domain \link[base]{numeric} vector of length-2, with default value c(0, 1)
  # @param extreme \link[base]{logical} scalar, whether to use the two extreme values of the domain. For now, we think \code{(x1, x2]} is not meaningful
  
  ### `Qp[,(nQ + 1L)/2]` a vector of the medians of markers from all patients in training data
  # `FRi` the integral from all patients in training data
  train_sign_wts <- sign(cor(train_pred[,(train_dm[2L] + 1L)/2], train_FRi)) # scalar
  
  ## ADJUST SIGN: multiply by SIGN OF CORRELATION of FRi and the median
  ## This does nothing if Q50 and FRi is positively correlated, otherwise it flips the sign of all FRi's, which reciprocates all HR for continuous FRi's.
  train_idx <- train_sign_wts * train_FRi # 500
  
  train_qs <- quantile(train_idx, probs = c(.25, .5, .75))
  train_idx_std <- (train_idx - train_qs[2L]) / (train_qs[3L] - train_qs[1L])
  
  train_idx_wts <- train_sign_wts * train_beta # 99
  test_idx_wts <- train_sign_wts * test_beta # 19
  
  tmp <- test_pred * tcrossprod(rep(1, times = test_dm[1L]), test_idx_wts) # 122 * 19
  
  #test_idx <- ((1-0)/(test_dm[2L]+1L)) * (test_dm[2L] * rowMeans(tmp) - tmp[,1L]/2 - tmp[,test_dm[2L]]/2) # Newton integral
  test_idx <- ((max(testArgs)-min(testArgs))/(test_dm[2L]+1L)) * (test_dm[2L] * rowMeans(tmp) - tmp[,1L]/2 - tmp[,test_dm[2L]]/2) # Newton integral
  
  test_qs <- quantile(test_idx, probs = c(.25, .5, .75))
  test_idx_std <- (test_idx - test_qs[2L]) / (test_qs[3L] - test_qs[1L])
  
  return(list(
    trainData = data.frame(trainData, FRindex = train_idx, FRindex_std = train_idx_std), 
    testData = if (!missing(testData)) data.frame(testData, FRindex = test_idx, FRindex_std = test_idx_std), # else NULL
    trainArgs = trainArgs,
    testArgs = if (!missing(testData)) testArgs, # else NULL
    train_sign_wts = train_sign_wts, 
    train_idx_wts = train_idx_wts, 
    #test_idx_wts = if (!missing(testData) && (train_dm[2L] != test_dm[2L])) test_idx_wts, # else NULL # KEEP THIS LINE
    test_idx_wts = if (!missing(testData)) test_idx_wts, # else NULL
    gam = train_gam
  ))
}


