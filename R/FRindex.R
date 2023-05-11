
setOldClass('gam')

#' @title \linkS4class{FRindex} Object: Functional Regression Index Model
#' 
#' @description
#' An S4 class for the functional regression index model.
#' 
#' @slot gam \link[mgcv]{gam} object
#' 
#' @slot sign \link[base]{numeric} scalar
#' 
#' @slot FRi \link[base]{numeric} \link[base]{vector}, functional regression index
#' 
#' @slot FRwt \link[base]{numeric} \link[base]{vector}, functional regression weights
#' 
#' @slot xarg \link[base]{numeric} \link[base]{vector}, \eqn{x}-domain
#' 
#' @slot X \link[base]{numeric} \link[base]{matrix}, 
#' each row representing the tabulated, on a common grid, 
#' functional predictor values for each subject.
#' 
#' @export
setClass(Class = 'FRindex', slots = c(
  gam = 'gam',
  sign = 'numeric', # scalar
  FRi = 'numeric', # vector
  FRwt = 'numeric', # vector
  xarg = 'numeric', # vector
  X = 'matrix'
))


#' @title Functional Regression Index Model
#' 
#' @description 
#' Fits functional regression model and computes functional indices based on functional predictors in the fitted model
#' 
#' @param formula \link[stats]{formula}.
#' Left-hand-side is the name of the response variable. 
#' Right-hand-side is the name of the tabulated \link[base]{matrix} functional predictor
#' 
#' @param data \link[base]{data.frame}, including 
#' a response variable and a \link[base]{matrix} of tabulated functional predictor.
#' Each row of such \link[base]{matrix} represents functional predictor values for each subject tabulated on a fine grid 
#' and each column corresponds to the same point on the same fine grid common for all subjects
#' If the functional predictor is the \link[stats]{quantile} function,
#' then \code{data} is preferably the returned object of \link{clusterQp}.
#' 
#' @param xarg \link[base]{numeric} \link[base]{vector},
#' \eqn{x}-domain to be passed to \link{FR_gam}
#' corresponding to the common grid on which the functional predictor is tabulated
#' 
#' @param sign_prob \link[base]{numeric} scalar,
#' the argument for selecting the nearest-even quantile in \code{xarg}, 
#' which is used to define the \link[base]{sign} of the weight function \code{@@FRwt} for \code{@@FRi}.
#' Default is \code{.5} corresponding to the nearest-even \link[stats]{median}
#' 
#' @param ... additional parameters of \link{FR_gam}
#' 
#' @details 
#' 
#' Functional regression index for each subject is defined as a functional predictor in the functional regression model
#' (integral of subject-specific functional predictor multiplied by the weight function common for all subjects).
#' The weight function \code{@@FRwt} is either equal or negated estimated functional coefficient in the functional regression model. 
#' The \link[base]{sign} of the weight function \code{@@FRwt} is selected so that the resulting \code{@@FRi} 
#' are positively associated with the values of functional predictor when \code{xarg} is equal to \code{sign_prob}.
#' 
#' @returns 
#' \link{FRindex} returns an \linkS4class{FRindex} object.
#' 
#' @examples 
#' 
#' library(survival)
#' 
#' pt = unique(Ki67$PATIENT_ID)
#' length(pt) # 622
#' # set.seed if needed
#' train_pt = sample(pt, size = 500L)
#' Ki67q = clusterQp(data = Ki67, exclude = c('tissueID','inner_x','inner_y'), contX = 'Marker')
#' train_q = subset(Ki67q, PATIENT_ID %in% train_pt)
#' test_q = subset(Ki67q, !(PATIENT_ID %in% train_pt))
#' train_q$Marker = log1p(train_q$Marker)
#' test_q$Marker = log1p(test_q$Marker)
#' 
#' FRi = FRindex(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = train_q)
#' FRi@@FRi # functional regression index
#' FRi@@FRwt # functional regression weights
#' 
#' (FRi_test = predict.FRindex(FRi, newdata = test_q$Marker))
#' 
#' FRi_train = predict.FRindex(FRi)
#' stopifnot(identical(FRi@@FRi, c(FRi_train)), 
#'   identical(FRi@@FRwt, attr(FRi_train, 'FRwt')))
#' 
#' # set.seed if needed
#' Ki67bbc_v2 = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, 
#'   data = data.frame(train_q, FRi_std = std_IQR(FRi_train)), 
#'   contX = 'FRi_std', R = 100L)
#' summary(Ki67bbc_v2)
#' 
#' @importFrom pracma cumtrapz
#' @importFrom stats cor quantile
#' @export
FRindex <- function(
    formula,
    data,
    xarg = as.numeric(colnames(X)),
    sign_prob = .5,
    ...
) {
  
  contX <- formula[[3L]]
  X <- data[[contX]]
  
  N <- length(xarg)
  if (!N || !is.numeric(xarg) || anyNA(xarg)) stop('`xarg` needs to be a numeric sequence')
  
  gam_obj <- FR_gam(y = eval(formula[[2L]], envir = data), X = X, xarg = xarg, ...)
  beta_ <- gam2beta(gam_obj, n = N)
  
  # integral from all patients
  intg_ <- cumtrapz(x = xarg, y = t.default(X) * beta_)[N, ]
  
  # medians of markers from all patients
  med_id <- which(xarg == quantile(xarg, probs = sign_prob, type = 3L))[1L]
  sign_ <- sign(cor(X[, med_id], intg_)) # scalar
  
  return(new(
    Class = 'FRindex', 
    gam = gam_obj,
    sign = sign_,
    FRi = sign_ * intg_, # functional regression index
    FRwt = sign_ * beta_, # functional regression weight
    xarg = xarg,
    X = X
  ))
  
}


#' @title Functional Regression Index Prediction
#' 
#' @description 
#' Computes functional regression indices for new data using the \linkS4class{FRindex} object
#' 
#' @param object an \linkS4class{FRindex} object
#' 
#' @param newdata \link[base]{matrix}, 
#' each row representing the tabulated, on a common grid, 
#' functional predictor values for each subject.
#' The number of values in the grid of \code{newdata}
#' does not need to be the same as that number used to obtain \code{object}.
#' If the functional predictor is the \link[stats]{quantile} function,
#' then \code{newdata} is preferably the \link[base]{matrix} column of 
#' the returned object of \link{clusterQp}.
#' 
#' @param new_xarg \link[base]{numeric} \link[base]{vector}, \eqn{x}-domain
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' Functional regression index for each subject in \code{newdata} is defined as 
#' an integral of subject-specific functional predictor multiplied by the weight function common for all subjects.
#' The weight function \code{@@FRwt} is obtained from \code{object@@gam} using function \link{gam2beta}.
#' 
#' @returns 
#' \link{predict.FRindex} returns a \link[base]{numeric} \link[base]{vector}, 
#' which is the functional regression index of \code{newdata}.
#' 
#' @examples 
#' # see ?`Qindex-package`
#' 
#' @importFrom pracma cumtrapz
#' @export predict.FRindex
#' @export
predict.FRindex <- function(
    object, newdata = object@X,
    new_xarg = as.numeric(colnames(newdata)),
    ...
) {
  
  obj_xarg <- object@xarg
  
  new_N <- length(new_xarg)
  if (!new_N || !is.numeric(new_xarg) || anyNA(new_xarg)) stop('`new_xarg` needs to be a numeric sequence')
  if (!isTRUE(all.equal.numeric(min(obj_xarg), min(new_xarg))) ||
      !isTRUE(all.equal.numeric(max(obj_xarg), max(new_xarg)))) stop('New and old data should share the same domain')
  
  new_beta <- gam2beta(object@gam, n = new_N)
  
  new_intg <- cumtrapz(x = new_xarg, y = t.default(newdata) * new_beta)[new_N, ]
  
  new_FRi <- object@sign * new_intg
  attr(new_FRi, which = 'FRwt') <- object@sign * new_beta
  
  return(new_FRi)
  
}




