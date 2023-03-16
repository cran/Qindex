

#' @title Qindex: Continuous and dichotomized index predictors based on distribution quantiles
#' 
#' @description
#' 
#' Select optimal functional or dichotomized quantile-based predictors for survival/logistic/numeric outcome and 
#' perform optimal dichotomization with optimistic bias correction for any continuous predictors 
#'
#' @details
#' 
#' The package provides tools to
#' 
#' \enumerate{
#' 
#' \item use \link{sampleQp} to calculate user-selected sample quantiles in each independent 
#' cluster of observations.  This function is simply a wrapper of \link[stats]{aggregate}
#' for \link[stats]{quantile} function.
#' 
#' \item use \link{eval_dichotom} to estimate the effect size for dichotomized user-selected predictors of interest (e.g. sample quantiles)
#' (absolute value of the corresponding hazard ratio, odds ratio or regression coefficient) using repeated random split train/test sampling.
#' In this function, we first use \link[rpart]{rpart} to identify optimal cutoff in each training set 
#' and use this cutoff to dichotomize each predictor of interest in the corresponding independent test set.
#' The effect size for dichotomized predictor is estimated in the test set by fitting 
#' \link[survival]{coxph}, \link[stats]{glm} or \link[stats]{lm} to fit a Cox proportional hazard model,
#' logistic regression, or linear regression
#' for \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint.
#' 
#' \item use \link{optQp} to select the set of optimal quantiles that 
#' has the largest effect size (absolute value of the corresponding hazard ratio, odds ratio or regression coefficient)
#' for a given \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint.  
#' \link{optQp} is a wrapper of \link{summary.eval_dichotom}.
#' 
#' \item use \link{BBC_dichotom} to dichotomize predictors of interest and to obtain bootstrap-based optimism corrected effect size
#' from Cox model, logistic regression, or linear regression.  
#' Internally, \link{BBC_dichotom} calls \link{dichotom_int} to dichotomize each predictor in \code{contX} based on univariate model setting
#' and \code{model_dichotom} to fit Cox proportional hazard model, logistic regression, or linear regression 
#' for \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint
#' with the dichotomized predictors from \link{dichotom_int}.
#' 
#' \item use \link{FRindex} to derive a scalar functional regression index as a predictor in the functional regression model
#' with any response supported by \CRANpkg{gam} package. 
#' Function \link{FRindex} calls \link[mgcv]{gam} to fit a generalized additive model(GAM) to the training set
#' and makes use of \link[mgcv]{plot.gam} to extract functional coefficient tabulated on the same grid as functional predictor(s)
#' in the training and test set (if the test set is provided).
#' }
#'
#' @references 
#' Selection of optimal quantile protein biomarkers based on cell-level immunohistochemistry data,
#' Misung Yi, Tingting Zhan , Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva, under review
#' 
#' Quantile index biomarkers based on single-cell expression data,
#' Misung Yi, Tingting Zhan , Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva, under review
#' 
#' @import stats utils
#' @importFrom boot boot
#' @importFrom grDevices png dev.off
#' @importFrom matrixStats colMedians
#' @importFrom mgcv gam cox.ph plot.gam
#' @importFrom rpart rpart rpart.control
#' @importFrom survival coxph Surv
#' 
#' @docType package
#' @keywords package
#' @name Qindex-package
NULL



# @importFrom mgcv s  
