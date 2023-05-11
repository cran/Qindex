

#' @title Continuous and Dichotomized Index Predictors Based on Distribution Quantiles
#' 
#' @description
#' 
#' Select optimal functional regression or dichotomized quantile predictors for survival/logistic/numeric outcome and 
#' perform optimistic bias correction for any optimally dichotomized \link[base]{numeric} predictor(s)
#'
#' @details
#' 
#' Primary functions in this package are
#' 
#' \describe{
#' 
#' \item{\link{clusterQp}}{calculate user-selected sample quantiles in each cluster of observations.}
#' 
#' \item{\link{optim_splitSample_dichotom}}{optimal predictor selection via dichotomizing split sample}
#' 
#' \item{\link{BBC_dichotom}}{Bootstrap-based optimism correction for dichotomizing selected \link[base]{numeric} predictor(s)}
#' 
#' \item{\link{FRindex}}{Functional regression index as a predictor in the functional regression model}
#' 
#' }
#'
#' @references 
#' Selection of optimal quantile protein biomarkers based on cell-level immunohistochemistry data,
#' Misung Yi, Tingting Zhan , Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva.  Under revision
#' 
#' Quantile index biomarkers based on single-cell expression data,
#' Misung Yi, Tingting Zhan , Amy P. Peck, Jeffrey A. Hooke, Albert J. Kovatich, Craig D. Shriver, 
#' Hai Hu, Yunguang Sun, Hallgeir Rui and Inna Chervoneva. 
#' Lab Invest, 2023 Apr 21;100158. \doi{10.1016/j.labinv.2023.100158}. Online ahead of print.
#' 
#' 
#' @examples
#' 
#' library(survival)
#' 
#' Ki67q = clusterQp(data = Ki67, exclude = c('tissueID','inner_x','inner_y'), contX = 'Marker')
#' stopifnot(!anyDuplicated.default(Ki67q$subjID))
#' head(Ki67q)
#' sapply(Ki67q, FUN = class)
#'
#' # set.seed if needed
#' Ki67c_toy = optim_splitSample_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = Ki67q, 
#'   nsplit = 5L, include = (highX > .15 & highX < .85), top = 2L)
#' head(Ki67c_toy)   
#' attr(Ki67c_toy, 'top') # under-the-hood statistics for the optimal predictors 
#' \donttest{
#' # slow
#' Ki67c = optim_splitSample_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = Ki67q, 
#'   nsplit = 20L, include = (highX > .15 & highX < .85), top = 2L) 
#' }
#' 
#' # set.seed if needed
#' Ki67bbc_toy = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, 
#'   data = Ki67c_toy, contX = 'Marker', R = 100L)
#' summary(Ki67bbc_toy)
#' attr(Ki67bbc_toy, 'median_optimism')
#' attr(attr(Ki67bbc_toy, 'median_optimism'), 'boot_branch')
#' attr(Ki67bbc_toy, 'apparent_branch')
#' \donttest{
#' Ki67bbc = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ NodeSt + Tstage, 
#'   data = Ki67c, contX = 'Marker', R = 100L)
#' }
#' 
#' # see more examples in ?FRindex
#' 
#' @import methods
#' 
#' @docType package
#' @keywords package
#' @name Qindex-package
NULL



# @importFrom mgcv s  


# if (FALSE) { # for developers for now
# Ki67a = optim_splitSample_dichotom(RECFREESURV_MO ~ Marker, data = Ki67q, 
#  nsplit = 20L, top = 2L)
# set.seed(1)
# mod_a = BBC_dichotom(RECFREESURV_MO ~ NodeSt + Tstage, data = Ki67a, 
#   contX = 'Marker', R = 100L)
# summary(mod_a)
# 
# Ki67b = optim_splitSample_dichotom(RECURRENCE ~ Marker, data = Ki67q, nsplit = 20L, top = 2L)
# set.seed(1)
# mod_b = BBC_dichotom(RECURRENCE ~ NodeSt + Tstage, data = Ki67b, 
#   contX = 'Marker', R = 100L)
# summary(mod_b)
# }
# 



# if (FALSE) { # for developers for now
# FRindex(trainData = train_q, testData = test_q, response = RECURRENCE)
# FRindex(trainData = train_q, testData = test_q, response = RECFREESURV_MO)
# 
# set.seed(1)
# mod = BBC_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ 1, data = FRQI$trainData, 
#   contX = 'FRindex_std', R = 200L)
# summary(mod)
# names(attributes(mod))
# attr(mod, 'median_optimism')
# }





#' @title Ki67 Data
#' 
#' @description 
#' Ki67 cell data containing 622 patients
#' 
#' @format
#' \describe{
#'   \item{\code{PATIENT_ID}}{\link[base]{factor}, unique patient identifier}
#'   \item{\code{tissueID}}{\link[base]{factor}, TMA core identifier}
#'   \item{\code{RECURRENCE}}{\link[base]{integer}, recurrence indicator, 1 = Recurred, 0 = not Recurred}
#'   \item{\code{RECFREESURV_MO}}{\link[base]{integer}, recurrence-free survival time in months}
#'   \item{\code{Marker}}{\link[base]{numeric}, cell signal intensity of the protein immunofloerscence signal}
#'   \item{\code{inner_x}}{\link[base]{integer}, \eqn{x}-coordinate in the cell centroid in the TMA core}
#'   \item{\code{inner_y}}{\link[base]{integer}, \eqn{y}-coordinate in the cell centroid in the TMA core}
#'   \item{\code{AGE_AT_DX}}{\link[base]{integer}, age at diagnosis}
#'   \item{\code{Tstage}}{\link[base]{integer}, tumor stage}
#'   \item{\code{NodeSt}}{\link[base]{integer}, node stage, -1 = unknown, 0 = Node Negative, 1 = Node Positive}
#'   \item{\code{HRpos}}{\link[base]{integer}, indicator of hormone positive status (ER+ or PR+), 1 = positive, 0 = negative}
#'   \item{\code{HistologicalGrade}}{\link[base]{integer}, histology grade}
#'   \item{\code{Her2_path_qIF}}{\link[base]{integer}, Her2 status, 1 = positive, 0 = negative}
#'   \item{\code{RACE}}{\link[base]{character}, race, White, Black, Asian, Native Hawaiian or Other Pacific Islander, American Indian or Alaska Native, Unknown}
#'   \item{\code{RadjCHEMO}}{\link[base]{integer}, adjuvant chemo treatment, 0 = unknown,  1 = done, 2 = NOT done}
#'   \item{\code{RadjRAD}}{\link[base]{integer}, adjuvant radiation treatment, 0 = unknown,  1 = done, 2 = NOT done}
#'   \item{\code{HORM_4cat}}{\link[base]{integer}, hormone treatment, 0 = unknown, 1 = not indicated, 2 = done, 3 = recommended, but not done}
#'   \item{\code{MSI}}{\link[base]{numeric}, mean signal intensity (mean over all cells in the TMA core)}
#' }
#' 
#' @name celldata
'Ki67'

