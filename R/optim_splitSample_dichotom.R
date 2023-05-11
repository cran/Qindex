

#' @title Selection of an Optimal Dichotomous Predictor via Repeated Sample Splits
#' 
#' @description
#' 
#' Repeated sample splits are used to (i) determine the optimal cutoff for dichotomizing 
#' each predictor of interest in the training split part and (ii) compute the corresponding
#' effect size in the test split part. 
#' For each predictor, the effect sizes in test splits are computed for 100-200 splits,
#' and the optimal predictor is selected with the highest median effect size.
#' 
#' @param formula \link[stats]{formula}. 
#' Left-hand-side is the name of a \link[survival]{Surv}, \link[base]{logical}, or \link[base]{numeric} endpoint.
#' Right-hand-side is the name of a \link[base]{matrix} predictor. 
#' See details of parameter \code{data}
#' 
#' @param data \link[base]{data.frame}, with at least 
#' \itemize{
#' \item {
#' two columns including the \link[base]{numeric} time-to-event and the \link[base]{logical} event indicator,
#' or one \link[survival]{Surv} column,
#' or one \link[base]{logical} column, 
#' or one \link[base]{numeric} column, as the endpoint}
#' \item {one \link[base]{numeric} \link[base]{matrix} column, 
#' the columns of which are a collection of predictors to be selected from}
#' }
#' 
#' @param include \link[base]{language} object, 
#' inclusion criteria for restricting the predictor(s) to guarantee
#' a user-desired range of proportions in \code{highX}.
#' A suggested choice is \code{(highX > .15 & highX < .85)}.
#' See explanation of \code{highX} in \link{splitSample_dichotom}.
#' 
#' @param top \link[base]{integer} scalar, number of optimal predictors, default \code{3L}
#' 
#' @param nsplit,... additional parameters for \link{stratifiedSplitSample}
#' 
#' @details 
#' 
#' To select an optimal dichotomous predictor via repeated sample splits,
#' 
#' \enumerate{
#' 
#' \item generates \code{nsplit} different sample splits
#' (via \link{stratifiedSplitSample}) and 
#' stores them in variable \code{splitIDs}
#' 
#' \item finds the median effect size for multiple sample splits
#' (via \link{median_splitSample_dichotom}) for each predictor, 
#' using the same sample splits \code{splitIDs}.
#' In this step, we obtain one \link{splitSample_dichotom} object for each predictor
#' 
#' \item (optional) selects a subset of 
#' \link{splitSample_dichotom} objects obtained from Step 2 (one for each predictor),
#' based on the inclusion criteria specified in parameter \code{include}.
#' Typically, we restrict the predictor(s) to guarantee
#' a user-desired range of proportions for \code{highX} (see *Value* section of \link{splitSample_dichotom}).
#' A suggested choice is \code{(highX > .15 & highX < .85)}.
#' 
#' \item ranks the \link{splitSample_dichotom} objects, 
#' obtained from Step 2 or Step 3,
#' by the decreasing order of the absolute values (\link[base]{abs}) of 
#' the regression coefficients of dichotomized predictors \code{attr(, 'coef')}.
#' The predictor(s) corresponding to the \link{splitSample_dichotom} object(s) 
#' with largest absolute value(s) of the regression coefficients of dichotomized predictors,
#' are referred to as the optimal predictor(s).
#' 
#' }
#' 
#' @returns 
#' \link{optim_splitSample_dichotom} returns a \link[base]{data.frame}, 
#' containing only the optimal predictors, with \link[base]{attributes}
#' \describe{
#' \item{\code{attr(,'top')}}{a \link[base]{data.frame} with 
#' regression coefficients of dichotomized predictors, 
#' branches, 
#' and proportion \code{highX}}
#' }
#' 
#' 
#' @examples 
#' # see ?`Qindex-package`
#' 
#' @export
optim_splitSample_dichotom <- function(
    formula, data,
    #include = (highX > .15 & highX < .85), # ?devtools::check warning
    include,
    top = 3L,
    nsplit,
    ...
) {
  
  yval <- eval(formula[[2L]], envir = data)
  splitIDs <- stratifiedSplitSample(yval, nsplit = nsplit, ...) # using same split for all predictors
  
  xval <- eval(formula[[3L]], envir = data) # 'matrix' of predictors
  if (anyNA(xval)) stop('do not allow missingness in the \'', deparse1(formula[[3L]]), '\' matrix; for now')
  
  tmp <- lapply(seq_len(dim(xval)[2L]), FUN = function(p) {
    median_splitSample_dichotom(y = yval, x = xval[,p], splitIDs = splitIDs)
  })
  mssd <- do.call(what = Map, args = c(list(f = c), lapply(tmp, FUN = function(i) attributes(i)[c('branch', 'highX', 'coef')])))
  
  if (!missing(include)) {
    id_excl <- !eval(call(name = 'with.default', data = quote(mssd), expr = substitute(include)))
    mssd$coef[id_excl] <- mssd$branch[id_excl] <- mssd$highX[id_excl] <- NA_real_
  }
  
  #stop('use model type')
  anymod <- tmp[[which(lengths(tmp) > 0L)[1L]]]
  
  if (inherits(anymod, what = 'coxph')) {
    mssd$HazardsRatio <- exp(mssd$coef)
  } else if (inherits(anymod, what = 'glm')) {
    mssd$OddsRatio <- exp(mssd$coef)
  } else if (inherits(anymod, what = 'lm')) {
    # do nothing
  }
  
  id_top <- order(abs(mssd$coef), decreasing = TRUE)[seq_len(top)]
  if (anyNA(mssd$coef[id_top])) stop('Decrease `top` (containing coef\'s which do not satisfy `include`)')
  
  data_top <- data
  xval_top <- xval[, id_top]
  data_top[[formula[[3L]]]] <- xval_top
  
  id_top_named <- id_top
  names(id_top_named) <- colnames(xval)[id_top]
  data_top[[paste0('d_', deparse1(formula[[3L]]))]] <- do.call(cbind, args = lapply(id_top_named, FUN = function(i) {
    branch_ <- mssd$branch[[i]]
    branch_[[2L]] <- quote(xval[, i])
    eval(branch_)
  }))
  
  mssd$branch <- vapply(mssd$branch, FUN = deparse1, FUN.VALUE = '')
  attr(data_top, which = 'top') <- as.data.frame.list(mssd, row.names = colnames(xval))[id_top, ]
  return(data_top)
  
}




# only for developers

# if (FALSE) {
# # optim_splitSample_dichotom(TStage ~ Marker, data = Ki67q, nsplit = 20) # next task
# }
# 
# if (FALSE) {
# colnames(Ki67_opt$Marker) = paste0('Ki67_', colnames(Ki67_opt$Marker))
#
# PRq = clusterQp(data = PR, contX = 'Marker')
# PR_eval = optim_splitSample_dichotom(Surv(RECFREESURV_MO, RECURRENCE) ~ Marker, data = PRq, nsplit = 20)
# PR_opt = optim_splitSample_dichotom(PR_eval, n = 3L)
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
#   contX = 'Marker', R = 100L)
# names(mod0)
# }
# 
