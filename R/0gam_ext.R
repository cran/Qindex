
#' @title Using \link[mgcv]{gam} for Functional Regression
#' 
#' @description 
#' Estimate the functional coefficient by fitting functional regression model (via \link[mgcv]{gam}). 
#' 
#' @param y either a \link[base]{numeric} \link[base]{vector}, 
#' a \link[base]{logical} \link[base]{vector},
#' or a \link[survival]{Surv} object.
#' 
#' @param X \link[base]{numeric} \link[base]{matrix} of tabulated functional predictor.
#' Each row represents the tabulated, on a common grid, 
#' functional predictor values for each subject.
#' Each column corresponds to the same point on the common grid.
#' 
#' @param xarg \link[base]{numeric} \link[base]{vector}, default to 
#' \code{as.numeric(colnames(X))}
#' 
#' @param family \link[stats]{family} object specifying the distribution 
#' and link to use in \link[mgcv]{gam}
#' 
#' @param knot_pct \link[base]{numeric} scalar, 
#' percentage of the number of columns of \code{X},
#' to be used as \code{knot.value}.  
#' Default is \eqn{40\%}.
#' If \code{knot.value} is provided by the end-user, then \code{knot_pct} is ignored.
#' 
#' @param knot.value \link[base]{integer} scalar, number of knots 
#' (i.e., parameter \code{k} in the spline smooth function \link[mgcv]{s})
#' used in \link[mgcv]{gam}.
#' Default is the \link[base]{ceiling} of \code{knot_pct} of
#' the column dimension of the \link[base]{matrix} column \code{contX}.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#'
#' Using \link[mgcv]{gam} to estimate the functional coefficient by fitting functional regression model. 
#'
#' 
#' @returns
#' \link{FR_gam} returns a \link[mgcv]{gam} object
#' 
#' @references 
#' 
#' Cui, E., Crainiceanu, C. M., & Leroux, A. (2021). Additive Functional Cox Model. Journal of Computational and Graphical Statistics. 2021;30(3):780-793.
#' \doi{10.1080/10618600.2020.1853550}
#' 
#' Gellar, J. E., Colantuoni, E., Needham, D. M., & Crainiceanu, C. M. (2015). Cox regression models with functional covariates for survival data. Statistical Modelling, 15(3), 256-278.
#' \doi{10.1177/1471082X14565526}
#' 
#' @examples 
#' Ki67q = clusterQp(data = Ki67, 
#'   exclude = c('tissueID','inner_x','inner_y'), contX = 'Marker')
#' Ki67q$Marker = log1p(Ki67q$Marker)
#' 
#' yy = FR_gam(
#'  y = with(Ki67q, expr = survival::Surv(RECFREESURV_MO, RECURRENCE)),
#'  X = Ki67q[['Marker']])
#' 
#' @importFrom mgcv gam cox.ph s
#' @importFrom stats binomial gaussian
#' @export
FR_gam <- function(
    y, 
    X,
    xarg = as.numeric(colnames(X)),
    family,
    knot_pct = .4,
    knot.value = ceiling(length(xarg) * knot_pct), 
    ...
) {
  
  dm <- dim(X)
  # stopifnot(dm[2L] == length(xarg))
  
  ### `Tm`: time indices
  Tm <- tcrossprod(rep(1, times = dm[1L]), xarg)
  # identical to `t.default(array(xarg, dim = c(dm[2L], dm[1L])))`, which is much slower
  
  L <- array(1/dm[2L], dim = dm)
  s_term <- quote(s(Tm, by = L * X, bs = 'cr', k = knot.value))
  # original
  
  #s_term <- quote(s(Tm, by = X / dm[2L], bs = 'cr', k = knot.value)) 
  # error. why?
  
  #by_ <- array(1/dm[2L], dim = dm) * X
  #s_term <- quote(s(Tm, by = by_, bs = 'cr', k = knot.value))
  # not identical to original.  why?
  
  gam_cl <- if (inherits(y, what = 'Surv')) {
    call('gam', 
         formula = call('~', quote(y[,1L]), s_term),
         weights = y[,2L], 
         family = if (missing(family)) quote(cox.ph()) else substitute(family))
  } else {
    call('gam', formula = call('~', quote(y), s_term), 
         family = if (!missing(family)) {
           substitute(family) 
         } else if (is.logical(y) || all(y %in% c(0, 1))) {
           quote(binomial(link = 'logit'))
         } else if (is.numeric(y)) {
           quote(gaussian(link = 'identity'))
         } else stop('not supported yet'))
  }
  
  return(eval(gam_cl))
  
}



# future: nonlinear functional regression gam
# nFR_gam



#' @title Extract Functional Coefficient from \link[mgcv]{gam} Return
#' 
#' @description
#' Extract the function coefficient from \link[mgcv]{gam} return
#' 
#' @param object a \link[mgcv]{gam} object
#' 
#' @param ... additional parameters of \link[mgcv]{plot.gam}, most importantly \code{n}
#' 
#' @details
#' Extract beta function from \link[mgcv]{gam} return, using \link[mgcv]{plot.gam}.
#' 
#' @note
#' Suppress the figure from printing; inspired by \code{oddsratio::no_plot}
#' 
#' @returns
#' \link{gam2beta} returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @examples
#' library(mgcv)
#' # ?mgcv::gam examples
#' set.seed(0)
#' f1 = function(x) {exp(2 * x)}
#' f2 = function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10 
#' f3 = function(x) {x*0}
#' n = 200
#' sig2 = 4
#' x0 = rep(1:4,50)
#' x1 = runif(n, 0, 1)
#' x2 = runif(n, 0, 1)
#' x3 = runif(n, 0, 1)
#' e = rnorm(n, 0, sqrt(sig2))
#' y = 2*x0 + f1(x1) + f2(x2) + f3(x3) + e
#' x0 = factor(x0)
#' b = gam(y~x0+s(x1)+s(x2)+s(x3))
#' plot(b)
#' gam2beta(b, n = 100)
#' 
#' @importFrom grDevices png dev.off
#' @importFrom mgcv plot.gam
#' @export
gam2beta <- function(object, ...) {
  png(tmp_figure <- file.path(tempdir(), 'tmp.png'))
  ret <- c(plot.gam(object, ...)[[1]]$fit)
  dev.off()
  file.remove(tmp_figure)
  return(ret)
}




if (FALSE) { # only for developer
  # why no element carries `length(xarg)` ??
  #id_mat = vapply(xx, inherits, what = 'matrix', FUN.VALUE = NA)
  #lapply(xx[id_mat], dim)
  #id_vec = vapply(xx, is.vector, FUN.VALUE = NA)
  #lengths(xx[id_vec])
} # only for developer