
#' @title Integrand Surface(s) of Sign-Adjusted Quantile Indices \linkS4class{Qindex}
#' 
#' @description
#' An interactive \CRANpkg{htmlwidgets} of the \link[graphics]{persp}ective plot for 
#' \linkS4class{Qindex} model(s)
#' using package \CRANpkg{plotly}.
#' 
#' @param ... one or more \linkS4class{Qindex} models
#' based on *a same training set*.
#' 
#' @param newdata \link[base]{data.frame}, with at least 
#' the response \eqn{y^{\text{new}}} and
#' the \link[base]{double} \link[base]{matrix} of 
#' functional predictor values \eqn{X^{\text{new}}}
#' of the *test set*. 
#' The predictor \eqn{X^{\text{new}}} are 
#' tabulated on the same \eqn{p}-grid as 
#' the training functional predictor values \eqn{X}.
#' If missing, the training set will be used.
#' 
#' @param proj_Q_p \link[base]{logical} scalar, whether to show 
#' the projection of \eqn{\hat{S}\big(p, Q_i(p)\big)}
#' (see sections **Details** and **Value**)
#' to the \eqn{(p,q)}-plain, default `TRUE`
#' 
#' @param proj_S_p \link[base]{logical} scalar, whether to show
#' the projection of \eqn{\hat{S}\big(p, Q_i(p)\big)} to the \eqn{(p,s)}-plain, default `TRUE`
#' 
#' @param proj_beta \link[base]{logical} scalar, whether to show
#' \eqn{\hat{\beta}(p)} on the \eqn{(p,s)}-plain when applicable, default `TRUE`
#' 
#' @param n \link[base]{integer} scalar, fineness of visualization,
#' default `501L`. See parameter `n.grid` of function \link[mgcv]{vis.gam}.
#' 
#' @param newid \link[base]{integer} scalar or \link[base]{vector},
#' row indices of `newdata` to be visualized. 
#' Default `1:2`, i.e., the first two test subjects.
#' Use `newid = NULL` to disable visualization of `newdata`.
#' 
#' @param qlim \link[base]{length}-2 \link[base]{double} \link[base]{vector},
#' range on \eqn{q}-axis. Default is the range of \eqn{X} and \eqn{X^{\text{new}}} combined.
#' 
#' @param axis_col \link[base]{length}-3 \link[base]{character} \link[base]{vector},
#' colors of the \eqn{(p,q,s)} axes
#' 
#' @param beta_col \link[base]{character} scalar, color 
#' of \eqn{\hat{\beta(p)}}
#' 
#' @param surface_col \link[base]{length}-2 \link[base]{character} \link[base]{vector},
#' color of the integrand surface(s), for lowest and highest surface values
#' 
#' @section Integrand Surface:
#' 
#' The quantile index (QI), 
#' \deqn{\text{QI}=\displaystyle\int_0^1\beta(p)\cdot Q(p)\,dp}
#' with a linear functional coefficient \eqn{\beta(p)}
#' can be estimated by fitting a functional generalized linear model (FGLM, James, 2002) to exponential-family outcomes, 
#' or by fitting a linear functional Cox model (LFCM, Gellar et al., 2015) to survival outcomes. 
#' More flexible non-linear quantile index (nlQI)
#' \deqn{
#' \text{nlQI}=\displaystyle\int_0^1 F\big(p, Q(p)\big)\,dp
#' }
#' with a bivariate twice differentiable function \eqn{F(\cdot,\cdot)}
#' can be estimated by fitting a functional generalized additive model (FGAM, McLean et al., 2014) to exponential-family outcomes, 
#' or by fitting an additive functional Cox model (AFCM, Cui et al., 2021) to survival outcomes. 
#' 
#' The estimated **integrand surface** of quantile indices and non-linear quantile indices, defined on 
#' \eqn{p\in[0,1]} and 
#' \eqn{q\in\text{range}\big(Q_i(p)\big)} for all training subjects \eqn{i=1,\cdots,n}, 
#' is
#' \deqn{
#' \hat{S}_0(p,q) = 
#' \begin{cases}
#' \hat{\beta}(p)\cdot q & \text{for QI}\\
#' \hat{F}(p,q) & \text{for nlQI}
#' \end{cases}
#' }
#' 
#' @section Sign-Adjustment:
#' 
#' Ideally, we would wish that, *in the training set*, the estimated linear and/or non-linear quantile indices
#' \deqn{
#' \widehat{\text{QI}}_i = \displaystyle\int_0^1 \hat{S}_0\big(p, Q_i(p)\big)dp
#' }
#' be *positively correlated* with a more intuitive quantity, e.g., quantiles \eqn{Q_i(\tilde{p})} at a user-specified \eqn{\tilde{p}}, for the interpretation of downstream analysis, 
#' Therefore, we define the sign-adjustment term
#' \deqn{
#' \hat{c} = \text{sign}\left(\text{corr}\left(Q_i(\tilde{p}), \widehat{\text{QI}}_i\right)\right),\quad i =1,\cdots,n
#' }
#' as the \link[base]{sign} of the \link[stats]{cor}relation between 
#' the estimated quantile index \eqn{\widehat{\text{QI}}_i}
#' and the quantile \eqn{Q_i(\tilde{p})},
#' for training subjects \eqn{i=1,\cdots,n}.
#' 
#' The estimated **sign-adjusted integrand surface** is
#' \eqn{\hat{S}(p,q) = \hat{c} \cdot \hat{S}_0(p,q)}.
#' 
#' The estimated **sign-adjusted quantile indices**
#' \eqn{\int_0^1 \hat{S}\big(p, Q_i(p)\big)dp}
#' are positively correlated with subject-specific sample medians
#' (default \eqn{\tilde{p} = .5}) in the training set.
#' 
#' 
#' @returns 
#' Function [integrandSurface] returns a pretty \CRANpkg{htmlwidgets} created by **R** package \CRANpkg{plotly}
#' to showcase the \link[graphics]{persp}ective plot of the
#' estimated sign-adjusted integrand surface \eqn{\hat{S}(p,q)}.
#' 
#' If a set of training/test subjects is selected (via parameter `newid`), then 
#' \itemize{
#' \item {the estimated **sign-adjusted line integrand curve** \eqn{\hat{S}\big(p, Q_i(p)\big)} 
#' of subject \eqn{i} 
#' is displayed on the surface \eqn{\hat{S}(p,q)};}
#' \item {the quantile curve \eqn{Q_i(p)} 
#' is projected on the \eqn{(p,q)}-plain of the 3-dimensional \eqn{(p,q,s)} cube, 
#' if `proj_Q_p=TRUE` (default);}
#' \item {the user-specified \eqn{\tilde{p}} is marked on the \eqn{(p,q)}-plain of the 3D cube, 
#' if `proj_Q_p=TRUE` (default);}
#' \item {\eqn{\hat{S}\big(p, Q_i(p)\big)}
#' is projected on the \eqn{(p,s)}-plain of the 3-dimensional \eqn{(p,q,s)} cube, 
#' if one and only one \linkS4class{Qindex} model is provided in in
#' put argument `...` and `proj_S_p=TRUE` (default);}
#' \item {the estimated *linear functional coefficient* \eqn{\hat{\beta}(p)} is shown on the \eqn{(p,s)}-plain of the 3D cube, 
#' if one and only one *linear* \linkS4class{Qindex} model is provided in input argument `...` and `proj_beta=TRUE` (default).}
#' }
#' 
#' @note
#' The maintainer is not aware of any functionality of projection of arbitrary curves in package \CRANpkg{plotly}.
#' Currently, the projection to \eqn{(p,q)}-plain is hard coded on \eqn{(p,q,s=\text{min}(s))}-plain.
#' 
#' @references 
#' James, G. M. (2002). *Generalized Linear Models with Functional Predictors*,
#' \doi{10.1111/1467-9868.00342} 
#' 
#' Gellar, J. E., et al. (2015). *Cox regression models with functional covariates for survival data*,
#' \doi{10.1177/1471082X14565526}
#' 
#' Mathew W. M., et al. (2014) *Functional Generalized Additive Models*,
#' \doi{10.1080/10618600.2012.729985}
#' 
#' Cui, E., et al. (2021). *Additive Functional Cox Model*,
#' \doi{10.1080/10618600.2020.1853550}
#' 
#' 
#' @examples 
#' # see ?`Qindex-package`
#' @importFrom mgcv predict.gam
#' @importFrom plotly plot_ly add_paths add_surface layout
#' @importFrom stats asOneSidedFormula predict
#' @export
integrandSurface <- function(
    ...,
    newdata = data,
    proj_Q_p = TRUE, 
    proj_S_p = TRUE,
    proj_beta = TRUE,
    n = 501L,
    newid = seq_len(min(50L, .row_names_info(newdata, type = 2L))), 
    qlim = range(X, newX),
    axis_col = c('dodgerblue', 'deeppink', 'darkolivegreen'),
    beta_col = 'purple',
    surface_col = 
      # c('lightyellow', 'lightpink') # nice
      # c('beige', 'lightpink') # nice
      # c('white', 'deeppink') # not good!
      # c('white', 'magenta') # not good!
      c('white', 'lightgreen') # nice
    # c('white', 'darkgreen') # not good!
    # c('white', 'lightgoldenrod') # my R do not recognize
    # c('white', 'lightslateblue') # my R do not recognize
    # c('white', 'yellow') # nice
) {
  
  xs <- list(...)
  if (!all(vapply(xs, FUN = inherits, what = 'Qindex', FUN.VALUE = NA))) stop('all input needs to be `Qindex`')
  
  rhs_ <- unique(lapply(xs, FUN = function(i) i@formula[[3L]]))
  if (length(rhs_) > 1L) stop('endpoints not same?')
  rhs <- rhs_[[1L]]
  
  data_ <- unique(lapply(xs, FUN = function(i) i@gam$data))
  if (length(data_) > 1L) stop('data not same')
  data <- data_[[1L]]
  
  sign_prob_ <- unique(lapply(xs, FUN = function(i) i@sign_prob))
  if (length(sign_prob_) > 1L) stop('@sign_prob not same')
  sign_prob <- sign_prob_[[1L]]
  
  X <- data[[rhs]]
  
  xgrid <- as.double(colnames(X))
  nxgrid <- length(xgrid)
  
  newX <- newdata[[rhs]]
  if (!is.matrix(newX)) stop('`newdata` does not contain a matrix column of functional predictor values')
  new_xgrid <- as.double(colnames(newX))
  if (!all.equal.numeric(new_xgrid, xgrid)) stop('grid of training and test data must be exactly the same')
  
  # plot!!
  # *surface* based on training model
  x_ <- seq.int(from = min(xgrid), to = max(xgrid), length.out = n)
  y_ <- seq.int(from = qlim[1L], to = qlim[2L], length.out = n)
  d_surface <- data.frame(
    expand.grid(xgrid_ = x_, y_), # span `x_` first, then span `y_`
    L = 1/nxgrid
  )
  names(d_surface)[2] <- c(as.character(rhs))
  zs <- lapply(xs, FUN = function(i) {
    y0 <- i@sign * predict.gam(i@gam, newdata = d_surface, se.fit = FALSE, type = 'link')
    t.default(array(y0, dim = c(n, n), dimnames = NULL))
    # ?base::t.default important here!!!
    # plot_ly(, type = 'surface') lay out `z` differently from ?graphics::persp !!!
  })
  
  zmin <- min(unlist(zs))
  zmax <- max(unlist(zs))
  
  p <- plot_ly(x = x_, y = y_)
  for (z_ in zs) {
    p <- add_surface(
      p = p, 
      z = z_, cmin = zmin, cmax = zmax, 
      # contours = list(z = list(show = TRUE)), # works :)
      colorscale = list(c(0, 1), surface_col), 
      showscale = FALSE)
  }
  p <- layout(p = p, scene = list(
    xaxis = list(title = 'Probability (p)', tickformat = '.0%', color = axis_col[1L]), 
    yaxis = list(title = 'Quantile (q)', color = axis_col[2L]),
    zaxis = list(title = 'Integrand (s)', color = axis_col[3L])
  ))
  
  if (!length(newid)) return(p)
  
  if (!is.integer(newid) || anyNA(newid) || any(newid > nrow(newX))) stop('illegal `newid`')
  
  d_subj <- data.frame(
    xgrid_ = xgrid,
    y = c(t.default(newX[newid, , drop = FALSE])),
    id = rep(newid, each = nxgrid),
    L = 1/nxgrid
  )
  names(d_subj)[2] <- as.character(rhs)
  z_subj <- lapply(xs, FUN = function(i) {
    i@sign * predict.gam(i@gam, newdata = d_subj, se.fit = FALSE, type = 'link')
  })
  
  if (proj_Q_p) {
    
    p <- add_paths(
    #p <- add_trace(
      p = p, data = d_subj, 
      x = ~ xgrid_, y = asOneSidedFormula(rhs), z = zmin, 
      name = ~ id, color = ~ id,
      showlegend = FALSE, # (length(newid) <= 10L),
      line = list(
        width = if (length(newid) <= 5L) 5/5 else 2/5
      ))
    
    p <- add_paths(
      p = p, data = d_subj, 
      x = sign_prob, y = asOneSidedFormula(rhs), 
      z = zmin, 
      showlegend = FALSE,
      line = list(
        width = 1,
        color = axis_col[1L] #,
        # dash = 'dot' # 'dash', 'dot'; # only works for ?plotly::add_trace
      ))
    
  } # projection on x-y plain: Q(p) curve
  
  if (proj_S_p) {
    # projection on x-z plain, F(p, Q(p)) curve
    # this is only done if (length(xs) == 1L); otherwise too messy
    if (length(xs) == 1L) {
      for (i in seq_along(xs)) {
        p <- add_paths(
          p = p, data = d_subj, 
          x = ~ xgrid_, y = qlim[2L], z = z_subj[[i]], 
          name = ~ id, color = ~ id,
          showlegend = FALSE, # (length(newid) <= 10L),
          line = list(
            width = if (length(newid) <= 5L) 5/5 else 2/5
          ))
      }
    }
  } # projection on x-z plain, F(p, Q(p)) curve
  
  if (proj_beta) {
    if (length(xs) == 1L) {
      for (i in seq_along(xs)) {
        if (xs[[i]]@gpf$formula[[3L]][[1L]] == 's') {
          
          d_beta <- data.frame(
            xgrid_ = xgrid,
            y = rep(1, times = nxgrid),
            id = rep(1, times = nxgrid),
            L = 1/nxgrid
          )
          names(d_beta)[2] <- as.character(rhs)
          z_beta <- lapply(xs, FUN = function(i) {
            i@sign * predict.gam(i@gam, newdata = d_beta, se.fit = FALSE, type = 'link')
          })
          
          p <- add_paths(
            p = p, data = d_beta, 
            x = ~ xgrid_, y = qlim[2L], z = z_beta[[i]], 
            showlegend = FALSE,
            line = list(
              color = beta_col,
              width = 4
            ))  
        }
      }
    }
  } # projection of beta(p), for linear models
  
  for (i in seq_along(xs)) {
    p <- add_paths(
      p = p, data = d_subj, 
      x = ~ xgrid_, y = asOneSidedFormula(rhs), 
      z = z_subj[[i]], 
      name = ~ id, color = ~ id,
      showlegend = FALSE, # (length(newid) <= 10L),
      line = list(
        width = if (length(newid) <= 5L) 5 else 2
      ))
  }
  
  
  
  #p <- layout(p = p, legend = list(
  #  title = list(text = if (identical(newdata, data)) 'Training Subj' else 'Test Subj'))
  #)
  return(p)
  
}






