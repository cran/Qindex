

#' @title Visualize \linkS4class{Qindex} object using R package \pkg{graphics}
#' 
#' @description
#' Create \link[graphics]{persp}ective and \link[graphics]{contour}
#' plots of FR-index integrand using R package \pkg{graphics}.
#' 
#' End users are encouraged to use function [integrandSurface]
#' with \CRANpkg{plotly} work horse.
#' 
#' @param x \linkS4class{Qindex} object
#' 
#' @param n \link[base]{integer} scalar, fineness of visualization,
#' default `501L`. See parameter `n.grid` of function \link[mgcv]{vis.gam}.
#' 
#' @param xlab,ylab \link[base]{character} scalars
#' 
#' @param zlab \link[base]{character} scalar, for function [persp.Qindex]
#' 
#' @param image_col argument `col` of \link[graphics]{image.default}
#' 
#' @param ... ..
#' 
#' @returns
#' Function [persp.Qindex], 
#' a method dispatch of S3 generic \link[graphics]{persp},
#' does not have a return value.
#' 
#' @keywords internal
#' @name Qindex_graphics
#' @importFrom graphics persp
#' @export persp.Qindex
#' @export
persp.Qindex <- function(
    x, 
    n = 31L, 
    xlab = 'Percentages',
    ylab = 'Quantiles',
    zlab = 'Integrand of FR-index',
    ...
) {
  
  z <- z_Qindex(x, n = n)
  # ?graphics:::persp.default
  persp(x = attr(z, which = 'xy', exact = TRUE),
        z = z, 
        xlab = xlab, ylab = ylab, zlab = zlab,
        ...)
  
  return(invisible()) # ?graphics:::persp.default has an invisible return!
}





#' @returns
#' Function [contour.Qindex],
#' a method dispatch of S3 generic \link[graphics]{contour},
#' does not have a return value
#' 
#' @rdname Qindex_graphics
#' @importFrom graphics contour contour.default image.default
#' @importFrom grDevices topo.colors
#' @export contour.Qindex
#' @export
contour.Qindex <- function(
    x, 
    n = 501L,
    image_col = topo.colors(20L),
    xlab = 'Percentages',
    ylab = 'Quantiles',
    ...
) {
  z <- z_Qindex(x, n = n)
  xy <- attr(z, which = 'xy', exact = TRUE)
  
  image.default(
    x = xy, z = z, 
    col = image_col, xlab = xlab, ylab = ylab, ...
  )
  
  contour.default(x = xy, z = z, add = TRUE, ...)
  
  return(invisible())
}



if (FALSE) {
  # mgcv::vis.gam
  debug(mgcv::vis.gam)
  undebug(mgcv::vis.gam)
  #vis.gam(nlfr@gam, view = c("xgrid", "v"), n.grid = 11L, plot.type = "contour", color = "topo")
  mgcv::vis.gam(nlfr@gam, view = c("xgrid", "marker"), n.grid = 101L, plot.type = "contour", color = "topo")
  mgcv::vis.gam(nlfr@gam, view = c("xgrid", "marker"), n.grid = 101L, plot.type = "persp", color = "topo")
  persp(nlfr)
}




#' @importFrom mgcv predict.gam
z_Qindex <- function(
    x, # returned object from [Qindex]
    n = 501L, 
    ...
) {
  
  rhs <- x@formula[[3L]]
  X <- x@gam$data[[rhs]] 
  xgrid <- as.double(colnames(X))
  nxgrid <- length(xgrid)
  
  # inspired by ?mgcv::vis.gam
  xy <- list(
    x = seq.int(from = min(xgrid), to = max(xgrid), length.out = n),
    y = seq.int(from = min(X), to = max(X), length.out = n)
  )
  d_surface <- data.frame(
    expand.grid(xy), # span `x` first, then span `y`
    L = 1/nxgrid
  )
  names(d_surface)[1:2] <- c('xgrid_', as.character(rhs))
  
  z0 <- x@sign * predict.gam(x@gam, newdata = d_surface, se.fit = FALSE, type = 'link')
  z <- array(z0, dim = c(n, n), dimnames = NULL)
  attr(z, which = 'xy') <- xy
  return(z)
  
}


