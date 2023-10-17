
#' @title Extract Functional Coefficient from \link[mgcv]{gam} Return
#' 
#' @description
#' Extract the function coefficient from \link[mgcv]{gam} return
#' 
#' @param object a \link[mgcv]{gam} object
#' 
#' @param ... additional parameters of \link[mgcv]{plot.gam}, 
#' most importantly `n`
#' 
#' @details
#' Extract beta function from \link[mgcv]{gam} return, using \link[mgcv]{plot.gam}.
#' 
#' @note
#' Suppress the figure from printing; inspired by `oddsratio::no_plot`
#' 
#' @returns
#' Function [gam2beta()] returns a \link[base]{double} \link[base]{vector}.
#' 
#' @examples
#' library(mgcv)
#' 
#' 
#' # one dimensional
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
#' vis.gam(b)
#' gam2beta(b)
#' 
#' 
#' # two dimensional
#' # ?mgcv::ti examples
#' test1 <- function(x,z,sx=0.3,sz=0.4) { 
#' x <- x*20
#' (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
#' }
#' n <- 500
#' x <- runif(n)/20;z <- runif(n);
#' xs <- seq(0,1,length=30)/20;zs <- seq(0,1,length=30)
#' pr <- data.frame(x=rep(xs,30),z=rep(zs,rep(30,30)))
#' truth <- matrix(test1(pr$x,pr$z),30,30)
#' f <- test1(x,z)
#' y <- f + rnorm(n)*0.2
#' 
#' b2 <- gam(y ~ te(x,z))
#' vis.gam(b2)
#' dim(gam2beta(b2)) # works
#' 
#' 
#' b3 <- gam(y ~ ti(x) + ti(z) + ti(x,z))
#' vis.gam(b3)
#' plot(b3)
#' # debug(gam2beta); gam2beta(b3) # not working yet!
#' 
#' 
#' 
#' ## now illustrate partial ANOVA decomp...
#' b4 <- gam(y~ ti(x) + ti(x,z,mc=c(0,1))) ## note z constrained!
#' vis.gam(b4)
#' # debug(gam2beta); gam2beta(b4) # not working yet!
#' 
#' 
#' 
#' @importFrom grDevices png dev.off
#' @importFrom mgcv plot.gam
#' @keywords internal
#' @export
gam2beta <- function(object, ...) {
  png(tmp_figure <- file.path(tempdir(), 'tmp.png'))
  fig <- plot.gam(object, ...)
  fig1 <- fig[[1]] # not sure why ?mgcv::plot.gam organize the return this way, yet
  
  ret0 <- fig1$fit
  if (!is.matrix(ret0)) stop('mgcv package update?')
  if (dim(ret0)[2L] != 1L) stop('mgcv package update?')

  nx <- length(fig1[['x']])
  ny <- length(fig1[['y']])
  ret <- if (nx & !ny) {
    # one-dimensional
    c(ret0) # ncol-1 matrix to vector
  } else if (nx & ny) {
    # two-dimensional
    array(c(ret0), dim = c(nx, ny)) # ncol-1 matrix to `nx`-by-`ny` matrix
  }
  
  dev.off()
  file.remove(tmp_figure)
  return(ret)
}

