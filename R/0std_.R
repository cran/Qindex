
#' @title Standardization via \link[stats]{median} and \link[stats]{IQR}
#' 
#' @description
#' Standardization via \link[stats]{median} and \link[stats]{IQR}
#' 
#' @param x \link[base]{numeric} \link[base]{vector}
#' 
#' @return 
#' \link{std_IQR} returns a \link[base]{numeric} \link[base]{vector} of the same length as \code{x}
#' 
#' @examples
#' std_IQR(rnorm(20))
#' 
#' @importFrom stats quantile
#' @export
std_IQR <- function(x) {
  qs <- quantile(x, probs = c(.25, .5, .75), na.rm = TRUE)
  (x - qs[2L]) / (qs[3L] - qs[1L])
}


#' @title Standardization via \link[stats]{median} and \link[stats]{mad}
#' 
#' @description
#' Standardization via \link[stats]{median} and \link[stats]{mad}.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}
#' 
#' @return 
#' \link{std_mad} returns a \link[base]{numeric} \link[base]{vector} of the same length as \code{x}
#' 
#' @examples
#' std_mad(rnorm(20))
#' 
#' @importFrom stats median.default mad
#' @export
std_mad <- function(x) {
  m <- median.default(x, na.rm = TRUE)
  mad_ <- mad(x, center = m, na.rm = TRUE)
  (x - m) / mad_
}




