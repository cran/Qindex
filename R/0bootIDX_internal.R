
#' @title Generate Bootstrap Indices
#' 
#' @description
#' Generate a series of \link[boot]{boot}strap indices.
#' 
#' @param n positive \link[base]{integer} scalar, sample size
#' 
#' @param R positive \link[base]{integer} scalar, number of bootstrap replicates
#' 
#' @returns 
#' Function [bootIDX()] returns a \link[base]{list} of positive \link[base]{integer} \link[base]{vector}s.
#' Each element is the indices of one bootstrap sample.
#' 
#' @details
#' Function [bootIDX()] is designed to generate the same bootstrap indices as 
#' from the default options of function \link[boot]{boot}, 
#' given the same \link[base]{Random} seed.  
#' 
#' See details in `boot:::index.array()` and `boot:::ordinary.array()`.
#' 
#' @examples
#' set.seed(1345); boot::boot(data = 1:10, statistic = function(data, ind) ind, R = 3L)[['t']]
#' set.seed(1345); bootIDX(10L, R = 3L) # same copies of indices
#' 
#' @keywords internal
#' @export
bootIDX <- function(n, R) {
  
  # wont be the same !!
  # probably equivalent to `dim(ret) <- c(n, R)`
  #replicate(n = R, expr = {
  #  sample.int(n = n, replace = TRUE)
  #}, simplify = simplify)
  
  ret <- sample.int(n = n, size = n * R, replace = TRUE)
  dim(ret) <- c(R, n) # Tingting will prefer c(n, R)
  lapply(seq_len(R), FUN = function(i) ret[i,])
  
}



