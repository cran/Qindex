
#' @title Bootstrap Indexes
#' 
#' @description
#' Generate a *series of* bootstrap indexes.
#' 
#' @param N positive \link[base]{integer} scalar, sample size
#' 
#' @param ... additional parameters of \link[boot]{boot}, most importantly the number of bootstrap replicates \code{R}
#' 
#' @returns 
#' \link{bootIDX} returns a \link[base]{list} of \link[base]{integer} \link[base]{vector}s.
#' Each element is the indexes of a boot strap sample.
#' 
#' @examples
#' bootIDX(10L, R = 5L)
#' 
#' @importFrom boot boot
#' @export
bootIDX <- function(N, ...) {
  m <- boot(data = seq_len(N), statistic = function(data, ind) ind, ...)[['t']]
  seqR <- seq_len(dim(m)[1L]) # seq_len(R); `R` as in ?boot::boot 
  lapply(seqR, FUN = function(i) m[i,])
}