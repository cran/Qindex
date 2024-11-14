

#' @title Back Compatibility
#' 
#' @description
#' Functions that have been \link[base]{.Defunct}.
#' 
#' @param ... parameters that have been \link[base]{.Defunct}.
#' 
#' @keywords internal
#' @name defunct
#' @export 
FRindex <- function(...) Qindex_defunct(new = Qindex)
  
#' @rdname defunct
#' @export 
predict.FRindex <- function(...) Qindex_defunct(new = predict.Qindex)
  
#' @rdname defunct
#' @export 
optim_splitSample_dichotom <- function(...) Qindex_defunct(new = optimSplit_dichotom)






Qindex_defunct <- function(new) {
  if (!is.function(new)) stop('`new` must be an existing function')
  .Defunct(msg = sprintf(fmt = '** check out updated examples at ?Qindex::`Qindex-package` **\n  this function is replaced by %s()', deparse1(substitute(new))))
}
